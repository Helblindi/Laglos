#include "laglos_solver.hpp"
#include <cassert>


namespace mfem
{

namespace hydroLO
{

void VisualizeField(socketstream &sock, const char *vishost, int visport,
                    ParGridFunction &gf, const char *title,
                    int x, int y, int w, int h, bool vec)
{
   gf.HostRead();
   ParMesh &pmesh = *gf.ParFESpace()->GetParMesh();
   MPI_Comm comm = pmesh.GetComm();

   int num_procs, myid;
   MPI_Comm_size(comm, &num_procs);
   MPI_Comm_rank(comm, &myid);

   bool newly_opened = false;
   int connection_failed;

   do
   {
      if (myid == 0)
      {
         if (!sock.is_open() || !sock)
         {
            sock.open(vishost, visport);
            sock.precision(8);
            newly_opened = true;
         }
         sock << "solution\n";
      }

      pmesh.PrintAsOne(sock);
      gf.SaveAsOne(sock);

      if (myid == 0 && newly_opened)
      {
         const char* keys = (gf.FESpace()->GetMesh()->Dimension() == 2)
                            ? "mAcRjl" : "mmaaAcl";

         sock << "window_title '" << title << "'\n"
              << "window_geometry "
              << x << " " << y << " " << w << " " << h << "\n"
              << "keys " << keys;
         if ( vec ) { sock << "vvv"; }
         sock << std::endl;
      }

      if (myid == 0)
      {
         connection_failed = !sock && !newly_opened;
      }
      MPI_Bcast(&connection_failed, 1, MPI_INT, 0, comm);
   }
   while (connection_failed);
}


/****************************************************************************************************
* Function: Constructor
* Parameters:
*  h1            - ParFiniteElementSpace to be used for the mesh motion.  The node corresponding to
*                  the cell center is an averaged valued computed in FillCenterVelocitiesWithAvg
*  l2            - ParFiniteElementSpace to be used for the hydrodynamics variables.  We use a DG0
*                  approximation in this code.
*  l2v           - ParFiniteElementSpace for the velocity component.
*  cr            - ParFiniteElementSpace to assist with the mesh velocity reconstruction.  This object
*                  holds the intermediate face velocity object computed in equation (5.7) of the paper.
*  m             - ParLinearForm containing the initial density of each cell to be used to compute the
*                  initial cell mass, which is guaranteed to be conserved locally.
*  _pb           - ProblemBase class containing problem specific functions and quantities.
* Options:
*  use_viscosity - Boolean for the optional addition of viscosity to ensure our method is IDP.
*  mm            - Boolean to allow for mesh motion to guarantee local conservation of mass.
*  CFL           - Double representing our time step restriction.
*
* Purpose: Instantiate LagrangianLOOperator class.
****************************************************************************************************/
LagrangianLOOperator::LagrangianLOOperator(const int &_dim,
                                           const int size,
                                           ParFiniteElementSpace &h1,
                                           ParFiniteElementSpace &h1_l,
                                           ParFiniteElementSpace &l2,
                                           ParFiniteElementSpace &l2v,
                                           ParFiniteElementSpace &cr,
                                           const ParGridFunction &rho0_gf,
                                           ParLinearForm *m,
                                           ProblemBase *_pb,
                                           Array<int> offset,
                                           bool use_viscosity,
                                           int elastic_eos,
                                           bool mm, 
                                           double CFL) :
   dim(_dim),
   TimeDependentOperator(size),
   H1(h1),
   H1_L(h1_l),
   H1Lc(H1_L.GetParMesh(), H1_L.FEColl(), 1),
   L2(l2),
   L2V(l2v),
   CR(cr),
   CRc(CR.GetParMesh(), CR.FEColl(), 1),
   rho0_gf(rho0_gf),
   x_gf(&H1),
   mv_gf(&H1),
   v_CR_gf(&CR),
   v_CR_gf_corrected(&CR), 
   v_CR_gf_fluxes(&CR),
   v_geo_gf(&H1_L),
   cell_bdr_flag_gf(&L2),
   LagrangeMultipliers(&L2),
   pmesh(H1.GetParMesh()),
   m_lf(m),
   pb(_pb),
   block_offsets(offset),
   geom(_dim, offset, H1, L2),
   Vsize_H1(H1.GetVSize()),
   TVSize_H1(H1.TrueVSize()),
   GTVSize_H1(H1.GlobalTrueVSize()),
   NDofs_H1(H1.GetNDofs()),
   NDofs_H1L(H1_L.GetNDofs()),
   Vsize_H1L(H1_L.GetVSize()),
   NVDofs_H1(H1.GetNVDofs()), // Scalar vertex dofs
   Vsize_L2(L2.GetVSize()),
   TVSize_L2(L2.TrueVSize()),
   GTVSize_L2(L2.GlobalTrueVSize()),
   NDofs_L2(L2.GetNDofs()),
   Vsize_L2V(L2V.GetVSize()),
   TVSize_L2V(L2V.TrueVSize()),
   GTVSize_L2V(L2V.GlobalTrueVSize()),
   NDofs_L2V(L2V.GetNDofs()),
   NE(pmesh->GetNE()),
   NBE(pmesh->GetNBE()),
   BdrElementIndexingArray(pmesh->GetNumFaces()),
   BdrVertexIndexingArray(pmesh->GetNV()),
   num_elements(L2.GetNE()),
   num_faces(L2.GetNF()),
   num_vertices(pmesh->GetNV()),
   num_edges(pmesh->GetNEdges()),
   vertex_element(pmesh->GetVertexToElementTable()),
   face_element(pmesh->GetFaceToElementTable()),
   edge_vertex(pmesh->GetEdgeVertexTable()),
   ess_bdr(pmesh->bdr_attributes.Max()),
   // Options
   use_viscosity(use_viscosity),
   use_elasticity(elastic_eos),
   mm(mm),
   CFL(CFL),
   ir(IntRules.Get(pmesh->GetElementBaseGeometry(0), 3 * H1.GetOrder(0) + L2.GetOrder(0) - 1)) //NF//MS
{
   // Transpose face_element to get element_face
   Transpose(*face_element, element_face);
   Transpose(*edge_vertex, vertex_edge);

   cout << "Instantiating hydro op\n";
   cout << "block offsets: ";
   block_offsets.Print(cout);
   // block_offsets[0] = 0;
   // block_offsets[1] = block_offsets[0] + Vsize_H1;
   // block_offsets[2] = block_offsets[1] + Vsize_L2;
   // block_offsets[3] = block_offsets[2] + Vsize_L2V;
   // block_offsets[4] = block_offsets[3] + Vsize_L2;

   // This way all face indices that do not correspond to bdr elements will have val -1.
   // I.e. if global face index k corresponds to an interior face, then BdrElementIndexingArray[k] = -1.
   BdrElementIndexingArray = -1.;
   CreateBdrElementIndexingArray();

   // Similar to BdrElementIndexingArray, set default to -1.
   BdrVertexIndexingArray = -1.;
   CreateBdrVertexIndexingArray();

   // Fill cell boundaries
   cell_bdr_flag_gf = -1.;
   FillCellBdrFlag();

   // Initialize LagrangeMultipliers
   LagrangeMultipliers = 0; // TODO: Play with different value here for initialization???

   // Build mass vector
   m_hpv = m_lf->ParallelAssemble();

   // resize v_CR_gf to correspond to the number of faces
   if (dim == 1)
   {
      v_CR_gf.SetSize(num_faces);
      lambda_max_vec.SetSize(num_faces);
   }
   
   // Initialize values of intermediate face velocities
   v_geo_gf = 0.;
   v_CR_gf = 0.;
   v_CR_gf_corrected = 0.;
   v_CR_gf_fluxes = 0.;
   // assert(v_CR_gf.Size() == dim * num_faces);
   // cout << "v_CR_gf.Size after: " << v_CR_gf.Size() << endl;

   // Initialize Dij sparse
   InitializeDijMatrix();   

   // Initialize arrays of tdofs
   // Boundary conditions: all tests use v.n = 0 on the boundary, and we assume
   // that the boundaries are straight.
   if (ess_bdr.Size() > 0)
   {
      for (int d = 0; d < dim; d++)
      {
         // Attributes 1/2/3 correspond to fixed-x/y/z boundaries,
         // i.e., we must enforce v_x/y/z = 0 for the velocity components.
         ess_bdr = 0; ess_bdr[d] = 1;
         H1_L.GetEssentialTrueDofs(ess_bdr, dofs_list, d);
         ess_tdofs.Append(dofs_list);
      }
      /* Set bdr_vals for dofs that should be 0 */
      ess_tdofs_cart_size = ess_tdofs.Size();
      bdr_vals.SetSize(ess_tdofs_cart_size);
      bdr_vals = 0.;

      if (ess_bdr.Size() > 4)
      {
         MFEM_WARNING("May need to enforce additional BCs.\n");
         if (pb->has_mv_boundary_conditions())
         {
            pb->get_additional_BCs(H1_L, ess_bdr, add_ess_tdofs, add_bdr_vals, &geom);
         }

         ess_tdofs.Append(add_ess_tdofs);
         bdr_vals.Append(add_bdr_vals);
      }
   }

   // Set integration rule for Rannacher-Turek space
   IntegrationRules IntRulesLo(0, Quadrature1D::GaussLobatto);
   RT_ir = IntRules.Get(CR.GetFE(0)->GetGeomType(), RT_ir_order);

   // Reset chronos
   chrono_mm.Clear();
   chrono_mm_lin.Clear();
   chrono_state.Clear();
   chrono_dij.Clear();
   chrono_hiop.Clear();

   /* Initialize elasticity object */
   if (use_elasticity)
   {
      elastic = new Elastic(dim, elastic_eos, H1, L2, rho0_gf, ir);
      SetShearModulus(pb->get_shear_modulus());
   }

   // Print some dimension information
   cout << "Vsize_H1: " << Vsize_H1 << endl;
   cout << "Vsize_H1L: " << Vsize_H1L << endl;
   cout << "TVSize_H1: " << TVSize_H1 << endl;
   cout << "GTVSize_H1: " << GTVSize_H1 << endl;
   cout << "NDofs_H1: " << NDofs_H1 << endl;
   cout << "NVDofs_H1: " << NVDofs_H1 << endl;
   cout << "NDofs_H1L: " << NDofs_H1L << endl;
   cout << "NDofs_L2: " << NDofs_L2 << endl;
   cout << "NDofs_L2V: " << NDofs_L2V << endl;
   cout << "Vsize_L2V: " << Vsize_L2V << endl;

   // More readable format
   cout << "num_elements: " << num_elements << endl;
   cout << "num_faces: " << num_faces << endl;
   cout << "num_vertices: " << num_vertices << endl;
   cout << "num_edges: " << num_edges << endl;
}

LagrangianLOOperator::~LagrangianLOOperator()
{
   delete dij_sparse;
   dij_sparse = nullptr;
   if (use_elasticity)
   {
      delete elastic;
      elastic = nullptr;
   }
}


/**
 * @brief Computes the sigma component for a given element.
 *
 * This function calculates the sigma component for a specific element in the mesh.
 * It retrieves the state vector for the given element, computes the density and 
 * specific internal energy, and evaluates the stress tensor and flux based on 
 * the elasticity model.
 *
 * @tparam dim The dimension of the problem (1D, 2D, or 3D).
 * @param S The input state vector containing the solution variables.
 * @param e The index of the element for which the sigma component is computed.
 * @return The computed sigma component for the specified element.
 *
 * @note This function assumes that the elasticity model is being used and that
 *       the necessary data structures (e.g., elastic object) are properly initialized.
 */
void LagrangianLOOperator::ComputeSigmaDComp(const Vector &S, const int &e, DenseMatrix &sigmaD_e) const
{
   assert(use_elasticity);
   sigmaD_e.SetSize(3);
   sigmaD_e = 0.;

   if (pmesh->GetAttribute(e) != 50)
   {
      return;
   }

   Array<int> idx(dim), idy(dim);
   DenseMatrix flux(dim+2,dim);
   Vector U(dim+2);

   /* Fill sigma_e arrays */
   for (int i = 0; i < dim; i++)
   {
      idx[i] = i+1;
      idy[i] = i;
   }
   GetCellStateVector(S,e,U);
   double rho = 1./U[0], es = 0.;
   elastic->ComputeS(e, rho, sigmaD_e);
   // es = elastic->e_sheer(e);
   // flux = pb->ElasticFlux(sigmaD, es, U, pmesh->GetAttribute(e));
   // flux.GetSubMatrix(idx, idy, sigma_e);
   // sigma_e *= -1.;
}


/**
 * @brief Computes the sigma grid function (sigma_gf) based on the input vector S.
 *
 * This function computes the stress tensor grid function (sigma_gf) using the provided vector S.
 * It assumes that the elasticity is being used and that the size of sigma_gf matches the
 * number of degrees of freedom in the L2 space (NDofs_L2).
 *
 * @tparam dim The dimension of the problem.
 * @param S The input vector containing the state information.
 * @param sigma_gf The output parameter that will hold the computed sigma grid function.
 *
 * @note This function will need to be modified for 2D and 3D.
 */
void LagrangianLOOperator::ComputeSigmaGF(const Vector &S, ParGridFunction &sigma_gf) const
{
   assert(this->use_elasticity);
   assert(sigma_gf.Size() == NDofs_L2);
   sigma_gf = 0.;
   DenseMatrix sigmaD_e(3);

   for (int e = 0; e < NDofs_L2; e++)
   {
      if (pmesh->GetAttribute(e) == 50)
      {
         ComputeSigmaDComp(S,e,sigmaD_e);
         sigma_gf[e] = sigmaD_e.FNorm();
      }
   }
}

/**
 * @brief Computes the Jacobian grid function (f_gf) based on the elasticity tensor.
 *
 * This function computes the Jacobian grid function (f_gf) using the elasticity tensor.
 * It assumes that the elasticity is being used and that the size of f_gf matches the
 * number of degrees of freedom in the L2 space (NDofs_L2).
 *
 * @tparam dim The dimension of the problem.
 * @param f_gf The output parameter that will hold the computed F grid function.
 *
 * @note This function will need to be modified for 2D and 3D.
 */
void LagrangianLOOperator::ComputeFGF(ParGridFunction &f_gf) const
{
   assert(this->use_elasticity);
   assert(f_gf.Size() == NDofs_L2);

   DenseMatrix F(3);

   for (int e = 0; e < NDofs_L2; e++)
   {
      elastic->ComputeAvgF(e, F);
      f_gf[e] = F(0,0);
   }
}

void LagrangianLOOperator::ComputeESheerGF(ParGridFunction &e_sheer_gf) const
{
   assert(this->use_elasticity);
   assert(e_sheer_gf.Size() == NDofs_L2);

   for (int e = 0; e < NDofs_L2; e++)
   {
      e_sheer_gf[e] = elastic->e_sheer(e);
   }
}


/* This Mult method is not mass conservative by itself */
void LagrangianLOOperator::Mult(const Vector &S, Vector &dS_dt) const
{
   // Make sure that the mesh positions correspond to the ones in S. This is
   // needed only because some mfem time integrators don't update the solution
   // vector at every intermediate stage (hence they don't change the mesh).
   UpdateMesh(S);

   dS_dt = 0.0;

   /* Update sv_gf, v_gf, ste_gf */
   SolveHydro(S, dS_dt);

   /* Update dx_dt */
   if (this->compute_mv)
   {
      SolveMeshVelocities(S, dS_dt);
   }
   else
   {
      /* This assumes that before Mult has been called, this->mv_gf has been set externally */
      if (this->mv_gf.Size() != Vsize_H1)
      {
         cout << "NVDofs_H1: " << NVDofs_H1 << ", mv_gf.Size(): " << mv_gf.Size() << endl;
         MFEM_ABORT("Invalid mesh velocity size.");
      }
      ParGridFunction dxdt_gf;
      dxdt_gf.MakeRef(&H1, dS_dt, block_offsets[0]);
      dxdt_gf = this->mv_gf;
   }
}


/**
 * @brief Solves the hydrodynamic equations for the Lagrangian low-order operator.
 *
 * This function computes the time derivative of the state vector for each cell in the mesh.
 * It handles both interior and boundary faces, applying appropriate boundary conditions
 * and flux calculations.
 *
 * @tparam dim The dimension of the problem (1, 2, or 3).
 * @param S The input state vector.
 * @param dS_dt The output time derivative of the state vector.
 */
void LagrangianLOOperator::SolveHydro(const Vector &S, Vector &dS_dt) const
{
   chrono_state.Start();

   Vector rhs(dim+2), c(dim), n_int(dim), U_i(dim+2), U_j(dim+2);
   Array<int> fids, oris; 

   int cj = 0;
   double d;

   mfem::Mesh::FaceInformation FI;
   H1.ExchangeFaceNbrData();

   bool is_boundary_cell = false;

   Vector sum_validation(dim);

    /* Add in e_source for taylor green */
    LinearForm *e_source = nullptr;
    if (pb->get_indicator() == "TaylorGreen")
    {
        // Needed since the Assemble() defaults to PA.
       L2.GetMesh()->DeleteGeometricFactors();
       e_source = new LinearForm(&L2);
       TaylorCoefficient coeff;
       IntegrationRule ir = IntRules.Get(L2.GetFE(0)->GetGeomType(), 2);
       DomainLFIntegrator *d = new DomainLFIntegrator(coeff, &ir);
       e_source->AddDomainIntegrator(d);
       e_source->Assemble();
    }

   for (int ci = 0; ci < NDofs_L2; ci++) // Cell iterator
   {
      is_boundary_cell = false;
      Array<int> cell_bdr_arr;
      sum_validation = 0.;
      GetCellStateVector(S, ci, U_i);

      switch (dim)
      {
         case 1: 
         {
            pmesh->GetElementVertices(ci, fids);
            break;
         }
         case 2:
         {
            pmesh->GetElementEdges(ci, fids, oris);
            break;
         }
         case 3:
         {
            pmesh->GetElementFaces(ci, fids, oris);
            break;
         }
         default:
         {
            MFEM_ABORT("Incorrect dimension provided.\n");
         }
      }
      DenseMatrix F_i;
      //NF//MS
      int ci_attr = pmesh->GetAttribute(ci);
      if (use_elasticity && ci_attr == 50)
      {
         const double _rho = 1./U_i[0];
         double _es = 0.;
         DenseMatrix _sigmaD(3);

         // _sigmaD is 3 x 3
         elastic->ComputeS(ci, _rho, _sigmaD);
         _es = elastic->e_sheer(ci);

         F_i = pb->ElasticFlux(_sigmaD, _es, U_i, ci_attr);
      }
      else
      {
         F_i = pb->flux(U_i, ci_attr);
      }

      rhs = 0.;

      for (int j=0; j < fids.Size(); j++) // Face iterator
      {
         CalcOutwardNormalInt(S, ci, fids[j], n_int);
         c = n_int;
         c /= 2.;
         sum_validation.Add(1., c);

         FI = pmesh->GetFaceInformation(fids[j]);

         if (FI.IsInterior())
         {
            // Get index information/state vector for second cell
            if (ci == FI.element[0].index) { 
               cj = FI.element[1].index; 
            }
            else { 
               assert(ci == FI.element[1].index);
               cj = FI.element[0].index; 
            }
            GetCellStateVector(S, cj, U_j); 

            // flux contribution
            DenseMatrix dm;
            //NF//MS
            int cj_attr = pmesh->GetAttribute(cj);
            if (use_elasticity && cj_attr == 50)
            {
               const double _rho = 1./U_j[0];
               double _es = 0.;
               DenseMatrix _sigmaD(3);

               // _sigmaD is 3 x 3
               elastic->ComputeS(cj, _rho, _sigmaD);
               _es = elastic->e_sheer(cj);

               dm = pb->ElasticFlux(_sigmaD, _es, U_j, cj_attr);
            }
            else
            {
               dm = pb->flux(U_j, cj_attr);
            }

            dm += F_i; 
            Vector y(dim+2);
            dm.Mult(c, y);

            rhs -= y;

            /* viscosity contribution */
            if (use_viscosity)
            {
               d = dij_sparse->Elem(ci, cj); 

               Vector z = U_j;
               z -= U_i;
               rhs.Add(d, z);

               /* REMOVE, checking visocity */
               // Array<int> tmp_dofs(2);
               // tmp_dofs[0] = 1, tmp_dofs[1]=2;
               // Vector tmp_vel(dim);
               // z.GetSubVector(tmp_dofs, tmp_vel);
               // if (tmp_vel[0] < 0.)
               // {
               //    cout << "\tcells ci " << ci << " and cj " << cj << " have a neg vel contribution:\n";
               //    tmp_vel.Print(cout);
               //    cout << "corresponding d: " << d << endl;
               // }
            }
         }
         else
         {
            assert(FI.IsBoundary());
            is_boundary_cell = true;
            int bdr_attr = BdrElementIndexingArray[fids[j]];
            // cout << "boundary attribute for face " << fids[j] << ": " << bdr_attr << endl;
            cell_bdr_arr.Append(bdr_attr);
            Vector y_temp(dim+2), y_temp_bdry(dim+2), U_i_bdry(dim+2);
            F_i.Mult(c, y_temp);
            U_i_bdry = U_i;

            /* Enforce Boundary Conditions */
            if (pb->get_indicator() == "saltzmann")
            {
               // Check bdry flag
               int bdr_attribute = BdrElementIndexingArray[fids[j]]; 

               switch (bdr_attribute)
               {
               case 5:
               // Left wall
               {
                  double _rho = 3.9992502342988532;
                  double _p = 1.3334833281256551;
                  Vector _v(dim);
                  _v = 0.;
                  _v[0] = 1.; //e_x
                  Vector _v_neg = _v, _vp = _v;
                  _v_neg *= -1.;
                  _vp *= _p;

                  DenseMatrix F_i_bdry(dim+2, dim);
                  F_i_bdry.SetRow(0, _v_neg);
                  for (int i = 0; i < dim; i++)
                  {
                     F_i_bdry(i+1, i) = _p;
                  }
                  F_i_bdry.SetRow(dim+1, _vp);

                  F_i_bdry.Mult(c, y_temp_bdry);
                  y_temp += y_temp_bdry;
                  break;
               }
               // case 1:
               // case 2:
               // case 3:
               // {
               //    // Negate velocity
               //    for (int _it = 0; _it < dim; _it++)
               //    {
               //       U_i_bdry[_it + 1] = U_i_bdry[_it + 1] * -1;
               //    }
               //    DenseMatrix F_i_slip = pb->flux(U_i_bdry, pmesh->GetAttribute(i));
               //    F_i_slip.Mult(c, y_temp_bdry);
               //    // y_temp *= 2.;
               //    break;
               // }
               
               default:
               {
                  // cout << "invalid boundary attribute: " << bdr_attribute << endl;
                  // cout << "cell : " << ci << endl;
                  // cout << "face: " << fids[j] << endl;  
                  y_temp *= 2.; 
                  break;
               }
               } // switch (bdr_attr)
            } // if saltzmann
            else if (pb->get_indicator() == "Kidder")
            {
               int bdr_attribute = BdrElementIndexingArray[fids[j]];

               switch (bdr_attribute)
               {
               case 4: // inner radius
               case 5: // outer radius
               {
                  // Need face_x for boundary velocity
                  Array<int> face_dofs;
                  Vector face_x(dim), cell_v(dim), cell_vp(dim);
                  H1.GetFaceDofs(fids[j], face_dofs);
                  int face_dof = face_dofs[2];
                  geom.GetNodePositionFromBV(S,face_dof, face_x);
                  /* 
                  Need boundary state for Kidder
                  Note that t should have been set using update_additional_BCs
                  */
                  pb->GetBoundaryState(face_x, bdr_attribute, U_i_bdry);
                  // MFEM_ABORT("Need to replace 't'\n");

                  /* 
                  * Instead of calling the flux on the ghost state, we want the flux where only boundary pressure is enforced
                  * So instead of this call:
                  *    DenseMatrix F_i_bdry = pb->flux(U_i, pmesh->GetAttribute(ci));
                  * We modify just the pressure in the flux F_i
                  */
                  DenseMatrix F_i_bdry = F_i;
                  F_i_bdry.GetRow(0, cell_v);
                  cell_v *= -1.; // Negate the velocity
                  double _rho = 1. / U_i_bdry[0];
                  double _esi = 0.;
                  if (use_elasticity) { MFEM_WARNING("What is the elastic energy on a ghost cell??\n"); }
                  double _sie = pb->specific_internal_energy(U_i_bdry, _esi);
                  double _press = pb->pressure(_rho, _sie, bdr_attribute);
                  for (int i = 0; i < dim; i++)
                  {
                     F_i_bdry(i+1, i) = _press;
                  }
                  cell_vp = cell_v;
                  cell_vp *= _press;
                  F_i_bdry.SetRow(dim+1, cell_vp);

                  F_i_bdry.Mult(c, y_temp_bdry);
                  y_temp += y_temp_bdry;
                  break;
               }
               default:
               {
                  MFEM_ABORT("Invalid boundary attribute, only full ring currently implemented.\n");
                  y_temp *= 2.; 
                  break;
               }
               }
            }
            // else if (use_elasticity && ci_attr == 50)
            // {
            //    /* Negate sigma */
            //    DenseMatrix F_i_bdry = F_i;

            //    for (int i = 0; i < dim; i++)
            //    {
            //       for (int j = 0; j < dim; j++)
            //       {
            //          F_i_bdry(i+1, j) *= -1.;
            //       }
            //    }

            //    Vector sigmap;
            //    F_i_bdry.GetRow(dim+1, sigmap);
            //    sigmap *= -1.;
            //    F_i_bdry.SetRow(dim+1, sigmap);

            //    F_i_bdry.Mult(c, y_temp_bdry);
            //    y_temp += y_temp_bdry;
            // }
            else
            {
               y_temp *= 2.;
            }
            
            // Add in boundary contribution
            rhs -= y_temp;
         } // End boundary face        

      } // End Face iterator

      /* Add in e_source for taylor green */
      if (pb->get_indicator() == "TaylorGreen")
      {
         // cout << "LF size: " << e_source->Size() << endl;
         rhs[dim+1] += e_source->Elem(ci);
      }

      /* 
      Compute current mass, rather than assume mass is conserved
      In several methods we've attempted, mass has not been conserved
      */
      double k = pmesh->GetElementVolume(ci);
      double _mass = k / U_i[0];
      rhs /= _mass;

      // In either case, update dS_dt
      SetCellStateVector(dS_dt, ci, rhs);
   } // End cell iterator
   chrono_state.Stop();
}


/**
 * @brief Enforces L2 boundary conditions on the state vector.
 *
 * This function modifies the state vector to enforce boundary conditions on the velocity components.
 * It iterates over all boundary elements and applies the appropriate boundary conditions based on
 * the boundary attribute.
 *
 * @tparam dim The dimension of the problem (1, 2, or 3).
 * @param S The state vector to be modified.
 * @param t The current time.
 * @param dt The time step.
 *
 * @note This function should be run after the LagrangianLOOperator::Mult` function.
 */
void LagrangianLOOperator::EnforceL2BC(Vector &S, const double &t, const double &dt)
{
   int el, info;
   Array<int> vel_dofs(dim);
   for (int i = 0; i < dim; i++)
   {
      vel_dofs[i] = i + 1;
   }
   Vector vel(dim), Ui(dim+2);

   // Post processing modify computed values to enforce BCs
   if (pb->has_th_boundary_conditions())
   {
      for (int bdr_ci = 0; bdr_ci < NBE; bdr_ci++)
      {
         /* Get bdr attr and adjacent cell to enforce conditions on */
         int bdr_attr = pmesh->GetBdrAttribute(bdr_ci);
         pmesh->GetBdrElementAdjacentElement(bdr_ci, el, info);

         // cout << "bdr el " << bdr_ci << " has attr " << bdr_attr << " and adj cell " << el << endl;
         GetCellStateVector(S, el, Ui);
         
         Ui.GetSubVector(vel_dofs, vel);

         switch (bdr_attr)
         {
         case 1: // Vx = 0
            vel[0] = 0.;
            break;
         case 2: // Vy = 0
            vel[1] = 0.;
            break;
         case 3: // Vz = 0
            vel[2] = 0.;
            MFEM_ABORT("3D not implemented.\n");
            break;
         case 4: // Radial mesh boundary
            assert(pb->get_indicator() == "SodRadial" || pb->get_indicator() == "Sedov");
            MFEM_WARNING("Radial mesh boundary attributes not implemented.\n");
            break;
         default: // All other cases
         {
            if (pb->get_indicator() == "saltzmann")
            {
               if (bdr_attr == 5)
               {
                  // left
                  double _psi = bdr_vals.Last();
                  vel[0] = 1. * _psi;
               }
               else { MFEM_ABORT("Invalid boundary attribute provided.\n"); }
            }
            else if (pb->get_indicator() == "IsentropicVortex")
            {
               if (bdr_attr == 5) { EnforceExactBCOnCell(S, el, t, dt, Ui); }
               else { MFEM_ABORT("Invalid boundary attribute provided.\n"); }
            }
            break;
         } // end default
         } // end switch(bdr_attr)

         Ui.SetSubVector(vel_dofs, vel);
         SetCellStateVector(S, el, Ui);
      }
   } // End pb->has_th_boundary_conditions()
}


/**
 * @brief Updates the mesh velocity boundary conditions.
 *
 * This function updates the boundary conditions for the mesh velocity based on the current time and time step.
 * It checks if the boundary conditions need updating and, if so, resets the additional boundary conditions
 * from the previous iteration and updates them for the current time step.
 *
 * @tparam dim The dimension of the problem (1, 2, or 3).
 * @param t The current time.
 * @param dt The time step.
 *
 * @note This function should be run before the `ODESolver::Step` function.
 */
void LagrangianLOOperator::UpdateMeshVelocityBCs(const double &t, const double &dt)
{
   if (pb->get_mv_bcs_need_updating())
   {
      if (timestep_first == 0.)
      {
         timestep_first = GetTimestep();
      }

      /* Reset additional boundary conditions from previous iteration */
      bdr_vals.SetSize(ess_tdofs_cart_size);
      pb->update_additional_BCs(t, timestep_first, add_bdr_vals, &geom, &x_gf);
      bdr_vals.Append(add_bdr_vals);
   }
}



/**
 * @brief Solves for the mesh velocities.
 *
 * This function computes the mesh velocities based on the current state vector. It handles both 1D and higher dimensions,
 * applying different methods for computing the mesh velocities depending on the specified options.
 *
 * @tparam dim The dimension of the problem (1, 2, or 3).
 * @param S The input state vector.
 * @param dS_dt The output time derivative of the state vector.
 */
void LagrangianLOOperator::SolveMeshVelocities(const Vector &S, Vector &dS_dt) const
{
   // cout << "========================================\n"
   //      << "          SolveMeshVelocities           \n"
   //      << "========================================\n";
   chrono_mm.Start();
   ParGridFunction dxdt_gf, dxdt_gf_l(&H1_L);
   dxdt_gf.MakeRef(&H1, dS_dt, block_offsets[0]);

   ComputeIntermediateFaceVelocities(S);

   if (dim == 1)
   {
      /* Use intermediate face velocities to move the mesh */
      for (int face = 0; face < num_faces; face++)
      {
         Vector _vel(dim);
         GetIntermediateFaceVelocity(face, _vel);
         geom.UpdateNodeVelocity(dxdt_gf, face, _vel);
      }
      // Since we do not project here, must fill cell center velocities
      FillCenterVelocitiesWithAvg(dS_dt);
   }
   else
   {
      switch (mv_option)
      {
      /* 
      0-9   - no lagrange multipliers 
      10-19 - Dense LM
      20-29 - Sparse LM
      30    - Sparse Viscous LM
      */

      case 0: // arithmetic average of adjacent cells with distributed viscosity
      {
         MFEM_ABORT("Current mesh velocity calculation broken\n");
         // bool is_weighted = false;
         // int td_flag = 0;
         // CalcCellAveragedCornerVelocityVector(S, S_old, is_weighted, td_flag, mv_gf_l);
         break;
      }
      case 1: // arithmetic average of adjacent cells, no viscosity
      {
         MFEM_ABORT("Current mesh velocity calculation broken\n");
         // bool is_weighted = false;
         // int td_flag = 0;
         // CalcCellAveragedCornerVelocityVector(S, S_old, is_weighted, td_flag, mv_gf_l);
         // DistributeFaceViscosityToVelocity(S_old, mv_gf_l);
         break;
      }

      case 2:
      {
         ComputeGeoVNormal(S, dxdt_gf_l);
         break;
      }
      case 3: // weighted average of adjacent cell midpoint in time
      {
         MFEM_ABORT("Current mesh velocity calculation broken\n");
         /* Petrov-Galerkin Justified weighted average */
         // bool is_weighted = true;
         // int td_flag = 2;
         // CalcCellAveragedCornerVelocityVector(S, S_old, is_weighted, td_flag, mv_gf_l);
         break;
      }
      case 4: // weighted average of adjacent cell midpoint in time with distributed viscosity
      {
         MFEM_ABORT("Current mesh velocity calculation broken\n");
         /* Petrov-Galerkin Justified weighted average */
         // bool is_weighted = true;
         // int td_flag = 2;
         // CalcCellAveragedCornerVelocityVector(S, S_old, is_weighted, td_flag, mv_gf_l);
         // DistributeFaceViscosityToVelocity(S_old, mv_gf_l);
         break;
      }
      
      // case 3:
      //    ComputeGeoVCellFaceNormal(S);
      //    break;

      // case 1:
      //    ComputeGeoVRaviart(S);
      //    break;

      // case 4:
      //    ComputeGeoVCAVEAT(S);
      //    break;
      
      // case 5:
      //    ComputeGeoVCAVEATCellFace(S);
      //    break;
      
      // case 6:
      //    ComputeGeoVCAVEATCellFaceWeighted(S);
      //    break;

      // case 8: // HiOp Dense Solver for optimization problem
      // {
      //    // MFEM_ABORT("Hiop Dense not implemented.\n");
      //    SolveHiOpDense(S, S_old, t, dt, mv_gf_l);
      //    break;
      // }
      // case 9: // HiOp Sparse Solver for optimization problem
      // {
      //    SolveHiOp(S, S_old, t, dt, mv_gf_l);
      //    break;
      // }
      // case 10: // mv2 but with viscosity distribution
      // {
      //    ComputeGeoVNormalDistributedViscosity(S);
      //    // MFEM_ABORT("Viscosity distribution mesh movement not yet implemented.\n");
      //    break;
      // }
      // case 11: // arithmetic average of adjacent cells with distributed viscosity
      // {
      //    bool is_weighted = false;
      //    CalcCellAveragedCornerVelocityVector(S_old, is_weighted, mv_gf_l);
      //    DistributeFaceViscosityToVelocity(S_old, mv_gf_l);

      //    break;
      // }
      // case 12: // arithmetic average of adjacent cells, no viscosity
      // {
      //    bool is_weighted = false;
      //    CalcCellAveragedCornerVelocityVector(S_old, is_weighted, mv_gf_l);
      //    break;
      // }

      case 10:
      case 11: 
      case 12:
      case 13:
      case 14:
      case 15:
      case 16:
      case 17:
      case 18:
      case 19:
      {
         MFEM_ABORT("Current mesh velocity calculation broken\n");
         // int target_option = mv_option % 10;
         // int lm_option = 1; // Indicates Target
         // SolveHiOp(S, S_old, lm_option, target_option, t, dt, mv_gf_l);
         break;
      }
      case 20:
      case 21: 
      case 22:
      case 23:
      case 24:
      case 25:
      case 26:
      case 27:
      case 28:
      case 29:
      {
         MFEM_ABORT("Current mesh velocity calculation broken\n");
         // int target_option = mv_option % 10;
         // SolveHiOpDense(S, S_old, target_option, t, dt, mv_gf_l);
         break;
      }
      case 30:
      {
         MFEM_ABORT("Current mesh velocity calculation broken\n");
         // int target_option = -1;
         // int lm_option = 2; // Indicates Viscous
         // SolveHiOp(S, S_old, lm_option, target_option, t, dt, mv_gf_l);
         break;
      }

      default:
         MFEM_ABORT("Invalid mesh velocity option.\n");
      }

      /* mesh velocity only needs to be linearized and boundary condition enforced 
      only when not solving for them via a constrained optimization problem*/
      if (mv_option < 10)
      {
         /* Optionally, Linearize velocity */
         if (this->do_mv_linearization)
         {
            MFEM_ABORT("Mesh velocity linearization is not currently available\n");
            // ParGridFunction dxdt_gf_l_linearized(&H1_L);
            // ComputeLinearizedNodeVelocities(dxdt_gf_l, dxdt_gf_l_linearized, t, dt);
            // mv_gf_l = mv_gf_l_linearized;
         }

         /* Optionally, enforce boundary conditions */
         if (pb->has_mv_boundary_conditions())
         {
            for (int i = 0; i < ess_tdofs.Size(); i++) { 
               // cout << "ess_tdof: " << ess_tdofs[i] << ", bdr val: " << bdr_vals[i] <<endl;
               dxdt_gf_l(ess_tdofs[i]) = bdr_vals[i]; 
            }

         }
      }
      
      /* 
      Project back onto mv_gf 
      No need to fill center velocities as this is handled by the projection
      */
      dxdt_gf.ProjectGridFunction(dxdt_gf_l);

      if (NDofs_H1 != NDofs_H1L)
      {
         /* Must compute mesh velocities on cell centers and faces */
         /* Choose behavior for face velocity */
         switch (fv_option)
         {
         case 1:
            MFEM_ABORT("Current face mesh velocity calculation broken\n");
            // ComputeCorrectiveFaceVelocities(S, t, dt);
            break;
         
         case 2:
            FillFaceVelocitiesWithAvg(dxdt_gf);
            break;

         case 3:
            FillFaceVelocitiesWithButterfly(dxdt_gf);
            break;
         
         default:
            /* Do nothing */
            break;
         } // End face velocity switch case

      }
   
   } // End dim > 1

   dxdt_gf.SyncAliasMemory(dS_dt);
   chrono_mm.Stop();
}

/****************************************************************************************************
* Function: InitializeDijMatrix
*
* Purpose:
*  The purpose of this function is to initialize a SparseMatrix with the sparsity pattern
*  that corresponds to the adjacency matrix of the mesh.  I.E. this functions constructs a matrix
*  of size num_elements x num_elements where the nonzero elements are for each element row
*  are the cell indices in the stencil. This SparseMatrix will encode the max wave speed of
*  each corresponding local Riemann problem.  An example of cells in the stencil of cell c is given
*  below.
*
*                  -----------
*                 | 0 | 1 | 0 |
*                 |-----------|
*                 | 1 | c | 1 |
*                 |-----------|
*                 | 0 | 1 | 0 |
*                  -----------
* 
* This is accomplished in the following steps:
*  1) Create ParBilinearForm based on L2FESpace
*  2) Create DGTraceIntegrator(DG state coeff, CG coeff (ones), double, double)
*  3) AddInteriorFaceIntegrator using ^
*  4) Assemble, Finalize, ParallelAssemble calls
*  5) HypreParMatrix::MergeDiagAndOffd --> SparseMatrix
****************************************************************************************************/
void LagrangianLOOperator::InitializeDijMatrix()
{
   switch (dim)
   {
      case 1: 
      {
         el_num_faces = pmesh->GetElement(0)->GetNVertices();
         break;
      }
      case 2: 
      {
         el_num_faces = pmesh->GetElement(0)->GetNEdges();
         break;
      }
      case 3: 
      {
         el_num_faces = pmesh->GetElement(0)->GetNFaces();
         break;
      }
      default: MFEM_ABORT("Wrong dimension given.\n");
   }

   // Initialize SparseMatrix
   dij_sparse = new SparseMatrix(num_elements, num_elements, el_num_faces+1);

   // Create dummy coefficients
   using namespace std::placeholders;
   std::function<double(const Vector &,const double)> rho0_static = 
      std::bind(&ProblemBase::rho0, pb, std::placeholders::_1, std::placeholders::_2);

   FunctionCoefficient rho_coeff(rho0_static);

   // Assemble SparseMatrix object
   ParBilinearForm k(&L2);
   k.AddInteriorFaceIntegrator(new DGDiffusionIntegrator(rho_coeff, 1.0, 1.0));

   k.Assemble();
   k.Finalize();

   HypreParMatrix * k_hpm = k.ParallelAssemble();
   k_hpm->MergeDiagAndOffd(*dij_sparse);

   if (k_hpm) { delete k_hpm; } // Clear memory to avoid leak
   // From here, we can modify the sparse matrix according to the sparsity pattern
}


/****************************************************************************************************
* Function: BuildDijMatrix
* Parameters:
*  S - BlockVector that stores mesh information, mesh velocity, and state variables.
*
* Purpose:
****************************************************************************************************/
void LagrangianLOOperator::BuildDijMatrix(const Vector &S)
{
   // cout << "=======================================\n"
   //      << "           Build Dij Matrix            \n"
   //      << "=======================================\n";
   mfem::Mesh::FaceInformation FI;
   int c, cp;
   Vector Uc(dim+2), Ucp(dim+2), n_int(dim), c_vec(dim), n_vec(dim);
   double F, lambda_max, d, dij_avg = 0.;

   for (int face = 0; face < num_faces; face++) // face iterator
   {
      // cout << "face: " << face << endl;
      FI = pmesh->GetFaceInformation(face);
      c = FI.element[0].index;
      cp = FI.element[1].index;

      GetCellStateVector(S, c, Uc);
      if (FI.IsInterior())
      {
         GetCellStateVector(S, cp, Ucp);

         // Get normal, d, and |F|
         CalcOutwardNormalInt(S, c, face, n_int);
         n_vec = n_int;
         double F = n_vec.Norml2();
         n_vec /= F;

         if (1. - n_vec.Norml2() > 1e-12)
         {
            cout << "n_vec: ";
            n_vec.Print(cout);
            MFEM_ABORT("Invalid normal vector.\n");
         }
         c_vec = n_int;
         c_vec /= 2.;
         double c_norm = c_vec.Norml2();

         /* Compute max wave speed */ 
         // Some cases require an update to b_covolume at every interface.  This can be done through
         // the function ProblemBase::lm_update(), which is overridden only in the functions it is used.
         // double b_covolume = .1 / (max(1./Uc[0], 1./Ucp[0]));
         // pb->lm_update(b_covolume);

         // Compute sheer energy, if applicable
         double esl = 0., esr = 0.;
         if (use_elasticity)
         {
            esl = elastic->e_sheer(c);
            esr = elastic->e_sheer(cp);
         }

         // Compute pressure with given EOS
         double rhoL = 1. / Uc[0];
         double sieL = pb->specific_internal_energy(Uc, esl);
         double pl = pb->pressure(rhoL, sieL, pmesh->GetAttribute(c));

         double rhoR = 1. / Ucp[0];
         double sieR = pb->specific_internal_energy(Ucp, esr);
         double pr = pb->pressure(rhoR, sieR, pmesh->GetAttribute(cp));

         // Finally compute lambda max
         if (use_elasticity)
         {
            lambda_max = pb->compute_lambda_max(Uc, Ucp, n_vec, esl, esr, pl, pr, this->use_greedy_viscosity);
            lambda_max = sqrt(pow(lambda_max, 2) + std::max(rhoL,rhoR) * 4./3. * pb->get_shear_modulus());
         }
         else
         {
            lambda_max = pb->compute_lambda_max(Uc, Ucp, n_vec, esl, esr, pl, pr, this->use_greedy_viscosity);
         }
         
         d = lambda_max * c_norm;
         dij_avg += d;

         dij_sparse->Elem(c,cp) = d;
         dij_sparse->Elem(cp,c) = d;

         if (dim == 1)
         {
            lambda_max_vec[face] = lambda_max; // TODO: remove, only temporary
         }
      }
   } // End face iterator
   dij_avg = dij_avg / num_faces;
   ts_dijavg.Append(dij_avg);
}


/****************************************************************************************************
* Function: CalculateTimestep
* Parameters:
*  S - BlockVector that stores mesh information, mesh velocity, and state variables.
*
* Purpose:
*  Computes the maximum time step that can be taken under the CFL condition given in section 4.4:
*           dt * Sum_{c' \in J^*(c)} \frac{d_cc'}{mcrho} \leq \frac{CFL}{2} \forall c.
*
*  This is accomplished by iterative over each cell and computing the corresponding max dt, then 
*  a minimum of the current dt, and the computed restriction.
****************************************************************************************************/
void LagrangianLOOperator::CalculateTimestep(const Vector &S)
{
   // cout << "CalculateTimestep\n";
   double t_min = 1.;
   double t_temp = 0;
   double mi = 0;
   int cj = 0;

   Array<int> fids, oris;
   Vector U_i(dim+2), U_j(dim+2);
   double d=0., temp_sum = 0.;

   mfem::Mesh::FaceInformation FI;

   for (int ci = 0; ci < NDofs_L2; ci++) // Cell iterator
   {
      // cout << "\tcell: " << ci << endl;
      temp_sum = 0.;

      GetCellStateVector(S, ci, U_i);

      // Compute mass at time tn
      const double k = pmesh->GetElementVolume(ci);
      mi = k / U_i[0];

      if (mi <= 0.)
      {
         cout <<  "Invalid mass at cell " << ci << ": " << mi << endl;
         MFEM_ABORT("Invalid mass.\n");
      }

      H1.ExchangeFaceNbrData();

      switch (dim)
      {
         case 1: 
         {
            pmesh->GetElementVertices(ci, fids);
            break;
         }
         case 2:
         {
            pmesh->GetElementEdges(ci, fids, oris);
            break;
         }
         case 3:
         {
            pmesh->GetElementFaces(ci, fids, oris);
         }
      }

      for (int j=0; j < fids.Size(); j++) // Face iterator
      {
         FI = pmesh->GetFaceInformation(fids[j]);

         if (FI.IsInterior())
         {
            // Get index information/state vector for second cell
            if (ci == FI.element[0].index) { 
               cj = FI.element[1].index; 
            }
            else { 
               cj = FI.element[0].index; 
            }

            d = dij_sparse->Elem(ci, cj); 

            // cout << "face: " << fids[j] << endl;
            // cout << "d for ci " << ci << " and cj " << cj << ": " << d << endl;

            temp_sum += d;
         }
      }
      
      t_temp = 0.5 * ((CFL * mi) / temp_sum );

      if (t_temp < t_min && t_temp > 1.e-12) { 
         // cout << "timestep reduced\n";
         t_min = t_temp;
      }
   } // End cell iterator

   this->timestep = t_min;
}


/****************************************************************************************************
* Function: GetEntityDof
* Parameters:
*     GDof   - Global degree of freedom with a max of H1.GetNDofs()
*     Entity - enum representing what type of node GDof corresponds to
*     EDof   - Entity index of the GDof
*
* Purpose:
*     This function is a helper function to identify which entity a node corresponds to
*     and to convert from the global numbering to the entity's numbering.
****************************************************************************************************/
void LagrangianLOOperator::GetEntityDof(const int GDof, DofEntity & entity, int & EDof)
{
   // cout << "================== \n GetEntityDof \n================== \n";
   
   switch (dim)
   {
      case 1:
      {
         if (GDof < num_vertices)
         {
            // In 1d a face is a vertex, however our method to compute RT velocities
            // should return the intermediate face velocity.  Thus we shall treat
            // them as faces.
            entity = DofEntity::face;
            EDof = GDof;
         }
         else if (GDof < Vsize_H1) 
         {
            entity = DofEntity::cell;
            EDof = GDof - num_vertices;
         }
         break;
      }
      case 2:
      {
         if (GDof < num_vertices)
         {
            entity = DofEntity::corner;
            EDof = GDof;
         }
         else if (GDof < num_vertices + num_faces)
         {
            entity = DofEntity::face;
            EDof = GDof - num_vertices;
         }
         // TODO: Add edge implementation
         else
         {
            entity = DofEntity::cell;
            EDof = GDof - num_vertices - num_faces;
         }
         break;
      }
      case 3:
      {
         MFEM_ABORT("3D not implemented yet.");
         break;
      }
      default:
      {
         MFEM_ABORT("Incorrect dim provided.");
      }
   }
}


/****************************************************************************************************
* Function: CreateBdrElementIndexingArray
*  
* Purpose: 
*  To fill the Vector BdrElementIndexingArray of size NumFaces. If a face is a boundary face, then
*  BdrElementIndexingArray[face] = bdr_attr, and if the face is an interior face, then we will have
*  BdrElementIndexingArray[face] = -1.
****************************************************************************************************/
void LagrangianLOOperator::CreateBdrElementIndexingArray()
{
   cout << "Constructing BdrElementIndexingArray:\n";
   for (int i = 0; i < pmesh->GetNBE(); i++)
   {
      int bdr_attr = pmesh->GetBdrAttribute(i);
      int index = pmesh->GetBdrElementFaceIndex(i);
      BdrElementIndexingArray[index] = bdr_attr;
   }
}


/****************************************************************************************************
* Function: CreateBdrVertexIndexingArray
*
* Purpose: 
*  To fill the Vector BdrVertexIndexingArray of size Numvertices. If a vertex is a boundary vertex, 
*  then BdrVertexIndexingArray[vertex] = 1, and if the vertex is an interior vertex, then we will 
*  have BdrVertexIndexingArray[vertex] = -1. This is done by iterating over boundary elements, 
*  grab edges, and fill in the corresponding boundary attribute for the adjacent vertices
*
*  Note that this function will label a vertex with 5 if that vertex lies on the corner, indicating 
*  that the vertex velocity should be set to 0 to preserve the slip BCs on both of its faces.
****************************************************************************************************/
void LagrangianLOOperator::CreateBdrVertexIndexingArray()
{
   // 3DTODO: Will need to modify this for faces instead of edges
   Array<int> fids, oris;
   Array<int> verts;

   for (int i = 0; i < pmesh->GetNBE(); i++)
   {
      int bdr_attr = pmesh->GetBdrAttribute(i);
      int face = pmesh->GetBdrElementFaceIndex(i);

      if (dim == 1)
      {
         BdrVertexIndexingArray[face] = bdr_attr;
      }
      else
      {
         assert(dim == 2);

         pmesh->GetEdgeVertices(face, verts);
         for (int k = 0; k < verts.Size(); k++)
         {
            int index = verts[k];
            // TODO: Get rid of all get_indicator() funcalls, replace with use of class int problem.
            if (pb->get_indicator() == "saltzmann")
            {
               // Replace the bdr attribute in the array as long as it is not
               // the dirichlet condition (For Saltzmann Problem)
               // This ensures the left wall vertices have the proper indicator
               if (BdrVertexIndexingArray[index] != 5)
               {
                  BdrVertexIndexingArray[index] = bdr_attr;
               }
            }
            else if (pb->get_indicator() == "Sod" ||
                     pb->get_indicator() == "TriplePoint" || 
                     pb->get_indicator() == "riemann" ||
                     pb->get_indicator() == "Sedov" || 
                     pb->get_indicator() == "SodRadial" || 
                     pb->get_indicator() == "TaylorGreen")
            {
               // Mark corner vertices as 5
               // These nodes should not move at all during the simulation
               // Identify these corner vertices as those that already have
               // a value non negative
               if (BdrVertexIndexingArray[index] != -1 && BdrVertexIndexingArray[index] != bdr_attr)
               {
                  BdrVertexIndexingArray[index] = 5;
               }
               else 
               {
                  BdrVertexIndexingArray[index] = bdr_attr;
               }
            }
            else
            {
               BdrVertexIndexingArray[index] = bdr_attr;
            }
         } // end vertex iterator
      }

   } // end boundary elements
}


/****************************************************************************************************
* Function: FillCellBdrFlag
*
* Purpose: 
*  To fill an L2 gridfunction of size NDofs_L2 (num cells) with a value indicating if the cell
*  is on the boundary or not.  This ParGridFunction can be retrieved with GetCellBdrFlagGF().
*
*  Note that this function will label a cell with 5 if that cell lies on the corner, indicating 
*  that the cell velocity should be set to 0 to preserve the slip BCs on both of its faces.
****************************************************************************************************/
void LagrangianLOOperator::FillCellBdrFlag()
{
   for (int i = 0; i < pmesh->GetNBE(); i++)
   {
      int bdr_attr = pmesh->GetBdrAttribute(i);
      int face = pmesh->GetBdrElementFaceIndex(i);

      Array<int> row;
      // Get cell
      face_element->GetRow(face, row);
      
      // Each face should only have 1 adjacent cell
      assert(row.Size() == 1);

      if (cell_bdr_flag_gf[row[0]] != -1 && 
          (pb->get_indicator() == "Sod" || 
           pb->get_indicator() == "TriplePoint" || 
           pb->get_indicator() == "riemann" ||
           pb->get_indicator() == "Sedov" || 
           pb->get_indicator() == "SodRadial"))
      {
         // Corner cell
         cell_bdr_flag_gf[row[0]] = 5;
      }
      else
      {
         cell_bdr_flag_gf[row[0]] = bdr_attr;
      }
      
      // cout << "cell: " << row[0] << ", bdr_attr: " << bdr_attr << endl;
   }
}


void LagrangianLOOperator::ComputeKidderAvgIntExtRadii(const Vector &S, double &avg_rad_int, double &avg_rad_ext)
{
   // cout << "ComputeKidderAvgIntExtRadii\n";
   Vector face_x(dim);
   Array<int> face_dofs_row;
   int num_int_faces = 0, num_ext_faces = 0;
   avg_rad_int = 0., avg_rad_ext = 0.;

   /* Iterate over boundary faces */
   for (int i = 0; i < pmesh->GetNBE(); i++)
   {
      // Get Face coords
      int bdr_attr = pmesh->GetBdrAttribute(i);
      int face = pmesh->GetBdrElementFaceIndex(i);
      H1.GetFaceDofs(face, face_dofs_row);
      int face_dof = face_dofs_row[2];
      geom.GetNodePositionFromBV(S, face_dof, face_x);

      if (bdr_attr == 4)
      {
         // cout << "interior face\n";
         double val = face_x.Norml2();
         avg_rad_int += val;
         num_int_faces++;
      }
      else if (bdr_attr == 5)
      {
         // cout << "exterior face\n";
         double val = face_x.Norml2();
         avg_rad_ext += val;
         num_ext_faces++;
      }
      else
      {
         // cout << "other boundary face\n";
         MFEM_ABORT("Only full ring is implemented in Kidder problem.\n");
      }
   }
   
   /* By definition */
   assert(num_int_faces == num_ext_faces);

   /* Compute average */
   avg_rad_int /= num_int_faces;
   avg_rad_ext /= num_ext_faces;
}


void LagrangianLOOperator::ComputeKidderAvgDensityAndEntropy(const Vector &S, double &avg_density, double &avg_entropy)
{
   ParGridFunction sv_gf;
   Vector* sptr = const_cast<Vector*>(&S);
   sv_gf.MakeRef(&L2, *sptr, block_offsets[1]);
   Vector U(dim+2);
   avg_density = 0., avg_entropy = 0.;

   for (int cell_it = 0; cell_it < NDofs_L2; cell_it++)
   {
      GetCellStateVector(S, cell_it, U);
      double density = 1. / U[0]; 
      double e_sheer = 0.;
      if (use_elasticity) { e_sheer = elastic->e_sheer(cell_it); }
      double sie = pb->specific_internal_energy(U, e_sheer);
      double pressure = pb->pressure(density, sie, pmesh->GetAttribute(cell_it));
      double entropy = pressure / pow(density, pb->get_gamma());

      avg_density += density;
      avg_entropy += entropy;
   }

   avg_density /= NDofs_L2;
   avg_entropy /= NDofs_L2;
}


/****************************************************************************************************
* Function: ComputeMinDetJ
* Parameters:
*  minDetJ - Minimum value for the determinant of the cell transformation over the mesh
*  cell    - coresponding cell
*
* Purpose:
*     Compute the minimum determinant of the jacobian of a cell transformation for all cells in the mesh.
*     This value can be used to verify if the mesh has collapsed. If minDetJ < 0., then the mesh has collapsed.
*     The Jacobian is checked at the following quadrature points:
*        (.211, .211)
*        (.789, .211)
*        (.789, .789)
*        (.211, .789)
*     If the determinant of the Jacobian at any of the above points
*     is negative, this means that our element transformation has 
*     reversed orientation, or twisted, and thus the mesh has 
*     collapsed.
****************************************************************************************************/
void LagrangianLOOperator::ComputeMinDetJ(int &cell, double &minDetJ)
{
   ElementTransformation * trans;
   minDetJ = 100.;

   for (int el = 0; el < L2.GetNE(); el++)
   {
      trans = pmesh->GetElementTransformation(el);
      /* Iterate over 4 integration points */
      for (int i = 0; i < RT_ir.GetNPoints(); i++)
      {
         const IntegrationPoint &ip = RT_ir.IntPoint(i);
         trans->SetIntPoint(&ip);
         double detJ = trans->Weight();
         if (detJ < minDetJ)
         {
            minDetJ = detJ;
            cell = el;
         }
      }
   }
}


/****************************************************************************************************
* Function: ComputeTimeSeriesData
* Parameters:
*  t   - Current time
*  dt  - Current timestep
*
* Purpose:
*  This function computes various quantities that are saved at each time iteration
*        ts_dijmax
*        td_t
*        td_timestep
*        ts_min_detJ
*        ts_min_detJ_cell
*        ts_kidder_avg_rad_int
*        ts_kidder_avg_rad_ext
*        ts_kidder_avg_rad_density
*        ts_kidder_avg_rad_entropy
*
*  In addition, if ts_min_detJ which is calculated in the routine ComputeMinDetJ returns negative,
*  this function returns true to indicate that the mesh has collapsed.
****************************************************************************************************/
bool LagrangianLOOperator::ComputeTimeSeriesData(const Vector & S, const double &t, const double &dt)
{
   ts_dijmax.Append(dij_sparse->MaxNorm());
   ts_t.Append(t);
   ts_timestep.Append(dt);
   int minDetJCell;
   double minDetJ;
   ComputeMinDetJ(minDetJCell, minDetJ); 
   ts_min_detJ.Append(minDetJ);
   ts_min_detJ_cell.Append(minDetJCell);
   if (pb->get_indicator() == "Kidder")
   {
      // Compute average exterior and interior radius
      double avg_rad_int, avg_rad_ext;
      ComputeKidderAvgIntExtRadii(S, avg_rad_int, avg_rad_ext);
      ts_kidder_avg_rad_int.Append(avg_rad_int);
      ts_kidder_avg_rad_ext.Append(avg_rad_ext);
      double avg_density, avg_entropy;
      ComputeKidderAvgDensityAndEntropy(S, avg_density, avg_entropy);
      ts_kidder_avg_density.Append(avg_density);
      ts_kidder_avg_entropy.Append(avg_entropy);
   }

   /*
   If mesh has tangled, modify isCollapsed parameter to stop computation and print final data
   */
   if (minDetJ < 0.)
   {
      return true;
   }

   return false;
}

/****************************************************************************************************
* Function: IsMeshCollapsedGeom
*
* Purpose:
*     Verify if the mesh has collapsed. If collapsed, return true.
*     This method checks each cell to see if any corner node has 
*     passed the opposite face. This is accomplished by iterating
*     over the faces of a cell and checking the two nodes that are
*     not adjacent to the face. The vector 1F and 2F are dotted 
*     with the outward normal at F. If either of [these quantities 
*     is positive, we know the mesh has inverted.
*
*           1----------2
*           |          |
*           |          |
*           |          |
*           |          |
*           3----F-----4
*
****************************************************************************************************/
bool LagrangianLOOperator::IsMeshCollapsedGeom(const Vector &S)
{
   mfem::Mesh::FaceInformation FI;
   Array<int> fids, oris, verts, face_dofs;
   Vector n_vec(dim), temp_vec(dim), face_x(dim), vert_x(dim);
   int face_dof;

   for (int ci = 0; ci < NDofs_L2; ci++)
   {
      if (dim != 2)
      {
         MFEM_ABORT("Function IsMeshCollapsed is hardcoded for dim == 2.")
      }

      pmesh->GetElementVertices(ci, verts);
      pmesh->GetElementEdges(ci, fids, oris);

      // cout << "cell: " << ci << endl;

      for (int j=0; j < fids.Size(); j++) // Face iterator
      {
         FI = pmesh->GetFaceInformation(fids[j]);
         H1.GetFaceDofs(fids[j], face_dofs);
         // cout << "\tface: " << fids[j];
         face_dof = face_dofs[2];

         // Get outward normal
         CalcOutwardNormalInt(S, ci, fids[j], n_vec);

         // Get nodal loc of face dof
         geom.GetNodePositionFromBV(S, face_dof, face_x);

         // Iterate over corner vertices of cell
         for (int k=0; k < verts.Size(); k++)
         {
            // We only care about verts opposite to the face
            if (face_dofs.Find(verts[k]) == -1)
            {
               // cout << ", opp node: " << verts[k];
               geom.GetNodePositionFromBV(S, verts[k], vert_x);
               subtract(vert_x, face_x, temp_vec);
               double indicator = temp_vec * n_vec;
               if (indicator > 0.)
               {
                  cout << "mesh should collapse\n";
                  return true;
               }
            }
         }
         // cout << endl;
      }

   }
   // MFEM_ABORT("End check.");
   return false;
}


/****************************************************************************************************
* Function: EnforceExactBCOnCell
* Parameters:
*  cell      - Index representing the cell
*  state_val - Returned vector containing the exacthydrodynamic state variables
*
* Purpose:
*  Given a cell, return the exact hydrodynamic state variables of that cell.
*  This function is used to enforce the exact solution on the boundary at a
*  particular cell.
****************************************************************************************************/
void LagrangianLOOperator::EnforceExactBCOnCell(const Vector &S, const int & cell, const double &t, 
                                                     const double &dt, Vector & state_val)
{
   // cout << "========================================\n"
   //      << "         EnforceExactBCOnCell           \n"
   //      << "========================================\n";

   int cell_vdof;
   Vector cell_x(dim);

   switch (dim)
   {
      case 1:
      {
         cell_vdof = num_faces + cell;
         break;
      }
      case 2:
      {
         cell_vdof = NVDofs_H1 + num_faces + cell;
         break;
      }
      case 3:
      {
         MFEM_ABORT("3D not implemented.\n");
      }
      default:
      {
         MFEM_ABORT("Incorrect dim value provided.\n");
      }
   }

   geom.GetNodePositionFromBV(S, cell_vdof, cell_x);

   // Fill in val with exact solution
   state_val[0] = pb->sv0(cell_x, t+dt);
   Vector v_exact(dim);
   pb->v0(cell_x, t+dt, v_exact);
   for (int i = 0; i < dim; i++)
   {
      state_val[1 + i] = v_exact[i];
   }
   state_val[dim + 1] = pb->ste0(cell_x, t+dt);
}


/****************************************************************************************************
* Function: SetMassConservativeDensity
* Parameters:
*  S     - BlockVector that stores mesh information, mesh velocity, and state variables.
*
* Purpose:
*  Postprocess the density to be exactly mass conservative.  This function must be called
*  after the mesh motion has been calculated.
****************************************************************************************************/
void LagrangianLOOperator::SetMassConservativeDensity(Vector &S, double &pct_corrected, double &rel_mass_corrected)
{
   // cout << "========================================\n"
   //      << "       SetMassConservativeDensity       \n"
   //      << "========================================\n";
   UpdateMesh(S);
   int num_corrected_cells = 0;
   ParGridFunction sv_gf;
   sv_gf.MakeRef(&L2, S, block_offsets[1]);
   rel_mass_corrected = 0.;

   for (int cell_it = 0; cell_it < NDofs_L2; cell_it++)
   {
      // Get cell mass
      const double m = m_hpv->Elem(cell_it);

      // Get new cell volume
      const double k_new = pmesh->GetElementVolume(cell_it);
      const double sv_new= sv_gf.Elem(cell_it);

      // Compute sv that gives exact mass conservation
      const double sv_new_mc = k_new / m;

      // If needed, replace sv with mass conservative sv
      double val = abs(sv_new - sv_new_mc);
      if (val > 1.E-12)
      {
         num_corrected_cells += 1;
         sv_gf.Elem(cell_it) = sv_new_mc;
         rel_mass_corrected += val / sv_new_mc;
         if (sv_new_mc <= 0.)
         {
            cout << "cell " << cell_it << ", sv: " << sv_new_mc << endl;
            MFEM_ABORT("Invalid value for the specific volume\n");
         }
      }
   }

   // Set pct corrected to be appended to timeseries data
   pct_corrected = double(num_corrected_cells)/NDofs_L2;

   /* Append onto timeseries data */
   ts_ppd_pct_cells.Append(pct_corrected);
   ts_ppd_rel_mag.Append(rel_mass_corrected);
   sv_gf.SyncAliasMemory(S);
}


/****************************************************************************************************
* Function: ComputeDensity
* Parameters:
*  S      - BlockVector that stores mesh information, mesh velocity, and state variables.
*  rho_gf - ParGridFunction that will be set to the inverse of the specific volume
*
* Purpose:
*  Simply invert the specific volume
****************************************************************************************************/
void LagrangianLOOperator::ComputeDensity(const Vector &S, ParGridFunction &rho_gf) const
{
   Vector* sptr = const_cast<Vector*>(&S);
   ParGridFunction sv_gf;
   sv_gf.MakeRef(&L2, *sptr, block_offsets[1]);
   rho_gf.SetSize(NDofs_L2);
   for (int i = 0; i < NDofs_L2; i++)
   {
      rho_gf[i] = 1. / sv_gf.Elem(i);
   }
}


/****************************************************************************************************
* Function: GetCellStateVector
* Parameters:
*  S    - BlockVector that stores mesh information, mesh velocity, and state variables.
*  cell - Index representing the cell
*  U    - Returned vector contained the hydrodynamic state variables
*
* Purpose:
*  Given a cell, return the hydrodynamic state variables of that cell, i.e. in vector for get
*  (specific volume, velocity, specific total energy)
****************************************************************************************************/
void LagrangianLOOperator::GetCellStateVector(const Vector &S, const int cell, Vector &U) const
{
   U.SetSize(dim + 2);

   // This happens when a face is a boundary face, one of the cells index is -1q
   if (cell == -1) {
      U = 0.;
      return;
   }

   // Retrieve information from BlockVector S
   Vector* sptr = const_cast<Vector*>(&S);
   ParGridFunction sv_gf, v_gf, ste_gf;
   sv_gf.MakeRef(&L2, *sptr, block_offsets[1]);
   v_gf.MakeRef(&L2V, *sptr, block_offsets[2]);
   ste_gf.MakeRef(&L2, *sptr, block_offsets[3]);

   // Retrieve specific volume on cell
   U[0] = sv_gf.Elem(cell);

   // Retrieve cell velocity 
   // Here the corresponding gridfunction is a stacked vector:
   // Ex: x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4
   for (int i = 0; i < dim; i++)
   {
      int index = cell + i*NDofs_L2;
      assert(index < v_gf.Size());
      U[i+1] = v_gf.Elem(index);
   }

   // Retrieve cell specific total energy
   U[dim+1] = ste_gf.Elem(cell);
}


/****************************************************************************************************
* Function: SetCellStateVector
* Parameters:
*  S    - BlockVector that stores mesh information, mesh velocity, and state variables.
*  cell - Index representing the cell
*  U    - Vector containing state variables to be put into BlockVector object S.
*
* Purpose:
*  This function is used to update the BlockVector S with the newly computed hydrodynamic state
*  variables on the given cell.
****************************************************************************************************/
void LagrangianLOOperator::SetCellStateVector(Vector &S, const int cell, const Vector &U) const
{
   Array<int> dofs;
   Array<int> sub_dofs;
   dofs.Append(cell);
   sub_dofs.Append(1);

   // Get current state grid functions
   ParGridFunction sv_gf, v_gf, ste_gf;
   sv_gf.MakeRef(&L2, S, block_offsets[1]);
   v_gf.MakeRef(&L2V, S, block_offsets[2]);
   ste_gf.MakeRef(&L2, S, block_offsets[3]);

   // Set subvectors according to U
   sv_gf.SetSubVector(dofs, U[0]);
   ste_gf.SetSubVector(dofs, U[dim+1]);

   for (int i = 1; i < dim; i++)
   {
      dofs.Append(cell + i*NDofs_L2);
      sub_dofs.Append(i+1);
   }

   Vector vel;
   U.GetSubVector(sub_dofs, vel);
   v_gf.SetSubVector(dofs, vel);

   // Sync update gridfunctions
   sv_gf.SyncAliasMemory(S);
   v_gf.SyncAliasMemory(S);
   ste_gf.SyncAliasMemory(S);
}


/****************************************************************************************************
* Function: CalcOutwardNormalInt
* Parameters:
*     S    - BlockVector corresponding to nth timestep that stores mesh information, 
*            mesh velocity, and state variables.
*     cell - Inside integer cell index
*     face - Face index to compute normal on
*     res  - Resulting vector passed by reference and modified.
*
* Purpose:
*     This function returns the integral over the given face of the normal pointing away from
*     the cell c and storing this value into the vector res.
*
*        /\                                 /\
*        |                                  |
*        | n ds   = [orthogonal(s) / ||s||] | 1 ds   = orthogonal(s).
*        |                                  |
*       \/ (face)                          \/ (face)
****************************************************************************************************/
void LagrangianLOOperator::CalcOutwardNormalInt(const Vector &S, const int cell, const int face, Vector & res) const
{
   // cout << "=======================================\n"
   //      << "  Calculating outward normal integral  \n"
   //      << "=======================================\n";
   res = 0.;

   mfem::Mesh::FaceInformation FI;
   Array<int> row;
   FI = pmesh->GetFaceInformation(face);

   int c, cp, face_dof;
   c = FI.element[0].index;
   cp = FI.element[1].index;

   switch(dim)
   {
      case 1:
      {
         Vector cell_center_x(dim), face_x(dim);

         vertex_element->GetRow(face, row);

         int cell_gdof = cell + num_faces;
         geom.GetNodePositionFromBV(S, cell_gdof, cell_center_x);
         geom.GetNodePositionFromBV(S, face, face_x);

         subtract(face_x, cell_center_x, res);
         res *= 1. / res.Norml2();
         break;
      }
      case 2:
      {
         H1.GetFaceDofs(face, row);
         int face_vdof1 = row[1], face_vdof2 = row[0];
         face_dof = row[2];

         Vector face_x(dim), vdof1_x(dim), vdof2_x(dim), vdof1_v(dim), vdof2_v(dim);

         geom.GetNodePositionFromBV(S, face_dof, face_x);
         geom.GetNodePositionFromBV(S, face_vdof1, vdof1_x);
         geom.GetNodePositionFromBV(S, face_vdof2, vdof2_x);

         Vector secant(dim);
         subtract(vdof2_x, vdof1_x, secant);

         res = secant;
         geom.Orthogonal(res);

         /* Ensure orientation of normal */
         if (FI.IsInterior()) // Interior face
         {
            // Orientation of the normal vector depends on the indices of
            // the cells that share that face.  The normal points towards
            // the cell with the lower global index.  This normal must be
            // flipped if the cell we are working with is the one of the
            // lower index.
            if (cell == max(c, cp))
            {
               res *= -1.;
            }
         }

         break;
      }
      case 3:
      {
         MFEM_ABORT("3D not implemented.\n");
      }
      default:
      {
         MFEM_ABORT("Invalid dim provided.\n");
      }
   }
}



/****************************************************************************************************
* Function: ComputeIntermediateFaceVelocities
* Parameters:
*  S        - BlockVector representing FiniteElement information
*
* Purpose:
*  This function computes the intermediate face velocities as outlined in equation (5.7).  If
*  the flag is set to "testing" and a test_vel function is passed in, the velocity at the faces
*  will be the pointwise evaluation of the test_vel function.
*
*  The definition of the intermediate face velocity is given by eq (5.7).
****************************************************************************************************/
void LagrangianLOOperator::ComputeIntermediateFaceVelocities(const Vector &S) const
{
   // cout << "=======================================\n"
   //      << "   ComputeIntermediateFaceVelocities   \n"
   //      << "=======================================\n";
   
   mfem::Mesh::FaceInformation FI;
   int c, cp;
   Vector Uc(dim+2), Ucp(dim+2), n_int(dim), c_vec(dim), Vf(dim), Vf_flux(dim); 
   Vector n_vec(dim);
   double d, c_norm, F;

   Array<int> row;

   for (int face = 0; face < num_faces; face++) // face iterator
   {
      // cout << "face: " << face << endl;
      Vf = 0.;
      FI = pmesh->GetFaceInformation(face);
      c = FI.element[0].index;
      cp = FI.element[1].index;

      GetCellStateVector(S, c, Uc);

      // Get normal, d, and |F|
      CalcOutwardNormalInt(S, c, face, n_int);
      n_vec = n_int;
      F = n_vec.Norml2();
      n_vec /= F;
      if (1. - n_vec.Norml2() > 1e-12)
      {
         cout << "n_vec: ";
         n_vec.Print(cout);
         MFEM_ABORT("Invalid normal vector.\n");
      }

      if (FI.IsInterior())
      {
         GetCellStateVector(S, cp, Ucp);

         // Get max wave speed
         d = dij_sparse->Elem(c, cp);

         // Compute intermediate face velocity
         pb->velocity(Uc, Vf);
         Vector vcp; pb->velocity(Ucp, vcp);
         Vf.Add(1., vcp);
         Vf *= 0.5;
         if (use_viscosity)
         {
            double coeff = d * (Ucp[0] - Uc[0]) / F; // This is how 5.7b is defined.
            Vf.Add(coeff, n_vec);
         }
      }

      else 
      {
         assert(FI.IsBoundary());
         pb->velocity(Uc, Vf);             
      }

      // Put face velocity into object (for state var update)
      for (int i = 0; i < dim; i++)
      {
         int index = face + i * num_faces;
         v_CR_gf[index] = Vf[i];
      }
   } // End face iterator
}


/****************************************************************************************************
* Function: SetMVTargetViscCoeff
* Parameters:
*  coeff - double parameter indicating the coefficient for the viscous term to be used in
*          the target mesh velocity optimization with HiOp.  See the function SolveHiOp().
*
* Purpose:
*  To set the viscosity coefficient to be used with the TargetOptimizedMeshVelocityProblem
*  in the SolveHiOp function.
****************************************************************************************************/
void LagrangianLOOperator::SetMVTargetViscCoeff(const double & coeff)
{
   if (coeff < 0.)
   {
      MFEM_ABORT("Viscosity coefficient for TargetOptimizedMeshVelocityProblem cannot be negative.\n")
   }
   this->mv_target_visc_coeff = coeff;
}

/****************************************************************************************************
* Function: SetMVOption
* Parameters:
*  option - integer parameter indicating which mesh velocity to use.
*
* Purpose:
*  To set the mesh velocity motion computation type:
*     1) Raviart-Thomas reconstruction method.
*     2) Normal vector type
*     3) Adjacent cell face based mesh movement (LS based)
*     4) CAVEAT Weighted LS
****************************************************************************************************/
void LagrangianLOOperator::SetMVOption(const int & option)
{
   this->mv_option = option;
}


/****************************************************************************************************
* Function: SetFVOption
* Parameters:
*  option - integer parameter indicating which face velocity to use.
*
* Purpose:
*  To set the mesh velocity motion computation type:
*     1) Mass corrective.
*     2) Average of corners, Q1 type.
*     Default) Use computed velocity, RT or Vf
****************************************************************************************************/
void LagrangianLOOperator::SetFVOption(const int & option)
{
   this->fv_option = option;
}


/****************************************************************************************************
* Function: SetMVIterationOption
* Parameters:
*  option - integer parameter indicating which mesh velocity iteration type to use.
*
* Purpose:
*  To set the mesh velocity motion computation type:
*     1) Option 1 - Assume V3nperp = 0. Regular V3perp in corrective face velocity.
*                   c0 is evaluated implicitly.
*     2) Option 2 - In addition to the above, assume V3nperp is the perpendicular 
*                   component of the average of the two adjacent corner velocities. 
*                   c0 is still evaluated implicitly, c1 and V3nperp are evaluated 
*                   explicitly. Regular V3perp in corrective face velocity. 
*     3) Option 3 - The same as option 2 except the V3nperp in the corrective face 
*                   velocity is set to the average of the two adjacent corner 
*                   velocities.
*     4) Default  - No mesh iteration.
****************************************************************************************************/
void LagrangianLOOperator::SetMVIterationOption(const int &option)
{
   this->mv_it_option = option;
}


/****************************************************************************************************
* Function: SetMVIteration
* Parameters:
*  num_iterations - integer parameter indicating the number of times to iterate
*                   on the mesh velocities
*
* Purpose:
*  To enable the use of and set the number of iterations for the 
*  IterativeCornerVelocityMC() function.
****************************************************************************************************/
void LagrangianLOOperator::SetMVIteration(const int num_iterations) { 
   this->use_corner_velocity_MC_iteration = true;
   this->corner_velocity_MC_num_iterations = num_iterations; 
}


/****************************************************************************************************
* Function: SetMMViscFace
* Parameters:
*  mm_visc - Sets the amount of viscosity to add to the mesh velocity iteration
*
* Purpose:
*  To enable the use of and set the amount of viscosity to add to the 
*  cell-area based corner node mesh velocity iteration.
****************************************************************************************************/
void LagrangianLOOperator::SetMMViscFace(const double mm_visc)
{
   this->mm_visc_face = mm_visc;
}


/****************************************************************************************************
* Function: SetMMCell
* Parameters:
*  mm_visc - Sets the coefficient of the consistency term defined on the cell.
*
* Purpose:
*  To enable the use of and set the amount of viscosity to add to the 
*  cell-area based corner node mesh velocity iteration.
****************************************************************************************************/
void LagrangianLOOperator::SetMMCell(const double mm_consistency)
{
   this->mm_cell = mm_consistency;
}


/****************************************************************************************************
* Function: GetIntermediateFaceVelocity
* Parameters:
*  face - Face index
*  vel  - V_{F_m^n} corresponding to the given face
*
* Purpose:
*  Return intermediate face velocity corresponding to the given face.
*  NOTE: This function requires that ComputeIntermediateFaceVelocities has been run.
****************************************************************************************************/
void LagrangianLOOperator::GetIntermediateFaceVelocity(const int & face, Vector & vel) const
{
   assert(face < num_faces);
   // Retrieve face velocity from object
   for (int i = 0; i < dim; i++)
   {
      int index = face + i * num_faces;
      vel[i] = v_CR_gf[index];
   }
}


/****************************************************************************************************
* Function: SetCorrectedFaceVelocity
* Parameters:
*  face - Face index
*  vel  - V_{b}^3 corresponding to the given face
*
* Purpose:
*  Sets corrected face velocity corresponding to the given face.
*  NOTE: This function requires that ComputeIntermediateFaceVelocities has been run,
*        though this value will change after the corrective velocities have been computed.
****************************************************************************************************/
void LagrangianLOOperator::SetCorrectedFaceVelocity(const int & face, const Vector & vel)
{
   assert(face < num_faces);
   // Retrieve face velocity from object
   for (int i = 0; i < dim; i++)
   {
      int index = face + i * num_faces;
      v_CR_gf_corrected[index] = vel[i];
   }
}


/****************************************************************************************************
* Function: GetCorrectedFaceVelocity
* Parameters:
*  face - Face index
*  vel  - V_{b}^3 corresponding to the given face
*
* Purpose:
*  Return corrected face velocity corresponding to the given face.
*  NOTE: This function requires that ComputeIntermediateFaceVelocities has been run,
*        though this value will change after the corrective velocities have been computed.
****************************************************************************************************/
void LagrangianLOOperator::GetCorrectedFaceVelocity(const int & face, Vector & vel)
{
   assert(face < num_faces);
   // Retrieve face velocity from object
   for (int i = 0; i < dim; i++)
   {
      int index = face + i * num_faces;
      vel[i] = v_CR_gf_corrected[index];
   }
}


/****************************************************************************************************
* Function: SetCorrectedFaceFlux
* Parameters:
*  face - Face index
*  flux - V_{b}^3 corresponding to the given face
*
* Purpose:
*  Sets corrected face flux corresponding to the given face.
*  NOTE: This function requires that ComputeIntermediateFaceVelocities has been run, and that 
*        that v_CR_gf_corrected has been filled (i.e. that ComputeCorrectiveFaceVelocities() 
*        has been run).
****************************************************************************************************/
void LagrangianLOOperator::SetCorrectedFaceFlux(const int & face, const Vector & flux)
{
   assert(face < num_faces);
   // Retrieve face velocity from object
   for (int i = 0; i < dim; i++)
   {
      int index = face + i * num_faces;
      // double new_flux_comp = (v_CR_gf_fluxes[index] + flux[i]) / 2.;
      // v_CR_gf_fluxes[index] = new_flux_comp;
      v_CR_gf_fluxes[index] = flux[i];
   }
}


/****************************************************************************************************
* Function: GetCorrectedFaceFlux
* Parameters:
*  face - Face index
*  flux - V_{b}^3 corresponding to the given face
*
* Purpose:
*  Return corrected face flux corresponding to the corrected face velocity on the given face.
*  NOTE: This function requires that ComputeIntermediateFaceVelocities has been run, and that 
*        that v_CR_gf_corrected has been filled (i.e. that ComputeCorrectiveFaceVelocities() 
*        has been run).
****************************************************************************************************/
void LagrangianLOOperator::GetCorrectedFaceFlux(const int & face, Vector & flux)
{
   assert(face < num_faces);
   // Retrieve face velocity from object
   for (int i = 0; i < dim; i++)
   {
      int index = face + i * num_faces;
      flux[i] = v_CR_gf_fluxes[index];
   }
}


/****************************************************************************************************
* Function: CalcMassLoss
* Parameters:
*  S - BlockVector corresponding to nth timestep that stores mesh information, mesh velocity, 
*      and state variables.
*
* Purpose:
*  This function is used when computing values for the convergence tables.  It does this by computing
*  the sum of the local mass loss at each element.
****************************************************************************************************/
double LagrangianLOOperator::CalcMassLoss(const Vector &S)
{
   // cout << "=======================================\n"
   //      << "             CalcMassLoss              \n"
   //      << "=======================================\n";
   Vector U_i(dim + 2);

   double num = 0., denom = 0.;

   for (int ci = 0; ci < NDofs_L2; ci++)
   {
      if ((pb->get_indicator() == "ElasticNoh" || pb->get_indicator() == "Noh") && cell_bdr_flag_gf[ci] != -1)
      {
         // Skip boundary cells
         continue;
      }
      const double m = m_hpv->Elem(ci);
      const double k = pmesh->GetElementVolume(ci);
      GetCellStateVector(S, ci, U_i);
      num += abs((k / U_i[0]) - m);
      denom += abs(m);
   }

   return num / denom;
}


/****************************************************************************************************
* Function: CheckMassConservation
* Parameters:
*  S - BlockVector corresponding to nth timestep that stores mesh information, mesh velocity, 
*      and state variables.
*
* Purpose:
*  This function calculates the percentage of cells in the mesh where mass has not been conserved.
****************************************************************************************************/
void LagrangianLOOperator::CheckMassConservation(const Vector &S, ParGridFunction & mc_gf)
{
   // cout << "=======================================\n"
   //      << "         CheckMassConservation         \n"
   //      << "=======================================\n";
   Vector U_i(dim + 2);
   int counter = 0;

   double current_mass = 0., current_mass_sum = 0.;
   double num = 0., denom = 0.;
   double temp_num = 0., temp_denom = 0., val = 0.;
   double interior_num = 0., interior_denom = 0.;
   
   for (int ci = 0; ci < NDofs_L2; ci++)
   {
      const double m = m_hpv->Elem(ci);
      const double k = pmesh->GetElementVolume(ci);
      GetCellStateVector(S, ci, U_i);

      current_mass = k / U_i[0];
      temp_num = abs(current_mass - m);
      temp_denom = abs(m);

      num += temp_num;
      denom += temp_denom;
      current_mass_sum += current_mass;

      // Increment internal mass loss values
      if (cell_bdr_flag_gf[ci] == -1)
      {
         // Have interior node
         interior_num += temp_num;
         interior_denom += temp_denom;
      }

      val = temp_num / temp_denom;
      if (val > pow(10, -12))
      {
         counter++;
         // cout << "cell: " << ci << endl;
         // cout << "MASS CONSERVATION BROKEN!!\n";
         // cout << "k / U_i[0] - m = " << k / U_i[0] - m << endl;
         // cout << "K: " << k << endl;
         // cout << "T: " << U_i[0] << endl;
         // cout << "m: " << m << endl;
         // cout << "K/T: " << k / U_i[0] << endl;
         // cout << endl;
         // mc_gf[ci] = val;
      }
      // Fill corresponding cell to indicate graphically the local change in mass, if any
      mc_gf[ci] = val;
   }

   double cell_ratio = (double)counter / (double)NDofs_L2;
   double _mass_error = num / denom;
   double _interior_mass_error = interior_num / interior_denom;

   if (_mass_error > 1.E-12 || _interior_mass_error > 1.E-12 || cell_ratio > 1.E-12)
   {
      cout << "--------------------------------------\n"
           << "Percentage of cells where mass conservation was broken: " << cell_ratio << endl
           << "Initial mass sum: " << denom 
           << ", Current mass sum: " << current_mass_sum << endl
           << "Mass Error: " << _mass_error << endl
           << "Interior Mass Error: " << _interior_mass_error << endl
           << "--------------------------------------\n";
   }
}


/****************************************************************************************************
* Function: ComputeDeterminant
* Parameters:
*  C         - Dense matrix representing C_i^geo from section 5.
*  dt        - Timestep (assumes full timestep is given, not half step)
*  alpha     - Smallest postivie value satisfying equation (5.12) in [1]
*  obj_index - Nodal index, mostly to report where failure occurs
*
* Purpose:
*  This function computes the quantity given by equation (5.12) in [1].  Specifically, this function
*  computes the largest positive value alpha_i satisfying
*        det(alpha_i\mathbb{I} - \frac{dt}{2}C_i) = alpha_i^{d-1}
****************************************************************************************************/
void LagrangianLOOperator::ComputeDeterminant(const DenseMatrix &C, const double &dt, double & alpha, int obj_index)
{
   // cout << "=====================\n";
   // cout << "Computing determinant\n";
   // cout << "=====================\n";
   double trace = C.Trace();
   double det = C.Det();

   double a = 1.;
   double b = -1. * (1. + (dt / 2.) * trace);
   double c = (pow(dt, 2) / 4.) * det;

   double pm = sqrt(pow(b,2) - 4. * a * c);

   double alpha1 = -1. * b + pm;
   alpha1 /= (2. * a);

   double alpha2 = -1. * b - pm;
   alpha2 /= (2. * a);

   alpha = std::max(alpha1, alpha2);

   if (alpha <= 0.)
   {
      cout << "Failure in ComputeDeterminant at node: " << obj_index << endl;
      cout << "Dense Matrix:\n";
      C.Print(cout);
      cout << "trace: " << trace << ", det: " << det << endl;
      cout << "dt: " << dt << endl;
      cout << "a: " << a << ", b: " << b << ", c: " << c << endl;
      cout << "d1: " << alpha1 << ", d2: " << alpha2 << endl;
      MFEM_ABORT("Alpha_i should be positive.\n");
   }
}


/***********************************************************************************************************
* Function: ComputeCorrectiveFaceVelocities
* Parameters:
*   S        - BlockVector representing FiniteElement information
*   t        - Current time
*   dt       - Timestep
*   flag     - Flag used to indicate testing
*   test_vel - Velocity used for testing
* Purpose:
*   This function computes the node velocities on the faces which
*   are designed to bubble in the direction of the normal vector 
*   to conserve mass locally.
***********************************************************************************************************/
void LagrangianLOOperator::
   ComputeCorrectiveFaceVelocities(Vector &S, const double & t, const double & dt,
                                    const string flag, // Default NA
                                    void (*test_vel)(const Vector&, const double&, Vector&)) // Default NULL
{
   // cout << "==============================================\n";
   // cout << "Computing corrective interior face velocities.\n";
   // cout << "==============================================\n";

   assert(dim > 1); // No need to correct face velocities in dim=1

   /* Parameters needed for face velocity calculations */
   mfem::Mesh::FaceInformation FI;
   Vector Vf(dim), n_int(dim), n_vec(dim), face_velocity(dim);
   Vector Uc(dim+2), Ucp(dim+2);
   Vector cell_c_v(dim), cell_cp_v(dim), cell_center_v(dim); // For computation of bmn
   Array<int> row;
   int face_dof = 0, face_vdof2 = 0, c = 0, cp = 0;

   double V3n, V3nperp;

   // Iterate over faces
   for (int face = 0; face < num_faces; face++) // face iterator
   {  
      // Get intermediate face velocity, face information, and face normal
      face_velocity = 0.;
      Vf = 0.;

      FI = pmesh->GetFaceInformation(face);

      /* adjacent corner indices */
      H1.GetFaceDofs(face, row);
      int face_vdof1 = row[1], face_vdof2 = row[0], face_dof = row[2]; // preserve node orientation discussed in appendix A

      // retrieve old face and corner locations
      Vector face_x(dim), vdof1_x(dim), vdof2_x(dim), vdof1_v(dim), vdof2_v(dim);
      geom.GetNodePositionFromBV(S, face_dof, face_x);
      geom.GetNodePositionFromBV(S, face_vdof1, vdof1_x);
      geom.GetNodePositionFromBV(S, face_vdof2, vdof2_x);

      // retrieve corner velocities
      geom.GetNodeVelocity(S, face_vdof1, vdof1_v);
      geom.GetNodeVelocity(S, face_vdof2, vdof2_v);

      if (FI.IsInterior())
      {
         // Get adjacent cell information
         c = FI.element[0].index;
         cp = FI.element[1].index;
         GetCellStateVector(S, c, Uc);
         GetCellStateVector(S, cp, Ucp);
         pb->velocity(Uc, cell_c_v);
         pb->velocity(Ucp, cell_cp_v);

         // Calculate outer normal
         CalcOutwardNormalInt(S, c, face, n_int);
         n_vec = n_int;
         double F = n_vec.Norml2();
         n_vec /= F;

         if (1. - n_vec.Norml2() > 1e-12)
         {
            cout << "n_vec: ";
            n_vec.Print(cout);
            MFEM_ABORT("Invalid normal vector.\n");
         }

         // Calculate new corner locations and half  step locations
         Vector vdof1_x_new(dim), vdof2_x_new(dim), vdof1_x_half(dim), vdof2_x_half(dim);
         vdof1_x_new = vdof1_x;
         vdof1_x_new.Add(dt, vdof1_v);
         vdof1_x_half = vdof1_x;
         vdof1_x_half.Add(dt/2., vdof1_v);

         vdof2_x_new = vdof2_x;
         vdof2_x_new.Add(dt, vdof2_v);
         vdof2_x_half = vdof2_x;
         vdof2_x_half.Add(dt/2., vdof2_v);

         // calculate a_{12}^{n+1} (new tangent midpoint)
         Vector vdof12_x_new(dim);
         vdof12_x_new = vdof1_x_new;
         vdof12_x_new += vdof2_x_new;
         vdof12_x_new /= 2.;

         // Compute D (A.4c)
         Vector n_vec_R(dim), temp_vec(dim), temp_vec_2(dim);
         n_vec_R = n_vec;
         geom.Orthogonal(n_vec_R);
         subtract(vdof1_v, vdof2_v, temp_vec); // V1 - V2 = temp_vec
         subtract(vdof2_x_half, vdof1_x_half, temp_vec_2); // A2-A1
         geom.Orthogonal(temp_vec_2);
         double D = dt * (temp_vec * n_vec_R) + 2. * (n_vec * temp_vec_2);

         // Compute c1 (A.4a)
         subtract(vdof2_v, vdof1_v, temp_vec); // only change temp_vec, since temp_vec_2 is same from D calculation (half step representation)
         double c1 = ( dt * (temp_vec * n_vec) + 2. * (temp_vec_2 * n_vec_R) ) / D; // TRYING SOMETHING HERE. WILL NEED CORRECTED

         // Compute c0 (A.4b)
         Vector n_vec_half(dim);
         subtract(vdof2_x_half, vdof1_x_half, n_vec_half);
         geom.Orthogonal(n_vec_half);
         GetIntermediateFaceVelocity(face, Vf);
         
         double bmn = Vf * n_vec;
         bmn *= F;

         temp_vec = vdof1_x_half;
         geom.Orthogonal(temp_vec); // A1R
         temp_vec_2 = vdof2_x_half;
         geom.Orthogonal(temp_vec_2); // A2R
         double const1 = vdof1_v * temp_vec - vdof2_v * temp_vec_2; // V1*A1R - V2*A2R
         double const2 = vdof1_v * temp_vec_2 - vdof2_v * temp_vec; // V1*A2R - V2*A1R
         
         temp_vec = face_x;
         geom.Orthogonal(temp_vec);
         subtract(vdof2_v, vdof1_v, temp_vec_2);
         double const3 = temp_vec_2 * temp_vec; // (V2 - V1) * a3nR
         double c2 = (3. / D) * (const1 / 2. + const2 / 6. + 2. * const3 / 3.);
         
         /**
          * Compute V3n_perp
          * This value we choose:
          *    Option 1) new node is places in the normal direction
          *              of the midpoint of the secant line.
          *    Option 2) Choose perpendicular component of the 
          *              Raviart-Thomas velocity.
          *    Option 3) Take to be average of tangential components
          *              of V1 and V2.
          * 
          * When it comes to the mesh velocity iteration algorithm,
          * mv_it_options 1 and 2 use V3nperp option 1. mv_it_option 3
          * used V3nperp option 3.
          * */ 

         // Necessary quantity for either calculation
         Vector n_vec_perp(dim);
         n_vec_perp = n_vec_R;
         n_vec_perp *= -1.;

         V3nperp = 0.;
         if (mv_it_option == 3)
         {
            /* Option 3*/
            add(0.5, vdof2_v, 0.5, vdof1_v, temp_vec);
            V3nperp = temp_vec * n_vec_perp;
         }
         else 
         {
            /* Option 1 */ 
            subtract(vdof2_x_new, vdof1_x_new, temp_vec);
            subtract(face_x, vdof12_x_new, temp_vec_2);
            temp_vec_2.Add((c2+ (3. * bmn / D))*dt, n_vec);
            const1 = temp_vec * temp_vec_2; // numerator
            temp_vec_2 = n_vec_R;
            temp_vec_2 *= -1.;
            temp_vec_2.Add(c1, n_vec);
            const2 = temp_vec * temp_vec_2;
            const2 *= dt; // denominator
            V3nperp = -1. * const1 / const2;
         }

         /* Option 2 */
         // Vector rt_face_vel(dim);
         // geom.GetNodeVelocity(S, face_dof, rt_face_vel);
         // V3nperp = rt_face_vel * n_vec_perp;

         // Compute V3n (5.11)
         V3n = c1 * V3nperp + c2 + 3. * bmn / D;

         // Add in normal and tangentional contributions
         face_velocity = 0.;
         face_velocity.Add(V3n, n_vec);
         face_velocity.Add(V3nperp, n_vec_perp); // 10/2

         // Check perpendicular only in the case where new 
         // nodal location is exactly above tangent location
         // Vector face_x_new(dim);
         // face_x_new = face_x;
         // face_x_new.Add(dt, face_velocity);
         // subtract(face_x_new, vdof12_x_new, temp_vec);
         // subtract(vdof2_x_new, vdof1_x_new, temp_vec_2);
         // if (abs(temp_vec * temp_vec_2) > pow(10, -12))
         // {
         //    MFEM_ABORT("vectors are not orthogonal!\n");
         // }
      } // End interior face

      // On the boundary, we just need the average of the adjacent 
      // vertices to account for the corrective corner velocities
      // NOTE that even in the test case, the corner velocities will
      // be modified and so prescribing the exact velocities on the
      // boundary faces will not preserve mass.
      else
      {
         assert(FI.IsBoundary());
         // TODO: Add in boundary conditions similar to corner vertex bcs
         for (int j = 0; j < dim; j++)
         {
            face_velocity[j] = (vdof1_v[j] + vdof2_v[j]) / 2.;
         }
        
      } // End boundary face

      // Lastly, put face velocity into gridfunction object
      if (face_velocity[0] != face_velocity[0] || face_velocity[1] != face_velocity[1])
      {
         MFEM_ABORT("NaN values encountered in corrective face velocity calculation.\n");
      }

      // cout << "compute corrective face velocity face (" << face << ") dof: " << face_dof << ": ";
      // face_velocity.Print(cout);
      geom.UpdateNodeVelocity(S, face_dof, face_velocity);
      SetCorrectedFaceVelocity(face, face_velocity);
   }
}


/***********************************************************************************************************
* Function: ComputeCorrectiveFaceFluxes
* Parameters:
*   S        - BlockVector representing FiniteElement information
*   t        - Current time
*   dt       - Timestep
* Purpose:
*   This function computes the geometric face node velocities which
*   are used in the iterative procedure to recompute corner node 
*   velocities.
***********************************************************************************************************/
void LagrangianLOOperator::
   ComputeCorrectiveFaceFluxes(Vector &S, const double & t, const double & dt)
{
   assert(dim > 1); // There is no need for this function in dim=1
   // cout << "========================================================\n";
   // cout << "Computing corrective interior geometric face velocities.\n";
   // cout << "========================================================\n";
   /* Parameters needed for face velocity calculations */
   mfem::Mesh::FaceInformation FI;
   Vector Vf(dim), n_int(dim), n_vec(dim), face_velocity(dim);
   Vector cell_c_v(dim), cell_cp_v(dim), cell_center_v(dim); // For computation of bmn
   Array<int> row;
   int face_dof = 0, face_vdof2 = 0, c = 0, cp = 0;

   double V3n, V3nperp;

   // Iterate over faces
   for (int face = 0; face < num_faces; face++) // face iterator
   {  
      // Get intermediate face velocity, face information, and face normal
      face_velocity = 0.;
      Vf = 0.;

      FI = pmesh->GetFaceInformation(face);
      GetCorrectedFaceVelocity(face, face_velocity);

      /* adjacent corner indices */
      H1.GetFaceDofs(face, row);
      int face_dof = row[2];

      if (FI.IsInterior())
      {
         // Get adjacent cell information
         c = FI.element[0].index;
         cp = FI.element[1].index;

         // Calculate outer normal
         CalcOutwardNormalInt(S, c, face, n_int);
         n_vec = n_int;
         double F = n_vec.Norml2();
         n_vec /= F;

         if (1. - n_vec.Norml2() > 1e-12)
         {
            cout << "n_vec: ";
            n_vec.Print(cout);
            MFEM_ABORT("Invalid normal vector.\n");
         }

         /*** Compute Face Flux ***/

         // Compute C_i
         DenseMatrix Ci;
         // mfem::Array<int> elements;
         // face_element->GetRow(face_dof, elements);
         ComputeCiGeo(v_geo_gf, face_dof, Ci);

         // Compute alpha_i
         double alpha_i;
         ComputeDeterminant(Ci, dt, alpha_i, face);

         // Compute flux
         Vector face_flux(dim);
         DenseMatrix _mat(Ci);
         _mat *= - dt / 2.;
         for (int i = 0; i < dim; i++)
         {
            _mat(i,i) += alpha_i;
         }

         _mat.Mult(face_velocity, face_flux);
         SetCorrectedFaceFlux(face, face_flux); 
      }
   }
}

/****************************************************************************************************
* Function: SetCellCenterAsCenter
* Parameters:
*   S        - BlockVector representing FiniteElement information
*
* Purpose:
*  This function is only needed since we are not able to implement the Serendipity finite elements.
*  Since the mesh motion does not depends on a value for the node corresponding to the center of the 
*  element, we fill the actual grid function with the average of the moved adjacent corners
****************************************************************************************************/
void LagrangianLOOperator::SetCellCenterAsCenter(Vector &S)
{
   x_gf.MakeRef(&H1, S, block_offsets[0]);

   // Since we cannot use Serendipity elements, we must update cell center velocities
   Vector cell_x(dim), node_x(dim);

   Array<int> verts;
   for (int ci = 0; ci < NDofs_L2; ci++)
   {
      // Get center node dof
      int cell_vdof;
      switch (dim)
      {
         case 1:
         {
            cell_vdof = num_faces + ci;
            break;
         }
         case 2:
         {
            cell_vdof = NVDofs_H1 + num_faces + ci;
            break;
         }
         case 3:
         {
            MFEM_ABORT("3D not implemented.\n");
         }
         default:
         {
            MFEM_ABORT("Incorrect dim value provided.\n");
         }
      }

      cell_x = 0.;
      pmesh->GetElementVertices(ci, verts);

      for (int j = 0; j < verts.Size(); j++)
      {
         geom.GetNodePositionFromBV(S, verts[j], node_x);
         cell_x += node_x;
      }

      cell_x /= verts.Size();

      geom.UpdateNodePosition(S, cell_vdof, cell_x);
   }
}


/****************************************************************************************************
* Function: FillCenterVelocitiesWithL2
* Parameters:
*   S       - BlockVector containing cell L2 velocity
*   dSdt    - BlockVector containing derivative information,
*             containing mesh velocities
*
* Purpose:
*  This function is only needed since we are not able to implement 
*  the Serendipity finite elements. Since the mesh motion does not 
*  depends on a value for the node corresponding to the center of 
*  the element, we need to fill this with a velocity which will 
*  ensure the center node remains inside the cell.  This is done 
*  by taking the hydrodynamic velocity at the cell.
****************************************************************************************************/
void LagrangianLOOperator::FillCenterVelocitiesWithL2(const Vector &S, Vector &dSdt) const
{
   // Since we cannot use Serendipity elements, we must update cell center velocities
   ParGridFunction dxdt;
   dxdt.MakeRef(&H1, dSdt, block_offsets[0]);

   Vector Uc(dim+2), node_v(dim);

   for (int ci = 0; ci < NDofs_L2; ci++)
   {
      // Get center node dof
      int cell_vdof;
      switch (dim)
      {
         case 1:
         {
            cell_vdof = num_faces + ci;
            break;
         }
         case 2:
         {
            cell_vdof = NVDofs_H1 + num_faces + ci;
            break;
         }
         case 3:
         {
            MFEM_ABORT("3D not implemented.\n");
         }
         default:
         {
            MFEM_ABORT("Incorrect dim value provided.\n");
         }
      }

      GetCellStateVector(S, ci, Uc);
      pb->velocity(Uc, node_v);
      
      // Get corresponding velocity from ParGridFunction
      geom.UpdateNodeVelocity(dSdt, cell_vdof, node_v);
   }
}


/****************************************************************************************************
* Function: FillCenterVelocitiesWithAvg
* Parameters:
*   S        - BlockVector representing FiniteElement information
*
* Purpose:
*  This function is only needed since we are not able to implement 
*  the Serendipity finite elements. Since the mesh motion does not 
*  depends on a value for the node corresponding to the center of 
*  the element, we need to fill this with a velocity which will 
*  ensure the center node remains inside the cell.  This is done 
*  by average the values at the four corner nodes of the cell.
****************************************************************************************************/
void LagrangianLOOperator::FillCenterVelocitiesWithAvg(Vector &dxdt) const
{
   // Since we cannot use Serendipity elements, we must update cell center velocities
   Vector Vc(dim), Uc(dim+2), node_v(dim);

   Array<int> verts;
   for (int ci = 0; ci < NDofs_L2; ci++)
   {
      // Get center node dof
      int cell_vdof;
      switch (dim)
      {
         case 1:
         {
            cell_vdof = num_faces + ci;
            break;
         }
         case 2:
         {
            cell_vdof = NVDofs_H1 + num_faces + ci;
            break;
         }
         case 3:
         {
            MFEM_ABORT("3D not implemented.\n");
         }
         default:
         {
            MFEM_ABORT("Incorrect dim value provided.\n");
         }
      }
      
      Vc = 0.;

      if (dim == 1)
      {
         pmesh->GetElementVertices(ci, verts);
         for (int j = 0; j < verts.Size(); j++)
         {
            geom.GetNodeVelocity(dxdt, verts[j], node_v);
            Vc += node_v;
         }
         Vc /= verts.Size();
      }
      // Average velocity from corner vertices
      else if (dim == 2)
      {
         pmesh->GetElementVertices(ci, verts);
         for (int j = 0; j < verts.Size(); j++)
         {
            geom.GetNodeVelocity(dxdt, verts[j], node_v);
            Vc += node_v;
         }

         Vc /= verts.Size();
      }

      geom.UpdateNodeVelocity(dxdt, cell_vdof, Vc);
   }
}


/****************************************************************************************************
* Function: FillFaceVelocitiesWithAvg
* Parameters:
*   S        - BlockVector representing FiniteElement information
*
* Purpose:
*  This function is merely for testing purposes and should not be implemented in the full program
*  as the purpose of moving the face nodes is to ensure local mass conservation.
****************************************************************************************************/
void LagrangianLOOperator::FillFaceVelocitiesWithAvg(ParGridFunction &dxdt) const
{
   /* Parameters needed for face velocity calculations */
   mfem::Mesh::FaceInformation FI;
   Vector Uc(dim+2), Ucp(dim+2), face_velocity(dim);
   Vector cell_c_v(dim), cell_cp_v(dim); // For computation of bmn
   Array<int> row;
   int face_dof = 0, face_vdof1 = 0, face_vdof2 = 0;

   // Iterate over faces
   for (int face = 0; face < num_faces; face++) // face iterator
   {  
      face_velocity = 0.;
      FI = pmesh->GetFaceInformation(face);

      /* adjacent corner indices */
      H1.GetFaceDofs(face, row);
      int face_vdof1 = row[1], face_vdof2 = row[0], face_dof = row[2]; // preserve node orientation discussed in appendix A

      // retrieve old face and corner locations
      Vector vdof1_v(dim), vdof2_v(dim);

      // retrieve corner velocities
      geom.GetNodeVelocity(dxdt, face_vdof1, vdof1_v);
      geom.GetNodeVelocity(dxdt, face_vdof2, vdof2_v);

      // Average nodal velocities
      for (int j = 0; j < dim; j++)
      {
         face_velocity[j] = (vdof1_v[j] + vdof2_v[j]) / 2;
      }

      // Lastly, put face velocity into gridfunction object
      geom.UpdateNodeVelocity(dxdt, face_dof, face_velocity);
   }
}


/****************************************************************************************************
* Function: FillFaceVelocitiesWithButterfly
* Parameters:
*   S        - BlockVector representing FiniteElement information
*
* Purpose:
****************************************************************************************************/
void LagrangianLOOperator::FillFaceVelocitiesWithButterfly(ParGridFunction &dxdt) const
{
   mfem::Mesh::FaceInformation FI;
   Vector face_velocity(dim), tmp_vel(dim), vdof1_v(dim), vdof2_v(dim);
   Vector macro_v1(dim), macro_v2(dim), sub_vel(dim);
   Array<int> row, edge_row, cell_faces, side_faces, tmp_row;
   bool is_singular_macro = false;
   int sub_index;

   for (int face = 0; face < num_faces; face++)
   {
      // cout << "======================\n";
      // cout << "face: " << face << endl;
      face_velocity = 0.;
      FI = pmesh->GetFaceInformation(face);

      /* adjacent corner indices */
      H1.GetFaceDofs(face, row);
      int face_vdof1 = row[0], face_vdof2 = row[1], face_dof = row[2];

      // retrieve corner velocities
      geom.GetNodeVelocity(dxdt, face_vdof1, vdof1_v);
      geom.GetNodeVelocity(dxdt, face_vdof2, vdof2_v);

      /* If boundary cell, take average of adjacent vertices */
      if (FI.IsBoundary())
      {
         for (int j = 0; j < dim; j++)
         {
            face_velocity[j] = (vdof1_v[j] + vdof2_v[j]) / 2.;
         }
      }
      else {
         assert(FI.IsInterior());
         is_singular_macro = false;
         /* Get faces of adjacent cells ---> side_faces */
         int c = FI.element[0].index, cp = FI.element[1].index;
         element_face.GetRow(c, side_faces);
         element_face.GetRow(cp, cell_faces);
         side_faces.Append(cell_faces);
         // cout << "side faces: ";
         // side_faces.Print(cout);

         /* adj node 1, comput macro_v1 */
         vertex_edge.GetRow(face_vdof1, edge_row);
         if (edge_row.Size() == 4)
         {
            add(0.75, vdof1_v, 0.375, vdof2_v, macro_v1);
            // cout << "node 1: " << face_vdof1 << "\n"
            //      << "\tadj edges: ";
            // edge_row.Print(cout);
            /* Find adj face not in side faces */
            for (int vert_adj_face_it = 0; vert_adj_face_it < edge_row.Size(); vert_adj_face_it++)
            {
               int vert_adj_face = edge_row[vert_adj_face_it];
               if (side_faces.Find(vert_adj_face) == -1)
               {
                  /* face not found */
                  // cout << "face " << vert_adj_face << " is the face you seek!!!\n";

                  /* Get sub vel */
                  H1.GetFaceDofs(vert_adj_face, tmp_row);
                  sub_index = (tmp_row[0] == face_vdof1) ? tmp_row[1] : tmp_row[0];
                  // cout << "face dofs we seek. Face: " << vert_adj_face << ". vdof1: " << tmp_row[0] << ". vdof2: " << tmp_row[1] << ". Subindex: " << sub_index << endl;
                  geom.GetNodeVelocity(dxdt, sub_index, sub_vel);
                  // cout << "sub index: " << sub_index << ", vel: ";
                  // sub_vel.Print(cout);

                  macro_v1.Add(-.125, sub_vel);
                  // cout << "macro_v1: ";
                  // macro_v1.Print(cout);
               }
            }
         }
         else
         {
            // cout << "node 1 is not an interior node.\n";
            is_singular_macro = true;
            macro_v1 = 0.;
         }
         

         /* adj node 2, comput macro_v2 */
         vertex_edge.GetRow(face_vdof2, edge_row);
         if (edge_row.Size() == 4)
         {
            add(0.75, vdof2_v, 0.375, vdof1_v, macro_v2);
            // cout << "node 1: " << face_vdof2 << "\n"
            //      << "\tadj edges: ";
            // edge_row.Print(cout);
            /* Find adj face not in side faces */
            for (int vert_adj_face_it = 0; vert_adj_face_it < edge_row.Size(); vert_adj_face_it++)
            {
               int vert_adj_face = edge_row[vert_adj_face_it];
               if (side_faces.Find(vert_adj_face) == -1)
               {
                  /* face not found */
                  // cout << "face " << vert_adj_face << " is the face you seek!!!\n";

                  /* Get sub vel */
                  H1.GetFaceDofs(vert_adj_face, tmp_row);
                  sub_index = (tmp_row[0] == face_vdof2) ? tmp_row[1] : tmp_row[0];
                  // cout << "face dofs we seek. Face: " << vert_adj_face << ". vdof1: " << tmp_row[0] << ". vdof2: " << tmp_row[1] << ". Subindex: " << sub_index << endl;
                  geom.GetNodeVelocity(dxdt, sub_index, sub_vel);
                  // cout << "sub index: " << sub_index << ", vel: ";
                  // sub_vel.Print(cout);

                  macro_v2.Add(-.125, sub_vel);
                  // cout << "macro_v1: ";
                  // macro_v2.Print(cout);
               }
            }
         }
         else
         {
            // cout << "node 2 is not an interior node.\n";
            is_singular_macro = true;
            macro_v2 = 0.;
         }

         /* Lastly put macro_v1 and macro_v2 together */
         if (!is_singular_macro)
         {
            add(0.5, macro_v1, 0.5, macro_v2, face_velocity);
         }
         else
         {
            // cout << "---had a singular macro---\n";
            add(1., macro_v1, 1., macro_v2, face_velocity);
         }
      }
      

      // Lastly, put face velocity into gridfunction object
      if (face_velocity[0] != face_velocity[0] || face_velocity[1] != face_velocity[1])
      {
         MFEM_ABORT("NaN values encountered in corrective face velocity calculation.\n");
      }

      geom.UpdateNodeVelocity(dxdt, face_dof, face_velocity);
   }
}


/****************************************************************************************************
* Function: SetViGeo
* Parameters:
*   node - Global index of node in question.
*   vel  - Velocity Vi^geo at given node
*
* Purpose:
*  This function sets ViGeo.
****************************************************************************************************/
void LagrangianLOOperator::SetViGeo(const int &node, const Vector &vel)
{
   MFEM_ABORT("Function deprecated.\n");
   for (int i = 0; i < dim; i++)
   {
      int index = node + i * NVDofs_H1;
      v_geo_gf[index] = vel[i];
   }
}


/****************************************************************************************************
* Function: GetViGeo
* Parameters:
*   node - Global index of node in question.
*   vel  - Velocity Vi^geo at given node
*
* Purpose:
*  This function returns Vi_geo.
*  NOTE: Function LagrangianLOOperator::ComputeGeoVRaviart must be called first.
****************************************************************************************************/
void LagrangianLOOperator::GetViGeo(const int & node, Vector & vel)
{
   for (int i = 0; i < dim; i++)
   {
      int index = node + i * NVDofs_H1;
      vel[i] = v_geo_gf[index];
   }
}


/****************************************************************************************************
* Function: UpdateMesh
* Parameters:
*   S    - BlockVector representing FiniteElement information
*
* Purpose:
*  This function returns the cartesian location corresponding to a global node.
****************************************************************************************************/
void LagrangianLOOperator::UpdateMesh(const Vector & S) const 
{
   Vector* sptr = const_cast<Vector*>(&S);
   x_gf.MakeRef(&H1, *sptr, block_offsets[0]);
   H1.GetMesh()->NewNodes(x_gf, false);
   L2.GetMesh()->NewNodes(x_gf, false);
   pmesh->NewNodes(x_gf, false);
}


/****************************************************************************************************
* Function: SaveStateVecsToFile
* Parameters:
*  S                  - BlockVector representing FiniteElement information
*  output_file_prefix - Parameters used to define output file location
*  output_file_suffix - See above.
*
* Purpose:
*  The results of this program can be used to plot the same problem at multiple refinements.
****************************************************************************************************/
void LagrangianLOOperator::SaveStateVecsToFile(const Vector &S, 
                                                    const string &output_file_prefix, 
                                                    const string &output_file_suffix)
{
   Vector center(dim), U(dim+2), vel(dim);
   double pressure=0., ss=0.;

   // Form filenames and ofstream objects
   std::string sv_file = output_file_prefix + output_file_suffix;
   std::ofstream fstream_sv(sv_file.c_str());
   // Header
   fstream_sv << "x,y,z,rho,vx,vy,vz,ste,sie,p,ss,cell_type";
   if (use_elasticity)
   {
      fstream_sv << ",sd11,sd12,sd13,sd21,sd22,sd23,sd31,sd32,sd33,s11,s22,s33,es";
   }
   fstream_sv << "\n";

   for (int i = 0; i < NDofs_L2; i++)
   {
      // compute pressure and sound speed on the fly
      GetCellStateVector(S, i, U);
      pb->velocity(U, vel);
      double rho = 1. / U[0];
      double e_sheer = 0.;
      if (use_elasticity) { e_sheer = elastic->e_sheer(i); }
      double sie = pb->specific_internal_energy(U, e_sheer);
      double attr = pmesh->GetAttribute(i);
      pressure = pb->pressure(rho, sie, attr);
      ss = pb->sound_speed(rho, pressure, attr);
      pmesh->GetElementCenter(i, center);

      // Get relevant quantities
      double _x = center[0], _y = 0., _z = 0;
      double _vx = vel[0], _vy = 0., _vz = 0.;
      if (dim > 1) {
         _y = center[1];
         _vy = vel[1];
         if (dim > 2) {
            _z = center[2];
            _vz = vel[2];
         }
      }
      
      // Output to file
      fstream_sv << _x << "," << _y << "," << _z << ","    // x,y,z,
                 << 1./U[0] << ","                         // rho
                 << _vx << "," << _vy << "," << _vz << "," // vx,vy,vz
                 << U[dim+1] << ","                        // ste
                 << sie << ","                             // sie
                 << pressure << ","                        // pressure
                 << ss << ",";                             // sound speed

      // Print flag if interior or bdr
      if (cell_bdr_flag_gf[i] == -1.)
      {
         fstream_sv << "int";                             // cell_type
      }
      else 
      {
         fstream_sv << "bdr";
      }
 
      if (use_elasticity)
      {
         DenseMatrix sigmaD(3);
         ComputeSigmaDComp(S, i, sigmaD); // sigma
         // for (int i = 0; i < dim; i++) { fstream_sv << "," << sigma(0,i); }
         fstream_sv << "," << sigmaD(0,0) << "," << sigmaD(0,1) << "," << sigmaD(0,2)
                    << "," << sigmaD(1,0) << "," << sigmaD(1,1) << "," << sigmaD(1,2)
                    << "," << sigmaD(2,0) << "," << sigmaD(2,1) << "," << sigmaD(2,2)
                    << "," << sigmaD(0,0) - pressure
                    << "," << sigmaD(1,1) - pressure
                    << "," << sigmaD(2,2) - pressure
                    << "," << e_sheer;
      }
      fstream_sv << "\n";
   }
}

/****************************************************************************************************
* Function: SaveTimeSeriesArraysToFile
* Parameters:
*  output_file_prefix - Parameters used to define output file location
*  output_file_suffix - See above.
*
* Purpose:
*  Save timeseries arrays to file.
****************************************************************************************************/
void LagrangianLOOperator::SaveTimeSeriesArraysToFile(const string &output_file_prefix, const string &output_file_suffix)
{
   const int num_timesteps = ts_timestep.Size();

   // Verify all timeseries arrays are the proper size
   assert(ts_min_detJ.Size() == num_timesteps);
   assert(ts_min_detJ_cell.Size() == num_timesteps);

   if (post_process_density)
   {
      assert(ts_ppd_pct_cells.Size() == num_timesteps);
      assert(ts_ppd_rel_mag.Size() == num_timesteps);
   }

   // Form filenames and ofstream objects
   std::string sv_file = output_file_prefix + output_file_suffix;
   std::ofstream fstream_sv(sv_file.c_str());
   fstream_sv << "timestep,t,dt,dij_max,dij_avg,minDetJ,minDetJCell";

   /* Add optional entries dependent on command line options and problem */
   if (post_process_density)
   {
      fstream_sv << ",ppdPctCells,ppdRelMag";
   }
   if (pb->get_indicator() == "Kidder")
   {
      fstream_sv << ",avgInteriorRadius,avgExteriorRadius,avgDensity,avgEntropy";
   }
   fstream_sv << "\n";


   for (int i = 0; i < num_timesteps; i++)
   {
      // Output to file
      fstream_sv << i << ","
                 << ts_t[i] << ","
                 << ts_timestep[i] << ","
                 << ts_dijmax[i] << ","
                 << ts_dijavg[i] << ","
                 << ts_min_detJ[i] << ","
                 << ts_min_detJ_cell[i];

      // Save ppd values if enabled
      if (post_process_density)
      {
         fstream_sv << ","
                    << ts_ppd_pct_cells[i] << ","
                    << ts_ppd_rel_mag[i];
      }
      // Save interior and exterior radius in case of Kidder problem
      if (pb->get_indicator() == "Kidder")
      {
         fstream_sv << ","
                    << ts_kidder_avg_rad_int[i] << ","
                    << ts_kidder_avg_rad_ext[i] << ","
                    << ts_kidder_avg_density[i] << ","
                    << ts_kidder_avg_entropy[i];
      }
      fstream_sv << "\n";
   }
}


/****************************************************************************************************
* Function: Normal direction based Velocity Functions
****************************************************************************************************/

/****************************************************************************************************
* Function: tensor
* Parameters:
*  v1     - Vector object
*  v2     - Vector object
*  dm     - Resultant DenseMatrix
*
* Purpose:
*  Compute the tensor product of two vectors.
*        
****************************************************************************************************/
void LagrangianLOOperator::tensor(const Vector & v1, const Vector & v2, DenseMatrix & dm) const
{
   const int v1_len = v1.Size(), v2_len = v2.Size();
   for (int i = 0; i < v1_len; i++)
   {
      for (int j = 0; j < v2_len; j++)
      {
         dm.Elem(i,j) = v1[i] * v2[j];
      }
   }
}


/****************************************************************************************************
* Function: ComputeGeoVNormal
* Parameters:
*  S    - BlockVector corresponding to nth timestep that stores mesh information, 
*         mesh velocity, and state variables.
*
* Purpose:
*  This function computes the geometric velocity from the adjacent face normals
*  and face velocities as follows:
*
*  Define Mi    = Sum (nf x nf)
*         Ri    = Sum [(nf x nf) Vf]
*         Vigeo = Mi^-1 Ri
*
*  Note: Vigeo is only the geometric velocity, and is not the velocity used 
*        to move the mesh.  For the mesh velocity, one must solve for the 
*        alpha_i, which process is found in ComputeNodeVelocitiesFromVgeo().
*
*  Note: This function fills the Q1 velocity field v_geo_gf and also the mv_gf
*        in the BlockVector S.
****************************************************************************************************/
void LagrangianLOOperator::ComputeGeoVNormal(const Vector &S, ParGridFunction &mv_gf_l) const
{
   mfem::Mesh::FaceInformation FI;
   Vector node_v(dim), Ri(dim), n_vec(dim), Vf(dim), y(dim);
   DenseMatrix Mi(dim), dm_tmp(dim);

   Array<int> faces_row, oris;
   int faces_length;
   double c_norm;

   // Iterate over all corner mesh velocity dofs
   for (int node = 0; node < NVDofs_H1; node++) // Vertex iterator
   {
      /* Reset Mi, Ri, Vi */
      Mi = 0., Ri = 0., node_v = 0.;

      /* Get cell faces */
      vertex_edge.GetRow(node, faces_row);
      faces_length = faces_row.Size();

      /* iterate over adjacent faces */
      for (int face_it = 0; face_it < faces_length; face_it++) // Adjacent face iterator
      {
         int face = faces_row[face_it];
         GetIntermediateFaceVelocity(face, Vf);
         FI = pmesh->GetFaceInformation(face);

         // Retrieve corresponding normal (orientation doesn't matter)
         CalcOutwardNormalInt(S, FI.element[0].index, face, n_vec);
         c_norm = n_vec.Norml2();
         n_vec /= c_norm;

         tensor(n_vec, n_vec, dm_tmp);
         Mi += dm_tmp;
         dm_tmp.Mult(Vf, y);
         Ri += y;
      }

      /* Solve for node_v */
      Mi.Invert();
      Mi.Mult(Ri, node_v);

      /* In every case, update the corresponding nodal velocity in S */
      geom.UpdateNodeVelocityVecL(mv_gf_l, node, node_v);

   } // End node iterator
}


/****************************************************************************************************
* Function: ComputeGeoVNormalDistributedViscosity
* Parameters:
*  S    - BlockVector corresponding to nth timestep that stores mesh information, 
*         mesh velocity, and state variables.
*
* Purpose:
*  This function computes the geometric velocity from the adjacent face normals
*  and face velocities as follows:
*
*  Define Mi    = Sum (nf x nf)
*         Ri    = Sum [(nf x nf) Vf]
*         Vigeo = Mi^-1 Ri
*
*  Note: Vigeo is only the geometric velocity, and is not the velocity used 
*        to move the mesh.  For the mesh velocity, one must solve for the 
*        alpha_i, which process is found in ComputeNodeVelocitiesFromVgeo().
*
*  Note: This function fills the Q1 velocity field v_geo_gf and also the mv_gf
*        in the BlockVector S.
****************************************************************************************************/
void LagrangianLOOperator::ComputeGeoVNormalDistributedViscosity(Vector &S)
{
   mfem::Mesh::FaceInformation FI;
   Vector node_v(dim), Ri(dim), n_vec(dim), Vf(dim), y(dim);
   DenseMatrix Mi(dim), dm_tmp(dim);

   Array<int> faces_row, oris;
   int faces_length;
   double c_norm, d, pcp, pc, coeff;

   DofEntity entity;
   int EDof, c, cp;

   Vector Uc(dim+2), Ucp(dim+2), vcp(dim), visc_contr(dim);

   // Iterate over all mesh velocity dofs
   for (int node = 0; node < NDofs_H1; node++) // Vertex iterator
   {
      node_v = 0.; // Reset node velocity
      GetEntityDof(node, entity, EDof);

      switch (entity)
      {
         case 0: // corner
         {
            // Reset Mi, Ri, Vi
            Mi = 0., Ri = 0., node_v = 0., visc_contr = 0.;

            // Get cell faces
            vertex_edge.GetRow(node, faces_row);
            faces_length = faces_row.Size();

            // iterate over adjacent faces
            for (int face_it = 0; face_it < faces_length; face_it++) // Adjacent face iterator
            {
               /* Just the velocity */
               Vf = 0.; 
               int face = faces_row[face_it];
               FI = pmesh->GetFaceInformation(face);

               c = FI.element[0].index;
               cp = FI.element[1].index;
               GetCellStateVector(S,c,Uc);
               // Retrieve corresponding normal (orientation doesn't matter)
               CalcOutwardNormalInt(S, c, face, n_vec);
               c_norm = n_vec.Norml2();
               n_vec /= c_norm;

               if (FI.IsInterior())
               {
                  GetCellStateVector(S,cp,Ucp);
                  d = dij_sparse->Elem(c,cp);

                  /* Compute average velocity across face */
                  pb->velocity(Uc, Vf);
                  pb->velocity(Ucp, vcp);
                  Vf.Add(1., vcp);
                  Vf *= 0.5;
               }
               else
               {
                  assert(FI.IsBoundary());
                  pb->velocity(Uc, Vf);
               }
               
               tensor(n_vec, n_vec, dm_tmp);
               Mi += dm_tmp;
               dm_tmp.Mult(Vf, y);
               Ri += y;

               /* Get viscosity contribution */
               coeff = 0.5 * d * (Ucp[0] - Uc[0]) / c_norm; // This is how 5.7b is defined.
               if (BdrVertexIndexingArray[node] == -1)
               {
                  visc_contr.Add(0.5 * coeff, n_vec);
               }
               else 
               {
                  visc_contr.Add(coeff, n_vec);
               }
            }

            Mi.Invert();
            Mi.Mult(Ri, node_v);
            node_v.Add(1., visc_contr);

            // Only update the v_geo_gf in the case of the corners
            // Using Q1 for this velocity field
            SetViGeo(node, node_v);

            break;
         }
         case 1: // face
         {
            // face nodes just set to ivf
            assert(face >= 0);
            GetIntermediateFaceVelocity(EDof, node_v);

            break;
         }
         case 2: // Cell Center
         {
            Vector Uc(dim+2);
            GetCellStateVector(S, EDof, Uc);
            pb->velocity(Uc, node_v);

            break;
         }
         default:
         {
            MFEM_ABORT("Invalid entity\n");
         }
      }
      // cout << "node velocity from mv2 for node (" << node << "): ";
      // node_v.Print(cout);
      // In every case, update the corresponding nodal velocity in S
      geom.UpdateNodeVelocity(S, node, node_v);

   } // End node iterator
}


/****************************************************************************************************
* Function: ComputeGeoVCellFaceNormal
* Parameters:
*  S    - BlockVector corresponding to nth timestep that stores mesh information, 
*         mesh velocity, and state variables.
*
* Purpose:
*  This function computes the geometric velocity from the face normals and face 
*  velocities from all adjacent cells as follows:
*
*  CellFaceNormalTODO:
*        | num cells | num faces |
*        |-----------|-----------|
*        |     1     |     4     |
*        |     2     |     7     |
*        |     4     |    12     |
*
*  Note: Vigeo is only the geometric velocity, and is not the velocity used 
*        to move the mesh.  For the mesh velocity, one must solve for the 
*        alpha_i, which process is found in ComputeNodeVelocitiesFromVgeo().
*
*  Note: This function fills the Q1 velocity field v_geo_gf and also the mv_gf
*        in the BlockVector S.
****************************************************************************************************/
void LagrangianLOOperator::ComputeGeoVCellFaceNormal(Vector &S)
{
   cout << "ComputeGeoVCellFaceNormal call\n";
   mfem::Mesh::FaceInformation FI;
   
   Vector node_v(dim), node_x(dim), n_vec(dim), Vf(dim);
   Vector vdof1_x(dim), vdof2_x(dim), face_x(dim), vdof_avg(dim);
   double F;

   Array<int> element_row, cell_faces_row, oris, face_dofs;
   int row_index, length_element_row, length_cell_faces, num_rows_A, num_cols_A = 6;

   DenseMatrix A, AT, ATA(num_cols_A);
   Vector RHS, RHS_mod(num_cols_A), res(num_cols_A), temp_vec;

   DofEntity entity;
   int EDof;

   std::unordered_set<int> uniqueFaceIndices;

   // Iterate over all mesh velocity dofs
   // cout << "ndofs_h1: " << NDofs_H1 << endl;
   for (int node = 0; node < NDofs_H1; node++) // Vertex iterator
   {
      // cout << "NODE: " << node << endl;
      node_v = 0.; // Reset node velocity
      GetEntityDof(node, entity, EDof);
      // cout << "EDof: " << EDof << endl;

      switch (entity)
      {
         case 0: // corner
         {
            // Reset vals
            uniqueFaceIndices.clear();
            F = 0.;
            geom.GetNodePositionFromBV(S, node, node_x);

            // Get adjacent cells
            vertex_element->GetRow(node, element_row);
            // cout << "element_row: ";
            // element_row.Print(cout);
            length_element_row = element_row.Size();

            // Set size of overdetermined matrix A depending on number of adjacent cells
            switch (length_element_row)
            {
               case 1:
               {
                  // TODO: What is to be done at the corners of our domain?
                  node_v = 1.;
                  continue;
               }
               case 2:
               {
                  // TODO: This results in A^TA being a singular matrix
                  num_rows_A = 7;
                  cout << "===== boundary node: " << node << "\n";
                  continue;
               }
               case 4:
               {
                  num_rows_A = 12;
                  break;
               }
               default:
               {
                  MFEM_ABORT("Invalid number of adjacent ceells.\n");
                  break;
               }
            }
            A.SetSize(num_rows_A, num_cols_A);
            RHS.SetSize(num_rows_A);
            temp_vec.SetSize(num_rows_A);

            // cout << "===== new corner node ===== [at: " << node << "]\n";
            /*** Fill matrix representing overdetermined system and RHS vector ***/
            // Iterate over adjacent cells
            for (int element_it = 0; element_it < length_element_row; element_it++)
            {
               int el_index = element_row[element_it];

               pmesh->GetElementEdges(el_index, cell_faces_row, oris);
               length_cell_faces = cell_faces_row.Size();

               // Iterate over cell faces
               for (int face_it = 0; face_it < length_cell_faces; face_it++)
               {
                  int face_index = cell_faces_row[face_it];

                  // Check if cell face is unique
                  if (uniqueFaceIndices.find(face_index) == uniqueFaceIndices.end())
                  {
                     // New face, fill corresponding row in A and RHS
                     uniqueFaceIndices.insert(face_index);
                     row_index = uniqueFaceIndices.size() - 1;
                     // cout << "row_index: " << row_index << endl;
                     
                     // Get face velocity, outward normal, and face position
                     // n_vec is the outward normal times |F|
                     GetIntermediateFaceVelocity(face_index, Vf);
                     CalcOutwardNormalInt(S, el_index, face_index, n_vec);

                     double indicator = Vf * n_vec;
                     // if (indicator < 0.)
                     // {
                     //    n_vec *= -1.;
                     // }

                     // Get average of face edge nodes
                     H1.GetFaceDofs(face_index, face_dofs);
                     geom.GetNodePositionFromBV(S, face_dofs[1], vdof1_x);
                     geom.GetNodePositionFromBV(S, face_dofs[0], vdof2_x);
                     face_x = vdof1_x;
                     face_x.Add(1., vdof2_x);
                     face_x *= 0.5;
                     
                     // Cout info
                     // cout << "Unique face corresponds to el: " << el_index << ", face: " << face_index << endl;
                     // cout << "Vf: ";
                     // Vf.Print(cout);
                     // cout << "nf: ";
                     // n_vec.Print(cout);
                     // cout << "face_x: ";
                     // face_x.Print(cout);
                     // geom.GetNodePositionFromBV(S, face_dofs[2], face_x);

                     // Fill rows of A
                     A(row_index, 0) = n_vec(0);
                     A(row_index, 1) = n_vec(1);
                     A(row_index, 2) = n_vec(0)*face_x(0);
                     A(row_index, 3) = n_vec(1)*face_x(0);
                     A(row_index, 4) = n_vec(0)*face_x(1);
                     A(row_index, 5) = n_vec(1)*face_x(1);

                     // Fill rows of RHS
                     RHS(row_index) = Vf * n_vec;
                     // cout << endl;
                  } // End unique new face
               } // End face iterator
            } // End cell iterator

            // cout << "RHS: ";
            // RHS.Print(cout);

            // Solve for vgeo on this node
            AT = A;
            AT.Transpose();
            mfem::Mult(AT, A, ATA);
            AT.Mult(RHS, RHS_mod);
            // cout << "A: " << endl;
            // A.Print(cout);
            // cout << "ATA: " << endl;
            // ATA.Print(cout);
            ATA.Invert();
            ATA.Mult(RHS_mod,res);

            // Compute residual
            A.Mult(res, temp_vec);
            // cout << "temp_vec: ";
            // temp_vec.Print(cout);

            temp_vec -= RHS;
            // cout << "residual: " << temp_vec.Norml2() << endl;

            // cout << "coefficient vector: ";
            // res.Print(cout);

            node_v(0) = res[0] + res[2]*node_x[0] + res[4]*node_x[1];
            node_v(1) = res[1] + res[3]*node_x[0] + res[5]*node_x[1];

            // cout << "node_x: ";
            // node_x.Print(cout);
            // cout << "node_v: ";
            // node_v.Print(cout);

            // Only update the v_geo_gf in the case of the corners
            // Using Q1 for this velocity field
            SetViGeo(node, node_v);

            break;
         } // End corner
         case 1: // face
         {
            // face nodes just set to ivf
            assert(face >= 0);
            GetIntermediateFaceVelocity(EDof, node_v);

            break;
         }
         case 2: // Cell Center
         {
            Vector Uc(dim+2);
            GetCellStateVector(S, EDof, Uc);
            pb->velocity(Uc, node_v);

            break;
         }
         default:
         {
            MFEM_ABORT("Invalid entity\n");
         }
      }
      // if (node_v.Norml2() > 0.0001)
      // {
      //    cout << "node: " << node << endl       
      //         << "nodev: ";
      //    node_v.Print(cout);
      // }

      // In every case, update the corresponding nodal velocity in S
      geom.UpdateNodeVelocity(S, node, node_v);

   } // End node iterator
}


/****************************************************************************************************
* Function: ComputeGeoVCAVEAT
* Parameters:
*  S    - BlockVector corresponding to nth timestep that stores mesh information, 
*         mesh velocity, and state variables.
*
* Purpose:
*  This function computes the geometric velocity from the face normals and face 
*  velocities from all adjacent faces following the Caveat weighted least squares
*  algorithm described in https://www.osti.gov/biblio/5366795.
*
*  The weights used for the weighted least squares method is simply the sum of 
*  the densities across any given face.
*
*
*  Note: Vigeo is only the geometric velocity, and is not the velocity used 
*        to move the mesh.  For the mesh velocity, one must solve for the 
*        alpha_i, which process is found in ComputeNodeVelocitiesFromVgeo().
*
*  Note: This function fills the Q1 velocity field v_geo_gf and also the mv_gf
*        in the BlockVector S.
****************************************************************************************************/
void LagrangianLOOperator::ComputeGeoVCAVEAT(Vector &S)
{
   // cout << "==========ComputeGeoVCaveat==========\n";
   mfem::Mesh::FaceInformation FI;
   Vector node_v(dim), Ri(dim), n_vec(dim), Vf(dim), Uc(dim+2), Ucp(dim+2);
   DenseMatrix A, AT, WA, ATWA(dim,dim), W;
   Vector RHS, RHS_mod;

   Array<int> faces_row, oris;
   int faces_length, c, cp;
   double c_norm, weight;

   DofEntity entity;
   int EDof;

   // Iterate over all mesh velocity dofs
   for (int node = 0; node < NDofs_H1; node++) // Vertex iterator
   {
      // cout << "\n\t $$$node: " << node << endl;
      node_v = 0.; // Reset node velocity
      GetEntityDof(node, entity, EDof);

      switch (entity)
      {
         case 0: // corner
         {
            node_v = 0.;
            // Get cell faces
            vertex_edge.GetRow(node, faces_row);
            faces_length = faces_row.Size();

            // Set sizes of nodal matrices
            A.SetSize(faces_length, dim);
            AT.SetSize(dim, faces_length);
            W.SetSize(faces_length);
            WA.SetSize(faces_length, dim);
            RHS.SetSize(faces_length);
            RHS_mod.SetSize(dim);
            A = 0., AT = 0., ATWA = 0., W = 0., WA = 0., RHS = 0., RHS_mod = 0.;

            // iterate over adjacent faces
            for (int face_it = 0; face_it < faces_length; face_it++) // Adjacent face iterator
            {
               int face = faces_row[face_it];
               GetIntermediateFaceVelocity(face, Vf);
               FI = pmesh->GetFaceInformation(face);
               c = FI.element[0].index;
               cp = FI.element[1].index;

               // Retrieve corresponding normal (orientation doesn't matter)
               CalcOutwardNormalInt(S, FI.element[0].index, face, n_vec);
               c_norm = n_vec.Norml2();
               n_vec /= c_norm;

               // Compute weight
               double rhoc = 0., rhocp = 0.;
               GetCellStateVector(S, c, Uc);
               
               rhoc = 1. / Uc[0];
               if (FI.IsInterior())
               {
                  GetCellStateVector(S, cp, Ucp);
                  rhocp = 1 / Ucp[0];
               }
               weight = rhoc + rhocp;

               // Add in contribution to system
               A.SetRow(face_it, n_vec);
               W(face_it, face_it) = weight;
               double val = n_vec * Vf;
               RHS[face_it] = weight * val;

               // Print stuff
               // cout << "face: " << face << ", vf: ";
               // Vf.Print(cout);
               // cout << "normal: ";
               // n_vec.Print(cout);
               // cout << "weight: " << weight << endl;
            }
            // Matrix setup
            AT = A;
            AT.Transpose();
            mfem::Mult(W, A, WA);
            // cout << "WA: ";
            // WA.Print(cout);
            mfem::Mult(AT,WA,ATWA);
            // cout << "ATWA:\n";
            // ATWA.Print(cout);
            AT.Mult(RHS, RHS_mod);

            // Solve for best fit
            ATWA.Invert();
            ATWA.Mult(RHS_mod, node_v);

             // Print stuff
            // cout << "A:\n";
            // A.Print(cout);
            // cout << "AT:\n";
            // AT.Print(cout);
            // cout << "W:\n";
            // W.Print(cout);
            // cout << "RHS: ";
            // RHS.Print(cout);
            // cout << "RHS_mod: ";
            // RHS_mod.Print(cout);
            // cout << "node_v: ";
            // node_v.Print(cout); 

            // Only update the v_geo_gf in the case of the corners
            // Using Q1 for this velocity field
            SetViGeo(node, node_v);

            break;
         }
         case 1: // face
         {
            // face nodes just set to ivf
            assert(face >= 0);
            GetIntermediateFaceVelocity(EDof, node_v);

            break;
         }
         case 2: // Cell Center
         {
            Vector Uc(dim+2);
            GetCellStateVector(S, EDof, Uc);
            pb->velocity(Uc, node_v);

            break;
         }
         default:
         {
            MFEM_ABORT("Invalid entity\n");
         }
      }

      // In every case, update the corresponding nodal velocity in S
      geom.UpdateNodeVelocity(S, node, node_v);

   } // End node iterator
}


/****************************************************************************************************
* Function: ComputeGeoVCAVEATCellFace
* Parameters:
*  S    - BlockVector corresponding to nth timestep that stores mesh information, 
*         mesh velocity, and state variables.
*
* Purpose:
*  This function computes the geometric velocity from the face normals and face 
*  velocities from all adjacent faces.
*
*  On the boundary vertices, the above CAVEAT algorithm is used, which is 
*  described in https://www.osti.gov/biblio/5366795.
*
*  On the interior vertices, the above adjacent Cell face based method is used.
*
*  Note: Vigeo is only the geometric velocity, and is not the velocity used 
*        to move the mesh.  For the mesh velocity, one must solve for the 
*        alpha_i, which process is found in ComputeNodeVelocitiesFromVgeo().
*
*  Note: This function fills the Q1 velocity field v_geo_gf and also the mv_gf
*        in the BlockVector S.
****************************************************************************************************/
void LagrangianLOOperator::ComputeGeoVCAVEATCellFace(Vector &S)
{
   // cout << "ComputeGeoVCAVEATCellFace call\n";
   mfem::Mesh::FaceInformation FI;
   
   Vector node_v(dim), node_x(dim), n_vec(dim), Vf(dim);
   Vector vdof1_x(dim), vdof2_x(dim), face_x(dim), vdof_avg(dim);
   double F, c_norm;

   Array<int> element_row, cell_faces_row, oris, face_dofs;
   int row_index, length_element_row, length_cell_faces, num_rows_A, num_cols_A = 6;
   int c, cp;
   Vector Uc(dim + 2), Ucp(dim+2);

   DenseMatrix A, AT, ATA(num_cols_A), W, WA, ATWA(dim);
   Vector RHS, RHS_mod(num_cols_A), res(num_cols_A), temp_vec;

   DofEntity entity;
   int EDof;

   std::unordered_set<int> uniqueFaceIndices;

   // Iterate over all mesh velocity dofs
   for (int node = 0; node < NDofs_H1; node++) // Vertex iterator
   {
      node_v = 0.; // Reset node velocity
      GetEntityDof(node, entity, EDof);

      switch (entity)
      {
         case 0: // corner
         {
            // Reset vals
            uniqueFaceIndices.clear();
            F = 0.;
            node_v = 0.;
            geom.GetNodePositionFromBV(S, node, node_x);

            // Get adjacent cells
            vertex_element->GetRow(node, element_row);
            length_element_row = element_row.Size();

            // Set size of overdetermined matrix A depending on number of adjacent cells
            switch (length_element_row)
            {
               case 1: // Domain corner vertex
               case 2: // Boundary vertex
               {
                  // Since in the Cell Face based method we get a singular matrix for
                  // the domain corner and boundary vertices, we use the CAVEAT method
                  // of projecting the vertex velocity onto the adjacent cell faces
                  // and solve the system using weighted least squares.
                  if (length_element_row == 1) {
                     num_rows_A = 4;
                  } else {
                     num_rows_A = 7;
                  }

                  // Set sizes of nodal matrices
                  A.SetSize(num_rows_A, dim);
                  AT.SetSize(dim, num_rows_A);
                  W.SetSize(num_rows_A);
                  WA.SetSize(num_rows_A, dim);
                  RHS.SetSize(num_rows_A);
                  RHS_mod.SetSize(dim);
                  A = 0., AT = 0., ATWA = 0., W = 0., WA = 0., RHS = 0., RHS_mod = 0.;

                  /*** Fill matrix representing overdetermined system and RHS vector ***/
                  // Iterate over adjacent cells
                  for (int element_it = 0; element_it < length_element_row; element_it++)
                  {
                     int el_index = element_row[element_it];

                     pmesh->GetElementEdges(el_index, cell_faces_row, oris);
                     length_cell_faces = cell_faces_row.Size();

                     // Iterate over cell faces
                     for (int face_it = 0; face_it < length_cell_faces; face_it++)
                     {
                        int face_index = cell_faces_row[face_it];

                        // Check if cell face is unique
                        if (uniqueFaceIndices.find(face_index) == uniqueFaceIndices.end())
                        {
                           // New face, fill corresponding row in A and RHS
                           uniqueFaceIndices.insert(face_index);
                           row_index = uniqueFaceIndices.size() - 1;

                           // Get adjacent cell information for the weighting
                           FI = pmesh->GetFaceInformation(face_index);
                           c = FI.element[0].index;
                           cp = FI.element[1].index;
                           
                           // Get face velocity, outward normal, and face position
                           // n_vec is the outward normal times |F|
                           GetIntermediateFaceVelocity(face_index, Vf);
                           CalcOutwardNormalInt(S, el_index, face_index, n_vec);

                           // Compute weight
                           double rhoc = 0., rhocp = 0.;
                           GetCellStateVector(S, c, Uc);
                           
                           // rhoc = 1. / Uc[0];
                           // if (FI.IsInterior())
                           // {
                           //    GetCellStateVector(S, cp, Ucp);
                           //    rhocp = 1 / Ucp[0];
                           // }
                           // double weight = rhoc + rhocp;
                           double weight = 1.;

                           // Fill rows of A
                           A.SetRow(row_index, n_vec);
                           // Fill weight matrix
                           W(row_index, row_index) = weight;
                           // Fill rows of RHS
                           double val = n_vec * Vf;
                           RHS(row_index) = weight * val;
                           // cout << endl;
                        } // End unique new face
                     } // End face iterator
                  } // End cell iterator

                  // Matrix setup
                  AT = A;
                  AT.Transpose();
                  mfem::Mult(W, A, WA);
                  mfem::Mult(AT,WA,ATWA);
                  AT.Mult(RHS, RHS_mod);

                  // Solve for best fit
                  ATWA.Invert();
                  ATWA.Mult(RHS_mod, node_v);

                  // Compute residual
                  // temp_vec.SetSize(num_rows_A);
                  // A.Mult(node_v, temp_vec);
                  // temp_vec -= RHS;
                  // cout << "Boundary/Corner Node: " << node << ", residual: " << temp_vec.Norml2() << endl;

                  break;
               } // End corner and boundary vertex
               case 4:
               {
                  num_rows_A = 12;

                  A.SetSize(num_rows_A, num_cols_A);
                  RHS.SetSize(num_rows_A);
                  temp_vec.SetSize(num_rows_A);

                  // cout << "===== new corner node ===== [at: " << node << "]\n";
                  /*** Fill matrix representing overdetermined system and RHS vector ***/
                  // Iterate over adjacent cells
                  for (int element_it = 0; element_it < length_element_row; element_it++)
                  {
                     int el_index = element_row[element_it];

                     pmesh->GetElementEdges(el_index, cell_faces_row, oris);
                     length_cell_faces = cell_faces_row.Size();

                     // Iterate over cell faces
                     for (int face_it = 0; face_it < length_cell_faces; face_it++)
                     {
                        int face_index = cell_faces_row[face_it];

                        // Check if cell face is unique
                        if (uniqueFaceIndices.find(face_index) == uniqueFaceIndices.end())
                        {
                           // New face, fill corresponding row in A and RHS
                           uniqueFaceIndices.insert(face_index);
                           row_index = uniqueFaceIndices.size() - 1;
                           
                           // Get face velocity, outward normal, and face position
                           // n_vec is the outward normal times |F|
                           GetIntermediateFaceVelocity(face_index, Vf);
                           CalcOutwardNormalInt(S, el_index, face_index, n_vec);

                           // Get average of face edge nodes
                           H1.GetFaceDofs(face_index, face_dofs);
                           geom.GetNodePositionFromBV(S, face_dofs[1], vdof1_x);
                           geom.GetNodePositionFromBV(S, face_dofs[0], vdof2_x);
                           face_x = vdof1_x;
                           face_x.Add(1., vdof2_x);
                           face_x *= 0.5;

                           // Fill rows of A
                           A(row_index, 0) = n_vec(0);
                           A(row_index, 1) = n_vec(1);
                           A(row_index, 2) = n_vec(0)*face_x(0);
                           A(row_index, 3) = n_vec(1)*face_x(0);
                           A(row_index, 4) = n_vec(0)*face_x(1);
                           A(row_index, 5) = n_vec(1)*face_x(1);

                           // Fill rows of RHS
                           RHS(row_index) = Vf * n_vec;
                        } // End unique new face
                     } // End face iterator
                  } // End cell iterator

                  // Solve for vgeo on this node
                  AT = A;
                  AT.Transpose();
                  mfem::Mult(AT, A, ATA);
                  AT.Mult(RHS, RHS_mod);
                  ATA.Invert();
                  ATA.Mult(RHS_mod,res);

                  // Compute residual
                  // A.Mult(res, temp_vec);
                  // temp_vec -= RHS;
                  // cout << "Interior Node: " << node << ", residual: " << temp_vec.Norml2() << endl;

                  // Set node velocity from coefficient vector
                  node_v(0) = res[0] + res[2]*node_x[0] + res[4]*node_x[1];
                  node_v(1) = res[1] + res[3]*node_x[0] + res[5]*node_x[1];
                  break;
               }
               default:
               {
                  MFEM_ABORT("Invalid number of adjacent ceells.\n");
                  break;
               }
            } // end switch between corner, boundary, and interior vertex

            // Only update the v_geo_gf in the case of the corners
            // Using Q1 for this velocity field
            SetViGeo(node, node_v);

            break;
         } // End corner
         case 1: // face
         {
            // face nodes just set to ivf
            assert(face >= 0);
            GetIntermediateFaceVelocity(EDof, node_v);

            break;
         } // end face
         case 2: // Cell Center
         {
            Vector Uc(dim+2);
            GetCellStateVector(S, EDof, Uc);
            pb->velocity(Uc, node_v);

            break;
         } // End cell center
         default:
         {
            MFEM_ABORT("Invalid entity\n");
         }
      } // end entity switch

      // In every case, update the corresponding nodal velocity in S
      geom.UpdateNodeVelocity(S, node, node_v);

   } // End node iterator
}


/****************************************************************************************************
* Function: ComputeGeoVCAVEATCellFaceWeighted
* Parameters:
*  S    - BlockVector corresponding to nth timestep that stores mesh information, 
*         mesh velocity, and state variables.
*
* Purpose:
*  This function computes the geometric velocity from the face normals and face 
*  velocities from all adjacent faces.
*
*  On the boundary vertices, the above CAVEAT algorithm is used, which is 
*  described in https://www.osti.gov/biblio/5366795.
*
*  On the interior vertices, the above adjacent Cell face based method is used.
*
*  Note: Vigeo is only the geometric velocity, and is not the velocity used 
*        to move the mesh.  For the mesh velocity, one must solve for the 
*        alpha_i, which process is found in ComputeNodeVelocitiesFromVgeo().
*
*  Note: This function fills the Q1 velocity field v_geo_gf and also the mv_gf
*        in the BlockVector S.
****************************************************************************************************/
void LagrangianLOOperator::ComputeGeoVCAVEATCellFaceWeighted(Vector &S)
{
   // cout << "ComputeGeoVCAVEATCellFaceWeighted call\n";
   mfem::Mesh::FaceInformation FI;
   
   Vector node_v(dim), node_x(dim), n_vec(dim), Vf(dim);
   Vector vdof1_x(dim), vdof2_x(dim), face_x(dim), vdof_avg(dim);
   double F, c_norm;

   Array<int> element_row, cell_faces_row, oris, face_dofs;
   int row_index, length_element_row, length_cell_faces, num_rows_A, num_cols_A = 6;
   int c, cp;
   Vector Uc(dim + 2), Ucp(dim+2);

   DenseMatrix A, AT, ATA(num_cols_A), W, WA, ATWA(dim);
   Vector RHS, RHS_mod(num_cols_A), res(num_cols_A), temp_vec;

   DofEntity entity;
   int EDof;

   std::unordered_set<int> uniqueFaceIndices;

   // Iterate over all mesh velocity dofs
   for (int node = 0; node < NDofs_H1; node++) // Vertex iterator
   {
      node_v = 0.; // Reset node velocity
      GetEntityDof(node, entity, EDof);

      switch (entity)
      {
         case 0: // corner
         {
            // Reset vals
            uniqueFaceIndices.clear();
            F = 0.;
            node_v = 0.;
            geom.GetNodePositionFromBV(S, node, node_x);

            // Get adjacent cells
            vertex_element->GetRow(node, element_row);
            length_element_row = element_row.Size();

            // Set size of overdetermined matrix A depending on number of adjacent cells
            switch (length_element_row)
            {
               case 1: // Domain corner vertex
               case 2: // Boundary vertex
               {
                  // Since in the Cell Face based method we get a singular matrix for
                  // the domain corner and boundary vertices, we use the CAVEAT method
                  // of projecting the vertex velocity onto the adjacent cell faces
                  // and solve the system using weighted least squares.
                  if (length_element_row == 1) {
                     num_rows_A = 4;
                  } else {
                     num_rows_A = 7;
                  }

                  // Set sizes of nodal matrices
                  A.SetSize(num_rows_A, dim);
                  AT.SetSize(dim, num_rows_A);
                  W.SetSize(num_rows_A);
                  WA.SetSize(num_rows_A, dim);
                  RHS.SetSize(num_rows_A);
                  RHS_mod.SetSize(dim);
                  A = 0., AT = 0., ATWA = 0., W = 0., WA = 0., RHS = 0., RHS_mod = 0.;

                  /*** Fill matrix representing overdetermined system and RHS vector ***/
                  // Iterate over adjacent cells
                  for (int element_it = 0; element_it < length_element_row; element_it++)
                  {
                     int el_index = element_row[element_it];

                     pmesh->GetElementEdges(el_index, cell_faces_row, oris);
                     length_cell_faces = cell_faces_row.Size();

                     // Iterate over cell faces
                     for (int face_it = 0; face_it < length_cell_faces; face_it++)
                     {
                        int face_index = cell_faces_row[face_it];

                        // Check if cell face is unique
                        if (uniqueFaceIndices.find(face_index) == uniqueFaceIndices.end())
                        {
                           // New face, fill corresponding row in A and RHS
                           uniqueFaceIndices.insert(face_index);
                           row_index = uniqueFaceIndices.size() - 1;

                           // Get adjacent cell information for the weighting
                           FI = pmesh->GetFaceInformation(face_index);
                           c = FI.element[0].index;
                           cp = FI.element[1].index;
                           
                           // Get face velocity, outward normal, and face position
                           // n_vec is the outward normal times |F|
                           GetIntermediateFaceVelocity(face_index, Vf);
                           CalcOutwardNormalInt(S, el_index, face_index, n_vec);

                           // Compute weight
                           double rhoc = 0., rhocp = 0.;
                           GetCellStateVector(S, c, Uc);
                           
                           rhoc = 1. / Uc[0];
                           if (FI.IsInterior())
                           {
                              GetCellStateVector(S, cp, Ucp);
                              rhocp = 1. / Ucp[0];
                           }
                           double weight = rhoc + rhocp;
                           // double weight = 1.;

                           // Fill rows of A
                           A.SetRow(row_index, n_vec);
                           // Fill weight matrix
                           W(row_index, row_index) = weight;
                           // Fill rows of RHS
                           double val = n_vec * Vf;
                           RHS(row_index) = weight * val;
                           // cout << endl;

                           // Print stuff
                           // cout << "face: " << face << ", vf: ";
                           // Vf.Print(cout);
                           // cout << "normal: ";
                           // n_vec.Print(cout);
                           // cout << "weight: " << weight << endl;
                        } // End unique new face
                     } // End face iterator
                  } // End cell iterator

                  // Matrix setup
                  AT = A;
                  AT.Transpose();
                  mfem::Mult(W, A, WA);
                  mfem::Mult(AT,WA,ATWA);
                  AT.Mult(RHS, RHS_mod);

                  // Solve for best fit
                  ATWA.Invert();
                  ATWA.Mult(RHS_mod, node_v);

                  // Compute residual
                  // temp_vec.SetSize(num_rows_A);
                  // A.Mult(node_v, temp_vec);
                  // temp_vec -= RHS;
                  // cout << "Boundary/Corner Node: " << node << ", residual: " << temp_vec.Norml2() << endl;

                  break;
               } // End corner and boundary vertex
               case 4:
               {
                  num_rows_A = 12;

                  A.SetSize(num_rows_A, num_cols_A);
                  RHS.SetSize(num_rows_A);
                  temp_vec.SetSize(num_rows_A);

                  // cout << "===== new corner node ===== [at: " << node << "]\n";
                  /*** Fill matrix representing overdetermined system and RHS vector ***/
                  // Iterate over adjacent cells
                  for (int element_it = 0; element_it < length_element_row; element_it++)
                  {
                     int el_index = element_row[element_it];

                     pmesh->GetElementEdges(el_index, cell_faces_row, oris);
                     length_cell_faces = cell_faces_row.Size();

                     // Iterate over cell faces
                     for (int face_it = 0; face_it < length_cell_faces; face_it++)
                     {
                        int face_index = cell_faces_row[face_it];

                        // Check if cell face is unique
                        if (uniqueFaceIndices.find(face_index) == uniqueFaceIndices.end())
                        {
                           // New face, fill corresponding row in A and RHS
                           uniqueFaceIndices.insert(face_index);
                           row_index = uniqueFaceIndices.size() - 1;
                           
                           // Get face velocity, outward normal, and face position
                           // n_vec is the outward normal times |F|
                           GetIntermediateFaceVelocity(face_index, Vf);
                           CalcOutwardNormalInt(S, el_index, face_index, n_vec);

                           // Compute weight
                           double rhoc = 0., rhocp = 0., weight = 0.;
                           GetCellStateVector(S, c, Uc);
                           
                           rhoc = 1. / Uc[0];
                           if (FI.IsInterior())
                           {
                              GetCellStateVector(S, cp, Ucp);
                              rhocp = 1. / Ucp[0];
                              weight = rhoc + rhocp;
                           }
                           else 
                           {
                              weight = rhoc;
                           }
                           
                           // double weight = 1.;

                           // Get average of face edge nodes
                           H1.GetFaceDofs(face_index, face_dofs);
                           geom.GetNodePositionFromBV(S, face_dofs[1], vdof1_x);
                           geom.GetNodePositionFromBV(S, face_dofs[0], vdof2_x);
                           face_x = vdof1_x;
                           face_x.Add(1., vdof2_x);
                           face_x *= 0.5;

                           // Fill rows of A
                           A(row_index, 0) = n_vec(0);
                           A(row_index, 1) = n_vec(1);
                           A(row_index, 2) = n_vec(0)*face_x(0);
                           A(row_index, 3) = n_vec(1)*face_x(0);
                           A(row_index, 4) = n_vec(0)*face_x(1);
                           A(row_index, 5) = n_vec(1)*face_x(1);

                           // Fill weight matrix
                           W(row_index, row_index) = weight;

                           // Fill rows of RHS
                           double rhs_val = Vf * n_vec;
                           rhs_val *= weight;
                           RHS(row_index) = rhs_val;

                           // Print stuff
                           // cout << "face: " << face << ", vf: ";
                           // Vf.Print(cout);
                           // cout << "normal: ";
                           // n_vec.Print(cout);
                           // cout << "weight: " << weight << endl;
                        } // End unique new face
                     } // End face iterator
                  } // End cell iterator

                  // Matrix setup
                  AT = A;
                  AT.Transpose();
                  mfem::Mult(W, A, WA);
                  mfem::Mult(AT,WA,ATWA);
                  AT.Mult(RHS, RHS_mod);

                  // Solve for best fit
                  ATWA.Invert();
                  ATWA.Mult(RHS_mod, res);

                  // Compute residual
                  // A.Mult(res, temp_vec);
                  // temp_vec -= RHS;
                  // cout << "Interior Node: " << node << ", residual: " << temp_vec.Norml2() << endl;

                  // Set node velocity from coefficient vector
                  node_v(0) = res[0] + res[2]*node_x[0] + res[4]*node_x[1];
                  node_v(1) = res[1] + res[3]*node_x[0] + res[5]*node_x[1];
                  break;
               }
               default:
               {
                  MFEM_ABORT("Invalid number of adjacent ceells.\n");
                  break;
               }
            } // end switch between corner, boundary, and interior vertex

            // Only update the v_geo_gf in the case of the corners
            // Using Q1 for this velocity field
            SetViGeo(node, node_v);

            break;
         } // End corner
         case 1: // face
         {
            // face nodes just set to ivf
            assert(face >= 0);
            GetIntermediateFaceVelocity(EDof, node_v);

            break;
         } // end face
         case 2: // Cell Center
         {
            Vector Uc(dim+2);
            GetCellStateVector(S, EDof, Uc);
            pb->velocity(Uc, node_v);

            break;
         } // End cell center
         default:
         {
            MFEM_ABORT("Invalid entity\n");
         }
      } // end entity switch

      // In every case, update the corresponding nodal velocity in S
      geom.UpdateNodeVelocity(S, node, node_v);

   } // End node iterator
}


/****************************************************************************************************
* Function: IterativeCornerVelocityMC
* Parameters:
*  S    - BlockVector corresponding to nth timestep that stores mesh information, 
*         mesh velocity, and state variables.
*  dt   - Current timestep size
*
* Purpose:
*  This function constitutes one iteration on the previously computed corner node 
*  velocity to reduce total face motion while still preserving mass conservation.
*  The idea is that we correct the velocity at the vertices to ensure mass 
*  conservation and dampen the parachuting effect seen at the faces. This method 
*  was proposed on 04/03/2024.
*
*  Note: This function relies on the mesh velocities having already been computed
*  and linearized.  One can use whichever mesh velocity computation before 
*  iteration.
*
*  Note: This function will modify mv_gf in the BlockVector S
****************************************************************************************************/
void LagrangianLOOperator::IterativeCornerVelocityMC(Vector &S, const double & dt)
{
   // cout << "=====IterativeCornerVelocityMC=====\n";
   // Optional run time parameters
   bool do_theta_averaging = false;
   double theta = 1.;

   // Values needed during iteration
   Array<int> faces_row, face_dofs_row;
   int faces_length;
   mfem::Mesh::FaceInformation FI;
   H1.ExchangeFaceNbrData();
   int Vadj_index, c, bdr_ind = 0;
   double F = 0.;
   double Vnode_n_comp = 0.; // The quantity we solve for in each face iteration
   double Vnode_nR_comp = 0.;
   double Vnode_prev_it_nR_comp = 0., Vadj_n_comp = 0., Vadj_nR_comp = 0.;

   Vector predicted_node_v(dim), Vf(dim), Vnode_np1(dim), Vnode_prev_it(dim);
   Vector n_int(dim), n_vec(dim), n_vec_R(dim);
   Vector Vnode_x(dim), Vadj_x(dim), face_x(dim);
   Vector Vnode_xnp1(dim), Vadj_xnp1(dim), face_xnp1(dim), a12_xnp1(dim);
   Vector Badj(dim), Bnode(dim);

   Vector Vnode_n(dim), Vadj_n(dim), Vface_n(dim);
   Vector Vnode_half(dim), Vadj_half(dim), n_vec_half(dim);
   Vector temp_vec(dim), temp_vec_2(dim);
   double D = 0., c1 = 0., c0 = 0., V3nperp = 0., Dc1 = 0.;
   int num_faces_for_avg;

   // Averaging
   DenseMatrix Mi(dim), dm_tmp(dim);
   Vector Ri(dim), v_tmp(dim);

   /* Iterate over corner nodes */
   for (int node = 0; node < NVDofs_H1; node++) // TODO: Is NDofs_H1L == NVDofs_H1?
   {
      // Reset new nodal velocity, and averaging objects
      Mi = 0., Ri = 0.;
      predicted_node_v = 0.;
      num_faces_for_avg = 0;

      // Set corresponding boundary indicator
      // If node is boundary node, will have to add corresponding correction
      // node is interior node if bdr_ind == -1
      bdr_ind = BdrVertexIndexingArray[node];

      // Get current nodal velocity and position 
      geom.GetNodeVelocity(S, node, Vnode_n);
      geom.GetNodePositionFromBV(S, node, Vnode_x);

      // Compute An and an+1 for node
      add(Vnode_x, dt/2., Vnode_n, Vnode_half);
      add(Vnode_x, dt, Vnode_n, Vnode_xnp1);

      // Get cell faces
      vertex_edge.GetRow(node, faces_row);
      faces_length = faces_row.Size();

      /* Iterate over cell faces */
      for (int face_it = 0; face_it < faces_length; face_it++) // Adjacent face iterator
      {
         // Increment averaging denom 
         num_faces_for_avg += 1;
         // Get face information
         int face = faces_row[face_it];
         GetIntermediateFaceVelocity(face, Vf);
         FI = pmesh->GetFaceInformation(face);

         // Calculate outer normal
         c = FI.element[0].index;
         CalcOutwardNormalInt(S, c, face, n_int);
         n_vec = n_int;
         F = n_vec.Norml2();
         n_vec /= F;

         /* adjacent corner indices */
         // preserve node orientation where cell to right 
         // of face is with lower cell index
         H1.GetFaceDofs(face, face_dofs_row);
         int face_vdof1 = face_dofs_row[1], 
             face_vdof2 = face_dofs_row[0], 
             face_dof = face_dofs_row[2]; 

         // Grab corresponding vertex velocity from S
         // We are always solving for V2. 
         // If the index for Vnode does not match that
         // for V2, we must flip the normal.
         if (node == face_vdof1) {
            // Get adj index
            Vadj_index = face_vdof2;
            // Flip normal vector
            n_vec *= -1.;
         } else {
            // Get adj index
            Vadj_index = face_vdof1;
            // Node matches v2
         }

         /* Get face and adj node velocity and position */
         geom.GetNodeVelocity(S, Vadj_index, Vadj_n);
         geom.GetNodeVelocity(S, face_dof, Vface_n);
         geom.GetNodePositionFromBV(S, Vadj_index, Vadj_x);
         geom.GetNodePositionFromBV(S, face_dof, face_x);

         // Check to make sure orientation is correct
         // This part is mostly for assertion
         Vector ttvec(dim);
         subtract(Vnode_x, Vadj_x, ttvec);
         geom.Orthogonal(ttvec);
         double _val = ttvec * n_vec;
         // Vector should be pointing in same direction
         if (_val < 1E-6)
         {
            cout << "interior analysis:\n";
            cout << "node: " << node << endl;
            cout << "nodex: ";
            Vnode_x.Print();
            cout << "Vadj_x: ";
            Vadj_x.Print(cout);
            cout << "nvec: ";
            n_vec.Print(cout);
            assert(false);
         }

         // Compute future face and adjacent node locations
         add(Vadj_x, dt/2., Vadj_n, Vadj_half);
         add(Vadj_x, dt, Vadj_n, Vadj_xnp1);
         add(0.5, Vadj_xnp1, 0.5, Vnode_xnp1, a12_xnp1);

         // Perpendicular vector
         n_vec_R = n_vec;
         geom.Orthogonal(n_vec_R);

         // Calculate bmn
         double bmn = Vf * n_vec;
         bmn *= F;

         // Get Vnode_prev_it component in tangent direction from previous iteration
         Vnode_prev_it_nR_comp = Vnode_n * n_vec_R;

         // Get normal and rotated components of Vadj_n
         Vadj_n_comp = Vadj_n * n_vec;
         Vadj_nR_comp = Vadj_n * n_vec_R;
         
         // Compute geometrical vectors
         Bnode = face_x;
         Bnode *= -2.;
         Bnode.Add(-.5, Vadj_x);
         Bnode.Add(2.5, Vnode_x);
         geom.Orthogonal(Bnode);

         Badj = face_x;
         Badj *= -2.;
         Badj.Add(-.5, Vnode_x);
         Badj.Add(2.5, Vadj_x);
         geom.Orthogonal(Badj);

         // Compute dot products
         double badjn = Badj * n_vec;
         double badjnr = Badj * n_vec_R;
         double bnoden = Bnode * n_vec;
         double bnodenr = Bnode * n_vec_R;

         // Evaluate numerator, which depends on geometric orientation
         // The center of the cell should be on the right when 
         // traversing from adj node to node.
         double numer = 3. * bmn;

         // Add in Dc1V3nper contribution (explicit)
         // V3nperp is set to 0 in option 1.
         if (mv_it_option == 2 || mv_it_option == 3)
         {
            // Calculate D and c1 (perturbations only need to be handled explicitly)
            subtract(Vadj_n, Vnode_n, temp_vec); // V1 - V2 = temp_vec
            subtract(Vnode_half, Vadj_half, temp_vec_2); // A2-A1

            geom.Orthogonal(temp_vec_2);
            D = dt * (temp_vec * n_vec_R) + 2. * (n_vec * temp_vec_2);

            // Compute c1 (A.4a)
            subtract(Vnode_n, Vadj_n, temp_vec); // only change temp_vec, since temp_vec_2 is same from D calculation (half step representation)
            Dc1 = dt * (temp_vec * n_vec) + 2. * (temp_vec_2 * n_vec_R); 

            /* Compute V3nperp using previous iteration */
            add(0.5, Vnode_n, 0.5, Vadj_n, temp_vec);
            V3nperp = -1. * (temp_vec * n_vec_R);
            numer += Dc1*V3nperp;
         }
         numer += (badjn - dt * Vadj_nR_comp) * Vadj_n_comp;
         numer += badjnr * Vadj_nR_comp;
         numer += ((dt/2) * Vadj_n_comp - bnodenr)*Vnode_prev_it_nR_comp;

         double denom = -dt*Vnode_prev_it_nR_comp + (dt/2)*Vadj_nR_comp + bnoden;
   
         Vnode_n_comp = numer / denom;

         // Add in contribution to predicted nodal velocity for this face
         Vnode_np1 = 0.;
         Vnode_np1.Add(Vnode_n_comp, n_vec);
         Vnode_np1.Add(Vnode_prev_it_nR_comp, n_vec_R);

         // Add contribution to the averaging objects
         tensor(n_vec, n_vec, dm_tmp);
         dm_tmp.Mult(Vnode_np1, v_tmp);
         Mi.Add(1., dm_tmp);
         Ri.Add(1., v_tmp);

         // In the case of boundary nodes, check if this is an interior face
         if (bdr_ind != -1 && FI.IsInterior())
         {
            // Increment averaging denom 
            num_faces_for_avg += 1;
            // We must add the ghost node contribution 
            // We will be employing reflective boundary conditions

            // Compute vector from node to face 
            Vector face_dx(dim);
            subtract(Vnode_x, face_x, face_dx);

            // Compute vector from node to adjacent node
            Vector Vadj_dx(dim);
            subtract(Vnode_x, Vadj_x, Vadj_dx);

            // Sod
            switch (bdr_ind)
            {
            case 1: // bottom
            case 3: // top
               n_vec[0] *= -1.; // preserves orientation of face
               face_dx[0] *= -1.;
               Vadj_dx[0] *= -1.;
               Vadj_n[1] *= -1.;
               Vface_n[1] *= -1.;
               break;
            
            case 2: // right
            case 4: // left
               n_vec[1] *= -1.; // preserves orientation of face
               face_dx[1] *= -1.;
               Vadj_dx[1] *= -1.;
               Vadj_n[0] *= -1.;
               Vface_n[0] *= -1.;
               break;
            
            case -1:
               // Not a boundary vertex
               continue;

            default:
               MFEM_ABORT("Incorrect bdr attribute encountered while enforcing mesh velocity BCs.\n");
               break;
            }

            // Recompute bmn for face
            double bmn = Vf * n_vec;
            bmn *= F;

            n_vec_R = n_vec;
            geom.Orthogonal(n_vec_R);
            Vnode_prev_it_nR_comp = Vnode_n * n_vec_R;

            // Get normal and rotated components of Vadj_n
            Vadj_n_comp = Vadj_n * n_vec;
            Vadj_nR_comp = Vadj_n * n_vec_R;
            
            /* Get flipped face_x and vadj_x */
            add(Vnode_x, face_dx, face_x);
            add(Vnode_x, Vadj_dx, Vadj_x);
            // and half location
            add(Vadj_x, dt/2., Vadj_n, Vadj_half);

            Vector ttvec(dim);
            subtract(Vnode_x, Vadj_x, ttvec);
            geom.Orthogonal(ttvec);
            double _val = ttvec * n_vec;
            if (_val < 1E-6)
            {
               cout << "boundary analysis:\n";
               cout << "node: " << node << endl;
               cout << "nodex: ";
               Vnode_x.Print();
               cout << "Vadj_x: ";
               Vadj_x.Print(cout);
               cout << "nvec: ";
               n_vec.Print(cout);
               assert(false);
            }

            // Compute vector constants for ghost node
            Bnode = face_x;
            Bnode *= -2.;
            Bnode.Add(-.5, Vadj_x);
            Bnode.Add(2.5, Vnode_x);
            geom.Orthogonal(Bnode);

            Badj = face_x;
            Badj *= -2.;
            Badj.Add(-.5, Vnode_x);
            Badj.Add(2.5, Vadj_x);
            geom.Orthogonal(Badj);

            double badjn = Badj * n_vec;
            double badjnr = Badj * n_vec_R;
            double bnoden = Bnode * n_vec;
            double bnodenr = Bnode * n_vec_R;

            /* Solve for Vnode_n_comp */
            numer = 3. * bmn;
            // Add in Dc1V3nper contribution (explicit)
            // V3nperp is set to 0 in option 1.
            if (mv_it_option == 2 || mv_it_option == 3)
            {
               add(Vadj_x, dt/2., Vadj_n, Vadj_half);

               subtract(Vadj_n, Vnode_n, temp_vec); // V1 - V2 = temp_vec
               subtract(Vnode_half, Vadj_half, temp_vec_2); // A2-A1
               
               geom.Orthogonal(temp_vec_2);
               D = dt * (temp_vec * n_vec_R) + 2. * (n_vec * temp_vec_2);

               // Compute c1 (A.4a)
               subtract(Vnode_n, Vadj_n, temp_vec); // only change temp_vec, since temp_vec_2 is same from D calculation (half step representation)
               
               Dc1 = dt * (temp_vec * n_vec) + 2. * (temp_vec_2 * n_vec_R); 

               /* Compute V3nperp using previous iteration */
               add(0.5, Vnode_n, 0.5, Vadj_n, temp_vec);
               V3nperp = -1. * (temp_vec * n_vec_R);
               numer += Dc1*V3nperp;
            }
            numer += (badjn - dt * Vadj_nR_comp) * Vadj_n_comp;
            numer += badjnr * Vadj_nR_comp;
            numer += ((dt/2) * Vadj_n_comp - bnodenr)*Vnode_prev_it_nR_comp;

            double denom = -dt*Vnode_prev_it_nR_comp + (dt/2)*Vadj_nR_comp + bnoden;
      
            Vnode_n_comp = numer / denom;

            // Add in contribution to predicted nodal velocity for this face
            Vnode_np1 = 0.;
            Vnode_np1.Add(Vnode_n_comp, n_vec);
            Vnode_np1.Add(Vnode_prev_it_nR_comp, n_vec_R);

            // Add contribution to the averaging objects
            tensor(n_vec, n_vec, dm_tmp);
            dm_tmp.Mult(Vnode_np1, v_tmp);
            Mi.Add(1., dm_tmp);
            Ri.Add(1., v_tmp);
         } // End ghost node
      }

      // Average node_v
      Mi.Invert();
      Mi.Mult(Ri, predicted_node_v);

      // Theta average with previous velocity
      if (do_theta_averaging)
      {
         predicted_node_v *= theta; 
         predicted_node_v.Add(1. - theta, Vnode_n);
      }

      // Put velocity in S
      // Gauss-Seidel > Jacobi
      geom.UpdateNodeVelocity(S, node, predicted_node_v);
   }
}


/****************************************************************************************************
* Function: ComputeIterationNormMC
* Parameters:
*  S        - BlockVector representing FiniteElement information
*  dt       - Current time step
*
* Purpose:
* 
*  NOTE: Interior faces. 
*  
****************************************************************************************************/
double LagrangianLOOperator::ComputeIterationNormMC(Vector &S, const double & dt)
{
   double val = 0., denom_val = 0.;
   int num_broken = 0, total_num = 0;
   
   /* Parameters needed for face velocity calculations */
   mfem::Mesh::FaceInformation FI;
   Vector Vf(dim), n_int(dim), n_vec(dim), face_velocity(dim);
   Vector Uc(dim+2), Ucp(dim+2);
   Vector cell_c_v(dim), cell_cp_v(dim), cell_center_v(dim); // For computation of bmn
   Array<int> row;
   int face_dof = 0, face_vdof2 = 0, c = 0, cp = 0;

   double V3n, V3nperp;

   // Iterate over faces
   for (int face = 0; face < num_faces; face++) // face iterator
   {  
      // cout << "face: " << face << endl;
      // Get intermediate face velocity, face information, and face normal
      face_velocity = 0.;
      Vf = 0.;

      FI = pmesh->GetFaceInformation(face);
      if (FI.IsInterior())
      {
         /* adjacent corner indices */
         H1.GetFaceDofs(face, row);
         int face_vdof1 = row[1], face_vdof2 = row[0], face_dof = row[2]; // preserve node orientation discussed in appendix A

         // retrieve old face and corner locations
         Vector face_x(dim), vdof1_x(dim), vdof2_x(dim), vdof1_v(dim), vdof2_v(dim);
         geom.GetNodePositionFromBV(S, face_dof, face_x);
         geom.GetNodePositionFromBV(S, face_vdof1, vdof1_x);
         geom.GetNodePositionFromBV(S, face_vdof2, vdof2_x);

         // retrieve corner velocities
         geom.GetNodeVelocity(S, face_vdof1, vdof1_v);
         geom.GetNodeVelocity(S, face_vdof2, vdof2_v);

         // cout << "node 1 velocity (" << face_vdof1 << "): ";
         // vdof1_v.Print(cout);
         // cout << "node 2 velocity (" << face_vdof2 << "): ";
         // vdof2_v.Print(cout);

         // if (FI.IsInterior())
         // {
         // Get adjacent cell information
         c = FI.element[0].index;
         cp = FI.element[1].index;
         GetCellStateVector(S, c, Uc);
         GetCellStateVector(S, cp, Ucp);
         pb->velocity(Uc, cell_c_v);
         pb->velocity(Ucp, cell_cp_v);

         // Calculate outer normal
         CalcOutwardNormalInt(S, c, face, n_int);
         n_vec = n_int;
         double F = n_vec.Norml2();
         n_vec /= F;

         if (1. - n_vec.Norml2() > 1e-12)
         {
            cout << "n_vec: ";
            n_vec.Print(cout);
            MFEM_ABORT("Invalid normal vector.\n");
         }

         // Calculate new corner locations and half  step locations
         Vector vdof1_x_new(dim), vdof2_x_new(dim), vdof1_x_half(dim), vdof2_x_half(dim);
         vdof1_x_new = vdof1_x;
         vdof1_x_new.Add(dt, vdof1_v);
         vdof1_x_half = vdof1_x;
         vdof1_x_half.Add(dt/2., vdof1_v);

         vdof2_x_new = vdof2_x;
         vdof2_x_new.Add(dt, vdof2_v);
         vdof2_x_half = vdof2_x;
         vdof2_x_half.Add(dt/2., vdof2_v);

         // calculate a_{12}^{n+1} (new tangent midpoint)
         Vector vdof12_x_new(dim);
         vdof12_x_new = vdof1_x_new;
         vdof12_x_new += vdof2_x_new;
         vdof12_x_new /= 2.;

         // Compute D (A.4c)
         Vector n_vec_R(dim), temp_vec(dim), temp_vec_2(dim);
         n_vec_R = n_vec;
         geom.Orthogonal(n_vec_R);
         subtract(vdof1_v, vdof2_v, temp_vec); // V1 - V2 = temp_vec
         subtract(vdof2_x_half, vdof1_x_half, temp_vec_2); // A2-A1
         geom.Orthogonal(temp_vec_2);
         double D = dt * (temp_vec * n_vec_R) + 2. * (n_vec * temp_vec_2);

         // Compute c1 (A.4a)
         subtract(vdof2_v, vdof1_v, temp_vec); // only change temp_vec, since temp_vec_2 is same from D calculation (half step representation)
         double c1 = (1./D) * (dt * (temp_vec * n_vec) + 2. * (temp_vec_2 * n_vec_R)); 

         /* Compute V3nperp using previous iteration */
         add(0.5, vdof2_v, 0.5, vdof1_v, temp_vec);
         V3nperp = -1. * (temp_vec * n_vec_R);

         // Compute c0 (A.4b)
         // Vector n_vec_half(dim);
         // subtract(vdof2_x_half, vdof1_x_half, n_vec_half);
         // geom.Orthogonal(n_vec_half);
         GetIntermediateFaceVelocity(face, Vf);
         
         double bmn = Vf * n_vec;
         bmn *= F;

         temp_vec = vdof1_x_half;
         geom.Orthogonal(temp_vec); // A1R
         temp_vec_2 = vdof2_x_half;
         geom.Orthogonal(temp_vec_2); // A2R
         double const1 = vdof1_v * temp_vec - vdof2_v * temp_vec_2; // V1*A1R - V2*A2R
         double const2 = vdof1_v * temp_vec_2 - vdof2_v * temp_vec; // V1*A2R - V2*A1R
         
         temp_vec = face_x;
         geom.Orthogonal(temp_vec);
         subtract(vdof2_v, vdof1_v, temp_vec_2);
         double const3 = temp_vec_2 * temp_vec; // (V2 - V1) * a3nR
         double c0 = (3. / D) * (bmn + const1 / 2. + const2 / 6. + 2. * const3 / 3.);
         
         // Add to val for interior faces only
         add(vdof1_v, vdof2_v, temp_vec);
         double av_nc = (temp_vec * n_vec) / 2.;
         
         double rhs = c0;
         // RHS is more than just c0 in mv_it_option 2/3
         if (mv_it_option == 2 || mv_it_option == 3)
         {
            rhs += c1 * V3nperp;
         }
         double fval = av_nc - rhs;

         if (abs(fval) > 1.e-8)
         {
            num_broken++;
         }
         // cout << "denom val: " << denom_val << endl;

         denom_val += abs(av_nc);
         val += abs(fval);
         total_num++;
      }
   }
   double perc_broken = double(num_broken)/double(total_num);
   // cout << "percentage of broken faces: " << perc_broken << endl;
   return val / denom_val;
}


/****************************************************************************************************
* Function: IterativeCornerVelocityLS
* Parameters:
*  S    - BlockVector corresponding to nth timestep that stores mesh information, 
*         mesh velocity, and state variables.
*  dt   - Current timestep size
*
* Purpose:
*  This function constitutes one iteration on the previously computed corner node 
*  velocity to reduce total face motion while still preserving mass conservation.
*  The idea is an addition to the previous iteration method, where now we solve a 
*  least squares system where the unknowns are the normal components of both nodal
*  velocities adjacent to the given face. This method was initially discussed on 
*  05/20/2024.
*
*  Note: This function relies on the mesh velocities having already been computed
*  and possibly linearized.  One can use whichever mesh velocity computation before 
*  iteration.
*
*  Note: This function will modify mv_gf in the BlockVector S
****************************************************************************************************/
void LagrangianLOOperator::IterativeCornerVelocityLS(Vector &S, const double & dt)
{
   // cout << "=====IterativeCornerVelocityLS=====\n";
   // Optional run time parameters
   bool do_theta_averaging = false;
   double theta = 0.5;
   bool _add_viscosity = true;
   double visc_weight = 1.;
   double visc_wn = 1., visc_wt = 0.;

   // Values needed during iteration
   Vector S_new = S;
   Array<int> faces_row, face_dofs_row;
   int faces_length;
   mfem::Mesh::FaceInformation FI;
   H1.ExchangeFaceNbrData();
   int Vadj_index, c, bdr_ind = 0;
   double F = 0.;
   double Vnode_prev_it_nR_comp = 0., Vadj_n_comp = 0., Vadj_nR_comp = 0.;

   Vector predicted_node_v(dim), Vf(dim), Vnode_np1(dim), Vnode_prev_it(dim);
   Vector tau_vec(dim), n_int(dim), n_vec(dim), n_vec_R(dim);
   Vector Vnode_x(dim), Vadj_x(dim), face_x(dim);
   Vector Vnode_xnp1(dim), Vadj_xnp1(dim), face_xnp1(dim), a12_xnp1(dim);
   Vector Badj(dim), Bnode(dim);

   Vector Vnode_n(dim), Vadj_n(dim), Vface_n(dim);
   Vector Vnode_half(dim), Vadj_half(dim);
   Vector temp_vec(dim), temp_vec_2(dim);
   double D = 0., c1 = 0., c0 = 0., V3nperp = 0., Dc1 = 0.;

   // Averaging
   DenseMatrix Mi(dim), dm_tmp(dim), dm_tmp2(dim);
   Vector Ri(dim);

   /* Iterate over corner nodes */
   for (int node = 0; node < NVDofs_H1; node++) // TODO: Is NDofs_H1L == NVDofs_H1?
   {
      // Reset new nodal velocity, and averaging objects
      Mi = 0., Ri = 0.;
      predicted_node_v = 0.;

      // Set corresponding boundary indicator
      // If node is boundary node, will have to add corresponding correction
      // node is interior node if bdr_ind == -1
      bdr_ind = BdrVertexIndexingArray[node];

      // Get current nodal velocity and position 
      geom.GetNodeVelocity(S, node, Vnode_n);
      geom.GetNodePositionFromBV(S, node, Vnode_x);

      // Compute An and an+1 for node
      add(Vnode_x, dt/2., Vnode_n, Vnode_half);
      add(Vnode_x, dt, Vnode_n, Vnode_xnp1);

      // Get cell faces
      vertex_edge.GetRow(node, faces_row);
      faces_length = faces_row.Size();

      /* Iterate over cell faces */
      for (int face_it = 0; face_it < faces_length; face_it++) // Adjacent face iterator
      {
         // Get face information
         int face = faces_row[face_it];
         GetIntermediateFaceVelocity(face, Vf);
         FI = pmesh->GetFaceInformation(face);

         /* adjacent corner indices */
         // preserve node orientation where cell to right 
         // of face is with lower cell index
         H1.GetFaceDofs(face, face_dofs_row);
         int face_vdof1 = face_dofs_row[1], 
             face_vdof2 = face_dofs_row[0], 
             face_dof = face_dofs_row[2]; 

         // Grab corresponding vertex velocity from S
         // We are always solving for V2. 
         // If the index for Vnode does not match that
         // for V2, we must flip the normal.
         if (node == face_vdof1) {
            // Get adj index
            Vadj_index = face_vdof2;
         } else {
            // Get adj index
            Vadj_index = face_vdof1;
            // Node matches v2
         }

         /* Get face and adj node velocity and position */
         geom.GetNodeVelocity(S, Vadj_index, Vadj_n);
         geom.GetNodeVelocity(S, face_dof, Vface_n);
         geom.GetNodePositionFromBV(S, Vadj_index, Vadj_x);
         geom.GetNodePositionFromBV(S, face_dof, face_x);

         // Get n, tau, and F
         subtract(Vnode_x, Vadj_x, tau_vec);
         F = tau_vec.Norml2();
         tau_vec /= F;
         n_vec = tau_vec;
         geom.Orthogonal(n_vec);

         // Compute future face and adjacent node locations
         add(Vadj_x, dt/2., Vadj_n, Vadj_half);
         add(Vadj_x, dt, Vadj_n, Vadj_xnp1);
         // add(0.5, Vadj_xnp1, 0.5, Vnode_xnp1, a12_xnp1);

         // Perpendicular vector
         n_vec_R = n_vec;
         geom.Orthogonal(n_vec_R);

         // Calculate bmn
         double bmn = Vf * n_vec;
         bmn *= F;

         // Get Vnode_prev_it component in tangent direction from previous iteration
         Vnode_prev_it_nR_comp = Vnode_n * n_vec_R;

         // Get normal and rotated components of Vadj_n
         Vadj_n_comp = Vadj_n * n_vec;
         Vadj_nR_comp = Vadj_n * n_vec_R;
         
         // Compute geometrical vectors
         Bnode = face_x;
         Bnode *= -2.;
         Bnode.Add(-.5, Vadj_x);
         Bnode.Add(2.5, Vnode_x);
         geom.Orthogonal(Bnode);

         Badj = face_x;
         Badj *= -2.;
         Badj.Add(-.5, Vnode_x);
         Badj.Add(2.5, Vadj_x);
         geom.Orthogonal(Badj);

         // Compute dot products
         double badjn = Badj * n_vec;
         double badjnr = Badj * n_vec_R;
         double bnoden = Bnode * n_vec;
         double bnodenr = Bnode * n_vec_R;

         // Compute aFi, aFj, bF
         double aFi = 0.5 * dt * Vadj_nR_comp;
         aFi -= dt*Vnode_prev_it_nR_comp;
         aFi += bnoden;
         double aFj = dt*Vadj_nR_comp;
         aFj -= 0.5*dt*Vnode_prev_it_nR_comp;
         aFj -= badjn;
         double bF = 3 * bmn;
         bF += badjnr * Vadj_nR_comp;
         bF -= bnodenr * Vnode_prev_it_nR_comp;

         // Add in Dc1V3nper contribution (explicit)
         // V3nperp is set to 0 in option 1.
         if (mv_it_option == 2 || mv_it_option == 3)
         {
            // Calculate D and c1 (perturbations only need to be handled explicitly)
            subtract(Vadj_n, Vnode_n, temp_vec); // V1 - V2 = temp_vec
            subtract(Vnode_half, Vadj_half, temp_vec_2); // A2-A1

            geom.Orthogonal(temp_vec_2);
            D = dt * (temp_vec * n_vec_R) + 2. * (n_vec * temp_vec_2);

            // Compute c1 (A.4a)
            subtract(Vnode_n, Vadj_n, temp_vec); // only change temp_vec, since temp_vec_2 is same from D calculation (half step representation)
            Dc1 = dt * (temp_vec * n_vec) + 2. * (temp_vec_2 * n_vec_R); 

            /* Compute V3nperp using previous iteration */
            add(0.5, Vnode_n, 0.5, Vadj_n, temp_vec);
            V3nperp = -1. * (temp_vec * n_vec_R);

            /* Add in Dc1 explicitly at k to bF */
            bF += Dc1*V3nperp;
         }

         // Add contribution to the averaging objects
         double Mcoeff = 2. * pow(aFi, 2);
         double Rcoeff = aFj*Vadj_n_comp - bF;
         Rcoeff *= -2. * aFi;
         tensor(n_vec, n_vec, dm_tmp);

         Mi.Add(Mcoeff, dm_tmp);
         Ri.Add(Rcoeff, n_vec);

         /* Optionally, add viscosity to the system */
         if (_add_viscosity)
         {
            // Add to Mi
            double _w = 2 * visc_weight * pow(F, 2);
            Mi.Add(_w * visc_wn, dm_tmp);
            tensor(tau_vec, tau_vec, dm_tmp2);
            Mi.Add(_w * visc_wt, dm_tmp2);

            // Add to Ri
            dm_tmp.Mult(Vadj_n, temp_vec);
            dm_tmp2.Mult(Vadj_n, temp_vec_2);
            Ri.Add(_w * visc_wn, temp_vec);
            Ri.Add(_w * visc_wt, temp_vec_2);
         }

         // In the case of boundary nodes, check if this is an interior face
         if (bdr_ind != -1 && FI.IsInterior())
         {
            // We must add the ghost node contribution 
            // We will be employing reflective boundary conditions

            // Compute vector from node to face 
            Vector face_dx(dim);
            subtract(Vnode_x, face_x, face_dx);

            // Compute vector from node to adjacent node
            Vector Vadj_dx(dim);
            subtract(Vnode_x, Vadj_x, Vadj_dx);

            // Sod
            switch (bdr_ind)
            {
            case 1: // bottom
            case 3: // top
               n_vec[0] *= -1.; // preserves orientation of face
               tau_vec[1] *= -1.;
               face_dx[0] *= -1.;
               Vadj_dx[0] *= -1.;
               Vadj_n[1] *= -1.;
               Vface_n[1] *= -1.;
               break;
            
            case 2: // right
            case 4: // left
               n_vec[1] *= -1.; // preserves orientation of face
               tau_vec[0] *= -1.;
               face_dx[1] *= -1.;
               Vadj_dx[1] *= -1.;
               Vadj_n[0] *= -1.;
               Vface_n[0] *= -1.;
               break;
            
            case -1:
               // Not a boundary vertex
               continue;

            default:
               MFEM_ABORT("Incorrect bdr attribute encountered while enforcing mesh velocity BCs.\n");
               break;
            }

            // Recompute bmn for face
            double bmn = Vf * n_vec;
            bmn *= F;

            n_vec_R = n_vec;
            geom.Orthogonal(n_vec_R);
            Vnode_prev_it_nR_comp = Vnode_n * n_vec_R;

            // Get normal and rotated components of Vadj_n
            Vadj_n_comp = Vadj_n * n_vec;
            Vadj_nR_comp = Vadj_n * n_vec_R;
            
            /* Get flipped face_x and vadj_x */
            add(Vnode_x, face_dx, face_x);
            add(Vnode_x, Vadj_dx, Vadj_x);
            // and half location
            add(Vadj_x, dt/2., Vadj_n, Vadj_half);

            // Compute vector constants for ghost node
            Bnode = face_x;
            Bnode *= -2.;
            Bnode.Add(-.5, Vadj_x);
            Bnode.Add(2.5, Vnode_x);
            geom.Orthogonal(Bnode);

            Badj = face_x;
            Badj *= -2.;
            Badj.Add(-.5, Vnode_x);
            Badj.Add(2.5, Vadj_x);
            geom.Orthogonal(Badj);

            // Compute dot products
            double badjn = Badj * n_vec;
            double badjnr = Badj * n_vec_R;
            double bnoden = Bnode * n_vec;
            double bnodenr = Bnode * n_vec_R;

            // Compute aFi, aFj, bF
            double aFi = 0.5 * dt * Vadj_nR_comp;
            aFi -= dt*Vnode_prev_it_nR_comp;
            aFi += bnoden;
            double aFj = dt*Vadj_nR_comp;
            aFj -= 0.5*dt*Vnode_prev_it_nR_comp;
            aFj -= badjn;
            double bF = 3 * bmn;
            bF += badjnr * Vadj_nR_comp;
            bF -= bnodenr * Vnode_prev_it_nR_comp;

            // Add in Dc1V3nper contribution (explicit)
            // V3nperp is set to 0 in option 1.
            if (mv_it_option == 2 || mv_it_option == 3)
            {
               add(Vadj_x, dt/2., Vadj_n, Vadj_half);

               subtract(Vadj_n, Vnode_n, temp_vec); // V1 - V2 = temp_vec
               subtract(Vnode_half, Vadj_half, temp_vec_2); // A2-A1
               
               geom.Orthogonal(temp_vec_2);
               D = dt * (temp_vec * n_vec_R) + 2. * (n_vec * temp_vec_2);

               // Compute c1 (A.4a)
               subtract(Vnode_n, Vadj_n, temp_vec); // only change temp_vec, since temp_vec_2 is same from D calculation (half step representation)
               
               Dc1 = dt * (temp_vec * n_vec) + 2. * (temp_vec_2 * n_vec_R); 

               /* Compute V3nperp using previous iteration */
               add(0.5, Vnode_n, 0.5, Vadj_n, temp_vec);
               V3nperp = -1. * (temp_vec * n_vec_R);

               /* Add in Dc1 explicitly at k to bF */
               bF += Dc1*V3nperp;
            }

            // Add contribution to the averaging objects
            double Mcoeff = 2. * pow(aFi, 2);
            double Rcoeff = aFj*Vadj_n_comp - bF;
            Rcoeff *= -2. * aFi;
            tensor(n_vec, n_vec, dm_tmp);

            Mi.Add(Mcoeff, dm_tmp);
            Ri.Add(Rcoeff, n_vec);

            /* Optionally, add viscosity to the system */
            if (_add_viscosity)
            {
               // Add to Mi
               double _w = 2 * visc_weight * pow(F, 2);
               Mi.Add(_w * visc_wn, dm_tmp);
               tensor(tau_vec, tau_vec, dm_tmp2);
               Mi.Add(_w * visc_wt, dm_tmp2);

               // Add to Ri
               dm_tmp.Mult(Vadj_n, temp_vec);
               dm_tmp2.Mult(Vadj_n, temp_vec_2);
               Ri.Add(_w * visc_wn, temp_vec);
               Ri.Add(_w * visc_wt, temp_vec_2);
            }
         } // End ghost node
      }

      // Average node_v
      Mi.Invert();
      Mi.Mult(Ri, predicted_node_v);

      // Theta average with previous velocity
      if (do_theta_averaging)
      {
         predicted_node_v *= theta; 
         predicted_node_v.Add(1. - theta, Vnode_n);
      }

      // Put velocity in S
      // Gauss-Seidel > Jacobi
      // geom.UpdateNodeVelocity(S, node, predicted_node_v);
      // Try Jacobi
      geom.UpdateNodeVelocity(S_new, node, predicted_node_v);
   } // End node iterator
   S = S_new;
}


/****************************************************************************************************
* Function: compute_alpha
* Parameters:
*
* Purpose:
****************************************************************************************************/
double LagrangianLOOperator::compute_alpha(const Vector &Vadj, const Vector &AadjR, 
                                                const Vector &Vnode, const Vector &AnodeR, 
                                                const Vector &V3n, const Vector &a3nR, 
                                                const Vector &n_vec, const Vector &tau_vec,
                                                const double &dt, const double &F)
{
   double ret_val = 0.;

   Vector tmp_vec(dim);
   double tmp1, tmp2, tmp3;
   tmp1 = Vnode * AnodeR - Vadj * AadjR;
   tmp2 = Vadj * AnodeR - Vnode * AadjR; 
   subtract(Vadj, Vnode, tmp_vec); // Va - Vn
   tmp3 = tmp_vec * a3nR;
   ret_val = 0.5 * tmp1 - (1./6.) * tmp2 + (2./3.) * tmp3;

   double tau_comp = tmp_vec * n_vec;
   tau_comp *= 2. * dt / 3.;
   double V3n_tau = V3n * tau_vec;
   ret_val += tau_comp * V3n_tau;

   tmp_vec *= -1; // Vn - Va
   double n_comp = tmp_vec * tau_vec;
   n_comp *= dt;
   n_comp += F;
   n_comp *= 2./3.;
   double V3n_norm = V3n * n_vec;
   ret_val += n_comp * V3n_norm;

   return ret_val;
}


/****************************************************************************************************
* Function: ComputeIterativeLSGamma
* Parameters:
*  S        - BlockVector representing FiniteElement information
*  dt       - Current time step
*
* Purpose:
*  Compute object we are trying to minimize, to include viscosity.
*  NOTE: Interior vertices. 
****************************************************************************************************/
double LagrangianLOOperator::ComputeIterativeLSGamma(Vector &S, const double & dt)
{
   double ret_val = 0.;
   double vw = 1., wn = 1., wt = 0.;
   Vector Vnode(dim), Vadj(dim), Vface(dim);
   Vector anode(dim), aadj(dim), a3n(dim);
   Vector AadjR(dim), AnodeR(dim), a3nR(dim);
   Vector n_vec(dim), tau_vec(dim);
   mfem::Mesh::FaceInformation FI;
   Vector Vf(dim);
   double bmn = 0., F = 0.;
   Array<int> faces_row, face_dofs_row;

   int Vadj_index;
   /* Iterate over corner nodes */
   for (int node = 0; node < NDofs_H1L; node++)
   { 
      // Check if node is on boundary
      int bdr_ind = BdrVertexIndexingArray[node];

      // Get current nodal velocity and position 
      geom.GetNodeVelocity(S, node, Vnode);
      geom.GetNodePositionFromBV(S, node, anode);

      // Get cell faces
      vertex_edge.GetRow(node, faces_row);
      int faces_length = faces_row.Size();

      if (bdr_ind == -1)
      {
         /* Iterate over cell faces */
         for (int face_it = 0; face_it < faces_length; face_it++) // Adjacent face iterator
         {
            // Get face information
            int face = faces_row[face_it];
            GetIntermediateFaceVelocity(face, Vf);
            FI = pmesh->GetFaceInformation(face);

            /* adjacent corner indices */
            // preserve node orientation where cell to right 
            // of face is with lower cell index
            H1.GetFaceDofs(face, face_dofs_row);
            int face_vdof1 = face_dofs_row[1], 
               face_vdof2 = face_dofs_row[0], 
               face_dof = face_dofs_row[2]; 

            // Grab corresponding vertex velocity from S
            // We are always solving for V2. 
            // If the index for Vnode does not match that
            // for V2, we must flip the normal.
            if (node == face_vdof1) {
               // Get adj index
               Vadj_index = face_vdof2;
            } else {
               // Get adj index
               Vadj_index = face_vdof1;
               // Node matches v2
            }

            /* Get face and adj node velocity and position */
            geom.GetNodeVelocity(S, Vadj_index, Vadj);
            // geom.GetNodeVelocity(S, face_dof, Vface);
            add(0.5, Vnode, 0.5, Vadj, Vface);
            geom.GetNodePositionFromBV(S, Vadj_index, aadj);
            geom.GetNodePositionFromBV(S, face_dof, a3n);

            // Get n, tau, and F
            subtract(anode, aadj, tau_vec);
            F = tau_vec.Norml2();
            tau_vec /= F;
            n_vec = tau_vec;
            geom.Orthogonal(n_vec);

            // Compute AadjR, AnodeR, a3nR
            add(1., aadj, dt / 2., Vadj, AadjR);
            geom.Orthogonal(AadjR);
            add(1., anode, dt/2., Vnode, AnodeR);
            geom.Orthogonal(AnodeR);
            a3nR = a3n;
            geom.Orthogonal(a3nR);

            double alpha = compute_alpha(Vadj, AadjR, Vnode, AnodeR, Vface, a3nR, n_vec, tau_vec, dt, F);

            double Vnode_n = Vnode * n_vec;
            double Vnode_t = Vnode * tau_vec;
            double Vadj_n = Vadj * n_vec;
            double Vadj_t = Vadj * tau_vec;
            ret_val += pow(alpha - bmn, 2) + vw * pow(F,2) * (wn * pow(Vnode_n - Vadj_n, 2) + wt * pow(Vnode_t - Vadj_t, 2));
         }
      }
   }

   return ret_val;
}


/****************************************************************************************************
* Function: ComputeIterationNorm
* Parameters:
*  S             - BlockVector representing FiniteElement information
*  mv_gf_prev_it - mesh velocity grid function from the previous iteration 
*  dt            - Current time step
*
* Purpose:
*  Checks how well the iterative corner velocity methods are as converging
*  the corner velocities.  Computes
*                        | V_i^(k+1) - V_i^(k) |
*                 SUM_i  -----------------------
*                              | V_i^(k) |
*  NOTE: Interior vertices. 
****************************************************************************************************/
double LagrangianLOOperator::ComputeIterationNorm(const Vector &S, const ParGridFunction &mv_gf_prev_it, const double & dt)
{
   Vector Vi_prev(dim), Vi_next(dim), temp_vec(dim);
   double numer = 0., denom = 0.;
   /* Iterate over corner nodes */
   for (int node = 0; node < NDofs_H1L; node++)
   { 
      // Check if node is on boundary
      int bdr_ind = BdrVertexIndexingArray[node];

      if (bdr_ind == -1)
      {
         // Get updated iterated velocity and old velocity
         geom.GetNodeVelocity(S, node, Vi_next);
         geom.GetNodeVelocity(mv_gf_prev_it, node, Vi_prev);
         subtract(Vi_next, Vi_prev, temp_vec);
         numer += temp_vec.Norml2();
         denom += Vi_prev.Norml2();
      }
   }
   return numer / denom;
}


/****************************************************************************************************
* Function: ComputeFaceSecantNorm
* Parameters:
*  S        - BlockVector representing FiniteElement information
*  dt       - Current time step
*
* Purpose:
*  Checks how close the corrective face velocity is to the average of the adjacent 
*  corner node velocities.

*  NOTE: Interior vertices. 
****************************************************************************************************/
double LagrangianLOOperator::ComputeFaceSecantNorm(Vector &S, const double & D)
{
   double num = 0., denom = 0.;
   int num_broken = 0, total_num = 0;
   
   /* Parameters needed for face velocity calculations */
   mfem::Mesh::FaceInformation FI;
   Vector Vdof1_n(dim), Vdof2_n(dim), Vface_n(dim), temp_vec(dim);
   Array<int> row;

   // Iterate over faces
   for (int face = 0; face < num_faces; face++) // face iterator
   {
      FI = pmesh->GetFaceInformation(face);
      if (FI.IsInterior())
      {
         /* adjacent corner indices */
         H1.GetFaceDofs(face, row);
         int face_vdof1 = row[1], face_vdof2 = row[0], face_dof = row[2]; // preserve node orientation discussed in appendix A

         // retrieve corner velocities
         geom.GetNodeVelocity(S, face_vdof1, Vdof1_n);
         geom.GetNodeVelocity(S, face_vdof2, Vdof2_n);
         geom.GetNodeVelocity(S, face_dof, Vface_n);

         add(0.5, Vdof1_n, 0.5, Vdof2_n, temp_vec);
         denom += temp_vec.Norml2();
         temp_vec -= Vface_n;
         double val = temp_vec.Norml2();
         num += val;
         if (val > 1e-12)
         {
            num_broken++;
         }

         total_num++;
      }
   }
   double perc_broken = double(num_broken)/double(total_num);
   // cout << "percentage of broken faces: " << perc_broken << endl;
   return num / denom;
}


/****************************************************************************************************
* Function: IterativeCornerVelocityTNLSnoncart
* Parameters:
*  S    - BlockVector corresponding to nth timestep that stores mesh information, 
*         mesh velocity, and state variables.
*  dt   - Current timestep size
*
* Purpose:
*  This function constitutes one iteration on the previously computed corner node 
*  velocity to reduce total face motion while still preserving mass conservation.
*  The idea is an addition to the LS iteration method, where now we solve a 
*  least squares system composed of both the normal and tangential components.
*  This method was initially discussed on 06/05/2024.
*
*  Note: This function relies on the mesh velocities having already been computed
*  and possibly linearized.  One can use whichever mesh velocity computation before 
*  iteration.
*
*  Note: This function will modify mv_gf in the BlockVector S
*  Note: The basis taken for the l2 norm is that formed by the normal and tangent
*  vectors at the half step.
****************************************************************************************************/
void LagrangianLOOperator::IterativeCornerVelocityTNLSnoncart(Vector &S, const double & dt)
{
   // cout << "=====IterativeCornerVelocityTNLSnoncart=====\n";
   // Optional run time parameters
   bool do_theta_averaging = true;
   double theta = 0.5; // Works best for distorted sod.
   Vector S_new = S;

   // Optionally weight the normal and the tangent contributions
   double wn = 1.;
   double wt = 1.;

   // Values needed during iteration
   Array<int> faces_row, face_dofs_row;
   int faces_length;
   mfem::Mesh::FaceInformation FI;
   H1.ExchangeFaceNbrData();
   int Vadj_index, c, bdr_ind = 0;
   double F = 0.;
   double Vadj_n_comp = 0., Vadj_tau_comp = 0., vnode_tau_comp_prev = 0.;

   Vector predicted_node_v(dim), Vf(dim);
   Vector n_int(dim), n_vec(dim), tau_vec(dim);
   Vector anode_n(dim), aadj_n(dim), a3n(dim);
   Vector anode_np1(dim), aadj_np1(dim);
   Vector Badj(dim), Bnode(dim);

   Vector Vnode_n(dim), Vadj_n(dim);
   Vector Anode(dim), Aadj(dim);
   Vector n_vec_half(dim), tau_vec_half(dim);
   Vector temp_vec(dim);

   // Averaging
   DenseMatrix Mi(dim), dm_tmp(dim), dm_tmp2(dim);
   Vector Ri(dim);

   /* Iterate over corner nodes */
   for (int node = 0; node < NDofs_H1L; node++)
   {
      // Reset new nodal velocity, and averaging objects
      Mi = 0., Ri = 0.;
      predicted_node_v = 0.;

      // Set corresponding boundary indicator
      // If node is boundary node, will have to add corresponding correction
      // node is interior node if bdr_ind == -1
      bdr_ind = BdrVertexIndexingArray[node];

      // Get current nodal velocity and position 
      geom.GetNodeVelocity(S, node, Vnode_n);
      geom.GetNodePositionFromBV(S, node, anode_n);

      // Compute An for node
      add(anode_n, dt/2., Vnode_n, Anode);

      // Get cell faces
      vertex_edge.GetRow(node, faces_row);
      faces_length = faces_row.Size();

      /* Iterate over cell faces */
      for (int face_it = 0; face_it < faces_length; face_it++) // Adjacent face iterator
      {
         // Get face information
         int face = faces_row[face_it];
         GetIntermediateFaceVelocity(face, Vf);
         FI = pmesh->GetFaceInformation(face);

         /* adjacent corner indices */
         // preserve node orientation where cell to right 
         // of face is with lower cell index
         H1.GetFaceDofs(face, face_dofs_row);
         int face_vdof1 = face_dofs_row[1], 
             face_vdof2 = face_dofs_row[0], 
             face_dof = face_dofs_row[2]; 

         // Grab corresponding vertex velocity from S
         // We are always solving for V2. 
         if (node == face_vdof1) {
            // Get adj index
            Vadj_index = face_vdof2;
         } else {
            // Get adj index
            Vadj_index = face_vdof1;
            // Node matches v2
         }

         /* Get face and adj node velocity and position */
         geom.GetNodeVelocity(S, Vadj_index, Vadj_n);
         geom.GetNodePositionFromBV(S, Vadj_index, aadj_n);
         geom.GetNodePositionFromBV(S, face_dof, a3n);

         // Get n, tau, and F
         subtract(anode_n, aadj_n, tau_vec);
         F = tau_vec.Norml2();
         tau_vec /= F;
         n_vec = tau_vec;
         geom.Orthogonal(n_vec);

         // Compute future adjacent node locations
         add(aadj_n, dt/2., Vadj_n, Aadj);

         // Compute half
         subtract(Anode, Aadj, tau_vec_half);
         double Fhalf = tau_vec_half.Norml2();
         tau_vec_half /= Fhalf;
         n_vec_half = tau_vec_half;
         geom.Orthogonal(n_vec_half);

         // Get normal and tangential components of Vadj_n
         Vadj_n_comp = Vadj_n * n_vec;
         Vadj_tau_comp = Vadj_n * tau_vec;
         
         // Compute geometrical vectors
         Bnode = 0.;
         Bnode.Add(1./6., aadj_n);
         Bnode.Add(0.5, anode_n);
         Bnode.Add(-2./3., a3n);
         geom.Orthogonal(Bnode);

         Badj = 0.;
         Badj.Add(0.5, aadj_n);
         Badj.Add(1./6., anode_n);
         Badj.Add(-2./3., a3n);
         geom.Orthogonal(Badj);

         // Compute dot products
         double badjn = Badj * n_vec;
         double badjtau = Badj * tau_vec;
         double bnoden = Bnode * n_vec;
         double bnodetau = Bnode * tau_vec;
         double vadjbadj = Vadj_n * Badj;

         // Compute alphas and betas
         double alpha1 = bnoden - (dt/2.) * Vadj_tau_comp + (1./3.) * F;
         alpha1 /= Fhalf;
         double alpha2 = bnodetau + (dt/2.) * Vadj_n_comp;
         alpha2 /= Fhalf;
         double alpha3 = -1. * vadjbadj + (1./3.) * F * Vadj_n_comp;
         alpha3 /= Fhalf;

         temp_vec = 0.;
         temp_vec.Add(0.5, Anode);
         temp_vec.Add(1./6., aadj_n);
         temp_vec.Add(-2./3., a3n);
         double beta1 = n_vec * temp_vec;
         beta1 /= Fhalf;
         double beta2 = tau_vec * temp_vec + (1./3.) * F;
         beta2 /= Fhalf;
         temp_vec = 0.;
         temp_vec.Add(2./3., a3n);
         temp_vec.Add(F/3., tau_vec);
         temp_vec.Add(-1./2., Aadj);
         temp_vec.Add(-1./6., anode_n);
         double beta3 = Vadj_n * temp_vec;
         beta3 /= Fhalf;

         double cn = alpha3 - Vf * n_vec_half;
         double ctau = beta3 - Vf * tau_vec_half;

         // Add contribution to the averaging objects
         tensor(n_vec, n_vec, dm_tmp);
         Mi.Add(wn*pow(alpha1,2) + wt*pow(beta1,2), dm_tmp);

         tensor(tau_vec, tau_vec, dm_tmp);
         Mi.Add(wn*pow(alpha2,2) + wt*pow(beta2,2), dm_tmp);

         tensor(n_vec, tau_vec, dm_tmp);
         tensor(tau_vec, n_vec, dm_tmp2);
         dm_tmp.Add(1., dm_tmp2);
         Mi.Add(wn*alpha1*alpha2 + wt*beta1*beta2, dm_tmp);

         Ri.Add(wn*alpha1*cn + wt*beta1*ctau, n_vec);
         Ri.Add(wn*alpha2*cn + wt*beta2*ctau, tau_vec);

         // In the case of boundary nodes, check if this is an interior face
         if (bdr_ind != -1 && FI.IsInterior())
         {
            // We must add the ghost node contribution 
            // We will be employing reflective boundary conditions

            // Compute vector from face to node 
            Vector face_dx(dim);
            subtract(anode_n, a3n, face_dx);

            // Compute vector from adjacent node to node
            Vector Vadj_dx(dim);
            subtract(anode_n, aadj_n, Vadj_dx);

            // Sod
            switch (bdr_ind)
            {
            case 1: // bottom
            case 3: // top
               n_vec[0] *= -1.; // preserves orientation of face
               tau_vec[1] *= -1.;
               face_dx[0] *= -1.;
               Vadj_dx[0] *= -1.;
               Vadj_n[1] *= -1.;
               break;
            
            case 2: // right
            case 4: // left
               n_vec[1] *= -1.; // preserves orientation of face
               tau_vec[0] *= -1.;
               face_dx[1] *= -1.;
               Vadj_dx[1] *= -1.;
               Vadj_n[0] *= -1.;
               break;
            
            case -1:
               // Not a boundary vertex
               continue;

            default:
               MFEM_ABORT("Incorrect bdr attribute encountered while enforcing mesh velocity BCs.\n");
               break;
            }

            // Get normal and rotated components of Vadj_n
            Vadj_n_comp = Vadj_n * n_vec;
            Vadj_tau_comp = Vadj_n * tau_vec;
            
            /* Get flipped a3n and vadj_x */
            add(anode_n, face_dx, a3n);
            add(anode_n, Vadj_dx, aadj_n);
            // and half location
            add(aadj_n, dt/2., Vadj_n, Aadj);

            /* Compute n, tau, F at half */
            subtract(Anode, Aadj, tau_vec_half);
            double Fhalf = tau_vec_half.Norml2();
            tau_vec_half /= Fhalf;
            n_vec_half = tau_vec_half;
            geom.Orthogonal(n_vec_half);

            // Compute geometrical vectors
            Bnode = 0.;
            Bnode.Add(1./6., aadj_n);
            Bnode.Add(0.5, anode_n);
            Bnode.Add(-2./3., a3n);
            geom.Orthogonal(Bnode);

            Badj = 0.;
            Badj.Add(0.5, aadj_n);
            Badj.Add(1./6., anode_n);
            Badj.Add(-2./3., a3n);
            geom.Orthogonal(Badj);

            // Compute dot products
            double badjn = Badj * n_vec;
            double badjtau = Badj * tau_vec;
            double bnoden = Bnode * n_vec;
            double bnodetau = Bnode * tau_vec;
            double vadjbadj = Vadj_n * Badj;

            // Compute alphas and betas
            double alpha1 = bnoden - (dt/2.) * Vadj_tau_comp + (1./3.) * F;
            alpha1 /= Fhalf;
            double alpha2 = bnodetau + (dt/2.) * Vadj_n_comp;
            alpha2 /= Fhalf;
            double alpha3 = -1. * vadjbadj + (1./3.) * F * Vadj_n_comp;
            alpha3 /= Fhalf;

            temp_vec = 0.;
            temp_vec.Add(0.5, Anode);
            temp_vec.Add(1./6., aadj_n);
            temp_vec.Add(-2./3., a3n);
            double beta1 = n_vec * temp_vec;
            beta1 /= Fhalf;
            double beta2 = tau_vec * temp_vec + (1./3.) * F;
            beta2 /= Fhalf;
            temp_vec = 0.;
            temp_vec.Add(2./3., a3n);
            temp_vec.Add(F/3., tau_vec);
            temp_vec.Add(-1./2., Aadj);
            temp_vec.Add(-1./6., anode_n);
            double beta3 = Vadj_n * temp_vec;
            beta3 /= Fhalf;

            double cn = alpha3 - Vf * n_vec_half;
            double ctau = beta3 - Vf * tau_vec_half;

            // Add contribution to the averaging objects
            tensor(n_vec, n_vec, dm_tmp);
            Mi.Add(wn*pow(alpha1,2) + wt*pow(beta1,2), dm_tmp);

            tensor(tau_vec, tau_vec, dm_tmp);
            Mi.Add(wn*pow(alpha2,2) + wt*pow(beta2,2), dm_tmp);

            tensor(n_vec, tau_vec, dm_tmp);
            tensor(tau_vec, n_vec, dm_tmp2);
            dm_tmp.Add(1., dm_tmp2);
            Mi.Add(wn*alpha1*alpha2 + wt*beta1*beta2, dm_tmp);

            Ri.Add(wn*alpha1*cn + wt*beta1*ctau, n_vec);
            Ri.Add(wn*alpha2*cn + wt*beta2*ctau, tau_vec);
         } // End ghost node
      }

      // Average node_v
      Mi.Invert();
      Ri *= -1.;
      Mi.Mult(Ri, predicted_node_v);

      // Theta average with previous velocity
      if (do_theta_averaging)
      {
         predicted_node_v *= theta; 
         predicted_node_v.Add(1. - theta, Vnode_n);
      }

      // cout << "Predicted node velocity for node " << node << ": ";
      // predicted_node_v.Print(cout);
      // Put velocity in S
      // Gauss-Seidel > Jacobi
      // geom.UpdateNodeVelocity(S, node, predicted_node_v);
      // Try Jacobi
      geom.UpdateNodeVelocity(S_new, node, predicted_node_v);
   } // End node iterator
   S = S_new;
}


/****************************************************************************************************
* Function: IterativeCornerVelocityLSCellVolumeFaceVisc
* Parameters:
*  S    - BlockVector corresponding to nth timestep that stores mesh information, 
*         mesh velocity, and state variables.
*  dt   - Current timestep size
*
* Purpose:
*  This function constitutes one iteration on the previously computed corner node 
*  velocity to reduce total face motion while still preserving mass conservation.
*  The idea is to minimize the local change in mass for each adjacent cell to a 
*  given vertex, with some optional viscosity defined on all the adjacent faces.
*
*  This method was initially discussed on 06/12/2024.
*
*  Note: This function relies on the mesh velocities having already been computed
*  and possibly linearized.  One can use whichever mesh velocity computation before 
*  iteration.
*  Note: This function will modify mv_gf in the BlockVector S
****************************************************************************************************/
void LagrangianLOOperator::IterativeCornerVelocityLSCellVolumeFaceVisc(Vector &S, const Vector &S_old, const double &dt)
{
   // cout << "=====IterativeCornerVelocityLSCellVolumeFaceVisc=====\n";
   /* Optional run time parameters */
   bool do_theta_averaging = true;
   double theta = 0.5;

   /* Objects needed for function */
   mfem::Mesh::FaceInformation FI;
   Vector* sptr_old = const_cast<Vector*>(&S_old);
   ParGridFunction sv_gf, sv_old_gf;
   sv_old_gf.MakeRef(&L2, *sptr_old, block_offsets[1]);
   sv_gf.MakeRef(&L2, S, block_offsets[1]);

   Vector S_new = S;
   Vector predicted_node_v(dim);
   Vector anode_n(dim), al1_n(dim), al2_n(dim), al3_n(dim);
   Vector Vnode_n(dim), Vl1(dim), Vl2(dim), Vl3(dim);
   Vector al1_np1(dim), al2_np1(dim), al3_np1(dim);
   Vector temp_vec(dim), temp_vec2(dim), bj_vec(dim);

   double Bj, cj;

   Array<int> faces_row, face_dofs_row, cells_row, verts;
   int faces_length, cells_length;
   int bdr_ind;
   int Vadj_index;
   Vector aadj_n(dim), Vadj_n(dim);

   DenseMatrix Mi(dim), dm_temp(dim);
   Vector Ri(dim), visci_vec(dim);
   double visc_contr = 0.;

   /* Iterate over interior nodes */
   for (int node = 0; node < NDofs_H1L; node++)
   {
      /* Reset values for new node */
      Mi = 0., Ri = 0., predicted_node_v = 0.;

      /* bdr_ind = -1 for interior nodes */
      bdr_ind = BdrVertexIndexingArray[node];
      // cout << "dt: " << dt << ", node: " << node << ", bdr_ind: " << bdr_ind << endl;

      if (bdr_ind == -1)
      {
         /* Get adjacent cells and faces */
         vertex_element->GetRow(node, cells_row);
         vertex_edge.GetRow(node, faces_row);
         cells_length = cells_row.Size();
         faces_length = faces_row.Size();

         /* Get nodal position */
         geom.GetNodePositionFromBV(S, node, anode_n);
         geom.GetNodeVelocity(S, node, Vnode_n);

         /* Iterate over adjacent cells */ 
         for (int cell_it = 0; cell_it < cells_length; cell_it++)
         {
            int el_index = cells_row[cell_it];
            // cout << "cell: " << el_index << endl;
            // cout << "cell: " << el_index << endl;
            /* Get remaining node locations and velocities */
            pmesh->GetElementVertices(el_index, verts);
            int verts_length = verts.Size();
            int node_index = verts.Find(node);
            int l2_vi = verts[(node_index + 1) % verts_length],
               l3_vi = verts[(node_index + 2) % verts_length],
               l1_vi = verts[(node_index + 3) % verts_length];

            geom.GetNodePositionFromBV(S, l1_vi, al1_n);
            geom.GetNodePositionFromBV(S, l2_vi, al2_n);
            geom.GetNodePositionFromBV(S, l3_vi, al3_n);

            geom.GetNodeVelocity(S, l1_vi, Vl1);
            geom.GetNodeVelocity(S, l2_vi, Vl2);
            geom.GetNodeVelocity(S, l3_vi, Vl3);
            
            // cout << "cell: " << el_index << ", node: " << node << ", index: " << node_index << endl;
            // cout << "\tl1_i: " << l1_vi << endl;
            // cout << "\tl2_i: " << l2_vi << endl;
            // cout << "\tl3_i: " << l3_vi << endl;

            /* Compute nodal locations at next time step */
            add(1., al1_n, dt, Vl1, al1_np1);
            add(1., al2_n, dt, Vl2, al2_np1);
            add(1., al3_n, dt, Vl3, al3_np1);

            /* Compute Bj, c, and b vector */
            // Bj
            // cout << "Elem: " << el_index << endl;
            // cout << "Vol: " << pmesh->GetElementVolume(el_index);
            // cout << ", old sv: " << sv_old_gf.Elem(el_index);
            // cout << ", new sv: " << sv_gf.Elem(el_index) << endl;
            Bj = sv_gf.Elem(el_index) * pmesh->GetElementVolume(el_index) / sv_old_gf.Elem(el_index);

            // c and bvec
            subtract(al2_np1, al1_np1, temp_vec); // a_2^np1 - a_1^np1
            geom.Perpendicular(temp_vec);
            subtract(anode_n, al3_np1, temp_vec2); // a_r^n - a_3^np1

            // cout << "al1_n: ";
            // al1_n.Print(cout);
            // cout << "Vl1: ";
            // Vl1.Print(cout);
            // cout << "al1_np1: ";
            // al1_np1.Print(cout);

            // cout << "al2_n: ";
            // al2_n.Print(cout);
            // cout << "Vl2: ";
            // Vl2.Print(cout);
            // cout << "al2_np1: ";
            // al2_np1.Print(cout);

            // cout << "(al2_np1 - al1_np1)^perp: ";
            // temp_vec.Print(cout);
         
            cj = temp_vec * temp_vec2;
            cj *= 0.5;

            bj_vec = temp_vec;
            bj_vec *= dt / 2.;

            /* Add contribution to Mi and Ri */
            tensor(bj_vec, bj_vec, dm_temp);
            Mi.Add(1., dm_temp);

            Ri.Add(Bj - cj, bj_vec);
         } // End adjacent cells

         if (mm_visc_face != 0.)
         {
            visc_contr = 0.;
            visci_vec = 0.;
            /* Iterate over adjacent faces */
            for (int face_it = 0; face_it < faces_length; face_it++)
            {
               int face = faces_row[face_it];
               /* Get face normal and adjacent velocity */
               FI = pmesh->GetFaceInformation(face);

               H1.GetFaceDofs(face, face_dofs_row);
               int face_vdof1 = face_dofs_row[1],
                  face_vdof2 = face_dofs_row[0],
                  face_dof = face_dofs_row[2];

               // Grab corresponding vertex velocity from S
               if (node == face_vdof1) {
                  // Get adj index
                  Vadj_index = face_vdof2;
               } else {
                  // Get adj index
                  Vadj_index = face_vdof1;
               }
               geom.GetNodePositionFromBV(S, Vadj_index, aadj_n);
               geom.GetNodeVelocity(S, Vadj_index, Vadj_n);

               // Compute F
               subtract(anode_n, aadj_n, temp_vec);
               double F = temp_vec.Norml2();
               double F2 = pow(F,2);

               /* Collect contibution from viscosity */
               visc_contr += F2;
               visci_vec.Add(F2, Vadj_n);
            } // End adjacent faces iteration
            // Add contribution from viscosity to Mi and Ri
            for (int i = 0; i < dim; i++)
            {
               Mi.Elem(i,i) += visc_contr * mm_visc_face * pow(dt,2);
            }
            Ri.Add(mm_visc_face*pow(dt,2), visci_vec);
         } // End add viscosity

         /* Solve for Vnode */
         // cout << "Mi for node: " << node << endl;
         // Mi.Print(cout);
         Mi.Invert();
         Mi.Mult(Ri, predicted_node_v);

         /* Theta average with previous velocity */
         if (do_theta_averaging)
         {
            predicted_node_v *= theta; 
            predicted_node_v.Add(1. - theta, Vnode_n);
         }

         // Put velocity in S
         // Gauss-Seidel > Jacobi
         // geom.UpdateNodeVelocity(S, node, predicted_node_v);
         // Try Jacobi
         geom.UpdateNodeVelocity(S_new, node, predicted_node_v);
      } // End interior node

   } // End node iterator
   S = S_new;
}


/****************************************************************************************************
* Function: IterativeCornerVelocityLSCellVolumeCellVisc
* Parameters:
*  S    - BlockVector corresponding to nth timestep that stores mesh information, 
*         mesh velocity, and state variables.
*  dt   - Current timestep size
*
* Purpose:
*  This function constitutes one iteration on the previously computed corner node 
*  velocity to reduce total face motion while still preserving mass conservation.
*  The idea is to minimize the local change in mass for each adjacent cell to a 
*  given vertex, with some optional viscosity defined on all the adjacent cells.
*
*  This method was initially discussed on 06/26/2024.
*
*  Note: This function relies on the mesh velocities having already been computed
*  and possibly linearized.  One can use whichever mesh velocity computation before 
*  iteration.
*  Note: This function will modify mv_gf in the BlockVector S
****************************************************************************************************/
void LagrangianLOOperator::IterativeCornerVelocityLSCellVolumeCellVisc(Vector &S, const Vector &S_old, const double &dt)
{
   // cout << "=====IterativeCornerVelocityLSCellVolumeCellVisc=====\n";
   /* Optional run time parameters */
   bool do_theta_averaging = true;
   double theta = 0.5;

   /* Objects needed for function */
   mfem::Mesh::FaceInformation FI;
   Vector* sptr_old = const_cast<Vector*>(&S_old);
   ParGridFunction sv_gf, sv_old_gf;
   sv_old_gf.MakeRef(&L2, *sptr_old, block_offsets[1]);
   sv_gf.MakeRef(&L2, S, block_offsets[1]);

   Vector S_new = S;
   Vector predicted_node_v(dim);
   Vector anode_n(dim), al1_n(dim), al2_n(dim), al3_n(dim);
   Vector Vnode_n(dim), Vl1(dim), Vl2(dim), Vl3(dim);
   Vector al1_np1(dim), al2_np1(dim), al3_np1(dim);
   Vector temp_vec(dim), temp_vec2(dim), bj_vec(dim);

   double Bj, cj, Kcn;

   Array<int> cells_row, verts;
   int cells_length, bdr_ind;
   Vector Vcell_np1(dim), Uc(dim+2);

   DenseMatrix Mi(dim), dm_temp(dim);
   Vector Ri(dim), sum_vk(dim);
   double sum_k = 0.;

   /* Iterate over interior nodes */
   for (int node = 0; node < NDofs_H1L; node++)
   {
      /* Reset values for new node */
      Mi = 0., Ri = 0., predicted_node_v = 0.;
      sum_k = 0., sum_vk = 0.;
      bdr_ind = BdrVertexIndexingArray[node];

      /* bdr_ind = -1 for interior nodes */
      // if (bdr_ind == -1)
      if (bdr_ind != 5)
      {
         /* Get adjacent cells and faces */
         vertex_element->GetRow(node, cells_row);
         cells_length = cells_row.Size();

         /* Get nodal position */
         geom.GetNodePositionFromBV(S, node, anode_n);
         geom.GetNodeVelocity(S, node, Vnode_n);

         /* Iterate over adjacent cells */ 
         for (int cell_it = 0; cell_it < cells_length; cell_it++)
         {
            int el_index = cells_row[cell_it];
            Kcn = pmesh->GetElementVolume(el_index);

            /* Get remaining node locations and velocities */
            pmesh->GetElementVertices(el_index, verts);
            int verts_length = verts.Size();
            int node_index = verts.Find(node);
            int l2_vi = verts[(node_index + 1) % verts_length],
                  l3_vi = verts[(node_index + 2) % verts_length],
                  l1_vi = verts[(node_index + 3) % verts_length];

            geom.GetNodePositionFromBV(S, l1_vi, al1_n);
            geom.GetNodePositionFromBV(S, l2_vi, al2_n);
            geom.GetNodePositionFromBV(S, l3_vi, al3_n);

            geom.GetNodeVelocity(S, l1_vi, Vl1);
            geom.GetNodeVelocity(S, l2_vi, Vl2);
            geom.GetNodeVelocity(S, l3_vi, Vl3);

            /* Compute nodal locations at next time step */
            add(1., al1_n, dt, Vl1, al1_np1);
            add(1., al2_n, dt, Vl2, al2_np1);
            add(1., al3_n, dt, Vl3, al3_np1);

            /* Compute Bj, c, and b vector */
            // Bj
            Bj = sv_gf.Elem(el_index) * Kcn / sv_old_gf.Elem(el_index);

            // c and bvec
            subtract(al2_np1, al1_np1, temp_vec); // a_2^np1 - a_1^np1
            geom.Perpendicular(temp_vec);
            subtract(anode_n, al3_np1, temp_vec2); // a_r^n - a_3^np1
         
            cj = temp_vec * temp_vec2;
            cj *= 0.5;

            bj_vec = temp_vec;
            bj_vec *= dt / 2.;

            /* Add contribution to Mi and Ri */
            tensor(bj_vec, bj_vec, dm_temp);
            // cout << "bj_vec: ";
            // bj_vec.Print(cout);
            Mi.Add(1., dm_temp);

            Ri.Add(Bj - cj, bj_vec);

            // Since the viscosity is defined according to adjacent cells, 
            // we compute its contribution in the cell loop
            if (mm_visc_face != 0.)
            {
               // Retrieve cell velocity
               GetCellStateVector(S, el_index, Uc);
               pb->velocity(Uc, Vcell_np1);
               // cout << "cell: " << el_index << "vel: ";
               // Vcell_np1.Print(cout);
               sum_k += Kcn;
               sum_vk.Add(1., Vcell_np1);
            } // End add viscosity contribution from cell
         } // End adjacent cells

         // Add in contribution of viscosity
         for (int i = 0; i < dim; i++)
         {
            Mi.Elem(i,i) += mm_visc_face * pow(dt,2) * sum_k * pow(cells_length, 2);
         }
         Ri.Add(mm_visc_face*pow(dt,2)*sum_k*cells_length, sum_vk);

         // Solve
         Mi.Invert();
         Mi.Mult(Ri, predicted_node_v);

         /* Theta average with previous velocity */
         if (do_theta_averaging)
         {
            predicted_node_v *= theta; 
            predicted_node_v.Add(1. - theta, Vnode_n);
         }

         // Put velocity in S
         // Gauss-Seidel > Jacobi
         // geom.UpdateNodeVelocity(S, node, predicted_node_v);
         // Try Jacobi
         geom.UpdateNodeVelocity(S_new, node, predicted_node_v);
      } // End interior node

   } // End node iterator
   S = S_new;
}


/****************************************************************************************************
* Function: IterativeCornerVelocityLSCellVolumeMv2Visc
* Parameters:
*  S     - BlockVector corresponding to (n+1)th timestep that stores mesh information, 
*          mesh velocity, and state variables.
*  S_old - BlockVector corresponding to nth timestep that stores mesh information, 
*          mesh velocity, and state variables.
*  dt    - Current timestep size
*
* Purpose:
*  This function constitutes one iteration on the previously computed corner node 
*  velocity to reduce total face motion while still preserving mass conservation.
*  The idea is to minimize the local change in mass for each adjacent cell to a 
*  given vertex, with some optional viscosity defined on all the adjacent cells.
*
*  This method was initially discussed on 06/28/2024.
*
*  Note: This function relies on the mesh velocities having already been computed
*  and possibly linearized.  One can use whichever mesh velocity computation before 
*  iteration.
*  Note: This function will modify mv_gf in the BlockVector S
****************************************************************************************************/
void LagrangianLOOperator::IterativeCornerVelocityLSCellVolumeMv2Visc(Vector &S, const Vector &S_old, const ParGridFunction &mv2_gf, const double &dt)
{
   // cout << "=====IterativeCornerVelocityLSCellVolumeMv2Visc=====\n";
   /* Optional run time parameters */
   bool do_theta_averaging = true;
   double theta = 0.5;

   /* Objects needed for function */
   mfem::Mesh::FaceInformation FI;
   Vector* sptr_old = const_cast<Vector*>(&S_old);
   ParGridFunction sv_gf, sv_old_gf;
   sv_old_gf.MakeRef(&L2, *sptr_old, block_offsets[1]);
   sv_gf.MakeRef(&L2, S, block_offsets[1]);

   Vector S_new = S;
   Vector predicted_node_v(dim), Vnode_mv2(dim);
   Vector anode_n(dim), al1_n(dim), al2_n(dim), al3_n(dim);
   Vector Vnode_n(dim), Vl1(dim), Vl2(dim), Vl3(dim);
   Vector al1_np1(dim), al2_np1(dim), al3_np1(dim);
   Vector temp_vec(dim), temp_vec2(dim), bj_vec(dim);

   double Bj, cj, Kcn;

   Array<int> cells_row, verts;
   int cells_length, bdr_ind;

   DenseMatrix Mi(dim), dm_temp(dim);
   Vector Ri(dim);
   double sum_k = 0.;

   /* Iterate over interior nodes */
   for (int node = 0; node < NDofs_H1L; node++)
   {
      /* Reset values for new node */
      Mi = 0., Ri = 0., predicted_node_v = 0.;
      sum_k = 0.;
      bdr_ind = BdrVertexIndexingArray[node];

      /* bdr_ind = -1 for interior nodes */
      // if (bdr_ind == -1)
      if (bdr_ind != 5)
      {
         /* Get adjacent cells and faces */
         vertex_element->GetRow(node, cells_row);
         cells_length = cells_row.Size();

         /* Get nodal position */
         geom.GetNodePositionFromBV(S, node, anode_n);
         geom.GetNodeVelocity(S, node, Vnode_n);

         /* Iterate over adjacent cells */ 
         for (int cell_it = 0; cell_it < cells_length; cell_it++)
         {
            int el_index = cells_row[cell_it];
            Kcn = pmesh->GetElementVolume(el_index);

            /* Get remaining node locations and velocities */
            pmesh->GetElementVertices(el_index, verts);
            int verts_length = verts.Size();
            int node_index = verts.Find(node);
            int l2_vi = verts[(node_index + 1) % verts_length],
                  l3_vi = verts[(node_index + 2) % verts_length],
                  l1_vi = verts[(node_index + 3) % verts_length];

            geom.GetNodePositionFromBV(S, l1_vi, al1_n);
            geom.GetNodePositionFromBV(S, l2_vi, al2_n);
            geom.GetNodePositionFromBV(S, l3_vi, al3_n);

            geom.GetNodeVelocity(S, l1_vi, Vl1);
            geom.GetNodeVelocity(S, l2_vi, Vl2);
            geom.GetNodeVelocity(S, l3_vi, Vl3);

            /* Compute nodal locations at next time step */
            add(1., al1_n, dt, Vl1, al1_np1);
            add(1., al2_n, dt, Vl2, al2_np1);
            add(1., al3_n, dt, Vl3, al3_np1);

            /* Compute Bj, c, and b vector */
            // Bj
            Bj = sv_gf.Elem(el_index) * Kcn / sv_old_gf.Elem(el_index);

            // c and bvec
            subtract(al2_np1, al1_np1, temp_vec); // a_2^np1 - a_1^np1
            geom.Perpendicular(temp_vec);
            subtract(anode_n, al3_np1, temp_vec2); // a_r^n - a_3^np1
         
            cj = temp_vec * temp_vec2;
            cj *= 0.5;

            bj_vec = temp_vec;
            bj_vec *= dt / 2.;

            /* Add contribution to Mi and Ri */
            tensor(bj_vec, bj_vec, dm_temp);
            // cout << "bj_vec: ";
            // bj_vec.Print(cout);
            Mi.Add(1., dm_temp);

            Ri.Add(Bj - cj, bj_vec);

            // Add up cell volumes for viscosity contribution
            if (mm_cell != 0.)
            {
               sum_k += Kcn;
            } 
         } // End adjacent cells

         // Add in contribution of viscosity
         if (mm_cell != 0.)
         {
            // Get mv2 velocity for node i 
            geom.GetNodeVelocity(mv2_gf, node, Vnode_mv2);
            for (int i = 0; i < dim; i++)
            {
               Mi.Elem(i,i) += mm_cell * pow(dt,2) * sum_k;
            }
            Ri.Add(mm_cell*pow(dt,2)*sum_k, Vnode_mv2);
         }

         // Solve
         Mi.Invert();
         Mi.Mult(Ri, predicted_node_v);

         /* Theta average with previous velocity */
         if (do_theta_averaging)
         {
            predicted_node_v *= theta; 
            predicted_node_v.Add(1. - theta, Vnode_n);
         }

         // Put velocity in S
         // Gauss-Seidel > Jacobi
         // geom.UpdateNodeVelocity(S, node, predicted_node_v);
         // Try Jacobi
         geom.UpdateNodeVelocity(S_new, node, predicted_node_v);
      } // End interior node

   } // End node iterator
   S = S_new;
}


/****************************************************************************************************
* Function: IterativeCornerVelocityLSCellVolumeMv2FaceVisc
* Parameters:
*  S      - BlockVector corresponding to (n+1)th timestep that stores mesh information, 
*           mesh velocity, and state variables.
*  S_old  - BlockVector corresponding to nth timestep that stores mesh information, 
*           mesh velocity, and state variables.
*  mv2_gf - 
*  dt     - Current timestep size
*
* Purpose:
*  This function constitutes one iteration on the previously computed corner node 
*  velocity to reduce total face motion while still preserving mass conservation.
*  The idea is to minimize the local change in mass for each adjacent cell to a 
*  given vertex, with some optional viscosity defined on all the adjacent cells.
*
*  This method was initially discussed on 06/28/2024. Face flux and mv2 flux with
*  different omega values.
*
*  Note: This function relies on the mesh velocities having already been computed
*  and possibly linearized.  One can use whichever mesh velocity computation before 
*  iteration.
*  Note: This function will modify mv_gf in the BlockVector S
****************************************************************************************************/
void LagrangianLOOperator::IterativeCornerVelocityLSCellVolumeMv2FaceVisc(Vector &S, const Vector &S_old, const ParGridFunction &mv2_gf, const double &dt)
{
   // cout << "=====IterativeCornerVelocityLSCellVolumeMv2FaceVisc=====\n";
   /* Optional run time parameters */
   bool do_theta_averaging = true;
   double theta = 0.5;

   /* Objects needed for function */
   mfem::Mesh::FaceInformation FI;
   Vector* sptr_old = const_cast<Vector*>(&S_old);
   ParGridFunction sv_gf, sv_old_gf;
   sv_old_gf.MakeRef(&L2, *sptr_old, block_offsets[1]);
   sv_gf.MakeRef(&L2, S, block_offsets[1]);

   Vector S_new = S;
   Vector predicted_node_v(dim), Vnode_mv2(dim);
   Vector anode_n(dim), al1_n(dim), al2_n(dim), al3_n(dim);
   Vector Vnode_n(dim), Vl1(dim), Vl2(dim), Vl3(dim);
   Vector al1_np1(dim), al2_np1(dim), al3_np1(dim);
   Vector temp_vec(dim), temp_vec2(dim), bj_vec(dim);

   double Bj, cj, Kcn;

   Array<int> cells_row, verts;
   int cells_length, bdr_ind;

   DenseMatrix Mi(dim), dm_temp(dim);
   Vector Ri(dim);
   double sum_k = 0.;

   /* Iterate over all nodes except corners */
   for (int node = 0; node < NDofs_H1L; node++)
   {
      /* Reset values for new node */
      Mi = 0., Ri = 0., predicted_node_v = 0.;
      sum_k = 0.;
      bdr_ind = BdrVertexIndexingArray[node];

      /* bdr_ind = -1 for interior nodes */
      // if (bdr_ind == -1)
      if (bdr_ind != 5)
      {
         /* Get adjacent cells */
         vertex_element->GetRow(node, cells_row);
         cells_length = cells_row.Size();

         /* Get nodal position */
         geom.GetNodePositionFromBV(S, node, anode_n);
         geom.GetNodeVelocity(S, node, Vnode_n);

         /* Get mv2 velocity if needed in cell visc contribution */
         geom.GetNodeVelocity(mv2_gf, node, Vnode_mv2);

         /* Iterate over adjacent cells */ 
         for (int cell_it = 0; cell_it < cells_length; cell_it++)
         {
            int el_index = cells_row[cell_it];
            Kcn = pmesh->GetElementVolume(el_index);

            /* Get remaining node locations and velocities */
            pmesh->GetElementVertices(el_index, verts);
            int verts_length = verts.Size();
            int node_index = verts.Find(node);
            int l2_vi = verts[(node_index + 1) % verts_length],
                  l3_vi = verts[(node_index + 2) % verts_length],
                  l1_vi = verts[(node_index + 3) % verts_length];

            geom.GetNodePositionFromBV(S, l1_vi, al1_n);
            geom.GetNodePositionFromBV(S, l2_vi, al2_n);
            geom.GetNodePositionFromBV(S, l3_vi, al3_n);

            geom.GetNodeVelocity(S, l1_vi, Vl1);
            geom.GetNodeVelocity(S, l2_vi, Vl2);
            geom.GetNodeVelocity(S, l3_vi, Vl3);

            /* Compute nodal locations at next time step */
            add(1., al1_n, dt, Vl1, al1_np1);
            add(1., al2_n, dt, Vl2, al2_np1);
            add(1., al3_n, dt, Vl3, al3_np1);

            /* Compute Bj, c, and b vector */
            // Bj
            Bj = sv_gf.Elem(el_index) * Kcn / sv_old_gf.Elem(el_index);

            // c and bvec
            subtract(al2_np1, al1_np1, temp_vec); // a_2^np1 - a_1^np1
            geom.Perpendicular(temp_vec);
            subtract(anode_n, al3_np1, temp_vec2); // a_r^n - a_3^np1
         
            cj = temp_vec * temp_vec2;
            cj *= 0.5;

            bj_vec = temp_vec;
            bj_vec *= dt / 2.;

            /* Add contribution to Mi and Ri */
            tensor(bj_vec, bj_vec, dm_temp);
            Mi.Add(1., dm_temp);
            Ri.Add(Bj - cj, bj_vec);

            // Add up cell volumes for cell viscosity contribution
            if (mm_cell != 0.)
            {
               sum_k += Kcn;
            } 
         } // End adjacent cells

         // Add in contribution from cell viscosity
         if (mm_cell != 0.)
         {
            for (int i = 0; i < dim; i++)
            {
               Mi.Elem(i,i) += mm_cell * pow(dt,2) * sum_k;
            }
            Ri.Add(mm_cell*pow(dt,2)*sum_k, Vnode_mv2);
         }

         // Add in contribution of viscosity on faces
         if (mm_visc_face != 0.)
         {
            double sum_F2 = 0.;
            int Vadj_index, faces_length;
            Vector aadj_n(dim), Vadj_n(dim), visci_vec(dim);
            visci_vec = 0.;

            /* Get adjacent faces */
            Array<int> faces_row, face_dofs_row;
            vertex_edge.GetRow(node, faces_row);
            faces_length = faces_row.Size();

            /* Iterate over adjacent faces */
            for (int face_it = 0; face_it < faces_length; face_it++)
            {
               int face = faces_row[face_it];
               /* Get face normal and adjacent velocity */
               FI = pmesh->GetFaceInformation(face);

               H1.GetFaceDofs(face, face_dofs_row);
               int face_vdof1 = face_dofs_row[1],
                  face_vdof2 = face_dofs_row[0],
                  face_dof = face_dofs_row[2];

               // Grab corresponding adjacent vertex velocity from S
               if (node == face_vdof1) {
                  // Get adj index
                  Vadj_index = face_vdof2;
               } else {
                  // Get adj index
                  Vadj_index = face_vdof1;
               }
               geom.GetNodePositionFromBV(S, Vadj_index, aadj_n);
               geom.GetNodeVelocity(S, Vadj_index, Vadj_n);

               // Compute F
               subtract(anode_n, aadj_n, temp_vec);
               double F = temp_vec.Norml2();
               double F2 = pow(F,2);

               /* Collect contibution from face viscosity */
               sum_F2 += F2;
               visci_vec.Add(F2, Vadj_n);
            } // End adjacent faces iteration

            // Add contribution from face viscosity to Mi and Ri
            for (int i = 0; i < dim; i++)
            {
               Mi.Elem(i,i) += sum_F2 * mm_visc_face * pow(dt,2);
            }
            Ri.Add(mm_visc_face*pow(dt,2), visci_vec);
         }

         /* Solve */
         Mi.Invert();
         Mi.Mult(Ri, predicted_node_v);

         /* Theta average with previous velocity */
         if (do_theta_averaging)
         {
            predicted_node_v *= theta; 
            predicted_node_v.Add(1. - theta, Vnode_n);
         }

         // Put velocity in S
         // Gauss-Seidel > Jacobi
         // geom.UpdateNodeVelocity(S, node, predicted_node_v);
         // Try Jacobi
         geom.UpdateNodeVelocity(S_new, node, predicted_node_v);
      } // End interior node

   } // End node iterator
   S = S_new;
}


/****************************************************************************************************
* Function: VerifyContributions
* Parameters:
*  S      - BlockVector corresponding to (n+1)th timestep that stores mesh information, 
*           mesh velocity, and state variables.
*  S_old  - BlockVector corresponding to nth timestep that stores mesh information, 
*           mesh velocity, and state variables.
*  mv2_gf - 
*  dt     - Current timestep size
*
* Purpose:
*  This function is to be run following the function
*           IterativeCornerVelocityLSCellVolumeMv2FaceVisc
*  to compare the contribution
****************************************************************************************************/
void LagrangianLOOperator::VerifyContributions(const Vector &S, const Vector &S_old, const ParGridFunction &mv2_gf, const double &dt, const int &it)
{
   MFEM_ABORT("Function deprecated\n");
   Vector* sptr_old = const_cast<Vector*>(&S_old);
   double gamma_numer = 0., gamma_c_numer = 0., gamma_f_numer = 0., denom = 0.;

   ParGridFunction sv_gf, sv_old_gf, mv_gf;
   sv_old_gf.MakeRef(&L2, *sptr_old, block_offsets[1]);

   // Move vertices of S_temp to compute cell volume norm
   Vector S_temp = S;
   Vector *stemp_ptr = const_cast<Vector*>(&S_temp);
   x_gf.MakeRef(&H1, *stemp_ptr, block_offsets[0]);
   MFEM_ABORT("Need to remove mv_gf implementation\n");
   mv_gf.MakeRef(&H1, *stemp_ptr, block_offsets[1]);
   sv_gf.MakeRef(&L2, *stemp_ptr, block_offsets[1]);
   add(x_gf, dt, mv_gf, x_gf);

   Array<int> cells_row;
   double sum_k = 0., Kcn = 0.;
   Vector Uc(dim+2), Vnode(dim), Vnode_mv2(dim), temp_vec(dim);
   int bdr_ind, cells_length;

   Vector V1(dim), V2(dim);
   Vector a1(dim), a2(dim);
   mfem::Mesh::FaceInformation FI;
   double F;
   Array<int> face_dofs_row;

   // Sum over all cells
   for (int ci = 0; ci < NDofs_L2; ci++)
   {
      if ((pb->get_indicator() == "ElasticNoh" || pb->get_indicator() == "Noh") && cell_bdr_flag_gf[ci] != -1)
      {
         // Skip boundary cells
         continue;
      }

      double Kn = ComputeCellVolume(S_old, ci);
      double Knp1 = ComputeCellVolume(S_temp, ci);
      double Tn = sv_old_gf.Elem(ci);
      double Tnp1 = sv_gf.Elem(ci);

      double tmp_val = pow((Kn / Tn) - (Knp1 / Tnp1), 2);
      gamma_numer += tmp_val;
      denom += pow(Kn/Tn, 2);
   }

   // Iterate over all geometric nodes except corners
   for (int node = 0; node < NDofs_H1L; node++)
   {
      sum_k = 0.;
      bdr_ind = BdrVertexIndexingArray[node];

      // All nodes minus corner nodes will contribute
      if (bdr_ind != 5)
      {
         vertex_element->GetRow(node, cells_row);
         cells_length = cells_row.Size();

         geom.GetNodeVelocity(S, node, Vnode);
         geom.GetNodeVelocity(mv2_gf, node, Vnode_mv2);

         for (int cell_it = 0; cell_it < cells_length; cell_it++)
         {
            int el_index = cells_row[cell_it];
            Kcn = pmesh->GetElementVolume(el_index);

            sum_k += Kcn;
         }
         subtract(Vnode, Vnode_mv2, temp_vec);
         gamma_c_numer += sum_k * pow(temp_vec.Norml2(), 2);
      }
   }

   // Iterate over all faces
   for (int face = 0; face < num_faces; face++)
   {
      FI = pmesh->GetFaceInformation(face);
      H1.GetFaceDofs(face,face_dofs_row);
      int face_dof1 = face_dofs_row[0], face_dof2 = face_dofs_row[1];

      geom.GetNodePositionFromBV(S, face_dof1, a1);
      geom.GetNodePositionFromBV(S, face_dof2, a2);
      geom.GetNodeVelocity(S, face_dof1, V1);
      geom.GetNodeVelocity(S, face_dof2, V2);

      subtract(a2, a1, temp_vec);
      F = temp_vec.Norml2();

      subtract(V2, V1, temp_vec);
      gamma_f_numer += pow(F,2) * pow(temp_vec.Norml2(), 2);
   }

   double gamma = sqrt(gamma_numer / denom);
   double gamma_c = sqrt(mm_cell * pow(dt,2) * gamma_c_numer / denom);
   double gamma_f = sqrt(mm_visc_face * pow(dt,2) * gamma_f_numer / denom);
   cout << "it: " << it 
        <<  ", minimized quantity: " << gamma
        << ", viscosity on cells: " << gamma_c 
        << ", viscosity on faces: " << gamma_f << endl;
}


/****************************************************************************************************
* Function: ComputeCellVolume
* Parameters:
*  S        - BlockVector representing FiniteElement information
*  cell     - Cell to comput volume at
*
* Purpose: 
*  To compute the volume of a given cell according to Despres paper.  This is a Q1
*  computation and assumes the faces do not move for conservation.
****************************************************************************************************/
double LagrangianLOOperator::ComputeCellVolume(const Vector &S, const int &cell)
{
   Array<int> verts;
   switch (dim)
   {
      case 1:
      {
         MFEM_ABORT("Dimension 1 not implemented.\n");
         return 0.;
      }
      case 2:
      {
         // Get element nodes and positions
         pmesh->GetElementVertices(cell, verts);
         int verts_length = verts.Size();
         Vector a0(dim), a1(dim), a2(dim), a3(dim);
         geom.GetNodePositionFromBV(S, verts[0], a0);
         geom.GetNodePositionFromBV(S, verts[1], a1);
         geom.GetNodePositionFromBV(S, verts[2], a2);
         geom.GetNodePositionFromBV(S, verts[3], a3);

         // Compute area of square
         Vector temp_vec1(dim), temp_vec2(dim);
         subtract(a3, a1, temp_vec1);
         geom.Perpendicular(temp_vec1);
         subtract(a2, a0, temp_vec2);
         double area = temp_vec1 * temp_vec2;
         area *= 0.5;
         return area;
      }
      case 3:
      {
         MFEM_ABORT("Dimension 3 not implemented.\n");
         return 0.;
      }
      default:
      {
         MFEM_ABORT("Invalid dimension value provided.\n");
         return 0.;
      }
   }
}


/****************************************************************************************************
* Function: ComputeCellVolumeNorm
* Parameters:
*  S        - BlockVector representing FiniteElement information at tn+1
*  S_old    - BlockVector representing FiniteElement information at tn
*  dt       - Timestep
*
* Purpose: 
*  Compute the sum of the change in mass over all cells in the mesh.
****************************************************************************************************/
double LagrangianLOOperator::ComputeCellVolumeNorm(const Vector &S, const Vector &S_old, const double &dt)
{
   MFEM_ABORT("Function deprecated\n");
   Vector* sptr_old = const_cast<Vector*>(&S_old);

   ParGridFunction sv_gf, sv_old_gf, mv_gf;
   sv_old_gf.MakeRef(&L2, *sptr_old, block_offsets[1]);

   // Move vertices of S_temp to compute cell volume norm
   Vector S_temp = S;
   Vector *stemp_ptr = const_cast<Vector*>(&S_temp);
   x_gf.MakeRef(&H1, *stemp_ptr, block_offsets[0]);
   MFEM_ABORT("Need to remove mv_gf implementation\n");
   mv_gf.MakeRef(&H1, *stemp_ptr, block_offsets[1]);
   sv_gf.MakeRef(&L2, *stemp_ptr, block_offsets[1]);
   add(x_gf, dt, mv_gf, x_gf);

   double num = 0., denom = 0.;
   double cell_val = 0., face_val = 0.;

   // Sum over all cells
   for (int ci = 0; ci < NDofs_L2; ci++)
   {
      if ((pb->get_indicator() == "ElasticNoh" || pb->get_indicator() == "Noh") && cell_bdr_flag_gf[ci] != -1)
      {
         // Skip boundary cells
         continue;
      }

      double Kn = ComputeCellVolume(S_old, ci);
      double Knp1 = ComputeCellVolume(S_temp, ci);
      double Tn = sv_old_gf.Elem(ci);
      double Tnp1 = sv_gf.Elem(ci);

      double tmp_val = pow((Kn / Tn) - (Knp1 / Tnp1), 2);
      num += tmp_val;
      denom += pow(Kn / Tn, 2);
   }

   return sqrt(num/denom);
}


/****************************************************************************************************
* Function: compare_gamma
* Parameters:
*  S        - BlockVector representing FiniteElement information at tn+1
*  S_old    - BlockVector representing FiniteElement information at tn
*  dt       - Timestep
*  it       - Iteration number
*
* Purpose: 
*  Compare the value being minimized and the amount of viscosity being added
*  based on the average velocity of the adjacent cells.
* 
*  NOTE: We must use Jacobi method in stead of Gauss-Seidel.
****************************************************************************************************/
void LagrangianLOOperator::compare_gamma2(const Vector &S, const Vector &S_old, const double &dt, const int &it)
{
   MFEM_ABORT("Function deprecated\n");
   Vector* sptr_old = const_cast<Vector*>(&S_old);
   double alphaF = 0., alphaV = 0., denom = 0.;

   ParGridFunction sv_gf, sv_old_gf, mv_gf;
   sv_old_gf.MakeRef(&L2, *sptr_old, block_offsets[1]);

   // Move vertices of S_temp to compute cell volume norm
   Vector S_temp = S;
   Vector *stemp_ptr = const_cast<Vector*>(&S_temp);
   x_gf.MakeRef(&H1, *stemp_ptr, block_offsets[0]);
   MFEM_ABORT("Need to remove mv_gf implementation\n");
   mv_gf.MakeRef(&H1, *stemp_ptr, block_offsets[1]);
   sv_gf.MakeRef(&L2, *stemp_ptr, block_offsets[1]);
   add(x_gf, dt, mv_gf, x_gf);

   Array<int> cells_row;
   double sum_k = 0., Kcn = 0.;
   Vector Uc(dim+2), Vcell(dim), Vnode(dim), sum_vk(dim), temp_vec(dim);
   int bdr_ind, cells_length;

   // Sum over all cells
   for (int ci = 0; ci < NDofs_L2; ci++)
   {
      if ((pb->get_indicator() == "ElasticNoh" || pb->get_indicator() == "Noh") && cell_bdr_flag_gf[ci] != -1)
      {
         // Skip boundary cells
         continue;
      }

      double Kn = ComputeCellVolume(S_old, ci);
      double Knp1 = ComputeCellVolume(S_temp, ci);
      double Tn = sv_old_gf.Elem(ci);
      double Tnp1 = sv_gf.Elem(ci);

      double tmp_val = pow((Kn / Tn) - (Knp1 / Tnp1), 2);
      alphaF += tmp_val;
      denom += pow(Kn/Tn, 2);
   }

   // Iterate over all geometric nodes
   for (int node = 0; node < NDofs_H1L; node++)
   {
      sum_k = 0.;
      sum_vk = 0.;
      bdr_ind = BdrVertexIndexingArray[node];

      // All nodes minus corner nodes will contribute
      if (bdr_ind != 5)
      {
         vertex_element->GetRow(node, cells_row);
         cells_length = cells_row.Size();

         geom.GetNodeVelocity(S, node, Vnode);

         for (int cell_it = 0; cell_it < cells_length; cell_it++)
         {
            int el_index = cells_row[cell_it];
            Kcn = pmesh->GetElementVolume(el_index);

            GetCellStateVector(S, el_index, Uc);
            pb->velocity(Uc, Vcell);

            sum_k += Kcn;
            sum_vk.Add(1, Vcell);
         }
         add(cells_length, Vnode, -1., sum_vk, temp_vec);
         alphaV += sum_k * pow(temp_vec.Norml2(), 2);
      }
   }

   cout << "it: " << it 
        <<  ", minimized quantity: " << sqrt(alphaF/denom) 
        << ", viscosity added: " << sqrt(mm_visc_face*pow(dt,2)*alphaV/denom) << endl;
}


/****************************************************************************************************
* Function: ComputeAverageVelocities
* Parameters:
*  S        - BlockVector representing FiniteElement information
*
* Purpose: Compute average of adjacent corner velocities for interior nodes and average this with 
*          current velocity at each node.  This addition restriction on the nodal velocities was 
*          discussed on 06/03/2024.
* 
*  NOTE: We must use Jacobi method in stead of Gauss-Seidel.
****************************************************************************************************/
void LagrangianLOOperator::ComputeAverageVelocities(Vector &S)
{
   cout << "=====ComputeAverageVelocities=====\n";
   double theta = 0.01;
   Vector S_new = S;
   Array<int> faces_row, face_dofs_row;
   int faces_length;
   mfem::Mesh::FaceInformation FI;
   Vector Vnode(dim), Vadj(dim), Vpred(dim);
   Vector n_int(dim), n_vec(dim);
   int Vadj_index;

   for (int node = 0; node < NDofs_H1L; node++)
   {
      Vpred = 0.;
      geom.GetNodeVelocity(S, node, Vnode);

      // Only average for interior nodes
      int bdr_ind = BdrVertexIndexingArray[node];
      if (bdr_ind == -1)
      {
         // Get cell faces
         vertex_edge.GetRow(node, faces_row);
         faces_length = faces_row.Size();
         assert(faces_length == 4);
         /* Iterate over cell faces */
         for (int face_it = 0; face_it < faces_length; face_it++) // Adjacent face iterator
         {
            // Get face information
            int face = faces_row[face_it];
            FI = pmesh->GetFaceInformation(face);

            // Weight with F
            int c = FI.element[0].index;
            CalcOutwardNormalInt(S, c, face, n_int);
            n_vec = n_int;
            double F = n_vec.Norml2();

            /* adjacent corner indices */
            // preserve node orientation where cell to right 
            // of face is with lower cell index
            H1.GetFaceDofs(face, face_dofs_row);
            int face_vdof1 = face_dofs_row[1], 
               face_vdof2 = face_dofs_row[0], 
               face_dof = face_dofs_row[2]; 

            // Grab corresponding vertex velocity from S
            // We are always solving for V2. 
            // If the index for Vnode does not match that
            // for V2, we must flip the normal.
            if (node == face_vdof1) {
               // Get adj index
               Vadj_index = face_vdof2;
            } else {
               // Get adj index
               Vadj_index = face_vdof1;
            }

            /* Get face and adj node velocity and position */
            geom.GetNodeVelocity(S, Vadj_index, Vadj);

            // Vpred += Vadj;
            Vpred.Add(F, Vadj);
         }

         // Average 
         Vpred /=  faces_length;

         // Add in viscosity
         add(1. - theta, Vnode, theta, Vpred, Vnode);
         geom.UpdateNodeVelocity(S_new, node, Vnode);
      } // End interior node
   } // End node iteration
   S = S_new;
}


/*
*
*/
void LagrangianLOOperator::IterativeCornerVelocityFLUXLS(Vector &S, const double & dt)
{
   // cout << "=====IterativeCornerVelocityFLUXLS=====\n";
   // Optional run time parameters
   bool do_theta_averaging = false;
   double theta = 0.5;
   Vector S_new = S;

   // Values needed during iteration
   Array<int> faces_row, face_dofs_row;
   int faces_length;
   mfem::Mesh::FaceInformation FI;
   H1.ExchangeFaceNbrData();
   int Vadj_index, c, bdr_ind = 0;
   double F = 0.;
   double Vnode_prev_it_nR_comp = 0., Vadj_n_comp = 0., Vadj_nR_comp = 0.;

   Vector predicted_node_v(dim), Vf(dim), Vnode_np1(dim), Vnode_prev_it(dim);
   Vector n_int(dim), n_vec(dim), n_vec_R(dim);
   Vector Vnode_x(dim), Vadj_x(dim), face_x(dim);
   Vector Vnode_xnp1(dim), Vadj_xnp1(dim), face_xnp1(dim), a12_xnp1(dim);
   Vector Badj(dim), Bnode(dim);

   Vector Vnode_n(dim), Vadj_n(dim), Vface_n(dim);
   Vector Vnode_half(dim), Vadj_half(dim);
   Vector temp_vec(dim), temp_vec_2(dim);
   double D = 0., c1 = 0., c0 = 0., V3nperp = 0., Dc1 = 0.;

   // Averaging
   DenseMatrix Mi(dim), dm_tmp(dim);
   Vector Ri(dim);

   /* Iterate over corner nodes */
   for (int node = 0; node < NVDofs_H1; node++) // TODO: Is NDofs_H1L == NVDofs_H1?
   {
      // Reset new nodal velocity, and averaging objects
      Mi = 0., Ri = 0.;
      predicted_node_v = 0.;

      // Set corresponding boundary indicator
      // If node is boundary node, will have to add corresponding correction
      // node is interior node if bdr_ind == -1
      bdr_ind = BdrVertexIndexingArray[node];

      // Get current nodal velocity and position 
      geom.GetNodeVelocity(S, node, Vnode_n);
      geom.GetNodePositionFromBV(S, node, Vnode_x);

      // Compute An and an+1 for node
      add(Vnode_x, dt/2., Vnode_n, Vnode_half);
      add(Vnode_x, dt, Vnode_n, Vnode_xnp1);

      // Get cell faces
      vertex_edge.GetRow(node, faces_row);
      faces_length = faces_row.Size();

      /* Iterate over cell faces */
      for (int face_it = 0; face_it < faces_length; face_it++) // Adjacent face iterator
      {
         // Get face information
         int face = faces_row[face_it];
         GetIntermediateFaceVelocity(face, Vf);
         FI = pmesh->GetFaceInformation(face);

         // Calculate outer normal
         c = FI.element[0].index;
         CalcOutwardNormalInt(S, c, face, n_int);
         n_vec = n_int;
         F = n_vec.Norml2();
         n_vec /= F;

         /* adjacent corner indices */
         // preserve node orientation where cell to right 
         // of face is with lower cell index
         H1.GetFaceDofs(face, face_dofs_row);
         int face_vdof1 = face_dofs_row[1], 
             face_vdof2 = face_dofs_row[0], 
             face_dof = face_dofs_row[2]; 

         // Grab corresponding vertex velocity from S
         // We are always solving for V2. 
         // If the index for Vnode does not match that
         // for V2, we must flip the normal.
         if (node == face_vdof1) {
            // Get adj index
            Vadj_index = face_vdof2;
            // Flip normal vector
            n_vec *= -1.;
         } else {
            // Get adj index
            Vadj_index = face_vdof1;
            // Node matches v2
         }

         /* Get face and adj node velocity and position */
         geom.GetNodeVelocity(S, Vadj_index, Vadj_n);
         geom.GetNodeVelocity(S, face_dof, Vface_n);
         geom.GetNodePositionFromBV(S, Vadj_index, Vadj_x);
         geom.GetNodePositionFromBV(S, face_dof, face_x);

         // Check to make sure orientation is correct
         // This part is mostly for assertion
         Vector ttvec(dim);
         subtract(Vnode_x, Vadj_x, ttvec);
         geom.Orthogonal(ttvec);
         double _val = ttvec * n_vec;
         // Vector should be pointing in same direction
         if (_val < 1E-6)
         {
            cout << "interior analysis:\n";
            cout << "node: " << node << endl;
            cout << "nodex: ";
            Vnode_x.Print();
            cout << "Vadj_x: ";
            Vadj_x.Print(cout);
            cout << "nvec: ";
            n_vec.Print(cout);
            assert(false);
         }

         // Compute future face and adjacent node locations
         add(Vadj_x, dt/2., Vadj_n, Vadj_half);
         add(Vadj_x, dt, Vadj_n, Vadj_xnp1);
         // add(0.5, Vadj_xnp1, 0.5, Vnode_xnp1, a12_xnp1);

         // Perpendicular vector
         n_vec_R = n_vec;
         geom.Orthogonal(n_vec_R);

         // Calculate bmn
         double bmn = Vf * n_vec;
         bmn *= F;

         // Get Vnode_prev_it component in tangent direction from previous iteration
         Vnode_prev_it_nR_comp = Vnode_n * n_vec_R;

         // Get normal and rotated components of Vadj_n
         Vadj_n_comp = Vadj_n * n_vec;
         Vadj_nR_comp = Vadj_n * n_vec_R;
         
         // Compute geometrical vectors
         Bnode = face_x;
         Bnode *= -2.;
         Bnode.Add(.5, Vadj_x);
         Bnode.Add(1.5, Vnode_x);
         geom.Orthogonal(Bnode);

         Badj = face_x;
         Badj *= -2.;
         Badj.Add(.5, Vnode_x);
         Badj.Add(1.5, Vadj_x);
         geom.Orthogonal(Badj);

         // Compute dot products
         double badjn = Badj * n_vec;
         double badjnr = Badj * n_vec_R;
         double bnoden = Bnode * n_vec;
         double bnodenr = Bnode * n_vec_R;

         // Compute aFi, aFj
         double aFi = - 0.5 * dt * Vadj_nR_comp;
         aFi += bnoden;
         double aFj = 0.5*dt*Vnode_prev_it_nR_comp;
         aFj -= badjn;

         // Calculate D and c1 (perturbations only need to be handled explicitly)
         subtract(Vadj_n, Vnode_n, temp_vec); // V1 - V2 = temp_vec
         subtract(Vnode_half, Vadj_half, temp_vec_2); // A2-A1

         geom.Orthogonal(temp_vec_2);
         D = dt * (temp_vec * n_vec_R) + 2. * (n_vec * temp_vec_2);

         // Compute c1 (A.4a)
         subtract(Vnode_n, Vadj_n, temp_vec); // only change temp_vec, since temp_vec_2 is same from D calculation (half step representation)
         Dc1 = 2 * dt * (temp_vec * n_vec);
         
         // bF
         double bF = (3. - D / F) * bmn;
         bF += badjnr * Vadj_nR_comp;
         bF -= bnodenr * Vnode_prev_it_nR_comp;

         // Add in Dc1V3nper contribution (explicit)
         // V3nperp is set to 0 in option 1.
         if (mv_it_option == 2 || mv_it_option == 3)
         {

            /* Compute V3nperp using previous iteration */
            // add(0.5, Vnode_n, 0.5, Vadj_n, temp_vec);
            // V3nperp = -1. * (temp_vec * n_vec_R);
            // Compute V3nperp from flux
            V3nperp = -1. * (Vf * n_vec_R);

            /* Add in Dc1 explicitly at k to bF */
            bF += Dc1*V3nperp;
         }

         // Add contribution to the averaging objects
         double Mcoeff = 2. * pow(aFi, 2);
         double Rcoeff = aFj*Vadj_n_comp - bF;
         Rcoeff *= -2. * aFi;
         tensor(n_vec, n_vec, dm_tmp);

         Mi.Add(Mcoeff, dm_tmp);
         Ri.Add(Rcoeff, n_vec);

         // In the case of boundary nodes, check if this is an interior face
         if (bdr_ind != -1 && FI.IsInterior())
         {
            // We must add the ghost node contribution 
            // We will be employing reflective boundary conditions

            // Compute vector from node to face 
            Vector face_dx(dim);
            subtract(Vnode_x, face_x, face_dx);

            // Compute vector from node to adjacent node
            Vector Vadj_dx(dim);
            subtract(Vnode_x, Vadj_x, Vadj_dx);

            // Sod
            switch (bdr_ind)
            {
            case 1: // bottom
            case 3: // top
               n_vec[0] *= -1.; // preserves orientation of face
               face_dx[0] *= -1.;
               Vadj_dx[0] *= -1.;
               Vadj_n[1] *= -1.;
               Vface_n[1] *= -1.;
               break;
            
            case 2: // right
            case 4: // left
               n_vec[1] *= -1.; // preserves orientation of face
               face_dx[1] *= -1.;
               Vadj_dx[1] *= -1.;
               Vadj_n[0] *= -1.;
               Vface_n[0] *= -1.;
               break;
            
            case -1:
               // Not a boundary vertex
               continue;

            default:
               MFEM_ABORT("Incorrect bdr attribute encountered while enforcing mesh velocity BCs.\n");
               break;
            }

            // Recompute bmn for face
            double bmn = Vf * n_vec;
            bmn *= F;

            n_vec_R = n_vec;
            geom.Orthogonal(n_vec_R);
            Vnode_prev_it_nR_comp = Vnode_n * n_vec_R;

            // Get normal and rotated components of Vadj_n
            Vadj_n_comp = Vadj_n * n_vec;
            Vadj_nR_comp = Vadj_n * n_vec_R;
            
            /* Get flipped face_x and vadj_x */
            add(Vnode_x, face_dx, face_x);
            add(Vnode_x, Vadj_dx, Vadj_x);
            // and half location
            add(Vadj_x, dt/2., Vadj_n, Vadj_half);

            Vector ttvec(dim);
            subtract(Vnode_x, Vadj_x, ttvec);
            geom.Orthogonal(ttvec);
            double _val = ttvec * n_vec;
            if (_val < 1E-6)
            {
               cout << "boundary analysis:\n";
               cout << "node: " << node << endl;
               cout << "nodex: ";
               Vnode_x.Print();
               cout << "Vadj_x: ";
               Vadj_x.Print(cout);
               cout << "nvec: ";
               n_vec.Print(cout);
               assert(false);
            }

            // Compute vector constants for ghost node
            Bnode = face_x;
            Bnode *= -2.;
            Bnode.Add(.5, Vadj_x);
            Bnode.Add(1.5, Vnode_x);
            geom.Orthogonal(Bnode);

            Badj = face_x;
            Badj *= -2.;
            Badj.Add(.5, Vnode_x);
            Badj.Add(1.5, Vadj_x);
            geom.Orthogonal(Badj);

            // Compute dot products
            double badjn = Badj * n_vec;
            double badjnr = Badj * n_vec_R;
            double bnoden = Bnode * n_vec;
            double bnodenr = Bnode * n_vec_R;

            // Compute aFi, aFj
            double aFi = - 0.5 * dt * Vadj_nR_comp;
            aFi += bnoden;
            double aFj = 0.5*dt*Vnode_prev_it_nR_comp;
            aFj -= badjn;

            // Calculate D and c1 (perturbations only need to be handled explicitly)
            subtract(Vadj_n, Vnode_n, temp_vec); // V1 - V2 = temp_vec
            subtract(Vnode_half, Vadj_half, temp_vec_2); // A2-A1

            geom.Orthogonal(temp_vec_2);
            D = dt * (temp_vec * n_vec_R) + 2. * (n_vec * temp_vec_2);

            // Compute c1 (A.4a)
            subtract(Vnode_n, Vadj_n, temp_vec); // only change temp_vec, since temp_vec_2 is same from D calculation (half step representation)
            Dc1 = 2 * dt * (temp_vec * n_vec);
            
            // bF
            double bF = (3. - D / F) * bmn;
            bF += badjnr * Vadj_nR_comp;
            bF -= bnodenr * Vnode_prev_it_nR_comp;

            // Add in Dc1V3nper contribution (explicit)
            // V3nperp is set to 0 in option 1.
            if (mv_it_option == 2 || mv_it_option == 3)
            {

               /* Compute V3nperp using previous iteration */
               // add(0.5, Vnode_n, 0.5, Vadj_n, temp_vec);
               // V3nperp = -1. * (temp_vec * n_vec_R);
               // Compute V3nperp from flux
               V3nperp = -1. * (Vf * n_vec_R);

               /* Add in Dc1 explicitly at k to bF */
               bF += Dc1*V3nperp;
            }

            // Add contribution to the averaging objects
            double Mcoeff = 2. * pow(aFi, 2);
            double Rcoeff = aFj*Vadj_n_comp - bF;
            Rcoeff *= -2. * aFi;
            tensor(n_vec, n_vec, dm_tmp);

            Mi.Add(Mcoeff, dm_tmp);
            Ri.Add(Rcoeff, n_vec);
         } // End ghost node
      }

      // Average node_v
      Mi.Invert();
      Mi.Mult(Ri, predicted_node_v);

      // Theta average with previous velocity
      if (do_theta_averaging)
      {
         predicted_node_v *= theta; 
         predicted_node_v.Add(1. - theta, Vnode_n);
      }

      // if (node == 99)
      // {
      //    cout << "node: " << node << endl;
      //    cout << "Mi inverse: ";
      //    Mi.Print(cout);
      //    cout << "Ri: ";
      //    Ri.Print(cout);
      //    cout << "predicted v: ";
      //    predicted_node_v.Print(cout);
      // }

      // Put velocity in S
      // Gauss-Seidel > Jacobi
      // geom.UpdateNodeVelocity(S, node, predicted_node_v);
      // Try Jacobi
      geom.UpdateNodeVelocity(S_new, node, predicted_node_v);
   } // End node iterator
   S = S_new;
}


/****************************************************************************************************
* Function: ComputeCellAverageVelocityAtNode
* Parameters:
*  S           - BlockVector representing FiniteElement information tnp1
*  S_old       - BlockVector representing FiniteElement information tn
*  node        - Node on which to compute the weighted average cell velocities
*  is_weighted - boolean representing if the average should be a weighted 
*                average or not.
*  td_flag     - Which time discretized value should be used for Vcell:
*                    0 - tn
*                    1 - tnp1
*                    2 - tnp1/2
*  mv_gf       - Corresponding mesh velocities
*  node_v      - Averaged velocity on node
*
* Purpose:
*  To compute the weighted average of the adjacent cell velocities.
****************************************************************************************************/
void LagrangianLOOperator::ComputeCellAverageVelocityAtNode(const Vector &S, const Vector &S_old, const int node, const bool &is_weighted, int td_flag, Vector &node_v)
{
   // cout << "========================================\n"
   //      << "ComputeCellAverageVelocityAtNode\n"
   //      << "========================================\n";
   Array<int> cells_row;
   Vector Uc(dim+2), Ucnp1(dim+2), Vcell(dim), Vcellnp1(dim);
   double sum_k = 0.;
   node_v = 0.;

   /* Get adjacent cells */
   vertex_element->GetRow(node, cells_row);
   int cells_length = cells_row.Size();

   /* Iterate over adjacent cells */
   for (int cell_it = 0; cell_it < cells_length; cell_it++)
   {
      /* Get cell index and volume */
      int el_index = cells_row[cell_it];

      /* Retrieve cell velocity */
      switch (td_flag)
      {
      case 0: // S_OLD
      {
         GetCellStateVector(S_old, el_index, Uc);
         pb->velocity(Uc, Vcell);
         break;
      }
      case 1: // S_NEW 
      {
         GetCellStateVector(S, el_index, Uc);
         pb->velocity(Uc, Vcell);
         break;
      }
      case 2: // MIDPOINT
      {
         GetCellStateVector(S_old, el_index, Uc);
         GetCellStateVector(S, el_index, Ucnp1);

         pb->velocity(Uc, Vcell);
         pb->velocity(Ucnp1, Vcellnp1);

         Vcell.Add(1., Vcellnp1);
         Vcell *= .5;
         break;
      }
      default:
         MFEM_ABORT("Should not have reached this point.\n");
         break;
      }

      /* Add in contribution of cell */
      if (is_weighted) // weighted
      {
         double Kcn = pmesh->GetElementVolume(el_index);
         sum_k += Kcn;
         node_v.Add(Kcn, Vcell);
      }
      else // arithmetic
      {
         sum_k += 1.;
         node_v.Add(1., Vcell);
      }
   }

   /* Finally, take the average */
   node_v /= sum_k;
}


/****************************************************************************************************
* Function:  CalcMassVolumeVector
* Parameters:
*  S        - BlockVector representing FiniteElement information at tn+1
*  S_old    - BlockVector representing state at time tn
*  dt       - Current time step
*  massvec  - Vector of length num_cells representing mass of each cell 
*             at time tn multiplied by the specific volume of the cell
*             at time tn+1.
*
* Purpose:
*  This function computes the massvec that is used as the constraint in the 
*  HiOp implementation of the Lagrange Multiplier method to solve for the 
*  mesh velocity that guarantees mass conservation cell-wise.
*
*  Note: This function assumes that m_hpv has been populated and contains
*  the initial mass of all cells in the mesh.
****************************************************************************************************/
void LagrangianLOOperator::CalcMassVolumeVector(const Vector &S, const double &dt, Vector &massvec)
{
   // cout << "=======================================\n"
   //      << "          CalcMassVolumeVector         \n"
   //      << "=======================================\n";

   Vector *sptr = const_cast<Vector*>(&S);
   ParGridFunction sv_gf;
   sv_gf.MakeRef(&L2, *sptr, block_offsets[1]);

   /* Iterate over all cells in mesh and compute constraint */
   for (int cell_it = 0; cell_it < NDofs_L2; cell_it++)
   {
      double Tnp1 = sv_gf.Elem(cell_it);
      double m = m_hpv->Elem(cell_it);

      /* Compute constraint on cell */
      double val = Tnp1 * m;

      massvec[cell_it] = val;
      if (val <= 0.)
      {
         MFEM_ABORT("massvec should be positive.\n");
      }
   }
}


/****************************************************************************************************
* Function:  CalcCellAveragedVelocityVector
* Parameters:
*  S           - BlockVector representing FiniteElement information tnp1
*  S_old       - BlockVector representing FiniteElement information tn
*  is_weighted - boolean representing if the average should be a weighted 
*                average or not.
*  td_flag     - Which time discretized value should be used for Vcell:
*                    0 - tn
*                    1 - tnp1
*                    2 - tnp1/2
*  mv_gf       - Corresponding mesh velocities
*  mv_gf       - Corresponding mesh velocities
*
* Purpose:
*  
****************************************************************************************************/
void LagrangianLOOperator::CalcCellAveragedCornerVelocityVector(const Vector &S, const Vector &S_old, const bool &is_weighted, int td_flag, ParGridFunction &mv_gf_l)
{
   // cout << "=======================================\n"
   //      << "    CalcCellAveragedVelocityVector     \n"
   //      << "=======================================\n";
   Vector node_v(dim);

   /* Iterate over corner node */
   for (int node = 0; node < NVDofs_H1; node++)
   {
      ComputeCellAverageVelocityAtNode(S, S_old, node, is_weighted, td_flag, node_v);

      /* Put node_v in place */
      geom.UpdateNodeVelocityVecL(mv_gf_l, node, node_v);
   }
}


/****************************************************************************************************
* Function:  DistributeFaceViscosityToVelocity
* Parameters:
*  S        - BlockVector representing FiniteElement information at tn+1
*  mv_gf    - Gridfunction representing velocity of all geometric corner nodes
*
* Purpose:
* Note: This function assumes that mv_gf is defined on all geometric corner nodes of the mesh.  
****************************************************************************************************/
void LagrangianLOOperator::DistributeFaceViscosityToVelocity(const Vector &S, Vector &mv_gf)
{
   MFEM_ABORT("Function deprecated\n");
   assert(mv_gf.Size() == dim * NVDofs_H1);

   mfem::Mesh::FaceInformation FI;
   int c, cp, node_i, node_j;
   double F, d, coeff;
   Vector n_int(dim), Vi(dim), Vj(dim), Uc(dim+2), Ucp(dim+2);
   Array<int> face_dofs_row;

   for (int face = 0; face < num_faces; face++) // face iterator
   {
      FI = pmesh->GetFaceInformation(face);
      if (FI.IsInterior())
      {
         c = FI.element[0].index; cp = FI.element[1].index;
         CalcOutwardNormalInt(S, c, face, n_int);
         F = n_int.Norml2();
         n_int /= F;
         d = dij_sparse->Elem(c,cp);
         GetCellStateVector(S, c, Uc);
         GetCellStateVector(S, cp, Ucp);

         /* Get corner indices */
         H1.GetFaceDofs(face, face_dofs_row);
         node_i = face_dofs_row[0];
         node_j = face_dofs_row[1];
         geom.GetNodeVelocityVecL(mv_gf, node_i, Vi);
         geom.GetNodeVelocityVecL(mv_gf, node_j, Vj);
         coeff = 0.5 * d * (Ucp[0] - Uc[0]) / F;

         /* Handle ghost face on boundary */
         if (BdrVertexIndexingArray[node_i] == -1) // interior
         {
            Vi.Add(coeff, n_int);
         } else {
            Vi.Add(2*coeff, n_int);
         }
         if (BdrVertexIndexingArray[node_j] == -1) // interior
         {
            Vj.Add(coeff, n_int);
         } else {
            Vj.Add(2*coeff, n_int);
         }
         
         /* Updated nodal velocity */
         geom.UpdateNodeVelocityVecL(mv_gf, node_i, Vi);
         geom.UpdateNodeVelocityVecL(mv_gf, node_j, Vj);
      }
   }
   // cout << "DistributeFaceViscosityToVelocity - DONE\n";
}


/************
 * Function: SolveHiopDense
 * TODO: combine with SolveHiop function with a dense/sparse flag
 */
void LagrangianLOOperator::SolveHiOpDense(const Vector &S, const Vector &S_old, const int & target_option, const double &t, const double &dt, ParGridFunction &mv_gf_l)
{
   // cout << "=======================================\n"
   //      << "            SolveHiOpDense             \n"
   //      << "=======================================\n";
   assert(target_option < 10);
   assert(mv_gf_l.Size() == dim*NVDofs_H1);

   MFEM_ABORT("Why are you using the dense implementation when Sparse is available?\n");

   Vector massvec(NDofs_L2);
   ParGridFunction V_target(&H1_L), mv_gf_l_out(mv_gf_l);
   Vector* sptr = const_cast<Vector*>(&S);
   x_gf.MakeRef(&H1, *sptr, block_offsets[0]);

   OptimizationSolver *optsolver = NULL;
   const int optimizer_type = 2; // TODO: change this to be a set param
   if (optimizer_type == 2)
   {
#ifdef MFEM_USE_HIOP
   HiopNlpOptimizer *tmp_opt_ptr = new HiopNlpOptimizer();
   optsolver = tmp_opt_ptr;
#else
      MFEM_ABORT("MFEM is not built with HiOp support!");
#endif
   }

   /* Set min/max velocities */
   Vector xmin(V_target.Size()), xmax(V_target.Size());
   xmin = -1.E12;
   xmax = 1.E12;

   /* Compute equality restrictions */
   CalcMassVolumeVector(S, dt, massvec);

   /* Options */
   bool is_weighted = false;

   /* Compute targeted velocity */
   switch (target_option)
   {
   case 0: // arithmetic average of adjacent cells with distributed viscosity
   {
      int td_flag = 0;
      CalcCellAveragedCornerVelocityVector(S, S_old, is_weighted, td_flag, V_target);
      break;
   }
   case 1: // arithmetic average of adjacent cells with viscosity
   {
      int td_flag = 0;
      CalcCellAveragedCornerVelocityVector(S, S_old, is_weighted, td_flag, V_target);
      DistributeFaceViscosityToVelocity(S_old, V_target);
      break;
   }
   case 2:
   {
      ComputeGeoVNormal(S, V_target);
      break;
   }
   case 3: // weighted average of adjacent cell midpoint in time
   {
      /* Petrov-Galerkin Justified weighted average */
      bool is_weighted = true;
      int td_flag = 2;
      CalcCellAveragedCornerVelocityVector(S, S_old, is_weighted, td_flag, V_target);
      break;
   }
    case 4: // weighted average of adjacent cell midpoint in time with distributed viscosity
   {
      /* Petrov-Galerkin Justified weighted average */
      bool is_weighted = true;
      int td_flag = 2;
      CalcCellAveragedCornerVelocityVector(S, S_old, is_weighted, td_flag, V_target);
      DistributeFaceViscosityToVelocity(S_old, V_target);
      break;
   }
   default:
      MFEM_ABORT("Invalid mesh velocity target in Hiop solver.\n");
   }

   /* Linearize target velocity */
   if (this->do_mv_linearization)
   {
      ParGridFunction V_target_lin(V_target);
      ComputeLinearizedNodeVelocities(V_target, V_target_lin, t, dt);
      V_target = V_target_lin;
   }

   /* Solve for the corner node velocities */
   OptimizedMeshVelocityProblemDense omv_problem(dim, geom, V_target, massvec, x_gf, NDofs_L2, dt, xmin, xmax);
   optsolver->SetOptimizationProblem(omv_problem);
   optsolver->SetMaxIter(this->corner_velocity_MC_num_iterations);
   optsolver->SetPrintLevel(0);
   optsolver->SetRelTol(1E-6);
   optsolver->SetAbsTol(1E-8);
   mv_gf_l = 0.;
   optsolver->Mult(mv_gf_l, mv_gf_l_out);

   mv_gf_l = mv_gf_l_out;

   // delete optsolver;
   // optsolver = nullptr;
}

/****************************************************************************************************
* Function:  ComputeVelocityLumpedMass
* Parameters:
*  row_sums - Vector representing the diagonal from the lumped mass matrix, to be changed
*             in this function
*
* Purpose:
*  This function computes the lumped mass matrix for the velocity space.
****************************************************************************************************/
void LagrangianLOOperator::ComputeVelocityLumpedMass(Vector & row_sums)
{
   // cout << "=======================================\n"
   //      << "       ComputeVelocityLumpedMass       \n"
   //      << "=======================================\n";

   /* Compute the Mis for the velocity discretization */
   row_sums.SetSize(H1Lc.GetNDofs());
   BilinearForm pbl(&H1Lc);
   Coefficient *sigma = new ConstantCoefficient(1.0);
   pbl.AddDomainIntegrator(new MassIntegrator());
   pbl.Assemble();
   SparseMatrix mass;
   pbl.FormSystemMatrix(mfem::Array<int>(), mass);
   // Vector row_sums(81);
   mass.GetRowSums(row_sums);

   // cout << "row_sums: " << endl;
   // row_sums.Print(cout);
   // cout << "mass matrix: \n";
   // mass.Print(cout);
   // assert(false);
   // delete sigma;
   // assert(false);
}


/****************************************************************************************************
* Function:  ComputeVelocityLumpedMass
* Parameters:
*  row_sums - Vector representing the diagonal from the lumped mass matrix, to be changed
*             in this function
*
* Purpose:
*  This function computes the lumped mass matrix for the velocity space.
****************************************************************************************************/
void LagrangianLOOperator::ComputeVelocityLumpedMassByHand(const Vector &S, Vector & row_sums)
{
   // cout << "=======================================\n"
   //      << "     ComputeVelocityLumpedMassHand     \n"
   //      << "=======================================\n";
   row_sums.SetSize(H1Lc.GetNDofs());
   Array<int> cells_row;

   /* Compute the Mis for the velocity discretization */
   for (int i = 0; i < H1Lc.GetNDofs(); i++)
   {
      // cout << "--- i: " << i << endl;
      vertex_element->GetRow(i, cells_row);
      int cells_length = cells_row.Size();
      double sum = 0.;
      for (int cell_it = 0; cell_it < cells_length; cell_it++)
      {
         /* Get cell index and volume */
         int el_index = cells_row[cell_it];
         double Kcn = pmesh->GetElementVolume(el_index);
         // cout << "\tel: " << el_index << ", vol: " << Kcn << endl;
         sum += Kcn;

      }
      sum *= 1./3.;
      row_sums[i] = sum;
   }
}


/****************************************************************************************************
* Function:  ComputeAverageMassAtGeo
* Parameters:
*  S           - BlockVector representing FiniteElement information at tn+1
*  vec_weights - Vector representing average mass of the adjacent cells to a given node
*
* Purpose:
*  Compute the average mass of the adjacent cells to a given node.
****************************************************************************************************/
void LagrangianLOOperator::ComputeAverageMassAtGeo(const Vector &S, Vector &vec_weights)
{
   vec_weights.SetSize(H1Lc.GetNDofs());
   Array<int> cells_row;

   /* Iterate over geometric corner nodes */
   for (int i = 0; i < H1Lc.GetNDofs(); i++)
   {
      /* Get adjacent cells*/
      vertex_element->GetRow(i, cells_row);
      int cells_length = cells_row.Size();

      /* Compute average mass of adjacent cells*/
      double avg = 0.;
      for (int cell_it = 0; cell_it < cells_length; cell_it++)
      {
         int el_index = cells_row[cell_it];
         const double m = m_hpv->Elem(cell_it);
         avg += m;
      }
      avg /= cells_length;

      /* Store that value */
      vec_weights[i] = avg;
   }
}


/****************************************************************************************************
* Function:  SolveHiOp
* Parameters:
*  S        - BlockVector representing FiniteElement information at tn+1
*  S_old    - BlockVector representing state at time tn
*  dt       - Current time step
*  massvec  - Vector of length num_cells representing mass of each cell 
*             at time tn multiplied by the specific volume of the cell
*             at time tn+1.
*
* Purpose:
*  This function solves for the velocity vector at the corner vertices by using
*  an OptimizationSolver from MFEM.
****************************************************************************************************/
void LagrangianLOOperator::SolveHiOp(const Vector &S, const Vector &S_old, const int &lm_option, const int & target_option, const double &t, const double &dt, ParGridFunction &mv_gf_l)
{
   // cout << "=======================================\n"
   //      << "               SolveHiOp               \n"
   //      << "=======================================\n";
   assert(target_option < 10);
   assert(mv_gf_l.Size() == Vsize_H1L);

   Vector massvec(NDofs_L2);
   ParGridFunction V_target(&H1_L), mv_gf_l_out(mv_gf_l);
   Vector* sptr = const_cast<Vector*>(&S);
   x_gf.MakeRef(&H1, *sptr, block_offsets[0]);

   /* Set min/max velocities */
   Vector xmin(Vsize_H1L), xmax(Vsize_H1L);
   xmin = -1.E12;
   xmax = 1.E12;

   /* Adjust the above for boundary conditions */
   // if (pb->has_mv_boundary_conditions())
   // {
   //    double bdr_tol = 1.E-12;
   //    for (int i = 0; i < ess_tdofs.Size(); i++) { 
   //       xmin(ess_tdofs[i]) = bdr_vals[i] - bdr_tol; 
   //       xmax(ess_tdofs[i]) = bdr_vals[i] + bdr_tol;
   //    }
   // }
   // cout << "xmin: ";
   // xmin.Print(cout);
   // cout << "xmax: ";
   // xmax.Print(cout);

   /* Compute equality restrictions */
   CalcMassVolumeVector(S, dt, massvec);

   /* Options */
   bool is_weighted = false;

   /* 
   Calculate sparsity patterns for Gradient of the constraint vector as it is the same
   whether or not the target or viscous option is used
   */
   SetHiopConstraintGradSparsityPattern(pmesh, num_elements, NVDofs_H1, HiopCGradIArr, HiopCGradJArr);

   OptimizationProblem * omv_problem = NULL;
   switch (lm_option)
   {
      case 1: // Have a target velocity
      {
         /* Choose whether to add viscosity to the objective function */
         bool use_obj_visc = false;
         if (this->mv_target_visc_coeff > 0.)
         {
            use_obj_visc = true;
         }

         /* Calculate sparsity pattern Hessian which depends on if we choose to use viscosity */
         if (use_obj_visc)
         {
            SetHiopHessianSparsityPatternViscous(pmesh, geom, H1, NVDofs_H1, HiopHessIArr, HiopHessJArr);
         }
         else 
         {
            SetHiopHessianSparsityPattern(pmesh, H1, NVDofs_H1, HiopHessIArr, HiopHessJArr);
         }

         // Sparsity pattern for Boundary Constraint implementation which did not yield good results
         // SetHiopBoundaryConstraintGradSparsityPattern(ess_tdofs, HiopDGradIArr, HiopDGradJArr, HiopDGradData);

         /* Calculate target velocity */
         switch (target_option)
         {
         case 0: // arithmetic average of adjacent cells with distributed viscosity
         {
            int td_flag = 0;
            CalcCellAveragedCornerVelocityVector(S, S_old, is_weighted, td_flag, V_target);
            break;
         }
         case 1: // arithmetic average of adjacent cells with viscosity
         {
            int td_flag = 0;
            CalcCellAveragedCornerVelocityVector(S, S_old, is_weighted, td_flag, V_target);
            DistributeFaceViscosityToVelocity(S_old, V_target);
            break;
         }
         case 2:
         {
            ComputeGeoVNormal(S, V_target);
            break;
         }
         case 3: // weighted average of adjacent cell midpoint in time
         {
            /* Petrov-Galerkin Justified weighted average */
            is_weighted = true;
            int td_flag = 2;
            CalcCellAveragedCornerVelocityVector(S, S_old, is_weighted, td_flag, V_target);
            break;
         }
         case 4: // weighted average of adjacent cell midpoint in time with distributed viscosity
         {
            /* Petrov-Galerkin Justified weighted average */
            is_weighted = true;
            int td_flag = 2;
            CalcCellAveragedCornerVelocityVector(S, S_old, is_weighted, td_flag, V_target);
            DistributeFaceViscosityToVelocity(S_old, V_target);
            break;
         }

         default:
            MFEM_ABORT("Invalid mesh velocity target in Hiop solver.\n");
         }

         if (this->do_mv_linearization)
         {
            ParGridFunction V_target_lin(V_target);
            ComputeLinearizedNodeVelocities(V_target, V_target_lin, t, dt);
            V_target = V_target_lin;
         }

         /* Enforce BCs on target velocity */
         if (pb->has_mv_boundary_conditions())
         {
            for (int i = 0; i < ess_tdofs.Size(); i++) { V_target(ess_tdofs[i]) = bdr_vals[i]; }
         }

         /* Compute lumped mass in velocity space */
         Vector _vecWeights;
         // ComputeVelocityLumpedMass(_vecWeights);
         // ComputeVelocityLumpedMassByHand(S, _vecWeights);
         ComputeAverageMassAtGeo(S, _vecWeights);

         /* Instantiate problem */
         omv_problem = new TargetOptimizedMeshVelocityProblem(
            dim, geom, V_target, massvec, x_gf, NDofs_L2, dt, xmin, xmax, 
            HiopHessIArr, HiopHessJArr, HiopCGradIArr, HiopCGradJArr, 
            HiopDGradIArr, HiopDGradJArr, HiopDGradData, bdr_vals,
            ess_tdofs, BdrVertexIndexingArray, this->mv_target_visc_coeff,
            _vecWeights);
         break;
      }
      case 2: // Viscous objective function
      {
         /* Need to modify the solution bound restrictions to accomodate boundary conditions */
         if (pb->has_mv_boundary_conditions())
         {
            double bdr_tol = 1.E-12;
            for (int i = 0; i < ess_tdofs.Size(); i++) {
               xmin(ess_tdofs[i]) = bdr_vals[i] - bdr_tol;
               xmax(ess_tdofs[i]) = bdr_vals[i] + bdr_tol;
            }
         }

         /* Hessian sparsity pattern is slightly different due to the different objective function */
         SetHiopHessianSparsityPatternViscous(pmesh, geom, H1, NVDofs_H1, HiopHessIArr, HiopHessJArr);

         /* Instantiate problem */
         omv_problem = new ViscousOptimizedMeshVelocityProblem(dim, geom, massvec, x_gf, NDofs_L2, dt, xmin, xmax, HiopHessIArr, HiopHessJArr, HiopCGradIArr, HiopCGradJArr, ess_tdofs, BdrVertexIndexingArray);
         break;
      }
      default:
      {
         MFEM_ABORT("Invalid Lagrange Multiplier option.\n");
         break;
      }
   }

   chrono_hiop.Start();
   OptimizationSolver *optsolver = NULL;
#ifdef MFEM_USE_HIOP
   HiopNlpSparseOptimizer *tmp_opt_ptr = new HiopNlpSparseOptimizer();
   tmp_opt_ptr->SetNNZSparse(8); // FIXME: adjust hardcoded parameter
   optsolver = tmp_opt_ptr;
#else
      MFEM_ABORT("MFEM is not built with HiOp support!");
#endif

   /* Solve for corner node velocities */
   // optsolver->SetOptimizationProblem(*omv_problem);
   optsolver->SetOptimizationProblem(*omv_problem);
   optsolver->SetMaxIter(this->corner_velocity_MC_num_iterations);
   optsolver->SetPrintLevel(0);
   optsolver->SetRelTol(1E-1);
   optsolver->SetAbsTol(1E-1);
   mv_gf_l = 0.;
   // ComputeGeoVNormal(S, mv_gf_l);
   // is_weighted = true;
   // int td_flag = 2;
   // CalcCellAveragedCornerVelocityVector(S, S_old, is_weighted, td_flag, mv_gf_l);
   // DistributeFaceViscosityToVelocity(S_old, mv_gf_l);
   optsolver->Mult(mv_gf_l, mv_gf_l_out); 
   chrono_hiop.Stop();

   mv_gf_l = mv_gf_l_out;

   /* Properly dispose of allocated memory */
   #ifdef MFEM_USE_HIOP
      delete tmp_opt_ptr;
      tmp_opt_ptr = nullptr;
   #endif
   
   delete omv_problem;
   optsolver = nullptr;
   omv_problem = nullptr;
}


/**********************************************************
 * Linearization functions where v_geo_gf is not used but
 * is passed in as an argument.
 * Also, the BlockVector S is not updated in these functions,
 * but the linearized velocity is returned by way of a reference
 * argument.  This is necessary to utilize a target velocity 
 * function in the Lagrange Multiplier application
 **********************************************************/
/****************************************************************************************************
* Function: ComputeCiGeo
* Parameters:
*  node  - Global DOF corresponding to the node the CiGeo should be evaluated at.
*          DOF should correspond to a corner or face node, as there is no reason 
*          evaluate CiGeo at a cell center.
*  res   - DenseMatrix representing CiGeo evaluated at the provided node.
*
* Purpose:
*  This function serves to evaluate CiGeo at a given global dof.  Recall that
*  CiGeo is defined as
*     CiGeo =  [1 / Sum |K|] Sum int grad v_geo_gf
*
*  This function is essential to go from ViGeo to Vi.
****************************************************************************************************/
void LagrangianLOOperator::ComputeCiGeo(const ParGridFunction &mv_gf_l, const int & node, DenseMatrix & res)
{
   // cout << "=======================================\n"
   //      << "             ComputeCiGeo              \n"
   //      << "=======================================\n";
   assert(node < NDofs_H1); // "Invalid nodal index"
   res.SetSize(dim);
   res = 0.;

   // Check if node corresponds to a face or to a 
   // vertex and get adjacent cell indices
   double denom = 0.;
   DenseMatrix dm_temp(dim);
   mfem::Array<int> row;

   DofEntity entity;
   int EDof;
   GetEntityDof(node, entity, EDof);

   switch (entity)
   {
      case 0: // corner
      {
         vertex_element->GetRow(node, row);
         break;
      }
      case 1: // face
      {
         face_element->GetRow(EDof, row);
         break;
      }
      case 2: // cell
      {
         MFEM_ABORT("No need to compute C_i at cell centers.\n");
      }
      default:
      {
         MFEM_ABORT("Invalid entity\n");
      }
   }

   int row_length = row.Size();
   for (int row_it = 0; row_it < row_length; row_it++)
   {
      int row_el = row[row_it];
      double cell_vol = pmesh->GetElementVolume(row_el);
      denom += cell_vol;

      dm_temp = 0.;
      IntGrad(mv_gf_l, row_el, dm_temp);
      res.Add(cell_vol, dm_temp);
   }
   res *= 1./denom;
}


/****************************************************************************************************
* Function: IntGrad
* Parameters:
*    cell - index corrseponding to the cell (K_c)
*    res  - (dim x dim) DenseMatrix representing the outer product of Vfm with 
*           the gradient of its corresponding scalar RT shape function.
* Purpose:
*    This function calculates the integral of the gradient of the v_geo velocity function
*    on a given cell.

* NOTE: This function assumes that dim > 1.
* NOTE: This function assumes that the functions:
*          - LagrangianLOOperator::ComputeIntermediateFaceVelocities()
*          - LagrangianLOOperator::ComputeGeoVRaviart() 
*       have already been called.  If the functions have not been called, then the 
*       returned velocity will be 0.
****************************************************************************************************/
void LagrangianLOOperator::IntGrad(const ParGridFunction &mv_gf_l, const int cell, DenseMatrix & res)
{
   // cout << "=======================================\n"
   //      << "             IntGrad            \n"
   //      << "=======================================\n";

   ParGridFunction H1Lc_gf(&H1Lc), mv_gf_copy(mv_gf_l);
   const int _size = H1Lc.GetVSize();

   ElementTransformation * trans = pmesh->GetElementTransformation(cell);
   Vector grad(dim), row(dim);
   res.SetSize(dim);
   res = 0., row = 0.;

   // Iterate over quadrature
   for (int i = 0; i < RT_ir.GetNPoints(); i++)
   {
      const IntegrationPoint &ip = RT_ir.IntPoint(i);
      trans->SetIntPoint(&ip);

      for (int j = 0; j < dim; j++)
      {
         H1Lc_gf.MakeRef(&H1Lc, mv_gf_copy, j*_size);
         H1Lc_gf.GetGradient(*trans, grad);

         // Put information into Dense Matrix
         res.GetRow(j, row);
         row.Add(ip.weight, grad);
         res.SetRow(j, row);
      }
   }
}


/****************************************************************************************************
* Function: ComputeNodeVelocitiesFromVgeo
* Parameters:
*  S        - BlockVector representing FiniteElement information
*  t        - Current time
*  dt       - Current time step
*  flag     - Flag used to indicate testing, can be "testing" or "NA"
*  test_vel - Velocity used for testing
*
* Purpose:
*  This function iterates over all Serendipity nodes and computes the 
*  velocity with which to move the mesh from the geometric velocity.
*  
****************************************************************************************************/
void LagrangianLOOperator::ComputeLinearizedNodeVelocities(
   const ParGridFunction &mv_gf_l, ParGridFunction &mv_gf_linearized,
   const double &t, const double &dt, const string, 
   void (*test_vel)(const Vector&, const double&, Vector&))
{
   // cout << "ComputeLinearizedNodeVelocities\n";
   chrono_mm_lin.Start();
   assert(mv_gf_l.Size() == dim * NVDofs_H1);
   Vector node_v(dim);
   bool is_dt_changed = false;

   // Iterate over vertices and faces
   for (int node = 0; node < NVDofs_H1; node++) // Vertex and face iterator
   {
      ComputeLinearizedNodeVelocity(mv_gf_l, node, dt, node_v, is_dt_changed);

      // if (node_v[0] != node_v[0] || node_v[1] != node_v[1])
      // {
      //    cout << "NaN velocity encountered in ComputeNodeVelocities at node: " << node << endl;
      //    cout << "node_v in violation: ";
      //    node_v.Print(cout);
      //    MFEM_ABORT("Aborting due to NaNs.\n");
      // }

      geom.UpdateNodeVelocityVecL(mv_gf_linearized, node, node_v);

      // If we restricted the timestep, we must recompute the vertex velocities that were computed previously
      // if (is_dt_changed)
      // {
      //    vertex = -1;
      //    is_dt_changed = false;
      //    cout << "Restarting vertex iterator\n";
      // }
   } // End Vertex iterator
   chrono_mm_lin.Stop();
}


/****************************************************************************************************
* Function: ComputeLinearizedNodeVelocity
* Parameters:
*  S             - BlockVector representing FiniteElement information
*  node          - Global DOF to compute velocity at
*  dt            - Current time step
*  node_v        - The computed velocity to be returned
*  is_dt_changed - Boolean representing if the timestep was modified
*                  due to a timestep restriction from the computation
*                  of alpha_i
*
* Purpose:
*  This function is responsible for computing the nodal velocity
*  from the geometric representation of the velocity through the
*  algorithm outlined in the paper:
*     [alpha_i I - (dt/2) C_i]^-1 Vigeo.
*  
****************************************************************************************************/
void LagrangianLOOperator::ComputeLinearizedNodeVelocity(
   const ParGridFunction &mv_gf_l, const int & node, const double & dt, Vector &node_v, bool &is_dt_changed)
{
   // cout << "=======================================\n"
   //      << "      ComputeNodeVelocityFromVgeo      \n"
   //      << "=======================================\n";
   switch (dim)
   {
      case 1:
      {
         // Remark 5.3 in the paper states that since each geometric node belongs to one
         // face only, the RT velocity is equal to the IFV previously computed.
         assert(node < num_faces);
         GetIntermediateFaceVelocity(node, node_v);
         break;
      }
      default:
      {
         DenseMatrix Ci(dim);
         Vector Vgeo(dim);
         double d = 0.;
         
         Ci = 0., Vgeo = 0., node_v = 0.;

         ComputeCiGeo(mv_gf_l, node, Ci);

         // Enforce time restriction imposed by calculation of alpha_i
         // if (dim == 2)
         // {
         //    double trace = Ci.Trace();
         //    double det = Ci.Det();
         //    double val = 2. * sqrt(det);
         //    cout << "trace: " << trace << ", det: " << det << endl;
         //
         //    // alpha_i must be positive restriction
         //    if (det <= 1.e-12 && trace < 0.)
         //    {
         //       if (dt > 2. / abs(trace))
         //       {
         //          dt = 2. / (abs(trace) + 1.);
         //          is_dt_changed = true;
         //          cout << "timestep restriction from velocity computation, det = 0 case.\n";
         //       }
         //    }
         //  
         //    // alpha_i must be real
         //    else if (det > 0.)
         //    {
         //       if (trace < val)
         //       {
         //          // Enforce timestep restriction
         //          double zero2 = 1. / (-trace + val);
         //          if (dt > 2 * zero2)
         //          {
         //             cout << "timestep restriction from velocity computation\n";
         //             cout << "dt: " << dt << ", z2: " << zero2 << endl;
         //             dt = 2. * zero2;
         //             is_dt_changed = true;
         //          }
         //       }
         //       else 
         //       {
         //          cout << "positive determinant but no timestep restriction needed.\n";
         //       }
         //    }
         //    else
         //    {
         //       cout << "negative determinant, no timestep restriction needed.\n";
         //    }
         // } // End time restriction from velocity computation

         geom.GetNodeVelocityVecL(mv_gf_l, node, Vgeo);
         
         ComputeDeterminant(Ci, dt, d, node);

         // Compute V_i^n
         DenseMatrix _mat(Ci);
         _mat *= - dt / 2.;
         for (int i = 0; i < dim; i++)
         {
            _mat(i,i) += d;
         }

         _mat.Invert();
         _mat.Mult(Vgeo, node_v);

         break;
      }
   }
}


/****************************************************************************************************
* Function: Raviart-Thomas Velocity Functions
****************************************************************************************************/

/****************************************************************************************************
* Function: ComputeGeoVRaviart
* Parameters:
*  S    - BlockVector corresponding to nth timestep that stores mesh information, 
*         mesh velocity, and state variables.
*
* Purpose:
*  This function computes the geometric velocity using the Raviart-Thomas elements,
*  which are face-based vector finite elements.  This function must take into account
*  the Piola transformations, which lucky for us MFEM will handle by use of the 
*  ElementTransformation class.
*
*  Note: Vigeo is only the geometric velocity, and is not the velocity used 
*        to move the mesh.  For the mesh velocity, one must solve for the 
*        alpha_i, which process is found in ComputeNodeVelocitiesFromVgeo().
*
*  Note: This function fills the Q1 velocity field v_geo_gf and also the mv_gf
*        in the BlockVector S.
****************************************************************************************************/
void LagrangianLOOperator::ComputeGeoVRaviart(Vector &S)
{
   // cout << "=======================================\n"
   //      << "          ComputeGeoVRaviart           \n"
   //      << "=======================================\n";

   Vector node_v(dim), n_vec(dim);
   Vector ivf(dim), v_cell(dim);
   Array<int> element_row, cell_faces_row;
   int row_length, cell_faces_length;

   DofEntity entity;
   int EDof;

   for (int node = 0; node < NDofs_H1 - NDofs_L2; node++) // Vertex iterator
   {
      node_v = 0.; // Reset node velocity

      // Check if node corresponds to a face or to a 
      // vertex and get adjacent cell indices and edge indices
      GetEntityDof(node, entity, EDof);

      switch (entity)
      {
         case 0: // corner
         {
            // cout << "ComputeGeoVRaviart::corner node\n";
            vertex_element->GetRow(node, element_row);
            row_length = element_row.Size();

            break;
         }
         case 1: // face
         {
            // cout << "ComputeGeoVRaviart::face node for face: " << EDof << endl;
            face_element->GetRow(EDof, element_row);
            row_length = element_row.Size();
            // cout << "row length: " << row_length << endl;

            break;
         }
         case 2: // Cell Center
         {
            Vector Uc(dim+2);
            GetCellStateVector(S, EDof, Uc);
            pb->velocity(Uc, node_v);

            /* Put this velocity into mv_gf */
            geom.UpdateNodeVelocity(S,node,node_v);
            
            continue;
         }
         default:
         {
            MFEM_ABORT("Invalid entity\n");
         }
      }

      // Once you have the adjacent elements to a node, the computation is the same
      // Iterate over cells that share this node and average their contribution to vertex velocity
      for (int element_it = 0; element_it < row_length; element_it++)
      {
         v_cell = 0.;
         int el_index = element_row[element_it];
         // cout << "\t node: " << node << endl;

         Array<int> oris;
         pmesh->GetElementEdges(el_index, cell_faces_row, oris);
         cell_faces_length = cell_faces_row.Size();

         // Get Transformation and IntegrationPoint
         ElementTransformation * trans = pmesh->GetElementTransformation(el_index);
         Vector node_x(dim);
         geom.GetNodePositionFromBV(S, node, node_x);
         IntegrationPoint ip;
         trans->TransformBack(node_x, ip);
         trans->SetIntPoint(&ip);

         // Get vector shape function evaluations
         DenseMatrix shape_funcs(4,2);
         const FiniteElement * CR_fe = CR.GetFE(el_index);
         CR_fe->CalcPhysVShape(*trans, shape_funcs);

         // Iterate over faces of this cell that touch vertex
         for (int face_it = 0; face_it < cell_faces_length; face_it++)
         {
            ivf = 0.;
            int face_index = cell_faces_row[face_it];

            // GetCorrectedFaceFlux(face_index, ivf); 
            GetIntermediateFaceVelocity(face_index, ivf);
            CalcOutwardNormalInt(S, el_index, face_index, n_vec);

            double coeff = n_vec * ivf;

            Vector _shape_fun(dim);
            shape_funcs.GetRow(face_it, _shape_fun);

            v_cell.Add(coeff, _shape_fun);
         }

         node_v.Add(1., v_cell);
      }
      // Average contributions
      node_v /= row_length;

      /* Put this velocity into v_geo_gf*/
      if (entity == 0) // Only for corners
      {
         SetViGeo(node, node_v);
      }

      geom.UpdateNodeVelocity(S, node, node_v);
      
   } // End vertex iterator
}


/****************************************************************************************************
* Function: ComputeGeoVRaviart2
* Parameters:
*  S    - BlockVector corresponding to nth timestep that stores mesh information, 
*         mesh velocity, and state variables.
*
* Purpose:
*  Iterate over cells and average contribution to nodes. Alternate to above
*  function.
****************************************************************************************************/
void LagrangianLOOperator::ComputeGeoVRaviart2(const Vector &S)
{
   // cout << "=======================================\n"
   //      << "          ComputeGeoVRaviart2          \n"
   //      << "=======================================\n";
   Vector n_vec(dim), ivf(dim), v_cell(dim), node_x(dim), _shape_fun(dim);
   Array<int> cell_faces_row, oris, zones_per_dof, vdofs;
   int row_length, cell_faces_length;

   // reset Vgeo
   v_geo_gf = 0.;

   // Element transformation specific
   ElementTransformation * trans;
   IntegrationPoint ip;
   int ldof;
   DenseMatrix shape_funcs(4,2);

   zones_per_dof.SetSize(NDofs_H1);
   zones_per_dof = 0;

   // Iterate over all cells and distribute its contribution to the nodes
   for (int el = 0; el < L2.GetNE(); el++)
   {
      pmesh->GetElementEdges(el, cell_faces_row, oris);
      cell_faces_length = cell_faces_row.Size();
      
      trans = pmesh->GetElementTransformation(el);
      H1.GetElementDofs(el, vdofs);

      for (int i = 0; i < vdofs.Size(); i++)
      {
         v_cell = 0.;
         ldof = vdofs[i];

         // Get integration point corresponding to dof
         geom.GetNodePositionFromBV(S, ldof, node_x);
         trans->TransformBack(node_x, ip);
         trans->SetIntPoint(&ip);

         // Get vector shape function evaluations
         const FiniteElement * CR_fe = CR.GetFE(el);
         CR_fe->CalcPhysVShape(*trans, shape_funcs);

         // Iterate over faces of this cell
         for (int face_it = 0; face_it < cell_faces_length; face_it++)
         {
            ivf = 0.;
            int face_index = cell_faces_row[face_it];

            GetCorrectedFaceFlux(face_index, ivf); 
            CalcOutwardNormalInt(S, el, face_index, n_vec);

            double coeff = n_vec * ivf;
            shape_funcs.GetRow(face_it, _shape_fun);
            v_cell.Add(coeff, _shape_fun);
         }

         /* Put this velocity into v_geo_gf*/
         for (int i = 0; i < dim; i++)
         {
            int index = ldof + i * NVDofs_H1;
            v_geo_gf[index] += v_cell[i];
         }
         
         // Keep a tally per dof of contributions for later averaging
         zones_per_dof[ldof]++;
      }
   }

   // Average contributions from all cells
   for (int i = 0; i < zones_per_dof.Size(); i++)
   {
      const int nz = zones_per_dof[i];
      if (nz) { v_geo_gf(i) /= nz; }
   }
}

} // end ns hydroLO

} // end ns mfem
