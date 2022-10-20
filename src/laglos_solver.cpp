#include "laglos_solver.hpp"
#include <cassert>

extern "C" {
   void __arbitrary_eos_lagrangian_lambda_module_MOD_lagrangian_lambda_arbitrary_eos(
      double *in_rhol, double *in_ul, double *in_el, double *in_pl,
      double *in_rhor, double *in_ur, double *in_er, double *in_pr,
      double *in_tol, bool *no_iter,double *lambda_maxl_out,
      double *lambda_maxr_out, double *pstar, int *k, double *b_covolume);
}

namespace mfem
{

namespace hydrodynamics
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

template<int dim>
LagrangianLOOperator<dim>::LagrangianLOOperator(ParFiniteElementSpace &h1,
                                           ParFiniteElementSpace &l2,
                                           ParFiniteElementSpace &l2v,
                                           ParLinearForm *m) :
   H1(h1),
   L2(l2),
   L2V(l2v),
   pmesh(H1.GetParMesh()),
   m_lf(m),
   Vsize_H1(H1.GetVSize()),
   TVSize_H1(H1.TrueVSize()),
   GTVSize_H1(H1.GlobalTrueVSize()),
   NDofs_H1(H1.GetNDofs()),
   Vsize_L2(L2.GetVSize()),
   TVSize_L2(L2.TrueVSize()),
   GTVSize_L2(L2.GlobalTrueVSize()),
   NDofs_L2(L2.GetNDofs()),
   Vsize_L2V(L2V.GetVSize()),
   TVSize_L2V(L2V.TrueVSize()),
   GTVSize_L2V(L2V.GlobalTrueVSize()),
   NDofs_L2V(L2V.GetNDofs()),
   block_offsets(6),
   x_gf(&H1),
   // ess_tdofs(this->ess_tdofs),
   NE(pmesh->GetNE()),
   l2dofs_cnt(L2.GetFE(0)->GetDof()),
   h1dofs_cnt(H1.GetFE(0)->GetDof()),
   l2vdofs_cnt(L2V.GetFE(0)->GetDof())

{
   block_offsets[0] = 0;
   block_offsets[1] = block_offsets[0] + Vsize_H1;
   block_offsets[2] = block_offsets[1] + Vsize_H1;
   block_offsets[3] = block_offsets[2] + Vsize_L2;
   block_offsets[4] = block_offsets[3] + Vsize_L2V;
   block_offsets[5] = block_offsets[4] + Vsize_L2;

   // Build mass vector
   m_hpv = m_lf->ParallelAssemble();
}

template<int dim>
LagrangianLOOperator<dim>::~LagrangianLOOperator()
{

}

template<int dim>
void LagrangianLOOperator<dim>::CalculateTimestep(const Vector &S)
{
   /* TODO: Fix, issue with both timestep dependence on d and d dependence on timestep */
   // int n = m_hpv->Size(); // NDofs_L2
   double t_min = 1.;
   double t_temp = 0;
   double mi = 0;
   int ci = 0, cj = 0;

   Array<int> fids, oris;
   Vector c(dim), n(dim);
   Vector U_i(dim+2), U_j(dim+2);
   double c_norm = 0., d=0., temp_sum = 0.;

   mfem::Mesh::FaceInformation FI;

   for (ci; ci < L2.GetNE(); ci++) // Cell iterator
   {
      temp_sum = 0.;
      mi = m_hpv->Elem(ci);
      // cout << "(-1-) mi: " << mi << endl;
      assert(mi > 0); // Assumption, equation (3.6)

      GetCellStateVector(S, ci, U_i);

      H1.ExchangeFaceNbrData();

      pmesh->GetElementEdges(ci, fids, oris); // 3DTODO: GetElementFaces()

      for (int j=0; j < fids.Size(); j++) // Face iterator
      {
         CalcOutwardNormalInt(S, ci, fids[j], c);

         c_norm = c.Norml2();
         n = c;
         n /= c_norm;

         FI = pmesh->GetFaceInformation(fids[j]);

         if (FI.IsInterior())
         {
            // cout << "Interior face\n";
            // Get index information/state vector for second cell
            if (ci == FI.element[0].index) { 
               cj = FI.element[1].index; 
            }
            else { 
               cj = FI.element[0].index; 
            }

            GetCellStateVector(S, cj, U_j);

            // viscosity contribution
            // cout << "Ui: \n";
            // U_i.Print(cout);
            // cout << "Uj: \n";
            // U_j.Print(cout);
            // cout << "n: \n";
            // n.Print(cout);
            // cout << "c_norm: " << c_norm << endl;
            // cout << "Computing lambda max at ci: " << ci << " and cj: " << cj << endl;   
            d = compute_lambda_max(U_i, U_j, n) * c_norm; 

            temp_sum += d/mi;
            // cout << "(-2-) d: " << d << endl;
            // cout << endl;
         }
      }
      // cout << "(-3-) temp_sum: " << temp_sum << endl << endl;
      
      t_temp = 0.5 / temp_sum;

      if (t_temp < t_min && t_temp > 1e-12) { t_min = t_temp; }
   }

   this->timestep = t_min;
}

template<int dim>
double LagrangianLOOperator<dim>::GetTimestep()
{
   return timestep;
}

template<int dim>
double LagrangianLOOperator<dim>::GetTimeStepEstimate(const Vector &S) const
{
   // Todo: implement parallelized time step estimator here.
   return 0.001;
}

template<int dim>
void LagrangianLOOperator<dim>::MakeTimeStep(Vector &S, double & t, double dt)
{
   // cout << "Carrying out a timestep\n";
   // Retrieve data from Monolithic block vector S
   Vector* sptr = const_cast<Vector*>(&S);
   ParGridFunction mv_gf;
   mv_gf.MakeRef(&H1, *sptr, block_offsets[1]);
   
   // Variables needed for iteration
   Vector val(dim+2), c(dim), n(dim), U_i(dim+2), U_j(dim+2), sums(dim+2);
   Array<int> fids, oris;

   int ci = 0, cj = 0;
   double d, c_norm;

   mfem::Mesh::FaceInformation FI;

   for (ci; ci < L2.GetNE(); ci++) // Cell iterator
   {
      H1.ExchangeFaceNbrData();
      GetCellStateVector(S, ci, U_i);
      // cout << "m: " << m_hpv->Elem(ci) << ", dt: " << dt << endl;
      val = U_i; // TODO: Must mult by m and divide by dt

      pmesh->GetElementEdges(ci, fids, oris); // 3DTODO: GetElementFaces()
      DenseMatrix F_i = flux(U_i);
      sums = 0.;
      for (int j=0; j < fids.Size(); j++) // Face iterator
      {
         CalcOutwardNormalInt(S, ci, fids[j], c);

         c_norm = c.Norml2();
         n = c;
         n /= c_norm;

         FI = pmesh->GetFaceInformation(fids[j]);

         if (FI.IsInterior())
         {
            // cout << "Interior face\n";
            // Get index information/state vector for second cell
            if (ci == FI.element[0].index) { 
               cj = FI.element[1].index; 
            }
            else { 
               cj = FI.element[0].index; 
            }

            GetCellStateVector(S, cj, U_j); // issue here

            // flux contribution
            DenseMatrix dm = flux(U_j);
            dm += F_i;
            Vector y(dim+2);
            dm.Mult(c, y);
            sums -= y;

            // viscosity contribution
            // cout << "Carrying out a timestep at ci: " << ci << " and cj: " << cj << endl;
            d = compute_lambda_max(U_i, U_j, n) * c_norm; 
            Vector z = U_j;
            z -= U_i;
            sums.Add(d, z);
         }
         else
         {
            // cout << "boundary face\n";
            assert(FI.IsBoundary());
            // TODO: Add bdry summation stuff
            Vector y_temp(dim+2);
            F_i.Mult(c, y_temp);
            y_temp *= 2;
            sums -= y_temp;
         }         

      } // End Face iterator

      sums *= dt;
      sums /= m_hpv->Elem(ci);
      val += sums;
      SetCellStateVector(S, ci, val);
      
   } // End cell iterator

   t += dt;
}

template<int dim>
void LagrangianLOOperator<dim>::GetCellStateVector(const Vector &S, 
                                              const int cell, 
                                              Vector &U)
{
   U.SetSize(dim + 2);

   // Retrieve information from BlockVector S
   Vector* sptr = const_cast<Vector*>(&S);
   ParGridFunction sv_gf, v_gf, ste_gf;
   sv_gf.MakeRef(&L2, *sptr, block_offsets[2]);
   v_gf.MakeRef(&L2, *sptr, block_offsets[3]);
   ste_gf.MakeRef(&L2, *sptr, block_offsets[4]);

   // Retrieve specific volume on cell
   U[0] = sv_gf.Elem(cell);
   // Retrieve cell velocity 
   // Here the corresponding gridfunction is a stacked vector:
   // Ex: x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4
   for (int i = 0; i < dim; i++)
   {
      U[i+1] = v_gf.Elem(cell + i*NDofs_L2V);
   }

   // Retrieve cell specific total energy
   U[dim+1] = ste_gf.Elem(cell);
}

template<int dim>
void LagrangianLOOperator<dim>::SetCellStateVector(Vector &S_new,
                                              const int cell,
                                              const Vector &U)
{
   // cout << "Setting cell state vector\n";
   Array<int> dofs;
   Array<int> sub_dofs;
   dofs.Append(cell);
   sub_dofs.Append(1);

   // Get current state grid functions
   Vector* sptr = const_cast<Vector*>(&S_new);
   ParGridFunction sv_gf, v_gf, ste_gf;
   sv_gf.MakeRef(&L2, *sptr, block_offsets[2]);
   v_gf.MakeRef(&L2, *sptr, block_offsets[3]);
   ste_gf.MakeRef(&L2, *sptr, block_offsets[4]);

   // Set subvectors according to U
   sv_gf.SetSubVector(dofs, U[0]);
   ste_gf.SetSubVector(dofs, U[dim+1]);

   for (int i = 1; i < dim; i++)
   {
      dofs.Append(cell + i*NDofs_L2V);
      sub_dofs.Append(i+1);
   }

   Vector vel;
   U.GetSubVector(sub_dofs, vel);
   v_gf.SetSubVector(dofs, vel);

   // Sync update gridfunctions
   sv_gf.SyncAliasMemory(S_new);
   v_gf.SyncAliasMemory(S_new);
   ste_gf.SyncAliasMemory(S_new);
}

/* cij computation */
/*
*
* Build normal matrix
*
*/
template<int dim>
void LagrangianLOOperator<dim>::build_C(const Vector &S, double dt)
{
   /*
   iterate over faces
   TODO: 
      - Need to do specific things for interface vs boundary faces.
      - Incorporate direction of normal somehow.
   */
   Vector* sptr = const_cast<Vector*>(&S);
   ParGridFunction mv_gf;
   mv_gf.MakeRef(&H1, *sptr, block_offsets[1]);
   Vector res(dim);
   Array<int> fids, oris;

   for (int i = 0; i < H1.GetNE(); i++) // Cell iterator
   {
      H1.ExchangeFaceNbrData();
      pmesh->GetElementEdges(i, fids, oris); // 3DTODO: GetElementFaces()
      for (int j=0; j < fids.Size(); j++) // Face iterator
      {
         CalcOutwardNormalInt(S, i, fids[j], res);
      } // End Face iterator
   } // End cell iterator
}

template<int dim>
void LagrangianLOOperator<dim>::CalcOutwardNormalInt(const Vector &S, const int cell, const int face, Vector & res)
{
   Vector* sptr = const_cast<Vector*>(&S);
   ParGridFunction x_gf, mv_gf;
   x_gf.MakeRef(&H1, *sptr, block_offsets[0]);
   // Fill nodes with mesh nodes for retrieval
   const int nNodes = x_gf.Size() / dim; // TODO: Implement class variables for this

   Array<int> face_dofs;
   Vector a(dim);
   Vector face_node(dim); // Temp
   res = 0.;

   // Get int der shape functions
   Vector shapes = GetIntDerRefShapeFunctions();
   assert (dim == 2); // "This code only works for dim==2."

   H1.GetFaceDofs(face, face_dofs);

   // Compute C_{cc'} = 1/2 \Sum_{i\in {1:3}} [a_i^n \int_0^1 \theta_i(\xi) d\xi]
   for (int d = 0; d < face_dofs.Size(); d++)
   {
      a = 0.;
      for (int _dim = 0; _dim < dim; _dim++)
      {
         int index = _dim * nNodes + face_dofs[d];
         a[_dim] = x_gf(index);
      }
      if (d == face_dofs.Size() - 1)
      {
         face_node = a;
      }

      // Now that vectors are retrieved, sum
      a *= shapes[d];
      add(res, a, res);
   }

   Orthogonal(res);
   res *= 0.5;

   /* Ensure orientation of normal */
   const auto FI = pmesh->GetFaceInformation(face);

   if (FI.IsInterior()) // Interior face
   {
      // Orientation of the normal vector depends on the indices of
      // the cells that share that face.  The normal points towards
      // the cell with the lower global index.  This normal must be
      // flipped if the cell we are working with is the one of the
      // lower index.
      if (cell == min(FI.element[0].index, FI.element[1].index))
      {
         res *= -1.;
      }
   }
   else // Boundary face
   {
      // By default, boundary normals point toward cell.  This always
      // needs to be reversed to compute outward normal.
      assert(FI.IsBoundary());
      res *= -1.;
   }

   // Output normal
   // {
   //    cout << "Face node coords: ";
   //    face_node.Print(cout);
   //    cout << "Normal vector: ";
   //    res.Print(cout);
   // }
}

template<int dim>
void LagrangianLOOperator<dim>::Orthogonal(Vector &v)
{
   if (dim == 2)
   {
      double x = v(0), y = v(1);
      v(0) = -1 * y;
      v(1) = x;
   }
}

template<int dim>
Vector LagrangianLOOperator<dim>::GetIntDerRefShapeFunctions() 
{
   /* 
   Test stuff 
   Can use this to evaluate shape functions on the reference face.
   */
   // const FiniteElement * test_fe_p = H1FEC.FiniteElementForGeometry(Geometry::SEGMENT); 
   // cout << "FE derive type: " << test_fe_p->GetDerivType() << endl;
   // const IntegrationRule ir = test_fe_p->GetNodes();
   // cout << "num integration points: " << ir.GetNPoints() << endl;

   // for (int p = 0; p < ir.GetNPoints(); p++)
   // {
   //    const IntegrationPoint &ip = ir.IntPoint(p);
   //    cout << "ip x: " << ip.x << ", ip y: " << ip.y << ", ip weight: " << ip.weight << endl;
   //    DenseMatrix dm(3,2);
   //    test_fe_p->CalcDShape(ip, dm);
   //    dm.Print(cout);
   // }
   assert (dim == 2); // "This code only works for dim==2."

   Vector r(3);
   r[0] = -1.;
   r[1] = 1.;
   r[2] = 0.;

   return r;
}

/* ProblemDescription Functions */
template<int dim>
double LagrangianLOOperator<dim>::internal_energy(const Vector & U)
{
   const double &rho = 1./U[0];
   const double &e = specific_internal_energy(U);
   return rho * e;
}

template<int dim>
double LagrangianLOOperator<dim>::specific_internal_energy(const Vector & U)
{
   const Vector  v   = velocity(U);
   const double &E   = U[dim + 1]; // specific total energy
   return E - 0.5 * pow(v.Norml2(), 2);
}

template<int dim>
double LagrangianLOOperator<dim>::pressure(const Vector & U)
{
   return (gamma - 1.) * internal_energy(U);
}

template<int dim>
double LagrangianLOOperator<dim>::compute_lambda_max(const Vector & U_i,
                                                     const Vector & U_j,
                                                     const Vector & n_ij,
                                                     const string flag)
{
   double in_rhol, in_ul, in_el, in_pl, in_rhor, in_ur, in_er, in_pr;
   if (flag == "testing")
   {
      in_rhol = U_i[0];
      in_ul = U_i[1];
      in_pl = U_i[2];
      in_el = U_i[3];

      in_rhor = U_j[0]; 
      in_ur = U_j[1];
      in_pr = U_j[2];
      in_er = U_j[3];
   }
   else 
   {
      assert(flag == "NA");
      in_rhol = 1. / U_i[0];
      in_ul = U_i * n_ij; // TODO: Validate.
      in_el = specific_internal_energy(U_i);
      in_pl = pressure(U_i);

      in_rhor = 1. / U_j[0]; 
      in_ur = U_j * n_ij; // TODO: Check that this shouldn't be n_ji.
      in_er = specific_internal_energy(U_j);
      in_pr = pressure(U_j);
   }
   // std::cout << "Computing Viscosity\n";
   double b_covolume = 0.1/max(in_rhol,in_rhor);

   double in_tol,lambda_maxl_out,lambda_maxr_out,pstar;
   bool no_iter = true;
   int k = 10;
   __arbitrary_eos_lagrangian_lambda_module_MOD_lagrangian_lambda_arbitrary_eos(
      &in_rhol,&in_ul,&in_el,&in_pl,&in_rhor,&in_ur,&in_er,&in_pr,&in_tol,
      &no_iter,&lambda_maxl_out,&lambda_maxr_out,&pstar,&k,&b_covolume);

   return std::max(std::abs(lambda_maxl_out), std::abs(lambda_maxr_out));
}

template<int dim>
Vector LagrangianLOOperator<dim>::velocity(const Vector & U)
{
   Vector v;
   const double &rho = 1./U[0];
   v.SetSize(dim);
   Array<int> dofs;
   for (int i = 0; i < dim; i++)
   {
      dofs.Append(i + 1);
   }
   U.GetSubVector(dofs, v);

   return v;
}

template<int dim>
DenseMatrix LagrangianLOOperator<dim>::flux(const Vector &U)
{
   DenseMatrix result(dim+2, dim);

   const Vector v = velocity(U);
   const double p = pressure(U);

   // * is not overridden for Vector class, but *= is
   Vector v_neg = v, vp = v;
   v_neg *= -1.;
   vp *= p;

   // Set f(U) according to (2.1c)
   result.SetRow(0,v_neg);

   for (int i = 0; i < dim; i++)
   {
      result(i+1, i) = p;
   }
   result.SetRow(dim+1, vp);

   return result;
}

/* Explicit Instantiation */
template class LagrangianLOOperator<1>;
template class LagrangianLOOperator<2>;
template class LagrangianLOOperator<3>;

} // end ns hydrodynamics


} // end ns mfem