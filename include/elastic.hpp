//NF//MS
/*
 * This is the elastic class, it is used to compute the elastic energy 
 * and the deviatoric part of the stress tensor.
 */
 #ifndef ELASTIC
 #define ELASTIC
 
 #include "mfem.hpp"
 #include "problem_base.h"
 
 using namespace std;
 
 namespace mfem
 {
 
 namespace hydroLO
 {
 
 enum ShearEOS {
    NEO_HOOKEAN,
    MOONEY_RIVLIN,
    AORTIC         // anisotropic model
 };
 
 enum ShearEnergyMethod {
    AVERAGE_F,
    AVERAGE_C,
    AVERAGE_ENERGY
 };
 
 class Elastic
 {
 private:
    const int dim;
    ShearEnergyMethod shear_method;
    ShearEOS shear_eos;
    ParFiniteElementSpace &H1, &L2;
    const ParGridFunction &rho0_gf;
    const IntegrationRule &ir;
    const int NE, NDofs_L2, nqp;
 
    // Reference to physical Jacobian for the initial mesh.
    // These are computed only at time zero and stored here.
    DenseTensor Jac0inv;
    double mu = -1.; // Shear modulus
 
    /* 
    If using the aortic eos, a fiber direction must be defined. 
    This direction is used to compute the invariants:
       m = (cos(phi), sin(phi), 0) 
    */
    double phi = 0.; // angle of fiber direction
    Vector mi;
    DenseMatrix Gi;
    double A1 = 771.8, B1 = 21.2, D1 = 3.8, w1 = 0.4971;
 
 public:
    Elastic(const int &_dim,
            ParFiniteElementSpace &h1_fes,
            ParFiniteElementSpace &l2_fes,
            const ParGridFunction &rho0_gf,
            const IntegrationRule &ir,
            ShearEOS _shear_eos = ShearEOS::MOONEY_RIVLIN,
            ShearEnergyMethod method = ShearEnergyMethod::AVERAGE_C) : 
       dim(_dim),
       H1(h1_fes), 
       L2(l2_fes),
       rho0_gf(rho0_gf),
       ir(ir),
       NE(H1.GetParMesh()->GetNE()),
       NDofs_L2(L2.GetNDofs()), 
       nqp(ir.GetNPoints()),
       Jac0inv(dim, dim, NE * nqp),
       shear_eos(_shear_eos),
       shear_method(method),
       mi(3),
       Gi(3)
    {
       /*************** CONFIGURATION OUTPUT ***************/
       cout << "\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
       << "@         Elastic Class Configuration      @\n"
       << "@------------------------------------------@\n"
       << "@ shear_eos         : " << std::setw(20) << std::left;
       switch (shear_eos)
       {
          case NEO_HOOKEAN:
             cout << "NEO_HOOKEAN";
             break;
          case MOONEY_RIVLIN:
             cout << "MOONEY_RIVLIN";
             break;
          case AORTIC:
             cout << "AORTIC";
             break;
          default:
             MFEM_ABORT("Invalid value for shear_eos.");
       }
       cout << " @\n";
       if (shear_eos == ShearEOS::AORTIC)
       {
          cout << "@ Aortic parameters : " << std::setw(20) << std::left
               << "@ phi               : " << std::setw(20) << std::left << phi << " @\n"
               << "@ m                : " << std::setw(20) << std::left;
          for (int i = 0; i < dim; i++)
          {
             cout << mi(i) << " ";
          }
          cout << " @\n";
       }
 
       
       cout << "@ ShearEnergyMethod : " << std::setw(20) << std::left;
       switch (shear_method)
       {
          case ShearEnergyMethod::AVERAGE_F:
             cout << "AVERAGE_F";
             break;
          case ShearEnergyMethod::AVERAGE_C:
             cout << "AVERAGE_C";
             break;
          default:
             MFEM_ABORT("Unknown shear energy method");
       }
       cout << " @\n";
       
       cout << "@ nqp               : " << std::setw(20) << std::left << nqp      << " @\n"
             << "@ NDofs_L2          : " << std::setw(20) << std::left << NDofs_L2 << " @\n"
             << "@ NE                : " << std::setw(20) << std::left << NE       << " @\n"
             << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n\n";
       /***************** END CONFIGURATION *****************/
 
       if (NDofs_L2 != NE)
       {
          MFEM_ABORT("Error: Number of L2 dofs must be equal to the number of elements.\n");
       }
 
       /* Construct initial inverse jacobian transformations*/
       Jac0inv = 0.;
       for (int e = 0; e < NE; e++)
       {
          ElementTransformation *T = H1.GetElementTransformation(e);
          for (int q = 0; q < nqp; q++)
          {
             const IntegrationPoint &ip = ir.IntPoint(q);
             T->SetIntPoint(&ip);
             DenseMatrixInverse Jinv(T->Jacobian());
             Jinv.GetInverseMatrix(Jac0inv(e*nqp + q));
          }
       }
 
       /* Compute anisotropic invariants */
       if (shear_eos == ShearEOS::AORTIC)
       {
          /* Set fiber direction */
          mi(0) = cos(phi);
          mi(1) = sin(phi);
          mi(2) = 0.;
          tensor(mi, mi, Gi);
       }
    }
 
    void tensor(const Vector & v1, const Vector & v2, DenseMatrix & dm) const
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
    void set_shear_modulus(const double &_mu) { this->mu = _mu; }
 
    void setShearEnergyMethod(ShearEnergyMethod method) { shear_method = method; }
 
    double e_sheer(const int &e) const
    {
       DenseMatrix c(3), c2(3);
       DenseMatrix cF(3), cc(3);
       const double rho0 = rho0_gf(e);
       switch (shear_method) {
          case ShearEnergyMethod::AVERAGE_F:
          {
             Compute_cFromAvgF(e, c);
             break;
          }
          case ShearEnergyMethod::AVERAGE_C:
          {
             Compute_cAvg(e,c);
             break;
          }
          case ShearEnergyMethod::AVERAGE_ENERGY:
             MFEM_ABORT("Not implemented");
          default:
             MFEM_ABORT("Unknown shear energy method");
       }
       return e_sheer(c, rho0);
    }
    
    double e_sheer(const DenseMatrix &c, const double &rho0) const
    {
       /* Compute trace values */
       double trc = c.Trace();
       DenseMatrix c2(3);
       mfem::Mult(c,c,c2);
       double trc2 = c2.Trace();
 
       if (mu == -1.)
       {
          MFEM_ABORT("Must set shear modulus.\n");
       }
 
       switch (shear_eos)
       {
       case NEO_HOOKEAN:
          return mu / 2. * (trc - 3.) / rho0;
       case MOONEY_RIVLIN:
          return mu / 32. / rho0 * (pow(trc,4) - 2*trc2*pow(trc,2) + pow(trc2,2) - 8 * trc - 12.);
          // return 0.;
       case AORTIC:
       {
          /* compute i4 and i5 */
          DenseMatrix _res(3);
          Mult(c, Gi, _res);
          double i4 = _res.Trace();
          Mult(c, c, _res);
          Mult(_res, Gi, _res);
          double i5 = _res.Trace();
 
          /* Compute shear energy */
          double val = (1. - 2 * w1) * ((A1 / 2.) * (trc - 3.) + (B1 / 4.) * (trc*trc - trc2 - 6.));
          val += 2. * w1 * D1 * ((1. / A1) * (exp(A1 * (i4 - 1.)) - 1.) + (1. / B1) * (exp(B1 * (i5 - trc * i4 + (trc*trc - trc2) / 2. - 1.)) - 1.));
          return val;
       }
       default:
          MFEM_ABORT("Invalid value for shear_eos.");
       }
        
       /* alternate eos from paper */
       // double j2 = pow(trc2 - pow(trc,2),2) - 2. * trc;
       // return mu * (j2 - 3.) / 8. / rho0;
       /* favrie 2014 */
       // return mu * (trc2 - 3.) / 8 / rho0;
 
       return -1.;
    }
 
    double des_dtrc(const DenseMatrix &c, const double &rho0) const
    {
       /* Compute trace values */
       const double trc = c.Trace();
       DenseMatrix c2(3);
       mfem::Mult(c,c,c2);
       const double trc2 = c2.Trace();
 
       if (mu == -1.)
       {
          MFEM_ABORT("Must set shear modulus.\n");
       }
 
       switch (shear_eos)
       {
       case NEO_HOOKEAN:
          return mu / 2. / rho0;
       case MOONEY_RIVLIN:
          return mu / 8. / rho0 * (pow(trc,3) - trc*trc2 - 2.);
       case AORTIC:
       {
          /* compute i4 and i5 */
          DenseMatrix _res(3);
          Mult(c, Gi, _res);
          double i4 = _res.Trace();
          Mult(c, c, _res);
          Mult(_res, Gi, _res);
          double i5 = _res.Trace();
 
          double val = (1. - 2. * w1) * (A1 / 2. + (B1 / 2.) * trc);
          val += 2. * w1 * D1 * exp(B1 * (i5 - trc * i4 + (trc*trc - trc2) / 2. - 1.)) * (-i4 + trc);
          return val;
       }
       default:
          MFEM_ABORT("Invalid value for shear_eos.");
       }
 
       /* alternate eos from paper */
       // return mu * ( -4. * trc2 * trc + 4. * pow(trc,3) - 2. ) / 8. / rho0;
       /* favrie 2014 */
 
       return -1.;
    }
 
    double des_dtrc2(const DenseMatrix &c, const double &rho0) const
    {
       /* Compute trace values */
       const double trc = c.Trace();
       DenseMatrix c2(3);
       mfem::Mult(c,c,c2);
       const double trc2 = c2.Trace();
 
       if (mu == -1.)
       {
          MFEM_ABORT("Must set shear modulus.\n");
       }
 
       switch (shear_eos)
       {
       case NEO_HOOKEAN:
          return 0.;
       case MOONEY_RIVLIN:
          return mu / 16. / rho0 * (-pow(trc,2) + trc2);
          case AORTIC:
          {
             /* compute i4 and i5 */
             DenseMatrix _res(3);
             Mult(c, Gi, _res);
             double i4 = _res.Trace();
             Mult(c, c, _res);
             Mult(_res, Gi, _res);
             double i5 = _res.Trace();
    
             double val = (2. * w1 - 1.) * (B1 / 4.);
             val -= w1 * D1 * exp(B1 * (i5 - trc * i4 + (trc*trc - trc2) / 2. - 1.));
             return val;
          }
       default:
          MFEM_ABORT("Invalid value for shear_eos.");
       }
 
       /* alternate eos from paper */
       // return mu * (2. * trc2 - 2. * pow(trc,2)) / 8. / rho0;
       // return mu / 8. / rho0;
 
       return -1.;
    }
 
    double des_di4(const DenseMatrix &c) const
    {
       /* Compute trace values */
       const double trc = c.Trace();
       DenseMatrix c2(3);
       mfem::Mult(c,c,c2);
       const double trc2 = c2.Trace();
 
       assert(shear_eos == ShearEOS::AORTIC);
 
       if (mu == -1.)
       {
          MFEM_ABORT("Must set shear modulus.\n");
       }
 
       /* compute i4 and i5 */
       DenseMatrix _res(3);
       Mult(c, Gi, _res);
       double i4 = _res.Trace();
       Mult(c, c, _res);
       Mult(_res, Gi, _res);
       double i5 = _res.Trace();
 
       double val = 2. * w1 * D1 * (exp(A1 * (i4 - 1.)) - trc * exp(B1 * (i5 - trc * i4 + (trc*trc - trc2) / 2. - 1.)));
       return val;
    }
 
    double des_di5(const DenseMatrix &c) const
    {
       /* Compute trace values */
       const double trc = c.Trace();
       DenseMatrix c2(3);
       mfem::Mult(c,c,c2);
       const double trc2 = c2.Trace();
 
       assert(shear_eos == ShearEOS::AORTIC);
 
       if (mu == -1.)
       {
          MFEM_ABORT("Must set shear modulus.\n");
       }
 
       /* compute i4 and i5 */
       DenseMatrix _res(3);
       Mult(c, Gi, _res);
       double i4 = _res.Trace();
       Mult(c, c, _res);
       Mult(_res, Gi, _res);
       double i5 = _res.Trace();
 
       double val = 2. * w1 * D1 * exp(B1 * (i5 - trc * i4 + (trc*trc - trc2) / 2. - 1.));
       return val;
    }
 
    void ComputeAvgF(const int e, DenseMatrix &F) const
    {
       DenseMatrix F_dim(dim); F_dim = 0.;
       ElementTransformation *T = H1.GetElementTransformation(e);
       for (int q = 0; q < nqp; q++)
       {
          const IntegrationPoint &ip = ir.IntPoint(q);
          T->SetIntPoint(&ip);
          DenseMatrix Jr = T->Jacobian();
          /* Instead of a mapping from the reference element,
             we need a mapping from the initial configuration
             of the given element */
          DenseMatrix Ji(dim);
          mfem::Mult(Jr, Jac0inv(e*nqp + q), Ji);
          F_dim.Add(ip.weight, Ji);
          // TODO: Try this averaging on F F^T instead of F.
       }
       /* If dim != 3, will need to augment F with 1s on the diagonal */
       if (dim == 3)
       {
          F = F_dim;
       }
       else {
          /* Set submatrix to F_dim */
          F.SetSize(3);
          F = 0.;
          Array<int> idx(dim);
          for (int i = 0; i < dim; i++)
          {
             idx[i] = i;
          }
          F.SetSubMatrix(idx, idx, F_dim);
 
          /* Fill remaining entries with a 1 */
          for (int i = dim; i < 3; i++)
          {
             F(i,i) = 1.;
          }
       }
    }
 
    void Compute_cFromAvgF(const int &e, DenseMatrix &c) const
    {
       DenseMatrix F(3), FT(3), C(3);
       c.SetSize(3);
       ComputeAvgF(e,F);
       FT.Transpose(F);
       mfem::Mult(F, FT, C);
       c = C;
       c *= std::pow(C.Det(), -1./3.);
    }
 
    void Compute_cAvg(const int &e, DenseMatrix &c) const
    {
       // cout << "!!!!!!!!!!!!!!!!!Elastic::Compute_cAvg!!!!!!!!!!!!!!!!!\n";
       DenseMatrix F_dim(dim);
       DenseMatrix F(3), FT(3), FTF(3);
       DenseMatrix C(3), cqp(3);
       ElementTransformation *T = H1.GetElementTransformation(e);
 
       c.SetSize(3);
       c = 0.;
       for (int q = 0; q < nqp; q++)
       {
          F_dim = 0.;
          F = 0.; FT = 0.; FTF = 0.; C = 0.;
          const IntegrationPoint &ip = ir.IntPoint(q);
          T->SetIntPoint(&ip);
          DenseMatrix Jr = T->Jacobian();
          /* Instead of a mapping from the reference element,
             we need a mapping from the initial configuration
             of the given element */
          DenseMatrix Ji(dim);
          mfem::Mult(Jr, Jac0inv(e*nqp + q), F_dim);
          if (dim == 3)
          {
             F = F_dim;
          }
          else { /* If dim != 3, will need to augment F with 1s on the diagonal */
             /* Set submatrix to F_dim */
             Array<int> idx(dim);
             for (int i = 0; i < dim; i++)
             {
                idx[i] = i;
             }
             F.SetSubMatrix(idx, idx, F_dim);
 
             /* Fill remaining entries with a 1 */
             for (int i = dim; i < 3; i++)
             {
                F(i,i) = 1.;
             }
          }
 
          /* Compute FTF at quadrature point, and normalize */
          FT.Transpose(F);
          mfem::Mult(F, FT, cqp);
          cqp *= std::pow(cqp.Det(), -1./3.);
 
          /* Add in contribution from quadrature point */
          c.Add(ip.weight, cqp);
       }
    }
 
    void ComputeS(const int &e, const double &rho, DenseMatrix &S) const
    {
       /* Objects needed in function */
       DenseMatrix c(3), c2(3), I(3);
       for (int i = 0; i < 3; i++) { I(i,i) = 1.; }
 
       S.SetSize(3);
       S = 0.;
 
       /* Compute c */
       switch (shear_method) {
          case ShearEnergyMethod::AVERAGE_F:
          {
             Compute_cFromAvgF(e, c);
             break;
          }
          case ShearEnergyMethod::AVERAGE_C:
          {
             Compute_cAvg(e,c);
             break;
          }
          case ShearEnergyMethod::AVERAGE_ENERGY:
             MFEM_ABORT("Not implemented");
          default:
             MFEM_ABORT("Unknown shear energy method");
       }
 
       /* Compute sheer energy, and save in class member */
       mfem::Mult(c, c, c2);
       const double rho0 = rho0_gf(e);
       
       /* Compute deviatoric part of stess tensor */
       DenseMatrix c_tf(3), c2_tf(3); // 'trace-free objects'
       // remove multiple of trace from diagonal
       mfem::Add(c, I, -1./3. * c.Trace(), c_tf);
       mfem::Add(c2, I, -1./3. * c2.Trace(), c2_tf);
 
       double c_coeff = des_dtrc(c, rho0);
       double c2_coeff = 2.*des_dtrc2(c, rho0);
       mfem::Add(c_coeff, c_tf, c2_coeff, c2_tf, S);
 
       /* Handle anistropic contribution */
       if (shear_eos == ShearEOS::AORTIC)
       {
          double c4_coeff = des_di4(c);
          double c5_coeff = des_di5(c);
          DenseMatrix c4_tf(3), c5_tf(3);
          MFEM_ABORT("Function must be modified to take in F instead of c??");
       }
 
       /* Finally, multiply by 2rho */
       S *= 2. * rho;
    }
 
 };
 } // end ns hydroLO
 
 } // end ns mfem
 
 #endif // ELASTIC