#include "laglos_solver.hpp"

// Fortran subroutine
extern "C" {
   void __arbitrary_eos_lagrangian_lambda_module_MOD_lagrangian_lambda_arbitrary_eos(
      double *in_rhol, double *in_ul, double *in_el, double *in_pl,
      double *in_rhor, double *in_ur, double *in_er, double *in_pr,
      double *in_tol, bool *no_iter,double *lambda_maxl_out,
      double *lambda_maxr_out, double *pstar, int *k);
}

namespace mfem
{

LagrangianLOOperator::LagrangianLOOperator(const int size,
                                           ParFiniteElementSpace &h1,
                                           ParFiniteElementSpace &l2,
                                           ParFiniteElementSpace &l2v) :
   TimeDependentOperator(size),
   H1(h1),
   L2(l2),
   L2V(l2v),
   pmesh(H1.GetParMesh()),
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
   dim(pmesh->Dimension()),
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
}

LagrangianLOOperator::~LagrangianLOOperator()
{

}

void LagrangianLOOperator::Mult(const Vector &S, Vector &S_new) const
{

}

void LagrangianLOOperator::GetCellStateVector(const Vector &S, 
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
   U[dim+2] = ste_gf.Elem(cell);
}

void LagrangianLOOperator::SetCellStateVector(Vector &S_new,
                                              const int cell,
                                              const Vector &U)
{
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

double LagrangianLOOperator::ComputeViscosity()
{
   std::cout << "Computing Viscosity\n";
   double in_rhol,in_ul,in_el,in_pl,in_rhor,in_ur,in_er,in_pr,in_tol,lambda_maxl_out,lambda_maxr_out,pstar;
   bool no_iter;
   int k;
   __arbitrary_eos_lagrangian_lambda_module_MOD_lagrangian_lambda_arbitrary_eos(
      &in_rhol,&in_ul,&in_el,&in_pl,&in_rhor,&in_ur,&in_er,&in_pr,
      &in_tol,&no_iter,&lambda_maxl_out,&lambda_maxr_out,&pstar,&k);

   return 1.0;
}


} // end ns mfem