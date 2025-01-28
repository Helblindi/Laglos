# Lagrangian
> ## Optional runtime parameters
> - -mv - mesh velocity option
> - -fv - face velocity option 
> - -mv-iter-n - number of iterations that should be used for mesh velocity calculations
> - -m - mesh file to use 
> - -p - problem to run
> - -rs - number of serial refinements
> - -rp - number of parallel refinements (parallel not implemented)
> - -of - designate the output folder where results should be saved to
> - -print -
> - -mm
> - -visc
> - -greedy - Use greedy viscosity described in GMPST2024.
> - -ggn - Number of GMV steps to take before using Greedy lambda
> - -ms - maximum steps that the solver should take 
> - -vs - the number of steps that should pass between visualizing the solution
> - -vis - enable visualization
> - -visit - Enable outputting visit files
> - -gfprint
> - -tf - designate final time 
> - -ti - initial time (does not work for all test cases)
> - -cfl - Courant-Friedrichs-Lewy condition
> - -cm - check mesh
> - -do-mv-lin - enable mv linearization
> - -tvc - Target velocity optimization viscosity coefficient parameter
> - -ppd - Post process the density to give mass conservative update


## Mass Error Calculation
All calculations related to the local mass error are calculated in the function LagrangianLOOperator<dim>::CheckMassConservation.
The mass error is a relative quantity and is defined as 
$$ \text{error}_{\text{mass}} = \frac{\sum_{c \in \eta^{\text{Cel}}} \frac{\left|K_c^n\right|}{T_c^n} - m_c^{\rho}}{\sum_{c \in \eta^{\text{Cel}}}m_c^{\rho}} $$ 
where $m_c^{\rho}$ is the initial mass of a cell, $\left|K_c^n\right|$ is the measure of the cell at time $t^n$, and $T_c^n$ is the specific volume of a cell at time $t^n$.

The percentage of cells in which mass loss is broken is simply the number of cells in which $\frac{\left|K_c^n\right|}{T_c^n} - m_c^{\rho} > 10^{-12}$ divided by the total number of cells in the mesh.

# Select HiOp Notes
> ## SolveStatus
> - 0 - Solve successful,
> - 1 - Solve successful due to RelTol
> - 2 - Solve acceptable level
> - 5 - Infeasible problem
> - 10 - Max iteration exceeded 
> - -16 - Error in feasiblity restoration
>
> A more complete list of hiopSolveStatus values is available in the file src/Interface/hiopInterface.hpp where the enum hiopSolveStatus is defined.
