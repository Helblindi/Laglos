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

# Select script insructions [located in ${source}/scripts directory]
## Creating convergence tables
> To show proper convergence of the utilized numerical method, one can construct a table of multiple error rates from a given test problem.  To valide the numerical method, a convergence order can then be computed from two different refinements.


## Organize output pvtu files
> When running a test problem with the flags `-pview -print -vs #`, pvtu and vtu files
will be generated and placed in individual directories according to their timesteps.
In an effort to make importing all these files faster in a visualization software 
such as ParaView, it is helpful for all pvtu files to be in the same directory.  
The following command (executed from the ${source}/scripts directory) will handle 
renaming all pvtu files and replacing references to their corresponding vtu files
for each processor.  
```
$ python3 reorganize_paraview.py <input_dir> <output_dir>
````

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


# Select Elastic notes
Since our exact solution here is given in terms of approximation data from Dr. Favrie,
we cannot compute errors the conventional way.  To this end, I have created the file
compute_errors to handle both the computation of the approximation errors as well 
as convergence tables. To execute this, run the following
```
$ python3 ../scripts/compute_error.py results/testing/st1/state_vectors/ ../exact_sol/elastic/shocktube/all_results.txt 7 16
```
The third and fourth arguments correspond to matching columns of the respective files that are being compared. We list here the following pairs and their corresponding 
representative value

| parameter | arg 3 | arg 4 |
|:----------|:-----:|:-----:|
| $\rho$    |   1   |   13  |
| $E$       | 3     | 15    |
| $\sigma$  | 7     | 16    |

## Meshes
Laglos has the capability to handle multimaterial test problems, and this is implemented through the 
cell attribute values defined in the mesh. For each element that should be treated as a solid,
the cell attribute in the mesh should be set to 50. When the use-elasticity ['-ue'] is used in 
a Laglos execution, the elastic flux and elastic sheer with be computed only if the cell
attribute value is set to 50. Otherwise, a non-elastic flux will be used.

To see this implemented, see the execution of the multi-material isentropic vortex problem
```
./Laglos -m data/elastic/square-vortex-mz.mesh -p 53 -tf 10 -cfl 0.5 -ue 
```