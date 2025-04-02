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

# Handling boundary conditions
## Boundary attributes
The user has the option in the problem file to impose boundary conditions on
the thermodynamics and/or on the mesh velocity.
Both of these types of boundary conditions have corresponding flags ```_thbcs```
and ```_mvbcs``` which can be set to ```true``` in the problem.h file.

Boundary conditions are implemented via a flag in the mesh file.  Common boundaries 
that can be used are

1. Used to enforce $v_x = 0$.
2. Used to enforce $v_y = 0$.
3. Used to enforce $v_z = 0$.
4. Used to enforce $v_r = 0$, or in other words 0 radial movement.
5. Used to enforce arbitrary bcs, to be handled in the problem.h file. [5+]

If one chooses to implement BC that are not some version of $v\cdot n = 0$, then
the functions ProblemBase::get_additional_BCs and ProblemBase::update_additions_BCs
must be overridden. See for example Dirichlet conditions enforced on the left wall 
of the Saltzman problem.


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

## Examples

### Elastic shocktube
```
./Laglos -m ../data/elastic/ref-segment.mesh -p 50 -tf 0.00005 -cfl 0.5 -ue -rs 8
```
### Elastic impact
```
./Laglos -m ../data/elastic/ref-segment.mesh -p 51 -tf 0.00005 -cfl 0.5 -ue -rs 8
```
### Elastic shear
```
./Laglos -m ../data/elastic/tube-100x1y.mesh -p 52 -tf 0.00005 -cfl 0.5 -ue -ppd -rs 0
./Laglos -m ../data/elastic/distube-100x1y.mesh -p 52 -tf 0.00005 -cfl 0.5 -ue -ppd -rs 0
./Laglos -m ../data/elastic/tube-2x100y.mesh -p 52 -tf 0.00005 -cfl 0.5 -ue -ppd -rs 0
```
### Elastic Isentropic Vortex
```
./Laglos -m ../data/elastic/square-vortex-mz.mesh -p 53 -tf 1 -cfl 0.5 -ue -rs 1
```
### Elastic Projectile Plate
The final time depends on the shear modulus used, which is modified in the problem file.
The final times reported in the referenced paper corresponding to their shear moduli are as follows
#### $\mu = 9.2\times 10^{10}$ Pa
1. $t = 3.6\times 10^{-5}$
2. $t = 1.06\times 10^{-4}$
3. $t = 3.2\times 10^{-4}$
4. $t = 6.09\times 10^{-4}$

#### $\mu = 1\times 10^{9}$ Pa
1. $t = 3.5\times 10^{-5}$
2. $t = 1.4\times 10^{-4}$
3. $t = 4.2\times 10^{-4}$
4. $t = 7.1\times 10^{-4}$

#### $\mu = 0$ Pa
1. $t = 7.5\times 10^{-5}$
2. $t = 1.87\times 10^{-4}$
3. $t = 6\times 10^{-4}$
4. $t = 1.03\times 10^{-4}$
```
./Laglos -m ../data/elastic/projectile-plate.mesh -p 54 -tf .000103 -cfl 0.5 -ue -rs 1
```