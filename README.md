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
> - -print
> - -mm
> - -visc - type of viscosity to use
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
All calculations related to the local mass error are calculated in the function LagrangianLOOperator::ValidateMassConservation.
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


## Building Laglos

Laglos has the following external dependencies:

<!-- - *hypre*, used for parallel linear algebra, we recommend version 2.11.2<br>
  https://github.com/hypre-space/hypre/releases/tag/v2.11.2 

- METIS, used for parallel domain decomposition (optional), we recommend [version 4.0.3](https://github.com/mfem/tpls/blob/gh-pages/metis-4.0.3.tar.gz) <br>
  https://github.com/mfem/tpls -->

- COINHSL, collection of linear algebra libraries bundled for use with IPOPT and other applications <br>
  http://hsl.rl.ac.uk/ipopt

- HIOP, HPC solver for nonlinear optimization problems, its GitHub master branch <br>
  https://github.com/LLNL/hiop

- MFEM, used for (high-order) finite element discretization, its GitHub master branch <br>
  https://github.com/mfem/mfem

To build the miniapp, first download *CoinHSL*, *HIOP*, and *MFEM* from the links above
and put everything on the same level as the `Laglos` directory:
```sh
~> ls
Laglos/  coinhsl-x.y.z.tar.gz  hiop/ mfem/
```

<!-- Build *hypre*:
```sh
~> tar -zxvf v2.11.2.tar.gz
~> cd hypre-2.11.2/src/
~/hypre-2.11.2/src> ./configure --disable-fortran
~/hypre-2.11.2/src> make -j
~/hypre-2.11.2/src> cd ../..
~> ln -s hypre-2.11.2 hypre
```
For large runs (problem size above 2 billion unknowns), add the
`--enable-bigint` option to the above `configure` line.

Build METIS:
```sh
~> tar -zxvf metis-4.0.3.tar.gz
~> cd metis-4.0.3
~/metis-4.0.3> make
~/metis-4.0.3> cd ..
~> ln -s metis-4.0.3 metis-4.0
```
This build is optional, as MFEM can be build without METIS by specifying
`MFEM_USE_METIS = NO` below. -->

Build COINHSL:

1. Go to http://hsl.rl.ac.uk/ipopt
2. Submit a registration form (Should respond within one working business day). Reponse email should contain a link to download CoinHSL
3. Download the tarball and unpack into workspace directory, then establish symlink

```
~/Workspace> gunzip coinhsl-x.y.z.tar.gz
~/Workspace> tar xf coinhsl-x.y.z.tar
~/Workspace> ln- -s coinhsl-x.y.z coinhsl
~/Workspace/coinhsl> cd coinhsl

~/Workspace/coinhsl> meson setup build -Dcpp_args="-fPIC" -Dc_args="-fPIC" -Dfortran_args="-fPIC" --buildtype=release --prefix={absolute dir to install in coinhsl}
~/Workspace/coinhsl> meson compile -C build
~/Workspace/coinhsl> meson install -C build
```

Build HIOP:
```
hiop/build> cmake -DCMAKE_INSTALL_PREFIX=../install
-DHIOP_COINHSL_DIR=../../coinhsl/install 
-DCMAKE_POSITION_INDEPENDENT_CODE=ON 
-S ../
hiop/build> make -j 8
hiop/build> make test -j 8
hiop/build> make install
```

Clone and build the parallel version of MFEM:
```sh
~> git clone https://github.com/mfem/mfem.git ./mfem
~> cd mfem/
~/mfem> git checkout master
~/mfem> mkdir build
~/mfem> cd build
~/mfem/build> cmake -DCMAKE_INSTALL_PREFIX=./install -DCMAKE_POSITION_INDEPENDENT_CODE=ON -S ../
~/mfem/build> make parallel -j
```
The above uses the `master` branch of MFEM.
See the [MFEM building page](http://mfem.org/building/) for additional details.

(Optional) Clone and build GLVis:
```sh
~> git clone https://github.com/GLVis/glvis.git ./glvis
~> cd glvis/
~/glvis> make
~/glvis> cd ..
```
The easiest way to visualize Laglos results is to have GLVis running in a
separate terminal. Then the `-vis` option in Laglos will stream results directly
to the GLVis socket.

Build Laglos
```sh
~> cd Laglos/
~/Laglos> mkdir build
~/Laglos> cd build
~/Laglos/build> cmake -DCMAKE_INSTALL_PREFIX=./install -DCMAKE_POSITION_INDEPENDENT_CODE=ON ../
~/Laglos/build> make
~/Laglos/build> make install
```
This can be followed by `make test` and `make install` to check and install the
build respectively. See `make help` for additional options.

See also the `make setup` target that can be used to automated the
download and building of hypre, METIS and MFEM.


### Notes on build
```sh
mfem/build$ cmake -DCMAKE_INSTALL_PREFIX=./install -DCMAKE_POSITION_INDEPENDENT_CODE=ON ..
```

The flag `-DCMAKE_POSITION_INDEPENDENT_CODE=ON` is important on Linux systems as it handles shared library linking and position-independent code (PIC) different from mac.

On Linux, static libraries (.a) do not include position-independent code (-fPIC) unless explicitly compiled with it. When you try to link a static library (libmfem.a) that was not built with -fPIC into a shared library (libLaglos.so), the linker (ld) throws an error because it cannot use relocatable code in a shared object.

TODO: Add section of optional build if optimization is to be utilized.

## Build Laghos to limit

```sh
~/Workspace> cd Laghos
~/Workspace/Laghos> mkdir build
~/Workspace/Laghos> cd build
~/Workspace/Laghos/build> cmake -DLaglos_DIR=../../Laglos -S .. 
~/Workspace/Laghos/build> make -j 8
```

# Handling boundary conditions
## Boundary attributes
The user has the option in the problem file to impose boundary conditions on
the thermodynamics and/or on the mesh velocity.
Both of these types of boundary conditions have corresponding flags ```_thbcs```
and ```_mvbcs``` which can be set to ```true``` in the problem.h file.

Boundary conditions are implemented via a flag in the mesh file.  Common boundaries 
that can be used are

> - 0  - Free boundary condition
> - 1  - Used to enforce $v_x = 0$.
> - 2  - Used to enforce $v_y = 0$.
> - 3  - Used to enforce $v_z = 0$.
> - 4  - Used to enforce $v_r = 0$, or in other words 0 radial movement.
> - 5  - Used to enforce arbitrary bcs, to be handled in the problem.h file. [5+]
> - 99 - Off limits as this is used in BdrVertexIndexingArray to indicate corner vertices that should not move at all.

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

## Examples

### Elastic shocktube
```
./Laglos -m ../data/elastic/ref-segment.mesh -p 50 -tf 0.00005 -cfl 0.5 -ue 1 -rs 8
```
### Elastic impact
```
./Laglos -m ../data/elastic/ref-segment.mesh -p 51 -tf 0.00005 -cfl 0.5 -ue 1 -rs 8
```
### Elastic shear
```
./Laglos -m ../data/elastic/tube-100x1y.mesh -p 52 -tf 0.00005 -cfl 0.5 -ue 2 -ppd -rs 0
./Laglos -m ../data/elastic/distube-100x1y.mesh -p 52 -tf 0.00005 -cfl 0.5 -ue 2 -ppd -rs 0
./Laglos -m ../data/elastic/tube-2x100y.mesh -p 52 -tf 0.00005 -cfl 0.5 -ue 2 -ppd -rs 0
```
### Elastic shear y direction (rotated x)
```
./Laglos -m ../data/elastic/tube-1x100y.mesh -p 55 -tf 0.00005 -cfl 0.5 -ue 2 -ppd -rs 0
```
### Elastic impact + shear
```
./Laglos -m ../data/elastic/tube-100x1y.mesh -p 56 -tf 0.00005 -cfl 0.5 -ue 2 -ppd -rs 0
```
### Elastic 2D Rotation
```
./Laglos -m ../data/elastic/ref-square-c0-p15.mesh -p 57 -tf 0.00001 -cfl 0.5 -vis -vs 1 -ue 1 -rs 6 -ppd
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
./Laglos -m ../data/elastic/projectile-plate.mesh -p 54 -tf .000103 -cfl 0.5 -ue 1 -ppd -rs 1
```

### Elastic Noh
```
./Laglos -m ../data/elastic/ref-square-c0.mesh -p 58 -tf 0.000002 -cfl 0.5 -ue 1 -ppd -rs 6
./Laglos -m ../data/elastic/noh-nonuniform.mesh -p 58 -tf 0.000002 -cfl 0.5 -ue 1 -ppd -rs 1
```

### Elastic twist
```
./Laglos -m ../data/elastic/ref-square-c0.mesh -p 57 -tf 0.00005 -cfl 0.5 -ue 1 -ppd -rs 6
```

### Elastic projectile impact
This is a test case outlined in vilar-main-shu-2d. Currently
not yielding great results, perhaps due to the mesh velocity 
computation we employ.  Have tried it on a cartesian and
distorted mesh. 
```
./Laglos -m ../data/elastic/test-distorted-nonsymmetric.mesh -p 59 -tf 0.005 -cfl 0.5 -ue 1 -rs 2 -vis -vs 1

./Laglos -m ../data/elastic/test-distorted-nonsymmetric.mesh -p 59 -tf 0.005 -cfl 0.5 -ue 1 -rs 2 -vis -vs 1
```