# Laglos - Lagrangian Low-Order Solver

Laglos (LAGrangian Low-Order Solver) is a high-performance computational fluid dynamics miniapp for solving the time-dependent Euler equations of compressible gas dynamics in a moving Lagrangian frame. Built on the MFEM finite element library, Laglos employs unstructured low-order finite element spatial discretization with forward Euler time-stepping to simulate complex hydrodynamic phenomena.

## Table of Contents
- [Key Features](#key-features)
- [Applications](#applications)
- [Building Laglos](#building-laglos)
  - [Dependencies](#dependencies)
  - [Build Instructions](#build-instructions)
  - [Building Laghos Extension](#building-laghos-extension)
- [Quick Start](#quick-start)
  - [Hydrodynamics Examples](#hydrodynamics-examples)
  - [Elasticity Examples](#elasticity-examples)
- [Runtime Parameters](#runtime-parameters)
- [Verification and Validation](#verification-and-validation)
  - [Mass Conservation](#mass-conservation)
  - [Convergence Analysis](#convergence-analysis)
- [Advanced Topics](#advanced-topics)
  - [Boundary Conditions](#boundary-conditions)
  - [Multi-Material Meshes](#multi-material-meshes)
  - [Elastic Problems](#elastic-problems)
- [Contributing](#contributing)
- [References](#references)

## Key Features

- **Lagrangian Framework**: Solves hydrodynamic equations on a moving mesh that follows material flow
- **Low-Order Finite Elements**: Uses robust DG0/Q1 discretizations for improved stability
- **Flexible Problem Suite**: Includes classic test problems (Sod shock tube, Sedov blast wave, Noh implosion, etc.)
- **Advanced Physics**: Supports various equations of state including Van der Waals gas models
- **Multi-Material Support**: Handle elastic-plastic materials and multi-phase flows
- **Viscosity Options**: Multiple artificial viscosity formulations for shock capturing
- **Mesh Velocity Algorithms**: Configurable mesh motion schemes for optimal Lagrangian evolution
- **HiOp Integration**: Advanced optimization capabilities for mesh movement

## Applications

Laglos is designed for:
- Compressible hydrodynamics research
- Shock physics simulations
- Elastic-plastic material modeling
- Inertial confinement fusion modeling
- Algorithm development for Lagrangian methods

---

## Building Laglos

### Dependencies

Laglos requires the following external libraries:

| Library | Version | Purpose | Required |
|---------|---------|---------|----------|
| [MFEM](https://github.com/mfem/mfem) | master branch | Finite element discretization | Yes |
| [HiOp](https://github.com/LLNL/hiop) | master branch | Nonlinear optimization solver | Optional* |
| [CoinHSL](http://hsl.rl.ac.uk/ipopt) | latest | Linear algebra for optimization | Optional* |

**Optional:**
- [GLVis](https://github.com/GLVis/glvis) - Real-time visualization

> **Note on Optional Dependencies**: CoinHSL and HiOp are only required for certain advanced mesh velocity calculations involving optimization. These are research-specific dependencies and are **not necessary for out-of-the-box use cases** of Laglos. If you're not using optimization-based mesh velocity algorithms (e.g., `-mv` options that require HiOp), you can skip building these libraries and only build MFEM.
>
> **MFEM branch for optimization builds**: If you choose to use HiOp/CoinHSL, you must build MFEM from the `hiop-sparse` branch:
> `git checkout hiop-sparse`. Other branches may not include the required sparse/HiOp interfaces.

### Build Instructions

#### 1. Download and Setup Dependencies

Create a workspace directory and organize libraries:

```sh
~> mkdir Workspace && cd Workspace
~/Workspace> ls
Laglos/  coinhsl-x.y.z.tar.gz  hiop/  mfem/
```

#### 2. Build CoinHSL (Optional - for optimization)

1. Register at http://hsl.rl.ac.uk/ipopt to obtain download link
2. Download and extract:

```sh
~/Workspace> gunzip coinhsl-x.y.z.tar.gz
~/Workspace> tar xf coinhsl-x.y.z.tar
~/Workspace> ln -s coinhsl-x.y.z coinhsl
~/Workspace> cd coinhsl
~/Workspace/coinhsl> meson setup build \
    -Dcpp_args="-fPIC" \
    -Dc_args="-fPIC" \
    -Dfortran_args="-fPIC" \
    --buildtype=release \
    --prefix=$(pwd)/install
~/Workspace/coinhsl> meson compile -C build
~/Workspace/coinhsl> meson install -C build
```

#### 3. Build HiOp (Optional - for optimization)

```sh
~/Workspace> cd hiop
~/Workspace/hiop> mkdir build && cd build
~/Workspace/hiop/build> cmake \
    -DCMAKE_INSTALL_PREFIX=../install \
    -DHIOP_COINHSL_DIR=../../coinhsl/install \
    -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
    -S ../
~/Workspace/hiop/build> make -j8
~/Workspace/hiop/build> make test
~/Workspace/hiop/build> make install
```

#### 4. Build MFEM

```sh
~/Workspace> git clone https://github.com/mfem/mfem.git ./mfem
~/Workspace> cd mfem/
~/Workspace/mfem> git checkout master
~/Workspace/mfem> mkdir build && cd build
~/Workspace/mfem/build> cmake \
    -DCMAKE_INSTALL_PREFIX=$(pwd)/install \
    -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
    -S ../
~/Workspace/mfem/build> make -j8
~/Workspace/mfem/build> make install
```

**Note on `-DCMAKE_POSITION_INDEPENDENT_CODE=ON`:**  
This flag is critical on Linux systems for proper shared library linking. It ensures static libraries are compiled with position-independent code (PIC), allowing them to be linked into shared objects. On Linux, static libraries (.a) do not include PIC unless explicitly compiled with it.

#### 5. Build Laglos

```sh
~/Workspace> cd Laglos
~/Workspace/Laglos> mkdir build && cd build
~/Workspace/Laglos/build> cmake \
    -DCMAKE_INSTALL_PREFIX=$(pwd)/install \
    -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
    ../
~/Workspace/Laglos/build> make -j8
~/Workspace/Laglos/build> make test
~/Workspace/Laglos/build> make install
```

#### 6. Build GLVis (Optional)

For visualization:

```sh
~/Workspace> git clone https://github.com/GLVis/glvis.git ./glvis
~/Workspace> cd glvis
~/Workspace/glvis> make
```

Run GLVis in a separate terminal, then use `-vis` flag in Laglos to stream results.

### Building Laghos Extension

The Laghos miniapp can be built against Laglos:

```sh
~/Workspace> cd Laghos
~/Workspace/Laghos> mkdir build && cd build
~/Workspace/Laghos/build> cmake -DLaglos_DIR=../../Laglos/build -S ..
~/Workspace/Laghos/build> make -j8
```

---

## Quick Start

### Hydrodynamics Examples

#### Sod Shock Tube (1D)
```sh
~/Laglos/build> ./Laglos -m ../data/ref-segment.mesh -p 2 -tf 0.225 -cfl 0.5 -rs 8 -vis
```

#### Sedov Blast Wave (2D)
```sh
~/Laglos/build> ./Laglos -m ../data/ref-square-N15.mesh -p 1 -tf 0.9 -cfl 1 -rs 5 -vis
```

#### Triple Point Problem
```sh
~/Laglos/build> ./Laglos -m ../data/ref-square.mesh -p 3 -tf 5.0 -cfl 0.5 -rs 6 -vis
```

### Elasticity Examples

Enable elastic flux with `-ue <flag>` option.

#### Elastic Shock Tube
```sh
~/Laglos/build> ./Laglos -m ../data/elastic/ref-segment.mesh -p 50 -tf 0.00005 -cfl 0.5 -ue 1 -rs 8
```

#### Elastic Impact
```sh
~/Laglos/build> ./Laglos -m ../data/elastic/ref-segment.mesh -p 51 -tf 0.00005 -cfl 0.5 -ue 1 -rs 8
```

#### Elastic Shear
```sh
~/Laglos/build> ./Laglos -m ../data/elastic/tube-100x1y.mesh -p 52 -tf 0.00005 -cfl 0.5 -ue 2 -ppd -rs 0
~/Laglos/build> ./Laglos -m ../data/elastic/distube-100x1y.mesh -p 52 -tf 0.00005 -cfl 0.5 -ue 2 -ppd -rs 0
~/Laglos/build> ./Laglos -m ../data/elastic/tube-2x100y.mesh -p 52 -tf 0.00005 -cfl 0.5 -ue 2 -ppd -rs 0
```

#### Elastic Shear (Y Direction)
```sh
~/Laglos/build> ./Laglos -m ../data/elastic/tube-1x100y.mesh -p 55 -tf 0.00005 -cfl 0.5 -ue 2 -ppd -rs 0
```

#### Elastic Impact + Shear
```sh
~/Laglos/build> ./Laglos -m ../data/elastic/tube-100x1y.mesh -p 56 -tf 0.00005 -cfl 0.5 -ue 2 -ppd -rs 0
```

#### Elastic 2D Rotation
```sh
~/Laglos/build> ./Laglos -m ../data/elastic/ref-square-c0-p15.mesh -p 57 -tf 0.00001 -cfl 0.5 -vis -vs 1 -ue 1 -rs 6 -ppd
```

#### Elastic Projectile Plate

The final time depends on the shear modulus used (modified in the problem file). Final times from the referenced paper:

**For μ = 9.2×10¹⁰ Pa:**
```sh
~/Laglos/build> ./Laglos -m ../data/elastic/projectile-plate.mesh -p 54 -tf 0.000036 -cfl 0.5 -ue 1 -ppd -rs 1  # t = 3.6×10⁻⁵
~/Laglos/build> ./Laglos -m ../data/elastic/projectile-plate.mesh -p 54 -tf 0.000106 -cfl 0.5 -ue 1 -ppd -rs 1  # t = 1.06×10⁻⁴
~/Laglos/build> ./Laglos -m ../data/elastic/projectile-plate.mesh -p 54 -tf 0.00032 -cfl 0.5 -ue 1 -ppd -rs 1   # t = 3.2×10⁻⁴
~/Laglos/build> ./Laglos -m ../data/elastic/projectile-plate.mesh -p 54 -tf 0.000609 -cfl 0.5 -ue 1 -ppd -rs 1  # t = 6.09×10⁻⁴
```

**For μ = 1×10⁹ Pa:**
```sh
~/Laglos/build> ./Laglos -m ../data/elastic/projectile-plate.mesh -p 54 -tf 0.000035 -cfl 0.5 -ue 1 -ppd -rs 1  # t = 3.5×10⁻⁵
~/Laglos/build> ./Laglos -m ../data/elastic/projectile-plate.mesh -p 54 -tf 0.00014 -cfl 0.5 -ue 1 -ppd -rs 1   # t = 1.4×10⁻⁴
~/Laglos/build> ./Laglos -m ../data/elastic/projectile-plate.mesh -p 54 -tf 0.00042 -cfl 0.5 -ue 1 -ppd -rs 1   # t = 4.2×10⁻⁴
~/Laglos/build> ./Laglos -m ../data/elastic/projectile-plate.mesh -p 54 -tf 0.00071 -cfl 0.5 -ue 1 -ppd -rs 1   # t = 7.1×10⁻⁴
```

**For μ = 0 Pa:**
```sh
~/Laglos/build> ./Laglos -m ../data/elastic/projectile-plate.mesh -p 54 -tf 0.000075 -cfl 0.5 -ue 1 -ppd -rs 1  # t = 7.5×10⁻⁵
~/Laglos/build> ./Laglos -m ../data/elastic/projectile-plate.mesh -p 54 -tf 0.000187 -cfl 0.5 -ue 1 -ppd -rs 1  # t = 1.87×10⁻⁴
~/Laglos/build> ./Laglos -m ../data/elastic/projectile-plate.mesh -p 54 -tf 0.0006 -cfl 0.5 -ue 1 -ppd -rs 1    # t = 6×10⁻⁴
~/Laglos/build> ./Laglos -m ../data/elastic/projectile-plate.mesh -p 54 -tf 0.000103 -cfl 0.5 -ue 1 -ppd -rs 1  # t = 1.03×10⁻⁴
```

#### Elastic Noh
```sh
~/Laglos/build> ./Laglos -m ../data/elastic/ref-square-c0.mesh -p 58 -tf 0.000002 -cfl 0.5 -ue 1 -ppd -rs 6
~/Laglos/build> ./Laglos -m ../data/elastic/noh-nonuniform.mesh -p 58 -tf 0.000002 -cfl 0.5 -ue 1 -ppd -rs 1
```

#### Elastic Twist
```sh
~/Laglos/build> ./Laglos -m ../data/elastic/ref-square-c0.mesh -p 57 -tf 0.00005 -cfl 0.5 -ue 1 -ppd -rs 6
```

#### Elastic Projectile Impact
Test case from Vilar-Main-Shu 2D. Note: Currently yields mixed results on distorted meshes.
```sh
~/Laglos/build> ./Laglos -m ../data/elastic/test-distorted-nonsymmetric.mesh -p 59 -tf 0.005 -cfl 0.5 -ue 1 -rs 2 -vis -vs 1
```

---

## Runtime Parameters

> ⚠️ **Experimental Features**: Some options (marked with **[Experimental]**) are under active development. Results should be carefully validated.

| Parameter | Description |
|-----------|-------------|
| `-m <mesh>` | Mesh file to use |
| `-p <int>` | Problem number to run |
| `-rs <int>` | Number of serial mesh refinements |
| `-rp <int>` | Number of parallel refinements (not implemented) |
| `-tf <double>` | Final simulation time |
| `-ti <double>` | Initial time (not available for all problems) |
| `-cfl <double>` | Courant-Friedrichs-Lewy condition for time step |
| `-visc <int>` | Artificial viscosity type (0=none, 1-4=various formulations) |
| `-greedy` | Enable greedy viscosity (see GMPST2024) |
| `-ggn <int>` | Number of GMV steps before using greedy lambda |
| `-ue <int>` | Use elasticity (1=full elastic, 2=elastic shear only) |
| `-ppd` | Post-process density for mass-conservative update |
| `-vis` | Enable GLVis visualization |
| `-visit` | Enable VisIt output files |
| `-vs <int>` | Visualization snapshot frequency (every N steps) |
| `-of <path>` | Output folder for results |
| `-print` | Print state vectors to disk |
| `-gfprint` | Print grid functions |
| `-mv <int>` | Mesh velocity algorithm option |
| `-fv <int>` | Face velocity computation option |
| `-mv-iter-n <int>` | Number of mesh velocity solver iterations |
| `-do-mv-lin` | Enable mesh velocity linearization |
| `-tvc <double>` | Target velocity optimization viscosity coefficient |
| `-ms <int>` | Maximum number of time steps |
| `-cm` | Check mesh quality |
| `-mm` | Mesh motion flag |

For complete parameter list, run:
```sh
~/Laglos/build> ./Laglos --help
```

### Details

- ODE solvers (`-s`):
  - 1: Forward Euler **(DEFAULT)**
  - 2: RK2 SSP
  - 3: RK3 SSP
  - 4: RK4
  - 6: RK6
  - 7: RK2Avg
  - Note: IDP variants are also available internally (11,12,13,14,16), but these are only needed in the limiting Laghos case.

- Viscosity (`-visc`):
  - 0: None
  - 1: GMS-GV (Guaranteed Maximum Speed Graph Viscosity) **(DEFAULT)**
  - 2: Greedy viscosity **[Experimental]**
  - 3: Artificial graph viscosity (HO, Binder) **[Experimental]**
  - 4: Artificial graph viscosity (HO, Non-binder) **[Experimental]**

- Elasticity (`-ue`):
  - 0: No shear energy (fluid) **(DEFAULT)**
  - 1: Neo-Hookean
  - 2: Mooney-Rivlin
  - 3: Aortic
  - 4: Transversely isotropic
  - 5: Multiple shear EOS

- Mesh velocity (`-mv`):
  - 0: Arithmetic average of adjacent cells
  - 1: Arithmetic average with distributed viscosity
  - 2: Cell-face-normal with viscosity **(DEFAULT)**
  - 1*: Sparse HiOp LM with “*” target (requires HiOp/CoinHSL and MFEM branch hiop-sparse) **[Experimental]**

- Face velocity (`-fv`):
  - 0: Do nothing **(DEFAULT)**
  - 1: Mass-conservative bubble (Q2)
  - 2: Average (Q1-type)
  - 3: Butterfly (Q2) **[Experimental]**

---

## Verification and Validation

### Mass Conservation

Laglos monitors local and global mass conservation at each time step. All calculations related to mass error are performed in `LagrangianLOOperator::ValidateMassConservation()`. The relative mass error is computed as:

$$
\text{error}_{\text{mass}} = \frac{\sum_{c \in \eta^{\text{Cel}}} \left|\frac{|K_c^n|}{T_c^n} - m_c^{\rho}\right|}{\sum_{c \in \eta^{\text{Cel}}}m_c^{\rho}}
$$

where:
- $m_c^{\rho}$ = initial cell mass
- $|K_c^n|$ = cell volume at time $t^n$
- $T_c^n$ = specific volume at time $t^n$

A cell is flagged as violating mass conservation if:
$$
\left|\frac{|K_c^n|}{T_c^n} - m_c^{\rho}\right| > 10^{-12}
$$

The percentage of cells violating mass conservation is reported as the number of flagged cells divided by the total number of cells in the mesh.

### Convergence Analysis

For problems with known exact solutions, Laglos can compute convergence rates.

**Step 1:** Configure test parameters in `scripts/convergence_test_script.sh`

```bash
PROBLEM=40                           # Problem number
MESH="../data/ref-segment.mesh"      # Base mesh
FINAL_TIME=0.6                       # Simulation end time
CFL=0.5                              # CFL number
rs_levels=(4 5 6 7)                  # Refinement levels to test
OUTPUT_DIR="convergence_results"     # Output directory
```

**Step 2:** Run convergence suite

```sh
~/Laglos/scripts> ./convergence_test_script.sh
```

This will execute Laglos at each specified refinement level and save results.

**Step 3:** Generate convergence table

```sh
~/Laglos/scripts> python3 compute_convergence.py --results-dir ../build/convergence_results/
```

**Expected Output:**

**TODO:** Add an example convergence table here for the Sod case. Edit default convergence analysis script to correspond to this output.

The composite relative $L^1$ error is computed as:

$$
\delta^1(t) := \frac{\|\tau_h(t)-\tau(t)\|_{L^1(D)}}{\|\tau(t)\|_{L^1(D)}} + \frac{\|\mathbf{v}_h(t)-\mathbf{v}(t)\|_{L^1(D)}}{\|\mathbf{v}(t)\|_{L^1(D)}} + \frac{\|E_h(t)-E(t)\|_{L^1(D)}}{\|E(t)\|_{L^1(D)}}
$$

The convergence order between two refinement levels is computed by comparing error norms:
$$
\text{rate} = \frac{\log(e_1/e_2)}{\log(h_1/h_2)} \approx p
$$

where $e_i$ is the error at refinement level $i$, $h_i$ is the corresponding mesh size, and $p$ is the expected order of convergence.

**Expected convergence rates:**
- **DG0 (piecewise constant)**: Rate ≈ 1.0
- **Q1 (bilinear elements)**: Rate ≈ 2.0

For detailed instructions on other utility scripts, see [scripts/README.md](scripts/README.md).

---

## Advanced Topics

### Boundary Conditions

Laglos supports thermodynamic and mesh velocity boundary conditions via mesh attributes. Both types can be enabled in the problem file using flags:

```cpp
_thbcs = true;  // Enable thermodynamic BCs
_mvbcs = true;  // Enable mesh velocity BCs
```

#### Standard Boundary Attributes

| Attribute | Meaning |
|-----------|---------|
| 1 | Enforce $v_x = 0$ |
| 2 | Enforce $v_y = 0$ |
| 3 | Enforce $v_z = 0$ |
| 4 | Enforce $v_r = 0$ (no radial movement) |
| 5-8 | User-defined BCs (implement in problem file) |
| 9 | Free boundary (outflow) |
| 99 | Corner vertices (no movement) |

#### Custom Boundary Conditions

To implement non-standard BCs beyond $v \cdot n = 0$, override these functions in your problem class:

```cpp
void get_additional_BCs(const FiniteElementSpace &fes, 
                        Array<int> ess_bdr,
                        Array<int> &add_ess_tdofs,
                        Array<double> &add_bdr_vals);

void update_additional_BCs(const Vector &S,
                           Array<int> &add_ess_tdofs,
                           Array<double> &add_bdr_vals);
```

See `SaltzmannProblem` for an example with Dirichlet conditions on the left wall.

### Multi-Material Meshes

Laglos handles multi-material problems through element attributes in the mesh file:

- **Attribute < 50**: Fluid/gas elements (standard hydrodynamics)
- **Attribute ≥ 50**: Solid/elastic elements (uses elastic constitutive model)

When using `-ue` flag, only elements with attribute ≥ 50 will use elastic flux and stress calculations.

**Example multi-material mesh:**
```
# MFEM mesh file
...
elements
100
1 3 0 1 4 3     # Fluid region (attribute 1)
1 3 1 2 5 4     # Fluid region (attribute 1)
50 3 6 7 10 9   # Solid region (attribute 50)
50 3 7 8 11 10  # Solid region (attribute 50)
...
```

### Elastic Problems

#### Computing Errors with Tabulated Data

Since exact solutions for elastic problems are often given as tabulated approximation data, use the `compute_error.py` script:

```sh
~/Laglos/scripts> python3 compute_error.py \
    results/elastic/st1/state_vectors/ \
    ../exact_sol/elastic/shocktube/all_results.txt \
    7 16
```

The third and fourth arguments correspond to matching columns of the respective files being compared.

**Column mappings:**

| Variable | Laglos Column | Exact Solution Column |
|:---------|:-------------:|:---------------------:|
| ρ (density) | 1 | 13 |
| E (energy) | 3 | 15 |
| σ (stress) | 7 | 16 |

#### HiOp Solve Status

When using optimization for mesh movement, HiOp returns:

| Code | Meaning |
|------|---------|
| 0 | Solve successful |
| 1 | Solve successful (RelTol) |
| 2 | Acceptable solution |
| 5 | Infeasible problem |
| 10 | Max iterations exceeded |
| -16 | Error in feasibility restoration |

For a complete list of `hiopSolveStatus` values, see `hiop/src/Interface/hiopInterface.hpp`.

---

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Submit a pull request with clear description of changes

---

## References

For the underlying algorithms and theory, see:

- **Lagrangian Methods**: 
   - Dobrev et al., "High-Order Curvilinear Finite Element Methods for Lagrangian Hydrodynamics," SIAM J. Sci. Comput., 2012
   - Guermond et al., "Invariant-Domain Preserving and Locally Mass Conservative Approximation of the Lagrangian Hydrodynamics Equations," CMAME, 2025
- **Graph Viscosity**: 
   - Guermond et al., "Invariant Domains and First-Order Continuous Finite Element Approximation for Hyperbolic Systems, SIAM J. Numer. Anal. 54, 2016
- **Elastic-Plastic Models**: 
   - Favrie et al., "A Thermodynamically Compatible Splitting Procedure in Hyperelasticity," JCP, 2014.
   - Chaimoon et al., "An Anisotropic Hyperelastic Model with an Application to Soft Tissues," European Journal of Mechanics-A/Solids, 2019

---

## Contact

Madison Sheridan \
madison.sheridan94@gmail.com