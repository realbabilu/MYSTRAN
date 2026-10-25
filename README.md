#MYSTRAN 18.a
============
Active Branch is https://github.com/realbabilu/MYSTRAN/tree/v18.00.a

MYSTRAN is an acronym for “My Structural Analysis” (https://www.mystran.com)
This is a beta repository using added some features to test
---
This Mystran 18.a is fork of MYSTRANSOLVER https://github.com/MystranSolver/MYSTRANSolver from version 18 edition.
The different is this MYSTRAN are using some decks that may not existed in Nastran, but useful in the Civil Engineering Purpose.
For Example: PBEAMZ, PBeamZ is PBeamL but using more variables like auto-Tapered Cubic/Linear/Parabolic, and
Stiffness Modifier that made Area Stiffness different, but real area can be used for grav selfweight density,
More easier dimension parameter for section, etc.

[Build Instructions](#Build-Instructions) |
[Introduction](#Introduction) |
[Features](#Features) |
[Get EXE or Make Binary](#Get-EXE-or-Make-Binary) |
[Documentation](#Documentation) |
[Four Repositories](#Four-Repositories) |
[Developmental Goals](#Developmental-Goals) |
[Ways You Can Help](#ways-you-can-help) |
[Community](#community)

---
## October 25 update: linear shell elements and benchmark suite

The Windows packages now contain the executable used for the October 25 linear-shell benchmark: source `21540039` plus the supplied element-local patch, published with benchmark commit `86a28855`. The exact source delta, executable hash and numerical limits are included in the package; this is a tested snapshot with the candidate status stated below.

- CTRIA3 MITC3+ uses Lee + consistent-spin drilling. CTRIAR KMITC3+ is the Nguyen-based successor to retired MITC3+HB; their benchmark results remain separate.
- DSQK uses covariant shear and consistent spin. DSG3 includes consistent spin and the corrected transverse local-z geometric-stiffness block. DKMQ24 elastic, thermal and buckling work is represented by the tested source snapshot.
- DKMT18 uses the SF shear adaptation in this executable. **SF remains a candidate:** pipe-stress trade-offs, inherited rigid-rotation stabilization and shear-patch/output gates remain open; the edge-covariant alternative is retained in research history.
- The linear publication covers 22 selectors across CQUAD4, CQUADR, CTRIA3 and CTRIAR: 1,694 MYSTRAN configurations, 308 corresponding Nastran configurations and eight additional rank decks. It includes 5,096 metric rows, 17 reports, 99 figures and the tested DAT archive. NATIVE means Nastran for the corresponding element card. Quadratic reports remain separate.
- 1,693 MYSTRAN configurations completed; corrected DSG3 2-016c BEND still fails six-mode residual acceptance. All 308 Nastran configurations completed. Completion is not accuracy acceptance. Buckling decks explicitly request ARP_TOL=1e-10; default DARPACK and the 1e-6 physical residual limit are unchanged. No shared iterative-refinement/Jacobi prototype is included.
- Added 2-002 two-axis torsion and 2-007 hemisphere refinement series. Hemisphere GRID-CD results are transformed to basic coordinates; torsion uses the two physical tip corners. DKMT18 SF reaches a 0.99610044 torsion/reference ratio on 192x16.

[Linear comparisons](https://github.com/realbabilu/MYSTRAN/blob/v18.00.a/update/shell_benchmark/comparison_all_linear.md) · [CSV](https://github.com/realbabilu/MYSTRAN/blob/v18.00.a/update/shell_benchmark/comparison_all_linear.zip) · [Publication QA](https://github.com/realbabilu/MYSTRAN/blob/v18.00.a/update/shell_benchmark/PUBLICATION_QA_linear.md)

Extract the full Windows package, keeping its matching DLL beside the executable. Standalone users must update both mystran.exe.zip and libopenblas.dll.zip. Earlier dated package validation below is historical evidence and does not assert a fresh solid-element regression for this shell snapshot.

## CHEXA8 TIED update — 18 October 2026

- NEWSOLID CHEXA8 defaults to `PARAM,HEXATYP,TIED`, using mean-gradient-corrected physical shear tying with EAS9 condensation. `PARAM,HEXATYP,EAS9` retains EAS9v2.
- The [updated graph](update/solid_benchmark/comparison_all.md#chexa8) compares TIED and EAS9v2 on the existing beam/dome benchmark meshes against the previously published Nastran results.
- On the 24x24x1 dome, point A error is **0.759% for TIED**, **7.333% for EAS9v2**, and **0.682% for Nastran**. Default TIED uses ZETA; `HEXATAX` also provides explicit XI/ETA and opt-in AUTO controls.

[Windows binary and runtime](update/binaries/mystran18a-chexa8-tied-win64.zip) · [CHEXA8 formulation](update/Formulation/Solid/NEWSOLID_CONTINUUM_FORMULATIONS.md#35-chexa8-tied-default-and-retained-eas9-option) · [Benchmark comparison](update/solid_benchmark/comparison_all.md#chexa8)

## Solid recovery update - 18 October 2026

The [5-001 recovery correction](update/solid_element_recovery_5001_2026-10-18.md) restores ordinary OP2 element stress and strain for all 462 elements in the eight-solid patch. CHEXA8 corner recovery uses valid interior states; CPYRAM5/13 now supply center and corner tensors. The original patch displacements and GPSTRESS VOLUME payloads are unchanged. See the [recovery audit](update/solid_element_recovery_5001_2026-10-18.md) and audit for validation scope.

## Solid update - 17 October 2026

The [October release](https://github.com/realbabilu/MYSTRAN/releases/tag/mystran18a-october2026) includes the latest MYSTRAN 18.0.A Windows binary from `v18.00.a`.

- `NEWSOLID` is now the default; select `PARAM,SOLIDTYP,LEGACY` for historical solid paths.
- Enhanced CHEXA8 EAS9, CPENTA6 ABE/EAS9, CTETRA4 smoothing and Bedrosian CPYRAM13 are integrated. CPENTA6 defaults to ABE; `PARAM,PENTATYP,EAS9` selects EAS9. `EAS6` is not supported. CPYRAM accepts 5 or 13 nodes.
- CPENTA15 and CHEXA20 selective shear integration is corrected. Both pass the 5-001 six constant-strain patches and rigid-body rotation audit below 1e-7 relative error, with no extra global DOFs.
- Solid GPSTRESS VOLUME writes OP2 direct tensors, principal stresses and directions, pressure and von Mises for all eight 5-001 solid groups. Corner-grid rows and metadata match the tested Nastran reference. In the tested 5-001 case, seven volumes produce OP2 float payloads identical to Nastran. CPENTA6 EAS9 MYSTRAN is closer to the analytical solution than the Nastran reference. Thus, the tested GPSTRESS results are equal to or better than that reference for this case; this comparison does not extend to all solid benchmarks.
- Solid OP2 recovery indexing, mixed linear/quadratic records and engineering-shear principal strains are corrected. CTETRA4 NEWSOLID geometric stiffness matches the audited Python buckling operator and is independent of output coordinates.

Validation remains scoped: ABE on the distorted 5-001 patch and some CPENTA6 buckling benchmarks remain outside the established analytical acceptance tolerance. The ordinary CHEXA8 recovery and CPYRAM13 OES/OSTR limitations recorded on 17 October are superseded by the 18 October recovery correction above.

See the [5-001 OP2 audit](update/op2_solid_5001_audit.md), [quadratic patch audit](update/quadratic_patch_5001_audit.md), [solid defaults and limitations](update/newsolid_op2_default_2026-10-17.md), and [CTETRA4 buckling audit](update/ctetra4_python_kg_parity_2026-10-17.md).

### Solid benchmark comparisons

- [GPSTRESS VOLUME and solid-element patch-test example](update/solid_benchmark/5-001_8_SOLIDS_GPSTRESS_PENTA_REMESH_FIXED_QUADRATIC_MYSTRAN_EAS9.dat): validated 5-001 MYSTRAN input with eight solid groups, reconciled quadratic node ordering and CPENTA6 EAS9. See the [example notes](update/solid_benchmark/README.md#5-001-gpstress-volume-and-solid-element-patch-test-example) for validation scope.
- [MYSTRAN solid comparison vs Nastran](update/solid_benchmark/mystran_solid_comparison.md): patch, beam, dome, buckling and thermal benchmarks.
- [All solid elements - SAP2000 benchmarks](update/solid_benchmark/comparison_all.md): eight solid types, matched Nastran/MYSTRAN meshes and refinement curves, with independent and published SAP2000 references.
- [Download benchmark input decks only](update/solid_benchmark/solid_benchmark_inputs.zip) - no raw solver outputs.
- [Benchmark index and numerical data](update/solid_benchmark/README.md).

## Quadratic CQUAD8 and CTRIA6 Shell Formulations

Research formulations for quadratic shell elements with 6 DOF per node: five **CQUAD8** and four **CTRIA6** variants. The theory pages and CSI benchmark results distinguish independent-reference validation from pending investigations.

- [Quadratic shell formulation overview](update/Formulation/Shell/CQUAD8_CTRIA6_QUADRATIC_SHELL_FORMULATIONS.md)
- [Individual Q8/T6 formulations](update/Formulation/Shell/README.md)
- [Shell benchmark overview](update/shell_benchmark/README.md)
- [All-case reference-based comparison](update/shell_benchmark/comparison_all.md)
- [Browser comparison](update/shell_benchmark/comparison_all.md)
- [Benchmark DAT archive](update/shell_benchmark/DAT/quadratic_shell_benchmark_DAT.zip)

Note: CSV/Python files are research-internal and not part of the public repository.


## Formulation

Mathematical formulations, nodal numbering, integration rules, figures and existing benchmark comparisons are available in the reports below.

### Solid elements

[Solid formulation index](update/Formulation/Solid/README.md) · [Continuum overview: all eight solid topologies](update/Formulation/Solid/NEWSOLID_CONTINUUM_FORMULATIONS.md)

The continuum overview also covers CTETRA10, CPENTA15, CHEXA20 and CPYRAM5, including the applicable selective-integration choices. Detailed reports for the special formulations are:

- [TETRA4 NEWSOLID — nodal strain smoothing](update/Formulation/Solid/TETRA4_Formulation.md)
- [CPYRAM13 Bedrosian — rational quadratic pyramid](update/Formulation/Solid/CPYRAM13_Bedrosian_Formulation.md)
- [CHEXA8 EAS — covariant enhancement and consistent shear tying](update/Formulation/Solid/CHEXA8_EAS_Formulation.md)
- [CPENTA6 ABE — assumed shear and condensed enhanced strains](update/Formulation/Solid/CPENTA6_ABE_Formulation.md)

### Shell elements

[Shell formulation directory](update/Formulation/Shell) · [Quadratic Q8/T6 overview](update/Formulation/Shell/CQUAD8_CTRIA6_QUADRATIC_SHELL_FORMULATIONS.md)

#### Four-node quadrilateral shells

- [DSQK](update/Formulation/Shell/DSQK_formulation.md)
- [DKMQ20](update/Formulation/Shell/DKMQ20_Formulation.md)
- [DKMQ24](update/Formulation/Shell/DKMQ24_formulation.md)
- [DKMQ24 EAS](update/Formulation/Shell/DKMQ24EAS_Formulation_and_Benchmarks.md)
- [DKMQ24 AU](update/Formulation/Shell/DKMQ24AU_Formulation.md)
- [Simo 1989](update/Formulation/Shell/Simo1989_Formulation.md)
- [Simo 1993 Q4](update/Formulation/Shell/Simo1993_Q4_Formulation.md)
- [Q4EASANS](update/Formulation/Shell/Q4EASANS_Formulation.md)
- [MBP1C0](update/Formulation/Shell/MBP1C0_Formulation.md)
- [Q4RS](update/Formulation/Shell/Q4RS_Formulation.md)

#### Three-node triangular shells

- [DKMT18](update/Formulation/Shell/DKMT18_Formulation.md)
- [DSG3](update/Formulation/Shell/DSG3_Formulation.md)
- [MITC3+ Lee](update/Formulation/Shell/MITC3Plus_Lee_Formulation.md)
- [kMITC3+ Nguyen 2024](update/Formulation/Shell/KMITC3Plus_Nguyen2024_Formulation.md)
- [T3FF / T3FFD](update/Formulation/Shell/T3FF_Formulation.md)

#### Quadratic quadrilateral shells — Q8

- [SIMOQ8](update/Formulation/Shell/SIMOQ8_v13.md)
- [MITC8 / MITC8D](update/Formulation/Shell/MITC8_v4.md)
- [ANS8BDG6](update/Formulation/Shell/ANS8BDG6_v4.md)
- [HBQ8](update/Formulation/Shell/HBQ8_v3.md)
- [MacNeal Q8 with drilling](update/Formulation/Shell/MACQ8D_v3_drill.md)

#### Quadratic triangular shells — T6

- [SIMOT6](update/Formulation/Shell/SIMOT6_v2.md)
- [MITC6](update/Formulation/Shell/MITC6_v5e.md)
- [MacNeal MH6T](update/Formulation/Shell/MH6T_v4.md)
- [Rezaiee 2017](update/Formulation/Shell/REZAIEE_v3.md)

# Build Instructions

See [BUILD.md](BUILD.md) for both Windows and Linux build (compiling) instructions.

Since it added with several features, a minimum need OpenBlas library, use -DTPL_BLAS_LIBRARIES: to pointing where BLAS is, or it will automatically find openblas.
create build folder and run it from there
mkdir build
cd build

Minimum build :

cmake -G "MinGW Makefiles"  -DCMAKE_BUILD_TYPE=Release  -DCMAKE_C_COMPILER=gcc.exe  -DCMAKE_CXX_COMPILER=g++.exe  -DCMAKE_Fortran_COMPILER=gfortran.exe .. -DTPL_BLAS_LIBRARIES="..\libopenblas.dll.a"  -DMYSTRAN_USE_EXTERNAL_FEAST=OFF  -DMYSTRAN_USE_DMUMPS_SOLVER=OFF

Maximum build : download feast and mumps in each folder.

cmake -G "MinGW Makefiles"  -DCMAKE_BUILD_TYPE=Release  -DCMAKE_C_COMPILER=gcc.exe  -DCMAKE_CXX_COMPILER=g++.exe  -DCMAKE_Fortran_COMPILER=gfortran.exe .. -DTPL_BLAS_LIBRARIES="..\libopenblas.dll.a"  -DMYSTRAN_USE_DMUMPS_SOLVER=ON -DMYSTRAN_USE_EXTERNAL_FEAST=ON

Maximum build with already built FEAST and MUMPS-non MPI Libraries

cmake -G "MinGW Makefiles"  -DCMAKE_BUILD_TYPE=Release  -DCMAKE_C_COMPILER=gcc.exe  -DCMAKE_CXX_COMPILER=g++.exe  -DCMAKE_Fortran_COMPILER=gfortran.exe .. -DTPL_BLAS_LIBRARIES="..\libopenblas.dll.a"  -DMYSTRAN_USE_DMUMPS_SOLVER=ON -DMYSTRAN_USE_EXTERNAL_FEAST=ON -DMYSTRAN_EXTERNAL_SUPERLU_LIB="../superlu/src/libsuperlu.a" -DMYSTRAN_EXTERNAL_SUPERLU_INCLUDE_DIR="../superlu/src" -DMYSTRAN_EXTERNAL_SUPERLU_CONFIG_DIR="../superlu/src" -DMYSTRAN_EXTERNAL_SUPERLU_DRIVER="../superlu/FORTRAN/c_fortran_dgssv.c" -DMYSTRAN_USE_EXTERNAL_FEAST=ON  -DMYSTRAN_FEAST_EXTRA_LIBS="C../feast/4.0/lib/x64/libfeast.a" -DMYSTRAN_USE_DMUMPS_SOLVER=ON -DMYSTRAN_DMUMPS_INCLUDE_DIR="../mumps/include"   -DMYSTRAN_DMUMPS_EXTRA_LIBS="../mumps/lib/libdmumps.a;../mumps/lib/libmpiseq.a;../mumps/lib/libmumps_common.a;../mumps/lib/libpord.a;../mumps/lib/libsmumps.a"

Use mumps with no-mpi, but using multi-thread optimized like OPENBLAS

For Windows intel OneAPI use MKL LP64 and -DMYSTRAN_DISABLE_NDEBUG=ON for Release version, and edit slu_cnames.h SuperLU to UPCASE only.

# Introduction

MYSTRAN is a general purpose finite element analysis computer program for
structures that can be modeled as linear (i.e. displacements, forces and
stresses proportional to applied load). MYSTRAN is an acronym for
“My Structural Analysis”, to indicate its usefulness in solving a wide variety
of finite element analysis problems.

For anyone familiar with the popular NASTRAN computer program developed by NASA
(National Aeronautics and Space Administration) in the 1970’s and popularized
in several commercial versions since, the input to MYSTRAN will look quite
familiar. Many structural analyses modeled for execution in NASTRAN will
execute in MYSTRAN with little, or no, modification. MYSTRAN, however, is not
NASTRAN. It is an independent program written in modern Fortran 95.

# Features

- Nastran compatibility
- Linear Static Analysis
- Modal analysis
- Linear Elastic Buckling Analysis
- Full Suite of 1D, 2D, and 3D elements
- Selectable shell formulations through PARAM,QUAD8TYP ; PARAM,TRIA6TYP ; PARAM,QUAD4TYP ; PARAM,TRIA3TYP ; PARAM,QUADRTYP ; PARAM,TRIARTYP.
- CQUAD4 shell formulations: MIN4, MIN4T, MITC4, MITC4+, DSQK, DKMQ20, and Simo1989 with Hughes–Brezzi stabilization.
- CQUADR shell formulations: Simo1993, DKMQ24, DKMQ24 with EAS membrane enhancement, DKMQ24 AU, Q4EASANS, Müller–Bischoff P1C0, Krysl Q4RS, and MITC4+/D with Hughes–Brezzi stabilization.
- CTRIA3 shell formulations: MIN3, MITC3+, Krysl T3FF, and DSG3.
- CTRIAR shell formulations: DKMT18, Krysl T3FFD, and MITC3+ with Hughes–Brezzi stabilization.
- CQUAD8 shell formulations: SIMOQ8, MITC8, MITC8D, ANS8BDG6, HBQ8, and MacNeal Q8.
- CTRIA6 shell formulations: SIMOT6, MITC6, MacNeal MH6T, and Rezaiee2017
- Enhanced assumed-strain shell formulations: automatic shear blending and membrane tying for selected Q8 elements, with configurable ANS field, shear, and membrane options.
- Extended GPSTRESS nodal averaging for quadratic/linear shells, including midside nodes.
- CBEAM with PBEAML Nastran
- Faster Solver MUMPS for alternative SUPERLU
- Some Eigen Solver: FEAST, Subspace, DSYEV
- Faster RCM for Banded Optimization
- New Solid with EAS linear and quadratic including CPYRAM, CHEXA, CTETRA, CPENTA
- Consistent linear and quadratic-shell pressure loads: surface-normal and specified-direction loading, including midside-node contributions and curved-surface integration.
- Structural thermal loading for shells: uniform temperature changes and through-thickness temperature gradients.
- Consistent and HRZ lumped shell mass, including nonstructural mass; rotary inertia in selected upgraded formulations.
- Geometric stiffness and linear buckling support for shells, validated on flat plate and column benchmarks.
- Expanded shell recovery: membrane forces, bending moments, transverse shear, and top/bottom fiber stresses at center and nodal locations.
- Native linear and quadratic-shell OP2 stress and force output, with improved SURFACE/GPSTRESS handling and F06 consistency.
- Shell FORCE output selectors: CENTER, CORNER, and combined locations for native Q8/T6 shells

# Get EXE or Make Binary

Windows EXE (executable) for can be found in the "Releases" section of this page (right hand pane).

Static Linux binaries have been built, but releases are in work.
For now, it is better to build it yourself -- it's really
straightforward.

# Shell and Solid Elements
| Elements card | PARAM selector | Values |
|---|---|---|
| CQUAD4 | `QUAD4TYP` | `MIN4`, `MIN4T`, `MITC4`, `MITC4+`, `DSQK`, `DKMQ20`, `SIMO` |
| CQUADR | `QUADRTYP` | `SIMO`, `DKMQ24`, `DKM24EA`, `DKM24AU`, `Q4EASANS`, `MBP1C0`, `Q4RS`, `MITC4PD` |
| CTRIA3 | `TRIA3TYP` | `MIN3`, `MITC3+`, `T3FF`, `DSG3` |
| CTRIAR | `TRIARTYP` | `DKMT18`, `T3FFD`, `KMITC3+` (experimental) |
| CQUAD8 | `QUAD8TYP` | `SIMOQ8`, `MITC8`, `MITC8D`, `ANS8BDG6`, `HBQ8`, `MACQ8D` |
| CTRIA6 | `TRIA6TYP` | `SIMOT6`, `MITC6`, `MH6T`, `REZAIEE` |
| SOLID  | `SOLIDTYP` | `NEWSOLID` (default), `LEGACY` |
| CHEXA8 NEWSOLID | `HEXATYP` | `TIED` (default), `EAS9` |
| CHEXA8 TIED axis | `HEXATAX` | `ZETA` (default), `XI`, `ETA`, `AUTO` (opt-in) |
| CPENTA6 NEWSOLID | `PENTATYP` | `ABE` (default), `EAS9` |

# Documentation

The end user documentation is located the [MYSTRAN_Documentation](https://github.com/MYSTRANsolver/MYSTRAN_Documentation) repository.
This includes a Quick Setup Guide, User Manual, and Theory Manual.

# Four Repositories

The MYSTRAN project consists of five repositories.

1 - This repository contains the source code and build instructions.

2 - The [MYSTRAN_Documentation](https://github.com/MYSTRANsolver/MYSTRAN_Documentation) repository contains various documents related to the MYSTRAN program.

3 - The [MYSTRAN_Resources](https://github.com/MYSTRANsolver/MYSTRAN_Resources) repository consists of files for MYSTRAN developers.
It also contains information and files related to pre- and post-processors relevant to MYSTRAN.

4 - The [MYSTRAN_Benchmark](https://github.com/MYSTRANsolver/MYSTRAN_Benchmark) repository contains the test cases and utilities used to verify that a new build produces results consistent with prior builds and models that have been verified.


# Developmental Goals

- Continue the implementation of the MITC shell elements and shell element buckling capability
- Validation effort (hundreds/thousands of test cases)
- Discover and resolve bugs
- Improve performance

# Ways You Can Help

- Join the MYSTRAN Discord Channel and/or Forum (links below)
- Report bugs and inconsistencies
- Report issues with Documentation
- A large validation  effort is underway. This will require the assistance of the community. Any help would be greatly appreciated.

# Community

- [Join our Discord Channel](https://discord.gg/9k76SkHpHM) - Very active.
