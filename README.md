MYSTRAN 18.a
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
