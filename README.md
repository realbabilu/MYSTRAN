MYSTRAN
=======

Patched MYSTRAN 18.0.0 A 23-05-2026

MYSTRAN is an acronym for “My Structural Analysis” (https://www.mystran.com)


---

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
- Support for Classical Laminated Plate Theory
- OP2 Support

# Patched Distribution Notes

This branch is maintained as the active `MYSTRAN 18.0.0` patched distribution.
It tracks the working source stack used by the current Windows build and keeps
internal Codex notes, validation packs, and scratch artifacts out of the public
repository tree.

Major patch groups currently included in the source tree:

- Response spectrum modernization:
  `SDAMP`, `TABDMP1`, `DTI,SPECSEL`, `PARAM,OPTION` modal combination control,
  effective mass participation reporting, and initial `SOL SEMODES` alias
  support with explicit "Not supported yet" guards for unsupported NX-style
  features.
- LAPACK / ARPACK peel-off and helper split:
  helper entry points for dense, banded, and tridiagonal kernels, plus Lanczos
  and eigen extraction updates used by the current optimized solver path.
- Shell formulation work:
  ongoing MITC and shell buckling improvements, shell material/orientation
  cleanup, and related FEMAP/output updates.
- Composite shell support:
  active source support for laminated shell routing and benchmarked
  `CQUADR` / `CTRIAR` composite paths.
- Solid element work:
  active source updates related to the current solid element stack, including
  the newer pyramid/solid integration work present in this distribution.
- Parser and load-processing hardening:
  fixes around continuation handling, load processing, mass assembly, and
  supporting utility routines required by the current benchmarked workflows.

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
- [Join our Forums](https://mystran.com/forums) - Little activity. Mostly for archive purposes.
