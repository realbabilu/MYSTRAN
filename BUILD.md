# Building MYSTRAN 18 (public patched distribution)

This branch is published as a whole patched MYSTRAN 18 source tree. It does **not** bundle third-party binary libraries. Build scripts and wrapper sources are included, but the external `.a` / `.dll` artifacts must be provided by the builder.

## Supported Windows workflow

The primary Windows workflow for this distribution is:

- MinGW toolchain from equation.com (or another compatible MinGW GCC/GFortran toolchain)
- CMake
- out-of-tree build directory
- explicit external library paths

The root example script is:

- `example_nometis.bat`

That script creates a `build` directory, changes into it, configures CMake, and builds `mystran.exe`.

## Build requirements

1. A C/C++ and Fortran toolchain:
   - GCC/GFortran from Equation.com, or
   - Intel oneAPI HPC C++ and Fortran Compiler for Windows, or
   - another compatible GCC/GFortran toolchain on other platforms
2. A build driver:
   - `make`, `mingw32-make`, `nmake`, or `ninja`
3. `cmake`
4. A BLAS implementation:
   - OpenBLAS, MKL, AOCL, or BLIS/FLAME
5. Pre-compiled SuperLU:
   - with or without METIS, depending on your chosen build
6. Optional pre-compiled libraries:
   - MUMPS for the optional sparse solver path
   - FEAST 4.0 for the optional eigensolver path

Notes:

- `mystran.exe` in this branch is documented around an external SuperLU path.
- BLAS is required.
- MUMPS and FEAST are optional.
- CHASE is not part of the documented public build path here.

## BLAS is mandatory

This distribution must be linked with a BLAS implementation. Use one of:

- OpenBLAS
- Intel MKL
- AMD AOCL
- BLIS/FLAME

For the MinGW example in this repository, BLAS is passed explicitly with:

```bat
-DTPL_BLAS_LIBRARIES="C:/gcc/openblas32/lib/libopenblas.dll.a"
```

At runtime you must also make sure the corresponding BLAS DLL is visible on `PATH` when your BLAS package uses shared libraries.

## External libraries used in this branch

### 1. SuperLU

This public branch is centered on the external SuperLU path.

Relevant CMake options:

```bat
-DMYSTRAN_USE_EXTERNAL_SUPERLU=ON
-DMYSTRAN_EXTERNAL_SUPERLU_LIB="C:/gcc/libsuperlu/libsuperlu.a"
-DMYSTRAN_EXTERNAL_SUPERLU_INCLUDE_DIR="C:/gcc/libsuperlu/include"
-DMYSTRAN_EXTERNAL_SUPERLU_CONFIG_DIR="C:/gcc/libsuperlu"
-DMYSTRAN_EXTERNAL_SUPERLU_DRIVER="%SRC%/superlu/FORTRAN/c_fortran_dgssv.c"
```

Required pieces:

- SuperLU static or import library, typically `libsuperlu.a`
- SuperLU headers
- SuperLU config directory if your build exports one
- the Fortran/C wrapper source `c_fortran_dgssv.c`

If your SuperLU build depends on METIS, enable the METIS options shown below.

### 2. METIS (optional)

Only enable this if your SuperLU or MUMPS build actually requires external METIS/GKlib.

```bat
-DTPL_ENABLE_METISLIB=ON
-DTPL_METIS_INCLUDE_DIRS="C:/gcc/libmetis/include"
-DTPL_METIS_LIBRARIES="C:/gcc/libmetis/libmetis.a;C:/gcc/libmetis/libGKlib.a"
```

The default sample script keeps this off:

```bat
-DTPL_ENABLE_METISLIB=OFF
```

### 3. FEAST (optional)

```bat
-DMYSTRAN_USE_EXTERNAL_FEAST=ON
-DMYSTRAN_FEAST_EXTRA_LIBS="C:/gcc/feast32/libfeast.a"
```

Add any extra BLAS/LAPACK-related libraries required by your FEAST package.

### 4. MUMPS / DMUMPS (optional)

```bat
-DMYSTRAN_USE_DMUMPS_SOLVER=ON
-DMYSTRAN_DMUMPS_INCLUDE_DIR="C:/gcc/libmumps/include"
-DMYSTRAN_DMUMPS_EXTRA_LIBS="C:/gcc/mumps32_nonmpi/libdmumps.a;C:/gcc/mumps32_nonmpi/libmpiseq.a;C:/gcc/mumps32_nonmpi/libmumps_common.a;C:/gcc/mumps32_nonmpi/libpord.a;C:/gcc/mumps32_nonmpi/libsmumps.a"
```

Typical pieces needed for a non-MPI MinGW build:

- `libdmumps.a`
- `libmumps_common.a`
- `libmpiseq.a`
- `libpord.a`
- any precision-specific archives required by your package
- MUMPS headers such as `dmumps_struc.h`

### 5. FEAST wrapper note

FEAST does **not** use a small standalone wrapper source file in the same style as the external SuperLU path.

In this tree, FEAST is wired directly inside the MYSTRAN source, mainly through:

- `Source/LK4/EIGRL_EXTRACT_SOLVERS.F90`

So for FEAST, the public branch documents:

- compile-time enable with `-DMYSTRAN_USE_EXTERNAL_FEAST=ON`
- link-time library list with `-DMYSTRAN_FEAST_EXTRA_LIBS=...`
- no separate `wrapper/feast/*.c` shim is required in the current public path

### 6. Intel oneAPI on Windows

If you build third-party libraries with Intel oneAPI on Windows and then link them into a MinGW MYSTRAN build, pay close attention to symbol naming and runtime compatibility.

For the external SuperLU path, the practical convention for this branch is:

- keep `slu_Cnames.h` in an `UPCASE`-only configuration for the active build path
- remove or disable the other naming convention branches used in your local third-party package if they conflict with your Fortran/C symbol binding choice

In other words, for this public MYSTRAN branch the SuperLU interface should be treated as a single naming-convention build, not a many-convention package.

The relevant files are typically:

- `superlu/SRC/slu_Cnames.h`
- `superlu/CBLAS/slu_Cnames.h`

If your external SuperLU package was prepared separately, make sure it matches the naming convention expected by the wrapper source and your selected Fortran compiler.

## Example MinGW / equation.com build

From a shell where `gcc`, `g++`, `gfortran`, and `cmake` are available:

```bat
@echo off
setlocal

if not exist build mkdir build
cd /d build

set "SRC=..\MYSTRANSolver-18.0.0"

cmake -G "MinGW Makefiles" ^
  -DCMAKE_BUILD_TYPE=Release ^
  -DCMAKE_C_COMPILER=gcc.exe ^
  -DCMAKE_CXX_COMPILER=g++.exe ^
  -DCMAKE_Fortran_COMPILER=gfortran.exe ^
  -DCMAKE_C_FLAGS_RELEASE="-O3" ^
  -DCMAKE_CXX_FLAGS_RELEASE="-O3" ^
  -DCMAKE_Fortran_FLAGS_RELEASE="-O3 -ffree-line-length-none" ^
  -DMYSTRAN_DISABLE_NDEBUG=ON ^
  -DTPL_BLAS_LIBRARIES="C:/gcc/openblas32/lib/libopenblas.dll.a" ^
  -DMYSTRAN_USE_EXTERNAL_SUPERLU=ON ^
  -DMYSTRAN_EXTERNAL_SUPERLU_LIB="C:/gcc/libsuperlu/libsuperlu.a" ^
  -DMYSTRAN_EXTERNAL_SUPERLU_INCLUDE_DIR="C:/gcc/libsuperlu/include" ^
  -DMYSTRAN_EXTERNAL_SUPERLU_CONFIG_DIR="C:/gcc/libsuperlu" ^
  -DMYSTRAN_EXTERNAL_SUPERLU_DRIVER="%SRC%/superlu/FORTRAN/c_fortran_dgssv.c" ^
  -DTPL_ENABLE_METISLIB=OFF ^
  -DMYSTRAN_USE_EXTERNAL_FEAST=ON ^
  -DMYSTRAN_FEAST_EXTRA_LIBS="C:/gcc/feast32/libfeast.a" ^
  -DMYSTRAN_USE_DMUMPS_SOLVER=ON ^
  -DMYSTRAN_DMUMPS_INCLUDE_DIR="C:/gcc/libmumps/include" ^
  -DMYSTRAN_DMUMPS_EXTRA_LIBS="C:/gcc/mumps32_nonmpi/libdmumps.a;C:/gcc/mumps32_nonmpi/libmpiseq.a;C:/gcc/mumps32_nonmpi/libmumps_common.a;C:/gcc/mumps32_nonmpi/libpord.a;C:/gcc/mumps32_nonmpi/libsmumps.a" ^
  "%SRC%"

cmake --build . --config Release
```

The executable is produced at:

- `Binaries\mystran.exe`

## Runtime DLL note

If you linked against import libraries such as `libopenblas.dll.a`, the corresponding runtime DLLs must be available when launching `mystran.exe`.

Typical examples include:

- OpenBLAS DLL
- MinGW runtime DLLs from the selected GCC toolchain

## Wrapper sources

This repository may include source-only wrappers under `wrapper/` for public build reference. These are provided as source shims only; the actual third-party libraries are not bundled in this branch.

## Licensing note

Third-party license texts and notices used by this public distribution are collected in:

- `LICENSE/`

