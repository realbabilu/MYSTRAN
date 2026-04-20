# f2c Build Note

This repository was successfully built on Windows with:

- MinGW `gcc/gfortran`
- OpenBLAS supplied through `TPL_BLAS_LIBRARIES`
- regular `SuperLU`

In that configuration, the project built and `mystran.exe` ran correctly after removing the `f2c` block from `CMakeLists.txt`.

## Practical conclusion

For the tested Windows build configuration:

- `libf2c.a` was not required
- `f2c.h` was not required

The working build path was:

- native MYSTRAN Fortran sources
- SuperLU
- vendor BLAS / OpenBLAS

without building or linking the legacy `f2c` runtime.

## Why `f2c` is not needed here

The current repository is overwhelmingly native Fortran. In the tested configuration:

- MYSTRAN core is compiled directly with `gfortran`
- SuperLU is compiled as C code
- BLAS is provided externally by OpenBLAS

That means there is no confirmed need for the old `libf2c` runtime library in this build path.

`f2c` is historically associated with:

- translated Fortran-to-C runtime support
- generated headers such as `arith.h`
- legacy support functions used by `f2c`-style C code

But those were not required for this working Windows build.

## What changed in this repository

The previous `CMakeLists.txt` included logic to:

- build `arithchk`
- generate `arith.h`
- compile the `f2c` C sources into `libf2c`
- link `mystran` against `f2c`

That logic was removed for the tested Windows configuration, and the build still worked.

## Upload/repository note

For this OpenBLAS-based workflow, `f2c` support was removed from the main `CMakeLists.txt`.

Practical implication for this repository upload:

- no need to upload `f2c` sources/runtime artifacts for this build path
- no need to upload bundled `superlu`/`f2c` source trees when your target environment already provides the needed dependencies

This note is specifically about the tested OpenBLAS + normal SuperLU build path.

## Scope note

This note applies to the tested configuration only:

- Windows
- MinGW
- external BLAS / OpenBLAS
- normal SuperLU path

A different future configuration may still require `f2c.h` or related compatibility headers, especially if an internal fallback numerical path is enabled. But for the working build recorded here, neither `f2c.h` nor `libf2c.a` was needed.
