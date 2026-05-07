External prebuilt SuperLU packaging note.

Files in this folder:

- `CMakeLists.txt`
- `build.bat`

Purpose:

- avoid rebuilding bundled regular SuperLU on every MYSTRAN build
- allow MYSTRAN to link a prebuilt archive such as:
  - `C:/gcc/libsuperlu/superlu-metis/libsuperlu.a`
- keep CHASE / FEAST / DMUMPS / METIS wiring in one repeatable MinGW build recipe

Live source change:

- top-level `CMakeLists.txt` now supports:
  - `MYSTRAN_USE_EXTERNAL_SUPERLU=ON`
  - `MYSTRAN_EXTERNAL_SUPERLU_LIB`
  - `MYSTRAN_EXTERNAL_SUPERLU_INCLUDE_DIR`
  - `MYSTRAN_EXTERNAL_SUPERLU_CONFIG_DIR`
  - `MYSTRAN_EXTERNAL_SUPERLU_DRIVER`

Batch file:

- `build.bat` is the prepared MinGW release configure/build command that uses:
  - external SuperLU+METIS archive
  - external CHASE
  - external FEAST
  - non-MPI MUMPS libraries
