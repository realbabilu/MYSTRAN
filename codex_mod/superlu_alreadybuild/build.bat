@echo off
setlocal

rem MinGW Release build for MYSTRAN with:
rem - external prebuilt SuperLU+METIS archive
rem - external CHASE
rem - external FEAST
rem - direct DMUMPS sparse path using non-MPI MUMPS libs

set "SRC=..\MYSTRANSolver-18.0.0"

cmake -G "MinGW Makefiles" ^
  -DCMAKE_BUILD_TYPE=Release ^
  -DCMAKE_C_COMPILER=gcc.exe ^
  -DCMAKE_CXX_COMPILER=g++.exe ^
  -DCMAKE_Fortran_COMPILER=gfortran.exe ^
  -DCMAKE_C_FLAGS_RELEASE="-O3 -march=znver4" ^
  -DCMAKE_CXX_FLAGS_RELEASE="-O3 -march=znver4" ^
  -DCMAKE_Fortran_FLAGS_RELEASE="-O3 -march=znver4" ^
  -DMYSTRAN_DISABLE_NDEBUG=ON ^
  -DTPL_BLAS_LIBRARIES="C:/gcc/openblas32/lib/libopenblas.dll.a" ^
  -DMYSTRAN_USE_EXTERNAL_SUPERLU=ON ^
  -DMYSTRAN_EXTERNAL_SUPERLU_LIB="C:/gcc/libsuperlu/superlu-metis/libsuperlu.a" ^
  -DMYSTRAN_EXTERNAL_SUPERLU_INCLUDE_DIR="%SRC%/superlu/SRC" ^
  -DMYSTRAN_EXTERNAL_SUPERLU_CONFIG_DIR="C:/gcc/libsuperlu" ^
  -DMYSTRAN_EXTERNAL_SUPERLU_DRIVER="%SRC%/superlu/FORTRAN/c_fortran_dgssv.c" ^
  -DTPL_ENABLE_METISLIB=ON ^
  -DTPL_METIS_INCLUDE_DIRS="C:/gcc/libmetis/include" ^
  -DTPL_METIS_LIBRARIES="C:/gcc/libmetis/shared/libgk.dll.a;C:/gcc/libmetis/shared/libmetis.dll.a" ^
  -DMYSTRAN_USE_EXTERNAL_CHASE=ON ^
  -DMYSTRAN_CHASE_EXTRA_LIBS="C:/gcc/chase32/nompi/libchase_c.a;C:/gcc/chase32/nompi/libchase_f.a;C:/gcc/chase32/libmyfort.a;C:/gcc/chase32/libsymbols.a;stdc++;gomp" ^
  -DMYSTRAN_USE_EXTERNAL_FEAST=ON ^
  -DMYSTRAN_FEAST_EXTRA_LIBS="C:/gcc/feast32/libfeast.a" ^
  -DMYSTRAN_USE_DMUMPS_SOLVER=ON ^
  -DMYSTRAN_DMUMPS_INCLUDE_DIR="C:/gcc/libmumps/include" ^
  -DMYSTRAN_DMUMPS_EXTRA_LIBS="C:/gcc/mumps32_nonmpi/libdmumps.a;C:/gcc/mumps32_nonmpi/libmpiseq.a;C:/gcc/mumps32_nonmpi/libmumps_common.a;C:/gcc/mumps32_nonmpi/libpord.a;C:/gcc/mumps32_nonmpi/libsmumps.a" ^
  "%SRC%"

if errorlevel 1 (
  echo.
  echo CMake configure failed.
  exit /b 1
)

cmake --build . --config Release --target mystran -- -j8

if errorlevel 1 (
  echo.
  echo Build failed.
  exit /b 1
)

echo.
echo Build finished.

endlocal
