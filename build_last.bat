@echo off
setlocal EnableExtensions

rem Latest local MYSTRAN build recipe.
rem Run this file from anywhere; paths are resolved from this .bat location.
rem
rem Toolchain:
rem   C:/gcc/bin/gfortran.exe
rem
rem Sparse libraries:
rem   SuperLU built with METIS:
rem     C:/gcc/libsuperlu/superlu-metis/libsuperlu.a
rem   METIS/GK shared import libraries:
rem     C:/gcc/libmetis/shared/libgk.dll.a
rem     C:/gcc/libmetis/shared/libmetis.dll.a
rem   Runtime DLLs copied to Binaries:
rem     C:/gcc/libmetis/shared/libgk.dll
rem     C:/gcc/libmetis/shared/libmetis.dll

set "ROOT=%~dp0"
set "ROOT=%ROOT:~0,-1%"
set "BUILD=%ROOT%\build_mingw"
set "BIN=%ROOT%\Binaries"
set "SRC=..\MYSTRANSolver-18.0.0"

set "GCC=C:\gcc"
set "GCC_FWD=C:/gcc"
set "SUPERLU_LIB=%GCC_FWD%/libsuperlu/libsuperlu.a"
set "SUPERLU_INC=%ROOT%\superlu\SRC"
set "SUPERLU_CONFIG=%GCC_FWD%/libsuperlu"
set "SUPERLU_DRIVER=%ROOT%\superlu\FORTRAN\c_fortran_dgssv.c"
set "OPENBLAS=%GCC_FWD%/openblas32/lib/libopenblas.dll.a"
set "FEAST_LIBS=%GCC_FWD%/feast32/libfeast.a"
set "MUMPS_INC=%GCC_FWD%/libmumps/include"
set "MUMPS_LIBS=%GCC_FWD%/mumps32_nonmpi/libdmumps.a;%GCC_FWD%/mumps32_nonmpi/libmpiseq.a;%GCC_FWD%/mumps32_nonmpi/libmumps_common.a;%GCC_FWD%/mumps32_nonmpi/libpord.a;%GCC_FWD%/mumps32_nonmpi/libsmumps.a"
set "OPT_FLAGS=-O3 -march=znver4"
set "FOPT_FLAGS=%OPT_FLAGS% -ffree-line-length-none"

if exist "C:\Program Files\Git\usr\bin" set "PATH=C:\Program Files\Git\usr\bin;%PATH%"
if exist "C:\Program Files\Git\mingw64\bin" set "PATH=C:\Program Files\Git\mingw64\bin;%PATH%"
if exist "C:\Program Files\Git\cmd" set "PATH=C:\Program Files\Git\cmd;%PATH%"
set "PATH=%PATH%;%GCC%\bin"

echo.
echo === MYSTRAN latest local build ===
echo ROOT  = %ROOT%
echo BUILD = %BUILD%
echo.

if not exist "%BUILD%" mkdir "%BUILD%"

if exist "%BUILD%\CMakeCache.txt" if /I not "%~1"=="configure" goto build_only

cmake -S "%SRC%" -B "%BUILD%" -G "MinGW Makefiles" ^
  -DCMAKE_BUILD_TYPE=Release ^
  -DCMAKE_C_COMPILER=gcc.exe ^
  -DCMAKE_CXX_COMPILER=g++.exe ^
  -DCMAKE_Fortran_COMPILER=gfortran.exe ^
  -DCMAKE_C_FLAGS_RELEASE="%OPT_FLAGS%" ^
  -DCMAKE_CXX_FLAGS_RELEASE="%OPT_FLAGS%" ^
  -DCMAKE_Fortran_FLAGS_RELEASE="%FOPT_FLAGS%" ^
  -DMYSTRAN_DISABLE_NDEBUG=ON ^
  -DTPL_BLAS_LIBRARIES="%OPENBLAS%" ^
  -DMYSTRAN_USE_EXTERNAL_SUPERLU=ON ^
  -DMYSTRAN_EXTERNAL_SUPERLU_LIB="%SUPERLU_LIB%" ^
  -DMYSTRAN_EXTERNAL_SUPERLU_INCLUDE_DIR="%SUPERLU_INC%" ^
  -DMYSTRAN_EXTERNAL_SUPERLU_CONFIG_DIR="%SUPERLU_CONFIG%" ^
  -DMYSTRAN_EXTERNAL_SUPERLU_DRIVER="%SUPERLU_DRIVER%" ^
  -DTPL_ENABLE_METISLIB=OFF ^
  -DMYSTRAN_USE_EXTERNAL_FEAST=ON ^
  -DMYSTRAN_FEAST_EXTRA_LIBS="%FEAST_LIBS%" ^
  -DMYSTRAN_USE_DMUMPS_SOLVER=ON ^
  -DMYSTRAN_DMUMPS_INCLUDE_DIR="%MUMPS_INC%" ^
  -DMYSTRAN_DMUMPS_EXTRA_LIBS="%MUMPS_LIBS%"

if errorlevel 1 (
  echo.
  echo CMake configure failed.
  exit /b 1
)

:build_only
echo Using existing CMake cache in build_mingw.
echo Pass "configure" to this .bat to regenerate the cache.
echo.

cmake --build "%BUILD%" --target mystran --config Release -- -j8

if errorlevel 1 (
  echo.
  echo Build failed.
  exit /b 1
)

if not exist "%BIN%" mkdir "%BIN%"

echo.
echo Build finished:
echo   %BIN%\mystran.exe

endlocal
