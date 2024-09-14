MYSTRAN
=======

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

a. Build using MSYS2 for Windows
   -  Get and install msys2 https://www.msys2.org/ 
   - Launch msys2 msys base from startmenu
   - Get latest repository by type
     pacman ‐Syu
   - Get requirement
     pacman ‐S mingw‐w64‐x86_64‐gcc‐fortran 
     pacman ‐S mingw‐w64‐x86_64‐cmake
	   pacman -S mingw-w64-x86_64-make
     pacman -S mingw-w64-x86_64-f2c
	   pacman -S git
	   pacman -S mingw-w64-x86_64-openblas64 
   - Get environment setting
     export PATH="/mingw64/bin:$PATH"
     Get git clone latest mystran
     git clone https://github.com/MYSTRANsolver/MYSTRAN
  - Use existing Cmakelists.txt or modify it to use Cmakelists.pc
    mv Cmakelists.txt Cmakelists.original
    mv Cmakelists.pc Cmakelists.txt
  - Make directory for build and get into it
    mkdir build
    cd build
  - Build cmake
    OPTION 1: For internal blas build with
      cmake -G "MinGW Makefiles" ..
      mingw32-make
    OPTION 2: Using Openblas (check where it installed for example = /mingw64/lib/libopenblas_64.dll.a   
      cmake -G "MinGW Makefiles" -D"TPL_ENABLE_BLAS=TRUE" -D"USE_XSDK_DEFAULTS_DEFAULT=TRUE" -D"BLAS_LIBRARIES:PATH=/mingw64/lib/libopenblas_64.dll.a" -D"XSDK_ENABLE_Fortran=TRUE"  -DBLA_VENDOR=OpenBlas -DBUILD_SHARED_LIBS=NO  ..  -D"TPL_ENABLE_BLAS=TRUE" -D"TPL_BLAS_LIBRARIES:PATH=/mingw64/lib/libopenblas_64.dll.a"
      mingw32-make
   -For run please copy the libopenblass.dll.a and other libopenblas files to run it.  
    
b. Build using native linux gcc and gfortran.
   - Get gcc and gfortran, cmake, make, git, openblas using package installer
   - Use existing Cmakelists.txt or modify it to use Cmakelists.pc
    mv Cmakelists.txt Cmakelists.original
    mv Cmakelists.pc Cmakelists.txt
  - Make directory for build and get into it
    mkdir build
    cd build
  - Build cmake
    OPTION 1: For internal blas build with
      cmake -G "Unix Makefiles" ..
      make
    OPTION 2: Using Openblas (check where it installed for example = /usr/lib/x86_64-linux-gnu/libopenblas.so   
      cmake -G "Unix Makefiles" -D"TPL_ENABLE_BLAS=TRUE" -D"USE_XSDK_DEFAULTS_DEFAULT=TRUE" -D"BLAS_LIBRARIES:PATH=/usr/lib/x86_64-linux-gnu/libopenblas.so" -D"XSDK_ENABLE_Fortran=TRUE"  -DBLA_VENDOR=OpenBlas -DBUILD_SHARED_LIBS=NO  ..  -D"TPL_ENABLE_BLAS=TRUE" -D"TPL_BLAS_LIBRARIES:PATH=/usr/lib/x86_64-linux-gnu/libopenblas.so"
      make

c. Build using native Windows Intel fortran via command line.
  - Get Microsoft Visual Studio latest and Intel® oneAPI Base Toolkit and Intel® HPC Toolkit latest, cmake
  - Use Cmakelists.pc, rename the original
    mv Cmakelists.txt Cmakelists.original
    mv Cmakelists.pc Cmakelists.txt
  - Download zipped superlu (without metis and gklib) , f2c, and mystran latest source.
  - Extract superlu in superlu folder of extracted mystran, and f2c in f2c folder of mystran source
  - Recommendeded using IntelMKL Blas and LAPACK for compiling, instead of internal BLAS
  - Make directory for build and get into it
    mkdir build
    cd build
  - Several file need to be changed, patcher via python available soon.
    The file that changed is GET_MYSTRAN_DIR.f90
      ! INTRINSIC                                    :: GETENV ! this one commented
      CALL get_environment_variable  ( 'MYSTRAN_directory', MYSTRAN_DIR ) ! CALL GETENV ( 'MYSTRAN_directory', MYSTRAN_DIR ) ! this one commended
  - Use cmake
    cmake -G "NMake Makefiles JOM"  -D"CMAKE_FORTRAN_COMPILER=ifort.exe" -D"CMAKE_C_COMPILER=icx.exe" -D"CMAKE_CXX_COMPILER=icx.exe"  -D"TPL_ENABLE_BLAS=TRUE" -DBLA_VENDOR=Intel10_64lp -D"CMAKE_BUILD_TYPE=RELEASE" ..
    nmake

  d. Build using native Windows Intel fortran via VS STUDIO, cmake
  - Build library superlu first, extract go to SRC
    mkdir lib
    icx *.c /c /O3   /DUpcase /DXSDK_INDEX_SIZE=32 /tune:skylake /arch:skylake /Qm64  -DF77_CALL_C=UPCASE
    xilib *.obj /OUT:libsuperlu.lib
    del *.obj
    icx *.c /c /O3   /DUpcase /DXSDK_INDEX_SIZE=64 /tune:skylake /arch:skylake /Qm64   -DF77_CALL_C=UPCASE
    xilib *.obj /OUT:libsuperlu.lib
    copy *.h lib
    move *.lib lib
 -  go to fortran folder of superlu
    icx *.c /c /O3   /DUpcase /DXSDK_INDEX_SIZE=64 /tune:skylake /arch:skylake /Qm64   -DF77_CALL_C=UPCASE
  - Create console Fortran project file, drag all WITHOUT \Source\BLAS in latest mystran or \Source\Modules\LAPACK\Unresolved_Externals_Problem into source list files
  - Drag also c_fortran_dgssv.obj and libsuperlu.lib created before into source list files.
  - add configuration - linker - command line:  /heap-arrays:0 /F200000000
    or configuration - linker - command line - stack reserve : 200000000
  - activated intelmkl - configuration - Fortran - Libraries - Use Intel MKL Libraries - Parallel (preffered) or Sequential
  - Several file need to be changed, patcher via python available soon.
    The file that changed is GET_MYSTRAN_DIR.f90
      ! INTRINSIC                                    :: GETENV ! this one commented
      CALL get_environment_variable  ( 'MYSTRAN_directory', MYSTRAN_DIR ) ! CALL GETENV ( 'MYSTRAN_directory', MYSTRAN_DIR ) ! this one commended
  - Build debug setting for debug, release for faster with debugging.

  e. Build using native Windows equation.com GCC
  - get equation.com gcc install at c:\gcc , openblas binary extract at c:\gcc\openblas, extract cmake
  - put superlu in superlu folder
  - Run cmake and make
    cd MYSTRAN-15.2.1
    mkdir build
    cd build
    cmake -G "MinGW Makefiles" -D"CMAKE_MAKE_PROGRAM:PATH=c:\gcc\bin\make.exe" -DWIN32=TRUE -D"CMAKE_Fortran_COMPILER:PATH=c:\gcc\bin\gfortran.exe" -D"TPL_ENABLE_BLAS=TRUE" -D"TPL_BLAS_LIBRARIES:PATH=c:\gcc\openblas\libopenblas.DLL" -D"USE_XSDK_DEFAULTS_DEFAULT=TRUE" -D"BLAS_LIBRARIES:PATH=c:\gcc\openblas\libopenblas.DLL" -D"XSDK_ENABLE_Fortran=TRUE" -D"CMAKE_BUILD_TYPE=RELEASE" -DBLA_VENDOR=OpenBlas  ..
 


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

- NASTRAN compatibility
- OP2 Support
- Modal analysis
- Linear Static Analysis
- Linear Elastic Buckling Analysis (All But Shell Elements)
- Support for True Classical Laminate Plate Theory

# Get EXE or Make Binary

Windows EXE (executable) for can be found in the 
[MYSTRAN Releases](https://github.com/MYSTRANsolver/MYSTRAN_Releases) repository.

Static Linux binaries have been built, but releases are in work.
For now, it is better to build it yourself -- it's really
straightforward.

# Documentation

The end user documentation is located in the "User_Documents" folder of this repository.

# Four Repositories

The MYSTRAN project consits of four repositoires.

1 - This repository contains the source code, build instructions, and end user documentation.

2 - The [MYSTRAN_Resources](https://github.com/MYSTRANsolver/MYSTRAN_Resources) repository consists of files for MYSTRAN developers.
It also contains informationa and files realated to pre- and post-processors relavant to MYSTRAN.

3 - The [MYSTRAN_Releases](https://github.com/MYSTRANsolver/MYSTRAN_Releases) repository contains the most current release as well as prior releases.

4 - The [MYSTRAN_Benchmark](https://github.com/MYSTRANsolver/MYSTRAN_Benchmark) repository contains the test cases and utilities used to verify that a new build produces results consistent with prior builds and models that have been verified.


# Developmental Goals

- Implement the MITC shell elements and shell element buckling capability
- Improve OP2 output
- Creating easier ways to acquire MYSTRAN would be nice. This would include, but
  is not limited to, entry into the Arch Linux User Repository (AUR), the
  Debian Advanced Package Manager (apt), the snapcraft store (snap), the
  chocolatey package manager for Windows, an AppImage, or flatpak.
- As a longer term goal, geometric nonlinear support is desirable.

# Ways You Can Help

- Join the MYSTRAN forum and/or Discord Channel (links below)
- Contribute your MYSTRAN runs to the list of demonstration problems
- Report bugs and inconsistencies

# Community

- [Join our Forums](https://mystran.com/forums)
- [Join our Discord Channel](https://discord.gg/9k76SkHpHM)
