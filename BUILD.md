# Building MYSTRAN from source

###### Last updated 2024-09-26.

## Setting up a build environment

In order to build (compile) MYSTRAN using CMake, you first have to set up a
proper build environment (i.e. toolchain and required programs/libraries).

You can skip this part if you've done it already (or if you really know what
you're doing).

### Steps for Windows (x86_64) using Intelone CMAKE using IntelMKL Oneapi 
First, download and install MS Visual Studio Community
(https://visualstudio.microsoft.com/downloads/).
Second, download IntelOneAPI Basekit for C++ compiling with MKL and IntelOneAPI HPCkit for fortran compiling
(https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit-download.html)
(https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit-download.html)
Better to have all compiler and VS installed in drive c: otherwise some nerd setvars editing needed.
Download jom.exe put in c:\windows folder for faster compiling 
(https://github.com/qt-labs/jom)

Download superlu, preferably 5.3
https://github.com/xiaoyeli/superlu/archive/refs/tags/v5.3.0.zip
The latest version should be done but need some arguments for not compiling Metis and GTK for simplicity.
https://github.com/xiaoyeli/superlu

Extract the main mystran zipped to folder like c:\mystran, 
Make sure Source and Cmakelists.txt folder in c:\mystran\Source 
Extract the superlu in c:\mystran\Superlu, 
Make sure SRC folder and Fortran folder and Cmakelists.txt, Makefiles in the c:\mystran\Superlu

Open Command Line Intel OneApi x64 Enviroment.
1. Go to c:\mystran folder by type **`cd c:\mystran`** enter
2. just type run.bat enter for automatically built
3. Alternatively, make build folder by type **`mkdir build`** enter
4. Go to build folder by **`Cd build`** enter
5. For using DSS PARDISO and Intel MKL BLAS type  **`cmake -G "NMake Makefiles JOM"  -D"CMAKE_FORTRAN_COMPILER=ifort.exe" -D"CMAKE_C_COMPILER=icx.exe" -D"CMAKE_CXX_COMPILER=icx.exe"  -D"TPL_ENABLE_BLAS=TRUE" -DBLA_VENDOR=Intel10_64lp -D"USE_XSDK_DEFAULTS_DEFAULT=TRUE" -D"XSDK_ENABLE_Fortran=TRUE" -D"CMAKE_BUILD_TYPE=RELEASE" -DMKLDSS=TRUE  ..`** in single line
6. Alternatively,  **`cmake -G "NMake Makefiles JOM"  -D"CMAKE_FORTRAN_COMPILER=ifort.exe" -D"CMAKE_C_COMPILER=icx.exe" -D"CMAKE_CXX_COMPILER=icx.exe"  -D"TPL_ENABLE_BLAS=TRUE" -DBLA_VENDOR=Intel10_64lp -D"CMAKE_BUILD_TYPE=RELEASE" ..`** for using just Superlu
7. run with **`jom /j 12`** for 12 threads or whatever your cpu capable, alternatively just **`nmake`** for single threads compiling
8. The binaries will be at **`c:\mystran\binaries`**

### Steps for Windows (x86_64) using Intelone CMAKE for Visual Studio Enviromental Debug\Release
Install Intel OneAPI and VS Studio as above, Download also superlu first extract at coressponding folder as above.
pen Command Line Intel OneApi x64 Enviroment.
1. Go to c:\mystran folder by type **`cd c:\mystran`** enter
3. make build folder by type **`mkdir build`** enter
4. make build folder **binaries** by type **`mkdir binaries`** enter
5. Go to binaries folder by **`Cd binaries`** enter
6. Create folder Release by type **`mkdir Release`** and debug **`mkdir Debug`**
7. Go back to folder build **`Cd ..\build`** enter
8. type **`cmake -G "Visual Studio 2022"  -D"CMAKE_FORTRAN_COMPILER=ifort.exe" -D"CMAKE_C_COMPILER=icx.exe" -D"CMAKE_CXX_COMPILER=icx.exe"  -D"TPL_ENABLE_BLAS=TRUE" -DBLA_VENDOR=Intel10_64lp -D"USE_XSDK_DEFAULTS_DEFAULT=TRUE" -D"XSDK_ENABLE_Fortran=TRUE" -D"CMAKE_BUILD_TYPE=RELEASE" -DMKLDSS=TRUE  ..`** in single line for VS 2022 for using Pardiso and DSS
9. alternatively, type **`cmake -G "Visual Studio 2022"  -D"CMAKE_FORTRAN_COMPILER=ifort.exe" -D"CMAKE_C_COMPILER=icx.exe" -D"CMAKE_CXX_COMPILER=icx.exe"  -D"TPL_ENABLE_BLAS=TRUE" -DBLA_VENDOR=Intel10_64lp -D"USE_XSDK_DEFAULTS_DEFAULT=TRUE" -D"XSDK_ENABLE_Fortran=TRUE" -D"CMAKE_BUILD_TYPE=RELEASE" -DMKLDSS=false  ..`** in single line for VS 2022 for not using DSS Pardiso


### Steps for Windows (x86_64) using MSYS - CBLAS and internal BLAS

First, download and install MSYS2 from the
[official site](https://www.msys2.org/).

Open the MSYS2 terminal and run the following commands:

  1. **`pacman -Syu`**
This updates repository information and installed packages, and might require
you close and reopen MSYS2 terminals.
  2. **`pacman -S mingw-w64-x86_64-gcc-fortran mingw-w64-x86_64-cmake mingw-w64-x86_64-make git`**
This installs the required compilers (the GNU C and Fortran compilers), CMake
itself, and `git`.
  3. **`export PATH="/mingw64/bin:$PATH"`**
This makes the MinGW toolchain programs (such as `make` and the compilers)
visible so CMake can find them more easily. Note that this command's effects
are lost when you reopen the terminal, so you might want to append it to your
`~/.bashrc` to save time.

### Steps for Linux (any)

Follow your distribution's steps to install the following programs/libraries:
  - **`gcc`**
  - **`g++`**
  - **`gfortran`**
  - **`make`**
  - **`cmake`**
  - **`git`**

All of those are fairly common, so get in touch in the MYSTRAN Forums or
MYSTRAN Discord if you have trouble installing any of them. Also, note that
most distros have a "base" package group for developers (e.g. Arch's
`base-devel` or Ubuntu's `build-essential`) that includes necessary tooling
such as `gcc` and `make`. If that's the case, install it!

If your distribution doesn't ship CMake 3.18+ yet, check if your distro has a
some sort of testing/unstable channel before attempting to
[install it manually](https://cmake.org/install/).


For WSL (Linux for Windows)
===========================
Mystran won't work with Ubuntu 20.04, hasn't been tested on 22.04 and should work on 24.04 (what we're testing).

If you're upgrading your WSL, open PowerShell as Administrator and run:
```
wsl --update
wsl --install --distribution Ubuntu-24.04
```

Now that you've got into a modern version of Ubuntu
```
sudo apt update
sudo apt upgrade
apt install gcc g++ gfortran make cmake git
```

---
### Steps for Linux in WSL or any
A. Using AOCL Clang and Flang with OpenBLAS
1. `sudo apt-get update` to get updated apt library
2. `sudo apt install ./aocc-compiler-4.2.0_1_amd64.deb` to get aocc compiler optimized for Ryzen
3. `sudo apt install ./aocl-linux-aocc-4.2.0_1_amd64.deb` to BLIS a CBLAS optimized for Ryzen
4. Get CMAKE with `sudo apt install cmake g++ make`
5. Get Git by `sudo apt install git`
6. Prepare enviroment for AOCC AOCL
   `cd /opt/AMD`
   `cd aocc-compiler-4.2.0/`
   `sudo ./install.sh`
   `source /opt/AMD/setenv_AOCC.sh`
   `which clang`
   `clang -v`
   `cd /opt/AMD/aocl/aocl-linux-aocc-4.2.0/aocc`
   `sudo ./set_aocl_interface_symlink.sh`
7. Get OpenBLAS compiled   
   `sudo apt-cache search openblas`
   `sudo apt -y install libopenblas0-pthread` for x86 Openblas
   `sudo apt -y install libopenblas64-0-pthread` for x64 Openblas
8.  create mystran folder and `git clone https://github.com/realbabilu/MYSTRAN` or any version mystran you want
9.  create build folder and go to build folder, create mystran without DSS Pardiso using Openblas
   'cmake -G "Unix Makefiles" -DBLA_VENDOR=OpenBLAS  .. -D"TPL_ENABLE_BLAS=TRUE"  -D"BLAS_LIBRARIES:PATH=/usr/lib/x86_64-linux-gnu/openblas64-pthread/libopenblas64.so.0" -D"CMAKE_FORTRAN_COMPILER=flang" -D"CMAKE_C_COMPILER=clang" -D"CMAKE_CXX_COMPILER=clang"' 

B. Using gcc and gfortran with OpenBLAS 
1. Do point 1 above, to update apt-get
2. get gcc`sudo apt install gcc` 
3. get gfortran `sudo apt install gfortran`
4. get cmake like above
5. get git like above
6. get Openblas like point 7.
7. create mystran folder then `git clone https://github.com/realbabilu/MYSTRAN` or any version mystran you want
8. create build folder and go to build folder, create mystran without DSS Pardiso using Openblas
    `cmake -G "Unix Makefiles" -DBLA_VENDOR=OpenBLAS  .. -DMKLDSS=off -D"TPL_ENABLE_BLAS=TRUE"  -D"BLAS_LIBRARIES:PATH=/usr/lib/x86_64-linux-gnu/openblas64-pthread/libopenblas64.so.0"`
 


## Building MYSTRAN

If your build environment is already set up, building MYSTRAN is quite
straightforward.

### Steps for Windows (any)

  1. Open the MSYS2 shell.
  2. Re-run step #3 of the previous section if needed.
  3. Fetch the source code if you haven't already. If you're using Git, you can
  clone the repo with
  **`git clone https://github.com/MYSTRANsolver/MYSTRAN.git`**.
  4. Move the terminal to the MYSTRAN folder. If you've just run `git clone`,
     just do a **`cd MYSTRAN`**.
  5. Generate the build scripts by running **`cmake -G "MinGW Makefiles" .`**.
  6. Compile with **`mingw32-make`**. If you have an N-core processor, running
  **`mingw32-make -Oline -jN`** will probably be much faster. A good choice of N is
  printed in the previous step, right before the end. The `-Oline` argument prevents garbled output when `N` > 1.
  7. The executable will reside at **`Binaries/mystran.exe`**.

### Steps for Linux (any)

  1. Open a terminal.
  2. Fetch the source code if you haven't already. If you're using Git, you can
  clone the repo with
  **`git clone https://github.com/MYSTRANsolver/MYSTRAN.git`**.
  3. Move the terminal to the MYSTRAN folder. If you've just run `git clone`,
  just do a **`cd MYSTRAN`**.
  1. Generate the build scripts by running **`cmake .`**.
  2. Compile with **`make`**. If you have an N-core processor, running
  **`make -jN`** will probably me much faster. A good choice of N is printed in
  the previous step, right before the end. You can also find the number of
  cores/threads with the `nproc` command (not all distros ship it
  out-of-the-box though).
  1. The executable will reside at **`Binaries/mystran`**.

---

## Troubleshooting

While this process is meant to be straightforward, here is a list of some of
the more common issues that can arise. Other issues users find might be added
here if they're not too specific.

If your issue isn't here, you can always ask for help at the
[MYSTRAN forums](https://www.mystran.com/forums/) or the
[Discord server](https://discord.gg/9k76SkHpHM)

---

### "I'm getting "file not found" errors when running the step #2 setup command!"

Run a **`pacman -Syyu`** (note the two 'y's) and try again.

---

### "CMake is complaining about not being able to find the toolchain or the Fortran compiler or the "make" command!"

Try running the commands `make`/`mingw32-make`, `gcc`, and `gfortran`. If any
of these comes up as a "command not found", make sure they've been installed.
If you're **sure** they are, they might not be in the PATH.

Windows users, have a look at step #3 of the setup. Linux users, check out your
distro documentation, because whatever's happening should not be happening at
all.

---

### "CMake complains about `ARCHIVE_EXTRACT`!"

Check out the output of `cmake --version`. You must have version 3.18 or newer.
If you don't, first ensure it's up to date -- perform a system-wide update.
Windows users should not find this issue relevant -- MSYS2 ships CMake 3.27.1
as of this writing. Linux users should use their own package manager.

If your system is up to date and you still run into this issue, that means your
distro ships CMake 3.17 or older. Bad luck there. Here's what you can do:

  1. Enable a testing/unstable package channel (not all distros have one)
  2. Install the latest CMake [manually](https://cmake.org/install/)
  (might piss off your package manager)
  1. Download and extract `libf2c.zip` yourself, and comment out the
  `ARCHIVE_EXTRACT` stuff in `CMakeLists.txt`.

---

### "I'm getting random SuperLU build errors!"

SuperLU is included as a submodule. A recent update to the submodule might
require a clean build. Run `make clean` and delete the `superlu` subdirectory
and run the appropriate `cmake` command again.

---

### "I'm getting cryptic linker errors related to BLAS!"

SuperLU requires BLAS. Its build script can look for and link against your
system's installed BLAS implementation (we recommend OpenBLAS). However, your
install might be lacking the appropriate static (`.a`) library files.

If you don't know how to fix that and just want to build, you can use the
integrated BLAS subroutines bundled with the SuperLU source. To do that, run
the appropriate `cmake` command with the extra option
`-Denable_internal_blaslib=YES` *before* the `.` argument.

Please be aware that the bundled CBLAS might be slow when compared to a proper
BLAS install. That might have an impact on the time it takes to run larger
models.

---

### "I want to build offline, but the CMake script attempts to download stuff!"

Download the `superlu` submodule and `libf2c.zip` beforehand, and you should be
fine.

---

### "The terminal output is garbled during compilation!"

Multiple threads are printing to standard output simultaneously. That issue can
sometimes happen as a result of running `make` instead of `mingw32-make` on
Windows, but it can affect both. It's annoying, but harmless.

However, if you *really* need compiler output to be readable, ensure `make`
only runs with one thread by passing the option `-j1`. This will make
compilation slower, but at least you'll be able to read the output.

And if it's errors you're looking for, you can build fast with `-j[number]`,
and then `-j1` just to see the error again.

---

If your issue isn't here, you can always ask for help at the
[MYSTRAN forums](https://www.mystran.com/forums/) or the
[Discord server](https://discord.gg/9k76SkHpHM)
