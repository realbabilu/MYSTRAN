rmdir build /s
mkdir build
cd build
cmake -G "NMake Makefiles JOM"  -D"CMAKE_FORTRAN_COMPILER=ifort.exe" -D"CMAKE_C_COMPILER=icx.exe" -D"CMAKE_CXX_COMPILER=icx.exe"  -D"TPL_ENABLE_BLAS=TRUE" -DBLA_VENDOR=Intel10_64lp -D"CMAKE_BUILD_TYPE=RELEASE" ..   
jom /j 16
