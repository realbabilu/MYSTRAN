# CMakeLists.txt for compiling MYSTRAN with gfortran in a GNU environment
# based on an older CMakeLists.txt by ceanwang@gmail.com
# made to work again by Bruno Borges Paschoalinoto (2020)

# set up basic project info
cmake_minimum_required(VERSION 3.18)
include(CheckFunctionExists)
include(CheckFortranSourceCompiles)

project(Mystran)
enable_language(Fortran)

if(WIN32)
  enable_language(RC)
endif()

find_program(MYSTRAN_MAKE_PROGRAM NAMES mingw32-make make)

# ! --- hijack mode for extern blas begin--- !
check_fortran_source_compiles(
  "
      use iso_c_binding, only: c_int
      integer, parameter :: ok = selected_int_kind(9)
      integer(c_int) :: i
      integer :: assert_size(merge(1, -1, storage_size(i) == 32))
      end
  "
  MYSTRAN_FORTRAN_C_INT_IS_32BIT
)

if(NOT MYSTRAN_FORTRAN_C_INT_IS_32BIT)
  message(FATAL_ERROR "OpenBLAS hijack requires 32-bit Fortran INTEGER/C int ABI.")
endif()




###

set(MYSTRAN_OPENBLAS_IMPORT_LIB "" CACHE FILEPATH
    "OpenBLAS import library for BLAS hijack mode")

if(NOT MYSTRAN_OPENBLAS_IMPORT_LIB AND TPL_BLAS_LIBRARIES)
  set(MYSTRAN_OPENBLAS_IMPORT_LIB "${TPL_BLAS_LIBRARIES}")
endif()

if(NOT MYSTRAN_OPENBLAS_IMPORT_LIB)
  find_library(MYSTRAN_OPENBLAS_IMPORT_LIB
    NAMES openblas libopenblas libopenblas.dll.a openblas.dll.a
  )
endif()

if(NOT MYSTRAN_OPENBLAS_IMPORT_LIB)
  message(FATAL_ERROR
    "Could not find OpenBLAS import library. "
    "Set TPL_BLAS_LIBRARIES or MYSTRAN_OPENBLAS_IMPORT_LIB manually."
  )
endif()

set(TPL_BLAS_LIBRARIES "${MYSTRAN_OPENBLAS_IMPORT_LIB}" CACHE FILEPATH
    "BLAS library used by MYSTRAN" FORCE)
set(MYSTRAN_USE_OPENBLAS ON CACHE BOOL "Use OpenBLAS DLL import library" FORCE)
set(enable_internal_blaslib OFF CACHE BOOL "Disable internal BLAS/f2c path" FORCE)
set(BLAS_FOUND TRUE)
message(STATUS "Using OpenBLAS for BLAS calls (${TPL_BLAS_LIBRARIES})")



# ! --- hijack mode for extern blas end--- !

# basic compiler and output options
set(CMAKE_SOURCE_DIR "${PROJECT_SOURCE_DIR}/Source")
set(BLAS_SOURCE_DIR "${PROJECT_SOURCE_DIR}/BLAS")
set(PROJECT_BINARY_DIR "${PROJECT_SOURCE_DIR}/Binaries")
set(CMAKE_FLAGS "-c -fbacktrace")

# set(CMAKE_BUILD_TYPE Debug)

# set some dirs
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib")
set(CMAKE_Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/mod")
include_directories("${CMAKE_SOURCE_DIR}/INCLUDE")

# uncomment this to debug
# set(CMAKE_VERBOSE_MAKEFILE true)

# register the Profile build type, just in case
set(CMAKE_CONFIGURATION_TYPES
  Debug
  Release
  RelWithDebInfo
  MinSizeRel
  Profiling
  Deterministic
  AVX2
  CACHE STRING "" FORCE
)

# suppress cmake warnings for superlu
if(NOT DEFINED CMAKE_SUPPRESS_DEVELOPER_WARNINGS)
  set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS 1 CACHE INTERNAL "No dev warnings")
endif()

# recommend the appropriate make command
if(WIN32)
  set(RECOMMENDED_MAKE "mingw32-make")
else()
  set(RECOMMENDED_MAKE "make")
endif()

# submodules (i.e. SuperLU)
find_package(Git QUIET)

if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
  option(GIT_SUBMODULE "Check submodules during build" ON)

  if(GIT_SUBMODULE)
    message(STATUS "Updating submodules")
    execute_process(
      COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive --force
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      RESULT_VARIABLE GIT_SUBMOD_RESULT
    )

    if(NOT GIT_SUBMOD_RESULT EQUAL "0")
      message(
        FATAL_ERROR
        "git submodule update --init --recursive --force failed"
      )
    endif()
  endif()
endif()

set(SUPERLU_DIR "${PROJECT_SOURCE_DIR}/superlu")
set(SUPERLU_MT_DIR "${PROJECT_SOURCE_DIR}/superlu_mt")
set(SUPERLU_PATCHES_DIR "${PROJECT_SOURCE_DIR}/superlu_mt_patches")

set(enable_examples OFF FORCE)
set(enable_tests OFF FORCE)
set(enable_doc OFF FORCE)

if(USE_SUPERLU_MT)
  message(STATUS "Will link against faster SuperLU_MT.")
  set(PLAT "_OPENMP" CACHE STRING "threading flavor _PTHREAD/_OPENMP" FORCE)
  message(STATUS "Applying patches to SuperLU_MT...")
  file(GLOB_RECURSE PATCH_FILES RELATIVE "${SUPERLU_PATCHES_DIR}" "${SUPERLU_PATCHES_DIR}/*")

  foreach(file IN LISTS PATCH_FILES)
    set(src "${SUPERLU_PATCHES_DIR}/${file}")
    set(dst "${SUPERLU_MT_DIR}/${file}")

    # ensure the destination directory exists
    get_filename_component(dst_dir "${dst}" DIRECTORY)
    file(MAKE_DIRECTORY "${dst_dir}")

    # print the file being copied
    message(STATUS "  Patching: ${file}")

    # copy the file
    file(COPY "${src}" DESTINATION "${dst_dir}")
  endforeach()

  message(STATUS "SuperLU_MT patched.")

  if(NOT EXISTS "${SUPERLU_MT_DIR}/CMakeLists.txt")
    message(
      FATAL_ERROR
      "The SuperLU_MT submodule was not downloaded. You have to do it manually!"
    )
  endif()

  if(WIN32)
    # force static linking to OpenMP library
    # temporarily set CMAKE_FIND_LIBRARY_SUFFIXES to prefer .a over .dll.a
    set(CMAKE_FIND_LIBRARY_SUFFIXES_BACKUP ${CMAKE_FIND_LIBRARY_SUFFIXES})
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")

    set(OpenMP_C_FLAGS "-fopenmp" CACHE STRING "" FORCE)
    set(OpenMP_C_LIB_NAMES "gomp" CACHE STRING "" FORCE)

    # find the static gomp library explicitly
    find_library(GOMP_STATIC_LIB NAMES libgomp.a gomp PATHS /mingw64/lib NO_DEFAULT_PATH)

    if(GOMP_STATIC_LIB)
      set(OpenMP_gomp_LIBRARY "${GOMP_STATIC_LIB}" CACHE FILEPATH "" FORCE)
    endif()

    find_package(OpenMP REQUIRED)

    # restore original library suffixes
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES_BACKUP})
  endif()

  include_directories("${SUPERLU_MT_DIR}/SRC")
  add_subdirectory(${SUPERLU_MT_DIR})
  set(SLU_DRIVER "${SUPERLU_MT_DIR}/SRC/c_fortran_pdgssv.c")
else()
  message(STATUS "Will link against regular SuperLU.")

  if(NOT EXISTS "${SUPERLU_DIR}/CMakeLists.txt")
    message(
      FATAL_ERROR
      "The SuperLU submodule was not downloaded. You have to do it manually!"
    )
  endif()

  include_directories("${SUPERLU_DIR}/SRC")
  add_subdirectory(${SUPERLU_DIR})
  set(SLU_DRIVER "${SUPERLU_DIR}/FORTRAN/c_fortran_dgssv.c")
endif()

if(NOT MYSTRAN_USE_OPENBLAS)
  # f2c stuff
  set(F2C_DIR "${PROJECT_SOURCE_DIR}/f2c")
  set(F2C_INCLUDE_DIR "${F2C_DIR}/include")
  set(F2C_FN "${F2C_DIR}/libf2c.zip")
  set(F2C_URL "https://www.netlib.org/f2c/libf2c.zip")

  # download f2c
  if(NOT EXISTS ${F2C_DIR})
    message(STATUS "Downloading libf2c source from ${F2C_URL}...")
    make_directory("${F2C_DIR}")
    file(DOWNLOAD ${F2C_URL} ${F2C_FN} TIMEOUT 60 STATUS DOWNLOAD_STATUS)

    # Check if download was successful.
    list(GET DOWNLOAD_STATUS 0 STATUS_CODE)
    list(GET DOWNLOAD_STATUS 1 ERROR_MESSAGE)

    if(${STATUS_CODE} EQUAL 0)
      message(STATUS "Done downloading libf2c.")
    else()
      # Exit CMake if the download failed, printing the error message.
      file(REMOVE_RECURSE ${F2C_DIR})
      message(FATAL_ERROR "Error downloading libf2c: ${ERROR_MESSAGE}")
    endif()
  endif()

  # extract libf2c source
  file(ARCHIVE_EXTRACT INPUT ${F2C_FN} DESTINATION ${F2C_DIR})

  # prepare libf2c header files
  file(GLOB_RECURSE F2C_PREHEADERS "${F2C_DIR}/*.h0")

  foreach(H0 ${F2C_PREHEADERS})
    string(REGEX REPLACE "[.]h0$" ".h" H0_R ${H0})
    file(RENAME "${H0}" "${H0_R}")
    file(COPY "${H0_R}" DESTINATION "${F2C_INCLUDE_DIR}")
  endforeach()

  # get a load of this: f2c generates its own "arith.h" on the fly
  # so we gotta compile arithchk and run it
  set(F2C_ARITHCHK_SRC "${F2C_DIR}/arithchk.c")
  set(F2C_ARITHCHK_BIN "${F2C_DIR}/arithchk")

  if(WIN32)
    set(F2C_ARITHCHK_BIN "${F2C_ARITHCHK_BIN}.exe")
  endif()

  set(F2C_ARITH_H "${F2C_INCLUDE_DIR}/arith.h")
  set_source_files_properties(
    ${F2C_ARITHCHK_SRC} PROPERTIES COMPILE_FLAGS "-DNO_LONG_LONG -DNO_FPINIT"
  )
  add_executable(arithchk ${F2C_ARITHCHK_SRC})
  target_link_libraries(arithchk m)
  set_target_properties(
    arithchk PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${F2C_DIR}
  )
  add_custom_command(
    OUTPUT ${F2C_ARITH_H}
    COMMAND ${F2C_ARITHCHK_BIN} > ${F2C_ARITH_H}
    DEPENDS ${F2C_ARITHCHK_BIN}
  )

  # add libf2c to the compilation procedures
  include_directories(${F2C_INCLUDE_DIR})
  set(F2C_CFILES
    abort_.c backspac.c c_abs.c c_cos.c c_div.c c_exp.c c_log.c c_sin.c
    c_sqrt.c cabs.c close.c ctype.c d_abs.c d_acos.c d_asin.c d_atan.c
    d_atn2.c d_cnjg.c d_cos.c d_cosh.c d_dim.c d_exp.c d_imag.c d_int.c
    d_lg10.c d_log.c d_mod.c d_nint.c d_prod.c d_sign.c d_sin.c d_sinh.c
    d_sqrt.c d_tan.c d_tanh.c derf_.c derfc_.c dfe.c dolio.c dtime_.c
    due.c ef1asc_.c ef1cmc_.c endfile.c erf_.c erfc_.c err.c etime_.c
    exit_.c f77_aloc.c f77vers.c fmt.c fmtlib.c ftell_.c ftell64_.c
    getarg_.c getenv_.c h_abs.c h_dim.c h_dnnt.c h_indx.c h_len.c h_mod.c
    h_nint.c h_sign.c hl_ge.c hl_gt.c hl_le.c hl_lt.c i77vers.c i_abs.c
    i_dim.c i_dnnt.c i_indx.c i_len.c i_mod.c i_nint.c i_sign.c iargc_.c
    iio.c ilnw.c inquire.c l_ge.c l_gt.c l_le.c l_lt.c lbitbits.c lbitshft.c
    lread.c lwrite.c main.c open.c pow_ci.c pow_dd.c pow_di.c pow_hh.c
    pow_ii.c pow_qq.c pow_ri.c pow_zi.c pow_zz.c qbitbits.c qbitshft.c
    r_abs.c r_acos.c r_asin.c r_atan.c r_atn2.c r_cnjg.c r_cos.c r_cosh.c
    r_dim.c r_exp.c r_imag.c r_int.c r_lg10.c r_log.c r_mod.c r_nint.c
    r_sign.c r_sin.c r_sinh.c r_sqrt.c r_tan.c r_tanh.c rdfmt.c rewind.c
    rsfe.c rsli.c rsne.c s_cat.c s_cmp.c s_copy.c s_paus.c s_rnge.c
    s_stop.c sfe.c sig_die.c signal_.c sue.c system_.c typesize.c uio.c
    uninit.c util.c wref.c wrtfmt.c wsfe.c wsle.c wsne.c xwsne.c z_abs.c
    z_cos.c z_div.c z_exp.c z_log.c z_sin.c z_sqrt.c
  )
  list(TRANSFORM F2C_CFILES PREPEND "${F2C_DIR}/")
  add_definitions(-DINTEGER_STAR_8)
  add_library(f2c ${F2C_CFILES} ${F2C_ARITH_H})

  # add some extra win32 flags for libf2c
  if(WIN32)
    add_definitions(-DUSE_CLOCK -DMSDOS)
  endif()
endif()

# set some extra vars for MSYS builds to make the binary portable
if(WIN32)
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static -static-libgcc -static-libstdc++")
endif()

# collect all fortran source files
file(GLOB_RECURSE ALL_FORTRAN_FILES
  "${CMAKE_SOURCE_DIR}/*.f"
  "${CMAKE_SOURCE_DIR}/*.F"
  "${CMAKE_SOURCE_DIR}/*.f90"
  "${CMAKE_SOURCE_DIR}/*.F90"
  "${CMAKE_SOURCE_DIR}/*.f95"
  "${CMAKE_SOURCE_DIR}/*.F95"
  "${CMAKE_SOURCE_DIR}/*.f03"
  "${CMAKE_SOURCE_DIR}/*.F03"
)

# ! --- MUMPS_COO add begin --- !

if(EXISTS "${PROJECT_SOURCE_DIR}/Source/Modules/LAPACK/LAPACK_LIN_EQN_DPB.f")
  set_source_files_properties(
    "${PROJECT_SOURCE_DIR}/Source/Modules/LAPACK/LAPACK_LIN_EQN_DPB.f"
    PROPERTIES Fortran_FORMAT FIXED
  )
endif()

if(EXISTS "${PROJECT_SOURCE_DIR}/Source/Modules/LAPACK/LAPACK_LANCZOS_EIG.f")
  set_source_files_properties(
    "${PROJECT_SOURCE_DIR}/Source/Modules/LAPACK/LAPACK_LANCZOS_EIG.f"
    PROPERTIES Fortran_FORMAT FIXED
  )
endif()

option(MYSTRAN_USE_DMUMPS_SOLVER "Enable sequential non-MPI MUMPS sparse solver path" OFF)
set(MYSTRAN_DMUMPS_INCLUDE_DIR "C:/gcc/mumps32_nonmpi/include" CACHE PATH
  "Directory containing dmumps_struc.h and mpif.h for sequential MUMPS.")
set(MYSTRAN_DMUMPS_EXTRA_LIBS
  "C:/gcc/mumps32_nonmpi/libdmumps.a;C:/gcc/mumps32_nonmpi/libmumps_common.a;C:/gcc/mumps32_nonmpi/libpord.a;C:/gcc/mumps32_nonmpi/libmpiseq.a"
  CACHE STRING "Static libraries for sequential MUMPS sparse solver path.")

foreach(MYSTRAN_CPP_F90
    "${PROJECT_SOURCE_DIR}/Source/LK3/LINK3.f90"
    "${PROJECT_SOURCE_DIR}/Source/LK4/EIGRL_EXTRACT_SOLVERS.F90"
    "${PROJECT_SOURCE_DIR}/Source/Modules/DMUMPS_STUF.f90"
)
  if(EXISTS "${MYSTRAN_CPP_F90}")
    set_source_files_properties(
      "${MYSTRAN_CPP_F90}"
      PROPERTIES Fortran_PREPROCESS ON COMPILE_OPTIONS "-cpp"
    )
  endif()
endforeach()

# ! --- MUMPS_COO add end --- !

# ! --- FEAST V4.0 external hook --- begin !
option(MYSTRAN_USE_EXTERNAL_FEAST "Enable native FEAST backend when FEAST v4.0 is provided locally" OFF)
set(MYSTRAN_FEAST_ROOT "${PROJECT_SOURCE_DIR}/feast" CACHE PATH
  "Root directory where the user unpacked FEAST v4.0 (expects src/ and lib/<arch>/).")
set(MYSTRAN_FEAST_ARCH "x64" CACHE STRING
  "FEAST library architecture subdirectory under lib/.")
set(MYSTRAN_FEAST_F90 "gfortran" CACHE STRING
  "Fortran compiler name passed to FEAST's own makefile.")
set(MYSTRAN_FEAST_BUILD_FLAGS "MKL=no F90FLAGS=-O3 -fopenmp -ffree-line-length-none -ffixed-line-length-none -cpp -fallow-argument-mismatch" CACHE STRING
  "Additional flags passed to FEAST's own makefile when building libfeast.a.")

set(MYSTRAN_FEAST_REAL_ROOT "${MYSTRAN_FEAST_ROOT}")
if(EXISTS "${MYSTRAN_FEAST_ROOT}/FEAST/4.0/src")
  set(MYSTRAN_FEAST_REAL_ROOT "${MYSTRAN_FEAST_ROOT}/FEAST/4.0")
elseif(EXISTS "${MYSTRAN_FEAST_ROOT}/4.0/src")
  set(MYSTRAN_FEAST_REAL_ROOT "${MYSTRAN_FEAST_ROOT}/4.0")
endif()

set(MYSTRAN_FEAST_LIB "${MYSTRAN_FEAST_REAL_ROOT}/lib/${MYSTRAN_FEAST_ARCH}/libfeast.a" CACHE FILEPATH
  "Expected FEAST static library path after FEAST v4.0 is built.")
set(MYSTRAN_FEAST_LIB_FALLBACK "${MYSTRAN_FEAST_REAL_ROOT}/src/libfeast.a")

if(EXISTS "${PROJECT_SOURCE_DIR}/Source/LK4/EIGRL_EXTRACT_SOLVERS.F90")
  set_source_files_properties(
    "${PROJECT_SOURCE_DIR}/Source/LK4/EIGRL_EXTRACT_SOLVERS.F90"
    PROPERTIES Fortran_PREPROCESS ON COMPILE_OPTIONS "-cpp"
  )
endif()

if(EXISTS "${PROJECT_SOURCE_DIR}/Source/LK4/EIG_LANCZOS_FEAST.f90")
  set_source_files_properties(
    "${PROJECT_SOURCE_DIR}/Source/LK4/EIG_LANCZOS_FEAST.f90"
    PROPERTIES Fortran_PREPROCESS ON COMPILE_OPTIONS "-cpp"
  )
endif()
# ! --- FEAST V4.0 external hook --- end !

# same BLAS-finding subroutine as SuperLU
if(MYSTRAN_USE_OPENBLAS)
  message(STATUS "Using OpenBLAS for BLAS calls (${TPL_BLAS_LIBRARIES})")
  set(missing_blas_src "${BLAS_SOURCE_DIR}/XERBLA.f")
elseif(NOT enable_internal_blaslib)
  if(TPL_BLAS_LIBRARIES)
    set(BLAS_FOUND TRUE)
  else()
    find_package(BLAS)

    if(BLAS_FOUND)
      set(TPL_BLAS_LIBRARIES "${BLAS_LIBRARIES}" CACHE FILEPATH
        "Set from FindBLAS.cmake BLAS_LIBRARIES." FORCE)
    endif()
  endif()
endif()

if(NOT BLAS_FOUND)
  message(STATUS "BLAS not found, building local BLAS")
  list(
    APPEND blas_fns dgemm dgemv dlamch dlanst dscal dsteqr dsterf dswap dtrsm
    dtrtri ilaenv lsame xerbla
  )

  foreach(fname IN LISTS blas_fns)
    # message(STATUS "Checking for BLAS subr ${fname}")
    check_function_exists("${fname}" BLAS_FN_EXISTS)
    set(BLAS_FN_EXISTS CACHE "1" STRING)

    if(NOT BLAS_FN_EXISTS)
      # message(STATUS "BLAS subr ${fname} not found")
      string(TOUPPER ${fname} fname_upper)
      list(APPEND missing_blas_src "${BLAS_SOURCE_DIR}/${fname_upper}.f")
      list(APPEND missing_blas_fns ${fname})
    endif()

    unset(BLAS_FN_EXISTS CACHE)
  endforeach()

  # if any subroutines have bene found, create an inner blas library
  list(LENGTH missing_blas_fns MISSING_FNS_TOTAL)

  if(MISSING_FNS_TOTAL GREATER 0)
    if(MISSING_FNS_TOTAL GREATER 1)
      message(
        STATUS
        "BLAS subrs (${missing_blas_fns}) are absent and will be built locally."
      )
    else()
      message(
        STATUS
        "BLAS subr ${missing_blas_fns} is absent and will be built locally."
      )
    endif()
  endif()
else()
  message(STATUS "Using system BLAS (${TPL_BLAS_LIBRARIES})")
endif()

# prepare the main executable, linked against the specifics and the m
# it appears utils used to be a module, but that is no longer the case?
# file(GLOB UTIL_FILES "${CMAKE_SOURCE_DIR}/UTIL/*.f*")
# file(GLOB MAIN_FILES "${CMAKE_SOURCE_DIR}/MAIN/*.[fF]*")
add_executable(
  mystran
  ${ALL_FORTRAN_FILES}
  ${SLU_DRIVER}
  ${missing_blas_src}
)

# ! --- FEAST V4.0 external hook --- begin !
if(MYSTRAN_USE_EXTERNAL_FEAST)
  if(NOT EXISTS "${MYSTRAN_FEAST_REAL_ROOT}/src")
    message(FATAL_ERROR
      "MYSTRAN_USE_EXTERNAL_FEAST=ON but FEAST src/ was not found at "
      "${MYSTRAN_FEAST_ROOT}. Download FEAST v4.0, unpack it under a local "
      "'feast' folder, then point MYSTRAN_FEAST_ROOT there.")
  endif()

  if(NOT MYSTRAN_MAKE_PROGRAM)
    message(FATAL_ERROR
      "MYSTRAN_USE_EXTERNAL_FEAST=ON but no make program was found. "
      "Install make/mingw32-make so CMake can build FEAST v4.0 from feast/src.")
  endif()

  set(MYSTRAN_FEAST_EFFECTIVE_LIB "${MYSTRAN_FEAST_LIB}")
  if(NOT EXISTS "${MYSTRAN_FEAST_EFFECTIVE_LIB}" AND EXISTS "${MYSTRAN_FEAST_LIB_FALLBACK}")
    set(MYSTRAN_FEAST_EFFECTIVE_LIB "${MYSTRAN_FEAST_LIB_FALLBACK}")
  endif()

  if(EXISTS "${MYSTRAN_FEAST_EFFECTIVE_LIB}")
    add_custom_target(mystran_feast_build)
  else()
    add_custom_command(
      OUTPUT "${MYSTRAN_FEAST_LIB}"
      COMMAND ${CMAKE_COMMAND} -E env PWD=${MYSTRAN_FEAST_REAL_ROOT}/src
              ${MYSTRAN_MAKE_PROGRAM} ARCH=${MYSTRAN_FEAST_ARCH} F90=${MYSTRAN_FEAST_F90} ${MYSTRAN_FEAST_BUILD_FLAGS} feast
      COMMAND ${CMAKE_COMMAND} -E make_directory "${MYSTRAN_FEAST_REAL_ROOT}/lib/${MYSTRAN_FEAST_ARCH}"
      COMMAND ${CMAKE_COMMAND} -E copy_if_different "${MYSTRAN_FEAST_LIB_FALLBACK}" "${MYSTRAN_FEAST_LIB}"
      WORKING_DIRECTORY "${MYSTRAN_FEAST_REAL_ROOT}/src"
      COMMENT "Building FEAST v4.0 static library at ${MYSTRAN_FEAST_LIB}"
      VERBATIM
    )
    add_custom_target(mystran_feast_build DEPENDS "${MYSTRAN_FEAST_LIB}")
    set(MYSTRAN_FEAST_EFFECTIVE_LIB "${MYSTRAN_FEAST_LIB}")
  endif()

  add_library(mystran_feast STATIC IMPORTED GLOBAL)
  set_target_properties(mystran_feast PROPERTIES IMPORTED_LOCATION "${MYSTRAN_FEAST_EFFECTIVE_LIB}")
  add_dependencies(mystran_feast mystran_feast_build)

  target_compile_definitions(mystran PRIVATE MYSTRAN_HAVE_EXTERNAL_FEAST)
  target_include_directories(mystran PRIVATE "${MYSTRAN_FEAST_REAL_ROOT}/include")
  target_link_libraries(mystran mystran_feast)
endif()
# ! --- FEAST V4.0 external hook --- end !

# ! --- MUMPS_COO add begin --- !
if(MYSTRAN_USE_DMUMPS_SOLVER)
  target_compile_definitions(mystran PRIVATE DMUMPS_Solver)
  if(NOT EXISTS "${MYSTRAN_DMUMPS_INCLUDE_DIR}/dmumps_struc.h"
     OR NOT EXISTS "${MYSTRAN_DMUMPS_INCLUDE_DIR}/mpif.h")
    message(FATAL_ERROR
      "MYSTRAN_USE_DMUMPS_SOLVER=ON but dmumps_struc.h or mpif.h was not found in "
      "${MYSTRAN_DMUMPS_INCLUDE_DIR}")
  endif()
  target_include_directories(mystran PRIVATE "${MYSTRAN_DMUMPS_INCLUDE_DIR}")

  if(NOT MYSTRAN_DMUMPS_EXTRA_LIBS)
    message(FATAL_ERROR
      "MYSTRAN_USE_DMUMPS_SOLVER=ON requires MYSTRAN_DMUMPS_EXTRA_LIBS.")
  endif()
  foreach(MYSTRAN_DMUMPS_LIB IN LISTS MYSTRAN_DMUMPS_EXTRA_LIBS)
    if(NOT EXISTS "${MYSTRAN_DMUMPS_LIB}")
      message(FATAL_ERROR "Required static MUMPS library was not found: ${MYSTRAN_DMUMPS_LIB}")
    endif()
  endforeach()
  target_link_libraries(mystran ${MYSTRAN_DMUMPS_EXTRA_LIBS})

  if(WIN32)
    target_link_libraries(mystran pthread)
  endif()
endif()
# ! --- MUMPS_COO add end --- !

# determine which superlu variant to link against
if(USE_SUPERLU_MT)
  if(PLAT STREQUAL "_PTHREAD")
    message(WARNING "We recommend using OpenMP, not pthread!")
    if(MYSTRAN_USE_OPENBLAS)
      target_link_libraries(mystran superlu_mt_PTHREAD ${TPL_BLAS_LIBRARIES})
    else()
      target_link_libraries(mystran superlu_mt_PTHREAD f2c)
    endif()
  elseif(PLAT STREQUAL "_OPENMP")
    if(MYSTRAN_USE_OPENBLAS)
      target_link_libraries(mystran superlu_mt_OPENMP ${TPL_BLAS_LIBRARIES})
    else()
      target_link_libraries(mystran superlu_mt_OPENMP f2c)
    endif()
  endif()

  target_compile_definitions(mystran PRIVATE USE_SUPERLU_MT)
else()
  if(MYSTRAN_USE_OPENBLAS)
    target_link_libraries(mystran superlu ${TPL_BLAS_LIBRARIES})
  else()
    target_link_libraries(mystran superlu f2c)
  endif()
endif()

set_target_properties(
  mystran PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
)

if(WIN32)
  target_sources(mystran PRIVATE mystran.rc)
endif()

if(CMAKE_COMPILER_IS_GNUCC)
  set(MYSTRAN_GNU_WARN_FLAGS "-Wall -Wextra \
    -Wno-unused-variable -Wno-unused-label -Wno-unused-parameter -Wno-tabs \
    -Wno-compare-reals -Wno-character-truncation -Wno-unused-dummy-argument \
    -Wmaybe-uninitialized -Wrealloc-lhs \
    -ffree-line-length-none -ffixed-line-length-none")
  if(CMAKE_BUILD_TYPE STREQUAL "AVX2")
    set(MYSTRAN_AVX2_FLAGS "-O3 -funroll-loops -march=core-avx2 -mtune=core-avx2 -fno-fast-math -ffp-contract=off")
    set(CMAKE_Fortran_FLAGS "${MYSTRAN_GNU_WARN_FLAGS} ${MYSTRAN_AVX2_FLAGS}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MYSTRAN_AVX2_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MYSTRAN_AVX2_FLAGS}")
    foreach(AVX2_TARGET mystran f2c superlu blas)
      if(TARGET ${AVX2_TARGET})
        target_compile_options(${AVX2_TARGET} PRIVATE -O3 -funroll-loops -march=core-avx2 -mtune=core-avx2 -fno-fast-math -ffp-contract=off)
      endif()
    endforeach()
  else()
    set(CMAKE_Fortran_FLAGS "${MYSTRAN_GNU_WARN_FLAGS} -fbounds-check -fcheck=all")
  endif()
  set(CMAKE_Fortran_FLAGS_DEBUG
    "-g -O0 -fcheck=all -fbacktrace -fbounds-check \
    -fno-inline -fno-ipa-sra -fno-ipa-cp -fno-optimize-sibling-calls")
  set(CMAKE_Fortran_FLAGS_DETERMINISTIC
    "-O0 -g \
     -fno-fast-math \
     -ffp-contract=off \
     -fno-unsafe-math-optimizations \
     -fno-associative-math \
     -fno-reciprocal-math \
     -frounding-math"
  )
  set(CMAKE_C_FLAGS_DETERMINISTIC
    "-O0 -g \
     -fno-fast-math \
     -ffp-contract=off \
     -frounding-math"
  )
  set(CMAKE_CXX_FLAGS_DETERMINISTIC
    "${CMAKE_C_FLAGS_DETERMINISTIC}"
  )
  set(CMAKE_PROFILING_FLAGS "-O2 -g -fno-inline -fno-omit-frame-pointer -fno-inline-functions")
  set(CMAKE_C_FLAGS_PROFILING "${CMAKE_PROFILING_FLAGS}")
  set(CMAKE_CXX_FLAGS_PROFILING "${CMAKE_PROFILING_FLAGS}")
  set(CMAKE_Fortran_FLAGS_PROFILING "${CMAKE_PROFILING_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS_PROFILING "")
  set(CMAKE_SHARED_LINKER_FLAGS_PROFILING "")
endif()

# issue a couple messages about compilation
include(ProcessorCount)
ProcessorCount(NCPU)
message(STATUS "You can now compile MYSTRAN with `cmake --build .`")

if(NOT NCPU EQUAL 0)
  message(STATUS "Compile faster by passing the -j${NCPU} flag.")
endif()

