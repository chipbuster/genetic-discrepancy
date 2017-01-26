cmake_minimum_required(VERSION 3.3.0)
project (genetic-star-dirscrepancy)

include(ExternalProject)

# Set output directories so that the build directory is NEAT UND TIDY
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Set OpenMP
FIND_PACKAGE( OpenMP REQUIRED)
if(OPENMP_FOUND)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif(OPENMP_FOUND)

# Require C++11 to be used
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANARD_REQUIRED ON)

# Use new IN_LIST behavior
cmake_policy(SET CMP0057 NEW)

################################
## CHECK COMPILER INFORMATION ##
################################

# Require that icc is used for release--important for AVX512 support!
# This may be relaxed in the future if gcc supports Knight's Landing
# Relax this requirement for non-release builds.

if("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
  if(NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    message(FATAL_ERROR "In Release mode, the compiler for this project must be icc.\
    Please either load the appropriate module or manually specify CMAKE_CXX_COMPILER\
    on the command line.")
  endif(NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")

  # Require icc > 17.0.0 for KNL support

  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 17.0.0)
    message(ERROR "Your compiler version older than what is supported by this\
    project. We will attempt to continue, but bizzare things may happen.\
    For best results, please use an icc newer than version 17.0.0. You are\
    currently using ${CMAKE_CXX_COMPILER_VERSION}")
  endif(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 17.0.0)
endif("${CMAKE_BUILD_TYPE}" STREQUAL "Release")

########################
## BUILD BLIS PROJECT ##
########################

option(BUILD_BLIS "Build the BLIS project" ON)

if(BUILD_BLIS)

  set(BLIS_SOURCE_DIR "${PROJECT_SOURCE_DIR}/blis")
  set(BLIS_PREFIX "${CMAKE_BINARY_DIR}/thirdparty")

  # Set/check architectures. We need to be on some x86 system, so limit valid values
  set(BLIS_VALID_ARCHITECTURES "auto" "haswell" "knl" "mic" "piledriver" "bulldozer" "sandybridge" "reference")
  if(NOT BLIS_ARCH)
    set(BLIS_ARCH "auto") #Set a default value
  endif(NOT BLIS_ARCH)
  if(NOT BLIS_ARCH IN_LIST BLIS_VALID_ARCHITECTURES)
    message(SEND_ERROR "Your specified BLIS_ARCH is not valid. Please choose a valid BLIS_ARCH.

  Note that this project is restricted to x86 architectures, so not all valid BLIS architectures
  are valid here. Your choices for valid architectures are:
    - auto
    - haswell
    - knl
    - mic
    - piledriver
    - bulldozer
    - sandybridge
    - reference

  I am resetting your architecture to \"auto\".")
    set(BLIS_ARCH "auto")
  endif(NOT BLIS_ARCH IN_LIST BLIS_VALID_ARCHITECTURES)

  ExternalProject_Add(BLIS
    ## Configure phase commands
    SOURCE_DIR "${BLIS_SOURCE_DIR}"
    CONFIGURE_COMMAND "${BLIS_SOURCE_DIR}/configure" "--prefix=${BLIS_PREFIX}"
    --disable-cblas --disable-blas
    --enable-verbose-make
    --enable-shared --enable-static
    ${BLIS_ARCH}

    ## Build phase commands
    BUILD_COMMAND $(MAKE)      # This is not portable...oh whale

    ## Install phase commands
    INSTALL_COMMAND $(MAKE) install

    TEST_COMMAND make test

    LOG_BUILD 1
    )

endif(BUILD_BLIS)
####################
## BUILD OUR CODE ##
####################

include_directories("${PROJECT_SOURCE_DIR}/include")

### Build base code
add_executable(stardiscrep-noaccel.x "${PROJECT_SOURCE_DIR}/base/star-discrepancy.cpp")
add_executable(stardiscrep-noaccel-refact.x "${PROJECT_SOURCE_DIR}/base/star-discrepancy-refact.cpp")

####################
## RUN SOME TESTS ##
####################

enable_testing()

add_test(base_correctness "${PROJECT_SOURCE_DIR}/../tests/simple_correctness/run_simple_correctness_test.sh"
  "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/stardiscrep-noaccel-refact.x"
  "${PROJECT_SOURCE_DIR}/../tests/simple_correctness/input_data")
