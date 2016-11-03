cmake_minimum_required(VERSION 3.6.0)
project (genetic-star-dirscrepancy)

include(ExternalProject)

# Set output directories so that the build directory is NEAT UND TIDY
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

########################
## BUILD BLIS PROJECT ##
########################

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
  BUILD_COMMAND make

  ## Install phase commands
  INSTALL_COMMAND make install

  TEST_COMMAND make test
)

####################
## BUILD OUR CODE ##
####################

#include_directories("${PROJECT_SOURCE_DIR}/include")

### Build base code
#add_executable("${PROJECT_SOURCE_DIR}/base/star-discrepancy.cpp")
