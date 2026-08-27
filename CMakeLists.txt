cmake_minimum_required(VERSION 3.21)

project(FoldedXXZ LANGUAGES CXX)

set(ITENSOR_ROOT "" CACHE PATH
    "Path to the root of the built ITensor C++ source tree")

if(NOT ITENSOR_ROOT)
  message(FATAL_ERROR
    "ITENSOR_ROOT is required. Configure with "
    "-DITENSOR_ROOT=/path/to/itensor")
endif()

find_path(ITENSOR_INCLUDE_DIR
  NAMES itensor/all.h
  HINTS "${ITENSOR_ROOT}"
  NO_DEFAULT_PATH
  REQUIRED
)

find_library(ITENSOR_LIBRARY
  NAMES itensor
  HINTS "${ITENSOR_ROOT}/lib"
  NO_DEFAULT_PATH
  REQUIRED
)

find_package(Threads REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# ITensor C++ v3 supplies a static library through its Makefile build.
add_library(ITensor::itensor STATIC IMPORTED GLOBAL)
set_target_properties(ITensor::itensor PROPERTIES
  IMPORTED_LOCATION "${ITENSOR_LIBRARY}"
  INTERFACE_INCLUDE_DIRECTORIES "${ITENSOR_INCLUDE_DIR}"
  INTERFACE_LINK_LIBRARIES "Threads::Threads;BLAS::BLAS;LAPACK::LAPACK"
)

add_executable(3siteHam
  3siteHam.cc
  observables.cc
)

target_compile_features(3siteHam PRIVATE cxx_std_17)
target_link_libraries(3siteHam PRIVATE ITensor::itensor)
