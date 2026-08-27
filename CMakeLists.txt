cmake_minimum_required(VERSION 3.21)

project(FoldedXXZ LANGUAGES CXX)

set(ITENSOR_ROOT "" CACHE PATH
    "Path to the root of the built ITensor C++ source tree")

if(NOT EXISTS "${ITENSOR_ROOT}/itensor/all.h" OR
   NOT EXISTS "${ITENSOR_ROOT}/lib/libitensor.a")
  message(FATAL_ERROR
    "ITENSOR_ROOT must contain itensor/all.h and lib/libitensor.a. "
    "Configure with -DITENSOR_ROOT=/path/to/itensor")
endif()

find_package(Threads REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

add_library(folded_xxz_observables
  src/observables.cc
)

add_library(folded_xxz::observables ALIAS folded_xxz_observables)

target_compile_features(folded_xxz_observables PUBLIC cxx_std_17)
target_include_directories(folded_xxz_observables PUBLIC
  "${CMAKE_CURRENT_SOURCE_DIR}/include"
  "${ITENSOR_ROOT}"
)
target_link_libraries(folded_xxz_observables PUBLIC
  "${ITENSOR_ROOT}/lib/libitensor.a"
  Threads::Threads
  BLAS::BLAS
  LAPACK::LAPACK
)

add_executable(3siteHam src/main.cc)
target_link_libraries(3siteHam PRIVATE folded_xxz::observables)
