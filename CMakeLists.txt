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

add_executable(3siteHam
  src/main.cc
  src/parameters.cc
  src/simulation.cc
)
target_link_libraries(3siteHam PRIVATE folded_xxz::observables)

include(CTest)

if(BUILD_TESTING)
  add_executable(observables_test tests/observables_test.cc)
  target_link_libraries(observables_test PRIVATE folded_xxz::observables)
  add_test(NAME observables_test COMMAND observables_test)

  add_executable(simulation_schedule_test tests/simulation_schedule_test.cc)
  target_link_libraries(simulation_schedule_test PRIVATE folded_xxz::observables)
  add_test(NAME simulation_schedule_test COMMAND simulation_schedule_test)
endif()

find_program(CLANG_FORMAT_EXECUTABLE clang-format)

if(CLANG_FORMAT_EXECUTABLE)
  file(GLOB_RECURSE FORMAT_SOURCES CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc"
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cc"
  )

  add_custom_target(format
    COMMAND "${CLANG_FORMAT_EXECUTABLE}" -i ${FORMAT_SOURCES}
    COMMENT "Formatting C++ sources"
  )

  add_custom_target(format-check
    COMMAND "${CLANG_FORMAT_EXECUTABLE}" --dry-run --Werror ${FORMAT_SOURCES}
    COMMENT "Checking C++ formatting"
  )
endif()
