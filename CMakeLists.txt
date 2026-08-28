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

add_library(folded_xxz_physics
  src/initial_state.cc
  src/observables.cc
  src/three_site_hamiltonian.cc
  src/trotter_evolution.cc
)

add_library(folded_xxz::physics ALIAS folded_xxz_physics)

target_compile_features(folded_xxz_physics PUBLIC cxx_std_17)
target_include_directories(folded_xxz_physics PUBLIC
  "${CMAKE_CURRENT_SOURCE_DIR}/include"
  "${ITENSOR_ROOT}"
)
target_link_libraries(folded_xxz_physics PUBLIC
  "${ITENSOR_ROOT}/lib/libitensor.a"
  Threads::Threads
  BLAS::BLAS
  LAPACK::LAPACK
)

add_library(folded_xxz_runtime
  src/output.cc
  src/parameters.cc
  src/simulation.cc
  src/simulation_config.cc
  src/simulation_runner.cc
)
add_library(folded_xxz::runtime ALIAS folded_xxz_runtime)
target_compile_features(folded_xxz_runtime PUBLIC cxx_std_17)
target_link_libraries(folded_xxz_runtime PUBLIC folded_xxz::physics)

add_executable(3siteHam src/main.cc)
target_link_libraries(3siteHam PRIVATE folded_xxz::runtime)

include(CTest)

if(BUILD_TESTING)
  add_executable(observables_test tests/observables_test.cc)
  target_link_libraries(observables_test PRIVATE folded_xxz::physics)
  add_test(NAME observables_test COMMAND observables_test)

  add_executable(simulation_schedule_test tests/simulation_schedule_test.cc)
  target_compile_features(simulation_schedule_test PRIVATE cxx_std_17)
  target_include_directories(simulation_schedule_test PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")
  add_test(NAME simulation_schedule_test COMMAND simulation_schedule_test)

  add_executable(trotter_evolution_test tests/trotter_evolution_test.cc)
  target_link_libraries(trotter_evolution_test PRIVATE folded_xxz::physics)
  add_test(NAME trotter_evolution_test COMMAND trotter_evolution_test)

  add_executable(initial_state_test tests/initial_state_test.cc)
  target_link_libraries(initial_state_test PRIVATE folded_xxz::physics)
  add_test(NAME initial_state_test COMMAND initial_state_test)

  add_executable(simulation_config_test tests/simulation_config_test.cc)
  target_link_libraries(simulation_config_test PRIVATE folded_xxz::runtime)
  add_test(NAME simulation_config_test COMMAND simulation_config_test)
endif()

find_program(CLANG_FORMAT_EXECUTABLE clang-format)

if(CLANG_FORMAT_EXECUTABLE)
  file(GLOB_RECURSE FORMAT_SOURCES CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.h"
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
