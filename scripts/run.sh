#!/usr/bin/env bash

set -euo pipefail

project_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
build_dir=${BUILD_DIR:-"${project_dir}/build"}
itensor_root=${ITENSOR_ROOT:-"${HOME}/Programming/itensor"}
build_type=${CMAKE_BUILD_TYPE:-Release}

cmake -S "${project_dir}" -B "${build_dir}" \
  -DITENSOR_ROOT="${itensor_root}" \
  -DCMAKE_BUILD_TYPE="${build_type}"

cmake --build "${build_dir}" --parallel
cmake --build "${build_dir}" --target format-check
ctest --test-dir "${build_dir}" --output-on-failure

exec "${build_dir}/3siteHam" "$@"
