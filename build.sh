#!/bin/env bash
n_build_tasks=8

echo "Building release"
build_dir=builds/release
nek_dir=
rm -rf $build_dir
cmake -B $build_dir -DNektar++_DIR="$nek_dir"
cmake --build $build_dir -j $n_build_tasks

echo "Building debug"
debug_build_dir=builds/debug
debug_nek_dir=
rm -rf $debug_build_dir
cmake -B $debug_build_dir -DCMAKE_BUILD_TYPE=DEBUG -DNektar++_DIR="$debug_nek_dir"
cmake --build $debug_build_dir -j $n_build_tasks
echo "Done"