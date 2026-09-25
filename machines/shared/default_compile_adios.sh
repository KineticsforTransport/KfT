#!/bin/bash

set -e

SYSTEM_MPI=$(../../bin/julia --project ../shared/get_mk_preference.jl use_system_mpi)
USE_ADIOS=$(../../bin/julia --project ../shared/get_mk_preference.jl use_adios)
BUILD_ADIOS=$(../../bin/julia --project ../shared/get_mk_preference.jl build_adios)

if [[ $BUILD_ADIOS == "y" && -d adios-build ]]; then
  echo "ADIOS2 appears to have been downloaded, compiled and installed already."
  echo "Do you want to download, compile and install again, overwriting the existing "
  echo "version? y/[n]"
  read -p "> " input
  while [[ ! -z $input && !( $input == "y" || $input == "n" ) ]]; do
    echo
    echo "$input is not a valid response: y/[n]"
    read -p "> " input
  done
  if [[ -z $input || $input == "n" ]]; then
    BUILD_ADIOS="n"
  else
    # Remove the install directory if it exists already
    if [[ -d adios-build ]]; then
      rm -r adios-build
    fi
  fi
fi

if [[ $BUILD_ADIOS == "y" ]]; then
  # Clone the source code.
  # The version probably needs to match the one that would be provided by
  # ADIOS2_jll, as this will be compatible with ADIOS2.jl
  git clone https://github.com/ornladios/ADIOS2.git -b v2.12.1 --depth=1

  cd ADIOS2/

  # Build using CMake.
  # Exclude some optional features that sometimes cause linking errors when
  # system versions of libraries conflict with the Julia package manager's
  # versions.
  cmake -B build -DCMAKE_INSTALL_PREFIX=$ARTIFACT_DIR/adios-build -DADIOS2_USE_SZ=OFF -DADIOS2_USE_UCX=OFF | tee cmake_config.log
  cmake --build build --parallel 8 | tee cmake_build.log
  cmake --install build | tee cmake_install.log
fi

exit 0
