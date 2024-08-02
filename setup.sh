#!/bin/bash

export CURRENT_PATH=$(pwd)
export SRC_AREA=$(dirname $BASH_SOURCE)
export WAVEANA_BUILD_DIR=$SRC_AREA/../WaveformAna_build

echo "Found source directory: ${SRC_AREA}"
echo "Setting up build directory: ${WAVEANA_BUILD_DIR}"

mkdir -p $WAVEANA_BUILD_DIR
cd $WAVEANA_BUILD_DIR
cmake -B $WAVEANA_BUILD_DIR -S $SRC_AREA -DCMAKE_BUILD_TYPE=RELEASE
cd $CURRENT_PATH

alias build='cmake --build $WAVEANA_BUILD_DIR'
# there must be a easier way to set alias for executable
alias run_Ana='$WAVEANA_BUILD_DIR/analysisDriver/run_Ana'
alias ls_Ana='$WAVEANA_BUILD_DIR/analysisDriver/ls_Ana'
alias check_Ana='$WAVEANA_BUILD_DIR/analysisDriver/check_Ana'
