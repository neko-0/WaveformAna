#!/bin/bash

export TestArea=$(pwd)
export WAVEANA_BUILD_DIR=$TestArea/../WaveformAna_build

mkdir -p $WAVEANA_BUILD_DIR
cd $WAVEANA_BUILD_DIR
cmake -B $WAVEANA_BUILD_DIR -S $TestArea

alias build='cmake --build $WAVEANA_BUILD_DIR'
# there must be a easier way to set alias for executable
alias run_Ana='$WAVEANA_BUILD_DIR/analysisDriver/run_Ana'
alias ls_Ana='$WAVEANA_BUILD_DIR/analysisDriver/ls_Ana'
alias check_Ana='$WAVEANA_BUILD_DIR/analysisDriver/check_Ana'
