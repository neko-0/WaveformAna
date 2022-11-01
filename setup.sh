#!/bin/bash

export TestArea=$(pwd)

# alias build='cmake --build $TestArea/../build; source $TestArea/../build/*/setup.sh'
alias build='cmake --build $TestArea/../build'

cd $TestArea/../build

cmake $TestArea

# there must be a easier way to set alias for executable
alias run_Ana='$TestArea/../build/analysisDriver/run_Ana'
