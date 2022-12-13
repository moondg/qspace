#!/bin/bash

# minimal centralized startup script to setup Matlab environment
# using BASH_SOURCE[0] rather than $0 so this also works for `source'

  if command -v realpath >/dev/null; then
       P0=`dirname $(realpath "${BASH_SOURCE[0]}")`;  # linux
  else P0=`dirname $(readlink "${BASH_SOURCE[0]}")`;  # macOS
  fi

# -------------------------------------------------------------------- #
# NB! matlab requires very specific gcc version for mex/mcc compilation

  if [ -z $MATLAB_ROOT ]; then
   # this assumes that the module system is set up
   # otherwise ensure that `matlab' and `gcc' can be found on the PATH
   # or that MATLAB_ROOT and GCC are set (the latter is required
   # for mex-compilation only)
     if ! `type module >/dev/null 2>&1`; then
        printf "\n  ERR $0 : module environment not available\n";
        exit 1
     fi

   # module load matlab/2016a # NB! matlab/R2016 pairs with gcc/4.7.4
   # module load gcc/4.7.4

   # module load matlab/2018b # NB! matlab/R2018 pairs with gcc/6.3.x
   # module load gcc/6.3.0

     module load matlab/2020b
     module load gcc/8.4.0

   # on macOS
   # matlab/2020b -> xcode/11.x
   # matlab/2022b -> xcode/13.x
  fi

  if [ -z "$MATLAB_ROOT" ]; then
     MATLAB_ROOT=`which matlab`;
     export MATLAB_ROOT=${MATLAB_ROOT//bin*matlab*/}
  fi

# -------------------------------------------------------------------- #
# directory setting (please adapt)

# your matlab root directory in HOME
# script 'ml' starts in here -> ensure that $MYMATLAB/startup.m is present
# export MYMATLAB=$HOME/Matlab
  export MYMATLAB="${P0//system/}"
  export MSLOTS=12

# LMA = local matlab data directory
  if [ ! $LMA ]; then
     export LMA=/data/$USER/Data
  fi
  if [ ! $RC_STORE ]; then
     export RC_STORE=$LMA/RCStore  # for non-abelian symmetries
     export RC_SYNC=$LMA/RCSync    # just for lock files
  fi
  
# this overwrites the defaults chosen by matlab_setup.pl below
# export QSP_NUM_THREADS=2  #> QSpace
# export OMP_NUM_THREADS=4  #> Intel MKL library
# export MKL_NUM_THREADS=4  #> Intel MKL library

# -------------------------------------------------------------------- #
# double check and complete parameter setting

  if [[ "$*" == '-t' ]]; then
       $P0/matlab_setup.pl -t  # use this to check/test output
  else eval "$($P0/matlab_setup.pl "$@")"
  fi

# -------------------------------------------------------------------- #

  unset P0

