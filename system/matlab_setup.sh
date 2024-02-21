#!/bin/bash

# -------------------------------------------------------------------- #
# Minimal centralized startup script to setup Matlab environment.      #
# -------------------------------------------------------------------- #
# Please no longer modify this file, since it is tracked in the        #
# git repository. Rather use system/matlab_setup_user.sh               #
# or set the suggested environmental variable QS_CONFIG_ML_SH          #
# pointing to your configuration file, instead. // Wb,Feb15,24         #
# A template is provided in system/matlab_setup_user.sh-template       #
# -------------------------------------------------------------------- #

  if command -v realpath >/dev/null; then
     # using BASH_SOURCE[0] rather than $0 since $0 is not set
     # when using `source *this'
       P0=`realpath    "${BASH_SOURCE[0]}"`;  # linux
  else P0=`readlink -f "${BASH_SOURCE[0]}"`;  # macOS
  fi

  P0=`dirname "$P0"`

# script 'ml' starts form MYMATLAB
# to ensure that $MYMATLAB/startup.m is called.
# You may change the default here based on the location
# of this script in matlab_setup_user.sh called below
  export MYMATLAB="${P0//system/}"

# You need to put your own setup into a separate file that is
# not part of this git repository, so it does not get overwritten
# by git updates; you may do this by defining an environmental
# variable QS_CONFIG_ML_SH that points to script; otherwise,
# by default, this looks for matlab_setup_user.sh in the same
# directory as this file // Wb,Feb14,24
  if [ -z "$QS_CONFIG_ML_SH" ]; then
     msh="$P0/matlab_setup_user.sh"
     if [ -f "$msh" ]; then
        export QS_CONFIG_ML_SH="$msh"
   # else
     # issue error in matlab_setup.pl script below
     # printf "\n \e[31m matlab_setup.sh:\n  %s\n  %s\n  %s\e[0m\n\n" \
     #  "Please setup the system environment for QSpace first" \
     #  "e.g., set environmental variable QS_CONFIG_ML_SH or use file" "$msh"
     fi
  fi

  if [ ! -z "$QS_CONFIG_ML_SH" ]; then
     source $QS_CONFIG_ML_SH
  fi

# double check and complete parameter setting
  msh="$P0/matlab_setup.pl"
  if [[ "$*" == '-t' ]]; then
   # use "matlab_setup.sh -t" to check/test printed output
     "$msh" -t 
  else
     sout="$($msh "$@")"
     if [ $? -eq 0 ]; then
          eval "$sout"
     else echo "$sout"; fi
  fi

  unset P0 msh sout

# -------------------------------------------------------------------- #

