#!/bin/bash

# usage: rnrg.sh <matfile>

# export OMP_NUM_THREADS=$[NCPU-2] (set in .matlab_setup)
. ~/.matlab_setup -ld -f
# unset OMP_NUM_THREADS

  NRG=NRGWilsonCG 
  FDM=fdmNRG_CG
  TAG=FDM-NRG

  perl -we '
     if ($#ARGV<0) { exit 1; }
     if ($#ARGV>0) { exit 2; }
     $_=shift(@ARGV); if (! -f $_ || !/\.mat$/) { exit 2; }
     exit 0;
  ' -- "$@" ; r=$?

  if [ $r -eq 1 ]; then
     printf "\n  usage: rnrg.sh <matfile>\n"; exit $r
  elif [ $r -eq 2 ]; then
     $MEX/$NRG "$@" ; exit $r
  fi

  MAT=$1; shift 1

# echo [$r] run ... ; exit 1

  ODIR=$LMA/../../bin/NRG_$$/
  if [ -d "$ODIR" ]; then
     printf "\n  ERR output directory already exists!\n  $ODIR\n"
     exit $?
  fi

  mkdir $ODIR
  if [ $? -gt 0 ]; then
     printf "\n  ERR failed to setup output directory!\n  $ODIR\n"
     exit $?
  fi

  cp -pf $MEX/$NRG $MEX/$FDM $MAT $ODIR
  cd $ODIR ; mkdir NRG

  LOG=Wb$(date +"%y%m%d")_${TAG}_$$.log
  MAT=$(/bin/ls *.mat)

  print_status() {
     printf "\n$(hostname)//$(pwd):\n\n"
     ll --color=none ; printf "\n"
     du -h NRG/ ; printf "\n"
     df -h | egrep 'Avail|data|/(project|home|tmp)' | egrep -v filer
  }

  printf "$(date '+%d-%b-%Y %R:%S')  %-10s %s  $MAT\n" \
  $HOST $(pwdprompt $DIR) >> $MLBLOG

{ EMAIL=$USER@theorie.physik.uni-muenchen.de
  { printf "\n  DIR: $(pwd)\n" ; /bin/ls -ort
  } 2>&1 | mail -s "Job $TAG started on $(hostname)" $EMAIL

  { set | egrep 'OMP_' ; printf "\n"; # OMP_NUM_THREADS

    printf "   started: $(date)\n"
    printf "  location: $(hostname)//$(pwd) :: $NRG\n"
  ####################################################
    nice -19 ./$NRG $MAT ; r=$?                      #
  ####################################################
    printf "\nfinished: $(date)\n" ; print_status
    if [ ! $r -eq 0 ]; then exit $r; fi

    printf "   started: $(date)\n"
    printf "  location: $(hostname)//$(pwd) :: $FDM\n"
  ####################################################
    nice -19 ./$FDM NRG/NRG $MAT op1 op2 ; r=$?      #
  ####################################################
    printf "\nfinished: $(date)\n\n" ; print_status
    if [ ! $r -eq 0 ]; then exit $r; fi

    mv -v ./NRG/*info* .  # 2>/dev/null
    rm -rf ./NRG/ # may still be required for subsequent FDM run!

  } > $LOG 2>&1

  { print_status
  } 2>&1 | mail -s "Job $TAG finished on $(hostname) [$r]" $EMAIL
} &

  printf "\n  $NRG + $FDM started ...\n"
  printf "  dir: $(pwd)\n  log: $LOG\n  mat: $MAT\n"


