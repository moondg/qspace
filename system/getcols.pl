#!/usr/bin/perl -w
# Usage: NAME
# Wb,Nov11,11

  my $n=90; my $vflag=1;
  if (@ARGV && $ARGV[0] eq '-q') { $vflag=0; }

  my @ll=`resize 2>/dev/null`; # if ($?) { die "\n  ERR resize returned $? !?\n"; }
  foreach (@ll) {
     if (/COL.*=(\d+)/) { if ($1>30 && $1<$n) { $n=$1; }}
  }

  if ($vflag) { print $n; }
  exit $n;

