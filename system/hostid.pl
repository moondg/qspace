#!/usr/bin/perl -w
# Usage: NAME
#
#    get few-letter host-id
#    e.g. used for log files (see mbatch).
#
# Wb,Feb24,10

  require "$ENV{HOME}/bin/plib.pl";
  my $pflag=0; my $jflag=0;

  while (@ARGV) { $_=$ARGV[0];
     if (/^-[h\?]$/) { die usage($0); }
     elsif ($_ eq '-pid') { $pflag=1; shift; }
     elsif ($_ eq '-jid') { $jflag=1; shift; }
     else { last; }
  }

  if (@ARGV) {
     print STDERR "\n  WRN input arguments will be ignored! (".
     join(', ',@ARGV).")\n";
  }

  my $hid=hostid();

  if ($jflag) { my $j=jobid();
     if ($hid=~/\d$/) { $hid.=('_'.$j); }
     else { $hid.=$j; }
  }
  elsif ($pflag) { $hid.='.'.getppid(); }

  print $hid;

