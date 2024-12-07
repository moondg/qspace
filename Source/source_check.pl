#!/usr/bin/perl -w
# Usage: eval $(source_check.pl [tag] ...)
#
#    check for specific tags set / defined within given C source file.
#    NB! exit 0 in case of success (tag found), and non-zero otherwise
#    This then allows simple usage with bash scripts, as in
#    if NAME [-tag] source.cc; then .. else .. fi;
#
# Options #1: check wether LD_CLEBSCH_QS is defined in any source file
#    --cgc mexfun1.cc ...
#
# Options #2: find location of binary of each of the source files
#    --loc [default-output-directory] mexfun1.cc ...
#
# Wb,Oct10,14 ; Wb,Oct03,18

  my $me=$0; $me=~s/.*\///;

  if (!@ARGV) { die "\n  ERR $me: invalid usage\n"; }

  my $tag=shift; my ($f,@ll,@lx,$pat,@opts,@ff,@dd);

  foreach (@ARGV) {
     if    (-f $_) { push(@ff,$_); }
     elsif (-d $_) { push(@dd,$_); } else { push(@opts,$_); }
  }

  if ($tag eq '--cgc') {
     $pat='LD_CLEBSCH_QS'; # looking for "#define LD_CLEBSCH_QS"
     my $q=0;

     if (!@ff && @opts) { die "\n  ERR $me: ".
        "invalid file(s) with $tag: '".join("' '",@opts)."'\n";
     }
     if (@dd || @opts) { die "\n  ERR $me: invalid usage with $tag\n"; }
     foreach (@ff) { $f=$_;
        if (!/[^\.]*\.[chm].*$/) { die "\n  ERR invalid C-file $_\n\n"; }

        open(FH,'<',$f); 
        foreach (<FH>) {
           if (/^\s*#(define|undef)\s+(LD_CLEBSCH_QS|QS_SKIP_MPFR)/) {
              my $a=$1; my $b=$2; $b=($b=~/LD/ ? 1 : 2);
              if ($a eq 'define')
                   { $q |=  $b; }
              else { $q &= ~$b; }
           }
        }; close(FH);
        if ($q==1) { exit 0; }
     }

     exit 1;
  }
  elsif ($tag eq '--loc') {
     if (!@dd) {
        foreach ('./bin','../bin') {
        if (-d $_) { push(@dd,$_); last; }}
     }
     if (!@dd) { die
        "\n  ERR $me: invalid usage (missing default directory)\n"; }
     if (@dd>1) { die "\n  ERR $me: invalid usage with $tag\n"; }

     foreach $f (@ff,@opts) { my (%D,@d1); $_=$f; s/\.[^\.]+$//;
        @ll=`find .. -name '*$_*'`;
        @lx=grep(/\.mex/,@ll);
        foreach my $m (@lx) { $_=$m; s/\/[^\/]+$//; $D{$_}++; }
        @d1=sort(keys %D);

        if (@d1>1) { @d1=grep(!/Archive/i,@d1); }
        if (@d1>1) { 
           printf STDERR
           "\n  WRN got mex-files in multiple directories for $f:\n";
              for (my $i=0; $i<@lx; ++$i) { $_=$lx[$i]; chomp;
               printf "  %5d)  %s\n",$i+1,$_; }
           print STDERR "\n";

           $pat=`pwd`; $pat=~s/[\/\s]*$//; $pat=~s/.*\///; $pat="\b$pat\b";
           foreach (@d1) { if (/$pat/) { @d1=grep(/$pat/,@d1); last; }}
        }

        $_=(@d1 ? $d1[0] : $dd[0]);
        print $_;
     }
     exit 0;
  }
  else { die "\n  ERR $me: invalid usage (invalid tag=$tag)\n"; }

