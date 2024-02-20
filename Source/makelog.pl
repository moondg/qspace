#!/usr/bin/perl -w
# Usage: [make ... ] | NAME
#
#    This acts like "| head -20" to mex compiler output
#    while also color coding warnings and errors.
#
#    The default of 20 lines can be changed by setting
#    the environmental variable NMAX_MAKELOG.
#
# Wb,Apr22,10 ; Wb,Feb18,24

  use strict; use warnings;

  my ($P,$me); $P=$0; $P=~s/\/?([^\/]*)$//; $me=$1;
  if ($P eq '.') { $P=".."; } else { $P.="/.."; }
  require "$P/system/plib.pl";

  my ($lflag); my $tflag=0; $_=$ENV{LOG_FLAG};
  my $vflag=1; if ($_ && $_ eq '-q') { $vflag=0; }

  foreach (@ARGV) {
     if (/^-[h\?]$/) { die usage($0); }
     elsif (/^-(v)$/i) { $vflag+=($1 eq 'v' ? 1 : 2); }
     elsif (/^-(t)$/i) { $tflag+=($1 eq 't' ? 1 : 2); }
     elsif (/^--test$/i) { $tflag|=16; }
     else { die usage(__FILE__,__LINE__,"invalid usage:",@ARGV); }
  }

  my ($n,$n1,$lcont,@l2,$L); my $ncols='80'; my $e=0;
  my $H=$ENV{HOME};
  my $N=$ENV{NMAX_MAKELOG}; if (!$N) { $N=20; }
  my $MML=$ENV{MYMATLAB};

  if (-t STDIN) { die usage(); }

  if ($tflag<8)
       { parse_makelog_regular(); }
  else { parse_makelog_test(); }

sub parse_makelog_regular {

  foreach(<STDIN>) {
     if (/lib(mex|mat|mx|mwm|mws).*\.so: undefined reference/) { ++$n1; next; }
     if (/^\/usr\/bin\/ld: warning: lib[\w]+.so[.\d]*, needed by.*try using/) {
        ++$n1; next;
     }
     if (/^\/usr\/local.*matlab.*libmwlapack.*undefined reference/) {
        ++$n1; next;
     }

     if (s/Building with '(.*)'/built with $1/){ chomp; s/^\s*//; 
        print "\r  ... $_ ... \r";  next;
     }
     if (/MEX completed successfully/) { next; }

     if (++$n>$N && $e>1) { push(@l2,$_); next; };

     if (++$lflag==1) {
         $L=('=' x $ncols)."\n"; print "\n\n".$L;
         $L=('-' x $ncols)."\n";
     }

     if ($MML) { s/$MML/\$MLL/g; }
     if ($H)   { s/$H/~/g; }

     s/( no matching function [\s\w]+)/$1\n/g;

     if (/note:/) {
        if (!s/(.*candidate function not viable:)\s*/\e[38;5;8m$1\e[0m\n/) {
        if (!s/(candidate function)/\e[32m$1\e[0m/) {
             s/^(.* note: )/\n$1\n/g;
        }}
     }

     if (/\swarning:/) { s/^(.*warning:)/\e[32m$1\e[0m/i; if ($e) {++$e;}}
     elsif (/\serror:/) { s/^(.*error:)/\e[31m$1\e[0m/i; ++$e; }
     elsif (/undefined.*main.*/) { s/(reference to)/ref to/; if ($e) {++$e;} }
     elsif (!/main/ && /\sundefined/) { $_="\n$_"; ++$lcont;
         s/\s*(undefined reference to)\s*/\n\e[33m$1\e[0m\n/i; ++$e; }
     elsif (/ld returned/) { if ($lcont) { $_="\n$_"; }}

     print $_;
  }

  if (@l2 && @l2<5) { print @l2; @l2=(); }
  if ($n1 || @l2) {
     printf("%s\n=> skipping %d+%d lines.\n",$L,$n1,$#l2+1);
  }
};

# -------------------------------------------------------------------- #
# Wb,Feb18,24

sub parse_makelog_test {

  my $e1="\e[38;5;8m";
  my $er="\e[31m";
  my $ew="\e[35m";
  my $em="\e[0m";

  my ($q,@q,%cmp,%cxx,%xml,$ml,$gcc); my $agreed='';
  my $OK=0;

  $_='LIBRARY_PATH';      if ($ENV{$_}) { push(@q,$_); }
  $_='LD_LIBRARY_PATH';   if ($ENV{$_}) { push(@q,$_); }
  if (@q==2 && $q[0] eq $q[1]) { pop(@q); }
  $_='DYLD_LIBRARY_PATH'; if ($ENV{$_}) { push(@q,$_); }

  foreach (@q) { my @p=split(/:/,$ENV{$_});
     if (@p>1) 
          { print "  $_ = ",join("\n     ",'',@p),"\n\n"; }
     else { print "  $_ = ",@p,"\n\n"; }
  }

  foreach(<STDIN>) { s/^\s*//; chomp;
     if (/^(Found.*compiler|Building with)[:\s]*/) { my $i=$1;
        chomp($_=$'); s/^[\s'"]*//; s/[\.\s'"]*$//; ++$cmp{$_}; 
        $q=($i=~/Found/); $OK |= ($q? 1 : 2);
        printf "  $e1$i$em: $_%s\n", $q && $agreed ? " (agreed \@ $agreed)":'';

     }
     elsif (/^(Options file):\s*/) { my $i=$1; $OK|=8;
        chomp($_=$'); s/^[\s'"]*//; s/[\.\s'"]*$//; ++$xml{$_}; 
        s/$ENV{HOME}/~/;
        printf "  $e1$i$em: $_\n";
     }
     elsif (/^(CXX)\s*:\s*/) { my $i=$1;
        chomp($_=$'); s/^[\s'"]*//; s/[\.\s'"]*$//;
        @q=`$_ --version`;
        if ($? || @q<1) { die "\n  ERR failed to get '$_ --version'\n"; }
        chomp($q=$q[0]); $q=~s/\s*\(.*\)$//;
        printf "  $e1$i$em: $_ (%s)\n",$q; ++$cxx{$q}; 
     }
     elsif (/^agreed[:=\s*]*([\w\.]+)/) { chomp($agreed=$1); }
     elsif (/Yes .*'$agreed'/) { $OK|=4; }
  }
  print "\n";

  $_='matlab'; chomp($ml=`which $_ 2>/dev/null`);
  if ($ml ) {
     my $x=$ml; my $r='release??';
     if ($x=~s/bin\/matlab/VersionInfo.xml/ && -f $x) { 
        open(FH,'<',$x);
        foreach (<FH>) { if (/\<release\>(\w+)/) { $r=$1; last; }}
        close FH;
     };
     printf "  ${e1}which %-6s ->${em} %s (%s)\n",$_,$ml,$r;
  }

  $_='CXX'; $gcc=$ENV{$_}; if ($gcc) { $_="env \$$_" } else {
  $_='GCC'; $gcc=$ENV{$_}; if ($gcc) { $_="env \$$_" } else {
  $_='which gcc'; chomp($gcc=`$_ 2>/dev/null`); }}

  if ($gcc) { my ($wrn);
     if (!-f $gcc) { print " $ew WRN invalid $_ ($gcc)$em\n"; }
     else {
        @q=`$gcc --version`; if ($? || @q<1) {
            die "\n $er ERR failed to get '$q --version'$em\n"; }
        chomp($q=$q[0]); $q=~s/\s*\(.*\)$//;
        if (%cxx) { my @cxx=keys %cxx;
           if ($q ne $cxx[0]) {
              $wrn='Makefile uses different C++ compiler than system default';
           }
        }
        if (!$wrn) {
           if (ismac()) { if ($q!~/clang|apple|xcode/i) {
              $wrn='Makefile expects Xcode on macOS'; }
           }
           elsif (islinux()) { if ($q!~/gcc/i) {
              $wrn='Makefile expects gcc on linux'; }
           }
        }

        if (!$wrn) { $OK|=16; } else { $q="$ew$q$em"; }
        printf " $e1 %-12s ->$em $gcc ($q)\n",$_;
        if ($wrn) { print "\n $ew WRN $wrn$em\n"; }
     }
  }

  $q=$ENV{MATLAB_ROOT};
  if (!$q) { die "  ERR $q not defined\n"; }
  if (!$ml || $ml!~/$q/) {
     die "  ERR got mismatch of $q with 'which matlab'\n"; }

  if (($OK&7)==7) {
     $_='checks for matlab/compiler setup passed';
     if ($OK==31)
          { printf "\n \e[32m [OK] general $_.$em\n\n"; }
     else { printf "\n \e[32m --> other $_ @ q=$OK$em\n\n"; }

     if (!($OK&8)) {
        print "  (do make sure, though, to run 'mex -setup C++')\n\n";
     }
  }
  else {
     print "\n $er ERR please check matlab/compiler setup (%s)$em\n",
     "OK=$OK/7, agreed=$agreed";
     if (!($OK&12)) { print "  hint: did you run 'mex -setup C++'?\n"; }
  }
};

# -------------------------------------------------------------------- #

