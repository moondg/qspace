#!/usr/bin/perl -w
# NAME: acts like | head -20 but color code on warning and error
# Wb,Apr22,10

  my ($n,$n1,$flag,@l2,$L); my $ncols='80'; my $e=0;
  my $H=$ENV{HOME};
  my $N=$ENV{NMAX_MAKELOG}; if (!$N) { $N=20; }
  my $MML=$ENV{MYMATLAB};

  my $LFLAG=$ENV{'LOG_FLAG'}; my $lflag=1;
  if (defined($LFLAG) && $LFLAG eq '-q') { $lflag=0; }

  foreach(<>) {
   # skip default errors when compiling mex file from linux command line prompt
     if (/lib(mex|mat|mx|mwm|mws).*\.so: undefined reference/) { ++$n1; next; }
     if (/^\/usr\/bin\/ld: warning: lib[\w]+.so[.\d]*, needed by.*try using/) {
        ++$n1; next;
     }
     if (/^\/usr\/local.*matlab.*libmwlapack.*undefined reference/) {
        ++$n1; next;
     }

     if (/MEX completed successfully/) { next; }
     if (s/Building with '(.*)'/built with $1/
     ){ chomp; s/^\s*//; print "\r  ... $_ ... \r";  next; }

   # show max N (default 20) lines total
   # but skip forward to show at least one full error block
     if (++$n>$N && $e>1) { push(@l2,$_); next; };

     if (++$lflag==2) {
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
     elsif (/\serror:/) { s/^(.*error:)/\e[31m$1\e[0m/i; ++$e; } # ;47
     elsif (/undefined.*main.*/) { s/(reference to)/ref to/; if ($e) {++$e;} }
     elsif (!/main/ && /\sundefined/) { $_="\n$_"; $flag++;
         s/\s*(undefined reference to)\s*/\n\e[33m$1\e[0m\n/i; ++$e; }
     elsif (/ld returned/) { if ($flag) { $_="\n$_"; }}

     print $_;
  }

  if (@l2 && @l2<5) { print(@l2); @l2=(); }

  if ($n1 || @l2) {
     printf("%s\n=> skipping %d+%d lines.\n",$L,$n1,$#l2+1);
  }

