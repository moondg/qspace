#!/usr/bin/perl -w
# Usage: NAME <prog>
#
#    get output directory of given progs
#
# Wb,May18,20

  my ($f,$F,$d,$e,$p,$n1,$n2,$nx,@q,@qq,@ff,@fx,@mopts,%FD);
  my ($vflag,$mflag); $vflag=0;

  my @odir=('../bin','../util');
  my $odef=$ENV{MXC_DIR};
  if (!$odef) { $odef=$odir[0]; }

  my $me=$0; $me=~s/.*\///;

  foreach (@ARGV) {
     if (/^-[h\?]$/) { die usage($0); }
     elsif (/^-(v)$/i) { $vflag+=($1 eq 'v' ? 1 : 2); }
     elsif (/^--odir$/) {
         $mflag+=1;
     }
     elsif (/^--make-/) { my $x=$';
        if    ($x=~/info/) { disp_info(@ARGV); exit 0; }
        elsif ($x=~/defs/) { disp_defs(@ARGV); exit 0; }
        else { die("invalid option $_"); }
     }
     elsif (/^-/) { die("invalid option $_"); }
     else { push(@ff,$_); }
  }

  if (!@ff) { $vflag+=2; }

  if ($mflag) {
     if (@ff>1) { @mopts=split(/[:;,\|]/,' '.$ff[1].' '); }
     if (@ff!=2 || @mopts!=4) {
        die sprintf("invalid usage %s -> 4 args expected: '%s' (%d)",
        $ff[0], join("', '",@ff[1 ..  $#ff]),$#mopts+1);
     }
     @ff=($ff[0]);

     foreach (@mopts) { s/^\s+//; s/\s+$//; }
  }

  if ($vflag>1) { print STDERR "\n"; }

  foreach $d (@odir) {
     opendir(DH,$d) || die("invalid directory ($d)");
     @fx=grep(/\w\.m\w*$/,readdir(DH));
     close(DH); $e=0; $n1=$n2=$nx=0;

     foreach $f (@fx) { $F="$d/$f";
        if (-d $F) { ++$nx; next; }
        if (!-f $F) {
           printf STDERR "  WRN %s not a file !?\n", $F;
           ++$n2; next;
        }

        if ($f=~/\.mex\w*$/) { ++$n1;
           if (!-x $F && !$mflag) {
           printf STDERR "  WRN %s not executable\n", $F; }
        }
        elsif ($f=~/\.m/) { ++$n2; }
        else { ++$nx; next; }

        if ($FD{$f}) {
           printf STDERR "  ERR multiple occurance of %s also in %s\n",
           $F,$FD{$f}; ++$e;
        }
        else { $FD{$f}=$d; }
     }
     if ($e || $vflag>1) {
        printf STDERR "  %3s $me %-12s %g/%g executables\n",
        $e ? 'WRN':'ok.',$d,$n1,$n1+$n2;
     }
  }

  @fx = sort keys %FD;

  if ($vflag>1) {
     printf STDERR "\n  --> %g files total\n",$#fx+1;
     if (@ff) { print "\n"; }
  }

  foreach $f (@ff) {
     if (-d $f) { next; }
     $f=~s/.*\///;

     $p=$f; $p=~s/\..*$//; $p=~s/.*\///;
     @q=grep(/^$p\b/,@fx);

     if (@q>1) {
        my $i=1; my $D0=$FD{$q[0]};
        for (; $i<=$#q; ++$i) { if ($FD{$q[$i]} ne $D0) { last; }}
        if ($i<@q) { printf STDERR
           "  WRN $me got %g entries for $f at different paths:\n",$#q+1;
           foreach (@q) { printf STDERR "   %s / %s\n",$FD{$_}, $_; }
           print "\n";
        }
        else { @q=($q[0]); }
     }

     if (@q==1)  { 
        if ($vflag>1)
             { printf "  %-30s -> %s\n",$f,$FD{$q[0]}; }
        else { push(@qq,$FD{$q[0]}); }
     }
     elsif (!@q) { printf STDERR
     "  WRN $me failed to match %s -> using default odir %s\n",$f,$odir[0];
        push(@qq,$odir[0]);
     }
  }

  if (@ff!=1 || !$mflag) { exit $#ff+1; }

  my ($q,$l,@lb,@gg,@mo);

  $l=' '.$mopts[2];
  $l=~s/\s{2,}/ /g;

  if ($l=~/blas|lapack/) { my $m=0; $_=$l;
     if (s/\s(-l|-lmw|\/[^\s]+\/lib)blas\b[^\s]*//  ) { ++$m; }
     if (s/\s(-l|-lmw|\/[^\s]+\/lib)lapack\b[^\s]*//) { ++$m; }
     if ($m==2) { push(@l,'-lblas/lapack'); $l=$_; }
  }
  if ($l=~/mpfr|gmp/) { my $m=0; $_=$l;
     if (s/\s(-l|-lmw|\/[^\s]+\/lib)gmp\b[^\s]*// ) { ++$m; }
     if (s/\s(-l|-lmw|\/[^\s]+\/lib)mpfr\b[^\s]*//) { ++$m; }
     if ($m==2) { push(@l,'-lmpfr'); $l=$_; }
  }

  chomp($_=$l);
  if (!/^\s*$/) { s/^\s*//; push(@l,$_); }

  push(@ff,@qq); $vflag=0;

  if (($q=$ENV{DEBUG})) {
     if    ($q>=3) { push(@gg,'-V'); $vflag+=2; }
     elsif ($q>=2) { push(@gg,'-v'); $vflag+=1; }
  }
  if ($ENV{ML_DEBUG}) {
     if (!@gg) { push(@gg,'-g'); }
     else { print(STDERR "\n  WRN using DEBUG (ignoring ML_DEBUG)\n\n"); }
  }

  if ($_=$mopts[0]) { s/.*\///;
     s/\.app$//i; # osx/mac has trailin \.app
     push(@mo,substr($_,-5));
  }

  if ($_=$mopts[1]) {
     push(@mo,substr($_,3));
  }

  if ($_=$mopts[3]) {
     s/^\s*//; s/\s*$//; s/\s+/ /g;
     push(@mo,$_);
  };

  push(@mo,join(' ',@l),@gg);
  if (@ff<2) { $ff[1]='!?'; }

  printf STDERR "  mex %-36s %-16s => %s\n",join('; ',@mo),@ff;
  if ($vflag) { printf STDERR "\n" x $vflag; }

  print join(' ',@qq,@gg);

# ---------------------------------------------------------- #
# Wb,Dec16,24

sub get_host_name {
   my ($q,$H);
   foreach $H (@_,$ENV{HOSTNAME},$ENV{HOST}) {
      if ($H) { $H=~s/\..*//; }
      if ($H && $H!~/^\d[\d-]*(\.|$)/) { last; }
   }
   if (!$H) { my $ismac=0;
      foreach (`uname -s`) { if (/Darwin/i) { ++$ismac; }}
      if ($ismac)
           { foreach $q (`scutil --get ComputerName`) { chomp($H=$q); }}
      else { $H='(hostname)'; }
   }

   $q=$ENV{HOSTTYPE}; if (!$q) { chomp($q=`arch`); }
   if ($H=~s/-Air//) {
      my @q=`sysctl -n sysctl.proc_translated 2>/dev/null`;
      if (!$? && $q[0]) { $q.="/rosetta"; }
   }
   if ($q) { $H="$H ($q)"; }

   return $H;
};

# ---------------------------------------------------------- #
# Wb,Apr27,20

sub disp_info {
   my (@q,%q,$tag,@tag,@ll);

   $_=join(' ',@_); s/\s*=\s*/= /g;
   @q=split(/\s+/,$_); $tag=''; push(@tag,$tag);

   select STDERR;

   foreach (@q) {
      if (/^--make-info/) { next; }
      if (/^\/[^\s]*\/(matlab\/[^\s]*)/i) { $_=$1;
         s/\.app$//;
      }
      if (/^bin(\w*)$/) { $_=$1; }
      if (/([A-Za-z]\w+)=$/) { $tag=$1;
         if (!defined $q{$tag}) { push(@tag,$tag); }
      }
      else { push(@{$q{$tag}},$_); }
   }

   foreach $tag (@tag) { if (!defined $q{$tag}) { next; }
      @q=@{$q{$tag}}; if (@q) {
      push(@ll, ($tag ? "\e[01;30m[".$tag."]\e[0m ":'').join(' ',@q)); }
   }

   if (@ll) { $_=join(' ',@ll);
      if (length($_)<128) { print "  \e[01;30musing\e[0m $_\n"; }
      else {
         foreach (@ll) {
            if (/(\[\w+\])\s*\s*(.*)/)
                 { printf "  %-8s   %s\n",$1,$2; }
            else { printf "  %s\n",$_; }
         }
      }
   }
};

# ---------------------------------------------------------- #
# usage / input format
# VAR needs to consist of uppercase letters and underscore
#  
#    <VAR>=def     => -D<VAR>
#
#    <VAR>=        => -U<VAR>
#        WRN! this does not work as input to make, though, since
#       `make VAR=' keeps VAR undefined inside Makefile (ifdef returns 0)
#    <VAR>=undef   => -U<VAR>
#
# Wb,May13,20

sub disp_defs {
   my (@ll,$e,@args,$vflag);
   foreach (@_) {
      if (/^--make-defs/) { next; }
      if (/^-(v)$/i) { $vflag += ($1 eq 'v' ? 1:2); }
      else {
         s/=undef$/=/; s/=def$//;
         push(@args,$_);
      }
   }

   if (!grep(/HOST_NAME/,@args)) {
      push(@args,'-DHOST_NAME='.get_host_name());
   }
   if ($ENV{DBSTOP}) {
      if (!grep(/DBSTOP/,@args)) { push(@args,'DBSTOP'); }
   }

   foreach (@args) {
      if (/^([A-Z_]{4,})(=|$)/) {
         my $a=$1; my $x=$2;
         my $v=$'; 
            if ($a=~/CXX_FLAGS/) { $v=~s/(^|\s)-Wall\b//; }
            $v=~s/^\s*//; $v=~s/[\s\x0A-\x0D]*$//g; $v=~s/ +/ /g;

         if ($v) { if (!$x) { $e=1; } else {
            if ($a=~/HOST_NAME/) { $v=get_host_name($v); }
            push(@ll, $v=~/[,;\s]/ ? "-D$a='$v'" : "-D$a=$v");
         }}
         else {
            push(@ll, $x? "-U$a" : "-D$a");
         }
      }
      elsif (/^-[UD]/) { push(@ll,$_); }
      else { $e=2; }

      if ($e) {
         print STDERR "  ERR $me: ignoring $_ (e=$e)\n";
         $e=0;
      }
   }

   if ($vflag && @ll) {
      if ($vflag>1)
           { print STDERR join("\n   ",'',@ll,'',''); }
      else {
         my @q=grep(!/DATE|MLVER|HOST_NAME/,@ll); if (@q) {
         print STDERR "  \e[01;30musing [DEFS]\e[0m ",join(' ',@q,"\n\n"); }
      }
   }

   if (@ll) { print join(' ',@ll); }
};

# ---------------------------------------------------------- #

