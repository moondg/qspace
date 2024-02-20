
# to load use: require "$ENV{HOME}/bin/plib.pl";

  use strict;
  use warnings;

  use File::stat;
  use Cwd 'realpath';

# -------------------------------------------------------------------- #
# usage: getenv($vname [,default_value])
  sub getenv {
     if (@_<1 || @_>2) { die "\n  ERR getenv() invalid usage\n  ".
        "ERR ".join(', ',@_).sprintf(' (nargs=%g)\n',$#_+1);
     }
     if ($ENV{$_[0]})
          { return $ENV{$_[0]}; }
     else { return (@_>1 ? $_[1] : ''); }
  }

# -------------------------------------------------------------------- #
# usage: set_limit_aux( ulimit_opt, new_val) // Wb,Jul01,19
# e.g. set_limit_aux('-n',4096);

  sub set_limit_aux {
     if (@_!=2 || $_[1]!~/^\d/ || $_[0]!~/^-(\w)$/) { wbdie("invalid usage"); }
     my $o=$1; my $x=$_[1];
     my @ll=`ulimit -$o ; ulimit -H$o`;

     if ($? || @ll!=2) { wblog("ERR ulimit call returned error !?"); }
     else { chomp(@ll);
        if ($ll[1]=~/unlimited/i) { $ll[1]=1E99; }
        if ($x>$ll[1]) { $x=$ll[1]; }
        if ($ll[0]<$x) { printf "  ulimit -S -$o $x   # old: %g\n",$ll[0]; }
        else { printf "# ulimit -$o having $ll[0] / $x - ok.\n"; }
     }
  }

# -------------------------------------------------------------------- #
# reduce hostname to few-character ID

sub hostid() {
  my $id='xxx';
  $_=`hostname`; chomp; if ($?) { die; }

  if (/th-(cl|ws)-(.*)/) { $id=$1; }
  else {
     s/[\s_.-]+//g;
     $id=substr($_,0,3);
  }
  return $id;
};

# -------------------------------------------------------------------- #

sub islinux {
   foreach (`uname -s`) { if (/Linux/i) { return 1; }}
   return 0;
};

sub ismac {
   foreach (`uname -s`) { if (/Darwin/i) { return 1; }}
   return 0;
};

sub getARCH {
   my @ll=`uname -sm`;
   foreach (@ll) {
      if (/Linux/i) {
         if (/x86_64/) { return 'glnxa64'; }
         if (/i686/  ) { return 'glnx86';  }
      }
      elsif (/Darwin/i) {
         if (/x86_64|arm64/) { return 'maci64'; }
      }
   }
   wblog("ERR %F() failed to determine ARCH");
   return '';
};

sub getMEXEXT {
   my @ll=`uname -sm`;
   foreach (@ll) {
      if (/Linux/i) {
         if (/x86_64/) { return 'mexa64'; }
         if (/i686/  ) { return 'mexglx';  }
      }
      elsif (/Darwin/i) {
         if (/x86_64|arm64/) { return 'maci64'; }
      }
   }
   wblog("ERR %F() failed to determine MEXEXT");
   return '';
};

# -------------------------------------------------------------------- #
# usage: wblog([k,] fmt, args);
#
#    where k specifies how far up in the caller stack to go
#    (default: k=0); fmt+args are handed over to and thus
#    must be in the format of a call to printf(fmt,args)
#
#    if fmt starts with ERR, WRN, TST, ' * ', etc., this is
#    interpreted, and results e.g. in color coding for ERR and WRN.
#    also: %F is replaced by the caller's function name.
#
# Wb,May10,18

sub wblog {

  my ($k,@S,@S_,$fmt,$f,$s,$h,$x1,$w,$e); my $nh=20;
  my $in='';  # leadig newlines

  my $nl="\n";

  $k=0; if (@_ && "$_[0]"=~/^\d+$/) { $k=shift; }
  @S = caller($k); @S=($S[1],'',$S[2]);
  @S_= caller($k+1);

  $fmt=shift;
  if ($fmt=~s/^(\s*\n)//) { my $q=$1; $q=~s/ //g;  $in=$q; }

  if ($fmt=~s/\\+$//) { $nl=''; }
  elsif ($fmt=~s/(\s+)$//) {
   # trailing newlines (incl. space after new line)
     my $q=$1; $q=~s/^ +//; $nl.=$q;
  }

  if ($fmt=~s/^\s*([\w<>!]{3})( |$)//) { $w=$1;
     if    ($w eq 'ERR') { $x1="\e[31m"; $e=1; }
     elsif ($w eq 'WRN') { $x1="\e[35m"; $e=1; }
  }
  elsif ($fmt=~s/^( *\*+) +//) { $w=sprintf('%-3s',$1); } # e.g. ' * '
  elsif (@ARGV && length($fmt)<=3 && $fmt!~/\%/) {
      $w=sprintf('%-3s',$fmt); $fmt=shift; }
  else { $w='   '; }

  $S[0]=~s/.*\///;
  $S[2]=sprintf(":%-3d",$S[2]);

  if (@S_ && $S_[3]) {
     my $f=$S_[3]; $f=~s/.*:://;
     if ($fmt=~/$f/) { print STDERR
        "\n !? wblog() may want to replace sub-name $f by \%F".
        "\n !? in fmt=$fmt)\n\n";
     }
     elsif ($fmt!~s/\%F/$f/g) { my $n=$nh-length($S[2]);
        my $l0=length($S[0]);
        my $l1=length($f);
        if    ($l0>$n) { $S[0]=substr($S[0],0,$n-1)."'"; $S[1]=''; }
        elsif ($l0>$n-2) { $S[1]=''; }
        elsif ($l0+$l1>=$n) { $S[1]='/'.substr($f,0,$n-$l0-2)."'"; }
     }
  }
  $h=sprintf("%-${nh}s ",join('',@S));

  if ($x1) { my $x2="\e[0m"; $h.=$x1.$w.' ';
     $fmt=~s/\n/$x2\n$h/g;
     $fmt=$x1.$fmt.$x2;
  }
  else { $h.=$w.' '; $fmt=~s/\n/\n$h/g; }

  if (@_) { $s=sprintf($fmt,@_) } else { $s=$fmt; }
  $s=$in.$h.$s.$nl;

  if ($e) { print STDERR $s; }
  else    { print STDOUT $s; }
};

# -------------------------------------------------------------------- #
# print caller stack // Wb,Jun06,18

sub print_stack {

  my ($i,@q,@ll,$FL,$sflag); my $k=0;
  my $is='     '; my $h='';

  while (@_) {
     if ($_[0]=~/^\d+$/) { $k=shift; }
     elsif ($_[0]=~/^-s/) { $sflag=shift; }
     else {
        die "ERR invalid usage ($_[0])\n ";
     }
  }

  for ($i=0; $i<99; ++$i) { @q=caller($i);
     if (!@q) {
        push(@ll, sprintf($is."[%g] %-20s -\n",$i-$k,$FL));
        last;
     }

     $q[3]=~s/.*:://;
     if ($i>$k) { push(@ll,sprintf($is."[%g] %-20s %s()\n",$i-$k,$FL,$q[3])); }
     if ($i==1) { $h=sprintf("%s -> %s()",$FL,$q[3]); }
     $q[1]=~s/.*\///; $FL=$q[1].':'.$q[2];
  }

  @ll=reverse(@ll);
  if (!$sflag) { print "\n",$h,": current caller stack\n\n",@ll,"\n"; }

  return @ll;
};

# -------------------------------------------------------------------- #
# usage: wbdie([k,] fmt, args);
#
#    where k specifies how far up in the caller stack to go
#    (default: k=0); fmt+args are handed over to and thus
#    must be in the format of a call to printf(fmt,args)
#
#    once the message is written, a regular `die' is issued.
#
# Wb,May12,18

sub wbdie {

  my ($q,$FL,@ll,@args,$qflag);

  foreach $q (@_) {
     if ($q eq '-q') { $qflag=1; }
     else { push(@args,$q); }
  }

  my $k=0; if (@args && "$args[0]"=~/^\d+$/) { $k=shift(@args); }

  my @S=caller($k);
  my @S_=caller($k+1);

  @S=($S[1],$S[2],'');

  my $fmt; if (@args)
       { chomp($fmt=shift(@args)); $fmt=~s/\bERR\b//g; }
  else { $fmt='!?'; }

     $fmt=~s/^\s*/\n/;
     $fmt=~s/\s*\n\s*/\n  ERR /g;

  $S[0]=~s/.*\///;
  if (@S_ && $S_[3]) {
     my $f=$S_[3]; $f=~s/.*:://;
     if ($fmt=~/$f/) { print STDERR
        "\n !? wbdie() may want to replace sub-name $f by \%F".
        "\n !? in fmt=$fmt)\n\n";
     }
     elsif ($fmt!~s/\%F/$f/g) { $S[2]=' -> '.$f.'()'; }
  }

  @ll=print_stack(1,'-s');

  if (!$qflag) { printf STDERR "\n%s:%d%s\n",@S; }
  if (@ll>1) { print STDERR "\n",@ll; }

  if (@args) { $q=sprintf($fmt,@args); } else { $q=$fmt; }
  print STDERR "\e[31m".$q."\e[0m";

  die "\n";
};

# -------------------------------------------------------------------- #
# example usage
#    die usage($0);
#    die usage(__FILE__,__LINE__);
#    die usage(__FILE__,__LINE__,"Invalid usage: ",@ARGV);

sub usage() {
  if (!@_ || ! -f $_[0]) { die "\n  Usage: usage(\$0)\n "; }

# read header info from specified script
# which includes all initial comments that start with ^#
  open(FID, $_[0]) or die "\n  ERR failed to open file $_[0]\n. ";
  my $name=$_[0]; $name=~s/.*\///g;
  my @header=("\n");

  foreach my $l (<FID>) {
     if ($l=~/^#!/) { next; }
     if ($l=~/^#/) { $l=~s/^#/ /; $l=~s/NAME/$name/g; push(@header,$l); }
     else { last; }
  }
  close(FID);

  if (@_>1 && $_[1]=~/^\d+$/) {
     my $f=shift; repHome($f);
     push(@header,"\n  See $f:".(shift)."\n");
  }

  if (@_) { my @ll=@_;
     push(@header,"\n");
     foreach (@ll) { s/^\s+//; push(@header,"  $_\n"); }
  }
  return join('',@header);
};

# -------------------------------------------------------------------- #
# replace HOME and other standard PATHs by shortcuts
# switched from in-place to copy // Wb,Oct19,16

  sub repHome { my $q; my @DD=@_;
     foreach (@DD) {
        $q='HOME';     if ($ENV{$q}) { s/$ENV{$q}/~/g; }
        $q='LMA';      if ($ENV{$q}) { s/$ENV{$q}/\$LMA/g; }
        $q='RC_STORE'; if ($ENV{$q}) { s/$ENV{$q}/\$RCS/g; }
        $q='RC_SYNC';  if ($ENV{$q}) { s/$ENV{$q}/\$RCS/g; }
     }

     if (@DD==1)
          { return $DD[0]; }
     else { return @DD; }
  }

# -------------------------------------------------------------------- #
# run system command in safe mode
# -q  quiet (do not show stdout, unless error occurred)
# -c  continue in case of errror

sub rsys {

  my $qflag=0; my $vflag=0; my $cflag=0;
  while (@_) { my $q=$_[0];
     if ($q eq '-q') { $qflag=1; shift; }
     elsif ($q eq '-v') { $vflag=1; shift; }
     elsif ($q eq '-c') { $cflag=1; shift; }
     else { last; }
  }
  if (@_!=1) { die "\n  ERR invalid usage: rsys [opts] command\n"; }

  my $cmd=shift; my @ll=();
  if ($vflag) { print "  cmd: $cmd\n"; }
  if ($qflag) { @ll=`$cmd`; } else { system($cmd); }

  if ($? && !$cflag) {
     die join('',@ll)."\n  ERR invalid command '$cmd' ($?)\n";
  }
  return $?;
};

# -------------------------------------------------------------------- #

sub is_cluster_job {

  my $q=0;
  foreach ('PBS_JOBID','PBS_O_WORKDIR') { if ($ENV{$_}) { ++$q; }}
  if ($q>1) { return 1; }

  $q='';
  foreach ('SLURM_ARRAY_JOB_ID','SLURM_JOB_ID','SLURM_JOBID') {
     if ($ENV{$_} && $ENV{$_}=~/^\s*([\d\._-]+)\s*$/) { $q=$1; last; }
  }; if (!$q) { return 0; }

  $q='';
  foreach ('SLURM_ARRAY_TASK_ID','SLURM_TASK_PID') {
   # tid (if no array is set, use process ID)
     if ($ENV{$_} && $ENV{$_}=~/^\s*(\d+)\s*$/) { $q=$1; last; }
  }; if (!$q) { return 0; }

# $q='';
# foreach (
#   'SLURM_JOB_NODELIST','SLURM_NODELIST', # executing host
#   'SLURM_NPROCS','SLURM_NTASKS'  # see get_nslots() below
# ){
#    if ($ENV{$_}) { $q=$ENV{$_}; last; }
# }; if (!$q) { return 0; }

  return 1;
};

# -------------------------------------------------------------------- #

sub env_num_threads {

  my $np=-1;
  if (is_cluster_job()) {
    foreach ('SLURM_NPROCS','SLURM_NTASKS','PBS_NP','PBS_NUM_PPN') {
       if ($ENV{$_} && $ENV{$_}=~/^\s*(\d+)\s*$/) {
          if ($np<$1) { $np=$1; if ($np>0) { return $np; }}
       }
    }
  }
  else {
    foreach ('NSLOTS','OMP_NUM_THREADS') {
       if ($ENV{$_} && $ENV{$_}=~/^\s*(\d+)\s*$/) {
          $np=$1; if ($np>0) { return $np; }
       }
    }
  }

  if ($np<0) {
     my $f='/proc/cpuinfo';
     if (-f $f) { 
        my @ll=grep(/processor/,`cat $f`);
        $np=int($#ll+1)/2; if ($np>4) { $np=4; }
     }
  }

  return $np;
};

# -------------------------------------------------------------------- #
# Wb,Apr12,18

sub remove_paths {

   my ($pat,$x,$P,@P,$q,@p,@q,@out,$bflag); my $vflag=1;

   foreach (@_) {
      if (/^-/) {
         if ($_ eq '--bash') { $bflag=1; }
         elsif ($_ eq '-q') { $vflag=0; }
         elsif ($_ eq '-v') { $vflag+=1; }
         else { die "\n".'  invalid usage: remove_paths ($0)'."\n"; }
      }
      else { push(@P,$_); }
   }

   if (@P<2) { die "\n".'  invalid usage: remove_paths ($0)'."\n"; }
   $pat=shift(@P);

   foreach (@P) {
      if ($bflag) {  # assume got names of environmental variables
         if (!exists $ENV{$_}) {
            if ($vflag) { print STDERR 
            "  WRN ignoring non-existing environemental variable $_\n";}
            next;
         }
         $P=$ENV{$_};
      }
      else { $P=$_; }

      if ($P=~/([|:;])/) 
           { $x=$1; @p=split(/$x/,$P); }
      else { $x=''; @p=($P); }

      @q=grep(!/$pat/,@p);
      $q=join($x,@q);

      if ($bflag) { if ($#q!=$#p) {
         if (@q)  # only export command if path actually changed
              { push(@out,"export $_=".$q); }
         else { push(@out,"unset $_"); }
      }}
      else { push(@out,$q); }
   }

   return @out;
}

# -------------------------------------------------------------------- #
# Wb,Apr12,18

sub add_paths {
   my ($N,$x,$P,%P,@p,@x,@X,$bflag,$prepend); my $vflag=1;

   foreach (@_) {
      if (/^-/) {
         if ($_ eq '--bash') { $bflag=1; }
         elsif ($_ eq '-q') { $vflag=0; }
         elsif ($_ eq '-v') { $vflag+=1; }
         elsif ($_ eq '--pre') { $prepend=1; }
         else { die "\n".'  invalid usage: add_paths ($0)'."\n"; }
      }
      else { push(@X,$_); }
   }

   if (@X<2) { die "\n".'  invalid usage: add_paths ($0)'."\n"; }
   $P=shift(@X);

   if ($bflag) {  # assume got name of environmental variables
      if (!exists $ENV{$P}) { print STDERR 
      "  ERR ignoring non-existing environemental variable $P\n"; }
      $N=$P; $P=$ENV{$P};
   }

   if ($P=~/([|:;])/) 
        { $x=$1; @p=split(/$x/,$P); }
   else { $x=':'; @p=($P); }

   foreach (@p) {
      if (!-d $_) { die "  ERR add_paths() invalid path '$_' !?\n"; }
      s/[\/\s]*$//; $P{$_}++;
   }

   foreach (@X) {
      if (!-d $_) { die "  ERR add_paths() invalid path '$_' !?\n"; }
      s/[\/\s]*$//;
      if (++$P{$_}==1) { push(@x,$_); }
   }

   if (!@x) { 
    # only export command if path actually changed
      if ($bflag) { $P=''; }; return $P;
   }

   if ($prepend)
        { push(@x,@p); $P=join($x,@x); }
   else { push(@p,@x); $P=join($x,@p); }

   if ($bflag)
        { return "export $N=".$P; }
   else { return $P; }
};

# -------------------------------------------------------------------- #
1;  # keep this
