#!/usr/bin/perl
# Usage: eval "$(matlab_setup.pl [opts])"
# 
#    This script double checks current matlab bash environment.
#    It generates further relevant bash settings in text form
#    used by matlab_setup.sh required for running matlab.
# 
#    Note that the matlab and C++ compiler version must already
#    be setup and available.
# 
# Options
# 
#    -t   test flag (prints status info to stderr)
#    -q   quiet mode, i.e., suppress certain warnings
#         (intended for calling *this in Makefiles)
#
#    -LD  enforces buildup of full LD_LIBRARY_PATH
#    -ld  using matlab's default LD_LIBRARY_PATH (once it starts)
# 
# Wb,Jan09,19

# [02/19/2024] removed --pre-run
# [02/19/2024] removed gcc/GCC checks and gver/ML_GCC_VERSION
# since 'mex -setup' may use other C++ compiler anyway
# => rather use cd Source && make test, instead

  use strict; use warnings;
  use Cwd 'realpath';

  $|=1;  # auto flush after every printf

  my $P=realpath($0); $P=~s/\/[^\/]*$//;
  require "$P/plib.pl";

  my $nerr=0; my $nwrn=0; my $vflag=1;
  my ($q,$mver,$tflag, $ldflag,@task);

  while (defined ($_=shift)) {
     if (/^-[h\?]$/) { die usage($0); }
     elsif (/^-(t)$/i) { $tflag+=($1 eq 'v'?1:2);  }
     elsif (/^-q$/) { $vflag=0; }
     elsif (/^-(ld)$/i) { $ldflag+=($1 eq 'ld'?1:2);  }
     elsif (/^-/) { push(@task,$_); }
     else { 
        die usage(__FILE__,__LINE__,"got extra arguments: ",$_,@ARGV);
     }
  }

# -------------------------------------------------------------------- #
# check MATLAB_ROOT and consistency with `which matlab`

  my ($matlab,$mex);
  my $mlr='MATLAB_ROOT';
  my $MLR=$ENV{$mlr};

  chomp($matlab=`which matlab 2>/dev/null`);
  chomp($mex=`which mex 2>/dev/null`);

  if (!$matlab || ($MLR && $matlab!~/$MLR/)) {
  if (!$ENV{QS_CONFIG_ML_SH}) { $_="$P/matlab_setup_user.sh";
     if (!-f $_) { $q='MYMATLAB';
        if ($ENV{$q} && s/$ENV{$q}\/*/\$$q\//) 
             { $q=" (having $q = ".repHome($ENV{$q}).')'; } 
        else { $q=''; s/$ENV{HOME}/~/; }

        print STDERR join('  ',"\n",
        "system/matlab_setup.sh$q:\n\n\e[31m",
        "    Please setup the QSpace system environment first\n",
        "    i.e., create/edit file $_\n",
        "    e.g., based on $_-template\n",
        "    or set environmental variable QS_CONFIG_ML_SH.\e[0m\n\n",
        "See QSpace documentation for more detailed information [App. B2].\n");
        exit 1;
     }
  }}

  if (@task>1) {
     wbdie("invalid usage (multiple tasks: %s)",join(' ',@task));
  }

  foreach (@task) {
     if (/^--clear$/) {
        my $ML=$ENV{MYMATLAB};
        if ($MLR) { print join("\n", remove_paths('--bash','-q',$ML,
           'PATH','LD_LIBRARY_PATH','XAPPLRESDIR'));
        }
        exit 0;
     }
     elsif (/^--check$/) { exit check_ENV(); }
  }

# insist that `which matlab' exists and agrees with MATLAB_ROOT
  if (!$matlab) { ++$nerr; wblog(
     "ERR command 'matlab' not available on PATH (see matlab_setup*.sh)"); }
  if (!$mex) { ++$nwrn; wblog(
     "WRN command 'mex' not available on PATH (see matlab_setup*.sh)"); }

  if (!$MLR) {
     if ($matlab) {
        $MLR=$matlab; $MLR=~s/\/bin.*matlab.*//;
        export_plain($MLR,$matlab);
     }
     else { ++$nerr; wblog(
       "ERR env $mlr not defined (see matlab_setup*.sh)");
     }
  }

  if (!$nerr) {
     if (!-d $MLR.'/bin/') { ++$nerr; wblog(
        "ERR invalid env $mlr=%s (missing subdirectory ./bin",$MLR); }
     elsif ($matlab!~/$MLR/) { ++$nerr; wblog(
        "ERR invalid env $mlr=%s (%s)",$MLR,$matlab); }
     elsif ($MLR!~/matlab/i) { ++$nwrn; wblog(
        "WRN unexpected value for env $mlr=%s",$MLR); }
     if ($matlab=~/[_\/rR]?(\d{4}\w?)[_\.\/]/) { $mver=$1; }
     if ($mex!~/$MLR/) { ++$nwrn; wblog(
        "WRN unexpected mex command (outside matlab folder?)");
     }
  }

  if (!$mver && $MLR) { 
   # checkout $(MATLAB_ROOT)/VersionInfo.xml // Wb,Feb15,24
     my $v="$MLR/VersionInfo.xml";
     if (-f $v) { my (%R,@R);
        open(FH,$v);
        foreach (<FH>) {
           if (/<release>R?(\w+)<\/release>/) { $mver=$1; ++$R{$1}; }
           if (/R(\d{4}\w?)/) { ++$R{$1}; }
        }; close FH;

        if (%R) { @R=keys %R;
           if (@R==1) { if (!$mver) { $mver=$R[0]; }}
           else {
              wblog("WRN failed to auto-determine matlab release (%s)",
              join(', ',sort @R));
           }
        }
     }
  }

  if (!$nerr) {
  if (!$mver) { ++$nwrn; wblog(
     "ERR failed to derive matlab version from env\n$mlr = $MLR\n%s",
     "this is required for MEX compilation, etc.");
  }}

  if ($nerr) { exit 1; }

  if ($nwrn) { exit 2; }

# -------------------------------------------------------------------- #
# check MYMATLAB, MEX, MCC
# MYMATLAB is my Matlab home directory (e.g. points to startup.m)

  my ($MM,$ML, $MC,$MCC, $MX,$MEX);
  $MM='MYMATLAB'; $ML =$ENV{$MM};
  $MC='MCC';      $MCC=$ENV{$MC};
  $MX='MEX';      $MEX=$ENV{$MX};

  if (!$ML) { ++$nerr; $ML=''; wblog("ERR env $MM not defined\n".
     "(this is the local matlab directory e.g. that contains startup.m)"); }
  elsif (!-d $ML) { ++$nerr; wblog("ERR invalid directory $MM=$ML"); }

  $q=$ML."/$MC";
  if (-d $q) { $MCC=$q; }
  elsif ($MCC) { if (!-d $MCC)  { ++$nerr; wblog(
     "ERR invalid directory env $MC=$MCC"); }
  }
  else { $MCC=''; if (!$nerr) { ++$nerr; wblog(
     "ERR env $MC not defined (nor does ./MCC exist)");
  }}

  $q=$ML.'/Source';
  if (-d $q) { $MEX=$q; }
  elsif (!$nerr && $MEX && !-d $MEX) { ++$nerr; wblog(
     "ERR invalid directory env $MX=%s\n",$MEX); 
  }

  my $mar=getARCH(); my $mxt=getMEXEXT();
  if (!$mxt) { ++$nerr; }

  if ($tflag) {
     printf STDERR "\n  %-20s %s\n",$MM,repHome($ML);
     if ($mver) { printf STDERR
        "  %-20s %-14s (%s)\n",'matlab version',$mver,$matlab;
     }
  }

  if ($mver) {
     $q=$mver; $q=~s/^20//;
     if ($q=~s/([a-z])$//i) {
        my $a=lc($1); $q+=(ord($a)-ord('a')+1)/10;
     }

     export_plain('MATLAB_VERNUM',$q);
     export_plain('MATLAB_VERSION',$mver);

     if ($q>=18) {
      # required for MEX files to avoid compiler errors with functions
      # not found such as: mxComplexDouble, mxGetComplexDoubles, etc.
      # NB! matlab>=2018 has interleaved complex format
      # without this, mex compiles with the old non-iterleaved mex API.
      # which may be explicitly specified by the flag '-R2017b'
        export_plain('MEX_R2018','-R2018a');  # used with Makefile
     }
  }

# NB! needed by MatLab runtime library
  $q='ARCH';    if (!$ENV{$q}) { export_plain($q,$mar); }
  $q='MCC_TAG'; if (!$ENV{$q}) {
    # MCC output directory (not mex-extension)
      export_plain($q,'bin'.substr($mxt,3));
  }

  print "\n";

  if ($MEX) {
     export_f('MEX',  $MEX);
     export_f('MYMEX',$MEX);
  }
  export_f('MCC',$MCC);

  $q=export_chk('MCC_BIN',"$MCC/$mar");

  if ($nerr) { exit 1; }

# -------------------------------------------------------------------- #
# check LD_LIBRARY_PATH

  my (@ll,$LD,$cmd);

  if ($ldflag) {
     my $ld='LD_LIBRARY_PATH';
     $cmd="$matlab -nodesktop -nodisplay -n";
     my @ML=`$cmd`; my $e=$?; @ll=grep(/$ld/,@ML);
     if ($e || @ll!=1) { wbdie(
        "failed to determine matlab default env $ld\ncmd: $cmd (e=$e)"); }

   # reduce LD_LIBRARY_PATH to just matlab paths, prepended by '.'
     $_=shift(@ll); s/.*$ld[\s=]*//; chomp; @ll=('.');
     foreach (split(/:/,$_)) {
      # NB! -LD option is relevant in the deployment of mcc sources only
        if (/$MLR.*$mar/ && !/cef|java/) { push(@ll,$_);  }
     }; $LD=join(':',sort @ll);

     if ($ldflag<=1) {
        s/:/\n     : /g; s/$MLR/\$MLR/g;
        print STDERR "\n  # using matlab's default $ld\n\n       ",$_,"\n\n";
     }
     else { print("\n");
        export_plain($ld,$LD);
      # Preload glibc_shim in case of RHEL7 variants
        my $ldd='/usr/bin/ldd';
        if (-f $ldd) {
           foreach (`ldd --version`) {
              if (/GNU libc.*(2\.17)/) { pre_load_shim($1); last; }
           }
        }; print("\n");
     }
  }

# -------------------------------------------------------------------- #
# check LMA, RC_STORE, RC_SYNC

  my $LMA=$ENV{LMA}; if (!$LMA) { $LMA=''; }

  foreach $q ('LMA','RC_STORE','RC_SYNC') { my $D=$ENV{$q};
     if (!$D) { ++$nerr; wbdie("missing env directory $q"); }
     elsif (!-d $D) {
        if ($D=~/(.*)\/(RCStore|RCSync|Data)\/?$/) { my $D1=$1;
           if (-d $D1 || $D1=~$LMA) {
              wblog("WRN mkdir %-8s (%s)",$q, $q ne 'LMA'? repHome($D):$D);
              if (!$tflag) { rsys("mkdir -p '$D'"); }
              next;
           }
        }
        ++$nerr; wblog("ERR invalid env directory $q=%s",repHome($ENV{$q}));
     }
  }

  if ($nerr) { exit 1; }

# -------------------------------------------------------------------- #
# OMP_NUM_THREADS // WRN! matlab is not always listening to this!
# -> used maxNumCompThreads() within matlab [instead]

  if (($q=env_num_threads())>0) {
     export_plain('NSLOTS',$q);
     my $x=$ENV{ML_USING_PARPOOL}; if ($x && $x>1) { $q/=$x; }
     my $n=$ENV{QSP_NUM_THREADS};
     if ($n) { 
        my $q_=int($q/$n);
        if ($q_>1) { $q=$q_; } else { $n=$q; $q=1; }
     }
     else { $n=1; }

     $n=int($n); $q=int($q);

     export_plain('QSP_NUM_THREADS',$n);

     export_plain('OMP_NUM_THREADS',$q>1? $q:'-unset');
     export_plain('MKL_NUM_THREADS',$q>1? $q:'-unset');

if (0) {
     if ($n>1 && $q>1) {
      # export_plain('OMP_NESTED','true'); // deprecated for OpenMP v5
        my $na=10*($n>1? $n : 1)*($q>1? $q : 1);
      # NB! MAX_ACTIVE_LEVELS includes nesting!
        $na*=32; if ($na<4096) { $na=4096; }
        export_plain('OMP_MAX_ACTIVE_LEVELS',$na);
        export_plain('OMP_DYNAMIC','false');

      # default: true => MKL may use fewer threads than indicated
      # in particular, when within an openMP parallel region
      # MKL uses only 1 thread by default (i.e., no nested parallelism)
      # for nested parallelism, accoding to MKL developer guide (v2021.4)
      # MKL requires: OMP_NESTED=true, OMP_DYNAMIC=MKL_DYNAMIC=false
        export_plain('MKL_DYNAMIC','false');
     }
}
  }

  export_plain('CG_VERBOSE',3);

# if (@DEFS) {
#    $_=join(' -D','',@DEFS); s/\s+//g;
#    export_plain('WBDEFS',$_);
#  # print STDERR "WBDEFS = $_\n";
# }

  set_limit_aux('-n', 4096);

# -------------------------------------------------------------------- #
# check basic environment assumed defined by earlier call to *this
# Wb,Feb11,26

sub check_ENV {

   my ($v,$x,@ll);

   $v='MATLAB_ROOT';
   if (!$ENV{$v}) { push(@ll,"undefined ENV ".sprintf('%-12s (',$v).
      ($matlab? "having $matlab" : "also command matlab undefined").")"); }

   $v='MEX_R2018';
   if (!$ENV{$v}) { push(@ll, "undefined ENV ".sprintf('%-12s ',$v).
      "(auto-defined in system/matlab_setup.pl)"); }

   $v='ARCH';
   if (!$ENV{$v}) { push(@ll, "undefined ENV ".sprintf('%-12s ',$v).
      "(auto-defined in system/matlab_setup.pl)"); }

   if (@ll) {
      select STDERR; my $Lsep='-' x70;
      print "\e[31m  ",join("\n  ",'',@ll),"\e[0m\n";
      $_=$0; s/.*(system)/[QSpace-repository\/] $1/;
      die "  $Lsep\n  Please check QSpace environment (see Docs).\n".
          "  Did you source $_?\n".
          "  $Lsep\n\n";
   };
};

# -------------------------------------------------------------------- #
# LD_PRELOAD shim.so as suggested by mathworks support
# to avoid errors such as: undefined symbol: __cxa_thread_atexit_impl
# outsourced // Wb,May13,22

sub pre_load_shim {

   my $v=shift;

# $MLR/bin/glnxa64/glibc-2.17_shim.so
   my @ll=`ls $MLR/bin/glnxa64/glibc-*_shim.so 2>/dev/null`;
   foreach (@ll) {
      if (/$v/) { chomp;
         export_plain('LD_PRELOAD', $_);
         return;
      }
   }
   if (@ll) {
      print STDERR join("\n  WRN ", '',$0,
     "got mismatch in glibc version (ldd --version => $v)", @ll),"\n";
   }
};

# -------------------------------------------------------------------- #
# plain export_plain(vname,value [,doit])

sub export_plain {

   if (@_<2 || @_>3) { wbdie("export() invalid usage",@_);  }
   my $cflag=(@_==3 && !$_[2] ? 1 : 0);
   my $x="$_[1]"; 

   if ($x=~/^-unset/i) { $x='unset '.$_[0]; }
   elsif ($tflag) {
      if ($x=~/[a-z]/i) { $x=repHome($x); }
      $x=sprintf("%-18s %s\n",$_[0],$x);
   }
   else { $x="export $_[0]=$x"; }

   if ($tflag) {
      print '  ', $cflag ? '#':'-', " $x\n";
      return 0;
   }
   else {
      print $cflag ? '#':' ', " $x\n";
      return ($cflag ? 0 : 1);
   }
};

# -------------------------------------------------------------------- #
# export_chk(vname,value) - only export, if variable is not yet set

sub export_chk {
   if (@_!=2) { wbdie("export_chk() invalid usage",@_);  }
   my $v=$_[0]; my $q=$ENV{$v}; if ($q) 
        { return export_plain($v,$q, 0); }
   else { return export_plain($v,$_[1]); }
};

# -------------------------------------------------------------------- #
# export_f(vname,value) - force export (possibly overwriting ENV)

sub export_f {
   if (@_!=2) { wbdie("export_f() invalid usage",@_);  }
   my $v=$_[0];
   my $q=$ENV{$v}; if ($q && $q eq $_[1]) 
        { return export_plain($v,$q, 0); }
   else { return export_plain($v,$_[1]); }
};

# -------------------------------------------------------------------- #
1;  # keep this

