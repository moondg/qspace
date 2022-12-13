#!/usr/bin/perl
# Usage: eval "$(matlab_setup.pl [opts])"
# 
#    this script double checks current matlab bash environment
#    and generates further relevant bash settings in text form
#    if required for running matlab.
# 
#    the correct matlab and gcc version must already be available
#    e.g. loaded via the module system.
# 
# Options
# 
#   -t   test flag (prints status info to stderr)
#   -q   quiet mode, i.e., suppress certain warnings
#        intended for calling *this in Makefiles
#
#   -LD  enforces buildup of full LD_LIBRARY_PATH
#   -ld  using matlab's default LD_LIBRARY_PATH (once it starts)
# 
# Wb,Jan09,19

  use Cwd 'realpath';

  my $P=realpath($0); $P=~s/\/[^\/]*$//;
  require "$P/plib.pl";

  my $nerr=0; my $nwrn=0; my $vflag=1;
  my ($q,$mver,$gver, $tflag, $ldflag,$clear,$pre,@pre);

  while (defined ($_=shift)) {
     if (/^-[h\?]$/) { die usage($0); }
     elsif (/^-(t)$/i) { $tflag+=($1 eq 'v'?1:2);  }
   # elsif (/^-(v)$/i) { $vflag+=($1 eq 'v'?1:2);  }
     elsif (/^-q$/) { $vflag=0; }
     elsif (/^-(ld)$/i) { $ldflag+=($1 eq 'ld'?1:2);  } # -ld -LD
     elsif (/^--clear$/) { ++$clear; }
     elsif (/^--pre-run$/) { ++$pre; } # Wb,May05,19May04,19
     elsif (/^-\d+[a-z]$/i) { push(@pre,$_); }
     else { 
        die usage(__FILE__,__LINE__,"got extra arguments: ",$_,@ARGV);
     }
  }

  if ($pre) { $_=$pre[0];
     if (@pre==1 && /^-(\d+)([a-z])$/i) { # Wb,May05,19
        my ($y,$v,@ll); $y=$1; $v=lc($2);
        if (length($1)<=2) { $y+=2000; }

        push(@ll,"module load matlab/$y$v");

        if ($y<2014 || $y>2020) { wbdie("invalid matlab version $y$v ($_)"); }
        if ($y<2018)
             { push(@ll,"module load gcc/4.7.4"); }
        else { push(@ll,"module load gcc/6.3.0"); }

        $v=($tflag ? stderr : stdout);
        foreach (@ll) { printf($v "  $_\n"); }
     }
     elsif (@pre) { wbdie("invalid usage (".join(', ',@pre).')'); }
   # elsif ($tflag) { print(stderr "\n  TST got empty --pre-run\n\n"); }

     exit 0;
  }
# else ignore @pre

# $\="\n  ";
  $|=1; # auto flush after every printf

# -------------------------------------------------------------------- #
# check MATLAB_ROOT and consistency with `which matlab`

  my $mlr='MATLAB_ROOT';
  my $MLR=$ENV{$mlr};
  my $matlab=`which matlab 2>/dev/null`;
  my ($gcc);

  if ($clear) {
     my $ML=$ENV{MYMATLAB};
     if ($MLR) { print join("\n", remove_paths('--bash','-q',$ML,
        'PATH','LD_LIBRARY_PATH','XAPPLRESDIR'));
     }
     exit 0;
  }

  if (!$MLR || !$matlab) { ++$nerr; wblog(
     "ERR env $mlr not defined (hint: load module for matlab?)"); }
  else { chomp($matlab); 
     if (!-d $MLR.'/bin/') { ++$nerr; wblog(
        "ERR invalid env $mlr=%s (missing subdirectory ./bin",$MLR); }
     elsif ($matlab!~/$MLR/) { ++$nerr; wblog(
        "ERR invalid env $mlr=%s (%s)",$MLR,$matlab); }
     elsif ($MLR!~/matlab/i) { ++$nwrn; wblog(
        "WRN unexpected value for env $mlr=%s",$MLR); }
     if ($matlab=~/[_\/rR]?(\d{4}\w?)[_\.\/]/) { $mver=$1; }
     # e.g. mac path /Applications/MATLAB_R2018a.app
     # e.g. LRZ: matlab/R2018a_Update6
  }

# check gcc version and consistency with MATLAB_ROOT

  if (islinux()) {
     $gcc=`which gcc 2>/dev/null`;

     if ($gcc) { chomp($gcc);
        foreach (`$gcc --version 2>/dev/null`) {
        if (/gcc[^\d]* (\d+\.\d+.\d+)[^\w]/) { $gver=$1; }}
     }
     else { ++$nerr;
     wblog("ERR gcc not found (hint: load module for gcc?)"); }
  }
  elsif (ismac()) {
     $gcc=`which gcc 2>/dev/null`; # mac also writes to stderr here

     if ($gcc) { chomp($gcc);
        foreach (`$gcc --version 2>/dev/null`) {
           if (/Apple.*clang/) {
              foreach (`xcodebuild -version`) {
              if (/Xcode *(\d+[\.\d]+)/) { $gver='clang/'.$1; last; }}
              last;
           }
        }
     }
     else { ++$nerr;
     wblog("ERR gcc not found (hint: load module for gcc?)"); }
  }
  else { wbdie("got unknown OS %s !?",`uname -s`); }

  if (!$nerr) {
  if (!$mver) { ++$nwrn; wblog(
     "WRN failed to derive matlab version from env $mlr"); }
  elsif (!$gver) { ++$nerr; wblog(
     "ERR failed to derive gcc version (having $gcc)"); }
  elsif (islinux()) {
     if ($mver=~/2016a/) { if ($gver!~/^4\.7\./) { ++$nerr; wblog(
        "ERR unsupported gcc version $gver for matlab/$mver"); }}
     elsif ($mver=~/2016b/) { if ($gver!~/^4\.8\./) { ++$nerr; wblog(
        "ERR unsupported gcc version $gver for matlab/$mver"); }}
     elsif ($mver=~/201[89]|2020a/) { if ($gver!~/^6\.3\./) { ++$nerr; wblog(
        "ERR unsupported gcc version $gver for matlab/$mver"); }}
     elsif ($mver=~/2020b/) { if ($gver!~/^8\.[3-5]\./) { ++$nerr; wblog(
        "ERR unsupported gcc version $gver for matlab/$mver"); }}
     elsif ($mver=~/2022b/) { if ($gver!~/^([789]|10)\./) { ++$nerr; wblog(
        "ERR unsupported gcc version $gver for matlab/$mver"); }}
     else { ++$nwrn; if ($vflag) { wblog(
        "WRN check valid gcc version ($gver) for matlab/$mver"); 
     }}
  }
  elsif (ismac()) { my $wrn=0;
     if ($mver=~/201[89]|2020/) {
      # matlab/2020b supports Xcode 10.x and 11.x
        if ($gver!~/\/1[01]\.[0-7]/) { ++$wrn; }
     }
     elsif ($mver=~/2022/) {
      # matlab/2022a supports Xcode 12.x and 13.x
        if ($gver!~/\/1[23]\.[0-9]/) { ++$wrn; }
     }
     else { ++$wrn; }

     if ($wrn) { $nwrn+=$wrn; if ($vflag) {
        wblog("WRN unsupported C++ compiler version $gver for matlab/$mver");
     }}
  }}

  if ($nerr) { exit 1; }

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

  my $mar=getARCH(); my $mxt=getMEXEXT(); # plib.pl
  if (!$mxt) { ++$nerr; }

  if ($tflag) {
     printf(STDERR "\n  %-20s %s\n",$MM,repHome($ML));
     if ($mver) { printf(STDERR "  %-20s %-12s # %s\n",'matlab version',$mver,$matlab); }
     if ($gver) { printf(STDERR "  %-20s %-12s # %s\n\n",'gcc version',$gver,$gcc); }
  }

  if ($mver) {
     $q=$mver; $q=~s/^20//;
     if ($q=~s/([a-z])$//i) {
        my $a=lc($1); $q+=(ord($a)-ord('a')+1)/10;
     }

     export_plain('MATLAB_VERNUM',$q);     # numeric value for conditionals
     export_plain('MATLAB_VERSION',$mver); # to be used with 'module load ...'

     if ($q>=18) {
      # required for MEX files to avoid compiler errors with functions
      # not found such as: mxComplexDouble, mxGetComplexDoubles, etc.
      # NB! matlab>=2018 has interleaved complex format
      # without this, mex compiles with the old non-iterleaved mex API.
      # which may be explicitly specified by the flag '-R2017b'
        export_plain('MEX_R2018','-R2018a'); # used with Makefile
     }
  }

  if ($gver) {
   # just a safeguard to tell Makefile that gcc version has been checked
     export_plain('ML_GCC_VERSION',$gver);
  }

# NB! needed by MatLab runtime library
  $q='ARCH';    if (!$ENV{$q}) { export_plain($q,$mar); }
  $q='MCC_TAG'; if (!$ENV{$q}) {
    # MCC output directory (not mex-extension)
      export_plain($q,'bin'.substr($mxt,3));
  }

  print "\n";

  if ($MEX) {
   # MYMEX is my MEX root directory (used by MEX/make)
     export_frc('MEX',  $MEX);
     export_frc('MYMEX',$MEX);
  }
  export_frc('MCC',$MCC);

  $q=export_chk('MCC_BIN',"$MCC/$mcd");
  # if (!$q) { $q='MCC_BIN'; if (!-d $ENV{$q}) { ++$nerr; wblog(
  # "ERR invalid directory env $q=$ENV{$q}"); }}

  if ($nerr) { exit 1; }

# -------------------------------------------------------------------- #
# check LD_LIBRARY_PATH

  my (@ll,$LD,$cmd);

  if ($ldflag) {
     my $ld='LD_LIBRARY_PATH';
   # NB! this prepends to the current LD_LIBRARY_PATH -> unset
   # $cmd="unset $ld ARCH; $matlab -nodesktop -nodisplay -n";
   # inherit existing LD_LIBRARY_PATH e.g. with entries on gcc 
     $cmd="$matlab -nodesktop -nodisplay -n";
     my @ML=`$cmd`; my $e=$?; @ll=grep(/$ld/,@ML);
     if ($e || @ll!=1) { wbdie(
        "failed to determine matlab default env $ld\ncmd: $cmd (e=$e)"); }

   # reduce LD_LIBRARY_PATH to just matlab paths, prepended by '.'
     $_=shift(@ll); s/.*$ld[\s=]*//; chomp; @ll=('.');
     foreach (split(/:/,$_)) {
      # NB! -LD option is relevant in the deployment of mcc sources only
      # e.g. as also generated by mcc in run_*.sh // Wb,Nov24,21
        if (/$MLR.*$mar/ && !/cef|java/) { push(@ll,$_);  }
     }; $LD=join(':',sort @ll); # mcc run_*.sh also has /sys/ last

     if ($ldflag<=1) {
        s/:/\n     : /g; s/$MLR/\$MLR/g;
        print STDERR "\n  # using matlab's default $ld\n\n       ",$_,"\n\n";
     }
     else { print("\n");
        export_plain($ld,$LD);
      # adapted from mcc compiled run_*.sh script // Wb,Nov24,21
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
        if ($q_>1) { $q=$q_; } else { $n=$q; $q=1; } # $n=1;
     }
     else { $n=1; }

     $n=int($n); $q=int($q);

     export_plain('QSP_NUM_THREADS',$n);

   # NB! MKL automatically goes serial if in parallel environement
     export_plain('OMP_NUM_THREADS',$q>1? $q:'-unset');
     export_plain('MKL_NUM_THREADS',$q>1? $q:'-unset');

if (0) { # Wb,Nov25,21
     if ($n>1 && $q>1) {
      # export_plain('OMP_NESTED','true'); // deprecated for OpenMP v5
        my $na=10*($n>1? $n : 1)*($q>1? $q : 1);         # Wb,Nov24,21
      # NB! MAX_ACTIVE_LEVELS includes nesting!
      # --> significantly increase $na, otherwise job runs mostly serial!
        $na*=32; if ($na<4096) { $na=4096; }    # default is MAX_INT!
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

  if (@DEFS) {
     $_=join(' -D','',@DEFS); s/\s+//g;
     export_plain('WBDEFS',$_);
   # print STDERR "WBDEFS = $_\n";
  }

# NB! increase `ulimit -n' to at least 4096 // Wb,Jul01,19
  set_limit_aux('-n', 4096);
# set_limit_aux('-s',65536);

# -------------------------------------------------------------------- #
# LD_PRELOAD shim.so as suggested by mathworks support
# to avoid errors such as: undefined symbol: __cxa_thread_atexit_impl
# outsourced // Wb,May13,22

sub pre_load_shim {
  
   my $v=shift; # ldd glibc version

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
# export_frc(vname,value) - force export (possibly overwriting ENV)

sub export_frc {
   if (@_!=2) { wbdie("export_frc() invalid usage",@_);  }
   my $v=$_[0]; my $q=$ENV{$v}; if ($q && $q eq $_[1]) 
        { return export_plain($v,$q, 0); } # comment out
   else { return export_plain($v,$_[1]); }
};

# -------------------------------------------------------------------- #
1; # keep this

