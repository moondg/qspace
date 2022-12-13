#!/usr/bin/perl
# NAME [opts] [explicit files to include]
#
#    Synchronize this repo (QSpace-3.2) with MEX sources.
#
#    May explicitly also specify files to be included
#    e.g. because skipped since time step was preserved
#    while nevertheless changing a file. E.g. search for
#    find . -maxdepth 1 -type f -mtime +100 -ctime -30
#
# Options
#
#    --fix     fix time stamps if files are up to date otherwise
#    --FIX     check ALL files whether they are up-to-date
#              before and after decommenting
# Examples
#
#    rsync_v4s.pl -F -V --fpat '\.[ch]'   # selects C++ files only
#    rsync_v4s.pl -F -V --fpat '\.m'      # selects matlab files only
#
#  # just fixes time stamps if files are uptodate otherwise
#    rsync_v4s.pl -t --Fix 
#
# Wb,Aug19,18

  require $ENV{HOME}.'/bin/plib.pl';
  use File::stat;

  my ($tflag,$vflag,$force,$fpat,@FX,$FXp,@task); $force=0; $vflag=1;
  my ($fixme,$fixed1,$fixed2,$fix1,$fix2,$fix3,$fixed3,@FX0,@FX2);

  while (defined ($_=shift)) {
     if (/^-[h\?]$/) { die usage($0); }
     elsif (/^-(v)$/i) { $vflag+=($1 eq 'v' ? 1 : 2); } # -v,-V
     elsif (/^-q$/i) { $vflag=0; }
     elsif (/^-(t)$/i) { $tflag+=($1 eq 't' ? 1 : 2); } # -t,-T
     elsif (/^-(f)$/i) { $force+=($1 eq 'f' ? 1 : 2); } # -f,-F
     elsif ($_ eq '--fpat') { $fpat=shift; }
     elsif (/^--(fix|check)/i) { push(@task,$_); }
     elsif (!/^-$/) { push(@FX,$_); }
     else { die usage(__FILE__,__LINE__,$_,@ARGV); }
  }

  my $me=myname();

# my $DST=$ENV{HOME}.'/BitBucket/qspace-3.2'; # Wb,Jan12,19
  my $DST=realpath($0); $DST=~s/\/[^\/]*$//;  # qs40

  my $ML=$ENV{HOME}.'/Matlab';
  my $MEX=$ML.'/MEX'; # $SRC

  if (!-d $DST) { wbdie("invalid DST directory %s",repHome($DST)); }
  if (!-d $MEX) { wbdie("invalid MEX directory %s",repHome($MEX)); }

  foreach (@task) {
     if (/^--(fix\w*)$/i) {
        my $q=$1; my @q=split(//,$q); $fixme|=1; # NB! used as bit pattern
        if ($q[0]         eq 'F') { $fixme|= 2; }
        if ($q[1]         eq 'I') { $fixme|= 4; }
        if ($q[2]         eq 'X') { $fixme|= 8; }
        if (@q>3 && $q[3] eq 'M') { $fixme|=16; }
        if (@q>4 && $q[4] eq 'E') { $fixme|=32; }
     }
     elsif (/^--check-git-files$/) { check_git_files(); exit $?; }
     else { wbdie("invalid task $_"); }
  }

# my $cp='/bin/cp --preserve=timestamps'; # -u  // unix!
  my $cp=cp_sys();

  chdir($DST);
# printf("\n  Current directory: %s",`pwd`);

  printf "\n  %s\n      %-12s : %s\n  --> %-12s : %s (pwd)\n\n",
     $me,'source',repHome($MEX),'destination',repHome($DST);

  if ($fpat) { print "  using file pattern '$fpat'\n\n"; }

  if (@FX) { $FXp=join('|',@FX);
     printf "  Using: $FXp\n\n";
  }

# if ($tflag) { exit 1; }

# -------------------------------------------------------------------- #
# usage: $nf = rsync_dir(target_dir, source_dir);

  my %GFS; # global file listing of source (MEX)

# NB! decomment_source.pl checks for comments of the type
#     /* LICENSE/license.txt [; Class: ...] */
# in the very header lines of each file, and then inserts the
# corresponding license boiler plate, with the path specified
# relative to source file as is or relative to the source
# file's parent directory // Wb,Aug31,22

  my $Lver='v4.0 pre-release';

# NB! this preserves time stamps
  my $decomment="decomment_source.pl --Lver='$Lver'";

  if (!$tflag) {
     rsys("$decomment --check-Lver LICENSE.txt"); }

  $decomment.=' --nlog 256'; # default: 96

# -------------------------------------------------------------------- #

# my @FF=LS_files('.');
# print join("\n",@FF); printf("\n  %g files\n",$#FF+1);

  my $nf=0;

  $nf+=rsync_dir('Source',$MEX);
  $nf+=rsync_dir('NRG',   $ML.'/NRG',$MEX);
  $nf+=rsync_dir('DMRG',  $MEX,$ML.'/MPS',$ML.'/util');
  $nf+=rsync_dir('lib',   $ML.'/lib', $ML.'/libx');
  $nf+=rsync_dir('bin',   $MEX);
  $nf+=rsync_dir('tensor',$MEX,$ML.'/MPS',$ML.'/QCI',$ML.'/NRG');
  $nf+=rsync_dir('Class', $ML.'/Class');

  $nf+=rsync_dir('setup', $MEX,$ML.'/NRG',$ML.'/Phys'); # Wb,Nov28,18

  $nf+=rsync_fls('bin' ,'-h', $MEX.'/bin/*.m' );  # Wb,Mar27,19
  $nf+=rsync_fls('util','-h', $MEX.'/util/*.m');  # Wb,Mar27,19

# $nf+=rsync_fls('Docu',$ML.'/qspace-v4.0/README*');
# $nf+=rsync_fls('Docu',$ML.'/qspace-v4.0/Docu/*pdf');
# $nf+=rsync_fls('MCC', $ML.'/qspace-v4.0/MCC/[rM]*'); # runMCC, Makefile

  if ($vflag>1 || $fixme) { printf "\n  Total of %g files checked.\n",$nf; }

  check_readable($DST);
  # { my @ll=`cto.pl -q -d mp32x`;
  #   if (@ll) { check_readable($ll[0]); }
  # }

  if ($fixme) { $fix1+=0; $fix2+=0; $fix3+=0;
     printf "\n  check_and_fix_diff() checked %g files\n".
     "  %4g files actually up to date (yet with different time stamps)\n".
     "  %4g actually not different (BEFORE decomment)\n".
     "  %4g files really different\n",
     $fix1+$fix2+$fix3, $fix1, $fix2, $fix3;

     if ($fixed1 ||  $fixed2 ||  $fixed3) {
         $fixed1+=0; $fixed2+=0; $fixed3+=0; printf "\n".
     "  %4g mtime stamps fixed\n".
     "  %4g mtime stamps fixed (while also using decomment)\n",
     "  %4g files fully fixed\n", $fixed1, $fixed2, $fixed3;
     }

     if (@FX0) {
        print "\nFX0: ",join(' ',@FX0),"\n";
        print "\nFX2: ",join(' ',@FX2),"\n";
     }
  }

# end of main // done

# -------------------------------------------------------------------- #
# usage: LS_files($vflag,$dir) // like `ls'
# where vflag indicates verbose mode (0 = quiet)

sub LS_files { # recursive

  if (@_!=2 || !-d $_[1]) { wbdie("invalid usage (%s)",join(', ',@_)); }
  my (@ff,@FF,@fd,$F); my ($v,$D)=@_; # $_[0], $_[1]

  opendir(DH,$D); @fd=readdir(DH); closedir(DH);

  foreach my $f (@fd) { $F=$D.'/'.$f;
      if ($f=~/^\./) { next; } # e.g. '.' and '..'

      if (-d $F) {
         if ($f=~/(Archive|master|dev|tmp)/i || $f=~/\d{6}/) { next; }}
      elsif ($f=~/(README|test|html)/i || $f=~/\d{6}/) { next; }

      if ($f=~/_R20\d\d\w?$/) { next; }
      if ($f=~/\b(cto|dir2|repHome|hostid)\b/i) { next; }
       # WRN! wbvector shares wbve'cto'r
      if ($f=~/(startup|finish|mlinfo|smaxis)/i) { 
         if ($f!~/startup_(info|numthr)/) { next; }
      }
      if ($f=~/(Makefile|\.mex|\.o|\.pl$)/i) { next; }
       # startup.m startup_aux.m
       # libx/hostid -> calls hostid.pl
      if ($f=~/(make.m|vi.m|setRCStore)/i) { next; }
      if ($f=~/rnrg|(plot|run)DMRG|tst_(openWilson|Hubbard)/) { next; }
       # WRN! include fgrNRG.cc -> case sensitive! // Wb,Feb14,19
       # routines based on (tst|run)_Hamilton1D
      if    (-f $F) { push(@ff,$F); }
      elsif (-d $F) { push(@FF,LS_files($v,$F)); }
      elsif ($v) { printf(STDERR
        "  WRN skipping non-regular file '%s'\n",repHome($F));
      }
  }

  if ($fpat) { @ff=grep(/$fpat/,@ff); }

  @ff=sort(@ff); # let files come first (higher priority)
  push(@ff,@FF);

  return @ff;
};

# -------------------------------------------------------------------- #

sub rsync_dir {

  if (@_<2) { wbdie("invalid usage (%g args)",$#_+1); }
  foreach (@_) { if (!-d $_) {
     wbdie("invalid directory '%s'",repHome($_)); }
  }

  my ($D,$f,@f,$ns,%FS,@FF,$ok,$nup,@fx,$xp);

  $D=shift;
  if ($D=~/Class/)
       { $xp='Class\/(.*)$'; } # $xp='\/\@(.*)$';
  else { $xp='\/([^\/]+)$'; }

  foreach (@_) { 
     if (!exists $GFS{$_}) {
        if ($vflag>1) { printf "  LS: listing files in %s\n",repHome($_); }
        push(@{$GFS{$_}},LS_files(0,$_));
     }
     $ns+=($#{$GFS{$_}}+1);
   # generate hash for faster search
     foreach (@{$GFS{$_}}) {
        if (/$xp/) { push(@{$FS{$1}},$_); }
        else { wblog("WRN $_ !?"); }
     }
  }

  @FF=LS_files(1,$D);
  if (!@FF && $fpat) { return $#FF+1; }

  if ($vflag>1) {
     my $s=repHome($D); my $l=64-length($s);
     printf "\n# %s %s ---\n", '-' x $l,$s;
     printf("  checking %g files\n",$#FF+1);
     foreach (@_) {
        printf "   *  %s (%g files)\n",repHome($_),$#{$GFS{$_}}+1;
     }
     if (@_>1) { printf "  ... %g files total\n",$ns; }
  }

  foreach (@FF) {
     if (/$xp/) { $f=$1;
        if (defined $FS{$f}) { @f=@{$FS{$f}};
           if (@f==1) { ++$ok; my $upd=0;
              my ($f0,$f2)=($f[0],$_); 
              if (!-f $f0){ wbdie("invalid file %s",repHome($f0)); }
              if (!-f $f2){ wbdie("invalid file %s",repHome($f2)); }

              my $s0=stat($f0); 
              my $s2=stat($f2);

              if ($fixme) { if ($fixme>14 || $s2->mtime != $s0->mtime) {
                 if (check_and_fix_diff($f0,$f2)<0) { # --FIX => fixme=15
                    $s0=stat($f0); $s2=stat($f2);
                 }
              }}

              if ($s2->mtime < $s0->mtime || (@FX && /$FXp/)) { ++$nup; ++$upd; 
               # NB! may have updated file, but preserved time stamp!
               # -> rather specify manually at command line when calling *this
               #  $s2->mtime < $s0->ctime ||
                 if ($vflag>2 || $tflag) { printf
                    "\e[32;1m   +  %-30s -> %s\e[0m\n",
                    repHome($f0),repHome($f2);
                 }
              }
              elsif ($s2->mtime > $s0->mtime) { ++$upd;
                 wblog("WRN got newer file at destination !?");
                 print repHome(ls($f0,$f2));
              }

              if (($upd || $force>1) && !$tflag) {
                 rsys("$cp '$f0' '$f2'"); if ($f2=~/\.[chm]+$/) {
                 rsys($decomment.($vflag>2 ? '':' -q')." '$f2'"); }
              }
           }
           else {
              wblog("WRN got %g matches for '$f'",$#f+1);
              print join("\n",@f,'');
           }
        }
        else {
         # wblog("WRN no entry for '$f' !?");
           push(@fx,$_); 
        }
     }
     else { wblog("WRN $_ !?"); }
  }

  if ($vflag>1) {
     printf "  ok. %g/%g files uniquely identified\n",$ok,$#FF+1;
     if ($nup) { printf "  UPD %g/%g to update\n",$nup,$#FF+1; }
  }
  if (@fx) {
     wblog("ERR failed to identify %g file(s)",$#fx+1);
     foreach (@fx) { print '     ',$_,"\n"; }
  }

  return $#FF+1;
}

# -------------------------------------------------------------------- #
# usage: $nf = rsync_fls(target_dir, source_fls);
# with source_fls patterns written in rsync format
# Wb,Nov24,18

sub rsync_fls {
  if (@_<2) { wbdie("invalid usage (%g args)",$#_+1); }

  my $TAR=shift(@_);
  if (!-d $TAR) { wbdie("invalid TAR directory '%s'",repHome($TAR)); }

  my $rsync="/usr/bin/rsync -rptv";
  if ($tflag) { $rsync.=" --dry-run"; } # -n

  my ($cmd,@ll,@ff,$SRC,$excl); my $excl='';

  foreach my $fls (@_) {
     if ($fls eq '-h') { 
      # exclude certain help files in user space MEX // Wb,Mar09,20
        if (!$excl) { $excl='--exclude tstmex* --exclude wbtic*'; }
        else { wbdie("invalidag usage (already got $excl)"); }
      # wblog("TST $excl"); 
        next;
     }
     $SRC=$fls; $SRC=~s/[^\/]*$//;
     if (!-d $SRC) { wbdie("invalid SRC directory '%s'",repHome($SRC)); }

     $cmd="$rsync $excl $fls $TAR/";
     @ll=`$cmd`; @ff=();
     if ($?) { printf "\n  ERR cmd: %s\n\n",$cmd; wbdie(""); }
     foreach (@ll) {
         if (/^sending incremental.*list/) { next; }
         if (/^building file list.*done/) { next; }
         if (/^sent.*bytes.*received/) { next; }
         if (/^total.*size.*speedup/) { next; }
         if (/^\s*$/) { next; }
         push(@ff,$_);
     }
     if (@ff) {
        printf "\n   %s%g file(s) %-30s --> %s\n", $tflag ? 'TST ':'',
        $#ff+1, repHome($SRC), repHome($TAR);
        print join('      ','',@ff);
     }
  }
}

# -------------------------------------------------------------------- #
# fix mp32 files (originally written when switching LMU->BNL repo)
# Wb,Mar09,20

sub check_and_fix_diff {

   my $isdiff=0; if (@_!=2) { wbdie("invalid usage"); }

   my ($f0,$f2)=@_; my $fx=$f2;
   if (!-f $f0) { wbdie("invalid file $f0"); }
   if (!-f $f2) { wbdie("invalid file $f2"); }
   if ($fx=~s/([^\.])\.(\w{1,3})$/$1-tmp.$2/) { my $x=$2;
      if ($x!~/[chm]+$/) {
         if ($f0=~/ClebschGordan_Arne/) 
              { wblog("TST skipping $f0"); return $isdiff; }
         else { wbdie("unexpected file name $f2"); }
      }
   }
   else {  wbdie("unexpected file name $f2"); }
   if (-f $fx) { wbdie("file $fx already exists !?"); }

   my $s0=stat($f0);
   my $s2=stat($f2);

   rsys("$cp '$f0' '$fx'");
   rsys($decomment.($vflag>2 ? '':' -Q')." '$fx'");
   system("cmp -s '$fx' '$f2'");

   if (!$?) {
    # file f2 already up to date -> just fix time stamp if necessary
      if ($s0->mtime != $s2->mtime) { print "\n"; ++$fix1;
         wblog("OK! file up to date, but differs by time stamp only");
         if (!$tflag || $fixme>1) { $isdiff=-1; ++$fixed1; # --Fix
            if ($s0->mtime < $s2->mtime) # set to older time stamp
                 { rsys("touch -r '$f0' '$f2'"); }
            else { rsys("touch -r '$f2' '$f0'"); }
         }
         ll('-p',$f0,$f2,$fx);
      }
      rsys("rm '$fx'"); return $isdiff;
   }; print "\n";

   system("vimdiff -c 'syntax off' '$fx' '$f2'");
   sleep(1); # print " press key to continue ... "; my $i=<STDIN>;

 # keep comments with certain files
   my $keepcmts=grep(/tst_tdDMRG|tst_Hamilton1D/,$f0);

   system("cmp -s '$f0' '$f2'");
   if (!$?) {
      if ($keepcmts) { rsys("rm '$fx'");
         if ($s0->mtime!=$s2->mtime) { ++$fix1;
            wblog("WRN! files are actually the same (keeping comments)");
            ll('-p',$f0,$f2);

            if (!$tflag || $fixme>1) { $isdiff=-3; ++$fixed1;
               if ($s0->mtime < $s2->mtime) # set to older time stamp
                    { rsys("touch -r '$f0' '$f2'"); }
               else { rsys("touch -r '$f2' '$f0'"); }
            }
         }
      }
      else {
         wblog("WRN! files are actually the same (but needs decomment)");
         ll('-p',$f0,$f2,$fx); ++$fix2;

         if (($fixme&6)==6 && ($s0->mtime != $s2->mtime)) {
         #  --FIx => bits 2 & 3 => 2+4 = 6
            if ($s0->mtime < $s2->mtime) # set to older time stamp
                 { rsys("touch -r '$f0' '$fx'"); }
            else { rsys("touch -r '$f2' '$f0'"); }
            rsys("mv '$fx' '$f2'"); $isdiff=-2; ++$fixed2;
         }
         else { rsys("rm '$fx'"); }
      }
      return $isdiff;
   }

   wblog("NB! files are actually different"); ++$fix3; $isdiff=1; 
   ll('-p',$f0,$f2,$fx);

   if (($fixme&48)==48) { # --fixME => bits 4 & 5 => 16+32=48
    # update file in mp32 (time stamp in fx is preserved / up to date)
      if ($keepcmts) { rsys("rm '$fx'");
         system("vimdiff -c 'syntax off' '$f0' '$f2'");
         printf "\n  update $f2 <- $f0 ? [y]n  ";
         chomp (my $q=<STDIN>); $q=~s/\s+//g; if (!length($q)) { $q='y'; }
         if ($q!~/^[n0]/i) { print "  "; 
            rsys("$cp -v '$f0' '$f2'"); $isdiff=-4; ++$fixed3;
         }
      }
      else {
         system("vimdiff -c 'syntax off' '$fx' '$f2'");
         printf "\n  overwrite $fx -> $f2 ? [y]n  ";
         chomp (my $q=<STDIN>); $q=~s/\s+//g; if (!length($q)) { $q='y'; }
         if ($q!~/^[n0]/i) { print "  "; 
                rsys("mv -v '$fx' '$f2'"); $isdiff=-4; ++$fixed3; }
         else { rsys("rm '$fx'"); }
      }
   }
   else { rsys("rm '$fx'"); }

   if ($isdiff>0) { push(@FX0,$f0); push(@FX2,$f2); }

   return $isdiff;
};

# -------------------------------------------------------------------- #
# check whether all files in DST are readable // safeguard
# Wb,Oct06,21

sub check_readable {
   
   my (@dd,$cmd,$fargs,@ll);
   foreach (@_) {
      if (-d $_) { push(@dd,$_); }
      else { wbdie("invalid directory '$_'"); }
   }
   if (!@dd) { push(@dd,'.'); }

   $fargs='! -perm '.(ismac()? '+':'/').'044';

   foreach (@dd) {
      if ($vflag>1) { wblog(" *  %s  %s",$_,$fargs); }
      push(@ll,`find '$_' $fargs`);
   }

   if (@ll) { $_='not readable by group';
      if (@ll!=1)
           { wblog("WRN found %g files that are $_",$#ll+1); }
      else { wblog("WRN one file $_"); }

      foreach (@ll) { s/$ENV{HOME}/~/; $_='   '.$_; }
      if (@ll<12) { print "\n",@ll; }
      else {
         print "\n",@ll[0 .. 3],"\n    : \n",$ll[$#ll-3 .. $#ll];
      }
   }
}

# -------------------------------------------------------------------- #
# Wb,Sep30,22

sub check_git_files {

   my ($q,@x,@f0,@f2,%C0,%C2,@C0,@C2);

   chdir($ML) || wbdie("invalid directory $ML"); 
   @f0=`git ls-files`;
   foreach (@f0) { if (/(^|\/)(Class\/.*)/) { ++$C0{$2}; }}
   @C0=sort keys %C0;

   chdir($DST) || wbdie("invalid directory $DST"); 
   @f2=`git ls-files`;
   foreach (@f2) { if (/(^|\/)(Class\/.*)/) { ++$C2{$2}; }}
   @C2=sort keys %C2;

   if ($vflag) {
      printf "\n%6d / %4d files at SRC/Class (%s)\n",
         $#C0+1,$#f0+1,repHome($ML); 
      printf   "%6d / %4d files at DST/Class (%s)\n\n",
         $#C2+1,$#f2+1,repHome($DST); 
   }

   foreach (@C0) { if (!$C2{$_}) {
       if (/\/(tst_|Lab\/)/) { next; }
       if (/\.jpg$/i) { next; }
       push(@x,$_);
   }}
   if (@x) {
      printf "  %d Class git-files not in DST:\n",$#x+1;
      print join("\n     ",'',@x),"\n\n";
   }
   else { print "  all Class git-files in DST\n"; }

   @x=();
   foreach (@C2) { if (!$C0{$_}) { push(@x,$_); }}
   if (@x) {
      printf "\n  %d Class git-files not in SRC:\n",$#x+1;
   }
   else { print "  all Class git-files in SRC\n"; }

}

# -------------------------------------------------------------------- #

