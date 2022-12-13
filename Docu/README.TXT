
# NRGWilsonCG + fdmNRG_QS update
# Wb,Feb09,21

  Added support for complex NRG as well as possibly complex ff
  cross-couplings, e.g. in the presence of a non-diagonal hybridization
  function. Subsequently, this also results in complex (off-diagona)
  dynamical correlation functions in fdmNRG_QS, such that the returned
  a0 therein is complex. fdmNRG_QS listens to the Inrg.paras.complex
  flag as set by NRGWilsonQS based on its imput data.

# ==================================================================== #
# MPSPACK / QSpace v3.2.0 (matlab/2018b) released Wb,May18,20          #
# ==================================================================== #

  QSP_NUM_THREADS Parallelization // Wb, Mar-May 2020

     Frequently when dealing with higher-rang tensors, 
     contractions etc. need to work through (very) long lists,
     while individual contractions actually deal with smallish
     tensors. This makes QSP_NUM_THREADS an attractive option.

     Therefore I have put significant effort into improving
     on the QSP_NUM_THREADS parallelization [note, however,
     that due the MKL bug described below [MKL_BUG_OMP] either
     QSP_NUM_THREADS>1 *xor* MKL_NUM_THREADS>1 shall be used,
     i.e., they need to be used exclusively.

     A significant number of race conditions and dead locks
     due to QSP_NUM_THREADS have been removed making extensive
     use of internal locks for shared memory applications.
     These are based on standard OMP omp_nest_lock_t (mutex).
     They force threads to wait until certain RC_STORE related
     objects again become unlocked after updates.

  Data format for the outer multiplicity coefficients cgw and cgt
  in info.cgr has been changed from MPFR to plain double // Wb,May 2020

     Since for any QSpace (say A) in practice, A.info.cgr.cgw
     and A.info.cgr.cgt are used together with A.data which maintains
     double precision (or complex double), throughout, it made
     little sense to keep higher precision there.

     Therefore, while the X-Store (for the X-symbols)
     remains in higher precision MPFR format, once computed
     the X-symbols as well as the cgw and cgt coefficients
     are loaded into memory in standard double format
     as this can significantly speed up numerical operations
     (it had turned out, that MPFR in A.info.cgr.cgw could
     become the numerical bottle neck if one deals with
     smaller A.data blocks but very many of them).

  Other changes

   - every QSpace mex routine now supports the options --version
     which shows matlab, compiler, and QSpace versions, etc.

   - eigQS and svdQS now also parallelize using QSP_NUM_THREADS

   - bug fixes related to parallelization and race conditions

# ==================================================================== #
# MPSPACK / QSpace v3.1 (matlab/2018b)                                 #
# this version is officially retired with matlab/2018b                 #
# effective 05/18/2020 // AW                                           #
# ==================================================================== #

# see follow up / preliminary conclusions drawn out of discussion
# with Matan Lotem in discussion on nested parallelization below;
# search for UPDATE_NESTED // Wb,Jul24,19

# -------------------------------------------------------------------- #
# [LMU specific] Fix on MPFR library version (mismatch)
# Wb,Jul15,19

  Workstations and older cluster machines are configured with
  Ubuntu 16.04, whereas th-cl-naples* nodes now have Ubuntu 18.04

  QSpace is linked against -lmpfr which as it turns out
  points to /usr/lib/x86_64-linux-gnu/libmpfr.so.4
  the latter, however, are no longer available e.g. on the
  th-cl-naples-* nodes.

  Two possible workarounds:

# Work around #1 (quick and dirty)
  In your cluster submission script, since '.' is on the
  LD_LIBRARY_PATH by the mcc-generated run_*.sh scripts,
  the following provides (a valid instantation of) libmpfr.so.4
  for the mex binaries:

  $ rsync -l $SRC:/usr/lib/x86_64-linux-gnu/libmpfr.* $DST

    SRC = your favorite workstation with Ubuntu 16.04
    DST = local working dirctory from which you start your job

# Work around #2
  I recompiled the matlab binaries explicitly to use the mpfr libraries
  shipped with matlab (libmpfr.so.1); so if you update your binaries,
  these should now also work on the th-cl-naples-* machines.

# -------------------------------------------------------------------- #
# Note concerning a new matlab/2018[ab] -> intel/mkl library bug
# when parallelizing based on QSP_NUM_THREADS>1 and MKL_NUM_THREADS>1 *)
# that has not been present in matlab/2016 // MKL_BUG_OMP
# Wb,Jul12,19

  Unfortunately there is a bug in the Intel MKL / intel OMP
  libraries shipped with matlab>=2018 that can result in a deadlock
  for nested parallelism. This was also confirmed by mathworks
  support based on a minimal sample code from my side. This bug
  occurs completely independently of, and hence has nothing
  to do with QSpace per se.
     
  Within QSpace, QSP_NUM_THREADS represents an outer level of
  parallism e.g. used when working through a possibly long list
  of tensor operations such as contractions. In contrast,
  MKL_NUM_THREADS sets the parallelization at the level of
  the MKL library, e.g. for an individual matrix multiplcation.

  Now having QSP_NUM_THREADS>1 *AND* MKL_NUM_THREADS>1 *)
  with OMP_NESTED=true, this nested parallelism can result
  in a deadlock such that an entire job just stops:
  meaning it is still alive, but falls asleep at 0 CPU load.
  What one sees in the debugger is that one of the threads
  in the outer level of parallelization, even though finished,
  gets stuck for unknown reasons in a wait state that
  never terminates.

  Hence, contrary to the earlier matlab/2016 setting,
  nested parallelism needs to be turned off in matlab>=2018.
  That means in practice, that if you encounter jobs that stall
  in the above sense, you need to make sure that only one
  of the paraellization levels is active: eiher
  QSP_NUM_THREADS>1 xor MKL_NUM_THREADS>1 *), but not both.

# preliminary conclusions drawn out of further discussion
# with Matan Lotem // Wb,Jul24,19 // UPDATE_NESTED

  Within QSpace, if lists are long, e.g., in the presence
  of multiple / larger symmetries or higher rank tensors,
  it should not really slow down a calculation if with
  QSP_NUM_THREADS*MKL_NUM_THREADS = nthreads = const,
  the balance is fully tilted towards
  QSP_NUM_THREADS = nthreads, whereas MKL_NUM_THREADS = 1.

  Now while not all blas/lapack routines parallelize well,
  i.e., can make full use of a large pool of parallel threads,
  working in parallel fashion through longer lists in principle
  should parallelize perfectly. In this sense, the following
  is the recommended environmental setting up to further notice:

    * set QSP_NUM_THREADS to the maximal number of threads
      you intend to use.

    * do not set OMP_NESTED=true (default is false)
      by which MKL automatically goes serial when called
      within a parallel section; meaning the value of
      MKL_NUM_THREADS if set, is ignored.

    * do not set any MKL_* environmental variables either
      [except in the sense of setNumThreads() below *)]

  [only for the rather uncommon situation that one uses no symmetries
  at all within QSpace, i.e., where every tensor consists of one
  full/large block, it appears desirable to set MKL_NUM_THREADS maximal,
  while taking QSP_NUM_THREADS->1].

  *) MKL_NUM_THREADS is actually unset e.g. by the matlab desktop;
     so in order to set the number of threads in matlab (and hence MKL)
     you may rather use my routine setNumThreads() which uses
     feature('numthreads',nthreads).

# -------------------------------------------------------------------- #
# increase soft-limit on max. number of open files (default: 1024)
# Wb,Jun08,19

  the soft limit on open files (ulimit -n) set by default to 1024
  can be too low for matlab!

  Background: when matlab starts, it typically opens about 50(!)
  threads right away without having done any calculation.
  e.g. see /proc/<matlab-process-id>/

  Moreover, mex files are effectively dealt with / compiled as
  shared libaries; now as it turns out, when a mex file is called,
  each thread keeps an open link to the mex file; so when calling
  e.g. some 20+ different mex routines, this quickly amasses 
  >1000 open files, without even having explicitly opened a
  single file! e.g. see lsof +D /path

  As a result, this can lead to errors when explicitly opening
  files within a program, since the soft limit is reached
  (resulting in seemingly unjustified errors that matlab
  cannot open clearly existing files)

  The soft limit on open files can be increased by each user,
  however, all the way up to the hard limit. Therefor you (also)
  may want to write something like "ulimit -S -n 4096" into your
  ~/.bashrc file. Ralph Simmler [LMU, 07/22/2019] recommends
  to use prlimit instead, i.e. "prlimit -p $$ --nofile=4096:"
  where the trailing column is important to leave the hard-limit
  unchanged.

# -------------------------------------------------------------------- #
# Wb,Jun01,19

  Minor Changes

  permuteQS and contractQS

     The optional permutation previously must have precisely matched
     its length to the rank of the affected QSpace. This has been relaxed:
     now the permutation can be shorter than the QSpace it is
     applied to: it then simply acts on the initial group of indices,
     leaving the order of the trailing remainder of indices intact.

  C=getIdentityQS(A,..,B..)

     added "#i:j.." options to inherit itags from A or B in C
     see help getIdentityQS for more details

  orthoQS

     option 'itag' on newly added intermediate index was broken; fixed.

# ==================================================================== #
# QSpace v3.1 (matlab>=2018)                                          #
# ==================================================================== #

  QSpace v3.1 moved on form matlab/2016b to matlab/2018[ab].
  This was partly motivated by new matlab installs outside LMU
  as well as the os-x10-mojave release which was already
  also effecively v3.1 and which requres matlab/2018b.

  With this upgrade, the binary QSpace v3.1 routines have
  undergone some significant changes (Wb,May16,19):

# Matlab>=2018 supports interleaved complex data format

  Matlab finally (after decades) also switched to the standard
  interleaved format for complex numbers (previously real and
  imaginary part of matrices were stored separtely).

  The interleaved format was already inside the QSpace mex
  files all along, yet the data had to be converted forth and
  back when interacting with matlab.  Essentially, this resulted
  in extra copies for complex data for every input and output.
  All of these have been removed in QSpace v3.1. If there is no
  explicit need to copy the input data, it is now referenced,
  not just for real, but also for complex data.

  That is, all mex-binaries in v3.1 are compiled with the
  mex flag -R2018a.

# Avoiding copies where possible

  Along with the above, I also changed some of the internal
  data handling for the output QSpaces, in that they can be
  allocated and constructed already in matlab format, such that
  no further copy is required when returning the data back to
  matlab. In particular, this concerns plusQS and contractQS
  in the absence of no outer multiplicity [1].
  In these case, the output is generated in-place, and no
  other copies are made.

# Improved parallelization

  plusQS was mostly operating in serial mode previously, also
  making redunant copies; while plusQS is typically not a
  performance bottle neck, it can become so, if contractions
  are relatively cheap in comparison, e.g. when moving towards
  tree tensor networks. The parallelization of plusQS listens
  to NMAX := max(QSP_NUM_THREADS,OMP_NUM_THREADS).

  Copying of QSpaces, where necessary, as well as data I/O
  for higher-precision MPFR data I/O for has been parallelized
  in the same spirit, also listing to NMAX (the latter is
  relevant, in practice, only for generalized Clebsch-Gordan
  tensors for larger non-abelian symmetries).

  For comparison, the parallelization of contraction can
  utilize up to NPROD := QSP_NUM_THREADS*OMP_NUM_THREADS (!)
  thread, namely parallelizing loops when working through
  the contractions in the QSpace.data sector based on
  QSP_NUM_THREADS, where each individual contraction itself
  further parallelizes by calling threaded blas/lapack
  library routines which listen to OMP_NUM_THREADS.

# Improved thread safety

  The mex API is not thread safe, as mathworks keeps emphasizing.
  And as it turns out, matlab>=2018 is even more touchy in that
  respect, e.g. as compared to matlab/2016b. Mex API
  routines called within a parallel section can crash entire
  matlab sessions which gives rise to rather cryptic `assertion'
  errors. Importantly, for matlab matlab>=2018 this even includes
  a plain simple printf() call, since printf is redefined
  by the mex API into mexPrintf(). The latter is necessary to
  display output in the matlab desktop (if one uses it), and
  not in the shell/terminal from which matlab was started.

  Hence log output from auxilliary threads in a parallel
  section is now buffered, and only flushed with the next
  log-output when back to serial mode (or at the end
  of function call, which ever comes first).


  [1] Outer multiplicity can only be present with non-abelian
      symmetries; for SU(2) this is the case for rank>3 tensors,
      and for general non-abelian symmetries for rank>2. Rank-2
      is always multiplicity free, and so are fully-abelian
      calculations for any rank.

# ==================================================================== #
# MPSPACK / QSpace v3.0 (matlab/2016b)                                 #
# this version is officially retired with the old matlab/2016b         #
# effective 04/20/2019 // AW                                           #
# ==================================================================== #
  
  MPSPACK_v3.0pre (prelease) synced to $WBX : Feb 2015
  MPSPACK_v3 synced to $WBX : Feb 2017

  added DMRG package based on @Hamilton1D [(C) AW 2014-2018]
    * synced to $WBX : Wb,Feb06,18
    * see ./DMRG/runDMRG.m

  For other more recent changes,
  ==> see end of this readme file!

# -------------------------------------------------------------------- #
# new environmental variables with MPSPACK v3.0

  RC_STORE

    path to somewhere in /data where you would like to
    store all Clebsch Gordan related data
    (if your path ends with .../RCStore/ this last directory
    will be created automatically if it does not exist).

    Update 09/2016: by now, this path may contain multiple paths,
    e.g: export RC_STORE=PATH1:PATH2:...:PATHN with the following
    procedural convention:

       Symmetry related data is looked up in the *last* PATHN
       first (because it is considered the most recent local path!);
       if not found, the search proceeds to PATH1, PATH2, etc.
       if not found all the way to PATH_(N-1), the requested data
       is computed and store in PATHN. Therefore paths PATH1 to
       PATH(N-1) are effectively read-only.

    Background: RC_STORE gets large for larger SU(N) symmetries;
    therefore when working on a cluster, one does not want to
    copy the entire RC_STORE for every single job. Therefore
    when using large symmetries such as SU(N>3) [for SU(N<=3)
    one may as well just synchronize or regenerate symmetry
    related data], the recommended setting is:

       RC_STORE=NETWORK_PATH:LOCAL_PATH

    where NETWORK_PATH e.g. is some globally visible
    /project/theorie/.. network directory, whereas LOCAL_PATH
    e.g. is some job specific directory on /data/$USER/myjob/..
    This realizes differential storage, in that only new or
    updated data is stored locally. In cases where the RC_STORE
    on the NETWORK_PATH was already complete, the LOCAL_PATH
    will then remain empty.  Conversely, if significant new
    data was generated locally, one can update i.e. rsync the
    RCStore in the NETWORK_PATH from some LOCAL_PATH.
    It is recommended to do this after all jobs are finished
    to avoid interference with other running jobs.

    Note that instead of using /project/*, the NETWORK_PATH
    can also point to some LOCAL machine's /data/$USER/ ...
    This directory can then be mounted in the job directory using

        cd /data/$USER/myjob/...
        mkdir ./RCStore
        sshfs $USER@some_machine:/data/$USER/etc. ./RCStore
        run_my_program

    which should be unmounted once the job is finished

        fusermount -u ./RCStore

  CG_VERBOSE

    level of log output;
    general guideline on the value of CG_VERBOSE:

     = 0  basically log nothing concerning CG activity
    <= 4  log mainly information regarding newly generated RCStore data
     > 5  also log information on RCStore data that is (newly)
          read from RCStore.

# -------------------------------------------------------------------- #
# comments on itags (index tags)

  QSpace v3. introduced 'itags' i.e. string labels for indices
  with up to 8chars (since, internally, these are converted long
  int format), and which  are specified in QSpace.info.itags.
  For example, this offers automated  contraction ('auto-contraction')
  of pairs of tensors solely based on matching itags! See 'contractQS -h'
  for more detailed information.

  Itags must always also contain individual conjugate flags     
  (this represent the bare minimum that must be specified with  
  each QSpacein v3): the conjugate flags on individual indices  
  (legs) of a tensor determine  whether that index (leg) is     
  in- or out-going, with the convention that                    
																
	  all out-going indices have a trailing * in their itags!   
																
  For example, an A-tensor with L(eft), R(ight), s(=local)      
  indices may have itags A.info.itags={'L','R*',s'}             
  assuming (L,R,s) index order.                                 

# -------------------------------------------------------------------- #
# Notes regarding Clebsch Gordan tensors

# In general: do not even think about tinkering with the new
# info.cgr references to Clebsch Gordan spaces!

  e.g. must be careful about signs for degenerate q-labels
  which are left unchanged through permutation! e.g. for SU2:
  permute((1,1;0),[2 1 3]) acquires a minus sign!

# only exception for (small) CGC data during setup:
# set basic CData info fields (type,qset,qdir,cgw),
# while keeping the rest empty:

  c=emptystruct(A_someQSpace.info.cgr);
    c.type=qtype;
    c.qset=repmat(q,1,numel(Ak.Q));
    c.qdir=Ik.cgr(1).qdir;
    c.cgw=['1' 0];

  A_someother.info.cgr(:,j)=c;

# -------------------------------------------------------------------- #
# further comments regarding specific (new) routines

  @QSpace/getvac.m 

    allows to obtain vacuum state
    (through identity QSpace) for given symmetry setting,
    e.g. see setupSIAM_SU2x2.m: q=getvac(Z);

# -------------------------------------------------------------------- #
# removed environmental variable WB_CONTRACT_QFLAG
  and with it the contractQS( '-[lL]' options! // -l, -q
  => using automated lenient flag, instead // Wb,Nov22,14

# renamed routines for better naming consistency
  see $MEX/rename_qs.sh and archived versions therefrom

     for historical reasons
        MPSPACK_v1  => got QS label
        MPSPACK_v2  => got CG label

     mpsMakeUnique   => makeUniqueQS   # Wb,Nov22,14
     mpsNormQS       => normQS
     mpsPermuteQS    => permuteQS
     mpsPlusQS       => plusQS
     mpsSkipZerosQS  => skipZerosQS
     mpsGetDim       => getDimQS
     mpsGetQDim      => getQDimQS
     mpsIsIdentityCG => isIdentityCG
     mpsIsIdentityQS => isIdentityQS
     mpsTimesEl      => timesElQS
     mpsIsDiagQS     => isDiagQS

   * mpsMaxDiffQS    => maxDiffQS      # (mostly) unused
   * mpsGetDRange    => getDRangeQS    # (mostly) unused

     NRGWilsonCG     => NRGWilsonQS    # Wb,Nov24,14
     fdmNRG_CG       => fdmNRG_QS

     mpsOrthoQS      => orthoQS        # Wb,Aug25,15
     mpsGetSVD       => svdQS          # Wb,Aug25,15
     mpsEigenSymmQS  => eigQS          # Wb,Aug25,15

# -------------------------------------------------------------------- #
# Wb,May01,14
  outsourcing CGC's in QSpace into CGC library (single file)
  => not very flexible, thus replaced by the following:

# Wb,Dec15,14
  outsourcing CGC's in QSpace into CGC library (file system)
  => requires new environmental variable RC_STORE
  => since MANY files are generated, DO NOT store
     in backed up file systems such as HOME! => use /data

# -------------------------------------------------------------------- #
# Wb,Aug25,15 :: update / renamed routine

  mpsEigenSymmQS  => eigQS

  the behavior of this routine changed, in that in the presence
  of non-abelian symmetries the first return argument now contains
  a 2nd column which contains degeneracy (i.e. combined multiplet
  dimension)

# -------------------------------------------------------------------- #
# Wb,Aug24,15 :: reactivated / renamed routine

  mpsOrthoQS => orthoQS
  mpsGetSVD  => svdQS

  Here svdQS() is a specialization of orthoQS
  in that instead of returning [U,S,V] it returns [U,(SV')'=VS];
  Both routines accept truncation arguments (Nkeep, stol)

  Note that orthonormalization only makes sense in the presence
  of orthonormal input spaces! Hence the above routines now
  require in most cases that all indices are of the same type
  (i.e. all-in).

  The input index (idx) indicates which indices shall be 
  combined. Since non-abelian symmetries require proper multiplet
  fusion, here the routines orthoQS and svdQS have been 
  constrained such that idx eiter represents a single index,
  or (r-1) indices with r the rank (=dimensionality) of the 
  QSpace tensor.

# -------------------------------------------------------------------- #
# switched from matlab/2013a to matlab/2013b
# since it supports gcc-4.7.4 (-std=c++11)
# which now supports 'explicit' type conversion
# Wb,Jan04,16

  WRN! starting with matlab/2013b, flip() is also a matlab routine now!
  it has pretty much the same in syntax as my previous lib/flip.m which
  led to infinite recursive calls [e.g. my flip() calls flipud -> calls
  flip(x,1), which however again reroutes to my flip.m routine]
  Solution: simply removed my lib/flip.m routine.
  Wb,Jan04,16

  Drawback: 'clear all' appears to unload shared libraries
  => error message: shared object not open
  => requires to restart Matlab to return to regular operation

# -------------------------------------------------------------------- #
# altered treatment of decomposition of tensor-product of ireps
# Wb,Jan04,16

  further parallelized sparse tensor-product decomposition
  e.g. to speedup SU(6) calculations; this is mostly irrelevant
  in terms of speed for SU(N<=4) calculations.

  orthogonalization of inner multiplicity spaces is now 
  based on QR decomposition (the two-fold Gram-Schmidt used so far
  for numerical stability is numerically more expensive, yet still
  not as accurate anyway)

  => NB! this can lead to rotation in inner multiplicity spaces
     for SU(N>=3); therefore old RC_STORES need to be recomputed!

# -------------------------------------------------------------------- #
# minor adjustments for qlabel representation of non-abelian symmetries:
# Wb,Jan04,16

  For a more readable / compact representation of multiplet labels
  (q-labels) of non-abelian irreps, the numbers 10-35 are now replaced
  by the letters A-Z. E.g. the SU(4) multiplet (8 10 12) is 
  written in the compact form (8AC); similarly, this also affects
  SU(2) in that for example the multiplet q=2S=12 may be simply
  written as q=(C).

  This only affects log-outputs and file names in the RC_STORE.

# -------------------------------------------------------------------- #
# change in log behavior # Wb,Jan09,16

  log output referring to computation of generalized CGCs has
  been removed from STDOUT and are now written centralized
  with the RC_STORE; for larger CG_VERBOSE, the output to
  the log-file eventually is also replicated to STDOUT.

# -------------------------------------------------------------------- #
# I changed my default Matlab version from 2013b to to R2016a
# Major reason: severe SVD bug in Matlab 2013b .. 2015b
# Wb,Aug30,16

  Since in the meantime matlab also changed its graphics handles
  to actual class objects, many minor changes in scripts that
  produce graphics were incurred. Several plot scripts that
  I have used in the meantime, are updated. More to follow.

  I kept the old version in MPSPACK_v3.0pre/
  and moved the new version to MPSPACK_v3.0/

  Mostly affected directories

    bin/           please update all mex files
    util/          please update all mex files, etc.
    lib/           please change as far as necessary.
    Class/@QSpace
    Class/@SymOp

# -------------------------------------------------------------------- #
# Added Z_n symmetry
# this essentially acts like U(1) with charge q taken module n.
# Wb,Aug31,16

  New abelian symmetry labels that can be used to replace 'Acharge'
  (already also with getLocalSpace)

   * 'ZNcharge' with N=2,3,...
      total U(1) charge symmery reduced to a Z_N symmetry
      with q-labels = [U(1) charge modulo N] = {0,1,...,N-1}

   * 'Pcharge'  charge parity
      this is equivalent to 'Z2charge', except that symmetry
      labels are {-1,1} instead of {0,1}).

# -------------------------------------------------------------------- #
# Parallelization in QSpace // Wb,Feb07,17

  All of the parallelization in QSpace concerns threading
  in shared memory environments and is therefore based
  on the OpenMP standard. It is controlled via a set of
  environmental variables.

  Relevant environmental variables:

  OMP_NUM_THREADS
   * generic variable that specifies the default maximum number of
     threads when parallelizing
   * affects shared public libraries with internal parallelization
     such as lapack/blas
   * also used with QSpace via my wb::sparse() class to deal
     with generic sparse tensors (e.g. all Clebsch-Gordan data)

  MKL_NUM_THREADS
   * affecs the shared intel MKL libraries (similar to lapack/blas)
     note that these are delivered with the matlab API, and hence
     are the ones that are included in mex-files.

  QSP_NUM_THREADS

   * introduces a higher level parallelization within the QSpace
     tensor library. For example, consider the following piece
     of matlab-like pseudo-code

     | ee=cell(1,numel(A.data));
     |
     | parallel_for i=1:numel(A.data) <-- distribute this to QSP_NUM_THREADS threads
     |    ee{i}=parallel_eig(A.data{i}); <-- calculate using MKL_NUM_THREADS threads
     | end

     > each of the QSP_NUM_THREADS threads that gets assigned a
       unique value i in the loop therefore has a task [here eig()]
       which itself parallelizes, say using MKL_NUM_THREADS threads.
       This is referred to as `nested parallelization'.
  
     > hence the actual maximal number of threads is given by
       QSP_NUM_THREADS * max(OMP_NUM_THREADS, MKL_NUM_THREADS) !!

     > for cluster jobs: make sure that the number or requested
       cores matches, i.e. never exceeds (NTASK * QSP_NUM_THREADS).

   * the QSpace parallelization based on QSP_NUM_THREADS affects both,
     A.data{} as well as A.info.cgr(), e.g. parallelizes when working
     through reduced sets of reduced matrix elements as well as their
     corresponding generalized Clebsch Gordan coefficient data.

   * so far, QSP_NUM_THREADS mostly affects the contractQS() routine
     where by matching the symmtery sectors of two higher-rank tensors,
     say A and B, can quickly result in a huge list to work through.

   * a change in QSP_NUM_THREADS will be seen by mex routines
     if they are called after the changed (i.e. a 'clear functions'.
     may be required within matlab).

   * threads get locked (i.e. need to wait) when other threads
     e.g. already happen to work on the same Clebsch-Gordan tensor;
     this is based on standard OpenMP locking. Therefore large
     QSP_NUM_THREADS will results in overhead with less overall
     speedup.

   * I have tested QSP_NUM_THREADS with (NRG + fdmNRG) both, with
     abelian and non-abelian symmetries. If you see strange kinds of
     crashes try turning off QSP_NUM_THREADS (the reason being that
     different threads may get into each others way, despite that
     I tried to eliminate as much of possible of this scenario).

  The default value for the above variables is 1, i.e. do not parallelize
  (if you experience something different for mex-compiled routines,
  please let me know).

  Suggested usage:

   * set OMP_NUM_THREADS and MKL_NUM_THREADS the same, e.g.
     with NTASK the number of threads assigned to a given task,

        export OMP_NUM_THREADS $NTASK
        export MKL_NUM_THREADS $NTASK

        and within matlab: >> num_threads(NTASK);  % see MPS_PACK/lib

   * if one expects large A.data{} sets [which carry (reduced)
     matrix elements as in typical NRG calculations], you may simply
     ignore QSP_NUM_THREADS (default value for QSP_NUM_THREADS is 1).

   * if you are dealing with larger-rank tensors such as in PEPS
     that can have QSpace with lists of lengths hundreds to thousands,
     you man want to turn on QSP_NUM_THREADS while reducing NTASK
     accordingly.

  When using interactive matlab:
  it appears matlab completely ignores MKL_NUM_THREADS etc;
  Mathworks support 08/02/2016

  "MATLAB does not respond to MKL_NUM_THREADS or OMP_NUM_THREADS,
   as these settings are overridden. While MKL_DYNAMIC and
   OMP_DYNAMIC are not overridden, the default MKL_DYNAMIC and
   OMP_DYNAMIC settings are true. When they are set to true, MKL
   will use the up to the number of threads specified with
   maxNumCompThreads(). [...] But unfortunately, as of MATLAB/R2009b,
   the use of maxNumCompThreads() is now deprecated. So the only
   possible way of doing this is by setting the CPU affinity for
   MATLAB process at the OS level. Given that you are working
   on MATLAB R2016a, the last option i.e. 'Setting CPU affinity
   for a MATLAB process at the operating system level' is the
   only option available to you."

   e.g. use 'taskset 0xff ml' to start matlab constrained
   to 8 cores.

   On the cluster, the number of cores may be enforced by the
   cluster scheduler to ensure that a job does go in excess
   of the specified number of cores.

# -------------------------------------------------------------------- #
# fixed memory leak in QSpace (data sector)

  thanks to Benedikt for pointing this out (e.g. within svdQS)
  Wb,Jun27,17

# -------------------------------------------------------------------- #
# altered behavior of Nkeep/stol in svdQS and orthoQS
# in response to comments by Wei Li // Wb,Sep17,17

  Nkeep is interpreted as Nkeep_max
   - if Nkeep is specified, but stol is not,
     the default value of stol (1E-8) changes to 0
   - added option 'Nkmin' to specify minimum number
     of states/multiplets to keep

# -------------------------------------------------------------------- #
# added vector behavior to 'NC' in getLocalSpace
# in response to discussion with Fabian Kugler and Seung-Sub Lee
# Wb,Feb09,18

# example
  [FF,Z,SS,IS]=getLocalSpace('FermionS',...
  'SU2spin,Acharge,SUNchannel','NC',[2 1],'-v');

  {IS.SOP.info}' # for more info on symmetry setup

  [FF,Z,IS]=getLocalSpace('Fermion','Acharge,SUNchannel','NC',[2 1],'-v');

# -------------------------------------------------------------------- #
# Wb,Nov24,18

  Minor, but likely, convenient update with contractQS:
  (also updated help section with 'contractQS -h'; please have a look!)

  the option "-op:.." now accepts (posix) regular expressions.
  These provide a hint on where to contract without giving the
  precise label (conj-flags are ignored).
  
  For example, if A.info.itags={'K10','K11*','s11'} is some
  A-tensor with already well-defined itags, whereas 
  B.info.itags={'','*','*'} is some local operator that may act
  on any site, hence with only the (required) conj flags set, then

     contractQS(A,nloc,'-op:s\d+')

  searches for the a matching itag in A, and in given case
  finds the regex match 's11';  similarly, contractQS(A,nloc,'-op:^s')
  looks for an itag that `starts with s'; if the result is unique,
  this will, for the purpose of the contraction, assume
  the itags {'s11','s11*'[,'op*']} for B (the third tag,
  only if rank-3); this ensures that the itag of the local
  index is also carried on after the contraction.

  Special characters such as the back-slash or ranges [..]
  are taken as indicator to interpret given string as a
  regexp; plain text strings after '-op:' still must find
  and exactly matching itag.

  Technical note: posix regexp actually does not understand \d,
  which in posix need to be written as  [0-9] or [[:digit:]];
  in simple examples, as in the above, QSpace translates
  \d to the posix equivalents.

# -------------------------------------------------------------------- #
# Wb,Dec12,19

  As a follow-up to discussions with Mathan Lotem:
  matlab changed its default eig() solver from dsyev() to dsyevd()
  which is about a factor x2 faster for all matrix sizes!
  [the switch to dysevd() was indicated by mathworks support].
  Hence the following changes were made to QSpace:

     dsyev  -> dsyevd()
     zheev  -> zheevd()
     dgesvd -> dgesvdd()
     zgesvd -> zgesvdd()

# follow-up by bug-report by Seung-Sup // Wb,Jan12,19

  dsyevd() can return `internal error' info>0 (!!)

  This may be related to the current lapack version used by matlab;
  since such an error was actually reported and apparently fixed 
  in LAPACK; as a temporary solution: I changed UPLO='U'->'L'
  which in given case solved the problem.

# -------------------------------------------------------------------- #
# Wb,Dec16,18

  Changes to my mex logging routine

  With threading in mind, where individual threads can write
  simultaneously to stadout, my mex-log routine now completes
  each entry line as a string first, and then flushes the
  string, i.e. the entire log entry. This is intended
  to avoid mingled log output from  different threads.

  Moreover, in the past it occured at times that the matlab
  desktop wrote its mex text output to the calling terminal (!)
  rather than the matlab desktop. By switching (back to)
  matlab's mexPrintf(), this should have changed by now,
  in that mex output should occour within the matlab desktop.

  One persisting nuisance: matlab does not understand \r
  (carriage return), but interprets it the same way as \n (newline).
  Morever, the backspace \b erases the current line, hence
  is also not as useful. Thereforee, while in command-line mode
  in the terminal, mex output can overwrite previous lines
  simply by printing \r (e.g. for showing progress report),
  in the matlab desktop mode this leads to newlines,
  hence to more extended `running' output.

# -------------------------------------------------------------------- #
# Wb,Jan12,19

  Updated the way help works for many mex-files by outsourcing
  the help header for QSpace mex-files into a separate m-file.
  Note that the sole purpose of a m-file myprog.m is that it
  provides the help for myprog.mexa64; hence the m-file only
  contains a (long) matlab comment. 

  >> myprog --help
  >> myprog -h
  >> myprog -?

  are then equivalent to the standard matlab usage >> help myprog
  For the QSpace mex-files, the helper m-files are newly generated
  and all copied into the bin folder (same as for the 
  corresponding mex-files; note that in order for this work,
  the m-files need to be in the same folder as the mex-file of
  later in the matlab path.

  For some of the (non-QSpace) utility routines, the help header
  is still left with the mex-file. Therefore the old usages above
  e.g. >> myprog -h
  still work perfectly well for all mex files. If the help m-file is
  missing [because not with the mex-file] a proper warning is issued.

# -------------------------------------------------------------------- #

