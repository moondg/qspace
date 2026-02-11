
# Makefile for mex files to be compiled on the UNIX prompt, e.g.
# > make -f Makefile getIdentityQS
# > make -B -f Makefile mpsset

  SHELL = /bin/bash

# allow for parallel compilation
# WRN! may be ignored when set in here; then set with calling make
# or define MAKEFLAGS in the environment
  ifdef MSLOTS
     MAKEFLAGS += -j $(MSLOTS)
  endif

#! make sure that system/matlab_setup.pl was run as part of
#! the environmental setup -> enforce definitions of ENV variables:
#! MATLAB_ROOT, ARCH, MEX_R2018
  q := $(shell ../system/matlab_setup.pl --check; echo $$?)
  ifneq ($(q),0)
     $(error Please check QSpace environment)
  endif

  ifndef MXC_DIR
     MXC_DIR = ../bin
  endif

  ifndef CLIBS
   # NB! -liomp5 moved to mexopts.sh file // Wb,Jan14,16
     CLIBS = -lmwblas -lmwlapack # -lmysqlclient # -liomp5
  endif

  OPTS = # actual (direct) mex compiler options
  XFLAGS = -Wall # goes into CXXFLAGS
  LFLAGS = # goes into LDFLAGS
  DEFS =

  I = -I.
  L = -L.

# the following options assume gcc on linux and Xcode/clang for macOS;
# this is checked with `make test' in makelog.pl
# NB! removed ML_GCC_VERSION since mex -setup may use C++ compiler
# other than what is on the path via `which gcc' // Wb,Feb19,24

  ifeq (maca,$(findstring maca, $(ARCH))) # mac on ARM // Wb,Jun11,25
     GMP = /usr/local/lib/libgmp.dylib
     MPF = $(MATLAB_ROOT)/bin/$(ARCH)/libmpfr.dylib
     OMP = $(MATLAB_ROOT)/bin/$(ARCH)/libomp.dylib
     ismac = 1

  else ifeq (maci,$(findstring maci, $(ARCH))) # mac on Intel
     GMP = /usr/local/lib/libgmp.dylib  # Wb,Oct19,20 / matlab/2020b
     MPF = $(MATLAB_ROOT)/bin/maci64/libmpfr.1.dylib # Wb,Sep19,19
     OMP = $(MATLAB_ROOT)/sys/os/$(ARCH)/libiomp5.dylib
     ismac = 2

  else # assume linux environment
     GMP = -lgmp
     MPF = $(MATLAB_ROOT)/bin/glnxa64/libmpfr.so.1
     OMP = $(MATLAB_ROOT)/sys/os/$(ARCH)/libiomp5.so
     ismac = 0

   # WRN! no comment with leading space for macos
   # otherwise the space is part of the value!
     CLIBS += -lrt
     # required by clock_gettime() starting gcc6 [matlab/2018a]
     # (this option is gcc/linux specific, i.e., is not present for
     # clang on osx/mac at all; including it would result in error)
  endif

  ifndef MPF_LIB
   # MPF_LIB = -lmpfr -lgmp
     MPF_LIB = $(MPF) $(GMP)
  endif

  ifdef QS_USING_OMP # Wb,Nov29,22
     ifeq ($(shell expr $(QS_USING_OMP) \<= 0), 1)
         DEFS += -DQS_USING_OMP=0 # unset
     endif
  endif

# QSpace now permits to link against the high-performance
# tensor transpose (HPTT) library -lhptt
# This library needs to be installed separatly.
# In order to use it, simply define the environmental
# variable HPTT_ROOT to point to the root directory 
# of its installation prior to calling this Makefile
# Wb,Dec06,24
  ifdef HPTT_ROOT
      DEFS += -DQS_USING_HPTT
      I += -I$(HPTT_ROOT)/include
      L += -L$(HPTT_ROOT)/lib
      CLIBS += -lhptt
      LFLAGS += -Wl,-rpath,$(HPTT_ROOT)/lib # required on linux
  endif

# using GFLAG=1 with make for debugging purposes
# e.g. more detaild stack backtrace in case of error // Wb,Jun26,20
  ifdef GFLAG
     OPTS += -g  # keep with OPTS!
     XFLAGS += -O1
  else
     XFLAGS += -O3
  endif

  ifeq ($(ismac),0)
     XFLAGS += -Wno-format-truncation
     CFLAGS =
  else
     CFLAGS = -Xclang
  endif

  CFLAGS += -fopenmp -pthread

  XOPT = CFLAGS="\$$CFLAGS $(CFLAGS)"

# avoid compiler warning on (invalid) target macOS version 10.14
  XOPT += MACOSX_DEPLOYMENT_TARGET=12.0

  ifdef XFLAGS
     XOPT += CXXFLAGS="\$$CXXFLAGS $(XFLAGS)"  # works :)
  endif
  ifdef LFLAGS
     XOPT += LDFLAGS="\$$LDFLAGS $(LFLAGS)"    # Wb,Apr24,20
  endif

  XOPT = CXXLIBS="\$$CXXLIBS $(OMP)"

  ifdef VFLAG
     OPTS += '-v' # e.g. include by running: make ... VFLAG=1
  endif

  ifdef MEX_R2018
     OPTS += $(MEX_R2018)
  else
   # NB! marlab2018s -R2017b mex option
   # is not / no longer compatible with -largeArrayDims // Wb,Mar16,19
   # "Calling MEX with both 'R2017b' and 'largeArrayDims' not supported"
     OPTS += '-largeArrayDims'
  endif

# print general options used once and for all // tags: using
# OPTS=$(OPTS) // shown with each `make.pl --make ..' later anyway
  $(shell make.pl --make-info $(MATLAB_ROOT) $(MCC_TAG) \
    OPTS=$(OPTS) XFLAGS=$(XFLAGS)  LFLAGS=$(LFLAGS) );

  DEFS += "__DATE_COMP__=$(shell date '+%m/%d/%Y %Z%H%M%S')" # __TIME__ TSTAMP
  DEFS += HOST_NAME=$(shell hostname -s) # Wb,Nov05,21
  DEFS += MEX_EXT=$(shell mexext) # Wb,Jun27,23
  ifdef MATLAB_VERSION
     DEFS += MEX_MLVER=$(MATLAB_VERSION)
  endif
  ifdef DBSTOP
     DEFS += -DDBSTOP
  endif

# NB! make.pl checks whether DBSTOP is defined
  DEFS := $(shell make.pl -v --make-defs $(DEFS) "MEX_CXX_FLAGS=$(XFLAGS)")

  MEXC = $(MATLAB_ROOT)/bin/mex # see mex -h[elp] on unix prompt

# NB! default C++ compiler is chosen via 'mex -setup cpp'
# this is the recommended way as this automatically looks for compilers
# that are compatible with particular matlab version.
# mex -setup may be called from a system shell or within matlab
  ifdef CXX
     MEXC += "CXX=$(CXX)"
  else
     ifdef GCC
        MEXC += "GCC=$(GCC)"
     endif
  endif

  I += -I$(MATLAB_ROOT)/extern/include
# I += -I/usr/include/x86_64-linux-gnu -I$(MATLAB_ROOT)/extern/include
  L += -L$(MATLAB_ROOT)/bin/$(ARCH)

  ifeq ($(ismac), 1)
     I += -I/usr/local/include # gmp.h, omp.h
  endif

.SUFFIXES: .c .cc

  STD_DEP =

# --------------------------------------------------------------------- #
# MAKE SURE MATLAB_ROOT PATHS ARE SET PROPERLY (RUN ~/.MATLAB_SETUP -f -ld)
# --------------------------------------------------------------------- #
# implicit rules (basic recipe to compile), default dependency

# ./make.pl --odir $< "$$MATLAB_ROOT;$$MCC_TAG;$$CLIBS;$(OPTS)"` \

.cc: $(STD_DEP)
	@( eval "$$(../system/matlab_setup.pl -LD)"; XLIBS=; \
       if source_check.pl --cgc $<; then XLIBS="$(MPF_LIB)"; fi ; \
       OD="`./make.pl --odir $< ";;$(CLIBS) $$XLIBS;"`"; \
       $(MEXC) $< $(OPTS) $(L) $(I) $(CLIBS) $$XLIBS $(XOPT) \
      -outdir "$$OD" $(DEFS)  2>&1 | makelog.pl \
    )
#	NB! -largeArrayDims required by mxGet(Ir|Jc)

# -------------------------------------------------------------------------- #
# -------------------------------------------------------------------------- #

# default: tstmex
default: mpsset nrg FDM TDM fgr mpsother rcdata 

all: default utils

FDM: fdmNRG_QS getSmoothSpec
fdm: fdmNRG_QS
fgr: fgrNRG

TDM: tdmNRG    getSmoothTDM
tdm: tdmNRG

nrg: NRGWilsonQS
fdmnrg: fdm nrg

svd: svdQS
eig: eigQS
ortho: svdQS orthoQS

mpsset: getSymStates compactQS contractQS plusQS 

rcdata: setupRCStore getRC

mpsother: \
    skipZerosQS isIdentityCG isIdentityQS getIdentityQS \
    traceQS diagQS normQS getDimQS getQDimQS permuteQS \
    eigQS svdQS orthoQS makeUniqueQS isHConjQS

# traceQS: Wb,Aug17,16

mpsall: mpsset mpsother rcdata

# NB! util als tries to compile util.cc
# dealstruct cell2matc mat2cellc // Wb,May07,19
utils: helloworld mexworld printfc getchar filestat \
       kramers matchIndex mpfr2dec mpfrplus sparse2wb \
       uniquerows wbeig wbhash wbhist wbrat wbsvd # wbtrace getpid 

mpfr: mpfr2dec mpfrplus

world: helloworld mexworld

# -------------------------------------------------------------------------- #
# Matlab requires very specific C++ compilers for mex/mcc compilation;
# for supported matlab/C++ compiler combinations, see
# https://www.mathworks.com/support/requirements/supported-compilers-linux.html
# https://www.mathworks.com/support/requirements/previous-releases.html
#
#    R2016a       -> gcc-4.7.x                        # macOS
#    R2016b       -> gcc-4.8.x            R2018-20   -> Xcode [10-11].x
#    R2018-R2020a -> gcc-6.3.x            R2021      -> Xcode [11-12].x
#    R2020b       -> gcc-8.[3-5].x        R2022a     -> Xcode [12-13].x
#    R2022b       -> gcc-[7-10].x         R2022b-23  -> Xcode [12-14].x
#
# -------------------------------------------------------------------------- #

# Make sure to have successfully run 'mex -setup C++':
# based on a test compile here (without output due to -n flag)
# check whether mex runs correctly, including C++ compiler imployed.
# Formerly, this was part of system/matlab_setup.sh with explicit
# additional compiler check in system/matlab_setup.pl // Wb,Feb19,24

env: test

test:
	@(eval "$$(../system/matlab_setup.pl -LD)" ; \
      mex -n -v helloworld.cc 2>&1 | makelog.pl --test \
    )

# -------------------------------------------------------------------------- #
# Wb,Jan23,18

# NB! compiles, but by now also returns an empty array!?

feature: feature_tst # NB! keep feature('s') as matlab's command!

feature_tst: feature_tst.cc
	@$(MEXC) $< $(L) -lmwservices -outdir $(MXC_DIR)

# keep this (otherwise this tries to make util.cc into a binary)
util:
	@echo call "'make utils'" instead!

clean:
	@rm -vf ../bin/*.o ../util/*.o

# -------------------------------------------------------------------------- #

