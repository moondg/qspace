/* ---------------------------------------------------------------------
 * Project : QSpace tensor library (v4.0)
 *
 * Copyright 2024 Andreas Weichselbaum
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * --------------------------------------------------------------------- */

#ifndef __WB_LIB_MEX_HCC__
#define __WB_LIB_MEX_HCC__

/* ------------------------------------------------------------------ //
   header file for mex programs
   A. Weichselbaum (C) Jan 2006
// ------------------------------------------------------------------ */

#define PP_STR__(a) #a
#define PP_STRFY(a) PP_STR__(a)

#ifndef QS_USING_OMP
   #define QS_USING_OMP  1
#elif QS_USING_OMP<=0
   #undef QS_USING_OMP
#endif

#if __linux__
#elif __unix__
   #warning got unix
#elif __APPLE__
#else
   #error GOT UNSUPPORTED OPERATING SYSTEM
#endif

#include <assert.h> 
#include <execinfo.h> 

#include <time.h>
#include <signal.h>
#include <unistd.h>
#include <fcntl.h> 
#include <regex.h> 
#include <regex>   

#include <cstdio>
#include <cstdlib>
#include <chrono>    

#include <exception> 

#include <cmath>

#include <ctype.h>
#include <climits>
#include <stdarg.h> 
#include <string.h>

#include <string>
#include <iostream>
#include <typeinfo>

#include <sys/ipc.h>
#include <sys/msg.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <map>
#include <queue>         
#include <unordered_map> 
#include <algorithm>     
#include <vector>        

#ifdef QS_USING_OMP
#include <omp.h>   
#else
   int omp_get_max_threads() { return 1; }
   int omp_get_num_threads() { return 1; }
   int omp_get_thread_num()  { return 0; }
   int omp_in_parallel()     { return 0; }
   int omp_get_level()       { return 0; } 
#endif

#if __APPLE__
  #ifndef exp10
     inline double exp10(double x) { return pow(double(10),x); }
  #endif
  #ifndef sighandler_t
     #define sighandler_t sig_t
  #endif

  #include <mach/mach_host.h> 
  #include <mach/task.h>      

  #include <sys/sysctl.h>     
  #undef setbit 
  #undef isset  

#endif

   unsigned ARG_CHECK=1;

#define STRLEN  1023
   char str[STRLEN+1];

#define sprintf_str(...) snprintf(str,STRLEN,__VA_ARGS__)

#define MX_CHECK_HELPER_NARGS(n1,n2,m) \
   if (nargin>0 && \
       checkHelpVersion(argin[0], nargout ? argout : NULL)) { return; }; \
   if (nargin<n1 || (n2>=0 && nargin>n2)) { \
      strcpy(str,"ERR invalid number of input args"); unsigned l=strlen(str); \
      if (nargin<n1) \
           { snprintf(str+l,32," (%d / %d required)",nargin,n1); } \
      else { snprintf(str+l,32," (%d / %d %s)",nargin,n2, \
         n1==n2? "expected":"at most"); }; \
      if (nargin || nargout) wblog(FL,str); else usage(FL,str); \
   }; \
   if (m>=0 && nargout>m) { \
      sprintf_str("ERR invalid number of output args (%d)",nargout); \
      if (nargin || nargout) wblog(FL,str); else usage(FL,str); \
   }

template <class T> class wbvector;
template <class T> class wbMatrix;
template <class T> class wbarray;
template <class TD> class wbsparray;
template <class TQ, class TD> class QSpace;
template <class T> class sparseIndex2D;

class wbcomplex;
class wbstring;
class bitset;
class wbperm;
class wbindex;

class iTags;
class ctrIdx;

namespace CG {
   class FileLock;
}

namespace Wb {
   class Clock;
   class ClockSet;
}

using namespace std;

#ifdef MATLAB_MEX_FILE
class MXPut;
#endif

#define ENABLE_IF_COMPLEX(T__) typename \
   std::enable_if<  WbUtil<T__>::isComplex(), T__ >::type* = nullptr

#define ENABLE_IF_isPOD(T__) typename \
   std::enable_if<  WbUtil<T__>::isPOD(), T__ >::type* = nullptr

#define ENABLE_IF_noPOD(T__) typename \
   std::enable_if< !WbUtil<T__>::isPOD(), T__ >::type* = nullptr

#define ENABLE_IF_isINT(T__) typename \
   std::enable_if<std::is_integral<T__>::value>* = nullptr

#if 1
   #define   widx_t size_t  
   #define  swidx_t long    
   #define  wperm_t size_t
   #define swperm_t long   
#else
   #define   widx_t unsigned
   #define  swidx_t int
   #define  wperm_t unsigned
   #define swperm_t int
#endif

#if 1
   #define  SPIDX_T  size_t
   #define sSPIDX_T  long
#else
   #define  SPIDX_T  unsigned
   #define sSPIDX_T  int
#endif

#define WBINDEX  wbvector<widx_t>
#define WBPERM   wbvector<wperm_t>
#define WBIDXMAT wbMatrix<widx_t>

#define   UVEC        wbvector<unsigned>
#define C_UVEC  const wbvector<unsigned>
#define   DMAT        wbMatrix<double>
#define C_DMAT  const wbMatrix<double>
#define   UMAT        wbMatrix<unsigned>
#define C_UMAT  const wbMatrix<unsigned>
#define   IMAT        wbMatrix<int>
#define C_IMAT  const wbMatrix<int>
#define C_TMAT  const wbMatrix<T>
#define C_MEX   const mxArray
#define C_UINT  const unsigned
#define C_CHAR  const char

#define DMAT0   wbMatrix<double>()

#define __FL__ (file ? file : __FILE__), (line ? line : __LINE__)

#define F_L   (F ? F : __FILE__), (F ? L : __LINE__)
#define F_LF  (F ? F : __FILE__), (F ? L : __LINE__), __FUNCTION__
#define F_L_F (F ? F : __FILE__), (F ? L : __LINE__), (fct ? fct : __FUNCTION__)

#define FSTR  __FILE__ ": "
#define FLINE __FILE__, __LINE__
#define FL    __FILE__, __LINE__
#define FLF   __FILE__, __LINE__, __FUNCTION__
#define FCT   __FUNCTION__
#define FCL   __FUNCTION__, __LINE__

std::string getName(const std::type_info &type_id, char vflag=1);

#define PFL  shortFL(__FILE__,__LINE__,-1,""), WBL_GOT_SHORTFL
#define PF_L shortFL(F?F:__FILE__,F?L:__LINE__,-1,""), WBL_GOT_SHORTFL

#ifndef TIME
#define TIME Wb::TimeStamp('T').data
#endif

#ifndef shortFLT
#define shortFLT shortFL(__FILE__,__LINE__), Wb::TimeStamp('T').data
#define SHORT_FL shortFL(__FILE__,__LINE__) 
#endif

#ifndef PSTR
  #ifdef QS_USING_OMP
    #define PSTR Wb::ompID2Str().data
  #else
    #define PSTR "" 
  #endif
#endif

#ifndef Inf
#define Inf FP_INFINITE
#endif

#if !defined(MAX)
   template <class T> inline
   const T& MAX(const T&a, const T&b) { return ((b)>(a) ? (b) : (a)); }
#endif

   template <class T> inline
   const T& MAX3(const T&a, const T&b, const T&c) {
      if ((b)>(a))
           { return ((c)>(b) ? (c) : (b)); }
      else { return ((c)>(a) ? (c) : (a)); }
   }

#if !defined(MIN)
   template <class T> inline
   const T& MIN(const T&a, const T&b) { return ((b)<(a) ? (b) : (a)); };
#endif

   template <class T> inline
   const T& MIN3(const T&a, const T&b, const T&c) {
      if ((b)<(a))
           { return ((c)<(b) ? (c) : (b)); }
      else { return ((c)<(a) ? (c) : (a)); }
   };

#if !defined(SWAP)
   template <class T>
   inline void SWAP(T &a, T &b) { T x=a; a=b; b=x; }
#endif

#if !defined(NSTR) 
#define NSTR(a)    Wb::num2Str(a  ).data        
#endif
#if !defined(NSTRf)
#define NSTRf(a,f) Wb::num2Str(a,f).data        
#endif

#if !defined(RSTR)                              
#define RSTR(a,n) wbstring().initrep(a,n).data
#endif

#if !defined(SEC2STR)                            
#define SEC2STR(a) Wb::sec2Str(a).data
#endif

#if !defined(TSTR)
#define TSTR(a) getName(typeid(a),1).c_str()
#endif

#if !defined(sTSTR)
#define sTSTR(a) getName(typeid(a),0).c_str() 
#endif

#if !defined(SSTR)
#define SSTR(a) (a).sizeStr().data
#endif

#if !defined(SSTR_)
#define SSTR_(a) (a)->sizeStr().data
#endif

#if !defined(SSTRM)
#define SSTRM(S,r) (S).toStrf("%d","x",r,"_").data
#endif

#if !defined(SSTRM_)
#define SSTRM_(r,sd,m,sm) Wb::sizeStrM(r,sd,m,sm).data
#endif

#if !defined(QSTR)
#define QSTR(a) (a).qStr().data
#endif

#if !defined(QSTR_)
#define QSTR_(a) (a)->qStr().data
#endif

#if !defined(STR)
#define STR(a) (a).toStr().data
#endif

#if !defined(STR_)
#define STR_(a) (a)->toStr().data
#endif

#if !defined(STR2)
#define STR2(a,b) (a).toStr(b).data
#endif

#if !defined(STR2_)
#define STR2_(a,b) (a)->toStr(b).data
#endif

#if !defined(cSTR)
#define cSTR(a) Wb::char2Str(a).data
#endif

#if !defined(I2STR)
#define I2STR(a) Wb::int2Str(a).data
#endif

#if !defined(RATS)
#define RATS(a) Wb::rat2Str(__FILE__,__LINE__,a).data
#endif

#if !defined(STRREP)
#define STRREP(str,n) wbstring().initrep(str,n).data
#endif

#if !defined(BITS)
#define BITS(a) Wb::bits(a).data
#endif

#if !defined(POW)

   template <class T>
   inline T POW(const T &a, const unsigned &n) {
      T x = n ? a : 1;
      for (unsigned i=1; i<n; i++) x*=n;
      return x;
   }

#endif
#if !defined(SGN)

   template <class T>
   inline int SGN(const T &a) {
      return ((T(0)<a)-(a<T(0))); 
   }

#endif
#if !defined(NUMCMP)

   template <class T>
   inline char NUMCMP(const T &a, const T &b) {
      if (a<b) return -1; else
      if (a>b) return +1;
      return 0;
   }

   template <class Ta, class Tb>
   inline char NUMCMP(const Ta &a, const Tb &b) {
      if (a<b) return -1; else
      if (a>b) return +1;
      return 0;
   }

#endif

#define ISCOMPLX_(a) WbUtil<a>::isComplex()
#define ISREAL(a) WbUtil<a>::isReal()

#ifdef TST_MAC
   #include "Include/nomex.hh"

   #include <Accelerate/Accelerate.h>

#elif defined(NOMEX)
   #include "nomex.hh"
#else

   #include <mex.h>
   #include <mat.h>

   #define QS_VERSION 4.0
   #define QS_VERSION_SUB   0
   #define QS_VERSION_SUB_ ""     

   #ifdef MEX_MLVER
      #define MEX_MLVER_STR PP_STRFY(MEX_MLVER)
   #else
      #define MEX_MLVER_STR "(unknown)"
   #endif

   #ifdef HOST_NAME
      #define HOST_NAME_STR PP_STRFY(HOST_NAME)
   #else
      #define HOST_NAME_STR "(host)"
   #endif

#define mxIsNumChar(a) (mxIsNumeric(a) || mxIsChar(a) || mxIsLogical(a))

#endif

#ifdef MAIN 
   #undef printf  
#endif

#ifndef PROG
#ifdef MATLAB_MEX_FILE
   #define PROG mexFunctionName()
#else
   #define PROG "(?prog?)" 
#endif
#endif

   enum WBL_COLOR_SCHEME { WLC_OFF,
       WLC_DARK,
   NUM_WBL_COLOR_SCHEME };

namespace Wb {
   class myname__ {
     public:
       myname__() { unsigned i=0, k=0;
          #ifdef PROG
             snprintf(data,60,"%s",PROG);
          #elif defined MATLAB_MEX_FILE
             snprintf(data,60,"%s",mexFunctionName());
          #else
             snprintf(data,60,"(?myname?)");
          #endif

           for (i=0; data[i]; ++i) { if (data[i]=='/') k=i+1; }
           if (k && data[k]) {
              for (i=0; data[k]; ++k, ++i) { data[i]=data[k]; }
              data[i]=0;
           }

          #ifdef PROG_TAG 
             snprintf(data+60,4,"%.3s",PROG_TAG);
          #else
             snprintf(data+60,4,"%.3s",data);
          #endif
       };

       char data[64];

     private:
   };

   myname__ myName; 

#define myname Wb::myName.data
#define mytag  Wb::myName.data+60 

   WBL_COLOR_SCHEME useCol = WLC_DARK; 

   int get_WB_VERBOSE(const char *F=0, int L=0);
   int got_DBSTOP(const char *F=0, int L=0);
   int got_DESKTOP();
   int is_DEPLOYED(const char *F=0, int L=0);

   int envVRB= 0; 
   int envDKT=-1; 

   int envDBG= 0; 

   int my_caller_tid=0; 

}; 

void dbstop(const char* file, int line);
void ExitMsg(const char* s="", char xflag=0);
void doflush();

#ifdef LOAD_CGC_QSPACE
#ifndef QS_SKIP_MPFR 
   #define QS_USING_MPFR
#endif
#endif

   #define QS_ITAG_ "__"  

#ifdef  __WBDEBUG__
#define       __WB_MEM_CHECK__
#undef WB_SKIP_ASSERT
#elif defined __WB_MEM_CHECK__
#define __WBDEBUG__
#endif

namespace Wb {
class ARGV { 
  public:
    ARGV(va_list *vl=NULL) : args(vl) {};
   ~ARGV() { if (args) { va_end(*args); args=NULL; }};

  private:
    va_list *args;
};
};

namespace wblog { 
   const int SLEN=256;  
};

#ifdef MATLAB_MEX_FILE
#define PRINTF(...) \
 { if (omp_get_level()==0) \
        { printf(__VA_ARGS__); } \
   else { fprintf(stdout,__VA_ARGS__); } \
 }
#else
#define PRINTF printf
#endif

#include "wblibx.h" 

#include "wblog.h"

#ifdef QS_USING_OMP
#include "wbomp.hh" 
#endif

#include "memtrack.hh"

#include "util.hh"
#include "memlib.hh"
#include "mxlib.hh"

#ifdef QS_USING_MPFR
#include <gmp.h>
#include <mpfr.h> 
#include "wbmpfr.hh" 
#endif

#include "wbcomplex.hh"
#include "gathered.hh"
#include "wbsort.hh"

#include "wbvector.hh"
#include "wbstring.hh"
#include "wbbitset.hh"
#include "wbio.h"     

#include "wbperm.hh"  
#include "wbindex.hh" 

#include "wbclock.hh"
#include "wbMatrix.hh"

#include "mexlib.hh"

#include "wbMatIO.hh"

#include "wbarray.hh"
#include "wbsparray.hh"

#include "wbutil.hh"

#include "wbblas.hh"
#include "wbMatrix_blas.hh"
#include "wbarray_blas.hh"

#include "wbopts.hh"

#ifdef DBG_GCX_LOCKS
   Wb::tmpFile LKF(FL,"-lks"); 
#endif

#ifdef DBG_QSX_BUF
   Wb::tmpFile BFF(FL,"-buf"); 
#endif

#ifdef LOAD_CGC_QSPACE
#include "clebsch.hh"
#include "QSpace_aux.hh"
#include "QSpace.hh"
#include "mpsortho.hh" 
#endif

#include "wblog.c"

#ifdef QS_USING_MPFR
#include "wbmpfr.cc"
#endif

#ifdef QS_USING_OMP
#include "wbomp.cc"   
#endif

#include "wbsort.cc"
#include "wbstring.cc"
#include "wbbitset.cc"
#include "memlib.cc"
#include "wbcomplex.cc"
#include "gathered.cc"
#include "util.cc"

#ifdef LOAD_CGC_QSPACE
#include "QSpace_aux.cc"
#include "QSpace.cc"
#endif

#include "wbvector.cc"
#include "wbperm.cc"
#include "wbindex.cc"
#include "wbMatrix.cc"
#include "wbarray.cc"
#include "wbsparray.cc"

#include "wbclock.cc"
#include "mexlib.cc"
#include "matlib.cc"  
#include "wbMatrix_blas.cc"
#include "wbarray_blas.cc"

#include "wbio.c"
#include "wbmsg.hh"

#ifdef LOAD_CGC_QSPACE
#include "clebsch.cc"
#include "clebsch_aux.cc"
#include "clebsch_io.cc"  
#include "clebsch_old.cc"
#include "mpsortho.cc"    
#endif

namespace Wb {

   class gpara__ {
     public:
       gpara__() : x(0) {
         #ifdef __WBDEBUG__
           wblog(FL,"TST %s (%s)",FCT,myname);
         #endif
           wb_srand(); 

         #ifdef __WB_MPFR_HH__
           gmp_randinit_default(Wb::wb_rstate);
         #endif

          envDKT=got_DESKTOP();
          my_caller_tid=omp_get_thread_num(); 

          init();
       };

       void init() {
          memset(str,0,STRLEN+1);   
          envVRB=get_WB_VERBOSE();  
          envDBG=got_DBSTOP();

       #ifdef LOAD_CGC_QSPACE
       #ifdef QS_USING_OMP
          Wb::GetNumThreads(FL,OMP_NUM_THREADS,"OMP_NUM_THREADS");
          Wb::GetNumThreads(FL,QSP_NUM_THREADS,"QSP_NUM_THREADS");

          sp_num_threads=MAX(OMP_NUM_THREADS,QSP_NUM_THREADS);
       #endif
       #endif

          my_caller_tid=omp_get_thread_num(); 

       };

       char x;

    private:
   };

   gpara__ gpara;

class CleanUp {
  public:
    CleanUp() : i(0) {
       gpara.init(); 
    };

   void Check() {
      wblog::check_ERR_pending();
   };

   ~CleanUp() {
       wblog::myIO.clear();
       if (gwb_Profs.size()) { Wb::save_and_clear_Profiling(); }
    };

    int i;
};

}; 

#ifdef MAIN
   #ifdef MATLAB_MEX_FILE
      #warning got MAIN with MATLAB_MEX_FILE
      #pragma message "\n\n"  \
      "   Don't define MAIN with MEX files, this kills the whole matlab seesion!\n" \
      "   See ExitMsg() in wblib.h (!)\n" \
""
   #endif

   void ExitMsg(const char* s, char xflag) {
      if (s && s[0])
           { printf("\n%s !? %s\n\n",shortFL(FL),s); }
      else { printf("\n"); }

      if (Wb::envDBG) {
         dbstop(FL); if (xflag) { exit(-1); }
      }
      else {

         if (xflag) { 
            exit(-1); 
         }
         else {
            mexErrMsgIdAndTxt("Wb:MEX:wblog",s);
         }
      }
   }

   void doflush() { fflush(0); };

#else

   #ifdef DBSTOP
      void ExitMsg(const char* s, char xflag) {
         if (s && s[0]) printf("%s >%s<\n",shortFL(FL),s ? s:"");
         dbstop(FL);
         printf("\n%s dbstop() and on we go ...\n\n", shortFL(FL));
         mexErrMsgIdAndTxt("Wb:ERR:mex",s);
      }
   #else
      void ExitMsg(const char* s, char xflag __attribute__ ((unused))) {
         printf("\n");
         if (Wb::is_DEPLOYED()) { 
             Wb::print_backtrace(FL,"from deployed code"); printf("\n"); }

         mexErrMsgIdAndTxt("Wb:ERR:mex",s);
         fprintf(stdout,"%s TST %s()\n",shortFL(FL),FCT); 
      }
   #endif

   void doflush() {
       fflush(0); 

       #if 0

       if (!omp_in_parallel()) {
          mexEvalString("pause(0);");
       }
       #endif

   };

#endif

#endif

