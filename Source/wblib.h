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

#ifndef __WB_LIB_H__
#define __WB_LIB_H__

// #define QS_FERM_ACTIVE

#define PP_STR__(a) #a
#define PP_STRFY(a) PP_STR__(a)

#ifndef QS_USING_OMP
#  ifdef _OPENMP
#     define QS_USING_OMP 1
#  else
#     warning compiling QSpace without openMP
#  endif
#elif QS_USING_OMP<=0
   #undef QS_USING_OMP
#endif

#if __linux__
   #define QS_UNUSED_VAR __attribute__ ((unused))
#elif __unix__
   #warning got unix
#elif __APPLE__
   #define QS_UNUSED_VAR __attribute__ ((unused))
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

#include <locale.h> 
#include <wchar.h>

#include <typeinfo>

#include <sys/ipc.h>
#include <sys/msg.h>
#include <sys/types.h>
#include <sys/stat.h>

#if __linux__
#include <sys/sysinfo.h>     
#include <sys/syscall.h>     
#elif __APPLE__
#include <pthread.h>         
#endif

#include <map>
#include <queue>             
#include <unordered_map>     
#include <algorithm>         
#include <vector>            
#include <initializer_list>  

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

#ifdef QS_USING_HPTT
#include "hptt.h" 
#undef I  
#endif

   unsigned ARG_CHECK=1;

#define STRLEN  1023
   char str[STRLEN+1];

#define sprintf_str(...) snprintf(str,STRLEN,__VA_ARGS__)

#if 1
#   define  widx_t   size_t  
#   define  widx_ts  long    
#   define  wperm_t  size_t
#   define  wperm_ts long    
#else
#   define  widx_t   unsigned
#   define  widx_ts  int
#   define  wperm_t  unsigned
#   define  wperm_ts int
#endif

#if 1
#   define  SPIDX_T  size_t
#   define sSPIDX_T  long
#else
#   define  SPIDX_T  unsigned
#   define sSPIDX_T  int
#endif

#define WBINDEX  wbvector<widx_t>
#define WBPERM   wbvector<wperm_t>
#define WBIDXMAT wbMatrix<widx_t>

#define  MVEC        wbvector<int8_t> 
#define cMVEC  const wbvector<int8_t> 

#define  UVEC        wbvector<unsigned>
#define cUVEC  const wbvector<unsigned>
#define  DMAT        wbMatrix<double>
#define cDMAT  const wbMatrix<double>
#define  UMAT        wbMatrix<unsigned>
#define cUMAT  const wbMatrix<unsigned>
#define  IMAT        wbMatrix<int>
#define cIMAT  const wbMatrix<int>
#define cTMAT  const wbMatrix<T>
#define cMEX   const mxArray
#define cUINT  const unsigned
#define cCHAR  const char

#define DMAT0   wbMatrix<double>()

#define __FL__ (file ? file : __FILE__), (line ? line : __LINE__)

#define F_L   (F ? F : __FILE__), (F ? L : __LINE__)
#define FL_   (F ? F : NULL), (F ? L : 0)  
#define FL_q(x) (x?__FILE__:NULL), (x?__LINE__:0)  
#define F_LF  (F ? F : __FILE__), (F ? L : __LINE__), __FUNCTION__
#define F_L_F (F ? F : __FILE__), (F ? L : __LINE__), (fct ? fct : __FUNCTION__)

#define FSTR  __FILE__ ": "
#define FLINE __FILE__, __LINE__
#define FL    __FILE__, __LINE__
#define FLF   __FILE__, __LINE__, __FUNCTION__
#define FCT   __FUNCTION__
#define FCL   __FUNCTION__, __LINE__

#define C_FCT  Wb::class_function(__PRETTY_FUNCTION__,1) 
#define pFCT  Wb::class_function(__PRETTY_FUNCTION__,0)

std::string getName(const std::type_info &type_id, char vflag=1);

#define PFL  shortFL(    __FILE__,    __LINE__,-1,""), WBL_GOT_SHORTFL
#define PF_L shortFL(F?F:__FILE__,F?L:__LINE__,-1,""), WBL_GOT_SHORTFL

#ifndef TIME
#define TIME Wb::TimeStamp('T').data
#endif

#ifndef shortFLT
#define shortFLT shortFL(__FILE__,__LINE__), Wb::TimeStamp('T').data
#define SHORT_FL shortFL(__FILE__,__LINE__) 
#endif

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

#ifndef PSTR
#  ifdef QS_USING_OMP
#    define PSTR Wb::ompID2Str().data
#  else
#    define PSTR "" 
#  endif
#endif

#ifndef Inf
#  define Inf FP_INFINITE
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
#define  TSTR(a) getName(typeid(a),1).c_str() 
#endif

#if !defined(sTSTR)
#define sTSTR(a) getName(typeid(a),0).c_str() 
#endif

#if !defined(SSTR)
#define SSTR(a) (a).sizeStr().data
#endif

#if !defined(SSTR_)
#define SSTR_(a) (a ? (a)->sizeStr().data : "(null)")
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
#define QSTR_(a) (a ? (a)->qStr().data : "(null)")
#endif

#if !defined(STR)
#define STR(a) (a).toStr().data
#endif

#if !defined(STR_)
#define STR_(a) (a ? (a)->toStr().data : "(null)")
#endif

#if !defined(STR2)
#define STR2(a,b) (a).toStr(b).data
#endif

#if !defined(STR2_)
#define STR2_(a,b) (a ? (a)->toStr(b).data : "(null)")
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

#if !defined(BITS_) 
#define BITS_(a,b) Wb::bits(a,b).data
#endif

#if !defined(UNSET_BINARY) 
#define UNSET_BINARY(a,b) ((a) &= ~(b))
#endif

#if !defined(ABS)
#define ABS(a) Wb::abs(a)
#endif

#if !defined(NORM)
#define NORM(a) Wb::norm(a)
#endif

#if !defined(NORM2)
#define NORM2(a) Wb::norm2(a)
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
      return ( (a>T(0)) - (a<T(0)) ); 
   }
#endif

#if !defined(NUMCMP)
   template <class T>
   inline char NUMCMP(const T &a, const T &b) {
      if (a<b) { return -1; } else
      if (a>b) { return +1; } else { return 0; }
   };

   template <class Ta, class Tb>
   inline char NUMCMP(const Ta &a, const Tb &b) {
      if (a<b) { return -1; } else
      if (a>b) { return +1; } else { return 0; }
   };
#endif

#define ISREAL(x)    WbUtil<x>::isReal()
#define ISCOMPLX_(x) WbUtil<x>::isComplex()

#ifdef TST_MAC
   #include "Include/nomex.hh"

   #include <Accelerate/Accelerate.h>

#elif defined(NOMEX)

   #include "nomex.hh"

#else

   #include <mex.h>
   #include <mat.h>

   #define QS_VERSION 4.1 

   #define QS_VERSION_SUB   0
   #define QS_VERSION_SUB_ ""     

#   ifdef MEX_MLVER
#      define MEX_MLVER_STR PP_STRFY(MEX_MLVER)
#   else
#      define MEX_MLVER_STR "(unknown)"
#   endif

#   ifdef HOST_NAME
#      define HOST_NAME_STR PP_STRFY(HOST_NAME)
#   else
#      define HOST_NAME_STR "(host)"
#   endif

#if 1
#  include <blas.h>    
#  include <lapack.h>  
#endif

#define mxIsNumChar(a) (mxIsNumeric(a) || mxIsChar(a) || mxIsLogical(a))

#endif

#ifdef MAIN 
#  undef printf  
#endif

#ifndef PROG
#  ifdef MATLAB_MEX_FILE
#    define PROG mexFunctionName()
#  else
#    define PROG "(?prog?)" 
#  endif
#endif

   enum WBL_COLOR_SCHEME { WLC_OFF,
       WLC_DARK,
   NUM_WBL_COLOR_SCHEME };

#ifdef LD_CLEBSCH_QS
#  ifndef QS_SKIP_MPFR 
#    define QS_USING_MPFR
#  endif
#endif

   #define QS_ITAG_ "__"  

#ifdef __WBDEBUG__
#  define  __WB_MEM_CHECK__
#  undef WB_SKIP_ASSERT
#elif defined __WB_MEM_CHECK__
#  define __WBDEBUG__
#endif

#ifdef MATLAB_MEX_FILE
#define PRINTF(...) \
 { if (omp_get_level()==0) \
        { printf(__VA_ARGS__); } \
   else { fprintf(stdout,__VA_ARGS__); } \
 }
#else
#  define PRINTF printf
#endif

#include "wblib.hh"

#include "wblog.h"

#ifdef QS_USING_OMP
#include "wbomp.hh" 
#endif

#include "memtrack.hh"

#include "util.hh"
#include "memlib.hh"
#include "mxlib.hh"

#include "wbcomplex.hh" 

#ifdef QS_USING_MPFR
#include <gmp.h>
#include <mpfr.h> 
#include "wbmpfr.hh" 
#endif

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

#ifndef refblas_h 
#  include "wbblas.hh"
#endif

#include "wbMatrix_blas.hh"
#include "wbarray_blas.hh"

#include "wbopts.hh"

#ifdef DBG_GCX_LOCKS
   Wb::tmpFile LKF(FL,"-lks"); 
#endif

#ifdef DBG_QSX_BUF
   Wb::tmpFile BFF(FL,"-buf"); 
#endif

#ifdef LD_CLEBSCH_QS
#include "clebsch.hh"
#include "QSpace_aux.hh"
#include "QSpace.hh"
#include "mpsortho.hh"    
#endif

#include "wblib.cc" 
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

#ifdef LD_CLEBSCH_QS
#include "QSpace_aux.cc"
#include "QSpace.cc"
#endif

#include "wbvector.cc"
#include "wbperm.cc"
#include "wbindex.cc"
#include "wbMatrix.cc"
#include "wbarray.cc"
#include "wbsparray.cc"

#include "wbopts.cc"
#include "wbclock.cc"

#include "mexlib.cc"
#include "matlib.cc"  
#include "wbMatrix_blas.cc"
#include "wbarray_blas.cc"

#include "wbio.c"
#include "wbmsg.hh"

#ifdef LD_CLEBSCH_QS
#include "clebsch.cc"
#include "clebsch_aux.cc"
#include "clebsch_io.cc"  
#include "clebsch_old.cc"
#include "mpsortho.cc"    
#endif

#ifdef MAIN
   #ifdef MATLAB_MEX_FILE
      #warning got MAIN with MATLAB_MEX_FILE
      #pragma message "\n\n"  \
      "   Don't define MAIN with MEX files, this kills the whole matlab seesion!\n" \
      "   See ExitMsg() in wblib.h (!)\n\n" \
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

#elif defined(MATLAB_MEX_FILE)

   #ifdef DBSTOP
      void ExitMsg(const char* s, char xflag) {
         if (s && s[0]) printf("%s %s() :: %s\n",shortFL(FL),FCT,s ? s:"");
         dbstop(FL);
         printf("\n%s %s() and on we go ...\n\n",shortFL(FL),FCT);

         mexErrMsgIdAndTxt("Wb:ERR:mex",s);
      }
   #else
      void ExitMsg(const char* s, char xflag QS_UNUSED_VAR) {
         printf("\n");
         if (Wb::is_DEPLOYED()) { 
             Wb::print_backtrace(FL,"from deployed code"); printf("\n"); }

         mexErrMsgIdAndTxt("Wb:ERR:mex",s);
         fprintf(stdout,"%s TST %s()\n",shortFL(FL),FCT); 
      }
   #endif

   void doflush() {
     #if defined(__OMP_H) || defined(QS_USING_OMP)
      if (!Wb::omp_parallel())
     #endif
      { fflush(0); } 

      #if 0

      if (!omp_in_parallel()) {
         mexEvalString("pause(0);");                 
      }
      #endif

   };

#endif

#endif

