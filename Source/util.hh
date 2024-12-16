/* ---------------------------------------------------------------------
 * Project : QSpace tensor library (v4.0)
 * Class   : QSpace utility routines
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

#ifndef __WB_UTIL_HH__
#define __WB_UTIL_HH__

/* -------------------------------------------------------------------- */
namespace Wb {
/* -------------------------------------------------------------------- */

   int Rational(
      double &x,          
      long &P, long &Q,   
      double *dx=NULL,    
      unsigned *nmax=NULL,
      wbvector<double>* aa=NULL, 
      unsigned niter=6,   
      long pqmax=-1,      
      double eps=1E-8,    
      double eps2=1E-14,  
      char vflag=0
   );

   template<class T>
   double FixRational(  
      const char *F, int L, T *d, unsigned n, char rflag=0,
      unsigned niter=6,
      long pqmax=-1,    
      double eps=1E-8,  
      double eps2=1E-14,
      double *rz=NULL, double *ra=NULL, char vflag=0
   );

   wbstring rat2Str( 
      const char *F, int L, double d, unsigned niter=6,
      long pqmax=99999, double eps=1E-8, double eps2=1E-14, char vflag=0
   );

   template<class T>
   double SkipZeros(T *d, unsigned n, double eps=1E-14){ return 0; };

   void ResSummary(const char *F, int L, const char *istr=NULL);

#ifndef __WB_MEM_TRACK_HH__
   void MemStat(const char *F=0, int L=0, char lflag=0);
#else
   void MemStat(const char *F, int L, char lflag);
#endif

#ifndef __APPLE__
   long   getCpuInfo (const char *tag, const char *f="/proc/cpuinfo");
   size_t getProcSize(const char *tag, const char *f="/proc/meminfo");
#endif
   size_t getProcSize(const char *tag, const pid_t &p);
   size_t getProcSize(const char *tag, const char *f);

   wbstring hostname(unsigned len=16);
   wbstring hostid(char pflag=0);

   wbstring size2Str(double n, char p=4);
   const char* filename(const char *s);

   size_t getFileSize(const char* f); 
   double getFileTime(const char* f, char w='m'); 

   double GET_mtime(const char *F, int L, const char *f, const struct stat &S);
   double GET_atime(const char *F, int L, const char *f, const struct stat &S);
   double GET_ctime(const char *F, int L, const char *f, const struct stat &S);

   bool fexist(const char *f, char type='*');
   bool isFile(const char *f);
   bool isDir(const char *f);

   int mkdir( 
      const char *F, int L, const char *path,
      char priv=0 
   ){
      mode_t p=(S_IRWXU 
              | S_IXGRP 
              | S_IXOTH 
      ); 
      if (!priv) p |= (
                S_IRGRP 
              | S_IROTH 
      ); 

      int e;

     #pragma omp critical (manage_WB_FILE_IO) 
      { e=::mkdir(path,p); }

      if (e && F) wblog(F_L,
         "ERR %s() failed to create directory (e=%d)\n%s",FCT,e,path);
      return e;
   };

   template <class TI> 
   TI sub2ind(const TI *s, const TI *I, unsigned n);

   template <class T1, class T2> inline
   void ind2sub(T1 k, const T1 *s, T2 *I, unsigned n);

   int atoi(const char *s, int &iout, unsigned n, char white=0);

   int isUIntString(const char *s, char white=1); 

   template <class T>
   size_t findfirst_sorted(const char *F, int L,
      const T* d0, const T* dd, size_t m, size_t N, size_t M,
      char lex=1 
   );

   template <class T>
   size_t findlast_sorted(const char *F, int L,
      const T* d0, const T* dd, size_t m, size_t N, size_t M,
      char lex=1 
   );

   template <class T0, class T> inline
   size_t findfirst_sorted(const char *F, int L,
      const T0* d0, const T* dd, size_t m, size_t N, size_t M, char lex=1
   ){
      T d2[m];
      for (unsigned i=0; i<m; ++i) d2[i]=T(d0[i]);
      return findfirst_sorted(F,L,d2,dd,m,N,M,lex);
   };

   template <class T0, class T> inline
   size_t findlast_sorted(const char *F, int L,
      const T0* d0, const T* dd, size_t m, size_t N, size_t M, char lex=1
   ){
      T d2[m];
      for (unsigned i=0; i<m; ++i) d2[i]=T(d0[i]);
      return findlast_sorted(F,L,d2,dd,m,N,M,lex);
   };

   int dstrlen_utf8(const char *s);

   template <class T>
   wbstring bits(const T &x, char compact=1); 

class IOstat { 
  public:

    IOstat()
     : nread(0), nwrite(0), sread(0), swrite(0), aux1(0), aux2(0) {};

    IOstat& init() {
       nread=nwrite=0; sread=swrite=aux1=aux2=0;
       return *this;
    }

    void gotread (size_t s=0) { nread +=1; sread +=s; }
    void gotwrite(size_t s=0) { nwrite+=1; swrite+=s; }

    bool isEmpty() {
       return (!nread && !nwrite && !sread && !swrite && !aux1 && !aux2);
    };

    void print(const char *istr,
       const char *hstr=0, const char *xstr1=0, const char *xstr2=0);

    unsigned nread, nwrite;
    size_t sread, swrite, aux1, aux2;

  protected:
  private:
};

template <class T, class T2>
class iterator { 
  public:

    iterator(T* p_) : i(0), p(p_) {
       if (!p_) wblog(FL,"ERR %s() got p=NULL !?",FCT);
       n=p->numel();
    };

    iterator& operator++() { ++i; return *this; }
    bool end() { return i>=n; }

    size_t iter () const { return i; }
    size_t numel() const { return n; }

    T2& operator[](size_t k) {
       if (!p) wblog(FL,"ERR %s() object not yet initialized !?",FCT);
       return (*p)[k];
    }
    T2& at() { 
       if (i>=n) wblog(FL,
          "ERR %s() index out of bounds (%d/%d) !?",FCT,i,n);
       return (*p)[i];
    };

    size_t i,n;
    T* p;

  protected:
  private:
};

template <class T, class T2>
class iterator<T,T2*> { 
  public:

    iterator(T* p_) : i(0), p(p_) {
       if (!p_) wblog(FL,"ERR %s() got p=NULL !?",FCT);
       n=p->numel();
    };

    iterator& operator++() { ++i; return *this; }
    bool end() { return i>=n; }

    size_t iter () const { return i; }
    size_t numel() const { return n; }

    T2* operator[](size_t k) {
       if (!p) wblog(FL,"ERR %s() object not yet initialized !?",FCT);
       return (*p)[k];
    }
    T2* at() { 
       if (i>=n) wblog(FL,
          "ERR %s() index out of bounds (%d/%d) !?",FCT,i,n);
       return (*p)[i];
    };

    size_t i,n;
    T* p;

  protected:
  private:
};

template <class T, class T2>
class citerator { 
  public:

    citerator(const T* p_) : i(0), p(p_) {
       if (!p_) wblog(FL,"ERR %s() got p=NULL !?",FCT);
       n=p->numel();
    };

    citerator& operator++() { ++i; return *this; }
    bool end() { return i>=n; }

    size_t iter () const { return i; }
    size_t numel() const { return n; }

    const T2& operator[](size_t k) {
       if (!p) wblog(FL,"ERR %s() object not yet initialized !?",FCT);
       return (*p)[k];
    }
    const T2& at() { 
       if (i>=n) wblog(FL,
          "ERR %s() index out of bounds (%d/%d) !?",FCT,i,n);
       return (*p)[i];
    };

    size_t i,n;
    const T* p;

  protected:
  private:
};

template <class T, class T2> 
class citerator<T,T2*> { 
  public:

    citerator(const T* p_) : i(0), p(p_) {
       if (!p_) wblog(FL,"ERR %s() got p=NULL !?",FCT);
       n=p->numel();
    };

    citerator& operator++() { ++i; return *this; }
    bool end() { return i>=n; }

    size_t iter () const { return i; }
    size_t numel() const { return n; }

    const T2* operator[](size_t k) {
       if (!p) wblog(FL,"ERR %s() object not yet initialized !?",FCT);
       return (*p)[k];
    }
    const T2* at() { 
       if (i>=n) wblog(FL,
          "ERR %s() index out of bounds (%d/%d) !?",FCT,i,n);
       return (*p)[i];
    };

    size_t i,n;
    const T* p;

  protected:
  private:
};

class rand { 

  public:

    rand() {};

    double get_dfac() { return dfac; }

    template <class T>
    T& rand_(T &x); 

    template <class T>
    T& rand_(T &x, double fac); 

    template <class T>
    T& rands(T &x); 

    template <class T>
    T& randb(T &x); 

    template <class T>
    void scale_rand(T *x, size_t n, double fac) {
       for (size_t i=0; i<n; ++i) { x[i]*=fac; }
    }

  protected:
  private:

    static double dfac;
};

double rand::dfac = (1.-1./double(long(1)<<52)) / double(RAND_MAX);

template<> inline 
double& rand::rand_(double &x) { return (x=dfac*::rand()); };

template<> inline 
int& rand::rand_(int &x) { return (x=(100*dfac)*::rand()); };

template<> inline 
unsigned& rand::rand_(unsigned&x) { return (x=(100*dfac)*::rand()); };

template<> inline 
long& rand::rand_(long &x) { return (x=(100*dfac)*::rand()); };

template<> inline 
size_t& rand::rand_(size_t &x) { return (x=(100*dfac)*::rand()); };

template<> inline
int& rand::rand_(int&x, double fac) {
    return (x=(fac*dfac)*::rand());
};

template<> inline
unsigned& rand::rand_(unsigned &x, double fac) {
    return (x=(fac*dfac)*::rand());
};

template<> inline
long& rand::rand_(long &x, double fac) {
    return (x=(fac*dfac)*::rand());
};

template<> inline
size_t& rand::rand_(size_t &x, double fac) {
    return (x=(fac*dfac)*::rand());
};

template<> inline
double& rand::randb(double &x) { return (x=(2*dfac)*::rand()-1.); }

template<> inline
int& rand::randb(int &x) { return (x=(199*dfac)*::rand()-99); }

template<> inline
double& rand::rands(double &x) {
   x=(2*dfac*::rand()-1);
   if (fabs(x)>=1) wblog(FL,"ERR %s() got x=%.4g (1%+.4g) !?",FCT,x,x-1);
   return (x=atanh(x));
};

template<> inline
int& rand::rands(int &x) {
   double d=dfac*::rand();
   if (fabs(d)>=1) wblog(FL,"ERR %s() got d=%.4g (1%+.4g) !?",FCT,d,d-1);
   return (x=(100*atanh(d)+0.5)); 
};

template<> inline
void rand::scale_rand(int *x, size_t n, double fac) {
   fac/=double(99);
   for (size_t i=0; i<n; ++i) { x[i]=((double(x[i]))*fac); }
};

template <class T,
typename std::enable_if< std::is_floating_point<T>::value , T>::type* = nullptr >
bool is_finite(const T *x, size_t n) {
   for (size_t i=0; i<n; ++i) { if (!::isfinite(x[i])) { return 0; }}
   return 1;
};

bool is_finite(const wbcomplex *x, size_t n);

template <class T, 
typename std::enable_if< is_integral<T>::value , T>::type* = nullptr >
bool is_finite(const T *x, size_t n) { return 1; }

}; 

#ifdef __APPLE__

namespace Mac {
   double getMemSize(const char *tag);
   size_t getProcMemSize(
      pid_t p=0,   
      char res=1   
   );
};

#endif

   void wb_srand();

template <class T>
class WbUtil { 
  public:

     static void adjust_tnorm(char &tnorm);

     static constexpr bool isPOD() {
        return std::is_trivially_copyable<T>::value;
     };

     static bool isFloat();  
     static bool isInt();    

     static constexpr char isComplex() { return 0; };

     static char hasConj() { return isComplex(); };
     static char isReal() {
        char q=isComplex(); if (q>=0) { q=(q? 0 : 1);  }
        return q;
     };

     T eps(); 

  protected:
  private:

};

template <class T> inline
void WbUtil<T>::adjust_tnorm(char &tnorm){ tnorm=0; };

template <> inline
void WbUtil<wbcomplex>::adjust_tnorm(char &tnorm __attribute__ ((unused))){};

template <> inline
constexpr char WbUtil<wbcomplex>::isComplex() { return 1; };

template <> inline
constexpr char WbUtil<void>::isComplex() { return -1; }; 

#ifdef MATLAB_MEX_FILE
template <> inline
constexpr char WbUtil<mxComplexDouble>::isComplex() { return 2; };
#endif

#ifdef _COMPLEX_H  
template <> inline
constexpr char WbUtil<complex<double> >::isComplex() { return 4; };
template <> inline
constexpr char WbUtil<complex<float > >::isComplex() { return 8; };
#endif

template <class T> inline
bool WbUtil<T>::isFloat(){ return 0; };

template <> inline bool WbUtil<double>::isFloat(){ return 1; };
template <> inline bool WbUtil<long double>::isFloat(){ return 1; };
template <> inline bool WbUtil<float>::isFloat(){ return 1; };
template <> inline bool WbUtil<wbcomplex>::isFloat(){ return 1; };

template <class T> inline
bool WbUtil<T>::isInt(){ return 0; }

template <> inline bool WbUtil<int>::isInt(){ return 1; }
template <> inline bool WbUtil<char>::isInt(){ return 1; }
template <> inline bool WbUtil<long>::isInt(){ return 1; }

template <> inline bool WbUtil<unsigned>::isInt(){ return 1; }
template <> inline bool WbUtil<unsigned char>::isInt(){ return 1; }
template <> inline bool WbUtil<unsigned long>::isInt(){ return 1; } 

template <> inline double WbUtil<double>::eps(){ return DBL_EPSILON; }
template <> inline float  WbUtil<float >::eps(){ return FLT_EPSILON; }

template <> inline
long double WbUtil<long double>::eps(){ return LDBL_EPSILON; }

#include <sys/time.h>
#include <sys/resource.h>

class cpu_time {

  public:

    cpu_time() { tcpu=time('c'); tsys=time('s'); };

    void init() { tcpu=time('c'); tsys=time('s'); };

    cpu_time since() {
       cpu_time dt; dt.tcpu-=tcpu; dt.tsys-=tsys;
       return dt;
    };

    static double time(char flag) {
       struct rusage result;
       getrusage(RUSAGE_SELF,&result);
       if (flag=='c') { 
          return (
             double(result.ru_utime.tv_sec)
           + 1e-6 * result.ru_utime.tv_usec
          );
       } else { 
          return (
             double(result.ru_stime.tv_sec)
           + 1e-6 * result.ru_stime.tv_usec
          );
       }
    };

    wbstring toStr(char flag='c');
    double tcpu, tsys; 

  protected:
  private:
};

  cpu_time gCPUTime;

namespace wbsys {

    wbstring MemTot2Str();
    wbstring MemFree2Str();
    wbstring SwapTot2Str();
    wbstring SwapFree2Str();

    size_t getMemTot();
    size_t getMemFree();
    size_t getSwapTot();
    size_t getSwapFree();

    char checkSwapSpace(const char *F=NULL, int L=0);

    int getNumCores();      
    int getCacheLineSize(); 
};

class wbtop { 

  public:

    wbtop() {
       pid=getpid(); ppid=getppid();
       sysmem=wbsys::getMemTot(); 
    };

    wbstring VmSize2Str() const;
    size_t getVmSize() const { return Wb::getProcSize("VmSize",pid); };

    wbstring VmPeak2Str() const;
    size_t getVmPeak() const { return Wb::getProcSize("VmPeak",pid); };

    int runningLarge(const char *F=NULL, int L=0,
       double th0=0.80, 
       double fac=1.2   
    );

    void MemStat(const char *F=0, int L=0, char lflag=0) const;

    pid_t pid, ppid;
    size_t sysmem;

  protected:
  private:
};

#endif

