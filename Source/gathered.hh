/* ---------------------------------------------------------------------
 * Project : QSpace tensor library (v4.0 pre-release)
 * Class   : QSpace additional routines
 *
 * Copyright 2022 Andreas Weichselbaum
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

#ifndef __WB_GATHERED_HH__
#define __WB_GATHERED_HH__

/* -------------------------------------------------------------------- */
/* -------------------------------------------------------------------- */

namespace Wb {
   const char* basename(const char *data, char c='/');

   void print_backtrace(
      const char *F=NULL, int L=0, const char *istr=NULL);

   int GetNumThreads(const char *F, int L, int &n, const char *name);

   template<class T>
   int GetEnv (const char *F, int L, const char *name, T& val);

   int EnvIsSet(const char *F, int L, const char *name);

   wbstring repHome(const char *file);

   inline int INT(const double a) { 
      int i=::round(a);
      if (std::fabs(i-a)>1E-10) wblog(FL,"ERR %s(%g) !?",FCT,a);
      return i;
   }

   const char* strstri (const char *s1, const char *s2);
   const char* strstrw (const char *s1, const char *s2);
   const char* strstrwi(const char *s1, const char *s2);
   const char* strchri (const char *s, char c);

   const char* rstrstr (const char *s1, const char *s2);

   void shift(char *s, unsigned n, int p, char erase=0);

   void shift(char *s, unsigned n, int p, const char *srep);

   bool isword(const char &c) { return (
        (c>='0' && c<='9') ||
        (c>='a' && c<='z') ||
        (c>='A' && c<='Z') || c=='_'
   ); };

   template <class T>
   double mod(double a, T b) {
       if (b<=0) wblog(FL,"ERR %s() got b=%g",FCT,b);
       a=::fmod(a,double(b)); if (a<0) { a+=b; }
       return a;
   };

   template <class T, ENABLE_IF_isINT(T)>
   int mod(int a, T b) {
       if (b<=0) wblog(FL,"ERR %s() got b=%d",FCT,b);
       a %= int(b); if (a<0) { a+=b; }
       return a;
   };

   unsigned strnlen(const char *s, unsigned n) {
      for (unsigned i=0; i<n; ++i) if (!s[i]) return i;
      return n;
   };

   template <class T>
   bool isLower(const T* a, const T* b, const size_t n);

   template <class T>
   bool isEqual(const T* a, const T* b, const size_t n);

   template <class T>
   bool allEqual(const T* a, const size_t n, const T &x);

   template <class T>
   bool anyEqual(const T* a, const size_t n, const T &x);

   template <class T>
   bool anyUnqual(const T* a, const size_t n, const T &x);

   template <class T>
   T getdscale(const T* d, size_t n);

   template <class T> inline
   void scale_eps(T &eps, const T* d, size_t n);

} 

template <class T, class T2> inline
char GOT_EPS(const T &a, const T2 &eps) {
   if (a==0 || (eps!=0 && ABS(a)<eps)) return 1;
   return 0;
}

void setSigHandler(char i); 
void wbSigHandler(int sig);

char *gsh_F=NULL; int gsh_L=0;

namespace Wb {
class SigHandler {

  public:
    SigHandler(const char *F, int L) {
       if (F && F[0]) {
          const char *s=strrchr(F,'/'); if (!s) { s=F; }
          gsh_F = new char[strlen(s)+1]; strcpy(gsh_F,s);
          gsh_L=L;
       }
       else {
          const char s0[]=__FILE__, *s=strrchr(s0,'/'); if (!s) { s=s0; }
          gsh_F = new char[strlen(s)+1]; strcpy(gsh_F,s);
          gsh_L=__LINE__;
       }
    #ifdef DBSTOP
    #else
       setSigHandler('i');
    #endif
    };

   ~SigHandler() {
    #ifdef DBSTOP
    #else
       setSigHandler('r');
    #endif
       gsh_L=0; if (gsh_F) { delete [] gsh_F; gsh_F=NULL; }
    };

    int check911(char tflag=0) { doflush();
       if (tflag<=0) {
          #ifndef DBSTOP
             wbSigHandler(911);
          #endif
          return 0;
       }
       else { return check911_(); }
    };

    static int check911_(const char *F=NULL, int L=0) {
       int i=SigHandler::icount,
           a=SigHandler::acount, n=SigHandler::NC0;

       if (i>=n) { i=0; } else { i=n-i; if (i>9) { i=9; }}
       if (a>=n) { a=0; } else { a=n-a; if (a>9) { a=9; }}
       if (F) { wblog(F_L,"ERR %s() ",FCT); }
       return (i + 10*a);
    };

    static unsigned NC0, icount, acount;

  protected:
  private:

};

   unsigned SigHandler::NC0=3;

   unsigned SigHandler::icount=SigHandler::NC0;
   unsigned SigHandler::acount=SigHandler::NC0;

}; 

template <class TX>
class wbREF {

  public:

    wbREF (const TX &A, char flag) : orig(&A), ref(0) {
       if (flag) ref=new TX(A); else ref=&A;
    };

    TX& Ref() { return *ref; }

   ~wbREF() { if (ref!=orig) delete ref; }

  protected:
  private:
    const TX *orig, *ref;
};

namespace Wb {

template <class T>
class num2Fmt { 

  public:
   num2Fmt(const char *F, int L,
      int m_=-2, int p_=-2, const char *t_="") : m(m_), p(p_) {

      if (t_ && t_[0]) {
         if (strlen(t_)>2) wblog(F_L,"ERR %s() invalid t='%s'",FCT,t_);
         strcpy(t,t_);
      }
      else { memset(t,0,3); }

      check_init();
   };

   int check_init();

   explicit operator wbstring() const;

  private: 
   char m,p;
   char t[3]; 

   void get_fmt(wbstring &s) const;
};

   template<class T> inline
   int num2int(const char *F, int L, const T &x);

   template<class T> inline
   unsigned checkInt(
      const char *F, int L, const T* x, size_t n,
      double eps=1E-14
   );

   char* strpad(char *s, char c, unsigned n, unsigned w=1) {
      unsigned i,l=strlen(s);
      for (i=0; i<w; ++i) s[l++]=' '; 
      for (i=l; i<n; ++i) s[l++]=c;
      return s;
   };

   unsigned strpad_(char *s, unsigned n, char c=' ', int z=1) {
      if (int(n)>0) {
         unsigned i=0; for (; i<n; ++i) { s[i]=c; }
         if (z) s[i]=0;
         return i;
      }
      return 0;
   };

   template <class T>
   int charGetNumber(const char *F, int L, const char *s, T &x);

   bool isBaseType(const std::type_info &type_id);
   bool isComplexType(const std::type_info &type_id);

   wbstring TimeStamp(char type=0);

   void inl(unsigned n) { for (unsigned i=0; i<n; i++) printf("\n"); }

   void quietErrLog(const char *F, int L, const char *msg, char qflag) {
       if (qflag) snprintf(str,256,"%s [%s:%d]",msg,Wb::basename(F),L);
       else wblog(F,L,"ERR %s",msg);
   }

   template <class T1, class T2> inline
   void safeConvert(const char *F, int L, const T1 &x1, T2 &x2);

   template <class T1, class T2> inline
   int checkUpdate(T1& a, const T2 &r) {
       if (a!=T1(r)) { a=T1(r); return 1; } else return 0;
   };

   template <class T>
   bool uniformRange(const T* d, const size_t n);

   template <class T>
   T maxRange(const T* d, const size_t n);

   template <class T>
   T minRange(const T* d, const size_t n);

   template <class T>
   T addRange2(const T* d, size_t n, const size_t stride=1); 

   template <class T>
   T addRange(const T* d, size_t n, const size_t stride=1); 

   template <class T>
   void addRange( 
      const T* a, const T* b, T* c, size_t n, 
      const size_t stride=1
   );

   template <class T>
   void addRange(const T* a, T* c, size_t n, const size_t stride=1); 

   template <class T>
   void minusRange( 
      const T* a, const T* b, T* c, const size_t n, 
      size_t stride=1
   );

   template <class T>
   void diffRange(
      const T* a, const T* b, T* c, size_t n, 
      size_t stride=1
   );

   template <class T> 
   void setRange2avg(T* a, size_t m);

   template <class T> 
   size_t nnzRange(const T* a, size_t m);

   template <class T, class T2> inline
   void cpyStride(
      T2* d,          
      const T* d0,    
      unsigned n,     
      const unsigned *idx, unsigned m,
      unsigned D2=-1, 
      unsigned D0=-1, 
      char add_flag=0
   );

   template <class T, class T2> inline
   void cpyStride(
      T2* d,          
      const T* d0,    
      unsigned n,     
      unsigned m,     
      unsigned D2=-1, 
      unsigned D0=-1, 
      char add_flag=0
   ){
      return cpyStride(d,d0,n,(unsigned*)NULL,m,D2,D0,add_flag);
   };

   template <class T>        
   void cpyRange(T* a, const T* b, size_t n); 

   template <class T, class Tb>
   void cpyRange(T* a, const Tb* b, size_t n, char check_type=1);

   template <class T>
   void cpyRangeR(T* a, const wbcomplex* z, size_t n);
   template <class T>
   void cpyRangeI(T* a, const wbcomplex* z, size_t n);
   template <class T>
   void cpyRangeA(T* a, const wbcomplex* z, size_t n);

   template <class T>
   T prodRange(const T* d, const size_t n);
   template <class T>
   T prodRange(const T* d, const size_t n, T x0);

   template <class T> 
   T dotProd(const T* a, const T* b, size_t n, char conj=0);
   template <class T> 
   T dotProd(const T* a, const T* b, size_t n, T x, char conj=0);

   template <class T> inline
   void TimesElRange(T* a, const T* b, size_t n, char conj=0);

   template <class T> inline
   void timesElRange(T* c, const T* a, const T* b, size_t n, char conj=0);

   template <class T> inline
   char check_conj_flag(const T* a, char &conj);

   template <class T> inline
   void timesElRange_OM( 
      T* c, 
      const T* a, const T* b, 
      size_t n, unsigned m,     
      char conj=0);

   template <class T>
   void timesElRange_add(
      T* c, const T* a, const T* b, size_t n, T bfac=T(1), char conj=0);

   template <class TD, class TX> inline
   void timesRange(TD* d, TX x, const size_t n, size_t stride=1);

   template <class T> 
   T rangeNorm2(const T* d, size_t n, size_t *k=NULL);

   template <class T> 
   T rangeNormDiff2(const T* d1, const T* d2, size_t n,
      const T &fac, size_t *k=NULL
   );

   template <class T> 
   T rangeNormDiff2(const T* d1, const T* d2, size_t n);

   template <class T>
   T rangeMaxDiff(const T* d1, const T* d2, size_t n, size_t *k=NULL);

   template<class T>
   T overlap(
      const T *v1, const T *v2, size_t n,
      size_t stride=1, char tnorm=0
   );

   template<class T>
   T gs_project_range(
      T *v1, const T *v2, 
      size_t n, size_t stride=1, char isnorm=0, char tnorm=0
   );

   template<class T>
   void chopTiny_float(T *d, size_t n, T dref=-1) { return; };

   template<class T>  
   void chopTiny_z(T *d, size_t n, double eps=1E-14) { return; };

   template<class T> 
   double chopTiny_imag(T *d, size_t n, double eps=1E-14) { return 0; };

   template <class T> inline
   int recCompare(
      const T* a, const T* b, const size_t n, char lex=1
   ){
      if (lex>0) {
         for (size_t i=0; i<n; ++i) {
            if (a[i]<b[i]) return -1; else
            if (a[i]>b[i]) return +1;
         }
      }
      else {
         for (size_t i=n-1; i<n; --i) {
            if (a[i]<b[i]) return -1; else
            if (a[i]>b[i]) return +1;
         }
      }
      return 0;
   };

   template <class T> inline
   int recCompare( 
      const T* a, const T* b, const size_t n, char lex,
      T eps, int *nwrn_=NULL
   ){
      if (eps<=T(0)) return recCompare(a,b,n,lex);

      int nwrn0=0, *nwrn=(nwrn_ ? nwrn_ : &nwrn0);
      T x=0.2*eps;

      if (lex>0) {
         for (size_t i=0; i<n; ++i) {
            if (a[i]<b[i]-eps) return -1; else
            if (a[i]>b[i]+eps) return +1; else
            if (!nwrn[0] && std::fabs(a[i]-b[i])>x) { 
               wblog(FL,"WRN %s() got |a-b|/eps=%.3g (%d/%d)",
               FCT,(a[i]-b[i])/eps,i+1,n); ++nwrn[0];
            }
         }
      }
      else {
         for (size_t i=n-1; i<n; --i) {
            if (a[i]<b[i]-eps) return -1; else
            if (a[i]>b[i]+eps) return +1; else
            if (!nwrn[0] && std::fabs(a[i]-b[i])>x) { 
               wblog(FL,"WRN %s() got |a-b|/eps=%.3g (%d/%d)",
               FCT,(a[i]-b[i])/eps,i+1,n); ++nwrn[0];
            }
         }
      }
      return 0;
   };

   template <class T>
   int cmpRange( 
      const T *a, const T *b, const size_t n, char lex=+1
   ){
      if (lex>0) {
         for (size_t i=0; i<n; ++i) { if (a[i]!=b[i]) {
            return (a[i]<b[i] ? -1 : +1);
         }}
         return 0;
      }
      else {
         if (!lex) wblog(FL,"WRN %s() got lex=%d",FCT,lex);
         for (size_t i=n-1; i<n; --i) { if (a[i]!=b[i]) {
            return (a[i]<b[i] ? -1 : +1);
         }}
         return 0;
      }
   };

   template <class T> void cpyZRange(
   const double *R, const double *I, T *Z, const size_t n);

   template <class T> void splitZRange(
   const T *Z, double *R, double *I, const size_t n);

   template<class T>
   void invertIndex(
      const wbvector<T> &i1, unsigned N,
      wbvector<T> &i2, char lflag=1
   );

   void setRand(wbvector<wbcomplex> &zz, double fac=1., double shift=0.);

   template<class T>
   void set2avg(T *dd, const size_t n);

   template<class T>
   void getDiff(const wbvector<T> &xx, wbvector<T> &dx);

   template<class TX, class TY>
   TY IntTrapez(const TX *xx, const TY *yy, const size_t n);

   inline bool contains(const char* s, const char &x, const unsigned &m=16) {
      if (s)
      for (unsigned i=0; s[i] && i<m; i++) { if (s[i]==x) return 1; };

      return 0;
   };

}; 

   template<class T>
   void markSet(
      wbvector< wbvector<T>* > &E0,
      wbvector< wbvector<char> > &mark,
      unsigned Nkmin,     
      int &Nkeep,         
      double Etrunc,      
      wbvector<T> &E,     
      double eps=0.,      
      double b=0.,        
      int dmax=-1,        
      const char *dir=0   
   );

   template<class T>
   void markSet(const char *F, int L,
      const wbvector<T> &E,
      wbvector<char> &mark,
      unsigned Nkmin,     
      int &Nkeep,         
      double Etrunc=0,    
      double b=0.,        
      int dmax=-1,        
      const char *dir=0   
   );

#endif

