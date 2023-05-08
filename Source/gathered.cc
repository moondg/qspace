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

#ifndef __WB_GATHERED_CC__
#define __WB_GATHERED_CC__

// auxilliary routine for mexlib.* -> Mx::Array
// to deal with type conversions in a controlled way // Wb,May11,19

void setSigHandler(char w) {

   static sighandler_t prev_SIGINT;
   static sighandler_t prev_SIGABRT;

   static int iflag=0; 

   if (w=='i') {    
      if (++iflag==1) { 
         prev_SIGINT  = signal(SIGINT,  wbSigHandler);
         prev_SIGABRT = signal(SIGABRT, wbSigHandler);
      }
   }
   else if (w=='r') { 
      if (--iflag==0) {
         signal(SIGINT,  prev_SIGINT );
         signal(SIGABRT, prev_SIGABRT);
      }
      else if (iflag<0) { wblog(gsh_F, gsh_L,
        "ERR SigHandler not yet initialized (iflag=%d)",iflag);
      }
   }
   else wblog(gsh_F, gsh_L,"ERR invalid flag=%s",cSTR(w));
};

void wbSigHandler(int sig) {

   unsigned i=Wb::SigHandler::icount,
            a=Wb::SigHandler::acount, n=Wb::SigHandler::NC0;
   fflush(0);

   if (sig==911) {
      if (i!=n || a!=n) {
         Wb::SigHandler::icount=Wb::SigHandler::acount=n;
         wblog(gsh_F, gsh_L,
           "ERR terminating due to %s (%d/%d) --------------------",
            i<a ? "SIGINT":"SIGABRT", n-MIN(i,a), n
         );
      }
      return;
   }

   unsigned l=128; char istr[l], s[16]; strcpy(s,"(signal)");

   if (sig==SIGINT) { 
      strcpy(s,"SIGINT"); --Wb::SigHandler::icount;
   }
   else if (sig==SIGABRT) { 
      strcpy(s,"SIGABRT"); --Wb::SigHandler::acount;
      i=a; 
   }
   else { printf("\n");
      snprintf(istr,l,"%s ERR terminating due to uncaught signal %d",
         shortFL(gsh_F,gsh_L), sig);
      ExitMsg(istr); fflush(0); exit(1); 
   }

   if (i<=0) { 
      snprintf(istr,l,"%s ERR terminating due to %s (%d/%d) %s",
         shortFL(gsh_F,gsh_L),s,n,n,"--------------------");
      ExitMsg(istr);
      fflush(0); exit(-1); 
   }

   if (i>=3) { 
      snprintf(istr,l,"--- received %s (%d) %s",
      s,i,"--------------------------------"); }
   else if (i==2) {
      snprintf(istr,l,"--- received %s (%s)",
      s,"waiting for natural termination pt"); }
   else {
      snprintf(istr,l,"WRN next %s will result in immediate stop %s",
      s,"----------");
   }

#ifdef __APPLE__
   time_t curr_time=time(NULL);
   struct tm *tblock=localtime(&curr_time);
   char t[16]; strftime(t,16,"%T",tblock); 

   fprintf(stdout,"%-20s %s  %s\n",shortFL(gsh_F,gsh_L),t,istr);
#else
   wblog(gsh_F,gsh_L,istr);
#endif
};

void dbstop(const char* F, int L) {

   unsigned ith=omp_get_thread_num(), nth=omp_get_num_threads();
   char sx[8]; sx[0]=0;
   if (nth>1) { snprintf(sx,8," (%d/%d)",ith,nth); }

   fprintf(stdout,"%s %s %s() raising SIGTRAP%s ...\n\n",
      shortFL(F_L),TIME,FCT,sx);
   fflush(0);

   raise(SIGTRAP);
};

template<class T> inline
int num2int(const char *F, int L, const T &x){
   double d=x; int i=round(d);
   if (F) {
      double e=fabs(d-i); if (e>1E-12) {
      wblog(F,L,"ERR %s() got non-int value %g !?",FCT,d); }
   }
   return i;
};

template<> inline
int num2int(const char *F __attribute__ ((unused)),
   int L __attribute__ ((unused)), const char &x){ return x; }

template<> inline
int num2int(const char *F __attribute__ ((unused)),
   int L __attribute__ ((unused)), const unsigned char &x){ return x; }

template<> inline
int num2int( const char *F __attribute__ ((unused)),
   int L  __attribute__ ((unused)), const int &x){ return x; }

template<> inline
int num2int( const char *F __attribute__ ((unused)),
   int L  __attribute__ ((unused)), const long &x){ return x; }

template<> inline
int num2int(const char *F __attribute__ ((unused)),
   int L __attribute__ ((unused)), const unsigned &x){ return x; }

template<> inline
int num2int(const char *F __attribute__ ((unused)),
   int L __attribute__ ((unused)), const unsigned long &x){ return x; }

template<class T> inline
unsigned checkInt(const char *F, int L,
   const T* x, size_t n, double eps
){
    for (size_t i=0; i<n; ++i) {
        if (fabs(x[i]-round(x[i]))>eps) {
           if (F) wblog(F_L,
              "ERR %s() got data[%d]=%g (having %s) !?",FCT,
               i+1, double(x[i]), TSTR(T));
           return (i+1);
        }
    }
    return 0;
};

template<> inline
unsigned checkInt(
   const char *F     __attribute__ ((unused)),
   int L             __attribute__ ((unused)),
   const int* x      __attribute__ ((unused)),
   size_t n          __attribute__ ((unused)),
   double eps        __attribute__ ((unused))
){ return 0; }

template<> inline
unsigned checkInt(
   const char *F     __attribute__ ((unused)),
   int L             __attribute__ ((unused)),
   const unsigned* x __attribute__ ((unused)),
   size_t n          __attribute__ ((unused)),
   double eps        __attribute__ ((unused))
){ return 0; }

template<> inline
unsigned checkInt(
   const char *F     __attribute__ ((unused)),
   int L             __attribute__ ((unused)),
   const char* x     __attribute__ ((unused)),
   size_t n          __attribute__ ((unused)),
   double eps        __attribute__ ((unused))
){ return 0; }

const char* Wb::basename(const char *data, char c) {

    int i, n=strlen(data);                
    for (i=n-1; i>=0; i--) if (data[i]==c) { i++; break; }

    return (data + (i>0 ? i : 0));
};

inline const char* Wb::strchri(const char *s, char c) {

   if (!s) wblog(FL,"ERR %s() got invalid null string",FCT);

   for (c=tolower(c); *s; ++s) {
      if (tolower(*s)==c) return s;
   }

   return NULL;
};

inline const char* Wb::rstrstr(const char *s1, const char *s2) {
   if (!s1 || !s2) wblog(FL,"ERR %s() got null strings (%p/%p)",FCT,s1,s2);
   if (!*s2) { return s1; } 
   if ( *s1) { unsigned n=strlen(s1), m=strlen(s2);
      if (n>=m) { unsigned i,j,il=n-1;
         for (; il<n; --il) {
            for (i=il, j=m-1; i<n && j<m && s1[i]==s2[j]; --i, --j) {};
            if (j>m) { ++i; return s1+i; }
         }
      }
   }
   return NULL;
};

inline const char* Wb::strstri(const char *s1, const char *s2) {

   if (!s1 || !s2) wblog(FL,
      "ERR %s() got invalid null strings (%p, %p)",FCT,s1,s2);

   const char *a, *b;

   for(; *s1; ++s1) {
      for(a=s1, b=s2; *a && tolower(*a)==tolower(*b); ++a, ++b) {};
      if(!*b) return s1;
   }

   return NULL;
};

inline const char* Wb::strstrw(const char *s1, const char *s2) {

   if (!s1 || !s2) wblog(FL,
      "ERR %s() got invalid null strings (%p, %p)",FCT,s1,s2);

   unsigned i=0, j=0, k=0;

   for(; s1[i]; ++i) {
      for (j=i, k=0; s1[j] && s1[j]==s2[k]; ++j, ++k) {};
      if (!s2[k]) {
         if ((!i || !Wb::isword(s1[i-1])) && !Wb::isword(s1[j]))
         return s1+i;
      } 
   }

   return NULL;
};

inline const char* Wb::strstrwi(const char *s1, const char *s2) {

   if (!s1 || !s2) wblog(FL,
      "ERR %s() got invalid null strings (%p, %p)",FCT,s1,s2);

   unsigned i=0, j=0, k=0;

   for(; s1[i]; ++i) {
      for (j=i, k=0; s1[j] && tolower(s1[j])==tolower(s2[k]); ++j, ++k);
      if (!s2[k]) {
         if ((!i || !Wb::isword(s1[i-1])) && !Wb::isword(s1[j]))
         return s1+i;
      } 
   }

   return NULL;
};

void Wb::shift(char *s, unsigned n, int p, char erase) {

   if (unsigned(abs(p))>=n) {
      if (erase && n) { unsigned i=0;
         if (p>0) { for (; i<n; ++i) { s[i]=' '; }}
         else {
            for (; i<n && s[i]; ++i) { s[i]=' '; }
            if (i<n) { s[0]=0; }
         }
      }
   }
   else if (p>0) {
      unsigned i=n-p-1; for (; i<n; --i) { s[i+p]=s[i]; }
      if (erase) { for (i=0; int(i)<p; ++i) { s[i]=' '; }}
   }
   else if (p<0) {
      unsigned i=-p; for (; i<n && s[i]; ++i) { s[i+p]=s[i]; } 
      if (erase) {
         if (n>i) { n=i; } 
         for (i=i+p+1; i<n && s[i]; ++i) { s[i]=' '; }
      }
   }
};

void Wb::shift(char *s, unsigned n, int p, const char *s_) {

   if (!p) return;
   if (!s || !s_) wblog(FL,"ERR %s() invalid input %p / %p",FCT,s,s_);

   if (unsigned(abs(p))>=n) { if (n) {
      if (p>0) { strncpy(s,s_,n); }
      else {
         unsigned i=0; for (; i<n && s[i]; ++i) { }
         if (i<n) { s[0]=0; }
         else { strncpy(s,s_,n); }
      }
   }}
   else if (p>0) {
      unsigned i=n-p-1; for (; i<n; --i) { s[i+p]=s[i]; }
      for (i=0; int(i)<p ; ++i) { s[i]=s_[i]; if (!s[i]) break; }
   }
   else if (p<0) { unsigned i=-p;
      for (; i<n && s[i]; ++i) { s[i+p]=s[i]; } 
      for (i=n+p; i<n; ++i) { s[i]=s_[i]; if (!s[i]) break; }
   }
};

template <class T>
int Wb::charGetNumber(const char *F, int L, const char *s, T &x) {

    float fl=0; 
    int i;

    if (!s) {
       sprintf(str,"cannot read number from empty string");
       return 1;
    }

    i=sscanf (s,"%f",&fl);
    x=(T)fl;

    if (!i || fl!=(float)x) {
       sprintf(str,"invalid input for type `%s' (%g, %d)",
       TSTR(T), fl, i); return 1;
    }

    return 0;
}

template<>
int Wb::GetEnv (const char *F, int L, const char *name, double &val) {

    char *s, *s2;
    double dbl;

    s=getenv(name);    if (!s || !s[0]) return -1;
    dbl=strtod(s,&s2); if (!s2) return -2;
    val=dbl;

    if (F) wblog1(F,L," *  getenv %-6s = %g", name, val);

    if (*s2) {
       wblog1(FL,"WRN %s() got invalid %s='%s' !?",FCT,name,s);
       if (isspace(*s2)) {
          ++s2; while (*s2 && isspace(*s2)) ++s2;
       }
    }

    return (*s2 ? -3 : 0);
};

template<>
int Wb::GetEnv (const char *F, int L, const char *name, int &val) {

    double dbl=0;

    int i=GetEnv(F,L,name,dbl); if (i) return i;

    i=int(dbl); if (double(i)!=dbl) return -2;
    val=i;

    return 0;
}

template<>
int Wb::GetEnv (const char *F, int L, const char *name, unsigned &val) {

    double dbl=0;
    int i=GetEnv(F,L,name,dbl); if (i) return i;

    i=unsigned(dbl); if (double(i)!=dbl) return -2;
    val=i;

    return 0;
};

template<>
int Wb::GetEnv (const char *F, int L, const char *name, char &val) {

    char *s=getenv(name); if (!s || !s[0]) return -1;
    if (s[1]) { return -2; }
    val=s[0];

    if (F) wblog1(F,L," *  getenv %-6s = '%c'(%d)", name, val,val);
    return 0;
};

template<>
int Wb::GetEnv (const char *F, int L, const char *name, wbstring &val) {

   if (!name || !name[0]) wblog1(F_L,
      "ERR %s() invalid env '%s' !?",FCT,name);

   char *s=getenv(name);

   if (!s) {
      val.init(); return -1;
   }
   else {
      val=s;
      if (F) wblog1(F,L," *  getenv %-6s = '%s'",name,val.data);
      return 0;
   }
};

int Wb::EnvIsSet(const char *F, int L, const char *name) {
   int x=0, e=GetEnv(F,L,name,x);
   return (!e && x>0 ? 1 : 0);
};

int Wb::get_WB_VERBOSE(const char *F, int L) {

   const char *vname="WB_VERBOSE";
   int i, k=0xF1; 

   if ((i=GetEnv(0,0,vname,k))==0) {
      if (k<0) wblog1(F_L, 
         "WRN %s() invalid %s=%d/255",FCT,vname,k);
      else if (k!=Wb::envVRB) {
         if ((F && k) || k&8) { unsigned n=32; char s[n];
            i=snprintf(s,n,"%s()",myname);
            if (i<0) wblog(FL,"WRN %s() i=%d",FCT,i);
            wblog1(PF_L,"ENV %-14s %s = %X->%X",s,vname,Wb::envVRB,k); 
         }
         Wb::envVRB=k;
      }
   }
   else if (int(i)!=-1) {
      wblog1(F_L,"WRN %s() invalid %s (e=%d)",FCT,vname,i);
   }

   { int q=0;
     i=GetEnv(0,0,"WB_LOG_COLOR",q);
     if (i) { q=(Wb::envDKT>0 ? WLC_DARK : WLC_OFF); } else
     if (q>=NUM_WBL_COLOR_SCHEME) { q=WLC_DARK; } else
     if (q<0) { q=WLC_OFF; }
     Wb::useCol=(WBL_COLOR_SCHEME)q;
   }

   return k;
};

int Wb::got_DBSTOP(const char *F, int L) {

   int i=0,x=0,e,q;
   const char *ss[]={"DBSTOP","ML_DEBUG","DEBUG"};

   for (; i<3; ++i) { q=0;
      e=GetEnv(0,0,ss[i],q);
      if (!e) { if (q) { x|=(1<<i); }} else
      if (e!=-1 && F) wblog1(F_L,
         "WRN %s() invalid value for ENV %s (e=%d)",FCT,ss[i],e);
   }

   if (x!=Wb::envDBG) {
      if ((F && x && Wb::envVRB&15) ||
          (Wb::envDBG && Wb::envVRB&8)) wblog1(PF_L," *  "
         "using %-20s = %d -> %d","Wb::envDBG",Wb::envDBG,x);
      Wb::envDBG=x;             
   }

   return x;
};

#ifdef MATLAB_MEX_FILE
int Wb::got_DESKTOP() { 

   int q=(
       Wb::CallMatlab(0,0,"isdeployed") ? 0 :
       Wb::CallMatlab(0,0,"usejava","desktop") 
   );

   if ((Wb::envDKT>=0 && Wb::envDKT!=q) || (Wb::envVRB&8)) wblog1(FL,
      "ENV %s() Wb::envDKT = %d (%s)",myname,q,q?"desktop":"terminal");
   Wb::envDKT=q;

   return q;
};

int Wb::is_DEPLOYED(const char *F, int L) { 
   return (int)Wb::CallMatlab(F,L,"isdeployed");
};

#else
int Wb::got_DESKTOP() { return 0; } 
int Wb::is_DEPLOYED(
   const char *F=0 __attribute__ ((unused),
   int L=0 __attribute__ ((unused)) { return 1; } 
#endif

int Wb::GetNumThreads(const char *F, int L, int &n, const char *name) {

#ifndef QS_USING_OMP
    return 0; 
#endif

    double q=0;
    int i=GetEnv(0,0,name,q), ncpu=wbsys::getNumCores();

    if (i<0) {
       if (i<-1 && F) wblog1(F_L,
          "ERR %s() invalid value for env %s (i=%d)",FCT,name,i);
       return q;
    }
    if (q<0 || double(int(q))!=q) wblog1(F_L,
       "ERR %s() invalid value for env %s (%g)",FCT,name,q);
    if (q>ncpu) wblog1(F_L,"ERR %s() "
       "env %s=%g exceeds number of cores (%d)",FCT,name,q,ncpu);

    if (n!=q) { i=0;

       if (!GetEnv(0,0,"ML_DEBUG",i) && i<1) { q=0; } else
       if (!GetEnv(0,0,"DEBUG",   i) && i<1) { q=0; } else
       if (F && ((Wb::envVRB &12) 
          #ifdef LOAD_CGC_QSPACE
           || CG_VERBOSE>6
          #endif        
       )) wblog1(PF_L," *  using %-20s = %d -> %g",name,n,q);
       n=q;
    }

    return q;
};

int Wb::strrep(
   const char *S0, 
   const char *t0, 
   const char *t2, 
   char *Sout,     
   size_t N,       
   char gflag      
){
   if (!S0 || !S0[0] || !t0 || !t0[0]) return 0;

   int nrep=0;
   size_t k, l=0, m0=strlen(t0), m=(t2 ? strlen(t2):0);
   const char *s, *s0=S0; char *s2=Sout;

   while ((s=strstr(s0,t0))) { k=s-s0; ++nrep;
      if (k) {
         if (l+k>N) wblog(FL,
           "ERR string out of bounds (%d+%d / %d)",l,k,N);
         memcpy(s2,s0,k); l+=k; s2+=k;
      }
      if (m) {
         if (l+m>N) wblog(FL,
           "ERR string out of bounds (%d+%d / %d)",l,m,N);
         memcpy(s2,t2,m); l+=m; s2+=m;
      }
      s0+=(k+m0); if (!gflag) break;
   }

   if ((k=strlen(s0))) { 
      if (l+k>N) wblog(FL,
        "ERR string out of bounds (%d+%d / %d)",l,k,N);
      memcpy(s2,s0,k); l+=k; s2+=k;
   }
   s2[0]=0; 

   return nrep;
};

wbstring Wb::repHome(const char *file) {

   size_t n=2*strlen(file);
   const char *s; char F0[n+1], F[n+1];
   char *f0=F0, *f=F; strcpy(f,file);

   if ((s=getenv("MEX" ))) {
      SWAP(f,f0); if (Wb::strrep(f0,s,"$MEX", f,n)) return f; }
   if ((s=getenv("HOME"))) {
      SWAP(f,f0); if (Wb::strrep(f0,s,"$HOME",f,n)) return f; }
   if ((s=getenv("LMA" ))) {
      SWAP(f,f0); if (Wb::strrep(f0,s,"$LMA", f,n)) return f; }
   if ((s=getenv("PROJ"))) {
      SWAP(f,f0); if (Wb::strrep(f0,s,"$PROJ",f,n)) return f; }
   if ((s=getenv("PROJ_LEO"))) { 
      SWAP(f,f0); if (Wb::strrep(f0,s,"$PROJ_LEO",f,n)) return f; }
   if ((s=getenv("USER"))) {
      SWAP(f,f0); if (Wb::strrep(f0,s,"$USER",f,n)) return f; }

   return f;
};

void Wb::print_backtrace(const char *F, int L, const char *istr) {

   unsigned i=0, n=4;
   GetEnv(0,0,"QS_NUM_BTRACE",i); if (n<i) { n=i; }

   void *array[n];

   n=backtrace(array,n); 

   char **ss=backtrace_symbols(array,n); 

   printf("\n# %s: stack backtrace",shortFL(F_L)); 
      if (istr && istr[0]) { printf(" %s",istr); }
      printf(" (%d frame%s)\n",n, n!=1? "s":"");
   printf("# hint: use addr2line to decode this\n\n");

   if (ss) {
      for (i=0; i<n; ++i) {
         printf("  [%*d]  %s\n", n>10? 2:1, i,ss[i]);
      }; printf("\n");
      free(ss); 
   }
};

std::string getName(
   const std::type_info &t, 
   char v 
){

   std::string s; 

   if (t==typeid(float        )) s=(v? "float"       : "f"  ); else
   if (t==typeid(double       )) s=(v? "double"      : "d"  ); else
   if (t==typeid(bool         )) s=(v? "bool"        : "b"  ); else
   if (t==typeid(char         )) s=(v? "char"        : "c"  ); else
   if (t==typeid(unsigned char)) s=(v? "uchar"       : "uc" ); else
   if (t==typeid(int          )) s=(v? "int"         : "i"  ); else
   if (t==typeid(unsigned     )) s=(v? "unsigned"    : "u"  ); else
   if (t==typeid(long         )) s=(v? "long"        : "l"  ); else
   if (t==typeid(size_t       )) s=(v? "size_t"      : "st" ); else
   if (t==typeid(unsigned long long)) s=(v? "uLong"  : "uL" ); else 

   if (t==typeid(double*      )) s=(v? "double*"     : "d*" ); else
   if (t==typeid(float*       )) s=(v? "float*"      : "f*" ); else
   if (t==typeid(bool*        )) s=(v? "bool*"       : "b*" ); else
   if (t==typeid(char*        )) s=(v? "char*"       : "c*" ); else
   if (t==typeid(int*         )) s=(v? "int*"        : "i*" ); else
   if (t==typeid(unsigned*    )) s=(v? "unsigned*"   : "u*" ); else
   if (t==typeid(long*        )) s=(v? "long*"       : "l*" ); else
   if (t==typeid(size_t*      )) s=(v? "size_t*"     : "st*"); else
   if (t==typeid(unsigned long long*))s=(v? "uLong*" : "uL*"); else 

   if (t==typeid(const double*)) s=(v? "c-double*"   : "Kd*"); else
   if (t==typeid(const float* )) s=(v? "c-float*"    : "Kf*"); else
   if (t==typeid(const bool*  )) s=(v? "c-bool*"     : "Kb*"); else
   if (t==typeid(const char*  )) s=(v? "c-char*"     : "Kc*"); else
   if (t==typeid(const int*   )) s=(v? "c-int*"      : "Ki*"); else
   if (t==typeid(const unsigned*)) s=(v? "c-unsigned*":"Ku*"); else
   if (t==typeid(const long*  )) s=(v? "c-long*"     : "Kl*"); else
   if (t==typeid(const size_t*)) s=(v? "c-size_t*"   : "Kst*"); else

   if (t==typeid(wbcomplex    )) s=(v? "wbcomplex"   : "z"  ); else
   if (t==typeid(wbcomplex*   )) s=(v? "wbcomplex*"  : "z*" ); else

   if (t==typeid(long double  )) s=(v? "long double" : "dl" ); else
   if (t==typeid(long double* )) s=(v? "long double*": "dl*"); else

#ifdef __WB_MPFR_HH__
   if (t==typeid(Wb::quad     )) s=(v? "Wb::quad"    : "q"  ); else
   if (t==typeid(Wb::qquad    )) s=(v? "Wb::qquad"   : "qq" ); else
   if (t==typeid(Wb::quad2    )) s=(v? "Wb::quad2"   : "q2" ); else
   if (t==typeid(Wb::quad3    )) s=(v? "Wb::quad3"   : "q2" ); else

   if (t==typeid(Wb::quad*    )) s=(v? "Wb::quad"    : "q*" ); else
   if (t==typeid(Wb::qquad*   )) s=(v? "Wb::qquad"   : "qq*"); else
   if (t==typeid(Wb::quad2*   )) s=(v? "Wb::quad2"   : "q2*"); else
   if (t==typeid(Wb::quad3*   )) s=(v? "Wb::quad3"   : "q2*"); else
#endif

   s=(char*)(t.name());

   return s;
};

bool Wb::isBaseType(const std::type_info &t) { 

   if ( t==typeid(unsigned)
     || t==typeid(double)
     || t==typeid(float)
     || t==typeid(int)
     || t==typeid(char)
     || t==typeid(unsigned char)
     || t==typeid(long long)
     || t==typeid(size_t)
     || t==typeid(wbcomplex)
   ) return 1;

   return 0;
};

bool isComplexType(const std::type_info &t) {

   if ( t==typeid(wbcomplex) ) return 1;

   return 0;
};

template<class T>
Wb::num2Fmt<T>::operator wbstring() const {
   wbstring s; 
   get_fmt(s); return s;
};

template<class T>
void Wb::num2Fmt<T>::get_fmt(wbstring &s) const {

   unsigned l=0; s.init(12); 

   if (m<-1 || p<-1 || t[0]<=0) wblog(FL,
      "ERR %s<%s> invalid type %d.%d%s",FCT,TSTR(T),m,p,t);

   if (m>0) {
      if (p>=0)
           { l=snprintf(s.data,s.len,"%%%d.%d%s",m,p,t); }
      else { l=snprintf(s.data,s.len,"%%%d%s",m,t); }
   }
   else {
     if (p>=0)
          { l=snprintf(s.data,s.len,"%%.%d%s",p,t); }
     else { l=snprintf(s.data,s.len,"%%%s",t); }
   }

   if (l>=s.len) wblog(FL,
      "ERR %s<%s> string out of bounds (%%%d.%d%s; %d/%d)",
      FCT,TSTR(T),m,p,t,l,s.len
   );
};

template<>
void Wb::num2Fmt<wbcomplex>::get_fmt(wbstring &s) const {

   unsigned l=0; s.init(16); 

   if (m<-1 || p<-1 || t[0]<=0) wblog(FL,
      "ERR %s<%s> invalid type %d.%d%s",FCT,TSTR(wbcomplex),m,p,t);

   for (int i=0; i<2; ++i) {
      l+=snprintf(s.data+l,s.len-l,"%%%s",i?"+":"");
      if (m>0) {
         if (p>=0)
              { l+=snprintf(s.data+l,s.len-l,"%d.%d%s",m,p,t); }
         else { l+=snprintf(s.data+l,s.len-l,"%d%s",m,t); }
      }
      else {
        if (p>=0)
             { l+=snprintf(s.data+l,s.len-l,".%d%s",p,t); }
        else { l+=snprintf(s.data+l,s.len-l,"%s",t); }
      }
   }

   if (l+2>=s.len) wblog(FL,
      "ERR %s<%s> string out of bounds (%%%d.%d%s; %d/%d)",
      FCT,TSTR(wbcomplex),m,p,t,l,s.len
   );

   snprintf(s.data+l,s.len-l,"i");
};

template<class T>
int Wb::num2Fmt<T>::check_init() {
   wblog(FL,"ERR %s<> invalid data type `%s'",FCT,TSTR(T));
   return 1;
};

template<> inline
int Wb::num2Fmt<double>::check_init() {
   if (m<-1) m=8;
   if (p<-1) p=5;
   if (t[0])
        { if (!Wb::strchri("ge",t[0])) return 't'; }
   else { strcpy(t,"g"); }
   return 0;
};
template<> inline
int Wb::num2Fmt<float>::check_init() {
   return ((num2Fmt<double>*)this)->check_init();
};
template<> inline
int Wb::num2Fmt<wbcomplex>::check_init() {
   return ((num2Fmt<double>*)this)->check_init();
};

template<> inline
int Wb::num2Fmt<int>::check_init() {
   if (m<-1) m=4;
   if (p>=0) return 'p'; 
   if (t[0])
        { if (!strchr("ducxX",t[0])) return 't'; } 
   else { strcpy(t,"d"); }
   return 0;
};
template<> inline
int Wb::num2Fmt<unsigned>::check_init() {
   return ((num2Fmt<int>*)this)->check_init();
};

template<> inline
int Wb::num2Fmt<long>::check_init() {
   if (m<-1) m=4;
   if (p>=0) return 'p'; 
   if (t[0])
        { if (!strchr("lducxX",t[0])) return 't'; } 
   else { strcpy(t,"ld"); }
   return 0;
};

template<> inline
int Wb::num2Fmt<unsigned long>::check_init() {
   return ((num2Fmt<long>*)this)->check_init();
};

template<> inline
int Wb::num2Fmt<char>::check_init() {
   if (m<-1) m=4;
   if (p>=0) return 'p'; 
   if (t[0])
        { if (!strchr("ducxX",t[0])) return 't'; } 
   else { strcpy(t,"d"); }
   return 0;
};

template<> inline
int Wb::num2Fmt<char*>::check_init() {
   if (m<-1) m=-1;
   if (p>=0) return 'p'; 
   if (t[0])
        { if (t[0]!='s') return 't'; } 
   else { strcpy(t,"s"); }
   return 0;
};

template <class T> inline
char* defaultFmt(char *fmt,
const T &x __attribute__ ((unused)), int m, int p, char t) {

   const size_t n=16; size_t l=0;
   char dd[n]; dd[0]=0; fmt[0]=0;

   if (m>0) {
        if (p>=0)
             l=snprintf(dd,n,"%d.%d",m,p);
        else l=snprintf(dd,n,"%d",m);
   }
   else if (p>=0) l=snprintf(dd,n,".%d",p);

   if (l>=n) wblog(FL,"ERR %s() string out of bounds (%d/%d)",FCT,l,n);

   if (typeid(T)==typeid(double) || typeid(T)==typeid(float)) {
      if (!t) t='g';
      sprintf(fmt,"%%%s%c",dd,t);
   }
   else if (typeid(T)==typeid(wbcomplex)) {
      if (!t) t='g';
      sprintf(fmt,"%%%s%c%%+%s%ci",dd,t,dd,t);
   }
   else if (
      typeid(T)==typeid(int) || typeid(T)==typeid(char) ||
      typeid(T)==typeid(unsigned)
    ){
      if (p>=0) wblog(FL,"WRN %s() got format *.%d* for type '%s'",
          FCT, p, TSTR(T));
      if (!t) t='d';
      sprintf(fmt,"%%%s%c",dd,t);
   }
   else if (typeid(T)==typeid(char*)) {
       if (t && t!='s') wblog(FL,
          "ERR %s() got format type '%c'<%d> for type '%s'",
           FCT,t,t,TSTR(T)
       );
       sprintf(fmt,"%%%ss",dd);
   }
   else wblog(FL,
     "ERR %s() unsupported type '%s'",FCT,TSTR(T));

   return fmt;
};

wbstring Wb::TimeStamp(char type) {

   wbstring s; 
   unsigned l,n;

  #pragma omp critical (got_CPTR_TIME)
   { time_t t=time(NULL);
     struct tm *tb=localtime(&t);
     if (type=='t' || type=='T') {
        s.init(32); l=strftime(s.data,s.len,"%T",tb);

        if (type=='T') { 
           using namespace std::chrono;
           time_point<system_clock> t=system_clock::now();

           l+=snprintf(s.data+l,s.len-l,".%03ld", 
           size_t(duration_cast<milliseconds>(t.time_since_epoch()).count())%1000);
        }
     }
     else if (type=='D') {
        s.init(32); l=strftime(s.data,s.len,"%d-%b-%Y %Z",tb);
     }
     else { s=asctime(tb); } 
   }

   n=strlen(s.data);
   if (n && s.data[n-1]<32) s.data[n-1]=0; 

   return s;
}

template <class T1, class T2> inline
void Wb::safeConvert(const char *F, int L, const T1 &x1, T2 &x2) {
   x2=T2(x1); if (T1(x2)!=x1) wblog(F_L,
      "ERR converting %s -> %s changes value\n%g, %g",
       TSTR(T1), TSTR(T2),
       double(x1), double(x2)
   );
}

template <> inline
void Wb::safeConvert(const char *F, int L, const wbcomplex &x1, double &x2) {
   if (x1.i!=0.) wblog(F_L,
      "ERR cannot cast complex number to real (%.4g%+.4gi)",x1.r,x1.i);
   x2=x1.r;
}

template <> inline 
void Wb::safeConvert(
   const char *F __attribute__ ((unused)), int L __attribute__ ((unused)),
   const double &x1, wbcomplex &x2
){ x2=wbcomplex(x1,0); }

template <class T> inline
bool Wb::isLower(const T* a, const T* b, const size_t n) {
   for (size_t i=0; i<n; i++) {
       if (a[i]<b[i]) return 1; else
       if (a[i]>b[i]) return 0;
   }
   return 0;
}

template <class T> inline
bool Wb::isEqual(const T* a, const T* b, const size_t n) {

   for (size_t i=0; i<n; i++) { if (a[i]!=b[i]) return 0; }
   return 1;
}

template <class T> inline
bool Wb::allEqual(const T* a, const size_t n, const T &x) {

   for (size_t i=0; i<n; i++) { if (a[i]!=x) return 0; }
   return 1;
}

template <class T> inline
bool Wb::anyUnequal(const T* a, const size_t n, const T &x) {
   for (size_t i=0; i<n; i++) { if (a[i]!=x) return 1; }
   return 0;
};

template <class T> inline
bool Wb::anyEqual(const T* a, const size_t n, const T &x) {

   for (size_t i=0; i<n; i++) { if (a[i]==x) return 1; }
   return 0;
};

template <class T> inline
void Wb::scale_eps(T &eps, const T* d, size_t n) {
   if (double(eps)>0) {
      if (double(eps)>1E-8) wblog(FL,
         "WRN %s() got eps=%g (ignore)",FCT,double(eps));
      else {
         T x=getdscale(d,n);
         if (x>T(1E-8)) eps*=x; 
      }
   }
};

template <class T>
T Wb::getdscale(const T* data, size_t n) {

   if (n) { T d,x=0;
      for (size_t i=0; i<n; ++i) { d=ABS(data[i]); if (x<d) x=d; }
      return x;
   }
   return 0;
};

inline 
void Wb::update_range_stride(size_t &n, const size_t stride) {
   if (stride!=1) {
      if (int(stride)<=0 && n) wblog(FL,
         "ERR %s() got stride=%ld (n=%ld) !?",FCT,stride,n);
      n*=stride;
   }
};

template <class T> inline
T Wb::addRange2(const T* d, size_t n, const size_t stride) {
   T s=0; Wb::update_range_stride(n,stride);
   for (size_t i=0; i<n; i+=stride) { s+=ABS2(d[i]); }
   return s;
};

template <class T> inline
T Wb::addRange(const T* d, size_t n, const size_t stride) {
   T s=0; Wb::update_range_stride(n,stride);
   for (size_t i=0; i<n; i+=stride) { s+=d[i]; }
   return s;
};

template <class T> inline
void Wb::addRange(const T* a, T* c, size_t n, const size_t stride) {
   Wb::update_range_stride(n,stride);
   for (size_t i=0; i<n; i+=stride) { c[i]+=a[i]; }
};

template <class T> inline
void Wb::addRange( 
   const T* a, const T* b, T* c, size_t n, size_t stride
){
   Wb::update_range_stride(n,stride);
   for (size_t i=0; i<n; i+=stride) { c[i]=a[i]+b[i]; }
};

template <class T> inline
void Wb::minusRange( 
   const T* a, const T* b, T* c, size_t n, size_t stride
){
   Wb::update_range_stride(n,stride);
   for (size_t i=0; i<n; i+=stride) { c[i]=a[i]-b[i]; }
};

template <class T> inline
void Wb::diffRange( 
   const T* a, const T* b, T* c, size_t n, size_t stride
){
   Wb::update_range_stride(n,stride);
   for (size_t i=0; i<n; i+=stride) { c[i]=a[i]-b[i]; }
};

template <class T> inline
void Wb::setRange2avg(T* a, size_t n) {
   T x=0; size_t i;
   for (i=0; i<n; ++i) { x+=a[i]; }; x*=(T(1)/n);
   for (i=0; i<n; ++i) { a[i]=x;  }
};

template <class T> inline 
size_t Wb::nnzRange(const T* a, size_t n) {
   size_t m=0, i=0;
   for (; i<n; ++i) { if (a[i]!=0) m++; }
   return m;
};

template <class T> inline
T Wb::rangeNormDiff2(
   const T* d1, const T* d2, size_t n
){
   if (!n || d1==d2) return 0;

   T x2, x2sum=0;

   for (size_t i=0; i<n; ++i) {
      x2=d1[i]-d2[i]; x2*=CONJ(x2);
      x2sum+=x2;
   }

   return x2sum;
};

template <class T> inline
T Wb::rangeNormDiff2(
   const T* d1, const T* d2, size_t n, const T &fac, size_t *k
){
   if (fac==1 && !k) { return rangeNormDiff2(d1,d2,n); }
   if (!n || (d1==d2 && fac==1)) { if (k) { (*k)=0; }; return 0; }

   T x2, x2sum=0;

   if (k) { T x2max=0;
      for (size_t i=0; i<n; ++i) {
         x2=d1[i]-fac*d2[i]; x2*=CONJ(x2);
         x2sum+=x2; if (x2>x2max) { x2max=x2; (*k)=i; }
      }
   }
   else {
      for (size_t i=0; i<n; ++i) {
         x2=d1[i]-fac*d2[i]; x2*=CONJ(x2);
         x2sum+=x2;
      }
   }

   return x2sum;
};

template <class T> inline
T Wb::rangeNorm2(const T* d, size_t n, size_t *k
){
   T x2sum=0;
   if (!n) { if (k) { (*k)=0; }; return x2sum; }

   if (k) {
      T x2, x2max=0; (*k)=0;

      for (size_t i=0; i<n; ++i) {
         x2=CONJ(d[i])*d[i]; if (x2>x2max) { x2max=x2; (*k)=i; }
         x2sum+=x2;
      }
   }
   else {
      for (size_t i=0; i<n; ++i) { x2sum+=CONJ(d[i])*d[i]; }
   }

   return x2sum;
};

template <class T> inline 
T Wb::rangeMaxDiff(const T* d1, const T* d2, size_t n, size_t *k) {
   size_t imax=0, i=0; T x,xmax=0;

   if (!n || d1==d2) { if (k) (*k)=0;
      if (!n) wblog(FL,"WRN %s() got empty input",FCT);
      return 0;
   }

   for (; i<n; i++) {
      x=ABS(d1[i]-d2[i]); if (x>xmax) { xmax=x; imax=i; }
   }

   if (k) (*k)=imax;
   return xmax;
};

template <class T, class T2>
void Wb::cpyStride(
   T2* d2,         
   const T* d0,    
   unsigned n,     
   const unsigned *idx, unsigned m, 
   unsigned D2,    
   unsigned D0,    
   char add
){
   if (!n || !m) { return; }
   if (!idx) { return Wb::cpyStride(d2,d0,n,m,D2,D0, T2(add? 1:0), T(1)); }

   unsigned i=0, j; const T *d0_=d0;

   if (int(D0)<0) { D0=n; } else if (D0<n) { ++i; }
   if (int(D2)<0) { D2=n; } else if (D2<n) { ++i; }

   if (i) wblog(FL,
      "ERR %s() stride too small (%d,%d/%d)",FCT,D0,D2,n);
   if (gotMemOverlap(d0, (m-1)*D0+n, d2, (m-1)*D2+n)) wblog(FL,
      "ERR %s() must not copy onto itself",FCT); 

   for (j=0; j<m; ++j, d2+=D2) { d0 = d0_+ idx[j] * D0;
      if (!add) { for (i=0; i<n; ++i) { d2[i]  = d0[i]; }}
      else      { for (i=0; i<n; ++i) { d2[i] += d0[i]; }}
   }
};

template <class T, class T2>
void Wb::cpyStride(
   T2* d2,       
   const T* d0,  
   size_t n,     
   size_t m,     
   size_t D2,    
   size_t D0,    
   T2 fac2, T fac0
){
   if (!n || !m) { return; }

   unsigned j=0;
   if (int(D0)<0) { D0=n; } else if (D0<n) { ++j; }
   if (int(D2)<0) { D2=n; } else if (D2<n) { ++j; }

   if (j) wblog(FL,
      "ERR %s() stride too small (%d,%d/%d)",FCT,D0,D2,n);
   if (gotMemOverlap(d0, (m-1)*D0+n, d2, (m-1)*D2+n)) wblog(FL,
      "ERR %s() must not copy onto itself",FCT);

   for (j=0; j<m; ++j, d2+=D2, d0+=D0) {
      Wb::cpyRange(d2,d0,n,fac2,fac0);
   }
};

template <class T, class TB>
void Wb::cpyRange(T* a, const TB* b, size_t n, T afac, TB bfac) {

   size_t i=0;

   if (bfac) {
      if (!afac) {
         if (bfac==TB( 1)) { for (; i<n; ++i) { a[i] = b[i]; }} else
         if (bfac==TB(-1)) { for (; i<n; ++i) { a[i] =-b[i]; }}
         else              { for (; i<n; ++i) { a[i] = b[i]*bfac; }}
      }
      else if (afac==T(1)) {
         if (bfac==TB( 1)) { for (; i<n; ++i) { a[i]+= b[i]; }} else
         if (bfac==TB(-1)) { for (; i<n; ++i) { a[i]-= b[i]; }}
         else              { for (; i<n; ++i) { a[i]+= b[i]*bfac; }}
      }
      else {
         if (bfac==TB( 1)) { for (; i<n; ++i) { (a[i]*=afac) += b[i]; }} else
         if (bfac==TB(-1)) { for (; i<n; ++i) { (a[i]*=afac) -= b[i]; }}
         else              { for (; i<n; ++i) { (a[i]*=afac) += b[i]*bfac; }}
      }
   }
   else if (afac!=T(1)) {
      if (afac) {    for (; i<n; ++i) { a[i]*=afac; }}
      else { T z=0;  for (; i<n; ++i) { a[i]=z;     }}
   }
};

template <class T>
void Wb::cpyRange(T* a, const T* b, size_t n) { 
   for (size_t i=0; i<n; ++i) { a[i]=b[i]; }
};

template <class T, class Tb>
void Wb::cpyRange(T* a, const Tb* b, size_t n, char tcheck) { 

   if (!tcheck) {
      for (size_t i=0; i<n; ++i) { a[i]=T(b[i]); }} 
   else {
      for (size_t i=0; i<n; ++i) { a[i]=T(b[i]);
         if (Tb(a[i])!=b[i]) wblog(FL,
            "ERR %s() rounding error (%g/%g)",FCT,double(a[i]),double(b[i])
         );
      }
   }
};

template <>
void Wb::cpyRange(double* a, const double* b, size_t n) { 
   memcpy(a,b,n*sizeof(double));
};

template <>
void Wb::cpyRange(wbcomplex* a, const wbcomplex* b, size_t n) { 
   memcpy(a,b,n*sizeof(wbcomplex));
};

template <>
void Wb::cpyRange(int* a, const int* b, size_t n) { 
   memcpy(a,b,n*sizeof(int));
};

template <> 
void Wb::cpyRange(double* a, const wbcomplex* b, size_t n, char tcheck) {

   double i2=0;
   for (size_t i=0; i<n; ++i) {
       a[i]=b[i].r; if (b[i].i) { i2+=(b[i].i*b[i].i); }
   }

   if (i2) {
      sprintf(str,"%s() got |Im(b)|^2 = %.3g",FCT,i2);
      wblog(FL,tcheck? "ERR %s":"WRN %s",str);
   }
};

template <> 
void Wb::cpyRange(wbcomplex* a, const double* b, size_t n, char tcheck) {
   for (size_t i=0; i<n; ++i) { a[i]=b[i]; }
};

template <class T>
void Wb::cpyRangeR(T* a, const wbcomplex* z, size_t n) { 
   for (size_t i=0; i<n; ++i) { a[i]=z[i].r; }
};

template <class T>
void Wb::cpyRangeI(T* a, const wbcomplex* z, size_t n) { 
   for (size_t i=0; i<n; ++i) { a[i]=z[i].i; }
};

template <class T>
void Wb::cpyRangeA(T* a, const wbcomplex* z, size_t n) { 
   for (size_t i=0; i<n; ++i) { a[i]=z[i].abs(); }
};

template <class T> inline
bool Wb::uniformRange(const T* d, size_t n) {
   for (size_t i=1; i<n; ++i) if (d[i]!=d[0]) return 0;
   return 1;
};

template <class T> inline
T Wb::maxRange(const T* d, size_t n) {

   if (!n) wblog(FL,"ERR %s() got empty range",FCT);

   T x=d[0];
   for (size_t i=1; i<n; ++i) { if (x<d[i]) { x=d[i]; }}
   return x;
};

template <class T> inline
T Wb::minRange(const T* d, size_t n) {

   if (!n) wblog(FL,"ERR %s() got empty range",FCT);

   T x=d[0];
   for (size_t i=1; i<n; ++i) { if (x>d[i]) { x=d[i]; }}
   return x;
};

template <class T> inline
T Wb::prodRange(const T* d, size_t n) {
   if (!n) wblog(FL,"WRN %s() got empty range (returning 0)",FCT);
   return prodRange(d,n,T(0));
};

template <class T> inline
T Wb::prodRange(const T* d, size_t n, T x) {
   if (n) { size_t i=1;
      if (!d) wblog(FL,"ERR %s() got null space (n=%d)",FCT,n);
      for (x=d[0]; i<n; ++i) x*=d[i];
   }
   return x;
};

template <class TD, class TX> inline
void Wb::timesRange(TD* d, TX x, size_t n, size_t stride) {
   if (!n) return;
   if (!d) wblog(FL,"ERR %s() got null space (n=%d)",FCT,n);

   Wb::update_range_stride(n,stride);
   for (size_t i=0; i<n; i+=stride) { d[i]*=x; }
};

template <class T> inline
char Wb::check_conj_flag(const T* a, char &conj) {
   if (conj) {
      if (!ISCOMPLX_(T)) { conj=0; } else {
      if (conj>2) { 
         if (strchr("2bBcC*",conj)) { conj=2; } else
         if (strchr("1aA",   conj)) { conj=1; }
      }
      if (conj<0 || conj>2) wblog(FL,
         "WRN %s() invalid conj=%s",FCT,cSTR(conj));
      }
   }
   return conj;
};

template <class T> inline
T Wb::dotProd(const T* a, const T* b, size_t n, T x, char conj) {
   if (n) {
      size_t i=0; x=0;
      if (!a || !b) wblog(FL,"ERR %s() got null (%p, %p)",FCT,a,b);

      if (check_conj_flag(b,conj)) { 
         if (conj==1)
              { for (; i<n; ++i) x+=CONJ(a[i])*b[i]; }
         else { for (; i<n; ++i) x+=a[i]*CONJ(b[i]); }
      }
      else { for (; i<n; ++i) x+=(a[i]*b[i]); }
   }
   return x;
};

template <class T> inline
T Wb::dotProd(const T* a, const T* b, size_t n, char conj) {
   T x=T(0); 

   if (n) { x=dotProd(a,b,n,x,conj); }
   else {
      wblog(FL,"WRN %s() got empty range",FCT);
   }

   return x;
};

template <class T> inline
void Wb::TimesElRange(T* a, const T* b, size_t n, char conj) {

   if (!n) { return; }

   if (!a || !b) wblog(FL,
      "ERR %s() got null space (%p, %p)",FCT,a,b);

   if (check_conj_flag(b,conj)) {
      if (conj==1)
           { for (size_t i=0; i<n; ++i) a[i]=CONJ(a[i])*b[i]; }
      else { for (size_t i=0; i<n; ++i) a[i]*=CONJ(b[i]); }
   }
   else { for (size_t i=0; i<n; ++i) a[i]*=b[i]; }
};

template <class T> inline
void Wb::timesElRange(T* c, const T* a, const T* b, size_t n, char conj) {

   if (!n) { return; }
   if (!a || !b || !c) wblog(FL,
      "ERR %s() got null space (%p, %p, %p)",FCT,a,b,c);
   size_t i=0;

   if (check_conj_flag(b,conj)) {
      if (conj==1)
           { for (; i<n; ++i) { c[i] = CONJ(a[i]) *b [i]; }}
      else { for (; i<n; ++i) { c[i] = a[i] * CONJ(b[i]); }}
   }
   else { for (; i<n; ++i) { c[i] = a[i]*b[i]; }}
};

template <class T> inline
void Wb::timesElRange_add(
   T* c, const T* a, const T* b, size_t n, T bfac, char conj
){
   size_t i=0;

   if (!n || !bfac) { return; }
   if (!a || !b || !c) wblog(FL,
      "ERR %s() got null space (%p, %p, %p)",FCT,a,b,c);

   if (check_conj_flag(b,conj)) {
      if (conj==1) {
         if (bfac==T(+1)) { for (; i<n; ++i) c[i]+=CONJ(a[i])*b[i]; } else
         if (bfac==T(-1)) { for (; i<n; ++i) c[i]-=CONJ(a[i])*b[i]; }
         else             { for (; i<n; ++i) c[i]+=CONJ(a[i])*b[i]*bfac; }
      }
      else { 
         if (bfac==T(+1)) { for (; i<n; ++i) c[i]+=a[i]*CONJ(b[i]); } else
         if (bfac==T(-1)) { for (; i<n; ++i) c[i]-=a[i]*CONJ(b[i]); }
         else             { for (; i<n; ++i) c[i]+=a[i]*CONJ(b[i])*bfac; }
      }
   }
   else {if (bfac==T(+1)) { for (; i<n; ++i) c[i]+=a[i]*b[i]; } else
         if (bfac==T(-1)) { for (; i<n; ++i) c[i]-=a[i]*b[i]; }
         else             { for (; i<n; ++i) c[i]+=a[i]*b[i]*bfac; }
   }
};

template <class T> inline
void Wb::timesElRange_OM(T* c, 
   const T* a, const T* b, 
   size_t N, unsigned M, char conj
){
   if (!N || !M) { wblog(FL,
      "WRN %s() got empty NxM=%dx%d - return",FCT,N,M); return; }
   if (!a || !b || !c) wblog(FL,
      "ERR %s() got null (%p, %p, %p)",FCT,a,b,c);

   unsigned m=0; size_t i=0;
   check_conj_flag(b,conj);

   for (; m<M; ++m) {
      if (m) { i=0; a+=N; b+=N;
         if (!conj)
              { for (; i<N; ++i) { c[i] += a[i]*b[i]; }} else
         if (conj==1)
              { for (; i<N; ++i) { c[i] += CONJ(a[i]) *b[i];  }}
         else { for (; i<N; ++i) { c[i] += a[i] * CONJ(b[i]); }}
      }
      else {
         if (!conj)
              { for (; i<N; ++i) { c[i]  = a[i]*b[i]; }} else
         if (conj==1)
              { for (; i<N; ++i) { c[i]  = CONJ(a[i]) *b[i];  }}
         else { for (; i<N; ++i) { c[i]  = a[i] * CONJ(b[i]); }}
      }
   }
};

template<class T> inline
T Wb::overlap(
   const T* a, const T* b, size_t n,
   size_t stride, 
   char tnorm __attribute__ ((unused))
){
   T x2=0; 
   if (n) {
      Wb::update_range_stride(n,stride);
      if (!a || !b) wblog(FL,"ERR %s() got null (%p, %p)",FCT,a,b);
      for (size_t i=0; i<n; i+=stride) { x2+=CONJ(a[i])*b[i]; }
   }
   return x2;
};

template<> inline
wbcomplex Wb::overlap(
   const wbcomplex* a, const wbcomplex* b, size_t n,
   size_t stride, char tnorm 
){
   wbcomplex z2=0;
   if (n) { size_t i=0;
      Wb::update_range_stride(n,stride);
      if (!a || !b) wblog(FL,"ERR %s() got null (%p, %p)",FCT,a,b);
      if (tnorm) { for (; i<n; i+=stride) { z2+=a[i]*b[i]; }} else
      if (a!=b ) { for (; i<n; i+=stride) { z2+=a[i].conj()*b[i]; }}
      else { double nrm2=0;
         for (; i<n; i+=stride) { nrm2+=a[i].abs2(); }
         z2=nrm2;
      }
   }
   return z2;
};

template<class T>
T Wb::gs_project_range(
   T *a, const T *b, 
   size_t n, size_t stride, char isnorm, char tnorm
){
   T x=overlap(b,a,n,stride,tnorm);

   if (!isnorm) { 
      T nrm2=overlap(b,b,n,stride,tnorm);
      if (std::fabs(double(nrm2))<1E-12) wblog(FL,
         "WRN %s() got norm=%.3g !?",FCT,double(nrm2));
      x/=nrm2;
   }

   Wb::update_range_stride(n,stride);
   for (size_t i=0; i<n; i+=stride) { a[i]-=x*b[i]; }

   return x;
};

template <class T> inline 
size_t Wb::replRange(T* a, size_t n, T x, T v) {

   size_t i,m=0;

   if (!ISNAN(x))
        { for (i=0; i<n; ++i) if (a[i]==x)     { a[i]=v; ++m; }}
   else { for (i=0; i<n; ++i) if (ISNAN(a[i])) { a[i]=v; ++m; }}

   return m;
}

inline size_t Wb::countNaN(double* a, size_t m) {
   size_t i,n=0;
   for (i=0; i<m; i++) { if (std::isnan(a[i])) ++n; }
   return n;
};

   template <class T>
   void Wb::cpyZRange(
      const double *R, const double *I, T *Z, size_t n
   ){
      wblog(FL,"ERR %s not applicable for type '%s'\n(%p,%p,%p,%d)",
      FCT, TSTR(T),R,I,Z,n);
   };

   template <>
   void Wb::cpyZRange(
      const double *R, const double *I, wbcomplex *Z, size_t n
   ){
      for (size_t i=0; i<n; i++)
      Z[i].set(R[i], I ? I[i] : 0);
   }

   template <class T>
   void Wb::splitZRange(
      const T *Z, double *R, double *I, size_t n
   ){
      wblog(FL,"ERR %s not applicable for type `%s'",
      FCT, TSTR(T));
   }

   template <>
   void Wb::splitZRange(
      const wbcomplex *Z, double *R, double *I, size_t n
   ){
      for (size_t i=0; i<n; i++) {
      R[i]=Z[i].r; if (I) I[i]=Z[i].i; }
   }

template<> 
void Wb::chopTiny_float(double *d, size_t n, double ref) {

   size_t i; if (!n) return;

   if (ref<0) { 
      for (ref=0, i=0; i<n; ++i) ref=MAX(ref,::fabs(d[i]));
      if (ref==0) return;
   }

   if (ref!=0) {
      int k=int(std::log(ref)/std::log(2)); 
      if (k>=0) ref=(1<<k); else ref=1/(1<<-k);

      for (i=0; i<n; ++i) { 
         if (d[i]>0.) d[i]=float(d[i]+ref)-ref; else
         if (d[i]<0.) d[i]=float(d[i]-ref)+ref;
      }
   }
   else {
      for (i=0; i<n; ++i)
      if (d[i]!=0.) d[i]=double(float(d[i]));
   }
};

template<> 
double Wb::chopTiny_imag(wbcomplex *z, size_t n, double eps) {

   double i2=0, r2=0; size_t k=0; 

   for (; k<n; ++k) { r2+=(z[k].r*z[k].r); i2+=(z[k].i*z[k].i); }
   if (i2) {
      if (r2>1) { eps*=std::sqrt(r2); }
      if ((i2=std::sqrt(i2))>eps) { i2=-i2; } 
      else {
         for (k=0; k<n; ++k) { z[k].i=0; }
      }
   }
   return i2;
};

template<> 
void Wb::chopTiny_z(wbcomplex *z, size_t n, double eps) {

   size_t k=0; double a2, r2=0, i2=0;

   for (; k<n; ++k) {
      r2+=(z[k].r*z[k].r);
      i2+=(z[k].i*z[k].i);
   }
   a2=r2+i2; if (a2>1) { eps*=std::sqrt(a2); }

   if (std::sqrt(i2)<=eps) { for (k=0; k<n; ++k) { z[k].i=0; }} else
   if (std::sqrt(r2)<=eps) { for (k=0; k<n; ++k) { z[k].r=0; }}
};

template<class T>
void Wb::invertIndex(
    const wbvector<T> &i1, unsigned N,
    wbvector<T> &i2, char lflag
){
   wbvector<char> flag(N);
   unsigned i,k,n;

   for (i=0; i<i1.len; i++) {
       k=i1[i]; if (k<N) flag[k]=1; else
       wblog(FL, "ERR index out of bounds (%d/%d)", k, N);
   }
   n=flag.isum();

   if (n!=i1.len && lflag)
   wblog(FL, "WRN index set not unique!? (%d/%d)", n, i1.len);

   i2.init(N-n);
   for (k=i=0; i<N; i++) if (flag[i]==0) i2[k++]=i;
}

void Wb::setRand(wbvector<wbcomplex> &zz, double fac, double shift) {
   size_t i=0; fac/=(double)RAND_MAX;

   if (shift==0.)
   for (; i<zz.len; ++i)
        zz.data[i].set(fac*std::rand(), fac*std::rand());
   else
   for (; i<zz.len; ++i)
        zz.data[i].set(fac*std::rand()+shift,fac*std::rand()+shift);
}

template<class T>
void Wb::set2avg(T *dd, const size_t n) {
   size_t i;  T dbl=0; if (!n) return;
   for (i=0; i<n; ++i) { dbl+=dd[i]; }; dbl/= n;
   for (i=0; i<n; ++i) { dd[i]=dbl;  }
};

template<class T>
void Wb::getDiff(const wbvector<T> &xx, wbvector<T> &dx) {
   size_t i, n=xx.len; dx.init(n);
   for (i=1; i<n; i++) {
       dx[i-1] += (dx[i] = xx[i]-xx[i-1]);
       if (i>1) dx[i-1] *= 0.5;
   }
}

template<class TX, class TY>
TY Wb::IntTrapez(const TX *xx, const TY *yy, const size_t n) {
   TY s=0;

   for (size_t i=1; i<n; ++i) {
      s += 0.5 * (xx[i] - xx[i-1]) * (yy[i] + yy[i-1]);
   }

   return s;
};

template<class T>
void markSet(
   wbvector< wbvector<T>* > &E0, 
   wbvector< wbvector<char> > &mark, 
   unsigned Nkmin,   
   int &Nkeep,       
   double Etrunc,    
   wbvector<T> &E,   
   double eps,       
   double b,         
   int dmax,         
   const char *dir   
){
   unsigned i,d,l,N,n=E0.len;
   wbvector<char> mm;
   wbperm P,iP;

   for (N=i=0; i<n; i++) N+=(E0[i]->len);

   E.init(N); mark.initDef(n);
   if (!N) return;

   for (l=i=0; i<n; i++, l+=d) {
      d=(E0[i]->len);  mark[i].init(d);
      memcpy(E.data+l, E0[i]->data, d*sizeof(T));
   }

   E.Sort(P); P.invert(iP);

   if (eps>0) {
      unsigned i0=0;

      eps *= ABS(E.last()-E[0]); 
      for (i=1; i<=N; i++) {
         if (i==N || ABS(E[i]-E[i-1])>eps) {
            d=i-i0; if (d>1)
            Wb::setRange2avg(E.data+i0,d);
            i0=i;
         }
      }

      wbvector<T> x; E.permute(x,iP);
      for (l=i=0; i<n; ++i, l+=d) {
         d=E0[i]->len;
         memcpy(E0[i]->data, x.data+l, d*sizeof(T));
      }
   }

   if (Nkeep<0) Nkeep=N; 
   if (Nkeep==0) return;
   if (Nkeep>=(int)N && Etrunc<=0) {
      for (i=0; i<n; i++) mark[i].set(1);
      Nkeep=N; return;
   }

   markSet(FL,E,mm,Nkmin,Nkeep,Etrunc,b,dmax,dir);

   mm.Permute(iP);

   for (l=i=0; i<n; i++, l+=d) {
      d=mark[i].len;
      memcpy(mark[i].data, mm.data+l, d*sizeof(char));
   }

};

template<>
void markSet(
   wbvector< wbvector<wbcomplex>* > &E0, 
   wbvector< wbvector<char> > &mark,
   unsigned Nkmin,   
   int &Nkeep,       
   double Etrunc,    
   wbvector<wbcomplex> &E,   
   double eps,       
   double b,         
   int dmax,         
   const char *dir_  
){
   unsigned i,d,l,N, n=E0.len;
   const char *s=NULL;

   wbvector<char> mm;
   wbvector<double> ee;
   wbperm P,iP;

   if (dir_ && dir_[0]) {
      if (!strncmp(dir_,"asc", 3)) { s=dir_+3; } else
      if (!strncmp(dir_,"desc",4)) { s=dir_+4; } else
      wblog(FL,"ERR %s() invalid direction >%s<",FCT,dir_);
   }

   for (N=i=0; i<n; i++) N+=(E0[i]->len);
   if (!N) return;

   E.init(N); ee.init(N); mark.initDef(n);

   for (l=i=0; i<n; i++, l+=d) {
      d=(E0[i]->len);  mark[i].init(d);
      memcpy(E.data+l, E0[i]->data, d*sizeof(wbcomplex));
   }

   if (s==NULL || !strcmp(s,"R")) {
      Wb::cpyRangeR(ee.data,E.data,N);
   }
   else if (!strcmp(s,"I")) {
      Wb::cpyRangeI(ee.data,E.data,N);
   }
   else if (!strcmp(s,"A")) {
      wbcomplex z=E.min();
      wbvector<wbcomplex> x(E); z.i=0; x-=z;
      Wb::cpyRangeA(ee.data,x.data,N);
   }
   else wblog(FL,"ERR markSet() invalid direction >%s<", dir_);

   ee.Sort(P); E.Select(P); P.invert(iP);

   if (Nkeep<0) Nkeep=N; 
   if (Nkeep==0) return;

   if (Nkeep>=(int)N && Etrunc<=0) {
      for (i=0; i<n; ++i) mark[i].set(1);
      Nkeep=N; return;
   }

   markSet(FL,ee,mm,Nkmin,Nkeep,Etrunc,MAX(b,eps),dmax,dir_);

   mm.Permute(iP);

   for (l=i=0; i<n; i++, l+=d) {
      d=mark[i].len;
      memcpy(mark[i].data, mm.data+l, d*sizeof(char));
   }
}

template<class T>
void markSet(const char *F, int L,
   const wbvector<T> &E, 
   wbvector<char> &mark,
   unsigned Nkmin,   
   int &Nkeep,       
   double Etrunc,    
   double b,         
   int dmax,         
   const char *dir_  
){

   unsigned i,imax,l;
   unsigned N=E.len, nx;
   double x, xmax=0;
   char dir=+1; 

   if (dir_ && dir_[0]) {
      if (!strncmp(dir_,"asc", 3)) dir=+1; else  
      if (!strncmp(dir_,"desc",4)) dir=-1; else
      wblog(F,L,"ERR %s() invalid direction >%s<",FCT,dir_);
   }

   mark.init(N);

   if (Nkeep<0 || Nkeep>int(N)) Nkeep=N; 
   if (Nkeep==0) return;

   if (Nkeep>=(int)N && Etrunc<=0) { mark.set(1); Nkeep=N; return; }
   if (N<1) wblog(FL,"ERR %s() N=%d !?",FCT,N); 

   if (int(Nkmin)<0) Nkmin=0; 

   if (dir>0) {
      if (Etrunc>0) { 
         x=Etrunc+E[0]; 
         for (l=1; l<N; ++l) {
            if (E[l-1]>E[l]) wblog(FL,"ERR %s() data not sorted!",FCT);
            if (E[l]>=x && l>Nkmin) break;
         }
         if (int(l)>Nkeep) l=Nkeep; 
      }
      else { l=unsigned(Nkeep); }

      if (l==N) { mark.set(1); Nkeep=N; return; }
      if (l==0) { Nkeep=0; return; }

      nx=(dmax>=0 ? (unsigned)dmax : MAX(16U,l/5));

      if (b>0) {
         for (imax=i=l; i<N && i<l+nx; ++i) {
            x=fabs(E[i]-E[i-1]);
            if (x<=b) { if (xmax<x) { xmax=x; imax=i; }}
            else break;
         }
         if (i>=l+nx) wblog(F,L,
            "WRN %s() reached Nkeep+dmax=%d+%d states\n"
            "(db=%.3g too small !? => %d @ %.3g)", FCT,l,nx,b,imax,xmax);

         Nkeep=l=imax;
      }
      else if (b<0 && nx && l<N) {

         double q=10/(-b*nx);
         unsigned NX=5/q, i2=MIN(N,l+NX);
         i=l-NX; if (int(i)<1) i=1; 

         for (imax=l; i<i2; ++i) { x=q*int(i-l);
            x=exp(-x*x)*fabs(E[i]-E[i-1]);
            if (x>xmax) { xmax=x; imax=i; }
         }
         if (i==N) { x=q*int(i-l);
            x=0.1*exp(-x*x)*fabs(E[N-1]-E[0]);
            if (x>xmax) { xmax=x; imax=i; }
         }

         Nkeep=l=imax;
      }

      for (i=0; i<l; ++i) mark[i]=1;
   }
   else {
      if (Etrunc>0) { 
         unsigned lx=(Nkmin<N ? N-Nkmin : 0); 

         x=Etrunc; 

         for (l=N-1; l<N; --l) {
            if (l && E[l-1]>E[l]) wblog(FL,"ERR %s() data not sorted!",FCT);
            if (E[l]<=x && l<lx) break;
         }
         if (int(N-1-l)>Nkeep) { l = N-1-(unsigned)Nkeep; }
      }
      else { l=N-1-(unsigned)Nkeep; }

      if (l>N ) { mark.set(1); Nkeep=N; return; }
      if (l==N-1) { Nkeep=0; return; }

      nx=(dmax>=0 ? (unsigned)dmax : MAX(16U,(N-l)/5));

      if (b>0) { 
         for (imax=i=l; i<N && (l-i)<nx; --i) {
            x=fabs(E[i+1]-E[i]);
            if (x<=b) { if (xmax<x) { xmax=x; imax=i; }}
            else break;
         }
         if ((l-i)>=nx) wblog(F,L,
            "WRN %s() reached Nkeep+dmax=%d+%d states (%d/%d)\n"
            "(db=%.3g too small !? => %d @ %.3g)",FCT,
            Nkeep,nx, Nkeep+nx,N, b,N-imax,xmax);
         l=imax;
      }
      else if (b<0 && nx && l<N ) {

         double q=10/(-b*nx);
         unsigned NX=5/q, i1=l-NX; if (int(i1)<0) i1=0;
         i=MIN(N-1,l+NX); 

         for (imax=l; i>i1; --i) { x=q*int(i-l);
            x=exp(-x*x)*fabs(E[i+1]-E[i]);
            if (xmax<x) { xmax=x; imax=i; }
         }
         if (i==0) { x=q*int(i-l);
            x=0.1*exp(-x*x)*fabs(E[N-1]-E[0]);
            if (x>xmax) { xmax=x; imax=i-1; }
         }
         l=imax;
      }

      if (l>N) i=0; 
      else { i=l+1; }; Nkeep=(l<N ? N-1-l : 0);

      for (; i<N; ++i) mark[i]=1;
   }
};

#endif

