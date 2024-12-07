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

#ifndef __WB_UTIL_CC__
#define __WB_UTIL_CC__

/* -------------------------------------------------------------------- */
// NB! this requires
//    fabs() => [std]::fabs
//    sqrt() => [std]::sqrt(), etc.
// since Wb::fabs() is only defined for MPFR (!)
// and so compiler complains that MPFR routine got double // Wb,Jul05,22
/* -------------------------------------------------------------------- */

int Wb::Rational( 
   double &x0,          
   long &P, long &Q,    
   double *dx,          
   unsigned *nmax,      
   wbvector<double>* aa,
   unsigned niter,      
   long pqmax,

   double eps,  

   double eps2, 

   char vflag
){
   unsigned i; int err=0;
   long p[niter+2], q[niter+2];
   double e,xi,r=0, x=x0, a[niter+1];

   if (eps<eps2 || eps2<=0) wblog(FL,
      "ERR %s() invalid eps=%g, eps2=%g",FCT,eps,eps2);

   if (pqmax<=0) { pqmax=long(MIN(1E12, 1/eps)); }
   if (pqmax<=0 || pqmax>(1L<<50)) {
      wblog(FL,"ERR %s() pqmax=%ld !?",FCT,pqmax);
   }

   if (::fabs(x0)>pqmax || (x!=0 && ::fabs(1/x0)>pqmax)) {
      if (vflag) wblog(FL,"WRN %s() double out of bounds long int",FCT);
      if (dx) { (*dx)=0; }
      if (nmax) { (*nmax)=0; }; if (aa) { aa->init(); }
      return -1; 
   }

   p[0]=1; p[1]=a[0]=xi=::round(x); 
   q[0]=0; q[1]=1;

   if (vflag) printf(
      "\n  %s() finding rational approximation for %.16g ...\n"
      "\n%6d: %4g   | %8ld / %8ld\n",FCT,x, 0,a[0],p[1],q[1]);

   if (::fabs(xi-x)<eps) {
      P=p[1]; Q=q[1];
      if (dx) { (*dx)=::fabs(xi-x); }
      if (nmax) { (*nmax)=0; }; if (aa) { aa->init(1,a); }
      x0=xi; return 0;
   }

   if (::fabs(x)>1) {
      eps *=::fabs(x);
   }

   r=::fabs(xi-x0);

   for (i=0; i<niter; ++i) {
      if ((e=::fabs(x-a[i]))<eps || r<eps2) break;

      x=1./(x-a[i]); a[i+1]=::round(x); 
      p[i+2]=a[i+1]*p[i+1]+p[i];
      q[i+2]=a[i+1]*q[i+1]+q[i];

      xi=double(p[i+2])/q[i+2]; 
      r=::fabs(xi-x0);

      if (vflag) {
         printf("%6d: %4g   | %8ld / %8ld = %20.16g  @  %8.3g (%.3g)\n",
         i+1, a[i+1], p[i+2],q[i+2], xi,r, x-a[i+1]);
      }
      if (labs(p[i+2])>pqmax || labs(q[i+2])>pqmax) {
         ++i; err=1; break;
      }
   }
   P=p[i+1]; Q=q[i+1]; if (Q<0) { P=-P; Q=-Q; }
   if (i>=niter) err=2;

   if (vflag) {
      if (i<niter) printf("\n  "
         "Rational representation: %ld/%ld\n\n",P,Q);
      else printf("\n  "
         "Failed to find rational representation (%d @ %g)\n\n",
          niter,eps
      );
   }

   if (dx) (*dx)=r;
   if (aa) aa->init(i,a);
   if (nmax) (*nmax)=i;

   if (!err) x0=xi;
   return err;
};

template<class T>
double Wb::FixRational(
   const char *F  __attribute__ ((unused)),
   int L          __attribute__ ((unused)),
   T *d           __attribute__ ((unused)),
   unsigned n     __attribute__ ((unused)),
   char rflag     __attribute__ ((unused)),
   unsigned niter __attribute__ ((unused)),
   long pqmax     __attribute__ ((unused)),
   double eps     __attribute__ ((unused)),
   double eps2    __attribute__ ((unused)),
   double *rz     __attribute__ ((unused)),
   double *ra     __attribute__ ((unused)),
   char vflag     __attribute__ ((unused))
) {
   return 0;
};

template<> 
double Wb::FixRational( 
   const char *F, int L, double *A, unsigned n,
   char rflag, unsigned niter,
   long pqmax,     
   double eps,     
   double eps2,    
   double *r2z,    
   double *r2m,    
   char vflag
){
   unsigned i,m,m2, nmax=0, nmissed=0; int e1,e2; long p,q, p2,q2;
   double x,Ai,Ai2,a2,d2,a=0,r=0,r2=0,r2max=0,z2=0,err2=0, epsi=eps;

   for (i=0; i<n; ++i) {
      x=::fabs(A[i]); if (a<x) { a=x; }

      x=::fabs(A[i]-::round(A[i])); if (r<x) { r=x; }
      if (x<eps) { d2=(x*x);
         A[i]=::round(A[i]);
         if (A[i]) { err2+=d2; } else { z2+=d2; }
         if (r2max<d2) { r2max=d2; }
      }
   }
   if (r2z) { (*r2z)=z2; } else { err2+=z2; }

   if (r<eps || int(niter)<1) {
      if (r2m) { (*r2m)=r2max; }
      return err2;
   }
   if (a>1) { epsi*=a; }

   for (i=0; i<n; ++i) { a=Ai=A[i]; e2=-1;
       e1=Rational(a,p,q,&r,&m,NULL,niter,pqmax,eps,eps2,vflag>1);

       if (rflag && (e1 || abs(p)>999 || abs(q)>999)) {
          if ( (a2=Ai2=Ai*Ai) > eps) {
		  e2=Rational(a2,p2,q2,&r2,&m2,NULL,niter,pqmax,
			 a2<0.01 ? eps *::fabs(Ai) : eps,
			 a2<0.01 ? eps2*::fabs(Ai) : eps2, vflag>1);
		  }
       }

       if (!e1 && !e2) {
          double q1=::fabs(A[i]*A[i]-Ai2), q2=::fabs(a2-Ai2);
          if (q1>q2) { if (vflag>1) {
             PRINTF("  ==> picking sqrt() over plain (e=%.3g / %.3g)\n",
                ::sqrt(q2), ::sqrt(q1)); }
             e1=-2;
          }
       }

       if (!e1) {
          if (nmax<m) { nmax=m; }
          if (r>epsi) { if (q<0) { p=-p; q=-q; }; wblog(F_L,
             "WRN fixing number %.16g = %ld/%ld (%d; %.3g)",Ai,p,q,m,r);
          }
          A[i]=a; 
          d2=r*r; err2+=d2; if (r2max<d2) { r2max=d2; }
          continue;
       }
       if (e2) {  
          nmissed++; continue;
       }

       if (nmax<m2) { nmax=m2; }
       a=::sqrt(a2); if (Ai<0) { a=-a; }
       r=::fabs(a-Ai);

       if (r>epsi) {
          wblog(F_L,"WRN fixing %.16g = %ssqrt(%ld/%ld) (%.2g)",
            A[i], A[i]<0 ? "-":"", p2,q2,r);
          wblog(FL,"TST %s() pqmax=%ld, eps=%g, eps2=%g,niter=%d",
            FCT,pqmax,eps,eps2,niter);
          double a2=Ai*Ai; 
          Rational(a2,p2,q2,&r2,&m2,NULL,niter,pqmax,eps,eps2,'v');
       }
       d2=r*r; err2+=d2; if (r2max<d2) { r2max=d2; }
       A[i]=a;
   }

   if (r2m) { (*r2m)=r2max; }

   if (nmissed && vflag) wblog(F_L,
      "WRN %s() %d/%d values unaltered (%d)\r\\",FCT,nmissed,n,nmax);

   return err2;
};

wbstring Wb::rat2Str( 
   const char *F, int L, double d, unsigned niter,
   long pqmax,    
   double eps,    
   double eps2,   
   char vflag
){
   wbstring sout(32); 
   long p;

   if (::fabs(d-(p=::round(d)))<eps) {
      if (sout.printf(FL,"%.16g",d)>8) { 
          sout.printf(FL,"%s%ld", p || d>=0 ? "":"-", p);
      }; return sout;
   }

   if (int(niter)<1) {
      sout.printf(FL,"%.5g",d);
      return sout;
   }

   int e1, e2=-1; unsigned m,m2; long q, p2,q2;
   double r,r2, d_=d, d2, d2_, a=::fabs(d); if (a>1) { eps*=a; }

   e1=Rational(d,p,q,&r,&m,NULL,niter,pqmax,eps,eps2,vflag>1);

   if (e1 || abs(p)>999 || abs(q)>999) { 
      if ((d2_=d_*d_)>eps) { d2=d2_;
      e2=Rational(d2,p2,q2,&r2,&m2,NULL,niter,pqmax,
         d2<0.01 ? eps *a : eps,
         d2<0.01 ? eps2*a : eps2, vflag>1
      );
   }}

   if (!e1 && !e2) {
      double q1=::fabs(d*d-d2_), q2=::fabs(d2-d2_);
      if (q1>q2) { if (vflag>1) {
         PRINTF("  ==> picking sqrt() over plain (e=%.3g / %.3g)\n",
            ::sqrt(q2), ::sqrt(q1)); }
         e1=-2;
      }
   }

   if (!e1) {
      sout.printf(FL,"%ld/%ld",p,q);
      return sout;
   }
   if (e2) {
      sout.printf(FL,"%.5g",d);
      return sout;
   }

   wbstring sx(24); 

   char sgn[]="-"; if (d>=0) { sgn[0]=0; }
   sout.data[0]=0;

   if (q2==1) {
      sout.printf(FL,"%s%s",  sgn,Wb::surdStrf(sx,"%ld",p2));
   }
   else if (p2==1) {
      sout.printf(FL,"%s1/%s",sgn,Wb::surdStrf(sx,"%ld",q2));
   }
   else { r=::sqrt(double(q2));
      if (::fabs(r-::round(r))<eps) {
         sout.printf(FL,"%s%s/%g",sgn,Wb::surdStrf(sx,"%ld",p2),r);
      }
   }
   if (!sout.data[0]) { r=::sqrt(double(p2));
      if (::fabs(r-::round(r))<eps)
           { sout.printf(FL,"%s%g/%s",sgn,r,Wb::surdStrf(sx,"%ld",q2)); }
      else { sout.printf(FL,"%s%s",sgn,Wb::surdStrf(sx,"(%ld/%ld)",p2,q2)); }
   }
   return sout;
};

template<>
double Wb::SkipZeros(double *d, unsigned n, double eps){

   unsigned i=0; double s2=0;
   for (; i<n; ++i) {
       if (::fabs(d[i])<eps) { s2+=(d[i]*d[i]); d[i]=0; }
   }

   return s2;
};

template<>
double Wb::SkipZeros(long double *d, unsigned n, double eps){

   unsigned i=0; double s2=0;
   for (; i<n; ++i) {
       if (::fabs(d[i])<eps) { s2+=(d[i]*d[i]); d[i]=0; }
   }

   return s2;
};

template<>
double Wb::SkipZeros(wbcomplex *z, unsigned n, double eps){

   unsigned i=0; double s2=0;
   for (; i<n; ++i) {
       if (::fabs(z[i].r)<eps) { s2+=(z[i].r*z[i].r); z[i].r=0; }
       if (::fabs(z[i].i)<eps) { s2+=(z[i].i*z[i].i); z[i].i=0; }
   }

   return s2;
};

inline wbstring Wb::size2Str(
   double n,   
   char p  
){
   if (!n) { return "0"; } 

   unsigned l=0, m=16; size_t x;
   char s[m], ustr[8], fmt[m];

   x=size_t(1)<<30; if (n>x) { strcpy(ustr,"GB"); } else {
   x=size_t(1)<<20; if (n>x) { strcpy(ustr,"MB"); } else {
   x=size_t(1)<<10; if (n>x) { strcpy(ustr,"kB"); } else {
   x=size_t(1);                strcpy(ustr,"bytes"); }}}

   if (p<1 || p>9) {
      if (p=='c') { p=2; } 
      else p=3;
   }

   if (p>5) { snprintf(fmt,16,"%%.%dg %%s",p); }
   else { snprintf(fmt,16,"%%.%dg%%1.1s",p); } 
   l=snprintf(s,m,fmt,n/x,ustr);

   if (l>=m) wblog(FL,
      "ERR %s() string out of bounds (%d/%d)",FCT,l,m);
   return s;
};

inline const char* Wb::filename(const char *s) { 

   if (s==NULL) { wblog(FL,"WRN %s() got NULL string",FCT); }
   else {
      for (int i=strlen(s)-1; i>=0; --i) {
      if (s[i]=='/' || s[i]=='\\') return (s+i+1); }
   }
   return s;
};

size_t Wb::getFileSize(const char* f) {
   struct stat S;
   if (!f || !f[0] || stat(f, &S)!=0) { return 0; }
   return S.st_size;
};

#ifdef __APPLE__

double Wb::GET_mtime(const char *F, int L, const char *f, const struct stat &S) {

   if (S.st_mtimespec.tv_nsec>1E9) wblog(F,L,
      "WRN %s%s got mtime %ld @ %.3g", f?f:"", f?"()":"",
      S.st_mtimespec.tv_sec, 1E-9*S.st_mtimespec.tv_nsec
   );
   return (S.st_mtimespec.tv_sec + 1E-9*double(S.st_mtimespec.tv_nsec));
};

#else

double Wb::GET_mtime(const char *F, int L, const char *f, const struct stat &S) {

   if (S.st_mtim.tv_sec!=S.st_mtime || S.st_mtim.tv_nsec>1E9)
      wblog(F,L,"WRN %s() got difference in mtime (%ld @ %.3g)",
      f?f:"",S.st_mtim.tv_sec, S.st_mtime, 1E-9*S.st_mtim.tv_nsec
   );
   return S.st_mtim.tv_sec + 1E-9*double(S.st_mtim.tv_nsec);
};

#endif

#ifdef __APPLE__

double Wb::GET_atime(const char *F, int L, const char *f, const struct stat &S) {

   if (S.st_atimespec.tv_nsec>1E9) wblog(F,L,
      "WRN %s%s got atime %ld @ %.3g", f?f:"", f?"()":"",
      S.st_atimespec.tv_sec, 1E-9*S.st_atimespec.tv_nsec
   );
   return S.st_atimespec.tv_sec + 1E-9*double(S.st_atimespec.tv_nsec);
};

#else

double Wb::GET_atime(const char *F, int L, const char *f, const struct stat &S) {

   if (S.st_atim.tv_sec!=S.st_atime || S.st_atim.tv_nsec>1E9)
      wblog(F,L,"WRN %s() got difference in atime (%ld @ %.3g)",
      f?f:"",S.st_atim.tv_sec, S.st_atime, 1E-9*S.st_atim.tv_nsec
   );
   return S.st_atim.tv_sec + 1E-9*double(S.st_atim.tv_nsec);
};

#endif

#ifdef __APPLE__

double Wb::GET_ctime(const char *F, int L, const char *f, const struct stat &S) {

   if (S.st_ctimespec.tv_nsec>1E9) wblog(F,L,
      "WRN %s%s got ctime %ld @ %.3g", f?f:"", f?"()":"",
      S.st_ctimespec.tv_sec, 1E-9*S.st_ctimespec.tv_nsec
   );
   return S.st_ctimespec.tv_sec + 1E-9*double(S.st_ctimespec.tv_nsec);
};

#else

double Wb::GET_ctime(const char *F, int L, const char *f, const struct stat &S) {

   if (S.st_ctim.tv_sec!=S.st_ctime || S.st_ctim.tv_nsec>1E9)
      wblog(F,L,"WRN %s() got difference in ctime (%ld @ %.3g)",
      f?f:"",S.st_ctim.tv_sec, S.st_ctime, 1E-9*S.st_ctim.tv_nsec
   );
   return S.st_ctim.tv_sec + 1E-9*double(S.st_ctim.tv_nsec);
};

#endif

double Wb::getFileTime(const char* f, char w) {

   struct stat S;
   if (!f || !f[0] || stat(f, &S)!=0) { return 0; }

   switch (w) {
      case 'm': return GET_mtime(FLF,S);
      case 'a': return GET_atime(FLF,S);
      case 'c': return GET_ctime(FLF,S);

      case 'M': return S.st_mtime;
      case 'A': return S.st_atime;
      case 'C': return S.st_ctime;

      default: wblog(FL,"ERR %s() invalid w=%c<%d>",FCT,w,w);
   }

   return 0;
};

bool Wb::fexist(const char *f, char type) {

   if (!f || !f[0]) wblog(FL,"ERR %s() invalid (empty) file name",FCT);

   struct stat S;
   if (stat(f, &S)!=0) { return 0; }

   if (type=='f') return (S.st_mode & S_IFREG); else
   if (type=='d') return (S.st_mode & S_IFDIR); else
   return 1;
};

bool Wb::isFile(const char *f) { 
   if (f && f[0]) {
      struct stat S;
      if (stat(f,&S)==0) { return (S.st_mode & S_IFREG); }
   }
   return 0;
};

bool Wb::isDir(const char *f) { 
   if (f && f[0]) {
      struct stat S;
      if (stat(f,&S)==0) { return (S.st_mode & S_IFDIR); }
   }
   return 0;
};

wbstring Wb::hostname(unsigned len) {

   wbstring s(MAX(16U,len)); 
   int e=gethostname(s.data,s.len);

   if (e) { unsigned i=0;
      if (errno==ENAMETOOLONG) {
         for (; i<s.len; ++i) {
            if (s.data[i]=='.') { s.data[i]=0; break; }
         }
      }
      else {
         for (; i<s.len; ++i) {
            if (!isprint(s.data[i])) { s.data[i]=0; break; }
         }
         if (i<4) { strcpy(s.data,"(host!?)"); } else
         if (i>=s.len) {
            s.data[s.len-1]=0; 
         }
      }
      wblog(FL,"WRN %s() received e=%d (%s)",FCT,errno,s.data);
   }
   return s;
};

wbstring Wb::hostid(char pflag) {

   wbstring hid=hostname(pflag ? 24 : 16); 

   unsigned i=0, k=0, n=hid.len; char *s=hid.data;
   if (!n || !s[0] || !isalpha(s[0])) {
      wblog(FL,"WRN %s() received hostname `%s'",FCT,s);
      return hid;
   }

   for (; i<n; ++i) {
      if (s[i]=='-') { k=i+1; } else
      if (s[i]=='.') { if (i>4 && k+3<i) { s[i]=0; n=i; break; }} else
      if (!s[i]) { n=i; break; }
   }
   if (i<n) { n=i; } else
   if (i==hid.len) { s[i-1]=0;  
      wblog(FL,"WRN %s() received `%s'",FCT,s);
      return hid;
   }
   if (k && k+1<n) {
      for (i=k; i<n; ++i) { s[i-k]=s[i]; }
      s[i-k]=0; i-=k;
   }

   if (pflag && i+2<hid.len) { char fmt[]="/%d";
      if (pflag=='x' || pflag=='X') { fmt[2]=pflag; } 
      snprintf(s+i,hid.len-i,fmt,getpid());
   }
   return hid;
};

template <class TI> inline
TI Wb::sub2ind(const TI *s, const TI *I, unsigned n){
   if (n) {
      unsigned j=n-1; TI k=I[j];
      for (--j; j<n; --j) k=k*s[j]+I[j];
      return k;
   }
   else { wblog(FL,"ERR %s() got null set !?", FCT); return -1; }
};

template <class T1, class T2> inline
void Wb::ind2sub(T1 k, const T1 *s, T2 *I, unsigned n) {

   unsigned x,i=0;
   for (; i<n; ++i) { 
      if (s[i]) { x=k/s[i]; I[i]=k-x*s[i]; k=x; }
      else {
         wblog(FL,"ERR %s() got null-size %s",FCT,
         SSTR(wbvector<T1>(n,s,'r')));
      }
   }
   if (k) { wblog(FL,
      "ERR %s() index out of bounds\n(%d => [%s; %s])",FCT,
      (wbvector<T2>(n,I,'r')+1).toStr().data,
      SSTR(wbvector<T1>(n,s,'r')));
   }
};

int Wb::atoi(const char *s, int &iout,
   unsigned n, 
   char white  
){
   if (!s || !s[0]) { return -1; }
   unsigned i=0; if (int(n)<=0) n=9;

   if (white) {
      for (; s[i]; ++i) { if (s[i]!=' ') break; }
      if (i) {
         if (!s[i]) { return -2; }
         s+=i; i=0;
      }
   }

   if (white<0) {
      if (s[0]=='+' || s[0]=='-') { if (!s[++i]) return -3; }
   }

   for (; s[i] && i<n; ++i) {
      if (!isdigit(s[i])) return -4;
   }
   if (i==n && isdigit(s[i])) return -5;

   iout=::atoi(s);
   return i;
};

int Wb::isUIntString(const char *s, char white) {
   if (!s) { return -1; }

   if (white) {
      unsigned i=0, n=0, m=0; 

      for (; s[i]; ++i) {
         if (isdigit(s[i])) ++m; else
         if (isspace(s[i])) { if (m) { ++n; m=0; }} 
         else return -(100+n);
      }
      return (m ? n+1 : n);
   }
   else {
      unsigned i=0, n=0;

      while (isspace(*s)) ++s;      
      while (isdigit(s[n])) ++n;
      for (i=n; isspace(s[i]); ++i) {}; 

      if (!s[i])
           return n; 
      else return -(10+n);
   }
};

int Wb::dstrlen_utf8(const char *s) {
   int q=0; 
   if (!s || !s[0]) { return q; }

   for (int i=0, l=strlen(s); i<l && s[i]; ++i) {
      if (s[i]<0 || s[i]>127) {
         if ((s[i] & 0xE0) == 0xC0) { q+=1; } else 
         if ((s[i] & 0xF0) == 0xE0) { q+=2; } else 
         if ((s[i] & 0xF8) == 0xF0) { q+=3; } else 
         { q=0; break; } 
      }
   }
   return q;
};

template <class T>
wbstring Wb::bits(const T &x, char compact) {
   wbstring s_; 
   int i, j=0, n=sizeof(T), l=9*n; { s_.init(l); if (!l) return s_; }
   char *s; unsigned char c;
   const char *sx = (char*)(&x);

   if (compact) {
      for (i=0; i<n; ++i) { if (sx[i]) j=i; }
      n=j+1; l=9*n;
   }
   s=s_.data+l-1; s[0]=0;

   for (i=0; i<n; ++i, s-=9) { c=sx[i];
      for (j=1; j<=8; ++j, c>>=1) { s[-j]=(c&1 ? '1' : '-'); }
      if (i+1<n) s[-j]=' ';
   }

   return s_;
};

void Wb::IOstat::print(const char *istr,
   const char *hstr, const char *xstr1, const char *xstr2
){
   if (hstr && hstr[0]) printf("\n  I/O stats for %-9s "
      "#read%-6s size_read   #write%-5s size_write\n",
      hstr, xstr1?xstr1:"", xstr2?xstr2:"");

   char s[64];
   snprintf(s,   30,"%4ld",aux1);
   snprintf(s+30,30,"%4ld",aux2);

   printf("     %-19s %6d %4s %10s   %6d %4s %10s\n",istr ? istr:"",
      nread,  aux1 ? s    : "", size2Str(sread ).data,
      nwrite, aux2 ? s+30 : "", size2Str(swrite).data
   );
};

bool Wb::is_finite(const wbcomplex *x, size_t n) {
   for (size_t i=0; i<n; ++i) { if (!x[i].isfinite()) { return 0; }}
   return 1;
};

wbstring cpu_time::toStr(char flag) {

   size_t l=0, n=32; char tstr[n];

   cpu_time dt=since();
   double t=(flag=='c' ? dt.tcpu : dt.tsys);

   if (t<0) {
      wblog(FL,"WRN %s() got negative time (%ld)",FCT,long(t)); t=-t; }

   if (t<60) {
      l=snprintf(tstr,n,"%.4g sec",t); }
   else if (t<1E14) {
      int s,m,h,d; long x=t,
      fac=24*3600; d=t/fac; x-=d*fac;
      fac=3600;    h=x/fac; x-=h*fac;
      fac=60;      m=x/fac; x-=m*fac; s=x; x-=s;

      if (d)
           l=snprintf(tstr,n,"%d-%02d:%02d:%02d",d,h,m,s);
      else if (h || m>9)
           l=snprintf(tstr,n,   "%02d:%02d:%02d",  h,m,s);
      else l=snprintf(tstr,n,        "%02d:%02d",    m,s);

      if (t>2E9) l=-l; 
   }
   else {
      l=snprintf(tstr,n,"(WRN %.3g yrs (%ld sec))",
      t/(3600*24*365), long(t));
      l=-l; 
   }

   if (l>=n) {
     long tl=(flag=='c' ? dt.tcpu : dt.tsys);
     double yrs=double(tl)/(3600*24*365);
     wblog(FL, 
        "WRN cpu_time() %s `%s'\n'%c' dt=%ld => %.12g yrs (%d) !?\n"
        "%15.12g %15.12g %% now [CPU/SYS]\n"
        "%15.12g %15.12g %% then",
        int(l)<0 ? "possibly invalid CPU time":"string out of bounds\n",
        tstr, flag,tl,yrs,l, dt.time('c'), dt.time('s'), tcpu, tsys
     );
   }

   return tstr;
};

#ifdef __APPLE__

double Mac::getMemSize(const char *tag) {

   mach_msg_type_number_t count = HOST_VM_INFO_COUNT;
   vm_statistics_data_t S;

   if (host_statistics(  
       mach_host_self(), 
       HOST_VM_INFO,     
       (host_info_t)&S,  
       &count            
   )) wblog(FL,"ERR %s() %s",FCT,strerror(errno));

   long n = getpagesize(),
        sused = S.wire_count + S.active_count,
        sfree = S.inactive_count + S.free_count;

   if (!tag || !tag[0]) wblog(FL,"ERR %s() got empty tag !?",FCT);

   if (!strcmp(tag,"total")) { return n*(sused+sfree); }
   if (!strcmp(tag,"free" )) { return n*(sfree); }
   if (!strcmp(tag,"usage")) {
      return 100*double(sused)/(sfree+sused); 
   }

   wblog(FL,"ERR %s() tag '%s' not yet defined",FCT,tag);
   return 0;
};

size_t Mac::getProcMemSize(pid_t p, char res) {

   struct task_basic_info S;
   mach_msg_type_number_t len = TASK_BASIC_INFO_COUNT;

   task_t tid = MACH_PORT_NULL;
   task_for_pid(current_task(), p ? p : getpid(), &tid);

   if (task_info( 
      tid, TASK_BASIC_INFO, (task_info_t)&S, &len)
   ){ wblog(FL,"ERR %s() %s",FCT,strerror(errno)); }

   return (res ? S.resident_size : S.virtual_size);
};

size_t Wb::getProcSize(const char *tag, const pid_t &p) {
   return Mac::getProcMemSize(p);
};

#else 

size_t Wb::getProcSize(const char *tag, const pid_t &p) {
   const int n=32; char f[n];
   int i=snprintf(f,n,"/proc/%d/status",p);
   if (i>=n) wblog(FL,
      "ERR %s() string out of bounds (%s; %d)",FCT,f,n);
   return Wb::getProcSize(tag,f);
};

size_t Wb::getProcSize(const char *tag, const char *f) {

   char *s, *line=NULL, u[8];
   size_t val; unsigned i=0, m=0; size_t len=0; ssize_t n;

   static unsigned ncall=0;

   FILE *fid; ++ncall;
   if (!f || !f[0] || !(fid=fopen(f,"r"))) {
      wblog(FL,"WRN failed to access %s '%s' (%d)",f?f:"",tag?tag:"",ncall);
      return 0;
   }

   while ((n=getline(&line,&len,fid))>=0) { ++i;
      if (n) line[n-1]=0; 
      if ((s=strstr(line,tag))) { u[0]=0; ++m;
         sscanf(s,"%*s %lu %s",&val,u);
         if (strlen(u)==2 && u[1]=='B') {
             if (u[0]=='k' || u[0]=='K') { val*=(1<<10); } else
             if (u[0]=='M')              { val*=(1<<20); } else
             if (u[0]=='G')              { val*=(1<<30); }
             else wblog(FL,"WRN %s() got unit `%s'",FCT,u);
         }
         else wblog(FL,
            "WRN %s() invalid unit >%s< (%d)",FCT,u,strlen(u));
         break;
      }
   }

   if (line) { free(line); line=NULL; }
   fclose(fid); fid=NULL;

   if (!m) wblog(FL,"WRN %s() '%s' not found in %s",FCT,tag,f);

   return val;
};

long Wb::getCpuInfo(const char *tag, const char *f) {

   char *s, *line=NULL;
   int m=0, val=0, i=0, vmax=-1; size_t len=0; ssize_t n;

   FILE *fid=fopen(f,"r");
   if (fid==NULL) {
      wblog(FL,"WRN failed to access %s '%s'",f?f:"",tag?tag:"");
      return 0;
   }

   if (!strcmp(tag,"processor")) {
      while ((n=getline(&line,&len,fid))>=0) { ++i;
         if (n) line[n-1]=0; 
         if ((s=strstr(line,tag))) { ++m;
            sscanf(s,"%*s : %d",&val); if (vmax<val) vmax=val;
         }
      }

      if (vmax+1!=m) wblog(FL,
         "WRN %s() inconsistent ncores = %d/%d !?",FCT,vmax,m);
      val=m;
   }
   else wblog(FL,"ERR %s() invalid tag %s -> %s",FCT,f,tag);

   if (line) { free(line); line=NULL; }
   fclose(fid); fid=NULL;

   if (!m) wblog(FL,"WRN %s() '%s' not found in %s",FCT,tag,f);

   return val;
};

#endif 

#ifdef __APPLE__

size_t wbsys::getMemTot()  { return Mac::getMemSize("total"); };
size_t wbsys::getMemFree() { return Mac::getMemSize("free" ); };

int wbsys::getNumCores() {
   int n=-1; 

   size_t l=4; 
   int e=sysctlbyname("hw.physicalcpu",&n,&l,NULL,0);

   if (e || n<1) { wblog(FL,
      "ERR %s() sysctl returned e=%d (n=%d)\n%s",
      FCT,e,n, e? strerror(errno):"''");
   }
   return n;
};

int wbsys::getCacheLineSize() {
   size_t n=-1; 
   size_t l=sizeof(n); 
   int e=sysctlbyname("hw.cachelinesize",&n,&l,NULL,0);

   if (e || n<1) { wblog(FL,
      "ERR %s() sysctl returned e=%d (n=%d)\n%s",
      FCT,e,n, e? strerror(errno):"''");
   }
   return n;
};

#else

size_t wbsys::getMemTot()   { return Wb::getProcSize("MemTotal" ); };
size_t wbsys::getMemFree()  { return Wb::getProcSize("MemFree"  ); };
size_t wbsys::getSwapTot()  { return Wb::getProcSize("SwapTotal"); };
size_t wbsys::getSwapFree() { return Wb::getProcSize("SwapFree" ); };
int    wbsys::getNumCores() { return Wb::getCpuInfo("processor" ); };

int wbsys::getCacheLineSize() { 
   int n=sysconf(_SC_LEVEL1_DCACHE_LINESIZE); 
   if (n<8) { wblog(FL,"ERR %s() got L1d cache-linesize ls=%d !?",FCT,n); }
   return n;
}

#endif

void Wb::ResSummary(const char *F, int L, const char *istr) {

   static int do_cleanup=1;
   if (F==NULL && L) { do_cleanup=L; return; }

  #ifndef __APPLE__
   if (do_cleanup) { do_cleanup=0; wbtop PS;
      wbstring s1(PS.VmSize2Str()), s2(PS.VmPeak2Str());
      if (istr) { sprintf_str("%s() ",istr); } else { str[0]=0; }

      wblog(F_L," *  %sVMEM size: %s / %s", str, s1.data, s2.data);
      wblog(F_L," *  %scputime  : %s",str, STR2(gCPUTime,'c')); 
   }
  #else
   if (do_cleanup) { do_cleanup=0;
      if (istr) { sprintf_str("%s() ",istr); } else { str[0]=0; }
      wblog(F_L," *  %sVMEM size / cputime (skipped for OSX)", str);
   }
  #endif
};

void Wb::MemStat(const char *F, int L, char lflag) {
   wbtop().MemStat(F,L,lflag);
};

inline wbstring wbsys::MemTot2Str() {
   return Wb::size2Str(getMemTot()); };

inline wbstring wbsys::MemFree2Str() {
   return Wb::size2Str(getMemFree()); };

inline wbstring wbsys::SwapTot2Str() {
   return Wb::size2Str(getSwapTot()); };

inline wbstring wbsys::SwapFree2Str() {
   return Wb::size2Str(getSwapFree()); };

#ifdef __APPLE__

char wbsys::checkSwapSpace(const char *F, int L) {

   static double xref=0.25; 

   int e, mib[2] = { CTL_VM, VM_SWAPUSAGE };
   struct xsw_usage S; size_t l=sizeof(S);

   if ((e=sysctl(mib,2,&S,&l,NULL,0))) wblog(FL,
      "ERR %s() sysctl returned %d\n(%s)",FCT,e,strerror(errno));

   if (S.xsu_total<=0 || double(S.xsu_avail)/S.xsu_total>xref) { return 0; }

   wblog(F_L,"WRN free swap space: %s / %s @ %.3g",
      Wb::size2Str(S.xsu_avail).data,
      Wb::size2Str(S.xsu_total).data, xref); xref-=0.05;
   return 1;
};

#else

char wbsys::checkSwapSpace(const char *F, int L) {
   static double xref=0.25;
   double x=getSwapTot(); if (x<=0) return 0;
   x=getSwapFree()/x; if (x>xref) { return 0; }

   xref-=0.05; wblog(F_L,"WRN free swap space: %s / %s",
   SwapFree2Str().data, SwapTot2Str().data);
   return 1;
};
#endif

inline wbstring wbtop::VmSize2Str() const {
   return Wb::size2Str(getVmSize()); };

inline wbstring wbtop::VmPeak2Str() const {
   return Wb::size2Str(getVmPeak()); };

int wbtop::runningLarge(const char *F, int L,
   double th0, 
   double fac  
){
   static double thr=th0; 
   int r=0; size_t s0, M;

   s0=wbsys::getMemTot();
   M=getVmPeak();

   if (M<thr*s0) { wbsys::checkSwapSpace(F_L); }
   else {
      wblog(F_L,"%s() %s is using %s",FCT,myname,VmPeak2Str().data);
      thr*=fac; r=1;
   }
   doflush(); wblogBuf.flush();

   return r;
};

void wbtop::MemStat(const char *F, int L, char lflag) const {
   static size_t vm_size_max=0; 
   static size_t vm_size_last=0; 

   if (lflag=='m' || lflag=='M') {
      size_t q=getVmSize();
      if (vm_size_max*(lflag=='m' ? 1.10 : 1)<q)
           { vm_size_max=q; }
      else { lflag=0; }
   }
   else if (!lflag && getVmPeak()>0.66*sysmem) {
      size_t q=getVmSize();
      if (!vm_size_last ||
         fabs((q-vm_size_last)/double(vm_size_last))>0.05
      ){ vm_size_last=q; lflag=2; }
   }

   if (lflag) {
      double gfac=1/double(1<<30);
      wblog(F_L,"SYS %s() %.1fG (%.1f/%.0fG)",FCT,
         gfac*getVmSize(), gfac*getVmPeak(), gfac*sysmem);
      doflush();
   }
};

template <class T>
size_t Wb::findfirst_sorted(const char *F, int L,
   const T* d0, const T* dd, size_t m, size_t N, size_t M, char lex
){
   if (!N || !M || m>=M) { if (F) wblog(F,L,
      "ERR %s() got empty/invalid input (%dx%d; %d",FCT,N,M,m);
      return -1;
   }

   int q=Wb::cmpRange(dd,d0,m,lex);
   if (q>=0) {
      if (q>0) { if (F) wblog(F,L, 
         "ERR %s() empty match (<first)",FCT);
         return -1; 
      }
      return 0; 
   }

   size_t i1=0, i2=N-1, i=(i2+i1)/2;
   int q2=Wb::cmpRange(dd+i2*M,d0,m,lex);

   if (q2<0) { if (F) wblog(FL, 
      "ERR %s() empty match (>last)",FCT);
      return N; 
   }

   while (i>i1) {
      q=Wb::cmpRange(dd+i*M,d0,m,lex);
      if (q<0) { i1=i; }
      else {
         if (q>0 && q!=q2) wblog(F_L,"ERR %s() "
            "got unsorted data (upper i=%d: %d/%d) !?",FCT,i,q,q2);
         q2=q; i2=i;
      }
      i=(i1+i2)/2;
   }

   if (q2) { if (F) wblog(FL,
      "ERR %s() empty match (i=%d/%d @ %d,%d)",FCT,i,N,q,q2);
      return -1;
   }
   return i2;
};

template <class T>
size_t Wb::findlast_sorted(const char *F, int L,
   const T* d0, const T* dd, size_t m, size_t N, size_t M, char lex
){
   if (!N || !M || m>=M) { if (F) wblog(F,L,
      "ERR %s() got empty/invalid input (%dx%d; %d",FCT,N,M,m);
      return -1;
   }

   int q1=Wb::cmpRange(dd,d0,m,lex);

   if (q1>0) { if (F) wblog(F,L, 
      "ERR %s() empty match (<first)",FCT);
      return -1;
   }

   size_t i1=0, i2=N-1, i=(i2+i1)/2;
   int q=Wb::cmpRange(dd+i2*M,d0,m,lex);

   if (q<=0) {
      if (q<0) { if (F) wblog(F,L,
         "ERR %s() empty match (>last)",FCT);
         return N;
      }
      return i2;  
   }

   while (i>i1) {
      q=Wb::cmpRange(dd+i*M,d0,m,lex);
      if (q>0) { i2=i; }
      else {
         if (q<0 && q!=q1) wblog(F_L,"ERR %s() "
            "got unsorted data (lower i=%d: %d/%d) !?",FCT,i,q,q1);
         q1=q; i1=i;
      }
      i=(i1+i2)/2;
   }

   if (q1) { if (F) wblog(FL,
      "ERR %s() empty match (i=%d/%d @ %d,%d)",FCT,i,N,q1,q);
      return -1;
   }
   return i1;
};

void wb_srand() {

   struct timespec tnow; 
   clock_gettime(CLOCK_REALTIME,&tnow); 

   unsigned long low, high; 
   __asm__ __volatile__("rdtsc" : "=a" (low), "=d" (high) : );

   ::srand(::rand() ^ unsigned(getpid()) ^ high ^ low);
   ::srand(::rand() ^ tnow.tv_sec ^ tnow.tv_nsec);
};

#ifdef QS_USING_MPFR

namespace Wb {

template<> inline 
quad& rand::rand_(quad &x) { return (x.Rand()); };

template<> inline 
quad& rand::randb(quad &x) {
   x.Rand(); x*=2; x-=1; return x;
};

template<> inline 
quad& rand::rands(quad &x) {
    x.Rand(); x*=2; x-=1;
    return (x.Atanh());
};

}; 

template <> inline
Wb::quad WbUtil<Wb::quad>::eps() {
    Wb::quad x;
    return x.eps();
};

#endif

template <> inline
wbcomplex WbUtil<wbcomplex>::eps() { return DBL_EPSILON; };

#endif

