/* ---------------------------------------------------------------------
 * Project : QSpace tensor library (v4.0 pre-release)
 * Class   : QSpace NRG routines
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

#ifndef __WB_SPECTRAL_HCC__
#define __WB_SPECTRAL_HCC__

// default smoothing: norm-preserving log-gauss broadening
   #define LOG_GAUSS_BRD  2

namespace Wb {

template <class TS>
void Fourier(
   const wbvector<double> &om,
   const wbMatrix<TS> &A,
   const wbvector<double> &tt,
   wbMatrix<wbcomplex> &az,
   char domflag=0,
   double alpha=0., 
   char vflag=1
);

template <class TS>
void dispSumRule(
   const wbvector<double> &om,
   const wbMatrix<TS> &A
);

}; 

void KKimag(
   const wbMatrix<double> &Rp,
   const wbMatrix<double> &Rn,
   const wbvector<double> &OM,
   wbMatrix<double> &Ip,
   wbMatrix<double> &In
);

template <class TS> 
class Spectral {
 public:

    Spectral()
     : sigma(0.6), sigma2(0.5), alpha(0.25), Emin(1E-8), Emax(10.),
       eps(1E-6), nlog(128), l0(0), fac(1.), Delta(0.) {};

    void init() {
       Ap.init(); An.init(); Ap_buf.init(); An_buf.init();
       Ar.init(); A0.init(); mspec.init();
       istr.init();
    };

    void init(const char *istr, unsigned Nops,
       double delta=0., double T_=0, unsigned Niter=0
    );

    void initBuf() {
       if (Ap.dim1!=An.dim1) wblog(FL,
       "ERR severe dimension mismatch (%d,%d).", Ap.dim1, An.dim1);

       Ap_buf.init(Ap.dim1, Ap.dim2);
       An_buf.init(An.dim1, An.dim2);
    };

    void calcSpectralMoments(const char *F, int L, unsigned n) {
       if (n) {
          if (n>4) wblog(F,L,"WRN %s with n=%d !?",FCT,n);
          else wblog(F,L,"<i> calculate %d spectral moment(s)",n);
          mspec.init(Ap.dim1,n+1);
       }
    };

    void Add(
       const wbarray<double> &Om, const wbarray<TS> &CC, unsigned s,
       char dblock=0, 
       char mflag=0   
    );

    void Add2buf(const wbarray<double> &Om, const wbarray<TS> &CC, unsigned s);

    void crossAddBuf(const char minor=1, char vflag=1);

    void saveIter(const char *F, int L, unsigned iter);

    void getOmega(wbvector<double> &om, unsigned &linfit) const;
    void getOmega(wbvector<double> &om) const {
       unsigned linfit=0; getOmega(om,linfit);
    };

    void applyIROPfac(const char *F, int L, const wbvector<double> &fac);

    void getRawData(
       wbvector<double> &om, wbMatrix<TS> &A,
       unsigned linfit=0, const char *isbuf=""
    ) const;

    void getPARTIAL(
       wbarray<TS> &A1, wbarray<TS> &A2) const;

    void pairRawData();

    void detailedBalance(
    const double T, wbvector<double> &om, wbMatrix<TS> &aa) const;

    void dispSumRule();
    void dispSumRule(const wbvector<double> &om, const wbMatrix<TS> &A);

    void getSmoothSpec(
       wbvector<double> &om, wbMatrix<TS> &A, char vflag=1);

    void getSmoothSpec(
       wbMatrix<TS> &A, const wbvector<double> &om, char vflag=1);

    void KKreal(
       wbMatrix<wbcomplex> &GZ,
       wbvector<double> &OM,
       wbvector<double> &DOM
    ) const;

    wbMatrix<TS> Ap, An; 
    wbMatrix<TS> Ap_buf, An_buf;
    wbvector<TS> Ar;     
    wbMatrix<TS> A0;     
    wbMatrix<TS> mspec;  

    wbarray<TS> AP, AN;

    wbstring istr;
    double Temp; 

    double sigma, sigma2, alpha, Emin, Emax, eps;
    unsigned nlog;

 protected:
 private:

    double l0, fac, Delta;

    double OM2IDX(const double &om) const {
       return ((log(fabs(om))-l0)*fac);
    }

    inline unsigned om2idx(double om) const { 
       if (Delta==0) { if (om==0.) return 0; }
       else {
          if (om> Delta) om-=Delta; else
          if (om<-Delta) om+=Delta; else return 0;
       }

       int k = (int)floor(OM2IDX(om));
       if (k<1) k=1; else
       if (k>=(int)Ap.dim2) { k=Ap.dim2-1; }

       return (unsigned)k;
    }

    inline double idx2om(const unsigned &k) const {
       return exp( (double(k)+0.5)/fac + l0 );
    }

    int crossAddBuf_aux(
       TS *a, const TS *b, const unsigned M,
       const char minor=1, char vflag=1
    );

    TS foldLogGaussian(
       const double om, const TS *aa, const unsigned n
    );

    TS foldGaussian(
       const double om, const double b, 
       const double *oo, const TS *aa, const unsigned n
    );
};

template <class TS> 
class RSpecData {   
  public:
     RSpecData() { };

     RSpecData& Append(const double *op, const TS *cp, unsigned n) {
        om.Append(n,op);
        cc.Append(n,cp); return *this;
     };

     mxArray* mxCreateStruct(unsigned m, unsigned n) const {
        const char *flds[]={"om","cc"};
        return mxCreateStructMatrix(m,n,2,flds);
     };
     void add2MxStruct(mxArray *S, unsigned i, char tflag=0) const {
         mxSetFieldByNumber(S,i,0,om.toMx(tflag));
         mxSetFieldByNumber(S,i,1,cc.toMx(tflag));
     };

     wbvector<double> om;
     wbvector<TS> cc;

 protected:
 private:

};

template <class TS> 
class TDSpectral {  

 public:

    TDSpectral() : raw(0), l0(0), fac(1.), check_avg_om(0) {};

    void init(const char *F, int L,
       const char *istr0, const wbvector<double> &t0,
       unsigned Nnrg, unsigned Nop,
       double Emin=1E-8, double Emax=10., unsigned nl=512,
       char vflag=1
    );

    void checkAvgOm() {
       check_avg_om=1;
       AW.init(2,AP.dim(2),AP.dim(3));
    };
    wbarray<TS> getAvgOm();

    void Add(
       const wbarray<double> &Om, const wbarray<TS> &CC,
       unsigned iter, unsigned s
    );

    void addPartialFourier();

    void getOmega(wbvector<double> &om, unsigned &linfit) const;
    void getOmega(wbvector<double> &om) const {
       unsigned linfit=0; getOmega(om,linfit);
    }

    void applyIROPfac(
       const char *F, int L, const wbvector<double> &fac, char rixs=0);

    void getSpecData(
       wbvector<double> &om, wbarray<TS> &A,
       unsigned linfit=0
    ) const;

    void getSpecData(
       wbvector<double> &om, wbMatrix<TS> &A,
       unsigned linfit=0
    );

    void mergeNRGData();

    void getSmoothSpec_t(double sigma, double eps,
       char bswitch = LOG_GAUSS_BRD, char vflag=1);

    void Fourier(
       wbvector<double> &om, wbMatrix<wbcomplex> &ST,
       double alpha, double Lambda, char vflag=1
    );

    void Fourier(
       wbvector<double> &om, wbMatrix<TS> &A,
       const wbvector<double> &tt,
       wbMatrix<wbcomplex> &at,
       double alpha=0., char vflag=1
    );

    void dispSumRule();

    wbarray<TS> AP,AN;  
    wbMatrix<TS> Ap,An; 
    wbMatrix<wbcomplex> At;
    wbvector<double> tt;    
    wbstring istr;

    mxArray *toMx() const;
    int init(
       const mxArray *S,
       double &emin, double &emax, unsigned &nlog
    );

    unsigned raw;
    wbvector< RSpecData<TS> > RAW;

 protected:
 private:

    double l0,fac; 
    unsigned N;

    char check_avg_om;
    wbarray<TS> AW;

    double OM2IDX(const double &om) const {
       return ((log(fabs(om))-l0)*fac);
    }

    inline unsigned om2idx(const double &om) const {
       if (om==0.) return 0;

       int k = (int)floor(OM2IDX(om));
       if (k<1) k=1; else if (k>=(int)N) k=N-1;

       return (unsigned)k;
    }

    inline double idx2om(const unsigned &k) const {
       return exp( ((double)k+0.5)/fac + l0 );
    }

    void initLG(
       wbMatrix<double> &wlg, unsigned n, double sigma,
       char bswitch = LOG_GAUSS_BRD, char vflag=1
    );

    void foldLogGaussian(
       const wbMatrix<double> &wlg,
       TS *aa, TS a0, unsigned i0
    );
};

template <class TS>
void Spectral<TS>::init (
   const char *istr0,
   unsigned Nops, 
   double d,      
   double T,      
   unsigned Niter 
){
   unsigned Nom;
   double l2;

   if (Emin<=0. || Emax<=0. || Emin>=0.1*Emax || nlog==0) wblog(FL,
      "ERR %s() invalid range E = [%g .. %g @ %d]%s", myname, Emin, Emax, nlog,
      (Emin>0 && Emax>0 && Emin<Emax) ? " (require Emin<Emax/10)":"");
   if (d<0) wblog(FL,"ERR invalid Delta (%g neg!)",d);

   if (eps<Emin) wblog(FL,"NB! eps=%g < min(|om|)=%g", eps, Emin);
   if (alpha!=0.25) { 
      if (alpha==0.5) wblog(FL,
         "WRN %s() using symmetric log-Gaussian broadening",FCT);
      else wblog(FL,"WRN %s() using "
         "log-Gaussian with shift alpha=%g !?",FCT,alpha);
   }

   Delta=d;

   l0 =log(Emin);
   l2 =log(Emax);
   fac=log10(M_E)*double(nlog); 
   Nom=(unsigned)((l2-l0)*fac)+1;

   Ap.init(Nops,Nom);
   An.init(Nops,Nom); Ar.init(Nops);

   if (int(Niter)>0) {
      AP.init(Nom,Nops,Niter); 
      AN.init(Nom,Nops,Niter);
   }

   istr=istr0; Temp=T;
};

template <class TS>
void Spectral<TS>::Add(
   const wbarray<double> &Om,
   const wbarray<TS> &CC,
   unsigned s,
   char dblock, 
   char mflag   
){
   unsigned i,k, n=Om.SIZE.prod(); 
   double x, emin=0.1*idx2om(0);
   TS q, r=0;

   double bfac=(!mflag && Temp>0 ? 10/Temp : 0); 

   double zfac=(Temp>0 ? 1/Temp : 0);

   if (!Om.sameSize(CC)) {
      Om.info("Om"); CC.info("CC");
      wblog(FL,"ERR size mismatch between Om and CC!");
   }
   if (dblock) {
      if (Om.SIZE.len!=2 || Om.SIZE[0]!=Om.SIZE[1]) wblog(FL,
         "ERR %s() quadratic block expected (%s)",
         FCT,Om.sizeStr().data);
   }

   if (s>=Ap.dim1 || s>=Ar.len) wblog(FL,
      "ERR index out of bounds (%d/%d;%d)",s, Ap.dim1,Ar.len);

   for (i=0; i<n; ++i) { x=Om.data[i]; k=om2idx(x);
       if (x<0. || (x<=emin && s%2)) 
            An(s,k)+=CC.data[i];
       else Ap(s,k)+=CC.data[i];

       if (x!=0) {
          q=CC.data[i]/x;
          if (bfac>0) { x=bfac*fabs(x); if (x<10) {
             double w=exp(-x*x);  
             q*=(1-w);
          }}
          r+=q;
       }
       else if (mflag && !(s%2)) { 
          r+=zfac*CC.data[i]; 
       }
   }

   Ar[s]+=r;

   if (mspec.dim2) {
      double dbl; unsigned j;

      if (mspec.dim1!=Ap.dim1) {
         wblog(FL, 
           "WRN mspec must be initialized AFTER ASpec!\n"
           "resize (%d->%d x %d)", mspec.dim1, Ap.dim1, mspec.dim2
         );
         mspec.Resize(Ap.dim1,mspec.dim2);
      }

      for (i=0; i<n; i++)
      for (dbl=CC.data[i], j=0; j<mspec.dim2; j++) {
         if (j) dbl*=Om.data[i];
         mspec(s,j)+=dbl;
      }
   }
};

template <class TS>
void Spectral<TS>::Add2buf(
   const wbarray<double> &Om, const wbarray<TS> &CC,
   unsigned s
){
   unsigned i,k,n=Om.SIZE.prod();

   if (!Om.sameSize(CC)) {
      Om.info("Om"); CC.info("CC");
      wblog(FL,"ERR Size mismatch between Om and CC!");
   }

   if (s>Ap_buf.dim1)
   wblog(FL, "ERR Index out of bounds (%d,%d)",s, Ap_buf.dim1);

   for (i=0; i<n; i++) {
       k=om2idx(Om.data[i]);
       if (Om.data[i]<0.)
            An_buf(s,k)+=CC.data[i];
       else Ap_buf(s,k)+=CC.data[i];
   }
};

template <class TS>
void Spectral<TS>::crossAddBuf(
   const char minor, 
   char vflag
){
   unsigned s, Nom=Ap.dim2;
   int i;

   if (!Ap.sameSize(Ap_buf) || !An.sameSize(An_buf) ||
       !Ap.sameSize(An)
   ){
       Ap.info("Ap"); Ap_buf.info("Ap_buf");
       An.info("An"); An_buf.info("An_buf");
       wblog(FL,"ERR Spectral - severe size mismatch.");
   }

   for (s=0; s<Ap.dim1; s++) { i=0;

      i+=crossAddBuf_aux(Ap.rec(s), Ap_buf.rec(s), Nom,minor,vflag);
      i+=crossAddBuf_aux(An.rec(s), An_buf.rec(s), Nom,minor,vflag);

      if (i>1) {
         wblog(FL,"WRN empty buffer (%d/%d) !?",s,Ap.dim1);
         continue;
      }
   }

};

template <class TS>
int Spectral<TS>::crossAddBuf_aux(
   TS *a,        
   const TS *b,  
   const unsigned Nom,
   const char minor, 
   char vflag
){
   unsigned i, i1, i2, amin, amax, bmin, bmax; 
   double t, tfac;

   if (minor<1 || minor>2) wblog(FL,
      "ERR invalid minor flag %d", minor);

   for (i=0;     i<Nom;   ++i) { if (a[i]!=0.) break; }; amin=i;
   for (i=Nom-1; i>=amin; --i) { if (a[i]!=0.) break; }; amax=i;
   for (i=0;     i<Nom;   ++i) { if (b[i]!=0.) break; }; bmin=i;
   for (i=Nom-1; i>=bmin; --i) { if (b[i]!=0.) break; }; bmax=i;

   if (bmin>bmax) {
      return vflag;
   }

   if (amin>amax) {
      for (i=bmin; i<=bmax; i++) a[i]=b[i];
      return 0;
   }

   i1=MAX(amin,bmin); i2=MIN(amax,bmax);

   if (i1==i2) { a[i1]=0.5*(a[i1]+b[i2]); return 0; }

   tfac = 1. / double(i2-i1);

   for (i=0; i<Nom; ++i) {
      if (i<i1 || i>i2) a[i]+=b[i];
      else {
         t = tfac * double(i-i1); 
         if (minor==1)
              a[i] = (1.-t)*a[i] + t*b[i];
         else a[i] = t*a[i] + (1.-t)*b[i];
      }
   }

   return 0;
};

template <class TS>
void Spectral<TS>::saveIter(const char *F, int L, unsigned iter){

   unsigned n=Ap.numel();

   if (AP.SIZE.len!=3 || AP.SIZE!=AN.SIZE || n!=AP.SIZE[0]*AP.SIZE[1])
   wblog(F_L,
     "ERR %s() invalid initialization of AP/AN [%s; %s; %dx%d]",
      FCT,AP.sizeStr().data, AN.sizeStr().data, Ap.dim1, Ap.dim2
   );
   if (iter>=AP.SIZE[2]) wblog(F_L,
     "ERR %s() iter out of bounds AP/AN [%s; %d/%d]",
      FCT,AP.sizeStr().data, iter+1, AP.SIZE[2]
   );

   memcpy(AP.data+iter*n, Ap.data, n*sizeof(Ap.data[0]));
   memcpy(AN.data+iter*n, An.data, n*sizeof(An.data[0]));
};

template <class TS>
void Spectral<TS>::getOmega(wbvector<double> &om, unsigned &linfit) const {

   unsigned i, Nom=Ap.dim2;
   double dbl;

   if (linfit) linfit=unsigned(fac);
   if (linfit>Nom) wblog(FL,"ERR Nom<linfit (%d,%d) !?",Nom,linfit);

   om.init(2*Nom);

   dbl = idx2om(linfit) / (linfit+0.5);

   for (i=0; i<linfit; ++i) om[Nom-1-i] = -(om[Nom+i] = dbl*(i+0.5));
   for (   ; i<Nom;    ++i) om[Nom-1-i] = -(om[Nom+i] = idx2om(i)  );
}

template <class TS>
void Spectral<TS>::applyIROPfac(
   const char *F, int L, const wbvector<double> &fac) {

   unsigned i1,i2,s,m, n=fac.len;
   double x;

   if (2*n!=Ap.dim1 || 2*n!=An.dim1 || 2*n!=Ar.len) wblog(FL,
      "ERR %s() size mismatch (%d,%d,%d/2*%d)",
       FCT,Ap.dim1,An.dim1,Ar.len,n);

   for (s=0; s<n; ++s) { if (fac[s]!=1) {
      if (fac[s]<1 || fac[s]!=round(fac[s])) wblog(FL,
         "ERR %s() got fac[%d]=%g !?",FCT,s+1,fac[s]);
      wblog(FL," *  applying IROP factor 1/%g to op(%d)",fac[s],s+1);

      x=1/fac[s]; i1=2*s; i2=2*s+1;
      Ap.rowTimes(i1,x); An.rowTimes(i1,x); Ar[i1]*=x;
      Ap.rowTimes(i2,x); An.rowTimes(i2,x); Ar[i2]*=x;
   }}

   if (!Ap_buf.isEmpty() || !An_buf.isEmpty()) {
      if (!Ap_buf.sameSize(Ap) || !An_buf.sameSize(Ap)) wblog(FL,
         "ERR %s() severe size mismatch %dx%d <> %dx%d <> %dx%d !?", FCT,
          Ap_buf.dim1,Ap_buf.dim2,An_buf.dim1,An_buf.dim2,Ap.dim1,Ap.dim2
      );
      for (s=0; s<n; ++s) { if (fac[s]>1) { x=1/fac[s]; i1=2*s; i2=2*s+1;
         Ap_buf.rowTimes(i1,x); An_buf.rowTimes(i1,x);
         Ap_buf.rowTimes(i2,x); An_buf.rowTimes(i2,x);
      }}
   }

   if (!AP.isEmpty() || !AN.isEmpty()) {
      if (AP.SIZE.len!=3 || !AP.sameSize(AN) || AP.SIZE[1]!=2*n)
      wblog(FL,
         "ERR %s() size mismatch (%s; %s; %dx%d) !?", FCT,
          AP.sizeStr().data, AN.sizeStr().data, Ap.dim1,Ap.dim2
      );

      TS *ap=AP.data, *an=AN.data;
      unsigned n2=2*AP.SIZE[0]; m=n*AP.SIZE[2];

      for (i1=0; i1<m; ++i1, ap+=n2, an+=n2) { s=i1%n;
         if (fac[s]>1) { x=1/fac[s];
            Wb::timesRange(ap, x, n2);
            Wb::timesRange(an, x, n2);
         }
      }
   }

   if (!A0.isEmpty() || !mspec.isEmpty()) { wblog(FL,
        "WRN IROPfac yet not applied to A0 and mspec (%d,%d)",
         A0.isEmpty(), mspec.isEmpty()); 
      wblog(FL," *  A0: %dx%d",A0.dim1,A0.dim2);
      wblog(FL," *  sp: %dx%d",mspec.dim1,mspec.dim2);
   }
};

template <class TS>
void Spectral<TS>::getRawData(
   wbvector<double> &om, wbMatrix<TS> &A, unsigned linfit,
   const char *isbuf
) const {
   unsigned i,s,Nom=Ap.dim2;

   getOmega(om, linfit); 

   A.init(Ap.dim1, 2*Nom);

   if (A.dim2!=om.len || !Ap.sameSize(An)) wblog(FL,
   "ERR severe dimension mismatch (%d,%d).", A.dim2, om.len);

   if (isbuf && isbuf[0]) {
      if (strcmp(isbuf,"buf"))
      wblog(FL, "ERR invalid flag `%s'", isbuf);

      if (!Ap.sameSize(Ap_buf)) {
         if (Ap_buf.isEmpty())
         wblog(FL, "ERR buffer not set !?"); else
         wblog(FL, "ERR severe size mismatch !?");
      }

      for (s=0; s<Ap.dim1; s++) { A(s,Nom)=NAN;
         for (i=0; i<Nom; ++i) {
            A(s,Nom-i-1) = An_buf(s,i);
            A(s,Nom+i  ) = Ap_buf(s,i);
         }
      }
   }
   else {
      for (s=0; s<Ap.dim1; s++) { A(s,Nom)=NAN;
         for (i=0; i<Nom; ++i) {
            A(s,Nom-i-1) = An(s,i);
            A(s,Nom+i  ) = Ap(s,i);
         }
      }
   }

   if (linfit) {
      for (s=0; s<A.dim1; s++) {
          Wb::set2avg(A.rec(s)+Nom,        linfit);
          Wb::set2avg(A.rec(s)+Nom-linfit, linfit);
      }
   }
};

template <class TS>
void Spectral<TS>::pairRawData() {
   unsigned i,j,k, m=Ap.dim1, n=Ap.dim2;

   if (!Ap.sameSize(An)) wblog(FL, "ERR severe size mismatch!");
   if (Ap.dim1%2) wblog(FL,
   "ERR Ap has odd number of rows !? (%d)", Ap.dim1);

   for (i=0; i<m; i+=2)
   for (k=i/2, j=0; j<n; j++) {
      Ap(k,j)=Ap(i,j)+Ap(i+1,j);
      An(k,j)=An(i,j)+An(i+1,j);
   }

   Ap.Resize(m/2,n);
   An.Resize(m/2,n);

   if (mspec.dim1==0) return;

   m=mspec.dim1; n=mspec.dim2;

   if (m%2) wblog(FL,
   "ERR mspec has odd number of rows !? (%d)", m);

   for (i=0; i<m; i+=2)
   for (k=i/2, j=0; j<n; j++) mspec(k,j)=mspec(i,j)+mspec(i+1,j);

   mspec.Resize(m/2,n);
}

template <class TS>
void Spectral<TS>::getPARTIAL(
   wbarray<TS> &A1, wbarray<TS> &A2
 ) const {

   if (AP.isEmpty() && AN.isEmpty()) {
      A1.init(); A2.init(); return;
   }

   if (AP.SIZE.len!=3 || AP.SIZE!=AN.SIZE || AP.SIZE[1]%2) wblog(FL,
     "ERR %s() invalid initialization of AP/AN [%s; %s; %dx%d]",
      FCT,AP.sizeStr().data, AN.sizeStr().data, Ap.dim1, Ap.dim2
   ); 

   int i,j, m=AP.SIZE[1]/2, n=AP.SIZE[0], K=AP.SIZE[2], n1=n-1, n2=2*n*m;
   const TS *ap=AP.data, *an=AN.data;
   TS *a1, *a2;

   A1.init(2*n,m,K); a1=A1.data; 
   A2.init(2*n,m,K); a2=A2.data;

   for (K*=m, j=0; j<K; ++j) {
      if (j>=m) {
         for (i=0; i<n; ++i) { a1[i]=an[n1-i] - an[n1-i-n2]; }; a1+=n; an+=n;
         for (i=0; i<n; ++i) { a1[i]=ap[   i] - ap[   i-n2]; }; a1+=n; ap+=n;
         for (i=0; i<n; ++i) { a2[i]=an[n1-i] - an[n1-i-n2]; }; a2+=n; an+=n;
         for (i=0; i<n; ++i) { a2[i]=ap[   i] - ap[   i-n2]; }; a2+=n; ap+=n;
      }
      else {
         for (i=0; i<n; ++i) { a1[i]=an[n1-i]; }; a1+=n; an+=n;
         for (i=0; i<n; ++i) { a1[i]=ap[   i]; }; a1+=n; ap+=n;
         for (i=0; i<n; ++i) { a2[i]=an[n1-i]; }; a2+=n; an+=n;
         for (i=0; i<n; ++i) { a2[i]=ap[   i]; }; a2+=n; ap+=n;
      }
   }
};

template <class TS>
void Spectral<TS>::detailedBalance(
   const double T,
   wbvector<double> &om, wbMatrix<TS> &aa
) const {

   unsigned i,j,k1,k2,m,n;
   wbvector<double> f1,f2;
   wbMatrix<TS> ar;

   if (T<0) wblog(FL,"ERR T=%g !?", T);

   getRawData(om,ar); 

   n=Ap.dim2; f1.init(2*n);

   if (T==0) {
      for (i=0; i<n; i++) f1[i]=1;
   }
   else {
      double beta=1./T;
      for (i=0; i<om.len; i++) f1[i]=1./(1.+exp(-beta*om[i]));
   }

   f1.flip(f2);

   m=Ap.dim1; n=ar.dim2;
   aa.init(2*m,n);

   for (i=0; i<m; i++) {
      k1=2*i; k2=k1+1;

      for (j=0; j<n; j++) {
         aa(k1,j)=ar(i,j)*f2[j]; 
         aa(k2,j)=ar(i,j)*f1[j]; 
      }
   }
}

template <class TS>
void Spectral<TS>::dispSumRule() {

   unsigned i; double w=0;
   wbvector<TS> sp=Ap.recSum(),sn=An.recSum();
   Wb::termcolor Tc(Wb::TCOLS::BLUE);

   wblog(FL,"\r%60s%N"
       "   sum %sp_raw  : %s%N"      
       "    +  %sn_raw  : %s ","",   
       istr.data, sp.toStrf("%8.5f").data,
       istr.data, sn.toStrf("%8.5f").data
   );
   sp+=sn;

   if (!A0.isEmpty()) {
      if (A0.dim1!=sp.len) wblog(FL,
      "ERR size mismatch (%d/%dx%d)",sp.len,A0.dim1,A0.dim2);
      sp+=A0.recSum();
   }

   if (sp.len) { w=REAL(sp[0]); }
   for (i=0; i<sp.len; ++i) { if (ABS(sp[i]-w)>1E-2) break; }

   if (i==sp.len && fabs(w)>1E-4 && fabs(w-round(w))<1E-4) { sp-=w;
      wblog(FL,"\r%60s\r   %s\\","",Tc.e1);
      printf("total%-+7g: %s",-w,sp.toStrf("%8.2g").data);
      printf("%s\n",Tc.em);
   }
   else {
      printf("       %stotal   : ",Tc.e1);
      for (i=0; i<sp.len; ++i) { w=round(REAL(sp[i]));
         if (ABS(sp[i]-w)<1E-5) {
            if (w!=0)
                 printf("%2g%+5.0E ",w,REAL(sp[i])-w);
            else printf("%8.1E ",REAL(sp[i]));
         }
         else { printf("%8.5f ",REAL(sp[i])); }
      }; printf("%s\n",Tc.em);
   }
};

template <class TS>
void Spectral<TS>::dispSumRule(
   const wbvector<double> &om, const wbMatrix<TS> &A
){
   wbvector<TS> s(A.dim1);
   unsigned i; double m=0;

   if (om.len!=A.dim2) wblog(FL,
      "ERR Dimension mismatch (%d,%d)", om.len, A.dim2);

   for (unsigned i=0; i<A.dim1; ++i)
   s[i] = Wb::IntTrapez(om.data, A.rec(i), om.len);

   wblog(FL,"\r%60s\r   "
     "int AA(om)  : %s","",s.toStrf("%8.5f").data);

   if (!A0.isEmpty()) {
      if (A0.dim1!=s.len) wblog(FL,
      "ERR size mismatch (%d/%dx%d)",s.len,A0.dim1,A0.dim2);
      s+=A0.recSum();
   }

   if (s.len) {
      for (m=round(REAL(s[0])), i=0; i<s.len; ++i) {
        if (ABS(s[i]-m)>1E-2) break;
      }
      if (m && i==s.len) { s-=m;
         char ms[8]; Wb::termcolor Tc(Wb::TCOLS::BLUE);
         if (m) {
            i=snprintf(ms,8,"%+3g",-m);
            if (int(i)<0) wblog(FL,"WRN %s() i=%d",FCT,i);
         }
         else { *ms=0; }
         wblog(FL,"\r%60s\r    %s*  total%-3s: %s%s",
           "",Tc.e1, ms, s.toStrf("%8.2g").data,Tc.em);
      }
   }
};

template <class TS>
void Spectral<TS>::getSmoothSpec(
   wbvector<double> &om, wbMatrix<TS> &A, char vflag
){
   unsigned linfit=0;

   getOmega(om,linfit);
   getSmoothSpec(A,om,vflag);
};

template <class TS>
void Spectral<TS>::getSmoothSpec(
    wbMatrix<TS> &A, const wbvector<double> &omega,
    char vflag
){
    unsigned i,s,Nom=Ap.dim2;
    wbvector<double> OM(Nom), I;
    double x, dbl, om, aom;
    TS w1=0., w2=0.;

    if (fac<=0 || sigma<=0 || sigma2<=0 || eps<=0 || Nom<4)
    wblog(FL,"ERR %s()\n"
       "fac=%.3g, sigma=%.3g, sigma2=%.3g, eps=%.3g, Nom=%d !?",
        FCT,fac,sigma,sigma2, eps,Nom
    );

    for (i=0; i<OM.len; i++) OM[i]=idx2om(i);

    if (vflag) {
       double w1=omega.min(), w2=omega.max(), w0=omega.aMin();
       if (w1==-w2)
            sprintf(str,        "%.3g .. %.3g",   w0,w2);
       else sprintf(str,"%.3g .. %.3g .. %.3g",w1,w0,w2);

       wblog(FL," *  "
         "%s() %g bins/dec (len=%d)\n"
         "sigma=%g [eps=%.3g @ %.3g]\n" 
         "|omega| range: %s",
          myname,fac/log10(M_E),Nom,sigma,eps,sigma2,str);
       if (fabs((w0-idx2om(0    ))/w0)>0.01 ||
           fabs((w2-idx2om(Nom-1))/w2)>0.01) wblog(FL,
         "|emin...emax|= %.3g .. %.3g (bins)",idx2om(0),idx2om(Nom-1));
       if (vflag && (Ap.dim1<8 || vflag>1)) { dispSumRule(); }
    }  

    A.init(Ap.dim1, omega.len);

    for (s=0; s<Ap.dim1;   ++s)
    for (i=0; i<omega.len; ++i) { om=omega[i]; aom=fabs(om);

        if (aom>eps) x=1.;
        else {
            if (aom==0.) x=0.;
            else {
               dbl = (log(aom)-log(eps))/sigma;
               x = exp(-dbl*dbl);
               if (x<1E-10) x=0.;
            }
        }

        if (x!=0.)
        w1 = foldLogGaussian(aom, om>0 ? Ap.rec(s) : An.rec(s), Nom);

        if (x!=1.) {
           double b=MAX(sigma2*eps,om*exp(-sigma*sigma/4));

           w2 = foldGaussian( om, b, OM.data, Ap.rec(s), Nom)
              + foldGaussian(-om, b, OM.data, An.rec(s), Nom);
        }

        A(s,i) = x*w1 + (1.-x)*w2;
    }

    if (vflag && (Ap.dim1<8 || vflag>1)) dispSumRule(omega,A);
};

template <class TS>
TS Spectral<TS>::foldLogGaussian(
   const double om, const TS *aa, const unsigned n
){
   double k, dbl, efac;
   TS xa=0.;

   if (sigma<=0. || fac<=0.) wblog(FL,
      "ERR sigma=%g, fac=%g !?", sigma, fac);

   k=OM2IDX(om)-0.5; 
   efac=1./(fac*sigma);

   double dx=alpha*sigma; 

   for (unsigned i=0; i<n; ++i) {
       dbl=efac*(k-double(i))-dx;
       xa += aa[i]*exp(-dbl*dbl);
   }
   xa /= (sigma*sqrt(M_PI)*om);

   if (alpha!=0.25) {
      xa *= exp(+sigma*sigma*(alpha-0.25));
   }

   return xa;
};

template <class TS>
TS Spectral<TS>::foldGaussian(
   const double om, const double b,
   const double *oo, const TS *aa, const unsigned n
){
   double dbl, efac;
   TS xa=0.;

   if (b<=0.) wblog(FL, "ERR b=%g !?", b);

   efac = 1./b;

   for (unsigned i=0; i<n; i++) {
       dbl=efac*(om-oo[i]);
       xa += aa[i]*exp(-dbl*dbl);
   }

   xa *= 1. / (b*sqrt(M_PI));

   return xa;
};

template <class TS>
void Spectral<TS>::KKreal(
   wbMatrix<wbcomplex> &GZ, wbvector<double> &OM, wbvector<double> &DOM
) const {

   unsigned i,j,s, d=Ap.dim1, Nom;
   wbMatrix<double> II, RR;
   double doi, doj;

   getRawData(OM, II); Nom=OM.len;
   Wb::getDiff(OM,DOM);

#if 0
   wblog(FL,"TST Get smooth data ...");

   getSmoothSpec(OM, II, 1E-6, 0.3, 0); Nom=OM.len;
   Wb::getDiff(OM,DOM);

   for (s=0; s<d; s++)
   for (i=0; i<Nom; i++) II(s,i) *= DOM[i];
#endif

   RR.init(d,Nom);

   for (i=0;   i<Nom; ++i)
   for (j=i+1; j<Nom; ++j) { 

       doi=DOM[i]/(OM[j]-OM[i]);
       doj=DOM[j]/(OM[i]-OM[j]);

       for (s=0; s<d; s++) {
           RR(s,i) += II(s,j)*doi;
           RR(s,j) += II(s,i)*doj;
       }
   }

#if 0
   wbMatrix<double> Ip2, In2;
   wblog(FL,"TST backtransform real->imag");

   RR*=(1/M_PI); RR.put("RR","base",'r');

   KKimag(Rp,Rn,OM,Ip2,In2);

   Ip2.put("Ip","base",'r');
   In2.put("In","base",'r');

   wblog(FL,"ERR"); 
#endif

   II *= M_PI; 

   GZ.set(RR,II);
}

void KKimag(
   const wbMatrix<double> &RR,
   const wbvector<double> &OM,
   wbMatrix<double> &II
){
   unsigned i,j,s, d=RR.dim1, Nom=RR.dim2;
   double doi, doj;
   wbvector<double> DOM;

   Wb::getDiff(OM, DOM); II.init(d,Nom);

   if (OM.len!=Nom) wblog(FL,
      "ERR severe dimension mismatch (%d,%d) !?", OM.len, Nom);

   for (i=0;   i<Nom; ++i)
   for (j=i+1; j<Nom; ++j) { 
       doi=DOM[i]/(OM[i]-OM[j]); doj=DOM[j]/(OM[j]-OM[i]);
       for (s=0; s<d; s++) {
           II(s,i) += RR(s,j)*doi;
           II(s,j) += RR(s,i)*doj;
       }
   }

   II *= (1./M_PI);
}

template <class TS>
void TDSpectral<TS>::init(const char *F, int L, const char *istr0,
   const wbvector<double> &t0,
   unsigned Nnrg, unsigned Nop,
   double Emin, double Emax, unsigned nl, char vflag
){
   double l2;

   if (Emin<=0. || Emax<=0. || Emin>=0.1*Emax || nl==0) wblog(F,L,
      "ERR invalid parameter set [%g %g %d]%s", Emin, Emax, nl,
      (Emin>0 && Emax>0 && Emin<Emax)? " (require Emin<Emax/10)":"");

   l0 = log(Emin);
   l2 = log(Emax);

   fac= log10(M_E)*double(nl); 

   N = (unsigned)((l2-l0)*fac)+1;

   AP.init(N,Nop,Nnrg);
   AN.init(N,Nop,Nnrg); if (check_avg_om) {
   AW.init(2,Nop,Nnrg); }

   if (vflag) {
      wblog(F_L,"%3s TDSpectral: %.2g .. %.2g (%d/dec; 2x%d)",
      istr0 ? istr0 : "<i>", Emin, Emax, nl, N);
   }

   istr=istr0; tt=t0;
};

template <class TS>
wbarray<TS> TDSpectral<TS>::getAvgOm() {

   if (AW.rank()!=3 || AW.dim(1)!=2) wblog(FL,
      "ERR invalid AW data structure [%s] !?",AW.sizeStr().data);

   unsigned i=0, n=AW.numel();
   wbarray<TS> aw(AW); if (aw.isEmpty()) return aw;
   TS *a=AW.data, *b=aw.data;

   for (; i<n; i+=2) { b[i] = ((a[i+1]!=0.) ? a[i]/a[i+1] : a[i]); }

   return aw;
};

template <class TS>
void TDSpectral<TS>::Add(
   const wbarray<double> &Om,
   const wbarray<TS> &CC, unsigned iter, unsigned s
){
   static int iter_last=-1;
   unsigned i,k,n=Om.numel();
   TS *ap, *an;

   if (!Om.sameSize(CC)) {
      Om.info("Om"); CC.info("CC");
      wblog(FL,"ERR size mismatch between Om and CC!");
   }

   if (s>=AP.dim(2) || iter>=AP.dim(3)) wblog(FL,
      "ERR index out of bounds (op=%d/%d; k=%d/%d)",
       s+1,AP.dim(2),iter+1,AP.dim(3));

   ap=&AP(0,s,iter);
   an=&AN(0,s,iter);

   for (i=0; i<n; ++i) {
       k=om2idx(Om.data[i]);
       if (Om.data[i]<0.)
            an[k]+=CC.data[i];
       else ap[k]+=CC.data[i];
   }

   if (raw) {
      if (iter_last!=(int)iter) {
wblog(FL,"%NTST iter=%d",iter);
         if (RAW.len) addPartialFourier();
         RAW.initDef(AP.dim(2));
wblog(FL,"TST done");
      }

      RAW[s].Append(Om.data, CC.data, n);

      iter_last=iter;
   }

   if (check_avg_om) {
      if (AW.rank()!=3 || AW.dim(1)!=2 ||
          AW.dim(2)!=AP.dim(2) || AW.dim(3)!=AP.dim(3)) wblog(FL,
         "ERR size mismatch for AW [%s] ? [%s]",
          AW.sizeStr().data, AP.sizeStr().data
      );
      AW(0,s,iter)+=CC.dotProd(Om); 
      AW(1,s,iter)+=CC.sum();
   }
};

template <class TS> inline
void TDSpectral<TS>::addPartialFourier() {

   if (!RAW.len) return;

   unsigned s=0, m=RAW.len;
   wbMatrix<wbcomplex> at;
   wbMatrix<TS> aa;

   if (m!=AP.dim(2)) wblog(FL,
      "ERR severe size mismatch (%d/%d)",m,AP.dim(2));

   if (At.isEmpty()) At.init(AP.dim(2),tt.len);

   for (s=0; s<m; ++s) {
      const RSpecData<TS> &R=RAW[s];
      if (R.om.len!=R.cc.len) wblog(FL,
         "ERR RAW[%d].len=%d/%d !?",s,R.om.len,R.cc.len);

      aa.init2ref(1,R.cc.len,R.cc.data);

      Wb::Fourier(R.om,aa,tt,at); 
      At.recAddP(s,at.data);
   }

   RAW.init();
};

template <class TS>
void TDSpectral<TS>::getOmega(
   wbvector<double> &om, unsigned &linfit
) const {

   unsigned i; 
   double dbl;

   if (linfit) linfit=unsigned(fac);
   if (linfit>N) wblog(FL,"ERR linfit>N (%d,%d) !?", N, linfit);

   om.init(2*N);

   dbl = idx2om(linfit) / (linfit+0.5);

   for (i=0; i<linfit; i++) om[N-1-i] = -(om[N+i] = dbl*(i+0.5));
   for (   ; i<N;      i++) om[N-1-i] = -(om[N+i] = idx2om(i)  );
}

template <class TS>
void TDSpectral<TS>::applyIROPfac(
   const char *F, int L, const wbvector<double> &fac, char rixs) {

   unsigned i,k=0, n=fac.len;
   double x;

   for (i=0; i<n; ++i) { if (fac[i]!=1) {
      if (fac[i]<1 || fac[i]!=round(fac[i])) wblog(FL,
         "ERR %s() got fac[%d]=%g !?",FCT,i+1,fac[i]);
      else {
         wblog(FL," *  applying IROP factor 1/%g%s to op(%d)",
            fac[i], rixs? "^2":"",i+1);
         ++k;
      }
   }}
   if (!k) return;

   if (!Ap.isEmpty() || !An.isEmpty()) {
      if (n!=Ap.dim1 || n!=An.dim1) wblog(FL,
         "ERR %s() size mismatch (%d,%d/2*%d)",FCT,Ap.dim1,An.dim1,n);

      for (i=0; i<n; ++i) {
         if (fac[i]!=1) { x=1/fac[i]; if (rixs) { x*=x; }
            Ap.rowTimes(i,x);
            An.rowTimes(i,x);
         }
      }
   }

   if (!AP.isEmpty() || !AN.isEmpty()) {
      if (AP.SIZE.len!=3 || !AP.sameSize(AN) || AP.SIZE[1]!=n)
         wblog(FL,"ERR %s() size mismatch (%s; %s /%d) !?",
         FCT, SSTR(AP), SSTR(AN), n);

      TS *ap=AP.data, *an=AN.data;
      unsigned n0=AP.SIZE[0], m=n*AP.SIZE[2];

      for (i=0; i<m; ++i, ap+=n0, an+=n0) { k=i%n;
         if (fac[k]>1) { x=1/fac[k]; if (rixs) { x*=x; }
            Wb::timesRange(ap, x, n0);
            Wb::timesRange(an, x, n0);
         }
      }
   }
};

template <class TS>
void TDSpectral<TS>::getSpecData(
   wbvector<double> &om, wbarray<TS> &A,
   unsigned linfit
) const {
   unsigned i,r,s, d1=AP.dim(2), d2=AP.dim(3); 

   getOmega(om,linfit); 
   A.init(2*N,d1,d2);

   if (om.len!=2*N || !AP.sameSize(AN)) wblog(FL,
      "ERR severe dimension mismatch (%d,%d).",2*N,om.len);

   for (r=0; r<d1; r++)
   for (s=0; s<d2; s++) {
      for (i=0; i<N; i++) {
         A(N-i-1,r,s) = AN(i,r,s);
         A(N+i  ,r,s) = AP(i,r,s);
      }
   }

   if (linfit) {
      for (r=0; r<d1; r++)
      for (s=0; s<d2; s++) {
          Wb::set2avg(A.ref(r,s)+N,        linfit);
          Wb::set2avg(A.ref(r,s)+N-linfit, linfit);
      }
   }
};

template <class TS>
void TDSpectral<TS>::getSpecData(
   wbvector<double> &om, wbMatrix<TS> &A,
   unsigned linfit
){

   unsigned i,r; 

   getOmega(om,linfit); 

   if (Ap.isEmpty()) mergeNRGData();
   A.init(Ap.dim1,2*N);

   if (A.dim2!=om.len || !Ap.sameSize(An)) wblog(FL,
      "ERR Severe dimension mismatch (%d,%d).",A.dim2,om.len);

   for (r=0; r<A.dim1; r++)
   for (i=0; i<N; i++) {
      A(r,N-i-1) = An(r,i);
      A(r,N+i  ) = Ap(r,i);
   }

   if (linfit) {
      for (r=0; r<A.dim1; r++) {
          Wb::set2avg(A.rec(r)+N,        linfit);
          Wb::set2avg(A.rec(r)+N-linfit, linfit);
      }
   }
}

template <class TS>
void TDSpectral<TS>::Fourier(
   wbvector<double> &om,
   wbMatrix<wbcomplex> &ST,
   double alpha, double Lambda, char vflag
){
   unsigned k,it,ik,is,io,l=0,no, nt=tt.len, ns=AP.dim(2), nk=AP.dim(3);

   wbarray<wbcomplex> TKS(nt,nk,ns), KST;
   wbMatrix<double>   aTK(nt,nk); 
   wbarray<TS> AR;
   double ek,oi,*at; TS kso,*ap;
   wbcomplex z,*st;

   if (vflag) wblog(FL,
      "<i> %s (damping alpha=%g, Lambda=%g)",FCT,alpha,Lambda);

   getSpecData(om,AR);   
   ap=AR.data; no=om.len;

   if (alpha==0) aTK.set(1); else
   for (ik=0; ik<nk; ++ik) {
      ek=pow(Lambda, -double(ik)/2.) * 0.5*(Lambda+1);
      for (it=0; it<nt; it++)
      aTK(it,ik)=exp(-fabs(alpha*ek*tt[it]));
   }

   for (ik=0; ik<nk; ik++)
   for (is=0; is<ns; is++)
   for (io=0; io<no; io++, l++) if (ap[l]!=0.) { kso=ap[l]; oi=om[io];
       for (it=0; it<nt; it++)
       TKS(it,ik,is)+=kso*expi(oi*tt[it]);
   }

   TKS.permute(KST,"2 3 1");
   ST.init(ns,nt);

   for (is=0; is<ns; ++is)
   for (it=0; it<nt; ++it) {
      st=KST.ref(is,it); z=0;
      if (alpha==0)
           for (k=0; k<nk; ++k) z+=st[k];
      else for (at=aTK.rec(it), k=0; k<nk; ++k) z+=at[k]*st[k];
      ST(is,it)=z;
   }

};

template <class TS>
void TDSpectral<TS>::Fourier(
   wbvector<double> &om, wbMatrix<TS> &A,
   const wbvector<double> &tt,
   wbMatrix<wbcomplex> &at,
   double alpha, char vflag
){
   if (Ap.isEmpty()) { mergeNRGData(); dispSumRule(); }
   getSpecData(om,A);

   Wb::Fourier(om,A,tt,at,0,alpha,vflag); 
};

template <class TS>
void Wb::Fourier(
   const wbvector<double> &om, const wbMatrix<TS> &A,
   const wbvector<double> &tt,
   wbMatrix<wbcomplex> &az,
   char domflag,  
   double alpha,
   char vflag
){
   unsigned i,j,l,m,n,b1,b2, nom, nt=tt.len, D=512, stop=0;
   wbMatrix<wbcomplex> ez;
   wbcomplex z;
   double f,ot;

   wbvector< wbMatrix<wbcomplex> > AZ;
   wbvector<double> dom(om.len);

   m=A.dim1; nom=om.len;

   if (nom!=A.dim2) wblog(FL,
      "ERR %s - size mismatch (%d/%d)",nom,A.dim2);
   if (!tt.len || !m || !nom) {
      wblog(FL,"WRN Empty data [%d %d %d] - return",tt.len,m,nom);
      return;
   }
   if (nom<2) wblog(FL,"ERR %s - omega too short (%d/2)",nom);

   n=om.len-1;
   for (i=1; i<n; i++) dom[i]=0.5*(om[i+1]-om[i-1]);
   dom[0]=om[1]-om[0]; dom[n]=om[n]-om[n-1];

   for (i=0; i<dom.len; i++) if (dom[i]<=0) { wblog(FL,
      "WRN omega not strictly ascending !?"); break; } 

   for (i=0; i<tt.len; i++) if (tt[i]<0) { wblog(FL, 
      "WRN negative tt values !?"); break; }

   if (domflag) {
      wbMatrix<double> B(A);

      for (i=0; i<B.dim1; ++i)
      for (j=0; j<B.dim2; ++j) B(i,j)*=dom[j];

      Wb::Fourier(om,B,tt,az,0,alpha,vflag);
      return;
   }

   if (vflag) wblog(FL," *  %s() spectral data "
      "(alpha=%.4g, %s)",FCT,alpha, domflag? "func":"disc");

   AZ.initDef(nt<D ? 1 : nt/D); 

   for (l=b1=0; !stop; l++,b1+=D) {
      b2=b1+D-1; if (b2+D>=nt) { b2=nt-1; stop=1; }
      n=b2-b1+1;

      if (n*nom*sizeof(double)/(1<<20)>900) { 
         wblog(FL,"WRN trying to allocate %.6g MB !?", 
         n*nom*sizeof(double)/double(1<<20));
      }

      ez.init(nom,n);
      for (i=0; i<nom; i++)
      for (j=0; j<n; j++) { ot=om[i]*tt[b1+j];
         ez(i,j).init_expi(ot); 

         if (alpha>0) {    
            f=::exp(-::fabs(alpha*dom[i]*tt[b1+j]));
            ez(i,j)*=f;
         }
      }
      Wb::MatProd(A,ez,AZ[l]);
   }

   az.CAT(2,AZ);
};

template <class TS>
void TDSpectral<TS>::mergeNRGData() {

   unsigned i,k, d1=AP.dim(1), d2=AP.dim(2), d3=AP.dim(3);
   unsigned s=d1*d2, S=s*d3; 
   TS *Dp,*Dn,*dp,*dn;

   if (!AP.sameSize(AN)) wblog(FL, "ERR severe size mismatch!");

   Ap.init(d2,d1); dp=Ap.data; Dp=AP.data;
   An.init(d2,d1); dn=An.data; Dn=AN.data;

   for (i=0; i<S; i++) { k=i%s;
      dp[k]+=Dp[i];
      dn[k]+=Dn[i];
   }
};

template <class TS>
void TDSpectral<TS>::dispSumRule() {

   unsigned i=0; double m=0;
   wbvector<TS> sp,sn;
   wbarray<TS> x;

   AP.sum("1,3",x); sp.init(x.numel(),x.data);
   AN.sum("1,3",x); sn.init(x.numel(),x.data);

   wblog(FL," *  om>0: %s\nom<0: %s",
      sp.toStrf("%8.5f").data, sn.toStrf("%8.5f").data);

   sp+=sn;

   if (sp.len) {
      m=ROUND(REAL(sp[0]));
      for (i=0; i<sp.len; i++) if (ABS(sp[i]-m)>1E-2) break;
   }

   Wb::termcolor Tc(Wb::TCOLS::BLUE);
   if (i==sp.len) { sp-=TS(m); wblog(FL,
   " *  %stotal %s%s (%g)",Tc.e1,sp.toStrf("%8.2g").data,Tc.em,m); }
   else wblog(FL,
   " *  %stotal %s%s",     Tc.e1,sp.toStrf("%8.5f").data,Tc.em);
};

template <class TS>
void Wb::dispSumRule(
   const wbvector<double> &om,
   const wbMatrix<TS> &A
){
   unsigned i; double m=0;
   wbvector<double> s(A.dim1);
   Wb::termcolor Tc(Wb::TCOLS::BLUE);

   if (om.len!=A.dim2) wblog(FL,
      "ERR Dimension mismatch (%d,%d)", om.len, A.dim2);

   for (unsigned i=0; i<A.dim1; i++)
   s[i] = Wb::IntTrapez(om.data, A.rec(i), om.len);

   wblog(FL,"SUM %s", s.toStrf("%8.5f").data);

   if (s.len) {
      for (m=ROUND(s[0]), i=1; i<s.len; ++i) {
         if (ABS(s[i]-m)>1E-2) { break; }
      }
      if (i==s.len) { s-=m;
         wblog(FL,"%s *  total-%g     = %s%s",
         Tc.e1,m,s.toStrf("%8.2g").data,Tc.em);
      }
   }
}

template <class TS>
void TDSpectral<TS>::getSmoothSpec_t(
    double sigma, double eps, char bswitch, char vflag
) {
    unsigned i,s; 
    wbMatrix<double> Ap0(Ap), An0(An), wlg;

    if (fac<=0 || sigma<=0 || N<4) wblog(FL,
       "ERR invalid paramters (fac=%g, sigma=%g, N=%d)",fac,sigma,N);
    if (fac*sigma<1) wblog(FL,
       "WRN sigma < discrete log spacing (%.3g / %.3g)",sigma,1/fac);
    if (vflag) wblog(FL,
       "<i> %s() %dx%d bins (sigma=%g)",myname,AP.dim(2),N,sigma);

    if (Ap.isEmpty()) { mergeNRGData(); };
    if (vflag) { dispSumRule(); }

    Ap.set(0); An.set(0);
    initLG(wlg,N,sigma,bswitch,vflag);

    for (s=0; s<Ap.dim1; ++s)
    for (i=0; i<Ap.dim2; ++i) {
       foldLogGaussian(wlg, Ap.rec(s), Ap0(s,i), i);
       foldLogGaussian(wlg, An.rec(s), An0(s,i), i);
    }

    if (eps) {
       if (vflag>1) {
       wblog(FL," *  %s() eps=%.3g/%.3g",myname,eps,exp(l0)); }

       unsigned j,l, i0=0, i1=0, i2=0;
       unsigned nE=(1<<12); 

       double eps0=0.05*eps, eps1=1.5*eps, x1=5., xfac=x1/double(nE);
       double x=1./sqrt(M_PI), dbl, lfac=1/xfac;

       wbvector<double> om(Ap.dim2), dom(om.len);
       wbvector<double> gauss_lin(2*nE+1), erf_om(om.len);

       if (om.len<10) wblog(FL,"ERR %s() om.len = %d !?",FCT,om.len);
       for (i=0; i<om.len; ++i) {
          om[i]=idx2om(i);

          if (!i0 && i && om[i]>eps0) { i0=i; }

          if (!i1 && i && (om[i]-om[i-1])>eps1) { i1=i; }

          if (i>1) {
          dom[i-1] = 0.5*( om[i] - om[i-2] ); }
       }
       dom[i-1] = 2*dom[i-2]; 
       dom[0] = 0.5*( om[0] + om[1] ); 

       if (i1) { if (vflag) {
          sprintf(str,"(%.1g <) |om| < %.1g (%.0f%%)",om[i0],om[i1],
             100*((i1-i0)/double(om.len)));
          wblog(FL,"--> Gaussian broadening for %s",str);
       }}
       else {
          i1=om.len-1; sprintf(str,"(%.3g; eps=%.3g)",om[i1],eps);
          wblog(FL,"WRN Gaussian broadening for all omega %s",str);
       }

       i2=i1+10; if (i2>=om.len) { i2=om.len-1; }

       for (i=0; i<erf_om.len; ++i) {
          erf_om[i] = (0.5*( std::erf( om[i]/eps ) + 1 ));
       }

       for (i=0; i<gauss_lin.len; ++i) {
          dbl = i*xfac - x1;        
          gauss_lin[i] = x*exp(-dbl*dbl);
       }
       dom*=(1/eps);                             

       Ap0=Ap; An0=An;
       for (s=0; s<Ap.dim1; ++s)
       for (i=0; i<i1; ++i) { Ap(s,i) = An(s,i) = 0; }

       for (s=0; s<Ap.dim1; ++s) {
          for (x=0., i=0; i<i0; ++i) { x += (Ap0(s,i) + An0(s,i)); }
          for (i=0; i<i2; ++i) {
              Ap(s,i) += ( dbl = x*(erf_om[i] - (i? erf_om[i-1] : 0.5)) );
              An(s,i) += dbl;
          }
       }

       for (s=0; s<Ap.dim1; ++s)
       for (i=i0; i<i1; ++i) {
          for (j=0; j<i2; ++j) {
             dbl = (om[j]-om[i])/eps;
             if (fabs(dbl)<x1 ) { l = (dbl + x1) * lfac + 0.5;
             if (l<gauss_lin.len) {
                Ap(s,j) += Ap0(s,i) * (dbl = gauss_lin[l] * dom[j]);
                An(s,j) += An0(s,i) *  dbl;
             }}

             dbl = (om[j]+om[i])/eps;
             if (fabs(dbl)<x1 ) { l = (dbl + x1) * lfac + 0.5;
             if (l<gauss_lin.len) {
                An(s,j) += Ap0(s,i) * (dbl = gauss_lin[l] * dom[j]);
                Ap(s,j) += An0(s,i) *  dbl;
             }}
          }
       }
    }

};

template <class TS>
void TDSpectral<TS>::initLG(
   wbMatrix<double> &wlg,
   unsigned n, double sigma, char bswitch, char vflag
){
   unsigned i,j;
   double dbl,efac,s=0; 

   if (sigma<=0. || fac<=0.) wblog(FL,"ERR sigma=%g, fac=%g",sigma,fac);

   efac=1./(fac*sigma); wlg.init(2,n);

   switch (bswitch) {
      case 0:
          strcpy(str,"symmetric in log. space");
          s=0; break;
      case 1:
          strcpy(str,"symmetric log-gauss");
          s=sigma/2.; 
          break;
      case 2: 
          strcpy(str,"norm preserving log-gauss");
          s=sigma/4.; 
          break;
      default:
      wblog(FL,"ERR invalid bswitch=%d",bswitch);
   }
   if (vflag) wblog(FL,"--> %s() using %s",FCT,str);

   for (i=0; i<n; i++) {
       dbl=efac*double(i)+s; wlg(0,i)=exp(-dbl*dbl); 
       dbl=efac*double(i)-s; wlg(1,i)=exp(-dbl*dbl); 
   }

   wlg*=(1./wlg.sum());

   for (i=0; i<2; ++i)
   for (j=0; j<n; ++j) if (fabs(wlg(i,j))<1E-20) wlg(i,j)=0;
};

template <class TS> inline
void TDSpectral<TS>::foldLogGaussian(
   const wbMatrix<double> &wlg,
   TS *aa, TS a0, unsigned i0
){
   unsigned i, n=wlg.dim2;
   const double *wp=wlg.rec(0), *wn=wlg.rec(1);
   double wi, nrm=0.;

   for (i=0; i<n; ++i) {
      nrm += (i0>i ? wp[i0-i] : wn[i-i0]); }
   nrm=a0/nrm;

   for (i=0; i<n; ++i) {
      wi=( i0>i ? wp[i0-i] : wn[i-i0] );
      if (wi) aa[i]+=nrm*wi;
   }
};

template <class TS>
mxArray* TDSpectral<TS>::toMx() const {
   mxArray *S=mxCreateStructMatrix(1,1,0,NULL);

   mxAddField2Scalar(FL,S,"name",istr.toMx());
   mxAddField2Scalar(FL,S,"AP", AP.toMx());
   mxAddField2Scalar(FL,S,"AN", AN.toMx());
   mxAddField2Scalar(FL,S,"fac",numtoMx(fac));
   mxAddField2Scalar(FL,S,"l0", numtoMx(l0 ));

   if (raw) {
   mxAddField2Scalar(FL,S,"RAW",RAW.toMx()); }

   return S;
};

template <class TS>
int TDSpectral<TS>::init(const mxArray *S,
   double &emin, double &emax, unsigned &nlog
){
   if (!mxIsStruct(S) || !Mx::IsScalar(S) ||
        mxGetFieldNumber(S,"name")<0 ||
        mxGetFieldNumber(S,"AP"  )<0 ||
        mxGetFieldNumber(S,"AN"  )<0 ||
        mxGetFieldNumber(S,"fac" )<0 ||
        mxGetFieldNumber(S,"l0"  )<0
   ) return 1;

   istr.init(mxGetField(S,0,"name"  ));
   AP  .init(mxGetField(S,0,"AP"));
   AN  .init(mxGetField(S,0,"AN"));
   mxGetNumber(mxGetField(S,0,"fac"),fac);
   mxGetNumber(mxGetField(S,0,"l0"), l0 );

   if (AP.SIZE!=AN.SIZE || AP.isEmpty() || AN.isEmpty()) return 2;
   if (fac<=0) return 3;

   N=AP.dim(1);

   emin=exp(l0);
   emax=idx2om(N);
   nlog=(unsigned)(fac/log10(M_E) + 0.5);

   wblog(FL,
     "<i> TDSpectral `%s' from RAW data: %dx%d bins [x%d]",
      istr.data, AP.dim(1), AP.dim(2), AP.dim(3));

   return 0;
};

#endif
