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

#ifndef __WB_NRGDATA_HCC__
#define __WB_NRGDATA_HCC__

/* ------------------------------------------------------------------ //
// change log

   Wb,Feb22,19: changed internally K,T labels to K,D labels,
   besides that T is also temperature, TK is Kondo tempearture
   TD is a template data type name, etc.

// ------------------------------------------------------------------ */

   const double DEPS2=1E-6*DBL_EPSILON;

   template <class TQ, class TD>
   class NRGData;

  #define NRG_FLEN 128
   const char NRG_EMPTY_STR[]="(NRG)";

   unsigned NRG_N, NRG_ITER;
   char NRG_FILE[NRG_FLEN];

   unsigned STRICT_ITER0=1;

   Wb::ClockSet fdmClocks;

   wbvector<double> gES;

inline double gEScale(int iter) {
   if (iter<0 || iter>=(int)gES.len) wblog(FL,
      "ERR %s() index out of bounds (%d/%d)",FCT,iter,gES.len);
   return gES[iter];
};

void wblog_ERR_updateOp(
   const char *F, int L, const char *istr, int i, int nop, int iter) {

   wblog(F_L,
     "ERR failed to update operator (%s-op %d/%d @ iter=%d)",
      istr && istr[0] ? istr : "(fdm)", i+1,nop,iter
   );
}

template <class TQ, class TD, class TZ>
void WeightE0_dim1(QSpace<TQ,TZ> &X, const QSpace<TQ,TD> &H, double w);

double nrgGetE0(
   const wbvector<double> &dE,
   const wbvector<double> &ES,
   wbvector<double> &E0 
){
   if (E0.len!=dE.len) { E0.init(dE.len); }
   if (!E0.len) { wblog(FL,"WRN %s() got empty E0 !?",FCT); return 0; }
   if (dE.len!=ES.len) wblog(FL,
      "ERR %s() length mismatch (%d/%d)",FCT,dE.len,ES.len);

   unsigned k=dE.len-1;
   double x=0, xl=dE[k]*ES[k]; E0[k]=0;

   for (; k; --k) {
      x=dE[k-1]*ES[k-1]; 
      E0[k-1]=E0[k]-xl;  
      xl=x;
   }

   return x-E0[0]; 
};

template <class TQ, class TD>
int getSpecDimOps(
   const wbvector<QSpace<TQ,TD> > &FF, wbvector<double> &fac
){
   unsigned i=0, r,m, n=FF.len, gotfac=0; widx_t nq;
   fac.init2val(n,1);

   for (; i<n; ++i) {
      r=FF[i].rank(FL); if (r<=2) continue;
      if (r>3) wblog(FL,"ERR %s() invalid rank-%g !?",FCT,r);

      m=FF[i].getDim(2,&nq); 
      if (int(nq)<1) wblog(FL,"ERR %s() %d/%d !?",FCT,m,nq);
      if (nq>1) { fac[i]=nq; ++gotfac; }
   }
   return gotfac;
};

class RIXS_para { 

  public:

    RIXS_para()
     : Eoffset(0.), sigma(0.), afac(-1.), wrxs(0),
       E0(0.), E2(0.), dloc(0), iter(0) {};

    void init() {
       Einc=0; Eoffset=sigma=E0=E2=0; wrxs=0; dloc=iter=0;
       afac=-1; 
    };

    void setE0(int it, double e0, double e2) {
       iter=it; E0=e0+Eoffset; E2=e2;  
    }

    void check(const char *F=0, int L=0, char init=0);

    Wb::String toStr() const;
    mxArray* toMx() const;

    wbcomplex Einc; 
    double Eoffset; 
    double sigma;   
    double afac;    
    char wrxs;      

    double E0, E2;  
    unsigned dloc;
    int iter;       

  private:
};

   RIXS_para gRX;

void RIXS_para::check(const char *F, int L, char init) {

   if (Einc.i<=0) {
      if (!Einc.i && init) {
         Einc.i=gES.last(); 
      }
      else wblog(F_L,
         "ERR %s() got imag(Einc)=%g < 0 (RIXS) !?",FCT,Einc.i
      );
   }
   else if (F && init) {
      if (Einc.i<1E-2*Einc.r || Einc.i>1E2*Einc.r) wblog(FL,
         "WRN %s() got Einc.i / Einc_r = %.3g / %.3g = %.3g",
         PROG, Einc.i, Einc.r, Einc.i/Einc.r
      );
   }

   if (sigma<0) wblog(F_L,
      "ERR RIXS::%s() invalid sigma=%g<0",FCT,sigma);
   if (int(dloc)<=1 || dloc>(1<<10)) wblog(F_L,
      "ERR RIXS::%s() invalid dloc=%d",FCT,dloc);

   if (init) {
      if (!sigma) {
         if (Einc.r) { sigma = Einc.i/Einc.r; }
         else wblog(FL,"ERR %s() got bot sigma and real(Einc) zero !?",FCT);
      }
      else if (sigma<1E-3 && F) {
         wblog(FL,"WRN %s() got sigma = %g",PROG,sigma);
      }
   }

};

Wb::String RIXS_para::toStr() const {

   Wb::String s(64); 
   s.pushf("Einc = %.3g",Einc.r);

   if (Einc.r) {
      double x=Einc.i/Einc.r;
      if (x && fabs((x-sigma)/x)<1E-4)
           s.pushf(" @ sigma=%.3g",sigma);
      else s.pushf(" @ %.3gi, sigma=%g",x,sigma);
   }
   else { s.pushf(" + %.3gi, sigma=%g",Einc.i,sigma); }

   if (wrxs   ) { s.pushf(", mode=%d",wrxs); }
   if (afac>=0) { s.pushf(", afac=%g",afac); }

   return s;
};

mxArray* RIXS_para::toMx() const {
    const char* fields[]={
    "Einc","Eoffset","sigma","wrxs","dloc"}; 

    mxArray* a=mxCreateStructMatrix(1,1,5,fields);
    mxSetFieldByNumber(a,0,0, numtoMx(Einc ));
    mxSetFieldByNumber(a,0,1, numtoMx(Eoffset));
    mxSetFieldByNumber(a,0,2, numtoMx(sigma));
    mxSetFieldByNumber(a,0,3, numtoMx(wrxs ));
    mxSetFieldByNumber(a,0,4, numtoMx(dloc ));

    return a;
};

template <class TQ, class TD>
class NRGIndex { 

  public:

    NRGIndex() : idx(-1), nrg(NULL) {};
    NRGIndex(const NRGIndex &I) : idx(I.idx), nrg(I.nrg), tag(I.tag) {};

    void init() { idx=-1; nrg=NULL; };
    void init(int i, NRGData<TQ,TD> *d) { idx=i; nrg=d; }
    void init(int i, const char *t)     { idx=i; tag=t; }

    mxArray* toMx() const;
    mxArray* mxCreateStruct(unsigned m, unsigned n) const;
    void add2MxStruct(mxArray *S, unsigned i, char tst=0) const;

    void put (const char *vname, const char *ws="base") const {
       mxArray *S=toMx();
       mxPutAndDestroy(FL,S,vname,ws);
    };

    int idx;              
    NRGData<TQ,TD> *nrg;  
    wbstring tag;         

  private:
};

template <class TQ, class TD>
mxArray* NRGIndex<TQ,TD>::mxCreateStruct(unsigned m, unsigned n) const {
   const char *fields[] = {"idx","nrg","tag"};
   return  mxCreateStructMatrix(m,n,3,fields);
};

template <class TQ, class TD>
void NRGIndex<TQ,TD>::add2MxStruct(
   mxArray *S, unsigned i, char tst __attribute__ ((unused))
 ) const {

   mxSetFieldByNumber(S,i,0, numtoMx(idx));
   mxSetFieldByNumber(S,i,1, numtoMx((unsigned long)(nrg)));
   mxSetFieldByNumber(S,i,2, mxCreateString(tag.isEmpty() ? tag.data:""));
};

namespace Wb {

   template <class TQ, class TD>
   int updateOp(const char *F, int L,
      const NRGData<TQ,TD> &NRG,
      wbvector< QSpace<TQ,TD> > &F12,
      const QSpace<TQ,TD> &A1,
      const QSpace<TQ,TD> &A2,
      const wbvector< QSpace<TQ,TD> > &FKK,
      unsigned iter
   );

   template<class TT>
   void initLocalOp(const char *F, int L,
      const mxArray* a, TT &C,
      int k=0, char force=0, char *name=NULL
   );

   template<class TQ, class TD>
   inline void initLocalOp_aux(const char *F, int L,
      const mxArray* a, wbvector< QSpace<TQ,TD> > &C
   ){ mxInitQSpaceVec(F,L,a,C); }

   template<class TQ, class TD>
   inline void initLocalOp_aux(const char *F, int L,
      const mxArray* a, wbMatrix< QSpace<TQ,TD> > &C
   ){ mxInitQSpaceMat(F,L,a,C); }

   template<class TQ, class TD>
   inline void initLocalOp_aux(const char *F, int L,
      const mxArray* a, wbvector< wbvector< QSpace<TQ,TD> > > &C
   ){ mxInitQSpaceVecVec(F,L,a,C); }
};

template <class TQ, class TD>
class NRGData { 
  public:

    NRGData(const char *s=NULL)
     : name(s), calc(0), store(0), kloc(0), locMX(0), MX(NULL), aux(NULL) {};

    NRGData(const NRGData &A) 
     : name(A.name), nrgIdx(A.nrgIdx), calc(0), store(0), kloc(A.kloc),
       locMX(0), MX(NULL), aux(NULL)
    { 
       wblog(FL,"ERR generating copy of NRGData !?");
    };

   ~NRGData() {
       if (MX && !NAME.isEmpty() && !locMX) wblog(FL,
          "WRN ** data return `%s' not implemented yet **\n(%lX; %d)",
           NAME.data, MX, locMX);
       init();
    };

    void initName(const char *n) { name=n; };

    void init(const char *F, int L, const mxArray* a, 
       wbvector< QSpace<TQ,TD> > &C0,
       int k=0, char force=0, const char *vname0=0
    );

    void init(const char *F, int L, const char *XX, int iter=-1);

    void initXX(const char *F, int L, int m=4, int iter=-1) { 
        if (m==1 || m==4) { init(F,L,"KK",iter); } else
        if (m==2 || m>4) wblog(FL,"ERR %s() invalid m=%d !?",FCT,m);
        if (m>2) {
           init(F,L,"KD",iter);
           init(F,L,"DK",iter);
           init(F,L,"DD",iter);
        }
    };

    void initX(const char *F, int L, int m=2, int iter=-1) { 
        if (m==2) { init(F,L,"D",iter); }
        else if (m!=1) wblog(FL,"ERR %s() invalid m=%d !?",FCT,m);
        init(F,L,"K",iter);
    };

    void init(char flag=0) {
       if (!flag) { name.init(); NAME.init(); }

       nrgIdx.init(); calc=0; store=0;
       if (aux) {
          mxDestroyArray(aux);
          aux=NULL; 
       }

       if (MX && locMX) { mxDestroyArray(MX); } 
       MX=NULL; locMX=0;

       KK.init(); DK.init(); K.clearQSpace();
       KD.init(); DD.init(); D.clearQSpace();
    };

    void setupMX(
       const char *F, int L, int N, char store_=1, char locMX_=1
    ){
       if (MX || locMX) wblog(F_L,
          "ERR %s() already got MX=%p (%d,%d) !?",FCT,MX,locMX,store);
       MX=mxCreateStructMatrix(1,N,0,NULL);
       store=store_; locMX=locMX_;
    };

    void PUT( 
       const char *F, int L, const char *vn, const char *ws="caller"
    ){
       if (MX && locMX) {
          if (!vn || !vn[0]) wblog(FL,"ERR %s() invalid usage",FCT);
          int i=mexPutVariable(ws,vn,MX); 
          if (!i) wblog(FL,"I/O putting %s to %s",vn,ws);
          else wblog(FL,"ERR %s() mexPutVariable() returned i=%d !?",FCT,i);
       }
       else wblog(FL,"WRN %s() failed having %p (%d)",FCT,MX,locMX);
    };

    void check_dloc(
       const char *F, int L, const char *fun,
       unsigned dloc, unsigned iter) const;

    template<class TD_>
    int checkOpConsistency(
        const char *F, int L, const NRGData<TQ,TD_> &B) const;

    bool isQSpaceVec(const char *nm, int iter=-1);

    bool checkVec(const char *F, int L,
       const char *XX, unsigned R, const char *istr="");

    bool checkVecVec(const char *XX, unsigned r);

    void forceCalc();
    void getRelevantOps(wbindex &I) const;

    wbvector< QSpace<TQ,TD> >& getOpsXX(const char *t) {
       if (t && strlen(t)==2)
       switch (t[0]) {
          case 'K': return (t[1]=='K' ? KK : KD);
          case 'D': return (t[1]=='K' ? DK : DD);
       }
       wblog(FL,"ERR invalid tag `%s'", t ? t : "(null)");
       return KK; 
    };

    const QSpace<TQ,TD>& getQSpace(const char *t, unsigned i=0) const;
    QSpace<TQ,TD>& getQSpace(const char *t, unsigned i=0) {
       return const_cast<QSpace<TQ,TD>&>(
       static_cast<const NRGData<TQ,TD>*>(this)->getQSpace(t,i));
    }

    const QSpace<TQ,TD>& getQSpace(char t1, char t2, unsigned i=0) const {
       char t[3]={t1,t2,0}; return getQSpace(t,i);
    };
    QSpace<TQ,TD>& getQSpace(char t1, char t2, unsigned i=0) {
       char t[3]={t1,t2,0}; return getQSpace(t,i);
    };

    const QSpace<TQ,TD>& getQSpace(char t) const {
       switch (t) {
          case 'K': return K;
          case 'D': return D;
       }
       wblog(FL,"ERR invalid tag `%c'(%d)",t,t);
       return K; 
    };
    QSpace<TQ,TD>& getQSpace(char t) {
       return const_cast<QSpace<TQ,TD>&>(
       static_cast<const NRGData<TQ,TD>*>(this)->getQSpace(t));
    };

    void info(const char *vstr=NULL) const;
    void dispIdx(const char *F, int L) const;
    void info(const char *F, int L, const char *istr=" * ") const;
    void disp_otype(const char *F=0, int L=0) const;

    const mxArray* getMxInfo(const char *n); 

    template <class TX>
    void getMxInfo(const char *F, int L, const char *v, TX &X) {
       const mxArray *a=getMxInfo(v);
       if (a) {
          X.init(F_L,a);
       }
       else wblog(F_L,
         "ERR failed to read '%s' from %s_info.mat",v,wnrg());
    };

    mxArray* getMxData(   
    const char *nm, int iter=-1, char ionly=0);

    const char* wnrg() { 
       if (NAME.data) return NAME.data;
       return (name.data ? name.data : NRG_EMPTY_STR);
    };

    mxArray* toMx() const;

    void put(const char *vname, const char* ws="base") {
       mxPutAndDestroy(FL, toMx(), vname);
    };

    void put(const char *F, int L, const char *vname, const char* ws="base"
    ){ wblog(F,L,"I/O putting '%s' to %s", vname, ws);
       mxPutAndDestroy(FL, toMx(), vname);
    };

    template<class T>
    int getPara(const char *name, T &x, int iter=-1, char i=0) {
wblog(FL,"ERR %s() check!",FCT);  
       mxArray *a=getMxData(iter);
       return mxGetNumber(a,x);
    }

    bool opExists(char tflag=1);

    void initOp(unsigned n=1, char calc_=1, char store_=0) {
       KK.initDef(n); DK.initDef(n); nrgIdx.initDef(n);
       KD.initDef(n); DD.initDef(n); calc=calc_; store=store_;
    };

    bool initOp(const char *F, int L,
       const wbvector< QSpace<TQ,TD> > &C0,
       char tflag=1,
       NRGData<TQ,TD>* CR=NULL, 
       const char *istr=""
    );

    bool initOp(const char *F, int L,
       const wbvector< wbvector< QSpace<TQ,TD> > > &CC,
       NRGData<TQ,TD>* CR=NULL, 
       const char *istr="",
       char tflag=1
    );

    void initSym(const char *tag);

    void updatePara(
       const char *F, int L,
       const char *tag, mxArray *, int iter=-1
    ) const;

    void updatePara(const char *F, int L,
       const char *tag, int iter=-1,
       const char *tag2=NULL, char all=-1
    ) const;

    void updateOp(const char *F, int L,
       const NRGData<TQ,TD> &A,
       const NRGData<TQ,TD> &B, unsigned nop, int iter=-1
    );

    template<class TC>
    void updateOpProd(
       const char *F, int L, unsigned iter,
       const NRGData<TQ,TD> &B,
       const NRGData<TQ,TC> &C  
    );

    void updateOp(const char *F, int L,
       const NRGData<TQ,TD> &A, unsigned nop, int iter=-1
    ) { updateOp(F,L,A,A,nop,iter); };

    template <class TA> 
    void backPropagateOp(const char *F, int L, 
       const NRGData<TQ,TA> &A0,
       const NRGData<TQ,TA> &A2
    );

    template<class TD_>
    void initBG(
       const char *F, int L,
       const NRGData<TQ,TD_> &B,
       const NRGData<TQ,TD_> &H1, const NRGData<TQ,TD_> &H2,
       int m=3  
    );

    void unsetAllRefs();
    unsigned skipOpZeros(double eps=DEPS, char bflag='b');

    void saveOp(const char*, int, int iter=-1, unsigned m=4) const;

    void times(TD fac, const char *flag);

    void setupIO(const char*, int,
    const char*, const mxArray*, const char *istr=0);

    void updateInfo(
       const char *F, int L,
       const char*, mxArray *, char keep=-1
    ) const;

    void updateInfop( 
       const char *F, int L,
       const char *vtag, mxArray *, char keep=-1
    ) const;

    const char* getFileName(unsigned k, char lenient=0) const {
       if (!lenient && k>=NRG_N) wblog(FL,
       "ERR %s - index out of bounds (%d/%d)", __FUNCTION__, k, NRG_N);
       return getFileName_aux( k);
    };

    wbstring name;

    wbvector< NRGIndex<TQ,TD> > nrgIdx;

    char calc, store; 
    int  kloc;  
    char locMX; 

    wbvector< wbvector< QSpace<TQ,TD> > > CI;

    wbvector< QSpace<TQ,TD> > KK, KD, DK, DD;

    QSpace<TQ,TD> K, D;

    wbvector<char> issym;

    mxArray* MX;   
    wbstring NAME; 

    mxArray *aux;

  private:

    const char* getFileNameI() const {
       return getFileName_aux(-1);
    };

    const char* getFileName_aux(int k) const {
       int i;
       if (k>=0)
            i=snprintf(NRG_FILE,NRG_FLEN,"%s_%02d.mat", NAME.data, k);
       else i=snprintf(NRG_FILE,NRG_FLEN,"%s_info.mat", NAME.data);
       if (i>=NRG_FLEN) wblog(FL,"ERR File name exceeds length (%d/127)",i);
       return NRG_FILE;
    }
};

template<class TQ, class TD>
void NRGData<TQ,TD>::setupIO(
   const char *F, int L, const char *s, const mxArray *a,
   const char *istr 
){
   if (!s || !s[0]) {
      if (!a || (!mxIsCell(a) && !mxIsStruct(a))) wblog(F,L,
         "ERR invalid input NRG data (%s)", a ? mxGetClassName(a):"NULL");
      NAME.init(); NRG_N=0; NRG_ITER=0;
      MX=(mxArray*)a;
   }
   else {
      Wb::matFile f; unsigned n=strlen(s);
      Wb::Clock clio("data:I/O",1,0,&fdmClocks); 

      if (!n || n+16>NRG_FLEN) wblog(FL,
      "ERR invalid file/variable name `%s' (%d/%d)", s, n, NRG_FLEN);

      NRG_N=0; NRG_ITER=0; MX=NULL;
      NAME=s;

#ifdef MATLAB_MEX_FILE
      mxArray *a=mexGetVariable("caller",s); 
      if (a && (mxIsStruct(a) || mxIsCell(a))) {
         if (istr) wblog(F,L,"<i> I/O NRG%s structure variable: %s",istr,s);

         MX=a; 
         return;
      }
#endif

      if (!f.open(F,L, getFileName(0,'l'),"r"))
           wblog(F,L,"ERR cannot open files of type %s_##.mat", s);
      else f.close();

      if (istr) {
         if (strlen(s)<40) {
            wblog(F,L,"<i> NRG%s is %s_##.mat",istr,s);
         }
         else { char s[256];
            if (getcwd(s,256)==NULL) wblog(FL,"ERR failed read cwd");
            wblog(F,L,"<i> using NRG%s data ...",istr);
            printf("\n   pwd: %s\n   I/O: %s_##.mat\n\n",
            Wb::repHome(s).data, Wb::repHome(s).data);
         }
      }
   }
};

template<class TQ, class TD> 
void NRGData<TQ,TD>::init(const char *F, int L, const char *XX, int iter) {

   char s[64]; unsigned n=strlen(XX), e=0;
   mxArray *a;

   if (!n || n>2) wblog(FL,"ERR invalid operator tag `%s'",XX);
   if (iter<0) iter=NRG_ITER;

   strncpy(s,name.data,60); s[60]=0; strcat(s,XX);
   a=getMxData(s,iter);
   if (!a) wblog(FL,"ERR failed to read NRGData `%s' @ iter=%d",s,iter);

   if (n==1) {
      if (!strcmp(XX,"K")) K.init(F,L,a); else
      if (!strcmp(XX,"D")) D.init(F,L,a); else ++e;
   }
   else {
      if (!strcmp(XX,"KK")) mxInitQSpaceVec(F,L,a,KK); else
      if (!strcmp(XX,"KD")) mxInitQSpaceVec(F,L,a,KD); else
      if (!strcmp(XX,"DK")) mxInitQSpaceVec(F,L,a,DK); else
      if (!strcmp(XX,"DD")) mxInitQSpaceVec(F,L,a,DD); else ++e;
   }
   if (e) wblog(FL,"ERR invalid operator tag `%s'", XX);
};

template<class TQ, class TD>
void NRGData<TQ,TD>::init( 
   const char *F, int L, const mxArray* a,
   wbvector< QSpace<TQ,TD> > &C0, 
   int k,             
   char force,        
   const char *vname0 
){
   char nstr[32]; nstr[0]=0;

   Wb::initLocalOp(F,L,a,C0,k,force,nstr);

   if (vname0 && vname0[0]) { name=vname0; }
   else if (nstr[0]) { name=nstr; }

}

template<class TT>
void Wb::initLocalOp(const char *F, int L,
   const mxArray* a, 
   TT &C,            
   int k,            
   char force,       
   char *name        
){
   size_t n=12; char istr[n];
   if (k)
        snprintf(istr,n,"arg #%d",k);
   else strcpy(istr,"operator");

   if (!a || mxIsEmpty(a)) {
      if (force) wblog(F,L,"ERR invalid %s%s",istr,a ? " (empty)":"");
      return;
   }

   if (mxIsChar(a)) {
      const unsigned n=32; wbstring aux;

      if (name==NULL) { aux.init(n+1); name=aux.data; }
      name[0]=0;

      if (mxGetString(a,name,n) || name[0]==0) wblog(F,L,
         "ERR invalid string of %s%s",istr,name[0] ? "":" (empty)");

      a=mexGetVariablePtr("caller",name);

      if (!a) wblog(F,L,
      "ERR Can not read QSpace `%s' (%s) from workspace",name,istr);
   }

   Wb::initLocalOp_aux(F,L,a,C);

   if (force && C.isEmpty()) wblog(F,L,
   "ERR invalid QSpace (%s is empty)",istr);
};

template<class TQ, class TD>
void NRGData<TQ,TD>::forceCalc() {

   wbvector< NRGIndex<TQ,TD> > &I=nrgIdx;
   if (I.len==0) return;

   for (unsigned i=0; i<I.len; ++i) {
      if (I[i].idx<0) continue;
      if (I[i].nrg==NULL) I[i].idx=-1; 
      if (!I[i].tag.isEmpty()) {
         I[i].init(); I[i].idx=-1; 
      }
   }

   calc=store=1;
}

template<class TQ, class TD>
void NRGData<TQ,TD>::updateOp(
   const char *F, int L,
   const NRGData<TQ,TD> &A, const NRGData<TQ,TD> &B,
   unsigned nop, 
   int iter
){
   unsigned i=0,j,nrefs=0, m=nrgIdx.len;
   const char* tags[4] = { "KK","KD","DK","DD" };
   int e=0,k;

   wbvector< QSpace<TQ,TD> > X, X2;

   if (iter<0) iter=NRG_ITER;
   if (nrgIdx.len==0) return;

   for (; i<m; ++i) { if (nrgIdx[i].idx>=0) ++nrefs; }

   unsetAllRefs();

   if (nrefs<nrgIdx.len) { 
   if (nop==4) {
         Wb::updateOp(F,L,*this, DD, A.D, B.D, KK, iter);
         Wb::updateOp(F,L,*this, DK, A.D, B.K, KK, iter);
         Wb::updateOp(F,L,*this, KD, A.K, B.D, KK, iter);
      e+=Wb::updateOp(F,L,*this, KK, A.K, B.K, KK, iter); 
   }
   else if (nop==3) { 
         Wb::updateOp(F,L,*this, DD, A.D, B.D, KK, iter);
         Wb::updateOp(F,L,*this, DK, A.D, B.K, KK, iter);
      e+=Wb::updateOp(F,L,*this, KK, A.K, B.K, KK, iter); 

      if (KD.len!=DK.len) KD.initDef(DK.len);
      for (i=0; i<DK.len; ++i) DK[i].transp(KD[i]);  
   }
   else if (nop==1) {
      e+=Wb::updateOp(F,L,*this, KK, A.K, B.K, KK, iter); 
   }
   else wblog(FL,"ERR %s() invalid nop (%d)",FCT,nop);
   }

   if (e) wblog(FL,
      "ERR %s() operator %s[%d] updates to empty QSpace!",
       FCT,name.data,iter
   );

   if (nrefs) {
   for (i=0; i<nop; ++i) {

      wbvector< QSpace<TQ,TD> > &C = getOpsXX(tags[i]);
      if (C.len!=m) C.initDef(m);

      X.init(); X2.init(); 

      for (j=0; j<m; ++j) {
         k=nrgIdx[j].idx; if (k<0) continue;

         if (nrgIdx[j].nrg) {
            wbvector< QSpace<TQ,TD> > &R=nrgIdx[j].nrg->getOpsXX(tags[i]);

            if (k>=(int)R.len) wblog(FL,
               "ERR index out of bounds (%d/%d)",k,R.len);

            C[j].init2ref(R[k]);
         }
         else if (nrgIdx[j].tag.isEmpty()) {
            if (!X.len) {
               sprintf(str,"%s%s", name.data, tags[i]);
               mxInitQSpaceVecR23(F,L,getMxData(str,iter),X);

               if (k>=int(X.len)) { wblog(FL,
                  "ERR index out of bounds (%d/%d; %d)\n"
                  "(hint: is stored NRG operator data up to date?)",
                  k,X.len,nrgIdx.len);
               }
            }
            if (k>=int(X.len)) wblog(FL, 
               "ERR index out of bounds (%d/%d; %d)",k,X.len,nrgIdx.len);

            X[k].save2(C[j]);
         }
         else {
            if (!X2.len) {
               sprintf(str,"%s%s", nrgIdx[j].tag.data, tags[i]);
               mxInitQSpaceVecR23(F,L, getMxData(str,iter),X2);
            }
            if (k>=(int)X2.len) wblog(FL, 
               "ERR index out of bounds (%d/%d)",k,X2.len);

            X2[k].save2(C[j]);
         }
      }
   }}
};

template <class TQ, class TS>
template<class TD>
void NRGData<TQ,TS>::updateOpProd(
   const char *F, int L, unsigned iter,
   const NRGData<TQ,TS> &B,
   const NRGData<TQ,TD> &C  
){
   unsigned i=0,r,n;

   QSpace<TQ,TS> Xd;

   char ics[13] = "2  2* 23 23*", *ic=0, *ic_=0;
   for (; ics[i]; ++i) { if (ics[i]==' ') ics[i]=0; }

   n=B.checkOpConsistency(0,1,C);
   if (int(n)<0) wblog(FL,
      "ERR got dim/itag inconsistency between B and C !?");
   if (!n || B.KK.len!=n || B.DD.len!=n) wblog(FL,
      "ERR %s() got empty operators (%d/%d/%d) !?",FCT,B.KK.len,B.DD.len,n);
   KK.init(n); KD.init(n); DK.init(n); DD.init(n);

   for (i=0; i<n; ++i) {
      r=B.KK[i].rank(FL); if (!r) r=B.DD[i].rank(FL);
      if (r==2) { ic=ics+0; ic_=ics+3; } else
      if (r==3) { ic=ics+6; ic_=ics+9; } else
      wblog(FL,"ERR %s() invalid rank r=%d (i=%d)",FCT,r,i);

      B.KK[i].contract(F,L,ic,C.KK[i],ic_,KK[i]);
      B.KD[i].contract(F,L,ic,C.KD[i],ic_,Xd); KK[i]+=Xd;

      B.KK[i].contract(F,L,ic,C.DK[i],ic_,KD[i]);
      B.KD[i].contract(F,L,ic,C.DD[i],ic_,Xd); KD[i]+=Xd;

      B.DK[i].contract(F,L,ic,C.KK[i],ic_,DK[i]);
      B.DD[i].contract(F,L,ic,C.KD[i],ic_,Xd); DK[i]+=Xd;

      B.DK[i].contract(F,L,ic,C.DK[i],ic_,DD[i]);
      B.DD[i].contract(F,L,ic,C.DD[i],ic_,Xd); DD[i]+=Xd;

      if (KK[i].itags.len>2) wblog(FL, 
         "ERR %s() KK[%d] '%s' !?",FCT,i+1,IT2STR(KK[i]));
      if (DD[i].itags.len>2) wblog(FL,
         "ERR %s() DD[%d] '%s' !?",FCT,i+1,IT2STR(DD[i]));
   }
};

template <class TQ, class TS, class TD>
void applyBG(const char *F, int L,
   QSpace<TQ,TS> &B,        
   const QSpace<TQ,TD> &H0,
   const QSpace<TQ,TD> &H2  
){
   unsigned i,j,k,m,n;

   wbvector<unsigned> I1,I2;
   wbvector<double> E1,E2;
   wbindex i1h,i2h;

   wbarray<TD> Omega;
   TD *e1,*e2;

   wbcomplex z, Einc = gRX.Einc + gRX.E0, dEinc = Einc - gRX.E2; 

   gRX.check(); 

   double e12, sgm=gRX.sigma/2.;

   if (B.isEmpty()) { return; }
   matchH12(FL,B,H0,H2,gRX.iter, i1h,E1,I1, i2h,E2,I2); 

   for (k=0; k<i1h.len; ++k) {
      wbarray<TS> &Bk = *B.DATA[k];
      TS *xk=Bk.data;

      i=i1h[k]; m=I1[i+1]-I1[i];
      j=i2h[k]; n=I2[j+1]-I2[j];

      if (!Bk.sameSize(m,n)) { Bk.info(FL,"B",k); 
         MXPut(FL,"i").add(Bk,"B").add(m,"m").add(n,"n").add(k+1,"k");
         wblog(FL,"ERR %s() dimension mismatch (%d) !?",FCT,k+1);
      }

      e1=E1.data+I1[i]; e2=E2.data+I2[j];

      for (j=0; j<n; ++j) {
      for (i=0; i<m; ++i) { e12 = e2[j] - e1[i];
         z.r = dEinc.r - e12;
         z.i =-Einc.i;

         if  (!gRX.wrxs) { 
            z/=sqrt(abs(z.i)); 
         }
         else {
wblog(FL,"ERR %s() deprecated using BHeisenberg!!",FCT); 

            z.i -= sgm*(gRX.E0+gRX.E2);
            if (z.i>=0) wblog(FL,"ERR %s() got z=%g%+gi",FCT,z.r,z.i);

            z/=sqrt(abs(z.i)); 
         }

         xk[i+m*j]/=z;  
      }}
   }
};

template <class TQ, class TS>
template <class TD>
void NRGData<TQ,TS>::initBG(
   const char *F, int L,
   const NRGData<TQ,TD> &B,  
   const NRGData<TQ,TD> &H0,
   const NRGData<TQ,TD> &H2, 
   int m
){
   unsigned i;

   KD.initT_(B.KD); for (i=0; i<KD.len; ++i) applyBG(F_L,KD[i], H0.K, H2.D);
   DK.initT_(B.DK); for (i=0; i<DK.len; ++i) applyBG(F_L,DK[i], H0.D, H2.K);
   DD.initT_(B.DD); for (i=0; i<DD.len; ++i) applyBG(F_L,DD[i], H0.D, H2.D);

   if (m==3) return; else
   if (m!=4) wblog(FL,"ERR %s() invalid m=%d !?",FCT,m);

   KK.initT_(B.KK); for (i=0; i<KK.len; ++i) applyBG(F_L,KK[i], H0.K, H2.K);
};

template <class TQ, class TD, class TZ>
void WeightE0_dim1(QSpace<TQ,TZ> &X, const QSpace<TQ,TD> &H, double a) {

   unsigned i,j,k; widx_t m1,m2,d,D;
   double efac,w=1;

   wbMatrix<TQ> Qh,Qx;
   wbindex ih,ix;

   if (X.isEmpty()) { return; }
   if ((i=X.rank())<2 || i>3 || H.rank()!=2) wblog(FL,
      "ERR %s() got rank-%d data (%d)",FCT,i,H.rank());

   if (a<0) { w=0; a=-a; } else
   if (!a) wblog(FL,"ERR %s() invalid w=%g",FCT,a);

   X.getQsub(0,Qx);
   H.getQsub(0,Qh); 

   matchIndex(Qh,Qx,ih,ix,1,&m1,&m2); 
   if (m1 || ix.len!=Qx.dim1) wblog(FL,"ERR %s() "
      "need unique match (%d/%d; %d/%d)",FCT,m1,m2,ix.len,Qx.dim1);
   if (ix.len!=X.DATA.len) wblog(FL,"ERR %s() "
      "severe QSpace inconsistency (%d/%d)!",FCT,ix.len,X.DATA.len);

   for (k=0; k<ix.len; ++k) {
      wbarray<TZ> &xk=*X.DATA[ix[k]];
      const wbarray<TD> &ek=*H.DATA[ih[k]];
      xk.getMatSize(d,D,1);

      if (!ek.isVector()) wblog(FL,
         "ERR %s() non-diagonal H (%d: %s)",FCT,k+1,SSTR(ek));
      if (d!=ek.dim1()) wblog(FL,"ERR %s() "
         "size mismatch (%d: d=%ld/%ld x %ld)",FCT,k+1,d,ek.dim1(),D);

      for (i=0; i<d; ++i) {
         efac=exp(-a*ek[i]*ek[i]); if (!w) { efac=1-efac; }
         efac=sqrt(efac); 
         for (j=0; j<D; ++j) { xk[i+d*j]*=efac; } 
      }
   }
};

template <class TQ, class TS, class TD>
void matchH12(
   const char *F, int L,
   QSpace<TQ,TS> &X, 
   const QSpace<TQ,TD> &H1, const QSpace<TQ,TD> &H2, unsigned iter,
   wbindex &i1h, wbvector<double> &E1, wbvector<unsigned> &I1,
   wbindex &i2h, wbvector<double> &E2, wbvector<unsigned> &I2,
   char redCGR=0 
){
   unsigned i,j,e=0; widx_t m1,m2;
   wbindex i1,i2;
   wbperm P;

   wbvector<unsigned> D1,D2;
   wbMatrix<TQ> Q1h,Q2h,Q1,Q2;

   double Escale = gEScale(iter);

   if (X.isEmpty()) {
      wblog(F_L,"ERR %N empty matrix elements !?");
   }

   if (X.rank(F_L)==3 || !X.CGR.isEmpty()) { X.isConsistent(F_L);
      if (X.CGR.isEmpty()) {
         if (!X.allAbelian()) { wblog(F_L, 
            "ERR %s() expecting CG data for rank-3 data (%s)",
            FCT, X.qtype.toStr().data);
         }
      }
      else {
      for (i=0; i<X.DATA.len; ++i) {
         wbarray<TS> &D=(*X.DATA[i]);

         if (D.SIZE.len==3) {
            if (D.SIZE[2]!=1) wblog(F_L,
               "ERR expecting operator to have dim3=1 (%s)",SSTR(D));
            if (redCGR) { D.SIZE.len=2; } 
         }
         else if (D.SIZE.len!=2) wblog(F_L,
            "ERR invalid rank operator (%s)",SSTR(D));

         if (redCGR) { 
            for (j=0; j<X.CGR.dim2; ++j) {
               if (!X.CGR(i,j).isScalar()) wblog(F_L,
                  "ERR %s() got non-scalar cref (%d,%d): %s",
                  FCT,i+1,j+1,X.CGR(i,j).toStr().data
               );
            }
         }
      }}
   }

   getEData(H1,E1,D1); E1*=Escale; D1.cumsum_(I1,'x'); 
   getEData(H2,E2,D2); E2*=Escale; D2.cumsum_(I2,'x');

   H1.getQsub(0,Q1h); X.getQsub(0,Q1);
   H2.getQsub(0,Q2h); X.getQsub(1,Q2);

   i=matchIndex(Q1h,Q1,i1h,i1,1,&m1,&m2); if (m1 || i1.len!=Q1.dim1) ++e;
   i=matchIndex(Q2h,Q2,i2h,i2,1,&m1,&m2); if (m1 || i2.len!=Q2.dim1) ++e;

   if (e) wblog(F_L,
      "ERR %s() need unique and exact match (%d)!",FCT,e);
   if (i1.len!=X.DATA.len) wblog(F_L,
      "ERR %s() severe QSpace inconsistency (%d %d)!",
       FCT, i1.len, X.DATA.len
   );

   i1.toPerm(P); P.Invert(); i1h.Select(P);
   i2.toPerm(P); P.Invert(); i2h.Select(P);
};

template<class TQ, class TD>
void NRGData<TQ,TD>::unsetAllRefs() {

   unsigned i,j,m=0;
   const char* tags[4] = { "KK","KD","DK","DD" };

   for (i=0; i<4; ++i) {
      wbvector< QSpace<TQ,TD> > &C = getOpsXX(tags[i]);

      for (j=0; j<C.len; ++j) {
         if (nrgIdx[j].idx<0 || !C[j].isref) continue;
         C[j].init(); ++m;
      }
   }

}

template<class TQ, class TD>
unsigned NRGData<TQ,TD>::skipOpZeros(double eps, char bflag) {

   wbvector< QSpace<TQ,TD> > *XX[] = { &KK, &KD, &DK, &DD };
   unsigned i,j,n,r=0;

   for (i=0; i<4; ++i) {
      n=XX[i]->len;
      for (j=0; j<n; ++j)
      r += (*XX[i])[j].SkipZeroData(eps,bflag);
   }

   return r;
}

template<class TQ, class TD>
void NRGData<TQ,TD>::times(TD fac, const char *flag) {

   unsigned i,j,n=(flag ? strlen(flag) : 0);

   if (n==2) {
      if (!strcmp(flag,"XX")) {
         const char* tags[4] = { "KK","KD","DK","DD" };
         for (i=0; i<4; ++i) {
            wbvector< QSpace<TQ,TD> > &X=getOpsXX(tags[i]);
            for (j=0; j<X.len; ++j) X[j]*=fac;
         }
      }
      else {
         wbvector< QSpace<TQ,TD> > &X=getOpsXX(flag);
         for (j=0; j<X.len; ++j) X[j]*=fac;
      }
   }
   else if (n==1) {
      if (flag[0]=='X') {
         getQSpace('K')*=fac;
         getQSpace('D')*=fac;
      }
      else {
         getQSpace(flag[0])*=fac;
      }
   }
   else wblog(FL, "ERR invalid flag `%s'", flag ? flag : "(null)");
};

template<class TQ, class TD>
void NRGData<TQ,TD>::saveOp(
   const char *F, int L, int iter, unsigned m) const {

   unsigned i=0, n=nrgIdx.len, nrefs=0;
   const char* tags[4]={"KK","KD","DK","DD"};

   if (!m || m==2 || m>4) wblog(FL,
      "ERR value out of bounds (m=%d/%d)",m,4);
   if (iter<0) iter=NRG_ITER;

   for (; i<n; ++i) { if (nrgIdx[i].idx>=0) { ++nrefs; }}

   if (store && (!n || nrefs<n)) { i=0;  
      if (m==3) { i=1; m=4; } 
     #if 0
      if (WBLOG_IO) { sprintf(str,
           "%s() iter=%2d: store=%d, nrefs=%d/%ld, i=%d/%d ",
            FCT, iter, store, nrefs,nrgIdx.len, i,m);
         wblog(FL, MX? "MX: %s":"MAT %s",str);
      }
     #endif
      for (; i<m; ++i) { updatePara(F,L,tags[i],iter); }
   }
};

template<class TQ, class TD>
void NRGData<TQ,TD>::updatePara(
   const char *F, int L, const char *tag, mxArray *a, int iter) const {

   sprintf(str,"%s%s", name.data, tag);
   if (!str[0]) wblog(F,L,"ERR invalid name (empty)");

   if (iter<0) iter=NRG_ITER;
   if (MX) {
      if (!NAME.isEmpty())
      mxUpdateField(F,L, MX,str,iter,a); 
   }
   else {
      Wb::Clock clio("data:I/O",1,0,&fdmClocks); 
      matPutVariable(F,L,getFileName(iter),str,a); 
   }
};

template<class TQ, class TD>
void NRGData<TQ,TD>::updatePara(const char *F, int L,
   const char *tag,  
   int iter,
   const char *tag2, 
   char all
) const {

   unsigned n=strlen(tag), e=0;
   mxArray *a=0;

   if (n==1) {
      if (tag[0]=='K') { a=K.toMx(); } else
      if (tag[0]=='D') { a=D.toMx(); } else e|=1;
   }
   else if (n==2) {
      const wbvector< QSpace<TQ,TD> > *XX=0;
      if (!strcmp(tag,"KK")) { XX=&KK; } else
      if (!strcmp(tag,"KD")) { XX=&KD; } else
      if (!strcmp(tag,"DK")) { XX=&DK; } else
      if (!strcmp(tag,"DD")) { XX=&DD; } else e|=2;

      if (!e) {
         if (all>0 || (all<0 && store>1)) { a=XX->toMx(); }
         else {
            wbindex I; getRelevantOps(I);
            a=QSpaceVec2Mx(*XX,I);
         }
      }
   }
   else { e|=4; }
   if (e) wblog(FL,"ERR invalid tag `%s'",tag);

   sprintf(str,"%s%s", name.data, tag2 && tag2[0] ? tag2 : tag);
   if (!str[0] || !name.isName()) wblog(F,L,"ERR invalid name `%s'",str);

   if (iter<0) iter=NRG_ITER;
   if (MX) {
      mxUpdateField(F,L, MX,str,iter,a); 
   }
   else {
      Wb::Clock clio("data:I/O",1,0,&fdmClocks); 
      matPutVariable(F,L,getFileName(iter),str,a); 
   }
};

template<class TQ, class TD>
inline void NRGData<TQ,TD>::getRelevantOps(wbindex &I) const {

   unsigned i,k;
   I.init(nrgIdx.len);

   for (k=i=0; i<nrgIdx.len; ++i) {
      if (nrgIdx[i].idx<0) I[k++]=i;
   }

   I.len=k;
}

template<class TQ, class TD>
const QSpace<TQ,TD>& NRGData<TQ,TD>::getQSpace(
   const char *t, unsigned i
 ) const {

   int d=-1;

   if (!t || !index("KD",t[0]) || !index("KD",t[1]) || t[2])
   wblog(FL,"ERR invalid tag `%s'", t ? t : "(null)");

   if (t[0]=='K') {
      if (t[1]=='K')
           { if (i<KK.len) return KK[i]; else d=KK.len; }
      else { if (i<KD.len) return KD[i]; else d=KD.len; }
   }
   else {
      if (t[1]=='K')
           { if (i<DK.len) return DK[i]; else d=DK.len; }
      else { if (i<DD.len) return DD[i]; else d=DD.len; }
   }

   wblog(FL,"ERR index out of bounds (%s: %d/%d)",t,i,d);
   return K; 
};

template<class TQ, class TD>
void NRGData<TQ,TD>::updateInfo(const char *F, int L,
  const char *vname, mxArray *a, char keep
) const {

   if (!a) wblog(F,L,
   "WRN update info `%s' with NULL mxArray", vname);

   if (MX) {
      mexPutVariable("caller", vname, a);
      if (!keep) mxDestroyArray(a);
   }
   else {
      Wb::Clock clio("data:I/O",1,0,&fdmClocks); 
      matPutVariable(F,L,getFileNameI(),vname,a,keep);
   }
}

template<class TQ, class TD>
void NRGData<TQ,TD>::updateInfop(const char *F, int L,
  const char *vtag, mxArray *a, char keep
) const {

   if (!name.data || !name.data[0]) wblog(F,L,
   "WRN update info `%s' of unnamed operator", vtag);

   sprintf(str,"%s%s",
      name.data ? name.data : "",
      vtag ? vtag : ""
   );

   updateInfo(F,L,str,a,keep);
}

template<class TQ, class TD>
bool NRGData<TQ,TD>::opExists(char tflag) {

   if (tflag) {
      unsigned i, k=(NRG_N>2 ? NRG_N-2 : 0);
      const char* tag[3] = { "KK","KD","DK" };

      for (i=0; i<3; ++i) {
         sprintf(str,"%s%s",name.data,tag[i]);
         if (!isQSpaceVec(str,k)) return 0;
      }
   }
   else {
      sprintf(str,"%sKK",name.data);
      return isQSpaceVec(str,0);
   }

   return 1;
}

template<class TQ, class TD>
bool NRGData<TQ,TD>::initOp(const char *F, int L,
   const wbvector< QSpace<TQ,TD> > &C0,
   char tflag,
   NRGData<TQ,TD>* CR, 
   const char *istr
){
   unsigned i,j,m,iter,e,r, nrefs=0; 
   const char* tag[3] = { "KK","KD","DK" };
   unsigned ntag=sizeof(tag)/sizeof(char*);
   wbvector< QSpace<TQ,TD> > NRG; 
   mxArray *a;

   nrgIdx.initDef(C0.len); kloc=0;

   if (C0.len==0) {
      init(1); return 0;
   }

   if (istr && istr[0]) sprintf(str,"%s ",istr); else str[0]=0;

   if (!name.isName(99)) wblog(F,L, 
      "ERR invalid name for %soperator (%d,`%s')",
       str, nrgIdx.len, name.data
   );

   for (i=0; i<C0.len; ++i) {
      r=-1; 
      if (!C0[i].isConsistent(r) || r<2 || r>4) wblog(F,L,
         "ERR invalid operator (got rank-%d QSpace) !?\n"
         "%s (%s)",r,SHORT_FL,str);
      if (C0[i].QDIM!=C0[0].QDIM) wblog(FL,
         "ERR invalid operator (QDIM=%d, %d) !?\n%s:%d %s",
          C0[i].QDIM, C0[0].QDIM,FL,str);
   }

   KK.initDef(C0);

   for (i=0; i<C0.len; ++i) {
      if (nrgIdx[i].idx>=0) continue;
      for (j=i+1; j<C0.len; ++j)
      if (C0[j]==C0[i]) {
         nrgIdx[j].init(i,this);
         ++nrefs; 
      }
   }

   if (CR) {
      m=CR->KK.len;
      for (i=0; i<C0.len; ++i) {
         if (nrgIdx[i].idx>=0) continue; 
         for (j=0; j<m; ++j) if (C0[i]==CR->KK[j]) {
            nrgIdx[i].init(j,CR);
            ++nrefs; 
            break;
         }
      }
      sprintf(str,"%s%s",CR->name.data,"KK_");
      a=getMxData(str,0);
      if (a && mxIsQSpaceVecR23(FL,a)) { 
         mxInitQSpaceVecR23(F,L,a,NRG); m=NRG.len;
         for (i=0; i<C0.len; ++i) {
            if (nrgIdx[i].idx>=0) continue; 
            for (j=0; j<m; ++j) if (C0[i]==NRG[j]) {
               nrgIdx[i].init(j,CR->name.data);
               ++nrefs; 
               break;
            }
         }
      }

      if (nrefs==nrgIdx.len) {
         return nrefs; 
      }
   }

   e=0;

   if (tflag) {  
      iter=(NRG_N>2 ? NRG_N-2 : 0);
      for (i=0; i<ntag; ++i) {
         sprintf(str,"%s%s",name.data,tag[i]);
         if (!isQSpaceVec(str,iter)) { ++e; break; }
      }
   }

   iter=0; 

   sprintf(str,"%s%s",name.data,"KK");
   if (!isQSpaceVec(str,iter)) ++e;

   sprintf(str,"%s%s",name.data,"KK_");
   a=getMxData(str,iter); if (!a || !mxIsQSpaceVecR23(FL,a)) ++e;

   if (e) {
      calc=store=1;
      return nrefs;
   }

   mxInitQSpaceVecR23(F,L,a,NRG);

   for (i=0; i<C0.len; ++i) {
      if (nrgIdx[i].idx>=0 && nrgIdx[i].tag.isEmpty()) continue;
      for (j=0; j<NRG.len; ++j) if (C0[i]==NRG[j]) {
         nrgIdx[i].tag.init();
         nrgIdx[i].idx=j; ++nrefs; 
         break;
      }
   }

   calc= (nrefs==C0.len ? 0 : 1);
   store=(nrefs==C0.len ? 0 : 1);

   return nrefs;
};

template<class TQ, class TD>
bool NRGData<TQ,TD>::initOp(const char *F, int L,
   const wbvector< wbvector< QSpace<TQ,TD> > > &CC,
   NRGData<TQ,TD>* CR, 
   const char *istr,
   char tflag          
){
   unsigned i,j,l,m,iter,q=0,e,r, nrefs=0;
   const char* tag[3]={ "KK","KD","DK" }; 
   unsigned ntag=sizeof(tag)/sizeof(char*);
   const mxArray *a;

   wbvector< wbvector< QSpace<TQ,TD> > > NRG;

   nrgIdx.initDef(CC.len);

   if (CC.len==0) {
      init(1); return 0;
   }

   if (istr && istr[0])
   sprintf(str,"%s ",istr); else str[0]=0;

   if (!name.isName(99)) wblog(F,L, 
      "ERR invalid name for %soperator (%d,`%s')",
       str, nrgIdx.len, name.data
   );

   kloc=-1;

   for (i=0; i<CC.len; ++i) {
   for (l=j=0; j<CC[i].len; ++j) {
      const QSpace<TQ,TD> &x = CC[i][j]; if (x.isEmpty()) continue;
      r=-1; l=j; if (!q) q=x.QDIM;

      if (!x.isConsistent(r) || r<2 || r>3) wblog(F,L,
         "ERR invalid operator (rank-%d QSpace) !?\n"
         "--> %s %s",r, SHORT_FL, str);
      if (x.QDIM!=q) wblog(FL,
         "ERR invalid operator (QDIM=%d, %d) !?\n--> %s %s",
          x.QDIM, q, SHORT_FL, str
      );
   } kloc=MAX(kloc,(int)l); }

   if (kloc<0) wblog(FL,"ERR all empty operators !?");
   kloc=MAX(0,kloc-1); 

   KK.initDef(CC.len);
   CI.initDef(CC.len);
   for (i=0; i<CC.len; ++i) {
      if (!CC[i].len) wblog(FL,"ERR empty operator set !?");
      CI[i].initDef(CC[i]);
      KK[i]=CC[i][0]; 
   }

   for (i=0; i<CC.len; ++i) { if (nrgIdx[i].idx>=0) continue;
      for (j=i+1; j<CC.len; ++j)
      if (CC[i].isEqual(CC[j])) { ++nrefs;
         nrgIdx[j].init(i,this);
      }
   }

   if (CR && CR->CI.len) { m=CR->CI.len;
      for (i=0; i<CC.len; ++i) {
         if (nrgIdx[i].idx>=0) continue; 

         for (j=0; j<m; ++j)
         if (CC[i].isEqual(CR->CI[j])) {
            nrgIdx[i].init(j,CR);
            ++nrefs; break; 
         }
      }

      sprintf(str,"%s%s",CR->name.data,"CI_");
      a=getMxInfo(str);
      if (a) { if (mxIsQSpaceVecVec(a,2)) {
         mxInitQSpaceVecVec(F,L,a,NRG); m=NRG.len;
         for (i=0; i<CC.len; ++i) {
            if (nrgIdx[i].idx>=0) continue; 
            for (j=0; j<m; ++j) if (CC[i]==NRG[j]) {
               nrgIdx[i].init(j,CR->name.data);
               ++nrefs; break;
            }
         }} else wblog(FL,
         "WRN %s::CI_ not a QSpace cell vector !?", CR->name.data);
      }

      if (nrefs==nrgIdx.len) { return nrefs; } 
   }

   e=0;

   if (tflag && NRG_N>2) {
      iter=NRG_N-2;
      for (i=0; i<ntag; ++i) {
         sprintf(str,"%s%s",name.data,tag[i]);
         if (!isQSpaceVec(str,iter)) { ++e; break; }
      }
   }

   iter=0;
   sprintf(str,"%s%s",name.data,"KK");
   if (!isQSpaceVec(str,iter)) ++e;

   sprintf(str,"%s%s",name.data,"CI_");
   a=getMxData(str,iter); if (!a || !mxIsQSpaceVec(FL,a)) ++e;

   if (e) { calc=store=1;
      wblog(FL," *  calculate all operators %s (%d)",
      name.data, nrgIdx.len); return nrefs;
   }

   mxInitQSpaceVecVec(F,L,a,NRG); m=NRG.len;

   for (i=0; i<CC.len; ++i) {
      if (nrgIdx[i].idx>=0 && nrgIdx[i].tag.isEmpty()) continue;
      for (j=0; j<m; ++j) if (CC[i].isEqual(NRG[j])) {
         nrgIdx[i].tag.init();
         nrgIdx[i].idx=j; ++nrefs;
         break;
      }
   }

   calc  = (nrefs==CC.len) ? 0 : 1;
   store = (nrefs==CC.len) ? 0 : 1;

   return nrefs;
};

template<class TQ, class TD>
void NRGData<TQ,TD>::initSym(const char *tag) {

  if (strlen(tag)==2) {
     const wbvector< QSpace<TQ,TD> > &K=getOpsXX(tag);
     unsigned i,n=K.len;
     issym.init(n);

     for (i=0; i<n; ++i) {
        if (K[i].isHConj()) issym[i]=+1; else
        if (K[i].isAHerm()) issym[i]=-1;
     }
  }
  else if (strlen(tag)==1) {
     const QSpace<TQ,TD> &K=getQSpace(tag);
     issym.init(1);
     if (K.isHConj()) issym[0]=+1; else
     if (K.isAHerm()) issym[0]=-1;
  }
  else wblog(FL,"ERR invalid tag `%s'", tag ? tag : "(null)");
}

template<class TQ, class TD>
void NRGData<TQ,TD>::check_dloc(
   const char *F, int L, const char *fun, unsigned dloc, unsigned iter) const {

   widx_t d;
   K.getDim(2,&d);  

   if (iter+1>=NRG_N) {
      if (d) wblog(F_L,
         "ERR %s() got dloc=%li/%d at k=%d/%d !?",fun,d,dloc,iter,NRG_N);
      return;
   }

   if (d!=dloc && (iter || !D.isEmpty())) {
      char s[64]; snprintf(s,64,
         "got inconsistent dloc=%li/%d (k=%d/%d)",d,dloc,iter,NRG_N);

      if (D.isEmpty())
           wblog(F_L,"WRN %s() %s",fun,s);
      else wblog(F_L,"ERR %s() %s",fun,s);
   }
};

template<class TQ, class TD>
template<class TD_>
int NRGData<TQ,TD>::checkOpConsistency(
   const char *F, int L, const NRGData<TQ,TD_> &B) const {

   unsigned e=0, n = (KK.len ? KK.len : B.KK.len);
   if (!n) n = (DD.len ? DD.len : B.DD.len);
   if (!n) { return n; }

   if (KK.len && B.KK.len) {
     if (KK.len!=B.KK.len) { ++e;
        if (L) wblog(F_L,"WRN ... KK.len=%d/%d",KK.len,B.KK.len); }
     if (KD.len!=B.KD.len || KD.len!=n) { ++e;
        if (L) wblog(F_L,"WRN ... KD.len=%d/%d/%d",KD.len,B.KD.len,n); }
     if (DK.len!=B.DK.len || DK.len!=n) { ++e;
        if (L) wblog(F_L,"WRN ... DK.len=%d/%d/%d",DK.len,B.DK.len,n); }
     if (DD.len!=B.DD.len || DD.len!=n) { ++e;
        if (L) wblog(F_L,"WRN ... DD.len=%d/%d/%d",DD.len,B.DD.len,n); }
   }
   else if (KK.len) {
     if (KD.len!=n) { ++e;
        if (L) wblog(F_L,"WRN ... KD.len=%d/%d/%d",KD.len,B.KD.len,n); }
     if (DK.len!=n) { ++e;
        if (L) wblog(F_L,"WRN ... DK.len=%d/%d/%d",DK.len,B.DK.len,n); }
     if (DD.len!=n) { ++e;
        if (L) wblog(F_L,"WRN ... DD.len=%d/%d/%d",DD.len,B.DD.len,n); }
   }
   else {
     if (B.KD.len!=n) { ++e;
        if (L) wblog(F_L,"WRN ... KD.len=%d/%d/%d",KD.len,B.KD.len,n); }
     if (B.DK.len!=n) { ++e;
        if (L) wblog(F_L,"WRN ... DK.len=%d/%d/%d",DK.len,B.DK.len,n); }
     if (B.DD.len!=n) { ++e;
        if (L) wblog(F_L,"WRN ... DD.len=%d/%d/%d",DD.len,B.DD.len,n); }
   }
   if (e) {
      if (F) wblog(FL,"ERR %s() got size consistency",FCT);
      return -e;
   }

   return n;
};

template<class TQ, class TD>
bool NRGData<TQ,TD>::isQSpaceVec(const char *nm, int iter) {

    if (iter<0) iter=NRG_ITER;
    if (iter>=(int)NRG_N) wblog(FL,
       "ERR Index out of bounds (%d/%d)", iter, NRG_N);

    const mxArray *a=getMxData(nm,iter,'i');   
    if (a && mxIsQSpaceVecR23(FL,a)) return 1; 

    return 0;

}

template<class TQ, class TD>
bool NRGData<TQ,TD>::checkVec(
   const char *F, int L,
   const char *XX, unsigned R, const char *istr
){
   unsigned l,n=128;
   char nm[n]; int fid;

   char cflag=WbUtil<TD>::isComplex();

   l=snprintf(nm,n,"%s%s", name.data, XX); 
   if (l>=n) wblog(FL,"ERR strlen out of bounds (%d/%d)",l,n);

   if (MX) {
      const size_t *s=mxGetDimensions(MX);
      const int n=mxGetNumberOfDimensions(MX);

      fid=mxGetFieldNumber(MX, nm);
      if (fid<0) {
         if (istr) {
            wblog(F,L,"WRN invalid %s", istr);
            wblog(FL, "ERR field `%s' could not be found",nm);
         }
         return 1;
      }

      if (n!=2 || (s[0]!=1 && s[1]!=1)) {
         if (istr) {
            wblog(F,L,"WRN invalid %s (%s)",istr,nm);
            wblog(FL, "ERR NRG data must be row structure vector (%d,%d)",
            (int)s[0],(int)s[1]);
         }
         return 1;
      }

      if (!NRG_N) NRG_N=unsigned(s[0]*s[1]); else
      if (NRG_N!=unsigned(s[0]*s[1])) {
         if (istr) {
            wblog(F,L,"WRN invalid %s (%s)",istr,nm);
            wblog(FL, "ERR size mismatch (%d/%d)", NRG_N, int(s[0]*s[1]));
         }
         return 1;
      }

      return !mxsIsQSpaceVec(MX,fid,R,cflag);
   }
   else {
      int i=mxIsQSpaceVec(F_L,NAME.data,nm,R,0,0,cflag,NRG_N);
      if (i<0) {
         wblog(F,L,"ERR invalid %s (%s)%N%N%s%N",istr? istr:"!?",nm,str);
         return 1;
      }

      if (!NRG_N) NRG_N=i; else
      if (i!=(int)NRG_N) {
         wblog(F,L,"WRN invalid %s (%s)", istr?istr:"!?", nm);
         wblog(FL,"ERR size mismatch (%d,%d)", i, NRG_N);
         return 1;
      }
   }

   return 0;
};

template<class TQ, class TD>
bool NRGData<TQ,TD>::checkVecVec(const char *XX, unsigned R) {

   unsigned l,n=128;
   char nm[n]; int fid;

   l=snprintf(nm,n,"%s%s",name.data,XX); 
   if (l>=n) wblog(FL,"ERR strlen out of bounds (%d/%d)",l,n);

   if (MX) {
      fid=mxGetFieldNumber(MX,nm);
      if (fid<0) {
         sprintf(str, "%s:%d field `%s' could not be found.",
         FL,nm); return 1;
      }

      if (!NRG_N) {
         const size_t *s=mxGetDimensions(MX);
         const int n=mxGetNumberOfDimensions(MX);
         if (n!=2 || (s[0]!=1 && s[1]!=1)) { sprintf(str,
            "%s:%d NRG data must be row structure vector (%d,%d)",
             FL, (int)s[0], (int)s[1]); return 1;
         }
         NRG_N=unsigned(s[0]*s[1]);
      }
      return !mxsIsQSpaceVEC(MX, fid, R);
   }
   else {

      int i=mxIsQSpaceVEC(FL,NAME.data, nm, R);
      if (i<0) return 1;
      if (!NRG_N) NRG_N=i;
      if (i!=(int)NRG_N) {
         sprintf(str,"%s:%d size mismatch (%d,%d)",
         FL, i, NRG_N); return 1;
      }
   }

   return 0;
}

template<class TQ, class TD>
mxArray* NRGData<TQ,TD>::getMxData(const char *vs, int iter, char ionly) {
    if (iter<0) iter=NRG_ITER;
    if (iter>=(int)NRG_N) wblog(FL,
       "ERR %s() index out of bound (%d/%d)",FCT,iter,NRG_N);

    if (MX) {
       int fid=mxGetFieldNumber(MX, vs);
       if (fid<0) return NULL;
       return mxGetFieldByNumber(MX,iter,fid);
    }
    else {
       Wb::matFile F(FL,getFileName(iter),"r");
       Wb::Clock clio("data:I/O",1,0,&fdmClocks); 

       if (aux) { mxDestroyArray(aux); aux=NULL; }
       if (ionly)
            { aux = matGetVariableInfo(F.mfp,vs); }
       else { aux = matGetVariable    (F.mfp,vs); }

       F.close(); return aux;
    }
};

template<class TQ, class TD>
const mxArray* NRGData<TQ,TD>::getMxInfo(const char *v) {

    if (!v || !v[0]) wblog(FL,
       "ERR %s() invalid usage (got %s)",FCT,v? "empty":"null");

    if (aux) { mxDestroyArray(aux); aux=NULL; }

    if (MX) {
       return mexGetVariablePtr("caller",v); 
    }
    else {
       Wb::matFile F(FL,getFileNameI(),"r");
       Wb::Clock clio("data:I/O",1,0,&fdmClocks); 

       aux=matGetVariable(F.mfp, v);  

       return aux;
    }
};

template<class TQ, class TD>
mxArray* NRGData<TQ,TD>::toMx() const {

    MXPut S(0,0);
    S.add(name,"name")
     .add(calc,"calc").add(store,"store")
     .add(kloc,"kloc").add(issym,"issym")
     .addP(MXPut(0,0)
        .add(CI,"CI").add(nrgIdx,"nrgIdx").add(NAME,"NAME")
        .add(Wb::num2Str(MX,"0x%lX"),"MX").add(locMX,"locMX")
        .add(NRG_ITER,"ITER")
        .add(NRG_N,"N").toMx(),"Idx"
     );

    if (KK.len || KD.len || DK.len || DD.len) {
    S.add(KK,"KK").add(KD,"KD").add(DK,"DK").add(DD,"DD"); }

    if (!K.isEmpty() || !D.isEmpty()) {
    S.add(K,"K").add(D,"D"); }

    return S.toMx();
};

template<class TQ, class TD>
void NRGData<TQ,TD>::info(const char *vstr) const {

    unsigned i,nl=0;
    char sx[nrgIdx.len+1], s_[128];

    for (i=0; i<nrgIdx.len; ++i) {
       if ( nrgIdx[i].idx<0) sx[i]='C';
       if (!nrgIdx[i].nrg) sx[i]='R'; 
       else if ( nrgIdx[i].nrg==this)
            sx[i]='1'; 
       else sx[i]='2'; 
    }; sx[i]=0;

    snprintf(s_,128,"%d/%d; calc=%d (%s; %d); store=%d",
    NRG_ITER+1, NRG_N, calc, sx, nrgIdx.len, store);

    if (MX) {
       printf("\n%s Structure `%s' (0x%lX), %s\n",
       vstr && vstr[0] ? vstr : (name.data ? name.data : "ans"),
       NAME.data ? NAME.data : "",MX,s_);
    }
    else {
       printf("\n%s File structure `%s', %s\n",
       vstr && vstr[0] ? vstr : (name.data ? name.data : "ans"),
       NAME.data ? NAME.data : "",s_);
    }

    if (!KK.isEmpty()) {
       if (!nl) { printf("\n"); ++nl; }
       for (i=0; i<KK.len; ++i) {
          sprintf(s_,"KK[%d]", i);
          KK[i].info(s_);
       }
    }

    if (!KD.isEmpty()) {
       if (!nl) { printf("\n"); ++nl; }
       for (i=0; i<KD.len; ++i) {
          sprintf(s_,"KD[%d]", i);
          KD[i].info(s_);
       }
    }
    if (!DK.isEmpty()) {
       if (!nl) { printf("\n"); ++nl; }
       for (i=0; i<DK.len; ++i) {
          sprintf(s_,"DK[%d]", i);
          DK[i].info(s_);
       }
    }
    if (!DD.isEmpty()) {
       if (!nl) { printf("\n"); ++nl; }
       for (i=0; i<DD.len; ++i) {
          sprintf(s_,"DD[%d]", i);
          DD[i].info(s_);
       }
    }

    nl=0;

    if (!K.isEmpty()) { if (!nl) { printf("\n"); ++nl; }; K.info("K"); }
    if (!D.isEmpty()) { if (!nl) { printf("\n"); ++nl; }; D.info("D"); }

    if (nl) printf("\n");
};

template <class TQ, class TD>
void NRGData<TQ,TD>::dispIdx(const char *F, int L) const {

   char s[32];
   wblog(F,L,"%NNRGIndex %s (calc=%d, store=%d)%N", name.data, calc, store);

   for (unsigned i=0; i<nrgIdx.len; ++i) {
      if (!nrgIdx[i].tag.isEmpty())
           strncpy(s, nrgIdx[i].tag.data,32);
      else if (nrgIdx[i].nrg)
           strcpy(s, nrgIdx[i].nrg==this ? "(this)" : "(that)");
      else if (nrgIdx[i].idx>=0)
           strcpy(s, "(file)");
      else strcpy(s, "(calc)");

      printf("%6d: %4d  %s:%lX\n", i, nrgIdx[i].idx, s, nrgIdx[i].nrg);
   }

   if (nrgIdx.len) printf("\n");
}

template <class TQ, class TD>
void NRGData<TQ,TD>::info(const char *F, int L, const char *istr) const {

   unsigned i, nt1=0, nt2=0, nx1=0, nx2=0, nc=0;
   wbstring mark(nrgIdx.len+1);
   wbstring S(128);

   for (i=0; i<nrgIdx.len; ++i) {

      if (nrgIdx[i].nrg && nrgIdx[i].nrg!=this) { ++nx1; mark[i]='x'; } else
      if (!nrgIdx[i].tag.isEmpty()) { ++nx2; mark[i]='X'; } else
      if (nrgIdx[i].nrg==this) { ++nt1; mark[i]='r'; } else
      if (nrgIdx[i].idx>=0) { ++nt2; mark[i]='R'; }
      else { ++nc; mark[i]='*'; }
   }

   if (nc) {
      S.pushf(FL,", %s%d/%d op%s",
      calc ? "calc ":"", nc, KK.len, KK.len!=1 ? "s":"");
   }
   else {
      S.pushf(FL,", %d op%s%s", KK.len, KK.len!=1 ? "s":"",
      calc ? ", calc!?":", all from store or referenced");
   }

   if (store) S.push(FL," and store");
   else if (nc || calc) {
      S.push(FL," (not stored)");
   }

   if (kloc      ) S.pushf(FL,", nloc=%d", kloc+1);

   if (strcmp(istr,name.data)) {
      wblog(F,L," *  %s (%3s): [%s] %s",
      name.data, istr, mark.data, S.data[0] ? S.data+2 : "");
   }
   else {
      wblog(F,L," *  %s [%s] %s",
      name.data, mark.data, S.data[0] ? S.data+2 : "");
   }
};

template <class TQ, class TD>
void NRGData<TQ,TD>::disp_otype(const char *F, int L) const {

   wbstring s(128); unsigned i=0;

   if (name.data) s.cpy(FL,name.data);

   if (KK.len) {
      s.pushf(FL,"%sKK: ",i?"; ":"");
      for (i=0; i<KK.len; ++i) s.pushf(FL,"%d",KK[i].otype);
   }
   if (DK.len) {
      s.pushf(FL,"%sTK: ",i?"; ":"");
      for (i=0; i<DK.len; ++i) s.pushf(FL,"%d",DK[i].otype);
   }
   if (KD.len) {
      s.pushf(FL,"%sKT: ",i?"; ":"");
      for (i=0; i<KD.len; ++i) s.pushf(FL,"%d",KD[i].otype);
   }
   if (DD.len) {
      s.pushf(FL,"%sTT: ",i?"; ":"");
      for (i=0; i<DD.len; ++i) s.pushf(FL,"%d",DD[i].otype);
   }

   wblog(F_L," *  NRGData of type %s() %s ",FCT,s);
};

template <class TQ, class TD, class TZ>
void backPropagateOps(const char *F, int L,
   wbvector< QSpace<TQ,TZ> > &BK,
   const wbvector< QSpace<TQ,TZ> > &BX,
   const QSpace<TQ,TD> &A0,
   const QSpace<TQ,TD> &A2
){
   if (!BX.len) return;
   if (!BK.len) BK.init(BX.len); else
   if (BK.len!=BX.len) wblog(F_L,
      "ERR %s() length mismatch (%d/%d)",FCT,BK.len,BX.len);

   unsigned i=0,r;
   QSpace<TQ,TZ> Q1,Q2;

   for (; i<BX.len; ++i) { r=BX[i].rank(F_L);
      if (!r && !i) {
         for (; i<BX.len; ++i) { if (BX[i].rank(F_L)) break; }
         if (i==BX.len) return; 
      }
      if (r<2 || r>3) wblog(F_L,
         "ERR %s() invalid rank r=%d (%d/%d) !?",FCT,r,i+1,BX.len);

      A2.contract(F_L,"2*",BX[i],2,Q1); 
      A0.contract(F_L,"23",Q1,"32",Q2); 
      BK[i]+=Q2;
   }
};

template <class TQ, class TZ>
template <class TD>
void NRGData<TQ,TZ>::backPropagateOp(const char *F, int L,
   const NRGData<TQ,TD> &A0,
   const NRGData<TQ,TD> &A2
){
   unsigned i=0, e=0;
   wbvector< QSpace<TQ,TZ> > BK;

   A2.check_dloc(F_LF,gRX.dloc,gRX.iter); 

   backPropagateOps(F,L,BK, KK, A0.K, A2.K);
   backPropagateOps(F,L,BK, KD, A0.K, A2.D);
   backPropagateOps(F,L,BK, DK, A0.D, A2.K);
   backPropagateOps(F,L,BK, DD, A0.D, A2.D);

   for (; i<BK.len; ++i) { if (BK[i].isEmpty()) { ++e; }
      BK[i]*=(1./gRX.dloc); 
      BK[i].SkipZeroData(DEPS2,'b');
   }
   if (e && gRX.iter>kloc) {
      MXPut(FL,"Ie").add(*this,"B02").add(gRX.iter,"iter").add(BK,"BK")
       .add(A0,"A0").add(A2,"A2");
      wblog(FL,"ERR %s() got %d/%d empty BK at iter=%d/%d !?",
        FCT,e,BK.len,gRX.iter,kloc);
   }

   BK.save2(KK);
};

template <class TQ, class TD> inline
void updateOp_L(const char *F, int L,
   const QSpace<TQ,TD> &Xin, 
   const QSpace<TQ,TD> &A1,  
   const QSpace<TQ,TD> &A2,
   QSpace<TQ,TD> &Xout       
){
   unsigned rk, r_=-1;
   QSpace<TQ,TD> Xk;

   if (Xin.isEmpty()) { wblog(F_L,
      "ERR %s() got empty L-operator",FCT); return; }
   if (int(rk=Xin.isOperator(&r_))<=0) { Xin.info("L-op");
      wblog(F_L,"ERR %s() got invalid L-operator (%d/%d)",FCT,rk,r_);
   }

   Xin.contract(F,L,2,A2,1,Xk); 
   A1.contract(F,L,"1,3;*",Xk, (rk!=3 ? "1,3":"1,4"), Xout); 

   if (rk==3) Xout.Permute("1,3,2"); 

   Xout.otype=Xin.otype;
   Xout.SkipZeroData();
};

template <class TQ, class TD> inline
void updateOp_s(const char *F, int L,
   const QSpace<TQ,TD> &x,  
   const QSpace<TQ,TD> &A1, 
   const QSpace<TQ,TD> &A2,
   QSpace<TQ,TD> &Xout      
){
   unsigned rk, r_=-1;
   QSpace<TQ,TD> Xk;

   if (x.isEmpty()) { wblog(FL,
      "WRN %s() got empty L-operator",FCT); return; }
   if (int(rk=x.isOperator(&r_))<=0) { x.info("local-op");
      wblog(FL,"ERR %s() got invalid L-operator (%d/%d)",FCT,rk,r_);
   }

   A2.contract(F,L,3,x,2, Xk); 
   A1.contract(F,L,"1,3;*", Xk,"1,3", Xout); 

   Xout.otype=x.otype;
   Xout.SkipZeroData();
};

template <class TQ, class TD>
inline void updateOp_Ls(
   const QSpace<TQ,TD> &Xin, 
   const QSpace<TQ,TD> &A1,  
   const QSpace<TQ,TD> &A2,
   QSpace<TQ,TD> &Xout       
){
   unsigned rk, r_=-1;
   QSpace<TQ,TD> Xk;

   if (Xin.isEmpty()) { wblog(FL,
      "WRN %s() got empty L-operator",FCT); return; }
   if ((rk=Xin.isOperator(&r_,'x'))!=4) { Xin.info("Ls-op");
      wblog(FL,"ERR %s() got invalid Ls-operator (%d/%d)",FCT,rk,r_);
   }

   Xin.contract("3,4",A2,"1,3", Xk);  
   A1.contract("1,3;*",Xk,"1,2",Xout);

   Xout.otype=Xin.otype;
   Xout.SkipZeroData();
};

template <class TQ, class TD>
inline void updateOp_sL(
   const QSpace<TQ,TD> &Xin, 
   const QSpace<TQ,TD> &A1,  
   const QSpace<TQ,TD> &A2,
   QSpace<TQ,TD> &Xout
){
   unsigned rk, r_=-1;
   QSpace<TQ,TD> Xk;

   if (Xin.isEmpty()) { wblog(FL,
      "WRN %s() got empty L-operator",FCT); return; }
   if ((rk=Xin.isOperator(&r_))!=4) { Xin.info("sL-op");
      wblog(FL,"ERR %s() got invalid sL-operator (%d)",FCT,rk,r_);
   }

   Xin.contract("4,3",A2,"1,3", Xk);  
   A1.contract("1,3;*",Xk,"2,1",Xout);

   Xout.otype=Xin.otype;
   Xout.SkipZeroData();
};

template <class TQ, class TD> inline
void updateOp_loc(const char *F, int L,
   const NRGData<TQ,TD> &NRG, unsigned i, 
   const QSpace<TQ,TD> &XL, 
   const QSpace<TQ,TD> &x,  
   const QSpace<TQ,TD> &A1,
   const QSpace<TQ,TD> &A2,
   QSpace<TQ,TD> &Xout,     
   unsigned iter
){
   bool o1=0, o2=0;
   unsigned rk, r_=-1;
   QSpace<TQ,TD> Ax,XA; 

   if (!x.isEmpty()) {
      if (int(rk=x.isOperator(&r_))<=0) { x.info("local_op"); wblog(FL,
         "ERR %s() got invalid local operator (%d/%d)",FCT,rk,r_); }
      if (XL.otype!=x.otype && !XL.isEmpty()) wblog(FL,
         "WRN %s() got different operator types (%s; %s)",
          FCT, QS_STR[XL.otype],QS_STR[x.otype]);
      o2=(rk==3);

      wblog(FL," *  contracting %s[%d] onto A(%d).%s",
         NRG.name.data,i, iter, A2.itags.toStrk(str,8,2,"s"));

      A2.contract(F,L,3,x,"-op",Ax); 
   }
   else Ax.init2ref(A2);

   if (iter || !XL.isEmpty()) { r_=-1;
      if (int(rk=XL.isOperator(&r_))<=0) { x.info("L-op"); wblog(FL,
         "ERR %s() got invalid L-operator (%d/%d)",FCT,rk,r_); }
      o1=(rk==3);

      if (iter<1) wblog(FL,"TST contracting onto A[%d].%s",
         iter, A2.itags.toStrk(str,8,0,"L"));

      if (o1 && o2) { o1=o2=0; 
           XL.contract(F,L,"23",Ax,"14",XA); } 
      else {
         XL.contract(F,L,"-op",Ax,1,XA); 
      }
   }
   else XA.init2ref(Ax);

   if (o1)
        A1.contract(F,L,"13*",XA,"14", Xout,"132"); 
   else A1.contract(F,L,"13*",XA,"13", Xout); 

   if (o1) Xout.otype=XL.otype; else
   if (o2) Xout.otype=x .otype;

   Xout.SkipZeroData();
};

template <class TQ, class TD>
int Wb::updateOp(const char *F, int L,
   const NRGData<TQ,TD> &NRG,
   wbvector< QSpace<TQ,TD> > &F12, 
   const QSpace<TQ,TD> &A1,
   const QSpace<TQ,TD> &A2,
   const wbvector< QSpace<TQ,TD> > &FKK, 
   unsigned iter
){
   unsigned i, j=0, m=0, r_=-1, rk, nop=FKK.len;
   int e=0, i1=0;

   const wbvector< wbvector<QSpace<TQ,TD> > > &CI=NRG.CI;
   const wbvector< NRGIndex<TQ,TD> > &nrgIdx=NRG.nrgIdx;

   if (A1.isEmpty() || A2.isEmpty()) { F12.initDef(nop); return 0; }

   if (nrgIdx.len!=nop) wblog(FL,
      "ERR %s() severe size mismatch (nrgIdx: %d,%d)",FCT,nrgIdx.len,nop);

   if (iter>0) {
      for (i=0; i<nop; ++i) { 
         if (nrgIdx[i].idx>=0) continue;
         if (!(i1=FKK[i].isConsistent(FL))
         || int(rk=FKK[i].isOperator(&r_))<=0) break;
      }
   }
   else if (CI.isEmpty()) {
      for (i=0; i<nop; ++i) {
         if (!(i1=FKK[i].isConsistent(FL)) ||
            int(rk=FKK[i].isOperator(&r_,'x'))<=0
         ){ wblog(FL,
            "WRN FKK[%d] -> (%d; %d/%d/%d) ",i+1,i1,rk,r_,FKK[i].rank());
            break;
         }
      }
   }
   else {
      if (nop!=CI.len) wblog(FL,
         "ERR size inconsistency (%d/%d)",nop,CI.len);

      for (i=0; i<nop; ++i) { m=CI[i].len;
      for (j=0; j<m; ++j) { if (!CI[i][j].isEmpty()) {
         if (!(i1=CI[i][j].isConsistent(FL)) ||
            int(rk=CI[i][j].isOperator(&r_))<=0) break;
      }}}
   }

   if (i<nop || j<m) wblog(F,L,
      "ERR invalid rank-%d operator (iter=%d, op=%d/%d)",rk,iter,i+1,nop);
   if (A1.rank(F,L)!=3 || A2.rank(F,L)!=3 || A1.QDIM!=A2.QDIM) wblog(F,L,
      "ERR %s() severe rank inconsistency (%s; %s)",
       FCT, A1.sizeStrQ().data, A2.sizeStrQ().data);

   if (F12.len!=nop) 
   F12.initDef(nop);

   if ((int)iter>NRG.kloc) {
      for (i=0; i<nop; ++i) {
         if (nrgIdx[i].idx>=0) continue;
         updateOp_L(F,L, FKK[i], A1, A2, F12[i] );
         wbtop().runningLarge(FL);

         if (F12[i].isEmpty() && !FKK[i].isEmpty()) ++e;
      }
   }
   else if (!CI.isEmpty()) {
      QSpace<TQ,TD> E; 

      for (i=0; i<nop; ++i) { if (nrgIdx[i].idx>=0) continue;
         try { updateOp_loc(F,L,NRG,i,
            FKK[i], iter+2>CI[i].len ? E : CI[i][iter+1],
            A1, A2, F12[i],iter);
         }
         catch (...) { wblog_ERR_updateOp(FL,"loc",i,nop,iter); }

         if (F12[i].isEmpty() && !FKK[i].isEmpty()) ++e;
      }
   }
   else {
      char Lflag = A1.getDim(1)>1 || A2.getDim(1)>1;
      unsigned l=0, tlen=16; char tag[tlen];

      wbstring mark(nop,'*');

      static time_t last_call_iter0=0; 
      time_t tnow=time(0);

      for (i=0; i<nop; ++i) {
         if (nrgIdx[i].idx>=0) continue;
         rk=FKK[i].rank(FL);

         if (rk==2 || (rk==3 && FKK[i].otype==QS_OPERATOR)) {
            if (Lflag && 
                FKK[i].hasQOverlap(1,A1,1,'<')>0 &&
                FKK[i].hasQOverlap(2,A1,1,'<')>0
            ){
               try { updateOp_L(F,L,FKK[i],A1,A2, F12[i]); }
               catch (...) { wblog_ERR_updateOp(FL,"L",i,nop,iter); }
               mark[i]='L';
            }
            else
            if (FKK[i].hasQOverlap(1,A1,3,'<')>0 &&
                FKK[i].hasQOverlap(2,A1,3,'<')>0
            ){
               try { updateOp_s(F,L,FKK[i], A1, A2, F12[i]); }
               catch (...) { wblog_ERR_updateOp(FL,"s",i,nop,iter); }
               mark[i]='s';
            }
            else
            if (FKK[i].hasQOverlap(1,A1,2,'<')>0 &&
                FKK[i].hasQOverlap(2,A1,2,'<')>0
            ){
               mark[i]='R'; 
            } 
            else wblog(FL,
              "ERR %s() failed to deal input operator with A%d",FCT,iter);
         }
         else if (rk==4) {

            wbvector<int> m1,m2;
            int mm[]= {
                FKK[i].hasQOverlap(1,A1,1,'<')>0, 
                FKK[i].hasQOverlap(2,A2,3,'<')>0,
                FKK[i].hasQOverlap(3,A2,1,'<')>0,
                FKK[i].hasQOverlap(4,A2,3,'<')>0,

                FKK[i].hasQOverlap(1,A1,3,'<')>0, 
                FKK[i].hasQOverlap(2,A2,1,'<')>0,
                FKK[i].hasQOverlap(3,A2,3,'<')>0,
                FKK[i].hasQOverlap(4,A2,1,'<')>0
            };
            m1.init(4,mm); m2.init(4,mm+4);

            if (m1.allUnequal(0)) {
               try { updateOp_Ls(FKK[i], A1, A2, F12[i]); }
               catch (...) { wblog_ERR_updateOp(FL,"Ls",i,nop,iter); }
               mark[i]='D';
            }
            else
            if (m2.allUnequal(0)) {
               try { updateOp_sL(FKK[i], A1, A2, F12[i]); }
               catch (...) { wblog_ERR_updateOp(FL,"sL",i,nop,iter); }
               mark[i]='d';
            }
            else {
               if (STRICT_ITER0) { mark[i]='?'; 
               if (mark.count('?')==1 && last_call_iter0+5<tnow){
               wblog(FL,
                 "WRN %NFailed to fully associate rank-4 tensor product\n"
                 "operator with first site. [ %s; %s ]\n"
                 "A1: [ %s ];  A2: [ %s ]%N%N"
                 "NB! Hint: A0 may not be complete due to truncation "
                 "in NRGWilsonQS already at iter=0 !?%N",
                  m1.toStr().data, m2.toStr().data,
                  SSTR(A1.getDim()), SSTR(A2.getDim()));
               }}

               try { updateOp_Ls(FKK[i], A1, A2, F12[i]); }
               catch (...) { wblog_ERR_updateOp(FL,"Ls",i,nop,iter); }
            }
         }
         else wblog(FL,
        "ERR %s() invalid rank-%d operator (iter=%d)",FCT,rk,iter);

         if (F12[i].isEmpty() && !FKK[i].isEmpty()) ++e;
      }

      l=snprintf(tag,tlen,"FDM %2d/%d",iter,NRG_N); {
         if (NRG.name.data && NRG.name.data[0] && l<tlen)
            l+=snprintf(tag+l,tlen-l," %s:",NRG.name.data);
         if (l>=tlen) wblog(FL,
        "ERR %s() string out of bounds (%d/%d)",FCT,l,tlen);
      }

      if (!nop)
         wblog(FL,"%s got empty spectral ops",tag);
      else {
         unsigned i1=0, m=mark.count('*');
         for (i=0; i<nop; ++i) { if (mark[i]!='*') { i1=i; break; }}
         for (; i<nop; ++i) { if (mark[i]!='*') {
            if (mark[i]!=mark[i1]) break;
         }}

         if (m==nop) { wblog(FL,
            "%s no spectral ops to contract (all referenced)",
            tag,mark.data);
         }
         else if (i==nop && !m) { 
            char ostr[32];
            if (nop>1)
                 sprintf(ostr,"all %d spectral ops",nop);
            else sprintf(ostr,"spectral operator");

            if (mark[0]!='R') wblog(FL,
               "%s contracting %s onto %c",tag,ostr,mark[0]);
            else wblog(FL,
               "%s %s already assumed in %c",tag,ostr,mark[0]
            );
         }
         else if (i==nop && m) { 
            if (mark[i1]!='R') wblog(FL,"%s contracting %d/%d "
               "spectral ops onto [%s]",tag,nop-m,nop,mark.data);
            else wblog(FL,"%s all %d/%d spectral "
               "ops already assumed in %s",tag,nop-m,nop,mark.data);
         }
         else { 
            wblog(FL,"%s contracting "
            "%d/%d spectral ops onto [%s]",tag,nop-m,nop,mark.data);
         }
      }

      last_call_iter0=tnow;
   }

   if (iter==0) {
      for (i=0; i<F12.len; ++i) {
         if (F12[i].itags.len>2 && !F12[i].itags[2]) {
            F12[i].itags[2].init("op");
         }
      }
   }

   wbtop().runningLarge(FL);
   return e;
};

template <class TQ, class TD>
unsigned applyZ0(const QSpace<TQ,TD> &Z0, QSpace<TQ,TD> &CK){

   unsigned r=CK.rank(FL);
   const QS_TYPES otype=CK.otype;

   if (r==2 || r==3) {
      if (r==3 && CK.otype!=QS_OPERATOR) wblog(FL,
         "ERR %s() expecting otype==%s for rank-%d op (%s)",
          FCT,QS_STR[QS_OPERATOR],r,QS_STR[CK.otype]);
      if (Z0.rank(FL)!=2) wblog(FL,
         "ERR %s() Z0 of rank-2 required for rank-%d op (%d)",
          FCT,r,Z0.rank(FL)
      );

      if (Z0.hasQOverlap(2, CK,1,'>')<=0) {
         wblog(FL,"ERR %s() Z0 does not match operator space",FCT);
         return 1;
      }
      Z0.contract(2,CK,1,CK); 
      CK.otype=otype;
   }
   else if (r==4) {
      if (Z0.rank(FL)!=4) wblog(FL,
         "ERR tensor Z0 required for tensor operator (rank=%d/%d)",
          Z0.rank(FL), r
      );

      if (Z0.hasQOverlap(3, CK,1,'>')<=0 ||
          Z0.hasQOverlap(4, CK,2,'>')<=0) {
         wblog(FL,"ERR %s() Z0 does not match operator space",FCT);
         return 2;
      }
      Z0.contract("3 4",CK,"1 2",CK); 
      CK.otype=otype;
   }
   else wblog(FL,"ERR %s() invalid rank=%d.",FCT,r);

   return 0;
};

template <class TQ, class TD>
unsigned getDLoc(
   NRGData<TQ,TD> &A, Wb::SigHandler *sig=NULL, unsigned *Dk=NULL,
   char fflag=0
){
   unsigned iter, N=NRG_N; widx_t d0=0, d=0;
   char cgflag=-1;

   if (Dk) { (*Dk)=0; }; printf("\r");

   for (iter=0; iter<N; ++iter) { if (sig) sig->check911();
      printf("  %s: %2d/%d loading AK ...   \r",FCT,iter,N);
      fflush(0);

      A.init(FL,"K",iter);

      if (iter==0) cgflag=A.K.gotCGS(FL);

      if (A.K.isEmpty()) {
         if (iter==N-1) continue; else wblog(FL,
         "ERR %s() empty AK at intermediate iteration %d/%d",FCT,iter,N);
      }
      else {
         unsigned r=-1;
         A.K.isConsistent(FL,r);
         if (cgflag<=1 && !A.K.QIDX.isUnique()) wblog(FL,
            "ERR invalid AK at iteration %d/%d",iter,N);
      }

    #ifndef WB_SKIP_ASSERT
      if (fflag) { QSpace<TQ,TD> E;
         A.K.contract("1,3;*",A.K,"1,3",E); 
         if (!E.isIdentityMatrix(1E-12)) {
            MXPut(FL).add(A,"A").add(E,"E");
            wblog(FL,"ERR NRG[%d].AK not in LRs index order",iter,N);
         }
      }
    #endif

      A.K.getDim(2,&d); 

      if (iter>2 || d0) {
         if (d!=d0) {
            MXPut(FL,"a").add(A.K,"K").add(d,"d").add(d0,"d0")
              .add(iter,"iter").add(N,"N");
            wblog(FL,"ERR %s() inconsistent local dimension "
              "(d=%d/%d @ k=%d/%d)",FCT,d,d0,iter,N
            );
         }
      }
      else if (iter==2 || iter+2==N) { d0=d; }
      else {
          A.init(FL,"D",iter);
          if (!A.D.isEmpty()) { d0=d; }
      }

      if (Dk) { 
         d=A.K.getDim(1);  
         if (*Dk<d) { *Dk=d; continue; }
      }
      if (!fflag && iter>2) break; 
   }

   if (!d0) wblog(FL,"ERR %s() "
      "failed to get local dimension (d=%d; %d/%d)",FCT,d0,iter,N);
   return d0;
};

#endif

