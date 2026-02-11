/* ---------------------------------------------------------------------
 * Project : QSpace tensor library (v4.0)
 * Class   : QSpace
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

#ifndef __WB_QSPACE_COL_MAJOR_HH__
#define __WB_QSPACE_COL_MAJOR_HH__

/* ----------------------------------------------------------------- //
   QSpace - set of objects (possibly complex matrices)
   organized by generalized index set storeed as row vectors (recs)
   in index matrix. // Wb,Aug08,05  Apr07,06

   Introduction of non-abelian QSpace through variable 'qtype'.
   see also clebsch.hh // Wb,Sep25,09

   NB! QSpace currently only implemented on linux system
   (matlab+gcc)

   Trying to get it run on a macbook // Wb,Jan18,19
   Problems encountered:
    * OSX, by default, has case-*IN*sensitive files system (APFS)
      encountered troubel with wbcmat.hh + wbCMat.hh
    * default "gcc" compiler actually is with Xcode
      which is much more picky concerning header files!

// ----------------------------------------------------------------- */

#ifdef WB_CLOCK
#define WBC_QSPACE_IO
#endif

enum QS_TYPES {
   QS_NONE, QS_OPERATOR, QS_AMATRIX,
   QS_NUM_TYPES
};

const char* QS_STR[QS_NUM_TYPES] = {
  "", "operator", "A-matrix",
};

template <class TQ, class TD>
class QSpace { 

  public:

    QSpace(unsigned n=0)
     : QDIM(n), otype(QS_NONE), mt(Wb::MEM_DEF), isref(0), ctime(0.) {};

    QSpace(unsigned m, unsigned n, unsigned nq)
     : QSpace(nq) { init(m,n,nq); };

    QSpace(const QSpace &B)
     : QIDX(B.QIDX), CGR(B.CGR), qtype(B.qtype), QDIM(B.QDIM),
       otype(QS_NONE), itags(B.itags), fdir(B.fdir), mt(Wb::MEM_DEF),
       isref(0), ctime(B.ctime)
    {
       if (this!=&B) {
          unsigned i=0; setupDATA();
          for (; i<B.DATA.len; ++i) { DATA[i]->init(*(B.DATA[i])); }
       }
    };

    template<class TQ_, class TD_>
    QSpace(const QSpace<TQ_,TD_> &B,
       char xflag = WbUtil<TD>::isComplex() ? 1 : 0 
    ): QIDX(B.QIDX), CGR(B.CGR), qtype(B.qtype), QDIM(B.QDIM),
       otype(QS_NONE), itags(B.itags), fdir(B.fdir), mt(Wb::MEM_DEF),
       isref(0), ctime(B.ctime)
    {
       if (!xflag) wblog(FL,"ERR %s() verify type cast %s -> %s",
          FCT,TSTR(TD_),TSTR(TD));
       setupDATA();
       for (unsigned i=0; i<B.DATA.len; ++i) {
           DATA[i]->initT(*(B.DATA[i]));
       }
    };

    QSpace( const char *F, int L, const mxArray *a,
       char ref=0, char skip_empty=1, char check_cgw_norm=1)
     : otype(QS_NONE), mt(Wb::MEM_DEF), isref(0), ctime(0.)
     { init(F,L,a,ref,skip_empty,check_cgw_norm); };

    QSpace(const mxArray *a,
       char ref=0, char skip_empty=1, char check_cgw_norm=1)
     : otype(QS_NONE), mt(Wb::MEM_DEF), isref(0), ctime(0.)
     { init(FL,a,ref,skip_empty,check_cgw_norm); };

   ~QSpace() { clearQSpace(); };

    QSpace& clearQSpace(char mflag=0) {
       if (QIDX.data || DATA.data) { mt=Wb::MEM_DEF; } 
       QIDX.init(); qtype.init(); itags.init(); fdir.init(); 
       clearDATA(mflag); CGR.init(); isref=0; ctime=0.;
       return *this;
    };

    void clearDATA(char mflag=0) { 
       if (DATA.len) {
          if (!isref && !DATA.isref) { 
             for (unsigned i=0; i<DATA.len; ++i) { if (DATA[i]) {
                if (mflag) { DATA[i]->initX(); } 
                WB_DELETE_1(DATA[i]);
             }}
          }
          DATA.init();
       }
    };

    unsigned rank(const char *F=nullptr, int L=0) const {
       unsigned r=0; 
       if (QDIM) { r=QIDX.dim2/QDIM;
          if (QIDX.dim2!=r*QDIM) wblog(F_L,
          "ERR %s() QSpace inconsistency (QIDX: %dx%d @ %d)",
          FCT,QIDX.dim1,QIDX.dim2,QDIM);
       }
       else if (QIDX.dim2) wblog(F_L,
          "ERR %s() QSpace inconsistency (QIDX: %dx%d @ %d)",
           FCT,QIDX.dim1,QIDX.dim2,QDIM);
       return r;
    };

    unsigned len() const {
       if (QIDX.dim1!=DATA.len) wblog(FL, 
          "ERR qspace inconsistency (QIDX.dim1=%d, len=%d).",
           QIDX.dim1, DATA.len);
       return QIDX.dim1;
    };

    unsigned nsym() const { 
       if (CGR.dim2 && qtype.len!=CGR.dim2) wblog(FL,
          "ERR qspace inconsistency (nsym=%d '%s'; %s !?)",
           qtype.len, STR(qtype), SSTR(CGR));
       return (qtype.len ? qtype.len : QDIM);
    };

    TQ* Qref(unsigned r, unsigned k) { unsigned l=k*QDIM;
        if (!QDIM || r>=QIDX.dim1 || (l+QDIM)>QIDX.dim2) wblog(FL,
           "ERR %s() index out of bounds (%d/%d; %d*%d/%d)",
            FCT,r,QIDX.dim1,k,QDIM,QIDX.dim2
        );
        return (QIDX.data + r*QIDX.dim2 + l);
    };

    const TQ* Qref(unsigned r, unsigned k) const { 
        unsigned l=k*QDIM;
        if (!QDIM || r>=QIDX.dim1 || (l+QDIM)>QIDX.dim2) wblog(FL,
           "ERR %s() index out of bounds (%d/%d; %d*%d/%d)",
            FCT,r,QIDX.dim1,k,QDIM,QIDX.dim2
        );
        return (QIDX.data + r*QIDX.dim2 + l);
    };

    void init() { clearQSpace(); };

    void init(unsigned m, unsigned n, unsigned nq) { 
       clearQSpace(); 
       QDIM=nq; QIDX.init(m,nq*n);
       setupDATA();
    };

    void init(const wbMatrix<TQ> &Q, unsigned qdim, const QVec *qt=nullptr) {
       clearQSpace(); 
       QIDX=Q; QDIM=qdim; if (qt) qtype=(*qt); else qtype.init();
       setupDATA(FL);
    };

    QSpace& init(const QSpace &B, char ref=0) { 
        initQT_safe(B);
        for (unsigned i=0; i<B.DATA.len; ++i) {
           DATA[i]->init(*(B.DATA[i]),ref); 
        }
        return *this;
    };

    template<class TB> 
    QSpace& initT(const QSpace<TQ,TB> &B) { 
        initQT_safe(B);
        for (unsigned i=0; i<B.DATA.len; ++i) {
           DATA[i]->initT(*(B.DATA[i]));
        }; return *this;
    };

    QSpace& init(const CDATA_TQ &Cb, wbperm *P=nullptr, char inv=0);

    template<class TR> 
    QSpace& getReal(QSpace<TQ,TR> &R) const { 
        unsigned i=0; R.initQT_safe(*this);
        for (; i<DATA.len; ++i) { DATA[i]->getReal(*(R.DATA[i])); }
        return R;
    };

    QSpace<TQ,double> getReal() const { 
        QSpace<TQ,double> R; unsigned i=0; R.initQT_safe(*this); 
        for (; i<DATA.len; ++i) { DATA[i]->getReal(*(R.DATA[i])); }
        return R;
    };

    template<class TI> 
    QSpace& getImag(QSpace<TQ,TI> &I) const { 
        unsigned i=0; I.initQT_safe(*this);
        for (; i<DATA.len; ++i) { DATA[i]->getImag(*(I.DATA[i])); }
        return I;
    };

    QSpace<TQ,double> getImag() const { 
        QSpace<TQ,double> I; unsigned i=0; I.initQT_safe(*this); 
        for (; i<DATA.len; ++i) { DATA[i]->getImag(*(I.DATA[i])); }
        return I;
    };

    void init(const char *F, int L,
       const mxArray *S, char ref, unsigned k,
       char skip_empty=1,
       char check_cgw_norm=1);

    void init(const char *F, int L, const mxArray *S,
         char ref=0, char skip_empty=1, char check_cgw_norm=1
     ) {
        unsigned n = (S ? mxGetNumberOfElements(S) : 0);
        if (n!=1) {
           if (!n) wblog(FL,"ERR %s() got empty array",FCT);
           else    wblog(FL,"ERR %s() got QSpace array (%d entries)",FCT,n);
        }
        init(F,L,S,ref,0,skip_empty,check_cgw_norm);
    };

    QSpace& setupDATA(const char *F=nullptr, int L=0){

       if (!QIDX.dim1) { clearDATA(); return *this; }

       if (QIDX.dim2 && (!QDIM || QIDX.dim2%QDIM)) wblog(F_L,
          "ERR invalid QIDX (%dx%d; %d)",QIDX.dim1,QIDX.dim2,QDIM);
       clearDATA();

       if (isref>1) wblog(FL,"ERR %s() got isref=%d",FCT,isref);

       DATA.init(QIDX.dim1);

       if (!isref) {
          if (mt) { if (mt==Wb::MEX_RETURN) {
             if (WBLOG_MMEX) wblog(FL,
             "TST QSpace::init() to %s",Wb::MTYPE_STR[mt]); }
             else wblog(FL,
             "ERR %s() unexpected mtype=%s",FCT,Wb::MTYPE_STR[mt]);
          }

          for (unsigned i=0; i<DATA.len; ++i) {
              WB_NEW_1(DATA[i]);  
              DATA[i]->mtype=mt;  
          }
       }

       return *this;
    };

    int gotCGS(const char *F=nullptr, int L=0, char xflag=0) const;

    int gotCGX(const char *F=nullptr, int L=0) const;

    QSpace& Enlarge(unsigned m) {
       if (!isConsistent(FL) || !QIDX.dim1) wblog(FL,
          "ERR %s() can't enlarge empty QSpace (%d)",FCT,QIDX.dim1);
       if (!m) return *this;

       unsigned i, n0=QIDX.dim1, n=n0+m;
       gotCGS(FL);

       QIDX.Resize(n, QIDX.dim2);
       DATA.Resize(n);
       for (i=n0; i<n; ++i) { WB_NEW_1(DATA[i]); }

       CGR.Resize(n,CGR.dim2);
       return *this;
    };

    QSpace& Trim(unsigned n) {
       unsigned i, r=-1, n2=QIDX.dim1;
       if (!isConsistent(FL,r,0) || n>QIDX.dim1) wblog(FL,
          "ERR %s() can't trim to LARGER QSpace (%d/%d)",FCT,n,QIDX.dim1);
       if (n==QIDX.dim1) { return *this; }

       gotCGS(FL);

       QIDX.Resize(n, QIDX.dim2);

       for (i=n; i<n2; ++i) { WB_DELETE_1(DATA[i]); }
       DATA.Resize(n);

       CGR.Resize(n,CGR.dim2);

       return *this;
    };

    QSpace& initScalar(const TD &x) { 
        clearQSpace(); DATA.init(1);
        QIDX.init(1,0); 

        WB_NEW_2(DATA[0],wbarray<TD>(1,1));
        DATA[0]->data[0]=x;

        return *this;
    };

    void init2ref(const mxArray *S){
       init(FL,S,'r',0);
    };

    void init2ref(const QSpace &B); 

    void initOpZ_WET( 
       const char *F, int L, const QVec &qvec,
       const wbMatrix<TQ> &Q1, const wbMatrix<TQ> &Q2,
       const wbMatrix<TQ> &Q, const wbarray<TD> &D3,
       const double eps1=1e-10,
       const double eps2=1e-13
    );

    void reduceMatEl(const char *F, int L,
       const wbMatrix<TQ> &Q1, const wbvector<TD> &dd,
       const wbMatrix<TQ> &Q2, const wbvector<TD> &dc,
       wbvector<unsigned> &I1, wbvector<unsigned> &I2,
       wbvector<TD> &dr
    ) const;

    void unRef() { 

        if (!isref && !qtype.isref && !itags.isref && !fdir.isref &&
            !QIDX.isref && !DATA.isref && !CGR.isref
        ) return;

        wblog(FL,"WRN %s() check this [%d%d%d%d%d%d%d]",FCT, isref,
        QIDX.isref, DATA.isref, CGR.isref, qtype.isref,itags.isref,fdir.isref);

        QIDX.unRef(); qtype.unRef(); itags.unRef(); fdir.unRef();

        if (DATA.isref) DATA.unRef(); else
        for (unsigned i=0; i<DATA.len; ++i) DATA[i]->unRef();

        if (CGR.isref) CGR.unRef();
    };

    QSpace& initIdentity(const QSpace &A, char dflag=0); 
    QSpace& initIdentity(const QType &t, const TQ *qs);
    QSpace& initIdentity_bare(const wbMatrix<TQ> &Q);

    QSpace& FusionTree(
       const char *F, int L, const QType &t, const wbMatrix<TQ> &Q,
       QSpace *A=nullptr);

    QSpace& setupCGR(const char *F=nullptr, int L=0); 

    QSpace& initIdentityCGS(const char *F=nullptr, int L=0);
    QSpace& initIdentityCGS(const char *F, int L, unsigned k, unsigned r=0);

    template<class TD2>
    QSpace& initIdentityCGS(const char *F, int L, const QSpace<TQ,TD2>& H);

    template<class TA>
    QSpace<TQ,TD>& initIdentityCG(
       wbvector< const QSpace<TQ,TA>* > A, 
       const wbindex &ia, char zflag=0 
    );

    template<class TA, class TB>
    QSpace<TQ,TD>& initIdentityCG( 
       wbvector< const QSpace<TQ,TA>* > A, const wbindex &ia,
       wbvector< const QSpace<TQ,TB>* > B, const wbindex &ib,
       char vflag='v'
    );

    template<class TA, class TB>
    QSpace<TQ,TD>& initIdentityCG(
       wbvector< const QSpace<TQ,TA>* > A, 
       wbvector< const QSpace<TQ,TB>* > B,
       const char *pstr=nullptr 
    ){ wbindex i;
       if (!cPVEC_nnz(FLF,A,'w') || !cPVEC_nnz(FLF,B,'w')) wblog(FL,
          "ERR %s() got empty input",FCT); 

       initIdentityCG(A,i,B,i); 
       if (pstr && pstr[0]) Permute(wbperm(pstr));
       return *this;
    };

    QSpace& initDiagonal(
        const wbMatrix<TQ> &Q,
        const wbvector<unsigned> &S,
        const wbvector<TD> &D,
        char rcflag
    );

    QSpace& init2DiffOp(
       const QSpace &A, unsigned ia, 
       const QSpace &B, unsigned ib
    );

    template<class TB>
    QSpace& initQT_safe(const QSpace<TQ,TB> &B) {
        clearQSpace(); 

        QIDX=B.QIDX; QDIM=B.QDIM; qtype=B.qtype; otype=B.otype;
        setupDATA();

        itags=B.itags; fdir=B.fdir; CGR=B.CGR;

        if (CGR.dim2) {
           if (qtype.len!=CGR.dim2) wblog(FL,
              "ERR %s() got CGR size mismatch (%dx%d <> %dx%d)",
              FCT,CGR.dim1,CGR.dim2,QIDX.dim1,qtype.len
           );
        }
        else if (QIDX.data) { 
           if (DATA.len && !qtype.allAbelian()) wblog(FL,
              "ERR %s() missing CGR data for %s (%dx%d <> %dx%d)",
              FCT,STR(qtype),CGR.dim1,CGR.dim2,QIDX.dim1,qtype.len
           );
        }
        return *this;
    };

    template <class TA>
    QSpace& initQ(const QSpace<TQ,TA> &A, char full=1){
       if (isref) wblog(FL,"WRN %s() got QSpace ref!",FCT);
       qtype=A.qtype; QDIM=A.QDIM; if (full) {
       otype=A.otype; itags=A.itags; fdir=A.fdir; }
       return *this;
    };

    template <class TA, class TB>
    QSpace& initQ(const char *F, int L,
       const QSpace<TQ,TA> &A, char full=1, const QSpace<TQ,TB> *B=nullptr
    ){
       if (B) A.checkQ(F,L,*B);
       return initQ(A,full);
    };

    void initQ(unsigned n, unsigned r, const QSpace &A){
       clearQSpace(); 
       QDIM=A.QDIM; qtype=A.qtype; QIDX.init(n,r*QDIM);
       setupDATA(FL);
       if (!qtype.allAbelian()) setupCGR(FL); 
    };

    void initQ(const wbMatrix<TQ> &Qk, const QSpace &A){ 
       clearQSpace(); 
       QDIM=A.QDIM; qtype=A.qtype; QIDX=Qk;
       setupDATA(FL);
       if (A.CGR.data) setupCGR(FL);
    };

    void init_itags( const char *F, int L, const char *s,
       unsigned r=-1, char lflag=-1);

    QSpace& SetTag(unsigned k, const itag_ &b) { 
       if ((--k)>=itags.len) wblog(FL,
          "ERR %s() index out of bounds (%d/%d)",FCT,k+1,itags.len);
       itags.data[k].tSet(b);
       return *this;
    };

    QSpace& SetTags(const itag_ &b) { 
       unsigned k=rank(FL);
       if (b.t && k>2) wblog(FL,
          "ERR %s() got '%s' for rank-%d QSpace !?",FCT,STR(b),k);
       for (k=0; k<itags.len; ++k) { itags.data[k].tSet(b); }

       return *this;
    };

    QSpace& SetConjTags(unsigned k1, unsigned k2) { 
       if (!k1 || !k2) wblog(FL,
          "ERR %s(%d) expects 1-based index",FCT,k1,k2);
       itags.SetConj(0,0,k1-1);
       itags.SetConj(0,0,k2-1); return *this;
    };

    QSpace& SetConjTag(unsigned k) { 
       if (!k) wblog(FL,"ERR %s(%d) expects 1-based index",FCT,k);
       itags.SetConj(0,0,k-1); return *this;
    };

    QSpace& SetConjTag(const char* F, int L, unsigned k) { 
       if (!k) wblog(F_L,"ERR %s(%d) expects 1-based index",FCT,k);
       itags.SetConj(F,L,k-1); return *this;
    };

    QSpace& SetFlag(const char *F, int L, unsigned k, unsigned l) {
       if (!k) wblog(F_L,"ERR %s(%d) expects 1-based index",FCT,k);
       itags.SetFlag(F_L,k-1,l);
       return *this;
    };
    QSpace& UnsetFlag(const char *F, int L, unsigned k, unsigned l) {
       if (!k) wblog(F_L,"ERR %s(%d) expects 1-based index",FCT,k);
       itags.UnsetFlag(F_L,k-1,l);
       return *this;
    };
    QSpace& UnsetFlags() { 
       itags.UnsetFlags();
       return *this;
    };
    QSpace& SetFlags(const char *F, int L, unsigned k, unsigned l=-1) {
       if (!k) wblog(F_L,"ERR %s(%d) expects 1-based index",FCT,k);
       itags.SetFlags(F_L,k-1,l);
       return *this;
    };

    QSpace& init_tags() { itags.init_tags(); return *this; };

    template<class TB>
    unsigned matchITags(
       const char *F, int L, const QSpace<TQ,TB> &B,
       ctrIdx &ia, ctrIdx &ib) const;

    QSpace& SetFDir(const QSpace &A) { fdir=A.fdir; return *this; };
    QSpace& SetFDir(const char *F, int L, const QSpace &A, const QSpace &B);

    QSpace& SetFDir_op( 
       const char *F, int L, const iFerm &fdr, const wbindex &I);

    QSpace& SetFDir_Atensor(const char *F, int L,
       const iFerm &fdrA, const wbindex &Ia, 
       const iFerm &fdrB, const wbindex &Ib);

    unsigned RemoveZLabels(
       const char *F, int L, wbMatrix<TQ> &QQ, wbMatrix<TQ> *Z=nullptr) const;

    void RemoveZLabels(const char *F, int L, wbMatrix<TQ> *Z){
       if (QDIM==0 && qtype.len) { QDIM=qtype.Qlenz(); }
       QDIM=RemoveZLabels(F_L,QIDX,Z);
    };

    void ExpandDiagonal(unsigned i1=1, unsigned i2=2); 

    size_t getDataSize(wbMatrix<size_t> &S) const; 
    size_t getDataSize(wbvector<size_t> &S) const;
    size_t getDataSize() const;
    size_t getCGSSize() const;
    wbstring totSize2Str() const;

    wbstring itags2Str() const {  
       return itags.toStr().data; 
    };

    unsigned map2Vec (const char *F, int L, wbvector<TD> &V,
       const wbvector<unsigned> &S=wbvector<unsigned>());

    unsigned map2Vec(TD *v, char Iflag=0);

    void ExpandQ( 
    CPAT<TQ,TD> &CP, double nrm, wbvector<unsigned> &xflag, char disp=0);

    void ResetRec(unsigned k) {
       if (k>=QIDX.dim1) wblog(FL,
          "ERR %s() index out of bounds (%d/%d)",FCT,k,QIDX.dim1);
       isConsistent(FL);

       QIDX.setRec(k,0); DATA[k]->init();
       for (unsigned j=0; j<CGR.dim2; ++j) { CGR(k,j).init(); }
    };

    void setDATA(const TD c) {
       for (unsigned i=0; i<DATA.len; ++i) { DATA[i]->set(c); }
    };

    void setRand(double nrm=0);

    void markUnequalZero(const TD m=1, const TD eps=0) {
        for (unsigned i=0; i<DATA.len; ++i)
        DATA[i]->markUnequalZero(m,eps);
    };

    void PrependSingletons(unsigned R);

    void swap(QSpace &B) {
       if (isref || QIDX.isref || DATA.isref || CGR.isref) wblog(FL,
          "ERR must not save referenced data (isref=%d,%d,%d)",
          isref, QIDX.isref,DATA.isref);
       if (B.isref || B.QIDX.isref || B.DATA.isref || B.CGR.isref) wblog(FL,
          "ERR must not save referenced data (isref=%d,%d,%d)",
          B.isref, B.QIDX.isref, B.DATA.isref);

        DATA .swap(B.DATA );
        QIDX .swap(B.QIDX );  SWAP(QDIM, B.QDIM );
        CGR  .swap(B.CGR  );  SWAP(otype,B.otype); 
        qtype.swap(B.qtype);
        itags.swap(B.itags);  SWAP(isref,B.isref); SWAP(mt,B.mt);
    };

    QSpace& save2(QSpace &B) {
        if (isref || QIDX.isref || DATA.isref || CGR.isref) wblog(FL,
           "ERR must not save referenced data (isref=%d,%d,%d)",
           isref, QIDX.isref,DATA.isref);
        for (unsigned i=0; i<DATA.len; ++i) {
           if (DATA[i]->isRef()) wblog(FL,"%NERR using save2 with "
              "referenced data\n(data[%d] has isref=%d)",
              i,DATA[i]->isRef()
           );
        }
        B.clearQSpace();    

        DATA .save2(B.DATA ); 
        QIDX .save2(B.QIDX ); B.QDIM=QDIM;
        CGR  .save2(B.CGR  );

        qtype.save2(B.qtype);
        itags.save2(B.itags);
        B.otype=otype; otype=QS_NONE; 

        B.isref=isref; isref=0;
        B.ctime=ctime;
        B.mt=mt;

        return B;
    };

    unsigned Reduce2AbelianOp(const char *F=nullptr, int L=0);

    char allAbelian() const {
       return qtype.allAbelian();
    };

    char isNonAbelian() const {
       return qtype.isNonAbelian();
    };

    bool permitsOM() const { 
       bool q=0;
       if (CGR || !qtype.allAbelian()) {
          q=qtype.permitsOM(rank(FL)); }
       return q;
    };

    bool hasOM() const {
       bool om=0; 

      #ifndef WB_SKIP_ASSERT
       int q=QIDX.isSorted(), pm=permitsOM();
       if (! q   ) wblog(FL,"WRN %s() got non-sorted QIDX",FCT);
       if (!(q&2)) wblog(FL,"WRN %s() got non-unique QIDX "
          "(%s @ r=%d; %d, %d)",FCT,STR(qtype),rank(),q,pm);
      #endif

       unsigned i=0, n=CGR.numel();
       for (; i<n; ++i) { if (CGR[i].wnumel()>1) { om=1; break; }}
       return om;
    };

    bool allAdditive() const { 
       return qtype.allAdditive();
    };

    bool isEmpty() const;

    int isFerm(const char *F=nullptr, int L=0, const iFerm *f=nullptr) const {
       int q=fdir.check(itags.len, f);
       if (q<0 && F) wblog(F,L,
          "ERR %s() invalid data (%s /%d)", C_FCT,STR(fdir),itags.len);
       return q; 
    };

    explicit operator bool() const { return !isEmpty(); } 
    bool operator!() const { return isEmpty(); } 

    bool isConsistent(  
       const char *F, int L, unsigned &r, char level=3) const; 

    bool isConsistent(
      const char *F=nullptr, int L=0 
    ) const { unsigned r=-1; return isConsistent(F,L,r,3); };

    bool isConsistent(unsigned &r, char level=3) const {
        return isConsistent(nullptr,0,r,level); };

    bool isConsistent_r(unsigned r, char level=3) const {
        return isConsistent(nullptr,0,r,level);
    };

    bool isScalar() const { 
       unsigned r=rank(FL);
       if (r || DATA.len!=1) { return 0; }

       if (QIDX.dim1!=1 || DATA[0]->numel()!=1) { 
          wblog(FL,"ERR %s() inconsistency (%d; S=%s)",
          FCT,QIDX.dim1,SSTR_(DATA[0]));
       }
       return 1;
    };

    bool isComplex() const {
       unsigned i=0, n=DATA.len; isConsistent(FL);
       for (; i<n; ++i) { if (DATA[i]->isComplex()) return 1; }
       return 0;
    };

    bool isDiagBlock(unsigned k) const;

    bool isBlockDiagMatrix(const char *F, int L, char dflag=0) const;
    bool isBlockDiagMatrix(char dflag=0) const {
         return isBlockDiagMatrix(nullptr,0,dflag); }

    bool isDiagMatrix(const double eps=1e-14) const;
    bool isIdentityMatrix(const TD eps=1e-14) const;
    char hasIdentityCGS(const TD eps=1e-14) const;

    bool isR3Op() const {
       return (rank(FL)==3 && otype==QS_OPERATOR);
    };

    bool isOp(int r=-1, char lflag=0) const { 
        return itags.isOp(r,lflag); };

    int isOperator(unsigned *r_, char xflag=0) const;

    bool isOperator() const { return (isOperator(nullptr,0)>0 ? 1 : 0); };

    bool isAtensor(unsigned r=-1) const { return itags.isAtensor(r); }

    int  checkAbelianOp() const; 
    bool isAbelianOp() const { return (checkAbelianOp()>0); };

    bool isQSym(const char *F=0, int L=0, char dflag=0) const;

    bool isHConj(
      const char *F=0, int L=0, double eps=1e-14, char vflag=0
    ) const;

    bool isAHerm(
      const char *F=0, int L=0, double eps=1e-14, char vflag=0
    ) const;

    int isSingleton(unsigned k) const;

    template <class TB>
    bool hasSameQ(const QSpace<TQ,TB> &B) const;

    bool hasSameQ( 
       unsigned r1, unsigned k1, 
       unsigned r2, unsigned k2  
    ) const {
       if (memcmp( Qref(r1,k1), Qref(r2,k2), QDIM*sizeof(TQ) )) return 0;
       if (DATA[r1]->dim0(k1) != DATA[r2]->dim0(k2)) {
          wblog(FL, 
             "WRN %s()=1, yet DATA size mismatch (%s <> %s; %d,%d)",
              FCT, SSTR_(DATA[r1]), SSTR_(DATA[r2]), k1+1, k2+1);
          return 0;
       }
       return 1;
    };

    template <class T2>
    bool hasSameRank(const QSpace<TQ,T2> &B) const {
        return (QIDX.dim2==B.QIDX.dim2 && QDIM==B.QDIM);
    };

    char checkQ(const char *F=nullptr, int L=0) const;

    template <class TB>
    void checkQ(
       const char *F, int L, const QSpace<TQ,TB> &B,
       char aflag=0
     ) const;

    char checkQ_CGR(const char *F, int L, cUVEC *I=nullptr) const;

    int check_CGR_cgw(const char *F=0, int L=0) const;

    void checkScalarCGS(const char *F=0, int L=0) const;

    template <class TB>
    bool sameType(
       const QSpace<TQ,TB> &B, int r=-1, const char *istr=nullptr
    ) const; 

    QDir& getQDir(QDir &qdir) const;
    char getQDir(unsigned i) const; 

    bool sameQDir(const char *F, int L, const QSpace &B) const;

    bool hasQOverlap(const QSpace &A) const;

    template <class TB>
    int hasQOverlap(
       unsigned ica, const QSpace<TQ,TB> &B, unsigned icb,
       const char type=0
    ) const;

    int getQOverlapU(
       const char *F, int L, const QSpace &B,
       wbMatrix<TQ> &QQ, wbMatrix<int> &IQ,
       char olonly=0
    ) const;

    int getQOverlap(
       const char *F, int L, const QSpace &B,
       wbMatrix<TQ> &QQ, wbMatrix<int> &IQ,
       wbvector<widx_t> &DQ 
    ) const;

    bool findDimQ(const TQ* q, unsigned &d) const;
    bool findDimQ(const TQ* q, unsigned k, unsigned &d) const;

    bool operator==(const QSpace &B) const {
       if (this!=&B) {
          if (!(QIDX ==B.QIDX )) { return 0; }
          if (!(qtype==B.qtype)) { return 0; }
          if (!(CGR  ==B.CGR)  ) { return 0; }
          if (!DATA.deepEqualP(B.DATA)) { return 0; }
       }
       return 1;
    };

    bool operator!=(const QSpace &B) const { return (!((*this)==B)); }

    template<class T>
    void operator*= (const T fac);

    template<class T>
    QSpace& times(const T fac, QSpace &B) const {
       B=*this; B*=fac; return B;
    };

    QSpace& operator= (const QSpace &A) { 
       if (this!=&A) { init(A); }
       return *this;
    };

    void operator+= (const QSpace &Q);
    void operator-= (const QSpace &Q);

    QSpace& Cat(const char *F, int L, 
       const QSpace<TQ,TD> &A, const QSpace<TQ,TD> &B,
       TD afac=1, TD bfac=1, char uflag=1
    );

    QSpace& plus( 
       const QSpace &B, QSpace &C, TD bfac=1, char vflag=0) const;

    QSpace& minus(
       const QSpace &B, QSpace &C, TD bfac=1, char vflag=0) const;

    QSpace& Plus (const QSpace &A, TD bfac=1);
    QSpace& Minus(const QSpace &A, TD bfac=1);

    QSpace& plus_plain(const QSpace &B, QSpace &C, TD bfac=1) const;

    int NormCGW(char full=0, char skipzeros=1, char rcpy=0);
    int isNormCGW(char full=0) const;

    TD norm2(char checks=1) const;
    TD norm (char checks=1) const { return sqrt(norm2(checks)); }

    TD normDiff2(const QSpace &B) const;
    TD scalarProd(const QSpace &B) const;

    void checkNorm(
       const char *F, int L, double nrm, char vflag=1, double eps=1e-12
     ) const;

    QSpace& TimesEl(const QSpace &B, char conj=0, unsigned r2=-1);

    TD sumData() const;

    template <class TM>
    wbMatrix<TM>& diag(wbMatrix<TM> &xd) const;

    QSpace& TensorProd(const QSpace &B, QSpace &C,
       const char aflag='N', const char bflag='N') const;

    template <class TC, ENABLE_IF_isComplex(TC) >
    QSpace<TQ,TC>& TensorProd(const QSpace &B, QSpace<TQ,TC> &C,
       const char aflag='N', const char bflag='N'
     ) const { QSpace<TQ,TD> C_;
       this->TensorProd(B,C_,aflag,bflag);
       return C.initT(C_);
    };

    template <class TC, ENABLE_IF_isComplex(TC) >
    QSpace<TQ,TC>& TensorProd(const QSpace<TQ,TC> &B, QSpace<TQ,TC> &C,
       const char aflag='N', const char bflag='N'
     ) const { QSpace<TQ,TC> A; A.initT(*this);
       return A.TensorProd(B,C,aflag,bflag);
    };

    void TensorProdUnity(const wbMatrix<TQ> &Q, QSpace &C) const;

    void EigenSymmetric(wbMatrix<double> &Etot) const {
       QSpace Ak,Ad;
       QSpace<TQ,double> Ek,Et;
       wbMatrix<unsigned> DB; int Nkeep=-1; 

       EigenSymmetric(Ak,Ad,Ek,Et,Etot,DB,Nkeep);
    };

    void EigenSymmetric(
       QSpace &Ak, QSpace<TQ,double> &Ek, wbMatrix<double> &Etot
    ) const {
       QSpace Ad;
       QSpace<TQ,double> Et;
       wbMatrix<unsigned> DB; int Nkeep=-1; 

       EigenSymmetric(Ak,Ad,Ek,Et,Etot,DB,Nkeep);
    };

    void EigenSymmetric(QSpace &Ak, QSpace<TQ,double> &Ek) const {
       QSpace Ad;
       QSpace<TQ,double> Et;
       wbMatrix<unsigned> DB; int Nkeep=-1; 
       wbMatrix<double> Etot;

       EigenSymmetric(Ak,Ad,Ek,Et,Etot,DB,Nkeep);
    };

    void EigenSymmetric( 
       QSpace &Ak, QSpace &Ad,
       QSpace<TQ,double> &Ek, QSpace<TQ,double> &Et, wbMatrix<double> &Etot,
       wbMatrix<unsigned> &DD, int &Nkeep,
       double Etrunc=0, double *E0=nullptr,
       char mKD=0, 
       const wbperm &P=wbperm(), double eps=0., double b=0., int dmax=-1,
       const char *sort=nullptr
    ) const;

    void GroupIndizes(unsigned K);

    int ExpandOM(
       const char *F, int L, wbvector<unsigned> &iOM, char fflag=0);

    int ExpandOM_( 
       const char *F, int L, unsigned ia, 
       QSpace<TQ,TD> &B, unsigned ib);

    int ExpandOM_(const char *F, int L, unsigned ia, unsigned ib) {
       return ExpandOM_(F,L,ia,*this,ib);
    };

    void Append(unsigned k); 
    void Append(const wbvector<TQ> &, const wbarray<TD> &D);

    QSpace& Append(const char *F, int L, const QSpace &A, char uflag=1);

    int Append2AndDestroy(
       const char *F, int L, QSpace &A, char unique=0);

    template <class TA, ENABLE_IF_isComplex(TA) >
    int Append2AndDestroy(
       const char *F, int L, QSpace<TQ,TA> &A, char unique=0
    ){
        QSpace<TQ,TA> X; X.initT(*this);
        return X.Append2AndDestroy(F,L,A,unique);
    };

    void TransferSpace(const char *F, int L,
       QSpace &B, unsigned k,
       const QSpace<TQ,double> &W, double eps
    );

    void TransferSpace(const char *F, int L,
       QSpace &B, unsigned k
    );

    unsigned MakeUnique(); 

    bool isZero(const double eps=0., char flag=0) const;
    unsigned SkipZeroData( 
       double eps=0., char bflag=0, char cgflag=1,
       char all=0 
    );

    unsigned skipZeroOffDiag(double eps=0., char bflag=0);

    unsigned SkipEmptyData(const char *F=0, int L=0);

    bool gotZeroData(unsigned i, const double eps=0.) const;

    void Set(unsigned k, const wbvector<TQ> &Q, const wbarray<TD> &D);
    void SetQ(unsigned i, unsigned k, const TQ *d);

    QSpace& Select(const wbindex &I, char data_only=0);

    template <class TM>
    QSpace& SelectMarked(const wbvector<TM> &mark){
       if (mark.len!=QIDX.dim1) wblog(FL,
          "ERR %s() size inconsistency (%d/%d)",FCT,mark.le,QIDX.dim1);
       wbindex I; mark.find(I);
       return Select(I);
    };

    void Sort(const QSpace &A); 
    void Sort() {
       wbperm P; QIDX.SortRecs(P); DATA.Select(P);
       if (qtype.allAbelian()) { 
          if (!CGR.isEmpty()) wblog(FL,"ERR %s() "
         "invalid all-abelian CGR (%dx%d)",FCT,CGR.dim1,CGR.dim2);
       }
       else CGR.recPermute(P);
    };

    unsigned getOM(
       wbMatrix<unsigned> &OM, wbvector<unsigned> *om1=nullptr,
       char wflag=0) const; 

    unsigned cgsDimScalar(unsigned k) const;
    unsigned cgsDim(unsigned k, unsigned r) const;
    unsigned getDim(unsigned k, widx_t *D=nullptr) const; 

    unsigned getDIM(unsigned k) const { 
       widx_t D=0; getDim(k,&D);
       return D;
    };

    void getDim(wbvector<widx_t> &D, wbvector<widx_t> *DD=nullptr) const;

    wbvector<widx_t> getDim() const {
       wbvector<widx_t> D; getDim(D); return D;
    }

    void getQDim(unsigned k, 
       wbMatrix<TQ> &Q,
       wbvector<widx_t>&S,
       wbMatrix<widx_t>*SC=nullptr  
    ) const;

    void getQDim(wbMatrix<TQ> &Q,
       wbvector<widx_t> &S,
       wbMatrix<widx_t>*SC=nullptr  
    ) const;

    void getDRange(TD &dmin, TD &dmax) const; 

    void getDistQtot(wbMatrix<TQ> &Qtot, wbvector<double> &w) const;
    widx_t getDQtot() { 
       wbMatrix<TQ> Q; QIDX.blockSum(QDIM,Q); Q.makeUnique();
       return Q.dim1;
    };

    wbMatrix<TQ>& getQsub(const wbindex &I, wbMatrix<TQ> &QI) const;
    wbMatrix<TQ>& getQsub(
    const unsigned k, wbMatrix<TQ> &Qk, char iflag=0) const;

    void getQsub(const wbindex &I, wbMatrix<TQ> &QI, wbMatrix<TQ> &Q2) const;

    wbMatrix<TQ> getQsub(const unsigned k) const {
       wbMatrix<TQ> Qk; getQsub(k,Qk); return Qk;
    };

    wbMatrix<TQ>& getQtot(wbMatrix<TQ> &Qtot, char uflag=0) const;
    wbMatrix<TQ>& getQsum(const wbindex &J, wbMatrix<TQ> &Q) const;
    wbMatrix<TQ>& getQsum(const wbMatrix<TQ> &Q0, wbMatrix<TQ> &Q) const;

    wbMatrix<TQ>& getQfinal(const wbindex &J, wbMatrix<TQ> &Q) const;

    QSpace& getSub(const wbindex &I, QSpace &A, char ref=0) const;

    QSpace& Set2Sub(const wbindex &I, char ref=0) { 
       QSpace X; getSub(I,X,ref);
       return X.save2(*this);
    };

    QSpace& getSubInit(const wbindex &I, QSpace &A) const;

    void saveSub2(QSpace &A, const wbindex &I);

    void toFull(
       wbarray<TD> &A,                    
       wbvector< wbMatrix<TQ> > &Q,       
       wbvector< wbvector<unsigned> > &S  
    ) const;

    void toFull2(
       wbarray<TD> &A,        
       wbMatrix<TQ> &Q,       
       wbMatrix<widx_t> &S   
    ) const;

    void ExtendQ(const QSpace &B);

    void takeMainQS(char disp=1);

    void setQ(unsigned k, const QSpace &B, unsigned k0);
    void setQRec(unsigned r, unsigned j, const TQ *q0, const TQ *q2=nullptr);

    bool hasEqualQ(unsigned i1, unsigned i2) const;

    double maxDiff(const QSpace &B) const;

    int getFermZ2(wbMatrix<unsigned> &Z2, wbvector<unsigned> *z2=nullptr) const;
    int getFermZ2(wbperm &P, MVEC &z2,
       const ctrIdx *ic=nullptr, const wbindex *zx=nullptr) const;

    int getFermSigns(wbperm &P, MVEC &fS, 
       const ctrIdx *ic=nullptr, const wbindex *zx=nullptr) const {

       int q=getFermZ2(P,fS,ic,zx);
           Wb::z2ferm2sign(fS);
       return q;
    };

    int getFermZ2_conj(MVEC &z2) const;

    void PermuteQ(const wbperm &P);
    void permuteQ(wbMatrix<TQ> &Q, const wbperm &P) const;

    QSpace& Permute(wbperm P, char rcpy=0);         
    QSpace& permute(QSpace &B, wbperm P) const;     

    void permute(
       const char *pstr, QSpace &B, unsigned offset=1   
     ) const { permute(B, wbperm(pstr,offset) ); };

    void Permute(const char* pstr, unsigned offset=1    
     ) { Permute(wbperm(pstr,offset)); };

    QSpace& PermuteFirstTo(unsigned k);
    QSpace& permuteFirstTo(unsigned k, QSpace &C) const;

    QSpace& PermuteLastTo (unsigned k);
    QSpace& permuteLastTo (unsigned k, QSpace &C) const;

    char permute_to(char pos, 
       ctrIdx &ic, QSpace<TQ,TD> &B, char ref=0) const;

    void FlipQ();

    QSpace& trace(
       const char *F, int L, ctrIdx &i1, ctrIdx &i2, QSpace &C) const;

    TD trace() const;

    QSpace& revertLeg( 
       const char *F, int L, unsigned ic, QSpace &C) const;

    QSpace& RevertLeg(const char *F, int L, unsigned ic) {
       QSpace C; revertLeg(F,L,ic,C);
       return C.save2(*this);
    };

    QSpace& contractMat(
       unsigned ia, const QSpace &B, unsigned ib, QSpace &C) const;
    QSpace& contractMat(
       const ctrIdx &ia, const QSpace<TQ,TD> &B, const ctrIdx &ib,
       QSpace<TQ,TD> &C) const;

    template <class TB, class TC>
    int contract_getIdxSet(const char *F, int L,
       const ctrIdx &ica, const QSpace<TQ,TB> &B, const ctrIdx &icb,
       wbindex &Ia, wbindex &Ib, wbvector<widx_t> &Dc,
       QSpace<TQ,TC> &C, MVEC &fA, MVEC &fB
    ) const;

    template<class TB, class TC>
    double contract( 
       const char *F, int L,   ctrIdx ica, 
       const QSpace<TQ,TB> &B, ctrIdx icb, QSpace<TQ,TC> &C,
       const wbperm &P=wbperm(), char preview=0 
     ) const;

    template<class TB, class TC>
    double contract(const char *F, int L,
       unsigned ia, const QSpace<TQ,TB> &B, unsigned ib, 
       QSpace<TQ,TC> &C, const wbperm &P =wbperm(), char preview=0
     ) const {

       unsigned ra=rank(FL), rb=B.rank(FL);
       if (!ia || (ia>ra && (QIDX.dim1 || QIDX.dim2))) wblog(FL,
          "ERR %s() index out of bounds (ia=%d/%d; %s)",
          FCT,ia,ra,SSTR(QIDX));
       if (!ib || (ib>rb && (B.QIDX.dim1 || B.QIDX.dim2))) wblog(FL,
          "ERR %s() index out of bounds (ib=%d/%d; %s)",
          FCT,ib,rb,SSTR(B.QIDX)
       );

       ctrIdx ica(1,&(--ia)), icb(1,&(--ib));
       return contract(F,L,ica,B,icb,C,P,preview);
    };

    template<class TB, class TC>
    int contract(const char *F, int L,
       unsigned ia,  
       const QSpace<TQ,TB> &B, const char *sb, 
       QSpace<TQ,TC> &C, const wbperm &P =wbperm(), char preview=0
    ) const {

       ctrIdx ica(1,&(--ia)), icb;

       if (isEmpty() || B.isEmpty()) { C.init(); return 1; }
       if (ia>=itags.len) {
          wblog(FL,"ERR %s() index out of bounds (ia=%d/%d; %s)",
          FCT,ia+1,itags.len,STR(itags));
       }

       if (sb && isdigit(sb[0]))
            { icb.init(F_L,sb); }
       else { icb.initOp(FL,itags[ia],B,sb); } 

       return contract(F,L,ica,B,icb,C,P,preview);
    };

    template<class TB, class TC>
    int contract(const char *F, int L,
       const char *sa, 
       const QSpace<TQ,TB> &B, unsigned ib, 
       QSpace<TQ,TC> &C, const wbperm &P =wbperm(), char preview=0
     ) const {

       ctrIdx ica, icb(1,&(--ib));

       if (isEmpty() || B.isEmpty()) { C.init(); return 1; }

       if (ib>=B.itags.len) wblog(FL,
          "ERR %s() index out of bounds (ib=%d/%d; %s)",
          FCT,ib+1,B.itags.len,STR(B.itags));

       if (sa && isdigit(sa[0]))
            { ica.init(F_L,sa); }
       else { ica.initOp(FL,B.itags[ib],*this,sa); } 

       return contract(F,L,ica,B,icb,C,P,preview);
    };

    template<class TB, class TC> inline
    int contract(const char *F, int L,
       const char *sa, const QSpace<TQ,TB> &B, const char *sb,
       QSpace<TQ,TC> &C, const wbperm &P =wbperm(), char preview=0
    ) const {

       ctrIdx ica(F_L,sa), icb(F_L,sb); 
       return contract(F,L,ica,B,icb,C,P,preview);
    };

    template<class TB, class TC>
    int contract(
       const char *sa, const QSpace<TQ,TB> &B, const char *sb,
       QSpace<TQ,TC> &C, const wbperm &P =wbperm(), char preview=0
    ) const { return contract(FL,sa,B,sb,C,P,preview); };

    template<class TB, class TC>
    int contract(
       const ctrIdx &ica, const QSpace<TQ,TB> &B, const ctrIdx &icb,
       QSpace<TQ,TC> &C, const wbperm &P =wbperm(), char preview=0
    ) const { return contract(FL,ica,B,icb,C,P,preview); };

    template<class TB, class TC>
    int contract(
       unsigned ia, const QSpace<TQ,TB> &B, unsigned ib, 
       QSpace<TQ,TC> &C, const wbperm &P =wbperm(), char preview=0
    ) const { return contract(FL,ia,B,ib,C,P,preview); };

    mxArray* toMx(bool setcls=0) const;      

    mxArray* save2Mx(char vflag=0); 

    mxArray* QIDX_toMx() const; 
    mxArray* INFO_toMx() const;
    mxArray* DATA_toMx() const;
    mxArray* DATA_save2Mx(char vflag=0);

    mxArray* mxCreateStruct(unsigned m, unsigned n) const;
    void add2MxStruct(mxArray *S, unsigned i, bool setcls=0, char chk=0) const;

    void save2MxStruct(mxArray *S, unsigned i, char tst=0, char vflag=0);

    mxArray* mxCreateCell(
       unsigned m QS_UNUSED_VAR,
       unsigned n QS_UNUSED_VAR) const { wblog(FL,"ERR %s()",FCT); return 0; };
    void add2MxCell(
       mxArray *S QS_UNUSED_VAR,
       unsigned i QS_UNUSED_VAR,
       char tst   QS_UNUSED_VAR = 0) const { wblog(FL,"ERR %s()",FCT); };

    void put (
       const char *F, int L,
       const char *vname, const char *ws="caller"
    ) const {
       wblog(F,L,"I/O putting '%s' to %s", vname, ws);
       put(vname,ws);
    };

    void put (const char *vname, const char *ws="caller") const {
       mxArray *S=toMx(); 
       mxPutAndDestroy(0,0,S,vname,ws);
    };

    wbstring sizeStr(char vflag=0) const;

    wbstring sizeStrQ(const char *F=nullptr, int L=0) const;

    wbstring sizeStrD(unsigned i) const { 
       if (i>=DATA.len) wblog(FL,
          "ERR %s() index out of bounds (%d/%d)",FCT,i,DATA.len);
       return DATA[i]->sizeStr(rank(FL),"_");
    };

    void print_SIZE(const char *F=nullptr, int L=0) const;
    void print_rankDATA(const char *F=nullptr, int L=0) const;

    wbstring qStr() const { return qtype.toStr(); };

    wbstring qrec2Str(unsigned k) const {
        return QIDX.rec2Str(k,""," ",  QDIM," ;");
    };

    void disp_cgs(const char *F, int L, const char *istr) const;

    void info(const char *vname=nullptr,
       char nl=0,    
       char nind=2   
    ) const;

    void print() const { print("",3); }; 
    void print(const char *vname, char vflag=3) const;

    void print_tst_( 
       const char *F, int L, const char *istr, unsigned k=0) const {

       if (k>=DATA.len) { wblog(FL,
          "WRN %s() index out of bounds (%d/%d)",FCT,k,DATA.len);
          return;
       }
       unsigned i=0;
       const TD *dk=DATA[k]->data;
       wbvec<char> s(64); s.catf(FL,"%8s",NSTR(dk[0]));
       if (CGR.dim2) {
          for (; i<CGR.dim2; ++i) {
          s.catf(0,0," %c %4s",i?'*':'x',RATS(CGR(k,i).wel(0))); }
       }
       wblog(F_L,"TST %s @%2d/%d: %-32s %s",istr,k+1,DATA.len,s.data,STR(itags));
    };                         

    void initOType(const char *F, int L, const mxArray* a);

    wbstring otype2Str(const char *fmt="%s") const;

    int ConjOpScalar() { 
       int e=0; 
       unsigned r=rank(FL);
       if (r!=3) {
          if (r!=2) { wblog(FL,"WRN %s() got rank-%d QSpace !?",FCT,r); }
          return e;
       }
       for (unsigned n=CGR.numel(), i=0; i<n; ++i) {
           if ((e=CGR[i].HConjOpScalar())) { return e; }
       }

       if (itags.len!=3) wblog(FL,
          "ERR %s() got itags.len=%d/3 !?",FCT,itags.len);
       itags[2].Conj();

       return e;
    };

    QSpace& SortDegQ(const char *F=0, int L=0) {
       if (CGR) { unsigned i=0, n=CGR.numel();
          for (; i<n; ++i) { CGR.data[i].SortDegQ(F_L); }
       }
       return *this;
    };

    QSpace& Conj(char check_ferm=1) {
       for (unsigned i=0; i<itags.len; ++i) { itags.data[i].Conj(); }
       fdir.Conj();

       DConj(check_ferm);
       RConj();

       return *this;
    };

    QSpace& DConj(char check_ferm=1) { 
       unsigned i=0, n=DATA.len;
       if (check_ferm && Wb::envFERM) {
          MVEC z2; getFermZ2_conj(z2);
          if (z2) {
             if (z2.len!=n) wblog(FL,"ERR %s() %d/%d",FCT,z2.len,n);
             if (WbUtil<TD>::isComplex())
                  { for (; i<n; ++i) { DATA[i]->ConjTimes(z2[i]%2 ? -1:1); }}
             else { for (; i<n; ++i) { DATA[i]->    Times(z2[i]%2 ? -1:1); }}
             return *this;
          }
       }

       if (WbUtil<TD>::isComplex()) {
          for (; i<n; ++i) { DATA[i]->Conj(); }
       }
       return *this;
    };

    QSpace& RConj() { 
       if (CGR) {
          unsigned i=0, n=CGR.numel();
          for (; i<n; ++i) { CGR[i].Conj(); } 
          SortDegQ(FL); 
       }
       return *this;
    };

    QSpace& HConj() {
       QSpace B; this->save2(B);
       return B.hconj(*this);
    };

    QSpace& hconj(QSpace &B) const {
       if (isEmpty()) { B.clearQSpace(); }
       else { unsigned r=rank(FL);
          if (r%2 && (r!=3 || otype!=QS_OPERATOR)) wblog(FL,"ERR %s() "
             "requires even-rank object or operator (%dx%d/%d; %s)",
              FCT, QIDX.dim1, QIDX.dim2, QDIM, QS_STR[otype]);
          transp(B); B.Conj(); B.otype=otype;
       }
       return B;
    };

    QSpace& transp(QSpace &B) const {
       if (isEmpty()) { B.clearQSpace(); }
       else {
          unsigned r=rank(FL); wbperm Pt;
             if ((r%2)==0) { Pt.initTranspose(r); } else
             if (r==3 && otype==QS_OPERATOR) { Pt.init("213"); }
             else wblog(FL,"ERR %s got rank-%d object",FCT,r);
          permute(B,Pt);
       }
       return B;
    };

    const QSpace& opA(char tflag, QSpace& X) const {
       if (!strchr("NTC",tflag))
       wblog(FL,"ERR %s() invalid flag %c<%d>",FCT,tflag,tflag);

       if (tflag=='N' || isEmpty()) return (*this);
       else {
          wbperm P2("21"); unsigned r=rank(FL);

          if (r!=2) wblog(FL,
          "ERR %s() only applicable to rank-2 tensors (%d)",FCT,r);

          permute(X,P2); if (tflag=='C') X.Conj();
          return X;
       }
    };

    wbMatrix<TQ> QIDX;

    wbvector<wbarray<TD>*> DATA;
    wbMatrix< CRef<TQ> > CGR;

    QVec qtype; 

    unsigned QDIM;

    QS_TYPES otype; 
    iTags itags;

    iFerm fdir;

    Wb::MTYPE mt;

    char isref;

    double ctime;

  protected: 
  private:

    bool isSym_aux(
      const char *F, int L, const char *fct,
      RTD eps=1e-14, char symflag='s', char vflag=0
    ) const;

    void recSave2(unsigned k, unsigned i);  

}; 

template <class TQ, class TD>
void getQDimGen(
   wbvector< const QSpace<TQ,TD>* > A,
   wbMatrix<TQ> &Q, wbvector<widx_t> &S,
   wbindex I, 
   wbMatrix<widx_t> *SC_=nullptr 
);

template <class TD>
class IterOM_DATA { 

  public:

    template <class TQ>
    IterOM_DATA(const QSpace<TQ,TD> &A, unsigned i)
     : iter(0), niter(1), nS(0), data0(nullptr) {

       unsigned r, rA=A.rank(FL);
       if (i>=A.DATA.len) wblog(FL,
          "ERR %s() index out of bounds (%d/%d)",FCT,i,A.DATA.len);
       data0=A.DATA[i]->data;
       X.SIZE.init2ref(A.DATA[i]->SIZE); r=X.SIZE.len;
       if (!rA || r<rA || r>rA+1) wblog(FL,
          "ERR %s() rank inconsistency (r=%d/%d)",FCT,r,rA);

       if (r>rA) { niter=X.SIZE[rA]; --X.SIZE.len; }
       nS=X.SIZE.prod();
    };

    void begin() { iter=0; }
    bool end() { return iter>=niter; }

    IterOM_DATA& operator++() { ++iter; return *this; };

    wbarray<TD> DATA() {
       if (iter>=niter) wblog(FL,
          "ERR %s() index out of bounds (%d/%d)",FCT,iter,niter);
       return X.init2ref(data0+iter*nS);
    };

   unsigned iter, niter;
   size_t nS;
   TD *data0;
   wbarray<TD> X;

  protected: 
  private:
};

template <class TQ, class TD>
template <class T> inline
void QSpace<TQ,TD>::operator*= (const T fac) {
   if (fac==T(1)) return;
   for (unsigned i=0; i<DATA.len; ++i)
   (*DATA[i]) *= fac;
};

template <class TQ, class TD> inline
void QSpace<TQ,TD>::operator+= (const QSpace &B) {
   QSpace<TQ,TD> X;
   this->plus(B,X); X.save2(*this);
};

template <class TQ, class TD> inline
void QSpace<TQ,TD>::operator-=(const QSpace &B) {
   QSpace<TQ,TD> X;
   this->plus(B,X,-1); X.save2(*this);
};

template <class TQ, class TD> inline
QSpace<TQ,TD>& QSpace<TQ,TD>::Plus(const QSpace &B, TD bfac) {
   QSpace<TQ,TD> X;
   this->plus(B,X,bfac); return X.save2(*this);
};

template <class TQ, class TD> inline
QSpace<TQ,TD>& QSpace<TQ,TD>::Minus(const QSpace &B, TD bfac) {
   QSpace<TQ,TD> X;
   this->plus(B,X,-bfac); return X.save2(*this);
};

template <class TQ, class TD> inline
QSpace<TQ,TD>& QSpace<TQ,TD>::minus(
   const QSpace &B, QSpace &C, TD bfac, char vflag
) const {
   plus(B, C,-bfac,vflag); return C;
}

template <class TQ, class TD>
QSpace<TQ,TD>& QSpace<TQ,TD>::plus( 
   const QSpace &B, QSpace &C, TD bfac, char vflag
 ) const {

   if (B.isEmpty() || !bfac) { C=(*this); } else
   if (this->isEmpty()) { C=B; if (bfac!=1.) { C*=bfac; }} else
   if (gotCGS(FL)<2) { plus_plain(B,C,bfac); }
   else {
      C=B; if (bfac!=1.) { C*=bfac; }
      C.Append(FL,*this, vflag<1 ? 'q':vflag);
   }

   return C;
};

template <class TQ, class TD>
TD QSpace<TQ,TD>::scalarProd(const QSpace<TQ,TD> &B) const {

   unsigned i; int i1,i2;
   wbMatrix<int> IQ;
   wbMatrix<TQ> QQ;
   TD x=0;

   if (this==&B) return norm2();

#ifdef WB_CLOCK
   Wb::Clock clk("QS:scalarProd",1); 
#endif

   if (isEmpty() || B.isEmpty()) {
      wblog(FL,"WRN Overlap with empty object (%d,%d) !?",
      isEmpty(), B.isEmpty()); return x;
   }

   getQOverlapU(FL,B,QQ,IQ,1); 

   for (i=0; i<IQ.dim1; ++i) {
       i1=IQ(i,0); i2=IQ(i,1);
       if (i1<0 || i2<0) wblog(FL, "ERR i1=%d, i2=%d !?", i1, i2);

       x+=(DATA[i1]->scalarProd(*B.DATA[i2]));
   }

   return x;
};

template <class TQ, class TD>
QSpace<TQ,TD>& QSpace<TQ,TD>::TensorProd(
   const QSpace<TQ,TD> &B, QSpace<TQ,TD> &C,
   const char aflag, const char bflag
 ) const {

   if (isEmpty() || B.isEmpty()) { C.clearQSpace(); return C; }

   if (gotCGS(FL)>0 || B.gotCGS(FL)>0) wblog(FL,
      "ERR %s() not implemented for non-abelian symmetries (got %s)",
      FCT,STR(qtype));

   if (aflag!='N' || bflag!='N') {
      QSpace<TQ,TD> AX, BX;
      const QSpace<TQ,TD> *Ap=this, *Bp=&B;

      if (!strchr("NTC",aflag) || rank(FL)!=2) wblog(FL,
         "ERR invalid aflag=%c<%d> (rank-%d)",aflag,aflag,rank(FL));
      if (!strchr("NTC",bflag) || B.rank(FL)!=2) wblog(FL,
         "ERR invalid bflag=%c<%d> (rank-%d)",bflag,bflag,B.rank(FL));

      if (aflag!='N') {
          permute("21",AX); if (aflag=='C') AX.Conj();
          Ap=&AX;
      }
      if (bflag!='N') {
         B.permute("21",BX); if (bflag=='C') BX.Conj();
         Bp=&BX;
      }

      return Ap->TensorProd(*Bp,C,'N','N');
   }

   unsigned i,ia,ib,l,r;
   const unsigned n1=QIDX.dim1, n2=B.QIDX.dim1, R=rank(FL), s=QDIM*sizeof(TQ);
   const TQ *qa, *qb; TQ *q;

   if (R!=B.rank(FL) || QDIM!=B.QDIM) wblog(FL,
      "ERR %s() size mismatch (%d,%d)",FCT,R,B.rank(FL));

   C.clearQSpace();
   C.QIDX.init(n1*n2, 2*QIDX.dim2); C.QDIM=QDIM; C.setupDATA();

   C.initQ(FL,*this,'f',&B);

   if (itags.len || B.itags.len) { C.itags.init(2*R);
      if (!B.itags.len) {
         if (itags.len!=R) wblog(FL,
            "ERR %s() invalid itag %s",FCT, STR(itags));
         for (i=0; i<R; ++i) {
            C.itags[2*i]=C.itags[2*i+1]=itags[i];
         }
      }
      else if (!itags.len) {
         if (B.itags.len!=R) wblog(FL,
            "ERR %s() invalid itag %s",FCT, STR(B.itags));
         for (i=0; i<R; ++i) {
            C.itags[2*i]=C.itags[2*i+1]=B.itags[i];
         }
      }
      else {
         if (itags.len!=B.itags.len || itags.len!=R || !sameQDir(FL,B))
            wblog(FL,"ERR %s() itag/qdir mismatch (%s <> %s)",FCT,
            STR(itags),STR(B.itags));
         for (i=0; i<R; ++i) {
            C.itags[2*i  ]=  itags[i];
            C.itags[2*i+1]=B.itags[i];
         }
      }
   }

   for (i=ib=0; ib<n2; ++ib) {
      qb=B.QIDX.rec(ib);

      for (ia=0; ia<n1; ++ia, ++i) {
         qa=QIDX.rec(ia); q=C.QIDX.rec(i);

         for (l=r=0; r<R; ++r, l+=QDIM) {
            memcpy(q,qa+l,s); q+=QDIM;
            memcpy(q,qb+l,s); q+=QDIM;
         }
         DATA[ia]->tensorProd(*B.DATA[ib], *C.DATA[i]);
      }
   }

   return C;
};

template <class TQ, class TD> inline
void QSpace<TQ,TD>::GroupIndizes(unsigned K) {
   unsigned i, r=-1;

   if (!isConsistent(r)) wbdie(FL,str);

   if (!K || unsigned(r)%K) wblog(FL,
   "ERR %s() cannot block QSpace (%d/%d)",FCT,r,K);

   QDIM=K*QDIM;

   for (i=0; i<DATA.len; ++i)
   DATA[i]->GroupIndizes(K);
}

template <class TQ, class TD>
size_t QSpace<TQ,TD>::getDataSize(wbMatrix<size_t> &S) const {

   size_t s=0; unsigned i,m,n=DATA.len;

   if (n!=QIDX.dim1) wblog(FL, 
      "ERR size mismatch (%d/%d)",QIDX.dim1,DATA.len);

   if (isEmpty()) { S.init(); return s; }

   m=DATA[0]->SIZE.len; S.init(n,m);

   for (i=0; i<n; ++i) {
      const wbvector<size_t> &si=DATA[i]->SIZE;
      if (si.len!=m) wblog(FL,"ERR size mismatch (%d/%d)",si.len,m);
      S.recSetP(i,si.data);
      s+=si.prod(0);
   }
   return s;
};

template <class TQ, class TD> inline
size_t QSpace<TQ,TD>::getDataSize(wbvector<size_t> &S) const {
   size_t s=0; S.init(DATA.len);
   for (unsigned i=0; i<DATA.len; ++i) { s+=(S[i]=DATA[i]->SIZE.prod()); }
   return s;
};

template <class TQ, class TD> inline
size_t QSpace<TQ,TD>::getDataSize() const {
   size_t s=0;
   for (unsigned i=0; i<DATA.len; ++i) s+=DATA[i]->SIZE.prod();
   return s;
};

template <class TQ, class TD> inline
size_t QSpace<TQ,TD>::getCGSSize() const {
   size_t s=0; 
   if (gotCGS(FL)>0) {
      unsigned i,j;
      for (i=0; i<CGR.dim1; ++i) {
      for (j=0; j<CGR.dim2; ++j) s+=CGR(i,j).numel(); }
   }
   return s;
};

template <class TQ, class TD> inline
wbstring QSpace<TQ,TD>::totSize2Str() const {

   unsigned long ss[3]={
       sizeof(TD)*getDataSize(),
       sizeof(TQ)*(QIDX.dim1*QIDX.dim2),
       sizeof(double)*getCGSSize()
   };
   wbvec<char> s_(32);
   Wb::memsize2Str(ss[0]+ss[1]+ss[2],s_.data,s_.len);
   return s_.data;
};

template <class TQ, class TD> inline
bool QSpace<TQ,TD>::isZero(const double eps, char flag) const {

   for (unsigned i=0; i<DATA.len; ++i)
   if (!(DATA[i]->isZero(eps,flag))) return 0;

   return 1;
}

template <class TQ, class TD> inline
void QSpace<TQ,TD>::Set(
    unsigned k, const wbvector<TQ> &Q, const wbarray<TD> &D
){

   QIDX.recSet(k,Q); 

   if (k>=DATA.len)
   wblog(FL,"ERR %s() index out of bounds (%d/%d)",FCT,k,DATA.len);

   DATA[k]->init(D);
}

template <class TQ, class TD> inline
void QSpace<TQ,TD>::SetQ(unsigned i, unsigned k, const TQ *d){

   if (i>=QIDX.dim1 || k>=rank(FL)) wblog(FL,
   "ERR index out of bounds (%d/%d; %d;%d)",i,QIDX.dim1,k,rank(FL));

   memcpy(QIDX.data+QIDX.dim2*i+QDIM*k, d, QDIM*sizeof(TQ));
};

template <class TQ, class TD> inline
QSpace<TQ,TD>& QSpace<TQ,TD>::SetFDir(
   const char *F, int L, const QSpace<TQ,TD> &A, const QSpace<TQ,TD> &B) {

   if (A.fdir) { 
      if (B.fdir && A.fdir!=B.fdir) wblog(F_L,"ERR %s() " 
         "fdir mismatch %s <> %s",FCT,STR(A.fdir),STR(B.fdir));
      fdir=A.fdir;
   }
   else { fdir=B.fdir; }
   return *this;
};

template <class TQ, class TD> inline
QSpace<TQ,TD>& QSpace<TQ,TD>::SetFDir_op(
   const char *F, int L, const iFerm &fdr, const wbindex &I) {

   if (!fdr) {
      if (fdir.data) { fdir.init(); }
   }
   else {
      unsigned r=rank(FL);
      char f=fdr.get(FL,I); 
      if (r==2) { fdir.init_op(F_L,"pm", f,fdr); } else
      if (r==3) { fdir.init_op(F_L,"pm-",f,fdr); }
      else wblog(FL,"ERR %s() got rank r=%d",FCT,r);
   }
   return *this;
};

template <class TQ, class TD> inline
QSpace<TQ,TD>& QSpace<TQ,TD>::SetFDir_Atensor(const char *F, int L,
   const iFerm &fdrA, const wbindex &Ia,
   const iFerm &fdrB, const wbindex &Ib) {

   if (!fdrA || !fdrB) {
      if (fdrA.data || fdrB.data) wblog(FL,
         "ERR %s() fdir inconsistency (%s / %s)",FCT,STR(fdrA),STR(fdrB));
      if (fdir.data) { fdir.init(); }
   }
   else {
      char fa=fdrA.get(FL,Ia), fb=fdrB.get(FL,Ib);

      unsigned r=rank(FL); 
      if (r!=3) wblog(FL,"ERR %s() got rank-%d tensor",FCT,r);
      fdir.init_A(F_L,"pp-",fa,fb,fdrA);
   }
   return *this;
};

template <class TQ, class TD> inline
void QSpace<TQ,TD>::Sort(const QSpace<TQ,TD> &R) {
   int e;
   wbindex P,Ir; wbperm Pr;

   if (QIDX.dim1==0) return;

   e=matchIndex(QIDX,R.QIDX, P, Ir);

   if (e || P.len!=QIDX.dim1) wblog(FL,
      "ERR QSpace R not complete to be used as sorting reference "
      "(%d/%d)",P.len, QIDX.dim1);

   Ir.Sort(Pr); P.Select(Pr);

   QIDX.recPermute(P); if (!CGR.isEmpty())
   CGR .recPermute(P);
   DATA.Select(P);
};

template <class TQ, class TD> inline
int QSpace<TQ,TD>::isSingleton(unsigned k) const {
   unsigned i;

   if (k>=rank(FL)) wblog(FL,
   "ERR %s() index out of bounds (%d,%d)",FCT,k,rank(FL));

   for (i=0; i<DATA.len; ++i)
   if (DATA[i]->SIZE[k]!=1) break;

   return (i==DATA.len);
}

template <class TQ, class TD> inline
void QSpace<TQ,TD>::getDRange(TD &dmin, TD &dmax) const {
    unsigned i,n=DATA.len;

    if (!isConsistent()) { info("this"); wbdie(FL,str); }

    for (i=0; i<n; ++i) if (DATA[i] && !DATA[i]->isEmpty()) {
       dmin=dmax=DATA[i]->data[0];
       break;
    }

    for (i=0; i<n; ++i) if (DATA[i] && !DATA[i]->isEmpty()) {
       dmin=MIN(dmin, DATA[i]->min());
       dmax=MAX(dmax, DATA[i]->max());
    }
}

template <class TQ, class TD> inline
QSpace<TQ,TD>& QSpace<TQ,TD>::getSub(
   const wbindex &I, QSpace<TQ,TD> &A, char ref) const {

   if ((void*)this==(void*)&A) { QSpace<TQ,TD> X;
      if (  ref) wblog(FL,"ERR %s() requesting ref=%s on self",FCT,cSTR(ref));
      if (isref) wblog(FL,"ERR %s() got ref QSpace",FCT);
      return getSub(I,X).save2(A); 
   }

   unsigned i,j;
   char cgflag=gotCGS(FL);

   A.init();
   QIDX.getRecs(I,A.QIDX); A.QDIM=QDIM; A.qtype=qtype; A.itags=itags;

   A.isref=(ref ? 1 : 0);

   A.setupDATA(); 

   if (cgflag>0) {
      if (CGR.dim1!=DATA.len) wblog(FL,
         "ERR %s() size mismatch (CGR: %d/%dx%d)",
          FCT,DATA.len,CGR.dim1,CGR.dim2);
      A.setupCGR(); 
   }

   for (i=0; i<I.len; ++i) {
      if (I[i]>=DATA.len) wblog(FL, 
         "ERR %s() index out of bounds (%d,%d)",FCT,I[i],DATA.len);

      if (A.isref)
           { A.DATA[i]=DATA[I[i]]; } 
      else { A.DATA[i]->init(*DATA[I[i]]); }

      if (cgflag>0) {
         for (j=0; j<CGR.dim2; ++j) {
            A.CGR(i,j)=CGR(I[i],j);
         }
      }
   }

   return A;
};

template <class TQ, class TD> inline
QSpace<TQ,TD>& QSpace<TQ,TD>::getSubInit(
    const wbindex &I,
    QSpace<TQ,TD> &A
) const {

    if (this==&A) { 
       QSpace<TQ,TD> X; A.save2(X);
       return X.getSubInit(I,A);
    }

    unsigned i,j;
    char cgflag=gotCGS(FL);

    A.init();
       QIDX.getRecs(I,A.QIDX); A.QDIM=QDIM; A.qtype=qtype;
    A.setupDATA();

    if (cgflag) { A.setupCGR();

       for (i=0; i<I.len; ++i) {
          if (I[i]>=QIDX.dim1) wblog(FL, 
             "ERR %s() index out of bounds (%d/%d)",FCT,I[i]+1,QIDX.dim1);
          for (j=0; j<CGR.dim2; ++j) A.CGR(i,j).init(CGR(I[i],j));
       }
    }

    return A;
};

template <class TQ, class TD> inline
void QSpace<TQ,TD>::saveSub2(QSpace<TQ,TD> &A, const wbindex &I) {

    if (this==&A) wblog(FL,
       "ERR %s() got overlapping input and output space!",FCT);
    if (isref) wblog(FL,"ERR %s() called for ISREF!",FCT);

    if (!I.len) { A.init(); return; }

    unsigned i,j;
    char cgflag=gotCGS(FL);

    A.init();
    QIDX.getRecs(I,A.QIDX); A.QDIM=QDIM; A.qtype=qtype;
    A.setupDATA();
    A.setupCGR();

    for (i=0; i<I.len; ++i) {
       if (I[i]>=DATA.len) wblog(FL, 
          "ERR %s() index out of bounds (%d,%d)",FCT,I[i],DATA.len);

       DATA[I[i]]->save2(*A.DATA[i]);
       WB_DELETE_1(DATA[I[i]]);

       if (cgflag) {
          for (j=0; j<CGR.dim2; ++j) {
          CGR(I[i],j).save2(A.CGR(i,j)); }
       }
    }

    wbindex I2; I.invert(QIDX.dim1,I2);

    QIDX.Set2Recs(I2);
    DATA.Select(I2); if (cgflag) { CGR.Set2Recs(I2); }
};

template <class TQ, class TD>
void QSpace<TQ,TD>::toFull(
   wbarray<TD> &A,                    
   wbvector< wbMatrix<TQ> > &QB,      
   wbvector< wbvector<unsigned> > &SB 
) const {

   if (gotCGS(FL)>0 || gotCGX()) wblog(FL,
      "ERR %s() not implemented yet for %d %d",FCT,
      gotCGS(FL), gotCGX(FL)
   );

   widx_t i,j,k,M=QIDX.dim1; unsigned r=-1;
   wbvector<widx_t> N,SA,D;
   wbarray<unsigned> mark;
   wbMatrix<widx_t> G;
   wbindex I,J;
   wbperm P,iP;

   if (!isConsistent(r)) { wblog(FL,"%s",str); }
   if (isEmpty()) {
      A.init(); QB.init(); SB.init();
      return;
   }

   QB.init(r); SB.init(r); D.init(r);
   I.init(r); P.init(r); iP.init(r); G.init(M,r); N.init(r);

   for (k=0; k<r; ++k) {
      getQsub(k,QB[k]).groupRecs(P,D); P.invert(iP);

      I.BlockIndex(D).Permute(iP); 
      G.setCol(k,I.data);

      N[k]=D.len; SB[k].init(N[k]);
   }

   mark.init(N+1); 

   for (i=0; i<M; ++i) {
      G.getRec(i,I);

      if (mark(I)++) wblog(FL,
         "ERR %s() block (%s) not unique (%d)",
          FCT,(I+1).toStrf("",",").data, mark(I)
      );

      J.wbvector<widx_t>::operator=(N);

      for (k=0; k<r; ++k) {
         wbvector<unsigned> &SBk=SB[k];

         J[k]=j=I[k]; if (k) J[k-1]=N[k-1];
         if (!mark(J)++) {
            SBk[j]=DATA[i]->SIZE[k];
         }
         else if (SBk[j]!=DATA[i]->SIZE[k]) wblog(FL,
           "ERR %s - block size mismatch (%d: %d,%d)",
            FCT,i+1, SBk[j], DATA[i]->SIZE[k]
         );
      }
   }

   SA.init(r); for (k=0; k<r; ++k) SA[k]=SB[k].sum();
   A.init(SA);

   for (i=0; i<DATA.len; ++i) {
      G.getRec(i,I);
      A.setBlock(SB, I, *DATA[i]);
   }
}

template <class TQ, class TD>
void QSpace<TQ,TD>::toFull2(
   wbarray<TD> &A,        
   wbMatrix<TQ> &QB,      
   wbMatrix<widx_t> &SB   
) const {

   if (gotCGS(FL)>0 || gotCGX()) wblog(FL,
      "ERR %s() not implemented yet for %d %d",FCT,
      gotCGS(FL), gotCGX(FL)
   );

   unsigned i,j,k,r,K,N,M=QIDX.dim1, r_=-1;
   wbMatrix<TQ> Q1,Q2;
   wbvector< wbvector<widx_t> > S2;
   wbvector<widx_t> s,D,SA;
   wbMatrix<unsigned> mark; wbMatrix<widx_t> IJ;
   wbindex I,I1,I2;
   wbarray<TD> a;
   wbperm P,iP;

   if (!isConsistent(r_)) wblog(FL,"%s",str);
   if (isEmpty()) {
      A.init(); QB.init(); SB.init();
      return;
   }

   if (r_==3 && getDim(2)==1) { --r_; }

   if (r_%2) wblog(FL,"ERR %s even-rank object required (r=%d)", FCT,r_);
   r=unsigned(r_);
   K=r/2;

   I.Index(0,K-1); getQsub(I,Q1);
   I.Index(K,r-1); getQsub(I,Q2);

   QB.Cat(1,Q1,Q2);

   QB.groupRecs(P,D); P.invert(iP);

   I.BlockIndex(D); I.Select(iP);

   I1.init(M, I.data);   
   I2.init(M, I.data+M); 

   N=D.len;
   mark.init(N,N+1);     
   SB.init(N,K);         
   IJ.init(M,2);

   for (k=0; k<M; ++k) {
      i=IJ(k,0)=I1[k]; j=IJ(k,1)=I2[k];

      if (mark(i,j)++) wblog(FL, 
      "WRN %s - block (%d,%d) not unique (%d)", FCT, i,j, mark(i,j));

      s.init2ref(K, DATA[k]->SIZE.data);
      if (mark(i,N)++) {
         if (!SB.recEqual(i,s.data)) wblog(FL,
         "ERR %s - block size mismatch (%d: %s != %s)",FCT,k+1,
         SB.rec2Str(i,"","x").data, SSTR(s));
      }
      else SB.recSet(i,s);

      s.init2ref(K, DATA[k]->SIZE.data+K);
      if (mark(j,N)++) {
         if (!SB.recEqual(j,s.data)) wblog(FL,
         "ERR %s - block size mismatch (%d: %s != %s)",FCT,k+1,
         SB.rec2Str(j,"","x").data, SSTR(s));
      }
      else SB.recSet(j,s);
   }

   S2.init(2); SB.recProd(S2[0]); S2[1]=S2[0];
   SA.init(2); SA[0]=SA[1]=S2[0].sum();

   A.init(SA);

   for (i=0; i<DATA.len; ++i) {
      I.init2ref(2,IJ.rec(i));

      a.init2ref(*DATA[i]); 
      a.GroupIndizes(K);

      A.setBlock(S2,I,a);
   }
}

template <class TQ, class TD>
void QSpace<TQ,TD>::ExtendQ(const QSpace<TQ,TD> &B) {

   unsigned i,k,m,n=B.QIDX.dim1;
   wbvector<char> mark;
   wbindex Ia,Ib;

   if (B.isEmpty()) return;
   if (isEmpty()) {
      (*this)=B;
      for (i=0; i<DATA.len; ++i) DATA[i]->set(TD(0));
      return;
   }

   if (QIDX.dim2!=B.QIDX.dim2) wblog(FL,
   "ERR %s - size mismatch (%d/%d)", FCT, QIDX.dim2, B.QIDX.dim2);

   matchIndex(QIDX,B.QIDX,Ia,Ib);

   mark.init(n);
   for (i=0; i<Ib.len; ++i) mark[Ib[i]]++;

   m=mark.count(0);
   if (m==0) return; 

   k=QIDX.dim1;
   Append(m);

   for (i=0; i<n; ++i) {
      if (mark[i]) continue;
      QIDX.recSetP(k,B.QIDX.rec(i));
      DATA[k]->init(B.DATA[i]->SIZE);
      ++k;
   }
};

template <class TQ, class TD>
void QSpace<TQ,TD>::takeMainQS(char disp) {

   unsigned i,j,n,i0,imax=0;
   wbvector<widx_t> d, dc, I;
   wbperm P;

   wbMatrix<TQ> QT;
   wbvector<TD> Z;
   TD z, zmax=0, ztot=0, eps=1e-3;

   getQtot(QT);
   QT.groupRecs(P,d);

   d.cumsum_(dc); Z.init(d.len);

   for (i=0; i<d.len; ++i) {
      for (z=0, n=d[i], i0=dc[i], j=0; j<n; ++j)
      z+=DATA[P[i0+j]]->norm2();

      Z[i]=z;

      if (z>zmax) { zmax=z; imax=i; }
      ztot+=z;
   }

   if (disp) { char s_[128];
      z=fabs((ztot-zmax)/ztot);
      snprintf(s_,128,"%s skipping weight %s of total %s (%s%s)\n"
      "    Q=[%s]", (ABS(z)>eps ? "WRN" : " * "),
           NSTR(ztot-zmax), NSTR(ztot), NSTR(z*100.),
           ABS(z)>0.1 ? "% !?" : "%", QT.rec2Str(imax).data
      );

      wblog(FL,s_);
   }

   I.init(d[imax], P.data+dc[imax]);

   QIDX.Set2Recs(I);
   DATA.Select(I);
}

template <class TQ, class TD> inline
wbMatrix<TQ>& QSpace<TQ,TD>::getQsub(
    const unsigned k, wbMatrix<TQ> &Qk, char iflag
) const {

    unsigned i, r=rank(FL), s=QDIM*sizeof(TQ), offset=k*QDIM;

    if (iflag) { 
       wbindex I; I.Index_ex(r,k); getQsub(I,Qk);
       return Qk;
    }

    if (k>=r) wblog(FL,
       "ERR %s() index out of bounds (%dx%d/%d; %d)",
        FCT, QIDX.dim1, QIDX.dim2, QDIM, k);

    Qk.init(QIDX.dim1, QDIM);

    for (i=0; i<QIDX.dim1; ++i)
    memcpy(Qk.rec(i), QIDX.rec(i)+offset, s);

    return Qk;
};

template <class TQ, class TD> inline
wbMatrix<TQ>& QSpace<TQ,TD>::getQsub(
    const wbindex &I, 
    wbMatrix<TQ> &QI  
) const {
    unsigned i,k, r=rank(FL);
    TQ *q;

    if (&QI==&QIDX) wblog(FL,
    "ERR %s() QIN and QOUT space are the same!",FCT);
    for (i=0; i<I.len; ++i) if (I[i]>=r) wblog(FL,
    "ERR %s() index out of bounds (%d/%d)",FCT,I[i],r);

    QI.init(QIDX.dim1, I.len*QDIM);

    for (i=0; i<QIDX.dim1; ++i) {
        const TQ *q0=QIDX.rec(i); q=QI.rec(i);

        for (k=0; k<I.len; ++k, q+=QDIM) {
           Wb::cpyRange(q, q0+I[k]*QDIM, QDIM);
        }
    }

    return QI;
}

template <class TQ, class TD> inline
void QSpace<TQ,TD>::getQsub(
    const wbindex &I,  
    wbMatrix<TQ> &QI,  
    wbMatrix<TQ> &Q2   
) const {

    unsigned i, k, r=rank(FL), s=QDIM*sizeof(TQ);
    wbindex I2;
    TQ *q1,*q2; const TQ *q0;

    I.invert(r,I2); 

    if (&QI==&QIDX || &Q2==&QIDX)
    wblog(FL,"ERR getQsub requires distinct QOUT space!");

    QI.init(QIDX.dim1, I .len*QDIM);
    Q2.init(QIDX.dim1, I2.len*QDIM);

    for (i=0; i<QIDX.dim1; ++i) {
        q0=QIDX.rec(i); q1=QI.rec(i); q2=Q2.rec(i);

        for (k=0; k<I .len; ++k) memcpy(q1+k*QDIM, q0+I [k]*QDIM, s);
        for (k=0; k<I2.len; ++k) memcpy(q2+k*QDIM, q0+I2[k]*QDIM, s);
    }
}

template <class TQ, class TD> inline
wbMatrix<TQ>& QSpace<TQ,TD>::getQtot(
    wbMatrix<TQ> &Qtot,
    char uflag 
) const {

    unsigned i,k,l,r, r0=rank(FL);
    TQ *q; const TQ *q0;

    if (&Qtot==&QIDX)
    wblog(FL,"ERR getQtot requires distinct QOUT space!");

    Qtot.init(QIDX.dim1, QDIM);

    for (i=0; i<QIDX.dim1; ++i) {
        q0=QIDX.rec(i); q=Qtot.rec(i);
        for (l=r=0; r<r0; ++r)
        for (k=0; k<QDIM; ++k, ++l) q[k]+=q0[l];
    }

    if (uflag) Qtot.makeUnique();

    return Qtot;
}

template <class TQ, class TD> inline
void QSpace<TQ,TD>::getDistQtot(
   wbMatrix<TQ> &Q, wbvector<double> &w
) const {
   unsigned i,j,l,d;
   wbperm P; wbvector<widx_t> D;
   double n2;

   QIDX.blockSum(QDIM,Q); Q.groupRecs(P,D);
   w.init(D.len);

   for (l=i=0; i<D.len; ++i, l+=d) { d=D[i];
      for (n2=0, j=0; j<d; ++j) n2 += (DATA[P[l+j]]->norm2());
      w[i]=sqrt((double)n2);
   }
}

template <class TQ, class TD> inline
wbMatrix<TQ>& QSpace<TQ,TD>::getQsum(
  const wbindex &J, wbMatrix<TQ> &Q
) const {

    unsigned i,j, r=rank(FL);
    const TQ *q0; TQ *q;

    if (J.isEmpty()) { Q.init(); return Q; }

    for (i=0; i<J.len; ++i) if (J[i]>=r) wblog(FL, 
       "ERR block index out of range (%d/%d)",J[i],r);

    Q.init(QIDX.dim1,QDIM);

    if (allAdditive() || J.len==1) {
       for (i=0; i<QIDX.dim1; ++i) {
          q=Q.rec(i); q0=QIDX.rec(i);
          for (j=0; j<J.len; ++j) Wb::addRange(q0+J[j]*QDIM, q, QDIM);
       }
    }
    else if (allAbelian()) {
       if (qtype.len!=QDIM) wblog(FL,
          "ERR %s() got inconsistent qtype (len=%d/%d; %s)",
          FCT,qtype.len,QDIM,qStr().data
       );

       for (i=0; i<QIDX.dim1; ++i) { q=Q.rec(i); q0=QIDX.rec(i);
          for (j=0; j<QDIM; ++j) {
             q[j]=get_qtot_abelian(qtype[j],q0+j,J.len,J.data,QDIM);
          }
       }
    }
    else wblog(FL,
       "ERR %s() got non-abelian symmetries (%s)",FCT,qStr().data);

    return Q;
};

template <class TQ, class TD> inline
wbMatrix<TQ>& QSpace<TQ,TD>::getQsum(
  const wbMatrix<TQ> &Q0, wbMatrix<TQ> &Q
) const {

    if (!Q0) { Q.init(); return Q; }
    if (Q0.dim2==QDIM) { Q=Q0; return Q; }

    if (!QDIM || (Q0.dim2 % QDIM)) wblog(FL, 
       "ERR size mismatch Q0 (%dx%d /%d)",Q0.dim1,Q0.dim2,QDIM);

    unsigned i,j, m=Q0.dim2/QDIM;
    const TQ *q0; TQ *q;

    Q.init(Q0.dim1,QDIM);

    if (allAdditive()) {
       for (i=0; i<Q0.dim1; ++i) { q=Q.rec(i); q0=Q0.rec(i);
          for (j=0; j<m; ++j) Wb::addRange(q0+j*QDIM, q, QDIM);
       }
    }
    else if (allAbelian()) {
       if (qtype.len!=QDIM) wblog(FL, 
          "ERR %s() got inconsistent qtype (len=%d/%d; %s)",
          FCT,qtype.len,QDIM,qStr().data
       );

       for (i=0; i<Q0.dim1; ++i) { q=Q.rec(i); q0=Q0.rec(i);
          for (j=0; j<QDIM; ++j) {
             q[j]=get_qtot_abelian(qtype[j],q0+j,m,nullptr,QDIM);
          }
       }
    }
    else wblog(FL,"ERR %s() "
        "got non-abelian symmetries (%s, m=%d)",FCT,qStr().data,m);

    return Q;
};

template <class TQ, class TD> inline
int QSpace<TQ,TD>::getFermZ2(
   wbMatrix<unsigned> &Z2, wbvector<unsigned> *z2_) const {

   int nferm=0; 
   if (!fdir || !DATA.len || isScalar()) {
      Z2.init(); if (z2_) { z2_->init(); }
      return nferm;
   }

   unsigned i,j,k, n, r=rank(FL); widx_t *iF; int z2;
   const TQ *qidx=QIDX.data; QType *qf;
   wbvector<QType> Qf;
   wbindex If;

   Z2.init2val(DATA.len,r,+1);

   n=fdir.iSym(Qf,If,qtype, qtype? -1 : QDIM);

   if (!r || r!=fdir.r) wblog(FL,
      "ERR %s() rank inconsistency r=%d/%d",FCT,r,fdir.r);
   if (!n || n!=QDIM) wblog(FL,
      "ERR %s() got QDIM=%d/%d",FCT,n,QDIM);
   n=Qf.len; iF=If.data; qf=Qf.data;

   for (i=0; i<DATA.len; ++i) {
      for (j=0; j<r; ++j, qidx+=QDIM) { z2=0;
         for (k=0; k<n; ++k) { z2 += qf[k].getZ2(qidx[iF[k]]); }
         if (( Z2(i,j) = Wb::mod(z2,2) )) { ++nferm; }
      }
   }

   if (!nferm) {
      Z2.init(); if (z2_) { z2_->init(); }
   }
   else if (z2_) {
      unsigned nf=0, *z2;     
      Z2.recSum(*z2_); n=z2_->len; z2=z2_->data;
      for (i=0; i<n; ++i) {
         if (z2[i]%2) {
            MXPut(FL,"Iz2","base").add(Z2,"Z2").add(*z2_,"z2").add(*this,"A");
            wblog(FL,"ERR %s() tensor not parity preserving (i=%d/%d)",FCT,i,n);
         }
         if (( z2[i] = (z2[i]/2)%2 )) { ++nf; }
      }
      if (!nf) { z2_->init(); } 
   }

   #if 0
    MXPut Idbg(FL,"Idbg","base"); Idbg.add(QIDX,"QIDX").add(QDIM,"QDIM")
       .add(Qf,"Qf").add(If,"If").add(fdir,"fdir").add(Z2,"Z2")
       .add(nferm,"nferm");
    if (z2_) { Idbg.add(*z2_,"z2"); }
   #endif

   return nferm;
};

template <class TQ, class TD> inline
int QSpace<TQ,TD>::getFermZ2(
   wbperm &P, MVEC &z2, const ctrIdx *ic_, const wbindex *iz_) const {

   unsigned i,j,r=rank(FL); int nferm=-1;
   wbMatrix<unsigned> Z2,T2;
   wbvector<unsigned> zx;

   if (!P && (!ic_ || !*ic_)  && (!iz_ || !*iz_)) {
      z2.init(); return 0;
   }
   nferm=getFermZ2(Z2, P.conj ? &zx : nullptr);

   if (nferm<0) wblog(FL,"ERR %s() got nferm=%d",FCT,nferm); else
   if (!nferm) { z2.init();
if (Wb::envFERM) wblog(FL,"TST P=%s (DATA.len=%d) => nferm=%d",STR(P),DATA.len,nferm);
      return 0; }

   z2.init(DATA.len);

   if (P.len>Z2.dim2) wblog(FL,"ERR %s() permutation out of range "
      "(len=%d/%d/%d)",FCT,P.len,Z2.dim2,rank(FL));

   P.getTranspositionsNN(T2);

if (Wb::envFERM) wblog(FL,"TST P=%s (DATA.len=%d)",STR(P),DATA.len);
   if (T2) { unsigned *t2, *z2i=Z2.data;
      for (i=0; i<DATA.len; ++i, z2i+=Z2.dim2) { t2=T2.data;
         for (j=0; j<T2.dim1; ++j, t2+=T2.dim2) {
            z2[i] += (z2i[t2[0]] * z2i[t2[1]]);
         }
      }
if (Wb::envFERM) wblog(FL,"TST 1) zf=[%8s]  T2=[%s]",STR(z2),STR(T2));
   }

   if (iz_ && *iz_) {
      unsigned n=iz_->len, *z2i=Z2.data;
      size_t *iz = iz_->data;
      if (!iz_->isUnique(r)) wblog(FL,
         "ERR %s() invalid index %s",FCT,STR_(iz_));

      for (i=0; i<DATA.len; ++i, z2i+=Z2.dim2) {
         for (j=0; j<n; ++j) { z2[i] += z2i[iz[j]]; }
      }
if (Wb::envFERM) wblog(FL,"TST 2) zf=[%8s]  iz=%s",STR(z2),STR_(iz_));
   }

   if (ic_ && ic_->len) { 
      unsigned *ic=ic_->data, n=ic_->len, *z2i=Z2.data;
      int z, z2_, zj; ic_->checkUnique(FL,r);

      for (i=0; i<DATA.len; ++i, z2i+=Z2.dim2) {
         for (z=z2_=0, j=0; j<n; ++j) { zj=z2i[ic[j]]; z+=zj;  z2_+=zj*zj; }
         z=z*z-z2_; 
         if (z%2) wblog(FL,"ERR %s() unexpected parity (%d^2-%d)/2",FCT,z,z2_);
         z2[i]+=(z/2);
      }
if (Wb::envFERM) wblog(FL,"TST 3) zf=[%8s]  ic=%s",STR(z2),STR_(ic_));
   }

   if (P.conj && zx) { 
      if (zx.len!=z2.len) wblog(FL,
         "ERR %s() size inconsistency (len=%d/%d)",FCT,zx.len,z2.len);

      for (i=0; i<zx.len; ++i) { z2[i]+=(zx[i]%2); }
if (Wb::envFERM) wblog(FL,"TST 4) zf=[%8s]  zx=[%s] since P.conj",STR(z2),STR(zx));
   }

   for (nferm=i=0; i<z2.len; ++i) { if ((z2[i]%=2)) { ++nferm; }}
   if (!nferm) { z2.init(); }

if (Wb::envFERM) { wblog(FL,"TST -> zf=[%8s]  nferm=%d%N",STR(z2),nferm); }
   return nferm;
};

template <class TQ, class TD> inline
int QSpace<TQ,TD>::getFermZ2_conj(MVEC &z2) const {

   int nferm=0; 
   wbMatrix<unsigned> Z2;
   wbvector<unsigned> zx;

   if ((nferm=getFermZ2(Z2,&zx))<=0 || !zx) {
      z2.init(); return nferm;
   }

   if (zx.len!=DATA.len) wblog(FL,
      "ERR %s() size inconsistency (len=%d/%d)",FCT,zx.len,z2.len);
   z2.init(DATA.len);

   nferm=0;
   for (unsigned i=0; i<zx.len; ++i) { 
      if ((z2[i]=(zx[i]%2))) { ++nferm; }
   }
   if (!nferm) { z2.init(); }

   return nferm;
};

template <class TQ, class TD> inline
void QSpace<TQ,TD>::FlipQ() {
   wbperm P(rank(FL)); P.Flip();
   PermuteQ(P);
}

template <class TQ, class TD> inline
void QSpace<TQ,TD>::PermuteQ(const wbperm &P) {

   unsigned i,k,r=rank(FL), s=QDIM*sizeof(TQ);
   wbMatrix<TQ> Q0(QIDX);
   TQ *q0, *qp;

   if (!QDIM) wblog(FL,"WRN QDIM=0!");
   if (P.isEmpty()) { return; }
   if (isEmpty() || !QIDX.dim2) { return; }

   if (int(i=P.isValidPerm(r,1))<=0) wblog(FL, 
      "ERR %s() invalid permutation [%s] (r=%d, e=%d)",FCT,STR(P),r,i);

   for (i=0; i<QIDX.dim1; ++i) {
       q0=Q0.rec(i); qp=QIDX.rec(i);
       for (k=0; k<P.len; ++k) { memcpy(qp+k*QDIM, q0+P[k]*QDIM, s); }
   }
};

template <class TQ, class TD> inline
void QSpace<TQ,TD>::permuteQ(wbMatrix<TQ> &Q, const wbperm &P) const {

   if (&Q==&QIDX) {
      wblog(FL,"ERR do not permuteQ onto itself");
   }
   else {
      unsigned i, k, r=rank(FL), s=QDIM*sizeof(TQ);
      const TQ *q0; TQ *qp;

      if (P.isEmpty()) { Q=QIDX; return; }

      if (int(i=P.isValidPerm(r,1))<=0) wblog(FL,
         "ERR invalid permutation '%s' (e=%d)",STR(P),i);

      Q.init(QIDX.dim1, QIDX.dim2);

      for (i=0; i<QIDX.dim1; ++i) {
          q0=QIDX.rec(i); qp=Q.rec(i);
          for (k=0; k<P.len; ++k) memcpy(qp+k*QDIM, q0+P[k]*QDIM, s);
      }
   }
}

template <class TQ, class TD>
QSpace<TQ,TD>& QSpace<TQ,TD>::Permute( 
   wbperm P, 
   char rcpy 
 ) {         

   char q=P.relevant();
   unsigned r=rank(FL);

   if (!q || (!r && isEmpty())) { return *this; }

   if (!(q&1)) { 
      if (P.conj  ) { Conj(); }
      if (P.fac!=1) { (*this)*=(P.inv ? (1./P.fac) : P.fac); }
      return *this;
   }

   if (!r) wblog(PFL,
      "ERR %s() invalid permutation %s for scalar (r=0)",FCT,STR(P));
   if (P.len>r) wblog(PFL, 
      "ERR %s() permutation %s out of bounds (r=%d)",FCT,STR(P),r);

   unsigned i,n;
   double P_fac=P.fac;
   MVEC z2;

   P.flatten();

   this->getFermZ2(P,z2);

   PermuteQ(P); 
   if (!P_fac) wblog(FL,"WRN %s() got fac = %g",FCT,P.fac); 

   for (i=0; i<DATA.len; ++i) {
      if (Wb::envFERM && z2) { P.fac = P_fac*(z2[i]%2 ? -1 : 1); }
      DATA[i]->Permute(P,rcpy);
   }

   if (fdir ) { fdir.Permute(P);  } 
   if (itags) { itags.Permute(P); } 

   if (CGR) {
      for (n=CGR.numel(), i=0; i<n; ++i) {
         CGR[i].Permute(P);
      }
      NormCGW(0,0,rcpy); 
   }

   return *this;
};

template <class TQ, class TD> inline
QSpace<TQ,TD>& QSpace<TQ,TD>::permute( 
   QSpace<TQ,TD> &B, wbperm P) const { 

   char q=P.relevant();
   unsigned r=rank(FL);

   if (!q) { return B.init(*this); }

   if (!(q&1)) { B.init(*this); 
      if (P.conj  ) { B.Conj(); }
      if (P.fac!=1) { B*=(P.inv ? (1./P.fac) : P.fac); }
      return B;
   }

   if (!r) wblog(PFL,
      "ERR %s() invalid permutation %s for scalar (r=0)",FCT,STR(P));
   if (P.len>r) wblog(PFL, 
      "ERR %s() permutation %s out of bounds (r=%d)",FCT,STR(P),r);

   unsigned i,n;
   double P_fac=P.fac;
   MVEC z2;

   P.flatten(); P.fac=1; B.init();

   this->getFermZ2(P,z2);

   fdir.permute(B.fdir,P); 

   permuteQ(B.QIDX,P);
   B.QDIM=QDIM; B.qtype=qtype;

   itags.permute(B.itags,P); 

   if (otype && P.isOpTranspose()) {
      B.otype=otype;
   }

   B.setupDATA();
   if (!P_fac) wblog(FL,"WRN %s() got fac = %g",FCT,P_fac);

   P.fac=P_fac; 
   for (n=DATA.len, i=0; i<n; ++i) {
      if (Wb::envFERM && z2) { P.fac=P_fac*(z2[i]%2 ? -1 : 1); }
      DATA[i]->permute(*B.DATA[i],P);
   }

   if (CGR.isEmpty()) { B.CGR.init(); }
   else {
      B.setupCGR(); n=CGR.numel(); P.fac=1;
      if (n!=B.CGR.numel()) wblog(FL,"ERR %s() got CGR size mismatch "
         "(%s <> %s)",FCT,CGR.sizeStr().data,B.CGR.sizeStr().data);

      for (i=0; i<n; ++i) {
         CGR[i].permute(B.CGR[i],P); 
      }
   }

   return B;
};

template <class TQ, class TD> inline
QSpace<TQ,TD>& QSpace<TQ,TD>::PermuteFirstTo(unsigned k) {

   unsigned r=rank(FL);
   if (k>=r) wblog(FL,"ERR %s() index out of range (%d/%d)",FCT,k,r);
   wbperm P; P.initFirstTo(k,r);
   return Permute(P);
};

template <class TQ, class TD> inline
QSpace<TQ,TD>& QSpace<TQ,TD>::permuteFirstTo(unsigned k, QSpace &B) const {

   if (this==&B) { return B.PermuteFirstTo(k); }
   else {
      unsigned r=rank(FL);
      if (k>=r) wblog(FL,"ERR %s() index out of range (%d/%d)",FCT,k,r);
      wbperm P; P.initFirstTo(k,r);
      return permute(B,P);
   }
};

template <class TQ, class TD> inline
QSpace<TQ,TD>& QSpace<TQ,TD>::PermuteLastTo(unsigned k) {
   unsigned r=rank(FL);
   if (k>=r) wblog(FL,"ERR %s() index out of range (%d,%d)",FCT,k,r);
   wbperm P; P.initLastTo(k,r);
   return Permute(P);
};

template <class TQ, class TD> inline
QSpace<TQ,TD>& QSpace<TQ,TD>::permuteLastTo(unsigned k, QSpace &B) const {

   if (this==&B) { return B.PermuteLastTo(k); }
   else {
      unsigned r=rank(FL);
      if (k>=r) wblog(FL,"ERR %s() index out of range (%d,%d)",FCT,k,r);
      wbperm P; P.initLastTo(k,r);
      return permute(B,P);
   }
};

template <class TQ, class TD> inline
char QSpace<TQ,TD>::permute_to(
   char pos, 
   ctrIdx &ic, QSpace<TQ,TD> &B,
   char ref  
 ) const {

   char q; 
   unsigned i=0, r=rank(FL);

   unsigned l=(pos==2 ? 0 : pos);  
   wbperm P;

   for (; i<2; ++i) {
      if (l==0) { P.init2Front(ic,r); } else
      if (l==1) { P.init2End  (ic,r); }
      else wblog(FL,"ERR %s() invalid pos=%s (i=%d)",FCT,cSTR(pos),i);

      q=P.relevant(); 
      if (!q) {
         pos=1+l; if (!i) { q|=8; }
         break;
      }
      else if (q>0) { l=(l+1)%2; } 
      else wblog(FL,"ERR %s() q=%d",FCT,q);
   }

   if ((q&7)) { unsigned n=ic.len;
      if (ic.conj) { P.conj=1; ic.conj=0; } 

if (Wb::envFERM) wblog(PFL,"TST %s() P=%-6s(ic=%2s, pos=%d)",pFCT,STR(P),STR(ic),pos);
      permute(B,P);

      if (pos==1) 
           { for (i=0;        i<n; ++i) { ic[i]=i;   }}
      else { for (i=0, l=r-n; i<n; ++i) { ic[i]=i+l; }}
   }
   else if (ref)
        { B.init2ref(*this);    }
   else { B.init(*this); q|=16; }

   return q;
};

template <class TQ, class TD> inline
void QSpace<TQ,TD>::setQ(unsigned k,
    const QSpace<TQ,TD> &B, unsigned k0
) {

    unsigned i, s=QDIM*sizeof(TQ), dd0=B.QIDX.dim2, dd2=QIDX.dim2;
    TQ *d0 = B.QIDX.data+k0*QDIM, *d2 = QIDX.data+k*QDIM;

    if (QIDX.dim1!=B.QIDX.dim1 || QDIM!=B.QDIM) wblog(FL,
       "ERR %s() incompatible objects (%d,%d; %d,%d)",
        FCT, QIDX.dim1, B.QIDX.dim1, QDIM, B.QDIM);

    if (rank(FL)<=k || B.rank(FL)<=k0) wblog(FL,
       "ERR %s() index out of bounds (%d/%d; %d/%d)",
        FCT, k, rank(FL), k0, B.rank(FL));

    for (i=0; i<QIDX.dim1; ++i, d0+=dd0, d2+=dd2)
    memcpy(d2,d0,s);
}

template <class TQ, class TD> inline
void QSpace<TQ,TD>::setQRec(
    unsigned r, unsigned j, const TQ *q0, const TQ *q2
){
    if (r>=QIDX.dim1 || QDIM*(j+1)>QIDX.dim2) wblog(FL,
       "ERR %s() index out of bounds (%d,%d; %d,%d)",
        FCT, r, j, QIDX.dim1, QDIM);

    TQ *q = QIDX.rec(r)+j*QDIM;
    if (q2==nullptr)
         for (unsigned i=0; i<QDIM; ++i) q[i]=q0[i];
    else for (unsigned i=0; i<QDIM; ++i) q[i]=q0[i]+q2[i];
}

template <class TQ, class TD> inline
bool QSpace<TQ,TD>::hasEqualQ(unsigned i1, unsigned i2) const {

    unsigned i, s=QDIM*sizeof(TQ), D=QIDX.dim2, r=rank(FL);
    TQ *d1 = QIDX.data+i1*QDIM, *d2 = QIDX.data+i2*QDIM;

    if (i1>=r || i2>=r) wblog(FL,
    "ERR %s() index out of bounds (%d,%d; %d)",FCT,i1,i2,r);

    for (i=0; i<QIDX.dim1; ++i, d1+=D, d2+=D)
    if (memcmp(d1,d2,s)) return 0;

    return 1;
};

template <class TQ, class TD>
bool QSpace<TQ,TD>::hasQOverlap(const QSpace<TQ,TD> &A) const {

    wbindex Ia, Ib;
    if (qtype!=A.qtype) {
       wblog(FL,"WRN got objects with different qtype");
       return 0;
    }

    matchIndex(QIDX, A.QIDX, Ia, Ib);
    return (Ia.len!=0);
};

template <class TQ, class TD>
template <class TB>
int QSpace<TQ,TD>::hasQOverlap(
   unsigned ica, const QSpace<TQ,TB> &B, unsigned icb,
   const char type
) const {

   unsigned ra=rank(FL), rb=B.rank(FL); int q=0;
   wbMatrix<TQ> QA,QB;
   wbindex Ia,Ib,Ja,Jb;

   if (!ica || !icb) wblog(FL,
      "ERR %s() index must be 1-based (%d,%d)",FCT,ica,icb);
   --ica; --icb;

   if (QDIM!=B.QDIM) wblog(FL,
      "ERR %s() incompatible objects\nQDIM=%d,%d",
       FCT, QDIM, B.QDIM);
   if (ica>=ra || icb>=rb) wblog(FL,
      "ERR %s() index out of bounds\n%d/%d, %d/%d (QDIM=%d,%d)",
       FCT, ica,ra, icb,rb, QDIM, B.QDIM);

     getQsub(ica,QA).makeUnique(Ja);
   B.getQsub(icb,QB).makeUnique(Jb); matchIndex(QA,QB,Ia,Ib);

   if (Ia.len==QA.dim1) { q|=2; }
   if (Ib.len==QB.dim1) { q|=4; }

   if (type=='<') { if (!(q&2)) { q=0; }} else
   if (type=='>') { if (!(q&4)) { q=0; }} else
   if (type=='=') { if (  q!=6) { q=0; }} else {
      if (type) wblog(FL,"WRN %s() invalid type=%c<%d>",FCT,type,type);
      if (!Ia.len) { q=0; }
   }

   if (q && type) {
      for (unsigned i=0; i<Ia.len; ++i) {
         if (!DATA[Ja[Ia[i]]]->hasSameSize(ica, *B.DATA[Jb[Ib[i]]],icb)) {
            q=-q; break;
         }
      }
   }

   return q;
};

template <class TQ, class TD>
mxArray* QSpace<TQ,TD>::QIDX_toMx() const {

   unsigned k, r;
   mxArray *a;

   if (QIDX.isEmpty())
   return mxCreateCellMatrix(0,0);

   if (!QDIM || QIDX.dim2%QDIM) wblog(FL, 
      "ERR %s() invalid QSpace (QIDX: %dx%d; %d)",
       FCT,QIDX.dim1,QIDX.dim2,QDIM
   );

   r=QIDX.dim2/QDIM; a=mxCreateCellMatrix(1,r);
   if (!a) wblog(FL,"ERR could not allocate cell array (%d)", r);

   for (k=0; k<r; ++k)
   mxSetCell(a, k, getQsub(k).toMx());

   return a;
}

template <class TQ, class TD>
mxArray* QSpace<TQ,TD>::INFO_toMx() const {

   if ((qtype.isEmpty() && CGR.isEmpty() &&
       itags.isEmpty() && otype==QS_NONE) 
   ){ return mxCreateCellMatrix(0,0); }

   const char *fields[] = { "qtype","otype","itags","fdir","ctime","cgr" };

   mxArray *S=mxCreateStructMatrix(1,1,6,fields);

   mxSetFieldByNumber(S,0,0, qStr()     .toMx());
   mxSetFieldByNumber(S,0,1, otype2Str().toMx());
   mxSetFieldByNumber(S,0,2, itags      .toMx());
   mxSetFieldByNumber(S,0,3, fdir       .toMx()); 
   mxSetFieldByNumber(S,0,4, numtoMx(ctime)    ); 

   mxSetFieldByNumber(S,0,5,         CGR.toMx()); 

   return S;
};

template <class TQ, class TD>
mxArray* QSpace<TQ,TD>::DATA_toMx() const {

   mxArray *a=mxCreateCellMatrix(DATA.len, DATA.len ? 1 : 0);
   if (!a) wblog(FL,"ERR failed to allocate cell array (%d)",DATA.len);
   if (DATA.len) {
      unsigned i, n=DATA.len; int np=1;
      wbvec<mxArray*> x(n);

      if (n>1 && !omp_in_parallel()) { 
         np=MAX( QSP_NUM_THREADS, OMP_NUM_THREADS );
         if (np>int(n)) { np=n; }
      }

      for (i=0; i<n; ++i) { x[i]=0; }

      Wb::LogException ex; 

     #pragma omp parallel for num_threads(np)
      for (unsigned i=0; i<n; ++i) { if (!ex) {
         try { x[i]=DATA[i]->toMx(); }
         catch (Wb::LogException &e_) { ex+=e_; }
         catch (...) { ++ex; }
      }}
      ex.report(FLF);

      for (i=0; i<n; ++i) { if (x[i]) mxSetCell(a,i,x[i]); }
   }

   return a;
};

template <class TQ, class TD>
mxArray* QSpace<TQ,TD>::DATA_save2Mx(char vflag) {

   mxArray *a=mxCreateCellMatrix(DATA.len, DATA.len ? 1 : 0);
   if (!a) wblog(FL,"ERR failed to allocate cell array (%d)",DATA.len);

   if (vflag) { size_t sd=getDataSize();
      if (sd<(1<<30)) vflag=0; else wblog(FL,
      "TST %s() %d entries @ %.3gG ...",FCT,DATA.len,sd/double(1<<30));
   }

   for (unsigned k=0; k<DATA.len; ++k) {
      mxSetCell(a, k, DATA[k]->toMx());
      if (vflag && k%10==0) PRINTF("  %6d/%ld ...  \r",k+1,DATA.len);
      DATA[k]->init();
   }
   if (vflag) PRINTF("\r%60s\r","");

   return a;
};

template <class TQ, class TD>
mxArray* QSpace<TQ,TD>::toMx(bool setcls) const { 

   mxArray *S=mxCreateStruct(1,1);
   this->add2MxStruct(S,0,setcls);
   return S;
};

template <class TQ, class TD>
mxArray* QSpace<TQ,TD>::save2Mx(char vflag) {

   mxArray *S=mxCreateStruct(1,1);
   this->save2MxStruct(S,0,0,vflag);
   return S;
};

template <class TQ, class TD>
mxArray* QSpace<TQ,TD>::mxCreateStruct(unsigned m, unsigned n) const {

   const char *fields[] = { "Q","data","info" };
   return mxCreateStructMatrix(m,n,3,fields);
}

template <class TQ, class TD>
void QSpace<TQ,TD>::add2MxStruct(
   mxArray *S, unsigned i,
   bool setcls, 
   char chk
 ) const {

   mxArray *a; 

#ifdef WBC_QSPACE_IO
   Wb::Clock clk("QS:toMx",1); 
#endif

   if (chk) { unsigned s=0; int q=0; 
   if (S==nullptr || (s=mxGetNumberOfElements(S))<1 || i>=s ||
      (q=mxGetFieldNumber(S,"qtype"))<0) wblog(FL,
      "ERR %s() must follow mxCreateStruct()\n%lx, %d/%d, %d",
       FCT,S,i+1,s,q);
   }

   a=QIDX_toMx(); mxSetFieldByNumber(S,i,0,a);
   a=DATA_toMx(); mxSetFieldByNumber(S,i,1,a); 
   a=INFO_toMx(); mxSetFieldByNumber(S,i,2,a);

   if (setcls && Wb::my_caller_tid==omp_get_thread_num()) { 
      mxSetClassName(S,"QSpace");
   }
};

template <class TQ, class TD>
void QSpace<TQ,TD>::save2MxStruct(
   mxArray *S, unsigned i, char tst, char vflag
){
   mxArray *a; 

#ifdef WBC_QSPACE_IO
   Wb::Clock clk("QS:toMx",1); 
#endif

   if (tst) {
      unsigned s=0; int q=0; 
      if (S==nullptr || (s=mxGetNumberOfElements(S))<1 || i>=s
       || (q=mxGetFieldNumber(S,"qtype"))<0) wblog(FL,
      "ERR %s() must follow mxCreateStruct()\n%lx, %d/%d, %d",
       FCT,S,i+1,s,q);
   }

   a=QIDX_toMx();         mxSetFieldByNumber(S,i,0,a);
   a=DATA_save2Mx(vflag); mxSetFieldByNumber(S,i,1,a);
   a=INFO_toMx();         mxSetFieldByNumber(S,i,2,a);

   init();
};

template <class TQ, class TD> 
void QSpace<TQ,TD>::init(
   const char *F, int L, const mxArray *S,
   char ref, 
   unsigned k,
   char skip_empty, 
   char check_cgw_norm
){
   unsigned i,j,l,m,n,r,dim1=0,dim2=0; unsigned rk=-1;
   mxArray *aq,*ad,*ai;
   wbvector< wbMatrix<TQ>  > MQ;
   wbvector< wbMatrix<TQ>* > mq;
   const char cflag = (typeid(TD)==typeid(wbcomplex) ? 'C' : 0);
   char refD=0, refC=1, isa=0, cgr=0;

   #ifdef WBC_QSPACE_IO
      Wb::Clock clk("QS:mxInit",1); 
   #endif

   clearQSpace(); 

   if (!S || mxIsEmpty(S) || mxIsEmptyQSpace(S,k)) return;

   aq=mxGetField(S,k,"Q");
   ad=mxGetField(S,k,"data");
   ai=mxGetField(S,k,"info");

   if (mxIsQSpace(F_L,S,rk,cflag,k)<=0 || !aq || !ad) wblog(F_L,
      "ERR invalid QSpace (%d,%d)",aq==0,ad==0);

   if (ref) {
      if (ref=='r') { refD        = ref; } else
      wblog(FL,"ERR %s() got invalid ref=%c<%d>",FCT,ref,ref);

   }

   m=mxGetM(aq); r=mxGetN(aq);

   if ((m!=1 && r!=1 && m && r) || mxGetNumberOfDimensions(aq)>2) {
      wblog(F,L,"ERR %s() invalid QSpace\n"
      "Q must be blocked into row cell vector (%dx%d; %d)",
      FCT,m,r,mxGetNumberOfDimensions(aq));
   }

   r*=m; MQ.init(r); 

   for (i=0; i<r; ++i) {
      MQ[i].init(mxGetCell(aq,i));

      if (i==0) {
         dim1=MQ[i].dim1; dim2=MQ[i].dim2;
      }
      else if (MQ[i].dim1!=dim1 || MQ[i].dim2!=dim2) { wblog(F,L,
        "ERR %s() invalid QSpace\nsize mismatch in Q data (%dx%d; %dx%d)",
         FCT, MQ[i].dim1, MQ[i].dim2, dim1, dim2);
      }
   }

   mq.init(MQ.len); for (i=0; i<MQ.len; ++i) mq[i] = &MQ[i];

   QIDX.CAT(2, (const wbMatrix<TQ>**) mq.data, mq.len);
   QDIM=dim2;

   if (!QIDX.isUnique()) wblog(FL, 
      "ERR %s() invalid input QSpace\nrows in [Q{:}] are not unique",FCT);
   setupDATA();

   if (ai) {
      unsigned n=mxGetNumberOfElements(ai);
      if (n) {
         unsigned dq=0; 
         mxArray *a; QVec qv;

         if (!mxIsStruct(ai)) wblog(F,L,
            "ERR invalid QSpace.info (structure required)");
         if (n!=1) wblog(F,L,
            "ERR invalid QSpace.info (got %d elements)",n);

         a=mxGetField(ai,0,"qtype");
         if (a) { qtype.init(F,L,a);
          # ifdef QS_USING_OMP
            qtype.checkInit();
          # endif
         }
         else { qtype.init(); }

         isa=qtype.allAbelian();
         if (isa==1 && qtype.len) { qtype.ReduceU1(); }

         initOType(F_L,mxGetField(ai,0,"otype"));

         itags.init(F_L,mxGetField(ai,0,"itags"));
         if (itags.len && itags.len!=r) { wblog(FL,
            "ERR %s() invalid number of itags (%s; %d/%d)",
             FCT,IT2STR__,itags.len,r);
         }

         fdir.init(F_L,mxGetField(ai,0,"fdir"),r); 
         if (!Wb::envFERM) { fdir.init(); } 

         a=mxGetField(ai,0,"ctime");
         if (a && !mxIsEmpty(a)) { mxGetNumber(a,ctime); }

         a=mxGetField(ai,0,"cgr");
         if (a && !mxIsEmpty(a)) { cgr|=1;
            QSet<TQ> Q;
            char xflag=(isa ? 0 : 1); 
            wbvector<unsigned> qdc; qtype.Qpos(qdc);

            setupCGR(); 
            m=mxGetM(a); n=mxGetN(a);

            if (m==1 && !QIDX.dim1) { 
               if (!mxIsScalarQSpace(S) || m>1 || r) wblog(FL,
                  "ERR %s() invalid scalar QSpace",FCT);
               QIDX.init(1,0); QDIM=qtype.Qlen(); r=2; 
               setupDATA(); setupCGR(); 
            }
            else
            if (!QDIM || QIDX.dim2%QDIM || QIDX.dim2/QDIM!=r) wblog(FL,
               "ERR %s() empty QDIM while info.cgr is specified (%d/%d)",
                FCT,m,n);
            else
            if (mxGetNumberOfDimensions(ad)>2 ||
                m!=QIDX.dim1 || m!=CGR.dim1 ||
                n!=qtype.len || n!=CGR.dim2 || qtype.Qlen()!=QDIM)
            wblog(FL,"ERR %s() invalid dimensions for cell array "
               "info.cgr\n{%s} cgr: %dx%d / %dx%d; QIDX: %dx%d/%d",
                FCT, qStr().data, m,n, CGR.dim1, CGR.dim2,
                QIDX.dim1, QIDX.dim2, QDIM
            );

            if (!mxIsCell(a)&& !mxIsStruct(a)) wblog(FL,
               "ERR %s() cell array expected for info.cgr field (%s)",
                FCT,mxGetClassName(a));

            for (l=j=0; j<n; ++j) 
            for (  i=0; i<m; ++i, ++l) {
               Q.init(qtype[j], QIDX.rec(i)+qdc[j], r, QDIM, itags);
               CGR(i,j).init(F_L,a,l,refC,&Q,xflag); 
            }
         }

         if (qtype.len) { dq=qtype.Qlen();
            if (dq!=QDIM && (QIDX.dim2 || CGR.dim2)) wblog(FL,"ERR "
               "incompatible number of symmetry labels (%d/%d; %s; %s)",
               QDIM, dq, QIDX.sizeStr().data, CGR.sizeStr().data
            );
         }
      }

      fdir.check(FL, qtype, qtype.len ? -1 : QDIM);

      if (isa) {
         if (!CGR.isEmpty()) {
            for (i=0; i<CGR.dim1; ++i) 
            for (j=0; j<CGR.dim2; ++j) { if (CGR(i,j).cgb) {
               const cdata__ &c=(CGR(i,j).cgb->cgd);
               if (!c.isScalar()) wblog(FL,
                  "ERR %s() invalid abelian CGC space (%s)",FCT,SSTR(c));
               if (c.D.data && c.D.data[0]!=1) wblog(FL,
                  "ERR %s() invalid abelian CGC coefficient (%.4g)",
                   FCT,double(c.D.data[0])
               );
            }}
            checkQ_CGR(FL); CGR.init();
         }

      }
   }
   else qtype.init();

   m=mxGetM(ad); n=mxGetN(ad);

   if ((m>1 && n>1) || mxGetNumberOfDimensions(ad)>2) wblog(F,L,
      "ERR QSpace() cell VECTOR expected for data (%d,%d)",m,n);

   m*=n;

   if (DATA.len!=m) {
      if (m==1 && DATA.len==0 && QIDX.isEmpty() && CGR.isEmpty()) {
         QIDX.init(1,0); setupDATA();
         if (r<2) r=2; 
      }
      else { wblog(F,L,
        "ERR QSpace() dimension mismatch (QIDX: %dx%d, data: %d)",
         QIDX.dim1, QIDX.dim2, m);
      }
   }

   for (n=i=0; i<m; ++i) {
      DATA[i]->init(mxGetCell(ad,i), refD);
      DATA[i]->appendSingletons(FL,r,1); 
      if (!*DATA[i]) { ++n; }
   }

   if (n && skip_empty) {
      m=SkipEmptyData(FL); 
   }

   if (cgr && check_cgw_norm) { 
      int q=isNormCGW();
      if (ref) {
         if (q<=0) { wblog(PFL,
            "ERR %s() invalid CGW normalization of input QSpace ref",FCT);
         }
      }
      else {
         if (q<=0) wblog(PFL,
            "WRN %s() fixing CGW normalization of input QSpace",FCT);
         NormCGW();
      }
   }

};

template <class TQ, class TD> inline
void QSpace<TQ,TD>::init2ref(const QSpace<TQ,TD> &A) {

   clearQSpace(); 

   qtype.init2ref(A.qtype); otype=A.otype; QDIM=A.QDIM;
   itags.init2ref(A.itags);
   fdir .init2ref(A.fdir );

   QIDX .init2ref(A.QIDX );
   DATA .init2ref(A.DATA ); if (A.CGR) {
   CGR  .init2ref(A.CGR  ); }

   isref=1;
};

template <class TQ, class TD> inline
void QSpace<TQ,TD>::setRand(double nrm) {

   unsigned i;
   for (i=0; i<DATA.len; ++i) DATA[i]->setRand(); 

   if (nrm>0) {
      double fac=nrm/sqrt(double(norm2()));
      for (i=0; i<DATA.len; ++i) (*DATA[i]) *= fac;
   }
}

template <class TQ, class TD>
void QSpace<TQ,TD>::info(const char *vname, 
   char nl,   
   char nind  
 ) const {

   char cgflag=gotCGS(); 
   unsigned r=(QDIM? QIDX.dim2/QDIM : 0);

   unsigned lsz=2*r+2; if (lsz<8) { lsz=8; }

   wbvec<char> sout(256), ind; 
   wbvector<widx_t> D,DX;

   if (nind) {
      if (nind>0)
           { ind.init(nind+1); memset(ind.data,' ',nind); ind[nind]=0; }
      else { nind=0; }
   }

   if (QDIM && QIDX.dim2 && QIDX.dim2%QDIM) wblog(FL,
      "ERR %s() %d/%d = ?",FCT, QIDX.dim2, QDIM);
   getDim(D,&DX);

   sout.catf(0,0,"QSpace<%s,%s> %s",
      sTSTR(TQ), sTSTR(TD), qtype.len ? qStr().data : "");

   if (itags) { sout.catf(0,0," %s",STR(itags)); }
   if (fdir ) { sout.catf(0,1," %s",STR(fdir )); }

   if (nl&1) { PRINTF("\n"); }
   if (vname && *vname) {
      if (nind)
           { PRINTF("%s%s",ind.data,vname); }
      else { PRINTF("%-4s",         vname); }
   }
   else if (nind) { PRINTF("%s",ind.data); }

   PRINTF(" %-34s %2ldx {%d x%d}",sout.data,QIDX.dim1,r,QDIM);
   PRINTF("  %-*s",lsz,SSTR(D )); if (cgflag>0) {
   PRINTF("  %-*s",lsz,SSTR(DX)); }
   PRINTF("%s\n", isref? "  *REF*":"");

   if (nl&2) { PRINTF("\n"); }
};

template <class TQ, class TD> 
void QSpace<TQ,TD>::print(
   const char *istr, 
   char vflag) const {

   unsigned i=0,j,n, N=QIDX.dim1, r=rank(FL), m=(N<12? N:4); int q;
   const TQ *qs=QIDX.data;
   char cgflag=gotCGS(FL);

   wbvec<char> qstr(256);

   unsigned lsz=2*r+2;   
   if (lsz<8) { lsz=8; } 

   if (vflag & 'D') { vflag |=3; } 
   info(istr, vflag&3, 0); 

   QVec qt(qtype);
   if (!qt) { qt.init2val(QDIM,QType(QT_U1)); } 

   for (; i<N; ++i) { qstr.reset();
      if (i>=m) {
         if (N>2*m) { i=N-m; qs=QIDX.rec(i);
            PRINTF("    :   ...\n");
         }; m=N; 
      }

      if (i<N) {
         for (j=0; j<r; ++j, qs+=QDIM) {
            if (j) { qstr.catf(0,0," ; "); }
            if (qt.print_qset(FL,qs,qstr)<0) { break; }
         }
      } else { qstr.catf(0,0," !?"); }

      PRINTF("%5d. { %s }  %-*s ", i+1, qstr.data, lsz,
         i<DATA.len ? DATA[i]->sizeStr(r,"_").data : "!?");

      if (i>=DATA.len) { PRINTF("   !?\n"); continue; }
      if (vflag & 64) { 
         if (!(DATA[i]->printdata(istr))) { PRINTF("\n"); }
      }
      else {
         if (cgflag>0) { 
            for (j=0; j<CGR.dim2; ++j) { if (!CGR(i,j).isAbelian()) {
               PRINTF(" %-*s",lsz,SSTR(CGR(i,j))); }
            }
         }

         n=DATA[i]->numel();

         if (n==0) { PRINTF("   []\n"); } else
         if (n==1) { PRINTF(" %9s", Wb::num2Str(DATA[i]->data[0]).data); }
         else {
            n *= sizeof(TD);
            if (n<(1<<10)) { PRINTF("   %5d B", n); } else
            if (n<(1<<20)) { PRINTF("   %5.2f kB", n/double(1<<10)); }
            else           { PRINTF("   %5.2f MB", n/double(1<<20)); }
         }

         if ((q=DATA[i]->isRef())) {
            if (q<0)
                 { PRINTF("   *ref\n"); }
            else { PRINTF("   (%d ref%s)\n",q,q==1?"":"s"); }
         }
         else PRINTF("\n");
      }
   }

   if (vflag&2) { PRINTF("\n"); } 
};

#endif

