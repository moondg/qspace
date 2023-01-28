#ifndef __WB_ORTHO_HH__
#define __WB_ORTHO_HH__

#define WB_STOL_SVD 1E-8

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
// change log
//  * split off mpsortho.hh from mpsortho.cc
//  * transferred former QSpace member functions
//    toBlockMatrix() and initFromBlockMatrix()
//    to class QBlock(), as it fits more naturally there.
// Wb,Aug05,15

template <class TQ, class TD>
mxArray* orthoQS(
   const QSpace<TQ,TD> &PSI, 
   QSpace<TQ,TD> &A1,
   QSpace<TQ,TD> &A2, unsigned K,
   unsigned Nkmin, unsigned Nkeep,
   double stol,  
   itag_ tx,      
   const char *info
);

template <class TQ, class TD>
mxArray* getSVD(
   const QSpace<TQ,TD> &PSI, 
   QSpace<TQ,TD> &U,
   QSpace<TQ,TD> &S,
   QSpace<TQ,TD> &VT, unsigned K,
   unsigned Nkmin, unsigned Nkeep, double stol,
   itag_ tx,      
   const char *info
);

char mxGetLRDir(
   const char *F, int L, const mxArray *argin, unsigned k
);

template <class TQ, class T1, class T2, class T3>
void twoSiteInit(
   const QSpace<TQ,T1> &Psi1,
   const QSpace<TQ,T2> &Psi2,
   QSpace<TQ,T3> &PSI, unsigned &K,
   widx_t ic1, widx_t ic2,
   char twoSite, char ldir
);

template <class TQ, class T1, class T2, class T3>
void twoSiteFinal(
   const QSpace<TQ,T1> &Psi1,
   const QSpace<TQ,T2> &Psi2,
   QSpace<TQ,T3> &A1,
   QSpace<TQ,T3> &A2,
   unsigned ic1, unsigned ic2, char twoSite, char ldir
);

template <class TQ, class TD, class TE=double>
class QBlock { 
  public:

    QBlock() : qdim_tot(1) {};

    QBlock& init() { qdim_tot=1; qdim.init();
       I0.init(); Ik.init(); It.init(); idx.init();
       Q1.init(); S1.init(); D1.init(); A.init(); S.init();
       Q2.init(); S2.init(); D2.init(); U.init(); Vc.init();
       return *this;
    };

    QBlock& operator=(const QBlock &B) {
       if (this!=&B) { qdim_tot=B.qdim_tot; qdim=B.qdim;
          I0=B.I0; Ik=B.Ik; idx=B.idx;
          Q1=B.Q1; S1=B.S1; D1=B.D1; A=B.A; S=B.S;
          Q2=B.Q2; S2=B.S2; D2=B.D2; U=B.U; Vc=B.Vc;
       }
       return *this;
    };

    QBlock& init_bare_refA(const QBlock &B, unsigned w=0);

    void toBlockMatrix(
       wbarray<TD> &MM, 
       unsigned K,      
       const wbperm &P=wbperm()
    );

    void initFromBlockMatrix(const char *F, int L,
       QSpace<TQ,TD> &B,      
       const wbarray<TD> &MM, 
       const iTags &it,       
       const wbperm *cgp=NULL 
    ) const;

    void updateBlockDim(unsigned dim2);

    double get_s2t(const wbvector<char> &mark) {
       double s2t=0; unsigned i=0;
       if (mark.len!=S.len) wblog(FL,"ERR %s() "
          "length mismatch (mark.len=%d/%d)",FCT,mark.len,S.len);
       mark.find(Ik);
          for (; i<mark.len; ++i) { if (!mark[i]) s2t+=S[i]*S[i]; }
          if (qdim_tot!=1) { s2t*=qdim_tot; }
       return s2t;
    };

    void print(
       const char *F=NULL, int L=0, const iTags *t=NULL
     ) const {

       wblog(F_L,"%s() QBlock content ======================= %N",FCT);
       A.print("sub(PSI)"); if (t) {
       printf("=> itags: \"%s\"",STR_(t)); }
       printf("      qdim_tot=%d\n\n",qdim_tot);
       printf("      Q =[ %s; %s ]\n",STR(Q1),STR(Q2));
       printf("      S =[ %s; %s ]\n",STR(S1),STR(S2));
       printf("      D1=[ %s ]\n",STR(D1));
       printf("      D2=[ %s ]\n\n",STR(D2));
    };

    mxArray* toMx() const { 
       const char *f[]={ "A","qdim_tot","qdim",  
          "Q1","S1","D1", "Q2","S2","D2",        
          "I0","Ik","idx","U", "S", "Vc"         
       };
       unsigned i=0, n=1;
       mxArray *q=mxCreateStructMatrix(n,1,15,f);

       mxSetFieldByNumber(q,i, 0, A .toMx());
       mxSetFieldByNumber(q,i, 1, numtoMx(double(qdim_tot)));
       mxSetFieldByNumber(q,i, 2, qdim.toMx());
       mxSetFieldByNumber(q,i, 3, Q1.toMx());
       mxSetFieldByNumber(q,i, 4, S1.toMx());
       mxSetFieldByNumber(q,i, 5, D1.toMx());
       mxSetFieldByNumber(q,i, 6, Q2.toMx());
       mxSetFieldByNumber(q,i, 7, S2.toMx());
       mxSetFieldByNumber(q,i, 8, D2.toMx());
       mxSetFieldByNumber(q,i, 9, I0.toMx());
       mxSetFieldByNumber(q,i,10, Ik.toMx());
       mxSetFieldByNumber(q,i,11, idx.toMx());
       mxSetFieldByNumber(q,i,12, U .toMx());
       mxSetFieldByNumber(q,i,13, S .toMx());
       mxSetFieldByNumber(q,i,14, Vc.toMx());

       return q;
    };

    void put(const char *vname, const char *ws="caller") {
       mxArray *a=toMx();
       mxPutAndDestroy(FL,a,vname,ws);
    };

    unsigned qdim_tot;       
    wbvector<unsigned> qdim; 

    QSpace<TQ,TD> A;  
    wbindex I0;       
    wbindex Ik;       
    wbindex It;       

    wbMatrix<TQ>     Q1, Q2;  
    wbMatrix<widx_t> S1, S2;  
    wbMatrix<widx_t> D1, D2;  

    wbMatrix<int> idx; 

    wbarray<TD> U,Vc; 
    wbvector<TE> S;   

  protected:
  private:
};

template <class TQ, class TD>
class SVD_Data { 
  public:
    SVD_Data() {};

    void blockSVD(const QSpace<TQ,TD> &PSI, const unsigned K);

    unsigned dmrgTruncate(
       unsigned Nkmin, unsigned Nkeep, double stol,
       itag_ tx,            
       wbMatrix<double> &SV, 
       QSpace<TQ,TD> &UQ,  
       QSpace<TQ,TD> *SQ,  
       QSpace<TQ,TD> &VC,  
       mxArray *Sout=NULL, char vflag=0              
    );

    mxArray* toMx() const;

    void put(const char *vname, const char *ws="base") const {
       put(0,0,vname,ws);
    };

    void put(const char *F, int L,
       const char *vname, const char *ws="base"
    ) const {
       mxArray *a=toMx();
       int i=mexPutVariable(ws,vname,a);
       if (i) wblog(F,L,"ERR failed to write variable `%s' (%d)",vname,i);
       if (F) wblog(F,L,"I/O putting variable `%s' to `%s'",vname,ws);
       mxDestroyArray(a);
    };

    wbvector< QBlock<TQ,TD,double> > QB;

  protected:
  private:
};

#endif

