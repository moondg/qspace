#ifndef __WB_CPAT_HCC__
#define __WB_CPAT_HCC__

/* ------------------------------------------------------------------ *
 * class CPAT: class to store contraction pattern in QSpace structure
 * ------------------------------------------------------------------ *
 * contraction pattern (default: twoSite mode)
 *
 * deals with MatLab cell array of the form { A, i1, i2 [,opA] }
 * where A can be a single QSpace or a cell vector of these
 * i[12] specifies index onto which to contract onto Psi[12]
 * NB! either i1 or i2 is set, the other one must be empty, e.g. []
 *
 * single site mode: expects { A, ic, [,opA] } instead as
 * contraction pattern where A is always contracted onto "Psi1"
 *
 * ------------------------------------------------------------------ *
 * AW (C) Jun 2006
 * ------------------------------------------------------------------ */

#define CP1 "{ A, ic, [,opA] }"
#define CP2 "{ A, i1, i2 [,opA] }"
#define _CP (twoSite ? "{ A, i1, i2 [,opA] }" : "{ A, ic, [,opA] }")

template <class TQ, class TD>
class CPAT_Set;

int isCPAT(const mxArray *C, char twoSite=-1);
int isCPattern(const mxArray *C);

unsigned NR_ICTR_CALL=0;

template <class TQ, class TD>
class CPAT {

  public:

    CPAT() : z(0), sameCtrIdx(0), nr_check_ctr(0), HPsi_counter(0) {};

    CPAT(const CPAT<TQ,TD> &D) { 
       wblog(FL,"ERR do not copy CPAT object!");
    };

    void init(
      TD oz, const mxArray**, unsigned n,
      unsigned ic1, unsigned r1, unsigned ic2, unsigned r2,
      const char twoSite, const char ldir
    );

    int  getCtrIndex(const QSpace<TQ,TD> &PSI);

    QSpace<TQ,TD>& HZTimes(
       const QSpace<TQ,TD> &PSI, QSpace<TQ,TD> &HZPsi,
       const char ='N', char disp=1
    );

    QSpace<TQ,TD>& HZApplyDiag(
       const QSpace<TQ,TD> &PSI, QSpace<TQ,TD> &EPsi,
       TD a, TD b, char Iflag=0, char disp=1
    );

    void check_contraction(const QSpace<TQ,TD> &PSI, char disp=1);
    void initCounter() {
       HPsi_counter=0;
    };

    QSpace<TQ,TD>& HZTimes_tst( 
    const QSpace<TQ,TD> &PSI, QSpace<TQ,TD> &HZPsi,
    const char ='N', char disp=1);

    void getQRange(
       unsigned k, wbMatrix<TQ> &QQ, wbvector<unsigned> &SS,
       QSpace<TQ,TD> *B=NULL
    ) const;

    mxArray* toMx() const;

    void put(const char *vname, const char* ws="base") const {
       mxPutAndDestroy(FL,toMx(),vname, ws);
    };

    void PermIC(const wbperm &P);

    wbvector< wbvector< CPAT_Set<TQ,TD> > > CP;
    TD z;

    unsigned sameCtrIdx;
    unsigned nr_check_ctr; 
    unsigned HPsi_counter;

  protected:
  private:
};

template <class TQ, class TD>
class CPAT_Set {

  public:

    CPAT_Set () : opA(0), ic(0) {};

    CPAT_Set (const CPAT_Set<TQ,TD> &D) : opA(0), ic(0) {
       wblog(FL,"WRN Do not copy CPAT_Set!");
    };

    int init(
      const mxArray *a,
      unsigned ic1, unsigned r1, unsigned ic2, unsigned r2,
      const char twoSite, const char ldir
    );

    int getCtrIndex(const QSpace<TQ,TD> &PSI);
    int getCtrIndex(const wbvector< QSpace<TQ,TD> >&, const char ='N');

    int getCtrIndex(const unsigned i, const QSpace<TQ,TD>&, char ='N');

    void contract(
       const QSpace<TQ,TD> &PSI, wbvector< QSpace<TQ,TD> > &C,
       char='N') const;

    void contract(
       const wbvector< QSpace<TQ,TD> > &PSI1, wbvector< QSpace<TQ,TD> > &C,
       char='N') const;

    void contract(
       const unsigned k, const QSpace<TQ,TD> &PSI, QSpace<TQ,TD> &C,
       char='N') const;

    void contract_diag(
       const QSpace<TQ,TD> &PSI, wbvector< QSpace<TQ,TD> > &C,
       TD a, TD b, char Iflag) const;

    void contract_diag(
       const wbvector< QSpace<TQ,TD> > &PSI1, wbvector< QSpace<TQ,TD> > &C,
       TD a, TD b, char Iflag) const;

    void contract_diag(
       const unsigned k, const QSpace<TQ,TD> &PSI, QSpace<TQ,TD> &C,
       TD a, TD b, char Iflag) const;

    mxArray* toMx() const;

    wbvector< QSpace<TQ,TD> > A;
    char opA;    
    unsigned ic; 

    wbperm P;

    wbvector< QSpace<TQ,TD> >      CQn, CQt; 
    wbvector< wbMatrix<unsigned> > CIn, CIt; 

  protected:
  private:
};

template <class TQ, class TD>
void CPAT<TQ,TD>::init(
   TD oz, const mxArray** aa, unsigned N,
   unsigned ic1, unsigned r1, unsigned ic2, unsigned r2,
   const char twoSite, const char ldir
){
   unsigned i,j,e,k,p,m=0,n;
   wbvector<unsigned> dd;
   mxArray *a;

   if (ic1>=r1 || ic2>=r2) wblog(FL,
    "ERR Contraction index out of bounds (%d/%d; %d/%d (%d))",
     ic1+1, r1, ic2+1, r2, twoSite
   );

   z=oz;

   CP.initDef(N);

   for (p=k=0; k<N; k++) {
      if (!mxIsCell(aa[k]))
      wblog(FL, "ERR Invalid (pair) pattern %s", _CP);

      n=isCPAT(aa[k], twoSite);

      if (!n) wbdie(FL,str);

      CP[p].initDef(n);

      if (n==1) { 
         e=CP[p][0].init(aa[k],ic1,r1,ic2,r2,twoSite,ldir);
         if (e) wbdie(FL,str);
         m=CP[p][0].A.len;
      }
      else {
         for (i=0; i<n; i++) {
            a=mxGetCell(aa[k],i);
            e=CP[p][i].init(a,ic1,r1,ic2,r2,twoSite,ldir);

            if (e) {
               printf("\n%s:%d contraction pattern %d:%d\n", FL, k+1, i+1);
               printf("%s\n\n",str); ExitMsg("");
            }

            if (i>0) { if (m!=CP[p][i].A.len) wblog(FL,
              "ERR Pattern pair dimension mismatch in A (%d.%d: %d,%d)",
               k+1, i, CP[p][0].A.len, CP[p][i].A.len); }
            else m=CP[p][0].A.len;
         }
      }

      dd.init(m);
      for (i=0; i<n; i++) {
         const wbvector< QSpace<TQ,TD> > &Aki=CP[p][i].A;
         for (j=0; j<m; j++) if (!Aki[j].isEmpty()) dd[j]++;
      }

      for (j=0; j<m; j++) if (dd[j]==n) { p++; break; }
   }

   if (p<N) {
      if (p==0) wblog(FL,"WRN eigs: no relevant contraction pattern!");

      CP.len=p;
   }

   NR_ICTR_CALL=0;
}

template <class TQ, class TD>
int CPAT_Set<TQ,TD>::init(
   const mxArray *C,
   unsigned ic1, unsigned r1, unsigned ic2, unsigned r2,
   const char twoSite, const char ldir
){

   unsigned k,m,w, ra, ica, ic0, N=mxGetNumberOfElements(C);
   const char isC=(typeid(TD)==typeid(wbcomplex) ? 'C' : 0);
   int r=-1;

   mxArray *a; str[0]=0;

   if (ic1>=r1 || ic2>=r2) {
      sprintf(str, "%s:%d ERR Invalid input (%d/%d; %d/%d (%d))",
      FL, ic1+1, r1, ic2+1, r2, twoSite);
      return 1;
   }

   m=isCPattern(C);

   if ((char)m!=twoSite+1) wblog(FL,
   "ERR invalid input (not a %d-site CPattern; %d,%d)", twoSite+1,m);

   a=mxGetCell(C,0);

   if (mxIsQSpace(0,0,a,r,isC)>0) { 
      A.initDef(1);
      A[0].init(FL,a,'r');

      if (!A[0].isConsistent_r(2)) {
         sprintf(str,"%s:%d rank-2 object required (arg #1; %s)",FL,_CP);
         return 1;
      }
   }
   else { 
      if (mxIsQSpaceVec(a,2,isC)) { sprintf(str,
         "%s:%d rank-2 QSpace (vector) required (arg #1; %s)",FL,_CP);
          return 1;
      }
      mxInitQSpaceVec(FL,a,A,'r');
   }

   if (twoSite) {
      w = (mxIsEmpty(mxGetCell(C,1)) ? 0 : 1) +
          (mxIsEmpty(mxGetCell(C,2)) ? 0 : 2);
      k=w;

      if (w<1 || w>2) {
         sprintf(str,"%s:%d Either arg #2 or #3 in CPAT %s",FL,_CP);
         return 1;
      }
   }
   else {
      k=1; w=(ldir>0 ? 1 : 2);

      if (mxIsEmpty(mxGetCell(C,k))) {
         sprintf(str,"%s:%d Invalid arg #%d in CPAT %s",FL,k,CP1);
         return 1;
      }
   }

   ra  = (w==1 ?  r1 :  r2);
   ica = (w==1 ? ic1 : ic2);

   a=mxGetCell(C,k);

   if (mxGetNumber(a,ic)) { sprintf(str+strlen(str),
      "\n%s:%d Invalid arg #%d in CPAT %s", FL, k+1, _CP);
       return 1;
   }

   if (!ic || ic>ra) { sprintf(str+strlen(str),
      "\n%s contraction index out of bounds (%d/%d; %d,%d) in CPAT %s",
       SHORT_FL, ic, ra,r1,r2,_CP); return 1;
   }

   ic0=ic; 
   ic--;   

   if (twoSite && ic==ica) { sprintf(str+strlen(str),
      "\n%s Index %d overlaps with ctr index in CPAT %s",
       SHORT_FL, ic0, _CP); return 1;
   }

   if (twoSite) {
      if (ic>ica) ic--; 
      if (ldir>0 && w==2) ic+=(r1-1); else
      if (ldir<0 && w==1) ic+=(r2-1);
   }
   else {
      if (ic>ica) ic--; else
      if (ic==ica) ic=ra-1; 
   }

   opA=0; k=(twoSite ? 3 : 2);

   if (k<N) {
      a=mxGetCell(C,k);

      if (mxGetNumber(a,opA)) { sprintf(str,
         "%s:%d Invalid arg #%d in CPAT %s", FL, 3+twoSite, _CP);
          return 1;
      }
      if (N>k+1) { sprintf(str, 
         "%s:%d Too many arguments for CPAT %s", FL, _CP);
          return 1;
      }
   }

   opA=toupper(opA); if (!opA || opA==' ') opA='N';

   if (!Wb::strchri("NTC*",opA)) { sprintf(str,
     "%s:%d ERR invalid opA=%c<%d> in CPAT %s", FL, opA, opA, _CP);
      return 1;
   }

   if (opA!='N') {
      char tflag= (opA=='T' || opA=='C'),
           cflag= (opA=='C' || opA=='*');
      wbperm P("2 1");

      for (unsigned i=0; i<A.len; i++) {
         if (A[i].isEmpty()) continue; 
         if (tflag) A[i].Permute(P); 
         if (cflag) A[i].Conj();
      }
   }

   return 0;
}

template <class TQ, class TD>
void CPAT<TQ,TD>::PermIC(const wbperm &P) {
   unsigned i,j, m=CP.len, n;
   wbperm iP;

   P.invert(iP);

   for (i=0; i<m; i++) {
      for (n=CP[i].len, j=0; j<n; j++)
      CP[i][j].ic = iP[CP[i][j].ic];
   }
}

template <class TQ, class TD>
int CPAT<TQ,TD>::getCtrIndex(const QSpace<TQ,TD> &PSI) {

   unsigned i,j,n; int m=0;
   const QSpace<TQ,TD> PSI2;

   if ((++NR_ICTR_CALL)%20==0) wblog(FL,
   "%NWRN calling %s already %d times (?)", FCT, NR_ICTR_CALL);

   for (i=0; i<CP.len; i++) {
      m += CP[i][0].getCtrIndex(PSI); 

      for (n=CP[i].len, j=1; j<n; j++) {
         if (CP[i][0].A.len!=CP[i][j].A.len) wblog(FL,
            "ERR Severe dimension mismatch (%d: %d,%d)", j,
             CP[i][0].A.len, CP[i][j].A.len);

         m += CP[i][j].getCtrIndex(CP[i][j-1].CQn,'N');
         m += CP[i][j].getCtrIndex(CP[i][j-1].CQt,'T');
      }
   }

   return m;
}

template <class TQ, class TD>
int CPAT_Set<TQ,TD>::getCtrIndex(const QSpace<TQ,TD> &PSI) {

   unsigned i,m=0,n=A.len;

   if (CIn.len!=n) { CIn.initDef(n); m++; }
   if (CIt.len!=n) { CIt.initDef(n); m++; }
   if (CQn.len!=n) { CQn.initDef(n); m++; }
   if (CQt.len!=n) { CQt.initDef(n); m++; }

   if (A.len==0) return (m!=0);

   P.Index(PSI.rank());
   P.mvFirstTo(ic); 

   for (i=0; i<A.len; i++) {
      m += getCtrIndex(i,PSI,'N'); 
      m += getCtrIndex(i,PSI,'T'); 
   }

   return (m!=0);
}

template <class TQ, class TD>
int CPAT_Set<TQ,TD>::getCtrIndex(
   const wbvector< QSpace<TQ,TD> > &PSI1, 
   const char tflag
){
   unsigned i, m=0, n=A.len;

   if (!Wb::strchri("NT",tflag)) wblog(FL,
   "ERR Invalid flag %c<%d>",tflag,tflag);

   if (tflag=='N') {
      if (CIn.len!=n) { CIn.initDef(n); m++; }
      if (CQn.len!=n) { CQn.initDef(n); m++; }
   }
   else {
      if (CIt.len!=n) { CIt.initDef(n); m++; }
      if (CQt.len!=n) { CQt.initDef(n); m++; }
   }

   if (n==0) return (m!=0);

   if (n!=PSI1.len) wblog(FL,
   "ERR Severe dimension mismatch (%d,%d)", n, PSI1.len);

   for (i=1; i<n; i++) if (PSI1[i].QDIM!=PSI1[0].QDIM)
   wblog(FL, "ERR Severe QDIM inconsistency (%d,%d).",
   PSI1[i].QDIM, PSI1[0].QDIM);

   P.Index(PSI1[0].rank());
   P.mvFirstTo(ic); 

   for (i=0; i<n; i++)
   m += getCtrIndex(i, PSI1[i], tflag);

   return (m!=0);
}

template <class TQ, class TD>
int CPAT_Set<TQ,TD>::getCtrIndex(
   const unsigned i,
   const QSpace<TQ,TD> &PSI,
   char tflag
){
   int e;
   unsigned ica = tflag=='N' ? 2 : (tflag=='T' ? 1 : 0); 

   QSpace<TQ,TD>      &CQi = (tflag=='N' ? CQn[i] : CQt[i]);
   wbMatrix<unsigned> &CIi = (tflag=='N' ? CIn[i] : CIt[i]);

   QSpace<TQ,TD> CQ0;
   const wbMatrix<unsigned> CI0(CIi);

   if (!ica) wblog(FL, "ERR Invalid flag %c<%d>",tflag,tflag);

   CQ0.QDIM=CQi.QDIM; CQi.QDIM=PSI.QDIM;
   CQ0.QIDX=CQi.QIDX;

   e=A[i].contract_getIdxSet(ica, PSI, ic+1, CQi.QIDX, CIi);
   if (e<0) wbdie(FL,str); 

   CQi.PermuteQ(P);

   return (CQi.QDIM!=CQ0.QDIM || CQi.QIDX!=CQ0.QIDX || CIi!=CI0);
}

template <class TQ, class TD>
QSpace<TQ,TD>& CPAT<TQ,TD>::HZTimes(
   const QSpace<TQ,TD> &PSI,
   QSpace<TQ,TD> &HZPsi,
   const char tflag,
   char disp
){
   unsigned i,j,n;

#ifdef WB_CLOCK
   Wb::Clock clk("cpat:HZTimes",0); 
#endif

   wbvector< QSpace<TQ,TD> > X1,X2;

   if (!Wb::strchri("NTC",tflag))
   wblog(FL,"ERR Invalid flag %c<%d>", tflag, tflag);

   check_contraction(PSI,disp);
   HPsi_counter++;

   HZPsi=PSI;
   HZPsi *= (tflag=='C' ? -CONJ(z) : -z);

   for (i=0; i<CP.len; i++) {

      CP[i][0].contract(PSI, X1, tflag);

      for (n=CP[i].len, j=1; j<n; j++) {
         X1.save2(X2);
         CP[i][j].contract(X2,X1,tflag);
      }

      for (j=0; j<X1.len; j++)
      X1[j].Append2AndDestroy(FL,HZPsi);
   }

   HZPsi.makeUnique();

   return HZPsi;
};

template <class TQ, class TD>
QSpace<TQ,TD>& CPAT<TQ,TD>::HZApplyDiag(
   const QSpace<TQ,TD> &PSI, QSpace<TQ,TD> &EPsi,
   TD a, TD b, char Iflag, char disp
){
   unsigned i,j,n;

   wbvector< QSpace<TQ,TD> > X1,X2;

   check_contraction(PSI,disp);

   EPsi.init();

   for (i=0; i<CP.len; i++) {

      CP[i][0].contract_diag(PSI, X1, a,b,Iflag);

      for (n=CP[i].len, j=1; j<n; j++) {
         X1.save2(X2);
         CP[i][j].contract_diag(X2,X1, a,b,Iflag);
      }

      for (j=0; j<X1.len; j++) {
         X1[j].SkipEmptyData();
         X1[j].Append2AndDestroy(FL,EPsi);
      }
   }

   EPsi.makeUnique();

   if (PSI.QIDX!=EPsi.QIDX) { 
      wblog(FL,
        "WRN QIDX different for PSI and EPsi !? (%d/%d)",
         EPsi.QIDX.dim1,PSI.QIDX.dim1);
      wblog(FL,"TST a=%g, b=%g, '%c'",a,b,Iflag);
#ifdef MATLAB_MEX_FILE
      PSI.put("PSI"); EPsi.put("EPsi");
#endif
   }

   return EPsi;
};

template <class TQ, class TD>
void CPAT<TQ,TD>::check_contraction(
   const QSpace<TQ,TD> &PSI, char disp
){
   if (!nr_check_ctr) getCtrIndex(PSI); else
   if (sameCtrIdx<nr_check_ctr) {
      unsigned i=getCtrIndex(PSI);
      if (i==0) {
         sameCtrIdx++; 

         if (sameCtrIdx==nr_check_ctr && disp>10)
         wblog(FL,"<i> using sameCtrIdx=%d", sameCtrIdx);
      }
      else {
         if (sameCtrIdx && disp)
         wblog(FL,"WRN sameCtrIdx=%d reset to 0.", sameCtrIdx);

         sameCtrIdx=0;
      }
   }
   else ARG_CHECK=0;
};

template <class TQ, class TD>
inline void CPAT<TQ,TD>::getQRange(
   unsigned k,
   wbMatrix<TQ> &QQ,
   wbvector<unsigned> &SS,
   QSpace<TQ,TD> *B
) const {

   unsigned i,j,l,n,nA,m=CP.len; int e=0;
   wbvector<unsigned> S,D;
   wbMatrix<TQ> Q;
   wbperm P;

   for (i=0; i<m; i++) { n=CP[i].len;
   for (j=0; j<n; j++) { nA=CP[i][j].A.len; ;
       if (CP[i][j].ic!=k) continue;

       for (l=0; l<nA; l++) {
          CP[i][j].A[l].getQDim(Q,S); 
          QQ.cat(1,Q);
          SS.Append(S);
       }
   }}

   if (B) {
      B->getQDim(k,Q,S);
      QQ.cat(1,Q);
      SS.Append(S);
   }

   QQ.groupRecs(P,D);
   e=SS.set2Group(P,D,0,'i'); 

   if (e) wblog(FL,
   "ERR QSpace size inconsistency (dim=%d): %s (%d)",k+1,str,e);
}

template <class TQ, class TD>
inline void QSpace<TQ,TD>::ExpandQ(
   CPAT<TQ,TD> &CP, 
   double rnrm,
   wbvector<unsigned> &xflag,
   char disp
){
   unsigned i,j,k,l,n,nr,s,r=rank(FL), dmin=0;
   const char EXISTS=1, OLDSYM=1<<1, XPAND=1<<2, XPAND2=1<<3, rank2=(r<3);
   int e=0;

   wbvector< wbMatrix<TQ> > QQ;
   wbvector< wbvector<unsigned> > SS;
   wbvector<unsigned> D(r),S(r),dc;
   wbMatrix<unsigned> SA;
   wbvector<char> mark;
   wbMatrix<TQ> Q1,Q2,QK1,QK2;
   QSpace<TQ,TD> AA;
   wbindex I1,I2,d1,d2,I(r);
   wbperm p1,p2;
   TD n1,n2;

   QQ.initDef(r);
   SS.initDef(r);

   for (i=0; i<r; i++) {
      CP.getQRange(i,QQ[i],SS[i], this); 

      S[i]=SS[i].len;

      D[i]=SS[i].sum();
      if (D[i]>1) {
      if (dmin==0) dmin=D[i]; else dmin=MIN(dmin,D[i]); }
   }

   if (xflag.isEmpty()) {
      if (!rank2) {
         xflag.init(r);
         for (i=0; i<r; i++) if (D[i]<=dmin) xflag[i]=1;
      }
   }
   else if (xflag.len!=r) wblog(FL,
   "ERR %s - invalid xflag (dim. mismatch %d/%d)",FCT,xflag.len,r);

   if (disp>9)
   wblog(FL,"TST D=[%s], expand=[%s] (%d)",
   D.toStr().data, xflag.toStr().data, dmin);

   s=S.prod(); I.set(0);
   AA.init(s,r,QDIM); SA.init(s,r); mark.init(s);

   for (l=S.len-1, i=0; i<s; i++) {
       for (k=0; k<I.len; k++) {
          AA.setQRec(i, k, QQ[k].rec(I[k]));
          SA(i,k)=SS[k][I[k]];
       }

       I[l]++; k=l; 
       while(I[k]>=S[k] && k>0) { I[k]=0; ++I[--k]; }
   }

      QIDX.blockSum(QDIM,Q1); Q1.makeUnique();
   AA.QIDX.blockSum(QDIM,Q2); Q2.groupRecs(p2,d2);

   e=matchIndex(Q1,Q2,I1,I2);
   if (e || I1.len!=Q1.dim1) wblog(FL,
      "ERR Failed to locate all existing symmetry sectors\n"
      "(e=%d; len=%d/%d)", e, I1.len, Q1.dim1);

   d2.cumsum_(dc);
   for (i=0; i<I2.len; i++) {
      l=dc[I2[i]]; n=d2[I2[i]];
      for (j=0; j<n; j++) mark[p2[l+j]] |= OLDSYM;
   }

   if (disp>9) sprintf(str,
   "Psi having %d/%d symmetry sector(s)", Q1.dim1, Q2.dim1);

   xflag.find(0,I); if (!I.isEmpty()) {
         getQsub(I,QK1); QK1.groupRecs(p1,d1);
      AA.getQsub(I,QK2); QK2.groupRecs(p2,d2);

      e=matchIndex(QK1,QK2,I1,I2);
      if (e || I1.len!=QK1.dim1) wblog(FL,
         "ERR Failed to locate all existing symmetry sectors (%d/%d)",
          I1.len, QK1.dim1);

      d2.cumsum_(dc);
      for (i=0; i<I2.len; i++) {
         l=dc[I2[i]]; n=d2[I2[i]];
         for (j=0; j<n; j++) mark[p2[l+j]] |= XPAND;
      }
   }
   else {
      for (i=0; i<mark.len; i++) mark[i] |= XPAND;
   }

   e=matchIndex(QIDX, AA.QIDX, I1, I2);
   if (e || I1.len!=QIDX.dim1) wblog(FL,
      "ERR Failed to locate all existing data blocks (%d/%d)",
       I1.len, QIDX.dim1);

   for (i=0; i<I2.len; i++) {
      mark[I2[i]] |= EXISTS;
      DATA[I1[i]]->save2(*(AA.DATA[I2[i]]));

      if (!SA.recEqual(I2[i], AA.DATA[I2[i]]->SIZE.data)) {
         SA.recPrint(i,"SA"); 
         AA.DATA[I2[i]]->SIZE.print("SIZE"); 
         wblog(FL, "ERR Severe size inconsistency (%d)", i+1);
      }
   }

   if (rnrm<=0) {
      for (i=0; i<mark.len; i++)
      if (mark[i]>0 && !(mark[i] & OLDSYM)) mark[i]=-mark[i];
   }

   for (n=i=0; i<mark.len; i++) if (mark[i]>0) n++;

   init(n,r,QDIM);
   I1.init(n); 

   for (nr=l=i=0; i<mark.len; i++) {
      if (mark[i]<=0) continue;

      QIDX.recSetP(l, AA.QIDX.rec(i));

      if (AA.DATA[i]->isEmpty()) {
         S.init2ref(SA.dim2, SA.rec(i));
         DATA[l]->init(S);
         if (!(mark[i] & OLDSYM)) {
            DATA[l]->setRand('s'); 
            I1[l]=1; nr++;
         }
      }
      else AA.DATA[i]->save2(*DATA[l]);

      l++;
   }

   if (rnrm>0 && nr) {
      double dbl,fac=0,r2=rnrm*rnrm*1E-3;
      if (fabs(rnrm)>1) wblog(FL,"WRN rnrm = %g !?",rnrm);

      for (n1=n2=0,i=0; i<I1.len; i++) {
         if (I1[i]) n2+=(DATA[i]->norm2());
         else {
            dbl=(DATA[i]->norm2()); fac+=dbl; 
            if (dbl>r2) n1+=dbl;
            else {
               DATA[i]->setRand('s'); 
               I1[i]=2; n2+=(DATA[i]->norm2());
            }
         }
      }

      if (n1!=0 && n2!=0) {
         r2=rnrm*rnrm; fac/=(n1+r2);

         if (disp>9) wblog(FL,
         "TST n1=%g, n2=%g, r2=%g -> fac=%g",n1,n2,r2,fac);

         n1=(TD)sqrt(fac*n1);
         n2=(TD)sqrt(fac*(r2/n2));

         for (i=0; i<I1.len; i++)
         if (I1[i]) (*DATA[i])*=n2; else (*DATA[i])*=n1;
      }
      else wblog(FL,"WRN rand. norm = %.3g (%.3g) !?",n2,n1);
   }

   if (xflag.count(0)) {

      unsigned itry=0,iTRY=0, ntry=16, nQ=1, m=CP.CP.len, mx=0, m0=QIDX.dim1;
      QSpace<TQ,TD> AX;

      Q1.init(); Q2.init();

      for (iTRY=0; iTRY<ntry && nQ; iTRY++) {
      for (itry=0; itry<ntry; itry++) {
         CP.getCtrIndex(*this);

         for (i=0; i<m; i++) {
            j=CP.CP[i].len-1; 
            nQ=CP.CP[i][j].CQn.len;

            for (l=0; l<nQ; l++)
            Q1.cat(1,
               CP.CP[i][j].CQn[l].QIDX 
            );
         }

         Q1.makeUnique();
         if (Q1==Q2) break; else Q1.save2(Q2);
      }

      if (itry==ntry) wblog(FL,
      "ERR CP|psi> keeps generating new terms (%d).\n"
      "Symmetries not preserved !?", ntry);
      if (Q1.dim2!=QIDX.dim2) wblog(FL,
      "ERR Severe QIDX dimension mismatch (%d/%d)",Q1.dim2,QIDX.dim2);

      e=matchIndex(Q1,AA.QIDX,I1,I2);

      if (e || I1.len!=Q1.dim1) wblog(FL,
         "ERR Failed to locate all existing data blocks (%d/%d)",
          I1.len, Q1.dim1);

      for (nQ=i=0; i<I2.len; i++) if (mark[I2[i]]<=0) nQ++;
      if (nQ) {
         l=QIDX.dim1; if (l!=DATA.len) wblog(FL,"%d/%d ???",l,DATA.len);
         QIDX.Resize(l+nQ, QIDX.dim2);
         DATA.Resize(l+nQ);
         for (i=0; i<nQ; i++) DATA[l+i] = new wbarray<TD>;

         for (i=0; i<I2.len; i++) {
            j=I2[i]; if (mark[j]>0) continue; 
            if (mark[j]<0) mark[j]=-mark[j]; else
            mark[j] |= XPAND2; mx++;

            QIDX.recSetP(l, AA.QIDX.rec(j));

            if (AA.DATA[j]->isEmpty()) {
               S.init2ref(SA.dim2, SA.rec(j));
               DATA[l]->init(S);
            }
            else {
               wblog(FL,"WRN excuse moi !?"); 
               AA.DATA[j]->save2(*DATA[l]);
            }

            l++;
         }
      }} 

      if (disp>9) wblog(FL,"TST itry=%d, iTRY=%d (%d)",itry,iTRY,ntry);

      if (iTRY==ntry) wblog(FL,
      "ERR CP|psi> keeps generating new terms (%d).\n"
      "Symmetries not preserved !?", ntry);

      if (disp>9)
      wblog(FL,"TST D=[%s], %d/%d completing to %d+%d blocks",
      getDim().toStr().data, iTRY+1,itry+1,m0,mx);
   }

   if (disp>9)
   wblog(FL,"TST %s -> %d", str, getDQtot());
}

template <class TQ, class TD>
QSpace<TQ,TD>& CPAT<TQ,TD>::HZTimes_tst(
   const QSpace<TQ,TD> &PSI,
   QSpace<TQ,TD> &HZPsi,
   const char tflag,
   char disp
){
   unsigned i,j,it;
   static unsigned call=0;
   wbperm P;

   if (!(call++)) wblog(FL,"TST %N*** using HZTimes_tst ***%N");

   QSpace<TQ,TD> X1,X2,XX;

   if (tflag=='C' && ISREAL(TD)) tflag='T';

   if (!Wb::strchri("NTC",tflag)) wblog(FL,
   "ERR Invalid flag %c<%d>", tflag, tflag);

   it=( tflag=='N' ? 2 : ( tflag=='T' ? 1 : -1 )); 

   HZPsi=PSI;
   HZPsi *= (tflag=='C' ? -CONJ(z) : -z);

   for (i=0; i<CP.len; i++) {

      if (CP[i].len<=1) { 
         CPAT_Set<TQ,TD> &C1=CP[i][0];
         for (j=0; j<C1.A.len; j++) {
            P.initFirstTo(C1.ic, PSI.rank());
            if (tflag!='C') 
                 C1.A[j].              contract(it,PSI, C1.ic+1, X1, P);
            else C1.A[j].opA(tflag,XX).contract(2, PSI, C1.ic+1, X1, P);

            HZPsi+=X1;
         }
      }
      else { 
         CPAT_Set<TQ,TD> &C1=CP[i][0], &C2=CP[i][1];
         if (C1.A.len!=C2.A.len || CP[i].len>2) wblog(FL,"ERR ???");

         for (j=0; j<C1.A.len; j++) {
            P.initFirstTo(C1.ic, PSI.rank());
            if (tflag!='C') 
                 C1.A[j].              contract(it,PSI, C1.ic+1, X1, P);
            else C1.A[j].opA(tflag,XX).contract(2, PSI, C1.ic+1, X1, P);

            P.initFirstTo(C2.ic, PSI.rank());
            if (tflag!='C') 
                 C2.A[j].              contract(it,X1 , C2.ic+1, X2, P);
            else C2.A[j].opA(tflag,XX).contract(2, X1 , C2.ic+1, X2, P);

            HZPsi+=X2;
         }
      }
   }

   return HZPsi;
};

template <class TQ, class TD>
inline void CPAT_Set<TQ,TD>::contract(
   const QSpace<TQ,TD> &PSI,
   wbvector< QSpace<TQ,TD> > &C,
   char tflag
) const {

   unsigned i,n=A.len;
   if (C.len!=n) C.initDef(n);

   for (i=0; i<n; i++)
   contract(i,PSI,C[i],tflag);
}

template <class TQ, class TD>
inline void CPAT_Set<TQ,TD>::contract_diag(
   const QSpace<TQ,TD> &PSI,
   wbvector< QSpace<TQ,TD> > &C,
   TD a, TD b, char Iflag
) const {

   unsigned i,n=A.len;
   if (C.len!=n) C.initDef(n);

   for (i=0; i<n; i++)
   contract_diag(i,PSI,C[i],a,b,Iflag);
}

template <class TQ, class TD>
inline void CPAT_Set<TQ,TD>::contract(
   const wbvector< QSpace<TQ,TD> > &PSI1,
   wbvector< QSpace<TQ,TD> > &C,
   char tflag
) const {

   unsigned i,n=A.len;
   if (C.len!=n) C.initDef(n);

   for (i=0; i<n; i++)
   contract(i,PSI1[i],C[i],tflag);
}

template <class TQ, class TD>
inline void CPAT_Set<TQ,TD>::contract_diag(
   const wbvector< QSpace<TQ,TD> > &PSI1,
   wbvector< QSpace<TQ,TD> > &C,
   TD a, TD b, char Iflag
) const {

   unsigned i,n=A.len;
   if (C.len!=n) C.initDef(n);

   for (i=0; i<n; i++)
   contract_diag(i,PSI1[i],C[i],a,b,Iflag);
}

template <class TQ, class TD>
void CPAT_Set<TQ,TD>::contract(
   const unsigned k,
   const QSpace<TQ,TD> &PSI,
   QSpace<TQ,TD> &C,
   char tflag
) const {

   unsigned i,n,ia,ib,iu,it;
   wbMatrix<unsigned> &Iabu = (tflag=='N' ? CIn[k] : CIt[k]);
   wbarray<TD> XX;

#ifdef WB_CLOCK
   Wb::Clock clk("cpat:contract",0); 
#endif

   if (tflag=='C' && ISREAL(TD)) tflag='T';

   if (!Wb::strchri("NTC",tflag)) wblog(FL,
   "ERR Invalid flag %c<%d>",tflag,tflag);

   C.QDIM=PSI.QDIM;
   C.QIDX=(tflag=='N' ? CQn[k].QIDX : CQt[k].QIDX);
   C.setupDATA();

   it=( tflag=='N' ? 2 : ( tflag=='T' ? 1 : -1 )); 

   for (n=Iabu.dim1, i=0; i<n; i++) {
      ia=Iabu(i,0); ib=Iabu(i,1); iu=Iabu(i,2);

      if (tflag!='C') {
         A[k].DATA[ia]->
         contract(it,*(PSI.DATA[ib]), ic+1, *(C.DATA[iu]), P); }
      else {
         A[k].DATA[ia]->opA(tflag,XX).
         contract(2, *(PSI.DATA[ib]), ic+1, *(C.DATA[iu]), P);
      }
   }
}

template <class TQ, class TD>
void CPAT_Set<TQ,TD>::contract_diag(
   const unsigned k,
   const QSpace<TQ,TD> &PSI,
   QSpace<TQ,TD> &C,
   TD a, TD b, char Iflag
) const {

   unsigned i,n,ia,ib,iu;
   wbMatrix<unsigned> &Iabu = CIn[k];
   wbarray<TD> XX;

#ifdef WB_CLOCK
   Wb::Clock clk("cpat:contract",0); 
#endif

   C.QDIM=PSI.QDIM; C.QIDX=CQn[k].QIDX;
   C.setupDATA();

   for (n=Iabu.dim1, i=0; i<n; i++) {
      ia=Iabu(i,0); ib=Iabu(i,1); iu=Iabu(i,2);

      if (A[k].isDiagBlock(ia))
         PSI.DATA[ib]->contractDiag(
         ic+1,*A[k].DATA[ia], a,b, *C.DATA[iu], Iflag
      );
   }
}

template <class TQ, class TD>
mxArray* CPAT<TQ,TD>::toMx() const {

   unsigned i,j,l,k, m=CP.len, n;
   mxArray *S, *C, *a;

   C=mxCreateCellMatrix(m,1);

   for (i=0; i<m; i++) {
      n=CP[i].len;
      S=mxCreateStructMatrix(1,n,0,NULL);

      for (j=0; j<n; j++) {
         a=CP[i][j].toMx();
         k=mxGetNumberOfFields(a);

         if (j==0)
         for (l=0; l<k; l++) {
            if (mxAddField(S,mxGetFieldNameByNumber(a,l))<0) wblog(FL,
            "ERR Could not add field `%s'", mxGetFieldNameByNumber(a,l));
         }

         for (l=0; l<k; l++) { 
            mxSetFieldByNumber(S, j, l,
            mxDuplicateArray(mxGetFieldByNumber(a,0,l)));
         }

         mxDestroyArray(a);
      }
      mxSetCell(C,i,S);
   }

   S=mxCreateStructMatrix(1,1,0,NULL);
   mxAddField2Scalar(FL,S,"CP",C);
   mxAddField2Scalar(FL,S,"z",numtoMx(z));

   return S;
}

template <class TQ, class TD>
mxArray* CPAT_Set<TQ,TD>::toMx() const {

   mxArray *S, *a;
   unsigned i;

   S=mxCreateStructMatrix(1,1,0,NULL);

   mxAddField2Scalar(FL,S,"A",    A.toMx()     );
   mxAddField2Scalar(FL,S,"ic",   numtoMx(ic+1));
   mxAddField2Scalar(FL,S,"opA",  numtoMx(opA) );
   mxAddField2Scalar(FL,S,"P",    (P+1).toMx() );

   for (i=0; i<CQn.len; ++i) CQn[i].setupDATA();
   for (i=0; i<CQt.len; ++i) CQt[i].setupDATA();

   mxAddField2Scalar(FL,S,"CQn", CQn.toMx());
   mxAddField2Scalar(FL,S,"CQt", CQt.toMx());

   a=mxCreateCellMatrix(CIn.len,1);
   for (unsigned i=0; i<CIn.len; i++)
   mxSetCell(a, i, CIn[i].toMx());
   mxAddField2Scalar(FL,S,"CIn",a);

   a=mxCreateCellMatrix(CIt.len,1);
   for (unsigned i=0; i<CIt.len; i++)
   mxSetCell(a, i, CIt[i].toMx());
   mxAddField2Scalar(FL,S,"CIt",a);

   return S;
}

int isCPAT(const mxArray *CC, char twoSite) {
   unsigned i,i0,k,n;

   if (!CC || !mxIsCell(CC)) return 0;
   if (isCPattern(CC)) return 1;

   n=mxGetNumberOfElements(CC);
   if (!n) return 0;

   i0=isCPattern(mxGetCell(CC,0));
   if (!i0) return 0;

   for (k=1; k<n; k++) {
      i=isCPattern(mxGetCell(CC,k));
      if (!i) return 0;

      if (i!=i0) wblog(FL, 
      "ERR CPAT mismatch twoSite=%d,%d ???", i0, i);
   }

   if (twoSite>=0) if (1+twoSite!=(char)i0) wblog(FL,
     "ERR CPAT is not of type %d-site (%d)", 1+twoSite, i0);

   return n; 
}

int isCPattern(const mxArray *C) {
   unsigned k,n,w;
   char twoSite, hasOp;
   mxArray *a;

   if (!C || !mxIsCell(C)) {
      sprintf(str,"%s:%d invalid CP",FL);
      return 0;
   }

   n=mxGetNumberOfElements(C);
   if (n<2 || n>4) {
      sprintf(str,"%s:%d invalid CP (%d entries)",FL,n);
      return 0;
   }

   a=mxGetCell(C,0);
   if (!mxIsStruct(a) || mxIsCell(a)) {
      sprintf(str,"%s:%d invalid arg #1",FL);
      return 0;
   }

   a=mxGetCell(C,n-1);
   hasOp=mxIsChar(a);

   if (hasOp) if (mxGetNumberOfElements(a)!=1) {
      sprintf(str,"%s:%d arg #%d not a flag (char)",FL,n);
      return 0;
   }

   k=n-hasOp;
   if (k==2) twoSite=0; else
   if (k==3) twoSite=1; else {
      sprintf(str,"%s:%d invalid CP (k=%d)",FL,k);
      return 0;
   }

   if (twoSite) {
      w = (mxIsEmpty(mxGetCell(C,1)) ? 0 : 1) +
          (mxIsEmpty(mxGetCell(C,2)) ? 0 : 2);
      if (w<1 || w>2) {
         sprintf(str, "%s:%d either arg #2 or #3",FL);
         return 0;
      }
   }
   else w=1;

   a=mxGetCell(C,w);
   if (!Mx::IsNumber(a)) {
      sprintf(str,"%s:%d arg #%d not a number",FL,w+1);
      return 0;
   }

   return (1+twoSite);
}

#endif

