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

#ifndef __WB_QSPACE_COL_MAJOR_CC__
#define __WB_QSPACE_COL_MAJOR_CC__

// NB! OpenMP requires that a loop construct processes EACH iteration
// => breaking out of loops through intermediate return, break, throw,
// error, wblog(ERR,...) statements, etc. is NOT allowed! // Wb,Aug04,14
// e.g. see http://bisqwit.iki.fi/story/howto/openmp/
// => using try/catch statement, query shared variable (e.g. ex below)
// to quickly flush remainder of loop, immediately issuing an error
// once outside the parallel section.  // OMP_LOOP_CATCH
// see also ~/Source/MPI/tst_openmp.cc // Wb,Jan24,17

template<class TQ, class TD>
void QSpace<TQ,TD>::init_itags(
   const char *F, int L, const char *tag, unsigned r_, char lflag
){
   if (QIDX.dim2 && QDIM) {
      unsigned l, r=QIDX.dim2/QDIM; if (int(r_)<0 && r) { r_=-r; }

      l=itags.Init(F,L,tag,r_,lflag);

      if (QDIM) { if (r && r!=l) wblog(F_L,
         "ERR %s() mismatch in number of itags (%d/%d)",
          FCT,r,QIDX.dim2/QDIM);
      }
      else if (r || QIDX.dim1 || QIDX.dim2) wblog(FL,
         "ERR %s() got invalid labels `%s' having (%dx%d;%d)",
          FCT,IT2STR__,QIDX.dim1,QIDX.dim2,QDIM
      );
   }
};

template<class TQ, class TD>
void QSpace<TQ,TD>::initOType(const char *F, int L, const mxArray* a){

   if (!a) { otype=QS_NONE; return; }
   else {
      unsigned i=0; wbstring s(F,L,a);
      if (!s.data || !s.data[0]) { otype=QS_NONE; return; }

      for (i=1; i<QS_NUM_TYPES; ++i)
      if (!strcmp(s.data,QS_STR[i])) { otype=QS_TYPES(i); return; }

      if (i==QS_NUM_TYPES) wblog(FL,
      "ERR invalid object type '%s'",s.data);
   }
};

template<class TQ, class TD>
wbstring QSpace<TQ,TD>::otype2Str(const char *fmt) const {

   const char *q=0;

   switch (otype) {
     case QS_NONE     : 
     case QS_AMATRIX  : 
     case QS_OPERATOR : q=QS_STR[otype]; break;
     default :
        char s[16]; snprintf(s,16,"(invalid: %d)",otype);
        return s;
   }
   if (fmt && fmt[0] && q && q[0]) {
      char s[16]; snprintf(s,16,fmt,q);
      return s;
   }
   return q;
};

template<class TQ, class TD>
template<class TB>
unsigned QSpace<TQ,TD>::matchITags(
   const char *F, int L,
   const QSpace<TQ,TB> &B, ctrIdx &ia, ctrIdx &ib
 ) const {

   unsigned n=0, ra=rank(FL), rb=B.rank(FL);

   if (itags.len!=ra || B.itags.len!=rb) wblog(F_L,
      "ERR %s() got empty or invalid set\n"
      "'%s' (%d/%d) <> '%s' (%d/%d)",FCT,IT2STR__,itags.len,ra,
      IT2STR(B), B.itags.len, rb
   );

   n=itags.getCtr(F,L,B.itags,ia,ib);

   if (!n) { if (F) wblog(F,L,
     "ERR %s() got empty match\n"
     "'%s' <> '%s'", FCT,IT2STR__, IT2STR(B)
   ); }

   return ia.len;
};

template<class TQ, class TD>
void QSpace<TQ,TD>::PrependSingletons(unsigned R) {

    unsigned r=-1;

    if (!isConsistent(r)) { wbdie(FL,str); }
    if (r<R) {
       unsigned i,s=QIDX.dim2*sizeof(TQ);
       wbMatrix<TQ> Q0(QIDX);
       TQ *q, *q0=Q0.data;

       if (itags.len) wblog(FL, 
          "ERR %s() also need to adjust itags!  (%d/%d-%d)",
          FCT,itags.len,r,R);

       QIDX.init(QIDX.dim1, R*QDIM); q=QIDX.data+(R-r)*QDIM;
       for (i=0; i<QIDX.dim1; ++i) {
          memcpy(q,q0,s); q+=QIDX.dim2; q0+=Q0.dim2;
          DATA[i]->prependSingletons(R);
       }
    }
    else if (r>R) {
    wblog(FL,"ERR %s() index reduction (%d->%d) !?",FCT,r,R); }
};

template <class TQ, class TD> inline
QSpace<TQ,TD>& QSpace<TQ,TD>::setupCGR(const char *F, int L) {

    if (QIDX.dim1 && QIDX.dim2) { 
    if (!QDIM || QIDX.dim2%QDIM) wblog(F_L,
       "ERR invalid QIDX (%dx%d; %d)",QIDX.dim1,QIDX.dim2,QDIM);
    if (qtype.len && qtype.Qlen()!=QDIM) wblog(F_L,
       "ERR qtype inconsistency (%d/%d)",QDIM,qtype.Qlen());
    }

    if (QIDX.dim1 && qtype.len) {
       CGR.init(QIDX.dim1,qtype.len);
    }
    else CGR.init();

    return *this;
};

template <class TQ, class TD> inline
QSpace<TQ,TD>& QSpace<TQ,TD>::initIdentity(
   const QSpace<TQ,TD> &A, char dflag) {

   wbMatrix<TQ> Q1;
   wbvector<widx_t> D; wbperm P;
   unsigned i,d,l;

   if (dflag && !strchr("rc",dflag)) wblog(FL,
      "ERR %s() invalid dflag=%c<%d>",FCT,dflag,dflag);

   clearQSpace();

   if (A.isEmpty()) { return *this; }
   if (!A.isConsistent_r(2)) { wbdie(FL,str); }

   if (!A.isQSym(0,0,dflag)) {
      if (dflag)
           wblog(FL,"WRN %s() got non-compact QSpace !?\n%s",FCT,str);
      else wblog(FL,"WRN %s() got non-symmetric QIDX !?\n%s",FCT,str);
   }

   A.QIDX.getBlock(0,A.QDIM,Q1);
   Q1.groupRecs(P,D);

   qtype=A.qtype;
   QDIM=Q1.dim2; QIDX.Cat(2,Q1,Q1); setupDATA();

   for (l=i=0; i<DATA.len; ++i, l+=D[i]) {
      d=A.DATA[P[l]]->SIZE.max(); 
      DATA[i]->initIdentity(d,dflag);
   }

   return initIdentityCGS();
};

template <class TQ, class TD> inline
QSpace<TQ,TD>& QSpace<TQ,TD>::initIdentity_bare(const wbMatrix<TQ> &Q) {

   clearQSpace();
   QDIM=Q.dim2; QIDX.Cat(2,Q,Q); setupDATA();

   for (unsigned i=0; i<DATA.len; ++i) {
      DATA[i]->init(1,1);
      DATA[i]->data[0]=1;
   }
   return *this;
};

template <class TQ, class TD> inline
QSpace<TQ,TD>& QSpace<TQ,TD>::initIdentity(const QType &t, const TQ *qs) {

   unsigned n=t.qlen();
   if (!n || !t || !qs) wblog(FL,
      "ERR %s() invalid input '%s' (n=%d, %p)",FCT,STR(t),n,qs);

   clearQSpace();
   QDIM=n; QIDX.init(1,2*n).recSetP(0,qs,n,qs,n);

   setupDATA();
      DATA[0]->init(1,1);
      DATA[0]->data[0]=1;

   qtype.init(t);
   if (!t.isAbelian()) { initIdentityCGS(); }

   return *this;
};

template <class TQ, class TD> inline
QSpace<TQ,TD>& QSpace<TQ,TD>::init(
   const CDATA_TQ &Cb, wbperm *P, char inv) {

   if (!Cb.valid()) { 
      wblog(FL,"ERR %s() got %s CGT %s",
      FCT, Cb.NP_zero()? "non-permissible":"invalid", STR(Cb));
   }

   unsigned r=Cb.rank(FL), m=0; clearQSpace(); 
   if (!Cb.t || !Cb.cstat.gotuser_BUF()) wblog(FL,
      "ERR %s() invalid input: %s",FCT,STR2(Cb,2));
   QDIM=Cb.t.qlen();
   QIDX.init(1,Cb.qs.len,Cb.qs.data);
   itags.init(Cb.qdir);

   if (Cb.t.isU1())
        { qtype.init(); }
   else { qtype.init(Cb.t); }

   setupDATA(FL);
   if (!Cb.t.isAbelian()) {
      setupCGR(FL); CGR[0].initBase(&Cb); 
      m=Cb.numOM(FL); 
   }
   DATA[0]->initSingleton(r,1,m);

   if (P && *P) { Permute(*P); }

   return *this;
};

template <class TQ, class TD> inline
QSpace<TQ,TD>& QSpace<TQ,TD>::FusionTree(
   const char *F, int L,
   const QType &t, const wbMatrix<TQ> &Q, QSpace<TQ,TD> *AK) {

   if (!t) { wblog(F_L,"ERR %s() invalid symmetry '%s'",FCT,STR(t)); }
   if (Q) {
      if (Q.dim2!=t.qlen()) wblog(F_L,
         "ERR %s() size mismatch (%dx%d / %d)",FCT,Q.dim1,Q.dim2,t.qlen());
      initIdentity(t,Q.data);
   }
   else { init(); qtype.init(t); }

   if (AK) { AK->init(*this); }
   if (Q.dim1<=1) { return *this; }

   QSpace<TQ,TD> X;
   wbindex Ia(1), Ib(1), i3(1);
   wbMatrix<TQ> Qa,Qb;
   wbvector<widx_t> Sa,Sb;
   QMap<TQ> M;

   QSpace<TQ,TD> A3, Ak(*this), Ek;

   i3[0]=2; 

   Wb::iterLevel<int> wd(&CG::thread_fOM.me(),256);

   for (unsigned k=1; k<Q.dim1; ++k) {
      getQDimGen( cPVEC1_(Ak), Qa,Sa,Ia); Ek.initIdentity(t,Q.rec(k));
      getQDimGen( cPVEC1_(Ek), Qb,Sb,Ib);
      gCS.getQfinal(FL,qtype,Qa,Qb,M,'c'); 
      M.getIdentityQ(FL,Sa,Sb,A3);

      if (AK) {
      AK->contract(FL,k+1,A3,1,X); X.save2(*AK); }
      Ak=initIdentityCG( cPVEC1_(A3), i3);
   }

   return Ak.save2(*this);
};

template <class TQ, class TD> inline
QSpace<TQ,TD>& QSpace<TQ,TD>::initDiagonal(
   const wbMatrix<TQ> &Q,
   const wbvector<unsigned> &D, 
   const wbvector<TD> &E,       
   char rc                      
){
   unsigned i,j,d,l;

   if (D.len!=Q.dim1 || D.sum()!=E.len) wblog(FL,
      "ERR %s() size mismatch (D.len=%d/%d; sum(D)=%d/%d)",
       FCT, D.len, Q.dim1, D.sum(), E.len);

   if (rc && !strchr("rc",rc)) wblog(FL,
      "ERR %s() invalid rcflag=%s",FCT,cSTR(rc));

   clearQSpace(); if (Q.isEmpty()) { return *this; }

   QIDX.Cat(2,Q,Q); QDIM=Q.dim2; setupDATA();

   for (l=i=0; i<D.len; ++i, l+=d) {
      wbarray<TD> &A=(*DATA[i]); d=D[i];

      if (rc) {
         if (rc=='r') A.init(1,d); else A.init(d,1);
         memcpy(A.data, E.data+l, d*sizeof(TD));
      }
      else {
         A.init(d,d);
         for (j=0; j<d; ++j) A.data[j+j*d]=E[l+j];
      }
   }
   return *this;
};

template <class TQ, class TD> inline
char QSpace<TQ,TD>::checkQ_CGR(const char *F, int L, cUVEC *I) const {

    char q=0; 

    if (!QIDX.dim1 || !QIDX.dim2) { return q; }
    if (!CGR.data) {
       if (!qtype.allAbelian()) wblog(FL,
          "ERR %s() got missing CGR data (%s)",FCT,STR(qtype));
       return q;
    }

    if (CGR.dim1!=QIDX.dim1 || CGR.dim2!=qtype.len) wblog(FL,
       "ERR %s() got CGR size mismatch (%s <> %s; %s)",
       FCT,SSTR(QIDX),SSTR(CGR),STR(qtype));
    if (!QDIM || QIDX.dim2%QDIM) wblog(F_L,
       "ERR %s() invalid QIDX (%dx%d; %d)",FCT,QIDX.dim1,QIDX.dim2,QDIM);

    unsigned i,j, k=0, n=(I? I->len : CGR.dim1);
    unsigned r=QIDX.dim2/QDIM; 
    wbvector<unsigned> qdc;
    QSet<TQ> Q;

    if (qtype.Qpos(qdc)!=QDIM) {  
       wblog(F_L,"ERR %s() qtype inconsistency (%s: %d/%d)",
       FCT, qStr().data,qtype.Qlen(),QDIM);
    }

    for (; k<n; ++k) { i=(I? I->data[k] : k);
       const TQ *qi=QIDX.rec(i); 
       for (j=0; j<CGR.dim2; ++j) { if (CGR(i,j).cgb) {
          Q.init(qtype[j],qi+qdc[j],r,QDIM,itags);
          q=CGR(i,j).checkQ(F,L,Q);
          if (q) {
             return q;
          }
       }}
    }
    return q;
};

template <class TQ, class TD> inline
QSpace<TQ,TD>& QSpace<TQ,TD>::initIdentityCGS(const char *F, int L) {
    unsigned i,j, r=2;

    if (!QIDX) { CGR.init(0,qtype.len); return *this; }

    if (!QDIM || QIDX.dim2%QDIM) wblog(F_L,
       "ERR %s() invalid QIDX (%dx%d; %d)",FCT,QIDX.dim1,QIDX.dim2,QDIM);
    if (isref) wblog(FL,"ERR %s() got isref (bailing out)",FCT);
    if (QIDX.dim2/QDIM!=r) wblog(F_L,
       "ERR %s() for rank-%d objects only (%dx%d/%d; %s)",
        FCT,r,QIDX.dim1,QIDX.dim2,QDIM,qStr().data
    );

    if (qtype.len) {
       CGR.init(QIDX.dim1,qtype.len);
       wbvector<unsigned> qdc;

       if (qtype.Qpos(qdc)!=QDIM) 
          wblog(F_L,"ERR %s() qtype inconsistency (%s: %d/%d)",
          FCT, qStr().data,qtype.Qlen(),QDIM
       );

       for (i=0; i<CGR.dim1; ++i) {
          const TQ *qi=QIDX.rec(i);
          if (memcmp(qi,qi+QDIM,QDIM*sizeof(*qi))) wblog(FL,
             "ERR %s() got non-block-diagonal rank-2 QSpace !?\n"
             "%d: %s",FCT,i+1,QIDX.rec2Str(i).data);
          for (j=0; j<CGR.dim2; ++j) {
             CGR(i,j).initIdentityR(FL, qtype[j], qi+qdc[j]);
          }
       }
    }
    else { CGR.init(); }

    if (!itags.len) { itags.init_qdir("+-"); }
    else if (itags.len==2 && itags[0].isConj(itags[1],0)) {
       if (itags[0].isConj()) { 
          wbperm p("21");
          for (unsigned n=CGR.numel(), i=0; i<n; ++i) {
             CGR[i].Permute(p);
          }
       }
    }
    else wblog(FL,
      "ERR %s() unexpected itags %s (len=%d/%d)",FCT,STR(itags),itags.len,r);
    return *this;
};

template <class TQ, class TD> inline
QSpace<TQ,TD>& QSpace<TQ,TD>::initIdentityCGS(
    const char *F, int L, unsigned i, unsigned r
){
    unsigned j;

    if (isref) wblog(FL,"ERR %s() got isref (bailing out)",FCT);
    if (QDIM==0 || QIDX.dim2%QDIM) wblog(F_L,
       "ERR %s() invalid QIDX (%dx%d; %d)",FCT,QIDX.dim1,QIDX.dim2,QDIM);
    if (r>=(j=QIDX.dim2/QDIM)) wblog(F_L,
       "ERR %s() index out of bounds (%d/%d; %dx%d/%d; %s)",
        FCT,r,j, QIDX.dim1, QIDX.dim2, QDIM, qStr().data);

    if (i>=CGR.dim1) {
       if (!CGR.dim1 && !qtype.len) return *this;
       else wblog(FL,
          "ERR %s() index out of bounds (%d/%dx%d; %s)",
          FCT,i, CGR.dim1, CGR.dim2, qStr().data
       );
    }
    else {
       const TQ *qi=QIDX.ref(i,r*QDIM);
       wbvector<unsigned> qdc;

       if (qtype.Qpos(qdc)!=QDIM) 
          wblog(F_L,"ERR %s() qtype inconsistency (%s: %d/%d)",
          FCT, qStr().data, qtype.Qlen(), QDIM
       );

       for (j=0; j<CGR.dim2; ++j) {
          CGR(i,j).initIdentityR(FL, qtype[j], qi+qdc[j]);
       }
    }

    return *this;
};

template <class TQ, class TD>
template<class TD2> inline
QSpace<TQ,TD>& QSpace<TQ,TD>::initIdentityCGS(
   const char *F, int L, const QSpace<TQ,TD2> &H
){
    unsigned i,j, r=2;
    wbindex I1,I2;

    if (!QIDX.dim1 || !QIDX.dim2) {
       CGR.init(0,qtype.len); return *this;
    }

    if (QDIM==0 || QIDX.dim2%QDIM) wblog(F_L,
       "ERR %s() invalid QIDX (%dx%d; %d)",FCT,QIDX.dim1,QIDX.dim2,QDIM);
    if (qtype.isEmpty()) wblog(F_L,"ERR %s() got empty qtype",FCT);
    if (isref) wblog(FL,"ERR %s() got isref (bailing out)",FCT);
    if (qtype.Qlen()!=QDIM) wblog(F_L,
       "ERR %s() qtype inconsistency (%s: %d/%d)",
        FCT, qStr().data,qtype.Qlen(),QDIM);
    if (QIDX.dim2/QDIM!=r) wblog(F_L,
       "ERR %s() for rank-%d objects only (%dx%d/%d; %s)",
        FCT,r,QIDX.dim1,QIDX.dim2,QDIM,qStr().data);
    if (qtype!=H.qtype) wblog(F_L,"ERR %s() mismatch in qtype (%s <> %s)",
        FCT, qStr().data, H.qStr().data);

    i=(unsigned)matchIndex(QIDX, H.QIDX, I1, I2);
    if (i || QIDX.dim1!=I1.len) wblog(FL,
       "ERR %s() QIDX not fully contained in reference space (%d,%d)",
       FCT, I1.len, QIDX.dim1
    );

    CGR.init(QIDX.dim1,qtype.len);

    for (i=0; i<CGR.dim1; ++i) {
    for (j=0; j<CGR.dim2; ++j) { CGR(I1[i],j)=H.CGR(I2[i],j); }}

    return *this;
};

template <class TQ, class TD> inline
void QSpace<TQ,TD>::ExpandDiagonal(
   unsigned i1, 
   unsigned i2  
){
   unsigned i,j1,j2, n=QDIM*sizeof(TQ), dflag=0, r=-1;

   isConsistent(FL,r);
   if (!i1 || !i2 || i1>r || i2>r || i1==i2) wblog(FL,
      "ERR %s() invalid index (%d %d; %d)",FCT,i1,i2,r);

   i1--; i2--; 

   j1=i1*QDIM;
   j2=i2*QDIM;

   for (i=0; i<DATA.len; ++i) {
      if (memcmp(QIDX.ref(i,j1), QIDX.ref(i,j2),n)) { 
         continue;
      }

      const wbvector<unsigned> &S=DATA[i]->SIZE;

      if (S[i1]==S[i2]) continue;

      if (S[i1]!=1) {
         if (dflag) { if (dflag!='c') wblog(FL,
            "ERR %s() expecting `col'\n%s (%d,%d)",FCT, SSTR(S),i1+1,i2+1);
         } else dflag='c';
      }
      else
      if (S[i2]!=1) {
         if (dflag) { if (dflag!='r') wblog(FL,
            "ERR %s() expecting `row'\n%s (%d,%d)",FCT,SSTR(S),i1+1,i2+1);
         } else dflag='r';
      }

      DATA[i]->ExpandDiagonal(i1,i2);
   }
};

template <class TQ, class TD>
QSpace<TQ,TD>& QSpace<TQ,TD>::init2DiffOp(
   const QSpace<TQ,TD> &A, unsigned ia, 
   const QSpace<TQ,TD> &B, unsigned ib
){
   wbvector<widx_t> d1,d2,D;
   wbMatrix<TQ> Qa,Qb,Qx;
   wbindex Ia,I0;
   wbperm P1,P2;

   A.getQsub(ia,Qa).groupRecs(P1,d1); I0.initGroup0(P1,d1);
   B.getQsub(ib,Qb).groupRecs(P2,d2);

   Qa.getDiffSorted(Qb,Ia);
   if (!Ia.len) { init(); return *this; }

   unsigned i,j,d; wbvector<TD> E;

   Qa.getRecs(Ia,Qx); I0.Select(Ia);

   D.init(Ia.len);
   for (i=0; i<I0.len; ++i) { D[i]=A.DATA[I0[i]]->SIZE.elx(ia,1); }
   E.init(D.sum()); 

   initDiagonal(Qx,D,E,0); 

   qtype=A.qtype; if (qtype.len) {
      wbvector<unsigned> qdc; qtype.Qpos(qdc);

      setupCGR();
      if (CGR.dim2!=A.CGR.dim2) wblog(FL,
         "ERR %s() size mismatch (%d/%d)",FCT,CGR.dim2,A.CGR.dim2);

      for (i=0; i<I0.len; ++i) {
         const CRef<TQ> *cg=A.CGR.rec(I0[i]);
         const TQ *qi = A.QIDX.ref(I0[i],ia*A.QDIM);
         for (j=0; j<CGR.dim2; ++j) {
            d=cg[j].Size(ia); 
            CGR(i,j).initIdentityR(FL, qtype[j], qi+qdc[j], d,
              cg[j].cgb ? '!': 0 
            );
         }
      }
   }
   else if (B.qtype.len) wblog(FL,
      "ERR %s() qtype mismatch (%s,%s)",FCT,
       A.qStr().data, B.qStr().data
   );

   itags.init_qdir("+-"); 

   return *this;
};

template <class TQ, class TD> inline
int QSpace<TQ,TD>::gotCGS(const char *F, int L,
   char xflag 
 ) const {

   int m=0; 

   if (CGR.isEmpty()) {
      if (qtype.allAbelian()
         || (!QIDX.dim1 && DATA.len==0) 
         || (!QIDX.dim2 && DATA.len==1) 
         ){ return m; }
      if (!QIDX.norm2()) { return m; }  
      wblog(F_L,"ERR %s() empty CGR for %s !?",FCT,qStr().data);
   }

   if (qtype.isEmpty() || CGR.dim2!=qtype.len) wblog(F_L,
      "ERR %s() got invalid qtype (%d/%d)",FCT,qtype.len,CGR.dim2);

   if (CGR.dim1!=QIDX.dim1 || 
       CGR.dim2!=qtype.len || DATA.len!=QIDX.dim1) wblog(FL,
      "ERR %s() QSpace inconsistency\n(CGR: %s <> %dx%d; %d)",
       FCT, SSTR(CGR), QIDX.dim1, qtype.len, DATA.len);

   if (!QDIM || QIDX.dim2%QDIM || QDIM!=qtype.Qlen()) wblog(F_L,
      "ERR %s() invalid QDIM (QIDX.dim2: %s @ %d/%d)",
      FCT, SSTR(QIDX), QDIM, qtype.Qlen());

   unsigned i,j,m_,nw=0; QDir qdir(itags);
   try {
      for (i=0; i<CGR.dim1; ++i)
      for (j=0; j<CGR.dim2; ++j) { const CRef<TQ> &Rij=CGR(i,j);
         Rij.check(FL,qtype[j],qdir);
         if (Rij.cgw || !qtype[j].isAbelian()) {
            m_=Rij.wdim1(FL);    
            if (m<int(m_)) { m=m_; }  
            if (m_==1 && Rij.cgw[0]!=1) { ++nw; }
         }
      }
   } catch (...) {
      wblog(F_L,"ERR %s() CGR(%d,%d)\nhaving %s [%dx%d]",
      FCT,i+1,j+1, qStr().data, CGR.dim1,CGR.dim2);
   }

   if (allAbelian()) {
      if (nw) wblog(FL,"ERR %s() "
         "non-trivial cgw for all-abelian %s (%d/%d)",FCT,STR(qtype),m,nw);
      if (m==1) { m=-1; } else
      if (m) wblog(FL,
        "ERR %s() got m=%d/%d for all-abelian %s",FCT,m,nw,STR(qtype));
   }
   else if (!m) { wblog(FL,"ERR %s() got m=%d for %s",FCT,m,STR(qtype)); }
   else if (xflag) {
      if (m==1 && !nw) { m=-2; }
   }
   else {
      m=1+qtype.permitsOM(rank(FL));
   }

   return m;
};

template <class TQ, class TD> inline
int QSpace<TQ,TD>::gotCGX(const char *F, int L) const {

   unsigned i,j;
   for (j=0; j<CGR.dim2; ++j)
   for (i=0; i<CGR.dim1; ++i) {
      if (CGR(i,j).cgp.len || CGR(i,j).conj) {
         if (!CGR(i,j).cgb) wblog(F_L,
            "ERR %s() unexpected CRef(%d,%d) data (%s)\n%s",
            FCT,i+1,j+1,qStr().data,STR(CGR(i,j)));
         if (!CGR(i,j).isScalar('d')) {
            return 1;
         }
      }
   }
   return 0;
};

template <class TQ, class TD>
template <class TB>
bool QSpace<TQ,TD>::sameType(
   const QSpace<TQ,TB> &B, int r, const char *istr
 ) const {

   if (istr) str[0]=0;
   if ((void*)this==(void*)&B) { if (r<0) return 1; }
   else {
      if (QIDX.dim2!=B.QIDX.dim2 || QDIM!=B.QDIM) {
         if (istr) sprintf_str("%s: %s QIDX: %d/%d, QDIM: %d/%d",
         SHORT_FL, istr, QIDX.dim2, B.QIDX.dim2, QDIM, B.QDIM);
         return 0;
      }
      if (!qtype.sameType(B.qtype)) { 
         if (istr) sprintf_str("%s: %s qtype: %s/%s",
         SHORT_FL,istr, qStr().data, B.qStr().data);
         return 0;
      }
      if (QDIM!=qtype.Qlen()) { 
         if (istr) { sprintf_str("%s: %s QDIM/qtype: %d/%d",
             SHORT_FL,istr,QDIM,qtype.Qlen());
             return 0;
         }
         else wblog(FL,"ERR qtype not compatible with QDIM (%d/%d)",
         QDIM,qtype.Qlen());
      }
   };

   if (r<0) return 1;

   if (r==0) {
      if (QIDX.dim2) { if (istr) sprintf_str(
         "%s: %s rank = %d/0",SHORT_FL,istr,QIDX.dim2);
          return 0;
      }
      else return 1; 
   }
   else if (QDIM==0 || QIDX.dim2/QDIM!=(unsigned)r) {
      if (istr) sprintf_str("%s: %s rank = %d/%d",
          SHORT_FL,istr,QIDX.dim2/QDIM,r);
      return 0;
   }

   return 1;
};

template <class TQ, class TD> inline
unsigned QSpace<TQ,TD>::getOM(
   wbMatrix<unsigned> &OM, wbvector<unsigned> *M_,
   char wflag 
 ) const {

   unsigned Mtot=0, nsym=(qtype.len ? qtype.len : QDIM);

   if (!QIDX.isUnique()) wblog(FL,
      "ERR %s() got non-unique QIDX",FCT);
   OM.init2val(QIDX.dim1,nsym,1);

   if (!CGR.dim1) { 
      if (M_) { M_->init2val(QIDX.dim1,1); }
      return Mtot;
   }
   if (CGR.dim2!=nsym) wblog(FL,
      "ERR %s() CGR/qtype inconsistency (%s/%d)",FCT,SSTR(CGR),nsym);

   unsigned i,j,m;

   for (i=0; i<CGR.dim1;  ++i) {
   for (m=1, j=0; j<nsym; ++j) { m*=(OM(i,j)=CGR(i,j).getOM(FL,wflag)); }
      Mtot+=(m-1);
   }

   if (wflag!=2) { if (M_) { OM.recProd(*M_); }}
   else {
      unsigned r=rank(FL);
      wbvector<unsigned> M; OM.recProd(M);

      if (M.len!=DATA.len) { wblog(FL,
         "ERR %s() size mismatch len=%d/%d",FCT,M.len,DATA.len); }
      for (i=0; i<DATA.len; ++i) { m=DATA[i]->getSizeM(r);
         if (m!=M[i]) wblog(FL,"ERR %s() got "
         "OM mismatch in DATA[%d] @ %d / %d",FCT,i+1,DATA.len,m,M[i]);
      }
      if (M_) { M.save2(*M_); }
   }

   return Mtot;
};

template <class TQ, class TD> inline
void QSpace<TQ,TD>::getDim( 
   wbvector<widx_t> &D, wbvector<widx_t> *DD
 ) const {

   unsigned k, r=-1;

   if (!isConsistent(r)) {
      wblog(FL,"ERR %s() %s",FCT,str);
   }
   if (r<0) wblog(FL,"ERR invalid rank r=%d",r);

   if (otype==QS_OPERATOR) { r=rank(FL); }

   D.init(r);
   if (DD) { DD->init(r);
          for (k=0; k<r; ++k) D[k]=getDim(k,DD->data+k); }
   else { for (k=0; k<r; ++k) D[k]=getDim(k); }
};

template <class TQ, class TD> inline
unsigned QSpace<TQ,TD>::getDim(unsigned k, widx_t *D_) const {

   widx_t i,s, l=0, Dm=0, Ds=0;
   wbvector<widx_t> dd;
   wbMatrix<TQ> Qk;
   wbperm pp;

   char cgflag=gotCGS(FL);

   if (QIDX.dim1==0) { if (D_) *D_=0; return 0; }

   getQsub(k,Qk).SkipTiny_float().groupRecs(pp,dd);

   if (cgflag==1 && qtype.permitsOM()>1) { cgflag=2; }
   if (cgflag>0) {
      unsigned j,ip,e=0,dc, r=rank(FL), m=0;
      for (i=0; i<qtype.len; ++i) { if (qtype[i].permitsOM(r)) { ++m; }}

      if (k>=r) wblog(FL,
         "ERR %s() index out of bounds (%d/%d)",FCT,k+1,r);

      for (l=i=0; i<dd.len; ++i, l+=dd[i-1]) { ip=pp[l];
         const wbarray<TD> &a=*DATA[ip];
         if (!a.isRankM(r,m)) wblog(FL, 
            "ERR %s() rank mismatch (%s / %d)",FCT,SSTR(a),r);
         s=a.SIZE[k]; Dm+=s;

         for (dc=1,j=0; j<CGR.dim2; ++j) {
            dc*=CGR(ip,j).Size(k);

         }
         if (cgflag<=1) { 
            unsigned d0=qtype.QDim(QIDX.ref(ip,QDIM*k)); if (d0!=dc)
            wblog(FL,"ERR %s() cgsDim inconsistency! (%d/%d)",FCT,dc,d0);
         }
         Ds+=(s*dc);
      }

      if (e) wblog(FL,
      "WRN %s() CGR data not yet initialized !? (%d)",FCT,e);
   }
   else {
      for (l=i=0; i<dd.len; ++i, l+=dd[i-1]) { 
         Dm+=(s=DATA[pp[l]]->SIZE[k]);
         if (qtype.len) Ds += s * qtype.QDim(Qk.rec(i));
      }
   }

   if (D_) { D_[0] = (Ds ? Ds : Dm); }

   return Dm;
};

template <class TQ, class TD> inline
unsigned QSpace<TQ,TD>::cgsDim(unsigned i, unsigned k) const {

   if (i>=QIDX.dim1) wblog(FL,
      "ERR %s() index out of bounds (%d/%d)",FCT,i,QIDX.dim1);

   unsigned d=1, R=rank(FL);
   if (k>=R) wblog(FL,"ERR %s() index out of bounds (%d/%d)",FCT,k,R);
   if (!CGR.dim2) return d;

   if (CGR.dim1!=QIDX.dim1) wblog(FL,
      "ERR %s() size mismatch (%dx%d <> %dx%d)",
       FCT,QIDX.dim1,QIDX.dim2,CGR.dim1,CGR.dim2);

   for (unsigned j=0; j<CGR.dim2; ++j) { d*=CGR(i,j).Size(k); }

   return d;
};

template <class TQ, class TD> inline
unsigned QSpace<TQ,TD>::cgsDimScalar(unsigned i) const {

   if (i>=QIDX.dim1) wblog(FL,
      "ERR %s() index out of bounds (%d/%d)",FCT,i,QIDX.dim1);
   unsigned d=1, r=rank(FL);

   if (r!=2) wblog(FL,"ERR %s() invalid operator rank-%d",FCT,r);
   if (CGR.dim2) {
      if (CGR.dim1!=QIDX.dim1) wblog(FL,
         "ERR %s() size mismatch (%dx%d <> %dx%d)",
         FCT,QIDX.dim1,qtype.len,CGR.dim1,CGR.dim2
      );
      for (unsigned j=0; j<CGR.dim2; ++j) {
         if (CGR(i,j).cgb) { d*=CGR(i,j).cgb->dim(); }
      }
   }

   if (!d) wblog(FL,"ERR %s() got dim(cgd)=%d !?",FCT,d); 

   return d;
};

template <class TQ, class TD>
void getQDimGen(
   wbvector< const QSpace<TQ,TD>* > A,
   wbMatrix<TQ> &Q, wbvector<widx_t> &S,
   wbindex I, 
   wbMatrix<widx_t> *SC_ 
){
   if (!cPVEC_nnz(FLF,A,'w')) wblog(FL,"ERR %s() got empty input",FCT);

   unsigned i,j,k,l,ic,d, d2=0, m=0, nq,p,i1,i2, n=0; widx_t *ip;
   unsigned r_=-1, rk=-1, QDIM=A[0]->QDIM;
   const QS_TYPES &otype=A[0]->otype;
   const QVec &qtype=A[0]->qtype;
   wbMatrix<widx_t> S0,SC,II;
   wbvector<widx_t> D;
   wbperm P;

   char cgflag=(SC_ ? 1 : 0), e=0;

   nq=A[0]->CGR.dim2;

   for (k=0; k<A.len; ++k) {
      const QSpace<TQ,TD> &Ak=(*A[k]);
      n+=Ak.QIDX.dim1;

      if (I.len) { r_=rk=Ak.rank(FL); }
      else {
         if (int(rk=Ak.isOperator(&r_,'x'))<=0) wblog(FL,
            "ERR %s() invalid rank-%d operator\n(%d: r=%d/%d, '%s')",
            FCT,Ak.itags.len,k+1,rk,r_,A[k]->otype2Str().data
         );
         if (Ak.otype!=otype) wblog(FL,
            "ERR %s() object type inconsistency (%s ; %s)",
            FCT, Ak.otype2Str().data, A[k]->otype2Str().data
         );
      }

      if (!k) d2=r_*QDIM;

      if (Ak.qtype!=qtype) wblog(FL,
         "ERR %s() Q-type inconsistency (%s ; %s)",
          FCT, Ak.qStr().data, STR(qtype));
      if (Ak.QDIM!=QDIM || Ak.QIDX.dim2<d2) wblog(FL,
         "ERR %s() rank inconsistency (%d: %d/%d %d/%d)",
          FCT, k+1, Ak.QDIM, QDIM, Ak.QIDX.dim2, d2);
      if (Ak.QIDX.dim1!=Ak.DATA.len)
          wblog(FL,"ERR %s() size mismatch (%dx%d <> %d)",
          FCT, Ak.QIDX.dim1, Ak.QIDX.dim2, Ak.DATA.len
      );

      if (Ak.gotCGS(FL)<=0) {
         if (!k) { cgflag=0; } else if (cgflag) { e|=1; }
      }
      else if (!k) { cgflag|=2; }
      else if (!(cgflag&2)) { e|=2; }

      if (cgflag && Ak.CGR.dim2!=nq) { e|=4; }
      if (e) { wblog(FL, 
         "ERR %s() CGR inconsistency (%d/%d: %d/%d)",
         FCT, k+1,A.len, Ak.CGR.dim2,nq);
      }
   }
   if (int(rk)<=0) wblog(FL,"ERR got empty QSpaces (r=%d)",rk);
   cgflag&=1;

   if (!I.len) I.Index(r_); 
   else if (I.anyGE(r_)) wblog(FL,
     "ERR %s() index out of bounds (%s; %d)",FCT,STR((I+1)),r_);

   S0.init(n,I.len);
   II.init(n,I.len*2);
   Q .init(n,I.len*QDIM);

   wbvector<unsigned> rm; 
   if (cgflag) {
      SC.init(n,nq*I.len).set(1);
      rm.init2val(nq,rk);
   }
   for (j=0; j<nq; ++j) {
      if (qtype[j].permitsOM(rk)) { m=1;
         if (cgflag) { ++rm[j]; } else break;
      }
   }

   for (l=k=0; k<A.len; ++k) {
      const QSpace<TQ,TD> &Ak=(*A[k]);
      const wbMatrix<TQ> &qq=Ak.QIDX; unsigned nd=Ak.DATA.len;

      for (i=0; i<nd; ++i, ++l) {
         const wbvector<widx_t> &s=Ak.DATA[i]->SIZE;
         if (!Ak.DATA[i]->isRankM(rk,m)) {
            wblog(FL,"ERR %s() rank mismatch (data[%d]: %s (r=%d)",
            FCT,i+1,SSTR_(Ak.DATA[i]), rk);
         }

         Q .recSetB(l,I,QDIM, qq.rec(i));
         S0.recSetB(l,I,1, s.data);

         if (cgflag) {
            for (j=0; j<nq; ++j) {
               wbvector<widx_t> s; Ak.CGR(i,j).getSize(s);
               if (!s.len || s.allEqual(1)) continue;
               if (s.len<rk || s.len>rm[j]) wblog(FL,
                  "ERR %s() rank mismatch (cgb(%d,%d): [%s](r=%d/%d)",
                   FCT,i+1,j+1, SSTR(Ak.CGR(i,j)), rk,r_
               );
               for (ic=j, i1=0; i1<I.len; ++i1, ic+=nq) {
                  SC(l,ic)=s[I.data[i1]];
               }
            }
         }

         ip=II.rec(l);  
         for (p=j=0; j<II.dim2; j+=2, ++p) { ip[j]=k; ip[j+1]=I[p]; }
      }
   }

   Q .Reshape(n*I.len,QDIM); 
   S0.Reshape(n*I.len,1);
   II.Reshape(n*I.len,2);

   if (cgflag)
   SC.Reshape(n*I.len,nq);

   Q.groupRecs(P,D); S.init(D.len);
   if (cgflag) { SC_->init(D.len,nq); }

   for (l=i=0; i<D.len; ++i) { d=D[i]; i1=P[l++];
      S[i]=S0(i1,0);
      if (cgflag) SC_->recSetP(i,SC.ref(i1));
      for (j=1; j<d; ++j, ++l) { i2=P[l];
         if (S0.recCompare(i1,i2)) {
            sprintf_str(
              "%s() data{} block size inconsistency (SC=%d)\n"
              "A(%ld) @%ld <> A(%ld) @%ld: Q=[%s] D = %s / %s !?",FCT,cgflag,
               II(i1,0)+1, II(i1,1)+1,
               II(i2,0)+1, II(i2,1)+1, Q.rec2Str(i).data, 
               S0.rec2Str(i1).data, 
               S0.rec2Str(i2).data  
            );
            if (II(i1,0)==II(i2,0)) wblog(FL,"ERR %s\n"
              "(specify explicit leg for non-hermitian objects?)",str);
            else wblog(FL,"ERR %s",str);
         }
         if (cgflag && SC.recCompare(i1,i2)) wblog(FL,
            "ERR %s() size inconsistency\nCGS: %d:%d <> %d:%d",
             FCT, II(i1,0)+1, II(i1,1)+1, II(i2,0)+1, II(i2,1)+1
         );
      }
   }
};

template <class TQ, class TD> inline
void QSpace<TQ,TD>::getQDim(wbMatrix<TQ> &Q,
   wbvector<widx_t> &S,
   wbMatrix<widx_t>*SC   
) const {

   wbindex I;
   wbvector< const QSpace* > A(1); A[0]=this;

   getQDimGen(A,Q,S,I,SC); 
};

template <class TQ, class TD> inline
void QSpace<TQ,TD>::getQDim(unsigned k,
   wbMatrix<TQ> &Q,
   wbvector<widx_t> &SD,
   wbMatrix<widx_t> *SC 
) const {

   unsigned i,j,q, r=rank(FL), r_=r, n=DATA.len;
   wbvector<widx_t> D;
   wbindex Ig;
   wbperm P;

   wbMatrix<TQ> QX;

   if (k>=r) wblog(FL,"ERR %s() dim out of bounds (%d/%d)",FCT,k+1,r);
   if (QIDX.dim1!=n) wblog(FL,
      "ERR %s() QSpace inconsistency (%d/%d)",FCT,QIDX.dim1,n);

   getQsub(k,Q).groupRecs(P,D,-1,1,&Ig,&QX);
   SD.init(D.len).set(-1);

   if (r>2 && permitsOM()) { ++r_; }

   for (i=0; i<n; ++i) {
      const wbvector<widx_t> &si=DATA[i]->SIZE;
      widx_t &s=SD.at(Ig[i]);
      if (k>=si.len || si.len>r_) wblog(FL,
         "ERR %s() index out of bounds (%d/%d,%d)",FCT,k+1,si.len,r);
      if (int(s)>=0) {
         if (s!=si.data[k]) wblog(FL,
            "ERR %s() size mismatch DATA(%d): %s [S(%d)=%d]",
            FCT,i+1,DATA[i]->sizeStr().data, k+1, s
         );
      }
      else s=si.data[k];
   }

   if (SC==nullptr) return;
   if (gotCGS(FL)<=0) { SC->init(); return; }

   SC->init(D.len,CGR.dim2).set(-1);

   for (i=0; i<n; ++i) {
      for (j=0; j<CGR.dim2; ++j) {
         widx_t &c=SC->at(Ig[i],j);
         q=CGR(i,j).qdim(k,&qtype[j]);

         if (int(c)>=0) {
            if (c!=q) wblog(FL,
               "ERR %s() size mismatch CGR(%d,%d): %s [S(%d)=%d]",
               FCT,i+1,j+1, CGR(i,j).sizeStr().data, k+1, q
            );
         } else c=q;
      }
   }
};

template <class TQ, class TD>
char QSpace<TQ,TD>::getQDir(unsigned i) const {
   if (i>=itags.len) wblog(FL,
      "ERR %s() index out of bounds (i=%d/%d)",FCT,i,itags.len);
   return (itags[i].isConj() ? -1 : 1);
};

template <class TQ, class TD>
QDir& QSpace<TQ,TD>::getQDir(QDir &qdir) const {

   qdir.init(itags);

   if (CGR.dim1) {
      for (unsigned j=0; j<CGR.dim2; ++j) {
         if (!CGR(0,j).sameQDir(FL,qdir)) wblog(FL,
            "ERR %s() got qdir inconsistency\nj=%d: %s <> %s",
            FCT,j+1,STR(CGR(0,j)), STR(qdir)
         );
      }
   }
   return qdir;
};

template <class TQ, class TD>
bool QSpace<TQ,TD>::sameQDir(
   const char *F, int L, const QSpace<TQ,TD> &B
 ) const {

   if (qtype!=B.qtype) wblog(F_L,
      "ERR %s() qtype mismatch (%s <> %s)", FCT,qStr().data,B.qStr().data);
   if (!QIDX && !B.QIDX) return 1; 

   if (!  CGR &&   qtype.allAbelian()) { return 1; } 
   if (!B.CGR && B.qtype.allAbelian()) { return 1; }

   if (qtype.len!=CGR.dim2 || qtype.len!=B.CGR.dim2) wblog(F_L,
      "ERR %s() qtype mismatch (%s <> %s %d/%d / %d/%d)", FCT,
      qStr().data,B.qStr().data, qtype.len, B.qtype.len, CGR.dim2, B.CGR.dim2
   );

   for (unsigned i=0, n=qtype.len; i<n; ++i) {
      if (!CGR(0,i).sameQDir(0,0,B.CGR(0,i))) {
         if (F) wblog(F_L,"ERR %s() got qdir mismatch\n"
            "(%s => %s <> %s => %s)", FCT,
              CGR(0,i).qdir2Str('V').data,   CGR(0,i).qdir2Str().data,
            B.CGR(0,i).qdir2Str('V').data, B.CGR(0,i).qdir2Str().data
         );
         return 0;
      }
   }

   return 1;
};

template <class TQ, class TD>
char QSpace<TQ,TD>::checkQ(const char *F, int L) const {

   if (!QDIM) {
      if (QIDX.dim2 || itags.len) {
         if (F) wblog(F,L,"ERR %s() QDIM not yet set "
            "(%dx%d,%d) !?",FCT, QIDX.dim1, QIDX.dim2, itags.len);
         else return 1;
      }
      return 0;
   }
   if (QIDX.dim2%QDIM) {
      if (F) wblog(F,L,
         "ERR %s() invalid QDIM (r=%d/%d)",FCT,QIDX.dim2,QDIM);
      else return 2;
   }

   if (qtype.allAbelian()) {
      if (CGR.dim1 || CGR.dim2) {
         if (F) wblog(F,L,"ERR %s() "
            "got CGR data (%dx%d) without qtype",FCT,CGR.dim1,CGR.dim2);
         else return 5;
      }
   }
   else {
      if (QDIM!=qtype.Qlen()) {
         if (F) wblog(F,L,
            "ERR %s() QDIM inconsistency (%d/%d; %s)",
            FCT, qtype.Qlen(), QDIM, qStr().data);
         else return 3;
      }
      if (CGR.dim2!=qtype.len || CGR.dim1!=QIDX.dim1) {
         if (F) wblog(F,L,
            "ERR %s() CGR inconsistency: %dx%d <> %dx%d",
            FCT, CGR.dim1, CGR.dim2, QIDX.dim1, qtype.len);
         else return 4;
      }
   }

   if (itags.len != QIDX.dim2/QDIM) {
      if (!itags.len) {
         if (F) wblog(F,L,"ERR %s() got empty itags !?",FCT);
         else return 6;
      }
      else {
         if (F) wblog(F,L,"ERR %s() got invalid itags "
            "(len=%d != (%d/%d)) !?",FCT,itags.len,QIDX.dim2,QDIM);
         else return 7;
      }
   }

   return 0;
};

template <class TQ, class TD>
template <class TB>
void QSpace<TQ,TD>::checkQ(
   const char *F, int L, const QSpace<TQ,TB> &B, char aflag
 ) const {

   if (QDIM!=B.QDIM) wblog(F_L,
      "ERR %s() QDIM inconsistency (%d/%d)",FCT,QDIM,B.QDIM);

   if (qtype.len && QDIM!=qtype.Qlen()) wblog(F_L,
      "ERR %s() length inconsistency in qtype (%d/%d)",
       FCT, qtype.Qlen(), QDIM);

   if (qtype!=B.qtype) {
     if (aflag && allAbelian()==1 && B.allAbelian()==1 &&
        (!qtype.len || !B.qtype.len)) { aflag='A'; } 
     if (aflag!='A') wblog(F_L, 
        "ERR %s() qtype inconsistency ('%s' <> '%s')",
         FCT,qStr().data,B.qStr().data
     );
   }

   if (aflag!='A') { 
      if (CGR.dim2!=B.CGR.dim2) wblog(F_L,
         "ERR %s() CGR inconsistency: %dx%d <> %dx%d",
         FCT, CGR.dim1, CGR.dim2, B.CGR.dim1, B.CGR.dim2
      );
   }
};

template <class TQ, class TD> inline
void QSpace<TQ,TD>::recSave2(unsigned k, unsigned i) { 

   if (DATA.len!=QIDX.dim1 || (CGR.dim2 && CGR.dim1!=QIDX.dim1)) {
      wblog(FL,"ERR QSpace inconsistency (%d/%d/%d)",
      CGR.dim1, DATA.len, CGR.dim1);
   }

   if (i==k) return;
   if (i>=QIDX.dim1 || k>=QIDX.dim1) wblog(FL,
      "ERR index out of bounds ({%d,%d}/%d)",i,k,QIDX.dim1);

   QIDX.recSet(k,i);

   WB_DELETE_1(DATA[k]);
   DATA[k]=DATA[i]; DATA[i]=nullptr;

   if (CGR.dim2) {
   for (unsigned j=0; j<CGR.dim2; ++j) {
      CGR(k,j)=CGR(i,j); CGR(i,j).init();
   }}
};

template <class TQ, class TD>
int QSpace<TQ,TD>::getQOverlapU(
   const char *F, int L,
   const QSpace<TQ,TD> &B,
   wbMatrix<TQ> &QQ,
   wbMatrix<int> &IQ,
   char olonly 
) const{

   unsigned i,i1,i2,d,k,l, m=QIDX.dim1; 
   wbvector<widx_t> D;
   wbperm P;

   if (!B.QIDX.dim1) { QQ=  QIDX; } else
   if (!  QIDX.dim1) { QQ=B.QIDX; }
   else {
      if (!QIDX.dim2 || QIDX.dim2!=B.QIDX.dim2 || rank(FL)!=B.rank(FL)) {
         if (F) wblog(F,L,"ERR QSpace mismatch (QDIM=%d/%d; r=%d/%d)",
            QDIM, B.QDIM, rank(FL), B.rank(FL));
         return 1;
      }
      QQ.Cat(1,QIDX,B.QIDX);
   }

   QQ.groupRecs(P,D);
   IQ.init(QQ.dim1,2);

   for (k=l=i=0; i<QQ.dim1; ++i, l+=d) { d=D[i];
      if (k<i) QQ.recSet(k,i);
      if (d==1) { if (!olonly) {
         i1=P[l]; if (i1<m)
              { IQ(k,0)=i1; IQ(k,1)=-1;   }
         else { IQ(k,0)=-1; IQ(k,1)=i1-m; }
         ++k;
      }}
      else if (d==2) { i1=P[l]; i2=P[l+1];
         if (i1<m && i2>=m) { i2-=m; } else
         if (i2<m && i1>=m) { i1-=m; SWAP(i1,i2); }
         else {
            if (F) wblog(F,L,
               "ERR %s() non-unique QIDX (%d,%d/%d)",FCT,i1,i2,m);
            return 2;
         }
         IQ(k,0)=i1; IQ(k,1)=i2;
         ++k;
      }
      else {
         if (F) wblog(F,L,
            "ERR %s() got non-unique QIDX (%d)",FCT,d);
         return 3;
      }
   }

   if (k<i) {
      QQ.Resize(k,QQ.dim2);
      IQ.Resize(k,IQ.dim2);
   }

   return 0;
};

template <class TQ, class TD>
int QSpace<TQ,TD>::getQOverlap(
   const char *F, int L, const QSpace<TQ,TD> &B,
   wbMatrix<TQ> &QQ, wbMatrix<int> &IQ,
   wbvector<widx_t> &DQ 
) const{

   if (!QIDX.dim1 || !B.QIDX.dim1) {
      QQ.init(); IQ.init(); DQ.init();
      if (F) wblog(F,L,"WRN got empty QSpace");
      return -1;
   }
   if (!QIDX.dim2 || QIDX.dim2!=B.QIDX.dim2 || rank(FL)!=B.rank(FL)) {
      if (F) wblog(F,L,"ERR QSpace mismatch (QDIM=%d/%d; r=%d/%d)",
         QDIM, B.QDIM, rank(FL), B.rank(FL));
      return 1;
   }

   unsigned i,l;
   wbindex I1,I2; wbperm P1,P2;
   wbMatrix<TQ> Q1(QIDX), Q2(B.QIDX);

   Q1.SortRecs(P1);
   Q2.SortRecs(P2); matchSortedIdx(Q1,Q2,I1,I2,DQ);

   IQ.init(I1.len,2);
   for (i=0; i<I1.len; ++i) {
      IQ(i,0)=P1[I1[i]];
      IQ(i,1)=P2[I2[i]];
   }

   QQ.init(DQ.len,Q1.dim2);
   for (l=i=0; i<DQ.len; ++i) { QQ.recSetP(i,Q1.rec(I1[l])); l+=DQ[i]; }

   return 0;
};

template <class TQ, class TD>
bool QSpace<TQ,TD>::findDimQ(const TQ* q0, unsigned &d) const {

   unsigned i,j,n, r=-1, s=QDIM*sizeof(TQ);
   const TQ *q=QIDX.data;

   if (!isConsistent(r)) wbdie(FL,str);
   if (QIDX.isEmpty()) return 0;

   r=rank(); 
   n=QIDX.dim1*r;

   for (i=0; i<n; ++i, q+=QDIM) { if (!memcmp(q,q0,s)) break; }
   if (i==n) {
      wbvector<TQ> x(QDIM,q0);
      wblog(FL,"TST %s() %d/%d (%d)",FCT,i,n,QDIM);
      MXPut(FL,"i").add(*this,"A").add(x,"q");
      return 0;
   }

   j=i%r; i=i/r; d=DATA[i]->SIZE[j];
   return 1;
};

template <class TQ, class TD>
bool QSpace<TQ,TD>::findDimQ(const TQ* q, unsigned k, unsigned &d) const {

   unsigned i, n=QIDX.dim1, r=-1, s=QDIM*sizeof(TQ);
   TQ *q0=QIDX.data+k*QDIM;

   if (!isConsistent(r)) wbdie(FL,str);
   if (QIDX.isEmpty()) return 0;

   for (i=0; i<n; ++i, q0+=QIDX.dim2) {
      if (!memcmp(q0,q,s)) break;
   }  if (i==n) return 0;

   d=DATA[i]->SIZE[k];
   return 1;
};

template <class TQ, class TD> inline
bool QSpace<TQ,TD>::isEmpty() const {

   if (QIDX.dim1!=DATA.len) { 
      if (QIDX.dim1 || DATA.len>1) wblog(FL, 
         "ERR QSpace inconsistency (%d,%d)",QIDX.dim1, DATA.len);
      return (!DATA.len); 
   }
   if (!QIDX.dim2) { 
      if (DATA.len>1) wblog(FL,
         "ERR QSpace inconsistency (%d,%d)",QIDX.dim1, DATA.len);
      return (!DATA.len);
   }

   if (itags.len && QIDX.data) {
      if (QIDX.dim2!=itags.len*QDIM) { wblog(FL,
         "ERR %s() QSpace inconsistency (QIDX %s / %dx%d)",
         FCT,SSTR(QIDX),itags.len,QDIM);
      }
   }

   if (CGR.data && (CGR.dim1!=QIDX.dim1 || CGR.dim2!=qtype.len)) {
      wblog(FL,"ERR QSpace size mismatch %s / %dx%d ",
      SSTR(CGR),QIDX.dim1,qtype.len);
   }

   return (!DATA.len); 
};

template <class TQ, class TD> inline
bool QSpace<TQ,TD>::isConsistent(
    const char *F, int L, unsigned &r0, char level) const {

    if (isEmpty()) {
       if (int(r0)<0) { r0=0; } 
       return 1;
    };

    unsigned r=rank(F,L); str[0]=0;

    if (int(r)<0) wblog(FL,"ERR %s() r=%d",FCT,r); 

    if (int(r0)<0) { r0=r; } else
    if (r0!=r && (r0!=2 || r!=1)) { sprintf_str(
       "rank mismatch (is %d, should be %d)",r,r0);
       if (F) wblog(F,L,"ERR %s",str);
       return 0;
    }

    if (itags.len && itags.len!=r) {
       sprintf_str("QSpace::itags inconsistency (%ldx%ld/%d, %ld/%d)",
          QIDX.dim1, QIDX.dim2, QDIM, itags.len,r);
       if (F) wblog(F,L,"ERR %s",str);
       return 0;
    }

    if (QIDX.dim1!=DATA.len) wblog(F_L, 
       "ERR %s() QSpace size mismatch %s / %d [x%d]",
       SSTR(QIDX), DATA.len, QDIM);
    if (level&1 && !QIDX.isUnique()) wblog(FL, 
       "ERR %s() QSpace got non-unique rows in [Q{:}]",FCT);

    gotCGS(F_L);

    if (isref || level<=1) { return 1; } 

    unsigned i,l,m=0;
    for (i=0; i<qtype.len; ++i) { if (qtype[i].permitsOM(r)) { ++m; }}

    for (i=0; i<DATA.len; ++i) {
       if (!DATA[i]) { sprintf_str("QSpace got null space DATA[%d]",i+1);
           if (F) wblog(F,L,"ERR %s",str);
           return 0;
       }

       if (!DATA[i]->isRankM(r,m)) {
          if ((l=DATA[i]->SIZE.len)==1 && !r) { continue; }

          sprintf_str("QSpace inconsistency (%s)\nDATA[%d] has rank "
             "%d / %d (%s)",SHORT_FL, i+1, l, r0, SSTR_(DATA[i]));
          if (F) wblog(F,L,"ERR %s",str);
          return 0;
       }
       if (r0>=2) { continue; }

       if ((r0==1 && (DATA[i]->SIZE[0]!=1 && DATA[i]->SIZE[1]!=1)) ||
           (r0==0 && (DATA[i]->SIZE[0]!=1 || DATA[i]->SIZE[1]!=1)) ) {
           sprintf_str("invalid size %s for rank-%d object",
           DATA[i]->sizeStr().data, r0);
           if (F) wblog(F,L,"ERR %s",str);
           return 0;
       }

       if (r0==0 && i) { sprintf_str(
          "rank-%d object can only have max 1 element (%d)", r0, i+1);
           if (F) wblog(F,L,"ERR %s",str);
           return 0;
       }
    }

    return 1;
};

template <class TQ, class TD>
int QSpace<TQ,TD>::isOperator(unsigned *r_, char xflag) const {

   unsigned rk=rank(FL), r=rk; 

   if (r_ && int(*r_)<0) { *r_=-rk; }

   if (otype==QS_OPERATOR) { char c=gotCGS(FL);
      if ((c>0 && r!=3) || (c<=0 && (r<2 || r>3))) { wblog(FL,
         "ERR %s() got rank-%d for %sabelian %s", FCT,r, c? "non-":"",
         otype2Str().data); }
      r=2;
   }
   else if (otype!=QS_NONE) { return -9; }

   if (int(rk)<2) {
      return -1;
   }

   if (r>2) { 
      if (r>3 && xflag) {
         if (rk%2) { return -3; }
         if (!itags.isOpX()) { return -4; }
      }
      else {
         if (!itags.isOp()) { return -2; }
         if (rk>3) { return 0; }
         r=2;
      }
   }

   if (r_) {
      if (int(*r_)>0)
           { if ((*r_)!=r) { return 0; }}
      else { (*r_)=r; }
   }

   return rk;
};

template <class TQ, class TD> 
int QSpace<TQ,TD>::checkAbelianOp() const {

   if (itags.len && !itags.isOp()) { return -1; }
   if (otype!=QS_NONE && otype!=QS_OPERATOR) { return -2; }
   if (gotCGS() || !CGR.isEmpty()) { return -3; }

   unsigned r=rank(FL);

   return (r==2 || r==3 ? 1 : 0);
};

template<class TQ, class TD> 
unsigned QSpace<TQ,TD>::Reduce2AbelianOp(const char *F, int L) {

   unsigned r=rank(F_L), i=0;

   if (r<=2) { return r; }
   if (r>3 || !isAbelianOp()) { info("",1); wblog(F_L,
      "ERR %s() got invalid operator (r=%d)",FCT,r); }

   if (!CGR.isEmpty() || QIDX.dim2!=3*QDIM) { info("",1); wblog(F_L,
      "ERR %s() unexpected / invalid operator (r=%d) !?",FCT,r); }
   if ((r=getDim(2))!=1) { info("",1); wblog(F_L,
      "ERR %s() got non-scalar abelian operator (d=%d) !?",FCT,r); }

   QIDX.Resize(QIDX.dim1,QIDX.dim2-QDIM);

   for (; i<DATA.len; ++i) { DATA[i]->skipSingletons(F_L,2); }

   if (otype==QS_OPERATOR) { otype=QS_NONE; }

   if (itags.len) {
      if (itags.len==3) { itags.SkipLast(); }
      else wblog(FL,"ERR %s() unexpected itags=%s !?",FCT,STR(itags));
   }

   return --r;
};

template <class TQ, class TD> inline
bool QSpace<TQ,TD>::isDiagBlock(unsigned k) const {

    unsigned s=QDIM*sizeof(TQ);
    const TQ *q=QIDX.rec(k);

    if (!isConsistent_r(2)) wbdie(FL,str);
    if (k>=QIDX.dim1) wblog(FL,
       "ERR index out of bounds (%d/%d)",k,QIDX.dim1);

    if (memcmp(q, q+QDIM, s)) return 0;
    if (DATA[k]->SIZE[0]!=DATA[k]->SIZE[1]) return 0;

    return 1;
};

template <class TQ, class TD> inline
bool QSpace<TQ,TD>::isBlockDiagMatrix(
    const char *F, int L, char dflag) const {

    unsigned r=rank(F_L); if (r%2) { return 0; }
    unsigned i, n=QIDX.dim1, r2=r/2, m=r2*QDIM,
       sq=m*sizeof(TQ), s2=r2*sizeof(TQ);
    TQ *q=QIDX.data;

    if (dflag && r>2) { wblog(FL,
       "WRN %s() dflag=%d will be ignored (r=%d)",FCT,dflag,r);
       dflag=0;
    }

    for (i=0; i<n; ++i, q+=QIDX.dim2) {
       if (memcmp(q, q+m, sq)) {
          if (F) wblog(F_L,"--> %s(): non-matching Q(%d,:)",FCT,i+1);
          return 0;
       }

       const wbvector<unsigned> &S=DATA[i]->SIZE;

       if (dflag) { 
          if (S.len!=2 || (S[0]!=1 && S[1]!=1)){
             if (F) wblog(F_L,"--> %s(): data{%d} got size %s [dflag]",
                FCT, i+1, SSTR_(DATA[i]),dflag);
             return 0;
          }
       }
       else {
          if (S.len<r || S.len>r+1) wblog(FL,"ERR %s() got "
             "rank mismatch (%s /%d)",FCT,STR(DATA[i]->SIZE),r);
          if (memcmp(S.data, S.data+r2, s2)!=0) {
             if (F) wblog(F_L,"--> %s(): data{%d} got size %s",
                FCT, i+1, SSTR_(DATA[i]));
             return 0;
          }
       }
    }

    return 1;
}

template <class TQ, class TD> inline
bool QSpace<TQ,TD>::isDiagMatrix(const double eps) const {

    if (!isBlockDiagMatrix()) return 0;

    for (unsigned i=0; i<DATA.len; ++i) {
       if (!(DATA[i]->isDiagMatrix(eps))) return 0;
    }

    if (gotCGS(FL)>0) {
       for (unsigned  n=CGR.numel(), i=0; i<n; ++i) {
          if (!CGR.data[i].isDiagCSC(eps)) return 0;
       }
    }

    return 1;
};

template <class TQ, class TD> inline
bool QSpace<TQ,TD>::isIdentityMatrix(const TD eps) const {

   if (!isBlockDiagMatrix()) return 0; 

   if (CGR.isEmpty() || qtype.allAbelian()) {
      for (unsigned i=0; i<DATA.len; ++i)
      if (!(DATA[i]->isIdentityMatrix(eps))) return 0;
   }
   else {
      if (CGR.dim1!=DATA.len) wblog(FL,"ERR %s() "
         "size mismatch (%d / %dx%d)",FCT,DATA.len,CGR.dim1,CGR.dim2);
      TD dval; wbvector<double> cfac; unsigned i=0, j=0, r=rank(FL);

      for (; i<DATA.len; ++i) { cfac.init();
         for (j=0; j<CGR.dim2; ++j) {
            if (!CGR(i,j).isIdentityCG(&cfac)) { return 0; }
         }
         if (cfac.norm()<1e-12) { return 0; }

         if (cfac.len==1 && r==DATA[i]->SIZE.len) {
            dval=TD(1)/cfac[0];
            if (!(DATA[i]->isProptoId(dval,eps))) { return 0; }
         }
         else {
            IterOM_DATA<TD> I(*this,i);
            if (I.niter!=cfac.len) wblog(FL,"WRN %s() IterOM_DATA "
               "length mismatch (%d/%d)",FCT,I.niter,cfac.len);

wblog(FL,"ERR %s() check IterOM_DATA",FCT);

            for (I.begin(); !I.end(); ++I) {
               dval=TD(1)/cfac[I.iter];
               if (!(I.DATA().isProptoId(dval,eps))) { return 0; }
            }
         }
      }
   }

   return 1;
};

template <class TQ, class TD> inline
char QSpace<TQ,TD>::hasIdentityCGS(const TD eps) const {

    for (unsigned n=CGR.numel(), i=0; i<n; ++i) {
       if (!CGR.data[i].isIdentityCG(nullptr,eps)) {
          return 0;
       }
    }

    return 1;
};

template <class TQ, class TD> inline
bool QSpace<TQ,TD>::isHConj(
   const char *F, int L, double eps, char vflag) const {

   bool cgflag=(gotCGS(FL)>0);

   if (cgflag) {
      unsigned r=rank(F,L);

      if (r<2 || (r!=3 && !isQSym(F,L)) || (r==3 && getDIM(2)>1)) {
         if (vflag || F) {
            if (r<2 || r==3) { sprintf_str("%s got rank-%d QSpace (%s)",
               SHORT_FL,r, sizeStr('v').data); }
            else { sprintf_str(
               "%s got non-symmetry Q data (r=%d)", SHORT_FL,r); }
            if (F) wblog(F_L,"ERR %s() %s",FCT,str);
         }
         return 0;
      }

      QSpace<TQ,TD> X; wbperm P;
      double e=0;

      if (r==3)  
           { P.initStr(FL,"2,1,3"); }
      else { P.initTranspose(r); } 
      permute(X,P).Conj();

      if (r==3) { if (X.ConjOpScalar()) return 0; } else
      if (!itags.isOp(r,'L')) { return 0; } 
      X.itags.tSet(itags); 

      try { X-=(*this); e=X.norm(); }
      catch (...) {
        if (r!=3) wblog(FL,
           "ERR %s() failed to subtract block-transpose QSpace\n"
           "ERR %s (p=[%s]; r=%d)",FCT,sizeStr('v').data,STR(P),r);
        else wblog(FL,
           "ERR %s() failed to subtract transpose IROP\n"
           "ERR %s (p=[%s]; r=%d)",FCT,sizeStr('v').data,STR(P),r
        );
      }
      if (e>eps) { double q=norm(); if (q>1) eps*=q;
         if (e>eps) {
            if (vflag || F) {
               sprintf_str("%s got non-symmetric data (%g)",SHORT_FL,e);
               if (F) wblog(F,L,"ERR %s() %s",FCT);
            }
            return 0;
         }
      }
      return 1;
   }
   else {
      return isSym_aux(F,L,FCT,eps,'s',vflag);
   }
};

template <class TQ, class TD> inline
bool QSpace<TQ,TD>::isAHerm(
   const char *F, int L, double eps, char vflag) const {

   bool cgflag=(gotCGS(FL)>0);

   if (cgflag) {
      if (!isQSym(F,L)) {
         if (vflag) wblog(FL,"ERR %s() got size mismatch",FCT);
         return 0;
      }
      unsigned r=rank(FL); QSpace<TQ,TD> X;
      wbperm P; P.initTranspose(r); permute(X,P); X+=(*this);
      double e=X.norm();
      if (e>eps) { if (vflag) wblog(FL,
         "ERR %s() got non-symmetric data (%g)",FCT,e);
         return 0;
      }
      else return 1;
   }
   else {
      return isSym_aux(F,L,FCT,eps,'s',vflag);
   }
};

template <class TQ, class TD> inline
bool QSpace<TQ,TD>::isSym_aux(
  const char *F, int L, const char *fct,
  RTD eps, char symflag, char vflag
) const {

    wbMatrix<TQ> Q1,Q2;
    wbperm P1,P2,iP2,P;

    unsigned i,j, r=-1, k=QDIM, K; int e=0;

    double deps=this->norm(); 
       if (deps<1) { deps=1; }; deps*=double(eps);

    gotCGS(FL);

    if (isEmpty()) { return 1; }
    str[0]=0; isConsistent(FL,r);

    if (r%2) {
       sprintf_str("%s() requires even-rank (%d)",fct,r);
       if (F) wblog(F,L,"ERR %s",str);
       return 0;
    }

    K=(unsigned)r/2;

    Q1=QIDX; Q1.SortRecs(P1); 

    P.init((unsigned)r).Rotate(K);
    QIDX.blockPermute(Q2,P); Q2.SortRecs(P2);

    if (Q1!=Q2) {
       sprintf_str("QIDX[:,1] does not match QIDX[:,2]");
       if (F || vflag) wblog(F_L,"ERR %s",str);
       return 0;
    }
    if (!Q1.isUniqueSorted()) {
       sprintf_str("WRN QIDX is not unique!");
       if (F || vflag) wblog(F_L,"ERR %s",str);
       return 0;
    }

    P2.invert(iP2);

    for (i=0; i<DATA.len; ++i) {
       j=P1[iP2[i]]; if (j<i) continue; 

       const wbarray<TD> &di=(*DATA[i]);
       const wbarray<TD> &dj=(*DATA[j]);

       if (symflag=='s') {
          if (i==j)
               { if (!di.isHConj(   deps)) e=1; }
          else { if (!di.isHConj(dj,deps)) e=2; }
       }
       else if (symflag=='a') {
          if (i==j)
               { if (!di.isAHerm(   deps)) e=3; }
          else { if (!di.isAHerm(dj,deps)) e=4; }
       }
       else wblog(F,L,"ERR %s() invalid flag %c<%d>",FCT,symflag,symflag);

       for (k=0; !e && k<CGR.dim2; ++k) {
          if (CGR(i,k).cgb && !CGR(i,k).cgb->cgd.isScalar()) {
          const cdata__ &ci=(CGR(i,k).cgb->cgd);
          const cdata__ &cj=(CGR(j,k).cgb->cgd);

          if (symflag=='s') {
             if (i==j)
                  { if (!ci.isHConj(   eps)) { e=11; break; }}
             else { if (!ci.isHConj(cj,eps)) { e=12; break; }}
          }
          else {
             if (i==j)
                  { if (!ci.isAHerm(   eps)) { e=13; break; }}
             else { if (!ci.isAHerm(cj,eps)) { e=14; break; }}
          }
       }}

       if (!e && !CGR.isEmpty()) { 
          cdata__ ci,cj; char gotC=0;

          for (k=0; k<CGR.dim2; ++k) { if (CGR(i,k).cgb) {
             if (gotC) {
                ci.Kron(CGR(i,k).cgb->cgd); if (i!=j) {
                cj.Kron(CGR(j,k).cgb->cgd); }
             }
             else { gotC=1;
                ci=(CGR(i,k).cgb->cgd); if (i!=j) {
                cj=(CGR(j,k).cgb->cgd); };
             }
          }}

          if (symflag=='s') {
              if (i==j)
                   { if (!ci.isHConj(   eps)) { e=21; break; }}
              else { if (!ci.isHConj(cj,eps)) { e=22; break; }}
          }
          else {
              if (i==j)
                   { if (!ci.isAHerm(   eps)) { e=23; break; }}
              else { if (!ci.isAHerm(cj,eps)) { e=24; break; }}
          }
       }

       if (e) { wbvec<char> istr(128);
          if (e<10) { istr.catf(FL,
             "%s %s() non-matching data{%d} <> data{%d} (e=%d)",
              SHORT_FL,fct, i+1, j+1,e);
          }
          else { istr.catf(FL,
             "%s %s() non-matching cgs{%d,%d} <> cgs{%d,%d} (e=%d)",
              SHORT_FL,fct, i+1, k+1, j+1, k+1,e);
          }

          if (str[0]) { strcat(str,"\n"); }
          strcat(str,istr.data);

          if (vflag) { 
             const char *s=strstr(istr.data,"non-match");
             if (s) wblog(F_L,"%s ",s);
          }
          else if (vflag=='!') { MXPut(FL,"qs")
             .add(*this,"A").add(i+1,"i").add(j+1,"j")
             .add(e<10 ? -1 : k+1,"k").add(e,"e");
          }

          return 0; 
       }
    }

    return 1; 
};

template <class TQ, class TD> inline
bool QSpace<TQ,TD>::isQSym(const char *F, int L, char dflag) const {

    wbMatrix<TQ> Q1,Q2;
    wbperm P1,P2,iP2,P;

    unsigned i,j,l,r2,r=-1, k=QDIM; int e=0;
    gotCGS(FL);

    if (isEmpty()) { return 1; }
    str[0]=0; isConsistent(FL,r);

    if (r%2) {
       sprintf_str("%s() requires even-rank (%d)",FCT,r);
       if (F) wblog(F,L,"ERR %s",str);
       return 0;
    }

    r2=(unsigned)r/2;

    Q1=QIDX; Q1.SortRecs(P1); 

    P.init((unsigned)r).Rotate(r2);
    QIDX.blockPermute(Q2,P); Q2.SortRecs(P2);

    if (Q1!=Q2) {
       sprintf_str("QIDX[:,1] does not match QIDX[:,2]");
       if (F) wblog(F,L,"ERR %s",str);
       return 0;
    }
    if (!Q1.isUniqueSorted()) {
       sprintf_str("WRN QIDX is not unique!");
       if (F) wblog(F,L,"ERR %s",str);
       return 0;
    }

    P2.invert(iP2);

    if (dflag && !strchr("rc",dflag)) wblog(FL,
       "ERR %s() invalid dflag=%s",FCT,cSTR(dflag));
    if (dflag && r!=2) wblog(FL,
       "ERR %s() expecting rank-2 for dflag=%s (%d)",FCT,cSTR(dflag),r);

    for (i=0; i<DATA.len; ++i) {
       j=P1[iP2[i]]; if (j<i) continue; 

       const wbvector<unsigned> &si=DATA[i]->SIZE, &sj=DATA[j]->SIZE;

       if (!dflag) {
          for (l=0; l<r2; ++l)
          if (si[l]!=sj[l+r2] || si[l+r2]!=sj[l]) {
              sprintf_str("data skew dimensional (%d,%d: %s vs. %s)",
                 i+1,j+1, SSTR(si), SSTR(sj));
              return 0;
          }
       }
       else
       if (i!=j || (dflag=='r' && si[0]!=1) || (dflag!='r' && si[1]!=1)) {
           sprintf_str("data expected to be diagonal (%s) "
           "\n(%d,%d: %s vs. %s)",cSTR(dflag),i+1,j+1,STR(si),STR(sj));
           return 0;
       }

       for (k=0; !e && k<CGR.dim2; ++k) { 

          if (!CGR(i,k).cgb || !CGR(j,k).cgb) {
             if (CGR(i,k).cgb || CGR(j,k).cgb || !qtype[k].isAbelian()) {
             wblog(FL,"ERR %s() (%d,%d; %d): 0x%lX 0x%lX (%s) !?",
                FCT,i,j,k, CGR(i,k).cgb, CGR(j,k).cgb, STR(qtype[k])
             ); }
             continue;
          }

          wbvector<unsigned> sa, sb;
             CGR(i,k).getSize(sa,'b');
             CGR(j,k).getSize(sb,'b');

          if (sa.len!=sb.len || (sa.len && sa.len!=r)) { wblog(FL,
             "ERR %s(%d,%d;%d) CGC size mismatch: %s <> %s (%d)\n"
             "%s\n%s", FCT,i+1,j+1,k+1, SSTR(sa), SSTR(sb), r,
             STR(CGR(i,k)), STR(CGR(j,k)));
          }
          if (!sa.len) continue;

          for (l=0; l<r2; ++l) {
          if (sa[l]!=sb[l+r2] || sa[l+r2]!=sb[l]) {
             sprintf_str("data skew dimensional (%d,%d; %d): %s / %s",
                i+1,j+1,k+1, SSTR(sa), SSTR(sb));
             if (F) wblog(F,L,"ERR %s",FCT,str);
             return 0;
          }}
       }
    }

    return 1;
};

template <class TQ, class TD>
template <class T2>
bool QSpace<TQ,TD>::hasSameQ(const QSpace<TQ,T2> &B) const {

    const QSpace<TQ,T2> &A=(*this);

    unsigned i,j,d,k,l,m, r=-1;
    wbperm pA, pB;
    wbvector<unsigned> dA, dB;
    wbMatrix<TQ> QAk, QBk;
    wbindex Ia, Ib;

    if (isEmpty() || B.isEmpty()) return 1;
    if (QIDX.dim2!=B.QIDX.dim2 || QDIM!=B.QDIM) return 0;

    if (!A.isConsistent(r)) { A.info("A"); wbdie(FL,str); }
    if (!B.isConsistent(r)) { B.info("B"); wbdie(FL,str); }

    for (k=0; k<r; ++k) {
       A.getQsub(k, QAk); QAk.groupRecs(pA,dA);
       B.getQsub(k, QBk); QBk.groupRecs(pB,dB);

       for (l=i=0; i<dA.len; ++i, l+=d) {
          d=dA[i]; m=dA[i]=DATA[pA[l]]->SIZE[k];

          for (j=1; j<d; ++j)  
          if (m!=(DATA[pA[l+j]]->SIZE[k])) { 
          wblog(FL,"ERR QSpace inconsistency (%d,%d) !?",
          m,DATA[pA[l+j]]->SIZE[k]); }
       }

       for (l=i=0; i<dB.len; ++i, l+=d) {
          d=dB[i]; m=dB[i]=B.DATA[pB[l]]->SIZE[k];

          for (j=1; j<d; ++j)  
          if (m!=(B.DATA[pB[l+j]]->SIZE[k])) { 
          wblog(FL,"ERR QSpace inconsistency (%d,%d) !?",
          m, B.DATA[pB[l+j]]->SIZE[k]); }
       }

       matchIndex(QAk, QBk, Ia, Ib);

       for (i=0; i<Ia.len; ++i)
       if (!QAk.recEqual(Ia[i], QBk.rec(Ib[i]))) return 0;
    }

    return 1;
}

template<class TQ, class TD>
template<class TA, class TB>
QSpace<TQ,TD>& QSpace<TQ,TD>::initIdentityCG(
   wbvector< const QSpace<TQ,TA>* > A, const wbindex &Ia,
   wbvector< const QSpace<TQ,TB>* > B, const wbindex &Ib,
   char vflag
){
   if (!cPVEC_nnz(FLF,A,'w') || !cPVEC_nnz(FLF,B,'w')) {
      wblog(FL,"ERR %s() got empty input",FCT); 
   }

   unsigned i, QDIM, dim2=A[0]->QIDX.dim2;
   wbMatrix<TQ> Qa,Qb;
   wbvector<widx_t> Sa,Sb;
   QVec qt0=A[0]->qtype; 
   QMap<TQ> M;

   itag_ itA, itB;
     itA.init(FL,A,Ia); 
     itB.init(FL,B,Ib);

   init(); QDIM=A[0]->QDIM;

   for (i=1; i<A.len; ++i) if (A[i]->qtype!=qt0) wblog(FL,
      "ERR qtype mismatch of input space 1 (%d: %s; %s)",
       i+1, A[i]->qStr().data, qStr().data);
   for (i=0; i<B.len; ++i) if (B[i]->qtype!=qt0) wblog(FL,
      "ERR qtype mismatch of input space 2 (%d: %s; %s)",
       i+1, B[i]->qStr().data, qStr().data);
   for (i=1; i<A.len; ++i)
   if (A[i]->QDIM!=QDIM || A[i]->QIDX.dim2!=dim2) wblog(FL,
      "ERR %s() size mismatch (%d: %d/%d %d/%d)",
       i+1, A[i]->QDIM, QDIM, A[i]->QIDX.dim2, dim2);
   for (i=1; i<B.len; ++i)
   if (B[i]->QDIM!=QDIM || B[i]->QIDX.dim2!=B[0]->QIDX.dim2) wblog(FL,
      "ERR %s() size mismatch (%d: %d/%d %d/%d)",
       FCT,i+1, B[i]->QDIM, QDIM, B[i]->QIDX.dim2, B[0]->QIDX.dim2);

   getQDimGen(A,Qa,Sa,Ia); 
   getQDimGen(B,Qb,Sb,Ib);

   gCS.getQfinal(FL,qt0,Qa,Qb,M,'c'); 

   M.getIdentityQ(FL,Sa,Sb, *this, vflag);

   if (itags.len!=3) wblog(FL,"ERR %s() itags.len=%d",FCT,itags.len);
   if (!itags[0].t) itags[0]=itA;
   if (!itags[1].t) itags[1]=itB;

   ctime=Wb::getTimeNow(); 

   return *this;
};

template<class TQ, class TD>
template<class TA>
QSpace<TQ,TD>& QSpace<TQ,TD>::initIdentityCG(
   wbvector< const QSpace<TQ,TA>* > A, const wbindex &ia, 
   char zflag 
){
   if (!cPVEC_nnz(FLF,A,'l')) { init(); return *this; }

   unsigned i,dim2=A[0]->QIDX.dim2;
   wbMatrix<TQ> Qa,Qb;
   wbvector<widx_t> Sa,Sb;
   wbMatrix<widx_t> Sc;

   init(); QDIM=A[0]->QDIM; qtype=A[0]->qtype;

   for (i=1; i<A.len; ++i) {
      if (A[i]->qtype!=qtype) wblog(FL,
         "ERR qtype mismatch of input space 1 (%d: %s; %s)",
          i+1, A[i]->qStr().data, qStr().data);
      if (A[i]->QDIM!=QDIM || A[i]->QIDX.dim2!=dim2) { wblog(FL,
         "ERR %s() size mismatch (%d: %d/%d %d/%d)",
          i+1, A[i]->QDIM, QDIM, A[i]->QIDX.dim2, dim2);
      }
   }

   getQDimGen(A,Qa,Sa,ia,&Sc);

   if (qtype.isEmpty()) {
      qtype.init(QDIM);
      for (i=0; i<QDIM; ++i) qtype[i]="A";
   }

   unsigned j=0; TQ *qij;
   char cgflag=qtype.isNonAbelian(); 
   wbvec<unsigned> d(qtype.len);

   itags.init(2);
   itags[1]=itags[0].init(FL,A,ia);
   itags[1].Conj();

   QIDX.Cat(2,Qa,Qa); 
   setupDATA(); if (cgflag) {
   setupCGR(); } 

   for (; j<qtype.len; ++j) d[j]=qtype[j].qlen();

   for (i=0; i<DATA.len; ++i) {
      DATA[i]->initIdentity(Sa[i]);

      if (cgflag)
      for (qij=Qa.rec(i), j=0; j<CGR.dim2; ++j) {
         CGR(i,j).initIdentityR(
            FL, qtype[j], qij, Sc.data ? Sc(i,j) : -1, cgflag);
         qij+=d[j];
      }
   }

   ctime=Wb::getTimeNow(); 
   if (!zflag) { return *this; }

   if (zflag>3) { 
      if (zflag=='z') { zflag=1; } else
      if (zflag=='Z') { zflag=2; } else
      wblog(FL,"WRN %s() unexpected zflag=%s",FCT,cSTR(zflag));
   }
   else if (zflag<-3) { 
      wblog(FL,"WRN %s() unexpected zflag=-%s",FCT,cSTR(-zflag));
   }

   itags[1].Conj(); 
   if (abs(zflag)==1) { itags[1].MarkDual(); } else 
   if (abs(zflag)==3) { itags[0].MarkDual(); } 

   for (i=0; i<DATA.len; ++i) {
      for (qij=QIDX.rec(i), j=0; j<qtype.len; ++j) {
         qtype[j].getDual(qij,qij+QDIM); if (cgflag) {
         CGR(i,j).initIdentity1J(
            FL, qtype[j], qij, Sc.data ? Sc(i,j) : -1, cgflag); }
         qij+=d[j];
      }
   }

   if (zflag<0) { Conj(); } 

   if (qtype.allU1()) {
      qtype.init();
   }

   return *this;
};

template <class TQ, class TD>
TD QSpace<TQ,TD>::norm2(char checks) const {

   if (gotCGS(FL)>0) {
      unsigned i=0, j=0; TD c2, x2=0;
      for (i=0; i<DATA.len; ++i) {
         for (c2=1,j=0; j<CGR.dim2; ++j) {
            if (CGR(i,j)) { c2*=TD(CGR(i,j).norm2(checks)); }
         }
         x2+=c2*(DATA[i]->norm2());
      }
      return x2;
   }
   else {
      TD x2=0;
      for (unsigned i=0; i<DATA.len; ++i) x2+=(DATA[i]->norm2());
      return x2;
   }
};

template <class TQ, class TD>
void QSpace<TQ,TD>::checkNorm(
   const char *F, int L, double nrm, char vflag,
   double eps 
 ) const {

   double dbl=norm();

   if (fabs(dbl-nrm)>eps) wblog(F_L,
      "ERR invalid norm %.4g/%.4g (%.3g)",dbl,nrm,(nrm-dbl)/nrm);
   else if (vflag) wblog(F,L,
      "TST norm = %.4g/%.4g (%.3g)",dbl,nrm,(nrm-dbl)/nrm
   );
};

template <class TQ, class TD>
unsigned QSpace<TQ,TD>::RemoveZLabels(
  const char *F, int L, wbMatrix<TQ> &QQ, wbMatrix<TQ> *Z
) const {

   if (QDIM<=0) wblog(FL,"ERR %s() missing QDIM (=%d)",FCT,QDIM);

   unsigned mz=qtype.Qrank(), m0=qtype.Qlen(), m=m0+mz;

   if (&QQ==&QIDX) { if (QDIM!=m) wblog(F_L,
      "ERR no zlabels included in *this (%d <> %d+%d)",QDIM,m0,mz); }
   else if (QDIM!=m0) wblog(F_L,
      "ERR QSpace inconsistency (%d/%d (+%d))",QDIM,m0,mz);
   if (QQ.dim2%m) wblog(F_L,
      "ERR incompatible QQ-label space (%dx%d/(%d+%d))",
       QQ.dim1,QQ.dim2,m0,mz);

   if (QQ.isEmpty()) {
      QQ.dim1=0; if (m && QQ.dim2)
      QQ.dim2=m0*(QQ.dim2/m);
      return m0;
   }

   int r=QQ.dim2/m;
   unsigned i,j,n, l=0, i0=0, iz=0, N=QQ.dim1*r;

   wbvector<unsigned> d0_,dz_;
   wbvector<char> isz_(m);
   char *isz=isz_.data;

   TQ *qz=nullptr, *qq=QQ.data;

   qtype.Qlen(d0_,dz_);
   if (m0!=d0_.sum() || mz!=dz_.sum()) wblog(FL,
      "ERR %s() qdim-inconsistency: %d/%d, %d/%d",
       FCT, m0, d0_.sum(), mz, dz_.sum());

   if (Z) { Z->init(QQ.dim1,mz*r); qz=(Z ? Z->data : nullptr); }
   if (!mz) return m0;

   for (l=i=0; i<qtype.len; ++i)
   for (l+=d0_[i], n=dz_[i], j=0; j<n; ++j) isz[l++]=1;

   for (l=i=0; i<N; ++i) {
      for (j=0; j<m; ++j, ++l) { 
         if (isz[j]==0) qq[i0++]=qq[l];
         else { if (qz) qz[iz]=qq[l]; ++iz; }
      }
   }

   if (!i0) wblog(F_L,"ERR removing z-labels leaves nothing !?");
   else if (i0 != QQ.dim1*m0*r || iz != QQ.dim1*mz*r) wblog(F_L,
     "ERR %s() size inconsistency\n%dx(%d*%d)=%d =>\n"
     "[ %dx(%d*%d)) = %d/%d ] + [ %dx(%d*%d)) = %d/%d ]",FCT,
       QQ.dim1,m, r, QQ.dim1*QQ.dim2,
       QQ.dim1,m0,r, QQ.dim1*m0*r, i0,
       QQ.dim1,mz,r, QQ.dim1*mz*r, iz
   );

   QQ.dim2=m0*r;
   return m0;
};

template <class TQ, class TD> 
void QSpace<TQ,TD>::initOpZ_WET(const char *F, int L,
   const QVec &qvec,
   const wbMatrix<TQ> &Q1, const wbMatrix<TQ> &Q2, const wbMatrix<TQ> &Q,
   const wbarray<TD> &D3,
   const double eps1, const double eps2
){
   unsigned d,i,j,k,l,id,nd, k1,k2, is,i0, r=0, im; 
   unsigned i1,i2,i3, j1,j2,j3, l1,l2,l3, N1,N2,N3, n1,n2,n3;
   unsigned qdimz, qdimz3, dz;

   wbvector<widx_t> d1,d2,d3, I1,I2,I3, M1,M2,M3, s1,s2,s3, D,I,S;
   wbvector<unsigned> m;
   wbvector<TD> d0,dc;
   wbvector<TD> cc;
   wbarray<char> mark;
   wbIndex m1;
   double x;
   TD w;

   wbMatrix<TQ> Q1r(Q1),Q2r(Q2),Q3r(Q), Q1u(Q1),Q2u(Q2),Q3u(Q);
   wbMatrix<TQ> QZ,QX,Z,q2,q3,qc; qset<TQ> q1;
   wbMatrix<double> x1,x2;
   const wbperm P3("312");
   wbperm p1,p2,p3, P;
   wbindex J1,J2,J3;
   QMap<TQ> M;

   clearQSpace(); qtype=qvec; QDIM=qvec.Qlen(); dz=qvec.Qrank();

   qdimz=QDIM+dz; 
   qdimz3=3*qdimz;

   if (Q1.dim2!=qdimz || Q2.dim2!=qdimz || Q.dim2!=qdimz) wblog(FL,
      "ERR inconsistency in number of symmetry labels\n"
      "(%d/%d/%d/%d)", Q1.dim2, Q2.dim2, Q.dim2, qdimz);
   if (D3.SIZE.len!=3) wblog(FL,
      "ERR expecting rank-3 array (%s)",D3.sizeStr().data);

   N1=D3.SIZE[0];
   N2=D3.SIZE[1];
   N3=D3.SIZE[2];

   if (Q1.dim1!=N1 || Q2.dim1!=N2 || Q.dim1!=N3) wblog(FL,
      "ERR dimension mismatch: symmetry labels <> data\n"
      "(%dx%dx%d <> %s)",Q1.dim1,Q2.dim1,Q.dim1,D3.sizeStr().data
   );

   RemoveZLabels(F_L,Q1r); Q1r.groupRecs(p1,d1,Q1u,-1,1,&J1);
   RemoveZLabels(F_L,Q2r); Q2r.groupRecs(p2,d2,Q2u,-1,1,&J2);
   RemoveZLabels(F_L,Q3r); Q3r.groupRecs(p3,d3,Q3u,-1,1,&J3);

   gCS.getQfinal_zdim(F,L,qvec,Q2u,Q3u,Q1u,s2,s3,s1,&m); 
      n1=Q1u.dim1;
      n2=Q2u.dim1;
      n3=Q3u.dim1;

   for (i=0; i<n1; ++i) {
      if (d1[i]%s1[i]==0) d1[i]/=s1[i];
      else wblog(FL,"ERR incomplete WET multiplet space "
        "(i1=%d/%d: %d/%d; %d) !?",i+1,n2,d1[i],s1[i],Q1.dim1
      );
   }
   for (i=0; i<n2; ++i) {
      if (d2[i]%s2[i]==0) d2[i]/=s2[i];
      else wblog(FL,"ERR incomplete WET multiplet space "
        "(i2=%d/%d: %d/%d; %d) !?",i+1,n2,d2[i],s2[i],Q2.dim1
      );
   }
   for (i=0; i<n3; ++i) {
      if (d3[i]%s3[i]==0) d3[i]/=s3[i];
      else wblog(FL,"ERR incomplete WET multiplet space "
        "(i3=%d/%d: %d/%d; %d) !?",i+1,n3,d3[i],s3[i],Q.dim1
      );
   }

   M1.init(N1); M2.init(N2); M3.init(N3);
   I1.init(n1); I2.init(n2); I3.init(n3); 

   for (i=0; i<N1; i+=d) { j=J1[i]; d=s1[j];
      for (l=1; l<d; ++l) {
         if (Q1r.recCompare(i,i+l)) wblog(FL,
          "ERR %s() multiplet not grouped in z-labels (1: %d/%d)",
           FCT,i,Q1.dim1);
      }
      for (id=I1[j]++, l=0; l<d; ++l) M1[i+l]=id;
   }
   if (I1!=d1) wblog(FL,"ERR %s() missed multiplets !?",FCT);

   for (i=0; i<N2; i+=d) { j=J2[i]; d=s2[j];
      for (l=1; l<d; ++l) {
         if (Q2r.recCompare(i,i+l)) wblog(FL,
          "ERR WET: multiplet not grouped in z-labels (2: %d/%d)",
           i,Q2.dim1);
      }
      for (id=I2[j]++, l=0; l<d; ++l) M2[i+l]=id;
   }
   if (I2!=d2) wblog(FL,"ERR %s() missed multiplets !?",FCT);

   for (i=0; i<N3; i+=d) { j=J3[i]; d=s3[j];
      for (l=1; l<d; ++l) {
         if (Q3r.recCompare(i,i+l)) wblog(FL,
          "ERR WET: multiplet not grouped in z-labels (3: %d/%d)",
           i,Q.dim1);
      }
      for (id=I3[j]++, l=0; l<d; ++l) M3[i+l]=id;
   }
   if (I3!=d3) wblog(FL,"ERR %s() missed multiplets",FCT);

   n3=D3.numel(); mark.init(D3.SIZE);
   for (nd=i=0; i<n3; ++i) { x=fabs(D3.data[i]);
      if (x>eps1) { ++nd; mark.data[i]=1; }
      else if (x>eps2) wblog(F,L,
      " *  ignoring small value data[%d]=%.3g @ %g",i+1,D3.data[i],eps2);
   }
   if (!nd) {
      wblog(FL,"WRN %s() got all-zero matrix elements !?",FCT);
      init(); return; 
   }

   QZ.init(nd,qdimz3); d0.init(nd);
   QIDX.init(nd,3*QDIM);
   setupDATA(FL); setupCGR(FL);

   for (i3=0; i3<N3; i3+=n3) { j3=J3[i3]; n3=s3[j3];
   for (i2=0; i2<N2; i2+=n2) { j2=J2[i2]; n2=s2[j2];
   for (i1=0; i1<N1; i1+=n1) { j1=J1[i1]; n1=s1[j1];
      d0.len=QZ.dim1=nd; l=0;

      for (k=0; k<n3; ++k) { l3=i3+k;
      for (j=0; j<n2; ++j) { l2=i2+j;
      for (i=0; i<n1; ++i) { l1=i1+i; if (!mark(l1,l2,l3)) continue;
         QZ.recSetB(l,0,qdimz,Q1.rec(l1));
         QZ.recSetB(l,1,qdimz,Q2.rec(l2));
         QZ.recSetB(l,2,qdimz,Q .rec(l3)); d0[l++]=D3(l1,l2,l3);
      }}}
      if (!l) continue;

      d0.len=QZ.dim1=l; 
      QX=QZ; QX.SkipTiny_float();
      RemoveZLabels(F_L,QX,&Z);

      q1.init(  QDIM, QX.ref(0,     0));
      q2.init(1,QDIM, QX.ref(0,  QDIM));
      q3.init(1,QDIM, QX.ref(0,2*QDIM));

      gCS.getQfinal(F,L,qvec,q2,q3,M,QF_LOADC);

      m1.init();
      M.getCGZlist(F,L,qc,dc,&I,nullptr,&q1,&m1); 

      if (!I.allEqual(I[0])) wblog(FL,"ERR %s() CData not unique",FCT);

      qc.BlockPermute(P3);  cc.init(1); cc[0]=1;
      QZ.isref=d0.isref=1; {
         QZ.SkipTiny_float(); QZ.SortRecs(P); d0.Permute(P);
         qc.SkipTiny_float(); qc.SortRecs(P); dc.Permute(P);
      }; QZ.isref=d0.isref=0;

      if (m1.numel()>1) {
         unsigned l, mq=m1.numel();
         wbvector< wbMatrix<TQ> > QC(mq+1);
         wbvector< wbvector<TD> > DC(mq+1);
         wbindex Ia,Ib;

         wblog(FL," *  %s() got outer multiplicity [%s]",FCT,STR(m1.SIZE));

         for (l=0; l<=mq; ++l) {
            if (l) {
               if (!M.getCGZlist(F,L,qc,dc,&I,nullptr,&q1,&m1)) break;
               qc.BlockPermute(P3);
               qc.SkipTiny_float(); qc.SortRecs(P); dc.Permute(P);
            }
            qc.findUnique(Ia);
            qc.getRecs(Ia,QC[l]);
            dc.select(Ia,DC[l]);
         }
         if (l!=mq) wblog(FL,"ERR %s() !? (%d/%d)",FCT,l,mq);

         QZ.findUnique(Ia);
         QZ.getRecs(Ia,QC[l]);
         d0.select(Ia,DC[l]); ++l;

         qc=QC[0];
         for (l=1; l<QC.len; ++l) {
            matchSortedIdx(qc,QC[l],Ia,Ib);
            QC[l]=QC[l].getRecs(Ib,qc); DC[l].Select(Ib);
            if (qc.isEmpty()) break;
         }
         if (qc.dim1<QC.len) wblog(FL,
            "ERR %s() insufficient number of unique recs (%d/%d)",
            FCT,qc.dim1,QC.len
         ); l=DC.len-2;
         if (DC[l].len!=DC[l+1].len) wblog(FL,
            "ERR %s() mismatch of unique recs (%d/%d)",
            FCT,DC[l].len,DC[l+1].len
         );
         for (l=0; l<QC.len; ++l) {
            matchSortedIdx(qc,QC[l],Ia,Ib);
            if (Ia.len!=qc.dim1) wblog(FL, 
               "ERR %s() %d/%d",FCT,Ia.len,qc.dim1);
            QC[l].Set2Recs(Ib); DC[l].Select(Ib);
         }

         wbarray<TD> b,c,X2,Xc,iX, X(qc.dim1, mq);
         for (l=0; l<mq; ++l) { X.setCol(l,DC[l].data); }

         Wb::MatProd(X,X,X2,'C');
         Wb::MatVProd(X,DC[l],b,'C');
         Wb::VMatProd(wbInverse(FL,X2,iX),b,c,'C');

         cc.init(c.numel(), c.data); m1.reset();
         Wb::VMatProd(X,c,Xc); x=Xc.normDiff(DC[l]);
         if (x>1e-12) wblog(FL,
            "ERR %s() failed to match OM (%.3g)",FCT,x);
         else wblog(FL,
            " *  %s() succeeded to match OM (@ %.2g)\n==> "
            "[%s] (%.5g)",FCT,x,STR(cc),cc.norm2());

         if (!M.getCGZlist(F,L,qc,dc,&I,nullptr,&q1,&m1,&cc)) wblog(FL,
            "ERR %s() failed to get CGZ list (%s)",FCT,STR(m1));

         qc.BlockPermute(P3);
         qc.SkipTiny_float(); qc.SortRecs(P); dc.Permute(P);

      }

      if (d0.len!=dc.len) {

         MXPut X(FL,"a"); wbvector<double> i(3);
         X.add(D3,"D3").add(Q1,"Q1").add(Q2,"Q2").add(Q,"Q")
          .add(Q1r,"Q1r").add(Q2r,"Q2r").add(Q3r,"Q3r")
          .add(Q1u,"Q1u").add(Q2u,"Q2u").add(Q3u,"Q3u")
          .add(QZ,"QZ").add(d0,"d0").add(qc,"qc").add(dc,"dc")
          .add(eps1,"eps1").add(eps2,"eps2")
          .add(q1,"q1").add(q2,"q2").add(q3,"q3").add(m,"m").add(M,"M");
         i[0]=i1; i[1]=i2; i[2]=i3; X.add(i+1,"i");
         i[0]=j1; i[1]=j2; i[2]=j3; X.add(i+1,"j");
         i[0]=n1; i[1]=n2; i[2]=n3; X.add(i,"n");
         i[0]=N1; i[1]=N2; i[2]=N3; X.add(i,"N");
         wblog(FL,"ERR WET data mismatch (len=%d/%d)",d0.len,dc.len);
      }

      x=Wb::sqrt(double(QZ.normDiff2(qc))); if (x>eps1) {
        MXPut(FL,"q").add(dc,"dc").add(qc,"qc").add(QZ,"QZ").add(d0,"d0");
        wblog(FL,"ERR %s() data mismatch (%.3g)",FCT,x);
      }

      k1=QZ.findUniqueRecSorted1(); if (int(k1)<0) {
         MXPut(FL,"a").add(QZ,"QZ").add(k1,"k1");
         wblog(FL,"ERR %s() failed to find unique Q3-record",FCT);
      }
      k2=qc.findRecSorted(QZ.rec(k1)); if (int(k2)<0)
         wblog(FL,"ERR %s() failed to find matching Q3-record",FCT);
      w=d0[k1]/dc[k2];

      x1.init(FL,QZ,d0).SortRecs();
      x2.init(FL,qc,dc); 
      x2.colTimes(x2.dim2-1,w).SortRecs();

      x=x1.normDiff(x2); if (x>eps1) {
         wbvector<double> i(3);
         MXPut X(FL,"a"); X.add(D3,"D3")
          .add(Q1,"Q1").add(Q2,"Q2").add(Q,"Q").add(QZ,"QZ").add(M,"M")
          .add(d0,"d0").add(qc,"qc").add(dc,"dc").add(w,"w");
         i[0]=i1; i[1]=i2; i[2]=i3; X.add(i+1,"i");
         i[0]=n1; i[1]=n2; i[2]=n3; X.add(i,"n");
         i[0]=N1; i[1]=N2; i[2]=N3; X.add(i,"N");
         i[0]=j1; i[1]=j2; i[2]=j3; X.add(i+1,"j");
         X.add(x1,"x1").add(x2,"x2");
        wblog(FL,"ERR WET data mismatch (%.3g)",x);
      }
      x=fabs(w); if (x>1E3 || x<1e-3)
        wblog(FL,"WRN %s() got WET=%.3g !?",FCT,w);

      i0=I[0]; I.init(3); m1.reset();
      for (im=0; im<cc.len; ++im, ++r) {
         QIDX.recSetB(r,0,QDIM,q1.data);  I[0]=M1[i1];
         QIDX.recSetB(r,1,QDIM,q2.data);  I[1]=M2[i2];
         QIDX.recSetB(r,2,QDIM,q3.data);  I[2]=M3[i3];

         if (!(++m1)) wblog(FL,
            "ERR %s() counter on outer multiplicity !?",FCT);

         for (is=0; is<qvec.len; ++is) {
            if (( d = M.cg3(i0,is)->getOM() )!=cc.len) wblog(FL, 
              "ERR %s() QMap mismatch (%d/%d)",FCT,d, cc.len);

            CGR(r,is).init_1( *M.cg3(i0,is), im); 
            CGR(r,is).Permute(P3); 
         }

         wbarray<TD> &a = DATA[r]->init(d1[j1], d2[j2], d3[j3]);
         a.element(I)=w*cc[im]; 

         if (fabs(w*cc[im])<1e-12) wblog(FL,"WRN %s() "
            "got w = %g * %g = %g (im=%d)",FCT,w,cc[im],w*cc[im],im);

      }
   }}}

   if (!r || r>QIDX.dim1) {
      MXPut(FL,"a").add(s1,"s1").add(s2,"s2").add(s3,"s3")
      .add(Q1,"Q1").add(Q2,"Q2").add(Q,"Q")
      .add(Q1u,"Q1u").add(Q2u,"Q2u").add(Q3u,"Qu").add(mark,"mark");
      wblog(FL,"ERR %s (r=%d)",FCT,r); 
   }

   QIDX.dim1 = DATA.len = CGR.dim1 = r;

   otype=QS_OPERATOR;

   itags.init_qdir("-++");

   Conj();

   MakeUnique(); 

};

template <class TQ, class TD>
void QSpace<TQ,TD>::reduceMatEl(
   const char *F, int L,
   const wbMatrix<TQ> &QQ, const wbvector<TD> &dd, 
   const wbMatrix<TQ> &Qc, const wbvector<TD> &dc, 
   wbvector<unsigned> &J,  
   wbvector<unsigned> &Jc, 
   wbvector<TD> &dr        
) const {

   unsigned i,j,k,m,n,r; char e=0;
   wbvector<size_t> d1,d2,D1,D2;
   wbperm pp,pc;
   wbindex I1,I2,i1,i2;
   wbvector<TD> xx,x1,x2;
   double x=0;

   wbMatrix<TQ> Z1,Z2,q2;
   QSpace<TQ,TD> A1,A2;

   if (QQ.dim1!=dd.len || Qc.dim1!=dc.len || QQ.dim2!=Qc.dim2) wblog(F,L,
      "ERR dimension mismatch (%d/%dx%d; %d/%dx%d)",
       dd.len, QQ.dim1, QQ.dim2, dc.len, Qc.dim1, Qc.dim2
   );

   A1.QIDX=QQ; A1.qtype=qtype;
   A1.RemoveZLabels(F,L,&Z1); A1.QIDX.groupRecs(pp,d1,Z1); 

   A2.QIDX=Qc; A2.qtype=qtype;
   A2.RemoveZLabels(F,L,&Z2); A2.QIDX.groupRecs(pc,d2,Z2); 

   i=matchIndex(A1.QIDX,A2.QIDX,I1,I2);
   if (i) wblog(F,L,"ERR got non-unique Q-set (%d)",i);
   if (I1.len!=A1.QIDX.dim1) {
       wbindex Ix; I1.invert(A1.QIDX.dim1,Ix);
       for (i=0; i<Ix.len; ++i) A1.QIDX.recPrint(Ix[i],"Q");
       wblog(F,L,
         "ERR got %d additional symmetry sectors (%d/%d matching)",
          A1.QIDX.dim1-I1.len, I1.len, A1.QIDX.dim1
       );
   }

   pc.BlockSelect(I2,d2); Qc.getRecs(pc,q2); 

   matchIndex(QQ,q2,i1,i2); 

   if (i1.len!=QQ.dim1) { ++e;
       wbindex Ix; i1.invert(QQ.dim1,Ix);
       for (i=0; i<Ix.len; ++i) QQ.recPrint(Ix[i],"QZ");
       wblog(F,L,"==> ERR got additional z-labels (%d/%d)",
       QQ.dim1-i1.len,QQ.dim1);
   }

   if (i2.len!=q2.dim1) { ++e;
       wbindex Ix; i1.invert(q2.dim1,Ix);
       for (unsigned i=0; i<Ix.len; ++i) q2.recPrint(Ix[i],"Qall");
       wblog(F,L,
          "==> ERR got missing z-labels (%d/%d) !?",
           q2.dim1-i2.len,q2.dim1
       );
   }

   if (e) ExitMsg("\n"); 

   if (dd.len<pp.len || dd.len%pp.len || pc.max()>=dc.len) wblog(F,L,
      "ERR index out of bounds (%d/%d, %d/%d)",
       dd.len, pp.len, pc.max(), dc.len
   ); 

   m=dd.len/pp.len; 

   d1.cumsum_(D1); d2.cumsum_(D2); dr.init(m*I1.len);
   J .init(I1.len);
   Jc.init(I1.len);

   for (i=0; i<I1.len; ++i) { n=d1[I1[i]];
      if (n!=d2[I2[i]]) wblog(F,L, 
      "ERR size mismatch (%d/%d)", d1[I1[i]], d2[I2[i]]);

      x1.init(n); x2.init(n); xx.init(n);
      k=D2[i]; Jc[i]=pc[k]; for (j=0; j<n; ++j) x2[j]=dc[pc[k+j]];
      k=D1[i]; J[i]=pp[k];

      for (r=0; r<m; ++r) { 
         for (j=0; j<n; ++j){  x1[j]=dd[pp[k+j]+r*pp.len];

             if (fabs(x1[j])<1e-12 || fabs(x2[j])<1e-12) wblog(F,L,
                "ERR got matrix element (%d,%d): %.4g, %.4g",
                 i+1,j+1,x1[j],x2[j] 
             );

             xx[j]=x1[j]/x2[j];
             if (j) { if (fabs(x-xx[j])>1e-12) {
#ifdef MATLAB_MEX_FILE
                 MXPut(FL).add(QQ,"QQ").add(dd,"dd").add(Qc,"qc")
                 .add(dc,"dc").add(pp+1,"pp").add(pc+1,"pc")
                 .add(D1+1,"d1").add(D2+1,"d2");
#endif
                 wblog(F,L,
                   "ERR inconsistent reduced matrix element (%d,%d/%d; %d,%d/%d):\n"
                   "%.4g/%.4g = %.4g (%.4g) @ %.4g",
                    pp[D1[i]]+r*pp.len+1, pp[D1[i]+j]+r*pp.len+1, dd.len,
                    pc[D2[i]]+1, pc[D2[i]+j]+1, dc.len,
                    x1[j], x2[j], xx[j], x, xx[j]-x
                 );
             }}
             else x=xx[j];
         }
         dr[i+r*I1.len]=xx.avg();
      }
   }
};

template <class TQ, class TD>
QSpace<TQ,TD>& QSpace<TQ,TD>::Select(const wbindex& I, char data_only) {

   widx_t *ix; wbindex Ix;
   char gotcgr=(!CGR.isEmpty());

   if (data_only) {
      if (QIDX.dim1!=I.len)
         wblog(FL,"ERR %s() dimension mismatch (%d/%d; %d/%d)",
         FCT, QIDX.dim1, I.len, CGR.dim1, DATA.len
      );
   }

   if (gotcgr && (CGR.dim1!=DATA.len || CGR.dim2!=qtype.len)) wblog(FL,
      "ERR %s() CGR size mismatch (%dx%d ; %dx%d)",
      FCT, CGR.dim1, CGR.dim2, DATA.len, qtype.len
   );

   if (isref) wblog(FL,"ERR %s() got QSpace isref=%d",FCT,isref);
   if (!I.isUnique(DATA.len)) wblog(FL,
      "ERR %s() got invalid index (%d/%d)",FCT,I.max()+1,DATA.len);

   I.invert(DATA.len,Ix); ix=Ix.data;
   if (!data_only) QIDX.Set2Recs(I);

   if (DATA.len) {
      wbvector< wbarray<TD>* > DD(DATA); DD.select(I,DATA);
      for (unsigned i=0; i<Ix.len; ++i) {
         wbarray<TD>* &d=DD[ix[i]]; WB_DELETE_1(d);
      }
   }

   if (gotcgr) { CGR.Set2Recs(I); }
   return *this;
};

template <class TQ, class TD> inline
void QSpace<TQ,TD>::Append(
   const wbvector<TQ> &Q,
   const wbarray<TD> &D
){

   wbarray<TD> *d=0; WB_NEW_1(d); d->init(D);

   QIDX.appendRow(Q); 
   DATA.Append(d);    
};

template <class TQ, class TD> inline
void QSpace<TQ,TD>::Append(unsigned k) {
   unsigned i, n=DATA.len;

   QIDX.Resize(QIDX.dim1+k, QIDX.dim2);
   DATA.Resize(n+k);

   for (i=0; i<k; ++i) {
   WB_NEW_1(DATA[n+i]); }
}

template <class TQ, class TD>
int QSpace<TQ,TD>::ExpandOM(
   const char *F, int L, wbvector<unsigned> &iOM,
   char full 
){
   int rval=0; 

   unsigned nsym=qtype.len, r=rank(F_L);
   iOM.init(qtype.len); if (!iOM || r<=2) { return rval; }

   if (CGR.dim2!=nsym || CGR.dim1!=DATA.len) wblog(F_L,
      "ERR %s() CGR size mismatch %s / %dx%d",FCT,SSTR(CGR),DATA.len,nsym);

   unsigned i,j=0, l=r; 
   for (; j<nsym; ++j) {
      if (qtype[j].permitsOM(r)) {
          for (i=0; i<CGR.dim1; ++i) {
             if (CGR(i,j).wOM()) { iOM[j]=l++; break; }
          }
      }
   }
   if (l==r) { 
      for (i=0; i<DATA.len; ++i) { wbarray<TD> &Ai=(*DATA[i]);
         if ((l=Ai.SIZE.len)>r) {
            if (Ai.skipSingletons(0,0,r)<0) wblog(FL,"ERR %s() "
               "unexpected DATA[%d/%d] size %s",FCT,i,DATA.len,SSTR(Ai));
            if (rval<int(l-=r)) { rval=l; }
         }
      }
      return (rval=-rval);  
   }
   rval=l-r; if (rval==1 && !full) { return rval; }

   if (full=='f') { full=1; } else
   if (full=='F') { full=2; } else
   if (full>2) wblog(FL,"WRN %s() unexpected full=%s",FCT,cSTR(full));

   unsigned m=l-r, M, ri;
   wbvector<unsigned> S(m), Sx; if (full>1) { Sx.init(m); }

   for (i=0; i<DATA.len; ++i) {
      wbarray<TD> &Ai=(*DATA[i]); ri=Ai.SIZE.len;

      for (l=j=0; j<CGR.dim2; ++j) { if (iOM[j]) {
          S [l]=CGR(i,j).wdim2(); if (full>1) {
          Sx[l]=CGR(i,j).wdim1(); }; ++l;
      }}
      M=S.prod(); 

      if (ri==r) {
         if (M!=1) wblog(F_L,
            "ERR %s(%d/%d) OM size mismatch %s @ %d (r=%d+%d)",
            FCT, i+1,DATA.len, SSTR(Ai), M,r,S.len
         );
         Ai.ExpandOM(FL,r,S, Sx.data); 
      }
      else if (ri==r+1) {
         if (M!=Ai.SIZE[r]) wblog(F_L,
            "ERR %s(%d/%d) OM size mismatch %s @ %d (r=%d+%d)",
            FCT, i+1,DATA.len, SSTR(Ai), M,r,S.len);

         if (S.len>1) { 
         Ai.ExpandOM(FL,r,S, Sx.data); } 
      }
      else if (ri==r+S.len) { 
         for (j=0; j<S.len; ++j) { if (Ai.SIZE[r+j]!=S[j]) break; }
         if (j<S.len) wblog(FL,
            "ERR %s(%d/%d) OM size mismatch %s / %s (r=%d+%d)",
            FCT,i+1,DATA.len, SSTR(Ai), SSTR(S),r,S.len);
      }
      else {
         wblog(FL,"ERR %s(%d/%d) unexpected rank %s / %s (r=%d+%d)",
         FCT,i+1,DATA.len, SSTR(Ai), SSTR(S),r,S.len);
      }
   }

   return rval;
};

template <class TQ, class TD>
int QSpace<TQ,TD>::ExpandOM_(const char *F, int L,
   unsigned ia, QSpace<TQ,TD> &B, unsigned ib 
) {
   unsigned r=rank(FL);

   if (ia>=   DATA.len) wblog(F_L,
      "ERR %s() index out of bounds (A: %d/%d)",FCT,ia+1,  DATA.len);
   if (ib>=B.DATA.len) wblog(F_L,
      "ERR %s() index out of bounds (B: %d/%d)",FCT,ib+1,B.DATA.len);
   if (!CGR) { return -1; }

   if (this!=&B) { 
      if (QIDX.dim2!=B.QIDX.dim2 || CGR.dim2!=B.CGR.dim2) wblog(FL,
      "ERR %s() QSpace mismatch (r=%d/%d)",FCT,r,B.rank(FL));
   }

   wbarray<TD> &Ai=(*DATA[ia]), &Bi=(*B.DATA[ib]);
   unsigned j,n1,n2, ma=Ai.numOM(r), mb=Bi.numOM(r), nsym=CGR.dim2;

   if (ma<=1 && mb<=1) {
      if (!Ai.sameSize(Bi)) wblog(FL,"ERR %s() got size difference "
         "(%d,%d): %s / %s", FCT,ia+1,ib+1,SSTR(Ai),SSTR(Bi));

      for (j=0; j<nsym; ++j) { 
         const CRef<TQ> &a=CGR(ia,j), &b=B.CGR(ib,j);
         n1=a.wnumel(); n2=b.wnumel();
         if (n1>1 || n2>1) { break; } 
         if (n1 && n2) { 
            if (fabs(a.cgw[0]-b.cgw[0])>1e-12) { wblog(FL,
               "ERR %s() cgw normalization mismatch %.5g / %.5g @ %.3g",
               FCT,a.cgw[0],b.cgw[0],fabs(a.cgw[0]-b.cgw[0]));
            }
         }
      }
      if (j>=nsym)
           { return 0; } 
   }

   wbvector<unsigned> Sx(nsym), jOM(nsym); 
   wbMatrix<unsigned> Sa(2,nsym), Sb(2,nsym);

   unsigned *sa=Sa.data, *sa_=Sa.data+nsym, *sb=Sb.data, *sb_=Sb.data+nsym;
   unsigned i, m=0, Ma=-1, Mb=-1;

   for (j=0; j<nsym; ++j) {
      const CRef<TQ> &a=CGR(ia,j), &b=B.CGR(ib,j);
      n1=a.wnumel(); n2=b.wnumel();
      if (n1>1 || n2>1) {
         if (!a.cgb || a.cgb!=b.cgb || !n1 || !n2) { wblog(FL,
            "ERR %s() cgb mismatch or null %p / %p\n   %s\n<> %s",
            FCT,a.cgb,b.cgb, STR(a),STR(b));
         }
         if (!a.wSame(b.cgw,'l')) { 
            wblog(FL,"ERR %s() cgw normalization mismatch %.5g / %.5g @ %.3g",
            FCT,a.cgw[0],b.cgw[0],fabs(a.cgw[0]-b.cgw[0]));
         }

         sa_[m] = ma = a.wdim12(FL,sa[m]);   
         sb_[m] = mb = b.wdim12(FL,sb[m]);   

         if (ma>1 && Ma>ma) { Ma=ma; } 
         if (mb>1 && Mb>mb) { Mb=mb; }

         Sx[m]=MAX(ma,mb);
         jOM[m]=j; ++m;
      }
   }

   if (!m) wblog(FL,"ERR %s() got m=%d",FCT,m);
   Sx.Shorten2(m);

   if (int(Ma)<0) { Ma=0; }
   if (int(Mb)<0) { Mb=0; }

   if (Sa==Sb) {
      #ifndef WB_SKIP_ASSERT
       if (!Ai.sameSize(Bi)) wblog(FL,
          "ERR %s() got size difference (%d,%d): %s / %s",
          FCT,ia+1,ib+1, SSTR(Ai),SSTR(Bi));
      #endif
      return 0;
   }

   Ai.ExpandOM(FL,r,m,sa); 
   Bi.ExpandOM(FL,r,m,sb); 

   if (Ma || Mb) {
      for (i=0; i<m; ++i) { if (Sx[i]>1) { 
         wbarray<double> &a=CGR(ia,jOM[i]).cgw;
         if (!a.isIdentityMatrix()) {
            Ai.ContractMat(FL,r+i+1,a,2); 
         }; a.initIdentity(Sx[i]); 
      }}
   }
   if (!Ai.sameSizeM(r,m,sa_)) { wblog(FL, 
      "ERR %s() size mismatch (%s: %s -> %s @ r=%d)",
      FCT,SSTR(Ai), STR(wbvector<unsigned>(m,sa_,'r')), STR(Sx),r);
   }
   Ai.ResizeM(r,Sx).FuseOM(FL,r);

   if (Mb || Ma) {
      for (i=0; i<m; ++i) { if (Sx[i]>1) { 
         wbarray<double> &b=B.CGR(ib,jOM[i]).cgw;
         if (!b.isIdentityMatrix()) {
            Bi.ContractMat(FL,r+i+1,b,2); 
         }; b.initIdentity(Sx[i]); 
      }}
   }
   if (!Bi.sameSizeM(r,m,sb_)) { wblog(FL, 
      "ERR %s() size mismatch (%s: %s -> %s @ r=%d)",
      FCT,SSTR(Bi), STR(wbvector<unsigned>(m,sb_,'r')), STR(Sx),r);
   }
   Bi.ResizeM(r,Sx).FuseOM(FL,r);

   return m;
};

template <class TQ, class TD>
QSpace<TQ,TD>& QSpace<TQ,TD>::Append(
   const char *F, int L, const QSpace<TQ,TD> &B, char uflag
){
   if (B.isEmpty()) return *this;
   if (this==&B) wblog(F_L,"ERR %s() requires distinct B",FCT);
   if (isEmpty()) { *this = B; return *this; }

   unsigned r=-1;

   if (!  isConsistent(r)) wbdie(F,L,str);
   if (!B.isConsistent(r)) wbdie(F,L,str);

   if (QDIM!=B.QDIM || QIDX.dim2!=B.QIDX.dim2) wblog(F,L,
      "ERR %s() rank mismatch (%d,%d; %d,%d)",
       FCT,QDIM, B.QDIM, QIDX.dim2, B.QIDX.dim2
   );
   if (B.QIDX.dim1!=B.DATA.len || QIDX.dim1!=DATA.len) wblog(FL,
      "ERR %s() QIDX/DATA size inconsistency\n%d/%d; %d/%d",
       FCT,B.QIDX.dim1,B.DATA.len,QIDX.dim1,DATA.len
   );
   if (qtype!=B.qtype) wblog(F,L,
      "ERR %s() qtype inconsistency (%s; %s)",
       FCT,B.qStr().data,qStr().data
   );
   if (itags!=B.itags) wblog(F_L, 
      "ERR %s() itag inconsistency ('%s' <> '%s')",
      FCT, STR(itags), STR(B.itags)
   );
   if (!CGR.isEmpty() || !B.CGR.isEmpty()) {
      if (qtype.isEmpty()) wblog(FL,
         "ERR %s() got empty qtype %s (CGR: %dx%d; %dx%d)",FCT,
          qStr().data,B.CGR.dim1,B.CGR.dim2,CGR.dim1,CGR.dim2
      );
      if (B.CGR.dim1!=B.QIDX.dim1 || CGR.dim1!=QIDX.dim1 ||
          B.CGR.dim2!=B.qtype.len || CGR.dim2!=B.CGR.dim2) wblog(FL,
          "ERR %s() QSpace inconsistency\nCGS: %dx%d; %dx%d/%d",
          FCT, B.CGR.dim1, B.CGR.dim2, CGR.dim1, CGR.dim2, qtype.len
      );
   }

   QIDX.appendRows(B.QIDX.dim1, B.QIDX.data);

   unsigned i=DATA.len, n;
   wbarray<TD> *a=nullptr;

   DATA.Append(B.DATA.len, B.DATA.data); 
   for (n=DATA.len; i<n; ++i) { 
       WB_NEW_2(a,wbarray<TD>(*DATA[i]));
       DATA[i]=a;
   }

   if (!CGR.isEmpty())
   CGR.appendRows(B.CGR.dim1,B.CGR.data);

   if (uflag) { MakeUnique(); }

   return *this;
};

template <class TQ, class TD>
QSpace<TQ,TD>& QSpace<TQ,TD>::Cat( 
   const char *F, int L,
   const QSpace<TQ,TD> &A, const QSpace<TQ,TD> &B,
   TD afac, TD bfac, char uflag
){
   if (A.isEmpty()) {
      if (this!=&B   ) { (*this)=B;     }
      if (bfac!=TD(1)) { (*this)*=bfac; }; return *this;
   }
   if (B.isEmpty()) {
      if (this!=&A   ) { (*this)=A;     }
      if (afac!=TD(1)) { (*this)*=afac; }; return *this;
   }
   if (&A==this || &B==this) {
      QSpace<TQ,TD> X; X.Cat(F,L,A,B,afac,bfac,uflag);
      return X.save2(*this);
   }

   unsigned r=-1;
   if (!A.isConsistent(r)) { wbdie(F,L,str); }
   if (!B.isConsistent(r)) { wbdie(F,L,str); }

   if (A.QDIM!=B.QDIM || A.QIDX.dim2!=B.QIDX.dim2) wblog(F,L,
      "ERR %s() rank mismatch (%d,%d; %d,%d)",
      FCT, A.QDIM, B.QDIM, A.QIDX.dim2, B.QIDX.dim2);

   if (A.QIDX.dim1!=A.DATA.len || B.QIDX.dim1!=B.DATA.len) wblog(FL,
      "ERR %s() QIDX/DATA size inconsistency\n%d/%d; %d/%d",
      FCT,A.QIDX.dim1,A.DATA.len,B.QIDX.dim1,B.DATA.len);

   if (A.qtype!=B.qtype) wblog(F,L,
      "ERR %s() qtype inconsistency (%s; %s)",
      FCT,A.qStr().data,B.qStr().data);

   if (A.itags!=B.itags) wblog(F_L,
      "ERR %s() itag inconsistency ('%s' <> '%s')",
      FCT, STR(A.itags), STR(B.itags));

   if (!A.CGR.isEmpty() || !B.CGR.isEmpty()) {
      if (A.qtype.isEmpty()) wblog(FL,
         "ERR %s() got empty qtype %s (CGR: %dx%d; %dx%d)",FCT,
         A.qStr().data,A.CGR.dim1,A.CGR.dim2,B.CGR.dim1,B.CGR.dim2
      );
      if (B.CGR.dim1!=B.QIDX.dim1 || A.CGR.dim1!=A.QIDX.dim1 ||
         B.CGR.dim2!=B.qtype.len || A.CGR.dim2!=B.CGR.dim2) wblog(FL,
         "ERR %s() QSpace inconsistency\nCGS: %dx%d; %dx%d/%d",
         FCT,A.CGR.dim1, A.CGR.dim2, A.CGR.dim1, A.CGR.dim2, A.qtype.len
      );
   }

   initQ(A);

   QIDX.Cat(1,A.QIDX,B.QIDX);
   CGR .Cat(1,A.CGR, B.CGR );

   unsigned i,n, na=A.DATA.len;
   wbarray<TD> *a;
   TD one=1;

   DATA.Cat(A.DATA, B.DATA); 
   for (n=DATA.len, i=0; i<n; ++i) {  
       WB_NEW_2(a,wbarray<TD>(*DATA[i]));
       DATA[i]=a;
   }

   if (afac!=one) { for (i=0; i<na; ++i) (*DATA[i])*=afac; }
   if (bfac!=one) { for (i=na; i<n; ++i) (*DATA[i])*=bfac; }

   if (uflag) { MakeUnique(); }

   return *this;
};

template <class TQ, class TD>
int QSpace<TQ,TD>::Append2AndDestroy(
   const char *F, int L, QSpace<TQ,TD> &A, char unique
){
   if (isEmpty()) return 0;
   if (this==&A) wblog(F,L,"ERR %s() requires distinct A",FCT);

   unsigned r=-1; int e=0;

   if (A.isEmpty()) { save2(A); return e; }

   if (!  isConsistent(r,3)) { wbdie(F,L,str); }
   if (!A.isConsistent(r,2)) { wbdie(F,L,str); }

   if (QDIM!=A.QDIM || QIDX.dim2!=A.QIDX.dim2) wblog(F,L,
      "ERR %s() rank mismatch (%d,%d; %d,%d)",
       FCT,QDIM, A.QDIM, QIDX.dim2, A.QIDX.dim2
   );
   if (A.QIDX.dim1!=A.DATA.len || QIDX.dim1!=DATA.len) wblog(FL,
      "ERR %s() QIDX/DATA size inconsistency\n%d/%d; %d/%d",
       FCT,A.QIDX.dim1,A.DATA.len,QIDX.dim1,DATA.len
   );
   if (A.qtype!=qtype) wblog(F,L,
      "ERR %s() qtype inconsistency (%s; %s)",
       FCT,A.qStr().data,qStr().data
   );
   if (!CGR.isEmpty() || !A.CGR.isEmpty()) {
      if (qtype.isEmpty()) wblog(FL,
         "ERR %s() got empty qtype %s (CGR: %dx%d; %dx%d)",FCT,
          qStr().data,A.CGR.dim1,A.CGR.dim2,CGR.dim1,CGR.dim2
      );
      if (A.CGR.dim1!=A.QIDX.dim1 || CGR.dim1!=QIDX.dim1 ||
          A.CGR.dim2!=A.qtype.len || CGR.dim2!=A.CGR.dim2) wblog(FL,
          "ERR %s() QSpace inconsistency\nCGS: %dx%d; %dx%d/%d",
          FCT, A.CGR.dim1, A.CGR.dim2, CGR.dim1, CGR.dim2, qtype.len
      );
   }

   if (unique && A.hasQOverlap(*this)) { ++e;
      sprintf_str("%s ERR objects have QIDX overlap", shortFL(F,L));
   }

   A.QIDX.appendRows(QIDX.dim1, QIDX.data);

   A.DATA.Append(DATA.len, DATA.data); if (!CGR.isEmpty()) {
   A.CGR.appendRows(CGR.dim1,CGR.data); }

   QIDX.init(0,QIDX.dim2); DATA.init(); CGR.init(); init();

   return e;
};

template <class TQ, class TD> inline
unsigned QSpace<TQ,TD>::SkipZeroData( 
   double eps, char bflag, char cgflag,
   char all 
){

   if (!QIDX.dim2 && QIDX.dim1==1 && DATA.len==1 && CGR.dim1==1) {
      unsigned j=0; double cfac=1;
      for (; j<CGR.dim2; ++j) {
         const CRef<TQ> &cj=CGR(0,j);
         if (!cj.isScalar()) wblog(FL,
            "ERR %s() invalid contraction to scalar\n%s",FCT,STR(cj));
         if (cj.wscalar()) { if (cj.cgw[0]!=1) wblog(FL,
            "ERR %s() invalid scalar normalization (%g)",FCT,cj.cgw[0]);
         }
         else if (cj.rtype!=CR_ABELIAN || cj.cgw) wblog(FL,
            "ERR %s() invalid scalar or empty\n%s",FCT,STR(cj)
         );
      }
      if (cfac<1e-12) wblog(FL,"ERR %s() got cfac=%g",FCT,cfac);
      (*DATA[0])*=cfac; CGR.init(1,0);
   }

   unsigned i,j, nz=0; char isz=0;
   double epsi=eps;

   if (!all) {
      if (rank(FL)!=2) { all|=2; } else
      if (!isBlockDiagMatrix('d') && !isHConj()) { all|=4; }
   }

   if (cgflag) { cgflag=(gotCGS(FL)>0); }

   if (isref) wblog(FL,"ERR must not change QSpace ref!");
   if (cgflag && (CGR.dim1!=DATA.len || CGR.dim2!=qtype.len)) wblog(FL,
      "ERR %s() CGR size mismatch (%dx%d ; %dx%d)",
       FCT, CGR.dim1, CGR.dim2, DATA.len, qtype.len
   );

   for (i=0; i<DATA.len; ++i) { wbarray<TD> &a=(*DATA[i]);

      epsi=eps; isz=0;

      if (cgflag) {
         double cfac, w2=1;
         for (j=0; j<CGR.dim2; ++j) { w2 *= double(CGR(i,j).norm2()); }
         cfac=sqrt(w2);

         if (cfac<eps) { isz=2;
            double x=a.aMax(); if (x>1e-6) { wblog(FL,
              "WRN %s() skipping |cgw|=%.3g having max(|data|)=%.3g !?",
              FCT,cfac,x
            ); }
         }
         else if (cfac>1) { epsi/=cfac; } 
      }

      if (!isz) {
         isz=(a.isZero(epsi,bflag) ? 1 : 0); 
      }

      if (isz)
           { if (all || isz>1) { a.initX(); ++nz; }}
      else { a.SkipTiny_imag(epsi); }
   }

   if (nz) { nz=SkipEmptyData(FL); }

   return nz;
};

template <class TQ, class TD> inline
unsigned QSpace<TQ,TD>::SkipEmptyData(const char *F, int L) {

   unsigned j, i=0, l=0, cgflag=(CGR.data ? 1:0);

   if (QIDX.dim1!=DATA.len || (cgflag && CGR.dim1!=DATA.len))
      wblog(FL,"ERR %s() unexpected QSpace (%d/%d/%d)",
      FCT,QIDX.dim1,DATA.len,CGR.dim1);
   if (isref) wblog(F_L,"ERR got QSpace reference");

   for (; i<DATA.len; ++i) {
      if (DATA[i]->isEmpty()) { 
         WB_DELETE_1(DATA[i]); if (cgflag) {
         for (j=0; j<CGR.dim2; ++j) { CGR(i,j).init(); }}
      }
      else {
         if (l<i) {
            QIDX.recSet(l,i);
            DATA[l]=DATA[i]; if (cgflag) {
            for (j=0; j<CGR.dim2; ++j) CGR(i,j).save2(CGR(l,j)); }
         }; ++l;
      }
   }

   if (!l) {
      DATA.init(); 
      QIDX.init(); if (cgflag) CGR.init();
      clearQSpace();
   }
   else if ((i=(DATA.len-l))) {
      QIDX.dim1=l; DATA.len=l; if (cgflag) CGR.dim1=l;
   }

   return i; 
};

template <class TQ, class TD> inline
bool QSpace<TQ,TD>::gotZeroData(unsigned i, const double eps) const {

   double x2, n2=1;

   if (i>=DATA.len) wblog(FL,
      "ERR %s() index out of bounds (%d/%d)",FCT,i,DATA.len);

   if (CGR.data) {
      if (i>=CGR.dim1) wblog(FL,
         "ERR %s() index out of bounds (%d/%dx%d)",FCT,i,CGR.dim1,CGR.dim2);
      for (unsigned j=0; j<CGR.dim2; ++j) { const CRef<TQ>& R=CGR(i,j);
         if (!R.isRefInit()) {
            n2*=(x2=double(R.norm2()));
            if (x2<1e-12) wblog(FL,"WRN %s() got small CGC (%g)",FCT,x2);
            if (!n2) return 1;
         }
      }
   }
   if (!DATA[i]) {
      wblog(FL,"WRN %s() got nullptr data[%d]",FCT,i);
      return 1;
   }
   n2 *= DATA[i]->norm2();

   return (std::sqrt(n2)<=eps);
};

template <class TQ, class TD> inline
unsigned QSpace<TQ,TD>::skipZeroOffDiag(double eps, char bflag){

   unsigned i,k,s, r=-1;
   wbindex I(DATA.len);
   TQ *qi;

   if (isEmpty()) { return 0; }
   isConsistent(FL,r);

   if (r%2) wblog(FL,
      "ERR %s() requires even-rank (%d)",FCT,r);

   s=QDIM*sizeof(TQ);

   for (k=i=0; i<DATA.len; ++i) { qi=QIDX.rec(i);
      if (memcmp(qi,qi+QDIM,s) && DATA[i]->isZero(eps,bflag)) {
         WB_DELETE_1(DATA[i]);
      }
      else {
         if (k<i) DATA[k]=DATA[i]; 
         I[k++]=i;
      }
   }

   i=DATA.len-k;

   if (k==DATA.len) return 0;
   if (k==0) {
      DATA.init(); 
      clearQSpace(); return i;
   }

   I.len=k; DATA.len=k;

   QIDX.Set2Recs(I);

   return i;    
}

template <class TQ, class TD>
void QSpace<TQ,TD>::TransferSpace(const char *F, int L,
   QSpace<TQ,TD> &B, unsigned dim, 
   const QSpace<TQ,double> &W, double eps
){
   unsigned i=0; int e=0;
   char gotpartial=0; 

   QSpace<TQ,TD> Aout, Bout, &A=*this;
   if ((void*)&A==(void*)&B) wblog(FL,"ERR got same space for A and B");

   if (this->isEmpty()) return;
   if (B.isEmpty()) { B.QDIM=A.QDIM; B.QIDX.init(0,A.QIDX.dim2); }

   if (A.QDIM!=B.QDIM || A.QIDX.dim2!=B.QIDX.dim2) {
      wblog(F_L,"ERR dimension mismatch (%d,%d; %d,%d)",
      A.QDIM, B.QDIM, A.QIDX.dim2, B.QIDX.dim2);
   }

   if (!W.isEmpty()) {
      for (i=0; i<W.DATA.len && gotpartial<3; ++i) {
         const wbarray<double> &w=(*W.DATA[i]);
         if (!w.isVector()) wblog(FL,
            "ERR %s() invalid weight space (%s)",FCT,w.sizeStr().data);
         for (unsigned n=w.numel(), j=0; j<n && gotpartial<3; ++j) {
            if (!(gotpartial & 1) && w.data[j]<eps) gotpartial|=1;
            if (!(gotpartial & 2) && w.data[j]>=eps) gotpartial|=2;
         }
      }
   }
   else gotpartial=1; 

   if (gotpartial==2) {
      return;
   }
   if (gotpartial==1) {
      TransferSpace(F_L,B,dim); return;
   }

   wbindex Ia,Ib,Ja,Jw,I2; --dim; 

   e=matchIndex(A.QIDX,B.QIDX,Ia,Ib);
   if (e) wblog(FL,"ERR %s() got non-unique QIDX !?",FCT);
   Ia.invert(QIDX.dim1,I2);

   wbMatrix<TQ> QA,QW;
   A.getQsub(dim,QA); 
   W.getQsub(1,QW);

   try {
      matchIndexU(FL,QA,QW,Ja,'f'); }
   catch (...) {
      wblog(F_L,"ERR %s() failed to match input weights uniquely",FCT);
   }

   if (QA.dim1!=Ja.len) wblog(FL,
      "ERR %s() failed to fully match QIDX with W",FCT);

   try {
      for (i=0; i<Ia.len; ++i) {
         DATA[Ia[i]]->Append2(FL,
          *B.DATA[Ib[i]], dim, *W.DATA[Ja[Ia[i]]], eps);
      }

      QSpace<TQ,TD> AX; getSubInit(I2,AX);
      for (i=0; i<I2.len; ++i) {
         DATA[I2[i]]->Append2(FL,
          *AX.DATA[i], dim, *W.DATA[Ja[I2[i]]], eps);

      }
      AX.Append2AndDestroy(FL,B);
   }
   catch (...) {
      MXPut(FL,"i").add(A,"A").add(B,"B").add(dim,"dim").add(W,"W")
        .add(QA,"QA").add(QW,"QW").add(Ia,"Ia").add(Ib,"Ib").add(Ja,"Ja");
      if (Ia.len && Ja.len)
         wblog(F_L,"ERR %s() %d/%d: %d,%d,%d",
         FCT, i+1,Ia.len, Ia[i]+1, Ib[i]+1, Ja[Ia[i]]+1);
      else wblog(F_L,"ERR %s() %d [%d %d]",FCT,i,Ia.len,Ja.len);
   }

   A.SkipEmptyData(F_L);
   B.SkipEmptyData(F_L);
};

template <class TQ, class TD>
void QSpace<TQ,TD>::TransferSpace(const char *F, int L,
   QSpace<TQ,TD> &B, unsigned dim
){
   unsigned i=0; int e=0, r=rank(FL);

   QSpace<TQ,TD> &A=*this;
   if ((void*)&A==(void*)&B) wblog(FL,"ERR got same space for A and B");

   if (A.QDIM!=B.QDIM || A.QIDX.dim2!=B.QIDX.dim2) {
      wblog(F_L,"ERR dimension mismatch (%d,%d; %d,%d)",
      A.QDIM, B.QDIM, A.QIDX.dim2, B.QIDX.dim2);
   }
   if (dim<1 || (int)dim>r) wblog(FL,
      "ERR %s() dimension out of bounds",FCT,dim,r);

   wbindex Ia,Ib,I2; --dim; 

   e=matchIndex(A.QIDX,B.QIDX,Ia,Ib);
   if (e) wblog(FL,"ERR %s() got non-unique QIDX !?",FCT);

   try {
      for (i=0; i<Ia.len; ++i) {
      DATA[Ia[i]]->Append2(*B.DATA[Ib[i]],dim); }
   }
   catch (...) { wblog(F_L,"ERR %s() %d/%d",FCT,i+1,Ia.len); }

   Ia.invert(QIDX.dim1,I2);
   if (I2.len) {
      QSpace<TQ,TD> AX;
      saveSub2(AX,I2); AX.Append2AndDestroy(F_L,B);
   }

   A.init();
};

template <class TQ, class TD>
QSpace<TQ,TD>& QSpace<TQ,TD>::TimesEl( 
   const QSpace<TQ,TD> &B, char conj,
   unsigned r2
){
   unsigned i,j,k; int i1,i2;

   if (isEmpty()) { return *this; }
   if (B.isEmpty()) { clearQSpace(); return *this; }
   checkQ(FL,B);

   unsigned r=rank(FL);
   wbMatrix<int> IQ;
   wbMatrix<TQ> QQ;
   wbstring mark(DATA.len);

   char cgflag=(CGR.isEmpty() ? 0 : 1);  
   if (cgflag && permitsOM() && (hasOM() || B.hasOM())) { cgflag|=2; }

   getQOverlapU(FL,B,QQ,IQ,'!'); 

   if (int(r2)>0 && r2!=r) {
      if (r2>r || r2<2) wblog(FL,
         "ERR %s() rank out of bound (r=%d/%d)",FCT,r2,r);
      if (r!=3 || r2!=2) wblog(FL,"TST %s() r = %d -> %d",FCT,r,r2);
   } else { r2=0; }

   for (i=0; i<IQ.dim1; ++i) { i1=IQ(i,0); i2=IQ(i,1);
      if (i1<0 || i2<0) wblog(FL,"ERR got (i1,i2)=(%d,%d)",i1,i2);
      if ((++mark[i1])>1) wblog(FL, 
         "ERR %s() got mark[%d]=%d", FCT, i1, mark[i1]);

      if (cgflag<=1) {
         if (r2)              
              { DATA[i1]->timesEl_OM(r2, *B.DATA[i2], *DATA[i1], conj); }
         else { DATA[i1]->TimesEl(*B.DATA[i2],conj); }

         if (cgflag) { double w2=1;
            for (j=0; j<CGR.dim2; ++j) {
               w2 *= CGR(i1,j).wget0() * B.CGR(i2,j).wget0();
               CGR(i1,j).initCtrScalar(); 
            }
            (*DATA[i1]) *= w2;
         }
      }
      else {
         unsigned m=CGR.dim2; double w2=1;
         wbvector<unsigned> Sa(m), Sb(m);
         wbarray<TD> Ai(*DATA[i1],'r'); wbarray<double> x2;

         for (j=0; j<B.CGR.dim2; ++j) {
            Sa[j]=  CGR(i1,j).wdim2(FL); 
            Sb[j]=B.CGR(i2,j).wdim2(FL);
         }

         Ai.ExpandOM(FL,r,Sa); 

         for (j=0; j<B.CGR.dim2; ++j) {
            CGR(i1,j).wProd(B.CGR(i2,j),x2);

            if (x2.numel()==1)
                 { w2*=x2[0]; }
            else { Ai.ContractMat(FL,r+j+1,x2); } 

            CGR(i1,j).initCtrScalar(); 
         }
         Ai.FuseOM(FL,r, Sb.prod(1));

         Ai.timesEl_OM(r2 ? r2 : r, *B.DATA[i2], *DATA[i1], conj);

         if (w2!=1) {
            if (fabs(w2)<CG_SKIP_DEPS1) wblog(FL,"WRN %s() got w2=%g",FCT,w2);
            Ai.init();         
            (*DATA[i1]) *= w2; 
         }
      }

      if (r2) { DATA[i1]->appendSingletons(r); }
   }

   for (k=i=0; i<DATA.len; ++i) {
      if (mark[i]) { if (k!=i) { recSave2(k,i); }; ++k; }
   }

   if (k && k<DATA.len) {
      QIDX.dim1=DATA.len=k; 

      unsigned n=CGR.numel();
      for (i=k*CGR.dim2; i<n; ++i) { CGR[i].init(); }
      CGR.dim1=k;
   }
   else if (!k) { clearQSpace(); }

   return *this;
};

template <class TQ, class TD>
QSpace<TQ,TD>& QSpace<TQ,TD>::plus_plain(
   const QSpace &B, QSpace &C, TD bfac) const {

   if (!QIDX || !B.QIDX) {
      if (B.isEmpty()) { C=*this;      return C; }
      if (  isEmpty()) { C=B; C*=bfac; return C; }
      if (!isScalar() || !B.isScalar()) {
         this->info("A"); B.info("B");
         wblog(FL,"ERR %s() invalid empty QSpaces !?",FCT);
      }
      C=*this; C.DATA[0]->data[0] += (bfac*B.DATA[0]->data[0]);
      return C;
   }

   unsigned i,j, n1=0, n2=0, n12=0; int i1,i2, np=1;
   unsigned cgflag=0; 
   wbvector<double> bfc;

   wbMatrix<int> IQ;
   wbMatrix<TQ> QQ;

   if (qtype.len && B.qtype.len) {
      if (qtype!=B.qtype) wblog(FL,"ERR %s() "
         "got inconsistent qtype %s/%s !?",FCT,qStr().data,B.qStr().data);
      if (!CGR.isEmpty() && !B.CGR.isEmpty()) cgflag=1;
   }
   if (!cgflag) {
      if (  !allAbelian() ||   !CGR.isEmpty() || 
          !B.allAbelian() || !B.CGR.isEmpty() || (qtype.len!=B.qtype.len))
      wblog(FL,"ERR %s() got qtype '%s' vs. '%s' !?",
      FCT,qStr().data,B.qStr().data);
   }

   getQOverlapU(FL,B,QQ,IQ);

   C.init(QQ, isEmpty() ? B.QDIM : QDIM, nullptr);
   C.initQ(FL,*this,1,&B);

   if (itags.len || B.itags.len) {
      if (!B.itags.len) C.itags=  itags; else
      if (!  itags.len) C.itags=B.itags;
      else {
         if (!itags.sameConj(B.itags)) wblog(FL,
            "ERR %s() got inconsistent itags\n%s <> %s",
            FCT,IT2STR__,IT2STR(B));
         C.itags=itags;
         if (itags!=B.itags) { C.itags.init_tags(); }
      }
   }

   checkQ(FL,B, cgflag? 0:'a');

   if (cgflag) { C.setupCGR(); str[0]=0; j=CGR.dim2;
      if (B.CGR.dim2!=j || C.CGR.dim2!=j)
         sprintf_str("CGR.dim2: %d/%ld/%ld",j,B.CGR.dim2,C.CGR.dim2); else
      if (CGR.dim1 && CGR.dim1!=QIDX.dim1)
         sprintf_str("A.CGR.dim1: %ld/%ld", CGR.dim1, QIDX.dim1); else
      if (B.CGR.dim1!=B.QIDX.dim1)
         sprintf_str("B.CGR.dim1: %ld/%ld",B.CGR.dim1,B.QIDX.dim1); else
      if (C.CGR.dim1!=C.QIDX.dim1)
         sprintf_str("C.CGR.dim1: %ld/%ld",C.CGR.dim1,C.QIDX.dim1);
      if (str[0]) wblog(FL,"ERR plus() CGR inconsistency (%s)",str);

      bfc.init2val(B.DATA.len,bfac);

      for (i=0; i<IQ.dim1; ++i) { i1=IQ(i,0), i2=IQ(i,1);
         if (i1>=0 && i2>=0) { 
            for (j=0; j<CGR.dim2; ++j)
            bfc[i2]*=double(CGR(i1,j).safeCpy(FL, B.CGR(i2,j), C.CGR(i,j)));
         }
         else if (i1>=0) {    
            for (j=0; j<CGR.dim2; ++j) C.CGR(i,j)=CGR(i1,j);
         }
         else if (i2>=0) {    
            for (j=0; j<CGR.dim2; ++j) C.CGR(i,j)=B.CGR(i2,j);
         }
         else wblog(FL,"ERR i1=%d, i2=%d !?",i1,i2);
      }
   }

   if (IQ.dim1>1 && !omp_in_parallel()) { 
      np=MAX( QSP_NUM_THREADS, OMP_NUM_THREADS );
      if (np>int(IQ.dim1)) { np=IQ.dim1; } else if (np<1) { np=1; }
   }

   Wb::LogException ex; 

  #pragma omp parallel for num_threads(np) reduction(+:n1,n2,n12)
   for (unsigned i=0; i<IQ.dim1; ++i) { if (!ex) { try {
      int i1=IQ(i,0), i2=IQ(i,1); 

      if (i1>=0 && i2>=0) { ++n12; 
         DATA[i1]->plus(*B.DATA[i2], *C.DATA[i], cgflag ? TD(bfc[i2]) : bfac);
      }
      else if (i1>=0) { ++n1; 
         C.DATA[i]->init(*DATA[i1]);
      }
      else if (i2>=0) { ++n2; 
         C.DATA[i]->set(*B.DATA[i2], bfac);
      }
      else wblog(FL,"ERR i1=%d, i2=%d !?",i1,i2);
   }
      catch (Wb::LogException &e_) { ex+=e_; }
      catch (...) { ++ex; }
   }}

   if (ex) { 
      C.clearQSpace('!'); 
      ex.report(FLF);
   }

   return C;
};

template <class TQ, class TD>
int QSpace<TQ,TD>::isNormCGW(char full) const {

   int rval=1;

   if (CGR) {
      unsigned i,j, nwrn=0;
      double cfac;

      for (i=0; i<CGR.dim1 && !nwrn; ++i) {
      for (j=0; j<CGR.dim2; ++j) {
         cfac = CGR(i,j).getNormSignW(0,0, full? 0:1);
         if (fabs(cfac-1)>CG_SKIP_DEPS1) { ++nwrn; }
      }}
      rval=(nwrn? 0 : 2);
   }

   return rval;
};

template <class TQ, class TD>
int QSpace<TQ,TD>::NormCGW(
    char full,
    char skipzeros, 
    char rcpy   
) {

   int rval=0; if (!CGR) { return rval; }
   unsigned i,j,m;
   double cfac;

   wbvector<unsigned> mark;

   m=gotCGS(FL,'x'); rval=m; if (int(m)<0) { m=0; }
   unRef(); 

   if (!full || m<=1) { 

      for (i=0; i<CGR.dim1; ++i) { cfac=1;
         for (j=0; j<CGR.dim2; ++j) {
            cfac *= CGR(i,j).NormSignW(0,0, full? 0:1);
         }

         if (cfac!=TD(1)) {
            if (fabs(cfac)>CG_SKIP_DEPS2 || !skipzeros) {
               DATA[i]->Times(cfac,rcpy); 
            }
            else {
               if (!mark.len) { mark.init(CGR.dim1); }
               mark[i]=1;
            }
         }
      }
   }
   else { 
      unsigned r=rank(FL); double cfacM;
      wbvector<unsigned> iOM;

      ExpandOM(FL,iOM,'f'); 

      for (i=0; i<CGR.dim1; ++i) { cfac=cfacM=1;
         for (j=0; j<CGR.dim2; ++j) { CRef<TQ> &Rij=CGR(i,j);
            if (!iOM[j] || Rij.wisId()>0) {
               cfac *= Rij.NormSignW(0,0, full? 0:1);
            }
            else {
               m=Rij.wdim1(); 
               cfacM *= Rij.cgw.norm(); 
               DATA[i]->ContractMat(FL,iOM[j]+1,Rij.cgw,2); 
               Rij.cgw.initIdentity(m);
            }
         }
         if (fabs(cfac*cfacM)>CG_SKIP_DEPS2 || !skipzeros) {
            DATA[i]->FuseOM(FL,r);
            if (cfac!=TD(1)) {
               DATA[i]->Times(cfac,rcpy); 
            }
         }
         else {
            if (!mark.len) { mark.init(CGR.dim1); }
            mark[i]=1;
         }
      }
   }

   if (mark) {
      wbindex I; mark.find(I,0);
      rval+=8*(DATA.len-I.len);
      Select(I);
   }

   if ((i=qtype.allAbelian())) {
      checkQ_CGR(FL); CGR.init();
      if (i==1) { qtype.ReduceU1(); } 
   }

   return rval;
};

template <class TQ, class TD>
unsigned QSpace<TQ,TD>::MakeUnique() {

   if (isEmpty()) { return 0; }

   wbMatrix<TQ> QS(QIDX);
   wbvector<widx_t> D;
   wbperm P;

   QS.groupRecs(P,D); if (D.allLE(1)) { return 0; } 

   unsigned i0,i2,ib,j,j0,k,l, mx=0, Di, r=rank(FL);
   double cfac0=1, cfac2=1;

   wperm_t const* p; char noM=(!permitsOM());
   wbvector<widx_t> I;

   if (QIDX.dim1!=DATA.len || !QDIM || QIDX.dim2%QDIM) wblog(FL,
      "ERR QSpace size mismatch (%d/%d,%d)",DATA.len,QIDX.dim1,QDIM);

   NormCGW( noM? 0 : 'f'); 

   if (QS.dim1!=QIDX.dim1)
        { QIDX.groupRecs(P,D); } 
   else { QS.save2(QIDX); }

   p=P.data; I.init(DATA.len);

   for (Di=k=ib=0; ib<D.len; ++ib, k+=Di) { Di=D[ib];
      wbvector<char> mark(Di,1); 
      for (j0=Di, j=0; j<Di; ++j) {
         if (gotZeroData(p[k+j])) { mark[j]=0; } else
         if (j0>j) { j0=j; } 
      }; if (j0==Di) { j0=0; mark[0]=-1; }

      i0=I[ib]=p[k+j0];
      wbarray<TD> &a0=(*DATA[i0]);
      CRef<TQ> *r0=nullptr, *r2=nullptr;

      if (CGR.data) { r0=CGR.ref(i0);
         if (noM) { cfac0=1;
            for (l=0; l<CGR.dim2; ++l) { cfac0 *= r0[l].NormSignW(); }
            if (fabs(cfac0)<1e-6) { wblog(FL,
               "ERR %s() got cfac = %g",FCT,cfac0);
            }
         }
      }

      for (j=0; j<Di; ++j) { if (j!=j0) {
         i2=p[k+j]; 
         wbarray<TD> &a2=(*DATA[i2]);

         if (r0) { char e; r2=CGR.ref(i2);
            if (noM) { cfac2=1;
               for (l=0; l<CGR.dim2; ++l) { cfac2 *= r2[l].NormSignW(); }
            }

            for (l=0; l<CGR.dim2; ++l) {
               if ((e=r0[l].sameAs_fix(0,0,r2[l])) && (e<11 || e>12)) {
                  MXPut(FL,"Iq","base")
                    .add(*this,"A").add(a0,"a0").add(a2,"a2")
                    .add(i0+1,"i1").add(i2+1,"i2").add(l+1,"j")
                    .add(r0[l],"r0").add(r2[l],"r2");
                  wblog(FL,"WRN %s() got CGR([%d %d], %d) mismatch",
                    FCT,i0+1,i2+1,l+1);
                  r0[l].sameAs_fix(FL,r2[l]); 
               }
            }
         }

         if (a0.numOM(r)>1 || a2.numOM(r)>1) {
            ExpandOM_(FL,i0,i2); 
         }
         if (!a2.sameSize(a0,0)) { 
            MXPut(FL,"Iq","base").add(*this,"A")
              .add(i0+1,"i1").add(i2+1,"i2");
            wblog(FL,"ERR size inconsistency (b=%d: %d,%d): %s / %s",
            ib+1, i0+1, i2+1, SSTR(a2), SSTR(a0));
         }

         if (mark[j]) {
            if (j<=j0) { wblog(FL,"ERR %s() %d/%d",FCT,j,j0); }
            if (cfac0!=1 || cfac2!=1) {
               if (fabs(cfac2)<1e-8) wblog(FL,
                  "ERR %s() got cfac = %g",FCT,cfac2);
               a0.Plus(a2,TD(cfac2),0,TD(cfac0)); cfac0=1;
            }
            else { a0+=a2; }
         }
         a2.init(); ++mx; 

         if (r0) { for (l=0; l<CGR.dim2; ++l) { r2[l].init(); }}

      }} 
   } 

   I.Shorten2(ib);

   this->Select((wbindex&)I,'d');
   return mx;
};

template <class TQ, class TD> inline
TD QSpace<TQ,TD>::trace() const {

    unsigned i,j,K,s, r=-1; gotCGS(FL);
    TQ *q=QIDX.data;
    TD tsum=0, t=0;

    if (isEmpty()) { return 0; }
    isConsistent(FL,r);

    if (r%2) wblog(FL,"ERR %s() requires even-rank (%d)",FCT,r);

    K=QDIM*(unsigned)r/2; s=K*sizeof(TQ);

    for (i=0; i<DATA.len; ++i, q+=QIDX.dim2)
    if (!memcmp(q, q+K, s)) {
       if (DATA[i]->SIZE.len==r) {
          t=DATA[i]->trace();
          for (j=0; j<CGR.dim2; ++j) { t*=TD(CGR(i,j).trace()); }
          tsum+=t;
       }
       else {
          wbvector<double> tc(1,1.);
          wbvector<TD> tx; DATA[i]->trace(r,tx);

          for (j=0; j<CGR.dim2; ++j) {
             tc.TensorProd(CGR(i,j).trace(FL)); 
          }
          if (tx.len!=tc.len) wblog(FL,
             "ERR %s() OM length mismatch (%d/%d)",FCT,tx.len,tc.len);
          tsum+=tx.dotProd(tc);
       }
    }

    return tsum;
};

template <class TQ, class TD>
template <class TM>
wbMatrix<TM>& QSpace<TQ,TD>::diag(wbMatrix<TM> &xd) const {

   unsigned r,j,d, i=0, n=DATA.len;
   wbvector<unsigned> dd(n);
   char q,dflag=-1; 
   TM *x, m;

   if (isEmpty()) { return xd.init(); }

   r=rank();
   if (r!=2) wblog(FL,
      "ERR %s() got invalid rank-%d QSpace !?",FCT,r);

   for (; i<n; ++i) { q=DATA[i]->isVector();
      if (!q && !DATA[i]->isSMatrix()) wblog(FL,
         "ERR %s() got non-square data (%d: %s) !?",FCT,i+1,SSTR_(DATA[i]));

      if (dflag>=0) {
         if (q>=0 && dflag!=q) { wblog(FL,"ERR %s() "
            "mismatch in diagonal representation (%d: %d/%d; %s) !?",
            FCT,i+1,q,dflag,SSTR_(DATA[i]));
         }
      }
      else { dflag=q; }

      dd[i] = (!dflag ? DATA[i]->dim2() : DATA[i]->numel());
   }

   xd.init(dd.sum(),2); x=xd.data;

   for (i=0; i<n; ++i) {
      const TD *a=DATA[i]->data; d=dd[i];
      m=TM(cgsDim(i,0));
      for (j=0; j<d; ++j, x+=2) {
         x[0]=TM(dflag<=0 ? a[j+d*j] : a[j]);
         x[1]=m;
      }
   }

   return xd;
};

template <class TQ, class TD>
double QSpace<TQ,TD>::maxDiff(const QSpace<TQ,TD> &B) const {

   unsigned i,m=0; int i1, i2;
   wbMatrix<int> IQ;
   wbMatrix<TQ> QQ;
   double dmax=0, d=0;

   if (QDIM!=B.QDIM || QIDX.dim2!=B.QIDX.dim2) {
      sprintf_str("%s:%d QDIM mismatch (%d,%d; %ld,%ld)", FL,
      QDIM, B.QDIM, QIDX.dim2, B.QIDX.dim2); return NAN;
   }

   getQOverlapU(FL,B,QQ,IQ);

   for (i=0; i<IQ.dim1; ++i) {
       i1=IQ(i,0); i2=IQ(i,1);

       if (i1>=0 && i2>=0) {
          d=DATA[i1]->maxDiff(*B.DATA[i2]); ++m;
          if (Wb::isnan(d)) return d;
       }
       else {
          if (i1>=0) d=  DATA[i1]->aMax(); else
          if (i2>=0) d=B.DATA[i2]->aMax(); else
          wblog(FL,"ERR i1=%d, i2=%d !?", i1, i2);
       }
       dmax=MAX(dmax,d);
   }

   if (!m) wblog(FL, "WRN no Q-overlap");

   return dmax;
}

template <class TQ, class TD>
TD QSpace<TQ,TD>::normDiff2(const QSpace<TQ,TD> &B) const {

   unsigned i;
   wbMatrix<int> IQ;
   wbMatrix<TQ> QQ;
   TD x=0;

   if (isEmpty() || B.isEmpty()) {
      wblog(FL,"WRN %s() with empty object (%d,%d) !?",
      FCT, isEmpty(), B.isEmpty()); return x;
   }

   if (this==&B) return x;

   if (QIDX!=B.QIDX) {
      int i1,i2;

      getQOverlapU(FL,B,QQ,IQ);

      for (i=0; i<IQ.dim1; ++i) {
          i1=IQ(i,0); i2=IQ(i,1);

          if (i1>=0 && i2>=0)
             x+=(DATA[i1]->normDiff2(*B.DATA[i2]));
          else {
             if (i1>=0) x+=(  DATA[i1]->norm2()); else
             if (i2>=0) x+=(B.DATA[i2]->norm2());
             else wblog(FL,"ERR i1=%d, i2=%d !?", i1, i2);
          }
      }
   }
   else {
      for (i=0; i<DATA.len; ++i)
      x+=DATA[i]->normDiff2(*B.DATA[i]);
   }

   return x;
};

template <class TQ, class TD>
TD QSpace<TQ,TD>::sumData() const {

    TD x=0;
    for (unsigned i=0; i<DATA.len; ++i) x+=(DATA[i]->sum());
    return x;
}

template <class TQ, class TD> inline
QSpace<TQ,TD>& QSpace<TQ,TD>::revertLeg(
   const char *F, int L, unsigned i1, 
   QSpace<TQ,TD> &C) const {

   unsigned r=rank(FL);
   if (i1>=r) wblog(FL,"ERR %s() index out of bounds (i1=%d/%d)",FCT,i1,r);

   ctrIdx I1(1,&i1), I2(0, getQDir(i1)>0 ? '*':' ');
   QSpace<TQ,TD> E;

   E.initIdentityCG( cPVEC1_(*this),I1,'z');
   contractMat(I1,E,I2,C); 

   return C;
};

template <class TQ, class TD> inline
QSpace<TQ,TD>& QSpace<TQ,TD>::contractMat(
   unsigned ia, const QSpace<TQ,TD> &B, unsigned ib, 
   QSpace<TQ,TD> &C
 ) const {

   unsigned ra=rank(FL), rb=B.rank(FL);

   if (!ia || !ib) wblog(FL,
      "ERR %s() invalid ctr-index (%d,%d; using 1-based)",FCT,ia,ib);
   if (ia>ra || ib>rb) wblog(FL,
      "ERR %s() index out of range (A: %d/%d; B: %d/%d)",FCT,ia,ra,ib,rb);
   if (rb!=2) wblog(FL,"ERR %s() requires rank-2 tensor B (r=%d)",FCT,rb);

   ctrIdx ica(1,&(--ia)), icb(1,&(--ib)); 
   wbperm P; P.initLastTo(ia,ra);

   contract(ica,B,icb,C,P);

   return C;
};

template <class TQ, class TD> inline
QSpace<TQ,TD>& QSpace<TQ,TD>::contractMat(
   const ctrIdx &ica, const QSpace<TQ,TD> &B, const ctrIdx &icb,
   QSpace<TQ,TD> &C) const {

   unsigned ra=rank(FL), rb=B.rank(FL);
   wbperm P;

   if (rb!=2) wblog(FL,"ERR %s() requires rank-2 tensor B (r=%d)",FCT,rb);
   if (ica.len!=1 || icb.len!=1 || ica[0]>=ra || icb[0]>=rb) { wblog(FL,
      "ERR %s() invalid ctr-index (A: %s/%d; B: %s/%d)",
      FCT,STR(ica),ra,STR(icb),rb);
   }

   P.initLastTo(ica[0],ra);
   contract(ica,B,icb,C,P);

   return C;
};

template <class TQ, class TD>
template <class TB, class TC>
int QSpace<TQ,TD>::contract_getIdxSet(const char *F, int L,
   const ctrIdx &ica, 
   const QSpace<TQ,TB> &B, const ctrIdx &icb,
   wbindex &Ia, wbindex &Ib, wbvector<widx_t> &Dc,
   QSpace<TQ,TC> &C,   
   MVEC &fA, MVEC &fB  
 ) const {

   unsigned i=0, l, ra=rank(F_L), rb=B.rank(F_L); int isf, e;

   const wbindex ica_(ica), icb_(icb);
   wbindex ika, ikb;
   wbMatrix<TQ> QAc, QAk, QBc, QBk;

   wbMatrix<TQ> &QC=C.QIDX;
   iTags &idC=C.itags;
   iFerm &fdC=C.fdir;

   wbMatrix<unsigned> Z2;

#ifdef WB_CLOCK
   Wb::Clock clk("QS:ctr:getIdx",1); 
#endif

   if (QDIM!=B.QDIM) wblog(F_L,"ERR %scontract() "
      "incompatible objects (QDIM=%d/%d)",F? "QSpace::":"",QDIM,B.QDIM);
   if (ica.len!=icb.len || !ica.len) wblog(F_L,
      "ERR invalid contraction [%s] <> [%s]",STR(ica), STR(icb));

   idC.init(); Ia.init(); Ib.init();

   if (!isUniqueIdxSet(ica,ra) || !isUniqueIdxSet(icb,rb) ||
       ica.len!=icb.len) wblog(F_L,"ERR "
      "invalid contraction indices\n[%s; %d], [%s; %d] (QDIM=%d/%d)",
       STR(ica), ra, STR(icb), rb, QDIM, B.QDIM
   );

   ica_.invert(ra,ika);
   icb_.invert(rb,ikb);

   if ((isf=isFerm(0,0,&B.fdir))<0) wblog(FL,"ERR %s() "
      "fdir inconsistency (%s / %s; e=%d)",FCT,STR(fdir),STR(B.fdir),isf);

   fdC.init();
   if (isf) {
      char a, b; unsigned l=0;
      wbindex iz(ica.len);
      wbperm pA,pB;

if (Wb::envFERM) {
wblog(FL,"TST %40R","-");
wblog(FL,"TST ica='%s' ;  icb='%s'",STR(ica),STR(icb));
}
      pA.init2End  (ica,ra); if (ica.conj) { pA.conj=1; }
      pB.init2Front(icb,rb); if (icb.conj) { pB.conj=1; }
if (Wb::envFERM) {
wblog(FL,"--> ica='%s' ;  icb='%s' => pA=%s ;  pB=%s",STR(ica),STR(icb),STR(pA),STR(pB));
}
      for (i=0; i<ica.len; ++i) {
         a=  fdir[ica[i]]; if (ica.conj) { a=-a; }
         b=B.fdir[icb[i]]; if (icb.conj) { b=-b; }
         if (!a || a!=-b) { wblog(FL,
            "ERR %s() iFerm direction mismatch\n  %-12s @ %s\n  %-12s @ %s",
            FCT, STR(fdir), STR(ica), STR(B.fdir), STR(icb));
         }
         if (a<0) {
            iz[l++] = ica[i]; 
         }
      }

      if (!l) { iz.init(); } else { iz.len=l; }

        getFermSigns(pA,fA,&ica,&iz);  
      B.getFermSigns(pB,fB);

wblog(FL,"TST %s() %d/%d %d/%d",FCT,pA.len,DATA.len,pB.len,B.DATA.len);

if (Wb::envFERM) {
   static unsigned i=0; snprintf(str,8,"Ix%d",++i);
   MXPut(FL,str,"base").add(pA,"pA").add(pB,"pB")
    .add(fA,"fA") .add(fB,"fB").add(iz,"iz");
   wblog(FL,"TST ica='%s' ;  fA=%s ",STR(ica),STR(fA));
   wblog(FL,"TST icb='%s' ;  fB=%s ",STR(icb),STR(fB));
}

      fdC.init(ika.len+ikb.len, fdir); 
      a=(ica.conj ? -1 : 1); 
      b=(icb.conj ? -1 : 1);

      for (     i=0; i<ika.len; ++i) { fdC[i  ] = a*  fdir.at(ika[i]); }
      for (l=i, i=0; i<ikb.len; ++i) { fdC[i+l] = b*B.fdir.at(ikb[i]); }
if (Wb::envFERM) {
wblog(FL,"TST %s | %s | %s",STR(fdir), STR(B.fdir), STR(fdC)); }

if (!fdC) wblog(FL,"ERR %s() ",FCT);
   }

   iTags ta, tb;

   if (itags.len) { ta=itags; }     
   if (ica.newtags.len) { ta.Update(ica.newtags); }

   if (B.itags.len) { tb=B.itags; } 
   if (icb.newtags.len) { tb.Update(icb.newtags); }

   if (ta.len!=ra && (CGR.data || ta.len)) { wblog(FL,
      "ERR ctr/getIdxSet() %s itags (A: len=%d/%d; %ld)",
      ta.len? "invalid":"missing", ta.len,ra, CGR.numel()); }
   if (tb.len!=rb && (B.CGR.data || tb.len)) { wblog(FL,
      "ERR ctr/getIdxSet() %s itags (B: len=%d/%d; %ld)",
      tb.len? "invalid":"missing", tb.len,rb, B.CGR.numel()); }

   idC.init(ika.len+ikb.len);
   for (     i=0; i<ika.len; ++i) { idC[i  ]=ta[ika[i]]; }
   for (l=i, i=0; i<ikb.len; ++i) { idC[i+l]=tb[ikb[i]]; }

   if (ta.len && tb.len) {  char q=0, zflag=0;
      for (i=0; i<ica.len; ++i)  {
         if ((q=ta[ica[i]].sameAs(tb[icb[i]],'l'))) { q=SGN(q);
            if (zflag) { if (zflag!=q) q=0; }
            else { zflag=q; }
         }
         if (!q) wblog(FL,"ERR contract() itag mismatch\n"
           "%s [%s] <> %s [%s] @ i=%d %N",
            IT2STR__, STR(ica+1),
            IT2STR(B),STR(icb+1), i+1
         );
      }

      if (zflag>0) {
         if (!(ica.conj^icb.conj)) { 
            wblog(FL,"WRN contract() check missing conj-flag!?");
            mexWRN("applying conj(A) flag");
            ((ctrIdx&)ica).Conj(); 
         }
      }
      if (ica.conj) { for (i=0;       i<ika.len; ++i) idC[i].Conj(); }
      if (icb.conj) { for (i=ika.len; i<idC.len; ++i) idC[i].Conj(); }
   }

     getQsub(ica_,QAc,QAk);
   B.getQsub(icb_,QBc,QBk);

   e=matchIndex(QAc,QBc,Ia,Ib);

   QC.Cat(2,QAk,Ia,QBk,Ib); 

   if (Ia.isEmpty()) {
      Dc.init(); 
   }
   else if (!ika.isEmpty() || !ikb.isEmpty()) { 
      wbperm P; QC.groupRecs(P,Dc);
      Ia.Permute(P); Ib.Permute(P);

   }
   else {
      Dc.init(1); Dc[0]=Ia.len;
      QC.init(1,0); 
   }

   return e;
};

template <class TQ, class TD>
template <class TB, class TC>
double QSpace<TQ,TD>::contract(const char *F, int L, 
   ctrIdx ica, const QSpace<TQ,TB> &B,
   ctrIdx icb, QSpace<TQ,TC> &C, const wbperm &P,
   char preview 
) const {

   if ((void*)this==(void*)&C || (void*)&B==(void*)&C) {
      QSpace<TQ,TC> X;
      int i=contract(F,L,ica,B,icb,X,P,preview); X.save2(C);
      return i;
   }

   C.clearQSpace(); 
   if (isEmpty() || B.isEmpty()) { return 1; }

#ifdef WB_CLOCK
   Wb::Clock clk("QS:ctr:actual",1); 
#endif

#ifdef LD_CLEBSCH_QS
#endif

   unsigned nc, Nx=0; int np=1;
   unsigned ra=rank(F_L), rb=B.rank(F_L), rc=(ra+rb)-ica.len-icb.len;
   char cgflag=gotCGS(FL);

   QSpace<TQ,TD> A2;
   QSpace<TQ,TB> B2;

   wbvector<widx_t> D,Dc;
   wbindex Ia,Ib;
   MVEC fA, fB;

   Wb::LogException ex;

   if (QDIM!=B.QDIM) wblog(F_L,
      "ERR QSpace() QDIM mismatch (%d/%d)",QDIM,B.QDIM);
   if (!P.isEmpty() && (P.isValidPerm()<=0 || P.len>rc)) wblog(FL,
      "ERR %s() invalid permutation [%s; %d]",FCT,STR(P),rc);

   if (qtype!=B.qtype) {
      if (!allAbelian() || !B.allAbelian() || (qtype.len && B.qtype.len))
      wblog(FL,"ERR %s() qtype inconsistency '%s' / '%s'",
      FCT, qStr().data, B.qStr().data);
   }

   if (!ica.isSorted() && !icb.isSorted()) {
      size_t sA=getDataSize(), sB=B.getDataSize();
      if (sA>sB)
           { ica.Sort_(icb); }
      else { icb.Sort_(ica); }
   }

     permute_to(2,ica,A2,'r'); 
   B.permute_to(1,icb,B2,'r'); 

   A2.contract_getIdxSet(FL,ica,B2,icb, Ia,Ib,D, C, fA, fB);

   C.initQ(F,L,A2,0, cgflag<=0? nullptr : &B2);
   C.setupDATA(); if (cgflag>0) {
   C.setupCGR();  }

   D.cumsum_(Dc,'x');

   nc=D.len;
   if (nc>1 && QSP_NUM_THREADS>1 && !omp_in_parallel()) {
      np=MIN( QSP_NUM_THREADS, int(nc) );
   }

   if (WBLOG_CTR || (WBLOG_RLARGE && Dc.last()>(1<<14))) {
      size_t Dsum=Dc.last();       
      wblog(FL,"CTR [A @ %s]*[B @ %s] %ld entries (avg. %.1f per QC block)",
      STR(ica),STR(icb),Dsum,Dsum/double(D.len));
   }

   itag_::Reset();

   #pragma omp parallel for num_threads(np)
    for (unsigned ic=0; ic<nc; ++ic) { if (!ex) {
     # if defined(DBSTOP) && defined(QS_USING_OMP) && defined(__OMP_H)
       char tid[8]; 
       snprintf(tid,8,"%d/%d",omp_get_thread_num(),omp_get_num_threads());
     # endif

       widx_t i0=Dc[ic], Di=D[ic];
       const wbvector<widx_t>
          Ia_(Di,Ia.data+i0,'r'),
          Ib_(Di,Ib.data+i0,'r');

       try { 
          Nx+=contractDATA_group(FL, 
              A2,fA,Ia_,ica, B2,fB,Ib_,icb, C,ic,preview);
       }
       catch (Wb::LogException &e_) { ex+=e_; }
       catch (...) { ++ex; }
    }}
    ex.report(FLF);

    if (preview) { return Nx; }

    if (WBLOG_CTR) wblog(FL,
       "--> contractDATA done (nq=%d/%d threads, %sabelian)",
       np,omp_get_num_threads(), cgflag>0 ? "non-" : cgflag<0 ? "[non-]":""
    );

   if (Nx)
        { C.SkipEmptyData(FL); }
   else { C.SkipZeroData(); }

   C.checkScalarCGS(FL);
   C.Permute(P);

   if (C.rank(FL)<=2) {
      C.NormCGW(); 
   }

   if (WBLOG_CTR) { 
      wblog(FL,"     %N");
      A2.info("A"); B2.info("B"); C.info("C"); printf("\n");
   }

   return 0; 
};

template <class TQ, class TD>
QSpace<TQ,TD>& QSpace<TQ,TD>::trace(
   const char *F, int L, ctrIdx &i1, ctrIdx &i2, QSpace<TQ,TD> &B
 ) const {

   unsigned i,j,d,k=0,l,m,n, r=rank(F_L), rB;
   int e=0; 

   MVEC mark(r);

   if (itags.len!=r) wblog(F_L,
      "ERR %s() itag length mismatch ('%s' / %d)",FCT,STR(itags),r);

   if (!i1.len && !i2.len) {
      i1.init(r); i2.init(r); l=0;
      for (i=0; i<r; ++i) { m=0;
         for (j=i+1; j<r; ++j) { if (itags[i].isConj(itags[j])) { ++m; k=j; }}
         if (m>1) wblog(FL,"ERR %s() "
            "auto-pairing of indices not unique (m=%d)",FCT,m); else
         if (m==1) { i1[l]=i; i2[l]=k; ++l; }
      }
      if (l) { i1.len=i2.len=l; }
      else { i1.init(); i2.init(); }
   }

   if (!i1.len || i1.len!=i2.len) { e|=1; } else
   if (i1.conj ^ i2.conj) { e|=2; } 
   else {
      for (i=0; i<i1.len; ++i) {
         j=i1[i]; if (j>=r || ++mark[j]>1) { e|= 4; }
         k=i2[i]; if (k>=r || ++mark[k]>1) { e|= 8; }
         if (!itags[j].isConj(itags[k]))   { e|=16; }
      }
   }
   if (e&15) wblog(F_L,
     "ERR %s() invalid indices i1=%s / i2=%s (r=%d, e=%d)",
      FCT,STR(i1),STR(i2),r,e);
   else if (e) wblog(F_L,
     "ERR %s() invalid indices (itag or conj mismatch)\n"
     "having '%s' @ %s <> %s !?",FCT,STR(itags),STR(i1),STR(i2)
   );

   rB=r-2*i1.len;

   if (!rB) {
      if (i1[0]<i2[0])
           { for (i=0; i<i1.len; ++i) { if (i2[i]!=i1[i]+i1.len) break; }}
      else { for (i=0; i<i1.len; ++i) { if (i1[i]!=i2[i]+i2.len) break; }}

      if (i==i1.len) { return B.initScalar(trace()); }
   }

   wbindex D, ik, ic(2*i2.len);
   wbMatrix<TQ> Qc, Qk;
   wbperm P;
   wbarray<TD> x2;

   ic.Cat(i1,i2); ic.invert(r,ik);

   getQsub(ic,Qc,Qk);
   Qk.groupRecs(P,D);

   {  unsigned i=0, k=0, il=0, kl=0; 
      unsigned j, m=i1.len*QDIM, msz=m*sizeof(TQ), m2=2*m;
      TQ *qc=Qc.data;

      for (; k<Qk.dim1; ++k) { n=0; d=D[k];
         for (j=0; j<d; ++j, ++i, qc+=m2) {
            if (!memcmp(qc,qc+m,msz)) {
               if (il<i) { Qc.recSet(il,i); P[il]=P[i]; }
               ++il; ++n;
            }
         }
         if (n) {
            if (kl<k) { Qk.recSet(kl,k); }
            D[kl++]=n;
         }
      }

      if (!il) { return B.initScalar(0); }
      if (!kl) { wblog(FL,"ERR %s() kl=%d having il=%d",FCT,kl,il); }

      P.len=il; Qc.dim1=il;
      D.len=kl; Qk.dim1=kl;
   }

   B.initQ(Qk,*this);
   itags.select(mark,B.itags,0);

   if (!B.CGR.data) {
      for (i=k=0; k<Qk.dim1; ++k) { d=D[k];
      for (  j=0; j<d ; ++j, ++i) {
         DATA[P[i]]->trace(i1,i2,*B.DATA[k]); 
      }}
   }
   else {
      wbarray<TD> Cl;
      for (i=k=0; k<Qk.dim1; ++k) { d=D[k];
      for (  j=0; j<d ; ++j, ++i) { l=P[i]; Cl.init();
         DATA[l]->trace(i1,i2,Cl);

         x2.init(1,1); x2[0]=1;
         for (m=0; m<CGR.dim2; ++m) {
            x2.Kron( CGR(l,m).trace(FL,i1,i2, j? nullptr : &B.CGR(k,m)) );
         }        

         if ((n=x2.numel())>1) { unsigned rl=Cl.SIZE.len;
            if (rl!=rB+1) wblog(FL,
               "ERR %s() unexpected rank %d+%d / %d",FCT,rl-1,1,rB);
            Cl.ContractMat(FL,rB+1,x2); 
            B.DATA[k]->Plus(Cl,TD(1),j==0); 
         }
         else if (n==1) {
            B.DATA[k]->Plus(Cl,x2[0],j==0);
         }
         else wblog(FL,"ERR %s() numel(x2)=%d",FCT,n);
      }}

      if (qtype.permitsOM(r) && !qtype.permitsOM(rB)) {
         for (k=0; k<B.DATA.len; ++k) {
            B.DATA[k]->skipSingletons(rB);
         }
      }
   }

   return B;
};

template <class TQ, class TD>
int QSpace<TQ,TD>::check_CGR_cgw(const char *F, int L) const {
   int e=0;

   if (CGR.data && qtype.len!=CGR.dim2) wblog(F_L,"ERR %s() "
      "size mismatch (%dx%d/%d)",FCT,CGR.dim1,CGR.dim2,qtype.len);
   if (!CGR.dim2) return e;

   unsigned i,j; double w2;

   for (i=0; i<CGR.dim1; ++i)
   for (j=0; j<CGR.dim2; ++j) { if (CGR(i,j).cgw.len) {
      w2=CGR(i,j).cgw.norm2();
      if (!isfinite(w2) || Wb::abs(w2-1)>1e-12) { ++e;
         if (F) wblog(F,L,"ERR %s() ",FCT);
      }
   }}
   return e;
};

template <class TQ, class TD>
void QSpace<TQ,TD>::checkScalarCGS(const char *F, int L) const {

   if (CGR.data && qtype.len!=CGR.dim2) wblog(F_L,"ERR %s() "
      "size mismatch (%dx%d/%d)",FCT,CGR.dim1,CGR.dim2,qtype.len);
   if (!CGR.dim2) return;

   unsigned i=0,j=0; 
   wbvec<char> isa(qtype.len);

   for (; j<CGR.dim2; ++j) {
      if ((isa[j]=qtype[j].isAbelian())!=0) ++i;
   }

   if (i) {
      for (i=0; i<CGR.dim1; ++i)
      for (j=0; j<CGR.dim2; ++j) { if (isa[j]) CGR(i,j).checkAbelian(FL); }
   }
};

template <class TQ, class TD>
void QSpace<TQ,TD>::disp_cgs(
   const char *F, int L, const char *istr
 ) const {

   unsigned i,j; wblog(F_L,"TST CGR data (%s)",istr);
   for (i=0; i<CGR.dim1; ++i) {
   for (j=0; j<CGR.dim2; ++j) {
      printf(" %8s",CGR(i,j).sizeStr().data);
   }; printf("\n"); }

};

template <class TQ, class TD>
wbstring QSpace<TQ,TD>::sizeStrQ(const char *F, int L) const {
   wbvec<char> s(32); isConsistent(F_L);
   if (QDIM)
        s.catf(FL,"%lix(%lix%d)",QIDX.dim1,QIDX.dim2/QDIM,QDIM);
   else s.catf(FL,"%lix%li/%d",QIDX.dim1,QIDX.dim2,QDIM);
   return s.data;
};

template <class TQ, class TD>
wbstring QSpace<TQ,TD>::sizeStr(char vflag) const {

   if (isEmpty()) { return wbstring(); }

   if (allAbelian()) {
      wbvector<widx_t> D; getDim(D);
      return D.sizeStr();
   }
   else {
      wbvector<widx_t> D,DD; getDim(D,&DD);
      if (vflag) { wbvec<char> s(64);
         s.catf(FL,"%s (%s)",SSTR(D),SSTR(DD));
         return s.data;
      }
      else { return DD.sizeStr(); }
   }
};

template <class TQ, class TD>
void QSpace<TQ,TD>::print_SIZE(const char *F, int L) const {

   unsigned i=0, r=rank(F_L);

   if (F)
        { wblog(F,L," *  data SIZE listing (rank-%d QSpace)",FCT,r); }
   else { PRINTF("\n"); }

   for (; i<DATA.len; ++i) {
      PRINTF("%6d.  %s\n", i,
      DATA[i]->SIZE.toStrf("%4d","",r,"  │").data); 
   }
};

template <class TQ, class TD>
void QSpace<TQ,TD>::print_rankDATA(const char *F, int L) const {

   unsigned i=0, r, r0=rank(F_L), rmin=r0+99, rmax=r0;
   size_t m, M=1;

   wbvec<char> sout(64);
   sout.catf(0,0,"QSpace::DATA rank @ %d",r0);

   for (; i<DATA.len; ++i) { r=DATA[i]->SIZE.len;
      if (rmin>r) { rmin=r; }
      if (rmax<r) { rmax=r; }
      if (r>r0) { m=DATA[i]->numOM(r0); if (M<m) { M=m; }}
   }

   if (rmin==rmax) { if (rmax!=r0) {
        sout.catf(0,0," + %d  [M<=%ld]",rmin-r0,M);
   }}
   else if (rmax==rmin+1)
        sout.catf(0,0," + (%d, %d)  [M<=%ld]",  rmin-r0,rmax-r0,M);
   else sout.catf(0,0," + (%d .. %d)  [M<=%ld]",rmin-r0,rmax-r0,M);

   sout.check_bounds(FL,0);

   wblog(F_L," *  %s",sout.data);
};

template <class TQ, class TD>
void QSpace<TQ,TD>::EigenSymmetric( 
   QSpace<TQ,TD>      &AK, QSpace<TQ,TD>     &AD,
   QSpace<TQ,double>  &EK, QSpace<TQ,double> &ED,
   wbMatrix<double> &EM,   
   wbMatrix<unsigned> &DD, 
   int &Nkeep,             
   double Etrunc,          
   double *E0,             
   char mKD,               
   const wbperm &PA,       
   double deps,            
   double db,              
   int dmax,               
   const char *sdir        
 ) const {

   unsigned e=0; {

      if ((void*)&AK==(void*)this) { ++e; }
      if ((void*)&AD==(void*)this) { ++e; }
      if ((void*)&EK==(void*)this) { ++e; }
      if ((void*)&ED==(void*)this) { ++e; }

      if (&AK==&AD) { e+=20; }
      if (&EK==&ED) { e+=20; }

      if (e) {
         if (e>10) wblog(FL,
            "ERR got variable overlap in output (e=%d) !?",e);
         if (e>1) wblog(FL,
            "ERR multiple overlap within output space (e=%d) !?",e);
         QSpace<TQ,TD> X(*this); 

         X.EigenSymmetric(
           AK,AD,EK,ED,EM,DD, Nkeep,Etrunc,E0,mKD,PA,deps,db,dmax,sdir);
         return;
      }
   }

   unsigned i,j,d,n,l,nk,nt,r2, r=-1; bool cgflag=(gotCGS(FL)>0);
   wbMatrix<TQ> Qb;
   wbarray<TD> UX;
   wbvector<widx_t> D,Dc;
   wbvector<double> Ek,Et,Etot;

   QSpace<TQ,TD> A1;
   QSpace<TQ,double> E1;
   double Eref=0;

   iTags tk,td,ts,ts_;
   wbindex Ib;
   wbperm P;

   wbvector< QBlock<TQ,TD,double> > QB;
   wbvector< wbvector<char> > mark;
   wbvector< wbvector<double>* > EE;

   QBlock<TQ,TD,double> x;

   if (!isConsistent(r)) wbdie(FL,str);
   if (r%2) wblog(FL,"ERR %s() requires even-rank (%d)",FCT,r);
   r2=unsigned(r)/2; Ib.Index(r2);

   EK.init(); AK.init(); ED.init(); AD.init();
   if (isEmpty()) return;

   if (itags.len || cgflag) { 
      if (!itags.got_op_labels(FL) || unsigned(r)!=itags.len)
         wblog(FL,"ERR %s() got invalid operator '%s' (%d)",
         FCT,STR(itags),r
      );
   }
   else {
      if (!CGR.isEmpty() || (qtype.len && !qtype.allAbelian()))
      wblog(FL,"ERR %s() got missing itags (empty)",FCT);
   }

   if (cgflag) {
      if (r2>1) wblog(FL,
         "ERR %s() not implemented yet for rank-%d (qtype '%s')",
          FCT,r,qStr().data);

      if (!isBlockDiagMatrix(FL)) {
         MXPut(FL).add(*this,"H");
         wblog(FL,"ERR %s() H is not block-diagonal (%d)",FCT,cgflag);
      }

      for (n=CGR.numel(), i=0; i<n; ++i) {
         if (CGR[i].isScalar()) continue;
         if (!CGR[i].isIdentityCG(nullptr,1e-12)) {
            MXPut(FL,"q").add(*this,"H").add(Nkeep,"Nkeep")
             .add(Etrunc,"Etrunc").add(i+1,"i").add(n,"n");
            wblog(FL,
             "ERR expecting identity as CGR coefficients (%s; %d/%d)",
              CGR[i].cgb ? CGR[i].cgb->sizeStr().data : "!?",i+1,n
           );
         }
      }
   }

   getQsum(Ib,Qb).groupRecs(P,D);

   EE.init(D.len);
   QB.init(D.len);
   DD.init(D.len, cgflag ? 2 : 1);

   Wb::LogException ex;

   int np=1;
   if (D.len>1 && QSP_NUM_THREADS>1 && !omp_in_parallel()) {
      np=MIN( QSP_NUM_THREADS, int(D.len) ); 
   }

   D.cumsum_(Dc);

  #pragma omp parallel for num_threads(np) 
   for (unsigned i=0; i<D.len; ++i) { if (!ex) { try {
      QBlock<TQ,TD> &b=QB[i];
      wbarray<TD> MM;

      getSub(wbindex(D[i],P.data+Dc[i]),b.A,'r');
      b.toBlockMatrix(MM,r2);

      if (b.Q1!=b.Q2 || b.S1!=b.S2) {

         if (!isHConj())
              wblog(FL,"ERR %s() requires symmetric QSpace!",FCT);
         else wblog(FL,"ERR %s() Q sectors couple to others\n"
         "=> quantum symmetries not preserved !?", FCT);
      }

      wbEigenS(MM, b.U, b.S); 

      EE[i]=&(b.S); DD(i,0)=(MM.SIZE.len>1 ? MM.SIZE[1]: 0);
   }
      catch (Wb::LogException &e_) { ex+=e_; }
      catch (...) { ++ex; }
   }}
   ex.report(FLF);

   markSet(EE,mark,0,Nkeep,Etrunc,Etot,deps,db,dmax,sdir);
   nk=nt=0;

   if (E0) { (*E0)=Eref=Etot[0]; Etot-=(*E0); }

   for (i=0; i<D.len; ++i) {
       QBlock<TQ,TD> &b=QB[i]; 
       mark[i].find(b.Ik,b.It); nk+=b.Ik.len; nt+=b.It.len;
   }

   if (nk+nt!=Etot.len) wblog(FL, 
      "ERR size inconsistency (%d+%d==%d)", nk,nt,Etot.len);

   if (cgflag)
        { EM.init(Etot.len,2); }
   else { EM.init(Etot.len,1, Etot.data); }

   if (QB.len && QB[0].A.itags.len) {
      QBlock<TQ,TD> &b=QB[0]; itag_ t;

      b.A.checkQ(FL); 
      if (b.A.itags.len!=b.S1.dim2+b.S2.dim2) wblog(FL,
         "ERR %s() got invalid/empty itags (%s; s=%d+%d !?)",
         FCT, STR(b.A.itags), b.S1.dim2+b.S2.dim2
      );

      if (r2==1) { t=b.A.itags.last(); } else
      if (!b.A.itags[0].isConj()) { t.Conj(); }

      tk.init(b.S1.dim2, b.A.itags.data, 1, &t); td=tk;
      ts.init(1,&t,1,&t); ts[0].Conj(); ts_=ts;

      if (nt && mKD) {
         tk.last().SetK(); 
         td.last().SetD(); 
         for (i=0; i<2; ++i) { ts[i].SetK(); ts_[i].SetD(); }
      }
   }

   for (l=i=0; i<D.len; ++i, l+=d) {
       QBlock<TQ,TD> &b=QB[i]; 
       d=b.S.len;

       if (b.U.rank()!=2) wblog(FL,"ERR invalid rank-%d for U",b.U.rank());
       if (!b.S) wblog(FL,"ERR got empty S (%s)",SSTR(b.S));

       x.init_bare_refA(b,'U');

       if (cgflag) {
          wbvector<unsigned> sc; qtype.QDim(x.Q2.data,sc);
          DD(i,1)=DD(i,0)*sc.prod(); 

          double *em=EM.rec(l);
          for (j=0; j<d; ++j, em+=2) {
             em[0]=b.S[j] - Eref; 
             em[1]=QB[i].qdim_tot;
          }
       }

       if (b.Ik.len) {
          const wbvector<widx_t> D2_(1,x.D2.data,'r');
          b.U.select0(b.Ik,1,UX); 
          b.S.select(b.Ik,Ek); if (E0) Ek-=(*E0);

          x.updateBlockDim(b.Ik.len);
          x.initFromBlockMatrix(FL,A1, UX, tk); 

          E1.initDiagonal(x.Q2, D2_,Ek,'r'); 

          A1.Append2AndDestroy(FL,AK);
          E1.Append2AndDestroy(FL,EK);
       }

       if (b.It.len) {
          const wbvector<widx_t> D2_(1,x.D2.data,'r');
          b.U.select0(b.It,1,UX); 
          b.S.select(b.It,Et); if (E0) Et-=(*E0);

          x.updateBlockDim(b.It.len);
          x.initFromBlockMatrix(FL,A1, UX, td); 

          E1.initDiagonal(x.Q2, D2_,Et, 'r');  

          A1.Append2AndDestroy(FL,AD);
          E1.Append2AndDestroy(FL,ED);
       }

      #ifndef WB_SKIP_ASSERT
       if (b.D1.colSum(0)!=mark[i].len || b.Ik.len+b.It.len!=mark[i].len)
          wblog(FL,"ERR QB[%d/%d]: %d/%d/%d+%d", i, D.len,
          b.D1.colSum(0), mark[i].len, b.Ik.len, b.It.len
       );
      #endif
   }

   EK.qtype=qtype; 
   ED.qtype=qtype; 

   if (cgflag) {
      EK.initIdentityCGS(FL,*this);
      ED.initIdentityCGS(FL,*this);

      EM.SortRecs(P);
   }

   EK.itags=ts; ED.itags=ts_; 

   if (!PA.isEmpty() && !PA.isIdentityPerm()) {
      if (AK) { AK.Permute(PA); }
      if (AD) { AD.Permute(PA); }
   }

   AK.Sort(); 
   AD.Sort(); 

#ifndef WB_SKIP_ASSERT
   if (!AK.QIDX.isUniqueSorted()) {
      MXPut(FL,"q").add(AK.QIDX,"IDX_K").add(AD.QIDX,"IDX_T");
      wblog(FL,"ERR AK does not have unique QIDX !?");
   }
   if (!AD.QIDX.isUniqueSorted()) {
      MXPut(FL,"q").add(AK.QIDX,"IDX_K").add(AD.QIDX,"IDX_T");
      wblog(FL,"ERR AD does not have unique QIDX !?");
   }
#endif
};

template <class TQ, class TD> inline
unsigned QSpace<TQ,TD>::map2Vec(const char *F, int L,
   wbvector<TD> &V, const wbvector<unsigned> &S
){
   if (!S.isEmpty()) {
      unsigned i,s, n=DATA.len, e=0, N=0;

      if (n!=S.len) ++e; else
      for (i=0; i<n && !e; ++i) {
         s=DATA[i]->SIZE.prod(); N+=s;
         if (S[i]!=s) ++e;
      }
      if (N!=V.len) ++e;

      if (e) wblog(F,L,
      "ERR %s() size changed (%d->%d; %d)",FCT,V.len,S.sum(),__LINE__);
   }
   else {
      unsigned N=getDataSize();
      if (V.len!=N) V.init(N);
   }

   return map2Vec(V.data);
};

template <class TQ, class TD> inline
unsigned QSpace<TQ,TD>::map2Vec(TD *v, char Iflag
){
   unsigned i,s, n=DATA.len, N=0;

   for (i=0; i<n; ++i, v+=s) { s=DATA[i]->SIZE.prod();
      if (Iflag)
           memcpy(DATA[i]->data, v, s*sizeof(TD));
      else memcpy(v, DATA[i]->data, s*sizeof(TD));
      N+=s;
   }
   return N;
};

#endif

