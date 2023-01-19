/* ---------------------------------------------------------------------
 * Project : QSpace tensor library (v4.0 pre-release)
 * Class   : wbsparray (sparse array of arbitary dimension)
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

#ifndef __WB_SPARRAY_CC__
#define __WB_SPARRAY_CC__

// ------------------------------------------------------------------ //
// check whether A=this fits into B for all dimensions
// return value
//   0  A=this does not fit into B (is larger in at least one dimension)
// else A fits, where
//   1  A  has same SIZE.len (strict)
//   2  B.SIZE.len is longer than in A (B has non-zero extra dimension)
//   3  A.SIZE.len is longer than in B (A has trailing singleton dimensions)
//   4  A and B have diagonal format
// Wb,Jun22,20

template <class TD>
char wbsparray<TD>::fitsSize(const wbsparray<TD> &B, char strict) const {

   bool i1=isDiag(), i2=B.isDiag();
   if (i1 || i2) {
      if (!i1 || !i2) wblog(FL,
         "ERR %s() got mixed diag format (%d/%d)",FCT,i1,i2);
      return (D.len<=B.D.len ? 4 : 0);
   }

   unsigned i=0, n=MIN(SIZE.len,B.SIZE.len);
   const SPIDX_T *s1=SIZE.data, *s2=B.SIZE.data;

   for (; i<n; ++i) { if (s1[i]>s2[i]) return 0; }

   if (SIZE.len==B.SIZE.len) { return 1; }
   if (strict) { return 0; }

   if (SIZE.len<B.SIZE.len) {
      for (n=B.SIZE.len; i<n; ++i) { if (!s2[i]) return 0; }
      return 2;
   }
   else {
      for (n=SIZE.len; i<n; ++i) { if (s1[i]>1) return 0; }
      return 3;
   }
};

template <class TD>
char wbsparray<TD>::sameSize(const wbsparray<TD> &B, char strict) const {

   if (SIZE.len!=B.SIZE.len && !strict) {
      if (  isDiag(FL)) return (B.isSMatrix(FL,  D.len) ? 11 : 0);
      if (B.isDiag(FL)) return (  isSMatrix(FL,B.D.len) ? 12 : 0);

      char s=sameSizeUp2Singletons(B.SIZE);
      if (s) {
         if (s>2)
            return 23; 
         else {
            return (SIZE.len>B.SIZE.len ? 21 : 22);
         }
      }
      else return 0;
   }

   if (SIZE==B.SIZE) {
      if (!SIZE.len && D.len)
           return (D.len==B.D.len ? 13 : 0); 
      else return 1;
   }

   return 0;
};

template <class TD> inline
char wbsparray<TD>::sameSizeUp2Singletons(
   const wbvector<SPIDX_T> &S, const char strict) const {

   if (SIZE==S) return 1;
   if (numel()!=S.prod(0)) return 0;

   if (isDiag()) {
      if (S.len<2 || S[0]!=S[1] || S[0]!=D.len) return 0;
      if (S.len==2) return 1;
      if (!strict) {
         for (unsigned i=2; i<S.len; ++i) if (S[i]!=1) return 0;
         return 2;
      }
      return 0;
   }

   {  char ok=1;
      unsigned n,N, i=0; const SPIDX_T *s;

      if (SIZE.len<S.len)
           { n=SIZE.len; N=   S.len; s=   S.data; }
      else { n=   S.len; N=SIZE.len; s=SIZE.data; }

      for (; i<n; ++i) { if (SIZE[i]!=S[i]) { ok=0; break; }}
      if (ok) {
      for (; i<N; ++i) if (s[i]!=1) { ok=0; break; }}

      if (ok) return 2;
   }

   if (!strict) {
      unsigned i=0, j=0, ra=SIZE.len, rb=S.len;
      const SPIDX_T *sa=SIZE.data, *sb=S.data;

      for (; i<ra; ++i) { if (sa[i]!=1) {
         for (; j<rb; ++j) { if (sb[j]!=1) break; }
         if (j<rb && sa[i]==sb[j]) { ++j; }
         else return 0;
      }}
      for (; j<rb; ++j) { if (sb[j]!=1) return 0; }
      return 3;
   }

   return 0;
};

template <class TD> inline
char wbsparray<TD>::matchWithSingletons(
   const char *F, int L, const wbvector<SPIDX_T> &S, wbvector<int> &M
 ) const {

   unsigned i=0, j=0, e=0, ra=SIZE.len, rb=S.len;
   const SPIDX_T *sa=SIZE.data, *sb=S.data;

   if (!S.len) {
      if (SIZE.len) wblog(FL,"ERR %s() got %s <> ()",FCT,SSTR_(this));
      M.init(); return 0;
   }

   M.init(S.len).set(-1); 

   for (; i<ra; ++i) { if (sa[i]!=1) {
      for (; j<rb; ++j) { if (sb[j]!=1) break; }
      if (j<rb && sa[i]==sb[j]) { M[j++]=i; }
      else { ++e; break; }
   }}

   if (!e)
   for (; j<rb; ++j) { if (sb[j]!=1) { ++e; break; }}

   if (e) {
      if (F) wblog(F,L,"ERR %s() got size mismatch (%s <> %s)",
      FCT, SSTR_(this), SSTR(S));
      return e;
   }

   return 0;
};

template <class TD> inline
bool wbsparray<TD>::hasSize(SPIDX_T d1) const {
   return (SIZE.len==1 && SIZE[0]==d1);
};

template <class TD> inline
bool wbsparray<TD>::hasSize(SPIDX_T d1, SPIDX_T d2) const {
   return (SIZE.len==2 && SIZE[0]==d1 && SIZE[1]==d2);
};

template <class TD> inline
bool wbsparray<TD>::hasSize(SPIDX_T d1, SPIDX_T d2, SPIDX_T d3) const {
   return (SIZE.len==3 && SIZE[0]==d1 && SIZE[1]==d2 && SIZE[2]==d3);
};

template <class TD> inline
wbsparray<TD>& wbsparray<TD>::setRand(double x, char sflag){

   SPIDX_T n=numel(); if (!n) return *this;
   n=MIN(n,SPIDX_T(x*pow(double(n),1./rank(FL))));

   init_nnz(n); if (!n) return *this;

   SPIDX_T i,j, m=SIZE.len, *I=IDX.data;
   wbindex I1; Wb::rand R;

   for (i=0; i<n; ++i, I+=m) {
      for (j=0; j<m; ++j) {
         R.rand_(I[j], SIZE[j]);
         if (I[j]>=SIZE[j]) wblog(FL,"ERR %s() index out of bounds "
            "(%d,%d: %d; %s)",FCT,i+1,j+1,I[j],SSTR_(this)
         );
      }
   }

   IDX.makeUnique(I1); D.Select(I1);
   setRand_data(x,sflag);

   Compress(); return *this;
};

template <class TD> inline
wbsparray<TD>& wbsparray<TD>::setRand_data(double x, char sflag) {

   if (D.len) {
      SPIDX_T i=0, n=D.len; Wb::rand R;

      if (sflag) {
         if (sflag=='d' || sflag=='D') {
            for (; i<n; ++i) { R.rands(D.data[i]); }
         }
         else {
            for (; i<n; ++i) { R.randb(D.data[i]); }
         }
      }
      else {
         for (; i<n; ++i) { R.rand_(D.data[i]); }
      }

      if (x!=1) { R.scale_rand(D.data,D.len, x); }
   }

   return *this;
};

template <class TD> inline
wbsparray<TD>& wbsparray<TD>::Resize(const wbvector<SPIDX_T> &S_) {

   if (SIZE.len!=S_.len) wblog(FL,
      "ERR %s() invalid size (%d/%d)",FCT,S_.len,SIZE.len);
   if (!D.len) { return init(S_); }
   if (SIZE.len!=IDX.dim2) wblog(FL,
      "ERR %s() got size %s !?",FCT,SSTR_(this));

   SPIDX_T i,j, l=0, N=IDX.dim1, m=IDX.dim2;
   const SPIDX_T *idx=IDX.data, *s=S_.data;
   wbsparray<TD> X; TD *x=D.data;

   for (j=0; j<m; ++j) { if (s[j]<SIZE[j]) break; }
   if (j==m) {
      SIZE=S_; return *this;
   }

   for (i=0; i<N; ++i, idx+=m) {
      for (j=0; j<m; ++j) { if (idx[j]>=s[j]) break; }
      if (j==m && x[i]!=0) ++l;
   }

   if (l) {
      save2(X); init(S_,l); l=0; idx=X.IDX.data; x=X.D.data;
      for (i=0; i<N; ++i, idx+=m) {
         for (j=0; j<m; ++j) { if (idx[j]>=s[j]) break; }
         if (j==m && x[i]!=0) {
            for (j=0; j<m; ++j) { IDX(l,j)=idx[j]; }
            D[l++]=x[i];
         }
      }
   }
   else init(S_);

   return *this;
};

template <class TD> inline
SPIDX_T wbsparray<TD>::findRecSortedP(
   SPIDX_T *idx, SPIDX_T n, char lex
){
   unsigned i=0;

   if (long(n)<0) n=IDX.dim2;
   else if (n!=SIZE.len) { 
      if (lex>0) {
         if (n<SIZE.len) {
            for (i=n; i<SIZE.len; ++i) { if (SIZE[i]!=1) { i=-1; break; }}
         }
         else if (n>SIZE.len) {
            for (i=SIZE.len; i<n; ++i) { if (idx[i]!=0) { i=-1; break; }}
            n=SIZE.len;
         }
      }
      else {
         wblog(FL,"ERR %s() got n=%d/%d (%d)",FCT,n,SIZE.len,lex);
      }
   }

   if (int(i)<0 || !n || !SIZE.len) wblog(FL,
      "ERR %s() invalid index [%s] (%s)",
      FCT, wbvector<SPIDX_T>(n,idx,'r'), SSTR_(this)
   );

   return IDX.findRecSorted(idx,n,lex);
};

template<class TD> inline
char wbsparray<TD>::sameUptoFac( 
  const wbsparray<TD> &B, TD *fac_, char lex, TD eps
) const {

   widx_t i=-1, n=D.len; TD a0,b0,fac=1;
   const TD *a=D.data, *b=B.D.data;
   wbindex Ia,Ib;
   char iA=0, iB=0, sab=sameSize(B);

   if (!sab) return 1;
   if (!n) {
      if (fac_) { (*fac_)=(B.D.len ? 0. : 1.); }
      return (B.D.len ? 0 : 2);
   }

   if (sab>=20) wblog(FL,"ERR %s() got singletons (%s <> %s)",
      FCT,SSTR_(this), SSTR(B));

   if (sab>=10) {        
      iA=((sab-10) & 1); 
      iB=((sab-10) & 2);
   }

   a0=D.aMax(&i);

   if (swidx_t(i)<0) wblog(FL,
      "ERR %s() got zero data (%g;%d)",FCT,double(a0),i
   );

   if (iB) {
      if (iA) b0=B.D[i];
      else {
         const SPIDX_T *I=IDX.ref(i);
         if (IDX.dim2!=2) wblog(FL,
            "ERR %s() %dx%d !?",FCT,IDX.dim1,IDX.dim2);
         b0=(I[0]==I[1] ? B.D[I[0]] : TD(0.));
      }

   }
   else {
      b0=B.value(IDX.rec(i),IDX.dim2); 
   }

   if (a0<=eps) {
      if (fac_) { (*fac_)=((a0==0 && b0) || ABS(b0)>eps ? 0. : 1.); }
      return (ABS(b0)>eps ? 0 : 3);
   }
   if (ABS(b0)<=eps) { return 4; }

   fac=a[i]/b0; if (fac_) (*fac_)=fac;

   if (sab>=10) { 
      if (iA && iB) {
         for (i=0; i<D.len; ++i) { if (ABS(a[i]-fac*b[i])>eps) return 11; }
         return 0;
      }
      else if (iA) {
         const SPIDX_T *I=B.IDX.data; wbvector<char> mark(D.len);
         if (B.IDX.dim2!=2) wblog(FL,
            "ERR %s() %dx%d !?",FCT,B.IDX.dim1,B.IDX.dim2);
         for (i=0; i<B.D.len; ++i, I+=2) {
            if (I[0]==I[1]) { mark.el(I[0])=1;
                 if (ABS(a[I[0]]-fac*b[i])>eps) return 12; }
            else if (ABS(b[i])>eps) return 13;
         }
         for (i=0; i<D.len; ++i) {
            if (!mark[i] && ABS(a[i])>eps) return 14;
         }
         return 0;
      }
      else { 
         const SPIDX_T *I=IDX.data; wbvector<char> mark(B.D.len);
         if (IDX.dim2!=2) wblog(FL,
            "ERR %s() %dx%d !?",FCT,IDX.dim1,IDX.dim2);
         for (i=0; i<D.len; ++i, I+=2) {
            if (I[0]==I[1]) { mark.el(I[0])=1;
                 if (ABS(a[i]-fac*b[I[0]])>eps) return 21; }
            else if (ABS(a[i])>eps) return 22;
         }
         for (i=0; i<B.D.len; ++i) {
            if (!mark[i] && ABS(b[i])>eps) return 23;
         }
         return 0;
      }
   }

   matchSortedIdxU(FL,IDX,B.IDX,Ia,Ib,-1,lex);
   for (n=Ia.len, i=0; i<n; ++i) {
      if (ABS(a[Ia[i]]-fac*b[Ib[i]])>eps) return 31;
   }

   Ia.Invert(  IDX.dim1,'u');
   for (i=0; i<Ia.len; ++i) { if (ABS(a[Ia[i]])>eps) return 32; }

   Ib.Invert(B.IDX.dim1,'u');
   for (i=0; i<Ib.len; ++i) { if (ABS(b[Ib[i]])>eps) return 33; }

   return 0;
};

template<class TD> inline
bool wbsparray<TD>::isZero(TD eps, char bflag) const {

   if (!D.len) return 1;
   Wb::scale_eps(eps,D.data,D.len); 

   if (eps==0) {
      for (SPIDX_T i=0; i<D.len; ++i) if (D[i]!=0) return 0;
   }
   else {
      if (bflag) { return (D.norm()<eps); }
      else {
         for (SPIDX_T i=0; i<D.len; ++i) {
         if (ABS(D[i])>eps) return 0; }
      }
   }

   return 1;
};

template <class TD>
bool wbsparray<TD>::isDiagMatrix(TD eps) const {

   if (!isDiag()) { unsigned r=rank(FL);

      if (r!=2) return 0;
      if (IDX.dim2==1) return 1; 

      SPIDX_T i=0; const SPIDX_T *id=IDX.data;

      if (IDX.dim2!=2) sperror_this(FLF);
      Wb::scale_eps(eps,D.data,D.len); 

      if (eps==0) {
         for (; i<IDX.dim1; ++i, id+=2) { if (id[0]!=id[1]) return 0; }
      }
      else {
         for (; i<IDX.dim1; ++i, id+=2) if (id[0]!=id[1]) {
            if (ABS(D[i])>eps) return 0;
         }
      }
   }

   return 1;
};

template <class TD>
wbvector<TD>& wbsparray<TD>::getDiag(
   const char *F, int L, wbvector<TD> &dd) const {

   if (isDiag(F_L)) { dd.init(D); }
   else {
      SPIDX_T i=0, n=-1; const SPIDX_T *I=IDX.data;

      if (!isSMatrix(F_L,&n) || IDX.dim2!=2) wblog(F_L,
         "ERR %s() requires rank-2 object (%s)",FCT,SSTR_(this));
      dd.init(n);

      for (; i<IDX.dim1; ++i, I+=2) {
         if (*I==I[1]) { dd.el(*I) += D[i]; }
      }
   }
   return dd;
};

template <class TD>
TD wbsparray<TD>::getDiag(const char *F, int L, SPIDX_T k) const {

   TD x=0; 

   if (isDiag(F_L)) {
      if (sSPIDX_T(k)<0) { k+=D.len; }
      if (k>=D.len) wblog(FL,
         "ERR %s() index out of bounds (%d/%d)",FCT,k,D.len);
      x=D[k];
   }
   else {
      SPIDX_T i=0, n=-1; const SPIDX_T *I=IDX.data;

      if (!isSMatrix(F_L,&n) || IDX.dim2!=2) wblog(F_L,
         "ERR %s() requires rank-2 object (%s)",FCT,SSTR_(this));
      if (sSPIDX_T(k)<0) { k+=n; }
      if (k>=n) wblog(FL,
         "ERR %s() index out of bounds (%d/%d)",FCT,k,n);

      for (; i<IDX.dim1; ++i, I+=2) {
         if (*I==I[1] && *I==k) { x+=D[i]; }
      }
   }
   return x;
};

template <class TD>
wbsparray<TD>& wbsparray<TD>::diag2reg(
  const char *F, int L, wbsparray<TD> &A, unsigned r
) const {

   if (int(r)<0 || r==2) { 
      if (!D.len) { checkSize(F_L); A=(*this); return A; }
      if (!isDiag(FL)) sperror_this(F_LF);

      SPIDX_T i=0, n=D.len, *I; const SPIDX_T s[]={n,n};
      wbvector<SPIDX_T> S(2,s);

      A.init(S,n); A.D=D; I=A.IDX.data;
      for (; i<n; ++i, I+=2) { I[0]=I[1]=i; }
   }
   else { 
      if (!r) wblog(FL,"ERR %s() invalid rank (%d)",FCT,r);
      if (!isScalar()) sperror_this(F_LF);
      A.SIZE.init(r).set(1);
      A.IDX.init(D.len,r); A.D=D; A.isref=0; 
   }

   return A;
};

template <class TD>
wbsparray<TD>& wbsparray<TD>::diag2reg(
  const char *F, int L, unsigned r
){
   if (isref) wblog(FL,"ERR %s() got reference (%d)",FCT,isref);

   if (int(r)<0 || r==2) { 
      if (!D.len) { checkSize(F_L); return *this; }
      if (!isDiag(FL)) sperror_this(F_LF);

      SPIDX_T i=0, n=D.len, *I;
      SIZE.init(2); SIZE[0]=SIZE[1]=n;
      IDX.init(n,SIZE.len);

      for (I=IDX.data; i<n; ++i, I+=2) { I[0]=I[1]=i; }
   }
   else { 
      if (!r) wblog(FL,"ERR %s() invalid rank (%d)",FCT,r);
      if (!isScalar()) sperror_this(F_LF);
      SIZE.init(r).set(1);
      IDX.init(D.len,r); 
   }

   return *this;
};

template <class TD>
bool wbsparray<TD>::isProptoId(TD &q, TD eps) const {

   if (!D.len) return 0;
   if (isDiag(FL)) { 
      if (q==0) q=D[0];
      for (SPIDX_T i=1; i<D.len; ++i) { if (ABS(D[i]-q)>eps) return 0; }
      return 1;
   }

   unsigned j, r=rank(FL), r2=r/2;
   SPIDX_T i=0,d=1; const SPIDX_T *id=IDX.data, *s=SIZE.data;

   if (IDX.dim2!=r) sperror_this(FLF);
   Wb::scale_eps(eps,D.data,D.len); 

   if (!r || r%2) return 0; 
   for (j=0; j<r2; ++j) { d*=s[j];
      if (!s[j] || s[j]!=s[j+r2]) return 0;
   }

   if (double(eps)<=0) {
      if (D.len!=d) return 0;
      for (; i<IDX.dim1; ++i, id+=r) {
         for (j=0; j<r2; ++j) { if (id[j]!=id[j+r2]) return 0; }
         if (i || q!=0) { if (ABS(D[i]-q)>eps) return 0; }
         else q=D[i];
      }
   }
   else {
      SPIDX_T n=0;
      for (; i<IDX.dim1; ++i, id+=r) {
         for (j=0; j<r2; ++j) { if (id[j]!=id[j+r2]) {
            if (ABS(D[i])>eps) return 0;
            else break;
         }}
         if (j==r2) {
            if ((++n)!=1 || q!=0) { if (ABS(D[i]-q)>eps) return 0; }
            else q=D[i];
         }
      }
      if (!n || n!=d) return 0;
   }

   return 1;
};

template <class TD>
char wbsparray<TD>::isIdentity(TD eps) const {

   if (!D.len) return 1;
   if (isDiag(FL)) { 
      for (SPIDX_T i=1; i<D.len; ++i) { if (ABS(D[i]-1)>eps) return 2; }
      return 0;
   }

   unsigned r=rank(FL), r2=r/2;
   SPIDX_T j,i=0,d=1; const SPIDX_T *id=IDX.data, *s=SIZE.data;

   if (IDX.dim2!=r) sperror_this(FLF);
   Wb::scale_eps(eps,D.data,D.len); 

   if (!r || r%2) return 3; 
   for (j=0; j<r2; ++j) { d*=s[j];
      if (!s[j] || s[j]!=s[j+r2]) return 4;
   }

   if (eps<=0) {
      if (D.len!=d) return 5;
      for (; i<IDX.dim1; ++i, id+=r) {
         for (j=0; j<r2; ++j) { if (id[j]!=id[j+r2]) return 6; }
         if (ABS(D[i]-1)>eps) return 7;
      }
   }
   else {
      SPIDX_T n=0;
      for (; i<IDX.dim1; ++i, id+=r) {
         for (j=0; j<r2; ++j) { if (id[j]!=id[j+r2]) {
            if (ABS(D[i])>eps) return 8;
            else break;
         }}
         if (j==r2) { ++n;
            if (ABS(D[i]-1)>eps) return 9;
         }
      }
      if (!n || n!=d) return 10;
   }

   return 0;
};

template <class TD>
char wbsparray<TD>::isIdentity(
   const wbperm &P, TD *dval_, TD eps) const {

   TD dval=1;

   if (!D.len) return 1;
   if (isDiag(FL)) { 
      if (dval_) { (*dval_)=dval=D.avg(); eps*=dval; }
      for (SPIDX_T i=1; i<D.len; ++i) { if (ABS(D[i]-dval)>eps) return 2; }
      return 0;
   }

   unsigned j1,j2, r=rank(FL);
   const SPIDX_T *id=IDX.data; 

   if (IDX.dim2!=r) sperror_this(FLF);

   if (P.len) {
      if (P.len!=r || r<2) wblog(FL,"ERR %s() "
         "invalid permutation P=[%s] (r=%d)",FCT,STR(P),r);
      j1=P[0]; j2=P[1];
   }
   else {
      unsigned i=0; j1=j2=SIZE.len;
      for (   ; i<SIZE.len; ++i) { if (SIZE[i]>1) { j1=i; break; }}
      for (++i; i<SIZE.len; ++i) { if (SIZE[i]>1) { j2=i; break; }}
      for (++i; i<SIZE.len; ++i) { if (SIZE[i]>1) {       break; }}
      if (j1==SIZE.len) {
         if (SIZE.allEqual(1)) {
            if (D.len!=1) return 3;
            dval=D[0]; if (dval_) { (*dval_)=dval; }
            return 0;
         }
         else return 4;
      }
      if (j2==SIZE.len || i<SIZE.len) return 5;
   }

   if (j1>=SIZE.len || j2>=SIZE.len || j1==j2) wblog(FL,
      "ERR %s() got (%d,%d/%d)",FCT,j1,j2,SIZE.len);
   if (SIZE[j1]!=SIZE[j2]) return 6;

   Wb::scale_eps(eps,D.data,D.len); 

   if (dval_) {
      SPIDX_T i=0, l=0; dval=0;
      for (; i<IDX.dim1; ++i, id+=r) {
         if (id[j1]==id[j2]) { dval+=D[i]; ++l; }
      }
      if (l!=SIZE[j1]) return 7;
      if (l) { dval/=l; eps*=dval; }; (*dval_)=dval;
      id=IDX.data;
   }

   if (eps<=TD(0)) { if (D.len!=SIZE[j1]) return 8;
      for (SPIDX_T i=0; i<IDX.dim1; ++i, id+=r) {
         if (id[j1]!=id[j2])     return 12; 
         if (ABS(D[i]-dval)>eps) return 11; 
      }
   }
   else {
      SPIDX_T i=0, l=0;
      for (; i<IDX.dim1; ++i, id+=r) {
         if (id[j1]==id[j2]) { ++l;
                if (ABS(D[i]-dval)>eps) return 11;  }
         else { if (ABS(D[i])     >eps) return 12; }
      }
      if (!l || l!=SIZE[j1]) return 7; 
   }

   return 0;
};

template <class TD>
bool wbsparray<TD>::isSym_aux(
  const char *F, int L, const char *fct,
  const wbsparray<TD> &B0, TD eps, TD* xref,
  const char symflag,
  const char lflag
) const {

   wbsparray<TD> B; B0.transpose(F_L,B);
   SPIDX_T i,r,s; widx_t m,ma,mb;
   wbindex Ia,Ib;
   wbvector<SPIDX_T> Ja(IDX.dim1), Jb(B.IDX.dim1);
   TD x;

   checkSize(FL,"A: "); if (&B!=this) checkSize(FL,"B: ");

   if (IDX.dim2!=B.IDX.dim2 || IDX.dim2%2) {
      if (lflag) sprintf(str,
         "%s %s() only applies to even-rank objects (%ld;%ld).",
          shortFL(F,L),fct, SIZE.len, B.SIZE.len);
      return 0;
   }
   if (isEmpty() && B.isEmpty()) return 1;

   m=matchIndex(IDX,B.IDX,Ia,Ib,1,&ma,&mb);
   if (m) wblog(FL,"ERR %s() input sparse matrices "
      "not compressed (%d,%d,%d) !?",FCT,m,ma,mb);

   if (xref) (*xref)=0;

   if (symflag=='s') {
      for (i=0; i<Ia.len; ++i) { r=Ia[i]; s=Ib[i];
         if ((++Ja[r])>1 || (++Jb[s])>1) { 
            wblog(FL,"ERR %s()",FCT);
         }
         x=ABS(D[r]-CONJ(B.D[s])); if (xref) *xref=MAX(*xref,x);
         if (x>eps) return 0;
      }
   }
   else if (symflag=='a') {
      for (i=0; i<Ia.len; ++i) { r=Ia[i]; s=Ib[i];
         if ((++Ja[r])>1 || (++Jb[s])>1) { 
            wblog(FL,"ERR %s()",FCT);
         }
         x=ABS(D[r]+CONJ(B.D[s])); if (xref) *xref=MAX(*xref,x);
         if (x>eps) return 0;
      }
   }
   else wblog(FL,"ERR invalid symflag=%c<%d>",symflag,symflag);

   for (i=0; i<Ja.len; ++i) { if (!Ja[i] && ABS(  D[i])>eps) return 0; }

   if (&B0!=this)
   for (i=0; i<Jb.len; ++i) { if (!Jb[i] && ABS(B.D[i])>eps) return 0; }

   return 1;
};

template <class TD>
TD wbsparray<TD>::dotProd(
  const char *F, int L, const wbsparray<TD> &B, char tnorm
) const {

   TD x=TD(); 
   char cflag = ((tnorm || WbUtil<TD>::hasNoConj()) ? 0 : 1);
   char iA=isDiag(), iB=B.isDiag();

   if (iA || iB) { 
      if ((!iA && SIZE.len!=2) || (!iB && B.SIZE.len!=2) || dim()!=B.dim())
         wblog(FL,"ERR %s() size mismatch (%s <> %s)",
         FCT,SSTR_(this),SSTR(B)
      );

      if (iA && iB) {    
         const TD *a=D.data, *b=B.D.data; SPIDX_T i=0;
         if (cflag)
              { for (; i<D.len; ++i) x+=CONJ(a[i])*b[i]; }
         else { for (; i<D.len; ++i) x+=     a[i] *b[i]; }
      }
      else if (iA)
           { x=B.dotProd_full_diag(*this, cflag ? 2:0); }
      else { x=  dotProd_full_diag(B,     cflag ? 1:0); } 

      return x;
   }

   if (SIZE!=B.SIZE) wblog(F_L,"ERR %s() "
      "got size mismatch (%s <> %s)",FCT,SSTR_(this),SSTR(B));
   if (!SIZE.len) { return x; }

   SPIDX_T ia=0, ib=0, na=IDX.dim1, nb=B.IDX.dim1, l=0, lnext=0; 
   const SPIDX_T *Ia=IDX.data, *Ib=B.IDX.data, Davg=SQRT(SIZE.prod());
   const TD *a=D.data, *b=B.D.data;
   unsigned m=SIZE.len;
   char q=0, lex=0; 

   while (ia<na && ib<nb) {
      q=Wb::recCompare(Ia,Ib,m,lex);

      if (q<0) { ++ia; Ia+=m; } else
      if (q>0) { ++ib; Ib+=m; } else {
         if (cflag)
              { x+=CONJ(a[ia])*b[ib]; }
         else { x+=     a[ia] *b[ib]; }
         ++ia; ++ib; Ia+=m; Ib+=m;

         if ((++l)>=lnext) { lnext+=Davg;
            if (ia<na) { q=Wb::recCompare(Ia-m,Ia,m,lex);
               if (q>=0) wblog(F_L,"ERR %s() A: "
               "input sparse IDX not sorted (%ld/%ld) !?",FCT,ia,na); }
            if (ib<nb) { q=Wb::recCompare(Ib-m,Ib,m,lex);
               if (q>=0) wblog(F_L,"ERR %s() B: "
               "input sparse IDX not sorted (%ld/%ld) !?",FCT,ib,nb);
            }
         }
      }
   }

   return x;
};

template <class TD> inline
TD wbsparray<TD>::dotProd_full_diag(const wbsparray<TD> &B, char cflag
) const {

   const TD *a=D.data, *b=B.D.data;
   const SPIDX_T *I=IDX.data;
   TD x=TD();

   if (IDX.dim2!=2 || B.IDX.dim1 || B.D.len!=SIZE[0]) wblog(FL,
      "ERR %s() %s <> %s !?",FCT, SSTR_(this), SSTR(B)
   );

   if (cflag==0) { 
      for (SPIDX_T i=0; i<D.len; ++i, I+=2) {
         if (I[0]==I[1]) { x += a[i] * b[I[0]]; }
      }
   }
   else if (cflag==1) {
      for (SPIDX_T i=0; i<D.len; ++i, I+=2) {
         if (I[0]==I[1]) { x += CONJ(a[i]) * b[I[0]]; }
      }
   }
   else if (cflag==2) {
      for (SPIDX_T i=0; i<D.len; ++i, I+=2) {
         if (I[0]==I[1]) { x += a[i] * CONJ(b[I[0]]); }
      }
   }
   else wblog(FL,"ERR %s() invalid cflag=%d",FCT,cflag);

   return x;
};

template <class TD>
wbsparray<TD>& wbsparray<TD>::TimesEl(
   const wbsparray<TD> &B, char tnorm
){
   char cflag = ((tnorm || WbUtil<TD>::hasNoConj()) ? 0 : 1);
   char iA=0, iB=0, sab=sameSize(B);

   if (!sab || sab>=20) wblog(FL,"ERR %s() "
      "size mismatch (%s <> %s)",FCT,SSTR_(this),SSTR(B));

   if (sab>=10) { wbsparray<TD> X;

      iA=((sab-10) & 1); 
      iB=((sab-10) & 2);

      if (iA && iB) {
         SPIDX_T i=0; const TD *a=D.data; TD *x;
         X=B; x=X.D.data;
         if (cflag)
              { for (; i<D.len; ++i) x[i]*=CONJ(a[i]); }
         else { for (; i<D.len; ++i) x[i]*=     a[i] ; }
      }
      else if (iA)
           { B.timesEl_full_diag(*this,X, cflag ? 2:0); }
      else {   timesEl_full_diag(B,    X, cflag ? 1:0); } 
   }
   else {
      wbindex Ia,Ib; SPIDX_T i=0, m,ma=0,mb=0;
      const TD *a=D.data, *b=B.D.data;

      m=matchIndex(IDX,B.IDX,Ia,Ib,1,&ma,&mb);
      if (m) wblog(FL,"ERR %s() input sparse "
        "objects not compressed (%d,%d,%d) !?",FCT,m,ma,mb);

      wbsparray<TD> X(SIZE,Ia.len);
      TD *x=X.D.data;

      for (i=0; i<Ia.len; ++i) { X.IDX.recSetP(i, IDX.rec(Ia[i])); }
      if (cflag)
           { for (i=0; i<Ia.len; ++i) { x[i]=CONJ(a[Ia[i]])*b[Ib[i]]; }}
      else { for (i=0; i<Ia.len; ++i) { x[i]=     a[Ia[i]] *b[Ib[i]]; }}

      X.save2(*this).Compress();
   }

   return *this;
};

template <class TD> inline
wbsparray<TD>& wbsparray<TD>::timesEl_full_diag(
  const wbsparray<TD> &B, wbsparray<TD> &X, char cflag
) const {

   const TD *a=D.data, *b=B.D.data;
   const SPIDX_T *I=IDX.data;

   X.initDiag(B.D.len); TD *x=X.D.data;

   if (IDX.dim2!=2 || B.IDX.dim1 || B.D.len!=IDX.dim1) wblog(FL,
      "ERR %s() (%d,%d)x%d <> (%d,%d)x%d !?",FCT,
      D.len, IDX.dim1,IDX.dim2, B.D.len, B.IDX.dim1,B.IDX.dim2
   );

   if (cflag==0) {
      for (SPIDX_T i=0; i<IDX.dim1; ++i, I+=2) { if (I[0]==I[1]) {
         x[I[0]]=a[i]*b[I[0]];
      }}
   }
   else if (cflag==1) {
      for (SPIDX_T i=0; i<IDX.dim1; ++i, I+=2) { if (I[0]==I[1]) {
         x[I[0]]=CONJ(a[i])*b[I[0]];
      }}
   }
   else if (cflag==2) {
      for (SPIDX_T i=0; i<IDX.dim1; ++i, I+=2) { if (I[0]==I[1]) {
         x[I[0]]=a[i]*CONJ(b[I[0]]);
      }}
   }
   else wblog(FL,"ERR %s() invalid cflag=%d",FCT,cflag);
   return X;
};

template <class TD> inline
wbsparray<TD>& wbsparray<TD>::setCol(
   const char *F, int L, SPIDX_T k, const wbsparray<TD> &v
){
   unsigned rv=v.rank(F_L);

   if (rank(F_L)!=2 || rv<1 || rv>2 || SIZE[0]!=v.SIZE[0] ||
      (rv==2 && v.SIZE[1]!=1)
    ) wblog(F_L,"ERR %s() invalid input (%s <> %s)",
      FCT, SSTR_(this), SSTR(v)
   );

   if (k>=SIZE[1]) wblog(F_L,
      "ERR %s() index out of bounds (%d/%d)",FCT,k,SIZE[1]);
   if ((IDX.dim1 && IDX.dim2!=2) || v.IDX.dim2!=rv) wblog(F_L,
      "ERR %s() %dx%d <> %dx%d/%d !?",FCT,IDX.dim1,IDX.dim2,
      v.IDX.dim2,v.IDX.dim2,rv);

   SPIDX_T i=0, l=0, l2, n=IDX.dim1, *I2=IDX.data+1;

   for (; i<n; ++i, I2+=2) { if ((*I2)!=k) {
       if (l<i) { setRec(l,i); }; ++l;
   }}

   if (l<i && F) wblog(F_L,
      "WRN %s() *this already contains data at col=%d/%d (%d/%d)",
      FCT,k,SIZE[1],l,SIZE[0]
   );
   if (!l) { 
      if (v.IDX.dim2==1)
           { v.IDX.resize(v.IDX.dim1, 2, IDX); }
      else { IDX=v.IDX; }
      IDX.setLastCol(k);

      D=v.D; return *this;
   }

   n=v.D.len; l2=l+n;
   if (l2) {
      if (l2<=IDX.dim1)
           { REDSIZE_NNZ(l2); }
      else { IDX.Resize(l2,IDX.dim2); D.Resize(l2); }
   }
   else { IDX.init(0,IDX.dim2); D.init(); }

   if (n) {
      const SPIDX_T *Iv=v.IDX.data;
      Wb::cpyRange(D.data+l, v.D.data, n); I2=IDX.rec(l);
      for (i=0; i<n; ++i, I2+=2, Iv+=rv) {
         I2[0]=(*Iv); 
         I2[1]=k;
      }
   }
   return *this;
};

template <class TD> inline
TD wbsparray<TD>::Normalize(char tnorm, char qflag){
   TD x = SQRT(Wb::overlap(D.data,D.data,D.len,1,tnorm));

   if (ABS(x)>TD(1E-15)) Wb::timesRange(D.data,TD(1)/x,D.len);
   else if (!qflag) {
      wblog(FL,"ERR %s() got vector with norm %.3g !?",FCT,double(x));
   }

   return x;
};

template <class TD>
TD wbsparray<TD>::NormalizeCol(SPIDX_T k, char tnorm, char qflag){

   if (rank()!=2 || IDX.dim2!=2) wblog(FL,
      "ERR %s() requires rank-2 (got %s)",FCT,SSTR_(this));
   if (k>=SIZE[1]) wblog(FL,
      "ERR %s() index out of bounds (%d/%d)",FCT,k,SIZE[1]);
   checkSize(FLF);

   SPIDX_T i=0, n=IDX.dim1, *I2=IDX.data+1;
   opFlags<TD> xflags(tnorm ? 'N' : 'C');
   TD x=0;

   if (xflags.conj())
        { for (; i<n; ++i, I2+=2) if ((*I2)==k) { x+=CONJ(D[i])*D[i]; }}
   else { for (; i<n; ++i, I2+=2) if ((*I2)==k) { x+=D[i]*D[i]; }}

   x=SQRT(x);
   if (ABS(x)>1E-15) { x=1/x;
      for (i=0; i<n; ++i, I2+=2) if ((*I2)==k) { D[i]*=x; }
   }
   else if (!qflag) {
      wblog(FL,"ERR %s() got vector with norm %.3g !?",FCT,x);
   }

   return x;
};

template <class TD>
wbsparray<TD>& wbsparray<TD>::QRdecomp(const char *F, int L,
   wbsparray<TD> &Q, 
   char posR, 
   TD eps     
){
   if (rank()!=2 || SIZE.len!=2 || IDX.dim2!=2) wblog(F_L,
      "ERR %s() requires rank-2 array (%s)",FCT,SSTR_(this));

   SPIDX_T d1=SIZE[0], d2=SIZE[1], *idx=IDX.data;
   if (!d1 || !d2) { Q.init(); return *this; }

   SPIDX_T i,l, k=0, k_, l1=0, l2=0, iu=0, N=IDX.dim1;
   TD x2, c2=0, e2=0, eps2=eps*eps;
   wbsparray<TD> u,X;
   int dk=0;

   eps2*=SQRT(TD(MAX(d1,d2)));

   if (eps2>TD(1E-8)) wblog(F_L,
      "ERR %s() got eps=%.3g !?",FCT,double(eps));

   if (!N) { SPIDX_T one=1;
      Q.initIdentity(d1,one); init(one,d2); 
      return *this;
   }

   if (N==1) { 
      wbvector<SPIDX_T> s(2); wbsparray<TD> R;
      s[0]=d1; s[1]=1; Q.init(s,1); Q.IDX(0,0)=IDX(0,0);
      s[0]=1; s[1]=d2; R.init(s,1); R.IDX(0,1)=IDX(0,1);

      if (posR && D[0]<0)
           { Q.D[0]=-1; R.D[0]=-D[0]; }
      else { Q.D[0]= 1; R.D[0]= D[0]; }

      return R.save2(*this);
   }

   wbvector< wbsparray<TD> > U(d2);

   wbvector<char> sR; 
   if (posR) { posR=1; sR.init2val(d2,1); }

   wbsparray<TD> A_(*this); 

   for (; k<=d2; ++k) {
      N=IDX.dim1; 

      if (k) {
         if (l1>=N) l1=N-1;
         if (IDX(l1,1)>=k_) { 
            while (l1>0) { if (IDX(--l1,1)<k_) { ++l1; break; }}
         }
         if (IDX(l1,1)!=k_) wblog(FL,
            "ERR %s() QR out of pace (%ld/%ld) !?",FCT,IDX(l1,1),k_);
         l2=l1;
      }

      idx=IDX.rec(l2); k_=idx[1]; dk=-1;

      while (l2<N) {
         while (l2<N && idx[0]<k && idx[1]==k_) { ++l2; idx+=2; }

         if ((++dk)==0 && posR && k) {
            const TD &rk=D[l2-1];

            if (idx[-2]!=k-1) {
               MXPut(FL,"q").add(A_,"A")
                 .add(*this,"R").add(IDX,"IDX").add(D,"D")
                 .add(l1,"l1").add(l2,"l2").add(k,"k").add(k_,"k_")
                 .add(c2,"c2").add(e2,"e2").add(eps,"eps");
               wblog(FL,
                 "ERR %s() got R(%ld/%ld,%ld)=%.4g (l1=%ld, l2=%ld)",
                  FCT,idx[-2],k-1,k_,double(rk),l1,l2-1
               );
            }
            if (NORM2(rk)<=eps2) {
               sprintf(str,"got R(%ld,%ld)=%.4g",k-1,k_,double(rk));
               if (rk) wblog(FL,"WRN %s() %s",FCT,str);
               else wblog(FL,"ERR %s() %s",FCT,str);
            }
            if (rk<0) { sR[k-1]=-1; posR=-1; }
         }

         if (l2>=N && idx[-2]<k) { l1=l2; break; }
         if (idx[1]!=k_) {
            k_=idx[1]; continue; 
         }

         l1=l2; c2=NORM2(D[l2]);      idx+=2; ++l2;
         while (l2<N && idx[1]==k_) { idx+=2; c2+=NORM2(D[l2++]); }

         if (c2<=eps2 || (l1 && dk==0)) { e2+=c2;
            if (c2>100*eps2) wblog(FL,"ERR %s() "
               "got c2=%.3g (%.3g) !?",FCT,double(c2),double(eps2));
            for (l=l1; l<l2; ++l) D[l]=0;
            if (l2<N) {
               k_=idx[1]; continue; 
            }
            else { l1=l2; break; }
         }
         else {
            break;
         }
      }; if (l1>=N) break;

      if (k_<k || k>=d2) {
         MXPut(FL,"q").add(A_,"A")
           .add(*this,"R").add(IDX,"IDX").add(D,"D")
           .add(l1,"l1").add(l2,"l2").add(k,"k").add(k_,"k_")
           .add(c2,"c2").add(e2,"e2").add(eps,"eps");
         wblog(FL,"ERR %s() k_=%ld/%ld/%ld !?",FCT,k_,k,d2);
      }

      x2=getHouseholderVec(FL,k,k_,u,eps2); 

      if (x2<eps2*100) { 
         if (x2>eps2) wblog(F_L,"WRN %s() got x2 @ %.3g / %.3g "
            "[%ldx%ld; %d]",FCT,double(x2),double(eps2),d1,d2,dk);
         continue;
      }

      LProject1(u,X,-2); 
      Plus(FL,X); 
      u.save2(U[iu++]);
   }

   if (!k || k>d2) wblog(FL,"ERR %s() got k=%ld/%ld !?",FCT,k,d2);

   if (k<d1) {
      N=IDX.dim1; idx=IDX.data; c2=0;
      for (l=0; l<N; ++l, idx+=2) {
         if (idx[0]>=k) { c2+=NORM2(D[l]); }
      }
      if (c2>eps2) { MXPut(FL,"q").add(*this,"R").add(c2,"c2");
         wblog(FL,"ERR %s() got non-triangular R matrix (%.3g) !?",
         FCT,double(c2));
      }
      SPIDX_T s[2]={k,d2};
      Resize(wbvector<SPIDX_T>(2,s));
   }
   e2+=SkipTiny(eps);

   Q.initIdentity(d1,MIN(d1,k)); 
   if (iu) {
      for (k=iu-1; k<iu; --k) {
         Q.LProject1(U[k],X,-2);  
         Q.Plus(FL,X); 
      }
   }
   e2+=Q.SkipTiny(eps);

   if (posR<0) {
      TD *x=Q.D.data; N=Q.IDX.dim1; idx=Q.IDX.data;
      for (i=0; i<N; ++i, idx+=2) { if (sR[idx[1]]<0) x[i]=-x[i]; }
      x=D.data; N=IDX.dim1; idx=IDX.data;
      for (i=0; i<N; ++i, idx+=2) { if (sR[idx[0]]<0) x[i]=-x[i]; }
   }

   return *this;
};

template <class TD>
wbsparray<TD>& wbsparray<TD>::OrthoNormalizeColsQR(
  const char *F, int L,
  TD eps 
){
   wbsparray<TD> Q;

   QRdecomp(F_L,Q,'+',eps);

   if (SIZE.len==2 || D.len) {
      SPIDX_T l, n=D.len, d2=SIZE[1]; char e=0; int lmax=-1; 
      wbvector<TD> dd(d2); TD eps=1E-12;
      wbvector<int> id; id.init2val(d2,-1);

      for (l=0; l<n; ++l) { if (NORM(D[l])>eps) {
         dd[IDX(l,1)]=D[l]; 
         id[IDX(l,1)]=IDX(l,0);
      }}

      for (l=0; l<d2 && !e; ++l) {
         if (id[l]>lmax) { 
            if (id[l]==lmax+1) { lmax=id[l];
               if (dd[l]<0) { e=1;
                  sprintf(str,"dd[%ld]=%.3g",l,double(dd[l]));
               }
            }
            else { e=2; sprintf(str,"%ld: %d -> %d",l,lmax,id[l]); }
         }
      }

      if (e) { MXPut(FL,"q") 
         .add(*this,"R").add(Q,"Q").add(SIZE,"S")
         .add(IDX,"IDX").add(D,"D").add(dd,"dd").add(id,"id");
         wblog(FL,"ERR %s() %s (eps=%g) !?",FCT,str,double(eps));
      }
   }

   return Q.save2(*this);
};

template <class TD>
wbsparray<TD>& wbsparray<TD>::OrthoNormalizeColsGS( 
   const char *F, int L,
   char tnorm,  
   char qflag,  
   TD eps,      
   unsigned np  
){
   if (rank()!=2) wblog(FL,
      "ERR %s() requires rank-2 array (%s)",FCT,SSTR_(this));

   SPIDX_T i,j,ip;
   char xflag=(qflag=='x' || qflag=='X');
   wbvector< wbsparray<TD> > X;
   wbvector<SPIDX_T> I0;
   TD a,z;

   splitSparseCM(F_L,X, xflag ? eps : TD(0));

   WbUtil<TD>().adjust_tnorm(tnorm);

   for (i=0; i<X.len; ++i) {
      for (ip=0; ip<np; ++ip) {
         for (j=0; j<i; ++j) { if (X[j].D.len) {
            z=X[i].dotProd(FL,X[j],tnorm); if (z!=0) {
            X[i].Plus(FL,X[j],-z); }
         }}
      }

      z=X[i].norm2(tnorm); a=ABS(z);

      if (a>eps) { X[i]*=(TD(1)/SQRT(a)); }
      else {
         if (qflag) {
            X[i].initz();
         }
         else if (!xflag) { 
            MXPut(FL,"x").add(*this,"V").add(X,"X")
              .add(qflag,"qflag").add(xflag,"xflag")
              .add(i+1,"i").add(a,"a").add(eps,"eps");
            wblog(F_L,"ERR |U(:,%d/%d)| = %.3g (%g) %s !?", i+1, X.len,
               double(a), double(eps), tnorm ? " (using t-norm)" : ""
            );
         }
      }
   }

   wbvector< wbsparray<TD>* > xp(X.len);
   if (xflag) {
      for (j=i=0; i<X.len; ++i) { if (X[i].D.len) {
         xp[j++]=(&X[i]);
      }}
      if (j<i) {
         if (!j) wblog(F_L,"ERR %s() got all null-vectors",FCT);
         xp.len=j;
      }
   }
   else {
      for (i=0; i<X.len; ++i) { xp[i]=(&X[i]); }
   }

   initCAT(FL,xp,&I0);
   if (!IDX.dim1) wblog(F_L,"ERR %s() got all null-vectors",FCT);

   return *this;
};

template <class TD>
template <class T2>
wbvector<T2>& wbsparray<TD>::norm2vec(unsigned l, 
   wbvector<T2> &a, char tnorm
 ) const {

   if (l>=SIZE.len) wblog(FL,
      "ERR %s() index out of bounds (l=%d/%d)",FCT,l,SIZE.len);
   a.init(SIZE[l]);

   unsigned m=IDX.dim2;
   SPIDX_T i=0, *ip=IDX.data+l;
   T2 *ad=a.data;

   if (tnorm)
        { for (; i<D.len; ++i, ip+=m) { ad[*ip] += D[i]*D[i]; }}
   else { for (; i<D.len; ++i, ip+=m) { ad[*ip] += CONJ(D[i])*D[i]; }}

   return a;
};

template <class TD>
TD wbsparray<TD>::trace() const {

   TD x=0; 

   unsigned r=rank(FL); if (!r) return x;

   if (r==2 && !IDX.dim1) return D.sum(); 
   if (r%2 || r!=SIZE.len || r!=IDX.dim2) sperror_this(FLF);

   SPIDX_T i=0; unsigned r2=r/2, n=r2*sizeof(SPIDX_T);
   const SPIDX_T *idx=IDX.data;

   for (; i<r2; ++i) {
      if (SIZE[i]!=SIZE[i+r2]) wblog(FL,
      "ERR %s() got non-symmetric tensor (%s)",FCT,SSTR_(this));
   }

   for (i=0; i<IDX.dim1; ++i, idx+=r) {
      if (memcmp(idx,idx+r2,n)==0) x+=D[i];
   }

   return x;
};

template <class TD>
template <class T2>
wbvector<T2>& wbsparray<TD>::trace(
   unsigned k, wbvector<T2> &A 
 ) const {

   SPIDX_T i=0;
   unsigned r=IDX.dim2, r2=(r-1)/2, n=r2*sizeof(SPIDX_T);
   const SPIDX_T *idx=IDX.data+(k?0:1), *sz=SIZE.data+(k?0:1);
   int l=(k ? r-1 : -1);

   if (k && k+1!=SIZE.len) wblog(FL,
      "ERR %s() only accepts k=1 or k=rank (%d/%d)",FCT,k+1,r);
   if (SIZE.len%2!=1) wblog(FL,"ERR %s() "
      "requires odd-rank tensor (r=%d; k=%d)",FCT,SIZE.len,k);
   if (SIZE.len!=IDX.dim2 || IDX.dim1!=D.len) sperror_this(FLF);

   A.init(SIZE[k]);
   T2 *a=A.data;

   for (; i<r2; ++i) {
      if (sz[i]!=sz[i+r2]) wblog(FL,"ERR %s() "
      "got non-symmetric tensor (%s; k=%d)",FCT,SSTR_(this),k+1);
   }

   for (i=0; i<IDX.dim1; ++i, idx+=r) {
      if (memcmp(idx,idx+r2,n)==0) {
         if (idx[l]>=A.len) wblog(FL,"ERR %s() "
            "index out of bounds (%d,%d: %d/%d)",FCT,i,l,idx[l],A.len);
         a[idx[l]]+=D[i];
      }
   }

   return A;
};

bool mxIsWbsparray(
   const char *F, int L, const mxArray *a, unsigned k
){

   if (!a || mxIsEmpty(a)) { return 0; }

   SPIDX_T n=mxGetNumberOfElements(a);
   if (k>=n) { if (F) wblog(F,L,
      "ERR sparse::%s() index exceeds dimension (%d/%d)",FCT,k,n);
      return 0;
   }

   int e[3] {
      mxGetFieldNumber(a,"S"),
      mxGetFieldNumber(a,"idx"), mxGetFieldNumber(a,"data")
   };

   if (e[0]<0 || e[1]<0 || e[2]<0) { if (F) wblog(FL,
      "ERR %s() invalid sparse array (missing fields%s%s%s)",FCT,
       e[0]?" S":"", e[1]?" idx":"", e[2]?" data":"");
      return 0;
   }

   return 1;
};

bool mxIsWbsparray(
   const char *F, int L, const mxArray *a, unsigned k,
   const mxArray **as_, const mxArray **ad_, const mxArray **ai_
){

   if (!a || mxIsEmpty(a)) { return 0; }

   SPIDX_T n=mxGetNumberOfElements(a);
   if (k>=n) { if (F) wblog(F,L,
      "ERR sparse::%s() index exceeds dimension (%d/%d)",FCT,k,n);
      return 0;
   }

   const mxArray
      *as=mxGetField(a,k,"S"),
      *ai=mxGetField(a,k,"idx"), 
      *ad=mxGetField(a,k,"data");

   if (!as || !ad) { if (F) {
      if (as && mxGetNumberOfElements(as)) wblog(F_L,
         "ERR %s() invalid sparse array (0x%lX, 0x%lX: S (%d)",
          FCT,as,ad, as ? mxGetNumberOfElements(as):-1);
      if (ad && mxGetNumberOfElements(ad)) wblog(F_L,
         "ERR %s() invalid sparse array (0x%lX, 0x%lX: data (%d)",
          FCT,as,ad, ad ? mxGetNumberOfElements(ad):-1);
      }
      return 0;
   }

   if (as_) (*as_)=as;
   if (ad_) (*ad_)=ad;
   if (ai_) (*ai_)=ai;

   return 1;
};

template <class TD>
wbsparray<TD>& wbsparray<TD>::toScalar(const char *F, int L, unsigned r) {

   if (D.len>1 || (SIZE.len && !SIZE.allEqual(1))) wblog(F_L,
      "ERR %s() invalid scalar %s",FCT,info2Str().data);

   if (!D.len) { D.init(1); D[0]=0; }

   if (int(r)<=0) {
      if (SIZE.data) SIZE.init();
      if (IDX.dim1 || IDX.dim2) IDX.init();
   }
   else {
      if (SIZE.len) { if (SIZE.len!=r) wblog(FL,
         "ERR %s() got rank mismatch (%d/%d)",FCT,SIZE.len,r); }
      else {
         SIZE.init(r).set(1);
         IDX.init(1,r);
      }
   }

   return *this;
};

template <class TD>
TD wbsparray<TD>::getScalar(const char *F, int L) const {

   if (D.len>1 || (SIZE.len && !SIZE.allEqual(1))) wblog(F_L,
      "ERR %s() invalid scalar %s",FCT,info2Str().data);

   return (D.len ? D[0] : TD(0));
};

template <class TD>
wbvector<TD>& wbsparray<TD>::getCol(SPIDX_T k, wbvector<TD> &v) const {

   if (isDiag(FL)) { 
      if (sSPIDX_T(k)<0) { k+=D.len; }
      if (k>=D.len) wblog(FL,
         "ERR %s() index out of bounds (%d/ %s)",FCT,k,SSTR_(this));
      v.init(D.len); v[k]=D[k];
   }
   else {
      if (SIZE.len!=2) wblog(FL,"ERR %s() invalid S=%d",FCT,SSTR_(this));
      if (sSPIDX_T(k)<0) { k+=SIZE[1]; } 
      if (k>=SIZE[1]) wblog(FL,
         "ERR %s() index out of bounds (%d/ %s)",FCT,k,SSTR_(this));

      SPIDX_T i=0, n=SIZE[0]; const SPIDX_T *I=IDX.data;
      v.init(n);

      for (; i<IDX.dim1; ++i, I+=2) { if (I[1]==k) {
         if (*I>=n) wblog(FL,"ERR %s() "
            "index out of bounds (%d/%d/ %s)",FCT,*I,n,SSTR_(this));
         v[*I]+=D[i];
      }}
   }
   return v;
};

template <class TD>
wbvector<TD>& wbsparray<TD>::getRow(SPIDX_T k, wbvector<TD> &v) const {

   if (isDiag(FL)) { 
      if (sSPIDX_T(k)<0) { k+=D.len; }
      if (k>=D.len) wblog(FL,
         "ERR %s() index out of bounds (%d/ %s)",FCT,k,SSTR_(this));
      v.init(D.len); v[k]=D[k];
   }
   else {
      if (SIZE.len!=2) wblog(FL,"ERR %s() invalid S=%d",FCT,SSTR_(this));
      if (sSPIDX_T(k)<0) { k+=SIZE[0]; } 
      if (k>=SIZE[0]) wblog(FL,
         "ERR %s() index out of bounds (%d/ %s)",FCT,k,SSTR_(this));

      SPIDX_T i=0, n=SIZE[1]; const SPIDX_T *I=IDX.data;
      v.init(n);

      for (; i<IDX.dim1; ++i, I+=2) { if (*I==k) {
         if (I[1]>=n) wblog(FL,"ERR %s() "
            "index out of bounds (%d/%d/ %s)",FCT,I[1],n,SSTR_(this));
         v[I[1]]+=D[i];
      }}
   }
   return v;
};

template <class TD>
int wbsparray<TD>::checkTrailingSingletons(
   const char *F, int L, unsigned *r
 ) const {

   if (SIZE.len!=IDX.dim2) wblog(F_L,
      "ERR %s() size mismatch (%d/%d)",FCT,SIZE.len,IDX.dim2);
   if (!SIZE.len) {
      if (r && int(*r)>=0) wblog(FL,
         "ERR %s() out of bounds (%d/%d)",FCT,*r,SIZE.len);
      return 0;
   }

   unsigned l=SIZE.len-1;
   for (; l<SIZE.len; --l) { if (SIZE[l]!=1) break; }; ++l;

   if (r) {
      if (int(*r)>=0) {
         if ((*r)>SIZE.len) wblog(FL,
            "ERR %s() out of bounds (%d/%d)",FCT,*r,SIZE.len);
         if (l>(*r)) wblog(FL,"ERR %s() "
            "got non-singletons for dim>=%d (%d,%d)",FCT,*r,l,SIZE.len);
         l=(*r);
      }
      else (*r)=l;
   }

   return (SIZE.len-l); 
};

template <class TD>
int wbsparray<TD>::skipTrailingSingletons(
   const char *F, int L, wbsparray<TD> &C, unsigned r
 ) const {

   C.init();

   if (SIZE.len && checkTrailingSingletons(F_L,&r)) {
      unsigned m=SIZE.len-r;

     #ifndef WB_SKIP_ASSERT
      SPIDX_T i=0, n=IDX.dim1, zero=0;
      for (; i<n; ++i) if (Wb::anyUnequal(IDX.ref(i,r),m,zero)) {
         wblog(F_L,"ERR %s() IDX out of bounds (%d/%dx%d: %d/%d) !?",
         FCT,i+1,IDX.dim1,IDX.dim2,r,SIZE.len);
      }
     #endif

      C.SIZE.init(r,SIZE.data);
      IDX.resize(IDX.dim1,r, C.IDX);
      C.D=D;

      return m; 
   }

   if (int(r)>=0) wblog(FL,
      "ERR %s() out of bounds (%d/%d)",FCT,r,SIZE.len);

   return 0;
};

template <class TD>
int wbsparray<TD>::SkipTrailingSingletons(
   const char *F, int L, unsigned r
){

   if (SIZE.len && checkTrailingSingletons(F_L,&r)) {
      unsigned m=SIZE.len-r;

     #ifndef WB_SKIP_ASSERT
      SPIDX_T i=0, n=IDX.dim1, zero=0;
      for (; i<n; ++i) if (Wb::anyUnequal(IDX.ref(i,r),m,zero)) {
         wblog(F_L,"ERR %s() IDX out of bounds (%d/%dx%d: %d/%d) !?",
         FCT,i+1,IDX.dim1,IDX.dim2,r,SIZE.len);
      }
     #endif

      if (r) SIZE.len=r; else SIZE.init();
      IDX.Resize(IDX.dim1,r);

      return m; 
   }

   if (int(r)>=0) wblog(FL,
      "ERR %s() out of bounds (%d/%d)",FCT,r,SIZE.len);

   return 0;
};

template <class TD> inline
bool wbsparray<TD>::hasSingletons(const wbindex &I) const {

   if (!I.len) wblog(FL,"ERR %s() got empty index set",FCT);
   for (unsigned i=0; i<I.len; ++i) {
      if (I[i]>=SIZE.len || SIZE[I[i]]!=1) return 0;
   }

   return 1;
};

template <class TD>
wbsparray<TD>& wbsparray<TD>::skipSingletons(
   const char *F, int L, const wbindex &I, wbsparray<TD> &A
 ) const {

   if (!I.len) { A=*this; return A; }
   if (&A==this) {
      wbsparray<TD> X; A.save2(X);
      return X.skipSingletons(F,L,I,A);
   }

   if (!hasSingletons(I)) wblog(FL,"ERR %s() got S=(%s) "
      "given I=[%s] !?", FCT,SSTR_(this), STR(I));

   wbindex I2; I.invert(IDX.dim2,I2,'u');

   SIZE.select(I2,A.SIZE); A.isref=0;
   IDX.getCols(I2,A.IDX); A.D=D;

   return A;
};

template <class TD>
wbsparray<TD>& wbsparray<TD>::AddTrailingSingletons(
   const char *F, int L, unsigned r
){
   unsigned rk=rank(FL); if (r==rk) { return *this; }

   if (int(r)<0) wblog(F_L,"ERR %s() invalid rank (%d/%d)",FCT,r,rk);
   if (r<rk) wblog(F_L,"ERR %s() cannot reduce rank (%d/%d)",FCT,r,rk);

   if (isDiag()) {
      if (IDX.data) wblog(F_L,
         "ERR %s() got IDX data (%s) !?",FCT,IDX.dim1,IDX.dim2);
      SIZE.init(r); SIZE[0]=SIZE[1]=D.len;
         for (unsigned i=2; i<r; ++i) SIZE[i]=1;
      IDX.init(D.len,r);
         for (SPIDX_T *I=IDX.data, i=0; i<D.len; ++i, I+=r) {
         I[0]=I[1]=i; 
      }
   }
   else {
      if (SIZE.len) {
         unsigned j=SIZE.len; SIZE.Resize(r);
         for (; j<r; ++j) SIZE[j]=1;
      }

      if (IDX.dim1) { 
         if (!SIZE.len) wblog(F_L,
            "ERR %s() got empty size (%d/%d)",FCT,SIZE.len,r);
         IDX.Resize(IDX.dim1,r); 
      }
   }

   return *this;
};

template <class TD>
char wbsparray<TD>::checkSize2(
  const char *F, int L, const wbsparray<TD> &B
) const {

    char sab=sameSize(B);
    this->checkSize(F_L,"A:"); B.checkSize(F_L,"B:");

    if (sab==1) {
       if (IDX.dim2!=B.IDX.dim2){ wblog(F_L,
             "ERR sparse::%s() IDX size mismatch (%dx%d <> %dx%d)",
             FCT,IDX.dim1,IDX.dim2,B.IDX.dim1,B.IDX.dim2);
          return -1;
       }
    }
    else if (sab<=0) {
       wblog(F_L,"ERR %s() mismatch %s <> %s (%d)",
          FCT,SSTR_(this),SSTR(B),sab);
    }

    return sab;
};

template <class TD> inline
char wbsparray<TD>::checkSize(
  const char *F, int L, const char *fct, const char *istr
) const {

   if (SIZE.len) {
      if (IDX.data) { 
         if (IDX.dim1!=D.len || IDX.dim2!=SIZE.len) {
            if (F) sperror_this(F_L_F,istr); else return 2;
         }
      }
      else if (D.len) {
         if (SIZE.len==2 && SIZE[0]==D.len && SIZE[1]==D.len) {
            wblog(FL,"WRN %s() got SIZE for diag-format (%s; d=%d)",
            FCT,SSTR_(this),D.len);
         }
         else { if (F) sperror_this(F_L_F,istr); else return 3; }
      }
   }
   else {
      if (!IDX.isEmpty()) {
         if (F) sperror_this(F_L_F,istr); else return 1;
      }
   }
   return 0;
};

template <class TD> inline
void wbsparray<TD>::check_IDX_range(const char *F, int L) const {

   if (!SIZE.len) {
      if (IDX.dim1 || IDX.dim2) sperror_this(F_LF);
      return;
   }

   if (IDX.dim1!=D.len ||
      (IDX.dim2!=SIZE.len && (IDX.dim1 || IDX.dim2))
   ) sperror_this(F_LF);

   if (IDX.dim1) {
      SPIDX_T i=0; unsigned j=0, m=IDX.dim2;
      const SPIDX_T *s=SIZE.data, *I=IDX.data;

      for (; i<IDX.dim1; ++i, I+=m) {
         for (j=0; j<m; ++j) {
            if (I[j]>=s[j]) wblog(FL,
               "ERR %s() index out of bounds (%d: %d <> %s)",
               FCT,j+1,I[j]+1, SSTR_(this)
            );
         }
      }
   }
};

template <class TD>
wbsparray<TD>& wbsparray<TD>::make2D(const wbindex &ic, wbsparray<TD> &B,
   char pos 
 ) const {

   if (B.D.len) { B.init(); }; B.D=D;
   IDX.toIndex2D(SIZE,ic,B.IDX,pos); 

   SPIDX_T s2, s1=SIZE.prod(ic.data,ic.len,s2);

   if (pos==1)
        { SPIDX_T s[2]={s1,s2}; B.SIZE.init(2,s); }
   else { SPIDX_T s[2]={s2,s1}; B.SIZE.init(2,s); }

   B.Compress(); return B;
};

template <class TD>
sparseIndex2D<SPIDX_T>& wbsparray<TD>::toIndex2D(
   const wbindex &ic, sparseIndex2D<SPIDX_T> &a,
   char pos, 
   char use_rmaj
 ) const {

   unsigned i,j;
   if (pos==1) { i=0; j=1; } else { i=1; j=0; }

   a.S.init(2);
   a.S[i]=SIZE.prod(ic.data,ic.len,a.S[j]);

#ifndef WB_SKIP_ASSERT
 { unsigned k; double s1=1, s2=1; wbvector<char> mark(SIZE.len);
   for (k=0; k<ic.len; ++k) { ++mark[ic[k]]; }
   for (k=0; k<SIZE.len; ++k) { (mark[k] ? s1 : s2)*=SIZE[k]; }
   if (double(a.S[i])!=s1 || double(a.S[j]!=s2)) wblog(FL,
      "ERR %s() type range out of bounds !?\n(%gx%g <> %ldx%ld)",
      FCT,s1,s2,a.S[0],a.S[1]
   );
 }
#endif

   IDX.toIndex2D(SIZE,ic,a.IJ,pos); 

   return a.Setup(use_rmaj);
};

template <class TD>
mxArray* wbsparray<TD>::toMxSp() const {

   if (!Wb::isBaseType(typeid(TD))) {
      double x;
      for (SPIDX_T i=0; i<D.len; ++i) { x=double(D[i]);
         if (fabs(double(TD(x)-D[i])/x)>1E-20) {
             MXPut(FL,"ans").add(*this,"data").add(i+1,"i").add(x,"x");
             wblog(FL,"ERR %s() got data type `%s' (d[%d]=%g @ %.3g)",
                FCT,TSTR(TD), i+1,x,double(TD(x)-D[i])
             );
         }
      }
   }

   if (isScalar()) {
      if (D.len) {
          if (D.len>1) sperror_this(FLF);
          return numtoMx(D[0]);
      }
      return numtoMx(0);
   }

   if (isDiag(FL)) {
      wbMatrix<SPIDX_T> IJ; SPIDX_T d1=1,d2=1;
      get2DIndex(IJ,d1,d2);

      return Wb::mxCreateSparse(FL,d1,d2,IJ, D);
   }

   if (!SIZE.len) {
      if (IDX.isEmpty() && D.isEmpty())
         return Wb::mxCreateSparse(FL,0,0,IDX,D);
      else wblog(FL,
     "ERR %s() invalid wbsparray (%s;%d)",FCT,SSTR_(this),D.len);
   }

   if (SIZE.len==2) {
      return Wb::mxCreateSparse(FL,SIZE[0],SIZE[1],IDX,D);
   }

   if (SIZE.len==1) {
      wbMatrix<SPIDX_T> IX(IDX.dim1,2);
      SPIDX_T i=0, *I=IX.data;

      if (IDX.dim2!=1) sperror_this(FLF);
      for (; i<IDX.dim1; ++i, I+=2) {
         I[0]=IDX.data[i]; I[1]=0;
      }

      return Wb::mxCreateSparse(FL,SIZE[0],1,IX,D);
   }

   wblog(FL,"ERR %s() requires rank-2 wbsparray (%s)",
   FCT,SSTR_(this)); return 0;
};

template <class TD>
mxArray* wbsparray<TD>::mxCreateStructX(
   unsigned m, unsigned n, int ma
) const {
   const char *fields[]={"S","idx","data","info","A"};
   return mxCreateStructMatrix(m,n,ma!=-99 ? 5:4,fields);
};

template <class TD>
void wbsparray<TD>::add2MxStructX(mxArray *a, unsigned i, int ma) const {

   char isId=isIdentityMatrix(), isc=isScalar();

   mxSetFieldByNumber(a,i,0, SIZE.toMx());

   if (isc) {
      mxSetFieldByNumber(a,i,1, IDX.toMx());    
      mxSetFieldByNumber(a,i,2, D.toMx('t'));   
      mxSetFieldByNumber(a,i,3, wbstring("scalar").toMx());
   }
   else if (!IDX.isEmpty() && IDX.allEqual(0) && D.allEqual(0)) {
      mxSetFieldByNumber(a,i,1,mxCreateSparse(IDX.dim1,IDX.dim2,0,mxREAL));
      mxSetFieldByNumber(a,i,2,mxCreateSparse(1,D.len,0,mxREAL));
      mxSetFieldByNumber(a,i,3, wbstring("(init)").toMx());
   }
   else if (isId) {
      mxSetFieldByNumber(a,i,2, numtoMx(D[0]));
      mxSetFieldByNumber(a,i,3, wbstring("identity").toMx()); 
   }
   else if (isDiag(FL)) { 
      mxSetFieldByNumber(a,i,2, D.toMx('t'));
      mxSetFieldByNumber(a,i,3, wbstring("diag").toMx()); 
   }
   else if (IDX.dim2!=1 && isDiagMatrix()) { 
      mxSetFieldByNumber(a,i,1, IDX.getCols(0,0).toMx());
      mxSetFieldByNumber(a,i,2, D.toMx('t'));
   }
   else {
      mxSetFieldByNumber(a,i,1, IDX.toMx());
      mxSetFieldByNumber(a,i,2, D.toMx('t'));
   };

   if (ma!=-99) {
      wbMatrix<SPIDX_T> IJ; SPIDX_T d1=1,d2=1;

      if (!SIZE.len) { d1=d2=0; }
      else ma=get2DIndex(IJ,d1,d2,ma);

      mxSetFieldByNumber(a,i,4,
         Wb::mxCreateSparse(FL,d1,d2,IJ,D)
      );
   }
};

template <class TD>
mxArray* wbsparray<TD>::mxCreateSTRUCT(
   unsigned m, unsigned n
) const {
   const char *fields[]={"S","idx","data","A"};
   return mxCreateStructMatrix(m,n,4,fields);
};

template <class TD>
mxArray* wbsparray<TD>::add2MxSTRUCT(mxArray *a, unsigned i) const {

   wbMatrix<SPIDX_T> IJ; SPIDX_T d1=1,d2=1; char tflag='t';
   get2DIndex(IJ,d1,d2);

   mxSetFieldByNumber(a,i,0, SIZE.toMx());
   mxSetFieldByNumber(a,i,1, IDX.toMx());

#ifdef QS_USING_MPFR
   if (typeid(TD)==typeid(Wb::quad)) {
      tflag=30; 
   }
#endif

   mxSetFieldByNumber(a,i,2, D.toMx(tflag));

   mxSetFieldByNumber(a,i,3,Wb::mxCreateSparse(FL,d1,d2,IJ,D));

   return a;
};

template <class TD> inline
mxArray* wbsparray<TD>::mxCreateStruct(unsigned m, unsigned n) const {
   const char *fields[]={"S","data"};
   return mxCreateStructMatrix(m,n,2,fields);
};

#ifdef QS_USING_MPFR

template <> inline mxArray*
wbsparray<Wb::quad>::mxCreateStruct(unsigned m, unsigned n) const {
   const char *fields[]={"S","idx","data"};
   return mxCreateStructMatrix(m,n,3,fields);
};

#endif

template <class TD>
void wbsparray<TD>::add2MxStruct(mxArray *a, unsigned i, int ma) const {

   wbMatrix<SPIDX_T> IJ;
   SPIDX_T d1,d2; get2DIndex(IJ,d1,d2,ma);
   int e=0;

   mxSetFieldByNumber(a,i,0, SIZE.toMx());
   mxSetFieldByNumber(a,i,1, Wb::mxCreateSparse(FL,d1,d2,IJ,D,NULL,&e));

   if (e && !isDiag()) {
      static unsigned ic=0; char s[8];
      snprintf(s,8,"i%02d",ic++);
      MXPut(FL,s).add(*this,"A").addP(toMX(),"A0").add(e,"e");
      wblog(FL,"ERR %s() got e=%d",FCT,e);
   }
};

#ifdef QS_USING_MPFR

template <>
void wbsparray<Wb::quad>::add2MxStruct(
   mxArray *a, unsigned i, int ma  __attribute__ ((unused))
 ) const {

   mxSetFieldByNumber(a,i,0, SIZE.toMx ());
   mxSetFieldByNumber(a,i,1, IDX .toMxT()); 
   mxSetFieldByNumber(a,i,2, D   .toMx ());
};

#endif

template <class TD> inline
mxArray* wbsparray<TD>::mxCreateCell(unsigned m, unsigned n) const {
   return mxCreateCellMatrix(m,n);
};

template <class TD>
void wbsparray<TD>::add2MxCell(mxArray *a, unsigned i) const {

   if (rank()>2) wblog(FL,
      "ERR %s() cell array expects rank <=2 (%s)",FCT,SSTR_(this));
   mxSetCell(a,i,toMxSp());

};

template <class TD>
mxArray* wbsparray<TD>::IDtoMx() const {

   const char *fields[]={"S","nnz","norm","stat"};
   mxArray *S=mxCreateStructMatrix(1,1,4,fields);

   const char *field2[]={"istat","idx","d"};
   mxArray *as=mxCreateStructMatrix(1,1,3,field2);

   if (isDiag()) {
      SPIDX_T s[2]={D.len,D.len};
      mxSetFieldByNumber(S,0,0, wbvector<SPIDX_T>(2,s).toMx());
   } else
   mxSetFieldByNumber(S,0,0, SIZE.toMx());

   mxSetFieldByNumber(S,0,1, numtoMx(nnz()));
   mxSetFieldByNumber(S,0,2, numtoMx(double(norm())));

   if (D.len>3) {
      unsigned i=0; const unsigned m=3;
      const TD *d=D.data; TD a;

      wbindex I(m); widx_t *k=I.data;
      wbvector<double> x_(m); double *x=x_.data;
         x[0]=d[0];
         x[1]=ABS(d[0]);
         x[2]=d[0];

      for (++i; i<D.len; ++i) { a=ABS(d[i]);
         if (x[0]>d[i]) { x[0]=d[i]; k[0]=i; } 
         if (x[1]>a   ) { x[1]=a;    k[1]=i; } 
         if (x[2]<d[i]) { x[2]=d[i]; k[2]=i; } 
      }

      mxSetFieldByNumber(as,0,0, I  .toMx());
      mxSetFieldByNumber(as,0,1, IDX.getRecs(I).toMx());
      mxSetFieldByNumber(as,0,2, x_ .toMx());
   }
   else if (D.len) {
      mxSetFieldByNumber(as,0,1, IDX.toMx());
      mxSetFieldByNumber(as,0,2, D  .toMx());
   }

   mxSetFieldByNumber(S,0,3,as);

   return S;
};

template <class TD>
wbsparray<TD>& wbsparray<TD>::initCAT(
   const char *F, int L, wbvector< wbsparray<TD>* > &X,
   wbvector<SPIDX_T> *I0, char del
){
   SPIDX_T i,id,l, r=0, i0=-1, i2=0, n=0, nz=0;
   wbvector<SPIDX_T> S;

   for (i=0; i<X.len; ++i) {
      if (!X[i] || X[i]->isEmpty()) continue;
      if (X[i]==this) {
         wbsparray<TD> A; A.initCAT(F,L,X,I0,del);
         return A.save2(*this);
      }
      X[i]->checkSize(F_LF);
      nz+=X[i]->nnz();

      if ((++n)==1) { i0=i; r=X[i]->SIZE.len+1;
         S.init(r,X[i]->SIZE); 
      }; i2=i;

      if (i!=i0 && !X[i0]->sameSize(*X[i])) wblog(F_L,
         "ERR %s() severe size mismatch (%d: %s <> %d: %s)",FCT,
         i0+1, SSTR_(X[i0]), i+1, SSTR_(X[i])
      );
   }
   if (long(i0)<0) { init(); return *this; } 

   char gotI0=(I0 && I0->len);
   if (gotI0 && I0->len!=X.len+1) wblog(F_L,
      "ERR %s() size mismatch (%d/%d+1)",FCT,I0->len,X.len);

   S.last()=(gotI0 ? I0->last() : n);
   init(S,nz);

   for (id=-1, l=0, i=i0; i<=i2; ++i) {
      if (X[i] && !X[i]->isEmpty()) {
         if (gotI0) { id=I0->data[i];
            if (id>=S.last()) { wblog(F_L,
               "ERR %s() index out of bounds (%d/%d)",FCT,id,S.last());
            }
         }
         else ++id;

         if (X[i]->nnz()) l+=catRecs(l,*X[i], id);
      }
   }

   if (l!=D.len) wblog(F_L,"ERR %s() %d/%d !?",FCT,l,D.len);

   if (del) for (i=0; i<X.len; ++i) {
      if (X[i]) { WB_DELETE_1(X[i]); }
      else break; 
   }

   Compress(); return *this;
};

template <class TD>
SPIDX_T wbsparray<TD>::catRecs(
   SPIDX_T i0, const wbsparray<TD> &a, SPIDX_T id
){
   SPIDX_T n=a.IDX.dim1;

   if (IDX.dim2!=a.IDX.dim2+1) wblog(FL,
      "ERR %s() size mismatch (%d+1/%d)",FCT,IDX.dim2,a.IDX.dim2);
   if (i0+n>IDX.dim1) wblog(FL,
      "ERR %s() index out of bounds (%d+%d = %d ?)",FCT,i0,n,IDX.dim1);

   if (!n && a.D.len==1) {
      if (IDX.dim2!=1 || IDX.dim1!=D.len) wblog(FL,"ERR %s() "
         "%s <> %dx%d (%d)",FCT,SSTR_(this),IDX.dim1,IDX.dim2,i0);
      if (a.D[0])
           { IDX(i0,0)=id; D[i0]=a.D[0]; return 1; }
      else { return 0; }
   }

   if (a.D.len!=n) wblog(FL,"ERR %s() %s <> %dx%d (%d)",
      FCT, SSTR(a), a.IDX.dim1, a.IDX.dim2, i0);

   Wb::cpyRange(D.data+i0, a.D.data, n);
   Wb::cpyStride(
      IDX.rec(i0), a.IDX.data, a.IDX.dim2, a.IDX.dim1,
      IDX.dim2
   );

   for (SPIDX_T l=IDX.dim2-1, i=0; i<n; ++i) {
      IDX(i0+i,l)=id;
   }

   return n;
};

template <class TD>
wbsparray<TD>& wbsparray<TD>::Cat(
   const char *F, int L,
   const wbsparray<TD> &B, unsigned k  
){
   if (this==&B) wblog(FL,"ERR %s() got same object !?",FCT);
   if (  isEmpty()) { (*this)=B; return *this; }
   if (B.isEmpty()) { return *this; }

   if (!SIZE.len || !B.SIZE.len) wblog(F_L,
      "ERR %s() got empty size (%s <> %s @ %d) !?",
      FCT, SSTR_(this), SSTR(B), k);
   checkSize(FLF); B.checkSize(FLF);

   if (!k || (k>SIZE.len && k>99)) wblog(F_L,
      "ERR %s() k out of bounds (%d/%d)",FCT,k,SIZE.len);
   if (k>SIZE.len) AddTrailingSingletons(FL,k);
   --k; 

   SPIDX_T i=0, n=MIN(SIZE.len, B.SIZE.len); int e=0;
   SPIDX_T D1=SIZE[k], d1=D.len, d2=B.D.len, d12=d1+d2;

   for (; i<n; ++i) { if (i!=k) {
      if (SIZE[i]!=B.SIZE[i]) { e=1; break; }
   }}
   if (e==0) { if (i==k) ++i;
      for (; i<  SIZE.len; ++i) { if (  SIZE[i]!=1) { e=2; break; }}
      for (; i<B.SIZE.len; ++i) { if (B.SIZE[i]!=1) { e=3; break; }}
   }
   if (e) wblog(F_L,
      "ERR %s() size mismatch (%s <> %s @ %d; e=%d)",
      FCT, SSTR_(this), SSTR(B), k, e
   );

   SIZE[k]+=(B.SIZE.len>k ? B.SIZE[k] : 1);
   D.Resize(d12, B.D.data);

   if (!IDX.dim1) {
      if (B.IDX.dim2==SIZE.len)
           { IDX=B.IDX; }
      else { B.IDX.resize(B.IDX.dim1,SIZE.len, IDX); }
   }
   else {
      if (B.IDX.dim2==IDX.dim2) {
         IDX.Resize(d12,IDX.dim2, B.IDX.data);
      }
      else {
         SPIDX_T j, *a, n=MIN(IDX.dim2,B.IDX.dim2);
         const SPIDX_T *b=B.IDX.data;

         IDX.Resize(d12,IDX.dim2); a=IDX.ref(d1,0);
         for (i=0; i<d2; ++i, a+=IDX.dim2, b+=B.IDX.dim2) {
         for (j=0; j<n; ++j) { a[j]=b[j]; } }
      }
      for (SPIDX_T i=d1; i<d12; ++i) IDX(i,k)+=D1;
   }

   Compress(); return *this;
};

template <class TD>
wbsparray<TD>& wbsparray<TD>::initCAT(
   const char *F, int L, const wbvector< const wbsparray<TD>* > &X,
   unsigned k 
){
   if (!X.len) { init(); return *this; }

   SPIDX_T i,j, i0=-1, i2=0, n=0;
   wbvector<SPIDX_T> S0, DD(X.len), dd(X.len);

   for (i=0; i<X.len; ++i) { if (X[i] && !X[i]->isEmpty()) {
      const wbvector<SPIDX_T> &S2=X[i]->SIZE;
      if (X[i]==this) {
         wbsparray<TD> A; A.initCAT(F,L,X,k);
         return A.save2(*this);
      }

      if (S2.len!=X[i]->IDX.dim2) X[i]->sperror_this(F_L);
      X[i]->checkSize(F_LF);

      if ((++n)==1) { i0=i; S0=S2;
         if (!k || k>S2.len) wblog(FL,"ERR %s() k out of bounds "
            "(%s; %d)",FCT,SSTR(S2),k);
         --k; 
      }
      else {
         if (S0.len!=S2.len) wblog(FL,"ERR %s() rank mismatch "
            "(%d: %s <> %d/%d: %s; %d)",FCT,i0+1,SSTR(S0),
            i+1,X.len, SSTR(S2), k+1);
         for (j=0; j<S0.len; ++j) {
            if (j!=k && S0[j]!=S2[j]) wblog(FL,"ERR %s() size mismatch "
            "(%d: %s <> %d/%d: %s; %d/%d)",FCT,i0+1,SSTR(S0),
            i+1,X.len, SSTR(S2), j+1,k+1);
         }
      }

      DD[i]=S2[k];
      dd[i]=X[i]->D.len; i2=i;
   }}

   SPIDX_T dall=dd.sum(), Dall=DD.sum();

   if (!n) { init(); return *this; }
   S0[k]=Dall; init(S0,dall); if (!dall) return *this;

   SPIDX_T *pi=IDX.data,N,s, d1=0, d2, D1=0, m=S0.len;
   TD *pd=D.data;

   for (i=i0; i<=i2; ++i, d1+=n, D1+=N) { n=dd[i]; N=DD[i];
      s=n*m; MEM_CPY<SPIDX_T>(pi,s,X[i]->IDX.data); pi+=s;
      s=n  ; MEM_CPY<TD>(pd,s,X[i]->D.data); pd+=s;
      for (d2=d1+n, j=d1; j<d2; ++j) IDX(j,k)+=D1;
   }
   if (d1!=IDX.dim1) wblog(FL,"ERR %s() %d/%d",FCT,d1,IDX.dim1);

   Compress(); return *this;
};

template <class TD>
wbsparray<TD>& wbsparray<TD>::init( 
   const char *F, int L, const mxArray *a, unsigned k
){

   init(); if (!a) {
      if (k) wblog(FL,"ERR %s() got k=%d for null mxArray !?",FCT,k);
      return *this;
   }

   if (mxIsEmpty(a) && mxIsDouble(a)) {
      SPIDX_T i,r=mxGetNumberOfDimensions(a);
      if (r) {
         const size_t *sp=mxGetDimensions(a);
         wbvector<SPIDX_T> S(r);

         for (i=0; i<r; ++i) S[i]=(SPIDX_T)sp[i];
         init(S,0);
      }
      return *this;
   }

   if (mxIsStruct(a)) {
      const mxArray *as,*ad,*ai;
      if (!mxIsWbsparray(F_L,a,k,&as,&ad,&ai))wblog(FL,
         "ERR %s() invalid sparse array",FCT);

      SIZE.init(F_L,as);
      IDX .init(F_L,ai);
      D   .init(F_L,ad);

      check_IDX_range(F_L);

      if (IDX.isEmpty()) {
         if (!SIZE.isEmpty() && !IDX.dim2) {
            if (IDX.dim1) wblog(FL,
               "ERR %s() got size %s !?",FCT,SSTR_(this));
            IDX.init(IDX.dim1,SIZE.len);
         }
      }
      if (!isDiag()) {
         if (IDX.dim1!=D.len || (IDX.dim1 && IDX.dim2!=SIZE.len))
         sperror_this(F_L);
      }

      Compress(); 
      return *this;
   }

   if (mxIsCell(a)) {
      a=mxGetCell(a,k); if (!a) return *this;
   }
   else if (k) wblog(FL,"ERR %s() got k=%d !?",FCT,k);

   if (mxIsSparse(a)) { 

      SPIDX_T j,k,l,m,n,d, nnz=mxGetNzmax(a);
      mwIndex *Ir,*Jc;

      if (!mxIsDouble(a) || mxGetNumberOfDimensions(a)!=2) wblog(FL,
         "ERR %s() invalid input (got %s)",FCT,mxTypeSize2Str(a).data);

      m=mxGetM(a); n=mxGetN(a); initz(m,n,nnz);
      Ir=mxGetIr(a);
      Jc=mxGetJc(a);

      if (mxIsComplex(a)) {
         if (typeid(TD)!=typeid(wbcomplex)) wblog(FL,
            "ERR %s() got complex sparse array for %s",FCT,TSTR(TD));
         wbcomplex *dz=(wbcomplex*)mxGetComplexDoubles(a);
         for (l=j=0; j<n; ++j    ) { d=SPIDX_T(Jc[j+1]-Jc[j]);
         for (  k=0; k<d; ++k,++l) { D[l]=dz[l];
            IDX(l,0)=SPIDX_T(Ir[l]); IDX(l,1)=j;
         }}
      }

      else {
         double *dr=mxGetDoubles(a);
         for (l=j=0; j<n; ++j    ) { d=SPIDX_T(Jc[j+1]-Jc[j]);
         for (  k=0; k<d; ++k,++l) { D[l]=dr[l];
            IDX(l,0)=SPIDX_T(Ir[l]); IDX(l,1)=j;
         }}
      }

      if (l && l!=nnz) wblog(FL, 
      "ERR severe sparse inconsistency (mex: %d,%d)",l,nnz);
   }
   else { 

      if (!mxIsDouble(a)) wblog(FL,
         "ERR %s() got `%s' data",FCT,mxGetClassName(a));
      if (typeid(TD)==typeid(wbcomplex)) wblog(FL,
         "WRN %s() ignores complex input data",FCT);

      SPIDX_T r=mxGetNumberOfDimensions(a);
      const size_t *sp=mxGetDimensions(a);
      const double *dr=mxGetDoubles(a);

      wbvector<SPIDX_T> S(r), I(r);
      SPIDX_T i,k,N, nz=0, iz=0, l=r-1, *s=S.data, *ip=I.data;

      for (i=0; i<r; ++i) S[i]=(SPIDX_T)sp[i];

      N=S.prod(0);
      for (i=0; i<N; ++i) { if (dr[i]!=0) ++nz; }
      init(S,nz);

      for (i=0; i<N; ++i) {
          if (dr[i]!=0) {
             if (iz>=nz) wblog(FL,
                "ERR %s() index out of bounds (%d/%d) !?",FCT,iz+1,nz);
             D[iz]=dr[i]; IDX.recSetP(iz,ip); ++iz;
          }

          k=0; ++ip[0]; 
          while (ip[k]>=s[k] && k<l) { ip[k]=0; ++ip[++k]; }
      }
   }

   return *this;
};

template <class TD>
wbsparray<TD>& wbsparray<TD>::init(
   const char *F, int L, const wbarray<TD> &A, TD eps
){
   unsigned k, r=A.rank(), l=r-1;
   wbvector<SPIDX_T> I(r);
   SPIDX_T i,N, nz=0, iz=0, *ip=I.data; const size_t *s=A.SIZE.data;

   N=A.numel();
   for (i=0; i<N; ++i) { if (ABS(A.data[i])>eps) ++nz; }
   init(A.SIZE,nz);

   for (i=0; i<N; ++i) {
       if (ABS(A.data[i])>eps) {
          if (iz>=nz) wblog(F_L,
             "ERR %s() index out of bounds (%d/%d) !?",FCT,iz+1,nz);
          D[iz]=A.data[i]; IDX.recSetP(iz,ip); ++iz;
       }

       k=0; ++ip[0]; 
       while (ip[k]>=s[k] && k<l) { ip[k]=0; ++ip[++k]; }
   }

   return *this;
};

template <class TD>
template <class TA>
wbarray<TA>& wbsparray<TD>::toFull(wbarray<TA> &A) const {

   if (isDiag()) { A.initDiag(D); }
   else {
      unsigned k, r=rank(FL), l=r-1;
      SPIDX_T i,j; const SPIDX_T *s=SIZE.data, *ip=IDX.data;

      A.init(SIZE);
      if (IDX.dim1 && A.data!=NULL) {
         for (i=0; i<IDX.dim1; ++i, ip+=IDX.dim2) {
             for (j=ip[l], k=l-1; k<r; --k) { j=j*s[k]+ip[k]; }
             A.data[j]=TA(D.data[i]);
         }
      }
   }
   return A;
};

template <class TD>
template <class T2>
wbvector<T2>& wbsparray<TD>::toFull(wbvector<T2> &A) const {

   if (SIZE.len) {
      unsigned m=IDX.dim2, r=rank(FL);
      SPIDX_T i=0; const SPIDX_T *s=SIZE.data, *ip=0;

      for (; i<r; ++i) { if (s[i]!=1) {
         if (!ip) { ip=IDX.data+i; }
         else wblog(FL,
            "ERR %s() invalid usage (got %s tensor)",
            FCT,SSTR_(this)
         );
      }}
      if (!ip) ip=IDX.data; 

      A.init(numel());
      for (i=0; i<D.len; ++i, ip+=m) {
         if ((*ip)>=A.len) wblog(FL, 
            "ERR %s() invalid sparse setting (%d/%d)",FCT,*ip,A.len);
         if (A[*ip]) wblog(FL,"ERR %s() " 
            "invalid sparse setting (%d: %g)",FCT,*ip,double(A[*ip]));
         A[*ip]=D.data[i];
      }
   }
   else { A.init(1);
      if (D.len) {
         if (D.len>1) wblog(FL,
            "ERR %s() invalid scalar (len=%d)",FCT,D.len);
         A.data[0]=D.data[0];
      }
      else { A.data[0]=0; }
   }
   return A;
};

template <class TD> inline
wbsparray<TD>& wbsparray<TD>::reshape(
   const wbvector<SPIDX_T> &S, wbsparray<TD> &A
 ) const {

   if (S==SIZE) { A=*this; return A; }

   if (!S.len && (SIZE.allEqual(1) && D.len<=1)) {
      A.init(); A.D=D;
      return A;
   }

   if (isDiag()) { 
      wbsparray<TD> X; diag2reg(FL,X);
      return X.reshape(S,A);
   }

   if (sameSizeUp2Singletons(S)) {
      unsigned j; SPIDX_T i=0, *I; const SPIDX_T *I0=IDX.data; int *ms;

      wbvector<int> Ms; matchWithSingletons(FL,S,Ms);
      A.init(S); A.D=D; A.IDX.init(IDX.dim1,S.len);
      ms=Ms.data; I=A.IDX.data;

      for (; i<IDX.dim1; ++i, I+=S.len, I0+=IDX.dim2) {
         for (j=0; j<S.len; ++j) {
            if (ms[j]>=0) { I[j]=I0[ms[j]]; }
         }
      }
      return A;
   }

   if (numel()!=S.prod(0)) wblog(FL, 
      "ERR %s() got size mismatch (%s => %s)",
      FCT,SSTR_(this), SSTR(S)
   );

   A.init(S); A.D=D; A.IDX.init(IDX.dim1,S.len);

   if (!IDX.isEmpty()) {
      SPIDX_T i,k, *I=A.IDX.data; unsigned r0=SIZE.len, r=S.len;
      const SPIDX_T *I0=IDX.data, *S0=SIZE.data, *S=A.SIZE.data;

      for (i=0; i<IDX.dim1; ++i, I0+=r0, I+=r) {
         k=Wb::sub2ind(S0,I0,r0);
         Wb::ind2sub(k,S,I,r);
      }
   }

   return A;
};

template <class TD>
template <class DB> inline
wbsparray<TD>& wbsparray<TD>::init(
   const char *F, int L,
   SPIDX_T d1, SPIDX_T d2, const wbsparray<DB> &B
){
   SPIDX_T s[]= {d1,d2};
   wbvector<SPIDX_T> S; S.init(2,s);

   if (B.numel()!=S.prod(0)) wblog(F_L,"ERR %s() size mismatch "
      "(%s <> %s)",FCT,SSTR(S),SSTR(B));
   (*this)=B; return Reshape(S);
};

template <class TD>
template <class DB> inline
wbsparray<TD>& wbsparray<TD>::init(
   const char *F, int L,
   SPIDX_T d1, SPIDX_T d2, SPIDX_T d3, const wbsparray<DB> &B
){
   SPIDX_T s[]= {d1,d2,d3};
   wbvector<SPIDX_T> S; S.init(3,s);

   if (B.numel()!=S.prod(0)) wblog(F_L,"ERR %s() size mismatch "
      "(%s <> %s)",FCT,SSTR(S),SSTR(B));
   (*this)=B; return Reshape(S);
};

template <class TD>
template <class DB> inline
wbsparray<TD>& wbsparray<TD>::init(
   const char *F, int L,
   SPIDX_T d1, SPIDX_T d2, SPIDX_T d3, SPIDX_T d4, const wbsparray<DB> &B
){
   SPIDX_T s[]= {d1,d2,d3,d4};
   wbvector<SPIDX_T> S; S.init(4,s);

   if (B.numel()!=S.prod(0)) wblog(F_L,"ERR %s() size mismatch "
      "(%s <> %s)",FCT,SSTR(S),SSTR(B));
   (*this)=B; return Reshape(S);
};

template <class TD>                    
wbsparray<TD>& wbsparray<TD>::splitLast(SPIDX_T k, wbsparray &a) const {

   if (&a==this) {
      wblog(FL,"ERR %s() output same as input space",FCT);
   }

   unsigned r=SIZE.len, l=r-1;
   SPIDX_T *I=IDX.data+l, i=0, n=0;

   if (r<2) {
      wblog(FL,"ERR %s() got S=%s (k=%d)",FCT,SSTR_(this),k);
   }

   if (k>=SIZE[l]) wblog(FL,
      "ERR %s() index out of bounds (%s /%d)",FCT,SSTR_(this),k);
   if (IDX.dim2!=r) wblog(FL,"ERR %s() %d/%d",FCT,IDX.dim2,r);

   for (; i<IDX.dim1; ++i, I+=r) { if (*I==k) ++n; }

   a.SIZE.init(l,SIZE.data);
   a.init(a.SIZE,n);

   if (n) {
   for (I=IDX.data+l, n=i=0; i<IDX.dim1; ++i, I+=r) {
      if (*I==k) {
         a.IDX.recSetP(n,IDX.rec(i));
         a.D[n++]=D[i];
      }
   }}

   return a;
};

template <class TD>
wbsparray<TD>& wbsparray<TD>::select0(
   const WBINDEX &P, unsigned dim, 
   wbsparray<TD> &X 
) const {

   if (&X==this) {
      wbsparray<TD> Q; X.save2(Q);
      return Q.select0(P,dim,X);
   }
   if (dim>=SIZE.len) wblog(FL,
      "ERR %s() dimension out of bounds (%d/%d)",FCT,dim+1,SIZE.len);

   SPIDX_T i,l, d0=SIZE[dim], *j;
   sSPIDX_T *I;

   wbvector<sSPIDX_T> I_(d0); I_.set(-1); I=I_.data;

   for (i=0; i<P.len; ++i) {
      if (P[i]>=d0) wblog(FL,
         "ERR %s() index out of bounds (%d/%d)",FCT,P[i],d0);
      if (I[P[i]]>=0) wblog(FL,
         "ERR %s() input index not unique (%d/%d)",FCT,P[i],d0);
      I[P[i]]=i;
   }

   for (j=IDX.data+dim, l=i=0; i<IDX.dim1; ++i, j+=IDX.dim2) {
      if (I[*j]>=0) ++l;
   }; X.init(SIZE,l); X.SIZE[dim]=P.len;

   for (j=IDX.data+dim, l=i=0; i<IDX.dim1; ++i, j+=IDX.dim2) {
      if (I[*j]>=0) {
         X.IDX.recSetP(l,IDX.rec(i)); X.IDX(l,dim)=I[*j];
         X.D[l++]=D[i];
      }
   }

   return X;
};

template <class TD>
wbsparray<TD>& wbsparray<TD>::getBlock(const char *F, int L,
   SPIDX_T i1, SPIDX_T i2, SPIDX_T j1, SPIDX_T j2, 
   wbsparray<TD> &B 
) const {

   SPIDX_T k1,k2,dim1,dim2, i=0, l=0, n=IDX.dim1;

   if (!isRank(2)) wblog(FL,
      "ERR %s() only applies to rank-2 objects (%d)",FCT,SIZE.len);

   if (!SIZE.len) { dim1=dim2=D.len; } 
   else {
      if (SIZE.len!=2 || IDX.dim2!=2) wblog(FL,"ERR %s() "
         "got invalid rank-2 object (%ld,%ld)!?",FCT,SIZE.len,IDX.dim2);
      dim1=SIZE[0]; dim2=SIZE[1];
   }

   if (sSPIDX_T(i2)<0) i2+=dim1;
   if (sSPIDX_T(j2)<0) j2+=dim2;

   if (sSPIDX_T(i1)<0 || i2>=dim1) wblog(FL,"ERR %s() \n"
      "index out of bounds (1: %ld ..  %ld / %ld)", FCT,i1,i2,dim1);
   if (sSPIDX_T(j1)<0 || j2>=dim2) wblog(FL,"ERR %s() \n"
      "index out of bounds (1: %ld ..  %ld / %ld)", FCT,j1,j2,dim2);

   k1=(i2>=i1 ? i2-i1+1 : 0);
   k2=(j2>=j1 ? j2-j1+1 : 0);

   B.init(k1,k2); if (!k1 || !k2) return B;

   if (!SIZE.len) { 
      k1=MAX(i1,j1); k2=MIN(i2,j2); l=(k2>k1 ? k2-k1 : 0);
      if (l) { B.init_nnz(l); SPIDX_T *idx=B.IDX.data;
         for (; i<l; ++i, ++k1, idx+=2) {
            idx[0]=k1-i1; 
            idx[1]=k1-j1; 
            B.D.data[i]=D.data[k1];
         }
      }
      return B;
   }

   const SPIDX_T *i0=IDX.data; SPIDX_T *idx;

   for (; i<n; ++i, i0+=2) {
      if (i0[1]>=j1) {
         if (i0[1]<=j2) { if (i0[0]>=i1 && i0[0]<=i2) ++l; }
         else break; 
      }
   }

   if (l) { B.init_nnz(l); idx=B.IDX.data; i0=IDX.data; l=0;
      for (i=0; i<n; ++i, i0+=2) {
         if (i0[1]>=j1) {
            if (i0[1]<=j2) { if (i0[0]>=i1 && i0[0]<=i2) {
               idx[0]=i0[0]-i1; 
               idx[1]=i0[1]-j1; 
               B.D.data[l]=D.data[i]; idx+=2; ++l;
            }} else break; 
         }
      }
   }

   return B;
};

template <class TD>
wbsparray<TD>& wbsparray<TD>::splitBlock(const char *F, int L,
   sSPIDX_T i1, sSPIDX_T i2, sSPIDX_T j1, sSPIDX_T j2, 
   wbsparray<TD> &B 
){

   SPIDX_T k1,k2,dim1,dim2, j,i_, i=0, l=0, n=IDX.dim1;

   if (!isRank(2) || SIZE.len!=2) wblog(FL,
      "ERR %s() only applies to rank-2 objects (%d)",FCT,SIZE.len);

   B.init(SIZE); if (!D.len) return B;
   dim1=SIZE[0]; dim2=SIZE[1];

   if (i2<0) i2+=dim1;
   if (j2<0) j2+=dim2;

   if (i1<0 || i2>=dim1) wblog(FL,"ERR %s() \n"
      "index out of bounds (1: %ld ..  %ld / %ld)", FCT,i1,i2,dim1);
   if (j1<0 || j2>=dim2) wblog(FL,"ERR %s() \n"
      "index out of bounds (1: %ld ..  %ld / %ld)", FCT,j1,j2,dim2);

   k1=(i2>=i1 ? i2-i1+1 : 0);
   k2=(j2>=j1 ? j2-j1+1 : 0); if (!k1 || !k2) return B;

   SPIDX_T *j0, *i0=IDX.data, *idx;

   for (; i<n; ++i, i0+=2) { if (i0[1]>=j1) break; }; i_=i;
   for (; i<n; ++i, i0+=2) { 
      if (i0[1]<=j2) { if (i0[0]>=i1 && i0[0]<=i2) ++l; }
      else break;
   }

   if (l) {
      B.init_nnz(l); idx=B.IDX.data; l=0; j=i=i_; i0=IDX.data+2*i;
      for (; i<n; ++i, i0+=2) {
         if (i0[1]<=j2 && i0[0]>=i1 && i0[0]<=i2) {
            idx[0]=i0[0];
            idx[1]=i0[1]; idx+=2;
            B.D.data[l]=D.data[i]; ++l;
         }
         else {
            if (j<i) { j0=IDX.data+2*j;
               j0[0]=i0[0];
               j0[1]=i0[1]; D[j]=D[i];
            }; ++j;
         }
      }
   }

   REDSIZE_NNZ(j);
   return B;
};

template <class TD>
TD wbsparray<TD>::getHouseholderVec(const char *F, int L,
   SPIDX_T k,        
   SPIDX_T j,        
   wbsparray<TD> &u, 
   TD eps2
) const {

   TD x2=0, b2=0; 

   if (!isRank(2) || SIZE.len!=2 || IDX.dim2!=2) wblog(F_L,
      "ERR %s() only applies to rank-2 objects (%s; %d,%d)",
      FCT,SSTR_(this),SIZE.len,IDX.dim2);

   SPIDX_T dim1=SIZE[0], dim2=SIZE[1], i1,i2, i=0, l=0, n=IDX.dim1;
   const SPIDX_T *i0=IDX.data;
   char gotk=0;

   if (k>=dim2) wblog(F_L,"ERR %s() "
      "col-index out of bounds (%d/%d)",FCT,k,dim2);
   u.init(dim1); if (!dim1) return x2;

   for (; i<n; ++i, i0+=2) { if (i0[1]>=j) break; }
   for (; i<n; ++i, i0+=2) { if (i0[0]>=k) break; }
   i1=i; if (i<n && i0[0]==k && i0[1]==j) gotk=1;
   for (; i<n; ++i, i0+=2) {
      if (i0[1]==j) { b2+=NORM2(D[i]); ++l; }
      else break;
   }
   i2=i;

   if (b2<=eps2) {
      if (b2==0) {
         wblog(FL,"WRN %s() got 0 vector (%s) !?",FCT,SSTR_(this));
         return x2;
      }
      wblog(FL,"WRN %s() got b2=%.3g !?",FCT,double(b2));
   }

   u.init_nnz(l+(gotk?0:1));

   SPIDX_T *idx=u.IDX.data; i0=IDX.rec(i1);

   if (gotk) { l=0; } else { u.D[0]=0; idx[0]=k; l=1; }

   for (i=i1; i<i2; ++i, i0+=2) {
      u.D[l]=D.data[i]; idx[l]=i0[0]; ++l;
   }

   x2=NORM2(u.D[0]);
   if (u.D[0]<0) { u.D[0]-=SQRT(b2); } else { u.D[0]+=SQRT(b2); }
   x2=b2-x2+NORM2(u.D[0]);

   u*=(TD(1)/SQRT(x2));

   return x2;
};

template <class TD>
wbsparray<TD>& wbsparray<TD>::IDX_Shift(sSPIDX_T i1, sSPIDX_T i2) {

   if (!isRank(2) || SIZE.len!=2 || IDX.dim2!=2) wblog(FL,
      "ERR %s() only applies to rank-2 objects (%s; %d,%d)",
      FCT,SSTR_(this),SIZE.len,IDX.dim2
   );

   SPIDX_T i=0, n=IDX.dim1, *I=IDX.data;
   SPIDX_T dim1=SIZE[0]+i1, dim2=SIZE[1]+i2; SIZE[0]=dim1; SIZE[1]=dim2;

   if (i1>=0 && i2>=0) {
      for (; i<n; ++i, I+=2) { I[0]+=i1; I[1]+=i2; }
   }
   else { 
      if (sSPIDX_T(dim1)<0) wblog(FL,
         "ERR %s() index out of bounds (%ld/%ld)",FCT,i1,SIZE[0]);
      if (sSPIDX_T(dim2)<0) wblog(FL,
         "ERR %s() index out of bounds (%ld/%ld)",FCT,i2,SIZE[1]);

      for (; i<n; ++i, I+=2) { I[0]+=i1; I[1]+=i2;
          if (I[0]>=dim1 || I[1]>=dim2) wblog(FL,
             "ERR %s() index out of bounds (%ld,%ld; %s)",
             FCT,I[0],I[1],SSTR_(this)
          );
      }
   }

   return *this;
};

template <class TD>
SPIDX_T wbsparray<TD>::LContractVec_omp(
   const wbsparray<TD> &b, wbvector<TD> &x, char tnorm) const {

   SPIDX_T i, k0=0, d2, l=0; 

   if (!isMatrix()) wblog(FL,
      "ERR %s() only applies to rank-2 objects (%d)",FCT,SIZE.len);

   if (!isEmpty() && !b.isEmpty()) {
      if (b.IDX.dim2!=1 || b.SIZE.len!=1) wblog(FL,
         "ERR %s() expecting sparse VECTOR\n(got b: IDX=%s; SIZE=[%s])",
         FCT,SSTR(b.IDX),STR(b.SIZE));
      if (dim0(0)!=b.SIZE[0]) wblog(FL,
         "ERR %s() size mismatch (%d/%d)",FCT,SSTR_(this),SSTR(b));
   }

   if (!D.len || !b.D.len) { x.init(); return l; }

   if (SIZE.len) {
      indexSparseRef<TD> ISP(FL,*this);
      d2=ISP.numCols(); x.init(d2);

     #ifndef WB_SKIP_ASSERT
      SPIDX_T d1=ISP.numRows();
      if (d1!=SIZE[0] || d2!=SIZE[1] || !d1 || !d2) wblog(FL,
         "ERR %s() got size inconsistency (%dx%d / %s)",
         FCT,d1,d2,SSTR_(this));
     #endif

      int np=1;
      if (sp_num_threads>1) {
         np=IDX.dim1/1024;
         if (np>sp_num_threads) np=sp_num_threads;
         if (np>int(d2)) np=d2;
      }

      for (; k0<d2; ++k0) { if (ISP.nnz(k0)) break; }

     #pragma omp parallel for num_threads(np) reduction(+:l)
      for (SPIDX_T k=k0; k<d2; ++k) {
         x[k]=ISP.LContractVec(b,k,tnorm);
         if (x[k]) l+=1;
      }

      if (l!=x.nnz()) {  
         MXPut(FL,"q").add(x,"x").add(l,"l"); SPIDX_T n=0;
         for (i=0; i<x.len; ++i) {
            if (x[i]) printf("%6ld: %3ld : %+.4g\n",i,++n,double(x[i])); else
            if (x[i]!=0) printf("%6ld: %3ld : %+.4g (!!)\n",i,++n,double(x[i]));
         }
         wblog(FL,"ERR %s() %ld/%ld !?",FCT,l,x.nnz());
      }
   }
   else {
      const SPIDX_T *idx=b.IDX.data;
      d2=D.len; 
      x.init(d2); k0=b.IDX[0]; l=b.IDX.dim1;
      for (i=0; i<l; ++i) {
         if (idx[i]>=d2 || (i && idx[i-1]>=idx[i])) wblog(FL,
            "ERR %s() b.IDX out of bounds (%ld/%ld) !?",FCT,idx[i],d2);
         x[idx[i]]=b.D[i]*D[idx[i]];
      }
   }

   return l;
};

template <class TD>
wbsparray<TD>& wbsparray<TD>::LProject1(
   const wbsparray<TD> &B, wbsparray<TD> &X, TD xfac, char tnorm) const {

   if (!D.len) { X=(*this); return X; }
   if (!B.D.len) { X.init(SIZE); return X; }

   wbvector<TD> x;
   SPIDX_T i, d1,d2;
   SPIDX_T l=LContractVec_omp(B,x,tnorm); 

   if (xfac!=1) {
      if (xfac) x*=xfac;
      else wblog(FL,"ERR %s() got xfac=0 !?",FCT);
   }

   d1=B.D.len; d2=x.len;
   if (SIZE.len) { X.init(SIZE, l*d1); }
   else {
      SPIDX_T s[2]={d2,d2};
      X.init(wbvector<SPIDX_T>(2,s), l*d1);
   }

   const SPIDX_T *I=B.IDX.data;
   SPIDX_T *idx=X.IDX.data, r=X.IDX.dim2; l=0;
   const TD *b=B.D.data; TD *xb=X.D.data;

   for (SPIDX_T k=0; k<d2; ++k) { if (x[k]) {
      for (i=0; i<d1; ++i, ++l, idx+=r) {
         xb[l]=x[k]*b[i]; idx[0]=I[i]; idx[1]=k;
      }
   }}

   if (l!=X.IDX.dim1) wblog(FL,"ERR %s() %ld/%ld !?",FCT,l,X.IDX.dim1);
   return X;
};

template <class TD> inline
wbsparray<TD>& wbsparray<TD>::Permute(const wbperm &P, char iflag){

   if (!SIZE.len && !IDX.dim1 && D.len<=1) return *this;
   if (!P.isEmpty()) {
      SIZE.Permute(P,iflag); IDX.colPermute(P,iflag);
      Sort();
   }
   return *this;
};

template <class TD> inline
wbsparray<TD>& wbsparray<TD>::permute(
   const wbperm &P, wbsparray<TD>&B, char iflag
) const {

   if (this==&B) return B.Permute(P,iflag);

   if (P.isEmpty() || (!SIZE.len && !IDX.dim1 && D.len<=1)) { B=(*this); }
   else {
      SIZE.permute(B.SIZE,P,iflag); B.D=D; B.isref=0;
      IDX.colPermute(P,B.IDX,iflag);
      B.Sort();
   }
   return B;
};

template <class TD> inline
wbsparray<TD>& wbsparray<TD>::Permute(const char* s, char iflag){
   return Permute(Str2Idx(s,1),iflag);
};

template <class TD> inline
wbsparray<TD>& wbsparray<TD>::permute(
   const char* s, wbsparray<TD>&Q, char iflag
) const { return permute(Q,Str2Idx(s,1),iflag); };

template <class TD> inline
wbsparray<TD>& wbsparray<TD>::ColPermute(const wbperm &P0, char iflag){

   if (isIdentityPerm(P0)) return *this;

   if (SIZE.len!=2 || IDX.dim2!=2) wblog(FL,
      "ERR %s() requires rank-2 (%s; %d)",FCT,SSTR_(this),IDX.dim2);
   if (P0.len!=SIZE[1]) wblog(FL,
      "ERR %s() size mismatch (%d/%d)",FCT,P0.len,SIZE[1]);

   SPIDX_T k,i=0, l=IDX.dim2-1; wbperm P;

   if (!iflag)
        P0.invert(P);
   else P.wbvector<wperm_t>::init2ref(P0);

   for (; i<IDX.dim1; ++i) { k=IDX(i,l);
      if (k>=P.len) wblog(FL,"ERR %s() %d/%d !?",FCT,k,P.len);
      IDX(i,l)=P[k];
   }
   return Sort();
};

template <class TD> inline
wbsparray<TD>& wbsparray<TD>::MatPermute(const wbperm &P0, char iflag){

   if (isIdentityPerm(P0)) return *this;

   SPIDX_T d=-1;
   if (!isSMatrix(FL,&d)) wblog(FL,
      "ERR %s() requires rank-2 (%s; %d)",FCT,SSTR_(this),IDX.dim2);
   if (P0.len!=d) wblog(FL,
      "ERR %s() size mismatch (%d/%d)",FCT,P0.len,d);

   if (isDiag()) { D.Permute(P0); }
   else {
      wbperm P;
      if (!iflag) 
           P0.invert(P);
      else P.wbvector<wperm_t>::init2ref(P0);

      SPIDX_T *I=IDX.data, n=IDX.numel(); wperm_t *p=P.data;
      for (SPIDX_T i=0; i<n; ++i) {
         if (I[i]>=P.len) wblog(FL,"ERR %s() %d/%d !?",FCT,I[i],P.len);
         I[i]=p[I[i]];
      }
   }

   return Sort();
};

template <class TD>
void wbsparray<TD>::setSIZE_kron( 
   const char *F, int L,
   const wbvector<SPIDX_T> &sa, const wbvector<SPIDX_T> &sb,
   char kflag
){
   wbvector<SPIDX_T> S;
   if (sa.len!=sb.len) wblog(F_L,
      "ERR %s() rank mismatch (%d/%d)",FCT,sa.len,sb.len);

   if (kflag) { S.init(sa.len); SPIDX_T *s=S.data;
      for (unsigned i=0; i<sa.len; ++i) {
         s[i]=sa[i]*sb[i]; 
      }
   }
   else { S.init(sa.len+sb.len); SPIDX_T *s=S.data;
      for (unsigned i=0, l=0; i<sa.len; ++i) {
         s[l++]=sa[i];
         s[l++]=sb[i];
      }
   }
   S.save2(SIZE);
};

template <class TD> 
void wbsparray<TD>::setSIZE_kron(const char *F, int L, char kflag) {

   if (kflag) {
      unsigned l=0, i=0, r=IDX.dim2; SPIDX_T *s=SIZE.data;
      if (SIZE.len!=2*r) wblog(F_L,
         "ERR %s() unexpected size (%dx%d <> %d/(%d/2))",
         FCT,IDX.dim1,IDX.dim2,D.len,SIZE.len
      );
      if (r) {
         for (; i<r; ++i,l+=2) { s[i]=s[l]*s[l+1]; }
         SIZE.len=r;
      }
   }
   else if (SIZE.len!=IDX.dim2) sperror_this(F_L);
};

template <class TD>
wbsparray<TD>& wbsparray<TD>::setRec_kron(
   SPIDX_T l, TD x, unsigned r, SPIDX_T *i1, SPIDX_T *i2,
   unsigned kflag
){
   if (l>IDX.dim1) wblog(FL,"ERR %s() index "
      "out of bounds (%d,%d <> %dx%d)",FCT,l,r,IDX.dim1,IDX.dim2);

   if (kflag) {
      if (l==0 && (SIZE.len!=2*r || IDX.dim2!=r))
         wblog(FL,"ERR %s() 2*%d/2*%d/%d",FCT,r,IDX.dim2,SIZE.len);
      SPIDX_T *q=IDX.rec(l);
      for (unsigned i=0; i<r; ++i) {
         q[i]=i1[i] + SIZE[2*i]*i2[i];
      }
   }
   else {
      if (l==0 && (SIZE.len!=2*r || SIZE.len!=IDX.dim2))
         wblog(FL,"ERR %s() index out of bounds (%d,%d <> %dx%d)",
         FCT,l,r,IDX.dim1,IDX.dim2
      );

      SPIDX_T *q=IDX.rec(l);
      for (unsigned j=0, i=0; i<r; ++i) {
          q[j++]=i1[i];
          q[j++]=i2[i];
      }
   }
   D[l]=x; return *this;
};

template <class TD>
wbsparray<TD>& wbsparray<TD>::transpose(
   const char *F, int L, wbsparray<TD> &A
 ) const {

   unsigned r=rank(); wbperm P;

   if (r%2) wblog(F_L,"ERR %s() "
      "applies to even-rank arrays only! (%s)",FCT,SSTR_(this));

   P.initTranspose(r);
   permute(P,A);

   return A.Sort();
};

template <class TD> inline
wbsparray<TD>& wbsparray<TD>::tensorProd(const char *F, int L,
   const wbsparray<TD>& B0, wbsparray<TD>& C,
   char aflag0, char bflag0,
   char kflag 
 ) const {

   char iA=isDiag(FL), iB=B0.isDiag(FL);

   if (iA || iB) {
      if (iA && iB) {
          if (kflag) {
             C.initDiag(D.len*B0.D.len); D.tensorProd(B0.D,C.D);
             return C;
          }
          wbsparray<TD> A2;    diag2reg(FL,A2);
          wbsparray<TD> B2; B0.diag2reg(FL,B2);
          return A2.tensorProd(F,L,B2,C,aflag0,bflag0,kflag);
      }
      else if (iA) {
          wbsparray<TD> A2; diag2reg(FL,A2, B0.rank());
          return A2.tensorProd(F,L,B0,C,aflag0,bflag0,kflag);
      }
      else { 
          wbsparray<TD> B2; B0.diag2reg(FL,B2, rank());
          return tensorProd(F,L,B2,C,aflag0,bflag0,kflag);
      }
   }

   if ((void*)this==(void*)(&C) || (void*)(&B0)==(void*)(&C)) {
       wbsparray<TD> X; tensorProd(F_L,B0,X,aflag0,bflag0,kflag);
       return X.save2(C);
   }

   unsigned ra=rank(), rb=B0.rank();

   if (!ra || ra%2 || ra!=rb) wblog(F_L,"ERR %s() invalid input to "
      "kron(%s,%s)",FCT,SSTR_(this),SSTR(B0));
   if (SIZE.len!=IDX.dim2 || B0.SIZE.len!=B0.IDX.dim2) wblog(FL,
      "ERR %s() data inconsistency (%s -> %dx%d, %s -> %dx%d)",FCT,
         SSTR_(this),   IDX.dim1,   IDX.dim2, 
      SSTR(B0),B0.IDX.dim1,B0.IDX.dim2);

   opFlags<TD> aflag(aflag0), bflag(bflag0);
   wbsparray<TD> A(*this),B(B0);

   if (!A.isCompressed()) wblog(FL,"ERR %s() A not sorted/unique",FCT); 
   if (!B.isCompressed()) wblog(FL,"ERR %s() B not sorted/unique",FCT);

   aflag.apply(F_L,A);
   bflag.apply(F_L,B);

   SPIDX_T i,j, l=0, na=A.IDX.dim1, nb=B.IDX.dim1;
   C.init_kron(SIZE,B.SIZE,na*nb,kflag);

   if (C.IDX.dim1) {
      for (i=0; i<na; ++i)
      for (j=0; j<nb; ++j,++l) {
         C.setRec_kron(
           l, A.D[i]*B.D[j], ra, A.IDX.rec(i), B.IDX.rec(j), kflag
         );
      }
   }

   C.setSIZE_kron(F_L,kflag); 
   C.Compress(FL); 

   return C;
};

template <class TD> inline
wbsparray<TD>& wbsparray<TD>::tensorProdX(const char *F, int L,
   const wbsparray& B, wbsparray& C,
   char conja, const TD &afac, char conjb, const TD &cfac,
   const wbperm *Pfinal
 ) const {

   if (cfac && cfac!=1) { C.D*=cfac; }
   if (afac==0) { return C; }

   char iA=isDiag(FL), iB=B.isDiag(FL), noconj=1;
   unsigned l=0, ra=(iA ? 2 : SIZE.len), rb=(iB ? 2 : B.SIZE.len);

   const SPIDX_T *Ia=IDX.data, *Ib=B.IDX.data;
   SPIDX_T ia=0, ib=0, ic=0, na=D.len, nb=B.D.len, *Ic;
   wbsparray<TD> X;

   if ((conja || conjb)) {
      if (ISCOMPLX_(TD)) { noconj=0; }
      else { conja=conjb=0; }
   }

   X.SIZE.init(ra+rb); Ic=X.SIZE.data;
      if (iA) { Ic[0]=Ic[1]=na; }
      else { for (; l<ra; ++l) { Ic[l]=SIZE[l]; }; l=0; }; Ic+=ra;
      if (iB) { Ic[0]=Ic[1]=nb; }
      else { for (; l<rb; ++l) { Ic[l]=B.SIZE[l]; }}
   X.init_nnz(na*nb); Ic=X.IDX.data;

   for (; ib<nb; ++ib) {
      for (ia=0; ia<na; ++ia, ++ic) { 
         if (iA) { Ic[0]=Ic[1]=ia; }
         else { for (l=0; l<ra; ++l) { Ic[l]=Ia[l]; }; Ia+=ra; }
         Ic+=ra;

         if (iB) { Ic[0]=Ic[1]=ib; }
         else { for (l=0; l<rb; ++l) { Ic[l]=Ib[l]; }; Ib+=rb; }
         Ic+=rb;

         if (noconj) { X.D[ic]=D[ia]*B.D[ib]; }
         else { X.D[ic] =
                (conja ? CONJ(  D[ia]) :   D[ia])
              * (conjb ? CONJ(B.D[ib]) : B.D[ib]);
         }
      }
   }

   if (Pfinal && !Pfinal->isIdentityPerm()) { X.Permute(*Pfinal); }

   if (afac!=1) { X*=afac; }
   if (cfac) {
      if (X.SIZE!=C.SIZE) wblog(FL,"ERR %s() "
         "got size mismatch (%s <> %s) !?",FCT,SSTR(X),SSTR(C));
      C.Plus(F_L,X);
   }
   else X.save2(C);

   return C;
};

template <class TD>
wbsparray<TD>& wbsparray<TD>::Plus(
   const char *F, int L,
   const wbsparray<TD> &B, TD bfac, TD afac,
   TD eps 
){
   if (isref) wblog(FL,"ERR %s() got isref=%d !?",FCT,isref);

   if ((void*)this==(void*)(&B)) {
      TD x=afac+bfac;
      if (x== 0) { D.init(); IDX.init(0,IDX.dim2); isref=0;    } else
      if (x==-1) { for (SPIDX_T i=0; i<D.len; ++i) D[i]=-D[i];} else
      if (x!=+1) { for (SPIDX_T i=0; i<D.len; ++i) D[i]*=x;   }
      return *this;
   }

   if (isEmpty()) { (*this)=B; (*this)*=bfac; return *this; }
   if (B.isEmpty()) { (*this)*=afac; return *this; }

   char ida=isDiag(), idb=B.isDiag();
   if (ida || idb) { 
      if (ida && idb) {
         if (D.len!=B.D.len) wblog(FL,
            "ERR %s() size mismatch (%ld/%dl)",FCT,D.len,B.D.len);
         for (SPIDX_T i=0; i<D.len; ++i) { D[i]=afac*D[i]+bfac*B.D[i]; }
         return *this;
      }
      else if (ida) { diag2reg();
         return Plus(F,L,B,bfac,afac,eps);
      }
      else if (idb) {
         wbsparray<TD> B2; B.diag2reg(FL,B2);
         return Plus(F,L,B2,bfac,afac,eps);
      }
   }

   if (SIZE!=B.SIZE || !SIZE.len) wblog(FL,
      "ERR %s() got size mismatch (%s <> %s)",FCT,SSTR_(this),SSTR(B));

   if (afac!=1) {
      if (afac== 0) { D.init(); IDX.init(0,IDX.dim2); isref=0; } else
      if (afac==-1) { for (SPIDX_T i=0; i<D.len; ++i) D[i]=-D[i];} else
      if (afac!=+1) { for (SPIDX_T i=0; i<D.len; ++i) D[i]*=afac;}
   }

   if (bfac==0) return *this;

   SPIDX_T i=0, ia=0, ib=0, l=0, nx, na=IDX.dim1, nb=B.IDX.dim1;
   const SPIDX_T *Ia=IDX.data, *Ib=B.IDX.data;
   unsigned m=SIZE.len;
   char lex=0; 

   wbsparray<TD> X(SIZE);
   const TD *a=D.data, *b=B.D.data; TD *x, e2=0;

   wbvector<char> qq_(na+nb); char q, *qq=qq_.data;

#ifndef WB_SKIP_ASSERT
   SPIDX_T Dmax=SIZE.max();
   for (ia=1; ia<na; ia+=Dmax) {
      if (Wb::recCompare(Ia+(ia-1)*m,Ia+ia*m,m,lex)>=0) wblog(FL,
      "ERR %s() A: input sparse IDX not sorted (%ld/%ld) !?",FCT,ia,na);
   }
   for (ib=1; ib<nb; ib+=Dmax) {
      if (Wb::recCompare(Ib+(ib-1)*m,Ib+ib*m,m,lex)>=0) wblog(FL,
      "ERR %s() B: input sparse IDX not sorted (%ld/%ld) !?",FCT,ib,nb);
   }
   ia=ib=0;
#endif

   while (ia<na && ib<nb) {
      qq[l++]=q=Wb::recCompare(Ia+ia*m,Ib+ib*m,m,lex);

      if (q<=0) { ++ia; }
      if (q>=0) { ++ib; }
   }

   if (ia<na)
        { for (; ia<na; ++ia) { qq[l++]=-1; } }
   else { for (; ib<nb; ++ib) { qq[l++]=+1; } }

   X.init_nnz(l); x=X.D.data; nx=X.IDX.dim1; ia=ib=l=0;

   for (; i<nx; ++i) {
      if (qq[i]<0) {
         x[l]=a[ia]; if (afac!=1) x[l]*=afac;
         if (ABS(x[l])>eps) X.IDX.recSetP(l++,Ia); else e2+=NORM2(x[l]);
         ++ia; Ia+=m;
      }
      else if (qq[i]>0) {
         x[l]=b[ib]; if (bfac!=1) x[l]*=bfac;
         if (ABS(x[l])>eps) X.IDX.recSetP(l++,Ib); else e2+=NORM2(x[l]);
         ++ib; Ib+=m;
      }
      else {
         x[l]=a[ia]; if (afac!=1) x[l]*=afac;
         x[l]+=(bfac!=1 ? bfac*b[ib] : b[ib]);
         if (ABS(x[l])>eps) X.IDX.recSetP(l++,Ia); else e2+=NORM2(x[l]);
         ++ia; ++ib; Ia+=m; Ib+=m;
      }
   }

   X.REDSIZE_NNZ(l);

   return X.save2(*this);
};

template <class TD>
double wbsparray<TD>::Compress( 
   const char *F, int L, TD eps, char lex, const wbperm *Pc
){
   wbperm P; WBINDEX dg; rank(F_L); 

   checkSize(F_L);
   if (!IDX.dim1) { 
      double e2=0;
      if (Pc && Pc->len!=2) wblog(FL,
         "ERR %s() invalid permutation (%s)",FCT,Pc->toStr().data);
      if (double(eps)>0) {
         for (SPIDX_T k=0; k<D.len; ++k) {
            TD &x=D.data[k];
            if (x && ABS(x)<=eps) { e2+=double(NORM2(x)); x=0; }
         }
      }
      return e2;
   }

   if (!Pc || Pc->isIdentityPerm()) {
      if (double(eps)<=0 && IDX.isUniqueSorted(+1,lex)) { return 0; }
      IDX.groupRecs(P,dg,-1,lex);
   }
   else {
      wbMatrix<SPIDX_T> X; IDX.colPermute(*Pc,X);
      X.groupRecs(P,dg,-1,lex);
   }

#ifdef WB_SPARSE_CLOCK
   Wb::Clock clk("sparse:cmpr",0); 
#endif

   SPIDX_T i,j,k,d; wperm_t *p=P.data;
   wbvector<TD> D0; D.save2(D0); D.init(dg.len);
   const TD *d0=D0.data;

   for (j=k=0; j<dg.len; ++j,++k) {
      TD &x=D.data[k]; x=d0[p[0]];
      for (d=dg.data[j], i=1; i<d; ++i) { x+=d0[p[i]]; }

      p+=d;
      if (x==0) { --k; } else if (k<j) { IDX.recSet(k,j); }
   }

   REDSIZE_NNZ(k);

   return (double)SkipTiny(eps); 
};

template <class TD>
int wbsparray<TD>::splitSparseCM(
   const char *F, int L, wbvector< wbsparray<TD> > &X,
   TD eps,     
   char ref,   
   char tnorm
){
   iterSparseCMref<TD> ISP(F_L,*this);
   SPIDX_T i=0, n=ISP.numIter(); if (!n) { X.init(); return X.len; }

   if (eps>0) {
      wbvector<char> mark(n); SPIDX_T l=0;
      for (; i<n; ++i) { if (ISP.norm2(i,tnorm)>=eps) { mark[i]=1; ++l; }}
      X.init(l); i=l=0;
      for (; i<n; ++i) { if (mark[i]) { ISP.getIter(i,X[l++],ref); }}
   }
   else {
      X.init(n);
      for (; i<n; ++i) { ISP.getIter(i,X[i],ref); }
   }

   return X.len;
};

template <class T>
inline int checkContract(const char* F, int L,
   const wbvector<T> &SA, const ctrIdx &ica,
   const wbvector<T> &SB, const ctrIdx &icb,
   const wbperm *pfinal, wbvector<T> *SC
){
   int e=0;

   if (ica.len!=icb.len) {
      if (F) wblog(F,L,
         "ERR %s() invalid index set (length mismatch, [%s], [%s])",
          FCT,STR(ica),STR(icb));
      return 1;
   }
   if (!ica.len) {
      if (F) wblog(F,L,"ERR contract() got empty ic[ab] !?");
      return 1;
   }

   for (unsigned i=0; i<ica.len; ++i) {
      if ((ica[i]<SA.len ? SA[ica[i]] : 1) !=
          (icb[i]<SB.len ? SB[icb[i]] : 1))
      wblog(FL, 
         "ERR invalid contraction @ i=%d/%d ... !?\n"
         "having: %s (%s) * %s (%s)",i+1,ica.len,
         SSTR(SA), STR(ica), SSTR(SB), STR(icb)
      );
   }

   if (SC) {
      wbvector<T> sC; sC.initI(SA,ica,SB,icb);

      if (SC->isEmpty()) { sC.save2(*SC); }
      else {
         if (pfinal && pfinal->len) sC.Permute(*pfinal);
         if ((*SC)!=sC) { e=6; if (F) wblog(FL,
           "ERR %s() severe size mismatch (C: %s <> %s)",
            FCT, SSTR(sC), SSTR_(SC));
         }
      }

      if (!SC->prod(0)) return -1;
   }

   return e;
};

template <class TD> inline
char wbsparray<TD>::contract_scalar(
   const char* F, int L,  const ctrIdx &ica,
   const wbsparray<TD> &B, const ctrIdx &icb, wbsparray<TD> &C,
   const wbperm &pfinal, const TD& afac, const TD& cfac
 ) const {

   char isa=isScalar(), isb=B.isScalar(), gotC=(cfac && !C.isEmpty());
   unsigned r=-1;
   if (!isa || !isb) return 0;

   TD x=(D.len && B.D.len ? afac*D[0]*B.D[0] : TD(0));
   if (gotC && C.D.len) { x+=cfac*C.D[0]; }

   if (ica.len!=icb.len) wblog(F_L,"ERR %s() length mismatch: "
      "%s <> %s",FCT, STR(ica), STR(icb));
   if (gotC && (!C.isScalar() || C.D.len>1)) wblog(F_L,
      "ERR %s() invalid scalar %s",FCT,C.info2Str("C").data);
   if (SIZE.len!=IDX.dim2 || B.SIZE.len!=B.IDX.dim2) wblog(F_L,
      "ERR %s() invalid scalars %s <> %s", FCT,
      info2Str("A").data, B.info2Str("B").data);

   if (SIZE.len) {
      if (B.SIZE.len) {
         checkContract(F_L, SIZE,ica,B.SIZE,icb,
            (gotC? &pfinal : NULL), (gotC? &C.SIZE : NULL));
         r=(SIZE.len + B.SIZE.len - 2*ica.len);
      }
      else {
         wbvector<SPIDX_T> Sb(2); Sb[0]=Sb[1]=B.D.len;
         checkContract(F_L, SIZE,ica,Sb,icb,
            (gotC? &pfinal : NULL), (gotC? &C.SIZE : NULL));
         r=(SIZE.len + 2 - 2*ica.len);
      }
   }
   else {
      if (B.SIZE.len) {
         wbvector<SPIDX_T> Sa(2); Sa[0]=Sa[1]=D.len;
         checkContract(F_L,Sa,ica,B.SIZE,icb,
            (gotC? &pfinal : NULL), (gotC? &C.SIZE : NULL));
         r=(2 + B.SIZE.len - 2*ica.len);
      }
   }

   C.initScalar(x,r);

   return 1;
};

template <class TD>
wbsparray<TD>& wbsparray<TD>::contract_diag_diag(
   const char* F, int L,  const ctrIdx &ica,
   const wbsparray<TD> &B, const ctrIdx &icb, wbsparray<TD> &C,
   const wbperm &pfinal, const TD& afac, const TD& cfac
 ) const {

   SPIDX_T i=0, ra=2, rb=2; int e=0;

   if (SIZE.len || B.SIZE.len) wblog(F_L,
      "ERR %s() (%s) <> (%s)",FCT,SSTR_(this),SSTR(B));
   if (D.len!=B.D.len) wblog(FL,
      "ERR %s() size mismatch (%d/%d)",FCT,D.len,B.D.len);

   if  (!ica.len) e=-1; else 
   if  (ica.len!=icb.len || ica.len>2) e=-2; else {
      for (; i<ica.len; ++i) if (ica[i]>=ra || icb[i]>=rb) e=3;
   }
   if (!e && ica.len==2 && (ica[0]==ica[1] || icb[0]==icb[1])) e=4;
   if (e) wblog(F_L,
      "ERR %s() invalid contraction %s (%s) <> %s (%s)",
      FCT,SSTR_(this),STR(ica),SSTR(B),STR(icb)
   );

   if (ISCOMPLX_(TD) && (ica.conj || icb.conj)) 
   wblog(FL,"ERR %s() complex not yet implemented",FCT);  

   wbvector<TD> x(D.len);
   if (afac) {
      for (i=0; i<D.len; ++i) { x[i]=D[i]*B.D[i]; }
   }

   if (ica.len==2) { TD q=afac*x.sum(); 
      if (cfac) {
         if (!C.isScalar()) wblog(F_L,
            "ERR %s() invalid scalar C (%s)",FCT,SSTR(C));
         C.D*=cfac; C[0]+=q;
      }
      else C.initScalar(q);

      if (!pfinal.isEmpty()) wblog(F_L,
         "ERR %s() got pfinal=[%s] for scalar !?",FCT,STR(pfinal)
      );
   }
   else { 
      if (cfac) {
         if (!C.isDiag() || C.D.len!=D.len) wblog(F_L,
            "ERR %s() invalid diag C (%s; %d) @ cfac=%.3g",
            FCT,SSTR(C),D.len,double(cfac));
         C.D*=cfac; C.D.Plus(x,afac);
      }
      else {
         if (afac) x*=afac;
         C.initDiag(D.len); x.save2(C.D);
      }

      if (!pfinal.isEmpty() && !pfinal.isValidPerm(0,0,2))
      wblog(F_L,"ERR %s() got pfinal=[%s] !?",FCT,STR(pfinal));
   }

   return C;
};

template <class TD>
char wbsparray<TD>::contract_check_2full(
   const char* F, int L,  const ctrIdx &ica,
   const wbsparray<TD> &B, const ctrIdx &icb, wbsparray<TD> &C,
   const wbperm &pfinal, const TD& afac, const TD& cfac,
   char flag, TD eps, SPIDX_T *nnz_C
) const {

   if (flag=='s' || flag=='S'
       || typeid(TD)!=typeid(double)  
   ) return 0; 

   if (flag && flag!='f' && flag!='F')   
   wblog(F_L,"ERR %s() invalid flag=%c<%d>",FCT,flag,flag);

   if (!flag || nnz_C!=NULL) {
      SPIDX_T NA=numel(), NB=B.numel(), M=SIZE.prod(ica.data,ica.len);
      SPIDX_T NC=(NA/M)*(NB/M); 
      double pa=nnz()/double(NA), pb=B.nnz()/double(NB);

      double pc=M*pa*pb; if (pc>1) pc=1; 

      if (!NA || !NB || !M) wblog(F_L,
         "ERR %s() got\n%s [%s] <> %s [%s] (%d)",FCT,SSTR_(this),
         ica.toStr().data, SSTR(B), STR(icb), M);

      if ((pc>0.50 && NC<(1<<24)) || 
          (pc>0.20 && NC<(1<<16)) || 
          (pc>0.10 && NC<(1<<12))    
       ) { flag='F';
         wbvector<SPIDX_T> sC; sC.initI(SIZE,ica,B.SIZE,icb);
      }

      if (nnz_C) {
         double nc=pc*NC;
         if (NC>128) (nc)*=1.20; else
         if (NC> 32) (nc)*=1.50; else
         if (NC>  4) (nc)*=2; else nc=NC;
         (*nnz_C)=MIN(SPIDX_T(nc),NC);
      }

   }

   if (!flag) return 0;

#ifdef WB_SPARSE_CLOCK
   Wb::Clock clk("sparse:ctr-full",0); 
#endif

   wbarray<TD> Af,Bf,ABf; this->toFull(Af); B.toFull(Bf);

   Af.contract(F_L,ica,Bf,icb,ABf,pfinal,afac);

   if (cfac!=0 && C.D.len) { C*=cfac;
      wbarray<TD> Cf; C.toFull(Cf);
      ABf+=Cf;
   }

   C.init(F_L,ABf,eps); return 1;
};

template <class TD>
wbsparray<TD>& wbsparray<TD>::contract( 
   const char* F, int L,  const ctrIdx &ica,
   const wbsparray<TD> &B, const ctrIdx &icb, wbsparray<TD> &C,
   const wbperm &pfinal, const TD& afac, const TD& cfac,
   char sflag, 
   TD eps
) const { 

   static int nvlog=0;

#ifdef WB_SPARSE_CLOCK
   Wb::Clock clk("sparse:ctr:all",0); 
#endif

   if (D.len<=1 && B.D.len<=1) { 
   if (contract_scalar(F_L,ica,B,icb,C,pfinal,afac,cfac)) {
      return C;
   }}

   if (!sflag && typeid(TD)!=typeid(double)) sflag='s';

   if (this==&C || &B==&C) { wbsparray<TD> X(C);
      contract(F,L,ica,B,icb,X,pfinal,afac,cfac,sflag,eps);
      return X.save2(C);
   }

   char iA=isDiag(FL), iB=B.isDiag(FL);
   if (iA || iB) {
      if (iA && iB) {
         if ((ica.len && ica.len<=2) && (icb.len && icb.len<=2))
            return contract_diag_diag(F_L,ica,B,icb,C,pfinal,afac,cfac);
         else {
            if (ica.len || icb.len) wblog(FL,
               "ERR %s() got %d/%d !?",FCT,ica.len,icb.len);
            return tensorProdX(FL,B,C,ica.conj,afac,icb.conj,cfac,&pfinal);
         }
      }
      else if (iA) {
          wbsparray<TD> A2; diag2reg(FL,A2);
          return A2.contract(F,L,ica,B,icb,C,pfinal,afac,cfac,sflag,eps);
      }
      else if (iB) {
          wbsparray<TD> B2; B.diag2reg(FL,B2);
          return contract(F,L,ica,B2,icb,C,pfinal,afac,cfac,sflag,eps);
      }
   }

   wbvector<SPIDX_T> sC; if (cfac) sC=C.SIZE; int e=0;

   if ((e=checkContract(F_L,SIZE,ica,B.SIZE,icb, &pfinal, &sC))) {
      if (e<0) { 
         if (sC.len) {
            C.init(sC); 
            return C;
         }
      }
      else wblog(FL,"ERR %s()",FCT);
   }

   if (!D.len || !B.D.len) {
      if (cfac) { C.D*=cfac; } else { C.init(sC); }
      if (e<0) { 
         if (C.SIZE.len || C.IDX.dim2) wblog(FL,"ERR %s() !?",FCT);
         if (!C.D.len) { C.D.init(1); C.D[0]=0; }
      }
      return C;
   }

   if (contract_check_2full(F_L,
         ica,B,icb,C,pfinal,afac,cfac,sflag,eps)) { return C; }

   wbsparray<TD> C0;

   if (cfac!=0) { C.D*=cfac; }
   if (afac==0) {
      if (!cfac) C.init(sC.Permute(pfinal));
      return C;
   }
   if (cfac!=0) C.save2(C0); 

   if (B.IDX.dim2==1) {
      if (B.SIZE.len!=1 || ica.len!=1) wblog(FL,
         "ERR %s() got S=[%s] !?",FCT,SSTR(B));
      wbvector<TD> b; B.toFull(b);
      contract(FL,ica[0],b,C); 

      C.Permute(pfinal); if (afac!=1) C*=afac;
      if (!C0.isEmpty()) { C+=C0; }

      return C;
   }

#ifdef WB_SPARSE_CLOCK
   Wb::Clock cl2("sparse:ctr:2",0); 
#endif

   sparseIndex2D<SPIDX_T> a2, b2;

     toIndex2D(ica,a2,2,'r'); 
   B.toIndex2D(icb,b2,1);     

   if (!a2.IJ.dim1 || a2.IJ.dim2!=2 || !a2.cIdx.len ||
       !b2.IJ.dim1 || b2.IJ.dim2!=2 || !b2.cIdx.len)
      wblog(FL,"ERR %s() unexpected IJ data (%s [%ld]; %s [%ld])",FCT,
      SSTR(a2.IJ), a2.cIdx.len,
      SSTR(b2.IJ), b2.cIdx.len
   );

   if (b2.rmaj) wblog(FL,
      "ERR %s() got unexpected sparse dimensions\n B: %s @ [%s]",
      FCT,SSTR(B), STR(icb)
   );

   int np=1; if (sp_num_threads>1) {
      np=(IDX.dim1+B.IDX.dim1)/1024; 
      if (np>sp_num_threads) { np=sp_num_threads; }
      if (unsigned(np)>b2.S[1]) { np=b2.S[1]; }
      if (np<1) { np=1; }
   }

  #ifdef LOAD_CGC_QSPACE
   if (np<sp_num_threads && CG_VERBOSE>5) {
      static time_t tlast=0;
      time_t tnow=time(NULL); 
      if (tnow>tlast+60) { size_t s1=(1<<26);
         if (numel()>s1 || B.numel()>s1 || sC.prod(0)>s1) {
            wblog(FL,"CTR %12s * %12s @ %d/%d threads",
               SSTR(a2), SSTR(b2), np, sp_num_threads
            );
            tlast=tnow;
         }
      }
   }
  #endif

   if (a2.rmaj) {
      size_t s1=a2.S[0], s2=b2.S[1], la=a2.IJ.dim1-1;
      const size_t nk=b2.cIdx.len-1;
      const SPIDX_T *pa=a2.P.data, *pb=b2.P.data;

      if (s1*s2>(1<<20)) wblog(FL,
         "ERR %s() got a.rmaj=%d having %dx%d * %dx%d !?!",FCT,
         a2.S[0], a2.S[1], b2.S[0], b2.S[1]
      );

      wbarray<TD> Cx(s1,s2);
      Wb::LogException ex; 

     #pragma omp parallel for num_threads(np)
      for (unsigned k=0; k<nk; ++k) { if (!ex) { try {
         size_t ib=b2.cidx(k), lb=b2.cidx(k+1);

      if (lb>ib) { --lb;
         size_t ia=0, ib_=ib;
         SPIDX_T *Ia=a2.IJ.data, *Ib=b2.IJ.data+2*ib;
         char c;

      while (ib<=lb) {
         c=NUMCMP(Ia[1],Ib[0]);

         if (c==0) {
            if (ia<la) {
               if (Ia[0]!=Ia[2]) c=-2; 
            }
            else {
               c=+3;
            }

            if (ib<lb) {
               if (!c && Ib[1]!=Ib[3]) c=+2; 
            }

            if (Ia[0]>=s1 || Ib[1]>=s2) wblog(FL,
               "ERR %s() index out of bounds (%ld/%ld; %ld/%ld)",
               FCT, Ia[0], s1, Ib[1], s2);

            Cx(Ia[0],Ib[1]) += D[pa[ia]] * B.D[pb[ib]];
         }

         if (ia<la && (c<0 || ib==lb)) {
            if (Ia[0]==Ia[2]) { if (Ia[1]>=Ia[3]) wblog(FL,
               "ERR %s() a.IJ not sorted or unique !?",FCT); }
            else {
               if (Ia[0]>Ia[2]) wblog(FL,
                  "ERR %s() a.IJ not sorted !?",FCT);
               Ib-=2*(ib-ib_); ib=ib_; 
            }
            ++ia; Ia+=2;
         }
         else if (c) { 
            if (ib<lb) {
               if (Ib[1]==Ib[3]) { if (Ib[0]>=Ib[2]) wblog(FL,
                  "ERR %s() b.IJ not sorted or unique !?",FCT); }
               else {
                  if (Ib[1]>Ib[3]) wblog(FL,
                     "ERR %s() b.IJ not sorted !?",FCT);
                  ib_=ib+1; Ia-=2*ia; ia=0; 
               }
            }
            ++ib; Ib+=2;
         }
         else {
            if (ia<la && Ia[1]>=Ia[3]) 
               wblog(FL,"ERR %s() a.IJ not sorted !?",FCT);
            if (ib<lb && Ib[0]>=Ib[2]) 
               wblog(FL,"ERR %s() b.IJ not sorted !?",FCT);

            ++ia; Ia+=2;
            ++ib; Ib+=2;
         }
      }}}
         catch (Wb::LogException &e_) { ex+=e_; }
         catch (...) { ++ex; }
      }}
      ex.report(FLF);

      C.init(FL,Cx).Reshape(sC);
      if (!C.SIZE.len && !C.D.len) {
         C.D.init(1); C.D[0]=0;
      }

      C.Permute(pfinal); if (afac!=1) C*=afac;
      if (!C0.isEmpty()) { C+=C0; }

      return C;
   }

   wbMatrix<SPIDX_T> MARK(np,a2.S[0]);
   groupIndex<widx_t> kc; { kc.init(b2.S[1]); }

   const SPIDX_T *pa=a2.P.data, *pb=b2.P.data;
   SPIDX_T k=0, ltot=0, nnzc=0;

  #pragma omp parallel for num_threads(np) reduction(+:ltot)
   for (k=0; k<b2.S[1]; ++k) { 
      int id=omp_get_thread_num();
      widx_t *mark=MARK.ref(id,0);
      SPIDX_T i, i_,i1,i2, j, j_, j1=b2.cidx(k), j2=b2.cidx(k+1), nnzc=0;
      for (j_=j1; j_<j2; ++j_) { j=b2.ridx(j_);
         i1=a2.cidx(j); i2=a2.cidx(j+1);
         for (i_=i1; i_<i2; ++i_) { i=a2.ridx(i_);
            if (mark[i]<k+1) { mark[i]=k+1; ++nnzc; }
         }
      }
      kc.setD(k,nnzc);
      ltot+=nnzc;
   }

   nnzc=kc.init_ic();
   if (nnzc!=ltot) wblog(FL,"ERR %s() "
      "nnz inconsistency (%ld/%ld) !?",FCT,ltot,nnzc);

#ifdef WB_SPARSE_CLOCK
   cl2.Switch("sparse::ctr:3"); 
#endif

   if (!nnzc) {
      if (C0.isEmpty()) C.init(sC); else C0.save2(C);
      return C;
   }

   wbMatrix<TD> W(np,a2.S[0]);
   wbMatrix<SPIDX_T> IC_(np,a2.S[0]), IC(np,a2.S[0]);
   MARK.set(0);

   ctrIdx ika, ikb;
   int largeD=0;

   ica.invert(  SIZE.len,ika);
   icb.invert(B.SIZE.len,ikb);

   {  size_t DX=(1<<26); 

      if (D.len>DX || B.D.len>DX || nnzc>DX) { largeD=2; nvlog=1; }
      else {
         DX=(1<<20); 
         if (D.len>DX || B.D.len>DX || nnzc>DX) { largeD=1;
            if ((++nvlog)%24==1) largeD|=2;
         }
      }

      if (largeD&2) {
         wblog(FL,"TST %s() large sparse arrays (%d/%d threads)\n"
            "  A: %-20s (%9.3g) @ %s\n  B: %-20s (%9.3g) @ %s\n"
            "> C: %-20s (%9.3g)",FCT, np, sp_num_threads,
         SSTR_(this), D.len/double(  numel()), STR(ica),
         SSTR(B),   B.D.len/double(B.numel()), STR(icb),
         STR(sC),      nnzc/double(sC.prod(1)));
      }
      else if (largeD) {
        wblog(FL,"CTR %s_%s * %s_%s = %s (np=%d/%d)", SSTR_(this),
        STR(ica), SSTR(B), STR(icb), SSTR(sC), np, sp_num_threads);
      }
      if (largeD) { Wb::MemStat(FL); doflush(); }
   }

   C.init(sC,nnzc);
   ltot=0; 

  #pragma omp parallel for num_threads(np) reduction(+:ltot)
   for (k=0; k<b2.S[1]; ++k) {
      int id=omp_get_thread_num();
      wbvector<SPIDX_T> ic(IC.dim2,IC.ref(id,0),'r');
      SPIDX_T *ic_=IC_.ref(id,0), is=0;
      widx_t *mark=MARK.ref(id,0);
      TD *w=W.ref(id,0);

      SPIDX_T i, i_,i1,i2, j, j_, j1=b2.cidx(k), j2=b2.cidx(k+1);
      bool sortflag=((j2-j1)<SPIDX_T(ceil(0.2*MARK.dim2)));

      for (j_=j1; j_<j2; ++j_) { j=b2.ridx(j_);
         const TD &Bjk=B.D[pb[j_]];
         i1=a2.cidx(j); i2=a2.cidx(j+1);
         for (i_=i1; i_<i2; ++i_) { i=a2.ridx(i_);
            if (mark[i]<k+1) { mark[i]=k+1; ic_[i]=pa[i_];
               if (sortflag) ic[is++]=i;
               w[i] = D[pa[i_]] * Bjk; 
            }
            else {
               w[i]+=(D[pa[i_]] * Bjk);
            }
         }
      }

      SPIDX_T l=kc.getGroup(k);

      if (sortflag) {
         if (is) {
            ic.len=is; ic.Sort();
            for (i_=0; i_<ic.len; ++i_) { i=ic[i_];
               C.IDX.recSetP(l,IDX.rec(ic_[i]),ika,B.IDX.rec(pb[j1]),ikb);
               C.D[l++]=w[i];
            }
            ic.len=MARK.dim2;
         }
         ltot+=is;
      }
      else {
         for (i=0; i<MARK.dim2; ++i) {
            if (mark[i]==k+1) { 
               C.IDX.recSetP(l,IDX.rec(ic_[i]),ika,B.IDX.rec(pb[j1]),ikb);
               C.D[l++]=w[i]; ++ltot;
            }
         }
      }
   }  

   if (ltot!=C.D.len) wblog(FL,
      "ERR %s() %d/%d/%d !?",FCT,ltot,C.D.len,nnzc);

   if (largeD) {
      size_t s=C.memSize();
      if (s>(1<<30)) wblog(FL,
         "TST %s() size(C) = %s",FCT,Wb::size2Str(s).data);
      Wb::MemStat(FL); 
   }
#ifdef WB_SPARSE_CLOCK
   cl2.stop();
#endif

   C.Permute(pfinal); if (afac!=1) C*=afac;
   if (!C0.isEmpty()) { C+=C0; }

   return C;
};

template <class TD>
wbsparray<TD>& wbsparray<TD>::contract( 
   const char *F, int L,
   unsigned k, const wbvector<TD> &B, wbsparray<TD> &C
) const {

   if (k>=SIZE.len) {
      if (B.len!=1) wblog(F_L,"ERR %s() "
         "size mismatch (%s ; %d @ %d)",FCT,SSTR_(this),B.len,k+1);
      C=*this; C*=B.data[0];
      if (ABS(B.data[0])<TD(1E-8)) C.Compress();
      return C;
   }

   if (SIZE[k]!=B.len) wblog(F_L,
      "ERR %s() size mismatch (%s ; %d @ %d)",FCT,SSTR_(this),B.len,k+1);
   if (IDX.dim1!=D.len || IDX.dim2!=SIZE.len) wblog(FL,
      "ERR %s() size inconsistency (%dx%d <> %dx%d) !?",
      FCT,IDX.dim1,IDX.dim2,D.len,SIZE.len);

   if (IDX.dim2>1) {
      const SPIDX_T *idx=IDX.data+k; 
      const TD *b=B.data; TD *c; wbindex Ik;

      Ik.Index_ex(IDX.dim2,k); C.init();
      IDX.getCols(Ik,C.IDX); C.D=D; c=C.D.data;
      SIZE.select(Ik,C.SIZE);

      for (SPIDX_T n=IDX.dim2, i=0; i<D.len; ++i, idx+=n) {
         if ((*idx)>=B.len) wblog(F_L,
            "ERR %s() index out of bounds (%d/%d)",FCT,(*idx)+1,B.len);
         c[i]*=b[*idx];
      }

      C.Compress();
   }
   else if (IDX.dim2==1) { 
      if (int(k)>=0 && k) wblog(FL,
         "ERR %s() ctrIdx out of bounds (%d/1)",FCT,k+1);
      C.initScalar();

      const SPIDX_T *idx=IDX.data;
      const TD *a=D.data, *b=B.data;
      TD &c=C.D[0]; c=0;

      for (SPIDX_T i=0; i<D.len; ++i) {
         if (idx[i]>=B.len) wblog(F_L,
            "ERR %s() index out of bounds (%d/%d)",FCT,idx[i]+1,B.len);
         c+=a[i]*b[idx[i]];
      }
   }
   else {
      info(FL); wblog(FL,"ERR %s()",FCT);
   }

   return C;
};

template <class TD>
wbarray<TD>& wbsparray<TD>::contract(
   const char* F, int L,  const ctrIdx &ica,
   const wbsparray<TD> &B, const ctrIdx &icb, wbarray<TD> &C,
   const wbperm &pfinal, const TD& afac, const TD& cfac
) const {

   wbarray<TD> Af,Bf,X;

   if (cfac && !C.isEmpty()) { C.save2(X); X*=cfac; }

   this->toFull(Af); B.toFull(Bf);
   Af.contract(FL,ica,Bf,icb,C,pfinal,afac);

   if (!X.isEmpty()) {
      if (!X.sameSize(C)) wblog(FL,
         "ERR %s() severe size mismatch (%s <> %s; %g)",
         FCT, SSTR(C), SSTR(X), cfac);
      C+=X;
   }
   return C;
};

template <class TD>
wbarray<TD>& contract(const char* F, int L,
   const wbarray<TD>   &A, const wbindex &ic1,
   const wbsparray<TD> &B, const wbindex &ic2, wbarray<TD> &C,
   const wbperm &pfinal, const TD& afac, const TD& cfac
){
   unsigned m=ic1.len, r1=A.rank(), r2=B.rank(FL); int e=0;
   wbarray<TD> X;

   if ((void*)(&A)==(void*)(&C)) wblog(F,L,
      "ERR sparse::%s() must have distinct target space",FCT);

   wbvector<SPIDX_T> sC; if (cfac) sC=C.SIZE;
   if ((e=checkContract(F_L,A.SIZE,ic1,B.SIZE,ic2, &pfinal, &sC))) {
      if (e<0) { C.init(sC); if (sC.len) return C; }
      else wblog(FL,"ERR %s()",FCT);
   }

   if (afac==0 || !B.D.len) {
      if (cfac!=0) { C*=cfac; } else { C.init(sC.Permute(pfinal)); }
      return C;
   }

   if (cfac!=0) { C*=cfac; C.save2(X); }
   C.init(sC);

   wbvector<SPIDX_T> s1,D2;
   wbMatrix<SPIDX_T> I2;
   wbperm Pa,P2;

   A.SIZE.getI(ic1,s1); Pa.init2End(ic1,r1);
   B.IDX.cols2End(ic2,I2).groupRecs(P2,D2,r2-m);

   wbindex Ia(r1),Ic(sC.len);
   SPIDX_T i, ig,d,k,l, k1=r1-ic1.len, k2=r2-ic2.len;
   SPIDX_T *i2, *ia=Ia.data, *ic=Ic.data, *pa=Pa.data;
   wbIndex K1(s1);

   for (l=ig=0; ig<D2.len; ++ig) { d=D2[ig];
   for (k=0; k<d; ++k,++l) {
       TD b=B.D.data[P2[l]]; if (b==0) continue; if (afac!=1) b*=afac;

       i2=I2.ref(l); if (k==0) {
       for (i=0; i<k2; ++i) ic[k1+i]=i2[i]; }     
       for (i=0; i<m; ++i) ia[pa[k1+i]]=i2[k2+i]; 

       for (K1.reset(); ++K1;) {
          for (i=0; i<k1; ++i) { ia[pa[i]]=ic[i]=K1.data[i]; }
          if (b== 1) { C(Ic)+=A(Ia); } else 
          if (b==-1) { C(Ic)-=A(Ia); } else { C(Ic)+=b*A(Ia); }
       }
   }}

   C.Permute(pfinal); if (!X.isEmpty()) C+=X;
   return C;
};

template <class TD>
wbarray<TD>& contract(const char* F, int L,
   const wbsparray<TD> &A, const wbindex &ic1,
   const wbarray<TD>   &B, const wbindex &ic2, wbarray<TD> &C,
   const wbperm &pfinal, const TD& afac, const TD& cfac
){
   unsigned m=ic1.len, r1=A.rank(FL), r2=B.rank(); int e=0;
   wbarray<TD> X;

   if ((void*)(&B)==(void*)(&C)) wblog(F,L,
      "ERR sparse::%s() must have distinct target space",FCT);

   wbvector<SPIDX_T> sC; if (cfac) sC=C.SIZE;
   if ((e=checkContract(F_L,A.SIZE,ic1,B.SIZE,ic2, &pfinal, &sC))) {
      if (e<0) { C.init(sC); if (sC.len) return C; }
      else wblog(FL,"ERR %s()",FCT);
   }

   if (afac==0 || !A.D.len) {
      if (cfac!=0) { C*=cfac; } else { C.init(sC.Permute(pfinal)); }
      return C;
   }

   if (cfac!=0) { C*=cfac; C.save2(X); }
   C.init(sC);

   wbvector<SPIDX_T> s2,D1;
   wbMatrix<SPIDX_T> I1;
   wbperm Pb,P1;

   A.IDX.cols2End(ic1,I1).groupRecs(P1,D1,r1-m);
   B.SIZE.getI(ic2,s2); Pb.init2End(ic2,r2);

   wbindex Ib(r2),Ic(sC.len);
   SPIDX_T i, ig,d,k,l, k1=r1-ic1.len, k2=r2-ic2.len, *i1,
      *ib=Ib.data, *ic=Ic.data, *pb=Pb.data;
   wbIndex K2(s2);

   for (l=ig=0; ig<D1.len; ++ig) { d=D1[ig];
   for (k=0; k<d; ++k,++l) {
       TD a=A.D.data[P1[l]]; if (a==0) continue; if (afac!=1) a*=afac;

       i1=I1.ref(l); if (k==0) {
       for (i=0; i<k1; ++i) ic[i]=i1[i]; }        
       for (i=0; i<m; ++i) ib[pb[k2+i]]=i1[k1+i]; 

       for (K2.reset(); ++K2;) {
          for (i=0; i<k2; ++i) { ib[pb[i]]=ic[k1+i]=K2.data[i]; }
          if (a== 1) { C(Ic)+=B(Ib); } else 
          if (a==-1) { C(Ic)-=B(Ib); } else { C(Ic)+=a*B(Ib); }
       }
   }}

   C.Permute(pfinal); if (!X.isEmpty()) C+=X;
   return C;
};

template <class TD>
TD Wb::VMatVprod(
  const wbsparray<TD> &v1, const wbsparray<TD> &M, const wbsparray<TD> &v2
){
   wbsparray<TD> x,X; wbindex i1(1), i2(1);
   i1[0]=0; i2[0]=1; 

   M.contract(FL,i2,v2,i1,X);
   v1.contract(FL,i1,X,i1,x);

   x.checkSize(FL);
   if (x.D.len && !x.isScalar()) {
      MXPut(FL,"q").add(M,"M").add(v1,"v1").add(v2,"v2")
        .add(X,"X").add(x,"x");
      wblog(FL,"ERR %s() got %s",FCT,SSTR(x));
   }

   return (x.D.len ? x[0] : 0.);
};

template <class TD>
wbsparray<TD>& Wb::MatProd(
  const wbsparray<TD> &A, const wbsparray<TD> &B, wbsparray<TD> &C,
  char aflag0, char bflag0, TD afac, TD cfac
){
   if (&C==&A || &C==&B) {
      wbsparray<TD> X; if (cfac) X=C;
      Wb::MatProd(A,B,X,aflag0,bflag0,afac,cfac);
      return X.save2(C);
   }

   int ra=A.rank(), rb=B.rank(), e=0;
   if (ra!=2) { if (ra || !A.isEmpty()) e+=1; }
   if (rb!=2) { if (rb || !B.isEmpty()) e+=2; }

   if (e) wblog(FL,
      "ERR %s() expects matrices (%s; %s)",FCT,SSTR(A),SSTR(B));
   if (!ra || !rb) { 
      if (!C.isEmpty() && cfac!=1) C*=cfac;
      return C;
   }

   opFlags<TD> aflag(aflag0), bflag(bflag0);
   const wbsparray<TD>
      &AX = aflag.applyConjOrRef(FL,A),
      &BX = bflag.applyConjOrRef(FL,B);
   ctrIdx ia(1),ib(1);

   ia[0]=(aflag.trans() ? 0:1);
   ib[0]=(bflag.trans() ? 1:0);

   AX.contract(FL,ia,BX,ib,C, wbperm(), afac, cfac);

   if (&AX!=&A) { delete &AX; }
   if (&BX!=&B) { delete &BX; }; return C;
};

template <class TD>
wbsparray<TD>& wbsparray<TD>::comm(const char *F, int L,
  const wbsparray<TD> &B, wbsparray<TD> &C, char aflag, char bflag
) const {

   if (&C==this || &C==&B) {
      wbsparray<TD> X; comm(F,L,B,X,aflag,bflag);
      return X.save2(C);
   }

   SPIDX_T d=-1;
   if (!isSMatrix(FL,&d) || !B.isSMatrix(FL,&d)) wblog(F_L,
      "ERR %s() invalid operators for [A,B] (%s; %s)", FCT,
      SSTR_(this), SSTR(B)
   );

   Wb::MatProd(*this,B,C,aflag,bflag);                
   Wb::MatProd(B,*this,C,bflag,aflag,TD(-1.),TD(1.)); 

   return C;
};

#endif

