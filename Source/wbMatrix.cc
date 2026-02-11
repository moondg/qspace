/* ---------------------------------------------------------------------
 * Project : QSpace tensor library (v4.0)
 * Class   : wbMatrix (matrix class, row-major)
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

#ifndef __WB_MATRIX_ROW_MAJOR_CC__
#define __WB_MATRIX_ROW_MAJOR_CC__

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
// quick integer data type // Wb,Dec02,11

template<class T> inline
bool wbMatrix<T>::thisIsInt() const { return 0; };

template <> inline
bool wbMatrix<int>::thisIsInt() const { return 1; };

template<> inline
bool wbMatrix<unsigned>::thisIsInt() const { return 1; };

template<> inline
bool wbMatrix<size_t>::thisIsInt() const { return 1; };

template <> inline
bool wbMatrix<char>::thisIsInt() const { return 1; };

template <> inline
bool wbMatrix<unsigned char>::thisIsInt() const { return 1; };

template <>
bool wbMatrix<double>::isNormal() const {
   for (size_t n=dim1*dim2, i=0; i<n; ++i)
       if (data[i] && !isnormal(data[i])) return 0;
   return 1;
};

template <>
bool wbMatrix<float>::isNormal() const {
   for (size_t n=dim1*dim2, i=0; i<n; ++i)
       if (data[i] && !isnormal(data[i])) return 0;
   return 1;
};

template <>
bool wbMatrix<wbcomplex>::isNormal() const { 
   for (size_t n=dim1*dim2, i=0; i<n; ++i) {
       if (data[i].r && !isnormal(data[i].r)) return 0;
       if (data[i].i && !isnormal(data[i].i)) return 0;
   }
   return 1;
};

template <> bool wbMatrix<unsigned>::isNormal() const { return 1; };
template <> bool wbMatrix<size_t  >::isNormal() const { return 1; };
template <> bool wbMatrix<int     >::isNormal() const { return 1; };
template <> bool wbMatrix<char    >::isNormal() const { return 1; };
template <> bool wbMatrix<long    >::isNormal() const { return 1; };

template <class T> inline
wbMatrix<T>& wbMatrix<T>::init( 
   const char *F, int L, const mxArray *a,
   char tflag, char ref, char tcheck
){
   if (!a || mxIsEmpty(a)) { init(); return *this; }

   unsigned r=mxGetNumberOfDimensions(a);
   const size_t *S=mxGetDimensions(a);

   if (r!=2 || !mxIsNumChar(a)) wblog(F,L,
      "ERR got rank-%d array `%s'",r,mxGetClassName(a));

   size_t d1=(tflag ? S[1]:S[0]), d2=(tflag ? S[0]:S[1]);

   Mx::Array<T> A(F_L,a,!tflag);
   if (!ref || !A.data) { if (ref) wblog(FL,"WRN invalid ref: %s",STR(A));
      RENEW(d1,d2);
      A.copyTo(data,tcheck);
   }
   else { init2ref(d1,d2,A.data); } 

   return *this;
};

template <class T> inline
wbMatrix<T>& wbMatrix<T>::init2ref(size_t r, size_t c, const T* d) {
   if (data && !isref) { init(); }; isref=1;
   dim1=r; dim2=c; data=(T*)d;
   return *this;
};

template <class T> inline
wbMatrix<T>& wbMatrix<T>::init2ref(const wbvector<T> &v, const char tflag) {
   if (data && !isref) { init(); }; isref=1; data=v.data;
   if (tflag)
        { dim1=v.len; dim2=1; }
   else { dim1=1; dim2=v.len; }
   return *this;
};

template <class T> inline
wbMatrix<T>& wbMatrix<T>::init2ref(const wbMatrix<T> &A) {
   if (data && !isref) { init(); }; isref=1;
   dim1=A.dim1; dim2=A.dim2; data=A.data;
   return *this;
};

template <class T>
template<class TA, class TB>
wbMatrix<T>& wbMatrix<T>::init(
   const char *F, int L,
   const wbMatrix<TA> &M, const wbvector<TB> &V, char dim
){
   size_t i,j,n1=M.dim1,n2=M.dim2;
   if (abs(dim)==2) {
      if (n1!=V.len) wblog(F_L,
         "ERR %s() size mismatch %d/%dx%d",FCT,V.len,n1,n2);
      RENEW(n1,n2+1,0,0); 
      if (dim>0) {
         for (i=0; i<n1; ++i) {
            for (j=0; j<n2; ++j) (*this)(i,j)=(T)M(i,j);
            (*this)(i,j)=(T)V[i];
         }
      }
      else {
         for (i=0; i<n1; ++i) {
            (*this)(i,0)=(T)V[i];
            for (j=0; j<n2; ++j) (*this)(i,j+1)=(T)M(i,j);
         }
      }
   }
   else if (abs(dim)==1) {
      if (n2!=V.len) wblog(F_L,
         "ERR %s() size mismatch %dx%d/%d",FCT,n1,n2,V.len);
      RENEW(n1+1,n2,0,0); 
      if (dim>0) {
         for (i=0; i<n1; ++i)
         for (j=0; j<n2; ++j) (*this)(i+1,j)=(T)M(i,j);
         for (j=0; j<n2; ++j) (*this)(i,j)=(T)V[j];
      }
      else {
         for (j=0; j<n2; ++j) (*this)(0,j)=(T)V[j];
         for (i=0; i<n1; ++i)
         for (j=0; j<n2; ++j) (*this)(i+1,j)=(T)M(i,j);
      }
   }
   else wblog(F_L,"ERR %s() invalid dim=%d",FCT,dim);

   return *this;
};

template <class T> inline
int matchSortedIdxU( 
    const char *F, int L,
    cTMAT &QA, cTMAT &QB, wbindex &Ia, wbindex &Ib,
    widx_t m,  
     char lex   
){
    widx_t m1=0, m2=0; int i=0;
    if (QA.dim2!=QB.dim2) wblog(FL,"ERR %s() dimension mismatch "
       "(%dx%d <> %d/%d)",FCT,QA.dim1,QA.dim2,QB.dim1,QB.dim2);

    try {
       i=Wb::matchSortedIdx(
          QA.data, QA.dim2, QA.dim1,
          QB.data, QB.dim2, QB.dim1, Ia,Ib, m, lex, &m1, &m2
       );
       if (m1 || m2) wblog(F_L,"ERR %s() "
          "got non-unique/sorted input records (%d,%d)",FCT,m1,m2);
    }
    catch (...) {
#ifdef MATLAB_MEX_FILE
       MXPut(FL,"q").add(QA,"A").add(QB,"B").add(int(m),"m").add(lex,"lex");
       wblog(FL,"ERR %s()",FCT);
#endif
    }

    return i;
};

template<class T> inline 
bool wbMatrix<T>::quickCheckSorted(
   char dir,             
   char lex,             
   const char *F, int L, 
   size_t n              
) const {

   char c; dir=(dir>0 ? +1 : -1);

   for (size_t i=1; i<dim1; i++) {
      c=recCompare(i,i-1,-1,lex);
      if (c) {
         if (c!=dir) {
            if (F) wblog(F,L,
              "ERR %s() input records not sorted %s (%d)",
               FCT, dir>0 ? "ascendingly":"descendingly", i);
            else return 0;
         }
         if (int(--n)<=0) break;
      }
   }

   return 1;
};

template<> 
WBINDEX& wbMatrix<size_t>::toIndex(
   WBINDEX &I, const WBINDEX *S) const {

   const size_t *idx=data;
   if (S==nullptr) {
      WBINDEX SX(dim2); widx_t *s=SX.data; size_t i=0,j;

      for (; i<dim1; ++i, idx+=dim2) {
         for (j=0; j<dim2; ++j) {
            if (s[j]<idx[j]) s[j]=idx[j];
         }
      }
      SX+=1; return toIndex(I,&SX);
   }

   I.init(dim1); if (!dim1) return I;
   if (!dim2) wblog(FL,
      "ERR %s() got empty index matrix (%dx%d)",FCT,dim1,dim2);

   size_t i=0,j, l=dim2-1;
   const widx_t *s=S->data;

   for (; i<dim1; ++i, idx+=dim2) { widx_t &q=I.data[i];
      for (j=l, q=idx[j--]; j<l; --j) { q = q*size_t(s[j]) + size_t(idx[j]); }
   }

   return I;
};

template<> 
WBIDXMAT& wbMatrix<size_t>::toIndex2D(
   const WBINDEX &S, const wbindex &Ic,
   WBIDXMAT &IJ, char pos 
 ) const {

   wbvector<char> M(S.len); 
   size_t i; char *m=M.data;

   if (S.len!=dim2) wblog(FL,
      "ERR %s() size mismatch (%dx%d/%d)",FCT,dim1,dim2,S.len);
   if (!dim2) wblog(FL,
      "ERR %s() got empty index matrix (%dx%d)",FCT,dim1,dim2);

   for (i=0; i<Ic.len; ++i) {
      if (Ic[i]>=S.len) wblog(FL,"ERR %s() index out of bounds "
         "(%s; %d)",FCT,STR(Ic+1),S.len);
      if ((++m[Ic[i]])>1) wblog(FL,"ERR %s() index not unique "
         "(%s; %d)",FCT,STR(Ic+1),S.len
      );
   }

   if ((--pos)>1) wblog(FL, 
      "ERR %s() index out of founds (pos=%d)",FCT,pos);

   if (dim2==2 && Ic.len==1) {
      if (int(Ic[0])==pos) { IJ=(*this); }
      else {
         IJ.init(dim1,dim2);
         size_t d12=dim1*dim2; widx_t *d=IJ.data;
         for (i=0; i<d12; i+=2) {
            d[i]=data[i+1]; 
            d[i+1]=data[i];
         }
      }
      return IJ;
   }

   IJ.init(dim1,2); if (!dim1) return IJ;

   wbindex Ik(S.len-Ic.len);
   unsigned j=0, k=0, lk=Ik.len-1, lc=Ic.len-1;
   const size_t *idx=data, *s=S.data;
   const widx_t *ic=Ic.data; widx_t *ik=Ik.data;

   for (i=0; i<S.len; ++i) { if (!m[i]) { ik[j++]=i; }}

   for (i=0; i<dim1; ++i, idx+=dim2, k+=2) {
      if (Ik.len) {
         widx_t &I=IJ.data[k+(1-pos)]; I=idx[ik[lk]]; 
         for (j=lk-1; j<lk; --j) { if (I) I*=s[ik[j]]; I+=idx[ik[j]]; }
      }
      if (Ic.len) {
         widx_t &J=IJ.data[k+pos]; J=idx[ic[lc]]; 
         for (j=lc-1; j<lc; --j) { if (J) J*=s[ic[j]]; J+=idx[ic[j]]; }
      }
   }

   return IJ;
};

template<> 
wbMatrix<double>& wbMatrix<double>::SkipTiny_float(double x) {
   Wb::chopTiny_float(data,dim1*dim2,x);
   return *this;
}

template<> 
wbMatrix<wbcomplex>& wbMatrix<wbcomplex>::SkipTiny_float(wbcomplex x) {
   Wb::chopTiny_float((double*)data,2*dim1*dim2,x.r);
   return *this;
}

template <class T> inline
size_t wbMatrix<T>::maxRec(char lex, size_t *m_) const {

   if (!dim1 || !dim2) { wblog(FL,
      "WRN %s() got empty matrix (%d,%d)",FCT,dim1,dim2);
      return 0;
   }

   size_t i=1, k=0, m=1; char c;
   for (; i<dim1; i++) { 
      c=recCompare(k,i,-1,lex);
      if (c<0) { k=i; m=1; }
      else if (c==0) { m++; }
   }

   if (m_) (*m_)=m;
   return k;
};

template <> inline
size_t wbMatrix<double>::maxRec_float(
   char lex, size_t *m_, double xref
 ) const {

   wbMatrix<double> M(*this);
   Wb::chopTiny_float(M.data,dim1*dim2,xref);
   return M.maxRec(lex,m_);
};

template <class T> inline
T wbMatrix<T>::recMax(size_t r, size_t *k) const {

   if (r>=dim1) wblog(FL,
      "ERR %s() index out of bounds (%d/%d)",FCT,r,dim1);
   if (!dim2) { wblog(FL,
      "WRN %s() from empty matrix (%dx%d)",FCT,dim1,dim2);
      return 0;
   }

   T *d=data+r*dim2, x=d[0]; if (k) (*k)=0; 

   for (size_t i=1; i<dim2; ++i) {
      if (x<d[i]) { x=d[i]; if (k) (*k)=i; }
   }

   return x;
};

template <class T> inline
T wbMatrix<T>::colMax(size_t c, size_t *k) const {

   if (c>=dim2) wblog(FL,
      "ERR %s() index out of bounds (%d/%d)",FCT,c,dim2);
   if (!dim1) { if (k) (*k)=0; wblog(FL,
      "WRN %s() from empty matrix (%dx%d)",FCT,dim1,dim2);
       return 0;
   }

   T *d=data+c, x=d[0]; d+=dim2; if (k) (*k)=0;

   for (size_t i=1; i<dim1; ++i, d+=dim2) { 
      if (x<d[0]) { x=d[0]; if (k) (*k)=i; }
   }

   return x;
};

template <class T> inline
double wbMatrix<T>::recMaxA(size_t r, size_t *k) const {
   if (k) (*k)=0;

   if (r>=dim1) wblog(FL,
      "ERR %s() index out of bounds (%d/%d)",FCT,r,dim1);
   if (!dim2) { wblog(FL,
      "WRN %s() from empty matrix (%dx%d)",FCT,dim1,dim2);
       return 0;
   }

   T *d=data+r*dim2; 
   double a=ABS(d[0]), amax=a;

   for (size_t i=1; i<dim2; i++) {
      a=ABS(d[i]); if (amax<a) { amax=a; if (k) (*k)=i; }
   }

   return amax;
};

template <class T>
double wbMatrix<T>::maxRelDiff(const wbMatrix<T> &M) const {

   size_t i, s=dim1*dim2;
   double maxdiff=0., maxval=0.;

   if (!sameSize(M)) wblog(FL,
   "ERR %s() incompatible matrix objects",FCT);

   for (i=0; i<s; i++) {
       maxdiff=MAX(maxdiff, fabs(double(M.data[i]-data[i])));
       maxval=MAX(maxval,
              MAX( fabs(double(data[i])), fabs(double(M.data[i])) ) );
   }

   return ((maxval==0.) ? 0. : maxdiff/maxval);
};

template <class T>
T wbMatrix<T>::normDiff2(const wbMatrix<T> &M, size_t *k) const {

   if (!sameSize(M)) wblog(FL,
      "ERR %s() incompatible matrix objects (%dx%d, %dx%d)",
       FCT,dim1,dim2,M.dim1,M.dim2);
   T fac=1;

   return Wb::rangeNormDiff2(M.data, data, dim1*dim2, fac, k);
};

template <class T>
template <class T2>
double wbMatrix<T>::normDiff(const wbMatrix<T2> &M) const {

   if (!sameSize(M)) wblog(FL,
      "ERR %s() incompatible matrix objects (%dx%d, %dx%d)",
       FCT,dim1,dim2,M.dim1,M.dim2);

   double x, x2=0;
   for (size_t n=numel(), i=0; i<n; ++i) {
      x=(data[i]-M.data[i]); x2+=(x*x);
   }

   return sqrt(x2);
};

template<class T>
bool wbMatrix<T>::isOrthoRows(T *x, char tnorm, T eps) const {

   if (!dim1 || !dim2 || dim1>dim2) return 0;

   wbarray<T> E, A(*this,'r');
   Wb::MatProd(A,A,E, tnorm? 'T':'C');

   if (x)
        return E.isProptoId(*x,eps); 
   else return E.isDiagMatrix();     
};

template<class T>
bool wbMatrix<T>::isOrthoRows( 
   const wbMatrix<T> &B_, T *x, char tnorm, T eps) const {

   if (dim1>dim2 || B_.dim1>B_.dim2) return 0;
   if (!(*this) && !B_) return 0;

   wbarray<T> E, A(*this,'r'), B(B_,'r');
   Wb::MatProd(A,B,E, tnorm? 'T':'C');

   if (x)
        return E.isProptoId(*x,eps); 
   else return E.isDiagMatrix();     
};

template <class T> inline
wbMatrix<T>& wbMatrix<T>::NormalizeCol(size_t k) {

   if (k>=dim2) wblog(FL,"ERR %s() index out of bounds (%d/%d)",FCT,k,dim2);
   if (!dim1) return *this;

   T *d=data+k, x2 = Wb::overlap(d,d,dim1,dim2);
   Wb::timesRange(d,1/Wb::sqrt(x2),dim1,dim2);

   return *this;
};

template <class T> inline
wbMatrix<T>& wbMatrix<T>::set2Cols(const size_t i1, const size_t i2){
   const wbMatrix<T> X; save2(X);
   X.getCols(i1,i2,*this);
   return *this;
};

template <class T> inline
wbMatrix<T>& wbMatrix<T>::getCols(
   size_t j1, size_t j2, 
   wbMatrix<T> &M 
 ) const {

   if (&M==this) {
      wbMatrix<T> X; getCols(j1,j2,X);
      return X.save2(M);
   }

   if (j1==0 && j2+1==dim2) return M.init(*this);

   if (j1>=dim2 || j2>=dim2) {
      if (int(j1-j2)==1 && (j1==0 || j1==dim2)) {
         return M.init(dim1,0);
      }
      else wblog(FL,
     "ERR index out of bounds (%dx%d: %d,%d)",dim1,dim2,j1,j2);
   }
   if (j1>j2) return M.init(dim1,0);

   size_t i=0, m=j2-j1+1; const T *p=data+j1;

   M.init(dim1,m);

   for (; i<dim1; ++i)
   MEM_CPY<T>(M.data+i*m, m, p+i*dim2); 

   return M;
};

template <class T> inline
wbvector<T>& wbMatrix<T>::getCol(
   const size_t k, wbvector<T> &v
 ) const {

   size_t i=0; const T *p=data+k;

   if (k>=dim2) wblog(FL,
      "ERR %s() index out of bounds (%dx%d/%d)",FCT,dim1,dim2,k);
   v.init(dim1); 
   for (; i<dim1; ++i, p+=dim2) { v[i]=p[0]; }

   return v;
};

template <class T>
template <class TI> inline
wbMatrix<T>& wbMatrix<T>::getCols(
   size_t n, const TI *I, wbMatrix<T> &M) const {

   wbMatrix<T> X(dim1,n); T *d0=data, *d=X.data;
   size_t i=0, j;

   for (j=0; j<n; ++j) if (I[j]>=dim2) wblog(FL,
      "ERR %s() index out of bounds (%d: %d/%d)",FCT,j,I[j],dim2);

   for (; i<dim1; ++i, d+=n, d0+=dim2) { 
   for (j=0; j<n; ++j) d[j]=d0[I[j]]; }

   X.save2(M); return M;
};

template <class T> inline
wbvector<T>& wbMatrix<T>::getRec(
   const size_t j0, wbvector<T> &v, 
   char ref
 ) const {

   if (j0>=dim1) wblog(FL,
      "ERR %s() index out of bounds (%d/%d)",FCT,j0,dim1);
   return v.init(dim2,rec(j0),ref);
};

template <class T>
wbMatrix<T>& wbMatrix<T>::getRecs(
   const WBINDEX &I, wbMatrix<T> &M 
 ) const {

   size_t i=0; M.init(I.len,dim2);
   for (; i<I.len; ++i) {
      if (I[i]>=dim1) wblog(FL,"ERR index out of bounds (%d/%d)",I[i],dim1);
      MEM_CPY<T>(M.data+i*dim2, dim2, data+I[i]*dim2);
   }
   return M;
};

template <class T>
wbMatrix<T>& wbMatrix<T>::getRecs(
   size_t i1, size_t i2, wbMatrix<T> &M,  
   size_t m  
 ) const {

    if (i1>=dim1 || i2>=dim1) wblog(FL,
       "ERR index out of bounds (%d,%d/%dx%d)",i1,i2,dim1,dim2);

    if (m==0) { 
       if (i1>i2) { M.init(0,dim2); return M; }
       M.init(i2-i1+1,dim2);
       MEM_CPY<T>(M.data, M.dim1*dim2, data+i1*dim2);
    }
    else {
       bool lflag=(int(m)<0); m=abs(int(m));
       if (m>dim2) wblog(FL,
          "ERR %s() m=%d/%d out of bounds",FCT,m,dim2);

       M.init(i1<=i2 ? i2-i1+1 : 0, m);
       if (M.dim1) {
          const T* d0=data+i1*dim2; if (lflag) d0+=(dim2-m);

          for (size_t i=0; i<M.dim1; ++i, d0+=dim2) {
             MEM_CPY<T>(M.data+i*m, m, d0);
          }
       }
    }

    return M;
};

template <class T> inline
wbMatrix<T>& wbMatrix<T>::getBlocks(  
    size_t j1, size_t j2, const size_t D,
    wbMatrix<T> &M
) const {

    size_t i,isnn, D2=2*D, s=D;
    T *p2,*p0;

    if (D*(j1+1)>dim2 || D*(j2+1)>dim2 || (D ? dim2%D : 0)) wblog(FL,
    "ERR invalid block parameters (%d,%d,%d,%d).",j1, j2, D, dim2);

    M.init(dim1,D2);

    isnn=(j2==(j1+1)); 
    if (isnn) s+=s;

    p0=data+D*j1; p2=M.data;
    for (i=0; i<dim1; i++) MEM_CPY<T>(p2+i*D2, s, p0+i*dim2);

    if (!isnn) {
       p0=data+D*j2; p2=M.data+D;
       for (i=0; i<dim1; i++) MEM_CPY<T>(p2+i*D2, s, p0+i*dim2);
    }

    return M;
}

template <class T>
void wbMatrix<T>::ColKron(
   const wbMatrix<T> &B,
   int d_ 
){
   if (&B==this) {
      wbMatrix<T> X(B); ColKron(X,d_);
      return;
   }

   size_t i,j, da=dim1, d2=dim2, db=B.dim1; 
   size_t d=( (d_<0) ? B.dim2 : (size_t)d_);

   if (d>B.dim2) wblog(FL,
      "ERR %s() index out of bounds (%d/%d)", FCT,d,B.dim2);

   if (B.isEmpty() || d==0) return;
   if (isEmpty()) { *this = B; return; }

   Repmat(db,1, 0,B.dim2); 

   T *dd, *d0;
   for (i=0; i<db; ++i) { 
      d0=B.data+i*B.dim2; dd = data + i*da*dim2 + d2;
      for (j=0; j<da; ++j, dd+=dim2) MEM_CPY<T>(dd,d,d0);
   }
};

template <class T>
wbMatrix<T>& wbMatrix<T>::repmat(
   size_t m, size_t n, wbMatrix<T> &Q,
   size_t pad1, size_t pad2
) const {

   if (this==&Q) {
      wbMatrix<T> X(*this);
      return X.repmat(m,n,Q,pad1,pad2);
   }

   size_t r,i,j;
   size_t d1=dim1+pad1, d2=dim2+pad2, D1=m*d1, D2=n*d2;

   Q.init(D1,D2);

   for (i=0; i<m; i++) 
   for (j=0; j<n; j++) {
      T *d = Q.data + (i*n*d1+j)*d2; 
      for (r=0; r<dim1; r++, d+=D2) MEM_CPY<T>(d, dim2, data+r*dim2);
   }

   return Q;
}

template <class T> inline
wbMatrix<T>& wbMatrix<T>::getBlock(  
    size_t j1, const size_t D, wbMatrix<T> &M
) const {
    T *d, *d0=data+D*j1;

    if (D*(j1+1)>dim2 || (D ? dim2%D : 0)) wblog(FL,
       "ERR %s() invalid block parameters (%d+%d/%d)",FCT,j1,D,dim2);
    M.init(dim1,D); d=M.data;

    for (size_t i=0; i<dim1; i++, d+=D, d0+=dim2)
    MEM_CPY<T>(d,D,d0); 

    return M;
};

template <class T> inline
wbMatrix<T>& wbMatrix<T>::getBlock(
   size_t i, size_t j, 
   size_t n, size_t m, 
   wbMatrix<T> &M
) const {

   if (i+n>dim1 || j+m>dim2) wblog(FL,
      "ERR %s() index out of bounds (%d..%d/%d; %d..%d/%d)",
       i+1,i+n,dim1,j+1,j+m,dim2
   );

   M.init(n,m); Wb::cpyStride(M.data,ref(i,j),m,n,-1,dim2);
   return M;
};

template <class T>
template <class T2>
wbMatrix<T>& wbMatrix<T>::getBlock(
    size_t k,      
    size_t D,      
    const wbMatrix<T2>& R, size_t k2, size_t D2, 
    wbMatrix<T> &M   
) const {

    size_t i,j;
    const T *d0; const T2 *d2; T *d;

    if (D*(k+1)>dim2 || !D || dim2%D
     || D2*(k2+1)>R.dim2 || !D2 || R.dim2%D2) wblog(FL,
       "ERR invalid block setting (%d,%d/%d; %d,%d/%d).",
        k,D,dim2, k2,D2,R.dim2);
    if (dim1!=R.dim1) wblog(FL,
       "ERR %s() dimension mismatch (%d/%d)",FCT,dim1,R.dim1);

    M.init(dim1,D+D2); d=M.data; d0=data+D*k; d2=R.data+D2*k2;

    for (i=0; i<dim1; i++) {
        MEM_CPY<T>(d,D,d0); d0+=dim2; d+=D;
        for (j=0; j<D2; j++) d[j]=T(d2[j]);
        d2+=R.dim2; d+=D2;
    }

    return M;
};

template <class T> inline
void wbMatrix<T>::setBlock(
   size_t k, const size_t D,
   const T* d0 
){
   T *d = data+k*D;

   if ((k+1)*D>dim2 || (D ? dim2%D : 0)) wblog(FL,
   "ERR invalid block parameters (%d,%d,%d).",k, D, dim2);

   for (size_t i=0; i<dim1; i++, d+=dim2)
   MEM_CPY<T>(d,D,d0); 
};

template <class T> inline
void wbMatrix<T>::SetBlock(
   size_t i0, size_t j0, 
   const wbMatrix<T> &M
){
   if (i0+M.dim1>dim1 || j0+M.dim2>dim2) wblog(FL,"ERR %s() index "
      "out of bounds (%d/%d)",FCT,i0+1,M.dim1,dim1,j0+1,M.dim2,dim2);
   const T* d0=M.data; T* d=ref(i0,j0);

   for (size_t i=0; i<M.dim1; i++, d+=dim2, d0+=M.dim2)
   MEM_CPY<T>(d,M.dim2,d0);
};

template <class T>
wbMatrix<T>& wbMatrix<T>::blockSum(size_t D, wbMatrix<T> &Q) const {

   size_t i,r,rk;
   const T *q0; T *q;

   if (&Q==this) wblog(FL,"ERR Output space same as input space!");
   if (!D || dim2%D) wblog(FL,"ERR Invalid block size (%d/%d)",D,dim2);

   rk=dim2/D; if (rk==1) { Q=(*this); return Q; }

   Q.init(dim1,D);

   for (i=0; i<dim1; i++) {
      q0=rec(i); q=Q.rec(i);
      for (r=0; r<rk; r++, q0+=D) Wb::addRange(q0,q,D);
   }

   return Q;
}

template <class T> inline
void wbMatrix<T>::blockSum( 
   size_t j1,
   size_t j2,
   const size_t D,    
   wbMatrix<double> &Q
) const {

   size_t i,k;
   double *q1, *q2, *q;

   if (D*(j1+1)>dim2 || D*(j2+1)>dim2 || (D ? dim2%D : 0)) wblog(FL,
   "ERR Invalid block parameters (%d,%d,%d,%d).",j1, j2, D, dim2);

   Q.init(dim1,D);
   for (k=0; k<dim1; k++) {
       q=rec(k); q1=q+D*j1; q2=q+D*j2;
       q=Q.rec(k);
       for (i=0; i<D; i++) q[i]=q1[i]+q2[i];
   }
}

template <class T> inline
void wbMatrix<T>::blockSum(
   const wbindex &J,
   const size_t D,    
   wbMatrix<double> &Q
) const {

   size_t i,j,r=(D ? dim2/D : 0);
   double *q;

   if (D ? dim2%D : 0) wblog(FL,"ERR Invalid D=%d (%d)", D, dim2);
   for (i=0; i<J.len; i++) if (J[i]>=r) 
   wblog(FL,"ERR Block index out of range (%d; %d,%d)",J[i],D,dim2);

   Q.init(dim1,D);

   for (i=0; i<dim1; i++) {
      q=Q.rec(i);
      for (j=0; j<J.len; j++) Wb::addRange(Q.ref(i,J[j]*D), q, D);
   }
}

template <class T> inline
void wbMatrix<T>::recSet(size_t k,
    const T* a, const wbindex &ia,
    const T* b, const wbindex &ib
){
    if (ia.len+ib.len!=dim2) wblog(FL,"ERR %s() size mismatch "
       "%dx%d/(%d: %d+%d)",FCT,dim1,dim2,k,ia.len,ib.len);
    if (k>=dim1) wblog(FL,
       "ERR %s() index out of bounds (%dx%d: %d)",FCT,dim1,dim2,k);

    if (dim2) { size_t i,l=0; T *d=rec(k);
       for (i=0; i<ia.len; ++i, ++l) { d[l]=a[ia[i]]; }
       for (i=0; i<ib.len; ++i, ++l) { d[l]=b[ib[i]]; }
    }
};

template <class T> inline
void wbMatrix<T>::recSetB(
    size_t r, size_t j, size_t D,
    const T *q0, const T *q2, const T a 
){
    T *q = data + r*dim2 + j*D;

    if (r>=dim1 || (j+1)*D>dim2 || (D ? dim2%D : 0)) wblog(FL,
       "ERR setRecB() index out of bounds (%d/%d; %d/%d %d)",
        r,dim1, j,dim2, D);

    if (q2==nullptr)
       for (size_t i=0; i<D; i++) q[i]=q0[i];
    else if (a) {
       if (a==+1) for (size_t i=0; i<D; i++) q[i]=q0[i]+q2[i]; else
       if (a==-1) for (size_t i=0; i<D; i++) q[i]=q0[i]-q2[i];
       else       for (size_t i=0; i<D; i++) q[i]=q0[i]+a*q2[i];
    }
}

template <class T> inline
void wbMatrix<T>::recSetB(
    size_t r, const wbindex &J, size_t D, const T *d0
){
    size_t i,j,k; const widx_t *jp=J.data;
    T *d = data + r*dim2;

    if (r>=dim1) wblog(FL,
    "ERR %s() index out of bounds (%d/%d)",FCT,r,dim1);
    if (J.len*D!=dim2) wblog(FL,
    "ERR %s() size mismatch ([%s] %d/%d)",FCT,STR(J+1),D,dim2);

    if (D==1) {
       for (j=0; j<J.len; j++) d[j]=d0[jp[j]];
    }
    else {
       for (j=0; j<J.len; j++, d+=D) { k=jp[j]*D;
          for (i=0; i<D; i++) { d[i]=d0[k+i]; }
       }
    }
}

template <class T> inline
void wbMatrix<T>::recSet(size_t i1, size_t i2) {

    if (i1>=dim1 || i2>=dim1) wblog(FL,
       "ERR index out of bounds (%d,%d; %d)", i1, i2, dim1);

    if (i1!=i2)
    MEM_CPY<T>(data+i1*dim2, dim2, data+i2*dim2);
}

template <class T> inline
void wbMatrix<T>::recSet(size_t i, const wbvector<T> &v) {

   if (i>=dim1) wblog(FL,"ERR %s() index out of bounds (%d/%d)",FCT,i,dim1);
   if (v.len!=dim2) wblog(FL,"ERR size mismatch (%d/%d)",v.len,dim2);
   MEM_CPY<T>(data+i*dim2, dim2, v.data);
};

template <class T>
template <class T2>
void wbMatrix<T>::recSetT(size_t i, const wbvector<T2> &v) {

   size_t j=0, n=(v.len<dim2 ? v.len : dim2); T *d = data + i*dim2;

   if (i>=dim1) wblog(FL,"ERR %s() index out of bounds (%d/%d)",FCT,i,dim1);
   if (v.len>dim2) wblog(FL,"ERR size mismatch (%d/%d)",v.len,dim2);

   for (; j<n; ++j) { d[j]=T(v.data[j]); }
   for (; j<dim2; ++j) { d[j]=T(0); }
};

template <class T> inline
void wbMatrix<T>::recSet (
    size_t i, const wbMatrix<T> &M, size_t j
){
    if (i>=dim1 || j>=M.dim1) wblog(FL,
       "ERR index out of bounds (%d/%d; %d/%d)", i, dim1, j, M.dim1);
    if (M.dim2!=dim2) wblog(FL,
       "ERR incompatible objects (%d/%d)", M.dim2, dim2);

    MEM_CPY<T>(data+i*dim2, dim2, M.data+j*M.dim2);
}

template <class T> inline 
void wbMatrix<T>::recSet (
    size_t i, const wbvector<T> &v1, const wbvector<T> &v2
){
    if (i>=dim1) wblog(FL,
       "ERR index out of bounds (%d/%d)", i, dim1);
    if (v1.len+v2.len!=dim2) wblog(FL,
       "ERR incompatible objects (%d+%d; %d)", v1.len, v2.len, dim2);

    MEM_CPY<T>(data+i*dim2, dim2, v1.len, v1.data, v2.data);
}

template <class T>
void wbMatrix<T>::recSet(size_t i, const T* v, T x) {
    if (i>=dim1) wblog(FL,
       "ERR record index out of bounds (%d/%s)",i,SSTR(*this));
    T *d=data+i*dim2; i=0; 
    if (d!=v) {
       if (x==T(+1)) { for (; i<dim2; ++i) d[i]= v[i]; } else
       if (x==T(-1)) { for (; i<dim2; ++i) d[i]=-v[i]; } else
       if (x==T( 0)) { for (; i<dim2; ++i) d[i]= T(0); }
       else          { for (; i<dim2; ++i) d[i]=x*v[i];}
    }
};

template <class T> inline
void wbMatrix<T>::recSetP(size_t i1, const T* v, size_t n) {
    if (i1>=dim1) wblog(FL,
       "ERR record index out of bounds (%d/%d)",i1,dim1);
    if (long(n)<0) n=dim2; else if (n>dim2) wblog(FL,
        "ERR %s() rec-length out of bounds (%d/%d)",FCT,n,dim2);
    T *p=data+i1*dim2; if (p!=v) MEM_CPY<T>(p,n,v);
};

template <class T> inline
void wbMatrix<T>::recSetP(size_t i1,
    const T *v1, size_t n1,
    const T *v2, size_t n2
){
    if (i1>=dim1) {
       if (n1 || n2) wblog(FL,
          "ERR record index out of bounds (%d/%d)",i1,dim1);
       return;
    }
    if (n1+n2>dim2) wblog(FL,
       "ERR %s() size out of bounds (%d+%d/%d)",FCT,n1,n2,dim2);

    size_t i; T *p=data+i1*dim2;
    for (i=0; i<n1; ++i) { p[i]=v1[i]; }; p+=n1;
    for (i=0; i<n2; ++i) { p[i]=v2[i]; };
};

template <class T> inline
void wbMatrix<T>::recSetP(size_t l,
    const T *v1, const wbindex &i1, const T *v2, const wbindex &i2
){
    if (l>=dim1) {
       if (i1.len || i2.len) wblog(FL,
          "ERR record index out of bounds (%d/%d)",l,dim1);
       return;
    }
    if (i1.len+i2.len>dim2) wblog(FL,
       "ERR %s() size out of bounds (%d+%d/%d)",FCT,i1.len,i2.len,dim2);

    size_t i=0; T *x=data+l*dim2; if (i1.len) {
    for (; i<i1.len; ++i) { x[i]=v1[i1.data[i]]; }; x+=i1.len; i=0; }
    for (; i<i2.len; ++i) { x[i]=v2[i2.data[i]]; };
};

template <class T> inline
void wbMatrix<T>::recAddP(size_t i1, const T* v) {
    T *d = data+i1*dim2;
    for (size_t i=0; i<dim2; i++) d[i]+=v[i];
}

template <class T>
template<class T2> inline
void wbMatrix<T>::setCol(size_t k, const wbvector<T2> &v) {
    if (v.len!=dim1) wblog(FL,
       "ERR %s() size mismatch (%d/%d)",FCT,v.len,dim1);
    setCol(k,v.data);
};

template <class T>
template<class T2> inline
void wbMatrix<T>::setCol(size_t k, const T2 *v, T fac, size_t stride) {

    if (k>=dim2) {
       if (dim2) wblog(FL,
          "ERR %s() index out of bounds (%d/%d)",FCT,k,dim2);
       else wblog(FL,
          "ERR %s() got empty matrix (%dx%d/%d)",FCT,dim1,dim2,k);
    }

    T *d=data+k; size_t i=0;
    if (fac==1)
         { for (; i<dim1; ++i, d+=dim2, v+=stride) { (*d)= T(*v);    }}
    else if (fac==-1)
         { for (; i<dim1; ++i, d+=dim2, v+=stride) { (*d)=-T(*v);    }}
    else { for (; i<dim1; ++i, d+=dim2, v+=stride) { (*d)=fac*T(*v); }}
};

template <class T> inline
void wbMatrix<T>::setCol(size_t k, const T &x) { 

    if (k>=dim2) {
       if (dim2) wblog(FL,
          "ERR %s() index out of bounds (%d/%d)",FCT,k,dim2);
       else wblog(FL,
          "ERR %s() got empty matrix (%dx%d/%d)",FCT,dim1,dim2,k);
    }

    T *d=data+k;
    for (size_t i=0; i<dim1; ++i, d+=dim2) (*d)=x;
};

template <class T>
template <class T2> inline
void wbMatrix<T>::setColT (size_t k, const wbvector<T2> &v) {
    if (v.len!=dim1) wblog(FL,
       "ERR %s() size mismatch (%d/%d)",FCT,v.len,dim1);
    setColT(k,v.data);
};

template <class T>
template <class T2> inline
void wbMatrix<T>::setColT(size_t k, const T2 *v, T fac, size_t stride) {

    T *d=data+k; size_t i=0;
    if (k>=dim2) wblog(FL,
       "ERR %s() index out of bounds (%d/%d)",FCT,k,dim2);

    if (fac==1)
         { for (; i<dim1; ++i, d+=dim2, v+=stride) { (*d)= T(*v); }}
    else if (fac==-1)
         { for (; i<dim1; ++i, d+=dim2, v+=stride) { (*d)=-T(*v); }}
    else { for (; i<dim1; ++i, d+=dim2, v+=stride) { (*d)=fac*T(*v); }}
};

template <class T> inline
wbMatrix<T>& wbMatrix<T>::colTimes(size_t i, T x) {

   if (i>=dim2) wblog(FL,
      "ERR %s() index out of bounds (%s/%d)",FCT,SSTR(*this),i);
   if (x!=T(1)) { T *d=data+i; i=0; 
      if (x==T( 0)) { for (; i<dim1; ++i, d+=dim2) d[0]= T(0); } else
      if (x==T(-1)) { for (; i<dim1; ++i, d+=dim2) d[0]=-d[0]; }
      else          { for (; i<dim1; ++i, d+=dim2) d[0]*=x;    }
   }
   return *this;
};

template <class T> inline
wbMatrix<T>& wbMatrix<T>::rowTimes(size_t i, T x) {

   if (i>=dim1) wblog(FL,
      "ERR %s() index out of bounds (%d /%s",FCT,i,SSTR(*this));
   if (x!=T(1)) { T *d=data+i*dim2; i=0; 
      if (x==T( 0)) { for (; i<dim2; ++i) d[i]= T(0); } else
      if (x==T(-1)) { for (; i<dim2; ++i) d[i]=-d[i]; }
      else          { for (; i<dim2; ++i) d[i]*=x;    }
   }
   return *this;
};

template <class T>
size_t wbMatrix<T>::skipNanRecs(wbindex &I QS_UNUSED_VAR) {
   wblog(FL,"WRN %s() irrelevant for type <%s>",FCT,TSTR(T));
   return 0;
}

template <> 
size_t wbMatrix<double>::skipNanRecs(wbindex &I) {

   size_t i,j,k=0,skipped=0;
   double *d=data;
   I.init(dim1);

   for (i=0; i<dim1; i++, d+=dim2) {
      for (j=0; j<dim2; j++) if (std::isnan(d[j])) break;
      if (j<dim2) { skipped++; continue; }

      if (k<i) recSet(k,i);
      I[k++]=i;
   }

   if (k==0)
        { Resize(0,dim2); I.init(); }
   else { dim1=k; I.len=k; } 

   return skipped;
};

template <class T> inline
wbMatrix<T>& wbMatrix<T>::recPermute(
   wbMatrix<T> &B, const wbperm &P) const {

   if (P.len!=dim1 || P.isValidPerm()<=0) wblog(FL,
      "ERR %s() invalid permutation [%s; %d]",FCT,STR(P),dim1);
   if (&B==this) wblog(FL,"ERR %s() got same object!",FCT);

   if (B.dim1!=dim1 || B.dim2!=dim2) { 
      B.init(dim1,dim2);
   }

   if (P.inv) { 
      const T *d=data;
      for (size_t i=0; i<P.len; ++i, d+=dim2)
      MEM_CPY<T>(B.data+P[i]*dim2, dim2, d); 
   }
   else {
      T *d=B.data;
      for (size_t i=0; i<P.len; ++i, d+=dim2)
      MEM_CPY<T>(d, dim2, data+P[i]*dim2); 
   }

   return B;
};

template <class T> inline
wbMatrix<T>& wbMatrix<T>::colPermute(
   wbMatrix<T> &B, const wbperm &P) const {

   char q=P.relevant(FL);

   if (q>1) wblog(FL,"ERR %s() got P=%s",FCT,STR(P));
   if (P.len>dim2) wblog(FL,
      "ERR %s() P.len=%d/%d out of bounds",FCT,P.len,dim2);

   if (!(q&1)) { B=*this; }
   else {
      size_t i,j;
      const wperm_t *p=P.data; const T* d0; T* d2;

      B.init(dim1,dim2); d2=B.data; d0=data;

      if (P.inv) {
         for (i=0; i<dim1; ++i, d0+=dim2, d2+=dim2) {
            for (j=0; j<P.len; ++j) { d2[p[j]]=d0[j]; }
            for (   ; j<dim2;  ++j) { d2[  j ]=d0[j]; }
         }
      }
      else {
         for (i=0; i<dim1; ++i, d0+=dim2, d2+=dim2) {
            for (j=0; j<P.len; ++j) { d2[j]=d0[p[j]]; }
            for (   ; j<dim2;  ++j) { d2[j]=d0[  j ]; }
         }
      }
   }

   return B;
};

template <class T> inline
void wbMatrix<T>::blockPermute(
   wbMatrix<T> &B, wbperm P) const { 

   if (this==&B) {
      wbMatrix<T> X(*this); X.blockPermute(B,P);
      return;
   }

   size_t i,j,D; const T *d0=data; T *d;

   if (!P.len || P.isValidPerm()<=0) wblog(FL,
      "ERR invalid permutation (%d;%d)", P.len, dim2);
   if (dim2%P.len) wblog(FL,
      "ERR block dimension mismatch ( %d = %d*?? )", dim2, P.len);

   B.init(dim1,dim2); d=B.data; D=dim2/P.len;
   P.flatten();

   for (i=0; i<dim1; i++, d0+=dim2)
   for (j=0; j<P.len; j++, d+=D) MEM_CPY<T>(d, D, d0+D*P[j]);
};

template <class T> inline
wbMatrix<T>& wbMatrix<T>::cols2Front(
  const WBINDEX &I1, wbMatrix<T> &B
) const {

   wbperm P; P.init2Front(I1,dim2);
   return colPermute(B,P);
};

template <class T> inline
wbMatrix<T>& wbMatrix<T>::cols2End(
  const WBINDEX &I1, wbMatrix<T> &B
) const {

   wbperm P; P.init2End(I1,dim2);
   return colPermute(B,P);
};

template <class T>
WBIDXMAT& wbMatrix<T>::toBlockIndex(
   widx_t D, WBIDXMAT &II_,
   wbvector< wbMatrix<T> > *QQ_  
) const {

   if (!D || dim2%D)
   wblog(FL,"ERR %s() invalid blocksize (%d/%d)",dim2,D);

   size_t i,m=dim2/D;

   wbvector< wbMatrix<T> > QQ;
   wbvector< wbindex > II;
   WBINDEX dd;
   wbperm P;

   wbvector<  WBINDEX const* > Ip;
   wbperm iP(P,'i');

   QQ.init(m);
   II.init(m); Ip.init(m);

   for (i=0; i<m; i++) {
      getBlock(i,D,QQ[i]).groupRecs(P,dd);
      II[i].BlockIndex(dd).Permute(iP);
      Ip[i]=(&II[i]);
   }

   II_.CAT(2,Ip);
   if (QQ_) QQ.save2(*QQ_);

   return II_;
};

template <class T>
template <class T2>
void wbMatrix<T>::toBlockIndex(widx_t D,
   wbMatrix<T2> &R, widx_t DR_, 
   WBIDXMAT *IB_,       
   WBIDXMAT *I2_,       
   WBIDXMAT *SS_,       
   wbvector< wbMatrix<T> > *QQ_   
) const {

   if (dim1!=R.dim1) wblog(FL,
      "ERR %s() size mismatch (%d/%d)",FCT,dim1,R.dim1);
   if (!D || dim2%D) wblog(FL,
      "ERR %s() invalid blocksize (%d/%d)",FCT,dim2,D);

   if (dim1==0 || dim2==0) {
      if (IB_) { IB_->init(); }
      if (I2_) { I2_->init(); }
      if (SS_) { SS_->init(); }
      if (QQ_) { QQ_->init(); }; return;
   }

   size_t i,DR, m=dim2/D;

   wbvector< wbMatrix<T> > QQ;
   wbvector< WBINDEX > IB,I2,SB;
   WBINDEX dd;
   wbMatrix<T2> Ri;
   wbperm pp;

   if (int(DR_)<=0) DR=R.dim2/m; else DR=DR_;
   if (m*DR_!=R.dim2) wblog(FL,
      "ERR block size mismatch (expecting %d*%d = %d / %d)",
       m,DR_,m*DR_,R.dim2
   );

   QQ.init(m);
   IB.init(m); I2.init(m); SB.init(m);

   for (i=0; i<m; i++) {
      getBlock(i,D, R,i,DR, QQ[i])
     .groupRecs(pp,dd, D,
        IB[i], I2[i], SB[i], 
        'p' 
      );

   }

   if (IB_) { IB_->CAT(2,IB); }
   if (I2_) { I2_->CAT(2,I2); }
   if (SS_) { SS_->CAT(2,SB); }; if (QQ_) QQ.save2(*QQ_);
};

template <class T>
void wbMatrix<T>::recPrint(
    size_t k, const char *istr0, char mflag) const {

    wbvec<char> s(64);
    wbvector<T> d;

    if (k>=dim1) wblog(FL,
       "ERR %s() index out of bounds (%d/%d)",FCT,k,dim1);
    s.catf(FL,"%.32s.rec(%ld)", istr0, k);

    d.init(dim2, (*this).rec(k));
    d.print(s.data,mflag); 

    return;
};

template <class T>
wbstring wbMatrix<T>::rec2Str(
   size_t k,
   const char *fmt,
   const char *sep,
   size_t stride,
   const char *sep2
 ) const {

   if (k>=dim1) wblog(FL,
      "ERR %s() index out of bounds (%d/%d)",FCT,k,dim1);
   if (!dim2) return "";

   wbvector<T> d; d.init2ref(dim2, data+k*dim2);
   return d.toStrf(fmt,sep,stride,sep2);
};

template <class T>
wbstring wbMatrix<T>::toStr(
    const char *fmt0, const char *sep, const char *rsep,
    size_t stride, const char *sep2
) const {

    size_t j, i=16+dim1*dim2*16;
    wbstring s,fmt;

    if (isComplex()) { i*=2; } 
    if (i<128) { i=128; }
    s.init(i);

    if (fmt0 && fmt0[0])
         { fmt=fmt0; }
    else { fmt.init2Fmt((T)0); }

    for (i=0; i<dim1; ++i) { if (i) s.push(FL,rsep);
    for (j=0; j<dim2; ++j) { if (j) s.push(FL,sep);
        s.pushf(FL,fmt.data, data[i*dim2+j]);
        if (stride && ((j+1)%stride)==0 && j+1<dim2) {
           s.push(FL,sep2);
        }
    }}

    return s;
};

template <>
wbstring wbMatrix<wbcomplex>::toStr(
   const char *fmt0, const char *sep, const char *rsep,
   size_t stride, const char *sep2
) const {

   wbstring s; 
   unsigned i,j; char fmt[16], zfmt[8]="%s", s1[32];

   snprintf(fmt,16,"%s",fmt0 && fmt0[0] ? fmt0 : "%.4g");

   for (i=0; fmt[i]; ++i) { if (fmt[i]=='%') {
      if (isdigit(fmt[++i])) { unsigned l=fmt[i]-'0';
         for (j=i+1; fmt[j]; ++j) {
            if (isdigit(fmt[j])) { l=10*l+(fmt[j]-'0'); } else break;
         }
         snprintf(zfmt,8,"%%%ds ",l); 
      }
      break;
   }}

   s.init(32+(4+i)*dim1*dim2); 

   for (i=0; i<dim1; ++i) { if (i) s.push(FL,rsep);
   for (j=0; j<dim2; ++j) { if (j) s.push(FL,sep);
       snprintf(s1,32,"%s",(*this)(i,j).toStr(fmt).data);
	   s.pushf(FL,zfmt,s1);
	   if (stride && ((j+1)%stride)==0 && j+1<dim2) {
		  s.push(FL,sep2);
	   }
   }}

   return s;
};

template <class T>
void wbMatrix<T>::toMxStruct(mxArray* S, const char *vname) const {

   mxArray *a=Mx::Array<double>(dim1,dim2,'t').copyFromTR(data);

   int fid = mxAddField2Scalar(FL,S,vname);
   mxSetFieldByNumber(S,0,fid, a);
};

template <class T>
void wbMatrix<T>::mat2mxs(mxArray* S, char fid) const {
   int n;

   if (!S) wblog(FL,"ERR %s() got null",FCT);
   if (!mxIsStruct(S)) wblog(FL,
      "ERR %s() got non-structure %s",FCT,mxGetClassName(S));
   if ((n=mxGetNumberOfElements(S))!=1) wblog(FL,
      "ERR %s() single struct required (%d)",FCT,n);
   if (fid>=(n=mxGetNumberOfFields(S))) wblog(FL,
      "ERR %s() field index out of bounds (%d/%d)",FCT,fid,n);

   mxArray *a=Mx::Array<double>(dim2,dim1).copyFromTR(data);
   mxSetFieldByNumber(S,0,fid, a);
};

#ifdef QS_USING_OMP

template <class T>
void getSortPerm_OMP(
   const wbMatrix<T> &A, wbperm &P, char dir, char lex
){
   if (!A.dim2) { P.init(); return; }

   wbRecs<T> R(A,P,dir,lex); 

 #ifdef WB_CLK_SPARSE
   Wb::Clock clk("mat:sort:recs-omp",0);
 #endif

   if (A.dim1<128) {
      sort(P.data, P.data+A.dim1, R); 
      return;
   }

   int id=0, i,l; 
   int np=1;
   if (!omp_in_parallel()) { 
      l=(1 << unsigned(floor(log2(double(A.dim1)))-6)); 
      np=MAX( QSP_NUM_THREADS, OMP_NUM_THREADS );
      np=MIN(l,np);
   }

   double nsub=double(A.dim1)/np;

   wbindex idx(np+1); wbperm PX(P.len);
   for (i=0; i<np; ++i) { idx[i]=size_t(i*nsub+0.5); }
   idx[i]=A.dim1;

  #pragma omp parallel for 
   for (int i=0; i<np; ++i) {
      size_t i1=idx.data[i], i2=idx.data[i+1];
      sort(P.data+i1, P.data+i2, R); 
      id=MAX(id,omp_get_thread_num());
   }

   while (idx.len>2) { int m2=np/2; P.swap(PX); 

     #pragma omp parallel for 
      for (int i=0; i<m2; ++i) { 
         size_t k=2*i, i1=idx.data[k], i2=idx.data[k+1], i3=idx.data[k+2];
         merge( 
            PX.data+i1, PX.data+i2,
            PX.data+i2, PX.data+i3, P.data+i1, R
         );
         id=MAX(id,omp_get_thread_num());
      }

      if (np%2) {
         size_t i1=idx.data[np-3], i2=idx.data[np-1], i3=idx.data[np];
         Wb::MemCpy(PX.data+i1,P.data+i1,(i2-i1)); 
         merge( 
            PX.data+i1, PX.data+i2,
            PX.data+i2, PX.data+i3, P.data+i1, R
         );
         idx[np-1]=idx[np]; idx.len=(np--);
      }

      for (l=1, i=2; i<np; i+=2, ++l) { idx[l]=idx[i]; }
      idx[l]=idx[np]; np=l; idx.len=l+1; 
   }

};

#endif

template <class T> inline
wbMatrix<T>& wbMatrix<T>::SortRecs(
   wbperm &P, char dir, char lex
){
   static int use_omp=-1;

   P.init(); if (!dim1) { return *this; }

#ifdef LD_CLEBSCH_QS
   if (use_omp<0) {
      use_omp=( 
         sp_num_threads > 0 ? 1 : 0 
      );

      if (use_omp && CG_VERBOSE>6)
      wblog(FL,"OMP %s() using openMP on top of STL",FCT);
   }

   char isLarge=(
      (dim1>(1<<30) && CG_VERBOSE>5) ||
      (dim1>(1<<28) && CG_VERBOSE>8) );

   if (isLarge) {
      if (use_omp>0) sprintf_str(
         "(parallel mode @ %d; %d)",omp_get_max_threads(),use_omp);
      else strcpy (str,"(serial mode)");
      wblog(FL,"TST %s() got %.3g x %d entries %s",
      FCT,double(dim1),int(dim2),str);
   }
#else
   if (use_omp<0) use_omp=0;
#endif

   if (use_omp<=0 || dim1<128) {
      Wb::hpsort(data,dim2,dim1, P,dir,lex);
      return *this; 
   }
   else {

     #ifdef QS_USING_OMP
      getSortPerm_OMP(*this,P,dir,lex);
     #else
      wblog(FL,"ERR %s()",FCT); 
     #endif

      return recPermute(P); 
   }

#ifdef LD_CLEBSCH_QS
   if (isLarge) wblog(FL,"TST %s() done",FCT);
#endif
   return *this;
};

template <class T>
template <class T2>
void wbMatrix<T>::groupRecs_rdeg(
   wbperm &P, WBINDEX &D, const wbMatrix<T2> &R,
   wbvector<widx_t> *Ib, 
   wbvector<widx_t> *I2, 
   wbvector<widx_t> *Sb, 
   char iflag 
){
   if (dim1==0) { P.init(); D.init(); return; }

   SortRecs(P);
   GroupSortedRecs(D);

   if (R.dim1!=P.len) wblog(FL,
      "ERR %s() dimension mismatch (%d/%d)",FCT,R.dim1,P.len);
   if (R.dim2==0) { wblog(FL,
      "WRN %s() got empty reference space",FCT); return;
   }

   size_t i,j,d,l=0, m=D.len;
   WBINDEX dd,v, I(D.max());
   wbvector< WBINDEX > SS;
   wbvector< wbindex > II;
   wbMatrix<T2> Q;
   wbperm pp, iP;

   if (I2 || Sb) { II.init(m); if (Sb) SS.init(m); }

   for (i=0; i<m; ++i, l+=d) { d=D[i];
       if (d>1) {
          I.len=d; for (j=0; j<d; j++) I[j]=P[l+j];
          R.getRecs(I,Q).groupRecs(pp,dd);
          for (j=0; j<d; j++) P[l+j]=I[pp[j]];

          if (I2 || Sb) {
             II[i].BlockIndex(dd); if (Sb)
             SS[i].init(pp.len).set(dd.len);
          }
       }
       else {
          if (d!=1) wblog(FL,"ERR d=%g !?",d); 
          if (I2 || Sb) {
             II[i].init(1).set(0); if (Sb) {
             SS[i].init(1).set(1); }
          }
          continue;
       }
   }

   if (iflag) { P.invert(iP); }

   if (Ib) {
      wbindex J; J.BlockIndex(D); J.save2(*Ib);
      if (Ib->len!=P.len) wblog(FL,"ERR size mismatch %d/%d",Ib->len,P.len);
      if (iflag) Ib->Permute(iP); 
   }

   if (I2) {
      wbvector< WBINDEX* > J(II.len);
      for (i=0; i<II.len; ++i) { J[i]=(WBINDEX*)(&II[i]); }; I2->Cat(J);
      if (I2->len!=P.len) wblog(FL,"ERR size mismatch %d/%d",I2->len,P.len);
      if (iflag) I2->Permute(iP); 
   }

   if (Sb) {
      Sb->Cat(SS);
      if (Sb->len!=P.len) wblog(FL,"ERR size mismatch %d/%d",Sb->len,P.len);
      if (iflag) Sb->Permute(iP); 
   }
};

template <class T>
void wbMatrix<T>::groupRecs(
   wbperm &P, WBINDEX &D, size_t nc, char lex,
   WBINDEX &Ib, 
   WBINDEX &I2, 
   WBINDEX &Sb, 
   char iflag   
){
   if (!dim1) { P.init(); D.init(); return; }

   SortRecs(P);

   if (nc>dim2) wblog(FL,"ERR dimension out of bounds (%d/%d)",nc,dim2);
   GroupSortedRecs(D,nc,lex,Ib,I2,Sb);

   if (iflag) { wbperm iP(P,'i');
      Ib.Permute(iP);
      I2.Permute(iP);
      Sb.Permute(iP);
   }
};

template <class T>
void wbMatrix<T>::groupRecs(
   wbperm &P, WBINDEX &D,
   size_t m,
   char lex,
   wbindex *Ig,
   wbMatrix<T> *X
){
   if (!dim1) { P.init(); D.init(); return; }
   if (!dim2) {
      P.init(dim1); D.init(1); D[0]=dim1;
      if (Ig) { Ig->init(dim1).set(0); }

      if (int(m)<0) dim1=1;
      else if (m) wblog(FL,
         "ERR %s() m=%d out of bounds (%dx%d)",FCT,m,dim1,dim2);
      return;
   }

 #ifdef WB_CLK_SPARSE
   Wb::Clock clk("mat:sort:recs",0);
 #endif

   if (X) {
     if (dim1>1) { X->init(*this); } else { X=nullptr; }
   }

   if (!isSorted(+1,lex)) 
        { SortRecs(P,+1,lex); } 
   else { P.init(dim1); }

 #ifdef WB_CLK_SPARSE
   clk.Switch("mat:group:recs"); 
 #endif

   if (int(m)<0) { GroupSortedRecs(D,-1,lex); } 
   else if (m) { groupSortedRecs(D,m,lex); }
   else { D.init(1); D[0]=dim1; }

   if (X && (dim1==X->dim1)) { X->save2(*this);
      P.init(dim1);
      D.init2val(dim1,1);
   }

   if (Ig) {
      size_t i,j,l=0,d=0; widx_t *ig; wperm_t *p=P.data;
      Ig->init(P.len); ig=Ig->data;
      for (i=0; i<D.len; ++i, l+=d) { d=D.data[i];
      for (j=0; j<d; ++j) ig[p[l+j]]=i; }
   }
};

template <class T>
void wbMatrix<T>::groupRecs(wbperm &P, WBINDEX &D,
   const WBINDEX *I, char lex
){
   if (dim1==0) { P.init(); D.init(); return; }

   SortRecs(P,+1,lex); 

   if (!I) { GroupSortedRecs(D,-1,lex); }
   else {
      wbMatrix<T> X; this->cols2Front(*I,X);
      X.groupSortedRecs(D,I->len,lex);
   }
};

template <class T>
wbMatrix<T>& wbMatrix<T>::ReduceBlocks(const WBINDEX &D) {

   if (!D || !*this) { return init(); }
   if (D.len!=dim1) {
      for (size_t i=0, ig=0; ig<D.len; i+=D[ig], ++ig) { 
         if (i>=dim1) wblog(FL,"ERR %s() "
           "dims out of bounds (ig=%d/%d: D=%d/%d)",FCT,ig,i,D.sum(),dim1);
         if (ig<i) {
         MEM_CPY<T>(data+ig*dim2, dim2, data+i*dim2); }
      }
      Resize(D.len,dim2);
   }
   else if (D.anyUnequal(1)) wblog(FL,
     "ERR %s() size mismatch (%d/%d)",FCT,D.sum(),dim1);
   return *this;
};

template <class T>
wbMatrix<T>& wbMatrix<T>::GroupSortedRecs(WBINDEX &D, size_t m, char lex) {
   groupSortedRecs(D,m,lex); 
   ReduceBlocks(D);
   return *this;
};

template <class T>
WBINDEX& wbMatrix<T>::groupSortedRecs( 
   WBINDEX &D, size_t m, char lex) const { 

   size_t i,ig;
   char c, cref=0;

   if (int(m)<0) { m=dim2; } else
   if (m>dim2) wblog(FL,"ERR size out of bounds (%d/%d)",m,dim2);
   if (!dim2 || !m) {
      if (!dim1) wblog(FL,"WRN %s() for %dx%d matrix",FCT,dim1,dim2);
      D.init(1); D[0]=dim1; return D;
   }

   D.init(dim1); if (!dim1) { return D; }

   ig=0; ++D[ig];
   for (i=1; i<dim1; ++i) { c=recCompare(i,i-1,m,lex);
      if (c) { ++ig;
         if (!cref) { cref=c; } else
         if (c!=cref) wblog(FL,
            "ERR recs not sorted (%d: %d/%d, m=%d/%d)",i,c,cref,m,dim2);
      }
      ++D[ig];
   }

   return D.Resize(ig+1);
};

template <class T>
WBINDEX& wbMatrix<T>::groupSortedRecs(
   WBINDEX &D, size_t m, char lex,  
   WBINDEX &Ig, 
   WBINDEX &I2, 
   WBINDEX *D2  
) const {
   size_t i,j=1,i2=0,i0=0, ig=0; widx_t *d2=nullptr;
   char c, c2=0, cref=0;

   if (int(m)<0) { m=dim2; }
   else if ((!m && dim2) || m>dim2) wblog(FL,
      "ERR number out of bounds (%d/%d)",m,dim2);

   if (!dim1 || !dim2) {
      D.init(); Ig.init(); I2.init(); if (D2) { D2->init(); }
      return D;
   }

   Ig.init(dim1); D.init(dim1);
   I2.init(dim1); 
   if (D2) { D2->init(dim1); d2=D2->data; }

   for (++D[ig], i=1; i<dim1; ++i, ++D[ig]) {
      c =recCompare(i,i-1, m,lex); if (!c && m<dim2) {
      c2=recCompare(j,j-1,-1,lex); } else { c2=c; }
      if (c2) {
         if (!cref) { cref=c2; } else
         if (c2!=cref) { wblog(FL, 
            "ERR input not sorted (%d: %d/%d/%d)",i,c,c2,cref);
         }
      }
      if (c) {
         if (d2) { for (j=i0; j<i; ++j) { d2[j]=i2; }}
         ++ig; i0=i; i2=0;
      }
      I2[i]=i2++;
      Ig[i]=ig;
   }
   if (d2) { for (j=i0; j<i; ++j) { d2[j]=i2; }}

   D.len=ig+1; 

   return D;
};

template <> inline
wbMatrix<double>& wbMatrix<double>::sortRecs_float(wbperm &P, char dir) {

   if (dim1==0) { P.init(); return *this; }

   char lex=1;
   wbMatrix<double> X(*this); X.SkipTiny_float();

   X.SortRecs(P,dir,lex); 

   return Set2Recs(P);
};

template <class T>
char wbMatrix<T>::findRecsInSet(
    const wbMatrix<T> &S, wbvector<int> &I) const {

    size_t j,l,k,d, e=0; int ns=(int)S.dim1, ifound;
    WBINDEX dg;
    wbMatrix<T> QQ(*this);
    wbperm is;

    QQ.groupRecs(is,dg); I.init(dim1);

    for (l=k=0; k<QQ.dim1; k++) {
        for (ifound=0; ifound<ns; ifound++) {
            if (QQ.recEqual(k, S.rec(ifound)))
            break;
        }
        if (ifound>=ns) { ifound=-1; e++; }

        for (d=dg[k], j=0; j<d; j++)
        I[is[l++]]=ifound;  
    }

    return e;
}

template <class T>
size_t wbMatrix<T>::findUniqueRecSorted1(size_t n, T eps) const {

   if (long(n)>=0 && (!n || n>dim2)) wblog(FL,
      "ERR %s() invalid n=%d/%d",FCT,n,dim2);
   if (dim1<=1) { return (dim1 ? 0:-1); }

   if (eps==0) {
      char c=0, c0=0;
      for (size_t l=dim1-1, i=0; i<l; ++i) {
         c=recCompare(i,i+1,n);
         if (!c || !c0) { c0=c; continue; }
         else { return i; }
      }
      if (dim1==2 && c) return 0;
   }
   else {
      T c=0, c0=0;
      for (size_t l=dim1-1, i=0; i<l; ++i) {
         c=recDiff2(i,i+1,n);
         if (c<eps || c0<eps) { c0=c; continue; }
         else { return i; }
      }
      if (dim1==2 && (c>eps)) return 0;
   }
   return -1;
};

template <class T>
int wbMatrix<T>::findRecSorted(const T* r, size_t n, char lex) const {

   if (int(n)>=0 && (!n || n>dim2)) wblog(FL,
      "ERR %s() invalid n=%d/%d",FCT,n,dim2);

   if (dim1<3) {
      for (size_t i=0; i<dim1; ++i) {
         if (!recCompareP(i,r,n,lex)) return i;
      }; return -1;
   }

   size_t k1=0, k2=dim1-1, k=(k2-k1)/2;
   char c, cref=recCompare(k1,k2,n,lex);

   if (cref==0) {
      wblog(FL,"WRN all recs the same (%dx%d) ???",dim1,dim2);
      return (recCompareP(0,r,n,lex) ? -1 : 0);
   }
   for (size_t m=MAX(size_t(1),(dim1-2)/4), i=1; i<dim1; i+=m)
   if ((c=recCompare(i-1,i,n,lex))==-cref) wblog(FL, 
      "ERR records not sorted (%d/%d) !?",c,cref);

   c=recCompareP(k1,r,n,lex); if (c!= cref) return (c ? -1 : 0);
   c=recCompareP(k2,r,n,lex); if (c!=-cref) return (c ? -1 : dim1-1);
   c=recCompareP(k ,r,n,lex);

   while (1) { 
      if (c==cref)
             { k1=k; k+=(k2-k1)/2; if (k==k1) return -1; } else
      if (c) { k2=k; k-=(k2-k1)/2; if (k==k2) return -1; } else
      return k; 

      if (recCompare(k1,k2,n,lex)==-cref) wblog(FL,
         "ERR records not sorted !?"); 
      c=recCompareP(k,r,n,lex);
   }

   return -1;
};

template <class T>
int wbMatrix<T>::findRec(const T* r, size_t n, char lex) const {

   for (size_t k=0; k<dim1; k++)
   if (!recCompareP(k,r,n,lex)) return k;

   return -1;
}

template <class T>
size_t wbMatrix<T>::UnionRecs(const wbMatrix<T> &B0, wbindex &Ia){
    size_t l=0,ia=0,ib=0, na=(dim2 ? dim1 : 0), nb=B0.dim1;
    char c;

    wbMatrix<T> A(*this), B(B0);
    wbperm P1,P2; wbindex I;

    if (A.isEmpty() || B.isEmpty()) { init(); return na; }

    if (A.dim2 != B.dim2) wblog(FL,
       "ERR %s() severe dimension mismatch (%d/%d)",
        FCT,A.dim2,B.dim2);

    A.SortRecs(P1); I.init(na);
    B.SortRecs(P2);

    while (ia<na && ib<nb) { c=A.recCompareP(ia, B.rec(ib));
       if (c<0) ia++; else
       if (c>0) ib++;
       else {
          do { I[l++]=ia++; } while (ia<na && A.recEqual(ia-1,ia));
          do {        ib++; } while (ib<nb && B.recEqual(ib-1,ib));
       }
    }

    if (l) I.len=l; else I.init();

    A.Set2Recs(I); P1.get(I,Ia);

    return (na-dim1);
}

template <class T>
bool wbMatrix<T>::gotRecOverlap(const wbMatrix<T> &B0) const {

    size_t ia=0,ib=0, na=(dim2 ? dim1 : 0), nb=B0.dim1;
    char c;

    wbMatrix<T> A(*this), B(B0);
    wbperm P1,P2;

    if (A.isEmpty() || B.isEmpty()) return 0;

    if (A.dim2!=B.dim2) wblog(FL,
    "ERR %s() severe dimension mismatch (%d/%d)",FCT,A.dim2,B.dim2);

    A.SortRecs(P1);
    B.SortRecs(P2);

    while (ia<na && ib<nb) { c=A.recCompareP(ia, B.rec(ib));
       if (c<0) ia++; else
       if (c>0) ib++; else return 1;
    }

    return 0;
};

template <class T>
bool wbMatrix<T>::gotRecOverlapSA(const wbMatrix<T> &B) const {

    size_t ia=0,ib=0, na=(dim2 ? dim1 : 0), nb=B.dim1;
    char c;

    if (isEmpty() || B.isEmpty()) return 0;

    if (dim2!=B.dim2) wblog(FL,
       "ERR %s() severe dimension mismatch (%dx%d/%dx%d)",
        FCT,dim1,dim2,B.dim1,B.dim2);

    while (ia<na && ib<nb) { c=recCompareP(ia, B.rec(ib));
       if (c<0) ia++; else
       if (c>0) ib++; else return 1;
    }

    return 0;
};

template <class T>
size_t wbMatrix<T>::findValsCol(
    const size_t k, const wbvector<T> &v, wbindex &Ia
) const {
    wbvector<T> vk; getCol(k,vk);
    vk.findValues(v,Ia);
    return Ia.len;
}

template <class T> inline
void matchIndexU(const char *F, int L,
    const wbMatrix<T> &QA,
    const wbMatrix<T> &QB,
    wbindex &IB, const char force
){
    size_t ma,mb; wbindex Ia, Ib;

    matchIndex(QA,QB,Ia,Ib,1,&ma,&mb);
    if (mb) wblog(F,L,"ERR index B not unique (%d)",mb);

    if (force) { wbperm iP;
       if (Ia.len!=QA.dim1) wblog(F,L,
       "ERR failed to find all records (%d/%d)",Ia.len,QA.dim1);

       getIPerm(Ia,iP); Ib.select(iP,IB);
    }
    else {
       IB.init(QA.dim1).set(size_t(-1));
       IB.Set(Ia,Ib);
    }
}

template <class T>
int matchIndex(
    const wbMatrix<T> &QA, const wbMatrix<T> &QB,
    wbindex &Ia, wbindex &Ib,
    char lex,
    widx_t *ma, widx_t *mb, 
    T eps
){
    size_t i;

    wbMatrix<T> Q1(QA), Q2(QB);
    wbperm P1, P2;
    wbindex ix1, ix2;

    if (QA.dim2!=QB.dim2 && QA.dim2 && QB.dim2) wblog(FL,
       "ERR %s severe dimension mismatch (%s <> %s)",
       FCT,SSTR(QA),SSTR(QB));

    if (QA.isEmpty() || QB.isEmpty()) {
       Ia.init(); Ib.init(); return 0;
    }

    Q1.SortRecs(P1);
    Q2.SortRecs(P2);

    i=matchSortedIdx(Q1,Q2,ix1,ix2,-1,lex,ma,mb,eps);

    P1.get(ix1,Ia);
    P2.get(ix2,Ib);

    return (int)i;
};

template <class T>
int wbMatrix<T>::getDiff(
    const wbMatrix<T> &B, wbindex &Ia, wbindex *Ib) const {

    wbMatrix<T> Q1(*this), Q2(B);
    wbperm P1,P2;
    wbindex ix1,ix2;

    if (dim2!=B.dim2) wblog(FL,
       "ERR %s() severe dimension mismatch (%d/%d)",
        FCT,dim2, B.dim2
    );

    if (isEmpty() || B.isEmpty()) {
       Ia.init(); if (Ib) Ib->init();
       return 0;
    }

    Q1.SortRecs(P1);
    Q2.SortRecs(P2);

    Q1.getDiffSorted(Q2,ix1, Ib ? &ix2 : nullptr);
    P1.get(ix1, Ia); if (Ib) {
    P2.get(ix2,*Ib); }

    return Ia.len;
};

template <class T>
int wbMatrix<T>::getDiffSorted(
   const wbMatrix<T> &B, wbindex &Ia, wbindex *IB
) const {

   size_t ia,ib,la,lb; widx_t *Ib=nullptr; char c;

   if (dim2!=B.dim2) wblog(FL,
      "ERR %s() dimension mismatch (%d/%d)",FCT,dim2,B.dim2);

   Ia.init(dim1);
   if (IB) { IB->init(B.dim1); Ib=IB->data; }

   for (la=lb=ia=ib=0; ia<dim1 && ib<B.dim1;) {
       c=recCompareP(ia, B.rec(ib));
       if (c<0) { Ia[la++]=ia++; } else
       if (c>0) { if (Ib) { Ib[lb++]=ib; }; ++ib; }
       else {
          while ((++ia)<dim1) {
             c=recCompare(ia-1,ia);
             if (c) { 
                if (c>0) wblog(FL,
                   "ERR %s() expecting ascending order (%d)",FCT,ia);
                break;
             }
          }
          while ((++ib)<B.dim1) {
             c=B.recCompare(ib-1,ib);
             if (c) { 
                if (c>0) wblog(FL,
                   "ERR %s() expecting ascending order (%d)",FCT,ib);
                break;
             }
          }
       }
   }

   while (ia<dim1) { Ia[la++]=ia++; }
   Ia.Resize(la);

   if (Ib) {
      while (ib<B.dim1) { Ib[lb++]=ib++; }
      IB->Resize(lb);
   }

   return Ia.len;
};

#endif

