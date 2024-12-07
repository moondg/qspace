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

#ifndef __WB_MATRIX_ROW_MAJOR_HH__
#define __WB_MATRIX_ROW_MAJOR_HH__

/* ---------------------------------------------------------------- *
 * wbMatrix - row-major matrix class
 * AW (C) Jan 2006
 * ---------------------------------------------------------------- */

template <class T>
class wbMatrix { 

  public:

    wbMatrix(size_t r=0, size_t c=0)
     : data(NULL), dim1(r), dim2(c), isdiag(0), isref(0) {
       if (r || c) NEW_DATA();
    };

    wbMatrix(size_t r, size_t c, T *d, const char ref=0)
     : data(NULL), dim1(r), dim2(c), isdiag(0), isref(0) {
       if (ref) { data=d; isref=1; }
       else if (r || c) NEW_DATA(d);
    };

    wbMatrix(const wbMatrix &M)
     : data(NULL), dim1(M.dim1), dim2(M.dim2), isdiag(M.isdiag), isref(0) {
       NEW_DATA(M.data);
    };

    template <class T2>
    wbMatrix(const wbMatrix<T2> &M)
     : data(NULL), dim1(0), dim2(0), isdiag(0), isref(0) {
       if (M.isdiag) wblog(FL,"ERR %s() got isdiag=%d",FCT,M.isdiag);
       initT(M);
    };

    wbMatrix(const mxArray *a, char tflag=0, char ref=0)
     : data(NULL), dim1(0), dim2(0), isdiag(0),isref(0) {
       init(FL,a,tflag,ref);
    };

    wbMatrix(const char *F, int L, const mxArray *a)
     : data(NULL), dim1(0), dim2(0), isdiag(0),isref(0) {
       init(F,L,a);
    };

    virtual ~wbMatrix() { if (data && !isref) {
       WB_DELETE(data);
    }};

    wbMatrix& init(size_t r=0, size_t c=0, const T* d=NULL) {
       if (isref) {
          isref=0; dim1=dim2=0; data=0; 
       }
       RENEW(r,c,d); return *this;
    };

    wbMatrix& init( 
       size_t r, size_t c, const T *d0, size_t stride
    ){
       unsigned i=0, j=0; T* d;

       if (int(stride)<0) wblog(FL,
          "ERR %s() got stride=%ld !?",FCT,stride);
       if (!d0) wblog(FL,"ERR %s() got null space !?",FCT);

       RENEW(r,c,NULL,0); d=data;
       for (; i<dim1; ++i, d+=dim2, d0+=stride) {
          for (j=0; j<dim2; ++j) d[j]=d0[j];
       }
       return *this;
    };

    wbMatrix& initIdentity(size_t n) { 
       size_t i=0, n2=n*n; RENEW(n,n,NULL,0);
       for (; i<n2; i+=(n+1)) { data[i]=T(1); }
    };

    wbMatrix& init2val(size_t r, size_t c, const T &d) {
       RENEW(r,c,NULL,0); set(d);
       return *this;
    };

    wbMatrix& initDef(size_t r=0, size_t c=0, T* d=NULL) {
       RENEW(r,c,NULL,0);
       if (d) for (size_t n=r*c, i=0; i<n; ++i) data[i]=d[i];
       return *this;
    };

    template<class T2>
    wbMatrix<T>& initT(size_t r, size_t c, T2* d) {
       RENEW(r,c,0,0);
       for (size_t n=r*c, i=0; i<n; ++i) data[i]=(T)d[i];
       return *this;
    };

    template<class T2>
    wbMatrix<T>& initT(const wbMatrix<T2> &B) {
       if (B.isdiag) wblog(FL,"ERR %s() got isdiag=%d",FCT,B.isdiag);
       RENEW(B.dim1,B.dim2,0,0);
       for (size_t n=dim1*dim2, i=0; i<n; ++i) data[i]=(T)B.data[i];
       return *this;
    };

    template<class T2>
    wbMatrix<T>& initT(const char *F, int L, const wbMatrix<T2> &B) {
       if (B.isdiag) wblog(FL,"ERR %s() got isdiag=%d",FCT,B.isdiag);
       RENEW(B.dim1,B.dim2,0,0);
       for (size_t n=dim1*dim2, i=0; i<n; ++i) { data[i]=(T)B.data[i];
          if (T2(data[i])!=B.data[i]) wblog(F,L,
             "ERR type conversion changes value (%g,%g)",
              double(B.data[i]), double(T2(data[i]))
          );
       }
       return *this;
    };

    template<class TA, class TB>
    wbMatrix<T>& init(const char *F, int L,
       const wbMatrix<TA> &B, const wbvector<TB> &V, char dim=2
    );

    template<class T2>
    wbMatrix<T>& initDim(const wbMatrix<T2> &B) { 
       RENEW(B.dim1, B.dim2);
       return *this;
    };

    wbMatrix& init(const wbMatrix &B) {
         RENEW(B.dim1, B.dim2, B.data); isdiag=B.isdiag;
         return *this;
    };

    wbMatrix& init( 
       const char *F, int L, const mxArray *a,
       char tflag=0, char ref=0, char tcheck=1); 

    void init(const mxArray *a) { init(FL,a); };

    void init0(const mxArray *a) {
       init(FL,a,1,0,1); 
    };

    void initTST() {           
        size_t i,j,k=0;      
        for (i=0; i<dim1; ++i) 
        for (j=0; j<dim2; ++j) data[k++]=10*(i+1)+(j+1);
    };

    wbMatrix& init2ref(size_t r, size_t c, const T* d);
    wbMatrix& init2ref(const wbvector<T> &v, const char tflag=0);
    wbMatrix& init2ref(const wbMatrix &M);

    wbMatrix& unRef();

    void no_ref_(const char *F, int L, const char *fct) {
       if (isref) { wblog(F_L,
          "ERR %s() got reference (isref=%d)",fct,isref);
       }
    };

    void swap(wbMatrix &B, char ref=0) { 
        if (this!=&B) {
           if (!ref && (isref || B.isref)) wblog(FL,
              "ERR %s() got arrays with isref=(%d,%d)",FCT,isref,B.isref);
           SWAP(data,   B.data  );
           SWAP(dim1,   B.dim1  );
           SWAP(dim2,   B.dim2  );
           SWAP(isref,  B.isref );
           SWAP(isdiag, B.isdiag);
        }
    };

    wbMatrix& save2(wbMatrix &A, char ref=0) {
       if (this!=&A) {
          if (!ref) { no_ref_(FLF); }
          this->swap(A,ref); if (data!=A.data) { init(); }
       }
       return A;
    };

    void copyStride(T* dd, size_t stride) const {
       T *d0=data;
       for (size_t i=0; i<dim1; ++i, d0+=dim2, dd+=stride)
       MEM_CPY<T>(dd, dim2, d0);
    };

    wbMatrix& operator= (const wbMatrix &B) {
       if (this!=&B) {
          RENEW(B.dim1, B.dim2, B.data);
          isdiag=B.isdiag;
       }
       return *this;
    };

    template<class TB>
    wbMatrix& operator= (const wbMatrix<TB> &B) {
       init(B.dim1, B.dim2);
       for (size_t i=0, n=B.numel(); i<n; ++i) { data[i]=T(B.data[i]); }
       isdiag=B.isdiag;
       return *this;
    };

    explicit operator bool() const {
       if (data)
            { if ( dim1 &&  dim2) { return 1; }}
       else { if (!dim1 || !dim2) { return 0; }}

       wblog(FL,"ERR unexpected matrix %gx%g (%p)",dim1,dim2,data);
       return 0; 
    };
    bool operator! () const { return !bool(*this); }

    bool operator==(const T &x) const {
       for (size_t n=dim1*dim2, i=0; i<n; ++i)
         { if (data[i]!=x) return 0; }
       return 1;
    };

    bool operator==(const wbMatrix &B) const {
       if (this!=&B) {
          if (dim1!=B.dim1 || dim2!=B.dim2) return 0;
          if (data!=B.data) {
             for (size_t n=dim1*dim2, i=0; i<n; ++i)
             if (data[i]!=B.data[i]) return 0;
          }
       }
       return 1;
    };

    bool operator!=(const wbMatrix &B) const {
       return !((*this)==B);
    };

    bool deepEqualP(const wbMatrix &B) const; 

    bool anyEqual(const T& x) const {
       for (size_t n=dim1*dim2, i=0; i<n; ++i) {
          if (data[i]==x) return 1; }
       return 0;
    };

    bool allEqual(const T& x) const {
       for (size_t n=dim1*dim2, i=0; i<n; ++i) {
          if (data[i]!=x) return 0; }
       return 1;
    };

    bool anyUnequal(const T& x) const { return !allEqual(x); };

    bool anyLT(const T& x) const {
       for (size_t n=dim1*dim2, i=0; i<n; ++i) {
           if (data[i]<x) return 1; }
       return 0;
    };

    bool anyGT(const T& x) const {
       for (size_t n=dim1*dim2, i=0; i<n; ++i) {
           if (data[i]>x) return 1; }
       return 0;
    };

    bool col_anyGT(size_t j, const T& x) const { 
       if (dim1) {
          if (j>=dim2) wblog(FL,
             "ERR %s() index out of bounds (%d/%d)",FCT,j,dim2);
          for (size_t i=0; i<dim1; ++i, j+=dim2) {
             if (data[j]>x) return 1; 
          }
       }
       return 0;
    };

    wbMatrix& operator+= (const wbMatrix &B) {
       size_t i,s=dim1*dim2;

       if (isdiag) if (!B.isdiag) isdiag=0;

       if (dim1!=B.dim1 || dim2!=B.dim2) wblog(FL,
       "ERR Dimension mismatch (%dx%d + %dx%d).",dim1,dim2,B.dim1,B.dim2);

       for (i=0; i<s; ++i) data[i]+=B.data[i];
       return *this;
    };

    wbMatrix operator+ (const wbMatrix &B) const {
       wbMatrix Mout(*this); Mout+=B; 
       return Mout;
    };

    wbMatrix operator+ (const T &x) const {
       wbMatrix Mout(*this); Mout+=x; 
       return Mout;
    };

    void operator-= (const wbMatrix &B) {
       if (isdiag) if (!B.isdiag) { isdiag=0; }
       if (dim1!=B.dim1 || dim2!=B.dim2) wblog(FL,
          "ERR size mismatch (%dx%d + %dx%d).",dim1,dim2,B.dim1,B.dim2);

       for (size_t i=0, n=dim1*dim2; i<n; ++i) {
           data[i]-=B.data[i];
       }
    };

    wbMatrix operator- (const wbMatrix &B) const {
       wbMatrix Mout(*this); Mout-=B;
       return Mout;
    };

    void operator*= (const T c) {
        size_t i, n=dim1*dim2;

        if (c==T(+1)) { return; } else
        if (c==T( 0)) { set((T)0); } else
        if (c==T(-1)) { for (i=0; i<n; ++i) data[i]=-data[i]; }
        else          { for (i=0; i<n; ++i) data[i]*=c; }
    };

    wbMatrix& operator+= (const T c) { if (c) {
       for (size_t n=dim1*dim2, i=0; i<n; ++i) data[i]+=c; }
       return *this;
    };

    void operator-= (const T c) { if (c) {
       for (size_t n=dim1*dim2, i=0; i<n; ++i) data[i]-=c; }
    };

    wbMatrix  operator* (const T x) const {
        wbMatrix X(*this); if (x!=1) X*=x; 
        return X;
    };

    wbMatrix& add  (const wbMatrix &B, const T& c=1, const char zflag=0);
    wbMatrix& minus(const wbMatrix &B, const T& c=1, const char zflag=0);
    wbMatrix& times(T fac, wbMatrix &X) const {
        X=*this; X*=fac; return X;
    };

    wbMatrix& rowTimes(size_t i, T fac); 
    wbMatrix& colTimes(size_t j, T fac); 

    T prod() const {
       size_t i,n=dim1*dim2; T x=0;

       if (n==0) {
          wblog(FL,"WRN %s() of null object !?",FCT);
          return x;
       }

       for (x=data[0], i=1; i<n; ++i)
       if (x<data[i]) x*=data[i];

       return x;
    };

    T min() const {
       size_t i,n=dim1*dim2; T x=0;

       if (n==0) { wblog(FL,
          "WRN %s() of null object !?",FCT); return x; }
       for (x=data[0], i=1; i<n; ++i) { if (x>data[i]) x=data[i]; }

       return x;
    };

    T max() const { 
       size_t i,n=dim1*dim2; T x=0;

       if (n==0) { wblog(FL,
          "WRN %s() of null object !?",FCT); return x; }
       for (x=data[0], i=1; i<n; ++i) { if (x<data[i]) x=data[i]; }

       return x;
    };

    T aMax() const { 
       T m=0; if (dim1 && dim2) {
          size_t i=0, n=dim1*dim2; T a;
          for (; i<n; ++i) { a=Wb::abs(data[i]); if (m<a) { m=a; }}
       }
       else wblog(FL,"WRN %s() got empty object",FCT);
       return m;
    };

    T aMax(unsigned &i_, unsigned &j_) const { 
       T m=0, a;
       if (!dim1 && !dim2) wblog(FL,"ERR %s() got empty object",FCT);

       size_t i,j,l=0;
       for (i=0; i<dim1; ++i)
       for (j=0; j<dim2; ++j, ++l) { if (m<(a=Wb::abs(data[i]))) {
           m=a; i_=i; j_=j;
       }}
       return m;
    };

    T sum() const { 
        T x=0; 
        size_t i, n=dim1*dim2;
        for (i=0; i<n; ++i) x+=data[i];
        return x;
    };

    T cumsum0prodRec(wbvector<T> &cs, char xflag=0) const {
       wbvector<T> x;
       return recProd(x).cumsum_(cs,xflag);
    };

    T cumsum0prodRec(
       wbvector<T> &cs, unsigned i1, unsigned i2, char xflag=0
     ) const {

       wbvector<T> x(dim1);
       if (i2>=dim2) wblog(FL,
          "ERR %s() index out of bounds (%d/%d)",FCT,i2+1,dim2);

       if (i1<=i2) {
          unsigned i=0, j, n=i2-i1+1; const T* d=data+i1;
          for (; i<dim1; ++i, d+=dim2) { 
             x[i]=d[0];
             for (j=1; j<n; ++j) { x[i]*=d[j]; }
          }
       }

       return x.cumsum_(cs,xflag);
    };

    T norm2() const {
       T x2=0; size_t i, n=dim1*dim2;
       for (i=0; i<n; ++i) x2+=Wb::norm2(data[i]);
       return x2;
    };
    T norm() const { return Wb::sqrt(norm2()); }

    double normReal() const;
    double normImag() const;

    T norm2Cols(long j1, long j2) const { 
       T x2=0; long i,j;

       if (j1<0) { j1+=dim2; }
       if (j2<0) { j2+=dim2; }; if (j1>j2 || !dim1) { return x2; }
       if (j1<0 || j2>=dim2) wblog(FL,
          "ERR %s() index out of bounds (%ld,%ld/%ld)",FCT,j1,j2,dim2);

       for (i=0; i<dim1; ++i)
       for (j=j1; j<=j2; ++j) { x2+=Wb::norm2(data[i*dim2+j]); } 
       return x2;
    };

    T norm2Recs(long i1, long i2) const { 
       T x2=0; long i,j;

       if (i1<0) { i1+=dim1; }
       if (i2<0) { i2+=dim1; }; if (i1>i2 || !dim2) { return x2; }
       if (i1<0 || i2>=dim1) wblog(FL,
          "ERR %s() index out of bounds (%ld,%ld/%ld)",FCT,i1,i2,dim1);

       for (i=i1; i<=i2; ++i)
       for (j=0; j<dim2; ++j) { x2+=Wb::norm2(data[i*dim2+j]); } 
       return x2;
    };

    wbvector<T>& norm2Cols(wbvector<T> &x2_) const { 
       x2_.init(dim2); T* x2=x2_.data; size_t i,j, l=0;
       for (i=0; i<dim1; ++i) 
       for (j=0; j<dim2; ++j, ++l) { x2[j]+=Wb::norm2(data[l]); }
       return x2_;
    };

    wbvector<T> norm2Cols() const {
       wbvector<T> x2; return norm2Cols(x2); 
    };

    wbvector<T> normCols() const {
       wbvector<T> x_; norm2Cols(x_); T *x=x_.data; 
       for (size_t i=0; i<dim2; ++i) { x[i]=Wb::sqrt(x[i]); }
       return x_;
    };

    size_t maxRec(char lex=1, size_t *m=NULL) const;
    size_t maxRec_float(char lex=1, size_t *m=NULL, T xref=-1) const {
       return maxRec(lex,m);
    };

    T recMax(size_t r, size_t *k=NULL) const;
    T colMax(size_t c, size_t *k=NULL) const;

    T colMax(size_t c, unsigned &k) const { T x; 
       size_t K; x=colMax(c,&K); k=K;
       if (size_t(k)!=K) wblog(FL,
          "ERR %s() unsigned out of bounds (%d/%ld)",FCT,k,K);
       return x;
    };

    wbMatrix& NormalizeCol(size_t k);

    bool isOrthoRows(T *x=NULL, char tnorm=0, T eps=1E-14) const;
    bool isOrthoRows(const wbMatrix &B, 
       T *x=NULL, char tnorm=0, T eps=1E-14) const;

    double recMaxA(size_t r, size_t *k=NULL) const;

    size_t skipNanRecs(wbindex &I);

    wbMatrix& ResizeXRecs(size_t r); 
    wbMatrix& Resize(size_t r, size_t c, const T* dx=NULL);
    wbMatrix& resize(size_t m, size_t n, wbMatrix &M) const;

    wbMatrix& Resize2Mult(const WBINDEX &M);

    wbMatrix& transpose(wbMatrix &M) const;
    wbMatrix& Transpose() {
        wbMatrix X; save2(X); return X.transpose(*this);
    };

    wbMatrix& Reshape(size_t r, size_t c);

    template <class T2>
    wbMatrix<T>& AddCols(size_t n, const T2* d0=NULL) {
       Resize(dim1,dim2+n);
       if (d0) Wb::cpyStride(data+dim2-n, d0, n, dim1, dim2);
       return *this;
    };

    template <class T2>
    wbMatrix<T>& AddRows(size_t n, const T2* d0=NULL) {
       Resize(dim1+n,dim2);
       if (d0) MEM_CPY<T>(data+(dim1-n)*dim2, n*dim2, d0);
       return *this;
    };

    bool isEmpty() const { return data==NULL; };
    bool isNormal() const;

    bool isFinite() const {
       return Wb::is_finite(data,dim1*dim2);
    };

    bool isSquare() const { return (dim1==dim2); };
    bool isSquare(size_t d) const { return (dim1==dim2 && d==dim1); };

    bool isScalar() const { return (dim1==1 && dim2==1); };
    bool isVector() const { return (dim1==1 || dim2==1); };

    bool isDiag() const;
    bool isDiagMatrix(T eps=1E-14) const;

    bool isIdentity(double eps, double &maxdiff) const;
    bool isProptoId(T &x, T eps=1E-14) const;

    bool isHConj( 
      const wbMatrix &B, double eps=1E-12, double *xref=NULL
    ) const { return isSym_aux(B, eps, xref, 's'); };

    bool isHConj(
      double eps=1E-12, double *xref=NULL
    ) const { return isSym_aux(*this,eps,xref,'s'); };

    bool isAHerm( 
      const wbMatrix &B, double eps=1E-12, double *xref=NULL
    ) const { return isSym_aux(B, eps, xref, 'a'); };

    bool isAHerm(
      double eps=1E-12, double *xref=NULL
    ) const { return isSym_aux(*this,eps,xref,'a'); };

    bool isComplex() const;

    bool isUnique() const;

    bool isUniqueSorted(char dir=0, char lex=1) const;

    int isSorted( 
       char dir=0, 
       char lex=1,
       size_t m=-1
    ) const;

    bool quickCheckSorted(
        char dir, 
        char lex=1,
        const char *F=NULL, int L=0, 
        size_t n=2 
    ) const;

    bool thisIsInt() const; 

    int getDiff(
    const wbMatrix &B, wbindex &Ia, wbindex *Ib=NULL) const;

    int getDiffSorted(
    const wbMatrix &B, wbindex &Ia, wbindex *Ib=NULL) const;

    size_t length() const {
       return dim1>dim2 ? (dim2 ? dim1 : 0) : (dim1 ? dim2 : 0);
    };

    size_t numel() const { return dim1*dim2; }; 

    void flipSign() {
        size_t i,s=dim1*dim2;
        for (i=0; i<s; ++i) data[i]=-data[i];
    };

    wbMatrix& Symmetrize() {
        size_t i,j,r,s;
        if (dim1!=dim2) wblog(FL,
        "ERR Symmetrize() called with %dx%d matrix !?",dim1,dim2);

        for (j=0; j<dim2; ++j)
        for (i=j+1; i<dim1; ++i) {
            r=i*dim2+j; s=j*dim2+i;
            data[r]=data[s]=0.5*(data[r]+data[s]);
        }

        return *this;
    };

    void set(const T &x) {
       for (size_t s=dim1*dim2, i=0; i<s; ++i) { data[i]=x; }
       isdiag=0;
    };

    wbMatrix& setRand(double fac=1., double shift=0.);
    void setDiagRand(double fac=1., double shift=0.);

    double maxRelDiff(const wbMatrix &M) const;

    template <class T2>
    double normDiff(const wbMatrix<T2> &M) const;

    T normDiff2(const wbMatrix &M, size_t *k=NULL) const;
    T normDiff (const wbMatrix &M, size_t *k=NULL) const {
       return Wb::sqrt(normDiff2(M,k));
    };

    bool sameSize(const wbMatrix B) const { 
       return (dim1==B.dim1 && dim2==B.dim2);
    };

    WBINDEX& toIndex(WBINDEX &I, const WBINDEX *S=NULL) const { 
       wblog(FL,"ERR %s() not defined for type '%s'",FCT,
       TSTR(T)); return I;
    };

    WBIDXMAT& toIndex2D(
      const WBINDEX &S, const wbindex &ic, WBIDXMAT &IJ,
      char pos=1 
    ) const { 
       wblog(FL,"ERR %s() not defined for type '%s'",FCT,
       TSTR(T)); return IJ;
    };

    wbMatrix& SkipTiny_float(T ref __attribute__ ((unused)) =-1) {
       return *this;
    };

    void reset() {
       if (dim1 && dim2) MEM_SET<T>(data,dim1*dim2);
       isdiag=0;
    }

    size_t nnz() const {
       size_t n=0, i=0,s=dim1*dim2;
       for (i=0; i<s; ++i) { if (data[i]!=T(0)) ++n; }
       return n;
    };

    size_t nnz(T eps) const {
       size_t n=0; 
       if (!eps) { n=nnz(); } else {
          size_t i,s=dim1*dim2;
          for (i=0; i<s; ++i) { if (Wb::abs(data[i])>eps) ++n; }
       }
       return n;
    };

    wbMatrix& CAT(C_UINT d12, const wbMatrix** M, C_UINT len);
    wbMatrix& CAT(C_UINT d12, const wbvector< wbMatrix > &M0);
    wbMatrix& CAT(C_UINT d12,       wbvector< wbMatrix const* >  M); 

    wbMatrix& CAT(C_UINT d12, const wbvector< wbvector<T> > &M);
    wbMatrix& CAT(C_UINT d12, const wbvector< wbvector<T> const* > &M);

    wbMatrix& cat(C_UINT d12, const wbMatrix&);
    wbMatrix& cat(C_UINT d12, const wbMatrix&, const wbMatrix&);
    wbMatrix& cat(C_UINT d12, const wbMatrix&, const wbMatrix&, const wbMatrix&);

    wbMatrix& Cat(C_UINT d12,
       const wbMatrix&, const wbMatrix&);
    wbMatrix& Cat(C_UINT d12,
       const wbMatrix&, const wbMatrix&, const wbMatrix&);
    wbMatrix& Cat(C_UINT d12,
       const wbMatrix&, const wbMatrix&, const wbMatrix&, const wbMatrix&);

    wbMatrix<T>& CAT(const unsigned dim,         
        wbvector< wbMatrix<T> const* > M,        
        wbvector< WBINDEX const* > I  
    );

    wbMatrix& Cat(const unsigned dim,
       const wbMatrix &M1, const wbindex &I1,
       const wbMatrix &M2, const wbindex &I2
    ){
        wbvector< wbMatrix<T> const* > M(2);
        wbvector< WBINDEX const* > I(2);
        M.data[0]=&M1; I.data[0]=&I1;
        M.data[1]=&M2; I.data[1]=&I2;
        return CAT(dim,M,I); 
    };

    void split(C_UINT d12, wbvector< wbMatrix* > &M) const;
    void split(C_UINT d12, wbMatrix** M, C_UINT len) const;
    void split(C_UINT d12, wbMatrix &M1, wbMatrix &M2, wbMatrix &M3) const;

    const T& operator() (size_t i, size_t j, const T&x) const {
       if (i<dim1 && j<dim2) return data[i*dim2+j];
       else return x;
    };
    T& operator() (size_t i, size_t j, const T&x) {
       if (i<dim1 && j<dim2) return data[i*dim2+j];
       else return x;
    };

    const T& operator() (size_t i, size_t j) const { return data[i*dim2+j]; };
          T& operator() (size_t i, size_t j)       { return data[i*dim2+j]; };

    const T& el(size_t i, size_t j) const { 
       if (i>=dim1 || j>=dim2) wblog(FL,"ERR %s() "
          "index out of bounds (%d,%d; %dx%d)",FCT,i+1,j+1,dim1,dim2);
       return data[i*dim2+j];
    };
    T& el(size_t i, size_t j) {
       if (i>=dim1 || j>=dim2) wblog(FL,"ERR %s() "
          "index out of bounds (%d,%d; %dx%d)",FCT,i+1,j+1,dim1,dim2);
       return data[i*dim2+j];
    };

    const T& els(long i, long j) const { 
       if (i<0) i+=dim1; 
       if (j<0) j+=dim2; 
       if (i>=dim1 || j>=dim2) wblog(FL,"ERR %s() "
          "index out of bounds (%d,%d; %dx%d)",FCT,i+1,j+1,dim1,dim2);
       return data[i*dim2+j];
    };
    T& els(long i, long j) {
       if (i>=dim1 || j>=dim2) wblog(FL,"ERR %s() "
          "index out of bounds (%d,%d; %dx%d)",FCT,i+1,j+1,dim1,dim2);
       return data[i*dim2+j];
    };

    const T* operator() (size_t i) const { return data+i*dim2; };
          T* operator() (size_t i)       { return data+i*dim2; };

    const T& operator[] (size_t i) const { return data[i]; };
          T& operator[] (size_t i)       { return data[i]; };

    const T* ref(size_t i, size_t j=0) const { 
       if (int(j)<0) j=dim2+j; 
       if (i>=dim1 || j>dim2) wblog(FL,"ERR %s() "
          "index out of bounds (%ld/%ld, %ld/%ld)",FCT,i,dim1,j,dim2);
       return data+(i*dim2+j);
    };
    T* ref(size_t i, size_t j=0) {
       if (int(j)<0) j=dim2+j;
       if (i>=dim1 || j>dim2) wblog(FL,"ERR %s() "
          "index out of bounds (%ld/%ld, %ld/%ld)",FCT,i,dim1,j,dim2);
       return data+(i*dim2+j);
    };

    const T* rec(size_t i) const {
       if (i>=dim1) wblog(FL,"ERR index out of bounds (%d/%d)",i,dim1);
       return data+i*dim2; 
    };

    T* rec(size_t i) {
       if (i>=dim1) wblog(FL,
          "ERR %s() index out of bounds (%d/%d)",FCT,i,dim1);
       return data+i*dim2;
    };

    void getDiag(wbvector<T> &D) const {
       size_t i, n=MIN(dim1,dim2); D.init(n);
       for (i=0; i<n; ++i) D[i]=data[i*dim2+i];
    };
    wbvector<T> getDiag() const {
       wbvector<T> D; getDiag(D); return D;
    };

    void appendRows (size_t n, const T* =NULL);
    void appendRow(const wbvector<T> &v);

    mxArray* toMx (char raw=0) const;  
    mxArray* toMxP(char flag=0) const; 
    mxArray* toMx_base(char raw=0) const;
    mxArray* toMx_base_d(char raw=0) const;

    mxArray* toMx_Struct() const;
    mxArray* toMxP_S() const; 
    mxArray* toMxP_C() const; 

    mxArray* toMxT(char raw=0) const { return toMx(); }

    mxArray* mxCreateStruct(unsigned m, unsigned n) const;
    void add2MxStruct(mxArray *S, unsigned i, char tst=0) const;

    void toMxStruct(mxArray* S, const char *vname) const;
    void mat2mx (mxArray* &a) const { a=toMx(); };
    void mat2mxs(mxArray* a, char field_nr) const;

    void info(const char *istr="ans") const;
    void print(const char *istr="", char mflag=0) const; 
    void Print(const char *istr="", char mflag=0) const; 

    void recPrint(size_t i, const char *istr="", char mflag=0) const;

    void printdata(
    const char *istr, const char *dbl_fmt=" %8.4g", const char *rsep="\n") const;

    wbstring toStr() const { return toStr(""," ","; ",0,""); };

    wbstring toStr(
       const char *fmt, const char *sep=" ", const char *rsep="; ",
       size_t stride=0, const char *sep2=""
    ) const;

    wbstring rec2Str(
       size_t k, const char *fmt="", const char *sep=" ",
       size_t stride=0, const char *sep2=""
    ) const;

    wbstring sizeStr() const { 
       wbstring s(16); 
       size_t l=snprintf(s.data,s.len,"%ldx%ld",dim1,dim2);
       if (l>=s.len) wblog(FL,
          "WRN %s() string out of bounds (%ld/%ld)",FCT,l,s.len);
       return s;
    };

    void put(const char *F, int L,
      const char *vname="ans", const char *ws="caller", const char raw=0
    ) const {
       mxArray *a=toMx(raw);
       mxPutAndDestroy(F_L,a,vname,ws);
    };

    void put(const char *vname="ans", const char *ws="caller",
    const char raw=0) const { put(0,0,vname,ws,raw); }

    void putg(const char *vname="ans", const char raw=0
    ) const { put(vname, "global",raw); return; };

    void putx(const char *vname="ans", const char raw=0
    ) const { put(vname, "caller", raw); return; };

    void put0(const char *vname="ans", const char *ws="caller") const {
       put(vname,ws,'r');
    }

    void getReal(wbMatrix<double> &R) const;
    void getImag(wbMatrix<double> &I) const;
    void set(const wbMatrix<double> &R, const wbMatrix<double> &I);
    wbMatrix& Conj();

    void setRec(size_t r, const T x) {
        T* d=data+r*dim2;
        if (r>=dim1) wblog(FL,"ERR index out of bounds (%d/%d)",r,dim1);
        for (size_t i=0; i<dim2; ++i) d[i]=x;
    };

    void setRec(size_t r, const T x1, const T x2) {
        size_t i=0, m=dim2%2, n=dim2-m; T* d=data+r*dim2;

        if (r>=dim1) wblog(FL,
           "ERR index out of bounds (%d/%d)",r,dim1);

        for (; i<n; i+=2) { d[i]=x1; d[i+1]=x2; }
        if (m) d[i]=x1;
    };

    void recSet(size_t idest, size_t isrc);

    void recSet(size_t, const wbvector<T>&);

    template<class T2> 
    void recSetT(size_t, const wbvector<T2>&);

    void recSet(size_t, const wbMatrix&, size_t);
    void recSet(size_t i, const wbvector<T>&, const wbvector<T>&);

    void recSet (size_t, const T*, T bfac=1); 
    void recSetP(size_t, const T*, size_t n=-1); 
    void recSetP(size_t, const T*, size_t, const T*, size_t);
    void recSetP(size_t,
       const T*, const wbindex &i1, const T*, const wbindex &i2);

    void recSetB(size_t, size_t, size_t D, const T*, 
                 const T* =NULL, const T bfac=1);
    void recSetB(size_t, const wbindex&, size_t D, const T*);

    void recSet(size_t k,
       const wbvector<T> &a, const wbindex &ia,
       const wbvector<T> &b, const wbindex &ib
    ){ return recSet(a.data,ia,b.data,ib); };

    void recSet(size_t k,
       const T* a, const wbindex &ia,
       const T* b, const wbindex &ib
    );

    template<class T2>
    void setCol(size_t k, const wbvector<T2> &v); 

    template<class T2>
    void setCol(size_t k, const T2 *v, T fac=1, size_t stride=1);

    void setCol(size_t k, const T& x);

    template <class T2>
    void setColT(size_t k, const wbvector<T2> &v);
    template <class T2>
    void setColT(size_t k, const T2* v, T fac=1, size_t stride=1);

    void setLastCol(const T& x) { setCol(dim2-1,x); };

    void recAddP(size_t, const T*);
    int  recLess (size_t, size_t, char lex=1) const;
    int  recLessE(size_t, size_t, char lex=1) const;
    int  recEqual(size_t, size_t) const;
    int  recEqual(size_t j1, const T *d) const;

    T recDiff2(size_t i1, size_t i2, size_t n=-1) const;
    T recDiff2P(size_t i1, const T* b, size_t n=-1) const;

    bool recAllEqual(const T *v) { 
       for (size_t i=0; i<dim1; ++i) { if (recCompareP(i,v)) return 0; }
       return 1;
    };

    bool recAllEqual() {
       for (size_t i=1; i<dim1; ++i) { if (recCompare(i,0)) return 0; }
       return 1;
    };

    bool recIsZero(size_t r) {
       T *d=data+r*dim2;
       if (r>=dim1) wblog(FL,"ERR index out of bounds (%d/%d)",r,dim1);
       for (size_t i=0; i<dim2; ++i) if (d[i]!=0) return 0;
       return 1;
    };

    char recCompare(
       size_t, size_t, size_t m=-1, char lex=1, T eps=0) const;
    char recCompare(
       size_t, const wbvector<T> &r2, char lex=1) const;
    char recCompareP(
       size_t, const T*, size_t m=-1, char lex=1, T eps=0) const;

    void set2Rec(size_t i) { recSet(0,i); Resize(1,dim2); };

    wbMatrix& Set2Recs(const WBINDEX &I) {
       wbMatrix X; save2(X);
       return X.getRecs(I,*this);
    };

    wbvector<T>& getRec(size_t j0, wbvector<T> &v, char ref=0) const;
    wbvector<T> getRec(size_t j0) const {
       wbvector<T> v; return getRec(j0,v);
    };

    size_t UnionRecs(const wbMatrix &B, wbindex &Ia);

    bool gotRecOverlap  (const wbMatrix &B) const;
    bool gotRecOverlapSA(const wbMatrix &B) const;

    size_t findValsCol(
       const size_t k, const wbvector<T> &v, wbindex &Ia) const;

    wbMatrix& getRecs(const WBINDEX &, wbMatrix &M) const;
    wbMatrix  getRecs(const WBINDEX &I) const {
       wbMatrix M; getRecs(I,M); 
       return M;
    };

    wbMatrix& getRecs(
       size_t j1, size_t j2, wbMatrix &M, size_t m=0) const;

    wbMatrix  getRecs(
       size_t i1, size_t i2, size_t m=0) const {
       wbMatrix M; getRecs(i1,i2,M,m); 
       return M;
    };

    wbvector<T>& getCol(const size_t c, wbvector<T> &v) const;
    wbvector<T>  getCol(const size_t c) const {
       wbvector<T> v; getCol(c,v); return v;
    };

    wbMatrix getCols(const size_t j1, const size_t j2) const {
       wbMatrix M; getCols(j1,j2,M); 
       return M;
    };

    wbMatrix& getCols(size_t j1, size_t j2, wbMatrix &M) const;

    template <class TI>
    wbMatrix& getCols(size_t n, const TI *I, wbMatrix<T> &M) const;

    wbMatrix& getCols(const wbperm &p, wbMatrix &M) const {
       return getCols(p.len,p.data,M);
    };

    wbMatrix& getCols(
       const wbindex &I, wbMatrix &M, char iflag=0) const {

       if (!iflag) { return getCols(I.len,I.data,M); }
       else {
          wbindex I2; I.invert(dim2,I2);
          return getCols(I2.len,I2.data,M);
       }
    };

    wbMatrix& set2Cols( 
       const size_t j1, const size_t j2);

    wbMatrix& set2Cols(const wbindex &J, char iflag=0) {
       wbMatrix X; this->save2(X); X.getCols(J,*this,iflag);
       return *this;
    };

    wbMatrix& recPermute(const wbperm &P, wbMatrix &M, char iflag=0) const;
    wbMatrix& recPermute(const wbperm &P, char iflag=0) {
       wbMatrix M(*this); 
       return M.recPermute(P,*this,iflag); 
    };

    wbMatrix& colPermute(const wbperm &P, wbMatrix &M, char iflag=0) const;
    wbMatrix& colPermute(const wbperm &P, char iflag=0){
       wbMatrix M; save2(M); return M.colPermute(P,*this,iflag);
    };

    wbMatrix& FlipRecs() { wbperm P(dim1,'r'); return recPermute(P); };
    wbMatrix& FlipCols() { wbperm P(dim2,'r'); return colPermute(P); };

    wbMatrix& cols2Front(
       const WBINDEX &I1, wbMatrix &M) const;
    wbMatrix& cols2End(
       const WBINDEX &I1, wbMatrix &M) const;

    void blockPermute(const wbperm &P, wbMatrix &M) const;
    void BlockPermute(const wbperm &P){
       wbMatrix X; save2(X); X.blockPermute(P,*this);
    };

    wbvector<T>& recProd(wbvector<T> &p) const; 
    T recProd(size_t) const;

    T recSum (size_t) const; 
    wbvector<T>& recSum(wbvector<T> &) const;
    wbvector<T>& recSumA(wbvector<T> &) const;
    wbvector<T> recSum() const { wbvector<T> s; return recSum(s); }; 
    T colSum (size_t) const;

    wbMatrix& repmat(
       size_t m, size_t n, wbMatrix &Q,
       size_t pad1=0, size_t pad2=0
    ) const;

    wbMatrix& Repmat(
       size_t m, size_t n,
       size_t pad1=0, size_t pad2=0
    ){ wbMatrix X(*this); return X.repmat(m,n,*this,pad1,pad2); };

    void ColKron(
       const wbMatrix &Q,
       int d=-1 
    );

    wbMatrix& getBlocks( 
       size_t j1, size_t j2, 
       const size_t D, wbMatrix &M
    ) const;

    wbMatrix& getBlock(
       size_t j1, 
       const size_t D, wbMatrix &M) const;

    wbMatrix getBlock(
       size_t j1, 
       const size_t D
    ) const {
       wbMatrix M; getBlock(j1,D,M); 
       return M;
    };

    template <class T2>
    wbMatrix<T>& getBlock(
       size_t k, size_t D, 
       const wbMatrix<T2>& R, size_t k2, size_t D2, 
       wbMatrix &M 
    ) const;

    wbMatrix& getBlock(
       size_t i, size_t j, 
       size_t n, size_t m, 
       wbMatrix &M
    ) const;

    wbMatrix& blockSum(size_t D, wbMatrix&) const;
    void blockSum(size_t, size_t, const size_t, wbMatrix<double>&) const;
    void blockSum(const wbindex &J, const size_t, wbMatrix<double>&) const;

    void setBlock(size_t k, const size_t D, const T*);
    void SetBlock(size_t i0, size_t j0, const wbMatrix &M);

    WBIDXMAT& toBlockIndex(
       widx_t D, WBIDXMAT &II,
       wbvector< wbMatrix > *QI=NULL
    ) const;

    template <class T2>
    void toBlockIndex(widx_t D,
       wbMatrix<T2> &R, widx_t DR, 
       WBIDXMAT *II=NULL, WBIDXMAT *IS=NULL,
       WBIDXMAT *SS=NULL, wbvector< wbMatrix > *QI=NULL
    ) const;

    void groupRecs(
       wbperm &P, WBINDEX &d,
       size_t m=-1,
       char lex=1,
       wbindex *Ig=NULL,
       wbMatrix<T> *X=NULL
    );

    void groupRecs(wbindex &Ig, widx_t m=-1, char lex=1) {
       wbperm P; wbvector<widx_t> d;
       groupRecs(P,d,m,lex,&Ig);
    };

    void groupRecs(
       groupIndex<widx_t> &IG, size_t m=-1, char lex=1,
       wbindex *Ig=NULL
    ){
       wbperm P; wbvector<widx_t> d;
       groupRecs(P,d,m,lex,Ig);
       IG.initX(P,d);
    };

    void groupRecs(
       groupIndex<widx_t> &IG, size_t m, char lex, const wperm_t **p
    ){
       wbperm P; wbvector<widx_t> d;
       groupRecs(P,d,m,lex);
       IG.initX(P,d,p);
    };

    void groupRecs(
       wbperm &P, WBINDEX &d, wbMatrix &B,
       size_t m=-1, char lex=1, wbindex *Ig=NULL
    ) const { B=*this; B.groupRecs(P,d,m,lex,Ig); };

    void groupRecs(
       wbperm &P, WBINDEX &d,
       const WBINDEX *I, char lex=1
    );

    template <class T2>
    void groupRecs_rdeg(wbperm &P, WBINDEX &D, const wbMatrix<T2> &R,
       WBINDEX *Ib=NULL, 
       WBINDEX *I2=NULL, 
       WBINDEX *Sb=NULL, 
       char iflag=0 
    );

    void groupRecs(
       wbperm &P, WBINDEX &D, size_t nc,
       WBINDEX &Ib, 
       WBINDEX &I2, 
       WBINDEX &Sb, 
       char iflag=0 
    );

    void groupSortedRecs(
       WBINDEX &d, char keepall=0, size_t m=-1, char lex=1);

    void groupSortedRecs(
       WBINDEX &d, const WBINDEX *I);

    void groupSortedRecs(
       WBINDEX &d,
       size_t mc, 
       char lex,
       WBINDEX &Ib, 
       WBINDEX &I2, 
       WBINDEX &Sb  
    );

    void makeUnique() { if (dim1<2) return;
       wbperm P; WBINDEX D;
       groupRecs(P,D);
    };
    void makeUnique(wbperm &P, WBINDEX &D) {
       if (dim1<2) {
          P.init(dim1); D.init(dim1); if (dim1) D[0]=1;
          return;
       }
       groupRecs(P,D);
    };

    void makeUnique(wbindex &I) { 
       if (dim1<2) { I.init(dim1); }
       else {
          size_t i=0, l=0, d; wbperm P; WBINDEX D;
          groupRecs(P,D); I.init(D.len);
          for (; i<D.len; ++i, l+=d) { d=D[i]; I[i]=P[l]; }
       }
    };

    void findUnique(wbindex &I) const {
       if (dim1<2) { I.init(dim1); return; }

       wbMatrix X(*this);
       wbperm P; WBINDEX D; X.groupRecs(P,D);
       I.init(D.numEQ(T(1)));
       for (size_t d, l=0, j=0, i=0; i<D.len; ++i, l+=d) {
          if ((d=D[i])==1) I[j++]=P[l];
       }
    };

    size_t findUniqueRecSorted1(size_t n=-1, T eps=T(1E-12)) const;

    void makeUnique_ig(wbindex &Iu){
       wbperm P; wbvector<widx_t> D; groupRecs(P,D);
       Iu.initGroup(P,D);
    };
    void makeUnique(groupIndex<widx_t> &IG){ groupRecs(IG); };

    wbMatrix& SortRecs(wbperm &P,
       char dir=+1, 
       char lex=1
    );

    wbMatrix& sortRecs(
       wbMatrix &B, wbperm &P, char dir=+1, char lex=1
    ){
       B=*this;
       return B.SortRecs(P,dir,lex);
    };

    wbMatrix& SortRecs(char dir=+1, char lex=1) { 
       wbperm P; return SortRecs(P,dir,lex);
    };

    wbMatrix& sortRecs_float(wbperm &P, char dir=+1) {
       P.init(dim1); return *this;
    };
    wbMatrix& sortRecs_float() { wbperm P; return sortRecs_float(P); };

    char findRecsInSet(const wbMatrix &S, wbvector<int> &I) const;
    int findRecSorted(const T* r, size_t n=-1, char lex=1) const;
    int findRec(const T* r, size_t n=-1, char lex=1) const;

    T *data;

    size_t dim1, dim2;

    char isdiag, isref;

  protected:

    void NEW_DATA(T* d=NULL) { size_t s=dim1*dim2;

        if (data) { WB_DELETE(data); }
        if (s) {
            WB_NEW(data,s);
            MEM_CPY<T>(data,s,d);
        }
        else data=NULL;
    };

    void RENEW(
        const size_t &d1, const size_t &d2, const T* dd=NULL, char iflag=1
    ){
        const size_t s=d1*d2; isdiag=0;

        if (isref) {
           if (d1==0 || d2==0) {
              if (d1 || d2) wblog(FL,
                 "WRN wbMatrix is declared as reference "
                 "(%dx%d => %dx%d; %d)",dim1,dim2,d1,d2,isref
              );

              dim1=d1; dim2=d2; isref=0;

              data=NULL;

              return;
           }
           else {
              wblog(FL,"ERR wbMatrix is declared as reference (%dx%d; %d)",
              d1,d2, isref);
           }
        }

        if (d1!=dim1 || d2!=dim2) {

            dim1=d1; dim2=d2; if (data) {
                if (dd==data) wblog(FL,
                   "ERR init space equals *this (use Resize instead)");
                WB_DELETE(data);
            }
            if (d1==0 || d2==0) return;
            WB_NEW(data,s);
        }

        if (!iflag && !dd) return;
        if (s) MEM_CPY<T>(data,s,dd);
    };

  private:

    bool isSym_aux(
       const wbMatrix &B, double eps, double *xref,
       const char symflag='s'
    ) const;
};

template <class T>
class wbRecs {

    const T* data;
    size_t n,lda; bool asc, lex;

  public:

    wbRecs(const wbMatrix<T> &A, wbperm &P, char dir_=1, char lex_=1)
     : data(A.data), n(A.dim1), lda(A.dim2),
       asc(dir_> 0 ? 1 : 0), lex(lex_> 0 ? 1 : 0)
    {
       P.init(n);
    };

   ~wbRecs() {}; 

    bool operator()(const size_t &pa, const size_t &pb) const {

       if (pa>=n || pb>=n) wblog(FL,
          "ERR %s() index out of bounds (%d,%d; %d)",FCT,pa,pb,n);

       const T *a=data+pa*lda, *b=data+pb*lda;
       if (lex) { size_t i=0;
          if (asc) {
             for (; i<lda; ++i) {
                 if (a[i]<b[i]) return 1;
                 if (a[i]>b[i]) return 0;
             }
          }
          else {
             for (; i<lda; ++i) {
                 if (a[i]<b[i]) return 0;
                 if (a[i]>b[i]) return 1;
             }
          }
       }
       else { size_t i=lda-1;
          if (asc) {
             for (; i<lda; --i) {
                 if (a[i]<b[i]) return 1;
                 if (a[i]>b[i]) return 0;
             }
          }
          else {
             for (; i<lda; --i) {
                 if (a[i]<b[i]) return 0;
                 if (a[i]>b[i]) return 1;
             }
          }
       }
       if (asc) return (pa<pb); 
       else return (pb<pa);
    };
};

template <class T>
void getSortPerm_OMP(
   const wbMatrix<T> &A, wbperm &P, char dir, char lex);

template <class T>
bool wbMatrix<T>::deepEqualP(const wbMatrix<T> &B) const {
    if (dim1!=B.dim1 || dim2!=B.dim2) return 0;

    if (data!=B.data) {
       size_t i=0, n=numel();
       for (; i<n; ++i) { if (data[i]==B.data[i]) continue;
          if ((!data[i] || !B.data[i])) return 0; else
          if (!((*data[i])==(*B.data[i]))) return 0;
       }
    }

    return 1;
};

template <class T> inline
bool wbMatrix<T>::isDiag() const {
   size_t i,j, l=0;

   for (i=0; i<dim1; ++i) 
   for (j=0; j<dim2; ++j, ++l) { if (i!=j && data[l]) return 0; }

   return 1;
};

template <class T>
bool wbMatrix<T>::isDiagMatrix(T eps) const { 

   if (!dim1 || dim1!=dim2) return 0;
   size_t i,j, l=0;

   for (i=0; i<dim1; ++i)  
   for (j=0; j<dim2; ++j, ++l) { if (i!=j && Wb::abs(data[l])>eps) return 0; }

   return 1;
};

template <class T>
bool wbMatrix<T>::isIdentity(double eps, double &maxdiff) const {

   if (!dim1 || dim1!=dim2) return 0;
   size_t i,j,l=0; double d; maxdiff=0;

   for (i=0; i<dim1; ++i)
   for (j=0; j<dim2; ++j, ++l) { 
       d = ABS( data[l] - (i!=j ? 0 : 1) );
       if (maxdiff<d) maxdiff=d;
   }

   return (maxdiff<eps);
};

template <class T>
bool wbMatrix<T>::isProptoId(T &x, T eps) const { 

   if (!dim1 || dim1!=dim2) return 0;
   size_t i,j,l=0;

   for (i=0; i<dim1; ++i)
   for (j=0; j<dim2; ++j, ++l) { 
       if (Wb::abs(data[l] - (i!=j ? 0 : data[0]))>eps) { return 0; }
   }
   x=data[0]; return 1;
};

template <class T> inline
bool wbMatrix<T>::isSym_aux(
   const wbMatrix<T> &B, double eps, double *xref,
   const char symflag
) const {

   const wbMatrix<T> &A=(*this);
   char issame = (&A==&B);
   size_t i,j,k;
   double x;

   if (A.dim1!=B.dim2 || A.dim2!=B.dim1) return 0;
   if (xref) *xref=0;

   if (symflag=='s') {
      if (eps==0) {
         for (i=0; i<dim1; ++i) { k = (issame ? i : 0);
         for (j=k; j<dim2; ++j) {
            if (A(i,j)!=Wb::CONJ(B(j,i))) {
               if (xref) {
                  x = Wb::abs(A(i,j) - Wb::CONJ(B(j,i)));
                  xref[0] = MAX(*xref,x);
               }
               else return 0;
            }
         }}
      }
      else {
         for (i=0; i<dim1; ++i) { k = (issame ? i : 0);
         for (j=k; j<dim2; ++j) {
            x = Wb::abs(A(i,j) - Wb::CONJ(B(j,i)));
            if (x>eps) { if (xref) xref[0]=MAX(*xref,x); else return 0; }
         }}
      }
   }
   else if (symflag=='a') {
      if (eps==0) {
         for (i=0; i<dim1; ++i) { k = (issame ? i : 0);
         for (j=k; j<dim2; ++j) {
            if (A(i,j)!=-Wb::CONJ(B(j,i))) {
               if (xref) {
                  x = Wb::abs(A(i,j) + Wb::CONJ(B(j,i)));
                  xref[0] = MAX(*xref,x);
               }
               else return 0;
            }
         }}
      }
      else {
         for (i=0; i<dim1; ++i) { k = (issame ? i : 0);
         for (j=k; j<dim2; ++j) {
            x = Wb::abs(A(i,j) + Wb::CONJ(B(j,i)));
            if (x>eps) { if (xref) xref[0]=MAX(*xref,x); else return 0; }
         }}
      }
   }
   else wblog(FL,"ERR invalid symflag=%c<%d>",symflag,symflag);

   if (xref && (*xref)>eps) return 0;

   return 1;
}

template <class T> inline
bool wbMatrix<T>::isComplex() const { return 0; };

template<> inline
bool wbMatrix<wbcomplex>::isComplex() const {
   for (size_t s=dim1*dim2, i=0; i<s; ++i)
   if (data[i].i!=0.) return 1;

   return 0;
}

template <class T> inline
bool wbMatrix<T>::isUnique() const {

   wbMatrix<T> Q(*this);
   wbperm P;

   if (!dim2) { return (dim1<=1); }

   Q.SortRecs(P);
   return Q.isUniqueSorted();
};

template <class T> inline
bool wbMatrix<T>::isUniqueSorted(char dir, char lex) const {

   if (dim1<=1 || !dim2) return 1;
   size_t i=1; char c=0;

   if (!dir) { if (dim1<=2) return 1;
      c=recCompare(1,0,-1,lex); if (!c) { return 0; }
      i=2;
   }
   else if (abs(int(dir))>1) { 
      if (dir=='A') c=+1; else 
      if (dir=='D') c=-1; else
      wblog(FL,"ERR %s() invalid c='%c'<%d> !?",FCT,c,c);
   }
   else { c=dir; } 

   for (; i<dim1; ++i) {
      if (recCompare(i,i-1,-1,lex)!=c) return 0;
   }

   return 1;
};

template <class T> inline
int wbMatrix<T>::isSorted(char dir, char lex, size_t m) const {

   if (!dim1 || !dim2) { return -1; }
   if (dim1==1) { return 3; } 

   size_t j=1, ndeg=0; char c=0;

   if (!dir) {
      for (; j<dim1; ++j, ++ndeg) {
      if ((c=recCompare(j,j-1,m,lex))) { dir=c; ++j; break; }}
   }
   else if (abs(int(dir))>1) {
      wblog(FL,"ERR %s() invalid direction s=%d",FCT,dir);
   }

   for (; j<dim1; ++j) {
      if ((c=recCompare(j,j-1,m,lex))) { if (c!=dir) return 0; }
      else { ++ndeg; }
   }

   return (1 | (ndeg ? 0 : 2) | (c<0 ? 4 : 0));
};

template <class T> inline
wbMatrix<T>& wbMatrix<T>::ResizeXRecs(size_t nr) {

    no_ref_(FLF);

    if (nr) { 
       if (nr<=dim1 && nr/double(dim1)>0.10) dim1=nr; 
       else Resize(nr,dim2);
    }
    else init(0,dim2);

    return *this;
};

template <class T> inline
wbMatrix<T>& wbMatrix<T>::resize(
   size_t m, size_t n, wbMatrix<T> &M
 ) const {

   if (this==&M) return M.Resize(m,n);

   M.init(m,n);
   m=MIN(m,dim1); n=MIN(n,dim2); if (m==0 || n==0) return M;

   for (size_t i=0; i<m; ++i)
   MEM_CPY<T>(M.data+i*M.dim2, n, data+i*dim2);

   return M;
};

template <class T> inline
wbMatrix<T>& wbMatrix<T>::Resize(size_t nr, size_t nc, const T* dx) {

   size_t d1,d2,s=nr*nc; T *d0=data;
   no_ref_(FLF);

   d1=MIN(dim1,nr); d2=MIN(dim2,nc);

   if (dim1==nr && dim2==nc) return *this; 
   if (d1==0 || d2==0) { RENEW(nr,nc); return *this; }

   WB_NEW(data,s);

   if (nc==dim2) {
      MEM_CPY<T>(data,d1*dim2,d0);
      if (nr>dim1) {
         if (dx)
              MEM_CPY<T>(data+dim1*dim2,(nr-dim1)*dim2,dx);
         else MEM_SET<T>(data+dim1*dim2,(nr-dim1)*dim2);
      }
   }
   else {
      size_t i; T *p=data, *p0=d0;
      for (i=0; i<d1; ++i, p+=nc, p0+=dim2) {
         MEM_CPY<T>(p,d2,p0); if (d2<nc) {
         MEM_SET<T>(p+d2,nc-d2); }
      }
      if (nr>dim1) {
         if (dx) wblog(FL,"ERR %s() "
            "allows no reference data if dim2=%d->%d",FCT,dim2,nc);
         MEM_SET<T>(data+dim1*dim2,(nr-dim1)*dim2);
      }
   }

   dim1=nr; dim2=nc;
   WB_DELETE(d0);

   return *this;
};

template <class T> inline
wbMatrix<T>& wbMatrix<T>::Resize2Mult(const WBINDEX &M) {

   size_t m=M.sum();

   if (M.len!=dim1) wblog(FL,
      "ERR %s() size mismatch (%d/%dx%d)",FCT,M.len,dim1,dim2);
   if (!m || !dim2) { return init(m,dim2); } else
   if (int(m)<0) wblog(FL,"ERR %s() got m=%d !?",FCT,m);

   size_t i,j;
   wbMatrix<T> X(m,dim2);
   const T *x0=data; T *x=X.data;

   for (i=0; i<M.len; ++i, x0+=dim2) { m=M.data[i];
   for (j=0; j<m; ++j, x+=dim2) MEM_CPY<T>(x, dim2, x0); } 

   return X.save2(*this); 
};

template <class T> inline
wbMatrix<T>& wbMatrix<T>::unRef() {

    if (!isref) return *this;
    if (!dim1 || !dim2) {
       if (data) wblog(FL,"ERR data=%lX (%dx%d) !?",data,dim1,dim2);
       isref=0; return *this;
    }

    size_t s=dim1*dim2;
    T const* const d0=data;

    WB_NEW(data,s);

    MEM_CPY<T>(data,s,d0);
    isref=0;

    return *this;
};

template <class T> inline
wbMatrix<T>& wbMatrix<T>::transpose(wbMatrix<T> &M) const {

    if (this==&M) return M.Transpose();

    size_t i,j,k=0; M.init(dim2,dim1);

    if (isdiag) {
       M.isdiag=1; k=MIN(dim1,dim2);
       for (i=0; i<k; ++i) M.data[i*dim1+i]=data[i*dim2+i];
       return M;
    }

    for (j=0; j<dim2; ++j)
    for (i=0; i<dim1; ++i) M.data[k++]=data[i*dim2+j]; 

    return M;
};

template <class T> inline
wbMatrix<T>& wbMatrix<T>::Reshape(size_t nr, size_t nc) {

    if (dim1!=nr || dim2!=nc) {
       if (dim1*dim2!=nr*nc) wblog(FL,
          "ERR %s() size not preserved (%dx%d => %dx%d)",
           FCT,dim1,dim2,nr,nc
       );
       dim1=nr; dim2=nc;
    }
    return *this;
};

template <class T> inline
wbMatrix<T>& wbMatrix<T>::add(
    const wbMatrix<T> &M,
    const T& c,
    const char zflag
){
    size_t i,s=dim1*dim2;

    if (data==NULL && zflag) {
       RENEW(M.dim1, M.dim2, M.data); (*this)*=c;
       isdiag=M.isdiag;
       return *this;
    }

    if (dim1!=M.dim1 || dim2!=M.dim2)
    wbdie(FL,"Dimension mismatch.");

    if (isdiag) if (!M.isdiag) isdiag=0;

    if (c==T(+1)) { for (i=0; i<s; ++i) data[i]+=M.data[i]; } else
    if (c==T( 0)) { return *this; } else
    if (c==T(-1)) { for (i=0; i<s; ++i) data[i]-=M.data[i]; }
    else          { for (i=0; i<s; ++i) data[i]+=(c * M.data[i]); }

    return *this;
}

template <class T> inline
wbMatrix<T>& wbMatrix<T>::minus(
    const wbMatrix<T> &M,
    const T& c,
    const char zflag
){
    size_t i,s=dim1*dim2;

    if (data==NULL && zflag) {
       RENEW(M.dim1, M.dim2, M.data); (*this)*=(-c);
       isdiag=M.isdiag;
       return *this;
    }

    if (dim1!=M.dim1 || dim2!=M.dim2)
    wbdie(FL,"Dimension mismatch.");

    if (isdiag) if (!M.isdiag) isdiag=0;

    if (c==T(+1)) { for (i=0; i<s; ++i) data[i]-=M.data[i]; } else
    if (c==T( 0)) { return *this; } else
    if (c==T(-1)) { for (i=0; i<s; ++i) data[i]+=M.data[i]; }
    else          { for (i=0; i<s; ++i) data[i]-=(c * M.data[i]); }

    return *this;
}

template <class T> inline
wbMatrix<T>& wbMatrix<T>::cat(const unsigned dim,
    const wbMatrix<T> &M2
){
    const wbMatrix<T>* M0[] = {this, &M2};
    wbvector< wbMatrix<T> const* > M; M.init(2,M0);
    return CAT(dim,M);
};

template <class T> inline
wbMatrix<T>& wbMatrix<T>::cat(const unsigned dim,
    const wbMatrix<T> &M2, const wbMatrix<T> &M3
){
    const wbMatrix<T>* M0[] = {this, &M2, &M3};
    wbvector< wbMatrix<T> const* > M; M.init(3,M0);
    return CAT(dim,M);
};

template <class T> inline
wbMatrix<T>& wbMatrix<T>::cat(const unsigned dim,
    const wbMatrix<T> &M2, const wbMatrix<T> &M3, const wbMatrix<T> &M4
){
    const wbMatrix<T>* M0[] = {this, &M2, &M3, &M4};
    wbvector< wbMatrix<T> const* > M; M.init(4,M0);
    return CAT(dim,M);
};

template <class T> inline
wbMatrix<T>& wbMatrix<T>::Cat(const unsigned dim,
    const wbMatrix<T> &M1, const wbMatrix<T> &M2
){
    wbMatrix<T> const* M0[] = {&M1, &M2};
    wbvector< wbMatrix<T> const* > M; M.init(2,M0);
    return CAT(dim,M);
};

template <class T> inline
wbMatrix<T>& wbMatrix<T>::Cat(const unsigned dim,
    const wbMatrix<T> &M1, const wbMatrix<T> &M2,
    const wbMatrix<T> &M3
){
    wbMatrix<T> const* M0[] = {&M1, &M2, &M3};
    wbvector< wbMatrix<T> const* > M; M.init(3,M0);
    return CAT(dim,M);
};

template <class T> inline
wbMatrix<T>& wbMatrix<T>::Cat(const unsigned dim,
    const wbMatrix<T> &M1, const wbMatrix<T> &M2,
    const wbMatrix<T> &M3, const wbMatrix<T> &M4
){
    wbMatrix<T> const* M0[] = {&M1, &M2, &M3, &M4};
    wbvector< wbMatrix<T> const* > M; M.init(4,M0);
    return CAT(dim,M);
};

template <class T> inline
wbMatrix<T>& wbMatrix<T>::CAT(const unsigned dim,
    wbMatrix<T> const* M[], const unsigned len
){
    wbvector< wbMatrix<T> const* > MM;
    MM.init(len,M);

    return CAT(dim,MM);
};

template <class T> inline
wbMatrix<T>& wbMatrix<T>::CAT(const unsigned dim,
    const wbvector< wbMatrix<T> > &M0
){
    wbvector< wbMatrix<T> const* > M(M0.len);
    for (unsigned i=0; i<M0.len; ++i) M[i]=&M0[i];

    return CAT(dim,M);
};

template <class T> inline
wbMatrix<T>& wbMatrix<T>::CAT(const unsigned dim,
    const wbvector< wbvector<T> > &V0
){
    wbvector< wbMatrix<T> > M0(V0.len); 
    wbvector< wbMatrix<T> const* > M(V0.len); 

    for (size_t i=0; i<V0.len; ++i) {
       if (dim==1)
            M0[i].init2ref(1,V0[i].len,V0[i].data);
       else M0[i].init2ref(V0[i].len,1,V0[i].data);

       M[i]=&M0[i];
    }

    return CAT(dim,M);
};

template <class T> inline
wbMatrix<T>& wbMatrix<T>::CAT(const unsigned dim,
    const wbvector< wbvector<T> const* > &V0
){
    wbvector< wbMatrix<T> > M0(V0.len);
    wbvector< wbMatrix<T> const* > M(V0.len);

    for (size_t i=0; i<V0.len; ++i) {
       if (dim==1)
            M0[i].init2ref(1,V0[i]->len,V0[i]->data);
       else M0[i].init2ref(V0[i]->len,1,V0[i]->data);

       M[i]=&M0[i];
    }

    return CAT(dim,M);
};

template <class T>
wbMatrix<T>& wbMatrix<T>::CAT(const unsigned dim,
    wbvector< wbMatrix<T> const* > M 
){
    wbMatrix<T> X;
    size_t i,k;

    if (isdiag) { wblog(FL,"WRN unsetting isdiag"); isdiag=0; }

    if (M.len==0) { this->init(); return *this; }

    for (k=i=0; i<M.len; ++i) {
       if (M[i]==this) {
          if (!(k++)) { this->save2(X); } 
          M[i]=&X;
       }
    }

    for (k=i=0; i<M.len; ++i) {
       if (M[i] && (M[i]->dim1 || M[i]->dim2)) {
          if (k<i) { M[k]=M[i]; }; ++k;
          if (M[i]->isdiag) wblog(FL,
             "ERR %s() got diagonal matrix (i=%d)",FCT, i+1
          );
       }
    }
    if (k<M.len) {
       if (k) M.len=k;
       else { this->init(); return *this; }
    }

    if (dim==1) {
       size_t m,d1,D1,d2=0;

       for (D1=k=0; k<M.len; ++k) { D1+=(M[k]->dim1);
          if (k==0) d2=(M[k]->dim2); else
          if (d2!=(M[k]->dim2)) wblog(FL,
             "ERR dimension mismatch (%d/%d)",d2,M[k]->dim2
          );
       }

       this->init(D1,d2);

       for (m=k=0; k<M.len; ++k) { d1=(M[k]->dim1);
          MEM_CPY<T>(data+d2*m, d1*d2, M[k]->data); 
          m+=d1;
       }
    }
    else if (dim==2) {
       size_t i,m,d1=0,d2,D2;
       T *p, *p2;

       for (D2=k=0; k<M.len; ++k) { D2+=(M[k]->dim2);
          if (k==0) d1=(M[k]->dim1); else
          if (d1!=(M[k]->dim1)) wblog(FL,
             "ERR dimension mismatch (%d/%d)",d1,M[k]->dim1
          );
       }

       this->init(d1,D2);

       for (m=k=0; k<M.len; ++k) { d2=(M[k]->dim2);
          p=data+m; p2=(M[k]->data);

          for (i=0; i<d1; ++i)
          MEM_CPY<T>(p+D2*i, d2, p2+d2*i); 

          m+=d2;
       }
    }
    else wblog(FL,"ERR invalid cat dim=%d",dim);

    return *this;
};

template <class T>
wbMatrix<T>& wbMatrix<T>::CAT( 
    const unsigned dim,
    wbvector< wbMatrix<T> const* > M,        
    wbvector< WBINDEX const* > I  
){
    wbMatrix<T> X;
    size_t i,k;

    if (M.len!=I.len) wblog(FL,
       "ERR %s() size mismatch %d/%d",FCT,M.len,I.len);
    if (!M.len) { this->init(); return *this; }

    for (k=i=0; i<M.len; ++i) {
       if (M[i]==this) {
          if (!(k++)) { this->save2(X); } 
          M[i]=&X;
       }
    }

    for (k=i=0; i<M.len; ++i) {
       if (M[i] && (M[i]->dim1 || M[i]->dim2 || I[i]->len)) {
          if (k<i) { M[k]=M[i]; I[k]=I[i]; }; ++k;
          if (M[i]->isdiag) wblog(FL,
             "ERR %s() got diagonal matrix (i=%d)",FCT, i+1
          );
       }
    }
    if (k<M.len) {
       if (k) { M.len=k; I.len=k; }
       else { this->init(); return *this; }
    }

    if (isdiag) isdiag=0;

    if (dim==1) {
       size_t d1,d2=0, D1=0;

       for (k=0; k<M.len; ++k) { D1+=(I[k]->len);
          if (k==0) d2=(M[k]->dim2); else
          if (d2!=(M[k]->dim2)) wblog(FL,
             "ERR dimension mismatch (vcat: dim2=%d/%d)",d2,M[k]->dim2
          );
       }

       this->init(D1,d2);
       T *x=data;

       for (k=0; k<M.len; ++k) {
          const T* x0=M[k]->data;
          const widx_t *Ik=I[k]->data; D1=M[k]->dim1; d1=I[k]->len;

          for (i=0; i<d1; ++i, x+=d2) {
              if (Ik[i]>=D1) wblog(FL,
                 "ERR %s() index out of bounds (%d/%d)",FCT,Ik[i]+1,D1);
              MEM_CPY<T>(x, d2, x0+Ik[i]*d2); 
          }
       }
    }
    else if (dim==2) {
       size_t i,d1=0,d2,D1,D2=0;

       for (k=0; k<M.len; ++k) { D2+=(d2=M[k]->dim2);
          if (k==0) d1=(I[k]->len); else
          if (d1!=(I[k]->len)) wblog(FL,
             "ERR dimension mismatch (%d/%d)",d1,I[k]->len
          );
       }

       this->init(d1,D2);
       T *x=data;

       for (k=0; k<M.len; ++k) {
          const T* x0=M[k]->data;
          const widx_t *Ik=I[k]->data; D1=M[k]->dim1; d2=(M[k]->dim2);

          for (i=0; i<d1; ++i) {
             if (Ik[i]>=D1) wblog(FL,
                "ERR %s() index out of bounds (%d/%d)",FCT,Ik[i]+1,D1);
             MEM_CPY<T>(x+i*D2, d2, x0+Ik[i]*d2); 
          }; x+=d2;
       }
    }
    else wblog(FL,"ERR invalid cat dim=%d",dim);

    return *this;
};

template <class T> inline
void wbMatrix<T>::split(
    const unsigned dim,
    wbvector< wbMatrix<T>* > &M) const {

    split(dim,M.data,M.len);
    return;
}

template <class T> inline
void wbMatrix<T>::split(
    const unsigned dim,
    wbMatrix<T> &M1,
    wbMatrix<T> &M2,
    wbMatrix<T> &M3
) const {

    wbMatrix<T>* M[] = {&M1,&M2,&M3};
    split(dim,M,3);
    return;
}

template <class T>
void wbMatrix<T>::split(
    const unsigned dim, wbMatrix<T>** M, const unsigned len
)const {

    unsigned k;

    for (k=0; k<len; ++k) {
       if (M[k]==this)
       wbdie(FL, "This matrix itself is included in split list!");
    }

    if (dim==1) {
       size_t m,d1,D1;

       for (D1=k=0; k<len; ++k) D1+=(*(M[k])).dim1;

       if (D1!=dim1)
       wbdie(FL, "Dimension mismatch.");

       for (k=0; k<len; ++k) {
          if (M[k]->dim2!=dim2) M[k]->init(M[k]->dim1, dim2);
       }

       for (m=k=0; k<len; ++k) { d1=M[k]->dim1;
          MEM_CPY<T>(M[k]->data, d1*dim2, data+dim2*m); 
          m+=d1;
       }
    }
    else if (dim==2) {
       size_t i,m,d2,D2;
       T *p, *p2;

       for (D2=k=0; k<len; ++k) D2+=M[k]->dim2;

       if (D2!=dim2)
       wblog(FL,"ERR dimension mismatch (%d/%d)",D2,dim2);

       for (k=0; k<len; ++k) {
          if (M[k]->dim1!=dim1) M[k]->init(dim1,M[k]->dim2);
       }

       for (m=k=0; k<len; ++k) {
          d2=M[k]->dim2; p=data+m; p2=M[k]->data;

          for (i=0; i<dim1; ++i)
          MEM_CPY<T>(p2+d2*i, d2, p+D2*i); 

          m+=d2;
       }
    }
    else wblog(FL,"ERR invalid cat index (%d).", dim);
};

template <class T> inline
int wbMatrix<T>::recLess(size_t i1, size_t i2, char lex) const {

    size_t i;
    T *r1=data+i1*dim2, *r2=data+i2*dim2;

    if (i1>=dim1 || i2>=dim1) wblog(FL,
       "ERR %s() index out of bounds (%d,%d/%d)",FCT,i1,i2,dim1);

    if (lex>0) {
       for (i=0; i<dim2; ++i)
       if (r1[i]<r2[i]) return 1; else
       if (r1[i]>r2[i]) return 0;
    }
    else {
       for (i=dim2-1; i<dim2; --i)
       if (r1[i]<r2[i]) return 1; else
       if (r1[i]>r2[i]) return 0;
    }

    return 0;
}

template <class T> inline
int wbMatrix<T>::recLessE(size_t i1, size_t i2, char lex) const {

    size_t i;
    T *r1=data+i1*dim2, *r2=data+i2*dim2;

    if (i1>=dim1 || i2>=dim1) wblog(FL,
       "ERR %s() index out of bounds (%d,%d/%d)",FCT,i1,i2,dim1);

    if (lex) {
       for (i=0; i<dim2; ++i)
       if (r1[i]<r2[i]) return 1; else
       if (r1[i]>r2[i]) return 0;
    }
    else {
       for (i=dim2-1; i<dim2; --i)
       if (r1[i]<r2[i]) return 1; else
       if (r1[i]>r2[i]) return 0;
    }

    return 1;
}

template <class T> inline
int wbMatrix<T>::recEqual(size_t i1, const T *d) const {
    if (i1>=dim1) wblog(FL,
       "ERR %s() index out of bounds (%d/%d)",FCT,i1,dim1);
    if (data==d) return 1;
    else return (!memcmp(data+i1*dim2, d, dim2*sizeof(T)));
};

template <class T> inline
int wbMatrix<T>::recEqual(size_t i1, size_t i2) const {
    if (i1>=dim1 || i2>=dim1) wblog(FL,
       "ERR %s() index out of bounds (%d,%d/%d)",FCT,i1,i2,dim1);
    if (i1==i2) return 1;
    else return (!memcmp(data+i1*dim2, data+i2*dim2, dim2*sizeof(T)));
};

template <class T> inline
T wbMatrix<T>::recDiff2(size_t i1, size_t i2, size_t n) const {
    if (i1>=dim1 || i2>=dim1) wblog(FL,
       "ERR %s() index out of bounds (%d,%d/%d)",i1,i2,dim1);
    if (i1==i2) return 0;
    return Wb::rangeNormDiff2( 
       data+i1*dim2, data+i2*dim2, int(n)<0 ? dim2 : n
    );
};

template <class T> inline
T wbMatrix<T>::recDiff2P(size_t i1, const T* b, size_t n) const {
    if (i1>=dim1) wblog(FL,
       "ERR %s() index out of bounds (%d/%d)",i1,dim1);
    return Wb::rangeNormDiff2( 
       data+i1*dim2, b, int(n)<0 ? dim2 : n
    );
};

template <class T> inline
char wbMatrix<T>::recCompare(
   size_t i1, size_t i2,
   size_t m, 
   char lex,   
   T eps       
 ) const {

   if (i1>=dim1 || i2>=dim1) wblog(FL,
      "ERR %s() index out of bounds (%d,%d/%d)",FCT,i1,i2,dim1);
   if (!dim2) wblog(FL,
      "ERR %s() got nothing to compare (dim2=%d)",FCT,dim2);

   if (int(m)<0) m=dim2; else
   if (!m || m>dim2) wblog(FL,"ERR size out of bounds (%d/%d)",m,dim2);

   i1*=dim2; i2*=dim2;
   if (lex<=0) {
      i1+=(dim2-m);
      i2+=(dim2-m);
   }

   if (eps<=T(0))
        return Wb::recCompare(data+i1, data+i2, m, lex);
   else return Wb::recCompare(data+i1, data+i2, m, lex, eps);
};

template <class T> inline
char wbMatrix<T>::recCompare(
   size_t i1, const wbvector<T> &v, char lex
 ) const {

   if (v.len!=dim2) wblog(FL,
      "ERR %s() size mismatch (%d,%d)",FCT,v.len,dim2);
   return Wb::recCompare(data+i1*dim2, v.data, dim2, lex);
};

template <class T> inline
char wbMatrix<T>::recCompareP(
  size_t i1, const T* v,
  size_t m,   
  char lex,   
  T eps       
) const {

    if (i1>=dim1) wblog(FL,
       "ERR %s() index out of bounds (%d/%d)",FCT,i1,dim1);
   if (!dim2) wblog(FL,
      "ERR %s() got nothing to compare (dim2=%d)",FCT,dim2);

    if (int(m)<0) m=dim2; else
    if (!m || m>dim2) wblog(FL,"ERR size out of bounds (%d/%d)",m,dim2);

    i1*=dim2;
    if (lex<=0) i1+=(dim2-m); 

    return Wb::recCompare(data+i1, v, m, lex, eps);
};

template <class T> inline
wbvector<T>& wbMatrix<T>::recProd(wbvector<T> &pp) const {
   size_t i=0, j; T *d=data, p;

   pp.init(dim1);

   if (!dim2) wblog(FL,"ERR %s() records of length 0",FCT);

   for (; i<dim1; ++i, d+=dim2) { 
   for (p=d[0], j=1; j<dim2; ++j) { p*=d[j]; }; pp[i]=p; }

   return pp;
};

template <class T> inline
T wbMatrix<T>::recProd(size_t i) const {

   if (i>=dim1 || dim2==0) { wblog(FL,
      "ERR %s() got empty record length (%dx%d; i)",FCT,dim1,dim2,i);
      return 0;
   }

   const T *d=data+i*dim2; T p=d[0];
   for (size_t j=1; j<dim2; ++j) { p*=d[j]; } 
   return p;
};

template <class T>
inline T wbMatrix<T>::colSum(size_t c) const {
   T *d=data+c, x=0;

   if (c>=dim2) wblog(FL,"ERR index out of bounds (%d/%d)", c, dim2);
   if (dim1==0) wblog(FL,
   "WRN Sum over column of length zero (%d/%d).", c, dim2);

   for (size_t i=0; i<dim1; ++i, d+=dim2) x+=(*d);
   return x;
}

template <class T> inline
T wbMatrix<T>::recSum(size_t i) const {

   if (i>=dim1) wblog(FL,
      "ERR %s() index out of bounds (%d,%dx%d)",FCT,i+1,dim1,dim2);
   if (!dim2) wblog(FL,
      "ERR %s() got empty object (%d,%dx%d)",i,dim1,dim2);

   return Wb::addRange(rec(i),dim2);
};

template <class T> inline
wbvector<T>& wbMatrix<T>::recSum(wbvector<T> &s) const {

   s.init(dim1);
   if (!dim2) { if (dim1) wblog(FL,
      "WRN %s() got empty object (%dx%d)",FCT,dim1,dim2);
      return s;
   }

   size_t i=0; const T* d=data;
   for (; i<dim1; ++i, d+=dim2) { s[i]=Wb::addRange(d,dim2); }

   return s;
};

template <class T> inline
wbvector<T>& wbMatrix<T>::recSumA(wbvector<T> &s) const {

   s.init(dim1);
   if (!dim2) { if (dim1) wblog(FL,
      "WRN %s() got empty object (%dx%d)",FCT,dim1,dim2);
      return s;
   }

   size_t i=0; const T* d=data; double q;
   for (; i<dim1; ++i, d+=dim2) {
      q=Wb::addRange2(d,dim2);
      s[i]=std::sqrt(q);
   }

   return s;
};

template <class T>
wbMatrix<T>& wbMatrix<T>::setRand(double fac, double shift) {
   size_t i,s=dim1*dim2;
   static char first_call=1;

   if (first_call) { wb_srand(); first_call=0; } 
   fac/=(double)RAND_MAX; isdiag=0;

   if (shift==0.)
        for (i=0; i<s; ++i) data[i]=(T)(fac*rand());
   else for (i=0; i<s; ++i) data[i]=(T)(fac*rand()+shift);

   return *this;
}

template<>
wbMatrix<wbcomplex>& wbMatrix<wbcomplex>::setRand(double fac, double shift) {
   size_t i,s=dim1*dim2;
   static char first_call=1;

   if (first_call) { wb_srand(); first_call=0; } 
   fac/=(double)RAND_MAX; isdiag=0;

   if (shift==0.) {
      for (i=0; i<s; ++i)
      data[i].set(fac*rand(), fac*rand());
   }
   else {
      for (i=0; i<s; ++i)
      data[i].set(fac*rand()+shift, fac*rand()+shift);
   }

   return *this;
}

template<>
wbMatrix<wbcomplex>& wbMatrix<wbcomplex>::Symmetrize() {
    size_t i,j,r,s;
    if (dim1!=dim2) wblog(FL,
    "ERR Symmetrize() called with %dx%d matrix !?",dim1,dim2);

    for (j=0; j<dim2; ++j)
    for (i=j+1; i<dim1; ++i) { r=i*dim2+j; s=j*dim2+i;
        data[r]=data[s]=wbcomplex(
            0.5*(data[r].r+data[s].r),
            0.5*(data[r].i-data[s].i)
        );
        data[s].Conj();
    }

    return *this;
};

template <class T>
inline void wbMatrix<T>::setDiagRand(double fac, double shift) {
   size_t i, s=MIN(dim1,dim2);
   static char first_call=1;

   if (data==NULL) return;

   if (first_call) { wb_srand(); first_call=0; } 
   fac/=(double)RAND_MAX;

   MEM_SET<T>(data,dim1*dim2);

   if (shift==0.)
        for (i=0; i<s; ++i) data[i*dim2+i]=(T)(fac*rand());
   else for (i=0; i<s; ++i) data[i*dim2+i]=(T)(fac*rand()+shift);

   isdiag=1;
}

template <class T>
mxArray* wbMatrix<T>::toMx_base_d(char raw) const {   
   size_t d1=(raw ? dim2:dim1), d2=(raw ? dim1:dim2);
   return Mx::Array<double>(d1,d2,!raw).copyFromTR(data);
}; 

template <class T>
mxArray* wbMatrix<T>::toMx_base(char raw) const {   
   size_t d1=(raw ? dim2:dim1), d2=(raw ? dim1:dim2);
   return Mx::Array<double>(d1,d2,!raw).copyFromTR(data);
}; 

template <> 
mxArray* wbMatrix<wbcomplex>::toMx_base(char raw) const {
   size_t d1=(raw ? dim2:dim1), d2=(raw ? dim1:dim2);
   return Mx::Array<wbcomplex>(d1,d2,!raw).copyFromTR(data);
}; 

template <>
mxArray* wbMatrix<char>::toMx_base(char raw) const {
   size_t d1=(raw ? dim2:dim1), d2=(raw ? dim1:dim2);
   return Mx::Array<mxChar>(d1,d2,!raw).copyFromTR(data);
}; 

template <> 
mxArray* wbMatrix<long>::toMx_base(char raw) const {
   size_t d1=(raw ? dim2:dim1), d2=(raw ? dim1:dim2);
   if (double(aMax())<1E15)  
        return Mx::Array<double >(d1,d2,!raw).copyFromTR(data);
   else return Mx::Array<int64_T>(d1,d2,!raw).copyFromTR(data);
};

template <> 
mxArray* wbMatrix<unsigned long>::toMx_base(char raw) const {
   size_t d1=(raw ? dim2:dim1), d2=(raw ? dim1:dim2);
   if (double(aMax())<1E15)  
        return Mx::Array<double >(d1,d2,!raw).copyFromTR(data);
   else return Mx::Array<uint64_T>(d1,d2,!raw).copyFromTR(data);
};

template <class T> inline
mxArray* wbMatrix<T>::toMx(char raw) const {
   if (raw) wblog(FL,
      "WRN raw flag will be ignored for typeid=%s",TSTR(T));
   return toMx_Struct();
};

template <class T> inline
mxArray* wbMatrix<T>::toMxP(char flag) const {

   if (!flag || flag=='s' || flag=='S') { return toMxP_S(); }
   if (         flag=='c' || flag=='C') { return toMxP_C(); }

   wblog(FL,"ERR %s() invalid flag (%s)",FCT,cSTR(flag));
   return 0;
};

template <> inline
mxArray* wbMatrix<unsigned>::toMx(char raw) const {
   return toMx_base(raw); 
};
template <> inline
mxArray* wbMatrix<unsigned>::toMxT(char raw) const {
   size_t d1=(raw ? dim2:dim1), d2=(raw ? dim1:dim2);
   return Mx::Array<uint32_T>(d1,d2,!raw).copyFromTR(data);
};

template <> inline
mxArray* wbMatrix<int>::toMx(char raw) const {
   return toMx_base(raw); 
};
template <> inline
mxArray* wbMatrix<int>::toMxT(char raw) const {
   size_t d1=(raw ? dim2:dim1), d2=(raw ? dim1:dim2);
   return Mx::Array<int32_T>(d1,d2,!raw).copyFromTR(data);
};

template <> inline
mxArray* wbMatrix<unsigned long>::toMx(char raw) const {
   return toMx_base(raw); 
};
template <> inline
mxArray* wbMatrix<unsigned long>::toMxT(char raw) const {
   size_t d1=(raw ? dim2:dim1), d2=(raw ? dim1:dim2);
   return Mx::Array<uint64_T>(d1,d2,!raw).copyFromTR(data);
};

template <> inline
mxArray* wbMatrix<long>::toMx(char raw) const {
   return toMx_base(raw); 
};
template <> inline
mxArray* wbMatrix<long>::toMxT(char raw) const {
   size_t d1=(raw ? dim2:dim1), d2=(raw ? dim1:dim2);
   return Mx::Array<int64_T>(d1,d2,!raw).copyFromTR(data);
};

template <> inline
mxArray* wbMatrix<unsigned char>::toMx(char raw) const {
   return toMx_base(raw); 
};
template <> inline
mxArray* wbMatrix<unsigned char>::toMxT(char raw) const {
   size_t d1=(raw ? dim2:dim1), d2=(raw ? dim1:dim2);
   return Mx::Array<uint8_T>(d1,d2,!raw).copyFromTR(data);
};

template <> inline
mxArray* wbMatrix<char>::toMx(char raw) const {
   return toMx_base(raw); 
};
template <> inline
mxArray* wbMatrix<char>::toMxT(char raw) const {
   size_t d1=(raw ? dim2:dim1), d2=(raw ? dim1:dim2);
   return Mx::Array<int8_T>(d1,d2,!raw).copyFromTR(data);
};

template <> inline
mxArray* wbMatrix<double>::toMx(char raw) const {
   return toMx_base(raw); 
};

template <> inline
mxArray* wbMatrix<wbcomplex>::toMx(char raw) const {
   return toMx_base(raw); 
};

template<class T>
mxArray* wbMatrix<T>::toMx_Struct() const {

   mxArray *S=data->mxCreateStruct(dim1,dim2); 

   size_t i,j;
   for (i=0; i<dim1; ++i)
   for (j=0; j<dim2; ++j) data[j+i*dim2].add2MxStruct(S,i+j*dim1);

   return S;
};

template<class T>
mxArray* wbMatrix<T>::toMxP_S() const {

   mxArray *S;
   if (data) { size_t i,j;
      S=(*data)->mxCreateStruct(dim1,dim2);

      for (i=0; i<dim1; ++i)
      for (j=0; j<dim2; ++j) data[j+i*dim2]->add2MxStruct(S,i+j*dim1);
   }
   else { S=mxCreateStructMatrix(dim1,dim2,0,NULL); }

   return S;
};

template<class T>
mxArray* wbMatrix<T>::toMxP_C() const {

   mxArray *S;

   if (data) { size_t i,j;
      S=(*data)->mxCreateCell(dim1,dim2);

      for (i=0; i<dim1; ++i)
      for (j=0; j<dim2; ++j) data[j+i*dim2]->add2MxCell(S,i+j*dim1);
   }
   else { S=mxCreateCellMatrix(dim1,dim2); }

   return S;
};

template<class T>
mxArray* wbMatrix<T>::mxCreateStruct(unsigned m, unsigned n) const {
   return mxCreateCellMatrix(m,n);
}

template<class T>
void wbMatrix<T>::add2MxStruct(mxArray *S, unsigned i, char tst) const {

   if (tst) {
      size_t s=0; 
      if (!S || (s=mxGetNumberOfElements(S))<1 || i>=s) wblog(FL,
      "ERR %s() must follow mxCreateCell()\n%lx, %d/%d",
       FCT,S,i+1,s);
   }

   mxSetCell(S,i,toMx());
}

template <class T>
void wbMatrix<T>::info(const char *istr) const {

    size_t l=0, n=16; char s[n];
    snprintf(s,n,"%lix%li",dim1,dim2);
    if (l>=n) wblog(FL,"ERR %s() string out of bounds (%d/%d)",FCT,l,n);
    fprintf(stdout,"  %-12s %-10s @ 0x%p  double array\n", istr, s, data);
};

template <class T>
void wbMatrix<T>::print(const char *istr, char mflag) const {

    mxArray *a=toMx_base(0);

    if (istr && istr[0]) {
       if (mflag) { wb_printf("\n%s = [\n",istr); }
       else {
          wb_printf("\n%s = [%dx%d double]\n\n",istr,dim1,dim2);
       }
    } else { wb_printf("\n"); }

    Wb::CallMatlab(0,NULL,1,&a,"disp");
    wb_printf("%s\n", mflag && istr && istr[0] ? "];":"");

    mxDestroyArray(a);
};

template <class T>
void wbMatrix<T>::Print(const char *istr, char mflag) const {

    size_t i,j;

    if (!mflag) {
        if (istr[0])
        printf("\n%s = [%ldx%ld double]\n\n", istr, dim1, dim2);
        else printf("\n");
    }
    else
    printf("\n%s = [\n", istr[0] ? istr : "ans");

    for (i=0; i<dim1; ++i) {
        for (j=0; j<dim2; ++j) printf(" %8.3g", (double)data[i*dim2+j]);
        printf("\n");
    }

    if (mflag) printf("];\n");
}

template <class T>
void wbMatrix<T>::printdata(
    const char *istr, const char *dbl_fmt, const char *rsep
) const {

    size_t i=0,j;
    if (istr && istr[0]) printf("%s = [\n ",istr); else printf("[");

    for (; i<dim1; ++i) { if (i) printf("%s ",rsep);
    for (j=0; j<dim2; ++j) { printf(dbl_fmt, double(data[i*dim2+j])); }}

    printf(" ]\n");
}

template <class T>
void wbMatrix<T>::appendRow(const wbvector<T> &v) {

    if (dim1) {
       if (dim2!=v.len) wblog(FL,
       "ERR Dimension mismatch (%dx%d : %d)", dim1,dim2, v.len);

       appendRows(1, v.data);
    }
    else RENEW(1,v.len,v.data);
};

template <class T>
void wbMatrix<T>::appendRows(size_t n, const T* v) {

    size_t s=dim1*dim2, ds=n*dim2;
    T *d0=data;

    dim1+=n; if (!n || !dim2) return;

    WB_NEW(data,s+ds);
    MEM_CPY<T>(data,s+ds,s,d0,v);
    WB_DELETE(d0);
};

template<class T>
void wbMatrix<T>::getReal(wbMatrix<double> &R) const {
   if (typeid(T)!=typeid(double)) {
      R.init(dim1,dim2);
      for (size_t n=dim1*dim2, i=0; i<n; ++i) { R.data[i]=double(data[i]); }
   }
   else { R.init(*this); }
   return;
};

template<class T>
void wbMatrix<T>::getImag(wbMatrix<double> &I) const {
   wblog(FL,"ERR wbMatrix::getImag not defined for type %s.",
   TSTR(T));
};

template<class T>
void wbMatrix<T>::set(const wbMatrix<double> &R, const wbMatrix<double> &I) {
   wblog(FL,"ERR wbMatrix::set(R,I) not defined for type %s.",
   TSTR(T));
};

template<class T> 
wbMatrix<T>& wbMatrix<T>::Conj() { return *this; };

template<class T>
double wbMatrix<T>::normReal() const { return Wb::sqrt(norm2()); };

template<class T>
double  wbMatrix<T>::normImag() const { return 0; };

template<>
double wbMatrix<wbcomplex>::normReal() const {
    double x=0;
    for (size_t n=dim1*dim2, i=0; i<n; ++i) x+=Wb::norm2(data[i].r);
    return Wb::sqrt(x);
};

template<>
double wbMatrix<wbcomplex>::normImag() const {
    double x=0;
    for (size_t n=dim1*dim2, i=0; i<n; ++i) x+=Wb::norm2(data[i].i);
    return Wb::sqrt(x);
};

template<>
wbMatrix<wbcomplex>& wbMatrix<wbcomplex>::Conj() {
   for (size_t s=dim1*dim2, i=0; i<s; ++i)
   data[i].i = -data[i].i;
   return *this;
}

template<>
void wbMatrix<wbcomplex>::getReal(wbMatrix<double> &R) const {
   R.init(dim1, dim2);

   if (isdiag) {
      size_t i,k,s=MIN(dim1,dim2);
      for (i=0; i<s; ++i) { k=i*dim2+i; R.data[k]=data[k].r; }
      R.isdiag=isdiag;
   }
   else {
      for (size_t s=dim1*dim2, i=0; i<s; ++i)
      R.data[i]=data[i].r;
   }
}

template<>
void wbMatrix<wbcomplex>::getImag(wbMatrix<double> &I) const {
   I.init(dim1, dim2);

   if (isdiag) {
      size_t i,k,s=MIN(dim1,dim2);
      for (i=0; i<s; ++i) { k=i*dim2+i; I.data[k]=data[k].i; }
      I.isdiag=isdiag;
   }
   else {
      for (size_t s=dim1*dim2, i=0; i<s; ++i)
      I.data[i]=data[i].i;
   }
}

template<>
void wbMatrix<wbcomplex>::set(
   const wbMatrix<double> &R,
   const wbMatrix<double> &I
){
   if (!R.sameSize(I)) wblog(FL,
   "ERR Dimension mismatch (%dx%d; %dx%d)!",R.dim1,R.dim2,I.dim1,I.dim2);

   init(R.dim1,R.dim2);

   for (size_t s=dim1*dim2, i=0; i<s; ++i)
   data[i].set(R.data[i], I.data[i]);
}

template <class T> inline
void matchIndexU(const char *F, int L,
    const wbMatrix<T> &QA, const wbMatrix<T> &QB,
    wbindex &Ib, const char force=1
);

template<class T>
int matchIndex(C_TMAT &QA, C_TMAT &QB, wbindex &Ia, wbindex &Ib,
    char lex=1, 
    widx_t *ma=NULL, widx_t *mb=NULL, 
    T eps=0
);

template <class T> inline
int matchIndex(
    const wbvector<T> &qA, const wbvector<T> &qB,
    wbindex &Ia, wbindex &Ib,
    char lex=1, 
    widx_t *ma=NULL, widx_t *mb=NULL, 
    T eps=0
){
    wbMatrix<T> QA, QB;  int r;
    QA.init2ref(qA,'t');
    QB.init2ref(qB,'t'); if (QA.dim2!=1 || QB.dim2!=1) wblog(FL,"ERR");

    r=matchIndex(QA,QB,Ia,Ib,lex,ma,mb,eps);
    return r;
};

template <class T> inline
int matchSortedIdx(
    const T *da, size_t lda, size_t na,
    const T *db, size_t ldb, size_t nb, wbindex &Ia, wbindex &Ib,
    widx_t m=-1,
    char lex=1
){
    widx_t m1=0,m2=0;
    return Wb::matchSortedIdx(da,lda,na, db,ldb,nb, Ia,Ib,m,lex,&m1,&m2);
};

template <class T> inline
int matchSortedIdx(C_TMAT &QA, C_TMAT &QB, wbindex &Ia, wbindex &Ib,
    widx_t m=-1,  
    char lex=1,    
    widx_t *ma=NULL, widx_t *mb=NULL,  
    T eps=0        
){
    return Wb::matchSortedIdx(
       QA.data, QA.dim2, QA.dim1,
       QB.data, QB.dim2, QB.dim1, Ia,Ib, m, lex, ma, mb, eps
    );
};

template <class T> inline
int matchSortedIdx(C_TMAT &QA, C_TMAT &QB,
    wbindex &Ia, wbindex &Ib, wbvector<widx_t> &D,
    widx_t m=-1,  
    char lex=1,    
    T eps=0        
){
    return Wb::matchSortedIdx(
       QA.data, QA.dim2, QA.dim1,
       QB.data, QB.dim2, QB.dim1, Ia,Ib,D, m, lex, eps
    );
};

template <class T> inline
int matchSortedIdxU( 
    const char *F, int L,
    C_TMAT &QA, C_TMAT &QB, wbindex &Ia, wbindex &Ib,
    widx_t m=-1, char lex=1
);

template <class T> inline
int matchSortedIdxU( 
    C_TMAT &QA, C_TMAT &QB, wbindex &Ia, wbindex &Ib, widx_t m=-1
){  return matchSortedIdxU(FL,QA,QB,Ia,Ib,m); };

#endif

