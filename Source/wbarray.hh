/* ---------------------------------------------------------------------
 * Project : QSpace tensor library (v4.0)
 * Class   : wbarray (array class, col-major)
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

#ifndef __WB_ARRAY_COL_MAJOR_HH__
#define __WB_ARRAY_COL_MAJOR_HH__

// #define CHECK_ELEMENT_RANGE

namespace Wb {

   template<class T>
   wbstring sizeStrM( 
      unsigned r, const T* sd, unsigned m, const T* sm,
      const char *sep1="x", const char *sepM="_"); 

   template<class T> inline
   void householder(const T *u, T *x, size_t n, const T& eps);

   template<class T> inline
   void Givens_eraseCol(
      T* a, unsigned i0, unsigned j0, unsigned dim1, unsigned dim2,
      T *x=nullptr, T eps=1e-14);
};

template <class T>
class wbarray { 

  public:

    wbarray() : mtype(Wb::MEM_DEF), sptr(nullptr), data(nullptr) { };

    wbarray(size_t d1) : wbarray() { init(d1); };
    wbarray(size_t d1, size_t d2) : wbarray() { init(d1,d2); };
    wbarray(size_t d1, size_t d2, size_t d3) : wbarray() { init(d1,d2,d3); };
    wbarray(size_t d1, size_t d2, size_t d3, size_t d4)
     : wbarray() { init(d1,d2,d3,d4); };

    wbarray(char *sidx) : wbarray() {
       WBINDEX S; Wb::Str2Idx(FL,sidx,S);
       init(S);
    };
    wbarray(const wbvector<unsigned> &S) : wbarray() { init(S); };

    wbarray(const wbvector<size_t> &S, const T *d0, char ref=0) : wbarray() {
       if (ref)  
            { INIT2REF(d0,&S); }
       else { NEW(S,d0); }
    };

    wbarray(const wbarray &a, char ref=0)
     : wbarray() { init(a,ref); };

    template<class T2>
    wbarray(const wbarray<T2> &a) : wbarray() {
       NEW(a.SIZE); init_T(a.data); 
    };

    wbarray(const mxArray* a, char ref=0, char vec=0)
     : wbarray() { init(a,ref,vec); };

    wbarray(const char *F, int L, const mxArray* a,
       char ref=0, 
       char vec=0)
     : wbarray() { init(F,L,a,ref,vec); };

    wbarray(const wbMatrix<T> &a, char ref=0)
     : wbarray() { init(a,ref); };

    virtual ~wbarray() {
       DELETE_DATA(); 
    };

    wbarray& initX() {
       DELETE_DATA('!'); SIZE.init();
       return *this;
    };

    wbarray& init() {
       DELETE_DATA();
       SIZE.init(); return *this;
    };

    wbarray& init(size_t d1, const T* d=nullptr) {
       if (d1) { NEW(wbvector<size_t>(1,&d1,'r'),d); } else { init(); }
       return *this;
    };

    wbarray& init(size_t d1, size_t d2, const T* d=nullptr) {
       size_t s[2]= {d1,d2}; NEW(wbvector<size_t>(2,s,'r'),d);
       return *this;
    };

    wbarray& init(size_t d1, size_t d2, size_t d3, const T* d=nullptr) {
       size_t s[3]= {d1,d2,d3}; NEW(wbvector<size_t>(3,s,'r'),d);
       return *this;
    };

    wbarray& init(size_t d1, size_t d2, size_t d3, size_t d4, const T* d=nullptr) {
       size_t s[4]= {d1,d2,d3,d4}; NEW(wbvector<size_t>(4,s,'r'),d);
       return *this;
    };

    wbarray& init(const wbvector<size_t> &S, const T *d0, char ref=0) {
       if (ref)
            { INIT2REF(d0,&S); }
       else { NEW(S,d0); }
       return *this;
    };

    wbarray& init2ref(const wbarray &a) {
       return INIT2REF(a.data,&a.SIZE);
    };

    wbarray& initScalar(T x=0, unsigned r=2) {
       if (r==2) { init(1,1,&x); } else
       if (r==1) { init(1); data[0]=x; } else
       if (!r) wblog(FL,"ERR %s() got r=%d (x=%g)",FCT,r,double(x));
       else { wbvector<size_t> S; NEW(S.init2val(r,1)); data[0]=x; }
       return *this;
    };

    wbarray& init(const wbarray &a, char ref=0) {
       if (this!=&a) {
          if (!ref) { NEW(a.SIZE, a.data); }
          else { DELETE_DATA(); SIZE=a.SIZE;
             a.check_consistency(FL);
             if (a.sptr) {
                (sptr=a.sptr)->add_dref();
                data=a.data;
             }
          }
       }
       return *this;
    };

    template<class T_>
    wbarray& set(const wbarray &A, T_ fac_) {
       if (fac_==T_(1) ) { init(A); } 
       else if (!A || !fac_) { NEW(A.SIZE); }
       else {
          size_t i=0, n=A.SIZE.prod(0);
          const T *a=A.data; double fac=fac_;
          SIZE=A.SIZE; NEW_DATA(n,nullptr,'\0','\0'); 

          if (fac==-1)
               { for (; i<n; ++i) data[i]=   -a[i]; }
          else { for (; i<n; ++i) data[i]=fac*a[i]; }
       }
       return *this;
    };

    wbarray& init(const wbMatrix<T> &a, char ref=0) { 
       if (ref) {
          size_t s[]={ a.dim2, a.dim1 }; SIZE.init(2,s);
          INIT2REF(a.data);
       }
       else {
          size_t i,j, l=0; init(a.dim1,a.dim2);
          for (i=0; i<a.dim1; ++i)
          for (j=0; j<a.dim2; ++j, ++l) { data[i+a.dim1*j]=a.data[l]; }
       }
       return *this;
    };

    wbindex& ind2sub(size_t k, wbindex &I) const; 

    void adjustMMat(
       const char* F, int L,
       char flag, const T fac, wbarray &A
    ) const {

       if (!isMatrix()) wblog(F,L,
          "ERR %d - rank-2 tensor required (%d)",FCT,SIZE.len);
       if (!strchr("NTC",flag)) wblog(F,L,
          "ERR %d - invalid flag %c<%d>",FCT,flag,flag);

       size_t s=numel();

       if (flag=='C' && typeid(T)!=typeid(wbcomplex)) flag='T';

       if (fac!=T(1) || flag!='N') { size_t i;
          if (flag!='N') permute(A,"21"); else A=(*this);
          if (flag=='C') for (i=0; i<s; i++) A[i]=Wb::CONJ(A[i]);
          if (fac!=T(1)) for (i=0; i<s; i++) A[i]*=fac;
       }
       else { A.init2ref(*this); } 
    }

    void adjustDMat(
       const char* F, int L,
       char flag, const T fac, wbarray &A
    ) const {

       if (SIZE.len!=1) wblog(F,L,
          "ERR %d - 1D representation required (%d)",FCT,SIZE.len);
       if (!strchr("NTC",flag)) wblog(F,L,
          "ERR %d - invalid flag %c<%d>",FCT,flag,flag);

       size_t n=SIZE[0];

       if (flag=='C' && typeid(T)!=typeid(wbcomplex)) flag='T';

       if (fac!=T(1) || flag=='C') { size_t i; A=(*this);
          if (flag=='C') for (i=0; i<n; i++) A[i]=Wb::CONJ(A[i]);
          if (fac!=T(1))    for (i=0; i<n; i++) A[i]*=fac;
       }
       else { A.init2ref(*this); } 
    }

    void adjustVec(
       const char* F, int L, char &flag, wbarray &A, char pos
     ) const { A=(*this); A.adjustVec(F,L,flag,pos); };

    void adjustVec(const char* F, int L, char &flag, char pos) {
       if (SIZE.len!=1) wblog(F,L,
          "ERR %s() vector expected (%d)",FCT,SIZE.len);
       if (!strchr("NTC",flag)) wblog(F,L,
          "ERR %s() invalid flag %c<%d>",FCT,flag,flag);

       if (flag=='C' && typeid(T)!=typeid(wbcomplex)) flag='T';
       if (flag=='C') {
          for (size_t s=numel(), i=0; i<s; ++i)
          data[i]=Wb::CONJ(data[i]);
       }
       flag='N'; 

       SIZE.Resize(2);
       if (pos==1) { SIZE[1]=SIZE[0]; SIZE[0]=1; } else
       if (pos==2) { SIZE[1]=1; } else
       wblog(F,L,"ERR %s() invalid pos=%d",FCT,pos);
    };

    void adjustCMat(const char* F, int L,
       T cfac, size_t s1, size_t s2, char cforce
    ){
       if (cfac!=T(0)) {
           if (data==nullptr) {
              if (cforce) wblog(F,L,      
             "WRN C = A*B + c*[] with c=%s !?", Wb::num2Str(cfac).data);
           }
           else if (!isMatrix() || SIZE[0]!=s1 || SIZE[1]!=s2) {
              wblog(F,L,"ERR %s() dimension mismatch: C=(%s) =? (%d,%d).",
              FCT, SSTR(*this), s1, s2); return;
           }
           else {
               if (cfac!=T(1)) { 
                  size_t i, s=SIZE.prod(0);
                  for (i=0; i<s; i++) data[i]*=cfac;
               }
               return;
           }
       }
       init(); 
    };

    void setBlock(
        const wbvector< WBINDEX > &D, 
        const wbindex &IB,   
        const wbarray &A
    );

    void addBlock(
        const WBINDEX &I, 
        const wbarray &A, const char dflag=0
    );

    wbarray& addBlock(  
       size_t i, size_t j, 
       size_t n, size_t m, 
       wbarray &a
    ) const;

    T norm2block(  
       size_t i, size_t j, 
       size_t n, size_t m  
    ) const;

    wbarray& getBlock(  
       size_t i, size_t j, 
       size_t n, size_t m, 
       wbarray &a
    ) const;

    wbarray& BlockDiag(const wbvector< wbarray > &D);

    wbarray& blockTrace( 
       size_t D, 
       wbarray &a, size_t D2=-1 
    ) const;

    void copyStride(T* dd, size_t stride, T afac=1) const;
    void addStride (T* dd, size_t stride, T afac=1) const;

    template<class T0>
    wbarray& init_T(const T0 *dd) {
       if (SIZE.len) {
       for (size_t n=SIZE.prod(), i=0; i<n; ++i) { data[i]=T(dd[i]); }}
       return *this;
    };

    template<class T0>
    wbarray& init_T(const char *F, int L, const T0 *dd) {
       for (size_t n=SIZE.prod(), i=0; i<n; ++i) { data[i]=T(dd[i]);
          if (T0(data[i])!=dd[i]) wblog(F_L,
             "ERR type conversion changes value (%g->%g) !?",
              double(dd[i]), double(data[i])
          );
       }
       return *this;
    };

    wbarray& initT(const wbarray<T> &a) { return init(a); };

    template<class T0>
    wbarray& initT(const wbarray<T0> &a) {
       NEW(a.SIZE); init_T(a.data);
       return *this;
    };

    template <class T0>
    wbarray& initT(const char *F, int L, const wbarray<T0> &a) {
       NEW(a.SIZE); init_T(F_L,a.data); return *this;
    };

    wbarray& init(char *s, T* d=nullptr) {
       WBINDEX S; Wb::Str2Idx(FL,s,S);
       init(S,d); return *this;
    };

    wbarray& init(const WBINDEX &S, T* d=nullptr) {
       return NEW(S,d); 
    };

    template <class IT>
    wbarray& init(const wbvector<IT> &S_, T* d=nullptr) {
       WBINDEX S(S_.len);
       for (size_t i=0; i<S.len; ++i) { S.data[i]=size_t(S_.data[i]); }
       return NEW(S,d); 
    };

    wbarray& init_bare(const WBINDEX &S) {
       size_t len=S.prod(0);
       if (len)
            { SIZE=S; NEW_DATA(len,nullptr,'\0','\0'); } 
       else { init(); }
       return *this;
    };

    wbarray& init_bare(size_t d1, size_t d2) { 
       wbvector<size_t> S(2); S[0]=d1; S[1]=d2;
       return init_bare(S); 
    };

    wbarray& init(const char *F, int L,
    const mxArray *a, char ref=0, char vec=0);  

    wbarray& init(const char *F, int L,
    const mxArray *a, wbperm &P0); 

    wbarray& init(const mxArray *a, char ref=0, char vec=0) {
       return init(FL,a,ref,vec); }
    wbarray& init(const mxArray *a, wbperm &P0) {
       return init(FL,a,P0); }

    wbarray& init2ref(const T *x, const WBINDEX &S) {
       SIZE=S; return INIT2REF(x);
    };
    wbarray& init2ref(const T *x) {
       return INIT2REF(x);
    };

    wbarray& init2ref(size_t d1, const T *x) {
       SIZE.init(1,&d1);
       return INIT2REF(x);
    };
    wbarray& init2ref(size_t d1, size_t d2, const T *x) {
       size_t s[]= {d1,d2}; SIZE.init(2,s);
       return INIT2REF(x);
    };
    wbarray& init2ref(size_t d1, size_t d2, size_t d3, const T *x) {
       size_t s[]= {d1,d2,d3}; SIZE.init(3,s);
       return INIT2REF(x);
    };

    int isRef() const { return sptr ? sptr->isRef() : 0; };

    wbarray& Instantiate(const char *F=nullptr, int L=0) {
	   if (isRef() && data) {
          if (!sptr) wblog(F_L,"ERR %s() got null sptr !?",FCT);
          T *d0=data; NEW_DATA(SIZE.prod(0),d0);
       }
       else if (F) wblog(F,L,
          "WRN %s() no need to instantiate %s",SSTR(*this));

       return *this;
    };

    wbarray& unRef() { return Instantiate(); };

    void init2Vec(const T *d0, size_t len, char ref=0) {
       SIZE.init(1,&len);
       if (data!=d0) { NEW_DATA(len,d0,ref); }
    };

    void init2ref(const wbarray &A, size_t s1, size_t s2) {
        if (s1*s2!=A.SIZE.prod(0)) wblog(FL,
          "ERR %s() cannot reshape array %s into %dx%d",
           SSTR(*this), s1, s2
        );

        SIZE.init(2); SIZE[0]=s1; SIZE[1]=s2;
        INIT2REF(A.data);
    };

    void initTst();

    void init2Vec(const wbvector<T> &a) {
       WBINDEX S(1); S[0]=a.len;
       NEW(S, a.data);
    };

    wbarray& initIdentity(size_t d, char dflag=0, T dval=1);
    wbarray& initIdentity(const WBINDEX &S, char dflag=0);
    wbarray& Reduce2Id(); 

    wbarray& initIdentityB(
        size_t d1, size_t d2, size_t k=0);

    wbarray& initIdentityB3(
        size_t d1, size_t d2, size_t D, size_t i0, T x=1);

    wbarray& initPermB(size_t d1, size_t d2, const wbperm &P);

    wbarray& Expand2Projector(const T eps=0);
    wbarray& ExpandDiagonal(unsigned i1, unsigned i2); 
    wbarray& ExpandDiagonal();

    void Diag2Vec();
    wbarray& Reduce2Diag(T eps=0);

    wbarray& initDiag(unsigned n, const T *d){
       init(n,n);
          for (size_t i=0; i<n; ++i) data[i+i*n]=d[i];
       return *this;
    };

    wbarray& initDiag(const wbvector<T> &d){
       return initDiag(d.len, d.data); };

    void initVec(size_t d, const T val=0) {
       init(d); if (val) set(val);
    };

    wbarray& resize(
       const wbvector<size_t> &S, wbarray &B) const;
    wbarray& resize( 
       const wbvector<size_t> &S, wbarray &B, double &dx) const;

    wbarray& resizeNumCols(size_t n, wbarray &B) const { 
       if (!n && isEmpty()) { return B.init(); }
       if (SIZE.len!=2) wblog(FL,
          "ERR %s() for matrices only (got %s)",FCT,SSTR(*this));
       wbvector<size_t> S(SIZE); S[1]=n;
       return resize(S,B);
    };

    wbarray& Resize(const wbvector<size_t> &S, const T* d0=nullptr) { 
       if (d0) { NEW(S,d0); } else
       if (S!=SIZE) {
          wbarray<T> X;
          resize(S,X).save2(*this);
       };
       return *this;
    };

    wbarray& Resize(size_t s1, size_t s2) {
       size_t s[2]={s1,s2}; wbarray X;
       if (SIZE.len!=2) wblog(FL,"ERR got rank %d (matrix required)",SIZE.len);
       resize(wbvector<size_t>(2,s),X).save2(*this);
       return *this;
    };
    wbarray& Resize(size_t s1, size_t s2, double &dx) {
       size_t s[2]={s1,s2}; wbarray X;
       if (SIZE.len!=2) wblog(FL,"ERR got rank %d (matrix required)",SIZE.len);
       resize(wbvector<size_t>(2,s),X,dx).save2(*this);
       return *this;
    };

    wbarray& Resize1(unsigned i, size_t Si) {
       if (i>=SIZE.len) wblog(FL,
          "ERR %s() index out of bounds (%d/%d)",FCT,i,SIZE.len);
       if (SIZE[i]!=Si) {
          wbvector<size_t> S(SIZE); wbarray X; S[i]=Si;
          resize(S,X).save2(*this);
       }
       return *this;
    };

    wbarray& ResizeM(unsigned r, const wbvector<size_t> &Sm);

    wbarray& Enlarge(
       size_t dim,     
       size_t D,       
       size_t istart=0 
    );

    void swap(wbarray &A) { 
       SWAP(sptr, A.sptr );
       SWAP(data, A.data );
       SWAP(mtype,A.mtype); SIZE.swap(A.SIZE);
    };

    wbarray& save2(wbarray &A) {
       if (this!=&A) {
          if (A.mtype) wblog(FL,"ERR %s() got %s / %s",
             FCT, Wb::MTYPE_STR[mtype], Wb::MTYPE_STR[A.mtype]);
          this->swap(A); if (data!=A.data) { init(); }
       }
       return A;
    };

    wbarray& save2(wbMatrix<T>&);

    wbarray& Append2(
       const char *F, int L, wbarray &B, unsigned dim); 
    wbarray& Append2(wbarray &B, unsigned dim){
       return Append2(0,0,B,dim);
    }

    void Append2(
       const char *F, int L,
       wbarray &B, unsigned dim, 
       wbarray<double> &W, double eps 
    );
    void Append2(wbarray &B, unsigned dim, wbarray &W, double eps){
       return Append2(0,0,B,dim,W,eps);
    };

    void rand(const UVEC &d) { init(d); setRand(); };
    void setRand(const char pnflag=0);

    wbarray& Reset(const T *dref) { 
       size_t n=SIZE.prod();
       if (dref)
            { Wb::MemCpy(data,dref,n); } 
       else { MEM_SET<T>(data,n); }
       return *this;
    };

    wbarray& set(const T& c) { 
       size_t n=SIZE.prod();
       if (c)
            { for (size_t i=0; i<n; ++i) data[i]=c; }
       else { MEM_SET<T>(data,n); }
       return *this;
    };

    wbarray& operator= (const wbarray &a) { return init(a); };

    template <class T0> 
    wbarray<T>& operator= (const wbarray<T0> &a) {
       if ((void*)this!=(void*)&a) {
          NEW(a.SIZE);
          for (size_t n=a.numel(), i=0; i<n; ++i) data[i]=T(a.data[i]);
       }
       return *this;
    };

    void SetRef(const widx_t *I, size_t len, const T* d0);
    void SetRef(size_t k, const T* d0) { SetRef(&k,1,d0); };

    const T* ref(const widx_t *I, size_t n) const {
       if (!n) return nullptr; 
       return data+serial_index(I,n);
    };
    T* ref(const widx_t *I, size_t n) {
       if (!n) return nullptr; 
       return data+serial_index(I,n);
    };

    T* col(size_t k) {
       if (SIZE.len!=2 || k>=SIZE[1]) wblog(FL,
          "ERR %s() index out of bounds (%s; %d)",FCT,SSTR(*this),k);
       return (data+k*SIZE[0]); 
    };
    const T* col(size_t k) const {
       if (SIZE.len!=2 || k>=SIZE[1]) wblog(FL,
          "ERR %s() index out of bounds (%s; %d)",FCT,SSTR(*this),k);
       return (data+k*SIZE[0]); 
    };

    T* row(size_t k) {
       if (SIZE.len!=2 || k>=SIZE[0]) wblog(FL,
          "ERR %s() index out of bounds (%d / %s)",FCT,k,SSTR(*this));
       return (data+k); 
    };
    const T* row(size_t k) const {
       if (SIZE.len!=2 || k>=SIZE[0]) wblog(FL,
          "ERR %s() index out of bounds (%d / %s)",FCT,k,SSTR(*this));
       return (data+k); 
    };

    wbvector<T>& getCol(size_t k, wbvector<T> &v) const;
    wbarray<T>&  getCol(size_t k, wbarray<T>  &a) const;

    wbvector<T>& getRow(size_t k, wbvector<T> &v) const;
    wbarray<T>&  getRow(size_t k, wbarray<T>  &a) const;

    wbvector<T> getDiag() const; 

    T* ref(const WBINDEX &I) { return ref(I.data,I.len); };
    const T* ref(const WBINDEX &I) const {
        return ref(I.data,I.len);
    };

    const T* ref(widx_t i) const { return ref(&i,1); }
    T* ref(widx_t i) { return ref(&i,1); }

    const T* ref(widx_t i1, widx_t i2) const {;
       widx_t I[]={i1,i2}; return ref(I,2);
    };
    T* ref(widx_t i1, widx_t i2) {;
       widx_t I[]={i1,i2}; return ref(I,2);
    };

    const T* ref(widx_t i1, widx_t i2, widx_t i3) const {
       widx_t I[]={i1,i2,i3}; return ref(I,3);
    };
    T* ref(widx_t i1, widx_t i2, widx_t i3) {
       widx_t I[]={i1,i2,i3}; return ref(I,3);
    };

    const T& operator[] (size_t i) const {
       return data[i];
    };
    T& operator[] (size_t i) {
       return data[i];
    }

    T& last(size_t i) {  
       size_t n=numel(); 
       if (!n) wblog(FL,
          "ERR %s() got empty %s array (%s)",FCT,TSTR(T),SSTR(*this));
       return data[n-1];
    };

    const T& operator() (size_t i) const { 
        if (SIZE.len!=1) wblog(FL,
           "ERR %s requires vector type (%s)",FCT,SSTR(*this));
        return data[i];
    };
    T& operator() (size_t i) {
        if (SIZE.len!=1) wblog(FL,
           "ERR %s requires vector type (%s)",FCT,SSTR(*this));
        return data[i];
    };

    const T& operator() (size_t i, size_t j) const {
        if (SIZE.len!=2) wblog(FL,
           "ERR %s(i,j) requires matrix type (%s)",FCT,SSTR(*this));
        return data[ i + j*SIZE[0] ]; 
    };
    T& operator() (size_t i, size_t j) {
        if (SIZE.len!=2) wblog(FL,
           "ERR %s(i,j) requires matrix type (%s)",FCT,SSTR(*this));
        return data[ i + j*SIZE[0] ]; 
    };

    const T& operator() (size_t i, size_t j, size_t k) const {
        if (SIZE.len!=3) wblog(FL,
           "ERR %s(i,j,k) requires 3D-object (%s)",FCT,SSTR(*this));
        return data[ i + SIZE[0] * (j + SIZE[1]*k) ];
    };
    T& operator() (size_t i, size_t j, size_t k) {
        if (SIZE.len!=3) wblog(FL,
           "ERR %s(i,j,k) requires 3D-object (%s)",FCT,SSTR(*this));
        return data[ i + SIZE[0] * (j + SIZE[1]*k) ];
    };

    const T& operator() (const WBINDEX &I) const {
       if (I.len!=SIZE.len) wblog(FL,
          "ERR %s() invalid index [%s]",FCT,STR(I));
       return data[serial_index(I.data,I.len)];
    };
    T& operator() (const WBINDEX &I) {
       if (I.len!=SIZE.len) wblog(FL,
       "ERR invalid index [%s] having %s",STR(I),SSTR(*this));
       return data[serial_index(I.data,I.len)];
    };

    const T& element(const WBINDEX &I) const {
       if (I.len!=SIZE.len) wblog(FL,"ERR %s() invalid index [%s] "
          "having %s",FCT,STR(I),SSTR(*this));
       for (size_t i=0; i<I.len; i++) if (I.data[i]>=SIZE.data[i])
           wblog(FL,"ERR %s() index out of bounds (%s; %s)",
           FCT,STR(I),SSTR(*this));
       return data[serial_index(I.data,I.len)];
    };
    T& element(const WBINDEX &I) {
       if (I.len!=SIZE.len) wblog(FL,"ERR %s() invalid index [%s] "
          "having %s",FCT,STR(I),SSTR(*this));
       for (size_t i=0; i<I.len; i++) if (I.data[i]>=SIZE.data[i])
           wblog(FL,"ERR %s() index out of bounds (%s; %s)",
           FCT,STR(I),SSTR(*this));
       return data[serial_index(I.data,I.len)];
    };

    size_t numel(unsigned r0=-1) const { 
       size_t N=0; 
       if (SIZE.len) {
          unsigned i=0, n=SIZE.len;
          const size_t *s=SIZE.data; N=1;

          if (int(r0)>=0) { if (r0>n) { r0=n; }
             for (; i<r0; ++i) { if (!s[i]) { N=0; break; }}
          }
          if (N && i<n) {
             for (N=s[i++]; i<n; ++i) { N*=s[i]; }
          }
       }
       return N;
    };

    size_t numOM(unsigned r0) const {
       size_t M=numel(r0); 
       if (int(r0)<0 || !M || (r0<=2 && M>1)) wblog(FL, 
          "ERR %s() got OM=%d for %s @ r0=%d",FCT,M,SSTR(*this),r0);
       return M;
    };

    wbstring sizeStr(unsigned stride=0, const char *sep2="") const {
       if (SIZE.len)
            { return SIZE.toStrf("","x",stride,sep2); }
       else { return wbstring("[]"); }
    };

    wbarray& Squeeze(); 

    wbarray& initSingleton(unsigned r, T x=1, unsigned m=0) {
       wbvector<size_t> S;
       if (!m) { S.init2val(r,1); }
       else {
           S.init2val(r+1,1); S[r]=m;
           if (m>1) { x/=sqrt(double(m)); }
       }
       init(S); set(x);
       return *this;
    };

    int SkipSingleton(const char *F, int L, unsigned i); 

    int skipSingletons(const char *F, int L, unsigned r=-1);
    int skipSingletons(unsigned r=-1) {
       return skipSingletons(FL,r); 
    };

    wbarray& appendSingletons(
       const char *F, int L,  unsigned r, unsigned m=0);
    wbarray& appendSingletons(unsigned r, unsigned m=0) {
       return appendSingletons(0,0,r,m); };

    void prependSingletons(unsigned r);

    wbarray& ExpandOM(const char *F, int L,  unsigned r0, 
       const wbvector<unsigned> &S,
       const unsigned *sx=nullptr);

    wbarray& ExpandOM(const char *F, int L, unsigned r0,
       unsigned m, const unsigned *s, const unsigned *sx=nullptr
    ){ return ExpandOM(F,L,r0, wbvector<unsigned>(m,s,'r'), sx); };

    int ExpandOM(const char *F, int L, wbarray &B, unsigned r0);

    wbarray& FuseOM(const char *F, int L, unsigned r, unsigned mc=-1) {
       if (!SIZE.len) { return *this; }

       unsigned i, m=1;
       for (i=r; i<SIZE.len; ++i) { m*=SIZE[i]; }

       if (r<=2 && m>1) wblog(F_L,
          "ERR %s() requesting OM m=%d for rank r=%d/%d",FCT,m,r,SIZE.len);
       if (r>SIZE.len) wblog(F_L,
          "ERR %s() rank out of bounds (r=%d/%d; m=%d)",FCT,r,SIZE.len,m);
       if (int(mc)>=0 && m!=mc && (m!=1 || mc>1)) wblog(F_L,"ERR %s() "
          "got OM mismatch %s @ %d/%d (r=%d)",FCT,SSTRM(SIZE,r),m,mc,r);

       if (m==1 && r)
            { SIZE.len=r; } 
       else { SIZE[r]=m; SIZE.len=r+1; } 

       return *this;
    };

    size_t length() const {
       size_t l=0; const size_t *s=SIZE.data;
       for (unsigned i=0; i<SIZE.len; ++i) {
          if (l<s[i]) { l=s[i]; } else
          if (!s[i]) { l=0; break; }
       }
       return l;
    }

    size_t size(unsigned i) const { return SIZE[i]; };

    size_t dim(unsigned i) const { 
       if (i==0 || i>SIZE.len) wblog(FL,
          "ERR index out of bounds (%d/%d)", i-1, SIZE.len);
       return SIZE[i-1];
    };
    size_t dim0(unsigned i) const { 
       if (i>=SIZE.len) wblog(FL,
          "ERR index out of bounds (%d/%d)", i, SIZE.len);
       return SIZE[i];
    };

    size_t dim1() const {
       size_t D=0; 
       if (SIZE.len==1) { D=SIZE[0]; } else
       if (SIZE.len) {
          unsigned i=0, n=SIZE.len, z=0, l=-1; const size_t *s=SIZE.data;
          for (; i<n; ++i) {
             if (s[i]>1) {
                if (l>i) { l=i; } else
                wblog(FL,"ERR %s() invalid vector %s",FCT,SSTR(*this));
             }
             else if (!s[i]) { ++z; } 
          }

          if (!z) { D=s[l<n ? l : 0]; } else
          if (n>2 || l<n) { 
             wblog(FL,"ERR %s() invalid vector %s",FCT,SSTR(*this));
          }
       }
       return D;
    };

    size_t dim2() const { 
       if (!SIZE.len) return 0;
       if (SIZE.len!=2) wblog(FL,
          "ERR %s() got rank-%d array (%s) !?",FCT,SIZE.len,SSTR(*this));
       if (SIZE[0]!=SIZE[1]) wblog(FL,
          "ERR %s() got non-square matrix (%s) !?",FCT,SSTR(*this));
       return SIZE[0];
    };

    void getMatSize(const char *F, int L, size_t &dim1, size_t &dim2) const;
    void getMatSize(size_t &dim1, size_t &dim2, unsigned m) const;
    void getMatSize(const ctrIdx &ic, size_t &dimc, size_t &dimk,
        unsigned r=-1, size_t *dom=nullptr) const;

    unsigned rank() const { return SIZE.len; }; 

    bool isRank(unsigned r) const { 
       bool q=1;
       if (r!=SIZE.len) { unsigned i=r;
          for (; i<SIZE.len; ++i) { if (SIZE[i]!=1) break; }
          if (i!=SIZE.len) { q=0; } 
       }
       return q;
    };

    bool isRankM(unsigned r, unsigned m) const {  
       bool q=1;
       if (SIZE.len==r) { return q; }
       if (r>2 && m) { if (SIZE.len==r+1 || SIZE.len==r+m ) { return q; }}
       return (q=0);
    };

    bool isOpS(size_t *n=nullptr) const; 

    wbarray& SkipTiny_float(const T eps QS_UNUSED_VAR =1e-14) { return *this; };
    double   SkipTiny_imag (double  eps QS_UNUSED_VAR =1e-14) { return 0; };

    double SkipTiny(double eps_=1e-14);

    wbarray<T>& SkipTrailingZeroSpace(
       unsigned k,  
       size_t n=0,  
       double eps=1e-14);

    bool isMatrix() const { return SIZE.len==2; }
    bool isMatrix(size_t d1, size_t d2) { 
       return (SIZE.len==2 && SIZE[0]==d1 && SIZE[1]==d2);
    };

    bool isSMatrix(unsigned d=-1) const {
       bool q=0;
       if (SIZE.len==2 && SIZE[0]==SIZE[1] &&
          (int(d)<0 || d==SIZE[0])) { q=1; }
       return q;
    };

    bool isSquare() const { 
       if (SIZE.len%2) return 0;
       else {
          for (unsigned i=0, n=SIZE.len/2; i<n; ++i) {
          if (SIZE[i]!=SIZE[n+i]) return 0; }
       }
       return 1;
    };

    bool isEmpty() const { return (data==nullptr); };
    bool isComplex() const;
    bool isZero(double eps=0., char flag=0) const;

    bool isFinite() const { 
       bool q=1; size_t i=0, n=numel();
       for (; i<n; ++i) { if (!isfinite(data[i])) { q=0; break; }
       }; return q;
    };

    explicit operator bool() const { return !isEmpty(); }
    bool operator! () const { return isEmpty(); }

    wbarray& Symmetrize(
       const char *F, int L,
       double *delta=nullptr, char cflag=1, char tflag=0, char fflag=0
    );

    wbarray& swapRows(size_t i1, size_t i2);
    wbarray& swapCols(size_t j1, size_t j2);

    wbarray& setCol(size_t k, size_t k0 );
    wbarray& setCol(size_t k, const T* d, T fac=1, size_t stride=1);
    wbarray& setRow(size_t k, const T* d);

    wbarray& setCol(size_t k, const wbvector<T> &v);
    template<class T2>
    wbarray<T>& setCol(size_t k, const wbvector<T2> &v);

    wbvector<T>& colNorm2(wbvector<T> &a) const;

    T colNorm2(size_t k) const;

    wbvector<T>& normDim2(unsigned k, wbvector<T> &xk) const;  
    wbvector<T>& normDim (unsigned k, wbvector<T> &xk) const {
       normDim2(k,xk);
       for (size_t n=xk.len, i=0; i<n; ++i) { xk[i]=Wb::sqrt(xk[i]); }
       return xk;
    };

    unsigned QRdecomp(
       const char *F, int L, wbarray<T> &Q, wbarray<T> &R,
       char useP=0, T eps=1e-15);

    wbarray& mldivide( const char *F, int L, 
       const wbarray<T> &B, wbarray<T> &X, T bfac=1, T eps=1e-15) const;

    int Householder(
       const char *F, int L,
       size_t k, T *u, 
       wbperm *P=nullptr, 
       T eps=1e-15);

    unsigned eigTriDiag(
       const char *F, int L, wbvector<T> &E,
       wbarray<T> *U=nullptr, T eps=1e-15) const;

    wbarray& ColProject(
       size_t k1, size_t k2, char nflag=0, char tnorm=0);

    wbarray& ColPermute(const wbperm &P); 

    bool isOrthogonalCol(size_t k, char tnorm=0) const;
    bool isOrthoCols(T *x2=nullptr, char tnorm=0, T eps=1e-14) const;

    wbarray& SignCol(size_t k, const T *d0, char tnorm=0);
    wbarray& FlipSignCol(size_t k);
    wbarray& SignConventionCol(
       size_t k=-1,    
       double eps=1e-12  
    );

    T NormalizeCol(size_t k, char tnorm=0, char qflag=0);
    wbarray& NormalizeCols(
       const char *F, int L,
       double *amin=nullptr, double *amax=nullptr, char tnorm=0
    );

    wbarray& OrthoNormalizeCols(
       const char *F=nullptr, int L=0, char tnorm=0,
       char qxflag=0, double eps=1e-14, unsigned np=1
    );

    wbarray& balanceOp(
       const char *F, int L, T &xref, double &xscale
    );

    bool isProptoId(T &x, T eps=1e-14) const;

    bool isProptoId(T *x_=nullptr, T eps=1e-14) const { 
       bool q=0; T x=T(0);
       if (isProptoId(x,eps)) { q=1; if (x_) { (*x_)=x; }}
       return q;
    };

    bool isDiagMatrix(T eps=1e-14) const {
       return isDiag_aux(eps,"isDiag");
    };

    bool isIdentityMatrix(T eps=1e-14) const { T one(1);
       if (data && Wb::abs(data[0]-one)<eps) { 
          if ((SIZE.len%2)==0) {
             if (SIZE.allEqual(1)) { return 1; } 
             else { 
                return isProptoId(one,eps);
             }
          }
       }
       return 0; 
    };

    bool isDiagMatrix(double *eps) const {
       return isDiag_aux(eps,"isDiag"); }

    bool isIdentityMatrix(double *eps) const {
       return isDiag_aux(eps,"isIdty"); }

    size_t nnz(const T eps=0, WBINDEX *I=nullptr) const;

    size_t ndiag() const { size_t n=0;
       if (SIZE.len==1) { if (SIZE[0]) { n=1; }} else
       if (SIZE.len>1) {
          size_t i=0,l; SIZE.max(&l);
          for (n=1; i<SIZE.len; ++i) { if (i!=l) { n*=SIZE[i]; }}
       }
       return n;
    };

    double sparsity(const T eps=0) const { 
       size_t n=numel();
       return (n ? nnz(eps)/double(n) : 0);
    };

    bool isHConj(
      const wbarray &B, double eps=1e-12, double* xref=nullptr
    ) const { return isSym_aux(FLF, B, eps, xref, 's'); };

    bool isHConj(
      double eps=1e-12, double* xref=nullptr
    ) const { return isSym_aux(FLF,*this,eps,xref,'s'); };

    bool isAHerm(
      const wbarray &B, double eps=1e-12, double* xref=nullptr
    ) const { return isSym_aux(FLF, B, eps, xref, 'a'); };

    bool isAHerm(
      double eps=1e-12, double* xref=nullptr
    ) const { return isSym_aux(FLF,*this,eps,xref,'a'); };

    char hasGroupSize(size_t s1, size_t s2) const;

    template<class TI>
    int aMax(T &x, TI &k) const; 

    T aMax() const {
      T x; size_t k; int e=aMax(x,k); 
      if (e) wblog(FL,"WRN %s() got empty array (e=%d)",FCT,e);
      return x;
    };

    template<class TI>
    T aMax(TI &i, TI &j, const wbarray *B=nullptr) const; 

    double maxDiff (const wbarray &B) const;
    T froNorm2(const wbarray &B) const;
    T norm2() const;
    T norm() const { return Wb::sqrt(norm2()); };

    T hasNorm2(const T &x2) const; 

    T normCol2(size_t k) const;

    T norm2Recs(long i1, long i2) const; 
    T norm2Cols(long j1, long j2) const; 

    wbvector<T>& norm2Cols(wbvector<T> &x2) const;
    wbvector<T> normCols() const;

    T scalarProd(const wbarray &B) const;
    T sum() const;

    wbarray& sum(const WBINDEX &I, wbarray &B) const;

    wbarray& sum(const char *sidx, wbarray &B) const {
       WBINDEX I; Wb::Str2Idx(FL,sidx,I, widx_t(1)); 
       sum(I,B); return B;
    };

    char fitsSize(const wbarray &, char strict=0) const;

    size_t getSize(unsigned i, char lflag=1) const { 
       size_t s=0; 
       if (i<SIZE.len) { s=SIZE[i]; } else
       if (lflag && SIZE.len) {
          for (i=0; i<SIZE.len; ++i) { if (!SIZE[i]) { break; }}
          if (i==SIZE.len) { s=1; }
       }
       return s;
    };

    wbvector<size_t>& getSize( 
       const wbvector<unsigned> &I, wbvector<size_t> &S
     ) const { return SIZE.select(I,S); };

    size_t getSizeM(unsigned r) { 
       size_t M=1; 
       if (r>2 && SIZE.len==r+1) { M=SIZE[r]; }
       else if (SIZE.len!=r) wblog(FL,
          "ERR %s() unexpected rank r=%d/%d",FCT,r,SIZE.len);
       return M;
    };

    bool sameSize(const size_t *s, unsigned n) const;

    bool sameSize(size_t d1, size_t d2) const { 
       size_t s[2]={d1,d2}; return sameSize(s,2);
    };

    template <class Ts>
    bool sameSizeM(unsigned r, unsigned m, const Ts *s) const {
       bool rval=0; unsigned n=r+m;
       if (n==SIZE.len) { unsigned i=0;
          for (; i<m; ++i) { if (SIZE[r+i]!=s[i]) { break; }}
          if (i>=m) { rval=1; } 
       } 
       return rval;
    };

    template<class TB>
    bool sameSize(const wbarray<TB> &, char lflag=2) const;

    bool sameSize1(const size_t *s, unsigned n) const;

    template<class TB> 
    bool hasSameSize(  
       unsigned ia, const wbarray<TB> &B, unsigned ib, char lflag=1
    ) const {
       size_t sa=getSize(ia,lflag), sb=B.getSize(ib,lflag);
       return (sa==sb);
    };

    char sameUptoFac(
      const wbarray &B, T *fac=nullptr, double eps=1e-12
    ) const;

    bool sameAs(const wbarray &B, double eps=1e-12) const;

    bool operator== (const wbarray &B) const { return  sameAs(B); };
    bool operator!= (const wbarray &B) const { return !sameAs(B); };

    T normDiff2(const wbarray &B, char sflag=0) const;
    T normDiff (const wbarray &B, char sflag=0) const {
        return Wb::sqrt(normDiff2(B,sflag)); };

    T normDiff2(const wbvector<T> &B) const;
    T normDiff (const wbvector<T> &B) const { return Wb::sqrt(normDiff2(B)); };

    void TimesEl(const wbarray &, char conj=0);
    wbarray& timesEl(
       const wbarray &B, wbarray &C, T bfac=1, T cfac=0, char conj=0) const;
    wbarray& timesEl_OM(unsigned r,
       const wbarray &B, wbarray &C, char conj=0) const;

    T dotProd(const wbarray &, char conj=0) const;
    T weightedAvg(const wbarray &) const;

    void TensorProd(const wbarray &B,
       const char aflag='N', const char bflag='N', const char kflag=0
    ){
       wbarray X(*this);
       X.tensorProd(B,*this,aflag,bflag,kflag);
    };

    wbarray& tensorProd(
       const wbarray &, wbarray &,
       const char aflag='N', const char bflag='N', const char kflag=0
    ) const;

    wbarray& kron( 
       const wbarray &B, wbarray &X, char aflag='N', char bflag='N'
     ) const { return tensorProd(B,X,aflag,bflag,'k'); };

    wbarray& Kron(const wbarray &B, char aflag='N', char bflag='N') {
       wbarray X; save2(X);
       return X.tensorProd(B,*this,aflag,bflag,'k');
    };

    wbarray& Cat(unsigned dim, wbvector< const wbarray* > &ap);

    wbarray& Cat(unsigned dim, const wbvector< wbarray > &aa) {
       wbvector< const wbarray<T>* > ap(aa.len);
       for (unsigned i=0; i<aa.len; ++i) { ap[i]=&aa[i]; }
       return Cat(dim,ap);
    };

    wbarray& BlockCat(
       const wbarray< wbarray > &aa,
       WBINDEX *D1=nullptr, WBINDEX *D2=nullptr
    ){
       wbarray< const wbarray* > ap(aa.SIZE);
       for (size_t n=aa.numel(), i=0; i<n; ++i) { ap[i]=&aa[i]; }
       return BlockCat(ap,D1,D2);
    };

    wbarray& BlockCat(
       const wbarray< const wbarray* > &ap,
       WBINDEX *D1=nullptr, WBINDEX *D2=nullptr
    );

    void operator+= (T x) {
       if (isRef()) wblog(FL,"ERR %s() reference is considered const",FCT);
       for (size_t n=numel(), i=0; i<n; ++i) { data[i]+=x; }
    };

    void operator-= (T x) {
       if (isRef()) wblog(FL,"ERR %s() reference is considered const",FCT);
       for (size_t n=numel(), i=0; i<n; ++i) { data[i]-=x; }
    };

    template <class T_>
    wbarray& Times(const T_ a, char rcpy=0) { 
       if (a!=T_(1)) {
          if (!rcpy || !isRef()) { return (*this)*=a; }
          else {
             wbarray<T> X; this->times(a,X);
             X.save2(*this);
          }
       }; return *this;
    };

    template <class T_>
    wbarray& times(const T_ a, wbarray<T> &X) const { 
       if (&X==this) { X*=a; } else
       if (a==T_(1)) { X.init(*this); } else
       if (!a) { X.init(SIZE); } else { X.set(*this,a); } 
       return X;
    };

    template <class T_>
    wbarray& operator*=(const T_ a) {
       if (a!=T_(1)) { size_t i=0, n=numel();
          if (isRef()) wblog(FL,"ERR %s() got ref (considered const)",FCT);
          if (!a) { MEM_SET<T>(data,n); } else
          if (a==T_(-1))
               { for (; i<n; ++i) { data[i] = -data[i]; }}
          else { for (; i<n; ++i) { data[i] *= a; }}
       }
       return *this;
    };

    wbarray& operator*=(const wbarray &b) { 
       if (SIZE.len!=2 || b.SIZE.len!=2) wblog(FL,
          "ERR %s() operator*= only supported for matrices (r=%d/%d)",
          FCT, SIZE.len, b.SIZE.len);
       ContractMat(FL,2,b,1);
       return *this;
    };

    template <class T_>
    wbarray& operator/=(const T_ &a) {
       if (SIZE && a!=T_(1)) { size_t i=0, n=numel();
          if (isRef()) wblog(FL,"ERR %s() reference is considered const",FCT);
          if (!a) wblog(FL,"ERR %s() got div/0",FCT);

          if (a==T_(-1)) for (; i<n; ++i) { data[i] = -data[i]; } else
          if (n<4 || WbUtil<T>::isInt()) { for (; i<n; ++i) { data[i]/=a; }}
          else {
             T x=1/T(a);
             if (fabs(double(x*a-1))>1e-12) 
                  { for (; i<n; ++i) { data[i]/=a; }}
             else { for (; i<n; ++i) { data[i]*=x; }}
          }
       }
       return *this;
    };

    void operator+= (const wbarray &B) { Plus(B); };
    void operator-= (const wbarray &B) { Plus(B,-1); };

    wbarray& plus( 
       const wbarray &B, wbarray &C, T bfac=1, char iflag=0) const;

    template<class TB>
    wbarray&  Plus(
       const wbarray<TB> &B, TB bfac=TB(1),
       char iflag=0, T afac=T(1));

    bool equal  (const wbarray&, T eps) const;
    bool unequal(const wbarray&, T eps) const;
    bool equal  (const wbarray&, double eps, double& maxdiff) const;
    bool unequal(const wbarray&, double eps, double& maxdiff) const;

    bool allEqual(const T &x) {
       for (size_t n=numel(), i=0; i<n; ++i) { if (data[i]!=x) return 0; }
       return 1;
    };
    bool allLE(const T &x) {
       for (size_t n=numel(), i=0; i<n; ++i) { if (data[i]>x) return 0; }
       return 1;
    };

    char isVector() const { 
        if (SIZE.len==1) return 1; 
        if (SIZE.len!=2 || (SIZE[0]!=1 && SIZE[1]!=1)) return 0;
        if (SIZE[0]==1 && SIZE[1]==1) return -1; 
        else return 2; 
    };

    unsigned isScalar() const {
       unsigned rval=0;
       if (SIZE.len && SIZE.allEqual(1)) { rval=SIZE.len; }
       return rval;
    };

    bool requiresDataPerm(const wbperm &P, char lflag=3) const;

    void markUnequalZero(const T m=1, const T eps=0) {
        size_t i=0, n=numel(); T z(0);
        if (eps==z)
             { for (; i<n; ++i) data[i] = (data[i]!=T(0)    ? m:z); }
        else { for (; i<n; ++i) data[i] = (Wb::abs(data[i])>eps ? m:z); }
    };

    void toMatrixRef(
        wbarray &A, 
        const UVEC &I, 
        int pos,       
        wbperm &P      
    ) const;

    bool toMatrixRef(const char *F, int L,
        wbarray &A, 
        const UVEC &I, 
        int pos,       
        char &tflag    
    ) const;

    void toMatrixRef(
      wbarray &A, unsigned K, const wbperm &P0 
    ) const;

    void toMatrixRefM( 
      wbarray &A, unsigned K, const wbperm &P0, unsigned R,
      wbperm pfin=wbperm() 
    ) const;

    void toMatrixRef(
      wbarray &A, unsigned i, int pos, wbperm &P
    ) const { return toMatrixRef( A, UVEC(1,&i), pos, P); };

    void toMatrixRef(
       wbarray &A, const char* sidx, int pos, wbperm &P) const {
       WBINDEX S; Wb::Str2Idx(FL,sidx,S,widx_t(1)); 
       return toMatrixRef(A,S,pos,P);
    };

    wbarray& Reshape(const wbvector<size_t> &, bool lflag=0);

    wbarray& Reshape(size_t s1, size_t s2, bool lflag=0) {
       WBINDEX S(2); S.data[0]=s1; S.data[1]=s2;
       return Reshape(S,lflag);
    };
    wbarray& Reshape(size_t s1, size_t s2, size_t s3, bool lflag=0) {
       WBINDEX S(3); S.data[0]=s1; S.data[1]=s2; S.data[2]=s3;
       return Reshape(S,lflag);
    };

    wbarray<T>& Reshape1(const size_t *s, unsigned n) {
       if (!sameSize1(s,n)) wblog(FL,
          "ERR %s() size mismatch (%s <> %s)",FCT,
          SSTR(*this), SSTR(wbvector<size_t>(n,s)));
       SIZE.init(n,s); return *this;
    };
    wbarray& Reshape1(size_t s1, size_t s2) {
       size_t s[2]={s1,s2}; return Reshape1(s,2);
    };
    wbarray& Reshape1(size_t s1, size_t s2, size_t s3) {
       size_t s[3]={s1,s2,s3}; return Reshape1(s,3);
    };

    void GroupIndizes(size_t K); 

    wbarray& groupIndizes_DREF(
       wbarray &A, unsigned k1, unsigned k2
    ) const {
       size_t i, s1=1, s2=1;

       if (k1+k2!=SIZE.len) wblog(FL,
       "ERR invalid index group %d+%d=%d !?",k1,k2,SIZE.len);

       k2=SIZE.len;

       for (i=0; i<k1; ++i) { s1*=SIZE[i]; }
       for (   ; i<k2; ++i) { s2*=SIZE[i]; }

       A.SIZE.init(2); A.SIZE[0]=s1; A.SIZE[1]=s2;

       if (sptr!=A.sptr) { 
          A.INIT2REF(data);
       }

       return A;
    };

    void GroupIndizes(unsigned k1, unsigned k2) {
       GroupIndizes_DREF(*this,k1,k2);
    };

    void groupIndizes_P(
       const WBINDEX&, int, wbperm&,
       size_t* =nullptr, size_t* =nullptr
    ) const;

    wbarray& transpose(const char *F, int L, wbarray &A) const;
    wbarray& Transpose(const char *F=0, int L=0) {
       wbarray X; transpose(F_L,X);
       return X.save2(*this);
    };

    wbarray& flipLR(const char *F, int L, wbarray &A) const;
    wbarray& FlipLR(const char *F=0, int L=0);

    wbarray& flipUD(const char *F, int L, wbarray &A) const;
    wbarray& FlipUD(const char *F=0, int L=0);

    wbarray& MatPermute(const wbperm &P);

    wbarray& Permute(const char *s, int offset=1, char rcpy=0);
    wbarray& Permute(wbperm, char rcpy=0);

    wbarray& permute(wbarray&, const char* s, int offset=1) const;
    wbarray& permute(wbarray&, wbperm) const; 

    wbarray permute(const char* s, int offset=1) const {
       wbarray x; permute(x,s,offset); 
       return x;
    };
    wbarray permute(wbperm &p) const {
       wbarray x; permute(x,p); 
       return x;
    };

    wbarray& select0( 
       const WBPERM &, unsigned dim, wbarray&) const;
    wbarray& Select0(const UVEC &I, unsigned dim) { 
       wbarray A(*this); A.select0(I,dim,*this);
       return *this;
    };

    wbarray& select0(
      size_t i1, size_t i2, unsigned dim, wbarray &Q
    ) const {
       if (dim>=SIZE.len) wblog(FL,
          "ERR dim out of bounds (%d/%d)",dim+1,SIZE.len);
       if (i1>=SIZE[dim] || i2>=SIZE[dim]) wblog(FL,
          "ERR index out of bounds (%d,%d/%d)",i1,i2,SIZE[dim]);

       if (i2<i1) {
          Q.init(); Q.SIZE=SIZE; Q.SIZE[dim]=0;
          return Q;
       }

       wbindex I; I.Index(i1,i2);
       return select0(I,dim,Q);
    };

    wbarray& selectSqueeze(
       wbarray &A, unsigned p,  
       unsigned dim) const; 

    T min() const;
    T max() const;

    T aMin(char zflag=0, size_t *k=nullptr) const;  

    void info(const char* ="") const;
    void info(const char *F, int L,
       const char* ="", unsigned k=-1, unsigned nlt=0, unsigned nlb=0) const;

    wbstring info_mem() const { 
       wbstring s(128);
       snprintf(s.data, s.len,"%s / %p", sptr? STR_(sptr):"(null)", data);
       return s;
    };

    void print(
       const char *F, int L,
       const char *istr="", const char *fmt=""
     ) const;

    void print(const char *istr="", const char *fmt="") const {
       print(0,0,istr,fmt);
    };

    int  printdata(const char *istr="ans",  const char *fmt="") const;
    void print_ref(const char *bstr=" ***", const char *nl="\n") const;

    wbstring toStr() const;

    void contractDiag(
       unsigned ic, const wbarray &B, 
       T a, T b, wbarray &C, char Iflag=0, char bflag='N'
    ) const;

    wbarray& contractVec( const char *F, int L,
       unsigned i1, const wbvector<T> &v, wbarray &B) const; 

    wbarray& ContractVec(
       const char *F, int L, unsigned i1, const wbvector<T>&v) {
       wbarray A; save2(A);
       return A.contractVec(F,L,i1,v,*this);
    };

    template<class TB, class TC>
    wbarray<TC>& contractMat(
       const char *F, int L, unsigned i1,
       const wbarray<TB>&M,  unsigned i2, 
       wbarray<TC>&) const;

    template<class TB, class TC>
    wbarray<TC>& contractMat( 
       const char *F, int L, unsigned i1, 
       const wbarray<TB>&M, wbarray<TC> &C
     ) const { return contractMat(F,L,i1,M,1,C); };

    wbarray& ContractMat(const char *F, int L,
       unsigned i1, const wbarray &M, unsigned i2=1 
    ){ wbarray A; save2(A);
       return A.contractMat(F,L,i1,M,i2,*this);
    };

    template<class TB, class TC>
    wbarray<TC>& contract(const char *F, int L, 
       const ctrIdx&, const wbarray<TB> &B,
       const ctrIdx&, wbarray<TC> &C,
       const wbperm &pfinal=wbperm(),
       T afac=1, TC cfac=1 
    ) const;

    template<class TB, class TC>
    wbarray<TC>& contract(
       const char *F, int L, const iTags &ta,
       const wbarray<TB> &B, const iTags &tb, wbarray<TC> &C, iTags &tc,
       const wbperm &pfinal=wbperm(),
       T afac=1, TC cfac=1 
     ) const;

    template<class TB>
    wbarray<T>& Contract(
       const wbvector<unsigned>&, const wbarray<TB>&,
       const wbvector<unsigned>&,
       const wbperm &pfinal=wbperm()
    );

    template<class TB>
    wbarray<T>& Contract(
       const char*, const wbarray<TB>&, const char*,
       const wbperm &pfinal=wbperm()
    );

    template<class TB, class TC>
    wbarray<TC>& contract(const char *F, const int L,
       const char*, const wbarray<TB>&, const char*, wbarray<TC>&,
       const wbperm &pfinal=wbperm(), T afac=1, TC cfac=1
    ) const;

    template<class TB, class TC>
    wbarray<TC>& contract( 
       char *i1, const wbarray<TB> &B, char *i2, wbarray<TC> &Q,
       const wbperm &pfinal=wbperm()
     ) const { return contract(FL,i1,B,i2,Q,pfinal); };

    template<class TB, class TC>
    wbarray<TC>& contract(
       const wbvector<unsigned> &i1, const wbarray<TB> &B0,
       const wbvector<unsigned> &i2,
       wbarray<TC> &Q0, const wbperm &P
     ) const { return contract(FL,i1,B0,i2,Q0,P); };

    template<class TB, class TC>
    wbarray<TC>& contract(const char *F, int L,
       unsigned i1, const wbarray<TB>&B, unsigned i2, 
       wbarray<TC>&C, const wbperm &pfinal=wbperm()
     ) const {

       if (!i1 || !i2) wblog(F_L,
          "ERR %s() expecting 1-based indices (%d,%d)",FCT,i1,i2);
       --i1; --i2;
       wbvector<unsigned> I1(1,&i1), I2(1,&i2);

       return contract(I1,B,I2,C, pfinal);
    };

    template<class TB, class TC>
    wbarray<TC>& contract(
       unsigned i1, const wbarray<TB>&B, unsigned i2,
       wbarray<TC>&C, const wbperm &pfinal=wbperm()
     ) const { return contract(FL,i1,B,i2,C,pfinal); }

    template<class TB, class TC>
    wbarray<TC>& comm(
      const wbarray<TB>&, wbarray<TC>&, char aflag='N', char bflag='N'
    ) const;

    template<class TB, class TC>
    wbarray<TC>& acomm(
      const wbarray<TB>&, wbarray<TC>&, char aflag='N', char bflag='N'
    ) const;

    wbarray operator+(const wbarray &B) const {
       wbarray C; 
       return plus(B,C);
    };

    wbarray operator-(const wbarray &B) const {
       wbarray C; 
       return plus(B,C,-1);
    };

    wbarray operator*(const wbarray &B) const {
       wbarray C; 
       wbvector<unsigned> ica(1), icb(1); ica[0]=1; icb[0]=0;
       return contract(FL,ica,B,icb,C);
    };

    T trace() const;
    wbarray& trace(unsigned i1, unsigned i2, wbarray &C) const;
    wbarray& trace(const ctrIdx &I1, const ctrIdx &I2, wbarray<T> &C) const;

    wbvector<T>& trace(unsigned k, wbvector<T> &t) const;

    mxArray* toMx() const { return toMx_Struct(); } 
    mxArray* toMx_Struct() const;
    mxArray* toMx_base() const; 

    mxArray* mxCreateStruct(unsigned m, unsigned n) const;
    void add2MxStruct(mxArray *S, unsigned i, char tflag=0) const;

    void put(const char *vname="ans", const char* ws="base") const {
       mxArray *a=toMx();
       mxPutAndDestroy(FL,a,vname,ws);
    };

    void put(const char *F, int L,
       const char *vname="ans", const char* ws="base"
    ) const {
       mxArray *a=toMx();
       wblog(F,L,"I/O putting array '%s' to %s.",vname,ws);
       mxPutAndDestroy(F,L,a,vname,ws);
    };

    void getReal(wbarray<double> &R) const;
    void getImag(wbarray<double> &I) const;

    wbarray& Conj();
    wbarray& ConjTimes(double fac);

    void set(const wbarray<double> &R, const wbarray<double> &I);

    const wbarray& opA(char tflag, wbarray& aux) {
       if (!isMatrix()) wblog(FL,
          "ERR cannot use opA() with rank-%d object.", SIZE.len);

       if (tflag=='N') { return *this; } else
       if (tflag=='T' || tflag=='C') {
          permute(aux, wbperm("21")); if (tflag=='C') aux.Conj();
          return aux;
       }
       wblog(FL, "ERR invalid flag %c<%d>", tflag, tflag);
       return *this;
    };

    bool isSym_aux(
      const char *file, int line, const char *fct,
      const wbarray &B, double eps, double* xref,
      const char symflag='s',
      const char lflag=1
    ) const;

    inline int check_consistency(const char *F=0, int L=0) const; 

    wbvector<size_t> SIZE;  

    Wb::MTYPE mtype;   

    Wb::sptr<T> *sptr; 

    T *data;           

  protected: 

    bool isDiag_aux(const T eps, const char* task) const;
    bool isDiag_aux(
       double *eps, 
       const char* task) const;

    size_t serial_index(const widx_t *I, size_t n) const;

  private:

    void print_rec (UVEC&, size_t, const char*, const char*) const;

    void make2D() {
       if (SIZE.len<2) { SIZE.Resize(2); if (SIZE[0]) SIZE[1]=1; } else
       if (SIZE.len>2) { wblog(FL,"ERR %s() got rank-%d",FCT,SIZE.len); }
    };

    void MemErrMsg(const char* file, int line, size_t s) {
        wblog(file,line, "ERR Out of memory? (%dx%d)", s, sizeof(T));
    };

    inline void DELETE_DATA(char mflag=0);

    inline wbarray& NEW_DATA( 
       size_t n, const T* d0=nullptr, char ref=0, char init=1);

    void RENEW_SPTR() {
       if (sptr) { DELETE_DATA(); } 
       WB_NEW_1(sptr);
    };

    wbarray& INIT2REF(const T *d0, const wbvector<size_t> *S=nullptr) {
       if (S) { SIZE=(*S); }
       size_t n=SIZE.prod(0);

       if (d0) { if (!sptr || !isRef() || d0!=data) {
          RENEW_SPTR();
          data=(sptr->init2ref(n,d0)); 
       }}
       else if (n) wblog(FL,
         "ERR data/size inconsistency: %p / %p => %s",data,d0,SSTR(*this));

       return *this;
    };

    wbarray& NEW(const wbvector<size_t> &S, const T* d0=nullptr) {
       SIZE=S;
       NEW_DATA(SIZE.prod(0),d0); 
       return *this;
    };

    wbarray& NEW(const wbvector<size_t> &S, size_t n, const T* d0=nullptr) {
       SIZE.init(n,S); NEW_DATA(SIZE.prod(0),d0);
       return *this;
    };

    T& dummy() const {
       static T d; memset(&d,0,sizeof(T));
       return d;
    };

}; 

template <class T>
int wbarray<T>::check_consistency(const char *F, int L) const {
   int e=0; 

   if (sptr) {
      if (!data || sptr->data!=data ) { e|=1; }
      if (sptr->check_consistency(0)) { e|=2; }
      if (!SIZE.len || SIZE.min()==0) { e|=4; }

      if (e && F) wblog(F,L,
         "ERR Wb::sptr = %s\n%p -> %p / %p (r=%ld, e=%d)",
          STR_(sptr), sptr, sptr? sptr->data: nullptr, data, SIZE.len, e
      );

      if (mtype!=sptr->mtype) {
      if (mtype || sptr->mtype!=Wb::MEM_REF) wblog(F_L,
         "WRN Wb:sptr mtype inconsistency (%s / %s)",
         Wb::MTYPE_STR[mtype], Wb::MTYPE_STR[sptr->mtype]);
      }
   }
   else {
      if (SIZE.len && SIZE.min()) { e|=8; }
      if (data) { e|=16; }

      if (e && F) wblog(F,L,
         "ERR Wb:sptr inconsistency: %p / %p",sptr,data
      );
   }
   return e;
};

template <class T>
void wbarray<T>::DELETE_DATA(char mflag) {

   if (sptr) { int q=0; check_consistency(FL);
      q=sptr->rm_dref(mflag);
      if (q) { 
         WB_DELETE_1(sptr);  
      }
      else { sptr=nullptr; }    
      data=nullptr;             

      mtype=Wb::MEM_DEF;
   }
   else if (data) wblog(FL,"ERR %s() %p / %p !?",FCT,data,sptr);
};

template <class T>
wbarray<T>& wbarray<T>::NEW_DATA(
   size_t n, const T* d0, char ref, char init) {

   if (!n) {
      DELETE_DATA(); return *this;
   }

   if (!ref) {
      if (!sptr) WB_NEW_1(sptr);

      if (mtype && WBLOG_MMEX) { wblog(FL,
         " |  %s     ... %p -> %s", Wb::MTYPE_STR[mtype],
         this, Wb::size2Str(n*sizeof(T)).data);
      }

      if (d0 || init)
           { data=sptr->malloc_data(n,d0,mtype); }
      else { data=sptr->malloc_base(n, 0,mtype); } 
   }
   else {
      RENEW_SPTR();
      data=(sptr->init2ref(n,d0)); 
   }
   return *this;
};

template <class T>
class wbperm_helper { 

public:

   wbperm_helper() 
    : fac(1), conj(0), rk_(0), rk(0), numel(0), sz(nullptr), stride(nullptr),
      l1(0), l2(0), m_blk(0) {};

   wbperm_helper(const wbvector<widx_t> &sz_, const wbperm &P)
    : fac(1), conj(0), rk_(sz_.len), rk(0), numel(0),
      sz( new widx_t [2*rk_] ), stride(sz+rk_), l1(0), l2(0), m_blk(0)
    {
      unsigned i,j;
      fac=P.fac;  

      if (P.conj && ISCOMPLX_(T)) {
         conj=Wb::conj2bool(P.conj); 
      }

      wbvec<widx_t> Stride_(rk_+1);
      widx_t *stride_=Stride_.data; stride_[0]=1;
      for (i=1; i<=rk_; ++i) { stride_[i] = stride_[i-1]*sz_[i-1]; }
      numel = stride_[i-1];

      if (rk_!=P.len) wblog(FL,
         "ERR %s() rank mismatch (r=%d/%d)",FCT,rk_,P.len);
      if (!numel) { rk=1; return; }

      for (i=0; i<rk_; ++i) { j=P[i];
         sz[i]     = sz_[j];
         stride[i] = stride_[j];
      }

      for (i=1; i<rk_; ++i) {
         if (sz[rk]==1) {
            sz[rk] *= sz[i];
            stride[rk] = stride[i]; 
         }
         else if (sz[i]==1 || stride[i] == stride[rk] * sz[rk]) {
            sz[rk] *= sz[i];
         }
         else if ((++rk)<i) { 
            sz[rk] = sz[i];
            stride[rk] = stride[i];
         }
      }

      if (++rk>=2) { widx_t s0;
         l1=(stride[0]==1 ? 1 : 0);

         s0=(l1 ? sz[0] : 1);
         for (i=0; i<rk; ++i) { if (stride[i]==s0) { l2=i; break;}}

         if (l2<=l1) { 
            wblog(FL,"ERR %s() failed to identify l2=%d "
            "(rk=%d, l1=%d, s0=%ld)",FCT,l2,rk,l1,s0);
         }
      }

      m_blk = wbsys::getCacheLineSize()/sizeof(T); 
   };

   ~wbperm_helper () {
       if (sz) { delete [] sz; sz=0; }
    };

   void permute(const T *src, T *dest, int np=1) const;

   wbperm_helper& save2(wbperm_helper &X) {
      X.rk_=rk_; X.rk=rk; X.numel=numel; rk_=rk=0;
      X.sz=sz;   X.stride=stride;        sz=stride=nullptr;
      X.l1=l1;   X.l2=l2; X.m_blk=m_blk;
      return X;
   };

   mxArray* toMx() const;

   double fac;      
   char conj;       

   unsigned rk_;    
   unsigned rk;     
   widx_t numel;    
   widx_t *sz;      
   widx_t *stride;  

   unsigned l1, l2; 

   unsigned m_blk;
};

template <class T>
void wbarray_permute__(
   T *B_data, const wbarray<T> &A, const wbperm &P, int np=1);

template <class T>
void wbarray_permute__(
   T *B_data, const wbarray<T> &A, const wbperm &P, int np) {

   wbperm_helper<T> PH(A.SIZE,P); 
   PH.permute(A.data,B_data,np);  
};

#ifdef QS_USING_HPTT

template <>
void wbarray_permute__(
   double *B_data, const wbarray<double> &A, const wbperm &P_, int np) {

   wbvector<int> Sz(A.SIZE), P(P_); 

   auto plan = hptt::create_plan(P.data, Sz.len,
       P.fac, A.data, Sz.data, nullptr,  
       0.,    B_data,          nullptr,  
       hptt::ESTIMATE,np>=1 ? np : 1);
   plan->execute();

   if (P.conj) { Conj(); } 
};

template <>
void wbarray_permute__(
   wbcomplex *B_data, const wbarray<wbcomplex> &A, const wbperm &P_, int np) {

   wbvector<int> Sz(A.SIZE), P(P_);
   auto plan = hptt::create_plan(P.data, Sz.len,
       P.fac, (hptt::DoubleComplex*)A.data, Sz.data, nullptr,
       0.,    (hptt::DoubleComplex*)B_data,          nullptr,
       hptt::ESTIMATE,np>=1 ? np : 1);
   plan->execute();

   if (P.conj) { Conj(); } 
};

#endif

class tensor_ref_ { 

  public:
    tensor_ref_() : Ar(nullptr), Az(nullptr) {};
    tensor_ref_(const tensor_ref_ &b) : Ar(b.Ar), Az(b.Az) {};

    tensor_ref_& init(const tensor_ref_ &b) {
       Ar=b.Ar; Az=b.Az;
       return *this;
    };
    tensor_ref_& operator=(const tensor_ref_ &b) { return init(b); };

    bool operator==(const wbarray<double> *A) const {
       return (Ar && Ar==A); };
    bool operator==(const wbarray<wbcomplex> *A) const {
       return (Az && Az==A); };

    tensor_ref_& init(const wbarray<double>    &A) { Ar=&A; Az=nullptr;
       return *this; };
    tensor_ref_& init(const wbarray<wbcomplex> &A) { Az=&A; Ar=nullptr;
       return *this; };

    explicit operator bool() const { return (Ar || Az); };
    bool operator! () const { return (!Ar && !Az); };

    int check(const char *F=nullptr, int L=0) const;
    int rank(unsigned r) const;

    bool sameSize(const wbvector<size_t> &S) const;

    double xflops() const { return (Az ? 2 : 1); };

    wbstring sizeStr() const;
    const wbvector<size_t>* getSIZE() const;

    template<class TC>
    void contract(
       const char *F, int L, const ctrIdx &ia,
       const tensor_ref_ &b, const ctrIdx &ib, wbarray<TC> &C) const;

    const wbarray<double>    *Ar; 
    const wbarray<wbcomplex> *Az; 

  protected:
  private:
};

class tensorRef_ { 
  public:
     tensorRef_ () : id(0), conj(0) {};

     tensorRef_ (const tensorRef_ &B) : id(B.id), conj(B.conj), R(B.R) {
        it=B.it; S=B.S; ID=B.ID; };

     template<class T>
     tensorRef_(const wbarray<T> &A, const iTags *t_=nullptr) {
        init(A,t_); }

     tensorRef_& operator=(const tensorRef_ &B) { return init(B); };

     tensorRef_& init(const tensorRef_ &B) {
        id=B.id; conj=B.conj; R=B.R; it=B.it; S=B.S; ID=B.ID;
        return *this;
     };

     tensorRef_& save2(tensorRef_ &B) {
        B.id=id; B.conj=conj; B.R=R;
        S.save2(B.S); it.save2(B.it); ID.save2(B.ID);
        return B;
     };

     template<class T>
     tensorRef_& init(
        const wbarray<T> &A, const iTags *t_=nullptr, char c=0) {

        R.init(A); 
        S.init2ref(A.SIZE); ID.init(); id=0; conj=c;

        if (!t_) { it.init(); }
        else {
           if (t_->len!=S.len) { unsigned l=t_->len;
              if (l>S.len || (l<=2 && A.numOM(l)!=1)) {
                 wblog(FL,"ERR %s() rank mismatch r=%d/%d (%s) %s",
                 FCT, l,S.len, SSTR(A), STR_(t_));
              }
           }
           it.init(*t_);
        }
        return *this;
     };

     size_t numel() const { return S.prod(0); };

     char check(const char* F=nullptr, int L=0) const;
     char Check(const char* F, int L, unsigned &id_);

     int UpdateItagsCtr(const char *F, int L,
		const ctrIdx &ia, tensorRef_ &b,
		const ctrIdx &ib, const char *tag 
     ) {
        int q=it.UpdateItagsCtr(F,L,ia,b.it,ib,tag);

          conj=ia.conj;
        b.conj=ib.conj;

        return q;
     };

     template<class T>
     int contract(const char *F, int L, 
        const tensorRef_ &b, tensorRef_ &c, double *flops,
        wbarray<T> *C) const;

     int contract(const char *F, int L,
        const tensorRef_ &b, tensorRef_ &c, double &flops
      ) const { return contract(F,L,b,c,&flops, (wbarray<double>*)nullptr); }

     template<class T>
     int contract(const char *F, int L,
        const tensorRef_ &b, tensorRef_ &c, wbarray<T> &C
      ) const { return contract(F,L,b,c,nullptr, &C); }

     wbstring toStr() const {
        wbstring s(ID.len ? 128 : 256); 
        unsigned l, n=s.len;

        l=snprintf(s.data,s.len,"%-16s", SSTR(S)); if (l<n) {
        l+=snprintf(s.data+l,n-l,ID.len ? "%-24s ":"%s",STR(it));  }
        for (unsigned i=0; i<ID.len && l<n; ++i) {
        l+=snprintf(s.data+l,n-l,"%s%d.%d",i?" ":"",ID[i]>>8,ID[i]&255); }

        if (l>=n) wblog(FL,"ERR %s() "
           "string out of bunds (%d/%d)%N%N%s%N%N",FCT,l,n,s.data);
        return s;
     };

     const tensorRef_& print(const char* istr=nullptr, char vflag=1) const {
        char q=(vflag ? check() : 0);
        PRINTF("%s %-50s%s\n", istr? istr:"      ", STR(*this),
           q ? (q&1 ? "empty data" : "length mismatch") : "");
        return *this;
     };

     unsigned id;

     char conj;

     tensor_ref_ R;

     iTags it;
     wbvector<size_t> S; 

     wbvector<unsigned> ID;

  protected:
  private:
};

class tensorRefs : public wbvector<tensorRef_> { 
  public:
     tensorRefs(unsigned l=0) : wbvector<tensorRef_>(l) {};

     char check(const char *F, int L) const {
        char q=0;
        for (unsigned i=0; i<len; ++i) { q|=data[i].check(F,L); }
        return q;
     };

     char Check(const char *F, int L) {
        char q=check(F,L);  
        if (!q) {
           unsigned i=0, id=0;
           for (   ; i<len; ++i) { if (id<data[i].id) { id=data[i].id; }}
           for (i=0; i<len; ++i) { data[i].Check(F_L,id); } 
        }
        return q;
     };

     void print() { print(nullptr,0); } 
     void print(const char *F, int L) const;

     double getOptimalCtrOrder( 
        const char *F, int L, wbperm &P,
        wbperm *pfin=nullptr, char useP=0, char vflag=0);

  protected:
  private:
};

template <class T>
class wbarrRef { 
  public:

    wbarrRef() : dc(0.), fac(1.), tflag(0), D(nullptr) {};

    wbarrRef(const wbarrRef<T> &r) {
       memcpy(this, &r, sizeof(r));
    };

    wbarrRef<T>& operator=(const wbarrRef<T> &r) {
       memcpy(this, &r, sizeof(r));
       return *this;
    };

    void print(const char *vname="ans") {
       printf("\nwbarrRef %s\n", vname);
       printf("... wbarray : 0x%lX (%s%s)\n", D,
          D && D->SIZE.len==1 ? "vec of length ":"", D? SSTR_(D):"");
       printf("... dc      : %g\n", dc);
       printf("... fac     : %g\n", fac);
       printf("... tflag   : %d\n\n", tflag);
    };

    double dc;   
    double fac;  
    bool tflag;  

    wbarray<T> *D;

  protected:
  private:
};

   namespace NR {

template<class T> inline
T pythag(T a, T b) { 

  if (a<0) a=-a;
  if (b<0) b=-b;

  if (a>b)
       { T x=b/a; return a*Wb::sqrt(T(1)+x*x); } else if (b==0) return 0;
  else { T x=a/b; return b*Wb::sqrt(T(1)+x*x); }
};

template<class T> inline
T sign(const T &a, const T &b) {
   return (b>=0 ? (a>=0 ? a : -a) : (a>=0 ? -a : a));
};

}; 

template<class T>
wbarray<T> TestContract(
    const wbarray<T> &A, cUVEC i1,
    const wbarray<T> &B, cUVEC i2, const wbperm &pfinal=wbperm());

template<class T>
void cell2mat(
   const wbMatrix< wbarrRef<T> > &C,
   wbarray<T> &M, WBINDEX &D1, WBINDEX &D2
);

template<class T> inline
size_t wbarray<T>::serial_index(const widx_t *I, size_t len) const {

   size_t i,k, r=SIZE.len, l=r-1, nz=r-len; 

   if (len==0) return -1;
   if (!r || len>r || (len && data==nullptr)) wblog(FL,
      "ERR %s[%s] having %s array",FCT,
      STR(WBINDEX(len,I)+1), SSTR(*this));

#ifdef CHECK_ELEMENT_RANGE
   for (i=0; i<len; i++) if (I[i]>=SIZE[i+nz]) {
       wblog(FL,"ERR %s[%s] index out of bounds (%s)",FCT,
      STR(WBINDEX(len,I)+1), SSTR(*this));
   }
#endif

   for (k=I[l-nz], i=l-1; i<r; i--) {
      k = k*SIZE[i] + ((i>=nz) ? I[i-nz] : 0);
   }

#ifdef CHECK_ELEMENT_RANGE
   if (k>=numel()) wblog(FL,
   "ERR %s() index out of bounds (%d; %s)",FCT,k,SSTR(*this));
#endif

   return k;
};

template<class T>
mxArray* wbarray<T>::toMx_Struct() const { 

   size_t i,dim1,dim2,r=SIZE.len, s=numel();
   WBINDEX I(r);
   T A; 

   if (r==2) { dim1=SIZE[0]; dim2=SIZE[1]; } else
   if (r==1) { dim1=SIZE[0]; dim2=(dim1 ? 1 : 0); } else
   if (r==0) { dim1=dim2=0; }
   else wblog(FL,"ERR %s() not implemented yet for rank>2 (%d)",SIZE.len);

   mxArray *S=A.mxCreateStruct(dim1,dim2);

   for (i=0; i<s; ++i) data[i].add2MxStruct(S,i);

   return S;
};

template<class T>
mxArray* wbarray<T>::mxCreateStruct(unsigned m, unsigned n) const {
   return mxCreateCellMatrix(m,n);
}

template<class T>
void wbarray<T>::add2MxStruct(mxArray *S, unsigned i, char tflag) const {

   if (tflag) { size_t s=0; 
      if (S==nullptr || (s=mxGetNumberOfElements(S))<1 || i>=s) wblog(FL,
      "ERR %s() must follow mxCreateStruct()\n%lx, %d/%d",FCT,S,i+1,s);
   }

   mxSetCell(S,i,toMx());
};

template <class T> inline
mxArray* wbarray<T>::toMx_base() const { 

   mxArray *a=nullptr;

   if (sptr && sptr->mtype) { check_consistency(FL);
      if (sptr->mtype==Wb::MEX_RETURN) {
         size_t l=gP2X.BUF.size();
         a=gP2X.Return(0,0,data); 

         if (a) {
            size_t s1=numel(), s2=mxGetNumberOfElements(a);
            if (s1==s2) {
              if (WBLOG_MMEX) wblog(FL,
                 "2Mx %s (PSX @ %d)",STR_(sptr),l);
               mxSetDimensions(a,SIZE.data,SIZE.len);
            }
            else wblog(FL,
               "TST %s() numel=%ld/%ld %p/%p (PSX @ %ld)",
               FCT,s1,s2,mxGetData(a),data,l
            );
         }
         else {
            wblog(FL,"WRN sptr=%s not / no longer in gP2X @ %ld",STR_(sptr),l);
         }
      }
      else {
         if (sptr->mtype!=Wb::MEM_REF)
         wblog(FL,"TST %s() got %s",FCT,Wb::MTYPE_STR[sptr->mtype]);
      }
   }

   if (!a) { a=cpyRange2Mx(data,SIZE); }

   return a;
};

template <> mxArray* wbarray<unsigned> ::toMx() const { return toMx_base(); }
template <> mxArray* wbarray<int>      ::toMx() const { return toMx_base(); }
template <> mxArray* wbarray<char>     ::toMx() const { return toMx_base(); }
template <> mxArray* wbarray<double>   ::toMx() const { return toMx_base(); }
template <> mxArray* wbarray<wbcomplex>::toMx() const { return toMx_base(); }

template<class T>
wbarray<T>& wbarray<T>::resize(
   const wbvector<size_t> &S, wbarray<T> &B) const {

   if (SIZE.len!=S.len) wblog(FL,
      "ERR length mismatch (%d/%d)",SIZE.len,S.len);
   if (data==B.data) wblog(FL,"WRN %s() got call onto self",FCT);

   if (S==SIZE) { B=(*this); return B; }
   else if (!S) { B.init();  return B; }
   else { B.NEW(S); }

   size_t i,j,k, r=S.len, l=r-1, *s1=SIZE.data, *s2=B.SIZE.data;
   T *b=B.data;
   WBINDEX I(r);

   wbvec<size_t> Smin(r); size_t *s=Smin.data;
   for (i=0; i<r; ++i) { s[i] = MIN( S[i], SIZE[i] ); }

   while (I[l]<s[l]) { 
      for (i=j=I[l], k=l-1; k<r; --k) {
          i = i*s1[k] + I[k];
          j = j*s2[k] + I[k];
      }

      b[j]=data[i];

      ++I[k=0];
      while (I[k]>=s[k] && k<l) { I[k]=0; ++I[++k]; }
   }

   return B;
};

template<class T>
wbarray<T>& wbarray<T>::resize(
   const wbvector<size_t> &S, wbarray<T> &B, double &dx) const {

   if (SIZE.len!=S.len) wblog(FL,
      "ERR length mismatch (%d/%d)",SIZE.len,S.len);
   if (data==B.data) wblog(FL,"WRN %s() got call onto self",FCT);

   if (S==SIZE) { B=(*this); return B; }
   else if (!S) { B.init();  return B; }
   else { B.NEW(S); }

   size_t i,j,k, r=S.len, l=r-1, *s1=SIZE.data, *s2=B.SIZE.data;
   int nx; T x2=0, *b=B.data;
   WBINDEX I(r);

   wbvec<size_t> Smin(r); size_t *s=Smin.data;
   for (i=0; i<r; ++i) { s[i] = MIN( S[i], SIZE[i] ); }

   i=0;

   while (I[l]<s1[l]) { 
      k=l; nx=(I[k]>=s[k] ? 1 : 0);
      for (j=I[k--]; k<r; --k) {
         if (I[k]>=s[k]) { ++nx; } else
         if (!nx) { j = j*s2[k] + I[k]; }
      }
      if (!nx)
           { b[j]=data[i]; }
      else { x2+=Wb::norm2(data[i]); }

      ++I[k=0]; ++i;
      while (I[k]>=s1[k] && k<l) { I[k]=0; ++I[++k]; }
   }

   dx=Wb::sqrt(double(x2));
   return B;
};

template<class T>
wbarray<T>& wbarray<T>::ResizeM(unsigned r, const wbvector<size_t> &Sm) {

   unsigned i=0, m=Sm.len;
   wbvector<size_t> S(SIZE);
   wbarray<T> X;

   if (r+m!=S.len) wblog(FL,
      "ERR %s() size mismatch %d+%d / %d",FCT,r,m,S.len);
   if (int(r)<=2 && Sm.anyUnequal(1)) wblog(FL, 
      "ERR %s() unexpected OM %s",FCT,SSTRM_(r,S.data,m,Sm.data));

   for (; i<m; ++i) { 
      if (S[r+i]>Sm[i]) wblog(FL,"WRN %s() "
         "reducing OM (%d/%d+%d): %d->%d",FCT,i+1,r,m,S[r+i],Sm[i]);
      S[r+i]=Sm[i];
   }
   return resize(S,X).save2(*this);
};

template<class T>
wbarray<T>& wbarray<T>::Enlarge(
   size_t dim,   
   size_t D,     
   size_t istart 
){
   WBINDEX S(SIZE);
   wbarray<T> X;
   wbIndex I(S);

   if (!dim || dim>S.len) wblog(FL,
      "ERR %s() dimension out of bounds (%d/%d)",FCT,dim,S.len);
   if (D<S[--dim]) wblog(FL,
      "ERR %s() envalid new size (%d: %d/%d)",FCT,dim+1,S[dim],D);
   if (istart+S[dim]>=D) wblog(FL,
      "ERR %s() index out of bounds (%d+%d/%d)",FCT,istart+1,S[dim],D);

   S[dim]=D; X.init(S);
   if (istart) {
      for (size_t i=0; ++I; ++i) { 
          I.data[dim]+=istart; X(I)=data[i];
          I.data[dim]-=istart;
      }
   } else { for (size_t i=0; ++I; ++i) X(I)=data[i]; }

   X.save2(*this);
   return *this;
};

template<class T>
wbarray<T>& wbarray<T>::Append2(
   const char *F, int L,
   wbarray<T> &B,
   unsigned dim 
){
   size_t i=0, D=0, ib=0;

   if (B.isEmpty()) { save2(B); return B; }

   if (SIZE.len!=B.SIZE.len) wblog(F_L,
      "ERR %s() incompatible objects (%s <> %s)",FCT,SSTR(*this),SSTR(B));
   if (dim>=SIZE.len) wblog(F_L,
      "ERR %s() dimension out of bounds (%d/%d)",FCT,dim+1,SIZE.len);

   for (i=0; i<SIZE.len; ++i) {
      if (i!=dim) { if (SIZE[i]!=B.SIZE[i]) wblog(F_L,
         "ERR %s() size mismatch (%s <> %s @ %d/%d; %d)",
          FCT,SSTR(*this),SSTR(B),i+1,SIZE.len,dim+1);
      }
      else { D=SIZE[i]+B.SIZE[i]; ib=B.SIZE[i]; }
   }

   WBINDEX S(SIZE); S[dim]=D;
   wbIndex Ia(SIZE), Ib(B.SIZE);
   size_t &ia=Ia.data[dim];
   wbarray<T> X(S);

   for (i=0; ++Ib; ++i) X(Ib)=B.data[i];
   for (i=0; ++Ia; ++i) { 
       ia+=ib; X(Ia)=data[i];
       ia-=ib;
   }

   X.save2(B); init();
   return B;
};

template<class T>
void wbarray<T>::Append2(
   const char *F, int L,
   wbarray<T> &B, unsigned dim, 
   wbarray<double> &W, double eps 
){
   size_t i,n, D=0, ib=0;
   const double *w=W.data;

   if (this->isEmpty()) return;
   if (dim>=SIZE.len) wblog(F_L, 
      "ERR %s() dimension out of bounds (%d/%d)",FCT,dim+1,SIZE.len);

   if (B.isEmpty()) {
      WBINDEX S(SIZE); S[dim]=0;
      B.init(S);
   }

   if (SIZE.len!=B.SIZE.len) wblog(F_L,
      "ERR %s() incompatible objects (%s <> %s)",
       FCT,SSTR(*this),SSTR(B));
   for (i=0; i<SIZE.len; ++i) { if (i!=dim) {
       if (SIZE[i]!=B.SIZE[i]) wblog(F_L,
      "ERR %s() size mismatch (%s <> %s; %d; %d/d)",
       FCT,SSTR(*this),SSTR(B),dim+1,i+1,SIZE.len);
   }}

   wbindex J(SIZE[dim]);
   size_t m=0, *j=J.data; n=W.numel();

   if (!W.isVector() || n!=SIZE[dim]) wblog(F_L,
      "ERR %s() invalid weigths W (%s <> %s @ %d)",
       FCT, SSTR(W), SSTR(*this), dim+1);
   for (i=0; i<n; ++i) { if (w[i]<eps) { j[m++]=i; }}

   if (!m) return; 
   if (m==n) { Append2(F_L,B,dim); return; } 

   J.len=m;
   D=m+B.SIZE[dim]; ib=B.SIZE[dim];

   WBINDEX S(SIZE); S[dim]=m;
   wbIndex Ix(S), Ib(B.SIZE);  S[dim]=D;
   wbindex Ia(S.len),J2;
   wbarray<T> X(S);

   size_t &ia=Ia.data[dim], &ix=Ix.data[dim];

   for (i=0; ++Ib; ++i) { X(Ib)=B.data[i]; }
   for (i=0; ++Ix; ++i) { Ia.setp(Ix.data); ia=j[ia];
       ix+=ib; X(Ix)=(*this)(Ia);
       ix-=ib;
   }

   X.save2(B);
   J.invert(SIZE[dim],J2); Select0(J2,dim);
};

template<class T>
wbarray<T>& wbarray<T>::init(  
   const char *F, int L, const mxArray *a, char ref, char vec) {

   if (!a) { init(); return *this; } 

   unsigned r=mxGetNumberOfDimensions(a);
   const size_t *sz=mxGetDimensions(a);

   if (!mxIsNumChar(a)) wblog(F,L,
      "ERR got rank-%d array %s",r,mxGetClassName(a));

   wbvector<size_t> S(r,sz); init();

   if (vec && r==2) {
      if (S[1]==1) { S.len=1; } else
      if (S[0]==1) { S[0]=S[1]; S.len=1; }
   }

   Mx::Array<T> A(a);

   if (ref && A.data) {
      INIT2REF((T*)A.data,&S);
   }
   else {
      init_bare(S); A.copy_to(data); 
      if (ref && WBLOG_TCAST) wblog(FL,
         "WRN ignoring ref (%s)",A.toStrT().data);
   }

   if (data && !Wb::is_finite(data,1)) wblog(FL,
      "WRN %s() encountered nan or inf in mex input data",FCT);

   return *this;
};

template<class T>
wbarray<T>& wbarray<T>::init(
   const char *F, int L, const mxArray *a, wbperm &P
){
   if (!a || mxIsEmpty(a)) { init(); return *this; }
   if (P.isIdentityPerm()) { return init(F,L,a); }

   unsigned i, r=mxGetNumberOfDimensions(a);
   const size_t *sz=mxGetDimensions(a);
   size_t len=1;

   if (!mxIsNumeric(a)) wblog(F,L,
      "ERR got rank-%d array %s",r,mxGetClassName(a));
   if (!r || !sz || P.len!=r) wblog(FL,
      "ERR %s() got %d/%d @ sz=%p !?",FCT,P.len,r,sz);

   SIZE.init(r);
   for (i=0; i<r; ++i) { SIZE[i]=sz[P[i]]; len*=sz[i]; }

   NEW_DATA(len,nullptr,'\0','\0'); 

   Mx::Array<T>(a).copyTo(data,P);

   return *this;
};

template<class T> inline
void wbarray<T>::initTst() {

   size_t i,j,k, r=SIZE.len, l=r-1, s=numel(), x;
   WBINDEX I(r);

   for (i=0; i<s; i++) {
      for (j=I[l], k=l-1; k<r; k--) j = j*SIZE[k] +I[k];
      for (x=I[0]+1, k=1; k<r; k++) x = x*10 + (I[k]+1);

      data[j]=(T)x;

      k=0; I[0]++;
      while (I[k]>=SIZE[k] && k<l) { I[k]=0; ++I[++k]; }
   }
}

template<class T> inline
wbarray<T>& wbarray<T>::initIdentity(size_t d, char dflag, T dval) {

   if (dflag) { init(d).set(dval); }
   else {
      init(d,d);
      for (size_t i=0; i<d; ++i) {
         data[i+d*i]=dval; 
      }
   }
   return *this;
};

template<class T> inline
wbarray<T>& wbarray<T>::initIdentity(
   const WBINDEX &S, char dflag
){
   size_t i,d,m=S.len/2;
   size_t const* const s=S.data;

   if (S.len%2) wblog(FL,
   "ERR %s() requires even rank object (%s)",FCT,SSTR(*this));
   for (i=0; i<m; i++) if (s[i]!=s[i+m]) wblog(FL,
   "ERR %s() requires symmetric object (%s)",FCT,SSTR(*this));
   if (m>1 && dflag) wblog(FL,
   "ERR %s() rank-2 tensor required with dflag (%s)",FCT,SSTR(*this));

   if (!m) { init(); return *this; }
   if (dflag) { init(s[0]).set(T(1)); return *this; }
   if (m==1 && s[0]==1) { init(S).set(T(1)); return *this; }

   init(S); for (d=s[0], i=1; i<m; i++) d*=s[i];
   for (i=0; i<d; i++) data[i+d*i]=1;

   return *this;
};

template<class T> inline
wbarray<T>& wbarray<T>::Reduce2Id() {
   if (SIZE.len==2) {
      size_t i, j, *s=SIZE.data; T *x=data;

      for (j=0; j<s[1]; ++j) { 
      for (i=0; i<s[0]; ++i, ++x) { *x = (i==j ? 1 : 0); }}
   }
   else if (SIZE.len) wblog(FL,
      "ERR %s() matrix expected (got %s)",FCT,SSTR(*this));

   return *this;
};

template<class T> inline
wbarray<T>& wbarray<T>::initIdentityB(
   size_t d1, size_t d2, size_t k 
){
   init(d1,d2); if (d1 && d2) {
   if (int(k)>=0) {
      if (k<d2) {
         size_t j=k, dx=d1+1; T *x=data+j*d1; 
         for (; j<d2 && (j-k)<d1; ++j, x+=dx) { x[0]=1; }
      }
      else wblog(FL,
      "WRN %s() got all-zero block (%dx%d w/ k=%d)",FCT,d1,d2, k);
   }
   else { k=-k;
      if (k<d1) {
         size_t i=k, dx=d1+1; T *x=data+i; 
         for (; i<d1 && (i-k)<d2; ++i, x+=dx) { x[0]=1; }
      }
      else wblog(FL,
      "WRN %s() got all-zero block (%dx%d w/ k=%d)",FCT,d1,d2,-k);
   }}

   return *this;
};

template<class T> inline
wbarray<T>& wbarray<T>::initIdentityB3(
   size_t d1, size_t d2, size_t D,
   size_t i0, 
   T one
){
   size_t i=0, d12=d1*d2, d12_=d12+1;

   if (i0+d12>D) wblog(FL,"ERR %s() "
      "index out of bounds (%ld*%ld) x %ld @ %ld",FCT,d1,d2,D,i0+1);

   init(d1,d2,D); 
   for (T *x = data + i0*d12; i<d12; ++i, x+=d12_) { (*x)=one; }

   return *this;
};

template<class T> inline
wbarray<T>& wbarray<T>::initPermB(size_t d1, size_t d2, const wbperm &P) {

   unsigned i, j=0; T *x;
   init(d1,d2); x=data;

   for (; j<d2; ++j, x+=d1) { i=P.at_(j);
      if (i>=d1) wblog(FL,"ERR %s() "
         "P out of bounds (j=%d/%d/%d: i=%d/%d)",FCT,j,d2,P.len,i,d1);
      x[i]=1;
   }
   return *this;
};

template<class T> inline
void wbarray<T>::Diag2Vec() {
   size_t i,m;
   WBINDEX S(1);

   if (!isSMatrix()) wblog(FL,
      "ERR Calling %s for rank-%d object (%s).",
       FCT, SIZE.len, SSTR(*this));

   m=S[0]=SIZE[0];
   for (i=1; i<m; i++) data[i]=data[i+i*m]; 

   Resize(S);
};

template<class T> inline
wbarray<T>& wbarray<T>::Reduce2Diag(T eps) {

   if (!SIZE.len) { return *this; }
   else if (!isSMatrix()) wblog(FL,
      "ERR %s() symmetric matrix required\n"
      "(got rank-%d object; %s)",FCT,rank(),SSTR(*this)
   );

   unsigned i, j=0, n=SIZE[0];
   wbarray<T> X(n,1); T *d0=data;

   X.init(n,1); 
   for (; j<n; ++j, d0+=n) { 
      X.data[j]=d0[j];
      for (i=0; i<n; ++i) { if (i!=j) {
         if (Wb::abs(d0[i])>eps) wblog(FL,
            "ERR %s() got off-diagonal matrix element (%g/%g)",
            FCT,double(d0[i]),double(eps)
         );
      }}
   }

   return X.save2(*this);
};

template<class T>
void wbarray<T>::contractDiag(
   unsigned ic, const wbarray<T> &B, 
   T a, T b, wbarray<T> &C, char Iflag, char bflag
) const {

   if (B.SIZE.len!=2) wblog(FL,"ERR invalid rank-%d",SIZE.len);
   if (ic<1 || ic>SIZE.len) wblog(FL,
      "ERR invalid contraction index (ic=%d)",ic);
   if (!strchr("NTCc",bflag)) wblog(FL,
   "ERR invalid flag %c<%d>",bflag,bflag);

   size_t i,k, s=numel(), r=SIZE.len, l=r-1;
   size_t n=B.SIZE.min(), m=B.SIZE[0]; 
   WBINDEX I_(r);
   wbvector<T> h_(n);
   wbarray<T> X;

   size_t *I=I_.data, *S=SIZE.data;
   T *h=h_.data;

   ic--; 
   if (n!=SIZE[ic]) wblog(FL,"ERR dimension mismatch (%d/%d)",n,SIZE[ic]);
   if (n==0) return;

   for (i=0; i<n; i++) h[i]=B.data[i+i*m];
   if (bflag=='C' || bflag=='c') h_.Conj();

   X=*this; if (a!=1) X*=a;

   for (i=0; i<s; i++) {
       if (Iflag)
            X.data[i]/=(b-h[I[ic]]+1e-33); 
       else X.data[i]*=(b-h[I[ic]]);

       k=0; I[0]++; 
       while (I[k]>=S[k] && k<l) { I[k]=0; ++I[++k]; }
   }

   if (C.isEmpty()) X.save2(C); else {
      if (C.SIZE!=SIZE) wblog(FL,"ERR size mismatch in C: %s <> %s",
          SSTR(*this), SSTR(C));
      C+=X;
   }
};

template<class T>
wbarray<T>& wbarray<T>::ExpandDiagonal() {

   size_t i,n=numel();
   wbarray<T> X;

   if (!isVector()) 
      wblog(FL,"ERR %s() got %s array",FCT,SSTR(*this));

   save2(X); init(n,n);
   for (i=0; i<n; ++i) data[i+n*i]=X.data[i];

   return *this;
};

template<class T>   
wbarray<T>& wbarray<T>::ExpandDiagonal(unsigned i1, unsigned i2) {

   if (i1>=SIZE.len || i2>=SIZE.len || i1==i2) wblog(FL,
      "ERR %s() invalid indices (%d %d; %d)",FCT,i1+1,i2+1,SIZE.len);
   if (SIZE[i1]!=1 && SIZE[i2]!=1) wblog(FL,
      "ERR %s() array should be diagonal in (%s; %d %d)",
       FCT, SSTR(*this), i1+1, i2+1);
   if (SIZE[i1]==1 && SIZE[i2]==1) { return *this; }

   if (SIZE[i1]!=1) SWAP(i1,i2); 

   size_t i,j,k,*s0,*s, r=SIZE.len, l=r-1;
   WBINDEX S=SIZE;
   wbarray<T> X;
   wbindex I(r);

   save2(X); s0=X.SIZE.data;
   S[i1]=S[i2]; init(S); s=SIZE.data;

   for (i=0; I[l]<s0[l]; i++) {
       for (j=0,k=l; k<r; k--) {
          if (j) { j*=s[k]; }
          j += (k!=i1 ? I[k] : I[i2]);
       }
       data[j] = X.data[i];

       k=0; I[0]++;
       while (I[k]>=s0[k] && k<l) { I[k]=0; ++I[++k]; }
   }
   return *this;
};

template<class T>
wbarray<T>& wbarray<T>::Expand2Projector(T eps) {

   size_t i,j,m,n=numel();
   WBINDEX S(2);
   wbvector<char> mark(n);

   if (n==0) { if (SIZE.len>1) SIZE[0]=0; return *this; }
   if (!isVector()) wblog(FL,
      "ERR %s() array must be 1-dim vector (got %s)",
       FCT, SSTR(*this)
   );

   for (m=i=0; i<n; i++) if (Wb::abs(data[i])>eps) { mark[i]=1; m++; }

   S[0]=n; S[1]=m; init(S);

   for (i=j=0; i<n; i++)
   if (mark[i]) { data[i+(j++)*n]=1; } 

   return *this;
};

template<class T> inline
T wbarray<T>::froNorm2(const wbarray<T> &B) const {
   T x=0;

   if (!sameSize(B)) wblog(FL,"ERR %s() size mismatch (%s/%s)",
   FCT, SSTR(*this), SSTR(B));

   for (size_t n=numel(), i=0; i<n; i++) x+=Wb::CONJ(data[i])*B.data[i];
   return x;
};

template<class T> inline
T wbarray<T>::norm2() const {
   T d=0;
   for (size_t n=numel(), i=0; i<n; ++i) d+=data[i]*data[i];
   return d;
};

template<> inline
wbcomplex wbarray<wbcomplex>::norm2() const {
   double d=0;

   for (size_t n=numel(), i=0; i<n; ++i)
   d+=data[i].abs2();

   return d;
};

template<class T> inline
T wbarray<T>::hasNorm2(const T &x2) const {
   T q2=0; size_t i=0, n=numel();
   for (; i<n; ++i) {
       if ((q2+=Wb::norm2(data[i]))>x2) { return q2; }
   }
   return -q2;
};

template<class T> inline
T wbarray<T>::normCol2(size_t k) const {

   if (SIZE.len!=2) wblog(FL,
      "ERR %s() invalid matrix %s",FCT,SSTR(*this));
   if (k>=SIZE[1]) wblog(FL,
      "ERR %s() index out of bounds (%d/%d)",FCT,k,SIZE[1]);

   size_t i=0, dim1=SIZE[0];
   const T *d=data+k*dim1; T d2=0; 

   for (i=0; i<dim1; ++i) { d2+=Wb::abs2(d[i]); }
   return d2;
};

template<class T> inline
wbvector<T>& wbarray<T>::norm2Cols(wbvector<T> &x2_) const {

   if (!SIZE.len) { if (x2_.len) { x2_.init(); }; return x2_; }
   if (SIZE.len!=2) wblog(FL,
      "ERR %s() invalid matrix %s",FCT,SSTR(*this));

   x2_.init(SIZE[1]); {
      T* x2=x2_.data; size_t i,j, l=0, dim1=SIZE[0], dim2=SIZE[1];
      for (j=0; j<dim2; ++j) 
      for (i=0; i<dim1; ++i, ++l) { x2[j]+=Wb::norm2(data[l]); }
   }
   return x2_;
};

template<class T> inline
wbvector<T> wbarray<T>::normCols() const {

   wbvector<T> x_; norm2Cols(x_); 
   if (x_) {
      size_t i=0; T *x=x_.data;
      for (; i<x_.len; ++i) { x[i]=Wb::sqrt(x[i]); }
   }
   return x_;
};

template<class T> inline
T wbarray<T>::scalarProd(const wbarray<T> &B) const {
   T x=0; size_t i=0, n=numel();  

   if (!sameSize(B)) wblog(FL,
      "ERR %s() size mismatch: %s <> %s",FCT,SSTR(*this),SSTR(B));
   for (; i<n; i++) { x+=data[i]*B.data[i]; }
   return x;
}

template<> inline
wbcomplex wbarray<wbcomplex>::scalarProd(
  const wbarray<wbcomplex> &B
) const {

   wbcomplex x=0; size_t i=0, n=numel();

   if (!sameSize(B)) wblog(FL,
      "ERR %s() size mismatch: %s <> %s",FCT,SSTR(*this),SSTR(B));
   for (; i<n; ++i) { x+=data[i]*B.data[i].conj(); }

   return x;
};

template<class T> inline
T wbarray<T>::sum() const {
   T x=0;

   for (size_t n=numel(), i=0; i<n; i++)
   x+=data[i];

   return x;
}

template<class T> inline
wbarray<T>& wbarray<T>::sum(const WBINDEX &K, wbarray<T>&A) const {

    size_t i,j,k,q, r=SIZE.len, l=r-1, s=numel();
    WBINDEX S=SIZE;
    wbvector<char> mark(r);
    char *m=mark.data;
    wbindex I;

    for (k=i=0; i<K.len; i++) {
       if (K[i]<SIZE.len) { mark[K[i]]++; k++;
          if (mark[K[i]]>1) wblog(FL,
           "ERR %s() index not unique (%s; %d)",
            FCT,STR(K+1),SIZE.len);
       }
       else wblog(FL,
        "ERR %s() index out of bounds (%s; %d/%d)",
         FCT, STR(K+1), K[i],SIZE.len);
    }

    for (k=i=0; i<mark.len; i++) if (!mark[i]) S[k++]=SIZE[i];
    if (k) S.len=k;
    else { A.init(1,1); A[0]=Wb::addRange(data,s); return A; }

    A.init(S); I.init(SIZE.len);
    for (q=l; q<r; q--) if (!m[q]) break;

    for (i=0; i<s; i++) {
        for (j=I[q], k=q-1; k<r; k--) if (!m[k]) j=j*SIZE[k]+I[k];

        A.data[j]+=data[i];

        k=0; I[0]++; 
        while (I[k]>=SIZE[k] && k<l) { I[k]=0; ++I[++k]; }
    }

    return A;
}

template<class T> inline
T wbarray<T>::trace() const {
   T x=0;

   if (SIZE.len%2) wblog(FL,
      "ERR %s() got object with odd rank %d",FCT,SIZE.len);
   if (!isSquare()) wblog(FL,
      "WRN %s() of non-square array %s",FCT,SSTR(*this));

   if (SIZE.len) {
      size_t j=0, N=SIZE[0], M, dN;
      T *d=data;

      if (SIZE.len==2) { M=SIZE[1]; }
      else {
         unsigned i=1, r2=SIZE.len/2; M=SIZE[r2];
         for (; i<r2; ++i) {
            N*=SIZE[   i];
            M*=SIZE[r2+i]; 
         }
      }

      dN=N+1; if (M>N) { M=N; } 

      for (; j<M; ++j, d+=dN) { x+=(*d); }
   }

   return x;
};

template<class T>
wbvector<T>& wbarray<T>::trace(unsigned k, wbvector <T> &tt) const {

   unsigned r=SIZE.len;

   if (r%2!=1) wblog(FL,
      "ERR %s() requires odd-rank tensor (r=%d; k=%d)",FCT,r,k);

   if (int(k)<0) { k=r-1; } else
   if (k && k+1!=r) wblog(FL,
     "ERR %s() only accepts to be first or last index (k=%d/%d)",FCT,k,r);

   unsigned r2=(r-1)/2;
   const size_t *sz=SIZE.data+(k? 0:1); const T *x=data; T *t;
   size_t i,j, n=SIZE[k], dN, N=1;

   tt.init(n); t=tt.data;

   for (i=0; i<r2; ++i) {
      if (sz[i]!=sz[i+r2]) wblog(FL,"ERR %s() "
         "got non-symmetric tensor (%s; k=%d/%d)",FCT,SSTR(*this),k,r);
      N*=sz[i];
   }

   if (k) { dN=N+1; 
      for (j=0; j<n; ++j, x+=(1-dN), ++t) { 
      for (i=0; i<N; ++i, x+=dN) { (*t)+=(*x); }}
   }
   else { dN=n*(N+1); 
      for (i=0; i<N; ++i, x+=dN) { 
      for (j=0; j<n; ++j) { t[j]+=x[j]; }}
   }

   return tt;
};

#endif

