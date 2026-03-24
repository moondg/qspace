/* ---------------------------------------------------------------------
 * Project : QSpace tensor library (v4.0)
 * Class   : QSpace index routines
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

#ifndef __WB_INDEX_HH__
#define __WB_INDEX_HH__

#define QS_SEP_IFERM  '/'  // 's' 'f' '/' '@' // Wb,Sep19,25

class wbindex : public wbvector<widx_t> { 

  public:

    wbindex(widx_t n=0, char iflag=0) {
       if (iflag) Index(n); else init(n);
    };

    wbindex(const ctrIdx& ic) { init(ic); };

    wbindex(widx_t n, const widx_t* d)
    : wbvector<widx_t>(n,d) {};

    wbindex(const char* s, widx_t offset=0) {
       int i=Wb::Str2Idx(0,0,s,*this,offset);
       if (i<=0) wblog(FL,"ERR %s() invalid index '%s' (%d)",FCT,s,i);
    };

    wbindex& init(widx_t n=0, const widx_t *d=NULL) {
       wbvector<widx_t>::init(n,d);
       return *this;
    };

    wbindex& init(const ctrIdx &ic);

    wbindex& init( 
       const char *F, int L, const mxArray *a, widx_t offset=0);

    widx_t init(
       const char *F, int L, const mxArray *a,
       const char *istr, char check_type=1
    ){ return wbvector<widx_t>::init(F,L,a,istr,check_type); };

    wbindex& init(const mxArray *a, widx_t offset=0) {
       return init(0,0,a,offset);
    };

    wbindex& initGroup(const wbvector<widx_t> &D) { 
       if (!D) { return init(); }
       else {
          widx_t i,d, l=0, ig=0; init(D.sum());
          for (; ig<D.len; ++ig) { d=D[ig];
             for (i=0; i<d; ++i, ++l) { data[l]=ig; }
          }
          return *this;
       }
    };

    wbindex& initGroup(
       const wbvector<widx_t> &P, const wbvector<widx_t> &D
    ){
       widx_t i,d, k=0, ig=0, m=D.sum();
       if (m!=P.len) wblog(FL,
          "ERR %s() invalid group index (%d/%d)",FCT,m,P.len);
       init(P.len+D.len);

       for (; ig<D.len; ++ig) {
          data[k+ig]=d=D[ig]; m=ig+1;
          for (i=0; i<d; ++i, ++k) data[k+m]=P[k];
       }

       return *this;
    };

    wbindex& initGroup0(
       const wbvector<widx_t> &P, const wbvector<widx_t> &D
    ){
       widx_t l=0, i=0; init(D.len);
       for (; i<D.len; ++i) {
          if (l>=P.len) wblog(FL,"ERR %s() index out of bounds",FCT,l,P.len);
          data[i]=P[l]; l+=D[i];
       }
       return *this;
    };

    void initStr(const char* s, widx_t offset=0) {
       int i=Wb::Str2Idx(0,0,s,*this,offset);
       if (i<=0) wblog(FL,"ERR %s() invalid index '%s' (%d)",FCT,s,i);
    };

    bool isIndex(widx_t n=0) { 
       if (n) {
          if (widx_ts(n)>0) {  
             for (widx_t i=0; i<len; ++i) {
             if (data[i]!=(n+i)) return 0; }
          }
          else { 
             n=(1-n-len); if (widx_ts(n)<0) return 0;
             for (widx_t i=0; i<len; ++i) {
             if (data[i]!=(n+i)) return 0; }
          }
       }
       else { 
          for (widx_t i=0; i<len; ++i) { if (data[i]!=i) return 0; }
       }
       return 1;
    };

    wbindex& Index(widx_t n) { 
       init(n); for (widx_t i=0; i<n; ++i) { data[i]=i; }
       return *this;
    };

    wbindex& Index(widx_t i1, widx_t i2) { 
       if (i1>i2) { init(); return *this; }
       init(i2-i1+1);
       for (widx_t i=0; i<len; ++i) data[i]=i+i1;

       return *this;
    };

    wbindex& Index_ex(widx_t n, widx_t ix) {
       widx_t i,k;
       if (ix>=n) wblog(FL,
          "WRN index to be excluded of no relevance (%d/%d)",ix,n);

       init(ix<n ? n-1 : n);
       if (n) { for (k=i=0; i<n; ++i) { if (i!=ix) { data[k++]=i; }}}
       return *this;
    };

    wbindex& Index(widx_t i1, widx_t i2, int p) {
       if (i1>i2) { init(); return *this; }

       widx_t n=i2-i1+1; p=(-p)%int(n);
          if (!p || n==1) { return Index(i1,i2); }
          if (p<0) p+=n; 

       init(n); {
          widx_t i, ix=i1+p, *d=data;
          for (i=ix; i<=i2; ++i, ++d) { (*d)=i; }
          for (i=i1; i< ix; ++i, ++d) { (*d)=i; }
       }
       return *this;
    };

    wbindex& BlockIndex(const wbvector<widx_t> &D, char sub=0) {
       widx_t i,j,d,l; init(D.sum());

       if (sub==0) {
          for (l=i=0; i<D.len; ++i)
          for (d=D[i], j=0; j<d; ++j,++l) data[l]=i; 
       }
       else {
          for (l=i=0; i<D.len; ++i)
          for (d=D[i], j=0; j<d; ++j,++l) data[l]=j; 
       }

       return *this;
    };

    wbindex& BlockIndex(widx_t n, widx_t D) {
       widx_t i,j,l; init(n*D);

       for (l=i=0; i<n; ++i)
       for (  j=0; j<D; ++j,++l) data[l]=i;

       return *this;
    };

    wbindex& BlockIndexI(widx_t n, widx_t D) {
       widx_t i,j,l; init(n*D);

       for (l=j=0; j<D; ++j)
       for (  i=0; i<n; ++i,++l) data[l]=i;

       return *this;
    };

    wbindex& IndexTranspose(widx_t m, widx_t n) {
       widx_t i,j,l;
       init(n*m); if (!n || !m) return *this;

       for (l=j=0; j<m; ++j) {
          data[l++]=j;
          for (i=1; i<n; ++i,++l) data[l]=data[l-1]+m;
       }

       return *this;
    };

    wbindex& operator= (const char* s) {
       int i=Wb::Str2Idx(0,0,s,*this,(widx_t)0);
       if (i<=0) wblog(FL,"ERR %s() invalid index '%s' (%d)",FCT,s,i);
       return *this;
    };

    int extend2Perm(widx_t N, wbperm &P, char toend=0) const;

    void toPerm(wbperm &P) const {
        P.initT(len, data);

        if (P.isValidPerm()<=0) { print("this");
        wblog(FL, "WRN No valid perm."); }
    };

    wbindex& invert(widx_t N, wbindex &I2, char uflag=1) const;

    wbindex& Invert(widx_t N, char uflag=1) {
       wbindex X; save2(X);
       return X.invert(N,*this,uflag);
    };

    wbindex flipIdx(widx_t ndim) const;
    char isUnique(widx_t imax, widx_t offset=0) const;

    mxArray* toMx(const char tflag=0) const {
       return toMx_offset(1,tflag); 
    };

  protected:
  private:

};

wbperm::wbperm (const wbindex &P, double fac_) : fac(fac_) {
   RENEW(P.len,P.data);
   isValidPerm(FL);
};

class wbIndex : public wbvector<widx_t> { 

  public:

    wbIndex() : wbvector<widx_t>() {}; 
    wbIndex(const wbvector<widx_t> &S) { init(S); };

    template<class TI>
    wbIndex(widx_t n, const TI *S) { init(n,S); };

    wbIndex(const wbIndex &I) : wbvector<widx_t>() {
       if (!I.isEmpty()) {
          wbvector<widx_t>::init(I.len); SIZE=I.SIZE;
          if (len) --data[0]; 
       }
    };

    wbIndex& init() {
       wbvector<widx_t>::init(); SIZE.init();
       return *this;
    };

    wbIndex& init(const wbvector<widx_t> &S, bool pp=1) {
       const widx_t* const &s=S.data;
       for (widx_t k=0; k<S.len; ++k) if (!s[k]) {
          wbvector<widx_t>::init(); SIZE.init();
          return *this;
       }

       wbvector<widx_t>::init(S.len); SIZE=S;

       if (len) {
          if (pp) 
               { --data[0]; }               
          else { data[len-1]=SIZE[len-1]; } 
       }
       return *this;
    };

    template<class TI>
    wbIndex& init(unsigned n, const TI *S, bool pp=1) {

       unsigned k=0; for (; k<n; ++k) { if (!S[k]) { break; }}
       if (k<n || (!n && !S)) {
          wbvector<widx_t>::init(); SIZE.init();
          return *this;
       }

       if (n) {
          wbvector<widx_t>::init(n); SIZE.initT(n,S);
          if (!S) { wblog(FL,"ERR %s() got n=%d with null S",FCT,n); }
       }
       else { 
          wbvector<widx_t>::init(1);  
          SIZE.init(1); SIZE[0]=1;    
       }

       if (pp)
            { --data[0]; }               
       else { data[len-1]=SIZE[len-1]; } 

       return *this;
    };

    wbIndex& reset(bool pp=1) {
       if (len!=SIZE.len) wblog(FL,
          "ERR %s() severe size mismatch (%d/%d)",FCT,len,SIZE.len);

       if (len) { wbvector<widx_t>::set(0);
          if (pp)
               { --data[0]; }               
          else { data[len-1]=SIZE[len-1]; } 
       }
       return *this;
    };

    wbIndex& set(wbvector<widx_t> &I) {
       widx_t k=0, *const &i=I.data, *const &s=SIZE.data;
       if (I.len!=len) wblog(FL,
          "ERR %s() size mismatch (%d/%d)",C_FCT,I.len,len);
       for (; k<len; k++) if (i[k]>=s[k]) { wblog(FL,
          "ERR %s() index out of bounds (%d: %d/%d)",C_FCT,k+1,i[k],s[k]);
           data[k]=i[k];
       }
       return *this;
    };

    bool isValid() const {
       if (len!=SIZE.len) return 0;
       for (widx_t i=0; i<len; ++i) { if (data[i]>=SIZE.data[i]) return 0; }
       return 1;
    };

    void checkValid(const char *F, int L) const {
       if (!isValid()) wblog(F_L,
          "ERR %s() invalid hyperindex\n%s",FCT,STR(*this));
    };

    bool operator++() {
       widx_t k=0, l=len-1, *const &s=SIZE.data;

       if (len!=SIZE.len) wblog(FL,
          "ERR wbIndex::++ size mismatch (len=%d/%d)",len,SIZE.len);
       if (!len) { return 0; }

       k=0; ++data[0]; 
       while (data[k]==s[k] && k<l) { data[k]=0; ++data[++k]; }

       if (data[k]>=s[k]) {
          if (k!=l) wblog(FL,"WRN wbIndex::++ "
             "premature stop (%d/%d: %s)",k+1,len,STR(*this));
          return 0;
       }
       return 1;
    };

    bool operator--() {
       widx_t k=0, l=len-1, *const &s=SIZE.data;

       if (len!=SIZE.len) wblog(FL,
          "ERR wbIndex::-- size mismatch (%d/%d)",len,SIZE.len);
       if (!len) { return 0; }

       k=0; data[0]--; 
       while (widx_ts(data[k])==-1 && k<l) { data[k]=s[k]-1; --data[++k]; }

       if (data[k]>=s[k]) {
          if (k!=l) wblog(FL,"WRN wbIndex::-- "
             "premature stop (%d/%d: %s)",k+1,len,STR(*this));
          return 0;
       }
       return 1;
    };

    bool isEmpty() const {
       if (len!=SIZE.len) wblog(FL,"ERR %d/%d",len,SIZE.len);
       if (!len) return 1;
       else {
          const widx_t *s=SIZE.data;
          for (widx_t i=0; i<len; i++) if (!s[i]) return 1;
       }
       return 0;
    };

    bool gotmore() const { 
       if (len!=SIZE.len) wblog(FL,"ERR %d/%d",len,SIZE.len);
       for (widx_t i=0; i<len; ++i) {
          if (data[i]+1<SIZE.data[i]) return 1; else
          if (data[i]>=SIZE.data[i]) wblog(FL,
             "WRN %s() index out of bounds (%s)",FCT,STR(*this)
          );
       }
       return 0;
    };

    bool isdiag() const {
       if (len%2) wblog(FL,"ERR %s() applies to "
          "even-rank objects only (%s)",FCT,SSTR(*this));
       for (widx_t r=len/2, i=0; i<r; ++i) {
          if (data[i]!=data[i+r]) { return 0; }
       }
       return 1;
    };

    widx_t numel() const { return  SIZE.prod(0); };

    widx_t serial() const {
       widx_t l, k=len-1;
       const widx_t *s=SIZE.data;

       if (!len) wblog(FL,"ERR %s() got empty object",C_FCT);
       for (l=data[k--]; k<len; --k) { l = l*s[k] + data[k]; } 
       return l;
    };

    widx_t serial(const widx_t *stride) const {
       widx_t l=0;
          if (!len   ) wblog(FL,"ERR %s() got empty object",C_FCT);
          if (!stride) wblog(FL,"ERR %s() got null strides",C_FCT);
       for (unsigned k=0; k<len; ++k) { l += data[k]*stride[k]; }
       return l;
    };

    wbstring sizeStr() const { return SIZE.toStrf("","x"); };
    wbstring toStr() const {
       char s[128]; snprintf(s,128,"[%s] %s",
          wbvector<widx_t>::toStr().data, SSTR(SIZE));
       return s;
    };

    mxArray* toMx() const {
       const char *fields[]={"S","idx"};
       mxArray *S=mxCreateStructMatrix(1,1,2,fields);

       mxSetFieldByNumber(S,0,0, SIZE.toMx());
       mxSetFieldByNumber(S,0,1,((*this)+1).toMx());
       return S;
    };

    wbvector<widx_t> SIZE;

  protected:
  private:

};

template <class T=widx_t> 
class groupIndex { 

  public:

    groupIndex() {};

    groupIndex(
       const wbvector<T> &P_, const wbvector<T> &D_,
       const wperm_t **p=NULL
    ){ P=P_; D=D_; setup(); if (p) (*p)=P.data; };

    groupIndex(const wbvector<T> &DI, T ng=-1) { init(DI,ng); }

    void init() { P.init(); D.init(); I0.init(); };

    groupIndex& init(
       const wbvector<T> &P_, const wbvector<T> &D_
    ){ P=P_; D=D_; setup(); return *this; };

    groupIndex& initX(
       wbvector<wperm_t> &P_, wbvector<T> &D_,
       const wperm_t **p=NULL
    ){
       P_.save2(P); D_.save2(D); setup(); if (p) (*p)=P.data;
       return *this;
    };

    groupIndex& init(const wbvector<T> &DI, T ng=-1) {
       if (DI.isEmpty()) { init(); return *this; }

       if (long(ng)<0) {
          T n=0,l=0;
          for (; l<DI.len; ++n, l+=(1+DI[l]));
          if (l!=DI.len) wblog(FL,
             "ERR %s() invalid group index (%d/%d)",FCT,l,DI.len);
          ng=n;
       }
       P.init(DI.len-ng); D.init(ng);

       T i,d,k=0,l=0; T *dp=P.data;
       while (l<DI.len) { d=DI.data[l++];
          for (i=0; i<d; ++i, ++k, ++l) { dp[k]=DI.data[l]; }
       }
       if (k!=P.len || l!=DI.len) wblog(FL,
          "ERR %s() %d/%d %d/%d",FCT,k,P.len,l,DI.len);
       setup(); return *this;
    };

    groupIndex& init(T ng) {
       D.init(ng); I0.init(ng+1); P.init();
       return *this;
    };

    groupIndex& init_g(T ig, T nz) {
       if (ig>=D.len || ig+1>=I0.len) wblog(FL,"ERR %s() "
          "index out of bounds (%d/%d,%d)",FCT,ig,D.len,I0.len-1);
       if (ig==0) I0[ig]=0;
       I0[ig+1]=nz; D[ig]=nz-I0[ig];
       return *this;
    };

    groupIndex& setD(T i, T nz) {
       if (i>=D.len) wblog(FL,"ERR %s() "
          "index out of bounds (%d/%d)",FCT,i,D.len);
       D[i]=nz; return *this;
    };

    T init_ic() {
       if (I0.len!=D.len+1) wblog(FL,"ERR %s() "
          "I0 not yet initialized !? (%d/%d+1)",FCT,I0.len,D.len+1);
       T i=0; I0[0]=0;
       for (; i<D.len; ++i) { I0[i+1]=I0[i]+D[i]; }
       return I0[i];
    };

    groupIndex& setup(char full=0) {
       T n=D.cumsum_(I0,'N'); if (n!=P.len) wblog(FL,
          "ERR %s() severe size mismatch (%d/%d)",FCT,n,P.len);
       if (full) { Ig.init(P.len); 
          for (T l=0, j=0, d=0, i=0; i<D.len; ++i, l+=d) {
          for (d=D[i], j=0; j<d; ++j) { Ig[P[l+j]]=i; }}
       }
       return *this;
    };

    T dim(T ig) {
       if (ig>=D.len) wblog(FL,
          "ERR %s() index out of bounds (%d/%d)",FCT,ig,D.len);
       return D.data[ig];
    };

    T numGroups() const { return D.len; };

    void getGroup(T ig, T* &p0, T &d) const {
       if (ig>=D.len || I0.len!=D.len+1) wblog(FL,
          "ERR %s() index out of bounds (%d/%d/%d)",FCT,ig,I0.len-1,D.len);
       p0=P.data+I0.data[ig];
       d=D.data[ig];
    };

    T getGroup(T ig, T &j, T &j2) const {
       if (long(ig)==-1) return 0;
       if (ig>=D.len || I0.len!=D.len+1) wblog(FL,
          "ERR %s() index out of bounds (%d/%d/%d)",FCT,ig,I0.len-1,D.len);
       j =I0.data[ig  ];
       j2=I0.data[ig+1]; return ig+1;
    };

    T getGroup(T ig) const { 
       if (long(ig)==-1) return 0;
       if (ig>=I0.len) wblog(FL,"ERR %s() "
          "index out of bounds (%d/%d/%d)",FCT,ig,I0.len-1,D.len);
       return I0.data[ig];
    };

    template <class T2>
    void getRecs0(const wbMatrix<T2> &, wbMatrix<T2> &, T m) const;

    T groupIter(T &ig, T &j, T &j2) const {
       if (long(ig)>=0) { ++ig; } else { ig=0;
          if (I0.len!=D.len+1) wblog(FL,"ERR %s() "
         "index out of bounds (%d/%d/%d)",FCT,ig,I0.len-1,D.len);
       }
       if (ig<D.len) {
          j =I0.data[ig  ];
          j2=I0.data[ig+1];
       }
       else {
          ig=-1;
       }

       return ig+1;
    };

    mxArray* toMx() const {
       const char *fields[]={"P","D","I0","pp"};
       mxArray *S=mxCreateStructMatrix(1,1,4,fields);

       mxSetFieldByNumber(S,0,0,(P+1).toMx());
       mxSetFieldByNumber(S,0,1, D.toMx());
       mxSetFieldByNumber(S,0,2,(I0+1).toMx());

       wbvector<T> p;
       mxArray *c=mxCreateCellMatrix(1,D.len);

       for (T i=0; i<D.len; ++i) {
          p.init(D[i],P.data+I0[i]);
          mxSetCell(c,i,p.toMx());
       }
       mxSetFieldByNumber(S,0,3,c);

       return S;
    };

     wbvector<T> D; 
     wbperm P;      

     wbvector<T> I0; 

     wbvector<T> Ig;

  protected:
  private:
};

#define IDT unsigned long 
#define ITAG_LEN sizeof(IDT)

#define CC_ITAG '*' 

#define QS_MARK_DUAL '\''
#define IS_CHAR_MARK(c) ((c)>=32 && (c)<= QS_MARK_DUAL)

#define IS_ITAGS_SEP(c) ((c)==',' || (c)==';' || (c)=='|')
#define IS_CHAR_ITAGS(c) ((c)>32 && (c)<127)
#define IS_CHAR_ITAG(c) (((c)>32 && (c)<127) && (c)!=CC_ITAG && !IS_ITAGS_SEP(c))

#define IT2STR(a) (a).itags.toStr().data
#define IT2STR__      itags.toStr().data

class itag_ { 

 public:

   itag_() : t(0) {};
   itag_(const char *F, int L, const char *s) { init(F,L,s); };
   itag_(const char *s) { init(FL,s); };
   itag_(const itag_ &x) : t(x.t) {};

   itag_& init(const char* F, int L, const char *s, unsigned n=-99);

   itag_& init(const char *s=NULL) { return init(FL,s); };

   itag_& init_tag() { 
      char *s = (char*)(&t), c128=char(128);
      unsigned k=0, n=ITAG_LEN; for (; k<n; ++k) { s[k]&=c128; }
      return *this;
   };

   itag_& init_qdir(char c) { t=0;
      if (c=='-' || c==CC_ITAG) Conj(); 
      else if (c!='+') wblog(FL,"ERR %s() invalid c=%c<%d>",FCT,c,c);
      return *this;
   };

   template<class TQ, class TD>
   itag_& init(const char *F, int L,
      wbvector< const QSpace<TQ,TD>* > A, wbindex ia);

#ifndef NOMEX
   itag_& init(const char* F, int L, const mxArray *a);
#endif

   itag_& operator=(unsigned n) { t=n; return *this; }
   itag_& operator=(const char *s) { return init(FL,s); }
   itag_& operator=(const itag_ &x) { t=x.t; return *this; } 

   bool operator==(const itag_ &x) const { return t==x.t; };
   bool operator!=(const itag_ &x) const { return t!=x.t; };

   char sameAs(const itag_ &x, char lflag=0) const; 

   bool isConj(const itag_ &x, char full=1) const { 
      IDT c128=128;
      return (full ? (t^c128)==x.t : (t&c128) ^ (x.t & c128) );
   };

   bool isConj_x(const itag_ &x, char lflag) const; 
   int isValid() const;

   bool isEmpty() const { 
      return (t & ~(IDT)(128) ? 0 : 1); 
   };

   bool isEmpty(char lenient) const;
   bool isEmpty_m() const; 

   char getMark() const { 
       char c=((char&)t) & 127; 
       if (c && !IS_CHAR_MARK(c)) { c=0; }
       return c;
   };

   explicit operator bool() const { return (t!=0); };
   bool operator! () const { return !t; }; 

   explicit operator const char*() const {
      return (const char*)(&t); 
   };

   itag_& Set(const char *s) { char c=isConj();
      init(FL,s);
         if (isConj()) wblog(FL, 
         "ERR %s() unexpected usage (%s; %d)",FCT,STR(*this),c);
      if (c) Conj(); 
      return *this;
   };

   itag_& tSet(const itag_ &x) { 
      if (isConj() ^ x.isConj()) { t=x.t; Conj(); } else { t=x.t; }
      return *this;
   };

   itag_& Set_(const itag_ &x, char keep=0) { 
      if (isEmpty() && !keep) { t=x.t; } else
      if (isConj() ^ x.isConj()) { Conj(); }
      return *this;
   };

   itag_& Conj() {
      ((char&)t) ^= char(128); 
      return *this;
   };

   itag_& setConj() { 
      ((char&)t) |= char(128);  
      return *this;
   };

   bool isConj() const {
      return (((char&)t)<0);
   };

   itag_& deConj() { 
      char &c=((char&)t); if (c<0) { c+=char(128); } 
      return *this;
   };

   itag_ conj()   const { itag_ x(*this); return x.Conj();   } 
   itag_ deconj() const { itag_ x(*this); return x.deConj(); } 

   itag_& SetFlag (unsigned l); 

   itag_& SetFlags(unsigned l=-1, char check=1); 

   itag_& UnsetFlag(unsigned l) { 
      if (l>=ITAG_LEN) wblog(FL,
         "ERR %s() index out of bounds (%d/%d)",FCT,l,ITAG_LEN);
      ((char*)(&t))[l] &= 127; 
      return *this;
   };

   itag_& UnsetFlags() {
      char *s = (char*)(&t);
      for (unsigned l=1; l<ITAG_LEN; ++l) { if (s[l]<0) {
         s[l]+=char(128);  
      }}
      return *this;
   };

   bool GotFlag(unsigned l) const { 
      if (l>=ITAG_LEN) wblog(FL,
         "ERR %s() index out of bounds (%d/%d)",FCT,l,ITAG_LEN);
      return ( ((char*)(&t))[l] < 0 ); 
   };

   unsigned GotFlags(unsigned l=1) const {
      unsigned q=0; 
      const char *s = (char*)(&t);
      for (; l<ITAG_LEN; ++l) { if (s[l]<0) { q|=(1<<l); }}
      return q;
   };

   itag_& SetSI( const char *F, int L, 
      const char *tag, unsigned k,
      char mflag=6, 
      char nfmt=1); 

   itag_& AppendChar(char q, const char *qs=NULL);
   itag_& PrependChar(char q, const char *qs=NULL);
   int CheckFirstChar(char q, const char *qs=NULL);

   itag_& MarkDual(char m=QS_MARK_DUAL, bool always=1);

   itag_& AddTilde() { return AppendChar('~',"KD~"); };
   itag_& SetM() { return AppendChar('~',"KD~"); }; 

   itag_& SetK() {
      if (CheckFirstChar('K',"AEX")>0)
           { return *this; }
      else { return AppendChar('K',"KD~"); }
   };
   itag_& SetD() {
      if (CheckFirstChar('D',"AEX")>0)
           { return *this; }
      else { return AppendChar('D',"KD~"); }
   };

   unsigned to_str(
      const char *F, int L, char *s, unsigned len, char cflag=1) const;

   wbstring toStr(char cflag=1) const { 
      wbstring s(ITAG_LEN+4); 
      to_str(FL,s.data,s.len,cflag);
      return s;
   };

   mxArray* toMx() const { return toStr().toMx(); };

   IDT t;

   static void Reset(unsigned char q=-1) { 
      if (char(q)<0) {
         if (itag_::flag_id & 2) { return; }  
         else { q=4; } 
      }
      if (q&1) wblog(FL,"ERR %s() invalid q=%d",FCT,q);
      if (q<4) { q+=4; } 
      itag_::flag_id=q;
   };

   static void iter_flag_id() {
      #pragma omp critical (iter_flag_id__)
      if ((itag_::flag_id+=4)<4) { itag_::flag_id+=4; }
   };

 protected: 
 private:

   static unsigned char flag_id;
};

   unsigned char itag_::flag_id = 4;

class iTags : public wbvector<itag_> { 

 public:

   iTags(unsigned n=0) : wbvector<itag_>(n) {};
   iTags(const char *F, int L, const char *s) { init(F,L,s); };
   iTags(const char *s) { init(FL,s); };
   iTags(const iTags &b) { wbvector<itag_>::init(b); }; 

   unsigned init(const char* F, int L, const char *s); 
   unsigned init(const char *s) { return init(FL,s); };
   unsigned init(const QDir &qdir);

   iTags& init(unsigned n=0) {
      return (iTags&)wbvector<itag_>::init(n);
   };
   iTags& init(const iTags &b) {
      return (iTags&)wbvector<itag_>::init(b);
   };
   iTags& operator=(const iTags &b) {
      return (iTags&)wbvector<itag_>::init(b);
   };

   iTags& Update(const iTags &b) { 
      if (len!=b.len) { wblog(FL,
         "ERR %s() itag length mismatch (%d/%d) '%s' / '%s'",
          FCT,len,b.len,STR(*this),STR(b));
      }
      for (unsigned i=0; i<len; ++i) {
         if (!data[i].isEmpty()) { if (data[i]!=b.data[i]) {
             wblog(FL,"ERR %s() itag mismatch '%s' / '%s'",
             FCT,STR(*this),STR(b),len,b.len);
         }}
         else if (data[i].isConj() ^ b.data[i].isConj()) {
             wblog(FL,"ERR %s() itag conj mismatch '%s' / '%s'",
             FCT,STR(*this),STR(b),len,b.len);
         }
      }
      return init(b);
   };

   iTags& tSet(const iTags &b) { 
      if (len!=b.len) { wblog(FL,
         "ERR %s() itag length mismatch (%d/%d) '%s' / '%s'",
          FCT,len,b.len,STR(*this),STR(b));
      }
      for (unsigned i=0; i<len; ++i) { data[i].tSet(b.data[i]); }
      return *this;
   }

#ifdef LD_CLEBSCH_QS
   unsigned init_qdir(const char *s);
#endif
#ifndef NOMEX
   unsigned init(const char* F, int L, const mxArray *a); 
#endif

   iTags& init3() {
      wbvector<itag_>::init(3); data[2].Conj();
      return *this;
   };

   iTags& init(size_t l1, const itag_* d1, size_t l2, const itag_* d2) {
      wbvector<itag_>::init(l1,d1,l2,d2);
      return *this;
   };

   unsigned Init(const char* F, int L,
       const char *s,  
       unsigned k=-1,  
       char lflag=-1);

   iTags& Init(unsigned r, int k, const char **st); 

   iTags& initOp(itag_ t, const iTags &b, const char *ot="op");

   iTags& init_alpha(
      unsigned r, unsigned s, const char *tag,
      char t, 
      unsigned k); 

   unsigned Set(const char* F, int L, const char *tag,
      unsigned r=-1,  
      const char *tom=NULL, 
      unsigned m=-1,  
      bool conjOM=1   
   );

   iTags& permute(iTags &B, const wbperm &P) const {
      wbvector<itag_>::permute(B,P('r'));
      if (P.conj%2) { B.Conj(); }
      return B;
   };

   iTags& Permute(const wbperm &P) {
      return iTags(*this).permute(*this,P);
   };

   iTags& SetConj(const char* F, int L, unsigned i) {
      if (len) {
         if (i>=len) wblog(F_L,
            "ERR %s() index out of bounds (%d/%d)",FCT,i,len);
         data[i].Conj();
      }; return *this;
   };

   iTags& Conj() { 
      for (unsigned i=0; i<len; ++i) { data[i].Conj(); }
      return *this;
   };

   iTags& SetFlag(const char* F, int L, unsigned k, unsigned l) {
      if (len) {
         if (k>=len) wblog(F_L,
            "ERR %s() index out of bounds (%d/%d)",FCT,k,len);
         data[k].SetFlag(l);
      }; return *this;
   };

   iTags& SetFlags( 
      const char* F, int L, unsigned k, unsigned l=-1, char check=1) {
      if (len) {
         if (k>=len) wblog(F_L,
            "ERR %s() index out of bounds (%d/%d)",FCT,k,len);
         data[k].SetFlags(l,check);
      }; return *this;
   };

   iTags& UnsetFlag(const char* F, int L, unsigned k, unsigned l) {
      if (len)  {
         if (k>=len) wblog(F_L,
            "ERR %s() index out of bounds (%d/%d)",FCT,k,len);
         data[k].UnsetFlag(l);
      }; return *this;
   };

   iTags& UnsetFlags() {
      for (unsigned k=0; k<len; ++k) {
         data[k].UnsetFlags(); 
      }; return *this;
   };

   iTags& init_tags() {  
      for (unsigned k=0; k<len; ++k) data[k].init_tag();
      return *this;
   };

   bool operator==(const iTags &b) const {
      if (len!=b.len) return 0;
      for (unsigned i=0; i<len; ++i) {
         if (data[i].sameAs(b.data[i],'l')<=0) return 0;
      }
      return 1;
   };
   bool operator!=(const iTags &b) const {
      return !((*this)==b);
   };

   int findRegEx(const char *r) const;
   int RegEx_replace(const char *F, int L, const char *rstr) const;

   int getCtr( 
      const char *F, int L, const iTags &b,
      ctrIdx &ia, ctrIdx &ib, iTags *c=NULL,
      wbvector<char> *ma_=NULL, wbvector<char> *mb_=NULL
    ) const;

   int UpdateItagsCtr(const char *F, int L, 
      const ctrIdx &ia, iTags &b,
      const ctrIdx &ib,
      const char *tag 
   );

   int FlagItagsCtr( 
      const char *F, int L, const ctrIdx &ia,
      iTags &B, const ctrIdx &ib,
      char all=1, 
      char sgn=0  
   );

   bool sameConj(const iTags &b) const { 
      if (len!=b.len) { return 0; }
      for (unsigned i=0; i<len; ++i) {
         if (data[i].isConj(b.data[i],0)) { return 0; }
      }
      return 1;
   };

   bool isOp(int r=-1, char lflag=0) const;

   bool isOpX(unsigned r=-1) const;

   bool isAtensor(unsigned r=-1) const; 

   bool got_qdir(const char *qdir, char fflag=0) const { 
      unsigned i=0, n=(qdir? strlen(qdir) : 0);
      if (n<2) wblog(FL,
         "ERR %s() invalid qdir=%s !?",FCT,qdir?qdir:"(null)");
      if (len<n || (fflag && len!=n)) return 0;

      for (; i<n; ++i) {
         if (qdir[i]=='+') { if ( data[i].isConj()) return 0; } else
         if (qdir[i]=='-') { if (!data[i].isConj()) return 0; }
         else wblog(FL,"ERR %s() invalid qdir=`%s' !?",FCT,qdir);
      }
      return 1;
   };

   unsigned SkipLast(unsigned n=1) { 
      if (n<len) { len-=n; }
      else {
         if (n>len) wblog(FL,
            "ERR %s() range out of bounds (%d/%d) !?",FCT,n,len);
         init();
      }
      return len;
   };

   wbstring toStr(char vflag=1) const;

   const char* toStrk( 
      char *sout, unsigned l, unsigned k,
      const char *st="(itag?)", char cflag=0
    ) const {
      if (k>=len || data[k].isEmpty())
           { snprintf(sout,l,"%s",st?st:""); }
      else { snprintf(sout,l,"%s",data[k].toStr(cflag).data); }
      return sout;
   };

   mxArray* toMx() const {
      if (len) {
         mxArray* a=mxCreateCellMatrix(1,len);
         for (unsigned i=0; i<len; ++i) {
            mxSetCell(a,i,data[i].toMx());
         }
         return a;
      }
      return mxCreateCellMatrix(0,0);
   };

   int got_op_labels(const char *F=NULL, int L=0) const {
      int q=0, x=0; 
      unsigned i=0, l=len/2;
      for (; i<l; ++i) {
         if (!(x=data[i].isValid())) wblog(F_L,
            "ERR %s() got invalid itags (%s)",FCT,STR(*this));
         if (data[i].isConj() || data[i].conj()!=data[i+l]) return 0;
         if (q>0) { if (q!=x) q=-1; } else
         if (!q) { q=x; }
      }
      return q;
   };

   unsigned gotLabels() const { 
      unsigned m=0; 
      for (int q, k=0; k<int(len); ++k) {
         q=data[k].isValid();
         if (q>1) { ++m; } else if (!q) return 0;
      }
      return m;
   };

 protected: 
 private:

   unsigned checkUnique(const char *F=NULL, int L=0) const {
      if (len) {
         unsigned i,j;
         for (i=0; i<len; ++i)
         for (j=i+1; j<len; ++j) if (data[i].t==data[j].t) {
            if (F) wblog(F,L,
               "WRN got non-unique itags (%s)",STR(*this));
            return 0;
         }
      }
      return len; 
   };

};

class iFerm : public MVEC { 
  public:
    iFerm() : r(0) {};
    iFerm(const iFerm &b) : r(0) { init(b); };
    iFerm(const char  *s) : r(0) { init(FL,s); };

    iFerm(unsigned r_, unsigned nsym) : r(r_) { MVEC::init(r_+nsym); };

    iFerm& operator=(const iFerm &b) { return init(b); }

    explicit operator bool () const { return (len || r ? 1 : 0); }

    bool operator==(const iFerm &b) const {
       return (r==b.r && MVEC::operator==(b));
    };
    bool operator!=(const iFerm &b) const { return !(*this==b); };

    iFerm& init() {
       if (data) { MVEC::init(); r=0; }
       return *this;
    };

    iFerm& init(const iFerm &b) {
       MVEC::init(b); r=b.r;
       return *this;
    };

    iFerm& init(unsigned r_, const iFerm &b) {
       unsigned i=0, ns=b.nsym();
       if (int(r_)<0 || ((!r_) ^ (!ns))) wblog(FL,
          "ERR %s() invalid usage (r=%d for nsym=%d)",FCT,r_,ns);

       MVEC::init(r_+ns); r=r_;
       for (; i<ns; ++i) { data[r+i]=b.data[b.r+i]; } 

       return *this;
    };

    iFerm& init(const char *F, int L, const char *s, unsigned r_=-1);
    iFerm& init(const char *F, int L, const mxArray *a, unsigned r_=-1);

    iFerm& init2ref(const iFerm &b) {
       MVEC::init2ref(b); r=b.r;
       return *this;
    };

    iFerm& init_op(const char *F, int L,
       const char *fdr, char f, const iFerm &B);

    iFerm& init_A(const char *F, int L,
       const char *fdr, char fA, char fB, const iFerm &B);

    iFerm& Conj() { 
       if (len && (!r || r>=len)) wblog(FL,
          "ERR %s() invalid fdir='%s'",FCT,STR(*this));
       for (unsigned i=0; i<r; ++i) { data[i]=-data[i]; }
       return *this;
    };

    iFerm& Permute(          const wbperm &P);
    iFerm& permute(iFerm &B, const wbperm &P) const;

    inline int check(
       const char *F, int L, unsigned r_, const iFerm *B=NULL) const;

    inline int check(unsigned r_, const iFerm *B=NULL) const {
       return check(0,0,r_,B); };

    char check(const char *F, int L, const QVec &qtype, unsigned QDIM=-1);

    unsigned nsym() const {
       int n=len-r;   
          if (data) { if (n<=0 || n>9) { n=-10; }}
          else if (r || len) { n=-20; }
          if (n<0) wblog(FL,"ERR %s() invalid `%s'",C_FCT,STR(*this));
       return n;
    };

    int8_t* iSym(unsigned *n_=NULL) const {
       int n=len-r; if (n_) { (*n_)=n; } 
       if ((data && n<=0) || (!data && (n || r))) wblog(FL,
          "ERR %s() invalid object (%s)",C_FCT,STR(*this));
       return (data+r);
    };

    unsigned iSym( 
       wbvector<QType> &qt, wbindex &I, const QVec &qtype, unsigned QDIM=-1
     ) const;

    char get(const char *F, int L, const wbindex &I) const;

    char* toStr(char *s, int n, char compact=-1) const;

    wbstring toStr(char compact=-1) const;

    mxArray* toMx() const {
       return mxCreateString(this->toStr().data);
    };

    unsigned r;

  protected:
  private:

    bool permits_compact() const {
       bool rval=1; unsigned i=0;
       for (; i<r; ++i) { if (abs(data[i])!=1) { return (rval=0); }}
       for (; i<len; ++i) {
          if (data[i]<0 || data[i]>9) wblog(FL,"ERR %s() "
          "invalid iFerm (r=%d/%d; data[%d]=%d)",FCT,r,len,i,data[i]);
       }
       return rval;
    };

   char adapt(char f, char sf) {

      if (sf=='+') { f= 1; } else
      if (sf=='-') { f=-1; } else
      if (sf=='m') { f=-f; } else
      if (sf!='p') { f=-127; } 

      return f;
   };

};

#define CTR_IS_NUM(c)      (((c)>='1' && (c)<='9') || ((c)>='a' && (c)<='z'))

#define CTR_IS_SEP_ANY(c)  (isspace(c) || (c)==',' || (c)==';')
#define CTR_IS_OTHER(c)    (CTR_IS_SEP_ANY(c) || CTR_IS_CONJ(c))

#define CTR_IS_VALID(c)    (CTR_IS_NUM(c) || CTR_IS_OTHER(c))

#define CTR_IS_CONJ(c)     ((c)=='*')
#define IDX_IS_CONJ(c)     ((c)=='*') 

class ctrIdx : public wbvector<unsigned> { 

 public:

    ctrIdx() : conj(0) {}; 

    ctrIdx(unsigned l, unsigned *d, char c=0)
     : wbvector<unsigned>(l,d), conj(0) {
       if (c) { conj=Wb::conj2bool(c); }
    };

    explicit ctrIdx(unsigned i, char c) 
     : wbvector<unsigned>(1), conj(0) { data[0]=i;
       if (c) { conj=Wb::conj2bool(c); }
    };

    ctrIdx(const char *F, int L, const char *s)
     : conj(0) { init(F,L,s); };

    ctrIdx(const char *F, int L, const mxArray *a)
     : conj(0) { init(F,L,a); };

    ctrIdx(const ctrIdx &I)
     : wbvector<unsigned>(I), conj(I.conj), newtags(I.newtags) {};

    ctrIdx(const wbvector<unsigned> &I, char c=0)
     : wbvector<unsigned>(I), conj(Wb::conj2bool(c)) { };

    ctrIdx& init(const ctrIdx &I, char c=0) {
       wbvector<unsigned>::init(I.len,I.data); conj=I.conj;
       if (c) { Wb::conj_add_z2(conj,c); }
       return *this;
    };

    ctrIdx& init(unsigned l=0, unsigned *d=NULL, char c=0) {
       wbvector<unsigned>::init(l,d);
       conj = (c? Wb::conj2bool(c) : 0);
       return *this;
    };

    template <class T>
    ctrIdx& initT(unsigned l=0, T* d=NULL, char c=0) {
       wbvector<unsigned>::initT(l,d);
       conj = (c? Wb::conj2bool(c) : 0);
       return *this;
    };

    ctrIdx& init(const wbvector<unsigned> &I, char c=0) {
       wbvector<unsigned>::init(I);
       conj = (c? Wb::conj2bool(c) : 0);
       return *this;
    };

    ctrIdx& init1(unsigned k, char c=0) { 
       wbvector<unsigned>::init(1); data[0]=k; 
       conj = (c? Wb::conj2bool(c) : 0);
       return *this;
    };

    ctrIdx& init2(unsigned k1, unsigned k2, char c=0) { 
       wbvector<unsigned>::init(2); 
       if (k1==k2) wblog(FL,"ERR %s() got k=(%d,%d) !?",FCT,k1,k2);
       data[0]=k1; data[1]=k2;
       conj = (c? Wb::conj2bool(c) : 0);
       return *this;
    };

   template<class TQ, class TD>
   ctrIdx& initOp(const char *F, int L,
      itag_ t1, const QSpace<TQ,TD> &B, const char *sopt);

    int init(const char *F, int L, const char *s, unsigned offset=1);

    ctrIdx& init(const char *F, int L, const mxArray *a);
    ctrIdx& init(const mxArray *a) { return init(FL,a); }

    ctrIdx& Index(unsigned l, char c=0) { 
       wbvector<unsigned>::init(l);
       conj = (c? Wb::conj2bool(c) : 0);
       for (unsigned i=0; i<l; ++i) { data[i]=i; }
       return *this;
    };

    ctrIdx& operator=(const ctrIdx &I) { return init(I); };

    explicit operator bool() const {
       if (len ) { return 1; }
       if (conj) { return (!conj || conj==1? conj : Wb::conj2bool(conj)); }
       return 0;
    };

    bool operator==(const ctrIdx &I) const {
       return (conj==I.conj && wbvector<unsigned>::operator==(I));
    };

    bool operator!=(const ctrIdx &I) const { 
       return (conj!=I.conj || wbvector<unsigned>::operator!=(I));
    };

    int checkUniqueS(const char *F=0, int L=0, unsigned r=-1) const;

    int checkUnique(
       const char *F=0, int L=0, unsigned r=-1, char sorted=0) const;

    ctrIdx& Sort(wbperm &P, unsigned r=-1) {
       if (len) {
          if (!isSorted())
               { wbvector<unsigned>::Sort(P); }
          else { P.init(); }
          if (int(r)>=0) checkUniqueS(FL,r);
       }
       else { P.init(); }
       return *this;
    };

    ctrIdx& sort(ctrIdx &I, wbperm &P, unsigned r=-1) const {
       if (len) {
          if (!isSorted())
               { wbvector<unsigned>::sort(I,P); }
          else { I.init(*this); P.init(); }
          if (int(r)>=0) I.checkUniqueS(FL,r);
       }
       else { I.init(); P.init(); }
       return I;
    };

    ctrIdx& Sort_(ctrIdx &J, unsigned r=-1) {
       if (len!=J.len) wblog(FL,
          "ERR %s() length mismatch (%s/%s)",FCT,STR(*this),STR(J));
       if (len) {
          if (!isSorted()) {
             wbperm P; UVEC::Sort(P);
             J.UVEC::Permute(P); 
          }
          if (int(r)>=0) { checkUniqueS(FL,r); }
       }
       return *this;
    };

    ctrIdx& invert(unsigned r, ctrIdx &I) const;
    ctrIdx& Invert(unsigned r) { 
       ctrIdx I; save2(I);
       return I.invert(r,*this);
    };

    ctrIdx& Shift_x(unsigned i1, unsigned i2=-1);

    ctrIdx& Conj() {
       Wb::conj_iter_z2(conj); 
       return *this;
    };

    int extend2Perm(widx_t N, wbperm &P, char toend=0) const;

    wbperm& extend2Perm(  unsigned ra, unsigned ma,
       const ctrIdx &icb, unsigned rb, unsigned mb,
       const wbperm &pab, unsigned mc,
       wbperm &Pbc, 
       ctrIdx &kb,  
       ctrIdx &kcb, 
       char Bflag=1 
    ) const;

    ctrIdx& save2(ctrIdx &J) { 
       wbvector<unsigned>::save2(J);
       J.conj=conj; return J;
    };

    mxArray* toMx() const { return toStr().toMx(); };

    wbstring toStr() const { return toStr(0); }; 
    wbstring toStr(char xflag) const;

    char conj;

    iTags newtags; 

 protected: 
 private:

};

bool isCtrIdx(const mxArray *a) {

   unsigned i, r=60; 

   if (Mx::IsDblVector(0,0,a)) { wbvector<double> x(FL,a);
      if (mxGetM(a)!=1 || x.len>r) return 0;

      for (i=0; i<x.len; ++i) {
         if (x[i]<0 || double(int(x[i]))!=x[i] || x[i]>r) return 0;
      }
      if (!x.wbvector<double>::isUnique()) return 0;
      return 1;
   }
   else if (mxIsChar(a)) { ctrIdx idx; wbstring s(a);
      if (idx.init(0,0,s.data)<=0) { return 0; }
      for (i=0; i<idx.len; ++i) { if (idx[i]>r) { return 0; }}
      if (!idx.wbvector<unsigned>::isUnique()) { return 0; }
      return 1;
   }
   return 0;
};

template <class T=SPIDX_T>
class sparseIndex2D { 

  public:

     sparseIndex2D() : rmaj(0) {};

     sparseIndex2D& Setup(char use_rmaj=0);

     T cidx(T j) const;
     T ridx(T ic) const;

     wbstring sizeStr() const { 
        wbstring s; 
        if (S.len) { s=S.toStrf("","x"); } else { s="[]"; }
        return s;
     };

     char rmaj; 
     wbMatrix<T> IJ;
     wbvector<T> S; 

     wbperm P;

     wbvector<T> cIdx;

  protected:
  private:

};

#endif

