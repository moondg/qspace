/* ---------------------------------------------------------------------
 * Project : QSpace tensor library (v4.0)
 * Class   : wbvector (template vector class)
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

#ifndef __WB_SIMPLE_VECTOR_HCC__
#define __WB_SIMPLE_VECTOR_HCC__

template <class T> class wbvector;
template <class T> class wbvec;

template<class T>
class wbvec { 

 public:

   wbvec(unsigned n=0)
    : data(nullptr), l(0), len(0), ref(0) { if (n) {init(n); }};

   wbvec(unsigned n, const T* d0)
    : data(nullptr), l(0), len(0), ref(0) { if (n) { init(n,d0,-1); }}

   wbvec(unsigned n, T* d0, char isr, unsigned l0=-1)
    : data(nullptr), l(0), len(0), ref(0) {
      if (n) {
         if (!isr)
              { init    (n,d0,l0); }
         else { init2ref(n,d0,l0); }
      }
   };

   wbvec(const wbvec &B)
    : data(nullptr), l(0), len(0), ref(0) { 
      if (B.len) { init(B.len,B.data); }
   };

  ~wbvec() {
      if (data && !ref) { free(data); } 
   };

   wbvec& init(unsigned n=0);
   wbvec& init(unsigned n, const T* d0, unsigned l=-1);
   wbvec& init2ref(unsigned n,   T* d0, unsigned l=-1);

   wbvec& Resize(unsigned n);

   const T& operator[] (unsigned i) const { return data[i]; }
         T& operator[] (unsigned i)       { return data[i]; }

   T* current() {
      if (l>len) wblog(FL,"ERR %s() invalid l=%d/%d",FCT,l,len);
      return (data+l);
   };
   const T* current() const {
      if (l>len) wblog(FL,"ERR %s() invalid l=%d/%d",FCT,l,len);
      return (data+l);
   };

   inline int got_left() { return (len-l); };

         T* current(unsigned &n)       { n=(len-l); return current(); };
   const T* current(unsigned &n) const { n=(len-l); return current(); };

   wbvec& operator=(const wbvec& B) { 
      if (B.len) { init(B.len,B.data); } else { init(); }
      return *this;
   };

   wbvec& operator=(const char* s) { 
      unsigned n=(s && s[0] ? strlen(s) : 0);
      if (!n) {
         l=0; if (len) { data[0]=0; }
         return *this;
      }
      if (!len) {
         unsigned q=n, l=0; while (q) { ++l; q>>=1; }
         init(1<<l);
      }

      if (n>=len) {
         if (!len)
              wblog(FL,"ERR %s() called on empty object (%d/%d)",FCT,n,len);
         else wblog(FL,"WRN %s() string out of bounds (%d/%d)",FCT,n,len);
         n=len-1;
      }
      strncat(data,s,n); l=n;
      return *this;
   };

   explicit operator bool() const { return (check()&1); }
   char check() const { char q=0; 
       if (l  ) { q|=1; }
       if (len) { q|=2; }
       if (l>len || ((q!=0) ^ (data!=nullptr))) wblog(FL,"ERR %s() "
          "invalid wbvec (l=%d/%d, %p",FCT,l,len,data);
       return q;
   };

   inline wbvec& reset() {
      l=0; if (len) { data[0]=0; }
      return *this;
   };

   inline char skip_last(T x) { 
      char q=0; if (l && l<=len && data[l-1]==x) { data[--l]=0; q=1; }
      return q;
   };

   inline char skip_last() {
      char q=0; if (l && l<=len) { data[--l]=0; q=1; }
      return q;
   };

   inline int append(T x) { return append(0,0,x); }
   inline int append(const char *F, int L, T x) {
      if (l<len) { data[l]=x; }
      return shift(F,L,1);
   };

   inline int shift(const char *F, int L, int n);

   int catf(const char *F, int L, const char *fmt, ...); 

   inline 
   int cat(const char *s, unsigned n=-1) { return cat(0,0,s,n); }
   int cat(const char *F, int L, const char *s, unsigned n=-1); 

   inline 
   int strcpy(const char *s, int n=-1) { return strcpy(0,0,s,n); }
   int strcpy(const char *F, int L, const char *s, int n=-1);

   int check_bounds(const char *F=nullptr, int L=0, char force=1, int l2=0);

   wbvec& set(T x) {
      if (x==0) { memset(data,0,len*sizeof(T)); }
      else {
         for (unsigned i=0; i<len; ++i) { data[i]=x; }
         if (len && typeid(T)==typeid(char)) { data[len-1]=0; }
      }
      return *this;
   };

   int pad(T x, unsigned n=-1) {
      unsigned l_=l;
      if (int(n)<0) { n=len; } else
      if (n>len) {
         wblog(FL,"WRN %s() length out of bounds (n=%d/%d)",FCT,n,len);
         n=len;
      }
      for (; l<n; ++l) { data[l]=x; }
      return (l-l_);
   };

   wbvec& tr(const char *from, const char *to, unsigned i0=0) {
      unsigned n=(from && *from ? strlen(from) : 0);
      if (n) {
         unsigned m=(to && *to ? strlen(to) : 0);
         if (m) { unsigned i=i0, j;
            for (--m; i<l; ++i) {
            for (j=0; j<n; ++j) {
                if (data[i]==from[j]) { data[i]=to[j<m?j:m]; break; }
            }}
         }
         else { 
            unsigned i=i0, i_=i0-1, j;
            for (   ; i<l; ++i) {
            for (j=0; j<n; ++j) { if (data[i]==from[j]) {
                if (++i_<i) { data[i_]=data[i]; }
                break;
            }}}
            l=(++i_);
         }
      }
      return *this;
   };

   wbvec& tr(const T& from, const T& to, unsigned i0=0) { 
      for (unsigned i=i0; i<l; ++i) { if (data[i]==from) { data[i]=to; }}
      return *this;
   };

   void print(const char *istr=nullptr) const {
      unsigned i=0;        printf("  %4s = [",istr?istr:"ans");
      for (; i<l; ++i) { printf(" %g",double(data[i])); };
      printf("\e[38;5;8m | ");
      for (; i<len; ++i) { printf(" %g",double(data[i])); fflush(0); }
      printf("\e[0m]  (%d / %d elements initialized)\n",l,len);
   };

   T* data;

   unsigned l;   
   unsigned len; 

   char ref;

 protected:
 private:

};

template <class T>
wbvec<T>& wbvec<T>::init(unsigned n) {

   if (ref) { ref=0; } else
   if (data && n!=len) { free(data); } 

   if (n) {
      data = (T*)malloc(n*sizeof(T));  
      if (!data) wblog(FL,
         "ERR %s() failed to allocate data (len=%d)",FCT,n);
      data[0]=0; 
      len=n; l=0;
   }
   else { data=nullptr; l=len=0; }
   return *this;
};

template <class T>
wbvec<T>& wbvec<T>::init(unsigned n, const T* d0, unsigned l0) {
   init(n);
   if (int(l0)<0) { if (n && d0) { l0=n; }} else
   if (l0>n) wblog(FL,
      "ERR %s() size out of bounds (len=%d/%d)",FCT,l0,n);

   if (int(l0)>0) {
      Wb::MemCpy(data,d0,l0); 
      l=l0;
   }
   return *this;
};

template <class T>
wbvec<T>& wbvec<T>::init2ref(unsigned n, T* d0, unsigned l0) {

   if (!d0) wblog(FL,"ERR %s() got null ref",FCT);

   if (data) { init(); }
   if (int(l0)<=0) { l= 0; } else
   if (    l0 < n) { l=l0; } else
   wblog(FL,"ERR %s() size out of bounds (len=%d/%d)",FCT,l0,n);

   len=n; data=d0; ref=1;

   return *this;
};

template <>
wbvec<char>& wbvec<char>::init(unsigned n, const char* d0, unsigned l0) {

   if (int(n)<0) {
      if (int(l0)>=0) { n=l0+1; } 
      else if (d0) {
         l0=strlen(d0); n=l0+1;
      }
      else { n=0; }
   }
   else if (int(l0)<0) {
      if (d0) {
         for (l0=0; l0<n && d0[l0]; ++l0) { }
         if (l0==n && n) { --l0; } 
      }
      else { l0=0; }
   }
   if (!n) { return init(); }

   init(n);

   l=(n<l0 ? n:l0);
   if (l) { Wb::MemCpy(data,d0,l); if (l<n) { data[l]=0; }} 

   return *this;
};

template <>
int wbvec<char>::check_bounds(const char *F, int L, char force, int l2) {
   int n=l;
   if (len) { data[len-1]=0; } 

   if (!l2 ) { l2=len; } else                    
   if (l2<0) { l2=len+l2; if (l2<0) l2=0; } else 
   if (size_t(l2)>len) { l2=len; }                         

   if (l>=size_t(l2)) {
      n=(n>0 ? -n : -1);  

      if (F || L) { char sx[64]; snprintf(sx,64,
         "%s() string out of bounds (l=%d/%d)",FCT,l,l2);
         if (F && force)
              { wblog(F_L,"ERR %s\n%s",sx,data?data:"(null)"); }
         else { wblog(F_L,"WRN %s\n%s",sx,data?data:"(null)"); }
      }
   }
   return n;
};

template <class T>
int wbvec<T>::shift(const char *F, int L, int n) {
   if (unsigned(l+n)<len) { l+=n; data[l]=0; } 
   else {
      if (len) { data[len-1]=0; } 
      if (F || L) { 
         char sx[64]; snprintf(sx,64,
            "%s() wbvec out of bounds (%d%+d/%d)",FCT,l,n,len);

         if (len<60) {
            if (F) wblog(F,L,"ERR %s\n%s%N",sx,data);
            else   wblog(FL, "WRN %s\n%s%N",sx,data);
         }
         else {
            if (F) wblog(F,L,"ERR %s%N%N%s%N",sx,data);
            else   wblog(FL, "WRN %s%N%N%s%N",sx,data);
         }
      }
      if (int(l+n)>=0) { l=len; } else { l=0; if (len) { data[0]=0; }}
      n= (int(n)>0 ? -n : -1); 
   }
   return n;
};

template <>
int wbvec<char>::catf(const char *F, int L, const char *fmt, ...) {
   int n=0;
   if (!l && len) { data[0]=0; }
   if (!fmt || !fmt[0]) { return n; }

   if (l<len) {
      va_list args;  Wb::ARGV wd(&args); 
      va_start(args,fmt);
      n=vsnprintf(data+l,len-l,fmt,args);
      n=shift(F,L,n); 
   }
   else { n=-1; }

   return n;
};

template <>
int wbvec<char>::cat(const char *F, int L, const char *s, unsigned n) {

   if (!l && len) { data[0]=0; }
   if (!s || !s[0] || !len) { return (n=0); }

   unsigned i=0; for (; s[i] && i<n && l<len; ++i, ++l) { data[l]=s[i]; };
   if (l==len && len) { data[len-1]=0; } 
   n=i; l-=n; 

   return shift(F,L,n); 
};

WBINDEX Index(widx_t i1, widx_t i2);

char isUniqueIdxSet(const WBINDEX &I, widx_t imax);

void getIPerm (const WBPERM &P, WBPERM &iP);

   bool mxIsWbvector(
      const char *F, int L, const mxArray *a,
      size_t *n_=nullptr, const char *istr=nullptr, char dflag='d');

   bool mxIsWbvector(const mxArray *a) {
      return mxIsWbvector(0,1,a); 
   };

template <class T>
class wbvector { 

  public:

    wbvector() : data(nullptr), len(0), isref(0) {};

    wbvector(size_t l, T x) : data(nullptr), len(0), isref(0) {
       RENEW(l,nullptr,0,0); 
       set(x);
    };

    wbvector(size_t l, const T* d0=nullptr, char ref=0)
     : data(nullptr), len(0), isref(0) {
       if (ref && d0 && l)
            { init2ref(l,d0); }
       else { RENEW(l,d0); }
    };

    wbvector(const wbvector &v)
     : data(nullptr), len(0), isref(0) { RENEW(v.len, v.data); };

    template <class T2>
    wbvector(const wbvector<T2> &v)
     : data(nullptr), len(0), isref(0) { initT(v.len,v.data); }

    wbvector(const wbvector &v1, const wbvector &v2)
     : data(nullptr), len(0), isref(0) {

       size_t i, k=0;
       RENEW(v1.len+v2.len);

       for (i=0; i<v1.len; ++i, ++k) data[k]=v1[i];
       for (i=0; i<v2.len; ++i, ++k) data[k]=v2[i];
    };

    wbvector(const char *F, int L, const mxArray *a, char check_type=1, char ref=0)
     : data(nullptr), len(0), isref(ref) {
       init(F,L,a,nullptr,check_type,ref);
    };

    wbvector(const mxArray *a, char check_type=1, char ref=0)
     : data(nullptr), len(0), isref(ref) {
       init(FL,a,nullptr,check_type,ref);
    };

    virtual ~wbvector() {
       if (!isref && data) { WB_DELETE(data); }
    };

    wbvector& init_bare(size_t l=0) { return RENEW(l,nullptr,0,0); };

    wbvector& init2val(size_t l, const T& x) {  
       return RENEW_VAL(l,x); };

    wbvector& Index(unsigned n) {
       if (int(n)<0) wblog(FL,"ERR %s() got n=%d",FCT,n); 
       init(n); for (unsigned i=0; i<n; ++i) { data[i]=i; }
       return *this;
    };

    wbvector& Index(int i0, int i2) {
       unsigned i=0, n=i2-i0+1; if (int(n)<0) { n=0; }
       init(n); for (; i<n; ++i) { data[i]=i0+i; }
       return *this;
    };

    wbvector& init(size_t l=0) { return RENEW(l); };
    wbvector& init(size_t l, const T* d0) { return RENEW(l,d0); };

    wbvector& init(size_t l, const T* d0, char ref) {
       if (ref && d0 && l)
            { return init2ref(l,d0); }
       else { return RENEW(l,d0); }
    };

    wbvector& initp(size_t l, const T* d0, size_t n, 
       char init=1 
    ){ return RENEW(l,d0,n,init); };

    wbvector& initStride(size_t l, const T* d0, size_t stride) {
       RENEW(l,nullptr,0,0); 
       for (size_t i=0; i<len; ++i, d0+=stride) { data[i]=*d0; }
       return *this;
    };

    wbvector& init_stride(
       const T* d0, size_t n, size_t m, size_t D0
    ){
       RENEW(n*m);
       Wb::cpyStride(data, d0, n, m, -1, D0);
       return *this;
    };

    wbvector& init(const wbvector &v) {
       RENEW(v.len);
       for (size_t i=0; i<len; ++i) data[i]=v.data[i];
       return *this;
    };

    wbvector& init(size_t n, const wbvector &v) {
       RENEW(n); if (n>v.len) { n=v.len; }
       for (size_t i=0; i<n; ++i) data[i]=v.data[i];
       return *this;
    };

    wbvector& init_data(const T *d0) { 
       for (size_t i=0; i<len; ++i) { data[i]=d0[i]; }
       return *this;
    };

    template<class T0>
    wbvector& init_data(const T0 *d0) {
       for (size_t i=0; i<len; ++i) { data[i]=T(d0[i]); }
       return *this;
    };

    template<class T0>
    wbvector& init_data(const char *F, int L, const T0 *d0) {
       for (size_t i=0; i<len; ++i) { data[i]=T(d0[i]);
          if (T0(data[i])!=d0[i]) wblog(F_L,
             "ERR type conversion changes value (%g->%g) !?",
              double(d0[i]), double(data[i])
          );
       }
       return *this;
    };

    wbvector<T>& init(const wbsparray<T> &S) { return initT(S); };

    template <class T2>
    wbvector<T>& initT(const wbsparray<T2> &S);

    template <class T2>
    wbvector<T>& initT(size_t l, const T2* v) {
       RENEW(l);
       for (size_t i=0; i<len; ++i) data[i]=T(v[i]);
       return *this;
    };

    template <class T2>
    wbvector<T>& initT(const wbvector<T2> &v) {
       RENEW(v.len);
       for (size_t i=0; i<len; ++i) data[i]=T(v.data[i]);
       return *this;
    };

    template <class T2>
    wbvector<T>& initT_(const wbvector<T2> &v) {
       RENEW(v.len);
       for (size_t i=0; i<len; ++i) data[i]=T(v.data[i],'!');
       return *this;
    };

    template <class T2> 
    wbvector<T>& initT(const char *F, int L, const wbvector<T2> &v);

    wbvector& init2ref(const wbvector &v) {
       return init2ref(v.len, v.data); 
    };

    wbvector& init2ref(size_t l, const T* d0) {
       if (l && !d0) wblog(FL,"ERR %s() got null ref (n=%d)",FCT,l);
       if (data && !isref) { RENEW(0); }
       isref=1; len=l; data=(len ? (T*)d0 : nullptr); 
       return *this;
    };

    wbvector& unRef();

    wbvector& init(size_t l1, const T* d1, size_t l2, const T* d2) {
       RENEW(l1+l2); MEM_CPY<T>(data,len,l1,d1,d2);
       return *this;
    };

    int init( 
       const char *F, int L, const mxArray *a, const char *istr=nullptr,
       char check_type QS_UNUSED_VAR =1,
       char ref=0
    ) { return init_Struct(F,L,a,istr,ref); };

    int init_Struct(
       const char *F, int L, const mxArray *a,
       const char *istr=nullptr, char ref=0
    );

    int init_base( 
       const char *F, int L, const mxArray *a,
       const char *istr=nullptr, char check_type=1, char ref=0
    );

    int init(
       const mxArray *a, const char *istr=nullptr, char check_type=1
    ){ return init(0,0,a,istr,check_type); };

    int init_mpfr( 
       const char *F, int L, const mxArray *a, char base=-1,
       const char *istr="", char tcheck=1
    );

    wbvector& operator=(const wbvector &v) {
       if (this!=&v) { RENEW(v.len, v.data); }
       return *this;
    };

    template<class TB>
    wbvector<T>& operator=(const wbvector<TB> &v) {
       init(v.len);
       for (size_t i=0; i<v.len; ++i) { data[i]=T(v.data[i]); }
       return *this;
    };

    wbvector& set(const T &x) { 
       for (size_t i=0; i<len; ++i) { data[i]=x; }
       return *this;
    };

    wbvector& setMin(const T &x) {
       for (size_t i=0; i<len; ++i) { if (data[i]<x) data[i]=x; }
       return *this;
    };
    wbvector& setMax(const T &x) {
       for (size_t i=0; i<len; ++i) { if (data[i]>x) data[i]=x; }
       return *this;
    };

    wbvector& set(const wbvector<size_t> &I, const T &x);
    wbvector& add(const wbvector<size_t> &I, const T &x);

    wbvector& Set(const wbvector<size_t> &I, const wbvector &v);

    void Sort(wbperm &p, char dir=+1) {  Wb::hpsort(*this,p,dir); };
    void sort(wbvector &v2, wbperm &p) const { v2=(*this); Wb::hpsort(v2,p); };

    void Sort(char dir=+1);
    void sort(wbvector &v2) const;

    wbperm& pSort(wbperm &p) const;

    bool isSorted() const { 
       for (size_t i=1; i<len; ++i) { if (data[i-1]> data[i]) return 0; }
       return 1;
    };
    bool isSortedU() const { 
       for (size_t i=1; i<len; ++i) { if (data[i-1]>=data[i]) return 0; }
       return 1;
    };

    bool isUnique() const;

    char isValidPerm(wperm_t n=-1, char f=3) const;

    int contains(const T& x) const { 
       for (size_t i=0; i<len; ++i) { if (data[i]==x) return 1; }
       return 0;
    };

    bool isNormal() const;

    bool isFinite() const {
       return Wb::is_finite(data,len);
    };

    size_t numel() const { return len; }

    size_t find1(const T& x) const {
       for (size_t i=0; i<len; ++i) { if (data[i]==x) return i; }
       return -1;
    };
    size_t find1() const {
       for (size_t i=0; i<len; ++i) { if (data[i]) return i; }
       return -1;
    };

    size_t count(const T& x) const {
       size_t i=0, n=0;
       for (; i<len; ++i) if (data[i]==x) ++n;
       return n;
    };

    size_t icount(const T& x) const { 
       size_t i=0,n=0;
       for (; i<len; ++i) if (data[i]!=x) ++n;
       return n;
    };

    size_t nnz() const { 
       size_t i=0, n=0;
       for (; i<len; ++i) if (data[i]) ++n;
       return n;
    };
    size_t nnz(T eps) const {
        if (eps==0) { return nnz();  }
        else { return len-numZeros(eps); } 
    };

    size_t numZeros(T eps=T(1e-14)) const;
    size_t numEQ(const T&) const;   
    size_t numGT(const T&) const;   
    size_t numLT(const T&) const;   

    wbindex& find(const T &x, wbindex &I, char iflag=0) const;
    wbindex& find(const T &x, wbindex &I, wbindex &Ix) const;
    wbindex& findGT(const T &x, wbindex &I, char iflag=0) const;
    wbindex& findGT(const T &x, wbindex &I, wbindex &Ix) const;

    wbindex& find(wbindex &I, wbindex &Ix) const {
       find(0,Ix,I); return I; }

    wbindex& find(wbindex &I, char nonzero=1) const;

    WBINDEX find() const {
       WBINDEX I; return find(I); }; 

    wbindex& findRange(
       const T &x1, const T &x2, wbindex &I, const char *w=nullptr
    ) const;

    int findClosestSorted(T r) const;
    int findClosest(T r) const;

    size_t findValues(const wbvector &v, wbindex &I) const;

    wbvector& setp(const T *x) { 
       for (size_t i=0; i<len; ++i) { data[i]=x[i]; }
       return *this;
    };

    void set(size_t offset, const T *x, size_t n) {
       if (offset+n>=len) wblog(FL,
          "ERR index out of bounds (%d+%d; %d)",offset,n,len);
       MEM_CPY<T>(data+offset,n,x);
    };

    void reset() { if (len) MEM_SET<T>(data,len); };

    void cpy(const T* v) {
       if (len) MEM_CPY<T>(data,len,v);
    };

    void cpy(const T *x, size_t m) {
       if (m>len) { wblog(FL,"ERR index out of bounds (%d/%d)",m,len); m=len; }
       if (m) MEM_CPY<T>(data,m,x);
    };

    wbvector& save2(wbvector &v) {
       swap(v); if (data!=v.data) { init(); }
       return v;
    };

    void swap(wbvector &v) { 
       if (this!=&v) {
          SWAP(len,   v.len  );
          SWAP(data,  v.data );
          SWAP(isref, v.isref);
       }
    };

    wbvector& swap(size_t i, size_t j) {
       if (i<len && j<len)
            { T x=data[i]; data[i]=data[j]; data[j]=x; }
       else { wblog(FL,
         "ERR %s() index out of bounds (%d,%d/%d)",FCT,i+1,j+1,len); }
       return *this;
    };

    wbvector& Collect(T& A1, T& A2) {
       init(2); A1.save2(data[0]); A2.save2(data[1]);
       return *this;
    };

    wbvector& Collect(T& A1, T& A2, T& A3) {
       init(3); A1.save2(data[0]); A2.save2(data[1]); A3.save2(data[2]);
       return *this;
    };

    wbvector& Skip(size_t);
    wbvector& Skip(const WBINDEX &, char rflag=0);

    wbvector& Move2_FE  (const WBINDEX &I, char FE_flag);
    wbvector& Move2front(const WBINDEX &I) { return Move2_FE(I,'F'); }
    wbvector& Move2end  (const WBINDEX &I) { return Move2_FE(I,'E'); }

    wbvector& Move2front (size_t k) {
       if (k>=len) wblog(FL,"ERR index out of bounds (%d,%d)", k, len);
       T xk=data[k];
       for (size_t i=k; i>0; i--) data[i]=data[i-1];
       data[0]=xk; return *this;
    }

    wbvector& Move2end (size_t k) {
       if (k>=len) wblog(FL,"ERR index out of bounds (%d,%d)", k, len);
       T xk=data[k];
       for (size_t i=k+1; i<len; i--) data[i-1]=data[i];
       data[len-1]=xk; return *this;
    }

    wbvector& mvLastTo (size_t k) { 
       if (k>=len) wblog(FL,"ERR index out of bounds (%d,%d)", k, len);
       T xk=data[len-1];
       for (size_t i=len-1; i>k; i--) data[i]=data[i-1];
       data[k]=xk; return *this;
    }

    wbvector& mvFirstTo (size_t k) { 
       if (k>=len) wblog(FL,"ERR index out of bounds (%d,%d)", k, len);
       T xk=data[0];
       for (size_t i=0; i<k; i++) data[i]=data[i+1];
       data[k]=xk; return *this;
    }

    bool checkSameLength(
       const wbvector &a,
       const char *F=nullptr, int L=0, const char *istr=""
    )  const;

    bool checkSameLength(
       const wbvector &a,
       const wbvector &b,
       const char *F=nullptr, int L=0, const char *istr=""
    )  const;

    template <class T_>
    void operator+= (const T_ &a) {
       if (a!=T(0)) { for (size_t i=0; i<len; ++i) { data[i]+=a; }}
    };

    template <class T_>
    void operator-= (const T_ &a) {
       if (a!=T(0)) { for (size_t i=0; i<len; ++i) { data[i]-=a; }}
    };

    template <class T_>
    wbvector& operator*= (const T_ &a) {
       if (len && a!=T_(1)) {
          if (isref) wblog(FL,"ERR %s() reference is considered const",FCT);
          if (!a) { MEM_SET<T>(data,len); }
          else { size_t i=0;
             if (a==T_(-1))
                  { for (; i<len; ++i) { data[i]=-data[i]; }}
             else { for (; i<len; ++i) { data[i]*=a; }}
          }
       }
       return *this;
    };

    wbvector& Times(double fac, char conj=0) {
       if (len) {
          Wb::rangeTimes(data,len,fac,conj);
       }
       return *this;
    };

    wbvector& times(double fac, wbvector &B, char conj=0) const {
       if (!len || !fac) { B.init(len); }
       else {
          B.RENEW(len,nullptr,0,0); 
          Wb::rangeTimes(B.data,data,len,fac,conj);
       }
       return B;
    };

    template <class T_>
    wbvector& operator/= (const T_ &a) {
       if (len && a!=T_(1)) { size_t i=0;
          if (isref) wblog(FL,"ERR %s() reference is considered const",FCT);
          if (!a) wblog(FL,"ERR got div/0");

          if (a==T_(-1)) { for (; i<len; ++i) data[i]=-data[i]; } else
          if (len<4 || WbUtil<T>::isInt()) { for (; i<len; ++i) data[i]/=a; }
          else {
             T x=1/T(a);
             if (fabs(double(x*a-1))>1e-12) 
                  { for (; i<len; ++i) { data[i]/=a; }}
             else { for (; i<len; ++i) { data[i]*=x; }}
          }
       }
       return *this;
    };

    template <class T_>
    wbvector operator* (const T_ &a) const {
       wbvector u(*this); if (a!=1) u*=a; 
       return u;
    };

    template <class T_>
    wbvector operator+ (const T_ &a) const {
       wbvector u(*this); 
       if (a!=0) for (size_t i=0; i<u.len; ++i) { u.data[i]+=a; }
       return u;
    };

    template <class T_>
    wbvector operator- (const T_ &a) const {
       wbvector u(*this); 
       if (a!=0) for (size_t i=0; i<u.len; ++i) { u.data[i]-=a; }
       return u;
    };

    wbvector operator+ (const wbvector &b) const {
       wbvector u(*this); checkSameLength(b,FLF); 
       for (size_t i=0; i<u.len; i++) u.data[i]+=b.data[i];
       return u;
    };

    wbvector operator- (const wbvector &b) const {
       wbvector u(*this); checkSameLength(b,FLF); 
       for (size_t i=0; i<u.len; i++) u.data[i]-=b.data[i];
       return u;
    };

    void operator+= (const wbvector &b) {
       checkSameLength(b,FLF);
       for (size_t i=0; i<len; ++i) data[i]+=b.data[i];
    };

    void operator-= (const wbvector &b) {
       checkSameLength(b,FLF);
       for (size_t i=0; i<len; ++i) data[i]-=b.data[i];
    };

    wbvector& Plus (const T* b) { 
       for (size_t i=0; i<len; ++i) data[i]+=b[i];
       return *this;
    };

    wbvector& Plus (const T* b, T bfac) {
       size_t i=0; 
       if (bfac==+1) for (; i<len; ++i) data[i]+=b[i]; else
       if (bfac==-1) for (; i<len; ++i) data[i]-=b[i]; else
       if (bfac!= 0) for (; i<len; ++i) data[i]+=(bfac*b[i]);
       return *this;
    };

    wbvector& Plus (const wbvector &b, T bfac=1) {
       checkSameLength(b,FLF);
       return Plus(b.data,bfac);
    };

    void plus (
       const wbvector &b, wbvector &c, T bfac=1, T cfac=0
    ) const {
       if (cfac) {
              checkSameLength(b,c,FLF); c*=cfac; c.Plus(*this); }
       else { checkSameLength(b,  FLF); c=*this; }

       c.Plus(b,bfac);
    };

    wbvector& RevertSigns();

    int isum() const {
       int s=0; if (len) {
          size_t i=1; s=data[0]; for (; i<len; ++i) { s+=data[i]; }
       }; return s;
    };

    T sum() const {
       T s=T(); if (len) {
          size_t i=1; s=data[0]; for (; i<len; ++i) { s+=data[i]; }
       }; return s;
    };

    T sum(size_t i1, const char *estr) const {
       if (!estr || strcmp(estr,"end")) wblog(FL,
          "ERR invalid usage (got `%s')", estr);
       return sum(i1,len-1);
    };

    T sum(size_t i1, size_t i2) const {
       T s=T();

       if (!len) { if (!i1 && !i2) return s; }
       if (i1>=len || i2>=len) wblog(FL,
          "ERR index out of bounds (%d,%d; %d)",i1,i2,len);
       if (i1<=i2) {
          size_t i=i1+1; s=data[i1]; for (; i<=i2; ++i) { s+=data[i]; }
       }
       return s;
    };

    T cumsum_(wbvector &S, char xflag=0) const;
    T Cumsum(char zflag=0) const;

    template <class T2> 
    T cumsum0prod(wbvector<T2> &b, wbvector<T> &cs, char xflag) const;

    wbvector& Eldiv( const wbvector &v) const;

    wbvector& tensorProd(const wbvector &v2, wbvector &vv) const;

    wbvector& TensorProd(const wbvector &v2) { 
       wbvector v1; save2(v1);
       return v1.tensorProd(v2,*this);
    };

    T dotProd(const wbvector &v) const { T x=T();
       if (len!=v.len) wblog(FL,
          "ERR %s() length mismatch (len=%d/%d)",FCT,len,v.len);
       for (size_t i=0; i<len; ++i) { x+=(data[i]*v[i]); }
       return x;
    };

    T scalarProd(const wbvector &v) const { 
       if (len!=v.len) wblog(FL,
          "ERR %s() length mismatch (len=%d/%d)",FCT,len,v.len);

       if (len) return scalarProd(v.data);
       else return 0;
    };

    T scalarProd(const T *v) const { T x=T();
       for (size_t i=0; i<len; ++i) { x+=(data[i]*Wb::CONJ(v[i])); }
       return x;
    };

    T norm() const { return Wb::sqrt(norm2()); }

    T norm2() const { 
       T x=0; if (!len) { return x; }
       x=Wb::norm2(data[0]);

       for (size_t i=1; i<len; ++i) { x+=Wb::norm2(data[i]); }
       return x;
    };

    T normDiff2() const { 
       if (!len) { return 0; }
       else {
          size_t i=1; T dx2=0;
          for (; i<len; ++i) { dx2+=Wb::norm2(data[i]-data[i-1]); }
          return (dx2/len);
       }
    };

    T normDiff() const { return Wb::sqrt(normDiff2()); }

    double normDiff(const wbvector &v, T fac=1) const {
       return std::sqrt(double(normDiff2(v,fac)));
    };

    T normDiff2(const wbvector &v, T fac=1) const;
    double normDiff2_(const wbvector &v, double fac=1) const;

    T GSProject(const T *v, char isnorm=0, char tflag=0) const {
       return gs_project_range(data,v,len,1,isnorm,tflag);
    };

    wbvector& Normalize(const char *F=nullptr, int L=0) {
       double n=std::sqrt(norm2());
       if (n==0) {
          if (F) wblog(F,L,"WRN normalize zero vector!? - randomize");
          setRand(1,-0.5); return Normalize(F,L);
       }
       (*this)*=(1./n); return *this;
    };

    wbvector& add(const T* d) {
       for (size_t i=0; i<len; ++i) data[i]+=d[i];
       return *this;
    };

    wbvector& subtract(const T* d) {
       for (size_t i=0; i<len; ++i) data[i]-=d[i];
       return *this;
    };

    void getReal(wbvector<double> &R) const;
    void getImag(wbvector<double> &I) const;
    void set(const wbvector<double> &R, const wbvector<double> &I);
    wbvector& Conj() { return *this; };

    wbvector& applyFlag(char &flag, const T &fac, wbvector &a);

    bool allEqual() const {
       for (size_t i=1; i<len; ++i) if (data[i]!=data[0]) return 0;
       return 1;
    }

    bool allZero() const { return allEqual(T(0)); };
    bool allOnes() const { return allEqual(T(1)); };

    bool allEqual(const T &a) const {
       bool q=1; size_t i=0;
       if (!len) wblog(FL,"WRN %s() got len=%d",FCT,len);
       for (; i<len; ++i) { if (data[i]!=a) { q=0; break; }}
       return q;
    }

    bool allEqual2(const T &a, T eps=1e-14) const {
       if (!eps) { return allEqual(a); } 
       else {
          T x, x2=0; if (a!=0)  eps*=a;
          if (!len) wblog(FL,"WRN %s() got len=%d",FCT,len);
          for (size_t i=0; i<len; ++i) { x=data[i]-a; x2+=(x*x); }
          return (Wb::abs(x2)<eps);
       }
    };

    bool allUnequal(const T &a) const {
       bool q=1; size_t i=0;
       if (!len) wblog(FL,"WRN %s() got len=%d",FCT,len);
       for (; i<len; ++i) { if (data[i]==a) { q=0; break; }}
       return q;
    }

    bool any() const { bool q=0; size_t i=0;
       for (; i<len; ++i) { if (data[i]) { q=1; break; }}
       return q;
    }

    bool anyUnequal(const T &a) const {
       bool q=0; size_t i=0;
       for (; i<len; ++i) { if (data[i]!=a) { q=1; break; }}
       return q;
    }

    bool anyEqual(const T &a, size_t *i1=nullptr) const {
       bool q=0; size_t i=0;
       for (; i<len; ++i) { if (data[i]==a) {
          if (i1) { (*i1)=i; }
          q=1; break;
       }}
       return q;
    };

    bool anyGT(const T &a) const {
       bool q=0; size_t i=0;
       for (; i<len; ++i) { if (data[i]> a) { q=1; break; }}
       return q;
    };
    bool anyGE(const T &a) const {
       bool q=0; size_t i=0;
       for (; i<len; ++i) { if (data[i]>=a) { q=1; break; }}
       return q;
    };

    bool anyLT(const T &a) const {
       bool q=0; size_t i=0;
       for (; i<len; ++i) { if (data[i]< a) { q=1; break; }}
       return q;
    };
    bool anyLE(const T &a) const {
       bool q=0; size_t i=0;
       for (; i<len; ++i) { if (data[i]<=a) { q=1; break; }}
       return q;
    };

    bool allGT(const T &a) const { return !anyLE(a); } 
    bool allGE(const T &a) const { return !anyLT(a); }
    bool allLE(const T &a) const { return !anyGT(a); }
    bool allLT(const T &a) const { return !anyGE(a); }

    bool allIn(const T &a1, const T &a2) const {
       bool q=1; size_t i=0;
       if (!len) wblog(FL,"WRN %s() got len=%d",FCT,len);
       for (; i<len; ++i) { if (data[i]<a1 || data[i]>a2) { q=0; break; }}
       return q;
    };

    bool anySmallVals(const T eps=1e-15) const {
       size_t i=0; if (eps>0)
            { for (; i<len; ++i) { if (Wb::abs(data[i])<eps) return 1; }}
       else { for (; i<len; ++i) { if (data[i]) return 1; }}
       return 0;
    };

    bool isEqual(const T* v) const {
       return (!memcmp(data, v, len*sizeof(T)));
    }

    bool anyEqual(const wbvector<size_t> &I, const T &x,
       size_t *k=nullptr) const;
    bool anyUnequal(const wbvector<size_t> &I, const T &x,
       size_t *k=nullptr) const;

    bool isEmpty() const {
       if ((!len) ^ (!data)) wblog(FL,
          "ERR severe data inconsistency (%ld,%p)",len,data);
       return (data ? 0 : 1);
    };

    explicit operator bool() const { return !isEmpty(); }
    bool operator! () const { return isEmpty(); }

    bool operator==(const wbvector &v) const { return  isEqual(v); };
    bool operator!=(const wbvector &v) const { return !isEqual(v); };

    bool isEqual(const wbvector &v) const { 
       bool q=1;
       if (this!=&v) {
          if (len!=v.len) { q=0; } else
          if (data!=v.data) {
             for (size_t i=0; i<len; ++i) {
                if (data[i]!=v.data[i]) { q=0; break; }
             }
          }
       }
       return q;
    };

    bool isEqualR(const wbvector &v) const {
       if (len!=v.len) return 0;
       for (size_t l=len-1, i=0; i<len; ++i) {
          if (data[i]!=v.data[l-i]) return 0;
       }
       return 1;
    };

    bool operator^(const wbvector &v) const {
       if (this==&v) { return (len ? 0 : 1); }
       if (len!=v.len) return 0;

       for (size_t i=0; i<len; ++i) {
       if (!((data[i]==0) ^ (v.data[i]==0))) return 0; }

       return 1; 
    };

    bool deepEqualP(const wbvector &B) const;

    char cmp(const wbvector &v) const;

    bool operator< (const wbvector &v) const { return cmp(v)< 0; };
    bool operator<=(const wbvector &v) const { return cmp(v)<=0; };
    bool operator> (const wbvector &v) const { return cmp(v)> 0; };
    bool operator>=(const wbvector &v) const { return cmp(v)>=0; };

    int operator< (const T &r) const {
       for (size_t i=0; i<len; ++i) if (data[i]>=r) return 0;
       return (len ? 1 : 0);
    };

    inline T& at(size_t i) const { 
       if (i>=len) wblog(FL,
          "ERR wbvector::%s() index out of bounds (i=%ld/%ld)",FCT,i,len);
       return data[i];
    };

    inline T& el(long i) const { 
       if (i>=long(len) || -i>long(len)) wblog(FL,
          "ERR wbvector::%s() index out of bounds (i=%ld/%ld)",FCT,i,len);
       return data[i>=0 ? i : long(len)+i];
    };

    T elx(size_t i, const T &x) const {
       return (i<len ? data[i] : x);
    };

    const T& operator[] (size_t i) const { return data[i]; };
          T& operator[] (size_t i)       { return data[i]; };

    T& last() const {
       if (!len) wblog(FL,"ERR %s() got empty %s vector",FCT,TSTR(T));
       return data[len-1];
    };

    T& last(size_t &i) const { 
       if (!len) wblog(FL,"ERR %s() got empty %s vector",FCT,TSTR(T));
       return data[i=len-1];
    };
    template<class TI> 
    T& last(TI &i_) const { size_t i=len-1; i_=i;
       if (!len) wblog(FL,"ERR %s() got empty %s vector",FCT,TSTR(T));
       return data[i];
    };

    T& p2last() const { 
       if (!len) wblog(FL,"ERR %s() got empty %s vector",FCT,TSTR(T));
       return data+(len-1);
    };

    T* ref(size_t i) const { return (data+i); };

    wbvector& setRand(double fac=1., double shift=0.);
    void setIdx (double shift=0., double fac=1.);

    int set2Group(const wbperm &P,
        const WBINDEX &D, const T* S0=nullptr, char iflag=0
    );

    wbvector& Resize_(size_t n) {
       if (isref) wblog(FL,"ERR %s() got reference (%d)",FCT,isref);

       if (n) { 
          if (n<=len && (len<32 || (4*n)>len))
               { len=n; }
          else { Resize(n); }
       }
       else init();

       return *this;
    };

    wbvector& Shorten2(size_t n) {
       if (!n) {
          if (isref) wblog(FL,"ERR %s() got reference (%d)",FCT,isref);
          init();
       }
       else if (n<len && (len<32 || (2*n)>len)) { len=n; }
       else if (n!=len) {
          if (n>len) wblog(FL,
             "WRN %s() got increasing size %ld -> %ld",FCT,len,n);
          Resize(n);
       }
       return *this;
    };

    wbvector& Resize(size_t n, const T* d=nullptr) { 
       if (isref) wblog(FL,"ERR %s() got reference (%d)",FCT,isref);

       if (n==len) { return *this; }
       if (!len) { return RENEW(n,d); }; 
       if (!n) { len=0; if (data) WB_DELETE(data); return *this; }

       T *d0=data;
       WB_NEW(data,n, n>len && !d ? 0 : 1); 
       if (d && n>len)
            { MEM_CPY<T>(data, n, len, d0, d); }
       else { MEM_CPY<T>(data, (n<len? n:len), d0); }

       len=n; WB_DELETE(d0);
       return *this;
    };

    wbvector& resize(size_t n, wbvector& X) const {
       if (!n) { X.init(); } else
       if (n==len) { X=(*this) ; }
       else {
          WB_NEW(X.data,n, n>len ? 0 : 1); X.len=n; 
          if (len) { MEM_CPY<T>(X.data, MIN(n,len), data); }
       }
       return X;
    };

    wbvector& cat(const wbvector &v2) {
       Resize(len+v2.len, v2.data);
       return *this;
    };
    wbvector& Cat(const wbvector &v1, const wbvector &v2) {
       RENEW(v1.len+v2.len);
       MEM_CPY<T>(data,len,v1.len,v1.data,v2.data);
       return *this;
    };
    wbvector& Cat(const wbvector<wbvector  > &vv);
    wbvector& Cat(const wbvector<wbvector* > &vv);

    template <class T1, class T2, class T3, class T4>
    wbvector& Cat(
       const T1* v1, size_t l1, const T2* v2, size_t l2,
       const T3* v3=0, size_t l3=0,
       const T4* v4=0, size_t l4=0
    );

    template <class T1, class T2, class T3, class T4>
    wbvector& Cat(
       const wbvector<T1> *v1,const wbvector<T2> *v2,
       const wbvector<T3> *v3=0,
       const wbvector<T3> *v4=0
    );

    wbvector& Append(const wbvector &v2) {
       Resize(len+v2.len, v2.data);
       return *this;
    };
    wbvector& Append(const T &x) { 
        Resize(len+1, &x);
        return *this;
    };
    wbvector& Append(size_t d2, const T *x) {
        Resize(len+d2, x);
        return *this;
    };

    wbvector append(const wbvector &v) const {
        wbvector v2; append(v.data,v.len,v2); 
        return v2;
    };

    wbvector& append(const wbvector &v2, wbvector &vout) const {
        append(v2.data, v2.len, vout);
        return vout;
    };
    void append(const T *x, size_t d2, wbvector &vout) const {
        vout.RENEW(len+d2);
        MEM_CPY<T>(vout.data,vout.len,len,data,x);
    };

    int recLess (const T* v) const;
    int recLessE(const T* v) const;
    int recCompare(const T* v) const;

    T max() const;
    T max(size_t *k) const; 
    T max(unsigned &k) const { T x; 
      size_t l=0; x=max(&l); Wb::safeConvert(FL,l,k);
      return x;
    };

    T max_(T x) const;

    T min() const;
    T min(size_t *k) const; 
    T min(unsigned &k) const { T x; 
      size_t l=0; x=min(&l); Wb::safeConvert(FL,l,k);
      return x;
    };
    T min_(T x) const;

    T aMax(size_t *k=nullptr) const; 

    T aMax(unsigned &k) const { T x; 
       size_t k_; x=aMax(&k_); Wb::safeConvert(FL,k_,k);
       return x;
    };

    T aMin(char zflag=0, size_t *k=nullptr) const;  

    T maxneg(size_t *k=nullptr) const;
    T minpos(size_t *k=nullptr) const;

    T avg() const; 
    T std() const; 

    void range(T &xmin, T&xmax) const;
    T range() const;

    wbvector& Permute(const wbperm &P); 
    wbvector& permute(wbvector &B, const wbperm &P) const;

    wbvector permute(const wbperm &P) const { 
       wbvector<T> X; permute(X,P);
       return X;
    };

    wbvector& blockPermute(wbvector &B, const wbperm &P) const;
    wbvector& BlockPermute(const wbperm &P) {
       wbvector<T> x(*this);
       return x.blockPermute(*this,P);
    };

    int blockCompare( 
       const wbperm &pa, const wbvector &b, const wbperm &pb) const;

    void get(const WBINDEX &I0, wbvector &v) const {
       const widx_t *I=I0.data;
       size_t i=0, n=I0.len; v.init(n);

       for (; i<n; i++) {
          if (I[i]>=len) wblog(FL,
          "ERR index out of bounds (%d/%d) !?",I[i],len);
          v.data[i]=data[I[i]];
       }
    };

    void get(const wbvector<size_t> &I, T* v) const {
        for (size_t i=0; i<I.len; i++) {
           if (I[i]>=len)
           wblog(FL,"ERR index out of bounds (%d/%d) !?", I[i],len);

           v[i]=data[I[i]];
        }
        return;
    };

    void get(size_t i1, size_t i2, wbvector &v) const {

        if (i1>=len || i2>=len) wbdie(FL,"index out of range.");
        if (i2<i1)
             v.init();
        else v.init(i2-i1, data+i1);
    };

    template<class TI>
    void selectSU(const wbvector<TI> &I);

    template<class TI>
    void select(const wbvector<TI> &I, T* r) const;

    template<class TI>
    wbvector& Select(const wbvector<TI> &I);

    template<class TI>
    wbvector& select(const wbvector<TI> &I, wbvector &B) const;

    template<class TI>
    wbvector  select(const wbvector<TI> &I) const {
       wbvector B; return select(I,B); 
    };

    wbvector& Select(cMVEC &mark, char nonzero=1);
    wbvector& select(cMVEC &mark, wbvector &B, char nonzero=1) const;
    wbvector  select(cMVEC &mark, char nonzero=1) const {
       wbvector B; return select(mark,B,nonzero); 
    };

    template<class TI>
    wbvector& BlockSelect(const wbvector<TI> &I, const wbvector<TI> &D);

    wbvector& getI(size_t k, wbvector &v) const;
    wbvector& getI(const wbindex &k, wbvector &v, char uflag=1) const;

    wbvector& getI(const WBINDEX &k, wbvector &v, char uflag=1) const {
       return getI((const wbindex&)k,v,uflag);
    };

    wbvector& initI(const wbvector &a, const wbindex &ia){
       return a.getI(ia,*this); };
    wbvector& initI(
       const wbvector &a, const wbindex &ia,
       const wbvector &b, const wbindex &ib,
       char uflag=1 
    );

    wbvector& Flip() { 
        size_t i, m2=len/2, l=len-1;
        for (i=0; i<m2; i++) SWAP(data[i],data[l-i]);
        return *this;
    };

    void flip(wbvector &a) const { 
        if (&a==this) { a.Flip(); return; }
        a.init(len);
        for (size_t l=len-1, i=0; i<len; ++i) a.data[i]=data[l-i];
    };

    wbvector flip() const {
       wbvector a; flip(a); 
       return a;
    };

    wbvector operator[] (const wbvector<size_t>& I) const {
       wbvector v(I.len); 
       size_t i=0, *const k=I.data;

       for (; i<I.len; ++i) {
          if (k[i]>=len) wblog(FL,
             "ERR index out of bounds (%d: %d/%d)",i,I.max(),len);
          v[i]=data[k[i]];
       }

       return v;
    };

    T prod() const {
       if (!len) wblog(FL,"WRN %s() empty data set (returning 0)",FCT);
       return prod(0);
    };

    T prod(T x) const {
       if (len) {
          size_t i=1; x=data[0]; for (; i<len; ++i) { x*=data[i]; }
       }; return x;
    };

    T prod(size_t i, size_t i2) const {
       T x=0; 
       if (i<len && i<=i2) {
          if ((++i2)>len) { i2=len; }
          for (x=data[i++]; i<i2; ++i) { x*=data[i]; }
       }
       else wblog(FL,
         "ERR empty index set [%d,%d] / %d (returning 0)",i,i2,len);
       return x;
    };

    T prod(size_t i, size_t i2, T x) const { 
       if (long(i )<0) { i =len+i; } 
       if (long(i2)<0) { i2=len+i2; }
       if (i<len && i<=i2) {
          if ((++i2)>len) { i2=len; }
          for (x=data[i++]; i<i2; ++i) { x*=data[i]; }
       }
       return x;
    };

    template <class TI>
    T prod(const TI *I, size_t n) const {
       T x=0; if (n) { size_t i=1;
         x=at(I[0]); for (; i<n; ++i) { x*=at(I[i]); }
       }
       return x;
    };

    template <class T1, class T2>
    T2 prod(const T1 *I, size_t n, T2 &s2) const {
       size_t i,k;
       std::vector<char> m(len); 

       for (i=0; i<n; ++i) { k=I[i];
          if (k>=len) wblog(FL,
             "ERR %s() index out of bounds (%d/%d)",FCT,k,len);
          if ((++m[k])>1) wblog(FL,"ERR %s() index not unqiue",FCT);
       }
       if (!len) { s2=0; return 0; }
       else {
          T2 s1=1; s2=1;
          for (i=0; i<len; ++i) { (m[i] ? s1 : s2) *= data[i]; }
          return s1;
       }
    };

    void prod2(T &p1, T &p2, const T &p0=0) const {
        size_t l=len/2; if (len%2) wblog(FL,
        "ERR cannot split vector into two equal parts (%d)",len);

        if (len==0) { p1=p0; p2=p0; return; }
        p1=data[0]; p2=data[l];

        for (size_t i=1; i<l; i++) { p1*=data[i]; p2*=data[i+l]; }
    }

    void info(const char *istr) const;
    void print(const char *istr, char mflag=0) const;
    void printdata(const char *istr=0) const {
       if (istr && istr[0]) printf("%s = ",istr);
       printf("[ %s ]\n",STR(*this));
    };

    void put(
      const char *F, int L,
      const char *vname="ans", char tflag=0, const char *ws="base"
    ) const;

    void put(const char *vname="ans", char tflag=0, const char *ws="base"
    ) const { put(0,0,vname,tflag,ws); }

    void dispWithMsg(const char* file, int line,
    const char *msg, const char* vname="") const {
        print(vname);
        wblog(file,line,msg);
    }

    void mat2mx(mxArray* &a, char dim=2) const;
    void mat2mxs(mxArray* a, char field_nr, char dim=2) const;

    mxArray* toMx_Struct (const char tflag=0) const;
    mxArray* toMxP_S(const char tflag=0) const; 

    mxArray* toMx (const char tflag=0) const { return toMx_Struct (tflag); }; 
    mxArray* toMxP(const char tflag=0) const { return toMxP_S(tflag); };

    mxArray* toMx_base(const char tflag=0) const {
       size_t m=1, n=len; if (tflag) { m=n; n=1; } 
       return Mx::Array<double>(m,n).copyFromTR(data);
    };

    mxArray* toMx_base_d(const char tflag=0) const {
       size_t m=1, n=len; if (tflag) { m=n; n=1; } 
       return Mx::Array<double>(m,n).copyFromTR(data);
    };

    mxArray* toMx_offset(const T& offset, const char tflag=0) const;

    mxArray* mxCreateStruct(unsigned m, unsigned n) const;
    void add2MxStruct(mxArray *S, unsigned i, char tflag=0) const;
    void add2MxStruct(mxArray *S, const char *vname, char tflag=0) const;

    mxArray* mxCreateCell(unsigned m, unsigned n) const;
    void add2MxCell(mxArray *S, unsigned i, char tflag=0) const;

    void putg(const char *vname="ans", char tflag=0) const {
        put(vname,tflag,"global");
        return;
    };

    void putx(const char *vname="ans", char tflag=0) const {
        put(vname,tflag,"base");
        return;
    };

    wbstring toStr(int n, const char *sep) const;
    wbstring toStr() const;  

    wbstring sizeStr() const;

    wbstring toStrf(
       const char *fmt0="", 
       const char *sep=" ",
       unsigned stride=0,   
       const char *sep2=""  
    ) const;

    T *data;

    size_t len;

    char isref;

  protected: 

    wbvector& RENEW(size_t n, const T* d=nullptr, size_t l=-1, char init=1) {
       if (data) {
          if (isref ) { data=nullptr; } else
          if (n!=len) { WB_DELETE(data); } 
       }
       isref=0; if (!n) { len=0; return *this; }

       if (!d ) { l=0; } else
       if (l>n) { l=n; } 

       if (init && l>=n) { init=0; }

       if (!data) { WB_NEW(data,n, init? 0 : 1); } else
       if (init ) { MEM_CPY<T>(data+l,n-l,nullptr); } 

       if (l) { MEM_CPY<T>(data,l,d); }
       len=n;

       return *this;
    };

    wbvector& RENEW_VAL(size_t n, const T& x) {
        if (isref) { data=nullptr; isref=0; } else
        if (n!=len && data) { WB_DELETE(data); } 

        if ((len=n)) {
           if (!data) { WB_NEW(data,n,1); } 
           for (size_t i=0; i<n; ++i) { data[i]=x; }
        }
        return *this;
    };

  private:

}; 

template<class T>
inline wbvector<const T*>cPVEC1_(const T &A) {
   return wbvector<const T*>(1,&A);
};

template<class T> 
inline wbvector<const T*> cPVEC__(const wbvector<T> &A) {
   wbvector<const T*> a(A.len);
   for (unsigned i=0; i<A.len; ++i) { a[i]=&A[i]; }
   return a;
};

template<class T> 
inline unsigned cPVEC_nnz(
   const  char *F, int L, const char *fct, wbvector<const T*> &A,
   char lflag=0) { 

   unsigned i=0, l=-1, nwrn=0;
   for (; i<A.len; ++i) {
       if (A[i]) {
          if (*A[i]) { if (++l<i) { A[l]=A[i]; }}
          else {
             if (lflag=='w') { ++nwrn; wblog(F_L,
             "WRN %s() got empty data[%d/%d]",fct?fct:FCT,i+1,A.len); }
          }
       }
       else {
          if (tolower(lflag)=='w') { ++nwrn; wblog(F_L, 
          "WRN %s() got null data[%d/%d]",fct?fct:FCT,i+1,A.len); }
       }
   }

   if (++l) { A.len=l; }
   else {
      if ((F || fct) && !nwrn) { char sout[128]; snprintf(sout,128,
         "%s() got empty data (l=%d/%ld)",fct?fct:FCT,l,A.len);
         if (lflag) 
              { wblog(F_L,"WRN %s",sout); }
         else { wblog(F_L,"ERR %s",sout); }
      }
      A.init();
   }
   return A.len;
};

template <class T>
wbvector<T>& wbvector<T>::unRef() {

    if (!isref) return *this;
    if (!len) {
       if (data) wblog(FL,"ERR data=%lX (%d)",data,len);
       isref=0; return *this;
    }

    T const* const d0=data;

    WB_NEW(data,len,1); 
    MEM_CPY<T>(data,len,d0);
    isref=0;

    return *this;
};

template <class T>
bool wbvector<T>::deepEqualP(const wbvector<T> &B) const {
    if (len!=B.len) return 0;

    if (data!=B.data) {
       for (size_t i=0; i<len; i++) { if (data[i]==B.data[i]) continue;
          if ((!data[i] || !B.data[i])) return 0; else
          if (!((*data[i])==(*B.data[i]))) return 0;
       }
    }

    return 1;
};

template <class T> inline 
char wbvector<T>::cmp(const wbvector<T> &v) const {

   size_t i=0, n=MIN(len,v.len); const T *b=v.data;
   for (; i<n; ++i) {
      if (data[i]<b[i]) return -1;
      if (data[i]>b[i]) return +1; 
   }

   if (len<v.len) return -1;
   if (len>v.len) return +1; else return 0;
};

template <class T> inline
bool wbvector<T>::checkSameLength(
   const wbvector<T> &a,
   const char *F, int L, const char *istr
) const {

   bool i=(len==a.len);
   if (!i && F) { 
      if (istr && istr[0]) sprintf_str("%s() - ",istr); else str[0]=0;
      wblog(F,L, "ERR %svector length mismatch (%d/%d)",str,len,a.len);
   }
   return i;
};

template <class T> inline
bool wbvector<T>::checkSameLength(
   const wbvector<T> &a,
   const wbvector<T> &b,
   const char *F, int L, const char *istr
) const {

   bool i=(len==a.len && len==b.len);
   if (!i && F) { 
      if (istr && istr[0]) sprintf_str("%s() - ",istr); else str[0]=0;
      wblog(F,L,"ERR %svector length mismatch (%d,%d/%d)",str,a.len,b.len,len);
   }
   return i;
};

template <class T> inline
void wbvector<T>::getReal(wbvector<double> &R) const {
   R.init(len);
   for (size_t i=0; i<len; ++i) { R.data[i]=double(data[i]); }
};

template <> inline
void wbvector<wbcomplex>::getReal(wbvector<double> &R) const {
   R.init(len);
   for (size_t i=0; i<len; ++i) { R.data[i]=data[i].r; }
};

template <class T> inline
void wbvector<T>::getImag(wbvector<double> &I) const {
   I.init(len);
};

template <> inline
void wbvector<wbcomplex>::getImag(wbvector<double> &I) const {
   I.init(len);
   for (size_t i=0; i<len; ++i) { I.data[i]=data[i].i; }
};

template <> inline
void wbvector<wbcomplex>::set( 
   const wbvector<double> &R, const wbvector<double> &I) {

   size_t i=0; const double *re=R.data, *im=I.data;

   init(R.len ? R.len : I.len);

   if (R.len) {
      if (I.len) {
         if (R.len!=I.len) wblog(FL,
            "ERR %s() got length inconsistency (%d/%d)",FCT,R.len,I.len);
         for (; i<len; ++i) { data[i].init(re[i],im[i]); }
      }
      else { for (; i<len; ++i) { data[i].init(re[i],0); }}
   }
   else if (I.len) { for (; i<len; ++i) { data[i].init(0,im[i]); }}
};

template <> inline
wbvector<wbcomplex>& wbvector<wbcomplex>::Conj() {
   for (size_t i=0; i<len; i++) data[i].Conj();
   return *this;
};

template <class T> inline
wbvector<T>& wbvector<T>::applyFlag(
   char &flag, const T &fac, wbvector<T> &a
){
   if (!strchr("NCTc",flag))
   wblog(FL,"ERR %s() invalid flag %c<%d>",__FUNCTION__,flag,flag);

   if (flag=='c') flag='C';
   if (flag=='C' && typeid(T)!=typeid(wbcomplex)) flag='N';

   if (fac!=1 || flag=='C') { a=(*this);
      if (flag=='C') a.Conj();
      if (fac!=1) a*=fac;
   }
   else a<<*this;
};

template <class T> inline
T wbvector<T>::min() const {

   if (len) { size_t i=1; T x=data[0];
      for (; i<len; ++i) { if (x>data[i]) { x=data[i]; }}
      return x;
   }
   wblog(FL,"ERR %s() got empty vector",FCT);
   return 0;
};

template <class T> inline
T wbvector<T>::min_(T x) const {
   if (len) { size_t i=1; x=data[0];
      for (; i<len; ++i) { if (x>data[i]) { x=data[i]; }}
   }
   return x;
};

template <class T> inline
T wbvector<T>::min(size_t *k) const {
   if (!k) { return min(); }
   if (len) { size_t i=1; T x=data[0]; *k=0;
      for (; i<len; ++i) { if (x>data[i]) { x=data[i]; *k=i; }}
      return x;
   }
   wblog(FL,"ERR %s() got empty vector",FCT);
   return 0;
};

template <class T> inline
T wbvector<T>::max() const {

   if (len) { unsigned i=1; T x=data[0];
      for (; i<len; ++i) { if (x<data[i]) { x=data[i]; }}
      return x;
   }
   wblog(FL,"ERR %s() got empty vector",FCT);
   return 0;
};

template <class T> inline
T wbvector<T>::max_(T x) const {
   if (len) { size_t i=1; x=data[0];
      for (; i<len; ++i) { if (x<data[i]) { x=data[i]; }}
   }
   return x;
};

template <class T> inline
T wbvector<T>::max(size_t *k) const {
   if (!k) { return max(); }
   if (len) { size_t i=1; T x=data[0]; *k=0;
      for (; i<len; ++i) { if (x<data[i]) { x=data[i]; *k=i; }}
      return x;
   }
   wblog(FL,"ERR %s() got empty vector",FCT);
   return 0;
};

template <class T> inline
T wbvector<T>::minpos(size_t *k) const {

   if (len) {
      T x=0; size_t i=0; if (k) { (*k)=-1; }
      for (; i<len; ++i) { if (data[i]>0) {
          x=data[i]; if (k) { (*k)=i; }; ++i; break;
      }}
      for (; i<len; ++i) { if (data[i]>0 && x>data[i]) {
          x=data[i]; if (k) { (*k)=i; }
      }}
      return x;
   }
   wblog(FL,"WRN %s() of empty vector",FCT);
   return 0;
};

template <class T> inline
T wbvector<T>::maxneg(size_t *k) const {

   if (len) {
      T x=0; size_t i=0; if (k) { (*k)=-1; }
      for (; i<len; ++i) { if (data[i]<0) {
         x=data[i]; if (k) { (*k)=i; }; ++i; break;
      }}
      for (; i<len; ++i) { if (data[i]<0 && x<data[i]) {
         x=data[i]; if (k) { (*k)=i; }
      }}
      return x;
   }
   wblog(FL,"WRN %s() of empty vector",FCT);
   return 0;
};

template <class T> inline
T wbvector<T>::aMax(size_t *k_) const {

   if (len) {
      size_t k=-1, i=0; T a,x=0;
      for (; i<len; ++i) {
         a=Wb::abs(data[i]); if (x<a) { x=a; k=i; }
      }
      if (k_) (*k_)=k;
      return x;
   }
   else {
      if (k_) { (*k_)=-1; wblog(FL,"WRN %s() got null vector",FCT); }
      return 0;
   }
};

template <class T> inline
T wbvector<T>::aMin(char zflag, size_t *k_) const {

   if (!len) {
      if (k_) { (*k_)=-1; wblog(FL,"WRN %s() got null vector",FCT); }
      return 0;
   }

   size_t i=0, k=-1; double a,x=0;

   if (zflag==0) { 
      x=Wb::abs(data[i++]); k=0; if (x!=0) {
         for (; i<len; ++i) { a=Wb::abs(data[i]); if (x>a) {
            x=a; k=i; if (x==0) break;
         }}
      }
   }
   else {
      for (; i<len; ++i) {
         a=Wb::abs(data[i]); if (a) { x=a; k=i; break; }
      }
      for (; i<len; ++i) { a=Wb::abs(data[i]);
         if (a && x>a) { x=a; k=i; }
      }
   }
   if (k_) (*k_)=k;
   return x;
};

template <class T> inline
void wbvector<T>::range(T &xmin, T&xmax) const {
   if (len==0) {
       xmin=xmax=0;
       wblog(FL,"WRN %s() of null vector",FCT);
   }
   else {
      xmin=xmax=data[0];
      for (size_t i=1; i<len; i++) {
         if (xmin>data[i]) xmin=data[i];
         if (xmax<data[i]) xmax=data[i];
      }
   }
};

template <class T> inline
T wbvector<T>::range() const {
   T xmin, xmax; range(xmin,xmax);
   return (xmax-xmin);
};

template <class T> inline
T wbvector<T>::avg() const {

   if (!len) {
      wblog(FL,"WRN %s() of null vector !?",FCT);
      return 0;
   }
   else {
      T x=data[0];
      for (size_t i=1; i<len; ++i) { x+=data[i]; }
      return (x/T(len));
   }
};

template <class T> inline
T wbvector<T>::std() const {

   if (len<2) {
      wblog(FL,"WRN %s() of vector of length %d !?",FCT,len);
      return 0;
   }
   else {
      T x=data[0], x2=data[0]*data[0];
      for (size_t i=1; i<len; ++i) { x+=data[i]; x2+=data[i]*data[i]; }
      return Wb::sqrt( ( x2 - (x*x)/T(len) ) / T(len-1) );
   }
};

template <class T> inline
int wbvector<T>::findClosest(T r) const {
   if (len==0) return -1;

   T dbl, eps=Wb::abs(data[0]-r); size_t i,k=0;

   for (i=1; i<len; i++) {
      dbl = Wb::abs(data[i]-r);
      if (dbl<eps) { k=i; eps=dbl; if (eps==0) break; }
   }

   return k;
}

template <class T>
int wbvector<T>::findClosestSorted(T r) const {
   if (len==0) return -1;

   size_t i, k, n1=0, n2=len, e=0;
   char c, cref;

   if (len==0) return -1;
   if (len==1) return 0;

   cref=NUMCMP(data[0], data[len-1]);

   if (cref) {
      if (len>6) {
        for (i=1; i<4; i++) if (NUMCMP(data[i-1],data[i])==-cref) e++;
        for (i=len-3; i<len; i++)
        if (NUMCMP(data[i-1],data[i])==-cref) e++;
      }
      else
      for (i=1; i<len; i++) if (NUMCMP(data[i-1],data[i])==-cref) e++;
   }
   else for (i=1; i<len; i++) if (NUMCMP(data[i-1],data[i])) e++;

   if (e) wblog(FL,"ERR Data not sorted! (%d)", e);
   if (cref==0) { wblog(FL,
      "WRN find closest: all data constant!?\n(%g; %g %dx)",
      (double)r, (double)data[0], len); return 0;
   }

   c=NUMCMP(data[len-1],r); if (c!=-cref) return len-1;
   c=NUMCMP(data[0],r);     if (c!= cref) return 0;

   for (k=0;;) {
      if (c==cref) {
          n1=n2/2; n2-=n1; if (n1) k+=n1; else {
             if (k+1<len && (Wb::abs(r-data[k]) > Wb::abs(r-data[k+1]))) k++;
             break;
          }
      }
      else if (c) {
          n2=n1/2; n1-=n2; if (n2) k-=n2; else {
             if (k) if (Wb::abs(r-data[k]) > Wb::abs(r-data[k-1])) k--;
             break;
          }
      }
      else break;

      c=NUMCMP(data[k],r);
   }

   return k;
}

template <class T>
wbvector<T>& wbvector<T>::setRand(double fac, double shift) {
   size_t i;
   static char first_call=1;

   if (first_call) { wb_srand(); first_call=0; } 
   fac/=(double)RAND_MAX;

   if (shift==0.)
        for (i=0; i<len; ++i) data[i]=(T)(fac*rand());
   else for (i=0; i<len; ++i) data[i]=(T)(fac*rand()+shift);

   return *this;
};

template<>
wbvector<wbcomplex>& wbvector<wbcomplex>::setRand(double fac, double shift) {
   size_t i;
   static char first_call=1;

   if (first_call) { wb_srand(); first_call=0; } 
   fac/=(double)RAND_MAX;

   if (shift==0.) {
      for (i=0; i<len; ++i)
      data[i].set(fac*rand(), fac*rand());
   }
   else {
      for (i=0; i<len; ++i)
      data[i].set(fac*rand()+shift, fac*rand()+shift);
   }

   return *this;
};

template <class T> inline
void wbvector<T>::setIdx(double shift, double fac) {
   size_t i;

   if (shift==0.) {
        if (fac==1.)
             for (i=0; i<len; i++) data[i]=(T)(i);
        else for (i=0; i<len; i++) data[i]=(T)(fac*i);
   }
   else {
        if (fac==1.)
             for (i=0; i<len; i++) data[i]=(T)(i+shift);
        else for (i=0; i<len; i++) data[i]=(T)(fac*i+shift);
   }
};

template <class T>
void wbvector<T>::mat2mx(mxArray* &a, char dim) const {

    if (dim==1) a=mxCreateDoubleMatrix(len,1,mxREAL); else
    if (dim==2) a=mxCreateDoubleMatrix(1,len,mxREAL);
    else wblog(FL,"ERR %s() invalid dim=%d",FCT,dim);

    double *d=mxGetDoubles(a);
    for (size_t i=0; i<len; ++i) d[i]=(double)data[i];
};

template <class T>
void wbvector<T>::mat2mxs(mxArray* a, char fid, char dim) const {

    size_t i=-1, m=1, n=len;

    if (!mxIsStruct(a) || (i=mxGetNumberOfElements(a))!=1) wblog(FL,
       "ERR %s() invalid input structure (n=%d)",FCT,i);
    if (fid>=int(i=mxGetNumberOfFields(a))) wblog(FL,
       "ERR %s() field index out of bounds (%d/%d)",FCT,fid,i);

    if (dim==1) { m=len; n=1; } else
    if (dim==2) { m=1; n=len; } else
       wblog(FL,"ERR %s() invalid dim=%d",FCT,dim);

    mxArray *b=mxCreateDoubleMatrix(m,n,mxREAL);
    double  *d=mxGetDoubles(b);

    for (i=0; i<len; ++i) { d[i]=(double)data[i]; }

    mxSetFieldByNumber(a,0,fid,b);
};

template <class T> inline
int wbvector<T>::recLess(const T* v) const {

    for (size_t i=0; i<len; i++)
    if (data[i]<v[i]) return 1; else
    if (data[i]>v[i]) return 0;

    return 0;
}

template <class T> inline
int wbvector<T>::recLessE(const T* v) const {

    for (size_t i=0; i<len; i++)
    if (data[i]<v[i]) return 1; else
    if (data[i]>v[i]) return 0;

    return 1;
}

template <class T> inline
int wbvector<T>::recCompare(const T* v) const {

    for (size_t i=0; i<len; i++)
    if (data[i]<v[i]) return -1; else
    if (data[i]>v[i]) return +1;

    return 0;
}

template<> 
mxArray* wbvector<long>::toMx_base(const char tflag) const {
    size_t m=1, n=len; if (tflag) { m=n; n=1; } 
    return Mx::Array<long>(m,n).copyFromTR(data);
};

template<> 
mxArray* wbvector<size_t>::toMx_base(const char tflag) const {
    size_t m=1, n=len; if (tflag) { m=n; n=1; } 
    return Mx::Array<size_t>(m,n).copyFromTR(data);
};

template<> 
mxArray* wbvector<wbcomplex>::toMx_base(const char tflag) const {
    size_t m=1, n=len; if (tflag) { m=n; n=1; } 
    return Mx::Array<wbcomplex>(m,n).copyFromTR(data);
};

template<> int wbvector<double>::init(
  const char *F, int L, const mxArray *a, const char *istr, char tcheck, char ref)
{ return init_base(F,L,a,istr,tcheck,ref); }

template<> int wbvector<wbcomplex>::init(
  const char *F, int L, const mxArray *a, const char *istr, char tcheck, char ref)
{ return init_base(F,L,a,istr,tcheck,ref); }

template<> int wbvector<char>::init(
  const char *F, int L, const mxArray *a, const char *istr, char tcheck, char ref)
{ return init_base(F,L,a,istr,tcheck,ref); }

template<> int wbvector<int>::init(
  const char *F, int L, const mxArray *a, const char *istr, char tcheck, char ref)
{ return init_base(F,L,a,istr,tcheck,ref); }

template<> int wbvector<unsigned>::init(
  const char *F, int L, const mxArray *a, const char *istr, char tcheck, char ref)
{ return init_base(F,L,a,istr,tcheck,ref); }

template<> int wbvector<size_t>::init(
  const char *F, int L, const mxArray *a, const char *istr, char tcheck, char ref)
{ return init_base(F,L,a,istr,tcheck,ref); }

template<> mxArray* wbvector<     double>::toMx(const char tflag) const {
  return toMx_base(tflag); }
template<> mxArray* wbvector<  wbcomplex>::toMx(const char tflag) const {
  return toMx_base(tflag); }

template<> mxArray* wbvector<     size_t>::toMx(const char tflag) const {
  return toMx_base(tflag); };
template<> mxArray* wbvector<       long>::toMx(const char tflag) const {
  return toMx_base(tflag); };
template<> mxArray* wbvector<   unsigned>::toMx(const char tflag) const {
  return toMx_base(tflag); };
template<> mxArray* wbvector<        int>::toMx(const char tflag) const {
  return toMx_base(tflag); }
template<> mxArray* wbvector<      short>::toMx(const char tflag) const {
  return toMx_base(tflag); }; 

template<> mxArray* wbvector<       char>::toMx(const char tflag) const {
  return toMx_base(tflag); }
template<> mxArray* wbvector<     int8_t>::toMx(const char tflag) const {
  return toMx_base(tflag); }
template<> mxArray* wbvector<    uint8_t>::toMx(const char tflag) const {
  return toMx_base(tflag); } 

template<> mxArray* wbvector<   uint16_t>::toMx(const char tflag) const {
  return toMx_base(tflag); };

template<class T>
mxArray* wbvector<T>::toMx_Struct(const char tflag) const { 

    mxArray *S;

    if (tflag)
         S=data->mxCreateStruct(len,1);
    else S=data->mxCreateStruct(1,len);

    for (size_t i=0; i<len; ++i) { data[i].add2MxStruct(S,i); }
    return S;
}

template<class T>
mxArray* wbvector<T>::toMxP_S(const char tflag) const {

    mxArray *S;

    if (tflag)
         S=(*data)->mxCreateStruct(len,1);
    else S=(*data)->mxCreateStruct(1,len);

    for (size_t i=0; i<len; ++i) { data[i]->add2MxStruct(S,i); }
    return S;
};

template<class T>
mxArray* wbvector<T>::toMx_offset(const T& offset, const char tflag)
const {
    if (!len) { return mxCreateDoubleMatrix(0,0,mxREAL); }

    mxArray *a=(tflag ?
       mxCreateDoubleMatrix(len,1,mxREAL) :
       mxCreateDoubleMatrix(1,len,mxREAL)
    );

    double *d=mxGetDoubles(a);
    for (size_t i=0; i<len; ++i) { d[i]=double(data[i]+offset); }

    return a;
};

template <class T>
void wbvector<T>::put(
    const char *file, int line, const char *vname, char tflag,
    const char* ws
) const {

    size_t i;
    mxArray *a=toMx(tflag);

    if (file) wblog(file,line,
    "I/O putting vector '%s' to %s%s",vname,ws,tflag ? " (t)" : ".");

    i=mexPutVariable(ws, vname, a);

    if (i==1) {
       info("this"); wblog(__FL__,
      "ERR failed to put variable %s into workspace %s.",vname,ws);
    }

    mxDestroyArray(a);
    return;
}

inline WBINDEX Index(size_t i1, size_t i2) {
    WBINDEX ii;

    if (i1>i2) return ii;

    ii.init(i2-i1+1);
    for (size_t i=0; i<ii.len; i++) ii[i]=i+i1;

    return ii;
};

inline char isUniqueIdxSet(const WBINDEX &I_, widx_t imax) {

   size_t i, n=I_.len;
   const widx_t *I=I_.data;
   wbvector<char> mark(imax); char *m=mark.data;

   for (i=0; i<n; ++i) {
      if (I[i]>=imax || ++m[I[i]]>1) { break; }
   }

   return (i==n);
};

inline void getIPerm(const WBPERM &P, WBPERM &iP) {

   size_t i,n=P.len;

 #ifdef __WBDEBUG__
   if (P.isValidPerm()<=0) { wbdie(FLINE,"invalid permutation"); }
  #endif

    iP.init(n);
    for (i=0; i<n; i++) iP[P[i]]=i;

    return;
}

char isIdentityPerm(const WBPERM &P) {

    for (size_t i=0; i<P.len; i++) if (P[i]!=i) return 0;
    return 1;
};

#endif

