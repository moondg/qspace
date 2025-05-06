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

#ifndef __WB_ARRAY_COL_MAJOR_CC__
#define __WB_ARRAY_COL_MAJOR_CC__

//-------------------------------------------------------------------//
// wbarray - N-D Array class
//
//    generic N-dimensional class that allows permutations
//    and reshaping; index order is column-major.
//
// Wb,Sep14,09 ;  Wb,Aug11,05
//-------------------------------------------------------------------//

template<class T>
wbindex& wbarray<T>::ind2sub(size_t k, wbindex &I) const { 
   const size_t *s=SIZE.data, n=SIZE.len;
   size_t i,r, e=0, *idx;

   I.init(n); idx=I.data; 
   for (i=0; i<n; ++i) {
       if (!s[i]) wblog(FL,
          "ERR wbarray::%s() got size %s",FCT,SSTR_(this));

       r=k/s[i]; idx[i]=k-r*s[i]; k=r;
       if (idx[i]>=s[i]) e++;
   }

   if (e || k) wblog(FL,
      "ERR wbarray::%s() out of bounds (%d => [%s; %s])",
       FCT, STR(I+1), SSTR_(this)
   );

   return I;
};

template<class T>
void wbarray<T>::setRand(const char pnflag) {
    static int firstcall=1; size_t i,s=numel(); double fac;

    if (WbUtil<T>::isFloat())
         fac = 1.  /(double)RAND_MAX;
    else fac = 100./(double)RAND_MAX;

    if (firstcall) { wb_srand(); firstcall=0; } 
    for (i=0; i<s; ++i) data[i] = (T)(fac*::rand()); 

    if (pnflag) {
       fac *= (0.5*(double)RAND_MAX);
       for (i=0; i<s; ++i) data[i] -= fac;
    }
}

template<class T>
double wbarray<T>::SkipTiny(double eps_) {

   T zero=T(0), x2=zero, eps=T(eps_);
   if (eps!=zero) { size_t i=0, n=numel(); T a;
      for (; i<n; ++i) { if (data[i]!=zero) {
         a=Wb::abs(data[i]); if (a<eps) { x2+=Wb::norm2(a); data[i]=zero; }
      }}
   }
   return double(x2); 
};

template<> inline 
double wbarray<wbcomplex>::SkipTiny(double eps) {
   double a,x2=0;
   for (size_t s=numel(), i=0; i<s; ++i) {
      a=fabs(data[i].r); if (a>0 && a<eps) { x2+=a*a; data[i].r=0; }
      a=fabs(data[i].i); if (a>0 && a<eps) { x2+=a*a; data[i].i=0; }
   }
   return std::sqrt(x2);
}

template<> inline 
wbarray<double>& wbarray<double>::SkipTiny_float(double ref) {
   Wb::chopTiny_float(data,numel(),ref);
   return *this;
};

template<> inline 
wbarray<wbcomplex>& wbarray<wbcomplex>::SkipTiny_float(wbcomplex ref) {
   Wb::chopTiny_float((double*)data,2*numel(),ref.r);
   return *this;
};

template<> inline 
double wbarray<wbcomplex>::SkipTiny_imag(double eps) {
   return Wb::chopTiny_imag(data,numel(),eps);
};

template<class T>
wbvector<T>& wbarray<T>::normDim2(unsigned id, wbvector<T> &xk) const {

   size_t i, n, N=numel(); if (!N) { return xk.init(); }
   unsigned k, r=SIZE.len, l=r-1;
   const size_t *s=SIZE.data; T *x2;
   WBINDEX I(r);

   if (id>=r) { wblog(FL,"ERR %s() index out of bounds (%d/%d)",FCT,id,r); }
   n=SIZE[id]; xk.init(n); x2=xk.data;

   for (i=0; i<N; ++i) {
      x2[I[id]] += Wb::norm2(data[i]);

      k=0; ++I[0];  
      while (I[k]>=s[k] && k<l) { I[k]=0; ++I[++k]; }
   }

   return xk;
};

template<class T>
wbarray<T>& wbarray<T>::SkipTrailingZeroSpace(
   unsigned k,  
   size_t n,    
   double eps
) {
   wbvector<T> x2; normDim2(k,x2); if (!x2) { return *this; }

   size_t N=SIZE[k]; 

   if (long(n)<0) { n=0; }
   if (n<N) {
      size_t i=N-1;
      double eps2=x2.max(); if (eps2<1) { eps2=1; }
      eps2*=(eps*eps);

      for (; i<N; --i) { if (x2[i]>eps2) { break; }}
      if ((++i)<N) { if (i<n) { i=n; }
         wbvector<size_t> S(SIZE); S[k]=i;
         Resize(S);
      }
   }
   return *this;
};

template<class T> inline
T wbarray<T>::normDiff2(const wbarray<T> &B, char sflag) const {

   T x2=0; 
   if (sameSize(B)) {
      x2=Wb::rangeNormDiff2(data, B.data, numel());
      return x2;
   }

   if (!sflag || SIZE.len!=2 || B.SIZE.len!=2) {
      if (sflag) wblog(FL,"ERR %s() size mismatch\n"
         "(sflag=%d for matrices only: %s / %s",FCT,sflag,SSTR_(this),SSTR(B));
      else wblog(FL,"ERR size mismatch (%s; %s)",SSTR_(this),SSTR(B));
   }

   char b1=1, b2=1; 
   size_t i,j, d1=SIZE[0], d2=SIZE[1], D1=B.SIZE[0], D2=B.SIZE[1];
   const T* a=data, *b=B.data;

   if (d1>D1) { SWAP(d1,D1); b1=0; } 
   if (d2>D2) { SWAP(d2,D2); b2=0; } 

   for (j=0; j<D2; ++j) { 
      if (j<d2) {
         for (i=0; i<d1; ++i) { x2+=Wb::norm2(a[i]-b[i]); }
         if (b1)
              { for (; i<D1; ++i) { x2+=Wb::norm2(b[i]); }}
         else { for (; i<D1; ++i) { x2+=Wb::norm2(a[i]); }}
         a+=SIZE[0]; b+=B.SIZE[0]; 
      }
      else if (b2) {
         if (j==d2) { D1=B.SIZE[0]; }
         for (i=0; i<D1; ++i) { x2+=Wb::norm2(b[i]); }; b+=D1;
      }
      else {
         if (j==d2) { D1=SIZE[0]; }
         for (i=0; i<D1; ++i) { x2+=Wb::norm2(a[i]); }; a+=D1;
      }
   }
   return x2;
};

template<class T> inline
T wbarray<T>::normDiff2(const wbvector<T> &b) const {
   if (!isVector() || numel()!=b.len) wblog(FL,
      "ERR size mismatch (%s; %d)",SSTR_(this),b.len);
   return Wb::rangeNormDiff2(data, b.data, numel());
};

template<class T>
wbarray<T>& wbarray<T>::plus( 
   const wbarray &B, wbarray &C, T bfac, char iflag
 ) const {

   if (isEmpty() && B.isEmpty()) { return C.init(); }

   if (C.isRef()) wblog(FL,"ERR %s() got ref (considered const)",FCT);

   if (SIZE==B.SIZE) { if (!bfac) { C=*this; } else {
      size_t i=0, n=numel(); const T *b=B.data; T* c;

      C.SIZE=SIZE;                   
      C.NEW_DATA(SIZE.prod(0),NULL,0,0); 
      c=C.data;

      if (bfac==T(+1)) for (; i<n; ++i) { c[i] = data[i] + b[i]; } else
      if (bfac==T(-1)) for (; i<n; ++i) { c[i] = data[i] - b[i]; }
      else             for (; i<n; ++i) { c[i] = data[i] + b[i]*bfac; }
   }}
   else if (!B && !bfac) { C.init(*this); } 
   else if (!iflag || (SIZE.len && B.SIZE.len)) {
      wblog(FL,"ERR %s() size mismatch %s / %s",FCT,SSTR_(this),SSTR(B));
   }
   else {
      if (SIZE.len) { C=*this; } 
      else {
         if (!bfac) { C.init(B.SIZE); } else
         if (bfac==T(1)) { C.init(B); }
         else {
            size_t i=0, n=B.numel(); const T* b=B.data; T* c;
            C.SIZE=B.SIZE;
            C.NEW_DATA(B.SIZE.prod(0),NULL,0,0); 
            c=C.data;

            if (bfac==T(-1))
                 { for (; i<n; ++i) { c[i] =-b[i]; }}
            else { for (; i<n; ++i) { c[i] = b[i]*bfac; }}
         }
      }
   }
   return C;
};

template<class T>
template<class TB>
wbarray<T>& wbarray<T>::Plus(
   const wbarray<TB> &B, TB bfac,
   char iflag, 
   T afac      
){
   if (isEmpty()) {
      if (B.isEmpty()) { return *this; }
      if (iflag || afac==T(0)) {
         if (iflag && afac && afac!=T(1)) { wblog(FL,
            "ERR %s() got empty A with afac=%s",FCT,NSTR(afac));
         }
         if (bfac==TB(1)) { return init(B); } else
         if (!bfac) {  return init(B.SIZE); } 
         else {
            SIZE=B.SIZE; 
            NEW_DATA(B.SIZE.prod(0),NULL,0,0); 
         }
      }
   }
   else { iflag=0; }

   if (isRef()) wblog(FL,"ERR %s() got ref (considered const)",FCT);
   if (SIZE!=B.SIZE) wblog(FL,
      "ERR %s() size mismatch %s / %s",FCT,SSTR_(this), SSTR(B));

   Wb::cpyRange(data, B.data, numel(), iflag? T(0) : afac, bfac);

   return *this;
};

template<class T> inline
char wbarray<T>::sameUptoFac(
  const wbarray<T> &B, T* fac_, double eps
) const {

   size_t i, n=numel(); T a,fac=1;

   if (SIZE!=B.SIZE) return 1;
   if (!n) return 0;

   a=this->aMax(&i);
   if (a<=eps) {
      if (a==0 || B.aMax(&i)>eps) { if (fac_) (*fac_)=0; return 0; }
      else return 1;
   }
   if (Wb::abs(B.data[i])<=eps) return 2;

   fac=data[i]/B.data[i]; if (fac_) (*fac_)=fac;

   for (i=0; i<n; ++i) {
      a=Wb::abs(data[i]-fac*B.data[i]);
      if (a>eps) return 3;
   }
   return 0;
};

template<class T> inline
bool wbarray<T>::sameAs(const wbarray &B, double eps) const {
   if (this!=&B) {
      if (SIZE!=B.SIZE) { return 0; }
      if (data!=B.data) {
         if (eps) {
            for (size_t i=0, n=numel(); i<n; ++i) {
            if (Wb::abs(data[i]-B.data[i])>eps) { return 0; }}
         }
         else { return (!memcmp(data, B.data, numel()*sizeof(T))); }
      }
   }
   return 1;
};

template<class T> inline
wbvector<T>& wbarray<T>::getCol(size_t k, wbvector<T> &v) const {
   const T* x=col(k); 
   return v.init(SIZE[0],x); 
};

template<class T> inline
wbvector<T>& wbarray<T>::getRow(size_t k, wbvector<T> &v) const {
   const T* x=row(k); 
   return v.initStride(SIZE[1],x,SIZE[0]); 
};

template<class T> inline
wbarray<T>& wbarray<T>::getCol(size_t k, wbarray<T> &a) const {
   const T* x=col(k); 
   size_t s[2]={ SIZE[0], 1 }; a.SIZE.init(2,s); a.NEW_DATA(s[0],x);
   return a;
};

template<class T> inline
wbarray<T>& wbarray<T>::getRow(size_t k, wbarray<T> &a) const {
   const T* x=row(k); 
   size_t s[2]={ 1, SIZE[1] }; a.SIZE.init(2,s);

   a.NEW_DATA(s[1],NULL,0,0); 
   Wb::cpyStride(a.data, x,1,NULL,s[1],-1,SIZE[0],0); 

   return a;
};

template<class T> inline
wbvector<T> wbarray<T>::getDiag() const {

   wbvector<T> x; 

   if (!SIZE.len) { return x; }
   if (SIZE.len!=2) wblog(FL,
      "ERR %s() invalid matrix %s",FCT,SSTR_(this));

   size_t i=0, j=0, m=SIZE[0]+1, n=MIN(SIZE[0],SIZE[1]);
   if (n) { x.init(n);
      for (; i<n; ++i, j+=m) { x[i]=data[j]; }
   }
   return x;
};

template<class T>
void wbarray<T>::getMatSize(
   const char *F, int L, size_t &dim1, size_t &dim2) const {

   size_t i, r=SIZE.len, r2=r/2;
   const size_t *const &s=SIZE.data;

   if (r%2) wblog(F_L,
      "ERR %s() even rank object required (%s)",FCT,SSTR_(this));
   if (r) {
      dim1=s[0]; dim2=s[r2];
      for (i=1; i<r2; i++) { dim1*=s[i]; dim2*=s[i+r2]; }
   } else { dim1=dim2=0; }
};

template<class T>
void wbarray<T>::getMatSize(size_t &dim1, size_t &dim2, unsigned m) const {

   unsigned i=0, r=SIZE.len;
   const size_t *const &s=SIZE.data;

   if (m>r) { m=r; }
   if (m)   { dim1=s[0]; for (++i; i<m; ++i) {  dim1*=s[i]; }}
   else     { dim1=(r?1:0); }

   if (i<r) { dim2=s[i]; for (++i; i<r; ++i) {  dim2*=s[i]; }}
   else     { dim2=(r?1:0); }
};

template<class T>
void wbarray<T>::getMatSize(const ctrIdx &ic,
   size_t &dimc,           
   size_t &dimk,           
   unsigned r, size_t *dom 
) const {

   unsigned i,j;

   if (!ic) wblog(FL,"ERR %s() got empty ctrIdx %s",FCT,STR(ic));
   if (int(r)<0) {
      if (dom) wblog(FL,"ERR %s() invalid usage (%d,%p)",FCT,r,dom);
      r=SIZE.len;
   }
   else if (r>SIZE.len) {
      sprintf_str("%s() rank out of bounds (r=%d/%ld)",FCT,r,SIZE.len);
      wblog(FL,dom ? "ERR %s":"WRN %s", str);
      r=SIZE.len;
   }

   wbvector<char> mm(r); 
   dimc=dimk=1;

   for (j=0; j<ic.len; ++j) { i=ic[j];
      if (i>=r) { char s[128];
         snprintf(s,127,"%s() ctrIdx out of bounds (%s/%d)",FCT,STR(ic),r);
         if (dom || i>99) 
              wblog(FL,"ERR %s",s);
         else wblog(FL,"WRN %s",s);
      }
      else {
         if (++mm[i]!=1) wblog(FL,
            "ERR %s() non-unique ctrIdx (%s/%d)",FCT,STR(ic),r);
         dimc *= SIZE[i];
      }
   }

   for (i=0; i<r; ++i) { if (!mm[i]) { dimk *= SIZE[i]; }}

   if (dom) { for (*dom=1; i<SIZE.len; ++i) { *dom *= SIZE[i]; }}
};

template<class T>
char wbarray<T>::fitsSize(const wbarray<T> &B, char strict) const {

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

template<class T>
template<class TB>
bool wbarray<T>::sameSize(const wbarray<TB> &B,
   char lflag 
 ) const {

   if (!lflag) { return (SIZE==B.SIZE); }
   if (numel()!=B.numel()) { return 0; } 
   if (!data) { return 1; } 

   unsigned i=0, j=0, ra=SIZE.len, rb=B.SIZE.len;
   const size_t *sa=SIZE.data, *sb=B.SIZE.data;

   if (!ra || !rb) wblog(FL,"ERR %s() got rank r=%d/%d",FCT,ra,rb);

   if (lflag=='l') { lflag=3; } else 
   if (lflag=='L') { lflag=7; } else
   if (lflag<1 || lflag>7) { wblog(FL,
      "ERR %s() invalid lflag=%s",FCT,cSTR(lflag)); }

   if (lflag&2) { 
      while (ra) { if (sa[ra-1]==1) { --ra; } else { break; }}
      while (rb) { if (sb[rb-1]==1) { --rb; } else { break; }}
   }

   if (lflag&1) { 
      for (; i<ra; ++i) { if (sa[i]!=1) { break; }}
      for (; j<rb; ++j) { if (sb[j]!=1) { break; }}
   }

   if (lflag&4) { 
      while (i<ra && j<rb) {
         if (sa[i]!=sb[j]) { return 0; }
         for (++i; i<ra; ++i) { if (sa[i]!=1) break; }
         for (++j; j<rb; ++j) { if (sb[j]!=1) break; }
      }
      if (i<ra || j<rb) { return 0; }
   }
   else {
      ra-=i; rb-=j; if (ra!=rb) { return 0; }
      sa+=i; sb+=j;
      for (i=0; i<ra; ++i) { if (sa[i]!=sb[i]) return 0; }
   }

   return 1;
};

template<class T> inline
bool wbarray<T>::sameSize(const size_t *s, unsigned n) const {

   unsigned i=0, m=MIN(n,unsigned(SIZE.len));
   for (i=0; i<m; ++i) { if (s[i]!=SIZE.data[i]) return 0; }
   for (   ; i<n; ++i) { if (s[i]!=1) return 0; }
   return 1;
};

template<class T> inline
bool wbarray<T>::sameSize1(const size_t *s2, unsigned n2) const {

   const size_t *s1=SIZE.data;
   unsigned i=0, j=0, n1=SIZE.len;

   while (1) {
      for (; i<n1 && s1[i]==1; ++i) {};
      for (; j<n2 && s2[j]==1; ++j) {};
      if (i==n1 || j==n2 || s1[i]!=s2[j]) break;
      ++i; ++j;
   }

   return ((i<n1 || j<n2 || (!n1)^(!n2)) ? 0 : 1);
};

template<class T> inline
bool wbarray<T>::equal (const wbarray<T> &B, T eps) const {
    size_t i,s=numel(); if (SIZE!=B.SIZE) { return 0; }
    if (eps)
         { for (i=0; i<s; ++i) { if (fabs(data[i]-B.data[i])>eps) return 0; }}
    else { if (memcmp(data,B.data,s*sizeof(T))) { return 0; }}
    return 1;
};

template<class T> inline
bool wbarray<T>::unequal (const wbarray<T> &B, T eps) const {
    return !equal(B,eps);
};

template<class T> inline
bool wbarray<T>::equal (
    const wbarray<T> &B, double eps, double &maxdiff) const {

    if (SIZE!=B.SIZE) return 0;

    size_t i,s=numel(); double d;

    for (maxdiff=0, i=0; i<s; i++) {
        d = Wb::abs( data[i] - B.data[i] );
        if (maxdiff<d) maxdiff=d;
    }

    return (maxdiff<eps);
};

template<class T> inline
bool wbarray<T>::unequal (
    const wbarray<T> &B, double eps, double &maxdiff) const {
    return !(*this).equal(B,eps,maxdiff);
};

template<class T>
bool wbarray<T>::isDiag_aux(const T eps, const char* task) const {

   size_t i=0;
   wbIndex I(SIZE); 

   char isDiag=!strcmp(task,"isDiag");
   char isIdty=!strcmp(task,"isIdty");

   if (!isDiag && !isIdty) wblog(FL,"ERR %s() invalid task '%s'",FCT,task);
   if (SIZE.len%2) wblog(FL,"ERR %s() for even-rank objects only (%s)",
     task,SSTR_(this));

   while (++I) {
       if (!I.isdiag()) {
          if (Wb::abs(data[i])>eps) return 0;
       }
       else {
          if (isIdty && (Wb::abs(data[i]-T(1)))>eps) return 0;
       }
       ++i;
   }

   return 1;
};

template<class T>
bool wbarray<T>::isDiag_aux(double *epsp, const char* task) const {
   if (!epsp) wblog(FL,"ERR %s() got null eps",FCT);
   if (!task) wblog(FL,"ERR %s() got null task",FCT);

   size_t i,k, s=numel(), r=SIZE.len, l=r-1;
   WBINDEX I(r); bool is=1; double a, eps=0, eref=*epsp;

   char isDiag=!strcmp(task,"isDiag");
   char isIdty=!strcmp(task,"isIdty");
   if (!isDiag && !isIdty) wblog(FL,"ERR invalid task '%s'", task);

   if (!isRank(2)) { info("this"); wblog(FL,
   "ERR %s only appies to rank-2 tensors (%s).",task,SSTR_(this)); }

   for (i=0; i<s; i++) {
       if (I[0]!=I[1]) { a=Wb::abs(data[i]);
          if (is && eref<a) is=0;
          if (eps<a) eps=a;
       }
       else if (isIdty) { a=Wb::abs(data[i]-T(1));
          if (is && eref<a) is=0;
          if (eps<a) eps=a;
       }

       k=0; I[0]++;  
       while (I[k]>=SIZE[k] && k<l) { I[k]=0; ++I[++k]; }
   }

   *epsp=eps; return is;
};

template<class T>
bool wbarray<T>::isProptoId(T &x, const T eps) const {

   if (!data) { return 0; } 

   size_t i=0;
   wbIndex I(SIZE); 

   if (SIZE.len%2) wblog(FL,
      "ERR %s() for even-rank objects only (%s)",FCT,SSTR_(this));
   if (!x && data) { x=data[0]; }

   while (++I) {
      if (!I.isdiag())
           { if (Wb::abs(data[i]  )>eps) return 0; }
      else { if (Wb::abs(data[i]-x)>eps) return 0; }; ++i;
   }

   if (i && !data) wblog(FL,
      "ERR %s() invalid empty array (len=%d)",FCT,i);

   return 1;
};

template<class T>
bool wbarray<T>::isOpS(size_t *n) const {

   size_t i,r=SIZE.len, r2=r/2;

   if (r==0) return 1;
   if (r%2) return 0;

   for (i=0; i<r2; i++) if (SIZE[i]!=SIZE[i+r2]) return 0;

   if (n) { (*n)=Wb::prodRange(SIZE.data,r2); }

   return 1;
}

template<class T>
wbarray<T>& wbarray<T>::balanceOp(
   const char *F, int L, T &dref, double &dscale
){
   size_t i,j,l,n=0;

   if (!isOpS(&n)) wblog(F_L,
   "ERR %s() requires operator (%s)",FCT,SSTR_(this));

   if (n<=1) { dref=0; dscale=1; return *this; }
   dref=0; dscale=0;

   for (l=i=0; i<n; i++, l+=(n+1)) dref+=data[l];
   dref /= n;

   for (l=i=0; i<n; i++)
   for (j=0; j<n; j++, l++) { if (i==j) data[l]-=dref;
       dscale += Wb::abs2(data[l]);
   }
   dscale = std::sqrt(dscale)/n;

   if (dscale!=0) (*this)*=(1./dscale);
   else dscale=1;

   return *this;
}

template<class T>
wbarray<T>& wbarray<T>::Symmetrize(
   const char *F, int L,
   double *delta,  
   char cflag,     
   char tflag,     
   char fflag
){
   size_t i,j,n=0; T d1,d2; double e,E=0,dmax=0;

   if (!isOpS(&n)) wblog(F,L,
      "ERR %s() (generalized) square matrix required (%s)",
       FCT,SSTR_(this));
   if (tflag && !delta) return *this; 

   if (cflag && typeid(T)!=typeid(wbcomplex)) cflag=0;

   for (i=0; i<n; i++)
   for (j=i; j<n; j++) { d1=data[i+j*n]; d2=data[j+i*n];
      dmax=MAX(MAX(dmax,Wb::abs(d1)),Wb::abs(d2));
      if (i!=j) {
         e=(cflag ? Wb::abs(d1-Wb::CONJ(d2)) : Wb::abs(d1-d2));
         E=MAX(E,e);
      }
   }

   dmax=MAX(1.,dmax);
   e=E; if (dmax!=0) e/=dmax;

   if (delta) {
      if (dmax!=0) {
         if (e>(*delta)) { tflag=1;
#ifdef MATLAB_MEX_FILE
         this->put(FL,"M_");
#endif
            sprintf_str("%s() non-symmetric operator (%d)\n"
              "eps=%.3g/%.3g @ dmax=%.3g",FCT,cflag,e,*delta,dmax);

            if (fflag) wblog(F,L,"ERR %s",str);
            else wblog(F,L,"WRN %s",str);
         }
         *delta=(e);
      }
      else if (e==0) *delta=0; else *delta=1;
   }

   if (tflag==0) {
      if (cflag) {
         for (i=0; i<n; i++)
         for (j=i+1; j<n; j++) {
            d1=0.5*(data[i+j*n]+Wb::CONJ(data[j+i*n]));
            data[i+j*n]=d1; data[j+i*n] = Wb::CONJ(d1);
         }
      }
      else {
         for (i=0; i<n; i++)
         for (j=i+1; j<n; j++) {
            d1=0.5*(data[i+j*n]+data[j+i*n]);
            data[i+j*n] = data[j+i*n] = d1;
         }
      }
   }

   return *this;
};

template<class T> 
wbarray<T>& wbarray<T>::swapRows(size_t i1, size_t i2) {

   if (rank()!=2) wblog(FL,
      "ERR %s() for matrices only (%s)",FCT,SSTR_(this));
   if (i1>=SIZE[0] || i2>=SIZE[0]) wblog(FL,"ERR %s() "
      "index out of bounds (%d,%d; %s)",FCT,i1+1,i2+1,SSTR_(this));

   if (i1!=i2) { size_t dim1=SIZE[0], dim2=SIZE[1]; if (dim2) {
      T x, *d1=data+i1, *d2=data+i2; 
      for (size_t j=0; j<dim2; ++j, d1+=dim1, d2+=dim1) {
         x=(*d1); *d1=(*d2); *d2=x;
      }
   }}
   return *this;
};

template<class T> 
wbarray<T>& wbarray<T>::swapCols(size_t j1, size_t j2) {

   if (rank()!=2) wblog(FL,
      "ERR %s() for matrices only (%s)",FCT,SSTR_(this));
   if (j1>=SIZE[1] || j2>=SIZE[1]) wblog(FL,"ERR %s() "
      "index out of bounds (%d,%d; %s)",FCT,j1+1,j2+1,SSTR_(this));

   if (j1!=j2) { size_t dim1=SIZE[0]; if (dim1) {
      T x, *d1=data+dim1*j1, *d2=data+dim1*j2; 
      for (size_t i=0; i<dim1; ++i) {
         x=d1[i]; d1[i]=d2[i]; d2[i]=x;
      }
   }}
   return *this;
};

template<class T>
wbarray<T>& wbarray<T>::setCol(size_t k, size_t k0) {
   if (rank()!=2) wblog(FL,
      "ERR %s() for matrices only (%s)",FCT,SSTR_(this));
   if (k>=SIZE[1] || k0>=SIZE[1]) wblog(FL,
      "ERR %s() index out of bounds (%d,%d; %s)",
       FCT,k+1,k0+1,SSTR_(this));

   size_t dim1=SIZE[0];
   if (k!=k0 && dim1) {
      MEM_CPY<T>(ref(0U,k),dim1,ref(0U,k0)); 
   }

   return *this;
};

template<class T>
wbarray<T>& wbarray<T>::setCol(
   size_t k, const T* v, T fac, size_t stride) {

   if (rank()!=2) wblog(FL,
      "ERR %s() for matrices only (%s)",FCT,SSTR_(this));
   if (k>=SIZE[1]) wblog(FL,
      "ERR %s() index out of bounds (%d; %s)",FCT,k+1,SSTR_(this));

   size_t i=0, dim1=SIZE[0]; T *d=data+k*dim1; 
   if (dim1 && d!=v) {
      if (fac== 1) { for (; i<dim1; ++i, v+=stride) { d[i]= (*v); }} else
      if (fac==-1) { for (; i<dim1; ++i, v+=stride) { d[i]=-(*v); }}
      else         { for (; i<dim1; ++i, v+=stride) { d[i]=fac*(*v); }}
   }
   return *this;
};

template<class T>
wbarray<T>& wbarray<T>::setCol(size_t k, const wbvector<T> &v) {

   if (rank()!=2) wblog(FL,
      "ERR %s() for matrices only (%s)",FCT,SSTR_(this));
   if (k>=SIZE[1]) wblog(FL,
      "ERR %s() index out of bounds (%d: %s)",FCT,k+1,SSTR_(this));

   size_t dim1=SIZE[0];
   if (v.len!=dim1) wblog(FL,"ERR size mismatch (%d/%d)",v.len,dim1);
   if (dim1) MEM_CPY<T>(data+k*dim1, dim1, v.data); 
   return *this;
};

template<class T>
template<class T2>
wbarray<T>& wbarray<T>::setCol(size_t k, const wbvector<T2> &v) {

   if (rank()!=2) wblog(FL,
      "ERR %s() for matrices only (%s)",FCT,SSTR_(this));
   if (k>=SIZE[1]) wblog(FL,
      "ERR %s() index out of bounds (%d: %s)",FCT,k+1,SSTR_(this));

   size_t dim1=SIZE[0]; T* d=data+k*dim1; 
   if (dim1) {
      size_t i=0, n=(v.len<dim1 ? v.len: dim1);
      for (; i<n;    ++i) { d[i]=T(v.data[i]); }
      for (; i<dim1; ++i) { d[i]=T(0); }
   }
   return *this;
};

template<class T>
wbarray<T>& wbarray<T>::setRow(size_t k, const T* d0) {

   if (rank()!=2) wblog(FL,
   "ERR %s() for matrices only (%s)",FCT,SSTR_(this));
   if (k>=SIZE[0]) wblog(FL,
   "ERR %s() index out of bounds (%d; %s)",FCT,k+1,SSTR_(this));

   size_t i=0, dim1=SIZE[0], dim2=SIZE[1];
   T* d=data+k; 
   for (; i<dim2; i++, d+=dim1) d[i]=d0[i];

   return *this;
};

template<class T>
wbvector<T>& wbarray<T>::colNorm2(wbvector<T> &a) const {

   if (rank()!=2) wblog(FL,
      "ERR %s() for matrices only (%s)",FCT,SSTR_(this));

   size_t i,j, dim1=SIZE[0], dim2=SIZE[1];
   const T* d0=data; T x2;

   a.init(dim2);
   for (j=0; j<dim2; j++, d0+=dim1) { 
      for (x2=T(), i=0; i<dim1; i++) x2+=(Wb::CONJ(d0[i])*d0[i]);
      a[j]=x2;
   }
   return a;
};

template<class T>
T wbarray<T>::colNorm2(size_t k) const {

   size_t j=0, dim1=SIZE[0], dim2=SIZE[1];
   const T* d=data+k*dim1; 
   T a=0;

   if (rank()!=2) wblog(FL,
      "ERR %s() for matrices only (%s)",FCT,SSTR_(this));
   if (k>=dim2) wblog(FL,
      "ERR %s() index out of bounds (%d/%d; %s)",FCT,k,dim2,SSTR_(this));

   for (; j<dim2; ++j) { a+=Wb::CONJ(d[j])*d[j]; }
   return a;
};

template<class T>
T wbarray<T>::NormalizeCol(size_t k, char tnorm, char qflag) {

   if (rank()!=2) wblog(FL,
      "ERR %s() matrix type required (%s)",FCT,SSTR_(this));

   const size_t dim1=SIZE[0], dim2=SIZE[1];
   if (k>=dim2) wblog(FL,
      "ERR %s() index out of bounds (%d/%d)",FCT,k,dim2);

   T *d=data+dim1*k, x2=Wb::overlap(d,d,dim1,1,tnorm); 
   x2=Wb::sqrt(x2);

   if (x2!=0) Wb::timesRange(d,1/x2,dim1);
   else if (!qflag) wblog(FL,"ERR %s() got vector with norm 0!",FCT);

   return x2;
};

template<class T>
wbarray<T>& wbarray<T>::ColProject(
   size_t k1, size_t k2, char nflag, char tnorm
){

   if (rank()!=2) wblog(FL,
      "ERR %s() matrix type required (%s)",FCT,SSTR_(this));

   const size_t dim1=SIZE[0], dim2=SIZE[1];
   if (k1>=dim2 || k2>=dim2) wblog(FL,
      "ERR %s() index out of bounds (%d,%d/%d)",FCT,k1,k2,dim2);
   if (k1==k2) wblog(FL,"WRN %s() got k1=k2=%d",FCT,k1+1);

   gs_project_range(
      data+dim1*k1, data+dim1*k2, dim1, 1, 
      nflag, tnorm
   );
   return *this;
};

template<class T>
bool wbarray<T>::isOrthogonalCol(size_t k0, char tnorm) const {

   if (rank()!=2) wblog(FL,
      "ERR %s() matrix type required (%s)",FCT,SSTR_(this));

   const size_t dim1=SIZE[0], dim2=SIZE[1];
   if (k0>=dim2) wblog(FL,
      "ERR %s() index out of bounds (%d/%d)",FCT,k0,dim2);
   if (!dim1) return 0;

   const T *d=data, *d0=data+dim1*k0;
   T nrm2=Wb::overlap(d0,d0,dim1,1,tnorm);

   for (size_t k=0; k<dim2; k++, d+=dim1) { 
      if (k!=k0 && Wb::abs(Wb::overlap(d,d0,dim1,1,tnorm)/nrm2)>1E-10)
      return 0;
   }

   return 1;
};

template<class T>
bool wbarray<T>::isOrthoCols(T *x2, char tnorm, T eps) const {

   bool q=0;
   const size_t *s=SIZE.data;

   if (SIZE.len==2 && s[1] && s[0]>=s[1]) {
      wbarray<T> E;
      Wb::MatProd(*this,*this,E, tnorm? 'T':'C');

      if (x2)
           { q=E.isProptoId(*x2,eps); }
      else { q=E.isDiagMatrix(); }
   }

   return q;
};

template<class T>
unsigned wbarray<T>::QRdecomp(
   const char *F, int L, wbarray<T> &Q, wbarray<T> &R, T eps
){
   if (isEmpty()) { Q.init(); R.init(); return 0; }
   if (rank()!=2 || !numel()) wblog(F_L,
      "ERR %s() matrix type required (%s)",FCT,SSTR_(this));

   size_t dim1=SIZE[0], dim2=SIZE[1], n=MIN(dim1,dim2);

   if (dim1<=1 || dim2<=1) {
      Q.init(dim1,n); R.init(n,dim2);

      if (!dim1 || !dim2) return n;
      if (dim1==1) { Q[0]=1; R=*this; return n; }
      if (dim2==1) { Q=*this;
         if ((R[0]=this->norm())!=0) { Q*=(T(1)/R[0]); }
         else { Q[0]=1; }
         return n;
      }
   }

   if (allEqual(0)) {
      Q.init(dim1,1); R.init(1,dim2); Q[0]=1;
      return 1;
   }

   size_t i,j,k=0;
   WBPERM P;  
   wbarray<T> U(dim1,n); 

   R=*this;

   for (k=0; k<n; ++k) {
      if (R.Householder(F_L,k,U.col(k),&P,eps)<0) break;
   }

if (!R.isFinite() || !U.isFinite()) wblog(FL,"ERR %s() ",FCT);

   if (k<n) {
      if (!k) wblog(FL,"ERR %s() resulted in k=%d/%d !?",FCT,k,n);
      R.Resize(k,dim2); n=k;
   }
   else if (n<dim1) {
      R.Resize(n,dim2);
   }

   Q.initIdentityB(dim1,n); 

   T zero=0, *x;
   const T *u=U.col(n-1);

   for (k=n; k>0; u-=dim1) { --k; x=Q.data;
      for (j=0; j<n; ++j, x+=dim1) {
         Wb::householder(u+k,x+k,dim1-k,eps);
      }
   }

   x=Q.data;
   for (j=0; j<n; ++j, x+=dim1) { 
      for (i=0; i<dim1; ++i) { if (Wb::abs(x[i])>eps) {
         if (x[i]<zero) { T *y=R.data+j;
            for (i=0; i<dim1; ++i) { x[i]=-x[i]; }
            for (k=0; k<dim2; ++k, y+=n) { y[0]=-y[0]; }
         }
         break;
      }}
   }

   return n;
};

template<class T>
int wbarray<T>::Householder( 
   const char *F, int L, size_t k, 
   T *u, WBPERM *P, 
   T eps
){
   if (rank()!=2) wblog(F_L,
      "ERR %s() matrix type required (%s)",FCT,SSTR_(this));

   size_t i=0, j=0, k2=k, dim1=SIZE[0], dim2=SIZE[1];
   if (k>=MIN(dim1,dim2)) wblog(F_L,"ERR %s() "
      "index out of bounds (%d; %s)",FCT,k+1,SSTR_(this));

   T nrm, x2=0., x0max=0., *x; 

   for (; i<k; ++i) u[i]=0.;

   if (P) {
      if (!k) P->init(dim2); else
      if (P->len!=dim2) wblog(FL,"ERR %s() "
         "got invalid P for (k=%d: %d/%d) !?",FCT,k+1,P->len,dim2);
      T x2max=0.;

      for (x=data, j=0; j<dim2; ++j, x+=dim1) { 
         if (!(P->data[j])) {
            for (x2=0., i=k; i<dim1; ++i) { x2 += (Wb::CONJ(x[i])*x[i]); }
            if (x2max<x2) { x2max=x2; k2=j; }
         }
         else {
            if (P->data[j]>k) wblog(FL, 
               "ERR %s() invalid P=[%s] (%d/%d/%d)",FCT,STR_(P),j,k,dim2);
            if (x0max<(x2=Wb::abs(x[P->data[j]-1]))) { x0max=x2; }
         }
      }; x2=x2max;

      x=data+k2*dim1;    

      for (i=k; i<dim1; ++i) { u[i]=x[i]; };
   }
   else {
      for (x=data, j=0; j<k; ++j, x+=dim1) {
         if (x[j]<0) wblog(FL,"ERR %s() x[%d]=%d",FCT,j,x[j]);
         if (x0max<(x2=Wb::abs(x[j]))) { x0max=x2; }
      }
      for (x2=0., i=k; i<dim1; ++i) { u[i]=x[i];
         x2+=(Wb::CONJ(x[i])*x[i]);
      }
   }

   eps*=(x0max<x2 ? x2 : x0max);

   nrm=Wb::sqrt(x2); if (nrm<=eps) { 
     if (!P && k+1<dim2) wblog(FL,
       "ERR %s() check cols (%d/%d)",FCT,k,dim2);
     return -1;
   }

   u[k] += (x[k]>0 ? +nrm : -nrm);

   nrm=T(1)/Wb::sqrt( x2 - Wb::CONJ(x[k])*x[k] + Wb::CONJ(u[k])*u[k] );
   if (!isfinite(nrm)) wblog(FL,"ERR %s() nrm=%g",FCT,nrm); 

   for (i=k; i<dim1; ++i) { u[i]*=nrm; }

   if (P) {
      for (x=data, j=0; j<dim2; ++j, x+=dim1) { 
         if (!(P->data[j])) {
         Wb::householder(u+k,x+k,dim1-k,eps); }
      }
      P->data[k2]=(k+1); 
   }
   else {
      for (x=data+k*dim1, j=k; j<dim2; ++j, x+=dim1) { 
         Wb::householder(u+k,x+k,dim1-k,eps);
      }
   }

   return k2;
};

template<class T> inline
void Wb::householder(const T *u, T *a, size_t n, const T& eps) {

    size_t i=0;
    T ua=0.; for (; i<n; ++i) { ua += (Wb::CONJ(u[i])*a[i]); }

    if (ua) { ua*=2.;
       for (i=0; i<n; ++i) { a[i] -= (u[i]*ua); } 
       if (eps>0.) { ua=0.; 
          for (i=0; i<n; ++i) ua += Wb::CONJ(a[i])*a[i];
          if (ua<eps) { for (i=0; i<n; ++i) a[i]=0.; }
       }
    }
};

template<class T>
unsigned wbarray<T>::eigTriDiag(
   const char *F, int L, wbvector<T> &E, wbarray<T> *U, T eps
 ) const {

   unsigned D;

   if (!isMatrix()) wblog(F_L,
      "ERR %s() matrix expected (%s)",FCT,SSTR_(this));

   D=dim(1); E.init(D); if (U) U->initIdentity(D);
   if (D<=1) {
      if (D) { E.data[0]=data[0]; }
      return 0;
   }

   wbvector<T> d_(D), e_(D);

   int i,k,m,l; unsigned iter=0;
   T s,r,p,g,f,da,ea,c,b, *d=d_.data, *e=e_.data;

   for (l=D-1, i=0; i<l; ++i) {
      if (Wb::abs(data[(i+1)+i*D] - data[i+(i+1)*D])>eps) wblog(F_L,
         "ERR %s() got non-symmetric matrix\nM(%d,%d+1): %.4g / %.4g !?",
         FCT,i+1,i+2, double(data[(i+1)+i*D]), double(data[i+(i+1)*D]));
      if (i>1 && (Wb::abs(data[i])>eps || Wb::abs(data[i*D])>eps)) wblog(F_L,
         "ERR %s() got non-triangular matrix\nM(%d,1): %.4g !?",
         FCT,i+1, double(data[i]));

      d[i]=data[i*D+i]; e[i]=0.5*(data[i*D+i+1] + data[(i+1)*D+i]);
   }; d[i]=data[i*D+i];

   for (l=0; l<int(D); ++l) { iter=0;
      do {
         for (m=l; m<int(D-1); ++m) {
            ea=Wb::abs(e[m]); da=Wb::abs(d[m]) + Wb::abs(d[m+1]);
            if (T(ea+da)==da) break;
         }
         if (m!=l) {
            if (++iter > 30) wblog(FL,
               "ERR %s() too many iterations (iter=%d) !?",FCT,iter);
            g=(d[l+1] - d[l]) / (T(2)*e[l]); 
            r=NR::pythag(g,T(1));
            g=d[m]-d[l]+e[l]/(g+NR::sign(r,g)); 
            s=c=1; p=0;

            for (i=m-1; i>=l; --i) {
               f=s*e[i];
               b=c*e[i]; e[i+1]=(r=NR::pythag(f,g));
               if (r==0) { d[i+1]-=p; e[m]=0; break; }

               s=f/r; 
               c=g/r; 
               g=d[i+1]-p;
               r=(d[i]-g)*s + 2*c*b;
               d[i+1]=g + (p=s*r);
               g=c*r-b;

               if (U) { wbarray<T> &z=(*U);
                  for (k=0; k<int(D); ++k) { f=z(k,i+1);
                     z(k,i+1) = s*z(k,i) + c*f;
                     z(k,i  ) = c*z(k,i) - s*f;
                  }
               }
            }

            if (r==0 && i>=l) continue;
            d[l]-=p; e[l]=g; e[m]=0;
         }
      } while (m!=l);
   }

   wbperm P;
   d_.sort(E,P);  
   if (U) U->ColPermute(P);

   return 0;
};

template<class T>
wbarray<T>& wbarray<T>::ColPermute(const wbperm &P, char iflag){

   if (rank()!=2) wblog(FL,
      "ERR %s() matrix type required (%s)",FCT,SSTR_(this));
   if (P.len!=SIZE[1]) wblog(FL,
      "ERR %s() size mismatch (%d/%d)",FCT,P.len,SIZE[1]);

   size_t i=0, dim1=SIZE[0], dim2=SIZE[1];
   wbarray<T> X(*this);

   for (; i<P.len; i++) {
      if (P[i]>=dim2) wblog(FL,
         "ERR %s() index out of bounds (%d/%d)",FCT,P[i],dim2);
      if (iflag==0) 
           MEM_CPY<T>(data+i*dim1, dim1, X.data+P[i]*dim1);
      else MEM_CPY<T>(data+P[i]*dim1, dim1, X.data+i*dim1);
   }

   return *this;
};

template<class T>
wbarray<T>& wbarray<T>::SignCol(size_t k, const T *d0, char tnorm) {

   if (rank()!=2) wblog(FL,
      "ERR %s() matrix type required (%s)",FCT,SSTR_(this));

   const size_t dim1=SIZE[0], dim2=SIZE[1];
   if (k>=dim2) wblog(FL,
      "ERR %s() index out of bounds (%d/%d)",FCT,k,dim2);

   T *d=data+dim1*k, x=Wb::overlap(d0,d,dim1,1,tnorm); 
   if (x<0) { 
      for (size_t i=0; i<dim1; i++) { d[i]=-d[i]; }
   }

   return *this;
};

template<class T>
wbarray<T>& wbarray<T>::FlipSignCol(size_t k) {

   if (rank()!=2) wblog(FL,
      "ERR %s() matrix type required (%s)",FCT,SSTR_(this));

   const size_t dim1=SIZE[0], dim2=SIZE[1];
   if (k>=dim2) wblog(FL,
      "ERR %s() index out of bounds (%d/%d)",FCT,k,dim2);

   T *d=data+dim1*k;
   for (size_t i=0; i<dim1; i++) { d[i]=-d[i]; }

   return *this;
};

template<class T>
wbarray<T>& wbarray<T>::SignConventionCol(
   size_t k0, 
   double deps  
){

   if (rank()!=2) wblog(FL,
      "ERR %s() matrix type required (%s)",FCT,SSTR_(this));

   const size_t dim1=SIZE[0], dim2=SIZE[1];
   size_t i,k,k1,k2; T *d, eps=
   ((typeid(T)!=typeid(double) && typeid(T)!=typeid(float)) ? 0 : deps);

   if (k>=dim2) wblog(FL,
      "ERR %s() index out of bounds (%d/%d)",FCT,k,dim2);

   if (int(k0)>=0) { k1=k0; k2=k0+1; }
   else { k1=0; k2=dim2; }

   for (k=k1; k<k2; k++) { d=data+dim1*k;
      for (i=0; i<dim1; i++) {
         if (d[i]> eps) break;
         if (d[i]<-eps) {
             for (i=0; i<dim1; i++) { d[i]=-d[i]; }
             break;
         }
      }
   }

   return *this;
};

template<class T>
wbarray<T>& wbarray<T>::NormalizeCols(
   const char *F, int L, double *amin, double *amax, char tnorm
){
   if (!SIZE.len) return *this;

   size_t i,j,m, n=SIZE[0];
   T *d=data, z=0; double a;
   wbvector<T> nn;

   if (SIZE.len>1) for (m=SIZE[1], i=2; i<SIZE.len; i++) m*=SIZE[i];
   else m=1;

   nn.init(m);

   for (j=0; j<m; j++, d+=n, z=0) {
      z=Wb::overlap(d,d,n,1,tnorm); z=Wb::sqrt(z);

      a=Wb::abs(z); if (a==0) { 
         wblog(F,L,"ERR |sum^2| of column returns %.3g!%s",a,
         tnorm && (typeid(T)==typeid(wbcomplex)) ? " (using t-norm)":"");
      }
      nn[j]=z;
   }

   if (amin) {
      a=nn.aMin(); if ((*amin)>0 && a<(*amin))
      wblog(F,L,"WRN %s() nMin=%.3g (%.3g)",FCT,a,*amin);
      (*amin)=a;
   }
   if (amax) {
      a=nn.aMax(); if ((*amax)>0 && a>(*amax))
      wblog(F,L,"WRN %s() nMax=%.3g (%.3g)",FCT,a,*amax);
      (*amax)=a;
   }

   for (d=data, j=0; j<m; j++, d+=n)
   for (z=T(1)/nn[j], i=0; i<n; i++) { d[i]*=z; }

   return *this;
};

template<class T>
wbarray<T>& wbarray<T>::OrthoNormalizeCols(
   const char *F, int L,
   char tnorm,  
   char qflag,  
   double eps,  
   unsigned np  
){
   if (!SIZE.len) return *this;

   size_t i,j,k,l,p,m, n=SIZE[0];
   T z, *d=0, *v=data; double a;

   char cflag=(typeid(T)==typeid(wbcomplex));
   char xflag=(qflag=='x' || qflag=='X');

   if (SIZE.len<2) m=1; else
   for (m=SIZE[1], i=2; i<SIZE.len; i++) m*=SIZE[i];

   if (!m || !n) return *this;

   for (l=k=0; k<m; ++k, ++l, v+=n) {
      for (p=0; p<np; ++p) {
         for (d=data, j=0; j<l; ++j, d+=n) {
            z=Wb::overlap(d,v,n,1,tnorm); if (z!=T(0)) {
            for (i=0; i<n; ++i) v[i]-=z*d[i]; } 
         }
      }

      z=Wb::overlap(v,v,n,1,tnorm);
      a=Wb::abs(z);

      if (a>eps) { a=Wb::sqrt(1/a);
         if (xflag && l<k)
              { for (i=0; i<n; ++i) d[i]=v[i]*a; }
         else { for (i=0; i<n; ++i) v[i]*=a; }
      }
      else {
         if (xflag) --l; 
         else if (qflag) {
            for (i=0; i<n; ++i) v[i]=0;
         }
         else { 
            wblog(F_L,"ERR |sum^2| of column returns %.3g!%s",
            a, tnorm  && cflag ? " (using t-norm)" : "");
         }
      }
   }

   if (xflag && l<k) {
      if (SIZE.len!=2) wblog(F_L,"ERR %s() qflag='%c' must be "
         "used with rank-2 tensors (%d)",FCT,qflag,qflag,SIZE.len);
      if (!l) wblog(F_L,"ERR %s() got all null-vectors (%g)",FCT,k);
      SIZE[1]=l;
   }

   return *this;
};

template<class T>
T wbarray<T>::norm2Cols(long j1, long j2) const {

   T x2=0; long i,j, dim1,dim2;
   if (SIZE.len!=2) wblog(FL,
      "ERR %s() got size %s array",FCT,SSTR_(this));
   dim1=SIZE[0]; dim2=SIZE[1];

   if (j1<0) { j1+=dim2; }
   if (j2<0) { j2+=dim2; }; if (j1>j2 || !dim1) { return x2; }
   if (j1<0 || j2>=dim2) wblog(FL,
      "ERR %s() index out of bounds (%ld,%ld/%ld)",FCT,j1,j2,dim2);

   for (j=j1; j<=j2; ++j)
   for (i=0; i<dim1; ++i) { x2+=Wb::norm2(data[i+j*dim1]); } 
   return x2;
};

template<class T>
T wbarray<T>::norm2Recs(long i1, long i2) const {

   T x2=0; long i,j,dim1,dim2;
   if (SIZE.len!=2) wblog(FL,
      "ERR %s() got size %s array",FCT,SSTR_(this));
   dim1=SIZE[0]; dim2=SIZE[1];

   if (i1<0) { i1+=dim1; }
   if (i2<0) { i2+=dim1; }; if (i1>i2 || !dim2) { return x2; }
   if (i1<0 || i2>=dim1) wblog(FL,
      "ERR %s() index out of bounds (%ld,%ld/%ld)",FCT,i1,i2,dim1);

   for (j=0; j<dim2; ++j)
   for (i=i1; i<=i2; ++i) { x2+=Wb::norm2(data[i+j*dim1]); } 
   return x2;
};

template<class T>
bool wbarray<T>::isSym_aux(
  const char *F, int L, const char *fct,
  const wbarray<T> &B, double eps, double* xref,
  const char symflag,
  const char lflag
) const {

   WBINDEX S1,S2;
   size_t i,j,k,dim1,dim2, r=SIZE.len, r2=r/2;
   double x;

   const wbarray<T> &A = (*this);
   bool issame = (&A==&B);
   wbarray<T> A2,B2;

   if (r%2 || r!=B.SIZE.len) {
      if (lflag) sprintf_str(
         "%s %s() only applies to even-rank objects (%ld;%ld)",
          shortFL(F,L),fct, r, B.SIZE.len);
      return 0;
   }
   if (isEmpty() && B.isEmpty()) return 1;

   dim1=dim2=1;
   for (i=0; i<r2; i++) { j=i+r2; dim1*=SIZE[i]; dim2*=SIZE[j];
      if (SIZE[i]!=B.SIZE[j] || SIZE[j]!=B.SIZE[i]) {
         if (lflag) sprintf_str( 
           "%s %s() size mismatch [%s; %s].", SHORT_FL, fct,
            SSTR_(this), SSTR(B));
         return 0;
      }
   }

   A.groupIndizes_DREF(A2,r2,r2);
   B.groupIndizes_DREF(B2,r2,r2);

   if (xref) *xref=0;
   if (eps<0) {
      eps=-eps; x=A2.aMax(); if (x>1) eps*=x; 
   }

   if (symflag=='s') { 
      if (eps==0) {
         for (i=0; i<dim1; i++) { k = (issame ? i : 0);
         for (j=k; j<dim2; j++) {
            if (A2(i,j)!=Wb::CONJ(B2(j,i))) {
               if (xref) {
                  x=Wb::abs(A2(i,j)-Wb::CONJ(B2(j,i)));
                  *xref = MAX(*xref,x);
               }
               else return 0;
            }
         }}
      }
      else {
         for (i=0; i<dim1; i++) { k = (issame ? i : 0);
         for (j=k; j<dim2; j++) {
            x=Wb::abs(A2(i,j)-Wb::CONJ(B2(j,i)));
            if (xref) *xref=MAX(*xref,x); else if (x>eps) return 0;
         }}
      }
   }
   else if (symflag=='a') { 
      if (eps==0) {
         for (i=0; i<dim1; i++) { k = (issame ? i : 0);
         for (j=k; j<dim2; j++) {
            if (A2(i,j)!=-Wb::CONJ(B2(j,i))) {
               if (xref) {
                  x = Wb::abs(A2(i,j) + Wb::CONJ(B2(j,i)));
                  *xref = MAX(*xref,x);
               }
               else return 0;
            }
         }}
      }
      else {
         for (i=0; i<dim1; i++) { k = (issame ? i : 0);
         for (j=k; j<dim2; j++) {
            x = Wb::abs(A2(i,j) + Wb::CONJ(B2(j,i)));
            if (x>eps) { if (xref) *xref=MAX(*xref,x); else return 0; }
         }}
      }
   }
   else wblog(FL,"ERR invalid symflag=%c<%d>",symflag,symflag);

   if (xref && (*xref)>eps) return 0;

   return 1;
}

template<class T>
char wbarray<T>::hasGroupSize(size_t s1, size_t s2) const {

   size_t i,s=1,s0=1;

   for (i=0; i<SIZE.len; i++) {
      s*=SIZE[i]; if (s0==s1 && (!s || s!=s1)) break;
      s0=s;
   }
   for (s=1; i<SIZE.len; i++) s*=SIZE[i];

   return ((s0!=s1 || s!=s2) ? 0 : 1);
};

template<class T> inline
bool wbarray<T>::isComplex() const { return 0; };

template<> inline
bool wbarray<wbcomplex>::isComplex() const {
   for (size_t s=numel(), i=0; i<s; i++)
   if (data[i].i!=0.) return 1;

   return 0;
};

template<class T> inline
bool wbarray<T>::isZero(double eps, char bflag) const {

   size_t i, s=numel(); if (!s) { return 1; }

   if (eps==0) {
      for (i=0; i<s; ++i) { if (data[i]!=0.) return 0; }
   }
   else {
      if (bflag) {
         double e2=eps*eps, n2=0.;
         for (i=0; i<s; ++i) {
            n2+=Wb::norm2(data[i]); if (n2>e2) { return 0; }
         }
      }
      else {
         for (i=0; i<s; ++i) { if (Wb::abs(data[i])>eps) { return 0; }}
      }
   }

   return 1;
}

template<class T> inline
size_t wbarray<T>::nnz(const T eps, WBINDEX *I) const {

   size_t i,n, s=numel(); T a;

   for (i=n=0; i<s; ++i) {
      a=data[i]; if (a<0) { a=-a; };
      if (a>eps) { ++n; }
   }

   if (I) {
      I->init(n);
      for (i=n=0; i<s; i++) {
         a=data[i]; if (a<0) a=-a;
         if (a>eps) I->data[n++]=i;
      }
   }

   return n;
}

template<class T> inline
void wbarray<T>::SetRef (
   const widx_t *I, size_t len, const T* d0
){
   size_t i,k,s, S=numel(), r=SIZE.len, l=len-1;
   size_t nz=r-len; 

   if (len>r || r==0 || (len && data==NULL)) wblog(FL,
      "ERR %s() index out of bounds ([%s]) for %s !?",FCT,
      (WBINDEX(len,I)+1).toStr().data, SSTR_(this));
   for (i=0; i<len; i++) if (I[i]>=SIZE[i+nz]) wblog(FL,
      "ERR %s() index out of bounds (%s; %s)",FCT,
      (WBINDEX(len,I)+1).toStr().data, SSTR_(this));

   if (!nz || !len) wblog(FL,"ERR invalid reference (len=%d/%d)",len,r);

   for (s=SIZE[0], i=1; i<nz; i++) s*=SIZE[i]; 
   for (k=I[l], i=l-1; i<l; i--) { k = k*SIZE[i+nz] + I[i]; }
   k*=s;

   if (k+s>S || d0==NULL) wblog(FL,
      "ERR index out of bounds ( %d+%d / %d; 0x%lX)\n"
      "having I=[%s] for %s array",k,s,S,d0,
      (WBINDEX(len,I)+1).toStr().data, SSTR_(this)
   );

   MEM_CPY<T>(data+k,s,d0);
};

template<class T> inline
wbarray<T>& wbarray<T>::Squeeze() { 

   size_t r=SIZE.len; if (r>1) {
      unsigned i=0, l=0; size_t *s=SIZE.data;
      for (; i<r; ++i) { if (s[i]!=1) {
         if (l!=i) { s[l]=s[i]; }
         ++l;
      }}
      if (l<r) {
         SIZE.len=(l? l:1); 
      }
   }
   return *this;
};

template<class T> inline
int wbarray<T>::skipSingletons(const char *F, int L,
   unsigned r2 
) {
   int m=0;
   unsigned i, r=SIZE.len;

   if (r2==r || (!r && r2<=0)) { return m; }
   if (!r) { return (m=-r2); } 

   if (int(r2)<=0) { r2=-r2; if (!r2) { r2=1; }
   if (r2<r) {
      for (i=r-1; i>=r2; --i) { if (SIZE[i]!=1) break; }

      m=r-1-i; 
      if (m) { SIZE.Shorten2(i+1); } 
   }}
   else if (r2<r) { int e=0;
      for (i=r2; i<r; ++i) { if (SIZE[i]!=1) { ++e; if (F) wblog(F_L,
         "ERR %s() non-singleton at %d/%d in [%s]",FCT,i+1,r2,STR(SIZE));
      }}
      if (e) { m=-e; }
      else   { m=r-r2; SIZE.Shorten2(r2); }
   }
   else { m=-r2; 
     if (F) wblog(F_L,"ERR %s() rank out of bounds (%d/%d)",FCT,r2,r);
   }

   return m;
};

template<class T> inline
wbarray<T>& wbarray<T>::appendSingletons(
   const char *F, int L, unsigned r, unsigned m) {

   unsigned i, l=SIZE.len;
   if (l<r) { unsigned s=1;
      if (!l) { s=0; 
      }
      SIZE.Resize(r); for (i=l; i<r; ++i) { SIZE[i]=s; }
   }
   else if (l>r) {
      size_t *s=SIZE.data;

      if (l==2 && r==1 && s[0]==1) {
         s[0]=s[1]; SIZE.len=1;
      }
      else {
         if (r>2) { r+=m; } else
         if (!r) wblog(FL,"ERR %s() got r=%d+%d",FCT,r,m);

         if (l>r) {
            for (i=r; i<l; ++i) { if (s[i]!=1) { wblog(F_L,
               "ERR %s() got r=%d->%d for %s",FCT,l,r,SSTR_(this));
            }}
            SIZE.len=r;
         }
      }
   }

   return *this;
};

template<class T> inline
void wbarray<T>::prependSingletons(unsigned r) {

   if (!SIZE.len) wblog(FL,"ERR %s() called on empty array",FCT);

   if (r>SIZE.len) {
      size_t i, l=SIZE.len-1, m=r-SIZE.len; 
      SIZE.Resize(r);
      for (i=l; i<r; --i) SIZE[i+m]=SIZE[i];
      for (i=0; i<m; ++i) SIZE[i]=1;
   }
   else if (r<SIZE.len) wblog(FL,
   "ERR %s() invalid extended rank !? (%d->%d)",FCT,SIZE.len,r);
};

template<class T> inline
wbarray<T>& wbarray<T>::ExpandOM(const char *F, int L,  unsigned r0,
    const wbvector<unsigned> &S,  
    const unsigned *sx  
){
    if (!SIZE) {
       if (!S && !sx) { return *this; }
       wblog(FL,"ERR %s() got empty vector with S=%s",FCT,STR(S));
    }
    if (!S) { if (!sx) {
       if (int(r0)<=0 || r0==SIZE.len) { return *this; }}
       wblog(FL,"ERR %s() got empty S=[%s] for %s",FCT,STR(S),SSTR_(this));
    }
    if (sx && S.isEqual(sx)) { sx=NULL; }

    unsigned i=0, l=SIZE.len;
    size_t M=1;

    if (int(r0)<=0) { r0=l-S.len;
       if (int(r0)<2) wblog(FL,
       "ERR %s() requesting %d/%d OM indices",FCT,S.len,l);
    }
    else if (r0+S.len<l) wblog(F_L,
       "ERR %s() rank out of bounds %d+%d / %d",FCT,r0,S.len,l);

    M=( r0<l ? numOM(r0) : 1 );

    if (M!=S.prod()) wblog(FL,
       "ERR %s() got OM size mismatch %s @ %d <> %s",FCT,
       SSTR_(this), r0, SSTR(S));

    if (r0<=2 && M>1) wblog(FL,"ERR %s() "
       "got OM for rank-%d array %s (l=%d)",FCT,r0,SSTR_(this), S.len);

    if (sx) { M=1;
       for (; i<S.len; ++i) { M*=sx[i]; if (sx[i]!=S[i]) {
          if (sx[i]<S[i]) { wblog(FL,"ERR %s() "
             "cannot zero-pad to smaller size\nS=[%s] -> [%s]",
             FCT,STR(S),STR(wbvector<unsigned>(S.len,sx,'r')));
          }
       }}
       if (M>1 && r0<=2) { wblog(FL,"ERR %s() " 
          "got OM for rank-%d array %s (M=%ld)",FCT,r0,SSTR_(this),M);
       }
    }

    if (l!=r0+1 && l!=r0+S.len) { unsigned j=0;
       for (i=r0; i<l && j<S.len; ++i, ++j) { if (SIZE[i]!=S[j]) {
          while (S[j]==1 && ++j<S.len) {
             if (SIZE[i]==S[j]) break; 
          }
       }}
       for (; i<l; ++i) {  if (SIZE[i]!=1) break; }
       for (; j<S.len; ++j) { if (S[j]!=1) break; }

       if (i<l || j<S.len) { wblog(FL,
          "WRN %s() got OM size reshape %s @ %d <> %s",
          FCT, SSTR_(this), r0, SSTR(S));
       }
    }

    wbvector<size_t> SX; SIZE.resize(r0+S.len,SX);
    for (i=0; i<S.len; ++i) { SX[r0+i]=S[i]; }
    Reshape(SX);

    if (sx) {
       for (i=0; i<S.len; ++i) { SX[r0+i]=sx[i]; }
       Resize(SX);
    }

    return *this;
};

template<class T> inline
int wbarray<T>::ExpandOM(
   const char *F, int L, wbarray<T> &B,
   unsigned r 
) {

   int rval=0, eq=0;

   if (isEmpty() ||  B.isEmpty()) { return (rval=-1); } 

   if (!r) { unsigned R=MAX(SIZE.len,B.SIZE.len);
      if (!isScalar() || !B.isScalar() || R>2) { wblog(FL,
         "ERR %s() unexpected empy/scalar QSpaces (%s / %s, r=%d)",
         FCT,SSTR_(this),SSTR(B),R);
      }
      if (  SIZE.len<R) {   appendSingletons(R); }
      if (B.SIZE.len<R) { B.appendSingletons(R); }
      return rval;
   }

   if ((eq=SIZE.isEqual(B.SIZE))) {
     #ifndef WB_SKIP_ASSERT
      if (r<=2 && r<SIZE.len) { unsigned i=r; 
         for (; i<SIZE.len; ++i) { if (SIZE[i]!=1) { wblog(FL,
            "ERR %s() got OM space %s for rank r=%d",FCT,SSTR_(this),r);
         }}
      }
     #endif
      return rval; 
   }

   if (SIZE.len<r || B.SIZE.len<r) wblog(F_L,"ERR %s() "
      "rank out of bounds (r=%d/%d/%d)",FCT,r,SIZE.len,B.SIZE.len);
   if (memcmp(SIZE.data,B.SIZE.data,r*sizeof(SIZE[0]))) wblog(F_L,
      "ERR %s() size mismatch %s / %s (r=%d)",FCT,SSTR_(this),SSTR(B),r);

   unsigned i, R=MAX(SIZE.len,B.SIZE.len);  

   if (SIZE.len!=B.SIZE.len) {
      if (  SIZE.len==r) {   appendSingletons(R); rval|=1; } else
      if (B.SIZE.len==r) { B.appendSingletons(R); rval|=2; }
      else {
          wblog(F_L,"ERR %s() rank mismatch (r=%d/%d)",
          FCT,SIZE.len,B.SIZE.len);
      }
      eq=SIZE.isEqual(B.SIZE); 
   }

   if (r<3) {  
      for (i=r; i<R; ++i) { if (  SIZE[i]!=1) { break; }}
      if (i<R || !eq) wblog(FL,"ERR %s() "
         "got OM space %s / %s for rank r=%d",FCT,SSTR_(this),SSTR(B),r);
   }
   else if (!eq) {
      wbvector<size_t> S(SIZE); unsigned q=0;
      for (i=r; i<SIZE.len; ++i) {
          if (S[i]<B.SIZE[i]) { S[i]=B.SIZE[i]; ++q; }}
      if (q) { this->Resize(S); rval|=4; }

      S=B.SIZE; q=0;
      for (i=r; i<SIZE.len; ++i) {
          if (S[i]<SIZE[i]) { S[i]=SIZE[i]; ++q; }}
      if (q) { B.Resize(S); rval|=8; }
   }

   return rval;
};

template<class T>
template<class TI> inline
int wbarray<T>::aMax(T &x, TI &k) const {
   int e=0; size_t n=numel(); x=0; k=0; 
   if (n) { size_t i=0; T a;
      for (; i<n; ++i) { if (x<(a=Wb::abs(data[i]))) { x=a; k=i; }}
   }
   else { k=-1; e=1; }
   return e;
};

template<class T>
template<class TI> inline
T wbarray<T>::aMax(TI &i_, TI &j_, const wbarray *B) const {
   T x=0; i_=j_=0; 
   if (SIZE.len!=2) wblog(FL,
      "ERR %s() invalid usage (s=%s)",FCT,SSTR_(this));

   if (!B) { 
      size_t k; int e=aMax(x,k);
      if (e)
           { wblog(FL, "WRN %s() got empty array (e=%d)",FCT,e); }
      else { i_=k%SIZE[0]; j_=k/SIZE[0]; } 
      return x;
   }

   if (B->SIZE.len!=2) wblog(FL,
      "ERR %s() invalid usage (s=%s)",FCT,SSTR_(B));

   TI na=SIZE[0], nb=B->SIZE[0],
      ma=SIZE[1], mb=B->SIZE[1], n=MIN(na,nb), m=MIN(ma,mb), i,j;

   if (n && m) { T a, b;
      for (i=0; i<n; ++i) {
      for (j=0; j<m; ++j) {
         a=Wb::abs(   data[i+j*na]); 
         b=Wb::abs(B->data[i+j*nb]); if (a<b) { a=b; }
         if (x<a) { x=a; i_=i; j_=j; }
      }}
   }
   else { wblog(FL, "WRN %s() got empty array (%dx%d)",FCT,n,m); }

   return x;
};

template<class T> inline
double wbarray<T>::maxDiff (const wbarray<T> &B) const {
   size_t i,s=numel();
   double dmax=0;

   if (SIZE!=B.SIZE) {
      sprintf_str("%s:%d data SIZE mismatch ([%s, %s])",
         FL,SSTR_(this),SSTR(B));
      return NAN;
   }

   for (i=0; i<s; i++)
   dmax=MAX(dmax, Wb::abs(data[i]-B.data[i]));

   return dmax;
};

template<class T> inline
void wbarray<T>::TimesEl(const wbarray<T> &B, char conj) {

   if (!sameSize(B)) wblog(FL,
      "ERR %s() size mismatch: %s <> %s",FCT,SSTR_(this),SSTR(B));

   Wb::TimesElRange(data,B.data,numel(),conj);
};

template<class T>
wbarray<T>& wbarray<T>::timesEl(
   const wbarray &B, wbarray &C, T bfac, T cfac, char conj) const {

   if (!sameSize(B)) wblog(FL,
      "ERR %s() size mismatch\n%s <> %s",FCT,SSTR_(this),SSTR(B));

   if (!cfac || C.isEmpty()) {
      if (bfac) { C=*this;
         Wb::TimesElRange(C.data,B.data,numel(),conj);
         if (bfac!=T(1)) { C*=bfac; }
      }
      else C.init(SIZE);
   }
   else if (!sameSize(C)) { wblog(FL,
     "ERR %s() size mismatch\n%s <> %s",FCT,SSTR_(this),SSTR(C)); }
   else {
      if (cfac!=T(1)) { C*=cfac; }
      if (bfac) {
         Wb::timesElRange_add(C.data, data,B.data,numel(),bfac,conj);
      }
   }
   return C;
};

template<class T>
wbarray<T>& wbarray<T>::timesEl_OM(
   unsigned r, const wbarray &B, wbarray &C, char conj
 ) const {

   if (!r || r>SIZE.len) wblog(FL,
      "ERR %s() invalid rank (r=%d/%d)",FCT,r,SIZE.len);
   if (!sameSize(B)) wblog(FL,
      "ERR %s() size mismatch %s / %s",FCT,SSTR_(this),SSTR(B));

   unsigned i=0; size_t N=1, M=1;
   const size_t *sa=SIZE.data, *sb=B.SIZE.data;

   bool self=( &C==this || &C==&B || C.data==data || C.data==B.data);

   for (; i<r; ++i) {
      if (sa[i]!=sb[i]) wblog(FL,"ERR %s() "
        "size mismatch %s / %s (i=%d/%d)",FCT,SSTR_(this),SSTR(B),i+1,r);
      N*=sa[i];
   }
   for (; i<SIZE.len; ++i) { M*=sa[i]; }

   if (!self) { 
      C.init_bare(wbvector<size_t>(r,SIZE.data,'r'));
   }

   Wb::timesElRange_OM(C.data, data,B.data,N,M,conj);

   if (self)
        { C.SIZE.len=r; } 

   return C;
};

template<class T> inline
T wbarray<T>::dotProd(const wbarray<T> &B, char conj) const {
   T x=T(0);
   if (!sameSize(B)) wblog(FL,
      "ERR %s() size mismatch: %s <> %s",FCT,SSTR_(this),SSTR(B));

   return Wb::dotProd(data,B.data,numel(),x,conj);
};

template<class T> inline
T wbarray<T>::weightedAvg(const wbarray<T> &B) const {
   size_t i,s=numel(); T w,x=0,n2=0;

   if (!sameSize(B)) wblog(FL,
      "ERR %s() size mismatch: %s <> %s",FCT,SSTR_(this),SSTR(B));

   for (i=0; i<s; i++) { w=B.data[i]; x+=(data[i]*w); n2+=w*Wb::CONJ(w); }

   if (n2<=0) wblog(FL,
      "ERR weighted avg with total weight %.4g",n2);
   return x/std::sqrt(n2);
};

template<class T>
wbarray<T>& wbarray<T>::tensorProd(
   const wbarray<T> &B_, wbarray<T> &C,
   const char aflag0, const char bflag0,
   const char kflag 
) const {

   if (this==&C || &B_==&C) {
      wbarray<T> X; tensorProd(B_,X,aflag0,bflag0,kflag);
      return X.save2(C);
   }

   wbarray<T> A,B;
   opFlags<T> aflag(aflag0), bflag(bflag0);
   size_t ra=SIZE.len, rb=B_.SIZE.len;

   aflag.applyOrRef(FL,*this,A); 
   bflag.applyOrRef(FL, B_,  B);

   if (ra!=rb) {
      if (ra>2 || rb>2) wblog(FL,
         "ERR %s() rank mismatch (r=%d/%d)",FCT,ra,rb);
      if (!ra || !rb) { return C.init(); }

      if (ra==1) { A.ExpandDiagonal(); ra=  SIZE.len; } else
      if (rb==1) { B.ExpandDiagonal(); rb=B.SIZE.len; }

      if (ra!=rb) wblog(FL, 
         "ERR %s() rank mismatch (%d/%d)",FCT,ra,rb);
   }

   const size_t r=A.SIZE.len, l=r-1, *sa=A.SIZE.data, *sb=B.SIZE.data;
   WBINDEX Ia(r), Ib(r), Ja(r), Jb(r), S;

   size_t i=0,k=0,ia,ib,N;

   if (kflag) { 
      S.init(r);
      for (; i<r; ++i) { S[i]=sa[i]*sb[i]; }
   }
   else {
      S.init(2*r);
      for (; i<r; ++i) { S[k++]=sa[i]; S[k++]=sb[i]; }
   }

   C.init(S); N=S.prod(); Ja[l]=Jb[l]=0;

   for (k=l,i=0; i<N; ++i) {
       for (--k; k<r; --k) { 
           Ja[k] = ( Ja[k+1] + Ia[k+1] ) * sa[k];
           Jb[k] = ( Jb[k+1] + Ib[k+1] ) * sb[k];
       }

       ia = Ja[0] + Ia[0];
       ib = Jb[0] + Ib[0];
       C.data[i] = A.data[ia] * B.data[ib];

       k=0;
       if ( (++Ia[k]) >= sa[k] ) { Ia[k]=0; ++Ib[k]; }
       while (Ib[k]>=sb[k] && k<l) {
          Ib[k]=Ia[k]=0; ++k; 
          if ( (++Ia[k]) >= sa[k] ) { Ia[k]=0; ++Ib[k]; }
       }
   }

   return C;
};

template<class T>
wbarray<T>& wbarray<T>::Cat(
   unsigned dim, 
   const wbvector< wbarray<T> > &aa
){
   wbvector< const wbarray<T>* > ap(aa.len);
   for (unsigned i=0; i<aa.len; i++) ap[i]=&aa[i];
   Cat(dim,ap); return *this;
};

template<class T>
wbarray<T>& wbarray<T>::Cat(
   unsigned dim, 
   const wbvector< wbarray<T>* > &ap
){
   if (ap.len==0) { init(); return *this; }
   if (ap.len==1) { *this=(*ap[0]); return *this; }

   size_t i,j=0,k,n, e=0, N, *s2=0;
   const size_t r=ap[0]->SIZE.len, l=r-1, *sref=ap[0]->SIZE.data;
   const WBINDEX &S0=ap[0]->SIZE;
   WBINDEX S(S0), I(r), I0(ap.len), idx;

   if (dim<1) wblog(FL,"ERR invalid index (not 1-based; %d/%d)",dim,r);
   dim--; 

   if (dim<r) N=sref[dim];

   for (i=1; i<ap.len; i++) {
      if (ap[i]==NULL || ap[i]->data==NULL) continue;

      if (r != ap[i]->SIZE.len) e++; else {
         s2=ap[i]->SIZE.data;
         for (j=0; j<r; j++) if (j!=dim && s2[j]!=sref[j]) e++;
      }
      if (e) wblog(FL,"ERR size mismatch %d/%d: %s <> %s (%d)",
         i+1,ap.len, SSTR_(ap[0]),SSTR_(ap[i]),dim+1);

      I0[i]=N; if (dim<r) N+=s2[dim];
   }

   if (dim>=r) { T *dd;
      S.Resize(dim+1); S[dim]=ap.len; init(S); dd=data; n=ap[0]->numel();
      for (i=0; i<ap.len; i++, dd+=n) {
         MEM_CPY<T>(dd,n,ap[i]->data);
      }
      return *this;
   }

   idx.init(N);
   for (k=i=0; i<ap.len; i++) {
      if (ap[i]==NULL || ap[i]->data==NULL) continue;
      for (n=ap[i]->SIZE.data[dim], j=0; j<n; j++, k++) idx[k]=i;
   }

   S[dim]=N; init(S); n=numel();

   for (i=0; i<n; i++) {
       widx_t q=idx[I[dim]];
       const size_t *s=ap[q]->SIZE.data;

       for (j=0,k=l; k<r; k--) {
          if (j) j*=s[k];
          j += (k!=dim ? I[k] : (I[k]-I0[q]));
       }

       data[i]=ap[q]->data[j];

       k=0; I[0]++; 
       while (I[k]>=SIZE[k] && k<l) { I[k]=0; ++I[++k]; }
   }

   return *this;
};

template<class T>
wbarray<T>& wbarray<T>::BlockCat(
   const wbarray< const wbarray<T>* > &ap,
   WBINDEX *D1_, WBINDEX *D2_
){
   size_t i,j,m,n,i0,j0,d1,d2,M,N;
   WBINDEX D1,D2;

   if (ap.isEmpty()) {
      init(); if (D1_) D1_->init(); if (D2_) D2_->init();
      return *this;
   }

   if (ap.rank()!=2) wblog(FL,
      "ERR %s() invalid input matrix (%s)",FCT,SSTR(ap));

   m=ap.SIZE[0]; n=ap.SIZE[1]; D1.init(m); D2.init(n);
   for (i=0; i<m; i++)
   for (j=0; j<n; j++) if (ap(i,j) && !ap(i,j)->isEmpty()) {
      const wbarray<T> &a=(*ap(i,j));
      if (a.rank()!=2) wblog(FL,
         "ERR %s() invalid block rank (%d,%d): %s",FCT,i+1,j+1,SSTR(a));

      if (D1[i]) {
         if (D1[i]!=a.SIZE[0]) wblog(FL,
         "ERR %s() block-size mismatch (%d,%d: %s <> %dx..",
          FCT, i+1, j+1, SSTR(a), D1[i]);
      }
      else D1[i]=a.SIZE[0];

      if (D2[j]) {
         if (D2[j]!=a.SIZE[0]) wblog(FL,
         "ERR %s() block-size mismatch (%d,%d: %s <> ..x%d",
          FCT, i+1, j+1, SSTR(a), D2[j]);
      }
      else D2[j]=a.SIZE[0];
   }

   if (D1_) D1_->init();
   if (D2_) D2_->init();

   M=D1.sum(); N=D2.sum();
   init(M,N); if (isEmpty()) return *this;

   for (j0=j=0; j<n; j++, j0+=d2) { d2=D2[j];
   for (i0=i=0; i<m; i++, i0+=d1) { d1=D1[i];
       ap(i,j)->copyStride(ref(i0,j0),M); 
   }}

   return *this;
};

template<class T>
wbarray<T>& wbarray<T>::blockTrace(
   size_t D1, 
   wbarray<T> &a, size_t D2 
) const {

   if (rank()!=2) wblog(FL,
      "ERR %s() only applies to rank-2 tensors (%d)",FCT,rank());

   size_t i, dim1=SIZE[0], dim2=SIZE[1], m=(D1 ? dim1/D1 : 0);
   if (int(D2)<0) { D2=D1; }

   if (!D1 || !D2 || dim1%D1 || dim2%D2 || (dim2/D2)!=m) wblog(FL,
      "ERR %s() invalid block specs (%d/%d; %d/%d; %d)",
       D1,dim1, D2,dim2, m);

   a.init(D1,D2);
   for (i=0; i<m; ++i) { addBlock(i*D1, i*D2, D1, D2, a); }

   return a;
};

template<class T>
wbarray<T>& wbarray<T>::BlockDiag(const wbvector< wbarray<T> > &D) {

   size_t i,d1=0,d2=0, i0=0, j0=0, D1=0, D2=0;

   for (i=0; i<D.len; i++) {
      if (!D[i].isMatrix()) wblog(FL,
         "ERR %s() matrices expected (%s)",FCT,SSTR(D[i]));
      D1+=D[i].SIZE[0];
      D2+=D[i].SIZE[1];
   }

   init(D1,D2);
   for (i=0; i<D.len; i++, i0+=d1,j0+=d2) {
      d1=D[i].SIZE[0];
      d2=D[i].SIZE[1];
      Wb::cpyStride(ref(i0,j0), D[i].data,d1,d2,D1);
   }

   return *this;
};

template<class T>
wbarray<T>& wbarray<T>::addBlock(
   size_t i, size_t j, 
   size_t n, size_t m, 
   wbarray<T> &a
) const {

   if (rank()!=2) wblog(FL,
      "ERR %s() only applies to rank-2 tensors (%s)",FCT,SSTR_(this));
   size_t dim1=SIZE[0], dim2=SIZE[1];

   if (i+n>dim1 || j+m>dim2) wblog(FL,
      "ERR %s() index out of bounds (%d:%d/%d; %d:%d/%d)",
       i+1,i+n,dim1,j+1,j+m,dim2);

   if (a.isEmpty()) a.init(n,m);
   Wb::cpyStride(a.data, ref(i,j), n,m, -1,dim1,'+'); 

   return a;
};

template<class T>
T wbarray<T>::norm2block(
   size_t i0, size_t j0, 
   size_t n,  size_t m   
) const {

   T x2=0; 

   if (rank()!=2) wblog(FL,
      "ERR %s() only applies to rank-2 tensors (%s)",FCT,SSTR_(this));
   size_t dim1=SIZE[0], dim2=SIZE[1];

   if (i0+n>dim1 || j0+m>dim2) wblog(FL,
      "ERR %s() index out of bounds (%d:%d/%d; %d:%d/%d)",
       i0+1,i0+n,dim1, j0+1,j0+m,dim2);

   if (m && n) {
      size_t i,j; const T* d = data + i0 + j0*dim1; 
      for (j=0; j<m; ++j, d+=dim1) {
      for (i=0; i<n; ++i) { x2+=Wb::norm2(d[i]); }}
   }

   return x2;
};

template<class T>
wbarray<T>& wbarray<T>::getBlock(
   size_t i, size_t j, 
   size_t n, size_t m, 
   wbarray<T> &a           
) const {

   if (rank()!=2) wblog(FL,
   "ERR %s() only applies to rank-2 tensors (%d)",FCT,rank());

   size_t dim1=SIZE[0], dim2=SIZE[1];

   if (i+n>dim1 || j+m>dim2) wblog(FL,
      "ERR %s() index out of bounds (%d..%d/%d; %d..%d/%d)",
       i+1,i+n,dim1,j+1,j+m,dim2);

   a.init(n,m); Wb::cpyStride(a.data, ref(i,j),n,m,-1,dim1);
   return a;
};

template<class T>
void wbarray<T>::copyStride(T* dd, size_t stride, T afac) const {

   if (SIZE.len!=2) wblog(FL,"ERR %s() "
      "only intended for rank-2 tensors (r=%d)",FCT,rank());

   size_t i, j=0, dim1=SIZE[0], dim2=SIZE[1];
   T *d0=data;

   if (afac==T(1)) { 
      for (; j<dim2; ++j, dd+=stride, d0+=dim1) {
         MEM_CPY<T>(dd,dim1,d0);
      }
   }
   else if (afac==T(-1)) {
      for (; j<dim2; ++j, dd+=stride, d0+=dim1) {
      for (i=0; i<dim1; ++i) { dd[i]= -d0[i]; }}
   }
   else if (afac) {
      for (; j<dim2; ++j, dd+=stride, d0+=dim1) {
      for (i=0; i<dim1; ++i) { dd[i]=afac*d0[i]; }}
   }
   else {
      for (; j<dim2; ++j, dd+=stride) {
      for (i=0; i<dim1; ++i) { dd[i]=0; }}
   }
};

template<class T>
void wbarray<T>::addStride(T* dd, size_t stride, const T afac) const {

   if (SIZE.len==2) wblog(FL,"ERR %s() "
      "only intended for rank-2 tensors (r=%d)",FCT,rank());

   size_t i=0, j=0, dim1=SIZE[0], dim2=SIZE[1];
   T *a=data; 

   if (afac==T(1)) {
      for (; j<dim2; ++j, a+=dim1, dd+=stride) {
      for (i=0; i<dim1; ++i) { dd[i] += a[i]; }} 
   }
   else if (afac==T(-1)) {
      for (; j<dim2; ++j, a+=dim1, dd+=stride) {
      for (i=0; i<dim1; ++i) { dd[i] -= a[i]; }} 
   }
   else if (afac) {
      for (; j<dim2; ++j, a+=dim1, dd+=stride) {
      for (i=0; i<dim1; ++i) { dd[i]+=(afac*a[i]); }} 
   }
};

template<class T> inline
wbarray<T>& wbarray<T>::Reshape(const wbvector<size_t> &S, bool lflag) {

   if (numel()!=S.prod(0)) wblog(FL,
      "ERR %s() length mismatch (%s => %s)",FCT, SSTR_(this), SSTR(S));

   if (!lflag) {
      unsigned i=0, j=0, e=0;
      const size_t *s1=SIZE.data, *s2=S.data; size_t sb;

      while (!e) {
         for (; i<SIZE.len; ++i) { if (SIZE[i]!=1) { break; }}
         for (; j<S.len;    ++j) { if (S[j]   !=1) { break; }}
         if (i>=SIZE.len || j>=S.len) { break; }

         if (s1[i]==s2[j]) { ++i; ++j; }
         else if (s1[i]<s2[j]) {
            for (sb=s1[i++]; i<SIZE.len; ++i) {
               if (( sb*=s1[i] ) >= s2[j]) {
                  if (sb==s2[j]) { ++i; ++j; } else { ++e; }
                  break;
               }
            }
         }
         else if (s1[i]>s2[j]) {
            for (sb=s2[j++]; j<S.len; ++j) {
               if (( sb*=s2[j] ) >= s1[i]) {
                  if (sb==s1[i]) { ++i; ++j; } else { ++e; }
                  break;
               }
            }
         }
      }
      if (e || i<SIZE.len || j<S.len) wblog(FL,"WRN %s() "
      "got non-trivial reshape %s -> %s",FCT,SSTR_(this),SSTR(S));
   }

   SIZE=S; return *this;
};

template<class T>
void wbarray<T>::GroupIndizes(size_t K) {
    size_t i,n;

    if (!K || SIZE.len%K) wblog(FL,
    "ERR %s() invalid block specs: mod(%d,%d)",FCT,SIZE.len,K);

    n=SIZE.len/K;

    for (i=0; i<n; ++i)
    SIZE[i]=Wb::prodRange(SIZE.data+i*K, K);

    SIZE.Resize(n);
}

template<class T> inline
void wbarray<T>::groupIndizes_P(
   const WBINDEX &I, int pos, wbperm &P,
   size_t *s1, size_t *s2
) const {

   size_t i,j,j1=0,j2=0,l, m=I.len, n=SIZE.len;
   wperm_t *p; widx_t *idx;
   char mark[n ? n : 1]; 

   if (m>n) wblog(FL,
      "ERR %s() invalid group index (len=%d/%d)",FCT,m,n);
   if (!n) {
      if (s1) { *s1=0; }
      if (s2) { *s2=0; }; return;
   }
   P.init(n); p=P.data; idx=I.data;

   if (pos==1) { j1=0;   j2=m; } else
   if (pos==2) { j1=n-m; j2=0; }
   else wblog(FL,"ERR %s() invalid pos=%d",FCT,pos);

   memset(mark,0,n*sizeof(char));

   for (j=j1, i=0; i<m; ++i, ++j) { l=p[j]=idx[i];
      if (l>=n) wblog(FL,"ERR %s() index not unique (%d/%d)",FCT,l,n);
      if (++mark[l]>1) wblog(FL,
         "ERR %s() index not unique (%d/%d)",FCT,l,n);
   }
   for (j=j2, i=0; i<n; ++i) { if (!mark[i]) { p[j++]=i; }}

   if (s1 || s2) {
      for (j1=j2=1, i=0; i<n; ++i) {
         if (mark[i])
              { j1*=SIZE[i]; }
         else { j2*=SIZE[i]; }
      }
      if (s1) { *s1 = (pos==1 ? j1 : j2); }
      if (s2) { *s2 = (pos==1 ? j2 : j1); }
   }
};

template<class T>
wbarray<T>& wbarray<T>::transpose(const char *F, int L, wbarray<T> &A) const {

    unsigned r=rank(); wbperm P;

    if (r%2) wblog(F_L,"ERR %s() "
       "applies to even-rank arrays only! (%s)",FCT,SSTR_(this));

    P.initTranspose(r);
    permute(A,P);

    return A;
};

template<class T>
wbarray<T>& wbarray<T>::flipLR(const char *F, int L, wbarray &B) const {

   if (rank()!=2) wblog(F_L,
      "ERR %s() requires matrix (%s)",FCT,SSTR_(this));
   if (&B==this || B.data==data) wblog(F_L,
      "ERR %s() requires different objects!",FCT);
   B.init(SIZE);

   if (SIZE[1]) {
      widx_t i=0, j=0, dim1=SIZE[0], dim2=SIZE[1];
      const T *a=data; T *b=B.data+dim1*(dim2-1);
      for (; j<dim2; ++j, a+=dim1, b-=dim1) { 
         for (i=0; i<dim1; ++i) b[i]=a[i];
      }
   }

   return B;
};

template<class T>
wbarray<T>& wbarray<T>::FlipLR(const char *F, int L) {

   if (rank()!=2) wblog(F_L,
      "ERR %s() requires matrix (%s)",FCT,SSTR_(this));

   widx_t i=0, j=0, dim1=SIZE[0], dim2=SIZE[1], d2=dim2/2;

   if (dim2>1) {
      T x, *a=data, *b=data+dim1*(dim2-1);
      for (; j<d2; ++j, a+=dim1, b-=dim1) { 
         for (i=0; i<dim1; ++i) {
            x=a[i]; a[i]=b[i]; b[i]=x; 
         }
      }
   }

   return *this;
};

template<class T>
wbarray<T>& wbarray<T>::flipUD(const char *F, int L, wbarray &B) const {

   if (rank()!=2) wblog(F_L,
      "ERR %s() requires matrix (%s)",FCT,SSTR_(this));
   if (&B==this || B.data==data) wblog(F_L,
      "ERR %s() requires different objects!",FCT);
   B.init(SIZE);

   if (SIZE[0]) {
      widx_t i=0, j=0, dim1=SIZE[0], dim2=SIZE[1], l=dim1-1;
      const T *a=data; T *b=B.data;
      for (; j<dim2; ++j, a+=dim1, b+=dim1) { 
         for (i=0; i<dim1; ++i) b[i]=a[l-i];
      }
   }

   return B;
};

template<class T>
wbarray<T>& wbarray<T>::FlipUD(const char *F, int L) {

   if (rank()!=2) wblog(F_L,
      "ERR %s() requires matrix (%s)",FCT,SSTR_(this));

   widx_t i=0, j=0, dim1=SIZE[0], dim2=SIZE[1], d2=dim1/2, l=dim1-1;

   if (dim1>1) {
      T x, *a=data;
      for (; j<dim2; ++j, a+=dim1) { 
         for (i=0; i<d2; ++i) {
            x=a[i]; a[i]=a[l-i]; a[l-i]=x; 
         }
      }
   }

   return *this;
};

template<class T>
wbarray<T>& wbarray<T>::MatPermute(const wbperm &P, char iflag){

   if (isIdentityPerm(P)) return *this;
   if (!isSMatrix()) wblog(FL,
      "ERR %s() expects square matrix (%s)",FCT,SSTR_(this));

   wbarray<T> X(*this);
   const wperm_t *p=P.data;
   size_t i,j, dim1=SIZE[0], dim2=SIZE[1];

   if (!iflag) {
      for (i=0; i<dim1; i++) 
      for (j=0; j<dim2; j++) { data[i+dim1*j]=X.data[p[i]+dim1*p[j]]; }
   }
   else {
      for (i=0; i<dim1; i++) 
      for (j=0; j<dim2; j++) { data[p[i]+dim1*p[j]]=X.data[i+dim1*j]; }
   }

   return *this;
};

template<class T>
wbarray<T>& wbarray<T>::Permute(
   const wbperm &P, char iflag, char rcpy) {

   if (!P) { return *this; }   

   if (P.len>SIZE.len) { appendSingletons(P.len); } 

   if (!requiresDataPerm(P)) { 
      SIZE.Permute(P,iflag);   
      return *this;
   }

   if (isRef()) { 
      if (rcpy=='c' || rcpy=='C' || rcpy=='i' || rcpy=='I') {
         wbarray<T> X;
         return permute(X,P,iflag).save2(*this);
      }
      else if (rcpy!='r' && rcpy!='R') {
         wblog(PFL,"ERR %s() got ref (rcpy=%s)",FCT,cSTR(rcpy));
      }
   }

   wbarray<T> X(*this);
   X.permute(*this,P,iflag);

   return *this;
};

template<class T>
wbarray<T>& wbarray<T>::Permute(const char* sidx, char iflag, char rcpy) {
    wbperm P(sidx); 
    if (P.isEmpty()) wblog(FL,"WRN got empty permutation '%s'", sidx);
    return Permute(P,iflag,rcpy);
};

template<class T>
wbarray<T>& wbarray<T>::permute(
    wbarray<T> &A,  
    const char* sidx,
    char iflag      
) const {

    wbperm P(sidx); 

    if (P.isEmpty()) {
       wblog (FL,"WRN got empty permutation '%s'", sidx);
       return A.init(*this);
    }
    return permute(A,P,iflag);
};

template<class T>
wbarray<T>& wbarray<T>::permute(   
    wbarray<T> &A,     
    const wbperm &P0,
    char iflag         
) const {

    if (!P0) { return A=(*this); } 
    if (&A==this) { 
       wbarray<T> X;
       return permute(X,P0,iflag).save2(A);
    }

    wbperm P;

#ifdef WB_CLOCK
   Wb::Clock clk("arr:permute",0); 
#endif

    if (iflag!=0 && iflag !='I') wblog(FL, 
       "ERR %s() invalid iflag '%c'<%d>",FCT,iflag,iflag);

    P.init(P0,iflag,SIZE.len); 

    if (!requiresDataPerm(P)) {
        A=(*this); A.SIZE.Select(P); 
        return A;
    }

#ifdef QS_USING_OMP
   int np=MAX(OMP_NUM_THREADS,QSP_NUM_THREADS);
#else
   int np=0;
#endif

    A.init_bare(SIZE); 
    SIZE.permute(A.SIZE,P);

    wbarray_permute__(A.data,*this,P,np);

    return A;
};

template<class T>
wbarray<T>& wbarray<T>::select0(
    const WBPERM &P,
    unsigned dim, 
    wbarray<T> &A 
) const {

    if (&A==this) {
       wbarray<T> X(A); X.select0(P,dim,A);
       return A;
    }

    const unsigned r=SIZE.len, l=r-1; size_t i,j,k,d;
    WBINDEX S,I;

    if (dim>=SIZE.len) wblog(FL,
       "ERR %s() dimension out of bounds (%d/%d)",FCT,dim+1,SIZE.len);

    S=SIZE; d=SIZE[dim]; S[dim]=P.len;

    for (i=0; i<P.len; i++) if (P[i]>=d) wblog(FL,
        "ERR %s() index out of bounds (%d/%d)",FCT,P[i],d);

    A.init(S); I.init(SIZE.len);
    widx_t *ip=I.data, *s=S.data; const wperm_t *p=P.data;

    for (d=A.numel(), i=0; i<d; i++) {
        for (j=0, k=l; k<r; k--) {
           if (j) j*=SIZE[k];
           j+=(k!=dim ? ip[k] : p[ip[k]]);
        }

        A.data[i]=data[j];

        k=0; ip[0]++; 
        while (ip[k]>=s[k] && k<l) { ip[k]=0; ++ip[++k]; }
    }

    return A;
}

template<class T>
wbarray<T>& wbarray<T>::selectSqueeze(
    wbarray<T> &A, unsigned p, 
    unsigned dim 
) const {

    if (&A==this) {
       wbarray<T> X(A); X.selectSqueeze(A,p,dim);
       return A;
    }

    const unsigned r=SIZE.len, l=r-1; size_t i,j,k,s;
    WBINDEX S,I;

    if (dim>=SIZE.len) wblog(FL,
    "ERR %s() dimension out of bounds (%d/%d)",FCT,dim+1,SIZE.len);

    if (p>=SIZE[dim]) wblog(FL,
    "ERR %s() index out of bounds (%d/%d)",FCT,p,SIZE[dim]);

    if (SIZE.len==1) {
       A.init(1); A.data[0]=data[p];
       return A;
    }

    S=SIZE; S[dim]=1; A.init(S); s=A.numel(); I.init(SIZE.len);

    for (i=0; i<s; i++) {
        for (j=0, k=l; k<r; k--) {
           if (j) j*=SIZE[k];
           j+=(k!=dim ? I[k] : p);
        }

        A.data[i]=data[j];

        k=0; I[0]++; 
        while (I[k]>=S[k] && k<l) { I[k]=0; ++I[++k]; }
    }

    for (i=dim+1; i<S.len; i++) A.SIZE[i-1]=A.SIZE[i];
    A.SIZE.len--; A.SIZE[A.SIZE.len]=0;

    return A;
}

template<class T>
void wbarray<T>::setBlock(
    const wbvector< WBINDEX > &DB, 
    const wbindex &IB,   
    const wbarray<T> &A
){
    size_t i,j,ik,k,s, r=SIZE.len, l=r-1;
    wbindex I,I1(r),I2(r); 

    if (DB.len!=r || IB.len!=r) wblog(FL,
       "ERR %s() size mismatch (r: %d %d / %d)",FCT,DB.len,IB.len,r);

    if (A.SIZE.len!=r) {
       for (i=r; i<A.SIZE.len; ++i) { if (A.SIZE[i]!=1) break; }
       if (i!=A.SIZE.len) wblog(FL,
       "ERR %s() size mismatch (r=%d / %d)",FCT,A.SIZE.len,r);
    }

    for (k=0; k<r; ++k) {
        const WBINDEX &Dk=DB[k]; ik=IB[k];

        if (ik>=Dk.len) wblog(FL,
           "ERR %s() block index out of bounds (%d: %d/%d)",
            FCT,k+1, ik+1, Dk.len);
        if (Dk[ik] != A.SIZE[k]) wblog(FL,
           "ERR %s() block size mismatch (%d: %d/%d)",
            FCT,k+1, Dk[ik], A.SIZE[k]);

        I1[k]=Wb::addRange(Dk.data,ik);
        I2[k]=I1[k]+Dk[ik];

        if (I2[k]>SIZE[k]) wblog(FL,
           "ERR %s() index out of bounds (%d: %d/%d)",
            FCT,k+1, I2[k], SIZE[k]);
    }

    I=I1; s=A.numel();

    for (i=0; i<s; ++i) {
        for (j=I[l],k=l-1; k<r; k--) j = j*SIZE[k] + I[k];

        data[j]=A.data[i];

        k=0; I[0]++;
        while (I[k]>=I2[k] && k<l) { I[k]=I1[k]; ++I[++k]; }
    }
}

template<class T>
void wbarray<T>::addBlock(
    const WBINDEX &I0, 
    const wbarray<T> &A,
    const char dflag 
){
    size_t i,j,k, s=A.numel(), r=SIZE.len, l=r-1;
    unsigned ra=(dflag ? r/2 : r), la=ra-1;
    const widx_t *S=A.SIZE.data;
    WBINDEX I(ra);

    if (I0.len!=r || A.SIZE.len!=ra || (dflag && r%2)) wblog(FL,
       "ERR %s() size mismatch (%d/%d, %s; %c<%d>)",
        FCT,I0.len, A.SIZE.len, SSTR_(this), dflag
    );
    for (i=0; i<r; i++) if (I0[i]+S[i>=ra ? i-ra : i]>SIZE[i])
         wblog(FL,"ERR %s() index out of bounds (%d: %s - %s %s)",
         FCT,i+1,I0.toStr().data,SSTR(A),SSTR_(this)
    );

    for (i=0; i<s; i++) {
       for (j=0,k=l; k<r; k--) { if (j) j*=SIZE[k];
          j+=(I0[k]+I[k<ra ? k : k-ra]);
       }

       data[j] += A.data[i];

       k=0; I[0]++;
       while (I[k]>=S[k] && k<la) { I[k]=0; ++I[++k]; }
    }
};

template<class T>
T wbarray<T>::min() const {
   size_t i,n=numel(); T x;

   if (n==0) {
       wblog(FLINE,"WRN min() of empty wbarray !?");
       memset(&x,0,sizeof(T)); return x;
   }

   for (x=data[0], i=1; i<n; i++)
   if (x>data[i]) x=data[i];

   return x;
}

template<class T>
T wbarray<T>::max() const {
   size_t i,n=numel(); T x;

   if (n==0) {
       wblog(FLINE,"WRN max() of empty wbarray !?");
       memset(&x,0,sizeof(T)); return x;
   }

   for (x=data[0], i=1; i<n; i++)
   if (x<data[i]) x=data[i];

   return x;
}

template<class T>
T wbarray<T>::aMin(char zflag, size_t *k_) const {

   size_t i=0, k=-1,n=numel(); T a,x=0;

   if (!n) {
      if (k_) { (*k_)=-1; wblog(FL,"WRN %s() got empty array",FCT); }
      return 0;
   }

   if (zflag==0) { 
      x=Wb::abs(data[i++]); k=0; if (x!=0) {
         for (; i<n; ++i) { a=Wb::abs(data[i]); if (x>a) {
            x=a; k=i; if (x==0) break;
         }}
      }
   }
   else {
      for (; i<n; ++i) {
         a=Wb::abs(data[i]); if (a) { x=a; k=i; break; }
      }
      for (; i<n; ++i) { a=Wb::abs(data[i]);
         if (a && x>a) { x=a; k=i; }
      }
   }
   if (k_) (*k_)=k;
   return x;
};

template<class T>
wbstring wbarray<T>::toStr() const {
   wbstring s; 
   unsigned l=2; size_t n=numel();
   if (n==1) { s.init(24+4*SIZE.len);
      l=snprintf(s.data,s.len, SIZE.len==2? "%s (%s)" : "%s [%s]",
      STR(wbvector<T>(numel(),data,'r')), SSTR_(this));
   }
   else if (n || SIZE.len) { s.init(5 + 10*(numel()+SIZE.len));
      l=snprintf(s.data,s.len,"[%s](%s)",
      STR(wbvector<T>(numel(),data,'r')), SSTR_(this));
   }
   else { s="[]"; }

   if (l>=s.len) wblog(FL,
      "WRN %s() string out of bounds (%d/%d)",FCT,l,s.len);
   return s;
};

template<class T>
wbstring Wb::sizeStrM( 
   unsigned r, const T* sd,
   unsigned m, const T* sm, const char *sep, const char *sepM
) {
   unsigned i=0, j, l=0, n=0, ndims=r+m;
   char *s; size_t x;
   wbstring sout;

   if (r>1 || m>1)
          { n+=(ndims-(r && m ? 2 : 1))*(sep ? strlen(sep) : 1); }
   if (m) { n+=(sepM ? strlen(sepM) : 1); }

   for (; i<ndims; ++i) {
      x=(i<r ? sd[i] : sm[i-r]); if (x<0) { x=-x; ++n; } 
      for (j=0; j<64; ++j) { if (!(x>>1)) break; }
      n+=(1+ceil(3*(double(i)/10))); 
   }; if (n<8) n=8;

   sout.init(n); s=sout.data;

   for (i=0; i<r && l<n; ++i) {
      l+=snprintf(s+l,n-l,"%s%ld",i? (sep? sep:" "):"", long(sd[i]));
   }
   for (i=0; i<m && l<n; ++i) {
      l+=snprintf(s+l,n-l,"%s%ld", i? (sep? sep:" ") : (sepM? sepM:"|"),
      long(sm[i]));
   }

   if (l>=n) wblog(FL,
      "WRN %s() string out of bounds (%d/%d)\n'%s'",FCT,l,n,s);
   return sout;
};

template<class T>
void wbarray<T>::print_ref (const char *sb, const char *nl) const {

   if (isRef() || mtype || (sptr && sptr->mtype)) {
      int l=0, n=32; char str[n]; str[0]=0;

      if (sptr) {
         l+=snprintf(str,n,"%s",STR_(sptr));
         if (mtype!=sptr->mtype && l<n)
         l+=snprintf(str+l,n-l," / %s !?",Wb::MTYPE_STR[mtype]);
      }
      else { l+=snprintf(str,n,"%s",Wb::MTYPE_STR[mtype]); }

      printf(" %s%s%s%s", sb?sb:"", str, sb?sb:"", nl?nl:"");

      if (l>=n) { printf("\n");
         wblog(FL,"WRN %s() string out of bounds (%d/%d)",FCT,l,n);
      }
   }
   else if ((nl && nl[0])) { printf("%s",nl); }

   if (sptr && data!=sptr->data) {
      if (!nl && !nl[0]) { PRINTF("\n"); }
      wblog(FL, "WRN %s() data=%p / %p !?",FCT,data,sptr->data);
   }
};

template<class T>
void wbarray<T>::info(const char *F, int L,
   const char* istr, unsigned k, 
   unsigned nlt, unsigned nlb 
 ) const {

   if (nlt) {
      printf("\n"); --nlt; if (nlt>1) {
      printf("\n"); --nlt; }

      if (nlt) { --nlt;
         printf("%-20s %-8s %-8s %s\n",
         "", 
         "name","size","data type");
      }
      if (nlt>4) nlt=4; 
      for (unsigned i=0; i<nlt; ++i) { printf("\n"); }
   }

   char s[32];
   if (int(k)>=0)
        { snprintf(s,32,"%s[%d]",istr?istr:"arr",k+1); }
   else { snprintf(s,32,"%s",istr?istr:""); }

   #pragma omp critical (__ensure_sequential_wblog__)
   {  printf("%-20s %-8s %-8s %s",shortFL(F_L),s,SSTR_(this),TSTR(T));
      print_ref();
   }

   if (nlb) {
      if (nlb>10) nlb=10; 
      for (unsigned i=0; i<nlb; ++i) { printf("\n"); }
   }
};

template<class T>
void wbarray<T>::info(const char* istr) const {

   unsigned l=0, n=16; char s[n];
   size_t b=numel()*sizeof(T); 

   if (b<(1<<10)) l=snprintf(s,n,"%ld bytes ",b); else
   if (b<(1<<20)) l=snprintf(s,n,"%.3g kB",b/double(1<<10)); else
                  l=snprintf(s,n,"%.3g MB",b/double(1<<20));
   if (l>=n) wblog(FL,"ERR %s() string out of bounds (%d/%d)",FCT,l,n);

   #pragma omp critical (__ensure_sequential_wblog__)
   {  printf("  %-12s %-10s %12s  @ %12p %s",
         istr? istr:"", SSTR_(this), s, data, TSTR(T));
      print_ref();
   }
};

template<class T>
void wbarray<T>::print(
   const char *F, int L, const char* istr, const char* fmt0
 ) const {

   if (F) {
      wbstring tstr;
      if (typeid(T)==typeid(double)) tstr="double";
      else tstr.init(((char*)(typeid(T).name())));

      wblog(F,L,"%s = [%s] %s array \\",
         istr[0]? istr:"ans", SSTR_(this), tstr.data);
      print_ref();
   }
   else info(istr);

   printdata(istr, fmt0);
};

template<class T>
int wbarray<T>::printdata( 
   const char *istr,
   const char *fmt0
) const {

   unsigned i,j,k, r=SIZE.len, l=r-1, s=numel();
   wbstring fmt;

   if (!fmt0[0]) {
       if (typeid(T)==typeid(double) || typeid(T)==typeid(float))
            fmt=" %12.5g";
       else fmt.init2Fmt((T)0, 6);
   } else fmt=fmt0;

   if (s>1) {
      if (SIZE.len<2) {
         printf("\n");
         for (i=0; i<s; i++) printf(fmt.data,data[i]);
         if (s) printf("\n\n");
      }
      else {
         WBINDEX I(r);

         wbarray<T> A; wbperm P(r); P[0]=1; P[1]=0;
         permute(A,P);

         for (k=2, i=0; i<s; i++) {
            if (k) { if (i) printf("\n");
               if (k>1) {
                  if (r>2) {
                     printf("\n  %s(:,:",istr);
                     for (j=2; j<r; j++) printf(",%ld",I[j]+1);
                     printf(") = [\n\n");
                  }
                  else printf("%s = [\n",istr);
               }
            }

            printf(fmt.data,A.data[i]);

            k=0; I[0]++;  
            while (I[k]>=A.SIZE[k] && k<l) { I[k]=0; ++I[++k]; }
         }
         printf("\n];\n");
      }
   }
   else if (s==1) {
      printf("  %s = ",istr);
      printf(fmt.data,data[0]); printf("\n");
   }
   else printf("\n");

   return (int)s;
};

template<>
int wbarray<wbcomplex>::printdata( 
   const char *istr, const char *fmt0
) const {

   unsigned i,j,k, r=SIZE.len, l=r-1, s=numel();
   wbstring fmt;

   fmt=((fmt0 && fmt0[0]) ? fmt0 : " %8.4g%+8.4gi");

   if (s>1) {
      if (SIZE.len<2) {
         printf("\n");
         for (i=0; i<s; i++) printf(fmt.data,data[i].r,data[i].i);
         if (s) printf("\n\n");
      }
      else {
         WBINDEX I(r);

         wbarray<wbcomplex> A; wbperm P(r); P[0]=1; P[1]=0;
         permute(A,P);

         for (k=2, i=0; i<s; i++) {
            if (k) { if (i) printf("\n");
               if (k>1) {
                  if (r>2) {
                     printf("\n  %s(:,:",istr);
                     for (j=2; j<r; j++) printf(",%ld",I[j]+1);
                     printf(") = [\n\n");
                  }
                  else printf("%s = [\n",istr);
               }
            }

            printf(fmt.data,A.data[i].r,A.data[i].i);

            k=0; I[0]++;  
            while (I[k]>=A.SIZE[k] && k<l) { I[k]=0; ++I[++k]; }
         }
         printf("\n];\n");
      }
   }
   else if (s==1) {
      printf("  %s = ",istr);
      printf(fmt.data,data[0].r,data[0].i); printf("\n");
   }
   else printf("\n");

   return (int)s;
};

template<class T>
template<class TB, class TC>
wbarray<TC>& wbarray<T>::comm(
  const wbarray<TB> &B, wbarray<TC> &C, char aflag, char bflag
) const {

   if ((void*)&C==(void*)this || (void*)&C==(void*)&B) {
      wbsparray<TC> X; comm(B,X,aflag,bflag);
      return X.save2(C);
   }

   if (!isMatrix() || !B.isMatrix() || SIZE[0]!=B.SIZE[1] ||
       SIZE[0]!=SIZE[1] || B.SIZE[0]!=B.SIZE[1]) wblog(FL,
      "ERR %s() invalid operators for [A,B] (%s; %s)",
       FCT,SSTR_(this), SSTR(B)
   );

   Wb::MatProd(*this,B,C,aflag,bflag);        
   Wb::MatProd(B,*this,C,bflag,aflag,-1.,1.); 

   return C;
};

template<class T>
template<class TB, class TC>
wbarray<TC>& wbarray<T>::acomm(
  const wbarray<TB> &B, wbarray<TC> &C, char aflag, char bflag
) const {

   if (!isMatrix() || !B.isMatrix() ||
       SIZE[0]!=SIZE[1] || B.SIZE[0]!=B.SIZE[1]) wblog(FL,
      "ERR %s() invalid operators for [A,B] (%s; %s)",
       FCT,SSTR_(this), SSTR(B)
   );

   Wb::MatProd(*this,B,C,aflag,bflag);        
   Wb::MatProd(B,*this,C,bflag,aflag,+1.,1.); 

   return C;
};

template<class T>
template<class TB, class TC>
wbarray<TC>& wbarray<T>::contractMat(
   const char *F, int L, unsigned i1,
   const wbarray<TB> &M, unsigned i2, 
   wbarray<TC> &X
 ) const {

   unsigned r=SIZE.len, l=r-1;
   wbvector<size_t> S; wbperm P;
   wbarray<T> A;
   wbarray<T> MA;

   if (!M.isMatrix()) wblog (F_L,
      "ERR %s() rank 2 object required (%d)",FCT,M.SIZE.len);
   if (i1<1 || i1>SIZE.len || i2<1 || i2>2) wblog(F_L,
      "ERR %s() invalid contraction indizes [%d %d, %d %d]", i1,i2,
       FCT,rank(), M.rank());

   if (i1==1) {
      toMatrixRef(MA,i1-1,1,P); SIZE.select(P,S); l=0;
      if (i2==2) 
           { Wb::MatProd(M,MA,X    ); S[l]=M.SIZE[0]; }
      else { Wb::MatProd(M,MA,X,'T'); S[l]=M.SIZE[1]; }
   }
   else {
      toMatrixRef(MA,i1-1,2,P); SIZE.select(P,S);
      if (i2==1) 
           { Wb::MatProd(MA,M,X        ); S[l]=M.SIZE[1]; }
      else { Wb::MatProd(MA,M,X,'N','T'); S[l]=M.SIZE[0]; }
   }

   X.Reshape(S).Permute(P,'I');

   return X;
};

template<class T>
template<class TB> inline
wbarray<T>& wbarray<T>::Contract(char *idx1,
   const wbarray<TB> &B, char *idx2,
   const wbperm &pfinal
){
   wbarray<T> A(*this);

   WBINDEX i1, i2;
   i1 = Str2Idx(idx1,1); 
   i2 = Str2Idx(idx2,1);

   A.contract(i1, B, i2, *this, pfinal);
   return *this;
}

template<class T>
template<class TB, class TC> inline
wbarray<TC>& wbarray<T>::contract(const char *F, int L,
   const char *idx1, const wbarray<TB> &B, const char *idx2,
   wbarray<TC> &Q, const wbperm &pfinal
 ) const {

   wbvector<unsigned> i1, i2;
   i1 = Str2Idx(idx1,1); 
   i2 = Str2Idx(idx2,1); 

   return contract(F,L, i1,B,i2,Q,pfinal);
}

template<class T>
template<class TB> inline
wbarray<T>& wbarray<T>::Contract(const wbvector<unsigned> &i1,
   const wbarray<TB> &B, const wbvector<unsigned> &i2,
   const wbperm &pfinal
) {
   wbarray<T> A(*this); A.contract(i1, B, i2, *this, pfinal);
   return *this;
}

template<class T>
template<class TB, class TC>
wbarray<TC>& wbarray<T>::contract( 
   const char *F, int L, const iTags &ta,
   const wbarray<TB> &B, const iTags &tb, wbarray<TC> &C, iTags &tc,
   const wbperm &pfinal, T afac, TC cfac
) const {

   ctrIdx ia,ib;

   ta.getCtr(F_L,tb,ia,ib,&tc);
   contract(F_L,ia,B,ib,C,pfinal,afac,cfac);

   tc.Permute(pfinal);
   return C;
};

template<class TA>
template<class TB, class TC>
wbarray<TC>& wbarray<TA>::contract( 
   const char *F, const int L,
   const ctrIdx &ica, const wbarray<TB> &B, const ctrIdx &icb,
   wbarray<TC> &C0, const wbperm &P, 
   TA afac, 
   TC cfac  
) const {

   if ((void*)(&C0)==(void*)this || (void*)(&C0)==(void*)(&B)) {
      wbarray<TC> Ci; 
      if (cfac) {
         bool q1=((void*)(&C0)==(void*)this),
              q2=((void*)(&C0)==(void*)(&B));
         wblog(FL,"ERR %s() got cfac=%s with self-reference (%s)",
         FCT, NSTR(cfac), q1 && q2 ? "A=B=C" : (q1 ? "A=C" : "B=C"));
      }
      contract(F,L,ica,B,icb,Ci,P,afac,cfac);
      return Ci.save2(C0);
   }

#ifdef WB_CLOCK
   Wb::Clock clk("arr:contract",0); 
#endif

   unsigned i,s; size_t one=1;
   char aflag, bflag;
   wbvector<size_t> S1,S2;
   wbarray<TA> MA;
   wbarray<TB> MB;
   wbarray<TC> Ci;
   wbperm p1,p2;

   wbvector<unsigned> i1(ica); 
   wbvector<unsigned> i2(icb);

   char gotC_=1; 
   if (C0.isEmpty()) {
      if (isEmpty() || B.isEmpty()) { return C0; }
      if (cfac && cfac!=TC(1)) wblog(FL,
         "ERR %s() got cfac=%s with empty C",FCT,NSTR(cfac));
      if (C0.mtype && C0.mtype!=Wb::MEX_RETURN) wblog(FL,
         "WRN %s() unexpected mt=%s",FCT,Wb::MTYPE_STR[C0.mtype]);
      SWAP(Ci.mtype,C0.mtype); 
      gotC_=0;
   }

   if (i1.isEmpty() || i2.len!=i1.len) wblog(F,L,
      "ERR %s() index set length mismatch [%s] / [%s]",FCT,STR(i1),STR(i2));
   for (s=SIZE.len, i=0; i<i1.len; ++i) { if (i1[i]>=s) wblog(F,L,
      "ERR %s() index out of bounds (A @ [%s]/%d)",FCT,STR(i1+1),s); }
   for (s=B.SIZE.len, i=0; i<i2.len; ++i) { if (i2[i]>=s) wblog(F,L,
      "ERR %s() index out of bounds (B @ [%s]/%d)",FCT,STR(i2+1),s); }
   for (i=0; i<i1.len; ++i) {
      if (SIZE[i1[i]]!=B.SIZE[i2[i]]) wblog(F,L,
      "ERR wbarray::%s() incompatible data size\n[%s] @ %d <> [%s] @ %d",
      FCT,SSTR_(this),i1[i]+1,SSTR(B),i2[i]+1); }

     toMatrixRef(FL,MA,i1,2,aflag);    SIZE.getI(i1,S1); 
   B.toMatrixRef(FL,MB,i2,1,bflag);  B.SIZE.getI(i2,S2);

   if ((i=S1.len+S2.len)<2) {
      if (!i)           { S2.init2val(2,one); } 
   }

   if (!afac) { 
      WBINDEX S(S1,S2); if (P) S.Permute(P);
      if (gotC_) {
         if (C0.SIZE!=S) wblog(F_L,
            "ERR %s() size mismatch %s / %s",FCT,SSTR(S),SSTR(C0));
         C0*=cfac;
      }
      else {
         Ci.init(S).save2(C0); 
      }
      return C0;
   }

   char isPtrans = P.isCyclic2F(S1.len,S2.len,'i');

   if (aflag<0 || aflag>1) wblog(FL,
      "ERR %s() invalid aflag=%d !?",FCT,aflag);
   if (bflag<0 || bflag>1) wblog(FL,
      "ERR %s() invalid bflag=%d !?",FCT,bflag);

   if (isPtrans) { aflag=!aflag; bflag=!bflag; }

   if (ica.conj && ISCOMPLX_(TA)) { 
      if (aflag)
           { aflag='C'; }
      else { aflag='N'; MA.Instantiate().Conj(); }
   }  else { aflag=(aflag ? 'T':'N'); }

   if (icb.conj && ISCOMPLX_(TB)) { 
      if (bflag)
           { bflag='C'; }
      else { bflag='N'; MB.Instantiate().Conj(); }
   }  else { bflag=(bflag ? 'T':'N'); }

   if (!isPtrans) {
      Wb::MatProd(MA,MB,Ci, aflag, bflag);
      Ci.Reshape(UVEC(S1,S2)).Permute(P);
   }
   else {
      Wb::MatProd(MB,MA,Ci, bflag, aflag);
      Ci.Reshape(UVEC(S2,S1));
   }

   if (gotC_) {
      if (!Ci.sameSize(C0)) wblog(F,L,
         "ERR %s() size mismatch\n%s cannot add %s / %s",
         FCT,SHORT_FL,SSTR(Ci),SSTR(C0));
      return C0.Plus(Ci, TC(afac), 0, cfac); 
   }
   else {
      if (cfac && cfac!=TC(1)) wblog(FL, 
         "WRN %s() got cfac=%s with empty input",FCT,NSTR(cfac));
      if (afac!=TA(1)) { Ci*=afac; }
      return Ci.save2(C0);
   }
};

template<class T>
wbarray<T>& wbarray<T>::trace(
   unsigned i1, unsigned i2, wbarray<T> &C) const {   

   size_t i,j,k,d,q,s,r,l;
   WBINDEX S,I0,I;

   if (&C==this) wblog (FL,
   "ERR %s() output space coincides with *this",FCT);

   if (i1==i2 || !i1 || !i2 || i1>SIZE.len || i2>SIZE.len ||
       SIZE[i1-1]!=SIZE[i2-1]) {
       wblog (FL,"ERR %s() invalid trace indizes (%d,%d; %s)",
       FCT,i1, i2, SSTR_(this)); return C;
   }

   i1--; i2--; 

   d=SIZE[i1];

   S.init(SIZE.len-2);
   for (k=i=0; i<SIZE.len; ++i) if (i!=i1 && i!=i2) S[k++]=SIZE[i];
   if (S.len==0) { S.init(1); S[0]=1; }

   if (C.isEmpty()) C.init(S); else {
      if (S!=C.SIZE) wblog(FL,
         "ERR %s() size mismatch of I/O array\n%s: [%s] vs. [%s]",
          FCT,SHORT_FL,SSTR(C), STR(S)
      );
   }

   I.init(S.len); s=S.prod(0); I0.init(SIZE.len);
   r=SIZE.len-1; l=S.len-1;

   for (i=0; i<s; ++i) {
      for (q=k=0; k<SIZE.len; ++k) if (k!=i1 && k!=i2) I0[k]=I[q++];
      for (q=0; q<d; ++q) {
          I0[i1]=I0[i2]=q;
          for (j=I0[r],k=r-1; k<r; k--) j = j*SIZE[k] + I0[k];

          C.data[i]+=data[j];
       }

       k=0; I[0]++;
       while (I[k]>=S[k] && k<l) { I[k]=0; ++I[++k]; }
   }

   return C;
}

template<class T>
wbarray<T>& wbarray<T>::trace( 
   const ctrIdx &I1, const ctrIdx &I2, wbarray<T> &C) const {

   size_t i,j,k,l, kc,kx, lc,lx,l0, Nc,Nx, *sx, *sc;
   const size_t *s0=SIZE.data;

   unsigned r0=SIZE.len, *mx;
   const unsigned *i1=I1.data, *i2=I2.data;

   WBINDEX Sx, Sc, I0, Ix, Ic; widx_t *i0, *ix, *ic;
   wbvector<unsigned > Mx;

   if (&C==this) wblog (FL,
      "ERR %s() output space coincides with *this !?",FCT);
   if (I1.len!=I2.len || 2*I1.len>r0) wblog (FL,
      "ERR %s() invalid index sets (%s <> %s)",FCT,STR(I1),STR(I2));

   Mx.init(r0);     mx=Mx.data;
   Sx.init(I1.len); sx=Sx.data; 

   for (i=0; i<I1.len; ++i) {
      j=i1[i]; if (j>=r0 || mx[j]++) { break; }
      j=i2[i]; if (j>=r0 || mx[j]++) { break; }

      sx[i]=s0[i1[i]]; if (sx[i]!=s0[i2[i]]) { break; }
   }
   if (i<I1.len) wblog (FL,
      "ERR %s() invalid index set\nhaving %s with %s <> %s",
      FCT,SSTR_(this),STR(I1),STR(I2)
   );

   if (!I1.len) {
      if (C.isEmpty()) C=(*this); else C+=(*this);
      return C;
   }

   k=r0-2*I1.len; Sc.init(k?k:1); sc=Sc.data;
   if (!k) { sc[0]=1; }
   else { for (l=i=0; i<r0; ++i) { if (!mx[i]) sc[l++]=s0[i]; }}

   Nx=Sx.prod(0); lx=Sx.len-1;
   Nc=Sc.prod(0); lc=Sc.len-1; l0=r0-1;

   if (C.isEmpty()) C.init(Sc); else {
      if (Sc!=C.SIZE) wblog(FL,
         "ERR %s() size mismatch of I/O array\n%s: [%s] vs. [%s]",
         FCT,SHORT_FL,SSTR(C), STR(Sc)
      );
   }

   I0.init(r0);     i0=I0.data;
   Ic.init(Sc.len); ic=Ic.data;
   Ix.init(Sx.len); ix=Ix.data;

   for (kc=0; kc<Nc; ++kc) {
       for (l=i=0; i<r0; ++i) if (!mx[i]) i0[i]=ic[l++];

       Ix.set(0);
       for (kx=0; kx<Nx; ++kx) {
          for (i=0; i<Ix.len; ++i) {
             i0[i1[i]] = i0[i2[i]] = ix[i]; 
          }
          for (k=i0[l0], i=l0-1; i<l0; i--) { k = k*s0[i] + i0[i]; }

          C.data[kc] += data[k];

          i=0; ++ix[0]; 
          while (ix[i]>=sx[i] && i<lx) { ix[i]=0; ++ix[++i]; }
       }

       i=0; ic[0]++; 
       while (ic[i]>=sc[i] && i<lc) { ic[i]=0; ++ic[++i]; }
   }

   return C;
};

template<class T>
wbarray<T> TestContract(
   const wbarray<T> &A, const WBINDEX i1,
   const wbarray<T> &B, const WBINDEX i2,
   const wbperm &pfinal
){
   WBINDEX Si,Sj,Sk, I,J,K, Ia,Ib,Ic, Iq; WBPERM Pa,Pb;
   size_t i,j,k, si,sj,sk, li,lj,lk, q, ra=A.SIZE.len, rb=B.SIZE.len;
   wbarray<T> C; 

   T dbl;

   if (i1.isEmpty() || i2.len!=i1.len) wblog (FL,
      "ERR %s() invalid index set {[%s], [%s]}",
       FCT,i1.toStr().data, i2.toStr().data);
   if (i1.anyGE(ra) || i2.anyGE(rb)) wblog (FL,
      "ERR %s() index out of range\n[%s; %d], [%s; %d]",
       FCT,i1.toStr().data,ra,i2.toStr().data,rb);

return A;

   for (k=0; k<i1.len; ++k) {
   if (A.SIZE[i1[k]]!=B.SIZE[i2[k]]) wblog (FL,
      "ERR %s() incompatible objects {[%s], [%s]}\n[%s] <> [%s] at %d",
       FCT,i1.toStr().data, i2.toStr().data, SSTR(A), SSTR(B), k);
   }

   Si=A.SIZE.Skip(i1); si=Si.len>0 ? Si.prod(0) : 1; li=Si.len-1;
   Sj=A.SIZE [i1];     sj=Sj.len>0 ? Sj.prod(0) : 1; lj=Sj.len-1;
   Sk=B.SIZE.Skip(i2); sk=Sk.len>0 ? Sk.prod(0) : 1; lk=Sk.len-1;

   Pa = Index(0,ra-1).Move2end  (i1); 
   Pb = Index(0,rb-1).Move2front(i2);

   Ia.init(ra); Ib.init(rb); C.init(Si.append(Sk));

   if (C.numel()==0) C.init(1);

   I.init(Si.len);
   for (i=0; i<si; ++i) {

       K.init(Sk.len);
       for (k=0; k<sk; ++k) {

           J.init(Sj.len);
           for (dbl=0, j=0; j<sj; ++j) {
               Iq=I.append(J); for (q=0; q<Iq.len; ++q) Ia[Pa[q]]=Iq[q];
               Iq=J.append(K); for (q=0; q<Iq.len; ++q) Ib[Pb[q]]=Iq[q];

               if (j==sj-1) Ic=I.append(K);

               dbl += A(Ia) * B(Ib);

               if (J.len>0) {
                   q=0; J[0]++; 
                   while (J[q]>=Sj[q] && q<lj) { J[q]=0; ++J[++q]; }
               }
           }
           C(Ic) = dbl;

           if (K.len>0) {
               q=0; K[0]++; 
               while (K[q]>=Sk[q] && q<lk) { K[q]=0; ++K[++q]; }
           }
       }
       if (I.len>0) {
           q=0; I[0]++; 
           while (I[q]>=Si[q] && q<li) { I[q]=0; ++I[++q]; }
       }
   }

   if (pfinal) C.Permute(pfinal);

   return C;
};

int tensor_ref_::check(const char *F, int L) const {
   int q=0;
   if (Ar) {
      if (!Az) { q=1; }
      else { char s[64]; q=-1; 
         snprintf(s,64,"%s() multiple refs set (%p/%p)",FCT,Ar,Az);
         if (F) wblog(F,L,"ERR %s",s);
         else   wblog(F_L,"WRN %s",s);
      }
   }
   else if (Az) { q=2; }

   if (q<=0 && F) wblog(FL,"ERR %s() invalid array ref (q=%d)",FCT,q);
   return q;
};

int tensor_ref_::rank(unsigned r) const {
   int q=-1;
   if (Ar) {
      if (Az) wblog(FL,"ERR %s() multiple refs set (%p/%p)",FCT,Ar,Az);
      q=(Ar->SIZE.len==r);
   }
   else if (Az) { q=(Az->SIZE.len==r); }
   return q;
};

bool tensor_ref_::sameSize(const wbvector<size_t> &S) const {
   bool q=0;
   if (Ar) {
      if (Az) wblog(FL,"ERR %s() multiple refs set (%p/%p)",FCT,Ar,Az);
      q=(Ar->SIZE.isEqual(S));
   }
   else if (Az) { q=(Az->SIZE.isEqual(S)); }
   return q;
};

wbstring tensor_ref_::sizeStr() const { 

   wbstring s; 
   if (Ar) {
      if (Az) wblog(FL,"ERR %s() multiple refs set (%p/%p)",FCT,Ar,Az);
      s=SSTR_(Ar);
   }
   else if (Az) { s=SSTR_(Az); }
   else { s="[]"; }

   return s;
};

const wbvector<size_t>* tensor_ref_::getSIZE() const {

   const wbvector<size_t> *s=NULL;
   if (Ar) {
      if (Az) wblog(FL,"ERR %s() multiple refs set (%p/%p)",FCT,Ar,Az);
      s=&(Ar->SIZE);
   }
   else if (Az) { s=&(Az->SIZE); }

   return s;
};

template<>
void tensor_ref_::contract(
   const char *F, int L, const ctrIdx &ia,
   const tensor_ref_ &b, const ctrIdx &ib, wbarray<double> &C
 ) const {

   if (Az || b.Az) wblog(FL,
      "ERR %s() got complex array refs (%p/%p)",FCT,Az,b.Az);
   if (!Ar || !b.Ar) wblog(FL,
      "ERR %s() got NULL array refs (%p/%p)",FCT,Ar,b.Ar);
   Ar->contract(F_L,ia,*b.Ar,ib,C,wbperm(),1.,0.);
};

template<>
void tensor_ref_::contract(
   const char *F, int L, const ctrIdx &ia,
   const tensor_ref_ &b, const ctrIdx &ib, wbarray<wbcomplex> &C
 ) const {

   wbperm P; wbcomplex zero(0);
   check(F_L); b.check(F_L);

   if (Ar) {
      if (b.Ar)
           { Ar->contract(F_L,ia,*b.Ar,ib,C, P,1.,zero); }
      else { Ar->contract(F_L,ia,*b.Az,ib,C, P,1.,zero); }
   }
   else {
      if (b.Ar)
           { Az->contract(F_L,ia,*b.Ar,ib,C, P,1.,zero); }
      else { Az->contract(F_L,ia,*b.Az,ib,C, P,1.,zero); }
   }
};

char tensorRef_::check(const char* F, int L) const {
   char q=0, r; 
   if ((!S && !R) || !it) { q|=1; if (F) wblog(F_L,
      "ERR %s() got uninitialized data (%d/%d)",FCT,S.len,it.len); }
   if (S.len && S.len!=it.len) { q|=2; if (F) wblog(F_L,
      "ERR %s() length mismatch (S: %d/%d)",FCT,S.len,it.len); }
   if (!(r=R.rank(it.len))) { q|=2; if (F) wblog(F_L, 
       "ERR %s() length mismatch (A: %d/%d)",FCT,r,it.len); }
   if (ID.len && ID.len!=it.len) { q|=4; if (F) wblog(F_L,
      "ERR %s() length mismatch (ID: %d/%d)",FCT,ID.len,it.len); }
   if (S.len && r>0 && !R.sameSize(S)) { q|=8; if (F) wblog(F_L,
      "ERR %s() SIZE mismatch %s / %s",FCT,SSTR(S),SSTR(R)); }
   return q;
};

char tensorRef_::Check(const char* F, int L, unsigned &k) {
   if (!id) { id=(++k); }
   if (!R) { if (F) wblog(FL,"ERR %s() got NULL array reference",FCT); }
   else if (!S.len) { S.init2ref(*R.getSIZE()); }

   return check(F_L);
};

template<class T>
int tensorRef_::contract(const char *F, int L,
   const tensorRef_ &b, tensorRef_ &c,
   double *flops,
   wbarray<T> *C 
 ) const {

   int nc=0; 

   if (this==&c || &b==&c) { tensorRef_ x;
      nc=contract(F,L,b,x,flops,C);
      x.save2(c); return nc;
   }

   unsigned i,k,l=0;

   ctrIdx ia, ib;
   wbvector<char> ma, mb;

   if ((i=(check() | b.check()))) wblog(F_L,
      "ERR %s() %s input\n%s\n%s", FCT,
      i&1 ? "got empty":"inconsistent", STR_(this),STR(b));

   nc=this->it.getCtr(
      flops? NULL : (F? F : __FILE__), 
      flops?    1 : (L? L : __LINE__),
      b.it, ia,ib,    
      &c.it, &ma, &mb 
   );

   if (nc<1 && !flops) wblog(F_L,
      "ERR %s() invalid contraction (empty match)",FCT);

   if (ma.len!=S.len || mb.len!=b.S.len) wblog(FL,"ERR %s() "
      "length mismatch %d/%d ; %d/%d",FCT,ma.len,S.len,mb.len,b.S.len);
   for (i=0; i<ia.len; ++i) { if (S[ia[i]]!=b.S[ib[i]]) { wblog(FL,
      "ERR %s() size mismatch\n%-40s @ %s\n%-40s @ %s",
      FCT, STR_(this), STR(ia), STR(b), STR(ib));
   }}

   l=c.it.len;
   c.S.init(l); c.ID.init(l);

   k=(id<<8);
   for (l=i=0; i<S.len; ++i) { if (!ma[i]) {
      if (c.it[l]!=  it[i]) wblog(FL,"ERR %s() %d/%d: %s / %s",
         FCT,i+1,  S.len,STR(c.it[l]),  STR(it[i]));
      c.ID[l]=(ID.len ? ID[i] : k+i);
      c.S[l]=  S[i]; ++l;
   }}

   k=(b.id<<8);
   for (i=0; i<b.S.len; ++i) { if (!mb[i]) {
      if (c.it[l]!=b.it[i]) wblog(FL,"ERR %s() %d/%d: %s / %s",
         FCT,i+1,b.S.len,STR(c.it[l]),STR(b.it[i]));
      c.ID[l]=k+i;
      c.S[l]=b.S[i]; ++l;
   }}
   if (l!=c.S.len) wblog(FL,"ERR %s() %d/%d",FCT,l,c.S.len);

   if (flops) { double f=R.xflops()*b.R.xflops();
      if (f==2) { f=4; }

      for (i=0; i< ia.len; ++i) { f *= S[ia[i]]; }
      for (i=0; i<c.S.len; ++i) { f *= c.S[i]; } 
      *flops = f;
   }

   if (C) {
      ia.conj=  conj;
      ib.conj=b.conj;

      R.contract(F_L,ia,b.R,ib,*C);
      c.R.init(*C);
   }

   return nc; 
};

double tensorRefs::getOptimalCtrOrder(const char *F, int L,
    wbperm &P,  
    wbperm *pfin,
    char useP,  
    char vflag
) { 

   double fmin=-1, flops=0, f;

   if (len<2 || len>8) wblog(FL,
      "ERR %s() got %d tensors (max. 8 supported)",FCT,len);

   Check(F_L); 

   unsigned ip=0,i0=0, j; int nc; wperm_t *p;
   wbMatrix<wperm_t> PP;
   tensorRef_ R, Rmin;

   if (vflag) {
      if (vflag=='v') vflag=1; else
      if (vflag=='V') vflag=3; else
      if (vflag>8) wblog(FL,"WRN %s() got vflag=%s",FCT,cSTR(vflag));
   }
   if (useP) { 
      if (P) { if (P.len!=len) wblog(FL,
         "ERR %s() size mismatch (%d/%d)",FCT,P.len,len);
      } else { useP=0; }
   }
   if (useP)
        { PP.init(1,len,P.data); }
   else { Wb::getPerms(PP,len); }
   p=PP.data;

   for (; ip<PP.dim1; ++ip, p+=PP.dim2) {
      if (p[0]>p[1]) { continue; }

      nc=data[p[0]].contract(FL,data[p[1]], R,f);
      if (!nc) { continue; }
      flops=f;

      for (j=2; j<len; ++j) { 
         nc=R.contract(FL,data[p[j]], R,f);
         if (!nc) { break; }
         flops+=f;
      }

      if (!nc) { continue; } 

      if (!ip || (fmin>flops)) { fmin=flops; R.save2(Rmin); i0=ip; }

      if (vflag>1) { PRINTF( 
         "%2d [%s] %-64s @%8ld flops\n",
         ip, STR(wbperm().init(PP.dim2,p,'r')), STR(R), size_t(flops));
      }
   }
   if (vflag) { print();
      wblog(FL,"--> %s() picking ip=%d/%ld\n\n",FCT,i0,PP.dim1);
   }

   P.init(PP.dim2,PP.rec(i0));

   if (pfin) {
      Rmin.ID.Sort(*pfin);
   }

   return fmin;
};

void tensorRefs::print(const char *F, int L) const {

   if (F || !len) wblog(F_L,
      "--- tensorRefs listing%s",len? "":" (empty)");

   if (len) { unsigned i=0; char istr[16];
      PRINTF("\n");
      for (; i<len; ++i) {
         snprintf(istr,16,"  %2d/%ld :",i+1,len);
         data[i].print(istr, F?'v':0); 
      }
      PRINTF("\n");
   }
};

template<class T> inline
bool wbarray<T>::requiresDataPerm(const wbperm &P, char lflag) const {

   if (!P.len) { return 0; }

   const wperm_t *p=P.data; const size_t *s=SIZE.data;
   unsigned i,j, j_=-1;

   if (P.len!=SIZE.len) {
      if (lflag<0 || lflag>3) wblog(FL,
         "WRN %s() got lflag=%s",FCT,cSTR(lflag));
      if (!( lflag & (P.len<SIZE.len ? 1 : 2) )) { wblog(FL,"ERR %s() "
         "invalid permutation (len=%d/%d; %d)",FCT,P.len,SIZE.len,lflag);
      }
   }
   if (numel()<=1) { return 0; }

   for (i=0; i<P.len; ++i) { if ((j=p[i])<SIZE.len && s[j]>1) {
      if (j_<j || int(j_)<0)
           { j_=j; } 
      else { return 1; }
   }}

   return 0;
};

template<class T> inline
void wbarray<T>::toMatrixRef(
   wbarray<T> &A,               
   const wbvector<unsigned> &I, 
   int pos,                     
   wbperm &P                    
) const {

   size_t s1,s2;
   groupIndizes_P(I,pos,P,&s1,&s2);

   if (!requiresDataPerm(P)) { A.init2ref(*this);  } 
   else { permute(A,P); }

   A.Reshape(s1,s2);
};

template<class T> inline
bool wbarray<T>::toMatrixRef(const char *F, int L,
   wbarray<T> &A,               
   const wbvector<unsigned> &I, 
   int pos,                     
   char &tflag
) const {

   size_t i,s1,s2; wbperm P;

   int k[2]={pos, pos==1 ? 2 : 1};

   for (i=0; i<2; ++i) {
      groupIndizes_P(I,k[i],P,&s1,&s2);

      if (!requiresDataPerm(P)) {
         tflag = (pos==k[i] ? 0 : 1);
         A.init2ref(*this);
         A.Reshape(s1,s2);  
         return 0;
      }
   }

   if (&A==this) wblog(F_L,"ERR %s() Wb::overlaping output space",FCT);

   toMatrixRef(A,I,pos,P); tflag=0;

   return 1;
};

template<class T>
void wbarray<T>::toMatrixRef(
   wbarray<T> &A,     
   unsigned K,        
   const wbperm &P0   
) const {

   size_t i,s1,s2; wbperm P;

   if (K>SIZE.len) wblog(FL,
      "ERR %s() index out of bounds (%d/%d)",FCT,K,SIZE.len);

   if (!P0.len) { P.Index(SIZE.len); }
   else {
      if (P0.len!=SIZE.len) wblog(FL,
         "ERR %s() rank mismatch (len=%d/%d)",FCT,P0.len,SIZE.len);
      P=P0;
   }

   for (s1=1,i=0; i<K; ++i) { s1*=SIZE[P[i]]; }
   for (s2=1; i<P.len; ++i) { s2*=SIZE[P[i]]; }

   if (!requiresDataPerm(P))
        { A.init2ref(*this);  } 
   else { permute(A,P); }

   A.Reshape(s1,s2);  
};

template<class T>
void wbarray<T>::toMatrixRefM(
   wbarray<T> &A,     
   unsigned K,        
   const wbperm &P0,  
   unsigned R,
   wbperm pfin
) const {

   size_t i,m=1,s1,s2; wbperm P;
   char use_ref=0;

   if (K>R) wblog(FL,
      "ERR %s() index out of bounds (K=%d/%d/%d)",FCT,K,R,SIZE.len);
   if (R<2 || R>SIZE.len) wblog(FL,
      "ERR %s() invalid rank (R=%d/%d; K=%d)",FCT,R,SIZE.len,K);
   if (SIZE.len>R+1) wblog(FL,"WRN %s() "
      "got multi-dimensional OM (%s @ %d; K=%d)",FCT,SSTR_(this),R,K);

   if (!P0.len) { P.Index(R); }
   else {
      if (P0.len!=R) wblog(FL,
         "ERR %s() rank mismatch (len=%d/%d/%d)",FCT,P0.len,R,SIZE.len);
      P=P0;
   }

   for (s1=1,i=0; i<K; ++i) { s1*=SIZE[P[i]]; }
   for (s2=1; i<P.len; ++i) { s2*=SIZE[P[i]]; }
   for ( ; i<SIZE.len; ++i) { m *=SIZE[  i ]; }

   use_ref=( requiresDataPerm(P) ? 0:1 );
   if (pfin.len) {
      if (pfin.len!=3) wblog(FL,"ERR %s() invalid pfin='%s'",FCT,STR(pfin));
      if (pfin[0]!=0 || m!=1) { use_ref=0; }
   }

   if (use_ref) 
        { A.init2ref(*this);  } 
   else { permute(A,P); }

   A.Reshape(s1,s2,m); 

   if (pfin.len) {
       A.Permute(pfin);
    }
};

template<class T> inline
wbarray<T>& wbarray<T>::save2(wbMatrix<T> &M) {

   if (M.data) { M.init(); }

   if (SIZE.len!=2) {
      if (SIZE.len) wblog(FL,"ERR invalid rank-%d",SIZE.len);
      return M;
   }

   check_consistency(FL);
   if (sptr) { M.data=sptr->save_dref(FL); } 

   M.dim1=SIZE[1]; M.dim2=SIZE[0]; SIZE.init();

   return M;
};

#ifdef QS_USING_MPFR

template<> inline
mxArray* wbarray<Wb::quad>::toMx() const {

   const char *fields[]={"S","data"};
   mxArray *a=mxCreateStructMatrix(1,1,2,fields);

   wbvector<Wb::quad> X(numel(),data,'r');

   mxSetFieldByNumber(a,0,0, SIZE.toMx());
   mxSetFieldByNumber(a,0,1, X.toMx());

   return a;
};

template<> inline
wbarray<Wb::quad>& wbarray<Wb::quad>::init(
   const char *F, int L, const mxArray *S,
   char ref QS_UNUSED_VAR,
   char vec QS_UNUSED_VAR
){
   if (!S || !mxIsStruct(S)) wblog(FL,"ERR %s() "
      "got invalid input `%s' / mpfr",FCT,S? mxGetClassName(S):"null");

   wbvector<size_t> sz(FL,mxGetFieldByNumber(S,0,0));
   wbvector<Wb::quad> X(mxGetFieldByNumber(S,0,1));

   if (X.len!=sz.prod(0)) wblog(F_L,
      "ERR %s() size mismatch %d/%d",FCT,(int)X.len,(int)sz.prod(0));
   if (X.len) init(sz,X.data);

   return *this;
};

#endif

template<class T>
void cell2mat(
   const wbMatrix< wbarrRef<T> > &C,
   wbarray<T> &M,  
   WBINDEX &D1, WBINDEX &D2
){

   size_t i,j,m,r,s,d1=0,d2=0,i0,j0;
   T *d0, *dd;
   double dc,fac;

   D1.init(C.dim1);
   D2.init(C.dim2);

   for (i=0; i<C.dim1; ++i)
   for (j=0; j<C.dim2; ++j) { if (C(i,j).D==NULL) continue;

       WBINDEX &S = C(i,j).D->SIZE;

       if (S.len==2) { d1=S[0]; d2=S[1]; } else
       if (S.len==1) { d1=d2=S[0]; }
       else wblog(FL,
         "ERR %s() requires 1/2-D arrays (got rank-%d).",FCT,S.len);

       if (C(i,j).tflag) SWAP(d1,d2);

       if (D2[j]==0) D2[j]=d2; else
       if (d2!=D2[j]) wblog(FL,
          "ERR %s() dimension mismatch (%d,%d: %d/%d).",FCT,i,j,d2,D2[j]);

       if (D1[i]==0) D1[i]=d1; else
       if (d1!=D1[i]) wblog(FL,
          "ERR %s() dimension mismatch (%d,%d: %d/%d).",FCT,i,j,d1,D1[i]);
   }

   M.init(D1.sum(),D2.sum());

   for (j0=j=0; j<C.dim2; ++j, j0+=d2) { d2=D2[j];
   for (i0=i=0; i<C.dim1; ++i, i0+=d1) { d1=D1[i];

      if (C(i,j).D==NULL) continue;

      d0=(C(i,j).D->data); dd=M.ref(i0,j0);
      fac=C(i,j).fac;
      dc =C(i,j).dc;

      if (C(i,j).D->SIZE.len==1) {
         if (C(i,j).tflag) wblog(FL,
            "WRN %s() tflag irrelevant for diag-matrix (%d,%d)",
             FCT,i+1,j+1);
         if (d1!=d2) wblog(FL,"ERR d1 != d2 !? (%,%d)",d1,d2);

         for (r=0; r<d1; ++r, dd+=(M.dim1+1))
         dd[r]+=(dc+fac*d0[r]);

         continue;
      }

      if (C(i,j).tflag==0) {

         if (fac==1.) {
            for (s=0; s<d2; ++s, dd+=M.dim1, d0+=d1)
            for (r=0; r<d1; ++r) dd[r]+=     d0[r]; 
         }
         else if (fac==-1.) {
            for (s=0; s<d2; ++s, dd+=M.dim1, d0+=d1)
            for (r=0; r<d1; ++r) dd[r]-=     d0[r];
         }
         else if (fac!=0.) {
            for (s=0; s<d2; ++s, dd+=M.dim1, d0+=d1)
            for (r=0; r<d1; ++r) dd[r]+= fac*d0[r];
         }
      }
      else {

         if (fac==1.) {
            for (s=0; s<d2; ++s, dd+=M.dim1)
            for (r=0; r<d1; ++r) dd[r]+=     d0[s+r*d2];
         }
         else if (fac==-1.) {
            for (s=0; s<d2; ++s, dd+=M.dim1)
            for (r=0; r<d1; ++r) dd[r]-=     d0[s+r*d2];
         }
         else if (fac!=0.) {
            for (s=0; s<d2; ++s, dd+=M.dim1)
            for (r=0; r<d1; ++r) dd[r]+= fac*d0[s+r*d2];
         }
      }
      if (dc!=0.) {
         if (d1!=d2) wblog(FL,
            "WRN %s() constant diagonal shift \n"
            "requires square block (%d,%d: %d,%d)",FCT,i+1,j+1,d1,d2);

         dd=M.ref(i0,j0); m=MIN(d1,d2);

         for (r=0; r<m; ++r, dd+=(M.dim1+1))
         dd[r]+=dc;
      }
   }}
};

template<class T> 
void wbarray<T>::getReal(wbarray<double> &R) const {
   R.initT(*this);
};

template<>
void wbarray<wbcomplex>::getReal(wbarray<double> &R) const {
   size_t i=0, n=numel();
   R.init(SIZE); for (; i<n; ++i) { R.data[i]=(double)data[i].r; }
};

template<class T>
void wbarray<T>::getImag(wbarray<double> &I QS_UNUSED_VAR) const {
   wblog(FL,"ERR wbarray::getImag not defined for type %s",TSTR(T));
};

template<> 
void wbarray<double>::getImag(wbarray<double> &I) const {
   I.init(SIZE);
};

template<>
void wbarray<wbcomplex>::getImag(wbarray<double> &I) const {
   size_t i=0, n=numel();
   I.init(SIZE); for (; i<n; ++i) { I.data[i]=data[i].i; }
};

template<class T>
void wbarray<T>::Conj() { return; }

template<>
void wbarray<wbcomplex>::Conj() {
   size_t i=0, n=numel();
   for (; i<n; ++i) { data[i].i = -data[i].i; }
};

template<class T>
void wbarray<T>::set(
   const wbarray<double> &R QS_UNUSED_VAR,
   const wbarray<double> &I QS_UNUSED_VAR
 ) { wblog(FL,"ERR wbarray::set(R,I) not defined for type %s",TSTR(T)); };

template<>
void wbarray<wbcomplex>::set(
   const wbarray<double> &R_, const wbarray<double> &I_) {

   size_t i=0, n;
   const double *R=R_.data, *I=I_.data;

   if (R) { init(R_.SIZE); n=numel();
      if (I) {
          if (!I_.sameSize(R_)) wblog(FL,"ERR %s() "
             "dimension mismatch %s / %s",FCT,SSTR(R_),SSTR(I_));
          for (; i<n; ++i) { data[i].set(R[i],I[i]); }
      }
      else { for (; i<n; ++i) { data[i].set(R[i],0.); }}
   }
   else if (I) {
      init(I_.SIZE); n=numel();
      for (; i<n; ++i) { data[i].set(0.,I[i]); }
   }
   else { init(); }
};

template <typename T> 
  void wbperm_helper<T>::permute(const T *src, T *dest, int np
) const {

   if (rk<2) {
      if (rk==1) { std::copy_n (src, sz[0], dest); }
      return;
   }

   widx_t m=2*m_blk, m2=m*m; 

   if (!src || !dest) wblog(FL,
      "ERR %s() got null input (s=%p, d=%p; rk=%d)",FCT,src,dest,rk);

   if (numel<8192) { 
      widx_t idest=-1;
      wbIndex I(rk, sz);
      while (++I) { dest[++idest]=src[I.serial(stride)]; }
      return;
   }

   if (l1 && (l1>1 || stride[0]!=1)) wblog(FL,
      "ERR %s() unexpected l1=%d with stride %d (l2=%d, rk=%d)",
      FCT,l1,stride[0],l2,rk);

   if (np<1) { np=1; }

   if (l1) {
      const unsigned k=2;  

      widx_t idest=-1, isrc=0, len=sz[k-1], step=stride[k-1];
      wbIndex I(rk-k, sz+k); 

      if (rk<3) wblog(FL, 
         "ERR %s() got rk=%d with l1=%d, l2=%d !?",FCT,rk,l1,l2);
      if (sz[0]<2) wblog(FL,"ERR %s() unexpected sz[0]=%d",FCT,sz[0]);

      if (np>1) { np=std::min(4,np); if ((int)len<4*np) { np=1; }}

      while (++I) { ++idest;
         isrc = I.serial(stride+k);
         #pragma omp parallel for num_threads(np)
         for (widx_t i=0; i<len; ++i) {
            std::copy_n( 
            src+isrc+i*step, sz[0], dest+(i+idest*len)*sz[0]);
         }
      }
      return;
   }

   unsigned i, l=0;
   int r_=rk-(l1+2); if (r_<1) { r_=1; }

   widx_t isrc, idest;
   widx_t sz_[3*r_], *stride_s = sz_+r_, *stride_d = sz_+2*r_;

   widx_t N_=0, M_;
   widx_t M=sz[l2], N=sz[l1]; 

   for (idest=1, i=0; i<=l1; ++i) { idest*=sz[i]; } 
   for (; i<rk; ++i) {
      if (i!=l2) { sz_[l] = sz[i];
         stride_s[l] = stride[i];
         stride_d[l] = idest;
         ++l;
      }
      else {
         N_=idest;  
      }
      idest*=sz[i]; 
   }

   M_=stride[l1];   
   if (!N_) wblog(FL,"ERR %s() failed to set N_",FCT);

   if (np>1) {  
      if (M*N > (1>>14) && M>=2*m) { 
         int np_=M/m; if (np>np_) { np=np_; }
      }
   }

   wbarray<T> Blk(m2,np);
   wbIndex I_(l,sz_);

   while (++I_) {
      isrc  = I_.serial(stride_s);
      idest = I_.serial(stride_d);

      #pragma omp parallel for num_threads(np)
      for (widx_t I=0; I<M; I+=m) {
         const T* s_; T* d_, *blk = Blk.data;
        #ifdef QS_USING_OMP
         blk += omp_get_thread_num()*m2;
        #endif

         widx_t i,j,J,jM_,n_,m_ = std::min(m,M-I); 
      for (J=0; J<N; J+=m) { n_ = std::min(m,N-J); 
         s_ = src + isrc + J*M_ + I;  
         for (j=0; j<n_; ++j) { jM_=j*M_;
         for (i=0; i<m_; ++i) { blk[j*m+i] = s_[jM_+i]; }}

         d_ = dest + idest + I*N_ + J;
         for (i=0; i<m_; ++i) { 
         for (j=0; j<n_; ++j) { d_[i*N_+j] = blk[j*m+i]; }}
      }}
   }
};

template <typename T>
mxArray* wbperm_helper<T>::toMx() const {

   const char *fields[] = { "rk","numel","sz","stride","blk" };
   mxArray *S=mxCreateStructMatrix(1,1,5,fields);
   wbvector<double> x(MAX(3U,rk));

   x.len=2; x[0]=rk; x[1]=rk_;
   mxSetFieldByNumber(S,0,0,x.toMx()); 

   x.len=1; x[0]=numel;
   mxSetFieldByNumber(S,0,1,x.toMx()); 

   if (rk) {
      x.len=rk;  
      Wb::cpyRange(x.data,sz,    rk); mxSetFieldByNumber(S,0,2,x.toMx());
      Wb::cpyRange(x.data,stride,rk); mxSetFieldByNumber(S,0,3,x.toMx());
   }

   x.len=3; x[0]=l1; x[1]=l2; x[2]=m_blk; 
   mxSetFieldByNumber(S,0,4,x.toMx());    

   return S;
};

#endif

