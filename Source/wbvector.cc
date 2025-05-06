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

#ifndef __WB_SIMPLE_VECTOR_CC__
#define __WB_SIMPLE_VECTOR_CC__

/* -------------------------------------------------------------------- */
/* -------------------------------------------------------------------- */

template <class T>
wbvector<T>& wbvector<T>::set(const wbvector<size_t> &I0, const T &x) {
   const size_t *I=I0.data;
   for (size_t i=0; i<I0.len; ++i) {
      if (I[i]>=len) wblog(FL,
         "ERR %s() index out of bounds (%d/%d)",FCT,I[i],len);
      data[I[i]]=x;
   }
   return *this;
};

template <class T>
wbvector<T>& wbvector<T>::add(const wbvector<size_t> &I0, const T &x) {
   const size_t *I=I0.data;
   for (size_t i=0; i<I0.len; ++i) {
      if (I[i]>=len) wblog(FL,
         "ERR %s() index out of bounds (%d/%d)",FCT,I[i],len);
      data[I[i]]+=x;
   }
   return *this;
};

template <class T>
bool wbvector<T>::anyUnequal(
   const wbvector<size_t> &I0, const T &x, size_t *k
 ) const {

   const size_t *I=I0.data;
   for (size_t i=0; i<I0.len; ++i) {
      if (I[i]>=len) wblog(FL,
         "ERR %s() index out of bounds (%d/%d)",FCT,I[i],len);
      if (data[I[i]]!=x) { if (k) (*k)=I[i]; return 1; }
   }
   return 0;
};

template <class T>
bool wbvector<T>::anyEqual(
   const wbvector<size_t> &I0, const T &x, size_t *k
 ) const {

   const size_t *I=I0.data;
   for (size_t i=0; i<I0.len; ++i) {
      if (I[i]>=len) wblog(FL,
         "ERR %s() index out of bounds (%d/%d)",FCT,I[i],len);
      if (data[I[i]]==x) { if (k) (*k)=I[i]; return 1; }
   }
   return 0;
};

template <class T>
wbvector<T>& wbvector<T>::Set(
   const wbvector<size_t> &I, const wbvector<T> &v
){
   if (!I.isUnique()) wblog(FL,"ERR index not unique");
   if (I.len!=v.len) wblog(FL,
      "ERR length mismatch (%d/%d)",I.len,v.len);

   for (size_t j,i=0; i<I.len; i++) { j=I[i];
      if (j>=len) wblog(FL,
      "ERR index out of bounds (%d: %d/%d)",i,j,len);

      data[j]=v[i];
   }
   return *this;
};

template <class T>
bool wbvector<T>::isUnique() const {
   wbvector v(*this); wbperm p; Wb::hpsort(v,p);
   for (size_t i=1; i<v.len; ++i) if (v[i]==v[i-1]) return 0;
   return 1;
};

template <>
bool wbvector<double>::isNormal() const {
   for (size_t i=0; i<len; ++i)
       if (data[i] && !isnormal(data[i])) return 0;
   return 1;
};

template <>
bool wbvector<float>::isNormal() const {
   for (size_t i=0; i<len; ++i)
       if (data[i] && !isnormal(data[i])) return 0;
   return 1;
};

template <> bool wbvector<unsigned>::isNormal() const { return 1; };
template <> bool wbvector<size_t  >::isNormal() const { return 1; };
template <> bool wbvector<int     >::isNormal() const { return 1; };
template <> bool wbvector<char    >::isNormal() const { return 1; };
template <> bool wbvector<long    >::isNormal() const { return 1; };

template <class T>
T wbvector<T>::normDiff2(const wbvector<T> &v, T fac) const {
   if (len!=v.len) wblog(FL,
      "ERR %s() vector length mismatch (%d/%d)",FCT,len,v.len);
   return Wb::rangeNormDiff2(data,v.data,len,fac);
};

template <class T>
double wbvector<T>::normDiff2_(const wbvector<T> &v, double fac) const {
   double n2=0; unsigned i=0;

   if (len!=v.len) wblog(FL,
      "ERR %s() vector length mismatch (%d/%d)",FCT,len,v.len);

   if (fac==1) {
      for (; i<len; ++i) { n2+=data[i].normDiff2(v.data[i]); }
   }
   else if (fac!=0) { T x;
      for (; i<len; ++i) {
         x=v.data[i]; x*=fac; n2+=data[i].normDiff2(x);
      }
   }
   return n2;
};

template <class T>
T wbvector<T>::cumsum_(wbvector<T> &S, char xflag) const {
   T Stot=0; 

   S.init( xflag? len+1 : len); 
   if (len>1) {
      size_t i=0, l=S.len-1; T *x=S.data;
      for (; i<l; ++i) { x[i+1]=x[i]+data[i]; }
      Stot=(xflag ? x[l] : x[l]+data[l]); 
   }
   else if (len) { 
      if (xflag) { S[1]=data[0]; }
      Stot=data[0];
   }
   return Stot;
};

template <class T>
T wbvector<T>::Cumsum(char zflag) const {
   if (!len) return 0;
   if (zflag) {
      T q,x=data[0]; data[0]=0;
      for (size_t i=1; i<len; ++i) {
         q=data[i]; data[i]=data[i-1]+x; x=q;
      }
      return data[len-1]+x;
   }
   else {
      for (size_t i=1; i<len; ++i) { data[i]+=data[i-1]; }
      return data[len-1];
   }
};

template <class T>
template <class T2> 
T wbvector<T>::cumsum0prod(wbvector<T2> &b, wbvector<T> &cs, char xflag) const {
   size_t i=0, l=(xflag ? len+1 : len); 
   if (b.len!=len) wblog(FL,
      "ERR %s() got length mismatch (%d/%d)",FCT,len,b.len);
   cs.init(l); if (!len) return 0;

   for (--l; i<l; ++i) { cs[i+1]=cs[i]+data[i]*b.data[i]; }
   return (xflag ? cs[i] : cs[i]+data[i]*b.data[i]); 
};

template <class T>
wbvector<T>& wbvector<T>::Eldiv(const wbvector<T> &v) const {
   if (len!=v.len) wblog(FL,
      "ERR %s() size mismatch (%d/%d)",FCT,len,v.len);
   for (size_t i=0; i<len; ++i) {
      data[i] /= v.data[i];
      if (v.data[i]==0) wblog(FL,
         "ERR %s() div/0 !?  (%g)",FCT,double(data[i])
      );
   }
   return *this;
};

template <class T>
wbvector<T>& wbvector<T>::tensorProd(
   const wbvector<T> &v2, wbvector<T> &vv) const {

   if (&vv==this || &vv==&v2) {
      wbvector<T> x;
      return this->tensorProd(v2,x).save2(vv);
   }

   vv.init(len*v2.len);
   if (vv) {
      size_t i, j=0, l=0;
      const T *d2=v2.data; T *d=vv.data;

      for (; j<v2.len; ++j) { 
      for (i=0; i<len; ++i, ++l) { d[l]=data[i]*d2[j]; }}
   }
   return vv;
};

template <class T>
wbvector<T>& wbvector<T>::RevertSigns() { 

   size_t i=0;
   for (; i<len; ++i) { if (data[i]) break; }
   if (i<len) {
      if (data[i]>0 && (-data[i])>0) wblog(FL,"ERR %s() "
         "got unsigned data type %s",FCT,TSTR(T));
      for (; i<len; ++i) data[i]=-data[i];
   }
   return *this;
};

template <class T>
void wbvector<T>::selectSU(const WBINDEX &I) {
    size_t i,n;
    WBINDEX mark(len);
    for (i=0; i<I.len; i++) {
       if (I[i]>=len)
       wblog(FL,"ERR index out of bounds (%d/%d).", i,len);
       mark[I[i]]++;
    }
    for (n=i=0; i<mark.len; i++)
    if (mark[i]) mark[n++]=i; 

    mark.Resize(n);
    Select(mark);
};

template <class T>
wbvector<T>& wbvector<T>::select(
    const WBINDEX &I, wbvector<T> &a) const {

    a.init(I.len); if (!I.len) { return a; }

    for (size_t i=0; i<a.len; ++i) {
       if (I[i]>=len) wblog(FL,
          "ERR %s() index out of bounds (%d/%d)",FCT,I[i],len);
       a.data[i]=data[I[i]];
    }
    return a;
};

template <class T>
void wbvector<T>::select(const WBINDEX &I, T* r) const {
    if (I.len==0) return;

    for (size_t i=0; i<I.len; i++) {
       if (I[i]>=len) wblog(FL,
       "ERR Index out of bounds (%d/%d).", I[i],len);

       r[i]=data[I[i]];
    }
};

template <class T>
wbvector<T>& wbvector<T>::Select(const WBINDEX &I) {

    if (isref) wblog(FL,"ERR must not resize vector reference!");
    size_t i; T* d0=data;

    if (I.len==0) { init(); return *this; }
    WB_NEW(data,I.len);

    for (i=0; i<I.len; ++i) {
       if (I[i]>=len) wblog(FL,
          "ERR index out of bounds (%d/%d)", I[i],len);
       data[i]=d0[I[i]];
    }

    len=I.len;
    WB_DELETE(d0); return *this;
};

template <class T>
wbvector<T>& wbvector<T>::BlockSelect(
    const WBINDEX &I, const WBINDEX &D
){
    WBINDEX dc;
    size_t i,j,l,n; T *d0=data, *d;

    if (isref) wblog(FL,"ERR shall not resize vector reference!");
    if (I.len==0) { init(); return *this; }

    for (n=i=0; i<I.len; i++) {
       j=I[i]; if (j>=D.len) wblog(FL,
         "ERR index out of bounds (%d/%d).",j,D.len);
       n+=D[j];
    }

    i=D.cumsum_(dc);
    if (i!=len) wblog(FL,"ERR severe size mismatch (%d/%d; %d)",i,len,n);

    WB_NEW(data,n); len=n;

    for (n=l=i=0; i<I.len; ++i, l+=n) { j=I[i];
       d=d0+dc[j]; n=D[j];
       for (j=0; j<n; j++) data[l+j]=d[j];
    }

    WB_DELETE(d0); return *this;
};

template <class T>
template <class T1, class T2, class T3, class T4>
wbvector<T>& wbvector<T>::Cat(
   const wbvector<T1> *v1, const wbvector<T2> *v2,
   const wbvector<T3> *v3, const wbvector<T3> *v4
){
   return Cat(
      v1 ? v1->data : (T1*)0, v1 ? v1->len : 0,
      v2 ? v2->data : (T2*)0, v2 ? v2->len : 0,
      v3 ? v3->data : (T3*)0, v3 ? v3->len : 0,
      v4 ? v4->data : (T4*)0, v4 ? v4->len : 0
   );
};

template <class T>
template <class T1, class T2, class T3, class T4>
wbvector<T>& wbvector<T>::Cat(
   const T1* v1, size_t l1, const T2* v2, size_t l2,
   const T3* v3, size_t l3, const T4* v4, size_t l4
){
   wbvector<T> X(l1+l2+l3+l4); 
   size_t i; T *x=X.data;

   for (i=0; i<l1; ++i) { x[i]=T(v1[i]); }; x+=l1;
   for (i=0; i<l2; ++i) { x[i]=T(v2[i]); }; x+=l2;
   for (i=0; i<l3; ++i) { x[i]=T(v3[i]); }; x+=l3;
   for (i=0; i<l4; ++i) { x[i]=T(v4[i]); }

   return X.save2(*this);
};

template <class T>
wbvector<T>& wbvector<T>::Cat(const wbvector<wbvector<T> > &vv) {
   size_t i,n=0; T *d;

   for (i=0; i<vv.len; i++) n+=vv[i].len;
   init(n); d=data;

   for (i=0; i<vv.len; i++) {
      memcpy(d,vv[i].data,sizeof(T)*vv[i].len); d+=vv[i].len;
   }
   return *this;
};

template <class T>
wbvector<T>& wbvector<T>::Cat(const wbvector<wbvector<T>* > &vv) {
   size_t i,n=0; T *d;

   for (i=0; i<vv.len; i++) n+=(vv[i] ? vv[i]->len : 0);
   init(n); d=data;

   for (i=0; i<vv.len; i++) { if (vv[i]) {
       const wbvector<T> &v=*vv[i];
       memcpy(d,v.data,sizeof(T)*v.len); d+=v.len;
   }}
   return *this;
};

template <class T>
void wbvector<T>::info(const char *istr) const {

   size_t l=0, n=32; char sstr[n];
   size_t b=len*sizeof(T);
   wbstring tstr;

   if (typeid(T)==typeid(double)) tstr="double";
   else tstr.init(TSTR(T));

   if (b<(1<<10)) l=snprintf(sstr,n,"%ld ",   b); else
   if (b<(1<<20)) l=snprintf(sstr,n,"%.3g kB",b/double(1<<10)); else
                  l=snprintf(sstr,n,"%.3g MB",b/double(1<<20));
   if (l>=n) wblog(FL,"ERR %s() string out of bounds (%d/%d)",FCT,l,n);

   printf("  %-12s %10ld %12s  @ %p  %s vector%s\n",
   istr, len, sstr, (void*)data, tstr.data, isref ? "  ***ISREF***" : "");
}

template <class T>
void wbvector<T>::print(const char *istr, char mflag) const {

    mxArray *a;

    if (mflag==1) { 
        char fmt[8];

        if (typeid(T)==typeid(double))
             strcpy(fmt," %8g");
        else strcpy(fmt," %4d");

        if (istr[0]) printf("%s:",istr);

        for (size_t i=0; i<len; i++) printf(fmt,data[i]);
        printf("%s\n", isref ? " ***ISREF***" : "");

        return;
    }

    a=mxCreateDoubleMatrix(1,len,mxREAL);

    double *d=mxGetDoubles(a);
    for (size_t i=0; i<len; i++) { d[i]=(double)data[i]; }

    if (!mflag) {
       if (istr[0])
            { wb_printf("\n%s = [1x%ld double]\n\n", istr, len); }
       else { wb_printf("\n"); }
    }
    else { 
       wb_printf("\n%s = [%s\n",
          istr[0]? istr:"ans", isref? " *isref*":"");
    }

    Wb::CallMatlab(0,NULL,1,&a,"disp");

    if (mflag>1) wb_printf("];\n");

    mxDestroyArray(a);
};

template <class T> inline
wbindex& wbvector<T>::find(const T &x, wbindex &I, char iflag) const {

   if (!len) { I.init(); }
   else {
      size_t i=0, l=0;
      for (; i<len; ++i) { if (data[i]==x) ++l; }

      if (!iflag) { I.init(l); if (I.len) { l=i=0;
         for (; i<len; ++i) { if (data[i]==x) I[l++]=i; }}
      }
      else { I.init(len-l); if (I.len) { l=i=0;
         for (; i<len; ++i) { if (data[i]!=x) I[l++]=i; }}
      }
   }
   return I;
};

template <class T> inline
wbindex& wbvector<T>::find(const T &x, wbindex &I, wbindex &Ix) const {

   if (!len) { I.init(); Ix.init(); }
   else {
      size_t i=0, l=0, k=0;
      for (; i<len; ++i) { if (data[i]==x) ++l; }
      I.init(l); Ix.init(len-l); l=i=0;
      for (; i<len; ++i) { if (data[i]==x) I[l++]=i; else Ix[k++]=i; }
   }
   return I;
};

template <class T> inline
wbindex& wbvector<T>::findGT(const T &x, wbindex &I, char iflag) const {

   if (!len) { I.init(); }
   else {
      size_t i=0, l=0;
      for (; i<len; ++i) { if (data[i]>x) ++l; }

      if (!iflag) { I.init(l); if (I.len) { l=i=0;
         for (; i<len; ++i) { if (data[i]>x) I[l++]=i; }}
      }
      else { I.init(len-l); if (I.len) { l=i=0;
         for (; i<len; ++i) { if (!(data[i]>x)) I[l++]=i; }}
      }
   }
   return I;
};

template <class T> inline
wbindex& wbvector<T>::findGT(const T &x, wbindex &I, wbindex &Ix) const {

   if (!len) { I.init(); Ix.init(); }
   else {
      size_t i=0, l=0, k=0;
      for (; i<len; ++i) { if (data[i]>x) ++l; };
      I.init(l); Ix.init(len-l); l=i=0;
      for (; i<len; ++i) { if (data[i]>x) I[l++]=i; else Ix[k++]=i; }
   }
   return I;
};

template <class T> inline
size_t wbvector<T>::numZeros(T eps) const {
   size_t n=0, i=0;
   if (!eps ) { 
      for (; i<len; ++i) { if (!data[i]) { ++n; }}
   }
   else if (eps>0) {
      for (; i<len; ++i) { if (Wb::abs(data[i])<eps) { ++n; }}
   }
   else {
       wblog(FL,"WRN %s() got eps=%.3g (-> 0) !?",FCT,double(eps));
       for (; i<len; ++i) { if (!data[i]) { ++n; }} 
   }
   return n;
};

template <class T> inline 
size_t wbvector<T>::numEQ(const T &x) const {
   size_t n=0, i=0;
   for (; i<len; ++i) { if (data[i]==x) ++n; }
   return n;
};

template <class T> inline 
size_t wbvector<T>::numGT(const T &x) const {
   size_t n=0, i=0;
   for (; i<len; ++i) { if (data[i]>x) ++n; }
   return n;
};

template <class T> inline 
size_t wbvector<T>::numLT(const T &x) const {
   size_t n=0, i=0;
   for (; i<len; ++i) { if (data[i]<x) ++n; }
   return n;
};

template <class T> inline
wbindex& wbvector<T>::findRange(
   const T &x1, const T &x2, wbindex &I, const char *w
) const {

   if (len==0) { I.init(); return I; }

   bool i1=1,i2=1; size_t k=0;
   if (w) {
      if ((w[0]!='[' && w[0]!=']') || (w[1]!='[' && w[1]!=']') || w[2])
      wblog(FL,"ERR invalid interval specification '%s'",w);
      i1=(w[0]=='[');
      i2=(w[1]==']');
   }

   I.init(len);
   for (size_t i=0; i<len; i++) { const T &x=data[i];
      if (x<x1 || x>x2 || (x==x1 && !i1) || (x==x2 && !i2)) continue;
      I[k++]=i;
   }
   if (k) I.len=k; else I.init();

   return I;
};

template <class T>
size_t wbvector<T>::findValues(const wbvector<T> &B0, wbindex &Ia) const {

   size_t l=0,i=0,j=0, m=len, n=B0.len;
   char c;

   wbvector<T> A(*this), B(B0);
   wbperm P1,P2; wbindex I;
   T *a, *b;

   if (A.isEmpty() || B.isEmpty()) { Ia.init(); return 0; }

   A.Sort(P1);
   B.Sort(P2); I.init(m); a=A.data; b=B.data;

   while (i<m && j<n) { c=NUMCMP(a[i],b[j]);
      if (c<0) i++; else
      if (c>0) j++;
      else {
         do { I[l++]=i++; } while (i<m && a[i-1]==a[i]);
         do {        j++; } while (j<n && b[j-1]==b[j]);
      }
   }

   if (l) I.len=l; else I.init();

   P1.get(I,Ia); 

   return Ia.len;
}

template <class T> inline
wbvector<T>& wbvector<T>::getI(size_t k, wbvector<T> &v) const {

    if (k>=len) wblog(FL,
       "ERR Index out of bounds (%d/%d) !?",k,len);

    v.init(len-1);
    for (size_t l=0, i=0; i<len; i++) if (i!=k) v.data[l++]=data[i];
    return v;
};

template <class T> inline
wbvector<T>& wbvector<T>::getI(
  const wbindex &I, wbvector<T> &v, char uflag
) const {

    wbindex J; I.invert(len,J,uflag);
    return select(J,v);
};

template <class T> inline
wbvector<T>& wbvector<T>::initI(
    const wbvector<T> &a, const wbindex &ia,
    const wbvector<T> &b, const wbindex &ib, char uflag
){
    if (!a.len && !b.len) {
       if (ia.len || ib.len) wblog(FL,
          "ERR %s() index out of bounds (%d+%d)",FCT,ia.len,ib.len);
       return init();
    }

    size_t i=0, j=0, na=a.len, nb=b.len, n=na+nb;
    char ma[n], *mb=ma+na; memset(ma,0,n*sizeof(char));

    if (uflag) {
       for (; i<ia.len; ++i) { j=ia[i];
          if (j>=na) wblog(FL,
             "ERR %s() index out of bounds (%d/%d)",FCT,j+1,a.len);
          if ((++ma[j])>1) wblog(FL,
             "ERR %s() index not unique (%d/%d)",FCT,j+1,a.len
          );
       }
       for (i=0; i<ib.len; ++i) { j=ib[i];
          if (j>=nb) wblog(FL,
             "ERR %s() index out of bounds (%d/%d)",FCT,j+1,b.len);
          if ((++mb[j])>1) wblog(FL,
             "ERR %s() index not unique (%d/%d)",FCT,j+1,b.len
          );
       }
       na-=ia.len; nb-=ib.len;
    }
    else {
       for (; i<ia.len; ++i) { j=ia[i];
          if (j>=a.len) wblog(FL,
             "ERR %s() index out of bounds (%d/%d)",FCT,j+1,a.len);
          if ((++ma[j])==1) --na; else ma[j]=1;
       }
       for (i=0; i<ib.len; ++i) { j=ib[i];
          if (j>=nb) wblog(FL,
             "ERR %s() index out of bounds (%d/%d)",FCT,j+1,b.len);
          if ((++mb[j])==1) --nb; else mb[j]=1;
       }
    }

    init(na+nb);

    if (len) { j=0;
       if (na) {
          for (i=0; i<a.len; ++i) { if (!ma[i]) { data[j]=a[i]; ++j; } }
       }
       if (nb) {
          for (i=0; i<b.len; ++i) { if (!mb[i]) { data[j]=b[i]; ++j; } }
       }
       if (j!=len) wblog(FL,"ERR %s() %d/%d !?",FCT,j,len); 
    }

    return *this;
};

template <class T>
template <class T2> 
wbvector<T>& wbvector<T>::initT(
   const char *F, int L, const wbvector<T2> &v) {

   RENEW(v.len);
   for (size_t i=0; i<len; ++i) {
      data[i]=T(v.data[i]);
      if (T2(data[i])!=v.data[i]) wblog(F_L,
         "ERR %s() got rounded value: %g (%s) -> %g (%s)",
         FCT, double(data[i]), TSTR(T2),
         double(v.data[i]), TSTR(T)
      );
   }
   return *this;
};

template <class T>
template <class TD> inline
wbvector<T>& wbvector<T>::initT(const wbsparray<TD> &S) {

   if (!S.isVector(FL)) wblog(FL,
      "ERR %s() invalid sparse vector (%s)",FCT,S.sizeStr().data);

   size_t i=0, j=0, n=S.IDX.dim1, m=S.IDX.dim2;
   const SPIDX_T *idx=S.IDX.data;

   for (; i<S.SIZE.len; ++i) { if (S.SIZE[i]>1) { idx+=i; break; }}

   RENEW(S.numel());

   for (i=0; i<n; ++i, idx+=m) { j=(*idx);
      if (j>=len) wblog(FL,
         "ERR %s() index out of bounds (%d/%d)",FCT,j+1,len);
      data[j]=T(S.D.data[i]);
   }

   return *this;
};

bool mxIsWbvector(
    const char *F, int L, const mxArray *a,
    size_t *n_, const char *istr, char dflag
){
    size_t n=0;

    if (a && (n=mxGetNumberOfElements(a))) { int i=0;
       if (istr && istr[0])
            i=Mx::IsDblVector(0,0,a,0,dflag);
       else i=Mx::IsDblVector(F,L,a,0,dflag);

       if (!i) {
          if (F) wblog(F,L,
            "ERR %s() %s%snumeric vector required\n%s   (%s; %d)", FCT,
             istr?istr:"", istr?" ":"", str, mxGetClassName(a), n);
          return 0;
       }
    }

    if (n_) { (*n_)=n; } 
    return n; 
};

template <class T>
int wbvector<T>::init_base( 
    const char *F, int L,
    const mxArray *a, const char *istr, char tcheck, char ref
){
    size_t n=0;

    if (!a || !mxIsWbvector(F_L,a,&n,istr,'*')) {
       init(); return 1; 
    }

    Mx::Array<T> A(a);
    if (ref && A.data)
         { init2ref(n,A.data); }
    else { init(n); A.copyTo(data,tcheck); }

    if (data && !Wb::is_finite(data,1)) wblog(PFL,
       "WRN %s() encountered nan or inf in mex input data",FCT);

    return 0;
};

#ifdef QS_USING_MPFR

template<>
int wbvector<Wb::quad>::init_mpfr( 
   const char *F, int L, const mxArray *a, char base,
   const char *istr, char tcheck
){
   if (!a) { init();
      if (F) wblog(FL,"ERR %s() got null mxArray",FCT);
      return -1;
   }
   if (!mxIsNumChar(a)) { init();
      if (F) wblog(FL,
         "ERR %s() invalid type `%s'",FCT,mxGetClassName(a));
      return -2;
   }
   if (mxIsEmpty(a)) { init(); return 0; }

   if (base<=0 && mxIsWbvector(0,0,a) && mxIsDouble(a)) {
	  wbvector<double> X; X.init(F_L,a,istr,tcheck);
	  this->initT(X); return 0;
   }

   if (base<=0) base=30; else 
   if (base==1 || base>62) wblog(FL,
      "ERR %s() invalid base=%d<%c>\n"
      "(interpreted as base for Wb::quad)",FCT,base,base);

   unsigned n=0, np=1; int e=0;

   unsigned r=mxGetNumberOfDimensions(a);
   const size_t *sz=mxGetDimensions(a);

   wbvector<char> ss; char *sd=NULL;

   if (r!=2) {
      if (F) wblog(FL,"ERR invalid rank-%d array",r);
      init(); return -3;
   }
   if (!sz[0] || sz[0]>62 || (sz[0]==1 && (!sz[1] || sz[1]>62))) {
      if (F) wblog(FL,"ERR invalid %ldx%ld array",sz[0],sz[1]);
      init(); return -4;
   }

   if (mxIsInt8(a)) {
      Mx::Array<char> A(FL,a); 
      sd=A.data;

      if (sz[0]>1)
           { n=sz[0]; init(sz[1]); } 
      else { n=sz[1]; init(sz[0]); } 
   }
   else if (mxIsChar(a) && (sz[0]==1 || sz[1]==1)) {
      Mx::Array<mxChar> A(FL,a);

      if (sz[0]>1)
           { n=sz[0]; init(sz[1]); }
      else { n=sz[1]; init(sz[0]); }  

      ss.initT(n,A.data); 
      sd=ss.data;
   }
   else { 
      MXPut(FL,"Idbg")
         .addP(mxDuplicateArray(a),"a").add(wbstring(istr),"istr");
      wblog(FL,"ERR %s() invalid input (%s: %d x %d)",
      FCT,mxGetClassName(a),sz[0],sz[1]);
   }

   if (len>=64) {
      np=MAX( QSP_NUM_THREADS, OMP_NUM_THREADS );
   }

  #pragma omp parallel for num_threads(np)
   for (size_t k=0; k<len; ++k) { if (!e) {

      const char *sk=sd+k*n; char cx=0;
      unsigned i=0, i0=0;

      try {
         while (i<n && sk[i]==' ') { ++i; }; i0=i;
         while (i<n && sk[i]!=' ' && sk[i]) { ++i; }

         if (i0<i && i<n) {
            cx=sk[i];                       ((char*)sk)[i]=0;
            data[k].init_s(F_L,sk+i0,base); ((char*)sk)[i]=cx;
         }
         else { char s[n+1]; ++e;  strncpy(s,sk,n); s[n]=0;
            wblog(FL,"... error i=%4ld/%ld: `%s' (%d..%d/%d)",
            k,len,s,i0,i,n);
         }
      }
      catch (...) {
         ++e; if (cx) { ((char*)sk)[i]=cx; }
      }
   }}
   if (e) wblog(FL,"ERR %s() e=%d",FCT,e);

   return 0;
};

template<>
int wbvector<Wb::quad>::init(
   const char *F, int L,
   const mxArray *a, const char *istr, char tcheck,
   char ref QS_UNUSED_VAR
){
   unsigned base=(tcheck<2 ? 30 : tcheck);

   if (base>30) wblog(FL,
      "WRN %s() using tcheck=%d as mpfr base !?",FCT,tcheck);

   return init_mpfr(F,L,a,base,istr,tcheck);
};

#endif

template <class T>
int wbvector<T>::init_Struct( 
    const char *F, int L, const mxArray *a, const char *istr, char ref
){
    if (!a || !mxIsStruct(a)) { if (F) wblog(F,L,
       "ERR %s() invalid input %s",FCT,istr?istr:"");
       init(); return 1;
    }

    unsigned i=0, n=mxGetNumberOfElements(a);
    init(n);

    for (; i<n; ++i) { data[i].init(F_L,a,ref,i); }

    return 0;
};

template <class T> inline
void wbvector<T>::Sort(char dir) {
   wbperm p; Wb::hpsort(*this,p,dir);
};

template <class T> inline
void wbvector<T>::sort(wbvector &v2) const {
   wbperm p; v2=(*this); Wb::hpsort(v2,p);
};

template <class T> inline
wbperm& wbvector<T>::pSort(wbperm &p) const { 
   if (!isSorted()) { wbvector x(*this); Wb::hpsort(x,p); }
   else p.init();
   return p;
};

template <class T>
inline wbvector<T>& wbvector<T>::Permute(const wbperm &P, char iflag) {
   const wperm_t *p=P.data; char q=0;

   if (!P.len || (P.len<=len && (q=P.isValidPerm())>1)) { return *this; }

   if (P.len!=len) {
      wbperm P_; P_.init(P,iflag,len);
      return Permute(P_);
   }

   if (q) {
      size_t i=0; wbvector<T> v(*this);
      if (!iflag)
           { for (; i<len; ++i) data[i]=v.data[p[i]]; }
      else { for (; i<len; ++i) data[p[i]]=v.data[i]; }
   }
   else if (P.len<12) wblog(FL,
       "ERR invalid permutation [%s] (len=%d/%d)",STR(P),P.len,len);
   else wblog(FL,"ERR invalid permutation (len=%d/%d)", P.len, len);

   return *this;
};

template <class T>
inline wbvector<T>& wbvector<T>::permute(
   wbvector<T> &v, const wbperm &P, char iflag
 ) const {

   if (P.isEmpty()) { v=*this; return v; }
   const wperm_t *p=P.data;

   if (len!=P.len || !validPerm(P)) {
      if (P.len<12) wblog(FL,
         "ERR invalid permutation [%s] (%d/%d)",STR(P),len,P.len);
      else wblog(FL,"ERR invalid permutation (%d/%d)",len,P.len);
   }

   if (!P.isIdentityPerm()) { v.init(len);
      if (iflag==0) {
             for (size_t i=0; i<len; i++) v.data[i]=data[p[i]]; }
      else { for (size_t i=0; i<len; i++) v.data[p[i]]=data[i]; }
   }
   else v=(*this);

   return v;
};

template <class T> inline
wbvector<T>& wbvector<T>::blockPermute(
   const wbperm &P, wbvector<T> &v, char iflag
) const {

   if (isEmpty() || P.isEmpty() || P.isIdentityPerm()) {
       if (P.len && len%P.len) wblog(FL,
          "ERR %s() data mismatch (%d mod %d !?)",FCT,len,P.len);
       v=*this; return v;
   }
   if (&v==this) {
      wbvector<T> x(*this);
      return x.blockPermute(P,v,iflag);
   }

   if (len%P.len || !validPerm(P)) wblog(FL,
      "ERR invalid permutation (%d/%d)", len, P.len);
   v.init(len);

   size_t i,j, m=len/P.len;
   const wperm_t *p=P.data; const T *d0=data; T *d=v.data;

   if (iflag==0) {
      for (i=0; i<P.len; ++i, d+=m) { d0=data+m*p[i];
      for (j=0; j<m; ++j) d[j]=d0[j]; }
   }
   else {
      for (i=0; i<P.len; ++i, d0+=m) { d=v.data+m*p[i];
      for (j=0; j<m; ++j) d[j]=d0[j]; }
   }

   return v;
};

template <class T> inline
int wbvector<T>::blockCompare(
   const wbperm &Pa, const wbvector &B, const wbperm &Pb
 ) const {

   if (len!=B.len ||
      (Pa.len &&   len%Pa.len) || (Pb.len && B.len%Pb.len))
      wblog(FL,"ERR %s() got length mismatch (%d@%d, %d@%d)",
      FCT,len,Pa.len,B.len,Pb.len
   );

   if ((!Pa.len && !Pb.len) || Pa.len==len || Pb.len==B.len) {
      if (Pa.len) {
         if (Pb.len) {
            if (Pa.len!=Pb.len) wblog(FL,
               "ERR %s() got length mismatch (%d@%d, %d@%d)",
               FCT,len,Pa.len,B.len,Pb.len
            );

            const wperm_t*pa=Pa.data, *pb=Pb.data;
            for (size_t i=0; i<len; ++i) {
               if (data[pa[i]]!=B.data[pb[i]]) {
                  return (data[pa[i]]<B.data[pb[i]] ? -1 : +1);
               }
            }
         }
         else {
            const wperm_t *pa=Pa.data;
            for (size_t i=0; i<len; ++i) {
               if (data[pa[i]]!=B.data[i]) {
                  return (data[pa[i]]<B.data[i] ? -1 : +1);
               }
            }
         }
      }
      else {
         if (Pb.len) {
            const wperm_t *pb=Pb.data;
            for (size_t i=0; i<len; ++i) {
               if (data[i]!=B.data[pb[i]]) {
                  return (data[i]<B.data[pb[i]] ? -1 : +1);
               }
            }
         }
         else {
            for (size_t i=0; i<len; ++i) {
               if (data[i]!=B.data[i]) {
                  return (data[i]<B.data[i] ? -1 : +1);
               }
            }
         }
      }
   }
   else {
      if (Pa.len) {
         if (Pb.len) {
            const wperm_t *pa=Pa.data, *pb=Pb.data;
            size_t j, i=0, m=len/Pa.len; 
            const T *a, *b;

            if (Pa.len!=Pb.len) wblog(FL,
               "ERR %s() got invalid permutations (%d@%d, %d@%d)",
               FCT,len,Pa.len,B.len,Pb.len
            );

            for (; i<Pa.len; ++i) { a=data+m*pa[i]; b=B.data+m*pb[i];
               for (j=0; j<m; ++j) {
                  if (a[j]!=b[j]) return (a[j]<b[j] ? -1 : +1);
               }
            }
         }
         else {
            const wperm_t *pa=Pa.data;
            size_t j, i=0, m=len/Pa.len;
            const T *a, *b=B.data;

            for (; i<Pa.len; ++i, b+=m) { a=data+m*pa[i];
               for (j=0; j<m; ++j) {
                  if (a[j]!=b[j]) { return (a[j]<b[j] ? -1 : +1); }
               }
            }
         }
      }
      else {
         if (Pb.len) {
            const wperm_t *pb=Pb.data;
            size_t j, i=0, m=B.len/Pb.len;
            const T *a=data, *b;

            for (; i<Pb.len; ++i, a+=m) { b=B.data+m*pb[i];
               for (j=0; j<m; ++j) {
                  if (a[j]!=b[j]) { return (a[j]<b[j] ? -1 : +1); }
               }
            }
         }
      }
   }

   return 0;
};

template <class T>
inline int wbvector<T>::set2Group(
   const wbperm &P, const WBINDEX &D,
   const T* S0,
   char iflag
){
   size_t i,j,l,d,n=D.len; int e=0;
   wbvector<T> X;

   if (!S0 && len!=P.len) wblog(FL,
   "ERR %s() severe size mismatch (%d/%d)", FCT, len, P.len);

   if (n==0) {
      if (P.len) wblog(FL,
      "ERR empty object inconsistency (%d,%d)", P.len, D.len);
      return 0;
   }

   if (S0==NULL || S0==data) { save2(X); S0=X.data; }

   init(n);

   for (l=i=0; i<n; i++, l+=d) { d=D[i];
      data[i]=S0[P[l]];

      for (j=1; j<d; j++) if (S0[P[l+j]]!=data[i]) {
         e++; if (iflag)
         sprintf_str("%g/%g", (double)data[i], (double)S0[P[l+j]]);
      }
   }

   return e;
}

template<class T>
mxArray* wbvector<T>::mxCreateStruct(unsigned m, unsigned n) const {
   return mxCreateCellMatrix(m,n);
}

template<class T> 
void wbvector<T>::add2MxStruct(mxArray *S, unsigned i, char tflag) const {
   mxSetCell(S,i,toMx(tflag));
}

template<class T>
mxArray* wbvector<T>::mxCreateCell(unsigned m, unsigned n) const {
   return mxCreateCellMatrix(m,n);
}

template<class T>
void wbvector<T>::add2MxCell(mxArray *S, unsigned i, char tflag) const {
   mxSetCell(S,i,toMx(tflag));
}

template<class T> inline
void wbvector<T>::add2MxStruct(
    mxArray *S, const char *vname, char tflag) const { 

    int fid; mxArray *a;

    fid = mxAddField2Scalar(FL,S,vname);

    if (!tflag)
         a=mxCreateDoubleMatrix(1,len,mxREAL);
    else a=mxCreateDoubleMatrix(len,1,mxREAL);

    double *d=mxGetDoubles(a);
    for (size_t i=0; i<len; ++i) d[i]=(double)data[i];

    mxSetFieldByNumber(S,0,fid, a);
};

template<class T> inline
wbstring wbvector<T>::sizeStr() const {
   wbstring s; 
   if (len) { s=toStrf("","x"); } else { s="[]"; }
   return s;
};

template<class T> inline
wbstring wbvector<T>::toStr() const { return toStr(-1," "); }

template<class T> inline
wbstring wbvector<T>::toStr(int n, const char *sep) const {

    wbstring s(MAX(size_t(32), len*MAX(16,int(n)+6))), fmt;
    fmt.init2Fmt((T)0,n);

    for (size_t i=0; i<len; ++i) {
        if (i>0) { s.push(FL,sep); }
        s.pushf(FL,fmt.data,data[i]);
    }
    return s;
};

template<> inline
wbstring wbvector<wbcomplex>::toStr(int n, const char *sep) const {

    wbstring s(MAX((size_t)32, len*MAX(16,abs(n)+6))), fmz, fmd;
    fmz.init2Fmt(wbcomplex(0),n); fmd.init2Fmt(double(0),n);

    for (size_t i=0; i<len; ++i) { if (i) { s.push(FL,sep); }
        if (data[i].r && data[i].i)
           s.pushf(FL,fmz.data,data[i].r,data[i].i);
        else if (data[i].r)
             { s.pushf(FL,fmd.data,data[i].r); }
        else { s.pushf(FL,fmd.data,data[i].i).push("i"); }
    }
    return s;
};

#ifdef QS_USING_MPFR
template<> inline
wbstring wbvector<Wb::quad>::toStr(int n, const char *sep) const {

    wbstring s(MAX((size_t)32, len*MAX(16,int(n)+6))), fmt;
    fmt.init2Fmt((double)0,n);

    for (size_t i=0; i<len; ++i) {
        if (i>0) { s.push(FL,sep); }
        s.pushf(FL,fmt.data,double(data[i]));
    }

    return s;
};
#endif

template<class T> inline
wbstring wbvector<T>::toStrf (
    const char *fmt_,  const char *sep, unsigned stride, const char *sep2
  ) const {

    wbstring s(MAX(32U, unsigned(len)*16)), fmt;  
    unsigned i_,isep=0;

    if (fmt_ && fmt_[0])
         { fmt=fmt_; }
    else { fmt.init2Fmt(T(0)); }

    for (unsigned i=0; i<len; ++i) { if (++isep>1) {
       s.push (FL,sep); }
       s.pushf(FL,fmt.data,data[i]); if (stride) { i_=i+1;
       if (!(i_%stride) && i_<len) { s.push(FL,sep2); isep=0; }}
    }

    return s;
};

template<> inline
wbstring wbvector<wbcomplex>::toStrf(
   const char *fmt_, const char *sep, unsigned stride, const char *sep2
 ) const {

   wbstring s; 
   unsigned w,i=0; char fmt[16], s1[32], flag=1;

   snprintf(fmt,16,"%s",fmt_ && fmt_[0] ? fmt_ : "%.4g");
   w=2 + sprintf_str(fmt, sqrt(2.))
       + sprintf_str(fmt,-sqrt(2.)*1e-20);

   s.init(32+(4+i)*len); 

   for (size_t i=0; i<len; ++i) {
      if (i>0) { if (flag) { s.push(FL,sep); } else { flag=1; }}

      snprintf(s1,32,"%s",data[i].toStr(fmt).data);
      s.pushf(FL,"%*s",w,s1);
      if (stride && ((i+1)%stride)==0 && i+1<len) {
         s.push(FL,sep2); flag=0;
      }
   }

   return s;
};

#ifdef QS_USING_MPFR

template<> inline
wbstring wbvector<Wb::quad>::toStrf(
   const char *fmt_, const char *sep, unsigned stride, const char *sep2
 ) const {

   wbvector<double> x; x.init(*this);
   return x.toStrf(fmt_,sep,stride,sep2);
};

#endif

#ifdef __WB_MPFR_HH__

template <> inline
mxArray* wbvector< Wb::quad >::toMx(char base) const { 

   if (base<=0) base=30; else 
   if (base==1 || base>62) wblog(FL,
      "ERR %s() invalid base=%d<%c>\n"
      "(interpreted as base for Wb::quad)",FCT,base,base);

   unsigned nmax=0, n=0; size_t S[2];
   mxArray *a;

   wbstring s_; char *s=NULL;
   if (len) { data[0].toStr(s_,base); n=s_.len; s=s_.data; }

   if (len!=1)
        { S[0]=n; S[1]=len; } 
   else { S[0]=len; S[1]=n; } 

   if (len==1 && data[0]==Wb::quad(int(double(data[0])))) {
      unsigned i=0, l=0; 

      while (i<n && s[i] && s[i]!='.') { ++i; }
      if (i<n && s[i]) { l=i;
         for (++i; i<n && s[i]; ++i) { if (s[i]!='0')
            wblog(FL,"ERR %s() expecting integer (%s) !?",FCT,s);
         }
         s[l]=0; S[1]=l+1; 

         return a=Mx::Array<mxChar>(2,S).copyFromTR(s);
      }
   }

   char *cdat=NULL;
   a=Mx::Array<char>(2,S).Return(cdat); 

   if (!len) { return a; }

   unsigned np=1; int e=0;
   if (len>=64) {
      np=MAX( QSP_NUM_THREADS, OMP_NUM_THREADS );
   }

  #pragma omp parallel for num_threads(np) 
   for (size_t k=0; k<len; ++k) { if (!e) {
      unsigned i=0, i0=0;
      char *sk=cdat+k*n;

      try { data[k].toStr(sk,n,base,-1,'a'); } 
      catch (...) { ++e; } 

      while (i<n && s[i]==' ') { ++i; }; i0=i;
      while (i<n && s[i]!=' ' && s[i]) { ++i; }

      if (i<2 || i0>=i || s[i] || (++i)>=n) { ++e;
         char s[n+1]; strncpy(s,sk,n); s[n]=0;
         wblog(FL,"... error i=%4ld/%ld: `%s' (%d..%d/%d)",
         k,len,s,i0,i,n);
      }

     #pragma omp critical (mpfr_collecting_nmax)
      { if (nmax<i) nmax=i; }
   }}

   if (e) wblog(FL,"ERR %s() e=%d",FCT,e);

   if (nmax==n) wblog(FL,
      "WRN %s<Wb::quad> got n=%d/%d (len=%ld)",FCT,nmax,n,len);
   if (nmax>n) wblog(FL,
      "ERR %s<Wb::quad> got n=%d/%d (len=%ld)",FCT,nmax,n,len);

   return a;
};

#endif 

template <class T>
inline wbvector<T>& wbvector<T>::Move2_FE(
    const WBINDEX &I,
    char iflag
){
    size_t i,k;
    wbvector<char> mark(len);   
    wbvector<T> x;

    if (I.isEmpty()) return *this;
    if (!isUniqueIdxSet(I,len)) {
        wblog(FL, "ERR Move2_FE Index out of bounds [%s; %d]",
        I.toStr().data, len); return *this;
    }

    for (i=0; i<I.len; i++) mark[I[i]]++;

    select(I,x);

    if (iflag=='E') {
        for (k=i=0; i<len; i++)  
        if (mark[i]==0) {
            if (k!=i) data[k]=data[i];
            k++;
        }
        memcpy(data+k, x.data, x.len*sizeof(T));  
    }
    else
    if (iflag=='F') { 
        for (k=i=len-1; int(i)>=0; i--) 
        if (mark[i]==0) {
            if (k!=i) data[k]=data[i];
            k--;
        }
        memcpy(data, x.data, x.len*sizeof(T));  
    }
    else wblog(FL,"ERR Move2_FE - Invalid flag %c<%d>", iflag, iflag);

    return *this;
}

template <class T>
inline wbvector<T>& wbvector<T>::Skip(size_t i1) {

    if (i1>=len) {
        wblog(FL, "ERR Skip - index out of range (%d,%d)\n[%s; %d]",
        i1, len, toStr().data, i1); return *this;
    }

    for (size_t i=i1+1; i<len; i++) data[i-1]=data[i];
    Resize(len-1);

    return *this;
}

template <class T>
inline wbvector<T>& wbvector<T>::Skip(
    const WBINDEX &I,
    char rflag 
){
    size_t i,k;
    wbvector<char> mark(len); char *m=mark.data;

    for (k=0; k<I.len; k++) {
        if (I[k]<len) m[I[k]]++;
        else {
            wblog(FL, "ERR Skip - index out of range (%d,%d)\n[%s; %d]",
            I[k]+1, len, toStr().data, k+1); return *this;
        }
    }

    for (k=i=0; i<len; ++i) {
       if (m[i]==0) { if (k!=i) { data[k]=data[i]; }
          ++k;
       }
    }

    if (!rflag) Resize(k);
    else {
       if (k) len=k; else init();
    }

    return *this;
};

#endif

