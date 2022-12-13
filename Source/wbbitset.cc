#ifndef __WB_BITSET_CC__
#define __WB_BITSET_CC__

// ---------------------------------------------------------------- //
   namespace Wb {
// ---------------------------------------------------------------- //

   template <class T>
   wbstring bitstr(const T &x) { // NB! i is 0-based
      Bitset B(8*sizeof(T)); memcpy(B.data,(void*)&x,B.len);
      return B.toStr();
   };

wbstring Bitset::toStr() const {
   wbstring s(nbits+len+1,' '); 
   unsigned i=0, l=0;
      s[l++]='['; s[s.len-2]=']';
   for (unsigned k=0; k<len; ++k) {
      const char &b=data[k]; char r=1;
      for (; r && i<nbits; r<<=1, ++i) { s[l++]=( b&r ? '1' : '-'); }
      if (!r && i<nbits) s[l++]='|';
   }
   if (l>s.len) wblog(FL,
      "ERR %s() index out of bounds (%d/%d)",FCT,l,s.len);

   return s;
};

mxArray* Bitset::toMx() const {
   mxArray *a=mxCreateStruct(1,1);
   add2MxStruct(a,0);
   return a;
};

mxArray* Bitset::mxCreateStruct(unsigned m, unsigned n) const {
   const char *fields[]={"nbits","bset","data"};
   return mxCreateStructMatrix(m,n,3,fields);
};

void Bitset::add2MxStruct(mxArray *S, unsigned i) const {
   mxSetFieldByNumber(S,i,0, numtoMx(nbits));
   mxSetFieldByNumber(S,i,1, toStr().toMx());
   wbvector<unsigned> v(len);
   for (unsigned j=0; j<len; ++j) v.data[j]=(unsigned char)data[j];
   mxSetFieldByNumber(S,i,2, v.toMx());
};

void Bitset::put(const char *file, int line,
   const char *vname, const char *ws
 ) const {

   mxArray *a=toMx();
   int i=mexPutVariable(ws,vname,a);

   if (i) wblog(__FL__,
      "ERR failed to write variable `%s' (%d)",vname,i);
   if (file) wblog(file,line,
      "I/O putting variable '%s' to %s",vname,ws);
   mxDestroyArray(a);
};

template <class T>
wbstring bitset<T>::toStr() const {
   unsigned len=sizeof(T);

   wbstring s(9*len+1,' '); 
   unsigned i=1, l=0; T r=1;

   s[l++]='['; s[s.len-2]=']';
   for (; r; r<<=1, ++i) {
       s[l++]=( data&r ? '1' : '-');
       if (i%8==0 && l+2<s.len) { s[l++]='|'; }
   }

   return s;
};

template <class T>
mxArray* bitset<T>::toMx() const {
   mxArray *a=mxCreateStruct(1,1);
   add2MxStruct(a,0);
   return a;
};

template <class T>
mxArray* bitset<T>::mxCreateStruct(unsigned m, unsigned n) const {
   const char *fields[]={"nbits","bset","data"};
   return mxCreateStructMatrix(m,n,3,fields);
};

template <class T>
void bitset<T>::add2MxStruct(mxArray *S, unsigned i) const {
   mxSetFieldByNumber(S,i,0, numtoMx(8*sizeof(T)));
   mxSetFieldByNumber(S,i,1, toStr().toMx());
   mxSetFieldByNumber(S,i,2, numtoMx(data));
};

template <class T>
void bitset<T>::put(const char *file, int line,
   const char *vname, const char *ws
 ) const {

   mxArray *a=toMx();
   int i=mexPutVariable(ws,vname,a);

   if (i) wblog(__FL__,
      "ERR failed to write variable `%s' (%d)",vname,i);
   if (file) wblog(file,line,
      "I/O putting variable '%s' to %s",vname,ws);
   mxDestroyArray(a);
};

}; 

#endif

