/* ---------------------------------------------------------------------
 * Project : QSpace tensor library (v4.0)
 * Class   : QSpace memory routines and smart pointers
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

#ifndef __WB_MEM_LIB_HH__
#define __WB_MEM_LIB_HH__

// ------------------------------------------------------------------ //
// ------------------------------------------------------------------ //

template <class T1, class T2>
bool gotMemOverlap(
   const T1* src, size_t len0,
   const T2* dest, size_t len2=-1 
);

template <class T>
inline void MEM_CPY_BASE(T* data, size_t n, T const* const d=NULL) {

#ifdef MLB_VFLAG
    wblog(FL,"TST base init");
#endif
    if (n) {
       if (!d)
            { memset((void*)data,0,n*sizeof(T)); }
       else { memcpy((void*)data,d,n*sizeof(T)); }
    }
};

template <class T>
inline void MEM_CPY_BASE(
    T* data, size_t n, size_t n1, T const* const d1, T const* const d2
){

#ifdef MLB_VFLAG
    wblog(FL,"TST base init %d->%d",n1,n);
#endif
    if (n1) {
       if (d1==NULL) wblog(FL,"ERR got (null) pointer"); else
       if (n1>n) { wblog(FL,"WRN length out of bounds (%d/%d)",n1,n); n1=n; }
       memcpy(data,d1,n1*sizeof(T));
    }
    if (n1<n) {
       if (d2) memcpy(        data+n1, d2, (n-n1)*sizeof(T));
       else    memset((void*)(data+n1), 0, (n-n1)*sizeof(T));
    }
};

namespace Wb {

enum MTYPE {
   MEM_DEF,  
   MEM_REF,      
   MEX_RETURN,   
   MEM_NUM_TYPES 
};

const char* MTYPE_STR[MEM_NUM_TYPES] {
   "MEM:default",
   "MEM:ref",
   "MEX:return"
};

template <class T>
class sptr { 

 public:

   sptr() : mtype(MEM_DEF), len(0), nref(0), data(NULL) {};

  ~sptr() {
      #pragma omp critical (manage_SPTR) 
      {  if (mtype || len || data || nref) wblog(FL,
            "WRN %s() got %s",FCT,STR_(this));
         init_def();
      }
   }

   void init_def() {
      mtype=MEM_DEF; len=nref=0; data=NULL;
   }

   T* malloc_base(size_t n, char init=1, MTYPE mt=MEM_DEF);
   T* malloc_data(size_t n, const T *d0=NULL, MTYPE mt=MEM_DEF);

   int rm_dref(char mflag=0); 

   unsigned add_dref();

   T* save_dref(const char *F=0, int L=0); 

   T* init2ref(size_t n, const T* d0);

   int isRef() const {
      return (nref || mtype==Wb::MEM_REF);
   };

   explicit operator bool() const {
      check_consistency(); return (len!=0); };
   bool operator! () const {
      check_consistency(); return (len==0); }

   int check_consistency(char force=1) const {
      int e=0; 

      #pragma omp critical (manage_SPTR) 
      { if (!len ^ !data) { e|=1; }      
        if (nref && !len) { e|=2; }
      }

      if (e && force) wblog(FL,
         "ERR Wb::sptr() inconsistency (%p @ n=%d, r=%d)",data,len,nref);
      return e;
   };

   wbstring toStr() const;

   MTYPE mtype; 

   size_t len;
   unsigned nref;

   T* data;

 protected:
 private:
};

}; 

template <class T, ENABLE_IF_isPOD(T) >
void MEM_CPY(T* data, size_t n, const T* d=NULL) {
   MEM_CPY_BASE(data,n,d);
};

template <class T, ENABLE_IF_noPOD(T) >
void MEM_CPY(T* data, size_t n, const T* d=NULL) {

  #ifdef MLB_VFLAG
   wblog(FL,"TST class init");
  #endif

   if (d) { for (size_t i=0; i<n; ++i) data[i]=d[i]; }
   else {
      const T z=T(); 
      for (size_t i=0; i<n; ++i) data[i]=z;
   }
};

template <class T, ENABLE_IF_isPOD(T) >
void MEM_CPY(T* data, size_t n, size_t n1, const T* d1, const T* d2) {
   MEM_CPY_BASE(data,n,n1,d1,d2);
};

template <class T, ENABLE_IF_noPOD(T) >
void MEM_CPY(T* data, size_t n, size_t n1, const T* d1, const T* d2) {
  #ifdef MLB_VFLAG
   wblog(FL,"TST class init %d->%d",n1,n);
  #endif

   size_t i=0;
   if (n1) {
      if (d1==NULL) wblog(FL,"ERR got (null) pointer"); else
      if (n1>n) { wblog(FL,
         "WRN length out of bounds (%d->%d)",n1,n); n1=n; }
      for (; i<n1; ++i) data[i]=d1[i];
   }
   if (d2) { for (; i<n; ++i) { data[i]=d2[i-n1]; }}
   else {
      const T z=T(); 
      for (; i<n; ++i) data[i]=z;
   }
};

template <class T>
void MEM_CPY(T** data, size_t n, T* const* const d=NULL) {
  #ifdef MLB_VFLAG
   wblog(FL,"TST got pointer space");
  #endif
   MEM_CPY_BASE(data,n,d);
};

template <class T>
void MEM_CPY(T** data, size_t n, size_t n1,
   T* const* const d1, T* const* const d2
){
  #ifdef MLB_VFLAG
   wblog(FL,"TST got pointer space");
  #endif
   MEM_CPY_BASE(data,n,n1,d1,d2);
};

template <class T> inline
void MEM_SET_BASE(T* data, size_t len, const T *a=NULL) {
#ifdef MLB_VFLAG
   wblog(FL,"TST base memset");
#endif
   if (!a)
        { memset((void*)data,0,len*sizeof(T)); }
   else { for (size_t i=0; i<len; ++i) data[i]=(*a); }
};

template <class T> inline  
void MEM_SET_BASE(T** data, size_t len, T* const *a=NULL) {
#ifdef MLB_VFLAG
   wblog(FL,"TST base memset (pointer space)");
#endif
   if (!a)
        { memset((void*)data,0,len*sizeof(T*)); }
   else { for (size_t i=0; i<len; ++i) data[i]=(*a); }
};

template <class T>
class MEM_SET {

  public:

     MEM_SET(
        T* data      QS_UNUSED_VAR,
        size_t len   QS_UNUSED_VAR,
        const T *a   QS_UNUSED_VAR = NULL
     ){ 
     };

  private:

};

template <class T>
class MEM_SET<T*> { 

  public:

     MEM_SET(T** data, size_t len, T* const *a=NULL) {
#ifdef MLB_VFLAG
        wblog(FL,"TST pointer memset");
#endif
        MEM_SET_BASE(data,len,a); 
     }

  private:

};

template <> inline MEM_SET<int>::MEM_SET(
   int* data, size_t len, const int *a
){ MEM_SET_BASE(data,len,a); }

template <> inline MEM_SET<unsigned>::MEM_SET(
   unsigned* data, size_t len, const unsigned *a
){ MEM_SET_BASE(data,len,a); }

template <> inline MEM_SET<unsigned long>::MEM_SET(
   unsigned long* data, size_t len, const unsigned long *a
){ MEM_SET_BASE(data,len,a); }

template <> inline MEM_SET<long>::MEM_SET(
   long* data, size_t len, const long *a
){ MEM_SET_BASE(data,len,a); }

template <> inline MEM_SET<char>::MEM_SET(
   char* data, size_t len, const char *a
){ MEM_SET_BASE(data,len,a); }

template <> inline MEM_SET<double>::MEM_SET(
   double* data, size_t len, const double *a
){ MEM_SET_BASE(data,len,a); }

template <> inline MEM_SET<float>::MEM_SET(
   float* data, size_t len, const float *a
){ MEM_SET_BASE(data,len,a); }

template <> inline MEM_SET<wbcomplex>::MEM_SET(
   wbcomplex* data, size_t len, const wbcomplex *a
){ MEM_SET_BASE(data,len,a); }

#endif

