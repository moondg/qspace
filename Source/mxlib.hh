/* ---------------------------------------------------------------------
 * Project : QSpace tensor library (v4.0)
 * Class   : QSpace MEX routines
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

#ifndef __WB_MXLIB_HH__
#define __WB_MXLIB_HH__

/* -------------------------------------------------------------------- */
// switching to matlab 2018a => got interleaved complex format
// => phasing out mxGetPr_17(FL,) and mxGetPi_17(FL,)

#ifdef MATLAB_MEX_FILE

   mxDouble* mxGetPr_17(const char *F, int L, const mxArray *pm) {
      wblog(F_L,"ERR mxGetPr() no longer supported in matlab>=2018");
      return 0;
   };

   mxDouble* mxGetPi_17(const char *F, int L, const mxArray *pm) {
      wblog(F_L,"ERR mxGetPi() no longer supported in matlab>=2018");
      return 0;
   };

#endif

int isHelpIndicator(const char *s) {
   if (s && s[0]=='-') {
      if (!strcmp(s,"-?") || !strcmp(s,"-h") || !strcmp(s,"--help"))
      return 1;
   }
   return 0;
};

#ifdef MATLAB_MEX_FILE

int isHelpIndicator(const mxArray *a) {
   char s[8];
   if (!a || mxGetString(a,s,8)) { return 0; }
   return isHelpIndicator(s);
};

int checkHelpVersion(const mxArray *a0, mxArray **argout) {

   char s[16];
   if (a0 && !mxGetString(a0,s,16)) {
      if (isHelpIndicator(s)) { usage(); return 1; }

      if (!strcmp(s,"--ping")) { return 2; }

      if (!strcmp(s,"--version")) { Wb::VersionInfo I;
         if (argout)
              { argout[0]=I.toMx(); } 
         else { I.print(); }
         return 1;
      }
   }
   return 0;
};

#else 

int checkHelpVersion(const char *s) {
   if (s) {
      if (isHelpIndicator(s)) { usage(); return 1; }
      if (!strcmp(s,"--ping")) { return 2; }

      if (!strcmp(s,"--version")) {
         Wb::VersionInfo().print();
         return 1;
      }
   }
   return 0;
};

#endif

inline mxArray* mx_get_field(mxArray *S, unsigned k, int i) {
   if (S && i>=0 && int(k)>=0) {
      mxArray *a=mxGetFieldByNumber(S,k,i);
      if (a && mxGetNumberOfElements(a)) return a;
   }
   return NULL;
};

   void mxPutArray(const char *file, int line,
      mxArray *a, const char *vname="ans", const char* ws="base",
      const char* dstr=""); 

   void mxPutAndDestroy(const char *file, int line,
      mxArray *a, const char *vname="ans", const char* ws="base");

namespace Mx {

int IsEqual(const mxArray *a, const char *s); 

template <class T> 
class Array { 

  public:

    Array() : ax(0), data(0),
    len(0), rank(0), trans(0), cmplx(0), mxref(0) { };

    Array(unsigned ndim, size_t *dims, char t_=0, char c_=-1)
     : ax(0), data(0), len(0), rank(ndim), trans(t_), cmplx(c_), mxref(0) {
       if (cmplx<0) { cmplx=ISCOMPLX_(T); }
       init(ndim,dims);
    };

    Array(size_t d1, size_t d2, char t_=0, char c_=-1)
     : ax(0), data(0), len(0), rank(2), trans(t_), cmplx(c_), mxref(0) {
       size_t dims[2]={d1,d2}; if (cmplx<0) { cmplx=ISCOMPLX_(T); }
       init(2,dims);
    };

    Array(const wbvector<size_t> &S, char t_=0, char c_=-1);

    Array(const mxArray* a, char t_=0, char c_=-1)
     : ax(0), data(0), len(0), rank(0), trans(t_), cmplx(c_), mxref(1) {
       if (cmplx<0) { cmplx=ISCOMPLX_(T); }
       init(a);
    };

    Array(const char *F, int L, const mxArray* a, char t_=0, char c_=-1)
     : ax(0), data(0), len(0), rank(0), trans(t_), cmplx(c_), mxref(1)
     { if (cmplx<0) { cmplx=ISCOMPLX_(T); }

       try { init(a); }
       catch (...) { wblog(F_L,"ERR %s() ",FCT); } 

       if (!ax || (!mxref && !data)) wblog(F_L,
       "ERR %s() invalid input\n%s %s",FCT,SHORT_FL,STR(*this));
    };

    wbstring toStr(char vflag=0, unsigned l=64) const;
    wbstring toStrT(char vflag=0) const;

    void init_xstr() {
       sprintf_str("Mx::Array<%s>",sTSTR(T));
       xstr=str;
    };

    const char *xStr() const {
       return xstr.c_str();
    };

   ~Array() {
       if ((ax || data) && mxref<=0) wblog(FL,
          "WRN %s() got ax=%p, data=%p !?",FCT,ax,data);
    }

    mxArray* Return(T*& dd) { dd=data; return Return(); };

    mxArray* Return() {  
       mxArray *a=ax;
          ax=NULL; data=NULL; len=0; 
          mxref=(mxref>0? 1:0);
       return a;
    };

    Array& init(unsigned ndim, size_t *dims, char t_=-1, char c_=-1);
    Array& init(const mxArray *a, char t_=-1, char c_=-1);

    template <class Tb> 
    mxArray* copyFromTR(const Tb* b, char tcheck=0) {

       if (!trans && typeid(T)==typeid(Tb) && data && b)
            { memcpy(data,b,len*sizeof(T)); }
       else { copy_from(b,tcheck); }

       return Return();
    };

    template <class Tb>
    int copy_from(const Tb* b, char tcheck=0);

    template <class Tb>
    size_t copyTo(Tb* b, char tcheck=0) const {
       if (!trans && data && b) { 
          if (mxref<=8) wblog(FL,"WRN %s() ref=%d !?",FCT,mxref);
          memcpy(b,data,len*sizeof(T));
       }
       else { copy_to(b,NULL,tcheck); }
       return len;
    };

    size_t copy_to(T* b, const wbperm *P=NULL, char tcheck=0) const;

    size_t rcopy_to(T* b, const wbperm *P=NULL, char tcheck=0) const;
    size_t zcopy_to(T* b, const wbperm *P=NULL, char tcheck=0) const;

    template <class Tx>
    size_t cpy_to_(T* b, const wbperm *P, char tcheck, const Tx *ad_) const;

    size_t ncopy_to(T* b, size_t n, char tcheck=0) {
       size_t l_=len; if (len>n) len=n;
       n=copy_to(b,NULL,tcheck); len=l_;
       return n;
    };

    mxArray *ax; 
    T *data;     
    size_t len;  
    unsigned rank;

    bool trans;  
    char cmplx;  
    char mxref;  

    std::string xstr;

  protected:
  private:

    void issue_type_error(
       const char *F, int L, unsigned ndim=0, const size_t *dims=NULL);

    int init_dptr();

    int safe_set_dptr(T *md) { data=md; return 0; }

    template <class TM>
    int safe_set_dptr(TM *md QS_UNUSED_VAR) { return 1; } 

    void MX_INIT_DATA(unsigned ndim, const size_t *dims);
};

template <> template <>
int Mx::Array<wbcomplex>::safe_set_dptr(mxComplexDouble *md) {
    data=(wbcomplex*)md; return 0;
};

template <> template <>
int Mx::Array<char>::safe_set_dptr(int8_T *md) {
    unsigned m=sizeof(char), n=sizeof(int8_T);
    if (m!=n) wblog(FL,"ERR %s() size mismatch (char: %d/%d)",FCT,m,n);
    data=(char*)md; return 0;
};

template <> template <>
int Mx::Array<size_t>::safe_set_dptr(uint64_T *md) {
    unsigned m=sizeof(size_t), n=sizeof(uint64_T);
    if (m!=n) wblog(FL,"ERR %s() size mismatch (size_t: %d/%d)",FCT,m,n);
    data=(size_t*)md; return 0;
};

template <> template <> 
int Mx::Array<unsigned long long>::safe_set_dptr(uint64_T *md) {
    unsigned m=sizeof(size_t), n=sizeof(uint64_T);
    if (m!=n) wblog(FL,"ERR %s() size mismatch (size_t: %d/%d)",FCT,m,n);
    data=(unsigned long long*)md; return 0;
};

template <> template <>
int Mx::Array<long long>::safe_set_dptr(int64_T *md) {
    unsigned m=sizeof(long long), n=sizeof(int64_T);
    if (m!=n) wblog(FL,"ERR %s() size mismatch (long: %d/%d)",FCT,m,n);
    data=(long long*)md; return 0;
};

template <>
size_t Mx::Array<double>::copy_to(double * b, const wbperm *P, char tcheck)
const { return rcopy_to(b,P,tcheck); };

template <>
size_t Mx::Array<int>::copy_to(int * b, const wbperm *P, char tcheck)
const { return rcopy_to(b,P,tcheck); };

template <>
size_t Mx::Array<unsigned>::copy_to(unsigned* b, const wbperm *P, char tcheck)
const { return rcopy_to(b,P,tcheck); };

template <>
size_t Mx::Array<size_t>::copy_to(size_t * b, const wbperm *P, char tcheck)
const { return rcopy_to(b,P,tcheck); };

template <>
size_t Mx::Array<char>::copy_to(char* b, const wbperm *P, char tcheck)
const { return rcopy_to(b,P,tcheck); };

template <> 
void Mx::Array<double>::MX_INIT_DATA(unsigned ndim, const size_t *dims) {
   data=0;
   #pragma omp critical (using_MEX_API)
   {  ax=mxCreateNumericArray(ndim,dims, mxDOUBLE_CLASS, mxREAL);
      if (ax) { safe_set_dptr(mxGetDoubles(ax)); }
   }; if (!data) issue_type_error(FL,ndim,dims);
};

template <> 
void Mx::Array<wbcomplex>::MX_INIT_DATA(unsigned ndim, const size_t *dims) {
   data=0;
   #pragma omp critical (using_MEX_API)
   {  ax=mxCreateNumericArray(ndim,dims, mxDOUBLE_CLASS, mxCOMPLEX);
      if (ax) { safe_set_dptr(mxGetComplexDoubles(ax)); }
   }; if (!data) issue_type_error(FL,ndim,dims);
};

template <>  
void Mx::Array<char>::MX_INIT_DATA(unsigned ndim, const size_t *dims) {
   data=0;
   #pragma omp critical (using_MEX_API)
   {  ax=mxCreateNumericArray(ndim,dims, mxINT8_CLASS, mxREAL);
      if (ax) { safe_set_dptr(mxGetInt8s(ax)); }
   }; if (!data) issue_type_error(FL,ndim,dims);
};

template <>  
void Mx::Array<signed char>::MX_INIT_DATA(unsigned ndim, const size_t *dims) {
   data=0;
   #pragma omp critical (using_MEX_API)
   {  ax=mxCreateNumericArray(ndim,dims, mxINT8_CLASS, mxREAL);
      if (ax) { safe_set_dptr(mxGetInt8s(ax)); }
   }; if (!data) issue_type_error(FL,ndim,dims);
};

template <> 
void Mx::Array<unsigned char>::MX_INIT_DATA(unsigned ndim, const size_t *dims) {
   data=0;
   #pragma omp critical (using_MEX_API)
   {  ax=mxCreateNumericArray(ndim,dims, mxUINT8_CLASS, mxREAL);
      if (ax) { safe_set_dptr(mxGetUint8s(ax)); }
   }; if (!data) issue_type_error(FL,ndim,dims);
};

template <> 
void Mx::Array<mxChar>::MX_INIT_DATA(unsigned ndim, const size_t *dims) {
   data=0;
   #pragma omp critical (using_MEX_API)
   {  ax=mxCreateCharArray(ndim,dims);
      if (ax) { safe_set_dptr(mxGetChars(ax)); }
   }; if (!data) issue_type_error(FL,ndim,dims);
};

template <>  
void Mx::Array<int>::MX_INIT_DATA(unsigned ndim, const size_t *dims) {
   data=0;
   #pragma omp critical (using_MEX_API)
   {  ax=mxCreateNumericArray(ndim,dims, mxINT32_CLASS, mxREAL);
      if (ax) { safe_set_dptr(mxGetInt32s(ax)); }
   }; if (!data) issue_type_error(FL,ndim,dims);
};

template <> 
void Mx::Array<unsigned>::MX_INIT_DATA(unsigned ndim, const size_t *dims) {
   data=0;
   #pragma omp critical (using_MEX_API)
   {  ax=mxCreateNumericArray(ndim,dims, mxUINT32_CLASS, mxREAL);
      if (ax) { safe_set_dptr(mxGetUint32s(ax)); }
   }; if (!data) issue_type_error(FL,ndim,dims);
};

template <> 
void Mx::Array<long>::MX_INIT_DATA(unsigned ndim, const size_t *dims) {
   data=0;
   #pragma omp critical (using_MEX_API)
   {  ax=mxCreateNumericArray(ndim,dims, mxINT64_CLASS, mxREAL);
      if (ax) { safe_set_dptr(mxGetInt64s(ax)); }
   }; if (!data) issue_type_error(FL,ndim,dims);
};

template <>
void Mx::Array<long long>::MX_INIT_DATA(unsigned ndim, const size_t *dims) {
   data=0;
   #pragma omp critical (using_MEX_API)
   {  ax=mxCreateNumericArray(ndim,dims, mxINT64_CLASS, mxREAL);
      if (ax) { safe_set_dptr(mxGetInt64s(ax)); }
   }; if (!data) issue_type_error(FL,ndim,dims);
};

template <> 
void Mx::Array<size_t>::MX_INIT_DATA(unsigned ndim, const size_t *dims) {
   data=0;
   #pragma omp critical (using_MEX_API)
   {  ax=mxCreateNumericArray(ndim,dims, mxUINT64_CLASS, mxREAL);
      if (ax) { safe_set_dptr(mxGetUint64s(ax)); }
   }; if (!data) issue_type_error(FL,ndim,dims);
};

template <>
void Mx::Array<unsigned long long>::MX_INIT_DATA(unsigned ndim, const size_t *dims) {
   data=0;
   #pragma omp critical (using_MEX_API)
   {  ax=mxCreateNumericArray(ndim,dims, mxUINT64_CLASS, mxREAL);
      if (ax) { safe_set_dptr(mxGetUint64s(ax)); }
   }; if (!data) issue_type_error(FL,ndim,dims);
};

int IsNumArray( 
   const char *F, int L,
   const mxArray* a, unsigned r=2, char cflag=0, char type='d'
);

bool IsDblMat(const char *F, int L, 
   const mxArray* a, char cflag=0, char dflag='d')
 { return Mx::IsNumArray(F,L,a,2,cflag,dflag)>0; }

bool IsDblArr(const char *F, int L, 
   const mxArray* a, char cflag=0, char dflag='d')
 { return Mx::IsNumArray(F,L,a,-1,cflag,dflag)>0; }

bool IsDblVector(const char *F, int L, 
   const mxArray *a, char cflag=0, char dflag='d')
 { return Mx::IsNumArray(F,L,a,1,cflag,dflag)>0; }

bool IsDblScalar(const char *F, int L, 
   const mxArray *a, char cflag=0, char dflag='d') {
   if (Mx::IsNumArray(F,L,a,1,cflag,dflag)<=0) return 0;
   return mxGetNumberOfElements(a)==1;
}

bool IsNumber(const char *F, int L, 
   const mxArray *a, char cflag=0, char dflag='*') {
   if (Mx::IsNumArray(F,L,a,1,cflag,dflag)<=0) return 0;
   return (mxGetNumberOfElements(a)==1);
}

bool IsScalar(const mxArray *a) { 
   return (mxGetNumberOfDimensions(a)==2 && mxGetNumberOfElements(a)==1);
};

char IsVector(const mxArray *a) { 
   int q=0; 
   unsigned n=mxGetNumberOfDimensions(a);
   if (n==2) {
      const size_t *s=mxGetDimensions(a);
      if (s[0]==1) { q|=1; }  
      if (s[1]==1) { q|=2; }  
      if (!q && !s[0] && !s[1]) { q=-1; }
   }
   else if (n<2) { q=-1; } 
   return q;
};

bool IsIndex(const mxArray *a, int base=0); 

template <class T>
bool IsIndex(const T *x, size_t n, int base=0) { 
   if (!WbUtil<T>().isPOD()) wblog(FL,
      "ERR %s() invalid data type `%s'",FCT,TSTR(T));

   if (WbUtil<T>().isInt()) {
      for (size_t i=0; i<n; ++i) { if (x[i]<base) return 0; }
   }
   else {
      for (size_t i=0; i<n; ++i) {
      if (x[i]<base || x[i]!=round(x[i])) return 0; }
   }
   return 1;
};

}; 

namespace Wb {
   int CallMatlab(
      int nargout, mxArray *argout[], int nargin, mxArray *argin[],
      const char *fctname);

   double CallMatlab(
      const char *F, int L, const char* cmd, const char *arg1=NULL);

   double CallMatlab(const char* cmd, const char *arg1=NULL) {
      return CallMatlab(0,0,cmd,arg1);
   };
};

#ifndef MATLAB_MEX_FILE

#undef  mexPutVariable
#define mexPutVariable mexPutVariable_wbx
int mexPutVariable_wbx(const char *ws, const char *vname, const mxArray *a){
    wblog(FL,"WRN mexPutVariable(%s,%lX,%s) not available outside MatLab",
    vname?vname:"",a,ws?ws:""); return 1;
};

#undef  mexCallMATLAB
#define mexCallMATLAB mexCallMATLAB_wbx
int mexCallMATLAB_wbx(
    int nargout, mxArray **argin, int nargin, mxArray **argout,
    const char *cmd
){  wblog(FL,"ERR mexCallMATLAB(%d,%lX,%d,%lX,%s) not available outside MatLab",
    nargin,argin,nargout,argout,cmd?cmd:""); return 1;
};

#undef  mexGetVariable
#define mexGetVariable mexGetVariable_wbx
int mexGetVariable_wbx(const char *ws, const char *vname) {
    wblog(FL,"WRN mexGetVariable(%s,%s) not available outside MatLab",
    vname?vname:"",ws?ws:""); return 1;
};

#undef  mexGetVariablePtr
#define mexGetVariablePtr mexGetVariablePtr_wbx
mxArray* mexGetVariablePtr_wbx(const char *ws, const char *vname) {
    wblog(FL,"WRN mexGetVariablePtr(%s,%s) not available outside MatLab",
    vname?vname:"",ws?ws:""); return NULL;
};

#undef  mexPrintf
#define mexPrintf mexPrintf_wbx
int mexPrintf_wbx(const char *fmt, ...){
    wblog(FL,"WRN mexPrintf(%s,...) not available outside MatLab",
    fmt?fmt:""); return 1;
};

#undef  mexFunctionName
#define mexFunctionName mexFunctionName_wbx

const char* mex_function_name="(program)";

const char* mexFunctionName_wbx(const char *fmt, ...){ 
    wblog(FL,"WRN mexFunctionNamef(%s,...) not available outside MatLab",
    fmt?fmt:""); return mex_function_name;
};

const char* mexFunctionName_wbx(){ 
    wblog(FL,"WRN mexFunctionNamef() not available outside MatLab");
    return mex_function_name;
};

#endif

#endif

