/* ---------------------------------------------------------------------
 * Project : QSpace tensor library (v4.0 pre-release)
 * Class   : QSpace MEX routines
 *
 * Copyright 2022 Andreas Weichselbaum
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

#ifndef __WB_MATLIB_CC__
#define __WB_MATLIB_CC__

/* -------------------------------------------------------------------- //
   NB! intended for mex routines that communicate with matlab desktop
   i.e. mex-files, but not programs e.g. with MAIN defined // Wb,Mar31,20
   NB! mexCallMATLAB(), mexDisp, mexPutVariable, mexGetVariablePtr, ...
// -------------------------------------------------------------------- */
#ifdef MATLAB_MEX_FILE
/* -------------------------------------------------------------------- //
   Safeguard: matlab `engine' may only be called from the same thread
   that also called the mex routint! see printf -> mexPrintf() error
   in tstmex_assert_error.cc // Wb,May06,19
// -------------------------------------------------------------------- */

int Wb::CallMatlab(
   int nargout, mxArray *argout[], int nargin, mxArray *argin[],
   const char *fctname
 ){
   if (Wb::my_caller_tid==omp_get_thread_num()) {
      return mexCallMATLAB(nargout,argout,nargin,argin,fctname);
   }
   else { fprintf(stdout,
     "mexCallMATLAB() to '%s' only allowed main thread (%d/%d)",
      fctname ? fctname : "(function?)", Wb::my_caller_tid,
      omp_get_thread_num());

      #ifdef DBSTOP
         dbstop(FL);
      #endif

      return -1;
   }
};

double Wb::CallMatlab(
   const char *F, int L, const char* cmd, const char *arg1) {

   double q=0; int e=0;
   mxArray *argout[1], *argin[1]={NULL};

   if (!cmd || !cmd[0]) wblog(F_L,"ERR %s() got empty cmd",FCT);

   if (Wb::my_caller_tid!=omp_get_thread_num()) {
      char s[16]=""; if (arg1) snprintf(s,16,"('%s')",arg1);

      fprintf(stdout, 
        "mexCallMATLAB() to %s%s only allowed main thread (%d/%d)",
         cmd, s, Wb::my_caller_tid, omp_get_thread_num()
      );
      #ifdef DBSTOP
         dbstop(FL);
      #endif

      return (q=-1);
   }

   if (arg1) {
      argin[0]=mxCreateString(arg1);
      try { e=mexCallMATLAB(1,argout,1,argin,cmd); }
      catch (...) { e=-1; }

      mxDestroyArray(argin[0]);
   }
   else {
      try { e=mexCallMATLAB(1,argout,0,argin,cmd); }
      catch (...) { e=-2; }
   }

   if (e) { int l=0, n=64; char s[n]; s[0]=0;
      if (e<0)
           l=snprintf(s,n,"invalid command >> %s",cmd);
      else l=snprintf(s,n,">> %s",cmd);

      if (l<n && arg1) l+=snprintf(s+l,n-l,"('%s')",arg1);
      if (l<n) l+=snprintf(s+l,n-l,
         e>0? "  %% returned error e=%d":" %% (e=%d)",e);

      if (F)
           wblog(F,L,"ERR %s %s",FCT,s);
      else wblog(FL, "WRN %s %s",FCT,s);

      return q=0;
   }

   mxGetNumber(argout[0],q);
   mxDestroyArray(argout[0]);

   if (F) { unsigned l, n=32; char s[n];
      l=snprintf(s,n,"%s",cmd);
      if (arg1 && l<n) l+=snprintf(s+l,n-l,"('%s')",arg1);
      wblog(F,L," *  %s >> %-24s %% %g",FCT,s,q);
   }

   return q;
};

void mxPutArray(const char *F, int L,
    mxArray *a, const char *vn, const char* ws, const char* dstr
 ){
    int i;

   #pragma omp critical (using_MEX_API)
    i=mexPutVariable(ws, vn, a);

    if (i) wblog(F_L,
       "ERR failed to put '%s' to '%s' (%d)",vn,ws,i);

    if (dstr!=NULL) {
       if (!strcmp(dstr,"destroy")) mxDestroyArray(a);
       else if (dstr[0]!=0)
       wblog(F_L,"ERR invalid dstr '%s'", dstr);
    }
}

void mexDisp(const mxArray *a, const char *vn) {
    if (vn && vn[0]) printf("\n%s = \n\n", vn);
    Wb::CallMatlab(0, NULL, 1, (mxArray**) &a,"display");
};

void mexIssueWRN(const char *s) { 
    char tag[32];
    snprintf(tag,32,"Wb:MEX:%.24s",myname); 
    mexWarnMsgIdAndTxt(tag, s && s[0] ? s : "");
};

wbstring sprint_info(const mxArray *a) {
   wbstring s(64); 
   unsigned l=0;

   if (!a) { strcpy(s.data,"(null)"); }
   else {
      unsigned ndim=mxGetNumberOfDimensions(a);
      const size_t *dims=mxGetDimensions(a);

      l+=snprintf(s.data,s.len,"%s array of type `%s'",
         SSTR(wbvector<size_t>(ndim,dims,'r')), mxGetClassName(a)
      );
   }

   if (l>=s.len) wblog(FL,
      "WRN %s() string out of bounds (%d/%d)",FCT,l,s.len);
   return s;
};

double wbtoc( 
   const char *F, int L, const char *istr, char restart
){
   double t=-1.; mxArray *a;

   if ((F==NULL) ? 1 : (F[0]==0 ? 1 : 0)) {
      Wb::CallMatlab(0, NULL, 0, NULL, "tic");
      return t;
   }

   Wb::CallMatlab(1,&a,0,NULL,"toc");
   Mx::Array<double>(a).ncopy_to(&t,1); 
   mxDestroyArray(a);

   if (istr) {
      if (!istr[0])
           wblog(F_L,"elapsed time: %.6f sec", t);
      else wblog(F_L,"elapsed time: %.6f sec (%s)", t, istr);
   }

   if (restart)
   Wb::CallMatlab(0, NULL, 0, NULL, "tic");

   return t;
};

template<class T>
int getNumGlobal(const char *vn, T &x, const char* ws){

    const mxArray *a=mexGetVariablePtr(ws, vn);
    if (!a || !Mx::IsNumber(0,0,a)) { return 1; } 
    else {
       double dbl=0;
       Mx::Array<double>(a).ncopy_to(&dbl,1); 

       x=T(dbl);
       if ((double)x!=dbl) wblog(FL,
         "WRN ()%s invalid value for type '%s' (%g)",FCT,TSTR(T),dbl);

       return 0;
    }
};

template<>
int getNumGlobal(const char *vn, wbstring &s, const char *ws) {
   const mxArray *a=mexGetVariablePtr(ws, vn);

   char istr[64];

   if (!a || !mxIsChar(a)) return 1;
   if (mxGetString(a,istr,63)) return 1;
   s=istr;

   return 0;
}

#else

int Wb::CallMatlab(
   int nargout, mxArray *argout[], int nargin, mxArray *argin[],
   const char *fctname
){
   wblog(FL,"ERR %s(%d,%p,%d,%p,'%s')\nnot available outsite MatLab",
   FCT,nargout,argout,nargin,argin,fctname?fctname:""); return 0;
};

double Wb::CallMatlab(
   const char *F, int L, const char* cmd, const char *arg1) {

   char s[16]=""; if (arg1) snprintf(s,16,"('%s')",arg1);
    wblog(FL,
     "ERR %s -> %s%s\nnot available outsite MatLab",
      FCT, cmd? cmd:"(cmd!?)",s);
   return 0;
};

double wbtoc(
   const char *F, int L, const char *istr,
   char restart __attribute__ ((unused))
){
   wblog(FL,"WRN %s(%s,%s) not available outsite MatLab",
   FCT,shortFL(F_L),istr ? istr:""); return 0;
}

void mxPutArray(const char *F, int L,
    mxArray *a, const char *vn, const char* ws,
    const char* dstr
){
   wblog(F_L,"ERR %s(%s,%s,%lX,%s) not available outsite MatLab",
   FCT,vn?vn:"",ws?ws:"",&a,dstr?dstr:"");
};

template<class T>
int getNumGlobal( 
   const char *vn  __attribute__ ((unused)),
   T &x            __attribute__ ((unused)),
   const char* ws  __attribute__ ((unused))
){
   wblog(FL,"ERR %s(%s,%s) not available outsite MatLab",
   FCT,vn?vn:"",ws?ws:""); return 0;
};

void mexDisp(const mxArray *a, const char *vn) {
   wblog(FL,"ERR %s(%s,%lX) not available outsite MatLab",
   FCT,vn ? vn:"",a);
};

#endif 
#endif 
