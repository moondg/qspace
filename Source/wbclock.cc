/* ---------------------------------------------------------------------
 * Project : QSpace tensor library (v4.0 pre-release)
 * Class   : QSpace clock routines
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

#ifndef __WB_CLOCK_CC__
#define __WB_CLOCK_CC__

/* -------------------------------------------------------------------- */
/* -------------------------------------------------------------------- */
// see also $JUL/WbModules::sec2str() // Wb,Sep22,16
// tags: sec2str, int2time, num2time, tstamp

inline wbstring Wb::sec2Str(double t) { 
   wbstring st(16); 
   unsigned l;

   if (t<100) { l=snprintf(st.data,st.len,"%.5g",t); } else
   if (t<3600) {
      double s=fmod(t,60), m=(t-s)/60;
      if (::fabs(s-::round(s))<1e-6)
           { l=snprintf(st.data,st.len,"%02.0f:%02f",  m,s); }
      else { l=snprintf(st.data,st.len,"%02.0f:%02.3f",m,s); }
   }
   else {
      double s=fmod(t,60);
      unsigned q=::round((t-s)/60), m,h,d;

      m=q%60; q=(q-m)/60;
      h=q%24; d=(q-h)/24;

      if (d)
           l=snprintf(st.data,st.len,"%d-%02d:%02d:%02.0f",d,h,m,s);
      else l=snprintf(st.data,st.len,   "%02d:%02d:%02.0f",  h,m,s);
   }

   if (l>=st.len) wblog(FL,"WRN %s() string out of bounds "
      "(len=%d/%d)\nstr = %s",FCT,l,st.len,st.data);
   return st;
};

void Wb::pause(const char *F, int L, double tsec) {
    if (tsec<=0) { return; }

    if (F) wblog(FL," *  %s() for %g secs",FCT,tsec); else
    if (tsec>99) wblog(F_L,"WRN %s() for %g secs !?",FCT,tsec);

    timespec t; double s=floor(tsec), n=floor((tsec-s)*1E9);
    t.tv_sec=s; t.tv_nsec=n;

    nanosleep(&t,0); 
};

mxArray* Wb::Clock::toMxS(const char *F, int L) {
   return MXPut(F,L).add(name,"istr")
     .add(sec2Str(gettime('c')),"cpu")
     .add(sec2Str(gettime(  )),"wall")
   .toMx();
};

mxArray* Wb::Clock::toMx(const char *F, int L) {
   return MXPut(F,L).add(name,"istr")
     .addP(numtoMx(gettime('c')),"cpu")
     .addP(numtoMx(gettime(  )),"wall")
   .toMx();
};

void Wb::save_and_clear_Profiling() {

   mxArray *S=mxCreateStructMatrix(1,1,0,NULL);

   for (auto I=gwb_Profs.begin(); I!=gwb_Profs.end(); ++I) {
      unsigned j=0, n=I->second.size(); double x[4];
      wbMatrix<double> M(n,4);
      for (auto J=I->second.begin(); J!=I->second.end(); ++J, ++j) {
         x[0]=J->first; 
         x[1]=J->second.ncall;
         x[2]=J->second.dt;
         x[3]=J->second.aux; if (x[1]>1) { x[2]/=x[1]; x[3]/=x[1]; }
         M.recSetP(j,x);
      }
      mxAddField2Scalar(FL,S,I->first.c_str(),M.SortRecs().toMx());
   }

   unsigned i=0, l=0,n=128; char fout[n]; fout[0]=fout[n-1]=0;
   const char *s;

   s=getenv("Wb_PROFILE_FOUT");
   if (s && s[0])
        { l=snprintf(fout,n,"%s",s); }
   else {
      const mxArray *a=mexGetVariablePtr("caller","profile_fout");
      if (a && mxGetString(a,fout,n-1)) { l=-1; }
   }
   if (l>=n) { fout[n-1]=0; wblog(FL,
      "WRN %s() string out of bounds (%d/%d) !?\ncheck %s",FCT,l,n,
      int(n)>0? "variable `profile_fout'":"env `WB_PROFILE_FOUT'");
      fprintf(stdout,"\n%s (%ld)\n",fout,strlen(fout)); fout[0]=0;
   }

   if (fout[0]) { Wb::matFile F;
      if (F.open(FL,fout,"w",0)) { 
         for (n=mxGetNumberOfFields(S), i=0; i<n; ++i) {
            F.put(FL, 
               mxGetFieldNameByNumber(S,i),
               mxGetFieldByNumber(S,0,i)
            );
         }
         F.close();
      }
      else {
         wblog(FL,"WRN %s() invalid file name\n%s",FCT,fout);
         fout[0]=0;
      }
   }

   static int ncall=0;

   if (1 || ++ncall==1) { strcpy(fout,"Iprof"); }
   else { sprintf(fout,"Iprof%02d",ncall); }

   mxPutAndDestroy(FL,S,fout,"base"); 
};

#endif

