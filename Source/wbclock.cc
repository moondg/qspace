/* ---------------------------------------------------------------------
 * Project : QSpace tensor library (v4.0)
 * Class   : QSpace clock routines
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

void Wb::get_Clock_name(wbstring &name, const char *s, char use_tag) {

   if (use_tag) { 
      unsigned l, n=16+(s ? strlen(s) : 0); char sx[n];
      l=snprintf(sx,n,"%.12s:",
        #ifdef PROG_TAG
          PROG_TAG
        #elif defined(myname)
          myname
        #else
          "???"
        #endif
      );
      if (s && l<n) {
         if (l && !strncmp(sx,s,l))
              { wblog(FL,"WRN %s() '%s :? %s'",FCT,sx,s); }
         else { l+=snprintf(sx+l,n-l,"%s",s); }
      }
      name=sx;
   }
   else { name=(s ? s : ""); }
};

Wb::Clock* Wb::ClockSet::insert(Wb::Clock* clk, char lflag) {

   if (!clk) { wblog(FL,"ERR %s() got null clock",FCT); }
   if (!clk->name) { wblog(FL,"WRN %s() got empty clock name",FCT); }

   auto q=buf.insert({ clk->name.data, clk });

   if (!q.second && !lflag) wblog(FL,
      "ERR %s() got existing clock entry \f'%s'",FCT,clk->name.data);
   if (q.first->second!=clk) wblog(FL,
      "ERR %s() clk inconsistency %p / %p",FCT,q.first->second,clk);
   return q.first->second;
};

int Wb::ClockSet::erase(Wb::Clock* clk) {
   if (!clk) { wblog(FL,"ERR %s() got null clock",FCT); }
   if (!clk->name) {
      wblog(FL,"ERR %s() got empty clock name",FCT);
   }

   return buf.erase(clk->name.data);
};

Wb::Clock* Wb::ClockSet::use(const char *istr, char mode) {
   Wb::Clock* clk; 

   if (!istr   ) wblog(FL,"ERR %s() got null clock string", FCT);
   if (!istr[0]) wblog(FL,"WRN %s() got empty clock string",FCT);

   auto i=buf.find(istr);
   if (i!=buf.end()) { clk=i->second; }
   else {
      clk = new Wb::Clock(istr,mode&2,0,NULL,NULL,0); 
      auto q=buf.insert({ istr, clk });         
      if (!q.second) wblog(FL,"WRN %s() failed to insert clock",FCT);
   }

   if (mode>1) {
      if (mode&1) { clk->resume(); } 
      if (mode>3) wblog(FL,"WRN %s() got mode = %s",FCT,cSTR(mode));
   }
   return clk;
};

void Wb::ClockSet::init(char vflag) {
   unsigned n=buf.size();

   if (vflag>3) {
      if (vflag=='v') { vflag=1; } else
      if (vflag=='V') { vflag=2; } else
      if (vflag=='q') { vflag=0; } else
      wblog(FL,"WRN %s() got vflag=%s",FCT,cSTR(vflag));
   }

   for (auto i=buf.begin(); i!=buf.end(); ++i) {
      if (!vflag) { i->second->init(); }
      delete i->second; 
   }
   if (n) { PRINTF("\n"); }
};

unsigned Wb::ClockSet::info(char vflag) const {
   for (auto i=buf.begin(); i!=buf.end(); ++i) {
      i->second->info(NULL,vflag);
   }
   return buf.size();
};

int Wb::ClockSet::info_u(unsigned u, char rflag) {
   unsigned m=0; 

   for (auto i=buf.begin(); i!=buf.end(); ++i) {
      Wb::Clock &clk = *(i->second);
      if ((u && (clk.user & u)) || u==clk.user) { ++m;
          clk.info();  if (rflag) {
          clk.reset(); }
      }
   }
   return m;
};

int Wb::ClockSet::reset_u(unsigned u, char iflag) {
   unsigned m=0; 

   for (auto i=buf.begin(); i!=buf.end(); ++i) {
      Wb::Clock &clk = *(i->second);
      if ((u && (clk.user & u)) || u==clk.user) { ++m;
          if (iflag) { clk.info(); }
          clk.reset();
      }
   }
   return m;
};

Wb::Clock* Wb::ClockSet::get(const char *istr, char use_tag) {
   wbstring s;
   get_Clock_name(s,istr,use_tag);

   auto i=buf.find(s.data);
   return ( i!=buf.end() ? i->second : NULL);
};

int Wb::ClockSet::reset(const char *istr, char use_tag) {

   int q=0; 
   Wb::Clock *clk=get(istr,use_tag);
   if (clk) { clk->reset(); ++q; }
   return q;
};

mxArray* Wb::ClockSet::toMx() const {
   unsigned i=0, n=buf.size();
   mxArray* C=mxCreateCellMatrix(n,3);

   for (auto it=buf.begin(); it!=buf.end(); ++it, ++i) {
      const Wb::Clock &ci=*(it->second);   
      mxSetCell(C,     i, ci.name.toMx()); 
      mxSetCell(C,   n+i, numtoMx(ci.gettime('c')));
      mxSetCell(C, 2*n+i, numtoMx(ci.gettime(   )));
   }
   return C;
};

mxArray* Wb::Clock::toMxS(const char *F, int L) const {
   return MXPut(F,L)
      .add(name,"istr")
      .add(sec2Str(gettime('c')),"cpu")
      .add(sec2Str(gettime(  )),"wall")
   .toMx();
};

mxArray* Wb::Clock::toMx(const char *F, int L) const {
   return MXPut(F,L).add(name,"istr")
     .addP(numtoMx(gettime('c')),"cpu" )
     .addP(numtoMx(gettime(   )),"wall")
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

   const unsigned flen=128;
   unsigned i=0, l=0; char fout[flen]; fout[0]=fout[flen-1]=0;
   const char *s;

   s=getenv("Wb_PROFILE_FOUT");
   if (s && s[0])
        { l=snprintf(fout,flen,"%s",s); }
   else {
      const mxArray *a=mexGetVariablePtr("caller","profile_fout");
      if (a && mxGetString(a,fout,flen-1)) { l=-1; }
   }
   if (l>=flen) { fout[flen-1]=0; wblog(FL,
      "WRN %s() string out of bounds (%d/%d) !?\ncheck %s",FCT,l,flen,
      int(l)>0? "variable `profile_fout'":"env `WB_PROFILE_FOUT'");
      fprintf(stdout,"\n%s (%ld)\n",fout,strlen(fout)); fout[0]=0;
   }

   if (fout[0]) { Wb::matFile F;
      if (F.open(FL,fout,"w",0)) { 
         unsigned n=mxGetNumberOfFields(S);
         for (i=0; i<n; ++i) {
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
   else { snprintf(fout,flen,"Iprof%02d",ncall); }

   mxPutAndDestroy(FL,S,fout,"base"); 
};

#endif

