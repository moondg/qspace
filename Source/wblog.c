/* ---------------------------------------------------------------------
 * Project : QSpace tensor library (v4.0)
 * Class   : wblog (logging routines)
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

#ifndef __WB_LOG_C__
#define __WB_LOG_C__

/* -------------------------------------------------------------------- */
/* -------------------------------------------------------------------- */

#ifndef LIBQS_SETUP

void usage(const char *F, int L, const char* estr) {

    int l=0;
    while (isspace(USAGE[l])) { ++l; }

    if (USAGE[l]) { printf("\n%s%s\n",l?"":"   ",USAGE); }
    else {
       const char *f=mexFunctionName();
       if (f && f[0]) {
          char cmd[256]; l=snprintf(cmd,256,
             "if exist('%s.m','file')==2, fprintf(1,'\\n'); help %s;\n"
             "else fprintf(1,'"
             "   ERR failed to find help %s.m file for mex function %s()\\n'); "
             "end", f,f,f,f);

          if (l<256) {
             if (F && F[0]) { mexPrintf("\n   usage %s\n",shortFL(F,L)); }
             mexEvalString(cmd);
          }
          else wblog(FL,"WRN %s() got f=%s (l=%d/256) !?",FCT,f,l);
       }
    }

    if (estr && estr[0]) { wbdie(F,L,estr); }
};

#endif

const char* shortFL(const char *F, int L,
   unsigned n, 
   const char *P, char sep 
) {

   static char fstr[512]="", *fs=fstr; 
   unsigned i,l, j=0, nx=(P? 40 : 32);
   char *sout=fstr, sx[nx]; 

   if (int(n)<0) { n=WBL_HLEN; } 
   else if (n<16 || n>=nx) { 
      if (!n && !F && !L)
           { memset(fstr,0,512); return sout; }
      else { wblog(FL,"ERR %s() got n=%d/%d\n\n",FCT,n+1,nx); }
   }

  #pragma omp critical (using_shortFL)
   { if (fs+n-fstr<512) 
          { sout=fs; fs+=(n+1); } 
     else { sout=fs=fstr; }
   }

   if (!F || !F[0]) { sout[0]=0; return sout; }
   if (F[0]=='/' || F[0]=='.') { l=strlen(F);
      for (i=l-1; i<l; --i) { if (F[i]=='.') { break; }}
      for (     ; i<l; --i) { if (F[i]=='/') { j=i+1; break; }}
   }

   if (P && !P[0]) { 
      P=mytag; if (!P || !P[0]) { P=myname; }
     #ifdef PROG
      if (P && P[0]) { l=strlen(P);
         if (!strncmp(F+j,P, l>3? 3:l)) { P=NULL; }
      }
     #endif
   }

   if (!P || !P[0]) {
      if (L>0)
           { l=snprintf(sx,nx,"%.24s:%d",F+j,L); }
      else { l=snprintf(sx,nx,"%.24s:",  F+j  ); }
   }
   else {
      if (L>0)
           { l=snprintf(sx,nx,"%.3s%c%.20s:%d",P,sep,F+j,L); }
      else { l=snprintf(sx,nx,"%.3s%c%.20s:",  P,sep,F+j  ); }
   }

   if (l>=nx) wblog(FL,
      "WRN %s() string out of bounds (%d/%d)\n%s",FCT,l,nx,sx);

   if (l<=n) { strcpy(sout,sx); } 
   else {
      for (i=l-1; i<l; --i) { if (sx[i]=='.') { break; }}
      for (j=i;   j<l; --j) { if (sx[j]=='_' || sx[j]=='-') { break; }}
      if (j<l && (l-j)<n-3) { 
         i=n+j-l; 
         strncpy(sout,sx,i);
         sout[i]='\''; 
         strcpy(sout+i+1,sx+j+1);
      }
      else if (i<l && (l-i)<n-8) { 
         j=n+i-l; 
         strncpy(sout,sx,j);
         sout[j]='\''; 
         strcpy(sout+j+1,sx+i+1);
      }
      else { i=5;
         strncpy(sout,sx,i); sout[i]='\'';
         strcpy(sout+i+1,sx+l-n+i+1);
      }
   }
   return sout;
};

const char* tmpmatFL(
   const char *F, int L, const char *vn, char *s, unsigned n) {

   if (!s || int(n)<8) wblog(FL,"ERR %s() invalid input (%p, %d)",FCT,s,n);
   unsigned l=snprintf(s,n,"./tmp");

   if (l<n && F && F[0]) { unsigned i=0, k=0;
	  for (; F[i]; ++i) { if (F[i]=='/') { k=i+1; } }
      i=l; l+=snprintf(s+l,n-l,"-%s",F+k);
	  for (; s[i] && i<n; ++i) { if (s[i]=='.') { s[i]=0; l=i; break; }}
   }

   if (l<n) { l+=snprintf(s+l,n-l,"-%d",L); } 
#ifdef QS_USING_OMP
   if (l<n) { 
      int j=omp_get_thread_num();
      if (j!=Wb::my_caller_tid) { l+=snprintf(s+l,n-l,".%d",j); }
   }
#endif

   if (l<n && vn && vn[0]) { l+=snprintf(s+l,n-l,"-%s",vn); }
   if (l<n) { l+=snprintf(s+l,n-l,".mat"); }

   if (l>=n) wblog(F_L,"WRN %s() "
      "string out of bounds (%d/%d)\n%s",FCT,l,n,s?s:"(null)");
   return s;
};

unsigned Wb::VersionInfo::str_cpy(char *s, unsigned n, const char *a) {
   unsigned i=0; s[0]=0;
   if (a) {
      for (; i<n; ++i) { if (!(s[i]=a[i])) break; }
      if (i>=n) { wblog(FL,
         "WRN %s() out of bounds (%d/%d)\n%s",FCT,strlen(a),n,a);
      }
   }
   return i;
};

Wb::VersionInfo& Wb::VersionInfo::init() {

   unsigned i,l,m,n; char *s;

   str_cpy(fctn,sizeof(fctn),myname);
  #ifdef PROG_TAG
   str_cpy(tag, sizeof(tag), PROG_TAG);
  #endif

   n=str_cpy(matlab,sizeof(matlab),MEX_MLVER_STR);

   s=matlab;
   for (l=i=0; i<n; ++i) { if (s[i]=='/') { l=i+1; }}
   if (l) {
      for (i=l; i<=n; ++i) { if (!(s[i-l]=s[l])) break; }
   }

   l=sizeof(os);
   #if __linux__
      snprintf(os,l,"linux / %s",  PP_STRFY(MEX_EXT));
   #elif __unix__
      snprintf(os,l,"unix / %s",   PP_STRFY(MEX_EXT));
   #elif __APPLE__
      snprintf(os,l,"macOS / %s",  PP_STRFY(MEX_EXT));
   #else
      snprintf(os,l,"unknown / %s",PP_STRFY(MEX_EXT));
   #endif

#ifdef QS_USING_MPFR
   str_cpy(mpfr,sizeof(mpfr),MPFR_VERSION_STRING);
#endif

   n=sizeof(compiler); l=6;

#if defined(__clang__)
   snprintf(compiler,n,"clang/%d.%d.%d",
    __clang_major__, __clang_minor__, __clang_patchlevel__);
#elif defined(__GNUC__)
   snprintf(compiler,n,"gcc/%d.%d.%d",
   __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
#endif

#ifdef LD_CLEBSCH_QS
   n=sizeof(qspace);
   l=snprintf(qspace,n,"%.1f",double(QS_VERSION));

   i=strlen(QS_VERSION_SUB_);
   if (QS_VERSION_SUB || !i)
          { l+=snprintf(qspace+l,n-l,".%g",double(QS_VERSION_SUB)); }
   if (i) { l+=snprintf(qspace+l,n-l,"-%s",QS_VERSION_SUB_); }
#endif

#ifdef __QS_GIT_BRANCH__
   str_cpy(git,sizeof(git),PP_STRFY(__QS_GIT_BRANCH__));
#endif

   n=sizeof(flags); l=0;

#ifdef MAIN
   l+=snprintf(flags+l,n-l,"%s","MAIN");
#endif
#ifdef MATLAB_MEX_FILE
   l+=snprintf(flags+l,n-l,"%s%s",l? ", ":"","MEX_FILE");
#endif

   i=l; m=0;
   l+=snprintf(flags+l,n-l," using");

#ifdef LD_CLEBSCH_QS
   if (l<n) { l+=snprintf(flags+l,n-l," QS"  ); ++m; }
#endif
#ifdef QS_USING_OMP
   if (l<n) { l+=snprintf(flags+l,n-l," OMP" ); ++m; }
#endif
#ifdef QS_USING_MPFR
   if (l<n) { l+=snprintf(flags+l,n-l," MPFR"); ++m; }
#endif
#ifdef QS_USING_HPTT  
   if (l<n) { l+=snprintf(flags+l,n-l," HPTT"); ++m; }
#endif

   if (m)
        { i=l; m=0; if (l+3<n) { strcpy(flags+l,", "); l+=2; }}
   else { l=i; m=1; } 

#ifdef WB_SKIP_ASSERT
   if (l<n) {
      l+=snprintf(flags+l,n-l,"%s%s",m? ", ":"","SKIP_ASSERT");
   }; ++m;
#endif
#ifdef MEX_CXX_FLAGS 
   if (l<n) {
      l+=snprintf(flags+l,n-l,"%s%s",m?", ":"",PP_STRFY(MEX_CXX_FLAGS));
   }; ++m;
#endif
#ifdef __WBDEBUG__
   if (l<n) {
      l+=snprintf(flags+l,n-l,"%s%s",m?", ":"","DEBUG");
   }; ++m;
#endif
#ifdef DBSTOP
   if (l<n) {
      l+=snprintf(flags+l,n-l,"%s%s",m?", ":"","DBSTOP");
   }; ++m;
#endif
#ifdef DBG_GCX_LOCKS
   if (l<n) { l+=snprintf(flags+l,n-l,
      "%s%s=%g",m?", ":"","DBG_GCX_LOCKS", double(DBG_GCX_LOCKS));
   }; ++m;
#endif
#ifdef DBG_QSX_BUF
   if (l<n) { l+=snprintf(flags+l,n-l,
      "%s%s=%g",m?", ":"","DBG_QSX_BUF", double(DBG_QSX_BUF));
   }; ++m;
#endif
#ifdef __WB_MEM_CHECK__
   if (l<n) {
      l+=snprintf(flags+l,n-l,"%s%s",m?", ":"","MEM_CHECK");
   }; ++m;
#endif

   if (m)
        { i=l; m=0; }
   else { l=i; }

   l=snprintf(compiled,sizeof(compiled),"%s on %s",
     #ifdef __DATE_COMP__
        PP_STRFY(__DATE_COMP__),
     #else
        __DATE__,
     #endif
     HOST_NAME_STR
   );

   return *this;
};

void Wb::VersionInfo::print() const { PRINTF("\n%s\n",toStr().data); };

wbstring Wb::VersionInfo::toStr() const {

   wbstring sout; unsigned i=0,l,n=512;
   sout.init(n);

   const unsigned m=8; char q[2], *sx=sout.data;

   const char *field[m] = {
      "Matlab","Compiler","System","MPFR","QSpace","Git-branch","Flags","Compiled"
   };
   const char *data[m] = {
      matlab,   compiler,  os,     mpfr,  qspace,  git,         flags,  compiled
   };

   l=snprintf(sx,n,"Version info for function `%s'",fctn);
   if (*tag)
        { l+=snprintf(sx+l,n-l," (%s)\n\n",tag); }
   else { l+=snprintf(sx+l,n-l,"\n\n"); }

   q[0]=q[1]=0;

   for (; i<m; ++i) { if (data[i] && l<n) {
      if (i==0) { *q='R'; } else
      if (i==3 || i==4) 
           { if (*data[i]) { *q='v'; } else continue; }
      else { *q=0; }

      l+=snprintf(sx+l,n-l,"   %-13s %s%s\n",field[i],q,data[i]);
   }}

   if (l<n)
        { if (l) { sout.len=l; }}
   else { wblog(FL,"WRN %s() string out of bounds (%d/%d)",FCT,l,n); }

   return sout;
};

mxArray* Wb::VersionInfo::toMx() const {

   unsigned i=0; const unsigned n=10;

   const char *fields[n] = {
       "function","tag","matlab","system","MPFR","compiler","QSpace",
       "git_branch", "flags", "compiled"
       };
   mxArray *S=mxCreateStructMatrix(1,1,n,fields);

   const char *data[n] = {
       fctn, tag, matlab, os, mpfr, compiler, qspace, git,
       flags, compiled
   };

   for (; i<n; ++i) {
      mxSetFieldByNumber(S,0,i,wbstring(data[i]).toMx());
   }

   return S;
};

void banner(unsigned n, const char *s, const char *istr, const char *fstr){

   unsigned i,r=0,l=0,m=0, s1=0, s2=0;

   if (s && s[0]) r=strlen(s);
   else wblog(FL,"ERR cannot repeat empty string!");

   if (istr && istr[0]) { s1=strlen(istr); printf("%s",istr); }
   if (fstr && fstr[0]) { s2=strlen(fstr); }

   m=s1+s2; if (n>=m) l=(n-m)/r;
   else wblog(FL,"ERR overall length too short (%d/%d)",n,m);

   for (i=0; i<l; i++) printf("%s",s);
   for (i=l*r+s1, l=n-s2; i<l; i++) printf("%s",s);

   if (fstr && fstr[0]) { printf("%s",fstr); }
   printf("\n");
};

void wbdie(const char *F, int L, const char* istr) {

   PRINTF("\n%s ERR\n\n%s\n\n", shortFL(F,L), istr);

  #ifdef QS_USING_OMP
   throw Wb::LogException(ERR); 
  #else
   ExitMsg();
  #endif
};

void wbdie(
   const char *F, int L,
   const char *istr1,
   const char *istr2
){
   wblog(F,L,"ERR %s%N%N%s%N", istr1, istr2);
  #ifdef QS_USING_OMP
   throw Wb::LogException(ERR); 
  #else
   ExitMsg();
  #endif
};

Wb::LogException::LogException(WBLOG_TYPE l, const char *s)
 : type(l), nrefs(NULL), count(0), ith(-1), nth(-1) {

   if (!s || !s[0]) { istr[0]=0; } else {
      if (snprintf(istr,64,"%s",s)>=64) {
         istr[61]=istr[62]='.'; istr[63]=0;
      }
   }

   if (l) { nrefs = new int; (*nrefs)=1; count=1;
      ith=omp_get_thread_num();
      nth=omp_get_num_threads();
     #ifdef DBSTOP
      fprintf(stdout,"%s %s dbstop() having Wb::LogException "
         "nrefs=%d, istr='%s'; tid=%d/%d\n",
         shortFLT,*nrefs,istr,ith,nth); fflush(0);
      dbstop(FL); 
     #endif
   }
};

void Wb::LogException::report(const char *F, int L, const char *func) {
   if (type || count) { 
      unsigned l, n=128; char s[n]; 

      l=snprintf(s,n,
         "%s ERR %s() at nerr=%d",shortFL(F_L),func?func:FCT,count);
      if (l<n) {
         if (ith<=1 && nth<=1)
              l+=snprintf(s+l,n-l," (serial @ %d/%d)",ith,nth);
         else l+=snprintf(s+l,n-l," (parallel section thread %d/%d)",ith,nth);
      }

      if (type)  { PRINTF("\n%s\n\n",s); throw(*this); } else
      if (count) { ExitMsg(s); }
   }
   init();
};

wbstring Wb::LogException::toStr() const {
   wbstring sout(128); 

   unsigned l, n=sout.len; char *s=sout.data;

   l=snprintf(s,n,"%s: ",istr[0] ? istr : "''");
   if (l<n) {
      if (type<NUM_TYPE)
           { l+=snprintf(s+l,n-l,"type=%s",WBLOG_TYPE_STR[type]); }
      else { l+=snprintf(s+l,n-l,"type=%d !?",unsigned(type)); }
   }
   if (l<n) {
     l+=snprintf(s+l,n-l,", %s (i=%d/%d; nerr=%d)",
     ith<=1 && nth<=1 ? "serial":"parallel section", ith,nth,count);
   }
   if (l<n) {
      if (nrefs)
           { l+=snprintf(s+l,n-l,", *nrefs=%d",*nrefs); }
      else { l+=snprintf(s+l,n-l, ", nrefs=null"); }
   }
   return sout;
};

void init_header(
   char *hstr, unsigned hlen, 
   const char *F, int L,
   const char *time_stamp, const char *tag, int xcol,
   unsigned lenFL=21 
){
   unsigned i,k=0,n; 

   if (xcol && hlen>12) { 
      i=Wb::termcolor::ID("hdr");
      if (int(i)>0) {
         if (i<8)
              { i=snprintf(hstr,hlen,"\e[3%dm",i); }
         else { i=snprintf(hstr,hlen,"\e[38;5;%dm",i); } 
      }
      else {
         i=snprintf(hstr,hlen,"\e[0m");
      }
      hstr+=i; hlen-=i;
   }

   if (F) {
      if (F[0]=='/' || (F[0]=='.' && F[1]=='/')) {
         for (i=0; F[i]; ++i) { if (F[i]=='/') k=i+1; }
      }

      if (L==WBL_GOT_SHORTFL) 
           { i=snprintf(hstr,hlen,"%.100s",   F+k  ); }
      else { i=snprintf(hstr,hlen,"%.100s:%d",F+k,L); }

      if (i>=hlen) wblog(FL,"ERR %s() "
         "string out of bounds (%d/%d)%N%N`%s'",FCT,i,hlen,hstr);

      n=strlen(hstr); if (n>lenFL) {
          i=lenFL/2; hstr[i]=hstr[i+1]='.'; i+=2;
          for (k=n-lenFL/2+2; k<n; ++k, ++i) hstr[i]=hstr[k];
          hstr[i]=0;
      }
   }
   else { i=snprintf(hstr,hlen,"(null):%d", L); }

   i=strlen(hstr); for (; i<lenFL; ++i) { hstr[i]=' '; }

   n=strlen(time_stamp+17); if (n>32) { n=32; } 
   if (i+n+12<hlen) { 
      memcpy(hstr+i,time_stamp+17,n); i+=n;
      hstr[i]=' '; hstr[++i]=0;
      hstr+=i; hlen-=i; i=0;
   }
   else wblog(FL,"ERR %s() "
     "string out of bounds (%d/%d)%N%N`%s'",FCT,i,hlen,hstr);

   if (xcol) {
      if (xcol>0) {
         if (xcol<8) 
              { i=snprintf(hstr,hlen,"\e[3%dm",xcol); }
         else { i=snprintf(hstr,hlen,"\e[38;5;%dm",xcol); }
      }
      else if (xcol<0) { i=snprintf(hstr,hlen,"\e[0m"); }
   }

   if (i<hlen) {
      if (tag && tag[0])
           { i+=snprintf(hstr+i,hlen-i," %s ",tag); }
      else { hstr[i]=' '; hstr[++i]=0; } 
   }
   if (i>=hlen) wblog(FL,"ERR %s() "
     "string out of bounds (%d/%d)%N%N`%s'",FCT,i,hlen,hstr);
};

int check_update_header(
   char *hstr, const char* fmt, const char *istr=NULL, int lenFL=21
){
   int n,k=0,i=0,q=0;

   for (; fmt[i]; ++i) {
      if (isalpha(fmt[i])) q|=1; else if (fmt[i]=='.') q|=2;  else
      if (isdigit(fmt[i])) q|=8; else if (fmt[i]=='_') q|=16; else
      if (fmt[i]=='/') { q|=4; k=i+1; } else break;
   }

   if ((q&3)!=3 || fmt[i]!=':') return 0;
   n=i++;

   for (; fmt[i]; ++i) if (!isdigit(fmt[i])) break;
   if (i==n+1 || isalpha(fmt[i])) return 0;

   if (istr && istr[0]) { 
      strcpy(hstr,istr);
      n=strlen(hstr); hstr+=n; lenFL-=n;
   }

   if (i-k<lenFL) {
      memcpy(hstr,fmt+k,i-k);
      for (k=i-k; k<lenFL; k++) hstr[k]=' ';
   }
   else {
      memcpy(hstr,fmt+k,lenFL/2); k=lenFL/2;
      hstr[k]=hstr[k+1]='.';
      memcpy(hstr+k+2, fmt+i-k+2,k-2); k=2*k;
   }

   if (k<=lenFL) {
      for (; hstr[k]; k++) if (hstr[k]!=' ') break;
      if (k>lenFL) {
         for (n=k-lenFL; hstr[k]; k++) hstr[k-n]=hstr[k];
         hstr[k-n]=0;
      }
   }

   for (; fmt[i]; ++i) {
      if (!isspace(fmt[i]) && fmt[i]!=':' && fmt[i]!='(' && fmt[i]!=')')
      break;
   }

   return i;
};

   void wbSetLogLevel(unsigned l) {
       wblogf(stdout,0,0,"log_level",l);
   };

int wblogf(FILE *fid,
   const char* file, int line, const char *fmt, va_list args
) {
   #ifdef QS_USING_OMP
      Wb::ompGuard myLK(wblog_lock,1); 
   #endif

   Wb::SBUF S; int l; Wb::LogException e;

   WBL_COLOR_SCHEME xcol=WLC_OFF;
   if (Wb::useCol && (fid==stdout || fid==stderr)) { xcol=WLC_DARK; }

   try { l=wblogs(S,xcol,file,line,fmt,args); }
   catch (Wb::LogException &e_) { e=e_; l=-11; }
   catch (...) { l=-12; }

    { S.flush(fid, l>1); }

   if (e.type) { throw(e); } else
   if (l<-10) { 
      sprintf_str("ERR wblog.h:%d encountered l=%d !?",__LINE__,l);
      ExitMsg(str);
   }

   return l;
};

int wblog1(const char* file, int line, const char *fmt, ...) {

   va_list args; Wb::ARGV wd(&args); 
   va_start(args,fmt);

   WBL_COLOR_SCHEME xcol=WLC_OFF;
   if (Wb::useCol) { xcol=WLC_DARK; }

   Wb::SBUF S; int l; Wb::LogException e;
   try { l=wblogs(S,xcol,file,line,fmt,args,1); } 
   catch (Wb::LogException &e_) { e=e_; l=-11; }
   catch (...) { l=-12; }

   S.flush(l>3 ? stderr : stdout,l>1); 

   if (e.type) { throw(e); } else
   if (l<-10) { 
      sprintf_str("ERR wblog.h:%d encountered l=%d !?",__LINE__,l);
      ExitMsg(str);
   }

   return l;
};

int wblogs( 
    Wb::SBUF &Sb, WBL_COLOR_SCHEME xcol_,
    const char* file, int line,
    const char *fmt, va_list args,
    char Hflag 
){
    static int log_level=1;
    static int then=0;

    int rval=0, xcol=0;

    unsigned i,j,k,l,m, nesc=0, hlen=128, flen=128;
    const unsigned nt=32;

    char isfmt=0, hflag=1, bflag=0, iflag=0, eflag=0, wflag=0, fwd=0;

    char fstr[flen];  
    char time_stamp[nt], tag[8];

    char c, *cp, cb, log_header[hlen];

    if (file==NULL && line<=0) {
       if (!strcmp(fmt,"log_level")) {
          log_level=(unsigned)va_arg(args,int);
          return 0;
       }
    }

    if (log_level<0) { return 0; }

    if (!fmt) wblog(FL,"ERR wblog got null fmt");

    if (fmt[0]=='\\') { hflag=0; fmt+=1; bflag=1; } 

    for (; *fmt; ++fmt) { if (*fmt=='\n') Sb.cat("\n"); else break; }

    tag[0]=0;

    eflag=wblog_check_tag(fmt,"ERR",tag); if (!eflag) {
    wflag=wblog_check_tag(fmt,"WRN",tag); if (!wflag) {
       if (fmt[0]=='\b') { fmt+=1; } 
       else {
           unsigned i=0, k,l;
           for (; fmt[i]; ++i) { if (fmt[i]!='\n') break; }
           for (k=i, l=k+3; i<l && fmt[i]; ++i) {
               if (fmt[i]=='%' || !isprint(fmt[i])) break; }

           if (i==l && fmt[i]==' ') {
              memcpy(tag,fmt+k,3); tag[3]=0; iflag=i+1;
           }

           if (xcol_ && iflag) {
              xcol=Wb::termcolor::ID(tag);
           }
        }
    }
    else { iflag=wflag; }}

    if ((log_level & 15)==0) { 
    if (!eflag && !iflag) return 0; }

    if (xcol_ && !xcol) { 
       xcol=Wb::termcolor::ID(tag); 
    }

    if (eflag || iflag) fmt+=(eflag+iflag);
    if (eflag) Sb.cat("\n");

   #pragma omp critical (got_CPTR_TIME)
   { time_t curr_time=time(NULL);
     struct tm *tblock=localtime(&curr_time);

     l=strftime(time_stamp,nt,"%a, %d %b %Y %T",tblock); 

     if (!l) {
        snprintf(time_stamp,nt,"Sat 22 Jul 1972 12:00:00%s",
        WBLOG_MSEC ? ".000" : "");
     }
     else if (WBLOG_MSEC && l+4<nt) {
        using namespace std::chrono;
        time_point<system_clock> t=system_clock::now();

        snprintf(time_stamp+l,nt-l,".%03ld", 
        size_t(duration_cast<milliseconds>(t.time_since_epoch()).count())%1000);

     }

     if (then && then!=tblock->tm_yday) { 
        char dstr[nt];
        l=strftime(dstr,nt,"%d-%b-%Y %T", tblock);  
        Sb.catf("\n>> TODAY %s (%s)\n\n",l? dstr:"??",myname); 

     }
     then=tblock->tm_yday; 
   } 

   if (!fmt[0] && !eflag && !iflag) {
       return 0;
   }

   init_header(log_header,hlen,file,line,time_stamp,tag,
     xcol ? xcol : (xcol_ ? -1 : 0));

   for (;;) {
      if (fmt[0]=='%' && fmt[1]=='N') { Sb.cat("\n"); fmt+=2; }
      else break;
   }

   if (hflag) { Sb.cat(log_header); hflag=0; }

   i=isfmt=0;

   while ((c=(*fmt++))) {
      if (c=='\f') { fwd|=1; continue; } 
      if (c=='\'') { fwd ^= 2; } else
      if (c=='\"') { fwd ^= 4; }

      fstr[i++]=c;

      if (i>124) { fstr[i-1]=0;
         snprintf(str,512,
           "ERR fstr out ouf bounds %d/128 (%s'%s)",i,fstr,fmt-1);
         ExitMsg(str);
      }

      if (c=='%') {
         if (fmt[0]=='%' && !isfmt) { fstr[i]=fmt[0]; ++fmt; ++i; } else
         if (i<2 || fstr[i-2]!='\\') { ++isfmt; }
         continue;
      }

      if (c=='\n') {
          fstr[i-1]=0; Sb.catf("%s\n",fstr);
          if (fmt[0]=='\n') {
             for (j=1; fmt[j]; ++j) { if (fmt[j]!='\n') break; }
             if (!fmt[j]) { Sb.cat(fmt); break; }
          }
          else {
             fmt+=check_update_header(log_header,fmt);
             Sb.cat(log_header);
             i=isfmt=0;
          }
      }
      else if (c=='\e' || (c=='\\' && *fmt=='e')) { ++nesc; }

      if (!isfmt || isdigit(c) || strchr("%-+. lL",c)) { continue; }

      if (c=='s') { 
         fstr[i]=0;
         if (i>2 && fstr[i-2]=='%') {
            fstr[i-2]=0; Sb.cat(fstr);
            strcpy(fstr,"%s"); i=2;
         }

         char *s0=va_arg(args, char*);

         if (!s0) { Sb.catf(fstr,"(null)"); } else {

         for (k=0; s0[k]>=30; ++k);

         if ((!s0[k] && k<4) || i>2) { Sb.catf(fstr,s0); }
         else {
            char ck=s0[k]; if (ck) { k+=strlen(s0+k); }
            char s1_[k+1]; char *s1=s1_; 
            strcpy(s1,s0);
            for (k=0; s1[k] || !k;) { 
               while (s1[k]) {
                  if (s1[k]=='\n') { 
                     Sb.catf("\n%s",log_header);
                     s1+=(k+1); k=0; 
                  }
                  else if (isspace(s1[k])) { ++k; }
                  else break;
               }

               k=(fwd ? 0 : check_update_header(log_header,s1,"> "));
               if (k) { s1+=k; k=0;
                  if (!ck && *fmt==' ') { ++fmt; }

                  Sb.catf("\r%s",log_header);
               }

               for (k=0; s1[k]; ++k) {
                  if (s1[k]=='\f') { s1[k]='\n'; } else
                  if (s1[k]=='\n') break;
               }

               if (s1[k]) { 
                  for (j=k+1; s1[j]; ++j) { if (s1[j]!='\n') break; }
                  if (!s1[j]) { Sb.catf("%s",s1); break; }

                  s1[j-1]=0; 
                  Sb.catf("%s\n",s1); 
                  Sb.cat(log_header);
                  s1+=j; k=0;
               }
               else { 
                  j=strlen(s1); if (j) { 
                  if (s1[j-1]=='\\') { bflag|=2; 
                     s1[j-1]=0;
                  }}
                  Sb.catf(fstr,s1); break;
               }
            }
         }}
         i=isfmt=0;
      }
      else if (strchr("dxXp",c)) {
         fstr[i]=0;
         m=0; if (i>1) if (fstr[i-2]=='l') m++;
         if (m) Sb.catf(fstr, va_arg(args, long));
         else   Sb.catf(fstr, va_arg(args, int ));

         i=isfmt=0;
      }
      else if (c=='D') {
         char s[nt]; int j,l,q;

         if (i>1 && fstr[i-2]=='l') {
            fstr[i-2]='s'; fstr[i-1]=0;
            l=snprintf(s,nt,"%ld",va_arg(args,long));
         }
         else {
            fstr[i-1]='s'; fstr[i]=0;
            l=snprintf(s,nt,"%d", va_arg(args,int ));
         }

         j=l+(l-(s[0]!='-' ? 1 : 2))/3;
         if (j>l && j<int(nt)) {
            s[j]=0; q=(l-1)%3; 
            for (--l, --j; l>0 && j>l; --j) {
               s[j]=s[l--]; if (l%3==q) { s[--j]=','; }
            }
         }
         Sb.catf(fstr,s);
         i=isfmt=0;
      }
      else if (c=='c') { 
         fstr[i]=0;
         Sb.catf(fstr, va_arg(args, int ));
         i=isfmt=0;
      }
      else if (strchr("eEgf",c)) { 
         fstr[i]=0;
         m=0; if (i>1) if (fstr[i-2]=='L') m++;
         if (m) Sb.catf(fstr, va_arg(args, long double));
         else   Sb.catf(fstr, va_arg(args, double));
         i=isfmt=0;
      }
      else if (c=='N') { 
         fstr[i-2]=0;
         Sb.catf("%s\n",fstr);
         i=isfmt=0;
      }
      else if (c=='Y') {
         fstr[i-2]=0; l=16;
         char c_; c_=time_stamp[l]; 
         time_stamp[l]=0; Sb.catf("%s%s",fstr,time_stamp);
         time_stamp[l]=c_;

         i=isfmt=0;
      }
      else if (c=='R') {
         k=i; fstr[--i]=0; 

         for (; i<flen; --i) { if (fstr[i]=='%') break; }

         if (i>=flen) {
            Sb.catf("%s:%d ERR invalid format '%s'\n", FLINE, fmt-k);
         }
         else {
            fstr[i]=0; Sb.catf("%s",fstr);

            if (sscanf(fstr+i+1,"%d",&k)<=0) {
               Sb.catf("%s:%d ERR invalid format '%s'\n", FLINE, fmt-k);
            }
            else {
               snprintf(fstr,flen,"%s",va_arg(args, char *));
               for (; k>0; --k) { Sb.catf("%s",fstr); }
            }
         }

         i=isfmt=0;
      }
      else if (c=='B') {

         k = i; fstr[--i]=0;  

         for (; int(i)>=0; i--)
         if (fstr[i]=='%') break;

         if (int(i)<0)
         Sb.catf("%s:%d ERR Invalid format >%s<", FLINE, fmt-k);
         else {
            fstr[i] = 0;
            Sb.catf("%s",fstr);  

            i = sscanf(fstr+i+1, "%d", &k);
            if (i<=0) k=8;    

            m = va_arg(args, int);   
            cp = (char*)&m;         

            if (k>0)
            for (i=l=0; i<((k-1)/8+1); i++) {
                for (cb=1, j=0; j<CHAR_BIT; ++j, cb<<=1) {
                    if (cp[i] & cb) Sb.catf("1");
                    else            Sb.catf("0");
                    if ((++l)>=k) break;
                }
            }
         }
         i=isfmt=0;
      }
      else {
         fstr[i]=0; i=isfmt=0;
         Sb.catf(" >>%s<< ERR-fmt [..%s]\n",fstr,fmt);
         if (xcol) { Sb.cat("\e[0m"); }
         return -1;
      }
   }
   fstr[i]=0;

   if (i && (fstr[i-1]=='\\')) { fstr[i-1]=0; bflag|=2; }
   if (!(bflag &2)) strcat(fstr,"\n");

   Sb.catf("%s",fstr);
   if (xcol) { Sb.cat("\e[0m"); }

   if (nesc && !xcol_) {
      Sb.skipEscCols();
   }

   if (eflag) { rval|=8; }
   if (wflag) { rval|=4; }
   if (iflag) { rval|=2; }; rval|=1;

   #ifdef DBG_GCX_LOCKS
    if (!Hflag && (eflag || wflag)) { LKF.flush(); }
   #endif

   if (eflag) { if (Hflag) { ++wblog::ERR_pending; } else {
      #ifdef DBSTOP
      {  unsigned n=wblog::myIO.size(); 
         if (n) {
            fprintf(stdout,"\n"
              "Got %d pending entr%s in I/O buffer:\n\n",n,n==1?"y":"ies");
            for (const wblog::stdio_buf &b : wblog::myIO) {
               b.print_stdout(stdout);
            }
         }
         Sb.print(FL,"DBSTOP"); fflush(stdout);
      }
      #endif
      throw Wb::LogException(ERR); 
   }}
#ifdef DBSTOP
   else if (wflag) {
      Sb.print(FL,"DBSTOP");
      dbstop(file ? file : __FILE__, file ? line : __LINE__);
   }
#endif

   if (!Hflag) { wblog::check_ERR_pending(); } 

   return rval;
};

void wblog::check_ERR_pending() {
   if (wblog::ERR_pending) { char s[32];
      int q=wblog::ERR_pending; wblog::ERR_pending=0;
      snprintf(s,32,"pending error%s (e=%d)",q>1? "s":"",q);
      ExitMsg(s); 
   }
};

unsigned wblog_check_tag(const char *fmt, const char *t0, char *tag) {

    const char *s;
    int i,k,l=strlen(t0); if (!l) { return 0; }

    s=strstr(fmt,t0); if (!s) { return 0; }
    k=s-fmt;

    if (isalnum(s[l])) { return 0; } 
    for (i=0; i<k; ++i) {            
       if (isspace(fmt[i])) { return 0; }
    }

    memcpy(tag,t0,l); tag[l]=0;

    l+=k; if (fmt[l]) { ++l; } 

    return l;
};

char wblog_findtoken(const char *istr, const char *tok, int maxoffset) {

    char f=1; const char *c=strstr(istr,tok);

    if (c==0) return 0;
    if (maxoffset>=0) if (c>istr+maxoffset) return 0;

    if (c>istr)
    if (isalnum(*(c-1))) f=0; 
    if (isalnum(*(c+3))) f=0; 

    return f;
}

void Wb::SBUF::flush(FILE *fid, char fflag) {

   if (!slen || !sbuf) { return; }

   if (l && sbuf[0]) {

     #ifdef MATLAB_MEX_FILE
      if (fid==stdout || fid==stderr) { int tid=omp_get_thread_num();

         if (tid!=Wb::my_caller_tid) {
            wblog::myIO.push_back(wblog::stdio_buf(*this,fid,tid)); 
         }
         else {
            wblog::myIO.clear(); 

            if (Wb::envDKT==1) { mexPrintf("%s",sbuf); }
            else {
               PRINTF("%s",sbuf);
               fflush(stdout);
            }
         }
      }
      else { fprintf(fid,"%s",sbuf); fflush(fid); } 
     #else
      { fprintf(fid,"%s",sbuf); fflush(fid); }
     #endif

      if (fflag && (fid==stdout || fid==stderr)) { doflush(); }
   }

   if (sbuf) { sbuf[0]=0; l=0; }
};

void Wb::SBUF::print(const char *F, int L, const char *istr) {

   if (!l) return;
   if (!slen || l>slen) { fprintf(stderr,
      "\n\n%s %s got l=%d/%d !?\n\n",shortFLT,l,slen);
      return;
   }

   unsigned n=l+64, j=(sbuf[0]=='\n' ? 1:0); char s[n];
   char *sj=sbuf+j, c=sbuf[l]; sbuf[l]=0; 

   if (istr) { snprintf(s,n, 
      "\n%s (%s in %s)\n%s", shortFL(F_L),istr,myname,sj);
   }
   else if (F) { snprintf(s,n,"\n%s\n%s",shortFL(F,L),sj); }
   else        { snprintf(s,n,"\rTST> `%s'",sj); }

   PRINTF("\e[3%dm%s\e[0m\n",j?1:5,s); 

   sbuf[l]=c;
};

void Wb::SBUF::increase_size(const char *F, int L, unsigned n) {

   if (!slen || !sbuf || l>=slen) {
      error_bounds(F_L,FCT,"");
      ExitMsg(""); 
   }
   else if (ref) { char s[64];
      snprintf(s,64,"\n\n%s "
        "ERR cannot increase size %d -> %d for ref=%d data\n\n",
         shortFL(F_L),slen,n,ref);
      ExitMsg(s);
   }
   else if (wblog::SLEN<=32) { char s[64];
      snprintf(s,64,"\n\n%s ERR got %p @ %d, L=%d !?\n\n",
         shortFL(F_L),sbuf,slen,wblog::SLEN);
      ExitMsg(s);
   }
   else if (n>=slen) {

      unsigned nx=(1+(n+1)/wblog::SLEN)*wblog::SLEN;
      char *sx = new char[nx];
      if (!sx) { char s[64]; snprintf(s,64,
         "\n\n%s ERR allocation error (n=%d)\n\n",shortFL(F_L),nx);
         ExitMsg(s);
      }
      sbuf[slen-1]=0; 
      strcpy(sx,sbuf); delete [] sbuf; sbuf=sx; slen=nx;

     #if 0
      nx=strlen(sbuf); 
      fprintf(stderr,"%s TST increasing sbuf[ %d -> %d ] having "
        "l=%d%+d (%d)\n",SHORT_FL,nx+1,slen, l,n-l,n);
      if (Wb::envVRB & 8) 
         fprintf(stderr,"%s' (%X)\n\n",sbuf,Wb::envVRB);
     #endif
   }
};

void Wb::SBUF::error_bounds(
   const char *F, int L, const char *fct, const char *fmt) {

   if (sbuf && slen) {
      sbuf[slen-1]=0; 
      fprintf(stderr,
        "\n\n%s  ERR %s() string out of bounds (l=%d%+ld/%d):\n\n%s\n\n",
         shortFL(F_L),fct?fct:"(fct)",l,fmt ? strlen(fmt):-1,slen,sbuf);
      if (fmt) fprintf(stderr,"  ERR fmt = '%s'\n\n",fmt);
      if (l>=slen) ExitMsg("");
   }
   else {
      char s[128]; snprintf(s,128,
        "\n\n%s ERR %s() got uninitialized %p @ l=%d%+ld/%d !?\n\n",
         shortFL(F_L),fct?fct:"(fct)",sbuf,l,fmt ? strlen(fmt):-1,slen);
      ExitMsg(s);
   }
};

int Wb::SBUF::cat(const char *s) {
   int n=0; 

   if (s && s[0]) {
      n=strlen(s); if (l+n>=slen) { increase_size(FL,l+n); }
      strcpy(sbuf+l,s); l+=n;
   }
   return n;
};

int Wb::SBUF::catf(const char *fmt, ...) {

   if (!fmt || !fmt[0]) { return 0; }

   va_list args; Wb::ARGV wd(&args); 
   va_start(args,fmt);

   unsigned iter=0, n=strlen(fmt)+32; 
   if (l+n>=slen) { increase_size(FL,l+n); }
   if (l>=slen) { 
      fprintf(stdout,"%s %s ERR got l=%d/%d (%+d)",shortFLT,l,slen,n);
      fflush(stdout);
   }

   while (1) {
      n=vsnprintf(sbuf+l,slen-l,fmt,args);
      if (l+n<slen) { break; } else increase_size(FL,l+n);

      va_end(args); va_start(args,fmt); 
      if (++iter>2) fprintf(stdout,"%s %s ERR got iter=%d",shortFLT,iter);
   }

   if (iter>1) {
      fprintf(stderr,"%s TST increasing sbuf[%d->%d] having "
        "l=%d%+d (%d)\n",SHORT_FL, wblog::SLEN, slen, l,n, l+n);
      if (Wb::envVRB & 8) 
         fprintf(stderr,"\n`%s' (%X)\n\n",sbuf,Wb::envVRB);
   }

   l+=n; return n;
};

int Wb::SBUF::skipEscCols() {
   int nesc=0, nskip=0; 
   if (!sbuf) { return (nesc=-1); }
   unsigned i=0, k=0, j=0;

   for (; sbuf[i]; ++i) {
      if (sbuf[i]=='\e' && sbuf[i+1]=='[') { j=i+2; } else
      if (sbuf[i]=='\\' && sbuf[i+1]=='e' && sbuf[i+2]=='[') { j=j+3; }
      if (!j) { if (k<i) { sbuf[k]=sbuf[i]; }; ++k; }
      else {
         for (; sbuf[j]; ++j) {
            if (!isdigit(sbuf[j]) && sbuf[j]!=';') { break; }
         }
         if (sbuf[j]=='m')
              { ++nesc; i=j; }
         else { ++nskip; if (k<i) { sbuf[k]=sbuf[i]; }; ++k; }
         j=0;
      }
   }

   if (k<i) { sbuf[k]=0; l=k; }

   if (!nesc) { nesc=-nskip; }
   return nesc;
};

char* Wb::surdStrf(wbstring &s, const char *fmt, ...) {
   va_list args; Wb::ARGV wd(&args); 
   va_start(args,fmt);
   if (!s) { s.init(32); }
   return surdStrf(s.data,s.len,fmt,args);
};

char* Wb::surdStrf(char *s, unsigned n, const char *fmt, va_list args) {

   if (!s || n<8) wblog(FL,"ERR %s() invalid usage (s=%g, n=%d)",FCT,s,n);
   if (!fmt || !fmt[0]) { s[0]=0; return s; }

   unsigned i,k,l;

   k=std::wcrtomb(s,u'\u221A',NULL); 

   i=l=k+1; 
   l+=vsnprintf(s+l,n-l,fmt,args);

   for (; i<n && s[i]; ++i) { if (s[i]<'0' || s[i]>'9') break; }

   if (l+(s[i] ? 1:0)>=n) wblog(FL,
      "WRN %s() string out of bounds (%d/%d)\n'%s'",FCT,l,n,s);

   if (s[i]) {
      s[k]='('; s[l]=')'; i=0; 
   }
   else {
      for (i=k; i; --i) { s[i]=s[i-1]; }
      s[i++]=' '; 
   }
   return s+i;
};

Wb::termcolor& Wb::termcolor::init(unsigned k) {

   if (em) { delete [] em; em=e1=NULL; }

   if (!Wb::useCol) { em=e1 = new char[1]; em[0]=0; }
   else {
      unsigned l, n=(k<8 ? 12 : 18);
      em = new char[n];
      l=snprintf(em,n,"\e[0m"); 
      e1 = em + (++l);

      if (k<  8) { l+=snprintf(e1,n-l,"\e[3%dm",k); } else 
      if (k<256) { l+=snprintf(e1,n-l,"\e[38;5;%dm",k); }  
      else wblog(FL,"ERR %s() color index out of bounds (%d)",FCT,k);

      if (l>=n) wblog(FL,
      "WRN %s() termcolor() string out of bounds (%d/%d)",FCT,l,n);
   }

   return *this;
};

int Wb::termcolor::ID(const char *tag) {

   int i=0; 

   if (Wb::useCol && tag && tag[0]) {
   switch (tag[0]) {
    case 'h':  
       if (!strcmp(tag,"hdr")) { i=240; } 
       break;
    case 'E':
       if (!strcmp(tag,"ERR")) { i=  1; } else 
       if (!strcmp(tag,"ENV")) { i=240; } 
       break;
    case 'W':
       if (!strcmp(tag,"WRN")) { i=  9; } 
       break;
    case ' ':
       if (!strcmp(tag,"  *")) { i=243; } else 
       if (!strcmp(tag," * ")) { i=246; }      
       break;
    case 'T':
       if (!strcmp(tag,"TST")) { i=246; } 
       break;
    case 'N':
       if (!strcmp(tag,"NB!")) { i= 34; } 
       break;
    case 'D':
        if (!strcmp(tag,"DBG")) { i=130; } 
        break;
    case 'X':
       if (tag[1]=='X') {
          if (tag[2]=='E') { i=1; } else 
          if (tag[2]=='W') { i=9; } else 
          if (tag[2]=='G') { i=2; } else 
          if (tag[2]=='Y') { i=3; } else 
          if (tag[2]=='B') { i=4; }      
       }
       break;
   }}

   return i;
};

#endif

