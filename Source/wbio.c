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

#ifndef __WB_WBIO_C__
#define __WB_WBIO_C__

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
// generate temporary file $TMPDIR/<p=myprog><f="">-<pid>.<x=log>
// in $TMPDIR (must be defined in environment)
// Wb,Apr30,20

wbstring Wb::tmpfile(const char *F, int L,
   const char *t, 
   const char *x, 
   const char *f  
) {

   wbstring fout(64); 

   unsigned l; const char *p=getenv("TMPDIR");
   if (!p) wblog(F_L,"ERR TMPFILE not defined (null)"); else
   if (p[0]!='/') wblog(F_L,"ERR invalid TMPFILE `%s'",p);

   l=snprintf(fout.data,fout.len,"%s/%s%s-%d%s", p,
      f && f[0] ? f : myname, t? t:"", getpid(), x? x : ".log");

   if (l>=fout.len) wblog(FL,
      "ERR %s() string out of bounds (%d/%d)",FCT,l,fout.len);

   return fout;
};

void Wb::tmpFile::init(const char *F, int L,
   const char *f,   
   const char *x,   
   const char *mode 
) {

   if (fid) wblog(FL,
      "ERR %s() alraedy got initialized fid=%d",FCT,fid);

   fout=tmpfile(F_L,f,x);
   fid=fopen(fout.data, mode && mode[0] ? mode : "a+");

   if (!fid) wblog(F_L,
      "ERR %s() failed to open file\n`%s'",FCT,fout.data);
   else {
      fprintf(fid,"\n%-18s %5s %s open %s\n\n", 
      shortFL(PF_L),"", Wb::TimeStamp('T').data, Wb::TimeStamp('D').data);
   }

   if (F) wblog(PF_L,"I/O opening tmpfile `%s'",fout.data);
};

int Wb::tmpFile::sepline(const char *x, unsigned n) {
   unsigned l=0;

   if (!n || !ncall) { return l; }
   if (int(n)< 0) { n=72; } else 
   if (int(n)>90) { n=90; } 

   unsigned m=(x ? strlen(x) : 0);
   wbvec<char> s_(n+1); char *s=s_.data;

   if (m<=1) { memset(s, m? x[0] :'-', n); }
   else {
      while (l+m<=n) { memcpy(s+l,x,m); l+=m; }
      if (l<n) { memcpy(s+l,x,n-l); }
   }

   if (!omp_in_parallel()) { this->flush(); }

   l=n; s[l]=0; {
      if (!fid) wblog(FL,"ERR %s() tmpFile not yet opened !?",FCT);
      #pragma omp critical (__tmpFile_IO__)
      { fprintf(fid,"%s\n",s); }
   }

   ncall=0;

   return l;
};

int Wb::tmpFile::blogf(const char *F, int L, const char *fmt, ...) {

   va_list args; Wb::ARGV wd(&args); 
   va_start(args,fmt); unsigned l=0;

   char sx[8];
   snprintf(sx,8,"%d/%d",omp_get_thread_num(), omp_get_num_threads());

   if ((++ncall)<1) { ncall=1; }

   if (fid) {
     #pragma omp critical (__tmpFile_IO__)
      { l+=fprintf(fid,"%-18s %5s %s ", 
           shortFL(F_L,18), sx, Wb::TimeStamp('T').data);
        l+=vfprintf(fid,fmt,args);
        l+=fprintf(fid,"\n");
      }
   }
   else {
      fprintf(stdout,"\n%s %s " 
      "ERR %s() tmpFile closed / not yet opened !?\n\n",shortFLT, FCT);
   }

   return l;
};

void WbPrintMatrixC (
    const wbMatrix<wbcomplex> &M, const char *istr, int space
){
    unsigned i,j;
    wbvec<char> dstr(64), sfmt(64);

    int rflag, iflag;

    PRINTF("%s%s", istr, istr[0] ? "\n" : "");

    if (space<1) sfmt.catf(FL," %%s");
    else         sfmt.catf(FL," %%%ds", space);

    for (i=0; i<M.dim1; i++) {
        for (j=0; j<M.dim2; j++) {

            rflag = std::isnan(M(i,j).r);
            if (!rflag) {
                if (M(i,j).r!=0.) {
                    rflag = (fabs(M(i,j).r)==1.) ? 1 : 2;
                    if (M(i,j).r<0) rflag = -rflag;
                }
            }

            iflag = std::isnan(M(i,j).i);
            if (!iflag) {
                if (M(i,j).i!=0.) {
                    iflag = (fabs(M(i,j).i)==1.) ? 1 : 2;
                    if (M(i,j).i<0) iflag = -iflag;
                }
            }

            if (rflag) {
                dstr.catf(0,0,"% g",M(i,j).r);
                if (iflag) {
                   if (abs(iflag)>1)
                        dstr.catf(0,0,"%+gi", M(i,j).i);
                   else dstr.catf(0,0,"%ci",  iflag==1? '+':'-');
                }
            }
            else {
                if (iflag)
                   if (abs(iflag)>1)
                        dstr.catf(0,0,"% gi", M(i,j).i);
                   else dstr.catf(0,0,"%ci",  iflag==1? ' ':'-');
                else {  dstr.catf(0,0,"%g", M(i,j).r); }
            }
            PRINTF(sfmt.data,dstr.data);
        }
        PRINTF("\n");
    }
    PRINTF("\n");
};

template<class T>
int Wb::Str2Idx(
   const char *F, int L, 
   const char* s_, wbvector<T> &idx,
   T offset,             
   char cmpct,           
   unsigned l_           
){
   if (!s_) {
      wblog(F_L,"ERR %s() got null input !?",FCT);
      return -1;
   }
   if (!s_[0]) { idx.init(); return  0; } 

   unsigned i,j, n=0,
      ma=0, 
      mx=0, 
      ms=1; 

   char c, sx=0, bflag=0;
   long x; int e=0;

   unsigned l=strlen(s_);
      if (l>l_) { l=l_; } 
      else      { l_=l; }

   wbvec<char> S_(l+1,s_); 
   char *s=S_.data, *s1=s, *s2=s;

   for (i=0; i<l && !e; ++i) { c=s[i]; if (!isspace(c)) {
      if (isalnum(c)) { ++ma; ms=0;
         if ((++n)>1 && !sx) { sx=','; }
         for (++i; i<l; ++i) {
            if (isalnum(s[i])) { ++ma; } else { --i; break; }}
         continue;
      }
      else { ++mx; s[i]=' '; }

      if (c==',' || c==';') { ++ms;
         if (n && ms==1)       
              { if (!sx) { sx=c; } else if (sx!=c) e=+__LINE__; }
         else { e=+__LINE__; } 
      }
      else if (c=='[') { if (!n && !bflag) { bflag|=1; } else { e=-__LINE__; }}
      else if (c==']') { if (bflag&2) { e=-__LINE__; } else {
         bflag|=2; ++i;
         for (; i<l && s[i]; ++i) { if (!isspace(s[i])) { e=__LINE__; break; }}
      }}
      else { e=__LINE__; } 
   }}

   if (!n && !mx) { idx.init(); return 0; }

   if (bflag && bflag!=3) {  
      if (e) { if (e>0) e=-e; }
      else { e=-__LINE__; }  
   }

   if (e<0) { if (F) wblog(F,L,  
      "ERR %s() invalid '%s' (i=%d/%d; b=%d)",FCT,s_,i,l,bflag);
      return -i;
   }

   if (e) { l=i; } 

   if (l>l_) wblog(F_L,"ERR %s() '%s' (e=%d, l=%d/%d)",FCT,s_,e,l,l_);

   if (n==1 && cmpct) { 
      n=ma; idx.init(n); j=-1;
      for (i=0; i<l; ++i) {
         if (isalnum(s[i])) { c=s[i];
            if ((++j)>=n) wblog(FL,"ERR %s() j=%d/%d !?",FCT,j,n);
            if (c>='0' && c<='9') { idx[j]=c-'0';    } else
            if (c>='A' && c<='Z') { idx[j]=c-'A'+10; } else
            if (c>='a' && c<='z') { idx[j]=c-'a'+36; }
            else { e=-__LINE__; ++i; break; } 
         }
      }
   }
   else { idx.init(n); 
      for (j=0; j<idx.len; ++j) { s1=s2;
         x=strtol(s1,&s2,0); if (s2==s1) { break; }
         idx[j]=T(x);
      }
      if (j!=idx.len || (s2-s)>l_) { e=-__LINE__; }
      if (!e) {
         for (i=s2-s; i<l_; ++i) {
             if (!isspace(s[i])) { l=i; e=__LINE__; }
         }
      }
   }

   if (e<0) { if (F) wblog(F,L,  
      "ERR %s() invalid '%s' (%d/%d;%d)",FCT,s_,i,l,bflag);
      return -l;
   }

   if (offset && e>=int(0)) {
      for (i=0; i<n; ++i) {
      if (idx[i]>=offset) { idx[i]-=offset; } else { e=-__LINE__; }}
   }

   if (e<0) {  
      if (F) wblog(F,L,"ERR %s() invalid '%s' "
         "(e=%d; %d/%d, %d @ %d)",FCT,s_,-e,idx.len,l,bflag,offset);
      return -l;
   }

   if (e) { 
      if (F) wblog(F,L,"WRN %s() trailing string in "
        "'%s' [e=%d, l=%d, %d/%d]",FCT,s_,e,l,bflag,cmpct);
      return -l; 
   }

   return idx.len;
};

template<class T>
int Wb::Str2Idx(
   const char* F, int L, 
   wbvector<T> &I, const mxArray *a, T offset, char cmpct
){
   if (!a) wblog(FL,"ERR %s() got null mxArray",FCT);

   if (mxIsChar(a)) { wbstring s(a);
      int n=Wb::Str2Idx(0,0,s.data,I,offset,cmpct);
      if (n<0) wblog(F_L,
         "ERR %s() '%s' => [%s] @ e=%d !?",FCT,s.data,STR(I),n);
      return I.len; 
   }
   else {
      I.init(F,L,a);
      if (offset) {
         if (I.anyLT(offset)) wblog(FL,
            "ERR %s() index incompatible with offset=%d",FCT,offset);
         I-=offset;
      }
      return 0;
   }
};

#endif

