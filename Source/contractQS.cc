
/* Additional options (undocumented) ================================= *

    --cg-preview // Wb,Oct01,16

      do not perform any contraction just print new contractions
      to be performed that are not in the current RCStore yet
      (together with usage #1 only)

   Usage #3 // Wb,Sep12,16

      similar usage as contractRC, yet with no & separating the
      contraction sequences

      contractQS --rcs sym A1 iai1 B1 ib1, A2 iai2 B2 ib2, ...

   Usage #4 // Wb,Sep24,16

      low-level CData contraction
      => return result of single contraction

      contractQS --rcd sym A1 iai1 B1 ib1

      Ex. contractQS('--rcd','SU4','(036,121,521)','123*','(036,121,521)','123')
      --> returns overlap in OM space [here eye(3)]

   Usage #5 // Wb,Sep26,16

      display current internal status of loaded/written RCStore data
      contractQS --stat

      return current internal data
      [ gCS, gXS, gRS ] = contractQS --getCG
*/

/* CHANGE LOG ======================================================== *

// Wb,Apr19,13 : increased functionality of cell-contraction
// also introducing the complex conjugation tags.
*/

char USAGE[]=""; 

#ifdef MATLAB_MEX_FILE
   #define PROG mexFunctionName()
#else
   #define PROG "contractQS"
#endif

#define PROG_TAG "ctr"

#define LD_CLEBSCH_QS
#include "wblib.h"

template<class TA, class TB, class TC>
int CONTRACT_QS(
    const QSpace<gTQ,TA> &A, const ctrIdx &ica,
    const QSpace<gTQ,TB> &B, const ctrIdx &icb,
    QSpace<gTQ,TC> &C, wbperm P,
    char cg_preview=0, 
    char vflag=0
);

int is_valid_cell_ctr(
    unsigned na, const mxArray *aa[],
    unsigned level,
    unsigned &lmax  
);

class icFlags { 
 public:

   icFlags(const mxArray* ain=NULL) : a(ain), conj(0), xflag(0) {
      checkQSpace(a); };

  ~icFlags() {}; 

   icFlags& checkQSpace(const mxArray *a);

   char check_arg(const mxArray *ain);

   icFlags& init(const mxArray *ain) {
      if (!ain) wblog(FL,"ERR %s() got NULL QSpace or cell !?",FCT);
      a=ain; checkQSpace(a);
      sidx.init(); otags.init(); regex_r.init(); xflag=conj=0;
      return *this;
   };

   char set(const mxArray *a, const char *istr="");

   template <class TQ, class TD>
   void apply(QSpace<TQ,TD> &X, char vflag=1);

   const mxArray *a;  

   char conj;         
   char xflag;        

   wbstring otags;    
   wbstring regex_r;  
   wbstring sidx;     

 protected:
 private:

};

char icFlags::check_arg(const mxArray *ain) {

   char q=0; 
   a=ain; 
   if (a) {
      if (mxIsChar(a) || mxIsNumeric(a)) {
         int i=set(a,0); 
         if (i>0) { q=1;
            if (i==2) { q|=2; }  
            else if (i!=1) wblog(FL,"ERR %s() invalid i=%d",FCT,i);
         }
         else if (!i) { q|=3; }  
      }
      else if (mxIsCell(a)) { q=12; }  
      else {
         int i=mxIsQSpace(0,0,a,'c');
         if (i>0) { q=4; } 
         else if (i<-1) {
            int n=mxGetNumberOfElements(a); 
            if (n>1) wblog(FL,
               "ERR invalid input QSpace (array with %d elements)",n);
            else wblog(FL,"ERR invalid input QSpace (n=%d, e=%d)",n,i);
         }
      }
   }
   return q;
};

icFlags& icFlags::checkQSpace(const mxArray *a) {

   if (a) {
      if (mxIsCell(a)) {
         if (mxGetNumberOfElements(a)<2) wblog(FL,
           "ERR %s() invalid cell-contraction\n"
           "(got cell with fewer than two entries (%d)!)",
            FCT, mxGetNumberOfElements(a)
         );
      }
      else {
         if (mxIsQSpace(FL,a,'c')<=0) wblog(FL,
           "ERR %s() invalid QSpace in cell-contraction (%s)",
            FCT, mxGetClassName(a));
         if (mxGetNumberOfElements(a)!=1) wblog(FL,
           "ERR %s() invalid QSpace(%d) in cell-contraction",
            FCT, mxGetNumberOfElements(a)
         );
      }
   }
   return *this;
};

char icFlags::set(const mxArray *a, const char *istr) {

   char q=0; 
   char tflag=(istr ? 0 : 1);

   if (mxIsChar(a)) {
      wbstring S(FL,a); const char *s=S.data;
      if (!s || !s[0]) {
         if (!tflag) wblog(FL,"ERR %s() %s\n"
           "invalid string option for cell-contraction (%s)",
            FCT, istr, s? "empty":"null");
         else { q=-1; return q; }
      }

      unsigned i=0; while (s[i]=='-') { ++i; }

      if ((i==1 || i==2) && !strncasecmp(s+i,"op:",3)) { if (!tflag) {
         if (otags) wblog(FL, 
            "ERR %s() otags already set ( %s / %s)",FCT,otags.data,s+1);
         otags=s+1;
      }}
      else if (!strncmp(s,"--itag:s/",9)) { if (!tflag) { 
         if (regex_r) wblog(FL,
            "ERR %s() regex_r already set (%s / %s)",FCT,regex_r.data,s);
         regex_r=s+7; 
      }}
      else { 
         unsigned l=0;
         if (s[0]=='!') { ++s; xflag=0; } else { xflag=1; }
         for (; s[l]; ++l) { if (!isdigit(s[l])) break; }

         if (s[l]) {
            if (s[l]==CC_ITAG && !s[l+1]) { if (!tflag) {
               if (conj) wblog(FL, 
                  "ERR %s() conj alrady set",FCT);
               conj=1;
            }}
            else if (tflag) { q=-2; return q; }
            else { wblog(FL,
               "ERR %s() %s\ninvalid cell-contraction option\n"
               "(expecting '[!##][*]': got %s)",FCT,istr,S.data);
            }
         }
         if (l && !tflag) {
            if (sidx) wblog(FL, 
               "ERR %s() sidx already set ( %s / %s)",FCT,sidx.data,s);
            sidx.init(s,l);
         }
      }
      q=1; 
   }
   else if (mxIsNumeric(a)) {
      if (Mx::IsVector(a)) {
         unsigned i=0, j;
         wbvector<double> x(FL,a); 
         if (x.len && !tflag) {
            if (sidx) wblog(FL, 
               "ERR ctr-idx already set ( %s / [%s])",sidx.data,STR(x));
            sidx.init(x.len+1);
            for (; i<x.len; ++i) { j=x[i];
               if (!j || j>9 || j!=x[i]) { wblog(FL,
                  "ERR invalid ctr-idx [%s]\nfailed to convert "
                  "to compact string (j=%d/%g)",STR(x),j,x[i]);
               }
               sidx[i]=j+'0'; 
            }
            sidx[i]=0; xflag=1;
         }
         q=(x.len ? 2 : 0);
      }
      else if (mxGetNumberOfElements(a)) { q=-3; }
   }
   else {
       q=(mxIsCell(a) || mxIsQSpace(FL,a,'c') ? -16 : -4);
       if (!tflag) { wblog(FL,
          "ERR %s() unexpected input %s",FCT,mxGetClassName(a));
       }
   }

   return q;
};

template <class TQ, class TD>
void icFlags::apply(QSpace<TQ,TD> &X, char vflag) {

   if (X.isEmpty()) { return; } 
   unsigned r=-1;

   if (otags) { 
      if (!strncmp(otags.data,"op:",3)) {  
         int q=X.isOperator(&r);
         if (q>0) { r=q; } 
         else { wblog(FL,  
            "ERR failed to set itags for rank%+d QSpace\n"
            "'%s' %s ('%s', q=%d/%d, e=%d)",
            q, X.itags2Str().data, X.otype2Str().data, otags.data,
            X.itags.isOp(), X.itags.got_op_labels(), q);
         }
      }
      X.init_itags(FL,otags.data); 
   }
   else { r=X.rank(FL); }

   if (regex_r && vflag) {
      if (!X.itags.RegEx_replace(FL,regex_r.data)) { wblog(FL,
         "WRN %s() regex_r = '%s' had no effect",FCT,regex_r.data);
         mexWRN("irrelevant regex_r"); 
      }
   }

   if (sidx) {
      unsigned i=0, j, r_=(r ? r : 9); const char *s=sidx.data;
      wbstring mark(r_); char *m=mark.data;

      for (; i<sidx.len && s[i]; ++i) { j=s[i]-'1';
         if (j>=r_) wblog(FL,
            "ERR %s() cell-ctr index out of bounds (%s; r=%d)",FCT,s,r_);
         if (++m[j]!=1) { wblog(FL,
            "ERR %s() invalid cell-ctr index %s (non-unique)",FCT,s);
         }
      }

      if (r!=X.itags.len) wblog(FL,"ERR %s() "
         "valid set of itags required (%s; r=%d)",FCT,IT2STR(X),r);

      for (i=0; i<r; ++i) {
         if ((!xflag) ^ (!m[i])) X.SetFlags(FL,i+1); 
      }
   }
};

mxArray* contract_plain(int nargin, const mxArray *argin[]);

template<class TD>
unsigned contractQS_itags(
   int nargin, const mxArray *argin[],
   unsigned level, unsigned &vflag, QSpace<gTQ,TD> &C,
   wbperm &P 
);

int contractRC(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[]
);

int contractCD(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[]
);

#ifdef DBG_GCX_LOCKS
   unsigned ncall_contract=0;
#endif

void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){ Wb::CleanUp aclu; try { 
    MX_CHECK_HELPER_NARGS(1,-1,3); 

    if (nargin==1 && !mxIsCell(argin[0])) {  
       if (Mx::IsEqual(argin[0],"--stat")) { 
          CG::MemStat(FL);
          if (nargout) wblog(FL,
             "ERR invalid usage (no output with option --stat)",FCT);
          return;
       }
       else if (Mx::IsEqual(argin[0],"--getCG")) {
          argout[0]=gCS.toMx(); if (nargout>1) {
          argout[1]=gXS.toMx(); if (nargout>2) {
          argout[2]=gRS.toMx(); }}
          return;
       }
       else wblog(FL,"ERR %s() invalid usage",FCT);
    }

    if (mxIsChar(argin[0])) {
       if (Mx::IsEqual(argin[0],"--rcs")) { 
          contractRC(nargout,argout,nargin-1,argin+1);
          return;
       }
       if (Mx::IsEqual(argin[0],"--rcd")) { 
          contractCD(nargout,argout,nargin-1,argin+1);
          return;
       }
    }

    get_CG_VERBOSE(FL);
    get_CG_FIXIT(FL);

#ifdef __WBDEBUG__
    Wb::MemCheck(FL,"start");
#endif
#ifdef DBG_GCX_LOCKS
    LKF.sepline(); 
    LKF.blogf(FL,"call #%d",++ncall_contract);
#endif
#if defined(DBG_QSX_BUF)
    BFF.sepline(); 
    BFF.blogf(FL,"--- %s() having BUF[%ld], XBUF[%ld]",
      FCT, gCS.BUF.size(), gXS.XBUF.size());
#endif

    if (nargin<1) wblog(FL,
       "ERR invalid number of input arguments (%d)",nargin);
    if (nargout>1) wblog(FL,"ERR invalid number of output arguments");

    mxArray *a=NULL;

    itag_::Reset(2); 

    if (nargin>=4 &&  
         !mxIsCell(argin[0]) && isCtrIdx(argin[1]) &&
         !mxIsCell(argin[3]) && isCtrIdx(argin[3])
       ) { 
       a=contract_plain(nargin,argin);
    }
    else {
       unsigned i,l=0, lmax=0, zflag=0,
          nc=nargin,
          level=0;  
       unsigned vflag=1; wbperm P;
       char q, nQS=0;

       icFlags Q;
       for (i=0; i<(unsigned)nargin; ++i) {
          q=Q.check_arg(argin[i]); if (q<=0) { nc=i; break; }
          if (q>3) { l=0; ++nQS; } 
          else if (q&2) { l=i; }   
       }

       if (!nQS) {
          if (i<(unsigned)nargin)
               { sprintf_str("arg #%d",i+1); }
          else { sprintf_str("QSpace/cell pair required; n=%d",nQS); }
          wblog(FL,"ERR %s() invalid usage (%s)",myname,str);
       }

       if (nc<(unsigned)nargin) {
          OPTS opts(argin+nc,nargin-nc);

          if (opts.getOpt("-v")) { vflag=2; } else
          if (opts.getOpt("-q")) { vflag=0; } 

          opts.checkAnyLeft(FL);
       }

       if (l+1==nc) { 
          int r=isValidPerm(argin[l]); 
          if (r>=0 && (r!=1 || l<2)) {
             P.init(FL,argin[l],1);  
             --nc; 
          }
       }

       if (int(i=is_valid_cell_ctr(nc,argin,0,lmax))<=0) wblog(FL,
          "ERR invalid cell contraction (e=%d)",-i);
       zflag=(i>>16); 

       if (vflag>1) {
          unsigned l,n=64, nQS=((i>>8)&255) + zflag; char s[n];
          l=snprintf(s,n,"%d level%s total",lmax+1, lmax?"s":"");
          if (zflag && l<n) { l+=snprintf(s+l,n-l,
             ", %d/%d QSpace%s complex",zflag, nQS, zflag!=1? "s":""); }
          wblog(FL,"=== cell contraction %34R","=");
          wblog(FL," *  %d contractions total\n%s",i&255,s);
       }

       if (zflag) {
          QSpace<gTQ,wbcomplex> C;
          contractQS_itags(nc,argin,level,vflag,C,P);
          a=C.toMx();
       }
       else {
          QSpace<gTQ,double> C;
          contractQS_itags(nc,argin,level,vflag,C,P);
          a=C.toMx();
       }
    }

    argout[0]=a; 

#ifdef DBG_GCX_LOCKS
    LKF.flush();
#endif
#if defined(DBG_QSX_BUF)
    BFF.flush();
#endif

#ifdef __WBDEBUG__
    Wb::MemCheck(FL,"stop");
#endif

}  catch (Wb::LogException &e) { ExitMsg(e.istr); }
   catch (...) { ExitMsg("caught exception in contractQS"); }
   aclu.Check();
};

int match_regex_itag(
   wbstring &otags, const iTags &itags,
   const char *Astr, const char *Bstr 
){
   int rval=0;

   if (!itags.len || !otags.strchr("[]\\^$")) {
      return rval;
   }

   const char *s=otags.data;

   if (!strncasecmp(s,"op:",3)) { 
      int i=-1, k=0;

      if (s[3]=='\\' && !s[5] && s[4]>='1' && s[4]<='9') {
         i=s[4]-'0';
      }
      else {
         for (unsigned j=4; s[j]; ++j) { if (s[j]==':') k=j; }
         if (k) { ((char*)s)[k]=0; } 
         i=itags.findRegEx(s+3);
         if (k && i<=0) { ((char*)s)[k]=':'; } 
      }

      if (unsigned(i)>itags.len) wblog(FL,
         "ERR %s() index out of bounds (%d/%d)",FCT,i,itags.len);

      if (i>0) { unsigned n=16; char sx[n]; 
         unsigned l=3; strncpy(sx,s,l);     

         l+=snprintf(sx+l,n-l,"%s",itags[i-1].toStr(0).data);
         if (k) { l+=snprintf(sx+l,n-l,":%s",s+k+1); }
         otags=sx; 
      }
      else if (!i) wblog(FL,"ERR auto-contract failed to match regex \n"
        "`%s' for %s in %s: `%s'",s+3,Astr,Bstr,STR(itags));
      else wblog(FL,"ERR got multiple (%d) itag matches for\n"
        "regex `%s' for %s in %s: `%s'",-i,s,Astr,Bstr,STR(itags));

      if (i>0) { rval=i; }
   }
   else wblog(FL,"WRN ignoring itag regex '%s' !?",s);

   return rval;
};

mxArray* contract_plain(int nargin, const mxArray *argin[]) {

   char isra=1, isrb=1, vflag=1, cg_preview=0;
   unsigned i=0, l=4, k=2;

   ctrIdx ica, icb;
   mxArray *a=NULL;

   OPTS opts;
   wbperm P;

   if (nargin<int(l)) wblog(FL,"ERR %s() "
      "invalid usage #1 (nargin=%d/%d)",PROG,nargin,l);
   if (!isCtrIdx(argin[1]) || !isCtrIdx(argin[3])) wblog(FL,"ERR %s() "
      "invalid contraction index\n(args #1 or #3 in usage #1)",PROG);

   ica.init(FL,argin[1]);
   icb.init(FL,argin[3]); 

   if (nargin>(int)l) {
      if (!mxIsChar(argin[l])) {
         P.init(FL,argin[l++],1); 
      }
      else { 
         char pstr[32];
         if (!mxGetString(argin[l],pstr,30)) {
            if (Wb::isUIntString(pstr)>0) { ++l;
               P.initStr(FL,pstr,1); 
            }
         }
      }
   }

   opts.init(argin+l,nargin-l);

   if (opts.getOpt("conjA")) {
      wblog(FL,"WRN %s() got deprecated option 'conjA'",FCT);
          mexWRN("deprecated conjA flag");
      ica.Conj(); ++i;
   }
   if (opts.getOpt("conjB")) {
      wblog(FL,"WRN %s() got deprecated option 'conjB'",FCT);
          mexWRN("deprecated conjB flag");
      icb.Conj(); ++i;
   }

   if (opts.getOpt("--cg-preview")) cg_preview=1; 
   if (opts.getOpt("-v")) vflag=2; else 
   if (opts.getOpt("-q")) vflag=0; 

   opts.checkAnyLeft(FL);

   for (i=0; i<2; ++i) { unsigned j=(i==1 ? k:i);
      try { mxIsQSpace(FL,argin[j],'c'); }
      catch (...) {
        wblog(FL,"ERR %s() invalid QSpace at arg #%d (%s)",
        FCT,i+1,mxGetClassName(argin[j]));
      }
   }

   if (!(isra=(mxIsQSpace(argin[0])>0))) { mxIsQSpace(FL,argin[0],'c'); }
   if (!(isrb=(mxIsQSpace(argin[k])>0))) { mxIsQSpace(FL,argin[k],'c'); }

   if (isra) {
      QSpace<gTQ,double> A(argin[0],'r');

      if (isrb) {
         QSpace<gTQ,double> B(argin[k],'r'); QSpace<gTQ,double> C;
         CONTRACT_QS(A,ica,B,icb,C,P,cg_preview,vflag); a=C.toMx();
      }
      else {
         QSpace<gTQ,wbcomplex> B(argin[k],'r'), C;
         CONTRACT_QS(A,ica,B,icb,C,P,cg_preview,vflag); a=C.toMx();
      }
   }
   else {
      QSpace<gTQ,wbcomplex> A(argin[0],'r');

      if (isrb) {
         QSpace<gTQ,double> B(argin[k],'r'); QSpace<gTQ,wbcomplex> C;
         CONTRACT_QS(A,ica,B,icb,C,P,cg_preview,vflag); a=C.toMx();
      }
      else {
         QSpace<gTQ,wbcomplex> B(argin[k],'r'), C;
         CONTRACT_QS(A,ica,B,icb,C,P,cg_preview,vflag); a=C.toMx();
      }
   }

   return a;
};

template<class TA, class TB, class TC>
int CONTRACT_QS(
    const QSpace<gTQ,TA> &A, const ctrIdx &ica,
    const QSpace<gTQ,TB> &B, const ctrIdx &icb,
    QSpace<gTQ,TC> &C, wbperm P,
    char cg_preview, 
    char vflag
){
    if (A.isEmpty() || B.isEmpty()) { C.init(); return 1; }

    int q=0;
    unsigned i, ra=A.rank(FL), rb=B.rank(FL); 

    for (i=0; i<ica.len; ++i) if (ica[i]>=ra) wblog(FL,
       "ERR index out of bounds ([%s]/%d)",STR(ica),ra);
    for (i=0; i<icb.len; ++i) if (icb[i]>=rb) wblog(FL,
       "ERR index out of bounds ([%s]/%d)",STR(icb),rb);

    if (!A.qtype.permitsOM(ra+rb-ica.len-icb.len)) {
       C.init(); C.mt=Wb::MEX_RETURN; 
    }

    if (vflag>1) { wblog(FL," *  plain pair-wise contraction\n"
    "   %-15s -> %s\n"
    "X  %-15s -> %s",IT2STR(A), STR(ica), IT2STR(B), STR(icb));
    }

    q=A.contract(FL,ica,B,icb,C,P,cg_preview); 
    if (cg_preview) { C.init(); return q; }

    C.SkipZeroData(); 
    C.NormCGW(); 

    C.ctime=Wb::getTimeNow();

    return q;
};

template<class TD>
unsigned contractQS_itags( 
   int nargin, const mxArray *argin[], unsigned level,
   unsigned &vflag,   
   QSpace<gTQ,TD> &C, 
   wbperm &P    
){
   unsigned i, l=-1, len=0; char mark[nargin];
   icFlags q;

   for (i=0; i<unsigned(nargin); ++i) {
      mark[i]=q.check_arg(argin[i]);
      if (mark[i]>=4) { ++len; } else
      if (mark[i]<=0 || !len) { 
         wblog(FL,"ERR invalid QSpace cell-structure (level %d @ %d)\n"
         "hint usage: ({QSpace} [,string opts or index])", FCT,level,i);
      }
   }

   wbvector<icFlags> Q(len);
   wbvector<QSpace<gTQ,TD> > X(len);

   for (i=0; i<(unsigned)nargin; ++i) {
      if (mark[i]>=4) { 
         icFlags &q=Q[++l];
         q.init(argin[i]); 

         if (mark[i]==4) { X[l].init(FL,argin[i],'r'); }
         else {
            if (!mxIsCell(q.a)) wblog(FL,"ERR %s() got non-cell !?",FCT);
            unsigned j=0, nc=mxGetNumberOfElements(q.a);
            const mxArray* ac[nc]; wbperm P_;
            for (; j<nc; ++j) ac[j]=mxGetCell(q.a,j);

            contractQS_itags(nc,ac,level+1,vflag,X[l],P_);
         }
      }
      else if (mark[i]&3) { 
         Q[l].set(argin[i]);
      }
      else wblog(FL,"ERR %s() mark[%d]=%d",FCT,i,mark[i]);
   }

   if (len<2) {
      if (len==1 && !level) { 
         X[0].save2(C);
         if (Q[0].conj) { C.Conj(); }
         if (Q[0].sidx) { wbperm P_;
            int i=P_.initStr(0,0,Q[0].sidx.data,1);
            if (P.len) { wblog(FL,
               i>0 ? "ERR got multiple permutations '%s' %s)"
                   : "ERR invalid usage (perm '%s' / %s)",
               STR(P),Q[0].sidx.data);
            }
            if (i<0) wblog(FL,
               "ERR invalid permutation '%s'",Q[0].sidx.data);
            P_.save2(P);
         }

         if (P) { unsigned r=C.rank(FL);
            if (P.len<=r) {
               if (vflag) { wblog(FL,
                  "--> applying permutation %s",STR(P)); }
               C.Permute(P);
            }
            else { wblog(FL,
               "ERR invalid input permutation '%s' (len=%d/%d)",
               STR(P),P.len,r);
            }
         }
      }
      else { wblog(FL,
         "ERR %s() invalid usage #2: at least\n"
         "two QSpaces required for each contraction (%d)",FCT,len);
      }
      return len;
   }

   unsigned k=len-1, ra,rb;
   ctrIdx ica,icb;

   QSpace<gTQ,TD> &B=X[k]; 

   if (k) match_regex_itag(Q[k].otags,X[k-1].itags,"opB","A");
   Q[k].apply(B,vflag); 
   if (Q[k].conj) { icb.conj=1; }

   for (--k; k<len; --k) {
      QSpace<gTQ,TD> &A=X[k];

      match_regex_itag(Q[k].otags,B.itags,"opA","B");
      Q[k].apply(A,vflag); 
      ica.conj=(Q[k].conj ? 1 : 0); 

      if (A.isEmpty() || B.isEmpty()) {
         if (vflag>1) {
            unsigned n=64; char s[n];
            char i[2]={ A.isEmpty(), B.isEmpty() };
            snprintf(s,n,"cell contraction (%d): %s empty",level,
              i[0] ? (i[1] ? "both QSpaces are": "1st QSpace is")
                   : (i[1] ? "2nd QSpace is" : "NEITHER (!?) QSpace is"));
            vflag += 256; 
            wblog(FL,"%2d) %s",vflag>>8,s);
         }
         C.init(); return len;
      }

      if (vflag>1) { unsigned l,n=16;
         char s[]="---------------------- ictr rank";
         if (level)
              { l=snprintf(s,n,"level %d",level); }
         else { l=snprintf(s,n,"base"); }

         if (len!=2 && l<n) {
            l+=snprintf(s+l,n-l," (%d/%d)",len-k-1,len-1);
         }; if (l<n) { s[l]= ' '; } else { s[n-1]=s[n]='?'; }

         vflag += 256; 
         wblog(FL,"%2d) cell contraction @ %s",vflag>>8,s);
      }

      ra=A.rank(FL); rb=B.rank(FL); 
      A.matchITags(FL,B,ica,icb);   

      for (i=0; i<ica.len; ++i) if (ica[i]>=ra) wblog(FL,
         "ERR index out of bounds (%d; %d)",ica[i],ra);
      for (i=0; i<icb.len; ++i) if (icb[i]>=rb) wblog(FL,
         "ERR %s() index out of bounds (%d; %d)",PROG,icb[i],rb);

      if (vflag>1) { wblog(FL,
         " *     %-37s @ %-5s %2d\nX  %-37s @ %-5s %2d",
         IT2STR(A), STR(ica), ra, IT2STR(B), STR(icb), rb);
      }
      if (level==0 && k+1==len) {
         if (!A.qtype.permitsOM(ra+rb-ica.len-icb.len)) {
            wblog(FL,"TST %s() ~~~~~~~~~~~~~~~~~~~~~~~~",FCT);
            C.init(); C.mt=Wb::MEX_RETURN; 
         }
      }

      try {
         A.contract(FL,ica,B,icb, C, P); 
      }
      catch (...) {
         A.print("A"); B.print("B"); 
         wblog(FL,"ERR invalid cell-contraction (l=%d: k=%d/%d)\nhaving "
         "ic = %s <> %s, P=[%s]",level,k,len,STR(ica),STR(icb),STR(P));
      }
      if (vflag>1) {
         wblog(FL, " *  =  %-45s %2d",IT2STR(C),ra+rb-2*ica.len);
      }

      C.UnsetFlags();

      C.SkipZeroData(); 
      if (level<=0) {
         C.NormCGW(); 
      }

      if (k) { C.save2(B); icb.conj=0; }
      else { C.ctime=Wb::getTimeNow(); }
   }

   return len;
};

int is_valid_cell_ctr(
   unsigned na, const mxArray *aa[],
   unsigned level, 
   unsigned &lmax  
){
   int rval=0, e=0, q; 
   unsigned i=0,n,nQS=0;
   icFlags Q; 

   if (lmax<level) { lmax=level; }

   for (; i<na; ++i) { const mxArray *a=aa[i];
      if (mxIsCell(a)) {
         unsigned i=0, nc=mxGetNumberOfElements(a);
         const mxArray* ac[nc]; ++nQS;
         for (; i<nc; ++i) { ac[i]=mxGetCell(a,i); }

         if ((q=is_valid_cell_ctr(nc,ac,level+1,lmax))>0)
              { rval+=q; }  
         else { rval=q; break; }
      }
      else if (mxIsChar(a) || mxIsNumeric(a)) {
         if (Q.set(a,0)<=0) { wbstring s(a); wblog(FL,
            "ERR %s() invalid QSpace cell-structure ('%s' @ l=%d)",
            FCT,s.data,level);
         }
      }
      else if ((n=mxGetNumberOfElements(a))!=1) {
         wblog(FL,"ERR invalid QSpace cell-structure "
         "(%s @ l=%d, %d elements)",mxGetClassName(a), level,n);
      }
      else {
         try { q=mxIsQSpace(FL,a,'c'); }
         catch (...) { ++e; q=0; }

         if (q>0 && !e) { ++nQS; 
            rval+=(q&4 ? 65536 : 256);   
         }
         else {
            wblog(FL,"ERR invalid QSpace cell-structure "
              "(%s @ l=%d; e=%d)",mxGetClassName(a),level,q);
            rval = (q<0 ? q : 0);
         }
      }
   }
   if (nQS>1) { rval+=(nQS-1); }
   else if (!nQS || level) { wblog(FL,
      "ERR invalid QSpace cell-structure\n"
      "(l=%d: only %d QSpace/cell objects encountered)",level,nQS);
   }

   return rval;
};

int contractRC(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[]
){

   int i,l, n=256, n1=256, n2=384, slen=2*n, i0, nc=nargin/4;
   char sbuf[slen];

   if ((nargin-1)%4) wblog(FL,"ERR %s() "
      "invalid usage [nargin: (%d-1)%%4] !?",FCT,nargin);

   for (i=0; i<nargin; ++i) {
      if (!mxIsChar(argin[i])) wblog(FL,"ERR %s() "
      "invalid usage (expecting strings; %d/%d) !?",FCT,i+1,nargin);
   }

   QSet<gTQ> Qa,Qb;
   CRef<gTQ> A,B,C;
   ctrIdx ia,ib;

   QType t(FL,argin[0]); 

   memset(sbuf,0,slen);
   mxGetString(argin[0],sbuf,n); l=strlen(sbuf); sbuf[l]=' '; sbuf[++l]=0;

   if (nc!=1)
        wblog(FL," *  %s() got %d contractions%N",myname,nc);
   else wblog(FL," *  %s() got single contraction%N",myname);

   for (i=0; i<nc; ++i) { i0=4*i+1;

      mxGetString(argin[i0  ],sbuf+l,n); Qa.init_str(FL,sbuf); ia.init(FL,argin[i0+1]);
      mxGetString(argin[i0+2],sbuf+l,n); Qb.init_str(FL,sbuf); ib.init(FL,argin[i0+3]);

      if (nc>1) {
         if (nc>9)
              { snprintf(sbuf+n,32,"%2d)", i+1   ); }
         else { snprintf(sbuf+n,32,"%d/%d",i+1,nc); }
      }
      else strcat(sbuf+n,"-->");

      wblog(FL,"%s %s @ %s | %s @ %s",sbuf+n,STR(Qa),STR(ia),STR(Qb),STR(ib));

      const CData<gTQ,RTD> &a = gCS.getBUF(FL,Qa,'l');
      const CData<gTQ,RTD> &b = gCS.getBUF(FL,Qb,'l');

      A.initBase(&a);
      B.initBase(&b);

      gXS.contractCGR(FL,A,ia,B,ib,C);

      snprintf(sbuf+n1,128,"%s @ %s",STR(Qa),STR(ia));
      snprintf(sbuf+n2,128,"%s @ %s",STR(Qb),STR(ib));

      printf("\n   ==> %-30s  %s\n",sbuf+n1,sbuf+n2);
      printf("   cgd %-30s  %s\n",STR(a),STR(b));
      printf("   cgb %-30s  %s\n",STR(A),STR(B));
      printf("   ==> %s\n\n   ==> %s\n\n",STR(C),STR_(C.cgb));
   }

   argout[0]=C.toMx(); if (nargout>1) {
   argout[1]=C.cgb->toMx(); }

   return i;
};

int contractCD(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[]
){

   int i,l,n=256,n1=256,n2=384;
   char sbuf[2*n];

   if (nargin!=5 || nargout>1) wblog(FL,"ERR %s() "
      "invalid usage [nargin: (%d-1)%%4] !?",FCT,nargin);

   for (i=0; i<nargin; ++i) {
      if (!mxIsChar(argin[i])) wblog(FL,"ERR %s() "
      "invalid usage (expecting strings; %d/%d) !?",FCT,i+1,nargin);
   }

   QSet<gTQ> Qa,Qb;
   ctrIdx ia,ib;

   QType t(FL,argin[0]); 

   memset(sbuf,0,2*n);
   mxGetString(argin[0],sbuf,n); l=strlen(sbuf); sbuf[l]=' '; sbuf[++l]=0;

   ++argin;

   mxGetString(argin[0],sbuf+l,n); Qa.init_str(FL,sbuf); ia.init(FL,argin[1]);
   mxGetString(argin[2],sbuf+l,n); Qb.init_str(FL,sbuf); ib.init(FL,argin[3]);

   wblog(FL,"--> %s @ %s | %s @ %s",STR(Qa),STR(ia),STR(Qb),STR(ib));

   const CData<gTQ,RTD> &a = gCS.getBUF(FL,Qa,'l');
   const CData<gTQ,RTD> &b = gCS.getBUF(FL,Qb,'l');
   cdata<RTD> c;

   a.cgd.contract(FL,ia,b.cgd,ib,c,0,0,0);

   snprintf(sbuf+n1,128,"%s @ %s",STR(Qa),STR(ia));
   snprintf(sbuf+n2,128,"%s @ %s",STR(Qb),STR(ib));

   printf("\n   ==> %-30s  %s\n",sbuf+n1,sbuf+n2);
   printf("   cgd %-30s  %s\n",STR(a),STR(b));
   printf("   ==> %s\n\n",SSTR(c));

   argout[0]=c.toMx();

   return 0;
};

