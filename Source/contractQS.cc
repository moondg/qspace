
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

#define LOAD_CGC_QSPACE
#include "wblib.h"

template<class TA, class TB, class TC>
int CONTRACT_QS(
    const QSpace<gTQ,TA> &A, const ctrIdx &ica,
    const QSpace<gTQ,TB> &B, const ctrIdx &icb,
    QSpace<gTQ,TC> &C, wbperm P,
    char cg_preview=0 
);

char isComplexQS(unsigned na, const mxArray *aa[], unsigned level=0);

class icFlags { 
 public:

   icFlags(const mxArray* ain=NULL) : a(ain), conj(0) {
      checkQSpace(a); };

  ~icFlags() {}; 

   icFlags& checkQSpace(const mxArray *a);

   char check_arg(const mxArray *ain) { a=ain;
      if (a) {
         if (mxIsCell(a)) { return 3; }
         if (mxIsChar(a)) { return (set(a,0) ? 0:1); }
         if (mxIsQSpace(0,0,a,'c')>0) return 2;
      }
      return 0;
   };

   icFlags& init(const mxArray *ain) {
      if (!ain) wblog(FL,"ERR %s() got NULL QSpace or cell !?",FCT);
      a=ain; checkQSpace(a);
      ktags.init(); otags.init(); conj=0;
      return *this;
   };

   char set(const mxArray *a, const char *istr="");

   template <class TQ, class TD>
   void apply(QSpace<TQ,TD> &X);

   const mxArray *a;  
   wbstring ktags;    
   wbstring otags;    
   char conj;         

 protected:
 private:

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

   wbstring S(FL,a); const char *s=S.data;

   if (!s || !s[0]) {
      if (istr) wblog(FL,"ERR %s() %s\n"
        "invalid string option for cell-contraction (%s)",
         FCT, istr, s? "empty":"null");
      else return 1; 
   }

   if (!strncasecmp(s,"-op:",4)) { otags=s+1; }
   else {
      unsigned l=0;
      if (s[0]=='!') {
         for (++s; s[l]; ++l) { if (!isdigit(s[l])) break; }
      }
      if (s[l]) {
         if (s[l]==CC_ITAG && !s[l+1]) { conj=1; }
         else {
            if (istr) wblog(FL,
               "ERR %s() %s\ninvalid cell-contraction option\n"
               "(expecting '[!##][*]': got %s)",FCT,istr,S.data);
            else return 1; 
         }
      }
      ktags=s; 
      if (conj) ktags.data[l]=0; 
   }

   return 0;
};

template <class TQ, class TD>
void icFlags::apply(QSpace<TQ,TD> &X) {

   if (X.isEmpty()) { return; } 
   unsigned r=-1;

   if (otags.len) { 
      if (!strncmp(otags.data,"op:",3)) {  
         int q=X.isOperator(&r);
         if (q<=0) { 
            wblog(FL,"ERR failed to set itags for rank%+d QSpace\n"
              "'%s' %s ('%s', q=%d/%d, e=%d)",
               r, X.itags2Str().data, X.otype2Str().data, otags.data,
               X.itags.isOp(), X.itags.got_op_labels(), q
            );
         }
      }
      X.init_itags(FL,otags.data); 
   }
   else { r=X.rank(FL); }

   if (ktags.len) {
      unsigned i=0, k;
      if (r!=X.itags.len) wblog(FL,
         "ERR %s() valid set of info.itags required (%s; %d)",
         FCT,IT2STR(X),r
      );

      if (r) { 
      for (; i<ktags.len; ++i) { k=ktags[i]; if (!k) break;
         k-='0'; 
         if (!k || k>r) wblog(FL,"ERR %s() index in cc-string "
            "out of bounds (%d: '%s'; %d/%d)",FCT,i,ktags.data,k,r);
         X.SetFlag(FL,k,1);
      }}
   }
};

mxArray* contract_plain(int nargin, const mxArray *argin[]);

template<class TD>
unsigned contractQS_itags(
   int nargin, const mxArray *argin[],
   unsigned level, unsigned vflag, QSpace<gTQ,TD> &C
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
    MX_CHECK_HELPER_NARGS(2,-1,3); 

    if (nargin==1) {
       if (Mx::IsEqual(argin[0],"--stat")) { 
          CG::MemStat(FL);
          return;
       }
       if (Mx::IsEqual(argin[0],"--getCG")) {
          argout[0]=gCS.toMx(); if (nargout>1) {
          argout[1]=gXS.toMx(); if (nargout>2) {
          argout[2]=gRS.toMx(); }}
          return;
       }
       wblog(FL,"ERR %s() invalid usage",FCT);
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

    if (nargin>1 && isCtrIdx(argin[1])) {
       a=contract_plain(nargin,argin);
    }
    else {
       unsigned i,
          nc=nargin,
          level=0;  
       char vflag=0; wbperm P;
       char q, nQS=0;

       icFlags Q;
       for (i=0; i<(unsigned)nargin; ++i) {
          if (!(q=Q.check_arg(argin[i]))) { nc=i; break; }
          if (q==2  || q==3 ) { ++nQS; }
       }

       if (nQS<2) wblog(FL,
          "ERR %s() invalid/incomplete QSpace/cell pair (%d)",myname,nQS);

       if (nc<(unsigned)nargin) {
          OPTS opts; unsigned l=nc;

          if (int(l)<nargin && !mxIsChar(argin[l])) {
             try { P.init(FL,argin[l],1); } 
             catch (...) { wblog(FL,"ERR %s() expecting permutation\n"
               "(got %s at input arg #%d/%d)",
               FCT,mxGetClassName(argin[l]),l+1,nargin);
             }; ++l;
          }
          opts.init(argin+l,nargin-l);

          vflag = opts.getOpt("-v");

          opts.checkAnyLeft(FL);
       }

       if (isComplexQS(nc,argin)) {
          QSpace<gTQ,wbcomplex> C;
          contractQS_itags(nc,argin,level,vflag,C);
          C.Permute(P); a=C.toMx();
       }
       else {
          QSpace<gTQ,double> C;
          contractQS_itags(nc,argin,level,vflag,C);
          C.Permute(P); a=C.toMx();
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

char isra=1, isrb=1, cg_preview=0;
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
          mexWarnMsgIdAndTxt("Wb:MEX:contractQS","check calling function");
      ica.Conj(); ++i;
   }
   if (opts.getOpt("conjB")) {
      wblog(FL,"WRN %s() got deprecated option 'conjB'",FCT);
          mexWarnMsgIdAndTxt("Wb:MEX:contractQS","check calling function");
      icb.Conj(); ++i;
   }

   if (opts.getOpt("--cg-preview")) cg_preview=1; 

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
         CONTRACT_QS(A,ica,B,icb,C,P,cg_preview); a=C.toMx();
      }
      else {
         QSpace<gTQ,wbcomplex> B(argin[k],'r'), C;
         CONTRACT_QS(A,ica,B,icb,C,P,cg_preview); a=C.toMx();
      }
   }
   else {
      QSpace<gTQ,wbcomplex> A(argin[0],'r');

      if (isrb) {
         QSpace<gTQ,double> B(argin[k],'r'); QSpace<gTQ,wbcomplex> C;
         CONTRACT_QS(A,ica,B,icb,C,P,cg_preview); a=C.toMx();
      }
      else {
         QSpace<gTQ,wbcomplex> B(argin[k],'r'), C;
         CONTRACT_QS(A,ica,B,icb,C,P,cg_preview); a=C.toMx();
      }
   }

   return a;
};

template<class TA, class TB, class TC>
int CONTRACT_QS(
    const QSpace<gTQ,TA> &A, const ctrIdx &ica,
    const QSpace<gTQ,TB> &B, const ctrIdx &icb,
    QSpace<gTQ,TC> &C, wbperm P,
    char cg_preview 
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

    q=A.contract(FL,ica,B,icb,C,P,cg_preview); 
    if (cg_preview) { C.init(); return q; }

    C.SkipZeroData(); 

    C.ctime=Wb::getTimeNow();

    return q;
};

template<class TD>
unsigned contractQS_itags( 
   int nargin, const mxArray *argin[],
   unsigned level, unsigned vflag,
   QSpace<gTQ,TD> &C 
){

   unsigned i, l=-1, len=0; char mark[nargin];
   icFlags q;

   for (i=0; i<(unsigned)nargin; ++i) { mark[i]=q.check_arg(argin[i]);
      if (mark[i]>1) ++len; else
      if (!mark[i] || (mark[i]==1 && !i)) wblog(FL,
         "ERR %s() invalid QSpace cell-structure (%d,%d)\n"
         "hint usage: ({QSpace} [,string opts])", FCT,level,i
      );
   }

   wbvector<icFlags> Q(len);
   wbvector<QSpace<gTQ,TD> > X(len);

   for (i=0; i<(unsigned)nargin; ++i) {
      if (mark[i]>1) {
         icFlags &q=Q[++l]; q.init(argin[i]);
         if (mark[i]==2) { X[l].init(FL,argin[i],'r'); }
         else {
            if (!mxIsCell(q.a)) wblog(FL,"ERR %s() got non-cell !?",FCT);
            unsigned j=0, nc=mxGetNumberOfElements(q.a);
               const mxArray* ac[nc];
               for (; j<nc; ++j) ac[j]=mxGetCell(q.a,j);
            contractQS_itags(nc,ac,level+1,vflag,X[l]);
         }
      }
      else if (mark[i]==1) { Q[l].set(argin[i]); }
      else if (mark[i]) {
         wblog(FL,"ERR %s() mark[%d]=%d !?",FCT,i,mark[i]);
      }
   }

   if (len<2) wblog(FL,"ERR %s() invalid usage #2: at least\n"
      "two QSpaces required for each contraction (%d)",FCT,len);

   unsigned k=len-1, ra,rb;
   ctrIdx ica,icb; wbperm P;

   QSpace<gTQ,TD> &B=X[k];

   if (k) match_regex_itag(Q[k].otags,X[k-1].itags,"opB","A");
   Q[k].apply(B); if (Q[k].conj) icb.conj=1;

   for (--k; k<len; --k) {
      QSpace<gTQ,TD> &A=X[k];

      match_regex_itag(Q[k].otags,B.itags,"opA","B");
      Q[k].apply(A); if (Q[k].conj) ica.conj=1;

      if (A.isEmpty() || B.isEmpty()) {
         if (vflag) {
            char i[2]={ A.isEmpty(), B.isEmpty() };
            char s[64]; snprintf(s,64,"%.24s(#%d) %s empty", myname, level,
            i[0] ? (i[1] ? "both QSpaces are": "QSpace 1 is")
                 : (i[1] ? "QSpace 2 is" : "NEITHER (!?) QSpace is"));
            wblog(FL,"%s",s);
         }
         C.init(); return nargin;
      }

      ra=A.rank(FL); rb=B.rank(FL); 
      A.matchITags(FL,B,ica,icb);   

      for (i=0; i<ica.len; ++i) if (ica[i]>=ra) wblog(FL,
         "ERR index out of bounds (%d; %d)",ica[i],ra);
      for (i=0; i<icb.len; ++i) if (icb[i]>=rb) wblog(FL,
         "ERR %s() index out of bounds (%d; %d)",PROG,icb[i],rb);

      if (vflag) wblog(FL,
         " *  %s(#%d)\n       %-15s -> %s\n    <> %-15s -> %s",
         myname, level, IT2STR(A), STR(ica), IT2STR(B), STR(icb)
      );

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
         wblog(FL,"ERR invalid cell-contraction (l=%d: k=%d/%d)\n"
        "%s <> %s, P=[%s]",level,k,len,STR(ica),STR(icb),STR(P));
      }

      C.UnsetFlags();
      C.SkipZeroData(); 

      if (k) { C.save2(B); icb.conj=0; }
      else { C.ctime=Wb::getTimeNow(); }
   }

   return nargin;
};

char isComplexQS(
   unsigned na, const mxArray *aa[],
   unsigned level 
){
   icFlags Q; 

   for (unsigned i=0; i<na; ++i) { const mxArray *a=aa[i];
      if (mxIsCell(a)) {
         unsigned i=0, nc=mxGetNumberOfElements(a);
            const mxArray* ac[nc];
            for (; i<nc; ++i) ac[i]=mxGetCell(a,i);
         if (isComplexQS(nc,ac,level+1)) { 
            return 1;
         }
      }
      else if (mxIsChar(a)) {
         if (Q.set(a,0)!=0) { wbstring s(a); wblog(FL,
            "ERR %s() invalid QSpace cell-structure ('%s'; l=%d)",
            FCT,s.data,level);
         }
      }
      else {
         try {
            if (mxIsQSpace(FL,a,'c')<=0 || mxGetNumberOfElements(a)!=1)
               throw "Wb: invalid QSpace";
            if (mxIsQSpace(0,0,a)<=0) return 1;
         }
         catch (...) {
            wblog(FL,"ERR %s() invalid QSpace cell-structure (%s; l=%d)",
            FCT,mxGetClassName(a),level);
         }
      }
   }

   return 0;
};

int contractRC(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[]
){

   int i,l,n=256,n1=256,n2=384,i0, nc=nargin/4;
   char sbuf[2*n];

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

   memset(sbuf,0,2*n);
   mxGetString(argin[0],sbuf,n); l=strlen(sbuf); sbuf[l]=' '; sbuf[++l]=0;

   if (nc!=1)
        wblog(FL," *  %s() got %d contractions%N",myname,nc);
   else wblog(FL," *  %s() got single contraction%N",myname);

   for (i=0; i<nc; ++i) { i0=4*i+1;

      mxGetString(argin[i0  ],sbuf+l,n); Qa.init_str(FL,sbuf); ia.init(FL,argin[i0+1]);
      mxGetString(argin[i0+2],sbuf+l,n); Qb.init_str(FL,sbuf); ib.init(FL,argin[i0+3]);

      if (nc>1) {
         if (nc>9) sprintf(sbuf+n,"%2d)",i+1);
         else  sprintf(sbuf+n,"%d/%d",i+1,nc);
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

