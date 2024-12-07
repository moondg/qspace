
/* Undocumented ====================================================== *

   Usage #2 // Wb,Sep26,16

      low-level CData contraction (adapted from tensorProdRC)
      getIdentityQS --rcd sym R1 R2

   Usage #3 // Wb,Sep30,16

      display current internal status of loaded/written RCStore data
      getIdentityQS --stat

   Also, a second output argument may return the contents of gCS.

 * =================================================================== */

char USAGE[]=""; // outsourced to getIdentityQS.m // Wb,Jan12,19

#define PROG_TAG "gId"

#ifdef MATLAB_MEX_FILE
   #define PROG mexFunctionName()
#else
   #define PROG "getIdentityQS"
#endif

#define LD_CLEBSCH_QS
#include "wblib.h"

int getIdentityCD(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[]
);

void check_itag_markerA(
   unsigned nargin, const mxArray** argin,
   int &l, 
   wbvector< wbstring > &t3s, 
   char id, const wbindex &ic
) {
   if (unsigned(l)>=nargin || !mxIsChar(argin[l])) return;

   unsigned j=2; char s[16]; s[0]=id;

   if (mxGetString(argin[l],s+2,14)) wblog(FL,"ERR %s() "
      "input string (arg#%d) out of bounds)",myname,l+1);
   if (strncmp(s+2,"-A",2)) return;

   if (ic.len!=1) wblog(FL,"ERR itags() "
      "invalid option `%s' (ic%c must have len=1 then)",s+2,tolower(id));
   if (ic[0]!=0 && ic[0]!=1) wblog(FL,"ERR itags() "
      "invalid option `%s' (ic%c=%d must be in [1,2])",s+2,tolower(id),ic[0]+1);
   s[1]=ic[0]+'1';

   if (t3s.len<=j || t3s[j]) wblog(FL,
      "ERR invalid usage (itag %d already set, having %s / %s)",
      j+1, s+2, t3s.len>j ? t3s[j].data : "(null)");

   t3s[j]=s; ++l; 
};

void check_itag_markers(
   const char *F, int L, unsigned nargin, const mxArray** argin,
   int &l, 
   wbvector< wbstring > &t3s, char id
){
   int i=0, j=0; char s[16]; s[0]=id;

   while (unsigned(l)<nargin && mxIsChar(argin[l])) {
      if (mxGetString(argin[l],s+1,15)) {
         if (!mxIsChar(argin[l])) wblog(F_L,
            "ERR %s() invalid input option of type %s (arg#%d)",
            myname, mxGetClassName(argin[l]), l+1);
         else { wblog(F_L,
            "ERR input itag string out of bounds (arg#%d; 16)",l+1);
         }
      }
      else if (!s[1] || !s[2] || !s[3]) {
          if (F) wblog(FL,
             "ERR invalid itags() option '%s' (arg#%d)",s+1,l+1);
          break;
      }
      i=s[2]-'1'; j=s[4]-'1'; 

      if (s[1]!='#' || s[3]!=':' || isdigit(s[5]) ||
         (i<0 || i>8 || j<0 || j>2)
      ){ if (F || s[1]=='#') wblog(FL,"ERR itags() "
         "invalid arg#%d `%s'\n(expected format '#n:m..')",l+1,s+1);
         break;
      }

      if (unsigned(j)>=t3s.len) wblog(F_L,
         "ERR itags() index out of bounds (%d/%d)",i,t3s.len);
      if (t3s[j]) wblog(F_L,"ERR itags() invalid usage "
         "(#n:%d.. already set)\nhaving %s / %s",j+1,s+1,t3s[j].data);
      t3s[j]=s; ++l;
   }
};

template <class TQ, class TA, class TB, class TC>
void set_ctags(
   const char *F, int L, QSpace<TQ,TC> &C,
   const wbvector< QSpace<TQ,TA> > &A,
   const wbvector< QSpace<TQ,TB> > &B,
   wbvector< wbstring > &t3s
){
   unsigned i=-1,j=-1,k=0,l;
   char c, *s; itag_ t;

   if (t3s.len!=3) wblog(F_L,"ERR %s() got len=%d !?",myname,t3s.len);

   for (; k<3; ++k) { if (t3s[k]) {
      s=t3s[k].data; c=s[0]; l=strlen(s);

      if (!l) wblog(FL,"ERR %s() got `%s' !?",myname,s);
      if (c=='A' && A.len!=1) wblog(F_L,"ERR %s() "
         "#n.$m.. options only with single QSpace (A.len=%d)",myname,A.len);
      if (c=='B' && B.len!=1) wblog(F_L,"ERR %s() "
         "#n.$m.. options only with single QSpace (B.len=%d)",myname,B.len);

      if (c=='C') { 
         if (s[1]=='-') {  
            if (s[2]=='m' && (s[3]==':' || s[3]=='!')) { 
               unsigned n=strlen(s+1); 
               char cflag=0, mflag=0;  

               for (  i=n; i<=n && s[i]==CC_ITAG; --i) { ++cflag; s[i]=0; }
               for (; i<=n && s[i]==QS_MARK_DUAL; --i) { ++mflag; s[i]=0; }

               if (((cflag % 2)==0) == (s[3]==':')) { ++mflag; } 

               t.init(F_L,s+4);  if (mflag % 2) { t.MarkDual(); }
            }
            else { wblog(FL,
              "ERR %s() invalid itag option '%s' (arg#%d)",myname,s+1,l+1);
            }
         }
         else {
            if (s[1] && !isalnum(s[1])) wblog(FL,
               "ERR %s() invalid itag specification '%s'",myname,s+1);
            t.init(F_L,s+1);
         }
         if (C.itags.len==3) { j=2; } 
      }
      else if (l>2 && !strncmp(s+2,"-A",2)) {
         itag_ tloc; i=s[1]-'1'; j=2;
         if (c=='A') {
            if (!A[0].isAtensor() || i>1) wblog(FL,
               "ERR itags() invalid A-tensor (%c: %s; %d)",c,s+2,i+1);
            if (i==0)        
                 { t=A[0].itags[1]; }
            else { t=A[0].itags[0]; }; tloc=A[0].itags[2];
         }
         else if (c=='B') {
            if (!B[0].isAtensor() || i>1) wblog(FL,
               "ERR itags() invalid A-tensor (%c: %s; %d)",c,s+2,i+1);
            if (i==0)        
                 { t=B[0].itags[1]; }
            else { t=B[0].itags[0]; }; tloc=B[0].itags[2];
         }
         else wblog(F_L,"ERR %s() invalid c=%c<%d>",myname,c,c);

         C.SetTag(2,tloc); 

         if (s[4]) {
            if (s[4]!=':') wblog(FL,"ERR %s() invalid option `%s'",myname,s+2);
            if (s[5]) { char sx[16]; 
               snprintf(sx,16,"%s%s",STR(t.deConj()),s+5);
               t.init(FL,sx);
            }
         }
      }
      else if (s[1]=='#') {
         if ((c!='A' && c!='B') || !s[2] || !s[3]) wblog(F_L,
            "ERR %s() invalid c=%c<%d> (%s)",myname,c,c,s);
         i=s[2]-'1'; j=s[4]-'1';

         if ((c=='A' && i>=A[0].itags.len) ||
             (c=='B' && i>=B[0].itags.len)) wblog(FL,
            "ERR itags() out of bounds (%c: %s @ %d)",c,s+1,i+1);

         if (c=='A')
              { t=A[0].itags[i]; }
         else { t=B[0].itags[i]; }

         if (s[5]) { char sx[16];
            snprintf(sx,16,"%s%s",STR(t.deConj()),s+5); 
            t.init(FL,sx);
         }
      }
      else wblog(FL,"ERR itags() invalid usage `%s' (%d)",s,k);

      if (C.itags.len==3) { 
         if (j>=C.itags.len) wblog(F_L,
            "ERR itags() out of bounds (%c: %s @ %d)",c,s+1,i+1);
         if (j!=k) wblog(F_L,
            "ERR %s() invalid usage (%s / k=%d)",myname,s+1,k);
         C.SetTag(j+1,t); 
      }
      else if (c=='C')
           { C.SetTags(t); } 
      else { wblog(FL,"WRN %s() unexpected usage (%s)",myname,s); }
   }}
};

template <class TQ, class TD>
void check_idx_bounds(const char *F, int L,
   const wbindex &I, const wbvector< QSpace<TQ,TD> > &A, const char *istr) {

   unsigned r=0, i=0, r2;
   for (i=0; i<A.len; ++i) { r2=A[i].rank(F_L); if (!r || r>r2) { r=r2; }}

   if (I.len>r) wblog(F_L,
      "ERR %s() invalid index %s=[%s]",PROG,istr,STR(I+1));
   if (!istr) wblog(F_L,"ERR %s() got null string !?",FCT);

   for (i=0; i<I.len; ++i) {
      if (I.data[i]>=r) wblog(F_L,
      "ERR %s() index out of bounds (%s=[%s])",PROG,istr,STR(I+1));
   }
};

void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){ Wb::CleanUp aclu; try { 

   int l=0, iB=0;
   char vflag=1, vflag_=0, zflag=0;
   char isra, isrb=-1; 

   wbindex Ia,Ib; wbperm P;
   mxArray *a=NULL;

   wbvector< wbstring > t3s(3);

   MX_CHECK_HELPER_NARGS(1,-1,2); 

   if (nargin==1 && mxIsChar(argin[0])) {
      if (Mx::IsEqual(argin[0],"--stat")) { 
         CG::MemStat(FL);
      }
      else if (Mx::IsEqual(argin[0],"--getCG")) {
         argout[0]=gCS.toMx();
         if (nargout>1) { argout[1]=gRS.toMx(); }
      }
      else if (Mx::IsEqual(argin[0],"--rcd")) {
         getIdentityCD(nargout,argout,nargin-1,argin+1);
      }
      else wblog(FL,"ERR %s() invalid usage",myname);
      return;
   }

   isra=mxIsQSpaceVec(0,0,argin[l++]);
   if (!isra && !mxIsQSpaceVec(0,1,argin[l-1],-1,0,0,'C')) wblog(FL,
      "ERR invalid input operator A (arg #%d)%N%N%s",l+1,str);

   if (l<nargin && Mx::IsIndex(argin[l])) {
      Ia.init(argin[l++],1); 
      check_itag_markerA(nargin,argin,l,t3s,'A',Ia);
   }

   check_itag_markers(0,0,nargin,argin,l,t3s,'A');

   if (l<nargin && mxIsStruct(argin[l])) { 
      isrb=mxIsQSpaceVec(0,0,argin[l]);
      if (!isrb && !mxIsQSpaceVec(0,1,argin[l],-1,0,0,'C')) wblog(FL,
         "ERR invalid input operator B (arg #%d)%N%N%s",l+1,str);
      iB=l++;

      if (l<nargin && Mx::IsIndex(argin[l])
          && mxGetNumberOfElements(argin[l])<3 
      ){ Ib.init(argin[l++],1);
         check_itag_markerA(nargin,argin,l,t3s,'B',Ib);
      }
   }

   check_itag_markers(0,0,nargin,argin,l,t3s,'B');

   for (; l<nargin; ++l) {
      if (mxIsChar(argin[l])) { char s[32];
         if (mxGetString(argin[l],s,31)) { s[31]=0;
            wblog(FL,"ERR %s() invalid usage (%d: %s)",myname,l+1,s);
         }

         if (!strcmp(s,"-v")) { ++vflag; } else 
         if (!strcmp(s,"-q")) { --vflag; } else
         if (!zflag) { 
            if (!strncmp(s,"-0",2)) { zflag=1; if (s[2]) {
               if (!strcmp(s,"-0fm")) 
                    { zflag|=2; }
               else { wblog(FL,"WRN invalid option '%s' (using -0)",s); }
            }} else
            if (!strcmp(s,"-z")) { zflag=2; } 
            if (zflag) { continue; }
         }
         if (!t3s[2]) {
            char sx[16]; sx[0]='C'; sx[1]=0;

            if (mxGetString(argin[l],sx+1,15)) { sx[15]=0; wblog(FL,
               "ERR invalid itag option '%s' (arg#%d)",sx,l+1); }
            t3s[2]=sx;
         }
         else if (!P.len && isValidPerm(s)>0) {
            P.initStr(FL,s);
         }
         else wblog(FL,"ERR %s() invalid option '%s'\n"
            "(itag already specified as '%s')",myname,s,t3s[2].data+1
         );
      }
      else if (mxIsNumeric(argin[l])) {
         if (!P.len) { P.init(FL,argin[l],1); }
         else wblog(FL,
            "ERR %s() invalid usage (multiple permutations !?)",myname);
      }
      else wblog(FL,
      "ERR %s() invalid usage (unreocgnized trailing options)",myname);
   }

   vflag_ = (vflag>0 ? vflag-1 : 0);

   if (vflag< 0) { vflag= 0;  } else
   if (vflag==2) { vflag='V'; } else
   if (vflag> 2) { wblog(FL,
      "WRN %s() got multiple options '-v' (%d)",myname,vflag); }

   if (isra) {
      wbvector< QSpace<gTQ,double> > A;
      mxInitQSpaceVec(FL,argin[0],A,'r');
      if (Ia.len) check_idx_bounds(FL,Ia,A,"i1");

      if (zflag && isrb>=0) wblog(FL,   
         "WRN ignoring input option %s",zflag&1 ? "-0":"-z");

      if (isrb>0) {
         wbvector< QSpace<gTQ,double> > B;
         QSpace<gTQ,double> C;

         mxInitQSpaceVec(FL,argin[iB],B,'r');
         if (Ib.len) check_idx_bounds(FL,Ib,B,"i2");

         C.initIdentityCG(A,Ia,B,Ib,vflag); 
            if (t3s.any()) set_ctags(FL,C,A,B,t3s);
            if (P.len) C.Permute(P);
         a=C.save2Mx(vflag_);
      }
      else if (!isrb) {
         wbvector< QSpace<gTQ,wbcomplex> > B;
         QSpace<gTQ,wbcomplex> C;

         mxInitQSpaceVec(FL,argin[iB],B);
         if (Ib.len) check_idx_bounds(FL,Ib,B,"i2");

         C.initIdentityCG(A,Ia,B,Ib,vflag); 
            if (t3s.any()) set_ctags(FL,C,A,B,t3s);
            if (P.len) C.Permute(P);
         a=C.save2Mx(vflag_);
      }
      else {
         QSpace<gTQ,double> C;
         if (!Ib.isEmpty()) wblog(FL,"ERR invalid usage (%d)",Ib.len);

         if (A.len==1 && A[0].isScalar()) { 
            if (t3s.any()) wblog(FL,
               "ERR %s() got itag specs with scalar QSpace !?",myname);
            if (P.len) wblog(FL,
               "ERR %s() got P=[%s] with scalar QSpace !?",myname,STR(P));
            C=A[0]; C.DATA[0]->data[0]=1;
         }
         else {
            C.initIdentityCG(A,Ia,zflag); 
            if (t3s.any()) set_ctags(FL,
               C, A, wbvector<QSpace<gTQ,double> >(), t3s);
            if (P.len) C.Permute(P);
         }
         a=C.save2Mx(vflag_);
      }
   }
   else {
      wbvector< QSpace<gTQ,wbcomplex> > A;
      mxInitQSpaceVec(FL,argin[0],A,'r');
      if (Ia.len) check_idx_bounds(FL,Ia,A,"i1");

      if (zflag && isrb>=0) wblog(FL,   
         "WRN ignoring input option %s",zflag&1 ? "-0":"-z");

      if (isrb>0) {
         wbvector< QSpace<gTQ,double> > B;
         QSpace<gTQ,double> C;

         mxInitQSpaceVec(FL,argin[iB],B,'r');
         if (Ib.len) check_idx_bounds(FL,Ib,B,"i2");

         C.initIdentityCG(A,Ia,B,Ib,vflag); 
            if (t3s.any()) set_ctags(FL,C,A,B,t3s);
            if (P.len) C.Permute(P);
         a=C.save2Mx(vflag_);
      }
      else if (!isrb) {
         wbvector< QSpace<gTQ,wbcomplex> > B;
         QSpace<gTQ,wbcomplex> C;

         mxInitQSpaceVec(FL,argin[iB],B);
         if (Ib.len) check_idx_bounds(FL,Ib,B,"i2");

         C.initIdentityCG(A,Ia,B,Ib,vflag); 
            if (t3s.any()) set_ctags(FL,C,A,B,t3s);
            if (P.len) C.Permute(P);
         a=C.save2Mx(vflag_);
      }
      else {
         QSpace<gTQ,wbcomplex> C;
         if (!Ib.isEmpty()) wblog(FL,"ERR invalid usage (%d)",Ib.len);

         C.initIdentityCG(A,Ia,zflag); 
            if (t3s.any()) set_ctags(FL,
               C, A, wbvector<QSpace<gTQ,double> > (), t3s);
            if (P.len) C.Permute(P);
         a=C.save2Mx(vflag_);
      }
   }

   argout[0]=a; 
   if (nargout>1) { argout[1]=gCS.toMx(); }

}  catch (Wb::LogException &e) { ExitMsg(e.istr); }
   catch (...) { ExitMsg("caught exception in getIdentityQS"); }
   aclu.Check();
};

int getIdentityCD(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[]
){
   unsigned i;

   if (nargin!=3) wblog(FL,"ERR %s() "
      "invalid usage [nargin=%d] !?",myname,nargin);

   for (i=0; int(i)<nargin; ++i) {
      if (!mxIsChar(argin[i])) wblog(FL,"ERR %s() "
      "invalid usage (expecting strings; %d/%d) !?",myname,i+1,nargin);
   }

   QType t(FL,argin[0]); 

   qset<gTQ> q1(FL,argin[1]),q2(FL,argin[2]);

   wblog(FL,"==> %s() %s: %s x %s",myname,STR(t),STR2(q1,t),STR2(q2,t));

   gCS.add_CData(t,q1,q2);

   return 0;
};

