/* ---------------------------------------------------------------------
 * Project : QSpace tensor library (v4.0)
 * Class   : clebsch (for abelian and non-abelian symmetries)
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

#ifndef __WB_CLEBSCH_IO_CC__
#define __WB_CLEBSCH_IO_CC__

// Wb,May17,20

int RCStore::get_rcs_path(
   const char *F, int L, char *path, unsigned len,
   const QType &t, const char *sub1, const char *sub2, const char *file,
   char mflag
 ) const {

   unsigned l,
      l1=(sub1 && sub1[0] ? strlen(sub1) : 0),
      l2=(sub2 && sub2[0] ? strlen(sub2) : 0),
      l3=(file && file[0] ? strlen(file) : 0);
   char rflag=(!l1 && !l2 && !l3); 

   char mdir=(mflag & 1); 
   char sync=(mflag & 2); 
   char last=(mflag & 4); 

   unsigned ir=0, nr=root.numel(); 
   int q=0;

   if (nr<2) wblog(F_L,
      "ERR %s() root not yet initialized (%d)",FCT,nr);
   if (mflag>>3) wblog(FL, 
      "ERR %s() invalid mflag=%d",FCT); 
   if (!l1 && l2) wblog(F_L,
      "ERR %s() invalid usage (l=%d/d/%d)",FCT,l1,l2,l3);

   if (sync) { ir=nr-1; } else  { --nr; } 

   if (!l1 || last) { 
      ir=nr-1; 
   }

   if (rflag) { 
      if (t.isUnknown()) {
         l=snprintf(path,len,"%s",root[ir]);
         if (!Wb::fexist(path,'d')) wblog(F_L,
            "ERR %s() RC directory does not exist\n'%s'",FCT,path);
         if (l>=len) wblog(FL,
            "ERR %s() string out of bounds (%d/%d)",FCT,l,len);
         return +30;
      }
   }

   if (!t.validType()) wblog(F_L,
      "ERR %s() got invalid symmetry (%s)",FCT,STR(t));
   if (!path) wblog(F_L,"ERR %s() got NULL buffer",FCT);

   const char *stag = sync ? "RCSync" : "RCStore";
   unsigned i0,i1,i2,i3, lx=128; char sub[lx]; sub[0]=0;

   if (sync && !l3) wblog(FL,
      "ERR %s() got empty file with %s",FCT,stag);

   l=snprintf(sub,lx,"%s",STR2(t,'t'));               i1=i2=i3=l;
   if (l1 && l<lx) { l+=snprintf(sub+l,lx-l,"/%s",sub1); i2=i3=l; }
   if (l2 && l<lx) { l+=snprintf(sub+l,lx-l,"/%s",sub2);    i3=l; }
   if (l3 && l<lx) { l+=snprintf(sub+l,lx-l,"/%s",file); }
   if ((l1 && l>=lx) || (!l1 && l+16>=lx)) wblog(FL,
      "ERR %s() string out of bounds (%d/%d)\n%s",FCT,l,lx,sub);

   {  unsigned k=0; char path_[len];

      wbindex Ir; Ir.Index(ir,nr-1,+1); 

      for (k=0; k<Ir.len; ++k) {
         char *p = (k ? path_ : path);
         l=snprintf(p,len,"%s/%s",root[Ir[k]],sub);
         if (l>=len) wblog(FL,
            "ERR %s() string out of bounds (%d/%d)\n%s",FCT,l,len,p);
         if (Wb::fexist(p)) { 
            if (k) strcpy(path,path_);
            break;
         }
      }

      if (k<Ir.len) {
         if (rflag)
              { path[l-1]=0; return +30; } 
         else { return +33; }
      }
   }

   if (rflag) {
      snprintf(path+l,len-l,"/%s","RStore"); 
   }

   int lflag=0;

   const char *dir=root[nr-1];

   i0=(dir && dir[0] ? strlen(dir) : 0);
   if (!i0) wblog(FL,
      "ERR %s() got empty root directory (%d)",FCT,nr);

   i1+=(i0+1); 
   i2+=(i0+1);
   i3+=(i0+1);

   path[i3]=0; q=Wb::fexist(path,'d'); if (!rflag) {
   path[i3]='/'; } 

   if (q) { return (l3 ? -30 : +30); }
   if (!mdir) { return (q ? -30 : -0); } 

   g_sym[t]++; 

#ifdef QS_USING_OMP
   Wb::ompGuard dLK(RCS_lock); 

   path[i3]=0; q=Wb::fexist(path,'d'); if (!rflag) { path[i3]='/'; }
   if (q) { return (l3 ? -30 : +30); }
#endif

   path[i1]=0; 
   if (!Wb::fexist(path,'d')) { 
      Wb::mkdir(FL,path); ++lflag;
      gStore.rclog(t,0,0,-1,"\n"
        "   NB! creating %s directories for symmetry %s\n"
        "   ... mkdir %s\n",
      stag, STR(t), Wb::repHome(path).data);
   }
   path[i1]='/';

   path[i2]=0; 

   if (!Wb::fexist(path,'d')) {

      char c=path[i1+1], rcx[]="RCX";

      if (!index(rcx,c)) wblog(FL,"ERR %s() "
         "unexpected path [c=%c; '%s']\n\%s",FCT,c,rcx,path);

      for (int i=0; i<3; ++i) { path[i1+1]=rcx[i];
         if (!Wb::fexist(path,'d')) {
            Wb::mkdir(FL,path); 
            if (++lflag==1) PRINTF("\n");
            gStore.rclog(t,0,0,-1,"   ... mkdir %s\n",Wb::repHome(path).data);
         }
      }

      if (lflag) {
         PRINTF("\n"); char sep[80]; memset(sep,'=',79); sep[79]=0;
         gStore.rclog(t,0,0,"\n"
         ">> This is a centralized log-file on Clebsch-Gordan tensors (CGTs)\n"
         ">> within QSpace for the symmetry %s setup on host %s.\n\n%s\n",
         STR(t), Wb::hostname().data, sep);
      }

      path[i1+1]=c; 

      if (!Wb::fexist(path,'d')) wblog(FL, 
         "ERR %s() invalid R+C store directory\n%s",FCT,path);
   }
   path[i2]='/';

   if (rflag) {
      path[i0]=0;  
      return +20;  
   }

   path[i3]=0;
   if (!Wb::fexist(path,'d')) Wb::mkdir(FL,path);
   path[i3]='/';

   return -20; 
};

template <class TQ>
int RCStore::get_file_name( 
   const char *F, int L, wbstring &file,
   const QSet<TQ> &Q, const char *ext, char mflag) const {

   unsigned l, l1=35, l2=219, n=l1+l2+2; 
   char sbuf[n], *s1=sbuf, *s2=s1+l1;
   char is3 = (ext && (!strcmp(ext,"mp3") || !strcmp(ext,"c1j"))) ? 1 : 0;
   int q;

   if (Q.qdir.len<2 || Q.qdir.len>l1) { 
      if (!Q.qdir.len) wblog(F_L,
        "ERR %s() got empty CData\n%s",FCT,STR(Q));
      else if (Q.qs.norm2()) wblog(FL,
        "ERR %s() got invalid rank-%d CData\n%s",FCT,Q.qdir.len,STR(Q));
   }

   l=snprintf(s1,l1,"%s", is3 ? "++-" : STR(Q.qdir));
   if (l>=l1) wblog(FL,
      "ERR %s() string out of bounds (%d/%d)\n%s",FCT,l,l1,s1);

   s2[0]='('; l=1;
   l+=snprintf(s2+l,l2-l,"%s",Q.QStrS(0,';').data);

   if (l<l2 && is3) {
      char *s=strstr(s2,";"); if (s) { s[0]=0; l=(s-s2); }
   }

   if (l<l2) l+=snprintf(s2+l,l2-l,").%s",ext?ext:"mat");
   else wblog(FL,
      "ERR %s() string out of bounds (%d/%d)\n%s",FCT,l,l2,s2);

   file.init(n);
   q=get_rcs_path(F_L,file.data,file.len,Q.t,"CStore",s1,s2,mflag);

   return (q>0 ? 1:0); 
};

template <class TQ, class TD>
int RCStore::get_file_name( 
   const char *F, int L, wbstring &file,
   const genRG_base<TQ,TD> &R, const char *ext, char mflag
 ) const {

   int i=0; 
   unsigned l, l1=64; char s1[l1];

   if (!R.J.len || R.J.len>99) wblog(F_L,"ERR %s() "
      "got empty symmetry labels (J=[%s])",FCT,STR(R.J));

   s1[0]='('; l=1;

   if (CG::gotQAlpha(R.J.data,R.J.len)) { 
      if (CG::qset2cstr(R.J.data,s1+l,R.J.len)) wblog(FL,
         "ERR %s() failed to obtain compact qset string",FCT);
      l+=R.J.len;
   }
   else {
      wbstring js=R.J.toStrf("", R.q.qlen()<2 ? "":" ");
      l+=snprintf(s1+l,l1-l,"%s",js.data);
   }

   if (l<l1) l+=snprintf(s1+l,l1-l,").%s",ext?ext:"mat");
   if (l>=l1) wblog(FL,
      "ERR %s() string out of bounds (%d/%d)\n%s",FCT,l,l1,s1);

   file.init(256);
   i=get_rcs_path(F_L,file.data,file.len,R.q,"RStore",NULL,s1,mflag);

   return (i>0 ? 1:0); 
};

int RCStore::get_file_name( 
   const char *F, int L, wbstring &file,
   const cgc_contract_id<MTI> &idc, const char *ext, char mflag,
   QType *t 
 ) const {

   unsigned l,rc, l1=35, l2=219, n=l1+l2+2;  
   char sbuf[n], *s1=sbuf, *s2=s1+l1;
   int q;

   ctrIdx ica,icb;
   CData<gTQ,RTD> a,b;

   idc.extract(a,ica,b,icb,'l');
   rc=a.qdir.len+b.qdir.len-ica.len-icb.len;

   if (rc>l1) wblog(FL, 
      "ERR %s() got unlikely rc=%d\na: %s\nb: %s",
      FCT, rc, STR(a), STR(b));
   if (!a.qs.len || a.qdir.len>l1 || !b.qs.len || b.qdir.len>l1)
      wblog(F_L,"ERR %s() got invalid symmetry labels\n" 
      "A: %s\nB: %s",FCT,STR(a),STR(b)
   );

   if (t) (*t)=a.t;

   QDir cdir; cdir.init2val(rc,+1);
   l=rc - a.qdir.nconj(ica) - b.qdir.nconj(icb);
   for (; l<rc; ++l) cdir[l]=-1;

   if (rc) 
        { l=snprintf(s1,l1,"%s",STR(cdir)); }
   else { strcpy(s1,"scalar"); l=6; }

   if (l>=l1) wblog(FL,
      "ERR %s() string out of bounds (%d/%d)\n%s",FCT,l,l1,s1);

   l=1; s2[0]='(';
   l+=snprintf(s2+l,l2-l,"%s)_%s%s(%s)_%s.%s", 
      a.QStr().data, STR(ica), ica.conj?"":" ",
      b.QStr().data, STR(icb), 
      ext?ext:"mat");
   if (l>=l2) wblog(FL,
      "ERR %s() string out of bounds (%d/%d)\n%s",FCT,l,l2,s2);

   file.init(n);
   q=get_rcs_path(F_L,file.data,file.len,a.t,"XStore",s1,s2,mflag);

   return (q>0 ? 1:0); 
};

template <class TQ, class TD>
int RCStore::save_CData(const char *F, int L, const CData<TQ,TD> &A) {

   wbstring fs;
   int q=get_file_name(F_L,fs,(const QSet<TQ>&)A,"cgd",RC_SAVE);

#ifndef WB_SKIP_ASSERT
   if (q>0) { int i; CDATA_TQ Cf;
      A.Load_CRef(FL,Cf,fs.data); i=A.cmpOM(Cf);
      if (i<0) wblog(FL,
         "ERR %s() got lower-OM CData (i=%d)\n%s: %s -> %s",
         FCT,i,STR(Cf),STR(A)
      );
   }
#endif

   char *sx=strstr(fs.data,"/CStore"), sm[8]="";
   unsigned dx=(sx ? sx-fs.data+1 : 0);

   int m=A.gotOM(F_L); if (!m) m=A.checkOM(F_L);
   if (m) { snprintf(sm,8," @ %d",m); }

   gStore.rclog(A.t, PFL, CG_VERBOSE>6 && F,
      "[%s] %s() #%05x %s%s", q>0?"W":"w",FCT,A.cstat.ID,fs.data+dx,sm);

   if (!A.cstat.ctime) wblog(FL,"ERR %s() got empty time stamp / ID\n"
      "%s\n%s", FCT, STR(A), A.cstat.toStr('V').data); 
   if (A.cstat==CGD_REF_INIT) wblog(FL,"ERR %s() got REF_INIT CData\n"
      "%s\n%s", FCT, STR(A), A.cstat.toStr('V').data); 

   gstatC.gotwrite(A.memSize());

   CData<TQ,TD> Ar; Ar.RefInit(FL,A);
   Ar.cstat.t=A.cstat.t; 

   mxArray *a=Ar.toMx();
   mxArray *c=A.cgd.toMx();

   Wb::matFile f(FL,fs.data, A.cgd.D.len<(1<<20) ? "w:std":"w");

   matPutVariable(f.mfp,"CRef",a);
   mxDestroyArray(a); 

   matPutVariable(f.mfp,"cdata",c);
   mxDestroyArray(c); 

   return q;
};

template <class TQ, class TD>
int RCStore::load_CData(
   const char *F, int L, const QSet<TQ> &Q, CData<TQ,TD> &A,
   char bflag
){
   int q=0; wbstring fs; mxArray *a;

   CData<TQ,TD> B;

   if (bflag && A.cstat.ID && Q!=(QSet<TQ>&)A) wblog(F_L,
      "WRN %s() got QSet mismatch\n%s\n%s",FCT,STR(Q),STR(A));
   if (!bflag ^ !A.gotuser_BUF()) wblog(FL,
      "WRN %s() got bflag=%c<%d>/%d",FCT,bflag,bflag,A.gotuser_BUF());

   if ((q=get_file_name(F_L,fs,Q,"cgd"))<=0) {
      if (!q && A.cstat.ID==CID_RANK1_Q0) {
         A.initScalar(Q); 
         return q;
      }

      if (bflag && A.cstat.ID) {
         wblog(F_L,"WRN %s() got missing CData file having (e=%d)\n"
         "Qin: %s\nBUF: %s\n%s",FCT,q,STR(Q),STR(A),STR2(A.cstat,'V'));
      }
      if (F) wblog(FL,
         "ERR %s() missing CData file\n%s",FCT,fs.data);
      A.init(); return q;
   }

   char *sx=strstr(fs.data,"/CStore");
   unsigned dx=(sx ? sx-fs.data+1 : 0);

   if (CG_VERBOSE>4 && F) wblog(PF_L,
      "[r] %s() %s %s",FCT,fs.data+dx,cSTR(bflag));

   Wb::matFile f(FL,fs.data,"r");
   if (!f.mfp) wblog(FL,"ERR %s() failed to open file\n%s",FCT,fs.data);

   a=matGetVariable(f.mfp,"CRef");
   if (!a) wblog(FL,"ERR %s() "
      "failed to read CRef from file\n%s",FCT,fs.data);

   q=B.init_mxCRef(FL,a,0);
   mxDestroyArray(a);

   if (B.qs && (B.qs!=Q.qs)) { str[0]=0;  
     #if __APPLE__
      if (B.qs.anyGE(36) || Q.qs.anyGE(36)) { sprintf_str("\n\n"
         "hint: copied RCStore from linux?\n"
         "macOS has issues with case sensitive files"); }
     #endif

      wblog(FL,"ERR %s() QSet mismatch\n%s / %s%s",FCT,STR(Q),STR(B),str);
   }

   if (bflag=='r') {
      B.cstat.t=CGD_REF_INIT; 
   }
   else {
      if (bflag && bflag!='b') wblog(FL,
         "WRN %s() got bflag=%s",FCT,cSTR(bflag));
      gstatC.aux1++;

      a=matGetVariable(f.mfp,"cdata");
      f.init(); 

      if (!a) {
         if (CG_FIXIT) { wblog(FL,"FIX %s() "
            "got corrupted 'cdata' -> clear data\n%s",FCT,fs.data);
            A.init(); return (q=-1);
         }
         wblog(FL,"ERR %s() failed to read 'cdata' from file\n%s",
         FCT,fs.data);
      }

      try { B.cgd.init(FL,a,0,B.t); }
      catch (...) {
          wblog(FL,"ERR %s() failed to load cdata\n%s",FCT,fs.data);
      }
      mxDestroyArray(a);

      B.checkNormSign(FL);
   }

   B.save2(A); 

   gstatC.gotread(A.memSize());
   return q;
};

template <class TQ>
int RCStore::save_Std3(
   const char *F, int L, const QType &t, const qset<TQ> &J12,
   const wbvector<double> *c2eps) const {

   wbstring fs;

   QSet<TQ> Q; { 
      Q.init3(t); if (J12.len) {
      memcpy(Q.qs.data,J12.data,J12.len*sizeof(TQ)); }
   }

   int q=get_file_name(F_L,fs,Q,"mp3",RC_SAVE);

   char *sx=strstr(fs.data,"/RCStore");
   unsigned dx=(sx ? sx-fs.data+1 : 0);

   gStore.rclog(t, PFL, CG_VERBOSE>6 && F,
      "[%s] %s() %s",q>0?"W":"w",FCT,fs.data+dx);

   Wb::matFile f(FL,fs.data,"w:std");

   mxArray *a=map3_CreateStructMatrix(1,1);
   map3_add2MxStruct(a,0,t,J12);

   matPutVariable(f.mfp,"map3",a); 
   mxDestroyArray(a);

   if (c2eps) { a=c2eps->toMx();
      matPutVariable(f.mfp,"c2eps",a);
      mxDestroyArray(a);
   }

   return q;
};

template <class TQ>
int RCStore::load_Std3(const char *F, int L,
   const QType &t, const qset<TQ> &J1, const qset<TQ> &J2,
   unsigned loadRC) const {

   if (!J1.len || J1.len!=t.qlen() || J1.len!=J2.len) wblog(FL,
      "ERR %s() invalid J12=[%s, %s] having %s",
      FCT,STR(J1),STR(J2),STR(t));

   qset<TQ> J12(J1,J2);
   wbstring fs;

   QSet<TQ> Q; { Q.init3(t);
      memcpy(Q.qs.data,       J1.data, J1.len*sizeof(TQ));
      memcpy(Q.qs.data+J1.len,J2.data, J2.len*sizeof(TQ));
   }

#ifdef QS_USING_OMP
   CG::Guard mLK;
#endif

   int n3=gCS.find_map3(t,J12); 
   if (n3<=0) { 
      #ifdef QS_USING_OMP
      mLK.acquire(FL,Q,"mp3");
      #endif
      n3=gCS.find_map3(t,J12); 
   }
   if (n3>0 && !loadRC) { return gCS.map3[t][J12].size(); }

   int q=get_file_name(F_L,fs,Q,"mp3");
   if (q<=0) {
      if (!F) return -1;

      genRG_struct<TQ,RTD> &B=gRS.Buf(t);

      q=B.getTensorProdReps_gen(J1,J2, F ? TP3_LDM : TP3_TST);
      if (q>0) {
         return gCS.map3[t][J12].size();
      }

      q=get_file_name(F_L,fs,Q,"mp3");
   }

   if (q<=0) {
      if (F) wblog(FL,"ERR %s() missing mp3 file\n%s",FCT,fs.data);
      return 0;
   }

   Wb::matFile f(FL,fs.data,"r");
   if (!f.mfp) wblog(FL,"ERR %s() failed to open file\n%s",FCT,fs.data);

   mxArray *a, *S=matGetVariable(f.mfp,"map3");
   if (!S) wblog(FL,"ERR %s() "
      "failed to read variable 'map3' from file\n%s",FCT,fs.data);

   unsigned n=mxGetNumberOfElements(S); int m;
   int i3=(mxIsStruct(S) ? mxGetFieldNumber(S,"cgr") : -1);

   char *sx=strstr(fs.data,"/CStore");
   unsigned dx=(sx ? sx-fs.data+1 : 0);

   if (n!=1 || i3<0) wblog(FL,
      "ERR %s() invalid mp3 data in RCStore (%d,%d)\n%s",
      FCT, n,i3, fs.data+dx);

   a=mxGetFieldByNumber(S,0,i3);
   m=gCS.add3(t,a, loadRC & MP3_LDC);

   Q.qs.len=2*J1.len; Q.qdir.len=2; 

   if (m>0) { if ((CG_VERBOSE>7 && L) || CG_VERBOSE>8) { wblog(PFL,
         "(+) map3[%03d] %s(%d)",gCS.map3[t].size(), STR(Q), m
      );
   }}
   else {
      if (m==0) { if (CG_VERBOSE>4)
           wblog(PFL,"TST %s() already got %s",FCT,STR(Q)); }
      else wblog(PFL,"ERR %s() m=%d",FCT,m);
   }

   if (loadRC & MP3_LDR) {
      const auto &M=gCS.map3[t][J12];
      for (auto it=M.begin(); it!=M.end(); ++it) {
         const qset<TQ> &J=it->first;
         const genRG_base<TQ,RTD> &R=gRS.Buf(t).RSet[J];

         if (R.isEmpty()) {
            int q=gStore.load_RSet(0,1,t,J);
            if (q<=0) wblog(FL,
               "ERR %s() %s irep (%s) not yet in RCStore (e=%d)",
               FCT, STR(t), STR(J), q
            );
         }
      }
   }

   mxDestroyArray(S);
   return m;
};

template <class TQ, class TD>
int RCStore::save_RSet(
   const char *F, int L, const genRG_base<TQ,TD> &R, char xflag) {

   wbstring fs;
   int q=get_file_name(F_L,fs,R,"rep",RC_SAVE);

   if (q>0 && xflag) {
      genRG_base<TQ,TD> X=R;
      gStore.load_RSet(0,0,X.q,X.J); 

      double e=X.normDiff(FL,R);
      if (e>1E-12) wblog(FL,"ERR %s() "
         "got RSet[%s] inconsistency (%.3g)",FCT,STR(X.J), e);
      return 2;
   }

   char *sx=strstr(fs.data,"/RCStore");
   unsigned dx=(sx ? sx-fs.data+1 : 0);

   gStore.rclog(R.q, PFL, CG_VERBOSE>6 && F,
      "[%s] %s() %s",q>0?"W":"w",FCT,fs.data+dx);

   mxArray *a=R.toMx();
   Wb::matFile f(FL,fs.data, R.Z.dim1<(1<<16)?"w:std":"w");

   matPutVariable(f.mfp,"RSet",a); mxDestroyArray(a);

   a=numtoMx(double(R.Z.dim1));
   matPutVariable(f.mfp,"qdim",a); mxDestroyArray(a);

   a=R.istr.toMx(); 
   matPutVariable(f.mfp,"istr",a); mxDestroyArray(a);
   a=numtoMx(R.err);
   matPutVariable(f.mfp,"err", a); mxDestroyArray(a);

   gstatR.gotwrite(R.memSize());
   return q;
};

template <class TQ>
int RCStore::load_RSet(
   const char *F, int L, const QType &t, const qset<TQ> &J, char rclog
){
   wbstring fs;
   genRG_base<TQ,RTD> R;

   QSet<TQ> Q; Q.init1(t,J.data);
#ifdef QS_USING_OMP
   CG::Guard qLK(FL,Q,"gRS"); 
#endif

   R.q=t; R.J=J;

   if (R.J.len!=t.qrank()) wblog(FL,"ERR %s() "
      "got invalid qlabels (%s) len=%d/%d",FCT,STR2(J,t),J.len,t.qrank());

   if (get_file_name(F_L,fs,R,"rep")<=0) {
      if (F) { wblog(FL,
         "ERR %s() got missing irrep (%s)\n%s\nhint: %s",
         FCT,STR(J),Wb::repHome(fs.data).data, strstr(myname,"compactQS") ?
         "try to generate using getIdentityQS" : "missing file");
      }
      return 0;
   }

   if (CG_VERBOSE>8 && L) { rclog |= RCL_V1; }

   if (rclog) { char sbuf[256];
      snprintf(sbuf,256,
         "(+) RBUF[%03d] %s",int(gRS.Buf(t).RSet.size()),STR(Q));
      if (rclog & RCL_V2)
           gStore.rclog(t, PF_L, rclog & RCL_V1, "%s",sbuf);
      else wblog(FL,"%s",sbuf);
   }

   Wb::matFile f(FL,fs.data,"r");
   if (!f.mfp) wblog(FL,"ERR %s() failed to open file\n%s",FCT,fs.data);

   unsigned n;
   mxArray *a=matGetVariable(f.mfp,"RSet");

   if (!a) wblog(FL,"ERR %s() "
      "failed to read variable 'RSet' from file\n%s",FCT,fs.data);
   if ((n=mxGetNumberOfElements(a))!=1) wblog(FL,
      "ERR %s() got invalid RSet (len=%d)",FCT,n);

   gRS.add(t,a,0); 
   mxDestroyArray(a);

   genRG_base<TQ,RTD> *Rp=gRS.find_RSet(t,J.data);
   if (!Rp) {
      char *f=rindex(fs.data,'/'); if (f) { ++f; } else { f=fs.data; }
      wblog(FL,"ERR %s() failed to load RSet %s %s\n"
      "%s [ %s ] r=%d/%d", FCT,STR(t),J.toStr(t).data, f,
      J.wbvector<TQ>::toStr().data, J.len, t.qrank());
   }

   a=matGetVariable(f.mfp,"istr");
   if (a) { Rp->istr.init(a); mxDestroyArray(a); }

   a=matGetVariable(f.mfp,"err");
   if (a) { mxGetNumber(a,Rp->err); mxDestroyArray(a); }

   return 1;
};

int RCStore::load_XMap(
   const char *F, int L, const cgc_contract_id<MTI> &idc 
 ) const {

   wbstring fs;
   if (get_file_name(F_L,fs,idc,"x3d")<=0) { if (F) wblog(FL,
      "ERR %s() got non-existing file\n%s",FCT,fs.data);
      return 0;
   }

   Wb::matFile f(FL,fs.data,"r");
   if (!f.mfp) wblog(FL,"ERR %s() failed to open file\n%s",FCT,fs.data);

   if (CG_VERBOSE>8) wblog(FL,"<-- %s",Wb::repHome(fs.data).data);

   unsigned n;
   mxArray *a=matGetVariable(f.mfp,"x3m");

   if (!a) {
   #ifdef QS_USING_OMP
      wbstring istr=CG::Guard::Status(idc);
      wblog(FL,"ERR %s() having %s\n"
      "failed to read `x3m' from file%N%N%s%N",FCT,istr.data,fs.data);
   #else
      wblog(FL,"ERR %s() failed to read `x3m' from file%N%N%s%N",
      FCT,fs.data); 
   #endif
   }
   if ((n=mxGetNumberOfElements(a))!=1) wblog(FL,
      "ERR %s() got invalid RSet (len=%d)",FCT,n);

   gXS.add(F_L,idc,a,0);
   mxDestroyArray(a);

   return 1;
};

int RCStore::save_XMap(
   const char *F, int L, const cgc_contract_id<MTI> &idc 
 ) const {

#ifdef QS_USING_OMP
   CG::Guard xLK(FL,idc);
#endif
   const x3map<gTQ,double> &M=gXS.getXBUF(FL,idc); 

   wbstring fs; QType t;
   int q=get_file_name(F_L,fs,idc,"x3d",RC_SAVE,&t); {
      unsigned l=0, i=0;
      for (; i<fs.len && fs[i]; ++i) { if (fs[i]=='/') l=i+1; }

      gStore.rclog(t, PF_L, CG_VERBOSE>6 && F,
      "[%s] %s", q>0?"W":"w", fs.data+l);
   }

   mxArray *a=M.toMx(&idc);

   Wb::matFile f(FL,fs.data,"w:std");

   matPutVariable(f.mfp,"x3m",a); 
   mxDestroyArray(a);

   gstatX.gotwrite(M.memSize());
   return q;
};

inline mxArray* QType::mxCreateStruct(unsigned m, unsigned n) const {
   return mxCreateCellMatrix(m,n);
};

void QType::add2MxStruct(mxArray *S, unsigned i) const {
   if (!S) wblog(FL,"ERR");
   mxSetCell(S,i,toStr('t').toMx());
};

template <class TQ>
mxArray* QMap<TQ>::mxCreateStruct(unsigned m, unsigned n) const {

   const char *fields[]={
      "qvec","Q1","Q2","Q", "D","I1","I2","II","CG"
   };
   return mxCreateStructMatrix(m,n,9,fields);
}

template <class TQ>
void QMap<TQ>::add2MxStruct(mxArray *S, unsigned i, char tst) const {

   if (tst) {
      unsigned s=0; int k=0; 
      if (S==NULL || (s=mxGetNumberOfElements(S))<1 || i>=s
       || (k=mxGetFieldNumber(S,"Q2"))<0) wblog(FL,
      "ERR %s must follow mxCreateStruct()\n%lx, %d/%d, %d",
       FCT,S,i+1,s,k);
   }

   mxSetFieldByNumber(S,i, 0, qvec.toStr().toMx());
   mxSetFieldByNumber(S,i, 1, Q1   .toMx());
   mxSetFieldByNumber(S,i, 2, Q2   .toMx());
   mxSetFieldByNumber(S,i, 3, Q    .toMx());
   mxSetFieldByNumber(S,i, 4, D    .toMx());
   mxSetFieldByNumber(S,i, 5,(I1+1).toMx());
   mxSetFieldByNumber(S,i, 6,(I2+1).toMx());
   mxSetFieldByNumber(S,i, 7,(II+1).toMx());
   mxSetFieldByNumber(S,i, 8, cg3  .toMxP_S()); 
};

template <class TQ>
QMap<TQ>& QMap<TQ>::init(const char *F, int L, const mxArray *S) {

   mxArray *a; if (!S || !mxIsStruct(S)) { init(); return *this; }

   a=mxGetField(S,0,"qvec"); if (a) qvec.init(F,L,a);
   a=mxGetField(S,0,"Q1");   if (a) Q1.init(F,L,a);
   a=mxGetField(S,0,"Q2");   if (a) Q2.init(F,L,a);
   a=mxGetField(S,0,"Q" );   if (a) Q .init(F,L,a);
   a=mxGetField(S,0,"D" );   if (a) D .init(F,L,a);
   a=mxGetField(S,0,"I1");   if (a) I1.init(F,L,a);
   a=mxGetField(S,0,"I2");   if (a) I2.init(F,L,a);
   a=mxGetField(S,0,"II");   if (a) II.init(F,L,a);

   return *this;
}

template <class T>
mxArray* pairPatternCG<T>::mxCreateStruct(unsigned m, unsigned n) const {
   const char *fields[]={"i","j","k","fac"};
   return mxCreateStructMatrix(m,n,4,fields);
};

template <class T>
void pairPatternCG<T>::add2MxStruct(mxArray *S, unsigned l) const {

   mxSetFieldByNumber(S,l, 0, numtoMx(i+1));
   mxSetFieldByNumber(S,l, 1, numtoMx(j+1));
   mxSetFieldByNumber(S,l, 2, (k+1).toMx());
   mxSetFieldByNumber(S,l, 3, fac.toMx());
};

namespace DY {

mxArray* Symmetry::toMx() const {

   const char *fields[]={"q","R","A","n2","M"};
   mxArray* a=mxCreateStructMatrix(1,1,5,fields);

   mxSetFieldByNumber(a,0,0, q .toMx());
   mxSetFieldByNumber(a,0,1, R .toMx());
   mxSetFieldByNumber(a,0,2, A .toMx());
   mxSetFieldByNumber(a,0,3, n2.toMx());
   mxSetFieldByNumber(a,0,4, M .toMx());

   return a;
};

template <class TQ>
mxArray* Weights<TQ>::toMx() const {

   const char *fields[]={"qq","dd"};
   mxArray *S=mxCreateStructMatrix(1,1,2,fields);

   size_t i=0, N=W.size(), n = (N ? W.begin()->first.len : 0);
   wbMatrix<TQ> qq(N,n);
   wbvector<TQ> dd(N);

   for (auto I=W.begin(); I!=W.end(); ++I, ++i) {
      qq.recSet(i,I->first);
      dd[i]=I->second.m;
   }

   mxSetFieldByNumber(S,0,0, qq.toMx());
   mxSetFieldByNumber(S,0,1, dd.toMx());

   return S;
};

}; 

template <class TQ> 
QSet<TQ>& QSet<TQ>::init(
   const char *F, int L, const mxArray *a, unsigned k
){
   if (!a) {
      init(); return *this;
   }
   if (!mxIsStruct(a)) wblog(FL,
      "ERR %s() invalid QSet (%s)",FCT,mxGetClassName(a));

   unsigned n=mxGetNumberOfElements(a);
   int itype = mxGetFieldNumber(a,"type"),
       iqset = mxGetFieldNumber(a,"qset"),
       iqdir = mxGetFieldNumber(a,"qdir");

   if (itype<0 || iqset<0 || iqdir<0) wblog(F_L,
      "ERR invalid CData structure (%d,%d,%d)",itype,iqset,iqdir);

   if (!n) wblog(FL,"ERR %s() got empty CData structure",FCT);
   if (int(k)>=0) {
      if (k>=n) wblog(FL,"ERR index out of bounds (%d/%d)",k+1,n); }
   else k=0;

   t   .init(FL,mxGetFieldByNumber(a,k,itype));
   qs  .init(FL,mxGetFieldByNumber(a,k,iqset),0,'!'); 
   qdir.init(FL,mxGetFieldByNumber(a,k,iqdir));

   return *this;
};

template <class TQ>
mxArray* QSet<TQ>::toMx() const {

   mxArray *S=mxCreateStruct(1,1);
   add2MxStruct(S,0); return S;

};

template <class TQ>
mxArray* QSet<TQ>::mxCreateStruct(unsigned m, unsigned n) const {

   const char *fields[]={"type","qset","qdir"};
   return mxCreateStructMatrix(m,n,3,fields);
};

template <class TQ>
void QSet<TQ>::add2MxStruct(mxArray *S, unsigned l) const {

   mxSetFieldByNumber(S,l,0, t.toMx());
   mxSetFieldByNumber(S,l,1, qs.toMx());
   mxSetFieldByNumber(S,l,2, qdir.toMx());
};

template <class TD>
cdata<TD>& cdata<TD>::init( 
   const char *F, int L, const mxArray *C, unsigned k,
   const QType &q, unsigned r
){
   if (C==NULL) {
       if (k) wblog(FL,"ERR %s() got k=%d for NULL mxArray",FCT,k);
       cgsparray::init();
       return *this;
   }

   try { cgsparray::init(F_L,C,k); }
   catch (...) { wblog(F_L,
      "ERR %s() invalid cdata (%s)",FCT,mxGetClassName(C));
   }

   const bool isa=q.isAbelian();

   if (isa && !isScalar()) wblog(F_L,
      "ERR %s() CG scalar inconsistency (%d/%d %s)",
       FCT, isScalar(), isa, this->sizeStr().data
   );

   if (int(r)>=0) {
   this->AddTrailingSingletons(FL,r); }

   return *this;
};

CGR_TYPE cgdStatus::init(const char *F, int L, const mxArray* a) {

   CGR_TYPE rt=CGR_DEFAULT;

   if (!a || mxIsEmpty(a)) {
      init(CGD_UNKNOWN); return rt;
   }

   if (mxIsDouble(a)) {
      wbvector<double> cid(F_L,a); int l=3;
      if (cid.len<4 || cid.len>5) wblog(F_L,
         "ERR %s() invalid CRef::cgb data (%d)",FCT,cid.len);

      ctime = cid[0];
      mtime = cid[1];
      ID    = cid[2];

      if (cid.len==5) { 
         rt=CGR_TYPE(cid[l++]);
         if (rt>=CGR_NUM_TYPES) wblog(F_L,
            "ERR %s() rtype out of bounds (%d/%d)",FCT,rt,CGR_NUM_TYPES
         );
      }

      t=CGD_TYPE(cid[l]);
      if (t>=CGD_NUM_TYPES) wblog(F_L,
         "ERR %s() invalid CData type %g",FCT,double(cid[4])
      );
   }
   else if (mxIsStruct(a)) {
      unsigned n=mxGetNumberOfElements(a);
      int ic=mxGetFieldNumber(a,"ctime"), id=mxGetFieldNumber(a,"ID"),
          im=mxGetFieldNumber(a,"mtime"), iflag=mxGetFieldNumber(a,"flag");

      if (ic<0 || im<0 || id<0 || iflag<0 || n!=1) wblog(F_L,
         "ERR %s() unexpected cgdStatus structure (%d %d %d %d; %d)",
         FCT,ic,im,id,iflag,n
      );

      mxGetNumber(mxGetFieldByNumber(a,0,ic), ctime);
      mxGetNumber(mxGetFieldByNumber(a,0,im), mtime);
      mxGetNumber(mxGetFieldByNumber(a,0,id),    ID);
      mxGetNumber(mxGetFieldByNumber(a,0,iflag),  t);

      if (t>=CGD_NUM_TYPES) wblog(F_L,
         "ERR %s() invalid CData type %g",FCT,t
      );
   }
   else wblog(F_L,
     "ERR %s() unexpected cgdStatus (%s)",FCT,mxGetClassName(a));

   return rt;
};

mxArray* cgdStatus::toMx() const {
   double cid[4]={ ctime, mtime, double(ID), double(t) };
   return wbvector<double>(4,cid).toMx();
};

mxArray* cgdStatus::toMx(CGR_TYPE rt) const {
   double cid[5]={ctime, mtime, double(ID), double(rt), double(t) };
   return wbvector<double>(5,cid).toMx();
};

mxArray* cgdStatus::mxCreateStruct(unsigned m, unsigned n) const {

   const char *fields[]={"ctime","mtime","ID","flag"};
   return mxCreateStructMatrix(m,n,4,fields);
};

mxArray* cgdStatus::add2MxStruct(mxArray *S, unsigned k) const {

   mxSetFieldByNumber(S,k,0, numtoMx(ctime));
   mxSetFieldByNumber(S,k,1, numtoMx(mtime));
   mxSetFieldByNumber(S,k,2, numtoMx(ID   ));
   mxSetFieldByNumber(S,k,3,
     t==CGD_UNKNOWN ? wbstring("").toMx() : toStr().toMx());

   return S;
};

int CG::isCData(const mxArray *S, unsigned k) {

   if (!S || !mxIsStruct(S)) { return 0; }

   if (mxGetFieldNumber(S,"type")<0) return -1;
   if (mxGetFieldNumber(S,"qset")<0) return -2;
   if (mxGetFieldNumber(S,"qdir")<0) return -3;
   if (mxGetFieldNumber(S,"cgd" )<0) return -4;
   if (mxGetFieldNumber(S,"cid" )<0) return -5;
   if (int(k)>=0 && k>=mxGetNumberOfElements(S)) return -99;

   return 1;
};

int CG::isQSet(const mxArray *S, unsigned k) {

   if (!S || !mxIsStruct(S)) { return 0; }

   if (mxGetFieldNumber(S,"type")<0) return -1;
   if (mxGetFieldNumber(S,"qset")<0) return -2;
   if (mxGetFieldNumber(S,"qdir")<0) return -3;
   if (int(k)>=0 && k>=mxGetNumberOfElements(S)) return -99;

   return 1;
};

int CG::isCRef(const mxArray *S, unsigned k) {

   if (!S || !mxIsStruct(S)) { return 0; }

   if (mxGetFieldNumber(S,"type")<0) return -1;
   if (mxGetFieldNumber(S,"qset")<0) return -2;
   if (mxGetFieldNumber(S,"qdir")<0) return -3;
   if (mxGetFieldNumber(S,"cid" )<0) return -5;
   if (mxGetFieldNumber(S,"cgw" )<0) return -5;
   if (mxGetFieldNumber(S,"size")<0) return -6;
   if (int(k)>=0 && k>=mxGetNumberOfElements(S)) return -99;

   return 1;
};

template <class TQ, class TD>
int CData<TQ,TD>::Load_CRef(
   const char *F, int L, CData<TQ,TD> &Cr, const char *cgd) const {

   if (!cgd || !cgd[0]) { 
      if (F) wblog(F,L,"ERR %s() got empty file",FCT);
      return -1;
   }

   Wb::matFile f(FL,cgd,"r");
   if (!f.mfp) {
      if (F) wblog(FL,"ERR %s() failed to open file\n%s",FCT,cgd);
      return -2;
   }

   mxArray *a=matGetVariable(f.mfp,"CRef");
   if (!a) wblog(F_L,"ERR %s() missing `CRef' in %s",FCT,cgd);

   Cr.init_mxCRef(F_L,a,0);
   mxDestroyArray(a);

   return 0;
};

template <class TQ, class TD>
int CData<TQ,TD>::init_mxCRef( 
   const char *F, int L, const mxArray *S, unsigned k
) {

   if (gotuser_BUF()) wblog(F_L,"ERR %s() ",FCT);

   if (!S) {
      init(); return 0; 
   }
   if (!mxIsStruct(S)) wblog(FL,
      "ERR %s() invalid CData (%s)",FCT,mxGetClassName(S));

   int itype = mxGetFieldNumber(S,"type"),
       iqset = mxGetFieldNumber(S,"qset"),
       iqdir = mxGetFieldNumber(S,"qdir"),
       icgd  = mxGetFieldNumber(S,"cgd" ),
       id    = mxGetFieldNumber(S,"cid" );
   unsigned n=mxGetNumberOfElements(S);

   if (itype<0 || iqset<0 || iqdir<0 || icgd<0 || id<0) wblog(F_L,
      "ERR invalid CData structure (%d,%d,%d,%d,%d)",
      itype,iqset,iqdir,icgd,id
   );

   if (!n) wblog(FL,"ERR %s() got empty CData structure",FCT);
   if (int(k)>=0) {
      if (k>=n) wblog(FL,"ERR index out of bounds (%d/%d)",k+1,n); }
   else k=0;

   this->t    .init(FL,mxGetFieldByNumber(S,k,itype));
   this->qs   .init(FL,mxGetFieldByNumber(S,k,iqset),0,'!'); 
   this->qdir .init(FL,mxGetFieldByNumber(S,k,iqdir));
   this->cgd  .init(FL,mxGetFieldByNumber(S,k,icgd),0,this->t);
   this->cstat.init(FL,mxGetFieldByNumber(S,k,id));

   return 1; 
};

template <class TQ, class TD>
mxArray* CData<TQ,TD>::toMx() const {

   mxArray *S=mxCreateStruct(1,1);
   add2MxStruct(S,0); return S;
};

template <class TQ, class TD>
mxArray* CData<TQ,TD>::mxCreateStruct(unsigned m, unsigned n) const {

#ifdef CG_CHECK_MW_PERM
   const char *fields[]={"type","qset","qdir","cid","cgd","P0"};
   return mxCreateStructMatrix(m,n,6,fields);
#else
   const char *fields[]={"type","qset","qdir","cid","cgd"};
   return mxCreateStructMatrix(m,n,5,fields);
#endif
};

template <class TQ, class TD>
void CData<TQ,TD>::add2MxStruct(mxArray *S, unsigned i, char tst) const {

   if (tst) {
      unsigned s=0; int k=0; 
      if (S==NULL || (s=mxGetNumberOfElements(S))<1 || i>=s
       || (k=mxGetFieldNumber(S,"qdir"))<0) wblog(FL,
      "ERR %s() missing mxCreateStruct()\n%lx, %d/%d, %d",
       FCT,S,i+1,s,k);
   }
   else if (!S) wblog(FL,"ERR %s() got non-initialized S",FCT);

   mxSetFieldByNumber(S,i,0, this->t.toMx());
   mxSetFieldByNumber(S,i,1, this->qs.toMx());
   mxSetFieldByNumber(S,i,2, this->qdir.toMx());
   mxSetFieldByNumber(S,i,3, cstat.toMx());
   mxSetFieldByNumber(S,i,4, cgd.toMx());

#ifdef CG_CHECK_MW_PERM
   if (!gRS.buf.isEmpty(this->t) && this->qdir=="++-" && rank(FL)==3) {
      unsigned n=this->t.qlen();

      const genRG_base<TQ,RTD> &R =
         gRS.getR(FL,this->t,this->qs.data+2*n);
      mxSetFieldByNumber(S,i,5, R.P0.toMx());
   }
#endif
};

template <class TQ>
CRef<TQ>& CRef<TQ>::init( 
   const char *F, int L, const mxArray *a, unsigned k,
   char refC, const QSet<TQ> *QS,
   char xflag 
){
   if (!a || mxIsEmpty(a)) {
      if (k) wblog(FL,
         "ERR CRef::%s() got k=%d for null mxArray\n%s",FCT,k,
         QS ? QS->toStr().data : "(no QS specified)");
      return init();
   }

   unsigned n=mxGetNumberOfElements(a);

   int lflag=0, e=0,
      ix = -1, 
      it = mxGetFieldNumber(a,"type"), 
      iq = mxGetFieldNumber(a,"qset"), 
      io = mxGetFieldNumber(a,"qdir"), 
      id = mxGetFieldNumber(a,"cid" ), 
      iw = mxGetFieldNumber(a,"cgw" ), 
      is = mxGetFieldNumber(a,"size"); 

   if (it<0 || iq<0 || io<0 || id<0 || iw<0 || is<0) wblog(FL,
      "ERR %s() invalid CRef (%d %d %d; %d %d %d)",FCT,it,iq,io,id,iw,is);
   if (k>=n) wblog(F,L,
      "ERR CRef::%s() index out of bounds (%d/%d)",FCT,k,n);

   const mxArray
      *as = (is>=0 ? mxGetFieldByNumber(a,k,is) : NULL),
      *ad = (id>=0 ? mxGetFieldByNumber(a,k,id) : NULL),
      *aw = (iw>=0 ? mxGetFieldByNumber(a,k,iw) : NULL),
      *at = (it>=0 ? mxGetFieldByNumber(a,k,it) : NULL),
      *aq = (iq>=0 ? mxGetFieldByNumber(a,k,iq) : NULL),
      *ao = (io>=0 ? mxGetFieldByNumber(a,k,io) : NULL);

   wbvector<SPIDX_T> S;
   cgdStatus st;
   QSet<TQ> Q;

   init();

   if (mxIsDouble(aw)) { cgw.init(F_L,aw); } else
   if (mxIsInt8(aw) || mxIsChar(aw)) {
      wbvector<RTD> w(F_L,aw);
      if (w.len!=1) wblog(FL,"ERR %s() got len=%d",FCT,w.len);
      cgw.initScalar(double(w[0]));
   }
   else wblog(FL,"ERR %s() ",FCT);

   if (!at || !aq || !ao) { 
      if (as && mxGetNumberOfElements(as)) e|= 1;
      if (ad && mxGetNumberOfElements(ad)) e|= 2;
      if (at && mxGetNumberOfElements(at)) e|= 4;
      if (aq && mxGetNumberOfElements(aq)) e|= 8;
      if (ao && mxGetNumberOfElements(ao)) e|=16;
      if (e) wblog(FL,"ERR %s() "
         "got invalid CRef (CGR_ABELIAN; k=%d, e=%s)",FCT,k+1,BITS(e));

      if (!QS || !xflag) {
         rtype=CGR_ABELIAN;
         return *this;
      }

      Q=(*QS); lflag=1;
      if (!cgw) { double x=1; cgw.init(1,1,&x); }
   }
   else {
      S     .init(F_L,as);
      st    .init(F_L,ad);
      Q.t   .init(F_L,at);
      Q.qs  .init(F_L,aq,0,'!'); 
      Q.qdir.init(F_L,ao);

      lflag=S.isEmpty();
      if (lflag ^ st.isEmpty()) wblog(FL,
         "ERR %s() got invalid minimal CRef data (%d,%d)\n"
         "(empty size and cid required to indicate initialization)!",
         FCT, lflag, STR2(st,'v')
      );

      if (QS) {
         if (Q!=(*QS)) wblog(FL, 
            "ERR %s() got CGR QSet mismatch\n%s <> %s",
            FCT,STR_(QS),STR(Q)
         );
      }
   }

   if (!Q.t.isAbelian()) {
      if (gStore.setupDirs(FL)) {   
      if (gStore.is_new(Q.t)) { 
         gStore.rclog(Q.t,FL,"");   
      }}
   }

   Q.Sort(&cgp,&conj,'i');
   CDATA_TQ &Cb=(
      refC && !lflag ? gCS.getBUF(0,0,Q,LB_REF) : gCS.getBUF(FL,Q,LB_UPD)
   );

   char gotC=(Cb.isEmpty() ? 0 : (Cb.isRefInit() ? -1 : 1));
   unsigned r=Q.rank(), m=Cb.getOM(), n1=wdim1();

   if (gotC>0 && n1>m) {
      if (refC) {
         if (Cb.Reduce2Ref(FL,'!')>0) { gotC=-1; }
      }
      else {
         wblog(FL,"ERR %s() OM out of bounds (%s /%d; %d,%d,%d)",
         FCT,SSTR(cgw),m, refC,lflag,gotC);
      }
   }

   if (lflag) {
      if (!gotC) wblog(FL,"ERR %s() got minimal CRef "
         "for non-existing CData\n%s",FCT,STR(Q));
      if (Cb!=Q) wblog(FL,"ERR %s() got BUF inconsistency\n"
         "(%s <> %s)",FCT,STR(Cb),STR(Q));
      if (!cgw || n1>Cb.getOM()) wblog(FL,
         "ERR %s() OM out of bounds (%s /%d)",FCT,SSTR(cgw),Cb.getOM());

      cgb=(&Cb); 

      if (CG_VERBOSE>6) wblog(PFL,
         " *  got minimal CRef %s",STR_(this));
      return *this;
   }

   if (st==CGD_UNKNOWN) wblog(FL,
      "ERR %s() got invalid CData status\n%s",FCT,st.toStr('v').data);

   if (Q.t.isAbelian()) {
      if (st.ctime) wblog(FL,
         "WRN %s() got cstat for abelian symmetry\n%s\n%s",
         FCT,STR(Q),st.toStr('V').data
      );
   }
   else if ((e=st.inValid())) wblog(FL,
      "WRN %s() inalid cstat for non-abelian symmetry (e=%d)"
      "\n%s\n%s", FCT,STR(Q),STR2(st,'V'), e
   );

   if (cgp.len && !cgp.isIdentityPerm()) {
      wbvector<SPIDX_T> Sx(S);
      if (cgp.len+1<S.len || cgp.len>S.len) wblog(FL,
         "ERR %s() invalid S.len=%d/%d",FCT,S.len,cgp.len);
      for (unsigned i=0; i<cgp.len; ++i) { Sx.el(cgp[i])=S[i]; }
      Sx.save2(S);
   }

   rtype=CGR_DEFAULT; 

   if (!gotC) { 
      if (!S.len) wblog(FL,"ERR %s() got empty cgb.size",FCT);

      Cb=Q; Cb.cstat=st;

      if (Cb.t.isAbelian()) {
         if (!S.allEqual(1)) wblog(FL,
            "ERR %s() invalid scalar size (%s)",FCT,Cb.cgd.sizeStr().data);
         Cb.cstat.init(CGD_ABELIAN);
         Cb.cgd.wbsparray<RTD>::init();
      }
      else if (refC) {
         Cb.cstat.t=CGD_REF_INIT; 
         ix=mxGetFieldNumber(a,"cgt");
         if (ix<0) Cb.cgd.init(S); 
      }
      else wblog(PFL,
         "ERR %s() got CRef to undefined CData (%s)\n%s",
         FCT,myname,STR(Cb)
      );
   }
   else if (gotC>0) {
      if (Cb.cstat.t==CGD_REF_INIT) wblog(FL, 
         "ERR %s() %s got cstat.t=%s",FCT,STR(Cb),Cb.cstat.tstr());
      if (Cb!=Q) wblog(FL,"ERR %s() got BUF inconsistency\n"
         "(%s <> %s)",FCT,STR(Cb),STR(Q));

      int q=Cb.sameSizeR(FL,S); 

      if (q<=0) wblog(FL,"ERR %s() "
         "CRef size mismatch (q=%d)\n   %s\n-> %s (cgw %s)\n%s",
         FCT,q, SSTR(Cb), SSTR(S), SSTR(cgw), STR(Cb));
      if (r==2 && q!=1) { wblog(FL, "ERR %s() " 
         "CRef OM mismatch (q=%d) for r=2\n   %s\n<> %s (cgw %s)\n%s",
         FCT, q, SSTR(Cb), SSTR(S), SSTR(cgw), STR(Cb));
      }

      if (q==2) {
         if (!refC) wblog(FL,"ERR %s() "
            "CRef with larger OM (q=%d)\n   %s\n-> %s (cgw %s)\n%s",
            FCT, q, SSTR(Cb), SSTR(S), SSTR(cgw), STR(Cb));

         q=Cb.cstat.cmp(0,0,st); 
         if (q!=-1) wblog(FL,"ERR %s() "
            "CRef unexpected cstat (q=%d)\n   %s\n<> %s (cgw %s)\n"
            "%s%N%N CData: %s%N CRef : %s%N",
             FCT, q, SSTR(Cb), SSTR(S), SSTR(cgw), STR(Cb),
             Cb.cstat.toStr('V').data, st.toStr('V').data);
         Cb.Reduce2Ref(FL,'!'); 
      }
      else { q=Cb.cstat.cmp(0,0,st); if (q) {
         if (q<-1 || q>1) e=1; 
         else if (r<=2) e=2;   
         else if (q< 0) {
            if (refC) Cb.Reduce2Ref(FL,'!'); 
            else e=3; 
         }

         if (e) { wblog(PFL,
            "WRN %s() unexpected CData status mismatch (e=%d)\n"
            "   %s\n<> %s%N%N   %s%N-> %s%N", FCT, e,
            STR(Q), STR(Cb), STR2(st,'V'), STR2(Cb.cstat,'V'));

            if (Q.qdir.len>3) wblog(FL,"ERR %s() ",FCT);
         }
      }}
   }
   else { 
      if (Cb.cstat.t!=CGD_REF_INIT || !refC) wblog(FL,
         "ERR %s() got isref CData.cstat='%s' (%d)",
         FCT,STR(Cb.cstat),refC);
      if (Cb!=Q) wblog(FL,"ERR %s() got BUF inconsistency\n"
         "(%d <> %s)",FCT,STR(Cb),STR(Q));

      int q=Cb.sameSizeR(FL,S); 
      if (q<=0) wblog(FL,
         "ERR %s() CRef size mismatch %s -> %s (%s)\n%s",FCT,
         SSTR(Cb), STR(S), SSTR(cgw), STR(Cb));

      if (!refC || st.isSet()) { char wrn=0;
         if (!Cb.cstat.sameAs(st,2)) { wrn=2; } else
         if (q==1) {
            if (!Cb.cstat.sameAs(st,1)) { wrn=1; }
         }

      if (wrn) {
         PRINTF("\n  CRef status mismatch [w=%d] for",wrn);
         if (Q==Cb)
              { PRINTF(" %s",STR(Q)); }
         else { PRINTF("\n     %s\n  <> %s\n", STR((QSet<TQ>&)Cb), STR(Q)); }
         PRINTF("\n     %s\n  <> %s\n", STR2(Cb.cstat,'V'), STR2(st,'V'));

         if (wrn>1 || (Cb.gotuser_BUF() && Cb.cstat.ID!=st.ID)) {
            PRINTF("  hint: got different RC_STORE?\n\n"); }

         wblog(FL, wrn<2? "WRN %s() %s":"ERR %s() %s",FCT,PSTR);
      }}

      if (q==2) {
         if ((q=Cb.cstat.cmp(0,0,st))!=-1) wblog(FL,"ERR %s() "
            "CRef unexpected cstat (q=%d)\n   %s\n<> %s (cgw %s)\n"
            "%s%N%N CData: %s%N CRef : %s%N",
             FCT, q, SSTR(Cb), SSTR(S), SSTR(cgw), STR(Cb),
             Cb.cstat.toStr('V').data, st.toStr('V').data
         );
         if (refC) 
            { if (Cb.cstat.mtime) Cb.cstat.mtime=st.mtime; }
         else wblog(FL,"ERR %s()",FCT); 

         ix=mxGetFieldNumber(a,"cgt");
         if (ix<0) Cb.cgd.init(S); 
      }
   }

   if (ix>=0) {
      const mxArray *ax=mxGetFieldByNumber(a,k,ix);

      if (Cb.cstat.t!=CGD_REF_INIT) wblog(FL, 
         "ERR %s() got non-ref CData\n%s",FCT,STR(Cb));
      if (r && S.len>r) {
         if (Cb.cgd.SIZE.len==S.len && Cb.cgd.SIZE[r]>S[r]) wblog(FL,
            "ERR %s() got reduced OM %d/%d",
            FCT,Cb.sizeStr().data,STR(S)
         );
      }

     #ifndef QS_SKIP_MPFR 
      if (mxIsDouble(ax) && mxGetNumberOfElements(ax)) { 
         static unsigned nlog=0; 

         wbvector<double> dd(F_L,ax);
         unsigned i=0, n=MIN(dd.len,Cb.cgd.D.len); double e,e2=0;
         for (; i<n; ++i) {
            e=dd[i]-double(Cb.cgd.D[i]); e2+=e*e;
         }
         if (e2>1e-24) { 
            MXPut X(FL,"Idbg"); X.add(i,"i").add(n,"n")
              .add(dd,"dd").add(Cb.cgd,"cgD").add(Q,"Q").add(Cb,"Cb");
               X.put("tmpfile"); 

            PRINTF("\ndata: input/cgt   loaded/internal\n");
            for (i=0; i<n; ++i) {
            PRINTF("  %2d: %9.5g %9.5g\n",i+1,dd[i],double(Cb.cgd.D[i])); }

            wblog(FL,"ERR %s() cgt[%d/%d] inconsisteny @ "
            "e2=%.3g\nQ: %s\nC: %s",FCT,i,n,SQRT(e2),STR(Q),STR(Cb));
         }

         if (dd.len>Cb.cgd.D.len) {
            wbvector<RTD> cgt(dd.len);
            for (i=0; i<n; ++i) { cgt[i]=Cb.cgd.D[i]; }
            for (; i<dd.len; ++i) {
               if (++nlog<8)
                    { cgt[i].init_d(FL, dd[i]); }
               else { cgt[i].init_d(0,0,dd[i]); }
            }
            Cb.RefInit_auxtr(FL,S,cgt);
         }
      } else
     #endif
      {  wbvector<RTD> cgt(F_L,ax);
         Cb.RefInit_auxtr(FL,S,cgt);
      }
   }

   cgb=(&Cb);

   n=wdim1();
   if (n>cgb->getOM() && !Cb.isRefInit()) wblog(FL, 
      "ERR %s() OM out of bounds (%s /%d; %d,%d,%d,%d) #%05x\n%s",
      FCT, SSTR(cgw), cgb->getOM(),
      refC, lflag, gotC, ix, cgb->cstat.ID, STR_(this)
   );

   return *this;
};

template <class TQ>
mxArray* CRef<TQ>::toMx(char flag) const {

   mxArray* a=mxCreateStruct(1,1,flag);
   this->add2MxStruct(a,0,flag);
   return a;
};

template <class TQ>
mxArray* CRef<TQ>::mxCreateStruct(
   unsigned m, unsigned n, char flag) const {

   if (m && n) {
      if (flag==0) { 
         const char *fields[]={
            "type","qset","qdir","cid","size","nnz","cgw","cgt"};
         return mxCreateStructMatrix(m,n,8,fields);
      }
      else {
         const char *fields[]={
            "type","qset","qdir","cgs","size","nnz","cgw"};
         return mxCreateStructMatrix(m,n,7,fields);
      }
   }
   else return NULL;
};

template <class TQ>
void CRef<TQ>::add2MxStruct(mxArray *S, unsigned k, char flag) const {

   if (cgb) {
      wbvector<size_t> S2; getSize(S2);
      if (!S2.len) wblog(FL,"WRN %s() got empty S2",FCT);

      if (!cgb->QSet<TQ>::isSorted()) wblog(FL,
         "ERR %s() got non-sorted CData\n%s",FCT,STR_(cgb));

      if (anyTrafo()) {
         QSet<TQ> QS; 
         if (!isSortedDegQ(&QS)) wblog(FL,
            "ERR %s() got non-standard CRef\n%s",FCT,STR_(this));

         mxSetFieldByNumber(S,k,0, QS.t   .toMx('t'));
         mxSetFieldByNumber(S,k,1, QS.qs  .toMx());
         mxSetFieldByNumber(S,k,2, QS.qdir.toMx());
         mxSetFieldByNumber(S,k,4, S2     .toMx());
      }
      else {
         const CDATA_TQ &C=*cgb;
         mxSetFieldByNumber(S,k,0, C.t    .toMx('t'));
         mxSetFieldByNumber(S,k,1, C.qs   .toMx());
         mxSetFieldByNumber(S,k,2, C.qdir .toMx());
         mxSetFieldByNumber(S,k,4, S2     .toMx());
      }

      if (!flag) { 
         wbvector<double> cgt; 
         if (rtype.t==CGR_CTR_ZERO) wblog(FL,
            "WRN CRef::toMx() got %s rtype=%s",STR_(this),STR(rtype));

         mxSetFieldByNumber(S,k,3, cgb->cstat.toMx(rtype.t));

         mxSetFieldByNumber(S,k,7, this->trace(0,0,cgt).toMx());
      }
      else {
         mxSetFieldByNumber(S,k,3, cgb->cgd.toMx());
      }
   }

   mxSetFieldByNumber(S,k,5, numtoMx(
      !isRefInit() ? (
        double(cgb ? cgb->cgd.nnz() : (cgw ? 1 : 0))
      ) : -1.
   ));

   mxSetFieldByNumber(S,k,6, cgw.toMx());

};

template <class TQ, class TD>
mxArray* genRG_base<TQ,TD>::mxCreateStruct(unsigned m, unsigned n) const {
#ifdef CG_CHECK_MW_PERM
   const char *fields[]={"type","J","Z","Sp","Sz","P"};
   return mxCreateStructMatrix(m,n,6,fields);
#else
   const char *fields[]={"type","J","Z","Sp","Sz"};
   return mxCreateStructMatrix(m,n,5,fields);
#endif
};

template <class TQ, class TD>
void genRG_base<TQ,TD>::add2MxStruct(mxArray *S, unsigned i) const {

   if (q.type==QTYPE_UNKNOWN) {
      if (J.isEmpty())
         mxSetFieldByNumber(S,i,0, wbstring().toMx());
      else {
         mxSetFieldByNumber(S,i,0, q.toMx());
      }
   }
   else mxSetFieldByNumber(S,i,0, q.toMx());

   mxSetFieldByNumber(S,i,1, J.toMx());
   mxSetFieldByNumber(S,i,2, Z.toMx());
   mxSetFieldByNumber(S,i,3, Sp.toMx());
   mxSetFieldByNumber(S,i,4, Sz.toMx());

#ifdef CG_CHECK_MW_PERM
   mxSetFieldByNumber(S,i,5, P0.toMx());
#endif
};

template <class TQ, class TD>
void genRG_struct<TQ,TD>::put(const char *F, int L,
   const char *vname, const char *ws
) const {
   mxArray *a=toMx(); int i=mexPutVariable(ws,vname,a);

   if (i) wblog(F_L,
      "ERR failed to write variable `%s' (%d)",vname,i);
   if (F) wblog(F,L,"I/O putting variable '%s' to %s",vname,ws);

   mxDestroyArray(a);
};

template <class TQ, class TD>
mxArray* genRG_struct<TQ,TD>::toMx(char bflag) const {

   if (!bflag) { 
      mxArray *a=MXPut().add(q,"type")
      .add(DZ,"DZ").add(CR,"CR").add(qdef,"qdef").toMx();
      return a;
   }

   unsigned i=0, n=RSet.size();
   const char* fields[]={"type","info","store"};
   mxArray *S, *a;

   if (!q.isKnown()) wblog(FL,"ERR %s() "
      "got invalid sym=%s (size=%d)",FCT,STR(q),RSet.size());

   S=mxCreateStructMatrix(1,1,3,fields);
   a=MXPut().add(DZ,"DZ").add(CR,"CR").add(qdef,"qdef").toMx();

   mxSetFieldByNumber(S,0, 0, q.toMx());
   mxSetFieldByNumber(S,0, 1, a);

   if (n) { i=0;
      for (auto I=RSet.begin(); I!=RSet.end(); ++I, ++i) {
         if (i==0) { a=I->second.mxCreateStruct(n,1); }
         I->second.add2MxStruct(a,i);
      }
      mxSetFieldByNumber(S,0, 2, a);
   }

   return S;
};

template <class TQ>
unsigned CStore<TQ>::add3(const QType &q, 
   const mxArray* S, 
   char loadC
){
   if (!S || !mxIsStruct(S)) wblog(FL,
      "ERR %s() invalid CStore.%s.std (%s <> cell)",
      FCT,q.toStr('t').data, S ? mxGetClassName(S) : "(null)");

   unsigned j=0, m=q.qlen(), n=mxGetNumberOfElements(S);
   double x=1.; wbarray<double> E;

   qset<TQ> J, J12;
   CRef<TQ> Rj;
   QDir qd;

   if (!n) wblog(FL,"ERR got empty CStore data for %s",STR2(q,'t'));

   int itype = mxGetFieldNumber(S,"type"),
       iqset = mxGetFieldNumber(S,"qset"),
       iqdir = mxGetFieldNumber(S,"qdir");

   if (itype<0 || iqset<0 || iqdir<0) wblog(FL,"ERR invalid "
      "CStore/%s/std (%d,%d,%d)",q.toStr('t').data, itype,iqset,iqdir);
   else {
      QSet<TQ> Q;
      Q.t.   init(FL,mxGetFieldByNumber(S,j,itype));
      Q.qs.  init(FL,mxGetFieldByNumber(S,j,iqset),0,'!'); 
      Q.qdir.init(FL,mxGetFieldByNumber(S,j,iqdir));

      if (Q.qs.len%3 || Q.qs.len/3!=m || Q.t!=q) wblog(FL,
         "ERR %s() invalid %s CRef data (%d/3*%d; %s)",FCT,
         q.toStr('t').data, Q.qs.len, m, Q.t.toStr('t').data);

      J12.init(2*m,Q.qs.data);
      qd=Q.qdir;
   }

   map < qset<TQ>, CRef<TQ> > M3, *M3x=NULL; 
   {  auto it1 = map3.find(q); 
      if (it1!=map3.end()) {
         auto it2 = it1->second.find(J12);
         if (it2!=it1->second.end()) { M3x=&(it2->second);
            if (n!=M3x->size()) wblog(FL,
            "ERR %s() size mismatch (n=%d/%d)",FCT,n,M3x->size());
         }
      }
   }

   for (; j<n; ++j) { 
      Rj.init(FL,S,j, !loadC); 
      if (!Rj.cgb) wblog(FL,   
         "ERR %s() missing cgb reference (%d/%d)",FCT,j+1,n);

      const QSet<TQ> &Q=(const QSet<TQ>&)(*Rj.cgb);
      J.init(m,Q.qs.data+2*m);

      if (j) {
         if (Q.t!=q) wblog(FL,
         "ERR %s() qtype mismatch (%s / %s)",FCT,STR2(Q.t,'t'),STR2(q,'t'));
         if (Q.qdir!=qd) wblog(FL,
         "ERR %s() qdir mismatch (%s / %s)",FCT,STR(Q.qdir),STR(qd));
      }

      if (!M3x) {  
         M3[J]=Rj; 

         Wb::MatProd(Rj.cgw,Rj.cgw,E,'T');
         if (!E.isProptoId(&x,CG_SKIP_DEPS1)) {
            wblog(FL,"ERR %s() got non-orthonormal CRef %s "
            "(x=%g)\n-> %s",FCT, STR(Rj), x, STR(Rj.cgb->cstat));
         }
      }
      else {
         const CRef<TQ> &Xj = (*M3x)[J];

         if (Xj.cgb!=Rj.cgb) wblog(FL,"ERR %s() "
            "cgb mismatch (%d/%d: 0x%lX, 0x%lX)",FCT,j+1,n,Xj.cgb,Rj.cgb);

         if (Rj!=Xj) { 
            unsigned e=0;
            Wb::MatProd(Rj.cgw,Xj.cgw,E,'T');
            if (!E.isProptoId(x,CG_SKIP_DEPS1)) { e|=1; }
            if (!Rj.isw2() || !Xj.isw2()) { e|=2; } else
            if (!Rj.cgw.sameSize(Xj.cgw,0)) { e|=4; }

            if (e) wblog(FL,"ERR %s() got CRef mismatch (e=%d; x=%g)\n"
              "(%s <> %s)\n   %s\n-> %s)",FCT,e,x, STR(Rj), STR(Xj),
               STR(Rj.cgb->cstat), STR(Xj.cgb->cstat)
            );
         }
      }
   }

   if (n!=M3.size()) wblog(FL,
      "ERR %s() size mismatch (n=%d/%d)",FCT,n,M3.size());

   if (!M3x) {
      map3[q][J12].swap(M3); 
      if ((j=M3.size())) wblog(FL,
         "ERR %s() got Std3 race condition (%s; %d/%d/%d)",
         FCT,STR(J12), j,n, map3[q][J12].size()
      );
   }
   else {
      n=0; 
   }

   return n;
};

template <class TQ>
mxArray* CStore<TQ>::toMx(const QType &q_) const {

   const char *field0[]={"std"};

   mxArray *c,*S, *C=mxCreateStructMatrix(1,1,0,NULL);

   CDATA_TQ X_; 

   unsigned i=0,j,j12;
   mxArray *a;

   for (auto I0=map3.begin(); I0!=map3.end(); ++I0) {
      const QType &q = I0->first;
      const auto &M12 = I0->second;
   if (!q_.isKnown() || q_==q) {

      S=mxCreateStructMatrix(1,1,1,field0);
      c=map3_CreateStructMatrix(1,M12.size()); j12=0;

      for (auto I12=M12.begin(); I12!=M12.end(); ++I12, ++j12) {
         map3_add2MxStruct(c,j12, q, I12->first); 
      }; mxSetFieldByNumber(S,0,0,c);

      map<unsigned,    
         map <qset<TQ>,    
            const CDATA_TQ* >    
      > M;

      for (auto it=BUF.begin(); it!=BUF.end(); ++it) {
         if (it->first.t==I0->first) {
            const qset<TQ> &qs=it->first.qs;
            if (qs!=it->second.qs) wblog(FL,
               "ERR %s() got inconsistent key [%s] <> [%s]",
               FCT,STR(qs),it->second.QStrS().data
            );
            M[it->second.rank(FL)][qs] = &(it->second);
         }
      }

      for (auto i1=M.begin(); i1!=M.end(); ++i1) { char s[8]; j=0;
         a=X_.mxCreateStruct(1,i1->second.size());
         for (auto i2=i1->second.begin(); i2!=i1->second.end(); ++i2) {
            i2->second->add2MxStruct(a,j++);
         }
         snprintf(s,8,"rank%d",i1->first);
         mxAddField2Scalar(FL,S,s,a);
      }
      mxAddField2Scalar(FL,C,q.toStr('t').data,S);
   }}

   wbvector<unsigned> nc(BUF.bucket_count());
   for (; i<nc.len; ++i) nc[i]=BUF.bucket_size(i);

   mxAddField2Scalar(FL,C,"QHash", MXPut(0,0)
      .add(BUF.max_load_factor(),"max_load_factor")
      .add(BUF.load_factor(),"load_factor")
      .add(nc,"filling_buckets")
      .add(BUF.size(),"nr_entries")
      .add(nc.max_(0),"max_collisions")
   .toMx());

   return C;
};

template <class TQ, class TD>
mxArray* RStore<TQ,TD>::toMx(const QType &q_, char bflag) const {

   char aflag=(q_.isKnown() ? 0 : 1); 
   mxArray *c, *C=NULL;

   if (aflag) C=mxCreateStructMatrix(1,1,0,NULL);

   for (auto it=buf.begin(); it!=buf.end(); ++it) {
      if (aflag || q_==it->first) {
         const QType &q = it->first;
         const genRG_struct<TQ,TD> &B = it->second;

         if (!B.q.isKnown()) wblog(FL,
            "ERR %s() got sym=%s / %s @ size=%d",
            FCT, STR(B.q), STR(it->first), B.RSet.size()
         );

         c=B.toMx(bflag);
         if (aflag)
              { mxAddField2Scalar(FL,C,q.toStr('t').data,c); }
         else { C=c; break; }
      }
   }
   return C;
};

template <class TQ, class TD>
unsigned RStore<TQ,TD>::add(
   const QType &q, const mxArray* S, unsigned k
){
   static int itype=-1, iJ=-1, iZ=-1, iSp=-1, iSz=-1;
   const mxArray *ap, *az;

   if (!S || !mxIsStruct(S)) wblog(FL,
      "ERR %s() invalid RStore.%s.store",FCT,q.toStr('t').data);
   if (k>=mxGetNumberOfElements(S)) wblog(FL,
      "ERR RStore.%s.store index out of bounds (%d/%d)",
      q.toStr('t').data, k+1, mxGetNumberOfElements(S));

   if (k==0) {
      int id[5]={ mxGetFieldNumber(S,"type"),
         mxGetFieldNumber(S,"J"),  mxGetFieldNumber(S,"Z"),
         mxGetFieldNumber(S,"Sp"), mxGetFieldNumber(S,"Sz")
      };
      if (itype<0) {
         if ((itype = id[0])<0 || (iJ = id[1])<0 || (iZ = id[2])<0 ||
             (iSp = id[3])<0 || (iSz = id[4])<0) wblog(FL,
            "ERR invalid RStore.%s.store (%d,%d,%d,%d,%d)",
             q.toStr('t').data, itype,iJ,iZ,iSp,iSz
         );
      }
      else {
         if (itype!=id[0] || iJ!=id[1] || iZ!=id[2] ||
             iSp!=id[3] || iSz!=id[4]) wblog(FL,
            "ERR invalid RStore.%s.store (%d/%d,%d/%d,%d/%d,%d/%d,%d/%d)",
             q.toStr('t').data, itype, id[0], iJ, id[1], iZ, id[2],
             iSp, id[3], iSz, id[4]
         );
      }
   }

   unsigned i;
   char xflag=0; 

   qset<TQ> qs(FL,mxGetFieldByNumber(S,k,iJ),0,'!'); 
   QType qk(FL,mxGetFieldByNumber(S,k,itype));

   genRG_base<TQ,TD> X; 

   if (qs.len!=q.qlen() || qk!=q) wblog(FL,
      "ERR %s() invalid %s RSet (%d/%d; %s)",
      FCT, q.toStr('t').data, qs.len, q.qlen(), qk.toStr('t').data
   );

   if (CG_VERBOSE>8) {
      QSet<TQ> Q; Q.init1(q,qs.data);
      wblog(PFL,"(+) RBUF[%03d] %s", gRS.buf[q].RSet.size(), STR(Q));
   }

   X.q=q;
   X.J=qs;
   X.Z.init(FL, mxGetFieldByNumber(S,k,iZ));

   if (!(ap=mxGetFieldByNumber(S,k,iSp)) || mxGetNumberOfElements(ap)!=qs.len)
      wblog(FL,"ERR %s() invalid %s RSet(%d)->Sp",FCT,q.toStr('t').data,k);
   if (!(az=mxGetFieldByNumber(S,k,iSz)) || mxGetNumberOfElements(az)!=qs.len)
      wblog(FL,"ERR %s() invalid %s RSet(%d)->Sz",FCT,q.toStr('t').data,k);

   X.Sp.init(qs.len);
   X.Sz.init(qs.len);

   for (i=0; i<qs.len; ++i) {
      X.Sp[i].init(FL,ap,i);
      X.Sz[i].init(FL,az,i);
   }

   genRG_base<TQ,TD> &R = Buf(q).RSet[qs];

   if (R.isEmpty()) { X.save2(R);
      gCS.getIdentityC(FL,q,qs.data); 
   }
   else { xflag=1;
      if (R.q!=X.q || R.J!=X.J || R.Z!=X.Z || R.Sp!=X.Sp || R.Sz!=X.Sz) {
         double e2=0, E2=0; char s[32];
         MXPut(FL,"Irx").add(R,"R").add(X,"X");

         if (R.q !=X.q ) { E2+=1; wblog(FL,"TST %s() "
            "q : [%s] <> [%s]",FCT, STR(R.q), STR(X.q)); }
         if (R.J !=X.J ) { E2+=R.J.normDiff2(X.J); wblog(FL,"TST %s() "
            "J : [%s] <> [%s]",FCT, STR(R.J), STR(X.J)); }
         if (R.Z !=X.Z ) { E2+=(e2=R.Z.normDiff2(X.Z));
             wblog(FL,"TST %s() Z : diff @ %g",FCT, sqrt(e2)); }
         if (R.Sp!=X.Sp) { E2+=(e2=R.Sp.normDiff2_(X.Sp));
             wblog(FL,"TST %s() Sp: diff @ %g",FCT, sqrt(e2)); }
         if (R.Sz!=X.Sz) { E2+=(e2=R.Sz.normDiff2_(X.Sz));
             wblog(FL,"TST %s() Sz: diff @ %g",FCT, sqrt(e2)); }

         E2=sqrt(E2);
         snprintf(s,32,"q=(%s) (@ %.3g)",STR2(R.J,q),E2);
         if (E2>1E-10)
              wblog(FL,"ERR %s() %s",FCT,s);
         else wblog(FL,"WRN %s() %s",FCT,s);
      }
   }

   gstatR.gotread(R.memSize());

   return (xflag ? 0 : 1);
};

template <class TQ, class TD>
mxArray* x3map<TQ,TD>::mxCreateStruct(unsigned m, unsigned n) const {

   const char *fields[]={ "idc", 
      "a","ica","b","icb","c", "x3",
      "Pab","conj","rtype"
   };
   return mxCreateStructMatrix(m,n,10,fields);
};

template <class TQ, class TD> 
void x3map<TQ,TD>::add2MxStruct(
   mxArray *S, unsigned l, const cgc_contract_id<MTI> *idc) const {

   if (!S) wblog(FL,"ERR %s() got null mxArray",FCT);
   if (l>=mxGetNumberOfElements(S)) wblog(FL,"ERR %s() "
      "index out of bounds (%d/%d)",FCT,l,mxGetNumberOfElements(S));

   ctrIdx ica,icb;

   wbvector<size_t> s3;
   if (cgb) cgb->getSize(s3,'b'); else s3.init();

   #ifndef WB_SKIP_ASSERT
   a.rank(FL); b.rank(FL); c.rank(FL);
   #endif

   if (idc) { CData<TQ,RTD> a_,b_;
   #ifdef QS_USING_OMP
      CG::Guard xLK(FL,*idc);
   #endif

      gXS.contains(FL,*this,*idc); 

      idc->extract(a_,ica,b_,icb, 0);

     #ifndef WB_SKIP_ASSERT
      if (a_!=(QSet<TQ>&)a || b_!=(QSet<TQ>&)b) wblog(FL,
         "ERR %s() QSet mismatch\n   %s | %s\n<> %s | %s",FCT,
         STR(a_), STR(b_), STR(a), STR(b)
      );
     #endif
   }

   mxSetFieldByNumber(S,l,1, a.toMx());
   mxSetFieldByNumber(S,l,3, b.toMx());

   if (idc) {
      mxSetFieldByNumber(S,l,0, idc->toMx());
      mxSetFieldByNumber(S,l,2, ica.toMx()); 
      mxSetFieldByNumber(S,l,4, icb.toMx());
   }

   if (rtype!=CGR_CTR_ZERO) {
      mxSetFieldByNumber(S,l,5, c .toMx());
   }

#ifdef QS_USING_MPFR
   if (x3.SIZE!=X3.SIZE) wblog(FL,
      "ERR %s() size mismatch x3: %s <> X3: %s",FCT,SSTR(x3),SSTR(X3));
   mxSetFieldByNumber(S,l,6, X3.toMx()); 
#else
   mxSetFieldByNumber(S,l,6, x3.toMx());
#endif

   mxSetFieldByNumber(S,l,7, pab.toMx());
   mxSetFieldByNumber(S,l,8, numtoMx(conj));
   mxSetFieldByNumber(S,l,9, numtoMx(rtype.t));
};

template <class TQ, class TD>
x3map<TQ,TD>& x3map<TQ,TD>::init( 
   const char *F, int L, const mxArray* S, unsigned k,
   ctrIdx *ica_, ctrIdx *icb_, char recalc
){
   int n, q=-99;
   ctrIdx ica, icb; 

   if (!S) wblog(FL,"ERR %s() got null mxArray",FCT);
   if ((n=mxGetNumberOfElements(S))!=1) wblog(FL,
      "ERR %s() got invalid number of elements (%d)",FCT,n);

   a  .init_mxCRef(FL, mxGetFieldByNumber(S,k,1), 0);
   ica.init(FL, mxGetFieldByNumber(S,k,2));
   b  .init_mxCRef(FL, mxGetFieldByNumber(S,k,3), 0);
   icb.init(FL, mxGetFieldByNumber(S,k,4));
   c  .init_mxCRef(FL, mxGetFieldByNumber(S,k,5), 0);

    X3.init(FL, mxGetFieldByNumber(S,k,6));
    x3.initT(X3);

   pab.init(FL, mxGetFieldByNumber(S,k,7), 1); 

   mxGetNumber(mxGetFieldByNumber(S,k,8), conj);
   mxGetNumber(mxGetFieldByNumber(S,k,9), rtype.t);

   if (rtype==CGR_CTR_ZERO || rtype==CGR_CTR_SCALAR) {
      if (!c.qdir.isEmpty() || !c.qs.isEmpty()) wblog(FL,
     "ERR %s() got %s (%s)",FCT,STR(c),rtype.tostr());
      if (c.t.isUnknown()) c.t=a.t;
   }
   else {
      if (rtype!=CGR_DEFAULT) wblog(F_L,
     "WRN %s() got rtype=%s",FCT,rtype.tostr());
   }

   gstatX.gotread(this->memSize());

   n=0;
   if (rtype==CGR_CTR_ZERO) {
      if (!x3.isEmpty()) wblog(FL,"ERR %s()",FCT); }
   else {
      if (x3.SIZE.len!=3)
         wblog(FL,"ERR %s() [%s]",FCT,SSTR(x3));
      if (rtype==CGR_DEFAULT) n=1; else
      if (rtype!=CGR_CTR_SCALAR) {
         wblog(FL,"ERR %s() [%s]",FCT,rtype.tostr());
      }
   }

   if (n) {
      CDATA_TQ &Cb = gCS.getBUF(0,0,(QSet<TQ>&)c,LB_REF); 

      rtype=CGR_DEFAULT; 
      cgb=&Cb;

      if (Cb.isEmpty()) {
         if (c.cstat!=CGD_REF_INIT) wblog(FL, 
            "ERR %s() unexpected cstat (%s)",FCT,STR(Cb.cstat));
         Cb=c; 
      }
      else if (Cb.cstat.sameID(c.cstat)) {
         q=Cb.cstat.cmp(FL,c.cstat);
         if (q<0) {
            if (Cb.cstat==CGD_REF_INIT) {
               Cb=c;
            }
            else {
              Cb=c; 
            }
         }
      }
      else {

         gStore.rclog(Cb.t, PF_L,1,
            "WRN X3Map::%s() got ID mismatch (#%05x <> #%05x)"
            "%N  a: %s @ %s%N  b: %s @ %s%N  c: %-26s %s%N"
            "?Cb? %-26s %s%N", FCT, c.cstat.ID, Cb.cstat.ID,
            STR(a),STR(ica), STR(b),STR(icb),
            STR(c),STR2(c.cstat,'V'), STR(Cb),STR2(Cb.cstat,'V'));

         if (!c.isRefInit()) wblog(FL,
            "ERR %s() got non-REF x3map->c data",FCT,STR(c));

         if (Cb.qdir.len==2) { 
            if (Cb.cstat.ctime!=c.cstat.ctime) {
               wblog(FL,"WRN %s() got cstat mismatch",FCT);
               c.init(); c=Cb;

               CRef<TQ> A(a,0), B(b,0);
               cgc_contract_id<MTI> idc(*A.cgb,ica, *B.cgb,icb);
               gStore.save_XMap(FL,idc);
               return *this;
            }
         }
         else if (c.qdir.len==3) { 
            if (recalc==1 && (CG_FIXIT & FIX_CID3)) {
               CRef<TQ> A(a,0), B(b,0);
               wbperm cgp;

               if (Cb.isRefInit()) { Cb.LoadRef(FL); }

               wblog(FL,
                 "TST X3Map::%s() %s [recalc=%d]",FCT,STR(Cb),recalc);
               c.cstat.init();

               this->initCtr(FL,a,ica,b,icb,cgp,'!'); 

               if (zflag) wblog(FL,
                  "ERR %s() got CTR_ZERO (%d)",FCT,zflag);

               cgc_contract_id<MTI> idc(*A.cgb,ica, *B.cgb,icb);

               gXS.contains(FL,*this,idc);

               gStore.save_XMap(FL,idc);
               return *this;
            }
            wblog(FL,"ERR %s() recalc=%d",FCT,recalc);
         }
         else if (c.qdir.len!=2) {
            q=Cb.cstat.cmp(FL,c.cstat);  
            wblog(FL,"ERR %s() q=%d",FCT,q); 
         }
         else wblog(FL,"ERR %s() q=%d",FCT,q); 
      }
   }

   if (ica_) { ica.save2(*ica_); }
   if (icb_) { icb.save2(*icb_); }

   if (CG_VERBOSE>8) wblog(PF_L, 
      "(+) XBUF[%03d] %s\n    %s @ %s\n    %s @ %s\n"
      "--> %s", gXS.XBUF.size(), x3Str().data,
      STR(a), STR(ica), STR(b), STR(icb), STR(c)
   );

   return *this;
};

template <class TQ, class TD>
mxArray* X3Map<TQ,TD>::toMx() const {

   mxArray *S1=NULL;

   map<QType,       
      map <wbvector<char>, 
         map <wbvector<unsigned>,  
            const cgc_contract_id<MTI>*    
   > > > X1;

   wbvector<char> r3(3);
   wbvector<unsigned> s3x; wbvector<SPIDX_T> sc;
   unsigned i=0; int ep=0;

#ifdef QS_USING_OMP
   Wb::ompGuard xLK(XS_buf); 
#endif

   try {

   for (auto it=XBUF.begin(); it!=XBUF.end(); ++it, ++i) {
      const x3map<TQ,TD> &M=it->second;  

      r3[0]=M.a.rank(FL); 
      r3[1]=M.b.rank(FL); 
      r3[2]=M.c.rank(FL); 

      if (M.cgb) M.cgb->getSize(sc,'b'); else sc.init();

      const wbvector<SPIDX_T> &sa=M.a.cgd.SIZE, &sb=M.b.cgd.SIZE;
      s3x.Cat(sa.data,sa.len, sb.data,sb.len, sc.data,sc.len, &i,1);

      if (M.c.t.isUnknown() && (!M.c.isEmpty() || M.rtype!=CGR_CTR_ZERO))
         wblog(FL,"WRN %s()\n%s\n%s\n%s",FCT,STR(M.a),STR(M.b),STR(M.c));
      X1[M.c.t][r3][s3x]=&(it->first);
   }

   S1=mxCreateStructMatrix(1,1,0,NULL);
   char s_[16];

   for (auto it1=X1.begin(); it1!=X1.end(); ++it1) {
      const auto &X2=it1->second;

      mxArray *S2=mxCreateStructMatrix(1,1,0,NULL);

   for (auto it2=X2.begin(); it2!=X2.end(); ++it2) {
      const auto &X3=it2->second;

      const wbvector<char> &r3=it2->first;
      mxArray *S3=x3map<TQ,TD>().mxCreateStruct(X3.size(),1);
      unsigned l=0;

      for (auto it3=X3.begin(); it3!=X3.end(); ++it3, ++l) {

         const cgc_contract_id<MTI> &idc=*(it3->second); 
         const auto itM=XBUF.find(idc);
         if (itM==XBUF.end()) wblog(FL,"ERR %s() key not found",FCT);
         const x3map<TQ,TD> &M=itM->second;          

         ctrIdx ica,icb;
         CData<TQ,RTD> a_,b_;

         idc.extract(a_,ica,b_,icb, 0);

         if (a_!=(QSet<TQ>&)M.a || b_!=(QSet<TQ>&)M.b) wblog(FL,
            "ERR %s() QSet mismatch\n   %s | %s\n<> %s | %s",FCT,
            STR(a_), STR(b_), STR(M.a), STR(M.b)
         );

         M.add2MxStruct(S3,l,&idc);
      }

      snprintf(s_,16,"X%d%d%d",r3[0],r3[1],r3[2]);
      mxAddField2Scalar(FL,S2,s_,S3);
   }
      snprintf(s_,16,"%s", 
         it1->first.isKnown() ? it1->first.toStr('t').data : "ctr_zero");
      mxAddField2Scalar(FL,S1,s_,S2);
   }

   wbvector<unsigned> nc(XBUF.bucket_count());
   for (i=0; i<nc.len; ++i) nc[i]=XBUF.bucket_size(i);

   mxAddField2Scalar(FL,S1,"MHash", MXPut(0,0)
      .add(XBUF.max_load_factor(),"max_load_factor")
      .add(XBUF.load_factor(),"load_factor")
      .add(nc,"filling_buckets")
      .add(XBUF.size(),"nr_entries")
      .add(nc.max_(0),"max_collisions")
   .toMx());

   } catch (...) { ++ep; } 

   if (ep) {
      throw Wb::LogException(ERR);
   }

   return S1;
};

template <class TQ, class TD>
unsigned X3Map<TQ,TD>::add(
   const char *F, int L, const cgc_contract_id<MTI> &idc,
   const mxArray* S, unsigned k
){

#ifdef QS_USING_OMP
   CG::Guard xLK(FL,idc); 
#endif
   x3map<TQ,TD> &Cm=getXBUF(0,0,idc); 

   if (!Cm.a.isEmpty() || !Cm.b.isEmpty() || !Cm.x3.isEmpty())
      wblog(FL,"ERR %s() already got non-empty XBUF entry",FCT);

   Cm.init(F_L,S,k);

#ifndef WB_SKIP_ASSERT
 { cgc_contract_id<MTI> id2;
   id2.init(FL, mxGetFieldByNumber(S,k,0));
   if (idc!=id2) {
      idc.wblog_(FL);
      id2.wblog_(FL);
      wblog(FL,"ERR %s() got contract ID mismatch",FCT);
   }
 }
#endif

   return 1;
};

mxArray* map3_CreateStructMatrix(unsigned m, unsigned n) {
   const char *fields[]={"J12","J","omult","cgr"};
   return mxCreateStructMatrix(m,n,4,fields);
};

template <class TQ>
void map3_add2MxStruct(
   mxArray *S, unsigned k, const QType &q, const qset<TQ> &J12) {

   const auto &M3=gCS.map3[q][J12];
   wbMatrix<TQ> JJ(M3.size(),J12.len/2);
   wbvector<unsigned> OM(M3.size());

   unsigned i=0;
   CRef<TQ> R_; 

   mxArray *a=R_.mxCreateStruct(OM.len,1);

   for (auto I3=M3.begin(); I3!=M3.end(); ++I3, ++i) {
      const qset<TQ> &J=I3->first;
      const CRef<TQ> &c3=I3->second;
      if (!c3.cgw) wblog(FL,"ERR %s() got empty cg3",FCT);

      if (c3.got3(FL,q,J12,J)!=0) wblog(FL, 
         "ERR %s() qset mismatch (%s: %s %s <> %s: %s; %d)",
         FCT, STR(q), STR(J12), STR(J),
         c3.qStr().data, c3.QStr().data, c3.got3(FL,q,J12,J)
      );

      OM[i]=c3.wdim2();
      JJ.recSet(i,(wbvector<TQ>&)J);
      c3.add2MxStruct(a,i);
   }

   mxSetFieldByNumber(S,k,0,J12.toMx());
   mxSetFieldByNumber(S,k,1,JJ.toMx());
   mxSetFieldByNumber(S,k,2,OM.toMx());
   mxSetFieldByNumber(S,k,3,a);
};

int CG::FileLock::tryGetLock(const char *F, int L, char flag_) {

   int i; 
   char flag=tolower(flag_), wait=(flag_<='Z' ? 1 : 0);
   struct flock fl;   

   if (!fname.data) {
      if (fid) {
           wblog(F_L,"ERR %s() got fname=null with fid=%d",FCT,fid); fid=0; }
      else wblog(F_L,"WRN %s() got fname=null with flag=%c",FCT,flag_);
      return fid;
   }

   memset(&fl,0,sizeof(fl));

      fl.l_whence = SEEK_SET; 
      fl.l_start = 0;         
      fl.l_len = 0;           
      fl.l_pid = getpid();    

   switch (flag) {
      case 'r': fl.l_type=F_RDLCK; break; 
      case 'w': fl.l_type=F_WRLCK; break; 
      case 'u': fl.l_type=F_UNLCK; break; 
      default : wblog(F_L,"ERR %s() invalid flag %s\n",FCT,cSTR(flag_));
   }

   if (flag=='u') { 
      if (fid>0) {  
         if ((i=fcntl(fid, F_SETLK,&fl))>=0) {
            gStore.rclog(q, PF_L, CG_VERBOSE>6 && F,
              "RCL flock released (#%X/%c, i=%d; %s)\n%s",fid,flag_,i,
               Wb::ompID2Str('v').data, 
               Wb::repHome(fname).data
            );
         }
         else { wblog(FL,
            "ERR %s() failed to unlock file (#%X/%c; e=%d)\n%s",
            FCT,fid,flag_,i,Wb::repHome(fname).data);
         }
         close(fid); fid=0;
      }
      else if (fid) { 
         wblog(F_L,"WRN cannot unlock file that is not locked "
           "(#%X/%c)\n%s",fid,flag_,fname.data);
         fid=0;
      }
      return fid;
   }

   char buf[64];

   if (!fid) { 
      mode_t m=umask(0);

      fid=open(fname.data, O_RDWR | O_CREAT | O_APPEND, 0664);
      umask(m); 

      if (fid<=0) { i=fid; fid=0;
         wblog(FL,"ERR %s() failed to open file (e=%d)\n'%s'",
         FCT,i,Wb::repHome(fname).data);
      }
      else if (fid<=2) {
         wblog(FL,"ERR %s() got fid=%d",FCT,fid);
      }
   }
   else if (fid>0) { wblog(FL,"WRN %s() "
     "got existing lock (#%X/%c)\n%s", FCT,fid,flag_,fname.data);
      return fid;
   }
   else  {
      wblog(FL,"WRN %s() got fid=%d/%c for %s", FCT,fid,flag_,fname.data);
      fid=-fid; 
   }

#ifdef QS_USING_OMP
   if (wait) { fLK.acquire(FL,*this); } 
#endif

   i=fcntl(fid, F_SETLK,&fl); 
   if (i<0) { if (wait) {
      gStore.rclog(q, PF_L, CG_VERBOSE>5 && F,
        "RCL waiting for flock (#%X/%c, i=%d; %s)\n%s",
         fid,flag_,i, Wb::hostid('p').data, Wb::repHome(fname).data
      );

      if ((i=fcntl(fid, F_SETLKW,&fl))<0) { 
         close(fid); fid=0; wblog(F_L,
         "ERR %s() failed to obtain flock (#%X/%c, i=%d)",FCT,fid,flag_,i);
      }
   }}
   else { wait=-wait; } 

   snprintf(buf,64,"(#%2d/%c, %d/%d/%.3f; %s)",
      fid, flag_, i,wait, fmod(Wb::getTimeNow(), 60),
      Wb::ompID2Str('x').data); 

   if (i<0) { 
      gStore.rclog(q, PF_L, CG_VERBOSE>5 && F,
      "WRN failed to obtain RCL lock %s\n%s",buf,Wb::repHome(fname).data);
      return (fid=-fid); 
   }

   gStore.rclog(q, PF_L, CG_VERBOSE>6 && F,
     "RCL flock acquired %s\n%s",buf,Wb::repHome(fname).data);

#ifdef DBG_RCS_FLOCKS
   if (flag=='w') { 
      static Wb::counter ctr(FL,"%d entries written to lock files");

      struct timespec tnow;
      unsigned l, i1=16, i2=32, n=64, nout=128;
      char s[n], sout[nout];

      gethostname(s,n); ++ctr;

      clock_gettime(CLOCK_REALTIME,&tnow); 
      strftime(s+i2,n-i2,"%a %D %T",localtime(&tnow.tv_sec));
      snprintf(s+i1,16," %.3f",1e-9*tnow.tv_nsec);
      if (s[i1+1]=='0' && s[i1+2]=='.') { i1+=2; }

      l=snprintf(sout,nout,"%s%-5s %-10s│ %-16s %8d │ %s\n",
         s+i2, s+i1, s,  
         Wb::ompID2Str('v').data, fl.l_pid,
         Wb::repHome(fname).data
      );
      i=write(fid,sout,l);
   }
#endif

   return fid;
};

#endif 
