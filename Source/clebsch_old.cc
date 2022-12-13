#ifndef __WB_CLEBSCH_OLD_CC__
#define __WB_CLEBSCH_OLD_CC__

// ------------------------------------------------------------------ //
// DEPRECATED parts of clebsch.cc
// but still kept for compatibility reasons but no longer intended
// for actual usage. See clebsch.hh for more information.
// Wb,Sep20,09 // Wb,Dec16,14
// ------------------------------------------------------------------ //

QType load_store_qtype(const char *F, int L, const Wb::matFile &M) {

   QType q0;

   mxArray *a=NULL, *ai=matGetVariable(M.mfp,"info");
   if (!ai) wblog(FL,"ERR LoadStore() missing info structure");
   if (!mxIsStruct(ai) || mxGetNumberOfElements(ai)!=1 ||
       !(a=mxGetField(ai,0,"qtype"))) wblog(FL,
      "ERR %s() missing field info.qtype",FCT);
   q0.init(F_L,a);

   mxDestroyArray(ai);
   return q0;
};

void load_RCStore(const char *F, int L, const QVec &qvec) {

   unsigned l,n=128;
   const char *sval, *stag="CG_STORE";
   char fname[n], reload=0;

   for (unsigned it=0; it<qvec.len; ++it) {
      if (qvec[it].isAbelian()) continue;

      if (!(sval=getenv(stag))) wblog(F_L,
         "ERR undefined environmental variable %s\n"
         "(=> expect CStore in <%s>_%s.mat)",stag,stag,qvec[it].toStr('t').data);

      l=snprintf(fname,n,"%s_%s.mat",sval,qvec[it].toStr('t').data);
      if (l>=n) wblog(FL,
         "ERR %s() env %s string out of bounds (%d/%d)",FCT,stag,l,n);

      if (!Wb::fexist(fname,'f')) wblog(F_L,
         "ERR invalid CStore / non-existing file:\n%s",fname);
      else {
         struct stat fs; stat(fname,&fs);
         time_t &ftime = load_cstore[qvec[it]];

         if (fs.st_mtime <= ftime) continue;

         reload=(ftime ? 1 : 0);
         ftime=fs.st_mtime;
      }

      if (CG_VERBOSE) {
         wblog(PF_L,"--- %46R\ni/o %sloading R+CStore (%s)\n *  %s @ %s", "-",
            reload?"re":"", myname, Wb::repHome(fname).data,
            Wb::size2Str(Wb::getFileSize(fname)).data
         );
      }

      gRS.LoadStore(F_L,fname);
      gCS.LoadStore(F_L,fname);
   }

#if 0
   wblog(FL, "TST %s() @ size=%d",FCT,load_cstore.size());
   for (auto it=load_cstore.begin(); it!=load_cstore.end(); ++it) {
      const QType &q=it->first;
      printf(" %2d> %-4s map3 size: %d\n",
      it->second, q.toStr().data, gCS.map3[q].size());
   }
   wblog(FL,"TST %s()",FCT);
#endif
};

template <class TQ>
unsigned CStore<TQ>::LoadStore(
   const char *F, int L, const char *fname
){
   Wb::matFile M;
   if (!M.open(F_L, fname,"r")) wblog(FL,
      "ERR failed to open file %s",fname ? fname : "(null)");

   mxArray *ai, *aj, *as, *ac=matGetVariable(M.mfp,"CStore");
   QType qt, q0=load_store_qtype(F_L,M);

   if (!ac) wblog(FL,"ERR %s() missing CStore structure",FCT);
   if (!mxIsStruct(ac) || mxGetNumberOfElements(ac)!=1) wblog(FL,
      "ERR %s() invalid CStore structure",FCT);

   unsigned j,n, mtot=0; 
   unsigned isym=0, nsym=mxGetNumberOfFields(ac);
   char istr[48];
   wbstring sq;

   for (; isym<nsym; ++isym) {
      sq.init(mxGetFieldNameByNumber(ac,isym));
      if (qt.init_s(sq.data) || qt!=q0) {
         wblog(FL,"TST %s() skipping field %s",FCT,sq.data);
         continue;
      }
      sq=qt.toStr('t');
      as=mxGetFieldByNumber(ac,0,isym);

      if (!as) wblog(FL,
         "ERR %s() invalid CStore.%s structure",FCT,sq.data);
      if (!mxIsStruct(as) || mxGetNumberOfElements(as)!=1) wblog(FL,
         "ERR %s() invalid CStore.%s structure",FCT,sq.data);

      unsigned ic,mc,Nc, nc=mxGetNumberOfFields(as), istd=-1;

      for (ic=0; ic<=nc; ++ic) {
         if (ic==nc) {
            if (int(istd)<0) wblog(FL,
               "ERR %s() missing field CStore.%s.std",FCT,sq.data);

            const char *fld = mxGetFieldNameByNumber(as,istd);
            ai=mxGetFieldByNumber(as,0,istd);
            n=mxGetNumberOfElements(ai);
            mc=Nc=0;

            int ir3=mxGetFieldNumber(ai,"cgr");
            if (ir3<0) wblog(FL,"ERR invalid CStore.%s.%s "
               "(missing field cgb)", qt.toStr('t').data,fld);

            for (j=0; j<n; ++j) {
               aj=mxGetFieldByNumber(ai,j,ir3);
               Nc+=mxGetNumberOfElements(aj);
               mc+=gCS.add3(qt,aj);
            }

            if (CG_VERBOSE && F && mc) {
               snprintf(istr,31,"CStore.%s.%s",sq.data,fld);
               wblog(PFL," *  %-16s %d/%d CGCs (%d tprods)",istr,mc,Nc,n);
            }

            break;
         }

         const char *fld = mxGetFieldNameByNumber(as,ic);

         ai=mxGetFieldByNumber(as,0,ic);
         if (!ai) wblog(FL,"ERR %s() missing CStore.%s.std",FCT,sq.data);
         if (!mxIsStruct(ai) || mxGetNumberOfElements(ai)<1) wblog(FL,
            "ERR %s() invalid cell array CStore.%s.%s",FCT,sq.data,fld);

         if (!strcmp(fld,"std")) { istd=ic; continue; }
         if (strncmp(fld,"rank",4) || !fld[4]) wblog(FL,
             "ERR %s() invalid field %s",FCT,fld);

         n=mxGetNumberOfElements(ai);
         mc=gCS.add2BUF(qt,ai);

         if (CG_VERBOSE && F && mc) {
            snprintf(istr,31,"CStore.%s.%s",sq.data,fld);
            wblog(PFL," *  %-16s %d/%d CGCs%s",istr,mc,n,
              fld[4]=='3'? " (standard CGCs)":"");
         }

         mtot+=mc;
      }

      if (qt==q0) break;
   }

   if (isym==nsym) wblog(FL,
       "WRN %s() failed to find CStore.%s in store",FCT,sq.data);

   mxDestroyArray(ac);

   return mtot;
};

template <class TQ>
unsigned CStore<TQ>::add2BUF(
   const QType &q, const mxArray* S, unsigned k
){
   if (!S || !mxIsStruct(S)) wblog(FL,"ERR %s() "
      "invalid CStore.%s.rank* structure",FCT,q.toStr('t').data);

   int itype = mxGetFieldNumber(S,"type"),
       iqset = mxGetFieldNumber(S,"qset"),
       iqdir = mxGetFieldNumber(S,"qdir"),
       icgd  = mxGetFieldNumber(S,"cgd" ),
       id    = mxGetFieldNumber(S,"cid" );

   if (itype<0 || iqset<0 || iqdir<0 || icgd<0 || id<0) wblog(FL,
      "ERR invalid CStore.%s.rank* (%d,%d,%d,%d,%d)",
      q.toStr('t').data, itype,iqset,iqdir,icgd,id
   );

   unsigned i=0, n=mxGetNumberOfElements(S), nk=n, l=0; 
   char xflag;
   CDATA_TQ X; QSet<TQ> Q;

   if (int(k)>=0) {
      if (k>=mxGetNumberOfElements(S)) wblog(FL,"ERR CStore.%s.rank* "
         "index out of bounds (%d/%d)", q.toStr('t').data, k+1, n);
      i=k; nk=k+1;
   }

   for (; i<nk; ++i) {
      Q.t   .init(FL,mxGetFieldByNumber(S,i,itype));
      Q.qs  .init(FL,mxGetFieldByNumber(S,i,iqset),0,'!'); 
      Q.qdir.init(FL,mxGetFieldByNumber(S,i,iqdir));

      if (Q.t!=q) wblog(FL,"ERR %s() qtype mismatch (%s <> %s)",
         FCT, Q.t.toStr('t').data, q.toStr('t').data
      );

      CDATA_TQ &C=getBUF(FL,Q);
      if ((xflag=!C.isEmpty())) { C.save2(X); }

      C=Q;
      C.cgd  .init(FL,mxGetFieldByNumber(S,i,icgd),0,q);
      C.cstat.init(FL,mxGetFieldByNumber(S,i,id));

      C.checkNormSign(FL); 

      if (xflag) {
         if (!C.sameAs(X,'l')) {
            MXPut(FL,"a").add(C,"C").add(X,"X");
            wblog(FL,"ERR %s() "
              "got QSet mismatch%N%N   %s  |  %s%N<> %s  |  %s%N", FCT,
               C.toStr().data, C.cstat.toStr('V').data,
               X.toStr().data, X.cstat.toStr('V').data
            );
         }

         if (C.cstat!=X.cstat){
            if (C.olderThan(X)) {
               if ((C.qdir.len>2 && CG_VERBOSE>1) || CG_VERBOSE>2)
               wblog(PFL,"{+} CBUF[%3d/%3d] %s",i+1,BUF.size(),C.toStr().data);
            }
            else {
               if ((C.qdir.len>2 && CG_VERBOSE) || CG_VERBOSE>1)
               wblog(PFL,"{-} CBUF[%3d/%3d] %s",i+1,BUF.size(),C.toStr().data);
               X.save2(C);
            }
         }
         else { 
            if (CG_VERBOSE>1) wblog(PFL,
               " ok CBUF[%3d/%3d] %s",i+1,BUF.size(),C.toStr().data);
         }
      }
      else { ++l; if (CG_VERBOSE>2) wblog(PFL,
        "(+) CBUF[%3d/%3d] %s",i+1,BUF.size(),C.toStr().data); 
      }
   }

   return l;
};

template <class TQ, class TD>
unsigned RStore<TQ,TD>::LoadStore(
   const char *F, int L, const char *fname
){
   Wb::matFile M;
   if (!M.open(F_L, fname,"r")) wblog(FL,
      "ERR failed to open file %s",fname ? fname : "(null)");

   mxArray *a, *as, *ar=matGetVariable(M.mfp,"RStore");
   QType qt, q0=load_store_qtype(F_L,M);

   if (!ar) wblog(FL,"ERR %s() missing RStore structure",FCT);
   if (!mxIsStruct(ar) || mxGetNumberOfElements(ar)!=1) wblog(FL,
      "ERR %s() invalid RStore structure",FCT);

   unsigned i,n, m=0, isym=0, nsym=mxGetNumberOfFields(ar);
   char istr[32];
   wbstring sq;

   for (; isym<nsym; ++isym) {
      qt.init(F_L,mxGetFieldNameByNumber(ar,isym));
      if (qt!=q0) continue;

      as=mxGetFieldByNumber(ar,0,isym);
      sq=qt.toStr('t');

      if (!as) wblog(FL,
         "ERR %s() invalid RStore.%s structure",FCT,sq.data);
      if (!mxIsStruct(as) || mxGetNumberOfElements(as)!=1) wblog(FL,
         "ERR %s() invalid RStore.%s structure",FCT,sq.data);

      as=mxGetField(as,0,"store");
      if (!as) wblog(FL,"ERR %s() missing RStore.%s.store",FCT,sq.data);
      if (!mxIsStruct(as) || mxGetNumberOfElements(as)<1 ||
          !(a=mxGetField(as,0,"Sp"))) wblog(FL,
         "ERR %s() invalid data RStore.%s.store",FCT,sq.data);

      gRS.buf[qt].checkInit(FL,qt); 

      i=0; n=mxGetNumberOfElements(as);
      for (; i<n; ++i) {
         m+=gRS.add(qt,as,i); 
      }

      if (CG_VERBOSE && F && m) {
         snprintf(istr,15,"RStore.%s",sq.data);
         if (m==n)
              snprintf(istr+16,15,"%d ireps",n);
         else snprintf(istr+16,15,"%d/%d ireps",m,n);
         wblog(PFL," *  %-16s %s and rank-2 CGCs",istr,istr+16);
      }

      if (qt==q0) break;
   }

   if (isym==nsym) wblog(FL,
       "WRN %s() failed to find RStore.%s in store",FCT,sq.data);

   mxDestroyArray(ar);

   return m;
};

#endif 
