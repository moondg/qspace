/* ---------------------------------------------------------------------
 * Project : QSpace tensor library (v4.0 pre-release)
 * Class   : QSpace
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

#ifndef __WB_QSPACE_AUX_CC__
#define __WB_QSPACE_AUX_CC__

//====================================================================//
// auxilliary map for QSpace => CGC contraction
//====================================================================//

template <class TQ, class TA, class TB, class TC>
double contractDATA_group(const char *F, int L,
   const QSpace<TQ,TA> &A, C_UVEC &Ia, const ctrIdx &ica,
   const QSpace<TQ,TB> &B, C_UVEC &Ib, const ctrIdx &icb,
   QSpace<TQ,TC> &C, unsigned ic,
   char preview 
){

   unsigned nCGR=0, Nx=0, nx=0, nsym=A.qtype.len, cgflag=A.qtype.maxRank();
   unsigned ra=A.rank(F_L), rb=B.rank(F_L), rc=ra+rb-2*ica.len;
   unsigned rmax=MAX3(ra,rb,rc);
   double x, rval=0, flops=0; 

   if (A.qtype!=B.qtype || A.qtype!=C.qtype) wblog(FL,
      "ERR %s() qtype mismatch '%s' / '%s' / '%s'",
      FCT,STR(A.qtype),STR(B.qtype),STR(C.qtype));
   if (ica.len!=icb.len || ica.len>ra || icb.len>rb) wblog(F_L,
      "ERR %s() invalid contraction indices (%s /%d; %s /%d)",
      FCT,STR(ica),ra,STR(icb),rb);

   if (!Ia.len || Ia.len!=Ib.len) wblog(F_L,
      "ERR %s() invalid index space (n=%d/%d)", FCT, Ia.len, Ib.len);
   if (rc!=C.rank(F_L)) wblog(FL,
      "ERR %s() inconsistent rank %d/%d for C",FCT,rc,C.rank(F_L));
   if (ic>=C.DATA.len) wblog(FL,
      "ERR %s() index out of bounds (ic=%d/%d)",FCT,ic+1,C.DATA.len);

   if (bool(A.CGR) ^ bool(B.CGR)) wblog(FL,
      "ERR %s() got non-intialized CGR data (%s <> %s)",
      FCT, SSTR(A.CGR), SSTR(B.CGR));
   if (!cgflag && (A.CGR || B.CGR || C.CGR)) wblog(FL,
      "ERR %s() got CGR data for all-abelian %s (%p,%p,%p)",
      FCT,STR(A.qtype), A.CGR.data, B.CGR.data, C.CGR.data);
   if (cgflag && (C.CGR.dim1!=C.DATA.len || C.CGR.dim2!=nsym))
      wblog(F_L,"ERR %s() invalid CGR data (C: %dx%d <> %dx%d)",
      FCT,C.CGR.dim1, C.CGR.dim2, C.DATA.len, nsym);

#ifndef WB_SKIP_ASSERT
   A.checkQ_CGR(FL,&Ia); 
   B.checkQ_CGR(FL,&Ib); 
   wbvector< QSet<TQ> > QC(nsym);
   QSet<TQ> Q;
#endif

   if (preview==1) { if (!cgflag) { return rval; }} 
   else if (preview<0 || preview>15) {
      wblog(FL,"ERR %s() invalid preview=%d",FCT);
   }

   unsigned i_,j,l,m=0,M=0, ia,ib; int q=0, xflag=0;
   double cfac=1;
   size_t sA, sB;

   wbvector<unsigned> iOM(nsym), iom(nsym); 
   wbarray<unsigned> Ma(nsym,2), Mb(nsym,2);

   CRef<TQ> Cr_;

   wbperm P, pfin; 
   tensorRef_ R;

   wbarray<TC> Ci, &Ck(*C.DATA.el(ic));

   Ck.init(); 

   if (cgflag) {
      for (M=j=0; j<nsym; ++j) { 
         C.CGR(ic,j).init();
         if (A.qtype[j].permitsOM(rmax)) { iOM[j]=(++M); }
      }            
   }

   #ifdef DBG_CONTRACT
      PRINTF("%s\n",STRREP("─",80));
      wblog(FL,"TST %s() ic=%d / %d",FCT,ic+1,Ia.len);
   #endif

   for (i_=0; i_<Ia.len; ++i_) { ia=Ia[i_]; ib=Ib[i_];
      wbarray<TA> Ai(*A.DATA.el(ia),'r'); 
      wbarray<TB> Bi(*B.DATA.el(ib),'r'); 

      if (cgflag) { iom.set(0);
         for (cfac=1., nx=m=0, j=0; j<nsym; ++j) {
            const CRef<TQ> &Ar=A.CGR(ia,j), &Br=B.CGR(ib,j);
            CRef<TQ> &Cr=C.CGR(ic,j); 

            if (Cr .cgb) { Cr_=Cr; } else
            if (Cr_.cgb) { Cr_.init(); }

            #ifndef WB_SKIP_ASSERT
               Q.init(F_L,Ar,ica,Br,icb);
               if (i_) { if (Q!=QC[j]) { wblog(FL,
                  "ERR %s() mixed QC sector (%d/%d)\n   %s\n<> %s",
                  FCT,i_+1,Ia.len,STR(Q),STR(QC[j]));
               }}
               else { QC[j]=Q; }
            #endif

            nCGR+=gXS.contractCGR(FL,Ar,ica,Br,icb, Cr, preview==1);

            if (preview==1) { continue; } 

            if (!Cr.cgb || !Cr.cgw) { q=0; 
               if (Cr.cgw) {
                  if (Cr.rtype==CGR_ABELIAN) {
                     if (Cr.wscalar1()) { Cr.cgw.init(); continue; }
                     else { q|=1; }  
                  }
                  else { size_t n=Cr.wnumel();
                     if (Cr.rtype==CGR_CTR_SCALAR) { if (!n) { q|=2; }} else
                     if (Cr.rtype==CGR_CTR_ZERO  ) { if ( n) { q|=4; }}
                     else { q|=8; }
                  }
                  if (q) {
                     wblog(FL,"ERR %s() invalid %sCRef\n%s",
                     FCT, q&1? "abelian ": (q&2 ? "scalar ":""), STR(Cr));
                  }
               }
               else if (Cr.isAbelian()) { 
                  continue; 
               }
               else if (Cr.cgb) { wblog(FL,
                  "ERR %s() invalid CRef (missing cgw)\n%s",FCT,STR(Cr)); }
               else if (Cr.rtype!=CGR_CTR_ZERO) wblog(FL,
                  "ERR %s() invalid CRef %s",FCT,STR(Cr));
            }

            if (Cr.cgw.norm()<1e-8) { ++nx; } 

            if (int(l=Cr.isw3())>1) {
               Ma(m,0)=Ar.wdim2();     
               Ma(m,1)=Cr.cgw.SIZE[0]; 

               Mb(m,0)=Br.wdim2();
               Mb(m,1)=Cr.cgw.SIZE[1]; ++m;

               if (!(iom[j]=iOM[j])) wblog(FL,
                  "ERR %s() missing iOM (l=%d, m=%d/%d)",FCT,l,m,M);
            }
            else if (Cr.rtype==CGR_CTR_ZERO) {
               cfac=0; 
               if (!Cr.cgw && Cr_.cgb && Cr_.cgw) { Cr_.save2(Cr); }
               break;                              
            }
            else if (l==1 || Cr.wscalar()) {
               cfac*=Cr.cgw[0];
            }
            else wblog(FL,"ERR %s() got empty x3 (%s)",FCT,SSTR(Cr.cgw));

           #ifdef DBG_CONTRACT
            PRINTF("%s [ic=%d] i=%d/%ld\n"
            "   %-60s@%4s\n   %-60s@%4s\n-> %-60s%5s  %5s\n\n",SHORT_FL,
            ic+1,i_+1,Ia.len, STR(Ar),STR(ica), STR(Br),STR(icb),
            STR2(Cr,'l'), RATS(cfac), RATS(Cr.cgw.norm()));
           #endif
         }
      }
      if (preview==1) { continue; } 

      if (fabs(cfac)<1e-8) { ++Nx;
         if (fabs(cfac)<CG_SKIP_EPS2) { continue; }
         wblog(FL,"WRN %s() got small cfac=%g (nx=%d/%d)",FCT,cfac,nx,nsym);
      }

      if (!m) {
         if (Ai.rank()!=ra || Bi.rank()!=rb) wblog(FL,
            "ERR %s() unexpected DATA rank %s\n%s /%d ; %s /%d",
            FCT, cgflag? "given no OM":"(abelian)",
            SSTR(Ai), ra, SSTR(Bi), rb);

         if (preview) { size_t k,M,K,N; 
            x=sizeof(double); x = sizeof(TA) * sizeof(TB) / (x*x);
            Ai.getMatSize(ica,K,M);
            Bi.getMatSize(icb,k,N); if (k!=K) wblog(FL,
               "ERR %s() K=%ld/%ld (flops=%g)",FCT,K,k,flops);
            flops += x*M*K*N; 
            continue;
         }

         if (!cgflag) {
            if (cfac!=1) wblog(FL,
               "WRN %s() got cfac=%g for all-abelian",FCT,cfac);
            Ai.contract(FL,ica,Bi,icb,Ck,wbperm(),cfac); 
         }
         else { Ci.init();
            Ai.contract(FL,ica,Bi,icb,Ci);

            if (Ck) { Ck.ExpandOM(FL,Ci,rc); }

            Ck.Plus(Ci,TC(cfac),'i',1); 
         }
         continue; 
      }

      if (!cgflag) wblog(FL,  
         "ERR %s() unexpected cgflag=%d",FCT,cgflag);

      tensorRefs TR(2+m); 

      Ai.ExpandOM(FL,ra, m, Ma.data, Ma.data+nsym); 
      Bi.ExpandOM(FL,rb, m, Mb.data, Mb.data+nsym);

      TR[0].init(Ai, &A.itags);
      TR[1].init(Bi, &B.itags);

      TR[0].UpdateItagsCtr(FL,ica,TR[1],icb, QS_ITAG_ "x");  

      TR[0].it.Set(FL, QS_ITAG_ "A",ra, QS_ITAG_ "a",m,'*'); 
      TR[1].it.Set(FL, QS_ITAG_ "B",rb, QS_ITAG_ "b",m,'*'); 

      TR[0].it.FlagItagsCtr(FL,ica,TR[1].it,icb,'!');

      for (l=1, j=0; j<nsym; ++j) { if (iom[j]) { ++l; 
         TR[l].init(C.CGR(ic,j).cgw); 
         TR[l].it.init_alpha(2,1, QS_ITAG_,'a',l-1);
      }}

      if ((sA=TR[0].numel())<1024 && (sB=TR[1].numel())<1024)
           { P.init2End(sA>sB ? 0 : 1,TR.len); }
      else { P.init(); }

      flops+=TR.getOptimalCtrOrder(FL,P,&pfin,'!'); 

      if (preview) {
         continue;
      }

      Ci.init();

      TR[P[0]].contract(FL,TR[P[1]],R,Ci);
      for (j=2; j<P.len; ++j) { R.contract(FL,TR[P[j]],R,Ci); }

      if (Ci.rank()!=rc+m) {
         if (rc || Ci.rank()!=2) wblog(FL,
         "ERR %s() rank mismatch %d+%d -> %d",FCT,rc,m,Ci.rank());
      }

      if (m<M) {
         unsigned i=0, l=0;
         wbvector<wperm_t> px; wbperm p(M);

         if (!pfin) { pfin.init(rc+m); } else
         if (pfin.len!=rc+m) wblog(FL,"ERR %s() "
            "length mismatch (perm.len=%d / %d+%d)",FCT,pfin.len,rc,m);

         pfin.Extend(rc+M);
         px.init(M,pfin.data+rc);

         for (; i<nsym; ++i) { if (iom[i]) { p[l++]=iOM[i]-1; }}
         p.Complete(l); 
         for (i=0; i<M; ++i) { pfin[rc+p[i]]=px[i]; }
      }

      Ci.Permute(pfin); 
      if (Ci.SIZE.len!=pfin.len) { 
         if (Ci.SIZE.len<pfin.len && !pfin) {
            Ci.appendSingletons(pfin.len);
         }
      }

      q=Ck.ExpandOM(FL,Ci,rc);
      xflag = (q<0 ? 0 : q&8);

     #ifdef DBG_CONTRACT
      static int iout=0; sprintf(str,"I%02d",++iout);
      MXPut Ix(FL,str,"base");
      Ix.add(Ck,"Ck_").add(cfac,"cfac").add(Ci,"Ci").add(ic,"ic").add(C,"C");
     #endif

      Ck.Plus(Ci,TC(cfac),'i',1); 

     #ifdef DBG_CONTRACT
      wbarray<TC> Ck_(Ck); Ix.add(Ck_,"Ck");
      if (cgflag) { Ix.add(C.CGR(ic,0),"Cr").add(*C.CGR(ic,0).cgb,"Cb"); }

      PRINTF("%-16s│ %-10s@%4s │ %-10s@%4s │ %d/%ld %d: %-10s %5s  %5s\n",
         SHORT_FL, SSTR(Ai),STR(ica),SSTR(Bi),STR(icb), i_+1, Ia.len, ic+1,
         SSTR(Ci), RATS(cfac), RATS(double(Ck.norm()))
      );
     #endif
   }

   if (cgflag && !preview) {
      unsigned d=-1; char l=Ck.isEmpty(); 
      if (l && (!Nx || !nx)) wblog(FL,"ERR %s() "
         "got empty C.DATA[%d]\nfor valid CGTs (%d/%d)",FCT,ic+1,nx,Nx);

      for (m=1, j=0; j<nsym; ++j) {
         if (xflag) { 
            if (iOM[j])
                 { d=Ck.SIZE.el(rc+iOM[j]-1); }
            else { d=1; }
         }
         m*=C.CGR(ic,j).Reduce_w3Id(d,l);
      }
      if (m) { Ck.FuseOM(FL,rc,m); } else
      if (!l) { wblog(FL,
         "ERR %s() got m=%d\nfor %s DATA w/norm %.3g (d=%d, l=%d)",
         FCT,m,SSTR(Ck),Ck.norm(),d,l);
      }
   }

   if (preview) {
      if (preview&1 && rval) wblog(FL,
         "TST CG_PREVIEW resulted in %d/%d new contractions",rval,Ia.len);
      rval = (preview>1 ? flops : nCGR);
   }
   else { rval=Nx; }

   return rval;
};

template <class TQ, class TA, class TB, class TC>
void contractDATA_plain(const char *F, int L,
   const QSpace<TQ,TA> &A, C_UVEC &Ia, const ctrIdx &ica,
   const QSpace<TQ,TB> &B, C_UVEC &Ib, const ctrIdx &icb,
   QSpace<TQ,TC> &C, unsigned ic
){

   wbvector<unsigned> S(A.qtype.len);
   unsigned i,j, ra=A.rank(F_L), rb=B.rank(F_L), rc=ra+rb-2*ica.len;

   if (A.qtype!=B.qtype) wblog(FL,
      "ERR %s() QSpace symmetry mismatch '%s' <> '%s'",myname,QSTR(A),QSTR(B));
   if (!A.allAbelian()) wblog(FL,
      "ERR %s() got non-abelian '%s' <> '%s'",myname,QSTR(A),QSTR(B));

   if (Ia.len!=Ib.len || !Ia.len) wblog(F_L,"ERR %s()\n"
      "invalid input data (%d/%d)",myname, Ia.len, Ib.len);
   if (C.QIDX.dim1!=C.DATA.len) wblog(F_L,
      "ERR %s() QSpace not yet setup (C: %dx%d <> %dx%d) !?",
      FCT,C.CGR.dim1, C.CGR.dim2, C.DATA.len,C.qtype.len);
   if (rc!=C.rank(F_L)) wblog(FL,
      "ERR %s() inconsistent rc=%d/%d !?",FCT,rc,C.rank(F_L));

   if (ic>=C.DATA.len) wblog(FL,
      "ERR %s() C.DATA out of bounds (%d/%d)",FCT,ic,C.DATA.len);
   if (!C.DATA[ic]) wblog(FL,
      "ERR %s() got null DATA (ic=%d/%d) !?", FCT,ic,C.DATA.len);

   if (C.CGR) {
      wblog(FL,"WRN %s() "
         "got non-empty CGR for qtype=%s",FCT,STR(C.qtype));
      if (C.CGR.dim1!=C.DATA.len || C.CGR.dim2!=C.qtype.len) wblog(F_L,
         "ERR %s() invalid CGR data (C: %dx%d <> %dx%d)",
         FCT,C.CGR.dim1, C.CGR.dim2, C.DATA.len,C.qtype.len
      );
      for (j=0; j<C.CGR.dim2; ++j) {
         if (!C.CGR(ic,j).isAbelian()) wblog(FL,
        "ERR %s() got non-abelian CRef data !?",FCT);
      }
   }

   C.DATA[ic]->init();

   if (C.mt && WBLOG_MMEX) {
      if (C.mt==C.DATA[ic]->mtype) wblog(FL,
         " *  %s ... %2d: %p",Wb::MTYPE_STR[C.mt],ic,C.DATA[ic]);
      else { wblog(FL,"WRN ... %2d: %s != %s !?", ic,
      Wb::MTYPE_STR[C.mt], Wb::MTYPE_STR[C.DATA[ic]->mtype]); }
   }

   try {
      for (i=0; i<Ia.len; ++i) {
         A.DATA.el(Ia[i])->contract(
            FL,ica, *B.DATA.el(Ib[i]),icb, *C.DATA[ic]
         );
      }
   }
   catch (...) {
    { unsigned ia=Ia[i], ib=Ib[i];
      const wbarray<TA> &a=(*A.DATA[ia]);
      const wbarray<TB> &b=(*B.DATA[ib]);
      const wbarray<TC> &c=(*C.DATA[ic]);

      wblog(F,L,"ERR %s::%s\n"
         "A(%d): %6s  [ %s ] @ [%s]\n"
         "B(%d): %6s  [ %s ] @ [%s]\n"
         "C(%d): %6s  [ %s ]", SHORT_FL, FCT,
      ia+1, SSTR(a), A.qrec2Str(ia).data, STR(ica+1),
      ib+1, SSTR(b), B.qrec2Str(ib).data, STR(icb+1),
      ic+1, SSTR(c), C.qrec2Str(ic).data);
   }}

   if (C.mt && WBLOG_MMEX) {
      if (C.mt!=C.DATA[ic]->mtype) { wblog(FL,
         "WRN ... %2d: %s != %s !?",ic,
         Wb::MTYPE_STR[C.mt], Wb::MTYPE_STR[C.DATA[ic]->mtype]);
      }
   }
};

template<class TQ, class TD>
void getQall(const wbvector< QSpace<TQ,TD> > &F, wbMatrix<TQ> &QA){

   unsigned i,i0,k=0,QDIM,dim2;

   wbvector < wbMatrix<TQ>* > pQ(F.len);

   for (i=0; i<F.len; ++i) if (F[i]) { break; }
   if (i>=F.len) { QA.init(); return; }

   pQ[k++]=&(F[i].QIDX); QDIM=F[i].QDIM; dim2=F[i].QIDX.dim2;

   if (QDIM==0 || dim2==0 || dim2 % QDIM) wblog(FL,
   "ERR invalid operator (empty; %d/%d)",dim2,QDIM);

   for (i0=i, i=i0+1; i<F.len; ++i) { if (!F[i]) continue;

      if (F[i].qtype!=F[i0].qtype) wblog(FL,
         "ERR qtype inconsistency %d,%d/%d: %s <> %s",i0+1,i+1,F.len,
         F[i].qStr().data, F[i0].qStr().data
      );

      if (F[i].QDIM!=QDIM || F[i].QIDX.dim2!=dim2)
      wblog(FL,"ERR dimensional inconsistency %d/%d: %d/%d <> %d/%d",
      i0+1, i+1, F.len, F[i].QIDX.dim2, F[i].QDIM, dim2, QDIM);

      pQ[k++]=&(F[i].QIDX);
   }

   pQ.len=k;

   QA.CAT(1,pQ);
   QA.Reshape(QA.dim1*QDIM, QA.dim2/QDIM);
   QA.makeUnique();
};

int mxIsQSpace( 
   const char *F, int L,
   const mxArray *A,
   unsigned &rank,       
   char cflag,           
   unsigned k,           
   const unsigned *rmin, 
   const unsigned *rmax, 
   const char *istr
){
   int q, isq=0, ise=0;
   unsigned l,nA,nA_, M=0, N=0, i=(istr ? strlen(istr) : 0);
   unsigned r, r1=(rmin ? *rmin : rank), r2=(rmax ? *rmax : r1);

   mxArray *Q, *D, *a;
   char *sp=str;

   int fidQ,fidD, vflag=(L || (i?1:0));

   if (F) { vflag+=2; }
   if (vflag) { str[0]=0; if ((vflag>1 && i) || int(k)>=0) {
      if (int(k)>=0)
           { sp+=snprintf(str,256,"%s(%d): ", i ? istr:"A", k+1); }
      else if (i>24) { sp+=snprintf(str,256,"%s\n",istr); }
      else { sp+=snprintf(str,256,"%s: ", i ? istr : "A"); }
   }}

   if (!A || !mxIsStruct(A) || !Mx::IsVector(A)) {
      if (vflag) { sprintf(sp,
         "invalid QSpace object (%s)",mxTypeSize2Str(A).data);
      if (vflag>1) wblog(F_L,"ERR %s",str); }
      return -1;
   }

   nA=nA_=mxGetNumberOfElements(A);
   if (int(k)>=int(nA)) {
      if (vflag) { sprintf(sp,"structure index out of bounds "
         "(%s <> %d)",mxTypeSize2Str(A).data,int(k)+1);
      if (vflag>1) wblog(F_L,"ERR %s",str); }
      return -2;
   }

   if (int(k)==-1 && nA>1) { 
      if (vflag) {
         sprintf(sp,"expecting single QSpace (got %d)",nA);
      if (vflag>1) wblog(F_L,"ERR %s",str); }
      return -3;
   }

   fidQ=mxGetFieldNumber(A,"Q");
   fidD=mxGetFieldNumber(A,"data");

   if (fidQ<0 || fidD<0) {
      if (vflag) {
         strcpy(sp,"mxArray is not of type {Q,data,...}");
      if (vflag>1) wblog(F_L,"ERR %s",str); }
      return -4;
   }

   if (int(k)<0) { k=0; } 
   else { nA_=k+1; }      

   for (; k<nA_; ++k) {
      Q=mxGetFieldByNumber(A,k,fidQ);
      D=mxGetFieldByNumber(A,k,fidD);

      if (!Q || !D) {
         if (!Q && !D) continue; 
         if ((Q && !mxIsEmpty(Q)) || (D && !mxIsEmpty(D))) {
            if (vflag) {
               strcpy(sp,"QSpace data only partially set !?");
            if (vflag>1) wblog(F_L,"ERR %s",str); }
            return -5;
         }
      }

      if ((!Q || mxIsEmpty(Q)) && (!D || mxIsEmpty(D))) {
         ise+=1; continue;
      }

      if (Q) { r=mxGetNumberOfElements(Q); } else { r=0; }
      if (r) {
         if (!Mx::IsVector(Q) || !mxIsCell(Q)) {
            if (vflag) { sprintf(sp,
               "Q-field not a vector cell (%s)",mxTypeSize2Str(Q).data);
            if (vflag>1) wblog(F_L,"ERR %s",str); }
            return -6;
         }
      }

      if (int(rank)<0) {
         rank=r; if (int(r1)<0) r1=rank; if (int(r2)<0) r2=r1;
         if (r1>r2) wblog(F_L,
         "ERR %s() invalid rank range [%d %d] !?)",FCT,r1,r2);
      }
      else if (r<r1 || r>r2) {
         if (vflag) { sprintf(sp,
            "rank r=%d not in range [%d %d])",r,r1,r2);
         if (vflag>1) wblog(F_L,"ERR %s",str); }
         return -7;
      }

      for (i=0; i<r; ++i) { a=mxGetCell(Q,i);
         if (Mx::IsNumArray(0,L,a)<=0) { 
            if (vflag) { unsigned l=0, n=512; char s[n];
               l+=snprintf(s,n,"invalid field %s",istr? istr:"A");
               if (nA>1) { l+=snprintf(s+l,n-l,"(%d/%d)",k+1,nA); }
               l+=snprintf(s+l,n-l,".Q{%d} (n=%d)\n%.256s",i+1,
                  a ? int(mxGetNumberOfElements(a)) : -1, str);
               if (vflag>1) wblog(F_L,"ERR %s",s); else strcpy(str,s);
            }
            return -8;
         }
         if (!i) { M=mxGetM(a); N=mxGetN(a); }
         if ((i && (M!=mxGetM(a) || N!=mxGetN(a))) ||
             mxGetNumberOfDimensions(a)>2
          ){
            if (vflag) { sprintf(sp,"dimension mismatch in Q{%d} !?",i+1);
            if (vflag>1) wblog(F_L,"ERR %s",str); }
            return -9;
         }
      }

      if (!D) {
         if (mxGetFieldNumber(A,"data")<0) {
            if (vflag) { strcpy(sp,"missing field 'data'");
            if (vflag>1) wblog(F_L,"ERR %s",str); }
            return -10;
         }
      }
      else {
         l=mxGetNumberOfElements(D);

         if (l!=M && (M || l!=1)) {
            if (vflag) { sprintf(sp,
               "QSpace inconsistency data: %d/%s",M,mxSize2Str(D).data);
            if (vflag>1) wblog(F_L,"ERR %s",str); }
            return -11;
         }
         if (!l) continue;

         if (!Mx::IsVector(D) || !mxIsCell(D)) {
            if (vflag) { sprintf(sp,
               "data field not a cell vector (%s)",mxTypeSize2Str(D).data);
            if (vflag>1) wblog(F_L,"ERR %s",str); }
            return -12;
         }

         for (i=0; i<l; ++i) { a=mxGetCell(D,i);
            if ((q=Mx::IsNumArray(0,L,a,r,cflag))<=0) {
               if (r>2 && unsigned(-q)==r+1) {
                  q=Mx::IsNumArray(0,L,a,r+1,cflag);  
               }
            }
            if (q<=0) {
               if (vflag) { unsigned l=0, n=512; char s[n];
                  l+=snprintf(s,n,"invalid field %s",istr? istr:"A");
                  if (nA>1) { l+=snprintf(s+l,n-l,"(%d/%d)",k+1,nA); }
                  l+=snprintf(s+l,n-l,".data{%d}\n%.256s",i+1,str);
                  if (vflag>1) wblog(F_L,"ERR %s",s); else strcpy(str,s);
               }
               return -13;
            }
            isq|=q; 
         }
      }
   }

   if (ise) { isq|=256; } 
   return isq;
};

bool mxIsEmptyQSpace(const mxArray *a, unsigned k) {

   mxArray *Q=NULL, *D=NULL;
   unsigned i1,i2;
   int fidQ, fidD;

   if (!a || mxIsEmpty(a)) return 1;

   if (k>=(unsigned)mxGetNumberOfElements(a)) wblog(FL,
      "ERR index out of bounds (%d/%d)",k,mxGetNumberOfElements(a));

   fidQ=mxGetFieldNumber(a,"Q");
   fidD=mxGetFieldNumber(a,"data");

   if (fidQ<0 || fidD<0) wblog(FL,
      "ERR mxArray is not of type {Q,data,...}", FL);

   Q=mxGetFieldByNumber(a,k,fidQ);  i1=(Q && !mxIsEmpty(Q));
   D=mxGetFieldByNumber(a,k,fidD);  i2=(D && !mxIsEmpty(D));

   if (i2 && !i1) { 
      if (mxGetNumberOfElements(D)==1) {
         const mxArray *a=mxGetCell(D,0);
         if (mxGetNumberOfElements(a)==1) { return 0; }
      }
   }
   if (i1 && !i2) { 
      unsigned i=0, n=mxGetNumberOfElements(Q);
      for (; i<n; ++i) {
         const mxArray *a=mxGetCell(Q,i);
         if (mxGetNumberOfElements(a)) { break; }
      }
      if (i==n) { return 1; }
   }

   if (i1 || i2) {
      if (i1 ^ i2) wblog(FL,
         "WRN Q or data set in QSpace, but not both !?");
      return 0;
   }

   return 1;
}

bool mxIsScalarQSpace(const mxArray *a, unsigned k) {

   int fidQ, fidD, fidI;
   mxArray *Q=NULL, *D=NULL, *I=NULL;

   if (!a || mxIsEmpty(a)) return 0;

   if (k>=(unsigned)mxGetNumberOfElements(a)) wblog(FL,
      "ERR index out of bounds (%d/%d)",k,mxGetNumberOfElements(a));

   fidQ=mxGetFieldNumber(a,"Q");
   fidD=mxGetFieldNumber(a,"data");
   fidI=mxGetFieldNumber(a,"info");

   if (fidQ<0 || fidD<0 || fidI<0) wblog(FL,
      "ERR mxArray is not of type {Q,data,info}", FL);

   Q=mxGetFieldByNumber(a,k,fidQ);
   if (Q && mxGetNumberOfElements(Q)) return 0;

   D=mxGetFieldByNumber(a,k,fidD);
   if (!mxIsCell(D) || mxGetNumberOfElements(D)!=1) return 0;
   else {
      const mxArray *a=mxGetCell(D,0);
      if (mxGetNumberOfElements(a)!=1) return 0;
   }

   I=mxGetFieldByNumber(a,k,fidI);
   if (I && (mxGetM(I)>1 || mxGetNumberOfDimensions(I)>2)) return 0;

   return 1;
};

bool mxIsQSpaceArr(
    const char *F, int L,
    const mxArray *S, unsigned rank, int arrdim, 
    const unsigned *rmin, 
    const unsigned *rmax, 
    char cflag
){
    unsigned e=0, m=mxGetM(S), n=mxGetN(S);
    if (L) str[0]=0;

    if (!m || !n) {
       if (L) { sprintf(str,"%s() received empty array (%dx%d)",FCT,m,n); }
       if (F) { wblog(F,L,"WRN %s",str); }
       return 1; 
    }

    if (!mxIsStruct(S) || mxGetNumberOfDimensions(S)!=2) {
       if (L) { sprintf(str,
          "invalid structure array (%s)", mxTypeSize2Str(S).data); }
       if (F) { wblog(F,L,"ERR %s",str); }
       return 0;
    }

    switch (arrdim) {
       case 0: if (m!=1 || n!=1) e++; break;
       case 1: if (m!=1 && n!=1) e++; break;
    }
    if (e) {
       if (L) sprintf(str,
          "invalid QSpace %s (%dx%d)",arrdim==0 ? "scalar":"vector",m,n);
       if (F) { wblog(F,L,"ERR %s",str); }
       return 0;
    }

    if (mxIsQSpace(F,L,S,rank,cflag,-2,rmin,rmax)<=0) {
        return 0;
    }

    return 1;
};

bool mxIsQSpaceVecVec(const mxArray *C, unsigned rank, char cflag) {

    unsigned i,m=0,n=0;
    const mxArray *a;

    if (C) { m=mxGetM(C); n=mxGetN(C); }

    if (!C || !m || !n) {
       wblog(FL,"WRN got empty object (%s; %dx%d)\n--> %s",
       mxGetClassName(C),m,n, SHORT_FL);
       return 1; 
    }

    if (!mxIsCell(C) || (m>1 && n>1)) {
       wblog(FL,"ERR got invalid object (%s)\n",mxTypeSize2Str(C).data);
       return 0;
    }

    n*=m;

    for (i=0; i<n; i++) { a=mxGetCell(C,i);
       mxIsQSpace(FL,a,rank,cflag,-2); 
       return 0;
    }

    return 1;
}

bool mxsIsQSpaceVec(const mxArray *S, int fid, unsigned rank, char cflag) {

    unsigned m=mxGetM(S), n=mxGetN(S);
    mxArray *a;

    if ((m!=1 && n!=1) || !mxIsStruct(S) || mxGetNumberOfDimensions(S)!=2) {
       sprintf(str, "%s:%d need vector structure (%dx%d)", FL,m,n);
       return 0;
    }

    n*=m;

    if (!n || fid<0) { sprintf(str,
      "%s:%d invalid or empty data set (%d;%d)", FL, n, fid);
       return 0;
    }

    for (unsigned k=0; k<n; k++) {
       a=mxGetFieldByNumber(S,k,fid); if (!a || mxIsEmpty(a)) continue;
       mxIsQSpace(FL,a,rank,cflag,0);
    }

    return 1;
};

bool mxsIsQSpaceVEC( 
    const mxArray *S, int fid,
    unsigned rank, char cflag
){
    unsigned i,k, d=0, m=mxGetM(S), n=mxGetN(S);
    mxArray *a;

    if ((m!=1 && n!=1) || !mxIsStruct(S) ||
        mxGetNumberOfDimensions(S)!=2) { sprintf(str,
      "%s:%d need vector structure of QSpace vectors \n"
      "containing {Q,data,...} structures (%d)", FL, m); return 0;
    }

    n*=m;

    if (!n || fid<0) { sprintf(str,
      "%s:%d invalid or empty data set (%d;%d)", FL, n, fid);
       return 0;
    }

    for (i=0; i<n; i++) {
       a=mxGetFieldByNumber(S,i,fid); if (!a || mxIsEmpty(a)) continue;

       m=mxGetM(a); k=mxGetN(a);

       if (d==0) d=m*k;

       if ((m!=1 && k!=1) || n!=d || mxGetNumberOfDimensions(a)!=2 ||
          !mxIsCell(a)) { sprintf(str,
          "%s:%d field must be %d-dim QSpace vector structure (%dx%d)",
           FL, d, m, k); return 0;
       }

       for (k=0; k<d; k++) {
           mxIsQSpace(FL,a,rank,cflag,k);
           return 0;
       }
    }

    return 1;
}

int mxIsQSpaceVec(
   const char *F, int L,
   const char *fname0, const char *vname,
   unsigned rank,
   const unsigned *rmin, 
   const unsigned *rmax, 
   char cflag,
   unsigned N 
){
    unsigned i, imax=1000, r1=-1,r2=-1, Nflag=0; int q;
    mxArray *a;

    if (int(N)>0) {
       if (N>imax) wblog(FL,"ERR %s() expected chain length N=%d !?",FCT,N);
       Nflag=2; 
    }

    wbstring fname(strlen(fname0)+12);
    char xstr[128]; xstr[127]=0;

    snprintf(xstr,126,"%s|%s|%s(%s%s)",
       shortFL(F_L), SHORT_FL, FCT,
       shortFL(fname0 ? fname0:"<data>"), (vname ? vname:"")
    );

    for (i=0; i<imax; ++i) { 
       if (Nflag && i>=Nflag) { i=N-Nflag; Nflag=imax; }
       snprintf(fname.data,fname.len,"%s_%02d.mat", fname0, i);

       Wb::matFile fid;
       if (!fid.open(FL,fname.data,"r")) {
          if ((i%10)!=1) {
             printf("\r     %s %2d ... \r",xstr,i-1); } 
          break;
       }

       printf("\r !X! "); 

       a=matGetVariableInfo(fid.mfp,vname);

       if ((i%10)==0) {
          printf("\r     %s %2d ... \r",xstr,i); } 

       if (!a || mxIsEmpty(a)) { continue; }
       if ((q=mxIsQSpace(F,L,a,rank,cflag))<=0) { return -1; }

       if (int(r1)<0 || int(r2)<0) {
          if (int(r1)<0) { r1=(rmin ? *rmin : rank); }
          if (int(r2)<0) { r2=(rmax ? *rmax : r1); }
          if (r1>r2) wblog(FL,
             "ERR %s() invalid rank-range [%d %d]",FCT,r1,r2);
       }

       if (rank<r1 || rank>r2) {
          if (F) wblog(F,L,
             "ERR invalid QSpaceVec (%d: rank=%d [%d..%d])\n%s '%s'",
              i,rank,r1,r2,SHORT_FL,str);
          return -2;
       }
       if (r1!=r2) rank=-1;

       mxDestroyArray(a);
    }
    printf("\r%85s\r",""); 

    if (Nflag && i!=N) wblog(FL,
       "ERR  inconsistent NRG data (chain length %d/%d) !?", i,N);
    if (i>=imax) wblog(FL,"ERR more than %d files !?", imax);

    return i;
};

int mxIsQSpaceVEC( 
   const char *F, int L,
   const char *fname0, const char *vname,
   unsigned rank,
   const unsigned *rmin, 
   const unsigned *rmax, 
   char cflag
){
    unsigned i,k,m,n,d=0,imax=1000, r1=-1,r2=-1;
    mxArray *a;

    wbstring fname(strlen(fname0)+12);

    for (i=0; i<imax; i++) {
       snprintf(fname.data,fname.len,"%s_%02d.mat", fname0, i);

       Wb::matFile fid;
       if (!fid.open(F_L,fname.data,"r")) break;

       a=matGetVariableInfo(fid.mfp,vname);

       if (!a || mxIsEmpty(a)) continue;

       m=mxGetM(a); n=mxGetN(a);

       if (d==0) d=m*n;

       if ((m!=1 && n!=1) || n!=d || !mxIsCell(a) ||
           mxGetNumberOfDimensions(a)!=2) {
           if (F) wblog(F,L,"ERR sequence of %d-dim vectors of QSpace\n"
              "structures required (%dx%d; %d)",d,m,n,i);
           return -1;
       }

       for (k=0; k<d; k++) {
          if (mxIsQSpace(F,L,a,rank,cflag,k)<=0) return -1;
          if (i==0 && k==0) {
             r1=(rmin ? *rmin : rank);
             r2=(rmax ? *rmax : r1);
             if (r1>r2) wblog(F,L,
                "ERR %s() invalid rank-range [%d %d]",FCT,r1,r2);
          }

          if (rank<r1 || rank>r2) {
             if (F) wblog(F,L,
                "ERR invalid QSpaceVec (%d: rank=%d [%d..%d])\n%s",
                 k+1,rank,r1,r2,str);
             return -1;
          }
          if (r1!=r2) rank=-1;
       }

       mxDestroyArray(a);
    }

    if (i>=imax) wblog(F,L,
    "WRN More than %d files - skip !?", imax);

    return i;
}

template<class TQ, class TD>
void mxInitQSpaceVec(
   const char *F, int L, const mxArray* S,
   wbvector< QSpace<TQ,TD> > &FF,
   const char ref,
   const unsigned *rmin, 
   const unsigned *rmax  
){
   unsigned k,m,n,l=0;
   char cflag=ISCOMPLX_(TD); 

   if (!S) { FF.init(); return; } 

   m=mxGetM(S); n=mxGetN(S);
   if (m==0 && n==0) { FF.init(); return; } 

   if (!mxIsQSpaceVec(F_L,S,-1,rmin,rmax,cflag)) {
      wblog(FL,"WRN %s() got invalid QSpaceVec !?",FCT); 
      return;
   }

   if (!mxIsStruct(S) || mxGetNumberOfDimensions(S)!=2)
       wblog(F,L,"ERR %s needs QSpace vector on input\n(%d, %dx%d)",
       FCT, mxIsStruct(S), m, n);

   n*=m;

   if (FF.len!=n) FF.initDef(n); 

   for (k=0; k<n; ++k) {
      FF[k].init(F,L,S,ref,k);
      if (l) FF[k].checkQ(F,L,FF[l-1]); else
      if (FF[k]) { l=k+1; }
   }
};

template<class TQ, class TD>
void mxInitQSpaceVecR23(const char *F, int L,
   const mxArray* S,
   wbvector< QSpace<TQ,TD> > &FF,
   const char ref=0
){
   unsigned rmin=2, rmax=3;
   mxInitQSpaceVec(F,L,S,FF,ref,&rmin,&rmax);
};

template<class TQ, class TD>
void mxcInitQSpaceVec(const char *F, int L,
   const mxArray* C,
   wbvector< QSpace<TQ,TD> > &FF,
   const char ref
){
   unsigned i,m,n;

   if (!C) { FF.init(); return; }

   m=mxGetM(C); n=mxGetN(C);

   if (!mxIsCell(C) || mxGetNumberOfDimensions(C)!=2 || (m!=1 && n!=1))
   wblog(F,L,"ERR %s needs QSpace vector on input",FCT);

   n*=m;

   if (FF.len!=n) FF.initDef(n); 

   for (i=0; i<n; i++)
   FF[i].init(F,L,mxGetCell(C,i),ref);
}

template<class TQ, class TD>
void mxInitQSpaceVecVec(const char *F, int L,
   const mxArray* C,
   wbvector< wbvector< QSpace<TQ,TD> > > &FF,
   const char ref
){
   unsigned i,m=0,n=0;
   const mxArray *a;

   if (C) { m=mxGetM(C); n=mxGetN(C); }

   if (!C || m==0 || n==0) {
   FF.init(); return; } 

   if (!mxIsCell(C) || (m!=1 && n!=1)) {
      wblog(F,L,"ERR invalid QSpace cell vector: got %s(%d,%d)\n"
      "--> %s:%d", mxGetClassName(C),m,n, Wb::basename(__FILE__),__LINE__);
   }

   n*=m; if (FF.len!=n) FF.initDef(n);

   for (i=0; i<n; i++) { a=mxGetCell(C,i);
       if (mxIsQSpaceVec(FL,a))
          mxInitQSpaceVec(F,L,a,FF[i],ref);
       else {
          wbstring istr(str); sprintf(str,
            "invalid QSpace cell vector (%d: %s)",i+1,mxGetClassName(C));
          wblog(F,L,"%s\n%s",str,istr.data);
       }
   }
}

template<class TQ, class TD>
void mxInitQSpaceMat(const char *F, int L,
   const mxArray* S,
   wbMatrix< QSpace<TQ,TD> > &FF,
   const char ref 
){
   unsigned i,j,m,n;

   if (!S) { FF.init(); return; } 

   if (mxIsCell(S)) { mxcInitQSpaceVec(F,L,S,FF,ref); return; }

   m=mxGetM(S); n=mxGetN(S);
   if (m==0 && n==0) { FF.init(); return; } 

   if (!mxIsQSpaceMat(FL,S)) { wbstring istr(str);
      sprintf(str,"invalid QSpace (%s)", mxGetClassName(S));
      wblog(F,L,"%s\n%s",str,istr.data);
   }

   if (!mxIsStruct(S) || mxGetNumberOfDimensions(S)!=2)
       wblog(F,L,"ERR %s needs QSpace array on input\n(%d, %dx%d)",
       FCT, mxIsStruct(S), m, n);

   if (FF.dim1!=m || FF.dim2!=n) FF.initDef(m,n);

   for (i=0; i<m; i++)
   for (j=0; j<n; j++) FF(i,j).init(F,L,S,ref,i+j*m);
}

template<class TQ, class TD>
void mxcInitQSpaceMat(const char *F, int L,
   const mxArray* C,
   wbMatrix< QSpace<TQ,TD> > &FF,
   const char ref
){
   unsigned i,j,m,n;
   const mxArray *a;

   if (!C) { FF.init(); return; }

   m=mxGetM(C); n=mxGetN(C);

   if (!mxIsCell(C) || mxGetNumberOfDimensions(C)!=2 || (m!=1 && n!=1))
   wblog(F,L,"ERR %s needs QSpace cell vector on input.",FCT);

   m*=n;

   for (i=0; i<m; i++) {
      a=mxGetCell(C,i); mxIsQSpaceVec(F,L,a);

      if (i>0) {
         if (n!=mxGetNumberOfElements(a)) { wbstring istr(str);
            sprintf(str,"cell QSpace dimensions mismatch (%d/%d)",
            n,mxGetNumberOfElements(a));
            wblog(F,L,"%s\n%s",str,istr.data);
            wblog(FL,"ERR");
         }
      }
      else {
         n=mxGetNumberOfElements(a);
         if (FF.dim1!=m || FF.dim2!=n) FF.initDef(m,n);
      }

      for (j=0; j<n; j++)
      FF(i,j).init(F,L,a,ref,j);
   }
}

template<class TQ, class TD>
mxArray* QSpaceVecVec2Mx(
   const wbvector< wbvector< QSpace<TQ,TD> > > &F
){
   mxArray *S;

   S=mxCreateCellMatrix(1,F.len);

   for (unsigned i=0; i<F.len; i++)
   mxSetCell(S,i,F[i].toMx());

   return S;
}

template<class TQ, class TD>
mxArray* QSpaceVec2Mx(
   const wbvector< QSpace<TQ,TD> > &F,
   wbindex &I
){
   if (!I.len) { QSpace<TQ,TD> X; return X.toMx(); }

   wbvector<const QSpace<TQ,TD>*> qq(I.len);
   for (unsigned i=0; i<I.len; i++) qq[i]=&(F[I[i]]);

   return qq.toMxP();
}

QVec qsGetSym(const char *F, int L, const mxArray *a, unsigned k) {

   QVec qtype;  
   int fidt, fidI; unsigned l,n;
   mxArray *I=NULL, *t;

   if (!a || mxIsEmpty(a)) return qtype;

   if (mxGetFieldNumber(a,"Q")<0 || mxGetFieldNumber(a,"data")<0)
      wblog(F_L,"ERR mxArray is not of type {Q,data,...}");
   if ((fidI=mxGetFieldNumber(a,"info"))<0) {
      return qtype; 
   }

   n=(unsigned)mxGetNumberOfElements(a);
   if (int(k)<0) { k=0; }
   else {
      if (k>=n) wblog(F_L,"ERR index out of bounds (%d/%d)",k,n);
      n=k+1;
   }

   for (unsigned k1=k; k<n; ++k) {
      I=mxGetFieldByNumber(a,k,fidI);
      if (I && !mxIsEmpty(I)) {
         if ((l=(unsigned)mxGetNumberOfElements(I))!=1) wblog(FL,
            "ERR %s() invalid Q.info structure (len=%d !?)",FCT,l);
         if ((fidt=mxGetFieldNumber(I,"qtype"))<0) wblog(FL,
            "ERR %s() missing Q.info.qtype field !?",FCT);
         t=mxGetFieldByNumber(I,0,fidt);

         if (k==k1) { qtype.init(F_L,t); }
         else {
            QVec q(F_L,t);
            if (q!=qtype) wblog(F_L,"ERR %s() "
               "got QSpace symmetry mismatch (%s/%s; %d/%d) !?",
               FCT,STR(q),STR(qtype),k1+1,k+1
            );
         }
      }
   }

   return qtype;
}

#endif
