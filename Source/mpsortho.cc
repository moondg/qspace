#ifndef __WB_ORTHO_CC__
#define __WB_ORTHO_CC__

// ================================================================== //
// QBlock :: class for block transformations of QSpace data{:}
// e.g. combining/blocking data into fused symmetry sectors
// ------------------------------------------------------------------ //
// toBlockMatrix() operates on this->A
// which is input and needs to be set prior to a call to *this
//
// => for all-abelian symmetries, *this can group multiple symmetry
//    sectors into single full matrix MM; typically, however, also
//    in the abelian case, either Q1 or Q2 represents a single
//    total (fused) symmetry sector, which needs to be taken care of
//    in the caller for all-abelian [e.g. using getQsum() on either
//    Q1 or Q2 then yields the same result].
//
// => NB! for non-abelian symmetries, this->A must represent single
//    well-defined symmetry sector to avoid conflicts with fusion
//    rules hence Q1 or Q2 must represent a single leg / index,
//    i.e., a single set of symmetry labels, and thus below either
//    n1 or n2 must equal 1! Then fusion rules are respected
//    automatically (since by construction, all symmetry labels
//    from Q1 must be fuseable into Q2 and vice versa, as already
//    described by the isometries (up to normalization) via
//    CRef -> CData.
//
// Wb,Aug23,08 // Wb,Aug04,15

template <class TQ, class TD, class TE>
void QBlock<TQ,TD,TE>::toBlockMatrix( // ::initBlockMatrix
   wbarray<TD> &MM,   // output matrix
   unsigned K,        // group first K dimension and also the rest
   const wbperm &P    // permutation on data *before* using K
){
   char cgflag=A.gotCGS(FL), wom=0;
   unsigned R=-1, nsym=A.CGR.dim2;
   widx_t i,j,k,d,l,l1=0, r,s, n1,n2, i1,i2, d1,d2, dim1,dim2;

   wbMatrix<unsigned> OM, OMD, mark;
   wbvector<unsigned> M1,M2, omk;
   wbarray<TD> Ak,Mk;
   wbarray<double> CGW;
   wbindex I1,I2;

   groupIndex<widx_t> g1, g2;
   wbvector<widx_t> s1,s2, cD1, cD2;

   if (!A.isConsistent(R)) { wblog(FL,"%s",str); }
   if (A.isEmpty()) {
      if (K || P) wblog(FL,
         "ERR %s() got empty QSpace with K=%d / P.len=%d",FCT,K,P.len);
      MM.init(); this->init();
      return;
   }

   if (int(R)<2) wblog(FL,"ERR %s() got invalid QSpace (rank r=%d)",FCT,R);
   if (int(K)<0) { K=R/2; } else
   if (K>R) wblog(FL,"ERR index out of bounds (k=%d/%d)",K,R);
   if (!K || K==R) wblog(FL,"ERR got empty group (k=%d/%d)",K,R);

   qdim_tot=1; qdim.init2val(nsym,1);

   if (cgflag>0) {
      if (K==   1 ) { wom=2; l1=0; } else 
      if (K==(R-1)) { wom=1; l1=K; }      
      else {       
         wblog(FL,"ERR %s() got non-abelian CGR with K=%d/%d (cg=%d)\n"
         "(hint: fuse tensor to rank-2 since non-abelian)",FCT,cgflag,K,R);
      }

      if (A.CGR) { qdim.init(A.CGR.dim2);
         for (i=0; i<A.CGR.dim1; ++i) {
         for (j=0; j<A.CGR.dim2; ++j) {
            A.CGR(i,j).norm2(); 
            if (i) {
               d = (A.CGR(i,j).cgb ? A.CGR(i,j).qdim(l1) : 1);
               if (d!=qdim[j]) wblog(FL,"ERR %s() "
                  "CGR(%d,%d) inconsisteny qdim=%d/%d",FCT,i+1,j+1,d,qdim[j]
               );
            }
            else if (A.CGR(i,j).cgb) {
               qdim_tot *= (qdim[j] = A.CGR(i,j).qdim(l1));
            }
         }}
      }
      else { qdim.init(); }
   }

   if (K!=1 && K!=(R-1) && (A.itags || A.CGR)) { 
      int e=0; QDir qdir; A.getQDir(qdir);

      if (qdir.len!=R) { e=1; }
      else {
         if (P) { qdir.Permute(P); }
         for (i=1; i<K; ++i) { if (qdir[i]!=qdir[0]) { e=2; break; }}
         for (++i; i<R; ++i) { if (qdir[i]!=qdir[K]) { e=3; break; }}
      }
      if (e) wblog(FL,
      "ERR %s() got mixed qdir=%s (k=%d/%d; e=%d)",FCT,STR(qdir), K,R,e);
   }

   if (P) {
      if (P.len!=R) wblog(FL,
         "ERR invalid permutation (len=%d/%d)",P.len,R);
      I1.init(  K, P.data  );
      I2.init(R-K, P.data+K);
   }
   else {
      I1.Index(0,K-1);
      I2.Index(K,R-1);
   }

   if (!A.QIDX.isUnique()) wblog(FL, 
      "ERR %s() got non-unique Q-labels",FCT);

   A.getQsub(I1,Q1).groupRecs(g1.P,g1.D); g1.setup('f');
   A.getQsub(I2,Q2).groupRecs(g2.P,g2.D); g2.setup('f');

   n1=g1.D.len; S1.init(n1,I1.len); 
   n2=g2.D.len; S2.init(n2,I2.len);

   A.getOM(OMD,&M2,2); 

   if (!A.getOM(OM,&M1,1)) {
      D1.init2val(n1,2, 1); 
      D2.init2val(n2,2, 1);
   }
   else if (n1==1 && !(n2==1 && wom==1)) {
      if (wom!=2) wblog(FL,
         "ERR %s() wom=%d/2 (n1=%d, n2=%d, K=%d/%d)",FCT,wom,n1,n2,K,R);
      D1.init2val(n1,2,         1);
      D2.init2val(n2,2+OM.dim2, 1);
      for (l=i=0; i<n2; ++i, l+=d) { d=g2.D[i];
         if ((D2(i,1)=M1[g2.P[l]])>1) { 
         Wb::cpyRange(D2.ref(i,2), OM.rec(g2.P[l]), OM.dim2); }
      }
   }
   else if (n2==1) {
      if (wom!=1) wblog(FL,
         "ERR %s() wom=%d/1 (n1=%d, n2=%d, K=%d/%d)",FCT,wom,n1,n2,K,R);
      D1.init2val(n1,2+OM.dim2, 1);
      D2.init2val(n2,2,         1);
      for (l=i=0; i<n1; ++i, l+=d) { d=g1.D[i];
         if ((D1(i,1)=M1[g1.P[l]])>1) { 
         Wb::cpyRange(D1.ref(i,2), OM.rec(g1.P[l]), OM.dim2); }
      }
   }
   else { wblog(FL,
      "ERR %s() got composite indices for all groups\n"
      "in the presence of outer-multiplicity (%s)\n"
      "(hint: fuse indices at least for one group since non-abelian)",
      FCT, A.qStr().data);
   }

   mark.init(n1+1,n2+1);   
   idx.init2val(n1,n2,-1); 

   for (i=0; i<A.DATA.len; ++i) {
      A.DATA[i]->getSize(I1,s1); r=g1.Ig[i];
      A.DATA[i]->getSize(I2,s2); s=g2.Ig[i]; 

      idx(r,s)=i; 

      if (++mark(r,s)>1 && !wom) { wblog(FL,
         "ERR %s() block %d not unique\n"
         "(got degenerate records in QSpace w/out OM)",FCT,i+1);
      }
      if (++mark(r,n2)==1) { S1.recSet(r,s1); D1(r,0)=s1.prod(0); }
      else if (S1.recCompare(r,s1)) { wblog(FL,
         "ERR %s() block size mismatch (%d, r: %d/%d)\n[%s, %s]",
         FCT, i, r+1, n1, S1.rec2Str(r,"","x").data, SSTR(s1));
      }
      if (++mark(n1,s)==1) { S2.recSet(s,s2); D2(s,0)=s2.prod(0); }
      else if (S1.recCompare(r,s1)) { wblog(FL,
         "ERR %s() block size mismatch (%d, c: %d/%d)\n[%s, %s]",
         FCT, i, s+1, n2, S1.rec2Str(r,"","x").data, SSTR(s1));
      }
   }

   dim1=D1.cumsum0prodRec(cD1,0,1); 
   dim2=D2.cumsum0prodRec(cD2,0,1); 
   MM.init(dim1,dim2); 

   for (k=0; k<A.DATA.len; ++k) {
      const CRef<TQ> *Rk=(nsym ? A.CGR.rec(k) : NULL);

      r  = g1.Ig[k];  s  = g2.Ig[k];  
      i1 = cD1[r];    i2 = cD2[s];    
      d1 = D1(r,0);   d2 = D2(s,0);   

      if (M1[k]==1) { double cfac=1;
         A.DATA[k]->toMatrixRef(Mk,K,P);

         for (i=0; i<nsym; ++i) {
            if ((l=Rk[i].wnumel())>1) wblog(FL,"ERR %s() "
               "unexpected CGR(%d,%d) %s",FCT,k+1,i+1,SSTR(Rk[i].cgw));
            if (l==1) { cfac*=Rk[i].wget0(); }
         }

         if (Mk.SIZE.len!=2 || Mk.SIZE[0]!=d1 || Mk.SIZE[1]!=d2) wblog(FL,
            "ERR %s() size inconsistency %s / %dx%d",FCT,SSTR(Mk),d1,d2);
         if (mark(r,s)!=1 || M2[k]!=1) { wblog(FL,"ERR %s() " 
            "expecting unique block [mark(%d/%d,%d/%d)=%d] / M2=%d",
            FCT,r+1,n1,s+1,n2,mark(r,s),M2[k]); }
         if (i2+Mk.SIZE[1]>dim2) { wblog(FL,
            "ERR %s() index out of bounds:\n" 
            "setting M(%d:%d,%d:%d) having size %s",FCT,
            i1+1, i1+Mk.SIZE[0], i2+1, i2+Mk.SIZE[1], SSTR(MM)); }
         cfac/=sqrt(double(qdim_tot)); 

         Mk.copyStride( MM.ref(i1,i2), dim1, cfac); 
      }
      else {
         if (mark(r,s)!=1) { wblog(FL,"ERR %s() "       
            "unexpected unique block [mark(%d/%d,%d/%d)=%d] @ M1=%d",
            FCT,r+1,n1,s+1,n2,mark(r,s),M1[k]);
         }
         if (M1[k]<2) { wblog(FL,"ERR %s() got M1[%d]=%d",FCT,k,M1[k]); }
         if (!wom   ) { wblog(FL,"ERR %s() got wom=%d",FCT,wom); }

         for (i=0; i<nsym; ++i) {
            if (i) { CGW.kron(Rk[i].cgw, CGW); }
            else   { CGW=Rk[i].cgw; }
         }

         if (!CGW.isMatrix(M1[k],M2[k])) wblog(FL,"ERR %s() "
            "size mismatch CGW %s / %dx%d",FCT,SSTR(CGW),M1[k],M2[k]);

         if (qdim_tot>1) {
            CGW *= sqrt(1./double(qdim_tot)); 
         }
         else if (qdim_tot<1) wblog(FL,
            "ERR %s() invalid qdim=%d",FCT,qdim_tot);

         A.DATA[k]->contractMat(FL,R+1,CGW,2,Ak); 

         Ak.toMatrixRefM(Mk,K,P,R, wom==1? "132":"");
         {  wbvector<size_t> S(Mk.SIZE);
            if (wom==1)
                 { S[0]*=S[1]; S[1]=S[2]; }
            else { S[1]*=S[2]; }
            S.len=2;
            Mk.Reshape(S);
         }

         if (i1+Mk.SIZE[0]>dim1 || i2+Mk.SIZE[1]>dim2) { wblog(FL,
            "ERR %s() index out of bounds:\n" 
            "setting M(%d:%d,%d:%d) having %s",FCT,
            i1+1,i1+Mk.SIZE[0], i2+1,i2+Mk.SIZE[1], SSTR(MM));
         }

         Mk.copyStride( MM.ref(i1,i2), dim1);  
      }
   }
};

template <class TQ, class TD, class TE>
QBlock<TQ,TD,TE>& QBlock<TQ,TD,TE>::init_bare_refA(
   const QBlock &b, unsigned flag) {

   if (this==&b) {
      if (flag) wblog(FL,"ERR %s() got same object (flag=%d)",FCT,flag);
      return *this;
   }

   if (!strchr("USVX",flag)) wblog(FL, 
      "ERR %s() invalid flag %c<%d>",FCT,flag,flag);

   I0.init(); Ik.init(); It.init(); U.init(); S.init(); Vc.init();

   A.init2ref(b.A);

   if (!flag) { 
      qdim_tot=b.qdim_tot; qdim=b.qdim;
      Q1=b.Q1; S1=b.S1; D1=b.D1;
      Q2=b.Q2; S2=b.S2; D2=b.D2; idx=b.idx;
      return *this;
   }

   int k=0, e=0;

   if (b.S2.dim2==1) { k=2;
      b.Q2.blockSum(b.A.QDIM,Q2);
      if (b.S2.dim2>1 && b.A.isNonAbelian()) { e=1; }
   }
   else { k=1; 
      b.A.getQsum(b.Q1,Q2); 
   }

   if (e || !Q2.recAllEqual() || !Q2.dim1) { 
      MXPut(FL,"q").add(b,"b").add(Q2,"Q2").add(k,"k");
      wblog(FL,
         "ERR %s() got inconsistency in block symmetry\n(%s; [%d+%d])",
         FCT,SSTR(b.idx),b.S1.dim2,b.S2.dim2
      );
   }
   Q2.dim1=1; 

   S2.init(1,1); 
   D2.init2val(1,2, 1);

   qdim_tot=b.qdim_tot;
   qdim=b.qdim;

   if (flag=='U') {
      Q1=b.Q1; S1=b.S1; D1=b.D1; D2[0]=S2[0]=b.U.dim(2);

      if (b.idx.dim2>1) {
         idx.init2val(b.idx.dim1,1, b.A.DATA.len);
      }
      else idx=b.idx;
   }
   else if (flag=='V' || flag=='X') {
      Q1=b.Q2; S1=b.S2; D1=b.D2; D2[0]=S2[0]=b.Vc.dim(2);

      if (b.idx.dim1>1) { 
         idx.init2val(b.idx.dim2,1, b.A.DATA.len); 
      }
      else b.idx.transpose(idx);
   }
   else if (flag=='S') { D2[0]=S2[0]=b.S.len;
      Q1=Q2; D1=D2; S1=S2;
      idx.init(1,1); 
   }
   else wblog(FL,"ERR %s() invalid flag=%d",FCT,flag);

#ifndef WB_SKIP_ASSERT
 { wbvector<widx_t> cD_; widx_t dim2_=S2[0],
      d1=b.D1.cumsum0prodRec(cD_,0,1),
      d2=b.D2.cumsum0prodRec(cD_,0,1);
   if (dim2_>d1 && dim2_>d2) wblog(FL, 
      "ERR %s() unexpected dS=%d (%dx%d)",FCT,dim2_,d1,d2
   ); }
#endif

   return *this;
};

template <class TQ, class TD, class TE>
void QBlock<TQ,TD,TE>::updateBlockDim(unsigned dim2_) {

   if (S2.numel()!=1 || D2.dim1!=1 || D2.dim2<2 || D2[1]!=1)
      wblog(FL,"ERR %s() QBlock not yet initialized\n"
      "S2: %s, D2: %s",FCT,STR(S2),STR(D2)
   );

   D2[0]=S2[0]=dim2_;
};

template <class TQ, class TD, class TE>
void QBlock<TQ,TD,TE>::initFromBlockMatrix(const char *F, int L,
   QSpace<TQ,TD> &B,      
   const wbarray<TD> &MM, 
   const iTags &it,       
   const wbperm *cgp      
) const {

   widx_t i,j,k, i0,i1,i2, d1,d2, r,s, dim1, dim2;
   wbvector<widx_t> cD1, cD2, S; 
   wbvector<TQ> Q;

  #ifndef WB_SKIP_ASSERT
   wbvector<unsigned> Ik(1);
  #endif

   char cgflag=A.gotCGS(F_L);
   char cpR=0; 
   double cfac=1;

   if (cgflag>0) {
      r=it.len; 
      s=S1.dim2+S2.dim2;
      if (s<2 || s!=r) wblog(FL,
         "ERR %s() got rank %d/%d",FCT,s,r);

      if (!cgp) {
         if (A.itags.sameConj(it)) { cpR=1; }
      }
      else if (A.itags.len==cgp->len) {
         iTags tp(A.itags); tp.Permute(*cgp);
         if (tp.sameConj(it)) { cpR=2; }
      }
      if (!cpR) {
         if (it.len!=2 || !it[0].isConj(it[1],0)) wblog(FL,
            "ERR %s() unexpected non-abelian setting\n"
            "r=%d: %s (A: %s); cgp=[%s]", FCT, it.len, STR(it),
            STR(A.itags), STR(cgp ? ((*cgp)) : wbperm())
         );
      }
   }

   if (Q1.dim1!=S1.dim1 || D1.dim1!=S1.dim1 || D1.dim2<2 ||
       Q2.dim1!=S2.dim1 || D2.dim1!=S2.dim1 || D2.dim2<2) { wblog(FL,
      "ERR block index mismatch (%d/%d/%dx%d; %d/%d/%dx%d)",
       Q1.dim1, S1.dim1, D1.dim1, D1.dim2,
       Q2.dim1, S2.dim1, D2.dim1, D2.dim2);
   }

   dim1=D1.cumsum0prodRec(cD1,0,1); 
   dim2=D2.cumsum0prodRec(cD2,0,1);

   if (MM.rank()!=2) wblog(FL,
      "ERR %s() got rank-%d object matrix",FCT,MM.rank());
   if (dim1!=MM.SIZE[0] || dim2!=MM.SIZE[1]) wblog(FL,
      "ERR block size mismatch (%dx%d; %s)",
       dim1, dim2, MM.sizeStr().data);
   if (!A.QDIM || Q1.dim2%A.QDIM || Q2.dim2%A.QDIM) wblog(FL,
      "ERR invalid QDIM (%d,%d / %d)",Q1.dim2,Q2.dim2,A.QDIM);

   B.init(D1.colSum(1)*D2.colSum(1), S1.dim2+S2.dim2, A.QDIM);
   B.qtype=A.qtype;
   B.itags=it;

   if (A.CGR.dim2) {
      B.setupCGR();
   }

   for (k=r=0; r<D1.dim1; ++r) {
   for (  s=0; s<D2.dim1; ++s, ++k) {
      i1 = cD1[r];  i2 = cD2[s];
      d1 = D1(r,0); d2 = D2(s,0); if (!d1 || !d2) {
         wblog(FL,"WRN %s() got empty %dx%d sub-block",FCT,d1,d2);
         continue;
      }

      if ((i0=idx(r,s))<0) {
         double x2=MM.norm2block(i1,i2,d1,d2); 
         if (x2>CG_SKIP_DEPS1 || !isfinite(x2)) { char s[64];
            snprintf(s,64,"block with norm x2=%.3g\nM(%ld:%ld,%ld:%ld) "
               "having %s",x2,i1+1,i1+d1,i2+1,i2+d2,SSTR(MM));
            if (x2>1e-8)
                 { wblog(FL,"ERR %s() got %s",FCT,s); }
            else { wblog(FL,"WRN %s() skipping %s",FCT,s); }
         }
         continue;
      }

      if (cpR && i0>=A.CGR.dim1) { 
         MXPut(FL,"q").add(A,"A").add(B,"B")
           .add(D1,"D1").add(cD1,"cD1").add(D2,"D2").add(cD2,"cD2")
           .add(idx,"idx").add(r+1,"r").add(s+1,"s").add(k+1,"k");
         wblog(FL,"ERR %s() index out of bounds:\nidx(%d,%d) = %d/%d",
         FCT,r+1,s+1,i0+1,A.CGR.dim1);
      }

      S.init(S1.dim2, S1.rec(r), S2.dim2, S2.rec(s));
      Q.init(Q1.dim2, Q1.rec(r), Q2.dim2, Q2.rec(s));

      B.QIDX.recSet(k,Q); cfac=1;
      if (cpR) {
         for (j=0; j<A.CGR.dim2; ++j) {
            B.CGR(k,j).init_wId( A.CGR(i0,j),1,-1,&qdim[j]);

            if (qdim[j]>1) { cfac*=B.CGR(k,j).NormSignW(FL); } 

            if (cgp) { B.CGR(k,j).Permute(*cgp); }
         }

         if (fabs(cfac)<1e-8) wblog(FL,
            "WRN %s() got small cfac=%g",FCT,cfac);

         #ifndef WB_SKIP_ASSERT
            Ik[0]=k;
            B.checkQ_CGR(FL,&Ik); 
         #endif
      }
      else if (cgflag) {
         B.initIdentityCGS(FL);
      }

      if (D1(r,1)==1 && D2(s,1)==1) { 
         if (i1+d1>dim1 || i2+d2>dim2) { wblog(FL,"ERR %s() "
            "index out of bounds (%d:%d; %d:%d) having %s",
            FCT, i1+1,i1+d1, i2+1,i2+d2, SSTR(MM));
         }

         wbarray<TD> Bk; Bk.init_bare(d1,d2);
         Wb::cpyStride(Bk.data, MM.ref(i1,i2), d1,d2,-1,dim1,TD(0),TD(cfac));

         Bk.Reshape(S).save2(*B.DATA[k]);
      }
      else {
         unsigned m=0, wom=0; int e=0;
         widx_t d1_=d1, d2_=d2;

         if (!cpR) wblog(FL,"ERR %s() " 
            "expecting cpR!=0 with OM=%dx%D",FCT,D1(r,1),D2(s,1));

         if (D1(r,1)==1) { wom=2; m=D2(s,1); d2_*=m;
            if (D2.dim2!=2+A.CGR.dim2) { e=1; }
         } else
         if (D2(s,1)==1) { wom=1; m=D1(r,1); d1_*=m;
            if (D1.dim2!=2+A.CGR.dim2) { e=2; }
         }
         else { e=3; }

         if (e) { char sx[80]; snprintf(sx,80,"block (%ld/%ld,%ld/%ld) "
            "at OM=%ldx%ld)\nhaving %s | %s | cpR=%d, '%s' (k=%ld)",
               r+1,D1.dim1, s+1,D2.dim1, D1(r,1), D2(s,1),
               D1.rec2Str(r).data, D2.rec2Str(s).data, cpR,STR(it),k);
            if (e<3)
                 wblog(FL,"ERR %s() invalid D1/D2 data (e=%d)\n%s",FCT,e,sx);
            else wblog(FL,"ERR %s() got multiple OM for\n%s",FCT,sx);
         }

         if (i1+d1_>dim1 || i2+d2_>dim2) { char sx[32];
            if (wom==1)
                 { snprintf(sx,32,"(%d*%ld=%ld) x %ld", m,d1,d1_, d2); }
            else { snprintf(sx,32,"%ld x (%d*%ld=%ld)", d1, m,d2,d2_); }

            wblog(FL,"ERR %s() index out of bounds "
               "(%d:%d/%d; %d:%d/%d)\nhaving block size %s",
               FCT, i1+1,i1+d1_,dim1, i2+1,i2+d2_,dim2, sx
            );
         }

         wbarray<TD> Bk; Bk.init_bare(d1_,d2_);
         Wb::cpyStride(Bk.data, MM.ref(i1,i2), d1_,d2_,-1,dim1,TD(0),TD(cfac));

         if (wom==1)
              { Bk.Reshape(d1,m,d2).Permute("132"); }

         S.Append(m);
         Bk.Reshape(S).save2(*B.DATA[k]);
      }
   }}

   for (k=i=0; i<B.DATA.len; ++i) {
      if (B.DATA[i]->isEmpty()) { continue; }
      if (k<i) {
         B.QIDX.recSet(k,i);
         for (j=0; j<B.CGR.dim2; ++j) { B.CGR(i,j).save2(B.CGR(k,j)); }
         B.DATA[k]=B.DATA[i];
      }; ++k;
   }

   if (k<B.DATA.len) {
      B.QIDX.Resize(k,B.QIDX.dim2);
      B.CGR .Resize(k,B.CGR .dim2);
      B.DATA.Resize(k);
   }
};

template <class TQ, class TD>
void SVD_Data<TQ,TD>::blockSVD(
   const QSpace<TQ,TD> &PSI,
   const unsigned K
){
   unsigned i,n, l=0, r=PSI.rank(FL);
   wbvector<widx_t> D,Dc;
   wbMatrix<TQ> Qb;
   wbindex Ib;
   wbperm P;

   if (PSI.isEmpty()) wblog(FL,"ERR %s() got empty PSI",FCT);
   if (!K || K>=r) wblog(FL,
      "ERR %s() index specs out of bounds (%d/%d)",FCT,K,r);

   if (K!=1) { QDir qdir; PSI.getQDir(qdir);
      if (qdir.len!=r) l=1;
      else {
         for (i=0; i<K; ++i) { if (qdir[i]<=0) { l=2; break; }}
         for (++i; i<r; ++i) { if (qdir[i]!=qdir[K]) { l=3; break; }}
      }
      if (l) wblog(FL,
         "ERR %s() got invalid orthonormalization\n"
         "['%s' => %s grouped into %d+%d, >>]\nhint: %s",
         FCT, STR(PSI.itags), STR(qdir), K,r-K, (l==2 ?
         "only ingoing indices can be fused into outgoing" : (l==3 ?
         "2nd/right group is mixed in+out" : "rank mismatch"))
      );
   }

   if (K==  1) { Ib.Index(0,K-1); } else
   if (K==r-1) { Ib.Index(K,r-1); }
   else {
      if (!PSI.allAbelian()) wblog(FL,
         "ERR %s() got non-abelian symmetries %s\nwith K=%d/%d",
         FCT,PSI.qStr().data,K,r);
      Ib.Index(0,K-1); 
   }

   PSI.getQsum(Ib,Qb).groupRecs(P,D);
   n=D.len; QB.init(n);

   Wb::LogException ex;

   int np=MIN(int(n),QSP_NUM_THREADS);
   if (np<1) np=1;

   D.cumsum_(Dc);

   #pragma omp parallel for num_threads(np) 
   for (i=0; i<n; ++i) { if (!ex) { try {
      QBlock<TQ,TD> &b=QB[i];
      wbarray<TD> MM;

      b.I0.init(D[i],P.data+Dc[i]);
      PSI.getSub(b.I0, b.A);

      b.toBlockMatrix(MM, K); 

      wbSVD(MM, b.U, b.S, b.Vc);  
   }
      catch (Wb::LogException &e_) { ex+=e_; }
      catch (...) { ++ex; }
   }}
   ex.report(FLF);
};

template <class TQ, class TD>
unsigned SVD_Data<TQ,TD>::dmrgTruncate( 
   unsigned Nkmin, unsigned Nkeep, double stol,
   itag_ tx,              
   wbMatrix<double> &SM, 
   QSpace<TQ,TD> &UQ,    
   QSpace<TQ,TD> *SQ,    
   QSpace<TQ,TD> &VC,    
   mxArray *Sout,        
   char vflag
){ 

   unsigned i,j,i0,d, nk=0, nt=0, nq=QB.len, flag=0;
   unsigned QDIM, N;
   int NK=Nkeep;

   int revert_new_idx=0; 

   iTags tu,ts,tv;
   wbperm P, cgp2;

   QBlock<TQ,TD,double> x;
   QSpace<TQ,TD> X;

   wbvector< wbvector<double>* > SV;
   wbvector< wbvector<char> > mark;
   wbvector<double> Stot;

   double snorm, snorm2=0, s2t=0, sfac=1,
        xtol2=(stol>0 ? stol*1E-14 : 0);
   char cgflag=0;

   unsigned l=0, n=128; char istr[n]; istr[0]=0;

   if (!nq) wblog(FL,"ERR %s() got empty QB data",FCT);
   else {
      QDIM=QB[0].A.QDIM;
         if (!QDIM) wblog(FL,"ERR QDIM=%d",QDIM);
      cgflag=QB[0].A.gotCGS(FL);
   }

   UQ.init(); if (SQ) { SQ->init(); }
   VC.init();

   for (N=i=0; i<nq; i++) { N+=QB[i].S.len; } 

   SV.init(nq);

   SM.init(N, cgflag? 2:1); 

   if (!N) wblog(FL,"ERR %s() got empty SVD data",FCT);
   else --N; 

   for (i0=i=0; i<nq; ++i, i0+=d) {
      const wbvector<double> &S=QB[i].S; d=S.len;

      SV[i]=&(QB[i].S);
      snorm2+=S.norm2()*QB[i].qdim_tot;

      if (!cgflag) {
         for (j=0; j<d; ++j) { SM[i0+j]=S[j]; }
      }
      else {
         double *sm=SM.rec(i0);
         for (j=0; j<d; ++j, sm+=2) {
            sm[0]=S[j]; 
            sm[1]=QB[i].qdim_tot;
         }
      }
   }

   SM.SortRecs(P,-1); 

   if (int(Nkmin)<0) { Nkmin=0; }

   snorm=SQRT(snorm2);
   markSet(SV,mark,Nkmin,NK,stol,Stot,0.,0.,-1,"desc"); 

   for (i=0; i<nq; ++i) {
      QBlock<TQ,TD> &b=QB[i];
      s2t+=b.get_s2t(mark[i]); 
      nk+=b.Ik.len; nt+=(mark[i].len-b.Ik.len);
   }

   if (nk!=unsigned(NK) || nk+nt!=SM.dim1) wblog(FL, 
      "WRN size inconsistency (nk=%d/%d, N=%d+%d / %d)",
      nk,NK, nk,nt,SM.dim1);

   if (SM(N,0)>stol) { if (SM.dim1>Nkeep) {
       if (stol) {
          l+=snprintf(istr,n,"Smin=%.3g > stol=%.3g (NK=%d/%d/%li)",
          SM(N,0), stol, Nkmin,NK,SM.dim1); flag++;
       }
       else {
          l+=snprintf(istr,n,"using Nkeep=%d (%d/%d/%li; "
          "stol=%g; Smin=%.3g)",NK,Nkmin,Nkeep,SM.dim1, stol, SM(N,0));
       }}
       else {
          l+=snprintf(istr,n,"keeping all (stol=%g; Smin=%.3g, "
          "NK=%d (%d/%d/%li)",stol, SM(N,0), NK,Nkmin,Nkeep,SM.dim1);
          xtol2=0;
       }
   }
   else {
      if (unsigned(NK)<Nkeep) { 
         if (Nkmin) {
            l+=snprintf(istr,n,"truncation by stol=%.3g or Nkmin=%d "
            "(Nkeep %d->%d; Smin=%.3g < %.3g)", stol,
            Nkmin,Nkeep,NK, SM(N,0), SM(nk? nk-1:0,0));
         }
         else {
            l+=snprintf(istr,n,"truncation by stol=%.3g "
            "(Nkeep %d->%d; Smin=%.3g < %.3g)", stol,
            Nkeep,NK, SM(N,0), SM(nk? nk-1:0,0));
         }
      }
      else if (unsigned(NK)<SM.dim1) {
         l+=snprintf(istr,n,"truncation by Nkeep ([%d,%d] -> %d/%li; "
            "Smin=%.3g < %.3g <= stol=%.3g)", Nkmin,Nkeep,NK,SM.dim1,
            SM(N,0), SM(MAX(0,int(NK)-1),0), stol);
         flag++;
      }
   }

   if (s2t) {
      if (l>=n) wblog(FL,"ERR %s() "
         "string out of bounds (%d/%d)%N%N%s%N",FCT,l,n,istr);
      else if (snorm2>s2t) {
         sfac=snorm/sqrt(snorm2-s2t);
         if (fabs(sfac-1)>0.01) { ++flag;
            snprintf(istr+l,n-l,
               "%sWRN readjusting norm by factor %.4g",l?"\n":"",sfac);
            if (fabs(sfac-1)>0.1 && (vflag || Wb::envVRB&8)) {
               wblog(FL,istr+(l?l+1:0)); 
            }
         }
         for (i=0; i<nq; ++i) QB[i].S*=sfac;
      }
      else if (snorm2==s2t) { snprintf(istr+l,n-l,
         "all states discarded (|SVD|=%.3g / %.3g; Nkeep=%d/%d/%li)",
         snorm,stol,NK,Nkeep,SM.dim1); 
      }
      else wblog(FL,"ERR %s() %g/%g",FCT,snorm2,s2t);
   }

   if (Sout) {
      mxAddField2Scalar(FL, Sout, "svd",    SM.toMx());
      mxAddField2Scalar(FL, Sout, "stol",   numtoMx(stol));
      mxAddField2Scalar(FL, Sout, "svd2tr", numtoMx(s2t)); 

      mxAddField2Scalar(FL, Sout, "sfac",   numtoMx(sfac));
      mxAddField2Scalar(FL, Sout, "Nkmin",  numtoMx(Nkmin));
      mxAddField2Scalar(FL, Sout, "Nkeep",  numtoMx(nk));
      mxAddField2Scalar(FL, Sout, "Ntot",   numtoMx(SM.dim1));
      mxAddField2Scalar(FL, Sout, "flag",   numtoMx(flag));
      mxAddField2Scalar(FL, Sout, "info",   mxCreateString(istr));
   }
   else if (vflag) {
      wblog(FL,"<i> sum(SVD) = %g = %g%+.4g",
         snorm2, round(snorm2), snorm2-round(snorm2));
      wblog(FL,"TST keep %d/%d states (%d)",nk,SM.dim1,Nkeep); if (nt) {
      wblog(FL,":x: truncate %d/%d states", nt,SM.dim1); }
   }

   tx.deConj();

   {  QBlock<TQ,TD> &b=QB[0];

      itag_ t(tx);

      b.A.checkQ(FL); 

      cgp2.init2End(0,b.S2.dim2+1);

      if (b.A.itags.len!=b.S1.dim2+b.S2.dim2) wblog(FL,
         "ERR %s() got invalid/empty itags (%s; s=%d+%d)",
         FCT, b.A.itags.toStr().data, b.S1.dim2+b.S2.dim2);

      if (b.S2.dim2==1) {
         t.Set_(b.A.itags.last());
         if (!t.isConj() && NK) {
            if (!b.A.isOp(2,'L')) { revert_new_idx=1; }
         }
      }
      else if (b.S1.dim2==1) {
         t.Set_(b.A.itags[0]); t.Conj();
      }
      else {
         if (!b.A.itags[0].isConj()) { t.Conj(); }
      }

      tu.init(b.S1.dim2, b.A.itags.data, 1, &t); 

      ts.init(1,&t,1,&t);
      ts[0].Conj(); 

      tv.init(b.S2.dim2, b.A.itags.data+b.S1.dim2, 1, &t);
      tv.last().Conj();  

   }

   for (i=0; i<nq; ++i) {
      QBlock<TQ,TD> &b=QB[i]; if (b.Ik.isEmpty()) continue;
      b.U .Select0(b.Ik,1); 
      b.S .Select(b.Ik);
      b.Vc.Select0(b.Ik,1); 

      if (b.U .rank()!=2) wblog(FL,"ERR U: invalid rank %s",SSTR(b.U));
      if (!b.S.len      ) wblog(FL,"ERR got empty S");
      if (b.Vc.rank()!=2) wblog(FL,"ERR V: invalid rank %s",SSTR(b.Vc));

      x.init_bare_refA(b,'U');     
      x.initFromBlockMatrix(FL,X,b.U,tu);

      if (X.Append2AndDestroy(FL,UQ,'u')) wblog(FL,
         "ERR %s() U: QIDX must not overlap\n%s",FCT,istr);

      if (SQ) {
         wbarray<double> Sf; 
         Sf.initDiag(b.S.len,b.S.data); 

         x.init_bare_refA(b,'S');    
         x.initFromBlockMatrix(FL,X,Sf,ts);

         if (X.DATA.len!=1) wblog(FL,
            "ERR %s() got S.DATA.len=%d",FCT,X.DATA.len);
         X.DATA[0]->Reduce2Diag();

         if (X.Append2AndDestroy(FL,*SQ,'u')) wblog(FL,
            "ERR %s() S: QIDX must not overlap\n%s",FCT,istr);

         x.init_bare_refA(b,'V');     
         x.initFromBlockMatrix(FL,X,b.Vc,tv,&cgp2);

         if (X.Append2AndDestroy(FL,VC,'u')) wblog(FL,
            "ERR %s() V: QIDX must not overlap\n%s",FCT,istr);
      }
      else {
         wbarray<TD> VS; Wb::DMatProd(b.Vc,b.S,VS); 

         if (VS.rank()!=2) wblog(FL,
            "ERR VS: invalid rank-%d",VS.rank());

         x.init_bare_refA(b,'X'); 
         x.initFromBlockMatrix(FL,X,VS,tv,&cgp2);

         X.SkipZeroData(xtol2,'b',0);

         if (X.Append2AndDestroy(FL,VC,'u')) wblog(FL,
            "ERR %s() VS: QIDX must not overlap\n%s",FCT,istr
         );
      }

   }

   UQ.checkQ(FL); if (SQ) SQ->checkQ(FL);
   VC.checkQ(FL);

   if (revert_new_idx) {
      unsigned r1=UQ.rank(FL); wbindex ia(1);
      wbvector< QSpace<TQ,TD> >X(1);
      QSpace<TQ,TD> U1J;

      if (UQ.qtype!=VC.qtype || (!UQ.qtype.allAbelian() &&
         (!UQ.CGR.data || !VC.CGR.data))) wblog(FL,"ERR %s() "
         "got uninitialzied CGR data\n(0x%lx, 0x%lx having %s)",
         FCT, UQ.CGR.data, VC.CGR.data, STR(VC.qtype));
      if (VC.itags.len!=2) wblog(FL,"ERR %s() "
         "unexpected rank-%d for %s",FCT,VC.itags.len, SQ ? "V":"VS");

      if (SQ) {
         SQ->save2(X[0]); U1J.initIdentityCG(X,ia,'z');
         U1J.contract(FL,1,X[0],ia[0]+1,*SQ); 
      }
      else {
         ia[0]=VC.rank(FL)-1;
         VC.save2(X[0]); U1J.initIdentityCG(X,ia,'z');
         X[0].contract(FL,ia[0]+1,U1J,1,VC); 
      }

      UQ.save2(X[0]); {
         X[0].contract(FL,r1,U1J,"1*",UQ); 
      }
   }

   return flag;
};

char mxGetLRDir( 
   const char *F, int L,
   const mxArray *a, unsigned k 
){
   char w=0, e=0; char s[8];

   if (mxIsChar(a)) {
      if (mxGetString(a,s,8)) {
         if (F) wblog(F,L,
            "ERR failed to read lrdir string (arg #%d)",k);
         return w;
      }

      if (!strcmp(s,">>")) w=+3; else
      if (!strcmp(s,"<<")) w=-3; else
      if (!strcmp(s,"LR")) w=+2; else
      if (!strcmp(s,"RL")) w=-2; else
      if (!strcmp(s,"+" )) w=+1; else   
      if (!strcmp(s,"-" )) w=-1; else { 
         if (F) wblog(F,L,
            "ERR invalid lrdir (arg #%d: '%s')",k,s);
         return w;
      }
   }
   else if (Mx::IsNumber(0,0,a)) {
      if (mxGetNumber(a,w)) e=1;
      else {
         if (!w || w<-3 || w>+3) { e=2; }
      }
      if (e) {
         if (F) wblog(F,L,
            "ERR invalid lrdir (arg #%d: %d; e=%d)",k,w,e);
         w=0; return w;
      }
   }
   else e=3;

   if (!w) e=4;
   if (e) {
      if (F) wblog(F,L,
         "ERR invalid lrdir (arg #%d; e=%d)",k,e);
      return w;
   }

   return w;
};

template <class TQ, class T1, class T2, class T3>
void twoSiteInit(
   const QSpace<TQ,T1> &A1,
   const QSpace<TQ,T2> &A2,
   QSpace<TQ,T3> &PSI, unsigned &K,
   widx_t ic1, widx_t ic2,
   char twoSite, char ldir
){
   const unsigned r1=A1.rank(), r2=A2.rank();
   wbindex I1, I2;
   wbperm P;

   if (!ldir) wblog(FL,"ERR ldir=%d",ldir);
   if ((ldir>0 && ic1>=r1) || (ldir<0 && ic2>=r2)) wblog(FL,
      "ERR contraction index out of bounds (%d/%d, %d/%d)",
      ic1,r1,ic2,r2);

   if (twoSite) {
      if (ldir>0) {
         A1.contract(ic1+1, A2, ic2+1, PSI); 
         I1.Index(r1-1); I2.Index(r2-1)+=(r1-1);
      }
      else {
         A2.contract(ic2+1, A1, ic1+1, PSI); 
         I1.Index(r2-1); I2.Index(r1-1)+=(r2-1);
      }
   }
   else {
      if (ldir>0)
           { PSI=A1; I1.Index(r1).Skip(ic1); I2.init(1,&ic1); }
      else { PSI=A2; I1.Index(r2).Skip(ic2); I2.init(1,&ic2); }
   }

   if (PSI.isEmpty()) wblog(FL,"WRN %s yields empty QSpace",FCT);

   P.Cat(I1,I2);
   PSI.Permute(P); K=I1.len;
};

template <class TQ, class T1, class T2, class T3>
void twoSiteFinal(
   const QSpace<TQ,T1> &Psi1,
   const QSpace<TQ,T2> &Psi2,
   QSpace<TQ,T3> &A1,
   QSpace<TQ,T3> &A2,
   unsigned ic1, unsigned ic2, char twoSite, char ldir
){
   QSpace<TQ,T3> AX;

   if (!ldir) wblog(FL,"ERR ldir=%d",ldir);

   if (twoSite) {
      if (ldir>0) {
         A1.PermuteLastTo (ic1);
         A2.PermuteFirstTo(ic2);
      }
      else {
         A2.permuteFirstTo(ic1, AX);
         A1.permuteLastTo (ic2, A2); AX.save2(A1);
      }
   }
   else {
      if (ldir>0) {
         A1.PermuteLastTo(ic1);
         A2.contract(2, Psi2, ic2+1, AX); 
         AX.permuteFirstTo(ic2, A2);
      }
      else {
         A2.contract(1, Psi1, ic1+1, AX); 
         A1.permuteLastTo(ic2, A2);
         AX.permuteFirstTo(ic1,A1);
      }
   }
}

template <class TQ, class TD>
mxArray* orthoQS(
   const QSpace<TQ,TD> &PSI, 
   QSpace<TQ,TD> &A1,
   QSpace<TQ,TD> &A2,
   unsigned K,     
   unsigned Nkmin, 
   unsigned Nkeep, 
   double stol,    
   itag_ tx,       
   const char *info
){
   mxArray *Sout=NULL; 
   unsigned i=0;
   char vflag=0, sflag=0;

   SVD_Data<TQ,TD> SVD;
   wbMatrix<double> SV;  

   if (info) {
      for (;; ++i) {
         if (info[i]=='d') { vflag|=1; } else
         if (info[i]=='D') { vflag|=2; } else
         if (info[i]=='s') { sflag|=1; } else
         if (info[i]=='S') { sflag|=2; } else break;
      }
      if (i>3 || info[i])
      wblog(FL,"ERR %s() invalid flags '%s' (i=%d)",FCT,info,i);
   }

   if (sflag)
   Sout=mxCreateStructMatrix(1,1,0,NULL);

   SVD.blockSVD(PSI,K); 

   i=SVD.dmrgTruncate(Nkmin,Nkeep,stol,tx,SV, A1,NULL,A2, Sout,vflag);

   if (Sout) {
      mxAddField2Scalar(FL,Sout,"K",  numtoMx(K));
      mxAddField2Scalar(FL,Sout,"D1", A1 .getDim().toMx());
      mxAddField2Scalar(FL,Sout,"DD", PSI.getDim().toMx());
   }
   else if (i) wblog(FL,str);

   if (sflag>1)
      mxAddField2Scalar(FL,Sout,"SVD", SVD.toMx());
   return Sout;
};

template <class TQ, class TD>
mxArray* getSVD(
   const QSpace<TQ,TD> &PSI, 
   QSpace<TQ,TD> &U,
   QSpace<TQ,TD> &S,
   QSpace<TQ,TD> &V,
   unsigned K,     
   unsigned Nkmin, 
   unsigned Nkeep, 
   double stol,    
   itag_ tx,       
   const char *info
){
   SVD_Data<TQ,TD> SVD;
   wbMatrix<double> SV;  

   mxArray *Sout=NULL;
   unsigned i=0; char sflag=0, vflag=0;

   if (info) {
      for (;; ++i) {
         if (info[i]=='d') { vflag|=1; } else
         if (info[i]=='D') { vflag|=2; } else
         if (info[i]=='s') { sflag|=1; } else
         if (info[i]=='S') { sflag|=2; } else break;
      }
      if (i>3 || info[i])
      wblog(FL,"ERR %s() invalid flags '%s' (i=%d)",FCT,info,i);
   }

   if (sflag)
      Sout=mxCreateStructMatrix(1,1,0,NULL);

   SVD.blockSVD(PSI,K); 

   i=SVD.dmrgTruncate(Nkmin,Nkeep,stol,tx,SV, U,&S,V, Sout,vflag);

   if (Sout) {
      mxAddField2Scalar(FL,Sout,"K", numtoMx(K));
      mxAddField2Scalar(FL,Sout,"D1",  U.getDim().toMx());
      mxAddField2Scalar(FL,Sout,"DD",PSI.getDim().toMx());
   }
   else if (i) wblog(FL,str);

   if (sflag>1)
      mxAddField2Scalar(FL,Sout,"SVD", SVD.toMx());
   return Sout;
};

template <class TQ, class TD>
mxArray* SVD_Data<TQ,TD>::toMx() const {

   const char *fn[]={"Q1","Q2","S1","S2","U","S","Vc","A","I0","Ik"};

   mxArray *S=mxCreateStructMatrix(QB.len,1,10,fn);

   for (unsigned i=0; i<QB.len; ++i) {
      mxSetFieldByNumber(S,i,0, QB[i].Q1.toMx() );
      mxSetFieldByNumber(S,i,1, QB[i].Q2.toMx() );
      mxSetFieldByNumber(S,i,2, QB[i].S1.toMx() );
      mxSetFieldByNumber(S,i,3, QB[i].S2.toMx() );
      mxSetFieldByNumber(S,i,4, QB[i].U .toMx() );
      mxSetFieldByNumber(S,i,5, QB[i].S .toMx() );
      mxSetFieldByNumber(S,i,6, QB[i].Vc.toMx() );
      mxSetFieldByNumber(S,i,7, QB[i].A .toMx() );
      mxSetFieldByNumber(S,i,8, QB[i].I0.toMx() );
      mxSetFieldByNumber(S,i,9, QB[i].Ik.toMx() );
   }

   return S;
}

#endif
