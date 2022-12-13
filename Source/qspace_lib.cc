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

#ifndef __WB_QSPACE_LIB_HCC__
#define __WB_QSPACE_LIB_HCC__

/* =======================================================================
 * #include "qspace_lib.cc"
 * -----------------------------------------------------------------------
 * contractIndex.cc (library file)
 *
 *    C++ program that gets index sets for generic contractQ()
 *
 * This is a MEX-file for MATLAB.
 * AW (C) Jan 2006
 * ======================================================================= */

void check_qdim     (C_MEX *C, unsigned &ma, wbvector<widx_t> &dq);
void check_cidx     (C_MEX *a, const unsigned &adim, UVEC &ica);
void getPerm        (C_MEX *a, const unsigned m, wbvector<unsigned> &P);
char isUInt         (const wbvector<double> &dd);

void getDataSize    (const mxArray *AD, const unsigned nd, const unsigned ma,
                     wbMatrix<unsigned> &usa);

void getQsub        (C_MEX *C, unsigned m, C_UVEC &dQ, C_UVEC &I, DMAT &Qc);
void getQall        (C_MEX *A, wbMatrix<double> &Q);
void getQall        (C_MEX *C, unsigned ma, C_UVEC &dq, wbMatrix<double> &Q);
void getQBlock      (C_MEX *C, C_UINT ic, DMAT &Qc);
void cell2mat       (C_MEX *C, unsigned *idx, unsigned n, DMAT &Qc);

char isUnique       (wbMatrix<double> Q);
char isUniqueSorted (C_DMAT &Q);

void makeUnique     (DMAT &Q, UVEC &Iu);
void makeUnique     (DMAT &Q);

int matchSortedIdx_old(
   C_DMAT &QA, C_DMAT &QB,
   wbvector<widx_t> &Ia, wbvector<widx_t> &Ib
);

int matchIndex_old(C_DMAT &QA, C_DMAT &QB,
   wbvector<widx_t> &Ia, wbvector<widx_t> &Ib
);

void QtoCell        (C_DMAT &Q, C_UVEC &dab, mxArray *&S, int fieldnr,
                     C_UVEC &P = wbvector<unsigned>());
void Data2vec       (mxArray *AD, C_IMAT sda, wbvector<double> &Avec);
void Sab2Cell       (C_IMAT &SA, C_IMAT &SB, C_UVEC &Ia, C_UVEC &Ib,
                     C_UVEC &Iu, C_UVEC &ika, C_UVEC &ikb,
                     mxArray* C, int fid1, int fid2);

void DtoCell        (C_IMAT &SA, C_IMAT &SB, C_UVEC &Ia, C_UVEC &Ib,
                     C_UMAT &IG, C_UVEC &ika, C_UVEC &ikb,
                     mxArray* C, int fid1, int fid2);

void check_qdim(
    const mxArray *C,
    unsigned &ma,
    wbvector<widx_t> &dq
){
    mxArray *a;
    unsigned k,da;

    if (!C || !mxIsCell(C) || mxIsEmpty(C))
    usage(FLINE, "Input arguments 1 and 2 must be cell arrays.");

    if (mxGetM(C)!=1)
    usage(FLINE, "Input arguments 1 and 2 must be cell arrays (row vectors).");

    da=mxGetN(C); dq.init(da);

    a=mxGetCell(C,0); dq[0]=mxGetN(a);
    ma=mxGetM(a);

    for (k=1; k<da; k++) {
        a=mxGetCell(C,k); dq[k]=mxGetN(a);

        if (ma!=(unsigned)mxGetM(a))
        usage(FLINE, "Q Space inconsistency in input arguments.");

        if (mxGetClassID(a)!=mxDOUBLE_CLASS)
        usage(FLINE, "Cell array of type double as input required.");
    }
}

void check_cidx(
    const mxArray *a, const unsigned &adim, wbvector<unsigned> &ica) {

    int e=0;
    try {
       Mx::Array<unsigned> A(a); ica.init(A.len);  
       if (A.len) {
          ica.init(A.len); A.copyTo(ica.data,'t'); 
       }
    }
    catch (...) { ++e; }
    if (e || !ica.len) usage(FLINE,"input arguments 3 and 4 must be vectors");

    unsigned i=0, n=ica.max(); char mark[n];

    for (; i<ica.len; ++i) {
       if (--ica[i]>n || ++mark[ica[i]]>1) wblog(FL,
       "ERR invalid contraction index (unique, 1-based)");
    }
};

void getPerm(
    const mxArray *a, const unsigned m0, wbvector<unsigned> &P) {

    Mx::Array<unsigned> A(a); 

    if (!A.len) { P.init(); return; }
    if (A.len!=m0) {
        sprintf(str,"permuation dimension mismatch (%ld,%d)",A.len,m0);
        usage(FLINE, str);
    }

    P.init(A.len); A.copyTo(P.data,'t'); 
    P-=1; 

    if (!validPerm(P)) {
       P.print("P"); wberror(FLINE,"invalid permuation");
    }
};

char isUInt(const wbvector<double> &dd) {

    double di;
    for (unsigned i=0; i<dd.len; i++) {
       di=dd[i];
       if (floor(di)!=di) return 0;
       if (di<0) return 0;
    }
    return 1;
}

void getDataSize(
    const mxArray *AD,
    const unsigned nd,
    const unsigned ma,
    wbMatrix<int> &sa
){
    unsigned i,n;
    mxArray *a;

    if (!mxIsCell(AD))
    wblog(FLINE,"ERR Data must be cell structure.");

    if (ma!=(n=(unsigned)mxGetNumberOfElements(AD)))
    wblog(FLINE,"ERR Dimension mismatch between Q and data (%d,%d)",ma,n);

    if (mxGetM(AD)!=1 && mxGetN(AD)!=1)
    wblog(FLINE,"ERR Q must be vectorized cell structure");

    sa.init(ma,nd); sa.set(1); 

    for (i=0; i<ma; i++) {
       a=mxGetCell(AD,i);

       n=mxGetNumberOfDimensions(a);
       if (n>nd) wblog(FLINE,"ERR Dimension mismatch in data (%d,%d)",n,nd);

       memcpy(sa.rec(i), mxGetDimensions(a), n*sizeof(int));
    }

    return;
}

void getQsub(
    const mxArray *C,
    unsigned ma,                  
    const wbvector<unsigned> &dQ, 
    const wbvector<unsigned> &I,  
    wbMatrix<double> &Qc          
){
    unsigned i,j,k,m,ib,nq,nC=mxGetNumberOfElements(C);
    double *qik;

    for (m=ib=0; ib<I.len; ib++) m += dQ[I[ib]];

    Qc.init(ma,m);

    for (k=ib=0; ib<I.len; ib++) {
        if (I[ib]>=nC) wberror(FLINE,"index out of bounds");

        Mx::Array<double> A(FL,mxGetCell(C,I[ib]));
        nq=dQ[I[ib]]; 

        for (i=0; i<ma; ++i) { qik=Qc.rec(i)+k;
        for (j=0; j<nq; ++j) { qik[j]=A.data[i+ma*j]; }} 

        k+=nq;
    }
};

void getQall(
    const mxArray *A,
    wbMatrix<double> &Q           
){
    unsigned ma;
    int iQ;

    mxArray *AQ;
    wbvector<widx_t> dqa;

    if (!mxIsStruct(A)) usage(FLINE,
    "Input argument must be structure {Q [,data,...]}.");

    if (mxGetNumberOfElements(A)!=1) usage(FLINE,
    "Input argument must be SINGLE structure {Q [,data,...]}.");

    iQ=mxGetFieldNumber(A,"Q");
    if (iQ<0) usage(FLINE,
    "Input argument must be structure of the type {Q [,data,...]}.");

    AQ=mxGetFieldByNumber(A,0,iQ);

    check_qdim(AQ, ma, dqa);

    getQall(AQ, ma, dqa, Q);

    return;
}

void getQall(
    const mxArray *C,
    unsigned ma,                  
    const wbvector<unsigned> &dq, 
    wbMatrix<double> &Q           
){
    unsigned i, j, k, D, ib, d2;
    double *q;

    D=dq.sum(); Q.init(ma,D);

    if (dq.len!=(unsigned)mxGetNumberOfElements(C))
    wberror(FLINE,"Dimension mismatch.");

    for (k=ib=0; ib<dq.len; ++ib) {
        Mx::Array<double> A(FL,mxGetCell(C,ib));
        d2=dq[ib]; 

        for (i=0; i<ma; ++i) { q=Q.rec(i)+k;
        for (j=0; j<d2; ++j) { q[j]=A.data[i+ma*j]; }} 

        k+=d2;
    }
};

void getQBlock(
    const mxArray *C,
    const unsigned ic,    
    wbMatrix<double> &Qc  
){
    unsigned i,j; double *qi;

    if (ic>=(j=mxGetNumberOfElements(C))) wblog(FL,
       "ERR %s() index out of bounds (%d/%d)",FCT,ic,j);

    Mx::Array<double> A(FL,mxGetCell(C,ic));
    if (A.rank!=2) wblog(FL,"ERR %s() invalid rank-%d object",FCT,A.rank);

    const size_t *sz=mxGetDimensions(A.ax);
    Qc.init(sz[0],sz[1]);

    for (i=0; i<sz[0]; ++i) { qi=Qc.rec(i);
    for (j=0; j<sz[1]; ++j) { qi[j]=A.data[i+sz[0]*j]; }} 

    return;
};

void cell2mat(
    const mxArray *C,
    const wbMatrix<unsigned> &I,
    wbMatrix<double> &M          
){
    unsigned i,j,k,d1,d2,N,r0,s0,r,s;
    double *d;
    mxArray *a;

    wbvector<unsigned> D1(I.dim1), D2(I.dim2);
    wbMatrix<double*> dd(I.dim1, I.dim2);

    N=(unsigned)mxGetNumberOfElements(C);

    for (i=0; i<I.dim1; i++)
    for (j=0; j<I.dim2; j++) {
        k=I(i,j);
        if (k>=N) wberror(FLINE, "Index out of bounds.");
        a=mxGetCell(C,k);

        if (i==0) D2[j]=(unsigned)mxGetN(a);
        else if (D2[j]!=(unsigned)mxGetN(a))
        wberror(FLINE, "Dimension mismatch.");

        if (j==0) D1[i]=(unsigned)mxGetM(a);
        else if (D1[i]!=(unsigned)mxGetM(a))
        wberror(FLINE, "Dimension mismatch.");

        dd(i,j)=Mx::Array<double>(FL,a).data;
    }

    M.init(D1.sum(), D2.sum());

    for (r0=i=0; i<I.dim1; i++) {
        d1=D1[i];
        for (s0=j=0; j<I.dim2; j++) {
             d2=D2[j]; d=dd(i,j);

             for (k=s=0; s<d2; s++)
             for (  r=0; r<d1; r++) M(r0+r,s0+s)=d[k++];

             s0+=d2;
        }
        r0+=d1;
    }
};

char isUnique(wbMatrix<double> Q) {

   wbperm iS;
   Q.SortRecs(iS);
   return isUniqueSorted(Q);
}

char isUniqueSorted(const wbMatrix<double> &Q) {

   unsigned j;

   for (j=1; j<Q.dim1; j++)
   if (Q.recEqual(j-1,j)) return 0;

   return 1;
}

void makeUnique(wbMatrix<double> &Q, wbvector<unsigned> &Iu) {

   unsigned i,k,m,ig,d, ng;
   wbvector<widx_t> dQ;
   wbperm iS;

   Q.groupRecs(iS,dQ); ng=dQ.len;

   Iu.init(iS.len+ng);

   for (k=ig=0; ig<ng; ig++) {
       d=dQ[ig];
       Iu[k+ig]=d; m=ig+1;
       for (i=0; i<d; i++, k++) Iu[k+m]=iS[k];
   }
}

void makeUnique(wbMatrix<double> &Q) {
   wbperm iS; wbvector<widx_t> dQ;
   Q.groupRecs(iS,dQ);
}

int matchIndex_old(
    const wbMatrix<double> &QA,
    const wbMatrix<double> &QB,
    wbvector<widx_t> &Ia,
    wbvector<widx_t> &Ib
){
    wbMatrix<double> Q1(QA), Q2(QB);
    wbvector<widx_t> ix1, ix2;
    wbperm is1, is2;

    Q1.SortRecs(is1);
    Q2.SortRecs(is2);

    matchSortedIdx_old(Q1,Q2,ix1,ix2);

    is1.get(ix1,Ia);
    is2.get(ix2,Ib);

    return 0;
}

int matchSortedIdx_old(
    const wbMatrix<double> &QA,
    const wbMatrix<double> &QB,
    wbvector<widx_t> &Ia,
    wbvector<widx_t> &Ib
){

    unsigned i,j,k, ng, ia,ib, ig, iga, igb, m1,m2, mtot=0;
    unsigned *p1, *p2, na=QA.dim1, nb=QB.dim1;
    char c;

    if (QA.dim2!=QB.dim2)
    wberror(FLINE, "Dimension mismatch");

    wbvector<unsigned> gra(2*na), grb(2*nb);

    for (ng=iga=igb=ia=ib=0; ia<na && ib<nb;) {
        c=QA.recCompareP(ia, QB.rec(ib));
        if (c<0) ia++; else
        if (c>0) ib++;
        else {
            p1=gra.data+(iga++); gra[iga++]=(ia++);
            p2=grb.data+(igb++); grb[igb++]=(ib++);

            m1=m2=1; ng++;
            while (ia<na) {
                if (QA.recEqual(ia-1, ia)) { gra[iga++]=(ia++); m1++; }
                else break;
            }
            while (ib<nb) {
                if (QB.recEqual(ib-1, ib)) { grb[igb++]=(ib++); m2++; }
                else break;
            }
            mtot+=m1*m2; (*p1)=m1; (*p2)=m2;
        }
    }

    Ia.init(mtot); Ib.init(mtot);

    for (k=ig=iga=igb=0; ig<ng; ig++) {
        m1=gra[iga]; p1=gra.ref(iga+1); iga+=(m1+1);
        m2=grb[igb]; p2=grb.ref(igb+1); igb+=(m2+1);
        for (i=0; i<m1; i++)
        for (j=0; j<m2; j++,k++) {
            Ia[k]=p1[i];
            Ib[k]=p2[j];
        }
    }

    return 0;
}

int matchUSortedIdx(
    const wbMatrix<double> &QA,
    const wbMatrix<double> &QB,
    wbvector<unsigned> &Ia,
    wbvector<unsigned> &Ib
){

    unsigned n,ia,ib; char c;

    if (QA.dim2!=QB.dim2)
    wberror(FLINE, "Dimension mismatch");

 #ifdef __WBDEBUG__
    if (!QA.isUniqueSorted('A') || !QB.isUniqueSorted('A'))
    wberror(FL,"ERR Input must be unique increasing set of records.");
 #endif

    n=MIN(QA.dim1,QB.dim1); Ia.init(n); Ib.init(n);

    for (n=ia=ib=0; ia<QA.dim1 && ib<QB.dim1 ;) {

        c=QA.recCompareP(ia, QB.rec(ib));

        if (c<0) ia++; else
        if (c>0) ib++;
        else {
            Ia[n]=ia;
            Ib[n]=ib; ia++; ib++; n++;
        }
    }

    Ia.Resize(n); Ib.Resize(n);

    return 0;
}

void Sab2Cell(
    const wbMatrix<int> &SA,
    const wbMatrix<int> &SB,
    const wbvector<unsigned> &Ia,
    const wbvector<unsigned> &Ib,
    const wbvector<unsigned> &Id,
    const wbvector<unsigned> &ika,
    const wbvector<unsigned> &ikb,
    mxArray* ST, int fid_SAB, int fid_Sab
){
    unsigned i,k,l,ia,ib,ng=Id.len,d2=ika.len+ikb.len;
    wbvector<int> ss, sk(d2);
    double *d, *SD, *sd;

    mxSetFieldByNumber(ST,0,fid_SAB, Mx::Array<double>( 2,ng).Return(SD));
    mxSetFieldByNumber(ST,0,fid_Sab, Mx::Array<double>(d2,ng).Return(sd));

    for (l=k=0; k<ng; k++) {

        ia=Ia[l]; ib=Ib[l]; l+=Id[k];

        if (ia>=SA.dim1 || ib>=SB.dim1)
        wberror(FLINE,"Index out of bounds.");

        ss.init(SA.dim2, SA.rec(ia)); ss.get(ika, sk.data);
        ss.init(SB.dim2, SB.rec(ib)); ss.get(ikb, sk.data+ika.len);

        i=k*2;
        SD[i  ]=(double)sk.prod(unsigned(0),ika.len-1);
        SD[i+1]=(double)sk.prod(ika.len,sk.len-1);

        d=sd+k*d2;

        for (i=0; i<d2; i++) d[i]=(double)sk[i];
    }

    return;
}

void QtoCell(
    const wbMatrix<double> &Q,
    const wbvector<unsigned> &dab,
    mxArray *&S,
    int field_nr,
    const wbvector<unsigned> &P
){
    unsigned i,j,k,m,d1=Q.dim1,d2,n=dab.len;
    double *dd;
    mxArray *C;

    wbvector<unsigned> D2(n);

    if (field_nr>=0) {
        mxSetFieldByNumber(S,0,field_nr, mxCreateCellMatrix(1,n));
        C=mxGetFieldByNumber(S,0,field_nr);
    } else { S=mxCreateCellMatrix(1,n); C=S; }

    if (Q.dim1==0 || dab.len==0) {
        for (i=0; i<n; i++)
        mxSetCell(C,i,mxCreateDoubleMatrix(0,dab[i],mxREAL));
        return;
    }

    D2[0]=0; for (j=1; j<n; j++) D2[j]=D2[j-1]+dab[j-1];

    if (D2[j-1]+dab[j-1]!=Q.dim2) wberror(FLINE, "Dimension mismatch.");
    if (P.len && P.len!=n)        wberror(FLINE, "Invalid permutation.");

    for (k=0; k<n; k++) {
        d2=dab[k];

        mxSetCell(C,k, Mx::Array<double>(d1,d2).Return(dd));
        m = P.len ? D2[P[k]] : D2[k];

        for (j=0; j<d2; ++j)
        for (i=0; i<d1; ++i) { dd[d1*j+i]=Q(i,m+j); } 
    }
};

void Data2vec(
    mxArray *AD, const wbMatrix<int> sda, wbvector<double> &Avec) {

    unsigned i,d,k,N;
    wbvector<int> SD;

    if (sda.dim1==0 || sda.dim2==0) { Avec.init(); return; }

    SD.init(sda.dim1);
    for (i=0; i<SD.len; i++) SD[i]=sda.recProd(i);

    N=(unsigned)SD.sum();

    Avec.init(N);
    for (k=i=0; i<SD.len; ++i) { d=(unsigned)SD[i];
        Mx::Array<double> A(FL,mxGetCell(AD,i));
        memcpy(Avec.data+k, A.data, d*sizeof(double));
        k+=d;
    }
};

void QI2Cell(
    const wbMatrix<double> &Q,
    const mxArray *AQ,
    mxArray* C,
    int field_nr,
    const wbvector<widx_t> &P
){
    unsigned i,m,M;
    wbvector<widx_t> iP, dqa, iv;
    wbindex iq,ir;

    wbMatrix<double> Qref;

    if (Q.data==NULL) {
        mxSetFieldByNumber(C, 0, field_nr,
        mxCreateDoubleMatrix(0,0,mxREAL)); return;
    }

    getIPerm(P,iP);

    check_qdim(AQ,m,dqa);
    getQsub   (AQ,m,dqa, iP, Qref);

    M=Q.dim1; iv.init(M);

    if (Q.dim2==Qref.dim2) {
        matchIndex(Q,Qref,iq,ir); m=iq.len;

        for (i=0; i<m; i++) iv[iq[i]]=ir[i]+1;
    }

    mxSetFieldByNumber(C,0,field_nr,
       Mx::Array<double>(M,1).copyFromTR(iv.data)
    );
};

void DtoCell(
    const wbMatrix<int> &SA,
    const wbMatrix<int> &SB,
    const wbvector<unsigned> &Ia,
    const wbvector<unsigned> &Ib,
    const wbMatrix<unsigned> &IG,
    const wbvector<unsigned> &ika,
    const wbvector<unsigned> &ikb,
    mxArray* ST, int fid1, int fid2
){
    unsigned k,l,ia,ib,ng=IG.dim1;
    wbvector<size_t> uu, sk(ika.len+ikb.len);

    mxArray *D=mxCreateCellMatrix(ng,1);
    double *sd;

    mxSetFieldByNumber(ST,0,fid1, D);
    mxSetFieldByNumber(ST,0,fid2, Mx::Array<double>(ng,2).Return(sd));

    for (k=0; k<ng; ++k) {

        l=IG(k,1); ia=Ia[l]; ib=Ib[l]; 

        if (ia>=SA.dim1 || ib>=SB.dim1)
        wberror(FLINE,"Index out of bounds.");

        uu.initT(SA.dim2, SA.rec(ia)); uu.get(ika, sk.data);
        uu.initT(SB.dim2, SB.rec(ib)); uu.get(ikb, sk.data+ika.len);

        sd[k   ]=(double)sk.prod(unsigned(0),ika.len-1);
        sd[k+ng]=(double)sk.prod(ika.len,sk.len-1);

        mxSetCell(D,k, Mx::Array<double>(sk).Return());
    };
};

void Sab2Site2Cell(
    const wbMatrix<int> &SA,
    const wbMatrix<int> &SB,
    const wbvector<unsigned> &Ia,
    const wbvector<unsigned> &Ib,
    const wbMatrix<unsigned> &IG,
    const wbvector<unsigned> &ika,
    const wbvector<unsigned> &ikb,
    const wbvector<unsigned> &ikP,
    mxArray* ST, int fid_Sab, int fid_SAB
){
    unsigned k,l,ia,ib, ng=IG.dim1, dk=ika.len+ikb.len; int *s1,*s2;
    wbMatrix<int> sd(ng,dk);
    wbvector<int> uu;
    double *SD;

    mxSetFieldByNumber(ST,0,fid_SAB, Mx::Array<double>(2,ng).Return(SD));

    for (k=0; k<ng; ++k) {

        l=IG(k,1); ia=Ia[l]; ib=Ib[l]; 

        if (ia>=SA.dim1 || ib>=SB.dim1)
        wberror(FLINE,"Index out of bounds.");

        s1=sd.rec(k); s2=s1+ika.len;

        uu.init(SA.dim2, SA.rec(ia)); uu.get(ika, s1);
        uu.init(SB.dim2, SB.rec(ib)); uu.get(ikb, s2);

        SD[2*k  ]=(double)Wb::prodRange(s1,ika.len);
        SD[2*k+1]=(double)Wb::prodRange(s2,ikb.len);
    }

    sd.colPermute(ikP);
    sd.mat2mxs(ST,fid_Sab);
};

#endif

