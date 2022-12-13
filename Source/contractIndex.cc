char USAGE[] =
/* ==================================================================
 * contractIndex.cc
 *
 *    C++ program that gets index sets for generic contractQ() */

"   Usage: CIDX = contractIndex(QA,QB,i1,ib [,perm, [Qconst]])      \n\
                                                                    \n\
   perm is an optional permutation of dimensions after contraction. \n\
   Output: CIDX is a structure with elements                        \n\
                                                                    \n\
      Q       remaining (new) Q data after contraction              \n\
      Sab     remaining size vector of A after contraction          \n\
      SAB     remaining summed upsize of A after contraction        \n\
      Ia      index into original QA.data for contraction           \n\
      Ib      index into original QB.data for contraction           \n\
      Id      dimension of consecutive groups in Ia and Ib          \n\
      pA      permutation index needed for dimenxions of A data     \n\
      pB      permutation index needed for dimensions of B data     \n\
                                                                    \n\
   AW (C) Jan 2006                                                   \n";

/* This is a MEX-file for MATLAB.
 * ================================================================== */

#include <mex.h>
#include <math.h>
#include <string.h>
#include "wblib.h"

#include "qspace_lib.cc"

void mexFunction(
    int nargout, mxArray *plhs[],
    int nargin, const mxArray *prhs[]
){
    unsigned i,j,k,n,m,ma,mb,qiflag;
    int ifield_Q, ifield_D;

    wbMatrix<double> QAk, QAc, QBc, QBk, QAB, QAUX;
    wbvector<unsigned> ica, icb, ika, ikb, dqa, dqb, dab,
                       Ia, Ib, Iu, Id, uu, P;
    wbindex Isa,Isb;
    wbperm isa,isb;
    mxArray *AQ, *ADATA, *BQ, *BDATA;
    wbMatrix<int> sda, sdb;

    if (nargin) if (isHelpIndicator(prhs[0])) { usage(); return; }

    if (nargin<4 || nargin>6)
                   usage(FLINE, "Invalid number of input arguments.");
    if (nargout>1) usage(FLINE, "Max one output argument available.");

    if (!mxIsStruct(prhs[0]) || !mxIsStruct(prhs[1])) usage(FLINE,
    "First two input arguments must be structures {Q,data,...}.");

    if (mxGetNumberOfElements(prhs[0])!=1 ||
        mxGetNumberOfElements(prhs[1])!=1) usage(FLINE,
    "First two input arguments must be SINGLE structures {Q,data,...}.");

    ifield_Q=mxGetFieldNumber(prhs[0], "Q");
    ifield_D=mxGetFieldNumber(prhs[0], "data");

    if (ifield_Q<0 || ifield_D<0) usage(FLINE,
    "First two input arguments must be structures {Q,data,...}.");

    AQ=mxGetFieldByNumber(prhs[0], 0, ifield_Q);
    ADATA=mxGetFieldByNumber(prhs[0], 0, ifield_D);

    ifield_Q=mxGetFieldNumber(prhs[1], "Q");
    ifield_D=mxGetFieldNumber(prhs[1], "data");

    if (ifield_Q<0 || ifield_D<0) usage(FLINE,
    "First two input arguments must be structures {Q,data,...}.");

    BQ=mxGetFieldByNumber(prhs[1], 0, ifield_Q);
    BDATA=mxGetFieldByNumber(prhs[1], 0, ifield_D);

    check_qdim(AQ, ma, dqa);
    check_qdim(BQ, mb, dqb);

    getDataSize(ADATA,dqa.len,ma, sda);
    getDataSize(BDATA,dqb.len,mb, sdb);

    check_cidx(prhs[2], dqa.len, ica);
    check_cidx(prhs[3], dqb.len, icb);

    if (ica.len!=icb.len)
    usage(FLINE, "Mismatch of number of indizes to contract.");

    for (i=0; i<ica.len; i++)
    if (dqa[ica[i]]!=dqb[icb[i]]) {
        mexPrintf("%d: %d <> %d\n", i, dqa[ica[i]], dqb[icb[i]]);
        dqa.print("dqa"); ica.print("dqa");
        dqb.print("dqa"); icb.print("dqa");
        usage(FLINE, "Q-space mismatch of indizes to contract.");
    }

    Wb::invertIndex(ica, dqa.len, ika);
    Wb::invertIndex(icb, dqb.len, ikb);

    if (nargin>4)
    getPerm(prhs[4],ika.len+ikb.len,P);

    qiflag = (dqb.len==2);

    if (nargin>5)
    qiflag=getFlag(prhs[5]);

    getQsub(AQ,ma,dqa,ica,QAc);
    getQsub(AQ,ma,dqa,ika,QAk);

    getQsub(BQ,mb,dqb,icb,QBc);
    getQsub(BQ,mb,dqb,ikb,QBk);

    QAc.SortRecs(isa);
    QBc.SortRecs(isb);

    matchSortedIdx(QAc,QBc,Isa,Isb);

    isa.get(Isa,Ia);
    isb.get(Isb,Ib);

    QAk.getRecs(Ia,QAB);
    QBk.getRecs(Ib,QAUX);

    QAB.cat(2,QAUX);

    makeUnique(QAB,Iu); 

    Id.init(Iu.len); uu.init(Ia.len);

    for (n=k=i=0; i<Iu.len;) {
        Id[n++]=m=Iu[i++];
        for (j=0; j<m; j++) uu[k++]=Iu[i++];
    }

    Id.Resize(n);
    Ia.Permute(uu); Ib.Permute(uu);

    if (nargout==0) return;

    const char *fnames[]= {
       "Q", "SAB", "Sab", "Ia", "Ib", "Id", "pA", "pB" };

    plhs[0]=mxCreateStructMatrix(1,1, 8, fnames);

    dqa.get(ika,dab); dqb.get(ikb,uu); dab.cat(uu);

    QtoCell(QAB, dab, plhs[0], 0, P);

    if (qiflag) {
        int fid=mxAddField(plhs[0], "iQA");
        if (fid<0) wberror(FLINE,"Could not add field to structure.");
        QI2Cell(QAB, AQ, plhs[0], fid, P);
    }

    Sab2Cell(
       sda, sdb, Ia, Ib, Id, ika, ikb,
       plhs[0], 1, 2 
    );

    uu=Ia+1; uu.mat2mxs(plhs[0], 3, 1);  
    uu=Ib+1; uu.mat2mxs(plhs[0], 4, 1);
             Id.mat2mxs(plhs[0], 5, 1);

    uu=ika; uu.cat(ica); uu+=1; uu.mat2mxs(plhs[0],6);
    uu=icb; uu.cat(ikb); uu+=1; uu.mat2mxs(plhs[0],7);

    return;
}

