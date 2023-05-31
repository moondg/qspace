
/* COMMENTS / CHANGE LOG ============================================= *
 * =================================================================== */

char USAGE[]=""; // outsourced to matchIndex.m // Wb,Feb14,19

#include "wblib.h"

#include "qspace_lib.cc"

void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){ Wb::CleanUp aclu; try { 

    wbMatrix<double> QA,QB;
    wbMatrix<wbcomplex> ZA,ZB; 
    wbindex In1,In2,ia,ib;
    unsigned nan1,nan2;
    double eps=0;

    char sflag=0, fflag=0, qflag=0;

    MX_CHECK_HELPER_NARGS(2,-1,3); 

    if (nargin>2) {
       OPTS opts; opts.init(argin+2, nargin-2);

       fflag=opts.getOpt("-f");
       qflag=opts.getOpt("-q");
       sflag=opts.getOpt("-s");
       opts.getOpt("eps", eps);

       opts.checkAnyLeft(); 
    }

    if (!mxIsComplex(argin[0]) && !mxIsComplex(argin[1])) {
       QA.init(FL,argin[0]);
       QB.init(FL,argin[1]);

       if (QA.dim1==1 && QB.dim1==1 && QA.dim2 != QB.dim2) {
           SWAP(QA.dim1,QA.dim2);
           SWAP(QB.dim1,QB.dim2);
       }
    }
    else {
       ZA.init(FL,argin[0]);
       ZB.init(FL,argin[1]);

       if (ZA.dim1==1 && ZB.dim1==1 && ZA.dim2 != ZB.dim2) {
           SWAP(ZA.dim1,ZA.dim2);
           SWAP(ZB.dim1,ZB.dim2);
       }

       QA.init2ref(ZA.dim1,2*ZA.dim2,(double*)ZA.data);
       QB.init2ref(ZB.dim1,2*ZB.dim2,(double*)ZB.data);
    }

    if (QA.dim2!=QB.dim2 && QA.dim2 && QB.dim2) wblog(FL,
       "ERR column size mismatch (%d/%d)", QA.dim2, QB.dim2);

    if (fflag) {
       QA.SkipTiny_float(); 
       QB.SkipTiny_float();
    }

    nan1=QA.skipNanRecs(In1);
    nan2=QB.skipNanRecs(In2);

    if (!qflag) if (nan1 || nan2) wblog(FL,
       "WRN skipping %d record(s) containing NaN",nan1+nan2);

    matchIndex(QA,QB,ia,ib,1,NULL,NULL,eps);

    if (sflag) {
       wbperm p; Wb::hpsort(ia,p); ib.Select(p);

       for (unsigned i=1,j=0; i<ia.len; i++) {
          if (ia[i-1]!=ia[i]) {
             if (j+1<i) {
                wbindex jj; jj.init2ref(i-j, ib.data+j);
                Wb::hpsort(jj,p);
             }
             j=i; continue;
          }
       }
    }

    if (nan1) { In1.Select(ia); In1.save2(ia); }
    if (nan2) { In2.Select(ib); In2.save2(ib); }

    if (nargout>=3) {
       mxArray *S=mxCreateStructMatrix(1,1,0,NULL);
       wbindex iax,ibx; wbvector<unsigned> sz(2);

       Wb::invertIndex(ia, QA.dim1, iax, 0);
       Wb::invertIndex(ib, QB.dim1, ibx, 0);

       sz[0]=QA.dim1;
       sz[1]=QB.dim1;

       mxAddField2Scalar(FL,S,"n",   numtoMx(ia.len));
       mxAddField2Scalar(FL,S,"dim", sz.toMx());
       mxAddField2Scalar(FL,S,"ix1", iax.toMx());
       mxAddField2Scalar(FL,S,"ix2", ibx.toMx());

       argout[2]=S;
    }

    argout[0]=ia.toMx(); 
    if (nargout>1) argout[1]=ib.toMx();

}  catch (Wb::LogException &e) { ExitMsg(e.istr); }
   catch (...) { ExitMsg("caught exception in matchIndex"); }
   aclu.Check();
};

