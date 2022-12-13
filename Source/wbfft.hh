
#ifndef __WB_FFTW_C__
#define __WB_FFTW_C__

#include <rfftw.h>

namespace Wb {

int FFT2D (
    wbMatrix<double> &X2,
    wbMatrix<double> &F2,
    char dir
){

    int info=-1, thesame=(&X2==&F2); 
    wbvector<int> dim(2);

    if (dir=='F') {
        if (X2.isEmpty()) wblog(FL,
        "ERR %s - invalid dim [%d,%d]",FCT,X2.dim1,X2.dim2);

        dim[0] =  X2.dim1;
        dim[1] =  X2.dim2;

        if (!thesame) F2=X2;

        F2.Resize(dim[0],dim[1]+2);

        info=FFT_ND(dim,F2[0],dir);
    }
    else
    if (dir=='B') {
        if (F2.isEmpty()) wblog(FL,
        "ERR %s - invalid dim [%d,%d]",FCT,F2.dim1,F2.dim2);

        dim[0] =  F2.dim1;
        dim[1] =  F2.dim2-2;

        if (!thesame) X2 = F2;

        info = FFT_ND (dim, X2[0], dir);

        X2.Resize(dim[0], dim[1]);
    }
    else
    wblog(FL,"ERR %s - invalid dir=%c<%d>",FCT,dir,dir);

    return info;
}

int FFT_ND (
    wbvector<int> &dim,
    double *zz,
    char dir
){

    unsigned i; int s=dim.prod(), info=-1;

    if (dim.len==0 || s<=0) wblog(FL,
    "ERR %s - invalid dimensions [%s]",FCT,dim.toStr().data);

    if (dim.last()%2) wblog(FL,
       "ERR %s - invalid dimensions [%s]\n"
       "leading/last dimension must be even!",FCT,dim.toStr().data);

    if (dir=='F') {

        rfftwnd_plan plan;

        plan = rfftwnd_create_plan (
        dim.len, dim.data, FFTW_FORWARD,  FFTW_ESTIMATE | FFTW_IN_PLACE);

        rfftwnd_one_real_to_complex (plan, zz, NULL);

        rfftwnd_destroy_plan(plan);
        info=0;
    }

    else if (dir=='B') {

        double fact; int ntot, nlast=dim.last();
        rfftwnd_plan plan;

        plan = rfftwnd_create_plan (
        dim.len, dim.data, FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);

        rfftwnd_one_complex_to_real (plan, (fftw_complex*)zz, NULL);

        for (ntot=1,i=0; i<dim.len; i++) ntot *= dim[i];
        fact = 1./ntot;
        ntot = (ntot / nlast) * (nlast+2);
        for (i=0; i<ntot; i++) zz[i] *= fact;

        rfftwnd_destroy_plan(plan);
        info = 0;
    }
    else wblog(FL,"ERR %s - invalid dir=%c<%d>",dir,dir);

    return info;
}

} 

#endif

