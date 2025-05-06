/* ---------------------------------------------------------------------
 * Project : QSpace tensor library (v4.0)
 * Class   : QSpace
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

#ifndef __WB_QSpace_AUX_HH__
#define __WB_QSpace_AUX_HH__

template <class TQ, class TD> class QSpace;
template <class TQ, class TD> class CPAT;

template <class TQ, class TA, class TB, class TC>
void contractDATA_plain(
   const char *F, int L,
   const QSpace<TQ,TA> &A, C_UVEC &Ia, const ctrIdx &ica,
   const QSpace<TQ,TB> &B, C_UVEC &Ib, const ctrIdx &icb,
   QSpace<TQ,TC> &C, unsigned ic
);

template <class TQ, class TA, class TB, class TC>
double contractDATA_group(
   const char *F, int L,
   const QSpace<TQ,TA> &A, C_UVEC &Ia, const ctrIdx &ica,
   const QSpace<TQ,TB> &B, C_UVEC &Ib, const ctrIdx &icb,
   QSpace<TQ,TC> &C, unsigned ic,
   char preview=0
);

template<class TQ, class TD> 
void mxInitQSpaceVec(
   const char *F, int L, const mxArray* C,
   wbvector< QSpace<TQ,TD> > &Fk,
   const char ref=0,
   const unsigned *rmin=NULL, 
   const unsigned *rmax=NULL  
);

template<class TQ, class TD> 
void mxInitQSpaceVecVec(
   const char *F, int L, const mxArray* C,
   wbvector< wbvector< QSpace<TQ,TD> > > &Fk, const char ref=0
);

template<class TQ, class TD> 
void mxcInitQSpaceVec(
   const char *F, int L, const mxArray* C,
   wbvector< QSpace<TQ,TD> > &Fk, const char ref=0
);

template<class TQ, class TD>
mxArray* QSpaceVec2Mx(
   const wbvector< QSpace<TQ,TD> > &Fk
);

template<class TQ, class TD>
mxArray* QSpaceVec2Mx(
   const wbvector< QSpace<TQ,TD> > &F, wbindex &I
);

template<class TQ, class TD>
mxArray* QSpaceVecVec2Mx(
   const wbvector< wbvector< QSpace<TQ,TD> > > &F,
   const char ref=0
);

template<class TQ, class TD>
mxArray* QSpaceVec2Mxc(
   const wbvector< QSpace<TQ,TD> > &F,
   const char ref=0
);

template<class TQ, class TD>
void putQSpaceVec(
   const wbvector< QSpace<TQ,TD> > &Fk,
   const char *vname, const char *ws="caller"
){
   mxArray *a=Fk.toMx();
   mxPutAndDestroy(FL,a,vname,ws);
};

template <class TQ, class TD>
void getDiffOp(
   const QSpace<TQ,TD> &A, unsigned ia,
   const QSpace<TQ,TD> &B, unsigned ib, const QSpace<TQ,TD> &C
);

template<class TQ, class TD>
void getQall(const wbvector< QSpace<TQ,TD> > &F, wbMatrix<TQ> &QA);

int mxIsQSpace( 
   const char *F, int L,
   const mxArray *a, unsigned &rank,
   char cflag=0,
   unsigned k=-1, 
   const unsigned *rmin=NULL, 
   const unsigned *rmax=NULL, 
   const char *istr=NULL
);

inline int mxIsQSpace(
   const char *F, int L, const mxArray *a, char cflag=0,
   const char *istr=NULL
){
   unsigned r=-1;
   return mxIsQSpace(F,L,a,r,cflag,-1,NULL,NULL,istr);
};

inline int mxIsQSpace(
   const mxArray *a, char cflag=0, const char *istr=NULL
){
   unsigned r=-1;
   return mxIsQSpace(NULL,0,a,r,cflag,-1,NULL,NULL,istr);
};

bool mxIsQSpaceArr(
   const char *F, int L,
   const mxArray *S,
   unsigned rank=-1,
   int arrdim=2, 
   const unsigned *rmin=NULL, 
   const unsigned *rmax=NULL, 
   char cflag=0
);

bool mxIsQSpaceVecVec(const mxArray *C, unsigned r=-1, char cflag=0);

inline bool mxIsQSpaceVec(
   const char *F, int L,
   const mxArray *S,
   unsigned rank=-1,
   const unsigned *rmin=NULL, 
   const unsigned *rmax=NULL, 
   char cflag=0
){ return mxIsQSpaceArr(F,L,S,rank,1,rmin,rmax,cflag); }; 

inline bool mxIsQSpaceVecR23( 
   const char *F, int L,
   const mxArray *S,
   char cflag=0
){
   unsigned rank=-1, rmin=2, rmax=3;
   return mxIsQSpaceArr(F,L,S,rank,1,&rmin,&rmax,cflag);
};

inline bool mxIsQSpaceMat(
   const char *F, int L,
   const mxArray *S,
   unsigned rank=-1,
   const unsigned *rmin=NULL, 
   const unsigned *rmax=NULL, 
   char cflag=0
){ return mxIsQSpaceArr(F,L,S,rank,2,rmin,rmax,cflag); };

bool mxIsEmptyQSpace(const mxArray *a, unsigned k=0);

bool mxIsScalarQSpace(const mxArray *a, unsigned k=0);

inline bool mxIsQSpaceScalar(
   const char *F, int L,
   const mxArray *S,
   unsigned rank=-1,
   const unsigned *rmin=NULL, 
   const unsigned *rmax=NULL, 
   char cflag=0
){ return mxIsQSpaceArr(F,L,S,rank,0,rmin,rmax,cflag); };

bool mxsIsQSpaceVec(const mxArray *S, int fid, unsigned r=-1, char cflag=0);
bool mxsIsQSpaceVEC(const mxArray *S, int fid, unsigned r=-1, char cflag=0);

int mxIsQSpaceVec(
   const char *F, int L,
   const char *fname, const char *vname,
   unsigned rank=-1,
   const unsigned *rmin=NULL, 
   const unsigned *rmax=NULL, 
   char cflag=0,
   unsigned N=-1 
);

int mxIsQSpaceVEC(
   const char *F, int L,
   const char *fname, const char *vname,
   unsigned rank=-1,
   const unsigned *rmin=NULL, 
   const unsigned *rmax=NULL, 
   char cflag=0
);

bool mxIsQSpaceVecOrEmpty(const mxArray *S, int rank, char cflag=0) {
   if (mxIsEmpty(S)) return 1;
   return mxIsQSpaceVec(FL,S,rank,NULL,NULL,cflag);
};

inline int mxIsQSpaceOrEmpty(
   const mxArray *a, unsigned rank=-1, char cflag=0, unsigned k=-1
){ if (mxIsEmpty(a)) return 1;
   else return mxIsQSpace(NULL,0,a,rank,cflag,k);
};

inline int mxIsQSpaceOrEmpty(const char *F, int L,
   const mxArray *a, unsigned rank=-1, char cflag=0, unsigned k=-1
){ if (mxIsEmpty(a)) return 1;
   else return mxIsQSpace(F,L,a,rank,cflag,k);
};

QVec qsGetSym( 
   const char *F, int L, const mxArray *a, unsigned k=-1);

char qsGotCGS( 
   const char *F, int L, const mxArray *a, unsigned k=-1
){
   QVec qtype=qsGetSym(F_L,a,k);
   return qtype.permitsOM();
};

#endif

