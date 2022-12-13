/* ---------------------------------------------------------------------
 * Project : QSpace tensor library (v4.0 pre-release)
 * Class   : wbMatrix (matrix class, row-major)
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

#ifndef __WB_MATRIX_BLAS_HH__
#define __WB_MATRIX_BLAS_HH__

// ----------------------------------------------------------------- //
// ----------------------------------------------------------------- //

void WbEigenSymmetric (const wbMatrix<double>&M,
     wbMatrix<double>&U, wbvector<double> &E);

namespace Wb {

template<class T> 
void MatProd(
    const wbMatrix<T> &A,
    const wbMatrix<T> &B, wbMatrix<T> &C,
    char aflag='N', char bflag='N',
    const T afac=1., const T cfac=0.,
    const char i00flag=0
);

template<class TA, class TB, class TC> 
void MatProd(
    const wbMatrix<TA> &A,
    const wbvector<TB> &B, wbvector<TC> &C,
    char aflag='N', const TA afac=1., const TC cfac=0.,
    const char i00flag=0
);

template<class T> 
void SVD(
    const wbMatrix<T> &A,
    wbMatrix<T> &U,
    wbMatrix<double> &S,
    wbMatrix<T> &Vt 
);

}; 

#endif

