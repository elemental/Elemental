/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>
#include <El/blas_like/level2.hpp>

#include "./QuasiTrsv/LN.hpp"
#include "./QuasiTrsv/LT.hpp"
#include "./QuasiTrsv/UN.hpp"
#include "./QuasiTrsv/UT.hpp"

namespace El {

template<typename F>
void QuasiTrsv
( UpperOrLower uplo, Orientation orientation, const Matrix<F>& A, Matrix<F>& x, 
  bool checkIfSingular )
{
    DEBUG_CSE
    if( uplo == LOWER )
    {
        if( orientation == NORMAL )
            quasitrsv::LN( A, x, checkIfSingular );
        else
            quasitrsv::LT( orientation, A, x, checkIfSingular );
    }
    else
    {
        if( orientation == NORMAL )
            quasitrsv::UN( A, x, checkIfSingular );
        else
            quasitrsv::UT( orientation, A, x, checkIfSingular );
    }
}

template<typename F>
void QuasiTrsv
( UpperOrLower uplo, Orientation orientation, 
  const AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& x, 
  bool checkIfSingular )
{
    DEBUG_CSE
    if( uplo == LOWER )
    {
        if( orientation == NORMAL )
            quasitrsv::LN( A, x, checkIfSingular );
        else
            quasitrsv::LT( orientation, A, x, checkIfSingular );
    }
    else
    {
        if( orientation == NORMAL )
            quasitrsv::UN( A, x, checkIfSingular );
        else
            quasitrsv::UT( orientation, A, x, checkIfSingular );
    }
}

#define PROTO(F) \
  template void QuasiTrsv \
  ( UpperOrLower uplo, Orientation orientation, \
    const Matrix<F>& A, Matrix<F>& x, bool checkIfSingular ); \
  template void QuasiTrsv \
  ( UpperOrLower uplo, Orientation orientation, \
    const AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& x, \
    bool checkIfSingular );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
