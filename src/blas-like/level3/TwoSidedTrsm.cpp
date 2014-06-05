/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

#include "./TwoSidedTrsm/Unblocked.hpp"
#include "./TwoSidedTrsm/LVar4.hpp"
#include "./TwoSidedTrsm/UVar4.hpp"

namespace El {

template<typename F> 
inline void
TwoSidedTrsm
( UpperOrLower uplo, UnitOrNonUnit diag, Matrix<F>& A, const Matrix<F>& B )
{
    DEBUG_ONLY(CallStackEntry cse("TwoSidedTrsm"))
    if( uplo == LOWER )
        twotrsm::LVar4( diag, A, B );
    else
        twotrsm::UVar4( diag, A, B );
}

template<typename F> 
inline void
TwoSidedTrsm
( UpperOrLower uplo, UnitOrNonUnit diag, 
  DistMatrix<F>& A, const DistMatrix<F>& B )
{
    DEBUG_ONLY(CallStackEntry cse("TwoSidedTrsm"))
    if( uplo == LOWER )
        twotrsm::LVar4( diag, A, B );
    else
        twotrsm::UVar4( diag, A, B );
}

#define PROTO(T) \
  template void TwoSidedTrsm \
  ( UpperOrLower uplo, UnitOrNonUnit diag, \
    Matrix<T>& A, const Matrix<T>& B ); \
  template void TwoSidedTrsm \
  ( UpperOrLower uplo, UnitOrNonUnit diag, \
    DistMatrix<T>& A, const DistMatrix<T>& B );

PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
