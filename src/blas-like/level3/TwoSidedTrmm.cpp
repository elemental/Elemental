/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

#include "./TwoSidedTrmm/Unblocked.hpp"
#include "./TwoSidedTrmm/LVar4.hpp"
#include "./TwoSidedTrmm/UVar4.hpp"

namespace El {

template<typename T> 
void TwoSidedTrmm
( UpperOrLower uplo, UnitOrNonUnit diag, Matrix<T>& A, const Matrix<T>& B )
{
    DEBUG_ONLY(CallStackEntry cse("TwoSidedTrmm"))
    if( uplo == LOWER )
        twotrmm::LVar4( diag, A, B );
    else
        twotrmm::UVar4( diag, A, B );
}

template<typename T> 
void TwoSidedTrmm
( UpperOrLower uplo, UnitOrNonUnit diag, 
  DistMatrix<T>& A, const DistMatrix<T>& B )
{
    DEBUG_ONLY(CallStackEntry cse("TwoSidedTrmm"))
    if( uplo == LOWER )
        twotrmm::LVar4( diag, A, B );
    else
        twotrmm::UVar4( diag, A, B );
}

#define PROTO(T) \
  template void TwoSidedTrmm \
  ( UpperOrLower uplo, UnitOrNonUnit diag, \
    Matrix<T>& A, const Matrix<T>& B ); \
  template void TwoSidedTrmm \
  ( UpperOrLower uplo, UnitOrNonUnit diag, \
    DistMatrix<T>& A, const DistMatrix<T>& B );

PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
