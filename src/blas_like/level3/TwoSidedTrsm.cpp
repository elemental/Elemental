/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#include "./TwoSidedTrsm/Unblocked.hpp"
#include "./TwoSidedTrsm/LVar4.hpp"
#include "./TwoSidedTrsm/UVar4.hpp"

namespace El {

template<typename F> 
void TwoSidedTrsm
( UpperOrLower uplo, UnitOrNonUnit diag, Matrix<F>& A, const Matrix<F>& B )
{
    DEBUG_ONLY(CallStackEntry cse("TwoSidedTrsm"))
    if( uplo == LOWER )
        twotrsm::LVar4( diag, A, B );
    else
        twotrsm::UVar4( diag, A, B );
}

template<typename F> 
void TwoSidedTrsm
( UpperOrLower uplo, UnitOrNonUnit diag, 
  AbstractDistMatrix<F>& A, const AbstractDistMatrix<F>& B )
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
    AbstractDistMatrix<T>& A, const AbstractDistMatrix<T>& B );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
