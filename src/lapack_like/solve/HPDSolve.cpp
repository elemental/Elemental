/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename F>
void HPDSolve
( UpperOrLower uplo, Orientation orientation, 
  Matrix<F>& A, Matrix<F>& B )
{
    DEBUG_ONLY(CallStackEntry cse("HPDSolve"))
    Cholesky( uplo, A );
    cholesky::SolveAfter( uplo, orientation, A, B );
}

template<typename F>
void HPDSolve
( UpperOrLower uplo, Orientation orientation, 
  AbstractDistMatrix<F>& APre, AbstractDistMatrix<F>& BPre )
{
    DEBUG_ONLY(CallStackEntry cse("HPDSolve"))

    auto APtr = ReadProxy<F,MC,MR>( &APre );  auto& A = *APtr;
    auto BPtr = WriteProxy<F,MC,MR>( &BPre ); auto& B = *BPtr;

    Cholesky( uplo, A );
    cholesky::SolveAfter( uplo, orientation, A, B );
}

#define PROTO(F) \
  template void HPDSolve \
  ( UpperOrLower uplo, Orientation orientation, \
    Matrix<F>& A, Matrix<F>& B ); \
  template void HPDSolve \
  ( UpperOrLower uplo, Orientation orientation, \
    AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& B );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
