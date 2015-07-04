/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename Ring>
void Hemm
( LeftOrRight side, UpperOrLower uplo,
  Ring alpha, const Matrix<Ring>& A, const Matrix<Ring>& B, 
  Ring beta,        Matrix<Ring>& C )
{
    DEBUG_ONLY(CSE cse("Hemm"))
    Symm( side, uplo, alpha, A, B, beta, C, true );
}

template<typename Ring>
void Hemm
( LeftOrRight side, UpperOrLower uplo,
  Ring alpha, const AbstractDistMatrix<Ring>& A, 
              const AbstractDistMatrix<Ring>& B,
  Ring beta,        AbstractDistMatrix<Ring>& C )
{
    DEBUG_ONLY(CSE cse("Hemm"))
    Symm( side, uplo, alpha, A, B, beta, C, true );
}

#define PROTO(Ring) \
  template void Hemm \
  ( LeftOrRight side, UpperOrLower uplo, \
    Ring alpha, const Matrix<Ring>& A, const Matrix<Ring>& B, \
    Ring beta,        Matrix<Ring>& C ); \
  template void Hemm \
  ( LeftOrRight side, UpperOrLower uplo, \
    Ring alpha, const AbstractDistMatrix<Ring>& A, \
                const AbstractDistMatrix<Ring>& B, \
    Ring beta,        AbstractDistMatrix<Ring>& C );

#define EL_ENABLE_QUAD
#include "El/macros/Instantiate.h"

} // namespace El
