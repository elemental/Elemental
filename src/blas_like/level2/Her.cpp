/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename T>
void Her( UpperOrLower uplo, Base<T> alpha, const Matrix<T>& x, Matrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("Her"))
    Syr( uplo, T(alpha), x, A, true );
}

template<typename T>
void Her
( UpperOrLower uplo, 
  Base<T> alpha, const AbstractDistMatrix<T>& x, 
                       AbstractDistMatrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("Her"))
    Syr( uplo, T(alpha), x, A, true );
}

#define PROTO(T) \
  template void Her \
  ( UpperOrLower uplo, Base<T> alpha, const Matrix<T>& x, Matrix<T>& A ); \
  template void Her \
  ( UpperOrLower uplo, \
    Base<T> alpha, const AbstractDistMatrix<T>& x, \
                         AbstractDistMatrix<T>& A );

// blas::Her not yet supported
#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
