/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#include "./Trtrmm/Unblocked.hpp"
#include "./Trtrmm/LVar1.hpp"
#include "./Trtrmm/UVar1.hpp"

namespace El {

template<typename T>
void Trtrmm( UpperOrLower uplo, Matrix<T>& A, bool conjugate )
{
    DEBUG_ONLY(
        CallStackEntry cse("Trtrmm");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
    )
    if( uplo == LOWER )
        trtrmm::LVar1( A, conjugate );
    else
        trtrmm::UVar1( A, conjugate );
}

template<typename T>
void Trtrmm( UpperOrLower uplo, AbstractDistMatrix<T>& A, bool conjugate )
{
    DEBUG_ONLY(
        CallStackEntry cse("Trtrmm");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
    )
    if( uplo == LOWER )
        trtrmm::LVar1( A, conjugate );
    else
        trtrmm::UVar1( A, conjugate );
}

#define PROTO(T) \
  template void Trtrmm( UpperOrLower uplo, Matrix<T>& A, bool conjugate ); \
  template void Trtrmm \
  ( UpperOrLower uplo, AbstractDistMatrix<T>& A, bool conjugate );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
