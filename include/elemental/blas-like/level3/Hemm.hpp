/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

#include "./Hemm/LL.hpp"
#include "./Hemm/LU.hpp"
#include "./Hemm/RL.hpp"
#include "./Hemm/RU.hpp"

namespace elem {

template<typename T>
inline void
Hemm
( LeftOrRight side, UpperOrLower uplo,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("Hemm");
#endif
    const char sideChar = LeftOrRightToChar( side );
    const char uploChar = UpperOrLowerToChar( uplo );
    blas::Hemm
    ( sideChar, uploChar, C.Height(), C.Width(),
      alpha, A.LockedBuffer(), A.LDim(),
             B.LockedBuffer(), B.LDim(),
      beta,  C.Buffer(),       C.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Hemm
( LeftOrRight side, UpperOrLower uplo,
  T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("Hemm");
#endif
    if( side == LEFT && uplo == LOWER )
    {
        internal::HemmLL( alpha, A, B, beta, C );
    }
    else if( side == LEFT )
    {
        internal::HemmLU( alpha, A, B, beta, C );
    }
    else if( uplo == LOWER )
    {
        internal::HemmRL( alpha, A, B, beta, C );
    }
    else
    {
        internal::HemmRU( alpha, A, B, beta, C );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
