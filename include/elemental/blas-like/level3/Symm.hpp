/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef BLAS_SYMM_HPP
#define BLAS_SYMM_HPP

#include "./Symm/LL.hpp"
#include "./Symm/LU.hpp"
#include "./Symm/RL.hpp"
#include "./Symm/RU.hpp"

namespace elem {

template<typename T>
inline void
Symm
( LeftOrRight side, UpperOrLower uplo,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C,
  bool conjugate=false )
{
#ifndef RELEASE
    PushCallStack("Symm");
#endif
    const char sideChar = LeftOrRightToChar( side );
    const char uploChar = UpperOrLowerToChar( uplo );
    if( conjugate )
    {
        blas::Hemm
        ( sideChar, uploChar, C.Height(), C.Width(),
          alpha, A.LockedBuffer(), A.LDim(),
                 B.LockedBuffer(), B.LDim(),
          beta,  C.Buffer(),       C.LDim() );
    }
    else
    {
        blas::Symm
        ( sideChar, uploChar, C.Height(), C.Width(),
          alpha, A.LockedBuffer(), A.LDim(),
                 B.LockedBuffer(), B.LDim(),
          beta,  C.Buffer(),       C.LDim() );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Symm
( LeftOrRight side, UpperOrLower uplo,
  T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C, 
  bool conjugate )
{
#ifndef RELEASE
    PushCallStack("Symm");
#endif
    if( side == LEFT && uplo == LOWER )
    {
        internal::SymmLL( alpha, A, B, beta, C, conjugate );
    }
    else if( side == LEFT )
    {
        internal::SymmLU( alpha, A, B, beta, C, conjugate );
    }
    else if( uplo == LOWER )
    {
        internal::SymmRL( alpha, A, B, beta, C, conjugate );
    }
    else
    {
        internal::SymmRU( alpha, A, B, beta, C, conjugate );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem

#endif // ifndef BLAS_SYMM_HPP
