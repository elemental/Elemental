/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef BLAS_SYR2K_HPP
#define BLAS_SYR2K_HPP

#include "./Syr2k/LN.hpp"
#include "./Syr2k/LT.hpp"
#include "./Syr2k/UN.hpp"
#include "./Syr2k/UT.hpp"

namespace elem {

template<typename T>
inline void
Syr2k
( UpperOrLower uplo, Orientation orientation,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C,
  bool conjugate )
{
#ifndef RELEASE
    PushCallStack("Syr2k");
    if( orientation == NORMAL )
    {
        if( A.Height() != C.Height() || A.Height() != C.Width() ||
            B.Height() != C.Height() ||B.Height() != C.Width()    )
            throw std::logic_error("Nonconformal Syr2k");
    }
    else 
    {
        if( A.Width() != C.Height() || A.Width() != C.Width() ||
            B.Width() != C.Height() || B.Width() != C.Width()   )
            throw std::logic_error("Nonconformal Syr2k");
    }
#endif
    const char uploChar = UpperOrLowerToChar( uplo );
    const char transChar = OrientationToChar( orientation );
    const int k = ( orientation == NORMAL ? A.Width() : A.Height() );
    if( conjugate )
    {
        blas::Her2k
        ( uploChar, transChar, C.Height(), k,
          alpha, A.LockedBuffer(), A.LDim(),
                 B.LockedBuffer(), B.LDim(),
          beta,  C.Buffer(),       C.LDim() );
    }
    else
    {
        blas::Syr2k
        ( uploChar, transChar, C.Height(), k,
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
Syr2k
( UpperOrLower uplo, Orientation orientation,
  T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C,
  bool conjugate )
{
#ifndef RELEASE
    PushCallStack("Syr2k");
#endif
    if( uplo == LOWER && orientation == NORMAL )
        internal::Syr2kLN( alpha, A, B, beta, C, conjugate );
    else if( uplo == LOWER )
        internal::Syr2kLT( alpha, A, B, beta, C, conjugate );
    else if( orientation == NORMAL )
        internal::Syr2kUN( alpha, A, B, beta, C, conjugate );
    else
        internal::Syr2kUT( alpha, A, B, beta, C, conjugate );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem

#endif // ifndef BLAS_SYR2K_HPP
