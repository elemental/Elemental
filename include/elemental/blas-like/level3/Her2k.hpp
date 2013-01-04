/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

#include "./Her2k/LC.hpp"
#include "./Her2k/LN.hpp"
#include "./Her2k/UC.hpp"
#include "./Her2k/UN.hpp"

namespace elem {

template<typename T>
inline void
Her2k
( UpperOrLower uplo, Orientation orientation,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("Her2k");
    if( orientation == NORMAL )
    {
        if( A.Height() != C.Height() || A.Height() != C.Width() ||
            B.Height() != C.Height() ||B.Height() != C.Width() )
            throw std::logic_error("Nonconformal Her2k");
    }
    else if( orientation == ADJOINT )
    {
        if( A.Width() != C.Height() || A.Width() != C.Width() ||
            B.Width() != C.Height() || B.Width() != C.Width() )
            throw std::logic_error("Nonconformal Her2k");
    }
    else
        throw std::logic_error
        ("Her2k only accepts NORMAL and ADJOINT options");
#endif
    const char uploChar = UpperOrLowerToChar( uplo );
    const char transChar = OrientationToChar( orientation );
    const int k = ( orientation == NORMAL ? A.Width() : A.Height() );
    blas::Her2k
    ( uploChar, transChar, C.Height(), k,
      alpha, A.LockedBuffer(), A.LDim(),
             B.LockedBuffer(), B.LDim(),
      beta,  C.Buffer(),       C.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Her2k
( UpperOrLower uplo, Orientation orientation,
  T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("Her2k");
    if( orientation == TRANSPOSE )
        throw std::logic_error("Her2k accepts NORMAL and ADJOINT options");
#endif
    if( uplo == LOWER && orientation == NORMAL )
        internal::Her2kLN( alpha, A, B, beta, C );
    else if( uplo == LOWER )
        internal::Her2kLC( alpha, A, B, beta, C );
    else if( orientation == NORMAL )
        internal::Her2kUN( alpha, A, B, beta, C );
    else
        internal::Her2kUC( alpha, A, B, beta, C );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
