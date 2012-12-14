/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

#include "./Syrk/LN.hpp"
#include "./Syrk/LT.hpp"
#include "./Syrk/UN.hpp"
#include "./Syrk/UT.hpp"

namespace elem {

template<typename T>
inline void
Syrk
( UpperOrLower uplo, Orientation orientation,
  T alpha, const Matrix<T>& A, T beta, Matrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("Syrk");
    if( orientation == NORMAL )
    {
        if( A.Height() != C.Height() || A.Height() != C.Width() )
            throw std::logic_error("Nonconformal Syrk");
    }
    else if( orientation == TRANSPOSE )
    {
        if( A.Width() != C.Height() || A.Width() != C.Width() )
            throw std::logic_error("Nonconformal Syrk");
    }
    else
        throw std::logic_error
        ("Syrk only accepts NORMAL and TRANSPOSE options");
#endif
    const char uploChar = UpperOrLowerToChar( uplo );
    const char transChar = OrientationToChar( orientation );
    const int k = ( orientation == NORMAL ? A.Width() : A.Height() );
    blas::Syrk
    ( uploChar, transChar, C.Height(), k,
      alpha, A.LockedBuffer(), A.LDim(),
      beta,  C.Buffer(),       C.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Syrk
( UpperOrLower uplo, 
  Orientation orientation,
  T alpha, const DistMatrix<T,MC,MR>& A,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("Syrk");
    if( orientation == ADJOINT )
        throw std::logic_error("Syrk accepts NORMAL and TRANSPOSE options");
#endif
    if( uplo == LOWER && orientation == NORMAL )
        internal::SyrkLN( alpha, A, beta, C );
    else if( uplo == LOWER )
        internal::SyrkLT( alpha, A, beta, C );
    else if( orientation == NORMAL )
        internal::SyrkUN( alpha, A, beta, C );
    else
        internal::SyrkUT( alpha, A, beta, C );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
