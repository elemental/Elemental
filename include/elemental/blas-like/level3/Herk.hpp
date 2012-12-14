/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

#include "./Herk/LC.hpp"
#include "./Herk/LN.hpp"
#include "./Herk/UC.hpp"
#include "./Herk/UN.hpp"

namespace elem {

template<typename T>
inline void
Herk
( UpperOrLower uplo, Orientation orientation,
  T alpha, const Matrix<T>& A, T beta, Matrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("Herk");
    if( orientation == NORMAL )
    {
        if( A.Height() != C.Height() || A.Height() != C.Width() )
            throw std::logic_error("Nonconformal Herk");
    }
    else if( orientation == ADJOINT )
    {
        if( A.Width() != C.Height() || A.Width() != C.Width() )
            throw std::logic_error("Nonconformal Herk");
    }
    else
        throw std::logic_error("Herk only accepts NORMAL and ADJOINT options.");
#endif
    const char uploChar = UpperOrLowerToChar( uplo );
    const char transChar = OrientationToChar( orientation );
    const int k = ( orientation == NORMAL ? A.Width() : A.Height() );
    blas::Herk
    ( uploChar, transChar, C.Height(), k,
      alpha, A.LockedBuffer(), A.LDim(),
      beta,  C.Buffer(),       C.LDim() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Herk
( UpperOrLower uplo, 
  Orientation orientation,
  T alpha, const DistMatrix<T,MC,MR>& A,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("Herk");
    if( A.Grid() != C.Grid() )
        throw std::logic_error
        ("A and C must be distributed over the same grid");
    if( orientation == TRANSPOSE )
        throw std::logic_error
        ("Herk accepts NORMAL and ADJOINT options");
#endif
    if( uplo == LOWER && orientation == NORMAL )
        internal::HerkLN( alpha, A, beta, C );
    else if( uplo == LOWER )
        internal::HerkLC( alpha, A, beta, C );
    else if( orientation == NORMAL )
        internal::HerkUN( alpha, A, beta, C );
    else
        internal::HerkUC( alpha, A, beta, C );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
