/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename Real>
void Conjugate( Matrix<Real>& A )
{ }

template<typename Real>
void Conjugate( Matrix<Complex<Real>>& A )
{
    DEBUG_ONLY(CallStackEntry cse("Conjugate (in-place)"))
    const Int m = A.Height();
    const Int n = A.Width();
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            A.Set(i,j,Conj(A.Get(i,j)));
}

template<typename T>
void Conjugate( const Matrix<T>& A, Matrix<T>& B )
{
    DEBUG_ONLY(CallStackEntry cse("Conjugate"))
    const Int m = A.Height();
    const Int n = A.Width();
    B.Resize( m, n );
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            B.Set(i,j,Conj(A.Get(i,j)));
}

template<typename T>
void Conjugate( AbstractDistMatrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("Conjugate (in-place)"))
    Conjugate( A.Matrix() );
}

template<typename T>
void Conjugate( const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B )
{
    DEBUG_ONLY(CallStackEntry cse("Conjugate"))
    Copy( A, B );
    Conjugate( B );
}

#define PROTO(T) \
  template void Conjugate( Matrix<T>& A ); \
  template void Conjugate( const Matrix<T>& A, Matrix<T>& B ); \
  template void Conjugate( AbstractDistMatrix<T>& A ); \
  template void Conjugate \
  ( const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B );

#include "El/macros/Instantiate.h"

} // namespace El
