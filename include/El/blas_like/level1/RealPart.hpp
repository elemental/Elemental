/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_REALPART_HPP
#define EL_BLAS_REALPART_HPP

namespace El {

template<typename T>
void RealPart( const Matrix<T>& A, Matrix<Base<T>>& AReal )
{
    EL_DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    AReal.Resize( m, n );
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            AReal(i,j) = RealPart(A(i,j));
}

template<typename T>
void RealPart
( const AbstractDistMatrix<T>& A, AbstractDistMatrix<Base<T>>& AReal )
{
    auto realPart = []( const T& alpha ) { return RealPart(alpha); };
    EntrywiseMap( A, AReal, MakeFunction(realPart) );
}

#ifdef EL_INSTANTIATE_BLAS_LEVEL1
# define EL_EXTERN
#else
# define EL_EXTERN extern
#endif

#define PROTO(T) \
  EL_EXTERN template void RealPart \
  ( const Matrix<T>& A, Matrix<Base<T>>& AReal ); \
  EL_EXTERN template void RealPart \
  ( const AbstractDistMatrix<T>& A, AbstractDistMatrix<Base<T>>& AReal );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

#undef EL_EXTERN

} // namespace El

#endif // ifndef EL_BLAS_REALPART_HPP
