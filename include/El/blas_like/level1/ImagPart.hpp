/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_IMAGPART_HPP
#define EL_BLAS_IMAGPART_HPP

namespace El {

template<typename T>
void ImagPart( const Matrix<T>& A, Matrix<Base<T>>& AImag )
{
    DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    AImag.Resize( m, n );
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            AImag(i,j) = ImagPart(A(i,j));
}

template<typename T>
void ImagPart
( const ElementalMatrix<T>& A, ElementalMatrix<Base<T>>& AImag )
{ 
    auto imagPart = []( T alpha ) { return ImagPart(alpha); };
    function<Base<T>(T)> realLambda( imagPart );
    EntrywiseMap( A, AImag, realLambda );
}

template<typename T>
void ImagPart
( const BlockMatrix<T>& A, BlockMatrix<Base<T>>& AImag )
{ 
    auto imagPart = []( T alpha ) { return ImagPart(alpha); };
    function<Base<T>(T)> realLambda( imagPart );
    EntrywiseMap( A, AImag, realLambda );
}

#ifdef EL_INSTANTIATE_BLAS_LEVEL1
# define EL_EXTERN
#else
# define EL_EXTERN extern
#endif

#define PROTO(T) \
  EL_EXTERN template void ImagPart \
  ( const Matrix<T>& A, Matrix<Base<T>>& AImag ); \
  EL_EXTERN template void ImagPart \
  ( const ElementalMatrix<T>& A, ElementalMatrix<Base<T>>& AImag ); \
  EL_EXTERN template void ImagPart \
  ( const BlockMatrix<T>& A, \
          BlockMatrix<Base<T>>& AImag );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

#undef EL_EXTERN

} // namespace El

#endif // ifndef EL_BLAS_IMAGPART_HPP
