/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_FILLDIAGONAL_HPP
#define EL_BLAS_FILLDIAGONAL_HPP

namespace El {

template<typename T>
void FillDiagonal( Matrix<T>& A, T alpha, Int offset )
{
    DEBUG_CSE
    const Int height = A.Height();
    const Int width = A.Width();
    for( Int j=0; j<width; ++j )
    {
        const Int i = j-offset;
        if( i >= 0 && i < height )
            A(i,j) = alpha;
    }
}

template<typename T>
void FillDiagonal( AbstractDistMatrix<T>& A, T alpha, Int offset )
{
    DEBUG_CSE
    const Int height = A.Height();
    const Int localWidth = A.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = A.GlobalCol(jLoc);
        const Int i = j-offset;
        if( i >= 0 && i < height )
            A.Set( i, j, alpha );
    }
}

#ifdef EL_INSTANTIATE_BLAS_LEVEL1
# define EL_EXTERN
#else
# define EL_EXTERN extern
#endif

#define PROTO(T) \
  EL_EXTERN template void FillDiagonal \
  ( Matrix<T>& A, T alpha, Int offset ); \
  EL_EXTERN template void FillDiagonal \
  ( AbstractDistMatrix<T>& A, T alpha, Int offset );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

#undef EL_EXTERN

} // namespace El

#endif // ifndef EL_BLAS_FILLDIAGONAL_HPP
