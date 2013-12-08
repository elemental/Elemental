/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_BLAS_HEMM_HPP
#define ELEM_BLAS_HEMM_HPP

#include "elemental/blas-like/level3/Symm.hpp"

namespace elem {

template<typename T>
inline void
Hemm
( LeftOrRight side, UpperOrLower uplo,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C )
{
#ifndef RELEASE
    CallStackEntry cse("Hemm");
#endif
    Symm( side, uplo, alpha, A, B, beta, C, true );
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
    CallStackEntry cse("Hemm");
#endif
    Symm( side, uplo, alpha, A, B, beta, C, true );
}

} // namespace elem

#endif // ifndef ELEM_BLAS_HEMM_HPP
