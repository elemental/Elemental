/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef BLAS_HEMV_HPP
#define BLAS_HEMV_HPP

#include "elemental/blas-like/level2/Symv.hpp"

namespace elem {

template<typename T>
inline void
Hemv
( UpperOrLower uplo,
  T alpha, const Matrix<T>& A, const Matrix<T>& x, T beta, Matrix<T>& y )
{
#ifndef RELEASE
    PushCallStack("Hemv");
#endif
    Symv( uplo, alpha, A, x, beta, y, true );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Hemv
( UpperOrLower uplo,
  T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& x,
  T beta,        DistMatrix<T>& y )
{
#ifndef RELEASE
    PushCallStack("Hemv");
#endif
    Symv( uplo, alpha, A, x, beta, y, true );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem

#endif // ifndef BLAS_HEMV_HPP
