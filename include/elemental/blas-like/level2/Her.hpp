/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef BLAS_HER_HPP
#define BLAS_HER_HPP

#include "elemental/blas-like/level2/Syr.hpp"

namespace elem {

template<typename T>
inline void
Her( UpperOrLower uplo, T alpha, const Matrix<T>& x, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("Her");
#endif
    Syr( uplo, alpha, x, A, true );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
Her
( UpperOrLower uplo,
  T alpha, const DistMatrix<T>& x,
                 DistMatrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("Her");
#endif
    Syr( uplo, alpha, x, A, true );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem

#endif // ifndef BLAS_HER_HPP
