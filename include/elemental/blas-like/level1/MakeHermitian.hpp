/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef BLAS_MAKEHERMITIAN_HPP
#define BLAS_MAKEHERMITIAN_HPP

#include "elemental/blas-like/level1/MakeSymmetric.hpp"

namespace elem {

template<typename T>
inline void
MakeHermitian( UpperOrLower uplo, Matrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("MakeHermitian");
#endif
    MakeSymmetric( uplo, A, true );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
MakeHermitian( UpperOrLower uplo, DistMatrix<T>& A )
{
#ifndef RELEASE
    PushCallStack("MakeHermitian");
#endif
    MakeSymmetric( uplo, A, true );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem

#endif // ifndef BLAS_MAKEHERMITIAN_HPP
