/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_BLAS_ADJOINT_HPP
#define ELEM_BLAS_ADJOINT_HPP

#include "elemental/blas-like/level1/Transpose.hpp"

namespace elem {

template<typename T>
inline void
Adjoint( const Matrix<T>& A, Matrix<T>& B )
{
#ifndef RELEASE
    CallStackEntry entry("Adjoint");
#endif
    Transpose( A, B, true );
}

template<typename T,Distribution U,Distribution V,
                    Distribution W,Distribution Z>
inline void
Adjoint( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B )
{
#ifndef RELEASE
    CallStackEntry entry("Adjoint");
#endif
    Transpose( A, B, true );
}

} // namespace elem

#endif // ifndef ELEM_BLAS_ADJOINT_HPP
