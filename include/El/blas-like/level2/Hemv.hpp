/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_HEMV_HPP
#define EL_HEMV_HPP

#include "./Symv.hpp"

namespace El {

template<typename T>
inline void
Hemv
( UpperOrLower uplo,
  T alpha, const Matrix<T>& A, const Matrix<T>& x, T beta, Matrix<T>& y )
{
    DEBUG_ONLY(CallStackEntry cse("Hemv"))
    Symv( uplo, alpha, A, x, beta, y, true );
}

template<typename T>
inline void
Hemv
( UpperOrLower uplo,
  T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& x,
  T beta,        DistMatrix<T>& y )
{
    DEBUG_ONLY(CallStackEntry cse("Hemv"))
    Symv( uplo, alpha, A, x, beta, y, true );
}

} // namespace El

#endif // ifndef EL_HEMV_HPP
