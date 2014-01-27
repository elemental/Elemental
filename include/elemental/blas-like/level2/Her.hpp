/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_HER_HPP
#define ELEM_HER_HPP

#include "./Syr.hpp"

namespace elem {

template<typename T>
inline void
Her( UpperOrLower uplo, T alpha, const Matrix<T>& x, Matrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("Her"))
    Syr( uplo, alpha, x, A, true );
}

template<typename T>
inline void
Her
( UpperOrLower uplo,
  T alpha, const DistMatrix<T>& x,
                 DistMatrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("Her"))
    Syr( uplo, alpha, x, A, true );
}

} // namespace elem

#endif // ifndef ELEM_HER_HPP
