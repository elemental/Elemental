/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_QUASITRSV_HPP
#define ELEM_QUASITRSV_HPP

#include "./QuasiTrsv/LN.hpp"
#include "./QuasiTrsv/LT.hpp"
#include "./QuasiTrsv/UN.hpp"
#include "./QuasiTrsv/UT.hpp"

namespace elem {

template<typename F>
inline void
QuasiTrsv
( UpperOrLower uplo, Orientation orientation, const Matrix<F>& A, Matrix<F>& x, 
  bool checkIfSingular=false )
{
    DEBUG_ONLY(CallStackEntry cse("QuasiTrsv"))
    if( uplo == LOWER )
    {
        if( orientation == NORMAL )
            internal::QuasiTrsvLN( A, x, checkIfSingular );
        else
            internal::QuasiTrsvLT( orientation, A, x, checkIfSingular );
    }
    else
    {
        if( orientation == NORMAL )
            internal::QuasiTrsvUN( A, x, checkIfSingular );
        else
            internal::QuasiTrsvUT( orientation, A, x, checkIfSingular );
    }
}

template<typename F>
inline void
QuasiTrsv
( UpperOrLower uplo, Orientation orientation, 
  const DistMatrix<F>& A, DistMatrix<F>& x, bool checkIfSingular=false )
{
    DEBUG_ONLY(CallStackEntry cse("QuasiTrsv"))
    if( uplo == LOWER )
    {
        if( orientation == NORMAL )
            internal::QuasiTrsvLN( A, x, checkIfSingular );
        else
            internal::QuasiTrsvLT( orientation, A, x, checkIfSingular );
    }
    else
    {
        if( orientation == NORMAL )
            internal::QuasiTrsvUN( A, x, checkIfSingular );
        else
            internal::QuasiTrsvUT( orientation, A, x, checkIfSingular );
    }
}

} // namespace elem

#endif // ifndef ELEM_QUASITRSV_HPP
