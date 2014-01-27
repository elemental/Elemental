/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_INVERSE_HPD_HPP
#define ELEM_INVERSE_HPD_HPP

#include "./HPD/CholeskyLVar2.hpp"
#include "./HPD/CholeskyUVar2.hpp"

namespace elem {

template<typename F>
inline void
HPDInverse( UpperOrLower uplo, Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("HPDInverse"))
    if( uplo == LOWER )
        hpd_inv::CholeskyLVar2( A );
    else
        hpd_inv::CholeskyUVar2( A );
}

template<typename F>
inline void
HPDInverse( UpperOrLower uplo, DistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("HPDInverse"))
    if( uplo == LOWER )
        hpd_inv::CholeskyLVar2( A );
    else
        hpd_inv::CholeskyUVar2( A );
}

template<typename F>
inline void
LocalHPDInverse( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("LocalHPDInverse"))
    HPDInverse( uplo, A.Matrix() );
}

} // namespace elem

#endif // ifndef ELEM_INVERSE_HPD_HPP
