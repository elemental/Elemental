/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_INVERSE_HPP
#define LAPACK_INVERSE_HPP

#include "elemental/lapack-like/Inverse/CholeskyLVar2.hpp"
#include "elemental/lapack-like/Inverse/CholeskyUVar2.hpp"
#include "elemental/lapack-like/Inverse/LUPartialPiv.hpp"

namespace elem {

template<typename F> 
inline void
Inverse( Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("Inverse");
#endif
    inverse::LUPartialPiv( A );
}

template<typename F>
inline void
HPDInverse( UpperOrLower uplo, Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("HPDInverse");
#endif
    if( uplo == LOWER )
        hpd_inverse::CholeskyLVar2( A );
    else
        hpd_inverse::CholeskyUVar2( A );
}

template<typename F> 
inline void
Inverse( DistMatrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("Inverse");
#endif
    inverse::LUPartialPiv( A );
}

template<typename F>
inline void
HPDInverse( UpperOrLower uplo, DistMatrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("HPDInverse");
#endif
    if( uplo == LOWER )
        hpd_inverse::CholeskyLVar2( A );
    else
        hpd_inverse::CholeskyUVar2( A );
}

template<typename F>
inline void
LocalInverse( DistMatrix<F,STAR,STAR>& A )
{
#ifndef RELEASE
    CallStackEntry entry("LocalInverse");
#endif
    Inverse( A.Matrix() );
}

template<typename F>
inline void
LocalHPDInverse( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A )
{
#ifndef RELEASE
    CallStackEntry entry("LocalHPDInverse");
#endif
    HPDInverse( uplo, A.Matrix() );
}

} // namespace elem

#endif // ifndef LAPACK_INVERSE_HPP
