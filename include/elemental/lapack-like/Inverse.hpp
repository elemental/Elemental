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
    PushCallStack("Inverse");
#endif
    inverse::LUPartialPiv( A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
HPDInverse( UpperOrLower uplo, Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("HPDInverse");
#endif
    if( uplo == LOWER )
        hpd_inverse::CholeskyLVar2( A );
    else
        hpd_inverse::CholeskyUVar2( A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> 
inline void
Inverse( DistMatrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("Inverse");
#endif
    inverse::LUPartialPiv( A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
HPDInverse( UpperOrLower uplo, DistMatrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("HPDInverse");
#endif
    if( uplo == LOWER )
        hpd_inverse::CholeskyLVar2( A );
    else
        hpd_inverse::CholeskyUVar2( A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
LocalInverse( DistMatrix<F,STAR,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("LocalInverse");
#endif
    Inverse( A.Matrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
LocalHPDInverse( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A )
{
#ifndef RELEASE
    PushCallStack("LocalHPDInverse");
#endif
    HPDInverse( uplo, A.Matrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem

#endif // ifndef LAPACK_INVERSE_HPP
