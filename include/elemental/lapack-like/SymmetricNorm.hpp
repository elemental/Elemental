/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_SYMMETRICNORM_HPP
#define LAPACK_SYMMETRICNORM_HPP

#include "elemental/lapack-like/SymmetricNorm/Frobenius.hpp"
#include "elemental/lapack-like/SymmetricNorm/Infinity.hpp"
#include "elemental/lapack-like/SymmetricNorm/Max.hpp"
#include "elemental/lapack-like/SymmetricNorm/Nuclear.hpp"
#include "elemental/lapack-like/SymmetricNorm/One.hpp"
#include "elemental/lapack-like/SymmetricNorm/Two.hpp"

namespace elem {

template<typename F>
inline typename Base<F>::type
SymmetricNorm( UpperOrLower uplo, const Matrix<F>& A, NormType type )
{
#ifndef RELEASE
    PushCallStack("SymmetricNorm");
#endif
    typename Base<F>::type norm = 0;
    switch( type )
    {
    case ONE_NORM:
        norm = SymmetricOneNorm( uplo, A );
        break;
    case INFINITY_NORM:
        norm = SymmetricInfinityNorm( uplo, A );
        break;
    case MAX_NORM:
        norm = SymmetricMaxNorm( uplo, A );
        break;
    case NUCLEAR_NORM:
        norm = SymmetricNuclearNorm( uplo, A );
        break;
    case FROBENIUS_NORM:
        norm = SymmetricFrobeniusNorm( uplo, A );
        break;
    case TWO_NORM:
        norm = SymmetricTwoNorm( uplo, A );
        break;
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

template<typename F>
inline typename Base<F>::type
SymmetricNorm( UpperOrLower uplo, const DistMatrix<F>& A, NormType type )
{
#ifndef RELEASE
    PushCallStack("SymmetricNorm");
#endif
    typename Base<F>::type norm = 0;
    switch( type )
    {
    case ONE_NORM:
        norm = SymmetricOneNorm( uplo, A );
        break;
    case INFINITY_NORM:
        norm = SymmetricInfinityNorm( uplo, A );
        break;
    case MAX_NORM:
        norm = SymmetricMaxNorm( uplo, A );
        break;
    case NUCLEAR_NORM:
        norm = SymmetricNuclearNorm( uplo, A );
        break;
    case FROBENIUS_NORM:
        norm = SymmetricFrobeniusNorm( uplo, A );
        break;
    case TWO_NORM:
        norm = SymmetricTwoNorm( uplo, A );
        break;
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

} // namespace elem

#endif // ifndef LAPACK_SYMMETRICNORM_HPP
