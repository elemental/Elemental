/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_HERMITIANNORM_HPP
#define LAPACK_HERMITIANNORM_HPP

#include "elemental/lapack-like/HermitianNorm/EntrywiseOne.hpp"
#include "elemental/lapack-like/HermitianNorm/Frobenius.hpp"
#include "elemental/lapack-like/HermitianNorm/Infinity.hpp"
#include "elemental/lapack-like/HermitianNorm/One.hpp"
#include "elemental/lapack-like/HermitianNorm/Max.hpp"

#include "elemental/lapack-like/HermitianNorm/Nuclear.hpp"
#include "elemental/lapack-like/HermitianNorm/Two.hpp"

namespace elem {

template<typename F>
inline typename Base<F>::type
HermitianNorm
( UpperOrLower uplo, const Matrix<F>& A, NormType type=FROBENIUS_NORM )
{
#ifndef RELEASE
    PushCallStack("HermitianNorm");
#endif
    typename Base<F>::type norm = 0;
    switch( type )
    {
    // The following norms are rather cheap to compute
    case FROBENIUS_NORM: 
        norm = HermitianFrobeniusNorm( uplo, A );
        break;
    case ENTRYWISE_ONE_NORM:
        norm = HermitianEntrywiseOneNorm( uplo, A );
        break;
    case INFINITY_NORM:
        norm = HermitianInfinityNorm( uplo, A );
        break;
    case MAX_NORM:
        norm = HermitianMaxNorm( uplo, A );
        break;
    case ONE_NORM:
        norm = HermitianOneNorm( uplo, A );
        break;
    // The following norms make use of an SVD
    case NUCLEAR_NORM:
        norm = HermitianNuclearNorm( uplo, A );
        break;
    case TWO_NORM:
        norm = HermitianTwoNorm( uplo, A );
        break;
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

template<typename F>
inline typename Base<F>::type
HermitianNorm
( UpperOrLower uplo, const DistMatrix<F>& A, NormType type=FROBENIUS_NORM )
{
#ifndef RELEASE
    PushCallStack("HermitianNorm");
#endif
    typename Base<F>::type norm = 0;
    switch( type )
    {
    // The following norms are rather cheap to compute
    case FROBENIUS_NORM:
        norm = HermitianFrobeniusNorm( uplo, A );
        break;
    case ENTRYWISE_ONE_NORM:
        norm = HermitianEntrywiseOneNorm( uplo, A );
        break;
    case INFINITY_NORM:
        norm = HermitianInfinityNorm( uplo, A );
        break;
    case MAX_NORM:
        norm = HermitianMaxNorm( uplo, A );
        break;
    case ONE_NORM:
        norm = HermitianOneNorm( uplo, A );
        break;
    // The following norms make use of an SVD
    case NUCLEAR_NORM:
        norm = HermitianNuclearNorm( uplo, A );
        break;
    case TWO_NORM:
        norm = HermitianTwoNorm( uplo, A );
        break; 
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

} // namespace elem

#endif // ifndef LAPACK_HERMITIANNORM_HPP
