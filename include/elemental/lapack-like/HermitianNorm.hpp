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

#include "elemental/lapack-like/HermitianSVD.hpp"
#include "elemental/lapack-like/Norm.hpp"

#include "./HermitianNorm/One.hpp"
#include "./HermitianNorm/Infinity.hpp"
#include "./HermitianNorm/Max.hpp"

#include "./HermitianNorm/Nuclear.hpp"
#include "./HermitianNorm/Frobenius.hpp"
#include "./HermitianNorm/Two.hpp"

namespace elem {

template<typename F>
inline typename Base<F>::type
HermitianNorm( UpperOrLower uplo, const Matrix<F>& A, NormType type )
{
#ifndef RELEASE
    PushCallStack("HermitianNorm");
#endif
    typename Base<F>::type norm = 0;
    switch( type )
    {
    case ONE_NORM:
        norm = internal::HermitianOneNorm( uplo, A );
        break;
    case INFINITY_NORM:
        norm = internal::HermitianInfinityNorm( uplo, A );
        break;
    case MAX_NORM:
        norm = internal::HermitianMaxNorm( uplo, A );
        break;
    case NUCLEAR_NORM:
        norm = internal::HermitianNuclearNorm( uplo, A );
        break;
    case FROBENIUS_NORM: 
        norm = internal::HermitianFrobeniusNorm( uplo, A );
        break;
    case TWO_NORM:
        norm = internal::HermitianTwoNorm( uplo, A );
        break;
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

template<typename F>
inline typename Base<F>::type
HermitianNorm( UpperOrLower uplo, const DistMatrix<F>& A, NormType type )
{
#ifndef RELEASE
    PushCallStack("HermitianNorm");
#endif
    typename Base<F>::type norm = 0;
    switch( type )
    {
    case ONE_NORM:
        norm = internal::HermitianOneNorm( uplo, A );
        break;
    case INFINITY_NORM:
        norm = internal::HermitianInfinityNorm( uplo, A );
        break;
    case MAX_NORM:
        norm = internal::HermitianMaxNorm( uplo, A );
        break;
    case NUCLEAR_NORM:
        norm = internal::HermitianNuclearNorm( uplo, A );
        break;
    case FROBENIUS_NORM: 
        norm = internal::HermitianFrobeniusNorm( uplo, A );
        break;
    case TWO_NORM:
        norm = internal::HermitianTwoNorm( uplo, A );
        break;
    }
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

} // namespace elem

#endif // ifndef LAPACK_HERMITIANNORM_HPP
