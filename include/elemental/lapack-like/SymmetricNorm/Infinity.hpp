/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_SYMMETRICNORM_INFINITY_HPP
#define LAPACK_SYMMETRICNORM_INFINITY_HPP

#include "elemental/lapack-like/HermitianNorm/Infinity.hpp"

namespace elem {

template<typename F>
inline typename Base<F>::type
SymmetricInfinityNorm( UpperOrLower uplo, const Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("SymmetricInfinityNorm");
#endif
    typedef typename Base<F>::type R;
    const R norm = HermitianInfinityNorm( uplo, A );
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

template<typename F>
inline typename Base<F>::type
SymmetricInfinityNorm( UpperOrLower uplo, const DistMatrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("SymmetricInfinityNorm");
#endif
    typedef typename Base<F>::type R;
    const R norm = HermitianInfinityNorm( uplo, A );
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

} // namespace elem

#endif // ifndef LAPACK_SYMMETRICNORM_INFINITY_HPP
