/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_SYMMETRICNORM_MAX_HPP
#define LAPACK_SYMMETRICNORM_MAX_HPP

#include "elemental/lapack-like/HermitianNorm/Max.hpp"

namespace elem {

template<typename F>
inline typename Base<F>::type
SymmetricMaxNorm( UpperOrLower uplo, const Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("SymmetricMaxNorm");
#endif
    typedef typename Base<F>::type R;
    const R norm = HermitianMaxNorm( uplo, A );
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

template<typename F>
inline typename Base<F>::type
SymmetricMaxNorm( UpperOrLower uplo, const DistMatrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("SymmetricMaxNorm");
#endif
    typedef typename Base<F>::type R;
    const R norm = HermitianMaxNorm( uplo, A );
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

} // namespace elem

#endif // ifndef LAPACK_SYMMETRICNORM_MAX_HPP
