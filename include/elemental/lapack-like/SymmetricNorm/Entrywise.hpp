/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_SYMMETRICNORM_ENTRYWISE_HPP
#define LAPACK_SYMMETRICNORM_ENTRYWISE_HPP

#include "elemental/lapack-like/HermitianNorm/Entrywise.hpp"

namespace elem {

template<typename F> 
inline typename Base<F>::type
SymmetricEntrywiseNorm
( UpperOrLower uplo, const Matrix<F>& A, typename Base<F>::type p )
{
#ifndef RELEASE
    PushCallStack("SymmetricEntrywiseNorm");
#endif
    typedef typename Base<F>::type R;
    const R norm = HermitianEntrywiseNorm( uplo, A, p );
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

template<typename F,Distribution U,Distribution V> 
inline typename Base<F>::type
SymmetricEntrywiseNorm
( UpperOrLower uplo, const DistMatrix<F,U,V>& A, typename Base<F>::type p )
{
#ifndef RELEASE
    PushCallStack("SymmetricEntrywiseNorm");
#endif
    typedef typename Base<F>::type R;
    const R norm = HermitianEntrywiseNorm( uplo, A, p );
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

} // namespace elem

#endif // ifndef LAPACK_SYMMETRICNORM_ENTRYWISE_HPP
