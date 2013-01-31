/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_HERMITIANNORM_NUCLEAR_HPP
#define LAPACK_HERMITIANNORM_NUCLEAR_HPP

#include "elemental/lapack-like/HermitianNorm/Schatten.hpp"

namespace elem {

template<typename F> 
inline typename Base<F>::type
HermitianNuclearNorm( UpperOrLower uplo, const Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("HermitianNuclearNorm");
#endif
    typedef typename Base<F>::type R;
    const R norm = HermitianSchattenNorm( uplo, A, R(1) ); 
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

template<typename F,Distribution U,Distribution V> 
inline typename Base<F>::type
HermitianNuclearNorm( UpperOrLower uplo, const DistMatrix<F,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("HermitianNuclearNorm");
#endif
    typedef typename Base<F>::type R;
    const R norm = HermitianSchattenNorm( uplo, A, R(1) ); 
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

} // namespace elem

#endif // ifndef LAPACK_HERMITIANNORM_NUCLEAR_HPP
