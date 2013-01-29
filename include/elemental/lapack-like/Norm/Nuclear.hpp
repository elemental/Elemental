/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_NORM_NUCLEAR_HPP
#define LAPACK_NORM_NUCLEAR_HPP

#include "elemental/lapack-like/Norm/One.hpp"
#include "elemental/lapack-like/SVD.hpp"

namespace elem {

template<typename F> 
inline typename Base<F>::type
NuclearNorm( const Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("NuclearNorm");
#endif
    typedef typename Base<F>::type R;
    Matrix<F> B( A );
    Matrix<R> s;
    SingularValues( B, s );
    const R norm = OneNorm( s );
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

template<typename F,Distribution U,Distribution V> 
inline typename Base<F>::type
NuclearNorm( const DistMatrix<F,U,V>& A )
{
#ifndef RELEASE
    PushCallStack("NuclearNorm");
#endif
    typedef typename Base<F>::type R;
    DistMatrix<F> B( A );
    DistMatrix<R,VR,STAR> s( A.Grid() );
    SingularValues( B, s );
    const R norm = OneNorm( s );
#ifndef RELEASE
    PopCallStack();
#endif
    return norm;
}

} // namespace elem

#endif // ifndef LAPACK_NORM_NUCLEAR_HPP
