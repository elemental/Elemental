/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MATRICES_HAAR_HPP
#define ELEM_MATRICES_HAAR_HPP

#include "elemental/lapack-like/QR.hpp"
#include "elemental/matrices/Gaussian.hpp"

namespace elem {

template<typename F>
inline void
Haar( Matrix<F>& A, Int n )
{
#ifndef RELEASE
    CallStackEntry cse("Haar");
#endif
    // TODO: Replace this with a quadratic scheme similar to Stewart's, which
    //       essentially generates random Householder reflectors
    Gaussian( A, n, n );
    qr::Explicit( A );
}

template<typename F>
inline Matrix<F>
Haar( Int n )
{
    auto A = Gaussian<F>( n, n );
    qr::Explicit( A );
    return A;
}

template<typename F>
inline void
ImplicitHaar( Matrix<F>& A, Matrix<F>& t, Int n )
{
#ifndef RELEASE
    CallStackEntry cse("ImplicitHaar");
#endif
    // TODO: Replace this with a quadratic scheme similar to Stewart's, which
    //       essentially generates random Householder reflectors
    Gaussian( A, n, n );
    QR( A, t );
}

template<typename F>
inline void
Haar( DistMatrix<F>& A, Int n )
{
#ifndef RELEASE
    CallStackEntry cse("Haar");
#endif
    // TODO: Replace this with a quadratic scheme similar to Stewart's, which
    //       essentially generates random Householder reflectors
    Gaussian( A, n, n );
    qr::Explicit( A );
}

template<typename F>
inline DistMatrix<F>
Haar( const Grid& g, Int n )
{
    auto A = Gaussian<F>( g, n, n );
    qr::Explicit( A );
    return A;
}

template<typename F>
inline void
ImplicitHaar( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& t, Int n )
{
#ifndef RELEASE
    CallStackEntry cse("Haar");
#endif
    // TODO: Replace this with a quadratic scheme similar to Stewart's, which
    //       essentially generates random Householder reflectors
    Gaussian( A, n, n );
    QR( A, t );
}

} // namespace elem

#endif // ifndef ELEM_MATRICES_HAAR_HPP
