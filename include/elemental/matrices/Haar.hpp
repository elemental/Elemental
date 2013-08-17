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

template<typename T>
inline void
Haar( Matrix<T>& A, Int n )
{
#ifndef RELEASE
    CallStackEntry entry("Haar");
#endif
    // TODO: Replace this with a quadratic scheme similar to Stewart's, which
    //       essentially generates random Householder reflectors
    Gaussian( A, n, n );
    qr::Explicit( A );
}

template<typename T>
inline void
Haar( DistMatrix<T>& A, Int n )
{
#ifndef RELEASE
    CallStackEntry entry("Haar");
#endif
    // TODO: Replace this with a quadratic scheme similar to Stewart's, which
    //       essentially generates random Householder reflectors
    Gaussian( A, n, n );
    qr::Explicit( A );
}

} // namespace elem

#endif // ifndef ELEM_MATRICES_HAAR_HPP
