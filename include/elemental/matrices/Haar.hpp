/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_HAAR_HPP
#define ELEM_HAAR_HPP

#include ELEM_QR_INC
#include ELEM_GAUSSIAN_INC

namespace elem {

template<typename F>
inline void
Haar( Matrix<F>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Haar"))
    // TODO: Replace this with a quadratic scheme similar to Stewart's, which
    //       essentially generates random Householder reflectors
    Gaussian( A, n, n );
    qr::Explicit( A );
}

template<typename F>
inline void
ImplicitHaar( Matrix<F>& A, Matrix<F>& t, Matrix<Base<F>>& d, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("ImplicitHaar"))
    // TODO: Replace this with a quadratic scheme similar to Stewart's, which
    //       essentially generates random Householder reflectors
    Gaussian( A, n, n );
    QR( A, t, d );
}

template<typename F>
inline void
Haar( DistMatrix<F>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Haar"))
    // TODO: Replace this with a quadratic scheme similar to Stewart's, which
    //       essentially generates random Householder reflectors
    Gaussian( A, n, n );
    qr::Explicit( A );
}

template<typename F>
inline void
ImplicitHaar
( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& t, DistMatrix<Base<F>,MD,STAR>& d, 
  Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Haar"))
    // TODO: Replace this with a quadratic scheme similar to Stewart's, which
    //       essentially generates random Householder reflectors
    Gaussian( A, n, n );
    QR( A, t, d );
}

} // namespace elem

#endif // ifndef ELEM_HAAR_HPP
