/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_QR_HPP
#define LAPACK_QR_HPP

#include "elemental/lapack-like/QR/BusingerGolub.hpp"
#include "elemental/lapack-like/QR/Cholesky.hpp"
#include "elemental/lapack-like/QR/Householder.hpp"
#include "elemental/lapack-like/QR/Explicit.hpp"

namespace elem {

// On exit, the upper triangle of A is overwritten by R, and the Householder
// transforms that determine Q are stored below the diagonal of A with an 
// implicit one on the diagonal. 
//
// In the complex case, the column-vector t stores the unit-magnitude complex 
// rotations that map the norms of the implicit Householder vectors to their
// coefficient:  
//                psi_j = 2 tau_j / ( u_j^H u_j ),
// where tau_j is the j'th entry of t and u_j is the j'th unscaled Householder
// reflector.

template<typename Real> 
inline void
QR( Matrix<Real>& A )
{
#ifndef RELEASE
    CallStackEntry entry("QR");
#endif
    qr::Householder( A );
}

template<typename Real> 
inline void
QR( DistMatrix<Real>& A )
{
#ifndef RELEASE
    CallStackEntry entry("QR");
#endif
    qr::Householder( A );
}

template<typename Real> 
inline void
QR( Matrix<Complex<Real> >& A, 
    Matrix<Complex<Real> >& t )
{
#ifndef RELEASE
    CallStackEntry entry("QR");
#endif
    qr::Householder( A, t );
}

template<typename Real> 
inline void
QR( DistMatrix<Complex<Real> >& A, 
    DistMatrix<Complex<Real>,MD,STAR>& t )
{
#ifndef RELEASE
    CallStackEntry entry("QR");
    if( A.Grid() != t.Grid() )
        throw std::logic_error("{A,s} must be distributed over the same grid");
#endif
    qr::Householder( A, t );
}

} // namespace elem

#endif // ifndef LAPACK_QR_HPP
