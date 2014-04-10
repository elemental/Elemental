/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_SCHUR_HPP
#define ELEM_SCHUR_HPP

namespace elem {
namespace schur {

template<typename Real>
inline void
CheckQuasiTriangular( const Matrix<Real>& U )
{
    DEBUG_ONLY(CallStackEntry cse("CheckQuasiTriangular")) 
    const Int n = U.Height();
    if( n < 3 )
        return;
    Real thisSub, nextSub=U.Get(1,0);
    for( Int j=0; j<n-2; ++j )
    {
        thisSub = nextSub;
        nextSub = U.Get(j+2,j+1);
        if( thisSub != Real(0) && nextSub != Real(0) )
            LogicError
            ("Quasi-triangular assumption broken at j=",j,
             ": subdiagonals were ",thisSub," and ",nextSub);
    }
}

template<typename Real>
inline void
CheckQuasiTriangular( const DistMatrix<Real>& U )
{
    DEBUG_ONLY(CallStackEntry cse("CheckQuasiTriangular")) 
    const Int n = U.Height();
    if( n < 3 )
        return;
    auto uSub = U.GetDiagonal( -1 );
    DistMatrix<Real,STAR,STAR> uSub_STAR_STAR( uSub );
    Real thisSub, nextSub=uSub_STAR_STAR.Get(0,0);
    for( Int j=0; j<n-2; ++j )
    {
        thisSub = nextSub;
        nextSub = uSub_STAR_STAR.Get(j+1,0);
        if( thisSub != Real(0) && nextSub != Real(0) )
            LogicError
            ("Quasi-triangular assumption broken at j=",j,
             ": subdiagonals were ",thisSub," and ",nextSub);
    }
}

} // namespace schur
} // namespace elem

#include "./Schur/QR.hpp"
#include "./Schur/SDC.hpp"
#include "./Schur/InverseFreeSDC.hpp"

namespace elem {

template<typename F>
inline void
Schur( Matrix<F>& A, Matrix<Complex<BASE(F)>>& w )
{
    DEBUG_ONLY(CallStackEntry cse("Schur"))
    schur::QR( A, w );
}

template<typename F>
inline void
Schur( Matrix<F>& A, Matrix<Complex<BASE(F)>>& w, Matrix<F>& Q )
{
    DEBUG_ONLY(CallStackEntry cse("Schur"))
    schur::QR( A, w, Q );
}

template<typename F>
inline void
Schur( DistMatrix<F>& A, DistMatrix<Complex<BASE(F)>,VR,STAR>& w )
{
    DEBUG_ONLY(CallStackEntry cse("Schur"))
#ifdef ELEM_HAVE_SCALAPACK
    schur::QR( A, w );
#else
    schur::SDC( A, w );
#endif
}

template<typename F>
inline void
Schur
( DistMatrix<F>& A, DistMatrix<Complex<BASE(F)>,VR,STAR>& w, DistMatrix<F>& Q )
{
    DEBUG_ONLY(CallStackEntry cse("Schur"))
#ifdef ELEM_HAVE_SCALAPACK
    schur::QR( A, w, Q );
#else
    schur::SDC( A, w, Q );
#endif
}

} // namespace elem

#endif // ifndef ELEM_SCHUR_HPP
