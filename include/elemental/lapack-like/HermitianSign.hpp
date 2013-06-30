/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_HERMITIANSIGN_HPP
#define LAPACK_HERMITIANSIGN_HPP

#include "elemental/lapack-like/HermitianFunction.hpp"

namespace elem {

// The Hermitian sign decomposition is equivalent to the Hermitian polar
// decomposition... A = (U sgn(Lambda) U') (U sgn(Lambda)Lambda U')
//                    = (U sgn(Lambda) U') (U |Lambda| U')

// Even though sgn(lambda) isn't well-defined when lambda=0, we will extend it
// from the right so that the sign decomposition of a singular Hermitian matrix
// is a polar decomposition (which always exists).

template<typename F>
inline void
HermitianSign( UpperOrLower uplo, Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("HermitianSign");
#endif
    typedef BASE(F) R;

    // Get the EVD of A
    Matrix<R> w;
    Matrix<F> Z;
    HermitianEig( uplo, A, w, Z );

    const int n = A.Height();
    for( int i=0; i<n; ++i )
    {
        const R omega = w.Get(i,0);
        if( omega >= 0 )
            w.Set(i,0,R(1));
        else 
            w.Set(i,0,R(-1));
    }

    // Reform the Hermitian matrix with the modified eigenvalues
    HermitianFromEVD( uplo, A, w, Z );
}

template<typename F>
inline void
HermitianSign( UpperOrLower uplo, Matrix<F>& A, Matrix<F>& N )
{
#ifndef RELEASE
    CallStackEntry entry("HermitianSign");
#endif
    typedef BASE(F) R;

    // Get the EVD of A
    Matrix<R> w;
    Matrix<F> Z;
    HermitianEig( uplo, A, w, Z );

    const int n = A.Height();
    Matrix<R> wSgn( n, 1 ), wAbs( n, 1 );
    for( int i=0; i<n; ++i )
    {
        const R omega = w.Get(i,0);
        if( omega >= 0 )
        {
            wSgn.Set(i,0,R(1));
            wAbs.Set(i,0,omega);
        }
        else
        {
            wSgn.Set(i,0,R(-1));
            wAbs.Set(i,0,-omega);
        }
    }

    // Form the Hermitian matrices with modified eigenvalues
    HermitianFromEVD( uplo, A, wSgn, Z );
    HermitianFromEVD( uplo, N, wAbs, Z );
}

template<typename F>
inline void
HermitianSign( UpperOrLower uplo, DistMatrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("HermitianSign");
#endif
    EnsurePMRRR();
    typedef BASE(F) R;

    // Get the EVD of A
    const Grid& g = A.Grid();
    DistMatrix<R,VR,STAR> w(g);
    DistMatrix<F> Z(g);
    HermitianEig( uplo, A, w, Z );

    const int n = A.Height();
    const int numLocalEigs = w.LocalHeight();
    for( int iLoc=0; iLoc<numLocalEigs; ++iLoc )
    {
        const R omega = w.GetLocal(iLoc,0);
        if( omega >= 0 )
            w.SetLocal(iLoc,0,R(1));
        else
            w.SetLocal(iLoc,0,R(-1));
    }

    // Reform the Hermitian matrix with the modified eigenvalues
    HermitianFromEVD( uplo, A, w, Z );
}

template<typename F>
inline void
HermitianSign( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<F>& N )
{
#ifndef RELEASE
    CallStackEntry entry("HermitianSign");
#endif
    EnsurePMRRR();
    typedef BASE(F) R;

    // Get the EVD of A
    const Grid& g = A.Grid();
    DistMatrix<R,VR,STAR> w(g);
    DistMatrix<F> Z(g);
    HermitianEig( uplo, A, w, Z );

    const int n = A.Height();
    const int numLocalEigs = w.LocalHeight();
    DistMatrix<R,VR,STAR> wSgn(g), wAbs(g);
    wSgn.AlignWith( w );
    wAbs.AlignWith( w );
    wSgn.ResizeTo( n, 1 );
    wAbs.ResizeTo( n, 1 );
    for( int iLoc=0; iLoc<numLocalEigs; ++iLoc )
    {
        const R omega = w.GetLocal(iLoc,0);
        if( omega >= 0 )
        {
            wSgn.SetLocal(iLoc,0,R(1));
            wAbs.SetLocal(iLoc,0,omega);
        }
        else
        {
            wSgn.SetLocal(iLoc,0,R(-1));
            wAbs.SetLocal(iLoc,0,-omega);
        }
    }

    // Form the Hermitian matrix with the modified eigenvalues
    HermitianFromEVD( uplo, A, wSgn, Z );
    HermitianFromEVD( uplo, N, wAbs, Z );
}

} // namespace elem

#endif // ifndef LAPACK_HERMITIANSIGN_HPP
