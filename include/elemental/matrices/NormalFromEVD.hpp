/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_NORMALFROMEVD_HPP
#define LAPACK_NORMALFROMEVD_HPP

#include "elemental/blas-like/level1/DiagonalScale.hpp"
#include "elemental/blas-like/level3/Gemm.hpp"
#include "elemental/matrices/Zeros.hpp"

namespace elem {

// A :=  Z Omega Z^T, where Omega is complex-valued and diagonal

template<typename R>
inline void
NormalFromEVD
(       Matrix<Complex<R> >& A,
  const Matrix<Complex<R> >& w,
  const Matrix<Complex<R> >& Z )
{
    DEBUG_ONLY(CallStackEntry cse("NormalFromEVD"))
    typedef Complex<R> C;

    Matrix<C> Y1, Z1Copy;

    const Int m = Z.Height();
    const Int n = Z.Width();
    Zeros( A, m, m );
    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);
        auto Z1 = LockedView( Z, 0, k, m,  nb );
        auto w1 = LockedView( w, k, 0, nb, 1  );

        Y1 = Z1Copy = Z1;
        DiagonalScale( RIGHT, ADJOINT, w1, Y1 );
        Gemm( NORMAL, NORMAL, C(1), Z1Copy, Y1, C(1), A );
    }
}

template<typename R>
inline Matrix<Complex<R> >
NormalFromEVD
( const Matrix<Complex<R> >& w,
  const Matrix<Complex<R> >& Z )
{
    Matrix<Complex<R> > A;
    NormalFromEVD( A, w, Z );
    return A;
}

template<typename R>
inline void
NormalFromEVD
(       DistMatrix<Complex<R> >& A,
  const DistMatrix<Complex<R>,VR,STAR>& w,
  const DistMatrix<Complex<R> >& Z )
{
    DEBUG_ONLY(CallStackEntry cse("NormalFromEVD"))
    typedef Complex<R> C;
    const Grid& g = A.Grid();
    DistMatrix<C,MC,  STAR> Z1_MC_STAR(g);
    DistMatrix<C,VR,  STAR> Z1_VR_STAR(g);
    DistMatrix<C,STAR,MR  > Z1Adj_STAR_MR(g);
    DistMatrix<C,STAR,STAR> w1_STAR_STAR(g);

    const Int m = Z.Height();
    const Int n = Z.Width();
    Zeros( A, m, m );
    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);
        auto Z1 = LockedView( Z, 0, k, m,  nb );
        auto w1 = LockedView( w, k, 0, nb, 1  );

        Z1_MC_STAR.AlignWith( A );
        Z1_MC_STAR = Z1;
        Z1_VR_STAR.AlignWith( A );
        Z1_VR_STAR = Z1_MC_STAR;
        w1_STAR_STAR = w1;

        DiagonalScale( RIGHT, ADJOINT, w1_STAR_STAR, Z1_VR_STAR );

        Z1Adj_STAR_MR.AlignWith( A );
        Z1Adj_STAR_MR.AdjointFrom( Z1_VR_STAR );
        LocalGemm( NORMAL, NORMAL, C(1), Z1_MC_STAR, Z1Adj_STAR_MR, C(1), A );
    }
}

template<typename R>
inline DistMatrix<Complex<R> >
NormalFromEVD
( const DistMatrix<Complex<R>,VR,STAR>& w,
  const DistMatrix<Complex<R> >& Z )
{
    DistMatrix<Complex<R>> A( w.Grid() );
    NormalFromEVD( A, w, Z );
    return A;
}

} // namespace elem

#endif // ifndef LAPACK_NORMALFROMEVD_HPP
