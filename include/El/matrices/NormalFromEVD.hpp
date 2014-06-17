/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_NORMALFROMEVD_HPP
#define EL_NORMALFROMEVD_HPP



namespace El {

// A :=  Z Omega Z^T, where Omega is complex-valued and diagonal

template<typename Real>
inline void
NormalFromEVD
(       Matrix<Complex<Real>>& A,
  const Matrix<Complex<Real>>& w,
  const Matrix<Complex<Real>>& Z )
{
    DEBUG_ONLY(CallStackEntry cse("NormalFromEVD"))
    typedef Complex<Real> C;

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

template<typename Real>
inline void
NormalFromEVD
(       DistMatrix<Complex<Real>>& A,
  const DistMatrix<Complex<Real>,VR,STAR>& w,
  const DistMatrix<Complex<Real>>& Z )
{
    DEBUG_ONLY(CallStackEntry cse("NormalFromEVD"))
    typedef Complex<Real> C;
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
        Z1_VR_STAR.AdjointPartialColAllGather( Z1Adj_STAR_MR );
        LocalGemm( NORMAL, NORMAL, C(1), Z1_MC_STAR, Z1Adj_STAR_MR, C(1), A );
    }
}

} // namespace El

#endif // ifndef EL_NORMALFROMEVD_HPP
