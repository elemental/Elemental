/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace El {
namespace mstrsm {

template<typename F>
inline void
LUT
( Orientation orientation, F alpha, 
  Matrix<F>& U, const Matrix<F>& shifts, Matrix<F>& X ) 
{
    DEBUG_ONLY(CallStackEntry cse("mstrsm::LUT"))
    Scale( alpha, X );

    const Int m = X.Height();
    const Int n = X.Width();
    const Int bsize = Blocksize();
    for( Int k=0; k<m; k+=bsize )
    {
        const Int nb = Min(bsize,m-k);

        auto U11 =       ViewRange( U, k, k,    k+nb, k+nb );
        auto U12 = LockedViewRange( U, k, k+nb, k+nb, m    );

        auto X1 = ViewRange( X, k,    0, k+nb, n );
        auto X2 = ViewRange( X, k+nb, 0, m,    n );

        LeftUnb( UPPER, orientation, F(1), U11, shifts, X1 );
        Gemm( orientation, NORMAL, F(-1), U12, X1, F(1), X2 );
    }
}

template<typename F>
inline void
LUT
( Orientation orientation, F alpha, 
  const DistMatrix<F>& U, const DistMatrix<F,VR,STAR>& shifts,
        DistMatrix<F>& X ) 
{
    DEBUG_ONLY(CallStackEntry cse("mstrsm::LUT"))
    Scale( alpha, X );

    const Grid& g = U.Grid();
    DistMatrix<F,STAR,STAR> U11_STAR_STAR(g);
    DistMatrix<F,STAR,MC  > U12_STAR_MC(g);
    DistMatrix<F,STAR,MR  > X1_STAR_MR(g);
    DistMatrix<F,STAR,VR  > X1_STAR_VR(g);

    const Int m = X.Height();
    const Int n = X.Width();
    const Int bsize = Blocksize();
    for( Int k=0; k<m; k+=bsize )
    {
        const Int nb = Min(bsize,m-k);

        auto U11 = LockedViewRange( U, k, k,    k+nb, k+nb );
        auto U12 = LockedViewRange( U, k, k+nb, k+nb, m    );

        auto X1 = ViewRange( X, k,    0, k+nb, n );
        auto X2 = ViewRange( X, k+nb, 0, m,    n );

        // X1[* ,VR] := U11^-'[*,*] X1[* ,VR]
        U11_STAR_STAR = U11; // U11[* ,* ] <- U11[MC,MR]
        X1_STAR_VR.AlignWith( shifts );
        X1_STAR_VR = X1;  // X1[* ,VR] <- X1[MC,MR]
        LUT
        ( orientation, F(1), 
          U11_STAR_STAR.Matrix(), shifts.LockedMatrix(), X1_STAR_VR.Matrix() );

        X1_STAR_MR.AlignWith( X2 );
        X1_STAR_MR  = X1_STAR_VR; // X1[* ,MR]  <- X1[* ,VR]
        X1          = X1_STAR_MR; // X1[MC,MR]  <- X1[* ,MR]

        // X2[MC,MR] -= (U12[* ,MC])' X1[* ,MR]
        //            = U12'[MC,*] X1[* ,MR]
        U12_STAR_MC.AlignWith( X2 );
        U12_STAR_MC = U12; // U12[* ,MC] <- U12[MC,MR]
        LocalGemm
        ( orientation, NORMAL, F(-1), U12_STAR_MC, X1_STAR_MR, F(1), X2 );
    }
}

} // namespace mstrsm
} // namespace El
