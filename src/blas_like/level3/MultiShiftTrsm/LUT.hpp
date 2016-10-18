/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace El {
namespace mstrsm {

template<typename F>
void LUT
( Orientation orientation,
        Matrix<F>& U,
  const Matrix<F>& shifts,
        Matrix<F>& X ) 
{
    DEBUG_CSE

    const Int m = X.Height();
    const Int bsize = Blocksize();

    for( Int k=0; k<m; k+=bsize )
    {
        const Int nb = Min(bsize,m-k);

        const Range<Int> ind1( k,    k+nb ),
                         ind2( k+nb, m    );

        auto U11 = U( ind1, ind1 );
        auto U12 = U( ind1, ind2 );

        auto X1 = X( ind1, ALL );
        auto X2 = X( ind2, ALL );

        LeftUnb( UPPER, orientation, U11, shifts, X1 );
        Gemm( orientation, NORMAL, F(-1), U12, X1, F(1), X2 );
    }
}

template<typename F>
void LUT
( Orientation orientation,
  const AbstractDistMatrix<F>& UPre,
  const AbstractDistMatrix<F>& shiftsPre,
        AbstractDistMatrix<F>& XPre ) 
{
    DEBUG_CSE

    DistMatrixReadProxy<F,F,MC,MR> UProx( UPre );
    DistMatrixReadProxy<F,F,VR,STAR> shiftsProx( shiftsPre );
    DistMatrixReadWriteProxy<F,F,MC,MR> XProx( XPre );
    auto& U = UProx.GetLocked();
    auto& shifts = shiftsProx.GetLocked();
    auto& X = XProx.Get();

    const Grid& g = U.Grid();
    DistMatrix<F,STAR,STAR> U11_STAR_STAR(g);
    DistMatrix<F,STAR,MC  > U12_STAR_MC(g);
    DistMatrix<F,STAR,MR  > X1_STAR_MR(g);
    DistMatrix<F,STAR,VR  > X1_STAR_VR(g);

    const Int m = X.Height();
    const Int bsize = Blocksize();

    for( Int k=0; k<m; k+=bsize )
    {
        const Int nb = Min(bsize,m-k);

        const Range<Int> ind1( k,    k+nb ),
                         ind2( k+nb, m    );

        auto U11 = U( ind1, ind1 );
        auto U12 = U( ind1, ind2 );

        auto X1 = X( ind1, ALL );
        auto X2 = X( ind2, ALL );

        // X1[* ,VR] := U11^-'[*,*] X1[* ,VR]
        U11_STAR_STAR = U11; // U11[* ,* ] <- U11[MC,MR]
        X1_STAR_VR.AlignWith( shifts );
        X1_STAR_VR = X1;  // X1[* ,VR] <- X1[MC,MR]
        LUT
        ( orientation, 
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
