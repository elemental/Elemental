/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace El {
namespace trdtrmm {

template<typename F>
void UVar1( Matrix<F>& U, bool conjugate=false )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( U.Height() != U.Width() )
          LogicError("U must be square");
    )
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );

    Matrix<F> S01;

    const Int n = U.Height();
    const Int bsize = Blocksize();
    const Int kLast = LastOffset( n, bsize );
    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const Range<Int> ind0( 0, k    ),
                         ind1( k, k+nb );

        auto U00 = U( ind0, ind0 );
        auto U01 = U( ind0, ind1 );
        auto U11 = U( ind1, ind1 );
        auto d1 = GetDiagonal(U11);

        S01 = U01;
        DiagonalSolve( LEFT, NORMAL, d1, U01, true );
        Trrk( UPPER, NORMAL, orientation, F(1), U01, S01, F(1), U00 );
        Trmm( RIGHT, UPPER, orientation, UNIT, F(1), U11, U01 );
        trdtrmm::UUnblocked( U11, conjugate );
    }
}

template<typename F>
void UVar1( AbstractDistMatrix<F>& UPre, bool conjugate=false )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( UPre.Height() != UPre.Width() )
          LogicError("U must be square");
    )
    const Int n = UPre.Height();
    const Int bsize = Blocksize();
    const Grid& g = UPre.Grid();
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );

    DistMatrixReadWriteProxy<F,F,MC,MR> UProx( UPre );
    auto& U = UProx.Get();

    DistMatrix<F,MC,  STAR> S01_MC_STAR(g);
    DistMatrix<F,VC,  STAR> S01_VC_STAR(g);
    DistMatrix<F,VR,  STAR> U01_VR_STAR(g);
    DistMatrix<F,STAR,MR  > U01Trans_STAR_MR(g);
    DistMatrix<F,STAR,STAR> U11_STAR_STAR(g);

    S01_MC_STAR.AlignWith( U );
    S01_VC_STAR.AlignWith( U );
    U01_VR_STAR.AlignWith( U );
    U01Trans_STAR_MR.AlignWith( U );

    const Int kLast = LastOffset( n, bsize );
    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const Range<Int> ind0( 0, k    ),
                         ind1( k, k+nb );

        auto U00 = U( ind0, ind0 );
        auto U01 = U( ind0, ind1 );
        auto U11 = U( ind1, ind1 );
        auto d1 = GetDiagonal(U11);

        S01_MC_STAR = U01;
        S01_VC_STAR = S01_MC_STAR;
        U01_VR_STAR = S01_VC_STAR;
        DiagonalSolve( RIGHT, NORMAL, d1, U01_VR_STAR );
        Transpose( U01_VR_STAR, U01Trans_STAR_MR, conjugate );
        LocalTrrk( UPPER, F(1), S01_MC_STAR, U01Trans_STAR_MR, F(1), U00 );

        U11_STAR_STAR = U11;
        LocalTrmm
        ( RIGHT, UPPER, orientation, UNIT, F(1), U11_STAR_STAR, U01_VR_STAR );
        U01 = U01_VR_STAR;

        Trdtrmm( UPPER, U11_STAR_STAR, conjugate );
        U11 = U11_STAR_STAR;
    }
}

} // namespace trdtrmm
} // namespace El
