/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace El {
namespace trdtrmm {

template<typename F>
inline void
LVar1( Matrix<F>& L, bool conjugate=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("trdtrmm::LVar1");
        if( L.Height() != L.Width() )
            LogicError("L must be square");
    )
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );
    Matrix<F> S10;

    const Int n = L.Height();
    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        auto L00 = ViewRange( L, 0, 0, k,    k    );
        auto L10 = ViewRange( L, k, 0, k+nb, k    );
        auto L11 = ViewRange( L, k, k, k+nb, k+nb );
        auto d1 = L11.GetDiagonal();
       
        S10 = L10;
        DiagonalSolve( LEFT, NORMAL, d1, L10, true );
        Trrk( LOWER, orientation, NORMAL, F(1), S10, L10, F(1), L00 );
        Trmm( LEFT, LOWER, orientation, UNIT, F(1), L11, L10 );
        trdtrmm::LUnblocked( L11, conjugate );
    }
}

template<typename F>
inline void
LVar1( Matrix<F>& L, const Matrix<F>& dSub, bool conjugate=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("trdtrmm::LVar1");
        if( L.Height() != L.Width() )
            LogicError("L must be square");
    )
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );
    Matrix<F> S10;

    const Int n = L.Height();
    const Int bsize = Blocksize();
    Int k=0;
    while( k < n )
    {
        const Int nbProp = Min(bsize,n-k);
        const bool in2x2 = ( k+nbProp<n && dSub.Get(k+nbProp-1,0) != F(0) );
        const Int nb = ( in2x2 ? nbProp+1 : nbProp );
        auto dSub1 = LockedViewRange( dSub, k, 0, k+nb-1, 1 );

        auto L00 = ViewRange( L, 0, 0, k,    k    );
        auto L10 = ViewRange( L, k, 0, k+nb, k    );
        auto L11 = ViewRange( L, k, k, k+nb, k+nb );
        auto d1 = L11.GetDiagonal();

        S10 = L10;
        QuasiDiagonalSolve( LEFT, LOWER, d1, dSub1, L10, conjugate );
        Trrk( LOWER, orientation, NORMAL, F(1), S10, L10, F(1), L00 );
        Trmm( LEFT, LOWER, orientation, UNIT, F(1), L11, L10 );
        trdtrmm::LUnblocked( L11, dSub1, conjugate );

        k += nb;
    }
}

template<typename F>
inline void
LVar1( DistMatrix<F>& L, bool conjugate=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("trdtrmm::LVar1");
        if( L.Height() != L.Width() )
            LogicError("L must be square");
    )
    const Grid& g = L.Grid();
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );

    DistMatrix<F,STAR,VR  > L10_STAR_VR(g);
    DistMatrix<F,STAR,VC  > S10_STAR_VC(g);
    DistMatrix<F,STAR,MC  > S10_STAR_MC(g);
    DistMatrix<F,STAR,MR  > L10_STAR_MR(g);
    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g);

    L10_STAR_VR.AlignWith( L );
    S10_STAR_VC.AlignWith( L );
    S10_STAR_MC.AlignWith( L );
    L10_STAR_MR.AlignWith( L );

    const Int n = L.Height();
    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        auto L00 = ViewRange( L, 0, 0, k,    k    );
        auto L10 = ViewRange( L, k, 0, k+nb, k    );
        auto L11 = ViewRange( L, k, k, k+nb, k+nb );
        auto d1 = L11.GetDiagonal();

        L10_STAR_VR = L10;
        S10_STAR_VC = L10_STAR_VR;
        S10_STAR_MC = S10_STAR_VC;
        DiagonalSolve( LEFT, NORMAL, d1, L10_STAR_VR, true );
        L10_STAR_MR = L10_STAR_VR;
        LocalTrrk
        ( LOWER, orientation, F(1), S10_STAR_MC, L10_STAR_MR, F(1), L00 );

        L11_STAR_STAR = L11;
        LocalTrmm
        ( LEFT, LOWER, orientation, UNIT, F(1), L11_STAR_STAR, L10_STAR_VR );
        L10 = L10_STAR_VR;

        LocalTrdtrmm( LOWER, L11_STAR_STAR, conjugate );
        L11 = L11_STAR_STAR;
    }
}

template<typename F>
inline void
LVar1
( DistMatrix<F>& L, const DistMatrix<F,MD,STAR>& dSub, bool conjugate=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("trdtrmm::LVar1");
        if( L.Height() != L.Width() )
            LogicError("L must be square");
    )
    const Grid& g = L.Grid();
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );

    DistMatrix<F,STAR,VR  > L10_STAR_VR(g);
    DistMatrix<F,STAR,VC  > S10_STAR_VC(g);
    DistMatrix<F,STAR,MC  > S10_STAR_MC(g);
    DistMatrix<F,STAR,MR  > L10_STAR_MR(g);
    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g), 
                            d1_STAR_STAR(g), dSub1_STAR_STAR(g);

    L10_STAR_VR.AlignWith( L );
    S10_STAR_VC.AlignWith( L );
    S10_STAR_MC.AlignWith( L );
    L10_STAR_MR.AlignWith( L );

    const Int n = L.Height();
    const Int bsize = Blocksize();
    Int k=0;
    while( k < n )
    {
        const Int nbProp = Min(bsize,n-k);
        const bool in2x2 = ( k+nbProp<n && dSub.Get(k+nbProp-1,0) != F(0) );
        const Int nb = ( in2x2 ? nbProp+1 : nbProp );
        auto dSub1 = LockedViewRange( dSub, k, 0, k+nb-1, 1 );

        auto L00 = ViewRange( L, 0, 0, k,    k    );
        auto L10 = ViewRange( L, k, 0, k+nb, k    );
        auto L11 = ViewRange( L, k, k, k+nb, k+nb );
        auto d1 = L11.GetDiagonal();

        L10_STAR_VR = L10;
        S10_STAR_VC = L10_STAR_VR;
        S10_STAR_MC = S10_STAR_VC;
        d1_STAR_STAR = d1;
        dSub1_STAR_STAR = dSub1;
        // TODO: LocalQuasiDiagonalSolve?
        QuasiDiagonalSolve
        ( LEFT, LOWER,
          d1_STAR_STAR.LockedMatrix(), dSub1_STAR_STAR.LockedMatrix(), 
          L10_STAR_VR.Matrix(), conjugate );
        L10_STAR_MR = L10_STAR_VR;
        LocalTrrk
        ( LOWER, orientation, F(1), S10_STAR_MC, L10_STAR_MR, F(1), L00 );

        L11_STAR_STAR = L11;
        LocalTrmm
        ( LEFT, LOWER, orientation, UNIT, F(1), L11_STAR_STAR, L10_STAR_VR );
        L10 = L10_STAR_VR;

        LocalTrdtrmm( LOWER, L11_STAR_STAR, dSub1_STAR_STAR, conjugate );
        L11 = L11_STAR_STAR;

        k += nb;
    }
}

} // namespace trdtrmm
} // namespace El
