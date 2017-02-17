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
void LVar1( Matrix<F>& L, bool conjugate=false )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
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

        const Range<Int> ind0( 0, k    ),
                         ind1( k, k+nb );

        auto L00 = L( ind0, ind0 );
        auto L10 = L( ind1, ind0 );
        auto L11 = L( ind1, ind1 );
        auto d1 = GetDiagonal(L11);
       
        S10 = L10;
        DiagonalSolve( LEFT, NORMAL, d1, L10, true );
        Trrk( LOWER, orientation, NORMAL, F(1), S10, L10, F(1), L00 );
        Trmm( LEFT, LOWER, orientation, UNIT, F(1), L11, L10 );
        trdtrmm::LUnblocked( L11, conjugate );
    }
}

template<typename F>
void LVar1( Matrix<F>& L, const Matrix<F>& dSub, bool conjugate=false )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
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
        auto dSub1 = dSub( IR(k,k+nb-1), ALL );

        const Range<Int> ind0( 0, k    ),
                         ind1( k, k+nb );

        auto L00 = L( ind0, ind0 );
        auto L10 = L( ind1, ind0 );
        auto L11 = L( ind1, ind1 );
        auto d1 = GetDiagonal(L11);

        S10 = L10;
        QuasiDiagonalSolve( LEFT, LOWER, d1, dSub1, L10, conjugate );
        Trrk( LOWER, orientation, NORMAL, F(1), S10, L10, F(1), L00 );
        Trmm( LEFT, LOWER, orientation, UNIT, F(1), L11, L10 );
        trdtrmm::LUnblocked( L11, dSub1, conjugate );

        k += nb;
    }
}

template<typename F>
void LVar1( AbstractDistMatrix<F>& LPre, bool conjugate=false )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( LPre.Height() != LPre.Width() )
          LogicError("L must be square");
    )
    const Int n = LPre.Height();
    const Int bsize = Blocksize();
    const Grid& g = LPre.Grid();
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );

    DistMatrixReadWriteProxy<F,F,MC,MR> LProx( LPre );
    auto& L = LProx.Get();

    DistMatrix<F,STAR,VR  > L10_STAR_VR(g);
    DistMatrix<F,STAR,VC  > S10_STAR_VC(g);
    DistMatrix<F,STAR,MC  > S10_STAR_MC(g);
    DistMatrix<F,STAR,MR  > L10_STAR_MR(g);
    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g);

    L10_STAR_VR.AlignWith( L );
    S10_STAR_VC.AlignWith( L );
    S10_STAR_MC.AlignWith( L );
    L10_STAR_MR.AlignWith( L );

    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const Range<Int> ind0( 0, k    ),
                         ind1( k, k+nb );

        auto L00 = L( ind0, ind0 );
        auto L10 = L( ind1, ind0 );
        auto L11 = L( ind1, ind1 );
        auto d1 = GetDiagonal(L11);

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

        Trdtrmm( LOWER, L11_STAR_STAR, conjugate );
        L11 = L11_STAR_STAR;
    }
}

template<typename F>
void LVar1
(       AbstractDistMatrix<F>& LPre,
  const AbstractDistMatrix<F>& dSubPre, 
  bool conjugate=false )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( LPre.Height() != LPre.Width() )
          LogicError("L must be square");
    )
    const Int n = LPre.Height();
    const Int bsize = Blocksize();
    const Grid& g = LPre.Grid();
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );

    DistMatrixReadWriteProxy<F,F,MC,MR> LProx( LPre );
    DistMatrixReadProxy<F,F,MD,STAR> dSubProx( dSubPre );
    auto& L = LProx.Get();
    auto& dSub = dSubProx.GetLocked();

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

    Int k=0;
    while( k < n )
    {
        const Int nbProp = Min(bsize,n-k);
        const bool in2x2 = ( k+nbProp<n && dSub.Get(k+nbProp-1,0) != F(0) );
        const Int nb = ( in2x2 ? nbProp+1 : nbProp );

        const Range<Int> ind0( 0, k    ),
                         ind1( k, k+nb );

        auto L00 = L( ind0, ind0 );
        auto L10 = L( ind1, ind0 );
        auto L11 = L( ind1, ind1 );
        auto d1 = GetDiagonal(L11);

        auto dSub1 = dSub( IR(k,k+nb-1), ALL );

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

        Trdtrmm( LOWER, L11_STAR_STAR, dSub1_STAR_STAR, conjugate );
        L11 = L11_STAR_STAR;

        k += nb;
    }
}

} // namespace trdtrmm
} // namespace El
