/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_TRTRMM_LVAR1_HPP
#define EL_TRTRMM_LVAR1_HPP

namespace El {
namespace trtrmm {

template<typename T>
void LVar1( Matrix<T>& L, bool conjugate=false )
{
    DEBUG_CSE
    const Int n = L.Height();
    const Int bsize = Blocksize();
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );
    
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);  

        const Range<Int> ind0( 0, k    ),
                         ind1( k, k+nb );

        auto L00 = L( ind0, ind0 );
        auto L10 = L( ind1, ind0 );
        auto L11 = L( ind1, ind1 );

        Trrk( LOWER, orientation, NORMAL, T(1), L10, L10, T(1), L00 );
        Trmm( LEFT, LOWER, orientation, NON_UNIT, T(1), L11, L10 );
        trtrmm::LUnblocked( L11, conjugate );
    }
}

template<typename T>
void LVar1( AbstractDistMatrix<T>& LPre, bool conjugate=false )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( LPre.Height() != LPre.Width() )
          LogicError("L must be square");
    )
    const Int n = LPre.Height();
    const Int bsize = Blocksize();
    const Grid& g = LPre.Grid();
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );

    DistMatrixReadWriteProxy<T,T,MC,MR> LProx( LPre );
    auto& L = LProx.Get();

    // Temporary distributions
    DistMatrix<T,STAR,VR  > L10_STAR_VR(g);
    DistMatrix<T,STAR,VC  > L10_STAR_VC(g);
    DistMatrix<T,STAR,MC  > L10_STAR_MC(g);
    DistMatrix<T,STAR,MR  > L10_STAR_MR(g);
    DistMatrix<T,STAR,STAR> L11_STAR_STAR(g);

    L10_STAR_VR.AlignWith( L );
    L10_STAR_VC.AlignWith( L );
    L10_STAR_MC.AlignWith( L );
    L10_STAR_MR.AlignWith( L );

    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const Range<Int> ind0( 0, k    ),
                         ind1( k, k+nb );

        auto L00 = L( ind0, ind0 );
        auto L10 = L( ind1, ind0 );
        auto L11 = L( ind1, ind1 );

        L10_STAR_VR = L10;
        L10_STAR_VC = L10_STAR_VR;
        L10_STAR_MC = L10_STAR_VC;
        L10_STAR_MR = L10_STAR_VR;
        LocalTrrk
        ( LOWER, orientation, T(1), L10_STAR_MC, L10_STAR_MR, T(1), L00 );

        L11_STAR_STAR = L11;
        LocalTrmm
        ( LEFT, LOWER, orientation, NON_UNIT, 
          T(1), L11_STAR_STAR, L10_STAR_VR );
        L10 = L10_STAR_VR;

        Trtrmm( LOWER, L11_STAR_STAR, conjugate );
        L11 = L11_STAR_STAR;
    }
}

} // namespace trtrmm
} // namespace El

#endif // ifndef EL_TRTRMM_LVAR1_HPP
