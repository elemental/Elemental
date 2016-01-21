/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_TRTRMM_UVAR1_HPP
#define EL_TRTRMM_UVAR1_HPP

namespace El {
namespace trtrmm {

template<typename T>
inline void
UVar1( Matrix<T>& U, bool conjugate=false )
{
    DEBUG_ONLY(CSE cse("trtrmm::UVar1"))
    const Int n = U.Height();
    const Int bsize = Blocksize();
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE ); 

    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const Range<Int> ind0( 0, k    ),
                         ind1( k, k+nb );

        auto U00 = U( ind0, ind0 );
        auto U01 = U( ind0, ind1 );
        auto U11 = U( ind1, ind1 );

        Trrk( UPPER, NORMAL, orientation, T(1), U01, U01, T(1), U00 );
        Trmm( RIGHT, UPPER, orientation, NON_UNIT, T(1), U11, U01 );
        trtrmm::UUnblocked( U11, conjugate );
    }
}

template<typename T>
inline void
UVar1( ElementalMatrix<T>& UPre, bool conjugate=false )
{
    DEBUG_ONLY(
      CSE cse("trtrmm::UVar1");
      if( UPre.Height() != UPre.Width() )
          LogicError("U must be square");
    )

    DistMatrixReadWriteProxy<T,T,MC,MR> UProx( UPre );
    auto& U = UProx.Get();

    const Grid& g = UPre.Grid();
    DistMatrix<T,MC,  STAR> U01_MC_STAR(g);
    DistMatrix<T,VC,  STAR> U01_VC_STAR(g);
    DistMatrix<T,VR,  STAR> U01_VR_STAR(g);
    DistMatrix<T,STAR,MR  > U01Trans_STAR_MR(g);
    DistMatrix<T,STAR,STAR> U11_STAR_STAR(g);

    U01_MC_STAR.AlignWith( U );
    U01_VC_STAR.AlignWith( U );
    U01_VR_STAR.AlignWith( U );
    U01Trans_STAR_MR.AlignWith( U );

    const Int n = UPre.Height();
    const Int bsize = Blocksize();
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const Range<Int> ind0( 0, k    ),
                         ind1( k, k+nb );

        auto U00 = U( ind0, ind0 );
        auto U01 = U( ind0, ind1 );
        auto U11 = U( ind1, ind1 );

        U01_MC_STAR = U01;
        U01_VC_STAR = U01_MC_STAR;
        U01_VR_STAR = U01_VC_STAR;
        Transpose( U01_VR_STAR, U01Trans_STAR_MR, conjugate );
        LocalTrrk( UPPER, T(1), U01_MC_STAR, U01Trans_STAR_MR, T(1), U00 );

        U11_STAR_STAR = U11;
        LocalTrmm
        ( RIGHT, UPPER, orientation, NON_UNIT, 
          T(1), U11_STAR_STAR, U01_VC_STAR );
        U01 = U01_VC_STAR;

        Trtrmm( UPPER, U11_STAR_STAR, conjugate );
        U11 = U11_STAR_STAR;
    }
}

} // namespace trtrmm
} // namespace El

#endif // ifndef EL_TRTRMM_UVAR1_HPP
