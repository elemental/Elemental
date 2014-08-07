/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace El {
namespace trsm {

// Right Upper (Conjugate)Transpose (Non)Unit Trsm
//   X := X triu(U)^-T, 
//   X := X triu(U)^-H,
//   X := X triuu(U)^-T, or
//   X := X triuu(U)^-H
template<typename F>
inline void
RUT
( Orientation orientation, UnitOrNonUnit diag,
  const AbstractDistMatrix<F>& UPre, AbstractDistMatrix<F>& XPre, 
  bool checkIfSingular )
{
    DEBUG_ONLY(
        CallStackEntry cse("trsm::RUT");
        if( orientation == NORMAL )
            LogicError("Expected (Conjugate)Transpose option");
    )
    const Int m = XPre.Height();
    const Int n = XPre.Width();
    const Int bsize = Blocksize();
    const Grid& g = UPre.Grid();

    DistMatrix<F> U(g), X(g);
    Copy( UPre, U, READ_PROXY );
    Copy( XPre, X, READ_WRITE_PROXY );

    DistMatrix<F,VR,  STAR> U01_VR_STAR(g);
    DistMatrix<F,STAR,MR  > U01Trans_STAR_MR(g);
    DistMatrix<F,STAR,STAR> U11_STAR_STAR(g);
    DistMatrix<F,VC,  STAR> X1_VC_STAR(g);
    DistMatrix<F,STAR,MC  > X1Trans_STAR_MC(g);
    
    const Int kLast = LastOffset( n, bsize );
    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,n-k);

        auto U01 = LockedViewRange( U, 0, k, k,    k+nb );
        auto U11 = LockedViewRange( U, k, k, k+nb, k+nb );

        auto X0 = ViewRange( X, 0, 0, m, k    );
        auto X1 = ViewRange( X, 0, k, m, k+nb );

        U11_STAR_STAR = U11;
        X1_VC_STAR.AlignWith( X0 );
        X1_VC_STAR = X1; 

        LocalTrsm
        ( RIGHT, UPPER, orientation, diag, 
          F(1), U11_STAR_STAR, X1_VC_STAR, checkIfSingular );

        X1Trans_STAR_MC.AlignWith( X0 );
        X1_VC_STAR.TransposePartialColAllGather( X1Trans_STAR_MC );
        X1.TransposeRowFilterFrom( X1Trans_STAR_MC );
        U01_VR_STAR.AlignWith( X0 );
        U01_VR_STAR = U01;
        U01Trans_STAR_MR.AlignWith( X0 );
        U01_VR_STAR.TransposePartialColAllGather
        ( U01Trans_STAR_MR, (orientation==ADJOINT) );

        // X0[MC,MR] -= X1[MC,* ] (U01[MR,* ])^(T/H)
        //            = X1^T[* ,MC] (U01^(T/H))[* ,MR]
        LocalGemm
        ( TRANSPOSE, NORMAL, 
          F(-1), X1Trans_STAR_MC, U01Trans_STAR_MR, F(1), X0 );
    }

    Copy( X, XPre, RESTORE_READ_WRITE_PROXY );
}

} // namespace trsm
} // namespace El
