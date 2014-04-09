/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_TRSM_LUT_HPP
#define ELEM_TRSM_LUT_HPP

#include ELEM_GEMM_INC

namespace elem {
namespace trsm {

// Left Upper (Conjugate)Transpose (Non)Unit Trsm
//   X := triu(U)^-T  X, 
//   X := triu(U)^-H  X,
//   X := triuu(U)^-T X, or
//   X := triuu(U)^-H X

// width(X) >> p
template<typename F>
inline void
LUTLarge
( Orientation orientation, UnitOrNonUnit diag,
  const DistMatrix<F>& U, DistMatrix<F>& X,
  bool checkIfSingular )
{
    DEBUG_ONLY(
        CallStackEntry cse("trsm::LUTLarge");
        if( orientation == NORMAL )
            LogicError("Expected (Conjugate)Transpose option");
    )
    const Int m = X.Height();
    const Int n = X.Width();
    const Int bsize = Blocksize();
    const Grid& g = U.Grid();

    DistMatrix<F,STAR,STAR> U11_STAR_STAR(g); 
    DistMatrix<F,STAR,MC  > U12_STAR_MC(g);
    DistMatrix<F,STAR,MR  > X1_STAR_MR(g);
    DistMatrix<F,STAR,VR  > X1_STAR_VR(g);

    for( Int k=0; k<m; k+=bsize )
    {
        const Int nb = Min(bsize,m-k);

        auto U11 = LockedViewRange( U, k, k,    k+nb, k+nb );
        auto U12 = LockedViewRange( U, k, k+nb, k+nb, m    );

        auto X1 = ViewRange( X, k,    0, k+nb, n );
        auto X2 = ViewRange( X, k+nb, 0, m,    n );

        U11_STAR_STAR = U11; // U11[* ,* ] <- U11[MC,MR]
        X1_STAR_VR    = X1;  // X1[* ,VR] <- X1[MC,MR]
        
        // X1[* ,VR] := U11^-[T/H][*,*] X1[* ,VR]
        LocalTrsm
        ( LEFT, UPPER, orientation, diag, F(1), U11_STAR_STAR, X1_STAR_VR,
          checkIfSingular );

        X1_STAR_MR.AlignWith( X2 );
        X1_STAR_MR  = X1_STAR_VR; // X1[* ,MR]  <- X1[* ,VR]
        X1          = X1_STAR_MR; // X1[MC,MR]  <- X1[* ,MR]
        U12_STAR_MC.AlignWith( X2 );
        U12_STAR_MC = U12;        // U12[* ,MC] <- U12[MC,MR]

        // X2[MC,MR] -= (U12[* ,MC])^(T/H) X1[* ,MR]
        //            = U12^(T/H)[MC,*] X1[* ,MR]
        LocalGemm
        ( orientation, NORMAL, F(-1), U12_STAR_MC, X1_STAR_MR, F(1), X2 );
    }
}

// width(X) ~= p
template<typename F>
inline void
LUTMedium
( Orientation orientation, UnitOrNonUnit diag,
  const DistMatrix<F>& U, DistMatrix<F>& X, bool checkIfSingular )
{
    DEBUG_ONLY(
        CallStackEntry cse("trsm::LUTMedium");
        if( orientation == NORMAL )
            LogicError("Expected (Conjugate)Transpose option");
    )
    const Int m = X.Height();
    const Int n = X.Width();
    const Int bsize = Blocksize();
    const Grid& g = U.Grid();

    DistMatrix<F,STAR,STAR> U11_STAR_STAR(g); 
    DistMatrix<F,STAR,MC  > U12_STAR_MC(g);
    DistMatrix<F,MR,  STAR> X1Trans_MR_STAR(g);

    for( Int k=0; k<m; k+=bsize )
    {
        const Int nb = Min(bsize,m-k);

        auto U11 = LockedViewRange( U, k, k,    k+nb, k+nb );
        auto U12 = LockedViewRange( U, k, k+nb, k+nb, m    );

        auto X1 = ViewRange( X, k,    0, k+nb, n );
        auto X2 = ViewRange( X, k+nb, 0, m,    n );

        U11_STAR_STAR = U11; // U11[* ,* ] <- U11[MC,MR]
        // X1[* ,VR] <- X1[MC,MR]
        X1Trans_MR_STAR.AlignWith( X2 );
        X1.TransposeColAllGather( X1Trans_MR_STAR, (orientation==ADJOINT) );
        
        // X1[* ,MR] := U11^-[T/H][*,*] X1[* ,MR]
        // X1^[T/H][MR,* ] := X1^[T/H][MR,* ] U11^-1[* ,* ]
        LocalTrsm
        ( RIGHT, UPPER, NORMAL, diag, 
          F(1), U11_STAR_STAR, X1Trans_MR_STAR, checkIfSingular );

        X1.TransposeColFilterFrom( X1Trans_MR_STAR, (orientation==ADJOINT) );
        U12_STAR_MC.AlignWith( X2 );
        U12_STAR_MC = U12; // U12[* ,MC] <- U12[MC,MR]

        // X2[MC,MR] -= (U12[* ,MC])^[T/H] X1[* ,MR]
        //            = U12^[T/H][MC,*] X1[* ,MR]
        LocalGemm
        ( orientation, orientation, 
          F(-1), U12_STAR_MC, X1Trans_MR_STAR, F(1), X2 );
    }
}

// width(X) << p
template<typename F,Dist rowDist>
inline void
LUTSmall
( Orientation orientation, UnitOrNonUnit diag,
  const DistMatrix<F,STAR,rowDist>& U, DistMatrix<F,rowDist,STAR>& X,
  bool checkIfSingular )
{
    DEBUG_ONLY(
        CallStackEntry cse("trsm::LUTSmall");
        if( U.Grid() != X.Grid() )
            LogicError("U and X must be distributed over the same grid");
        if( orientation == NORMAL )
            LogicError("Expected (Conjugate)Transpose option");
        if( U.Height() != U.Width() || U.Height() != X.Height() )
            LogicError
            ("Nonconformal: \n",
             "  U ~ ",U.Height()," x ",U.Width(),"\n",
             "  X ~ ",X.Height()," x ",X.Width(),"\n");
        if( U.RowAlign() != X.ColAlign() )
            LogicError("U and X are assumed to be aligned");
    )
    const Int m = X.Height();
    const Int n = X.Width();
    const Int bsize = Blocksize();
    const Grid& g = U.Grid();

    DistMatrix<F,STAR,STAR> U11_STAR_STAR(g), X1_STAR_STAR(g); 

    for( Int k=0; k<m; k+=bsize )
    {
        const Int nb = Min(bsize,m-k);

        auto U11 = LockedViewRange( U, k, k,    k+nb, k+nb );
        auto U12 = LockedViewRange( U, k, k+nb, k+nb, m    );

        auto X1 = ViewRange( X, k,    0, k+nb, n );
        auto X2 = ViewRange( X, k+nb, 0, m,    n );

        U11_STAR_STAR = U11; // U11[* ,* ] <- U11[* ,VR]
        X1_STAR_STAR = X1;   // X1[* ,* ] <- X1[VR,* ]
        
        // X1[* ,* ] := U11^-[T/H][* ,* ] X1[* ,* ]
        LocalTrsm
        ( LEFT, UPPER, orientation, diag,
          F(1), U11_STAR_STAR, X1_STAR_STAR, checkIfSingular );

        X1 = X1_STAR_STAR;

        // X2[VR,* ] -= U12[* ,VR]^[T/H] X1[* ,* ]
        LocalGemm( orientation, NORMAL, F(-1), U12, X1_STAR_STAR, F(1), X2 );
    }
}

} // namespace trsm
} // namespace elem

#endif // ifndef ELEM_TRSM_LUT_HPP
