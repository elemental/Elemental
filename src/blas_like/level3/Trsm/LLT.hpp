/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace El {
namespace trsm {

// Left Lower (Conjugate)Transpose (Non)Unit Trsm
//   X := tril(L)^-T,
//   X := tril(L)^-H,
//   X := trilu(L)^-T, or
//   X := trilu(L)^-H

// width(X) >> p
template<typename F>
inline void
LLTLarge
( Orientation orientation,
  UnitOrNonUnit diag,
  const AbstractDistMatrix<F>& LPre,
        AbstractDistMatrix<F>& XPre, 
  bool checkIfSingular )
{
    DEBUG_ONLY(
      CSE cse("trsm::LLTLarge");
      if( orientation == NORMAL )
          LogicError("Expected (Conjugate)Transpose option");
    )
    const Int m = XPre.Height();
    const Int bsize = Blocksize();
    const Grid& g = LPre.Grid();

    DistMatrixReadProxy<F,F,MC,MR> LProx( LPre );
    DistMatrixReadWriteProxy<F,F,MC,MR> XProx( XPre );
    auto& L = LProx.GetLocked();
    auto& X = XProx.Get();

    DistMatrix<F,STAR,MC  > L10_STAR_MC(g);
    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<F,STAR,MR  > X1_STAR_MR(g);
    DistMatrix<F,STAR,VR  > X1_STAR_VR(g);

    const Int kLast = LastOffset( m, bsize );
    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,m-k);

        const Range<Int> ind0( 0, k    ),
                         ind1( k, k+nb );

        auto L10 = L( ind1, ind0 );
        auto L11 = L( ind1, ind1 );

        auto X0 = X( ind0, ALL );
        auto X1 = X( ind1, ALL );

        L11_STAR_STAR = L11; // L11[* ,* ] <- L11[MC,MR]
        X1_STAR_VR    = X1;  // X1[* ,VR] <- X1[MC,MR]

        // X1[* ,VR] := L11^-[T/H][* ,* ] X1[* ,VR]
        LocalTrsm
        ( LEFT, LOWER, orientation, diag, F(1), L11_STAR_STAR, X1_STAR_VR,
          checkIfSingular );

        X1_STAR_MR.AlignWith( X0 );
        X1_STAR_MR  = X1_STAR_VR; // X1[* ,MR] <- X1[* ,VR]
        X1          = X1_STAR_MR; // X1[MC,MR] <- X1[* ,MR]
        L10_STAR_MC.AlignWith( X0 );
        L10_STAR_MC = L10;        // L10[* ,MC] <- L10[MC,MR]

        // X0[MC,MR] -= (L10[* ,MC])^(T/H) X1[* ,MR]
        //            = L10^[T/H][MC,* ] X1[* ,MR]
        LocalGemm
        ( orientation, NORMAL, F(-1), L10_STAR_MC, X1_STAR_MR, F(1), X0 );
    }
}

// width(X) ~= p
template<typename F>
inline void
LLTMedium
( Orientation orientation,
  UnitOrNonUnit diag,
  const AbstractDistMatrix<F>& LPre,
        AbstractDistMatrix<F>& XPre, 
  bool checkIfSingular )
{
    DEBUG_ONLY(
      CSE cse("trsm::LLTMedium");
      if( orientation == NORMAL )
          LogicError("Expected (Conjugate)Transpose option");
    )
    const Int m = XPre.Height();
    const Int bsize = Blocksize();
    const Grid& g = LPre.Grid();

    DistMatrixReadProxy<F,F,MC,MR> LProx( LPre );
    DistMatrixReadWriteProxy<F,F,MC,MR> XProx( XPre );
    auto& L = LProx.GetLocked();
    auto& X = XProx.Get();

    DistMatrix<F,STAR,MC  > L10_STAR_MC(g);
    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<F,MR,  STAR> X1Trans_MR_STAR(g);

    const Int kLast = LastOffset( m, bsize );
    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,m-k);

        const Range<Int> ind0( 0, k    ),
                         ind1( k, k+nb );

        auto L10 = L( ind1, ind0 );
        auto L11 = L( ind1, ind1 );

        auto X0 = X( ind0, ALL );
        auto X1 = X( ind1, ALL );

        L11_STAR_STAR = L11; // L11[* ,* ] <- L11[MC,MR]
        // X1[* ,MR] <- X1[MC,MR]
        X1Trans_MR_STAR.AlignWith( X0 );
        Transpose( X1, X1Trans_MR_STAR, (orientation==ADJOINT) );

        // X1[* ,MR] := L11^-[T/H][* ,* ] X1[* ,MR]
        // X1^[T/H][MR,* ] := X1^[T/H][MR,* ] L11^-1[* ,* ]
        LocalTrsm
        ( RIGHT, LOWER, NORMAL, diag, 
          F(1), L11_STAR_STAR, X1Trans_MR_STAR, checkIfSingular );

        Transpose( X1Trans_MR_STAR, X1, (orientation==ADJOINT) );
        L10_STAR_MC.AlignWith( X0 );
        L10_STAR_MC = L10; // L10[* ,MC] <- L10[MC,MR]

        // X0[MC,MR] -= (L10[* ,MC])^[T/H] X1[* ,MR]
        //            = L10^[T/H][MC,* ] X1[* ,MR]
        LocalGemm
        ( orientation, orientation, 
          F(-1), L10_STAR_MC, X1Trans_MR_STAR, F(1), X0 );
    }
}

// width(X) << p
template<typename F,Dist colDist>
inline void
LLTSmall
( Orientation orientation,
  UnitOrNonUnit diag,
  const DistMatrix<F,colDist,STAR>& L,
        DistMatrix<F,colDist,STAR>& X,
  bool checkIfSingular )
{
    DEBUG_ONLY(
      CSE cse("trsm::LLTSmall");
      AssertSameGrids( L, X );
      if( orientation == NORMAL )
          LogicError("Expected (Conjugate)Transpose option");
      if( L.Height() != L.Width() || L.Height() != X.Height() )
          LogicError
          ("Nonconformal: \n",DimsString(L,"L"),"\n",DimsString(X,"X"));
      if( L.ColAlign() != X.ColAlign() )
          LogicError("L and X must be aligned");
    )
    const Int m = X.Height();
    const Int bsize = Blocksize();
    const Grid& g = L.Grid();

    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g), Z1_STAR_STAR(g);

    const Int kLast = LastOffset( m, bsize );
    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,m-k);

        const Range<Int> ind1( k,    k+nb ),
                         ind2( k+nb, m    );

        auto L11 = L( ind1, ind1 );
        auto L21 = L( ind2, ind1 );

        auto X1 = X( ind1, ALL );
        auto X2 = X( ind2, ALL );

        // X1 -= L21' X2
        LocalGemm( orientation, NORMAL, F(-1), L21, X2, Z1_STAR_STAR );
        axpy::util::UpdateWithLocalData( F(1), X1, Z1_STAR_STAR );
        El::AllReduce( Z1_STAR_STAR, X1.DistComm() );

        // X1 := L11^-1 X1
        L11_STAR_STAR = L11;
        LocalTrsm
        ( LEFT, LOWER, orientation, diag, F(1), L11_STAR_STAR, Z1_STAR_STAR,
          checkIfSingular );
        X1 = Z1_STAR_STAR;
    }
}

template<typename F,Dist rowDist>
inline void
LLTSmall
( Orientation orientation,
  UnitOrNonUnit diag,
  const DistMatrix<F,STAR,rowDist>& L,
        DistMatrix<F,rowDist,STAR>& X,
  bool checkIfSingular )
{
    DEBUG_ONLY(
      CSE cse("trsm::LLTSmall");
      AssertSameGrids( L, X );
      if( orientation == NORMAL )
          LogicError("Expected (Conjugate)Transpose option");
      if( L.Height() != L.Width() || L.Height() != X.Height() )
          LogicError
          ("Nonconformal: \n",DimsString(L,"L"),"\n",DimsString(X,"X"));
      if( L.RowAlign() != X.ColAlign() )
          LogicError("L and X must be aligned");
    )
    const Int m = X.Height();
    const Int bsize = Blocksize();
    const Grid& g = L.Grid();

    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g), X1_STAR_STAR(g);

    const Int kLast = LastOffset( m, bsize );
    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,m-k);

        const Range<Int> ind0( 0, k    ),
                         ind1( k, k+nb );

        auto L10 = L( ind1, ind0 );
        auto L11 = L( ind1, ind1 );

        auto X0 = X( ind0, ALL );
        auto X1 = X( ind1, ALL );

        L11_STAR_STAR = L11; // L11[* ,* ] <- L11[* ,VR]
        X1_STAR_STAR = X1;   // X1[* ,* ] <- X1[VR,* ]

        // X1[* ,* ] := L11^-[T/H][* ,* ] X1[* ,* ]
        LocalTrsm
        ( LEFT, LOWER, orientation, diag,
          F(1), L11_STAR_STAR, X1_STAR_STAR, checkIfSingular );

        X1 = X1_STAR_STAR;

        // X0[VR,* ] -= L10[* ,VR]^(T/H) X1[* ,* ]
        LocalGemm( orientation, NORMAL, F(-1), L10, X1_STAR_STAR, F(1), X0 );
    }
}

} // namespace trsm
} // namespace El
