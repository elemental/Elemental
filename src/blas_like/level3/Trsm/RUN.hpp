/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace El {
namespace trsm {

// Right Upper Normal (Non)Unit Trsm
//   X := X triu(U)^-1, and
//   X := X triuu(U)^-1
template<typename F>
inline void
RUN
( UnitOrNonUnit diag, 
  const AbstractDistMatrix<F>& UPre,
        AbstractDistMatrix<F>& XPre,
  bool checkIfSingular )
{
    DEBUG_ONLY(CSE cse("trsm::RUN"))
    const Int n = XPre.Width();
    const Int bsize = Blocksize();
    const Grid& g = UPre.Grid();

    DistMatrixReadProxy<F,F,MC,MR> UProx( UPre );
    DistMatrixReadWriteProxy<F,F,MC,MR> XProx( XPre );
    auto& U = UProx.GetLocked();
    auto& X = XProx.Get();

    DistMatrix<F,STAR,STAR> U11_STAR_STAR(g); 
    DistMatrix<F,STAR,MR  > U12_STAR_MR(g);
    DistMatrix<F,VC,  STAR> X1_VC_STAR(g);    
    DistMatrix<F,STAR,MC  > X1Trans_STAR_MC(g);

    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const Range<Int> ind1( k,    k+nb ),
                         ind2( k+nb, n    );

        auto U11 = U( ind1, ind1 );
        auto U12 = U( ind1, ind2 );

        auto X1 = X( ALL, ind1 );
        auto X2 = X( ALL, ind2 );

        U11_STAR_STAR = U11; 
        X1_VC_STAR.AlignWith( X2 );
        X1_VC_STAR = X1;

        LocalTrsm
        ( RIGHT, UPPER, NORMAL, diag, F(1), U11_STAR_STAR, X1_VC_STAR,
          checkIfSingular );

        X1Trans_STAR_MC.AlignWith( X2 );
        Transpose( X1_VC_STAR, X1Trans_STAR_MC );
        Transpose( X1Trans_STAR_MC, X1 );
        U12_STAR_MR.AlignWith( X2 );
        U12_STAR_MR = U12; 

        // X2[MC,MR] -= X1[MC,* ] U12[* ,MR]
        //            = X1^T[* ,MC] U12[* ,MR]
        LocalGemm
        ( TRANSPOSE, NORMAL, F(-1), X1Trans_STAR_MC, U12_STAR_MR, F(1), X2 );
    }
}

} // namespace trsm
} // namespace El
