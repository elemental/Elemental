/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_INVERSE_TRIANGULAR_UVAR3_HPP
#define EL_INVERSE_TRIANGULAR_UVAR3_HPP

namespace El {
namespace triang_inv {

template<typename F>
void
UVar3Unb( UnitOrNonUnit diag, Matrix<F>& U )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( U.Height() != U.Width() )
          LogicError("Nonsquare matrices cannot be triangular");
    )
    const Int n = U.Height();
    const Int ldu = U.LDim();
    F* UBuffer = U.Buffer();
    for( Int j=n-1; j>=0; --j )
    {
        const F upsilon = ( diag==NON_UNIT ? UBuffer[j+j*ldu] : F(1) );
        for( Int k=0; k<j; ++k )
            UBuffer[k+j*ldu] /= -upsilon;
        blas::Geru
        ( j, n-(j+1), F(1),
          &UBuffer[j*ldu], 1, &UBuffer[j+(j+1)*ldu], ldu, 
          &UBuffer[(j+1)*ldu], ldu );
        if( diag == NON_UNIT )
        {
            for( Int k=j+1; k<n; ++k )
                UBuffer[j+k*ldu] /= upsilon;
            UBuffer[j+j*ldu] = F(1) / UBuffer[j+j*ldu];
        }
    }
}

template<typename F>
void
UVar3( UnitOrNonUnit diag, Matrix<F>& U )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( U.Height() != U.Width() )
          LogicError("Nonsquare matrices cannot be triangular");
    )
    const Int n = U.Height();
    const Int bsize = Blocksize();
    const Int kLast = LastOffset( n, bsize );
    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const Range<Int> ind0( 0,    k    ),
                         ind1( k,    k+nb ),
                         ind2( k+nb, n    );

        auto U01 = U( ind0, ind1 );
        auto U02 = U( ind0, ind2 );
        auto U11 = U( ind1, ind1 );
        auto U12 = U( ind1, ind2 );
        auto U22 = U( ind2, ind2 );

        Trsm( RIGHT, UPPER, NORMAL, diag, F(-1), U11, U01 );
        Gemm( NORMAL, NORMAL, F(1), U01, U12, F(1), U02 );
        Trsm( LEFT, UPPER, NORMAL, diag, F(1), U11, U12 );
        UVar3Unb( diag, U11 );
    }
}

template<typename F>
void
UVar3( UnitOrNonUnit diag, ElementalMatrix<F>& UPre )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( UPre.Height() != UPre.Width() )
          LogicError("Nonsquare matrices cannot be triangular");
    )

    DistMatrixReadWriteProxy<F,F,MC,MR> UProx( UPre );
    auto& U = UProx.Get();

    const Grid& g = U.Grid();
    DistMatrix<F,VC,  STAR> U01_VC_STAR(g);
    DistMatrix<F,STAR,STAR> U11_STAR_STAR(g);
    DistMatrix<F,STAR,VR  > U12_STAR_VR(g);
    DistMatrix<F,STAR,MC  > U01Trans_STAR_MC(g);
    DistMatrix<F,MR,  STAR> U12Trans_MR_STAR(g);

    const Int n = U.Height();
    const Int bsize = Blocksize();
    const Int kLast = LastOffset( n, bsize );
    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const Range<Int> ind0( 0,    k    ),
                         ind1( k,    k+nb ),
                         ind2( k+nb, n    );

        auto U01 = U( ind0, ind1 );
        auto U02 = U( ind0, ind2 );
        auto U11 = U( ind1, ind1 );
        auto U12 = U( ind1, ind2 );
        auto U22 = U( ind2, ind2 );

        U01_VC_STAR = U01;
        U11_STAR_STAR = U11;
        LocalTrsm
        ( RIGHT, UPPER, NORMAL, diag, F(-1), U11_STAR_STAR, U01_VC_STAR );

        // We transpose before the communication to avoid cache-thrashing
        // in the unpacking stage.
        U12Trans_MR_STAR.AlignWith( U02 );
        U01Trans_STAR_MC.AlignWith( U02 );
        Transpose( U12, U12Trans_MR_STAR );
        Transpose( U01_VC_STAR, U01Trans_STAR_MC );

        LocalGemm
        ( TRANSPOSE, TRANSPOSE, 
          F(1), U01Trans_STAR_MC, U12Trans_MR_STAR, F(1), U02 );
        Transpose( U01Trans_STAR_MC, U01 );

        Transpose( U12Trans_MR_STAR, U12_STAR_VR );
        LocalTrsm
        ( LEFT, UPPER, NORMAL, diag, F(1), U11_STAR_STAR, U12_STAR_VR );
        LocalTriangularInverse( UPPER, diag, U11_STAR_STAR );
        U11 = U11_STAR_STAR;
        U12 = U12_STAR_VR;
    }
}

} // namespace triang_inv
} // namespace El

#endif // ifndef EL_INVERSE_TRIANGULAR_UVAR3_HPP
