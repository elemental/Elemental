/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SCHUR_REALTOCOMPLEX_HPP
#define EL_SCHUR_REALTOCOMPLEX_HPP

namespace El {
namespace schur {

// TODO: Optimize these routines?

template<typename Real>
void RealToComplex( const Matrix<Real>& UQuasi, Matrix<Complex<Real>>& U )
{
    DEBUG_ONLY(CSE cse("schur::RealToComplex"))
    DEBUG_ONLY(CheckRealSchur(UQuasi))
    typedef Complex<Real> C;

    Copy( UQuasi, U );
    const Int n = U.Height();

    Matrix<C> V11(2,2), w1(2,1), U01Copy, U12Copy;
    const bool fullTriangle=true, multiplyQ=false;

    for( Int j=0; j<n-1; ++j )
    {
        if( U.Get(j+1,j) != Real(0) )
        {
            // Compute the Schur decomposition of the 2x2 block.
            // TODO: Switch to an analytical formula which exploits the fact
            //       that the block is in standard real Schur form
            lapack::HessenbergSchur
            ( 2, U.Buffer(j,j), U.LDim(), 
              w1.Buffer(), V11.Buffer(), V11.LDim(), fullTriangle, multiplyQ );
            U.Set(j+1,j,0);

            // Apply V11 from the right to U01
            auto U01 = U( IR(0,j), IR(j,j+2) );
            U01Copy = U01;
            Gemm( NORMAL, NORMAL, C(1), U01Copy, V11, C(0), U01 );

            // Apply V11^H from the left
            auto U12 = U( IR(j,j+2), IR(j+2,n) );
            U12Copy = U12;
            Gemm( ADJOINT, NORMAL, C(1), V11, U12Copy, C(0), U12 );
        }
    }
}

template<typename Real>
void RealToComplex
( const Matrix<Real>& UQuasi,
  const Matrix<Real>& QQuasi,
        Matrix<Complex<Real>>& U,
        Matrix<Complex<Real>>& Q )
{
    DEBUG_ONLY(CSE cse("schur::RealToComplex"))
    typedef Complex<Real> C;

    DEBUG_ONLY(CheckRealSchur(UQuasi))
    Copy( UQuasi, U );
    Copy( QQuasi, Q );
    const Int n = U.Height();

    Matrix<C> V11(2,2), w1(2,1), Q1Copy, U01Copy, U12Copy;
    const bool fullTriangle=true, multiplyV=false;

    for( Int j=0; j<n-1; ++j )
    {
        if( U.Get(j+1,j) != Real(0) )
        {
            // Compute the Schur decomposition of the 2x2 block.
            // TODO: Switch to an analytical formula which exploits the fact
            //       that the block is in standard real Schur form
            lapack::HessenbergSchur
            ( 2, U.Buffer(j,j), U.LDim(), 
              w1.Buffer(), V11.Buffer(), V11.LDim(), fullTriangle, multiplyV );
            U.Set(j+1,j,0);

            // Apply V11 from the right to Q1
            auto Q1 = Q( ALL, IR(j,j+2) );
            Q1Copy = Q1;
            Gemm( NORMAL, NORMAL, C(1), Q1Copy, V11, C(0), Q1 );

            // Apply V11 from the right to U01
            auto U01 = U( IR(0,j), IR(j,j+2) );
            U01Copy = U01;
            Gemm( NORMAL, NORMAL, C(1), U01Copy, V11, C(0), U01 );

            // Apply V11^H from the left to U12
            auto U12 = U( IR(j,j+2), IR(j+2,n) );
            U12Copy = U12;
            Gemm( ADJOINT, NORMAL, C(1), V11, U12Copy, C(0), U12 );
        }
    }
}

template<typename Real>
void RealToComplex
( const ElementalMatrix<        Real >& UQuasi, 
        ElementalMatrix<Complex<Real>>& UPre )
{
    DEBUG_ONLY(CSE cse("schur::RealToComplex"))
    typedef Complex<Real> C;

    DistMatrixWriteProxy<C,C,MC,MR> UProx( UPre );
    auto& U = UProx.Get();

    DEBUG_ONLY(CheckRealSchur(UQuasi))
    Copy( UQuasi, U );

    const Int n = U.Height();
    const Grid& g = U.Grid();

    DistMatrix<C> U01(g), U12(g);
    DistMatrix<C,STAR,STAR> U11_STAR_STAR(g), V11(2,2,g), w1(2,1,g);
    DistMatrix<C,VC,STAR> U01_VC_STAR(g), U01Copy_VC_STAR(g); 
    DistMatrix<C,STAR,VR> U12_STAR_VR(g), U12Copy_STAR_VR(g);
    const bool fullTriangle=true, multiplyQ=false;

    for( Int j=0; j<n-1; ++j )
    {
        if( U.Get(j+1,j) != Real(0) )
        {
            // Compute the Schur decomposition of the 2x2 block.
            // TODO: Switch to an analytical formula. Note that ScaLAPACK
            //       typically does not return the 2x2 blocks in standard form
            auto U11 = U( IR(j,j+2), IR(j,j+2) );
            U11_STAR_STAR = U11;
            lapack::HessenbergSchur
            ( 2, U11_STAR_STAR.Buffer(), U11_STAR_STAR.LDim(), 
              w1.Buffer(), V11.Buffer(), V11.LDim(), fullTriangle, multiplyQ );
            U11 = U11_STAR_STAR;
            U.Set(j+1,j,0);

            // Apply V11 from the right to U01
            auto U01 = U( IR(0,j), IR(j,j+2) );
            U01_VC_STAR = U01;
            U01Copy_VC_STAR = U01_VC_STAR;

            LocalGemm
            ( NORMAL, NORMAL, C(1), U01Copy_VC_STAR, V11, C(0), U01_VC_STAR );
            U01 = U01_VC_STAR;

            // Apply V11^H from the left to U12
            auto U12 = U( IR(j,j+2), IR(j+2,n) );
            U12_STAR_VR = U12;
            U12Copy_STAR_VR = U12_STAR_VR;
            LocalGemm
            ( ADJOINT, NORMAL, C(1), V11, U12Copy_STAR_VR, C(0), U12_STAR_VR );
            U12 = U12_STAR_VR;
        }
    }
}

template<typename Real>
void RealToComplex
( const ElementalMatrix<        Real >& UQuasi, 
  const ElementalMatrix<        Real >& QQuasi,
        ElementalMatrix<Complex<Real>>& UPre,
        ElementalMatrix<Complex<Real>>& QPre )
{
    DEBUG_ONLY(CSE cse("schur::RealToComplex"))
    typedef Complex<Real> C;

    DistMatrixWriteProxy<C,C,MC,MR> UProx( UPre ), QProx( QPre );
    auto& U = UProx.Get();
    auto& Q = QProx.Get();

    DEBUG_ONLY(CheckRealSchur(UQuasi))
    Copy( UQuasi, U );
    Copy( QQuasi, Q );

    const Int n = U.Height();
    const Grid& g = U.Grid();

    DistMatrix<C> U01(g), U12(g);
    DistMatrix<C,STAR,STAR> U11_STAR_STAR(g), V11(2,2,g), w1(2,1,g);
    DistMatrix<C,VC,STAR> Q1_VC_STAR(g), Q1Copy_VC_STAR(g),
                          U01_VC_STAR(g), U01Copy_VC_STAR(g); 
    DistMatrix<C,STAR,VR> U12_STAR_VR(g), U12Copy_STAR_VR(g);
    const bool fullTriangle=true, multiplyV=false;

    for( Int j=0; j<n-1; ++j )
    {
        if( U.Get(j+1,j) != Real(0) )
        {
            // Compute the Schur decomposition of the 2x2 block.
            // TODO: Switch to an analytical formula. Note that ScaLAPACK
            //       typically does not return the 2x2 blocks in standard form
            auto U11 = U( IR(j,j+2), IR(j,j+2) );
            U11_STAR_STAR = U11;
            lapack::HessenbergSchur
            ( 2, U11_STAR_STAR.Buffer(), U11_STAR_STAR.LDim(), 
              w1.Buffer(), V11.Buffer(), V11.LDim(), fullTriangle, multiplyV );
            U11 = U11_STAR_STAR;
            U.Set(j+1,j,0);

            // Apply V11 from the right to Q1
            auto Q1 = Q( ALL, IR(j,j+2) );
            Q1_VC_STAR = Q1;
            Q1Copy_VC_STAR = Q1_VC_STAR;
            LocalGemm
            ( NORMAL, NORMAL, C(1), Q1Copy_VC_STAR, V11, C(0), Q1_VC_STAR );
            Q1 = Q1_VC_STAR;

            // Apply V11 from the right to U01
            auto U01 = U( IR(0,j), IR(j,j+2) );
            U01_VC_STAR = U01;
            U01Copy_VC_STAR = U01_VC_STAR;
            LocalGemm
            ( NORMAL, NORMAL, C(1), U01Copy_VC_STAR, V11, C(0), U01_VC_STAR );
            U01 = U01_VC_STAR;

            // Apply V11^H from the left to U12
            auto U12 = U( IR(j,j+2), IR(j+2,n) );
            U12_STAR_VR = U12;
            U12Copy_STAR_VR = U12_STAR_VR;
            LocalGemm
            ( ADJOINT, NORMAL, C(1), V11, U12Copy_STAR_VR, C(0), U12_STAR_VR );
            U12 = U12_STAR_VR;
        }
    }
}

} // namespace schur
} // namespace El

#endif // ifndef EL_SCHUR_REALTOCOMPLEX_HPP
