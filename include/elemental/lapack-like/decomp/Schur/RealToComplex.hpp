/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_SCHUR_REALTOCOMPLEX_HPP
#define ELEM_SCHUR_REALTOCOMPLEX_HPP

#include ELEM_GEMM_INC

namespace elem {
namespace schur {

template<typename Real>
inline void
RealToComplex( const Matrix<Real>& UQuasi, Matrix<Complex<Real>>& U )
{
    DEBUG_ONLY(CallStackEntry cse("schur::RealToComplex"))
    DEBUG_ONLY(CheckRealSchur(UQuasi))
    typedef Complex<Real> C;

    Copy( UQuasi, U );
    const Int n = U.Height();

    Matrix<C> Q11(2,2), w1(2,1), U01Copy, U12Copy;
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
              w1.Buffer(), Q11.Buffer(), Q11.LDim(), fullTriangle, multiplyQ );
            U.Set(j+1,j,0);

            // Apply Q11 from the right
            auto U01 = ViewRange( U, 0, j, j, j+2 );
            U01Copy = U01;
            Gemm( NORMAL, NORMAL, C(1), U01Copy, Q11, C(0), U01 );

            // Apply Q11^H from the left
            auto U12 = ViewRange( U, j, j+2, j+2, n );
            U12Copy = U12;
            Gemm( ADJOINT, NORMAL, C(1), Q11, U12Copy, C(0), U12 );
        }
    }
}

template<typename Real>
inline void
RealToComplex( const DistMatrix<Real>& UQuasi, DistMatrix<Complex<Real>>& U )
{
    DEBUG_ONLY(CallStackEntry cse("schur::RealToComplex"))
    DEBUG_ONLY(CheckRealSchur(UQuasi))
    typedef Complex<Real> C;

    Copy( UQuasi, U );
    const Int n = U.Height();
    const Grid& g = U.Grid();

    DistMatrix<C> U01(g), U12(g);
    DistMatrix<C,STAR,STAR> U11_STAR_STAR(g), Q11(2,2,g), w1(2,1,g);
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
            auto U11 = ViewRange( U, j, j, j+2, j+2 );
            U11_STAR_STAR = U11;
            lapack::HessenbergSchur
            ( 2, U11_STAR_STAR.Buffer(), U11_STAR_STAR.LDim(), 
              w1.Buffer(), Q11.Buffer(), Q11.LDim(), fullTriangle, multiplyQ );
            U11 = U11_STAR_STAR;
            U.Set(j+1,j,0);

            // Apply Q11 from the right
            auto U01 = ViewRange( U, 0, j, j, j+2 );
            U01_VC_STAR = U01;
            U01Copy_VC_STAR = U01_VC_STAR;

            LocalGemm
            ( NORMAL, NORMAL, C(1), U01Copy_VC_STAR, Q11, C(0), U01_VC_STAR );
            U01 = U01_VC_STAR;

            // Apply Q11^H from the left
            auto U12 = ViewRange( U, j, j+2, j+2, n );
            U12_STAR_VR = U12;
            U12Copy_STAR_VR = U12_STAR_VR;
            LocalGemm
            ( ADJOINT, NORMAL, C(1), Q11, U12Copy_STAR_VR, C(0), U12_STAR_VR );
            U12 = U12_STAR_VR;
        }
    }
}

} // namespace schur
} // namespace elem

#endif // ifndef ELEM_SCHUR_REALTOCOMPLEX_HPP
