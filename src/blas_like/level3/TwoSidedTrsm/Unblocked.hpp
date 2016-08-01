/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_TWOSIDEDTRSM_UNBLOCKED_HPP
#define EL_TWOSIDEDTRSM_UNBLOCKED_HPP

namespace El {
namespace twotrsm {

template<typename F>
void LUnb( UnitOrNonUnit diag, Matrix<F>& A, const Matrix<F>& L )
{
    DEBUG_CSE
    // Use the Variant 4 algorithm
    const Int n = A.Height();
    const Int lda = A.LDim();
    const Int ldl = L.LDim();
    F* ABuffer = A.Buffer();
    const F* LBuffer = L.LockedBuffer();
    for( Int j=0; j<n; ++j )
    {
        const Int a21Height = n - (j+1);

        // Extract and store the diagonal value of L
        const F lambda11 = ( diag==UNIT ? 1 : LBuffer[j+j*ldl] );

        // a10 := a10 / lambda11
        F* a10 = &ABuffer[j];
        if( diag != UNIT )
            for( Int k=0; k<j; ++k )
                a10[k*lda] /= lambda11;

        // A20 := A20 - l21 a10
        F* A20 = &ABuffer[j+1];
        const F* l21 = &LBuffer[(j+1)+j*ldl];
        blas::Geru( a21Height, j, F(-1), l21, 1, a10, lda, A20, lda );

        // alpha11 := alpha11 / |lambda11|^2
        ABuffer[j+j*lda] /= lambda11*Conj(lambda11);
        const F alpha11 = ABuffer[j+j*lda];

        // a21 := a21 / conj(lambda11)
        F* a21 = &ABuffer[(j+1)+j*lda];
        if( diag != UNIT )
            for( Int k=0; k<a21Height; ++k )
                a21[k] /= Conj(lambda11);

        // a21 := a21 - (alpha11/2)l21
        for( Int k=0; k<a21Height; ++k )
            a21[k] -= (alpha11/F(2))*l21[k];

        // A22 := A22 - (l21 a21' + a21 l21')
        F* A22 = &ABuffer[(j+1)+(j+1)*lda];
        blas::Her2( 'L', a21Height, F(-1), l21, 1, a21, 1, A22, lda );

        // a21 := a21 - (alpha11/2)l21
        for( Int k=0; k<a21Height; ++k )
            a21[k] -= (alpha11/F(2))*l21[k];
    }
}

template<typename F>
void UUnb( UnitOrNonUnit diag, Matrix<F>& A, const Matrix<F>& U )
{
    DEBUG_CSE
    // Use the Variant 4 algorithm
    // (which annoyingly requires conjugations for the Her2)
    const Int n = A.Height();
    const Int lda = A.LDim();
    const Int ldu = U.LDim();
    F* ABuffer = A.Buffer();
    const F* UBuffer = U.LockedBuffer();
    vector<F> a12Conj( n ), u12Conj( n );
    for( Int j=0; j<n; ++j )
    {
        const Int a21Height = n - (j+1);

        // Extract and store the diagonal value of U
        const F upsilon11 = ( diag==UNIT ? 1 : UBuffer[j+j*ldu] );

        // a01 := a01 / upsilon11
        F* a01 = &ABuffer[j*lda];
        if( diag != UNIT )
            for( Int k=0; k<j; ++k )
                a01[k] /= upsilon11;

        // A02 := A02 - a01 u12
        F* A02 = &ABuffer[(j+1)*lda];
        const F* u12 = &UBuffer[j+(j+1)*ldu];
        blas::Geru( j, a21Height, F(-1), a01, 1, u12, ldu, A02, lda );

        // alpha11 := alpha11 / |upsilon11|^2
        ABuffer[j+j*lda] /= upsilon11*Conj(upsilon11);
        const F alpha11 = ABuffer[j+j*lda];

        // a12 := a12 / conj(upsilon11)
        F* a12 = &ABuffer[j+(j+1)*lda];
        if( diag != UNIT )
            for( Int k=0; k<a21Height; ++k )
                a12[k*lda] /= Conj(upsilon11);

        // a12 := a12 - (alpha11/2)u12
        for( Int k=0; k<a21Height; ++k )
            a12[k*lda] -= (alpha11/F(2))*u12[k*ldu];

        // A22 := A22 - (a12' u12 + u12' a12)
        F* A22 = &ABuffer[(j+1)+(j+1)*lda];
        for( Int k=0; k<a21Height; ++k )
            a12Conj[k] = Conj(a12[k*lda]);
        for( Int k=0; k<a21Height; ++k )
            u12Conj[k] = Conj(u12[k*ldu]);
        blas::Her2
        ( 'U', a21Height,
          F(-1), u12Conj.data(), 1, a12Conj.data(), 1, A22, lda );

        // a12 := a12 - (alpha11/2)u12
        for( Int k=0; k<a21Height; ++k )
            a12[k*lda] -= (alpha11/F(2))*u12[k*ldu];
    }
}

} // namespace twotrsm
} // namespace El

#endif // ifndef EL_TWOSIDEDTRSM_UNBLOCKED_HPP
