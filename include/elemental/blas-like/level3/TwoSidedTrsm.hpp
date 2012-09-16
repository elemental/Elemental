/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/

namespace elem {
namespace internal {

template<typename F>
inline void 
TwoSidedTrsmLUnb( UnitOrNonUnit diag, Matrix<F>& A, const Matrix<F>& L )
{
#ifndef RELEASE
    PushCallStack("internal::TwoSidedTrsmLUnb");
#endif
    // Use the Variant 4 algorithm
    const int n = A.Height();
    const int lda = A.LDim();
    const int ldl = L.LDim();
    F* ABuffer = A.Buffer();
    const F* LBuffer = L.LockedBuffer();
    for( int j=0; j<n; ++j )
    {
        const int a21Height = n - (j+1);

        // Extract and store the diagonal value of L
        const F lambda11 = ( diag==UNIT ? 1 : LBuffer[j+j*ldl] );

        // a10 := a10 / lambda11
        F* a10 = &ABuffer[j];
        if( diag != UNIT )
            for( int k=0; k<j; ++k )
                a10[k*lda] /= lambda11;

        // A20 := A20 - l21 a10
        F* A20 = &ABuffer[j+1];
        const F* l21 = &LBuffer[(j+1)+j*ldl];
        blas::Geru( a21Height, j, (F)-1, l21, 1, a10, lda, A20, lda );

        // alpha11 := alpha11 / |lambda11|^2
        ABuffer[j+j*lda] /= lambda11*Conj(lambda11);
        const F alpha11 = ABuffer[j+j*lda];

        // a21 := a21 / conj(lambda11)
        F* a21 = &ABuffer[(j+1)+j*lda];
        if( diag != UNIT )
            for( int k=0; k<a21Height; ++k )
                a21[k] /= Conj(lambda11);

        // a21 := a21 - (alpha11/2)l21
        for( int k=0; k<a21Height; ++k )
            a21[k] -= (alpha11/2)*l21[k];

        // A22 := A22 - (l21 a21' + a21 l21')
        F* A22 = &ABuffer[(j+1)+(j+1)*lda];
        blas::Her2( 'L', a21Height, (F)-1, l21, 1, a21, 1, A22, lda );

        // a21 := a21 - (alpha11/2)l21
        for( int k=0; k<a21Height; ++k )
            a21[k] -= (alpha11/2)*l21[k];
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void 
TwoSidedTrsmUUnb( UnitOrNonUnit diag, Matrix<F>& A, const Matrix<F>& U )
{
#ifndef RELEASE
    PushCallStack("internal::TwoSidedTrsmUUnb");
#endif
    // Use the Variant 4 algorithm
    // (which annoyingly requires conjugations for the Her2)
    const int n = A.Height();
    const int lda = A.LDim();
    const int ldu = U.LDim();
    F* ABuffer = A.Buffer();
    const F* UBuffer = U.LockedBuffer();
    std::vector<F> a12Conj( n ), u12Conj( n );
    for( int j=0; j<n; ++j )
    {
        const int a21Height = n - (j+1);

        // Extract and store the diagonal value of U
        const F upsilon11 = ( diag==UNIT ? 1 : UBuffer[j+j*ldu] );

        // a01 := a01 / upsilon11
        F* a01 = &ABuffer[j*lda];
        if( diag != UNIT )
            for( int k=0; k<j; ++k )
                a01[k] /= upsilon11;

        // A02 := A02 - a01 u12
        F* A02 = &ABuffer[(j+1)*lda];
        const F* u12 = &UBuffer[j+(j+1)*ldu];
        blas::Geru( j, a21Height, (F)-1, a01, 1, u12, ldu, A02, lda );

        // alpha11 := alpha11 / |upsilon11|^2
        ABuffer[j+j*lda] /= upsilon11*Conj(upsilon11);
        const F alpha11 = ABuffer[j+j*lda];

        // a12 := a12 / conj(upsilon11)
        F* a12 = &ABuffer[j+(j+1)*lda];
        if( diag != UNIT )
            for( int k=0; k<a21Height; ++k )
                a12[k*lda] /= Conj(upsilon11);

        // a12 := a12 - (alpha11/2)u12
        for( int k=0; k<a21Height; ++k )
            a12[k*lda] -= (alpha11/2)*u12[k*ldu];

        // A22 := A22 - (a12' u12 + u12' a12)
        F* A22 = &ABuffer[(j+1)+(j+1)*lda];
        for( int k=0; k<a21Height; ++k )
            a12Conj[k] = Conj(a12[k*lda]);
        for( int k=0; k<a21Height; ++k )
            u12Conj[k] = Conj(u12[k*ldu]);
        blas::Her2
        ( 'U', a21Height, (F)-1, &u12Conj[0], 1, &a12Conj[0], 1, A22, lda );

        // a12 := a12 - (alpha11/2)u12
        for( int k=0; k<a21Height; ++k )
            a12[k*lda] -= (alpha11/2)*u12[k*ldu];
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
LocalTwoSidedTrsm
( UpperOrLower uplo, UnitOrNonUnit diag, 
  DistMatrix<F,STAR,STAR>& A, const DistMatrix<F,STAR,STAR>& B )
{
#ifndef RELEASE
    PushCallStack("internal::LocalTwoSidedTrsm");
#endif
    TwoSidedTrsm( uplo, diag, A.LocalMatrix(), B.LockedLocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal
} // namespace elem

#include "./TwoSidedTrsm/LVar1.hpp"
#include "./TwoSidedTrsm/LVar2.hpp"
#include "./TwoSidedTrsm/LVar3.hpp"
#include "./TwoSidedTrsm/LVar4.hpp"
#include "./TwoSidedTrsm/LVar5.hpp"
#include "./TwoSidedTrsm/UVar1.hpp"
#include "./TwoSidedTrsm/UVar2.hpp"
#include "./TwoSidedTrsm/UVar3.hpp"
#include "./TwoSidedTrsm/UVar4.hpp"
#include "./TwoSidedTrsm/UVar5.hpp"

namespace elem {

template<typename F> 
inline void
TwoSidedTrsm
( UpperOrLower uplo, UnitOrNonUnit diag, Matrix<F>& A, const Matrix<F>& B )
{
#ifndef RELEASE
    PushCallStack("TwoSidedTrsm");
#endif
    if( uplo == LOWER )
        internal::TwoSidedTrsmLVar4( diag, A, B );
    else
        internal::TwoSidedTrsmUVar4( diag, A, B );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> 
inline void
TwoSidedTrsm
( UpperOrLower uplo, UnitOrNonUnit diag, 
  DistMatrix<F>& A, const DistMatrix<F>& B )
{
#ifndef RELEASE
    PushCallStack("TwoSidedTrsm");
#endif
    if( uplo == LOWER )
        internal::TwoSidedTrsmLVar4( diag, A, B );
    else
        internal::TwoSidedTrsmUVar4( diag, A, B );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
