/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_HESSENBERG_LOWER_UNBLOCKED_HPP
#define EL_HESSENBERG_LOWER_UNBLOCKED_HPP

namespace El {
namespace hessenberg {

template<typename F>
void LowerUnblocked( Matrix<F>& A, Matrix<F>& householderScalars )
{
    EL_DEBUG_CSE
    const Int n = A.Height();
    const Int householderScalarsHeight = Max(n-1,0);
    householderScalars.Resize( householderScalarsHeight, 1 );

    Matrix<F> z1, z21;

    for( Int k=0; k<n-1; ++k )
    {
        const Range<Int> ind1( k,   k+1 ),
                         ind2( k+1, n   );

        auto a12 = A( ind1, ind2    );
        auto A22 = A( ind2, ind2    );
        auto A2  = A( ind2, IR(0,n) );

        auto alpha12L = A( ind1, IR(k+1,k+2) );
        auto a12R     = A( ind1, IR(k+2,n)   );

        // Find tau and v such that
        //  |alpha12L a12R| /I - tauP | 1   | | 1 conj(v) |\ = |beta 0|
        //                  \         | v^T |              /
        const F tau = RightReflector( alpha12L, a12R );
        householderScalars(k) = tau;

        // Temporarily set a12 = | 1 v |
        const F beta = alpha12L(0);
        alpha12L(0) = F(1);

        // A2 := Hous(a12^T,tau)^H A2
        //     = (I - conj(tau) a12^T conj(a12)) A2
        //     = A2 - conj(tau) a12^T (A2^H a12^T)^H
        // -----------------------------------------
        // z1 := A2^H a12^T
        Zeros( z1, n, 1 );
        Gemv( ADJOINT, F(1), A2, a12, F(0), z1 );
        // A2 := A2 - conj(tau) a12^T z1^H
        Ger( -Conj(tau), a12, z1, A2 );

        // A22 := A22 Hous(a12^T,tau)
        //      = A22 (I - tau a12^T conj(a12))
        //      = A22 - tau (A22 a12^T) conj(a12)
        // --------------------------------------
        // z21 := A22 a12^T
        Zeros( z21, A22.Height(), 1 );
        Gemv( NORMAL, F(1), A22, a12, F(0), z21 );
        // A22 := A22 - tau z21 conj(a12)
        Ger( -tau, z21, a12, A22 );

        // Put beta back
        alpha12L(0) = beta;
    }
}

template<typename F>
void LowerUnblocked
( AbstractDistMatrix<F>& APre,
  AbstractDistMatrix<F>& householderScalarsPre )
{
    EL_DEBUG_CSE

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    DistMatrixWriteProxy<F,F,STAR,STAR>
      householderScalarsProx( householderScalarsPre );
    auto& A = AProx.Get();
    auto& householderScalars = householderScalarsProx.Get();

    const Grid& g = A.Grid();
    const Int n = A.Height();
    const Int householderScalarsHeight = Max(n-1,0);
    householderScalars.SetGrid( g );
    householderScalars.Resize( householderScalarsHeight, 1 );

    DistMatrix<F,MC,  STAR> z21_MC(g);
    DistMatrix<F,MR,  STAR> z1_MR(g);
    DistMatrix<F,STAR,MC  > a12_MC(g);
    DistMatrix<F,STAR,MR  > a12_MR(g);

    for( Int k=0; k<n-1; ++k )
    {
        const Range<Int> ind1( k,   k+1 ),
                         ind2( k+1, n   );

        auto a12 = A( ind1, ind2    );
        auto A22 = A( ind2, ind2    );
        auto A2  = A( ind2, IR(0,n) );

        auto alpha12L = A( ind1, IR(k+1,k+2) );
        auto a12R     = A( ind1, IR(k+2,n)   );

        // Find tau and v such that
        //  |alpha12L a12R| /I - tauP | 1   | | 1 conj(v) |\ = |beta 0|
        //                  \         | v^T |              /
        const F tau = RightReflector( alpha12L, a12R );
        householderScalars.Set(k,0,tau);

        // Temporarily set a12 = | 1 v |
        F beta = 0;
        if( alpha12L.IsLocal(0,0) )
            beta = alpha12L.GetLocal(0,0);
        alpha12L.Set(0,0,F(1));

        // A2 := Hous(a12^T,tau)^H A2
        //     = (I - conj(tau) a12^T conj(a12)) A2
        //     = A2 - conj(tau) a12^T (A2^H a12^T)^H
        // -----------------------------------------
        // z1 := A2^H a12^T
        a12_MC.AlignWith( A2 );
        a12_MC = a12;
        z1_MR.AlignWith( A2 );
        Zeros( z1_MR, n, 1 );
        LocalGemv( ADJOINT, F(1), A2, a12_MC, F(0), z1_MR );
        El::AllReduce( z1_MR, A2.ColComm() );
        // A2 := A2 - conj(tau) a12^T z1^H
        LocalGer( -Conj(tau), a12_MC, z1_MR, A2 );

        // A22 := A22 Hous(a12^T,tau)
        //      = A22 (I - tau a12^T conj(a12))
        //      = A22 - tau (A22 a12^T) conj(a12)
        // --------------------------------------
        // z21 := A22 a12^T
        a12_MR.AlignWith( A22 );
        a12_MR = a12;
        z21_MC.AlignWith( A22 );
        Zeros( z21_MC, A22.Height(), 1 );
        LocalGemv( NORMAL, F(1), A22, a12_MR, F(0), z21_MC );
        El::AllReduce( z21_MC, A22.RowComm() );
        // A22 := A22 - tau z21 conj(a12)
        LocalGer( -tau, z21_MC, a12_MR, A22 );

        // Put beta back
        if( alpha12L.IsLocal(0,0) )
            alpha12L.SetLocal(0,0,beta);
    }
}

} // namespace hessenberg
} // namespace El

#endif // ifndef EL_HESSENBERG_LOWER_UNBLOCKED_HPP
