/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BIDIAG_UPPER_UNBLOCKED_HPP
#define EL_BIDIAG_UPPER_UNBLOCKED_HPP

namespace El {
namespace bidiag {

template<typename F>
void UpperUnblocked
( Matrix<F>& A,
  Matrix<F>& householderScalarsP,
  Matrix<F>& householderScalarsQ )
{
    EL_DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    const Int householderScalarsPHeight = Max(n-1,0);
    const Int householderScalarsQHeight = n;
    EL_DEBUG_ONLY(
      if( m < n )
          LogicError("A must be at least as tall as it is wide");
    )
    householderScalarsP.Resize( householderScalarsPHeight, 1 );
    householderScalarsQ.Resize( householderScalarsQHeight, 1 );

    Matrix<F> x12Adj, w21;

    for( Int k=0; k<n; ++k )
    {
        auto alpha11 = A( IR(k),     IR(k)     );
        auto a12     = A( IR(k),     IR(k+1,n) );
        auto a21     = A( IR(k+1,m), IR(k)     );
        auto A22     = A( IR(k+1,m), IR(k+1,n) );
        auto aB1     = A( IR(k,m),   IR(k)     );
        auto AB2     = A( IR(k,m),   IR(k+1,n) );

        // Find tauQ and u such that
        //  / I - tauQ | 1 | | 1, u^H | \ | alpha11 | = | epsilonQ |
        //  \          | u |            / |     a21 | = |    0     |
        const F tauQ = LeftReflector( alpha11, a21 );
        householderScalarsQ(k) = tauQ;

        // Temporarily set aB1 = | 1 |
        //                       | u |
        const F epsilonQ = alpha11(0);
        alpha11(0) = F(1);

        // AB2 := Hous(aB1,tauQ) AB2
        //      = (I - tauQ aB1 aB1^H) AB2
        //      = AB2 - tauQ aB1 (AB2^H aB1)^H
        // -----------------------------------
        // x12^H := (aB1^H AB2)^H = AB2^H aB1
        Zeros( x12Adj, a12.Width(), 1 );
        Gemv( ADJOINT, F(1), AB2, aB1, F(0), x12Adj );
        // AB2 := AB2 - tauQ aB1 x12
        //      = (I - tauQ aB1 aB1^H) AB2
        Ger( -tauQ, aB1, x12Adj, AB2 );

        // Put epsilonQ back 
        alpha11(0) = epsilonQ;

        if( k+1 < n )
        {
            // Expose the subvector we seek to zero, a12R
            Matrix<F> alpha12L, a12R;
            PartitionRight( a12, alpha12L, a12R, 1 );

            // Find tauP and v such that
            //  |alpha12L a12R| /I - tauP |1  | |1, conj(v)|\ = |epsilonP 0|
            //                  \         |v^T|             /
            const F tauP = RightReflector( alpha12L, a12R );
            householderScalarsP(k) = tauP;

            // Temporarily set a12^T = | 1 | 
            //                         | v |
            const F epsilonP = alpha12L(0);
            alpha12L(0) = F(1);

            // A22 := A22 Hous(a12^T,tauP)
            //      = A22 (I - tauP a12^T conj(a12))
            //      = A22 - tauP (A22 a12^T) conj(a12)
            // ---------------------------------------
            // w21 := A22 a12^T = A22 | 1   |
            //                        | v^T |
            Zeros( w21, a21.Height(), 1 );
            Gemv( NORMAL, F(1), A22, a12, F(0), w21 );
            // A22 := A22 - tauP w21 conj(a12)
            Ger( -tauP, w21, a12, A22 );

            // Put epsilonP back 
            alpha12L(0) = epsilonP;
        }
    }
}

template<typename F> 
void UpperUnblocked
( AbstractDistMatrix<F>& APre, 
  AbstractDistMatrix<F>& householderScalarsPPre,
  AbstractDistMatrix<F>& householderScalarsQPre )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      AssertSameGrids( APre, householderScalarsPPre, householderScalarsQPre )
    )
    DistMatrixReadWriteProxy<F,F,MC,MR>
      AProx( APre );
    DistMatrixWriteProxy<F,F,STAR,STAR>
      householderScalarsPProx( householderScalarsPPre ),
      householderScalarsQProx( householderScalarsQPre );
    auto& A = AProx.Get();
    auto& householderScalarsP = householderScalarsPProx.Get();
    auto& householderScalarsQ = householderScalarsQProx.Get();

    const Int m = A.Height();
    const Int n = A.Width();
    EL_DEBUG_ONLY(
      if( m < n )
          LogicError("A must be at least as tall as it is wide");
    )
    const Grid& g = A.Grid();
    const Int householderScalarsPHeight = Max(n-1,0);
    const Int householderScalarsQHeight = n;
    householderScalarsP.Resize( householderScalarsPHeight, 1 );
    householderScalarsQ.Resize( householderScalarsQHeight, 1 );

    DistMatrix<F,STAR,MR  > a12_STAR_MR(g);
    DistMatrix<F,MC,  STAR> aB1_MC_STAR(g);
    DistMatrix<F,MR,  STAR> x12Adj_MR_STAR(g);
    DistMatrix<F,MC,  STAR> w21_MC_STAR(g);

    for( Int k=0; k<n; ++k )
    {
        auto alpha11 = A( IR(k),     IR(k)     );
        auto a12     = A( IR(k),     IR(k+1,n) );
        auto a21     = A( IR(k+1,m), IR(k)     );
        auto A22     = A( IR(k+1,m), IR(k+1,n) );
        auto aB1     = A( IR(k,m),   IR(k)     );
        auto AB2     = A( IR(k,m),   IR(k+1,n) );

        // Find tauQ and u such that
        //  / I - tauQ | 1 | | 1, u^H | \ | alpha11 | = | epsilonQ |
        //  \          | u |            / |     a21 | = |    0     |
        const F tauQ = LeftReflector( alpha11, a21 );
        householderScalarsQ.Set(k,0,tauQ );

        // Temporarily set aB1 = | 1 |
        //                       | u |
        F epsilonQ=0;
        if( alpha11.IsLocal(0,0) )
            epsilonQ = alpha11.GetLocal(0,0);
        alpha11.Set(0,0,F(1));

        // AB2 := Hous(aB1,tauQ) AB2
        //      = (I - tauQ aB1 aB1^H) AB2
        //      = AB2 - tauQ aB1 (AB2^H aB1)^H
        // ----------------------------------
        // x12^H := (aB1^H AB2)^H = AB2^H aB1
        aB1_MC_STAR.AlignWith( aB1 );
        aB1_MC_STAR = aB1;
        x12Adj_MR_STAR.AlignWith( AB2 );
        Zeros( x12Adj_MR_STAR, a12.Width(), 1 );
        LocalGemv( ADJOINT, F(1), AB2, aB1_MC_STAR, F(0), x12Adj_MR_STAR );
        El::AllReduce( x12Adj_MR_STAR, AB2.ColComm() );
        // AB2 := AB2 - tauQ aB1 x12
        LocalGer( -tauQ, aB1_MC_STAR, x12Adj_MR_STAR, AB2 );

        // Put epsilonQ back 
        if( alpha11.IsLocal(0,0) )
            alpha11.SetLocal(0,0,epsilonQ);

        if( k+1 < n )
        {
            // Expose the subvector we seek to zero, a12R
            DistMatrix<F> alpha12L(g), a12R(g);
            PartitionRight( a12, alpha12L, a12R, 1 );

            // Find tauP and v such that
            //  |alpha12L a12R| /I - tauP |1  | |1, conj(v)|\ = |epsilonP 0|
            //                  \         |v^T|             /
            const F tauP = RightReflector( alpha12L, a12R );
            householderScalarsP.Set(k,0,tauP);

            // Temporarily set a12^T = | 1   |
            //                         | v^T |
            F epsilonP=0;
            if( alpha12L.IsLocal(0,0) )
                epsilonP = alpha12L.GetLocal(0,0);
            alpha12L.Set(0,0,F(1));

            // A22 := A22 Hous(a12^T,tauP)
            //      = A22 (I - tauP a12^T conj(a12))
            //      = A22 - tauP (A22 a12^T) conj(a12)
            // -------------------------------------
            // w21 := A22 a12^T = A22 | 1   |
            //                        | v^T |
            a12_STAR_MR.AlignWith( a12 );
            a12_STAR_MR = a12;
            w21_MC_STAR.AlignWith( A22 );
            Zeros( w21_MC_STAR, a21.Height(), 1 );
            LocalGemv( NORMAL, F(1), A22, a12_STAR_MR, F(0), w21_MC_STAR );
            El::AllReduce( w21_MC_STAR, A22.RowComm() );
            // A22 := A22 - tauP w21 conj(a12)
            LocalGer( -tauP, w21_MC_STAR, a12_STAR_MR, A22 );

            // Put epsilonP back
            if( alpha12L.IsLocal(0,0) )
                alpha12L.SetLocal(0,0,epsilonP);
        }
    }
}

} // namespace bidiag
} // namespace El

#endif // ifndef EL_BIDIAG_UPPER_UNBLOCKED_HPP
