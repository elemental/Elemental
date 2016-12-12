/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_CHOLESKY_UPPER_MOD_HPP
#define EL_CHOLESKY_UPPER_MOD_HPP

// TODO: Blocked algorithms

namespace El {
namespace cholesky {

namespace mod {

template<typename F>
void UpperUpdate( Matrix<F>& U, Matrix<F>& V )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( U.Height() != U.Width() )
          LogicError("Cholesky factors must be square");
      if( V.Height() != U.Height() )
          LogicError("V is the wrong height");
    )
    const Int m = V.Height();

    Matrix<F> z12;

    for( Int k=0; k<m; ++k )
    {
        const IR ind1( k ), ind2( k+1, END );

        F& upsilon11 = U(k,k);
        auto u12 = U( ind1, ind2 );

        auto v1 = V( ind1, ALL );
        auto V2 = V( ind2, ALL );

        // Find tau and w such that
        //  /I - tau | 1 | | 1 w^H | \ | upsilon11 | = | -beta |
        //  \        | w |           / | v1^H      |   | 0     |
        // where beta >= 0
        Conjugate( v1 );
        const F tau = LeftReflector( upsilon11, v1 );

        // Apply the negative of the Householder reflector from the left:
        // | u12  | := -| u12  | + tau | 1 | | 1 w^H | | u12  |
        // | V2^H |     | V2^H |       | w |           | V2^H |
        //
        //           = -| u12  | + tau | 1 | (u12 + w^H V2^H)
        //              | V2^H |       | w |
        // Thus,
        //    u12 := -u12 + tau (u12 + w^H V2^H), and
        //    V2^H := -V2^H + tau w (u12 + w^H V2^H),
        //    V2   := -V2 + conj(tau) (u12^H + V2 w) w^H
        //
        upsilon11 = -upsilon11;
        Conjugate( u12, z12 );
        Gemv( NORMAL, F(1), V2, v1, F(1), z12 );
        V2 *= -1;
        Ger( Conj(tau), z12, v1, V2 );
        Conjugate( z12 );
        u12 *= -1;
        Axpy( tau, z12, u12 );
    }
}

template<typename F>
void UpperUpdate
( AbstractDistMatrix<F>& UPre,
  AbstractDistMatrix<F>& VPre )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( UPre.Height() != UPre.Width() )
          LogicError("Cholesky factors must be square");
      if( VPre.Height() != UPre.Height() )
          LogicError("V is the wrong height");
      AssertSameGrids( UPre, VPre );
    )

    DistMatrixReadWriteProxy<F,F,MC,MR> UProx( UPre ), VProx( VPre );
    auto& U = UProx.Get();
    auto& V = VProx.Get();

    const Int m = V.Height();
    const Grid& g = U.Grid();
    DistMatrix<F,STAR,MC> z12_STAR_MC(g), b12_STAR_MC(g);
    DistMatrix<F,STAR,MR> v1_STAR_MR(g);

    for( Int k=0; k<m; ++k )
    {
        const IR ind1( k ), ind2( k+1, END );

        F upsilon11 = U.Get(k,k);
        auto u12 = U( ind1, ind2 );

        auto v1 = V( ind1, ALL );
        auto V2 = V( ind2, ALL );

        // Find tau and w such that
        //  /I - tau | 1 | | 1 w^H | \ | upsilon11 | = | -beta |
        //  \        | w |           / | v1^H      |   | 0     |
        // where beta >= 0
        Conjugate( v1 );
        // LeftReflector assumes that v1 is a column-vector, so conjugate
        // the RightReflector result
        const F tau = RightReflector( upsilon11, v1 );
        Conjugate( v1 );

        // Apply the negative of the Householder reflector from the left:
        // | u12  | := -| u12  | + tau | 1 | | 1 w^H | | u12  |
        // | V2^H |     | V2^H |       | w |           | V2^H |
        //
        //           = -| u12  | + tau | 1 | (u12 + w^H V2^H)
        //              | V2^H |       | w |
        // Thus,
        //    u12 := -u12 + tau (u12 + w^H V2^H), and
        //    V2^H := -V2^H + tau w (u12 + w^H V2^H),
        //    V2   := -V2 + conj(tau) (u12^H + V2 w) w^H
        //
        U.Set( k, k, -upsilon11 );
        z12_STAR_MC.AlignWith( V2 );
        b12_STAR_MC.AlignWith( V2 );
        v1_STAR_MR.AlignWith( V2 );
        Conjugate( u12, z12_STAR_MC );
        v1_STAR_MR = v1;
        Zeros( b12_STAR_MC, 1, V2.Height() );
        LocalGemv( NORMAL, F(1), V2, v1_STAR_MR, F(0), b12_STAR_MC );
        El::AllReduce( b12_STAR_MC, V2.RowComm() );
        z12_STAR_MC += b12_STAR_MC;
        V2 *= -1;
        LocalGer( Conj(tau), z12_STAR_MC, v1_STAR_MR, V2 );
        Conjugate( z12_STAR_MC );
        u12 *= -1;
        Axpy( tau, z12_STAR_MC, u12 );
    }
}

template<typename F>
void UpperDowndate( Matrix<F>& U, Matrix<F>& V )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( U.Height() != U.Width() )
          LogicError("Cholesky factors must be square");
      if( V.Height() != U.Height() )
          LogicError("V is the wrong height");
    )
    const Int m = V.Height();

    Matrix<F> z12;

    for( Int k=0; k<m; ++k )
    {
        const IR ind1( k ), ind2( k+1, END );

        F& upsilon11 = U(k,k);
        auto u12 = U( ind1, ind2 );

        auto v1 = V( ind1, ALL );
        auto V2 = V( ind2, ALL );

        // Find tau and w such that
        //  /I - 1/tau | 1 | | 1 w^H | Sigma \ | upsilon11 | = | -beta |
        //  \          | w |                 / | v1^H      |   | 0     |
        // where Sigma = diag(+1,-1,...,-1) and beta >= 0
        Conjugate( v1 );
        const F tau = LeftHyperbolicReflector( upsilon11, v1 );

        // Apply the negative of the Householder reflector from the left:
        // | u12  | := -| u12  | + 1/tau | 1 | | 1 w^H | Sigma | u12  |
        // | V2^H |     | V2^H |         | w |                 | V2^H |
        //
        //           = -| u12  | + 1/tau | 1 | (u12 - w^H V2^H)
        //              | V2^H |         | w |
        // Thus,
        //    u12 := -u12 + 1/tau (u12 - w^H V2^H), and
        //    V2^H := -V2^H + 1/tau w (u12 - w^H V2^H),
        //    V2   := -V2 + 1/tau (u12^H - V2 w) w^H
        //
        upsilon11 = -upsilon11;
        Conjugate( u12, z12 );
        Gemv( NORMAL, F(-1), V2, v1, F(1), z12 );
        u12 *= -1;
        V2 *= -1;
        Ger( F(1)/tau, z12, v1, V2 );
        Conjugate( z12 );
        Axpy( F(1)/tau, z12, u12 );
    }
}

template<typename F>
void UpperDowndate
( AbstractDistMatrix<F>& UPre, AbstractDistMatrix<F>& VPre )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( UPre.Height() != UPre.Width() )
          LogicError("Cholesky factors must be square");
      if( VPre.Height() != UPre.Height() )
          LogicError("V is the wrong height");
      AssertSameGrids( UPre, VPre );
    )

    DistMatrixReadWriteProxy<F,F,MC,MR> UProx( UPre ), VProx( VPre );
    auto& U = UProx.Get();
    auto& V = VProx.Get();

    const Int m = V.Height();
    const Grid& g = U.Grid();
    DistMatrix<F,STAR,MC> z12_STAR_MC(g), b12_STAR_MC(g);
    DistMatrix<F,STAR,MR> v1_STAR_MR(g);

    for( Int k=0; k<m; ++k )
    {
        const IR ind1( k ), ind2( k+1, END );

        F upsilon11 = U.Get(k,k);
        auto u12 = U( ind1, ind2 );

        auto v1 = V( ind1, ALL );
        auto V2 = V( ind2, ALL );

        // Find tau and w such that
        //  /I - 1/tau | 1 | | 1 w^H | Sigma \ | upsilon11 | = | -beta |
        //  \          | w |                 / | v1^H      |   | 0     |
        // where Sigma = diag(+1,-1,...,-1) and beta >= 0
        // Since LeftHyperbolicReflector expects v1 to be a column vector,
        // use the Right variant and conjugate the result
        Conjugate( v1 );
        const F tau = RightHyperbolicReflector( upsilon11, v1 );
        Conjugate( v1 );

        // Apply the negative of the Householder reflector from the left:
        // | u12  | := -| u12  | + 1/tau | 1 | | 1 w^H | Sigma | u12  |
        // | V2^H |     | V2^H |         | w |                 | V2^H |
        //
        //           = -| u12  | + 1/tau | 1 | (u12 - w^H V2^H)
        //              | V2^H |         | w |
        // Thus,
        //    u12 := -u12 + 1/tau (u12 - w^H V2^H), and
        //    V2^H := -V2^H + 1/tau w (u12 - w^H V2^H),
        //    V2   := -V2 + 1/tau (u12^H - V2 w) w^H
        //
        U.Set( k, k, -upsilon11 );
        z12_STAR_MC.AlignWith( V2 );
        b12_STAR_MC.AlignWith( V2 );
        v1_STAR_MR.AlignWith( V2 );
        Conjugate( u12, z12_STAR_MC );
        v1_STAR_MR = v1;
        Zeros( b12_STAR_MC, 1, V2.Height() );
        LocalGemv( NORMAL, F(-1), V2, v1_STAR_MR, F(0), b12_STAR_MC );
        El::AllReduce( b12_STAR_MC, V2.RowComm() );
        z12_STAR_MC += b12_STAR_MC;
        u12 *= -1;
        V2 *= -1;
        LocalGer( F(1)/tau, z12_STAR_MC, v1_STAR_MR, V2 );
        Conjugate( z12_STAR_MC );
        Axpy( F(1)/tau, z12_STAR_MC, u12 );
    }
}

} // namespace mod

template<typename F>
void UpperMod( Matrix<F>& U, Base<F> alpha, Matrix<F>& V )
{
    EL_DEBUG_CSE
    typedef Base<F> Real;
    if( alpha == Real(0) )
        return;
    else if( alpha > Real(0) )
    {
        V *= Sqrt(alpha);
        mod::UpperUpdate( U, V );
    }
    else
    {
        V *= Sqrt(-alpha);
        mod::UpperDowndate( U, V );
    }
}

template<typename F>
void UpperMod
( AbstractDistMatrix<F>& U,
  Base<F> alpha,
  AbstractDistMatrix<F>& V )
{
    EL_DEBUG_CSE
    typedef Base<F> Real;
    if( alpha == Real(0) )
        return;
    else if( alpha > Real(0) )
    {
        V *= Sqrt(alpha);
        mod::UpperUpdate( U, V );
    }
    else
    {
        V *= Sqrt(-alpha);
        mod::UpperDowndate( U, V );
    }
}

} // namespace cholesky
} // namespace El

#endif // ifndef EL_CHOLESKY_UPPER_MOD_HPP
