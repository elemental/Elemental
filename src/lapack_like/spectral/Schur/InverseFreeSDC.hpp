/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SCHUR_INVERSEFREESDC_HPP
#define EL_SCHUR_INVERSEFREESDC_HPP

// See Z. Bai, J. Demmel, J. Dongarra, A. Petitet, H. Robinson, and K. Stanley's
// "The spectral decomposition of nonsymmetric matrices on distributed memory
// parallel computers". Currently available at:
// www.netlib.org/lapack/lawnspdf/lawn91.pdf

namespace El {
namespace schur {

//
// X = [B;A]
//

template<typename F>
int InverseFreeSign( Matrix<F>& X, Int maxIts=100, Base<F> tau=0 )
{
    EL_DEBUG_CSE
    typedef Base<F> Real;
    const Int n = X.Width();
    if( X.Height() != 2*n )
        LogicError("X must be 2n x n");
    // Compute the tolerance if it is unset
    if( tau == Real(0) )
        tau = n*limits::Epsilon<Real>();

    // Expose A and B in the original and temporary
    Matrix<F> XAlt( 2*n, n );
    auto B = X( IR(0,n  ), ALL );
    auto A = X( IR(n,2*n), ALL ); 
    auto BAlt = XAlt( IR(0,n  ), ALL );
    auto AAlt = XAlt( IR(n,2*n), ALL );

    // Flip the sign of A
    A *= -1;

    // Set up the space for explicitly computing the left half of Q
    Matrix<F> householderScalars;
    Matrix<Base<F>> signature;
    Matrix<F> Q( 2*n, n );
    auto Q12 = Q( IR(0,n  ), ALL );
    auto Q22 = Q( IR(n,2*n), ALL );

    // Run the iterative algorithm
    Int numIts=0;
    Matrix<F> R, RLast;
    while( numIts < maxIts )
    {
        XAlt = X;
        QR( XAlt, householderScalars, signature );
 
        // Form the left half of Q
        Zero( Q12 );
        MakeIdentity( Q22 );
        qr::ApplyQ( LEFT, NORMAL, XAlt, householderScalars, signature, Q );

        // Save a copy of R
        R = BAlt;
        MakeTrapezoidal( UPPER, R );
        
        // Form the new iterate
        Gemm( ADJOINT, NORMAL, F(1), Q12, A, F(0), AAlt );
        Gemm( ADJOINT, NORMAL, F(1), Q22, B, F(0), BAlt );
        X = XAlt;

        // Use the difference in the iterates to test for convergence
        ++numIts;
        if( numIts > 1 )
        {
            const Real oneRLast = OneNorm(RLast);     
            AxpyTrapezoid( UPPER, F(-1), R, RLast );
            const Real oneRDiff = OneNorm(RLast);
            if( oneRDiff <= tau*oneRLast )
                break;
        }
        RLast = R;
    }

    // Revert the sign of A and return
    A *= -1;
    return numIts;
}

template<typename F>
int InverseFreeSign
( AbstractDistMatrix<F>& XPre, Int maxIts=100, Base<F> tau=0 )
{
    EL_DEBUG_CSE

    DistMatrixReadWriteProxy<F,F,MC,MR> XProx( XPre );
    auto& X = XProx.Get();

    typedef Base<F> Real;
    const Grid& g = X.Grid();
    const Int n = X.Width();
    if( X.Height() != 2*n )
        LogicError("X must be 2n x n");
    // Compute the tolerance if it is unset
    if( tau == Real(0) )
        tau = n*limits::Epsilon<Real>();

    // Expose A and B in the original and temporary
    DistMatrix<F> XAlt( 2*n, n, g );
    auto B = X( IR(0,n  ), ALL );
    auto A = X( IR(n,2*n), ALL );
    auto BAlt = XAlt( IR(0,n  ), ALL );
    auto AAlt = XAlt( IR(n,2*n), ALL );

    // Flip the sign of A
    A *= -1;

    // Set up the space for explicitly computing the left half of Q
    DistMatrix<F,MD,STAR> householderScalars(g);
    DistMatrix<Base<F>,MD,STAR> signature(g);
    DistMatrix<F> Q( 2*n, n, g );
    auto Q12 = Q( IR(0,n  ), ALL );
    auto Q22 = Q( IR(n,2*n), ALL );

    // Run the iterative algorithm
    Int numIts=0;
    DistMatrix<F> R(g), RLast(g);
    while( numIts < maxIts )
    {
        XAlt = X;
        QR( XAlt, householderScalars, signature );
 
        // Form the left half of Q
        Zero( Q12 );
        MakeIdentity( Q22 );
        qr::ApplyQ( LEFT, NORMAL, XAlt, householderScalars, signature, Q );

        // Save a copy of R
        R = BAlt;
        MakeTrapezoidal( UPPER, R );
        
        // Form the new iterate
        Gemm( ADJOINT, NORMAL, F(1), Q12, A, F(0), AAlt );
        Gemm( ADJOINT, NORMAL, F(1), Q22, B, F(0), BAlt );
        X = XAlt;

        // Use the difference in the iterates to test for convergence
        ++numIts;
        if( numIts > 1 )
        {
            const Real oneRLast = OneNorm(RLast);     
            AxpyTrapezoid( UPPER, F(-1), R, RLast );
            const Real oneRDiff = OneNorm(RLast);
            if( oneRDiff <= tau*oneRLast )
                break;
        }
        RLast = R;
    }

    // Revert the sign of A and return
    A *= -1;
    return numIts;
}

template<typename F>
Base<F> InverseFreeSignDivide( Matrix<F>& X )
{
    EL_DEBUG_CSE
    typedef Base<F> Real;
    const Int n = X.Width();
    if( X.Height() != 2*n )
        LogicError("Matrix should be 2n x n");
   
    // Expose A and B, and then copy A
    auto B = X( IR(0,n  ), ALL );
    auto A = X( IR(n,2*n), ALL );
    Matrix<F> ACopy( A );

    // Run the inverse-free alternative to Sign
    InverseFreeSign( X );

    // Compute the pivoted QR decomp of inv(A + B) A [See LAWN91]
    // 1) B := A + B 
    // 2) [Q,R,Pi] := QRP(A)
    // 3) B := Q^H B
    // 4) [R,Q] := RQ(B)
    B += A;
    Matrix<F> householderScalars;
    Matrix<Base<F>> signature;
    Permutation perm;
    QR( A, householderScalars, signature, perm );
    qr::ApplyQ( LEFT, ADJOINT, A, householderScalars, signature, B );
    RQ( B, householderScalars, signature );

    // A := Q^H A Q
    A = ACopy;
    rq::ApplyQ( LEFT, ADJOINT, B, householderScalars, signature, A );
    rq::ApplyQ( RIGHT, NORMAL, B, householderScalars, signature, A );

    // Return || E21 ||1 / || A ||1
    ValueInt<Real> part = ComputePartition( A );
    part.value /= OneNorm(ACopy);
    return part;
}

template<typename F>
ValueInt<Base<F>> InverseFreeSignDivide( AbstractDistMatrix<F>& XPre )
{
    EL_DEBUG_CSE

    DistMatrixReadWriteProxy<F,F,MC,MR> XProx( XPre );
    auto& X = XProx.Get();

    typedef Base<F> Real;
    const Grid& g = X.Grid();
    const Int n = X.Width();
    if( X.Height() != 2*n )
        LogicError("Matrix should be 2n x n");
   
    // Expose A and B, and then copy A
    auto B = X( IR(0,n  ), ALL );
    auto A = X( IR(n,2*n), ALL );
    DistMatrix<F> ACopy( A );

    // Run the inverse-free alternative to Sign
    InverseFreeSign( X );

    // Compute the pivoted QR decomp of inv(A + B) A [See LAWN91]
    // 1) B := A + B 
    // 2) [Q,R,Pi] := QRP(A)
    // 3) B := Q^H B
    // 4) [R,Q] := RQ(B)
    B += A;
    DistMatrix<F,MD,STAR> householderScalars(g);
    DistMatrix<Base<F>,MD,STAR> signature(g);
    DistPermutation perm(g);
    QR( A, householderScalars, signature, perm );
    qr::ApplyQ( LEFT, ADJOINT, A, householderScalars, signature, B );
    RQ( B, householderScalars, signature );

    // A := Q^H A Q
    A = ACopy;
    rq::ApplyQ( LEFT, ADJOINT, B, householderScalars, signature, A );
    rq::ApplyQ( RIGHT, NORMAL, B, householderScalars, signature, A );

    // Return || E21 ||1 / || A ||1
    // Return || E21 ||1 / || A ||1
    ValueInt<Real> part = ComputePartition( A );
    part.value /= OneNorm(ACopy);
    return part;
}

// TODO: RandomizedInverseFreeSignDivide which uses high-level algorithm
//       of http://parlab.eecs.berkeley.edu/sites/all/parlab/files/Communication%20Avoiding%20Nonsymmetric.pdf

} // namespace schur
} // namespace El

#endif // ifndef EL_SCHUR_INVERSEFREESDC_HPP
