/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_POLAR_QDWH_HPP
#define EL_POLAR_QDWH_HPP

#include EL_IDENTITY_INC

namespace El {

// Based on Yuji Nakatsukasa's implementation of a QR-based dynamically 
// weighted Halley iteration for the polar decomposition. In particular, this
// implementation mirrors the routine 'qdwh', which is part of the zip-file
// available here:
//     http://www.mathworks.com/matlabcentral/fileexchange/36830
//
// No support for row-sorting yet.
//
// The careful calculation of the coefficients is due to a suggestion from
// Gregorio Quintana Orti.

namespace polar {

template<typename F>
inline Int 
QDWHInner( Matrix<F>& A, Base<F> sMinUpper, const PolarCtrl& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("polar::QDWHInner"))
    typedef Base<F> Real;
    typedef Complex<Real> Cpx;
    const Int m = A.Height();
    const Int n = A.Width();
    const Real oneThird = Real(1)/Real(3);
    if( m < n )
        LogicError("Height cannot be less than width");

    const Real eps = lapack::MachineEpsilon<Real>();
    const Real tol = 5*eps;
    const Real cubeRootTol = Pow(tol,oneThird);
    Real L = sMinUpper / Sqrt(Real(n));

    Real frobNormADiff;
    Matrix<F> ALast, ATemp, C;
    Matrix<F> Q( m+n, n );
    Matrix<F> QT, QB;
    PartitionDown( Q, QT, QB, m );
    Int numIts=0;
    while( numIts < ctrl.maxIts )
    {
        ALast = A;

        Real L2;
        Cpx dd, sqd;
        if( Abs(1-L) < tol )
        {
            L2 = 1;
            dd = 0;
            sqd = 1;
        }
        else
        {
            L2 = L*L;
            dd = Pow( 4*(1-L2)/(L2*L2), oneThird );
            sqd = Sqrt( Real(1)+dd );
        }
        const Cpx arg = Real(8) - Real(4)*dd + Real(8)*(2-L2)/(L2*sqd);
        const Real a = (sqd + Sqrt(arg)/Real(2)).real();
        const Real b = (a-1)*(a-1)/4;
        const Real c = a+b-1;
        const Real alpha = a-b/c;
        const Real beta = b/c;

        L = L*(a+b*L2)/(1+c*L2);

        if( c > 100 )
        {
            //
            // The standard QR-based algorithm
            //
            QT = A;
            Scale( Sqrt(c), QT );
            MakeIdentity( QB );
            qr::Explicit( Q, ctrl.colPiv );
            Gemm( NORMAL, ADJOINT, F(alpha/Sqrt(c)), QT, QB, F(beta), A );
        }
        else
        {
            //
            // Use faster Cholesky-based algorithm since A is well-conditioned
            //
            Identity( C, n, n );
            Herk( LOWER, ADJOINT, F(c), A, F(1), C );
            Cholesky( LOWER, C );
            ATemp = A;
            Trsm( RIGHT, LOWER, ADJOINT, NON_UNIT, F(1), C, ATemp );
            Trsm( RIGHT, LOWER, NORMAL, NON_UNIT, F(1), C, ATemp );
            Scale( beta, A );
            Axpy( alpha, ATemp, A );
        }

        ++numIts;
        Axpy( F(-1), A, ALast );
        frobNormADiff = FrobeniusNorm( ALast );
        if( frobNormADiff <= cubeRootTol && Abs(1-L) <= tol )
            break;
    }
    return numIts;
}

template<typename F>
inline Int 
QDWH( Matrix<F>& A, const PolarCtrl& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("polar::QDWH"))
    typedef Base<F> Real;
    const Real twoEst = TwoNormEstimate( A );
    Scale( F(1)/twoEst, A );

    // The one-norm of the inverse can be replaced with an estimate which is
    // a few times cheaper, e.g., via Higham and Tisseur's block algorithm
    // from "A Block Algorithm for Matrix 1-Norm Estimation, with an Application
    // to 1-Norm Pseudospectra".
    Real sMinUpper;
    Matrix<F> Y( A );
    if( A.Height() > A.Width() )
    {
        QR( Y );
        Y.Resize( A.Width(), A.Width() );
        MakeTriangular( UPPER, Y );
        try 
        {
            TriangularInverse( UPPER, NON_UNIT, Y );
            sMinUpper = Real(1) / OneNorm( Y );
        } catch( SingularMatrixException& e ) { sMinUpper = 0; }
    }
    else
    {
        try 
        {
            Inverse( Y );
            sMinUpper = Real(1) / OneNorm( Y );
        } catch( SingularMatrixException& e ) { sMinUpper = 0; }
    } 

    return QDWHInner( A, sMinUpper, ctrl );
}

template<typename F>
inline Int 
QDWH( Matrix<F>& A, Matrix<F>& P, const PolarCtrl& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("polar::QDWH"))
    Matrix<F> ACopy( A );
    const Int numIts = QDWH( A, ctrl );
    Zeros( P, A.Height(), A.Height() );
    Trrk( LOWER, NORMAL, NORMAL, F(1), A, ACopy, F(0), P );
    MakeHermitian( LOWER, P );
    return numIts;
}

template<typename F>
inline Int 
QDWHInner
( DistMatrix<F>& A, Base<F> sMinUpper, const PolarCtrl& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("polar::QDWHInner"))
    typedef Base<F> Real;
    typedef Complex<Real> Cpx;
    const Grid& g = A.Grid();
    const Int m = A.Height();
    const Int n = A.Width();
    const Real oneThird = Real(1)/Real(3);
    if( m < n )
        LogicError("Height cannot be less than width");

    const Real eps = lapack::MachineEpsilon<Real>();
    const Real tol = 5*eps;
    const Real cubeRootTol = Pow(tol,oneThird);
    Real L = sMinUpper / Sqrt(Real(n));

    Real frobNormADiff;
    DistMatrix<F> ALast(g), ATemp(g), C(g);
    DistMatrix<F> Q( m+n, n, g );
    DistMatrix<F> QT(g), QB(g);
    PartitionDown( Q, QT, QB, m );
    Int numIts=0;
    while( numIts < ctrl.maxIts )
    {
        ALast = A;

        Real L2;
        Cpx dd, sqd;
        if( Abs(1-L) < tol )
        {
            L2 = 1;
            dd = 0;
            sqd = 1;
        }
        else
        {
            L2 = L*L;
            dd = Pow( 4*(1-L2)/(L2*L2), oneThird );
            sqd = Sqrt( Real(1)+dd );
        }
        const Cpx arg = Real(8) - Real(4)*dd + Real(8)*(2-L2)/(L2*sqd);
        const Real a = (sqd + Sqrt(arg)/Real(2)).real();
        const Real b = (a-1)*(a-1)/4;
        const Real c = a+b-1;
        const Real alpha = a-b/c;
        const Real beta = b/c;

        L = L*(a+b*L2)/(1+c*L2);

        if( c > 100 )
        {
            //
            // The standard QR-based algorithm
            //
            QT = A;
            Scale( Sqrt(c), QT );
            MakeIdentity( QB );
            qr::Explicit( Q, ctrl.colPiv );
            Gemm( NORMAL, ADJOINT, F(alpha/Sqrt(c)), QT, QB, F(beta), A );
        }
        else
        {
            //
            // Use faster Cholesky-based algorithm since A is well-conditioned
            //
            Identity( C, n, n );
            Herk( LOWER, ADJOINT, F(c), A, F(1), C );
            Cholesky( LOWER, C );
            ATemp = A;
            Trsm( RIGHT, LOWER, ADJOINT, NON_UNIT, F(1), C, ATemp );
            Trsm( RIGHT, LOWER, NORMAL, NON_UNIT, F(1), C, ATemp );
            Scale( beta, A );
            Axpy( alpha, ATemp, A );
        }

        ++numIts;
        Axpy( F(-1), A, ALast );
        frobNormADiff = FrobeniusNorm( ALast );
        if( frobNormADiff <= cubeRootTol && Abs(1-L) <= tol )
            break;
    }
    return numIts;
}

template<typename F>
inline Int 
QDWH( DistMatrix<F>& A, const PolarCtrl& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("polar::QDWH"))
    typedef Base<F> Real;
    const Real twoEst = TwoNormEstimate( A );
    Scale( F(1)/twoEst, A );

    // The one-norm of the inverse can be replaced with an estimate which is
    // a few times cheaper, e.g., via Higham and Tisseur's block algorithm
    // from "A Block Algorithm for Matrix 1-Norm Estimation, with an Application
    // to 1-Norm Pseudospectra".
    Real sMinUpper;
    DistMatrix<F> Y( A );
    if( A.Height() > A.Width() )
    {
        QR( Y );
        Y.Resize( A.Width(), A.Width() );
        MakeTriangular( UPPER, Y );
        try
        {
            TriangularInverse( UPPER, NON_UNIT, Y );
            sMinUpper = Real(1) / OneNorm( Y );
        } catch( SingularMatrixException& e ) { sMinUpper = 0; }
    }
    else
    {
        try
        {
            Inverse( Y );
            sMinUpper = Real(1) / OneNorm( Y );
        } catch( SingularMatrixException& e ) { sMinUpper = 0; }
    }

    return QDWHInner( A, sMinUpper, ctrl );
}

template<typename F>
inline Int 
QDWH( DistMatrix<F>& A, DistMatrix<F>& P, const PolarCtrl& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("polar::QDWH"))
    DistMatrix<F> ACopy( A );
    const Int numIts = QDWH( A, ctrl );
    Zeros( P, A.Height(), A.Height() );
    Trrk( LOWER, NORMAL, NORMAL, F(1), A, ACopy, F(0), P );
    MakeHermitian( LOWER, P );
    return numIts;
}

} // namespace polar

namespace herm_polar {

template<typename F>
inline int
QDWHInner
( UpperOrLower uplo, Matrix<F>& A, Base<F> sMinUpper, const PolarCtrl& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("herm_polar::QDWH"))
    if( A.Height() != A.Width() )
        LogicError("Height must be same as width");

    typedef Base<F> Real;
    typedef Complex<Real> Cpx;
    const Int n = A.Height();
    const Real oneThird = Real(1)/Real(3);

    const Real eps = lapack::MachineEpsilon<Real>();
    const Real tol = 5*eps;
    const Real cubeRootTol = Pow(tol,oneThird);
    Real L = sMinUpper / Sqrt(Real(n));

    Real frobNormADiff;
    Matrix<F> ALast, ATemp, C;
    Matrix<F> Q( 2*n, n );
    Matrix<F> QT, QB;
    PartitionDown( Q, QT, QB, n );
    Int numIts=0;
    while( numIts < ctrl.maxIts )
    {
        ALast = A;

        Real L2;
        Cpx dd, sqd;
        if( Abs(1-L) < tol )
        {
            L2 = 1;
            dd = 0;
            sqd = 1;
        }
        else
        {
            L2 = L*L;
            dd = Pow( 4*(1-L2)/(L2*L2), oneThird );
            sqd = Sqrt( Real(1)+dd );
        }
        const Cpx arg = Real(8) - Real(4)*dd + Real(8)*(2-L2)/(L2*sqd);
        const Real a = (sqd + Sqrt(arg)/Real(2)).real();
        const Real b = (a-1)*(a-1)/4;
        const Real c = a+b-1;
        const Real alpha = a-b/c;
        const Real beta = b/c;

        L = L*(a+b*L2)/(1+c*L2);

        if( c > 100 )
        {
            //
            // The standard QR-based algorithm
            //
            MakeHermitian( uplo, A );
            QT = A;
            Scale( Sqrt(c), QT );
            MakeIdentity( QB );
            qr::Explicit( Q, ctrl.colPiv );
            Trrk( uplo, NORMAL, ADJOINT, F(alpha/Sqrt(c)), QT, QB, F(beta), A );
        }
        else
        {
            //
            // Use faster Cholesky-based algorithm since A is well-conditioned
            //
            // TODO: Think of how to better exploit the symmetry of A,
            //       e.g., by halving the work in the first Herk through 
            //       a custom routine for forming L^2, where L is strictly lower
            MakeHermitian( uplo, A );
            Identity( C, n, n );
            Herk( LOWER, ADJOINT, F(c), A, F(1), C );
            Cholesky( LOWER, C );
            ATemp = A;
            Trsm( RIGHT, LOWER, ADJOINT, NON_UNIT, F(1), C, ATemp );
            Trsm( RIGHT, LOWER, NORMAL, NON_UNIT, F(1), C, ATemp );
            Scale( beta, A );
            Axpy( alpha, ATemp, A );
        }

        Axpy( F(-1), A, ALast );
        frobNormADiff = HermitianFrobeniusNorm( uplo, ALast );

        ++numIts;
        if( frobNormADiff <= cubeRootTol && Abs(1-L) <= tol )
            break;
    }

    MakeHermitian( uplo, A );
    return numIts;
}

template<typename F>
inline Int 
QDWH( UpperOrLower uplo, Matrix<F>& A, const PolarCtrl& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("herm_polar::QDWH"))
    typedef Base<F> Real;
    MakeHermitian( uplo, A );
    const Real twoEst = TwoNormEstimate( A );
    Scale( F(1)/twoEst, A );

    // The one-norm of the inverse can be replaced with an estimate which is
    // a few times cheaper, e.g., via Higham and Tisseur's block algorithm
    // from "A Block Algorithm for Matrix 1-Norm Estimation, with an Application
    // to 1-Norm Pseudospectra".
    Real sMinUpper;
    Matrix<F> Y( A );
    try
    {
        Inverse( Y );
        sMinUpper = Real(1) / OneNorm( Y );
    } catch( SingularMatrixException& e ) { sMinUpper = 0; }

    return QDWHInner( uplo, A, sMinUpper, ctrl );
}

template<typename F>
inline Int
QDWH
( UpperOrLower uplo, Matrix<F>& A, Matrix<F>& P, const PolarCtrl& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("herm_polar::QDWH"))
    Matrix<F> ACopy( A );
    // NOTE: This might be avoidable
    MakeHermitian( uplo, ACopy );
    const Int numIts = QDWH( uplo, A, ctrl );
    Zeros( P, A.Height(), A.Height() );
    Trrk( uplo, NORMAL, NORMAL, F(1), A, ACopy, F(0), P );
    return numIts;
}

template<typename F>
inline int
QDWHInner
( UpperOrLower uplo, DistMatrix<F>& A, Base<F> sMinUpper, 
  const PolarCtrl& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("herm_polar::QDWH"))
    if( A.Height() != A.Width() )
        LogicError("Height must be same as width");

    typedef Base<F> Real;
    typedef Complex<Real> Cpx;
    const Grid& g = A.Grid();
    const Int n = A.Height();
    const Real oneThird = Real(1)/Real(3);

    const Real eps = lapack::MachineEpsilon<Real>();
    const Real tol = 5*eps;
    const Real cubeRootTol = Pow(tol,oneThird);
    Real L = sMinUpper / Sqrt(Real(n));

    Real frobNormADiff;
    DistMatrix<F> ALast(g), ATemp(g), C(g);
    DistMatrix<F> Q( 2*n, n, g );
    DistMatrix<F> QT(g), QB(g);
    PartitionDown( Q, QT, QB, n );
    Int numIts=0;
    while( numIts < ctrl.maxIts )
    {
        ALast = A;

        Real L2;
        Cpx dd, sqd;
        if( Abs(1-L) < tol )
        {
            L2 = 1;
            dd = 0;
            sqd = 1;
        }
        else
        {
            L2 = L*L;
            dd = Pow( 4*(1-L2)/(L2*L2), oneThird );
            sqd = Sqrt( Real(1)+dd );
        }
        const Cpx arg = Real(8) - Real(4)*dd + Real(8)*(2-L2)/(L2*sqd);
        const Real a = (sqd + Sqrt(arg)/Real(2)).real();
        const Real b = (a-1)*(a-1)/4;
        const Real c = a+b-1;
        const Real alpha = a-b/c;
        const Real beta = b/c;

        L = L*(a+b*L2)/(1+c*L2);

        if( c > 100 )
        {
            //
            // The standard QR-based algorithm
            //
            MakeHermitian( uplo, A );
            QT = A;
            Scale( Sqrt(c), QT );
            MakeIdentity( QB );
            qr::Explicit( Q, ctrl.colPiv );
            Trrk( uplo, NORMAL, ADJOINT, F(alpha/Sqrt(c)), QT, QB, F(beta), A );
        }
        else
        {
            //
            // Use faster Cholesky-based algorithm since A is well-conditioned
            //
            // TODO: Think of how to better exploit the symmetry of A,
            //       e.g., by halving the work in the first Herk through 
            //       a custom routine for forming L^2, where L is strictly lower
            MakeHermitian( uplo, A );
            Identity( C, n, n );
            Herk( LOWER, ADJOINT, F(c), A, F(1), C );
            Cholesky( LOWER, C );
            ATemp = A;
            Trsm( RIGHT, LOWER, ADJOINT, NON_UNIT, F(1), C, ATemp );
            Trsm( RIGHT, LOWER, NORMAL, NON_UNIT, F(1), C, ATemp );
            Scale( beta, A );
            Axpy( alpha, ATemp, A );
        }

        ++numIts;
        Axpy( F(-1), A, ALast );
        frobNormADiff = HermitianFrobeniusNorm( uplo, ALast );
        if( frobNormADiff <= cubeRootTol && Abs(1-L) <= tol )
            break;
    }
    MakeHermitian( uplo, A );
    return numIts;
}

template<typename F>
inline Int 
QDWH( UpperOrLower uplo, DistMatrix<F>& A, const PolarCtrl& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("herm_polar::QDWH"))
    typedef Base<F> Real;
    MakeHermitian( uplo, A );
    const Real twoEst = TwoNormEstimate( A );
    Scale( F(1)/twoEst, A );

    // The one-norm of the inverse can be replaced with an estimate which is
    // a few times cheaper, e.g., via Higham and Tisseur's block algorithm
    // from "A Block Algorithm for Matrix 1-Norm Estimation, with an Application
    // to 1-Norm Pseudospectra".
    Real sMinUpper;
    DistMatrix<F> Y( A );
    try 
    {   
        Inverse( Y );
        sMinUpper = Real(1) / OneNorm( Y );
    } catch( SingularMatrixException& e ) { sMinUpper = 0; }

    return QDWHInner( uplo, A, sMinUpper, ctrl );
}

template<typename F>
inline Int
QDWH
( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<F>& P, const PolarCtrl& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("herm_polar::QDWH"))
    DistMatrix<F> ACopy( A );
    // NOTE: This might be avoidable
    MakeHermitian( uplo, ACopy );
    const Int numIts = QDWH( uplo, A, ctrl );
    Zeros( P, A.Height(), A.Height() );
    Trrk( uplo, NORMAL, NORMAL, F(1), A, ACopy, F(0), P );
    return numIts;
}

} // namespace herm_polar

} // namespace El

#endif // ifndef EL_POLAR_QDWH_HPP
