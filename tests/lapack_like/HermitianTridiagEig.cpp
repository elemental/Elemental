/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
using namespace El;

template<typename Real,typename=EnableIf<IsReal<Real>>>
void TestGraded
( bool progress,
  bool useQR,
  const herm_tridiag_eig::QRCtrl& qrCtrl,
  bool print )
{
    DEBUG_CSE
    const Int n = 5;
    Output("Testing small graded matrix with ",TypeName<Real>());

    HermitianTridiagEigCtrl<Real> ctrl;
    ctrl.progress = progress;
    ctrl.useQR = useQR;
    ctrl.qrCtrl = qrCtrl;

    Matrix<Real> d(n,1), e(n-1,1);
    d(0) = Real(1);
    d(1) = Pow( Real(10), Real(3) );
    d(2) = Pow( Real(10), Real(6) );
    d(3) = Pow( Real(10), Real(9) ); 
    d(4) = Pow( Real(10), Real(12) );
    e(0) = Real(1);
    e(1) = Pow( Real(10), Real(2) );
    e(2) = Pow( Real(10), Real(5) );
    e(3) = Pow( Real(10), Real(8) );
    auto dOrig(d), eOrig(e);

    const Real TOne = HermitianTridiagOneNorm( d, e );
    Output("|| T ||_1 = ",TOne);
    if( print )
    {
        Print( d, "d" );
        Print( e, "e" );
    }

    Matrix<Real> w, Q;
    Timer timer;
    timer.Start();
    auto info = HermitianTridiagEig( d, e, w, Q, ctrl );
    Output("HermitianTridiagEig: ",timer.Stop()," seconds");
    if( ctrl.useQR )
        Output
        ("Convergence achieved after ",info.qrInfo.numIterations," iterations");
    if( print )
    {
        Print( w, "w" );
        Print( Q, "Q" );
    }

    // R := Q diag(w)
    Matrix<Real> R(Q);
    DiagonalScale( RIGHT, NORMAL, w, R );
    // R := R - T Q
    // TODO: Move into TridiagonalMultiply routine?
    for( Int j=0; j<n; ++j )
    {
        for( Int i=0; i<n; ++i )
        {
            if( i > 0 )
                R(i,j) -= eOrig(i-1)*Q(i-1,j);
            R(i,j) -= dOrig(i)*Q(i,j);
            if( i < n-1 )
                R(i,j) -= eOrig(i)*Q(i+1,j);
        }
    }
    const Real errFrob = FrobeniusNorm( R ); 
    // TODO(poulson): Throw exception if this is too large
    Output("|| T Q - Q diag(w) ||_F / || T ||_1 = ",errFrob/TOne);
    if( print )
        Print( R );
}

template<typename Real,typename=EnableIf<IsReal<Real>>>
void TestRandom
( Int n,
  bool progress,
  bool useQR,
  const herm_tridiag_eig::QRCtrl& qrCtrl,
  bool print )
{
    DEBUG_CSE
    Output("Testing random tridiagonal matrix with ",TypeName<Real>());

    HermitianTridiagEigCtrl<Real> ctrl;
    ctrl.progress = progress;
    ctrl.useQR = useQR;
    ctrl.qrCtrl = qrCtrl;

    Matrix<Real> d, e;
    Uniform( d, n, 1 );
    Uniform( e, n-1, 1 );
    auto dOrig(d), eOrig(e);

    const Real TOne = HermitianTridiagOneNorm( d, e );
    Output("|| T ||_1 = ",TOne);
    if( print )
    {
        Print( d, "d" );
        Print( e, "e" );
    }

    Matrix<Real> w, Q;
    Timer timer;
    timer.Start();
    auto info = HermitianTridiagEig( d, e, w, Q, ctrl );
    Output("HermitianTridiagEig: ",timer.Stop()," seconds");
    if( ctrl.useQR )
        Output
        ("Convergence achieved after ",info.qrInfo.numIterations," iterations");
    if( print )
    {
        Print( w, "w" );
        Print( Q, "Q" );
    }

    // R := Q diag(w)
    Matrix<Real> R(Q);
    DiagonalScale( RIGHT, NORMAL, w, R );
    // R := R - T Q
    // TODO: Move into TridiagonalMultiply routine?
    for( Int j=0; j<n; ++j )
    {
        for( Int i=0; i<n; ++i )
        {
            if( i > 0 )
                R(i,j) -= eOrig(i-1)*Q(i-1,j);
            R(i,j) -= dOrig(i)*Q(i,j);
            if( i < n-1 )
                R(i,j) -= eOrig(i)*Q(i+1,j);
        }
    }
    const Real errFrob = FrobeniusNorm( R ); 
    // TODO(poulson): Throw exception if this is too large
    Output("|| T Q - Q diag(w) ||_F / || T ||_1 = ",errFrob/TOne);
    if( print )
        Print( R );
}

int main( int argc, char* argv[] )
{
    Environment env( argc, argv );

    try
    {
        const Int n = Input("--n","random matrix size",60);
        const bool fullAccuracyTwoByTwo =
          Input
          ("--fullAccuracyTwoByTwo?","full accuracy 2x2 eigenvalues?",true);
        const bool progress = Input("--progress","print progress?",true);
        const bool print = Input("--print","print matrices?",false);
        const bool useQR = Input("--useQR","use QR algorithm?",true);
        ProcessInput();
        PrintInputReport();

        herm_tridiag_eig::QRCtrl qrCtrl;
        qrCtrl.fullAccuracyTwoByTwo = fullAccuracyTwoByTwo;

        TestGraded<float>( progress, useQR, qrCtrl, print );
        TestGraded<double>( progress, useQR, qrCtrl, print );
#ifdef EL_HAVE_QUAD
        TestGraded<Quad>( progress, useQR, qrCtrl, print );
#endif
#ifdef EL_HAVE_QD
        TestGraded<DoubleDouble>( progress, useQR, qrCtrl, print );
        TestGraded<QuadDouble>( progress, useQR, qrCtrl, print );
#endif
#ifdef EL_HAVE_MPC
        TestGraded<BigFloat>( progress, useQR, qrCtrl, print );
#endif

        TestRandom<float>( n, progress, useQR, qrCtrl, print );
        TestRandom<double>( n, progress, useQR, qrCtrl, print );
#ifdef EL_HAVE_QUAD
        TestRandom<Quad>( n, progress, useQR, qrCtrl, print );
#endif
#ifdef EL_HAVE_QD
        TestRandom<DoubleDouble>( n, progress, useQR, qrCtrl, print );
        TestRandom<QuadDouble>( n, progress, useQR, qrCtrl, print );
#endif
#ifdef EL_HAVE_MPC
        TestRandom<BigFloat>( n, progress, useQR, qrCtrl, print );
#endif
    }
    catch( std::exception& e ) { ReportException(e); }

    return 0;
}
