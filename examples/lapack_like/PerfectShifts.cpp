/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.  

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
using namespace El;

// The following example performs a bulge sweep using a pair of eigenvalues
// using both a given datatype and an equivalent substituting the base datatype
// for BigFloat. Note that the deviation is large for perfect shifts, but not
// for most random shifts, even when there are no bulge collapses or vigilant
// deflations.
//

#ifdef EL_HAVE_MPC

template<typename F>
void TestRandomHelper
( const Matrix<F>& H, 
  const HessenbergSchurCtrl& ctrl,
  bool print )
{
    DEBUG_CSE
    typedef Base<F> Real;
    const Int n = H.Height();
    const Real eps = limits::Epsilon<Real>();
    if( n < 2 )
        LogicError("Matrix is too small");

    Matrix<F> T;
    Matrix<Complex<Real>> w;
    Timer timer;

    // Compute just the eigenvalues of the Hessenberg matrix
    Matrix<Complex<Real>> wCopy;
    T = H;
    auto ctrlEig( ctrl );
    ctrlEig.progress = false;
    ctrlEig.fullTriangle = false;
    ctrlEig.wantSchurVecs = false;
    timer.Start();
    auto infoEig = HessenbergSchur( T, wCopy, ctrlEig );
    Output("HessenbergSchur (eigenvalues only): ",timer.Stop()," seconds");
    Output
    ("Convergence achieved after ",infoEig.numIterations," iterations");

    // Use a pair of eigenvalues as shifts
    Output("Testing perfect-shift eigensolver preconditioner");
    T = H;
    auto ctrlSweep( ctrl );
    ctrlSweep.wantSchurVecs = false;
    ctrlSweep.accumulateReflections = false;
    // Find a pair of shifts which are either both real or complex
    // conjugates
    w = wCopy;
    const Int numShifts = 2;
    Int shiftOffset = -1;
    if( IsComplex<F>::value )
    {
        shiftOffset = 0;
    }
    else
    { 
        for( Int i=0; i<n; ++i )
        {
            if( (w(i).imag() == Real(0) && w(i+1).imag() == Real(0)) || 
                w(i) == Conj(w(i+1)) )
            {
                shiftOffset = i;
                break;
            }
        }
    }
    if( shiftOffset == -1 )
        LogicError("Did not find an acceptable shift pair");
    auto wSweep = w( IR(0,numShifts)+shiftOffset, ALL );
    Print( wSweep, "wSweep" );
    timer.Start();
    Matrix<F> Z;
    hess_schur::Sweep( T, wSweep, Z, ctrlSweep );
    Output("hess_schur::Sweep<F>: ",timer.Stop()," seconds");
        
    // Run the same process using BigFloat rather than Base<F>
    typedef ConvertBase<F,BigFloat> FBig;
    Matrix<Complex<BigFloat>> wBig, wBigSweep;
    Copy( wSweep, wBigSweep );
    Matrix<FBig> TBig, ZBig;
    Copy( H, TBig ); 
    timer.Start();
    hess_schur::Sweep( TBig, wBigSweep, ZBig, ctrlSweep );
    Output("hess_schur::Sweep<FBig>: ",timer.Stop()," seconds");

    // Compare the results using different precisions
    Matrix<FBig> TDiff;
    Copy( T, TDiff );
    Axpy( FBig(-1), TBig, TDiff ); 
    const BigFloat TBigFrob = FrobeniusNorm( TBig );
    const BigFloat TDiffFrob = FrobeniusNorm( TDiff );
    Output("|| T - TBig ||_F / || TBig ||_F = ",TDiffFrob/TBigFrob);
    if( print )
    {
        Print( T, "T after sweep" );
        Print( TBig, "TBig after sweep" );
        Print( TDiff, "TDiff after sweep" );
    }
}

template<typename F>
void TestRandom( Int n, const HessenbergSchurCtrl& ctrl, bool print )
{
    DEBUG_CSE
    Output("Testing uniform Hessenberg with ",TypeName<F>());

    Matrix<F> H;
    Uniform( H, n, n );
    MakeTrapezoidal( UPPER, H, -1 );
    if( print )
        Print( H, "H" );

    TestRandomHelper( H, ctrl, print );
    Output("");
}

int main( int argc, char* argv[] )
{
    Environment env( argc, argv );

    try
    {
        const Int n = Input("--n","random matrix size",60);
        const Int algInt = Input("--alg","AED: 0, MultiBulge: 1, Simple: 2",0);
        const Int minMultiBulgeSize =
          Input
          ("--minMultiBulgeSize",
           "minimum size for using a multi-bulge algorithm",75);
        const bool accumulate =
          Input("--accumulate","accumulate reflections?",true);
        const bool sortShifts =
          Input("--sortShifts","sort shifts for AED?",true);
        const bool progress = Input("--progress","print progress?",true);
        const bool print = Input("--print","print matrices?",false);
        const mpfr_prec_t prec = Input("--prec","MPFR precision",256);
        ProcessInput();
        PrintInputReport();

        mpfr::SetPrecision(prec);

        HessenbergSchurCtrl ctrl;
        ctrl.alg = static_cast<HessenbergSchurAlg>(algInt);
        ctrl.minMultiBulgeSize = minMultiBulgeSize;
        ctrl.accumulateReflections = accumulate;
        ctrl.sortShifts = sortShifts;
        ctrl.progress = progress;

        TestRandom<float>( n, ctrl, print );
        TestRandom<double>( n, ctrl, print );
#ifdef EL_HAVE_QD
        TestRandom<DoubleDouble>( n, ctrl, print );
#endif
    }
    catch( std::exception& e ) { ReportException(e); }

    return 0;
}

#else // ifdef EL_HAVE_MPC

// We provide a trivial main program to prevent linker errors on some systems
int main( int argc, char* argv[] )
{
    Environment env( argc, argv );
    return 0;
}

#endif // ifdef EL_HAVE_MPC
