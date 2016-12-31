/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

// The following example performs a bulge sweep using a pair of eigenvalues
// using both a given datatype and an equivalent substituting the base datatype
// for BigFloat. Note that the deviation is large for perfect shifts, but not
// for most random shifts, even when there are no bulge collapses or vigilant
// deflations.
//

#ifdef EL_HAVE_MPC

template<typename Field>
void TestRandomHelper
( const El::Matrix<Field>& H,
  const El::HessenbergSchurCtrl& ctrl,
  bool print )
{
    EL_DEBUG_CSE
    typedef El::Base<Field> Real;
    const El::Int n = H.Height();
    const Real eps = El::limits::Epsilon<Real>();
    if( n < 2 )
        El::LogicError("Matrix is too small");

    El::Matrix<Field> T;
    El::Matrix<El::Complex<Real>> w;
    El::Timer timer;

    // Compute just the eigenvalues of the Hessenberg matrix
    El::Matrix<El::Complex<Real>> wCopy;
    T = H;
    auto ctrlEig( ctrl );
    ctrlEig.progress = false;
    ctrlEig.fullTriangle = false;
    ctrlEig.wantSchurVecs = false;
    timer.Start();
    auto infoEig = El::HessenbergSchur( T, wCopy, ctrlEig );
    El::Output("HessenbergSchur (eigenvalues only): ",timer.Stop()," seconds");
    El::Output
    ("Convergence achieved after ",infoEig.numIterations," iterations");

    // Use a pair of eigenvalues as shifts
    El::Output("Testing perfect-shift eigensolver preconditioner");
    T = H;
    auto ctrlSweep( ctrl );
    ctrlSweep.wantSchurVecs = false;
    ctrlSweep.accumulateReflections = false;
    // Find a pair of shifts which are either both real or complex
    // conjugates
    w = wCopy;
    const El::Int numShifts = 2;
    El::Int shiftOffset = -1;
    if( El::IsComplex<Field>::value )
    {
        shiftOffset = 0;
    }
    else
    {
        for( El::Int i=0; i<n; ++i )
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
        El::LogicError("Did not find an acceptable shift pair");
    auto wSweep = w( El::IR(0,numShifts)+shiftOffset, El::ALL );
    El::Print( wSweep, "wSweep" );
    timer.Start();
    El::Matrix<Field> Z;
    El::hess_schur::Sweep( T, wSweep, Z, ctrlSweep );
    El::Output("hess_schur::Sweep<Field>: ",timer.Stop()," seconds");

    // Run the same process using BigFloat rather than El::Base<Field>
    typedef El::ConvertBase<Field,El::BigFloat> FieldBig;
    El::Matrix<El::Complex<El::BigFloat>> wBig, wBigSweep;
    El::Copy( wSweep, wBigSweep );
    El::Matrix<FieldBig> TBig, ZBig;
    El::Copy( H, TBig );
    timer.Start();
    El::hess_schur::Sweep( TBig, wBigSweep, ZBig, ctrlSweep );
    El::Output("hess_schur::Sweep<FBig>: ",timer.Stop()," seconds");

    // Compare the results using different precisions
    El::Matrix<FieldBig> TDiff;
    El::Copy( T, TDiff );
    El::Axpy( FieldBig(-1), TBig, TDiff );
    const El::BigFloat TBigFrob = El::FrobeniusNorm( TBig );
    const El::BigFloat TDiffFrob = El::FrobeniusNorm( TDiff );
    El::Output("|| T - TBig ||_F / || TBig ||_F = ",TDiffFrob/TBigFrob);
    if( print )
    {
        El::Print( T, "T after sweep" );
        El::Print( TBig, "TBig after sweep" );
        El::Print( TDiff, "TDiff after sweep" );
    }
}

template<typename Field>
void TestRandom( El::Int n, const El::HessenbergSchurCtrl& ctrl, bool print )
{
    EL_DEBUG_CSE
    El::Output("Testing uniform Hessenberg with ",El::TypeName<Field>());

    El::Matrix<Field> H;
    El::Uniform( H, n, n );
    El::MakeTrapezoidal( El::UPPER, H, -1 );
    if( print )
        El::Print( H, "H" );

    TestRandomHelper( H, ctrl, print );
    El::Output("");
}

int main( int argc, char* argv[] )
{
    El::Environment env( argc, argv );

    try
    {
        const El::Int n = El::Input("--n","random matrix size",60);
        const El::Int algInt =
          El::Input("--alg","AED: 0, MultiBulge: 1, Simple: 2",0);
        const El::Int minMultiBulgeSize =
          El::Input
          ("--minMultiBulgeSize",
           "minimum size for using a multi-bulge algorithm",75);
        const bool accumulate =
          El::Input("--accumulate","accumulate reflections?",true);
        const bool sortShifts =
          El::Input("--sortShifts","sort shifts for AED?",true);
        const bool progress = El::Input("--progress","print progress?",true);
        const bool print = El::Input("--print","print matrices?",false);
        const mpfr_prec_t prec = El::Input("--prec","MPFR precision",256);
        El::ProcessInput();
        El::PrintInputReport();

        El::mpfr::SetPrecision(prec);

        El::HessenbergSchurCtrl ctrl;
        ctrl.alg = static_cast<El::HessenbergSchurAlg>(algInt);
        ctrl.minMultiBulgeSize = minMultiBulgeSize;
        ctrl.accumulateReflections = accumulate;
        ctrl.sortShifts = sortShifts;
        ctrl.progress = progress;

        TestRandom<float>( n, ctrl, print );
        TestRandom<double>( n, ctrl, print );
#ifdef EL_HAVE_QD
        TestRandom<El::DoubleDouble>( n, ctrl, print );
        TestRandom<El::QuadDouble>( n, ctrl, print );
#endif
    }
    catch( std::exception& e ) { El::ReportException(e); }

    return 0;
}

#else // ifdef EL_HAVE_MPC

// We provide a trivial main program to prevent linker errors on some systems
int main( int argc, char* argv[] )
{
    El::Environment env( argc, argv );
    return 0;
}

#endif // ifdef EL_HAVE_MPC
