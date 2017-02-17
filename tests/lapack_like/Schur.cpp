/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.  

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
using namespace El;

template<typename F>
void TestRandomHelper
( const Matrix<F>& A, 
  const HessenbergSchurCtrl& hessSchurCtrl,
  bool print )
{
    EL_DEBUG_CSE
    typedef Base<F> Real;
    const Int n = A.Height();
    const Real eps = limits::Epsilon<Real>();

    SchurCtrl<Base<F>> ctrl;
    ctrl.hessSchurCtrl = hessSchurCtrl;

    Matrix<F> T, Z;
    Matrix<Complex<Real>> w;
    Timer timer;

    T = A;
    timer.Start();
    Schur( T, w, Z, ctrl );
    Output("Schur: ",timer.Stop()," seconds");
    if( print )
    {
        Print( w, "w" );
        Print( Z, "Z" );
        Print( T, "T" );
    }

    Matrix<F> R;
    Gemm( NORMAL, NORMAL, F(1), Z, T, R );
    Gemm( NORMAL, NORMAL, F(1), A, Z, F(-1), R );
    const Real errFrob = FrobeniusNorm( R ); 
    const Real AFrob = FrobeniusNorm( A );
    const Real relErr = errFrob / (eps*n*AFrob);
    Output("|| A ||_F = ",AFrob);
    Output("|| A Z - Z T ||_F / (eps n || A ||_F) = ",relErr);
    if( print )
        Print( R );
    // TODO(poulson): A more refined failure condition
    if( relErr > Real(100) )
        LogicError("Relative error was unacceptably large");
    Output("Passed test");
    Output("");
}

template<typename F>
void TestRandomHelper
( const AbstractDistMatrix<F>& A, 
  const HessenbergSchurCtrl& hessSchurCtrl,
  bool print )
{
    EL_DEBUG_CSE
    typedef Base<F> Real;
    const Int n = A.Height();
    const Real eps = limits::Epsilon<Real>();
    const Grid& grid = A.Grid();

    SchurCtrl<Base<F>> ctrl;
    ctrl.hessSchurCtrl = hessSchurCtrl;

    DistMatrix<F,MC,MR,BLOCK> T(grid), Z(grid);
    DistMatrix<Complex<Real>,STAR,STAR> w(grid);
    Timer timer;

    T = A;
    timer.Start();
    Schur( T, w, Z, ctrl );
    if( grid.Rank() == 0 )
        Output("Schur: ",timer.Stop()," seconds");
    if( print )
    {
        Print( w, "w" );
        Print( Z, "Z" );
        Print( T, "T" );
    }

    DistMatrix<F> R(grid);
    Gemm( NORMAL, NORMAL, F(1), Z, T, R );
    Gemm( NORMAL, NORMAL, F(1), A, Z, F(-1), R );
    const Real errFrob = FrobeniusNorm( R ); 
    const Real AFrob = FrobeniusNorm( A );
    const Real relErr = errFrob / (eps*n*AFrob);
    if( grid.Rank() == 0 )
    {
        Output("|| A ||_F = ",AFrob);
        Output("|| A Z - Z T ||_F / (eps n || A ||_F) = ",relErr);
    }
    if( print )
        Print( R );
    // TODO(poulson): A more refined failure condition
    if( relErr > Real(100) )
        LogicError("Relative error was unacceptably large");
    if( grid.Rank() == 0 )
    {
        Output("Passed test");
        Output("");
    }
}

template<typename F>
void TestRandom( Int n, const HessenbergSchurCtrl& ctrl, bool print )
{
    EL_DEBUG_CSE
    Output("Testing uniform with ",TypeName<F>());
    Matrix<F> A;
    Uniform( A, n, n );
    if( print )
        Print( A, "A" );
    TestRandomHelper( A, ctrl, print );
}

template<typename F>
void TestRandom
( Int n, const Grid& grid, const HessenbergSchurCtrl& ctrl, bool print )
{
    EL_DEBUG_CSE
    Output("Testing uniform with ",TypeName<F>());
    DistMatrix<F> A(grid);
    Uniform( A, n, n );
    if( print )
        Print( A, "A" );
    TestRandomHelper( A, ctrl, print );
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
        const bool sequential = Input("--sequential","test sequential?",true);
        const bool distributed =
          Input("--distributed","test distributed?",true);
        const bool progress = Input("--progress","print progress?",true);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        HessenbergSchurCtrl ctrl;
        ctrl.alg = static_cast<HessenbergSchurAlg>(algInt);
        ctrl.minMultiBulgeSize = minMultiBulgeSize;
        ctrl.accumulateReflections = accumulate;
        ctrl.sortShifts = sortShifts;
        ctrl.progress = progress;

        // TODO(poulson): Allow the grid dimensions to be selected
        const Grid grid( mpi::COMM_WORLD );

        if( sequential && grid.Rank() == 0 )
        {
            TestRandom<float>( n, ctrl, print );
            TestRandom<Complex<float>>( n, ctrl, print );
            TestRandom<double>( n, ctrl, print );
            TestRandom<Complex<double>>( n, ctrl, print );
#ifdef EL_HAVE_QUAD
            TestRandom<Quad>( n, ctrl, print );
            TestRandom<Complex<Quad>>( n, ctrl, print );
#endif
#ifdef EL_HAVE_QD
            TestRandom<DoubleDouble>( n, ctrl, print );
            TestRandom<Complex<DoubleDouble>>( n, ctrl, print );
            TestRandom<QuadDouble>( n, ctrl, print );
            TestRandom<Complex<QuadDouble>>( n, ctrl, print );
#endif
#ifdef EL_HAVE_MPC
            TestRandom<BigFloat>( n, ctrl, print );
            TestRandom<Complex<BigFloat>>( n, ctrl, print );
#endif
        }
        if( distributed )
        {
            TestRandom<float>( n, grid, ctrl, print );
            TestRandom<Complex<float>>( n, grid, ctrl, print );
            TestRandom<double>( n, grid, ctrl, print );
            TestRandom<Complex<double>>( n, grid, ctrl, print );
#ifdef EL_HAVE_QUAD
            TestRandom<Quad>( n, grid, ctrl, print );
            TestRandom<Complex<Quad>>( n, grid, ctrl, print );
#endif
#ifdef EL_HAVE_QD
            TestRandom<DoubleDouble>( n, grid, ctrl, print );
            TestRandom<Complex<DoubleDouble>>( n, grid, ctrl, print );
            TestRandom<QuadDouble>( n, grid, ctrl, print );
            TestRandom<Complex<QuadDouble>>( n, grid, ctrl, print );
#endif
#ifdef EL_HAVE_MPC
            TestRandom<BigFloat>( n, grid, ctrl, print );
            TestRandom<Complex<BigFloat>>( n, grid, ctrl, print );
#endif
        }
    }
    catch( std::exception& e ) { ReportException(e); }

    return 0;
}
