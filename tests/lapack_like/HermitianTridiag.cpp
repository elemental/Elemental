/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

template<typename F> 
void TestCorrectness
( UpperOrLower uplo, 
  const DistMatrix<F>& A, 
  const DistMatrix<F,STAR,STAR>& t,
        DistMatrix<F>& AOrig,
  bool print,
  bool display )
{
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    const Int m = AOrig.Height();
    const Real infNormAOrig = HermitianInfinityNorm( uplo, AOrig );
    const Real frobNormAOrig = HermitianFrobeniusNorm( uplo, AOrig );
    if( g.Rank() == 0 )
        Output("Testing error...");

    // Grab the diagonal and subdiagonal of the symmetric tridiagonal matrix
    Int subdiagonal = ( uplo==LOWER ? -1 : +1 );
    auto d = GetRealPartOfDiagonal(A);
    auto e = GetRealPartOfDiagonal(A,subdiagonal);
     
    // Grab a full copy of e so that we may fill the opposite subdiagonal 
    DistMatrix<Real,STAR,STAR> e_STAR_STAR( e );
    DistMatrix<Real,MD,STAR> eOpposite(g);
    eOpposite.SetRoot( A.DiagonalRoot(-subdiagonal) );
    eOpposite.AlignCols( A.DiagonalAlign(-subdiagonal) );
    eOpposite = e_STAR_STAR;
    
    // Zero B and then fill its tridiagonal
    DistMatrix<F> B(g);
    B.AlignWith( A );
    Zeros( B, m, m );
    SetRealPartOfDiagonal( B, d );
    SetRealPartOfDiagonal( B, e, subdiagonal );
    SetRealPartOfDiagonal( B, eOpposite, -subdiagonal );
    if( print )
        Print( B, "Tridiagonal" );
    if( display )
        Display( B, "Tridiagonal" );

    // Reverse the accumulated Householder transforms, ignoring symmetry
    herm_tridiag::ApplyQ( LEFT, uplo, NORMAL, A, t, B );
    herm_tridiag::ApplyQ( RIGHT, uplo, ADJOINT, A, t, B );
    if( print )
        Print( B, "Rotated tridiagonal" );
    if( display )
        Display( B, "Rotated tridiagonal" );

    // Compare the appropriate triangle of AOrig and B
    MakeTrapezoidal( uplo, AOrig );
    MakeTrapezoidal( uplo, B );
    B -= AOrig;
    if( print )
        Print( B, "Error in rotated tridiagonal" );
    if( display )
        Display( B, "Error in rotated tridiagonal" );
    const Real infNormError = HermitianInfinityNorm( uplo, B );
    const Real frobNormError = HermitianFrobeniusNorm( uplo, B );

    // Compute || I - Q Q^H ||
    MakeIdentity( B );
    herm_tridiag::ApplyQ( RIGHT, uplo, ADJOINT, A, t, B );
    DistMatrix<F> QHAdj( g );
    Adjoint( B, QHAdj );
    MakeIdentity( B );
    herm_tridiag::ApplyQ( LEFT, uplo, NORMAL, A, t, B );
    QHAdj -= B;
    herm_tridiag::ApplyQ( RIGHT, uplo, ADJOINT, A, t, B );
    ShiftDiagonal( B, F(-1) );
    const Real infNormQError = InfinityNorm( B );
    const Real frobNormQError = FrobeniusNorm( B ); 

    if( g.Rank() == 0 )
        Output
        ("    ||A||_oo = ",infNormAOrig,"\n",
         "    ||A||_F  = ",frobNormAOrig,"\n",
         "    || I - Q^H Q ||_oo = ",infNormQError,"\n",
         "    || I - Q^H Q ||_F  = ",frobNormQError,"\n",
         "    ||A - Q T Q^H||_oo = ",infNormError,"\n",
         "    ||A - Q T Q^H||_F  = ",frobNormError);
}

template<typename F>
void InnerTestHermitianTridiag
( UpperOrLower uplo,
        DistMatrix<F>& A,
        DistMatrix<F,STAR,STAR>& t,
  const HermitianTridiagCtrl<F>& ctrl,
  bool testCorrectness,
  bool print,
  bool display )
{
    DistMatrix<F> AOrig( A ), ACopy( A );
    const Int m = A.Height();
    const Grid& g = A.Grid();

    if( g.Rank() == 0 )
        Output("  Starting tridiagonalization...");
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    HermitianTridiag( uplo, A, t, ctrl );
    mpi::Barrier( g.Comm() );
    const double runTime = mpi::Time() - startTime;
    const double realGFlops = 16./3.*Pow(double(m),3.)/(1.e9*runTime);
    const double gFlops = ( IsComplex<F>::value ? 4*realGFlops : realGFlops );
    if( g.Rank() == 0 )
        Output("  ",runTime," seconds (",gFlops," GFlop/s)");
    if( print )
    {
        Print( A, "A after HermitianTridiag" );
        Print( t, "t after HermitianTridiag" );
    }
    if( display )
    {
        Display( A, "A after HermitianTridiag" );
        Display( t, "t after HermitianTridiag" );
    }
    if( testCorrectness )
        TestCorrectness( uplo, A, t, AOrig, print, display );
    A = ACopy;
}

template<typename F>
void TestHermitianTridiag
( const Grid& g,
  UpperOrLower uplo,
  Int m,
  Int nbLocal,
  bool avoidTrmv,
  bool testCorrectness,
  bool print,
  bool display )
{
    DistMatrix<F> A(g), AOrig(g);
    DistMatrix<F,STAR,STAR> t(g);
    if( g.Rank() == 0 )
        Output("Testing with ",TypeName<F>());

    HermitianTridiagCtrl<F> ctrl;
    ctrl.symvCtrl.bsize = nbLocal;
    ctrl.symvCtrl.avoidTrmvBasedLocalSymv = avoidTrmv;

    Wigner( A, m );
    if( testCorrectness )
        AOrig = A;
    if( print )
        Print( A, "A" );
    if( display )
        Display( A, "A" );

    if( g.Rank() == 0 )
        Output("Normal algorithm:");
    ctrl.approach = HERMITIAN_TRIDIAG_NORMAL;
    InnerTestHermitianTridiag
    ( uplo, A, t, ctrl, testCorrectness, print, display );

    if( g.Rank() == 0 )
        Output("Square row-major algorithm:");
    ctrl.approach = HERMITIAN_TRIDIAG_SQUARE;
    ctrl.order = ROW_MAJOR;
    InnerTestHermitianTridiag
    ( uplo, A, t, ctrl, testCorrectness, print, display );

    if( g.Rank() == 0 )
        Output("Square column-major algorithm:");
    ctrl.approach = HERMITIAN_TRIDIAG_SQUARE;
    ctrl.order = COLUMN_MAJOR;
    InnerTestHermitianTridiag
    ( uplo, A, t, ctrl, testCorrectness, print, display );
}

int 
main( int argc, char* argv[] )
{
    Environment env( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::Rank( comm );
    const int commSize = mpi::Size( comm );

    try
    {
        Int r = Input("--gridHeight","height of process grid",0);
        const bool colMajor = Input("--colMajor","column-major ordering?",true);
        const char uploChar = Input("--uplo","upper or lower storage: L/U",'L');
        const Int m = Input("--height","height of matrix",100);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const Int nbLocal = Input("--nbLocal","local blocksize",32);
        const bool avoidTrmv = 
            Input("--avoidTrmv","avoid Trmv local Symv",true);
        const bool testCorrectness = Input
            ("--correctness","test correctness?",true);
        const bool print = Input("--print","print matrices?",false);
        const bool display = Input("--display","display matrices?",false);
        const bool testReal = Input("--testReal","test real matrices?",true);
        const bool testCpx = Input("--testCpx","test complex matrices?",true);
#ifdef EL_HAVE_MPC
        const mpfr_prec_t prec = Input("--prec","MPFR precision",256);
#endif
        ProcessInput();
        PrintInputReport();

#ifdef EL_HAVE_MPC
        mpc::SetPrecision( prec );
#endif

        if( r == 0 )
            r = Grid::FindFactor( commSize );
        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( comm, r, order );
        const UpperOrLower uplo = CharToUpperOrLower( uploChar );
        SetBlocksize( nb );

        ComplainIfDebug();
        if( commRank == 0 )
            Output("Will test HermitianTridiag",uploChar);

        if( testReal )
            TestHermitianTridiag<float>
            ( g, uplo, m, nbLocal, avoidTrmv, testCorrectness, print, display );
        if( testCpx )
            TestHermitianTridiag<Complex<float>>
            ( g, uplo, m, nbLocal, avoidTrmv, testCorrectness, print, display );

        if( testReal )
            TestHermitianTridiag<double>
            ( g, uplo, m, nbLocal, avoidTrmv, testCorrectness, print, display );
        if( testCpx )
            TestHermitianTridiag<Complex<double>>
            ( g, uplo, m, nbLocal, avoidTrmv, testCorrectness, print, display );

#ifdef EL_HAVE_QD
        if( testReal )
        {
            TestHermitianTridiag<DoubleDouble>
            ( g, uplo, m, nbLocal, avoidTrmv, testCorrectness, print, display );
            TestHermitianTridiag<QuadDouble>
            ( g, uplo, m, nbLocal, avoidTrmv, testCorrectness, print, display );
        }
#endif

#ifdef EL_HAVE_QUAD
        if( testReal )
            TestHermitianTridiag<Quad>
            ( g, uplo, m, nbLocal, avoidTrmv, testCorrectness, print, display );
        if( testCpx )
            TestHermitianTridiag<Complex<Quad>>
            ( g, uplo, m, nbLocal, avoidTrmv, testCorrectness, print, display );
#endif

#ifdef EL_HAVE_MPC
        if( testReal )
            TestHermitianTridiag<BigFloat>
            ( g, uplo, m, nbLocal, avoidTrmv, testCorrectness, print, display );
#endif
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
