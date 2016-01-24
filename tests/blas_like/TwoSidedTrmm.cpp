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
  UnitOrNonUnit diag,
  const DistMatrix<F>& A,
  const DistMatrix<F>& B,
  const DistMatrix<F>& AOrig,
  bool print )
{
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    const Int m = AOrig.Height();

    const Int k=100;
    DistMatrix<F> X(g), Y(g), Z(g);
    Uniform( X, m, k );
    Y = X;
    Zeros( Z, m, k );

    if( uplo == LOWER )
    {
        // Test correctness by comparing the application of A against a 
        // random set of k vectors to the application of 
        // tril(B)^H AOrig tril(B)
        Trmm( LEFT, LOWER, NORMAL, diag, F(1), B, Y );
        Hemm( LEFT, LOWER, F(1), AOrig, Y, F(0), Z );
        Trmm( LEFT, LOWER, ADJOINT, diag, F(1), B, Z );
        Hemm( LEFT, LOWER, F(-1), A, X, F(1), Z );
        Real infNormAOrig = HermitianInfinityNorm( uplo, AOrig );
        Real frobNormAOrig = HermitianFrobeniusNorm( uplo, AOrig );
        Real infNormA = HermitianInfinityNorm( uplo, A );
        Real frobNormA = HermitianFrobeniusNorm( uplo, A );
        Real oneNormError = OneNorm( Z );
        Real infNormError = InfinityNorm( Z );
        Real frobNormError = FrobeniusNorm( Z );
        if( g.Rank() == 0 )
        {
            Output("||AOrig||_1 = ||AOrig||_oo     = ",infNormAOrig);
            Output("||AOrig||_F                    = ",frobNormAOrig);
            Output("||A||_1 = ||A||_oo             = ",infNormA);
            Output("||A||_F                        = ",frobNormA);
            Output("||A X - L^H AOrig L X||_1  = ",oneNormError);
            Output("||A X - L^H AOrig L X||_oo = ",infNormError);
            Output("||A X - L^H AOrig L X||_F  = ",frobNormError);
        }
    }
    else
    {
        // Test correctness by comparing the application of A against a 
        // random set of k vectors to the application of 
        // triu(B) AOrig triu(B)^H
        Trmm( LEFT, UPPER, ADJOINT, diag, F(1), B, Y );
        Hemm( LEFT, UPPER, F(1), AOrig, Y, F(0), Z );
        Trmm( LEFT, UPPER, NORMAL, diag, F(1), B, Z );
        Hemm( LEFT, UPPER, F(-1), A, X, F(1), Z );
        Real infNormAOrig = HermitianInfinityNorm( uplo, AOrig );
        Real frobNormAOrig = HermitianFrobeniusNorm( uplo, AOrig );
        Real infNormA = HermitianInfinityNorm( uplo, A );
        Real frobNormA = HermitianFrobeniusNorm( uplo, A );
        Real oneNormError = OneNorm( Z );
        Real infNormError = InfinityNorm( Z );
        Real frobNormError = FrobeniusNorm( Z );
        if( g.Rank() == 0 )
        {
            Output("||AOrig||_1 = ||AOrig||_oo     = ",infNormAOrig);
            Output("||AOrig||_F                    = ",frobNormAOrig);
            Output("||A||_1 = ||A||_oo             = ",infNormA);
            Output("||A||_F                        = ",frobNormA);
            Output("||A X - U AOrig U^H X||_1  = ",oneNormError);
            Output("||A X - U AOrig U^H X||_oo = ",infNormError);
            Output("||A X - U AOrig U^H X||_F  = ",frobNormError);
        }
    }
}

template<typename F,typename=EnableIf<IsBlasScalar<F>>> 
void TestTwoSidedTrmm
( UpperOrLower uplo,
  UnitOrNonUnit diag, 
  Int m,
  const Grid& g,
  bool scalapack,
  bool print,
  bool correctness )
{
    if( g.Rank() == 0 )
        Output("Testing with ",TypeName<F>());
    DistMatrix<F> A(g), B(g), AOrig(g);
    HermitianUniformSpectrum( A, m, 1, 10 );
    HermitianUniformSpectrum( B, m, 1, 10 );
    MakeTrapezoidal( uplo, B );
    if( correctness )
        AOrig = A;
    if( print )
    {
        Print( A, "A" );
        Print( B, "B" );
    }

    if( scalapack )
    {
        DistMatrix<F,MC,MR,BLOCK> ABlock( A ), BBlock( B );
        if( g.Rank() == 0 )
            Output("  Starting ScaLAPACK TwoSidedTrmm");
        const double startTime = mpi::Time();
        TwoSidedTrmm( uplo, diag, ABlock, BBlock );
        const double runTime = mpi::Time() - startTime;
        double gFlops = Pow(double(m),3.)/(runTime*1.e9);
        if( IsComplex<F>::value )
            gFlops *= 4.;
        if( g.Rank() == 0 )
            Output("  Time = ",runTime," seconds. GFlops = ",gFlops);
    }

    if( g.Rank() == 0 )
        Output("  Starting TwoSidedTrmm");
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    TwoSidedTrmm( uplo, diag, A, B );
    mpi::Barrier( g.Comm() );
    const double runTime = mpi::Time() - startTime;
    double gFlops = Pow(double(m),3.)/(runTime*1.e9);
    if( IsComplex<F>::value )
        gFlops *= 4.;
    if( g.Rank() == 0 )
        Output("  Time = ",runTime," seconds. GFlops = ",gFlops);
    if( print )
        Print( A, "A after reduction" );
    if( correctness )
        TestCorrectness( uplo, diag, A, B, AOrig, print );
}

template<typename F> 
void TestTwoSidedTrmm
( UpperOrLower uplo,
  UnitOrNonUnit diag, 
  Int m,
  const Grid& g,
  bool print,
  bool correctness )
{
    if( g.Rank() == 0 )
        Output("Testing with ",TypeName<F>());
    DistMatrix<F> A(g), B(g), AOrig(g);
    HermitianUniformSpectrum( A, m, 1, 10 );
    HermitianUniformSpectrum( B, m, 1, 10 );
    MakeTrapezoidal( uplo, B );
    if( correctness )
        AOrig = A;
    if( print )
    {
        Print( A, "A" );
        Print( B, "B" );
    }

    if( g.Rank() == 0 )
        Output("  Starting TwoSidedTrmm");
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    TwoSidedTrmm( uplo, diag, A, B );
    mpi::Barrier( g.Comm() );
    const double runTime = mpi::Time() - startTime;
    double gFlops = Pow(double(m),3.)/(runTime*1.e9);
    if( IsComplex<F>::value )
        gFlops *= 4.;
    if( g.Rank() == 0 )
        Output("  Time = ",runTime," seconds. GFlops = ",gFlops);
    if( print )
        Print( A, "A after reduction" );
    if( correctness )
        TestCorrectness( uplo, diag, A, B, AOrig, print );
}

int 
main( int argc, char* argv[] )
{
    Environment env( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const Int commRank = mpi::Rank( comm );
    const Int commSize = mpi::Size( comm );

    try
    {
        Int r = Input("--r","height of process grid",0);
        const bool colMajor = Input("--colMajor","column-major ordering?",true);
        const char uploChar = Input
            ("--uplo","lower or upper triangular storage: L/U",'L');
        const char diagChar = Input("--unit","(non-)unit diagonal: N/U",'N');
        const Int m = Input("--m","height of matrix",100);
        const Int blocksize = Input("--blocksize","algorithmic blocksize",96);
#ifdef EL_HAVE_SCALAPACK
        const bool scalapack = Input("--scalapack","test ScaLAPACK?",true);
        const Int mb = Input("--mb","block height",32);
        const Int nb = Input("--nb","block width",32);
#else
        const bool scalapack = false;
        const Int mb = 32;
        const Int nb = 32;
#endif
        const bool correctness = Input
            ("--correctness","test correctness?",true);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        if( r == 0 )
            r = Grid::FindFactor( commSize );
        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( comm, r, order );
        const UpperOrLower uplo = CharToUpperOrLower( uploChar );
        const UnitOrNonUnit diag = CharToUnitOrNonUnit( diagChar );
        SetBlocksize( blocksize );
        SetDefaultBlockHeight( mb );
        SetDefaultBlockWidth( nb );

        ComplainIfDebug();
        if( commRank == 0 )
            Output("Will test TwoSidedTrmm",uploChar,diagChar);

        TestTwoSidedTrmm<float>
        ( uplo, diag, m, g, scalapack, print, correctness );
        TestTwoSidedTrmm<Complex<float>>
        ( uplo, diag, m, g, scalapack, print, correctness );

        TestTwoSidedTrmm<double>
        ( uplo, diag, m, g, scalapack, print, correctness );
        TestTwoSidedTrmm<Complex<double>>
        ( uplo, diag, m, g, scalapack, print, correctness );

#ifdef EL_HAVE_QD
        TestTwoSidedTrmm<DoubleDouble>
        ( uplo, diag, m, g, print, correctness );
        TestTwoSidedTrmm<QuadDouble>
        ( uplo, diag, m, g, print, correctness );
#endif

#ifdef EL_HAVE_QUAD
        TestTwoSidedTrmm<Quad>
        ( uplo, diag, m, g, print, correctness );
        TestTwoSidedTrmm<Complex<Quad>>
        ( uplo, diag, m, g, print, correctness );
#endif

#ifdef EL_HAVE_MPC
        TestTwoSidedTrmm<BigFloat>
        ( uplo, diag, m, g, print, correctness );
#endif
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
