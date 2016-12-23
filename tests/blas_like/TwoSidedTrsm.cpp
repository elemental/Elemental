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
        // tril(B)^-1 AOrig tril(B)^-H
        Trsm( LEFT, LOWER, ADJOINT, diag, F(1), B, Y );
        Hemm( LEFT, LOWER, F(1), AOrig, Y, F(0), Z );
        Trsm( LEFT, LOWER, NORMAL, diag, F(1), B, Z );
        Hemm( LEFT, LOWER, F(-1), A, X, F(1), Z );
        Real infNormAOrig = HermitianInfinityNorm( uplo, AOrig );
        Real frobNormAOrig = HermitianFrobeniusNorm( uplo, AOrig );
        Real infNormA = HermitianInfinityNorm( uplo, A );
        Real frobNormA = HermitianFrobeniusNorm( uplo, A );
        Real oneNormError = OneNorm( Z );
        Real infNormError = InfinityNorm( Z );
        Real frobNormError = FrobeniusNorm( Z );
        OutputFromRoot
        (g.Comm(),
         "||AOrig||_1 = ||AOrig||_oo     = ",infNormAOrig,"\n",Indent(),
         "||AOrig||_F                    = ",frobNormAOrig,"\n",Indent(),
         "||A||_1 = ||A||_oo             = ",infNormA,"\n",Indent(),
         "||A||_F                        = ",frobNormA,"\n",Indent(),
         "||A X - L^-1 AOrig L^-H X||_1  = ",oneNormError,"\n",Indent(),
         "||A X - L^-1 AOrig L^-H X||_oo = ",infNormError,"\n",Indent(),
         "||A X - L^-1 AOrig L^-H X||_F  = ",frobNormError);
    }
    else
    {
        // Test correctness by comparing the application of A against a
        // random set of k vectors to the application of
        // triu(B)^-H AOrig triu(B)^-1
        Trsm( LEFT, UPPER, NORMAL, diag, F(1), B, Y );
        Hemm( LEFT, UPPER, F(1), AOrig, Y, F(0), Z );
        Trsm( LEFT, UPPER, ADJOINT, diag, F(1), B, Z );
        Hemm( LEFT, UPPER, F(-1), A, X, F(1), Z );
        Real infNormAOrig = HermitianInfinityNorm( uplo, AOrig );
        Real frobNormAOrig = HermitianFrobeniusNorm( uplo, AOrig );
        Real infNormA = HermitianInfinityNorm( uplo, A );
        Real frobNormA = HermitianFrobeniusNorm( uplo, A );
        Real oneNormError = OneNorm( Z );
        Real infNormError = InfinityNorm( Z );
        Real frobNormError = FrobeniusNorm( Z );
        OutputFromRoot
        (g.Comm(),
         "||AOrig||_1 = ||AOrig||_oo     = ",infNormAOrig,"\n",Indent(),
         "||AOrig||_F                    = ",frobNormAOrig,"\n",Indent(),
         "||A||_1 = ||A||_oo             = ",infNormA,"\n",Indent(),
         "||A||_F                        = ",frobNormA,"\n",Indent(),
         "||A X - U^-H AOrig U^-1 X||_1  = ",oneNormError,"\n",Indent(),
         "||A X - U^-H AOrig U^-1 X||_oo = ",infNormError,"\n",Indent(),
         "||A X - U^-H AOrig U^-1 X||_F  = ",frobNormError);
    }
}

template<typename F,
         typename=EnableIf<IsBlasScalar<F>>>
void TestTwoSidedTrsm
( UpperOrLower uplo,
  UnitOrNonUnit diag,
  Int m,
  const Grid& g,
  bool scalapack,
  bool print,
  bool correctness )
{
    Output("Testing with ",TypeName<F>());
    PushIndent();

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

    Timer timer;
    if( scalapack )
    {
        OutputFromRoot(g.Comm(),"Starting ScaLAPACK TwoSidedTrsm");
        mpi::Barrier( g.Comm() );
        timer.Start();
        DistMatrix<F,MC,MR,BLOCK> ABlock( A ), BBlock( B );
        TwoSidedTrsm( uplo, diag, ABlock, BBlock );
        double runTime = timer.Stop();
        double gFlops = Pow(double(m),3.)/(runTime*1.e9);
        if( IsComplex<F>::value )
            gFlops *= 4.;
        OutputFromRoot(g.Comm(),"Time = ",runTime," seconds. GFlops = ",gFlops);
    }

    OutputFromRoot(g.Comm(),"Starting Elemental TwoSidedTrsm");
    mpi::Barrier( g.Comm() );
    timer.Start();
    TwoSidedTrsm( uplo, diag, A, B );
    mpi::Barrier( g.Comm() );
    double runTime = timer.Stop();
    double gFlops = Pow(double(m),3.)/(runTime*1.e9);
    if( IsComplex<F>::value )
        gFlops *= 4.;
    OutputFromRoot(g.Comm(),"Time = ",runTime," seconds. GFlops = ",gFlops);
    if( print )
        Print( A, "A after reduction" );
    if( correctness )
        TestCorrectness( uplo, diag, A, B, AOrig, print );

    PopIndent();
}

template<typename F>
void TestTwoSidedTrsm
( UpperOrLower uplo,
  UnitOrNonUnit diag,
  Int m,
  const Grid& g,
  bool print,
  bool correctness )
{
    OutputFromRoot(g.Comm(),"Testing with ",TypeName<F>());
    PushIndent();

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

    OutputFromRoot(g.Comm(),"Starting Elemental TwoSidedTrsm");
    mpi::Barrier( g.Comm() );
    Timer timer;
    timer.Start();
    TwoSidedTrsm( uplo, diag, A, B );
    mpi::Barrier( g.Comm() );
    double runTime = timer.Stop();
    double gFlops = Pow(double(m),3.)/(runTime*1.e9);
    if( IsComplex<F>::value )
        gFlops *= 4.;
    OutputFromRoot(g.Comm(),"Time = ",runTime," seconds. GFlops = ",gFlops);
    if( print )
        Print( A, "A after reduction" );
    if( correctness )
        TestCorrectness( uplo, diag, A, B, AOrig, print );

    PopIndent();
}

int
main( int argc, char* argv[] )
{
    Environment env( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;

    try
    {
        int gridHeight = Input("--gridHeight","height of process grid",0);
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

        if( gridHeight == 0 )
            gridHeight = Grid::DefaultHeight( mpi::Size(comm) );
        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( comm, gridHeight, order );
        const UpperOrLower uplo = CharToUpperOrLower( uploChar );
        const UnitOrNonUnit diag = CharToUnitOrNonUnit( diagChar );
        SetBlocksize( blocksize );
        SetDefaultBlockHeight( mb );
        SetDefaultBlockWidth( nb );

        ComplainIfDebug();
        OutputFromRoot(comm,"Will test TwoSidedTrsm",uploChar,diagChar);

        TestTwoSidedTrsm<float>
        ( uplo, diag, m, g, scalapack, print, correctness );
        TestTwoSidedTrsm<Complex<float>>
        ( uplo, diag, m, g, scalapack, print, correctness );

        TestTwoSidedTrsm<double>
        ( uplo, diag, m, g, scalapack, print, correctness );
        TestTwoSidedTrsm<Complex<double>>
        ( uplo, diag, m, g, scalapack, print, correctness );

#ifdef EL_HAVE_QD
       TestTwoSidedTrsm<DoubleDouble>
        ( uplo, diag, m, g, print, correctness );
       TestTwoSidedTrsm<QuadDouble>
        ( uplo, diag, m, g, print, correctness );

       TestTwoSidedTrsm<Complex<DoubleDouble>>
        ( uplo, diag, m, g, print, correctness );
       TestTwoSidedTrsm<Complex<QuadDouble>>
        ( uplo, diag, m, g, print, correctness );
#endif

#ifdef EL_HAVE_QUAD
        TestTwoSidedTrsm<Quad>
        ( uplo, diag, m, g, print, correctness );
        TestTwoSidedTrsm<Complex<Quad>>
        ( uplo, diag, m, g, print, correctness );
#endif

#ifdef EL_HAVE_MPC
        TestTwoSidedTrsm<BigFloat>
        ( uplo, diag, m, g, print, correctness );
        TestTwoSidedTrsm<Complex<BigFloat>>
        ( uplo, diag, m, g, print, correctness );
#endif
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
