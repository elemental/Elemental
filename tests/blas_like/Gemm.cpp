/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

template<typename T>
void TestCorrectness
( Orientation orientA, Orientation orientB,
  T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& B, 
  T beta,  const DistMatrix<T>& COrig, 
           const DistMatrix<T>& CFinal, 
  bool print )
{
    DEBUG_ONLY(CallStackEntry cse("TestCorrectness"))
    DistMatrix<T,CIRC,CIRC> ARoot( A ), BRoot( B ), 
                            COrigRoot( COrig ), CFinalRoot( CFinal );
    if( ARoot.Root() == ARoot.CrossRank() )
    {
        Matrix<T> CSeq( COrigRoot.Matrix() );
        Gemm
        ( orientA, orientB, 
          alpha, ARoot.Matrix(), BRoot.Matrix(),
          beta,  CSeq );
        const Base<T> CNrm = FrobeniusNorm( CFinalRoot.Matrix() );
        CFinalRoot.Matrix() -= CSeq;
        const Base<T> ENrm = FrobeniusNorm( CFinalRoot.Matrix() );
        Output(" || E ||_F = ",ENrm);
        Output(" || C ||_F = ",CNrm);
    }
}

template<typename T> 
void TestGemm
( Orientation orientA,
  Orientation orientB,
  Int m,
  Int n,
  Int k,
  T alpha,
  T beta,
  const Grid& g, 
  bool print,
  bool correctness,
  Int colAlignA=0, Int rowAlignA=0,
  Int colAlignB=0, Int rowAlignB=0,
  Int colAlignC=0, Int rowAlignC=0 )
{
    if( g.Rank() == 0 )
        Output("Testing with ",TypeName<T>());
    double startTime, runTime, realGFlops, gFlops;
    DistMatrix<T> A(g), B(g), COrig(g), C(g);

    A.Align( colAlignA, rowAlignA );
    B.Align( colAlignB, rowAlignB );
    C.Align( colAlignC, rowAlignC );

    if( orientA == NORMAL )
        Uniform( A, m, k );
    else
        Uniform( A, k, m );
    if( orientB == NORMAL )
        Uniform( B, k, n );
    else
        Uniform( B, n, k );
    Uniform( COrig, m, n );
    if( print )
    {
        Print( A, "A" );
        Print( B, "B" );
        Print( COrig, "COrig" );
    }

    // Test the variant of Gemm that keeps A stationary
    C = COrig;
    if( g.Rank() == 0 )
        Output("Stationary A algorithm:");
    mpi::Barrier( g.Comm() );
    startTime = mpi::Time();
    Gemm( orientA, orientB, alpha, A, B, beta, C, GEMM_SUMMA_A );
    mpi::Barrier( g.Comm() );
    runTime = mpi::Time() - startTime;
    realGFlops = 2.*double(m)*double(n)*double(k)/(1.e9*runTime);
    gFlops = ( IsComplex<T>::value ? 4*realGFlops : realGFlops );
    if( g.Rank() == 0 )
        Output("  Finished in ",runTime," seconds (",gFlops," GFlop/s)");
    if( print )
        Print( C, BuildString("C := ",alpha," A B + ",beta," C") );
    if( correctness )
        TestCorrectness( orientA, orientB, alpha, A, B, beta, COrig, C, print );

    // Test the variant of Gemm that keeps B stationary
    C = COrig;
    if( g.Rank() == 0 )
        Output("Stationary B Algorithm:");
    mpi::Barrier( g.Comm() );
    startTime = mpi::Time();
    Gemm( orientA, orientB, alpha, A, B, beta, C, GEMM_SUMMA_B );
    mpi::Barrier( g.Comm() );
    runTime = mpi::Time() - startTime;
    realGFlops = 2.*double(m)*double(n)*double(k)/(1.e9*runTime);
    gFlops = ( IsComplex<T>::value ? 4*realGFlops : realGFlops );
    if( g.Rank() == 0 )
        Output("  Finished in ",runTime," seconds (",gFlops," GFlop/s)");
    if( print )
        Print( C, BuildString("C := ",alpha," A B + ",beta," C") );
    if( correctness )
        TestCorrectness( orientA, orientB, alpha, A, B, beta, COrig, C, print );

    // Test the variant of Gemm that keeps C stationary
    C = COrig;
    if( g.Rank() == 0 )
        Output("Stationary C Algorithm:");
    mpi::Barrier( g.Comm() );
    startTime = mpi::Time();
    Gemm( orientA, orientB, alpha, A, B, beta, C, GEMM_SUMMA_C );
    mpi::Barrier( g.Comm() );
    runTime = mpi::Time() - startTime;
    realGFlops = 2.*double(m)*double(n)*double(k)/(1.e9*runTime);
    gFlops = ( IsComplex<T>::value ? 4*realGFlops : realGFlops );
    if( g.Rank() == 0 )
        Output("  Finished in ",runTime," seconds (",gFlops," GFlop/s)");
    if( print )
        Print( C, BuildString("C := ",alpha," A B + ",beta," C") );
    if( correctness )
        TestCorrectness( orientA, orientB, alpha, A, B, beta, COrig, C, print );
    
    if( orientA == NORMAL && orientB == NORMAL )
    {
        // Test the variant of Gemm for panel-panel dot products
        if( g.Rank() == 0 )
            Output("Dot Product Algorithm:");
        C = COrig;
        mpi::Barrier( g.Comm() );
        startTime = mpi::Time();
        Gemm( NORMAL, NORMAL, alpha, A, B, beta, C, GEMM_SUMMA_DOT );
        mpi::Barrier( g.Comm() );
        runTime = mpi::Time() - startTime;
        realGFlops = 2.*double(m)*double(n)*double(k)/(1.e9*runTime);
        gFlops = ( IsComplex<T>::value ? 4*realGFlops : realGFlops );
        if( g.Rank() == 0 )
            Output("  Finished in ",runTime," seconds (",gFlops," GFlop/s)");
        if( print )
            Print( C, BuildString("C := ",alpha," A B + ",beta," C") );
        if( correctness )
            TestCorrectness
            ( orientA, orientB, alpha, A, B, beta, COrig, C, print );
    }
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
        const bool colMajor = Input("--colMajor","column-major ordering?",true);
        Int r = Input("--r","height of process grid",0);
        const char transA = Input("--transA","orientation of A: N/T/C",'N');
        const char transB = Input("--transB","orientation of B: N/T/C",'N');
        const Int m = Input("--m","height of result",100);
        const Int n = Input("--n","width of result",100);
        const Int k = Input("--k","inner dimension",100);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const bool print = Input("--print","print matrices?",false);
        const bool correctness = Input("--correctness","correctness?",true);
        const Int colAlignA = Input("--colAlignA","column align of A",0);
        const Int colAlignB = Input("--colAlignB","column align of B",0);
        const Int colAlignC = Input("--colAlignC","column align of C",0);
        const Int rowAlignA = Input("--rowAlignA","row align of A",0);
        const Int rowAlignB = Input("--rowAlignB","row align of B",0); 
        const Int rowAlignC = Input("--rowAlignC","row align of C",0);
        ProcessInput();
        PrintInputReport();

        if( r == 0 )
            r = Grid::FindFactor( commSize );
        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( comm, r, order );
        const Orientation orientA = CharToOrientation( transA );
        const Orientation orientB = CharToOrientation( transB );
        SetBlocksize( nb );

        ComplainIfDebug();
        if( commRank == 0 )
            Output("Will test Gemm",transA,transB);

        TestGemm<float>
        ( orientA, orientB,
          m, n, k,
          float(3), float(4),
          g,
          print, correctness,
          colAlignA, rowAlignA,
          colAlignB, rowAlignB,
          colAlignC, rowAlignC );
        TestGemm<Complex<float>>
        ( orientA, orientB,
          m, n, k, 
          Complex<float>(3), Complex<float>(4),
          g,
          print, correctness,
          colAlignA, rowAlignA,
          colAlignB, rowAlignB,
          colAlignC, rowAlignC );

        TestGemm<double>
        ( orientA, orientB,
          m, n, k,
          double(3), double(4),
          g,
          print, correctness,
          colAlignA, rowAlignA,
          colAlignB, rowAlignB,
          colAlignC, rowAlignC );
        TestGemm<Complex<double>>
        ( orientA, orientB,
          m, n, k, 
          Complex<double>(3), Complex<double>(4),
          g,
          print, correctness,
          colAlignA, rowAlignA,
          colAlignB, rowAlignB,
          colAlignC, rowAlignC );

#ifdef EL_HAVE_QD
        TestGemm<DoubleDouble>
        ( orientA, orientB,
          m, n, k,
          DoubleDouble(3), DoubleDouble(4),
          g,
          print, correctness,
          colAlignA, rowAlignA,
          colAlignB, rowAlignB,
          colAlignC, rowAlignC );
        TestGemm<QuadDouble>
        ( orientA, orientB,
          m, n, k,
          QuadDouble(3), QuadDouble(4),
          g,
          print, correctness,
          colAlignA, rowAlignA,
          colAlignB, rowAlignB,
          colAlignC, rowAlignC );
#endif

#ifdef EL_HAVE_QUAD
        TestGemm<Quad>
        ( orientA, orientB,
          m, n, k,
          Quad(3), Quad(4),
          g,
          print, correctness,
          colAlignA, rowAlignA,
          colAlignB, rowAlignB,
          colAlignC, rowAlignC );
        TestGemm<Complex<Quad>>
        ( orientA, orientB,
          m, n, k, 
          Complex<Quad>(3), Complex<Quad>(4),
          g,
          print, correctness,
          colAlignA, rowAlignA,
          colAlignB, rowAlignB,
          colAlignC, rowAlignC );
#endif

#ifdef EL_HAVE_MPC
        TestGemm<BigFloat>
        ( orientA, orientB,
          m, n, k,
          BigFloat(3), BigFloat(4),
          g,
          print, correctness,
          colAlignA, rowAlignA,
          colAlignB, rowAlignB,
          colAlignC, rowAlignC );
#endif
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
