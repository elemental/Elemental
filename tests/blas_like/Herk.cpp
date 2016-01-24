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
( UpperOrLower uplo, Orientation orientation,
  Base<T> alpha, const DistMatrix<T>& A,
  Base<T> beta,  const DistMatrix<T>& COrig,
                 const DistMatrix<T>& C,
  bool print )
{
    DEBUG_ONLY(CSE cse("TestCorrectness"))
    DistMatrix<T,CIRC,CIRC> ARoot( A ), 
                            COrigRoot( COrig ), CRoot( C );
    if( ARoot.Root() == ARoot.CrossRank() )
    {
        Matrix<T> CSeq( COrigRoot.Matrix() );
        Herk( uplo, orientation, alpha, ARoot.Matrix(), beta,  CSeq );
        const Base<T> CNrm = FrobeniusNorm( CRoot.Matrix() );
        CRoot.Matrix() -= CSeq;
        const Base<T> ENrm = FrobeniusNorm( CRoot.Matrix() );
        Output(" || E ||_F = ",ENrm);
        Output(" || C ||_F = ",CNrm);
    }
}

template<typename T> 
void TestHerk
( UpperOrLower uplo,
  Orientation orientation,
  Int m,
  Int k,
  Base<T> alpha,
  Base<T> beta,
  const Grid& g,
  bool print,
  bool correctness,
  Int nbLocal,
  Int colAlignA=0,
  Int rowAlignA=0,
  Int colAlignC=0,
  Int rowAlignC=0 )
{
    if( g.Rank() == 0 )
        Output("Testing with ",TypeName<T>());
    SetLocalTrrkBlocksize<T>( nbLocal );

    DistMatrix<T> A(g), C(g);
    A.Align( colAlignA, rowAlignA );
    C.Align( colAlignC, rowAlignC );

    if( orientation == NORMAL )
        Uniform( A, m, k );
    else
        Uniform( A, k, m );
    HermitianUniformSpectrum( C, m, 1, 10 );
    auto COrig = C;
    if( print )
    {
        Print( A, "A" );
        Print( C, "C" );
    }

    if( g.Rank() == 0 )
        Output("  Starting Herk");
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    Herk( uplo, orientation, alpha, A, beta, C );
    mpi::Barrier( g.Comm() );
    const double runTime = mpi::Time() - startTime;
    const double realGFlops = double(m)*double(m)*double(k)/(1.e9*runTime);
    const double gFlops = ( IsComplex<T>::value ? 4*realGFlops : realGFlops );
    if( g.Rank() == 0 )
        Output("  Finished in ",runTime," seconds (",gFlops," GFlop/s)");
    if( print )
    {
        if( orientation == NORMAL )
            Print( C, BuildString("C := ",alpha," A A' + ",beta," C") );
        else
            Print( C, BuildString("C := ",alpha," A' A + ",beta," C") );
    }

    if( correctness )
        TestCorrectness( uplo, orientation, alpha, A, beta, COrig, C, print );
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
        const char uploChar = Input("--uplo","upper/lower storage: L/U",'L');
        const char transChar = Input("--trans","orientation: N/C",'N');
        const Int m = Input("--m","height of result",100);
        const Int k = Input("--k","inner dimension",100);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const Int nbLocal = Input("--nbLocal","local blocksize",32);
        const bool print = Input("--print","print matrices?",false);
        const bool correctness = Input("--correctness","test correct?",true);
        const Int colAlignA = Input("--colAlignA","col align of A",0);
        const Int colAlignC = Input("--colAlignC","col align of C",0);
        const Int rowAlignA = Input("--rowAlignA","row align of A",0);
        const Int rowAlignC = Input("--rowAlignC","row align of C",0);
        ProcessInput();
        PrintInputReport();

        if( r == 0 )
            r = Grid::FindFactor( commSize );
        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( comm, r, order );
        const UpperOrLower uplo = CharToUpperOrLower( uploChar );
        const Orientation orientation = CharToOrientation( transChar );
        SetBlocksize( nb );

        ComplainIfDebug();
        if( commRank == 0 )
            Output("Will test Herk ",uploChar,transChar);

        TestHerk<float>
        ( uplo, orientation, m, k,
          float(3), float(4),
          g, print, correctness,
          nbLocal,
          colAlignA, rowAlignA, colAlignC, rowAlignC );
        TestHerk<Complex<float>>
        ( uplo, orientation, m, k,
          float(3), float(4),
          g, print, correctness,
          nbLocal,
          colAlignA, rowAlignA, colAlignC, rowAlignC );

        TestHerk<double>
        ( uplo, orientation, m, k,
          double(3), double(4),
          g, print, correctness,
          nbLocal,
          colAlignA, rowAlignA, colAlignC, rowAlignC );
        TestHerk<Complex<double>>
        ( uplo, orientation, m, k,
          double(3), double(4),
          g, print, correctness,
          nbLocal,
          colAlignA, rowAlignA, colAlignC, rowAlignC );

#ifdef EL_HAVE_QD
        TestHerk<DoubleDouble>
        ( uplo, orientation, m, k,
          DoubleDouble(3), DoubleDouble(4),
          g, print, correctness,
          nbLocal,
          colAlignA, rowAlignA, colAlignC, rowAlignC );
        TestHerk<QuadDouble>
        ( uplo, orientation, m, k,
          QuadDouble(3), QuadDouble(4),
          g, print, correctness,
          nbLocal,
          colAlignA, rowAlignA, colAlignC, rowAlignC );
#endif

#ifdef EL_HAVE_QUAD
        TestHerk<Quad>
        ( uplo, orientation, m, k,
          Quad(3), Quad(4),
          g, print, correctness,
          nbLocal,
          colAlignA, rowAlignA, colAlignC, rowAlignC );
        TestHerk<Complex<Quad>>
        ( uplo, orientation, m, k,
          Quad(3), Quad(4),
          g, print, correctness,
          nbLocal,
          colAlignA, rowAlignA, colAlignC, rowAlignC );
#endif

#ifdef EL_HAVE_MPC
        TestHerk<BigFloat>
        ( uplo, orientation, m, k,
          BigFloat(3), BigFloat(4),
          g, print, correctness,
          nbLocal,
          colAlignA, rowAlignA, colAlignC, rowAlignC );
#endif
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
