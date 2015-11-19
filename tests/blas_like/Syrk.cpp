/*
   Copyright (c) 2009-2015, Jack Poulson
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
  T alpha, const DistMatrix<T>& A,
  T beta,  const DistMatrix<T>& COrig,
           const DistMatrix<T>& C,
  bool print )
{
    DEBUG_ONLY(CSE cse("TestCorrectness"))
    DistMatrix<T,CIRC,CIRC> ARoot( A ),
                            COrigRoot( COrig ), CRoot( C );
    if( ARoot.Root() == ARoot.CrossRank() )
    {
        Matrix<T> CSeq( COrigRoot.Matrix() );
        Syrk
        ( uplo, orientation, alpha, ARoot.Matrix(), beta,  CSeq );
        const Base<T> CNrm = FrobeniusNorm( CRoot.Matrix() );
        CRoot.Matrix() -= CSeq;
        const Base<T> ENrm = FrobeniusNorm( CRoot.Matrix() );
        Output(" || E ||_F = ",ENrm);
        Output(" || C ||_F = ",CNrm);
    }
}

template<typename T>
void TestSyrk
( UpperOrLower uplo, Orientation orientation,
  Int m, Int k, T alpha, T beta, const Grid& g, bool print, bool correctness,
  Int colAlignA=0, Int rowAlignA=0, Int colAlignC=0, Int rowAlignC=0,
  bool contigA=true, bool contigC=true )
{
    DistMatrix<T> A(g), C(g);
    A.Align( colAlignA, rowAlignA );
    C.Align( colAlignC, rowAlignC );

    if( orientation == NORMAL )
    {
        if( !contigA )
        {
            A.Resize( m, k );
            A.Resize( m, k, 2*A.LDim() );
        }
        Uniform( A, m, k );
    }
    else
    {
        if( !contigA )
        {
            A.Resize( k, m );
            A.Resize( k, m, 2*A.LDim() );
        }
        Uniform( A, k, m );
    }
    if( !contigC )
    {
        C.Resize( m, m );
        C.Resize( m, m, 2*C.LDim() );
    }
    Uniform( C, m, m );
    MakeTrapezoidal( uplo, C );
    auto COrig = C;
    if( print )
    {
        Print( A, "A" );
        Print( C, "C" );
    }

    if( g.Rank() == 0 )
        Output("  Starting Syrk");
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    Syrk( uplo, orientation, alpha, A, beta, C );
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
        const bool contigA = Input("--contigA","contiguous A?",true);
        const bool contigC = Input("--contigC","contiguous C?",true);
        ProcessInput();
        PrintInputReport();

        if( r == 0 )
            r = Grid::FindFactor( commSize );
        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( comm, r, order );
        const UpperOrLower uplo = CharToUpperOrLower( uploChar );
        const Orientation orientation = CharToOrientation( transChar );
        SetBlocksize( nb );
        SetLocalTrrkBlocksize<double>( nbLocal );
        SetLocalTrrkBlocksize<Complex<double>>( nbLocal );

        ComplainIfDebug();
        if( commRank == 0 )
            Output("Will test Syrk ",uploChar,transChar);

        if( commRank == 0 )
            Output("Testing with doubles");
        TestSyrk<double>
        ( uplo, orientation, m, k, 3., 4., g, print, correctness, 
          colAlignA, rowAlignA, colAlignC, rowAlignC, contigA, contigC );

        if( commRank == 0 )
            Output("Testing with Complex<double>");
        TestSyrk<Complex<double>>
        ( uplo, orientation, m, k, Complex<double>(3), Complex<double>(4), g,
          print, correctness, colAlignA, rowAlignA, colAlignC, rowAlignC,
          contigA, contigC );
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
