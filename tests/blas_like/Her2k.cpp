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
void TestHer2k
( UpperOrLower uplo,
  Orientation orientation,
  Int m,
  Int k,
  T alpha,
  Base<T> beta,
  const Grid& g,
  bool print,
  Int nbLocal )
{
    if( g.Rank() == 0 )
        Output("Testing with ",TypeName<T>());
    SetLocalTrr2kBlocksize<T>( nbLocal );

    DistMatrix<T> A(g), B(g), C(g);

    if( orientation == NORMAL )
    {
        Uniform( A, m, k );
        Uniform( B, m, k );
    }
    else
    {
        Uniform( A, k, m );
        Uniform( B, k, m );
    }
    HermitianUniformSpectrum( C, m, 1, 10 );
    if( print )
    {
        Print( A, "A" );
        Print( B, "B" );
        Print( C, "C" );
    }

    if( g.Rank() == 0 )
        Output("  Starting Her2k");
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    Her2k( uplo, orientation, alpha, A, B, beta, C );
    mpi::Barrier( g.Comm() );
    const double runTime = mpi::Time() - startTime;
    const double realGFlops = 2.*double(m)*double(m)*double(k)/(1.e9*runTime);
    const double gFlops = ( IsComplex<T>::value ? 4*realGFlops : realGFlops );
    if( g.Rank() == 0 )
        Output("  Finished in ",runTime," seconds (",gFlops," GFlop/s)");
    if( print )
    {
        if( orientation == NORMAL )
            Print( C, BuildString("C := ",alpha,"(A B' + B A') + ",beta," C") );
        else
            Print( C, BuildString("C := ",alpha,"(A' B + B' A) + ",beta," C") );
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
        Int r = Input("--r","height of process grid",0);
        const bool colMajor = Input("--colMajor","column-major ordering?",true);
        const char uploChar = Input("--uplo","upper/lower storage: L/U",'L');
        const char transChar = Input("--trans","orientation: N/C",'N');
        const Int m = Input("--m","height of result",100);
        const Int k = Input("--k","inner dimension",100);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const Int nbLocal = Input("--nbLocal","local blocksize",32);
        const bool print = Input("--print","print matrices?",false);
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
            Output("Will test Her2k ",uploChar,transChar);

        TestHer2k<float>
        ( uplo, orientation,
          m, k,
          float(3), float(4),
          g, print, nbLocal );
        TestHer2k<Complex<float>>
        ( uplo, orientation,
          m, k,
          Complex<float>(3), float(4),
          g, print, nbLocal );

        TestHer2k<double>
        ( uplo, orientation,
          m, k,
          double(3), double(4),
          g, print, nbLocal );
        TestHer2k<Complex<double>>
        ( uplo, orientation,
          m, k,
          Complex<double>(3), double(4),
          g, print, nbLocal );

#ifdef EL_HAVE_QD
        TestHer2k<DoubleDouble>
        ( uplo, orientation,
          m, k,
          DoubleDouble(3), DoubleDouble(4),
          g, print, nbLocal );
        TestHer2k<QuadDouble>
        ( uplo, orientation,
          m, k,
          QuadDouble(3), QuadDouble(4),
          g, print, nbLocal );
#endif

#ifdef EL_HAVE_QUAD
        TestHer2k<Quad>
        ( uplo, orientation,
          m, k,
          Quad(3), Quad(4),
          g, print, nbLocal );
        TestHer2k<Complex<Quad>>
        ( uplo, orientation,
          m, k,
          Complex<Quad>(3), Quad(4),
          g, print, nbLocal );
#endif

#ifdef EL_HAVE_MPC
        TestHer2k<BigFloat>
        ( uplo, orientation,
          m, k,
          BigFloat(3), BigFloat(4),
          g, print, nbLocal );
#endif
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
