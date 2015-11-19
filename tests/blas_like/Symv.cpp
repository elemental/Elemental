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
void TestSymv
( UpperOrLower uplo,
  Int m,
  T alpha,
  T beta, 
  bool print,
  const Grid& g )
{
    DistMatrix<T> A(g), x(g), y(g);

    Uniform( A, m, m );
    Uniform( x, m, 1 );
    Uniform( y, m, 1 ); 
    if( print )
    {
        Print( A, "A" );
        Print( x, "x" );
        Print( y, "y" );
    }

    // Test Symm
    if( g.Rank() == 0 )
        Output("  Starting Symv");
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    Symv( uplo, alpha, A, x, beta, y );
    mpi::Barrier( g.Comm() );
    const double runTime = mpi::Time() - startTime;
    const double realGFlops = 2.*double(m)*double(m)/(1.e9*runTime);
    const double gFlops = ( IsComplex<T>::value ? 4*realGFlops : realGFlops );
    if( g.Rank() == 0 )
        Output("  Finished in ",runTime," seconds (",gFlops," GFlop/s");
    if( print )
        Print( y, BuildString("y := ",alpha," Symm(A) x + ",beta," y") );
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
        const char uploChar = Input("--uplo","upper or lower storage: L/U",'L');
        const Int m = Input("--m","height of matrix",100);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const Int nbLocalDouble = Input
            ("--nbLocalDouble","local blocksize for real doubles",32);
        const Int nbLocalComplexDouble = Input
            ("--nbLocalComplexDouble","local blocksize for complex doubles",32);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        if( r == 0 )
            r = Grid::FindFactor( commSize );
        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( comm, r, order );
        const UpperOrLower uplo = CharToUpperOrLower( uploChar );
        SetBlocksize( nb );
        SetLocalSymvBlocksize<double>( nbLocalDouble );
        SetLocalSymvBlocksize<Complex<double>>( nbLocalComplexDouble );

        ComplainIfDebug();
        if( commRank == 0 )
            Output("Will test Symv ",uploChar);

        if( commRank == 0 )
            Output("Testing with doubles");
        TestSymv<double>( uplo, m, 3., 4., print, g );

        if( commRank == 0 )
            Output("Testing with Complex<double>");
        TestSymv<Complex<double>>
        ( uplo, m, Complex<double>(3), Complex<double>(4), print, g );
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
