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
void TestTrsv
( UpperOrLower uplo,
  Orientation orientation,
  UnitOrNonUnit diag,
  Int n,
  const Grid& g,
  bool print )
{
    if( g.Rank() == 0 )
        Output("Testing with ",TypeName<F>());
    typedef Base<F> Real;
    DistMatrix<F> A(g), x(g), y(g);

    // Generate random A and x
    HermitianUniformSpectrum( A, n, 1, 10 );
    Uniform( x, n, 1 );
    // Either y := op(L) x or y := op(U) x
    y = x;
    Trmm( LEFT, uplo, orientation, diag, F(1), A, y );

    if( print )
    {
        Print( A, "A" );
        Print( x, "x" );
        Print( y, "y" );
    }
    if( g.Rank() == 0 )
        Output("  Starting Trsv");
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    Trsv( uplo, orientation, diag, A, y );
    mpi::Barrier( g.Comm() );
    const double runTime = mpi::Time() - startTime;
    const double realGFlops = Pow(double(n),2.)/(1.e9*runTime);
    const double gFlops = ( IsComplex<F>::value ? 4*realGFlops : realGFlops );
    if( g.Rank() == 0 )
        Output("  Finished in ",runTime," seconds (",gFlops," GFlop/s)");
    if( print )
        Print( y, "y after solve" );

    y -= x;
    const Real xNorm = FrobeniusNorm( x );
    const Real yNorm = FrobeniusNorm( y );
    if( g.Rank() == 0 )
    {
        Output("  ||   x   ||_2 = ",xNorm);
        Output("  || x - y ||_2 = ",yNorm);
        Output("  || x - y ||_2 / || x ||_2 = ",yNorm/xNorm);
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
        const char uploChar = Input
            ("--uplo","upper or lower triangular: L/U",'L');
        const char transChar = Input
            ("--trans","orientation of triangular matrix: N/T/C",'N');
        const char diagChar = Input("--diag","(non-)unit diagonal: N/U",'N');
        const Int n = Input("--n","size of triangular matrix",100);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        if( r == 0 )
            r = Grid::FindFactor( commSize );
        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( comm, r, order );
        const UpperOrLower uplo = CharToUpperOrLower( uploChar );
        const Orientation orientation = CharToOrientation( transChar );
        const UnitOrNonUnit diag = CharToUnitOrNonUnit( diagChar );
        SetBlocksize( nb );

        ComplainIfDebug();
        if( commRank == 0 )
            Output("Will test Trsv ",uploChar,transChar,diagChar);

        TestTrsv<float>( uplo, orientation, diag, n, g, print );
        TestTrsv<Complex<float>>( uplo, orientation, diag, n, g, print );

        TestTrsv<double>( uplo, orientation, diag, n, g, print );
        TestTrsv<Complex<double>>( uplo, orientation, diag, n, g, print );

#ifdef EL_HAVE_QD
        TestTrsv<DoubleDouble>( uplo, orientation, diag, n, g, print );
        TestTrsv<QuadDouble>( uplo, orientation, diag, n, g, print );
#endif

#ifdef EL_HAVE_QUAD
        TestTrsv<Quad>( uplo, orientation, diag, n, g, print );
        TestTrsv<Complex<Quad>>( uplo, orientation, diag, n, g, print );
#endif

#ifdef EL_HAVE_MPC
        TestTrsv<BigFloat>( uplo, orientation, diag, n, g, print );
#endif
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
