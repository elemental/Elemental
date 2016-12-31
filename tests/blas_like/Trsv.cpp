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
void TestTrsv
( UpperOrLower uplo,
  Orientation orientation,
  UnitOrNonUnit diag,
  Int n,
  const Grid& g,
  bool print )
{
    OutputFromRoot(g.Comm(),"Testing with ",TypeName<F>());
    PushIndent();

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
    OutputFromRoot(g.Comm(),"Starting Trsv");
    mpi::Barrier( g.Comm() );
    Timer timer;
    timer.Start();
    Trsv( uplo, orientation, diag, A, y );
    mpi::Barrier( g.Comm() );
    const double runTime = timer.Stop();
    const double realGFlops = Pow(double(n),2.)/(1.e9*runTime);
    const double gFlops = ( IsComplex<F>::value ? 4*realGFlops : realGFlops );
    OutputFromRoot
    (g.Comm(),"Finished in ",runTime," seconds (",gFlops," GFlop/s)");
    if( print )
        Print( y, "y after solve" );

    y -= x;
    const Real xNorm = FrobeniusNorm( x );
    const Real yNorm = FrobeniusNorm( y );
    OutputFromRoot
     (g.Comm(),
      "||   x   ||_2 = ",xNorm,"\n",Indent(),
      "|| x - y ||_2 = ",yNorm,"\n",Indent(),
      "|| x - y ||_2 / || x ||_2 = ",yNorm/xNorm);

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
            ("--uplo","upper or lower triangular: L/U",'L');
        const char transChar = Input
            ("--trans","orientation of triangular matrix: N/T/C",'N');
        const char diagChar = Input("--diag","(non-)unit diagonal: N/U",'N');
        const Int n = Input("--n","size of triangular matrix",100);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        if( gridHeight == 0 )
            gridHeight = Grid::DefaultHeight( mpi::Size(comm) );
        const GridOrder order = colMajor ? COLUMN_MAJOR : ROW_MAJOR;
        const Grid g( comm, gridHeight, order );
        const UpperOrLower uplo = CharToUpperOrLower( uploChar );
        const Orientation orientation = CharToOrientation( transChar );
        const UnitOrNonUnit diag = CharToUnitOrNonUnit( diagChar );
        SetBlocksize( nb );

        ComplainIfDebug();
        OutputFromRoot(comm,"Will test Trsv ",uploChar,transChar,diagChar);

        TestTrsv<float>( uplo, orientation, diag, n, g, print );
        TestTrsv<Complex<float>>( uplo, orientation, diag, n, g, print );

        TestTrsv<double>( uplo, orientation, diag, n, g, print );
        TestTrsv<Complex<double>>( uplo, orientation, diag, n, g, print );

#ifdef EL_HAVE_QD
        TestTrsv<DoubleDouble>( uplo, orientation, diag, n, g, print );
        TestTrsv<QuadDouble>( uplo, orientation, diag, n, g, print );

        TestTrsv<Complex<DoubleDouble>>( uplo, orientation, diag, n, g, print );
        TestTrsv<Complex<QuadDouble>>( uplo, orientation, diag, n, g, print );
#endif

#ifdef EL_HAVE_QUAD
        TestTrsv<Quad>( uplo, orientation, diag, n, g, print );
        TestTrsv<Complex<Quad>>( uplo, orientation, diag, n, g, print );
#endif

#ifdef EL_HAVE_MPC
        TestTrsv<BigFloat>( uplo, orientation, diag, n, g, print );
        TestTrsv<Complex<BigFloat>>( uplo, orientation, diag, n, g, print );
#endif
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
