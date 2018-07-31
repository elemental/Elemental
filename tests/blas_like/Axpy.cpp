/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
using namespace El;

template<typename T>
void TestCorrectness
( T alpha, const DistMatrix<T>& X,
           const DistMatrix<T>& YOrig,
           const DistMatrix<T>& YFinal,
  bool print )
{
    EL_DEBUG_ONLY(CallStackEntry cse("TestCorrectness"))
      
    const Int m = X.Height();
    const Int n = X.Width();
    const Grid& g = X.Grid();
    const Base<T> YFrobNorm = FrobeniusNorm( YOrig );
    if( g.Rank() == 0 )
    {
        Matrix<T> E( m, n );
        for( Int j=0; j<n; ++j )
        {
            for( Int i=0; i<m; ++i )
            {
                T Eij = alpha * X.Get(i,j) + YOrig.Get(i,j);
                Eij -= YFinal.Get(i,j);
                E.Set( i, j, Eij );
            }
        }
        const Base<T> EFrobNorm = FrobeniusNorm( E );
        if( print )
            Print( E, "E" );
        Output
        ( "|| E ||_F / || Y ||_F = ",
          EFrobNorm, "/", YFrobNorm, "=", EFrobNorm/YFrobNorm );
    }
}

template<typename T>
void TestAxpy
( Int m,
  Int n,
  Int ldimX,
  Int ldimY,
  Int numThreads,
  const Grid& g,
  bool print,
  bool correctness )
{
    OutputFromRoot(g.Comm(),"Testing with ",TypeName<T>());
    PushIndent();

    double runTime, gflops;
    Timer timer;

#ifdef EL_HYBRID
    if( numThreads > 0 )
    {
        omp_set_num_threads(numThreads);
    }
#endif

    // Generate random matrices
    DistMatrix<T> X(g), Y(g), YOrig(g), XData(g), YData(g);
    DistMatrix<T,STAR,STAR> alpha(g);
    Gaussian( alpha, 1, 1 );
    Gaussian( XData, ldimX, n );
    View( X, XData, IR(0,m), ALL );
    Gaussian( YData, ldimY, n );
    View( Y, YData, IR(0,m), ALL );
    Zeros( YOrig, m, n );
    Copy( Y, YOrig );
    if( print )
    {
        Print( alpha, "alpha" );
        Print( X, "X" );
        Print( Y, "Y" );
    }

    // Clear L3 cache
    Matrix<float> temp;
    double cacheSize = 1e7;
    temp.Resize( Int(cacheSize/sizeof(float)), 1 );
    Scale( 2.f, temp );

    // Apply axpy
    if( g.Rank() == 0 )
    {
        Output("Starting Axpy");
    }
    mpi::Barrier( g.Comm() );
    timer.Start();
    Axpy( alpha.Get(0,0), X, Y );
    mpi::Barrier( g.Comm() );
    runTime = timer.Stop();

    // Print results
    gflops = 2e-9 * m * n / runTime;
    OutputFromRoot(g.Comm(),"Finished in ",runTime," seconds (",gflops," GFLOPS)");
    if( print )
        Print( Y, "Y" );
    if( correctness )
        TestCorrectness( alpha.Get(0,0), X, YOrig, Y, print );

    PopIndent();
}

int
main( int argc, char* argv[] )
{
    Environment env( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;

    try
    {
        // Get command-line arguments
        int gridHeight = Input("--gridHeight","height of process grid",0);
        const Int m = Input("--m","height of matrix",100);
        const Int n = Input("--n","width of matrix",100);
        Int ldimX = Input("--ldimX","leading dimension of matrix X",100);
        Int ldimY = Input("--ldimY","leading dimension of matrix Y",100);
        const Int numThreads = Input("--numThreads","number of OpenMP threads (-1=default)",-1);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const bool colMajor = Input("--colMajor","column-major ordering?",true);
        const bool print = Input("--print","print matrices?",false);
        const bool correctness = Input("--correctness","correctness?",true);
        ProcessInput();
        PrintInputReport();

        // Set Elemental parameters
        if( gridHeight == 0 )
            gridHeight = Grid::DefaultHeight( mpi::Size(comm) );
        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( comm, gridHeight, order );
        SetBlocksize( nb );
        ldimX = Max(m, ldimX);
        ldimY = Max(m, ldimY);
        ComplainIfDebug();

        // Message
        OutputFromRoot(comm,"Testing Axpy");

        // Run tests
        TestAxpy<float>( m, n, ldimX, ldimY, numThreads, g, print, correctness );
        TestAxpy<Complex<float>>( m, n, ldimX, ldimY, numThreads, g, print, correctness );
        TestAxpy<double>( m, n, ldimX, ldimY, numThreads, g, print, correctness );
        TestAxpy<Complex<double>>( m, n, ldimX, ldimY, numThreads, g, print, correctness );
#ifdef EL_HAVE_QD
        TestAxpy<DoubleDouble>( m, n, ldimX, ldimY, numThreads, g, print, correctness );
        TestAxpy<QuadDouble>( m, n, ldimX, ldimY, numThreads, g, print, correctness );
        TestAxpy<Complex<DoubleDouble>>( m, n, ldimX, ldimY, numThreads, g, print, correctness );
        TestAxpy<Complex<QuadDouble>>( m, n, ldimX, ldimY, numThreads, g, print, correctness );
#endif
#ifdef EL_HAVE_QUAD
        TestAxpy<Quad>( m, n, ldimX, ldimY, numThreads, g, print, correctness );
        TestAxpy<Complex<Quad>>( m, n, ldimX, ldimY, numThreads, g, print, correctness );      
#endif
#ifdef EL_HAVE_MPC
        TestAxpy<BigFloat>( m, n, ldimX, ldimY, numThreads, g, print, correctness );
        TestAxpy<Complex<BigFloat>>( m, n, ldimX, ldimY, numThreads, g, print, correctness );
#endif
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
