/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

template<typename F>
void TestCorrectness
( bool print,
  const ElementalMatrix<F>& A,
  const ElementalMatrix<F>& Z )
{
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    const Int n = Z.Height();
    const Int k = Z.Width();
    
    // X := AZ
    DistMatrix<F> X( Z );
    Trmm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), A, X );
    // Find the residual X-ZW = AZ-ZW
    DistMatrix<F> ZW( Z );
    DistMatrix<F> w( g );
    GetDiagonal( A, w );
    DiagonalScale( RIGHT, NORMAL, w, ZW );
    X -= ZW;
    // Find the Frobenius norms of A and AZ-ZW
    Real frobNormA = FrobeniusNorm( A );
    Real frobNormE = FrobeniusNorm( X );
    if( g.Rank() == 0 )
        Output("    ||A Z - Z W||_F / ||A||_F = ",frobNormE/frobNormA);
}

template<typename F,Dist U=MC,Dist V=MR,Dist S=MC>
void TestTriangEig
( bool testCorrectness,
  bool print,
  Int m, 
  const Grid& g,
  bool scalapack )
{
    typedef Base<F> Real;
    DistMatrix<F,U,V> A(g), AOrig(g), Z(g);
    DistMatrix<F,S,STAR> w(g);
    if( g.Rank() == 0 )
        Output("Testing with ",TypeName<F>());

    // TODO: option to construct matrix with repeated eigenvalues
    DistMatrix<F,U,V> B(g);
    Gaussian( B, m, m );
    LU( B );
    Transpose( B, A );
    MakeTrapezoidal( UPPER, A, 0 );
    if( testCorrectness ) {
        AOrig = A;
	GetDiagonal( A, w );
    }
    if( print )
        Print( A, "A" );

    if( scalapack && U == MC && V == MR )
    {
#if 0 // TODO: option to test ScaLAPACK. I am not sure why this does
      // not compile successfully.
      DistMatrix<F,MC,MR,BLOCK> ABlock( A ), ZBlock(g);
      const double startTime = mpi::Time();
      TriangEig( ABlock, ZBlock );
      const double runTime = mpi::Time() - startTime;
      if( g.Rank() == 0 )
	Output("  ScaLAPACK TriangEig time: ",runTime);
#endif
    }

    if( g.Rank() == 0 )
        Output("  Starting triangular eigensolver...");
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    TriangEig( A, Z );
    mpi::Barrier( g.Comm() );
    const double runTime = mpi::Time() - startTime;
    if( g.Rank() == 0 )
        Output("  Time = ",runTime," seconds");
    if( print )
    {
        Print( w, "eigenvalues:" );
	Print( Z, "eigenvectors:" );
    }
    if( testCorrectness )
        TestCorrectness( print, AOrig, Z );
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
        Int r = Input("--gridHeight","height of process grid",0);
        const bool colMajor = Input("--colMajor","column-major ordering?",true);
        const Int n = Input("--height","height of matrix",100);
        const Int nb = Input("--nb","algorithmic blocksize",96);
#ifdef EL_HAVE_SCALAPACK
        const bool scalapack = Input("--scalapack","test ScaLAPACK?",true);
#else
        const bool scalapack = false;
#endif
        const bool testCorrectness = Input
            ("--correctness","test correctness?",true);
        const bool print = Input("--print","print matrices?",false);
        const bool testReal = Input("--testReal","test real matrices?",true);
        const bool testCpx = Input("--testCpx","test complex matrices?",true);
        ProcessInput();
        PrintInputReport();

        if( r == 0 )
            r = Grid::FindFactor( commSize );
        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( comm, r, order );
        SetBlocksize( nb );
        ComplainIfDebug();

        if( commRank == 0 )
            Output("Normal algorithms:");
        if( testReal )
            TestTriangEig<double>
            ( testCorrectness, print,
              n, g, scalapack );
        if( testCpx )
            TestTriangEig<Complex<double>>
            ( testCorrectness, print,
              n, g, scalapack );

        // Also test with non-standard distributions
        if( commRank == 0 )
            Output("Nonstandard distributions:");
        if( testReal )
            TestTriangEig<double,MR,MC,MC>
            ( testCorrectness, print,
              n, g, scalapack );
        if( testCpx )
            TestTriangEig<Complex<double>,MR,MC,MC>
            ( testCorrectness, print,
              n, g, scalapack );
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
