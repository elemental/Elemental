/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

// Tests with local matrices
// TODO: replace this with tests for distributed matrices
#if 1

template<typename F>
void TestCorrectness
( bool print,
  const Matrix<F>& A,
  const Matrix<F>& X )
{
    typedef Base<F> Real;
    const Int n = X.Height();
    const Int k = X.Width();
    
    // Find the residual R = AX-XW
    Matrix<F> R( X );
    Trmm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), A, R );
    Matrix<F> XW( X );
    Matrix<F> w;
    GetDiagonal( A, w );
    DiagonalScale( RIGHT, NORMAL, w, XW );
    R -= XW;
    // Find the Frobenius norms of A and AX-XW
    Real frobNormA = FrobeniusNorm( A );
    Real frobNormR = FrobeniusNorm( R );
    // Find condition number
    Real condX = Condition( X );
    Output("    ||A X - X W||_F / ||A||_F = ",frobNormR/frobNormA);
    Output("    cond(X) = ", condX);
      
}

template<typename F>
void TestTriangEig
( bool testCorrectness,
  bool print,
  Int m, 
  string testMatrix )
{
    Matrix<F> A, AOrig, X;
    Matrix<F> w;
    Output("Testing with ",TypeName<F>());

    for( Int i=0; i<testMatrix.length(); ++i )
    {
        testMatrix[i] = std::tolower( testMatrix[i] );
    }
    
    // Generate test matrix
    if( testMatrix == "grcar" )
    {
        Matrix<Complex<Base<F>>> d;
	Grcar( A, m );
	Schur( A, d );
	MakeTrapezoidal( UPPER, A, 0 );
    }
    else if( testMatrix == "repeated" )
    {
        // LU factorization of Gaussian matrix
        Matrix<F> B;
	Gaussian( B, m, m );
	LU( B );
	Transpose( B, A );
	MakeTrapezoidal( UPPER, A, 0 );

	// 4/5 of eigenvalues are repeated
	F repeatList[4];
	for(Int i=0; i<Min(4,m); ++i)
	{
	    repeatList[i] = A.Get(i,i);
	}
	for(Int i=0; i<m; ++i)
	{
	    if( i%5 < 4 )
	    {
	        A.Set(i,i,repeatList[i%5]);
	    }
	}    
    }
    else
    {
        // LU factorization of Gaussian matrix
        Matrix<F> B;
	Gaussian( B, m, m );
	LU( B );
	Transpose( B, A );
	MakeTrapezoidal( UPPER, A, 0 );
    }
    
    if( testCorrectness )
    {
        AOrig = A;
	GetDiagonal( A, w );
    }
    if( print )
        Print( A, "A" );

    Output("  Starting triangular eigensolver...");
    const double startTime = mpi::Time();
    TriangEig( A, X );
    const double runTime = mpi::Time() - startTime;
    Output("  Time = ",runTime," seconds");
    if( print )
    {
        Print( w, "eigenvalues:" );
	Print( X, "eigenvectors:" );
    }
    if( testCorrectness )
        TestCorrectness( print, AOrig, X );
}

int 
main( int argc, char* argv[] )
{
    Environment env( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;

    try
    {
        const bool colMajor = Input("--colMajor","column-major ordering?",true);
        const Int n = Input("--height","height of matrix",100);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const bool testCorrectness = Input
            ("--correctness","test correctness?",true);
        const bool print = Input("--print","print matrices?",false);
        const bool testReal = Input("--testReal","test real matrices?",true);
        const bool testCpx = Input("--testCpx","test complex matrices?",true);
	std::string testMatrix = Input("--testMatrix","test matrix (gaussian, repeated, grcar)","gaussian");
        ProcessInput();
        PrintInputReport();

        SetBlocksize( nb );
        ComplainIfDebug();

        if( testReal )
            TestTriangEig<double>
            ( testCorrectness, print,
              n, testMatrix );
        if( testCpx )
            TestTriangEig<Complex<double>>
            ( testCorrectness, print,
              n, testMatrix );
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}

// Tests with distributed matrices
// TODO: restore this
#else

template<typename F>
void TestCorrectness
( bool print,
  const ElementalMatrix<F>& A,
  const ElementalMatrix<F>& X )
{
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    const Int n = X.Height();
    const Int k = X.Width();
    
    // Find the residual R = AX-XW
    DistMatrix<F> R( X );
    Trmm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), A, R );
    DistMatrix<F> XW( X );
    DistMatrix<F> w( g );
    GetDiagonal( A, w );
    DiagonalScale( RIGHT, NORMAL, w, XW );
    R -= XW;
    // Find the Frobenius norms of A and AX-XW
    Real frobNormA = FrobeniusNorm( A );
    Real frobNormR = FrobeniusNorm( R );
    // Find condition number
    Real condX = Condition( X );
    if( g.Rank() == 0 ) {
        Output("    ||A X - X W||_F / ||A||_F = ",frobNormR/frobNormA);
	Output("    cond(X) = ", condX);
    }
      
}

template<typename F,Dist U=MC,Dist V=MR,Dist S=MC>
void TestTriangEig
( bool testCorrectness,
  bool print,
  Int m, 
  const Grid& g,
  bool repeated )
{
    typedef Base<F> Real;
    DistMatrix<F,U,V> A(g), AOrig(g), X(g);
    DistMatrix<F,S,STAR> w(g);
    if( g.Rank() == 0 )
        Output("Testing with ",TypeName<F>());

    // Random triangular matrix is LU factor of Gaussian matrix
    DistMatrix<F,U,V> B(g);
    Gaussian( B, m, m );
    LU( B );
    Transpose( B, A );
    MakeTrapezoidal( UPPER, A, 0 );
    if( repeated )
    {
        F repeatList[4];
	for(Int i=0; i<4; ++i)
	    repeatList[i] = A.Get(i,i);
	for(Int i=0; i<m; ++i)
	{
	    if( i%5 < 4 )
	        A.Set(i,i,repeatList[i%5]);
	}    
    }
    if( testCorrectness )
    {
        AOrig = A;
	GetDiagonal( A, w );
    }
    if( print )
        Print( A, "A" );

    if( g.Rank() == 0 )
        Output("  Starting triangular eigensolver...");
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    TriangEig( A, X );
    mpi::Barrier( g.Comm() );
    const double runTime = mpi::Time() - startTime;
    if( g.Rank() == 0 )
        Output("  Time = ",runTime," seconds");
    if( print )
    {
        Print( w, "eigenvalues:" );
	Print( X, "eigenvectors:" );
    }
    if( testCorrectness )
        TestCorrectness( print, AOrig, X );
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
        const bool testCorrectness = Input
            ("--correctness","test correctness?",true);
        const bool print = Input("--print","print matrices?",false);
        const bool testReal = Input("--testReal","test real matrices?",true);
        const bool testCpx = Input("--testCpx","test complex matrices?",true);
        const bool repeated = Input("--repeated","test matrices with repeated eigenvalues?",false);
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
              n, g, repeated );
        if( testCpx )
            TestTriangEig<Complex<double>>
            ( testCorrectness, print,
              n, g, repeated );

        // Also test with non-standard distributions
        if( commRank == 0 )
            Output("Nonstandard distributions:");
        if( testReal )
            TestTriangEig<double,MR,MC,MC>
            ( testCorrectness, print,
              n, g, repeated );
        if( testCpx )
            TestTriangEig<Complex<double>,MR,MC,MC>
            ( testCorrectness, print,
              n, g, repeated );
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}

#endif
