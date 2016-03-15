/*
   Copyright (c) 2009-2016, Jack Poulson and Tim Moon
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

template<typename F>
void FoxLiSchurFactor( ElementalMatrix<F>& A, Int m )
{
    LogicError( "Fox-Li matrix is complex" );
}

template<typename Real>
void FoxLiSchurFactor( ElementalMatrix<Complex<Real>>& A, Int m )
{
    const Grid& g = A.Grid();
    DistMatrix<Complex<Real>,MC,MR> B(g);
    DistMatrix<Complex<Real>,VR,STAR> w(g);
    FoxLi( B, m, Real(-0.179) );
    Schur( B, w );
    MakeTrapezoidal( UPPER, B, 0 );
    A = B;
}

template<typename F>
void TestCorrectness
( bool print,
  const ElementalMatrix<F>& A,
  const ElementalMatrix<F>& X )
{
    typedef Base<F> Real;
    const Grid& g = A.Grid();

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
    if( g.Rank() == 0 )
    {
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
  Int testMatrix )
{
    typedef Base<F> Real;

    DistMatrix<F,U,V> A(g), AOrig(g), X(g);
    DistMatrix<F,S,STAR> w(g);
    if( g.Rank() == 0 )
        Output("Testing with ",TypeName<F>());

    // Generate test matrix
    switch( testMatrix )
    {

    case 0:
        {
            // LU factorization of Gaussian matrix
            DistMatrix<F> B(g);
            Gaussian( B, m, m );
            LU( B );
            Transpose( B, A );
            MakeTrapezoidal( UPPER, A, 0 );
            break;
        }

    case 1:
        {
            // Schur factorization of Fox-Li matrix
            FoxLiSchurFactor( A, m );
            break;
        }
        
    case 2:
        {
            // Schur factorization of Grcar matrix
            DistMatrix<Complex<Real>,VR,STAR> d(g);
            Grcar( A, m );
            Schur( A, d );
            MakeTrapezoidal( UPPER, A, 0 );
            break;
        }
        
    default:
        LogicError("Unknown test matrix");
        break;

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
    const double startTime = mpi::Time();
    TriangEig( A, X );
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
        const Int testMatrix = Input("--testMatrix","test matrix (0=Gaussian,1=Fox-Li,2=Grcar)",0);
        ProcessInput();
        PrintInputReport();

        if( r == 0 )
            r = Grid::FindFactor( commSize );
        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( comm, r, order );
        SetBlocksize( nb );
        ComplainIfDebug();

        // Test with default distributions
        if( commRank == 0 )
            Output("Normal algorithms:");
        if( testReal )
            TestTriangEig<double>
            ( testCorrectness, print,
              n, g, testMatrix );
        if( testCpx )
            TestTriangEig<Complex<double>>
            ( testCorrectness, print,
              n, g, testMatrix );

        // Test with non-standard distributions
        if( commRank == 0 )
            Output("Non-standard algorithms:");
        if( testReal )
            TestTriangEig<double,MR,MC,MC>
            ( testCorrectness, print,
              n, g, testMatrix );
        if( testCpx )
            TestTriangEig<Complex<double>,MR,MC,MC>
            ( testCorrectness, print,
              n, g, testMatrix );
        
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
