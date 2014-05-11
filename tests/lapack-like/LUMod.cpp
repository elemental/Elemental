/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include ELEM_GEMM_INC
#include ELEM_GER_INC
#include ELEM_GERU_INC
#include ELEM_LU_INC
#include ELEM_FROBENIUSNORM_INC
#include ELEM_INFINITYNORM_INC
#include ELEM_ONENORM_INC
#include ELEM_UNIFORM_INC
#include ELEM_EXPLICITPERMUTATION_INC
using namespace std;
using namespace elem;

template<typename F,Dist UPerm> 
void TestCorrectness
( bool print, 
  const DistMatrix<F>& A,
  const DistMatrix<Int,UPerm,STAR>& perm,
  const DistMatrix<F>& AOrig )
{
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    const Int n = AOrig.Width();

    if( g.Rank() == 0 )
        cout << "Testing error..." << endl;

    // Generate random right-hand sides
    DistMatrix<F> X(g);
    Uniform( X, n, 100 );
    auto Y( X );
    lu::SolveAfter( NORMAL, A, perm, Y );

    // Now investigate the residual, ||AOrig Y - X||_oo
    const Real oneNormOfX = OneNorm( X );
    const Real infNormOfX = InfinityNorm( X );
    const Real frobNormOfX = FrobeniusNorm( X );
    Gemm( NORMAL, NORMAL, F(-1), AOrig, Y, F(1), X );
    const Real oneNormOfError = OneNorm( X );
    const Real infNormOfError = InfinityNorm( X );
    const Real frobNormOfError = FrobeniusNorm( X );
    const Real oneNormOfA = OneNorm( AOrig );
    const Real infNormOfA = InfinityNorm( AOrig );
    const Real frobNormOfA = FrobeniusNorm( AOrig );

    if( g.Rank() == 0 )
    {
        cout << "||A||_1             = " << oneNormOfA << "\n"
             << "||A||_oo            = " << infNormOfA << "\n"
             << "||A||_F             = " << frobNormOfA << "\n"
             << "||X||_1             = " << oneNormOfX << "\n"
             << "||X||_oo            = " << infNormOfX << "\n"
             << "||X||_F             = " << frobNormOfX << "\n"
             << "||A A^-1 X - X||_1  = " << oneNormOfError << "\n"
             << "||A A^-1 X - X||_oo = " << infNormOfError << "\n"
             << "||A A^-1 X - X||_F  = " << frobNormOfError << endl;
    }
}

template<typename F,Dist UPerm> 
void TestLUMod
( bool conjugate, Base<F> tau,
  bool testCorrectness, bool print, Int m, const Grid& g )
{
    DistMatrix<F> A(g), AOrig(g);
    DistMatrix<Int,UPerm,STAR> perm(g);

    Uniform( A, m, m );
    if( testCorrectness )
    {
        if( g.Rank() == 0 )
        {
            cout << "  Making copy of original matrix...";
            cout.flush();
        }
        AOrig = A;
        if( g.Rank() == 0 )
            cout << "DONE" << endl;
    }
    if( print )
        Print( A, "A" );

    {
        if( g.Rank() == 0 )
        {
            cout << "  Starting full LU factorization...";
            cout.flush();
        }
        mpi::Barrier( g.Comm() );
        const double startTime = mpi::Time();
        LU( A, perm );
        mpi::Barrier( g.Comm() );
        const double runTime = mpi::Time() - startTime;
        const double realGFlops = 2./3.*Pow(double(m),3.)/(1.e9*runTime);
        const double gFlops = ( IsComplex<F>::val ? 4*realGFlops : realGFlops );
        if( g.Rank() == 0 )
            cout << "DONE.\n" << "  Time = " << runTime << " seconds. GFlops = "
                 << gFlops << endl;
    }

    if( print )
    {
        Print( A, "A after original factorization" );
        Print( perm, "perm after original factorization");
        DistMatrix<Int> P(g);
        ExplicitPermutation( perm, P );
        Print( P, "P" );
    }

    // Generate random vectors u and v
    DistMatrix<F> u(g), v(g);
    Uniform( u, m, 1 );
    Uniform( v, m, 1 );
    if( testCorrectness )
    {
        if( conjugate )
            Ger( F(1), u, v, AOrig );
        else
            Geru( F(1), u, v, AOrig );
    }

    { 
        if( g.Rank() == 0 )
        {
            cout << "  Starting rank-one LU modification...";
            cout.flush();
        }
        mpi::Barrier( g.Comm() );
        const double startTime = mpi::Time();
        LUMod( A, perm, u, v, conjugate, tau );
        mpi::Barrier( g.Comm() );
        const double runTime = mpi::Time() - startTime;
        if( g.Rank() == 0 )
            cout << "DONE.\n" << "  Time = " << runTime << " seconds." << endl;
    }

    if( print )
    {
        Print( A, "A after modification" );
        Print( perm, "perm after modification");
        DistMatrix<Int> P(g);
        ExplicitPermutation( perm, P );
        Print( P, "P" );
    }

    if( testCorrectness )
        TestCorrectness( print, A, perm, AOrig );
}

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const Int commRank = mpi::Rank( comm );
    const Int commSize = mpi::Size( comm );

    try
    {
        Int r = Input("--gridHeight","height of process grid",0);
        const bool colMajor = Input("--colMajor","column-major ordering?",true);
        const Int m = Input("--height","height of matrix",100);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const double tau = Input("--tau","pivot threshold",0.1);
        const bool conjugate = Input("--conjugate","conjugate v?",true);
        const bool testCorrectness = Input
            ("--correctness","test correctness?",true);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        if( r == 0 )
            r = Grid::FindFactor( commSize );
        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( comm, r, order );
        SetBlocksize( nb );
        ComplainIfDebug();

        if( commRank == 0 )
            cout << "Testing with doubles:" << endl;
        TestLUMod<double,MC>( conjugate, tau, testCorrectness, print, m, g );

        if( commRank == 0 )
            cout << "Testing with double-precision complex:" << endl;
        TestLUMod<Complex<double>,MC>
        ( conjugate, tau, testCorrectness, print, m, g );
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
