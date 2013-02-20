/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include "elemental/blas-like/level1/MakeTriangular.hpp"
#include "elemental/lapack-like/ApplyPackedReflectors.hpp"
#include "elemental/lapack-like/LQ.hpp"
#include "elemental/lapack-like/Norm/Frobenius.hpp"
#include "elemental/lapack-like/Norm/Infinity.hpp"
#include "elemental/lapack-like/Norm/One.hpp"
#include "elemental/matrices/Identity.hpp"
#include "elemental/matrices/Uniform.hpp"
using namespace std;
using namespace elem;

template<typename R> 
void TestCorrectness
( bool print, const DistMatrix<R>& A, DistMatrix<R>& AOrig )
{
    const Grid& g = A.Grid();
    const int m = A.Height();
    const int n = A.Width();
    const int minDim = std::min(m,n);

    if( g.Rank() == 0 )
        cout << "  Testing orthogonality of Q..." << endl;

    // Form Z := Q Q^H as an approximation to identity
    DistMatrix<R> Z(m,n,g);
    MakeIdentity( Z );
    ApplyPackedReflectors
    ( RIGHT, UPPER, HORIZONTAL, BACKWARD, 0, A, Z );
    ApplyPackedReflectors
    ( RIGHT, UPPER, HORIZONTAL, FORWARD, 0, A, Z );

    DistMatrix<R> ZUpper(g);
    View( ZUpper, Z, 0, 0, minDim, minDim );

    // Form Identity
    DistMatrix<R> X(minDim,minDim,g);
    MakeIdentity( X );

    // Form X := I - Q Q^H
    Axpy( R(-1), ZUpper, X );

    R oneNormOfError = OneNorm( X );
    R infNormOfError = InfinityNorm( X );
    R frobNormOfError = FrobeniusNorm( X );
    if( g.Rank() == 0 )
    {
        cout << "    ||Q Q^H - I||_1  = " << oneNormOfError << "\n"
             << "    ||Q Q^H - I||_oo = " << infNormOfError << "\n"
             << "    ||Q Q^H - I||_F  = " << frobNormOfError << endl;
    }

    if( g.Rank() == 0 )
        cout << "  Testing if A = LQ..." << endl;

    // Form L Q
    DistMatrix<R> L( A );
    MakeTriangular( LOWER, L );
    ApplyPackedReflectors
    ( RIGHT, UPPER, HORIZONTAL, BACKWARD, 0, A, L );

    // Form L Q - A
    Axpy( R(-1), AOrig, L );
    
    const R oneNormOfA = OneNorm( AOrig );
    const R infNormOfA = InfinityNorm( AOrig );
    const R frobNormOfA = FrobeniusNorm( AOrig );
    oneNormOfError = OneNorm( L );
    infNormOfError = InfinityNorm( L );
    frobNormOfError = FrobeniusNorm( L );
    if( g.Rank() == 0 )
    {
        cout << "    ||A||_1       = " << oneNormOfA << "\n"
             << "    ||A||_oo      = " << infNormOfA << "\n"
             << "    ||A||_F       = " << frobNormOfA << "\n"
             << "    ||A - LQ||_1  = " << oneNormOfError << "\n"
             << "    ||A - LQ||_oo = " << infNormOfError << "\n"
             << "    ||A - LQ||_F  = " << frobNormOfError << endl;
    }
}

template<typename R> 
void TestCorrectness
( bool print,
  const DistMatrix<Complex<R> >& A,
  const DistMatrix<Complex<R>,MD,STAR>& t,
        DistMatrix<Complex<R> >& AOrig )
{
    typedef Complex<R> C;

    const Grid& g = A.Grid();
    const int m = A.Height();
    const int n = A.Width();
    const int minDim = std::min(m,n);

    if( g.Rank() == 0 )
        cout << "  Testing orthogonality of Q..." << endl;

    // Form Z := Q Q^H as an approximation to identity
    DistMatrix<C> Z(m,n,g);
    MakeIdentity( Z );
    ApplyPackedReflectors
    ( RIGHT, UPPER, HORIZONTAL, BACKWARD, UNCONJUGATED, 0, A, t, Z );
    ApplyPackedReflectors
    ( RIGHT, UPPER, HORIZONTAL, FORWARD, CONJUGATED, 0, A, t, Z );
    
    DistMatrix<C> ZUpper(g);
    View( ZUpper, Z, 0, 0, minDim, minDim );

    // Form Identity
    DistMatrix<C> X(minDim,minDim,g);
    MakeIdentity( X );

    // Form X := I - Q Q^H
    Axpy( C(-1), ZUpper, X );

    R oneNormOfError = OneNorm( X );
    R infNormOfError = InfinityNorm( X );
    R frobNormOfError = FrobeniusNorm( X );
    if( g.Rank() == 0 )
    {
        cout << "    ||Q Q^H - I||_1  = " << oneNormOfError << "\n"
             << "    ||Q Q^H - I||_oo = " << infNormOfError << "\n"
             << "    ||Q Q^H - I||_F  = " << frobNormOfError << endl;
    }

    if( g.Rank() == 0 )
        cout << "  Testing if A = LQ..." << endl;

    // Form L Q
    DistMatrix<C> L( A );
    MakeTriangular( LOWER, L );
    ApplyPackedReflectors
    ( RIGHT, UPPER, HORIZONTAL, BACKWARD, UNCONJUGATED, 0, A, t, L );

    // Form L Q - A
    Axpy( C(-1), AOrig, L );
    
    const R oneNormOfA = OneNorm( AOrig );
    const R infNormOfA = InfinityNorm( AOrig );
    const R frobNormOfA = FrobeniusNorm( AOrig );
    oneNormOfError = OneNorm( L );
    infNormOfError = InfinityNorm( L );
    frobNormOfError = FrobeniusNorm( L );
    if( g.Rank() == 0 )
    {
        cout << "    ||A||_1       = " << oneNormOfA << "\n"
             << "    ||A||_oo      = " << infNormOfA << "\n"
             << "    ||A||_F       = " << frobNormOfA << "\n"
             << "    ||A - LQ||_1  = " << oneNormOfError << "\n"
             << "    ||A - LQ||_oo = " << infNormOfError << "\n"
             << "    ||A - LQ||_F  = " << frobNormOfError << endl;
    }
}

template<typename R>
void TestRealLQ
( bool testCorrectness, bool print, int m, int n, const Grid& g )
{
    DistMatrix<R> A(g), AOrig(g);
    Uniform( m, n, A );

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
        A.Print("A");

    if( g.Rank() == 0 )
    {
        cout << "  Starting LQ factorization...";
        cout.flush();
    }
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    LQ( A );
    mpi::Barrier( g.Comm() );
    const double runTime = mpi::Time() - startTime;
    const double mD = double(m);
    const double nD = double(n);
    const double gFlops = (2.*mD*mD*nD - 2./3.*mD*mD*mD)/(1.e9*runTime);
    if( g.Rank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
    if( print )
        A.Print("A after factorization");
    if( testCorrectness )
        TestCorrectness( print, A, AOrig );
}

template<typename R>
void TestComplexLQ
( bool testCorrectness, bool print, int m, int n, const Grid& g )
{
    typedef Complex<R> C;
    DistMatrix<C> A(g), AOrig(g);
    Uniform( m, n, A );

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
        A.Print("A");
    DistMatrix<C,MD,STAR> t(g);

    if( g.Rank() == 0 )
    {
        cout << "  Starting LQ factorization...";
        cout.flush();
    }
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    LQ( A, t );
    mpi::Barrier( g.Comm() );
    const double runTime = mpi::Time() - startTime;
    const double mD = double(m);
    const double nD = double(n);
    const double gFlops = (8.*mD*mD*nD - 8./3.*mD*mD*mD)/(1.e9*runTime);
    if( g.Rank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
    if( print )
        A.Print("A after factorization");
    if( testCorrectness )
        TestCorrectness( print, A, t, AOrig );
}

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::CommRank( comm );
    const int commSize = mpi::CommSize( comm );

    try
    {
        int r = Input("--gridHeight","height of process grid",0);
        const int m = Input("--height","height of matrix",100);
        const int n = Input("--width","width of matrix",100);
        const int nb = Input("--nb","algorithmic blocksize",96);
        const bool testCorrectness = Input
            ("--correctness","test correctness?",true);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        if( r == 0 )
            r = Grid::FindFactor( commSize );
        const int c = commSize / r;
        const Grid g( comm, r, c );
        SetBlocksize( nb );
#ifndef RELEASE
        if( commRank == 0 )
        {
            cout << "==========================================\n"
                 << " In debug mode! Performance will be poor! \n"
                 << "==========================================" << endl;
        }
#endif
        if( commRank == 0 )
            cout << "Will test LQ" << endl;

        if( commRank == 0 )
        {
            cout << "---------------------\n"
                 << "Testing with doubles:\n"
                 << "---------------------" << endl;
        }
        TestRealLQ<double>( testCorrectness, print, m, n, g );

        if( commRank == 0 )
        {
            cout << "--------------------------------------\n"
                 << "Testing with double-precision complex:\n"
                 << "--------------------------------------" << endl;
        }
        TestComplexLQ<double>( testCorrectness, print, m, n, g );
    }
    catch( ArgException& e ) { }
    catch( exception& e )
    {
        ostringstream os;
        os << "Process " << commRank << " caught error message:\n" << e.what()
           << endl;
        cerr << os.str();
#ifndef RELEASE
        DumpCallStack();
#endif
    }   
    Finalize();
    return 0;
}
