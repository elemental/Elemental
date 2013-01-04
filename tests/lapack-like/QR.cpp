/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <ctime>
#include "elemental.hpp"
using namespace std;
using namespace elem;

template<typename R> 
void TestCorrectness
( bool print,
  const DistMatrix<R>& A,
        DistMatrix<R>& AOrig )
{
    const Grid& g = A.Grid();
    const int m = A.Height();
    const int n = A.Width();
    const int minDim = std::min(m,n);

    if( g.Rank() == 0 )
        cout << "  Testing orthogonality of Q..." << endl;

    // Form Z := Q^H Q as an approximation to identity
    DistMatrix<R> Z(g);
    Identity( m, n, Z );
    ApplyPackedReflectors( LEFT, LOWER, VERTICAL, BACKWARD, 0, A, Z );
    ApplyPackedReflectors( LEFT, LOWER, VERTICAL, FORWARD, 0, A, Z );

    DistMatrix<R> ZUpper(g);
    View( ZUpper, Z, 0, 0, minDim, minDim );

    // Form Identity
    DistMatrix<R> X(g);
    Identity( minDim, minDim, X );

    // Form X := I - Q^H Q
    Axpy( R(-1), ZUpper, X );

    R oneNormOfError = Norm( X, ONE_NORM );
    R infNormOfError = Norm( X, INFINITY_NORM );
    R frobNormOfError = Norm( X, FROBENIUS_NORM );
    if( g.Rank() == 0 )
    {
        cout << "    ||Q^H Q - I||_1  = " << oneNormOfError << "\n"
             << "    ||Q^H Q - I||_oo = " << infNormOfError << "\n"
             << "    ||Q^H Q - I||_F  = " << frobNormOfError << endl;
    }

    if( g.Rank() == 0 )
        cout << "  Testing if A = QR..." << endl;

    // Form Q R
    DistMatrix<R> U( A );
    MakeTriangular( UPPER, U );
    ApplyPackedReflectors( LEFT, LOWER, VERTICAL, BACKWARD, 0, A, U );

    // Form Q R - A
    Axpy( R(-1), AOrig, U );
    
    const R oneNormOfA = Norm( AOrig, ONE_NORM );
    const R infNormOfA = Norm( AOrig, INFINITY_NORM );
    const R frobNormOfA = Norm( AOrig, FROBENIUS_NORM );
    oneNormOfError = Norm( U, ONE_NORM );
    infNormOfError = Norm( U, INFINITY_NORM );
    frobNormOfError = Norm( U, FROBENIUS_NORM );
    if( g.Rank() == 0 )
    {
        cout << "    ||A||_1       = " << oneNormOfA << "\n"
             << "    ||A||_oo      = " << infNormOfA << "\n"
             << "    ||A||_F       = " << frobNormOfA << "\n"
             << "    ||A - QR||_1  = " << oneNormOfError << "\n"
             << "    ||A - QR||_oo = " << infNormOfError << "\n"
             << "    ||A - QR||_F  = " << frobNormOfError << endl;
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

    // Form Z := Q^H Q as an approximation to identity
    DistMatrix<C> Z(g);
    Identity( m, n, Z );
    ApplyPackedReflectors
    ( LEFT, LOWER, VERTICAL, BACKWARD, UNCONJUGATED, 0, A, t, Z );
    ApplyPackedReflectors
    ( LEFT, LOWER, VERTICAL, FORWARD, CONJUGATED, 0, A, t, Z );
    
    DistMatrix<C> ZUpper(g);
    View( ZUpper, Z, 0, 0, minDim, minDim );

    // Form Identity
    DistMatrix<C> X(g);
    Identity( minDim, minDim, X );

    // Form X := I - Q^H Q
    Axpy( C(-1), ZUpper, X );

    R oneNormOfError = Norm( X, ONE_NORM );
    R infNormOfError = Norm( X, INFINITY_NORM );
    R frobNormOfError = Norm( X, FROBENIUS_NORM );
    if( g.Rank() == 0 )
    {
        cout << "    ||Q^H Q - I||_1  = " << oneNormOfError << "\n"
             << "    ||Q^H Q - I||_oo = " << infNormOfError << "\n"
             << "    ||Q^H Q - I||_F  = " << frobNormOfError << endl;
    }

    if( g.Rank() == 0 )
        cout << "  Testing if A = QR..." << endl;

    // Form Q R
    DistMatrix<C> U( A );
    MakeTriangular( UPPER, U );
    ApplyPackedReflectors
    ( LEFT, LOWER, VERTICAL, BACKWARD, UNCONJUGATED, 0, A, t, U );

    // Form Q R - A
    Axpy( C(-1), AOrig, U );
    
    const R oneNormOfA = Norm( AOrig, ONE_NORM );
    const R infNormOfA = Norm( AOrig, INFINITY_NORM );
    const R frobNormOfA = Norm( AOrig, FROBENIUS_NORM );
    oneNormOfError = Norm( U, ONE_NORM );
    infNormOfError = Norm( U, INFINITY_NORM );
    frobNormOfError = Norm( U, FROBENIUS_NORM );
    if( g.Rank() == 0 )
    {
        cout << "    ||A||_1       = " << oneNormOfA << "\n"
             << "    ||A||_oo      = " << infNormOfA << "\n"
             << "    ||A||_F       = " << frobNormOfA << "\n"
             << "    ||A - QR||_1  = " << oneNormOfError << "\n"
             << "    ||A - QR||_oo = " << infNormOfError << "\n"
             << "    ||A - QR||_F  = " << frobNormOfError << endl;
    }
}

template<typename R>
void TestRealQR
( bool testCorrectness, bool print,
  int m, int n, const Grid& g )
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
        cout << "  Starting QR factorization...";
        cout.flush();
    }
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    QR( A );
    mpi::Barrier( g.Comm() );
    const double runTime = mpi::Time() - startTime;
    const double mD = double(m);
    const double nD = double(n);
    const double gFlops = (2.*mD*nD*nD - 2./3.*nD*nD*nD)/(1.e9*runTime);
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
void TestComplexQR
( bool testCorrectness, bool print,
  int m, int n, const Grid& g )
{
    typedef Complex<R> C;
    DistMatrix<C> A(g), AOrig(g);
    DistMatrix<C,MD,STAR> t(g);

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
        cout << "  Starting QR factorization...";
        cout.flush();
    }
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    QR( A, t );
    mpi::Barrier( g.Comm() );
    const double runTime = mpi::Time() - startTime;
    const double mD = double(m);
    const double nD = double(n);
    const double gFlops = (8.*mD*nD*nD - 8./3.*nD*nD*nD)/(1.e9*runTime);
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
            cout << "Will test QR" << endl;

        if( commRank == 0 )
        {
            cout << "---------------------\n"
                 << "Testing with doubles:\n"
                 << "---------------------" << endl;
        }
        TestRealQR<double>( testCorrectness, print, m, n, g );

        if( commRank == 0 )
        {
            cout << "--------------------------------------\n"
                 << "Testing with double-precision complex:\n"
                 << "--------------------------------------" << endl;
        }
        TestComplexQR<double>( testCorrectness, print, m, n, g );
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
