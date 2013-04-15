/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include "elemental/lapack-like/Norm/Frobenius.hpp"
#include "elemental/lapack-like/Norm/Infinity.hpp"
#include "elemental/lapack-like/Norm/One.hpp"
#include "elemental/lapack-like/QR.hpp"
#include "elemental/matrices/Identity.hpp"
#include "elemental/matrices/Uniform.hpp"
using namespace std;
using namespace elem;

template<typename F> 
void TestCorrectness
( const DistMatrix<F,VC,  STAR>& Q,
  const DistMatrix<F,STAR,STAR>& R,
        DistMatrix<F,VC,  STAR>& A )
{
    typedef typename Base<F>::type Real;
    const Grid& g = A.Grid();
    const int m = A.Height();
    const int n = A.Width();

    // Form I - Q^H Q
    if( g.Rank() == 0 )
        cout << "  Testing orthogonality of Q..." << endl;
    DistMatrix<F> Z(g);
    Identity( n, n, Z );
    DistMatrix<F> Q_MC_MR( Q );
    Herk( UPPER, ADJOINT, F(-1), Q_MC_MR, F(1), Z );
    Real oneNormOfError = HermitianOneNorm( UPPER, Z );
    Real infNormOfError = HermitianInfinityNorm( UPPER, Z );
    Real frobNormOfError = HermitianFrobeniusNorm( UPPER, Z );
    if( g.Rank() == 0 )
    {
        cout << "    ||Q^H Q - I||_1  = " << oneNormOfError << "\n"
             << "    ||Q^H Q - I||_oo = " << infNormOfError << "\n"
             << "    ||Q^H Q - I||_F  = " << frobNormOfError << endl;
    }

    // Form A - Q R
    if( g.Rank() == 0 )
        cout << "  Testing if A = QR..." << endl;
    const Real oneNormOfA = OneNorm( A );
    const Real infNormOfA = InfinityNorm( A );
    const Real frobNormOfA = FrobeniusNorm( A );
    LocalGemm( NORMAL, NORMAL, F(-1), Q, R, F(1), A );
    oneNormOfError = OneNorm( A );
    infNormOfError = InfinityNorm( A );
    frobNormOfError = FrobeniusNorm( A );
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

template<typename F>
void TestQR
( bool testCorrectness, bool print,
  int m, int n, const Grid& g )
{
    DistMatrix<F,VC,STAR> A(g), Q(g);
    DistMatrix<F,STAR,STAR> R(g);

    Uniform( m, n, A );
    if( print )
        A.Print("A");
    Q = A;

    if( g.Rank() == 0 )
    {
        cout << "  Starting Cholesky QR factorization...";
        cout.flush();
    }
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    qr::Cholesky( Q, R );
    mpi::Barrier( g.Comm() );
    const double runTime = mpi::Time() - startTime;
    const double mD = double(m);
    const double nD = double(n);
    const double gFlops = (2.*mD*nD*nD + 1./3.*nD*nD*nD)/(1.e9*runTime);
    if( g.Rank() == 0 )
    {
        cout << "DONE. " << endl
             << "  Time = " << runTime << " seconds. GFlops = " 
             << gFlops << endl;
    }
    if( print )
    {
        Q.Print("Q");
        R.Print("R");
    }
    if( testCorrectness )
        TestCorrectness( Q, R, A );
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
            cout << "Will test CholeskyQR" << endl;

        if( commRank == 0 )
        {
            cout << "---------------------\n"
                 << "Testing with doubles:\n"
                 << "---------------------" << endl;
        }
        TestQR<double>( testCorrectness, print, m, n, g );

        if( commRank == 0 )
        {
            cout << "--------------------------------------\n"
                 << "Testing with double-precision complex:\n"
                 << "--------------------------------------" << endl;
        }
        TestQR<double>( testCorrectness, print, m, n, g );
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
