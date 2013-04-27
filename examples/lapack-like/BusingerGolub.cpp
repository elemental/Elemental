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
#include "elemental/lapack-like/ApplyColumnPivots.hpp"
#include "elemental/lapack-like/ApplyPackedReflectors.hpp"
#include "elemental/lapack-like/QR.hpp"
#include "elemental/lapack-like/Norm/Frobenius.hpp"
#include "elemental/matrices/Identity.hpp"
#include "elemental/matrices/Uniform.hpp"
using namespace std;
using namespace elem;

typedef double Real;
typedef Complex<Real> C;

int
main( int argc, char* argv[] )
{
    Initialize( argc, argv );

    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::CommRank( comm );

    try 
    {
        const int m = Input("--height","height of matrix",100);
        const int n = Input("--width","width of matrix",100);
        const bool alwaysRecompute = Input("--always","no norm updates?",false);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        DistMatrix<C> A;
        Uniform( A, m, n );
        const Real frobA = FrobeniusNorm( A );
        if( print )
            A.Print("A");

        // Compute the QR decomposition of A, but do not overwrite A
        DistMatrix<C> QRFact( A );
        DistMatrix<C,MD,STAR> t;
        DistMatrix<int,VR,STAR> p;
        qr::BusingerGolub( QRFact, t, p, alwaysRecompute );
        if( print )
        {
            QRFact.Print("QR");
            t.Print("t");
            p.Print("p");
        }

        // Check the error in the QR factorization, 
        // || A P - Q R ||_F / || A ||_F
        DistMatrix<C> E( QRFact );
        MakeTriangular( UPPER, E );
        ApplyPackedReflectors
        ( LEFT, LOWER, VERTICAL, BACKWARD, UNCONJUGATED, 0, QRFact, t, E );
        ApplyColumnPivots( A, p ); 
        Axpy( C(-1), A, E );
        const Real frobQR = FrobeniusNorm( E );
        if( print )
            E.Print("A P - Q R");

        // Check the numerical orthogonality of Q, || I - Q^H Q ||_F / || A ||_F
        Identity( E, m, n );
        ApplyPackedReflectors
        ( LEFT, LOWER, VERTICAL, BACKWARD, UNCONJUGATED, 0, QRFact, t, E );
        ApplyPackedReflectors
        ( LEFT, LOWER, VERTICAL, FORWARD, CONJUGATED, 0, QRFact, t, E );
        const int k = std::min(m,n);
        DistMatrix<C> EUpper;
        View( EUpper, E, 0, 0, k, k );
        DistMatrix<C> I;
        Identity( I, k, k );
        Axpy( C(-1), I, EUpper );
        const Real frobOrthog = FrobeniusNorm( EUpper ); 
        if( print )
            E.Print("I - Q^H Q");

        if( commRank == 0 )
        {
            std::cout << "|| A ||_F = " << frobA << "\n"
                      << "|| A P - Q R ||_F / || A ||_F   = " 
                      << frobQR/frobA << "\n"
                      << "|| I - Q^H Q ||_F / || A ||_F = "
                      << frobOrthog/frobA << "\n"
                      << std::endl;
        }
    }
    catch( ArgException& e )
    {
        // There is nothing to do
    }
    catch( exception& e )
    {
        ostringstream os;
        os << "Process " << commRank << " caught exception with message: "
           << e.what() << endl;
        cerr << os.str();
#ifndef RELEASE
        DumpCallStack();
#endif
    }

    Finalize();
    return 0;
}
