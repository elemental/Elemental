/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental.hpp"
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
        ProcessInput();
        PrintInputReport();

        const Grid g( comm );
        DistMatrix<C> A(g);
        Uniform( m, n, A );
        const Real frobA = Norm( A, FROBENIUS_NORM );

        // Compute the QR decomposition of A, but do not overwrite A
        DistMatrix<C> Q( A ), R(g);
        ExplicitQR( Q, R );

        // Check the error in the QR factorization, || A - Q R ||_F / || A ||_F
        DistMatrix<C> E( A );
        Gemm( NORMAL, NORMAL, C(-1), Q, R, C(1), E );
        const Real frobQR = Norm( E, FROBENIUS_NORM );

        // Check the numerical orthogonality of Q, || I - Q^H Q ||_F / || A ||_F
        const int k = std::min(m,n);
        Identity( k, k, E );
        Herk( LOWER, ADJOINT, C(-1), Q, C(1), E );
        const Real frobOrthog = HermitianNorm( LOWER, E, FROBENIUS_NORM ); 

        if( g.Rank() == 0 )
        {
            std::cout << "|| A ||_F = " << frobA << "\n"
                      << "|| A - Q R ||_F / || A ||_F   = " 
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
