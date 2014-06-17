/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "El.hpp" instead
#include "El-lite.hpp"
#include EL_IDENTITY_INC

using namespace std;
using namespace El;

// Typedef our real and complex types to 'Real' and 'C' for convenience
typedef double Real;
typedef Complex<Real> C;

int
main( int argc, char* argv[] )
{
    Initialize( argc, argv );

    try 
    {
        const Int m = Input("--height","height of matrix",100);
        const Int n = Input("--width","width of matrix",100);
        const bool colPiv = Input("--colPiv","QR with col pivoting?",false);
        ProcessInput();
        PrintInputReport();

        DistMatrix<C> A, Q, P;
        Uniform( A, m, n );
        const Real frobA = FrobeniusNorm( A );

        // Compute the polar decomp of A using a QR-based Dynamically Weighted
        // Halley (QDWH) iteration
        Q = A;
        PolarCtrl ctrl;
        ctrl.qdwh = true;
        ctrl.colPiv = colPiv;
        Polar( Q, ctrl );
        Zeros( P, n, n );
        Gemm( ADJOINT, NORMAL, C(1), Q, A, C(0), P );

        // Check and report overall and orthogonality error
        DistMatrix<C> B( A );
        Gemm( NORMAL, NORMAL, C(-1), Q, P, C(1), B );
        const Real frobQDWH = FrobeniusNorm( B );
        Identity( B, n, n );
        Herk( LOWER, ADJOINT, C(1), Q, C(-1), B );
        const Real frobQDWHOrthog = HermitianFrobeniusNorm( LOWER, B );
        if( mpi::WorldRank() == 0 )
        {
            std::cout << ctrl.numIts << " iterations of QDWH\n"
                      << "||A - QP||_F / ||A||_F = " 
                      << frobQDWH/frobA << "\n"
                      << "||I - QQ^H||_F / ||A||_F = " 
                      << frobQDWHOrthog/frobA << "\n"
                      << std::endl;
        }
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
