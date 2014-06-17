/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "El.hpp" instead
#include "El-lite.hpp"

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
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        SetBlocksize( nb );

        Grid g( mpi::COMM_WORLD );
        if( mpi::WorldRank() == 0 )
            std::cout << "Grid is " 
                      << g.Height() << " x " << g.Width() << std::endl;
        DistMatrix<C> A(g);
        Uniform( A, m, n );
        if( print )
            Print( A, "A" );

        // Compute just the singular values 
        DistMatrix<Real,VR,STAR> sOnly(g);
        auto U( A );
        SVD( U, sOnly );

        // Compute the SVD of A 
        DistMatrix<C> V(g);
        DistMatrix<Real,VR,STAR> s(g);
        U = A;
        SVD( U, s, V );
        if( print )
        {
            Print( U, "U" );
            Print( V, "V" );
            Print( s, "s" );
        }

        // Compare the singular values from both methods
        Axpy( Real(-1), s, sOnly );
        const Real singValDiff = FrobeniusNorm( sOnly );
        const Real twoNormOfA = MaxNorm( s );
        const Real maxNormOfA = MaxNorm( A );
        const Real oneNormOfA = OneNorm( A );
        const Real infNormOfA = InfinityNorm( A );
        const Real frobNormOfA = FrobeniusNorm( A );

        DiagonalScale( RIGHT, NORMAL, s, U );
        Gemm( NORMAL, ADJOINT, C(-1), U, V, C(1), A );
        const Real maxNormOfE = MaxNorm( A );
        const Real oneNormOfE = OneNorm( A );
        const Real infNormOfE = InfinityNorm( A );
        const Real frobNormOfE = FrobeniusNorm( A );
        const Real epsilon = lapack::MachineEpsilon<Real>();
        const Real scaledResidual = frobNormOfE / (max(m,n)*epsilon*twoNormOfA);

        if( mpi::WorldRank() == 0 )
        {
            cout << "||A||_max   = " << maxNormOfA << "\n"
                 << "||A||_1     = " << oneNormOfA << "\n"
                 << "||A||_oo    = " << infNormOfA << "\n"
                 << "||A||_F     = " << frobNormOfA << "\n"
                 << "||A||_2     = " << twoNormOfA << "\n"
                 << "\n"
                 << "||A - U Sigma V^H||_max = " << maxNormOfE << "\n"
                 << "||A - U Sigma V^H||_1   = " << oneNormOfE << "\n"
                 << "||A - U Sigma V^H||_oo  = " << infNormOfE << "\n"
                 << "||A - U Sigma V^H||_F   = " << frobNormOfE << "\n"
                 << "||A - U Sigma V_H||_F / (max(m,n) eps ||A||_2) = " 
                 << scaledResidual << "\n" 
                 << "\n"
                 << "|| sError ||_2 = " << singValDiff << std::endl;
        }
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
