/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental.hpp"
using namespace std;
using namespace elem;

// Typedef our real and complex types to 'R' and 'C' for convenience
typedef double R;
typedef Complex<R> C;

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
        const int nb = Input("--nb","algorithmic blocksize",96);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        SetBlocksize( nb );

        Grid g( comm );
        DistMatrix<C> A( g );
        Uniform( m, n, A );

        if( print )
            A.Print("A");

        // Compute just the singular values 
        DistMatrix<R,VR,STAR> sOnly( g );
        DistMatrix<C> U( A );
        SingularValues( U, sOnly );

        // Compute the SVD of A 
        DistMatrix<C> V( g );
        DistMatrix<R,VR,STAR> s( g );
        U = A;
        SVD( U, s, V );

        if( print )
        {
            U.Print("U");
            V.Print("V");
            s.Print("s");
        }

        // Compare the singular values from both methods
        Axpy( R(-1), s, sOnly );
        const R singValDiff = Norm( sOnly, FROBENIUS_NORM );

        const R twoNormOfA = Norm( s, MAX_NORM );
        const R maxNormOfA = Norm( A, MAX_NORM );
        const R oneNormOfA = Norm( A, ONE_NORM );
        const R infNormOfA = Norm( A, INFINITY_NORM );
        const R frobNormOfA = Norm( A, FROBENIUS_NORM );

        DiagonalScale( RIGHT, NORMAL, s, U );
        Gemm( NORMAL, ADJOINT, C(-1), U, V, C(1), A );
        const R maxNormOfE = Norm( A, MAX_NORM );
        const R oneNormOfE = Norm( A, ONE_NORM );
        const R infNormOfE = Norm( A, INFINITY_NORM );
        const R frobNormOfE = Norm( A, FROBENIUS_NORM );
        const R epsilon = lapack::MachineEpsilon<R>();
        const R scaledResidual = frobNormOfE / (max(m,n)*epsilon*twoNormOfA);

        if( commRank == 0 )
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
