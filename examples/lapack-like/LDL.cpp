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
        const int n = Input("--size","size of matrix to factor",100);
        const bool conjugate = Input("--conjugate","LDL^H?",false);
        ProcessInput();
        PrintInputReport();

        const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );
        Grid g( comm );
        DistMatrix<C> A( g );
        if( conjugate )
        {
            HermitianUniformSpectrum( n, A, -30, -20 );
        }
        else
        {
            Uniform( n, n, A );
            DistMatrix<C> ATrans( g );
            Transpose( A, ATrans );
            Axpy( C(1), ATrans, A );
        }

        // Make a copy of A and then overwrite it with its LDL factorization
        // WARNING: There is no pivoting here!
        DistMatrix<C> factA( A );
        DistMatrix<C,MC,STAR> d( g );
        if( conjugate )
            LDLH( factA, d );
        else
            LDLT( factA, d );

        DistMatrix<C> L( factA );
        MakeTriangular( LOWER, L );
        internal::SetDiagonalToOne( LEFT, 0, L );

        DistMatrix<C> LD( L );
        DiagonalScale( RIGHT, NORMAL, d, LD );
        Gemm( NORMAL, orientation, C(-1), LD, L, C(1), A );
        const R frobNormOfError = Norm( A, FROBENIUS_NORM );
        if( commRank == 0 )
            std::cout << "|| A - L D L^[T/H] ||_F = " << frobNormOfError << "\n"
                      << std::endl;
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

