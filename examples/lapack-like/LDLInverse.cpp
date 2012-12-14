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

        // Make a copy of A and then overwrite it with its inverse
        // WARNING: There is no pivoting here!
        DistMatrix<C> invA( A );
        if( conjugate )
            LDLH( invA );
        else
            LDLT( invA );
        TriangularInverse( LOWER, UNIT, invA );
        Trdtrmm( orientation, LOWER, invA );

        // Form I - invA*A and print the relevant norms
        DistMatrix<C> E( g );
        Identity( n, n, E );
        if( conjugate )
            Hemm( LEFT, LOWER, C(-1), invA, A, C(1), E );
        else
            Symm( LEFT, LOWER, C(-1), invA, A, C(1), E );

        const R frobNormA = Norm( A, FROBENIUS_NORM );
        const R frobNormInvA = 
            ( conjugate ? HermitianNorm( LOWER, invA, FROBENIUS_NORM )
                        : SymmetricNorm( LOWER, invA, FROBENIUS_NORM ) );
        const R frobNormError = Norm( E, FROBENIUS_NORM );
        if( g.Rank() == 0 )
        {
            std::cout << "|| A          ||_F = " << frobNormA << "\n"
                      << "|| invA       ||_F = " << frobNormInvA << "\n"
                      << "|| I - invA A ||_F = " << frobNormError << "\n"
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

