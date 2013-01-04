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
        const int n = Input("--size","size of HPD matrix",100);
        const bool upper = Input("--upper","upper storage?",false);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        Grid g( comm );
        DistMatrix<C> A( g );
        HermitianUniformSpectrum( n, A, R(1), R(20) );

        if( print )
            A.Print("A");

        // Make a copy of A and then overwrite it with its inverse
        const UpperOrLower uplo = ( upper ? UPPER : LOWER );
        DistMatrix<C> invA( A );
        HPDInverse( uplo, invA );

        if( print )
        {
            MakeHermitian( uplo, invA );
            invA.Print("inv(A)");
        }

        // Form I - invA*A and print the relevant norms
        DistMatrix<C> E( g );
        Identity( n, n, E );
        Hemm( LEFT, uplo, C(-1), invA, A, C(1), E );

        const R frobNormA = HermitianNorm( uplo, A, FROBENIUS_NORM );
        const R frobNormInvA = HermitianNorm( uplo, invA, FROBENIUS_NORM );
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

