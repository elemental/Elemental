/*
   Copyright (c) 2009-2013, Jack Poulson
   Copyright (c) 2011, The University of Texas at Austin
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include "elemental/matrices/Zeros.hpp"
using namespace elem;

int
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const Int commRank = mpi::CommRank( comm );
    const Int commSize = mpi::CommSize( comm );

    try 
    {
        const Int m = 3*commSize;
        const Int n = 2*commSize;

        Grid g( comm );

        for( Int k=0; k<50; ++k )
        {
            if( commRank == 0 )
                std::cout << "Iteration " << k << std::endl;

            DistMatrix<double> A(g);
            Zeros( A, m, n );

            AxpyInterface<double> interface;
            interface.Attach( LOCAL_TO_GLOBAL, A );
            Matrix<double> X( commSize, 1 );
            for( Int j=0; j<X.Width(); ++j )
                for( Int i=0; i<commSize; ++i )
                    X.Set(i,j,commRank+1);
            for( Int i=0; i<5; ++i )
            {
                interface.Axpy( 2, X, 2*commRank, commRank );
                interface.Axpy( 2, X, 2*commRank, commRank+1 );
            }
            interface.Detach();

            Print( A, "A" );

            interface.Attach( GLOBAL_TO_LOCAL, A );
            Matrix<double> Y;
            if( commRank == 0 )
            {
                Zeros( Y, m, n );
                interface.Axpy( 1.0, Y, 0, 0 );
            }
            interface.Detach();

            if( commRank == 0 )
                Print( Y, "Copy of global matrix on root process:" );

            // TODO: Check to ensure that the result is correct
        }
    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
