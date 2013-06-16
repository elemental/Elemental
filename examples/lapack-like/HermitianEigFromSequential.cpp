/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include "elemental/lapack-like/HermitianEig/Sort.hpp"
using namespace std;
using namespace elem;

// Typedef our real and complex types to 'R' and 'C' for convenience
typedef double R;
typedef Complex<R> C;

int
main( int argc, char* argv[] )
{
    // This detects whether or not you have already initialized MPI and 
    // does so if necessary. The full routine is elem::Initialize.
    Initialize( argc, argv );

    // Extract our MPI rank
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::CommRank( comm );

    // Surround the Elemental calls with try/catch statements in order to 
    // safely handle any exceptions that were thrown during execution.
    try 
    {
        const int n = Input("--size","size of matrix",100);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        // Create an n x n complex matrix residing on a single process.
        DistMatrix<C,CIRC,CIRC> HRoot( n, n );
        if( commRank == 0 )
        {
            // Set entry (i,j) to (i+j,i-j)
            for( int j=0; j<n; ++j )
                for( int i=0; i<n; ++i )
                    HRoot.SetLocal( i, j, C(i+j,i-j) );
        }
        if( print )
            Print( HRoot, "H on process 0" );

        // Redistribute into the usual matrix distribution
        DistMatrix<C> H( HRoot );
        if( print )
            Print( H, "H" );

        // Call the eigensolver. We first create an empty complex eigenvector 
        // matrix, X, and an eigenvalue column vector, w[VR,* ]
        DistMatrix<R,VR,STAR> w_VR_STAR;
        DistMatrix<C> X;
        // Optional: set blocksizes and algorithmic choices here. See the 
        //           'Tuning' section of the README for details.
        HermitianEig( LOWER, H, w_VR_STAR, X );

        // Sort the eigensolution
        hermitian_eig::Sort( w_VR_STAR, X );
        if( print )
        {
            Print( w_VR_STAR, "Eigenvalues of H" );
            Print( X, "Eigenvectors of H" );
        }

        // Store a complete copy of w and X on the root
        DistMatrix<R,CIRC,CIRC> wRoot( w_VR_STAR );
        DistMatrix<C,CIRC,CIRC> XRoot( X );
        if( print )
        {
            Print( wRoot, "Eigenvalues on root process" );
            Print( XRoot, "Eigenvectors on root process" );
        }
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}

