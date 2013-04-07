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

        // Create a 2d process grid from a communicator. In our case, it is
        // MPI_COMM_WORLD. There is another constructor that allows you to 
        // specify the grid dimensions, Grid g( comm, r, c ), which creates an 
        // r x c grid.
        Grid grid( comm );
    
        // Create an n x n complex matrix residing on a single process.
        Matrix<C> HRoot;
        if( commRank == 0 )
        {
            HRoot.ResizeTo( n, n, n ); // n x n, with a leading dimension of n
            // Set entry (i,j) to (i+j,i-j)
            for( int j=0; j<n; ++j )
                for( int i=0; i<n; ++i )
                    HRoot.Set( i, j, C(i+j,i-j) );
            HRoot.Print("H on process 0");
        }
        mpi::Barrier( comm );

        // Broadcast the local matrix to all processes so that they all have
        // a copy of H, say H[* ,* ]
        DistMatrix<C,STAR,STAR> H_STAR_STAR( n, n, n, grid );
        if( commRank == 0 ) 
        {
            mpi::Broadcast( HRoot.Buffer(), n*n, 0, comm );
            MemCopy( H_STAR_STAR.Buffer(), HRoot.Buffer(), n*n );
        }
        else
            mpi::Broadcast( H_STAR_STAR.Buffer(), n*n, 0, comm );
        if( print )
            H_STAR_STAR.Print("H[* ,* ]");

        // Now that we have a valid DistMatrix (in a [* ,* ] distribution), 
        // we can trivially redistribute into the usual matrix distribution
        DistMatrix<C> H( H_STAR_STAR );
        if( print )
            H.Print("H");

        // Call the eigensolver. We first create an empty complex eigenvector 
        // matrix, X, and an eigenvalue column vector, w[VR,* ]
        DistMatrix<R,VR,STAR> w_VR_STAR( grid );
        DistMatrix<C> X( grid );
        // Optional: set blocksizes and algorithmic choices here. See the 
        //           'Tuning' section of the README for details.
        HermitianEig( LOWER, H, w_VR_STAR, X );

        // Sort the eigensolution
        hermitian_eig::Sort( w_VR_STAR, X );
        if( print )
        {
            w_VR_STAR.Print("Eigenvalues of H");
            X.Print("Eigenvectors of H");
        }

        // Store a complete copy of w and X on the root
        Matrix<R> wLocal;
        Matrix<C> XLocal;
        {
            // Give every process a full copy of w and X
            DistMatrix<R,STAR,STAR> w_STAR_STAR( w_VR_STAR );
            DistMatrix<C,STAR,STAR> X_STAR_STAR( X );

            // Copy the data into a sequential matrix if we are the root matrix
            if( commRank == 0 )
            {
                wLocal = w_STAR_STAR.Matrix();
                XLocal = X_STAR_STAR.Matrix();
                if( print )
                {
                    wLocal.Print("Eigenvalues on root process");
                    XLocal.Print("Eigenvectors on root process");
                }
            }
        }
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

