/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include "elemental/blas-like/level1/DiagonalScale.hpp"
#include "elemental/blas-like/level1/Scale.hpp"
#include "elemental/blas-like/level3/Gemm.hpp"
#include "elemental/blas-like/level3/Herk.hpp"
#include "elemental/lapack-like/SkewHermitianEig.hpp"
#include "elemental/lapack-like/Norm/Frobenius.hpp"
#include "elemental/lapack-like/HermitianEig/Sort.hpp"
#include "elemental/matrices/Identity.hpp"
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
        Grid g( comm );
    
        // Create an n x n complex distributed matrix, 
        // We distribute the matrix using grid 'g'.
        //
        // There are quite a few available constructors, including ones that 
        // allow you to pass in your own local buffer and to specify the 
        // distribution alignments (i.e., which process row and column owns the
        // top-left element)
        DistMatrix<C> S( n, n, g );

        // Fill the matrix since we did not pass in a buffer. 
        //
        // We will fill entry (i,j) with the complex value (i-j,i+j) so that 
        // the global matrix is skew-Hermitian. However, only one triangle of 
        // the matrix actually needs to be filled, the symmetry can be implicit.
        //
        const int colShift = S.ColShift(); // first row we own
        const int rowShift = S.RowShift(); // first col we own
        const int colStride = S.ColStride();
        const int rowStride = S.RowStride();
        const int localHeight = S.LocalHeight();
        const int localWidth = S.LocalWidth();
        for( int jLocal=0; jLocal<localWidth; ++jLocal )
        {
            for( int iLocal=0; iLocal<localHeight; ++iLocal )
            {
                // Our process owns the rows colShift:colStride:n,
                //           and the columns rowShift:rowStride:n
                const int i = colShift + iLocal*colStride;
                const int j = rowShift + jLocal*rowStride;
                S.SetLocal( iLocal, jLocal, C(i-j,i+j) );
            }
        }

        // Make a backup of S before we overwrite it within the eigensolver
        DistMatrix<C> SCopy( S );

        // Call the eigensolver. We first create an empty complex eigenvector 
        // matrix, X[MC,MR], and an eigenvalue column vector, w[VR,* ]
        //
        // Optional: set blocksizes and algorithmic choices here. See the 
        //           'Tuning' section of the README for details.
        DistMatrix<R,VR,STAR> wImag( g );
        DistMatrix<C> X( g );
        SkewHermitianEig( LOWER, S, wImag, X ); // only use lower half of S

        // Optional: sort the eigenpairs
        hermitian_eig::Sort( wImag, X );

        if( print )
        {
            SCopy.Print("S");
            X.Print("X");
            wImag.Print("wImag");
        }

        // Check the residual, || S X - Omega X ||_F
        const R frobS = HermitianFrobeniusNorm( LOWER, SCopy );
        DistMatrix<C> E( X );
        Scale( C(0,1), E );
        DiagonalScale( RIGHT, NORMAL, wImag, E );
        Gemm( NORMAL, NORMAL, C(-1), SCopy, X, C(1), E );
        const R frobResid = FrobeniusNorm( E );

        // Check the orthogonality of X
        Identity( n, n, E );
        Herk( LOWER, NORMAL, C(-1), X, C(1), E );
        const R frobOrthog = HermitianFrobeniusNorm( LOWER, E );

        if( g.Rank() == 0 )
        {
            std::cout << "|| H ||_F = " << frobS << "\n"
                      << "|| H X - X Omega ||_F / || A ||_F = " 
                      << frobResid / frobS << "\n"
                      << "|| X X^H - I ||_F = " << frobOrthog / frobS
                      << "\n" << std::endl;
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
