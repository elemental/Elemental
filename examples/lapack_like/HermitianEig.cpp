/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace std;
using namespace El;

// Typedef our real and complex types to 'Real' and 'C' for convenience
typedef double Real;
typedef Complex<Real> C;

int
main( int argc, char* argv[] )
{
    // This detects whether or not you have already initialized MPI and 
    // does so if necessary. 
    Environment env( argc, argv );

    // Surround the Elemental calls with try/catch statements in order to 
    // safely handle any exceptions that were thrown during execution.
    try 
    {
        const Int n = Input("--size","size of matrix",100);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        // Create a 2d process grid from a communicator. In our case, it is
        // mpi::COMM_WORLD. There is another constructor that allows you to 
        // specify the grid dimensions, Grid g( comm, r ), which creates a
        // grid of height r.
        Grid g( mpi::COMM_WORLD );
    
        // Create an n x n complex distributed matrix, 
        // We distribute the matrix using grid 'g'.
        //
        // There are quite a few available constructors, including ones that 
        // allow you to pass in your own local buffer and to specify the 
        // distribution alignments (i.e., which process row and column owns the
        // top-left element)
        DistMatrix<C> H( n, n, g );

        // Fill the matrix since we did not pass in a buffer. 
        //
        // We will fill entry (i,j) with the complex value (i+j,i-j) so that 
        // the global matrix is Hermitian. However, only one triangle of the 
        // matrix actually needs to be filled, the symmetry can be implicit.
        //
        const Int localHeight = H.LocalHeight();
        const Int localWidth = H.LocalWidth();
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            // Our process owns the rows colShift:colStride:n,
            //           and the columns rowShift:rowStride:n
            const Int j = H.GlobalCol(jLoc);
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            {
                const Int i = H.GlobalRow(iLoc);
                H.SetLocal( iLoc, jLoc, C(i+j,i-j) );
            }
        }

        // Make a backup of H before we overwrite it within the eigensolver
        auto HCopy( H );

        // Call the eigensolver. We first create an empty complex eigenvector 
        // matrix, X[MC,MR], and an eigenvalue column vector, w[VR,* ]
        //
        // Optional: set blocksizes and algorithmic choices here. See the 
        //           'Tuning' section of the README for details.
        DistMatrix<Real,VR,STAR> w( g );
        DistMatrix<C> X( g );
        HermitianEig( LOWER, H, w, X, ASCENDING ); 

        if( print )
        {
            Print( HCopy, "H" );
            Print( X, "X" );
            Print( w, "w" );
        }

        // Check the residual, || H X - Omega X ||_F
        const Real frobH = HermitianFrobeniusNorm( LOWER, HCopy );
        auto E( X );
        DiagonalScale( RIGHT, NORMAL, w, E );
        Hemm( LEFT, LOWER, C(-1), HCopy, X, C(1), E );
        const Real frobResid = FrobeniusNorm( E );

        // Check the orthogonality of X
        Identity( E, n, n );
        Herk( LOWER, NORMAL, Real(-1), X, Real(1), E );
        const Real frobOrthog = HermitianFrobeniusNorm( LOWER, E );

        if( mpi::Rank() == 0 )
            Output
            ("|| H ||_F = ",frobH,"\n",
             "|| H X - X Omega ||_F / || A ||_F = ",frobResid/frobH,"\n",
             "|| X X^H - I ||_F = ",frobOrthog,"\n");
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
