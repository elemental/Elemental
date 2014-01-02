/*
   Copyright (c) 2009-2014, Jack Poulson
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
#include "elemental/matrices/Identity.hpp"
using namespace std;
using namespace elem;

// Typedef our real and complex types to 'Real' and 'C' for convenience
typedef double Real;
typedef Complex<Real> C;

int
main( int argc, char* argv[] )
{
    // This detects whether or not you have already initialized MPI and 
    // does so if necessary. The full routine is elem::Initialize.
    Initialize( argc, argv );

    // Surround the Elemental calls with try/catch statements in order to 
    // safely handle any exceptions that were thrown during execution.
    try 
    {
        const Int n = Input("--size","size of matrix",100);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        DistMatrix<C> S( n, n );

        // We will fill entry (i,j) with the complex value (i-j,i+j) so that 
        // the global matrix is skew-Hermitian. However, only one triangle of 
        // the matrix actually needs to be filled, the symmetry can be implicit.
        //
        const Int colShift = S.ColShift(); // first row we own
        const Int rowShift = S.RowShift(); // first col we own
        const Int colStride = S.ColStride();
        const Int rowStride = S.RowStride();
        const Int localHeight = S.LocalHeight();
        const Int localWidth = S.LocalWidth();
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            {
                // Our process owns the rows colShift:colStride:n,
                //           and the columns rowShift:rowStride:n
                const Int i = colShift + iLoc*colStride;
                const Int j = rowShift + jLoc*rowStride;
                S.SetLocal( iLoc, jLoc, C(i-j,i+j) );
            }
        }

        // Make a backup of S before we overwrite it within the eigensolver
        auto SCopy( S );

        // Call the eigensolver. We first create an empty complex eigenvector 
        // matrix, X[MC,MR], and an eigenvalue column vector, w[VR,* ]
        //
        // Optional: set blocksizes and algorithmic choices here. See the 
        //           'Tuning' section of the README for details.
        DistMatrix<Real,VR,STAR> wImag;
        DistMatrix<C> X;
        SkewHermitianEig( LOWER, S, wImag, X, ASCENDING );

        if( print )
        {
            Print( SCopy, "S" );
            Print( X, "X" );
            Print( wImag, "wImag" );
        }

        // Check the residual, || S X - Omega X ||_F
        const Real frobS = HermitianFrobeniusNorm( LOWER, SCopy );
        auto E( X );
        Scale( C(0,1), E );
        DiagonalScale( RIGHT, NORMAL, wImag, E );
        Gemm( NORMAL, NORMAL, C(-1), SCopy, X, C(1), E );
        const Real frobResid = FrobeniusNorm( E );

        // Check the orthogonality of X
        Identity( E, n, n );
        Herk( LOWER, NORMAL, C(-1), X, C(1), E );
        const Real frobOrthog = HermitianFrobeniusNorm( LOWER, E );

        if( mpi::WorldRank() == 0 )
        {
            std::cout << "|| H ||_F = " << frobS << "\n"
                      << "|| H X - X Omega ||_F / || A ||_F = " 
                      << frobResid / frobS << "\n"
                      << "|| X X^H - I ||_F = " << frobOrthog / frobS
                      << "\n" << std::endl;
        }
    }
    catch( exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
