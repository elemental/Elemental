/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

int
main( int argc, char* argv[] )
{
    // This detects whether or not you have already initialized MPI and
    // does so if necessary.
    El::Environment env( argc, argv );
    El::mpi::Comm comm = El::mpi::COMM_WORLD;

    // Surround the Elemental calls with try/catch statements in order to
    // safely handle any exceptions that were thrown during execution.
    try
    {
        typedef double Real;
        typedef El::Complex<Real> Scalar;

        const El::Int n = El::Input("--size","size of matrix",100);
        const bool print = El::Input("--print","print matrices?",false);
        El::ProcessInput();
        El::PrintInputReport();

        // Create a 2d process grid from a communicator. In our case, it is
        // mpi::COMM_WORLD. There is another constructor that allows you to
        // specify the grid dimensions, Grid g( comm, r ), which creates a
        // grid of height r.
        El::Grid grid( comm );

        // Create an n x n complex distributed matrix,
        // We distribute the matrix using grid 'grid'.
        El::DistMatrix<Scalar> H( n, n, grid );

        // Manually fill entry (i,j) with the complex value (i+j,i-j), which
        // results in a Hermitian matrix.
        const El::Int localHeight = H.LocalHeight();
        const El::Int localWidth = H.LocalWidth();
        for( El::Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const El::Int j = H.GlobalCol(jLoc);
            for( El::Int iLoc=0; iLoc<localHeight; ++iLoc )
            {
                const El::Int i = H.GlobalRow(iLoc);
                H.SetLocal( iLoc, jLoc, Scalar(i+j,i-j) );
            }
        }

        // Make a backup of H before we overwrite it within the eigensolver
        auto HCopy( H );

        // Call the eigensolver. We first create an empty complex eigenvector
        // matrix, Q[MC,MR], and an eigenvalue column vector, w[VR,* ]
        //
        // Optional: set blocksizes and algorithmic choices here. See the
        //           'Tuning' section of the README for details.
        El::DistMatrix<Real,El::VR,El::STAR> w( grid );
        El::DistMatrix<Scalar> Q( grid );
        El::HermitianEig( El::LOWER, H, w, Q );

        if( print )
        {
            El::Print( HCopy, "H" );
            El::Print( Q, "Q" );
            El::Print( w, "w" );
        }

        // Check the residual, || H Q - Omega Q ||_F
        const Real frobH = El::HermitianFrobeniusNorm( El::LOWER, HCopy );
        auto E( Q );
        El::DiagonalScale( El::RIGHT, El::NORMAL, w, E );
        El::Hemm( El::LEFT, El::LOWER, Scalar(-1), HCopy, Q, Scalar(1), E );
        const Real frobResid = El::FrobeniusNorm( E );

        // Check the orthogonality of Q
        El::Identity( E, n, n );
        El::Herk( El::LOWER, El::ADJOINT, Real(-1), Q, Real(1), E );
        const Real frobOrthog = El::HermitianFrobeniusNorm( El::LOWER, E );

        if( El::mpi::Rank(comm) == 0 )
            El::Output
            ("|| H ||_F = ",frobH,"\n",
             "|| H Q - Q Omega ||_F / || A ||_F = ",frobResid/frobH,"\n",
             "|| Q' Q - I ||_F = ",frobOrthog,"\n");
    }
    catch( std::exception& e ) { El::ReportException(e); }

    return 0;
}
