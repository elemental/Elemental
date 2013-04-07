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
#include "elemental/matrices/Legendre.hpp"
#include "elemental/matrices/Zeros.hpp"
using namespace elem;

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::CommRank( comm );

    try
    {
        const int n = Input("--size","size of matrix",10);
        const bool print = Input("--print","print matrix?",true);
        ProcessInput();
        PrintInputReport();

        DistMatrix<double> J;
        Legendre( n, J );

#ifndef HAVE_PMRRR
        if( print )
            J.Print("Jacobi matrix for Legendre polynomials");
#else
        // This will perform a lot of unnecessary work, but the code is simpler
        // than directly calling PMRRR
        //
        // We will compute Gaussian quadrature points and weights over [-1,+1]
        // using the eigenvalue decomposition of the Jacobi matrix for the 
        // Legendre polynomials.
        //
        DistMatrix<double,VR,STAR> points;
        DistMatrix<double> X;
        HermitianEig( LOWER, J, points, X );
        hermitian_eig::Sort( points, X );
        if( print )
            points.Print("points");
        DistMatrix<double> firstRow;
        View( firstRow, X, 0, 0, 1, n );
        DistMatrix<double,STAR,STAR> weights = firstRow;
        for( int j=0; j<n; ++j )
        {
            const double gamma = weights.Get( 0, j );
            weights.Set( 0, j, 2*gamma*gamma );
        }
        if( print )
            weights.Print("weights");
#endif
    }
    catch( ArgException& e )
    {
        // There is nothing to do
    }
    catch( std::exception& e )
    {
        std::ostringstream os;
        os << "Process " << commRank << " caught error message:\n" << e.what()
           << std::endl;
#ifndef RELEASE
        DumpCallStack();
#endif
    }

    Finalize();
    return 0;
}
