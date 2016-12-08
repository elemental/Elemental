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
    El::Environment env( argc, argv );
    El::mpi::Comm comm = El::mpi::COMM_WORLD;

    try
    {
        typedef double Real;
        typedef El::Complex<Real> Scalar;

        const El::Int m = El::Input("--height","matrix height",100);
        const El::Int n = El::Input("--width","matrix width",100);
        const bool print = El::Input("--print","print matrices?",false);
        El::ProcessInput();
        El::PrintInputReport();

        const El::Grid grid(comm);
        El::DistMatrix<Scalar> A(grid), Q(grid), P(grid);
        El::Uniform( A, m, n );

        // Compute the polar decomp of A (but do not overwrite A)
        Q = A;
        El::Polar( Q, P );
        if( print )
        {
            El::Print( A, "A" );
            El::Print( Q, "Q" );
            El::Print( P, "P" );
        }
    }
    catch( std::exception& e ) { El::ReportException(e); }

    return 0;
}
