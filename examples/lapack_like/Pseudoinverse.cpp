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

        const El::Int m = El::Input("--height","height of matrix",100);
        const El::Int n = El::Input("--width","width of matrix",100);
        const bool print = El::Input("--print","print matrices?",false);
        El::ProcessInput();
        El::PrintInputReport();

        El::DistMatrix<Scalar> A;
        El::Uniform( A, m, n );

        El::Timer timer;
        // Compute the pseudoinverseof A (but do not overwrite A)
        El::DistMatrix<Scalar> pinvA( A );
        if( El::mpi::Rank(comm) == 0 )
            timer.Start();
        El::Pseudoinverse( pinvA );
        if( El::mpi::Rank(comm) == 0 )
            timer.Stop();
        if( print )
        {
            El::Print( A, "A" );
            El::Print( pinvA, "pinv(A)" );
        }

        const Real frobA = El::FrobeniusNorm( A );
        const Real frobPinvA = El::FrobeniusNorm( pinvA );

        if( El::mpi::Rank(comm) == 0 )
        {
            El::Output("PseudoInverse time: ",timer.Total()," secs");
            El::Output
            ("||   A     ||_F = ",frobA,"\n",
             "|| pinv(A) ||_F = ",frobPinvA,"\n");
        }
    }
    catch( std::exception& e ) { El::ReportException(e); }

    return 0;
}
