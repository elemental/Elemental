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
        const El::Int k = El::Input("--order","generate 2^k x 2^k matrix",4);
        const bool binary = El::Input("--binary","binary data?",false);
        const bool display = El::Input("--display","display matrix?",false);
        const bool print = El::Input("--print","print matrix?",true);
        El::ProcessInput();
        El::PrintInputReport();

        const El::Grid grid( comm );

        // Generate a Walsh matrix of order k (a 2^k x 2^k matrix)
        El::DistMatrix<double> W(grid);
        El::Walsh( W, k, binary );
        if( display )
            El::Display( W, "Walsh matrix" );
        if( print )
            El::Print( W, "W(2^k)");

        if( !binary )
        {
            El::LDL( W, true );
            auto d = El::GetDiagonal(W);
            El::MakeTrapezoidal( El::LOWER, W );
            El::FillDiagonal( W, 1. );

            if( display )
            {
                El::Display( W, "Lower factor" );
                El::Display( d, "Diagonal factor" );
            }
            if( print )
            {
                El::Print( W, "L" );
                El::Print( d, "d" );
            }
        }
    }
    catch( std::exception& e ) { El::ReportException(e); }

    return 0;
}
