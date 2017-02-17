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
        const El::Int n = El::Input("--size","size of matrix",10);
        const bool display = El::Input("--display","display matrix?",false);
        const bool print = El::Input("--print","print matrix?",true);
        El::ProcessInput();
        El::PrintInputReport();

        const El::Grid grid( comm );

        const double& pi = El::Pi<double>();
        auto fourier = [&]( El::Int i, El::Int j ) { return (-2*pi*i*j)/n; };
        auto phase = [&]( El::Int i, El::Int j )
          {
            return fourier(i,j) +
              El::Sqrt(double(i)*double(i) + double(j)*double(j));
          };

        El::DistMatrix<El::Complex<double>> F(grid), G(grid);
        El::Egorov( F, El::MakeFunction(fourier), n );
        El::Egorov( G, El::MakeFunction(phase),   n );

        if( display )
        {
            El::Display( F, "Egorov with Fourier phase" );
            El::Display( G, "Egorov with more general phase" );
        }
        if( print )
        {
            El::Print( F, "Egorov with Fourier phase:" );
            El::Print( G, "Egorov with more general phase:" );
        }
    }
    catch( std::exception& e ) { El::ReportException(e); }

    return 0;
}
