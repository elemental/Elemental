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

        const El::Grid grid(comm);
        El::DistMatrix<double> J(grid);
        El::Legendre( J, n );
        if( display )
        {
            El::Display( J, "Jacobi matrix for Legendre polynomials" );
#ifdef EL_HAVE_QT5
            El::Spy( J, "Spy plot for Jacobi matrix" );
#endif
        }
        if( print )
            El::Print( J, "Jacobi matrix for Legendre polynomials" );

        // We will compute Gaussian quadrature points and weights over [-1,+1]
        // using the eigenvalue decomposition of the Jacobi matrix for the
        // Legendre polynomials.
        El::DistMatrix<double> points(grid), X(grid);
        El::HermitianTridiagEig
        ( GetDiagonal(J), GetDiagonal(J,-1), points, X );
        if( display )
            El::Display( points, "Quadrature points" );
        if( print )
            El::Print( points, "points" );
        auto firstRow = X( El::IR(0,1), El::ALL );
        El::DistMatrix<double> weights( firstRow );
        auto entryToWeight =
          []( const double& gamma ) { return 2*gamma*gamma; };
        El::EntrywiseMap( weights, El::MakeFunction(entryToWeight) );
        if( display )
            El::Display( weights, "Quadrature weights" );
        if( print )
            El::Print( weights, "weights" );
    }
    catch( std::exception& e ) { El::ReportException(e); }

    return 0;
}
