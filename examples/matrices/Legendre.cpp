/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

int 
main( int argc, char* argv[] )
{
    Environment env( argc, argv );

    try
    {
        const Int n = Input("--size","size of matrix",10);
        const bool display = Input("--display","display matrix?",true);
        const bool print = Input("--print","print matrix?",false);
        ProcessInput();
        PrintInputReport();

        DistMatrix<double> J;
        Legendre( J, n );
        if( display )
        {
            Display( J, "Jacobi matrix for Legendre polynomials" );
#ifdef EL_HAVE_QT5
            Spy( J, "Spy plot for Jacobi matrix" );
#endif
        }
        if( print )
            Print( J, "Jacobi matrix for Legendre polynomials" );

        // We will compute Gaussian quadrature points and weights over [-1,+1]
        // using the eigenvalue decomposition of the Jacobi matrix for the 
        // Legendre polynomials.
        DistMatrix<double,VR,  STAR> points;
        DistMatrix<double,STAR,VR  > X;
        HermitianTridiagEig
        ( GetDiagonal(J), GetDiagonal(J,-1), points, X, ASCENDING );
        if( display )
            Display( points, "Quadrature points" );
        if( print )
            Print( points, "points" );
        auto firstRow = View( X, 0, 0, 1, n );
        DistMatrix<double,STAR,STAR> weights = firstRow;
        for( Int j=0; j<n; ++j )
        {
            const double gamma = weights.Get( 0, j );
            weights.Set( 0, j, 2*gamma*gamma );
        }
        if( display )
            Display( weights, "Quadrature weights" );
        if( print )
            Print( weights, "weights" );
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
