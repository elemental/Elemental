/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "El.hpp" instead
#include "El-lite.hpp"

using namespace El;

typedef double Real;

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );

    try
    {
        const Int m = Input("--m","matrix height",200);
        const Int n = Input("--n","matrix width",200);
        const Int k = Input("--k","rank of approximation",10);
        const bool display = Input("--display","display matrices?",false);
        const bool print = Input("--print","print matrices",false);
        ProcessInput();
        PrintInputReport();

        DistMatrix<Real> A, X;
        Uniform( A, m, n );
        LowerClip( A );
        Uniform( X, m, k );
        LowerClip( X );
        if( print )
        {
            Print( A, "A" );
            Print( X, "X" );
        }
        if( display )
            Display( A, "A" );

        DistMatrix<Real> Y;
        NMF( A, X, Y );

        if( print )
        {
            Print( X, "X" );
            Print( Y, "Y" );
        }

        const double ANorm = FrobeniusNorm( A );
        Gemm( NORMAL, ADJOINT, Real(-1), X, Y, Real(1), A );
        const double ENorm = FrobeniusNorm( A );
        if( mpi::WorldRank() == 0 )
            std::cout << "|| A - X Y^H ||_F / || A ||_F = " << ENorm/ANorm
                      << std::endl;
    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
