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
        const Int m = Input("--m","height of C",1000);
        const Int n = Input("--n","width of C",1000);
        const Int k = Input("--k","inner dimension",1000);
        const double alpha = Input("--alpha","scale of A B",2.);
        const double beta = Input("--beta","scale of C",3.);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        Grid g( mpi::COMM_WORLD );
        if( g.Height() != g.Width() )
            LogicError("This routine requires a square grid");

        DistMatrix<double> A(g), B(g), C(g);
        Uniform( A, m, k );
        Uniform( B, k, n );
        Uniform( C, m, n );
        if( print )
        {
            Print( A, "A" );
            Print( B, "B" );
            Print( C, "C" );
        }

        Timer timer;
        timer.Start();
        Gemm( NORMAL, NORMAL, alpha, A, B, beta, C, GEMM_CANNON );
        const double gemmTime = timer.Stop();
        if( mpi::Rank() == 0 )
        {
            const double gFlops = (2.*m*n*k) / (1.e9*gemmTime);
            Output("Gemm took ",gemmTime," seconds (",gFlops," GFlop/s)");
        }
        if( print )
            Print( C, "C := alpha A B + beta C" );
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
