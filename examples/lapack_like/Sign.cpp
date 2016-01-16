/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

typedef double Real;
typedef Complex<Real> C;

int
main( int argc, char* argv[] )
{
    Environment env( argc, argv );

    try 
    {
        const Int m = Input("--height","height of matrix",100);
        const Int n = Input("--width","width of matrix",100);
        const SignScaling scaling = 
            static_cast<SignScaling>(Input("--scaling","scaling strategy",0));
        const Int maxIts = Input("--maxIts","max number of iter's",100);
        const double tol = Input("--tol","convergence tolerance",1e-6);
        const bool progress = Input("--progress","print sign progress?",true);
        const bool print = Input("--print","print matrix?",false);
        const bool display = Input("--display","display matrix?",false);
        ProcessInput();
        PrintInputReport();

        DistMatrix<C> A;
        Uniform( A, m, n );
        if( print )
            Print( A, "A" );
        if( display )
            Display( A, "A" );

        SignCtrl<Real> signCtrl;
        signCtrl.maxIts = maxIts;
        signCtrl.tol = tol;
        signCtrl.progress = progress;
        signCtrl.scaling = scaling;

        Timer timer;
        // Compute sgn(A)
        if( mpi::Rank() == 0 )
            timer.Start();
        Sign( A, signCtrl );
        if( mpi::Rank() == 0 )
            timer.Stop();
        if( print )
            Print( A, "A" );
        if( display )
            Display( A, "A" );
        if( mpi::Rank() == 0 )
            Output("Sign time: ",timer.Total()," secs");
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
