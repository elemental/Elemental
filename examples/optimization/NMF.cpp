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

int 
main( int argc, char* argv[] )
{
    Environment env( argc, argv );

    try
    {
        const Int m = Input("--m","matrix height",50);
        const Int n = Input("--n","matrix width",50);
        const Int k = Input("--k","rank of approximation",3);
        const Int maxIter = Input("--maxIter","max. iterations",20);
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


        NMFCtrl<Real> ctrl;
        ctrl.nnlsCtrl.approach = NNLS_QP;
        ctrl.nnlsCtrl.qpCtrl.mehrotraCtrl.print = false;
        ctrl.nnlsCtrl.qpCtrl.mehrotraCtrl.time = false;
        ctrl.nnlsCtrl.socpCtrl.mehrotraCtrl.print = false;
        ctrl.nnlsCtrl.socpCtrl.mehrotraCtrl.time = false;
        ctrl.maxIter = maxIter;

        Timer timer;
        DistMatrix<Real> Y;
        if( mpi::Rank() == 0 )
            timer.Start();
        NMF( A, X, Y, ctrl );
        if( mpi::Rank() == 0 )
            timer.Stop();

        if( print )
        {
            Print( X, "X" );
            Print( Y, "Y" );
        }

        const Real ANorm = FrobeniusNorm( A );
        Gemm( NORMAL, ADJOINT, Real(-1), X, Y, Real(1), A );
        const Real ENorm = FrobeniusNorm( A );
        if( mpi::Rank() == 0 )
        {
            Output("NMF time: ",timer.Total()," secs");
            Output("|| A - X Y^H ||_F / || A ||_F = ",ENorm/ANorm);
        }
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
