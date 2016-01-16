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
        const Int m = Input("--m","height of matrix",100);
        const Int n = Input("--n","width of matrix",200);
        const bool useIPM = Input("--useIPM","use Interior Point?",true);
        // TODO: Add options for controlling IPM
        const Int maxIter = Input("--maxIter","maximum # of iter's",500);
        const Real rho = Input("--rho","augmented Lagrangian param.",1.);
        const Real alpha = Input("--alpha","over-relaxation",1.2);
        const Real absTol = Input("--absTol","absolute tolerance",1e-6);
        const Real relTol = Input("--relTol","relative tolerance",1e-4);
        const Real pinvTol = Input("--pinvTol","pseudoinv tolerance",0.);
        const bool usePinv = Input("--usePinv","Directly compute pinv(A)",true);
        const bool progress = Input("--progress","print progress?",true);
        const bool display = Input("--display","display matrices?",false);
        const bool print = Input("--print","print matrices",false);
        ProcessInput();
        PrintInputReport();

        DistMatrix<Real> A, b, xTrue;
        Uniform( A, m, n );
        Uniform( b, m, 1 );
        if( print )
        {
            Print( A, "A" );
            Print( b, "b" );
        }
        if( display )
            Display( A, "A" );

        const bool sparse = false;
        BPCtrl<Real> ctrl(sparse);
        ctrl.useIPM = useIPM;
        ctrl.admmCtrl.rho = rho;
        ctrl.admmCtrl.alpha = alpha;
        ctrl.admmCtrl.maxIter = maxIter;
        ctrl.admmCtrl.absTol = absTol;
        ctrl.admmCtrl.relTol = relTol;
        ctrl.admmCtrl.usePinv = usePinv;
        ctrl.admmCtrl.pinvTol = pinvTol;
        ctrl.admmCtrl.progress = progress;

        DistMatrix<Real> x;
        Timer timer;
        if( mpi::Rank() == 0 )
            timer.Start();
        BP( A, b, x, ctrl );
        if( mpi::Rank() == 0 )
            timer.Stop();
        if( print )
            Print( x, "x" );
        const Int xZeroNorm = ZeroNorm( x );
        if( mpi::Rank() == 0 )
        {
            Output("Basis Pursuit time: ",timer.Total()," secs");
            Output("|| x ||_0 = ",xZeroNorm);
        }
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
