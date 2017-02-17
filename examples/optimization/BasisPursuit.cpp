/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

typedef double Real;

int
main( int argc, char* argv[] )
{
    El::Environment env( argc, argv );

    try
    {
        const El::Int m = El::Input("--m","height of matrix",100);
        const El::Int n = El::Input("--n","width of matrix",200);
        const bool useIPM = El::Input("--useIPM","use Interior Point?",true);
        // TODO(poulson): Add options for controlling IPM
        const El::Int maxIter =
          El::Input("--maxIter","maximum # of iter's",500);
        const Real rho = El::Input("--rho","augmented Lagrangian param.",1.);
        const Real alpha = El::Input("--alpha","over-relaxation",1.2);
        const Real absTol = El::Input("--absTol","absolute tolerance",1e-6);
        const Real relTol = El::Input("--relTol","relative tolerance",1e-4);
        const Real pinvTol = El::Input("--pinvTol","pseudoinv tolerance",0.);
        const bool usePinv =
          El::Input("--usePinv","Directly compute pinv(A)",true);
        const bool progress = El::Input("--progress","print progress?",true);
        const bool display = El::Input("--display","display matrices?",false);
        const bool print = El::Input("--print","print matrices",false);
        El::ProcessInput();
        El::PrintInputReport();

        El::DistMatrix<Real> A, b, xTrue;
        El::Uniform( A, m, n );
        El::Uniform( b, m, 1 );
        if( print )
        {
            El::Print( A, "A" );
            El::Print( b, "b" );
        }
        if( display )
            El::Display( A, "A" );

        const bool sparse = false;
        El::BPCtrl<Real> ctrl(sparse);
        ctrl.useIPM = useIPM;
        ctrl.admmCtrl.rho = rho;
        ctrl.admmCtrl.alpha = alpha;
        ctrl.admmCtrl.maxIter = maxIter;
        ctrl.admmCtrl.absTol = absTol;
        ctrl.admmCtrl.relTol = relTol;
        ctrl.admmCtrl.usePinv = usePinv;
        ctrl.admmCtrl.pinvTol = pinvTol;
        ctrl.admmCtrl.progress = progress;
        ctrl.lpIPMCtrl.mehrotraCtrl.print = true;

        El::DistMatrix<Real> x;
        El::Timer timer;
        if( El::mpi::Rank() == 0 )
            timer.Start();
        El::BP( A, b, x, ctrl );
        if( El::mpi::Rank() == 0 )
            timer.Stop();
        if( print )
            El::Print( x, "x" );
        const El::Int xZeroNorm = El::ZeroNorm( x );
        if( El::mpi::Rank() == 0 )
        {
            El::Output("Basis Pursuit time: ",timer.Total()," secs");
            El::Output("|| x ||_0 = ",xZeroNorm);
        }
        SoftThreshold( x, El::Sqrt(El::limits::Epsilon<Real>()) );
        if( print )
            El::Print( x, "xThresh" );
        const El::Int xZeroNormThresh = El::ZeroNorm( x );
        if( El::mpi::Rank() == 0 )
        {
            El::Output("|| xThresh ||_0 = ",xZeroNormThresh);
        }
    }
    catch( std::exception& e ) { El::ReportException(e); }

    return 0;
}
