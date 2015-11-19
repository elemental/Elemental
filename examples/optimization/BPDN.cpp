/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

// This driver calls an adaptation of the solver described at
//    http://www.stanford.edu/~boyd/papers/admm/lasso/lasso.html
//
// The Least Absolute Shrinkage and Selection Operator (LASSO)
//   minimizes || A x - b ||_2^2 + lambda || x ||_1,
// which is equivalent to minimizing || A x - b ||_2 subject to 
// || x ||_1 <= t for some t >= 0.

typedef double Real;

int 
main( int argc, char* argv[] )
{
    Environment env( argc, argv );

    try
    {
        // TODO: Extend to add options for IPM
        const Int m = Input("--m","height of matrix",100);
        const Int n = Input("--n","width of matrix",200);
        const Int maxIter = Input("--maxIter","maximum # of iter's",500);
        const Real lambda = Input("--lambda","one-norm coefficient",1.);
        const Real rho = Input("--rho","augmented Lagrangian param.",1.);
        const Real alpha = Input("--alpha","over-relaxation",1.2);
        const Real absTol = Input("--absTol","absolute tolerance",1e-6);
        const Real relTol = Input("--relTol","relative tolerance",1e-4);
        const bool inv = Input("--inv","use explicit inverse",true);
        const bool progress = Input("--progress","print progress?",true);
        const bool display = Input("--display","display matrices?",false);
        const bool useIPM = Input("--useIPM","use Interior Point?",true);
        const bool print = Input("--print","print matrices",false);
        ProcessInput();
        PrintInputReport();

        DistMatrix<Real> A, b;
        Uniform( A, m, n );
        Uniform( b, m, 1 );
        if( print )
        {
            Print( A, "A" );
            Print( b, "b" );
        }
        if( display )
            Display( A, "A" );

        BPDNCtrl<Real> ctrl;
        ctrl.useIPM = useIPM;
        ctrl.ipmCtrl.mehrotraCtrl.print = progress;
        ctrl.admmCtrl.rho = rho;
        ctrl.admmCtrl.alpha = alpha;
        ctrl.admmCtrl.maxIter = maxIter;
        ctrl.admmCtrl.absTol = absTol;
        ctrl.admmCtrl.relTol = relTol;
        ctrl.admmCtrl.inv = inv;
        ctrl.admmCtrl.progress = progress;

        DistMatrix<Real> z;
        BPDN( A, b, lambda, z, ctrl );
        if( print )
            Print( z, "z" );
        const Real zOneNorm = OneNorm( z );
        const Int  zZeroNorm = ZeroNorm( z );
        const Real bTwoNorm = FrobeniusNorm( b );
        Gemv( NORMAL, Real(-1), A, z, Real(1), b );
        const Real rTwoNorm = FrobeniusNorm( b );
        if( mpi::Rank() == 0 )
            Output
            ("|| A z - b ||_2 = ",rTwoNorm,"\n",
             "|| b ||_2 = ",bTwoNorm,"\n",
             "|| z ||_1 = ",zOneNorm,"\n",
             "|| z ||_0 = ",zZeroNorm,"\n");
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
