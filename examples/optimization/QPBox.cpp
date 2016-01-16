/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

// This driver is an adaptation of the solver described at
//    http://www.stanford.edu/~boyd/papers/admm/quadprog/quadprog.html
// which is derived from the distributed ADMM article of Boyd et al.
//
// This example attempts to solve the following convex optimization problem:
//     minimize    (1/2) x' Q x + c' x 
//     subject to  lb <= x <= ub
//

typedef double Real;

int main( int argc, char* argv[] )
{
    Environment env( argc, argv );

    try
    {
        const Int n = Input("--n","problem size",200);
        const Int maxIter = Input("--maxIter","maximum # of iter's",500);
        const Real lb = Input("--lb","lower bound for x",0.5);
        const Real ub = Input("--ub","upper bound for x",1.0);
        const Real lbEig = Input("--lbEig","spectral lower bound",1.);
        const Real ubEig = Input("--ubEig","spectral upper bound",2.);
        const Real rho = Input("--rho","augmented Lagrangian param.",1.);
        const Real alpha = Input("--alpha","over-relaxation",1.2);
        const Real absTol = Input("--absTol","absolute tolerance",1e-6);
        const Real relTol = Input("--relTol","relative tolerance",1e-4);
        const bool inv = Input("--inv","form inv(LU) to avoid trsv?",true);
        const bool progress = Input("--progress","print progress?",true);
        const bool display = Input("--display","display matrices?",false);
        const bool print = Input("--print","print matrices",false);
        ProcessInput();
        PrintInputReport();

        ADMMCtrl<Real> ctrl;
        ctrl.rho = rho;
        ctrl.alpha = alpha;
        ctrl.maxIter = maxIter;
        ctrl.absTol = absTol;
        ctrl.relTol = relTol;
        ctrl.inv = inv;
        ctrl.print = progress;

        DistMatrix<Real> Q, c, xTrue;
        HermitianUniformSpectrum( Q, n, lbEig, ubEig );
        // Alternate the entries of xTrue between ub and lb
        Zeros( xTrue, n, 1 );
        if( xTrue.LocalWidth() == 1 )
            for( Int iLoc=0; iLoc<xTrue.LocalHeight(); ++iLoc )
                xTrue.SetLocal( iLoc, 0, 
                    ( xTrue.GlobalRow(iLoc)%2==0 ? lb : ub ) );
        // Set c := - Q xTrue - du + dl
        Zeros( c, n, 1 );
        Hemv( LOWER, Real(-1), Q, xTrue, Real(0), c );
        if( c.LocalWidth() == 1 )
            for( Int iLoc=0; iLoc<c.LocalHeight(); ++iLoc )
                c.UpdateLocal( iLoc, 0,
                    ( c.GlobalRow(iLoc)%2==0 ? 0.5 : -0.5 ) );
        if( print )
        {
            Print( Q, "Q" );
            Print( c, "c" );
            Print( xTrue, "xTrue" );
        }
        if( display )
            Display( Q, "Q" );

        Timer timer;
        DistMatrix<Real> z;
        if( mpi::Rank() == 0 )
            timer.Start();
        qp::box::ADMM( Q, c, lb, ub, z, ctrl );
        if( mpi::Rank() == 0 )
            timer.Stop();

        if( print )
            Print( z, "z" );
        if( mpi::Rank() == 0 )
            Output("QPBox time: ",timer.Total(),"secs");
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
