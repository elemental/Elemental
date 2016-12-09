/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

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
    El::Environment env( argc, argv );

    try
    {
        const El::Int n = El::Input("--n","problem size",200);
        const El::Int maxIter =
          El::Input("--maxIter","maximum # of iter's",500);
        const Real lb = El::Input("--lb","lower bound for x",0.5);
        const Real ub = El::Input("--ub","upper bound for x",1.0);
        const Real lbEig = El::Input("--lbEig","spectral lower bound",1.);
        const Real ubEig = El::Input("--ubEig","spectral upper bound",2.);
        const Real rho = El::Input("--rho","augmented Lagrangian param.",1.);
        const Real alpha = El::Input("--alpha","over-relaxation",1.2);
        const Real absTol = El::Input("--absTol","absolute tolerance",1e-6);
        const Real relTol = El::Input("--relTol","relative tolerance",1e-4);
        const bool inv = El::Input("--inv","form inv(LU) to avoid trsv?",true);
        const bool progress = El::Input("--progress","print progress?",true);
        const bool display = El::Input("--display","display matrices?",false);
        const bool print = El::Input("--print","print matrices",false);
        El::ProcessInput();
        El::PrintInputReport();

        El::ADMMCtrl<Real> ctrl;
        ctrl.rho = rho;
        ctrl.alpha = alpha;
        ctrl.maxIter = maxIter;
        ctrl.absTol = absTol;
        ctrl.relTol = relTol;
        ctrl.inv = inv;
        ctrl.print = progress;

        El::DistMatrix<Real> Q, c, xTrue;
        El::HermitianUniformSpectrum( Q, n, lbEig, ubEig );
        // Alternate the entries of xTrue between ub and lb
        El::Zeros( xTrue, n, 1 );
        if( xTrue.LocalWidth() == 1 )
            for( El::Int iLoc=0; iLoc<xTrue.LocalHeight(); ++iLoc )
                xTrue.SetLocal( iLoc, 0,
                    ( xTrue.GlobalRow(iLoc)%2==0 ? lb : ub ) );
        // Set c := - Q xTrue - du + dl
        El::Zeros( c, n, 1 );
        El::Hemv( El::LOWER, Real(-1), Q, xTrue, Real(0), c );
        if( c.LocalWidth() == 1 )
            for( El::Int iLoc=0; iLoc<c.LocalHeight(); ++iLoc )
                c.UpdateLocal( iLoc, 0,
                    ( c.GlobalRow(iLoc)%2==0 ? 0.5 : -0.5 ) );
        if( print )
        {
            El::Print( Q, "Q" );
            El::Print( c, "c" );
            El::Print( xTrue, "xTrue" );
        }
        if( display )
            El::Display( Q, "Q" );

        El::Timer timer;
        El::DistMatrix<Real> z;
        if( El::mpi::Rank() == 0 )
            timer.Start();
        El::qp::box::ADMM( Q, c, lb, ub, z, ctrl );
        if( El::mpi::Rank() == 0 )
            timer.Stop();

        if( print )
            El::Print( z, "z" );
        if( El::mpi::Rank() == 0 )
            El::Output("QPBox time: ",timer.Total(),"secs");
    }
    catch( std::exception& e ) { El::ReportException(e); }

    return 0;
}
