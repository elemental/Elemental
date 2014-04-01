/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include ELEM_HEMV_INC
#include ELEM_QUADRATICPROGRAM_INC
#include ELEM_HERMITIANUNIFORMSPECTRUM_INC
#include ELEM_GAUSSIAN_INC
using namespace elem;

// This driver is an adaptation of the solver described at
//    http://www.stanford.edu/~boyd/papers/admm/quadprog/quadprog.html
// which is derived from the distributed ADMM article of Boyd et al.
//
// This example attempts to solve the following convex optimization problem:
//     minimize    (1/2) x' P x + q' x 
//     subject to  lb <= x <= ub
//

typedef double Real;

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );

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

        DistMatrix<Real> P, q, xTrue;
        HermitianUniformSpectrum( P, n, lbEig, ubEig );
        // Alternate the entries of xTrue between ub and lb
        Zeros( xTrue, n, 1 );
        if( xTrue.LocalWidth() == 1 )
            for( Int iLoc=0; iLoc<xTrue.LocalHeight(); ++iLoc )
                xTrue.SetLocal( iLoc, 0, 
                    ( xTrue.GlobalRow(iLoc)%2==0 ? lb : ub ) );
        // Set q := - P xTrue - du + dl
        Zeros( q, n, 1 );
        Hemv( LOWER, Real(-1), P, xTrue, Real(0), q );
        if( q.LocalWidth() == 1 )
            for( Int iLoc=0; iLoc<q.LocalHeight(); ++iLoc )
                q.UpdateLocal( iLoc, 0,
                    ( q.GlobalRow(iLoc)%2==0 ? 0.5 : -0.5 ) );
        if( print )
        {
            Print( P, "P" );
            Print( q, "q" );
            Print( xTrue, "xTrue" );
        }
        if( display )
            Display( P, "P" );

        DistMatrix<Real> x, z, u;
        QuadraticProgram
        ( P, q, lb, ub, x, z, u, rho, alpha, maxIter, absTol, relTol, inv, 
          progress );

        if( print )
        {
            Print( x, "x" );
            Print( z, "z" );
            Print( u, "u" );
        }
    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
