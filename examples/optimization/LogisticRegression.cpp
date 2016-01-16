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
    mpi::Comm comm = mpi::COMM_WORLD;

    try
    {
        const Int m = Input("--numExamples","number of examples",200);
        const Int n = Input("--numFeatures","number of features",100);
        const Int maxIter = Input("--maxIter","maximum # of iter's",500);
        const Real gamma = Input("--gamma","two-norm coefficient",1.);
        const Int penaltyInt = Input("--penalty","0: none, 1: l1, 2: l2",1);
        const Real rho = Input("--rho","augmented Lagrangian param.",1.);
        const bool inv = Input("--inv","use explicit inverse",true);
        const bool progress = Input("--progress","print progress?",true);
        const bool display = Input("--display","display matrices?",false);
        const bool print = Input("--print","print matrices",false);
        ProcessInput();
        PrintInputReport();

        Regularization penalty = static_cast<Regularization>(penaltyInt);

        // Define a random (affine) hyperplane 
        DistMatrix<Real> w;
        Gaussian( w, n, 1 );
        {
            const Real wNorm = FrobeniusNorm( w );
            w *= 1/wNorm;
        }
        Real offset = SampleNormal();
        mpi::Broadcast( offset, 0, comm );

        // Draw each example (row) from a Gaussian perturbation of a point
        // lying on the hyperplane
        DistMatrix<Real> G;
        Gaussian( G, m, n );
        DistMatrix<Real,MR,STAR> w_MR_STAR;
        w_MR_STAR.AlignWith( G );
        w_MR_STAR = w;
        for( Int jLoc=0; jLoc<G.LocalWidth(); ++jLoc )
            for( Int iLoc=0; iLoc<G.LocalHeight(); ++iLoc )
                G.UpdateLocal( iLoc, jLoc, w_MR_STAR.GetLocal(jLoc,0)*offset );
        
        // Label each example based upon its location relative to the hyperplane
        DistMatrix<Real> q;
        Ones( q, m, 1 );
        Gemv( NORMAL, Real(1), G, w, -offset, q );
        auto sgnMap = []( Real alpha ) 
                      { return alpha >= 0 ? Real(1) : Real(-1); }; 
        EntrywiseMap( q, function<Real(Real)>(sgnMap) );

        if( mpi::Rank(comm) == 0 )
            Output("offset=",offset);
        if( print )
        {
            Print( w, "w" );
            Print( G, "G" );
            Print( q, "q" );
        }
        if( display )
            Display( G, "G" );

        ModelFitCtrl<Real> ctrl;
        ctrl.rho = rho;
        ctrl.maxIter = maxIter;
        ctrl.inv = inv;
        ctrl.progress = progress;

        Timer timer;
        DistMatrix<Real> wHatLog;
        if( mpi::Rank() == 0 )
            timer.Start();
        LogisticRegression( G, q, wHatLog, gamma, penalty );
        if( mpi::Rank() == 0 )
            timer.Stop();
        auto wLog = View( wHatLog, 0, 0, n, 1 );
        const Real offsetLog = -wHatLog.Get(n,0);
        const Real wLogOneNorm = OneNorm( wLog );
        const Real wLogFrobNorm = FrobeniusNorm( wLog );
        if( mpi::Rank(comm) == 0 )
        {
            Output("Logistic Regression time: ",timer.Total()," secs");
            Output("|| wLog ||_1=",wLogOneNorm);
            Output("|| wLog ||_2=",wLogFrobNorm);
            Output("margin      =",Real(2)/wLogFrobNorm);
            Output("offsetLog=",offsetLog);
            Output("offsetLog / || wLog ||_2=",offsetLog/wLogFrobNorm);
            Output("");
        }
        if( print )
            Print( wLog, "wLog" );

        // Report the classification percentage using the returned hyperplane
        // TODO: Evaluate probabilities implied by logistic function
        DistMatrix<Real> qLog;
        Ones( qLog, m, 1 );
        Gemv( NORMAL, Real(1), G, wLog, -offsetLog, qLog );
        EntrywiseMap( qLog, function<Real(Real)>(sgnMap) );
        if( print )
            Print( qLog, "qLog" );
        qLog -= q;
        const Real numWrong = OneNorm(qLog) / Real(2);
        if( mpi::Rank(comm) == 0 )
            Output("ratio misclassified: ",numWrong,"/",m);
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
