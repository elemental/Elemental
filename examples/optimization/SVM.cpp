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
        const bool useIPM = Input("--useIPM","use Interior Point?",true);
        const Int maxIter = Input("--maxIter","maximum # of iter's",500);
        const Real gamma = Input("--gamma","two-norm coefficient",1.);
        const Real rho = Input("--rho","augmented Lagrangian param.",1.);
        const bool inv = Input("--inv","use explicit inverse",true);
        const bool progress = Input("--progress","print progress?",true);
        const bool display = Input("--display","display matrices?",false);
        const bool print = Input("--print","print matrices",false);
        ProcessInput();
        PrintInputReport();

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

        SVMCtrl<Real> ctrl;
        ctrl.useIPM = useIPM;
        ctrl.modelFitCtrl.rho = rho;
        ctrl.modelFitCtrl.maxIter = maxIter;
        ctrl.modelFitCtrl.inv = inv;
        ctrl.modelFitCtrl.progress = progress;

        Timer timer;
        DistMatrix<Real> wHatSVM;
        if( mpi::Rank() == 0 )
            timer.Start();
        SVM( G, q, gamma, wHatSVM, ctrl );
        if( mpi::Rank() == 0 )
            timer.Stop();
        auto wSVM = View( wHatSVM, 0, 0, n, 1 );
        const Real offsetSVM = -wHatSVM.Get(n,0);
        const Real wSVMNorm = FrobeniusNorm( wSVM );
        if( mpi::Rank() == 0 )
        {
            Output("SVM time: ",timer.Total()," secs");
            Output
            ("|| wSVM ||_2=",wSVMNorm,"\n",
             "margin      =",Real(2)/wSVMNorm,"\n",
             "offsetSVM   =",offsetSVM,"\n",
             "offsetSVM / || wSVM ||_2=",offsetSVM/wSVMNorm);
        }
        if( print )
            Print( wSVM, "wSVM" );

        // Report the classification percentage using the returned hyperplane
        DistMatrix<Real> qSVM;
        Ones( qSVM, m, 1 );
        Gemv( NORMAL, Real(1), G, wSVM, -offsetSVM, qSVM );
        EntrywiseMap( qSVM, function<Real(Real)>(sgnMap) );
        if( print )
            Print( qSVM, "qSVM" );
        qSVM -= q;
        const Real numWrong = OneNorm(qSVM) / Real(2);
        if( mpi::Rank(comm) == 0 )
            Output("ratio misclassified: ",numWrong,"/",m);
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
