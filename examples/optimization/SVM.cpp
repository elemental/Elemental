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
    El::mpi::Comm comm = El::mpi::COMM_WORLD;

    try
    {
        const El::Int m = El::Input("--numExamples","number of examples",200);
        const El::Int n = El::Input("--numFeatures","number of features",100);
        const double gamma = El::Input("--gamma","hinge-loss penalty",1.0);
        const bool display = El::Input("--display","display matrices?",false);
        const bool print = El::Input("--print","print matrices",false);
        El::ProcessInput();
        El::PrintInputReport();

        // Define a random (affine) hyperplane
        El::DistMatrix<Real> w;
        El::Gaussian( w, n, 1 );
        {
            const Real wNorm = El::FrobeniusNorm( w );
            w *= 1/wNorm;
        }
        Real offset = El::SampleNormal();
        El::mpi::Broadcast( offset, 0, comm );

        // Draw each example (row) from a Gaussian perturbation of a point
        // lying on the hyperplane
        El::DistMatrix<Real> G;
        El::Gaussian( G, m, n );
        El::DistMatrix<Real,El::MR,El::STAR> w_MR_STAR;
        w_MR_STAR.AlignWith( G );
        w_MR_STAR = w;
        for( El::Int jLoc=0; jLoc<G.LocalWidth(); ++jLoc )
            for( El::Int iLoc=0; iLoc<G.LocalHeight(); ++iLoc )
                G.UpdateLocal
                ( iLoc, jLoc, w_MR_STAR.GetLocal(jLoc,0)*offset );

        // Label each example
        El::DistMatrix<Real> q;
        El::Ones( q, m, 1 );
        El::Gemv( El::NORMAL, Real(1), G, w, -offset, q );
        auto sgnMap = []( const Real& alpha )
                      { return alpha >= 0 ? Real(1) : Real(-1); };
        El::EntrywiseMap( q, El::MakeFunction(sgnMap) );

        if( El::mpi::Rank(comm) == 0 )
            El::Output("offset=",offset);
        if( print )
        {
            El::Print( w, "w" );
            El::Print( G, "G" );
            El::Print( q, "q" );
        }
        if( display )
            El::Display( G, "G" );

        El::SVMCtrl<Real> ctrl;
        // TODO(poulson): Add support for configuring the IPM

        El::Timer timer;
        El::DistMatrix<Real> wHatSVM;
        if( El::mpi::Rank() == 0 )
            timer.Start();
        El::SVM( G, q, gamma, wHatSVM, ctrl );
        if( El::mpi::Rank() == 0 )
            timer.Stop();
        auto wSVM = wHatSVM( El::IR(0,n), El::IR(0,1) );
        const Real offsetSVM = -wHatSVM.Get(n,0);
        const Real wSVMNorm = El::FrobeniusNorm( wSVM );
        if( El::mpi::Rank() == 0 )
        {
            El::Output("SVM time: ",timer.Total()," secs");
            El::Output
            ("|| wSVM ||_2=",wSVMNorm,"\n",
             "margin      =",Real(2)/wSVMNorm,"\n",
             "offsetSVM   =",offsetSVM,"\n",
             "offsetSVM / || wSVM ||_2=",offsetSVM/wSVMNorm);
        }
        if( print )
            El::Print( wSVM, "wSVM" );

        // Report the classification percentage using the returned hyperplane
        El::DistMatrix<Real> qSVM;
        El::Ones( qSVM, m, 1 );
        El::Gemv( El::NORMAL, Real(1), G, wSVM, -offsetSVM, qSVM );
        El::EntrywiseMap( qSVM, El::MakeFunction(sgnMap) );
        if( print )
            El::Print( qSVM, "qSVM" );
        qSVM -= q;
        const Real numWrong = El::OneNorm(qSVM) / Real(2);
        if( El::mpi::Rank(comm) == 0 )
            El::Output("ratio misclassified: ",numWrong,"/",m);
    }
    catch( std::exception& e ) { El::ReportException(e); }

    return 0;
}
