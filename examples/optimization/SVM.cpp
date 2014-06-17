/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "El.hpp" instead
#include "El-lite.hpp"
#include EL_GAUSSIAN_INC


using namespace El;

typedef double Real;

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );

    try
    {
        const Int m = Input("--numExamples","number of examples",200);
        const Int n = Input("--numFeatures","number of features",100);
        const Int maxIter = Input("--maxIter","maximum # of iter's",500);
        const Real gamma = Input("--gamma","two-norm coefficient",0.01);
        const Real rho = Input("--rho","augmented Lagrangian param.",1.);
        const Real alpha = Input("--alpha","over-relaxation",1.2);
        const Real absTol = Input("--absTol","absolute tolerance",1e-6);
        const Real relTol = Input("--relTol","relative tolerance",1e-4);
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
            Scale( Real(1)/wNorm, w );
        }
        Real offset = SampleNormal();
        mpi::Broadcast( offset, 0, mpi::COMM_WORLD );

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
        EntrywiseMap
        ( q, []( Real alpha ) 
             { return ( alpha >=0 ? Real(1) : Real(-1) ); } );

        if( mpi::WorldRank() == 0 )
            std::cout << "offset=" << offset << std::endl;
        if( print )
        {
            Print( w, "w" );
            Print( G, "G" );
            Print( q, "q" );
        }
        if( display )
            Display( G, "G" );

        DistMatrix<Real> wHatSVM;
        SVM
        ( G, q, gamma, wHatSVM, rho, alpha, maxIter, absTol, relTol, inv, 
          progress );
        auto wSVM = View( wHatSVM, 0, 0, n, 1 );
        const Real offsetSVM = -wHatSVM.Get(n,0);
        const Real wSVMNorm = FrobeniusNorm( wSVM );
        if( mpi::WorldRank() == 0 )
            std::cout << "|| wSVM ||_2=" << wSVMNorm << "\n"
                      << "margin      =" << Real(2)/wSVMNorm << "\n"
                      << "offsetSVM=" << offsetSVM << "\n"
                      << "offsetSVM / || wSVM ||_2=" << offsetSVM/wSVMNorm 
                      << std::endl;
        if( print )
            Print( wSVM, "wSVM" );

        // Report the classification percentage using the returned hyperplane
        DistMatrix<Real> qSVM;
        Ones( qSVM, m, 1 );
        Gemv( NORMAL, Real(1), G, wSVM, -offsetSVM, qSVM );
        EntrywiseMap
        ( qSVM, []( Real alpha ) 
                { return ( alpha >=0 ? Real(1) : Real(-1) ); } );
        if( print )
            Print( qSVM, "qSVM" );
        Axpy( Real(-1), q, qSVM );
        const Real numWrong = OneNorm(qSVM) / Real(2);
        if( mpi::WorldRank() == 0 )
            std::cout << "ratio misclassified: " << numWrong << "/" << m 
                      << std::endl;
    }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
