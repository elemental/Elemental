/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

El::DistSparseMatrix<double> ConcatFD2D( El::Int n0, El::Int n1 )
{
    El::DistSparseMatrix<double> A;
    const El::Int height = n0*n1;
    const El::Int width = 2*n0*n1;
    A.Resize( height, width );
    const El::Int localHeight = A.LocalHeight();
    A.Reserve( 7*localHeight );

    for( El::Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const El::Int i = A.GlobalRow(iLoc);
        const El::Int x0 = i % n0;
        const El::Int x1 = i / n0;

        A.QueueLocalUpdate( iLoc, i, 15 );
        if( x0 > 0 )    A.QueueLocalUpdate( iLoc, i-1,  -1 );
        if( x0+1 < n0 ) A.QueueLocalUpdate( iLoc, i+1,   2 );
        if( x1 > 0 )    A.QueueLocalUpdate( iLoc, i-n0, -3 );
        if( x1+1 < n1 ) A.QueueLocalUpdate( iLoc, i+n0,  4 );

        const El::Int iRel = i + n0*n1;
        A.QueueLocalUpdate( iLoc, iRel, 1 );

        // For now, this is meant to use integer division (for reproducing)
        A.QueueLocalUpdate( iLoc, width-1, -10/height );
    }
    A.ProcessQueues();
    return A;
}

int main( int argc, char* argv[] )
{
    El::Environment env( argc, argv );
    El::mpi::Comm comm = El::mpi::COMM_WORLD;
    const int commRank = El::mpi::Rank( comm );

    try
    {
        const El::Int n0 = El::Input("--n0","height of 2D discretization",50);
        const El::Int n1 = El::Input("--n1","width of 2D discretization",50);
        const double lambda1 = El::Input("--lambda1","one-norm coefficient",3.);
        const double lambda2 = El::Input("--lambda2","two-norm coefficient",4.);
        const bool print = El::Input("--print","print matrices?",false);
        const bool prog = El::Input("--prog","print progress info?",true);
        const bool time = El::Input("--time","time each step of IPM?",true);
        const bool solveProg = El::Input("--solveProg","solver progress?",true);
        const bool solveTime = El::Input("--solveTime","solver timers?",true);
        El::ProcessInput();
        El::PrintInputReport();

        auto A = ConcatFD2D( n0, n1 );
        El::DistMultiVec<double> b;
        El::Gaussian( b, n0*n1, 1 );
        if( print )
        {
            El::Print( A, "A" );
            El::Print( b, "b" );
        }

        El::qp::affine::Ctrl<double> ctrl;
        ctrl.mehrotraCtrl.print = prog;
        ctrl.mehrotraCtrl.time = time;
        ctrl.mehrotraCtrl.solveCtrl.progress = solveProg;
        ctrl.mehrotraCtrl.solveCtrl.time = solveTime;

        El::DistMultiVec<double> x;
        El::Timer timer;
        El::mpi::Barrier( comm );
        if( commRank == 0 )
            timer.Start();
        El::EN( A, b, lambda1, lambda2, x, ctrl );
        if( commRank == 0 )
            El::Output("EN time: ",timer.Stop()," secs");
        if( print )
            El::Print( x, "x" );
    }
    catch( const std::exception& e ) { El::ReportException(e); }

    return 0;
}
