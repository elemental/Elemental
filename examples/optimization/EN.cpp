#include "El.hpp"
using namespace El;

DistSparseMatrix<double> ConcatFD2D( Int n0, Int n1 )
{
    DistSparseMatrix<double> A;
    const Int height = n0*n1;
    const Int width = 2*n0*n1;
    A.Resize( height, width );
    const Int localHeight = A.LocalHeight();
    A.Reserve( 7*localHeight );

    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = A.GlobalRow(iLoc);
        const Int x0 = i % n0;
        const Int x1 = i / n0;

        A.QueueLocalUpdate( iLoc, i, 15 );
        if( x0 > 0 )    A.QueueLocalUpdate( iLoc, i-1,  -1 );
        if( x0+1 < n0 ) A.QueueLocalUpdate( iLoc, i+1,   2 );
        if( x1 > 0 )    A.QueueLocalUpdate( iLoc, i-n0, -3 );
        if( x1+1 < n1 ) A.QueueLocalUpdate( iLoc, i+n0,  4 );

        const Int iRel = i + n0*n1;
        A.QueueLocalUpdate( iLoc, iRel, 1 );

        // For now, this is meant to use integer division (for reproducing)
        A.QueueLocalUpdate( iLoc, width-1, -10/height );
    }
    A.ProcessQueues();
    return A;
}

int main( int argc, char* argv[] )
{
    Environment env( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::Rank( comm );

    try
    {
        const Int n0 = Input("--n0","height of 2D discretization",50);
        const Int n1 = Input("--n1","width of 2D discretization",50);
        const double lambda1 = Input("--lambda1","one-norm coefficient",3.);
        const double lambda2 = Input("--lambda2","two-norm coefficient",4.);
        const bool print = Input("--print","print matrices?",false);
        const bool prog = Input("--prog","print progress info?",true);
        const bool time = Input("--time","time each step of IPM?",true);
        const bool solveProg = Input("--solveProg","solver progress?",true);
        const bool solveTime = Input("--solveTime","solver timers?",true);
        ProcessInput();
        PrintInputReport();

        auto A = ConcatFD2D( n0, n1 );
        DistMultiVec<double> b;
        Gaussian( b, n0*n1, 1 );
        if( print )
        {
            Print( A, "A" );
            Print( b, "b" );
        }
        
        qp::affine::Ctrl<double> ctrl; 
        ctrl.mehrotraCtrl.print = prog;
        ctrl.mehrotraCtrl.time = time;
        ctrl.mehrotraCtrl.solveCtrl.progress = solveProg;
        ctrl.mehrotraCtrl.solveCtrl.time = solveTime;
        
        DistMultiVec<double> x;
        Timer timer;
        mpi::Barrier( comm );
        if( commRank == 0 )
            timer.Start();
        EN( A, b, lambda1, lambda2, x, ctrl ); 
        if( commRank == 0 )
            Output("EN time: ",timer.Stop()," secs");
        if( print ) 
            Print( x, "x" );
    }
    catch( const exception& e ) { ReportException(e); }

    return 0;
}
