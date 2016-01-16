/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

//
// This driver generates a random low-rank matrix and then randomly corrupts
// a large percentage of the entries (by default 10%). Robust Principal 
// Component Analysis (RPCA) is then used to recover both the underlying 
// low-rank and sparse matrices.
//
// Please see <http://perception.csl.illinois.edu/matrix-rank/sample_code.html>
// for references and documentation on the Augmented Lagrange Multiplier (ALM) 
// and Alternating Direction Method of Multipliers (ADMM) for Robust PCA.
//

// Corrupt a portion of the entries with uniform samples from the unit ball
template<typename F>
int Corrupt( DistMatrix<F>& A, double probCorrupt )
{
    DEBUG_ONLY(CallStackEntry cse("Corrupt"))
    typedef Base<F> Real;

    Int numLocalCorrupt = 0;
    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            if( SampleUniform<Real>() <= probCorrupt )
            {
                ++numLocalCorrupt;
                const F perturb = SampleBall<F>();
                A.SetLocal( iLoc, jLoc, A.GetLocal(iLoc,jLoc)+perturb );
            }
        }
    }
    
    return mpi::AllReduce( numLocalCorrupt, A.DistComm() );
}

int 
main( int argc, char* argv[] )
{
    Environment env( argc, argv );
    const Int commRank = mpi::Rank();
    typedef Complex<double> C;

    try
    {
        const Int m = Input("--height","height of matrix",100);
        const Int n = Input("--width","width of matrix",100);
        const Int rank = Input("--rank","rank of structured matrix",10);
        const double probCorrupt = 
            Input("--probCorrupt","probability of corruption",0.1);
        const double tau = Input("--tau","sparse weighting factor",0.);
        const double beta = Input("--beta","step size",1.);
        const double rho = Input("--rho","stepsize multiple in ALM",6.);
        const Int maxIts = Input("--maxIts","maximum iterations",1000);
        const double tol = Input("--tol","tolerance",1.e-5);
        const bool usePivQR = 
            Input("--usePivQR","use pivoted QR approx?",false);
        const Int numPivSteps = 
            Input("--numPivSteps","number of steps of QR",75);
        const bool useALM = Input("--useALM","use ALM algorithm?",true);
        const bool display = Input("--display","display matrices",true);
        const bool print = Input("--print","print matrices",false);
        ProcessInput();
        PrintInputReport();

        DistMatrix<C> LTrue;
        {
            DistMatrix<C> U, V;
            Uniform( U, m, rank );
            Uniform( V, n, rank );
            Zeros( LTrue, m, n );
            // Since each entry of U and V is lies in the unit ball, every entry
            // of U V' will lie in the ball of radius 'rank', so scale this ball
            Gemm( NORMAL, ADJOINT, C(1./rank), U, V, C(0), LTrue );
        }
        const double frobLTrue = FrobeniusNorm( LTrue );
        const double maxLTrue = MaxNorm( LTrue );
        if( commRank == 0 )
        {
            Output("|| L ||_F = ",frobLTrue);
            Output("|| L ||_max = ",maxLTrue);
        }
        if( display )
            Display( LTrue, "True low-rank" );
        if( print )
            Print( LTrue, "True low-rank" );

        DistMatrix<C> STrue;
        Zeros( STrue, m, n );
        const Int numCorrupt = Corrupt( STrue, probCorrupt );
        const double frobSTrue = FrobeniusNorm( STrue );
        const double maxSTrue = MaxNorm( STrue );
        if( commRank == 0 )
        {
            Output("# of corrupted entries: ",numCorrupt);
            Output("|| S ||_F = ",frobSTrue);
            Output("|| S ||_max = ",maxSTrue);
        }
        if( display )
        {
            Display( STrue, "True sparse matrix" );
#ifdef EL_HAVE_QT5
            Spy( STrue, "True sparse spy plot" );
#endif
        }
        if( print )
            Print( STrue, "True sparse" );

        if( commRank == 0 )
            Output
            ("Using ",STrue.Grid().Height()," x ",STrue.Grid().Width(),
             " grid and blocksize of ",Blocksize());

        // M = LTrue + STrue
        DistMatrix<C> M( LTrue );
        M += STrue;
        if( display )
            Display( M, "Sum of low-rank and sparse");
        if( print )
            Print( M, "Sum of low-rank and sparse" );

        // Create a custom set of parameters for RPCA
        RPCACtrl<double> ctrl;
        ctrl.useALM = useALM;
        ctrl.usePivQR = usePivQR;
        ctrl.progress = print;
        ctrl.numPivSteps = numPivSteps;
        ctrl.maxIts = maxIts;
        ctrl.tau = tau;
        ctrl.beta = beta;
        ctrl.rho = rho;
        ctrl.tol = tol;

        Timer timer;
        DistMatrix<C> L, S;
        if( commRank == 0 )
            timer.Start();
        RPCA( M, L, S, ctrl );
        if( commRank == 0 )
            timer.Stop();

        if( display )
        {
            Display( L, "Estimated low-rank matrix" );
            Display( S, "Estimated sparse matrix" );
#ifdef EL_HAVE_QT5
            Spy( S, "Estimated sparse spy plot" );
#endif
        }
        if( print )
        {
            Print( L, "Estimated low-rank matrix" );
            Print( S, "Estimated sparse matrix" );
        }

        L -= LTrue;
        S -= STrue;
        const double frobLDiff = FrobeniusNorm( L );
        const double frobSDiff = FrobeniusNorm( S );
        if( commRank == 0 )
        {
            Output("RPCA time: ",timer.Total()," secs");
            Output
            ("Error in decomposition:\n",
             "  || L - LTrue ||_F / || LTrue ||_F = ",frobLDiff/frobLTrue,"\n",
             "  || S - STrue ||_F / || STrue ||_F = ",frobSDiff/frobSTrue,"\n");
        }

        if( display )
        {
            Display( L, "Error in low-rank estimate" );
            Display( S, "Error in sparse estimate" );
        }
        if( print )
        {
            Print( L, "Error in low-rank estimate" );
            Print( S, "Error in sparse estimate" );
        }
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
