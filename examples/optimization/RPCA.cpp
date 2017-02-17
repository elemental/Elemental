/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

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
template<typename Field>
int Corrupt( El::DistMatrix<Field>& A, double probCorrupt )
{
    EL_DEBUG_ONLY(El::CallStackEntry cse("Corrupt"))
    typedef El::Base<Field> Real;

    El::Int numLocalCorrupt = 0;
    const El::Int localHeight = A.LocalHeight();
    const El::Int localWidth = A.LocalWidth();
    for( El::Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        for( El::Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            if( El::SampleUniform<Real>() <= probCorrupt )
            {
                ++numLocalCorrupt;
                const Field perturb = El::SampleBall<Field>();
                A.SetLocal( iLoc, jLoc, A.GetLocal(iLoc,jLoc)+perturb );
            }
        }
    }

    return El::mpi::AllReduce( numLocalCorrupt, A.DistComm() );
}

int
main( int argc, char* argv[] )
{
    El::Environment env( argc, argv );
    const El::Int commRank = El::mpi::Rank();
    typedef El::Complex<double> Scalar;

    try
    {
        const El::Int m = El::Input("--height","height of matrix",100);
        const El::Int n = El::Input("--width","width of matrix",100);
        const El::Int rank = El::Input("--rank","rank of structured matrix",10);
        const double probCorrupt =
          El::Input("--probCorrupt","probability of corruption",0.1);
        const double tau = El::Input("--tau","sparse weighting factor",0.);
        const double beta = El::Input("--beta","step size",1.);
        const double rho = El::Input("--rho","stepsize multiple in ALM",6.);
        const El::Int maxIts = El::Input("--maxIts","maximum iterations",1000);
        const double tol = El::Input("--tol","tolerance",1.e-5);
        const bool usePivQR =
          El::Input("--usePivQR","use pivoted QR approx?",false);
        const El::Int numPivSteps =
          El::Input("--numPivSteps","number of steps of QR",75);
        const bool useALM = El::Input("--useALM","use ALM algorithm?",true);
        const bool display = El::Input("--display","display matrices",false);
        const bool print = El::Input("--print","print matrices",true);
        El::ProcessInput();
        El::PrintInputReport();

        El::DistMatrix<Scalar> LTrue;
        {
            El::DistMatrix<Scalar> U, V;
            El::Uniform( U, m, rank );
            El::Uniform( V, n, rank );
            El::Zeros( LTrue, m, n );
            // Since each entry of U and V is lies in the unit ball, every entry
            // of U V' will lie in the ball of radius 'rank', so scale this ball
            El::Gemm
            ( El::NORMAL, El::ADJOINT,
              Scalar(1./rank), U, V, Scalar(0), LTrue );
        }
        const double frobLTrue = El::FrobeniusNorm( LTrue );
        const double maxLTrue = El::MaxNorm( LTrue );
        if( commRank == 0 )
        {
            El::Output("|| L ||_F = ",frobLTrue);
            El::Output("|| L ||_max = ",maxLTrue);
        }
        if( display )
            El::Display( LTrue, "True low-rank" );
        if( print )
            El::Print( LTrue, "True low-rank" );

        El::DistMatrix<Scalar> STrue;
        El::Zeros( STrue, m, n );
        const El::Int numCorrupt = Corrupt( STrue, probCorrupt );
        const double frobSTrue = El::FrobeniusNorm( STrue );
        const double maxSTrue = El::MaxNorm( STrue );
        if( commRank == 0 )
        {
            El::Output("# of corrupted entries: ",numCorrupt);
            El::Output("|| S ||_F = ",frobSTrue);
            El::Output("|| S ||_max = ",maxSTrue);
        }
        if( display )
        {
            El::Display( STrue, "True sparse matrix" );
#ifdef EL_HAVE_QT5
            El::Spy( STrue, "True sparse spy plot" );
#endif
        }
        if( print )
            El::Print( STrue, "True sparse" );

        if( commRank == 0 )
            El::Output
            ("Using ",STrue.Grid().Height()," x ",STrue.Grid().Width(),
             " grid and blocksize of ",El::Blocksize());

        // M = LTrue + STrue
        El::DistMatrix<Scalar> M( LTrue );
        M += STrue;
        if( display )
            El::Display( M, "Sum of low-rank and sparse");
        if( print )
            El::Print( M, "Sum of low-rank and sparse" );

        // Create a custom set of parameters for RPCA
        El::RPCACtrl<double> ctrl;
        ctrl.useALM = useALM;
        ctrl.usePivQR = usePivQR;
        ctrl.progress = print;
        ctrl.numPivSteps = numPivSteps;
        ctrl.maxIts = maxIts;
        ctrl.tau = tau;
        ctrl.beta = beta;
        ctrl.rho = rho;
        ctrl.tol = tol;

        El::Timer timer;
        El::DistMatrix<Scalar> L, S;
        if( commRank == 0 )
            timer.Start();
        El::RPCA( M, L, S, ctrl );
        if( commRank == 0 )
            timer.Stop();

        if( display )
        {
            El::Display( L, "Estimated low-rank matrix" );
            El::Display( S, "Estimated sparse matrix" );
#ifdef EL_HAVE_QT5
            El::Spy( S, "Estimated sparse spy plot" );
#endif
        }
        if( print )
        {
            El::Print( L, "Estimated low-rank matrix" );
            El::Print( S, "Estimated sparse matrix" );
        }

        L -= LTrue;
        S -= STrue;
        const double frobLDiff = El::FrobeniusNorm( L );
        const double frobSDiff = El::FrobeniusNorm( S );
        if( commRank == 0 )
        {
            El::Output("RPCA time: ",timer.Total()," secs");
            El::Output
            ("Error in decomposition:\n",
             "  || L - LTrue ||_F / || LTrue ||_F = ",frobLDiff/frobLTrue,"\n",
             "  || S - STrue ||_F / || STrue ||_F = ",frobSDiff/frobSTrue,"\n");
        }

        if( display )
        {
            El::Display( L, "Error in low-rank estimate" );
            El::Display( S, "Error in sparse estimate" );
        }
        if( print )
        {
            El::Print( L, "Error in low-rank estimate" );
            El::Print( S, "Error in sparse estimate" );
        }
    }
    catch( std::exception& e ) { El::ReportException(e); }

    return 0;
}
