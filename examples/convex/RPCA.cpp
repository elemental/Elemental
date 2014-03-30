/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include ELEM_AXPY_INC
#include ELEM_SCALE_INC
#include ELEM_FROBENIUSNORM_INC
#include ELEM_ENTRYWISEONENORM_INC
#include ELEM_MAXNORM_INC
#include ELEM_TWONORM_INC
#include ELEM_ZERONORM_INC
#include ELEM_SVT_INC
#include ELEM_UNIFORM_INC
using namespace elem;

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

template<typename F,Distribution U,Distribution V>
void NormalizeEntries( DistMatrix<F,U,V>& A )
{
    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const F alpha = A.GetLocal( iLoc, jLoc );
            A.SetLocal( iLoc, jLoc, alpha/Abs(alpha) );
        }
    }
}

// If 'tau' is passed in as zero, it is set to 1/sqrt(max(m,n))
template<typename F>
void RPCA_ADMM
( const DistMatrix<F>& M, DistMatrix<F>& L, DistMatrix<F>& S, 
  Base<F> beta, 
  Base<F> tau, 
  Base<F> tol, 
  Int numStepsQR,
  Int maxIts,
  bool print )
{
    typedef Base<F> Real;
    const Int m = M.Height();
    const Int n = M.Width();
    const Int commRank = mpi::Rank( M.Grid().Comm() );

    // If tau is not specified, then set it to 1/sqrt(max(m,n))
    if( tau == Real(0) )
        tau = Real(1)/sqrt(Real(std::max(m,n)));

    if( beta <= Real(0) )
        LogicError("beta cannot be non-positive");
    if( tau <=  Real(0) )
        LogicError("tau cannot be non-positive");
    if( tol <= Real(0) )
        LogicError("tol cannot be non-positive");

    const double startTime = mpi::Time();
    if( commRank == 0 )
        std::cout << "Starting RPCA_ADMM" << std::endl;

    DistMatrix<F> E( M.Grid() ), Y( M.Grid() );
    Zeros( Y, m, n );

    const Real frobM = FrobeniusNorm( M );
    const Real maxM = MaxNorm( M );
    if( commRank == 0 )
        std::cout << "|| M ||_F = " << frobM << "\n"
                  << "|| M ||_max = " << maxM << std::endl;

    Int numIts = 0;
    while( true )
    {
        ++numIts;

        // ST_{tau/beta}(M - L + Y/beta)
        S = M;
        Axpy( F(-1), L, S );
        Axpy( F(1)/beta, Y, S );
        SoftThreshold( S, tau/beta );
        const Int numNonzeros = ZeroNorm( S );

        // SVT_{1/beta}(M - S + Y/beta)
        L = M;
        Axpy( F(-1), S, L );
        Axpy( F(1)/beta, Y, L );
        Int rank;
        if( numStepsQR == -1 )
            rank = SVT( L, Real(1)/beta );
        else
            rank = SVT( L, Real(1)/beta, numStepsQR );
      
        // E := M - (L + S)
        E = M;    
        Axpy( F(-1), L, E );
        Axpy( F(-1), S, E );
        const Real frobE = FrobeniusNorm( E );

        if( frobE/frobM <= tol )            
        {
            if( commRank == 0 )
                std::cout << "Converged after " << numIts << " iterations "
                          << " with rank=" << rank 
                          << ", numNonzeros=" << numNonzeros << " and "
                          << "|| E ||_F / || M ||_F = " << frobE/frobM
                          << ", and " << mpi::Time()-startTime << " total secs"
                          << std::endl;
            break;
        }
        else if( numIts >= maxIts )
        {
            if( commRank == 0 )
                std::cout << "Aborting after " << maxIts << " iterations and "
                          << mpi::Time()-startTime << " total secs" 
                          << std::endl;
            break;
        }
        else
        {
            if( commRank == 0 )
                std::cout << numIts << ": || E ||_F / || M ||_F = " 
                          << frobE/frobM << ", rank=" << rank 
                          << ", numNonzeros=" << numNonzeros 
                          << ", " << mpi::Time()-startTime << " total secs"
                          << std::endl;
        }
        
        // Y := Y + beta E
        Axpy( beta, E, Y );
    }
}

// If 'beta' or 'tau' is passed in as zero, then an estimate is used instead
template<typename F>
void RPCA_ALM
( const DistMatrix<F>& M, DistMatrix<F>& L, DistMatrix<F>& S, 
  Base<F> beta, Base<F> tau, Base<F> rho, Base<F> tol, 
  Int numStepsQR, Int maxIts, bool print )
{
    typedef Base<F> Real;

    const Int m = M.Height();
    const Int n = M.Width();
    const Int commRank = mpi::Rank( M.Grid().Comm() );

    // If tau is unspecified, set it to 1/sqrt(max(m,n))
    if( tau == Real(0) )
        tau = Real(1) / sqrt(Real(std::max(m,n)));

    if( tol <= Real(0) )
        LogicError("tol cannot be non-positive");
    if( tau <= Real(0) )
        LogicError("tau cannot be non-positive");

    const double startTime = mpi::Time();
    if( commRank == 0 )
        std::cout << "Starting RPCA_ALM" << std::endl;

    DistMatrix<F> Y( M );
    NormalizeEntries( Y );
    const Real twoNorm = TwoNorm( Y );
    const Real maxNorm = MaxNorm( Y );
    const Real infNorm = maxNorm / tau; 
    const Real dualNorm = std::max( twoNorm, infNorm );
    Scale( F(1)/dualNorm, Y );

    // If beta is unspecified, set it to 1 / 2 || sign(M) ||_2
    if( beta == Real(0) )
        beta = Real(1) / (2*twoNorm);

    if( beta <= Real(0) )
        LogicError("beta cannot be non-positive");

    const Real frobM = FrobeniusNorm( M );
    const Real maxM = MaxNorm( M );
    if( commRank == 0 )
        std::cout << "|| M ||_F = " << frobM << "\n"
                  << "|| M ||_max = " << maxM << std::endl;

    Int numIts=0, numPrimalIts=0;
    DistMatrix<F> LLast( M.Grid() ), SLast( M.Grid() ), E( M.Grid() );
    while( true )
    {
        ++numIts;
       
        Int rank, numNonzeros;
        while( true )
        {
            ++numPrimalIts;

            LLast = L;
            SLast = S;

            // ST_{tau/beta}(M - L + Y/beta)
            S = M;
            Axpy( F(-1), L, S );
            Axpy( F(1)/beta, Y, S );
            SoftThreshold( S, tau/beta );
            numNonzeros = ZeroNorm( S );

            // SVT_{1/beta}(M - S + Y/beta)
            L = M;
            Axpy( F(-1), S, L );
            Axpy( F(1)/beta, Y, L );
            if( commRank == 0 )
                std::cout << "beta=" << beta << std::endl;
            if( numStepsQR == -1 )
                rank = SVT( L, Real(1)/beta );
            else
                rank = SVT( L, Real(1)/beta, numStepsQR );

            Axpy( F(-1), L, LLast );
            Axpy( F(-1), S, SLast );
            const Real frobLDiff = FrobeniusNorm( LLast );
            const Real frobSDiff = FrobeniusNorm( SLast );

            if( frobLDiff/frobM < tol && frobSDiff/frobM < tol )
            {
                if( commRank == 0 )
                    std::cout << "Primal loop converged: " 
                              << mpi::Time()-startTime << " total secs"
                              << std::endl;
                break;
            }
            else 
            {
                if( commRank == 0 )
                    std::cout << "  " << numPrimalIts 
                              << ": \n"
                              << "   || Delta L ||_F / || M ||_F = " 
                              << frobLDiff/frobM << "\n"
                              << "   || Delta S ||_F / || M ||_F = "
                              << frobSDiff/frobM << "\n"
                              << "   rank=" << rank
                              << ", numNonzeros=" << numNonzeros 
                              << ", " << mpi::Time()-startTime << " total secs" 
                              << std::endl;
            } 
        }

        // E := M - (L + S)
        E = M;    
        Axpy( -1., L, E );
        Axpy( -1., S, E );
        const Real frobE = FrobeniusNorm( E );

        if( frobE/frobM <= tol )            
        {
            if( commRank == 0 )
                std::cout << "Converged after " << numIts << " iterations and "
                          << numPrimalIts << " primal iterations with rank=" 
                          << rank << ", numNonzeros=" << numNonzeros << " and "
                          << "|| E ||_F / || M ||_F = " << frobE/frobM
                          << ", " << mpi::Time()-startTime << " total secs"
                          << std::endl;
            break;
        }
        else if( numIts >= maxIts )
        {
            if( commRank == 0 )
                std::cout << "Aborting after " << maxIts << " iterations and "
                          << mpi::Time()-startTime << " total secs" 
                          << std::endl;
            break;
        }
        else
        {
            if( commRank == 0 )
                std::cout << numPrimalIts << ": || E ||_F / || M ||_F = " 
                          << frobE/frobM << ", rank=" << rank 
                          << ", numNonzeros=" << numNonzeros << ", "
                          << mpi::Time()-startTime << " total secs" 
                          << std::endl;
        }
        
        // Y := Y + beta E
        Axpy( beta, E, Y );
        beta *= rho;
    }
}

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const Int commRank = mpi::Rank( comm );
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
        const Int numStepsQR = Input("--numStepsQR","number of steps of QR",-1);
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
            std::cout << "|| L ||_F = " << frobLTrue << "\n"
                      << "|| L ||_max = " << maxLTrue << std::endl;
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
            std::cout << "number of corrupted entries: " << numCorrupt << "\n"
                      << "|| S ||_F = " << frobSTrue << "\n"
                      << "|| S ||_max = " << maxSTrue << std::endl;
        if( display )
        {
            Display( STrue, "True sparse matrix" );
#ifdef ELEM_HAVE_QT5
            Spy( STrue, "True sparse spy plot" );
#endif
        }
        if( print )
            Print( STrue, "True sparse" );

        if( commRank == 0 )
            std::cout << "Using " << STrue.Grid().Height() << " x " 
                      << STrue.Grid().Width() 
                      << " process grid and blocksize of " << Blocksize() 
                      << std::endl;

        // M = LTrue + STrue
        DistMatrix<C> M( LTrue );
        Axpy( C(1), STrue, M );
        if( display )
            Display( M, "Sum of low-rank and sparse");
        if( print )
            Print( M, "Sum of low-rank and sparse" );

        DistMatrix<C> L, S;
        Zeros( L, m, n );
        Zeros( S, m, n ); 

        if( useALM )
            RPCA_ALM( M, L, S, beta, tau, rho, tol, numStepsQR, maxIts, print );
        else
            RPCA_ADMM( M, L, S, beta, tau, tol, numStepsQR, maxIts, print );

        if( display )
        {
            Display( L, "Estimated low-rank matrix" );
            Display( S, "Estimated sparse matrix" );
#ifdef ELEM_HAVE_QT5
            Spy( S, "Estimated sparse spy plot" );
#endif
        }
        if( print )
        {
            Print( L, "Estimated low-rank matrix" );
            Print( S, "Estimated sparse matrix" );
        }

        Axpy( C(-1), LTrue, L );
        Axpy( C(-1), STrue, S );
        const double frobLDiff = FrobeniusNorm( L );
        const double frobSDiff = FrobeniusNorm( S );
        if( commRank == 0 )
            std::cout << "\n"
                      << "Error in computed decomposition:\n"
                      << "  || L - LTrue ||_F / || LTrue ||_F = " 
                      << frobLDiff/frobLTrue << "\n"
                      << "  || S - STrue ||_F / || STrue ||_F = " 
                      << frobSDiff/frobSTrue << "\n"
                      << std::endl;

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
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
