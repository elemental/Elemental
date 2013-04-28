/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include "elemental/blas-like/level1/Axpy.hpp"
#include "elemental/blas-like/level1/Scale.hpp"
#include "elemental/lapack-like/Norm/Frobenius.hpp"
#include "elemental/lapack-like/Norm/EntrywiseOne.hpp"
#include "elemental/lapack-like/Norm/Max.hpp"
#include "elemental/lapack-like/Norm/Two.hpp"
#include "elemental/lapack-like/Norm/Zero.hpp"
#include "elemental/convex/SingularValueSoftThreshold.hpp"
#include "elemental/matrices/Uniform.hpp"
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
#ifndef RELEASE
    CallStackEntry entry("Corrupt");
#endif
    typedef BASE(F) R;

    int numLocalCorrupt = 0;
    const int localHeight = A.LocalHeight();
    const int localWidth = A.LocalWidth();
    for( int jLocal=0; jLocal<localWidth; ++jLocal )
    {
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
        {
            if( Abs(SampleUnitBall<R>()) <= probCorrupt )
            {
                ++numLocalCorrupt;
                const F perturb = SampleUnitBall<F>();
                A.SetLocal( iLocal, jLocal, A.GetLocal(iLocal,jLocal)+perturb );
            }
        }
    }
    
    int numCorrupt;
    mpi::AllReduce
    ( &numLocalCorrupt, &numCorrupt, 1, mpi::SUM, A.Grid().VCComm() );
    return numCorrupt;
}

template<typename F,Distribution U,Distribution V>
void Sign( DistMatrix<F,U,V>& A )
{
    const int localHeight = A.LocalHeight();
    const int localWidth = A.LocalWidth();
    for( int jLocal=0; jLocal<localWidth; ++jLocal )
    {
        for( int iLocal=0; iLocal<localHeight; ++iLocal )
        {
            const F alpha = A.GetLocal( iLocal, jLocal );
            A.SetLocal( iLocal, jLocal, alpha/Abs(alpha) );
        }
    }
}

// If 'tau' is passed in as zero, it is set to 1/sqrt(max(m,n))
template<typename F>
void RPCA_ADMM
( const DistMatrix<F>& M, DistMatrix<F>& L, DistMatrix<F>& S, 
  BASE(F) beta, 
  BASE(F) tau, 
  BASE(F) tol, 
  int numStepsQR,
  int maxIts,
  bool print )
{
    typedef BASE(F) R;
    const int m = M.Height();
    const int n = M.Width();
    const int commRank = mpi::CommRank( M.Grid().Comm() );

    // If tau is not specified, then set it to 1/sqrt(max(m,n))
    if( tau == R(0) )
        tau = R(1)/sqrt(R(std::max(m,n)));

    if( beta <= R(0) )
        throw std::logic_error("beta cannot be non-positive");
    if( tau <=  R(0) )
        throw std::logic_error("tau cannot be non-positive");
    if( tol <= R(0) )
        throw std::logic_error("tol cannot be non-positive");

    const double startTime = mpi::Time();
    if( commRank == 0 )
        std::cout << "Starting RPCA_ADMM after " << startTime << " secs" 
                  << std::endl;

    DistMatrix<F> E( M.Grid() ), Y( M.Grid() );
    Zeros( Y, m, n );

    const R frobM = FrobeniusNorm( M );
    if( commRank == 0 )
        std::cout << "|| M ||_F = " << frobM << std::endl;

    int numIts = 0;
    while( true )
    {
        ++numIts;

        // ST_{tau/beta}(M - L + Y/beta)
        S = M;
        Axpy( F(-1), L, S );
        Axpy( F(1)/beta, Y, S );
        SoftThreshold( S, tau/beta );
        const int numNonzeros = ZeroNorm( S );

        // SVT_{1/beta}(M - S + Y/beta)
        L = M;
        Axpy( F(-1), S, L );
        Axpy( F(1)/beta, Y, L );
        int rank;
        if( numStepsQR == -1 )
            rank = SingularValueSoftThreshold( L, R(1)/beta );
        else
            rank = SingularValueSoftThreshold( L, R(1)/beta, numStepsQR );
      
        // E := M - (L + S)
        E = M;    
        Axpy( F(-1), L, E );
        Axpy( F(-1), S, E );
        const R frobE = FrobeniusNorm( E );

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
  BASE(F) beta, BASE(F) tau, BASE(F) rho, BASE(F) tol, 
  int numStepsQR, int maxIts, bool print )
{
    typedef BASE(F) R;

    const int m = M.Height();
    const int n = M.Width();
    const int commRank = mpi::CommRank( M.Grid().Comm() );

    // If tau is unspecified, set it to 1/sqrt(max(m,n))
    if( tau == R(0) )
        tau = R(1) / sqrt(R(std::max(m,n)));

    if( tol <= R(0) )
        throw std::logic_error("tol cannot be non-positive");
    if( tau <= R(0) )
        throw std::logic_error("tau cannot be non-positive");

    const double startTime = mpi::Time();
    if( commRank == 0 )
        std::cout << "Starting RPCA_ALM after " << startTime << " secs" 
                  << std::endl;

    DistMatrix<F> Y( M );
    Sign( Y );
    const R twoNorm = TwoNorm( Y );
    const R maxNorm = MaxNorm( Y );
    const R infNorm = maxNorm / tau; 
    const R dualNorm = std::max( twoNorm, infNorm );
    Scale( F(1)/dualNorm, Y );

    // If beta is unspecified, set it to 1 / 2 || sign(M) ||_2
    if( beta == R(0) )
        beta = R(1) / (2*twoNorm);

    if( beta <= R(0) )
        throw std::logic_error("beta cannot be non-positive");

    const R frobM = FrobeniusNorm( M );
    if( commRank == 0 )
        std::cout << "|| M ||_F = " << frobM << std::endl;

    int numIts=0, numPrimalIts=0;
    DistMatrix<F> LLast( M.Grid() ), SLast( M.Grid() ), E( M.Grid() );
    while( true )
    {
        ++numIts;
       
        int rank, numNonzeros;
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
                rank = SingularValueSoftThreshold( L, R(1)/beta );
            else
                rank = SingularValueSoftThreshold( L, R(1)/beta, numStepsQR );

            Axpy( F(-1), L, LLast );
            Axpy( F(-1), S, SLast );
            const R frobLDiff = FrobeniusNorm( LLast );
            const R frobSDiff = FrobeniusNorm( SLast );

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
        const R frobE = FrobeniusNorm( E );

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
    const int commRank = mpi::CommRank( comm );

    typedef Complex<double> C;

    try
    {
        const int m = Input("--height","height of matrix",100);
        const int n = Input("--width","width of matrix",100);
        const int rank = Input("--rank","rank of structured matrix",10);
        const double probCorrupt = 
            Input("--probCorrupt","probability of corruption",0.1);
        const double tau = Input("--tau","sparse weighting factor",0.);
        const double beta = Input("--beta","step size",1.);
        const double rho = Input("--rho","stepsize multiple in ALM",6.);
        const int maxIts = Input("--maxIts","maximum iterations",1000);
        const double tol = Input("--tol","tolerance",1.e-6);
        const int numStepsQR = Input("--numStepsQR","number of steps of QR",-1);
        const bool useALM = Input("--useALM","use ALM algorithm?",true);
        const bool print = Input("--print","print matrices",false);
        ProcessInput();
        PrintInputReport();

        DistMatrix<C> LTrue;
        {
            DistMatrix<C> U, V;
            Uniform( U, m, rank );
            Uniform( V, n, rank );
            Zeros( LTrue, m, n );
            Gemm( NORMAL, ADJOINT, C(1./std::max(m,n)), U, V, C(0), LTrue );
        }
        const double frobLTrue = FrobeniusNorm( LTrue );
        if( commRank == 0 )
            std::cout << "|| L ||_F = " << frobLTrue << std::endl;
        if( print )
            LTrue.Print("True L");

        DistMatrix<C> STrue;
        Zeros( STrue, m, n );
        const int numCorrupt = Corrupt( STrue, probCorrupt );
        const double frobSTrue = FrobeniusNorm( STrue );
        if( commRank == 0 )
            std::cout << "number of corrupted entries: " << numCorrupt << "\n"
                      << "|| S ||_F = " << frobSTrue << std::endl;
        if( print )
            STrue.Print("True S");

        if( commRank == 0 )
            std::cout << "Using " << STrue.Grid().Height() << " x " 
                      << STrue.Grid().Width() 
                      << " process grid and blocksize of " << Blocksize() 
                      << std::endl;

        // M = LTrue + STrue
        DistMatrix<C> M( LTrue );
        Axpy( C(1), STrue, M );

        DistMatrix<C> L, S;
        Zeros( L, m, n );
        Zeros( S, m, n ); 

        if( useALM )
            RPCA_ALM( M, L, S, beta, tau, rho, tol, numStepsQR, maxIts, print );
        else
            RPCA_ADMM( M, L, S, beta, tau, tol, numStepsQR, maxIts, print );

        if( print )
        {
            L.Print("L");
            S.Print("S"); 
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
        if( print )
        {
            L.Print("L - LTrue");
            S.Print("S - STrue");
        }
    }
    catch( ArgException& e )
    {
        // There is nothing to do
    }
    catch( std::exception& e )
    {
        std::ostringstream os;
        os << "Process " << commRank << " caught error message:\n" << e.what()
           << std::endl;
        std::cerr << os.str();
#ifndef RELEASE
        DumpCallStack();
#endif
    }

    Finalize();
    return 0;
}
