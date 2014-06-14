/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "El.hpp" instead
#include "El-lite.hpp"

//
// Please see <http://perception.csl.illinois.edu/matrix-rank/sample_code.html>
// for references and documentation on the Augmented Lagrange Multiplier (ALM) 
// and Alternating Direction Method of Multipliers (ADMM) for Robust PCA.
//

namespace El {

namespace rpca {

template<typename F,Distribution U,Distribution V>
inline void NormalizeEntries( DistMatrix<F,U,V>& A )
{ 
    EntrywiseMap
    ( A, []( F alpha ) { return alpha==F(0) ? F(1) : alpha/Abs(alpha); } ); 
}

// If 'tau' is passed in as zero, it is set to 1/sqrt(max(m,n))
template<typename F>
inline void ADMM
( const DistMatrix<F>& M, DistMatrix<F>& L, DistMatrix<F>& S, 
  const RpcaCtrl<Base<F>>& ctrl )
{
    typedef Base<F> Real;
    const Int m = M.Height();
    const Int n = M.Width();
    const Int commRank = mpi::Rank( M.Grid().Comm() );

    // If tau is not specified, then set it to 1/sqrt(max(m,n))
    const Base<F> tau = 
        ( ctrl.tau <= Real(0) ? 
          Real(1)/sqrt(Real(std::max(m,n))) : 
          ctrl.tau );
    if( ctrl.beta <= Real(0) )
        LogicError("beta cannot be non-positive");
    if( ctrl.tol <= Real(0) )
        LogicError("tol cannot be non-positive");
    const Base<F> beta = ctrl.beta;
    const Base<F> tol = ctrl.tol;

    const double startTime = mpi::Time();
    DistMatrix<F> E( M.Grid() ), Y( M.Grid() );
    Zeros( Y, m, n );

    const Real frobM = FrobeniusNorm( M );
    const Real maxM = MaxNorm( M );
    if( ctrl.progress && commRank == 0 )
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
        if( ctrl.usePivQR )
            rank = SVT( L, Real(1)/beta, ctrl.numPivSteps );
        else
            rank = SVT( L, Real(1)/beta );
      
        // E := M - (L + S)
        E = M;    
        Axpy( F(-1), L, E );
        Axpy( F(-1), S, E );
        const Real frobE = FrobeniusNorm( E );

        if( frobE/frobM <= tol )            
        {
            if( ctrl.progress && commRank == 0 )
                std::cout << "Converged after " << numIts << " iterations "
                          << " with rank=" << rank 
                          << ", numNonzeros=" << numNonzeros << " and "
                          << "|| E ||_F / || M ||_F = " << frobE/frobM
                          << ", and " << mpi::Time()-startTime << " total secs"
                          << std::endl;
            break;
        }
        else if( numIts >= ctrl.maxIts )
        {
            if( ctrl.progress && commRank == 0 )
                std::cout << "Aborting after " << numIts << " iterations and "
                          << mpi::Time()-startTime << " total secs" 
                          << std::endl;
            break;
        }
        else
        {
            if( ctrl.progress && commRank == 0 )
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
inline void ALM
( const DistMatrix<F>& M, DistMatrix<F>& L, DistMatrix<F>& S, 
  const RpcaCtrl<Base<F>>& ctrl )
{
    typedef Base<F> Real;

    const Int m = M.Height();
    const Int n = M.Width();
    const Int commRank = mpi::Rank( M.Grid().Comm() );

    // If tau is unspecified, set it to 1/sqrt(max(m,n))
    const Base<F> tau = 
      ( ctrl.tau <= Real(0) ?
        Real(1) / sqrt(Real(std::max(m,n))) :
        ctrl.tau );
    if( ctrl.tol <= Real(0) )
        LogicError("tol cannot be non-positive");
    const Base<F> tol = ctrl.tol;

    const double startTime = mpi::Time();

    DistMatrix<F> Y( M );
    NormalizeEntries( Y );
    const Real twoNorm = TwoNorm( Y );
    const Real maxNorm = MaxNorm( Y );
    const Real infNorm = maxNorm / tau; 
    const Real dualNorm = std::max( twoNorm, infNorm );
    Scale( F(1)/dualNorm, Y );

    // If beta is unspecified, set it to 1 / 2 || sign(M) ||_2
    Base<F> beta = 
      ( ctrl.beta <= Real(0) ? Real(1) / (2*twoNorm) : ctrl.beta );

    const Real frobM = FrobeniusNorm( M );
    const Real maxM = MaxNorm( M );
    if( ctrl.progress && commRank == 0 )
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
            // HERE
            if( ctrl.usePivQR )
                rank = SVT( L, Real(1)/beta, ctrl.numPivSteps );
            else
                rank = SVT( L, Real(1)/beta );

            Axpy( F(-1), L, LLast );
            Axpy( F(-1), S, SLast );
            const Real frobLDiff = FrobeniusNorm( LLast );
            const Real frobSDiff = FrobeniusNorm( SLast );

            if( frobLDiff/frobM < tol && frobSDiff/frobM < tol )
            {
                if( ctrl.progress && commRank == 0 )
                    std::cout << "Primal loop converged: " 
                              << mpi::Time()-startTime << " total secs"
                              << std::endl;
                break;
            }
            else 
            {
                if( ctrl.progress && commRank == 0 )
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
        Axpy( F(-1), L, E );
        Axpy( F(-1), S, E );
        const Real frobE = FrobeniusNorm( E );

        if( frobE/frobM <= tol )            
        {
            if( ctrl.progress && commRank == 0 )
                std::cout << "Converged after " << numIts << " iterations and "
                          << numPrimalIts << " primal iterations with rank=" 
                          << rank << ", numNonzeros=" << numNonzeros << " and "
                          << "|| E ||_F / || M ||_F = " << frobE/frobM
                          << ", " << mpi::Time()-startTime << " total secs"
                          << std::endl;
            break;
        }
        else if( numIts >= ctrl.maxIts )
        {
            if( ctrl.progress && commRank == 0 )
                std::cout << "Aborting after " << numIts << " iterations and "
                          << mpi::Time()-startTime << " total secs" 
                          << std::endl;
            break;
        }
        else
        {
            if( ctrl.progress && commRank == 0 )
                std::cout << numPrimalIts << ": || E ||_F / || M ||_F = " 
                          << frobE/frobM << ", rank=" << rank 
                          << ", numNonzeros=" << numNonzeros << ", "
                          << mpi::Time()-startTime << " total secs" 
                          << std::endl;
        }
        
        // Y := Y + beta E
        Axpy( beta, E, Y );
        beta *= ctrl.rho;
    }
}

} // namespace rpca

template<typename F>
void RPCA
( const DistMatrix<F>& M, DistMatrix<F>& L, DistMatrix<F>& S,
  const RpcaCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("RPCA"))
    if( ctrl.useALM )
        rpca::ALM( M, L, S, ctrl ); 
    else
        rpca::ADMM( M, L, S, ctrl );
}

#define PROTO(F) \
  template void RPCA \
  ( const DistMatrix<F>& M, DistMatrix<F>& L, DistMatrix<F>& S, \
    const RpcaCtrl<Base<F>>& ctrl );

PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
