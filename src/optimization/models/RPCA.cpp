/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

//
// Please see <http://perception.csl.illinois.edu/matrix-rank/sample_code.html>
// for references and documentation on the Augmented Lagrange Multiplier (ALM) 
// and Alternating Direction Method of Multipliers (ADMM) for Robust PCA.
//

namespace El {

namespace rpca {

template<typename F>
inline void NormalizeEntries( Matrix<F>& A )
{ 
    auto unitMap = []( F alpha ) 
                   { return alpha==F(0) ? F(1) : alpha/Abs(alpha); };
    EntrywiseMap( A, function<F(F)>(unitMap) );
}

template<typename F>
inline void NormalizeEntries( ElementalMatrix<F>& A )
{ 
    auto unitMap = []( F alpha ) 
                   { return alpha==F(0) ? F(1) : alpha/Abs(alpha); };
    EntrywiseMap( A, function<F(F)>(unitMap) );
}

// NOTE: If 'tau' is passed in as zero, it is set to 1/sqrt(max(m,n))

template<typename F>
inline void ADMM
( const Matrix<F>& M, 
        Matrix<F>& L,
        Matrix<F>& S, 
  const RPCACtrl<Base<F>>& ctrl )
{
    typedef Base<F> Real;
    const Int m = M.Height();
    const Int n = M.Width();

    // If tau is not specified, then set it to 1/sqrt(max(m,n))
    const Base<F> tau = 
        ( ctrl.tau <= Real(0) ? Real(1)/sqrt(Real(Max(m,n))) : ctrl.tau );
    if( ctrl.beta <= Real(0) )
        LogicError("beta cannot be non-positive");
    if( ctrl.tol <= Real(0) )
        LogicError("tol cannot be non-positive");
    const Base<F> beta = ctrl.beta;
    const Base<F> tol = ctrl.tol;

    const double startTime = mpi::Time();
    Matrix<F> E, Y;
    Zeros( Y, m, n );

    const Real frobM = FrobeniusNorm( M );
    const Real maxM = MaxNorm( M );
    if( ctrl.progress )
        cout << "|| M ||_F = " << frobM << "\n"
             << "|| M ||_max = " << maxM << endl;

    Zeros( L, m, n );
    Zeros( S, m, n );

    Int numIts = 0;
    while( true )
    {
        ++numIts;

        // ST_{tau/beta}(M - L + Y/beta)
        S = M;
        S -= L;
        Axpy( F(1)/beta, Y, S );
        SoftThreshold( S, tau/beta );
        const Int numNonzeros = ZeroNorm( S );

        // SVT_{1/beta}(M - S + Y/beta)
        L = M;
        L -= S;
        Axpy( F(1)/beta, Y, L );
        Int rank;
        if( ctrl.usePivQR )
            rank = SVT( L, Real(1)/beta, ctrl.numPivSteps );
        else
            rank = SVT( L, Real(1)/beta );
      
        // E := M - (L + S)
        E = M;    
        E -= L;
        E -= S;
        const Real frobE = FrobeniusNorm( E );

        if( frobE/frobM <= tol )            
        {
            if( ctrl.progress )
                cout << "Converged after " << numIts << " iterations "
                     << " with rank=" << rank 
                     << ", numNonzeros=" << numNonzeros << " and "
                     << "|| E ||_F / || M ||_F = " << frobE/frobM
                     << ", and " << mpi::Time()-startTime << " total secs"
                     << endl;
            break;
        }
        else if( numIts >= ctrl.maxIts )
        {
            if( ctrl.progress )
                cout << "Aborting after " << numIts << " iterations and "
                     << mpi::Time()-startTime << " total secs" 
                     << endl;
            break;
        }
        else
        {
            if( ctrl.progress )
                cout << numIts << ": || E ||_F / || M ||_F = " 
                     << frobE/frobM << ", rank=" << rank 
                     << ", numNonzeros=" << numNonzeros 
                     << ", " << mpi::Time()-startTime << " total secs"
                     << endl;
        }
        
        // Y := Y + beta E
        Axpy( beta, E, Y );
    }
}

template<typename F>
inline void ADMM
( const ElementalMatrix<F>& MPre,
        ElementalMatrix<F>& LPre, 
        ElementalMatrix<F>& SPre,
  const RPCACtrl<Base<F>>& ctrl )
{
    DistMatrixReadProxy<F,F,MC,MR>
      MProx( MPre );
    DistMatrixWriteProxy<F,F,MC,MR>
      LProx( LPre ),
      SProx( SPre );
    auto& M = MProx.GetLocked();
    auto& L = LProx.Get();
    auto& S = SProx.Get();

    typedef Base<F> Real;
    const Int m = M.Height();
    const Int n = M.Width();
    const int commRank = mpi::Rank( M.Grid().Comm() );

    // If tau is not specified, then set it to 1/sqrt(max(m,n))
    const Base<F> tau = 
        ( ctrl.tau <= Real(0) ? Real(1)/sqrt(Real(Max(m,n))) : ctrl.tau );
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
        cout << "|| M ||_F = " << frobM << "\n"
             << "|| M ||_max = " << maxM << endl;

    Zeros( L, m, n );
    Zeros( S, m, n );

    Int numIts = 0;
    while( true )
    {
        ++numIts;

        // ST_{tau/beta}(M - L + Y/beta)
        S = M;
        S -= L;
        Axpy( F(1)/beta, Y, S );
        SoftThreshold( S, tau/beta );
        const Int numNonzeros = ZeroNorm( S );

        // SVT_{1/beta}(M - S + Y/beta)
        L = M;
        L -= S;
        Axpy( F(1)/beta, Y, L );
        Int rank;
        if( ctrl.usePivQR )
            rank = SVT( L, Real(1)/beta, ctrl.numPivSteps );
        else
            rank = SVT( L, Real(1)/beta );
      
        // E := M - (L + S)
        E = M;    
        E -= L;
        E -= S;
        const Real frobE = FrobeniusNorm( E );

        if( frobE/frobM <= tol )            
        {
            if( ctrl.progress && commRank == 0 )
                cout << "Converged after " << numIts << " iterations "
                     << " with rank=" << rank 
                     << ", numNonzeros=" << numNonzeros << " and "
                     << "|| E ||_F / || M ||_F = " << frobE/frobM
                     << ", and " << mpi::Time()-startTime << " total secs"
                     << endl;
            break;
        }
        else if( numIts >= ctrl.maxIts )
        {
            if( ctrl.progress && commRank == 0 )
                cout << "Aborting after " << numIts << " iterations and "
                     << mpi::Time()-startTime << " total secs" 
                     << endl;
            break;
        }
        else
        {
            if( ctrl.progress && commRank == 0 )
                cout << numIts << ": || E ||_F / || M ||_F = " 
                     << frobE/frobM << ", rank=" << rank 
                     << ", numNonzeros=" << numNonzeros 
                     << ", " << mpi::Time()-startTime << " total secs"
                     << endl;
        }
        
        // Y := Y + beta E
        Axpy( beta, E, Y );
    }
}

// NOTE: If 'beta' or 'tau' is zero, then an estimate is used instead

template<typename F>
inline void ALM
( const Matrix<F>& M,
        Matrix<F>& L,
        Matrix<F>& S, 
  const RPCACtrl<Base<F>>& ctrl )
{
    typedef Base<F> Real;
    const Int m = M.Height();
    const Int n = M.Width();

    // If tau is unspecified, set it to 1/sqrt(max(m,n))
    const Base<F> tau = 
      ( ctrl.tau <= Real(0) ? Real(1) / sqrt(Real(Max(m,n))) :
        ctrl.tau );
    if( ctrl.tol <= Real(0) )
        LogicError("tol cannot be non-positive");
    const Base<F> tol = ctrl.tol;

    const double startTime = mpi::Time();

    Matrix<F> Y( M );
    NormalizeEntries( Y );
    const Real twoNorm = TwoNorm( Y );
    const Real maxNorm = MaxNorm( Y );
    const Real infNorm = maxNorm / tau; 
    const Real dualNorm = Max( twoNorm, infNorm );
    Y *= F(1)/dualNorm;

    // If beta is unspecified, set it to 1 / 2 || sign(M) ||_2
    Base<F> beta = 
      ( ctrl.beta <= Real(0) ? Real(1) / (2*twoNorm) : ctrl.beta );

    const Real frobM = FrobeniusNorm( M );
    const Real maxM = MaxNorm( M );
    if( ctrl.progress )
        cout << "|| M ||_F = " << frobM << "\n"
             << "|| M ||_max = " << maxM << endl;

    Zeros( L, m, n );
    Zeros( S, m, n );

    Int numIts=0, numPrimalIts=0;
    Matrix<F> LLast, SLast, E;
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
            S -= L;
            Axpy( F(1)/beta, Y, S );
            SoftThreshold( S, tau/beta );
            numNonzeros = ZeroNorm( S );

            // SVT_{1/beta}(M - S + Y/beta)
            L = M;
            L -= S;
            Axpy( F(1)/beta, Y, L );
            if( ctrl.usePivQR )
                rank = SVT( L, Real(1)/beta, ctrl.numPivSteps );
            else
                rank = SVT( L, Real(1)/beta );

            LLast -= L;
            SLast -= S;
            const Real frobLDiff = FrobeniusNorm( LLast );
            const Real frobSDiff = FrobeniusNorm( SLast );

            if( frobLDiff/frobM < tol && frobSDiff/frobM < tol )
            {
                if( ctrl.progress )
                    cout << "Primal loop converged: " 
                         << mpi::Time()-startTime << " total secs"
                         << endl;
                break;
            }
            else 
            {
                if( ctrl.progress )
                    cout << "  " << numPrimalIts 
                         << ": \n"
                         << "   || Delta L ||_F / || M ||_F = " 
                         << frobLDiff/frobM << "\n"
                         << "   || Delta S ||_F / || M ||_F = "
                         << frobSDiff/frobM << "\n"
                         << "   rank=" << rank
                         << ", numNonzeros=" << numNonzeros 
                         << ", " << mpi::Time()-startTime << " total secs" 
                         << endl;
            } 
        }

        // E := M - (L + S)
        E = M;    
        E -= L;
        E -= S;
        const Real frobE = FrobeniusNorm( E );

        if( frobE/frobM <= tol )            
        {
            if( ctrl.progress )
                cout << "Converged after " << numIts << " iterations and "
                     << numPrimalIts << " primal iterations with rank=" 
                     << rank << ", numNonzeros=" << numNonzeros << " and "
                     << "|| E ||_F / || M ||_F = " << frobE/frobM
                     << ", " << mpi::Time()-startTime << " total secs"
                     << endl;
            break;
        }
        else if( numIts >= ctrl.maxIts )
        {
            if( ctrl.progress )
                cout << "Aborting after " << numIts << " iterations and "
                     << mpi::Time()-startTime << " total secs" 
                     << endl;
            break;
        }
        else
        {
            if( ctrl.progress )
                cout << numPrimalIts << ": || E ||_F / || M ||_F = " 
                     << frobE/frobM << ", rank=" << rank 
                     << ", numNonzeros=" << numNonzeros << ", "
                     << mpi::Time()-startTime << " total secs" 
                     << endl;
        }
        
        // Y := Y + beta E
        Axpy( beta, E, Y );
        beta *= ctrl.rho;
    }
}

template<typename F>
inline void ALM
( const ElementalMatrix<F>& MPre, 
        ElementalMatrix<F>& LPre,
        ElementalMatrix<F>& SPre, 
  const RPCACtrl<Base<F>>& ctrl )
{
    DistMatrixReadProxy<F,F,MC,MR>
      MProx( MPre );
    DistMatrixWriteProxy<F,F,MC,MR>
      LProx( LPre ),
      SProx( SPre );
    auto& M = MProx.GetLocked();
    auto& L = LProx.Get();
    auto& S = SProx.Get();

    typedef Base<F> Real;
    const Int m = M.Height();
    const Int n = M.Width();
    const int commRank = mpi::Rank( M.Grid().Comm() );

    // If tau is unspecified, set it to 1/sqrt(max(m,n))
    const Base<F> tau = 
      ( ctrl.tau <= Real(0) ? Real(1) / sqrt(Real(Max(m,n))) : ctrl.tau );
    if( ctrl.tol <= Real(0) )
        LogicError("tol cannot be non-positive");
    const Base<F> tol = ctrl.tol;

    const double startTime = mpi::Time();

    DistMatrix<F> Y( M );
    NormalizeEntries( Y );
    const Real twoNorm = TwoNorm( Y );
    const Real maxNorm = MaxNorm( Y );
    const Real infNorm = maxNorm / tau; 
    const Real dualNorm = Max( twoNorm, infNorm );
    Y *= F(1)/dualNorm;

    // If beta is unspecified, set it to 1 / 2 || sign(M) ||_2
    Base<F> beta = 
      ( ctrl.beta <= Real(0) ? Real(1) / (2*twoNorm) : ctrl.beta );

    const Real frobM = FrobeniusNorm( M );
    const Real maxM = MaxNorm( M );
    if( ctrl.progress && commRank == 0 )
        cout << "|| M ||_F = " << frobM << "\n"
             << "|| M ||_max = " << maxM << endl;

    Zeros( L, m, n );
    Zeros( S, m, n );

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
            S -= L;
            Axpy( F(1)/beta, Y, S );
            SoftThreshold( S, tau/beta );
            numNonzeros = ZeroNorm( S );

            // SVT_{1/beta}(M - S + Y/beta)
            L = M;
            L -= S;
            Axpy( F(1)/beta, Y, L );
            if( ctrl.usePivQR )
                rank = SVT( L, Real(1)/beta, ctrl.numPivSteps );
            else
                rank = SVT( L, Real(1)/beta );

            LLast -= L;
            SLast -= S;
            const Real frobLDiff = FrobeniusNorm( LLast );
            const Real frobSDiff = FrobeniusNorm( SLast );

            if( frobLDiff/frobM < tol && frobSDiff/frobM < tol )
            {
                if( ctrl.progress && commRank == 0 )
                    cout << "Primal loop converged: " 
                         << mpi::Time()-startTime << " total secs" << endl;
                break;
            }
            else 
            {
                if( ctrl.progress && commRank == 0 )
                    cout << "  " << numPrimalIts 
                         << ": \n"
                         << "   || Delta L ||_F / || M ||_F = " 
                         << frobLDiff/frobM << "\n"
                         << "   || Delta S ||_F / || M ||_F = "
                         << frobSDiff/frobM << "\n"
                         << "   rank=" << rank
                         << ", numNonzeros=" << numNonzeros 
                         << ", " << mpi::Time()-startTime << " total secs" 
                         << endl;
            } 
        }

        // E := M - (L + S)
        E = M;    
        E -= L;
        E -= S;
        const Real frobE = FrobeniusNorm( E );

        if( frobE/frobM <= tol )            
        {
            if( ctrl.progress && commRank == 0 )
                cout << "Converged after " << numIts << " iterations and "
                     << numPrimalIts << " primal iterations with rank=" 
                     << rank << ", numNonzeros=" << numNonzeros << " and "
                     << "|| E ||_F / || M ||_F = " << frobE/frobM
                     << ", " << mpi::Time()-startTime << " total secs"
                     << endl;
            break;
        }
        else if( numIts >= ctrl.maxIts )
        {
            if( ctrl.progress && commRank == 0 )
                cout << "Aborting after " << numIts << " iterations and "
                     << mpi::Time()-startTime << " total secs" << endl;
            break;
        }
        else
        {
            if( ctrl.progress && commRank == 0 )
                cout << numPrimalIts << ": || E ||_F / || M ||_F = " 
                     << frobE/frobM << ", rank=" << rank 
                     << ", numNonzeros=" << numNonzeros << ", "
                     << mpi::Time()-startTime << " total secs" 
                     << endl;
        }
        
        // Y := Y + beta E
        Axpy( beta, E, Y );
        beta *= ctrl.rho;
    }
}

} // namespace rpca

template<typename F>
void RPCA
( const Matrix<F>& M,
        Matrix<F>& L,
        Matrix<F>& S,
  const RPCACtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("RPCA"))
    if( ctrl.useALM )
        rpca::ALM( M, L, S, ctrl ); 
    else
        rpca::ADMM( M, L, S, ctrl );
}

template<typename F>
void RPCA
( const ElementalMatrix<F>& M,
        ElementalMatrix<F>& L, 
        ElementalMatrix<F>& S,
  const RPCACtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("RPCA"))
    if( ctrl.useALM )
        rpca::ALM( M, L, S, ctrl ); 
    else
        rpca::ADMM( M, L, S, ctrl );
}

#define PROTO(F) \
  template void RPCA \
  ( const Matrix<F>& M, \
          Matrix<F>& L, \
          Matrix<F>& S, \
    const RPCACtrl<Base<F>>& ctrl ); \
  template void RPCA \
  ( const ElementalMatrix<F>& M, \
          ElementalMatrix<F>& L, \
          ElementalMatrix<F>& S, \
    const RPCACtrl<Base<F>>& ctrl );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
