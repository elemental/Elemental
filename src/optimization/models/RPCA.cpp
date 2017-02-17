/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

//
// Please see <http://perception.csl.illinois.edu/matrix-rank/sample_code.html>
// for references and documentation on the Augmented Lagrange Multiplier (ALM)
// and Alternating Direction Method of Multipliers (ADMM) for Robust PCA.
//

namespace El {

namespace rpca {

template<typename Field>
void NormalizeEntries( Matrix<Field>& A )
{
    auto unitMap = []( const Field& alpha )
      { return alpha==Field(0) ? Field(1) : alpha/Abs(alpha); };
    EntrywiseMap( A, MakeFunction(unitMap) );
}

template<typename Field>
void NormalizeEntries( AbstractDistMatrix<Field>& A )
{
    auto unitMap = []( const Field& alpha )
      { return alpha==Field(0) ? Field(1) : alpha/Abs(alpha); };
    EntrywiseMap( A, MakeFunction(unitMap) );
}

// NOTE: If 'tau' is passed in as zero, it is set to 1/sqrt(max(m,n))

template<typename Field>
void ADMM
( const Matrix<Field>& M,
        Matrix<Field>& L,
        Matrix<Field>& S,
  const RPCACtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Int m = M.Height();
    const Int n = M.Width();

    // If tau is not specified, then set it to 1/sqrt(max(m,n))
    const Real tau =
        ( ctrl.tau <= Real(0) ? Real(1)/sqrt(Real(Max(m,n))) : ctrl.tau );
    if( ctrl.beta <= Real(0) )
        LogicError("beta cannot be non-positive");
    if( ctrl.tol <= Real(0) )
        LogicError("tol cannot be non-positive");
    const Real beta = ctrl.beta;
    const Real tol = ctrl.tol;

    const double startTime = mpi::Time();
    Matrix<Field> E, Y;
    Zeros( Y, m, n );

    const Real frobM = FrobeniusNorm( M );
    const Real maxM = MaxNorm( M );
    if( ctrl.progress )
    {
        Output("|| M ||_F = ",frobM);
        Output("|| M ||_max = ",maxM);
    }

    Zeros( L, m, n );
    Zeros( S, m, n );

    Int numIts = 0;
    while( true )
    {
        ++numIts;

        // ST_{tau/beta}(M - L + Y/beta)
        S = M;
        S -= L;
        Axpy( Field(1)/beta, Y, S );
        SoftThreshold( S, tau/beta );
        const Int numNonzeros = ZeroNorm( S );

        // SVT_{1/beta}(M - S + Y/beta)
        L = M;
        L -= S;
        Axpy( Field(1)/beta, Y, L );
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
                Output
                ("Converged after ",numIts," iterations "," with rank=",rank,
                 ", numNonzeros=",numNonzeros," and || E ||_F / || M ||_F = ",
                 frobE/frobM,", and ",mpi::Time()-startTime," total secs");
            break;
        }
        else if( numIts >= ctrl.maxIts )
        {
            if( ctrl.progress )
                Output
                ("Aborting after ",numIts," iterations and ",
                 mpi::Time()-startTime," total secs");
            break;
        }
        else
        {
            if( ctrl.progress )
                Output
                (numIts,": || E ||_F / || M ||_F = ",frobE/frobM,", rank=",rank,
                 ", numNonzeros=",numNonzeros,", ",mpi::Time()-startTime,
                 " total secs");
        }

        // Y := Y + beta E
        Axpy( beta, E, Y );
    }
}

template<typename Field>
void ADMM
( const AbstractDistMatrix<Field>& MPre,
        AbstractDistMatrix<Field>& LPre,
        AbstractDistMatrix<Field>& SPre,
  const RPCACtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    DistMatrixReadProxy<Field,Field,MC,MR>
      MProx( MPre );
    DistMatrixWriteProxy<Field,Field,MC,MR>
      LProx( LPre ),
      SProx( SPre );
    auto& M = MProx.GetLocked();
    auto& L = LProx.Get();
    auto& S = SProx.Get();

    typedef Base<Field> Real;
    const Int m = M.Height();
    const Int n = M.Width();
    const int commRank = mpi::Rank( M.Grid().Comm() );

    // If tau is not specified, then set it to 1/sqrt(max(m,n))
    const Real tau =
        ( ctrl.tau <= Real(0) ? Real(1)/sqrt(Real(Max(m,n))) : ctrl.tau );
    if( ctrl.beta <= Real(0) )
        LogicError("beta cannot be non-positive");
    if( ctrl.tol <= Real(0) )
        LogicError("tol cannot be non-positive");
    const Real beta = ctrl.beta;
    const Real tol = ctrl.tol;

    const double startTime = mpi::Time();
    DistMatrix<Field> E( M.Grid() ), Y( M.Grid() );
    Zeros( Y, m, n );

    const Real frobM = FrobeniusNorm( M );
    const Real maxM = MaxNorm( M );
    if( ctrl.progress && commRank == 0 )
    {
        Output("|| M ||_F = ",frobM);
        Output("|| M ||_max = ",maxM);
    }

    Zeros( L, m, n );
    Zeros( S, m, n );

    Int numIts = 0;
    while( true )
    {
        ++numIts;

        // ST_{tau/beta}(M - L + Y/beta)
        S = M;
        S -= L;
        Axpy( Field(1)/beta, Y, S );
        SoftThreshold( S, tau/beta );
        const Int numNonzeros = ZeroNorm( S );

        // SVT_{1/beta}(M - S + Y/beta)
        L = M;
        L -= S;
        Axpy( Field(1)/beta, Y, L );
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
                Output
                ("Converged after ",numIts," iterations with rank=",rank,
                 ", numNonzeros=",numNonzeros," and || E ||_F / || M ||_F = ",
                 frobE/frobM,", and ",mpi::Time()-startTime," total secs");
            break;
        }
        else if( numIts >= ctrl.maxIts )
        {
            if( ctrl.progress && commRank == 0 )
                Output
                ("Aborting after ",numIts," iterations and ",
                 mpi::Time()-startTime," total secs");
            break;
        }
        else
        {
            if( ctrl.progress && commRank == 0 )
                Output
                (numIts,": || E ||_F / || M ||_F = ",frobE/frobM,", rank=",
                 rank,", numNonzeros=",numNonzeros,", ",mpi::Time()-startTime,
                 " total secs");
        }

        // Y := Y + beta E
        Axpy( beta, E, Y );
    }
}

// NOTE: If 'beta' or 'tau' is zero, then an estimate is used instead

template<typename Field>
void ALM
( const Matrix<Field>& M,
        Matrix<Field>& L,
        Matrix<Field>& S,
  const RPCACtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Int m = M.Height();
    const Int n = M.Width();

    // If tau is unspecified, set it to 1/sqrt(max(m,n))
    const Real tau =
      ( ctrl.tau <= Real(0) ? Real(1) / sqrt(Real(Max(m,n))) :
        ctrl.tau );
    if( ctrl.tol <= Real(0) )
        LogicError("tol cannot be non-positive");
    const Real tol = ctrl.tol;

    const double startTime = mpi::Time();

    Matrix<Field> Y( M );
    NormalizeEntries( Y );
    const Real twoNorm = TwoNorm( Y );
    const Real maxNorm = MaxNorm( Y );
    const Real infNorm = maxNorm / tau;
    const Real dualNorm = Max( twoNorm, infNorm );
    Y *= Field(1)/dualNorm;

    // If beta is unspecified, set it to 1 / 2 || sign(M) ||_2
    Real beta = ( ctrl.beta <= Real(0) ? Real(1) / (2*twoNorm) : ctrl.beta );

    const Real frobM = FrobeniusNorm( M );
    const Real maxM = MaxNorm( M );
    if( ctrl.progress )
    {
        Output("|| M ||_F = ",frobM);
        Output("|| M ||_max = ",maxM);
    }

    Zeros( L, m, n );
    Zeros( S, m, n );

    Int numIts=0, numPrimalIts=0;
    Matrix<Field> LLast, SLast, E;
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
            Axpy( Field(1)/beta, Y, S );
            SoftThreshold( S, tau/beta );
            numNonzeros = ZeroNorm( S );

            // SVT_{1/beta}(M - S + Y/beta)
            L = M;
            L -= S;
            Axpy( Field(1)/beta, Y, L );
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
                    Output
                    ("Primal loop converged: ",mpi::Time()-startTime,
                     " total secs");
                break;
            }
            else
            {
                if( ctrl.progress )
                {
                    Output("  ",numPrimalIts,": ");
                    Output("   || Delta L ||_F / || M ||_F = ",frobLDiff/frobM);
                    Output("   || Delta S ||_F / || M ||_F = ",frobSDiff/frobM);
                    Output
                    ("   rank=",rank,", numNonzeros=",numNonzeros,", ",
                     mpi::Time()-startTime," total secs");
                }
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
                Output
                ("Converged after ",numIts," iterations and ",numPrimalIts,
                 " primal iterations with rank=",rank,", numNonzeros=",
                 numNonzeros," and || E ||_F / || M ||_F = ",frobE/frobM,", ",
                 mpi::Time()-startTime," total secs");
            break;
        }
        else if( numIts >= ctrl.maxIts )
        {
            if( ctrl.progress )
                Output
                ("Aborting after ",numIts," iterations and ",
                 mpi::Time()-startTime," total secs");
            break;
        }
        else
        {
            if( ctrl.progress )
                Output
                (numPrimalIts,": || E ||_F / || M ||_F = ",frobE/frobM,
                 ", rank=",rank,", numNonzeros=",numNonzeros,", ",
                 mpi::Time()-startTime," total secs");
        }

        // Y := Y + beta E
        Axpy( beta, E, Y );
        beta *= ctrl.rho;
    }
}

template<typename Field>
void ALM
( const AbstractDistMatrix<Field>& MPre,
        AbstractDistMatrix<Field>& LPre,
        AbstractDistMatrix<Field>& SPre,
  const RPCACtrl<Base<Field>>& ctrl )
{
    DistMatrixReadProxy<Field,Field,MC,MR>
      MProx( MPre );
    DistMatrixWriteProxy<Field,Field,MC,MR>
      LProx( LPre ),
      SProx( SPre );
    auto& M = MProx.GetLocked();
    auto& L = LProx.Get();
    auto& S = SProx.Get();

    typedef Base<Field> Real;
    const Int m = M.Height();
    const Int n = M.Width();
    const int commRank = mpi::Rank( M.Grid().Comm() );

    // If tau is unspecified, set it to 1/sqrt(max(m,n))
    const Real tau =
      ( ctrl.tau <= Real(0) ? Real(1) / sqrt(Real(Max(m,n))) : ctrl.tau );
    if( ctrl.tol <= Real(0) )
        LogicError("tol cannot be non-positive");
    const Real tol = ctrl.tol;

    const double startTime = mpi::Time();

    DistMatrix<Field> Y( M );
    NormalizeEntries( Y );
    const Real twoNorm = TwoNorm( Y );
    const Real maxNorm = MaxNorm( Y );
    const Real infNorm = maxNorm / tau;
    const Real dualNorm = Max( twoNorm, infNorm );
    Y *= Field(1)/dualNorm;

    // If beta is unspecified, set it to 1 / 2 || sign(M) ||_2
    Real beta =
      ( ctrl.beta <= Real(0) ? Real(1) / (2*twoNorm) : ctrl.beta );

    const Real frobM = FrobeniusNorm( M );
    const Real maxM = MaxNorm( M );
    if( ctrl.progress && commRank == 0 )
    {
        Output("|| M ||_F = ",frobM);
        Output("|| M ||_max = ",maxM);
    }

    Zeros( L, m, n );
    Zeros( S, m, n );

    Int numIts=0, numPrimalIts=0;
    DistMatrix<Field> LLast( M.Grid() ), SLast( M.Grid() ), E( M.Grid() );
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
            Axpy( Field(1)/beta, Y, S );
            SoftThreshold( S, tau/beta );
            numNonzeros = ZeroNorm( S );

            // SVT_{1/beta}(M - S + Y/beta)
            L = M;
            L -= S;
            Axpy( Field(1)/beta, Y, L );
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
                    Output
                    ("Primal loop converged: ",mpi::Time()-startTime,
                     " total secs");
                break;
            }
            else
            {
                if( ctrl.progress && commRank == 0 )
                {
                    Output("  ",numPrimalIts,":");
                    Output("   || Delta L ||_F / || M ||_F = ",frobLDiff/frobM);
                    Output("   || Delta S ||_F / || M ||_F = ",frobSDiff/frobM);
                    Output
                    ("   rank=",rank,", numNonzeros=",numNonzeros,
                     ", ",mpi::Time()-startTime," total secs");
                }
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
                Output
                ("Converged after ",numIts," iterations and ",numPrimalIts,
                 " primal iterations with rank=",rank,", numNonzeros=",
                 numNonzeros," and || E ||_F / || M ||_F = ",frobE/frobM,
                 ", ",mpi::Time()-startTime," total secs");
            break;
        }
        else if( numIts >= ctrl.maxIts )
        {
            if( ctrl.progress && commRank == 0 )
                Output
                ("Aborting after ",numIts," iterations and ",
                 mpi::Time()-startTime," total secs");
            break;
        }
        else
        {
            if( ctrl.progress && commRank == 0 )
                Output
                (numPrimalIts,": || E ||_F / || M ||_F = ",frobE/frobM,
                 ", rank=",rank,", numNonzeros=",numNonzeros,", ",
                 mpi::Time()-startTime," total secs");
        }

        // Y := Y + beta E
        Axpy( beta, E, Y );
        beta *= ctrl.rho;
    }
}

} // namespace rpca

template<typename Field>
void RPCA
( const Matrix<Field>& M,
        Matrix<Field>& L,
        Matrix<Field>& S,
  const RPCACtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    if( ctrl.useALM )
        rpca::ALM( M, L, S, ctrl );
    else
        rpca::ADMM( M, L, S, ctrl );
}

template<typename Field>
void RPCA
( const AbstractDistMatrix<Field>& M,
        AbstractDistMatrix<Field>& L,
        AbstractDistMatrix<Field>& S,
  const RPCACtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    if( ctrl.useALM )
        rpca::ALM( M, L, S, ctrl );
    else
        rpca::ADMM( M, L, S, ctrl );
}

#define PROTO(Field) \
  template void RPCA \
  ( const Matrix<Field>& M, \
          Matrix<Field>& L, \
          Matrix<Field>& S, \
    const RPCACtrl<Base<Field>>& ctrl ); \
  template void RPCA \
  ( const AbstractDistMatrix<Field>& M, \
          AbstractDistMatrix<Field>& L, \
          AbstractDistMatrix<Field>& S, \
    const RPCACtrl<Base<Field>>& ctrl );

#define EL_NO_INT_PROTO
#include <El/macros/Instantiate.h>

} // namespace El
