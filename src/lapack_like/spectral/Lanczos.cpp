/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// Form 
//
//   A V = V T + v (beta e_{k-1})^H, 
//
// where A is explicitly Hermitian.
//
// T is returned instead of just the tridiagonal portion since later versions
// may account for T numerically having an upper-Hessenberg perturbation.
//
// The current version was implemented for the purposes of gaining rough
// estimates of the extremal singular values in order to keep the condition
// number of the augmented system
//
//      | alpha*I  A | | r/alpha | = | b |
//      |   A^H    0 | |    x    |   | 0 |
//
// as near to that of A as possible (e.g., by setting alpha to roughly the 
// minimum nonzero singular value of A), as well as to provide a means of
// scaling A down to roughly unit two-norm.
//

template<typename F>
void Lanczos
( const SparseMatrix<F>& A, Matrix<Base<F>>& T, Int basisSize )
{
    DEBUG_ONLY(CSE cse("Lanczos"))
    if( A.Height() != A.Width() )
        LogicError("A was not square");
    typedef Base<F> Real;
    const Real eps = Epsilon<Real>();
    const Int n = A.Height();

    Matrix<F> v_km1, v_k, v;
    basisSize = Min(n,basisSize);
    Zeros( v_km1, n, 1 );
    Zeros( v_k,   n, 1 );
    Zeros( T, basisSize, basisSize );

    // Create the initial unit-vector
    // ------------------------------
    {
        Uniform( v, n, 1 );
        Shift( v, SampleUniform<F>() );
        const Real beta = FrobeniusNorm( v ); 
        v *= 1/beta;
    }
        
    // TODO: Incorporate Frobenius norm of A
    const Real minBeta = eps;
    for( Int k=0; k<basisSize; ++k )
    {
        if( k > 0 )
            v_km1 = v_k;
        v_k = v;

        // w := A v
        // --------
        Multiply( NORMAL, F(1), A, v_k, F(0), v );

        // w := w - T(k-1,k) v_{k-1}
        // -------------------------
        if( k > 0 )
            Axpy( -T.Get(k-1,k), v_km1, v );

        // w := w - T(k,k) v_k
        // -------------------
        const Real tau = RealPart(Dot(v_k,v));
        T.Set( k, k, tau );
        Axpy( -tau, v_k, v );

        // v := w / || w ||_2
        // -----------------
        const Real beta = FrobeniusNorm( v );
        if( beta <= minBeta )
        {
            T.Resize( k+1, k+1 );
            break;
        }
        v *= 1/beta;

        // Expand the Rayleigh quotient as appropriate
        // -------------------------------------------
        if( k < basisSize-1 )
        {
            T.Set( k+1, k,   beta );
            T.Set( k,   k+1, beta );
        }
    }
}

template<typename F>
Base<F> LanczosDecomp
( const SparseMatrix<F>& A, Matrix<F>& V, 
        Matrix<Base<F>>& T, Matrix<F>& v,
  Int basisSize )
{
    DEBUG_ONLY(CSE cse("LanczosDecomp"))
    if( A.Height() != A.Width() )
        LogicError("A was not square");
    typedef Base<F> Real;
    const Real eps = Epsilon<Real>();
    const Int n = A.Height();

    basisSize = Min(n,basisSize);
    Zeros( V, n, basisSize );
    Zeros( T, basisSize, basisSize );

    // Initialize the first column of V 
    // --------------------------------
    {
        Uniform( v, n, 1 );
        Shift( v, SampleUniform<F>() );
        const Real beta = FrobeniusNorm( v ); 
        v *= 1/beta;
    }
    if( basisSize > 0 )
    {
        auto v0 = V( ALL, IR(0) );
        v0 = v;
    }
        
    Real beta = 0;
    // TODO: Incorporate Frobenius norm of A
    const Real minBeta = eps;
    for( Int k=0; k<basisSize; ++k )
    {
        // w := A v
        // --------
        auto v_k  = V( ALL, IR(k) );
        Multiply( NORMAL, F(1), A, v_k, F(0), v );

        // w := w - T(k-1,k) v_{k-1}
        // -------------------------
        if( k > 0 )
        {
            auto v_km1 = V( ALL, IR(k-1) );
            Axpy( -T.Get(k-1,k), v_km1, v );
        }

        // w := w - T(k,k) v_k
        // -------------------
        const Real tau = RealPart(Dot(v_k,v));
        T.Set( k, k, tau );
        Axpy( -tau, v_k, v );

        // v := w / || w ||_2
        // -----------------
        beta = FrobeniusNorm( v );
        if( beta <= minBeta )
        {
            T.Resize( k+1, k+1 );
            break;
        }
        v *= 1/beta;

        // Expand the Lanczos decomposition as appropriate
        // -----------------------------------------------
        if( k < basisSize-1 )
        {
            T.Set( k+1, k,   beta );
            T.Set( k,   k+1, beta );
            auto v_kp1 = V( ALL, IR(k+1) );
            v_kp1 = v;
        }
    }
    return beta;
}

template<typename F>
void Lanczos
( const DistSparseMatrix<F>& A, Matrix<Base<F>>& T, Int basisSize )
{
    DEBUG_ONLY(CSE cse("Lanczos"))
    if( A.Height() != A.Width() )
        LogicError("A was not square");
    typedef Base<F> Real;
    const Real eps = Epsilon<Real>();
    const Int n = A.Height();
    mpi::Comm comm = A.Comm();
    const int commRank = mpi::Rank( comm );

    DistMultiVec<F> v_km1(comm), v_k(comm), v(comm);
    basisSize = Min(n,basisSize);
    Zeros( v_km1, n, 1 );
    Zeros( v_k,   n, 1 );
    Zeros( T, basisSize, basisSize );

    // Create the initial unit-vector
    // ------------------------------
    {
        F shift;
        if( commRank == 0 )
            shift = SampleUniform<F>();
        mpi::Broadcast( shift, 0, comm ); 
        Uniform( v, n, 1 );
        Shift( v, shift );
        const Real beta = FrobeniusNorm( v );
        v *= 1/beta;
    }
        
    // TODO: Use the Frobenius norm of A
    const Real minBeta = eps;
    for( Int k=0; k<basisSize; ++k )
    {
        if( k > 0 )
            v_km1 = v_k;
        v_k = v;

        // w := A v
        // --------
        Multiply( NORMAL, F(1), A, v_k, F(0), v );

        // w := w - T(k-1,k) v_{k-1}
        // -------------------------
        if( k > 0 )
            Axpy( -T.Get(k-1,k), v_km1, v );

        // w := w - T(k,k) v_k
        // -------------------
        const Real tau = RealPart(Dot(v_k,v));
        T.Set( k, k, tau );
        Axpy( -tau, v_k, v );

        // v := w / || w ||_2
        // -----------------
        const Real beta = FrobeniusNorm( v );
        if( beta <= minBeta )
        {
            T.Resize( k+1, k+1 );
            break;
        }
        v *= 1/beta;

        // Expand the Rayleigh quotient as appropriate
        // -------------------------------------------
        if( k < basisSize-1 )
        {
            T.Set( k+1, k,   beta );
            T.Set( k,   k+1, beta );
        }
    }
}

template<typename F>
Base<F> LanczosDecomp
( const DistSparseMatrix<F>& A, DistMultiVec<F>& V, 
        Matrix<Base<F>>& T,     DistMultiVec<F>& v,
  Int basisSize )
{
    DEBUG_ONLY(CSE cse("LanczosDecomp"))
    if( A.Height() != A.Width() )
        LogicError("A was not square");
    typedef Base<F> Real;
    const Real eps = Epsilon<Real>();
    const Int n = A.Height();
    mpi::Comm comm = A.Comm();
    const int commRank = mpi::Rank( comm );

    DistMultiVec<F> v_km1(comm), v_k(comm);
    basisSize = Min(n,basisSize);
    Zeros( V, n, basisSize );
    Zeros( T, basisSize, basisSize );
    auto& VLoc = V.Matrix();

    // Choose the initial (unit-length) vector
    // ---------------------------------------
    {
        F shift;
        if( commRank == 0 )
            shift = SampleUniform<F>();
        mpi::Broadcast( shift, 0, comm ); 
        Uniform( v, n, 1 );
        Shift( v, shift );
        const Real beta = FrobeniusNorm( v );
        v *= 1/beta;
    }
    if( basisSize > 0 )
    {
        auto v0Loc = VLoc( ALL, IR(0) );
        v0Loc = v.Matrix();
    }

    Real beta = 0;
    // TODO: Use the Frobenius norm of A
    const Real minBeta = eps;
    for( Int k=0; k<basisSize; ++k )
    {
        // w := A v_k
        // ----------
        v_k  = V( ALL, IR(k) );
        Multiply( NORMAL, F(1), A, v_k, F(0), v );

        // w := w - T(k-1,k) v_{k-1}
        // -------------------------
        if( k > 0 )
        {
            v_km1 = V( ALL, IR(k-1) );
            Axpy( -T.Get(k-1,k), v_km1, v );
        }

        // w := w - T(k,k) v_k
        // -------------------
        const Real tau = RealPart(Dot(v_k,v));
        T.Set( k, k, tau );
        Axpy( -tau, v_k, v );

        // v := w / || w ||_2
        // -----------------
        beta = FrobeniusNorm( v );
        if( beta <= minBeta )
        {
            T.Resize( k+1, k+1 );
            break;
        }
        v *= 1/beta;

        // Expand the Lanczos decomposition as appropriate
        // -----------------------------------------------
        if( k < basisSize-1 )
        {
            T.Set( k+1, k,   beta );
            T.Set( k,   k+1, beta );
            auto v_kp1Loc = VLoc( ALL, IR(k+1) );
            v_kp1Loc = v.Matrix();
        }
    }
    return beta;
}

#define PROTO(F) \
  template void Lanczos \
  ( const SparseMatrix<F>& A, Matrix<Base<F>>& T, Int basisSize ); \
  template void Lanczos \
  ( const DistSparseMatrix<F>& A, Matrix<Base<F>>& T, Int basisSize ); \
  template Base<F> LanczosDecomp \
  ( const SparseMatrix<F>& A, Matrix<F>& V, \
          Matrix<Base<F>>& T, Matrix<F>& v, \
    Int basisSize ); \
  template Base<F> LanczosDecomp \
  ( const DistSparseMatrix<F>& A, DistMultiVec<F>& V, \
          Matrix<Base<F>>& T,     DistMultiVec<F>& v, \
    Int basisSize );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
