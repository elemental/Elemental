/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_SCHUR_SDC_HPP
#define ELEM_LAPACK_SCHUR_SDC_HPP

#include "elemental/blas-like/level1/Axpy.hpp"
#include "elemental/blas-like/level1/Scale.hpp"
#include "elemental/blas-like/level1/UpdateDiagonal.hpp"
#include "elemental/lapack-like/Median.hpp"
#include "elemental/lapack-like/Norm/Infinity.hpp"
#include "elemental/lapack-like/Norm/One.hpp"
#include "elemental/lapack-like/QR.hpp"
#include "elemental/lapack-like/Sign.hpp"
#include "elemental/matrices/Haar.hpp"

// See Z. Bai, J. Demmel, J. Dongarra, A. Petitet, H. Robinson, and K. Stanley's
// "The spectral decomposition of nonsymmetric matrices on distributed memory
// parallel computers". Currently available at:
// www.netlib.org/lapack/lawnspdf/lawn91.pdf
//
// as well as the improved version, which avoids pivoted QR, in J. Demmel, 
// I. Dumitriu, and O. Holtz, "Fast linear algebra is stable", 2007.
// www.netlib.org/lapack/lawnspdf/lawn186.pdf

namespace elem {
namespace schur {

template<typename F>
inline ValueInt<BASE(F)>
ComputePartition( Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry cse("schur::ComputePartition");
#endif
    typedef BASE(F) Real;
    const Int n = A.Height();
    if( n == 0 ) 
    {
        ValueInt<Real> part;
        part.value = -1;
        part.index = -1;
        return;
    }

    // Compute the sets of row and column sums
    std::vector<Real> colSums(n-1,0), rowSums(n-1,0);
    for( Int j=0; j<n-1; ++j )
        for( Int i=j+1; i<n; ++i )
            colSums[j] += Abs( A.Get(i,j) ); 
    for( Int i=1; i<n-1; ++i )
        for( Int j=0; j<i; ++j )
            rowSums[i-1] += Abs( A.Get(i,j) );

    // Compute the list of norms and its minimum value/index
    ValueInt<Real> part;
    std::vector<Real> norms(n-1);
    norms[0] = colSums[0];
    part.value = norms[0];
    part.index = 1;
    for( Int j=1; j<n-1; ++j )
    {
        norms[j] = norms[j-1] + colSums[j] - rowSums[j-1];
        if( norms[j] < part.value )
        {
            part.value = norms[j];
            part.index = j+1;
        }
    }

    return part;
}

// The current implementation requires O(n^2/p + n lg p) work. Since the
// matrix-matrix multiplication alone requires O(n^3/p) work, and n <= p for
// most practical computations, it is at least O(n^2) work, which should dwarf
// the O(n lg p) unparallelized component of this algorithm.
template<typename F>
inline ValueInt<BASE(F)>
ComputePartition( DistMatrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry cse("schur::ComputePartition");
#endif
    typedef BASE(F) Real;
    const Grid& g = A.Grid();
    const Int n = A.Height();
    if( n == 0 ) 
    {
        ValueInt<Real> part;
        part.value = -1;
        part.index = -1;
        return part;
    }

    // Compute the sets of row and column sums
    std::vector<Real> colSums(n-1,0), rowSums(n-1,0);
    const Int mLocal = A.LocalHeight();
    const Int nLocal = A.LocalWidth();
    const Int rowShift = A.RowShift();
    const Int colShift = A.ColShift();
    const Int rowStride = A.RowStride();
    const Int colStride = A.ColStride();
    for( Int jLoc=0; jLoc<nLocal; ++jLoc )
    {
        const Int j = rowShift + jLoc*rowStride;
        if( j < n-1 )
        {
            for( Int iLoc=0; iLoc<mLocal; ++iLoc )
            {
                const Int i = colShift + iLoc*colStride;
                if( i > j )
                {
                    colSums[j] += Abs( A.GetLocal(iLoc,jLoc) ); 
                    rowSums[i-1] += Abs( A.GetLocal(iLoc,jLoc) );
                }
            }
        }
    }
    mpi::AllReduce( &colSums[0], n-1, g.VCComm() );
    mpi::AllReduce( &rowSums[0], n-1, g.VCComm() );

    // Compute the list of norms and its minimum value/index
    // TODO: Think of the proper way to parallelize this if necessary
    ValueInt<Real> part;
    std::vector<Real> norms(n-1);
    norms[0] = colSums[0];
    part.value = norms[0];
    part.index = 1;
    for( Int j=1; j<n-1; ++j )
    {
        norms[j] = norms[j-1] + colSums[j] - rowSums[j-1];
        if( norms[j] < part.value )
        {
            part.value = norms[j];
            part.index = j+1;
        }
    }

    return part;
}

// G should be a rational function of A. If returnQ=true, G will be set to
// the computed unitary matrix upon exit.
template<typename F>
inline ValueInt<BASE(F)>
SignDivide( Matrix<F>& A, Matrix<F>& G, bool returnQ=false )
{
#ifndef RELEASE
    CallStackEntry cse("schur::SignDivide");
#endif
    typedef BASE(F) Real;
    const Int n = A.Height();

    // G := sgn(G)
    // G := 1/2 ( G + I )
    Sign( G );
    UpdateDiagonal( G, F(1) );
    Scale( F(1)/F(2), G );

    // Compute the pivoted QR decomposition of the spectral projection 
    Matrix<F> t;
    Matrix<Int> p;
    elem::QR( G, t, p );

    // A := Q^H A Q
    const Real oneA = OneNorm( A );
    if( returnQ )
    {
        ExpandPackedReflectors( LOWER, VERTICAL, UNCONJUGATED, 0, G, t );
        Matrix<F> B;
        Gemm( ADJOINT, NORMAL, F(1), G, A, B );
        Gemm( NORMAL, NORMAL, F(1), B, G, A );
    }
    else
    {
        qr::ApplyQ( LEFT, ADJOINT, G, t, A );
        qr::ApplyQ( RIGHT, NORMAL, G, t, A );
    }

    // Return || E21 ||1 / || A ||1 and the chosen rank
    auto part = ComputePartition( A );
    part.value /= oneA;
    return part;
}

template<typename F>
inline ValueInt<BASE(F)>
SignDivide( DistMatrix<F>& A, DistMatrix<F>& G, bool returnQ=false )
{
#ifndef RELEASE
    CallStackEntry cse("schur::SignDivide");
#endif
    typedef BASE(F) Real;
    const Grid& g = A.Grid();
    const Int n = A.Height();

    // G := sgn(G)
    // G := 1/2 ( G + I )
    Sign( G );
    UpdateDiagonal( G, F(1) );
    Scale( F(1)/F(2), G );

    // Compute the pivoted QR decomposition of the spectral projection 
    DistMatrix<F,MD,STAR> t(g);
    DistMatrix<Int,VR,STAR> p(g);
    elem::QR( G, t, p );

    // A := Q^H A Q
    const Real oneA = OneNorm( A );
    if( returnQ )
    {
        ExpandPackedReflectors( LOWER, VERTICAL, UNCONJUGATED, 0, G, t );
        DistMatrix<F> B(g);
        Gemm( ADJOINT, NORMAL, F(1), G, A, B );
        Gemm( NORMAL, NORMAL, F(1), B, G, A );
    }
    else
    {
        qr::ApplyQ( LEFT, ADJOINT, G, t, A );
        qr::ApplyQ( RIGHT, NORMAL, G, t, A );
    }

    // Return || E21 ||1 / || A ||1 and the chosen rank
    auto part = ComputePartition( A );
    part.value /= oneA;
    return part;
}

template<typename F>
inline ValueInt<BASE(F)>
RandomizedSignDivide
( Matrix<F>& A, Matrix<F>& G, 
  bool returnQ=false, Int maxIts=1, BASE(F) relTol=0 )
{
#ifndef RELEASE
    CallStackEntry cse("schur::RandomizedSignDivide");
#endif
    typedef BASE(F) Real;
    const Int n = A.Height();
    const Real oneA = OneNorm( A );
    if( relTol == Real(0) )
        relTol = 500*n*lapack::MachineEpsilon<Real>();

    // S := sgn(G)
    // S := 1/2 ( S + I )
    auto S( G );
    Sign( S );
    UpdateDiagonal( S, F(1) );
    Scale( F(1)/F(2), S );

    ValueInt<Real> part;
    Matrix<F> V, B, t;
    Int it=0;
    while( it < maxIts )
    {
        G = S;

        // Compute the RURV of the spectral projector
        ImplicitHaar( V, t, n );
        qr::ApplyQ( RIGHT, NORMAL, V, t, G );
        elem::QR( G, t );

        // A := Q^H A Q [and reuse space for V for keeping original A]
        V = A;
        if( returnQ )
        {
            ExpandPackedReflectors( LOWER, VERTICAL, UNCONJUGATED, 0, G, t );
            Gemm( ADJOINT, NORMAL, F(1), G, A, B );
            Gemm( NORMAL, NORMAL, F(1), B, G, A );
        }
        else
        {
            qr::ApplyQ( LEFT, ADJOINT, G, t, A );
            qr::ApplyQ( RIGHT, NORMAL, G, t, A );
        }

        // || E21 ||1 / || A ||1 and the chosen rank
        part = ComputePartition( A );
        part.value /= oneA;

        ++it;
        if( part.value <= relTol || it == maxIts )
            break;
        else
            A = V;
    }
    return part;
}

template<typename F>
inline ValueInt<BASE(F)>
RandomizedSignDivide
( DistMatrix<F>& A, DistMatrix<F>& G, 
  bool returnQ=false, Int maxIts=1, BASE(F) relTol=0 )
{
#ifndef RELEASE
    CallStackEntry cse("schur::RandomizedSignDivide");
#endif
    typedef BASE(F) Real;
    const Grid& g = A.Grid();
    const Int n = A.Height();
    const Real oneA = OneNorm( A );
    if( relTol == Real(0) )
        relTol = 500*n*lapack::MachineEpsilon<Real>();

    // S := sgn(G)
    // S := 1/2 ( S + I )
    auto S( G );
    Sign( S );
    UpdateDiagonal( S, F(1) );
    Scale( F(1)/F(2), S );

    ValueInt<Real> part;
    DistMatrix<F> V(g), B(g);
    DistMatrix<F,MD,STAR> t(g);
    Int it=0;
    while( it < maxIts )
    {
        G = S;

        // Compute the RURV of the spectral projector
        ImplicitHaar( V, t, n );
        qr::ApplyQ( RIGHT, NORMAL, V, t, G );
        elem::QR( G, t );

        // A := Q^H A Q [and reuse space for V for keeping original A]
        V = A;
        if( returnQ )
        {
            ExpandPackedReflectors( LOWER, VERTICAL, UNCONJUGATED, 0, G, t );
            Gemm( ADJOINT, NORMAL, F(1), G, A, B );
            Gemm( NORMAL, NORMAL, F(1), B, G, A );
        }
        else
        {
            qr::ApplyQ( LEFT, ADJOINT, G, t, A );
            qr::ApplyQ( RIGHT, NORMAL, G, t, A );
        }

        // || E21 ||1 / || A ||1 and the chosen rank
        part = ComputePartition( A );
        part.value /= oneA;

        ++it;
        if( part.value <= relTol || it == maxIts )
            break;
        else
            A = V;
    }
    return part;
}

template<typename Real>
inline ValueInt<Real>
SpectralDivide
( Matrix<Real>& A, Int maxInnerIts=1, Int maxOuterIts=10, Real relTol=0 )
{
#ifndef RELEASE
    CallStackEntry cse("schur::SpectralDivide");
#endif
    const Int n = A.Height();
    const ValueInt<Real> median = Median(A.GetDiagonal());
    const Real infNorm = InfinityNorm(A);
    const Real eps = lapack::MachineEpsilon<Real>();
    if( relTol == Real(0) )
        relTol = 500*n*eps;
    const Real spread = 100*eps*infNorm;

    Int it=0;
    ValueInt<Real> part;
    Matrix<Real> G, ACopy;
    if( maxOuterIts > 1 )
        ACopy = A;
    while( it < maxOuterIts )
    {
        const Real shift = SampleBall<Real>(-median.value,spread);

        G = A;
        UpdateDiagonal( G, shift );

        //part = SignDivide( A, G );
        part = RandomizedSignDivide( A, G, false, maxInnerIts, relTol );

        ++it;
        if( part.value <= relTol )
            break;
        else if( it != maxOuterIts )
            A = ACopy;
    }
    if( part.value > relTol )
    {
        std::ostringstream os;
        os << "Unable to split spectrum to specified accuracy: part.value="
           << part.value << ", relTol=" << relTol << std::endl;
        RuntimeError( os.str() );
    }

    return part;
}

template<typename Real>
inline ValueInt<Real>
SpectralDivide
( Matrix<Complex<Real> >& A, 
  Int maxInnerIts=1, Int maxOuterIts=10, Real relTol=0 )
{
#ifndef RELEASE
    CallStackEntry cse("schur::SpectralDivide");
#endif
    typedef Complex<Real> F;
    const Int n = A.Height();
    const Real infNorm = InfinityNorm(A);
    const Real eps = lapack::MachineEpsilon<Real>();
    if( relTol == Real(0) )
        relTol = 500*n*eps;
    const Real spread = 100*eps*infNorm;

    Int it=0;
    ValueInt<Real> part;
    Matrix<F> G, ACopy;
    if( maxOuterIts > 1 )
        ACopy = A;
    while( it < maxOuterIts )
    {
        const Real angle = Uniform<Real>(0,2*Pi);
        const F gamma = F(Cos(angle),Sin(angle));
        G = A;
        Scale( gamma, G );

        const auto median = Median(G.GetRealPartOfDiagonal());
        const F shift = SampleBall<F>(-median.value,spread);
        UpdateDiagonal( G, shift );

        //part = SignDivide( A, G );
        part = RandomizedSignDivide( A, G, false, maxInnerIts, relTol );

        ++it;
        if( part.value <= relTol )
            break;
        else if( it != maxOuterIts )
            A = ACopy;
    }
    if( part.value > relTol )
    {
        std::ostringstream os;
        os << "Unable to split spectrum to specified accuracy: part.value="
           << part.value << ", relTol=" << relTol << std::endl;
        RuntimeError( os.str() );
    }

    return part;
}

template<typename Real>
inline ValueInt<Real>
SpectralDivide
( Matrix<Real>& A, Matrix<Real>& Q, 
  Int maxInnerIts=1, Int maxOuterIts=10, Real relTol=0 )
{
#ifndef RELEASE
    CallStackEntry cse("schur::SpectralDivide");
#endif
    const Int n = A.Height();
    const auto median = Median(A.GetDiagonal());
    const Real infNorm = InfinityNorm(A);
    const Real eps = lapack::MachineEpsilon<Real>();
    if( relTol == Real(0) )
        relTol = 500*n*eps;
    const Real spread = 100*eps*infNorm;

    Int it=0;
    ValueInt<Real> part;
    Matrix<Real> ACopy;
    if( maxOuterIts > 1 )
        ACopy = A;
    while( it < maxOuterIts )
    {
        const Real shift = SampleBall<Real>(-median.value,spread);

        Q = A;
        UpdateDiagonal( Q, shift );

        //part = SignDivide( A, Q, true );
        part = RandomizedSignDivide( A, Q, true, maxInnerIts, relTol );

        ++it;
        if( part.value <= relTol )
            break;
        else if( it != maxOuterIts )
            A = ACopy;
    }
    if( part.value > relTol )
    {
        std::ostringstream os;
        os << "Unable to split spectrum to specified accuracy: part.value="
           << part.value << ", relTol=" << relTol << std::endl;
        RuntimeError( os.str() );
    }

    return part;
}

template<typename Real>
inline ValueInt<Real>
SpectralDivide
( Matrix<Complex<Real> >& A, Matrix<Complex<Real> >& Q, 
  Int maxInnerIts=1, Int maxOuterIts=10, Real relTol=0 )
{
#ifndef RELEASE
    CallStackEntry cse("schur::SpectralDivide");
#endif
    typedef Complex<Real> F;
    const Int n = A.Height();
    const Real infNorm = InfinityNorm(A);
    const Real eps = lapack::MachineEpsilon<Real>();
    if( relTol == Real(0) )
        relTol = 500*n*eps;
    const Real spread = 100*eps*infNorm;

    Int it=0;
    ValueInt<Real> part;
    Matrix<F> ACopy;
    if( maxOuterIts > 1 )
        ACopy = A;
    while( it < maxOuterIts )
    {
        const Real angle = Uniform<Real>(0,2*Pi);
        const F gamma = F(Cos(angle),Sin(angle));
        Q = A;
        Scale( gamma, Q );

        const auto median = Median(Q.GetRealPartOfDiagonal());
        const F shift = SampleBall<F>(-median.value,spread);
        UpdateDiagonal( Q, shift );

        //part = SignDivide( A, Q, true );
        part = RandomizedSignDivide( A, Q, true, maxInnerIts, relTol );

        ++it;
        if( part.value <= relTol )
            break;
        else if( it != maxOuterIts )
            A = ACopy;
    }
    if( part.value > relTol )
    {
        std::ostringstream os;
        os << "Unable to split spectrum to specified accuracy: part.value="
           << part.value << ", relTol=" << relTol << std::endl;
        RuntimeError( os.str() );
    }

    return part;
}

template<typename Real>
inline ValueInt<Real>
SpectralDivide
( DistMatrix<Real>& A, Int maxInnerIts=1, Int maxOuterIts=10, Real relTol=0 )
{
#ifndef RELEASE
    CallStackEntry cse("schur::SpectralDivide");
#endif
    const Int n = A.Height();
    const auto median = Median(A.GetDiagonal());
    const Real infNorm = InfinityNorm(A);
    const Real eps = lapack::MachineEpsilon<Real>();
    if( relTol == Real(0) )
        relTol = 500*n*eps;
    const Real spread = 100*eps*infNorm;

    Int it=0;
    ValueInt<Real> part;
    DistMatrix<Real> ACopy(A.Grid()), G(A.Grid());
    if( maxOuterIts > 1 )
        ACopy = A;
    while( it < maxOuterIts )
    {
        Real shift = SampleBall<Real>(-median.value,spread);
        mpi::Broadcast( shift, 0, A.Grid().VCComm() );

        G = A;
        UpdateDiagonal( G, shift );

        //part = SignDivide( A, G );
        part = RandomizedSignDivide( A, G, false, maxInnerIts, relTol );

        ++it;
        if( part.value <= relTol )
            break;
        else if( it != maxOuterIts )
            A = ACopy;
    }
    if( part.value > relTol )
    {
        std::ostringstream os;
        os << "Unable to split spectrum to specified accuracy: part.value="
           << part.value << ", relTol=" << relTol << std::endl;
        RuntimeError( os.str() );
    }

    return part;
}

template<typename Real>
inline ValueInt<Real>
SpectralDivide
( DistMatrix<Complex<Real> >& A, 
  Int maxInnerIts=1, Int maxOuterIts=10, Real relTol=0 )
{
#ifndef RELEASE
    CallStackEntry cse("schur::SpectralDivide");
#endif
    typedef Complex<Real> F;
    const Int n = A.Height();
    const Real infNorm = InfinityNorm(A);
    const Real eps = lapack::MachineEpsilon<Real>();
    if( relTol == Real(0) )
        relTol = 500*n*eps;
    const Real spread = 100*eps*infNorm;

    Int it=0;
    ValueInt<Real> part;
    DistMatrix<F> ACopy(A.Grid()), G(A.Grid());
    if( maxOuterIts > 1 )
        ACopy = A;
    while( it < maxOuterIts )
    {
        const Real angle = Uniform<Real>(0,2*Pi);
        F gamma = F(Cos(angle),Sin(angle));
        mpi::Broadcast( gamma, 0, A.Grid().VCComm() );
        G = A;
        Scale( gamma, G );

        const auto median = Median(G.GetRealPartOfDiagonal());
        F shift = SampleBall<F>(-median.value,spread);
        mpi::Broadcast( shift, 0, A.Grid().VCComm() );
        UpdateDiagonal( G, shift );

        //part = SignDivide( A, G );
        part = RandomizedSignDivide( A, G, false, maxInnerIts, relTol );

        ++it;
        if( part.value <= relTol )
            break;
        else if( it != maxOuterIts )
            A = ACopy;
    }
    if( part.value > relTol )
    {
        std::ostringstream os;
        os << "Unable to split spectrum to specified accuracy: part.value="
           << part.value << ", relTol=" << relTol << std::endl;
        RuntimeError( os.str() );
    }

    return part;
}

template<typename Real>
inline ValueInt<Real>
SpectralDivide
( DistMatrix<Real>& A, DistMatrix<Real>& Q, 
  Int maxInnerIts=1, Int maxOuterIts=10, Real relTol=0 )
{
#ifndef RELEASE
    CallStackEntry cse("schur::SpectralDivide");
#endif
    const Int n = A.Height();
    const Real infNorm = InfinityNorm(A);
    const auto median = Median(A.GetDiagonal());
    const Real eps = lapack::MachineEpsilon<Real>();
    if( relTol == Real(0) )
        relTol = 500*n*eps;
    const Real spread = 100*eps*infNorm;

    Int it=0;
    ValueInt<Real> part;
    DistMatrix<Real> ACopy(A.Grid());
    if( maxOuterIts > 1 )
        ACopy = A;
    while( it < maxOuterIts )
    {
        Real shift = SampleBall<Real>(-median.value,spread);
        mpi::Broadcast( shift, 0, A.Grid().VCComm() );

        Q = A;
        UpdateDiagonal( Q, shift );

        //part = SignDivide( A, Q, true );
        part = RandomizedSignDivide( A, Q, true, maxInnerIts, relTol );

        ++it;
        if( part.value <= relTol )
            break;
        else if( it != maxOuterIts )
            A = ACopy;
    }
    if( part.value > relTol )
    {
        std::ostringstream os;
        os << "Unable to split spectrum to specified accuracy: part.value="
           << part.value << ", relTol=" << relTol << std::endl;
        RuntimeError( os.str() );
    }

    return part;
}

template<typename Real>
inline ValueInt<Real>
SpectralDivide
( DistMatrix<Complex<Real> >& A, DistMatrix<Complex<Real> >& Q,
  Int maxInnerIts=1, Int maxOuterIts=10, Real relTol=0 )
{
#ifndef RELEASE
    CallStackEntry cse("schur::SpectralDivide");
#endif
    typedef Complex<Real> F;
    const Int n = A.Height();
    const Real infNorm = InfinityNorm(A);
    const Real eps = lapack::MachineEpsilon<Real>();
    if( relTol == Real(0) )
        relTol = 500*n*eps;
    const Real spread = 100*eps*infNorm;

    Int it=0;
    ValueInt<Real> part;
    DistMatrix<F> ACopy(A.Grid());
    if( maxOuterIts > 1 )
        ACopy = A;
    while( it < maxOuterIts )
    {
        const Real angle = Uniform<Real>(0,2*Pi);
        F gamma = F(Cos(angle),Sin(angle));
        mpi::Broadcast( gamma, 0, A.Grid().VCComm() );
        Q = A;
        Scale( gamma, Q );

        const auto median = Median(Q.GetRealPartOfDiagonal());
        F shift = SampleBall<F>(-median.value,spread);
        mpi::Broadcast( shift, 0, A.Grid().VCComm() );
        UpdateDiagonal( Q, shift );

        //part = SignDivide( A, Q, true );
        part = RandomizedSignDivide( A, Q, true, maxInnerIts, relTol );

        ++it;
        if( part.value <= relTol )
            break;
        else if( it != maxOuterIts )
            A = ACopy;
    }
    if( part.value > relTol )
    {
        std::ostringstream os;
        os << "Unable to split spectrum to specified accuracy: part.value="
           << part.value << ", relTol=" << relTol << std::endl;
        RuntimeError( os.str() );
    }

    return part;
}

template<typename F>
inline void
SDC
( Matrix<F>& A, Matrix<COMPLEX(F)>& w, Int cutoff=256, 
  Int maxInnerIts=1, Int maxOuterIts=10, BASE(F) relTol=0 )
{
#ifndef RELEASE
    CallStackEntry cse("schur::SDC");
#endif
    typedef BASE(F) Real;
    const Int n = A.Height();
    w.ResizeTo( n, 1 );
    if( n <= cutoff )
    {
        schur::QR( A, w );
        return;
    }

    // Perform this level's split
    const auto part = SpectralDivide( A, maxInnerIts, maxOuterIts, relTol );
    Matrix<F> ATL, ATR,
              ABL, ABR;
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, part.index );
    Matrix<COMPLEX(F)> wT, wB;
    PartitionDown( w, wT, wB, part.index );

    // Recurse on the two subproblems
    SDC( ATL, wT, cutoff, maxInnerIts, maxOuterIts, relTol );
    SDC( ABR, wB, cutoff, maxInnerIts, maxOuterIts, relTol );
}

template<typename F>
inline void
SDC
( Matrix<F>& A, Matrix<COMPLEX(F)>& w, Matrix<F>& Q, 
  bool formATR=true, Int cutoff=256, Int maxInnerIts=1, Int maxOuterIts=10, 
  BASE(F) relTol=0 )
{
#ifndef RELEASE
    CallStackEntry cse("schur::SDC");
#endif
    typedef BASE(F) Real;
    const Int n = A.Height();
    w.ResizeTo( n, 1 );
    Q.ResizeTo( n, n );
    if( n <= cutoff )
    {
        schur::QR( A, w, Q, formATR );
        return;
    }

    // Perform this level's split
    const auto part = SpectralDivide( A, Q, maxInnerIts, maxOuterIts, relTol );
    Matrix<F> ATL, ATR,
              ABL, ABR;
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, part.index );
    Matrix<COMPLEX(F)> wT, wB;
    PartitionDown( w, wT, wB, part.index );
    Matrix<F> QL, QR;
    PartitionRight( Q, QL, QR, part.index );

    // Recurse on the top-left quadrant and update Schur vectors and ATR
    Matrix<F> Z;
    SDC( ATL, wT, Z, formATR, cutoff, maxInnerIts, maxOuterIts, relTol );
    auto G( QL );
    Gemm( NORMAL, NORMAL, F(1), G, Z, QL );
    if( formATR )
        Gemm( ADJOINT, NORMAL, F(1), Z, ATR, G );

    // Recurse on the bottom-right quadrant and update Schur vectors and ATR
    SDC( ABR, wB, Z, formATR, cutoff, maxInnerIts, maxOuterIts, relTol );
    if( formATR )
        Gemm( NORMAL, NORMAL, F(1), G, Z, ATR ); 
    G = QR;
    Gemm( NORMAL, NORMAL, F(1), G, Z, QR );
}

template<typename F>
inline void PushSubproblems
( DistMatrix<F>& ATL,    DistMatrix<F>& ABR, 
  DistMatrix<F>& ATLSub, DistMatrix<F>& ABRSub,
  DistMatrix<COMPLEX(F),VR,STAR>& wT,    DistMatrix<COMPLEX(F),VR,STAR>& wB,
  DistMatrix<COMPLEX(F),VR,STAR>& wTSub, DistMatrix<COMPLEX(F),VR,STAR>& wBSub )
{
#ifndef RELEASE
    CallStackEntry cse("schur::PushSubproblems");
#endif
    // The trivial push
    /*
    ATLSub = View( ATL );
    ABRSub = View( ABR );
    wTSub = View( wT );
    wBSub = View( wB );
    */

    // Split based on the work estimates
    typedef BASE(F) Real;
    const Int nLeft = ATL.Height();
    const Int nRight = ABR.Height();
    const Real leftWork = Pow(Real(nLeft),Real(3));
    const Real rightWork = Pow(Real(nRight),Real(3));
    const Real ratio = leftWork / (leftWork+rightWork);
    const Grid& g = ATL.Grid();
    // It should be guaranteed that p >= 2
    const Int p = g.Size();
    const Int pLeftProp = round(ratio*p);
    const Int pLeft = std::max(1,std::min(p-1,pLeftProp));
    const Int pRight = g.Size() - pLeft;

    std::vector<int> leftRanks(pLeft), rightRanks(pRight);
    for( int j=0; j<pLeft; ++j )
        leftRanks[j] = j;
    for( int j=0; j<pRight; ++j )
        rightRanks[j] = j+pLeft;
    mpi::Group group = g.OwningGroup();
    mpi::Group leftGroup, rightGroup;
    mpi::GroupIncl( group, pLeft, &leftRanks[0], leftGroup );
    mpi::GroupIncl( group, pRight, &rightRanks[0], rightGroup );
    Grid *leftGrid, *rightGrid;
    leftGrid = new Grid( g.VCComm(), leftGroup, Grid::FindFactor(pLeft) );
    rightGrid = new Grid( g.VCComm(), rightGroup, Grid::FindFactor(pRight) );

    ATLSub.SetGrid( *leftGrid ); 
    ABRSub.SetGrid( *rightGrid );
    wTSub.SetGrid( *leftGrid );
    wBSub.SetGrid( *rightGrid );
    ATLSub = ATL;
    ABRSub = ABR;
}

template<typename F>
inline void PullSubproblems
( DistMatrix<F>& ATL,    DistMatrix<F>& ABR,
  DistMatrix<F>& ATLSub, DistMatrix<F>& ABRSub,
  DistMatrix<COMPLEX(F),VR,STAR>& wT,    DistMatrix<COMPLEX(F),VR,STAR>& wB,
  DistMatrix<COMPLEX(F),VR,STAR>& wTSub, DistMatrix<COMPLEX(F),VR,STAR>& wBSub )
{
#ifndef RELEASE
    CallStackEntry cse("schur::PullSubproblems");
#endif
    // The trivial pull is empty
    ATL = ATLSub;
    ABR = ABRSub;
    // This is a hack
    //wT = wTSub;
    //wB = wBSub;
    {
        DistMatrix<COMPLEX(F)> wTSub_MC_MR( wTSub.Grid() );
        if( wTSub.Participating() )
            wTSub_MC_MR = wTSub;
        wTSub_MC_MR.MakeConsistent();
        DistMatrix<COMPLEX(F)> wT_MC_MR(wT.Grid()); 
        wT_MC_MR = wTSub_MC_MR;
        wT = wT_MC_MR;
    }
    {
        DistMatrix<COMPLEX(F)> wBSub_MC_MR( wBSub.Grid() );
        if( wBSub.Participating() )
            wBSub_MC_MR = wBSub;
        wBSub_MC_MR.MakeConsistent();
        DistMatrix<COMPLEX(F)> wB_MC_MR(wB.Grid());
        wB_MC_MR = wBSub_MC_MR;
        wB = wB_MC_MR;
    }
    
    const Grid *leftGrid = &ATLSub.Grid();
    const Grid *rightGrid = &ABRSub.Grid();
    ATLSub.Empty();
    ABRSub.Empty();
    wTSub.Empty();
    wBSub.Empty();
    mpi::Group leftOwning = leftGrid->OwningGroup();
    mpi::Group rightOwning = rightGrid->OwningGroup();
    delete leftGrid;
    delete rightGrid;
    mpi::GroupFree( leftOwning );
    mpi::GroupFree( rightOwning );
}

template<typename F>
inline void
SDC
( DistMatrix<F>& A, DistMatrix<COMPLEX(F)>& w, Int cutoff=256, 
  Int maxInnerIts=1, Int maxOuterIts=10, BASE(F) relTol=0 )
{
#ifndef RELEASE
    CallStackEntry cse("schur::SDC");
#endif
    typedef BASE(F) Real;
    const Grid& g = A.Grid();
    const Int n = A.Height();
    w.ResizeTo( n, 1 );
    if( A.Grid().Size() == 1 )
    {
        schur::QR( A.Matrix(), w.Matrix() );
        return;
    }
    if( n <= cutoff )
    {
        DistMatrix<F,CIRC,CIRC> A_CIRC_CIRC( A );
        DistMatrix<COMPLEX(F),CIRC,CIRC> w_CIRC_CIRC( w );
        if( g.VCRank() == A_CIRC_CIRC.Root() )
            schur::QR( A_CIRC_CIRC.Matrix(), w_CIRC_CIRC.Matrix() );
        A = A_CIRC_CIRC;
        w = w_CIRC_CIRC;
        return;
    }

    // Perform this level's split
    const auto part = SpectralDivide( A, maxInnerIts, maxOuterIts, relTol );
    DistMatrix<F> ATL(g), ATR(g),
                  ABL(g), ABR(g);
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, part.index );
    DistMatrix<COMPLEX(F),VR,STAR> wT(g), wB(g);
    PartitionDown( w, wT, wB, part.index );

    DistMatrix<F> ATLSub, ABRSub;
    DistMatrix<COMPLEX(F),VR,STAR> wTSub, wBSub;
    PushSubproblems( ATL, ABR, ATLSub, ABRSub, wT, wB, wTSub, wBSub );
    if( ATLSub.Participating() )
        SDC( ATLSub, wTSub, cutoff, maxInnerIts, maxOuterIts, relTol );
    if( ABRSub.Participating() )
        SDC( ABRSub, wBSub, cutoff, maxInnerIts, maxOuterIts, relTol );
    PullSubproblems( ATL, ABR, ATLSub, ABRSub, wT, wB, wTSub, wBSub );
}

template<typename F>
inline void PushSubproblems
( DistMatrix<F>& ATL,    DistMatrix<F>& ABR, 
  DistMatrix<F>& ATLSub, DistMatrix<F>& ABRSub,
  DistMatrix<COMPLEX(F),VR,STAR>& wT,    DistMatrix<COMPLEX(F),VR,STAR>& wB,
  DistMatrix<COMPLEX(F),VR,STAR>& wTSub, DistMatrix<COMPLEX(F),VR,STAR>& wBSub,
  DistMatrix<F>& ZTSub,  DistMatrix<F>& ZBSub )
{
#ifndef RELEASE
    CallStackEntry cse("schur::PushSubproblems");
#endif
    // The trivial push
    /*
    ATLSub = View( ATL );
    ABRSub = View( ABR );
    wTSub = View( wT );
    wBSub = View( wB );
    ZTSub.SetGrid( ATL.Grid() );
    ZBSub.SetGrid( ABR.Grid() );
    */

    // Split based on the work estimates
    typedef BASE(F) Real;
    const Int nLeft = ATL.Height();
    const Int nRight = ABR.Height();
    const Real leftWork = Pow(Real(nLeft),Real(3));
    const Real rightWork = Pow(Real(nRight),Real(3));
    const Real ratio = leftWork / (leftWork+rightWork);
    const Grid& g = ATL.Grid();
    // It should be guaranteed that p >= 2
    const Int p = g.Size();
    const Int pLeftProp = round(ratio*p);
    const Int pLeft = std::max(1,std::min(p-1,pLeftProp));
    const Int pRight = g.Size() - pLeft;

    std::vector<int> leftRanks(pLeft), rightRanks(pRight);
    for( int j=0; j<pLeft; ++j )
        leftRanks[j] = j;
    for( int j=0; j<pRight; ++j )
        rightRanks[j] = j+pLeft;
    mpi::Group group = g.OwningGroup();
    mpi::Group leftGroup, rightGroup;
    mpi::GroupIncl( group, pLeft, &leftRanks[0], leftGroup );
    mpi::GroupIncl( group, pRight, &rightRanks[0], rightGroup );
    Grid *leftGrid, *rightGrid;
    leftGrid = new Grid( g.VCComm(), leftGroup, Grid::FindFactor(pLeft) );
    rightGrid = new Grid( g.VCComm(), rightGroup, Grid::FindFactor(pRight) );

    ATLSub.SetGrid( *leftGrid );
    ABRSub.SetGrid( *rightGrid );
    wTSub.SetGrid( *leftGrid );
    wBSub.SetGrid( *rightGrid );
    ZTSub.SetGrid( *leftGrid );
    ZBSub.SetGrid( *rightGrid );
    ATLSub = ATL;
    ABRSub = ABR;
}

template<typename F>
inline void PullSubproblems
( DistMatrix<F>& ATL,    DistMatrix<F>& ABR,
  DistMatrix<F>& ATLSub, DistMatrix<F>& ABRSub,
  DistMatrix<COMPLEX(F),VR,STAR>& wT,    DistMatrix<COMPLEX(F),VR,STAR>& wB,
  DistMatrix<COMPLEX(F),VR,STAR>& wTSub, DistMatrix<COMPLEX(F),VR,STAR>& wBSub,
  DistMatrix<F>& ZT,     DistMatrix<F>& ZB,
  DistMatrix<F>& ZTSub,  DistMatrix<F>& ZBSub )
{
#ifndef RELEASE
    CallStackEntry cse("schur::PullSubproblems");
#endif
    // The trivial pull
    /*
    ZT = View( ZTSub );
    ZB = View( ZBSub );
    */

    ATL = ATLSub;
    ABR = ABRSub;
    // This is a hack
    //wT = wTSub;
    //wB = wBSub;
    {
        DistMatrix<COMPLEX(F)> wTSub_MC_MR( wTSub.Grid() );
        if( wTSub.Participating() )
            wTSub_MC_MR = wTSub;
        wTSub_MC_MR.MakeConsistent();
        DistMatrix<COMPLEX(F)> wT_MC_MR(wT.Grid()); 
        wT_MC_MR = wTSub_MC_MR;
        wT = wT_MC_MR;
    }
    {
        DistMatrix<COMPLEX(F)> wBSub_MC_MR( wBSub.Grid() );
        if( wBSub.Participating() )
            wBSub_MC_MR = wBSub;
        wBSub_MC_MR.MakeConsistent();
        DistMatrix<COMPLEX(F)> wB_MC_MR(wB.Grid()); 
        wB_MC_MR = wBSub_MC_MR;
        wB = wB_MC_MR;
    }
    ZTSub.MakeConsistent();
    ZBSub.MakeConsistent();
    ZT = ZTSub;
    ZB = ZBSub;
    const Grid *leftGrid = &ATLSub.Grid();
    const Grid *rightGrid = &ABRSub.Grid();
    ATLSub.Empty();
    ABRSub.Empty();
    wTSub.Empty();
    wBSub.Empty();
    ZTSub.Empty();
    ZBSub.Empty();
    mpi::Group leftOwning = leftGrid->OwningGroup();
    mpi::Group rightOwning = rightGrid->OwningGroup();
    delete leftGrid;
    delete rightGrid;
    mpi::GroupFree( leftOwning );
    mpi::GroupFree( rightOwning );
}

template<typename F>
inline void
SDC
( DistMatrix<F>& A, DistMatrix<COMPLEX(F),VR,STAR>& w, DistMatrix<F>& Q, 
  bool formATR=true, Int cutoff=256, Int maxInnerIts=1, Int maxOuterIts=10, 
  BASE(F) relTol=0 )
{
#ifndef RELEASE
    CallStackEntry cse("schur::SDC");
#endif
    typedef BASE(F) Real;
    const Grid& g = A.Grid();
    const Int n = A.Height();
    w.ResizeTo( n, 1 );
    Q.ResizeTo( n, n );
    if( A.Grid().Size() == 1 )
    {
        schur::QR( A.Matrix(), w.Matrix(), Q.Matrix(), formATR );
        return;
    }
    if( n <= cutoff )
    {
        DistMatrix<F,CIRC,CIRC> A_CIRC_CIRC( A ), Q_CIRC_CIRC( n, n, g );
        DistMatrix<COMPLEX(F),CIRC,CIRC> w_CIRC_CIRC( n, 1, g );
        if( g.VCRank() == A_CIRC_CIRC.Root() )
            schur::QR
            ( A_CIRC_CIRC.Matrix(), w_CIRC_CIRC.Matrix(), Q_CIRC_CIRC.Matrix(),
              formATR );
        A = A_CIRC_CIRC;
        w = w_CIRC_CIRC;
        Q = Q_CIRC_CIRC;
        return;
    }

    // Perform this level's split
    const Real infNorm = InfinityNorm( A );
    const auto part = SpectralDivide( A, Q, maxInnerIts, maxOuterIts, relTol );
    DistMatrix<F> ATL(g), ATR(g),
                  ABL(g), ABR(g);
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, part.index );
    DistMatrix<COMPLEX(F),VR,STAR> wT(g), wB(g);
    PartitionDown( w, wT, wB, part.index );
    DistMatrix<F> QL(g), QR(g);
    PartitionRight( Q, QL, QR, part.index );

    // Recurse on the two subproblems
    DistMatrix<F> ATLSub, ABRSub, ZTSub, ZBSub;
    DistMatrix<COMPLEX(F),VR,STAR> wTSub, wBSub;
    PushSubproblems
    ( ATL, ABR, ATLSub, ABRSub, wT, wB, wTSub, wBSub, ZTSub, ZBSub );
    if( ATLSub.Participating() )
        SDC( ATLSub, wTSub, ZTSub, 
             formATR, cutoff, maxInnerIts, maxOuterIts, relTol );
    if( ABRSub.Participating() )
        SDC( ABRSub, wBSub, ZBSub, 
             formATR, cutoff, maxInnerIts, maxOuterIts, relTol );
    
    // Ensure that the results are back on this level's grid
    DistMatrix<F> ZT(g), ZB(g);
    PullSubproblems
    ( ATL, ABR, ATLSub, ABRSub, wT, wB, wTSub, wBSub, ZT, ZB, ZTSub, ZBSub );

    // Update the Schur vectors
    auto G( QL );
    Gemm( NORMAL, NORMAL, F(1), G, ZT, QL );
    G = QR;
    Gemm( NORMAL, NORMAL, F(1), G, ZB, QR );

    if( formATR )
    {
        // Update the top-right quadrant
        Gemm( ADJOINT, NORMAL, F(1), ZT, ATR, G );
        Gemm( NORMAL, NORMAL, F(1), G, ZB, ATR ); 
    }
}

} // namespace schur
} // namespace elem

#endif // ifndef ELEM_LAPACK_SCHUR_SDC_HPP
