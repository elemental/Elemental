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
#include "elemental/lapack-like/Norm/One.hpp"
#include "elemental/lapack-like/QR.hpp"
#include "elemental/lapack-like/Sign.hpp"
#include "elemental/matrices/Haar.hpp"
#include "elemental/matrices/Identity.hpp"

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
    PushCallStack("schur::ComputePartition");
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
    part.index = 0;
    for( Int j=1; j<n-1; ++j )
    {
        norms[j] = norms[j-1] + colSums[j] - rowSums[j-1];
        if( norms[j] < part.value )
        {
            part.value = norms[j];
            part.index = j;
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
    PushCallStack("schur::ComputePartition");
#endif
    typedef BASE(F) Real;
    const Grid& g = A.Grid();
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
    part.index = 0;
    for( Int j=1; j<n-1; ++j )
    {
        norms[j] = norms[j-1] + colSums[j] - rowSums[j-1];
        if( norms[j] < part.value )
        {
            part.value = norms[j];
            part.index = j;
        }
    }

    return part;
}

template<typename F>
inline ValueInt<BASE(F)>
SpectralDivide( Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry cse("schur::SpectralDivide");
#endif
    typedef BASE(F) Real;
    const Int n = A.Height();

    // S := sgn(A)
    Matrix<F> S( A );
    Sign( S );

    // Compute the spectral projector, B := 1/2 ( S + I ), and its trace
    Matrix<F> B; 
    Identity( B, n, n );
    Axpy( F(1), S, B );
    Scale( F(1)/F(2), B );

    // Compute the pivoted QR decomposition of the spectral projection
    Matrix<F> t;
    Matrix<Int> p;
    QR( B, t, p );

    // A := Q^H A Q
    const Real oneA = OneNorm( A );
    qr::ApplyQ( LEFT, ADJOINT, B, t, A );
    qr::ApplyQ( RIGHT, NORMAL, B, t, A );

    // Return || E21 ||1 / || A ||1 and the chosen rank
    ValueInt<Real> part = ComputePartition( A );
    part.value /= oneA;
    return part;
}

template<typename F>
inline ValueInt<BASE(F)>
SpectralDivide( DistMatrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry cse("schur::SpectralDivide");
#endif
    typedef BASE(F) Real;
    const Int n = A.Height();
    const Grid& g = A.Grid();

    // S := sgn(A)
    DistMatrix<F> S( A );
    Sign( S );

    // Compute the spectral projector, B := 1/2 ( S + I ), and its trace
    DistMatrix<F> B(g);
    Identity( B, n, n );
    Axpy( F(1), S, B );
    Scale( F(1)/F(2), B );

    // Compute the pivoted QR decomposition of the spectral projection
    DistMatrix<F,MD,STAR> t(g);
    DistMatrix<Int,VR,STAR> p(g);
    QR( B, t, p );

    // A := Q^H A Q
    const Real oneA = OneNorm( A );
    qr::ApplyQ( LEFT, ADJOINT, B, t, A );
    qr::ApplyQ( RIGHT, NORMAL, B, t, A );

    // Return || E21 ||1 / || A ||1 and the chosen rank
    ValueInt<Real> part = ComputePartition( A );
    part.value /= oneA;
    return part;
}

template<typename F>
inline ValueInt<BASE(F)>
RandomizedSpectralDivide( Matrix<F>& A, Int maxIts=10, BASE(F) relTol=0 )
{
#ifndef RELEASE
    CallStackEntry cse("schur::SpectralDivide");
#endif
    typedef BASE(F) Real;
    const Int n = A.Height();
    const Real oneA = OneNorm( A );
    if( relTol == Real(0) )
        relTol = n*lapack::MachineEpsilon<Real>();

    // S := sgn(A)
    Matrix<F> S( A );
    Sign( S );

    // Compute the spectral projector, B := 1/2 ( S + I ), and its trace
    Matrix<F> B; 
    Identity( B, n, n );
    Axpy( F(1), S, B );
    Scale( F(1)/F(2), B );

    ValueInt<Real> part;
    Matrix<F> V, t;
    for( Int it=0; it<maxIts; ++it )
    {
        // Compute the RURV of the spectral projector [and reuse S]
        ImplicitHaar( V, t, n );
        S = B;
        qr::ApplyQ( RIGHT, NORMAL, V, t, S );
        QR( S, t );

        // A := Q^H A Q [and reuse space for V for keeping original A]
        V = A;
        qr::ApplyQ( LEFT, ADJOINT, S, t, A );
        qr::ApplyQ( RIGHT, NORMAL, S, t, A );

        // || E21 ||1 / || A ||1 and the chosen rank
        part = ComputePartition( A );
        part.value /= oneA;

        if( part.value <= relTol )
            break;
        else
            A = V;
    }
    return part;
}

template<typename F>
inline ValueInt<BASE(F)>
RandomizedSpectralDivide( DistMatrix<F>& A, Int maxIts=10, BASE(F) relTol=0 )
{
#ifndef RELEASE
    CallStackEntry cse("schur::SpectralDivide");
#endif
    typedef BASE(F) Real;
    const Grid& g = A.Grid();
    const Int n = A.Height();
    const Real oneA = OneNorm( A );
    if( relTol == Real(0) )
        relTol = n*lapack::MachineEpsilon<Real>();

    // S := sgn(A)
    DistMatrix<F> S( A );
    Sign( S );

    // Compute the spectral projector, B := 1/2 ( S + I ), and its trace
    DistMatrix<F> B(g); 
    Identity( B, n, n );
    Axpy( F(1), S, B );
    Scale( F(1)/F(2), B );

    ValueInt<Real> part;
    DistMatrix<F> V(g);
    DistMatrix<F,MD,STAR> t(g);
    for( Int it=0; it<maxIts; ++it )
    {
        // Compute the RURV of the spectral projector [and reuse S]
        ImplicitHaar( V, t, n );
        S = B;
        qr::ApplyQ( RIGHT, NORMAL, V, t, S );
        QR( S, t );

        // A := Q^H A Q [and reuse space for V for keeping original A]
        V = A;
        qr::ApplyQ( LEFT, ADJOINT, S, t, A );
        qr::ApplyQ( RIGHT, NORMAL, S, t, A );

        // || E21 ||1 / || A ||1 and the chosen rank
        part = ComputePartition( A );
        part.value /= oneA;

        if( part.value <= relTol )
            break;
        else
            A = V;
    }
    return part;
}

} // namespace schur
} // namespace elem

#endif // ifndef ELEM_LAPACK_SCHUR_SDC_HPP
