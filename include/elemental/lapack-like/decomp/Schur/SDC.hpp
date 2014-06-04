/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_SCHUR_SDC_HPP
#define ELEM_SCHUR_SDC_HPP

#include ELEM_AXPY_INC
#include ELEM_SCALE_INC
#include ELEM_UPDATEDIAGONAL_INC

#include ELEM_MEDIAN_INC
#include ELEM_INFINITYNORM_INC
#include ELEM_ONENORM_INC
#include ELEM_QR_INC
#include ELEM_SIGN_INC

#include ELEM_HAAR_INC

// See Z. Bai, J. Demmel, J. Dongarra, A. Petitet, H. Robinson, and K. Stanley's
// "The spectral decomposition of nonsymmetric matrices on distributed memory
// parallel computers". Currently available at:
// www.netlib.org/lapack/lawnspdf/lawn91.pdf
//
// as well as the improved version, which avoids pivoted QR, in J. Demmel, 
// I. Dumitriu, and O. Holtz, "Fast linear algebra is stable", 2007.
// www.netlib.org/lapack/lawnspdf/lawn186.pdf

namespace elem {

// TODO: Create SignCtrl and use as a member
template<typename Real>
struct SdcCtrl {
    Int cutoff;
    Int maxInnerIts;
    Int maxOuterIts;
    Real tol;
    Real spreadFactor;
    bool random; 
    bool progress;

    SignCtrl<Real> signCtrl;

    SdcCtrl()
    : cutoff(256), maxInnerIts(2), maxOuterIts(10),
      tol(0), spreadFactor(1e-6),
      random(true), progress(false), signCtrl()
    { }
};

namespace schur {

template<typename F>
inline ValueInt<Base<F>>
ComputePartition( Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("schur::ComputePartition"))
    typedef Base<F> Real;
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
inline ValueInt<Base<F>>
ComputePartition( DistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("schur::ComputePartition"))
    typedef Base<F> Real;
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
    for( Int jLoc=0; jLoc<nLocal; ++jLoc )
    {
        const Int j = A.GlobalCol(jLoc);
        if( j < n-1 )
        {
            for( Int iLoc=0; iLoc<mLocal; ++iLoc )
            {
                const Int i = A.GlobalRow(iLoc);
                if( i > j )
                {
                    colSums[j] += Abs( A.GetLocal(iLoc,jLoc) ); 
                    rowSums[i-1] += Abs( A.GetLocal(iLoc,jLoc) );
                }
            }
        }
    }
    mpi::AllReduce( colSums.data(), n-1, g.VCComm() );
    mpi::AllReduce( rowSums.data(), n-1, g.VCComm() );

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
inline ValueInt<Base<F>>
SignDivide
( Matrix<F>& A, Matrix<F>& G, bool returnQ, const SdcCtrl<Base<F>>& sdcCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("schur::SignDivide"))

    // G := sgn(G)
    // G := 1/2 ( G + I )
    Sign( G, sdcCtrl.signCtrl );
    UpdateDiagonal( G, F(1) );
    Scale( F(1)/F(2), G );

    // Compute the pivoted QR decomposition of the spectral projection 
    Matrix<F> t;
    Matrix<Base<F>> d;
    Matrix<Int> p;
    elem::QR( G, t, d, p );

    // A := Q^H A Q
    const Base<F> oneA = OneNorm( A );
    if( returnQ )
    {
        ExpandPackedReflectors( LOWER, VERTICAL, CONJUGATED, 0, G, t );
        DiagonalScale( RIGHT, NORMAL, d, G );
        Matrix<F> B;
        Gemm( ADJOINT, NORMAL, F(1), G, A, B );
        Gemm( NORMAL, NORMAL, F(1), B, G, A );
    }
    else
    {
        qr::ApplyQ( LEFT, ADJOINT, G, t, d, A );
        qr::ApplyQ( RIGHT, NORMAL, G, t, d, A );
    }

    // Return || E21 ||1 / || A ||1 and the chosen rank
    auto part = ComputePartition( A );
    part.value /= oneA;
    return part;
}

template<typename F>
inline ValueInt<Base<F>>
SignDivide
( DistMatrix<F>& A, DistMatrix<F>& G, bool returnQ, 
  const SdcCtrl<Base<F>>& sdcCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("schur::SignDivide"))
    const Grid& g = A.Grid();

    // G := sgn(G)
    // G := 1/2 ( G + I )
    Sign( G, sdcCtrl.signCtrl );
    UpdateDiagonal( G, F(1) );
    Scale( F(1)/F(2), G );

    // Compute the pivoted QR decomposition of the spectral projection 
    DistMatrix<F,MD,STAR> t(g);
    DistMatrix<Base<F>,MD,STAR> d(g);
    DistMatrix<Int,VR,STAR> p(g);
    elem::QR( G, t, d, p );

    // A := Q^H A Q
    const Base<F> oneA = OneNorm( A );
    if( returnQ )
    {
        ExpandPackedReflectors( LOWER, VERTICAL, CONJUGATED, 0, G, t );
        DiagonalScale( RIGHT, NORMAL, d, G );
        DistMatrix<F> B(g);
        Gemm( ADJOINT, NORMAL, F(1), G, A, B );
        Gemm( NORMAL, NORMAL, F(1), B, G, A );
    }
    else
    {
        qr::ApplyQ( LEFT, ADJOINT, G, t, d, A );
        qr::ApplyQ( RIGHT, NORMAL, G, t, d, A );
    }

    // Return || E21 ||1 / || A ||1 and the chosen rank
    auto part = ComputePartition( A );
    part.value /= oneA;
    return part;
}

template<typename F>
inline ValueInt<Base<F>>
RandomizedSignDivide
( Matrix<F>& A, Matrix<F>& G, bool returnQ, const SdcCtrl<Base<F>>& sdcCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("schur::RandomizedSignDivide"))
    typedef Base<F> Real;
    const Int n = A.Height();
    const Real oneA = OneNorm( A );
    Real tol = sdcCtrl.tol;
    if( tol == Real(0) )
        tol = 500*n*lapack::MachineEpsilon<Real>();

    // S := sgn(G)
    // S := 1/2 ( S + I )
    auto S( G );
    Sign( S, sdcCtrl.signCtrl );
    UpdateDiagonal( S, F(1) );
    Scale( F(1)/F(2), S );

    ValueInt<Real> part;
    Matrix<F> V, B, t;
    Matrix<Base<F>> d;
    Int it=0;
    while( it < sdcCtrl.maxInnerIts )
    {
        G = S;

        // Compute the RURV of the spectral projector
        ImplicitHaar( V, t, d, n );
        qr::ApplyQ( RIGHT, NORMAL, V, t, d, G );
        elem::QR( G, t, d );

        // A := Q^H A Q [and reuse space for V for keeping original A]
        V = A;
        if( returnQ )
        {
            ExpandPackedReflectors( LOWER, VERTICAL, CONJUGATED, 0, G, t );
            DiagonalScale( RIGHT, NORMAL, d, G );
            Gemm( ADJOINT, NORMAL, F(1), G, A, B );
            Gemm( NORMAL, NORMAL, F(1), B, G, A );
        }
        else
        {
            qr::ApplyQ( LEFT, ADJOINT, G, t, d, A );
            qr::ApplyQ( RIGHT, NORMAL, G, t, d, A );
        }

        // || E21 ||1 / || A ||1 and the chosen rank
        part = ComputePartition( A );
        part.value /= oneA;

        ++it;
        if( part.value <= tol || it == sdcCtrl.maxInnerIts )
            break;
        else
            A = V;
    }
    return part;
}

template<typename F>
inline ValueInt<Base<F>>
RandomizedSignDivide
( DistMatrix<F>& A, DistMatrix<F>& G, bool returnQ,
  const SdcCtrl<Base<F>>& sdcCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("schur::RandomizedSignDivide"))
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    const Int n = A.Height();
    const Real oneA = OneNorm( A );
    Real tol = sdcCtrl.tol;
    if( tol == Real(0) )
        tol = 500*n*lapack::MachineEpsilon<Real>();

    // S := sgn(G)
    // S := 1/2 ( S + I )
    auto S( G );
    Sign( S, sdcCtrl.signCtrl );
    UpdateDiagonal( S, F(1) );
    Scale( F(1)/F(2), S );

    ValueInt<Real> part;
    DistMatrix<F> V(g), B(g);
    DistMatrix<F,MD,STAR> t(g);
    DistMatrix<Base<F>,MD,STAR> d(g);
    Int it=0;
    while( it < sdcCtrl.maxInnerIts )
    {
        G = S;

        // Compute the RURV of the spectral projector
        ImplicitHaar( V, t, d, n );
        qr::ApplyQ( RIGHT, NORMAL, V, t, d, G );
        elem::QR( G, t, d );

        // A := Q^H A Q [and reuse space for V for keeping original A]
        V = A;
        if( returnQ )
        {
            ExpandPackedReflectors( LOWER, VERTICAL, CONJUGATED, 0, G, t );
            DiagonalScale( RIGHT, NORMAL, d, G );
            Gemm( ADJOINT, NORMAL, F(1), G, A, B );
            Gemm( NORMAL, NORMAL, F(1), B, G, A );
        }
        else
        {
            qr::ApplyQ( LEFT, ADJOINT, G, t, d, A );
            qr::ApplyQ( RIGHT, NORMAL, G, t, d, A );
        }

        // || E21 ||1 / || A ||1 and the chosen rank
        part = ComputePartition( A );
        part.value /= oneA;

        ++it;
        if( part.value <= tol || it == sdcCtrl.maxInnerIts )
            break;
        else
            A = V;
    }
    return part;
}

template<typename Real>
inline ValueInt<Real>
SpectralDivide( Matrix<Real>& A, const SdcCtrl<Real>& sdcCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("schur::SpectralDivide"))
    const Int n = A.Height();
    const ValueInt<Real> median = Median(A.GetDiagonal());
    const Real infNorm = InfinityNorm(A);
    const Real eps = lapack::MachineEpsilon<Real>();
    Real tol = sdcCtrl.tol;
    if( tol == Real(0) )
        tol = 500*n*eps;
    const Real spread = sdcCtrl.spreadFactor*infNorm;

    Int it=0;
    ValueInt<Real> part;
    part.value = 2*tol; // initialize with unacceptable value
    Matrix<Real> G, ACopy;
    if( sdcCtrl.maxOuterIts > 1 )
        ACopy = A;
    while( it < sdcCtrl.maxOuterIts )
    {
        ++it;
        const Real shift = SampleBall<Real>(-median.value,spread);

        G = A;
        UpdateDiagonal( G, shift );

        if( sdcCtrl.progress )
            std::cout << "chose shift=" << shift << " using -median.value="
                      << -median.value << " and spread=" << spread << std::endl;

        try
        {
            if( sdcCtrl.random )
                part = RandomizedSignDivide( A, G, false, sdcCtrl );
            else
                part = SignDivide( A, G, false, sdcCtrl );

            if( part.value <= tol )
            {
                if( sdcCtrl.progress )
                    std::cout << "Converged during outer iter " << it-1
                              << std::endl;
                break;
            }
            else if( sdcCtrl.progress )
                std::cout << "part.value=" << part.value << " was greater than "
                          << tol << " during outer iter " << it-1 
                          << std::endl;
        } 
        catch( SingularMatrixException& e ) 
        { 
            if( sdcCtrl.progress )
                std::cout << "Caught singular matrix in outer iter " << it-1 
                          << std::endl;
        }
        if( it != sdcCtrl.maxOuterIts )
            A = ACopy;
    }
    if( part.value > tol )
        RuntimeError 
        ( "Unable to split spectrum to specified accuracy: part.value=",
          part.value, ", tol=", tol );

    return part;
}

template<typename Real>
inline ValueInt<Real>
SpectralDivide
( Matrix<Complex<Real>>& A, const SdcCtrl<Real>& sdcCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("schur::SpectralDivide"))
    typedef Complex<Real> F;
    const Int n = A.Height();
    const Real infNorm = InfinityNorm(A);
    const Real eps = lapack::MachineEpsilon<Real>();
    Real tol = sdcCtrl.tol;
    if( tol == Real(0) )
        tol = 500*n*eps;
    const Real spread = sdcCtrl.spreadFactor*infNorm;

    Int it=0;
    ValueInt<Real> part;
    part.value = 2*tol; // initialize with unacceptable value
    Matrix<F> G, ACopy;
    if( sdcCtrl.maxOuterIts > 1 )
        ACopy = A;
    while( it < sdcCtrl.maxOuterIts )
    {
        ++it;
        const Real angle = SampleUniform<Real>(0,2*Pi);
        const F gamma = F(Cos(angle),Sin(angle));
        G = A;
        Scale( gamma, G );

        const auto median = Median(G.GetRealPartOfDiagonal());
        const F shift = SampleBall<F>(-median.value,spread);
        UpdateDiagonal( G, shift );

        if( sdcCtrl.progress )
            std::cout << "chose gamma=" << gamma << " and shift=" << shift 
                      << " using -median.value=" << -median.value 
                      << " and spread=" << spread << std::endl;

        try
        {
            if( sdcCtrl.random )
                part = RandomizedSignDivide( A, G, false, sdcCtrl );
            else
                part = SignDivide( A, G, false, sdcCtrl );

            if( part.value <= tol )
            {
                if( sdcCtrl.progress )
                    std::cout << "Converged during outer iter " << it-1
                              << std::endl;
                break;
            }
            else if( sdcCtrl.progress )
                std::cout << "part.value=" << part.value << " was greater than "
                          << tol << " during outer iter " << it-1 
                          << std::endl;
        } 
        catch( SingularMatrixException& e ) 
        {
            if( sdcCtrl.progress )
                std::cout << "Caught singular matrix in outer iter " << it-1 
                          << std::endl;
        }
        if( it != sdcCtrl.maxOuterIts )
            A = ACopy;
    }
    if( part.value > tol )
        RuntimeError
        ( "Unable to split spectrum to specified accuracy: part.value=",
          part.value, ", tol=", tol );

    return part;
}

template<typename Real>
inline ValueInt<Real>
SpectralDivide
( Matrix<Real>& A, Matrix<Real>& Q, const SdcCtrl<Real>& sdcCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("schur::SpectralDivide"))
    const Int n = A.Height();
    const auto median = Median(A.GetDiagonal());
    const Real infNorm = InfinityNorm(A);
    const Real eps = lapack::MachineEpsilon<Real>();
    Real tol = sdcCtrl.tol;
    if( tol == Real(0) )
        tol = 500*n*eps;
    const Real spread = sdcCtrl.spreadFactor*infNorm;

    Int it=0;
    ValueInt<Real> part;
    part.value = 2*tol; // initialize with unacceptable value
    Matrix<Real> ACopy;
    if( sdcCtrl.maxOuterIts > 1 )
        ACopy = A;
    while( it < sdcCtrl.maxOuterIts )
    {
        ++it;
        const Real shift = SampleBall<Real>(-median.value,spread);

        Q = A;
        UpdateDiagonal( Q, shift );

        if( sdcCtrl.progress )
            std::cout << "chose shift=" << shift << " using -median.value="
                      << -median.value << " and spread=" << spread << std::endl;

        try
        {
            if( sdcCtrl.random )
                part = RandomizedSignDivide( A, Q, true, sdcCtrl );
            else
                part = SignDivide( A, Q, true, sdcCtrl );

            if( part.value <= tol )
            {
                if( sdcCtrl.progress )
                    std::cout << "Converged during outer iter " << it-1
                              << std::endl;
                break;
            }
            else if( sdcCtrl.progress )
                std::cout << "part.value=" << part.value << " was greater than "
                          << tol << " during outer iter " << it-1 
                          << std::endl;
        } 
        catch( SingularMatrixException& e ) 
        { 
            if( sdcCtrl.progress )
                std::cout << "Caught singular matrix in outer iter " << it-1 
                          << std::endl;
        }
        if( it != sdcCtrl.maxOuterIts )
            A = ACopy;
    }
    if( part.value > tol )
        RuntimeError
        ( "Unable to split spectrum to specified accuracy: part.value=",
          part.value, ", tol=", tol );

    return part;
}

template<typename Real>
inline ValueInt<Real>
SpectralDivide
( Matrix<Complex<Real>>& A, Matrix<Complex<Real>>& Q, 
  const SdcCtrl<Real>& sdcCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("schur::SpectralDivide"))
    typedef Complex<Real> F;
    const Int n = A.Height();
    const Real infNorm = InfinityNorm(A);
    const Real eps = lapack::MachineEpsilon<Real>();
    Real tol = sdcCtrl.tol;
    if( tol == Real(0) )
        tol = 500*n*eps;
    const Real spread = sdcCtrl.spreadFactor*infNorm;

    Int it=0;
    ValueInt<Real> part;
    part.value = 2*tol; // initialize with unacceptable value
    Matrix<F> ACopy;
    if( sdcCtrl.maxOuterIts > 1 )
        ACopy = A;
    while( it < sdcCtrl.maxOuterIts )
    {
        ++it;
        const Real angle = SampleUniform<Real>(0,2*Pi);
        const F gamma = F(Cos(angle),Sin(angle));
        Q = A;
        Scale( gamma, Q );

        const auto median = Median(Q.GetRealPartOfDiagonal());
        const F shift = SampleBall<F>(-median.value,spread);
        UpdateDiagonal( Q, shift );

        if( sdcCtrl.progress )
            std::cout << "chose gamma=" << gamma << " and shift=" << shift 
                      << " using -median.value=" << -median.value 
                      << " and spread=" << spread << std::endl;

        try
        {
            if( sdcCtrl.random )
                part = RandomizedSignDivide( A, Q, true, sdcCtrl );
            else
                part = SignDivide( A, Q, true, sdcCtrl );

            if( part.value <= tol )
            {
                if( sdcCtrl.progress )
                    std::cout << "Converged during outer iter " << it-1
                              << std::endl;
                break;
            }
            else if( sdcCtrl.progress )
                std::cout << "part.value=" << part.value << " was greater than "
                          << tol << " during outer iter " << it-1 
                          << std::endl;
        } 
        catch( SingularMatrixException& e ) 
        { 
            if( sdcCtrl.progress )
                std::cout << "Caught singular matrix in outer iter " << it-1 
                          << std::endl;
        }
        if( it != sdcCtrl.maxOuterIts )
            A = ACopy;
    }
    if( part.value > tol )
        RuntimeError
        ( "Unable to split spectrum to specified accuracy: part.value=",
          part.value, ", tol=", tol );

    return part;
}

template<typename Real>
inline ValueInt<Real>
SpectralDivide( DistMatrix<Real>& A, const SdcCtrl<Real>& sdcCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("schur::SpectralDivide"))
    const Int n = A.Height();
    const auto median = Median(A.GetDiagonal());
    const Real infNorm = InfinityNorm(A);
    const Real eps = lapack::MachineEpsilon<Real>();
    Real tol = sdcCtrl.tol;
    if( tol == Real(0) )
        tol = 500*n*eps;
    const Real spread = sdcCtrl.spreadFactor*infNorm;

    Int it=0;
    ValueInt<Real> part;
    part.value = 2*tol; // initialize with unacceptable value
    const Grid& g = A.Grid();
    DistMatrix<Real> ACopy(g), G(g);
    if( sdcCtrl.maxOuterIts > 1 )
        ACopy = A;
    while( it < sdcCtrl.maxOuterIts )
    {
        ++it;
        Real shift = SampleBall<Real>(-median.value,spread);
        mpi::Broadcast( shift, 0, g.VCComm() );

        G = A;
        UpdateDiagonal( G, shift );

        if( sdcCtrl.progress && g.Rank() == 0 )
            std::cout << "chose shift=" << shift << " using -median.value="
                      << -median.value << " and spread=" << spread << std::endl;

        try
        {
            if( sdcCtrl.random )
                part = RandomizedSignDivide( A, G, false, sdcCtrl );
            else
                part = SignDivide( A, G, false, sdcCtrl );

            if( part.value <= tol )
            {
                if( sdcCtrl.progress && g.Rank() == 0 )
                    std::cout << "Converged during outer iter " << it-1
                              << std::endl;
                break;
            }
            else if( sdcCtrl.progress && g.Rank() == 0 )
                std::cout << "part.value=" << part.value << " was greater than "
                          << tol << " during outer iter " << it-1 
                          << std::endl;
        } 
        catch( SingularMatrixException& e ) 
        { 
            if( sdcCtrl.progress && g.Rank() == 0 )
                std::cout << "Caught singular matrix in outer iter " << it-1 
                          << std::endl;
        }
        if( it != sdcCtrl.maxOuterIts )
            A = ACopy;
    }
    if( part.value > tol )
        RuntimeError
        ( "Unable to split spectrum to specified accuracy: part.value=",
          part.value, ", tol=", tol );

    return part;
}

template<typename Real>
inline ValueInt<Real>
SpectralDivide
( DistMatrix<Complex<Real>>& A, const SdcCtrl<Real>& sdcCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("schur::SpectralDivide"))
    typedef Complex<Real> F;
    const Int n = A.Height();
    const Real infNorm = InfinityNorm(A);
    const Real eps = lapack::MachineEpsilon<Real>();
    Real tol = sdcCtrl.tol;
    if( tol == Real(0) )
        tol = 500*n*eps;
    const Real spread = sdcCtrl.spreadFactor*infNorm;

    Int it=0;
    ValueInt<Real> part;
    part.value = 2*tol; // initialize with unacceptable value
    const Grid& g = A.Grid();
    DistMatrix<F> ACopy(g), G(g);
    if( sdcCtrl.maxOuterIts > 1 )
        ACopy = A;
    while( it < sdcCtrl.maxOuterIts )
    {
        ++it;
        const Real angle = SampleUniform<Real>(0,2*Pi);
        F gamma = F(Cos(angle),Sin(angle));
        mpi::Broadcast( gamma, 0, A.Grid().VCComm() );
        G = A;
        Scale( gamma, G );

        const auto median = Median(G.GetRealPartOfDiagonal());
        F shift = SampleBall<F>(-median.value,spread);
        mpi::Broadcast( shift, 0, g.VCComm() );
        UpdateDiagonal( G, shift );

        if( sdcCtrl.progress && g.Rank() == 0 )
            std::cout << "chose gamma=" << gamma << " and shift=" << shift 
                      << " using -median.value=" << -median.value 
                      << " and spread=" << spread << std::endl;

        try
        {
            if( sdcCtrl.random )
                part = RandomizedSignDivide( A, G, false, sdcCtrl );
            else
                part = SignDivide( A, G, false, sdcCtrl );

            if( part.value <= tol )
            {
                if( sdcCtrl.progress && g.Rank() == 0 )
                    std::cout << "Converged during outer iter " << it-1
                              << std::endl;
                break;
            }
            else if( sdcCtrl.progress && g.Rank() == 0 )
                std::cout << "part.value=" << part.value << " was greater than "
                          << tol << " during outer iter " << it-1 
                          << std::endl;
        } 
        catch( SingularMatrixException& e ) 
        { 
            if( sdcCtrl.progress && g.Rank() == 0 )
                std::cout << "Caught singular matrix in outer iter " << it-1 
                          << std::endl;
        }
        if( it != sdcCtrl.maxOuterIts )
            A = ACopy;
    }
    if( part.value > tol )
        RuntimeError
        ( "Unable to split spectrum to specified accuracy: part.value=",
          part.value, ", tol=", tol );

    return part;
}

template<typename Real>
inline ValueInt<Real>
SpectralDivide
( DistMatrix<Real>& A, DistMatrix<Real>& Q, const SdcCtrl<Real>& sdcCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("schur::SpectralDivide"))
    const Int n = A.Height();
    const Real infNorm = InfinityNorm(A);
    const auto median = Median(A.GetDiagonal());
    const Real eps = lapack::MachineEpsilon<Real>();
    Real tol = sdcCtrl.tol;
    if( tol == Real(0) )
        tol = 500*n*eps;
    const Real spread = sdcCtrl.spreadFactor*infNorm;

    Int it=0;
    ValueInt<Real> part;
    part.value = 2*tol; // initialize with unacceptable value
    const Grid& g = A.Grid();
    DistMatrix<Real> ACopy(g);
    if( sdcCtrl.maxOuterIts > 1 )
        ACopy = A;
    while( it < sdcCtrl.maxOuterIts )
    {
        ++it;
        Real shift = SampleBall<Real>(-median.value,spread);
        mpi::Broadcast( shift, 0, g.VCComm() );

        Q = A;
        UpdateDiagonal( Q, shift );

        if( sdcCtrl.progress && g.Rank() == 0 )
            std::cout << "chose shift=" << shift << " using -median.value=" 
                      << -median.value << " and spread=" << spread 
                      << std::endl;

        try
        {
            if( sdcCtrl.random )
                part = RandomizedSignDivide( A, Q, true, sdcCtrl );
            else
                part = SignDivide( A, Q, true, sdcCtrl );

            if( part.value <= tol )
            {
                if( sdcCtrl.progress && g.Rank() == 0 )
                    std::cout << "Converged during outer iter " << it-1
                              << std::endl;
                break;
            }
            else if( sdcCtrl.progress && g.Rank() == 0 )
                std::cout << "part.value=" << part.value << " was greater than "
                          << tol << " during outer iter " << it-1 
                          << std::endl;
        } 
        catch( SingularMatrixException& e ) 
        {
            if( sdcCtrl.progress && g.Rank() == 0 )
                std::cout << "Caught singular matrix in outer iter " << it-1 
                          << std::endl;
        }
        if( it != sdcCtrl.maxOuterIts )
            A = ACopy;
    }
    if( part.value > tol )
        RuntimeError
        ( "Unable to split spectrum to specified accuracy: part.value=",
          part.value, ", tol=", tol );

    return part;
}

template<typename Real>
inline ValueInt<Real>
SpectralDivide
( DistMatrix<Complex<Real>>& A, DistMatrix<Complex<Real>>& Q,
  const SdcCtrl<Real>& sdcCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("schur::SpectralDivide"))
    typedef Complex<Real> F;
    const Int n = A.Height();
    const Real infNorm = InfinityNorm(A);
    const Real eps = lapack::MachineEpsilon<Real>();
    Real tol = sdcCtrl.tol;
    if( tol == Real(0) )
        tol = 500*n*eps;
    const Real spread = sdcCtrl.spreadFactor*infNorm;

    Int it=0;
    ValueInt<Real> part;
    part.value = 2*tol; // initialize with unacceptable value
    const Grid& g = A.Grid();
    DistMatrix<F> ACopy(g);
    if( sdcCtrl.maxOuterIts > 1 )
        ACopy = A;
    while( it < sdcCtrl.maxOuterIts )
    {
        ++it;
        const Real angle = SampleUniform<Real>(0,2*Pi);
        F gamma = F(Cos(angle),Sin(angle));
        mpi::Broadcast( gamma, 0, g.VCComm() );
        Q = A;
        Scale( gamma, Q );

        const auto median = Median(Q.GetRealPartOfDiagonal());
        F shift = SampleBall<F>(-median.value,spread);
        mpi::Broadcast( shift, 0, g.VCComm() );
        UpdateDiagonal( Q, shift );

        if( sdcCtrl.progress && g.Rank() == 0 )
            std::cout << "chose gamma=" << gamma << " and shift=" << shift 
                      << " using -median.value=" << -median.value 
                      << " and spread=" << spread << std::endl;

        try
        {
            if( sdcCtrl.random )
                part = RandomizedSignDivide( A, Q, true, sdcCtrl );
            else
                part = SignDivide( A, Q, true, sdcCtrl );

            if( part.value <= tol )
            {
                if( sdcCtrl.progress && g.Rank() == 0 )
                    std::cout << "Converged during outer iter " << it-1
                              << std::endl;
                break;
            }
            else if( sdcCtrl.progress && g.Rank() == 0 )
                std::cout << "part.value=" << part.value << " was greater than "
                          << tol << " during outer iter " << it-1 
                          << std::endl;
        } 
        catch( SingularMatrixException& e ) 
        { 
            if( sdcCtrl.progress && g.Rank() == 0 )
                std::cout << "Caught singular matrix in outer iter " << it-1 
                          << std::endl;
        }
        if( it != sdcCtrl.maxOuterIts )
            A = ACopy;
    }
    if( part.value > tol )
        RuntimeError
        ( "Unable to split spectrum to specified accuracy: part.value=",
          part.value, ", tol=", tol );

    return part;
}

template<typename F>
inline void
SDC
( Matrix<F>& A, Matrix<Complex<Base<F>>>& w, 
  const SdcCtrl<Base<F>> sdcCtrl=SdcCtrl<Base<F>>() )
{
    DEBUG_ONLY(CallStackEntry cse("schur::SDC"))
    const Int n = A.Height();
    w.Resize( n, 1 );
    if( n <= sdcCtrl.cutoff )
    {
        if( sdcCtrl.progress )
            std::cout << n << " <= " << sdcCtrl.cutoff 
                      << ": switching to QR algorithm" << std::endl;
        schur::QR( A, w );
        return;
    }

    // Perform this level's split
    if( sdcCtrl.progress )
        std::cout << "Splitting " << n << " x " << n << " matrix" << std::endl;
    const auto part = SpectralDivide( A, sdcCtrl );
    Matrix<F> ATL, ATR,
              ABL, ABR;
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, part.index );
    MakeZeros( ABL );
    Matrix<Complex<Base<F>>> wT, wB;
    PartitionDown( w, wT, wB, part.index );

    // Recurse on the two subproblems
    if( sdcCtrl.progress )
        std::cout << "Recursing on " << ATL.Height() << " x " << ATL.Width() 
                  << " left subproblem" << std::endl;
    SDC( ATL, wT, sdcCtrl );
    if( sdcCtrl.progress )
        std::cout << "Recursing on " << ABR.Height() << " x " << ABR.Width() 
                  << " right subproblem" << std::endl;
    SDC( ABR, wB, sdcCtrl );
}

template<typename F>
inline void
SDC
( Matrix<F>& A, Matrix<Complex<Base<F>>>& w, Matrix<F>& Q, 
  bool fullTriangle=true, const SdcCtrl<Base<F>> sdcCtrl=SdcCtrl<Base<F>>() )
{
    DEBUG_ONLY(CallStackEntry cse("schur::SDC"))
    const Int n = A.Height();
    w.Resize( n, 1 );
    Q.Resize( n, n );
    if( n <= sdcCtrl.cutoff )
    {
        if( sdcCtrl.progress )
            std::cout << n << " <= " << sdcCtrl.cutoff 
                      << ": switching to QR algorithm" << std::endl;
        schur::QR( A, w, Q, fullTriangle );
        return;
    }

    // Perform this level's split
    if( sdcCtrl.progress )
        std::cout << "Splitting " << n << " x " << n << " matrix" << std::endl;
    const auto part = SpectralDivide( A, Q, sdcCtrl );
    Matrix<F> ATL, ATR,
              ABL, ABR;
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, part.index );
    MakeZeros( ABL );
    Matrix<Complex<Base<F>>> wT, wB;
    PartitionDown( w, wT, wB, part.index );
    Matrix<F> QL, QR;
    PartitionRight( Q, QL, QR, part.index );

    // Recurse on the top-left quadrant and update Schur vectors and ATR
    if( sdcCtrl.progress )
        std::cout << "Recursing on " << ATL.Height() << " x " << ATL.Width() 
                  << " left subproblem" << std::endl;
    Matrix<F> Z;
    SDC( ATL, wT, Z, fullTriangle, sdcCtrl );
    if( sdcCtrl.progress )
        std::cout << "Left subproblem update" << std::endl;
    auto G( QL );
    Gemm( NORMAL, NORMAL, F(1), G, Z, QL );
    if( fullTriangle )
        Gemm( ADJOINT, NORMAL, F(1), Z, ATR, G );

    // Recurse on the bottom-right quadrant and update Schur vectors and ATR
    if( sdcCtrl.progress )
        std::cout << "Recursing on " << ABR.Height() << " x " << ABR.Width() 
                  << " right subproblem" << std::endl;
    SDC( ABR, wB, Z, fullTriangle, sdcCtrl );
    if( sdcCtrl.progress )
        std::cout << "Right subproblem update" << std::endl;
    if( fullTriangle )
        Gemm( NORMAL, NORMAL, F(1), G, Z, ATR ); 
    G = QR;
    Gemm( NORMAL, NORMAL, F(1), G, Z, QR );
}

// This routine no longer attempts to evenly assign work/process between two
// teams since it was found to lead to horrendously non-square process grids
// in practice, even when the original number of processes was a large power
// of two. Instead, the grid is either split in half, or not split at all.
// The choice is made based upon whether or not one subproblem requires twice
// as much work as the other. There is a complicated calculus here that would
// require a much more sophisticated (machine- and problem-specific) model to
// make the 'best' splitting, but this approach should be a good compromise.
inline void SplitGrid
( int nLeft, int nRight, const Grid& grid, 
  const Grid*& leftGrid, const Grid*& rightGrid, bool progress=false )
{
    typedef double Real;
    const Real leftWork = Pow(Real(nLeft),Real(3));
    const Real rightWork = Pow(Real(nRight),Real(3));
    if( Max(leftWork,rightWork) > 2*Min(leftWork,rightWork) )
    {
        // Don't split the grid
        leftGrid = &grid;
        rightGrid = &grid;
        if( progress && grid.Rank() == 0 )
            std::cout << "leftWork/rightWork=" << leftWork/rightWork 
                      << ", so the grid was not split" << std::endl;
    }
    else
    {
        // Split the grid in half (powers-of-two remain so)
        const Int p = grid.Size();
        const Int pLeft = p/2;
        const Int pRight = p-pLeft;
        std::vector<int> leftRanks(pLeft), rightRanks(pRight);
        for( int j=0; j<pLeft; ++j )
            leftRanks[j] = j;
        for( int j=0; j<pRight; ++j )
            rightRanks[j] = j+pLeft;
        mpi::Group group = grid.OwningGroup();
        mpi::Group leftGroup, rightGroup;
        mpi::Incl( group, pLeft, leftRanks.data(), leftGroup );
        mpi::Incl( group, pRight, rightRanks.data(), rightGroup );
        const Int rLeft = Grid::FindFactor(pLeft);
        const Int rRight = Grid::FindFactor(pRight);
        if( progress && grid.Rank() == 0 )
            std::cout << "leftWork/rightWork=" << leftWork/rightWork 
                      << ", so split " << p << " processes into " 
                      << rLeft << " x " << pLeft/rLeft << " and "
                      << rRight << " x " << pRight/rRight << " grids" 
                      << std::endl;
        leftGrid = new Grid( grid.VCComm(), leftGroup, rLeft );
        rightGrid = new Grid( grid.VCComm(), rightGroup, rRight );
        mpi::Free( leftGroup );
        mpi::Free( rightGroup );
    }
}

template<typename F,typename EigType>
inline void PushSubproblems
( DistMatrix<F>& ATL,    DistMatrix<F>& ABR, 
  DistMatrix<F>& ATLSub, DistMatrix<F>& ABRSub,
  DistMatrix<EigType,VR,STAR>& wT,    
  DistMatrix<EigType,VR,STAR>& wB,
  DistMatrix<EigType,VR,STAR>& wTSub, 
  DistMatrix<EigType,VR,STAR>& wBSub,
  bool progress=false )
{
    DEBUG_ONLY(CallStackEntry cse("schur::PushSubproblems"))
    const Grid& grid = ATL.Grid();

    // Split based on the work estimates
    const Grid *leftGrid, *rightGrid;
    SplitGrid
    ( ATL.Height(), ABR.Height(), grid, leftGrid, rightGrid, progress );
    ATLSub.SetGrid( *leftGrid ); 
    ABRSub.SetGrid( *rightGrid );
    wTSub.SetGrid( *leftGrid );
    wBSub.SetGrid( *rightGrid );
    if( progress && grid.Rank() == 0 )
        std::cout << "Pushing ATL and ABR" << std::endl;
    ATLSub = ATL;
    ABRSub = ABR;
}

template<typename F,typename EigType>
inline void PullSubproblems
( DistMatrix<F>& ATL,    DistMatrix<F>& ABR,
  DistMatrix<F>& ATLSub, DistMatrix<F>& ABRSub,
  DistMatrix<EigType,VR,STAR>& wT,    
  DistMatrix<EigType,VR,STAR>& wB,
  DistMatrix<EigType,VR,STAR>& wTSub, 
  DistMatrix<EigType,VR,STAR>& wBSub,
  bool progress=false )
{
    DEBUG_ONLY(CallStackEntry cse("schur::PullSubproblems"))
    const Grid& grid = ATL.Grid();
    const bool sameGrid = ( wT.Grid() == wTSub.Grid() );

    if( progress && grid.Rank() == 0 )
        std::cout << "Pulling ATL and ABR" << std::endl;
    ATL = ATLSub;
    ABR = ABRSub;

    // This section is a hack since no inter-grid redistributions exist for 
    // [VR,* ] distributions yet
    if( progress && grid.Rank() == 0 )
        std::cout << "Pulling wT and wB" << std::endl;
    if( sameGrid )
    {
        wT = wTSub;
        wB = wBSub;
    }
    else
    {
        DistMatrix<EigType> wTSub_MC_MR( wTSub.Grid() );
        if( wTSub.Participating() )
            wTSub_MC_MR = wTSub;
        wTSub_MC_MR.MakeConsistent();
        DistMatrix<EigType> wT_MC_MR(wT.Grid()); 
        wT_MC_MR = wTSub_MC_MR;
        wT = wT_MC_MR;

        DistMatrix<EigType> wBSub_MC_MR( wBSub.Grid() );
        if( wBSub.Participating() )
            wBSub_MC_MR = wBSub;
        wBSub_MC_MR.MakeConsistent();
        DistMatrix<EigType> wB_MC_MR(wB.Grid()); 
        wB_MC_MR = wBSub_MC_MR;
        wB = wB_MC_MR;
    }
    
    const Grid *leftGrid = &ATLSub.Grid();
    const Grid *rightGrid = &ABRSub.Grid();
    ATLSub.Empty();
    ABRSub.Empty();
    wTSub.Empty();
    wBSub.Empty();
    if( !sameGrid )
    {
        delete leftGrid;
        delete rightGrid;
    }
}

template<typename F>
inline void
SDC
( DistMatrix<F>& A, DistMatrix<Complex<Base<F>>>& w, 
  const SdcCtrl<Base<F>> sdcCtrl=SdcCtrl<Base<F>>() )
{
    DEBUG_ONLY(CallStackEntry cse("schur::SDC"))
    const Grid& g = A.Grid();
    const Int n = A.Height();
    w.Resize( n, 1 );
    if( A.Grid().Size() == 1 )
    {
        if( sdcCtrl.progress && g.Rank() == 0 )
            std::cout << "One process: using QR algorithm" << std::endl;
        schur::QR( A.Matrix(), w.Matrix() );
        return;
    }
    if( n <= sdcCtrl.cutoff )
    {
        if( sdcCtrl.progress && g.Rank() == 0 )
            std::cout << n << " <= " << sdcCtrl.cutoff 
                      << ": using QR algorithm" << std::endl;
        DistMatrix<F,CIRC,CIRC> A_CIRC_CIRC( A );
        DistMatrix<Complex<Base<F>>,CIRC,CIRC> w_CIRC_CIRC( w );
        if( A_CIRC_CIRC.CrossRank() == A_CIRC_CIRC.Root() )
            schur::QR( A_CIRC_CIRC.Matrix(), w_CIRC_CIRC.Matrix() );
        A = A_CIRC_CIRC;
        w = w_CIRC_CIRC;
        return;
    }

    // Perform this level's split
    if( sdcCtrl.progress && g.Rank() == 0 )
        std::cout << "Splitting " << n << " x " << n << " matrix" << std::endl;
    const auto part = SpectralDivide( A, sdcCtrl );
    DistMatrix<F> ATL(g), ATR(g),
                  ABL(g), ABR(g);
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, part.index );
    MakeZeros( ABL );
    DistMatrix<Complex<Base<F>>,VR,STAR> wT(g), wB(g);
    PartitionDown( w, wT, wB, part.index );

    if( sdcCtrl.progress && g.Rank() == 0 )
        std::cout << "Pushing subproblems" << std::endl;
    DistMatrix<F> ATLSub, ABRSub;
    DistMatrix<Complex<Base<F>>,VR,STAR> wTSub, wBSub;
    PushSubproblems
    ( ATL, ABR, ATLSub, ABRSub, wT, wB, wTSub, wBSub, sdcCtrl.progress );
    if( ATLSub.Participating() )
        SDC( ATLSub, wTSub, sdcCtrl );
    if( ABRSub.Participating() )
        SDC( ABRSub, wBSub, sdcCtrl );
    if( sdcCtrl.progress && g.Rank() == 0 )
        std::cout << "Pulling subproblems" << std::endl;
    PullSubproblems
    ( ATL, ABR, ATLSub, ABRSub, wT, wB, wTSub, wBSub, sdcCtrl.progress );
}

template<typename F,typename EigType>
inline void PushSubproblems
( DistMatrix<F>& ATL,    DistMatrix<F>& ABR, 
  DistMatrix<F>& ATLSub, DistMatrix<F>& ABRSub,
  DistMatrix<EigType,VR,STAR>& wT,    
  DistMatrix<EigType,VR,STAR>& wB,
  DistMatrix<EigType,VR,STAR>& wTSub, 
  DistMatrix<EigType,VR,STAR>& wBSub,
  DistMatrix<F>& ZTSub,  DistMatrix<F>& ZBSub,
  bool progress=false )
{
    DEBUG_ONLY(CallStackEntry cse("schur::PushSubproblems"))
    const Grid& grid = ATL.Grid();

    // Split based on the work estimates
    const Grid *leftGrid, *rightGrid;
    SplitGrid
    ( ATL.Height(), ABR.Height(), grid, leftGrid, rightGrid, progress );
    ATLSub.SetGrid( *leftGrid );
    ABRSub.SetGrid( *rightGrid );
    wTSub.SetGrid( *leftGrid );
    wBSub.SetGrid( *rightGrid );
    ZTSub.SetGrid( *leftGrid );
    ZBSub.SetGrid( *rightGrid );
    if( progress && grid.Rank() == 0 )
        std::cout << "Pushing ATLSub" << std::endl;
    ATLSub = ATL;
    if( progress && grid.Rank() == 0 )
        std::cout << "Pushing ABRSub" << std::endl;
    ABRSub = ABR;
}

template<typename F,typename EigType>
inline void PullSubproblems
( DistMatrix<F>& ATL,    DistMatrix<F>& ABR,
  DistMatrix<F>& ATLSub, DistMatrix<F>& ABRSub,
  DistMatrix<EigType,VR,STAR>& wT,    
  DistMatrix<EigType,VR,STAR>& wB,
  DistMatrix<EigType,VR,STAR>& wTSub, 
  DistMatrix<EigType,VR,STAR>& wBSub,
  DistMatrix<F>& ZT,     DistMatrix<F>& ZB,
  DistMatrix<F>& ZTSub,  DistMatrix<F>& ZBSub,
  bool progress=false )
{
    DEBUG_ONLY(CallStackEntry cse("schur::PullSubproblems"))
    const Grid& grid = ATL.Grid();
    const bool sameGrid = ( wT.Grid() == wTSub.Grid() );

    if( progress && grid.Rank() == 0 )
        std::cout << "Pulling ATL and ABR" << std::endl;
    ATL = ATLSub;
    ABR = ABRSub;

    // This section is a hack for since no intergrid redistributions exist for 
    // [VR,* ] distributions yet
    if( progress && grid.Rank() == 0 )
        std::cout << "Pulling wT and wB" << std::endl;
    if( sameGrid )
    {
        wT = wTSub;
        wB = wBSub;
    }
    else
    {
        DistMatrix<EigType> wTSub_MC_MR( wTSub.Grid() );
        if( wTSub.Participating() )
            wTSub_MC_MR = wTSub;
        wTSub_MC_MR.MakeConsistent();
        DistMatrix<EigType> wT_MC_MR(wT.Grid()); 
        wT_MC_MR = wTSub_MC_MR;
        wT = wT_MC_MR;

        DistMatrix<EigType> wBSub_MC_MR( wBSub.Grid() );
        if( wBSub.Participating() )
            wBSub_MC_MR = wBSub;
        wBSub_MC_MR.MakeConsistent();
        DistMatrix<EigType> wB_MC_MR(wB.Grid()); 
        wB_MC_MR = wBSub_MC_MR;
        wB = wB_MC_MR;
    }

    if( progress && grid.Rank() == 0 )
        std::cout << "Pulling ZT and ZB" << std::endl;
    if( !sameGrid )
    {
        ZTSub.MakeConsistent();
        ZBSub.MakeConsistent();
    }
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
    if( !sameGrid )
    {
        delete leftGrid;
        delete rightGrid;
    }
}

template<typename F>
inline void
SDC
( DistMatrix<F>& A, DistMatrix<Complex<Base<F>>,VR,STAR>& w, DistMatrix<F>& Q, 
  bool fullTriangle=true, const SdcCtrl<Base<F>> sdcCtrl=SdcCtrl<Base<F>>() )
{
    DEBUG_ONLY(CallStackEntry cse("schur::SDC"))
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    const Int n = A.Height();
    w.Resize( n, 1 );
    Q.Resize( n, n );
    if( A.Grid().Size() == 1 )
    {
        if( sdcCtrl.progress && g.Rank() == 0 )
            std::cout << "One process: using QR algorithm" << std::endl;
        schur::QR( A.Matrix(), w.Matrix(), Q.Matrix(), fullTriangle );
        return;
    }
    if( n <= sdcCtrl.cutoff )
    {
        if( sdcCtrl.progress && g.Rank() == 0 )
            std::cout << n << " <= " << sdcCtrl.cutoff 
                      << ": using QR algorithm" << std::endl;
        DistMatrix<F,CIRC,CIRC> A_CIRC_CIRC( A ), Q_CIRC_CIRC( n, n, g );
        DistMatrix<Complex<Base<F>>,CIRC,CIRC> w_CIRC_CIRC( n, 1, g );
        if( A_CIRC_CIRC.CrossRank() == A_CIRC_CIRC.Root() )
            schur::QR
            ( A_CIRC_CIRC.Matrix(), w_CIRC_CIRC.Matrix(), Q_CIRC_CIRC.Matrix(),
              fullTriangle );
        A = A_CIRC_CIRC;
        w = w_CIRC_CIRC;
        Q = Q_CIRC_CIRC;
        return;
    }

    // Perform this level's split
    if( sdcCtrl.progress && g.Rank() == 0 )
        std::cout << "Splitting " << n << " x " << n << " matrix" << std::endl;
    const Real infNorm = InfinityNorm( A );
    const auto part = SpectralDivide( A, Q, sdcCtrl );
    DistMatrix<F> ATL(g), ATR(g),
                  ABL(g), ABR(g);
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, part.index );
    MakeZeros( ABL );
    DistMatrix<Complex<Base<F>>,VR,STAR> wT(g), wB(g);
    PartitionDown( w, wT, wB, part.index );
    DistMatrix<F> QL(g), QR(g);
    PartitionRight( Q, QL, QR, part.index );

    // Recurse on the two subproblems
    DistMatrix<F> ATLSub, ABRSub, ZTSub, ZBSub;
    DistMatrix<Complex<Base<F>>,VR,STAR> wTSub, wBSub;
    if( sdcCtrl.progress && g.Rank() == 0 )
        std::cout << "Pushing subproblems" << std::endl;
    PushSubproblems
    ( ATL, ABR, ATLSub, ABRSub, wT, wB, wTSub, wBSub, ZTSub, ZBSub, 
      sdcCtrl.progress );
    if( ATLSub.Participating() )
        SDC( ATLSub, wTSub, ZTSub, fullTriangle, sdcCtrl );
    if( ABRSub.Participating() )
        SDC( ABRSub, wBSub, ZBSub, fullTriangle, sdcCtrl );
    
    // Ensure that the results are back on this level's grid
    if( sdcCtrl.progress && g.Rank() == 0 )
        std::cout << "Pulling subproblems" << std::endl;
    DistMatrix<F> ZT(g), ZB(g);
    PullSubproblems
    ( ATL, ABR, ATLSub, ABRSub, wT, wB, wTSub, wBSub, ZT, ZB, ZTSub, ZBSub,
      sdcCtrl.progress );

    // Update the Schur vectors
    if( sdcCtrl.progress && g.Rank() == 0 )
        std::cout << "Updating Schur vectors" << std::endl;
    auto G( QL );
    Gemm( NORMAL, NORMAL, F(1), G, ZT, QL );
    G = QR;
    Gemm( NORMAL, NORMAL, F(1), G, ZB, QR );

    if( fullTriangle )
    {
        if( sdcCtrl.progress && g.Rank() == 0 )
            std::cout << "Updating top-right quadrant" << std::endl;
        // Update the top-right quadrant
        Gemm( ADJOINT, NORMAL, F(1), ZT, ATR, G );
        Gemm( NORMAL, NORMAL, F(1), G, ZB, ATR ); 
    }
}

} // namespace schur
} // namespace elem

#endif // ifndef ELEM_SCHUR_SDC_HPP
