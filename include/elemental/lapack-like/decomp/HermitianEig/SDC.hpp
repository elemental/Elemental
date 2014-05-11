/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_HERMITIANEIG_SDC_HPP
#define ELEM_HERMITIANEIG_SDC_HPP

#include ELEM_AXPY_INC
#include ELEM_SCALE_INC
#include ELEM_UPDATEDIAGONAL_INC

#include ELEM_QR_INC
#include ELEM_POLAR_QDWH_INC
#include ELEM_SCHUR_INC

#include ELEM_MEDIAN_INC
#include ELEM_INFINITYNORM_INC
#include ELEM_ONENORM_INC

#include ELEM_HAAR_INC

// TODO: Reference to Yuji's work

namespace elem {

namespace herm_eig {

using elem::schur::ComputePartition;
using elem::schur::SplitGrid;
using elem::schur::PushSubproblems;
using elem::schur::PullSubproblems;

// TODO: Exploit symmetry in A := Q^H A Q. Routine for A := X^H A X?

// G should be a rational function of A. If returnQ=true, G will be set to
// the computed unitary matrix upon exit.
template<typename F>
inline ValueInt<Base<F>>
QDWHDivide
( UpperOrLower uplo, Matrix<F>& A, Matrix<F>& G, bool returnQ=false )
{
    DEBUG_ONLY(CallStackEntry cse("herm_eig::QDWHDivide"))

    // G := sgn(G)
    // G := 1/2 ( G + I )
    herm_polar::QDWH( uplo, G ); 
    UpdateDiagonal( G, F(1) );
    Scale( F(1)/F(2), G );

    // Compute the pivoted QR decomposition of the spectral projection 
    Matrix<F> t;
    Matrix<Base<F>> d;
    Matrix<Int> p;
    elem::QR( G, t, d, p );

    // A := Q^H A Q
    MakeHermitian( uplo, A );
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
QDWHDivide
( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<F>& G, bool returnQ=false )
{
    DEBUG_ONLY(CallStackEntry cse("herm_eig::QDWHDivide"))

    // G := sgn(G)
    // G := 1/2 ( G + I )
    herm_polar::QDWH( uplo, G ); 
    UpdateDiagonal( G, F(1) );
    Scale( F(1)/F(2), G );

    // Compute the pivoted QR decomposition of the spectral projection 
    const Grid& g = A.Grid();
    DistMatrix<F,MD,STAR> t(g);
    DistMatrix<Base<F>,MD,STAR> d(g);
    DistMatrix<Int,VR,STAR> p(g);
    elem::QR( G, t, d, p );

    // A := Q^H A Q
    MakeHermitian( uplo, A );
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
( UpperOrLower uplo, Matrix<F>& A, Matrix<F>& G, bool returnQ, 
  const HermitianSdcCtrl<Base<F>>& sdcCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("herm_eig::RandomizedSignDivide"))

    typedef Base<F> Real;
    const Int n = A.Height();
    MakeHermitian( uplo, A );
    const Real oneA = OneNorm( A );
    Real tol = sdcCtrl.tol;
    if( tol == Real(0) )
        tol = 500*n*lapack::MachineEpsilon<Real>();

    // S := sgn(G)
    // S := 1/2 ( S + I )
    auto S( G );
    herm_polar::QDWH( uplo, S ); 
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
( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<F>& G, bool returnQ, 
  const HermitianSdcCtrl<Base<F>>& sdcCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("herm_eig::RandomizedSignDivide"))

    typedef Base<F> Real;
    const Grid& g = A.Grid();
    const Int n = A.Height();
    MakeHermitian( uplo, A );
    const Real oneA = OneNorm( A );
    Real tol = sdcCtrl.tol;
    if( tol == Real(0) )
        tol = 500*n*lapack::MachineEpsilon<Real>();

    // S := sgn(G)
    // S := 1/2 ( S + I )
    auto S( G );
    herm_polar::QDWH( uplo, S ); 
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

template<typename F>
inline ValueInt<Base<F>>
SpectralDivide
( UpperOrLower uplo, Matrix<F>& A, const HermitianSdcCtrl<Base<F>>& sdcCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("herm_eig::SpectralDivide"))

    typedef Base<F> Real;
    const Int n = A.Height();
    MakeHermitian( uplo, A );
    const auto median = Median(A.GetRealPartOfDiagonal());
    const Real infNorm = InfinityNorm(A);
    const Real eps = lapack::MachineEpsilon<Real>();
    Real tol = sdcCtrl.tol;
    if( tol == Real(0) )
        tol = 500*n*eps;
    const Real spread = sdcCtrl.spreadFactor*infNorm;

    Int it=0;
    ValueInt<Real> part;
    Matrix<F> G, ACopy;
    if( sdcCtrl.maxOuterIts > 1 )
        ACopy = A;
    while( it < sdcCtrl.maxOuterIts )
    {
        const Real shift = SampleBall<Real>(-median.value,spread);

        G = A;
        UpdateDiagonal( G, F(shift) );

        part = RandomizedSignDivide( uplo, A, G, false, sdcCtrl );

        ++it;
        if( part.value <= tol )
            break;
        else if( it != sdcCtrl.maxOuterIts )
            A = ACopy;
    }
    if( part.value > tol )
        RuntimeError
        ( "Unable to split spectrum to specified accuracy: part.value=",
          part.value, ", tol=", tol );

    return part;
}

template<typename F>
inline ValueInt<Base<F>>
SpectralDivide
( UpperOrLower uplo, Matrix<F>& A, Matrix<F>& Q, 
  const HermitianSdcCtrl<Base<F>>& sdcCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("herm_eig::SpectralDivide"))

    typedef Base<F> Real;
    const Int n = A.Height();
    MakeHermitian( uplo, A );
    const auto median = Median(A.GetRealPartOfDiagonal());
    const Real infNorm = InfinityNorm(A);
    const Real eps = lapack::MachineEpsilon<Real>();
    Real tol = sdcCtrl.tol;
    if( tol == Real(0) )
        tol = 500*n*eps;
    const Real spread = sdcCtrl.spreadFactor*infNorm;

    Int it=0;
    ValueInt<Real> part;
    Matrix<F> ACopy;
    if( sdcCtrl.maxOuterIts > 1 )
        ACopy = A;
    while( it < sdcCtrl.maxOuterIts )
    {
        const Real shift = SampleBall<Real>(-median.value,spread);

        Q = A;
        UpdateDiagonal( Q, F(shift) );

        part = RandomizedSignDivide( uplo, A, Q, true, sdcCtrl );

        ++it;
        if( part.value <= tol )
            break;
        else if( it != sdcCtrl.maxOuterIts )
            A = ACopy;
    }
    if( part.value > tol )
        RuntimeError
        ( "Unable to split spectrum to specified accuracy: part.value=",
          part.value, ", tol=", tol );

    return part;
}

template<typename F>
inline ValueInt<Base<F>>
SpectralDivide
( UpperOrLower uplo, DistMatrix<F>& A, 
  const HermitianSdcCtrl<Base<F>>& sdcCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("herm_eig::SpectralDivide"))

    typedef Base<F> Real;
    const Int n = A.Height();
    MakeHermitian( uplo, A );
    const auto median = Median(A.GetRealPartOfDiagonal());
    const Real infNorm = InfinityNorm(A);
    const Real eps = lapack::MachineEpsilon<Real>();
    Real tol = sdcCtrl.tol;
    if( tol == Real(0) )
        tol = 500*n*eps;
    const Real spread = sdcCtrl.spreadFactor*infNorm;

    Int it=0;
    ValueInt<Real> part;
    DistMatrix<F> ACopy(A.Grid()), G(A.Grid());
    if( sdcCtrl.maxOuterIts > 1 )
        ACopy = A;
    while( it < sdcCtrl.maxOuterIts )
    {
        Real shift = SampleBall<Real>(-median.value,spread);
        mpi::Broadcast( shift, 0, A.Grid().VCComm() );

        G = A;
        UpdateDiagonal( G, F(shift) );

        part = RandomizedSignDivide( uplo, A, G, false, sdcCtrl );

        ++it;
        if( part.value <= tol )
            break;
        else if( it != sdcCtrl.maxOuterIts )
            A = ACopy;
    }
    if( part.value > tol )
        RuntimeError
        ( "Unable to split spectrum to specified accuracy: part.value=",
          part.value, ", tol=", tol );

    return part;
}

template<typename F>
inline ValueInt<Base<F>>
SpectralDivide
( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<F>& Q, 
  const HermitianSdcCtrl<Base<F>>& sdcCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("herm_eig::SpectralDivide"))

    typedef Base<F> Real;
    const Int n = A.Height();
    MakeHermitian( uplo, A );
    const Real infNorm = InfinityNorm(A);
    const auto median = Median(A.GetRealPartOfDiagonal());
    const Real eps = lapack::MachineEpsilon<Real>();
    Real tol = sdcCtrl.tol;
    if( tol == Real(0) )
        tol = 500*n*eps;
    const Real spread = sdcCtrl.spreadFactor*infNorm;

    Int it=0;
    ValueInt<Real> part;
    DistMatrix<F> ACopy(A.Grid());
    if( sdcCtrl.maxOuterIts > 1 )
        ACopy = A;
    while( it < sdcCtrl.maxOuterIts )
    {
        Real shift = SampleBall<Real>(-median.value,spread);
        mpi::Broadcast( shift, 0, A.Grid().VCComm() );

        Q = A;
        UpdateDiagonal( Q, F(shift) );

        part = RandomizedSignDivide( uplo, A, Q, true, sdcCtrl );

        ++it;
        if( part.value <= tol )
            break;
        else if( it != sdcCtrl.maxOuterIts )
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
( UpperOrLower uplo, Matrix<F>& A, Matrix<Base<F>>& w, 
  const HermitianSdcCtrl<Base<F>> sdcCtrl=HermitianSdcCtrl<Base<F>>() )
{
    DEBUG_ONLY(CallStackEntry cse("herm_eig::SDC"))

    typedef Base<F> Real;
    const Int n = A.Height();
    w.Resize( n, 1 );
    if( n <= sdcCtrl.cutoff )
    {
        HermitianEig( uplo, A, w );
        return;
    }

    // Perform this level's split
    const auto part = SpectralDivide( uplo, A, sdcCtrl );
    Matrix<F> ATL, ATR,
              ABL, ABR;
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, part.index );
    if( uplo == LOWER )
        MakeZeros( ABL );
    else
        MakeZeros( ATR );
    Matrix<Real> wT, wB;
    PartitionDown( w, wT, wB, part.index );

    // Recurse on the two subproblems
    SDC( uplo, ATL, wT, sdcCtrl );
    SDC( uplo, ABR, wB, sdcCtrl );
}

template<typename F>
inline void
SDC
( UpperOrLower uplo, Matrix<F>& A, Matrix<Base<F>>& w, Matrix<F>& Q, 
  const HermitianSdcCtrl<Base<F>> sdcCtrl=HermitianSdcCtrl<Base<F>>() )
{
    DEBUG_ONLY(CallStackEntry cse("herm_eig::SDC"))

    typedef Base<F> Real;
    const Int n = A.Height();
    w.Resize( n, 1 );
    Q.Resize( n, n );
    if( n <= sdcCtrl.cutoff )
    {
        HermitianEig( uplo, A, w, Q );
        return;
    }

    // Perform this level's split
    const auto part = SpectralDivide( uplo, A, Q, sdcCtrl );
    Matrix<F> ATL, ATR,
              ABL, ABR;
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, part.index );
    if( uplo == LOWER )
        MakeZeros( ABL );
    else
        MakeZeros( ATR );
    Matrix<Real> wT, wB;
    PartitionDown( w, wT, wB, part.index );
    Matrix<F> QL, QR;
    PartitionRight( Q, QL, QR, part.index );

    // Recurse on the top-left quadrant and update eigenvectors
    Matrix<F> Z;
    SDC( uplo, ATL, wT, Z, sdcCtrl );
    auto G( QL );
    Gemm( NORMAL, NORMAL, F(1), G, Z, QL );

    // Recurse on the bottom-right quadrant and update eigenvectors
    SDC( uplo, ABR, wB, Z, sdcCtrl );
    G = QR;
    Gemm( NORMAL, NORMAL, F(1), G, Z, QR );
}

template<typename F>
inline void
SDC
( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<Base<F>,VR,STAR>& w, 
  const HermitianSdcCtrl<Base<F>> sdcCtrl=HermitianSdcCtrl<Base<F>>() )
{
    DEBUG_ONLY(CallStackEntry cse("herm_eig::SDC"))

    typedef Base<F> Real;
    const Grid& g = A.Grid();
    const Int n = A.Height();
    w.Resize( n, 1 );
    if( A.Grid().Size() == 1 )
    {
        HermitianEig( uplo, A.Matrix(), w.Matrix() );
        return;
    }
    if( n <= sdcCtrl.cutoff )
    {
        HermitianEig( uplo, A, w );
        return;
    }

    // Perform this level's split
    const auto part = SpectralDivide( uplo, A, sdcCtrl );
    DistMatrix<F> ATL(g), ATR(g),
                  ABL(g), ABR(g);
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, part.index );
    if( uplo == LOWER )
        MakeZeros( ABL );
    else
        MakeZeros( ATR );
    DistMatrix<Real,VR,STAR> wT(g), wB(g);
    PartitionDown( w, wT, wB, part.index );

    // Recurse on the two subproblems
    DistMatrix<F> ATLSub, ABRSub;
    DistMatrix<Real,VR,STAR> wTSub, wBSub;
    PushSubproblems
    ( ATL, ABR, ATLSub, ABRSub, wT, wB, wTSub, wBSub, sdcCtrl.progress );
    if( ATL.Participating() )
        SDC( uplo, ATLSub, wTSub, sdcCtrl );
    if( ABR.Participating() )
        SDC( uplo, ABRSub, wBSub, sdcCtrl );
    PullSubproblems( ATL, ABR, ATLSub, ABRSub, wT, wB, wTSub, wBSub );
}

template<typename F>
inline void
SDC
( UpperOrLower uplo, 
  DistMatrix<F>& A, DistMatrix<Base<F>,VR,STAR>& w, DistMatrix<F>& Q, 
  const HermitianSdcCtrl<Base<F>> sdcCtrl=HermitianSdcCtrl<Base<F>>() )
{
    DEBUG_ONLY(CallStackEntry cse("herm_eig::SDC"))

    typedef Base<F> Real;
    const Grid& g = A.Grid();
    const Int n = A.Height();
    w.Resize( n, 1 );
    Q.Resize( n, n );
    if( A.Grid().Size() == 1 )
    {
        HermitianEig( uplo, A.Matrix(), w.Matrix(), Q.Matrix() );
        return;
    }
    if( n <= sdcCtrl.cutoff )
    {
        HermitianEig( uplo, A, w, Q );
        return;
    }

    // Perform this level's split
    const auto part = SpectralDivide( uplo, A, Q, sdcCtrl );
    DistMatrix<F> ATL(g), ATR(g),
                  ABL(g), ABR(g);
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, part.index );
    if( uplo == LOWER )
        MakeZeros( ABL );
    else
        MakeZeros( ATR );
    DistMatrix<Real,VR,STAR> wT(g), wB(g);
    PartitionDown( w, wT, wB, part.index );
    DistMatrix<F> QL(g), QR(g);
    PartitionRight( Q, QL, QR, part.index );

    // Recurse on the two subproblems
    DistMatrix<F> ATLSub, ABRSub, ZTSub, ZBSub;
    DistMatrix<Real,VR,STAR> wTSub, wBSub;
    PushSubproblems
    ( ATL, ABR, ATLSub, ABRSub, wT, wB, wTSub, wBSub, ZTSub, ZBSub, 
      sdcCtrl.progress );
    if( ATLSub.Participating() )
        SDC( uplo, ATLSub, wTSub, ZTSub, sdcCtrl );
    if( ABRSub.Participating() )
        SDC( uplo, ABRSub, wBSub, ZBSub, sdcCtrl );

    // Pull the results back to this grid
    DistMatrix<F> ZT(g), ZB(g);
    PullSubproblems
    ( ATL, ABR, ATLSub, ABRSub, wT, wB, wTSub, wBSub, ZT, ZB, ZTSub, ZBSub );

    // Update the eigen vectors
    auto G( QL );
    Gemm( NORMAL, NORMAL, F(1), G, ZT, QL );
    G = QR;
    Gemm( NORMAL, NORMAL, F(1), G, ZB, QR );
}

} // namespace herm_eig
} // namespace elem

#endif // ifndef ELEM_HERMITIANEIG_SDC_HPP
