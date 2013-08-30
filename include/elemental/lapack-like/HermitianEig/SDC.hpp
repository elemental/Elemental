/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_HERMITIANEIG_SDC_HPP
#define ELEM_LAPACK_HERMITIANEIG_SDC_HPP

#include "elemental/blas-like/level1/Axpy.hpp"
#include "elemental/blas-like/level1/Scale.hpp"
#include "elemental/blas-like/level1/UpdateDiagonal.hpp"
#include "elemental/lapack-like/Median.hpp"
#include "elemental/lapack-like/Norm/Infinity.hpp"
#include "elemental/lapack-like/Norm/One.hpp"
#include "elemental/lapack-like/Polar/QDWH.hpp"
#include "elemental/lapack-like/QR.hpp"
#include "elemental/matrices/Haar.hpp"

#include "elemental/lapack-like/Schur.hpp"

// TODO: Reference to Yuji's work

namespace elem {
namespace hermitian_eig {

using elem::schur::ComputePartition;

// TODO: Exploit symmetry in A := Q^H A Q. Routine for A := X^H A X?

// G should be a rational function of A. If returnQ=true, G will be set to
// the computed unitary matrix upon exit.
template<typename F>
inline ValueInt<BASE(F)>
QDWHDivide( UpperOrLower uplo, Matrix<F>& A, Matrix<F>& G, bool returnQ=false )
{
#ifndef RELEASE
    CallStackEntry cse("hermitian_eig::QDWHDivide");
#endif
    typedef BASE(F) Real;
    const Int n = A.Height();

    // G := sgn(G)
    // G := 1/2 ( G + I )
    hermitian_polar::QDWH( uplo, G ); 
    UpdateDiagonal( G, F(1) );
    Scale( F(1)/F(2), G );

    // Compute the pivoted QR decomposition of the spectral projection 
    Matrix<F> t;
    Matrix<Int> p;
    elem::QR( G, t, p );

    // A := Q^H A Q
    MakeHermitian( uplo, A );
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
QDWHDivide
( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<F>& G, bool returnQ=false )
{
#ifndef RELEASE
    CallStackEntry cse("hermitian_eig::QDWHDivide");
#endif
    typedef BASE(F) Real;
    const Grid& g = A.Grid();
    const Int n = A.Height();

    // G := sgn(G)
    // G := 1/2 ( G + I )
    hermitian_polar::QDWH( uplo, G ); 
    UpdateDiagonal( G, F(1) );
    Scale( F(1)/F(2), G );

    // Compute the pivoted QR decomposition of the spectral projection 
    DistMatrix<F,MD,STAR> t(g);
    DistMatrix<Int,VR,STAR> p(g);
    elem::QR( G, t, p );

    // A := Q^H A Q
    MakeHermitian( uplo, A );
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
( UpperOrLower uplo, Matrix<F>& A, Matrix<F>& G, 
  bool returnQ=false, Int maxIts=1, BASE(F) relTol=0 )
{
#ifndef RELEASE
    CallStackEntry cse("hermitian_eig::RandomizedSignDivide");
#endif
    typedef BASE(F) Real;
    const Int n = A.Height();
    MakeHermitian( uplo, A );
    const Real oneA = OneNorm( A );
    if( relTol == Real(0) )
        relTol = 500*n*lapack::MachineEpsilon<Real>();

    // S := sgn(G)
    // S := 1/2 ( S + I )
    auto S( G );
    hermitian_polar::QDWH( uplo, S ); 
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
( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<F>& G, 
  bool returnQ=false, Int maxIts=1, BASE(F) relTol=0 )
{
#ifndef RELEASE
    CallStackEntry cse("hermitian_eig::RandomizedSignDivide");
#endif
    typedef BASE(F) Real;
    const Grid& g = A.Grid();
    const Int n = A.Height();
    MakeHermitian( uplo, A );
    const Real oneA = OneNorm( A );
    if( relTol == Real(0) )
        relTol = 500*n*lapack::MachineEpsilon<Real>();

    // S := sgn(G)
    // S := 1/2 ( S + I )
    auto S( G );
    hermitian_polar::QDWH( uplo, S ); 
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

template<typename F>
inline ValueInt<BASE(F)>
SpectralDivide
( UpperOrLower uplo, Matrix<F>& A, 
  Int maxInnerIts=1, Int maxOuterIts=10, BASE(F) relTol=0 )
{
#ifndef RELEASE
    CallStackEntry cse("hermitian_eig::SpectralDivide");
#endif
    typedef BASE(F) Real;
    const Int n = A.Height();
    MakeHermitian( uplo, A );
    const auto median = Median(A.GetRealPartOfDiagonal());
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
        const Real shift = SampleBall<Real>(-median.value,spread);

        G = A;
        UpdateDiagonal( G, shift );

        //part = SignDivide( uplo, A, G );
        part = RandomizedSignDivide( uplo, A, G, false, maxInnerIts, relTol );

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
inline ValueInt<BASE(F)>
SpectralDivide
( UpperOrLower uplo, Matrix<F>& A, Matrix<F>& Q, 
  Int maxInnerIts=1, Int maxOuterIts=10, BASE(F) relTol=0 )
{
#ifndef RELEASE
    CallStackEntry cse("hermitian_eig::SpectralDivide");
#endif
    typedef BASE(F) Real;
    const Int n = A.Height();
    MakeHermitian( uplo, A );
    const auto median = Median(A.GetRealPartOfDiagonal());
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
        const Real shift = SampleBall<Real>(-median.value,spread);

        Q = A;
        UpdateDiagonal( Q, shift );

        //part = SignDivide( uplo, A, Q, true );
        part = RandomizedSignDivide( uplo, A, Q, true, maxInnerIts, relTol );

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
inline ValueInt<BASE(F)>
SpectralDivide
( UpperOrLower uplo, DistMatrix<F>& A, 
  Int maxInnerIts=1, Int maxOuterIts=10, BASE(F) relTol=0 )
{
#ifndef RELEASE
    CallStackEntry cse("hermitian_eig::SpectralDivide");
#endif
    typedef BASE(F) Real;
    const Int n = A.Height();
    MakeHermitian( uplo, A );
    const auto median = Median(A.GetRealPartOfDiagonal());
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
        Real shift = SampleBall<Real>(-median.value,spread);
        mpi::Broadcast( shift, 0, A.Grid().VCComm() );

        G = A;
        UpdateDiagonal( G, shift );

        //part = SignDivide( uplo, A, G );
        part = RandomizedSignDivide( uplo, A, G, false, maxInnerIts, relTol );

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
inline ValueInt<BASE(F)>
SpectralDivide
( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<F>& Q, 
  Int maxInnerIts=1, Int maxOuterIts=10, BASE(F) relTol=0 )
{
#ifndef RELEASE
    CallStackEntry cse("hermitian_eig::SpectralDivide");
#endif
    typedef BASE(F) Real;
    const Int n = A.Height();
    MakeHermitian( uplo, A );
    const Real infNorm = InfinityNorm(A);
    const auto median = Median(A.GetRealPartOfDiagonal());
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
        Real shift = SampleBall<Real>(-median.value,spread);
        mpi::Broadcast( shift, 0, A.Grid().VCComm() );

        Q = A;
        UpdateDiagonal( Q, shift );

        //part = SignDivide( uplo, A, Q, true );
        part = RandomizedSignDivide( uplo, A, Q, true, maxInnerIts, relTol );

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
( UpperOrLower uplo, Matrix<F>& A, Int cutoff=256, 
  Int maxInnerIts=1, Int maxOuterIts=10, BASE(F) relTol=0 )
{
#ifndef RELEASE
    CallStackEntry cse("hermitian_eig::SDC");
#endif
    typedef BASE(F) Real;
    const Int n = A.Height();
    if( n <= cutoff )
    {
        Matrix<Real> w;
        HermitianEig( uplo, A, w );
        MakeZeros( A );
        A.SetRealPartOfDiagonal( w );
        return;
    }

    // Perform this level's split
    const auto part = 
        SpectralDivide( uplo, A, maxInnerIts, maxOuterIts, relTol );
    Matrix<F> ATL, ATR,
              ABL, ABR;
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, part.index );

    // Recurse on the two subproblems
    SDC( uplo, ATL, cutoff, maxInnerIts, maxOuterIts, relTol );
    SDC( uplo, ABR, cutoff, maxInnerIts, maxOuterIts, relTol );
}

template<typename F>
inline void
SDC
( UpperOrLower uplo, Matrix<F>& A, Matrix<F>& Q, 
  Int cutoff=256, Int maxInnerIts=1, Int maxOuterIts=10, BASE(F) relTol=0 )
{
#ifndef RELEASE
    CallStackEntry cse("hermitian_eig::SDC");
#endif
    typedef BASE(F) Real;
    const Int n = A.Height();
    if( n <= cutoff )
    {
        Matrix<Real> w;
        HermitianEig( uplo, A, w, Q );
        MakeZeros( A );
        A.SetRealPartOfDiagonal( w );
        return;
    }

    // Perform this level's split
    const auto part = 
        SpectralDivide( uplo, A, Q, maxInnerIts, maxOuterIts, relTol );
    Matrix<F> ATL, ATR,
              ABL, ABR;
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, part.index );
    Matrix<F> QL, QR;
    PartitionRight( Q, QL, QR, part.index );

    // Recurse on the top-left quadrant and update eigenvectors
    Matrix<F> Z;
    SDC( uplo, ATL, Z, cutoff, maxInnerIts, maxOuterIts, relTol );
    auto G( QL );
    Gemm( NORMAL, NORMAL, F(1), G, Z, QL );

    // Recurse on the bottom-right quadrant and update eigenvectors
    SDC( uplo, ABR, Z, cutoff, maxInnerIts, maxOuterIts, relTol );
    G = QR;
    Gemm( NORMAL, NORMAL, F(1), G, Z, QR );
}

template<typename F>
inline void
SDC
( UpperOrLower uplo, DistMatrix<F>& A, Int cutoff=256, 
  Int maxInnerIts=1, Int maxOuterIts=10, BASE(F) relTol=0 )
{
#ifndef RELEASE
    CallStackEntry cse("hermitian_eig::SDC");
#endif
    typedef BASE(F) Real;
    const Grid& g = A.Grid();
    const Int n = A.Height();
    if( n <= cutoff )
    {
        DistMatrix<Real,VR,STAR> w(A.Grid());
        HermitianEig( uplo, A, w );
        MakeZeros( A );
        DistMatrix<Real,MD,STAR> w_MD_STAR(A.Grid());
        w_MD_STAR.AlignWithDiag( A );
        w_MD_STAR = w;
        A.SetRealPartOfDiagonal( w_MD_STAR );
        return;
    }

    // Perform this level's split
    const auto part = 
        SpectralDivide( uplo, A, maxInnerIts, maxOuterIts, relTol );
    DistMatrix<F> ATL(g), ATR(g),
                  ABL(g), ABR(g);
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, part.index );

    // Recurse on the two subproblems
    SDC( uplo, ATL, cutoff, maxInnerIts, maxOuterIts, relTol );
    SDC( uplo, ABR, cutoff, maxInnerIts, maxOuterIts, relTol );
}

template<typename F>
inline void
SDC
( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<F>& Q, 
  Int cutoff=256, Int maxInnerIts=1, Int maxOuterIts=10, BASE(F) relTol=0 )
{
#ifndef RELEASE
    CallStackEntry cse("hermitian_eig::SDC");
#endif
    typedef BASE(F) Real;
    const Grid& g = A.Grid();
    const Int n = A.Height();
    if( n <= cutoff )
    {
        DistMatrix<Real,VR,STAR> w(A.Grid());
        HermitianEig( uplo, A, w, Q );
        MakeZeros( A );
        DistMatrix<Real,MD,STAR> w_MD_STAR(A.Grid());
        w_MD_STAR.AlignWithDiag( A );
        w_MD_STAR = w;
        A.SetRealPartOfDiagonal( w_MD_STAR );
        return;
    }

    // Perform this level's split
    const auto part = 
        SpectralDivide( uplo, A, Q, maxInnerIts, maxOuterIts, relTol );
    DistMatrix<F> ATL(g), ATR(g),
                  ABL(g), ABR(g);
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, part.index );
    DistMatrix<F> QL(g), QR(g);
    PartitionRight( Q, QL, QR, part.index );

    // Recurse on the two subproblems
    DistMatrix<F> ZT(g), ZB(g);
    SDC( uplo, ATL, ZT, cutoff, maxInnerIts, maxOuterIts, relTol );
    SDC( uplo, ABR, ZB, cutoff, maxInnerIts, maxOuterIts, relTol );

    // Update the eigen vectors
    auto G( QL );
    Gemm( NORMAL, NORMAL, F(1), G, ZT, QL );
    G = QR;
    Gemm( NORMAL, NORMAL, F(1), G, ZB, QR );
}

} // namespace hermitian_eig
} // namespace elem

#endif // ifndef ELEM_LAPACK_HERMITIANEIG_SDC_HPP
