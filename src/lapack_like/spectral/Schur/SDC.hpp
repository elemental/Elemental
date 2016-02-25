/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SCHUR_SDC_HPP
#define EL_SCHUR_SDC_HPP

// See Z. Bai, J. Demmel, J. Dongarra, A. Petitet, H. Robinson, and K. Stanley's
// "The spectral decomposition of nonsymmetric matrices on distributed memory
// parallel computers". Currently available at:
// www.netlib.org/lapack/lawnspdf/lawn91.pdf
//
// as well as the improved version, which avoids pivoted QR, in J. Demmel, 
// I. Dumitriu, and O. Holtz, "Fast linear algebra is stable", 2007.
// www.netlib.org/lapack/lawnspdf/lawn186.pdf

namespace El {

namespace schur {

template<typename F>
inline ValueInt<Base<F>>
ComputePartition( Matrix<F>& A )
{
    DEBUG_ONLY(CSE cse("schur::ComputePartition"))
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
    vector<Real> colSums(n-1,0), rowSums(n-1,0);
    for( Int j=0; j<n-1; ++j )
        for( Int i=j+1; i<n; ++i )
            colSums[j] += Abs( A.Get(i,j) ); 
    for( Int i=1; i<n-1; ++i )
        for( Int j=0; j<i; ++j )
            rowSums[i-1] += Abs( A.Get(i,j) );

    // Compute the list of norms and its minimum value/index
    ValueInt<Real> part;
    vector<Real> norms(n-1);
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
    DEBUG_ONLY(CSE cse("schur::ComputePartition"))
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
    vector<Real> colSums(n-1,0), rowSums(n-1,0);
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
    vector<Real> norms(n-1);
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
( Matrix<F>& A,
  Matrix<F>& G,
  bool returnQ,
  const SDCCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("schur::SignDivide"))

    // G := sgn(G)
    // G := 1/2 ( G + I )
    Sign( G, ctrl.signCtrl );
    ShiftDiagonal( G, F(1) );
    G *= F(1)/F(2);

    // Compute the pivoted QR decomposition of the spectral projection 
    Matrix<F> t;
    Matrix<Base<F>> d;
    Permutation Omega;
    El::QR( G, t, d, Omega );

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
( DistMatrix<F>& A,
  DistMatrix<F>& G,
  bool returnQ, 
  const SDCCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("schur::SignDivide"))
    const Grid& g = A.Grid();

    // G := sgn(G)
    // G := 1/2 ( G + I )
    Sign( G, ctrl.signCtrl );
    ShiftDiagonal( G, F(1) );
    G *= F(1)/F(2);

    // Compute the pivoted QR decomposition of the spectral projection 
    DistMatrix<F,MD,STAR> t(g);
    DistMatrix<Base<F>,MD,STAR> d(g);
    DistPermutation Omega(g);
    El::QR( G, t, d, Omega );

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
( Matrix<F>& A,
  Matrix<F>& G,
  bool returnQ,
  const SDCCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("schur::RandomizedSignDivide"))
    typedef Base<F> Real;
    const Int n = A.Height();
    const Real oneA = OneNorm( A );
    Real tol = ctrl.tol;
    if( tol == Real(0) )
        tol = 500*n*limits::Epsilon<Real>();

    // S := sgn(G)
    // S := 1/2 ( S + I )
    auto S( G );
    Sign( S, ctrl.signCtrl );
    ShiftDiagonal( S, F(1) );
    S *= F(1)/F(2);

    ValueInt<Real> part;
    Matrix<F> V, B, t;
    Matrix<Base<F>> d;
    Int it=0;
    while( it < ctrl.maxInnerIts )
    {
        G = S;

        // Compute the RURV of the spectral projector
        ImplicitHaar( V, t, d, n );
        qr::ApplyQ( RIGHT, NORMAL, V, t, d, G );
        El::QR( G, t, d );

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
        if( part.value <= tol || it == ctrl.maxInnerIts )
            break;
        else
            A = V;
    }
    return part;
}

template<typename F>
inline ValueInt<Base<F>>
RandomizedSignDivide
( DistMatrix<F>& A,
  DistMatrix<F>& G,
  bool returnQ,
  const SDCCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("schur::RandomizedSignDivide"))
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    const Int n = A.Height();
    const Real oneA = OneNorm( A );
    Real tol = ctrl.tol;
    if( tol == Real(0) )
        tol = 500*n*limits::Epsilon<Real>();

    // S := sgn(G)
    // S := 1/2 ( S + I )
    auto S( G );
    Sign( S, ctrl.signCtrl );
    ShiftDiagonal( S, F(1) );
    S *= F(1)/F(2);

    ValueInt<Real> part;
    DistMatrix<F> V(g), B(g);
    DistMatrix<F,MD,STAR> t(g);
    DistMatrix<Base<F>,MD,STAR> d(g);
    Int it=0;
    while( it < ctrl.maxInnerIts )
    {
        G = S;

        // Compute the RURV of the spectral projector
        ImplicitHaar( V, t, d, n );
        qr::ApplyQ( RIGHT, NORMAL, V, t, d, G );
        El::QR( G, t, d );

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
        if( part.value <= tol || it == ctrl.maxInnerIts )
            break;
        else
            A = V;
    }
    return part;
}

template<typename Real>
inline ValueInt<Real>
SpectralDivide( Matrix<Real>& A, const SDCCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("schur::SpectralDivide"))
    const Int n = A.Height();
    const ValueInt<Real> median = Median(GetDiagonal(A));
    const Real infNorm = InfinityNorm(A);
    const Real eps = limits::Epsilon<Real>();
    Real tol = ctrl.tol;
    if( tol == Real(0) )
        tol = 500*n*eps;
    const Real spread = ctrl.spreadFactor*infNorm;

    Int it=0;
    ValueInt<Real> part;
    part.value = 2*tol; // initialize with unacceptable value
    Matrix<Real> G, ACopy;
    if( ctrl.maxOuterIts > 1 )
        ACopy = A;
    while( it < ctrl.maxOuterIts )
    {
        ++it;
        const Real shift = SampleBall<Real>(-median.value,spread);

        G = A;
        ShiftDiagonal( G, shift );

        if( ctrl.progress )
            cout << "chose shift=" << shift << " using -median.value="
                 << -median.value << " and spread=" << spread << endl;

        try
        {
            if( ctrl.random )
                part = RandomizedSignDivide( A, G, false, ctrl );
            else
                part = SignDivide( A, G, false, ctrl );

            if( part.value <= tol )
            {
                if( ctrl.progress )
                    cout << "Converged during outer iter " << it-1 << endl;
                break;
            }
            else if( ctrl.progress )
                cout << "part.value=" << part.value << " was greater than "
                     << tol << " during outer iter " << it-1 << endl;
        } 
        catch( SingularMatrixException& e ) 
        { 
            if( ctrl.progress )
                cout << "Caught singular matrix in outer iter " << it-1 << endl;
        }
        if( it != ctrl.maxOuterIts )
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
( Matrix<Complex<Real>>& A, const SDCCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("schur::SpectralDivide"))
    typedef Complex<Real> F;
    const Int n = A.Height();
    const Real infNorm = InfinityNorm(A);
    const Real eps = limits::Epsilon<Real>();
    Real tol = ctrl.tol;
    if( tol == Real(0) )
        tol = 500*n*eps;
    const Real spread = ctrl.spreadFactor*infNorm;

    Int it=0;
    ValueInt<Real> part;
    part.value = 2*tol; // initialize with unacceptable value
    Matrix<F> G, ACopy;
    if( ctrl.maxOuterIts > 1 )
        ACopy = A;
    while( it < ctrl.maxOuterIts )
    {
        ++it;
        const Real angle = SampleUniform<Real>(0,2*Pi<Real>());
        const F gamma = F(Cos(angle),Sin(angle));
        G = A;
        G *= gamma;

        const auto median = Median(GetRealPartOfDiagonal(G));
        const F shift = SampleBall<F>(-median.value,spread);
        ShiftDiagonal( G, shift );

        if( ctrl.progress )
            cout << "chose gamma=" << gamma << " and shift=" << shift 
                 << " using -median.value=" << -median.value 
                 << " and spread=" << spread << endl;

        try
        {
            if( ctrl.random )
                part = RandomizedSignDivide( A, G, false, ctrl );
            else
                part = SignDivide( A, G, false, ctrl );

            if( part.value <= tol )
            {
                if( ctrl.progress )
                    cout << "Converged during outer iter " << it-1 << endl;
                break;
            }
            else if( ctrl.progress )
                cout << "part.value=" << part.value << " was greater than "
                     << tol << " during outer iter " << it-1 << endl;
        } 
        catch( SingularMatrixException& e ) 
        {
            if( ctrl.progress )
                cout << "Caught singular matrix in outer iter " << it-1 << endl;
        }
        if( it != ctrl.maxOuterIts )
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
( Matrix<Real>& A,
  Matrix<Real>& Q,
  const SDCCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("schur::SpectralDivide"))
    const Int n = A.Height();
    const auto median = Median(GetDiagonal(A));
    const Real infNorm = InfinityNorm(A);
    const Real eps = limits::Epsilon<Real>();
    Real tol = ctrl.tol;
    if( tol == Real(0) )
        tol = 500*n*eps;
    const Real spread = ctrl.spreadFactor*infNorm;

    Int it=0;
    ValueInt<Real> part;
    part.value = 2*tol; // initialize with unacceptable value
    Matrix<Real> ACopy;
    if( ctrl.maxOuterIts > 1 )
        ACopy = A;
    while( it < ctrl.maxOuterIts )
    {
        ++it;
        const Real shift = SampleBall<Real>(-median.value,spread);

        Q = A;
        ShiftDiagonal( Q, shift );

        if( ctrl.progress )
            cout << "chose shift=" << shift << " using -median.value="
                 << -median.value << " and spread=" << spread << endl;

        try
        {
            if( ctrl.random )
                part = RandomizedSignDivide( A, Q, true, ctrl );
            else
                part = SignDivide( A, Q, true, ctrl );

            if( part.value <= tol )
            {
                if( ctrl.progress )
                    cout << "Converged during outer iter " << it-1 << endl;
                break;
            }
            else if( ctrl.progress )
                cout << "part.value=" << part.value << " was greater than "
                     << tol << " during outer iter " << it-1 << endl;
        } 
        catch( SingularMatrixException& e ) 
        { 
            if( ctrl.progress )
                cout << "Caught singular matrix in outer iter " << it-1 << endl;
        }
        if( it != ctrl.maxOuterIts )
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
( Matrix<Complex<Real>>& A,
  Matrix<Complex<Real>>& Q, 
  const SDCCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("schur::SpectralDivide"))
    typedef Complex<Real> F;
    const Int n = A.Height();
    const Real infNorm = InfinityNorm(A);
    const Real eps = limits::Epsilon<Real>();
    Real tol = ctrl.tol;
    if( tol == Real(0) )
        tol = 500*n*eps;
    const Real spread = ctrl.spreadFactor*infNorm;

    Int it=0;
    ValueInt<Real> part;
    part.value = 2*tol; // initialize with unacceptable value
    Matrix<F> ACopy;
    if( ctrl.maxOuterIts > 1 )
        ACopy = A;
    while( it < ctrl.maxOuterIts )
    {
        ++it;
        const Real angle = SampleUniform<Real>(0,2*Pi<Real>());
        const F gamma = F(Cos(angle),Sin(angle));
        Q = A;
        Q *= gamma;

        const auto median = Median(GetRealPartOfDiagonal(Q));
        const F shift = SampleBall<F>(-median.value,spread);
        ShiftDiagonal( Q, shift );

        if( ctrl.progress )
            cout << "chose gamma=" << gamma << " and shift=" << shift 
                 << " using -median.value=" << -median.value 
                 << " and spread=" << spread << endl;

        try
        {
            if( ctrl.random )
                part = RandomizedSignDivide( A, Q, true, ctrl );
            else
                part = SignDivide( A, Q, true, ctrl );

            if( part.value <= tol )
            {
                if( ctrl.progress )
                    cout << "Converged during outer iter " << it-1 << endl;
                break;
            }
            else if( ctrl.progress )
                cout << "part.value=" << part.value << " was greater than "
                     << tol << " during outer iter " << it-1 << endl;
        } 
        catch( SingularMatrixException& e ) 
        { 
            if( ctrl.progress )
                cout << "Caught singular matrix in outer iter " << it-1 << endl;
        }
        if( it != ctrl.maxOuterIts )
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
SpectralDivide( DistMatrix<Real>& A, const SDCCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("schur::SpectralDivide"))
    const Int n = A.Height();
    const auto median = Median(GetDiagonal(A));
    const Real infNorm = InfinityNorm(A);
    const Real eps = limits::Epsilon<Real>();
    Real tol = ctrl.tol;
    if( tol == Real(0) )
        tol = 500*n*eps;
    const Real spread = ctrl.spreadFactor*infNorm;

    Int it=0;
    ValueInt<Real> part;
    part.value = 2*tol; // initialize with unacceptable value
    const Grid& g = A.Grid();
    DistMatrix<Real> ACopy(g), G(g);
    if( ctrl.maxOuterIts > 1 )
        ACopy = A;
    while( it < ctrl.maxOuterIts )
    {
        ++it;
        Real shift = SampleBall<Real>(-median.value,spread);
        mpi::Broadcast( shift, 0, g.VCComm() );

        G = A;
        ShiftDiagonal( G, shift );

        if( ctrl.progress && g.Rank() == 0 )
            cout << "chose shift=" << shift << " using -median.value="
                 << -median.value << " and spread=" << spread << endl;

        try
        {
            if( ctrl.random )
                part = RandomizedSignDivide( A, G, false, ctrl );
            else
                part = SignDivide( A, G, false, ctrl );

            if( part.value <= tol )
            {
                if( ctrl.progress && g.Rank() == 0 )
                    cout << "Converged during outer iter " << it-1 << endl;
                break;
            }
            else if( ctrl.progress && g.Rank() == 0 )
                cout << "part.value=" << part.value << " was greater than "
                     << tol << " during outer iter " << it-1 << endl;
        } 
        catch( SingularMatrixException& e ) 
        { 
            if( ctrl.progress && g.Rank() == 0 )
                cout << "Caught singular matrix in outer iter " << it-1 << endl;
        }
        if( it != ctrl.maxOuterIts )
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
( DistMatrix<Complex<Real>>& A, const SDCCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("schur::SpectralDivide"))
    typedef Complex<Real> F;
    const Int n = A.Height();
    const Real infNorm = InfinityNorm(A);
    const Real eps = limits::Epsilon<Real>();
    Real tol = ctrl.tol;
    if( tol == Real(0) )
        tol = 500*n*eps;
    const Real spread = ctrl.spreadFactor*infNorm;

    Int it=0;
    ValueInt<Real> part;
    part.value = 2*tol; // initialize with unacceptable value
    const Grid& g = A.Grid();
    DistMatrix<F> ACopy(g), G(g);
    if( ctrl.maxOuterIts > 1 )
        ACopy = A;
    while( it < ctrl.maxOuterIts )
    {
        ++it;
        const Real angle = SampleUniform<Real>(0,2*Pi<Real>());
        F gamma = F(Cos(angle),Sin(angle));
        mpi::Broadcast( gamma, 0, A.Grid().VCComm() );
        G = A;
        G *= gamma;

        const auto median = Median(GetRealPartOfDiagonal(G));
        F shift = SampleBall<F>(-median.value,spread);
        mpi::Broadcast( shift, 0, g.VCComm() );
        ShiftDiagonal( G, shift );

        if( ctrl.progress && g.Rank() == 0 )
            cout << "chose gamma=" << gamma << " and shift=" << shift 
                 << " using -median.value=" << -median.value 
                 << " and spread=" << spread << endl;

        try
        {
            if( ctrl.random )
                part = RandomizedSignDivide( A, G, false, ctrl );
            else
                part = SignDivide( A, G, false, ctrl );

            if( part.value <= tol )
            {
                if( ctrl.progress && g.Rank() == 0 )
                    cout << "Converged during outer iter " << it-1 << endl;
                break;
            }
            else if( ctrl.progress && g.Rank() == 0 )
                cout << "part.value=" << part.value << " was greater than "
                     << tol << " during outer iter " << it-1 << endl;
        } 
        catch( SingularMatrixException& e ) 
        { 
            if( ctrl.progress && g.Rank() == 0 )
                cout << "Caught singular matrix in outer iter " << it-1 << endl;
        }
        if( it != ctrl.maxOuterIts )
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
( DistMatrix<Real>& A,
  DistMatrix<Real>& Q,
  const SDCCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("schur::SpectralDivide"))
    const Int n = A.Height();
    const Real infNorm = InfinityNorm(A);
    const auto median = Median(GetDiagonal(A));
    const Real eps = limits::Epsilon<Real>();
    Real tol = ctrl.tol;
    if( tol == Real(0) )
        tol = 500*n*eps;
    const Real spread = ctrl.spreadFactor*infNorm;

    Int it=0;
    ValueInt<Real> part;
    part.value = 2*tol; // initialize with unacceptable value
    const Grid& g = A.Grid();
    DistMatrix<Real> ACopy(g);
    if( ctrl.maxOuterIts > 1 )
        ACopy = A;
    while( it < ctrl.maxOuterIts )
    {
        ++it;
        Real shift = SampleBall<Real>(-median.value,spread);
        mpi::Broadcast( shift, 0, g.VCComm() );

        Q = A;
        ShiftDiagonal( Q, shift );

        if( ctrl.progress && g.Rank() == 0 )
            cout << "chose shift=" << shift << " using -median.value=" 
                 << -median.value << " and spread=" << spread << endl;

        try
        {
            if( ctrl.random )
                part = RandomizedSignDivide( A, Q, true, ctrl );
            else
                part = SignDivide( A, Q, true, ctrl );

            if( part.value <= tol )
            {
                if( ctrl.progress && g.Rank() == 0 )
                    cout << "Converged during outer iter " << it-1 << endl;
                break;
            }
            else if( ctrl.progress && g.Rank() == 0 )
                cout << "part.value=" << part.value << " was greater than "
                     << tol << " during outer iter " << it-1 << endl;
        } 
        catch( SingularMatrixException& e ) 
        {
            if( ctrl.progress && g.Rank() == 0 )
                cout << "Caught singular matrix in outer iter " << it-1 << endl;
        }
        if( it != ctrl.maxOuterIts )
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
( DistMatrix<Complex<Real>>& A,
  DistMatrix<Complex<Real>>& Q,
  const SDCCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("schur::SpectralDivide"))
    typedef Complex<Real> F;
    const Int n = A.Height();
    const Real infNorm = InfinityNorm(A);
    const Real eps = limits::Epsilon<Real>();
    Real tol = ctrl.tol;
    if( tol == Real(0) )
        tol = 500*n*eps;
    const Real spread = ctrl.spreadFactor*infNorm;

    Int it=0;
    ValueInt<Real> part;
    part.value = 2*tol; // initialize with unacceptable value
    const Grid& g = A.Grid();
    DistMatrix<F> ACopy(g);
    if( ctrl.maxOuterIts > 1 )
        ACopy = A;
    while( it < ctrl.maxOuterIts )
    {
        ++it;
        const Real angle = SampleUniform<Real>(0,2*Pi<Real>());
        F gamma = F(Cos(angle),Sin(angle));
        mpi::Broadcast( gamma, 0, g.VCComm() );
        Q = A;
        Q *= gamma;

        const auto median = Median(GetRealPartOfDiagonal(Q));
        F shift = SampleBall<F>(-median.value,spread);
        mpi::Broadcast( shift, 0, g.VCComm() );
        ShiftDiagonal( Q, shift );

        if( ctrl.progress && g.Rank() == 0 )
            cout << "chose gamma=" << gamma << " and shift=" << shift 
                 << " using -median.value=" << -median.value 
                 << " and spread=" << spread << endl;

        try
        {
            if( ctrl.random )
                part = RandomizedSignDivide( A, Q, true, ctrl );
            else
                part = SignDivide( A, Q, true, ctrl );

            if( part.value <= tol )
            {
                if( ctrl.progress && g.Rank() == 0 )
                    cout << "Converged during outer iter " << it-1 << endl;
                break;
            }
            else if( ctrl.progress && g.Rank() == 0 )
                cout << "part.value=" << part.value << " was greater than "
                     << tol << " during outer iter " << it-1 << endl;
        } 
        catch( SingularMatrixException& e ) 
        { 
            if( ctrl.progress && g.Rank() == 0 )
                cout << "Caught singular matrix in outer iter " << it-1 << endl;
        }
        if( it != ctrl.maxOuterIts )
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
( Matrix<F>& A,
  Matrix<Complex<Base<F>>>& w, 
  const SDCCtrl<Base<F>> ctrl=SDCCtrl<Base<F>>() )
{
    DEBUG_ONLY(CSE cse("schur::SDC"))
    const Int n = A.Height();
    w.Resize( n, 1 );
    if( n <= ctrl.cutoff )
    {
        if( ctrl.progress )
            Output(n," <= ",ctrl.cutoff,": switching to QR algorithm");
        Schur( A, w, false );
        return;
    }

    // Perform this level's split
    if( ctrl.progress )
        Output("Splitting ",n," x ",n," matrix");
    const auto part = SpectralDivide( A, ctrl );
    auto ind1 = IR(0,part.index);
    auto ind2 = IR(part.index,n);

    auto ATL = A( ind1, ind1 );
    auto ATR = A( ind1, ind2 );
    auto ABL = A( ind2, ind1 );
    auto ABR = A( ind2, ind2 );

    auto wT = w( ind1, ALL );
    auto wB = w( ind2, ALL );

    Zero( ABL );

    // Recurse on the two subproblems
    if( ctrl.progress )
        Output("Recursing on ",ATL.Height()," x ",ATL.Width()," top-left");
    SDC( ATL, wT, ctrl );
    if( ctrl.progress )
        Output("Recursing on ",ABR.Height()," x ",ABR.Width()," bottom-right");
    SDC( ABR, wB, ctrl );
}

template<typename F>
inline void
SDC
( Matrix<F>& A,
  Matrix<Complex<Base<F>>>& w,
  Matrix<F>& Q, 
  bool fullTriangle=true,
  const SDCCtrl<Base<F>> ctrl=SDCCtrl<Base<F>>() )
{
    DEBUG_ONLY(CSE cse("schur::SDC"))
    const Int n = A.Height();
    w.Resize( n, 1 );
    Q.Resize( n, n );
    if( n <= ctrl.cutoff )
    {
        if( ctrl.progress )
            Output(n," <= ",ctrl.cutoff,": switching to QR algorithm");
        Schur( A, w, Q, fullTriangle );
        return;
    }

    // Perform this level's split
    if( ctrl.progress )
        Output("Splitting ",n," x ",n," matrix");
    const auto part = SpectralDivide( A, Q, ctrl );

    auto ind1 = IR(0,part.index);
    auto ind2 = IR(part.index,n);

    auto ATL = A( ind1, ind1 );
    auto ATR = A( ind1, ind2 );
    auto ABL = A( ind2, ind1 );
    auto ABR = A( ind2, ind2 );

    auto QL = Q( ALL, ind1 );
    auto QR = Q( ALL, ind2 );

    auto wT = w( ind1, ALL );
    auto wB = w( ind2, ALL );

    Zero( ABL );

    // Recurse on the top-left quadrant and update Schur vectors and ATR
    if( ctrl.progress )
        Output("Recursing on ",ATL.Height()," x ",ATL.Width()," top-left");
    Matrix<F> Z;
    SDC( ATL, wT, Z, fullTriangle, ctrl );
    if( ctrl.progress )
        Output("Left subproblem update");
    auto G( QL );
    Gemm( NORMAL, NORMAL, F(1), G, Z, QL );
    if( fullTriangle )
        Gemm( ADJOINT, NORMAL, F(1), Z, ATR, G );

    // Recurse on the bottom-right quadrant and update Schur vectors and ATR
    if( ctrl.progress )
        Output("Recursing on ",ABR.Height()," x ",ABR.Width()," bottom-right");
    SDC( ABR, wB, Z, fullTriangle, ctrl );
    if( ctrl.progress )
        Output("Right subproblem update");
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
( int nLeft,
  int nRight,
  const Grid& grid, 
  const Grid*& leftGrid,
  const Grid*& rightGrid,
  bool progress=false )
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
            cout << "leftWork/rightWork=" << leftWork/rightWork 
                 << ", so the grid was not split" << endl;
    }
    else
    {
        // Split the grid in half (powers-of-two remain so)
        const Int p = grid.Size();
        const Int pLeft = p/2;
        const Int pRight = p-pLeft;
        vector<int> leftRanks(pLeft), rightRanks(pRight);
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
            cout << "leftWork/rightWork=" << leftWork/rightWork 
                 << ", so split " << p << " processes into " 
                 << rLeft << " x " << pLeft/rLeft << " and "
                 << rRight << " x " << pRight/rRight << " grids" << endl;
        leftGrid = new Grid( grid.VCComm(), leftGroup, rLeft );
        rightGrid = new Grid( grid.VCComm(), rightGroup, rRight );
        mpi::Free( leftGroup );
        mpi::Free( rightGroup );
    }
}

template<typename F,typename EigType>
inline void PushSubproblems
( DistMatrix<F>& ATL,
  DistMatrix<F>& ABR, 
  DistMatrix<F>& ATLSub,
  DistMatrix<F>& ABRSub,
  DistMatrix<EigType,VR,STAR>& wT,    
  DistMatrix<EigType,VR,STAR>& wB,
  DistMatrix<EigType,VR,STAR>& wTSub, 
  DistMatrix<EigType,VR,STAR>& wBSub,
  bool progress=false )
{
    DEBUG_ONLY(CSE cse("schur::PushSubproblems"))
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
        Output("Pushing ATL and ABR");
    ATLSub = ATL;
    ABRSub = ABR;
}

template<typename F,typename EigType>
inline void PullSubproblems
( DistMatrix<F>& ATL,
  DistMatrix<F>& ABR,
  DistMatrix<F>& ATLSub,
  DistMatrix<F>& ABRSub,
  DistMatrix<EigType,VR,STAR>& wT,    
  DistMatrix<EigType,VR,STAR>& wB,
  DistMatrix<EigType,VR,STAR>& wTSub, 
  DistMatrix<EigType,VR,STAR>& wBSub,
  bool progress=false )
{
    DEBUG_ONLY(CSE cse("schur::PullSubproblems"))
    const Grid& grid = ATL.Grid();
    const bool sameGrid = ( wT.Grid() == wTSub.Grid() );

    if( progress && grid.Rank() == 0 )
        Output("Pulling ATL and ABR");
    ATL = ATLSub;
    ABR = ABRSub;

    // This section is a hack since no inter-grid redistributions exist for 
    // [VR,* ] distributions yet
    if( progress && grid.Rank() == 0 )
        Output("Pulling wT and wB");
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
( ElementalMatrix<F>& APre,
  ElementalMatrix<Complex<Base<F>>>& wPre, 
  const SDCCtrl<Base<F>> ctrl=SDCCtrl<Base<F>>() )
{
    DEBUG_ONLY(
      CSE cse("schur::SDC");
      AssertSameGrids( APre, wPre );
    )
    typedef Base<F> Real;
    typedef Complex<Real> C;

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    DistMatrixWriteProxy<C,C,VR,STAR> wProx( wPre );
    auto& A = AProx.Get();
    auto& w = wProx.Get();

    const Grid& g = A.Grid();
    const Int n = A.Height();
    w.Resize( n, 1 );
    if( A.Grid().Size() == 1 )
    {
        if( ctrl.progress && g.Rank() == 0 )
            Output("One process: using QR algorithm");
        Schur( A.Matrix(), w.Matrix(), false );
        return;
    }
    if( n <= ctrl.cutoff )
    {
        if( ctrl.progress && g.Rank() == 0 )
            Output(n," <= ",ctrl.cutoff,": using QR algorithm"); 
        DistMatrix<F,CIRC,CIRC> A_CIRC_CIRC( A );
        DistMatrix<Complex<Base<F>>,CIRC,CIRC> w_CIRC_CIRC( w );
        if( A_CIRC_CIRC.CrossRank() == A_CIRC_CIRC.Root() )
            Schur( A_CIRC_CIRC.Matrix(), w_CIRC_CIRC.Matrix(), false );
        A = A_CIRC_CIRC;
        w = w_CIRC_CIRC;
        return;
    }

    // Perform this level's split
    if( ctrl.progress && g.Rank() == 0 )
        Output("Splitting ",n," x ",n," matrix");
    const auto part = SpectralDivide( A, ctrl );
    auto ind1 = IR(0,part.index);
    auto ind2 = IR(part.index,n);

    auto ATL = A( ind1, ind1 );
    auto ATR = A( ind1, ind2 );
    auto ABL = A( ind2, ind1 );
    auto ABR = A( ind2, ind2 );

    auto wT = w( ind1, ALL );
    auto wB = w( ind2, ALL );

    Zero( ABL );

    if( ctrl.progress && g.Rank() == 0 )
        Output("Pushing subproblems");
    DistMatrix<F> ATLSub, ABRSub;
    DistMatrix<Complex<Base<F>>,VR,STAR> wTSub, wBSub;
    PushSubproblems
    ( ATL, ABR, ATLSub, ABRSub, wT, wB, wTSub, wBSub, ctrl.progress );
    if( ATLSub.Participating() )
        SDC( ATLSub, wTSub, ctrl );
    if( ABRSub.Participating() )
        SDC( ABRSub, wBSub, ctrl );
    if( ctrl.progress && g.Rank() == 0 )
        Output("Pulling subproblems");
    PullSubproblems
    ( ATL, ABR, ATLSub, ABRSub, wT, wB, wTSub, wBSub, ctrl.progress );
}

template<typename F,typename EigType>
inline void PushSubproblems
( DistMatrix<F>& ATL,
  DistMatrix<F>& ABR, 
  DistMatrix<F>& ATLSub,
  DistMatrix<F>& ABRSub,
  DistMatrix<EigType,VR,STAR>& wT,    
  DistMatrix<EigType,VR,STAR>& wB,
  DistMatrix<EigType,VR,STAR>& wTSub, 
  DistMatrix<EigType,VR,STAR>& wBSub,
  DistMatrix<F>& ZTSub,
  DistMatrix<F>& ZBSub,
  bool progress=false )
{
    DEBUG_ONLY(CSE cse("schur::PushSubproblems"))
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
        Output("Pushing ATLSub");
    ATLSub = ATL;
    if( progress && grid.Rank() == 0 )
        Output("Pushing ABRSub");
    ABRSub = ABR;
}

template<typename F,typename EigType>
inline void PullSubproblems
( DistMatrix<F>& ATL,
  DistMatrix<F>& ABR,
  DistMatrix<F>& ATLSub,
  DistMatrix<F>& ABRSub,
  DistMatrix<EigType,VR,STAR>& wT,    
  DistMatrix<EigType,VR,STAR>& wB,
  DistMatrix<EigType,VR,STAR>& wTSub, 
  DistMatrix<EigType,VR,STAR>& wBSub,
  DistMatrix<F>& ZT,
  DistMatrix<F>& ZB,
  DistMatrix<F>& ZTSub,
  DistMatrix<F>& ZBSub,
  bool progress=false )
{
    DEBUG_ONLY(CSE cse("schur::PullSubproblems"))
    const Grid& grid = ATL.Grid();
    const bool sameGrid = ( wT.Grid() == wTSub.Grid() );

    if( progress && grid.Rank() == 0 )
        Output("Pulling ATL and ABR");
    ATL = ATLSub;
    ABR = ABRSub;

    // This section is a hack since no inter-grid redistributions exist for 
    // [VR,* ] distributions yet
    if( progress && grid.Rank() == 0 )
        Output("Pulling wT and wB");
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
        Output("Pulling ZT and ZB");
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
( ElementalMatrix<F>& APre,
  ElementalMatrix<Complex<Base<F>>>& wPre, 
  ElementalMatrix<F>& QPre,
  bool fullTriangle=true, 
  const SDCCtrl<Base<F>> ctrl=SDCCtrl<Base<F>>() )
{
    DEBUG_ONLY(
      CSE cse("schur::SDC");
      AssertSameGrids( APre, wPre, QPre );
    )
    typedef Base<F> Real;
    typedef Complex<Real> C;

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    DistMatrixWriteProxy<C,C,VR,STAR> wProx( wPre );
    DistMatrixWriteProxy<F,F,MC,MR> QProx( QPre );
    auto& A = AProx.Get();
    auto& w = wProx.Get();
    auto& Q = QProx.Get();

    const Grid& g = A.Grid();
    const Int n = A.Height();
    w.Resize( n, 1 );
    Q.Resize( n, n );
    if( A.Grid().Size() == 1 )
    {
        if( ctrl.progress && g.Rank() == 0 )
            Output("One process: using QR algorithm");
        Schur( A.Matrix(), w.Matrix(), Q.Matrix(), fullTriangle );
        return;
    }
    if( n <= ctrl.cutoff )
    {
        if( ctrl.progress && g.Rank() == 0 )
            Output(n," <= ",ctrl.cutoff,": using QR algorithm");
        DistMatrix<F,CIRC,CIRC> A_CIRC_CIRC( A ), Q_CIRC_CIRC( n, n, g );
        DistMatrix<Complex<Base<F>>,CIRC,CIRC> w_CIRC_CIRC( n, 1, g );
        if( A_CIRC_CIRC.CrossRank() == A_CIRC_CIRC.Root() )
            Schur
            ( A_CIRC_CIRC.Matrix(), w_CIRC_CIRC.Matrix(), Q_CIRC_CIRC.Matrix(),
              fullTriangle );
        A = A_CIRC_CIRC;
        w = w_CIRC_CIRC;
        Q = Q_CIRC_CIRC;
        return;
    }

    // Perform this level's split
    if( ctrl.progress && g.Rank() == 0 )
        Output("Splitting ",n," x ",n," matrix");
    const auto part = SpectralDivide( A, Q, ctrl );
    auto ind1 = IR(0,part.index);
    auto ind2 = IR(part.index,n);

    auto ATL = A( ind1, ind1 );
    auto ATR = A( ind1, ind2 );
    auto ABL = A( ind2, ind1 );
    auto ABR = A( ind2, ind2 );

    auto wT = w( ind1, ALL );
    auto wB = w( ind2, ALL );

    auto QL = Q( ALL, ind1 );
    auto QR = Q( ALL, ind2 );

    Zero( ABL );

    // Recurse on the two subproblems
    DistMatrix<F> ATLSub, ABRSub, ZTSub, ZBSub;
    DistMatrix<Complex<Base<F>>,VR,STAR> wTSub, wBSub;
    if( ctrl.progress && g.Rank() == 0 )
        Output("Pushing subproblems");
    PushSubproblems
    ( ATL, ABR, ATLSub, ABRSub, wT, wB, wTSub, wBSub, ZTSub, ZBSub, 
      ctrl.progress );
    if( ATLSub.Participating() )
        SDC( ATLSub, wTSub, ZTSub, fullTriangle, ctrl );
    if( ABRSub.Participating() )
        SDC( ABRSub, wBSub, ZBSub, fullTriangle, ctrl );
    
    // Ensure that the results are back on this level's grid
    if( ctrl.progress && g.Rank() == 0 )
        Output("Pulling subproblems");
    DistMatrix<F> ZT(g), ZB(g);
    PullSubproblems
    ( ATL, ABR, ATLSub, ABRSub, wT, wB, wTSub, wBSub, ZT, ZB, ZTSub, ZBSub,
      ctrl.progress );

    // Update the Schur vectors
    if( ctrl.progress && g.Rank() == 0 )
        Output("Updating Schur vectors");
    auto G( QL );
    Gemm( NORMAL, NORMAL, F(1), G, ZT, QL );
    G = QR;
    Gemm( NORMAL, NORMAL, F(1), G, ZB, QR );

    if( fullTriangle )
    {
        if( ctrl.progress && g.Rank() == 0 )
            Output("Updating top-right quadrant");
        // Update the top-right quadrant
        Gemm( ADJOINT, NORMAL, F(1), ZT, ATR, G );
        Gemm( NORMAL, NORMAL, F(1), G, ZB, ATR ); 
    }
}

} // namespace schur
} // namespace El

#endif // ifndef EL_SCHUR_SDC_HPP
