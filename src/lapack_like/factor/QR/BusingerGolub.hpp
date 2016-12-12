/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_QR_BUSINGERGOLUB_HPP
#define EL_QR_BUSINGERGOLUB_HPP

namespace El {
namespace qr {

template<typename F>
Base<F> ColNorms( const Matrix<F>& A, vector<Base<F>>& norms )
{
    EL_DEBUG_CSE
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    Real maxNorm = 0;
    norms.resize( n );
    for( Int j=0; j<n; ++j )
    {
        norms[j] = blas::Nrm2( m, A.LockedBuffer(0,j), 1 );
        maxNorm = Max( maxNorm, norms[j] );
    }
    return maxNorm;
}

template<typename Real,class Compare=std::less<Real>>
ValueInt<Real> FindPivot
( const vector<Real>& norms,
        Int col,
        Compare compare=std::less<Real>() )
{
    EL_DEBUG_CSE
    if( norms.size()-col <= 0 )
    {
        ValueInt<Real> pivot;
        pivot.value = 0;
        pivot.index = -1;
        return pivot;
    }

    auto maxNorm = std::max_element( norms.begin()+col, norms.end(), compare );
    ValueInt<Real> pivot;
    pivot.value = *maxNorm;
    pivot.index = maxNorm - norms.begin();
    return pivot;
}

template<typename F>
void BusingerGolub
(       Matrix<F>& A,
        Matrix<F>& householderScalars,
        Matrix<Base<F>>& signature,
        Permutation& Omega,
  const QRCtrl<Base<F>> ctrl )
{
    EL_DEBUG_CSE
    typedef Base<F> Real;
    const Real zero(0), one(1);
    const F zeroF(0), oneF(1);
    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    const Int maxSteps = ( ctrl.boundRank ? Min(ctrl.maxRank,minDim) : minDim );
    householderScalars.Resize( maxSteps, 1 );
    signature.Resize( maxSteps, 1 );

    Matrix<F> z21;

    // Initialize two copies of the column norms, one will be consistently
    // updated, but the original copy will be kept to determine when the
    // updated quantities are no longer accurate.
    vector<Real> origNorms;
    const Real maxOrigNorm = ColNorms( A, origNorms );
    auto norms = origNorms;
    const Real updateTol = Sqrt(limits::Epsilon<Real>());

    Omega.MakeIdentity( n );
    Omega.ReserveSwaps( n );

    Int k=0;
    for( ; k<maxSteps; ++k )
    {
        const Range<Int> ind1( k ), ind2( k+1, END ), indB( k, END );

        auto alpha11 = A( ind1, ind1 );
        auto a12     = A( ind1, ind2 );
        auto a21     = A( ind2, ind1 );
        auto aB1     = A( indB, ind1 );
        auto AB2     = A( indB, ind2 );

        // Find the next column pivot
        ValueInt<Real> pivot;
        if( ctrl.smallestFirst )
        {
            // Adaptivity does not make sense for "rank-obscuring" factorization
            pivot = FindPivot( norms, k, std::greater<Real>() );
        }
        else
        {
            pivot = FindPivot( norms, k, std::less<Real>() );
            if( ctrl.adaptive && pivot.value <= ctrl.tol*maxOrigNorm )
                break;
        }
        Omega.Swap( k, pivot.index );

        // Perform the swap
        const Int jPiv = pivot.index;
        if( jPiv != k )
        {
            blas::Swap( m, &A(0,k), 1, &A(0,jPiv), 1 );
            norms[jPiv] = norms[k];
            origNorms[jPiv] = origNorms[k];
        }

        // Find tau and u such that
        //  / I - tau | 1 | | 1, u^H | \ | alpha11 | = | beta |
        //  \         | u |            / |     a21 | = |    0 |
        const F tau = LeftReflector( alpha11, a21 );
        householderScalars(k) = tau;

        // Temporarily set aB1 = | 1 |
        //                       | u |
        const F alpha = alpha11(0);
        alpha11(0) = 1;

        // AB2 := Hous(aB1,tau) AB2
        //      = (I - tau aB1 aB1^H) AB2
        //      = AB2 - tau aB1 (AB2^H aB1)^H
        Zeros( z21, AB2.Width(), 1 );
        Gemv( ADJOINT, oneF, AB2, aB1, zeroF, z21 );
        Ger( -tau, aB1, z21, AB2 );

        // Reset alpha11's value
        alpha11(0) = alpha;

        // Update the column norm estimates in the same manner as LAWN 176
        for( Int j=k+1; j<n; ++j )
        {
            if( norms[j] != zero )
            {
                Real gamma = Abs(A(k,j)) / norms[j];
                gamma = Max( zero, (one-gamma)*(one+gamma) );

                const Real ratio = norms[j] / origNorms[j];
                const Real phi = gamma*(ratio*ratio);
                if( phi <= updateTol || ctrl.alwaysRecomputeNorms )
                {
                    norms[j] = blas::Nrm2( m-(k+1), A.Buffer(k+1,j), 1 );
                    origNorms[j] = norms[j];
                }
                else
                    norms[j] *= Sqrt(gamma);
            }
        }
    }

    // Form d and rescale R
    auto R = A( IR(0,k), ALL );
    GetRealPartOfDiagonal(R,signature);
    auto sgn = [&]( const Real& delta ) { return delta >= zero ? one : -one; };
    EntrywiseMap( signature, MakeFunction(sgn) );
    DiagonalScaleTrapezoid( LEFT, UPPER, NORMAL, signature, R );

    // Ensure that t is the correct length
    householderScalars.Resize( k, 1 );
}

// TODO: Implement lambda op which is registered to MPI but can easily
//       be swapped out by Elemental?
//
//       In the mean time, we can simply make use of smallestFirst.
template<typename F>
ValueInt<Base<F>>
FindColPivot
( const DistMatrix<F>& A,
  const vector<Base<F>>& norms,
        Int col,
        bool smallestFirst=false )
{
    EL_DEBUG_CSE
    typedef Base<F> Real;
    const Int localColsBefore = A.LocalColOffset(col);
    ValueInt<Real> localPivot;
    if( smallestFirst )
        localPivot = FindPivot( norms, localColsBefore, std::greater<Real>() );
    else
        localPivot = FindPivot( norms, localColsBefore, std::less<Real>() );
    ValueInt<Real> pivot;
    pivot.value = localPivot.value;
    pivot.index = A.GlobalCol(localPivot.index);

    if( smallestFirst )
        return mpi::AllReduce
               ( pivot, mpi::MinLocOp<Real>(), A.Grid().RowComm() );
    else
        return mpi::AllReduce
               ( pivot, mpi::MaxLocOp<Real>(), A.Grid().RowComm() );
}

template<typename F>
Base<F> ColNorms( const DistMatrix<F>& A, vector<Base<F>>& norms )
{
    EL_DEBUG_CSE
    typedef Base<F> Real;
    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    const Matrix<F>& ALoc = A.LockedMatrix();
    mpi::Comm colComm = A.Grid().ColComm();
    mpi::Comm rowComm = A.Grid().RowComm();

    // Carefully perform the local portion of the computation
    vector<Real> localScales(localWidth,0),
                 localScaledSquares(localWidth,1);
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            UpdateScaledSquare
            ( ALoc(iLoc,jLoc), localScales[jLoc], localScaledSquares[jLoc] );

    // Find the maximum relative scales
    vector<Real> scales(localWidth);
    mpi::AllReduce
    ( localScales.data(), scales.data(), localWidth, mpi::MAX, colComm );

    // Equilibrate the local scaled sums to the maximum scale
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        if( scales[jLoc] != Real(0) )
        {
            const Real relScale = localScales[jLoc]/scales[jLoc];
            localScaledSquares[jLoc] *= relScale*relScale;
        }
    }

    // Now sum the local contributions (can ignore results where scale is 0)
    vector<Real> scaledSquares(localWidth);
    mpi::AllReduce
    ( localScaledSquares.data(), scaledSquares.data(), localWidth, colComm );

    // Finish the computation
    Real maxLocalNorm = 0;
    norms.resize( localWidth );
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        if( scales[jLoc] != Real(0) )
            norms[jLoc] = scales[jLoc]*Sqrt(scaledSquares[jLoc]);
        else
            norms[jLoc] = 0;
        maxLocalNorm = Max( maxLocalNorm, norms[jLoc] );
    }
    return mpi::AllReduce( maxLocalNorm, mpi::MAX, rowComm );
}

template<typename F>
void ReplaceColNorms
( const DistMatrix<F>& A,
        vector<Int>& inaccurateNorms,
        vector<Base<F>>& norms,
        vector<Base<F>>& origNorms )
{
    EL_DEBUG_CSE
    typedef Base<F> Real;
    const Int localHeight = A.LocalHeight();
    const Int numInaccurate = inaccurateNorms.size();
    const Matrix<F>& ALoc = A.LockedMatrix();
    mpi::Comm colComm = A.Grid().ColComm();

    // Carefully perform the local portion of the computation
    vector<Real> localScales(numInaccurate,0),
                 localScaledSquares(numInaccurate,1);
    for( Int s=0; s<numInaccurate; ++s )
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
            UpdateScaledSquare
            ( ALoc(iLoc,inaccurateNorms[s]),
              localScales[s], localScaledSquares[s] );

    // Find the maximum relative scales
    vector<Real> scales(numInaccurate);
    mpi::AllReduce
    ( localScales.data(), scales.data(), numInaccurate, mpi::MAX, colComm );

    // Equilibrate the local scaled sums to the maximum scale
    for( Int s=0; s<numInaccurate; ++s )
    {
        if( scales[s] != Real(0) )
        {
            const Real relScale = localScales[s]/scales[s];
            localScaledSquares[s] *= relScale*relScale;
        }
    }

    // Now sum the local contributions (can ignore results where scale is 0)
    vector<Real> scaledSquares(numInaccurate);
    mpi::AllReduce
    ( localScaledSquares.data(), scaledSquares.data(), numInaccurate, colComm );

    // Finish the computation
    for( Int s=0; s<numInaccurate; ++s )
    {
        const Int jLoc = inaccurateNorms[s];
        if( scales[s] != Real(0) )
            norms[jLoc] = scales[s]*Sqrt(scaledSquares[s]);
        else
            norms[jLoc] = 0;
        origNorms[jLoc] = norms[jLoc];
    }
}

template<typename F>
void BusingerGolub
( AbstractDistMatrix<F>& APre,
  AbstractDistMatrix<F>& householderScalars,
  AbstractDistMatrix<Base<F>>& signature,
  DistPermutation& Omega,
  const QRCtrl<Base<F>> ctrl )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(AssertSameGrids( APre, householderScalars, signature ))
    typedef Base<F> Real;
    const Real zero(0), one(1);

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.Get();

    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    const Int mLocal = A.LocalHeight();
    const Int maxSteps = ( ctrl.boundRank ? Min(ctrl.maxRank,minDim) : minDim );
    householderScalars.Resize( maxSteps, 1 );
    signature.Resize( maxSteps, 1 );

    // Initialize two copies of the column norms, one will be consistently
    // updated, but the original copy will be kept to determine when the
    // updated quantities are no longer accurate.
    vector<Real> origNorms( A.LocalWidth() );
    const Real maxOrigNorm = ColNorms( A, origNorms );
    auto norms = origNorms;
    const Real updateTol = Sqrt(limits::Epsilon<Real>());
    vector<Int> inaccurateNorms;

    Omega.MakeIdentity( n );
    Omega.ReserveSwaps( n );

    const Grid& g = A.Grid();
    DistMatrix<F> z21(g);
    DistMatrix<F,MC,STAR> aB1_MC(g);
    DistMatrix<F,MR,STAR> z21_MR(g);
    DistMatrix<F,STAR,MR> a12_MR(g);

    Int k=0;
    for( ; k<maxSteps; ++k )
    {
        const Range<Int> ind1( k ), ind2( k+1, END ), indB( k, END );

        auto alpha11 = A( ind1, ind1 );
        auto a12     = A( ind1, ind2 );
        auto a21     = A( ind2, ind1 );
        auto aB1     = A( indB, ind1 );
        auto AB2     = A( indB, ind2 );

        // Find the next column pivot
        ValueInt<Real> pivot;
        if( ctrl.smallestFirst )
        {
            // Adaptivity does not make sense for "rank-obscuring" factorization
            pivot = FindColPivot( A, norms, k, true );
        }
        else
        {
            pivot = FindColPivot( A, norms, k, false );
            if( ctrl.adaptive && pivot.value <= ctrl.tol*maxOrigNorm )
                break;
        }
        Omega.Swap( k, pivot.index );

        // Perform the swap
        const Int jPiv = pivot.index;
        const Int curOwner = A.ColOwner(k);
        const Int pivOwner = A.ColOwner(jPiv);
        const Int myCur = A.IsLocalCol(k);
        const Int myPiv = A.IsLocalCol(jPiv);
        if( jPiv != k )
        {
            if( myCur && myPiv )
            {
                const Int kLoc    = A.LocalCol(k);
                const Int jPivLoc = A.LocalCol(jPiv);
                blas::Swap
                ( mLocal, A.Buffer(0,kLoc), 1, A.Buffer(0,jPivLoc), 1 );
                norms[jPivLoc] = norms[kLoc];
                origNorms[jPivLoc] = origNorms[kLoc];
            }
            else if( myCur )
            {
                const Int kLoc = A.LocalCol(k);
                mpi::SendRecv
                ( A.Buffer(0,kLoc), mLocal, pivOwner, pivOwner, g.RowComm() );
                mpi::Send( norms[kLoc], pivOwner, g.RowComm() );
            }
            else if( myPiv )
            {
                const Int jPivLoc = A.LocalCol(jPiv);
                mpi::SendRecv
                ( A.Buffer(0,jPivLoc), mLocal,
                  curOwner, curOwner, g.RowComm() );
                norms[jPivLoc] = mpi::Recv<Real>( curOwner, g.RowComm() );
            }
        }

        // Find tau and u such that
        //  / I - tau | 1 | | 1, u^H | \ | alpha11 | = | beta |
        //  \         | u |            / |     a21 | = |    0 |
        const F tau = LeftReflector( alpha11, a21 );
        householderScalars.Set( k, 0, tau );

        // Temporarily set aB1 = | 1 |
        //                       | u |
        F alpha = 0;
        if( alpha11.IsLocal(0,0) )
        {
            alpha = alpha11.GetLocal(0,0);
            alpha11.SetLocal(0,0,1);
        }

        // AB2 := Hous(aB1,tau) AB2
        //      = (I - tau aB1 aB1^H) AB2
        //      = AB2 - tau aB1 (AB2^H aB1)^H
        aB1_MC.AlignWith( AB2 );
        aB1_MC = aB1;
        z21_MR.AlignWith( AB2 );
        Zeros( z21_MR, AB2.Width(), 1 );
        LocalGemv( ADJOINT, F(1), AB2, aB1_MC, F(0), z21_MR );
        El::AllReduce( z21_MR, AB2.ColComm() );
        Ger
        ( -tau, aB1_MC.LockedMatrix(), z21_MR.LockedMatrix(), AB2.Matrix() );

        // Reset alpha11's value
        if( alpha11.IsLocal(0,0) )
            alpha11.SetLocal(0,0,alpha);

        // Update the column norm estimates in the same manner as LAWN 176.
        // However, we do so in two steps in order to lower the communication
        // latency:
        //   1) Each process first computes which of its column norms are
        //      too inaccurate and need to be recomputed.
        //   2) Each process communicates within its process column in order
        //      to replace the inaccurate column norms.
        // Step 1: Perform all of the easy updates and mark inaccurate norms
        a12_MR = a12;
        inaccurateNorms.resize(0);
        const Int a12LocalWidth = a12_MR.LocalWidth();
        for( Int jLoc12=0; jLoc12<a12LocalWidth; ++jLoc12 )
        {
            const Int j = (k+1) + a12.GlobalCol(jLoc12);
            const Int jLoc = A.LocalCol(j);
            if( norms[jLoc] != zero )
            {
                const Real beta = Abs(a12_MR.GetLocal(0,jLoc12));
                Real gamma = beta / norms[jLoc];
                gamma = Max( zero, (one-gamma)*(one+gamma) );

                const Real ratio = norms[jLoc] / origNorms[jLoc];
                const Real phi = gamma*(ratio*ratio);
                if( phi <= updateTol || ctrl.alwaysRecomputeNorms )
                    inaccurateNorms.push_back( jLoc );
                else
                    norms[jLoc] *= Sqrt(gamma);
            }
        }
        // Step 2: Compute the replacement norms and also reset origNorms
        ReplaceColNorms( A(ind2,ALL), inaccurateNorms, norms, origNorms );
    }
    householderScalars.Resize( k, 1 );

    // Form d and rescale R
    auto R = A( IR(0,k), ALL );
    GetRealPartOfDiagonal(R,signature);
    auto sgn = [&]( const Real& delta ) { return delta >= zero ? one : -one; };
    EntrywiseMap( signature, MakeFunction(sgn) );
    DiagonalScaleTrapezoid( LEFT, UPPER, NORMAL, signature, R );
}

} // namespace qr
} // namespace El

#endif // ifndef EL_QR_BUSINGERGOLUB_HPP
