/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BIDIAG_SVD_DC_HPP
#define EL_BIDIAG_SVD_DC_HPP

namespace El {
namespace bidiag_svd {

// The following is analogous to LAPACK's {s,d}lasd{1,2,3} [CITATION] but does
// not accept initial sorting permutations for s0 and s1, nor does it enforce 
// any ordering on the resulting singular values. Several bugs in said LAPACK
// routines were found and reported to
// https://github.com/Reference-LAPACK/lapack/issues/34.
template<typename Real>
DCInfo
Merge
( Real alpha,
  // The right entry in the removed middle row of the bidiagonal matrix
  Real beta,
  // The non-deflated m0 (unsorted) singular values from B0.
  const Matrix<Real>& s0,
  // The non-deflated m1 (unsorted) singular values from B1.
  const Matrix<Real>& s1,
  // On entry, a packing of the left singular vectors from the two subproblems,
  //
  //   U = | U0, 0, 0  |,
  //       | 0,  1, 0  |
  //       | 0,  0, U1 |
  //
  // where U0 is m0 x (m0+1) and U1 is either m1 x m1 or m1 x (m1+1).
  //
  // On exit, the left singular vectors of the merged bidiagonal matrix.
  Matrix<Real>& U,
  // On exit, the (unsorted) singular values of the merged bidiagonal matrix
  Matrix<Real>& d,
  // On entry, a packing of the right singular vectors from the two subproblems,
  //
  //   V = | V0, 0  |,
  //       | 0,  V1 |
  //
  // where V0 is (m0+1) x (m0+1), with its last column lying in the null space
  // of B0, and V1 is either m1 x m1 or (m1+1) x (m1+1), where, in the latter
  // case, its last column must lie in the null space of B1.
  //
  // On exit, the right singular vectors of the merged bidiagonal matrix.
  Matrix<Real>& V,
  const DCCtrl<Real>& ctrl )
{
    DEBUG_CSE
    const Int m = U.Height();
    const Int n = V.Height();
    const Int m0 = s0.Height();
    const Int m1 = s1.Height();
    const Int n0 = m0 + 1;
    const Int n1 = n - n0;
    const bool square = ( m1 == n1 );
    DEBUG_ONLY(
      if( !square && n1 != m1+1 )
          LogicError("B1 has to be square or one column wider than tall");
    )
    DCInfo info;
    auto& secularInfo = info.secularInfo;
    auto& deflationInfo = info.deflationInfo;
    if( ctrl.progress )
        Output("m=",m,", n=",n,", m0=",m0,", n0=",n0,", m1=",m1,", n1=",n1);

    // U = | U0 0 0  |
    //     | 0  1 0  |
    //     | 0  0 U1 |
    auto U0 = U( IR(0,m0), IR(0,m0) );
    auto U1 = U( IR(n0,END), IR(n0,END) );

    // V = | V0 0  | 
    //     |  0 V1 |
    auto V0 = V( IR(0,n0), IR(0,n0) );
    auto V1 = V( IR(n0,END), IR(n0,END) );

    // Before permutation, 
    //
    //   r = [ alpha*V0(m0,:), beta*V1(0,:) ],
    //
    // but we reorder indices 0 and m0 to put r in the first position. We also
    // form d = [0; s0; s1]. Thus, d and r would provide a representation of
    //
    //    | r(0), r(1), ..., r(m-1) |,
    //    |       d(1),         .   |
    //    |              .      .   |
    //    |                  d(m-1) |
    //
    // or
    //
    //    | r(0), r(1), ..., r(m-1), rhoExtra |,
    //    |       d(1),         .        0    |
    //    |              .      .        0    |
    //    |                  d(m-1),     0    |
    //
    // depending upon whether B is m x m or m x (m+1). In the latter case, we
    // will rotate rhoExtra into r(0).
    //

    // Form d = [0; s0; s1].
    // This effectively cyclically shifts [0,m0] |-> [1,m0+1] mod (m0+1).
    d.Resize( m, 1 );
    d(0) = 0;
    for( Int j=0; j<m0; ++j )
    {
        d(j+1) = s0(j);
    }
    for( Int j=0; j<m1; ++j )
    {
        d(j+n0) = s1(j);
    }

    // Compute the scale of the problem and rescale {d,alpha,beta}. We will
    // rescale the singular values at the end of this routine.
    Real scale = Max( Abs(alpha), Abs(beta) );
    scale = Max( scale, MaxNorm(s0) );
    scale = Max( scale, MaxNorm(s1) );
    SafeScale( Real(1), scale, d );
    SafeScale( Real(1), scale, alpha );
    SafeScale( Real(1), scale, beta );

    // Now that the problem is rescaled, our deflation tolerance simplifies to
    //
    //   tol = deflationFudge eps max( || d ||_max, |alpha|, |beta| )
    //       = deflationFudge eps.
    //
    // Cf. LAPACK's {s,d}lasd2 [CITATION] for this tolerance.
    const Real eps = limits::Epsilon<Real>();
    const Real deflationTol = ctrl.deflationFudge*eps;

    Matrix<Real> r(m,1);
    Matrix<Int> columnTypes(m,1);
    // Form the reordered left portion
    r(0) = alpha*V0(m0,m0);
    columnTypes(0) = DENSE_COLUMN;
    for( Int j=0; j<m0; ++j )
    {
        r(j+1) = alpha*V0(m0,j);
        columnTypes(j+1) = COLUMN_NONZERO_IN_FIRST_BLOCK;
    }
    for( Int j=0; j<m1; ++j )
    {
        r(j+n0) = beta*V1(0,j);
        columnTypes(j+n0) = COLUMN_NONZERO_IN_SECOND_BLOCK;
    }
    // Form r(m) if B has one more column than row and then compute the cosine
    // and sine defining the Givens rotation for rotating it into r(0). Then
    // ensure that |r(0)| >= deflationTol. The Givens rotation is such that
    // 
    //   | r(0), rhoExtra | | cExtra,  -sExtra | = | gamma, 0 |.
    //                      | sExtra,   cExtra |
    Real rhoExtra=0, cExtra=1, sExtra=0;
    if( n == m+1 )
    {
        rhoExtra = beta*V1(0,m1);
        const Real gamma = SafeNorm( r(0), rhoExtra );
        if( gamma <= deflationTol )
        {
            r(0) = Sgn(r(0),false)*deflationTol;
        }
        else
        {
            cExtra = r(0) / gamma;
            sExtra = rhoExtra / gamma;
            r(0) = gamma;
        }
        if( cExtra != Real(1) || sExtra != Real(0) )
        {
            // Since V was originally block-diagonal and the two relevant
            // columns have not changed, the m0'th column is zero in the first
            // (m0+1) entries and the m'th column is nonzero in the last
            // (m1+1) entries. Thus, the rotation takes the form
            //
            //    | V(0:m0,m0),      0      | | cExtra, -sExtra |.
            //    |      0,     V(m0+1:m,m) | | sExtra,  cExtra |
            //
            for( Int i=0; i<m0+1; ++i )
            {
                const Real nu = V(i,m0);
                V(i,m0) =  cExtra*nu;
                V(i,m)  = -sExtra*nu;
            }
            for( Int i=m0+1; i<n; ++i )
            {
                const Real nu = V(i,m);
                V(i,m0) = sExtra*nu;
                V(i,m)  = cExtra*nu;
            }
            // V(:,m) should now lie in the null space of the inner matrix.
        }
    }
    else
    {
        if( Abs(r(0)) < deflationTol )
            r(0) = Sgn(r(0),false)*deflationTol;
    }

    // We could avoid sorting d(0)=0, but it should not significantly effect
    // performance. We force the sort to be stable to force the first entry of
    // d to remain in place.
    Permutation combineSortPerm;
    bool stableSort = true;
    SortingPermutation( d, combineSortPerm, ASCENDING, stableSort );
    combineSortPerm.PermuteRows( d );
    combineSortPerm.PermuteRows( r );
    combineSortPerm.PermuteRows( columnTypes );

    auto combinedToOrig = [&]( const Int& combinedIndex )
      {
          const Int preCombined = combineSortPerm.Preimage( combinedIndex );
          if( preCombined <= m0 )
              // Undo the cyclic shift [0,m0] |-> [1,m0+1] mod (m0+1) which
              // pushed the removed row into the first position.
              return Mod( preCombined-1, m0+1 );
          else
              return preCombined;
      };

    Permutation deflationPerm;
    deflationPerm.MakeIdentity( m );
    deflationPerm.MakeArbitrary();
    // Since we do not yet know how many undeflated entries there will be, we
    // must use the no-deflation case as our storage upper bound.
    Matrix<Real> dUndeflated(m,1), rUndeflated(m,1);
    dUndeflated(0) = 0;
    rUndeflated(0) = r(0);

    // Deflate all (off-diagonal) update entries sufficiently close to zero
    Int numDeflated = 0;
    Int numUndeflated = 1; // We do not deflate the first index
    // We will keep track of the last column that we encountered that was not
    // initially deflatable (but that could be deflated later due to close
    // diagonal entries if another undeflatable column is not encountered
    // first).
    Int revivalCandidate = 0;
    for( Int j=1; j<m; ++j )
    {
        if( Abs(r(j)) <= deflationTol )
        {
            // We can deflate due to the r component being sufficiently small
            deflationPerm.SetImage( j, (m-1)-numDeflated );
            if( ctrl.progress )
                Output
                ("Deflating via p(",j,")=",(m-1)-numDeflated,
                 " because |r(",j,")|=|",r(j),"| <= ",deflationTol);
            columnTypes(j) = DEFLATED_COLUMN;
            ++numDeflated;
            ++deflationInfo.numDeflations;
            ++deflationInfo.numSmallUpdateDeflations;
        }
        else if( d(j) <= deflationTol )
        {
            // We can deflate due to d(0)=0 being close to d(j). We rotate r(j)
            // into r(0) (Cf. the discussion surrounding Eq. (4.3) of 
            // Gu/Eisenstat's TR [CITATION]).
            //
            // In particular, we want
            //
            //   | r(0), r(j) | | c -s | = | gamma, 0 |,
            //                  | s  c |
            //
            // where gamma = || r(0); r(j) ||_2. Putting 
            //
            //   c = r(0) / gamma,
            //   s = r(j) / gamma,
            //
            // implies
            //
            //   |  c, s | | r(0) | = | gamma |.
            //   | -s, c | | r(j) |   |   0   |
            //
            const Real f = r(0);
            const Real g = r(j);
            const Real gamma = SafeNorm( f, g );
            const Real c = f / gamma;
            const Real s = g / gamma;
            r(0) = rUndeflated(0) = gamma;
            r(j) = 0;

            // Apply | c -s | from the right to V.
            //       | s  c |
            //
            // We are mixing nonzero structures in the first column of U,
            // so we might as well always treat the first column as dense.
            //
            // TODO(poulson): Exploit the nonzero structure of V?
            const Int jOrig = combinedToOrig( j );
            blas::Rot( n, &V(0,m0), 1, &V(0,jOrig), 1, c, s );

            deflationPerm.SetImage( j, (m-1)-numDeflated );
            if( ctrl.progress )
                Output
                ("Deflating via p(",j,")=",(m-1)-numDeflated,
                 " because d(",j,")=",d(j)," <= ",deflationTol);

            columnTypes(j) = DEFLATED_COLUMN;

            ++numDeflated;
            ++deflationInfo.numDeflations;
            ++deflationInfo.numCloseDiagonalDeflations;
        }
        else
        {
            revivalCandidate = j;
            if( ctrl.progress )
                Output("Breaking initial deflation loop at j=",j);
            break;
        }
    }
    for( Int j=revivalCandidate+1; j<m; ++j )
    {
        if( Abs(r(j)) <= deflationTol )
        {
            deflationPerm.SetImage( j, (m-1)-numDeflated );
            if( ctrl.progress )
                Output
                ("Deflating via p(",j,")=",(m-1)-numDeflated,
                 " because |r(",j,")|=|",r(j),"| <= ",deflationTol);
            columnTypes(j) = DEFLATED_COLUMN;
            ++numDeflated;
            ++deflationInfo.numDeflations;
            ++deflationInfo.numSmallUpdateDeflations;
        }
        else if( d(j)-d(revivalCandidate) <= deflationTol )
        {
            // Deflate the previously undeflatable index by rotating
            // r(revivalCandidate) into r(j) (Cf. the discussion
            // surrounding Eq. (4.4) of Gu/Eisenstat's TR [CITATION]
            // but recall that we are operating on the transposed system).
            //
            // In particular, we want
            //
            //   | r(j), r(revivalCandidate) | | c -s | = | gamma, 0 |,
            //                                 | s  c |
            //
            // where gamma = || r(revivalCandidate); r(j) ||_2. Putting 
            //
            //   c = r(j)                / gamma,
            //   s = r(revivalCandidate) / gamma,
            //
            // implies
            //
            //   |  c,  s | |        r(j)         | = | gamma |,
            //   | -s,  c | | r(revivalCandidate) |   |   0   |
            //
            const Real f = r(j);
            const Real g = r(revivalCandidate);
            const Real gamma = SafeNorm( f, g );
            const Real c = f / gamma;
            const Real s = g / gamma;
            r(j) = gamma;
            r(revivalCandidate) = 0;

            // Apply | c -s | from the right to U and V
            //       | s  c |
            //
            // TODO(poulson): Exploit the nonzero structure of U and V?
            const Int revivalOrig = combinedToOrig( revivalCandidate );
            const Int jOrig = combinedToOrig( j );
            blas::Rot( m, &U(0,jOrig), 1, &U(0,revivalOrig), 1, c, s );
            blas::Rot( n, &V(0,jOrig), 1, &V(0,revivalOrig), 1, c, s );

            deflationPerm.SetImage( revivalCandidate, (m-1)-numDeflated );
            if( ctrl.progress )
                Output
                ("Deflating via p(",revivalCandidate,")=",
                 (m-1)-numDeflated," because d(",j,")=",d(j),
                 " - d(",revivalCandidate,")=",d(revivalCandidate)," <= ",
                 deflationTol);

            if( columnTypes(revivalCandidate) != columnTypes(j) )
            {
                // We mixed top and bottom columns so the result is dense.
                columnTypes(j) = DENSE_COLUMN;
            }
            columnTypes(revivalCandidate) = DEFLATED_COLUMN;

            revivalCandidate = j;
            ++numDeflated;
            ++deflationInfo.numDeflations;
            ++deflationInfo.numCloseDiagonalDeflations;
        }
        else
        {
            // We cannot yet deflate index j, so we must give up on the previous
            // revival candidate and then set revivalCandidate = j.
            dUndeflated(numUndeflated) = d(revivalCandidate);
            rUndeflated(numUndeflated) = r(revivalCandidate);
            deflationPerm.SetImage( revivalCandidate, numUndeflated );
            if( ctrl.progress )
                Output
                ("Could not deflate with j=",j," and revivalCandidate=",
                 revivalCandidate,", so p(",revivalCandidate,")=",
                 numUndeflated);
            ++numUndeflated;

            revivalCandidate = j;
        }

        if( j == m-1 )
        {
            // This is the last iteration, so give up on the revival candidate
            dUndeflated(numUndeflated) = d(revivalCandidate);
            rUndeflated(numUndeflated) = r(revivalCandidate);
            deflationPerm.SetImage( revivalCandidate, numUndeflated );
            if( ctrl.progress )
                Output
                ("Final iteration, so p(",revivalCandidate,")=",numUndeflated);
            ++numUndeflated;
        }
    }
    // Now shrink dUndeflated and rUndeflated down to their proper size
    dUndeflated.Resize( numUndeflated, 1 );
    rUndeflated.Resize( numUndeflated, 1 );

    // Count the number of columns of U with each nonzero pattern
    std::vector<Int> packingCounts( NUM_DC_COMBINED_COLUMN_TYPES, 0 );
    for( Int j=0; j<m; ++j )
        ++packingCounts[columnTypes(j)];
    DEBUG_ONLY(
      if( packingCounts[DEFLATED_COLUMN] != numDeflated )
          LogicError
          ("Inconsistency between packingCounts[DEFLATED_COLUMN]=",
           packingCounts[DEFLATED_COLUMN],
           " and numDeflated=",numDeflated);
    )

    // Compute offsets for packing them
    std::vector<Int> packingOffsets( NUM_DC_COMBINED_COLUMN_TYPES, 0 );
    Int totalPacked = 0;
    for( Int columnType=0; columnType<NUM_DC_COMBINED_COLUMN_TYPES;
         ++columnType )
    {
        packingOffsets[columnType] = totalPacked;
        totalPacked += packingCounts[columnType];
        if( ctrl.progress )
            Output("packingCounts[",columnType,"]=",packingCounts[columnType]);
    }

    // Set up the index ranges of the three packed column subsets
    const Range<Int> packingInd0(packingOffsets[0],packingOffsets[1]),
                     packingInd1(packingOffsets[1],packingOffsets[2]),
                     packingInd2(packingOffsets[2],packingOffsets[3]);

    Matrix<Real> dPacked;
    Matrix<Real> UPacked, VPacked;
    dPacked.Resize( m, 1 );
    UPacked.Resize( m, m );
    VPacked.Resize( n, m );
    Permutation packingPerm;
    packingPerm.MakeIdentity( m );
    for( Int j=0; j<m; ++j )
    {
        // Recall that columnTypes maps the indices in the *undeflated* ordering
        // to their column type, whereas packingPerm maps the *deflated*
        // ordering into the packed ordering.
        //
        // It is important to notice that packingPerm will map entries from
        // [0,numUndeflated) back into [0,numUndeflated).
        const Int packingSource = deflationPerm.Image(j);
        const Int packingDest = packingOffsets[columnTypes(j)]++;
        packingPerm.SetImage( packingSource, packingDest );

        const Int jOrig = combinedToOrig( j );

        dPacked(packingDest) = d(j);
        // TODO(poulson): Exploit the nonzero structure of U and V?
        blas::Copy( m, &U(0,jOrig), 1, &UPacked(0,packingDest), 1 );
        blas::Copy( n, &V(0,jOrig), 1, &VPacked(0,packingDest), 1 );
    }

    // Put the deflated columns in their final destination and shrink UPacked
    // and VPacked back down to their final sizes
    //
    // TODO(poulson): Exploit the nonzero structure of U and V?
    if( numDeflated > 0 )
    {
        blas::Copy
        ( numDeflated, &dPacked(numUndeflated), 1, &d(numUndeflated), 1 );
        lapack::Copy
        ( 'A', m, numDeflated,
          &UPacked(0,numUndeflated), UPacked.LDim(),
          &U(0,numUndeflated), U.LDim() );
        lapack::Copy
        ( 'A', n, numDeflated,
          &VPacked(0,numUndeflated), VPacked.LDim(),
          &V(0,numUndeflated), V.LDim() );
    }
    UPacked.Resize( m, numUndeflated );
    VPacked.Resize( n, numUndeflated );

    // Now compute the updated singular vectors using UPacked/VPacked
    // ==============================================================
    auto undeflatedInd = IR(0,numUndeflated);
    const Real rUndeflatedNorm = FrobeniusNorm( rUndeflated );
    rUndeflated *= Real(1) / rUndeflatedNorm;
    const Real rho = rUndeflatedNorm*rUndeflatedNorm;

    if( ctrl.progress )
        Output("Computing corrected update vector");
    Matrix<Real> rCorrected;
    Ones( rCorrected, numUndeflated, 1 );
    auto vScratch = V(undeflatedInd,IR(0));
    for( Int j=0; j<numUndeflated; ++j )
    {
        auto u = U(undeflatedInd,IR(j));
        auto valueInfo =
          SecularSingularValue
          ( j, dUndeflated, rho, rUndeflated, d(j), u, vScratch,
            ctrl.secularCtrl );

        secularInfo.numIterations += valueInfo.numIterations;
        secularInfo.numAlternations += valueInfo.numAlternations;
        secularInfo.numCubicIterations += valueInfo.numCubicIterations;
        secularInfo.numCubicFailures += valueInfo.numCubicFailures;

        // u currently holds dUndeflated-d(j) and vScratch currently holds
        // dUndeflated+d(j). Overwrite u with their element-wise product since
        // that is all we require from here on out.
        for( Int k=0; k<numUndeflated; ++k )
            u(k) *= vScratch(k);

        rCorrected(j) *= u(j);
        for( Int k=0; k<numUndeflated; ++k )
        {
            if( k == j )
                continue;
            rCorrected(k) *= u(k) /
              ((dUndeflated(j)+dUndeflated(k))*(dUndeflated(j)-dUndeflated(k)));
        }
    }
    for( Int j=0; j<numUndeflated; ++j )
        rCorrected(j) = Sgn(rUndeflated(j),false) * Sqrt(Abs(rCorrected(j)));

    // Compute the unnormalized left and right singular vectors via Eqs. (3.4)
    // and (3.3), respectively, from Gu/Eisenstat [CITATION].
    for( Int j=0; j<numUndeflated; ++j )
    {
        auto u = U(undeflatedInd,IR(j));
        auto v = V(undeflatedInd,IR(j));
        {
            const Real deltaSqMinusShiftSq = u(0);
            u(0) = -1;
            v(0) = rCorrected(0) / deltaSqMinusShiftSq;
        }
        for( Int i=1; i<numUndeflated; ++i )
        {
            const Real deltaSqMinusShiftSq = u(i);
            v(i) = rCorrected(i) / deltaSqMinusShiftSq;
            u(i) = dUndeflated(i) * v(i);
        }
    }

    // Form the normalized left singular vectors with the rows permuted by
    // the inverse of the packing permutation in Q. This allows the product
    // of UPacked with Q to be equal to the unpacked U times the left singular
    // vectors from the secular equation.
    Matrix<Real> Q;
    Zeros( Q, numUndeflated, numUndeflated );
    for( Int j=0; j<numUndeflated; ++j )
    {
        auto u = U(undeflatedInd,IR(j));
        auto q = Q(undeflatedInd,IR(j));
        const Real uFrob = FrobeniusNorm( u );
        for( Int i=0; i<numUndeflated; ++i )
            q(i) = u(packingPerm.Preimage(i)) / uFrob;
    }
    // Overwrite the first 'numUndeflated' columns of U with the updated left
    // singular vectors by exploiting the partitioning of Z = UPacked as,
    //
    //   Z = | Z_{0,0} |    0    | Z_{0,2} |, 
    //       |---------|---------|---------|
    //       |    0    |    0    | z_{1,2} |
    //       |---------|---------|---------|
    //       |    0    | Z_{2,1} | Z_{2,2} |
    //
    // where the first, second, and third block rows are respectively of heights
    // m0, 1, and m1, and the first, second, and third block columns 
    // respectively have widths packingCounts[0], packingCounts[1], and
    // packingCounts[2].
    //
    // Conformally partitioning Q, we have
    //
    //  Z Q = Z_{:,2} Q2 + | Z_{0,0} Q_0 |.
    //                     |-------------|
    //                     |      0      |
    //                     |-------------|
    //                     | Z_{2,1} Q_1 |
    //
    if( ctrl.progress )
        Output("Overwriting left singular vectors");
    auto UUndeflated = U( ALL, undeflatedInd );
    if( ctrl.exploitStructure )
    {
        auto Z2 = UPacked( ALL, packingInd2 );
        auto Q2 = Q( packingInd2, ALL );
        Gemm( NORMAL, NORMAL, Real(1), Z2, Q2, UUndeflated );

        // Finish updating the first block row
        auto U0 = UUndeflated( IR(0,m0), ALL );
        auto Z00 = UPacked( IR(0,m0), packingInd0 );
        auto Q0 = Q( packingInd0, ALL );
        Gemm( NORMAL, NORMAL, Real(1), Z00, Q0, Real(1), U0 );

        // Finish updating the last block row
        auto U2 = UUndeflated( IR(n0,m), ALL );
        auto Z21 = UPacked( IR(n0,m), packingInd1 );
        auto Q1 = Q( packingInd1, ALL );
        Gemm( NORMAL, NORMAL, Real(1), Z21, Q1, Real(1), U2 );
    }
    else
    {
        Gemm( NORMAL, NORMAL, Real(1), UPacked, Q, UUndeflated );
    }

    // Form the normalized right singular vectors with the rows permuted by
    // the inverse of the packing permutation in Q. This allows the product
    // of VPacked with Q to be equal to the unpacked V times the right singular
    // vectors from the secular equation.
    for( Int j=0; j<numUndeflated; ++j )
    {
        auto v = V(undeflatedInd,IR(j));
        auto q = Q(undeflatedInd,IR(j));
        const Real vFrob = FrobeniusNorm( v );
        for( Int i=0; i<numUndeflated; ++i )
            q(i) = v(packingPerm.Preimage(i)) / vFrob;
    }
    // Overwrite the first 'numUndeflated' columns of V with the updated right
    // singular vectors by exploiting the partitioning of Z = VPacked as
    //
    //   Z = | Z_{0,0} |    0    | Z_{0,2} |,
    //       |---------|---------|---------|
    //       |    0    | Z_{1,1} | Z_{1,2} |
    //
    // where the first and second block rows have heights n0 and n1. The block
    // columns respectively have widths packingCounts[0], packingCounts[1], and
    // packingCounts[2].
    //
    // Conformally partitioning Q, we have
    //
    //   Z Q = Z_{:,2} Q2 + | Z_{0,0} Q_0 |,
    //                      |-------------|
    //                      | Z_{1,1} Q_1 |
    //
    if( ctrl.progress )
        Output("Overwriting right singular vectors");
    auto VUndeflated = V( ALL, undeflatedInd );
    if( ctrl.exploitStructure )
    {
        auto Z2 = VPacked( ALL, packingInd2 );
        auto Q2 = Q( packingInd2, ALL );
        Gemm( NORMAL, NORMAL, Real(1), Z2, Q2, VUndeflated );

        // Finish updating the first block row
        auto V0 = VUndeflated( IR(0,n0), ALL );
        auto Z00 = VPacked( IR(0,n0), packingInd0 );
        auto Q0 = Q( packingInd0, ALL );
        Gemm( NORMAL, NORMAL, Real(1), Z00, Q0, Real(1), V0 );

        // Finish updating the second block row
        auto V1 = VUndeflated( IR(n0,n), ALL );
        auto Z11 = VPacked( IR(n0,n), packingInd1 );
        auto Q1 = Q( packingInd1, ALL );
        Gemm( NORMAL, NORMAL, Real(1), Z11, Q1, Real(1), V1 );
    }
    else
    {
        Gemm( NORMAL, NORMAL, Real(1), VPacked, Q, VUndeflated );
    }

    // Rescale the singular values
    SafeScale( scale, Real(1), d );

    return info;
}

template<typename Real>
DCInfo
DivideAndConquer
( const Matrix<Real>& mainDiag,
  const Matrix<Real>& superDiag,
        Matrix<Real>& U,
        Matrix<Real>& s,
        Matrix<Real>& V,
  const DCCtrl<Real>& ctrl=DCCtrl<Real>() )
{
    DEBUG_CSE
    const Int m = mainDiag.Height();
    const Int n = superDiag.Height() + 1;
    DCInfo info;
    auto& secularInfo = info.secularInfo;
    auto& deflationInfo = info.deflationInfo;

    if( m <= ctrl.cutoff )
    {
        BidiagSVDCtrl<Real> bidiagSVDCtrl;
        bidiagSVDCtrl.approach = FULL_SVD; // We need any null space as well
        bidiagSVDCtrl.progress = ctrl.progress;
        bidiagSVDCtrl.useQR = true;
        BidiagSVD( UPPER, mainDiag, superDiag, U, s, V, bidiagSVDCtrl );
        return info;
    }

    // TODO(poulson): A more intelligent split point. Perhaps the row near 
    // m/2 with smallest norm should be chosen to encourage small entry
    // deflation (which tends to be vastly more common).
    const Int split = m/2;

    const Real& alpha = mainDiag(split);
    const Real& beta = superDiag(split);

    Identity( U, m, m );
    Zeros( V, n, n );

    auto mainDiag0 = mainDiag( IR(0,split), ALL );
    auto superDiag0 = superDiag( IR(0,split), ALL );
    auto U0 = U( IR(0,split), IR(0,split) );
    auto V0 = V( IR(0,split+1), IR(0,split+1) );
    Matrix<Real> s0;
    auto info0 = DivideAndConquer( mainDiag0, superDiag0, U0, s0, V0, ctrl );

    auto mainDiag1 = mainDiag( IR(split+1,END), ALL );
    auto superDiag1 = superDiag( IR(split+1,END), ALL );
    auto U1 = U( IR(split+1,END), IR(split+1,END) );
    auto V1 = V( IR(split+1,END), IR(split+1,END) );
    Matrix<Real> s1;
    auto info1 = DivideAndConquer( mainDiag1, superDiag1, U1, s1, V1, ctrl );

    info = Merge( alpha, beta, s0, s1, U, s, V, ctrl );

    secularInfo.numIterations += info0.secularInfo.numIterations;
    secularInfo.numIterations += info1.secularInfo.numIterations;
    secularInfo.numAlternations += info0.secularInfo.numAlternations;
    secularInfo.numAlternations += info1.secularInfo.numAlternations;
    secularInfo.numCubicIterations += info0.secularInfo.numCubicIterations;
    secularInfo.numCubicIterations += info1.secularInfo.numCubicIterations;
    secularInfo.numCubicFailures += info0.secularInfo.numCubicFailures;
    secularInfo.numCubicFailures += info1.secularInfo.numCubicFailures;
    deflationInfo.numDeflations += info0.deflationInfo.numDeflations;
    deflationInfo.numDeflations += info1.deflationInfo.numDeflations;
    deflationInfo.numSmallDiagonalDeflations +=
      info0.deflationInfo.numSmallDiagonalDeflations;
    deflationInfo.numSmallDiagonalDeflations +=
      info1.deflationInfo.numSmallDiagonalDeflations;
    deflationInfo.numCloseDiagonalDeflations +=
      info0.deflationInfo.numCloseDiagonalDeflations;
    deflationInfo.numCloseDiagonalDeflations +=
      info1.deflationInfo.numCloseDiagonalDeflations;
    deflationInfo.numSmallUpdateDeflations +=
      info0.deflationInfo.numSmallUpdateDeflations;
    deflationInfo.numSmallUpdateDeflations +=
      info1.deflationInfo.numSmallUpdateDeflations;

    return info;
}

} // namespace bidiag_svd
} // namespace El

#endif // ifndef EL_BIDIAG_SVD_DC_HPP
