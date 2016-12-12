/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BIDIAG_SVD_DC_HPP
#define EL_BIDIAG_SVD_DC_HPP

// TODO(poulson): Move said routine into a utility function
#include "../Schur/SDC.hpp"
using El::schur::SplitGrid;

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
( Int m,
  Int n,
  Real alpha,
  // The right entry in the removed middle row of the bidiagonal matrix
  Real beta,
  // The non-deflated m0 (unsorted) singular values from B0.
  const Matrix<Real>& s0,
  // The non-deflated m1 (unsorted) singular values from B1.
  const Matrix<Real>& s1,
  // If ctrl.wantU is true, then, on entry, U contains a packing of the left
  // singular vectors from the two subproblems,
  //
  //   U = | U0, 0, 0  |,
  //       | 0,  1, 0  |
  //       | 0,  0, U1 |
  //
  // where U0 is m0 x (m0+1) and U1 is either m1 x m1 or m1 x (m1+1).
  //
  // If ctrl.wantU is true, then, on exit, the left singular vectors of the
  // merged bidiagonal matrix.
  Matrix<Real>& U,
  // On exit, the (unsorted) singular values of the merged bidiagonal matrix
  Matrix<Real>& d,
  // If ctrl.wantV is true, then, on entry, a packing of the right singular
  // vectors from the two subproblems,
  //
  //   V = | V0, 0  |,
  //       | 0,  V1 |
  //
  // where V0 is (m0+1) x (m0+1), with its last column lying in the null space
  // of B0, and V1 is either m1 x m1 or (m1+1) x (m1+1), where, in the latter
  // case, its last column must lie in the null space of B1.
  //
  // If ctrl.wantV is false, then, on entry, V is the same as above, but with
  // only the row that goes through the last row of V0 and the row that goes
  // through the first row of V1 kept.
  //
  // If ctrl.wantV is true, on exit, V will contain the right singular vectors
  // of the merged bidiagonal matrix. If ctrl.wantV is false, then only the two
  // rows of the result mentioned above will be output.
  Matrix<Real>& V,
  const BidiagSVDCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Int m0 = s0.Height();
    const Int m1 = s1.Height();
    const Int n0 = m0 + 1;
    const Int n1 = n - n0;
    const bool square = ( m1 == n1 );
    EL_DEBUG_ONLY(
      if( !square && n1 != m1+1 )
          LogicError("B1 has to be square or one column wider than tall");
    )
    const auto& dcCtrl = ctrl.dcCtrl;

    DCInfo info;
    auto& secularInfo = info.secularInfo;
    if( ctrl.progress )
        Output("m=",m,", n=",n,", m0=",m0,", n0=",n0,", m1=",m1,", n1=",n1);

    Matrix<Real> V0, V1;
    if( ctrl.wantV )
    {
        // V = | V0 0  |
        //     |  0 V1 |
        View( V0, V, IR(0,n0), IR(0,n0) );
        View( V1, V, IR(n0,END), IR(n0,END) );
    }
    else
    {
        View( V0, V, IR(0), IR(0,n0) );
        View( V1, V, IR(1), IR(n0,END) );
    }

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
    const Real deflationTol = dcCtrl.deflationFudge*eps;

    Matrix<Real> r(m,1);
    Matrix<Int> columnTypes(m,1);
    // Form the reordered left portion
    const Int lastRowOfV0 = ( ctrl.wantV ? m0 : 0 );
    r(0) = alpha*V0(lastRowOfV0,m0);
    columnTypes(0) = DENSE_COLUMN;
    for( Int j=0; j<m0; ++j )
    {
        r(j+1) = alpha*V0(lastRowOfV0,j);
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
    if( n == m+1 )
    {
        Real cExtra=1, sExtra=0;
        Real rhoExtra = beta*V1(0,m1);
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
            if( ctrl.wantV )
            {
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
            }
            else
            {
                const Real nu0 = V(0,m0);
                V(0,m0) = cExtra*nu0;
                V(0,m) = -sExtra*nu0;
                const Real nu1 = V(1,m);
                V(1,m0) = sExtra*nu1;
                V(1,m) = cExtra*nu1;
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
    // and sufficiently small diagonal entries.
    Int numDeflated = 0;
    Int numUndeflated = 1; // We do not deflate the first index
    // We will keep track of the last column that we encountered that was not
    // initially deflatable (but that could be deflated later due to close
    // diagonal entries if another undeflatable column is not encountered
    // first).
    Int revivalCandidate = m;
    for( Int j=1; j<m; ++j )
    {
        if( Abs(r(j)) <= deflationTol )
        {
            // We can deflate due to the r component being sufficiently small
            const Int deflationDest = (m-1) - numDeflated;
            deflationPerm.SetImage( j, deflationDest );
            if( ctrl.progress )
                Output
                ("Deflating via p(",j,")=",deflationDest,
                 " because |r(",j,")|=|",r(j),"| <= ",deflationTol);
            columnTypes(j) = DEFLATED_COLUMN;
            ++numDeflated;
            ++secularInfo.numDeflations;
            ++secularInfo.numSmallUpdateDeflations;
            continue;
        }
        if( d(j) <= deflationTol )
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
            const Int jOrig = combinedToOrig( j );
            if( ctrl.wantV )
            {
                // TODO(poulson): Exploit the nonzero structure of V?
                blas::Rot( n, &V(0,m0), 1, &V(0,jOrig), 1, c, s );
            }
            else
            {
                blas::Rot( 2, &V(0,m0), 1, &V(0,jOrig), 1, c, s );
            }

            const Int deflationDest = (m-1) - numDeflated;
            deflationPerm.SetImage( j, deflationDest );
            if( ctrl.progress )
                Output
                ("Deflating via p(",j,")=",deflationDest,
                 " because d(",j,")=",d(j)," <= ",deflationTol);

            columnTypes(j) = DEFLATED_COLUMN;

            ++numDeflated;
            ++secularInfo.numDeflations;
            ++secularInfo.numCloseDiagonalDeflations;
            continue;
        }

        revivalCandidate = j;
        if( ctrl.progress )
            Output("Breaking initial deflation loop at j=",j);
        break;
    }
    // If we already fully deflated, then the following loop should be trivial
    const Int deflationRestart = revivalCandidate+1;
    for( Int j=deflationRestart; j<m; ++j )
    {
        if( Abs(r(j)) <= deflationTol )
        {
            const Int deflationDest = (m-1) - numDeflated;
            deflationPerm.SetImage( j, deflationDest );
            if( ctrl.progress )
                Output
                ("Deflating via p(",j,")=",deflationDest,
                 " because |r(",j,")|=|",r(j),"| <= ",deflationTol);
            columnTypes(j) = DEFLATED_COLUMN;
            ++numDeflated;
            ++secularInfo.numDeflations;
            ++secularInfo.numSmallUpdateDeflations;
            continue;
        }
        if( d(j)-d(revivalCandidate) <= deflationTol )
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
            //   |  c,  s | |        r(j)         | = | gamma |.
            //   | -s,  c | | r(revivalCandidate) |   |   0   |
            //
            // TODO(poulson): Switch to the more stringent deflation criterion
            // used for the symmetric tridiagonal D&C.
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
            const Int revivalOrig = combinedToOrig( revivalCandidate );
            const Int jOrig = combinedToOrig( j );
            if( ctrl.wantU )
            {
                // TODO(poulson): Exploit the nonzero structure of U?
                blas::Rot( m, &U(0,jOrig), 1, &U(0,revivalOrig), 1, c, s );
            }
            if( ctrl.wantV )
            {
                // TODO(poulson): Exploit the nonzero structure of V?
                blas::Rot( n, &V(0,jOrig), 1, &V(0,revivalOrig), 1, c, s );
            }
            else
            {
                blas::Rot( 2, &V(0,jOrig), 1, &V(0,revivalOrig), 1, c, s );
            }

            const Int deflationDest = (m-1) - numDeflated;
            deflationPerm.SetImage( revivalCandidate, deflationDest );
            if( ctrl.progress )
                Output
                ("Deflating via p(",revivalCandidate,")=",
                 deflationDest," because d(",j,")=",d(j),
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
            ++secularInfo.numDeflations;
            ++secularInfo.numCloseDiagonalDeflations;
            continue;
        }

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
    if( revivalCandidate < m )
    {
        // Give up on the revival candidate
        dUndeflated(numUndeflated) = d(revivalCandidate);
        rUndeflated(numUndeflated) = r(revivalCandidate);
        deflationPerm.SetImage( revivalCandidate, numUndeflated );
        if( ctrl.progress )
            Output
            ("Final iteration, so p(",revivalCandidate,")=",numUndeflated);
        ++numUndeflated;
    }
    // Now shrink dUndeflated and rUndeflated down to their proper size
    dUndeflated.Resize( numUndeflated, 1 );
    rUndeflated.Resize( numUndeflated, 1 );

    // Count the number of columns of U with each nonzero pattern
    std::vector<Int> packingCounts( NUM_DC_COMBINED_COLUMN_TYPES, 0 );
    for( Int j=0; j<m; ++j )
        ++packingCounts[columnTypes(j)];
    EL_DEBUG_ONLY(
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
    if( ctrl.wantU )
        UPacked.Resize( m, m );
    if( ctrl.wantV )
        VPacked.Resize( n, m );
    else
        VPacked.Resize( 2, m );
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
        if( ctrl.wantU )
        {
            // TODO(poulson): Exploit the nonzero structure of U?
            blas::Copy( m, &U(0,jOrig), 1, &UPacked(0,packingDest), 1 );
        }
        if( ctrl.wantV )
        {
            // TODO(poulson): Exploit the nonzero structure of V?
            blas::Copy( n, &V(0,jOrig), 1, &VPacked(0,packingDest), 1 );
        }
        else
        {
            blas::Copy( 2, &V(0,jOrig), 1, &VPacked(0,packingDest), 1 );
        }
    }

    // Put the deflated columns in their final destination and shrink UPacked
    // and VPacked back down to their final sizes
    //
    // TODO(poulson): Exploit the nonzero structure of U and V?
    if( numDeflated > 0 )
    {
        blas::Copy
        ( numDeflated, &dPacked(numUndeflated), 1, &d(numUndeflated), 1 );
        if( ctrl.wantU )
        {
            lapack::Copy
            ( 'A', m, numDeflated,
              &UPacked(0,numUndeflated), UPacked.LDim(),
              &U(0,numUndeflated), U.LDim() );
        }
        if( ctrl.wantV )
        {
            lapack::Copy
            ( 'A', n, numDeflated,
              &VPacked(0,numUndeflated), VPacked.LDim(),
              &V(0,numUndeflated), V.LDim() );
        }
        else
        {
            lapack::Copy
            ( 'A', 2, numDeflated,
              &VPacked(0,numUndeflated), VPacked.LDim(),
              &V(0,numUndeflated), V.LDim() );
        }
    }
    if( ctrl.wantU )
        UPacked.Resize( m, numUndeflated );
    if( ctrl.wantV )
        VPacked.Resize( n, numUndeflated );
    else
        VPacked.Resize( 2, numUndeflated );

    // Now compute the updated singular vectors using UPacked/VPacked
    // ==============================================================
    auto undeflatedInd = IR(0,numUndeflated);
    const Real rUndeflatedNorm = FrobeniusNorm( rUndeflated );
    rUndeflated *= Real(1) / rUndeflatedNorm;
    const Real rho = rUndeflatedNorm*rUndeflatedNorm;

    if( ctrl.progress )
        Output("Solving secular equation and correcting update vector");
    Matrix<Real> rCorrected;
    Ones( rCorrected, numUndeflated, 1 );

    // Ensure that there is sufficient space for storing the needed singular
    // vectors from the undeflated secular equation. Notice that we *always*
    // need to compute the right singular vectors of the undeflated secular
    // equation.
    Matrix<Real> USecular, VSecular;
    if( ctrl.wantU )
        View( USecular, U, undeflatedInd, undeflatedInd );
    if( ctrl.wantV )
        View( VSecular, V, undeflatedInd, undeflatedInd );
    else
        VSecular.Resize( numUndeflated, numUndeflated );

    // For temporarily storing dUndeflated + d(j)
    Matrix<Real> plusShift;
    if( ctrl.wantU )
        View( plusShift, USecular, ALL, IR(0) );
    else
        plusShift.Resize( numUndeflated, 1 );

    for( Int j=0; j<numUndeflated; ++j )
    {
        auto minusShift = VSecular( ALL, IR(j) );

        auto valueInfo =
          SecularSingularValue
          ( j, dUndeflated, rho, rUndeflated, d(j), minusShift, plusShift,
            dcCtrl.secularCtrl );
        if( ctrl.progress )
            Output("Secular singular value ",j," is ",d(j));

        secularInfo.numIterations += valueInfo.numIterations;
        secularInfo.numAlternations += valueInfo.numAlternations;
        secularInfo.numCubicIterations += valueInfo.numCubicIterations;
        secularInfo.numCubicFailures += valueInfo.numCubicFailures;

        // minusShift currently holds dUndeflated-d(j) and plusShift
        // holds dUndeflated+d(j). Overwrite minusShift with their
        // element-wise product since that is all we require from here on
        // out.
        for( Int k=0; k<numUndeflated; ++k )
            minusShift(k) *= plusShift(k);

        rCorrected(j) *= minusShift(j);
        for( Int k=0; k<numUndeflated; ++k )
        {
            if( k == j )
                continue;
            rCorrected(k) *= minusShift(k) /
              ((dUndeflated(j)+dUndeflated(k))*
               (dUndeflated(j)-dUndeflated(k)));
        }
    }
    for( Int j=0; j<numUndeflated; ++j )
        rCorrected(j) = Sgn(rUndeflated(j),false) * Sqrt(Abs(rCorrected(j)));

    // Compute the unnormalized left and right singular vectors via Eqs. (3.4)
    // and (3.3), respectively, from Gu/Eisenstat [CITATION].
    if( ctrl.progress )
        Output("Computing unnormalized singular vectors");
    if( ctrl.wantU )
    {
        for( Int j=0; j<numUndeflated; ++j )
        {
            auto u = USecular(ALL,IR(j));
            auto v = VSecular(ALL,IR(j));
            {
                const Real deltaSqMinusShiftSq = v(0);
                v(0) = rCorrected(0) / deltaSqMinusShiftSq;
                u(0) = -1;
            }
            for( Int i=1; i<numUndeflated; ++i )
            {
                const Real deltaSqMinusShiftSq = v(i);
                v(i) = rCorrected(i) / deltaSqMinusShiftSq;
                u(i) = dUndeflated(i) * v(i);
            }
        }
    }
    else
    {
        for( Int j=0; j<numUndeflated; ++j )
        {
            auto v = VSecular(ALL,IR(j));
            {
                const Real deltaSqMinusShiftSq = v(0);
                v(0) = rCorrected(0) / deltaSqMinusShiftSq;
            }
            for( Int i=1; i<numUndeflated; ++i )
            {
                const Real deltaSqMinusShiftSq = v(i);
                v(i) = rCorrected(i) / deltaSqMinusShiftSq;
            }
        }
    }

    // Form the normalized left singular vectors with the rows permuted by
    // the inverse of the packing permutation in Q. This allows the product
    // of UPacked with Q to be equal to the unpacked U times the left singular
    // vectors from the secular equation.
    if( ctrl.progress )
        Output("Forming undeflated left singular vectors");
    Matrix<Real> Q;
    if( ctrl.wantU )
    {
        Zeros( Q, numUndeflated, numUndeflated );
        for( Int j=0; j<numUndeflated; ++j )
        {
            auto u = USecular(ALL,IR(j));
            auto q = Q(ALL,IR(j));
            const Real uFrob = FrobeniusNorm( u );
            for( Int i=0; i<numUndeflated; ++i )
                q(i) = u(packingPerm.Preimage(i)) / uFrob;
        }
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
    if( ctrl.wantU )
    {
        if( ctrl.progress )
            Output("Overwriting left singular vectors");
        auto UUndeflated = U( ALL, undeflatedInd );
        if( dcCtrl.exploitStructure )
        {
            auto Z2 = UPacked( ALL, packingInd2 );
            auto Q2 = Q( packingInd2, ALL );
            Gemm( NORMAL, NORMAL, Real(1), Z2, Q2, UUndeflated );

            // Finish updating the first block row
            auto U0Undeflated = UUndeflated( IR(0,m0), ALL );
            auto Z00 = UPacked( IR(0,m0), packingInd0 );
            auto Q0 = Q( packingInd0, ALL );
            Gemm( NORMAL, NORMAL, Real(1), Z00, Q0, Real(1), U0Undeflated );

            // Finish updating the last block row
            auto U2Undeflated = UUndeflated( IR(n0,m), ALL );
            auto Z21 = UPacked( IR(n0,m), packingInd1 );
            auto Q1 = Q( packingInd1, ALL );
            Gemm( NORMAL, NORMAL, Real(1), Z21, Q1, Real(1), U2Undeflated );
        }
        else
        {
            Gemm( NORMAL, NORMAL, Real(1), UPacked, Q, UUndeflated );
        }
    }

    // Form the normalized right singular vectors with the rows permuted by
    // the inverse of the packing permutation in Q. This allows the product
    // of VPacked with Q to be equal to the unpacked V times the right singular
    // vectors from the secular equation.
    if( ctrl.progress )
        Output("Forming undeflated right singular vectors");
    Q.Resize( numUndeflated, numUndeflated );
    for( Int j=0; j<numUndeflated; ++j )
    {
        auto v = VSecular(ALL,IR(j));
        auto q = Q(ALL,IR(j));
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
    //   Z Q = Z_{:,2} Q2 + | Z_{0,0} Q_0 |.
    //                      |-------------|
    //                      | Z_{1,1} Q_1 |
    //
    if( ctrl.progress )
        Output("Overwriting right singular vectors");
    auto VUndeflated = V( ALL, undeflatedInd );
    if( ctrl.wantV )
    {
        if( dcCtrl.exploitStructure )
        {
            auto Z2 = VPacked( ALL, packingInd2 );
            auto Q2 = Q( packingInd2, ALL );
            Gemm( NORMAL, NORMAL, Real(1), Z2, Q2, VUndeflated );

            // Finish updating the first block row
            auto V0Undeflated = VUndeflated( IR(0,n0), ALL );
            auto Z00 = VPacked( IR(0,n0), packingInd0 );
            auto Q0 = Q( packingInd0, ALL );
            Gemm( NORMAL, NORMAL, Real(1), Z00, Q0, Real(1), V0Undeflated );

            // Finish updating the second block row
            auto V1Undeflated = VUndeflated( IR(n0,n), ALL );
            auto Z11 = VPacked( IR(n0,n), packingInd1 );
            auto Q1 = Q( packingInd1, ALL );
            Gemm( NORMAL, NORMAL, Real(1), Z11, Q1, Real(1), V1Undeflated );
        }
        else
        {
            Gemm( NORMAL, NORMAL, Real(1), VPacked, Q, VUndeflated );
        }
    }
    else
    {
        if( dcCtrl.exploitStructure )
        {
            auto Z2 = VPacked( ALL, packingInd2 );
            auto Q2 = Q( packingInd2, ALL );
            Gemm( NORMAL, NORMAL, Real(1), Z2, Q2, VUndeflated );

            // Finish updating the first block row
            auto V0Undeflated = VUndeflated( IR(0), ALL );
            auto Z00 = VPacked( IR(0), packingInd0 );
            auto Q0 = Q( packingInd0, ALL );
            Gemm( NORMAL, NORMAL, Real(1), Z00, Q0, Real(1), V0Undeflated );

            // Finish updating the second block row
            auto V1Undeflated = VUndeflated( IR(1), ALL );
            auto Z11 = VPacked( IR(1), packingInd1 );
            auto Q1 = Q( packingInd1, ALL );
            Gemm( NORMAL, NORMAL, Real(1), Z11, Q1, Real(1), V1Undeflated );
        }
        else
        {
            Gemm( NORMAL, NORMAL, Real(1), VPacked, Q, VUndeflated );
        }
    }

    // Rescale the singular values
    SafeScale( scale, Real(1), d );

    return info;
}

template<typename Real>
DCInfo
Merge
( Int m,
  Int n,
  Real alpha,
  // The right entry in the removed middle row of the bidiagonal matrix
  Real beta,
  // The non-deflated m0 (unsorted) singular values from B0.
  const DistMatrix<Real,STAR,STAR>& s0,
  // The non-deflated m1 (unsorted) singular values from B1.
  const DistMatrix<Real,STAR,STAR>& s1,
  // If ctrl.wantU is true, then, on entry, U contains a packing of the left
  // singular vectors from the two subproblems,
  //
  //   U = | U0, 0, 0  |,
  //       | 0,  1, 0  |
  //       | 0,  0, U1 |
  //
  // where U0 is m0 x (m0+1) and U1 is either m1 x m1 or m1 x (m1+1).
  //
  // If ctrl.wantU is true, then, on exit, the left singular vectors of the
  // merged bidiagonal matrix.
  DistMatrix<Real>& U,
  // On exit, the (unsorted) singular values of the merged bidiagonal matrix
  DistMatrix<Real,STAR,STAR>& d,
  // If ctrl.wantV is true, then, on entry, a packing of the right singular
  // vectors from the two subproblems,
  //
  //   V = | V0, 0  |,
  //       | 0,  V1 |
  //
  // where V0 is (m0+1) x (m0+1), with its last column lying in the null space
  // of B0, and V1 is either m1 x m1 or (m1+1) x (m1+1), where, in the latter
  // case, its last column must lie in the null space of B1.
  //
  // If ctrl.wantV is false, then, on entry, V is the same as above, but with
  // only the row that goes through the last row of V0 and the row that goes
  // through the first row of V1 kept.
  //
  // If ctrl.wantV is true, on exit, V will contain the right singular vectors
  // of the merged bidiagonal matrix. If ctrl.wantV is false, then only the two
  // rows of the result mentioned above will be output.
  DistMatrix<Real>& V,
  const BidiagSVDCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Grid& g = s0.Grid();
    const bool amRoot = ( g.Rank() == 0 );
    const Int m0 = s0.Height();
    const Int m1 = s1.Height();
    const Int n0 = m0 + 1;
    const Int n1 = n - n0;
    const bool square = ( m1 == n1 );
    EL_DEBUG_ONLY(
      if( !square && n1 != m1+1 )
          LogicError("B1 has to be square or one column wider than tall");
    )
    const auto& dcCtrl = ctrl.dcCtrl;

    DCInfo info;
    auto& secularInfo = info.secularInfo;
    // TODO(poulson): Switch to log files rather than Output due to the fact
    // that using a recursive process subdivision guarantees that the output
    // in a single stream would be garbled.
    if( ctrl.progress && amRoot )
        Output("m=",m,", n=",n,", m0=",m0,", n0=",n0,", m1=",m1,", n1=",n1);

    DistMatrix<Real> V0(g), V1(g);
    if( ctrl.wantV )
    {
        // V = | V0 0  |
        //     |  0 V1 |
        View( V0, V, IR(0,n0), IR(0,n0) );
        View( V1, V, IR(n0,END), IR(n0,END) );
    }
    else
    {
        View( V0, V, IR(0), IR(0,n0) );
        View( V1, V, IR(1), IR(n0,END) );
    }

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
    const auto& s0Loc = s0.LockedMatrix();
    const auto& s1Loc = s1.LockedMatrix();
    auto& dLoc = d.Matrix();
    dLoc(0) = 0;
    for( Int j=0; j<m0; ++j )
    {
        dLoc(j+1) = s0Loc(j);
    }
    for( Int j=0; j<m1; ++j )
    {
        dLoc(j+n0) = s1Loc(j);
    }

    // Compute the scale of the problem and rescale {d,alpha,beta}. We will
    // rescale the singular values at the end of this routine.
    Real scale = Max( Abs(alpha), Abs(beta) );
    scale = Max( scale, MaxNorm(s0Loc) );
    scale = Max( scale, MaxNorm(s1Loc) );
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
    const Real deflationTol = dcCtrl.deflationFudge*eps;

    // Get a full copy of the last row of V0 and the first row of V1.
    const Int lastRowOfV0 = ( ctrl.wantV ? m0 : 0 );
    DistMatrix<Real,STAR,STAR> v0Last( V0(IR(lastRowOfV0),ALL) ),
      v1First( V1(IR(0),ALL) );
    const auto& v0LastLoc = v0Last.LockedMatrix();
    const auto& v1FirstLoc = v1First.LockedMatrix();

    // Put U and V into [VC,STAR] distributions for Givens applications
    DistMatrix<Real,VC,STAR> U_VC_STAR(U), V_VC_STAR(V);
    auto& U_VC_STAR_Loc = U_VC_STAR.Matrix();
    auto& V_VC_STAR_Loc = V_VC_STAR.Matrix();

    Matrix<Real> r(m,1);
    Matrix<Int> columnTypes(m,1);
    // Form the reordered left portion
    r(0) = alpha*v0LastLoc(0,m0);
    columnTypes(0) = DENSE_COLUMN;
    for( Int j=0; j<m0; ++j )
    {
        r(j+1) = alpha*v0LastLoc(0,j);
        columnTypes(j+1) = COLUMN_NONZERO_IN_FIRST_BLOCK;
    }
    for( Int j=0; j<m1; ++j )
    {
        r(j+n0) = beta*v1FirstLoc(0,j);
        columnTypes(j+n0) = COLUMN_NONZERO_IN_SECOND_BLOCK;
    }

    // Form r(m) if B has one more column than row and then compute the cosine
    // and sine defining the Givens rotation for rotating it into r(0). Then
    // ensure that |r(0)| >= deflationTol. The Givens rotation is such that
    //
    //   | r(0), rhoExtra | | cExtra,  -sExtra | = | gamma, 0 |.
    //                      | sExtra,   cExtra |
    if( n == m+1 )
    {
        Real cExtra=1, sExtra=0;
        Real rhoExtra = beta*v1FirstLoc(0,m1);
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
            if( ctrl.wantV )
            {
                for( Int iLoc=0; iLoc<V_VC_STAR.LocalHeight(); ++iLoc )
                {
                    const Int i = V_VC_STAR.GlobalRow(iLoc);
                    if( i < m0+1 )
                    {
                        const Real nu = V_VC_STAR_Loc(iLoc,m0);
                        V_VC_STAR_Loc(iLoc,m0) =  cExtra*nu;
                        V_VC_STAR_Loc(iLoc,m)  = -sExtra*nu;
                    }
                    else
                    {
                        const Real nu = V_VC_STAR_Loc(iLoc,m);
                        V_VC_STAR_Loc(iLoc,m0) = sExtra*nu;
                        V_VC_STAR_Loc(iLoc,m)  = cExtra*nu;
                    }
                }
            }
            else
            {
                for( Int iLoc=0; iLoc<V_VC_STAR.LocalHeight(); ++iLoc )
                {
                    const Int i = V_VC_STAR.GlobalRow(iLoc);
                    if( i == 0 )
                    {
                        const Real nu = V_VC_STAR_Loc(iLoc,m0);
                        V_VC_STAR_Loc(iLoc,m0) =  cExtra*nu;
                        V_VC_STAR_Loc(iLoc,m)  = -sExtra*nu;
                    }
                    else
                    {
                        const Real nu = V_VC_STAR_Loc(iLoc,m);
                        V_VC_STAR_Loc(iLoc,m0) = sExtra*nu;
                        V_VC_STAR_Loc(iLoc,m)  = cExtra*nu;
                    }
                }
            }
            auto vNull_VC_STAR = V_VC_STAR( ALL, IR(m) );
            auto vNull = V( ALL, IR(m) );
            vNull = vNull_VC_STAR;
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
    SortingPermutation( dLoc, combineSortPerm, ASCENDING, stableSort );
    combineSortPerm.PermuteRows( dLoc );
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
    // and sufficiently small diagonal entries.
    Int numDeflated = 0;
    Int numUndeflated = 1; // We do not deflate the first index
    // We will keep track of the last column that we encountered that was not
    // initially deflatable (but that could be deflated later due to close
    // diagonal entries if another undeflatable column is not encountered
    // first).
    Int revivalCandidate = m;
    for( Int j=1; j<m; ++j )
    {
        if( Abs(r(j)) <= deflationTol )
        {
            // We can deflate due to the r component being sufficiently small
            const Int deflationDest = (m-1) - numDeflated;
            deflationPerm.SetImage( j, deflationDest );
            if( ctrl.progress && amRoot )
                Output
                ("Deflating via p(",j,")=",deflationDest,
                 " because |r(",j,")|=|",r(j),"| <= ",deflationTol);
            columnTypes(j) = DEFLATED_COLUMN;
            ++numDeflated;
            if( amRoot )
            {
                ++secularInfo.numDeflations;
                ++secularInfo.numSmallUpdateDeflations;
            }
            continue;
        }
        if( dLoc(j) <= deflationTol )
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
            const Int jOrig = combinedToOrig( j );
            // TODO(poulson): Exploit the nonzero structure of V?
            blas::Rot
            ( V_VC_STAR_Loc.Height(), V_VC_STAR_Loc.Buffer(0,m0), 1,
              V_VC_STAR_Loc.Buffer(0,jOrig), 1, c, s );

            const Int deflationDest = (m-1) - numDeflated;
            deflationPerm.SetImage( j, deflationDest );
            if( ctrl.progress && amRoot )
                Output
                ("Deflating via p(",j,")=",deflationDest,
                 " because d(",j,")=",dLoc(j)," <= ",deflationTol);

            columnTypes(j) = DEFLATED_COLUMN;

            ++numDeflated;
            if( amRoot )
            {
                ++secularInfo.numDeflations;
                ++secularInfo.numCloseDiagonalDeflations;
            }
            continue;
        }

        revivalCandidate = j;
        if( ctrl.progress && amRoot )
            Output("Breaking initial deflation loop at j=",j);
        break;
    }
    // If we already fully deflated, then the following loop should be trivial
    const Int deflationRestart = revivalCandidate+1;
    for( Int j=deflationRestart; j<m; ++j )
    {
        if( Abs(r(j)) <= deflationTol )
        {
            const Int deflationDest = (m-1) - numDeflated;
            deflationPerm.SetImage( j, deflationDest );
            if( ctrl.progress && amRoot )
                Output
                ("Deflating via p(",j,")=",deflationDest,
                 " because |r(",j,")|=|",r(j),"| <= ",deflationTol);
            columnTypes(j) = DEFLATED_COLUMN;
            ++numDeflated;
            if( amRoot )
            {
                ++secularInfo.numDeflations;
                ++secularInfo.numSmallUpdateDeflations;
            }
            continue;
        }
        if( dLoc(j)-dLoc(revivalCandidate) <= deflationTol )
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
            //   |  c,  s | |        r(j)         | = | gamma |.
            //   | -s,  c | | r(revivalCandidate) |   |   0   |
            //
            // TODO(poulson): Switch to the more stringent deflation criterion
            // used for the symmetric tridiagonal D&C.
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
            const Int revivalOrig = combinedToOrig( revivalCandidate );
            const Int jOrig = combinedToOrig( j );
            if( ctrl.wantU )
            {
                // TODO(poulson): Exploit the nonzero structure of U?
                blas::Rot
                ( U_VC_STAR_Loc.Height(), U_VC_STAR_Loc.Buffer(0,jOrig), 1,
                  U_VC_STAR_Loc.Buffer(0,revivalOrig), 1, c, s );
            }
            // TODO(poulson): Exploit the nonzero structure of V?
            blas::Rot
            ( V_VC_STAR_Loc.Height(), V_VC_STAR_Loc.Buffer(0,jOrig), 1,
              V_VC_STAR_Loc.Buffer(0,revivalOrig), 1, c, s );

            const Int deflationDest = (m-1) - numDeflated;
            deflationPerm.SetImage( revivalCandidate, deflationDest );
            if( ctrl.progress && amRoot )
                Output
                ("Deflating via p(",revivalCandidate,")=",
                 deflationDest," because d(",j,")=",dLoc(j),
                 " - d(",revivalCandidate,")=",dLoc(revivalCandidate)," <= ",
                 deflationTol);

            if( columnTypes(revivalCandidate) != columnTypes(j) )
            {
                // We mixed top and bottom columns so the result is dense.
                columnTypes(j) = DENSE_COLUMN;
            }
            columnTypes(revivalCandidate) = DEFLATED_COLUMN;

            revivalCandidate = j;
            ++numDeflated;
            if( amRoot )
            {
                ++secularInfo.numDeflations;
                ++secularInfo.numCloseDiagonalDeflations;
            }
            continue;
        }

        // We cannot yet deflate index j, so we must give up on the previous
        // revival candidate and then set revivalCandidate = j.
        dUndeflated(numUndeflated) = dLoc(revivalCandidate);
        rUndeflated(numUndeflated) = r(revivalCandidate);
        deflationPerm.SetImage( revivalCandidate, numUndeflated );
        if( ctrl.progress && amRoot )
            Output
            ("Could not deflate with j=",j," and revivalCandidate=",
             revivalCandidate,", so p(",revivalCandidate,")=",
             numUndeflated);
        ++numUndeflated;
        revivalCandidate = j;
    }
    if( revivalCandidate < m )
    {
        // Give up on the revival candidate
        dUndeflated(numUndeflated) = dLoc(revivalCandidate);
        rUndeflated(numUndeflated) = r(revivalCandidate);
        deflationPerm.SetImage( revivalCandidate, numUndeflated );
        if( ctrl.progress && amRoot )
            Output
            ("Final iteration, so p(",revivalCandidate,")=",numUndeflated);
        ++numUndeflated;
    }
    // Now shrink dUndeflated and rUndeflated down to their proper size
    dUndeflated.Resize( numUndeflated, 1 );
    rUndeflated.Resize( numUndeflated, 1 );

    // Count the number of columns of U with each nonzero pattern
    std::vector<Int> packingCounts( NUM_DC_COMBINED_COLUMN_TYPES, 0 );
    for( Int j=0; j<m; ++j )
        ++packingCounts[columnTypes(j)];
    EL_DEBUG_ONLY(
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
        if( ctrl.progress && amRoot )
            Output("packingCounts[",columnType,"]=",packingCounts[columnType]);
    }

    // Set up the index ranges of the three packed column subsets
    const Range<Int> packingInd0(packingOffsets[0],packingOffsets[1]),
                     packingInd1(packingOffsets[1],packingOffsets[2]),
                     packingInd2(packingOffsets[2],packingOffsets[3]);

    Matrix<Real> dPacked;
    DistMatrix<Real,VC,STAR> UPacked(g), VPacked(g);
    dPacked.Resize( m, 1 );
    if( ctrl.wantU )
        UPacked.Resize( m, m );
    if( ctrl.wantV )
        VPacked.Resize( n, m );
    else
        VPacked.Resize( 2, m );
    auto& UPackedLoc = UPacked.Matrix();
    auto& VPackedLoc = VPacked.Matrix();
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

        dPacked(packingDest) = dLoc(j);
        if( ctrl.wantU )
        {
            // TODO(poulson): Exploit the nonzero structure of U?
            blas::Copy
            ( U_VC_STAR_Loc.Height(), U_VC_STAR_Loc.Buffer(0,jOrig), 1,
              UPackedLoc.Buffer(0,packingDest), 1 );
        }
        // TODO(poulson): Exploit the nonzero structure of V?
        blas::Copy
        ( V_VC_STAR_Loc.Height(), V_VC_STAR_Loc.Buffer(0,jOrig), 1,
          VPackedLoc.Buffer(0,packingDest), 1 );
    }

    // Put the deflated columns in their final destination and shrink UPacked
    // and VPacked back down to their final sizes.
    //
    // TODO(poulson): Exploit the nonzero structure of U and V?
    if( numDeflated > 0 )
    {
        blas::Copy
        ( numDeflated, &dPacked(numUndeflated), 1, &dLoc(numUndeflated), 1 );
        if( ctrl.wantU )
        {
            auto UDeflated = U(ALL,IR(numUndeflated,END));
            auto UPackedDeflated = UPacked(ALL,IR(numUndeflated,END));
            UDeflated = UPackedDeflated;
        }
        auto VDeflated = V(ALL,IR(numUndeflated,END));
        auto VPackedDeflated = VPacked(ALL,IR(numUndeflated,END));
        VDeflated = VPackedDeflated;
    }
    if( ctrl.wantU )
        UPacked.Resize( m, numUndeflated );
    if( ctrl.wantV )
        VPacked.Resize( n, numUndeflated );
    else
        VPacked.Resize( 2, numUndeflated );

    // Now compute the updated singular vectors using UPacked/VPacked
    // ==============================================================
    auto undeflatedInd = IR(0,numUndeflated);
    const Real rUndeflatedNorm = FrobeniusNorm( rUndeflated );
    rUndeflated *= Real(1) / rUndeflatedNorm;
    const Real rho = rUndeflatedNorm*rUndeflatedNorm;

    if( ctrl.progress && amRoot )
        Output("Solving secular equation and correcting update vector");
    Matrix<Real> rCorrected;
    Ones( rCorrected, numUndeflated, 1 );

    // Ensure that there is sufficient space for storing the needed singular
    // vectors from the undeflated secular equation. Notice that we *always*
    // need to compute the right singular vectors of the undeflated secular
    // equation.
    DistMatrix<Real,VR,STAR> dSecular(numUndeflated,1,g);
    DistMatrix<Real,STAR,VR> USecular(g), VSecular(g);
    if( ctrl.wantU )
        USecular.Resize( numUndeflated, numUndeflated );
    VSecular.Resize( numUndeflated, numUndeflated );
    auto& dSecularLoc = dSecular.Matrix();
    auto& USecularLoc = USecular.Matrix();
    auto& VSecularLoc = VSecular.Matrix();

    // For temporarily storing dUndeflated + d(j)
    Matrix<Real> plusShift( numUndeflated, 1 );

    const Int numUndeflatedLoc = VSecularLoc.Width();
    for( Int jLoc=0; jLoc<numUndeflatedLoc; ++jLoc )
    {
        const Int j = VSecular.GlobalCol(jLoc);
        auto minusShift = VSecularLoc( ALL, IR(jLoc) );

        auto valueInfo =
          SecularSingularValue
          ( j, dUndeflated, rho, rUndeflated, dSecularLoc(jLoc),
            minusShift, plusShift, dcCtrl.secularCtrl );
        if( ctrl.progress && amRoot )
            Output("Secular singular value ",j," is ",dSecularLoc(jLoc));

        // We will sum these across all of the processors at the top-level
        secularInfo.numIterations += valueInfo.numIterations;
        secularInfo.numAlternations += valueInfo.numAlternations;
        secularInfo.numCubicIterations += valueInfo.numCubicIterations;
        secularInfo.numCubicFailures += valueInfo.numCubicFailures;

        // minusShift currently holds dUndeflated-d(j) and plusShift
        // holds dUndeflated+d(j). Overwrite minusShift with their
        // element-wise product since that is all we require from here on
        // out.
        for( Int k=0; k<numUndeflated; ++k )
            minusShift(k) *= plusShift(k);

        rCorrected(j) *= minusShift(j);
        for( Int k=0; k<numUndeflated; ++k )
        {
            if( k == j )
                continue;
            rCorrected(k) *= minusShift(k) /
              ((dUndeflated(j)+dUndeflated(k))*
               (dUndeflated(j)-dUndeflated(k)));
        }
    }
    AllReduce( rCorrected, g.VRComm(), mpi::PROD );
    for( Int j=0; j<numUndeflated; ++j )
        rCorrected(j) = Sgn(rUndeflated(j),false) * Sqrt(Abs(rCorrected(j)));
    // Get a full copy of the undeflated singular values now
    {
        auto sUndeflated = d( IR(0,numUndeflated), ALL );
        sUndeflated = dSecular;
    }

    // Compute the unnormalized left and right singular vectors via Eqs. (3.4)
    // and (3.3), respectively, from Gu/Eisenstat [CITATION].
    if( ctrl.progress && amRoot )
        Output("Computing unnormalized singular vectors");
    if( ctrl.wantU )
    {
        for( Int jLoc=0; jLoc<numUndeflatedLoc; ++jLoc )
        {
            auto u = USecularLoc(ALL,IR(jLoc));
            auto v = VSecularLoc(ALL,IR(jLoc));
            {
                const Real deltaSqMinusShiftSq = v(0);
                v(0) = rCorrected(0) / deltaSqMinusShiftSq;
                u(0) = -1;
            }
            for( Int i=1; i<numUndeflated; ++i )
            {
                const Real deltaSqMinusShiftSq = v(i);
                v(i) = rCorrected(i) / deltaSqMinusShiftSq;
                u(i) = dUndeflated(i) * v(i);
            }
        }
    }
    else
    {
        for( Int jLoc=0; jLoc<numUndeflatedLoc; ++jLoc )
        {
            auto v = VSecularLoc(ALL,IR(jLoc));
            {
                const Real deltaSqMinusShiftSq = v(0);
                v(0) = rCorrected(0) / deltaSqMinusShiftSq;
            }
            for( Int i=1; i<numUndeflated; ++i )
            {
                const Real deltaSqMinusShiftSq = v(i);
                v(i) = rCorrected(i) / deltaSqMinusShiftSq;
            }
        }
    }

    // Form the normalized left singular vectors with the rows permuted by
    // the inverse of the packing permutation in Q. This allows the product
    // of UPacked with Q to be equal to the unpacked U times the left singular
    // vectors from the secular equation.
    if( ctrl.progress && amRoot )
        Output("Forming undeflated left singular vectors");
    DistMatrix<Real,STAR,VR> Q(g);
    auto& QLoc = Q.Matrix();
    if( ctrl.wantU )
    {
        Zeros( Q, numUndeflated, numUndeflated );
        for( Int jLoc=0; jLoc<numUndeflatedLoc; ++jLoc )
        {
            auto u = USecularLoc(ALL,IR(jLoc));
            auto q = QLoc(ALL,IR(jLoc));
            const Real uFrob = FrobeniusNorm( u );
            for( Int i=0; i<numUndeflated; ++i )
                q(i) = u(packingPerm.Preimage(i)) / uFrob;
        }
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
    if( ctrl.wantU )
    {
        if( ctrl.progress && amRoot )
            Output("Overwriting left singular vectors");
        auto UUndeflated = U( ALL, undeflatedInd );
        if( dcCtrl.exploitStructure )
        {
            auto Z2 = UPacked( ALL, packingInd2 );
            auto Q2 = Q( packingInd2, ALL );
            Gemm( NORMAL, NORMAL, Real(1), Z2, Q2, UUndeflated );

            // Finish updating the first block row
            auto U0Undeflated = UUndeflated( IR(0,m0), ALL );
            auto Z00 = UPacked( IR(0,m0), packingInd0 );
            auto Q0 = Q( packingInd0, ALL );
            Gemm( NORMAL, NORMAL, Real(1), Z00, Q0, Real(1), U0Undeflated );

            // Finish updating the last block row
            auto U2Undeflated = UUndeflated( IR(n0,m), ALL );
            auto Z21 = UPacked( IR(n0,m), packingInd1 );
            auto Q1 = Q( packingInd1, ALL );
            Gemm( NORMAL, NORMAL, Real(1), Z21, Q1, Real(1), U2Undeflated );
        }
        else
        {
            Gemm( NORMAL, NORMAL, Real(1), UPacked, Q, UUndeflated );
        }
    }

    // Form the normalized right singular vectors with the rows permuted by
    // the inverse of the packing permutation in Q. This allows the product
    // of VPacked with Q to be equal to the unpacked V times the right singular
    // vectors from the secular equation.
    if( ctrl.progress && amRoot )
        Output("Forming undeflated right singular vectors");
    Q.Resize( numUndeflated, numUndeflated );
    for( Int jLoc=0; jLoc<numUndeflatedLoc; ++jLoc )
    {
        auto v = VSecularLoc(ALL,IR(jLoc));
        auto q = QLoc(ALL,IR(jLoc));
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
    //   Z Q = Z_{:,2} Q2 + | Z_{0,0} Q_0 |.
    //                      |-------------|
    //                      | Z_{1,1} Q_1 |
    //
    if( ctrl.progress && amRoot )
        Output("Overwriting right singular vectors");
    auto VUndeflated = V( ALL, undeflatedInd );
    if( ctrl.wantV )
    {
        if( dcCtrl.exploitStructure )
        {
            auto Z2 = VPacked( ALL, packingInd2 );
            auto Q2 = Q( packingInd2, ALL );
            Gemm( NORMAL, NORMAL, Real(1), Z2, Q2, VUndeflated );

            // Finish updating the first block row
            auto V0Undeflated = VUndeflated( IR(0,n0), ALL );
            auto Z00 = VPacked( IR(0,n0), packingInd0 );
            auto Q0 = Q( packingInd0, ALL );
            Gemm( NORMAL, NORMAL, Real(1), Z00, Q0, Real(1), V0Undeflated );

            // Finish updating the second block row
            auto V1Undeflated = VUndeflated( IR(n0,n), ALL );
            auto Z11 = VPacked( IR(n0,n), packingInd1 );
            auto Q1 = Q( packingInd1, ALL );
            Gemm( NORMAL, NORMAL, Real(1), Z11, Q1, Real(1), V1Undeflated );
        }
        else
        {
            Gemm( NORMAL, NORMAL, Real(1), VPacked, Q, VUndeflated );
        }
    }
    else
    {
        if( dcCtrl.exploitStructure )
        {
            auto Z2 = VPacked( ALL, packingInd2 );
            auto Q2 = Q( packingInd2, ALL );
            Gemm( NORMAL, NORMAL, Real(1), Z2, Q2, VUndeflated );

            // Finish updating the first block row
            auto V0Undeflated = VUndeflated( IR(0), ALL );
            auto Z00 = VPacked( IR(0), packingInd0 );
            auto Q0 = Q( packingInd0, ALL );
            Gemm( NORMAL, NORMAL, Real(1), Z00, Q0, Real(1), V0Undeflated );

            // Finish updating the second block row
            auto V1Undeflated = VUndeflated( IR(1), ALL );
            auto Z11 = VPacked( IR(1), packingInd1 );
            auto Q1 = Q( packingInd1, ALL );
            Gemm( NORMAL, NORMAL, Real(1), Z11, Q1, Real(1), V1Undeflated );
        }
        else
        {
            Gemm( NORMAL, NORMAL, Real(1), VPacked, Q, VUndeflated );
        }
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
  const BidiagSVDCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Int m = mainDiag.Height();
    const Int n = superDiag.Height() + 1;
    const auto& dcCtrl = ctrl.dcCtrl;

    DCInfo info;
    auto& secularInfo = info.secularInfo;

    if( m <= Max(dcCtrl.cutoff,3) )
    {
        auto ctrlMod( ctrl );
        ctrlMod.useQR = true;
        ctrlMod.approach = FULL_SVD; // We need any null space of V as well
        if( ctrl.wantV )
        {
            ctrlMod.accumulateU = false;
            ctrlMod.accumulateV = false;
            BidiagSVD( UPPER, mainDiag, superDiag, U, s, V, ctrlMod );
        }
        else
        {
            ctrlMod.accumulateU = false;
            ctrlMod.wantV = true;
            ctrlMod.accumulateV = true;
            Zeros( V, 2, n );
            V(0,0) = 1;
            V(1,n-1) = 1;
            BidiagSVD( UPPER, mainDiag, superDiag, U, s, V, ctrlMod );
        }
        return info;
    }

    // TODO(poulson): A more intelligent split point. Perhaps the row near
    // m/2 with smallest norm should be chosen to encourage small entry
    // deflation (which tends to be vastly more common).
    const Int split = m/2;

    const Real& alpha = mainDiag(split);
    const Real& beta = superDiag(split);

    auto mainDiag0 = mainDiag( IR(0,split), ALL );
    auto superDiag0 = superDiag( IR(0,split), ALL );

    auto mainDiag1 = mainDiag( IR(split+1,END), ALL );
    auto superDiag1 = superDiag( IR(split+1,END), ALL );

    Matrix<Real> U0, U1;
    if( ctrl.wantU )
    {
        Identity( U, m, m );
        View( U0, U, IR(0,split), IR(0,split) );
        View( U1, U, IR(split+1,END), IR(split+1,END) );
    }

    Matrix<Real> V0, V1;
    if( ctrl.wantV )
    {
        Zeros( V, n, n );
        View( V0, V, IR(0,split+1), IR(0,split+1) );
        View( V1, V, IR(split+1,END), IR(split+1,END) );
    }
    else
    {
        Zeros( V0, 2, split+1 );
        Zeros( V1, 2, n-(split+1) );
    }

    Matrix<Real> s0;
    auto info0 = DivideAndConquer( mainDiag0, superDiag0, U0, s0, V0, ctrl );

    Matrix<Real> s1;
    auto info1 = DivideAndConquer( mainDiag1, superDiag1, U1, s1, V1, ctrl );

    if( !ctrl.wantV )
    {
        // We must manually pack the last row of V0 and the first row of V1
        Zeros( V, 2, n );
        auto V0Last = V( IR(0), IR(0,split+1) );
        auto V1First = V( IR(1), IR(split+1,n) );
        V0Last = V0( IR(1), ALL );
        V1First = V1( IR(0), ALL );
    }
    info = Merge( m, n, alpha, beta, s0, s1, U, s, V, ctrl );

    secularInfo.numIterations += info0.secularInfo.numIterations;
    secularInfo.numAlternations += info0.secularInfo.numAlternations;
    secularInfo.numCubicIterations += info0.secularInfo.numCubicIterations;
    secularInfo.numCubicFailures += info0.secularInfo.numCubicFailures;
    secularInfo.numDeflations += info0.secularInfo.numDeflations;
    secularInfo.numSmallDiagonalDeflations +=
      info0.secularInfo.numSmallDiagonalDeflations;
    secularInfo.numCloseDiagonalDeflations +=
      info0.secularInfo.numCloseDiagonalDeflations;
    secularInfo.numSmallUpdateDeflations +=
      info0.secularInfo.numSmallUpdateDeflations;

    secularInfo.numIterations += info1.secularInfo.numIterations;
    secularInfo.numAlternations += info1.secularInfo.numAlternations;
    secularInfo.numCubicIterations += info1.secularInfo.numCubicIterations;
    secularInfo.numCubicFailures += info1.secularInfo.numCubicFailures;
    secularInfo.numDeflations += info1.secularInfo.numDeflations;
    secularInfo.numSmallDiagonalDeflations +=
      info1.secularInfo.numSmallDiagonalDeflations;
    secularInfo.numCloseDiagonalDeflations +=
      info1.secularInfo.numCloseDiagonalDeflations;
    secularInfo.numSmallUpdateDeflations +=
      info1.secularInfo.numSmallUpdateDeflations;

    return info;
}

// We pass in mainDiag and superDiag in sequential form to avoid potential
// confusion from avoiding unnecessarily creating separate copies distributed
// (trivially) over subtrees.
template<typename Real>
DCInfo
DivideAndConquer
( const Matrix<Real>& mainDiag,
  const Matrix<Real>& superDiag,
        DistMatrix<Real>& U,
        DistMatrix<Real,STAR,STAR>& s,
        DistMatrix<Real>& V,
  const BidiagSVDCtrl<Real>& ctrl,
  bool topLevel=true )
{
    EL_DEBUG_CSE
    const Grid& grid = U.Grid();
    const Int m = mainDiag.Height();
    const Int n = superDiag.Height() + 1;
    const auto& dcCtrl = ctrl.dcCtrl;
    DCInfo info;

    if( m <= Max(dcCtrl.cutoff,3) )
    {
        // Run the problem locally
        auto ctrlMod( ctrl );
        ctrlMod.useQR = true;
        ctrlMod.approach = FULL_SVD; // We need any null space of V as well
        Matrix<Real> ULoc, sLoc, VLoc;
        if( ctrl.wantV )
        {
            ctrlMod.accumulateU = false;
            ctrlMod.accumulateV = false;
            if( ctrl.qrCtrl.broadcast )
            {
                if( grid.VCRank() == 0 )
                {
                    BidiagSVD
                    ( UPPER, mainDiag, superDiag, ULoc, sLoc, VLoc, ctrlMod );
                }
                if( ctrl.wantU )
                {
                    ULoc.Resize( m, m );
                    El::Broadcast( ULoc, grid.VCComm(), 0 );
                }
                sLoc.Resize( m, 1 );
                El::Broadcast( sLoc, grid.VCComm(), 0 );
                VLoc.Resize( n, n );
                El::Broadcast( VLoc, grid.VCComm(), 0 );
            }
            else
            {
                // Let's cross our fingers and ignore the forward instability
                BidiagSVD
                ( UPPER, mainDiag, superDiag, ULoc, sLoc, VLoc, ctrlMod );
            }
        }
        else
        {
            ctrlMod.accumulateU = false;
            ctrlMod.wantV = true;
            ctrlMod.accumulateV = true;
            Zeros( VLoc, 2, n );
            VLoc(0,0) = 1;
            VLoc(1,n-1) = 1;
            if( ctrl.qrCtrl.broadcast )
            {
                if( grid.VCRank() == 0 )
                {
                    BidiagSVD
                    ( UPPER, mainDiag, superDiag, ULoc, sLoc, VLoc, ctrlMod );
                }
                if( ctrl.wantU )
                {
                    ULoc.Resize( m, m );
                    El::Broadcast( ULoc, grid.VCComm(), 0 );
                }
                sLoc.Resize( m, 1 );
                El::Broadcast( sLoc, grid.VCComm(), 0 );
                El::Broadcast( VLoc, grid.VCComm(), 0 );
            }
            else
            {
                // Let's cross our fingers and ignore the forward instability
                BidiagSVD
                ( UPPER, mainDiag, superDiag, ULoc, sLoc, VLoc, ctrlMod );
            }
        }

        s.Resize( m, 1 );
        for( Int i=0; i<m; ++i )
            s.SetLocal( i, 0, sLoc(i) );

        if( ctrl.wantU )
        {
            U.Resize( m, m );
            const Int ULocHeight = U.LocalHeight();
            const Int ULocWidth = U.LocalWidth();
            for( Int iLoc=0; iLoc<ULocHeight; ++iLoc )
            {
                const Int i = U.GlobalRow(iLoc);
                for( Int jLoc=0; jLoc<ULocWidth; ++jLoc )
                {
                    const Int j = U.GlobalCol(jLoc);
                    U.SetLocal( iLoc, jLoc, ULoc(i,j) );
                }
            }
        }

        if( ctrl.wantV )
            V.Resize( n, n );
        else
            V.Resize( 2, n );
        const Int VLocHeight = V.LocalHeight();
        const Int VLocWidth = V.LocalWidth();
        for( Int iLoc=0; iLoc<VLocHeight; ++iLoc )
        {
            const Int i = V.GlobalRow(iLoc);
            for( Int jLoc=0; jLoc<VLocWidth; ++jLoc )
            {
                const Int j = V.GlobalCol(jLoc);
                V.SetLocal( iLoc, jLoc, VLoc(i,j) );
            }
        }

        return info;
    }
    else if( grid.Size() == 1 )
    {
        // TODO(poulson): Avoid unnecessary copies
        Matrix<Real> ULoc, sLoc, VLoc;
        info = DivideAndConquer( mainDiag, superDiag, ULoc, sLoc, VLoc, ctrl );

        s.Resize( m, 1 );
        s.Matrix() = sLoc;

        U.Resize( ULoc.Height(), ULoc.Width() );
        U.Matrix() = ULoc;

        V.Resize( VLoc.Height(), VLoc.Width() );
        V.Matrix() = VLoc;

        return info;
    }

    // TODO(poulson): A more intelligent split point. Perhaps the row near
    // m/2 with smallest norm should be chosen to encourage small entry
    // deflation (which tends to be vastly more common).
    const Int split = m/2;
    const Grid *leftGrid, *rightGrid;
    SplitGrid( split, m-(split+1), grid, leftGrid, rightGrid );

    const Real alpha = mainDiag(split);
    const Real beta = superDiag(split);

    auto mainDiag0 = mainDiag( IR(0,split), ALL );
    auto superDiag0 = superDiag( IR(0,split), ALL );

    auto mainDiag1 = mainDiag( IR(split+1,END), ALL );
    auto superDiag1 = superDiag( IR(split+1,END), ALL );

    // Split U and V between the subtrees
    DistMatrix<Real> U0Sub(*leftGrid), U1Sub(*rightGrid);
    if( ctrl.wantU )
    {
        Zeros( U0Sub, split, split );
        Zeros( U1Sub, m-(split+1), m-(split+1) );
    }

    DistMatrix<Real> V0Sub(*leftGrid), V1Sub(*rightGrid);
    if( ctrl.wantV )
    {
        Zeros( V0Sub, split+1, split+1 );
        Zeros( V1Sub, n-(split+1), n-(split+1) );
    }
    else
    {
        Zeros( V0Sub, 2, split+1 );
        Zeros( V1Sub, 2, n-(split+1) );
    }

    DistMatrix<Real,STAR,STAR> s0Sub(*leftGrid), s1Sub(*rightGrid);
    DCInfo info0, info1;
    if( s0Sub.Participating() )
    {
        info0 =
          DivideAndConquer
          ( mainDiag0, superDiag0, U0Sub, s0Sub, V0Sub, ctrl, false );
    }
    if( s1Sub.Participating() )
    {
        info1 =
          DivideAndConquer
          ( mainDiag1, superDiag1, U1Sub, s1Sub, V1Sub, ctrl, false );
    }

    const bool includeViewers = true;
    U0Sub.MakeConsistent( includeViewers );
    s0Sub.MakeConsistent( includeViewers );
    V0Sub.MakeConsistent( includeViewers );
    U1Sub.MakeConsistent( includeViewers );
    s1Sub.MakeConsistent( includeViewers );
    V1Sub.MakeConsistent( includeViewers );

    DistMatrix<Real,STAR,STAR> s0(grid), s1(grid);
    s0 = s0Sub;
    s1 = s1Sub;
    if( ctrl.wantU )
    {
        Identity( U, m, m );
        auto U0 = U(IR(0,split),IR(0,split));
        auto U1 = U(IR(split+1,END),IR(split+1,END));
        U0 = U0Sub;
        U1 = U1Sub;
    }
    if( ctrl.wantV )
    {
        Zeros( V, n, n );
        auto V0 = V(IR(0,split+1),IR(0,split+1));
        auto V1 = V(IR(split+1,END),IR(split+1,END));
        V0 = V0Sub;
        V1 = V1Sub;
    }
    else
    {
        // We must manually pack the last row of V0 and the first row of V1
        Zeros( V, 2, n );
        auto V0Last = V( IR(0), IR(0,split+1) );
        auto V1First = V( IR(1), IR(split+1,n) );
        V0Last = V0Sub( IR(1), ALL );
        V1First = V1Sub( IR(0), ALL );
    }
    info = Merge( m, n, alpha, beta, s0, s1, U, s, V, ctrl );

    auto& secularInfo = info.secularInfo;
    if( s0Sub.Participating() )
    {
        secularInfo.numIterations += info0.secularInfo.numIterations;
        secularInfo.numAlternations += info0.secularInfo.numAlternations;
        secularInfo.numCubicIterations += info0.secularInfo.numCubicIterations;
        secularInfo.numCubicFailures += info0.secularInfo.numCubicFailures;
    }
    if( s1Sub.Participating() )
    {
        secularInfo.numIterations += info1.secularInfo.numIterations;
        secularInfo.numAlternations += info1.secularInfo.numAlternations;
        secularInfo.numCubicIterations += info1.secularInfo.numCubicIterations;
        secularInfo.numCubicFailures += info1.secularInfo.numCubicFailures;
    }

    if( leftGrid->Rank() == 0 )
    {
        secularInfo.numDeflations += info0.secularInfo.numDeflations;
        secularInfo.numSmallDiagonalDeflations +=
          info0.secularInfo.numSmallDiagonalDeflations;
        secularInfo.numCloseDiagonalDeflations +=
          info0.secularInfo.numCloseDiagonalDeflations;
        secularInfo.numSmallUpdateDeflations +=
          info0.secularInfo.numSmallUpdateDeflations;
    }
    if( rightGrid->Rank() == 0 )
    {
        secularInfo.numDeflations += info1.secularInfo.numDeflations;
        secularInfo.numSmallDiagonalDeflations +=
          info1.secularInfo.numSmallDiagonalDeflations;
        secularInfo.numCloseDiagonalDeflations +=
          info1.secularInfo.numCloseDiagonalDeflations;
        secularInfo.numSmallUpdateDeflations +=
          info1.secularInfo.numSmallUpdateDeflations;
    }

    if( topLevel )
    {
        // Sum all of the secular iteration/deflation information
        Matrix<Int> counts(7,1);
        counts(0) = secularInfo.numIterations;
        counts(1) = secularInfo.numAlternations;
        counts(2) = secularInfo.numCubicIterations;
        counts(3) = secularInfo.numCubicFailures;
        counts(4) = secularInfo.numSmallDiagonalDeflations;
        counts(5) = secularInfo.numCloseDiagonalDeflations;
        counts(6) = secularInfo.numSmallUpdateDeflations;

        AllReduce( counts, grid.Comm() );

        secularInfo.numIterations = counts(0);
        secularInfo.numAlternations = counts(1);
        secularInfo.numCubicIterations = counts(2);
        secularInfo.numCubicFailures = counts(3);
        secularInfo.numSmallDiagonalDeflations = counts(4);
        secularInfo.numCloseDiagonalDeflations = counts(5);
        secularInfo.numSmallUpdateDeflations = counts(6);
        secularInfo.numDeflations = counts(4) + counts(5) + counts(6);
    }

    return info;
}

} // namespace bidiag_svd
} // namespace El

#endif // ifndef EL_BIDIAG_SVD_DC_HPP
