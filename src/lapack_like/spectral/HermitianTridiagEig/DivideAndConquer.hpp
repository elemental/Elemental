/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_HERM_TRIDIAG_EIG_DC_HPP
#define EL_HERM_TRIDIAG_EIG_DC_HPP

// TODO(poulson): Move said routine into a utility function
#include "../Schur/SDC.hpp"
using El::schur::SplitGrid;

namespace El {
namespace herm_tridiag_eig {

// The following is analogous to LAPACK's {s,d}laed{1,2,3} [CITATION] but does
// not accept initial sorting permutations for w0 and w1, nor does it enforce 
// any ordering on the resulting eigenvalues.
template<typename Real>
DCInfo
Merge
( Real beta,
  // The n0 (unsorted) eigenvalues from T0.
  const Matrix<Real>& w0,
  // The n1 (unsorted) eigenvalues from T1.
  const Matrix<Real>& w1,
  // On exit, the (unsorted) eigenvalues of the merged tridiagonal matrix
  Matrix<Real>& d,
  // If ctrl.wantEigVecs is true, then, on entry, a packing of the eigenvectors
  // from the two subproblems,
  //
  //   Q = | Q0, 0  |,
  //       | 0,  Q1 |
  //
  // where Q0 is n0 x n0, and Q1 is n1 x n1.
  //
  // If ctrl.wantEigVecs is false, then, on entry, Q is the same as above, but
  // with only the row that goes through the last row of Q0 and the row that
  // goes through the first row of Q1 kept.
  //
  // If ctrl.wantEigVecs is true, on exit, Q will contain the eigenvectors of
  // the merged tridiagonal matrix. If ctrl.wantEigVecs is false, then only the
  // two rows of the result mentioned above will be output.
  Matrix<Real>& Q,
  const HermitianTridiagEigCtrl<Real>& ctrl )
{
    DEBUG_CSE
    const Int n0 = w0.Height();
    const Int n1 = w1.Height();
    const Int n = n0 + n1;
    const auto& dcCtrl = ctrl.dcCtrl;

    DCInfo info;
    auto& secularInfo = info.secularInfo;
    if( ctrl.progress )
        Output("n=",n,", n0=",n0,", n1=",n1);

    Matrix<Real> Q0, Q1;
    if( ctrl.wantEigVecs )
    {
        // Q = | Q0 0  |
        //     |  0 Q1 |
        View( Q0, Q, IR(0,n0), IR(0,n0) );
        View( Q1, Q, IR(n0,END), IR(n0,END) );
    }
    else
    {
        View( Q0, Q, IR(0), IR(0,n0) );
        View( Q1, Q, IR(1), IR(n0,END) );
    }

    // Before permutation, 
    //
    //   r = sqrt(2 |beta|) z,
    //
    // where
    //     
    //   z = [ sgn(beta)*Q0(n0-1,:), Q1(0,:) ] / sqrt(2).
    //
    // But we reorder indices 0 and n0-1 to put r in the first position. Thus,
    // we must form 
    //
    //   d = [w0(n0-1); w0(0:n0-2); w1]
    //
    // and consider the matrix
    //
    //   diag(d) + 2 |beta| z z'.
    //

    // Form d = [w0(n0-1); w0(0:n0-2); w1].
    // This effectively cyclically shifts [0,n0) |-> [1,n0+1) mod n0.
    d.Resize( n, 1 );
    d(0) = w0(n0-1);
    for( Int j=0; j<n0-1; ++j )
    {
        d(j+1) = w0(j);
    }
    for( Int j=0; j<n1; ++j )
    {
        d(j+n0) = w1(j);
    }

    // Compute the scale of the problem and rescale. We will rescale the 
    // eigenvalues at the end of this routine. Note that LAPACK's {s,d}laed2
    // [CITATION] uses max(|beta|,||z||_max), where || z ||_2 = sqrt(2),
    // which could be much too small if || r ||_2 is much larger than beta 
    // and sqrt(2).
    Real scale = Max( 2*Abs(beta), MaxNorm(d) );
    SafeScale( Real(1), scale, d );
    SafeScale( Real(1), scale, beta );

    // Now that the problem is rescaled, our deflation tolerance simplifies to
    //
    //   tol = deflationFudge eps max( || d ||_max, 2*|beta| )
    //       = deflationFudge eps.
    //
    // Cf. LAPACK's {s,d}lasd2 [CITATION] for this tolerance.
    const Real eps = limits::Epsilon<Real>();
    const Real deflationTol = dcCtrl.deflationFudge*eps;

    Matrix<Real> z(n,1);
    Matrix<Int> columnTypes(n,1);
    const Real betaSgn = Sgn( beta, false );
    const Int lastRowOfQ0 = ( ctrl.wantEigVecs ? n0-1 : 0 );
    const Real sqrtTwo = Sqrt( Real(2) );
    z(0) = betaSgn*Q0(lastRowOfQ0,n0-1) / sqrtTwo;
    columnTypes(0) = DENSE_COLUMN;
    for( Int j=0; j<n0-1; ++j )
    {
        z(j+1) = betaSgn*Q0(lastRowOfQ0,j) / sqrtTwo;
        columnTypes(j+1) = COLUMN_NONZERO_IN_FIRST_BLOCK;
    }
    for( Int j=0; j<n1; ++j )
    {
        z(j+n0) = Q1(0,j) / sqrtTwo;
        columnTypes(j+n0) = COLUMN_NONZERO_IN_SECOND_BLOCK;
    }

    Permutation combineSortPerm;
    SortingPermutation( d, combineSortPerm, ASCENDING );
    combineSortPerm.PermuteRows( d );
    combineSortPerm.PermuteRows( z );
    combineSortPerm.PermuteRows( columnTypes );

    auto combinedToOrig = [&]( const Int& combinedIndex )
      {
          const Int preCombined = combineSortPerm.Preimage( combinedIndex );
          if( preCombined < n0 )
              // Undo the cyclic shift [0,n0) |-> [1,n0+1) mod n0 which
              // pushed the removed row into the first position.
              return Mod( preCombined-1, n0 );
          else
              return preCombined;
      };

    Permutation deflationPerm;
    deflationPerm.MakeIdentity( n );
    deflationPerm.MakeArbitrary();
    // Since we do not yet know how many undeflated entries there will be, we
    // must use the no-deflation case as our storage upper bound.
    Matrix<Real> dUndeflated(n,1), zUndeflated(n,1);
    dUndeflated(0) = 0;
    zUndeflated(0) = z(0);

    // Deflate all (off-diagonal) update entries sufficiently close to zero
    Int numDeflated = 0;
    Int numUndeflated = 0;
    // We will keep track of the last column that we encountered that was not
    // initially deflatable (but that could be deflated later due to close
    // diagonal entries if another undeflatable column is not encountered
    // first).
    Int revivalCandidate = n;
    for( Int j=0; j<n; ++j )
    {
        if( Abs(2*beta*z(j)) <= deflationTol )
        {
            // We can deflate due to the r component being sufficiently small
            const Int deflationDest = (n-1) - numDeflated;
            deflationPerm.SetImage( j, deflationDest );
            if( ctrl.progress )
                Output
                ("Deflating via p(",j,")=",deflationDest,
                 " because |2*beta*z(",j,")|=|",2*beta*z(j),"| <= ",
                 deflationTol);
            columnTypes(j) = DEFLATED_COLUMN;
            ++numDeflated;
            ++secularInfo.numDeflations;
            ++secularInfo.numSmallUpdateDeflations;
        }
        else
        {
            revivalCandidate = j;
            if( ctrl.progress )
                Output("Breaking initial deflation loop at j=",j);
            break;
        }
    }
    // If we already fully deflated, then the following loop should be trivial
    const Int deflationRestart = revivalCandidate+1;
    for( Int j=deflationRestart; j<n; ++j )
    {
        if( Abs(2*beta*z(j)) <= deflationTol )
        {
            const Int deflationDest = (n-1) - numDeflated;
            deflationPerm.SetImage( j, deflationDest );
            if( ctrl.progress )
                Output
                ("Deflating via p(",j,")=",deflationDest,
                 " because |2*beta*z(",j,")|=|",2*beta*z(j),"| <= ",
                 deflationTol);
            columnTypes(j) = DEFLATED_COLUMN;
            ++numDeflated;
            ++secularInfo.numDeflations;
            ++secularInfo.numSmallUpdateDeflations;
            continue;
        }
        const Real gamma = SafeNorm( z(j), z(revivalCandidate) );
        const Real c = z(j) / gamma;
        const Real s = z(revivalCandidate) / gamma;
        const Real offDiagNew = c*s*(d(j)-d(revivalCandidate));
        if( Abs(offDiagNew) <= deflationTol )
        {
            // Deflate the previously undeflatable index by rotating
            // z(revivalCandidate) into z(j) (Cf. the discussion
            // surrounding Eq. (4.4) of Gu/Eisenstat's TR [CITATION]).
            //
            // In particular, we want
            //
            //   | z(j), z(revivalCandidate) | | c -s | = | gamma, 0 |,
            //                                 | s  c |
            //
            // where gamma = || z(revivalCandidate); z(j) ||_2. Putting 
            //
            //   c = z(j)                / gamma,
            //   s = z(revivalCandidate) / gamma,
            //
            // implies
            //
            //   |  c,  s | |        z(j)         | = | gamma |.
            //   | -s,  c | | z(revivalCandidate) |   |   0   |
            //
            z(j) = gamma;
            z(revivalCandidate) = 0;

            // Apply | c -s | to both sides of d
            //       | s  c |
            const Real deltaDeflate = d(revivalCandidate)*(c*c) + d(j)*(s*s);
            d(j) = d(j)*(c*c) + d(revivalCandidate)*(s*s);
            d(revivalCandidate) = deltaDeflate;

            // Apply | c -s | from the right to Q
            //       | s  c |
            //
            const Int revivalOrig = combinedToOrig( revivalCandidate );
            const Int jOrig = combinedToOrig( j );
            if( ctrl.wantEigVecs )
            {
                // TODO(poulson): Exploit the nonzero structure of Q?
                blas::Rot( n, &Q(0,jOrig), 1, &Q(0,revivalOrig), 1, c, s );
            }
            else
            {
                blas::Rot( 2, &Q(0,jOrig), 1, &Q(0,revivalOrig), 1, c, s );
            }

            const Int deflationDest = (n-1) - numDeflated;
            deflationPerm.SetImage( revivalCandidate, deflationDest );
            if( ctrl.progress )
                Output
                ("Deflating via p(",revivalCandidate,")=",
                 deflationDest," because |c*s*(d(",j,")-d(",revivalCandidate,
                 "))|=",offDiagNew," <= ",deflationTol);

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
        zUndeflated(numUndeflated) = z(revivalCandidate);
        deflationPerm.SetImage( revivalCandidate, numUndeflated );
        if( ctrl.progress )
            Output
            ("Could not deflate with j=",j," and revivalCandidate=",
             revivalCandidate,", so p(",revivalCandidate,")=",
             numUndeflated);
        ++numUndeflated;
        revivalCandidate = j;
    }
    if( revivalCandidate < n )
    {
        // Give up on the revival candidate
        dUndeflated(numUndeflated) = d(revivalCandidate);
        zUndeflated(numUndeflated) = z(revivalCandidate);
        deflationPerm.SetImage( revivalCandidate, numUndeflated );
        if( ctrl.progress )
            Output
            ("Final iteration, so p(",revivalCandidate,")=",numUndeflated);
        ++numUndeflated;
    }

    // Now shrink dUndeflated and zUndeflated down to their proper size
    dUndeflated.Resize( numUndeflated, 1 );
    zUndeflated.Resize( numUndeflated, 1 );

    // Count the number of columns of Q with each nonzero pattern
    std::vector<Int> packingCounts( NUM_DC_COMBINED_COLUMN_TYPES, 0 );
    for( Int j=0; j<n; ++j )
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
    Matrix<Real> QPacked;
    dPacked.Resize( n, 1 );
    if( ctrl.wantEigVecs )
        QPacked.Resize( n, n );
    else
        QPacked.Resize( 2, n );
    Permutation packingPerm;
    packingPerm.MakeIdentity( n );
    for( Int j=0; j<n; ++j )
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
        if( ctrl.wantEigVecs )
        {
            // TODO(poulson): Exploit the nonzero structure of Q?
            blas::Copy( n, &Q(0,jOrig), 1, &QPacked(0,packingDest), 1 );
        }
        else
        {
            blas::Copy( 2, &Q(0,jOrig), 1, &QPacked(0,packingDest), 1 );
        }
    }

    // Put the deflated columns in their final destination and shrink QPacked
    // down to its final size
    //
    // TODO(poulson): Exploit the nonzero structure of Q?
    if( numDeflated > 0 )
    {
        blas::Copy
        ( numDeflated, &dPacked(numUndeflated), 1, &d(numUndeflated), 1 );
        if( ctrl.wantEigVecs )
        {
            lapack::Copy
            ( 'A', n, numDeflated,
              &QPacked(0,numUndeflated), QPacked.LDim(),
              &Q(0,numUndeflated), Q.LDim() );
        }
        else
        {
            lapack::Copy
            ( 'A', 2, numDeflated,
              &QPacked(0,numUndeflated), QPacked.LDim(),
              &Q(0,numUndeflated), Q.LDim() );
        }
    }
    if( ctrl.wantEigVecs )
        QPacked.Resize( n, numUndeflated );
    else
        QPacked.Resize( 2, numUndeflated );

    // Now compute the updated eigenvectors using QPacked
    // ==================================================
    auto undeflatedInd = IR(0,numUndeflated);
    const Real zUndeflatedNorm = FrobeniusNorm( zUndeflated );
    zUndeflated *= Real(1) / zUndeflatedNorm;
    const Real rho = 2*Abs(beta)*zUndeflatedNorm*zUndeflatedNorm;

    if( ctrl.progress )
        Output("Solving secular equation and correcting update vector");
    Matrix<Real> rCorrected;
    Ones( rCorrected, numUndeflated, 1 );

    // Ensure that there is sufficient space for storing the needed 
    // eigenvectors from the undeflated secular equation. Notice that we
    // *always* need to compute the eigenvectors of the undeflated secular
    // equation.
    Matrix<Real> QSecular;
    if( ctrl.wantEigVecs )
        View( QSecular, Q, undeflatedInd, undeflatedInd );
    else
        QSecular.Resize( numUndeflated, numUndeflated );

    for( Int j=0; j<numUndeflated; ++j )
    {
        auto minusShift = QSecular( ALL, IR(j) );

        auto valueInfo =
          SecularEigenvalue
          ( j, dUndeflated, rho, zUndeflated, d(j), minusShift,
            dcCtrl.secularCtrl );
        if( ctrl.progress )
            Output("Secular eigenvalue ",j," is ",d(j));

        secularInfo.numIterations += valueInfo.numIterations;
        secularInfo.numAlternations += valueInfo.numAlternations;
        secularInfo.numCubicIterations += valueInfo.numCubicIterations;
        secularInfo.numCubicFailures += valueInfo.numCubicFailures;

        rCorrected(j) *= minusShift(j);
        for( Int k=0; k<numUndeflated; ++k )
        {
            if( k == j )
                continue;
            rCorrected(k) *= minusShift(k) / (dUndeflated(j)-dUndeflated(k));
        }
    }
    for( Int j=0; j<numUndeflated; ++j )
        rCorrected(j) = Sgn(zUndeflated(j),false) * Sqrt(Abs(rCorrected(j)));

    // Compute the unnormalized eigenvectors.
    if( ctrl.progress )
        Output("Computing unnormalized eigenvectors");
    for( Int j=0; j<numUndeflated; ++j )
    {
        auto q = QSecular(ALL,IR(j));
        for( Int i=0; i<numUndeflated; ++i )
            q(i) = rCorrected(i) / q(i);
    }

    // Form the normalized right singular vectors with the rows permuted by
    // the inverse of the packing permutation in U. This allows the product
    // of QPacked with U to be equal to the unpacked Q times the eigenvectors
    // from the secular equation.
    Matrix<Real> U;
    if( ctrl.progress )
        Output("Forming undeflated right singular vectors");
    U.Resize( numUndeflated, numUndeflated );
    for( Int j=0; j<numUndeflated; ++j )
    {
        auto q = QSecular(ALL,IR(j));
        auto u = U(ALL,IR(j));
        const Real qFrob = FrobeniusNorm( q );
        for( Int i=0; i<numUndeflated; ++i )
            u(i) = q(packingPerm.Preimage(i)) / qFrob;
    }
    // Overwrite the first 'numUndeflated' columns of Q with the updated 
    // eigenvectors by exploiting the partitioning of Z = QPacked as
    //
    //   Z = | Z_{0,0} |    0    | Z_{0,2} |,
    //       |---------|---------|---------|
    //       |    0    | Z_{1,1} | Z_{1,2} |
    //
    // where the first and second block rows have heights n0 and n1. The block
    // columns respectively have widths packingCounts[0], packingCounts[1], and
    // packingCounts[2].
    //
    // Conformally partitioning U, we have
    //
    //   Z U = Z_{:,2} U2 + | Z_{0,0} U_0 |.
    //                      |-------------|
    //                      | Z_{1,1} U_1 |
    //
    if( ctrl.progress )
        Output("Overwriting eigenvectors");
    auto QUndeflated = Q( ALL, undeflatedInd );
    if( ctrl.wantEigVecs )
    {
        if( dcCtrl.exploitStructure )
        {
            auto Z2 = QPacked( ALL, packingInd2 );
            auto U2 = U( packingInd2, ALL );
            Gemm( NORMAL, NORMAL, Real(1), Z2, U2, QUndeflated );

            // Finish updating the first block row
            auto Q0Undeflated = QUndeflated( IR(0,n0), ALL );
            auto Z00 = QPacked( IR(0,n0), packingInd0 );
            auto U0 = U( packingInd0, ALL );
            Gemm( NORMAL, NORMAL, Real(1), Z00, U0, Real(1), Q0Undeflated );

            // Finish updating the second block row
            auto Q1Undeflated = QUndeflated( IR(n0,n), ALL );
            auto Z11 = QPacked( IR(n0,n), packingInd1 );
            auto U1 = U( packingInd1, ALL );
            Gemm( NORMAL, NORMAL, Real(1), Z11, U1, Real(1), Q1Undeflated );
        }
        else
        {
            Gemm( NORMAL, NORMAL, Real(1), QPacked, U, QUndeflated );
        }
    }
    else
    {
        if( dcCtrl.exploitStructure )
        {
            auto Z2 = QPacked( ALL, packingInd2 );
            auto U2 = U( packingInd2, ALL );
            Gemm( NORMAL, NORMAL, Real(1), Z2, U2, QUndeflated );

            // Finish updating the first block row
            auto Q0Undeflated = QUndeflated( IR(0), ALL );
            auto Z00 = QPacked( IR(0), packingInd0 );
            auto U0 = U( packingInd0, ALL );
            Gemm( NORMAL, NORMAL, Real(1), Z00, U0, Real(1), Q0Undeflated );

            // Finish updating the second block row
            auto Q1Undeflated = QUndeflated( IR(1), ALL );
            auto Z11 = QPacked( IR(1), packingInd1 );
            auto U1 = U( packingInd1, ALL );
            Gemm( NORMAL, NORMAL, Real(1), Z11, U1, Real(1), Q1Undeflated );
        }
        else
        {
            Gemm( NORMAL, NORMAL, Real(1), QPacked, U, QUndeflated );
        }
    }

    // Rescale the eigenvalues
    SafeScale( scale, Real(1), d );

    return info;
}

template<typename Real>
DCInfo
Merge
( Real beta,
  // The n0 (unsorted) eigenvalues from T0.
  const DistMatrix<Real,STAR,STAR>& w0,
  // The n1 (unsorted) eigenvalues from T1.
  const DistMatrix<Real,STAR,STAR>& w1,
  // On exit, the (unsorted) eigenvalues of the merged tridiagonal matrix
  DistMatrix<Real,STAR,STAR>& d,
  // If ctrl.wantEigVecs is true, then, on entry, a packing of the eigenvectors
  // from the two subproblems,
  //
  //   Q = | Q0, 0  |,
  //       | 0,  Q1 |
  //
  // where Q0 is n0 x n0, and Q1 is n1 x n1.
  //
  // If ctrl.wantEigVecs is false, then, on entry, Q is the same as above, but
  // with only the row that goes through the last row of Q0 and the row that
  // goes through the first row of Q1 kept.
  //
  // If ctrl.wantEigVecs is true, on exit, Q will contain the eigenvectors of
  // the merged tridiagonal matrix. If ctrl.wantEigVecs is false, then only the
  // two rows of the result mentioned above will be output.
  DistMatrix<Real>& Q,
  const HermitianTridiagEigCtrl<Real>& ctrl )
{
    DEBUG_CSE
    const Grid& g = w0.Grid();
    const bool amRoot = ( g.Rank() == 0 );
    const Int n0 = w0.Height();
    const Int n1 = w1.Height();
    const Int n = n0 + n1;
    const auto& dcCtrl = ctrl.dcCtrl;

    DCInfo info;
    auto& secularInfo = info.secularInfo;
    // TODO(poulson): Switch to log files rather than Output due to the fact
    // that using a recursive process subdivision guarantees that the output
    // in a single stream would be garbled.
    if( ctrl.progress && amRoot )
        Output("n=",n,", n0=",n0,", n1=",n1);

    DistMatrix<Real> Q0(g), Q1(g);
    if( ctrl.wantEigVecs )
    {
        // Q = | Q0 0  |
        //     |  0 Q1 |
        View( Q0, Q, IR(0,n0), IR(0,n0) );
        View( Q1, Q, IR(n0,END), IR(n0,END) );
    }
    else
    {
        View( Q0, Q, IR(0), IR(0,n0) );
        View( Q1, Q, IR(1), IR(n0,END) );
    }

    // Before permutation, 
    //
    //   r = sqrt(2 |beta|) z,
    //
    // where
    //     
    //   z = [ sgn(beta)*Q0(n0-1,:), Q1(0,:) ] / sqrt(2).
    //
    // But we reorder indices 0 and n0-1 to put r in the first position. Thus,
    // we must form 
    //
    //   d = [w0(n0-1); w0(0:n0-2); w1]
    //
    // and consider the matrix
    //
    //   diag(d) + 2 |beta| z z'.
    //

    // Form d = [w0(n0-1); w0(0:n0-2); w1].
    // This effectively cyclically shifts [0,n0) |-> [1,n0+1) mod n0.
    d.Resize( n, 1 );
    const auto& w0Loc = w0.LockedMatrix();
    const auto& w1Loc = w1.LockedMatrix();
    auto& dLoc = d.Matrix();
    dLoc(0) = w0Loc(n0-1);
    for( Int j=0; j<n0-1; ++j )
    {
        dLoc(j+1) = w0Loc(j);
    }
    for( Int j=0; j<n1; ++j )
    {
        dLoc(j+n0) = w1Loc(j);
    }

    // Compute the scale of the problem and rescale. We will rescale the 
    // eigenvalues at the end of this routine. Note that LAPACK's {s,d}laed2
    // [CITATION] uses max(|beta|,||z||_max), where || z ||_2 = sqrt(2),
    // which could be much too small if || r ||_2 is much larger than beta 
    // and sqrt(2).
    Real scale = Max( 2*Abs(beta), MaxNorm(dLoc) );
    SafeScale( Real(1), scale, d );
    SafeScale( Real(1), scale, beta );

    // Now that the problem is rescaled, our deflation tolerance simplifies to
    //
    //   tol = deflationFudge eps max( || d ||_max, 2*|beta| )
    //       = deflationFudge eps.
    //
    // Cf. LAPACK's {s,d}lasd2 [CITATION] for this tolerance.
    const Real eps = limits::Epsilon<Real>();
    const Real deflationTol = dcCtrl.deflationFudge*eps;

    // Get a full copy of the last row of Q0 and the first row of Q1.
    const Int lastRowOfQ0 = ( ctrl.wantEigVecs ? n0-1 : 0 );
    DistMatrix<Real,STAR,STAR> q0Last( Q0(IR(lastRowOfQ0),ALL) ), 
      q1First( Q1(IR(0),ALL) );
    const auto& q0LastLoc = q0Last.LockedMatrix();
    const auto& q1FirstLoc = q1First.LockedMatrix();

    // Put Q into [VC,STAR] distribution for Givens applications
    DistMatrix<Real,VC,STAR> Q_VC_STAR(Q);
    auto& Q_VC_STAR_Loc = Q_VC_STAR.Matrix();

    Matrix<Real> z(n,1);
    Matrix<Int> columnTypes(n,1);
    const Real betaSgn = Sgn( beta, false );
    const Real sqrtTwo = Sqrt( Real(2) );
    z(0) = betaSgn*q0LastLoc(0,n0-1) / sqrtTwo;
    columnTypes(0) = DENSE_COLUMN;
    for( Int j=0; j<n0-1; ++j )
    {
        z(j+1) = betaSgn*q0LastLoc(0,j) / sqrtTwo;
        columnTypes(j+1) = COLUMN_NONZERO_IN_FIRST_BLOCK;
    }
    for( Int j=0; j<n1; ++j )
    {
        z(j+n0) = q1FirstLoc(0,j) / sqrtTwo;
        columnTypes(j+n0) = COLUMN_NONZERO_IN_SECOND_BLOCK;
    }

    Permutation combineSortPerm;
    SortingPermutation( dLoc, combineSortPerm, ASCENDING );
    combineSortPerm.PermuteRows( dLoc );
    combineSortPerm.PermuteRows( z );
    combineSortPerm.PermuteRows( columnTypes );

    auto combinedToOrig = [&]( const Int& combinedIndex )
      {
          const Int preCombined = combineSortPerm.Preimage( combinedIndex );
          if( preCombined < n0 )
              // Undo the cyclic shift [0,n0) |-> [1,n0+1) mod n0 which
              // pushed the removed row into the first position.
              return Mod( preCombined-1, n0 );
          else
              return preCombined;
      };

    Permutation deflationPerm;
    deflationPerm.MakeIdentity( n );
    deflationPerm.MakeArbitrary();
    // Since we do not yet know how many undeflated entries there will be, we
    // must use the no-deflation case as our storage upper bound.
    Matrix<Real> dUndeflated(n,1), zUndeflated(n,1);
    dUndeflated(0) = 0;
    zUndeflated(0) = z(0);

    // Deflate all (off-diagonal) update entries sufficiently close to zero
    Int numDeflated = 0;
    Int numUndeflated = 0;
    // We will keep track of the last column that we encountered that was not
    // initially deflatable (but that could be deflated later due to close
    // diagonal entries if another undeflatable column is not encountered
    // first).
    Int revivalCandidate = n;
    for( Int j=0; j<n; ++j )
    {
        if( Abs(2*beta*z(j)) <= deflationTol )
        {
            // We can deflate due to the r component being sufficiently small
            const Int deflationDest = (n-1) - numDeflated;
            deflationPerm.SetImage( j, deflationDest );
            if( ctrl.progress && amRoot )
                Output
                ("Deflating via p(",j,")=",deflationDest,
                 " because |2*beta*z(",j,")|=|",2*beta*z(j),"| <= ",
                 deflationTol);
            columnTypes(j) = DEFLATED_COLUMN;
            ++numDeflated;
            if( amRoot )
            {
                ++secularInfo.numDeflations;
                ++secularInfo.numSmallUpdateDeflations;
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
    for( Int j=deflationRestart; j<n; ++j )
    {
        if( Abs(2*beta*z(j)) <= deflationTol )
        {
            const Int deflationDest = (n-1) - numDeflated;
            deflationPerm.SetImage( j, deflationDest );
            if( ctrl.progress && amRoot )
                Output
                ("Deflating via p(",j,")=",deflationDest,
                 " because |2*beta*z(",j,")|=|",2*beta*z(j),"| <= ",
                 deflationTol);
            columnTypes(j) = DEFLATED_COLUMN;
            ++numDeflated;
            if( amRoot )
            {
                ++secularInfo.numDeflations;
                ++secularInfo.numSmallUpdateDeflations;
            }
            continue;
        }
        const Real gamma = SafeNorm( z(j), z(revivalCandidate) );
        const Real c = z(j) / gamma;
        const Real s = z(revivalCandidate) / gamma;
        const Real offDiagNew = c*s*(dLoc(j)-dLoc(revivalCandidate));
        if( Abs(offDiagNew) <= deflationTol )
        {
            // Deflate the previously undeflatable index by rotating
            // z(revivalCandidate) into z(j) (Cf. the discussion
            // surrounding Eq. (4.4) of Gu/Eisenstat's TR [CITATION]).
            //
            // In particular, we want
            //
            //   | z(j), z(revivalCandidate) | | c -s | = | gamma, 0 |,
            //                                 | s  c |
            //
            // where gamma = || z(revivalCandidate); z(j) ||_2. Putting 
            //
            //   c = z(j)                / gamma,
            //   s = z(revivalCandidate) / gamma,
            //
            // implies
            //
            //   |  c,  s | |        z(j)         | = | gamma |.
            //   | -s,  c | | z(revivalCandidate) |   |   0   |
            //
            z(j) = gamma;
            z(revivalCandidate) = 0;

            // Apply | c -s | to both sides of d
            //       | s  c |
            const Real deltaDeflate =
              dLoc(revivalCandidate)*(c*c) + dLoc(j)*(s*s);
            dLoc(j) = dLoc(j)*(c*c) + dLoc(revivalCandidate)*(s*s);
            dLoc(revivalCandidate) = deltaDeflate;

            // Apply | c -s | from the right to Q
            //       | s  c |
            //
            const Int revivalOrig = combinedToOrig( revivalCandidate );
            const Int jOrig = combinedToOrig( j );
            // TODO(poulson): Exploit the nonzero structure of Q?
            blas::Rot
            ( Q_VC_STAR_Loc.Height(), Q_VC_STAR_Loc.Buffer(0,jOrig), 1,
              Q_VC_STAR_Loc.Buffer(0,revivalOrig), 1, c, s );

            const Int deflationDest = (n-1) - numDeflated;
            deflationPerm.SetImage( revivalCandidate, deflationDest );
            if( ctrl.progress )
                Output
                ("Deflating via p(",revivalCandidate,")=",
                 deflationDest," because |c*s*(d(",j,")-d(",revivalCandidate,
                 "))|=",offDiagNew," <= ",deflationTol);
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
        dUndeflated(numUndeflated) = dLoc(revivalCandidate);
        zUndeflated(numUndeflated) = z(revivalCandidate);
        deflationPerm.SetImage( revivalCandidate, numUndeflated );
        if( ctrl.progress )
            Output
            ("Could not deflate with j=",j," and revivalCandidate=",
             revivalCandidate,", so p(",revivalCandidate,")=",
             numUndeflated);
        ++numUndeflated;
        revivalCandidate = j;
    }
    if( revivalCandidate < n )
    {
        // Give up on the revival candidate
        dUndeflated(numUndeflated) = dLoc(revivalCandidate);
        zUndeflated(numUndeflated) = z(revivalCandidate);
        deflationPerm.SetImage( revivalCandidate, numUndeflated );
        if( ctrl.progress )
            Output
            ("Final iteration, so p(",revivalCandidate,")=",numUndeflated);
        ++numUndeflated;
    }
    // Now shrink dUndeflated and zUndeflated down to their proper size
    dUndeflated.Resize( numUndeflated, 1 );
    zUndeflated.Resize( numUndeflated, 1 );

    // Count the number of columns of Q with each nonzero pattern
    std::vector<Int> packingCounts( NUM_DC_COMBINED_COLUMN_TYPES, 0 );
    for( Int j=0; j<n; ++j )
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
        if( ctrl.progress && amRoot )
            Output("packingCounts[",columnType,"]=",packingCounts[columnType]);
    }

    // Set up the index ranges of the three packed column subsets
    const Range<Int> packingInd0(packingOffsets[0],packingOffsets[1]),
                     packingInd1(packingOffsets[1],packingOffsets[2]),
                     packingInd2(packingOffsets[2],packingOffsets[3]);

    Matrix<Real> dPacked;
    DistMatrix<Real,VC,STAR> QPacked(g);
    dPacked.Resize( n, 1 );
    if( ctrl.wantEigVecs )
        QPacked.Resize( n, n );
    else
        QPacked.Resize( 2, n );
    auto& QPackedLoc = QPacked.Matrix();
    Permutation packingPerm;
    packingPerm.MakeIdentity( n );
    for( Int j=0; j<n; ++j )
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
        // TODO(poulson): Exploit the nonzero structure of Q?
        blas::Copy
        ( Q_VC_STAR_Loc.Height(), Q_VC_STAR_Loc.Buffer(0,jOrig), 1,
          QPackedLoc.Buffer(0,packingDest), 1 );
    }

    // Put the deflated columns in their final destination and shrink QPacked
    // down to its final size.
    //
    // TODO(poulson): Exploit the nonzero structure of Q?
    if( numDeflated > 0 )
    {
        blas::Copy
        ( numDeflated, &dPacked(numUndeflated), 1, &dLoc(numUndeflated), 1 );
        auto QDeflated = Q(ALL,IR(numUndeflated,END));
        auto QPackedDeflated = QPacked(ALL,IR(numUndeflated,END));
        QDeflated = QPackedDeflated;
    }
    if( ctrl.wantEigVecs )
        QPacked.Resize( n, numUndeflated );
    else
        QPacked.Resize( 2, numUndeflated );

    // Now compute the updated eigenvectors using QPacked
    // ==================================================
    auto undeflatedInd = IR(0,numUndeflated);
    const Real zUndeflatedNorm = FrobeniusNorm( zUndeflated );
    zUndeflated *= Real(1) / zUndeflatedNorm;
    const Real rho = 2*Abs(beta)*zUndeflatedNorm*zUndeflatedNorm;

    if( ctrl.progress && amRoot )
        Output("Solving secular equation and correcting update vector");
    Matrix<Real> rCorrected;
    Ones( rCorrected, numUndeflated, 1 );

    // Ensure that there is sufficient space for storing the needed eigenvectors
    // from the undeflated secular equation. Notice that we *always* need to
    // compute the eigenvectors of the undeflated secular equation.
    DistMatrix<Real,VR,STAR> dSecular(numUndeflated,1,g);
    DistMatrix<Real,STAR,VR> QSecular(g);
    QSecular.Resize( numUndeflated, numUndeflated );
    auto& dSecularLoc = dSecular.Matrix();
    auto& QSecularLoc = QSecular.Matrix();

    const Int numUndeflatedLoc = QSecularLoc.Width();
    for( Int jLoc=0; jLoc<numUndeflatedLoc; ++jLoc )
    {
        const Int j = QSecular.GlobalCol(jLoc);
        auto minusShift = QSecularLoc( ALL, IR(jLoc) );

        auto valueInfo =
          SecularEigenvalue
          ( j, dUndeflated, rho, zUndeflated, dSecularLoc(jLoc),
            minusShift, dcCtrl.secularCtrl );
        if( ctrl.progress && amRoot )
            Output("Secular eigenvalue ",j," is ",dSecularLoc(jLoc));

        // We will sum these across all of the processors at the top-level
        secularInfo.numIterations += valueInfo.numIterations;
        secularInfo.numAlternations += valueInfo.numAlternations;
        secularInfo.numCubicIterations += valueInfo.numCubicIterations;
        secularInfo.numCubicFailures += valueInfo.numCubicFailures;

        rCorrected(j) *= minusShift(j);
        for( Int k=0; k<numUndeflated; ++k )
        {
            if( k == j )
                continue;
            rCorrected(k) *= minusShift(k) / (dUndeflated(j)-dUndeflated(k));
        }
    }
    AllReduce( rCorrected, g.VRComm(), mpi::PROD );
    for( Int j=0; j<numUndeflated; ++j )
        rCorrected(j) = Sgn(zUndeflated(j),false) * Sqrt(Abs(rCorrected(j)));
    // Get a full copy of the undeflated eigenvalues now
    {
        auto wUndeflated = d( IR(0,numUndeflated), ALL );
        wUndeflated = dSecular;
    }

    // Compute the unnormalized eigenvectors.
    if( ctrl.progress && amRoot )
        Output("Computing unnormalized eigenvectors");
    for( Int jLoc=0; jLoc<numUndeflatedLoc; ++jLoc )
    {
        auto q = QSecularLoc(ALL,IR(jLoc));
        for( Int i=0; i<numUndeflated; ++i )
            q(i) = rCorrected(i) / q(i);
    }


    // Form the normalized eigenvectors with the rows permuted by
    // the inverse of the packing permutation in U. This allows the product
    // of QPacked with U to be equal to the unpacked Q times the eigenvectors
    // from the secular equation.
    if( ctrl.progress && amRoot )
        Output("Forming undeflated right singular vectors");
    DistMatrix<Real,STAR,VR> U(g);
    U.Resize( numUndeflated, numUndeflated );
    auto& ULoc = U.Matrix();
    for( Int jLoc=0; jLoc<numUndeflatedLoc; ++jLoc )
    {
        auto q = QSecularLoc(ALL,IR(jLoc));
        auto u = ULoc(ALL,IR(jLoc));
        const Real qFrob = FrobeniusNorm( q );
        for( Int i=0; i<numUndeflated; ++i )
            u(i) = q(packingPerm.Preimage(i)) / qFrob;
    }
    // Overwrite the first 'numUndeflated' columns of Q with the updated
    // eigenvectors by exploiting the partitioning of Z = QPacked as
    //
    //   Z = | Z_{0,0} |    0    | Z_{0,2} |,
    //       |---------|---------|---------|
    //       |    0    | Z_{1,1} | Z_{1,2} |
    //
    // where the first and second block rows have heights n0 and n1. The block
    // columns respectively have widths packingCounts[0], packingCounts[1], and
    // packingCounts[2].
    //
    // Conformally partitioning U, we have
    //
    //   Z U = Z_{:,2} U2 + | Z_{0,0} U_0 |,
    //                      |-------------|
    //                      | Z_{1,1} U_1 |
    //
    if( ctrl.progress && amRoot )
        Output("Overwriting eigenvectors");
    auto QUndeflated = Q( ALL, undeflatedInd );
    if( ctrl.wantEigVecs )
    {
        if( dcCtrl.exploitStructure )
        {
            auto Z2 = QPacked( ALL, packingInd2 );
            auto U2 = U( packingInd2, ALL );
            Gemm( NORMAL, NORMAL, Real(1), Z2, U2, QUndeflated );

            // Finish updating the first block row
            auto Q0Undeflated = QUndeflated( IR(0,n0), ALL );
            auto Z00 = QPacked( IR(0,n0), packingInd0 );
            auto U0 = U( packingInd0, ALL );
            Gemm( NORMAL, NORMAL, Real(1), Z00, U0, Real(1), Q0Undeflated );

            // Finish updating the second block row
            auto Q1Undeflated = QUndeflated( IR(n0,n), ALL );
            auto Z11 = QPacked( IR(n0,n), packingInd1 );
            auto U1 = U( packingInd1, ALL );
            Gemm( NORMAL, NORMAL, Real(1), Z11, U1, Real(1), Q1Undeflated );
        }
        else
        {
            Gemm( NORMAL, NORMAL, Real(1), QPacked, U, QUndeflated );
        }
    }
    else
    {
        if( dcCtrl.exploitStructure )
        {
            auto Z2 = QPacked( ALL, packingInd2 );
            auto U2 = U( packingInd2, ALL );
            Gemm( NORMAL, NORMAL, Real(1), Z2, U2, QUndeflated );

            // Finish updating the first block row
            auto Q0Undeflated = QUndeflated( IR(0), ALL );
            auto Z00 = QPacked( IR(0), packingInd0 );
            auto U0 = U( packingInd0, ALL );
            Gemm( NORMAL, NORMAL, Real(1), Z00, U0, Real(1), Q0Undeflated );

            // Finish updating the second block row
            auto Q1Undeflated = QUndeflated( IR(1), ALL );
            auto Z11 = QPacked( IR(1), packingInd1 );
            auto U1 = U( packingInd1, ALL );
            Gemm( NORMAL, NORMAL, Real(1), Z11, U1, Real(1), Q1Undeflated );
        }
        else
        {
            Gemm( NORMAL, NORMAL, Real(1), QPacked, U, QUndeflated );
        }
    }

    // Rescale the eigenvalues
    SafeScale( scale, Real(1), d );

    return info;
}

template<typename Real>
DCInfo
DivideAndConquer
(       Matrix<Real>& mainDiag,
        Matrix<Real>& superDiag,
        Matrix<Real>& w,
        Matrix<Real>& Q,
  const HermitianTridiagEigCtrl<Real>& ctrl )
{
    DEBUG_CSE
    const Int n = mainDiag.Height();
    const auto& dcCtrl = ctrl.dcCtrl;

    DCInfo info;
    auto& secularInfo = info.secularInfo;

    if( n <= Max(dcCtrl.cutoff,3) )
    {
        auto ctrlMod( ctrl );
        ctrlMod.alg = HERM_TRIDIAG_EIG_QR;
        if( ctrl.wantEigVecs )
        {
            ctrlMod.accumulateEigVecs = false;
            HermitianTridiagEig( mainDiag, superDiag, w, Q, ctrlMod );
        }
        else
        {
            ctrlMod.wantEigVecs = true;
            ctrlMod.accumulateEigVecs = true;
            Zeros( Q, 2, n );
            Q(0,0) = 1;
            Q(1,n-1) = 1;
            HermitianTridiagEig( mainDiag, superDiag, w, Q, ctrlMod );
        }
        return info;
    }

    // TODO(poulson): A more intelligent split point.
    const Int split = (n/2) + 1;
    const Real& beta = superDiag(split-1);

    auto mainDiag0 = mainDiag( IR(0,split), ALL );
    auto superDiag0 = superDiag( IR(0,split-1), ALL );
    mainDiag0(split-1) -= Abs(beta);

    auto mainDiag1 = mainDiag( IR(split,END), ALL );
    auto superDiag1 = superDiag( IR(split,END), ALL );
    mainDiag1(0) -= Abs(beta);

    Matrix<Real> Q0, Q1;
    if( ctrl.wantEigVecs )
    {
        Zeros( Q, n, n );
        View( Q0, Q, IR(0,split), IR(0,split) );
        View( Q1, Q, IR(split,END), IR(split,END) );
    }
    else
    {
        Zeros( Q0, 2, split );
        Zeros( Q1, 2, n-split );
    }

    Matrix<Real> w0;
    auto info0 = DivideAndConquer( mainDiag0, superDiag0, w0, Q0, ctrl );

    Matrix<Real> w1;
    auto info1 = DivideAndConquer( mainDiag1, superDiag1, w1, Q1, ctrl );

    if( !ctrl.wantEigVecs )
    {
        // We must manually pack the last row of Q0 and the first row of Q1
        Zeros( Q, 2, n );
        auto Q0Last = Q( IR(0), IR(0,split) );
        auto Q1First = Q( IR(1), IR(split,n) );
        Q0Last = Q0( IR(1), ALL );
        Q1First = Q1( IR(0), ALL );
    }
    info = Merge( beta, w0, w1, w, Q, ctrl );

    secularInfo.numIterations += info0.secularInfo.numIterations;
    secularInfo.numAlternations += info0.secularInfo.numAlternations;
    secularInfo.numCubicIterations += info0.secularInfo.numCubicIterations;
    secularInfo.numCubicFailures += info0.secularInfo.numCubicFailures;
    secularInfo.numDeflations += info0.secularInfo.numDeflations;
    secularInfo.numCloseDiagonalDeflations +=
      info0.secularInfo.numCloseDiagonalDeflations;
    secularInfo.numSmallUpdateDeflations +=
      info0.secularInfo.numSmallUpdateDeflations;

    secularInfo.numIterations += info1.secularInfo.numIterations;
    secularInfo.numAlternations += info1.secularInfo.numAlternations;
    secularInfo.numCubicIterations += info1.secularInfo.numCubicIterations;
    secularInfo.numCubicFailures += info1.secularInfo.numCubicFailures;
    secularInfo.numDeflations += info1.secularInfo.numDeflations;
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
(       Matrix<Real>& mainDiag,
        Matrix<Real>& superDiag,
        DistMatrix<Real,STAR,STAR>& w,
        DistMatrix<Real>& Q,
  const HermitianTridiagEigCtrl<Real>& ctrl,
  bool topLevel=true )
{
    DEBUG_CSE
    const Grid& grid = Q.Grid();
    const Int n = mainDiag.Height();
    const auto& dcCtrl = ctrl.dcCtrl;

    DCInfo info;

    if( n <= Max(dcCtrl.cutoff,3) )
    {
        // Run the problem redundantly locally
        auto ctrlMod( ctrl );
        ctrlMod.alg = HERM_TRIDIAG_EIG_QR;
        Matrix<Real> wLoc, QLoc;
        if( ctrl.wantEigVecs )
        {

            ctrlMod.accumulateEigVecs = false;
            HermitianTridiagEig( mainDiag, superDiag, wLoc, QLoc, ctrlMod );
        }
        else
        {
            ctrlMod.wantEigVecs = true;
            ctrlMod.accumulateEigVecs = true;
            Zeros( QLoc, 2, n );
            QLoc(0,0) = 1;
            QLoc(1,n-1) = 1;
            HermitianTridiagEig( mainDiag, superDiag, wLoc, QLoc, ctrlMod );
        }

        w.Resize( n, 1 );
        for( Int i=0; i<n; ++i )
            w.SetLocal( i, 0, wLoc(i) );

        if( ctrl.wantEigVecs )
            Q.Resize( n, n );
        else
            Q.Resize( 2, n );
        const Int QLocHeight = Q.LocalHeight();
        const Int QLocWidth = Q.LocalWidth();
        for( Int iLoc=0; iLoc<QLocHeight; ++iLoc )
        {
            const Int i = Q.GlobalRow(iLoc);
            for( Int jLoc=0; jLoc<QLocWidth; ++jLoc )
            {
                const Int j = Q.GlobalCol(jLoc);
                Q.SetLocal( iLoc, jLoc, QLoc(i,j) );
            }
        }

        return info;
    }
    else if( grid.Size() == 1 )
    {
        // TODO(poulson): Avoid unnecessary copies
        Matrix<Real> wLoc, QLoc;
        info = DivideAndConquer( mainDiag, superDiag, wLoc, QLoc, ctrl );

        w.Resize( n, 1 );
        w.Matrix() = wLoc;

        Q.Resize( QLoc.Height(), QLoc.Width() );
        Q.Matrix() = QLoc;

        return info;
    }

    // TODO(poulson): A more intelligent split point.
    const Int split = (n/2) + 1;
    const Grid *leftGrid, *rightGrid;
    SplitGrid( split, n-split, grid, leftGrid, rightGrid );
    const Real& beta = superDiag(split-1);

    auto mainDiag0 = mainDiag( IR(0,split), ALL );
    auto superDiag0 = superDiag( IR(0,split-1), ALL );
    mainDiag0(split-1) -= Abs(beta);

    auto mainDiag1 = mainDiag( IR(split,END), ALL );
    auto superDiag1 = superDiag( IR(split,END), ALL );
    mainDiag1(0) -= Abs(beta);

    // Split Q between the subtrees
    DistMatrix<Real> Q0Sub(*leftGrid), Q1Sub(*rightGrid);
    if( ctrl.wantEigVecs )
    {
        Zeros( Q0Sub, split, split );
        Zeros( Q1Sub, n-split, n-split );
    }
    else
    {
        Zeros( Q0Sub, 2, split );
        Zeros( Q1Sub, 2, n-split );
    }

    DistMatrix<Real,STAR,STAR> w0Sub(*leftGrid), w1Sub(*rightGrid);
    DCInfo info0, info1;
    if( w0Sub.Participating() )
    {
        info0 =
          DivideAndConquer( mainDiag0, superDiag0, w0Sub, Q0Sub, ctrl, false );
    }
    if( w1Sub.Participating() )
    {
        info1 =
          DivideAndConquer( mainDiag1, superDiag1, w1Sub, Q1Sub, ctrl, false );
    }

    const bool includeViewers = true;
    w0Sub.MakeConsistent( includeViewers );
    Q0Sub.MakeConsistent( includeViewers );
    w1Sub.MakeConsistent( includeViewers );
    Q1Sub.MakeConsistent( includeViewers );

    DistMatrix<Real,STAR,STAR> w0(grid), w1(grid);
    w0 = w0Sub;
    w1 = w1Sub;
    if( ctrl.wantEigVecs )
    {
        Zeros( Q, n, n );
        auto Q0 = Q(IR(0,split),IR(0,split));
        auto Q1 = Q(IR(split,n),IR(split,n));
        Q0 = Q0Sub;
        Q1 = Q1Sub;
    }
    else
    {
        // We must manually pack the last row of Q0 and the first row of Q1
        Zeros( Q, 2, n );
        auto Q0Last = Q( IR(0), IR(0,split) );
        auto Q1First = Q( IR(1), IR(split,n) );
        Q0Last = Q0Sub( IR(1), ALL );
        Q1First = Q1Sub( IR(0), ALL );
    }
    info = Merge( beta, w0, w1, w, Q, ctrl );

    auto& secularInfo = info.secularInfo;
    if( w0Sub.Participating() )
    {
        secularInfo.numIterations += info0.secularInfo.numIterations;
        secularInfo.numAlternations += info0.secularInfo.numAlternations;
        secularInfo.numCubicIterations += info0.secularInfo.numCubicIterations;
        secularInfo.numCubicFailures += info0.secularInfo.numCubicFailures;
    }
    if( w1Sub.Participating() )
    {
        secularInfo.numIterations += info1.secularInfo.numIterations;
        secularInfo.numAlternations += info1.secularInfo.numAlternations;
        secularInfo.numCubicIterations += info1.secularInfo.numCubicIterations;
        secularInfo.numCubicFailures += info1.secularInfo.numCubicFailures;
    }

    if( leftGrid->Rank() == 0 )
    {
        secularInfo.numDeflations += info0.secularInfo.numDeflations;
        secularInfo.numCloseDiagonalDeflations +=
          info0.secularInfo.numCloseDiagonalDeflations;
        secularInfo.numSmallUpdateDeflations +=
          info0.secularInfo.numSmallUpdateDeflations;
    }
    if( rightGrid->Rank() == 0 )
    {
        secularInfo.numDeflations += info1.secularInfo.numDeflations;
        secularInfo.numCloseDiagonalDeflations +=
          info1.secularInfo.numCloseDiagonalDeflations;
        secularInfo.numSmallUpdateDeflations +=
          info1.secularInfo.numSmallUpdateDeflations;
    }

    if( topLevel )
    {
        // Sum all of the secular iteration/deflation information
        Matrix<Int> counts(6,1);
        counts(0) = secularInfo.numIterations;
        counts(1) = secularInfo.numAlternations;
        counts(2) = secularInfo.numCubicIterations;
        counts(3) = secularInfo.numCubicFailures;
        counts(4) = secularInfo.numCloseDiagonalDeflations;
        counts(5) = secularInfo.numSmallUpdateDeflations;

        AllReduce( counts, grid.Comm() );

        secularInfo.numIterations = counts(0);
        secularInfo.numAlternations = counts(1);
        secularInfo.numCubicIterations = counts(2);
        secularInfo.numCubicFailures = counts(3);
        secularInfo.numCloseDiagonalDeflations = counts(4);
        secularInfo.numSmallUpdateDeflations = counts(5);
        secularInfo.numDeflations = counts(4) + counts(5);
    }

    return info;
}

} // namespace herm_tridiag_eig
} // namespace El

#endif // ifndef EL_HERM_TRIDIAG_EIG_DC_HPP
