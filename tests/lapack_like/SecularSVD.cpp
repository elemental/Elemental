/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
using namespace El;

extern "C" {

void EL_LAPACK(slasd4)
( const BlasInt* n,
  const BlasInt* i,
  const float* d,
  const float* z,
  float* dMinusShift,
  const float* rho,
  float* sigma,
  float* dPlusShift,
  BlasInt* info );

void EL_LAPACK(dlasd4)
( const BlasInt* n,
  const BlasInt* i,
  const double* d,
  const double* z,
  double* dMinusShift,
  const double* rho,
  double* sigma,
  double* dPlusShift,
  BlasInt* info );

} // extern "C"

template<typename Real>
void TestLAPACK
( const Matrix<Real>& d, const Real& rho, const Matrix<Real>& z );

void TestLAPACK
( const Matrix<float>& d, const float& rho, const Matrix<float>& z )
{
    typedef float Real;
    const Int n = d.Height();
    Timer timer;
    Matrix<Real> wLAPACK(n,1), dPlusShift(n,1), dMinusShift(n,1);
    BlasInt nBLAS = n;
    BlasInt infoLAPACK;
    timer.Start();
    for( Int i=0; i<n; ++i )
    {
        Real sigmaLAPACK;
        const BlasInt ip1 = i+1;
        EL_LAPACK(slasd4)
        ( &nBLAS, &ip1, d.LockedBuffer(), z.LockedBuffer(),
          dMinusShift.Buffer(), &rho, &sigmaLAPACK,
          dPlusShift.Buffer(), &infoLAPACK );
        wLAPACK(i) = sigmaLAPACK;
    }
    const Real lapackTime = timer.Stop(); 
    Output("LAPACK: ",lapackTime," seconds");
    Output("");
}

void TestLAPACK
( const Matrix<double>& d, const double& rho, const Matrix<double>& z )
{
    typedef double Real;
    const Int n = d.Height();
    Timer timer;
    Matrix<Real> wLAPACK(n,1), dPlusShift(n,1), dMinusShift(n,1);
    BlasInt nBLAS = n;
    BlasInt infoLAPACK;
    timer.Start();
    for( Int i=0; i<n; ++i )
    {
        Real sigmaLAPACK;
        const BlasInt ip1 = i+1;
        EL_LAPACK(dlasd4)
        ( &nBLAS, &ip1, d.LockedBuffer(), z.LockedBuffer(),
          dMinusShift.Buffer(), &rho, &sigmaLAPACK,
          dPlusShift.Buffer(), &infoLAPACK );
        wLAPACK(i) = sigmaLAPACK;
    }
    const Real lapackTime = timer.Stop(); 
    Output("LAPACK: ",lapackTime," seconds");
    Output("");
}

template<typename Real>
void GenerateData
( Int n, Matrix<Real>& d, Real& rho, Matrix<Real>& z, bool print )
{
    // Implicitly form a matrix
    //
    //   M = | sqrt(rho)*z(0), sqrt(rho)*z(1), ..., sqrt(rho)*z(n-1) |
    //       |                      d(1),                            |
    //       |                                 .                     |
    //       |                                                d(n-1) |
    //
    // where 0 = d(0) <= d(1) <= d(2) <= ... <= d(n-1).
    //
    Uniform( d, n, 1, Real(2), Real(2) );
    Sort( d );
    d(0) = 0;
    Gaussian( z, n, 1 );
    z *= Real(1) / FrobeniusNorm( z );
    rho = SampleUniform( Real(1), Real(1)/Real(2) );
    if( print )
    {
        Print( d, "d" );
        Output( "rho=", rho );
        Print( z, "z" );
    }
}

struct SecularDeflationInfo
{
    Int numDeflations=0;

    Int numSmallDiagonalDeflations=0;
    Int numCloseDiagonalDeflations=0;
    Int numSmallUpdateDeflations=0;
};
// TODO(poulson): Combine with SecularSingularValueInfo structure

// Cf. Section 4 of Gu and Eisenstat's "A Divide-and-Conquer Algorithm for the
// Bidiagonal SVD" [CITATION] and LAPACK's {s,d}lasd2 [CITATION].
//
// We begin with the decomposition
//
// B = | U_0, 0,  0  | |     diag(s_0),              0       | | V_0, 0   |^T,
//     |  0,  1,  0  | | alpha e_{m_0}^T V_0, beta*e_0^T*V_1 | |   0, V_1 |
//     |  0,  0, U_1 | |         0,              diag(s_1)   |
//
// where U_0 is m_0 x m_0, U_1 is m_1 x m_1, V_0 is (m0+1) x (m0+1), and V_1 is
// either m1 x m1 or (m1+1) x (m1+1). Thus, putting m = m_0 + 1 + m_1, B is
// either m x m or m x (m+1). On entry, U and V should be filled with their
// above depictions.
//
// We operationalize Gu and Eisenstat's [CITATION] deflation-tracking
// mechanism by initializing the tags for the nonzero structure of the
// columns of the singular vectors:
//
//   0: nonzero in first block
//   1: nonzero in second block
//   2: nonzero in both blocks
//   3: deflated
//
// Cf. LAPACK's {s,d}lasd2 [CITATION] for this mechanism. Note that LAPACK 
// currently ignores deflations of the form |d(0)-d(j)| <= deflationTol, 
// which results in the first column of U potentially becoming dense. We
// do not ignore such deflations and always mark the first column of U
// as dense for the sake of simplicity.
//
enum SecularCombinedColumnType {
  COLUMN_NONZERO_IN_FIRST_BLOCK = 0,
  COLUMN_NONZERO_IN_SECOND_BLOCK = 1,
  DENSE_COLUMN = 2,
  DEFLATED_COLUMN = 3
};
const Int NUM_SECULAR_COMBINED_COLUMN_TYPES = 4;

// The following is analogous to LAPACK's {s,d}lasd{1,2,3} [CITATION] but does
// not accept initial sorting permutations for s0 and s1, nor does it enforce 
// any ordering on the resulting singular values. Several bugs in said LAPACK
// routines were found and reported to
// https://github.com/Reference-LAPACK/lapack/issues/34.
template<typename Real>
SecularDeflationInfo SecularCombine
( const Real& alpha,
  // The right entry in the removed middle row of the bidiagonal matrix
  const Real& beta,
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
  const SecularSingularValueCtrl<Real>& ctrl )
{
    DEBUG_CSE
    SecularDeflationInfo info;
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
    if( ctrl.progress )
        Output("m=",m,", n=",n,", m0=",m0,", n0=",n0,", m1=",m1,", n1=",n1);

    // TODO(poulson): Turn this back on
    const bool exploitStructure = false;

    // TODO(poulson): Add scaling

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
    // will rotate rhoExtra into r(0) towards the end of the routine.
    //
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
    // We form r(m) and rotate it into r(0) after the deflation tolerance is
    // available.

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

    // Now that we can read off || diag(d) ||_2 = || d ||_max = d(m-1), we 
    // calculate the deflation tolerance,
    //
    //   tol = 8 eps max( || d ||_max, |alpha|, |beta| ).
    //
    // Cf. LAPACK's {s,d}lasd2 [CITATION] for this tolerance.
    const Real eps = limits::Epsilon<Real>();
    const Real deflationTolScale = Max( d(m-1), Max( Abs(alpha), Abs(beta) ) );
    const Real deflationTol = 8*eps*deflationTolScale;

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

    Permutation deflationPerm;
    deflationPerm.MakeIdentity( m );
    deflationPerm.MakeArbitrary();
    // Since we do not yet know how many undeflated entries there will be, we
    // must use the no-deflation case as our storage upper bound.
    Matrix<Real> dUndeflated(m,1), rUndeflated(m,1);
    dUndeflated(0) = 0;
    rUndeflated(0) = r(0);

    // Deflate all (off-diagonal) update entries sufficiently close to zero
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
            deflationPerm.SetImage( j, (m-1)-info.numDeflations );
            if( ctrl.progress )
                Output
                ("Deflating via p(",j,")=",(m-1)-info.numDeflations,
                 " because |r(",j,")|=|",r(j),"| <= ",deflationTol);
            columnTypes(j) = DEFLATED_COLUMN;
            ++info.numDeflations;
            ++info.numSmallUpdateDeflations; 
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

            // Apply | c -s | from the right to U and V
            //       | s  c |
            //
            // We are mixing nonzero structures in the first column of U,
            // so we might as well always treat the first column as dense.
            //
            // TODO(poulson): Exploit the nonzero structure of U and V?
            const Int jOrig = combinedToOrig( j );
            blas::Rot( m, &U(0,m0), 1, &U(0,jOrig), 1, c, s );
            blas::Rot( n, &V(0,m0), 1, &V(0,jOrig), 1, c, s );

            deflationPerm.SetImage( j, (m-1)-info.numDeflations ); 
            if( ctrl.progress )
                Output
                ("Deflating via p(",j,")=",(m-1)-info.numDeflations,
                 " because d(",j,")=",d(j)," <= ",deflationTol);

            columnTypes(j) = DEFLATED_COLUMN;

            ++info.numDeflations;
            ++info.numCloseDiagonalDeflations;
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
            deflationPerm.SetImage( j, (m-1)-info.numDeflations );
            if( ctrl.progress )
                Output
                ("Deflating via p(",j,")=",(m-1)-info.numDeflations,
                 " because |r(",j,")|=|",r(j),"| <= ",deflationTol);
            columnTypes(j) = DEFLATED_COLUMN;
            ++info.numDeflations;
            ++info.numSmallUpdateDeflations;
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

            deflationPerm.SetImage
            ( revivalCandidate, (m-1)-info.numDeflations ); 
            if( ctrl.progress )
                Output
                ("Deflating via p(",revivalCandidate,")=",
                 (m-1)-info.numDeflations," because d(",j,")=",d(j)," - d(",
                 revivalCandidate,")=",d(revivalCandidate)," <= ",deflationTol);

            if( columnTypes(revivalCandidate) != columnTypes(j) )
            {
                // We mixed top and bottom columns so the result is dense.
                columnTypes(j) = DENSE_COLUMN;
            }
            columnTypes(revivalCandidate) = DEFLATED_COLUMN;

            revivalCandidate = j;
            ++info.numDeflations;
            ++info.numCloseDiagonalDeflations;
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

    // While I am not aware of any documentation, LAPACK's {s,d}lasd2 [CITATION]
    // ensures that the first nonzero entry of d and dUndeflated are at least
    // deflationTol / 2. This would appear to be essentially impossible to 
    // activate given that the deflation condition would imply that d(1) would
    // have been deflated, and the fact that deflated values are stored in 
    // (roughly) descending order would imply that || d ||_2 ~= d(1), which
    // contradicts d(1) < deflationTol <= 8 eps || d ||_2 unless d(1) is vastly
    // different from || d ||_2. One possibility would perhaps be for a 3x3
    // matrix with an out-of-order deflation. Why it would be useful to enforce
    // a deflated singular value being at least roughly 4 eps || d ||_2 is 
    // beyond me.
    DEBUG_ONLY(
      if( d(1) < deflationTol/2 )
      {
          Output
          ("Exceptional d(1)=",d(1)," < deflationTol/2=",deflationTol/2,
           " case encountered");
          d(1) = deflationTol/2;
      }
    )

    // Count the number of columns of U with each nonzero pattern
    std::vector<Int> packingCounts( NUM_SECULAR_COMBINED_COLUMN_TYPES, 0 );
    for( Int j=0; j<m; ++j )
        ++packingCounts[columnTypes(j)];
    DEBUG_ONLY(
      if( packingCounts[DEFLATED_COLUMN] != info.numDeflations )
          LogicError
          ("Inconsistency between packingCounts[DEFLATED_COLUMN]=",
           packingCounts[DEFLATED_COLUMN],
           " and info.numDeflations=",info.numDeflations);
    )

    // Compute offsets for packing them
    std::vector<Int> packingOffsets( NUM_SECULAR_COMBINED_COLUMN_TYPES, 0 );
    Int totalPacked = 0;
    for( Int columnType=0; columnType<NUM_SECULAR_COMBINED_COLUMN_TYPES;
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
    if( info.numDeflations > 0 )
    {
        blas::Copy
        ( info.numDeflations, &dPacked(numUndeflated), 1,
          &d(numUndeflated), 1 );
        lapack::Copy
        ( 'A', m, info.numDeflations, &UPacked(0,numUndeflated), UPacked.LDim(),
          &U(0,numUndeflated), U.LDim() );
        lapack::Copy
        ( 'A', n, info.numDeflations, &VPacked(0,numUndeflated), VPacked.LDim(),
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
        auto info =
          SecularSingularValue
          ( j, dUndeflated, rho, rUndeflated, u, vScratch, ctrl );
        d(j) = info.singularValue;
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
    if( exploitStructure )
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
    if( exploitStructure )
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

    // TODO(poulson): Add scaling

    return info;
}

template<typename Real>
void TestBidiagonalCombine
( Int m,
  bool square,
  Int maxIter,  
  Int maxCubicIter,
  FlipOrClip negativeFix,
  bool progress,
  bool print )
{
    Output("Testing BidiagonalCombine(",square,") with ",TypeName<Real>());

    SecularSingularValueCtrl<Real> ctrl;
    ctrl.maxIterations = maxIter;
    ctrl.maxCubicIterations = maxCubicIter;
    ctrl.negativeFix = negativeFix;
    ctrl.progress = progress;

    const Int n = ( square ? m : m+1 );
    Matrix<Real> mainDiag, superDiag;
    Uniform( mainDiag, m, 1 );
    Uniform( superDiag, n-1, 1 );
    if( print )
    {
        Print( mainDiag, "mainDiag" );
        Print( superDiag, "superDiag" );
    }
    const Int split = m/2;

    const Real alpha = mainDiag(split);
    const Real beta = superDiag(split);
    if( print )
    {
        Output("alpha=",alpha,", beta=",beta);
    }

    auto mainDiag0 = mainDiag( IR(0,split), ALL );
    auto superDiag0 = superDiag( IR(0,split), ALL );
    if( print )
    {
        Print( mainDiag0, "mainDiag0" );
        Print( superDiag0, "superDiag0" );
    }

    auto mainDiag1 = mainDiag( IR(split+1,END), ALL );
    auto superDiag1 = superDiag( IR(split+1,END), ALL );
    if( print )
    {
        Print( mainDiag1, "mainDiag1" );
        Print( superDiag1, "superDiag1" );
    }

    BidiagSVDCtrl<Real> bidiagSVDCtrl;
    bidiagSVDCtrl.approach = FULL_SVD; // We need any null space as well
    bidiagSVDCtrl.progress = ctrl.progress;
   
    Matrix<Real> U, V;
    Identity( U, m, m );
    Identity( V, n, n );
    auto U0 = U( IR(0,split), IR(0,split) );
    auto U1 = U( IR(split+1,END), IR(split+1,END) );
    auto V0 = V( IR(0,split+1), IR(0,split+1) );
    auto V1 = V( IR(split+1,END), IR(split+1,END) );

    Matrix<Real> s0, s1;
    BidiagSVD( UPPER, mainDiag0, superDiag0, U0, s0, V0, bidiagSVDCtrl );
    BidiagSVD( UPPER, mainDiag1, superDiag1, U1, s1, V1, bidiagSVDCtrl );
    if( print )
    {
        Print( U0, "U0" );
        Print( s0, "s0" );
        Print( V0, "V0" );
        Print( U1, "U1" );
        Print( s1, "s1" );
        Print( V1, "V1" );
    }
   
    Matrix<Real> d;
    SecularCombine( alpha, beta, s0, s1, U, d, V, ctrl );
    if( print )
    {
        Print( U, "U" );
        Print( d, "d" );
        Print( V, "V" );
    }
}

template<typename Real>
void TestSecularHelper
( const Matrix<Real>& d,
  const Real& rho,
  const Matrix<Real>& z,
  Int maxIter,
  Int maxCubicIter,
  FlipOrClip negativeFix,
  bool progress,
  bool print,
  bool testFull )
{
    Output("Testing with ",TypeName<Real>());
    const Int n = d.Height();

    Timer timer;

    SecularSingularValueCtrl<Real> ctrl;
    ctrl.maxIterations = maxIter;
    ctrl.maxCubicIterations = maxCubicIter;
    ctrl.negativeFix = negativeFix;
    ctrl.progress = progress;

    Matrix<Real> s(n,1), wSecular(n,1);
    Int measMinIter=1e9, measMaxIter=0, measTotalIter=0,
        measMinCubicIter=1e9, measMaxCubicIter=0, measTotalCubicIter=0,
        measMinCubicFails=1e9, measMaxCubicFails=0, measTotalCubicFails=0;
    timer.Start();
    for( Int i=0; i<n; ++i )
    {
        auto info = SecularSingularValue( i, d, rho, z, ctrl );
        s(i) = info.singularValue;
        wSecular(i) = info.singularValue*info.singularValue;

        measMinIter = Min( measMinIter, info.numIterations );
        measMaxIter = Max( measMaxIter, info.numIterations );
        measTotalIter += info.numIterations;

        measMinCubicIter = Min( measMinCubicIter, info.numCubicIterations );
        measMaxCubicIter = Max( measMaxCubicIter, info.numCubicIterations );
        measTotalCubicIter += info.numCubicIterations;

        measMinCubicFails = Min( measMinCubicFails, info.numCubicFailures );
        measMaxCubicFails = Max( measMaxCubicFails, info.numCubicFailures );
        measTotalCubicFails += info.numCubicFailures;
    }
    const Real secularTime = timer.Stop();
    Output("Secular: ",secularTime," seconds");
    Output
    ("Iterations [min/max/total]: ",
     measMinIter,"/",measMaxIter,"/",measTotalIter);
    Output
    ("Cubic iter's [min/max/total]: ",
     measMinCubicIter,"/",measMaxCubicIter,"/",measTotalCubicIter);
    Output
    ("Cubic failures [min/max/total]: ",
     measMinCubicFails,"/",measMaxCubicFails,"/",measTotalCubicFails);
    Output("");

    // Now compute the singular values and vectors. We recompute the singular
    // values to avoid interfering with the timing experiment above.
    Matrix<Real> U, V;
    timer.Start();
    SecularSVD( d, rho, z, U, s, V, ctrl );
    const double secularSVDTime = timer.Stop();
    Output("Singular SVD: ",secularSVDTime," seconds");
    if( print )
    {
        Print( U, "U" );
        Print( s, "s" );
        Print( V, "V" );
    }

    // Explicitly form the matrix M
    Matrix<Real> M;
    Zeros( M, n, n );
    for( Int j=0; j<n; ++j )
        M(0,j) = z(j)*Sqrt(rho);
    for( Int j=1; j<n; ++j )
        M(j,j) = d(j);
    const Real MFrob = FrobeniusNorm( M );
    Output("|| M ||_F = ",MFrob);
    if( print )
        Print( M, "M" );

    // Test the Singular Value Decomposition of M
    Matrix<Real> UScaled( U );
    DiagonalScale( RIGHT, NORMAL, s, UScaled );
    Matrix<Real> E( M );
    Gemm( NORMAL, ADJOINT, Real(-1), UScaled, V, Real(1), E );
    const Real EFrob = FrobeniusNorm( E );
    Output("|| M - U Sigma V' ||_F = ",EFrob);

    // Test the orthonormality of U and V
    Identity( E, n, n );
    Gemm( NORMAL, ADJOINT, Real(-1), U, U, Real(1), E );
    const Real UOrthError = FrobeniusNorm( E );
    Output("|| I - U U' ||_F = ",UOrthError);
    Identity( E, n, n );
    Gemm( NORMAL, ADJOINT, Real(-1), V, V, Real(1), E );
    const Real VOrthError = FrobeniusNorm( E );
    Output("|| I - V V' ||_F = ",VOrthError);

    if( testFull )
    {
        Matrix<Real> A, w;
        Matrix<Real> dSquared;
        Hadamard( d, d, dSquared );
        Diagonal( A, dSquared );
        Syrk( LOWER, NORMAL, rho, z, Real(1), A );
        timer.Start();
        HermitianEig( LOWER, A, w );
        const Real fullTime = timer.Stop();
        Output("Full Hermitian: ",fullTime," seconds");
        if( print )
            Print( w, "w" );

        auto wDiff( w );
        wDiff -= wSecular;
        const Real diffNorm = FrobeniusNorm( wDiff );
        Output("|| w - wSecular ||_F = ", diffNorm);
        Output("");
    }

    TestBidiagonalCombine<Real>
    ( n, true, maxIter, maxCubicIter, negativeFix, progress, print );
    TestBidiagonalCombine<Real>
    ( n, false, maxIter, maxCubicIter, negativeFix, progress, print );
}

template<typename Real,typename=EnableIf<IsBlasScalar<Real>>>
void TestSecular
( Int n, 
  Int maxIter,
  Int maxCubicIter,
  FlipOrClip negativeFix,
  bool progress,
  bool print,
  bool testFull,
  bool lapack )
{
    Matrix<Real> d, z;
    Real rho;
    GenerateData( n, d, rho, z, print );

    TestSecularHelper<Real>
    ( d, rho, z, maxIter, maxCubicIter, negativeFix, progress, print,
      testFull );
    if( lapack )
    {
        TestLAPACK( d, rho, z );
    }
}

template<typename Real,typename=DisableIf<IsBlasScalar<Real>>,typename=void>
void TestSecular
( Int n, 
  Int maxIter,
  Int maxCubicIter,
  FlipOrClip negativeFix,
  bool progress,
  bool print,
  bool testFull )
{
    Matrix<Real> d, z;
    Real rho;
    GenerateData( n, d, rho, z, print );

    TestSecularHelper<Real>
    ( d, rho, z, maxIter, maxCubicIter, negativeFix, progress, print,
      testFull );
}

int main( int argc, char* argv[] )
{
    Environment env( argc, argv );

    try
    {
        const Int n = Input("--n","matrix size",100);
        const Int maxIter = Input("--maxIter","max iterations",400);
        const Int maxCubicIter = Input("--maxCubicIter","max cubic iter's",40);
        const Int flipOrClipInt = Input("--flipOrClip","0: flip, 1: clip",1);
        const bool progress = Input("--progress","print progress?",false);
        const bool testFull = Input("--testFull","test full eigensolver?",true);
        const bool lapack =
          Input("--lapack","test against LAPACK's secular solver?",true);
        const bool print = Input("--print","print matrices?",false);
#ifdef EL_HAVE_MPC
        const mpfr_prec_t prec = Input("--prec","MPFR precision",256);
#endif
        ProcessInput();

        FlipOrClip negativeFix = static_cast<FlipOrClip>(flipOrClipInt);

        TestSecular<float>
        ( n, maxIter, maxCubicIter, negativeFix, progress, print, testFull,
          lapack );
        TestSecular<double>
        ( n, maxIter, maxCubicIter, negativeFix, progress, print, testFull,
          lapack );

#ifdef EL_HAVE_QD
        TestSecular<DoubleDouble>
        ( n, maxIter, maxCubicIter, negativeFix, progress, print, testFull );
        TestSecular<QuadDouble>
        ( n, maxIter, maxCubicIter, negativeFix, progress, print, testFull );
#endif
#ifdef EL_HAVE_QUAD
        TestSecular<Quad>
        ( n, maxIter, maxCubicIter, negativeFix, progress, print, testFull );
#endif
#ifdef EL_HAVE_MPC
        mpfr::SetPrecision( prec );
        TestSecular<BigFloat>
        ( n, maxIter, maxCubicIter, negativeFix, progress, print, testFull );
#endif

    }
    catch( std::exception& e ) { ReportException(e); }

    return 0;
}
