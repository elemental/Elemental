/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

namespace triang_eig {

/* Determine machine dependent parameters to control overflow
 *   Note: LAPACK uses more complicated parameters to handle 
 *   issues that can happen on Cray machines.
 */
template<typename Real>
inline void
OverflowParameters( Real& smlnum, Real& bignum )
{
    const Real unfl = lapack::MachineSafeMin<Real>();
    const Real ovfl = lapack::MachineOverflowThreshold<Real>();
    const Real ulp  = lapack::MachinePrecision<Real>();
    smlnum = Max( unfl/ulp, 1/(ovfl*ulp) );
    bignum = 1/smlnum;
}

template<typename F>
void
DiagonalBlockSolve
(       Matrix<F>& U,
  const Matrix<F>& shifts,
        Matrix<F>& X,
	Matrix<F>& scales )
{

    typedef Base<F> Real;
  
    DEBUG_ONLY(
      CSE cse("triang_eig::DiagonalBlockSolve");
      if( shifts.Height() != X.Width() )
          LogicError("Incompatible number of shifts");
    )
    auto diag = GetDiagonal(U);
    const Int n = U.Height();
    const Int ldim = U.LDim();
    const Int numShifts = shifts.Height();
    Real smlnum, bignum;
    OverflowParameters<Real>( smlnum, bignum );
    
    // Default scale is 1
    Ones( scales, numShifts, 1 );

    // Compute infinity norms of columns of U (excluding diagonal)
    // TODO: scale cnorm if an entry is bigger than bignum
    Matrix<Real> cnorm( n, 1 );
    cnorm.Set( 0, 0, Real(0) );
    for( Int j=1; j<n; ++j )
    {
        cnorm.Set( j, 0, MaxNorm( U(IR(0,j),IR(j)) ) );
    }

    // Iterate through RHS's
    for( Int j=1; j<numShifts; ++j )
    {

        // Initialize triangular system
        ShiftDiagonal( U, -shifts.Get(j,0) );
	auto xj = X( IR(0,Min(n,j)), IR(j) );

	// Determine largest entry of RHS
	Real xjMax = MaxNorm( xj );
	if( xjMax >= bignum )
	{
	    const Real s = Real(0.5)*bignum/xjMax;
	    xj *= s;
	    xjMax *= s;
	    scales.Set( j, 0, s*scales.Get(j,0) );
	}
	xjMax = Max( xjMax, 2*smlnum );

	// Estimate growth of entries in triangular solve
	//   Note: See "Robust Triangular Solves for Use in Condition
	//   Estimation" by Edward Anderson for explanation of bounds.
	Real invGi = 1/xjMax;
	Real invMi = invGi;
	for( Int i=Min(n,j)-1; i>=0; --i )
	{
	    const Real absUii = SafeAbs( U.Get(i,i) );
	    if( invGi<=smlnum || invMi<=smlnum || absUii<=smlnum )
	    {
	        invGi = 0;
	        break;
	    }
	    invMi = Min( invMi, absUii*invGi );
	    if( i > 0 )
	    {
	        invGi *= absUii/(absUii+cnorm.Get(i,0));
	    }
	}
	invGi = Min( invGi, invMi );

	// Call TRSV if estimated growth is not too large
	if( invGi > smlnum )
	{
	    blas::Trsv
	    ( 'U', 'N', 'N', Min(n,j),
	      U.LockedBuffer(), ldim, X.Buffer(0,j), 1 );
	}

	// Perform backward substitution if estimated growth is large
	else
	{
	    for( Int i=Min(n,j)-1; i>=0; --i )
	    {

		// Perform division and check for overflow
		const Real absUii = SafeAbs( U.Get(i,i) );
		Real absXij = SafeAbs( xj.Get(i,0) );
		if( absUii > smlnum )
		{
		    if( absUii<=1 && absXij>=absUii*bignum )
		    {
			// Set overflowing entry to 0.5/U[i,i]
		        const Real s = Real(0.5)/absXij;
			scales.Set( j, 0, s*scales.Get(j,0) );
			xj *= s;
			xjMax *= s;
		    }
		    xj.Set( i, 0, xj.Get(i,0)/U.Get(i,i) );
		}
		else if( absUii > 0 )
		{
		    if( absXij >= absUii*bignum )
		    {
			// Set overflowing entry to bignum/2
		        const Real s = Real(0.5)*absUii*bignum/absXij;
			scales.Set( j, 0, s*scales.Get(j,0) );
			xj *= s;
			xjMax *= s;
		    }
		    xj.Set( i, 0, xj.Get(i,0)/U.Get(i,i) );
		}
		else
		{
		    // TODO: maybe this tolerance should be loosened to
		    //   | Xij | >= || A || * eps
		    if( absXij >= smlnum )
		    {
		        scales.Set( j, 0, F(0) );
			Zero( xj );
			xjMax = 0;
		    }
		    xj.Set( i, 0, F(1) );
		}

		if( i > 0 )
	        {

		    // Check for possible overflows in AXPY
		    // Note: G(i+1) <= G(i) + | Xij | * cnorm(i)
		    absXij = SafeAbs( xj.Get(i,0) );
		    if( absXij >= 1 &&
			cnorm.Get(i,0) >= (bignum-xjMax)/absXij )
		    {
		        const Real s = Real(0.25)/absXij;
			scales.Set( j, 0, s*scales.Get(j,0) );
			xj *= s;
			xjMax *= s;
			absXij *= s;
		    }
		    else if( absXij < 1 &&
			     absXij*cnorm.Get(i,0) >= bignum-xjMax )
		    {
		        const Real s = Real(0.25);
		        scales.Set( j, 0, s*scales.Get(j,0) );
			xj *= s;
			xjMax *= s;
			absXij *= s;
		    }
		    xjMax += absXij * cnorm.Get(i,0);

		    // AXPY
		    auto U01 = U( IR(0,i), IR(i) );
		    auto X1  = X( IR(0,i), IR(j) );
		    Axpy( -xj.Get(i,0), U01, X1 );

		}

	    }

	}

	// Reset matrix diagonal
        SetDiagonal( U, diag );
    }
}
  
} // namespace triang_eig

template<typename F>
inline void
TriangEig( Matrix<F>& U, Matrix<F>& X ) 
{
    typedef Base<F> Real;
  
    DEBUG_ONLY(CSE cse("TriangEig"))
    const Int m = U.Height();
    const Int bsize = Blocksize();
    const Int kLast = LastOffset( m, bsize );
    Real smlnum, bignum;
    triang_eig::OverflowParameters<Real>( smlnum, bignum );

    Matrix<F> shifts( m, 1 ),
              scales( m, 1 ),
              scalesUpdate( m, 1 );
    GetDiagonal( U, shifts );
    Ones( scales, m, 1 );
    
    // Make X the negative of the strictly upper triangle of  U
    X = U;
    MakeTrapezoidal( UPPER, X, 1 );
    Scale( F(-1), X );

    // Determine largest entry of each RHS
    Matrix<Real> XMax( m, 1 );
    XMax.Set( 0, 0, Real(0) );
    for( Int j=1; j<m; ++j )
    {
        auto xj = X( IR(0,j), IR(j) );
	Real xjMax = MaxNorm( xj );
	xjMax = Max( xjMax, 2*smlnum );
        XMax.Set( j, 0, xjMax );
    }

    // Block triangular solve
    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,m-k);

        const Range<Int> ind0( 0, k      ),
	                 ind1( k, k+nb   ),
	                 ind2( k+nb, END ),
	                 cols( k, END    );
	
        auto U01 = U( ind0, ind1 );
        auto U11 = U( ind1, ind1 );

        auto X0 = X( ind0, cols );
        auto X1 = X( ind1, cols );
	auto X2 = X( ind2, cols );

	auto scalesUpdateCurrent = scalesUpdate( cols, ALL );
	
	// Perform triangular solve on diagonal block
        triang_eig::DiagonalBlockSolve( U11, shifts(cols,ALL),
					X1, scalesUpdateCurrent );

	// Apply scalings on RHS
	for( Int j=k; j<m; ++j )
	{
	    const Real sj = scalesUpdate.GetRealPart(j,0);
	    if( sj < 1 )
	    {
	        scales.Set( j, 0, sj*scales.Get(j,0) );
	        auto X0j = X0( ALL, IR(j-k) );
	        auto X2j = X2( ALL, IR(j-k) );
		X0j *= sj;
		X2j *= sj;
		XMax.Set( j, 0, sj*XMax.Get(j,0) );
	    }
	}

	if( k > 0 )
	{
	    // Compute infinity norms of columns in U01
	    // Note: nb*cnorm is the sum of infinity norms
	    // TODO: scale cnorm if an entry is bigger than bignum
	    Real cnorm = 0;
	    for( Int j=0; j<nb; ++j )
	    {
	        cnorm += MaxNorm( U01(ALL,IR(j)) ) / nb;
	    }

	    // Check for possible overflows in GEMM
	    // Note: G(i+1) <= G(i) + nb*cnorm*|| X1[:,j] ||_infty
	    for( Int j=k; j<m; ++j )
	    {
	        auto xj = X( IR(0,j), IR(j) );
	        Real xjMax = XMax.Get(j,0);
		Real X1Max = MaxNorm( X1(ALL,IR(j-k)) );
		if( X1Max >= 1 &&
		    cnorm >= (bignum-xjMax)/X1Max/nb )
		{
		    const Real s = Real(0.5)/X1Max/nb;
		    scales.Set( j, 0, s*scales.Get(j,0) );
		    xj *= s;
		    xjMax *= s;
		    X1Max *= s;
		}
		else if( X1Max < 1 &&
			 cnorm*X1Max >= (bignum-xjMax)/nb )
		{
		    const Real s = Real(0.5)/nb;
		    scales.Set( j, 0, s*scales.Get(j,0) );
		    xj *= s;
		    xjMax *= s;
		    X1Max *= s;
		}
		xjMax += nb*cnorm*X1Max;
		XMax.Set( j, 0, xjMax );
	    }

	    // Update RHS with GEMM
	    Gemm( NORMAL, NORMAL, F(-1), U01, X1, F(1), X0 );
	}
    }
    SetDiagonal( X, scales );

    // Normalize eigenvectors
    for( Int j=0; j<m; ++j )
    {
        auto Xj = X( IR(0,j+1), IR(j) );
	Scale( 1/Nrm2(Xj), Xj );
    }

}
  
template<typename F>
inline void
TriangEig
( const ElementalMatrix<F>& UPre, 
        ElementalMatrix<F>& XPre ) 
{

    DEBUG_ONLY(CSE cse("TriangEig"))

// TODO: implement distributed version of TriangEig
#if 1
    LogicError("TriangEig is not yet implemented for distributed matrices");
#else
    DistMatrixReadProxy<F,F,MC,MR> UProx( UPre );
    DistMatrixWriteProxy<F,F,MC,MR> XProx( XPre );
    auto& U = UProx.GetLocked();
    auto& X = XProx.Get();

    const Grid& g = U.Grid();
    DistMatrix<F,VR,STAR> shifts(g);
    GetDiagonal( U, shifts );

    // TODO: Handle near and exact singularity

    // Make X the negative of the strictly upper triangle of  U
    X = U;
    MakeTrapezoidal( UPPER, X, 1 );
    Scale( F(-1), X );

    DistMatrix<F,MC,  STAR> U01_MC_STAR(g);
    DistMatrix<F,STAR,STAR> U11_STAR_STAR(g);
    DistMatrix<F,STAR,MR  > X1_STAR_MR(g);
    DistMatrix<F,STAR,VR  > X1_STAR_VR(g);

    const Int m = X.Height();
    const Int bsize = Blocksize();
    const Int kLast = LastOffset( m, bsize );

    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,m-k);

        const Range<Int> ind0( 0, k    ),
                         ind1( k, k+nb );

        auto U01 = U( ind0, ind1 );
        auto U11 = U( ind1, ind1 );

        auto X0 = X( ind0, IR(k,END) );
        auto X1 = X( ind1, IR(k,END) );

        // X1[* ,VR] := U11^-1[* ,* ] X1[* ,VR]
        U11_STAR_STAR = U11; // U11[* ,* ] <- U11[MC,MR]
        X1_STAR_VR.AlignWith( shifts );
        X1_STAR_VR = X1; // X1[* ,VR] <- X1[MC,MR]
        triang_eig::DiagonalBlockSolve
        ( U11_STAR_STAR.Matrix(), shifts(IR(k,END),ALL).LockedMatrix(), 
          X1_STAR_VR.Matrix() );

        X1_STAR_MR.AlignWith( X0 );
        X1_STAR_MR = X1_STAR_VR; // X1[* ,MR]  <- X1[* ,VR]
        X1 = X1_STAR_MR; // X1[MC,MR] <- X1[* ,MR]

        // X0[MC,MR] -= U01[MC,* ] X1[* ,MR]
        U01_MC_STAR.AlignWith( X0 );
        U01_MC_STAR = U01; // U01[MC,* ] <- U01[MC,MR]
        LocalGemm( NORMAL, NORMAL, F(-1), U01_MC_STAR, X1_STAR_MR, F(1), X0 );
    }
    FillDiagonal( X, F(1) );
    
    // Normalize eigenvectors
    for( Int k=1; k<m; ++k )
    {
        auto Xcol = View( X, IR(0,k+1), IR(k,k+1) );
	Scale( F(1)/Nrm2(Xcol), Xcol );
    }

#endif
    
}

#define PROTO(F) \
  template void TriangEig \
  (       Matrix<F>& T, \
          Matrix<F>& X ); \
  template void TriangEig \
  ( const ElementalMatrix<F>& T, \
          ElementalMatrix<F>& X );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
