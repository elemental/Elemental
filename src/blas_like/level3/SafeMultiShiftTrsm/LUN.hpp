/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace El {
namespace safemstrsm {

template<typename F>
inline void
LUNBlock
(       Matrix<F>& U,
  const Matrix<F>& shifts,
        Matrix<F>& X,
	Matrix<F>& scales )
{

    typedef Base<F> Real;
  
    DEBUG_ONLY(
      CSE cse("safemstrsm::LUNBlock");
      if( shifts.Height() != X.Width() )
          LogicError("Incompatible number of shifts");
    )
    const char uploChar = 'U';
    const char orientChar = 'N';
    auto diag = GetDiagonal(U);
    const Int n = U.Height();
    const Int ldim = U.LDim();
    const Int numShifts = shifts.Height();

    // Determine machine dependent parameters to control overflow
    //   Note: LAPACK uses more complicated parameters to handle 
    //   issues that can happen on Cray machines.
    const Real unfl = lapack::MachineSafeMin<Real>();
    const Real ovfl = lapack::MachineOverflowThreshold<Real>();
    const Real ulp  = lapack::MachinePrecision<Real>();
    const Real smlnum = std::max( unfl/ulp, 1/(ovfl*ulp) );
    const Real bignum = (1-ulp)/smlnum;
    
    // Default scale is 1
    Ones( scales, numShifts, 1 );

    // Iterate through RHS's
    for( Int j=0; j<numShifts; ++j )
    {
        ShiftDiagonal( U, -shifts.Get(j,0) );
	auto xj = X( ALL, IR(j) );

	// Backward substitution
	for( Int i=n-1; i>=0; --i)
	{

	    // Check if division causes overflow
	    const Real AbsUii = SafeAbs( U.Get(i,i) );
	    const Real AbsXij = SafeAbs( xj.Get(i,0) );
	    if( AbsUii > smlnum )
	    {
	        if( AbsUii<=1 && AbsXij>=AbsUii*bignum )
		{
		    // Set overflowing entry to 1/U[i,i]
		    Real s = 1/AbsXij;
		    scales.Set( j, 0, s*scales.Get(j,0) );
		    xj *= s;
		}
		xj.Set( i, 0, xj.Get(i,0)/U.Get(i,i) );
	    }
	    else if( AbsUii > 0 )
	    {
	        if( AbsXij >= AbsUii*bignum )
		{
		    // Set overflowing entry to bignum
		    Real s = AbsUii*bignum/AbsXij;
		    scales.Set( j, 0, s*scales.Get(j,0) );
		    xj *= s;
		}
		xj.Set( i, 0, xj.Get(i,0)/U.Get(i,i) );
	    }
	    else
	    {
	        // Set overflowing entry to 1
		scales.Set( j, 0, F(0) );
		Zero( xj );
		xj.Set( i, 0, F(1) );
	    }

	    // TODO: handle case to prevent overflow in axpy

	    if( i > 0 )
	    {
	        auto U01 = U( IR(0,i), IR(i) );
		auto X1  = X( IR(0,i), IR(j) );
		Axpy( -xj.Get(i,0), U01, X1 );
	    }

	}

        SetDiagonal( U, diag );
    }
}

template<typename F>
inline void
LUN( Matrix<F>& U, const Matrix<F>& shifts,
     Matrix<F>& X, Matrix<F>& scales ) 
{
    DEBUG_ONLY(CSE cse("safemstrsm::LUN"))
    const Int m = X.Height();
    const Int n = X.Width();
    const Int bsize = Blocksize();
    const Int kLast = LastOffset( m, bsize );

    Ones( scales, n, 1 );
    Matrix<F> scalesUpdate;
    
    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,m-k);

        const Range<Int> ind0( 0,    k    ),
	                 ind1( k,    k+nb ),
	                 ind2( k+nb, END  );

        auto U01 = U( ind0, ind1 );
        auto U11 = U( ind1, ind1 );

        auto X0 = X( ind0, ALL );
        auto X1 = X( ind1, ALL );
	auto X2 = X( ind2, ALL );

	LUNBlock( U11, shifts, X1, scalesUpdate );
	DiagonalScale( LEFT, NORMAL, scalesUpdate, scales );
	DiagonalScale( RIGHT, NORMAL, scalesUpdate, X0 );
	DiagonalScale( RIGHT, NORMAL, scalesUpdate, X2 );
        Gemm( NORMAL, NORMAL, F(-1), U01, X1, F(1), X0 );
    }
}

#if 0
  // TODO: template specialization for elemental matrices
template<typename F>
inline void
LUN
( const ElementalMatrix<F>& UPre, 
  const ElementalMatrix<F>& shiftsPre,
        ElementalMatrix<F>& XPre ) 
{
    DEBUG_ONLY(CSE cse("mstrsm::LUN"))

    DistMatrixReadProxy<F,F,MC,MR> UProx( UPre );
    DistMatrixReadProxy<F,F,VR,STAR> shiftsProx( shiftsPre );
    DistMatrixReadWriteProxy<F,F,MC,MR> XProx( XPre );
    auto& U = UProx.GetLocked();
    auto& shifts = shiftsProx.GetLocked();
    auto& X = XProx.Get();

    const Grid& g = U.Grid();
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

        auto X0 = X( ind0, ALL );
        auto X1 = X( ind1, ALL );

        // X1[* ,VR] := U11^-1[* ,* ] X1[* ,VR]
        U11_STAR_STAR = U11; // U11[* ,* ] <- U11[MC,MR]
        X1_STAR_VR.AlignWith( shifts );
        X1_STAR_VR = X1; // X1[* ,VR] <- X1[MC,MR]
        LUN
        ( U11_STAR_STAR.Matrix(), shifts.LockedMatrix(), 
          X1_STAR_VR.Matrix() );

        X1_STAR_MR.AlignWith( X0 );
        X1_STAR_MR = X1_STAR_VR; // X1[* ,MR]  <- X1[* ,VR]
        X1 = X1_STAR_MR; // X1[MC,MR] <- X1[* ,MR]

        // X0[MC,MR] -= U01[MC,* ] X1[* ,MR]
        U01_MC_STAR.AlignWith( X0 );
        U01_MC_STAR = U01; // U01[MC,* ] <- U01[MC,MR]
        LocalGemm( NORMAL, NORMAL, F(-1), U01_MC_STAR, X1_STAR_MR, F(1), X0 );
    }
}
#endif
  
} // namespace safemstrsm
} // namespace El
