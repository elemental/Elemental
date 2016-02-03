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

template<typename F>
inline void
DiagonalBlockSolve
(       Matrix<F>& T,
  const Matrix<F>& shifts,
        Matrix<F>& X ) 
{
    DEBUG_ONLY(
      CSE cse("triang_eig::DiagonalBlockSolve");
      if( shifts.Height() != X.Width() )
          LogicError("Incompatible number of shifts");
    )
    const char uploChar = 'U';
    const char orientChar = 'N';
    auto diag = GetDiagonal(T);
    const Int n = T.Height();
    const Int ldim = T.LDim();
    const Int numShifts = shifts.Height();

    const Base<F> eps = lapack::MachineEpsilon<Base<F>>(); 
    
    for( Int j=1; j<numShifts; ++j )
    {
        ShiftDiagonal( T, -shifts.Get(j,0) );

        // TODO: Handle small diagonal entries in the usual manner
	auto minDiag = 2*eps*Abs(shifts.Get(j,0));
	for( Int i=0; i<Min(n,j); ++i)
	{
	    if( Abs(T.Get(i,i))<minDiag )
		T.Set(i,i,F(minDiag));
	}
	
        blas::Trsv
        ( uploChar, orientChar, 'N', Min(n,j), 
          T.LockedBuffer(), ldim, X.Buffer(0,j), 1 );
        SetDiagonal( T, diag );
    }
}

} // namespace triang_eig

#if 0
  // OLD VERSION
  // TOOD: remove this
template<typename F>
inline void
TriangEig( Matrix<F>& U, Matrix<F>& X ) 
{
    DEBUG_ONLY(CSE cse("TriangEig"))
    const Int m = U.Height();
    const Int bsize = Blocksize();
    const Int kLast = LastOffset( m, bsize );

    Matrix<F> shifts;
    GetDiagonal( U, shifts );

    // TODO: Handle near and exact singularity

    // Make X the negative of the strictly upper triangle of  U
    X = U;
    MakeTrapezoidal( UPPER, X, 1 );
    Scale( F(-1), X );

    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,m-k);

        const Range<Int> ind0( 0, k    ),
                         ind1( k, k+nb );

        auto U01 = U( ind0, ind1 );
        auto U11 = U( ind1, ind1 );

        auto X0 = X( ind0, IR(k,END) );
        auto X1 = X( ind1, IR(k,END) );

        triang_eig::DiagonalBlockSolve( U11, shifts(IR(k,END),ALL), X1 );
        Gemm( NORMAL, NORMAL, F(-1), U01, X1, F(1), X0 );
    }
    FillDiagonal( X, F(1) ); 

    // Normalize eigenvectors
    for( Int k=1; k<m; ++k )
    {
        auto Xcol = View( X, IR(0,k+1), IR(k,k+1) );
	Scale( F(1)/Nrm2(Xcol), Xcol );
    }

}

#else
  // NEW VERSION
  // TOOD: make this work
template<typename F>
inline void
TriangEig( Matrix<F>& U, Matrix<F>& X ) 
{
    DEBUG_ONLY(CSE cse("TriangEig"))
    const Int n = U.Height();
    const Int bsize = Blocksize();
    const Int kLast = LastOffset( n, bsize );

    Matrix<F> shifts;
    GetDiagonal( U, shifts );

    Matrix<F> scalings;
    Ones( scalings, n, 1 );
    
    // TODO: Handle near and exact singularity

    // Make X the negative of the strictly upper triangle of  U
    X = U;
    MakeTrapezoidal( UPPER, X, 1 );
    Scale( F(-1), X );

    // TODO: blocked algorithm to get BLAS 3
    for( Int k=n-1; k>=0; --k )
    {
        const Range<Int> ind0( 0, k   ),
	                 ind1( k, k+1 );

	auto U01 = U( ind0, ind1 );
	auto U11 = U( ind1, ind1 );
	
    }
    
    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const Range<Int> ind0( 0, k    ),
                         ind1( k, k+nb );

        auto U01 = U( ind0, ind1 );
        auto U11 = U( ind1, ind1 );

        auto X0 = X( ind0, IR(k,END) );
        auto X1 = X( ind1, IR(k,END) );

        triang_eig::DiagonalBlockSolve( U11, shifts(IR(k,END),ALL), X1 );
        Gemm( NORMAL, NORMAL, F(-1), U01, X1, F(1), X0 );
    }
    FillDiagonal( X, F(1) ); 

    // Normalize eigenvectors
    for( Int k=1; k<n; ++k )
    {
        auto Xcol = View( X, IR(0,k+1), IR(k,k+1) );
	Scale( F(1)/Nrm2(Xcol), Xcol );
    }

}
#endif
  
template<typename F>
inline void
TriangEig
( const ElementalMatrix<F>& UPre, 
        ElementalMatrix<F>& XPre ) 
{

    DEBUG_ONLY(CSE cse("TriangEig"))

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
