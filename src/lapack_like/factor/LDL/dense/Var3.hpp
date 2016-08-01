/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_LDL_VAR3_HPP
#define EL_LDL_VAR3_HPP

namespace El {
namespace ldl {

// Unblocked serial LDL _without_ partial pivoting
//
// Since it is significantly more likely that this 
// (generally unstable) routine would fail without encountering
// an exactly zero pivot, it is likely not worth the overhead of
// exception handling to detect zero pivots
template<typename F> 
void Var3Unb( Matrix<F>& A, bool conjugate=false )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( A.Height() != A.Width() )
          LogicError("A must be square");
    )
    typedef Base<F> Real;
    const Int n = A.Height();
    const Int ALDim = A.LDim();

    if( conjugate )
    {
        for( Int j=0; j<n; ++j )
        {
            const Int a21Height = n - (j+1);

            const Real alpha11 = RealPart(A(j,j));

            DEBUG_ONLY(
              if( alpha11 == Real(0) )
                  throw ZeroPivotException();
            )
            const Real alpha11Inv = Real(1)/alpha11;
            F* a21 = A.Buffer(j+1,j  );
            F* A22 = A.Buffer(j+1,j+1);

            blas::Her( 'L', a21Height, -alpha11Inv, a21, 1, A22, ALDim );
            blas::Scal( a21Height, alpha11Inv, a21, 1 );
        }
    }
    else
    {
        for( Int j=0; j<n; ++j )
        {
            const Int a21Height = n - (j+1);

            const F alpha11 = A(j,j);
            DEBUG_ONLY(
              if( alpha11 == F(0) )
                  throw ZeroPivotException();
            )
            const F alpha11Inv = Real(1)/alpha11;
            F* a21 = A.Buffer(j+1,j  );
            F* A22 = A.Buffer(j+1,j+1);

            blas::Syr( 'L', a21Height, -alpha11Inv, a21, 1, A22, ALDim );
            blas::Scal( a21Height, alpha11Inv, a21, 1 );
        }
    }
}

// Blocked serial LDL _without_ partial pivoting
template<typename F>
void Var3( Matrix<F>& A, bool conjugate=false )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( A.Height() != A.Width() )
          LogicError("A must be square");
    )
    const Int n = A.Height();
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );

    Matrix<F> d1, S21;
    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const Range<Int> ind1( k,    k+nb ),
                         ind2( k+nb, n    );

        auto A11 = A( ind1, ind1 );
        auto A21 = A( ind2, ind1 );
        auto A22 = A( ind2, ind2 );

        ldl::Var3Unb( A11, conjugate );
        GetDiagonal( A11, d1 );
        Trsm( RIGHT, LOWER, orientation, UNIT, F(1), A11, A21 );
        S21 = A21;
        DiagonalSolve( RIGHT, NORMAL, d1, A21 );
        Trrk( LOWER, NORMAL, orientation, F(-1), S21, A21, F(1), A22 );
    }
}

template<typename F>
void Var3( ElementalMatrix<F>& APre, bool conjugate=false )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( APre.Height() != APre.Width() )
          LogicError("A must be square");
    )
    const Grid& g = APre.Grid();

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.Get();

    const Int n = A.Height();
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );

    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g), d1_STAR_STAR(g);
    DistMatrix<F,VC,  STAR> A21_VC_STAR(g);
    DistMatrix<F,VR,  STAR> A21_VR_STAR(g);
    DistMatrix<F,STAR,MC  > S21Trans_STAR_MC(g);
    DistMatrix<F,STAR,MR  > A21Trans_STAR_MR(g);

    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);
        
        const Range<Int> ind1( k,    k+nb ),
                         ind2( k+nb, n    );

        auto A11 = A( ind1, ind1 );
        auto A21 = A( ind2, ind1 );
        auto A22 = A( ind2, ind2 );

        A11_STAR_STAR = A11;
        LDL( A11_STAR_STAR, conjugate );
        GetDiagonal( A11_STAR_STAR, d1_STAR_STAR );
        A11 = A11_STAR_STAR;

        A21_VC_STAR.AlignWith( A22 );
        A21_VC_STAR = A21;
        LocalTrsm
        ( RIGHT, LOWER, orientation, UNIT,
          F(1), A11_STAR_STAR, A21_VC_STAR );

        S21Trans_STAR_MC.AlignWith( A22 );
        Transpose( A21_VC_STAR, S21Trans_STAR_MC );
        DiagonalSolve( RIGHT, NORMAL, d1_STAR_STAR, A21_VC_STAR );
        A21_VR_STAR.AlignWith( A22 );
        A21_VR_STAR = A21_VC_STAR;
        A21Trans_STAR_MR.AlignWith( A22 );
        Transpose( A21_VR_STAR, A21Trans_STAR_MR, conjugate );
        LocalTrrk
        ( LOWER, TRANSPOSE,
          F(-1), S21Trans_STAR_MC, A21Trans_STAR_MR, F(1), A22 );

        A21 = A21_VC_STAR;
    }
}

} // namespace ldl
} // namespace El

#endif // ifndef EL_LDL_VAR3_HPP
