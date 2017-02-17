/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_LQ_PANELHOUSEHOLDER_HPP
#define EL_LQ_PANELHOUSEHOLDER_HPP

namespace El {
namespace lq {

template<typename Field>
void
PanelHouseholder
( Matrix<Field>& A,
  Matrix<Field>& householderScalars,
  Matrix<Base<Field>>& signature )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    householderScalars.Resize( minDim, 1 );
    signature.Resize( minDim, 1 );

    Matrix<Field> z21;

    for( Int k=0; k<minDim; ++k )
    {
        const Range<Int> ind1( k ), ind2( k+1, END ), indR( k, n );

        Field& alpha11 = A(k,k);
        auto a12 = A( ind1, ind2 );
        auto a1R = A( ind1, indR );
        auto A2R = A( ind2, indR );

        // Find tau and v such that
        //  |alpha11 a12| /I - tau |1  | |1 conj(v)|\ = |beta 0|
        //                \        |v^T|            /
        const Field tau = RightReflector( alpha11, a12 );
        householderScalars(k) = tau;

        // Temporarily set a1R = | 1 v |
        const Field alpha = alpha11;
        alpha11 = 1;

        // A2R := A2R Hous(a1R^T,tau)
        //      = A2R (I - tau a1R^T conj(a1R))
        //      = A2R - tau (A2R a1R^T) conj(a1R)
        Zeros( z21, A2R.Height(), 1 );
        Gemv( NORMAL, Field(1), A2R, a1R, Field(0), z21 );
        Ger( -tau, z21, a1R, A2R );

        // Reset alpha11's value
        alpha11 = alpha;
    }
    // Form d and rescale L
    auto L = A( ALL, IR(0,minDim) );
    GetRealPartOfDiagonal(L,signature);
    auto sgn = []( const Real& delta )
               { return delta >= Real(0) ? Real(1) : Real(-1); };
    EntrywiseMap( signature, MakeFunction(sgn) );
    DiagonalScaleTrapezoid( RIGHT, LOWER, NORMAL, signature, L );
}

template<typename Field>
void
PanelHouseholder
( DistMatrix<Field>& A,
  AbstractDistMatrix<Field>& householderScalars,
  AbstractDistMatrix<Base<Field>>& signature )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(AssertSameGrids( A, householderScalars, signature ))
    typedef Base<Field> Real;
    const Grid& grid = A.Grid();
    DistMatrix<Field,STAR,MR  > a1R_STAR_MR(grid);
    DistMatrix<Field,MC,  STAR> z21_MC_STAR(grid);

    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    householderScalars.Resize( minDim, 1 );

    for( Int k=0; k<minDim; ++k )
    {
        const Range<Int> ind1( k ), ind2( k+1, END ), indR( k, n );

        auto alpha11 = A( ind1, ind1 );
        auto a12     = A( ind1, ind2 );
        auto a1R     = A( ind1, indR );
        auto A2R     = A( ind2, indR );

        // Find tau and v such that
        //  |alpha11 a12| /I - tau |1  | |1 conj(v)|\ = |beta 0|
        //                \        |v^T|            /
        const Field tau = RightReflector( alpha11, a12 );
        householderScalars.Set( k, 0, tau );

        // Temporarily set a1R = | 1 v |
        Field alpha = 0;
        if( alpha11.IsLocal(0,0) )
        {
            alpha = alpha11.GetLocal(0,0);
            alpha11.SetLocal(0,0,Field(1));
        }

        // A2R := A2R Hous(a1R^T,tau)
        //      = A2R (I - tau a1R^T conj(a1R))
        //      = A2R - tau (A2R a1R^T) conj(a1R)
        a1R_STAR_MR.AlignWith( A2R );
        a1R_STAR_MR = a1R;
        z21_MC_STAR.AlignWith( A2R );
        Zeros( z21_MC_STAR, A2R.Height(), 1 );
        LocalGemv( NORMAL, Field(1), A2R, a1R_STAR_MR, Field(0), z21_MC_STAR );
        El::AllReduce( z21_MC_STAR, A2R.RowComm() );
        Ger
        ( -tau, z21_MC_STAR.LockedMatrix(),
                a1R_STAR_MR.LockedMatrix(),
          A2R.Matrix() );

        // Reset alpha11's value
        if( alpha11.IsLocal(0,0) )
            alpha11.SetLocal(0,0,alpha);
    }
    // Form d and rescale L
    auto L = A( ALL, IR(0,minDim) );
    GetRealPartOfDiagonal(L,signature);
    auto sgn = []( const Real& delta )
               { return delta >= Real(0) ? Real(1) : Real(-1); };
    EntrywiseMap( signature, MakeFunction(sgn) );
    DiagonalScaleTrapezoid( RIGHT, LOWER, NORMAL, signature, L );
}

} // namespace lq
} // namespace El

#endif // ifndef EL_LQ_PANELHOUSEHOLDER_HPP
