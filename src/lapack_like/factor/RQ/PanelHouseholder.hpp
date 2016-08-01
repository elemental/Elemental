/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_RQ_PANEL_HPP
#define EL_RQ_PANEL_HPP

namespace El {
namespace rq {

template<typename F> 
void
PanelHouseholder
( Matrix<F>& A,
  Matrix<F>& phase,
  Matrix<Base<F>>& signature )
{
    DEBUG_CSE
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    const Int iOff = m-minDim;
    const Int jOff = n-minDim;

    phase.Resize( minDim, 1 );
    signature.Resize( minDim, 1 );

    Matrix<F> z01;
    for( Int k=minDim-1; k>=0; --k )
    {
        const Int ki = k + iOff;
        const Int kj = k + jOff;

        const Range<Int> ind0Vert( 0,  ki   ), ind0Horz( 0,  kj   ),
                         ind1Vert( ki, ki+1 ), ind1Horz( kj, kj+1 ),
                                               indL(     0,  kj+1 );

        auto a10     = A( ind1Vert, ind0Horz );
        auto alpha11 = A( ind1Vert, ind1Horz );
        auto A0L     = A( ind0Vert, indL     );
        auto a1L     = A( ind1Vert, indL     );

        // Find tau and v such that
        //  |a10 alpha11| /I - tau |v^T| |conj(v) 1|\ = |0 beta|
        //                \        |1  |            /
        const F tau = RightReflector( alpha11, a10 );
        phase(k) = tau;

        // Temporarily set a1L = | v 1 |
        const F alpha = alpha11(0);
        alpha11(0) = 1;

        // A2R := A2R Hous(a1L^T,tau)
        //      = A2R (I - tau a1L^T conj(a1L))
        //      = A2R - tau (A2R a1L^T) conj(a1L)
        Zeros( z01, ki, 1 );
        Gemv( NORMAL, F(1), A0L, a1L, F(0), z01 );
        Ger( -tau, z01, a1L, A0L );

        // Reset alpha11's value
        alpha11(0) = alpha;
    }
    // Form d and rescale R
    auto R = A( IR(0,m), IR(jOff,n) );
    GetRealPartOfDiagonal(R,signature);
    auto sgn = []( const Real& delta )
               { return delta >= Real(0) ? Real(1) : Real(-1); };
    EntrywiseMap( signature, function<Real(Real)>(sgn) );
    DiagonalScaleTrapezoid( RIGHT, UPPER, NORMAL, signature, R, -iOff );
}

template<typename F> 
void
PanelHouseholder
( DistMatrix<F>& A,
  ElementalMatrix<F>& phase,
  ElementalMatrix<Base<F>>& signature )
{
    DEBUG_CSE
    DEBUG_ONLY(AssertSameGrids( A, phase, signature ))
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    DistMatrix<F,STAR,MR  > a1L_STAR_MR(g);
    DistMatrix<F,MC,  STAR> z01_MC_STAR(g);

    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    const Int iOff = m-minDim;
    const Int jOff = n-minDim;
    phase.Resize( minDim, 1 );

    for( Int k=minDim-1; k>=0; --k )
    {
        const Int ki = k + iOff;
        const Int kj = k + jOff;

        const Range<Int> ind0Vert( 0,  ki   ), ind0Horz( 0,  kj   ),
                         ind1Vert( ki, ki+1 ), ind1Horz( kj, kj+1 ),
                                               indL(     0,  kj+1 );

        auto a10     = A( ind1Vert, ind0Horz );
        auto alpha11 = A( ind1Vert, ind1Horz );
        auto A0L     = A( ind0Vert, indL     );
        auto a1L     = A( ind1Vert, indL     );

        // Find tau and v such that
        //  |a10 alpha11| /I - tau |v^T| |conj(v) 1|\ = |0 beta|
        //                \        |1  |            /
        const F tau = RightReflector( alpha11, a10 );
        phase.Set( k, 0, tau );

        // Temporarily set a1L = | v 1 |
        F alpha = 0;
        if( alpha11.IsLocal(0,0) )
        {
            alpha = alpha11.GetLocal(0,0);
            alpha11.SetLocal(0,0,1);
        }

        // A2R := A2R Hous(a1L^T,tau)
        //      = A2R (I - tau a1L^T conj(a1L))
        //      = A2R - tau (A2R a1L^T) conj(a1L)
        a1L_STAR_MR.AlignWith( A0L );
        a1L_STAR_MR = a1L;
        z01_MC_STAR.AlignWith( A0L );
        Zeros( z01_MC_STAR, ki, 1 );
        LocalGemv( NORMAL, F(1), A0L, a1L_STAR_MR, F(0), z01_MC_STAR );
        El::AllReduce( z01_MC_STAR, A0L.RowComm() );
        Ger
        ( -tau, z01_MC_STAR.LockedMatrix(), a1L_STAR_MR.LockedMatrix(),
          A0L.Matrix() ); 

        // Reset alpha11's value
        if( alpha11.IsLocal(0,0) )
            alpha11.SetLocal(0,0,alpha);
    }
    // Form d and rescale R
    auto R = A( IR(0,m), IR(jOff,n) );
    GetRealPartOfDiagonal(R,signature);
    auto sgn = []( const Real& delta )
               { return delta >= Real(0) ? Real(1) : Real(-1); };
    EntrywiseMap( signature, function<Real(Real)>(sgn) );
    DiagonalScaleTrapezoid( RIGHT, UPPER, NORMAL, signature, R, -iOff );
}

} // namespace rq
} // namespace El

#endif // ifndef EL_RQ_PANEL_HPP
