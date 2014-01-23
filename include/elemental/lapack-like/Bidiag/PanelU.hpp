/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_BIDIAG_PANELU_HPP
#define ELEM_LAPACK_BIDIAG_PANELU_HPP

#include "elemental/blas-like/level1/Axpy.hpp"
#include "elemental/blas-like/level1/Conjugate.hpp"
#include "elemental/blas-like/level1/Scale.hpp"
#include "elemental/blas-like/level2/Gemv.hpp"
#include "elemental/lapack-like/Reflector.hpp"

namespace elem {
namespace bidiag {

template<typename F> 
inline void
PanelU
( DistMatrix<F>& A, 
  DistMatrix<F,MD,STAR>& tP,
  DistMatrix<F,MD,STAR>& tQ,
  DistMatrix<F>& X, 
  DistMatrix<F>& Y,
  DistMatrix<F,MC,STAR>& AColPan_MC_STAR,
  DistMatrix<F,STAR,MR>& ARowPan_STAR_MR )
{
    const Int mA = A.Height();
    const Int nA = A.Width();
    const Int nX = X.Width();
    DEBUG_ONLY(
        CallStackEntry cse("bidiag::PanelU");
        if( A.Grid() != tP.Grid() || tP.Grid() != tQ.Grid() || 
            tQ.Grid() != X.Grid() || X.Grid() != Y.Grid() ||
            Y.Grid() != AColPan_MC_STAR.Grid() || 
            Y.Grid() != ARowPan_STAR_MR.Grid() )
            LogicError("Grids must match");
        if( tP.Height() != nX || tP.Width() != 1 )
            LogicError("tP was not the right size");
        if( tQ.Height() != nX || tQ.Width() != 1 )
            LogicError("tQ was not the right size");
        if( mA < nA )
            LogicError("A must be at least as tall as it is wide");
        if( mA != X.Height() )
            LogicError("A and X must be the same height");
        if( nA != Y.Height() )
            LogicError("Y must be the same height as A's width");
        if( X.Height() < nX )
            LogicError("X must be a column panel");
        if( Y.Width() != nX )
            LogicError("Y is the wrong width");
        if( A.ColAlign() != X.ColAlign() || 
            A.RowAlign() != X.RowAlign() )
            LogicError("A and X must be aligned");
        if( A.ColAlign() != Y.ColAlign() ||
            A.RowAlign() != Y.RowAlign() )
            LogicError("A and Y must be aligned");
    )
    typedef Base<F> Real;
    const Grid& g = A.Grid();

    DistMatrix<F> q21(g), s01(g);
    DistMatrix<F,MC,  STAR> s01_MC_STAR(g), z01_MC_STAR(g),
                            q21_MC_STAR(g), s21_MC_STAR(g), z21_MC_STAR(g),
                            uB1_MC_STAR(g);
    DistMatrix<F,MR,  STAR> a01_MR_STAR(g), s01_MR_STAR(g), z01_MR_STAR(g),
                            q21_MR_STAR(g), z21_MR_STAR(g),
                            sB1_MR_STAR(g);
    DistMatrix<F,STAR,MR  > a10_STAR_MR(g), y10_STAR_MR(g);
    DistMatrix<F,STAR,MC  > x10_STAR_MC(g), a12_STAR_MC(g);
    DistMatrix<F,MR,  MC  > z01_MR_MC(g), q21_MR_MC(g), z21_MR_MC(g);

    DistMatrix<Real,MD,STAR> d(g), e(g);
    d.AlignWithDiagonal( A.DistData(), 0 );
    e.AlignWithDiagonal( A.DistData(), 1 );
    d.Resize( nX, 1 );
    e.Resize( nX, 1 );

    for( Int k=0; k<nX; ++k )
    {
        auto a01      = ViewRange( A, 0,   k,   k,   k+1 );
        auto A02      = ViewRange( A, 0,   k+1, k,   nA  );
        auto a10      = ViewRange( A, k,   0,   k+1, k   );
        auto alpha11  = ViewRange( A, k,   k,   k+1, k+1 );
        auto a12      = ViewRange( A, k,   k+1, k+1, nA  );
        auto alpha12L = ViewRange( A, k,   k+1, k+1, k+2 );
        auto a12R     = ViewRange( A, k,   k+2, k+1, nA  );
        auto a21      = ViewRange( A, k+1, k,   mA,  k+1 );
        auto A22      = ViewRange( A, k+1, k+1, mA,  nA  );
        auto ABL      = ViewRange( A, k,   0,   mA,  k   );
        auto ABR      = ViewRange( A, k,   k,   mA,  nA  );
        auto aB1      = ViewRange( A, k,   k,   mA,  k+1 );
        auto AB2      = ViewRange( A, k,   k+1, mA,  nA  );
        auto A2L      = ViewRange( A, k+1, 0,   mA,  k+1 );

        auto a12_STAR_MR = ViewRange( ARowPan_STAR_MR, k, k+1, k+1, nA  );
        auto aB1_MC_STAR = ViewRange( AColPan_MC_STAR, k, k,   mA,  k+1 );

        auto x10 = ViewRange( X, k,   0,   k+1, k   );
        auto X20 = ViewRange( X, k+1, 0,   mA,  k   );
        auto x21 = ViewRange( X, k+1, k,   mA,  k+1 );
        auto XBL = ViewRange( X, k,   0,   mA,  k   );

        auto y10 = ViewRange( Y, k,   0,   k+1, k   );
        auto Y20 = ViewRange( Y, k+1, 0,   nA,  k   );
        auto y21 = ViewRange( Y, k+1, k,   nA,  k+1 );
        auto Y2L = ViewRange( Y, k+1, 0,   nA,  k+1 );

        auto delta1   = View( d,  k, 0, 1, 1 );
        auto epsilon1 = View( e,  k, 0, 1, 1 );
        auto tauQ1    = View( tQ, k, 0, 1, 1 );
        auto tauP1    = View( tP, k, 0, 1, 1 );

        const bool thisIsMyRow = ( g.Row() == alpha11.ColAlign() );
        const bool thisIsMyCol = ( g.Col() == alpha11.RowAlign() );
        const bool nextIsMyCol = ( g.Col() == a12.RowAlign() ) ;

        // Update the current column of A:
        //   aB1 := aB1 - ABL y10^H - XBL a01
        if( k > 0 )
        {
            Conjugate( y10 );
            y10_STAR_MR.AlignWith( ABL );
            y10_STAR_MR = y10;
            // uB1[MC,* ] := ABL[MC,MR] y10^H[MR,* ]
            a01_MR_STAR.AlignWith( ABL );
            a01_MR_STAR = a01;
            uB1_MC_STAR.AlignWith( ABL );
            Zeros( uB1_MC_STAR, ABL.Height(), 1 );
            LocalGemv( NORMAL, F(1), ABL, y10_STAR_MR, F(0), uB1_MC_STAR );
            // uB1[MC,* ] := uB1[MC,* ] + XBL[MC,MR] a01[MR,* ]
            //             = ABL[MC,MR] y10^H[MR,* ] + XBL[MC,MR] a01[MR,* ]
            LocalGemv( NORMAL, F(1), XBL, a01_MR_STAR, F(1), uB1_MC_STAR );
            // Sum the partial contributions and subtract from aB1
            aB1.SumScatterUpdate( F(-1), uB1_MC_STAR );
        }

        // Find tauQ, u, and delta such that
        //     I - conj(tauQ) | 1 | | 1, u^H | | alpha11 | = | delta |
        //                    | u |            |   a21   | = |   0   |
        F tauQ = 0;
        if( thisIsMyCol )
        {
            tauQ = reflector::Col( alpha11, a21 );
            if( thisIsMyRow )
            {
                tauQ1.SetLocal(0,0,tauQ);
                // Store delta and force | alpha11 | = | 1 |
                //                       |   a21   |   | u |
                delta1.SetLocal(0,0,alpha11.GetLocalRealPart(0,0));
                alpha11.SetLocal(0,0,F(1));
            }
        }

        //
        // y21 := tauQ ( AB2^H aB1 - A02^H XBL^H aB1 - Y20 ABL^H aB1 )
        //
        aB1_MC_STAR = aB1;
        // z01[MR,* ] := ABL^H[MR,MC] aB1[MC,* ]
        z01_MR_STAR.AlignWith( ABL );
        Zeros( z01_MR_STAR, k, 1 );
        LocalGemv( ADJOINT, F(1), ABL, aB1_MC_STAR, F(0), z01_MR_STAR );
        // z21[MR,* ] := AB2^H[MR,MC] aB1[MC,* ]
        z21_MR_STAR.AlignWith( AB2 );
        Zeros( z21_MR_STAR, A22.Width(), 1 );
        LocalGemv( ADJOINT, F(1), AB2, aB1_MC_STAR, F(0), z21_MR_STAR );
        // Sum the partial contributions
        z01_MR_STAR.SumOverCol();
        // z21[MC,* ] := Y20[MC,MR] z01[MR,* ] = Y20[MC,MR] (ABL^H aB1)[MR,* ]
        z21_MC_STAR.AlignWith( Y20 );
        Zeros( z21_MC_STAR, A22.Width(), 1 );
        LocalGemv( NORMAL, F(1), Y20, z01_MR_STAR, F(0), z21_MC_STAR );
        // z01[MR,* ] := XBL^H[MR,MC] aB1[MC,* ]
        LocalGemv( ADJOINT, F(1), XBL, aB1_MC_STAR, F(0), z01_MR_STAR );
        // Sum the partial contributions to z01[MR,* ] and scatter the result
        z01_MR_MC.SumScatterFrom( z01_MR_STAR );
        // Redistribute the scattered summation
        z01_MC_STAR.AlignWith( A02 );
        z01_MC_STAR = z01_MR_MC;
        // z21[MR,* ] := z21[MR,* ] - A02^H[MR,MC] z01[MC,* ]
        //             = AB2^H[MR,MC] aB1[MC,* ] - 
        //               A02^H[MR,MC] (XBL^H aB1)[MC,* ]
        LocalGemv( ADJOINT, F(-1), A02, z01_MC_STAR, F(1), z21_MR_STAR );
        // Sum the partial contributions to z21[MR,* ] and scatter the result
        z21_MR_MC.SumScatterFrom( z21_MR_STAR );
        // Redistribute (and rename) the scattered summation
        y21 = z21_MR_MC;
        // Subtract z21 = Y20 ABL^H aB1 from y21
        y21.SumScatterUpdate( F(-1), z21_MC_STAR );
        if( thisIsMyCol )
            Scale( tauQ, y21 );

        //
        // y21 := y21 + Y20 a10^H
        //
        Conjugate( a10 );
        a10_STAR_MR.AlignWith( Y20 );
        a10_STAR_MR = a10;
        Conjugate( a10 );
        // q21[MC,* ] := Y20[MC,MR] a10^H[MR,* ]
        q21_MC_STAR.AlignWith( Y20 );
        Zeros( q21_MC_STAR, A22.Width(), 1 );
        LocalGemv( NORMAL, F(1), Y20, a10_STAR_MR, F(0), q21_MC_STAR );
        // Sum the partial contributions
        q21.AlignWith( y21 );
        q21.SumScatterFrom( q21_MC_STAR );
        if( thisIsMyCol )
            Axpy( F(1), y21, q21 );

        //
        // a12 := conj(a12 - a10 Y20^H - x10 A02)
        //
        Conjugate( x10 );
        x10_STAR_MC.AlignWith( A02 );
        x10_STAR_MC = x10;
        Conjugate( x10 );
        q21_MR_MC.AlignWith( a12 );
        q21_MR_MC = q21;
        // q21[MR,* ] := A02^H[MR,MC] x10^H[MC,* ]
        q21_MR_STAR.AlignWith( A02 );
        Zeros( q21_MR_STAR, A22.Width(), 1 );
        LocalGemv( ADJOINT, F(1), A02, x10_STAR_MC, F(0), q21_MR_STAR );
        // Sum the partial contributions onto q21[MR,MC] = (Y20 a10^H)[MR,MC]
        q21_MR_MC.SumScatterUpdate( F(1), q21_MR_STAR );
        // a12 := conj(a12) - q21^T = conj(a12 - a10 Y20^H - x10 A02)
        Conjugate( a12 );
        if( thisIsMyRow )
        {
            const Int localWidth = a12.LocalWidth();
            F* a12Buffer = a12.Buffer();
            const F* q21Buffer = q21_MR_MC.LockedBuffer();
            const Int a12LDim = a12.LDim();
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
                a12Buffer[jLoc*a12LDim] -= q21Buffer[jLoc];
        }

        // Find tauP, v, and epsilon such that
        //     I - conj(tauP) | 1 | | 1, v^H | | alpha12L | = | epsilon |
        //                    | v |            |  a12R^T  |   |    0    |
        F tauP = 0;
        if( thisIsMyRow )
        {
            tauP = reflector::Row( alpha12L, a12R );
            if( nextIsMyCol )
            {
                tauP1.SetLocal(0,0,tauP);
                // Store epsilon and force | alpha12L | = | 1 |
                //                         |  a12R^T  |   | v |
                epsilon1.SetLocal(0,0,alpha12L.GetLocalRealPart(0,0));
                alpha12L.SetLocal(0,0,F(1));
            }
        }
        mpi::Broadcast( &tauP, 1, alpha11.ColAlign(), g.ColComm() );

        //
        // (Keep in mind that a12 is currently overwritten with its conjugate.
        //  We will use the 'true' value in the following comments.)
        //
        // x21 := conj(tauP) ( A22 a12^H - A2L Y2L^H a12^H - X20 A02 a12^H )
        //
        a12_STAR_MC.AlignWith( Y2L );
        a12_STAR_MR = a12;
        a12_STAR_MC = a12;
        // s21[MC,* ] := A22[MC,MR] a12^H[MR,* ]
        s21_MC_STAR.AlignWith( A22 );
        Zeros( s21_MC_STAR, A22.Height(), 1 );
        LocalGemv( NORMAL, F(1), A22, a12_STAR_MR, F(0), s21_MC_STAR );
        // sB1[MR,* ] := Y2L^H[MR,MC] a12^H[MC,* ]
        sB1_MR_STAR.AlignWith( Y2L );
        Zeros( sB1_MR_STAR, Y2L.Width(), 1 );
        LocalGemv( ADJOINT, F(1), Y2L, a12_STAR_MC, F(0), sB1_MR_STAR );
        // Sum the partial contributions
        sB1_MR_STAR.SumOverCol(); 
        // s21[MC,* ] := s21[MC,* ] - A2L[MC,MR] sB1[MR,* ]
        //             = A22[MC,MR] a12^H[MR,* ] - A2L[MC,MR] sB1[MR,* ]
        // (still needs to be summed within each process row)
        LocalGemv( NORMAL, F(-1), A2L, sB1_MR_STAR, F(1), s21_MC_STAR );
        // s01[MC,* ] := A02[MC,MR] a12^H[MR,* ]
        s01_MC_STAR.AlignWith( A02 );
        Zeros( s01_MC_STAR, k, 1 );
        LocalGemv( NORMAL, F(1), A02, a12_STAR_MR, F(0), s01_MC_STAR );
        // Sum the partial contributions and then redistribute
        s01.SumScatterFrom( s01_MC_STAR ); // TODO: SumScatter to [VC,* ]?
        s01_MR_STAR.AlignWith( X20 );
        s01_MR_STAR = s01;
        // s21[MC,* ] := s21[MC,* ] - X20[MC,MR] s01[MR,* ]
        //             = A22[MC,MR] a12^H[MR,* ] - A2L[MC,MR] sB1[MR,* ]
        //                                       - X20[MC,MR] s01[MR,* ]
        LocalGemv( NORMAL, F(-1), X20, s01_MR_STAR, F(1), s21_MC_STAR );
        // Sum the partial contributions into x21
        x21.SumScatterFrom( s21_MC_STAR );
        Scale( tauP, x21.Matrix() );

        // Undo the in-place conjugation of a12
        Conjugate( a12 );
        Conjugate( a12_STAR_MR );
    }

    // Put back d and e
    auto ATL = View( A, 0, 0, nX, nX );
    auto ATLExpanded = View( A, 0, 0, nX, nX+1 );
    ATL.SetRealPartOfDiagonal( d, 0 );
    ATLExpanded.SetRealPartOfDiagonal( e, 1 );
}

} // namespace bidiag
} // namespace elem

#endif // ifndef ELEM_LAPACK_BIDIAG_PANELU_HPP
