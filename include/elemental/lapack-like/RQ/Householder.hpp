/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_RQ_HOUSEHOLDER_HPP
#define LAPACK_RQ_HOUSEHOLDER_HPP

#include "elemental/lapack-like/ApplyPackedReflectors.hpp"
#include "elemental/lapack-like/RQ/PanelHouseholder.hpp"

namespace elem {
namespace rq {

template<typename F> 
inline void
Householder( Matrix<F>& A, Matrix<F>& t )
{
#ifndef RELEASE
    CallStackEntry entry("rq::Householder");
#endif
    t.ResizeTo( std::min(A.Height(),A.Width()), 1 );

    // Matrix views
    Matrix<F>
        ATL, ATR,  A00, A01, A02,  ATopPan, 
        ABL, ABR,  A10, A11, A12,  ABottomPan,
                   A20, A21, A22;
    Matrix<F>
        tT,  t0,
        tB,  t1,
             t2;

    PartitionUpRightDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    PartitionUp
    ( t, tT,
         tB, 0 );
    while( ABR.Height() < A.Height() && ABR.Width() < A.Width() )
    {
        RepartitionUpDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );

        RepartitionUp
        ( tT,  t0,
               t1,
         /**/ /**/
          tB,  t2 );

        View1x2( ATopPan, A00, A01 );
        View1x2( ABottomPan, A10, A11 );

        //--------------------------------------------------------------------//
        PanelHouseholder( ABottomPan, t1 );
        ApplyPackedReflectors
        ( RIGHT, LOWER, HORIZONTAL, BACKWARD, CONJUGATED, 
          0, ABottomPan, t1, ATopPan ); 
        //--------------------------------------------------------------------//

        SlidePartitionUp
        ( tT,  t0,
         /**/ /**/
               t1,
          tB,  t2 );

        SlidePartitionUpDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );
    }
}

template<typename F> 
inline void
Householder( Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("rq::Householder");
#endif
    Matrix<F> t;
    Householder( A, t );
}

template<typename F> 
inline void
Householder( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& t )
{
#ifndef RELEASE
    CallStackEntry entry("rq::Householder");
    if( A.Grid() != t.Grid() )
        throw std::logic_error("{A,s} must be distributed over the same grid");
#endif
    if( t.Viewing() )
    {
        if( !t.AlignedWithDiagonal( A ) ) 
            throw std::logic_error("t was not aligned with A");
        if( t.Height() != std::min(A.Height(),A.Width()) || t.Width() != 1 )
            throw std::logic_error("t was not the appropriate shape");
    }
    else
    {
        t.AlignWithDiagonal( A );
        t.ResizeTo( std::min(A.Height(),A.Width()), 1 );
    }

    // Matrix views
    const Grid& g = A.Grid();
    DistMatrix<F>
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),  ATopPan(g), 
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),  ABottomPan(g),
                         A20(g), A21(g), A22(g);
    DistMatrix<F,MD,STAR>
        tT(g),  t0(g),
        tB(g),  t1(g),
                t2(g);

    PartitionUpLeftDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    PartitionUp
    ( t, tT,
         tB, 0 );
    while( ABR.Height() < A.Height() && ABR.Width() < A.Width() )
    {
        RepartitionUpDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );

        RepartitionUp
        ( tT,  t0,
               t1,
         /**/ /**/
          tB,  t2 );

        View1x2( ATopPan, A00, A01 );
        View1x2( ABottomPan, A10, A11 );

        //--------------------------------------------------------------------//
        PanelHouseholder( ABottomPan, t1 );
        ApplyPackedReflectors
        ( RIGHT, LOWER, HORIZONTAL, BACKWARD, CONJUGATED, 
          0, ABottomPan, t1, ATopPan );
        //--------------------------------------------------------------------//

        SlidePartitionUp
        ( tT,  t0,
         /**/ /**/
               t1,
          tB,  t2 );

        SlidePartitionUpDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );
    }
}

template<typename F> 
inline void
Householder( DistMatrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("rq::Householder");
#endif
    DistMatrix<F,MD,STAR> t(A.Grid());
    Householder( A, t );
}

} // namespace rq
} // namespace elem

#endif // ifndef LAPACK_RQ_HOUSEHOLDER_HPP
