/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

#include "./LQ/Panel.hpp"

namespace elem {

// On exit, the lower triangle of A is overwritten by L, and the Householder
// transforms that determine Q are stored above the diagonal of A with an 
// implicit one on the diagonal. 
//
// In the complex case, the column-vector t stores the unit-magnitude complex 
// rotations that map the norms of the implicit Householder vectors to their
// coefficient:  
//                psi_j = 2 tau_j / ( u_j^H u_j ),
// where tau_j is the j'th entry of t and u_j is the j'th unscaled Householder
// reflector.

template<typename Real> 
inline void
LQ( Matrix<Real>& A )
{
#ifndef RELEASE
    PushCallStack("LQ");
#endif
    if( IsComplex<Real>::val )
        throw std::logic_error("Called real routine with complex datatype");

    // Matrix views
    Matrix<Real>
        ATL, ATR,  A00, A01, A02,  ATopPan, ABottomPan,
        ABL, ABR,  A10, A11, A12,
                   A20, A21, A22;

    PartitionDownLeftDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    while( ATL.Height() < A.Height() && ATL.Width() < A.Width() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        View1x2( ATopPan, A11, A12 );
        View1x2( ABottomPan, A21, A22 );

        //--------------------------------------------------------------------//
        internal::PanelLQ( ATopPan );
        ApplyPackedReflectors
        ( RIGHT, UPPER, HORIZONTAL, FORWARD, 0, ATopPan, ABottomPan );
        //--------------------------------------------------------------------//

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Real> 
inline void
LQ( DistMatrix<Real>& A )
{
#ifndef RELEASE
    PushCallStack("LQ");
#endif
    if( IsComplex<Real>::val )
        throw std::logic_error("Called real routine with complex datatype");
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<Real>
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),  ATopPan(g), ABottomPan(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);

    PartitionDownLeftDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    while( ATL.Height() < A.Height() && ATL.Width() < A.Width() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        View1x2( ATopPan, A11, A12 );
        View1x2( ABottomPan, A21, A22 );

        //--------------------------------------------------------------------//
        internal::PanelLQ( ATopPan );
        ApplyPackedReflectors
        ( RIGHT, UPPER, HORIZONTAL, FORWARD, 0, ATopPan, ABottomPan );
        //--------------------------------------------------------------------//

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Real> 
inline void
LQ( Matrix<Complex<Real> >& A, 
    Matrix<Complex<Real> >& t )
{
#ifndef RELEASE
    PushCallStack("LQ");
#endif
    typedef Complex<Real> C;
    t.ResizeTo( std::min(A.Height(),A.Width()), 1 );

    // Matrix views
    Matrix<C>
        ATL, ATR,  A00, A01, A02,  ATopPan, ABottomPan,
        ABL, ABR,  A10, A11, A12,
                   A20, A21, A22;
    Matrix<C>
        tT,  t0,
        tB,  t1,
             t2;

    PartitionDownLeftDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    PartitionDown
    ( t, tT,
         tB, 0 );
    while( ATL.Height() < A.Height() && ATL.Width() < A.Width() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        RepartitionDown
        ( tT,  t0,
         /**/ /**/
               t1,
          tB,  t2 );

        View1x2( ATopPan, A11, A12 );
        View1x2( ABottomPan, A21, A22 );

        //--------------------------------------------------------------------//
        internal::PanelLQ( ATopPan, t1 );
        ApplyPackedReflectors
        ( RIGHT, UPPER, HORIZONTAL, FORWARD, CONJUGATED,
          0, ATopPan, t1, ABottomPan );
        //--------------------------------------------------------------------//

        SlidePartitionDown
        ( tT,  t0,
               t1,
         /**/ /**/
          tB,  t2 );

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Real> 
inline void
LQ( DistMatrix<Complex<Real> >& A, 
    DistMatrix<Complex<Real>,MD,STAR>& t )
{
#ifndef RELEASE
    PushCallStack("LQ");
    if( A.Grid() != t.Grid() )
        throw std::logic_error("{A,t} must be distributed over the same grid");
#endif
    typedef Complex<Real> C;
    const Grid& g = A.Grid();
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
    DistMatrix<C>
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),  ATopPan(g), ABottomPan(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);
    DistMatrix<C,MD,STAR>
        tT(g),  t0(g),
        tB(g),  t1(g),
                t2(g);

    PartitionDownLeftDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    PartitionDown
    ( t, tT,
         tB, 0 );
    while( ATL.Height() < A.Height() && ATL.Width() < A.Width() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        RepartitionDown
        ( tT,  t0,
         /**/ /**/
               t1,
          tB,  t2 );

        View1x2( ATopPan, A11, A12 );
        View1x2( ABottomPan, A21, A22 );

        //--------------------------------------------------------------------//
        internal::PanelLQ( ATopPan, t1 );
        ApplyPackedReflectors
        ( RIGHT, UPPER, HORIZONTAL, FORWARD, CONJUGATED,
          0, ATopPan, t1, ABottomPan );
        //--------------------------------------------------------------------//

        SlidePartitionDown
        ( tT,  t0,
               t1,
         /**/ /**/
          tB,  t2 );

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
