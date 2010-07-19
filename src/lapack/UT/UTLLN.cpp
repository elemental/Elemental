/*
   Copyright (C) 2009-2010 Jack Poulson <jack.poulson@gmail.com>

   This file is part of Elemental.

   Elemental is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   Elemental is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with Elemental.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "elemental/blas_internal.hpp"
#include "elemental/lapack_internal.hpp"
using namespace elemental;
using namespace std;

namespace {

template<typename T>
void
SetDiagonalToOne( int offset, DistMatrix<T,MC,MR>& H )
{
#ifndef RELEASE
    PushCallStack("SetDiagonalToOne");
#endif
    const int localHeight = H.LocalHeight();
    const int localWidth = H.LocalWidth();
    const int r = H.GetGrid().Height();
    const int c = H.GetGrid().Width();
    const int colShift = H.ColShift();
    const int rowShift = H.RowShift();

    for( int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const int i = colShift + iLoc*r + offset;
        if( i % c == rowShift )
        {
            const int jLoc = (i-rowShift) / c;
            if( jLoc < localWidth )
                H.LocalEntry(iLoc,jLoc) = (T)1;
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void 
HalveMainDiagonal( DistMatrix<R,Star,Star>& SInv )
{
#ifndef RELEASE
    PushCallStack("HalveMainDiagonal");
#endif
    for( int j=0; j<SInv.Height(); ++j )
        SInv.LocalEntry(j,j) /= (R)2;
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename R>
void 
FixDiagonal
( const DistMatrix<complex<R>,Star,Star>& t,
        DistMatrix<complex<R>,Star,Star>& SInv )
{
#ifndef RELEASE
    PushCallStack("FixDiagonal");
#endif
    for( int j=0; j<SInv.Height(); ++j )
        SInv.LocalEntry(j,j) = complex<R>(1)/t.LocalEntry(j,0);
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif

} // anonymous namespace

// This routine applies the accumulation of Householder transforms stored 
// in the portion of H below the diagonal marked by 'offset'. It is assumed 
// that the Householder transforms were accumulated left-to-right.

template<typename R>
void
elemental::lapack::internal::UTLLN
( int offset, 
  const DistMatrix<R,MC,MR>& H,
        DistMatrix<R,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::UTLLN");
    if( H.GetGrid() != A.GetGrid() )
        throw logic_error( "H and A must be distributed over the same grid." );
    if( offset > 0 )
        throw logic_error( "Transforms cannot extend above matrix." );
    if( offset < -H.Height() )
        throw logic_error( "Transforms cannot extend below matrix." );
    if( H.Height() != A.Height() )
        throw logic_error
              ( "Height of transforms must equal height of target matrix." );
#endif
    const Grid& g = H.GetGrid();

    // Matrix views    
    DistMatrix<R,MC,MR>
        HTL(g), HTR(g),  H00(g), H01(g), H02(g),  HPan(g), HPanCopy(g),
        HBL(g), HBR(g),  H10(g), H11(g), H12(g),
                         H20(g), H21(g), H22(g);
    DistMatrix<R,MC,MR>
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),  ABottom(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);

    DistMatrix<R,VC,  Star> HPan_VC_Star(g);
    DistMatrix<R,MC,  Star> HPan_MC_Star(g);
    DistMatrix<R,Star,Star> SInv_Star_Star(g);
    DistMatrix<R,Star,MR  > Z_Star_MR(g);
    DistMatrix<R,Star,VR  > Z_Star_VR(g);

    LockedPartitionDownDiagonal
    ( H, HTL, HTR,
         HBL, HBR, 0 );
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    while( HTL.Height() < H.Height() && HTL.Width() < H.Width() )
    {
        LockedRepartitionDownDiagonal
        ( HTL, /**/ HTR,  H00, /**/ H01, H02,
         /*************/ /******************/
               /**/       H10, /**/ H11, H12,
          HBL, /**/ HBR,  H20, /**/ H21, H22 );

        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        HPan.LockedView2x1
        ( H11,
          H21 );

        ABottom.View1x2( ABL, ABR );

        HPan_MC_Star.AlignWith( ABottom );
        Z_Star_MR.AlignWith( ABottom );
        Z_Star_VR.AlignWith( ABottom );
        Z_Star_MR.ResizeTo( HPan.Width(), ABottom.Width() );
        SInv_Star_Star.ResizeTo( HPan.Width(), HPan.Width() );
        //--------------------------------------------------------------------//
        HPanCopy = HPan;
        HPanCopy.MakeTrapezoidal( Left, Lower, offset );
        SetDiagonalToOne( offset, HPanCopy );

        HPan_VC_Star = HPanCopy;
        blas::Herk
        ( Upper, ConjugateTranspose, 
          (R)1, HPan_VC_Star.LockedLocalMatrix(),
          (R)0, SInv_Star_Star.LocalMatrix() );     
        SInv_Star_Star.AllSum();
        HalveMainDiagonal( SInv_Star_Star );

        HPan_MC_Star = HPanCopy;
        blas::internal::LocalGemm
        ( ConjugateTranspose, Normal, 
          (R)1, HPan_MC_Star, ABottom, (R)0, Z_Star_MR );
        Z_Star_VR.SumScatterFrom( Z_Star_MR );
        
        blas::internal::LocalTrsm
        ( Left, Upper, ConjugateTranspose, NonUnit, 
          (R)1, SInv_Star_Star, Z_Star_VR );

        Z_Star_MR = Z_Star_VR;
        blas::internal::LocalGemm
        ( Normal, Normal, (R)-1, HPan_MC_Star, Z_Star_MR, (R)1, ABottom );
        //--------------------------------------------------------------------//
        HPan_MC_Star.FreeAlignments();
        Z_Star_MR.FreeAlignments();
        Z_Star_VR.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( HTL, /**/ HTR,  H00, H01, /**/ H02,
               /**/       H10, H11, /**/ H12,
         /*************/ /******************/
          HBL, /**/ HBR,  H20, H21, /**/ H22 );

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
               /**/       A10, /**/ A11, A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, /**/ A21, A22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename R>
void
elemental::lapack::internal::UTLLN
( int offset, 
  const DistMatrix<complex<R>,MC,MR  >& H,
  const DistMatrix<complex<R>,MD,Star>& t,
        DistMatrix<complex<R>,MC,MR  >& A )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::UTLLN");
    if( H.GetGrid() != t.GetGrid() || t.GetGrid() != A.GetGrid() )
        throw logic_error
              ( "H, t, and A must be distributed over the same grid." );
    if( offset > 0 )
        throw logic_error( "Transforms cannot extend above matrix." );
    if( offset < -H.Height() )
        throw logic_error( "Transforms cannot extend below matrix." );
    if( H.Height() != A.Height() )
        throw logic_error
              ( "Height of transforms must equal height of target matrix." );
    if( t.Height() != H.DiagonalLength( offset ) )
        throw logic_error( "t must be the same length as H's 'offset' diag." );
    if( !t.AlignedWithDiag( H, offset ) )
        throw logic_error( "t must be aligned with H's 'offset' diagonal." );
#endif
    typedef complex<R> C;
    const Grid& g = H.GetGrid();

    // Matrix views    
    DistMatrix<C,MC,MR>
        HTL(g), HTR(g),  H00(g), H01(g), H02(g),  HPan(g), HPanCopy(g),
        HBL(g), HBR(g),  H10(g), H11(g), H12(g),
                         H20(g), H21(g), H22(g);
    DistMatrix<C,MC,MR>
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),  ABottom(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);
    DistMatrix<C,MD,Star>
        tT(g),  t0(g),
        tB(g),  t1(g),
                t2(g);

    DistMatrix<C,VC,  Star> HPan_VC_Star(g);
    DistMatrix<C,MC,  Star> HPan_MC_Star(g);
    DistMatrix<C,Star,Star> t1_Star_Star(g);
    DistMatrix<C,Star,Star> SInv_Star_Star(g);
    DistMatrix<C,Star,MR  > Z_Star_MR(g);
    DistMatrix<C,Star,VR  > Z_Star_VR(g);

    LockedPartitionDownDiagonal
    ( H, HTL, HTR,
         HBL, HBR, 0 );
    LockedPartitionDown
    ( t, tT,
         tB, 0 );
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    while( HTL.Height() < H.Height() && HTL.Width() < H.Width() )
    {
        LockedRepartitionDownDiagonal
        ( HTL, /**/ HTR,  H00, /**/ H01, H02,
         /*************/ /******************/
               /**/       H10, /**/ H11, H12,
          HBL, /**/ HBR,  H20, /**/ H21, H22 );

        LockedRepartitionDown
        ( tT,  t0,
         /**/ /**/
               t1,
          tB,  t2 );

        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        HPan.LockedView2x1
        ( H11,
          H21 );

        ABottom.View1x2( ABL, ABR );

        HPan_MC_Star.AlignWith( ABottom );
        Z_Star_MR.AlignWith( ABottom );
        Z_Star_VR.AlignWith( ABottom );
        Z_Star_MR.ResizeTo( HPan.Width(), ABottom.Width() );
        SInv_Star_Star.ResizeTo( HPan.Width(), HPan.Width() );
        //--------------------------------------------------------------------//
        HPanCopy = HPan;
        HPanCopy.MakeTrapezoidal( Left, Lower, offset );
        SetDiagonalToOne( offset, HPanCopy );

        HPan_VC_Star = HPanCopy;
        blas::Herk
        ( Upper, ConjugateTranspose, 
          (C)1, HPan_VC_Star.LockedLocalMatrix(),
          (C)0, SInv_Star_Star.LocalMatrix() );     
        SInv_Star_Star.AllSum();
        t1_Star_Star = t1;
        FixDiagonal( t1_Star_Star, SInv_Star_Star );

        HPan_MC_Star = HPanCopy;
        blas::internal::LocalGemm
        ( ConjugateTranspose, Normal, 
          (C)1, HPan_MC_Star, ABottom, (C)0, Z_Star_MR );
        Z_Star_VR.SumScatterFrom( Z_Star_MR );
        
        blas::internal::LocalTrsm
        ( Left, Upper, ConjugateTranspose, NonUnit, 
          (C)1, SInv_Star_Star, Z_Star_VR );

        Z_Star_MR = Z_Star_VR;
        blas::internal::LocalGemm
        ( Normal, Normal, (C)-1, HPan_MC_Star, Z_Star_MR, (C)1, ABottom );
        //--------------------------------------------------------------------//
        HPan_MC_Star.FreeAlignments();
        Z_Star_MR.FreeAlignments();
        Z_Star_VR.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( HTL, /**/ HTR,  H00, H01, /**/ H02,
               /**/       H10, H11, /**/ H12,
         /*************/ /******************/
          HBL, /**/ HBR,  H20, H21, /**/ H22 );

        SlideLockedPartitionDown
        ( tT,  t0,
               t1,
         /**/ /**/
          tB,  t2 ); 

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
               /**/       A10, /**/ A11, A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, /**/ A21, A22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif

template void elemental::lapack::internal::UTLLN
( int offset,
  const DistMatrix<float,MC,MR>& H,
        DistMatrix<float,MC,MR>& A );

template void elemental::lapack::internal::UTLLN
( int offset,
  const DistMatrix<double,MC,MR>& H,
        DistMatrix<double,MC,MR>& A );

#ifndef WITHOUT_COMPLEX
template void elemental::lapack::internal::UTLLN
( int offset,
  const DistMatrix<scomplex,MC,MR  >& H,
  const DistMatrix<scomplex,MD,Star>& t,
        DistMatrix<scomplex,MC,MR  >& A );

template void elemental::lapack::internal::UTLLN
( int offset,
  const DistMatrix<dcomplex,MC,MR  >& H,
  const DistMatrix<dcomplex,MD,Star>& t,
        DistMatrix<dcomplex,MC,MR  >& A );
#endif

