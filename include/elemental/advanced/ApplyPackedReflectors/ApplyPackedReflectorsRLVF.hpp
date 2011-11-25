/*
   Copyright (c) 2009-2011, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/

template<typename R> // representation of a real number
inline void
elemental::advanced::internal::ApplyPackedReflectorsRLVF
( int offset, 
  const DistMatrix<R,MC,MR>& H,
        DistMatrix<R,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("advanced::internal::ApplyPackedReflectorsRLVF");
    if( H.Grid() != A.Grid() )
        throw std::logic_error("{H,A} must be distributed over the same grid");
    if( offset > 0 )
        throw std::logic_error("Transforms cannot extend above matrix");
    if( offset < -H.Height() )
        throw std::logic_error("Transforms cannot extend below matrix");
    if( H.Height() != A.Width() )
        throw std::logic_error
              ("Height of transforms must equal width of target matrix");
#endif
    const Grid& g = H.Grid();

    // Matrix views    
    DistMatrix<R,MC,MR>
        HTL(g), HTR(g),  H00(g), H01(g), H02(g),  HPan(g), HPanCopy(g),
        HBL(g), HBR(g),  H10(g), H11(g), H12(g),
                         H20(g), H21(g), H22(g);
    DistMatrix<R,MC,MR>
        AL(g), AR(g),
        A0(g), A1(g), A2(g);

    DistMatrix<R,VC,  STAR> HPan_VC_STAR(g);
    DistMatrix<R,MR,  STAR> HPan_MR_STAR(g);
    DistMatrix<R,STAR,STAR> SInv_STAR_STAR(g);
    DistMatrix<R,STAR,MC  > Z_STAR_MC(g);
    DistMatrix<R,STAR,VC  > Z_STAR_VC(g);

    LockedPartitionDownDiagonal
    ( H, HTL, HTR,
         HBL, HBR, 0 );
    PartitionRight
    ( A, AL, AR, 0 );
    while( HTL.Height() < H.Height() && HTL.Width() < H.Width() )
    {
        LockedRepartitionDownDiagonal
        ( HTL, /**/ HTR,  H00, /**/ H01, H02,
         /*************/ /******************/
               /**/       H10, /**/ H11, H12,
          HBL, /**/ HBR,  H20, /**/ H21, H22 );

        RepartitionRight
        ( AL, /**/ AR,
          A0, /**/ A1, A2 );

        const int HPanHeight = H11.Height() + H21.Height();
        const int HPanWidth = 
            std::min( H11.Width(), std::max(HPanHeight+offset,0) );
        HPan.LockedView( H, H00.Height(), H00.Width(), HPanHeight, HPanWidth );

        HPan_MR_STAR.AlignWith( AR );
        Z_STAR_MC.AlignWith( AR );
        Z_STAR_VC.AlignWith( AR );
        Z_STAR_MC.ResizeTo( HPanWidth, AR.Height() );
        SInv_STAR_STAR.ResizeTo( HPanWidth, HPanWidth );
        //--------------------------------------------------------------------//
        HPanCopy = HPan;
        HPanCopy.MakeTrapezoidal( LEFT, LOWER, offset );
        SetDiagonalToOne( LEFT, offset, HPanCopy );

        HPan_VC_STAR = HPanCopy;
        basic::Syrk
        ( UPPER, TRANSPOSE, 
          (R)1, HPan_VC_STAR.LockedLocalMatrix(),
          (R)0, SInv_STAR_STAR.LocalMatrix() );     
        SInv_STAR_STAR.SumOverGrid();
        HalveMainDiagonal( SInv_STAR_STAR );

        HPan_MR_STAR = HPan_VC_STAR;
        basic::internal::LocalGemm
        ( ADJOINT, ADJOINT, (R)1, HPan_MR_STAR, AR, (R)0, Z_STAR_MC );
        Z_STAR_VC.SumScatterFrom( Z_STAR_MC );
        
        basic::internal::LocalTrsm
        ( LEFT, UPPER, ADJOINT, NON_UNIT, (R)1, SInv_STAR_STAR, Z_STAR_VC );

        Z_STAR_MC = Z_STAR_VC;
        basic::internal::LocalGemm
        ( ADJOINT, TRANSPOSE, (R)-1, Z_STAR_MC, HPan_MR_STAR, (R)1, AR );
        //--------------------------------------------------------------------//
        HPan_MR_STAR.FreeAlignments();
        Z_STAR_MC.FreeAlignments();
        Z_STAR_VC.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( HTL, /**/ HTR,  H00, H01, /**/ H02,
               /**/       H10, H11, /**/ H12,
         /*************/ /******************/
          HBL, /**/ HBR,  H20, H21, /**/ H22 );

        SlidePartitionRight
        ( AL,     /**/ AR,
          A0, A1, /**/ A2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename R> // representation of a real number
void
elemental::advanced::internal::ApplyPackedReflectorsRLVF
( Conjugation conjugation, int offset, 
  const DistMatrix<std::complex<R>,MC,MR  >& H,
  const DistMatrix<std::complex<R>,MD,STAR>& t,
        DistMatrix<std::complex<R>,MC,MR  >& A )
{
#ifndef RELEASE
    PushCallStack("advanced::internal::ApplyPackedReflectorsRLVF");
    if( H.Grid() != t.Grid() || t.Grid() != A.Grid() )
        throw std::logic_error
              ("{H,t,A} must be distributed over the same grid");
    if( offset > 0 )
        throw std::logic_error("Transforms cannot extend above matrix");
    if( offset < -H.Height() )
        throw std::logic_error("Transforms cannot extend below matrix");
    if( H.Height() != A.Width() )
        throw std::logic_error
              ("Height of transforms must equal width of target matrix");
    if( t.Height() != H.DiagonalLength( offset ) )
        throw std::logic_error("t must be the same length as H's offset diag");
    if( !t.AlignedWithDiagonal( H, offset ) )
        throw std::logic_error("t must be aligned with H's offset diagonal");
#endif
    typedef std::complex<R> C;
    const Grid& g = H.Grid();

    // Matrix views    
    DistMatrix<C,MC,MR>
        HTL(g), HTR(g),  H00(g), H01(g), H02(g),  HPan(g), HPanCopy(g),
        HBL(g), HBR(g),  H10(g), H11(g), H12(g),
                         H20(g), H21(g), H22(g);
    DistMatrix<C,MC,MR>
        AL(g), AR(g),
        A0(g), A1(g), A2(g);
    DistMatrix<C,MD,STAR>
        tT(g),  t0(g),
        tB(g),  t1(g),
                t2(g);

    DistMatrix<C,VC,  STAR> HPan_VC_STAR(g);
    DistMatrix<C,MR,  STAR> HPan_MR_STAR(g);
    DistMatrix<C,STAR,STAR> t1_STAR_STAR(g);
    DistMatrix<C,STAR,STAR> SInv_STAR_STAR(g);
    DistMatrix<C,STAR,MC  > Z_STAR_MC(g);
    DistMatrix<C,STAR,VC  > Z_STAR_VC(g);

    LockedPartitionDownDiagonal
    ( H, HTL, HTR,
         HBL, HBR, 0 );
    LockedPartitionDown
    ( t, tT,
         tB, 0 );
    PartitionRight
    ( A, AL, AR, 0 );
    while( HTL.Height() < H.Height() && HTL.Width() < H.Width() )
    {
        LockedRepartitionDownDiagonal
        ( HTL, /**/ HTR,  H00, /**/ H01, H02,
         /*************/ /******************/
               /**/       H10, /**/ H11, H12,
          HBL, /**/ HBR,  H20, /**/ H21, H22 );

        const int HPanHeight = H11.Height() + H21.Height();
        const int HPanWidth = 
            std::min( H11.Width(), std::max(HPanHeight+offset,0) );
        HPan.LockedView( H, H00.Height(), H00.Width(), HPanHeight, HPanWidth );

        LockedRepartitionDown
        ( tT,  t0,
         /**/ /**/
               t1,
          tB,  t2, HPanWidth );

        RepartitionRight
        ( AL, /**/ AR,
          A0, /**/ A1, A2 );

        HPan_MR_STAR.AlignWith( AR );
        Z_STAR_MC.AlignWith( AR );
        Z_STAR_VC.AlignWith( AR );
        Z_STAR_MC.ResizeTo( HPanWidth, AR.Height() );
        SInv_STAR_STAR.ResizeTo( HPanWidth, HPanWidth );
        //--------------------------------------------------------------------//
        HPanCopy = HPan;
        HPanCopy.MakeTrapezoidal( LEFT, LOWER, offset );
        SetDiagonalToOne( LEFT, offset, HPanCopy );

        HPan_VC_STAR = HPanCopy;
        basic::Herk
        ( UPPER, ADJOINT, 
          (C)1, HPan_VC_STAR.LockedLocalMatrix(),
          (C)0, SInv_STAR_STAR.LocalMatrix() );     
        SInv_STAR_STAR.SumOverGrid();
        t1_STAR_STAR = t1;
        FixDiagonal( conjugation, t1_STAR_STAR, SInv_STAR_STAR );

        HPan_MR_STAR = HPan_VC_STAR;
        basic::internal::LocalGemm
        ( ADJOINT, ADJOINT, (C)1, HPan_MR_STAR, AR, (C)0, Z_STAR_MC );
        Z_STAR_VC.SumScatterFrom( Z_STAR_MC );
        
        basic::internal::LocalTrsm
        ( LEFT, UPPER, ADJOINT, NON_UNIT, (C)1, SInv_STAR_STAR, Z_STAR_VC );

        Z_STAR_MC = Z_STAR_VC;
        basic::internal::LocalGemm
        ( ADJOINT, ADJOINT, (C)-1, Z_STAR_MC, HPan_MR_STAR, (C)1, AR );
        //--------------------------------------------------------------------//
        HPan_MR_STAR.FreeAlignments();
        Z_STAR_MC.FreeAlignments();
        Z_STAR_VC.FreeAlignments();

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

        SlidePartitionRight
        ( AL,     /**/ AR,
          A0, A1, /**/ A2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX
