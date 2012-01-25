/*
   Copyright (c) 2009-2012, Jack Poulson
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

namespace elem {

// Distributed E := alpha (A B^{T/H} + C^{T/H} D^{T/H}) + beta E
template<typename T>
inline void
internal::Trr2kNTTT
( UpperOrLower uplo,
  Orientation orientationOfB, 
  Orientation orientationOfC, Orientation orientationOfD,
  T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
           const DistMatrix<T,MC,MR>& C, const DistMatrix<T,MC,MR>& D,
  T beta,        DistMatrix<T,MC,MR>& E )
{
#ifndef RELEASE
    PushCallStack("internal::Trr2kNTTT");
    if( E.Height() != E.Width()  || A.Width()  != C.Height() ||
        A.Height() != E.Height() || C.Width()  != E.Height() ||
        B.Height() != E.Width()  || D.Height() != E.Width()  ||
        A.Width()  != B.Width()  || C.Height() != D.Width() )
        throw std::logic_error("Nonconformal Trr2kNTTT");
#endif
    const Grid& g = E.Grid();

    DistMatrix<T,MC,MR> AL(g), AR(g),
                        A0(g), A1(g), A2(g);
    DistMatrix<T,MC,MR> BL(g), BR(g),
                        B0(g), B1(g), B2(g);

    DistMatrix<T,MC,MR> CT(g),  C0(g),
                        CB(g),  C1(g),
                                C2(g);
    DistMatrix<T,MC,MR> DL(g), DR(g),
                        D0(g), D1(g), D2(g);

    DistMatrix<T,MC,  STAR> A1_MC_STAR(g);
    DistMatrix<T,VR,  STAR> B1_VR_STAR(g);
    DistMatrix<T,STAR,MR  > B1AdjOrTrans_STAR_MR(g);
    DistMatrix<T,STAR,MC  > C1_STAR_MC(g);
    DistMatrix<T,VR,  STAR> D1_VR_STAR(g);
    DistMatrix<T,STAR,MR  > D1AdjOrTrans_STAR_MR(g);

    LockedPartitionRight( A, AL, AR, 0 );
    LockedPartitionRight( B, BL, BR, 0 );
    LockedPartitionDown
    ( C, CT,
         CB, 0 );
    LockedPartitionRight( D, DL, DR, 0 );
    while( AL.Width() < A.Width() )
    {
        LockedRepartitionRight
        ( AL, /**/ AR,
          A0, /**/ A1, A2 );
        LockedRepartitionRight
        ( BL, /**/ BR,
          B0, /**/ B1, B2 );
        LockedRepartitionDown
        ( CT,  C0,
         /**/ /**/
               C1,
          CB,  C2 );
        LockedRepartitionRight
        ( DL, /**/ DR,
          D0, /**/ D1, D2 );

        A1_MC_STAR.AlignWith( E );
        B1_VR_STAR.AlignWith( E );
        B1AdjOrTrans_STAR_MR.AlignWith( E );
        C1_STAR_MC.AlignWith( E );
        D1_VR_STAR.AlignWith( E );
        D1AdjOrTrans_STAR_MR.AlignWith( E );
        //--------------------------------------------------------------------//
        A1_MC_STAR = A1;
        C1_STAR_MC = C1;
        B1_VR_STAR = B1;
        D1_VR_STAR = D1;
        if( orientationOfB == ADJOINT )
            B1AdjOrTrans_STAR_MR.AdjointFrom( B1_VR_STAR );
        else
            B1AdjOrTrans_STAR_MR.TransposeFrom( B1_VR_STAR );
        if( orientationOfD == ADJOINT )
            D1AdjOrTrans_STAR_MR.AdjointFrom( D1_VR_STAR );
        else
            D1AdjOrTrans_STAR_MR.TransposeFrom( D1_VR_STAR );
        internal::LocalTrr2k
        ( uplo, orientationOfC,
          alpha, A1_MC_STAR, B1AdjOrTrans_STAR_MR, 
                 C1_STAR_MC, D1AdjOrTrans_STAR_MR,
          beta,  E );
        //--------------------------------------------------------------------//
        D1AdjOrTrans_STAR_MR.FreeAlignments();
        D1_VR_STAR.FreeAlignments();
        C1_STAR_MC.FreeAlignments();
        B1AdjOrTrans_STAR_MR.FreeAlignments();
        B1_VR_STAR.FreeAlignments();
        A1_MC_STAR.FreeAlignments();

        SlideLockedPartitionRight
        ( DL,     /**/ DR,
          D0, D1, /**/ D2 );
        SlideLockedPartitionDown
        ( CT,  C0,
               C1,
         /**/ /**/
          CB,  C2 );
        SlideLockedPartitionRight
        ( BL,     /**/ BR,
          B0, B1, /**/ B2 );
        SlideLockedPartitionRight
        ( AL,     /**/ AR,
          A0, A1, /**/ A2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
