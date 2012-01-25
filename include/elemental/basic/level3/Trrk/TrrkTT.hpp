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

// Distributed C := alpha A^{T/H} B^{T/H} + beta C
template<typename T>
inline void
internal::TrrkTT
( UpperOrLower uplo,
  Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("internal::TrrkTN");
    if( C.Height() != C.Width() ||
        A.Width() != C.Height() || 
        B.Height() != C.Width() ||
        A.Height() != B.Width() )
        throw std::logic_error("Nonconformal TrrkTN");
    if( orientationOfA == NORMAL || orientationOfB == NORMAL )
        throw std::logic_error("Orientations must be TRANSPOSE or ADJOINT");
#endif
    const Grid& g = C.Grid();

    DistMatrix<T,MC,MR> AT(g),  A0(g),
                        AB(g),  A1(g),
                                A2(g);
    DistMatrix<T,MC,MR> BL(g), BR(g),
                        B0(g), B1(g), B2(g);

    DistMatrix<T,STAR,MC> A1_STAR_MC(g);
    DistMatrix<T,VR,STAR> B1_VR_STAR(g);
    DistMatrix<T,STAR,MR> B1AdjOrTrans_STAR_MR(g);

    LockedPartitionDown
    ( A, AT,
         AB, 0 );
    LockedPartitionRight( B, BL, BR, 0 );
    while( AT.Width() < A.Height() )
    {
        LockedRepartitionDown
        ( AT,  A0,
         /**/ /**/
               A1,
          AB,  A2 );
        LockedRepartitionRight
        ( BL, /**/ BR,
          B0, /**/ B1, B2 );

        A1_STAR_MC.AlignWith( C );
        B1_VR_STAR.AlignWith( C );
        B1AdjOrTrans_STAR_MR.AlignWith( C );
        //--------------------------------------------------------------------//
        A1_STAR_MC = A1;
        B1_VR_STAR = B1;
        if( orientationOfB == ADJOINT )
            B1AdjOrTrans_STAR_MR.AdjointFrom( B1_VR_STAR );
        else
            B1AdjOrTrans_STAR_MR.TransposeFrom( B1_VR_STAR );
        internal::LocalTrrk
        ( uplo, orientationOfA,
          alpha, A1_STAR_MC, B1AdjOrTrans_STAR_MR, beta, C );
        //--------------------------------------------------------------------//
        B1AdjOrTrans_STAR_MR.FreeAlignments();
        B1_VR_STAR.FreeAlignments();
        A1_STAR_MC.FreeAlignments();

        SlideLockedPartitionRight
        ( BL,     /**/ BR,
          B0, B1, /**/ B2 );
        SlideLockedPartitionDown
        ( AT,  A0,
               A1,
         /**/ /**/
          AB,  A2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
