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
namespace internal {

template<typename T>
inline void
Syr2kLN
( T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("internal::Syr2kLN");
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
        throw std::logic_error
        ("{A,B,C} must be distributed over the same grid");
    if( A.Height() != C.Height() || A.Height() != C.Width() ||
        B.Height() != C.Height() || B.Height() != C.Width() ||
        A.Width() != B.Width()                                 )
    {
        std::ostringstream msg;
        msg << "Nonconformal Syr2kLN:\n"
            << "  A ~ " << A.Height() << " x " << A.Width() << "\n"
            << "  B ~ " << B.Height() << " x " << B.Width() << "\n"
            << "  C ~ " << C.Height() << " x " << C.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
#endif
    const Grid& g = A.Grid();

    // Matrix views 
    DistMatrix<T,MC,MR> AL(g), AR(g),
                        A0(g), A1(g), A2(g);

    DistMatrix<T,MC,MR> BL(g), BR(g),
                        B0(g), B1(g), B2(g);

    // Temporary distributions
    DistMatrix<T,MC,  STAR> A1_MC_STAR(g);
    DistMatrix<T,MC,  STAR> B1_MC_STAR(g);
    DistMatrix<T,VR,  STAR> A1_VR_STAR(g);
    DistMatrix<T,VR,  STAR> B1_VR_STAR(g);
    DistMatrix<T,STAR,MR  > A1Trans_STAR_MR(g);
    DistMatrix<T,STAR,MR  > B1Trans_STAR_MR(g);

    // Start the algorithm
    ScaleTrapezoid( beta, LEFT, LOWER, 0, C );
    LockedPartitionRight( A, AL, AR, 0 );
    LockedPartitionRight( B, BL, BR, 0 );
    while( AR.Width() > 0 )
    {
        LockedRepartitionRight
        ( AL, /**/ AR,
          A0, /**/ A1, A2 );

        LockedRepartitionRight
        ( BL, /**/ BR,
          B0, /**/ B1, B2 );

        A1_MC_STAR.AlignWith( C );
        B1_MC_STAR.AlignWith( C );
        A1_VR_STAR.AlignWith( C );
        B1_VR_STAR.AlignWith( C );
        A1Trans_STAR_MR.AlignWith( C );
        B1Trans_STAR_MR.AlignWith( C );
        //--------------------------------------------------------------------//
        A1_VR_STAR = A1_MC_STAR = A1;
        A1Trans_STAR_MR.TransposeFrom( A1_VR_STAR );

        B1_VR_STAR = B1_MC_STAR = B1;
        B1Trans_STAR_MR.TransposeFrom( B1_VR_STAR );

        LocalTrr2k
        ( LOWER, 
          alpha, A1_MC_STAR, B1Trans_STAR_MR,
                 B1_MC_STAR, A1Trans_STAR_MR,
          (T)1,  C );
        //--------------------------------------------------------------------//
        A1_MC_STAR.FreeAlignments();
        B1_MC_STAR.FreeAlignments();
        A1_VR_STAR.FreeAlignments();
        B1_VR_STAR.FreeAlignments();
        A1Trans_STAR_MR.FreeAlignments();
        B1Trans_STAR_MR.FreeAlignments();

        SlideLockedPartitionRight
        ( AL,     /**/ AR,
          A0, A1, /**/ A2 );

        SlideLockedPartitionRight
        ( BL,     /**/ BR,
          B0, B1, /**/ B2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal
} // namespace elem
