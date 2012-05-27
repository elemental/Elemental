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

template<typename T>
inline void
internal::Syr2kUT
( T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("internal::Syr2kUT");
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
        throw std::logic_error
        ("{A,B,C} must be distributed over the same grid");
    if( A.Width() != C.Height() || 
        A.Width() != C.Width()  ||
        B.Width() != C.Height() ||
        B.Width() != C.Width()  ||
        A.Height() != B.Height()  )
    {
        std::ostringstream msg;
        msg << "Nonconformal Syr2kUT:\n"
            << "  A ~ " << A.Height() << " x " << A.Width() << "\n"
            << "  B ~ " << B.Height() << " x " << B.Width() << "\n"
            << "  C ~ " << C.Height() << " x " << C.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
#endif
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<T,MC,MR> AT(g),  A0(g),
                        AB(g),  A1(g),
                                A2(g);

    DistMatrix<T,MC,MR> BT(g),  B0(g),
                        BB(g),  B1(g),
                                B2(g);

    // Temporary distributions
    DistMatrix<T,MR,  STAR> A1Trans_MR_STAR(g);
    DistMatrix<T,MR,  STAR> B1Trans_MR_STAR(g);
    DistMatrix<T,STAR,VR  > A1_STAR_VR(g);
    DistMatrix<T,STAR,VR  > B1_STAR_VR(g);
    DistMatrix<T,STAR,MC  > A1_STAR_MC(g);
    DistMatrix<T,STAR,MC  > B1_STAR_MC(g);

    // Start the algorithm
    ScaleTrapezoid( beta, LEFT, UPPER, 0, C );
    LockedPartitionDown
    ( A, AT, 
         AB, 0 );
    LockedPartitionDown
    ( B, BT,
         BB, 0 );
    while( AB.Height() > 0 )
    {
        LockedRepartitionDown
        ( AT,  A0,
         /**/ /**/
               A1,
          AB,  A2 );

        LockedRepartitionDown
        ( BT,  B0,
         /**/ /**/
               B1,
          BB,  B2 );

        A1Trans_MR_STAR.AlignWith( C );
        B1Trans_MR_STAR.AlignWith( C );
        A1_STAR_MC.AlignWith( C );
        B1_STAR_MC.AlignWith( C );
        //--------------------------------------------------------------------//
        A1Trans_MR_STAR.TransposeFrom( A1 );
        A1_STAR_VR.TransposeFrom( A1Trans_MR_STAR );
        A1_STAR_MC = A1_STAR_VR;

        B1Trans_MR_STAR.TransposeFrom( B1 );
        B1_STAR_VR.TransposeFrom( B1Trans_MR_STAR );
        B1_STAR_MC = B1_STAR_VR;

        internal::LocalTrr2k
        ( UPPER, TRANSPOSE, TRANSPOSE, TRANSPOSE, TRANSPOSE, 
          alpha, A1_STAR_MC, B1Trans_MR_STAR,
                 B1_STAR_MC, A1Trans_MR_STAR,
          (T)1,  C );
        //--------------------------------------------------------------------//
        A1Trans_MR_STAR.FreeAlignments();
        B1Trans_MR_STAR.FreeAlignments();
        A1_STAR_MC.FreeAlignments();
        B1_STAR_MC.FreeAlignments();

        SlideLockedPartitionDown
        ( AT,  A0,
               A1,
         /**/ /**/
          AB,  A2 );

        SlideLockedPartitionDown
        ( BT,  B0,
               B1,
         /**/ /**/
          BB,  B2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
