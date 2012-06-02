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
SytrmmLUnblocked( Matrix<T>& L )
{
#ifndef RELEASE
    PushCallStack("internal::SytrmmLUnblocked");
#endif
    // Matrix views
    Matrix<T>
        LTL, LTR,  L00, l01,      L02,
        LBL, LBR,  l10, lambda11, l12,
                   L20, l21,      L22;

    PushBlocksizeStack( 1 );
    PartitionDownDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    while( LTL.Height() < L.Height() && LTL.Width() < L.Height() )
    {
        RepartitionDownDiagonal 
        ( LTL, /**/ LTR,  L00, /**/ l01,      L02,
         /*************/ /***********************/
               /**/       l10, /**/ lambda11, l12,
          LBL, /**/ LBR,  L20, /**/ l21,      L22 );

        //--------------------------------------------------------------------//
        const T lambda = lambda11.Get(0,0);

        if( L22.Height() > 0 )
        {
            const T gamma = Dotu( l21, l21 );
            lambda11.Set(0,0,lambda*lambda+gamma);
            Gemv( TRANSPOSE, (T)1, L20, l21, lambda, l10 );
        }
        else
        {
            Scale( lambda, l10 );
        }
        //--------------------------------------------------------------------//

        SlidePartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, l01,      /**/ L02,
               /**/       l10, lambda11, /**/ l12,
         /*************/ /***********************/
          LBL, /**/ LBR,  L20, l21,      /**/ L22 );
    }
    PopBlocksizeStack();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
SytrmmUUnblocked( Matrix<T>& U )
{
#ifndef RELEASE
    PushCallStack("internal::SytrmmUUnblocked");
#endif
    // Matrix views
    Matrix<T>
        UTL, UTR,  U00, u01,       U02,
        UBL, UBR,  u10, upsilon11, u12,
                   U20, u21,       U22;

    PushBlocksizeStack( 1 );
    PartitionDownDiagonal
    ( U, UTL, UTR,
         UBL, UBR, 0 );
    while( UTL.Height() < U.Height() && UTL.Width() < U.Height() )
    {
        RepartitionDownDiagonal 
        ( UTL, /**/ UTR,  U00, /**/ u01,       U02,
         /*************/ /************************/
               /**/       u10, /**/ upsilon11, u12,
          UBL, /**/ UBR,  U20, /**/ u21,       U22 );

        //--------------------------------------------------------------------//
        const T upsilon = upsilon11.Get(0,0);

        if( U22.Height() > 0 )
        {
            const T gamma = Dotu( u12, u12 );
            upsilon11.Set(0,0,upsilon*upsilon+gamma);
            Gemv( NORMAL, (T)1, U02, u12, upsilon, u01 );
        }
        else
        {
            Scale( upsilon, u01 );
        }
        //--------------------------------------------------------------------//

        SlidePartitionDownDiagonal
        ( UTL, /**/ UTR,  U00, u01,       /**/ U02,
               /**/       u10, upsilon11, /**/ u12,
         /*************/ /***********************/
          UBL, /**/ UBR,  U20, u21,      /**/ U22 );
    }
    PopBlocksizeStack();
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal
} // namespace elem
