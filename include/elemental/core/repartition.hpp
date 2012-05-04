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
#ifndef ELEMENTAL_REPARTITION_HPP
#define ELEMENTAL_REPARTITION_HPP 1

#include "elemental/core/dist_matrix.hpp"

namespace elem {

// To make our life easier. Undef'd at the bottom of the header
#define M  Matrix<T,Int>
#define DM DistMatrix<T,U,V,Int>

//
// RepartitionUp
//

template<typename T,typename Int>
void RepartitionUp
( M& AT, M& A0,
         M& A1,
  M& AB, M& A2, Int bsize=Blocksize() );

template<typename T,Distribution U,Distribution V,typename Int>
void RepartitionUp
( DM& AT, DM& A0,
          DM& A1,
  DM& AB, DM& A2, Int bsize=Blocksize() );

template<typename T,typename Int>
void LockedRepartitionUp
( const M& AT, M& A0,
               M& A1,
  const M& AB, M& A2, Int bsize=Blocksize() );

template<typename T,Distribution U,Distribution V,typename Int>
void LockedRepartitionUp
( const DM& AT, DM& A0,
                DM& A1,
  const DM& AB, DM& A2, Int bsize=Blocksize() );

//
// RepartitionDown
//

template<typename T,typename Int>
void RepartitionDown
( M& AT, M& A0,
         M& A1,
  M& AB, M& A2, Int bsize=Blocksize() );

template<typename T,Distribution U,Distribution V,typename Int>
void RepartitionDown
( DM& AT, DM& A0,
          DM& A1,
  DM& AB, DM& A2, Int bsize=Blocksize() );

template<typename T,typename Int>
void LockedRepartitionDown
( const M& AT, M& A0,
               M& A1,
  const M& AB, M& A2, Int bsize=Blocksize() );

template<typename T,Distribution U,Distribution V,typename Int>
void LockedRepartitionDown
( const DM& AT, DM& A0,
                DM& A1,
  const DM& AB, DM& A2, Int bsize=Blocksize() );

//
// RepartitionLeft
//

template<typename T,typename Int>
void RepartitionLeft
( M& AL, M& AR,
  M& A0, M& A1, M& A2, Int bsize=Blocksize() );

template<typename T,Distribution U,Distribution V,typename Int>
void RepartitionLeft
( DM& AL, DM& AR,
  DM& A0, DM& A1, DM& A2, Int bsize=Blocksize() );

template<typename T,typename Int>
void LockedRepartitionLeft
( const M& AL, const M& AR,
  M& A0, M& A1, M& A2, Int bsize=Blocksize() );

template<typename T,Distribution U,Distribution V,typename Int>
void LockedRepartitionLeft
( const DM& AL, const DM& AR,
  DM& A0, DM& A1, DM& A2, Int bsize=Blocksize() );

//
// RepartitionRight
//

template<typename T,typename Int>
void RepartitionRight
( M& AL, M& AR,
  M& A0, M& A1, M& A2, Int bsize=Blocksize() );

template<typename T,Distribution U,Distribution V,typename Int>
void RepartitionRight
( DM& AL, DM& AR,
  DM& A0, DM& A1, DM& A2, Int bsize=Blocksize() );

template<typename T,typename Int>
void LockedRepartitionRight
( const M& AL, const M& AR,
  M& A0, M& A1, M& A2, Int bsize=Blocksize() );

template<typename T,Distribution U,Distribution V,typename Int>
void LockedRepartitionRight
( const DM& AL, const DM& AR,
  DM& A0, DM& A1, DM& A2, Int bsize=Blocksize() );

//
// RepartitionUpDiagonal
//

template<typename T,typename Int>
void RepartitionUpDiagonal
( M& ATL, M& ATR, M& A00, M& A01, M& A02,
                  M& A10, M& A11, M& A12,
  M& ABL, M& ABR, M& A20, M& A21, M& A22, Int bsize=Blocksize() );

template<typename T,Distribution U,Distribution V,typename Int>
void RepartitionUpDiagonal
( DM& ATL, DM& ATR, DM& A00, DM& A01, DM& A02,
                    DM& A10, DM& A11, DM& A12,
  DM& ABL, DM& ABR, DM& A20, DM& A21, DM& A22, Int bsize=Blocksize() );

template<typename T,typename Int>
void LockedRepartitionUpDiagonal
( const M& ATL, const M& ATR, M& A00, M& A01, M& A02,
                              M& A10, M& A11, M& A12,
  const M& ABL, const M& ABR, M& A20, M& A21, M& A22, Int bsize=Blocksize() );

template<typename T,Distribution U,Distribution V,typename Int>
void LockedRepartitionUpDiagonal
( const DM& ATL, const DM& ATR, DM& A00, DM& A01, DM& A02,
                                DM& A10, DM& A11, DM& A12,
  const DM& ABL, const DM& ABR, DM& A20, DM& A21, DM& A22, 
  Int bsize=Blocksize() );

//
// RepartitionDownDiagonal
//

template<typename T,typename Int>
void RepartitionDownDiagonal
( M& ATL, M& ATR, M& A00, M& A01, M& A02,
                  M& A10, M& A11, M& A12,
  M& ABL, M& ABR, M& A20, M& A21, M& A22, Int bsize=Blocksize() );

template<typename T,Distribution U,Distribution V,typename Int>
void RepartitionDownDiagonal
( DM& ATL, DM& ATR, DM& A00, DM& A01, DM& A02,
                    DM& A10, DM& A11, DM& A12,
  DM& ABL, DM& ABR, DM& A20, DM& A21, DM& A22, Int bsize=Blocksize() );

template<typename T,typename Int>
void LockedRepartitionDownDiagonal
( const M& ATL, const M& ATR, M& A00, M& A01, M& A02,
                              M& A10, M& A11, M& A12,
  const M& ABL, const M& ABR, M& A20, M& A21, M& A22, Int bsize=Blocksize() );

template<typename T,Distribution U,Distribution V,typename Int>
void LockedRepartitionDownDiagonal
( const DM& ATL, const DM& ATR, DM& A00, DM& A01, DM& A02,
                                DM& A10, DM& A11, DM& A12,
  const DM& ABL, const DM& ABR, DM& A20, DM& A21, DM& A22, 
  Int bsize=Blocksize() );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

//
// RepartitionUp
//

template<typename T,typename Int>
inline void
RepartitionUp
( M& AT, M& A0,
         M& A1,
  M& AB, M& A2, Int A1Height )
{
#ifndef RELEASE
    PushCallStack("RepartitionUp [Matrix]");
    if( (AT.Buffer() + AT.Height()) != AB.Buffer() )
        throw std::logic_error("Noncontiguous 2x1 array of matrices");
#endif
    A1Height = std::min(A1Height,AT.Height());
    const Int offset = AT.Height()-A1Height; 
    A0.View( AT, 0,      0, offset,   AT.Width() );
    A1.View( AT, offset, 0, A1Height, AT.Width() );
    A2.View( AB );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V,typename Int>
inline void
RepartitionUp
( DM& AT, DM& A0,
          DM& A1,
  DM& AB, DM& A2, Int A1Height )
{
#ifndef RELEASE
    PushCallStack("RepartitionUp [DistMatrix]");
    if( (AT.LocalMatrix().Buffer() + AT.LocalHeight()) != 
         AB.LocalMatrix().Buffer() )
        throw std::logic_error
        ("Noncontiguous 2x1 array of distributed matrices");
#endif
    A1Height = std::min(A1Height,AT.Height());
    const Int offset = AT.Height()-A1Height; 
    A0.View( AT, 0,      0, offset,   AT.Width() );
    A1.View( AT, offset, 0, A1Height, AT.Width() );
    A2.View( AB );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
LockedRepartitionUp
( const M& AT, M& A0,
               M& A1,
  const M& AB, M& A2, Int A1Height )
{
#ifndef RELEASE
    PushCallStack("LockedRepartitionUp [Matrix]");
    if( (AT.LockedBuffer() + AT.Height()) != AB.LockedBuffer() )
        throw std::logic_error("Noncontiguous 2x1 array of matrices");
#endif
    A1Height = std::min(A1Height,AT.Height());
    const Int offset = AT.Height()-A1Height;
    A0.LockedView( AT, 0,      0, offset,   AT.Width() );
    A1.LockedView( AT, offset, 0, A1Height, AT.Width() );
    A2.LockedView( AB );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V,typename Int>
inline void
LockedRepartitionUp
( const DM& AT, DM& A0,
                DM& A1,
  const DM& AB, DM& A2, Int A1Height )
{
#ifndef RELEASE
    PushCallStack("LockedRepartitionUp [DistMatrix]");
    if( (AT.LockedLocalMatrix().LockedBuffer() + AT.LocalHeight()) != 
         AB.LockedLocalMatrix().LockedBuffer() )
        throw std::logic_error
        ("Noncontiguous 2x1 array of distributed matrices");
#endif
    A1Height = std::min(A1Height,AT.Height());
    const Int offset = AT.Height()-A1Height;
    A0.LockedView( AT, 0,      0, offset,   AT.Width() );
    A1.LockedView( AT, offset, 0, A1Height, AT.Width() );
    A2.LockedView( AB );
#ifndef RELEASE
    PopCallStack();
#endif
}

//
// RepartitionDown
//

template<typename T,typename Int>
inline void
RepartitionDown
( M& AT, M& A0,
         M& A1,
  M& AB, M& A2, Int A1Height )
{
#ifndef RELEASE
    PushCallStack("RepartitionDown [Matrix]");
    if( (AT.Buffer() + AT.Height()) != AB.Buffer() )
        throw std::logic_error("Noncontiguous 2x1 array of matrices");
#endif
    A1Height = std::min(A1Height,AB.Height());
    const Int offset = AB.Height()-A1Height; 
    A0.View( AT );
    A1.View( AB, 0,        0, A1Height, AB.Width() );
    A2.View( AB, A1Height, 0, offset,   AB.Width() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V,typename Int>
inline void
RepartitionDown
( DM& AT, DM& A0,
          DM& A1,
  DM& AB, DM& A2, Int A1Height )
{
#ifndef RELEASE
    PushCallStack("RepartitionDown [DistMatrix]");
    if( (AT.LocalMatrix().Buffer() + AT.LocalHeight()) != 
         AB.LocalMatrix().Buffer() )
        throw std::logic_error
        ("Noncontiguous 2x1 array of distributed matrices");
#endif
    A1Height = std::min(A1Height,AB.Height());
    const Int offset = AB.Height()-A1Height; 
    A0.View( AT );
    A1.View( AB, 0,        0, A1Height, AB.Width() );
    A2.View( AB, A1Height, 0, offset, AB.Width() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
LockedRepartitionDown
( const M& AT, M& A0,
               M& A1,
  const M& AB, M& A2, Int A1Height )
{
#ifndef RELEASE
    PushCallStack("LockedRepartitionDown [Matrix]");
    if( (AT.LockedBuffer() + AT.Height()) != AB.LockedBuffer() )
        throw std::logic_error("Noncontiguous 2x1 array of matrices");
#endif
    A1Height = std::min(A1Height,AB.Height());
    const Int offset = AB.Height()-A1Height;
    A0.LockedView( AT );
    A1.LockedView( AB, 0,        0, A1Height, AB.Width() );
    A2.LockedView( AB, A1Height, 0, offset,   AB.Width() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V,typename Int>
inline void
LockedRepartitionDown
( const DM& AT, DM& A0,
                DM& A1,
  const DM& AB, DM& A2, Int A1Height )
{
#ifndef RELEASE
    PushCallStack("LockedRepartitionDown [DistMatrix]");
    if( (AT.LockedLocalMatrix().LockedBuffer() + AT.LocalHeight()) != 
         AB.LockedLocalMatrix().LockedBuffer() )
        throw std::logic_error
        ("Noncontiguous 2x1 array of distributed matrices");
#endif
    A1Height = std::min(A1Height,AB.Height());
    const Int offset = AB.Height()-A1Height;
    A0.LockedView( AT );
    A1.LockedView( AB, 0,        0, A1Height, AB.Width() );
    A2.LockedView( AB, A1Height, 0, offset,   AB.Width() );
#ifndef RELEASE
    PopCallStack();
#endif
}

//
// RepartitionLeft
//

template<typename T,typename Int>
inline void
RepartitionLeft
( M& AL, M& AR,
  M& A0, M& A1, M& A2, Int A1Width )
{
#ifndef RELEASE
    PushCallStack("RepartitionLeft [Matrix]");
    if( (AL.Buffer() + AL.Width()*AL.LDim()) != AR.Buffer() )
        throw std::logic_error("Noncontiguous 1x2 array of matrices");
#endif
    A1Width = std::min(A1Width,AL.Width());
    const Int offset = AL.Width()-A1Width;
    A0.View( AL, 0, 0,      AL.Height(), offset   );
    A1.View( AL, 0, offset, AL.Height(), A1Width  );
    A2.View( AR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V,typename Int>
inline void
RepartitionLeft
( DM& AL, DM& AR,
  DM& A0, DM& A1, DM& A2, Int A1Width )
{
#ifndef RELEASE
    PushCallStack("RepartitionLeft [DistMatrix]");
    if( (AL.LocalMatrix().Buffer() + AL.LocalWidth()*AL.LocalLDim())
         != AR.LocalMatrix().Buffer() )
        throw std::logic_error
        ("Noncontiguous 1x2 array of distributed matrices");
#endif
    A1Width = std::min(A1Width,AL.Width());
    const Int offset = AL.Width()-A1Width;
    A0.View( AL, 0, 0,      AL.Height(), offset  );
    A1.View( AL, 0, offset, AL.Height(), A1Width );
    A2.View( AR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
LockedRepartitionLeft
( const M& AL, const M& AR,
  M& A0, M& A1, M& A2, Int A1Width )
{
#ifndef RELEASE
    PushCallStack("LockedRepartitionLeft [Matrix]");
    if( (AL.LockedBuffer() + AL.Width()*AL.LDim()) != AR.LockedBuffer() )
        throw std::logic_error("Noncontiguous 1x2 array of matrices");
#endif
    A1Width = std::min(A1Width,AL.Width());
    const Int offset = AL.Width()-A1Width;
    A0.LockedView( AL, 0, 0,      AL.Height(), offset  );
    A1.LockedView( AL, 0, offset, AL.Height(), A1Width );
    A2.LockedView( AR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V,typename Int>
inline void
LockedRepartitionLeft
( const DM& AL, const DM& AR,
  DM& A0, DM& A1, DM& A2, Int A1Width )
{
#ifndef RELEASE
    PushCallStack("LockedRepartitionLeft [DistMatrix]");
    if( (AL.LockedLocalMatrix().LockedBuffer() + AL.LocalWidth()*AL.LocalLDim())
         != AR.LockedLocalMatrix().LockedBuffer() )
        throw std::logic_error
        ("Noncontiguous 1x1 array of distributed matrices");
#endif
    A1Width = std::min(A1Width,AL.Width());
    const Int offset = AL.Width()-A1Width;
    A0.LockedView( AL, 0, 0,      AL.Height(), offset  );
    A1.LockedView( AL, 0, offset, AL.Height(), A1Width );
    A2.LockedView( AR );
#ifndef RELEASE
    PopCallStack();
#endif
}

//
// RepartitionRight
//

template<typename T,typename Int>
inline void
RepartitionRight
( M& AL, M& AR,
  M& A0, M& A1, M& A2, Int A1Width )
{
#ifndef RELEASE
    PushCallStack("RepartitionRight [Matrix]");
    if( (AL.Buffer() + AL.Width()*AL.LDim()) != AR.Buffer() )
        throw std::logic_error("Noncontiguous 1x2 array of matrices");
#endif
    A1Width = std::min(A1Width,AR.Width());
    const Int offset = AR.Width()-A1Width;
    A0.View( AL );
    A1.View( AR, 0, 0,       AR.Height(), A1Width );
    A2.View( AR, 0, A1Width, AR.Height(), offset  );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V,typename Int>
inline void
RepartitionRight
( DM& AL, DM& AR,
  DM& A0, DM& A1, DM& A2, Int A1Width )
{
#ifndef RELEASE
    PushCallStack("RepartitionRight [DistMatrix]");
    if( (AL.LocalMatrix().Buffer() + AL.LocalWidth()*AL.LocalLDim()) 
         != AR.LocalMatrix().Buffer() )
        throw std::logic_error
        ("Noncontiguous 1x2 array of distributed matrices");
#endif
    A1Width = std::min(A1Width,AR.Width());
    const Int offset = AR.Width()-A1Width;
    A0.View( AL );
    A1.View( AR, 0, 0,       AR.Height(), A1Width );
    A2.View( AR, 0, A1Width, AR.Height(), offset  );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
LockedRepartitionRight
( const M& AL, const M& AR,
  M& A0, M& A1, M& A2, Int A1Width )
{
#ifndef RELEASE
    PushCallStack("LockedRepartitionRight [Matrix]");
    if( (AL.LockedBuffer() + AL.Width()*AL.LDim()) != AR.LockedBuffer() )
        throw std::logic_error("Noncontiguous 1x2 array of matrices");
#endif
    A1Width = std::min(A1Width,AR.Width());
    const Int offset = AR.Width()-A1Width;
    A0.LockedView( AL );
    A1.LockedView( AR, 0, 0,       AR.Height(), A1Width );
    A2.LockedView( AR, 0, A1Width, AR.Height(), offset  );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V,typename Int>
inline void
LockedRepartitionRight
( const DM& AL, const DM& AR,
  DM& A0, DM& A1, DM& A2, Int A1Width )
{
#ifndef RELEASE
    PushCallStack("LockedRepartitionRight [DistMatrix]");
    if( (AL.LockedLocalMatrix().LockedBuffer() + AL.LocalWidth()*AL.LocalLDim())
         != AR.LockedLocalMatrix().LockedBuffer() )
        throw std::logic_error
        ("Noncontiguous 1x2 DistMatrices in LockedRepartitionRight");
#endif
    A1Width = std::min(A1Width,AR.Width());
    const Int offset = AR.Width()-A1Width;
    A0.LockedView( AL );
    A1.LockedView( AR, 0, 0,       AR.Height(), A1Width );
    A2.LockedView( AR, 0, A1Width, AR.Height(), offset  );
#ifndef RELEASE
    PopCallStack();
#endif
}

//
// RepartitionUpDiagonal
//

template<typename T,typename Int>
inline void
RepartitionUpDiagonal
( M& ATL, M& ATR, M& A00, M& A01, M& A02,
                  M& A10, M& A11, M& A12,
  M& ABL, M& ABR, M& A20, M& A21, M& A22, Int bsize )
{
#ifndef RELEASE
    PushCallStack("RepartitionUpDiagonal [Matrix]");
    if( (ATL.Buffer() + ATL.Height()) != ABL.Buffer() ||
        (ATR.Buffer() + ATR.Height()) != ABR.Buffer() ||
        (ATL.Buffer() + ATL.Width()*ATL.LDim()) != ATR.Buffer() ||
        (ABL.Buffer() + ABL.Width()*ABL.LDim()) != ABR.Buffer()    )
        throw std::logic_error("Noncontiguous 2x2 grid of matrices");
#endif
    bsize = std::min(bsize,std::min(ATL.Height(),ATL.Width()));
    const Int vOffset = ATL.Height()-bsize;
    const Int hOffset = ATL.Width()-bsize;
    A00.View( ATL, 0,       0,       vOffset,      hOffset     );
    A01.View( ATL, 0,       hOffset, vOffset,      bsize       );
    A02.View( ATR, 0,       0,       vOffset,      ATR.Width() );
    A10.View( ATL, vOffset, 0,       bsize,        hOffset     );
    A11.View( ATL, vOffset, hOffset, bsize,        bsize       );
    A12.View( ATR, vOffset, 0,       bsize,        ATR.Width() );
    A20.View( ABL, 0,       0,       ABL.Height(), hOffset     );
    A21.View( ABL, 0,       hOffset, ABL.Height(), bsize       );
    A22.View( ABR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V,typename Int>
inline void
RepartitionUpDiagonal
( DM& ATL, DM& ATR, DM& A00, DM& A01, DM& A02,
                    DM& A10, DM& A11, DM& A12,
  DM& ABL, DM& ABR, DM& A20, DM& A21, DM& A22, Int bsize )
{
#ifndef RELEASE
    PushCallStack("RepartitionUpDiagonal [DistMatrix]");
    if( (ATL.LocalMatrix().Buffer() + ATL.LocalHeight()) != 
         ABL.LocalMatrix().Buffer() ||
        (ATR.LocalMatrix().Buffer() + ATR.LocalHeight()) != 
         ABR.LocalMatrix().Buffer() ||
        (ATL.LocalMatrix().Buffer() + ATL.LocalWidth()*ATL.LocalLDim())
         != ATR.LocalMatrix().Buffer() ||
        (ABL.LocalMatrix().Buffer() + ABL.LocalWidth()*ABL.LocalLDim())
         != ABR.LocalMatrix().Buffer() )
        throw std::logic_error
        ("Noncontiguous 2x2 grid of distributed matrices");
#endif
    bsize = std::min(bsize,std::min(ATL.Height(),ATL.Width()));
    const Int vOffset = ATL.Height()-bsize;
    const Int hOffset = ATL.Width()-bsize;
    A00.View( ATL, 0,       0,       vOffset,      hOffset     );
    A01.View( ATL, 0,       hOffset, vOffset,      bsize       );
    A02.View( ATR, 0,       0,       vOffset,      ATR.Width() );
    A10.View( ATL, vOffset, 0,       bsize,        hOffset     );
    A11.View( ATL, vOffset, hOffset, bsize,        bsize       );
    A12.View( ATR, vOffset, 0,       bsize,        ATR.Width() );
    A20.View( ABL, 0,       0,       ABL.Height(), hOffset     );
    A21.View( ABL, 0,       hOffset, ABL.Height(), bsize       );
    A22.View( ABR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
LockedRepartitionUpDiagonal
( const M& ATL, const M& ATR, M& A00, M& A01, M& A02,
                              M& A10, M& A11, M& A12,
  const M& ABL, const M& ABR, M& A20, M& A21, M& A22, Int bsize )
{
#ifndef RELEASE
    PushCallStack("LockedRepartitionUpDiagonal [Matrix]");
    if( (ATL.LockedBuffer() + ATL.Height()) != ABL.LockedBuffer() ||
        (ATR.LockedBuffer() + ATR.Height()) != ABR.LockedBuffer() ||
        (ATL.LockedBuffer() + ATL.Width()*ATL.LDim()) != ATR.LockedBuffer() ||
        (ABL.LockedBuffer() + ABL.Width()*ABL.LDim()) != ABR.LockedBuffer() )
        throw std::logic_error("Noncontiguous 2x2 grid of matrices");
#endif
    bsize = std::min(bsize,std::min(ATL.Height(),ATL.Width()));
    const Int vOffset = ATL.Height()-bsize;
    const Int hOffset = ATL.Width()-bsize;
    A00.LockedView( ATL, 0,       0,       vOffset,      hOffset     );
    A01.LockedView( ATL, 0,       hOffset, vOffset,      bsize       );
    A02.LockedView( ATR, 0,       0,       vOffset,      ATR.Width() );
    A10.LockedView( ATL, vOffset, 0,       bsize,        hOffset     );
    A11.LockedView( ATL, vOffset, hOffset, bsize,        bsize       );
    A12.LockedView( ATR, vOffset, 0,       bsize,        ATR.Width() );
    A20.LockedView( ABL, 0,       0,       ABL.Height(), hOffset     );
    A21.LockedView( ABL, 0,       hOffset, ABL.Height(), bsize       );
    A22.LockedView( ABR );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V,typename Int>
inline void
LockedRepartitionUpDiagonal
( const DM& ATL, const DM& ATR, DM& A00, DM& A01, DM& A02,
                                DM& A10, DM& A11, DM& A12,
  const DM& ABL, const DM& ABR, DM& A20, DM& A21, DM& A22, Int bsize )
{
#ifndef RELEASE
    PushCallStack("LockedRepartitionUpDiagonal [DistMatrix]");
    if( (ATL.LockedLocalMatrix().LockedBuffer()+ATL.LocalHeight()) != 
         ABL.LockedLocalMatrix().LockedBuffer() ||
        (ATR.LockedLocalMatrix().LockedBuffer()+ATR.LocalHeight()) != 
         ABR.LockedLocalMatrix().LockedBuffer() ||
        (ATL.LockedLocalMatrix().LockedBuffer()+
         ATL.LocalWidth()*ATL.LocalLDim()) != 
         ATR.LockedLocalMatrix().LockedBuffer() ||
        (ABL.LockedLocalMatrix().LockedBuffer()+
         ABL.LocalWidth()*ABL.LocalLDim()) !=
         ABR.LockedLocalMatrix().LockedBuffer() )
        throw std::logic_error
        ("Noncontiguous 2x2 grid of distributed matrices");
#endif
    bsize = std::min(bsize,std::min(ATL.Height(),ATL.Width()));
    const Int vOffset = ATL.Height()-bsize;
    const Int hOffset = ATL.Width()-bsize;
    A00.LockedView( ATL, 0,       0,       vOffset,      hOffset     );
    A01.LockedView( ATL, 0,       hOffset, vOffset,      bsize       );
    A02.LockedView( ATR, 0,       0,       vOffset,      ATR.Width() );
    A10.LockedView( ATL, vOffset, 0,       bsize,        hOffset     );
    A11.LockedView( ATL, vOffset, hOffset, bsize,        bsize       );
    A12.LockedView( ATR, vOffset, 0,       bsize,        ATR.Width() );
    A20.LockedView( ABL, 0,       0,       ABL.Height(), hOffset     );
    A21.LockedView( ABL, 0,       hOffset, ABL.Height(), bsize       );
    A22.LockedView( ABR );
#ifndef RELEASE
    PopCallStack();
#endif
}

//
// RepartitionDownDiagonal
//

template<typename T,typename Int>
inline void
RepartitionDownDiagonal
( M& ATL, M& ATR, M& A00, M& A01, M& A02,
                  M& A10, M& A11, M& A12,
  M& ABL, M& ABR, M& A20, M& A21, M& A22, Int bsize )
{
#ifndef RELEASE
    PushCallStack("RepartitionDownDiagonal [Matrix]");
    if( (ATL.Buffer() + ATL.Height()) != ABL.Buffer() ||
        (ATR.Buffer() + ATR.Height()) != ABR.Buffer() ||
        (ATL.Buffer() + ATL.Width()*ATL.LDim()) != ATR.Buffer() ||
        (ABL.Buffer() + ABL.Width()*ABL.LDim()) != ABR.Buffer()    )
        throw std::logic_error("Noncontiguous 2x2 grid of matrices");
#endif
    bsize = std::min(bsize,std::min(ABR.Height(),ABR.Width()));
    const Int vOffset = ABR.Height()-bsize;
    const Int hOffset = ABR.Width()-bsize;
    A00.View( ATL );
    A01.View( ATR, 0,     0,     ATL.Height(), bsize       );
    A02.View( ATR, 0,     bsize, ATL.Height(), hOffset     );
    A10.View( ABL, 0,     0,     bsize,        ABL.Width() );
    A11.View( ABR, 0,     0,     bsize,        bsize       );
    A12.View( ABR, 0,     bsize, bsize,        hOffset     );
    A20.View( ABL, bsize, 0,     vOffset,      ABL.Width() );
    A21.View( ABR, bsize, 0,     vOffset,      bsize       );
    A22.View( ABR, bsize, bsize, vOffset,      hOffset     );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V,typename Int>
inline void
RepartitionDownDiagonal
( DM& ATL, DM& ATR, DM& A00, DM& A01, DM& A02,
                    DM& A10, DM& A11, DM& A12,
  DM& ABL, DM& ABR, DM& A20, DM& A21, DM& A22, Int bsize )
{
#ifndef RELEASE
    PushCallStack("RepartitionDownDiagonal [DistMatrix]");
    if( (ATL.LocalMatrix().Buffer() + ATL.LocalHeight()) != 
         ABL.LocalMatrix().Buffer() ||
        (ATR.LocalMatrix().Buffer() + ATR.LocalHeight()) != 
         ABR.LocalMatrix().Buffer() ||
        (ATL.LocalMatrix().Buffer() + ATL.LocalWidth()*ATL.LocalLDim()) !=
         ATR.LocalMatrix().Buffer() ||
        (ABL.LocalMatrix().Buffer() + ABL.LocalWidth()*ABL.LocalLDim()) != 
         ABR.LocalMatrix().Buffer() )
        throw std::logic_error
        ("Noncontiguous 2x2 grid of distributed matrices");
#endif
    bsize = std::min(bsize,std::min(ABR.Height(),ABR.Width()));
    const Int vOffset = ABR.Height()-bsize;
    const Int hOffset = ABR.Width()-bsize;
    A00.View( ATL );
    A01.View( ATR, 0,     0,     ATL.Height(), bsize       );
    A02.View( ATR, 0,     bsize, ATL.Height(), hOffset     );
    A10.View( ABL, 0,     0,     bsize,        ABL.Width() );
    A11.View( ABR, 0,     0,     bsize,        bsize       );
    A12.View( ABR, 0,     bsize, bsize,        hOffset     );
    A20.View( ABL, bsize, 0,     vOffset,      ABL.Width() );
    A21.View( ABR, bsize, 0,     vOffset,      bsize       );
    A22.View( ABR, bsize, bsize, vOffset,      hOffset     );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void
LockedRepartitionDownDiagonal
( const M& ATL, const M& ATR, M& A00, M& A01, M& A02,
                              M& A10, M& A11, M& A12,
  const M& ABL, const M& ABR, M& A20, M& A21, M& A22, Int bsize )
{
#ifndef RELEASE
    PushCallStack("LockedRepartitionDownDiagonal [Matrix]");
    if( (ATL.LockedBuffer() + ATL.Height()) != ABL.LockedBuffer() ||
        (ATR.LockedBuffer() + ATR.Height()) != ABR.LockedBuffer() ||
        (ATL.LockedBuffer() + ATL.Width()*ATL.LDim()) != ATR.LockedBuffer() ||
        (ABL.LockedBuffer() + ABL.Width()*ABL.LDim()) != ABR.LockedBuffer() )
        throw std::logic_error("Noncontiguous 2x2 grid of matrices");
#endif
    bsize = std::min(bsize,std::min(ABR.Height(),ABR.Width()));
    const Int vOffset = ABR.Height()-bsize;
    const Int hOffset = ABR.Width()-bsize;
    A00.LockedView( ATL );
    A01.LockedView( ATR, 0,     0,     ATL.Height(), bsize       );
    A02.LockedView( ATR, 0,     bsize, ATL.Height(), hOffset     ); 
    A10.LockedView( ABL, 0,     0,     bsize,        ABL.Width() );
    A11.LockedView( ABR, 0,     0,     bsize,        bsize       );
    A12.LockedView( ABR, 0,     bsize, bsize,        hOffset     );
    A20.LockedView( ABL, bsize, 0,     vOffset,      ABL.Width() );
    A21.LockedView( ABR, bsize, 0,     vOffset,      bsize       );
    A22.LockedView( ABR, bsize, bsize, vOffset,      hOffset     );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V,typename Int>
inline void
LockedRepartitionDownDiagonal
( const DM& ATL, const DM& ATR, DM& A00, DM& A01, DM& A02,
                                DM& A10, DM& A11, DM& A12,
  const DM& ABL, const DM& ABR, DM& A20, DM& A21, DM& A22, Int bsize )
{
#ifndef RELEASE
    PushCallStack("LockedRepartitionDownDiagonal [DistMatrix]");
    if( (ATL.LockedLocalMatrix().LockedBuffer()+ATL.LocalHeight()) != 
         ABL.LockedLocalMatrix().LockedBuffer() ||
        (ATR.LockedLocalMatrix().LockedBuffer()+ATR.LocalHeight()) != 
         ABR.LockedLocalMatrix().LockedBuffer() ||
        (ATL.LockedLocalMatrix().LockedBuffer()+
         ATL.LocalWidth()*ATL.LocalLDim()) !=
         ATR.LockedLocalMatrix().LockedBuffer() ||
        (ABL.LockedLocalMatrix().LockedBuffer()+
         ABL.LocalWidth()*ABL.LocalLDim()) !=
         ABR.LockedLocalMatrix().LockedBuffer() )
        throw std::logic_error
        ("Noncontiguous 2x2 grid of distributed matrices");
#endif
    bsize = std::min(bsize,std::min(ABR.Height(),ABR.Width()));
    const Int vOffset = ABR.Height()-bsize;
    const Int hOffset = ABR.Width()-bsize;
    A00.LockedView( ATL );
    A01.LockedView( ATR, 0,     0,     ATR.Height(), bsize  );
    A02.LockedView( ATR, 0,     bsize, ATR.Height(), hOffset     );
    A10.LockedView( ABL, 0,     0,     bsize,        ABL.Width() );
    A11.LockedView( ABR, 0,     0,     bsize,        bsize       );
    A12.LockedView( ABR, 0,     bsize, bsize,        hOffset     );
    A20.LockedView( ABL, bsize, 0,     vOffset,      ABL.Width() );
    A21.LockedView( ABR, bsize, 0,     vOffset,      bsize       );
    A22.LockedView( ABR, bsize, bsize, vOffset,      hOffset     );
#ifndef RELEASE
    PopCallStack();
#endif
}

#undef DM
#undef M

} // namespace elem

#endif /* ELEMENTAL_REPARTITION_HPP */

