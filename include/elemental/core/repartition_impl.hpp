/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef CORE_REPARTITION_IMPL_HPP
#define CORE_REPARTITION_IMPL_HPP

namespace elem {

// To make our life easier. Undef'd at the bottom of the header
#define M  Matrix<T,Int>
#define DM DistMatrix<T,U,V,Int>

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
    View( A0, AT, 0,      0, offset,   AT.Width() );
    View( A1, AT, offset, 0, A1Height, AT.Width() );
    View( A2, AB );
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
    View( A0, AT, 0,      0, offset,   AT.Width() );
    View( A1, AT, offset, 0, A1Height, AT.Width() );
    View( A2, AB );
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
    LockedView( A0, AT, 0,      0, offset,   AT.Width() );
    LockedView( A1, AT, offset, 0, A1Height, AT.Width() );
    LockedView( A2, AB );
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
    LockedView( A0, AT, 0,      0, offset,   AT.Width() );
    LockedView( A1, AT, offset, 0, A1Height, AT.Width() );
    LockedView( A2, AB );
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
    View( A0, AT );
    View( A1, AB, 0,        0, A1Height, AB.Width() );
    View( A2, AB, A1Height, 0, offset,   AB.Width() );
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
    View( A0, AT );
    View( A1, AB, 0,        0, A1Height, AB.Width() );
    View( A2, AB, A1Height, 0, offset, AB.Width() );
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
    LockedView( A0, AT );
    LockedView( A1, AB, 0,        0, A1Height, AB.Width() );
    LockedView( A2, AB, A1Height, 0, offset,   AB.Width() );
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
    LockedView( A0, AT );
    LockedView( A1, AB, 0,        0, A1Height, AB.Width() );
    LockedView( A2, AB, A1Height, 0, offset,   AB.Width() );
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
    View( A0, AL, 0, 0,      AL.Height(), offset   );
    View( A1, AL, 0, offset, AL.Height(), A1Width  );
    View( A2, AR );
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
    View( A0, AL, 0, 0,      AL.Height(), offset  );
    View( A1, AL, 0, offset, AL.Height(), A1Width );
    View( A2, AR );
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
    LockedView( A0, AL, 0, 0,      AL.Height(), offset  );
    LockedView( A1, AL, 0, offset, AL.Height(), A1Width );
    LockedView( A2, AR );
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
    LockedView( A0, AL, 0, 0,      AL.Height(), offset  );
    LockedView( A1, AL, 0, offset, AL.Height(), A1Width );
    LockedView( A2, AR );
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
    View( A0, AL );
    View( A1, AR, 0, 0,       AR.Height(), A1Width );
    View( A2, AR, 0, A1Width, AR.Height(), offset  );
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
    View( A0, AL );
    View( A1, AR, 0, 0,       AR.Height(), A1Width );
    View( A2, AR, 0, A1Width, AR.Height(), offset  );
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
    LockedView( A0, AL );
    LockedView( A1, AR, 0, 0,       AR.Height(), A1Width );
    LockedView( A2, AR, 0, A1Width, AR.Height(), offset  );
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
    LockedView( A0, AL );
    LockedView( A1, AR, 0, 0,       AR.Height(), A1Width );
    LockedView( A2, AR, 0, A1Width, AR.Height(), offset  );
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
    View( A00, ATL, 0,       0,       vOffset,      hOffset     );
    View( A01, ATL, 0,       hOffset, vOffset,      bsize       );
    View( A02, ATR, 0,       0,       vOffset,      ATR.Width() );
    View( A10, ATL, vOffset, 0,       bsize,        hOffset     );
    View( A11, ATL, vOffset, hOffset, bsize,        bsize       );
    View( A12, ATR, vOffset, 0,       bsize,        ATR.Width() );
    View( A20, ABL, 0,       0,       ABL.Height(), hOffset     );
    View( A21, ABL, 0,       hOffset, ABL.Height(), bsize       );
    View( A22, ABR );
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
    View( A00, ATL, 0,       0,       vOffset,      hOffset     );
    View( A01, ATL, 0,       hOffset, vOffset,      bsize       );
    View( A02, ATR, 0,       0,       vOffset,      ATR.Width() );
    View( A10, ATL, vOffset, 0,       bsize,        hOffset     );
    View( A11, ATL, vOffset, hOffset, bsize,        bsize       );
    View( A12, ATR, vOffset, 0,       bsize,        ATR.Width() );
    View( A20, ABL, 0,       0,       ABL.Height(), hOffset     );
    View( A21, ABL, 0,       hOffset, ABL.Height(), bsize       );
    View( A22, ABR );
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
    LockedView( A00, ATL, 0,       0,       vOffset,      hOffset     );
    LockedView( A01, ATL, 0,       hOffset, vOffset,      bsize       );
    LockedView( A02, ATR, 0,       0,       vOffset,      ATR.Width() );
    LockedView( A10, ATL, vOffset, 0,       bsize,        hOffset     );
    LockedView( A11, ATL, vOffset, hOffset, bsize,        bsize       );
    LockedView( A12, ATR, vOffset, 0,       bsize,        ATR.Width() );
    LockedView( A20, ABL, 0,       0,       ABL.Height(), hOffset     );
    LockedView( A21, ABL, 0,       hOffset, ABL.Height(), bsize       );
    LockedView( A22, ABR );
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
    LockedView( A00, ATL, 0,       0,       vOffset,      hOffset     );
    LockedView( A01, ATL, 0,       hOffset, vOffset,      bsize       );
    LockedView( A02, ATR, 0,       0,       vOffset,      ATR.Width() );
    LockedView( A10, ATL, vOffset, 0,       bsize,        hOffset     );
    LockedView( A11, ATL, vOffset, hOffset, bsize,        bsize       );
    LockedView( A12, ATR, vOffset, 0,       bsize,        ATR.Width() );
    LockedView( A20, ABL, 0,       0,       ABL.Height(), hOffset     );
    LockedView( A21, ABL, 0,       hOffset, ABL.Height(), bsize       );
    LockedView( A22, ABR );
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
    View( A00, ATL );
    View( A01, ATR, 0,     0,     ATL.Height(), bsize       );
    View( A02, ATR, 0,     bsize, ATL.Height(), hOffset     );
    View( A10, ABL, 0,     0,     bsize,        ABL.Width() );
    View( A11, ABR, 0,     0,     bsize,        bsize       );
    View( A12, ABR, 0,     bsize, bsize,        hOffset     );
    View( A20, ABL, bsize, 0,     vOffset,      ABL.Width() );
    View( A21, ABR, bsize, 0,     vOffset,      bsize       );
    View( A22, ABR, bsize, bsize, vOffset,      hOffset     );
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
    View( A00, ATL );
    View( A01, ATR, 0,     0,     ATL.Height(), bsize       );
    View( A02, ATR, 0,     bsize, ATL.Height(), hOffset     );
    View( A10, ABL, 0,     0,     bsize,        ABL.Width() );
    View( A11, ABR, 0,     0,     bsize,        bsize       );
    View( A12, ABR, 0,     bsize, bsize,        hOffset     );
    View( A20, ABL, bsize, 0,     vOffset,      ABL.Width() );
    View( A21, ABR, bsize, 0,     vOffset,      bsize       );
    View( A22, ABR, bsize, bsize, vOffset,      hOffset     );
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
    LockedView( A00, ATL );
    LockedView( A01, ATR, 0,     0,     ATL.Height(), bsize       );
    LockedView( A02, ATR, 0,     bsize, ATL.Height(), hOffset     ); 
    LockedView( A10, ABL, 0,     0,     bsize,        ABL.Width() );
    LockedView( A11, ABR, 0,     0,     bsize,        bsize       );
    LockedView( A12, ABR, 0,     bsize, bsize,        hOffset     );
    LockedView( A20, ABL, bsize, 0,     vOffset,      ABL.Width() );
    LockedView( A21, ABR, bsize, 0,     vOffset,      bsize       );
    LockedView( A22, ABR, bsize, bsize, vOffset,      hOffset     );
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
    LockedView( A00, ATL );
    LockedView( A01, ATR, 0,     0,     ATR.Height(), bsize  );
    LockedView( A02, ATR, 0,     bsize, ATR.Height(), hOffset     );
    LockedView( A10, ABL, 0,     0,     bsize,        ABL.Width() );
    LockedView( A11, ABR, 0,     0,     bsize,        bsize       );
    LockedView( A12, ABR, 0,     bsize, bsize,        hOffset     );
    LockedView( A20, ABL, bsize, 0,     vOffset,      ABL.Width() );
    LockedView( A21, ABR, bsize, 0,     vOffset,      bsize       );
    LockedView( A22, ABR, bsize, bsize, vOffset,      hOffset     );
#ifndef RELEASE
    PopCallStack();
#endif
}

#undef DM
#undef M

} // namespace elem

#endif // ifndef CORE_REPARTITION_IMPL_HPP
