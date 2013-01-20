/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef CORE_VIEW_IMPL_HPP
#define CORE_VIEW_IMPL_HPP

namespace elem {

template<typename T,Distribution U,Distribution V,typename Int>
inline void HandleDiagPath
( DistMatrix<T,U,V,Int>& A, const DistMatrix<T,U,V,Int>& B )
{ }

template<typename T,typename Int>
inline void HandleDiagPath
( DistMatrix<T,MD,STAR,Int>& A, const DistMatrix<T,MD,STAR,Int>& B )
{ A.diagPath_ = B.diagPath_; } 

template<typename T,typename Int>
inline void HandleDiagPath
( DistMatrix<T,STAR,MD,Int>& A, const DistMatrix<T,STAR,MD,Int>& B )
{ A.diagPath_ = B.diagPath_; } 

template<typename T,typename Int>
inline void View
( Matrix<T,Int>& A, Matrix<T,Int>& B )
{
#ifndef RELEASE
    PushCallStack("View");
#endif
    A.Empty();
    A.height_ = B.Height();
    A.width_ = B.Width();
    A.ldim_ = B.LDim();
    A.data_ = B.Buffer();
    A.viewing_ = true;
    A.lockedView_ = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V,typename Int>
inline void View
( DistMatrix<T,U,V,Int>& A, DistMatrix<T,U,V,Int>& B )
{
#ifndef RELEASE
    PushCallStack("View");
#endif
    A.Empty();
    A.grid_ = B.grid_;
    A.height_ = B.Height();
    A.width_ = B.Width();
    A.colAlignment_ = B.ColAlignment();
    A.rowAlignment_ = B.RowAlignment();
    HandleDiagPath( A, B );
    A.viewing_ = true;
    if( A.Participating() )
    {
        A.colShift_ = B.ColShift();
        A.rowShift_ = B.RowShift();
        View( A.LocalMatrix(), B.LocalMatrix() );
    }
    else
    {
        A.colShift_ = 0;
        A.rowShift_ = 0;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void LockedView
( Matrix<T,Int>& A, const Matrix<T,Int>& B )
{
#ifndef RELEASE
    PushCallStack("LockedView");
#endif
    A.Empty();
    A.height_ = B.Height();
    A.width_ = B.Width();
    A.ldim_ = B.LDim();
    A.lockedData_ = B.LockedBuffer();
    A.viewing_ = true;
    A.lockedView_ = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V,typename Int>
inline void LockedView
( DistMatrix<T,U,V,Int>& A, const DistMatrix<T,U,V,Int>& B )
{
#ifndef RELEASE
    PushCallStack("LockedView");
#endif
    A.Empty();
    A.grid_ = B.grid_;
    A.height_ = B.Height();
    A.width_ = B.Width();
    A.colAlignment_ = B.ColAlignment();
    A.rowAlignment_ = B.RowAlignment();
    HandleDiagPath( A, B );
    A.viewing_ = true;
    A.lockedView_ = true;
    if( A.Participating() )
    {
        A.colShift_ = B.ColShift();
        A.rowShift_ = B.RowShift();
        LockedView( A.LocalMatrix(), B.LockedLocalMatrix() );
    }
    else 
    {
        A.colShift_ = 0;
        A.rowShift_ = 0;
    }
#ifndef RELEASE 
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void View
( Matrix<T,Int>& A, Matrix<T,Int>& B,
  Int i, Int j, Int height, Int width )
{
#ifndef RELEASE
    PushCallStack("View");
    if( i < 0 || j < 0 )
        throw std::logic_error("Indices must be non-negative");
    if( height < 0 || width < 0 )
        throw std::logic_error("Height and width must be non-negative");
    if( (i+height) > B.Height() || (j+width) > B.Width() )
    {
        std::ostringstream msg;
        msg << "Trying to view outside of a Matrix: "
            << "(" << i << "," << j << ") up to (" 
            << i+height-1 << "," << j+width-1 << ") "
            << "of " << B.Height() << " x " << B.Width() << " Matrix.";
        throw std::logic_error( msg.str().c_str() );
    }
#endif
    A.Empty();
    A.height_ = height;
    A.width_ = width;
    A.ldim_ = B.LDim();
    A.data_ = B.Buffer(i,j);
    A.viewing_ = true;
    A.lockedView_ = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V,typename Int>
inline void View
( DistMatrix<T,U,V,Int>& A, DistMatrix<T,U,V,Int>& B,
  Int i, Int j, Int height, Int width )
{
#ifndef RELEASE
    PushCallStack("View");
    B.AssertValidSubmatrix( B, i, j, height, width );
#endif
    A.Empty();

    const elem::Grid& g = B.Grid();
    const Int colStride = B.ColStride();
    const Int rowStride = B.RowStride();

    A.grid_ = &g;
    A.height_ = height;
    A.width_  = width;

    A.colAlignment_ = (B.ColAlignment()+i) % colStride;
    A.rowAlignment_ = (B.RowAlignment()+j) % rowStride;
    HandleDiagPath( A, B );
    A.viewing_ = true;

    if( A.Participating() )
    {
        const Int colRank = B.ColRank();
        const Int rowRank = B.RowRank();
        A.colShift_ = Shift( colRank, A.ColAlignment(), colStride );
        A.rowShift_ = Shift( rowRank, A.RowAlignment(), rowStride );

        const Int localHeightBehind = LocalLength(i,B.ColShift(),colStride);
        const Int localWidthBehind  = LocalLength(j,B.RowShift(),rowStride);

        const Int localHeight = LocalLength( height, A.ColShift(), colStride );
        const Int localWidth  = LocalLength( width,  A.RowShift(), rowStride );

        View
        ( A.LocalMatrix(), B.LocalMatrix(), 
          localHeightBehind, localWidthBehind, localHeight, localWidth );
    }
    else
    {
        A.colShift_ = 0;
        A.rowShift_ = 0;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void LockedView
( Matrix<T,Int>& A, const Matrix<T,Int>& B,
  Int i, Int j, Int height, Int width )
{
#ifndef RELEASE
    PushCallStack("LockedView");
    if( i < 0 || j < 0 )
        throw std::logic_error("Indices must be non-negative");
    if( height < 0 || width < 0 )
        throw std::logic_error("Height and width must be non-negative");
    if( (i+height) > B.Height() || (j+width) > B.Width() )
    {
        std::ostringstream msg;
        msg << "Trying to view outside of a Matrix: "
            << "(" << i << "," << j << ") up to (" 
            << i+height-1 << "," << j+width-1 << ") "
            << "of " << B.Height() << " x " << B.Width() << " Matrix.";
        throw std::logic_error( msg.str().c_str() );
    }
#endif
    A.Empty();
    A.height_ = height;
    A.width_ = width;
    A.ldim_ = B.LDim();
    A.lockedData_ = B.LockedBuffer(i,j);
    A.viewing_ = true;
    A.lockedView_ = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V,typename Int>
inline void LockedView
( DistMatrix<T,U,V,Int>& A, const DistMatrix<T,U,V,Int>& B,
  Int i, Int j, Int height, Int width )
{
#ifndef RELEASE
    PushCallStack("LockedView");
    B.AssertValidSubmatrix( B, i, j, height, width );
#endif
    A.Empty();

    const elem::Grid& g = B.Grid();
    const Int colStride = B.ColStride();
    const Int rowStride = B.RowStride();
    const Int colRank = B.ColRank();
    const Int rowRank = B.RowRank();

    A.grid_ = &g;
    A.height_ = height;
    A.width_  = width;

    A.colAlignment_ = (B.ColAlignment()+i) % colStride;
    A.rowAlignment_ = (B.RowAlignment()+j) % rowStride;
    HandleDiagPath( A, B );
    A.viewing_ = true;
    A.lockedView_ = true;

    if( A.Participating() )
    {
        A.colShift_ = Shift( colRank, A.ColAlignment(), colStride );
        A.rowShift_ = Shift( rowRank, A.RowAlignment(), rowStride );

        const Int localHeightBehind = LocalLength(i,B.ColShift(),colStride);
        const Int localWidthBehind  = LocalLength(j,B.RowShift(),rowStride);

        const Int localHeight = LocalLength( height, A.ColShift(), colStride );
        const Int localWidth  = LocalLength( width,  A.RowShift(), rowStride );

        LockedView
        ( A.LocalMatrix(), B.LockedLocalMatrix(), 
          localHeightBehind, localWidthBehind, localHeight, localWidth );
    }
    else
    {
        A.colShift_ = 0;
        A.rowShift_ = 0;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void View1x2
( Matrix<T,Int>& A,
  Matrix<T,Int>& BL, Matrix<T,Int>& BR )
{
#ifndef RELEASE
    PushCallStack("View1x2");
    if( BL.Height() != BR.Height() )
        throw std::logic_error("1x2 must have consistent height to combine");
    if( BL.LDim() != BR.LDim() )
        throw std::logic_error("1x2 must have consistent ldims to combine");
    if( BR.Buffer() != (BL.Buffer()+BL.LDim()*BL.Width()) )
        throw std::logic_error("1x2 must have contiguous memory");
#endif
    A.Empty();
    A.height_ = BL.Height();
    A.width_ = BL.Width() + BR.Width();
    A.ldim_ = BL.LDim();
    A.data_ = BL.Buffer();
    A.viewing_ = true;
    A.lockedView_ = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V,typename Int>
inline void View1x2
( DistMatrix<T,U,V,Int>& A,
  DistMatrix<T,U,V,Int>& BL, DistMatrix<T,U,V,Int>& BR )
{
#ifndef RELEASE
    PushCallStack("View1x2");
    A.AssertConforming1x2( BL, BR );
    BL.AssertSameGrid( BR );
#endif
    A.Empty();
    A.grid_ = BL.grid_;
    A.height_ = BL.Height();
    A.width_ = BL.Width() + BR.Width();
    A.colAlignment_ = BL.ColAlignment();
    A.rowAlignment_ = BL.RowAlignment();
    HandleDiagPath( A, BL );
    A.viewing_ = true;
    if( A.Participating() )
    {
        A.colShift_ = BL.ColShift();
        A.rowShift_ = BL.RowShift();
        View1x2( A.LocalMatrix(), BL.LocalMatrix(), BR.LocalMatrix() );
    }
    else
    {
        A.colShift_ = 0;
        A.rowShift_ = 0;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void LockedView1x2
(       Matrix<T,Int>& A,
  const Matrix<T,Int>& BL,
  const Matrix<T,Int>& BR )
{
#ifndef RELEASE
    PushCallStack("LockedView1x2");
    if( BL.Height() != BR.Height() )
        throw std::logic_error("1x2 must have consistent height to combine");
    if( BL.LDim() != BR.LDim() )
        throw std::logic_error("1x2 must have consistent ldims to combine");
    if( BR.LockedBuffer() != (BL.LockedBuffer()+BL.LDim()*BL.Width()) )
        throw std::logic_error("1x2 must have contiguous memory");
#endif
    A.Empty();
    A.height_ = BL.Height();
    A.width_ = BL.Width() + BR.Width();
    A.ldim_ = BL.LDim();
    A.lockedData_ = BL.LockedBuffer();
    A.viewing_ = true;
    A.lockedView_ = true;
#ifndef RELEASE
    PopCallStack(); 
#endif
}

template<typename T,Distribution U,Distribution V,typename Int>
inline void LockedView1x2
(       DistMatrix<T,U,V,Int>& A,
  const DistMatrix<T,U,V,Int>& BL,
  const DistMatrix<T,U,V,Int>& BR )
{
#ifndef RELEASE
    PushCallStack("LockedView1x2");
    A.AssertConforming1x2( BL, BR );
    BL.AssertSameGrid( BR );
#endif
    A.Empty();
    A.grid_ = BL.grid_;
    A.height_ = BL.Height();
    A.width_ = BL.Width() + BR.Width();
    A.colAlignment_ = BL.ColAlignment();
    A.rowAlignment_ = BL.RowAlignment();
    HandleDiagPath( A, BL );
    A.viewing_ = true;
    A.lockedView_ = true;
    if( A.Participating() )
    {
        A.colShift_ = BL.ColShift();
        A.rowShift_ = BL.RowShift();
        LockedView1x2
        ( A.LocalMatrix(), BL.LockedLocalMatrix(), BR.LockedLocalMatrix() );
    }
    else
    {
        A.colShift_ = 0;
        A.rowShift_ = 0;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void View2x1
( Matrix<T,Int>& A,
  Matrix<T,Int>& BT,
  Matrix<T,Int>& BB )
{
#ifndef RELEASE
    PushCallStack("View2x1");
    if( BT.Width() != BB.Width() )
        throw std::logic_error("2x1 must have consistent width to combine");
    if( BT.LDim() != BB.LDim() )
        throw std::logic_error("2x1 must have consistent ldim to combine");
    if( BB.Buffer() != (BT.Buffer() + BT.Height()) )
        throw std::logic_error("2x1 must have contiguous memory");
#endif
    A.Empty();
    A.height_ = BT.Height() + BB.Height();
    A.width_ = BT.Width();
    A.ldim_ = BT.LDim();
    A.data_ = BT.Buffer();
    A.viewing_ = true;
    A.lockedView_ = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V,typename Int>
inline void View2x1
( DistMatrix<T,U,V,Int>& A,
  DistMatrix<T,U,V,Int>& BT,
  DistMatrix<T,U,V,Int>& BB )
{
#ifndef RELEASE
    PushCallStack("View2x1");
    A.AssertConforming2x1( BT, BB );
    BT.AssertSameGrid( BB );
#endif
    A.Empty();
    A.grid_ = BT.grid_;
    A.height_ = BT.Height() + BB.Height();
    A.width_ = BT.Width();
    A.colAlignment_ = BT.ColAlignment();
    A.rowAlignment_ = BT.RowAlignment();
    HandleDiagPath( A, BT );
    A.viewing_ = true;
    if( A.Participating() )
    {
        A.colShift_ = BT.ColShift();
        A.rowShift_ = BT.RowShift();
        View2x1( A.LocalMatrix(), BT.LocalMatrix(), BB.LocalMatrix() );
    }
    else
    {
        A.colShift_ = 0;
        A.rowShift_ = 0;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void LockedView2x1
(       Matrix<T,Int>& A,
  const Matrix<T,Int>& BT,
  const Matrix<T,Int>& BB )
{
#ifndef RELEASE
    PushCallStack("LockedView2x1");
    if( BT.Width() != BB.Width() )
        throw std::logic_error("2x1 must have consistent width to combine");
    if( BT.LDim() != BB.LDim() )
        throw std::logic_error("2x1 must have consistent ldim to combine");
    if( BB.LockedBuffer() != (BT.LockedBuffer() + BT.Height()) )
        throw std::logic_error("2x1 must have contiguous memory");
#endif
    A.Empty();
    A.height_ = BT.Height() + BB.Height();
    A.width_ = BT.Width();
    A.ldim_ = BT.LDim();
    A.lockedData_ = BT.LockedBuffer();
    A.viewing_ = true;
    A.lockedView_ = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V,typename Int>
inline void LockedView2x1
(       DistMatrix<T,U,V,Int>& A,
  const DistMatrix<T,U,V,Int>& BT,
  const DistMatrix<T,U,V,Int>& BB )
{
#ifndef RELEASE
    PushCallStack("LockedView2x1");
    A.AssertConforming2x1( BT, BB );
    BT.AssertSameGrid( BB );
#endif
    A.Empty();
    A.grid_ = BT.grid_;
    A.height_ = BT.Height() + BB.Height();
    A.width_ = BT.Width();
    A.colAlignment_ = BT.ColAlignment();
    A.rowAlignment_ = BT.RowAlignment();
    HandleDiagPath( A, BT );
    A.viewing_ = true;
    A.lockedView_ = true;
    if( A.Participating() )
    {
        A.colShift_ = BT.ColShift();
        A.rowShift_ = BT.RowShift();
        LockedView2x1
        ( A.LocalMatrix(), BT.LockedLocalMatrix(), BB.LockedLocalMatrix() );
    }
    else
    {
        A.colShift_ = 0;
        A.rowShift_ = 0;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void View2x2
( Matrix<T,Int>& A,
  Matrix<T,Int>& BTL, Matrix<T,Int>& BTR,
  Matrix<T,Int>& BBL, Matrix<T,Int>& BBR )
{
#ifndef RELEASE
    PushCallStack("View2x2");
    if( BTL.Width() != BBL.Width()   ||
        BTR.Width() != BBR.Width()   ||
        BTL.Height() != BTR.Height() ||
        BBL.Height() != BBR.Height()   )
        throw std::logic_error("2x2 must conform to combine");
    if( BTL.LDim() != BTR.LDim() ||
        BTR.LDim() != BBL.LDim() ||
        BBL.LDim() != BBR.LDim()   )
        throw std::logic_error("2x2 must have consistent ldims to combine");
    if( BBL.Buffer() != (BTL.Buffer() + BTL.Height()) ||
        BBR.Buffer() != (BTR.Buffer() + BTR.Height()) ||
        BTR.Buffer() != (BTL.Buffer() + BTL.LDim()*BTL.Width()) )
        throw std::logic_error("2x2 must have contiguous memory");
#endif
    A.Empty();
    A.height_ = BTL.Height() + BBL.Height();
    A.width_ = BTL.Width() + BTR.Width();
    A.ldim_ = BTL.LDim();
    A.data_ = BTL.Buffer();
    A.viewing_ = true;
    A.lockedView_ = false;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V,typename Int>
inline void View2x2
( DistMatrix<T,U,V,Int>& A,
  DistMatrix<T,U,V,Int>& BTL, DistMatrix<T,U,V,Int>& BTR,
  DistMatrix<T,U,V,Int>& BBL, DistMatrix<T,U,V,Int>& BBR )
{
#ifndef RELEASE
    PushCallStack("View2x2");
    A.AssertConforming2x2( BTL, BTR, BBL, BBR );
    BTL.AssertSameGrid( BTR );
    BTL.AssertSameGrid( BBL );
    BTL.AssertSameGrid( BBR );
#endif
    A.Empty();
    A.grid_ = BTL.grid_;
    A.height_ = BTL.Height() + BBL.Height();
    A.width_ = BTL.Width() + BTR.Width();
    A.colAlignment_ = BTL.ColAlignment();
    A.rowAlignment_ = BTL.RowAlignment();
    HandleDiagPath( A, BTL );
    A.viewing_ = true;
    if( A.Participating() )
    {
        A.colShift_ = BTL.ColShift();
        A.rowShift_ = BTL.RowShift();
        View2x2
        ( A.LocalMatrix(), BTL.LocalMatrix(), BTR.LocalMatrix(),
                           BBL.LocalMatrix(), BBR.LocalMatrix() ); 
    }
    else
    {
        A.colShift_ = 0;
        A.rowShift_ = 0;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,typename Int>
inline void LockedView2x2
(       Matrix<T,Int>& A,
  const Matrix<T,Int>& BTL,
  const Matrix<T,Int>& BTR,
  const Matrix<T,Int>& BBL,
  const Matrix<T,Int>& BBR )
{
#ifndef RELEASE
    PushCallStack("LockedView2x2");
    if( BTL.Width() != BBL.Width()   ||
        BTR.Width() != BBR.Width()   ||
        BTL.Height() != BTR.Height() ||
        BBL.Height() != BBR.Height()   )
        throw std::logic_error("2x2 must conform to combine");
    if( BTL.LDim() != BTR.LDim() ||
        BTR.LDim() != BBL.LDim() ||
        BBL.LDim() != BBR.LDim()   )
        throw std::logic_error("2x2 must have consistent ldims to combine");
    if( BBL.LockedBuffer() != (BTL.LockedBuffer() + BTL.Height()) ||
        BBR.LockedBuffer() != (BTR.LockedBuffer() + BTR.Height()) ||
        BTR.LockedBuffer() != (BTL.LockedBuffer() + BTL.LDim()*BTL.Width()) )
        throw std::logic_error("2x2 must have contiguous memory");
#endif
    A.Empty();
    A.height_ = BTL.Height() + BBL.Height();
    A.width_ = BTL.Width() + BTR.Width();
    A.ldim_ = BTL.LDim();
    A.lockedData_ = BTL.LockedBuffer();
    A.viewing_ = true;
    A.lockedView_ = true;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V,typename Int>
inline void LockedView2x2
(       DistMatrix<T,U,V,Int>& A,
  const DistMatrix<T,U,V,Int>& BTL,
  const DistMatrix<T,U,V,Int>& BTR,
  const DistMatrix<T,U,V,Int>& BBL,
  const DistMatrix<T,U,V,Int>& BBR )
{
#ifndef RELEASE
    PushCallStack("LockedView2x2");
    A.AssertConforming2x2( BTL, BTR, BBL, BBR );
    BTL.AssertSameGrid( BTR );
    BTL.AssertSameGrid( BBL );
    BTL.AssertSameGrid( BBR );
#endif
    A.Empty();
    A.grid_ = BTL.grid_;
    A.height_ = BTL.Height() + BBL.Height();
    A.width_ = BTL.Width() + BTR.Width();
    A.colAlignment_ = BTL.ColAlignment();
    A.rowAlignment_ = BTL.RowAlignment();
    HandleDiagPath( A, BTL );
    A.viewing_ = true;
    A.lockedView_ = true;
    if( A.Participating() )
    {
        A.colShift_ = BTL.ColShift();
        A.rowShift_ = BTL.RowShift();
        LockedView2x2
        ( A.LocalMatrix(), BTL.LockedLocalMatrix(), BTR.LockedLocalMatrix(),
                           BBL.LockedLocalMatrix(), BBR.LockedLocalMatrix() );
    }
    else
    {
        A.colShift_ = 0;
        A.rowShift_ = 0;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem

#endif // ifndef CORE_VIEW_IMPL_HPP
