/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_VIEWS_VIEW_HPP
#define EL_VIEWS_VIEW_HPP

namespace El {

template<typename T>
inline void View( Matrix<T>& A, Matrix<T>& B )
{
    DEBUG_ONLY(CallStackEntry cse("View"))
    A.Attach( B.Height(), B.Width(), B.Buffer(), B.LDim() );
}

template<typename T>
inline Matrix<T> View( Matrix<T>& B )
{
    Matrix<T> A;
    View( A, B );
    return A;
}

template<typename T>
inline void View( AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B )
{
    DEBUG_ONLY(
        CallStackEntry cse("View");
        AssertSameDist( A.DistData(), B.DistData() );
    )
    A.Attach
    ( B.Height(), B.Width(), B.Grid(), B.ColAlign(), B.RowAlign(), 
      B.Buffer(), B.LDim(), B.Root() );
}

template<typename T,Dist U,Dist V>
inline DistMatrix<T,U,V> View( DistMatrix<T,U,V>& B )
{
    DistMatrix<T,U,V> A(B.Grid());
    View( A, B );
    return A;
}

template<typename T>
inline void View( AbstractBlockDistMatrix<T>& A, AbstractDistMatrix<T>& B )
{
    DEBUG_ONLY(
        CallStackEntry cse("View");
        AssertSameDist( A.DistData(), B.DistData() );
    )
    A.Attach
    ( B.Height(), B.Width(), B.Grid(), 1, 1, B.ColAlign(), B.RowAlign(), 0, 0,
      B.Buffer(), B.LDim(), B.Root() );
}

template<typename T>
inline void View( AbstractDistMatrix<T>& A, AbstractBlockDistMatrix<T>& B )
{
    DEBUG_ONLY(
        CallStackEntry cse("View");
        AssertSameDist( A.DistData(), B.DistData() );
    )
    if( B.BlockHeight() != 1 || B.BlockWidth() != 1 )
        LogicError("Block size was ",B.BlockHeight()," x ",B.BlockWidth(),
                    "instead of 1x1");
    A.Attach
    ( B.Height(), B.Width(), B.Grid(), B.ColAlign(), B.RowAlign(), 
      B.Buffer(), B.LDim(), B.Root() );
}

template<typename T>
inline void LockedView( Matrix<T>& A, const Matrix<T>& B )
{
    DEBUG_ONLY(CallStackEntry cse("LockedView"))
    A.LockedAttach( B.Height(), B.Width(), B.LockedBuffer(), B.LDim() );
}

template<typename T>
inline Matrix<T> LockedView( const Matrix<T>& B )
{
    Matrix<T> A;
    LockedView( A, B );
    return A;
}

template<typename T>
inline void LockedView
( AbstractDistMatrix<T>& A, const AbstractDistMatrix<T>& B )
{
    DEBUG_ONLY(
        CallStackEntry cse("LockedView");
        AssertSameDist( A.DistData(), B.DistData() );
    )
    A.LockedAttach
    ( B.Height(), B.Width(), B.Grid(), B.ColAlign(), B.RowAlign(), 
      B.LockedBuffer(), B.LDim(), B.Root() );
}

template<typename T,Dist U,Dist V>
inline DistMatrix<T,U,V> LockedView( const DistMatrix<T,U,V>& B )
{
    DistMatrix<T,U,V> A(B.Grid());
    LockedView( A, B );
    return A;
}

template<typename T>
inline void LockedView
( AbstractBlockDistMatrix<T>& A, const AbstractDistMatrix<T>& B )
{
    DEBUG_ONLY(
        CallStackEntry cse("LockedView");
        AssertSameDist( A.DistData(), B.DistData() );
    )
    A.LockedAttach
    ( B.Height(), B.Width(), B.Grid(), 1, 1, B.ColAlign(), B.RowAlign(), 0, 0,
      B.LockedBuffer(), B.LDim(), B.Root() );
}

template<typename T>
inline void LockedView
( AbstractDistMatrix<T>& A, const AbstractBlockDistMatrix<T>& B )
{
    DEBUG_ONLY(
        CallStackEntry cse("LockedView");
        AssertSameDist( A.DistData(), B.DistData() );
    )
    if( B.BlockHeight() != 1 || B.BlockWidth() != 1 )
        LogicError("Block size was ",B.BlockHeight()," x ",B.BlockWidth(),
                    "instead of 1x1");
    A.LockedAttach
    ( B.Height(), B.Width(), B.Grid(), B.ColAlign(), B.RowAlign(), 
      B.LockedBuffer(), B.LDim(), B.Root() );
}

template<typename T>
inline void View
( Matrix<T>& A, Matrix<T>& B,
  Int i, Int j, Int height, Int width )
{
    DEBUG_ONLY(
        CallStackEntry cse("View");
        if( i < 0 || j < 0 )
            LogicError("Indices must be non-negative");
        if( height < 0 || width < 0 )
            LogicError("Height and width must be non-negative");
        if( (i+height) > B.Height() || (j+width) > B.Width() )
            LogicError
            ("Trying to view outside of a Matrix: (",i,",",j,") up to (", 
             i+height-1,",",j+width-1,") of ",B.Height()," x ",B.Width(),
             " Matrix");
    )
    A.Attach( height, width, B.Buffer(i,j), B.LDim() );
}

template<typename T>
inline Matrix<T> View( Matrix<T>& B, Int i, Int j, Int height, Int width )
{
    Matrix<T> A;
    View( A, B, i, j, height, width );
    return A;
}

template<typename T>
inline void View
( AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B,
  Int i, Int j, Int height, Int width )
{
    DEBUG_ONLY(
        CallStackEntry cse("View");
        AssertSameDist( A.DistData(), B.DistData() );
        B.AssertValidSubmatrix( i, j, height, width );
    )
    const Int colAlign = (B.ColAlign()+i) % B.ColStride();
    const Int rowAlign = (B.RowAlign()+j) % B.RowStride();
    if( B.Participating() )
    {
        const Int iLoc = Length( i, B.ColShift(), B.ColStride() );
        const Int jLoc = Length( j, B.RowShift(), B.RowStride() );
        A.Attach
        ( height, width, B.Grid(), colAlign, rowAlign, 
          B.Buffer(iLoc,jLoc), B.LDim(), B.Root() );
    }
    else
    {
        A.Attach
        ( height, width, B.Grid(), colAlign, rowAlign, 0, B.LDim(), B.Root() );
    }
}

template<typename T,Dist U,Dist V>
inline DistMatrix<T,U,V> View
( DistMatrix<T,U,V>& B, Int i, Int j, Int height, Int width )
{
    DistMatrix<T,U,V> A(B.Grid());
    View( A, B, i, j, height, width );
    return A;
}

template<typename T>
inline void LockedView
( Matrix<T>& A, const Matrix<T>& B,
  Int i, Int j, Int height, Int width )
{
    DEBUG_ONLY(
        CallStackEntry cse("LockedView");
        if( i < 0 || j < 0 )
            LogicError("Indices must be non-negative");
        if( height < 0 || width < 0 )
            LogicError("Height and width must be non-negative");
        if( (i+height) > B.Height() || (j+width) > B.Width() )
            LogicError
            ("Trying to view outside of a Matrix: (",i,",",j,") up to (",
             i+height-1,",",j+width-1,") of ",B.Height()," x ",B.Width(),
             " Matrix");
    )
    A.LockedAttach( height, width, B.LockedBuffer(i,j), B.LDim() );
}

template<typename T>
inline Matrix<T> LockedView
( const Matrix<T>& B, Int i, Int j, Int height, Int width )
{
    Matrix<T> A;
    LockedView( A, B, i, j, height, width );
    return A;
}

template<typename T>
inline void LockedView
( AbstractDistMatrix<T>& A, const AbstractDistMatrix<T>& B,
  Int i, Int j, Int height, Int width )
{
    DEBUG_ONLY(
        CallStackEntry cse("LockedView");
        AssertSameDist( A.DistData(), B.DistData() );
        B.AssertValidSubmatrix( i, j, height, width );
    )
    const Int colAlign = (B.ColAlign()+i) % B.ColStride();
    const Int rowAlign = (B.RowAlign()+j) % B.RowStride();
    if( B.Participating() )
    {
        const Int iLoc = Length( i, B.ColShift(), B.ColStride() );
        const Int jLoc = Length( j, B.RowShift(), B.RowStride() );
        A.LockedAttach
        ( height, width, B.Grid(), colAlign, rowAlign, 
          B.LockedBuffer(iLoc,jLoc), B.LDim(), B.Root() );
    }
    else
    {
        A.LockedAttach
        ( height, width, B.Grid(), colAlign, rowAlign, 0, B.LDim(), B.Root() );
    }
}

template<typename T,Dist U,Dist V>
inline DistMatrix<T,U,V> LockedView
( const DistMatrix<T,U,V>& B, Int i, Int j, Int height, Int width )
{
    DistMatrix<T,U,V> A(B.Grid());
    LockedView( A, B, i, j, height, width );
    return A;
}

template<typename T>
inline void ViewRange
( Matrix<T>& A, Matrix<T>& B, Int iBeg, Int jBeg, Int iEnd, Int jEnd )
{ 
    DEBUG_ONLY(CallStackEntry cse("ViewRange"))
    View( A, B, iBeg, jBeg, iEnd-iBeg, jEnd-jBeg ); 
}

template<typename T>
inline Matrix<T> ViewRange
( Matrix<T>& B, Int iBeg, Int jBeg, Int iEnd, Int jEnd )
{
    DEBUG_ONLY(CallStackEntry cse("ViewRange"))
    return View( B, iBeg, jBeg, iEnd-iBeg, jEnd-jBeg ); 
}

template<typename T>
inline void ViewRange
( AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B,
  Int iBeg, Int jBeg, Int iEnd, Int jEnd )
{
    DEBUG_ONLY(CallStackEntry cse("ViewRange"))
    View( A, B, iBeg, jBeg, iEnd-iBeg, jEnd-jBeg ); 
}

template<typename T,Dist U,Dist V>
inline DistMatrix<T,U,V> ViewRange
( DistMatrix<T,U,V>& B, Int iBeg, Int jBeg, Int iEnd, Int jEnd )
{
    DEBUG_ONLY(CallStackEntry cse("ViewRange"))
    return View( B, iBeg, jBeg, iEnd-iBeg, jEnd-jBeg ); 
} 

template<typename T>
inline void LockedViewRange
( Matrix<T>& A, const Matrix<T>& B, Int iBeg, Int jBeg, Int iEnd, Int jEnd )
{
    DEBUG_ONLY(CallStackEntry cse("LockedViewRange"))
    LockedView( A, B, iBeg, jBeg, iEnd-iBeg, jEnd-jBeg ); 
}

template<typename T>
inline Matrix<T> LockedViewRange
( const Matrix<T>& B, Int iBeg, Int jBeg, Int iEnd, Int jEnd )
{
    DEBUG_ONLY(CallStackEntry cse("LockedViewRange"))
    return LockedView( B, iBeg, jBeg, iEnd-iBeg, jEnd-jBeg ); 
}

template<typename T>
inline void LockedViewRange
( AbstractDistMatrix<T>& A, const AbstractDistMatrix<T>& B,
  Int iBeg, Int jBeg, Int iEnd, Int jEnd )
{
    DEBUG_ONLY(CallStackEntry cse("LockedViewRange"))
    LockedView( A, B, iBeg, jBeg, iEnd-iBeg, jEnd-jBeg ); 
}

template<typename T,Dist U,Dist V>
inline DistMatrix<T,U,V> LockedViewRange
( const DistMatrix<T,U,V>& B, Int iBeg, Int jBeg, Int iEnd, Int jEnd )
{
    DEBUG_ONLY(CallStackEntry cse("LockedViewRange"))
    return LockedView( B, iBeg, jBeg, iEnd-iBeg, jEnd-jBeg ); 
}

template<typename T>
inline void View1x2( Matrix<T>& A, Matrix<T>& BL, Matrix<T>& BR )
{
    DEBUG_ONLY(
        CallStackEntry cse("View1x2");
        if( BL.Locked() || BR.Locked() )
            LogicError("Cannot grab an unlocked view of a locked matrix");
        if( BL.Height() != BR.Height() )
            LogicError("1x2 must have consistent height to combine");
        if( BL.LDim() != BR.LDim() )
            LogicError("1x2 must have consistent ldims to combine");
        if( BR.Buffer() != (BL.Buffer()+BL.LDim()*BL.Width()) )
            LogicError("1x2 must have contiguous memory");
    )
    A.Attach( BL.Height(), BL.Width()+BR.Width(), BL.Buffer(), BL.LDim() );
}

template<typename T>
inline Matrix<T> View1x2( Matrix<T>& BL, Matrix<T>& BR )
{
    Matrix<T> A;
    View1x2( A, BL, BR );
    return A;
}

template<typename T>
inline void View1x2
( AbstractDistMatrix<T>& A, 
  AbstractDistMatrix<T>& BL, AbstractDistMatrix<T>& BR )
{
    DEBUG_ONLY(
        CallStackEntry cse("View1x2");
        AssertSameDist( A.DistData(), BL.DistData() );
        AssertSameDist( BL.DistData(), BR.DistData() );
        AssertConforming1x2( BL, BR );
        BL.AssertSameGrid( BR.Grid() );
    )
    A.Attach
    ( BL.Height(), BL.Width()+BR.Width(), BL.Grid(), 
      BL.ColAlign(), BL.RowAlign(), BL.Buffer(), BL.LDim(), BL.Root() );
}

template<typename T,Dist U,Dist V>
inline DistMatrix<T,U,V> View1x2( DistMatrix<T,U,V>& BL, DistMatrix<T,U,V>& BR )
{
    DistMatrix<T,U,V> A(BL.Grid());
    View1x2( A, BL, BR );
    return A;
}

template<typename T>
inline void LockedView1x2
( Matrix<T>& A, const Matrix<T>& BL, const Matrix<T>& BR )
{
    DEBUG_ONLY(
        CallStackEntry cse("LockedView1x2");
        if( BL.Height() != BR.Height() )
            LogicError("1x2 must have consistent height to combine");
        if( BL.LDim() != BR.LDim() )
            LogicError("1x2 must have consistent ldims to combine");
        if( BR.LockedBuffer() != (BL.LockedBuffer()+BL.LDim()*BL.Width()) )
            LogicError("1x2 must have contiguous memory");
    )
    A.LockedAttach
    ( BL.Height(), BL.Width()+BR.Width(), BL.LockedBuffer(), BL.LDim() );
}

template<typename T>
inline Matrix<T> LockedView1x2( const Matrix<T>& BL, const Matrix<T>& BR )
{
    Matrix<T> A;
    LockedView1x2( A, BL, BR );
    return A;
}

template<typename T>
inline void LockedView1x2
(       AbstractDistMatrix<T>& A,
  const AbstractDistMatrix<T>& BL, const AbstractDistMatrix<T>& BR )
{
    DEBUG_ONLY(
        CallStackEntry cse("LockedView1x2");
        AssertSameDist( A.DistData(), BL.DistData() );
        AssertSameDist( BL.DistData(), BR.DistData() );
        AssertConforming1x2( BL, BR );
        BL.AssertSameGrid( BR.Grid() );
    )
    A.LockedAttach
    ( BL.Height(), BL.Width()+BR.Width(), BL.Grid(), 
      BL.ColAlign(), BL.RowAlign(), BL.LockedBuffer(), BL.LDim(), BL.Root() );
}

template<typename T,Dist U,Dist V>
inline DistMatrix<T,U,V> LockedView1x2
( const DistMatrix<T,U,V>& BL, const DistMatrix<T,U,V>& BR )
{
    DistMatrix<T,U,V> A(BL.Grid());
    LockedView1x2( A, BL, BR );
    return A;
}

template<typename T>
inline void View2x1( Matrix<T>& A, Matrix<T>& BT, Matrix<T>& BB )
{
    DEBUG_ONLY(
        CallStackEntry cse("View2x1");
        if( BT.Locked() || BB.Locked() )
            LogicError("Cannot grab an unlocked view of a locked matrix");
        if( BT.Width() != BB.Width() )
            LogicError("2x1 must have consistent width to combine");
        if( BT.LDim() != BB.LDim() )
            LogicError("2x1 must have consistent ldim to combine");
        if( BB.Buffer() != (BT.Buffer() + BT.Height()) )
            LogicError("2x1 must have contiguous memory");
    )
    A.Attach( BT.Height()+BB.Height(), BT.Width(), BT.Buffer(), BT.LDim() );
}

template<typename T>
inline Matrix<T> View2x1( Matrix<T>& BT, Matrix<T>& BB )
{
    Matrix<T> A;
    View2x1( A, BT, BB );
    return A;
}

template<typename T>
inline void View2x1
( AbstractDistMatrix<T>& A, 
  AbstractDistMatrix<T>& BT, AbstractDistMatrix<T>& BB )
{
    DEBUG_ONLY(
        CallStackEntry cse("View2x1");
        AssertSameDist( A.DistData(), BT.DistData() );
        AssertSameDist( BT.DistData(), BB.DistData() );
        AssertConforming2x1( BT, BB );
        BT.AssertSameGrid( BB.Grid() );
    )
    A.Attach
    ( BT.Height()+BB.Height(), BT.Width(), BT.Grid(), 
      BT.ColAlign(), BT.RowAlign(), BT.Buffer(), BT.LDim(), BT.Root() );
}

template<typename T,Dist U,Dist V>
inline DistMatrix<T,U,V> View2x1( DistMatrix<T,U,V>& BT, DistMatrix<T,U,V>& BB )
{
    DistMatrix<T,U,V> A(BT.Grid());
    View2x1( A, BT, BB );
    return A;
}

template<typename T>
inline void LockedView2x1
( Matrix<T>& A, const Matrix<T>& BT, const Matrix<T>& BB )
{
    DEBUG_ONLY(
        CallStackEntry cse("LockedView2x1");
        if( BT.Width() != BB.Width() )
            LogicError("2x1 must have consistent width to combine");
        if( BT.LDim() != BB.LDim() )
            LogicError("2x1 must have consistent ldim to combine");
        if( BB.LockedBuffer() != (BT.LockedBuffer() + BT.Height()) )
            LogicError("2x1 must have contiguous memory");
    )
    A.LockedAttach
    ( BT.Height()+BB.Height(), BT.Width(), BT.LockedBuffer(), BT.LDim() );
}

template<typename T>
inline Matrix<T> LockedView2x1( const Matrix<T>& BT, const Matrix<T>& BB )
{
    Matrix<T> A;
    LockedView2x1( A, BT, BB );
    return A;
}

template<typename T>
inline void LockedView2x1
(       AbstractDistMatrix<T>& A,
  const AbstractDistMatrix<T>& BT,
  const AbstractDistMatrix<T>& BB )
{
    DEBUG_ONLY(
        CallStackEntry cse("LockedView2x1");
        AssertSameDist( A.DistData(), BT.DistData() );
        AssertSameDist( BT.DistData(), BB.DistData() );
        AssertConforming2x1( BT, BB );
        BT.AssertSameGrid( BB.Grid() );
    )
    A.LockedAttach
    ( BT.Height()+BB.Height(), BT.Width(), BT.Grid(), 
      BT.ColAlign(), BT.RowAlign(), BT.LockedBuffer(), BT.LDim(), BT.Root() );
}

template<typename T,Dist U,Dist V>
inline DistMatrix<T,U,V> LockedView2x1
( const DistMatrix<T,U,V>& BT, const DistMatrix<T,U,V>& BB )
{
    DistMatrix<T,U,V> A(BT.Grid());
    LockedView2x1( A, BT, BB );
    return A;
}

template<typename T>
inline void View2x2
( Matrix<T>& A,
  Matrix<T>& BTL, Matrix<T>& BTR,
  Matrix<T>& BBL, Matrix<T>& BBR )
{
    DEBUG_ONLY(
        CallStackEntry cse("View2x2");
        if( BTL.Locked() || BTR.Locked() || BBL.Locked() || BBR.Locked() )
            LogicError("Cannot grab an unlocked view of a locked matrix");
        if( BTL.Width() != BBL.Width()   ||
            BTR.Width() != BBR.Width()   ||
            BTL.Height() != BTR.Height() ||
            BBL.Height() != BBR.Height()   )
            LogicError("2x2 must conform to combine");
        if( BTL.LDim() != BTR.LDim() ||
            BTR.LDim() != BBL.LDim() ||
            BBL.LDim() != BBR.LDim()   )
            LogicError("2x2 must have consistent ldims to combine");
        if( BBL.Buffer() != (BTL.Buffer() + BTL.Height()) ||
            BBR.Buffer() != (BTR.Buffer() + BTR.Height()) ||
            BTR.Buffer() != (BTL.Buffer() + BTL.LDim()*BTL.Width()) )
            LogicError("2x2 must have contiguous memory");
    )
    A.Attach
    ( BTL.Height()+BBL.Height(), BTL.Width()+BTR.Width(), 
      BTL.Buffer(), BTL.LDim() );
}

template<typename T>
inline Matrix<T> View2x2
( Matrix<T>& BTL, Matrix<T>& BTR,
  Matrix<T>& BBL, Matrix<T>& BBR )
{
    Matrix<T> A;
    View2x2( A, BTL, BTR, BBL, BBR );
    return A;
}

template<typename T>
inline void View2x2
( AbstractDistMatrix<T>& A,
  AbstractDistMatrix<T>& BTL, AbstractDistMatrix<T>& BTR,
  AbstractDistMatrix<T>& BBL, AbstractDistMatrix<T>& BBR )
{
    DEBUG_ONLY(
        CallStackEntry cse("View2x2");
        AssertSameDist( A.DistData(), BTL.DistData() );
        AssertSameDist( BTL.DistData(), BTR.DistData() );
        AssertSameDist( BTR.DistData(), BBL.DistData() );
        AssertSameDist( BBL.DistData(), BBR.DistData() );
        AssertConforming2x2( BTL, BTR, BBL, BBR );
        BTL.AssertSameGrid( BTR.Grid() );
        BTL.AssertSameGrid( BBL.Grid() );
        BTL.AssertSameGrid( BBR.Grid() );
    )
    A.Attach
    ( BTL.Height()+BBL.Height(), BTL.Width()+BTR.Width(), BTL.Grid(),
      BTL.ColAlign(), BTL.RowAlign(), BTL.Buffer(), BTL.LDim(), BTL.Root() );
}

template<typename T,Dist U,Dist V>
inline DistMatrix<T,U,V> View2x2
( DistMatrix<T,U,V>& BTL, DistMatrix<T,U,V>& BTR,
  DistMatrix<T,U,V>& BBL, DistMatrix<T,U,V>& BBR )
{
    DistMatrix<T,U,V> A(BTL.Grid());
    View2x2( A, BTL, BTR, BBL, BBR );
    return A;
}

template<typename T>
inline void LockedView2x2
(       Matrix<T>& A,
  const Matrix<T>& BTL, const Matrix<T>& BTR,
  const Matrix<T>& BBL, const Matrix<T>& BBR )
{
    DEBUG_ONLY(
        CallStackEntry cse("LockedView2x2");
        if( BTL.Width() != BBL.Width()   ||
            BTR.Width() != BBR.Width()   ||
            BTL.Height() != BTR.Height() ||
            BBL.Height() != BBR.Height()   )
            LogicError("2x2 must conform to combine");
        if( BTL.LDim() != BTR.LDim() ||
            BTR.LDim() != BBL.LDim() ||
            BBL.LDim() != BBR.LDim()   )
            LogicError("2x2 must have consistent ldims to combine");
        if( BBL.LockedBuffer() != (BTL.LockedBuffer()+BTL.Height()) ||
            BBR.LockedBuffer() != (BTR.LockedBuffer()+BTR.Height()) ||
            BTR.LockedBuffer() != (BTL.LockedBuffer()+BTL.LDim()*BTL.Width()) )
            LogicError("2x2 must have contiguous memory");
    )
    A.LockedAttach
    ( BTL.Height()+BBL.Height(), BTL.Width()+BTR.Width(), 
      BTL.LockedBuffer(), BTL.LDim() );
}

template<typename T>
inline Matrix<T> LockedView2x2
( const Matrix<T>& BTL, const Matrix<T>& BTR,
  const Matrix<T>& BBL, const Matrix<T>& BBR )
{
    Matrix<T> A;
    LockedView2x2( A, BTL, BTR, BBL, BBR );
    return A;
}

template<typename T>
inline void LockedView2x2
(       AbstractDistMatrix<T>& A,
  const AbstractDistMatrix<T>& BTL, const AbstractDistMatrix<T>& BTR,
  const AbstractDistMatrix<T>& BBL, const AbstractDistMatrix<T>& BBR )
{
    DEBUG_ONLY(
        CallStackEntry cse("LockedView2x2");
        AssertSameDist( A.DistData(), BTL.DistData() );
        AssertSameDist( BTL.DistData(), BTR.DistData() );
        AssertSameDist( BTR.DistData(), BBL.DistData() );
        AssertSameDist( BBL.DistData(), BBR.DistData() );
        AssertConforming2x2( BTL, BTR, BBL, BBR );
        BTL.AssertSameGrid( BTR.Grid() );
        BTL.AssertSameGrid( BBL.Grid() );
        BTL.AssertSameGrid( BBR.Grid() );
    )
    A.LockedAttach
    ( BTL.Height()+BBL.Height(), BTL.Width()+BTR.Width(), BTL.Grid(),
      BTL.ColAlign(), BTL.RowAlign(), BTL.LockedBuffer(), BTL.LDim(),
      BTL.Root() );
}

template<typename T,Dist U,Dist V>
inline DistMatrix<T,U,V> LockedView2x2
( const DistMatrix<T,U,V>& BTL, const DistMatrix<T,U,V>& BTR,
  const DistMatrix<T,U,V>& BBL, const DistMatrix<T,U,V>& BBR )
{
    DistMatrix<T,U,V> A(BTL.Grid());
    LockedView2x2( A, BTL, BTR, BBL, BBR );
    return A;
}

} // namespace El

#endif // ifndef EL_VIEWS_VIEW_HPP
