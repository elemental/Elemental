/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_VIEW_HPP
#define EL_VIEW_HPP

namespace El {

// View an entire matrix
// =====================

// (Sequential) matrix
// -------------------

template<typename T>
inline void View( Matrix<T>& A, Matrix<T>& B )
{
    DEBUG_ONLY(CallStackEntry cse("View"))
    A.Attach( B.Height(), B.Width(), B.Buffer(), B.LDim() );
}

template<typename T>
inline void LockedView( Matrix<T>& A, const Matrix<T>& B )
{
    DEBUG_ONLY(CallStackEntry cse("LockedView"))
    A.LockedAttach( B.Height(), B.Width(), B.LockedBuffer(), B.LDim() );
}

template<typename T>
inline Matrix<T> View( Matrix<T>& B )
{
    Matrix<T> A;
    View( A, B );
    return A;
}

template<typename T>
inline Matrix<T> LockedView( const Matrix<T>& B )
{
    Matrix<T> A;
    LockedView( A, B );
    return A;
}

// DistMatrix
// ----------

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

// Return by value
// ^^^^^^^^^^^^^^^

template<typename T,Dist U,Dist V>
inline DistMatrix<T,U,V> View( DistMatrix<T,U,V>& B )
{
    DistMatrix<T,U,V> A(B.Grid());
    View( A, B );
    return A;
}

template<typename T,Dist U,Dist V>
inline DistMatrix<T,U,V> LockedView( const DistMatrix<T,U,V>& B )
{
    DistMatrix<T,U,V> A(B.Grid());
    LockedView( A, B );
    return A;
}

// BlockDistMatrix
// ---------------

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

// View a contiguous submatrix
// ===========================

// (Sequential) Matrix
// -------------------

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
inline void View
( Matrix<T>& A, Matrix<T>& B, Range<Int> I, Range<Int> J )
{ View( A, B, I.beg, J.beg, I.end-I.beg, J.end-J.beg ); }

template<typename T>
inline void LockedView
( Matrix<T>& A, const Matrix<T>& B, Range<Int> I, Range<Int> J )
{ LockedView( A, B, I.beg, J.beg, I.end-I.beg, J.end-J.beg ); }

// Return by value
// ^^^^^^^^^^^^^^^

template<typename T>
inline Matrix<T> View( Matrix<T>& B, Int i, Int j, Int height, Int width )
{
    Matrix<T> A;
    View( A, B, i, j, height, width );
    return A;
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
inline Matrix<T> View
( Matrix<T>& B, Range<Int> I, Range<Int> J )
{ return View( B, I.beg, J.beg, I.end-I.beg, J.end-J.beg ); }

template<typename T>
inline Matrix<T> LockedView
( const Matrix<T>& B, Range<Int> I, Range<Int> J )
{ return LockedView( B, I.beg, J.beg, I.end-I.beg, J.end-J.beg ); }

// DistMatrix
// ----------

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

template<typename T>
inline void View
( AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B, 
  Range<Int> I, Range<Int> J )
{ View( A, B, I.beg, J.beg, I.end-I.beg, J.end-J.beg ); }

template<typename T>
inline void LockedView
( AbstractDistMatrix<T>& A, const AbstractDistMatrix<T>& B, 
  Range<Int> I, Range<Int> J )
{ LockedView( A, B, I.beg, J.beg, I.end-I.beg, J.end-J.beg ); }

// Return by value
// ^^^^^^^^^^^^^^^

template<typename T,Dist U,Dist V>
inline DistMatrix<T,U,V> View
( DistMatrix<T,U,V>& B, Int i, Int j, Int height, Int width )
{
    DistMatrix<T,U,V> A(B.Grid());
    View( A, B, i, j, height, width );
    return A;
}

template<typename T,Dist U,Dist V>
inline DistMatrix<T,U,V> LockedView
( const DistMatrix<T,U,V>& B, Int i, Int j, Int height, Int width )
{
    DistMatrix<T,U,V> A(B.Grid());
    LockedView( A, B, i, j, height, width );
    return A;
}

template<typename T,Dist U,Dist V>
inline DistMatrix<T,U,V> View
( DistMatrix<T,U,V>& B, Range<Int> I, Range<Int> J )
{ return View( B, I.beg, J.beg, I.end-I.beg, J.end-J.beg ); }
 
template<typename T,Dist U,Dist V>
inline DistMatrix<T,U,V> LockedView
( const DistMatrix<T,U,V>& B, Range<Int> I, Range<Int> J )
{ return LockedView( B, I.beg, J.beg, I.end-I.beg, J.end-J.beg ); }

} // namespace El

#endif // ifndef EL_VIEW_HPP
