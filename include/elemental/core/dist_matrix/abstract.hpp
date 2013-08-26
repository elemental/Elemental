/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_CORE_DISTMATRIX_ABSTRACT_DECL_HPP
#define ELEM_CORE_DISTMATRIX_ABSTRACT_DECL_HPP

namespace elem {

#ifndef RELEASE
template<typename T>
void AssertConforming1x2
( const AbstractDistMatrix<T>& AL, const AbstractDistMatrix<T>& AR );

template<typename T>
void AssertConforming2x1
( const AbstractDistMatrix<T>& AT, const AbstractDistMatrix<T>& AB );

template<typename T>
void AssertConforming2x2
( const AbstractDistMatrix<T>& ATL, const AbstractDistMatrix<T>& ATR,
  const AbstractDistMatrix<T>& ABL, const AbstractDistMatrix<T>& ABR );
#endif // ifndef RELEASE

template<typename T> 
class AbstractDistMatrix
{
public:
    virtual ~AbstractDistMatrix();

    //-----------------------------------------------------------------------//
    // Routines that do NOT need to be implemented in derived classes        //
    //-----------------------------------------------------------------------//

    // Move constructor
    AbstractDistMatrix( AbstractDistMatrix<T>&& A );

    // Move operator=
    AbstractDistMatrix<T>& operator=( AbstractDistMatrix<T>&& A );

#ifndef RELEASE
    void AssertNotLocked() const;
    void AssertNotStoringData() const;
    void AssertValidEntry( Int i, Int j ) const;
    void AssertValidSubmatrix
    ( Int i, Int j, Int height, Int width ) const;
    void AssertSameGrid( const elem::Grid& grid ) const;
    void AssertSameSize( Int height, Int width ) const;
#endif // ifndef RELEASE

    //
    // Basic information
    //

    Int Height() const;
    Int Width() const;
    Int DiagonalLength( Int offset=0 ) const;
    Int LocalHeight() const;
    Int LocalWidth() const;
    Int LDim() const;
    size_t AllocatedMemory() const;

    const elem::Grid& Grid() const;

          T* Buffer( Int iLoc=0, Int jLoc=0 );
    const T* LockedBuffer( Int iLoc=0, Int jLoc=0 ) const;

          elem::Matrix<T>& Matrix();
    const elem::Matrix<T>& LockedMatrix() const;

    //
    // Alignments
    //

    void FreeAlignments();
    bool ConstrainedColAlignment() const;
    bool ConstrainedRowAlignment() const;
    Int ColAlignment() const;
    Int RowAlignment() const;
    Int ColShift() const;
    Int RowShift() const;

    void Align( Int colAlign, Int rowAlign );
    void AlignCols( Int colAlign );
    void AlignRows( Int rowAlign );

    //
    // Local entry manipulation
    //

    T GetLocal( Int iLoc, Int jLoc ) const;
    void SetLocal( Int iLoc, Int jLoc, T alpha );
    void UpdateLocal( Int iLoc, Int jLoc, T alpha );

    //
    // Though the following routines are meant for complex data, all but two
    // logically apply to real data.
    //

    BASE(T) GetRealPart( Int i, Int j ) const;
    BASE(T) GetImagPart( Int i, Int j ) const;
    BASE(T) GetLocalRealPart( Int iLoc, Int jLoc ) const;
    BASE(T) GetLocalImagPart( Int iLoc, Int jLoc ) const;
    void SetLocalRealPart( Int iLoc, Int jLoc, BASE(T) alpha );
    void UpdateLocalRealPart( Int iLoc, Int jLoc, BASE(T) alpha );
    // Only valid for complex data
    void SetLocalImagPart( Int iLoc, Int jLoc, BASE(T) alpha );
    void UpdateLocalImagPart( Int iLoc, Int jLoc, BASE(T) alpha );

    //
    // Viewing 
    //

    bool Viewing() const;
    bool Locked()  const;

    //
    // Utilities
    //

    void Empty();
    void EmptyData();
    void SetGrid( const elem::Grid& grid );

    //------------------------------------------------------------------------//
    // Routines that can be overridden in derived classes                     //
    //------------------------------------------------------------------------//

    virtual void Swap( AbstractDistMatrix<T>& A );

    virtual bool Participating() const;
    virtual void AlignWith( const elem::DistData& data );
    virtual void AlignWith( const AbstractDistMatrix<T>& A );
    virtual void AlignColsWith( const elem::DistData& data );
    virtual void AlignColsWith( const AbstractDistMatrix<T>& A );
    virtual void AlignRowsWith( const elem::DistData& data );
    virtual void AlignRowsWith( const AbstractDistMatrix<T>& A );

    virtual void MakeConsistent();

    //------------------------------------------------------------------------//
    // Routines that MUST be implemented in non-abstract derived classes      //
    //------------------------------------------------------------------------//

    //
    // Basic information
    //

    virtual elem::DistData DistData() const = 0;

    // So that the local row indices are given by
    //   A.ColShift():A.ColStride():A.Height()
    virtual Int ColStride() const = 0; 
    // So that the local column indices are given by
    //   A.RowShift():A.RowStride():A.Width()
    virtual Int RowStride() const = 0;
    virtual Int ColRank() const = 0;
    virtual Int RowRank() const = 0;

    //
    // Entry manipulation
    //

    virtual T Get( Int i, Int j ) const = 0;
    virtual void Set( Int i, Int j, T alpha ) = 0;
    virtual void Update( Int i, Int j, T alpha ) = 0;

    //
    // Though the following routines are meant for complex data, all but two
    // logically applies to real data.
    //

    virtual void SetRealPart( Int i, Int j, BASE(T) alpha ) = 0;
    // Only valid for complex data
    virtual void SetImagPart( Int i, Int j, BASE(T) alpha ) = 0;
    virtual void UpdateRealPart( Int i, Int j, BASE(T) alpha ) = 0;
    // Only valid for complex data
    virtual void UpdateImagPart( Int i, Int j, BASE(T) alpha ) = 0;

    //
    // Utilities
    //
    
    virtual void ResizeTo( Int height, Int width ) = 0;
    virtual void ResizeTo( Int height, Int width, Int ldim ) = 0;

protected:
    ViewType viewType_;
    Int height_, width_;
    Memory<T> auxMemory_;
    elem::Matrix<T> matrix_;
    
    bool constrainedColAlignment_, constrainedRowAlignment_;
    Int colAlignment_, rowAlignment_;
    Int colShift_, rowShift_;
    const elem::Grid* grid_;

    // Build around a particular grid
    AbstractDistMatrix( const elem::Grid& g );

    void SetShifts();
    void SetColShift();
    void SetRowShift();
    void SetGrid();

    void ComplainIfReal() const;

    void SetAlignmentsAndResize
    ( Int colAlign, Int rowAlign, Int height, Int width );
    void ForceAlignmentsAndResize
    ( Int colAlign, Int rowAlign, Int height, Int width );

    void SetColAlignmentAndResize
    ( Int colAlign, Int height, Int width );
    void ForceColAlignmentAndResize
    ( Int colAlign, Int height, Int width );

    void SetRowAlignmentAndResize
    ( Int rowAlign, Int height, Int width );
    void ForceRowAlignmentAndResize
    ( Int rowAlign, Int height, Int width );

#ifndef SWIG
    template<typename S,Distribution U,Distribution V> 
    friend void View( DistMatrix<S,U,V>& A, DistMatrix<S,U,V>& B );
    template<typename S,Distribution U,Distribution V> 
    friend void LockedView( DistMatrix<S,U,V>& A, const DistMatrix<S,U,V>& B );
    template<typename S,Distribution U,Distribution V> 
    friend void View
    ( DistMatrix<S,U,V>& A, DistMatrix<S,U,V>& B,
      Int i, Int j, Int height, Int width );
    template<typename S,Distribution U,Distribution V> 
    friend void LockedView
    ( DistMatrix<S,U,V>& A, const DistMatrix<S,U,V>& B,
      Int i, Int j, Int height, Int width );
    template<typename S,Distribution U,Distribution V> 
    friend void View1x2
    ( DistMatrix<S,U,V>& A, DistMatrix<S,U,V>& BL, DistMatrix<S,U,V>& BR );
    template<typename S,Distribution U,Distribution V> 
    friend void LockedView1x2
    (       DistMatrix<S,U,V>& A,
      const DistMatrix<S,U,V>& BL, const DistMatrix<S,U,V>& BR );
    template<typename S,Distribution U,Distribution V> 
    friend void View2x1
    ( DistMatrix<S,U,V>& A, DistMatrix<S,U,V>& BT, DistMatrix<S,U,V>& BB );
    template<typename S,Distribution U,Distribution V> 
    friend void LockedView2x1
    (       DistMatrix<S,U,V>& A,
      const DistMatrix<S,U,V>& BT, const DistMatrix<S,U,V>& BB );
    template<typename S,Distribution U,Distribution V> 
    friend void View2x2
    ( DistMatrix<S,U,V>& A,
      DistMatrix<S,U,V>& BTL, DistMatrix<S,U,V>& BTR,
      DistMatrix<S,U,V>& BBL, DistMatrix<S,U,V>& BBR );
    template<typename S,Distribution U,Distribution V> 
    friend void LockedView2x2
    (       DistMatrix<S,U,V>& A,
      const DistMatrix<S,U,V>& BTL, const DistMatrix<S,U,V>& BTR,
      const DistMatrix<S,U,V>& BBL, const DistMatrix<S,U,V>& BBR );

    template<typename S,Distribution U,Distribution V>
    friend class DistMatrix;
#endif // ifndef SWIG
};

} // namespace elem

#endif // ifndef ELEM_CORE_DISTMATRIX_ABSTRACT_DECL_HPP
