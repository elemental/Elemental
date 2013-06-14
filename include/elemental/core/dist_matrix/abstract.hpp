/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef CORE_DISTMATRIX_ABSTRACT_DECL_HPP
#define CORE_DISTMATRIX_ABSTRACT_DECL_HPP

namespace elem {

#ifndef RELEASE
template<typename T,typename Int>
void AssertConforming1x2
( const AbstractDistMatrix<T,Int>& AL, const AbstractDistMatrix<T,Int>& AR );

template<typename T,typename Int>
void AssertConforming2x1
( const AbstractDistMatrix<T,Int>& AT,
  const AbstractDistMatrix<T,Int>& AB );

template<typename T,typename Int>
void AssertConforming2x2
( const AbstractDistMatrix<T,Int>& ATL, const AbstractDistMatrix<T,Int>& ATR,
  const AbstractDistMatrix<T,Int>& ABL, const AbstractDistMatrix<T,Int>& ABR );
#endif // ifndef RELEASE

template<typename T,typename Int> 
class AbstractDistMatrix
{
public:
    virtual ~AbstractDistMatrix();

    //-----------------------------------------------------------------------//
    // Routines that do NOT need to be implemented in derived classes        //
    //-----------------------------------------------------------------------//

#ifndef RELEASE
    void AssertNotLocked() const;
    void AssertNotStoringData() const;
    void AssertValidEntry( Int i, Int j ) const;
    void AssertValidSubmatrix( Int i, Int j, Int height, Int width ) const;
    void AssertFreeColAlignment() const;
    void AssertFreeRowAlignment() const;
    void AssertSameGrid( const elem::Grid& grid ) const;
    void AssertSameSize( int height, int width ) const;
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

          T* Buffer( Int iLocal=0, Int jLocal=0 );
    const T* LockedBuffer( Int iLocal=0, Int jLocal=0 ) const;

          elem::Matrix<T,Int>& Matrix();
    const elem::Matrix<T,Int>& LockedMatrix() const;

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

    void Align( Int colAlignment, Int rowAlignment );
    void AlignCols( Int colAlignment );
    void AlignRows( Int rowAlignment );

    //
    // Local entry manipulation
    //

    T GetLocal( Int iLocal, Int jLocal ) const;
    void SetLocal( Int iLocal, Int jLocal, T alpha );
    void UpdateLocal( Int iLocal, Int jLocal, T alpha );

    //
    // Though the following routines are meant for complex data, all but two
    // logically applies to real data.
    //

    BASE(T) GetLocalRealPart( Int iLocal, Int jLocal ) const;
    BASE(T) GetLocalImagPart( Int iLocal, Int jLocal ) const;
    void SetLocalRealPart( Int iLocal, Int jLocal, BASE(T) alpha );
    void UpdateLocalRealPart( Int iLocal, Int jLocal, BASE(T) alpha );
    // Only valid for complex data
    void SetLocalImagPart( Int iLocal, Int jLocal, BASE(T) alpha );
    void UpdateLocalImagPart( Int iLocal, Int jLocal, BASE(T) alpha );

    //
    // Viewing 
    //

    bool Viewing() const;
    bool Locked() const;

    //
    // Utilities
    //

    void Empty();
    void EmptyData();
    void SetGrid( const elem::Grid& grid );

    //------------------------------------------------------------------------//
    // Routines that can be overridden in derived classes                     //
    //------------------------------------------------------------------------//

    virtual bool Participating() const;
    virtual void AlignWith( const elem::DistData<Int>& data );
    virtual void AlignWith( const AbstractDistMatrix<T,Int>& A );
    virtual void AlignColsWith( const elem::DistData<Int>& data );
    virtual void AlignColsWith( const AbstractDistMatrix<T,Int>& A );
    virtual void AlignRowsWith( const elem::DistData<Int>& data );
    virtual void AlignRowsWith( const AbstractDistMatrix<T,Int>& A );

    //------------------------------------------------------------------------//
    // Routines that MUST be implemented in non-abstract derived classes      //
    //------------------------------------------------------------------------//

    //
    // Basic information
    //

    virtual elem::DistData<Int> DistData() const = 0;

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

    virtual BASE(T) GetRealPart( Int i, Int j ) const = 0;
    virtual BASE(T) GetImagPart( Int i, Int j ) const = 0;
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

protected:
    bool viewing_, locked_;
    Int height_, width_;
    Memory<T> auxMemory_;
    elem::Matrix<T,Int> matrix_;
    
    bool constrainedColAlignment_, constrainedRowAlignment_;
    Int colAlignment_, rowAlignment_;
    Int colShift_, rowShift_;
    const elem::Grid* grid_;

    // Initialize with particular local dimensions
    AbstractDistMatrix
    ( Int height, Int width,
      bool constrainedColAlignment, bool constrainedRowAlignment,
      Int colAlignment, Int rowAlignment,
      Int colShift, Int rowShift,
      Int localHeight, Int localWidth,
      const elem::Grid& g );

    // Initialize with particular local dimensions and local leading dimensions
    AbstractDistMatrix
    ( Int height, Int width,
      bool constrainedColAlignment, bool constrainedRowAlignment,
      Int colAlignment, Int rowAlignment,
      Int colShift, Int rowShift,
      Int localHeight, Int localWidth,
      Int ldim,
      const elem::Grid& g );

    // View a constant distributed matrix's buffer
    AbstractDistMatrix
    ( Int height, Int width,
      Int colAlignment, Int rowAlignment,
      Int colShift, Int rowShift,
      Int localHeight, Int localWidth,
      const T* buffer,
      Int ldim,
      const elem::Grid& g );

    // View a mutable distributed matrix's buffer
    AbstractDistMatrix
    ( Int height, Int width,
      Int colAlignment, Int rowAlignment,
      Int colShift, Int rowShift,
      Int localHeight, Int localWidth,
      T* buffer,
      Int ldim,
      const elem::Grid& g );

protected:
    void SetShifts();
    void SetColShift();
    void SetRowShift();
    void SetGrid();

#ifndef SWIG
    template<typename S,Distribution U,Distribution V,typename Ord> 
    friend void View
    ( DistMatrix<S,U,V,Ord>& A, DistMatrix<S,U,V,Ord>& B );
    template<typename S,Distribution U,Distribution V,typename Ord> 
    friend void LockedView
    ( DistMatrix<S,U,V,Ord>& A, const DistMatrix<S,U,V,Ord>& B );
    template<typename S,Distribution U,Distribution V,typename Ord> 
    friend void View
    ( DistMatrix<S,U,V,Ord>& A, DistMatrix<S,U,V,Ord>& B,
      Ord i, Ord j, Ord height, Ord width );
    template<typename S,Distribution U,Distribution V,typename Ord> 
    friend void LockedView
    ( DistMatrix<S,U,V,Ord>& A, const DistMatrix<S,U,V,Ord>& B,
      Ord i, Ord j, Ord height, Ord width );
    template<typename S,Distribution U,Distribution V,typename Ord> 
    friend void View1x2
    ( DistMatrix<S,U,V,Ord>& A,
      DistMatrix<S,U,V,Ord>& BL, DistMatrix<S,U,V,Ord>& BR );
    template<typename S,Distribution U,Distribution V,typename Ord> 
    friend void LockedView1x2
    (       DistMatrix<S,U,V,Ord>& A,
      const DistMatrix<S,U,V,Ord>& BL,
      const DistMatrix<S,U,V,Ord>& BR );
    template<typename S,Distribution U,Distribution V,typename Ord> 
    friend void View2x1
    ( DistMatrix<S,U,V,Ord>& A,
      DistMatrix<S,U,V,Ord>& BT,
      DistMatrix<S,U,V,Ord>& BB );
    template<typename S,Distribution U,Distribution V,typename Ord> 
    friend void LockedView2x1
    (       DistMatrix<S,U,V,Ord>& A,
      const DistMatrix<S,U,V,Ord>& BT,
      const DistMatrix<S,U,V,Ord>& BB );
    template<typename S,Distribution U,Distribution V,typename Ord> 
    friend void View2x2
    ( DistMatrix<S,U,V,Ord>& A,
      DistMatrix<S,U,V,Ord>& BTL, DistMatrix<S,U,V,Ord>& BTR,
      DistMatrix<S,U,V,Ord>& BBL, DistMatrix<S,U,V,Ord>& BBR );
    template<typename S,Distribution U,Distribution V,typename Ord> 
    friend void LockedView2x2
    (       DistMatrix<S,U,V,Ord>& A,
      const DistMatrix<S,U,V,Ord>& BTL,
      const DistMatrix<S,U,V,Ord>& BTR,
      const DistMatrix<S,U,V,Ord>& BBL,
      const DistMatrix<S,U,V,Ord>& BBR );

    template<typename S,Distribution U,Distribution V,typename Ord>
    friend class DistMatrix;
#endif // ifndef SWIG
};

} // namespace elem

#endif // ifndef CORE_DISTMATRIX_ABSTRACT_DECL_HPP
