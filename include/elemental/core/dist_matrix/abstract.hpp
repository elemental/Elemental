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

template<typename Int>
void AssertConforming1x2
( const DistMatrix_Base<Int>& AL, const DistMatrix_Base<Int>& AR );

template<typename Int>
void AssertConforming2x1
( const DistMatrix_Base<Int>& AT,
  const DistMatrix_Base<Int>& AB );

template<typename Int>
void AssertConforming2x2
( const DistMatrix_Base<Int>& ATL, const DistMatrix_Base<Int>& ATR,
  const DistMatrix_Base<Int>& ABL, const DistMatrix_Base<Int>& ABR );

template<typename Int> 
class DistMatrix_Base
{
protected:
    // Build around a particular grid
    DistMatrix_Base( const elem::Grid& );

public:
    virtual ~DistMatrix_Base();

    //-----------------------------------------------------------------------//
    // Routines that do NOT need to be implemented in derived classes        //
    //-----------------------------------------------------------------------//

    void AssertNotLocked() const;
    void AssertNotStoringData() const;
    void AssertValidEntry( Int i, Int j ) const;
    void AssertValidDimensions( Int height, Int width ) const;
    void AssertValidDimensions( Int height, Int width, Int LDim ) const;
    void AssertValidSubmatrix( Int i, Int j, Int height, Int width ) const;
    void AssertSameGrid( const elem::Grid& grid ) const;
    void AssertSameSize( Int height, Int width ) const;

    //
    // Basic information
    //

    Int Height() const;
    Int Width() const;
    Int DiagonalLength( Int offset=0 ) const;
    virtual Int LocalHeight() const = 0;
    virtual Int LocalWidth() const = 0;
    virtual Int LDim() const = 0;
    virtual size_t DataSize() const = 0;
    virtual size_t AllocatedMemory() const = 0;
    
    const elem::Grid& Grid() const;
    void SetGrid( const elem::Grid& grid );

    void ResizeTo( Int height, Int width );
    void ResizeTo( Int height, Int width, Int ldim );
    
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
    // Viewing 
    //

    bool Viewing() const;
    bool Locked() const;

    //
    // Utilities
    //

    virtual void Empty();
    virtual void EmptyData();
    
    //------------------------------------------------------------------------//
    // Routines that can be overridden in derived classes                     //
    //------------------------------------------------------------------------//

    //
    // Basic information
    //

    virtual Int Root() const;
    virtual Int DiagPath() const;
    virtual bool Participating() const;
    
    virtual void AlignWith( const DistMatrix_Base<Int>& A );
    virtual void AlignColsWith( const DistMatrix_Base<Int>& A );
    virtual void AlignRowsWith( const DistMatrix_Base<Int>& A );
    
    void MakeConsistent();

    //
    // Utilities
    //

    virtual void Attach
    ( Int height, Int width, Int colAlign, Int rowAlign, 
      void* buffer, Int LDim, const elem::Grid& g );
    virtual void LockedAttach
    ( Int height, Int width, Int colAlign, Int rowAlign,
      const void* buffer, Int LDim, const elem::Grid& g );
 
    //------------------------------------------------------------------------//
    // Routines that MUST be implemented in non-abstract derived classes      //
    //------------------------------------------------------------------------//

    //
    // Basic information
    //

    virtual elem::Distribution RowDist() const = 0;
    virtual elem::Distribution ColDist() const = 0;

    // So that the local row indices are given by
    //   A.ColShift():A.ColStride():A.Height()
    virtual Int ColStride() const = 0; 
    // So that the local column indices are given by
    //   A.RowShift():A.RowStride():A.Width()
    virtual Int RowStride() const = 0;
    virtual Int ColRank() const = 0;
    virtual Int RowRank() const = 0;

    virtual bool Index( Int i, Int j, Int& iLocal, Int& jLocal, int& mpiSrc, mpi::Comm& mpiDst ) const = 0;

protected:
    ViewType viewType_;
    Int height_, width_;
    
    bool constrainedColAlignment_, constrainedRowAlignment_;
    Int colAlignment_, rowAlignment_;
    Int colShift_, rowShift_;
    const elem::Grid* grid_;

    virtual void LocalEmpty_() = 0;
    virtual void LocalResize_( Int h, Int w ) = 0;
    virtual void LocalResize_( Int h, Int w, Int ldim ) = 0;
    virtual void LocalAttach_( Int h, Int w, void* buffer, Int ldim ) = 0;
    virtual void LocalLockedAttach_( Int h, Int w, const void* buffer, Int ldim ) = 0;

    void SetColShift();
    void SetRowShift();
    void SetShifts();
    void ResizeTo_( Int height, Int width, Int LDim );

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

    template <Distribution U,Distribution V,typename Ord>
    friend class DistMatrix_Dist;
    template <typename S,typename Ord>
    friend class DistMatrix_Type;
    template<typename S,Distribution U,Distribution V,typename Ord>
    friend class DistMatrix;
#endif // ifndef SWIG
};

template <Distribution U,Distribution V,typename Int>
class DistMatrix_Dist : virtual public DistMatrix_Base<Int>
{
public:
    void AlignWith( const DistMatrix_Base<Int>& A );
    void AlignColsWith( const DistMatrix_Base<Int>& A );
    void AlignRowsWith( const DistMatrix_Base<Int>& A );
    elem::Distribution RowDist() const;
    elem::Distribution ColDist() const;
    Int ColStride() const; 
    Int RowStride() const;
    Int ColRank() const;
    Int RowRank() const;
    Int Root() const;
    bool Participating() const;
protected:
    DistMatrix_Dist( const elem::Grid& g, int root=0 );
    bool Index( Int i, Int j, Int& iLocal, Int& jLocal, Int& mpiSrc, mpi::Comm& comm ) const;
};

template<typename T,typename Int> 
class DistMatrix_Type : virtual public DistMatrix_Base<Int>
{
protected:
    DistMatrix_Type( const elem::Grid& g );
    ~DistMatrix_Type();
    
public:
    Int LocalHeight() const;
    Int LocalWidth() const;
    Int LDim() const;
    size_t AllocatedMemory() const;
    size_t DataSize() const;

          elem::Matrix<T,Int>& Matrix();
    const elem::Matrix<T,Int>& LockedMatrix() const;

          T* Buffer( Int iLocal=0, Int jLocal=0 );
    const T* LockedBuffer( Int iLocal=0, Int jLocal=0 ) const;
    
    //
    // Local entry manipulation
    //

    T GetLocal( Int iLocal, Int jLocal ) const;
    void SetLocal( Int iLocal, Int jLocal, T alpha );
    void UpdateLocal( Int iLocal, Int jLocal, T alpha );

    //
    // Though the following routines are meant for complex data, all but two
    // logically apply to real data.
    //

    BASE(T) GetRealPart( Int i, Int j ) const;
    BASE(T) GetImagPart( Int i, Int j ) const;
    BASE(T) GetLocalRealPart( Int iLocal, Int jLocal ) const;
    BASE(T) GetLocalImagPart( Int iLocal, Int jLocal ) const;
    void SetLocalRealPart( Int iLocal, Int jLocal, BASE(T) alpha );
    void UpdateLocalRealPart( Int iLocal, Int jLocal, BASE(T) alpha );
    // Only valid for complex data
    void SetLocalImagPart( Int iLocal, Int jLocal, BASE(T) alpha );
    void UpdateLocalImagPart( Int iLocal, Int jLocal, BASE(T) alpha );

    //
    // Entry manipulation
    //

    virtual T Get( Int i, Int j ) const;
    virtual void Set( Int i, Int j, T alpha );
    virtual void Update( Int i, Int j, T alpha );

    //
    // Though the following routines are meant for complex data, all but two
    // logically applies to real data.
    //

    virtual void SetRealPart( Int i, Int j, BASE(T) alpha );
    // Only valid for complex data
    virtual void SetImagPart( Int i, Int j, BASE(T) alpha );
    virtual void UpdateRealPart( Int i, Int j, BASE(T) alpha );
    // Only valid for complex data
    virtual void UpdateImagPart( Int i, Int j, BASE(T) alpha );

protected:
    Memory<T> auxMemory_;
    elem::Matrix<T,Int> matrix_;
    
    void LocalEmpty_();
    void LocalResize_( Int height, Int Width );
    void LocalResize_( Int height, Int Width, Int LDim );
    void LocalAttach_( Int height, Int width, void* buffer, Int ldim );
    void LocalLockedAttach_( Int height, Int width, const void* buffer, Int ldim );

    void ComplainIfReal() const;

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
