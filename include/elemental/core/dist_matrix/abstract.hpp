/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_CORE_DISTMATRIX_ABSTRACT_DECL_HPP
#define ELEM_CORE_DISTMATRIX_ABSTRACT_DECL_HPP

namespace elem {

template<typename T> 
class AbstractDistMatrix
{
public:
    // Typedefs
    // ========
    typedef AbstractDistMatrix<T> type;

    // Constructors and destructors
    // ============================
 #ifndef SWIG
    // Move constructor
    AbstractDistMatrix( type&& A );
#endif
    virtual ~AbstractDistMatrix();

    // Assignment and reconfiguration
    // ==============================
#ifndef SWIG
    // Move assignment
    type& operator=( type&& A );
#endif
    void Empty();
    void EmptyData();
    void SetGrid( const elem::Grid& grid );
    void Resize( Int height, Int width );
    void Resize( Int height, Int width, Int ldim );
    void MakeConsistent();
    // Realignment
    // -----------
    void Align( Int colAlign, Int rowAlign );
    void AlignCols( Int colAlign );
    void AlignRows( Int rowAlign );
    void FreeAlignments();
    void SetRoot( Int root );
    virtual void AlignWith( const elem::DistData& data );
    virtual void AlignColsWith( const elem::DistData& data );
    virtual void AlignRowsWith( const elem::DistData& data );

    // Basic queries
    // =============

    // Global matrix information
    // -------------------------
    Int Height() const;
    Int Width() const;
    Int DiagonalLength( Int offset=0 ) const;
    bool Viewing() const;
    bool Locked()  const;

    // Local matrix information
    // ------------------------
    Int LocalHeight() const;
    Int LocalWidth() const;
    Int LDim() const;
          elem::Matrix<T>& Matrix();
    const elem::Matrix<T>& LockedMatrix() const;
    size_t AllocatedMemory() const;
          T* Buffer();
          T* Buffer( Int iLoc, Int jLoc );
    const T* LockedBuffer() const;
    const T* LockedBuffer( Int iLoc, Int jLoc ) const;

    // Distribution information
    // ------------------------
    const elem::Grid& Grid() const;
    bool ColConstrained() const;
    bool RowConstrained() const;
    Int ColAlign() const;
    Int RowAlign() const;
    Int ColShift() const;
    Int RowShift() const;
    Int ColRank() const;
    Int RowRank() const;
    Int DistRank() const;
    Int CrossRank() const;
    Int RedundantRank() const;
    Int DistSize() const;
    Int CrossSize() const;
    Int RedundantSize() const;
    Int Root() const;
    bool Participating() const;
    Int RowOwner( Int i ) const;     // rank in ColComm
    Int ColOwner( Int j ) const;     // rank in RowComm
    Int Owner( Int i, Int j ) const; // rank in DistComm
    Int LocalRow( Int i ) const; // debug throws if row i is not locally owned
    Int LocalCol( Int j ) const; // debug throws if col j is not locally owned
    bool IsLocalRow( Int i ) const; 
    bool IsLocalCol( Int j ) const;
    bool IsLocal( Int i, Int j ) const;
    // Must be overridden
    // ^^^^^^^^^^^^^^^^^^
    virtual elem::DistData DistData() const = 0;
    virtual mpi::Comm DistComm() const = 0;
    virtual mpi::Comm CrossComm() const = 0;
    virtual mpi::Comm RedundantComm() const = 0;
    virtual mpi::Comm ColComm() const = 0;
    virtual mpi::Comm RowComm() const = 0;
    virtual Int ColStride() const = 0;
    virtual Int RowStride() const = 0;

    // Single-entry manipulation
    // =========================

    // Global entry manipulation 
    // -------------------------
    // NOTE: Local entry manipulation is often much faster and should be
    //       preferred in most circumstances where performance matters.
    T Get( Int i, Int j ) const;
    BASE(T) GetRealPart( Int i, Int j ) const;
    BASE(T) GetImagPart( Int i, Int j ) const;
    void Set( Int i, Int j, T alpha );
    void SetRealPart( Int i, Int j, BASE(T) alpha );
    void SetImagPart( Int i, Int j, BASE(T) alpha );
    void Update( Int i, Int j, T alpha );
    void UpdateRealPart( Int i, Int j, BASE(T) alpha );
    void UpdateImagPart( Int i, Int j, BASE(T) alpha );
    void MakeReal( Int i, Int j );
    void Conjugate( Int i, Int j );

    // Local entry manipulation
    // ------------------------
    T GetLocal( Int iLoc, Int jLoc ) const;
    BASE(T) GetLocalRealPart( Int iLoc, Int jLoc ) const;
    BASE(T) GetLocalImagPart( Int iLoc, Int jLoc ) const;
    void SetLocal( Int iLoc, Int jLoc, T alpha );
    void SetLocalRealPart( Int iLoc, Int jLoc, BASE(T) alpha );
    void SetLocalImagPart( Int iLoc, Int jLoc, BASE(T) alpha );
    void UpdateLocal( Int iLoc, Int jLoc, T alpha );
    void UpdateLocalRealPart( Int iLoc, Int jLoc, BASE(T) alpha );
    void UpdateLocalImagPart( Int iLoc, Int jLoc, BASE(T) alpha );
    void MakeRealLocal( Int i, Int j );
    void ConjugateLocal( Int i, Int j );

    // Diagonal manipulation
    // =====================
    // TODO

    // Arbitrary submatrix manipulation
    // ================================
    // TODO

protected:
    // Member variables
    // ================

    // Global and local matrix information 
    // -----------------------------------
    ViewType viewType_;
    Int height_, width_;
    Memory<T> auxMemory_;
    elem::Matrix<T> matrix_;
    
    // Process grid and distribution metadata
    // --------------------------------------
    bool colConstrained_, rowConstrained_;
    Int colAlign_, rowAlign_,
        colShift_, rowShift_;
    Int root_;
    const elem::Grid* grid_;

    // Construct using a particular process grid
    // =========================================
    AbstractDistMatrix( const elem::Grid& g );

    // Exchange metadata with another matrix
    // =====================================
    virtual void ShallowSwap( type& A );

    // Modify the distribution metadata
    // ================================
    void SetShifts();
    void SetColShift();
    void SetRowShift();
    void SetGrid();

    // Combined realignment and resize
    // ===============================
    // NOTE: Should these be public?
    void AlignAndResize
    ( Int colAlign, Int rowAlign, Int height, Int width, bool force=false );
    void AlignColsAndResize
    ( Int colAlign, Int height, Int width, bool force=false );
    void AlignRowsAndResize
    ( Int rowAlign, Int height, Int width, bool force=false );

    // Assertions
    // ==========
    void ComplainIfReal() const;
    void AssertNotLocked() const;
    void AssertNotStoringData() const;
    void AssertValidEntry( Int i, Int j ) const;
    void AssertValidSubmatrix( Int i, Int j, Int height, Int width ) const;
    void AssertSameGrid( const elem::Grid& grid ) const;
    void AssertSameSize( Int height, Int width ) const;

    // Friend declarations
    // ===================
#ifndef SWIG
    template<typename S,Dist U,Dist V> friend class DistMatrix;

    template<typename S,Dist U,Dist V> 
    friend void View( DistMatrix<S,U,V>& A, DistMatrix<S,U,V>& B );
    template<typename S,Dist U,Dist V> 
    friend void LockedView( DistMatrix<S,U,V>& A, const DistMatrix<S,U,V>& B );

    template<typename S,Dist U,Dist V> 
    friend void View
    ( DistMatrix<S,U,V>& A, DistMatrix<S,U,V>& B,
      Int i, Int j, Int height, Int width );
    template<typename S,Dist U,Dist V> 
    friend void LockedView
    ( DistMatrix<S,U,V>& A, const DistMatrix<S,U,V>& B,
      Int i, Int j, Int height, Int width );

    template<typename S,Dist U,Dist V> 
    friend void View1x2
    ( DistMatrix<S,U,V>& A, DistMatrix<S,U,V>& BL, DistMatrix<S,U,V>& BR );
    template<typename S,Dist U,Dist V> 
    friend void LockedView1x2
    (       DistMatrix<S,U,V>& A,
      const DistMatrix<S,U,V>& BL, const DistMatrix<S,U,V>& BR );

    template<typename S,Dist U,Dist V> 
    friend void View2x1
    ( DistMatrix<S,U,V>& A, DistMatrix<S,U,V>& BT, DistMatrix<S,U,V>& BB );
    template<typename S,Dist U,Dist V> 
    friend void LockedView2x1
    (       DistMatrix<S,U,V>& A,
      const DistMatrix<S,U,V>& BT, const DistMatrix<S,U,V>& BB );

    template<typename S,Dist U,Dist V> 
    friend void View2x2
    ( DistMatrix<S,U,V>& A,
      DistMatrix<S,U,V>& BTL, DistMatrix<S,U,V>& BTR,
      DistMatrix<S,U,V>& BBL, DistMatrix<S,U,V>& BBR );
    template<typename S,Dist U,Dist V> 
    friend void LockedView2x2
    (       DistMatrix<S,U,V>& A,
      const DistMatrix<S,U,V>& BTL, const DistMatrix<S,U,V>& BTR,
      const DistMatrix<S,U,V>& BBL, const DistMatrix<S,U,V>& BBR );
#endif // ifndef SWIG
};

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

} // namespace elem

#endif // ifndef ELEM_CORE_DISTMATRIX_ABSTRACT_DECL_HPP
