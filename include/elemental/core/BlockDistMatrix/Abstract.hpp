/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_BLOCKDISTMATRIX_ABSTRACT_DECL_HPP
#define ELEM_BLOCKDISTMATRIX_ABSTRACT_DECL_HPP

namespace elem {

template<typename T> 
class AbstractBlockDistMatrix
{
public:
    // Typedefs
    // ========
    typedef AbstractBlockDistMatrix<T> type;

    // Constructors and destructors
    // ============================
    // Move constructor
    AbstractBlockDistMatrix( type&& A ) ELEM_NOEXCEPT;

    virtual ~AbstractBlockDistMatrix();

    // Assignment and reconfiguration
    // ==============================
    // Move assignment
    type& operator=( type&& A );

    void Empty();
    void EmptyData();
    void SetGrid( const elem::Grid& grid );
    void Resize( Int height, Int width );
    void Resize( Int height, Int width, Int ldim );
    void MakeConsistent( bool includingViewers=false );
    void MakeSizeConsistent( bool includingViewers=false );

    // Realignment
    // -----------
    void Align
    ( Int blockHeight, Int blockWidth, 
      Int colAlign, Int rowAlign, Int colCut=0, Int rowCut=0, 
      bool constrain=true );
    void AlignCols
    ( Int blockHeight, Int colAlign, Int colCut=0, bool constrain=true );
    void AlignRows
    ( Int blockWidth, Int rowAlign, Int rowCut=0, bool constrain=true );
    void FreeAlignments();
    void SetRoot( Int root, bool constrain=true );
    void AlignWith( const elem::BlockDistData& data, bool constrain=true );
    virtual void AlignColsWith
    ( const elem::BlockDistData& data, bool constrain=true );
    virtual void AlignRowsWith
    ( const elem::BlockDistData& data, bool constrain=true );
    // TODO: The interface for these routines could be improved
    void AlignAndResize
    ( Int blockHeight, Int blockWidth, 
      Int colAlign, Int rowAlign, Int colCut, Int rowCut, 
      Int height, Int width, bool force=false, bool constrain=true );
    void AlignColsAndResize
    ( Int blockHeight, Int colAlign, Int colCut, Int height, Int width, 
      bool force=false, bool constrain=true );
    void AlignRowsAndResize
    ( Int blockWidth, Int rowAlign, Int rowCut, Int height, Int width, 
      bool force=false, bool constrain=true );

    // Buffer attachment
    // -----------------
    // (Immutable) view of a distributed matrix's buffer
    void Attach
    ( Int height, Int width, const elem::Grid& g, 
      Int blockHeight, Int blockWidth,
      Int colAlign, Int rowAlign, Int colCut, Int rowCut,
      T* buffer, Int ldim, Int root=0 );
    void LockedAttach
    ( Int height, Int width, const elem::Grid& g,
      Int blockHeight, Int blockWidth,
      Int colAlign, Int rowAlign, Int colCut, Int rowCut,
      const T* buffer, Int ldim, Int root=0 );
    void Attach
    ( Int height, Int width, const elem::Grid& g,
      Int blockHeight, Int blockWidth,
      Int colAlign, Int rowAlign, Int colCut, Int rowCut,
      elem::Matrix<T>& A, Int root=0 );
    void LockedAttach
    ( Int height, Int width, const elem::Grid& g,
      Int blockHeight, Int blockWidth,
      Int colAlign, Int rowAlign, Int colCut, Int rowCut,
      const elem::Matrix<T>& A, Int root=0 );

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
    bool RootConstrained() const;
    Int BlockHeight() const;
    Int BlockWidth() const;
    Int ColAlign() const;
    Int RowAlign() const;
    Int ColCut() const;
    Int RowCut() const;
    Int ColShift() const;
    Int RowShift() const;
    Int ColRank() const;
    Int RowRank() const;
    Int PartialColRank() const;
    Int PartialRowRank() const;
    Int PartialUnionColRank() const;
    Int PartialUnionRowRank() const; 
    Int DistRank() const;
    Int CrossRank() const;
    Int RedundantRank() const;
    Int Root() const;
    bool Participating() const;
    Int RowOwner( Int i ) const;     // rank in ColComm
    Int ColOwner( Int j ) const;     // rank in RowComm
    Int Owner( Int i, Int j ) const; // rank in DistComm
    Int LocalRow( Int i ) const; // debug throws if row i is not locally owned
    Int LocalCol( Int j ) const; // debug throws if col j is not locally owned
    Int LocalRowOffset( Int i ) const; // number of local rows before row i
    Int LocalColOffset( Int j ) const; // number of local cols before col j
    Int GlobalRow( Int iLoc ) const;
    Int GlobalCol( Int jLoc ) const;
    bool IsLocalRow( Int i ) const; 
    bool IsLocalCol( Int j ) const;
    bool IsLocal( Int i, Int j ) const;
    // Must be overridden
    // ^^^^^^^^^^^^^^^^^^
    virtual elem::BlockDistData DistData() const = 0;
    virtual mpi::Comm DistComm() const = 0;
    virtual mpi::Comm CrossComm() const = 0;
    virtual mpi::Comm RedundantComm() const = 0;
    virtual mpi::Comm ColComm() const = 0;
    virtual mpi::Comm RowComm() const = 0;
    virtual mpi::Comm PartialColComm() const;
    virtual mpi::Comm PartialRowComm() const;
    virtual mpi::Comm PartialUnionColComm() const;
    virtual mpi::Comm PartialUnionRowComm() const;
    virtual Int ColStride() const = 0;
    virtual Int RowStride() const = 0;
    virtual Int PartialColStride() const;
    virtual Int PartialRowStride() const;
    virtual Int PartialUnionColStride() const;
    virtual Int PartialUnionRowStride() const;
    virtual Int DistSize() const = 0;
    virtual Int CrossSize() const = 0;
    virtual Int RedundantSize() const = 0;

    // Single-entry manipulation
    // =========================

    // Global entry manipulation 
    // -------------------------
    // NOTE: Local entry manipulation is often much faster and should be
    //       preferred in most circumstances where performance matters.
    T Get( Int i, Int j ) const;
    Base<T> GetRealPart( Int i, Int j ) const;
    Base<T> GetImagPart( Int i, Int j ) const;
    void Set( Int i, Int j, T alpha );
    void SetRealPart( Int i, Int j, Base<T> alpha );
    void SetImagPart( Int i, Int j, Base<T> alpha );
    void Update( Int i, Int j, T alpha );
    void UpdateRealPart( Int i, Int j, Base<T> alpha );
    void UpdateImagPart( Int i, Int j, Base<T> alpha );
    void MakeReal( Int i, Int j );
    void Conjugate( Int i, Int j );

    // Local entry manipulation
    // ------------------------
    T GetLocal( Int iLoc, Int jLoc ) const;
    Base<T> GetLocalRealPart( Int iLoc, Int jLoc ) const;
    Base<T> GetLocalImagPart( Int iLoc, Int jLoc ) const;
    void SetLocal( Int iLoc, Int jLoc, T alpha );
    void SetLocalRealPart( Int iLoc, Int jLoc, Base<T> alpha );
    void SetLocalImagPart( Int iLoc, Int jLoc, Base<T> alpha );
    void UpdateLocal( Int iLoc, Int jLoc, T alpha );
    void UpdateLocalRealPart( Int iLoc, Int jLoc, Base<T> alpha );
    void UpdateLocalImagPart( Int iLoc, Int jLoc, Base<T> alpha );
    void MakeLocalReal( Int iLoc, Int jLoc );
    void ConjugateLocal( Int iLoc, Int jLoc );

    // Diagonal manipulation
    // =====================
    void MakeDiagonalReal( Int offset=0 );
    void ConjugateDiagonal( Int offset=0 );

    // Arbitrary-submatrix manipulation
    // ================================

    // Global submatrix manipulation
    // -----------------------------
    void GetSubmatrix
    ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
      DistMatrix<T,STAR,STAR>& ASub ) const;
    void GetRealPartOfSubmatrix
    ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
      DistMatrix<Base<T>,STAR,STAR>& ASub ) const;
    void GetImagPartOfSubmatrix
    ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
      DistMatrix<Base<T>,STAR,STAR>& ASub ) const;
    DistMatrix<T,STAR,STAR> GetSubmatrix
    ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd ) const;
    DistMatrix<Base<T>,STAR,STAR> GetRealPartOfSubmatrix
    ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd ) const;
    DistMatrix<Base<T>,STAR,STAR> GetImagPartOfSubmatrix
    ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd ) const;

    void SetSubmatrix
    ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
      const DistMatrix<T,STAR,STAR>& ASub );
    void SetRealPartOfSubmatrix
    ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
      const DistMatrix<Base<T>,STAR,STAR>& ASub );
    void SetImagPartOfSubmatrix
    ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
      const DistMatrix<Base<T>,STAR,STAR>& ASub );

    void UpdateSubmatrix
    ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
      T alpha, const DistMatrix<T,STAR,STAR>& ASub );
    void UpdateRealPartOfSubmatrix
    ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
      Base<T> alpha, const DistMatrix<Base<T>,STAR,STAR>& ASub );
    void UpdateImagPartOfSubmatrix
    ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd,
      Base<T> alpha, const DistMatrix<Base<T>,STAR,STAR>& ASub );

    void MakeSubmatrixReal
    ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd );
    void ConjugateSubmatrix
    ( const std::vector<Int>& rowInd, const std::vector<Int>& colInd );

    // Local submatrix manipulation
    // ----------------------------
    void GetLocalSubmatrix
    ( const std::vector<Int>& rowIndLoc, const std::vector<Int>& colIndLoc,
      elem::Matrix<T>& ASub ) const;
    void GetRealPartOfLocalSubmatrix
    ( const std::vector<Int>& rowIndLoc, const std::vector<Int>& colIndLoc,
      elem::Matrix<Base<T>>& ASub ) const;
    void GetImagPartOfLocalSubmatrix
    ( const std::vector<Int>& rowIndLoc, const std::vector<Int>& colIndLoc,
      elem::Matrix<Base<T>>& ASub ) const;
    elem::Matrix<T> GetLocalSubmatrix
    ( const std::vector<Int>& rowIndLoc, 
      const std::vector<Int>& colIndLoc ) const;
    elem::Matrix<Base<T>> GetRealPartOfLocalSubmatrix
    ( const std::vector<Int>& rowIndLoc, 
      const std::vector<Int>& colIndLoc ) const;
    elem::Matrix<Base<T>> GetImagPartOfLocalSubmatrix
    ( const std::vector<Int>& rowIndLoc, 
      const std::vector<Int>& colIndLoc ) const;

    void SetLocalSubmatrix
    ( const std::vector<Int>& rowIndLoc, const std::vector<Int>& colIndLoc,
      const elem::Matrix<T>& ASub );
    void SetRealPartOfLocalSubmatrix
    ( const std::vector<Int>& rowIndLoc, const std::vector<Int>& colIndLoc,
      const elem::Matrix<Base<T>>& ASub );
    void SetImagPartOfLocalSubmatrix
    ( const std::vector<Int>& rowIndLoc, const std::vector<Int>& colIndLoc,
      const elem::Matrix<Base<T>>& ASub );

    void UpdateLocalSubmatrix
    ( const std::vector<Int>& rowIndLoc, const std::vector<Int>& colIndLoc,
      T alpha, const elem::Matrix<T>& ASub );
    void UpdateRealPartOfLocalSubmatrix
    ( const std::vector<Int>& rowIndLoc, const std::vector<Int>& colIndLoc,
      Base<T> alpha, const elem::Matrix<Base<T>>& ASub );
    void UpdateImagPartOfLocalSubmatrix
    ( const std::vector<Int>& rowIndLoc, const std::vector<Int>& colIndLoc,
      Base<T> alpha, const elem::Matrix<Base<T>>& ASub );

    void MakeLocalSubmatrixReal
    ( const std::vector<Int>& rowIndLoc, const std::vector<Int>& colIndLoc );
    void ConjugateLocalSubmatrix
    ( const std::vector<Int>& rowIndLoc, const std::vector<Int>& colIndLoc );

    // Sum over a specified communicator
    // =================================
    void SumOver( mpi::Comm comm );

    // Assertions
    // ==========
    void ComplainIfReal() const;
    void AssertNotLocked() const;
    void AssertNotStoringData() const;
    void AssertValidEntry( Int i, Int j ) const;
    void AssertValidSubmatrix( Int i, Int j, Int height, Int width ) const;
    void AssertSameGrid( const elem::Grid& grid ) const;
    void AssertSameSize( Int height, Int width ) const;

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
    bool colConstrained_, rowConstrained_, rootConstrained_;
    Int blockHeight_, blockWidth_;
    Int colAlign_, rowAlign_,
        colCut_,   rowCut_,
        colShift_, rowShift_;
    Int root_;
    const elem::Grid* grid_;

    // Private constructors
    // ====================
    AbstractBlockDistMatrix( const elem::Grid& g=DefaultGrid(),  Int root=0 );
    AbstractBlockDistMatrix
    ( const elem::Grid& g, Int blockHeight, Int blockWidth, Int root=0 );

    // Exchange metadata with another matrix
    // =====================================
    virtual void ShallowSwap( type& A );

    // Modify the distribution metadata
    // ================================
    void SetShifts();
    void SetColShift();
    void SetRowShift();
    void SetGrid();

    // Friend declarations
    // ===================
    template<typename S,Dist J,Dist K> friend class GeneralDistMatrix;
    template<typename S,Dist J,Dist K> friend class DistMatrix;
    template<typename S,Dist J,Dist K> friend class GeneralBlockDistMatrix;
    template<typename S,Dist J,Dist K> friend class BlockDistMatrix;
};

template<typename T>
void AssertConforming1x2
( const AbstractBlockDistMatrix<T>& AL, const AbstractBlockDistMatrix<T>& AR );

template<typename T>
void AssertConforming2x1
( const AbstractBlockDistMatrix<T>& AT, const AbstractBlockDistMatrix<T>& AB );

template<typename T>
void AssertConforming2x2
( const AbstractBlockDistMatrix<T>& ATL, 
  const AbstractBlockDistMatrix<T>& ATR,
  const AbstractBlockDistMatrix<T>& ABL, 
  const AbstractBlockDistMatrix<T>& ABR );

} // namespace elem

#endif // ifndef ELEM_BLOCKDISTMATRIX_ABSTRACT_DECL_HPP
