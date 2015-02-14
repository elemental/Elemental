/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_DISTMATRIX_ABSTRACT_DECL_HPP
#define EL_DISTMATRIX_ABSTRACT_DECL_HPP

namespace El {

template<typename T> 
class AbstractDistMatrix
{
public:
    // Typedefs
    // ========
    typedef AbstractDistMatrix<T> type;

    // Constructors and destructors
    // ============================
    // Move constructor
    AbstractDistMatrix( type&& A ) EL_NOEXCEPT;

    virtual ~AbstractDistMatrix();

    virtual AbstractDistMatrix<T>* Construct
    ( const El::Grid& g, int root ) const = 0;
    virtual AbstractDistMatrix<T>* ConstructTranspose
    ( const El::Grid& g, int root ) const = 0;
    virtual AbstractDistMatrix<T>* ConstructDiagonal
    ( const El::Grid& g, int root ) const = 0;
    // TODO: ConstructPartialCol and friends?

    // Assignment and reconfiguration
    // ==============================
    const type& operator=( const DistMultiVec<T>& A );
    // Move assignment
    type& operator=( type&& A );

    void Empty();
    void EmptyData();
    void SetGrid( const El::Grid& grid );
    void Resize( Int height, Int width );
    void Resize( Int height, Int width, Int ldim );
    void MakeConsistent( bool includingViewers=false );
    void MakeSizeConsistent( bool includingViewers=false );

    // Realignment
    // -----------
    void Align( int colAlign, int rowAlign, bool constrain=true );
    void AlignCols( int colAlign, bool constrain=true );
    void AlignRows( int rowAlign, bool constrain=true );
    void FreeAlignments();
    void SetRoot( int root, bool constrain=true );
    void AlignWith
    ( const El::DistData& data, bool constrain=true, bool allowMismatch=false );
    void AlignColsWith
    ( const El::DistData& data, bool constrain=true, bool allowMismatch=false );
    void AlignRowsWith
    ( const El::DistData& data, bool constrain=true, bool allowMismatch=false );
    void AlignAndResize
    ( int colAlign, int rowAlign, Int height, Int width, 
      bool force=false, bool constrain=true );
    void AlignColsAndResize
    ( int colAlign, Int height, Int width, 
      bool force=false, bool constrain=true );
    void AlignRowsAndResize
    ( int rowAlign, Int height, Int width, 
      bool force=false, bool constrain=true );

    // Buffer attachment
    // -----------------
    // (Immutable) view of a distributed matrix's buffer
    void Attach
    ( Int height, Int width, const El::Grid& grid, 
      int colAlign, int rowAlign, T* buffer, Int ldim, int root=0 );
    void LockedAttach
    ( Int height, Int width, const El::Grid& grid,
      int colAlign, int rowAlign, const T* buffer, Int ldim, int root=0 );
    void Attach
    ( Int height, Int width, const El::Grid& grid,
      int colAlign, int rowAlign, El::Matrix<T>& A, int root=0 );
    void LockedAttach
    ( Int height, Int width, const El::Grid& grid,
      int colAlign, int rowAlign, const El::Matrix<T>& A, int root=0 );
    // (Immutable) view of a local matrix's buffer
    void Attach( const El::Grid& grid, El::Matrix<T>& A );
    void LockedAttach( const El::Grid& grid, const El::Matrix<T>& A );

    // Basic queries
    // =============

    // Global matrix information
    // -------------------------
    Int  Height()                       const;
    Int  Width()                        const;
    Int  DiagonalLength( Int offset=0 ) const;
    bool Viewing()                      const;
    bool Locked()                       const;

    // Local matrix information
    // ------------------------
          Int            LocalHeight()                      const;
          Int            LocalWidth()                       const;
          Int            LDim()                             const;
          El::Matrix<T>& Matrix();
    const El::Matrix<T>& LockedMatrix()                     const;
          size_t         AllocatedMemory()                  const;
          T*             Buffer();
          T*             Buffer( Int iLoc, Int jLoc );
    const T*             LockedBuffer()                     const;
    const T*             LockedBuffer( Int iLoc, Int jLoc ) const;

    // Distribution information
    // ------------------------
            const El::Grid&    Grid()                  const;
                  int          ColAlign()              const;
                  int          RowAlign()              const;
                  int          ColShift()              const;
                  int          RowShift()              const;
                  bool         ColConstrained()        const;
                  bool         RowConstrained()        const;
                  bool         RootConstrained()       const;
                  bool         Participating()         const;

                  int          RowOwner( Int i )       const;     
                  int          ColOwner( Int j )       const;     
                  int          Owner( Int i, Int j )   const; 
                  Int          LocalRow( Int i )       const;     
                  Int          LocalCol( Int j )       const; 
                  Int          LocalRowOffset( Int i ) const; 
                  Int          LocalColOffset( Int j ) const; 
                  Int          GlobalRow( Int iLoc )   const;
                  Int          GlobalCol( Int jLoc )   const;
                  bool         IsLocalRow( Int i )     const; 
                  bool         IsLocalCol( Int j )     const;
                  bool         IsLocal( Int i, Int j ) const;

    virtual       Dist         ColDist()               const = 0;
    virtual       Dist         RowDist()               const = 0;
    virtual       Dist         CollectedColDist()      const = 0;
    virtual       Dist         CollectedRowDist()      const = 0;
    virtual       Dist         PartialColDist()        const = 0;
    virtual       Dist         PartialRowDist()        const = 0;
    virtual       Dist         PartialUnionColDist()   const = 0;
    virtual       Dist         PartialUnionRowDist()   const = 0;

    virtual       mpi::Comm    ColComm()               const = 0;
    virtual       mpi::Comm    RowComm()               const = 0;
    virtual       mpi::Comm    PartialColComm()        const;
    virtual       mpi::Comm    PartialRowComm()        const;
    virtual       mpi::Comm    PartialUnionColComm()   const;
    virtual       mpi::Comm    PartialUnionRowComm()   const;
    virtual       mpi::Comm    DistComm()              const = 0;
    virtual       mpi::Comm    CrossComm()             const = 0;
    virtual       mpi::Comm    RedundantComm()         const = 0;

    // NOTE: While these would be equivalent to composing mpi::Size
    //       with ColComm(), RowComm(), etc., for processes *within*
    //       the communicator, this composition is incorrect for 
    //       processes *outside* of the communicator
    virtual       int          ColStride()             const = 0;
    virtual       int          RowStride()             const = 0;
    virtual       int          PartialColStride()      const;
    virtual       int          PartialRowStride()      const;
    virtual       int          PartialUnionColStride() const;
    virtual       int          PartialUnionRowStride() const;
    virtual       int          DistSize()              const = 0;
    virtual       int          CrossSize()             const = 0;
    virtual       int          RedundantSize()         const = 0;

    // NOTE: These are all clearly equivalent to composing mpi::Rank
    //       with ColComm(), RowComm(), etc., but it is not clear that
    //       they should be removed just yet.
                  int          ColRank()               const;
                  int          RowRank()               const;
                  int          PartialColRank()        const;
                  int          PartialRowRank()        const;
                  int          PartialUnionColRank()   const;
                  int          PartialUnionRowRank()   const; 
                  int          DistRank()              const;
                  int          CrossRank()             const;
                  int          RedundantRank()         const;

                  int          Root()                  const;
    virtual       El::DistData DistData()              const = 0;

    // Single-entry manipulation
    // =========================

    // Global entry manipulation 
    // -------------------------
    // NOTE: Local entry manipulation is often much faster and should be
    //       preferred in most circumstances where performance matters.
    T       Get( Int i, Int j )                            const;
    Base<T> GetRealPart( Int i, Int j )                    const;
    Base<T> GetImagPart( Int i, Int j )                    const;
    void    Set( Int i, Int j, T alpha );
    void    SetRealPart( Int i, Int j, Base<T> alpha );
    void    SetImagPart( Int i, Int j, Base<T> alpha );
    void    Update( Int i, Int j, T alpha );
    void    UpdateRealPart( Int i, Int j, Base<T> alpha );
    void    UpdateImagPart( Int i, Int j, Base<T> alpha );
    void    MakeReal( Int i, Int j );
    void    Conjugate( Int i, Int j );

    // Local entry manipulation
    // ------------------------
    // NOTE: Clearly each of the following routines could instead be performed
    //       via composing [Locked]Matrix() with the corresponding local
    //       routine, but a large amount of code might need to change if 
    //       these were removed.
    T       GetLocal( Int iLoc, Int jLoc )                            const;
    Base<T> GetLocalRealPart( Int iLoc, Int jLoc )                    const;
    Base<T> GetLocalImagPart( Int iLoc, Int jLoc )                    const;
    void    SetLocal( Int iLoc, Int jLoc, T alpha );
    void    SetLocalRealPart( Int iLoc, Int jLoc, Base<T> alpha );
    void    SetLocalImagPart( Int iLoc, Int jLoc, Base<T> alpha );
    void    UpdateLocal( Int iLoc, Int jLoc, T alpha );
    void    UpdateLocalRealPart( Int iLoc, Int jLoc, Base<T> alpha );
    void    UpdateLocalImagPart( Int iLoc, Int jLoc, Base<T> alpha );
    void    MakeLocalReal( Int iLoc, Int jLoc );
    void    ConjugateLocal( Int iLoc, Int jLoc );

    // Diagonal manipulation
    // =====================
    bool DiagonalAlignedWith( const El::DistData& d, Int offset=0 ) const;
    int DiagonalRoot( Int offset=0 ) const;
    int DiagonalAlign( Int offset=0 ) const;

    // Assertions
    // ==========
    void ComplainIfReal() const;
    void AssertNotLocked() const;
    void AssertNotStoringData() const;
    void AssertValidEntry( Int i, Int j ) const;
    void AssertValidSubmatrix( Int i, Int j, Int height, Int width ) const;
    void AssertSameSize( Int height, Int width ) const;

protected:
    // Member variables
    // ================

    // Global and local matrix information 
    // -----------------------------------
    ViewType viewType_;
    Int height_, width_;
    El::Matrix<T> matrix_;
    
    // Process grid and distribution metadata
    // --------------------------------------
    bool colConstrained_, rowConstrained_, rootConstrained_;
    int colAlign_, rowAlign_,
        colShift_, rowShift_;
    int root_;
    const El::Grid* grid_;

    // Private constructors
    // ====================
    // Create a 0 x 0 distributed matrix
    AbstractDistMatrix( const El::Grid& g=DefaultGrid(), int root=0 );

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
    template<typename S,Dist J,Dist K> friend class DistMatrix;

    template<typename S,Dist J,Dist K> friend class GeneralBlockDistMatrix;
    template<typename S,Dist J,Dist K> friend class BlockDistMatrix;
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

} // namespace El

#endif // ifndef EL_DISTMATRIX_ABSTRACT_DECL_HPP
