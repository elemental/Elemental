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
    AbstractDistMatrix( type&& A ) EL_NO_EXCEPT;

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
    void Empty();
    void EmptyData();
    void SetGrid( const El::Grid& grid );
    void Resize( Int height, Int width );
    void Resize( Int height, Int width, Int ldim );
    // NOTE: The 'consistent' in this name should not be confused with the
    //       'AssertLocallyConsistent()' function, which has to do with whether
    //       the local queues have been processed. Perhaps these routine names
    //       should be changed.
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

    // Operator overloading
    // ====================

    // Copy
    // ----
    const type& operator=( const type& A );
    const type& operator=( const DistMultiVec<T>& A );

    // Move assignment
    // ---------------
    type& operator=( type&& A );

    // Rescaling
    // ---------
    const type& operator*=( T alpha );

    // Addition/subtraction
    // --------------------
    const type& operator+=( const type& A );
    const type& operator-=( const type& A );

    // Basic queries
    // =============

    // Global matrix information
    // -------------------------
    Int  Height()                       const EL_NO_EXCEPT;
    Int  Width()                        const EL_NO_EXCEPT;
    Int  DiagonalLength( Int offset=0 ) const EL_NO_EXCEPT;
    bool Viewing()                      const EL_NO_EXCEPT;
    bool Locked()                       const EL_NO_EXCEPT;

    // Local matrix information
    // ------------------------
          Int            LocalHeight()                      const EL_NO_EXCEPT;
          Int            LocalWidth()                       const EL_NO_EXCEPT;
          Int            LDim()                             const EL_NO_EXCEPT;
          El::Matrix<T>& Matrix()                                 EL_NO_EXCEPT;
    const El::Matrix<T>& LockedMatrix()                     const EL_NO_EXCEPT;
          size_t         AllocatedMemory()                  const EL_NO_EXCEPT;
          T*             Buffer() EL_NO_RELEASE_EXCEPT;
          T*             Buffer( Int iLoc, Int jLoc ) EL_NO_RELEASE_EXCEPT;
    const T*             LockedBuffer()                     const EL_NO_EXCEPT;
    const T*             LockedBuffer( Int iLoc, Int jLoc ) const EL_NO_EXCEPT;

    // Distribution information
    // ------------------------
            const El::Grid&    Grid()                  const EL_NO_EXCEPT;
                  int          ColAlign()              const EL_NO_EXCEPT;
                  int          RowAlign()              const EL_NO_EXCEPT;
                  int          ColShift()              const EL_NO_EXCEPT;
                  int          RowShift()              const EL_NO_EXCEPT;
                  bool         ColConstrained()        const EL_NO_EXCEPT;
                  bool         RowConstrained()        const EL_NO_EXCEPT;
                  bool         RootConstrained()       const EL_NO_EXCEPT;
                  bool         Participating()         const EL_NO_RELEASE_EXCEPT;

                  int          RowOwner( Int i )       const EL_NO_EXCEPT;
                  int          ColOwner( Int j )       const EL_NO_EXCEPT;
                  int          Owner( Int i, Int j )   const EL_NO_EXCEPT;
                  Int          LocalRow( Int i )       const EL_NO_RELEASE_EXCEPT;
                  Int          LocalCol( Int j )       const EL_NO_RELEASE_EXCEPT;
                  Int          LocalRowOffset( Int i ) const EL_NO_EXCEPT;
                  Int          LocalColOffset( Int j ) const EL_NO_EXCEPT;
                  Int          GlobalRow( Int iLoc )   const EL_NO_EXCEPT;
                  Int          GlobalCol( Int jLoc )   const EL_NO_EXCEPT;
                  bool         IsLocalRow( Int i )     const EL_NO_RELEASE_EXCEPT;
                  bool         IsLocalCol( Int j )     const EL_NO_RELEASE_EXCEPT;
                  bool         IsLocal( Int i, Int j ) const EL_NO_RELEASE_EXCEPT;

    virtual       Dist         ColDist()               const EL_NO_EXCEPT = 0;
    virtual       Dist         RowDist()               const EL_NO_EXCEPT = 0;
    virtual       Dist         CollectedColDist()      const EL_NO_EXCEPT = 0;
    virtual       Dist         CollectedRowDist()      const EL_NO_EXCEPT = 0;
    virtual       Dist         PartialColDist()        const EL_NO_EXCEPT = 0;
    virtual       Dist         PartialRowDist()        const EL_NO_EXCEPT = 0;
    virtual       Dist         PartialUnionColDist()   const EL_NO_EXCEPT = 0;
    virtual       Dist         PartialUnionRowDist()   const EL_NO_EXCEPT = 0;

    virtual       mpi::Comm    ColComm()               const EL_NO_EXCEPT = 0;
    virtual       mpi::Comm    RowComm()               const EL_NO_EXCEPT = 0;
    virtual       mpi::Comm    PartialColComm()        const EL_NO_EXCEPT;
    virtual       mpi::Comm    PartialRowComm()        const EL_NO_EXCEPT;
    virtual       mpi::Comm    PartialUnionColComm()   const EL_NO_EXCEPT;
    virtual       mpi::Comm    PartialUnionRowComm()   const EL_NO_EXCEPT;
    virtual       mpi::Comm    DistComm()              const EL_NO_EXCEPT = 0;
    virtual       mpi::Comm    CrossComm()             const EL_NO_EXCEPT = 0;
    virtual       mpi::Comm    RedundantComm()         const EL_NO_EXCEPT = 0;

    // NOTE: While these would be equivalent to composing mpi::Size
    //       with ColComm(), RowComm(), etc., for processes *within*
    //       the communicator, this composition is incorrect for 
    //       processes *outside* of the communicator
    virtual       int          ColStride()             const EL_NO_EXCEPT= 0;
    virtual       int          RowStride()             const EL_NO_EXCEPT= 0;
    virtual       int          PartialColStride()      const EL_NO_EXCEPT;
    virtual       int          PartialRowStride()      const EL_NO_EXCEPT;
    virtual       int          PartialUnionColStride() const EL_NO_EXCEPT;
    virtual       int          PartialUnionRowStride() const EL_NO_EXCEPT;
    virtual       int          DistSize()              const EL_NO_EXCEPT= 0;
    virtual       int          CrossSize()             const EL_NO_EXCEPT= 0;
    virtual       int          RedundantSize()         const EL_NO_EXCEPT= 0;

    // NOTE: These are all clearly equivalent to composing mpi::Rank
    //       with ColComm(), RowComm(), etc., but it is not clear that
    //       they should be removed just yet.
                  int          ColRank()               const EL_NO_RELEASE_EXCEPT;
                  int          RowRank()               const EL_NO_RELEASE_EXCEPT;
                  int          PartialColRank()        const EL_NO_RELEASE_EXCEPT;
                  int          PartialRowRank()        const EL_NO_RELEASE_EXCEPT;
                  int          PartialUnionColRank()   const EL_NO_RELEASE_EXCEPT;
                  int          PartialUnionRowRank()   const EL_NO_RELEASE_EXCEPT;
                  int          DistRank()              const EL_NO_RELEASE_EXCEPT;
                  int          CrossRank()             const EL_NO_RELEASE_EXCEPT;
                  int          RedundantRank()         const EL_NO_RELEASE_EXCEPT;

                  int          Root()                  const EL_NO_EXCEPT;
    virtual       El::DistData DistData()              const = 0;

    // Single-entry manipulation
    // =========================

    // Global entry manipulation 
    // -------------------------
    // NOTE: Local entry manipulation is often much faster and should be
    //       preferred in most circumstances where performance matters.
    T       Get( Int i, Int j ) const EL_NO_RELEASE_EXCEPT;
    Base<T> GetRealPart( Int i, Int j ) const EL_NO_RELEASE_EXCEPT;
    Base<T> GetImagPart( Int i, Int j ) const EL_NO_RELEASE_EXCEPT;
    void    Set( Int i, Int j, T alpha ) EL_NO_RELEASE_EXCEPT;
    void    Set( const Entry<T>& entry ) EL_NO_RELEASE_EXCEPT;
    void    SetRealPart( Int i, Int j, Base<T> alpha ) EL_NO_RELEASE_EXCEPT;
    void    SetImagPart( Int i, Int j, Base<T> alpha ) EL_NO_RELEASE_EXCEPT;
    void    SetRealPart( const Entry<Base<T>>& entry ) EL_NO_RELEASE_EXCEPT;
    void    SetImagPart( const Entry<Base<T>>& entry ) EL_NO_RELEASE_EXCEPT;
    void    Update( Int i, Int j, T alpha ) EL_NO_RELEASE_EXCEPT;
    void    Update( const Entry<T>& entry ) EL_NO_RELEASE_EXCEPT;
    void    UpdateRealPart( Int i, Int j, Base<T> alpha ) EL_NO_RELEASE_EXCEPT;
    void    UpdateImagPart( Int i, Int j, Base<T> alpha ) EL_NO_RELEASE_EXCEPT;
    void    UpdateRealPart( const Entry<Base<T>>& entry ) EL_NO_RELEASE_EXCEPT;
    void    UpdateImagPart( const Entry<Base<T>>& entry ) EL_NO_RELEASE_EXCEPT;
    void    MakeReal( Int i, Int j ) EL_NO_RELEASE_EXCEPT;
    void    Conjugate( Int i, Int j ) EL_NO_RELEASE_EXCEPT;

    // Batch updating of remote entries
    // ---------------------------------
    void Reserve( Int numRemoteEntries );
    void QueueUpdate( const Entry<T>& entry ) EL_NO_RELEASE_EXCEPT;
    void QueueUpdate( Int i, Int j, T value ) EL_NO_RELEASE_EXCEPT;
    void ProcessQueues();

    // Local entry manipulation
    // ------------------------
    // NOTE: Clearly each of the following routines could instead be performed
    //       via composing [Locked]Matrix() with the corresponding local
    //       routine, but a large amount of code might need to change if 
    //       these were removed.
    T       GetLocal( Int iLoc, Int jLoc ) const EL_NO_RELEASE_EXCEPT;
    Base<T> GetLocalRealPart( Int iLoc, Int jLoc ) const EL_NO_RELEASE_EXCEPT;
    Base<T> GetLocalImagPart( Int iLoc, Int jLoc ) const EL_NO_RELEASE_EXCEPT;
    void    SetLocal( Int iLoc, Int jLoc, T alpha ) EL_NO_RELEASE_EXCEPT;
    void    SetLocal( const Entry<T>& localEntry ) EL_NO_RELEASE_EXCEPT;
    void    SetLocalRealPart( Int iLoc, Int jLoc, Base<T> alpha ) EL_NO_RELEASE_EXCEPT;
    void    SetLocalImagPart( Int iLoc, Int jLoc, Base<T> alpha ) EL_NO_RELEASE_EXCEPT;
    void    SetLocalRealPart( const Entry<Base<T>>& localEntry ) EL_NO_RELEASE_EXCEPT;
    void    SetLocalImagPart( const Entry<Base<T>>& localEntry ) EL_NO_RELEASE_EXCEPT;
    void    UpdateLocal( Int iLoc, Int jLoc, T alpha ) EL_NO_RELEASE_EXCEPT;
    void    UpdateLocal( const Entry<T>& localEntry ) EL_NO_RELEASE_EXCEPT;
    void    UpdateLocalRealPart( Int iLoc, Int jLoc, Base<T> alpha ) EL_NO_RELEASE_EXCEPT;
    void    UpdateLocalImagPart( Int iLoc, Int jLoc, Base<T> alpha ) EL_NO_RELEASE_EXCEPT;
    void    UpdateLocalRealPart( const Entry<Base<T>>& localEntry ) EL_NO_RELEASE_EXCEPT;
    void    UpdateLocalImagPart( const Entry<Base<T>>& localEntry ) EL_NO_RELEASE_EXCEPT;
    void    MakeLocalReal( Int iLoc, Int jLoc ) EL_NO_RELEASE_EXCEPT;
    void    ConjugateLocal( Int iLoc, Int jLoc ) EL_NO_RELEASE_EXCEPT;

    // Diagonal manipulation
    // =====================
    bool DiagonalAlignedWith( const El::DistData& d, Int offset=0 ) const EL_NO_EXCEPT;
    int DiagonalRoot( Int offset=0 ) const EL_NO_EXCEPT;
    int DiagonalAlign( Int offset=0 ) const EL_NO_EXCEPT;

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

    // Remote updates
    // --------------
    vector<Entry<T>> remoteUpdates_;

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
    template<typename S,Dist J,Dist K,DistWrap wrap> friend class DistMatrix;

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
