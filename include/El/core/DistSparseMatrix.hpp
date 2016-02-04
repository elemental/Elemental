/*
   Copyright (c) 2009-2016, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, Stanford University, and the
   Georgia Insitute of Technology.
   All rights reserved.
 
   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_CORE_DISTSPARSEMATRIX_DECL_HPP
#define EL_CORE_DISTSPARSEMATRIX_DECL_HPP

namespace El {

struct DistSparseMultMeta
{
    bool ready;
    // NOTE: The 'send' and 'recv' roles reverse for adjoint multiplication
    Int numRecvInds;
    vector<int> sendSizes, sendOffs,
                recvSizes, recvOffs;
    vector<Int> sendInds, colOffs;

    DistSparseMultMeta() : ready(false), numRecvInds(0) { }

    void Clear()
    {
        ready = false;
        numRecvInds = 0; 
        SwapClear( sendSizes );
        SwapClear( recvSizes );
        SwapClear( sendOffs );
        SwapClear( recvOffs );
        SwapClear( sendInds );
        SwapClear( colOffs );
    }

    const DistSparseMultMeta& operator=( const DistSparseMultMeta& meta )
    {
        ready = meta.ready;
        numRecvInds = meta.numRecvInds;
        sendSizes = meta.sendSizes;
        sendOffs = meta.sendOffs;
        recvSizes = meta.recvSizes;
        recvOffs = meta.recvOffs;
        sendInds = meta.sendInds;
        colOffs = meta.colOffs;
        return *this;
    }
};

// Use a simple 1d distribution where each process owns a fixed number of rows,
//     if last process,  height - (commSize-1)*floor(height/commSize)
//     otherwise,        floor(height/commSize)
template<typename T>
class DistSparseMatrix
{
public:
    // Constructors and destructors
    // ============================
    DistSparseMatrix( mpi::Comm comm=mpi::COMM_WORLD );
    DistSparseMatrix( Int height, Int width, mpi::Comm comm=mpi::COMM_WORLD );
    DistSparseMatrix( const DistSparseMatrix<T>& A );
    // TODO: Move constructor
    ~DistSparseMatrix();

    // Assignment and reconfiguration
    // ==============================

    // Change the size of the matrix
    // -----------------------------
    void Empty( bool freeMemory=true );
    void Resize( Int height, Int width );

    // Change the distribution
    // -----------------------
    void SetComm( mpi::Comm comm );

    // Assembly
    // --------
    void Reserve( Int numLocalEntries, Int numRemoteEntries=0 );

    void FreezeSparsity() EL_NO_EXCEPT;
    void UnfreezeSparsity() EL_NO_EXCEPT;
    bool FrozenSparsity() const EL_NO_EXCEPT;

    // Expensive independent updates and explicit zeroing
    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    void Update( const Entry<T>& entry );
    void Update( Int row, Int col, T value );
    void Zero( Int row, Int col );

    void UpdateLocal( const Entry<T>& localEntry );
    void UpdateLocal( Int localRow, Int col, T value );
    void ZeroLocal( Int localRow, Int col );

    // Batch updating and zeroing (recommended)
    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    void QueueUpdate( const Entry<T>& entry, bool passive=false )
    EL_NO_RELEASE_EXCEPT;
    void QueueUpdate( Int row, Int col, T value, bool passive=false )
    EL_NO_RELEASE_EXCEPT;
    void QueueZero( Int row, Int col, bool passive=false )
    EL_NO_RELEASE_EXCEPT;

    void QueueLocalUpdate( const Entry<T>& localEntry )
    EL_NO_RELEASE_EXCEPT;
    void QueueLocalUpdate( Int localRow, Int col, T value )
    EL_NO_RELEASE_EXCEPT;
    void QueueLocalZero( Int localRow, Int col )
    EL_NO_RELEASE_EXCEPT;

    void ProcessQueues();
    void ProcessLocalQueues();

    // Operator overloading
    // ====================

    // Make a copy
    // -----------
    const DistSparseMatrix<T>& operator=( const DistSparseMatrix<T>& A );
    // TODO: Move assignment

    // Make a copy of a submatrix
    // --------------------------
    DistSparseMatrix<T> operator()
    ( Range<Int> I, Range<Int> J ) const;
    DistSparseMatrix<T> operator()
    ( Range<Int> I, const vector<Int>& J ) const;
    DistSparseMatrix<T> operator()
    ( const vector<Int>& I, Range<Int> J ) const;
    DistSparseMatrix<T> operator()
    ( const vector<Int>& I, const vector<Int>& J ) const;

    // Rescaling
    // ---------
    const DistSparseMatrix<T>& operator*=( T alpha );

    // Addition/subtraction
    // --------------------
    const DistSparseMatrix<T>& operator+=( const DistSparseMatrix<T>& A );
    const DistSparseMatrix<T>& operator-=( const DistSparseMatrix<T>& A );

    // For manually accessing/modifying buffers
    // ----------------------------------------
    void ForceNumLocalEntries( Int numLocalEntries );
    void ForceConsistency( bool consistent=true ) EL_NO_EXCEPT;
    Int* SourceBuffer() EL_NO_EXCEPT;
    Int* TargetBuffer() EL_NO_EXCEPT;
    Int* OffsetBuffer() EL_NO_EXCEPT;
    T* ValueBuffer() EL_NO_EXCEPT;
    const Int* LockedSourceBuffer() const EL_NO_EXCEPT;
    const Int* LockedTargetBuffer() const EL_NO_EXCEPT;
    const Int* LockedOffsetBuffer() const EL_NO_EXCEPT;
    const T* LockedValueBuffer() const EL_NO_EXCEPT;

    // Queries
    // =======

    // High-level information
    // ----------------------
    Int Height() const EL_NO_EXCEPT;
    Int Width() const EL_NO_EXCEPT;
    Int NumEntries() const EL_NO_EXCEPT;
    El::DistGraph& DistGraph() EL_NO_EXCEPT;
    const El::DistGraph& LockedDistGraph() const EL_NO_EXCEPT;
    Int FirstLocalRow() const EL_NO_EXCEPT;
    Int LocalHeight() const EL_NO_EXCEPT;
    Int NumLocalEntries() const EL_NO_EXCEPT;
    Int Capacity() const EL_NO_EXCEPT;
    bool LocallyConsistent() const EL_NO_EXCEPT;

    // Distribution information
    // ------------------------
    mpi::Comm Comm() const EL_NO_EXCEPT;
    Int Blocksize() const EL_NO_EXCEPT;
    int RowOwner( Int i ) const EL_NO_RELEASE_EXCEPT;
    Int GlobalRow( Int iLoc ) const EL_NO_RELEASE_EXCEPT;
    Int LocalRow( Int i ) const EL_NO_RELEASE_EXCEPT;

    // Detailed local information
    // --------------------------
    Int Row( Int localInd ) const EL_NO_RELEASE_EXCEPT;
    Int Col( Int localInd ) const EL_NO_RELEASE_EXCEPT;
    T Value( Int localInd ) const EL_NO_RELEASE_EXCEPT;
    T GetLocal( Int row, Int col ) const EL_NO_RELEASE_EXCEPT;
    void Set( Int row, Int col, T val ) EL_NO_RELEASE_EXCEPT;
    Int RowOffset( Int localRow ) const EL_NO_RELEASE_EXCEPT;
    Int Offset( Int localRow, Int col ) const EL_NO_RELEASE_EXCEPT;
    Int NumConnections( Int localRow ) const EL_NO_RELEASE_EXCEPT;

    // Return the ratio of the maximum number of local nonzeros to the 
    // total number of nonzeros divided by the number of processes
    double Imbalance() const EL_NO_RELEASE_EXCEPT;

    mutable DistSparseMultMeta multMeta;
    DistSparseMultMeta InitializeMultMeta() const;

    void MappedSources
    ( const DistMap& reordering, vector<Int>& mappedSources ) const;
    void MappedTargets
    ( const DistMap& reordering, 
      vector<Int>& mappedTargets, vector<Int>& colOffs ) const;

    void AssertConsistent() const;
    void AssertLocallyConsistent() const;

private:
    El::DistGraph distGraph_;
    vector<T> vals_;
    vector<T> remoteVals_;

    void InitializeLocalData();

    static bool CompareEntries( const Entry<T>& a, const Entry<T>& b );

    template<typename U,typename V>
    friend void EntrywiseMap
    ( const DistSparseMatrix<U>& A, DistSparseMatrix<V>& B, 
      function<V(U)> func );

    template<typename U> friend class SparseMatrix;
};

} // namespace El

#endif // ifndef EL_CORE_DISTSPARSEMATRIX_DECL_HPP
