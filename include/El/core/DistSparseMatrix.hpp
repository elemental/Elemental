/*
   Copyright (c) 2009-2015, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, Stanford University, and the
   Georgia Insitute of Technology.
   All rights reserved.
 
   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_CORE_DISTSPARSEMATRIX_DECL_HPP
#define EL_CORE_DISTSPARSEMATRIX_DECL_HPP

namespace El {

template<typename T>
struct SparseMultMeta
{
    bool ready;
    // NOTE: The 'send' and 'recv' roles reverse for adjoint multiplication
    Int numRecvInds;
    vector<int> sendSizes, sendOffs,
                recvSizes, recvOffs;
    vector<Int> sendInds, colOffs;

    SparseMultMeta() : ready(false), numRecvInds(0) { }

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
    // TODO: Constructor for building from another DistSparseMatrix
    // TODO: Move constructor
    ~DistSparseMatrix();

    // Assignment and reconfiguration
    // ==============================

    // Make a copy
    // -----------
    // TODO: operator=
    // TODO: Move assignment

    // Make a copy of a submatrix
    // --------------------------
    DistSparseMatrix<T> operator()( Range<Int> I, Range<Int> J ) const;

    // Change the size of the matrix
    // -----------------------------
    void Empty( bool clearMemory=true );
    void Resize( Int height, Int width );

    // Change the distribution
    // -----------------------
    void SetComm( mpi::Comm comm );

    // Assembly
    // --------
    void Reserve( Int numLocalEntries, Int numRemoteEntries=0 );

    // A safe procedure for applying a local update or zeroing an entry
    void Update( const Entry<T>& entry, bool passive=true );
    void Update( Int row, Int col, T value, bool passive=true );
    void UpdateLocal( const Entry<T>& localEntry );
    void UpdateLocal( Int localRow, Int col, T value );
    void Zero( Int row, Int col, bool passive=true );
    void ZeroLocal( Int localRow, Int col );

    // For applying a sequence of updates and then forcing consistency
    void QueueUpdate( const Entry<T>& entry, bool passive=true );
    void QueueUpdate( Int row, Int col, T value, bool passive=true );
    void QueueLocalUpdate( const Entry<T>& localEntry );
    void QueueLocalUpdate( Int localRow, Int col, T value );
    void QueueZero( Int row, Int col, bool passive=true );
    void QueueLocalZero( Int localRow, Int col );
    void ProcessQueues();

    // Queries
    // =======

    // High-level information
    // ----------------------
    Int Height() const;
    Int Width() const;
    El::DistGraph& DistGraph();
    const El::DistGraph& LockedDistGraph() const;
    Int FirstLocalRow() const;
    Int LocalHeight() const;
    Int NumLocalEntries() const;
    Int Capacity() const;
    bool Consistent() const;

    // Distribution information
    // ------------------------
    mpi::Comm Comm() const;
    Int Blocksize() const;
    int RowOwner( Int i ) const;
    Int GlobalRow( Int iLoc ) const;

    // Detailed local information
    // --------------------------
    Int Row( Int localInd ) const;
    Int Col( Int localInd ) const;
    T Value( Int localInd ) const;
    Int EntryOffset( Int localRow ) const;
    Int NumConnections( Int localRow ) const;
    Int* SourceBuffer();
    Int* TargetBuffer();
    T* ValueBuffer();
    const Int* LockedSourceBuffer() const;
    const Int* LockedTargetBuffer() const;
    const T* LockedValueBuffer() const;

    mutable SparseMultMeta<T> multMeta;

private:
    El::DistGraph distGraph_;
    vector<T> vals_;
    vector<T> remoteVals_;

    void InitializeLocalData();

    static bool CompareEntries( const Entry<T>& a, const Entry<T>& b );

    void AssertConsistent() const;

    template<typename U> friend class SparseMatrix;
    template<typename U> friend struct DistSymmFront;

    template<typename U> friend void Syrk
    ( Orientation orientation, 
      U alpha, const DistSparseMatrix<U>& A, 
      U beta,        DistSparseMatrix<U>& C, bool conjugate );
};

} // namespace El

#endif // ifndef EL_CORE_DISTSPARSEMATRIX_DECL_HPP
