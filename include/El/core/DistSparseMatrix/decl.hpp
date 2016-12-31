/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   Copyright (c) 2012 Jack Poulson, Lexing Ying, and
   The University of Texas at Austin.
   All rights reserved.

   Copyright (c) 2013 Jack Poulson, Lexing Ying, and Stanford University.
   All rights reserved.

   Copyright (c) 2014 Jack Poulson and The Georgia Institute of Technology.
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_CORE_DISTSPARSEMATRIX_DECL_HPP
#define EL_CORE_DISTSPARSEMATRIX_DECL_HPP

namespace El  {


// Use a simple 1d distribution where each process owns a fixed number of rows,
//     if last process,  height - (commSize-1)*floor(height/commSize)
//     otherwise,        floor(height/commSize)
template<typename Ring>
class DistSparseMatrix
{
public:
    // Constructors and destructors
    // ============================
    DistSparseMatrix
    ( const El::Grid& grid=El::Grid::Default() );
    DistSparseMatrix
    ( Int height, Int width, const El::Grid& grid=El::Grid::Default() );
    DistSparseMatrix( const DistSparseMatrix<Ring>& A );
    // TODO(poulson): Move constructor
    ~DistSparseMatrix();

    // Assignment and reconfiguration
    // ==============================

    // Change the size of the matrix
    // -----------------------------
    void Empty( bool freeMemory=true );
    void Resize( Int height, Int width );

    // Change the distribution
    // -----------------------
    void SetGrid( const El::Grid& grid );

    // Assembly
    // --------
    void Reserve( Int numLocalEntries, Int numRemoteEntries=0 );

    void FreezeSparsity() EL_NO_EXCEPT;
    void UnfreezeSparsity() EL_NO_EXCEPT;
    bool FrozenSparsity() const EL_NO_EXCEPT;

    // Expensive independent updates and explicit zeroing
    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    void Update( const Entry<Ring>& entry );
    void Update( Int row, Int col, const Ring& value );
    void Zero( Int row, Int col );

    void UpdateLocal( const Entry<Ring>& localEntry );
    void UpdateLocal( Int localRow, Int col, const Ring& value );
    void ZeroLocal( Int localRow, Int col );

    // Batch updating and zeroing (recommended)
    // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    void QueueUpdate( const Entry<Ring>& entry, bool passive=false )
    EL_NO_RELEASE_EXCEPT;
    void QueueUpdate( Int row, Int col, const Ring& value, bool passive=false )
    EL_NO_RELEASE_EXCEPT;
    void QueueZero( Int row, Int col, bool passive=false )
    EL_NO_RELEASE_EXCEPT;

    void QueueLocalUpdate( const Entry<Ring>& localEntry )
    EL_NO_RELEASE_EXCEPT;
    void QueueLocalUpdate( Int localRow, Int col, const Ring& value )
    EL_NO_RELEASE_EXCEPT;
    void QueueLocalZero( Int localRow, Int col )
    EL_NO_RELEASE_EXCEPT;

    void ProcessQueues();
    void ProcessLocalQueues();

    // Operator overloading
    // ====================

    // Make a copy
    // -----------
    const DistSparseMatrix<Ring>& operator=( const DistSparseMatrix<Ring>& A );
    // TODO(poulson): Move assignment

    // Make a copy of a submatrix
    // --------------------------
    DistSparseMatrix<Ring> operator()
    ( Range<Int> I, Range<Int> J ) const;
    DistSparseMatrix<Ring> operator()
    ( Range<Int> I, const vector<Int>& J ) const;
    DistSparseMatrix<Ring> operator()
    ( const vector<Int>& I, Range<Int> J ) const;
    DistSparseMatrix<Ring> operator()
    ( const vector<Int>& I, const vector<Int>& J ) const;

    // Rescaling
    // ---------
    const DistSparseMatrix<Ring>& operator*=( const Ring& alpha );

    // Addition/subtraction
    // --------------------
    const DistSparseMatrix<Ring>& operator+=( const DistSparseMatrix<Ring>& A );
    const DistSparseMatrix<Ring>& operator-=( const DistSparseMatrix<Ring>& A );

    // For manually accessing/modifying buffers
    // ----------------------------------------
    void ForceNumLocalEntries( Int numLocalEntries );
    void ForceConsistency( bool consistent=true ) EL_NO_EXCEPT;
    Int* SourceBuffer() EL_NO_EXCEPT;
    Int* TargetBuffer() EL_NO_EXCEPT;
    Int* OffsetBuffer() EL_NO_EXCEPT;
    Ring* ValueBuffer() EL_NO_EXCEPT;
    const Int* LockedSourceBuffer() const EL_NO_EXCEPT;
    const Int* LockedTargetBuffer() const EL_NO_EXCEPT;
    const Int* LockedOffsetBuffer() const EL_NO_EXCEPT;
    const Ring* LockedValueBuffer() const EL_NO_EXCEPT;

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
    const El::Grid& Grid() const EL_NO_EXCEPT;
    Int Blocksize() const EL_NO_EXCEPT;
    int RowOwner( Int i ) const EL_NO_RELEASE_EXCEPT;
    Int GlobalRow( Int iLoc ) const EL_NO_RELEASE_EXCEPT;
    Int LocalRow( Int i ) const EL_NO_RELEASE_EXCEPT;

    // Detailed local information
    // --------------------------
    Int Row( Int localInd ) const EL_NO_RELEASE_EXCEPT;
    Int Col( Int localInd ) const EL_NO_RELEASE_EXCEPT;
    Ring Value( Int localInd ) const EL_NO_RELEASE_EXCEPT;
    Ring GetLocal( Int row, Int col ) const EL_NO_RELEASE_EXCEPT;
    void Set( Int row, Int col, const Ring& val ) EL_NO_RELEASE_EXCEPT;
    Int RowOffset( Int localRow ) const EL_NO_RELEASE_EXCEPT;
    Int Offset( Int localRow, Int col ) const EL_NO_RELEASE_EXCEPT;
    Int NumConnections( Int localRow ) const EL_NO_RELEASE_EXCEPT;

    // Return the ratio of the maximum number of local nonzeros to the
    // total number of nonzeros divided by the number of processes
    double Imbalance() const EL_NO_RELEASE_EXCEPT;

    DistGraphMultMeta InitializeMultMeta() const;

    void MappedSources
    ( const DistMap& reordering, vector<Int>& mappedSources ) const;
    void MappedTargets
    ( const DistMap& reordering,
      vector<Int>& mappedTargets, vector<Int>& colOffs ) const;

    void AssertConsistent() const;
    void AssertLocallyConsistent() const;

private:
    El::DistGraph distGraph_;
    vector<Ring> vals_;
    vector<Ring> remoteVals_;

    void InitializeLocalData();

    static bool CompareEntries( const Entry<Ring>& a, const Entry<Ring>& b );

    template<typename U> friend class SparseMatrix;
};

} // namespace El

#endif // ifndef EL_CORE_DISTSPARSEMATRIX_DECL_HPP
