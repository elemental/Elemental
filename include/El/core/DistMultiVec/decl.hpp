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
#ifndef EL_CORE_DISTMULTIVEC_DECL_HPP
#define EL_CORE_DISTMULTIVEC_DECL_HPP

namespace El {

// Use a simple 1d distribution where each process owns a fixed number of rows,
//     if last process,  height - (commSize-1)*floor(height/commSize)
//     otherwise,        floor(height/commSize)
template<typename Ring>
class DistMultiVec
{
public:
    // Constructors and destructors
    // ============================
    DistMultiVec( mpi::Comm comm=mpi::COMM_WORLD );
    DistMultiVec( Int height, Int width, mpi::Comm comm=mpi::COMM_WORLD );
    DistMultiVec( const DistMultiVec<Ring>& A );
    ~DistMultiVec();

    // Advanced
    // --------
    // The following are provided to aid duck-typing over Matrix, DistMatrix,
    // DistMultiVec, etc. and are functionally equivalent to the mpi::Comm
    // equivalents.
    explicit DistMultiVec( const Grid& grid );
    explicit DistMultiVec( Int height, Int width, const Grid& grid ); 

    // Assignment  and reconfiguration
    // ===============================

    // Change the matrix size
    // ----------------------
    void Empty( bool freeMemory=true );
    void Resize( Int height, Int width );

    // Change the distribution
    // -----------------------
    void SetComm( mpi::Comm comm );

    // Operator overloading
    // ====================

    // Make a copy of a submatrix
    // --------------------------
    DistMultiVec<Ring> operator()
    ( Range<Int> I, Range<Int> J ) const;
    DistMultiVec<Ring> operator()
    ( Range<Int> I, const vector<Int>& J ) const;
    DistMultiVec<Ring> operator()
    ( const vector<Int>& I, Range<Int> J ) const;
    DistMultiVec<Ring> operator()
    ( const vector<Int>& I, const vector<Int>& J ) const;

    // Assignment
    // ----------
    const DistMultiVec<Ring>& operator=( const DistMultiVec<Ring>& X );
    const DistMultiVec<Ring>& operator=( const AbstractDistMatrix<Ring>& X );

    // Rescaling
    // ---------
    const DistMultiVec<Ring>& operator*=( const Ring& alpha );

    // Addition/subtraction
    // --------------------
    const DistMultiVec<Ring>& operator+=( const DistMultiVec<Ring>& A );
    const DistMultiVec<Ring>& operator-=( const DistMultiVec<Ring>& A );

    // Queries
    // =======

    // High-level data
    // ---------------
    Int Height() const EL_NO_EXCEPT;
    Int Width() const EL_NO_EXCEPT;
    Int FirstLocalRow() const EL_NO_EXCEPT;
    Int LocalHeight() const EL_NO_EXCEPT;
          El::Matrix<Ring>& Matrix() EL_NO_EXCEPT;
    const El::Matrix<Ring>& LockedMatrix() const EL_NO_EXCEPT;

    // Distribution information
    // ------------------------
    mpi::Comm Comm() const EL_NO_EXCEPT;
    Int Blocksize() const EL_NO_EXCEPT;
    int RowOwner( Int i ) const EL_NO_EXCEPT;
    int Owner( Int i, Int j ) const EL_NO_EXCEPT;
    bool IsLocal( Int i, Int j ) const EL_NO_EXCEPT;
    bool IsLocalRow( Int i ) const EL_NO_EXCEPT;
    Int GlobalRow( Int iLoc ) const;
    Int LocalRow( Int i ) const;

    // Entrywise manipulation
    // ======================
    Ring Get( Int i, Int j ) const;
    Ring GetLocal( Int iLoc, Int j ) const;
    void Set( Int i, Int j, const Ring& value );
    void Set( const Entry<Ring>& entry );
    void SetLocal( Int iLoc, Int j, const Ring& value );
    void SetLocal( const Entry<Ring>& localEntry );
    void Update( Int i, Int j, const Ring& value );
    void Update( const Entry<Ring>& entry );
    void UpdateLocal( Int iLoc, Int j, const Ring& value );
    void UpdateLocal( const Entry<Ring>& entry );

    // Batch updating of remote entries
    // --------------------------------
    void Reserve( Int numRemoteEntries );
    void QueueUpdate( const Entry<Ring>& entry );
    void QueueUpdate( Int i, Int j, const Ring& value );
    void ProcessQueues();

private:
    Int height_=0, width_=0;

    mpi::Comm comm_;
    // Calling MPI_Comm_size within an inner loop is apparently a bad idea
    int commSize_;
    int commRank_;
    Int blocksize_;

    El::Matrix<Ring> multiVec_;

    // Remote updates
    // --------------
    vector<Entry<Ring>> remoteUpdates_;

    void InitializeLocalData();
};

} // namespace El

#endif // ifndef EL_CORE_DISTMULTIVEC_DECL_HPP
