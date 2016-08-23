/*
   Copyright (c) 2009-2016, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, Stanford University, and the
   Georgia Insitute of Technology.
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
template<typename T>
class DistMultiVec
{
public:
    // Constructors and destructors
    // ============================
    DistMultiVec( mpi::Comm comm=mpi::COMM_WORLD );
    DistMultiVec( Int height, Int width, mpi::Comm comm=mpi::COMM_WORLD );
    DistMultiVec( const DistMultiVec<T>& A );
    ~DistMultiVec();

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
    DistMultiVec<T> operator()
    ( Range<Int> I, Range<Int> J ) const;
    DistMultiVec<T> operator()
    ( Range<Int> I, const vector<Int>& J ) const;
    DistMultiVec<T> operator()
    ( const vector<Int>& I, Range<Int> J ) const;
    DistMultiVec<T> operator()
    ( const vector<Int>& I, const vector<Int>& J ) const;
   
    // Assignment
    // ----------
    const DistMultiVec<T>& operator=( const DistMultiVec<T>& X );
    const DistMultiVec<T>& operator=( const AbstractDistMatrix<T>& X );

    // Rescaling
    // ---------
    const DistMultiVec<T>& operator*=( T alpha );

    // Addition/subtraction
    // --------------------
    const DistMultiVec<T>& operator+=( const DistMultiVec<T>& A );
    const DistMultiVec<T>& operator-=( const DistMultiVec<T>& A );

    // Queries
    // =======

    // High-level data
    // ---------------
    Int Height() const EL_NO_EXCEPT;
    Int Width() const EL_NO_EXCEPT;
    Int FirstLocalRow() const EL_NO_EXCEPT;
    Int LocalHeight() const EL_NO_EXCEPT;
          El::Matrix<T>& Matrix() EL_NO_EXCEPT;
    const El::Matrix<T>& LockedMatrix() const EL_NO_EXCEPT;

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
    T Get( Int i, Int j ) const;
    T GetLocal( Int iLoc, Int j ) const;
    void Set( Int i, Int j, T value );
    void Set( const Entry<T>& entry );
    void SetLocal( Int iLoc, Int j, T value );
    void SetLocal( const Entry<T>& localEntry );
    void Update( Int i, Int j, T value );
    void Update( const Entry<T>& entry );
    void UpdateLocal( Int iLoc, Int j, T value );
    void UpdateLocal( const Entry<T>& entry );

    // Batch updating of remote entries
    // --------------------------------
    void Reserve( Int numRemoteEntries );
    void QueueUpdate( const Entry<T>& entry );
    void QueueUpdate( Int i, Int j, T value );
    void ProcessQueues();

private:
    Int height_, width_;

    mpi::Comm comm_;
    // Calling MPI_Comm_size within an inner loop is apparently a bad idea
    int commSize_;
    int commRank_;
    Int blocksize_;

    El::Matrix<T> multiVec_;

    // Remote updates
    // --------------
    vector<Entry<T>> remoteUpdates_;

    void InitializeLocalData();
};

} // namespace El

#endif // ifndef EL_CORE_DISTMULTIVEC_DECL_HPP
