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
    ~DistMultiVec();

    // Assignment  and reconfiguration
    // ===============================

    // Make a copy of a submatrix
    // --------------------------
    DistMultiVec<T> operator()( Range<Int> I, Range<Int> J ) const;
   
    // Assignment
    // ----------
    const DistMultiVec<T>& operator=( const DistMultiVec<T>& X );
    const DistMultiVec<T>& operator=( const AbstractDistMatrix<T>& X );

    // Change the matrix size
    // ----------------------
    void Empty();
    void Resize( Int height, Int width );

    // Change the distribution
    // -----------------------
    void SetComm( mpi::Comm comm );

    // Queries
    // =======

    // High-level data
    // ---------------
    Int Height() const;
    Int Width() const;
    Int FirstLocalRow() const;
    Int LocalHeight() const;
          El::Matrix<T>& Matrix();
    const El::Matrix<T>& LockedMatrix() const;

    // Distribution information
    // ------------------------
    mpi::Comm Comm() const;
    Int Blocksize() const;
    int RowOwner( Int i ) const;
    Int GlobalRow( Int iLoc ) const;

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

private:
    Int height_, width_;

    mpi::Comm comm_;

    Int blocksize_;
    Int firstLocalRow_;

    El::Matrix<T> multiVec_;

    void InitializeLocalData();
};

} // namespace El

#endif // ifndef EL_CORE_DISTMULTIVEC_DECL_HPP
