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
    DistMultiVec();
    DistMultiVec( mpi::Comm comm );
    DistMultiVec( Int height, Int width, mpi::Comm comm );
    ~DistMultiVec();

    // Assignment  and reconfiguration
    // ===============================
   
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
    T Get( Int row, Int col ) const;
    T GetLocal( Int localRow, Int col ) const;
    void SetLocal( Int localRow, Int col, T value );
    void UpdateLocal( Int localRow, Int col, T value );

private:
    Int height_, width_;

    mpi::Comm comm_;

    Int blocksize_;
    Int firstLocalRow_;

    El::Matrix<T> multiVec_;
};

} // namespace El

#endif // ifndef EL_CORE_DISTMULTIVEC_DECL_HPP
