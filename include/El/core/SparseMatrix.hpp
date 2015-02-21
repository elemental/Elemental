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
#ifndef EL_CORE_SPARSEMATRIX_DECL_HPP
#define EL_CORE_SPARSEMATRIX_DECL_HPP

namespace El {

// Forward declaration for constructor
template<typename T> class DistSparseMatrix;

template<typename T>
class SparseMatrix
{
public:
    // Constructors and destructors
    // ============================
    SparseMatrix();
    SparseMatrix( Int height );
    SparseMatrix( Int height, Int width );
    SparseMatrix( const SparseMatrix<T>& A );
    // NOTE: This requires A to be distributed over a single process
    SparseMatrix( const DistSparseMatrix<T>& A );
    // TODO: Move constructor
    ~SparseMatrix();

    // Assignment and reconfiguration
    // ==============================

    // Make a copy
    // -----------
    // For copying one matrix to another
    const SparseMatrix<T>& operator=( const SparseMatrix<T>& A );
    // NOTE: This requires A to be distributed over a single process
    const SparseMatrix<T>& operator=( const DistSparseMatrix<T>& A );
    // TODO: Move assignment

    // Change the size of the matrix
    // -----------------------------
    void Empty( bool clearMemory=true );
    void Resize( Int height, Int width );

    // Assembly
    // --------
    void Reserve( Int numEntries );

    // A safe update/zeroing mechanism
    void Update( Int row, Int col, T value );
    void Zero( Int row, Int col );

    // For performing many updates and forcing consistency at the end
    void QueueUpdate( Int row, Int col, T value );
    void QueueZero( Int row, Int col );
    void MakeConsistent();

    // Queries
    // =======

    // High-level information
    // ----------------------
    Int Height() const;
    Int Width() const;
    Int NumEntries() const;
    Int Capacity() const;
    bool Consistent() const;
    El::Graph& Graph();
    const El::Graph& LockedGraph() const;

    // Entrywise information
    // ---------------------
    Int Row( Int index ) const;
    Int Col( Int index ) const;
    T Value( Int index ) const;
    Int EntryOffset( Int row ) const;
    Int NumConnections( Int row ) const;
    Int* SourceBuffer();
    Int* TargetBuffer();
    T* ValueBuffer();
    const Int* LockedSourceBuffer() const;
    const Int* LockedTargetBuffer() const;
    const T* LockedValueBuffer() const;

private:
    El::Graph graph_;
    vector<T> vals_;

    static bool CompareEntries( const Entry<T>& a, const Entry<T>& b );

    void AssertConsistent() const;

    template<typename U> friend class DistSparseMatrix;
    template<typename U> 
    friend void CopyFromRoot
    ( const DistSparseMatrix<U>& ADist, SparseMatrix<U>& A );
};

} // namespace El

#endif // ifndef EL_CORE_SPARSEMATRIX_DECL_HPP
