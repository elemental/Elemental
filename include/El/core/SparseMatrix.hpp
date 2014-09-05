/*
   Copyright (c) 2009-2014, Jack Poulson, Lexing Ying,
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
    // Construction and destruction
    SparseMatrix();
    SparseMatrix( int height );
    SparseMatrix( int height, int width );
    SparseMatrix( const SparseMatrix<T>& A );
    // NOTE: This requires A to be distributed over a single process
    SparseMatrix( const DistSparseMatrix<T>& A );
    ~SparseMatrix();

    // High-level information
    int Height() const;
    int Width() const;
    El::Graph& Graph();
    const El::Graph& LockedGraph() const;

    // Assembly-related routines
    void StartAssembly();
    void StopAssembly();
    void Reserve( int numEntries );
    void Update( int row, int col, T value );
    int Capacity() const;

    // Data
    int Row( int index ) const;
    int Col( int index ) const;
    T Value( int index ) const;
    int NumEntries() const;
    int EntryOffset( int row ) const;
    int NumConnections( int row ) const;
    int* SourceBuffer();
    int* TargetBuffer();
    T* ValueBuffer();
    const int* LockedSourceBuffer() const;
    const int* LockedTargetBuffer() const;
    const T* LockedValueBuffer() const;

    // For modifying the size of the matrix
    void Empty();
    void Resize( int height, int width );

    // For copying one matrix to another
    const SparseMatrix<T>& operator=( const SparseMatrix<T>& A );
    // NOTE: This requires A to be distributed over a single process
    const SparseMatrix<T>& operator=( const DistSparseMatrix<T>& A );

private:
    El::Graph graph_;
    std::vector<T> vals_;

    static bool CompareEntries( const Entry<T>& a, const Entry<T>& b );

    void EnsureConsistentSizes() const;
    void EnsureConsistentCapacities() const;

    template<typename U> friend class DistSparseMatrix;
};

} // namespace El

#endif // ifndef EL_CORE_SPARSEMATRIX_DECL_HPP
