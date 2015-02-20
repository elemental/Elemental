/*
   Copyright (c) 2009-2015, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, Stanford University, and the
   Georgia Insitute of Technology.
   All rights reserved.
 
   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp" 
namespace El {

// Constructors and destructors
// ============================

template<typename T>
SparseMatrix<T>::SparseMatrix() { }

template<typename T>
SparseMatrix<T>::SparseMatrix( Int height )
: graph_(height)
{ }

template<typename T>
SparseMatrix<T>::SparseMatrix( Int height, Int width )
: graph_(height,width)
{ }

template<typename T>
SparseMatrix<T>::SparseMatrix( const SparseMatrix<T>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("SparseMatrix::SparseMatrix"))
    if( &A != this )
        *this = A;
    else
        LogicError("Tried to construct sparse matrix with itself");
}

template<typename T>
SparseMatrix<T>::SparseMatrix( const DistSparseMatrix<T>& A )
{ 
    DEBUG_ONLY(CallStackEntry cse("SparseMatrix::SparseMatrix"))
    *this = A;
}

template<typename T>
SparseMatrix<T>::~SparseMatrix() { }

// Assignment and reconfiguration
// ==============================

// Make a copy
// -----------
template<typename T>
const SparseMatrix<T>& SparseMatrix<T>::operator=( const SparseMatrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("SparseMatrix::operator="))
    graph_ = A.graph_;
    vals_ = A.vals_;
    return *this;
}

template<typename T>
const SparseMatrix<T>&
SparseMatrix<T>::operator=( const DistSparseMatrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("SparseMatrix::operator="))
    mpi::Comm comm = A.Comm();
    const int commSize = mpi::Size( comm );
    if( commSize != 1 )
        LogicError("Can not yet construct from distributed sparse matrix");

    graph_ = A.distGraph_;
    vals_ = A.vals_;
    return *this;
}

// Change the size of the matrix
// -----------------------------
template<typename T>
void SparseMatrix<T>::Empty( bool clearMemory )
{
    graph_.Empty( clearMemory );
    if( clearMemory )
        SwapClear( vals_ );
    else
        vals_.resize( 0 );
}

template<typename T>
void SparseMatrix<T>::Resize( Int height, Int width )
{
    graph_.Resize( height, width );
    vals_.resize( 0 );
}

// Assembly
// --------
template<typename T>
void SparseMatrix<T>::Reserve( Int numEntries )
{ 
    graph_.Reserve( numEntries );
    vals_.reserve( numEntries );
}

template<typename T>
void SparseMatrix<T>::Update( Int row, Int col, T value )
{
    DEBUG_ONLY(CallStackEntry cse("SparseMatrix::Update"))
    QueueUpdate( row, col, value );
    MakeConsistent();
}

template<typename T>
void SparseMatrix<T>::Zero( Int row, Int col )
{
    DEBUG_ONLY(CallStackEntry cse("SparseMatrix::Zero"))
    QueueZero( row, col );
    MakeConsistent();
}

template<typename T>
void SparseMatrix<T>::QueueUpdate( Int row, Int col, T value )
{
    DEBUG_ONLY(CallStackEntry cse("SparseMatrix::QueueUpdate"))
    graph_.QueueConnection( row, col );
    vals_.push_back( value );
}

template<typename T>
void SparseMatrix<T>::QueueZero( Int row, Int col )
{
    DEBUG_ONLY(CallStackEntry cse("SparseMatrix::QueueUpdate"))
    graph_.QueueDisconnection( row, col );
}

// Queries
// =======

// High-level information
// ----------------------
template<typename T>
Int SparseMatrix<T>::Height() const { return graph_.NumSources(); }
template<typename T>
Int SparseMatrix<T>::Width() const { return graph_.NumTargets(); }

template<typename T>
Int SparseMatrix<T>::NumEntries() const
{
    DEBUG_ONLY(CallStackEntry cse("SparseMatrix::NumEntries"))
    return graph_.NumEdges();
}

template<typename T>
Int SparseMatrix<T>::Capacity() const
{
    DEBUG_ONLY(CallStackEntry cse("SparseMatrix::Capacity"))
    return graph_.Capacity();
}

template<typename T>
bool SparseMatrix<T>::Consistent() const { return graph_.Consistent(); }

template<typename T>
El::Graph& SparseMatrix<T>::Graph() { return graph_; }
template<typename T>
const El::Graph& SparseMatrix<T>::LockedGraph() const { return graph_; }

// Entrywise information
// ---------------------
template<typename T>
Int SparseMatrix<T>::Row( Int index ) const
{ 
    DEBUG_ONLY(CallStackEntry cse("SparseMatrix::Row"))
    return graph_.Source( index );
}

template<typename T>
Int SparseMatrix<T>::Col( Int index ) const
{ 
    DEBUG_ONLY(CallStackEntry cse("SparseMatrix::Col"))
    return graph_.Target( index );
}

template<typename T>
Int SparseMatrix<T>::EntryOffset( Int row ) const
{
    DEBUG_ONLY(CallStackEntry cse("SparseMatrix::EntryOffset"))
    return graph_.EdgeOffset( row );
}

template<typename T>
Int SparseMatrix<T>::NumConnections( Int row ) const
{
    DEBUG_ONLY(CallStackEntry cse("SparseMatrix::NumConnections"))
    return graph_.NumConnections( row );
}

template<typename T>
T SparseMatrix<T>::Value( Int index ) const
{ 
    DEBUG_ONLY(
      CallStackEntry cse("SparseMatrix::Value");
      if( index < 0 || index >= vals_.size() )
          LogicError("Entry number out of bounds");
    )
    return vals_[index];
}

template<typename T>
Int* SparseMatrix<T>::SourceBuffer() { return graph_.SourceBuffer(); }
template<typename T>
Int* SparseMatrix<T>::TargetBuffer() { return graph_.TargetBuffer(); }
template<typename T>
T* SparseMatrix<T>::ValueBuffer() { return vals_.data(); }

template<typename T>
const Int* SparseMatrix<T>::LockedSourceBuffer() const
{ return graph_.LockedSourceBuffer(); }
template<typename T>
const Int* SparseMatrix<T>::LockedTargetBuffer() const
{ return graph_.LockedTargetBuffer(); }
template<typename T>
const T* SparseMatrix<T>::LockedValueBuffer() const
{ return vals_.data(); }

// Auxiliary routines
// ==================

template<typename T>
bool SparseMatrix<T>::CompareEntries( const Entry<T>& a, const Entry<T>& b )
{ return a.indices[0] < b.indices[0] || 
         (a.indices[0] == b.indices[0] && a.indices[1] < b.indices[1]); }

template<typename T>
void SparseMatrix<T>::MakeConsistent()
{
    DEBUG_ONLY(
      CallStackEntry cse("SparseMatrix::MakeConsistent");
      if( graph_.sources_.size() != graph_.targets_.size() || 
          graph_.targets_.size() != vals_.size() )
          LogicError("Inconsistent sparse matrix buffer sizes");
    )
    if( !graph_.consistent_ )
    {
        const Int numEntries = vals_.size();
        Int numRemoved = 0;
        vector<Entry<T>> entries( numEntries );
        for( Int s=0; s<numEntries; ++s )
        {
            pair<Int,Int> candidate(graph_.sources_[s],graph_.targets_[s]);
            if( graph_.markedForRemoval_.find(candidate) == 
                graph_.markedForRemoval_.end() )
            {
                entries[s-numRemoved].indices[0] = graph_.sources_[s];
                entries[s-numRemoved].indices[1] = graph_.targets_[s];
                entries[s-numRemoved].value = vals_[s];
            }
            else
            {
                ++numRemoved;
            }
        }
        graph_.markedForRemoval_.clear();
        entries.resize( numEntries-numRemoved );
        std::sort( entries.begin(), entries.end(), CompareEntries );

        // Compress out duplicates
        Int lastUnique=0;
        for( Int s=1; s<numEntries; ++s )
        {
            if( entries[s].indices[0] != entries[lastUnique].indices[0] ||
                entries[s].indices[1] != entries[lastUnique].indices[1] )
            {
                ++lastUnique;
                entries[lastUnique].indices[0] = entries[s].indices[0];
                entries[lastUnique].indices[1] = entries[s].indices[1];
                entries[lastUnique].value = entries[s].value;
            }
            else
                entries[lastUnique].value += entries[s].value;
        }
        const Int numUnique = lastUnique+1;
        entries.resize( numUnique );

        graph_.sources_.resize( numUnique );
        graph_.targets_.resize( numUnique );
        vals_.resize( numUnique );
        for( Int s=0; s<numUnique; ++s )
        {
            graph_.sources_[s] = entries[s].indices[0];
            graph_.targets_[s] = entries[s].indices[1];
            vals_[s] = entries[s].value;
        }

        graph_.ComputeEdgeOffsets();

        graph_.consistent_ = true;
    }
}

template<typename T>
void SparseMatrix<T>::AssertConsistent() const
{ graph_.AssertConsistent(); }

#define PROTO(T) template class SparseMatrix<T>;
#include "El/macros/Instantiate.h"

} // namespace El
