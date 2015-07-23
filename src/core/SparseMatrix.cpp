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
SparseMatrix<T>::SparseMatrix( Int height, Int width )
: graph_(height,width)
{ }

template<typename T>
SparseMatrix<T>::SparseMatrix( const SparseMatrix<T>& A )
{ 
    DEBUG_ONLY(CSE cse("SparseMatrix::SparseMatrix"))
    if( &A != this )
        *this = A;
    else
        LogicError("Tried to construct sparse matrix with itself");
}

template<typename T>
SparseMatrix<T>::SparseMatrix( const DistSparseMatrix<T>& A )
{ 
    DEBUG_ONLY(CSE cse("SparseMatrix::SparseMatrix"))
    *this = A;
}

template<typename T>
SparseMatrix<T>::~SparseMatrix() { }

// Assignment and reconfiguration
// ==============================

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
    DEBUG_ONLY(CSE cse("SparseMatrix::Resize"))
    if( Height() == height && Width() == width )
        return;
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
void SparseMatrix<T>::FreezeSparsity() { graph_.frozenSparsity_ = true; }
template<typename T>
void SparseMatrix<T>::UnfreezeSparsity() { graph_.frozenSparsity_ = false; }
template<typename T>
bool SparseMatrix<T>::FrozenSparsity() const { return graph_.frozenSparsity_; }

template<typename T>
void SparseMatrix<T>::Update( Int row, Int col, T value )
{
    DEBUG_ONLY(CSE cse("SparseMatrix::Update"))
    QueueUpdate( row, col, value );
    ProcessQueues();
}

template<typename T>
void SparseMatrix<T>::Update( const Entry<T>& entry )
{ Update( entry.i, entry.j, entry.value ); }

template<typename T>
void SparseMatrix<T>::Zero( Int row, Int col )
{
    DEBUG_ONLY(CSE cse("SparseMatrix::Zero"))
    QueueZero( row, col );
    ProcessQueues();
}

template<typename T>
void SparseMatrix<T>::QueueUpdate( Int row, Int col, T value )
{
    DEBUG_ONLY(CSE cse("SparseMatrix::QueueUpdate"))
    if( FrozenSparsity() )
    {
        const Int offset = Offset( row, col );
        vals_[offset] += value;
    }
    else
    {
        graph_.QueueConnection( row, col );
        vals_.push_back( value );
    }
}

template<typename T>
void SparseMatrix<T>::QueueUpdate( const Entry<T>& entry )
{ QueueUpdate( entry.i, entry.j, entry.value ); }

template<typename T>
void SparseMatrix<T>::QueueZero( Int row, Int col )
{
    DEBUG_ONLY(CSE cse("SparseMatrix::QueueUpdate"))
    if( FrozenSparsity() )
    {
        const Int offset = Offset( row, col );
        vals_[offset] = 0;
    }
    else
    {
        graph_.QueueDisconnection( row, col );
    }
}

// Operator overloading
// ====================

// Make a copy
// -----------
template<typename T>
const SparseMatrix<T>& SparseMatrix<T>::operator=( const SparseMatrix<T>& A )
{
    DEBUG_ONLY(CSE cse("SparseMatrix::operator="))
    graph_ = A.graph_;
    vals_ = A.vals_;
    return *this;
}

template<typename T>
const SparseMatrix<T>&
SparseMatrix<T>::operator=( const DistSparseMatrix<T>& A )
{
    DEBUG_ONLY(CSE cse("SparseMatrix::operator="))
    mpi::Comm comm = A.Comm();
    const int commSize = mpi::Size( comm );
    if( commSize != 1 )
        LogicError("Can not yet construct from distributed sparse matrix");

    graph_ = A.distGraph_;
    vals_ = A.vals_;
    return *this;
}

// Make a copy of a submatrix
// --------------------------
template<typename T>
SparseMatrix<T>
SparseMatrix<T>::operator()( Range<Int> I, Range<Int> J ) const
{
    DEBUG_ONLY(CSE cse("SparseMatrix::operator()"))
    return GetSubmatrix( *this, I, J );
}

// Rescaling
// ---------
template<typename T>
const SparseMatrix<T>& SparseMatrix<T>::operator*=( T alpha )
{
    DEBUG_ONLY(CSE cse("SparseMatrix::operator*=( T )"))
    Scale( alpha, *this );
    return *this;
}

// Addition/subtraction
// --------------------
template<typename T>
const SparseMatrix<T>& SparseMatrix<T>::operator+=( const SparseMatrix<T>& A )
{
    DEBUG_ONLY(CSE cse("SparseMatrix::operator+=( const SM& )"))
    Axpy( T(1), A, *this );
    return *this;
}

template<typename T>
const SparseMatrix<T>& SparseMatrix<T>::operator-=( const SparseMatrix<T>& A )
{
    DEBUG_ONLY(CSE cse("SparseMatrix::operator-=( const SM& )"))
    Axpy( T(-1), A, *this );
    return *this;
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
    DEBUG_ONLY(CSE cse("SparseMatrix::NumEntries"))
    return graph_.NumEdges();
}

template<typename T>
Int SparseMatrix<T>::Capacity() const
{
    DEBUG_ONLY(CSE cse("SparseMatrix::Capacity"))
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
    DEBUG_ONLY(CSE cse("SparseMatrix::Row"))
    return graph_.Source( index );
}

template<typename T>
Int SparseMatrix<T>::Col( Int index ) const
{ 
    DEBUG_ONLY(CSE cse("SparseMatrix::Col"))
    return graph_.Target( index );
}

template<typename T>
Int SparseMatrix<T>::RowOffset( Int row ) const
{
    DEBUG_ONLY(CSE cse("SparseMatrix::RowOffset"))
    return graph_.SourceOffset( row );
}

template<typename T>
Int SparseMatrix<T>::Offset( Int row, Int col ) const
{
    DEBUG_ONLY(CSE cse("SparseMatrix::Offset"))
    return graph_.Offset( row, col );
}

template<typename T>
Int SparseMatrix<T>::NumConnections( Int row ) const
{
    DEBUG_ONLY(CSE cse("SparseMatrix::NumConnections"))
    return graph_.NumConnections( row );
}

template<typename T>
T SparseMatrix<T>::Value( Int index ) const
{ 
    DEBUG_ONLY(
      CSE cse("SparseMatrix::Value");
      if( index < 0 || index >= Int(vals_.size()) )
          LogicError("Entry number out of bounds");
    )
    return vals_[index];
}

template<typename T>
Int* SparseMatrix<T>::SourceBuffer() { return graph_.SourceBuffer(); }
template<typename T>
Int* SparseMatrix<T>::TargetBuffer() { return graph_.TargetBuffer(); }
template<typename T>
Int* SparseMatrix<T>::OffsetBuffer() { return graph_.OffsetBuffer(); }
template<typename T>
T* SparseMatrix<T>::ValueBuffer() { return vals_.data(); }

template<typename T>
const Int* SparseMatrix<T>::LockedSourceBuffer() const
{ return graph_.LockedSourceBuffer(); }
template<typename T>
const Int* SparseMatrix<T>::LockedTargetBuffer() const
{ return graph_.LockedTargetBuffer(); }
template<typename T>
const Int* SparseMatrix<T>::LockedOffsetBuffer() const 
{ return graph_.LockedOffsetBuffer(); }
template<typename T>
const T* SparseMatrix<T>::LockedValueBuffer() const
{ return vals_.data(); }

// Auxiliary routines
// ==================

template<typename T>
bool SparseMatrix<T>::CompareEntries( const Entry<T>& a, const Entry<T>& b )
{ return a.i < b.i || (a.i == b.i && a.j < b.j); }

template<typename T>
void SparseMatrix<T>::ProcessQueues()
{
    DEBUG_ONLY(
      CSE cse("SparseMatrix::ProcessQueues");
      if( graph_.sources_.size() != graph_.targets_.size() || 
          graph_.targets_.size() != vals_.size() )
          LogicError("Inconsistent sparse matrix buffer sizes");
    )
    if( graph_.consistent_ )
        return;

    Int numRemoved = 0;
    const Int numEntries = vals_.size();
    vector<Entry<T>> entries( numEntries );
    if( graph_.markedForRemoval_.size() != 0 )
    {
        for( Int s=0; s<numEntries; ++s )
        {
            pair<Int,Int> candidate(graph_.sources_[s],graph_.targets_[s]);
            if( graph_.markedForRemoval_.find(candidate) == 
                graph_.markedForRemoval_.end() )
            {
                entries[s-numRemoved].i = graph_.sources_[s];
                entries[s-numRemoved].j = graph_.targets_[s];
                entries[s-numRemoved].value = vals_[s];
            }
            else
            {
                ++numRemoved;
            }
        }
        graph_.markedForRemoval_.clear();
        entries.resize( numEntries-numRemoved );
    }
    else
    {
        for( Int s=0; s<numEntries; ++s )
            entries[s] = 
              Entry<T>{graph_.sources_[s],graph_.targets_[s],vals_[s]};
    }
    std::sort( entries.begin(), entries.end(), CompareEntries );
    const Int numSorted = entries.size();

    // Compress out duplicates
    Int lastUnique=0;
    for( Int s=1; s<numSorted; ++s )
    {
        if( entries[s].i != entries[lastUnique].i ||
            entries[s].j != entries[lastUnique].j )
            entries[++lastUnique] = entries[s];
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
        graph_.sources_[s] = entries[s].i;
        graph_.targets_[s] = entries[s].j;
        vals_[s] = entries[s].value;
    }

    graph_.ComputeSourceOffsets();
    graph_.consistent_ = true;
}

template<typename T>
void SparseMatrix<T>::AssertConsistent() const
{ graph_.AssertConsistent(); }

#define PROTO(T) template class SparseMatrix<T>;
#define EL_ENABLE_QUAD
#include "El/macros/Instantiate.h"

} // namespace El
