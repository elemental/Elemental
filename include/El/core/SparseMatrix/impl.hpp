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
#ifndef EL_SPARSEMATRIX_IMPL_HPP
#define EL_SPARSEMATRIX_IMPL_HPP

#include <El/blas_like/level1/Axpy.hpp>
#include <El/blas_like/level1/Scale.hpp>

namespace El {

// Constructors and destructors
// ============================

template<typename Ring>
SparseMatrix<Ring>::SparseMatrix() { }

template<typename Ring>
SparseMatrix<Ring>::SparseMatrix( Int height, Int width )
: graph_(height,width)
{ }

template<typename Ring>
SparseMatrix<Ring>::SparseMatrix( const SparseMatrix<Ring>& A )
{
    EL_DEBUG_CSE
    if( &A != this )
        *this = A;
    else
        LogicError("Tried to construct sparse matrix with itself");
}

template<typename Ring>
SparseMatrix<Ring>::SparseMatrix( const DistSparseMatrix<Ring>& A )
{
    EL_DEBUG_CSE
    *this = A;
}

template<typename Ring>
SparseMatrix<Ring>::~SparseMatrix() { }

// Assignment and reconfiguration
// ==============================

// Change the size of the matrix
// -----------------------------
template<typename Ring>
void SparseMatrix<Ring>::Empty( bool clearMemory )
{
    graph_.Empty( clearMemory );
    if( clearMemory )
        SwapClear( vals_ );
    else
        vals_.resize( 0 );
}

template<typename Ring>
void SparseMatrix<Ring>::Resize( Int height, Int width )
{
    EL_DEBUG_CSE
    if( Height() == height && Width() == width )
        return;
    graph_.Resize( height, width );
    vals_.resize( 0 );
}

// Assembly
// --------
template<typename Ring>
void SparseMatrix<Ring>::Reserve( Int numEntries )
{
    const Int currSize = vals_.size();
    graph_.Reserve( numEntries );
    vals_.reserve( currSize+numEntries );
}

template<typename Ring>
void SparseMatrix<Ring>::FreezeSparsity() EL_NO_EXCEPT
{ graph_.frozenSparsity_ = true; }
template<typename Ring>
void SparseMatrix<Ring>::UnfreezeSparsity() EL_NO_EXCEPT
{ graph_.frozenSparsity_ = false; }
template<typename Ring>
bool SparseMatrix<Ring>::FrozenSparsity() const EL_NO_EXCEPT
{ return graph_.frozenSparsity_; }

template<typename Ring>
void SparseMatrix<Ring>::Update( Int row, Int col, const Ring& value )
{
    EL_DEBUG_CSE
    QueueUpdate( row, col, value );
    ProcessQueues();
}

template<typename Ring>
void SparseMatrix<Ring>::Update( const Entry<Ring>& entry )
{ Update( entry.i, entry.j, entry.value ); }

template<typename Ring>
void SparseMatrix<Ring>::Zero( Int row, Int col )
{
    EL_DEBUG_CSE
    QueueZero( row, col );
    ProcessQueues();
}

template<typename Ring>
void SparseMatrix<Ring>::QueueUpdate( Int row, Int col, const Ring& value )
EL_NO_RELEASE_EXCEPT
{
    EL_DEBUG_CSE
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

template<typename Ring>
void SparseMatrix<Ring>::QueueUpdate( const Entry<Ring>& entry )
EL_NO_RELEASE_EXCEPT
{ QueueUpdate( entry.i, entry.j, entry.value ); }

template<typename Ring>
void SparseMatrix<Ring>::QueueZero( Int row, Int col )
EL_NO_RELEASE_EXCEPT
{
    EL_DEBUG_CSE
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
template<typename Ring>
const SparseMatrix<Ring>& SparseMatrix<Ring>::operator=( const SparseMatrix<Ring>& A )
{
    EL_DEBUG_CSE
    graph_ = A.graph_;
    vals_ = A.vals_;
    return *this;
}

template<typename Ring>
const SparseMatrix<Ring>&
SparseMatrix<Ring>::operator=( const DistSparseMatrix<Ring>& A )
{
    EL_DEBUG_CSE
    if( A.Grid().Size() != 1 )
        LogicError("Can not yet construct from distributed sparse matrix");

    graph_ = A.distGraph_;
    vals_ = A.vals_;
    return *this;
}

// Make a copy of a submatrix
// --------------------------
template<typename Ring>
SparseMatrix<Ring>
SparseMatrix<Ring>::operator()( Range<Int> I, Range<Int> J ) const
{
    EL_DEBUG_CSE
    SparseMatrix<Ring> ASub;
    GetSubmatrix( *this, I, J, ASub );
    return ASub;
}

template<typename Ring>
SparseMatrix<Ring>
SparseMatrix<Ring>::operator()( const vector<Int>& I, Range<Int> J ) const
{
    EL_DEBUG_CSE
    SparseMatrix<Ring> ASub;
    GetSubmatrix( *this, I, J, ASub );
    return ASub;
}

template<typename Ring>
SparseMatrix<Ring>
SparseMatrix<Ring>::operator()( Range<Int> I, const vector<Int>& J ) const
{
    EL_DEBUG_CSE
    SparseMatrix<Ring> ASub;
    GetSubmatrix( *this, I, J, ASub );
    return ASub;
}

template<typename Ring>
SparseMatrix<Ring>
SparseMatrix<Ring>::operator()( const vector<Int>& I, const vector<Int>& J ) const
{
    EL_DEBUG_CSE
    SparseMatrix<Ring> ASub;
    GetSubmatrix( *this, I, J, ASub );
    return ASub;
}

// Rescaling
// ---------
template<typename Ring>
const SparseMatrix<Ring>& SparseMatrix<Ring>::operator*=( const Ring& alpha )
{
    EL_DEBUG_CSE
    Scale( alpha, *this );
    return *this;
}

// Addition/subtraction
// --------------------
template<typename Ring>
const SparseMatrix<Ring>& SparseMatrix<Ring>::operator+=
( const SparseMatrix<Ring>& A )
{
    EL_DEBUG_CSE
    Axpy( Ring(1), A, *this );
    return *this;
}

template<typename Ring>
const SparseMatrix<Ring>& SparseMatrix<Ring>::operator-=
( const SparseMatrix<Ring>& A )
{
    EL_DEBUG_CSE
    Axpy( Ring(-1), A, *this );
    return *this;
}

// Queries
// =======

// High-level information
// ----------------------
template<typename Ring>
Int SparseMatrix<Ring>::Height() const EL_NO_EXCEPT
{ return graph_.NumSources(); }
template<typename Ring>
Int SparseMatrix<Ring>::Width() const EL_NO_EXCEPT
{ return graph_.NumTargets(); }

template<typename Ring>
Int SparseMatrix<Ring>::NumEntries() const EL_NO_EXCEPT
{
    EL_DEBUG_CSE
    return graph_.NumEdges();
}

template<typename Ring>
Int SparseMatrix<Ring>::Capacity() const EL_NO_EXCEPT
{
    EL_DEBUG_CSE
    return graph_.Capacity();
}

template<typename Ring>
bool SparseMatrix<Ring>::Consistent() const EL_NO_EXCEPT
{ return graph_.Consistent(); }

template<typename Ring>
El::Graph& SparseMatrix<Ring>::Graph() EL_NO_EXCEPT
{ return graph_; }
template<typename Ring>
const El::Graph& SparseMatrix<Ring>::LockedGraph() const EL_NO_EXCEPT
{ return graph_; }

// Entrywise information
// ---------------------
template<typename Ring>
Int SparseMatrix<Ring>::Row( Int index ) const EL_NO_RELEASE_EXCEPT
{
    EL_DEBUG_CSE
    return graph_.Source( index );
}

template<typename Ring>
Int SparseMatrix<Ring>::Col( Int index ) const EL_NO_RELEASE_EXCEPT
{
    EL_DEBUG_CSE
    return graph_.Target( index );
}

template<typename Ring>
Int SparseMatrix<Ring>::RowOffset( Int row ) const EL_NO_RELEASE_EXCEPT
{
    EL_DEBUG_CSE
    return graph_.SourceOffset( row );
}

template<typename Ring>
Int SparseMatrix<Ring>::Offset( Int row, Int col ) const EL_NO_RELEASE_EXCEPT
{
    EL_DEBUG_CSE
    return graph_.Offset( row, col );
}

template<typename Ring>
Int SparseMatrix<Ring>::NumConnections( Int row ) const EL_NO_RELEASE_EXCEPT
{
    EL_DEBUG_CSE
    return graph_.NumConnections( row );
}

template<typename Ring>
Ring SparseMatrix<Ring>::Value( Int index ) const EL_NO_RELEASE_EXCEPT
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( index < 0 || index >= Int(vals_.size()) )
          LogicError("Entry number out of bounds");
    )
    return vals_[index];
}

template< typename Ring>
Ring SparseMatrix<Ring>::Get( Int row, Int col) const EL_NO_RELEASE_EXCEPT
{
    if( row == END ) row = graph_.numSources_ - 1;
    if( col == END ) col = graph_.numTargets_ - 1;
    Int index = Offset( row, col );
    if( Row(index) != row || Col(index) != col )
        return Ring(0);
    else
        return Value( index );
}

template< typename Ring>
void SparseMatrix<Ring>::Set
( Int row, Int col, const Ring& val) EL_NO_RELEASE_EXCEPT
{
    if( row == END ) row = graph_.numSources_ - 1;
    if( col == END ) col = graph_.numTargets_ - 1;
    Int index = Offset( row, col );
    if( Row(index) == row && Col(index) == col )
    {
        vals_[index] = val;
    }
    else
    {
        QueueUpdate( row, col, val );
        ProcessQueues();
    }
}

template<typename Ring>
Int* SparseMatrix<Ring>::SourceBuffer() EL_NO_EXCEPT
{ return graph_.SourceBuffer(); }
template<typename Ring>
Int* SparseMatrix<Ring>::TargetBuffer() EL_NO_EXCEPT
{ return graph_.TargetBuffer(); }
template<typename Ring>
Int* SparseMatrix<Ring>::OffsetBuffer() EL_NO_EXCEPT
{ return graph_.OffsetBuffer(); }
template<typename Ring>
Ring* SparseMatrix<Ring>::ValueBuffer() EL_NO_EXCEPT
{ return vals_.data(); }

template<typename Ring>
const Int* SparseMatrix<Ring>::LockedSourceBuffer() const EL_NO_EXCEPT
{ return graph_.LockedSourceBuffer(); }
template<typename Ring>
const Int* SparseMatrix<Ring>::LockedTargetBuffer() const EL_NO_EXCEPT
{ return graph_.LockedTargetBuffer(); }
template<typename Ring>
const Int* SparseMatrix<Ring>::LockedOffsetBuffer() const EL_NO_EXCEPT
{ return graph_.LockedOffsetBuffer(); }
template<typename Ring>
const Ring* SparseMatrix<Ring>::LockedValueBuffer() const EL_NO_EXCEPT
{ return vals_.data(); }

template<typename Ring>
void SparseMatrix<Ring>::ForceNumEntries( Int numEntries )
{
    EL_DEBUG_CSE
    graph_.ForceNumEdges( numEntries );
    vals_.resize( numEntries );
}

template<typename Ring>
void SparseMatrix<Ring>::ForceConsistency( bool consistent ) EL_NO_EXCEPT
{ graph_.ForceConsistency( consistent ); }

// Auxiliary routines
// ==================

template<typename Ring>
void SparseMatrix<Ring>::ProcessQueues()
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( graph_.sources_.size() != graph_.targets_.size() ||
          graph_.targets_.size() != vals_.size() )
          LogicError("Inconsistent sparse matrix buffer sizes");
    )
    if( graph_.consistent_ )
        return;

    Int numRemoved = 0;
    const Int numEntries = vals_.size();
    vector<Entry<Ring>> entries( numEntries );
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
              Entry<Ring>{graph_.sources_[s],graph_.targets_[s],vals_[s]};
    }
    CompareEntriesFunctor comparer;
    std::sort( entries.begin(), entries.end(), comparer );
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

template<typename Ring>
void SparseMatrix<Ring>::AssertConsistent() const
{ graph_.AssertConsistent(); }

#ifdef EL_INSTANTIATE_CORE
# define EL_EXTERN
#else
# define EL_EXTERN extern
#endif

#define PROTO(Ring) EL_EXTERN template class SparseMatrix<Ring>;
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

#undef EL_EXTERN

} // namespace El

#endif // ifndef EL_SPARSEMATRIX_IMPL_HPP
