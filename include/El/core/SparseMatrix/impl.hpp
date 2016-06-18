/*
   Copyright (c) 2009-2016, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, Stanford University, and the
   Georgia Insitute of Technology.
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

template<typename T>
SparseMatrix<T>::SparseMatrix() { }

template<typename T>
SparseMatrix<T>::SparseMatrix( Int height, Int width )
: graph_(height,width)
{ }

template<typename T>
SparseMatrix<T>::SparseMatrix( const SparseMatrix<T>& A )
{ 
    DEBUG_CSE
    if( &A != this )
        *this = A;
    else
        LogicError("Tried to construct sparse matrix with itself");
}

template<typename T>
SparseMatrix<T>::SparseMatrix( const DistSparseMatrix<T>& A )
{ 
    DEBUG_CSE
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
    DEBUG_CSE
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
    const Int currSize = vals_.size();
    graph_.Reserve( numEntries );
    vals_.reserve( currSize+numEntries );
}

template<typename T>
void SparseMatrix<T>::FreezeSparsity() EL_NO_EXCEPT
{ graph_.frozenSparsity_ = true; }
template<typename T>
void SparseMatrix<T>::UnfreezeSparsity() EL_NO_EXCEPT
{ graph_.frozenSparsity_ = false; }
template<typename T>
bool SparseMatrix<T>::FrozenSparsity() const EL_NO_EXCEPT
{ return graph_.frozenSparsity_; }

template<typename T>
void SparseMatrix<T>::Update( Int row, Int col, T value )
{
    DEBUG_CSE
    QueueUpdate( row, col, value );
    ProcessQueues();
}

template<typename T>
void SparseMatrix<T>::Update( const Entry<T>& entry )
{ Update( entry.i, entry.j, entry.value ); }

template<typename T>
void SparseMatrix<T>::Zero( Int row, Int col )
{
    DEBUG_CSE
    QueueZero( row, col );
    ProcessQueues();
}

template<typename T>
void SparseMatrix<T>::QueueUpdate( Int row, Int col, T value )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
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
EL_NO_RELEASE_EXCEPT
{ QueueUpdate( entry.i, entry.j, entry.value ); }

template<typename T>
void SparseMatrix<T>::QueueZero( Int row, Int col )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
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
    DEBUG_CSE
    graph_ = A.graph_;
    vals_ = A.vals_;
    return *this;
}

template<typename T>
const SparseMatrix<T>&
SparseMatrix<T>::operator=( const DistSparseMatrix<T>& A )
{
    DEBUG_CSE
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
    DEBUG_CSE
    SparseMatrix<T> ASub;
    GetSubmatrix( *this, I, J, ASub );
    return ASub;
}

template<typename T>
SparseMatrix<T>
SparseMatrix<T>::operator()( const vector<Int>& I, Range<Int> J ) const
{
    DEBUG_CSE
    SparseMatrix<T> ASub;
    GetSubmatrix( *this, I, J, ASub );
    return ASub;
}

template<typename T>
SparseMatrix<T>
SparseMatrix<T>::operator()( Range<Int> I, const vector<Int>& J ) const
{
    DEBUG_CSE
    SparseMatrix<T> ASub;
    GetSubmatrix( *this, I, J, ASub );
    return ASub;
}

template<typename T>
SparseMatrix<T>
SparseMatrix<T>::operator()( const vector<Int>& I, const vector<Int>& J ) const
{
    DEBUG_CSE
    SparseMatrix<T> ASub;
    GetSubmatrix( *this, I, J, ASub );
    return ASub;
}

// Rescaling
// ---------
template<typename T>
const SparseMatrix<T>& SparseMatrix<T>::operator*=( T alpha )
{
    DEBUG_CSE
    Scale( alpha, *this );
    return *this;
}

// Addition/subtraction
// --------------------
template<typename T>
const SparseMatrix<T>& SparseMatrix<T>::operator+=( const SparseMatrix<T>& A )
{
    DEBUG_CSE
    Axpy( T(1), A, *this );
    return *this;
}

template<typename T>
const SparseMatrix<T>& SparseMatrix<T>::operator-=( const SparseMatrix<T>& A )
{
    DEBUG_CSE
    Axpy( T(-1), A, *this );
    return *this;
}

// Queries
// =======

// High-level information
// ----------------------
template<typename T>
Int SparseMatrix<T>::Height() const EL_NO_EXCEPT
{ return graph_.NumSources(); }
template<typename T>
Int SparseMatrix<T>::Width() const EL_NO_EXCEPT
{ return graph_.NumTargets(); }

template<typename T>
Int SparseMatrix<T>::NumEntries() const EL_NO_EXCEPT
{
    DEBUG_CSE
    return graph_.NumEdges();
}

template<typename T>
Int SparseMatrix<T>::Capacity() const EL_NO_EXCEPT
{
    DEBUG_CSE
    return graph_.Capacity();
}

template<typename T>
bool SparseMatrix<T>::Consistent() const EL_NO_EXCEPT
{ return graph_.Consistent(); }

template<typename T>
El::Graph& SparseMatrix<T>::Graph() EL_NO_EXCEPT
{ return graph_; }
template<typename T>
const El::Graph& SparseMatrix<T>::LockedGraph() const EL_NO_EXCEPT
{ return graph_; }

// Entrywise information
// ---------------------
template<typename T>
Int SparseMatrix<T>::Row( Int index ) const EL_NO_RELEASE_EXCEPT
{ 
    DEBUG_CSE
    return graph_.Source( index );
}

template<typename T>
Int SparseMatrix<T>::Col( Int index ) const EL_NO_RELEASE_EXCEPT
{ 
    DEBUG_CSE
    return graph_.Target( index );
}

template<typename T>
Int SparseMatrix<T>::RowOffset( Int row ) const EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    return graph_.SourceOffset( row );
}

template<typename T>
Int SparseMatrix<T>::Offset( Int row, Int col ) const EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    return graph_.Offset( row, col );
}

template<typename T>
Int SparseMatrix<T>::NumConnections( Int row ) const EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    return graph_.NumConnections( row );
}

template<typename T>
T SparseMatrix<T>::Value( Int index ) const EL_NO_RELEASE_EXCEPT
{ 
    DEBUG_CSE
    DEBUG_ONLY(
      if( index < 0 || index >= Int(vals_.size()) )
          LogicError("Entry number out of bounds");
    )
    return vals_[index];
}

template< typename T>
T SparseMatrix<T>::Get( Int row, Int col) const EL_NO_RELEASE_EXCEPT
{
    if( row == END ) row = graph_.numSources_ - 1;
    if( col == END ) col = graph_.numTargets_ - 1; 
    Int index = Offset( row, col );
    if( Row(index) != row || Col(index) != col )
        return T(0); 
    else
        return Value( index );
}

template< typename T>
void SparseMatrix<T>::Set( Int row, Int col, T val) EL_NO_RELEASE_EXCEPT
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

template<typename T>
Int* SparseMatrix<T>::SourceBuffer() EL_NO_EXCEPT
{ return graph_.SourceBuffer(); }
template<typename T>
Int* SparseMatrix<T>::TargetBuffer() EL_NO_EXCEPT
{ return graph_.TargetBuffer(); }
template<typename T>
Int* SparseMatrix<T>::OffsetBuffer() EL_NO_EXCEPT
{ return graph_.OffsetBuffer(); }
template<typename T>
T* SparseMatrix<T>::ValueBuffer() EL_NO_EXCEPT
{ return vals_.data(); }

template<typename T>
const Int* SparseMatrix<T>::LockedSourceBuffer() const EL_NO_EXCEPT
{ return graph_.LockedSourceBuffer(); }
template<typename T>
const Int* SparseMatrix<T>::LockedTargetBuffer() const EL_NO_EXCEPT
{ return graph_.LockedTargetBuffer(); }
template<typename T>
const Int* SparseMatrix<T>::LockedOffsetBuffer() const EL_NO_EXCEPT
{ return graph_.LockedOffsetBuffer(); }
template<typename T>
const T* SparseMatrix<T>::LockedValueBuffer() const EL_NO_EXCEPT
{ return vals_.data(); }

template<typename T>
void SparseMatrix<T>::ForceNumEntries( Int numEntries )
{
    DEBUG_CSE
    graph_.ForceNumEdges( numEntries );
    vals_.resize( numEntries );
}

template<typename T>
void SparseMatrix<T>::ForceConsistency( bool consistent ) EL_NO_EXCEPT
{ graph_.ForceConsistency( consistent ); }

// Auxiliary routines
// ==================

template<typename T>
void SparseMatrix<T>::ProcessQueues()
{
    DEBUG_CSE
    DEBUG_ONLY(
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

template<typename T>
void SparseMatrix<T>::AssertConsistent() const
{ graph_.AssertConsistent(); }

#ifdef EL_INSTANTIATE_CORE
# define EL_EXTERN
#else
# define EL_EXTERN extern
#endif

#define PROTO(T) EL_EXTERN template class SparseMatrix<T>;
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

#undef EL_EXTERN

} // namespace El

#endif // ifndef EL_SPARSEMATRIX_IMPL_HPP
