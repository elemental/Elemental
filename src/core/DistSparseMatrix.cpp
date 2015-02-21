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
DistSparseMatrix<T>::DistSparseMatrix()
{ }

template<typename T>
DistSparseMatrix<T>::DistSparseMatrix( mpi::Comm comm )
: distGraph_(comm)
{ }

template<typename T>
DistSparseMatrix<T>::DistSparseMatrix( Int height, mpi::Comm comm )
: distGraph_(height,comm)
{ }

template<typename T>
DistSparseMatrix<T>::DistSparseMatrix( Int height, Int width, mpi::Comm comm )
: distGraph_(height,width,comm)
{ }

template<typename T>
DistSparseMatrix<T>::~DistSparseMatrix()
{ }

// Assignment and reconfiguration
// ==============================

// Make a copy
// -----------
// TODO

// Change the matrix size
// ----------------------
template<typename T>
void DistSparseMatrix<T>::Empty( bool clearMemory )
{
    distGraph_.Empty( clearMemory );
    if( clearMemory )
        SwapClear( vals_ );
    else
        vals_.resize( 0 );
    multMeta.Clear();
}

template<typename T>
void DistSparseMatrix<T>::Resize( Int height, Int width )
{
    distGraph_.Resize( height, width );
    vals_.resize( 0 );
}

// Change the distribution
// -----------------------
template<typename T>
void DistSparseMatrix<T>::SetComm( mpi::Comm comm )
{ 
    if( Comm() == comm )
        return;
    distGraph_.SetComm( comm ); 
    vals_.resize( 0 );
}

// Assembly
// --------
template<typename T>
void DistSparseMatrix<T>::Reserve( Int numLocalEntries )
{ 
    distGraph_.Reserve( numLocalEntries );
    vals_.reserve( numLocalEntries );
}

template<typename T>
void DistSparseMatrix<T>::Update( Int row, Int col, T value )
{
    DEBUG_ONLY(CallStackEntry cse("DistSparseMatrix::Update"))
    QueueUpdate( row, col, value );
    MakeConsistent();
}

template<typename T>
void DistSparseMatrix<T>::UpdateLocal( Int localRow, Int col, T value )
{
    DEBUG_ONLY(CallStackEntry cse("DistSparseMatrix::UpdateLocal"))
    QueueLocalUpdate( localRow, col, value );
    MakeConsistent();
}

template<typename T>
void DistSparseMatrix<T>::Zero( Int row, Int col )
{
    DEBUG_ONLY(CallStackEntry cse("DistSparseMatrix::Zero"))
    QueueZero( row, col );
    MakeConsistent();
}

template<typename T>
void DistSparseMatrix<T>::ZeroLocal( Int localRow, Int col )
{
    DEBUG_ONLY(CallStackEntry cse("DistSparseMatrix::ZeroLocal"))
    QueueLocalZero( localRow, col );
    MakeConsistent();
}

template<typename T>
void DistSparseMatrix<T>::QueueUpdate( Int row, Int col, T value )
{
    DEBUG_ONLY(CallStackEntry cse("DistSparseMatrix::QueueUpdate"))
    if( row >= FirstLocalRow() && row < FirstLocalRow()+LocalHeight() )
        QueueLocalUpdate( row-FirstLocalRow(), col, value );
}

template<typename T>
void DistSparseMatrix<T>::QueueLocalUpdate( Int localRow, Int col, T value )
{
    DEBUG_ONLY(CallStackEntry cse("DistSparseMatrix::QueueLocalUpdate"))
    distGraph_.QueueLocalConnection( localRow, col );
    vals_.push_back( value );
    multMeta.ready = false;
}

template<typename T>
void DistSparseMatrix<T>::QueueZero( Int row, Int col )
{
    DEBUG_ONLY(CallStackEntry cse("DistSparseMatrix::QueueZero"))
    if( row >= FirstLocalRow() && row < FirstLocalRow()+LocalHeight() )
        QueueLocalZero( row-FirstLocalRow(), col );
}

template<typename T>
void DistSparseMatrix<T>::QueueLocalZero( Int localRow, Int col )
{
    DEBUG_ONLY(CallStackEntry cse("DistSparseMatrix::QueueZero"))
    distGraph_.QueueLocalDisconnection( localRow, col );
    multMeta.ready = false;
}

template<typename T>
void DistSparseMatrix<T>::MakeConsistent()
{
    DEBUG_ONLY(
      CallStackEntry cse("DistSparseMatrix::MakeConsistent");
      if( distGraph_.sources_.size() != distGraph_.targets_.size() || 
          distGraph_.targets_.size() != vals_.size() )
          LogicError("Inconsistent sparse matrix buffer sizes");
    )

    if( !distGraph_.consistent_ )
    {
        const Int numLocalEntries = vals_.size();
        Int numRemoved = 0;
        vector<Entry<T>> entries( numLocalEntries );
        for( Int s=0; s<numLocalEntries; ++s )
        {
            pair<Int,Int> 
              candidate(distGraph_.sources_[s],distGraph_.targets_[s]);
            if( distGraph_.markedForRemoval_.find(candidate) ==
                distGraph_.markedForRemoval_.end() )
            {
                entries[s-numRemoved].indices[0] = distGraph_.sources_[s];
                entries[s-numRemoved].indices[1] = distGraph_.targets_[s];
                entries[s-numRemoved].value = vals_[s];
            }
            else
            {
                ++numRemoved;
            }
        }
        distGraph_.markedForRemoval_.clear();
        entries.resize( numLocalEntries-numRemoved );
        std::sort( entries.begin(), entries.end(), CompareEntries );

        // Compress out duplicates
        Int lastUnique=0;
        for( Int s=1; s<numLocalEntries; ++s )
        {
            if( entries[s].indices[0] != entries[lastUnique].indices[0] ||
                entries[s].indices[1] != entries[lastUnique].indices[1] )
            {
                ++lastUnique;
                if( s != lastUnique )
                {
                    entries[lastUnique].indices[0] = entries[s].indices[0];
                    entries[lastUnique].indices[1] = entries[s].indices[1];
                    entries[lastUnique].value = entries[s].value;
                }
            }
            else
                entries[lastUnique].value += entries[s].value;
        }
        const Int numUnique = lastUnique+1;
        entries.resize( numUnique );

        distGraph_.sources_.resize( numUnique );
        distGraph_.targets_.resize( numUnique );
        vals_.resize( numUnique );
        for( Int s=0; s<numUnique; ++s )
        {
            distGraph_.sources_[s] = entries[s].indices[0];
            distGraph_.targets_[s] = entries[s].indices[1];
            vals_[s] = entries[s].value;
        }

        distGraph_.ComputeEdgeOffsets();

        distGraph_.consistent_ = true;
    }
}

// Queries
// =======

// High-level information
// ----------------------
template<typename T>
Int DistSparseMatrix<T>::Height() const { return distGraph_.NumSources(); }
template<typename T>
Int DistSparseMatrix<T>::Width() const { return distGraph_.NumTargets(); }

template<typename T>
El::DistGraph& DistSparseMatrix<T>::DistGraph() { return distGraph_; }
template<typename T>
const El::DistGraph& DistSparseMatrix<T>::LockedDistGraph() const
{ return distGraph_; }

template<typename T>
Int DistSparseMatrix<T>::FirstLocalRow() const
{ return distGraph_.FirstLocalSource(); }

template<typename T>
Int DistSparseMatrix<T>::LocalHeight() const
{ return distGraph_.NumLocalSources(); }

template<typename T>
Int DistSparseMatrix<T>::NumLocalEntries() const
{
    DEBUG_ONLY(CallStackEntry cse("DistSparseMatrix::NumLocalEntries"))
    return distGraph_.NumLocalEdges();
}

template<typename T>
Int DistSparseMatrix<T>::Capacity() const
{
    DEBUG_ONLY(CallStackEntry cse("DistSparseMatrix::Capacity"))
    return distGraph_.Capacity();
}

template<typename T>
bool DistSparseMatrix<T>::Consistent() const
{ return distGraph_.Consistent(); }

// Distribution information
// ------------------------
template<typename T>
mpi::Comm DistSparseMatrix<T>::Comm() const { return distGraph_.Comm(); }
template<typename T>
Int DistSparseMatrix<T>::Blocksize() const { return distGraph_.Blocksize(); }
template<typename T>
int DistSparseMatrix<T>::RowOwner( Int i ) const 
{ return RowToProcess( i, Blocksize(), mpi::Size(Comm()) ); }

template<typename T>
Int DistSparseMatrix<T>::GlobalRow( Int iLoc ) const
{ 
    DEBUG_ONLY(CallStackEntry cse("DistSparseMatrix::GlobalRow"))
    if( iLoc < 0 || iLoc > LocalHeight() )
        LogicError("Invalid local row index");
    return iLoc + FirstLocalRow(); 
}

// Detailed local information
// --------------------------
template<typename T>
Int DistSparseMatrix<T>::Row( Int localInd ) const
{ 
    DEBUG_ONLY(CallStackEntry cse("DistSparseMatrix::Row"))
    return distGraph_.Source( localInd );
}

template<typename T>
Int DistSparseMatrix<T>::Col( Int localInd ) const
{ 
    DEBUG_ONLY(CallStackEntry cse("DistSparseMatrix::Col"))
    return distGraph_.Target( localInd );
}

template<typename T>
Int DistSparseMatrix<T>::EntryOffset( Int localRow ) const
{
    DEBUG_ONLY(CallStackEntry cse("DistSparseMatrix::EntryOffset"))
    return distGraph_.EdgeOffset( localRow );
}

template<typename T>
Int DistSparseMatrix<T>::NumConnections( Int localRow ) const
{
    DEBUG_ONLY(CallStackEntry cse("DistSparseMatrix::NumConnections"))
    return distGraph_.NumConnections( localRow );
}

template<typename T>
T DistSparseMatrix<T>::Value( Int localInd ) const
{ 
    DEBUG_ONLY(
        CallStackEntry cse("DistSparseMatrix::Value");
        if( localInd < 0 || localInd >= (Int)vals_.size() )
            LogicError("Entry number out of bounds");
        AssertConsistent();
    )
    return vals_[localInd];
}

template<typename T>
Int* DistSparseMatrix<T>::SourceBuffer() { return distGraph_.SourceBuffer(); }
template<typename T>
Int* DistSparseMatrix<T>::TargetBuffer() { return distGraph_.TargetBuffer(); }
template<typename T>
T* DistSparseMatrix<T>::ValueBuffer() { return vals_.data(); }

template<typename T>
const Int* DistSparseMatrix<T>::LockedSourceBuffer() const
{ return distGraph_.LockedSourceBuffer(); }

template<typename T>
const Int* DistSparseMatrix<T>::LockedTargetBuffer() const
{ return distGraph_.LockedTargetBuffer(); }

template<typename T>
const T* DistSparseMatrix<T>::LockedValueBuffer() const
{ return vals_.data(); }

// Auxiliary routines
// ==================

template<typename T>
bool DistSparseMatrix<T>::CompareEntries( const Entry<T>& a, const Entry<T>& b )
{ return a.indices[0] < b.indices[0] || 
         (a.indices[0] == b.indices[0] && a.indices[1] < b.indices[1]); }

template<typename T>
void DistSparseMatrix<T>::AssertConsistent() const
{ 
    if( !Consistent() )
        LogicError("Distributed sparse matrix must be consistent");
}

#define PROTO(T) template class DistSparseMatrix<T>;
#include "El/macros/Instantiate.h"

} // namespace El
