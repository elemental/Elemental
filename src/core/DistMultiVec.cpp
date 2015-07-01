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
DistMultiVec<T>::DistMultiVec( mpi::Comm comm )
: height_(0), width_(0)
{ 
    if( comm == mpi::COMM_WORLD )
        comm_ = comm;
    else
        mpi::Dup( comm, comm_ );
    InitializeLocalData();
}

template<typename T>
DistMultiVec<T>::DistMultiVec( Int height, Int width, mpi::Comm comm )
: height_(height), width_(width)
{ 
    if( comm == mpi::COMM_WORLD )
        comm_ = comm;
    else
        mpi::Dup( comm, comm_ );
    InitializeLocalData();
}

template<typename T>
DistMultiVec<T>::DistMultiVec( const DistMultiVec<T>& A )
{
    DEBUG_ONLY(CSE cse("DistMultiVec::DistMultiVec"))
    height_ = 0;
    width_ = 0;
    comm_ = mpi::COMM_WORLD;
    if( &A != this )
        *this = A;
    else
        LogicError("Tried to construct DistMultiVec via itself");
}

template<typename T>
DistMultiVec<T>::~DistMultiVec()
{ 
    if( !mpi::Finalized() )
        if( comm_ != mpi::COMM_WORLD )
            mpi::Free( comm_ );
}

// Assignment and reconfiguration
// ==============================

// Change the size of the matrix
// -----------------------------
template<typename T>
void DistMultiVec<T>::Empty()
{
    height_ = 0;
    width_ = 0;
    blocksize_ = 0;
    firstLocalRow_ = 0;
    multiVec_.Empty();
}

template<typename T>
void DistMultiVec<T>::InitializeLocalData()
{
    const Int commRank = mpi::Rank( comm_ );
    const Int commSize = mpi::Size( comm_ );
    blocksize_ = height_/commSize;
    firstLocalRow_ = commRank*blocksize_;
    const Int localHeight =
        ( commRank<commSize-1 ?
          blocksize_ :
          height_ - (commSize-1)*blocksize_ );
    multiVec_.Resize( localHeight, width_ );
}

template<typename T>
void DistMultiVec<T>::Resize( Int height, Int width )
{
    height_ = height;
    width_ = width;
    InitializeLocalData();
}

// Change the distribution
// -----------------------
template<typename T>
void DistMultiVec<T>::SetComm( mpi::Comm comm )
{ 
    if( comm == comm_ )
        return;

    if( comm_ != mpi::COMM_WORLD )
        mpi::Free( comm_ );
    if( comm == mpi::COMM_WORLD )
        comm_ = comm;
    else
        mpi::Dup( comm, comm_ );

    InitializeLocalData();
}

// Operator overloading
// ====================

// Make a copy
// -----------
template<typename T>
const DistMultiVec<T>& DistMultiVec<T>::operator=( const DistMultiVec<T>& A )
{
    DEBUG_ONLY(CSE cse("DistMultiVec::operator="))
    Copy( A, *this );
    return *this;
}

template<typename T>
const DistMultiVec<T>& 
DistMultiVec<T>::operator=( const AbstractDistMatrix<T>& A )
{
    DEBUG_ONLY(CSE cse("DistMultiVec::operator="))
    Copy( A, *this );
    return *this;
}

// Make a copy of a submatrix
// --------------------------
template<typename T>
DistMultiVec<T>
DistMultiVec<T>::operator()( Range<Int> I, Range<Int> J ) const
{
    DEBUG_ONLY(CSE cse("DistMultiVec::operator()"))
    return GetSubmatrix( *this, I, J );
}   

// Rescaling
// ---------
template<typename T>
const DistMultiVec<T>& DistMultiVec<T>::operator*=( T alpha )
{
    DEBUG_ONLY(CSE cse("DistMultiVec::operator*=( T )"))
    Scale( alpha, *this );
    return *this;
}

// Addition/subtraction
// --------------------
template<typename T>
const DistMultiVec<T>& DistMultiVec<T>::operator+=( const DistMultiVec<T>& A )
{
    DEBUG_ONLY(CSE cse("DistMultiVec::operator+=( const DMV& )"))
    Axpy( T(1), A, *this );
    return *this;
}

template<typename T>
const DistMultiVec<T>& DistMultiVec<T>::operator-=( const DistMultiVec<T>& A )
{
    DEBUG_ONLY(CSE cse("DistMultiVec::operator-=( const DMV& )"))
    Axpy( T(-1), A, *this );
    return *this;
}

// Queries
// =======

// High-level information
// ----------------------
template<typename T>
Int DistMultiVec<T>::Height() const { return height_; }
template<typename T>
Int DistMultiVec<T>::Width() const { return multiVec_.Width(); }
template<typename T>
Int DistMultiVec<T>::FirstLocalRow() const { return firstLocalRow_; }
template<typename T>
Int DistMultiVec<T>::LocalHeight() const { return multiVec_.Height(); }
template<typename T>
El::Matrix<T>& DistMultiVec<T>::Matrix() { return multiVec_; }
template<typename T>
const El::Matrix<T>& DistMultiVec<T>::LockedMatrix() const { return multiVec_; }

// Distribution information
// ------------------------
template<typename T>
mpi::Comm DistMultiVec<T>::Comm() const { return comm_; }

template<typename T>
Int DistMultiVec<T>::Blocksize() const { return blocksize_; }

template<typename T>
int DistMultiVec<T>::RowOwner( Int i ) const 
{ 
    if( i == END ) i = height_ - 1;
    return RowToProcess( i, blocksize_, mpi::Size(comm_) ); 
}

template<typename T>
int DistMultiVec<T>::Owner( Int i, Int j ) const
{ return RowOwner(i); }

template<typename T>
bool DistMultiVec<T>::IsLocalRow( Int i ) const
{ return RowOwner(i) == mpi::Rank(comm_); }

template<typename T>
bool DistMultiVec<T>::IsLocal( Int i, Int j ) const
{ return IsLocalRow(i); }

template<typename T>
Int DistMultiVec<T>::GlobalRow( Int iLoc ) const
{
    DEBUG_ONLY(CSE cse("DistMultiVec::GlobalRow"))
    if( iLoc == END ) iLoc = LocalHeight() - 1;
    if( iLoc < 0 || iLoc >= LocalHeight() )
        LogicError("Invalid local row index");
    return iLoc + FirstLocalRow();
}

template<typename T>
Int DistMultiVec<T>::LocalRow( Int i ) const
{
    DEBUG_ONLY(CSE cse("DistMultiVec::LocalRow")) 
    if( i == END ) i = Height() - 1;
    if( i < 0 || i >= Height() )
        LogicError("Invalid global row index");
    return i - FirstLocalRow();
}

// Detailed local information
// --------------------------
template<typename T>
T DistMultiVec<T>::Get( Int i, Int j ) const
{
    DEBUG_ONLY(CSE cse("DistMultiVec::Get"))
    if( i == END ) i = height_ - 1;
    const int rowOwner = RowOwner(i);
    T value;
    if( rowOwner == mpi::Rank(comm_) )
        value = GetLocal( i-FirstLocalRow(), j );
    mpi::Broadcast( value, rowOwner, comm_ );
    return value;
}

template<typename T>
void DistMultiVec<T>::Set( Int i, Int j, T value )
{
    DEBUG_ONLY(CSE cse("DistMultiVec::Set"))
    if( i == END ) i = height_ - 1;
    const Int firstLocalRow = FirstLocalRow();
    if( i >= firstLocalRow && i < firstLocalRow+LocalHeight() )
        SetLocal( i-firstLocalRow, j, value );
}

template<typename T>
void DistMultiVec<T>::Set( const Entry<T>& entry )
{ Set( entry.i, entry.j, entry.value ); }

template<typename T>
void DistMultiVec<T>::Update( Int i, Int j, T value )
{
    DEBUG_ONLY(CSE cse("DistMultiVec::Update"))
    if( i == END ) i = height_ - 1;
    const Int firstLocalRow = FirstLocalRow();
    if( i >= firstLocalRow && i < firstLocalRow+LocalHeight() )
        UpdateLocal( i-firstLocalRow, j, value );
}

template<typename T>
void DistMultiVec<T>::Update( const Entry<T>& entry )
{ Update( entry.i, entry.j, entry.value ); }

template<typename T>
T DistMultiVec<T>::GetLocal( Int iLoc, Int j ) const
{ 
    DEBUG_ONLY(CSE cse("DistMultiVec::GetLocal"))
    return multiVec_.Get(iLoc,j);
}

template<typename T>
void DistMultiVec<T>::SetLocal( Int iLoc, Int j, T value )
{
    DEBUG_ONLY(CSE cse("DistMultiVec::SetLocal"))
    multiVec_.Set(iLoc,j,value);
}

template<typename T>
void DistMultiVec<T>::SetLocal( const Entry<T>& localEntry )
{ SetLocal( localEntry.i, localEntry.j, localEntry.value ); }

template<typename T>
void DistMultiVec<T>::UpdateLocal( Int iLoc, Int j, T value )
{
    DEBUG_ONLY(CSE cse("DistMultiVec::UpdateLocal"))
    multiVec_.Update(iLoc,j,value);
}

template<typename T>
void DistMultiVec<T>::UpdateLocal( const Entry<T>& localEntry )
{ UpdateLocal( localEntry.i, localEntry.j, localEntry.value ); }

// Batch remote updates
// --------------------
template<typename T>
void DistMultiVec<T>::Reserve( Int numRemoteEntries )
{
    DEBUG_ONLY(CSE cse("DistMultiVec::Reserve"))
    remoteUpdates_.reserve( numRemoteEntries );
}

template<typename T>
void DistMultiVec<T>::QueueUpdate( const Entry<T>& entry )
{
    DEBUG_ONLY(CSE cse("DistMultiVec::QueueUpdate"))
    remoteUpdates_.push_back( entry );
}

template<typename T>
void DistMultiVec<T>::QueueUpdate( Int i, Int j, T value )
{ QueueUpdate( Entry<T>{i,j,value} ); }

template<typename T>
void DistMultiVec<T>::ProcessQueues()
{
    DEBUG_ONLY(CSE cse("DistMultiVec::ProcessQueues"))
    int commSize = mpi::Size( comm_ );
    // Compute the send counts
    // -----------------------
    vector<int> sendCounts(commSize);
    for( auto entry : remoteUpdates_ )
        ++sendCounts[Owner(entry.i,entry.j)];
    // Pack the send data
    // ------------------
    vector<int> sendOffs;
    const int totalSend = Scan( sendCounts, sendOffs );
    auto offs = sendOffs;
    vector<Entry<T>> sendEntries(totalSend);
    for( auto entry : remoteUpdates_ )
        sendEntries[offs[Owner(entry.i,entry.j)]++] = entry;
    SwapClear( remoteUpdates_ );
    // Exchange and unpack
    // -------------------
    auto recvEntries = 
      mpi::AllToAll( sendEntries, sendCounts, sendOffs, comm_ );
    for( auto entry : recvEntries )
        Update( entry );
}

#define PROTO(T) template class DistMultiVec<T>;

#define EL_ENABLE_QUAD
#include "El/macros/Instantiate.h"

} // namespace El
