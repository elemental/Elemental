/*
   Copyright (c) 2009-2016, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, Stanford University, and the
   Georgia Insitute of Technology.
   All rights reserved.
 
   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

#include <El/blas_like/level1/Axpy.hpp>
#include <El/blas_like/level1/Scale.hpp>

namespace El {

// Constructors and destructors
// ============================

template<typename T>
DistMultiVec<T>::DistMultiVec( mpi::Comm comm )
: height_(0), width_(0),
  commSize_(mpi::Size(comm)), commRank_(mpi::Rank(comm))
{ 
    if( comm == mpi::COMM_WORLD )
        comm_ = comm;
    else
        mpi::Dup( comm, comm_ );
    InitializeLocalData();
}

template<typename T>
DistMultiVec<T>::DistMultiVec( Int height, Int width, mpi::Comm comm )
: height_(height), width_(width), 
  commSize_(mpi::Size(comm)), commRank_(mpi::Rank(comm))
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
    DEBUG_CSE
    height_ = 0;
    width_ = 0;
    comm_ = mpi::COMM_WORLD;
    if( &A != this )
        *this = A;
    DEBUG_ONLY(
      else
          LogicError("Tried to construct DistMultiVec via itself");
    )
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
void DistMultiVec<T>::Empty( bool freeMemory )
{
    height_ = 0;
    width_ = 0;
    blocksize_ = 1;

    multiVec_.Empty( freeMemory );

    SwapClear( remoteUpdates_ );
}

template<typename T>
void DistMultiVec<T>::InitializeLocalData()
{
    blocksize_ = height_ / commSize_;
    if( blocksize_*commSize_ < height_ || height_ == 0 )
        ++blocksize_;
    const Int localHeight = Min(blocksize_,Max(0,height_-blocksize_*commRank_));
    multiVec_.Resize( localHeight, width_ );
}

template<typename T>
void DistMultiVec<T>::Resize( Int height, Int width )
{
    if( height_ == height && width == width_ )
        return;

    height_ = height;
    width_ = width;
    InitializeLocalData();

    SwapClear( remoteUpdates_ );
}

// Change the distribution
// -----------------------
template<typename T>
void DistMultiVec<T>::SetComm( mpi::Comm comm )
{ 
    commSize_ = mpi::Size(comm);
    commRank_ = mpi::Rank(comm);
    if( comm == comm_ )
        return;

    if( comm_ != mpi::COMM_WORLD )
        mpi::Free( comm_ );
    if( comm == mpi::COMM_WORLD )
        comm_ = comm;
    else
        mpi::Dup( comm, comm_ );

    Resize( 0, 0 );
}

// Operator overloading
// ====================

// Make a copy
// -----------
template<typename T>
const DistMultiVec<T>& DistMultiVec<T>::operator=( const DistMultiVec<T>& A )
{
    DEBUG_CSE
    Copy( A, *this );
    return *this;
}

template<typename T>
const DistMultiVec<T>& 
DistMultiVec<T>::operator=( const AbstractDistMatrix<T>& A )
{
    DEBUG_CSE
    Copy( A, *this );
    return *this;
}

// Make a copy of a submatrix
// --------------------------
template<typename T>
DistMultiVec<T>
DistMultiVec<T>::operator()( Range<Int> I, Range<Int> J ) const
{
    DEBUG_CSE
    DistMultiVec<T> ASub(this->Comm());
    GetSubmatrix( *this, I, J, ASub );
    return ASub;
}

template<typename T>
DistMultiVec<T>
DistMultiVec<T>::operator()( Range<Int> I, const vector<Int>& J ) const
{
    DEBUG_CSE
    DistMultiVec<T> ASub(this->Comm());
    GetSubmatrix( *this, I, J, ASub );
    return ASub;
}

template<typename T>
DistMultiVec<T>
DistMultiVec<T>::operator()( const vector<Int>& I, Range<Int> J ) const
{
    DEBUG_CSE
    DistMultiVec<T> ASub(this->Comm());
    GetSubmatrix( *this, I, J, ASub );
    return ASub;
}

template<typename T>
DistMultiVec<T>
DistMultiVec<T>::operator()( const vector<Int>& I, const vector<Int>& J ) const
{
    DEBUG_CSE
    DistMultiVec<T> ASub(this->Comm());
    GetSubmatrix( *this, I, J, ASub );
    return ASub;
}

// Rescaling
// ---------
template<typename T>
const DistMultiVec<T>& DistMultiVec<T>::operator*=( T alpha )
{
    DEBUG_CSE
    Scale( alpha, *this );
    return *this;
}

// Addition/subtraction
// --------------------
template<typename T>
const DistMultiVec<T>& DistMultiVec<T>::operator+=( const DistMultiVec<T>& A )
{
    DEBUG_CSE
    Axpy( T(1), A, *this );
    return *this;
}

template<typename T>
const DistMultiVec<T>& DistMultiVec<T>::operator-=( const DistMultiVec<T>& A )
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
Int DistMultiVec<T>::Height() const EL_NO_EXCEPT
{ return height_; }
template<typename T>
Int DistMultiVec<T>::Width() const EL_NO_EXCEPT
{ return multiVec_.Width(); }
template<typename T>
Int DistMultiVec<T>::FirstLocalRow() const EL_NO_EXCEPT
{ return blocksize_*commRank_; }
template<typename T>
Int DistMultiVec<T>::LocalHeight() const EL_NO_EXCEPT
{ return multiVec_.Height(); }
template<typename T>
El::Matrix<T>& DistMultiVec<T>::Matrix() EL_NO_EXCEPT
{ return multiVec_; }
template<typename T>
const El::Matrix<T>& DistMultiVec<T>::LockedMatrix() const EL_NO_EXCEPT
{ return multiVec_; }

// Distribution information
// ------------------------
template<typename T>
mpi::Comm DistMultiVec<T>::Comm() const EL_NO_EXCEPT { return comm_; }

template<typename T>
Int DistMultiVec<T>::Blocksize() const EL_NO_EXCEPT { return blocksize_; }

template<typename T>
int DistMultiVec<T>::RowOwner( Int i ) const EL_NO_EXCEPT
{ 
    if( i == END ) i = height_ - 1;
    return i / blocksize_;
}

template<typename T>
int DistMultiVec<T>::Owner( Int i, Int j ) const EL_NO_EXCEPT
{ return RowOwner(i); }

template<typename T>
bool DistMultiVec<T>::IsLocalRow( Int i ) const EL_NO_EXCEPT
{ return RowOwner(i) == commRank_; }

template<typename T>
bool DistMultiVec<T>::IsLocal( Int i, Int j ) const EL_NO_EXCEPT
{ return IsLocalRow(i); }

template<typename T>
Int DistMultiVec<T>::GlobalRow( Int iLoc ) const
{
    DEBUG_CSE
    if( iLoc == END ) iLoc = LocalHeight() - 1;
    DEBUG_ONLY(
      if( iLoc < 0 || iLoc >= LocalHeight() )
          LogicError("Invalid local row index");
    )
    return iLoc + FirstLocalRow();
}

template<typename T>
Int DistMultiVec<T>::LocalRow( Int i ) const
{
    DEBUG_CSE
    if( i == END ) i = Height() - 1;
    DEBUG_ONLY(
      if( i < 0 || i >= Height() )
          LogicError("Invalid global row index");
    )
    return i - FirstLocalRow();
}

// Detailed local information
// --------------------------
template<typename T>
T DistMultiVec<T>::Get( Int i, Int j ) const
{
    DEBUG_CSE
    if( i == END ) i = height_ - 1;
    const int rowOwner = RowOwner(i);
    T value;
    if( rowOwner == commRank_ )
        value = GetLocal( i-FirstLocalRow(), j );
    mpi::Broadcast( value, rowOwner, comm_ );
    return value;
}

template<typename T>
void DistMultiVec<T>::Set( Int i, Int j, T value )
{
    DEBUG_CSE
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
    DEBUG_CSE
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
    DEBUG_CSE
    return multiVec_.Get(iLoc,j);
}

template<typename T>
void DistMultiVec<T>::SetLocal( Int iLoc, Int j, T value )
{
    DEBUG_CSE
    multiVec_.Set(iLoc,j,value);
}

template<typename T>
void DistMultiVec<T>::SetLocal( const Entry<T>& localEntry )
{ SetLocal( localEntry.i, localEntry.j, localEntry.value ); }

template<typename T>
void DistMultiVec<T>::UpdateLocal( Int iLoc, Int j, T value )
{
    DEBUG_CSE
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
    DEBUG_CSE
    const Int currSize = remoteUpdates_.size();
    remoteUpdates_.reserve( currSize+numRemoteEntries );
}

template<typename T>
void DistMultiVec<T>::QueueUpdate( const Entry<T>& entry )
{
    DEBUG_CSE
    remoteUpdates_.push_back( entry );
}

template<typename T>
void DistMultiVec<T>::QueueUpdate( Int i, Int j, T value )
{ QueueUpdate( Entry<T>{i,j,value} ); }

template<typename T>
void DistMultiVec<T>::ProcessQueues()
{
    DEBUG_CSE
    // Compute the send counts
    // -----------------------
    vector<int> sendCounts(commSize_);
    for( const auto& entry : remoteUpdates_ )
        ++sendCounts[Owner(entry.i,entry.j)];
    // Pack the send data
    // ------------------
    vector<int> sendOffs;
    const int totalSend = Scan( sendCounts, sendOffs );
    auto offs = sendOffs;
    vector<Entry<T>> sendEntries(totalSend);
    for( const auto& entry : remoteUpdates_ )
        sendEntries[offs[Owner(entry.i,entry.j)]++] = entry;
    SwapClear( remoteUpdates_ );
    // Exchange and unpack
    // -------------------
    auto recvEntries = 
      mpi::AllToAll( sendEntries, sendCounts, sendOffs, comm_ );

    T* matBuf = multiVec_.Buffer();
    const Int matLDim = multiVec_.LDim();
    const Int firstLocalRow = FirstLocalRow();
    for( const auto& entry : recvEntries )
        matBuf[(entry.i-firstLocalRow)+entry.j*matLDim] += entry.value;
}

#ifdef EL_INSTANTIATE_CORE
# define EL_EXTERN
#else
# define EL_EXTERN extern
#endif

#define PROTO(T) EL_EXTERN template class DistMultiVec<T>;
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

#undef EL_EXTERN

} // namespace El
