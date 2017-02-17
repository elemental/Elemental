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
#include <El/blas_like/level1/Axpy.hpp>
#include <El/blas_like/level1/Scale.hpp>

namespace El {

// Constructors and destructors
// ============================

template<typename Ring>
DistMultiVec<Ring>::DistMultiVec( const El::Grid& grid )
: grid_(&grid)
{
    EL_DEBUG_CSE
    InitializeLocalData();
}

template<typename Ring>
DistMultiVec<Ring>::DistMultiVec( Int height, Int width, const El::Grid& grid )
: height_(height), width_(width), grid_(&grid)
{
    InitializeLocalData();
}

template<typename Ring>
DistMultiVec<Ring>::DistMultiVec( const DistMultiVec<Ring>& A )
{
    EL_DEBUG_CSE
    grid_ = &A.Grid();
    if( &A != this )
        *this = A;
    EL_DEBUG_ONLY(
      else
          LogicError("Tried to construct DistMultiVec via itself");
    )
}

template<typename Ring>
DistMultiVec<Ring>::~DistMultiVec() { }

// Assignment and reconfiguration
// ==============================

// Change the size of the matrix
// -----------------------------
template<typename Ring>
void DistMultiVec<Ring>::Empty( bool freeMemory )
{
    EL_DEBUG_CSE
    height_ = 0;
    width_ = 0;
    blocksize_ = 1;

    multiVec_.Empty( freeMemory );

    SwapClear( remoteUpdates_ );
}

template<typename Ring>
void DistMultiVec<Ring>::InitializeLocalData()
{
    EL_DEBUG_CSE
    const int gridRank = grid_->Rank();
    const int gridSize = grid_->Size();
    blocksize_ = height_ / gridSize;
    if( blocksize_*gridSize < height_ || height_ == 0 )
        ++blocksize_;
    const Int localHeight = Min(blocksize_,Max(0,height_-blocksize_*gridRank));
    multiVec_.Resize( localHeight, width_ );
}

template<typename Ring>
void DistMultiVec<Ring>::Resize( Int height, Int width )
{
    EL_DEBUG_CSE
    if( height_ == height && width == width_ )
        return;

    height_ = height;
    width_ = width;
    InitializeLocalData();

    SwapClear( remoteUpdates_ );
}

// Change the distribution
// -----------------------
template<typename Ring>
void DistMultiVec<Ring>::SetGrid( const El::Grid& grid )
{
    EL_DEBUG_CSE
    if( grid_ == &grid )
        return;
    grid_ = &grid;
    Resize( 0, 0 );
}

// Operator overloading
// ====================

// Make a copy
// -----------
template<typename Ring>
const DistMultiVec<Ring>& DistMultiVec<Ring>::operator=
( const DistMultiVec<Ring>& A )
{
    EL_DEBUG_CSE
    Copy( A, *this );
    return *this;
}

template<typename Ring>
const DistMultiVec<Ring>&
DistMultiVec<Ring>::operator=( const AbstractDistMatrix<Ring>& A )
{
    EL_DEBUG_CSE
    Copy( A, *this );
    return *this;
}

// Make a copy of a submatrix
// --------------------------
template<typename Ring>
DistMultiVec<Ring>
DistMultiVec<Ring>::operator()( Range<Int> I, Range<Int> J ) const
{
    EL_DEBUG_CSE
    DistMultiVec<Ring> ASub(this->Grid());
    GetSubmatrix( *this, I, J, ASub );
    return ASub;
}

template<typename Ring>
DistMultiVec<Ring>
DistMultiVec<Ring>::operator()( Range<Int> I, const vector<Int>& J ) const
{
    EL_DEBUG_CSE
    DistMultiVec<Ring> ASub(this->Grid());
    GetSubmatrix( *this, I, J, ASub );
    return ASub;
}

template<typename Ring>
DistMultiVec<Ring>
DistMultiVec<Ring>::operator()( const vector<Int>& I, Range<Int> J ) const
{
    EL_DEBUG_CSE
    DistMultiVec<Ring> ASub(this->Grid());
    GetSubmatrix( *this, I, J, ASub );
    return ASub;
}

template<typename Ring>
DistMultiVec<Ring>
DistMultiVec<Ring>::operator()
( const vector<Int>& I, const vector<Int>& J ) const
{
    EL_DEBUG_CSE
    DistMultiVec<Ring> ASub(this->Grid());
    GetSubmatrix( *this, I, J, ASub );
    return ASub;
}

// Rescaling
// ---------
template<typename Ring>
const DistMultiVec<Ring>& DistMultiVec<Ring>::operator*=( const Ring& alpha )
{
    EL_DEBUG_CSE
    Scale( alpha, *this );
    return *this;
}

// Addition/subtraction
// --------------------
template<typename Ring>
const DistMultiVec<Ring>& DistMultiVec<Ring>::operator+=
( const DistMultiVec<Ring>& A )
{
    EL_DEBUG_CSE
    Axpy( Ring(1), A, *this );
    return *this;
}

template<typename Ring>
const DistMultiVec<Ring>& DistMultiVec<Ring>::operator-=
( const DistMultiVec<Ring>& A )
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
Int DistMultiVec<Ring>::Height() const EL_NO_EXCEPT
{ return height_; }
template<typename Ring>
Int DistMultiVec<Ring>::Width() const EL_NO_EXCEPT
{ return multiVec_.Width(); }
template<typename Ring>
Int DistMultiVec<Ring>::FirstLocalRow() const EL_NO_EXCEPT
{ return blocksize_*grid_->Rank(); }
template<typename Ring>
Int DistMultiVec<Ring>::LocalHeight() const EL_NO_EXCEPT
{ return multiVec_.Height(); }
template<typename Ring>
El::Matrix<Ring>& DistMultiVec<Ring>::Matrix() EL_NO_EXCEPT
{ return multiVec_; }
template<typename Ring>
const El::Matrix<Ring>& DistMultiVec<Ring>::LockedMatrix() const EL_NO_EXCEPT
{ return multiVec_; }

// Distribution information
// ------------------------
template<typename Ring>
const El::Grid& DistMultiVec<Ring>::Grid() const EL_NO_EXCEPT { return *grid_; }

template<typename Ring>
Int DistMultiVec<Ring>::Blocksize() const EL_NO_EXCEPT { return blocksize_; }

template<typename Ring>
int DistMultiVec<Ring>::RowOwner( Int i ) const EL_NO_EXCEPT
{
    if( i == END ) i = height_ - 1;
    return i / blocksize_;
}

template<typename Ring>
int DistMultiVec<Ring>::Owner( Int i, Int j ) const EL_NO_EXCEPT
{ return RowOwner(i); }

template<typename Ring>
bool DistMultiVec<Ring>::IsLocalRow( Int i ) const EL_NO_EXCEPT
{ return RowOwner(i) == grid_->Rank(); }

template<typename Ring>
bool DistMultiVec<Ring>::IsLocal( Int i, Int j ) const EL_NO_EXCEPT
{ return IsLocalRow(i); }

template<typename Ring>
Int DistMultiVec<Ring>::GlobalRow( Int iLoc ) const
{
    EL_DEBUG_CSE
    if( iLoc == END ) iLoc = LocalHeight() - 1;
    EL_DEBUG_ONLY(
      if( iLoc < 0 || iLoc >= LocalHeight() )
          LogicError("Invalid local row index");
    )
    return iLoc + FirstLocalRow();
}

template<typename Ring>
Int DistMultiVec<Ring>::LocalRow( Int i ) const
{
    EL_DEBUG_CSE
    if( i == END ) i = Height() - 1;
    EL_DEBUG_ONLY(
      if( i < 0 || i >= Height() )
          LogicError("Invalid global row index");
    )
    return i - FirstLocalRow();
}

// Detailed local information
// --------------------------
template<typename Ring>
Ring DistMultiVec<Ring>::Get( Int i, Int j ) const
{
    EL_DEBUG_CSE
    if( i == END ) i = height_ - 1;
    const int rowOwner = RowOwner(i);
    Ring value;
    if( rowOwner == grid_->Rank() )
        value = GetLocal( i-FirstLocalRow(), j );
    mpi::Broadcast( value, rowOwner, grid_->Comm() );
    return value;
}

template<typename Ring>
void DistMultiVec<Ring>::Set( Int i, Int j, const Ring& value )
{
    EL_DEBUG_CSE
    if( i == END ) i = height_ - 1;
    const Int firstLocalRow = FirstLocalRow();
    if( i >= firstLocalRow && i < firstLocalRow+LocalHeight() )
        SetLocal( i-firstLocalRow, j, value );
}

template<typename Ring>
void DistMultiVec<Ring>::Set( const Entry<Ring>& entry )
{ Set( entry.i, entry.j, entry.value ); }

template<typename Ring>
void DistMultiVec<Ring>::Update( Int i, Int j, const Ring& value )
{
    EL_DEBUG_CSE
    if( i == END ) i = height_ - 1;
    const Int firstLocalRow = FirstLocalRow();
    if( i >= firstLocalRow && i < firstLocalRow+LocalHeight() )
        UpdateLocal( i-firstLocalRow, j, value );
}

template<typename Ring>
void DistMultiVec<Ring>::Update( const Entry<Ring>& entry )
{ Update( entry.i, entry.j, entry.value ); }

template<typename Ring>
Ring DistMultiVec<Ring>::GetLocal( Int iLoc, Int j ) const
{
    EL_DEBUG_CSE
    return multiVec_.Get(iLoc,j);
}

template<typename Ring>
void DistMultiVec<Ring>::SetLocal( Int iLoc, Int j, const Ring& value )
{
    EL_DEBUG_CSE
    multiVec_.Set(iLoc,j,value);
}

template<typename Ring>
void DistMultiVec<Ring>::SetLocal( const Entry<Ring>& localEntry )
{ SetLocal( localEntry.i, localEntry.j, localEntry.value ); }

template<typename Ring>
void DistMultiVec<Ring>::UpdateLocal( Int iLoc, Int j, const Ring& value )
{
    EL_DEBUG_CSE
    multiVec_.Update(iLoc,j,value);
}

template<typename Ring>
void DistMultiVec<Ring>::UpdateLocal( const Entry<Ring>& localEntry )
{ UpdateLocal( localEntry.i, localEntry.j, localEntry.value ); }

// Batch remote updates
// --------------------
template<typename Ring>
void DistMultiVec<Ring>::Reserve( Int numRemoteEntries )
{
    EL_DEBUG_CSE
    const Int currSize = remoteUpdates_.size();
    remoteUpdates_.reserve( currSize+numRemoteEntries );
}

template<typename Ring>
void DistMultiVec<Ring>::QueueUpdate( const Entry<Ring>& entry )
{
    EL_DEBUG_CSE
    remoteUpdates_.push_back( entry );
}

template<typename Ring>
void DistMultiVec<Ring>::QueueUpdate( Int i, Int j, const Ring& value )
{ QueueUpdate( Entry<Ring>{i,j,value} ); }

template<typename Ring>
void DistMultiVec<Ring>::ProcessQueues()
{
    EL_DEBUG_CSE
    // Compute the send counts
    // -----------------------
    vector<int> sendCounts(grid_->Size());
    for( const auto& entry : remoteUpdates_ )
        ++sendCounts[Owner(entry.i,entry.j)];
    // Pack the send data
    // ------------------
    vector<int> sendOffs;
    const int totalSend = Scan( sendCounts, sendOffs );
    auto offs = sendOffs;
    vector<Entry<Ring>> sendEntries(totalSend);
    for( const auto& entry : remoteUpdates_ )
        sendEntries[offs[Owner(entry.i,entry.j)]++] = entry;
    SwapClear( remoteUpdates_ );
    // Exchange and unpack
    // -------------------
    auto recvEntries =
      mpi::AllToAll( sendEntries, sendCounts, sendOffs, grid_->Comm() );

    Ring* matBuf = multiVec_.Buffer();
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

#define PROTO(Ring) EL_EXTERN template class DistMultiVec<Ring>;
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

#undef EL_EXTERN

} // namespace El
