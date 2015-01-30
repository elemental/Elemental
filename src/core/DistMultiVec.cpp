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
DistMultiVec<T>::DistMultiVec()
: height_(0), width_(0), comm_(mpi::COMM_WORLD), 
  blocksize_(0), firstLocalRow_(0)
{ }

template<typename T>
DistMultiVec<T>::DistMultiVec( mpi::Comm comm )
: height_(0), width_(0), comm_(mpi::COMM_WORLD), 
  blocksize_(0), firstLocalRow_(0)
{ 
    SetComm( comm );
}

template<typename T>
DistMultiVec<T>::DistMultiVec( Int height, Int width, mpi::Comm comm )
: height_(height), width_(width), comm_(mpi::COMM_WORLD)
{ 
    SetComm( comm );
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

// Make a copy
// -----------
template<typename T>
const DistMultiVec<T>& DistMultiVec<T>::operator=( const DistMultiVec<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("DistMultiVec::operator="))
    Copy( A, *this );
    return *this;
}

template<typename T>
const DistMultiVec<T>& 
DistMultiVec<T>::operator=( const AbstractDistMatrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("DistMultiVec::operator="))
    Copy( A, *this );
    return *this;
}

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
void DistMultiVec<T>::Resize( Int height, Int width )
{
    const Int commRank = mpi::Rank( comm_ );
    const Int commSize = mpi::Size( comm_ );
    height_ = height;
    width_ = width;
    blocksize_ = height/commSize;
    firstLocalRow_ = commRank*blocksize_;
    const Int localHeight =
        ( commRank<commSize-1 ?
          blocksize_ :
          height_ - (commSize-1)*blocksize_ );
    multiVec_.Resize( localHeight, width );
}

// Change the distribution
// -----------------------
template<typename T>
void DistMultiVec<T>::SetComm( mpi::Comm comm )
{ 
    if( comm_ != mpi::COMM_WORLD )
        mpi::Free( comm_ );

    if( comm != mpi::COMM_WORLD )
        mpi::Dup( comm, comm_ );
    else
        comm_ = comm;

    const Int commRank = mpi::Rank( comm );
    const Int commSize = mpi::Size( comm );
    blocksize_ = height_/commSize;
    firstLocalRow_ = blocksize_*commRank;
    const Int localHeight =
        ( commRank<commSize-1 ?
          blocksize_ :
          height_ - (commSize-1)*blocksize_ );
    multiVec_.Resize( localHeight, width_ );
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
{ return RowToProcess( i, blocksize_, mpi::Size(comm_) ); }

template<typename T>
Int DistMultiVec<T>::GlobalRow( Int iLoc ) const
{
    DEBUG_ONLY(CallStackEntry cse("DistMultiVec::GlobalRow"))
    if( iLoc < 0 || iLoc > LocalHeight() )
        LogicError("Invalid local row index");
    return iLoc + FirstLocalRow();
}

// Detailed local information
// --------------------------
template<typename T>
T DistMultiVec<T>::Get( Int row, Int col ) const
{
    DEBUG_ONLY(CallStackEntry cse("DistMultiVec::Get"))
    int rowOwner = RowOwner(row);
    T value;
    if( rowOwner == mpi::Rank(comm_) )
        value = GetLocal( row-FirstLocalRow(), col );
    mpi::Broadcast( value, rowOwner, comm_ );
    return value;
}

template<typename T>
T DistMultiVec<T>::GetLocal( Int localRow, Int col ) const
{ 
    DEBUG_ONLY(CallStackEntry cse("DistMultiVec::GetLocal"))
    return multiVec_.Get(localRow,col);
}

template<typename T>
void DistMultiVec<T>::SetLocal( Int localRow, Int col, T value )
{
    DEBUG_ONLY(CallStackEntry cse("DistMultiVec::SetLocal"))
    multiVec_.Set(localRow,col,value);
}

template<typename T>
void DistMultiVec<T>::UpdateLocal( Int localRow, Int col, T value )
{
    DEBUG_ONLY(CallStackEntry cse("DistMultiVec::UpdateLocal"))
    multiVec_.Update(localRow,col,value);
}

#define PROTO(T) template class DistMultiVec<T>;

#include "El/macros/Instantiate.h"

} // namespace El
