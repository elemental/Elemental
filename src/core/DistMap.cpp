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

DistMap::DistMap()
: numSources_(0), comm_(mpi::COMM_WORLD)
{ SetComm( mpi::COMM_WORLD ); } 

DistMap::DistMap( mpi::Comm comm )
: numSources_(0), comm_(mpi::COMM_WORLD)
{ SetComm( comm ); }

DistMap::DistMap( Int numSources, mpi::Comm comm )
: numSources_(numSources), comm_(mpi::COMM_WORLD)
{ SetComm( comm ); }

DistMap::~DistMap()
{ 
    if( !mpi::Finalized() )
        if( comm_ != mpi::COMM_WORLD )
            mpi::Free( comm_ ); 
}

void DistMap::StoreOwners
( Int numSources, std::vector<Int>& localInds, mpi::Comm comm )
{
    DEBUG_ONLY(CallStackEntry cse("DistMap::StoreOwners"))
    SetComm( comm );
    Resize( numSources );
    const int commSize = mpi::Size( comm );

    // Exchange via AllToAlls
    std::vector<int> sendSizes( commSize, 0 );
    const Int numLocalInds = localInds.size();
    for( Int s=0; s<numLocalInds; ++s )
        ++sendSizes[ RowOwner(localInds[s]) ]; 
    std::vector<int> recvSizes( commSize );
    mpi::AllToAll( sendSizes.data(), 1, recvSizes.data(), 1, comm );
    std::vector<int> sendOffs, recvOffs;
    const int numSends = Scan( sendSizes, sendOffs );
    const int numRecvs = Scan( recvSizes, recvOffs );
    DEBUG_ONLY(
        if( numRecvs != NumLocalSources() )
            LogicError("Incorrect number of recv indices");
    )
    auto offs = sendOffs;
    std::vector<Int> sendInds( numSends );
    for( Int s=0; s<numLocalInds; ++s )
    {
        const Int i = localInds[s];
        sendInds[offs[RowOwner(i)]++] = i;
    }
    std::vector<Int> recvInds( numRecvs );
    mpi::AllToAll
    ( sendInds.data(), sendSizes.data(), sendOffs.data(),
      recvInds.data(), recvSizes.data(), recvOffs.data(), comm );

    // Form map
    const Int firstLocalSource = FirstLocalSource();
    for( int q=0; q<commSize; ++q )
    {
        const int size = recvSizes[q];
        const int off = recvOffs[q];
        for( int s=0; s<size; ++s )
            SetLocal( recvInds[off+s]-firstLocalSource, q );
    }
}

void DistMap::Translate( std::vector<Int>& localInds ) const
{
    DEBUG_ONLY(CallStackEntry cse("DistMap::Translate"))
    const int commSize = mpi::Size( comm_ );
    const Int numLocalInds = localInds.size();

    // Count how many indices we need each process to map
    std::vector<int> requestSizes( commSize, 0 );
    for( Int s=0; s<numLocalInds; ++s )
    {
        const Int i = localInds[s];
        if( i < numSources_ )
            ++requestSizes[ RowOwner(i) ];
    }

    // Send our requests and find out what we need to fulfill
    std::vector<int> fulfillSizes( commSize );
    mpi::AllToAll( requestSizes.data(), 1, fulfillSizes.data(), 1, comm_ );

    // Prepare for the AllToAll to exchange request sizes
    std::vector<int> requestOffs, fulfillOffs;
    const int numRequests = Scan( requestSizes, requestOffs );
    const int numFulfills = Scan( fulfillSizes, fulfillOffs );

    // Pack the requested information 
    std::vector<int> requests( numRequests );
    auto offs = requestOffs;
    for( Int s=0; s<numLocalInds; ++s )
    {
        const Int i = localInds[s];
        if( i < numSources_ )
            requests[offs[RowOwner(i)]++] = i;
    }

    // Perform the first index exchange
    std::vector<int> fulfills( numFulfills );
    mpi::AllToAll
    ( requests.data(), requestSizes.data(), requestOffs.data(),
      fulfills.data(), fulfillSizes.data(), fulfillOffs.data(), comm_ );

    // Map all of the indices in 'fulfills'
    for( int s=0; s<numFulfills; ++s )
    {
        const Int i = fulfills[s];
        const Int iLocal = i - firstLocalSource_;
        DEBUG_ONLY(
            if( iLocal < 0 || iLocal >= (Int)map_.size() )
            {
                const int commRank = mpi::Rank( comm_ );
                LogicError
                ("invalid request: i=",i,", iLocal=",iLocal,
                 ", commRank=",commRank,", blocksize=",blocksize_);
            }
        )
        fulfills[s] = map_[iLocal];
    }

    // Send everything back
    mpi::AllToAll
    ( fulfills.data(), fulfillSizes.data(), fulfillOffs.data(),
      requests.data(), requestSizes.data(), requestOffs.data(), comm_ );

    // Unpack in the same way we originally packed
    offs = requestOffs;
    for( Int s=0; s<numLocalInds; ++s )
    {
        const Int i = localInds[s];
        if( i < numSources_ )
            localInds[s] = requests[offs[RowOwner(i)]++];
    }
}

void DistMap::FormInverse( DistMap& inverseMap ) const
{
    DEBUG_ONLY(CallStackEntry cse("DistMap::FormInverse"))
    const int commSize = mpi::Size( comm_ );
    const Int numLocalSources = map_.size();

    // How many pairs of original and mapped indices to send to each process
    std::vector<int> sendSizes( commSize, 0 );
    for( Int s=0; s<numLocalSources; ++s )
        sendSizes[RowOwner(map_[s])] += 2;

    // Coordinate all of the processes on their send sizes
    std::vector<int> recvSizes( commSize );
    mpi::AllToAll( sendSizes.data(), 1, recvSizes.data(), 1, comm_ );

    // Prepare for the AllToAll to exchange send sizes
    std::vector<int> sendOffs, recvOffs;
    const int numSends = Scan( sendSizes, sendOffs );
    const int numRecvs = Scan( recvSizes, recvOffs );
    DEBUG_ONLY(
        if( numSends != 2*numLocalSources )
            LogicError("Miscalculated numSends");
    )
    DEBUG_ONLY(
        if( numRecvs != 2*numLocalSources )
            LogicError("Mistake in number of receives");
    )

    // Pack our map information
    std::vector<int> sends( numSends );
    auto offs = sendOffs;
    for( Int s=0; s<numLocalSources; ++s )
    {
        const Int i = map_[s];
        const int q = RowOwner(i);
        sends[offs[q]++] = s+firstLocalSource_;
        sends[offs[q]++] = i;
    }

    // Send out the map information
    std::vector<int> recvs( numRecvs );
    mpi::AllToAll
    ( sends.data(), sendSizes.data(), sendOffs.data(),
      recvs.data(), recvSizes.data(), recvOffs.data(), comm_ );

    // Form our part of the inverse map
    inverseMap.numSources_ = numSources_;
    inverseMap.SetComm( comm_ );
    for( Int s=0; s<numRecvs; s+=2 )
    {
        const Int origInd = recvs[s];
        const Int mappedInd = recvs[s+1];
        inverseMap.SetLocal( mappedInd-firstLocalSource_, origInd );
    }
}

void DistMap::Extend( DistMap& firstMap ) const
{
    DEBUG_ONLY(
        CallStackEntry cse("DistMap::Extend");
        // TODO: Ensure that the communicators are congruent and that the maps
        //       are compatible sizes.
    )
    Translate( firstMap.map_ ); 
}

void DistMap::Extend( const DistMap& firstMap, DistMap& compositeMap ) const
{
    DEBUG_ONLY(CallStackEntry cse("DistMap::Extend"))
    compositeMap = firstMap;
    Extend( compositeMap );
}

Int DistMap::NumSources() const { return numSources_; }

void DistMap::SetComm( mpi::Comm comm )
{
    if( comm_ != mpi::COMM_WORLD )
        mpi::Free( comm_ );

    if( comm != mpi::COMM_WORLD )
        mpi::Dup( comm, comm_ );
    else
        comm_ = comm;

    const int commRank = mpi::Rank( comm );
    const int commSize = mpi::Size( comm );
    blocksize_ = numSources_/commSize;
    firstLocalSource_ = blocksize_*commRank;
    const Int numLocalSources =
        ( commRank<commSize-1 ?
          blocksize_ :
          numSources_ - (commSize-1)*blocksize_ );
    map_.resize( numLocalSources );
}

mpi::Comm DistMap::Comm() const { return comm_; }

Int DistMap::Blocksize() const { return blocksize_; }

Int DistMap::FirstLocalSource() const { return firstLocalSource_; }

Int DistMap::NumLocalSources() const { return map_.size(); }

int DistMap::RowOwner( Int i ) const 
{ return RowToProcess( i, blocksize_, mpi::Size(comm_) ); }

Int DistMap::GetLocal( Int localSource ) const
{ 
    DEBUG_ONLY(
        CallStackEntry cse("DistMap::GetLocal");
        if( localSource < 0 || localSource >= (Int)map_.size() )
            LogicError("local source is out of bounds");
    )
    return map_[localSource];
}

void DistMap::SetLocal( Int localSource, Int target )
{
    DEBUG_ONLY(
        CallStackEntry cse("DistMap::SetLocal");
        if( localSource < 0 || localSource >= (Int)map_.size() )
            LogicError("local source is out of bounds");
    )
    map_[localSource] = target; 
}

      std::vector<Int>& DistMap::Map()       { return map_; }
const std::vector<Int>& DistMap::Map() const { return map_; }

      Int* DistMap::Buffer()       { return map_.data(); }
const Int* DistMap::Buffer() const { return map_.data(); }

void DistMap::Empty()
{
    numSources_ = 0;
    blocksize_ = 0;
    firstLocalSource_ = 0;
    SwapClear( map_ );
}

void DistMap::Resize( Int numSources )
{
    const int commRank = mpi::Rank( comm_ );
    const int commSize = mpi::Size( comm_ );
    numSources_ = numSources;
    blocksize_ = numSources/commSize;
    firstLocalSource_ = commRank*blocksize_;
    const Int numLocalSources = 
        ( commRank<commSize-1 ? blocksize_
                              : numSources-blocksize_*(commSize-1) );
    map_.resize( numLocalSources );
}

const DistMap& DistMap::operator=( const DistMap& map )
{
    numSources_ = map.numSources_;
    SetComm( map.comm_ );
    map_ = map.map_;
    return *this;
}

} // namespace El
