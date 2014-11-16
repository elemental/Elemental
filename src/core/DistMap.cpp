/*
   Copyright (c) 2009-2014, Jack Poulson, Lexing Ying,
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

DistMap::DistMap( int numSources, mpi::Comm comm )
: numSources_(numSources), comm_(mpi::COMM_WORLD)
{ SetComm( comm ); }

DistMap::~DistMap()
{ 
    if( comm_ != mpi::COMM_WORLD )
        mpi::Free( comm_ ); 
}

void DistMap::StoreOwners
( int numSources, std::vector<int>& localInds, mpi::Comm comm )
{
    DEBUG_ONLY(CallStackEntry cse("DistMap::StoreOwners"))
    SetComm( comm );
    Resize( numSources );
    const int commSize = mpi::Size( comm );
    const int blocksize = Blocksize();
    const int firstLocalSource = FirstLocalSource();

    // Exchange via AllToAlls
    std::vector<int> sendSizes( commSize, 0 );
    const int numLocalInds = localInds.size();
    for( int s=0; s<numLocalInds; ++s )
    {
        const int i = localInds[s];
        const int q = RowToProcess( i, blocksize, commSize );
        ++sendSizes[q];
    }
    std::vector<int> recvSizes( commSize );
    mpi::AllToAll( &sendSizes[0], 1, &recvSizes[0], 1, comm );
    std::vector<int> sendOffs( commSize ), recvOffs( commSize );
    int numSends=0, numRecvs=0;
    for( int q=0; q<commSize; ++q )
    {
        sendOffs[q] = numSends;
        recvOffs[q] = numRecvs;
        numSends += sendSizes[q];
        numRecvs += recvSizes[q];
    }
    DEBUG_ONLY(
        if( numRecvs != NumLocalSources() )
            LogicError("Incorrect number of recv indices");
    )
    std::vector<int> offs = sendOffs;
    std::vector<int> sendInds( numSends );
    for( int s=0; s<numLocalInds; ++s )
    {
        const int i = localInds[s];
        const int q = RowToProcess( i, blocksize, commSize );
        sendInds[offs[q]++] = i;
    }
    std::vector<int> recvInds( numRecvs );
    mpi::AllToAll
    ( &sendInds[0], &sendSizes[0], &sendOffs[0],
      &recvInds[0], &recvSizes[0], &recvOffs[0], comm );

    // Form map
    for( int q=0; q<commSize; ++q )
    {
        const int size = recvSizes[q];
        const int off = recvOffs[q];
        for( int s=0; s<size; ++s )
        {
            const int i = recvInds[off+s];
            const int iLocal = i - firstLocalSource;
            SetLocal( iLocal, q );
        }
    }
}

void DistMap::Translate( std::vector<int>& localInds ) const
{
    DEBUG_ONLY(CallStackEntry cse("DistMap::Translate"))
    const int commSize = mpi::Size( comm_ );
    const int numLocalInds = localInds.size();

    // Count how many indices we need each process to map
    std::vector<int> requestSizes( commSize, 0 );
    for( int s=0; s<numLocalInds; ++s )
    {
        const int i = localInds[s];
        DEBUG_ONLY(
            if( i < 0 )
                LogicError("Index was negative");
        )
        if( i < numSources_ )
        {
            const int q = RowToProcess( i, blocksize_, commSize );
            ++requestSizes[q];
        }
    }

    // Send our requests and find out what we need to fulfill
    std::vector<int> fulfillSizes( commSize );
    mpi::AllToAll( &requestSizes[0], 1, &fulfillSizes[0], 1, comm_ );

    // Prepare for the AllToAll to exchange request sizes
    int numRequests=0;
    std::vector<int> requestOffs( commSize );
    for( int q=0; q<commSize; ++q )
    {
        requestOffs[q] = numRequests;
        numRequests += requestSizes[q];
    }
    int numFulfills=0;
    std::vector<int> fulfillOffs( commSize );
    for( int q=0; q<commSize; ++q )
    {
        fulfillOffs[q] = numFulfills;
        numFulfills += fulfillSizes[q];
    }

    // Pack the requested information 
    std::vector<int> requests( numRequests );
    std::vector<int> offs = requestOffs;
    for( int s=0; s<numLocalInds; ++s )
    {
        const int i = localInds[s];
        if( i < numSources_ )
        {
            const int q = RowToProcess( i, blocksize_, commSize );
            requests[offs[q]++] = i;
        }
    }

    // Perform the first index exchange
    std::vector<int> fulfills( numFulfills );
    mpi::AllToAll
    ( &requests[0], &requestSizes[0], &requestOffs[0],
      &fulfills[0], &fulfillSizes[0], &fulfillOffs[0], comm_ );

    // Map all of the indices in 'fulfills'
    for( int s=0; s<numFulfills; ++s )
    {
        const int i = fulfills[s];
        const int iLocal = i - firstLocalSource_;
        DEBUG_ONLY(
            if( iLocal < 0 || iLocal >= (int)map_.size() )
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
    ( &fulfills[0], &fulfillSizes[0], &fulfillOffs[0],
      &requests[0], &requestSizes[0], &requestOffs[0], comm_ );

    // Unpack in the same way we originally packed
    offs = requestOffs;
    for( int s=0; s<numLocalInds; ++s )
    {
        const int i = localInds[s];
        if( i < numSources_ )
        {
            const int q = RowToProcess( i, blocksize_, commSize );
            localInds[s] = requests[offs[q]++];
        }
    }
}

void DistMap::FormInverse( DistMap& inverseMap ) const
{
    DEBUG_ONLY(CallStackEntry cse("DistMap::FormInverse"))
    const int commSize = mpi::Size( comm_ );
    const int numLocalSources = map_.size();

    // How many pairs of original and mapped indices to send to each process
    std::vector<int> sendSizes( commSize, 0 );
    for( int s=0; s<numLocalSources; ++s )
    {
        const int i = map_[s];
        const int q = RowToProcess( i, blocksize_, commSize );
        sendSizes[q] += 2;
    }

    // Coordinate all of the processes on their send sizes
    std::vector<int> recvSizes( commSize );
    mpi::AllToAll( &sendSizes[0], 1, &recvSizes[0], 1, comm_ );

    // Prepare for the AllToAll to exchange send sizes
    int numSends=0;
    std::vector<int> sendOffs( commSize );
    for( int q=0; q<commSize; ++q )
    {
        sendOffs[q] = numSends;
        numSends += sendSizes[q];
    }
    DEBUG_ONLY(
        if( numSends != 2*numLocalSources )
            LogicError("Miscalculated numSends");
    )
    int numReceives=0;
    std::vector<int> recvOffs( commSize );
    for( int q=0; q<commSize; ++q )
    {
        recvOffs[q] = numReceives;
        numReceives += recvSizes[q];
    }
    DEBUG_ONLY(
        if( numReceives != 2*numLocalSources )
            LogicError("Mistake in number of receives");
    )

    // Pack our map information
    std::vector<int> sends( numSends );
    std::vector<int> offs = sendOffs;
    for( int s=0; s<numLocalSources; ++s )
    {
        const int i = map_[s];
        const int q = RowToProcess( i, blocksize_, commSize );
        sends[offs[q]++] = s+firstLocalSource_;
        sends[offs[q]++] = i;
    }

    // Send out the map information
    std::vector<int> recvs( numReceives );
    mpi::AllToAll
    ( &sends[0], &sendSizes[0], &sendOffs[0],
      &recvs[0], &recvSizes[0], &recvOffs[0], comm_ );

    // Form our part of the inverse map
    inverseMap.numSources_ = numSources_;
    inverseMap.SetComm( comm_ );
    for( int s=0; s<numReceives; s+=2 )
    {
        const int origInd = recvs[s];
        const int mappedInd = recvs[s+1];
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

int DistMap::NumSources() const { return numSources_; }

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
    const int numLocalSources =
        ( commRank<commSize-1 ?
          blocksize_ :
          numSources_ - (commSize-1)*blocksize_ );
    map_.resize( numLocalSources );
}

mpi::Comm DistMap::Comm() const { return comm_; }

int DistMap::Blocksize() const { return blocksize_; }

int DistMap::FirstLocalSource() const { return firstLocalSource_; }

int DistMap::NumLocalSources() const { return map_.size(); }

int DistMap::RowOwner( int i ) const 
{ return RowToProcess( i, blocksize_, mpi::Size(comm_) ); }

int DistMap::GetLocal( int localSource ) const
{ 
    DEBUG_ONLY(
        CallStackEntry cse("DistMap::GetLocal");
        if( localSource < 0 || localSource >= (int)map_.size() )
            LogicError("local source is out of bounds");
    )
    return map_[localSource];
}

void DistMap::SetLocal( int localSource, int target )
{
    DEBUG_ONLY(
        CallStackEntry cse("DistMap::SetLocal");
        if( localSource < 0 || localSource >= (int)map_.size() )
            LogicError("local source is out of bounds");
    )
    map_[localSource] = target; 
}

      std::vector<int>& DistMap::Map()       { return map_; }
const std::vector<int>& DistMap::Map() const { return map_; }

      int* DistMap::Buffer()       { return &map_[0]; }
const int* DistMap::Buffer() const { return &map_[0]; }

void DistMap::Empty()
{
    numSources_ = 0;
    blocksize_ = 0;
    firstLocalSource_ = 0;
    SwapClear( map_ );
}

void DistMap::Resize( int numSources )
{
    const int commRank = mpi::Rank( comm_ );
    const int commSize = mpi::Size( comm_ );
    numSources_ = numSources;
    blocksize_ = numSources/commSize;
    firstLocalSource_ = commRank*blocksize_;
    const int numLocalSources = 
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
