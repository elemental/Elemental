/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// TODO: Avoid sorting since the ordering can easily be preserved

void GetSubgraph
( const Graph& graph,
        Range<Int> I,
        Range<Int> J,
        Graph& subgraph )
{
    DEBUG_ONLY(CSE cse("GetSubgraph"))
    const Int* offsetBuf = graph.LockedOffsetBuffer();
    const Int* targetBuf = graph.LockedTargetBuffer();

    if( I.end == END )
        I.end = graph.NumSources();
    if( J.end == END )
        J.end = graph.NumTargets();
    const Int mSub = I.end-I.beg;
    const Int nSub = J.end-J.beg;
    subgraph.Empty();
    subgraph.Resize( mSub, nSub );

    // Reserve the number of edges that live within the subgraph
    Int numEdgesSub = 0;
    for( Int i=I.beg; i<I.end; ++i )
    {
        const Int offset = offsetBuf[i];
        const Int numConn = offsetBuf[i+1] - offset;
        for( Int e=offset; e<offset+numConn; ++e )
        {
            const Int j = targetBuf[e];
            if( j >= J.beg && j < J.end )
                ++numEdgesSub;
        }
    }
    subgraph.Reserve( numEdgesSub );

    // Insert the edges
    for( Int i=I.beg; i<I.end; ++i ) 
    {
        const Int offset = offsetBuf[i];
        const Int numConn = offsetBuf[i+1] - offset;
        for( Int e=offset; e<offset+numConn; ++e )
        {
            const Int j = targetBuf[e];
            if( j >= J.beg && j < J.end )
                subgraph.QueueConnection( i-I.beg, j-J.beg );
        }
    }
    subgraph.ProcessQueues();
}

void GetSubgraph
( const Graph& graph,
        Range<Int> I,
  const vector<Int>& J,
        Graph& subgraph )
{
    DEBUG_ONLY(CSE cse("GetSubgraph"))
    // TODO: Decide how to handle unsorted I and J with duplicates
    LogicError("This routine is not yet written");
}

void GetSubgraph
( const Graph& graph,
  const vector<Int>& I,
        Range<Int> J,
        Graph& subgraph )
{
    DEBUG_ONLY(CSE cse("GetSubgraph"))
    // TODO: Decide how to handle unsorted I and J with duplicates
    LogicError("This routine is not yet written");
}

void GetSubgraph
( const Graph& graph,
  const vector<Int>& I,
  const vector<Int>& J,
        Graph& subgraph )
{
    DEBUG_ONLY(CSE cse("GetSubgraph"))
    // TODO: Decide how to handle unsorted I and J with duplicates
    LogicError("This routine is not yet written");
}

void GetSubgraph
( const DistGraph& graph,
        Range<Int> I,
        Range<Int> J,
        DistGraph& subgraph )
{
    DEBUG_ONLY(CSE cse("GetSubgraph"))
    const Int* targetBuf = graph.LockedTargetBuffer();
    const Int* sourceBuf = graph.LockedSourceBuffer();
    if( I.end == END )
        I.end = graph.NumSources();
    if( J.end == END )
        J.end = graph.NumTargets();
    const Int mSub = I.end-I.beg;
    const Int nSub = J.end-J.beg;
    const Int numEdges = graph.NumLocalEdges();

    mpi::Comm comm = graph.Comm();
    const int commSize = mpi::Size( comm );
    subgraph.Empty();
    subgraph.SetComm( comm );
    subgraph.Resize( mSub, nSub );

    // Compute the metadata
    // ====================
    vector<int> sendCounts(commSize,0);
    for( Int e=0; e<numEdges; ++e )
    {
        const Int i = sourceBuf[e];
        const Int j = targetBuf[e];
        if( i >= I.end )
            break;
        else if( i >= I.beg && j >= J.beg && j < J.end )
            ++sendCounts[ subgraph.SourceOwner(i-I.beg) ];
    }

    // Pack the data
    // =============
    vector<int> sendOffs;
    const int totalSend = Scan( sendCounts, sendOffs );
    auto offs = sendOffs;
    vector<Int> sendSources(totalSend), sendTargets(totalSend);
    for( Int e=0; e<numEdges; ++e )
    {
        const Int i = sourceBuf[e];
        const Int j = targetBuf[e];
        if( i >= I.end )
            break;
        else if( i >= I.beg && j >= J.beg && j < J.end )
        {
            const Int iSub = i - I.beg;
            const Int jSub = j - J.beg;
            const int owner = subgraph.SourceOwner( iSub );
            sendSources[offs[owner]] = iSub;
            sendTargets[offs[owner]] = jSub;
            ++offs[owner];
        }
    }
    
    // Exchange and unpack the data
    // ============================
    auto recvSources = mpi::AllToAll( sendSources, sendCounts, sendOffs, comm );
    auto recvTargets = mpi::AllToAll( sendTargets, sendCounts, sendOffs, comm );
    const Int totalRecv = recvSources.size();
    subgraph.Reserve( totalRecv );
    for( Int i=0; i<totalRecv; ++i )
        subgraph.QueueConnection( recvSources[i], recvTargets[i] );
    subgraph.ProcessQueues();
}

void GetSubgraph
( const DistGraph& graph,
        Range<Int> I,
  const vector<Int>& J,
        DistGraph& subgraph )
{
    DEBUG_ONLY(CSE cse("GetSubgraph"))
    // TODO: Decide how to handle unsorted I and J with duplicates
    LogicError("This routine is not yet written");
}

void GetSubgraph
( const DistGraph& graph,
  const vector<Int>& I,
        Range<Int> J,
        DistGraph& subgraph )
{
    DEBUG_ONLY(CSE cse("GetSubgraph"))
    // TODO: Decide how to handle unsorted I and J with duplicates
    LogicError("This routine is not yet written");
}

void GetSubgraph
( const DistGraph& graph,
  const vector<Int>& I,
  const vector<Int>& J,
        DistGraph& subgraph )
{
    DEBUG_ONLY(CSE cse("GetSubgraph"))
    // TODO: Decide how to handle unsorted I and J with duplicates
    LogicError("This routine is not yet written");
}

} // namespace El
