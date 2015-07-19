/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

void GetSubgraph
( const Graph& graph, Range<Int> I, Range<Int> J,
        Graph& subgraph )
{
    DEBUG_ONLY(CSE cse("GetSubgraph"))
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
        const Int sourceOff = graph.SourceOffset(i);
        const Int numConn = graph.NumConnections(i);
        for( Int e=sourceOff; e<sourceOff+numConn; ++e )
        {
            const Int j = graph.Target(e);
            if( j >= J.beg && j < J.end )
                ++numEdgesSub;
        }
    }
    subgraph.Reserve( numEdgesSub );

    // Insert the edges
    for( Int i=I.beg; i<I.end; ++i ) 
    {
        const Int sourceOff = graph.SourceOffset(i);
        const Int numConn = graph.NumConnections(i);
        for( Int e=sourceOff; e<sourceOff+numConn; ++e )
        {
            const Int j = graph.Target(e);
            if( j >= J.beg && j < J.end )
                subgraph.QueueConnection( i-I.beg, j-J.beg );
        }
    }
    subgraph.ProcessQueues();
}

Graph GetSubgraph( const Graph& graph, Range<Int> I, Range<Int> J )
{
    Graph subgraph;
    GetSubgraph( graph, I, J, subgraph );
    return subgraph;
}

void GetSubgraph
( const DistGraph& graph, Range<Int> I, Range<Int> J,
        DistGraph& subgraph )
{
    DEBUG_ONLY(CSE cse("GetSubgraph"))
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
        const Int i = graph.Source(e);
        const Int j = graph.Target(e);
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
        const Int i = graph.Source(e);
        const Int j = graph.Target(e);
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
    subgraph.Reserve( recvSources.size() );
    for( Int i=0; i<recvSources.size(); ++i )
        subgraph.QueueConnection( recvSources[i], recvTargets[i] );
    subgraph.ProcessQueues();
}

DistGraph GetSubgraph( const DistGraph& graph, Range<Int> I, Range<Int> J )
{
    DistGraph subgraph(graph.Comm());
    GetSubgraph( graph, I, J, subgraph );
    return subgraph;
}

} // namespace El
