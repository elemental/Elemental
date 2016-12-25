/*
   Copyright (c) 2009-2016, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, Stanford University, and the
   Georgia Insitute of Technology.
   All rights reserved.
 
   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

#ifdef EL_HAVE_PARMETIS
# include "parmetis.h"
#else
# include "metis.h"
#endif

namespace El {

Int Bisect
( const Graph& graph,
  Graph& leftChild,
  Graph& rightChild,
  vector<Int>& perm,
  const BisectCtrl& ctrl )
{
    EL_DEBUG_CSE
#ifdef EL_HAVE_METIS
    // METIS assumes that there are no self-connections or connections 
    // outside the sources, so we must manually remove them from our graph
    const Int numSources = graph.NumSources();
    const Int numEdges = graph.NumEdges();
    const Int* sourceBuf = graph.LockedSourceBuffer();
    const Int* targetBuf = graph.LockedTargetBuffer();
    Int numValidEdges = 0;
    for( Int i=0; i<numEdges; ++i )
        if( sourceBuf[i] != targetBuf[i] && targetBuf[i] < numSources )
            ++numValidEdges;

    // Fill our connectivity (ignoring self and too-large connections)
    vector<idx_t> xAdj( numSources+1 );
    vector<idx_t> adjacency( Max(numValidEdges,1) );
    Int validCounter=0;
    Int sourceOff=0;
    Int prevSource=-1;
    for( Int edge=0; edge<numEdges; ++edge )
    {
        const Int source = sourceBuf[edge];
        const Int target = targetBuf[edge];
        EL_DEBUG_ONLY(
          if( source < prevSource )
              RuntimeError("sources were not properly sorted");
        )
        while( source != prevSource )
        {
            xAdj[sourceOff++] = validCounter;
            ++prevSource;
        }
        if( source != target && target < numSources )
            adjacency[validCounter++] = target;
    }
    while( sourceOff <= numSources)
    { xAdj[sourceOff++] = validCounter; }

    // Call METIS_ComputeVertexSeparator, which is meant to be used by ParMETIS
    idx_t nvtxs = numSources;
    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions( options );
    options[METIS_OPTION_NSEPS] = ctrl.numSeqSeps;
    vector<idx_t> part(numSources);
    idx_t sepSize;
    METIS_ComputeVertexSeparator
    ( &nvtxs, xAdj.data(), adjacency.data(), NULL, options, 
      &sepSize, part.data() );
    
    Int sizes[3] = { 0, 0, 0 };
    for( Int s=0; s<numSources; ++s ) 
        ++sizes[part[s]];
    Int offsets[3];
    offsets[0] = 0;
    offsets[1] = sizes[0];
    offsets[2] = sizes[1] + offsets[1];
    perm.resize( numSources ); 
    for( Int s=0; s<numSources; ++s )
        perm[s] = offsets[part[s]]++;
 
    EL_DEBUG_ONLY(EnsurePermutation( perm ))
    BuildChildrenFromPerm
    ( graph, perm, sizes[0], leftChild, sizes[1], rightChild );
    return sizes[2];
#else
    LogicError("METIS interface was not available");
    return -1;
#endif
}

Int Bisect
( const DistGraph& graph, 
        unique_ptr<Grid>& childGrid,
        DistGraph& child, 
        DistMap& perm,
        bool& onLeft, 
  const BisectCtrl& ctrl )
{
    EL_DEBUG_CSE
#ifdef EL_HAVE_METIS
    const Grid& grid = graph.Grid();
    const int commSize = grid.Size();
    const int commRank = grid.Rank();
    if( commSize == 1 )
        LogicError
        ("This routine assumes at least two processes are used, "
         "otherwise one child will be lost");

    // (Par)METIS assumes that there are no self-connections or connections 
    // outside the sources, so we must manually remove them from our graph
    const Int numSources = graph.NumSources();
    const Int numLocalEdges = graph.NumLocalEdges();
    const Int* sourceBuf = graph.LockedSourceBuffer();
    const Int* targetBuf = graph.LockedTargetBuffer();
    Int numLocalValidEdges = 0;
    for( Int i=0; i<numLocalEdges; ++i )
        if( sourceBuf[i] != targetBuf[i] && targetBuf[i] < numSources )
            ++numLocalValidEdges;

    // Fill our local connectivity (ignoring self and too-large connections)
    const Int blocksize = graph.Blocksize();
    const Int numLocalSources = graph.NumLocalSources();
    const Int firstLocalSource = graph.FirstLocalSource();
    vector<idx_t> xAdj( numLocalSources+1 );
    vector<idx_t> adjacency( Max(numLocalValidEdges,1) );
    Int validCounter=0;
    Int sourceOff=0;
    Int prevSource=firstLocalSource-1;
    for( Int localEdge=0; localEdge<numLocalEdges; ++localEdge )
    {
        const Int source = sourceBuf[localEdge];
        const Int target = targetBuf[localEdge];
        EL_DEBUG_ONLY(
          if( source < prevSource )
              RuntimeError("sources were not properly sorted");
        )
        while( source != prevSource )
        {
            xAdj[sourceOff++] = validCounter;
            ++prevSource;
        }
        if( source != target && target < numSources )
            adjacency[validCounter++] = target;
    }
    while( sourceOff <= numLocalSources)
    { xAdj[sourceOff++] = validCounter; }

    vector<idx_t> sizes(3);
    if( ctrl.sequential )
    {
        // Gather the number of local valid edges on the root process
        vector<Int> edgeSizes( commSize ), edgeOffs;
        mpi::AllGather
        ( &numLocalValidEdges, 1, edgeSizes.data(), 1, grid.Comm() );
        Int numEdges=0;
        if( commRank == 0 )
        {
            edgeOffs.resize( commSize );
            numEdges = Scan( edgeSizes, edgeOffs );
        }

        // Gather the edges on the root process (with padding)
        Int maxLocalValidEdges=0;
        for( int q=0; q<commSize; ++q )
            maxLocalValidEdges = Max( maxLocalValidEdges, edgeSizes[q] );
        adjacency.resize( Max(maxLocalValidEdges,1) );
        vector<idx_t> globalAdj;
        if( commRank == 0 )
            globalAdj.resize( maxLocalValidEdges*commSize, 0 );
        mpi::Gather
        ( adjacency.data(), maxLocalValidEdges, 
          globalAdj.data(), maxLocalValidEdges, 0, grid.Comm() );

        if( commRank == 0 )
        {
            // Remove the padding
            for( int q=1; q<commSize; ++q )
            {
                const Int edgeOff = q*maxLocalValidEdges;
                for( Int j=0; j<edgeSizes[q]; ++j )
                    globalAdj[edgeOffs[q]+j] = globalAdj[edgeOff+j];
            }
            globalAdj.resize( numEdges );
        }

        // Set up the global xAdj vector
        vector<idx_t> globalXAdj;
        if( commRank == 0 )
            globalXAdj.resize( numSources+1 );
        // For now, simply loop over the processes for the receives
        if( commRank == 0 )
            for( Int j=0; j<numLocalSources; ++j )
                globalXAdj[j] = xAdj[j];    
        for( int q=1; q<commSize; ++q )
        {
            const Int thisLocalSize =
              Min(blocksize,Max(numSources-q*blocksize,0));
            if( thisLocalSize == 0 )
                break;

            if( commRank == q )
            {
                mpi::Send( xAdj.data(), numLocalSources, 0, grid.Comm() );
            }
            else if( commRank == 0 )
            {
                mpi::Recv
                ( &globalXAdj[q*blocksize], thisLocalSize, q, grid.Comm() );
                for( Int j=0; j<thisLocalSize; ++j )
                    globalXAdj[q*blocksize+j] += edgeOffs[q];
            }
        }
        if( commRank == 0 )
            globalXAdj[numSources] = numEdges;

        vector<Int> seqPerm;
        if( commRank == 0 )
        {
            // Call METIS_ComputeVertexSeparator, which is meant for ParMETIS
            idx_t nvtxs = numSources;
            idx_t options[METIS_NOPTIONS];
            METIS_SetDefaultOptions( options );
            options[METIS_OPTION_NSEPS] = ctrl.numSeqSeps;
            vector<idx_t> part(numSources);
            idx_t sepSize;

            if( globalAdj.size() > 0 )
            {
                METIS_ComputeVertexSeparator
                ( &nvtxs, globalXAdj.data(), globalAdj.data(), NULL, options,
                  &sepSize, part.data() );
            } 
            else
            {
                sepSize = 0;
                for( Int i=0; i<numSources; ++i )
                {
                    if( i <= numSources/2 )
                        part[i] = 0;
                    else
                        part[i] = 1;
                } 
            }

            for( Int j=0; j<3; ++j )
                sizes[j] = 0;
            for( Int s=0; s<numSources; ++s )
                ++sizes[part[s]];
            Int offsets[3];
            offsets[0] = 0;
            offsets[1] = sizes[0];
            offsets[2] = sizes[1] + offsets[1];
            seqPerm.resize( numSources );
            for( Int s=0; s<numSources; ++s )
                seqPerm[s] = offsets[part[s]]++;
        }

        // Set up space for the distributed permutation
        perm.SetGrid( grid );
        perm.Resize( numSources );

        // For now, loop over the processes to send the data
        if( commRank == 0 )
            for( Int j=0; j<numLocalSources; ++j )
                perm.SetLocal(j,seqPerm[j]);
        for( int q=1; q<commSize; ++q )
        {
            const Int thisLocalSize =
              Min(blocksize,Max(numSources-q*blocksize,0));
            if( thisLocalSize == 0 )
                break;

            if( commRank == 0 )
                mpi::Send
                ( &seqPerm[q*blocksize], thisLocalSize, q, grid.Comm() );
            else if( commRank == q )
                mpi::Recv( perm.Buffer(), thisLocalSize, 0, grid.Comm() );
        }

        // Broadcast the sizes information from the root
        mpi::Broadcast
        ( reinterpret_cast<byte*>(sizes.data()), 3*sizeof(idx_t), 0,
          grid.Comm() );
    }
    else
    {
#ifdef EL_HAVE_PARMETIS
        // Describe the source distribution
        vector<idx_t> vtxDist( commSize+1 );
        for( int i=0; i<commSize; ++i )
            vtxDist[i] = i*blocksize;
        vtxDist[commSize] = graph.NumSources();

        // Create space for the result
        perm.SetGrid( grid );
        perm.Resize( numSources );

        vector<idx_t> perm_idx_t( perm.NumLocalSources() );
        // Use the custom ParMETIS interface
        idx_t nseqseps = ctrl.numSeqSeps;
        idx_t nparseps = ctrl.numDistSeps;
        real_t imbalance = 1.1;
        mpi::Comm comm = grid.Comm();
        ParMETIS_ComputeVertexSeparator
        ( vtxDist.data(), xAdj.data(), adjacency.data(), &nparseps, &nseqseps, 
          &imbalance, NULL, perm_idx_t.data(), sizes.data(), &comm.comm );

        // Since idx_t might be different than Int
        std::copy( perm_idx_t.begin(), perm_idx_t.end(), perm.Buffer() );
#else
        LogicError("ParMETIS was not available");
#endif
    }
    EL_DEBUG_ONLY(EnsurePermutation( perm ))
    BuildChildFromPerm
    ( graph, perm, sizes[0], sizes[1], onLeft, childGrid, child );
    return sizes[2];
#else
    LogicError("METIS was not available");
    return -1;
#endif
}

void EnsurePermutation( const vector<Int>& map )
{
    EL_DEBUG_CSE
    const Int numSources = map.size();
    vector<Int> timesMapped( numSources, 0 );
    for( Int i=0; i<numSources; ++i )
        ++timesMapped[map[i]];
    for( Int i=0; i<numSources; ++i )
        if( timesMapped[i] != 1 )
            LogicError
            (timesMapped[i]," vertices were relabeled as ",i,
             " in sequential map");
}

void EnsurePermutation( const DistMap& map )
{
    EL_DEBUG_CSE
    const Grid& grid = map.Grid();
    const int commRank = grid.Rank();
    const Int numSources = map.NumSources();
    const Int numLocalSources = map.NumLocalSources();
    vector<Int> timesMapped( numSources, 0 );
    for( Int iLocal=0; iLocal<numLocalSources; ++iLocal )
        ++timesMapped[map.GetLocal(iLocal)];
    mpi::Reduce( timesMapped.data(), numSources, MPI_SUM, 0, grid.Comm() );
    if( commRank == 0 )
        for( Int i=0; i<numSources; ++i )
            if( timesMapped[i] != 1 )
                LogicError
                (timesMapped[i]," vertices were relabeled as ",i,
                 " in parallel map");
}

void BuildChildrenFromPerm
( const Graph& graph,
  const vector<Int>& perm, 
  Int leftChildSize, Graph& leftChild,
  Int rightChildSize, Graph& rightChild )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      const Int sepSize = graph.NumSources() - leftChildSize - rightChildSize;
    )
    const Int numSources = graph.NumSources();
    const Int numTargets = graph.NumTargets();

    // Build the inverse permutation
    vector<Int> invPerm( numSources );
    for( Int i=0; i<numSources; ++i )
        invPerm[perm[i]] = i;

    // Get an upper bound on the number of edges in the child graphs
    Int leftChildUpperBound=0, rightChildUpperBound=0;
    for( Int s=0; s<leftChildSize; ++s )
        leftChildUpperBound += graph.NumConnections( invPerm[s] );
    for( Int s=0; s<rightChildSize; ++s )
        rightChildUpperBound += 
            graph.NumConnections( invPerm[s+leftChildSize] );

    // Build the left child's graph
    leftChild.Resize( leftChildSize, numTargets );
    leftChild.Reserve( leftChildUpperBound );
    for( Int s=0; s<leftChildSize; ++s )
    {
        const Int source = s;
        const Int invSource = invPerm[s];
        const Int off = graph.SourceOffset( invSource );
        const Int numConnections = graph.NumConnections( invSource );
        for( Int t=0; t<numConnections; ++t )
        {
            const Int invTarget = graph.Target( off+t );
            const Int target = ( invTarget < numSources ? 
                                 perm[invTarget] :
                                 invTarget );
            EL_DEBUG_ONLY(
              if( target >= leftChildSize && target < (numSources-sepSize) )
                  LogicError("Invalid bisection, left set touches right set at (",source,",",target,") since leftChildSize=",leftChildSize);
            )
            leftChild.QueueConnection( source, target );
        }
    }
    leftChild.ProcessQueues();

    // Build the right child's graph
    rightChild.Resize( rightChildSize, numTargets-leftChildSize );
    rightChild.Reserve( rightChildUpperBound );
    for( Int s=0; s<rightChildSize; ++s )
    {
        const Int source = s+leftChildSize;
        const Int invSource = invPerm[source];
        const Int off = graph.SourceOffset( invSource );
        const Int numConnections = graph.NumConnections( invSource );
        for( Int t=0; t<numConnections; ++t )
        {
            const Int invTarget = graph.Target( off+t );
            const Int target = ( invTarget < numSources ?
                                 perm[invTarget] :
                                 invTarget );
            EL_DEBUG_ONLY(
              if( target < leftChildSize )
                  LogicError("Invalid bisection, right set touches left set at (",source,",",target,") since leftChildSize=",leftChildSize);
            )
            // The targets that are in parent separators do not need to be
            rightChild.QueueConnection
            ( source-leftChildSize, target-leftChildSize );
        }
    }
    rightChild.ProcessQueues();
}

void BuildChildFromPerm
( const DistGraph& graph,
  const DistMap& perm,
        Int leftChildSize,
        Int rightChildSize,
        bool& onLeft,
        unique_ptr<Grid>& childGrid,
        DistGraph& child )
{
    EL_DEBUG_CSE
    const Int numTargets = graph.NumTargets();
    const Int numLocalSources = graph.NumLocalSources();
    EL_DEBUG_ONLY(
      const Int numSources = graph.NumSources();
      const Int sepSize = numSources - leftChildSize - rightChildSize;
    )

    const Grid& grid = graph.Grid();
    const int commRank = grid.Rank();
    const int commSize = grid.Size();

    // Build the child graph from the partitioned parent
    const int smallTeamSize = commSize/2;
    const int largeTeamSize = commSize - smallTeamSize;
    const bool inSmallTeam = ( commRank < smallTeamSize );
    const bool smallOnLeft = ( leftChildSize <= rightChildSize );
    const int leftTeamSize = ( smallOnLeft ? smallTeamSize : largeTeamSize );
    const int rightTeamSize = ( smallOnLeft ? largeTeamSize : smallTeamSize );
    const int leftTeamOff = ( smallOnLeft ? 0 : smallTeamSize );
    const int rightTeamOff = ( smallOnLeft ? smallTeamSize : 0 );
    onLeft = ( inSmallTeam == smallOnLeft );

    // TODO(poulson): Generalize to 2D distributions?
    Int leftTeamBlocksize = leftChildSize / leftTeamSize;
    if( leftTeamBlocksize*leftTeamSize < leftChildSize )
        ++leftTeamBlocksize;
    // TODO(poulson): Generalize to 2D distributions?
    Int rightTeamBlocksize = rightChildSize / rightTeamSize;
    if( rightTeamBlocksize*rightTeamSize < rightChildSize )
        ++rightTeamBlocksize;

    // Count how many rows we must send to each process 
    vector<int> rowSendSizes( commSize, 0 );
    for( Int s=0; s<numLocalSources; ++s )
    {
        const Int i = perm.GetLocal(s);
        if( i < leftChildSize )
        {
            const int q = leftTeamOff + i / leftTeamBlocksize;
            ++rowSendSizes[q];
        }
        else if( i < leftChildSize+rightChildSize )
        {
            const int q = rightTeamOff + (i-leftChildSize) / rightTeamBlocksize;
            ++rowSendSizes[q];
        }
    }

    // Exchange the number of rows
    vector<int> rowRecvSizes( commSize );
    mpi::AllToAll
    ( rowSendSizes.data(), 1, rowRecvSizes.data(), 1, grid.Comm() );

    // Prepare for the AllToAll to exchange the row indices and 
    // the number of column indices per row
    vector<int> rowSendOffs, rowRecvOffs;
    const int numSendRows = Scan( rowSendSizes, rowSendOffs );
    const int numRecvRows = Scan( rowRecvSizes, rowRecvOffs );

    // Pack the row indices and how many column entries there will be per row
    vector<Int> rowSendLengths( numSendRows );
    vector<Int> rowSendInds( numSendRows );
    auto offs = rowSendOffs;
    for( Int s=0; s<numLocalSources; ++s )
    {
        const Int i = perm.GetLocal(s);
        if( i < leftChildSize )
        {
            const int q = leftTeamOff + i / leftTeamBlocksize;
            rowSendInds[offs[q]] = i;
            rowSendLengths[offs[q]] = graph.NumConnections( s );
            ++offs[q];
        }
        else if( i < leftChildSize+rightChildSize )
        {
            const int q = rightTeamOff + (i-leftChildSize) / rightTeamBlocksize;
            rowSendInds[offs[q]] = i;
            rowSendLengths[offs[q]] = graph.NumConnections( s );
            ++offs[q];
        }
    }

    // Perform the row lengths exchange
    vector<Int> rowRecvLengths( numRecvRows );
    mpi::AllToAll
    ( rowSendLengths.data(), rowSendSizes.data(), rowSendOffs.data(),
      rowRecvLengths.data(), rowRecvSizes.data(), rowRecvOffs.data(),
      grid.Comm() );

    // Perform the row indices exchange
    vector<Int> rowRecvInds( numRecvRows );
    mpi::AllToAll
    ( rowSendInds.data(), rowSendSizes.data(), rowSendOffs.data(),
      rowRecvInds.data(), rowRecvSizes.data(), rowRecvOffs.data(),
      grid.Comm() );
    SwapClear( rowSendInds );

    // Set up for sending the column indices
    int numSendInds=0;
    vector<int> indSendSizes( commSize, 0 ), indSendOffs(commSize);
    for( int q=0; q<commSize; ++q )
    {
        const int numRows = rowSendSizes[q];
        const int off = rowSendOffs[q];
        for( int s=0; s<numRows; ++s )
            indSendSizes[q] += rowSendLengths[off+s];

        indSendOffs[q] = numSendInds;
        numSendInds += indSendSizes[q];
    }
    SwapClear( rowSendLengths );
    int numRecvInds=0;
    vector<int> indRecvSizes( commSize, 0 ), indRecvOffs(commSize);
    for( int q=0; q<commSize; ++q )
    {
        const int numRows = rowRecvSizes[q];
        const int off = rowRecvOffs[q];
        for( int s=0; s<numRows; ++s )
            indRecvSizes[q] += rowRecvLengths[off+s];

        indRecvOffs[q] = numRecvInds;
        numRecvInds += indRecvSizes[q];
    }
    SwapClear( rowSendSizes );
    SwapClear( rowSendOffs );

    // Pack the indices
    vector<Int> sendInds( numSendInds );
    offs = indSendOffs;
    for( Int s=0; s<numLocalSources; ++s )
    {
        const Int i = perm.GetLocal(s);
        if( i < leftChildSize )
        {
            const int q = leftTeamOff + i / leftTeamBlocksize;

            const Int numConnections = graph.NumConnections( s );
            const Int localEdgeOff = graph.SourceOffset( s );
            for( Int j=0; j<numConnections; ++j )
                sendInds[offs[q]++] = graph.Target( localEdgeOff+j );
        }
        else if( i < leftChildSize+rightChildSize )
        {
            const int q = rightTeamOff + (i-leftChildSize) / rightTeamBlocksize;
               
            const Int numConnections = graph.NumConnections( s );
            const Int localEdgeOff = graph.SourceOffset( s );
            for( Int j=0; j<numConnections; ++j )
                sendInds[offs[q]++] = graph.Target( localEdgeOff+j );
        }
    }

    // Send/recv the column indices
    vector<Int> recvInds( numRecvInds );
    mpi::AllToAll
    ( sendInds.data(), indSendSizes.data(), indSendOffs.data(),
      recvInds.data(), indRecvSizes.data(), indRecvOffs.data(), grid.Comm() );
    SwapClear( sendInds );
    SwapClear( indSendSizes );
    SwapClear( indSendOffs );

    // Get the indices after reordering
    perm.Translate( recvInds );

    // Put the connections into our new graph
    const int childTeamRank = 
      onLeft ? commRank-leftTeamOff : commRank-rightTeamOff;
    mpi::Comm childComm;
    mpi::Split( grid.Comm(), onLeft, childTeamRank, childComm );
    // TODO(poulson): Decide on the grid dimensions.
    childGrid.reset( new Grid(childComm) );
    child.SetGrid( *childGrid );
    mpi::Free( childComm );
    if( onLeft )
        child.Resize( leftChildSize, numTargets );
    else
        child.Resize( rightChildSize, numTargets-leftChildSize );

    child.Reserve( recvInds.size() );
    Int off=0;
    for( Int s=0; s<numRecvRows; ++s )
    {
        const Int source = rowRecvInds[s];
        const Int numConnections = rowRecvLengths[s];
        const Int childFirstLocalSource = child.FirstLocalSource();
        for( Int t=0; t<numConnections; ++t )
        {
            const Int target = recvInds[off++];
            if( onLeft )
            {
                EL_DEBUG_ONLY(
                  if( target >= leftChildSize && 
                      target < (numSources-sepSize) )
                      LogicError
                      ("Invalid bisection, left set touches right:\n",
                       "  ",source," touches ",target," and leftChildSize=",
                       leftChildSize);
                )
                child.QueueLocalConnection
                ( source-childFirstLocalSource, target );
            }
            else
            {
                EL_DEBUG_ONLY(
                  if( target < leftChildSize )
                      LogicError
                      ("Invalid bisection, right set touches left set");
                )
                child.QueueLocalConnection
                ( source-leftChildSize-childFirstLocalSource, 
                  target-leftChildSize );
            }
        }
    }
    child.ProcessLocalQueues();
}

} // namespace El
