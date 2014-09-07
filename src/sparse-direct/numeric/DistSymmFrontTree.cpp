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

template<typename T>
DistSymmFrontTree<T>::DistSymmFrontTree()
{ }

template<typename T>
void DistSymmFrontTree<T>::Initialize
( const DistSparseMatrix<T>& A, 
  const DistMap& reordering,
  const DistSeparatorTree& sepTree, 
  const DistSymmInfo& info,
  bool conjugate )
{
    DEBUG_ONLY(
        CallStackEntry cse("DistSymmFrontTree::Initialize");
        if( A.LocalHeight() != reordering.NumLocalSources() )
            LogicError("Local mapping was not the right size");
    )
    frontType = SYMM_2D;
    isHermitian = conjugate;
    
    mpi::Comm comm = A.Comm();
    const DistGraph& graph = A.LockedDistGraph();
    const int blocksize = A.Blocksize();
    const int commSize = mpi::Size( comm );
    const int numLocal = sepTree.localSepsAndLeaves.size();
    const int numDist = sepTree.distSeps.size();
    DEBUG_ONLY(const int numSources = graph.NumSources())

    // Get the reordered indices of the targets of our portion of the 
    // distributed sparse matrix
    std::set<int> targetSet( graph.targets_.begin(), graph.targets_.end() );
    std::vector<int> targets( targetSet.size() );
    std::copy( targetSet.begin(), targetSet.end(), targets.begin() );
    std::vector<int> mappedTargets = targets;
    reordering.Translate( mappedTargets );

    // Set up the indices for the rows we need from each process
    std::vector<int> recvRowSizes( commSize, 0 );
    for( int s=0; s<numLocal; ++s )
    {
        const SepOrLeaf& sepOrLeaf = *sepTree.localSepsAndLeaves[s];
        const int numInds = sepOrLeaf.inds.size();
        for( int t=0; t<numInds; ++t )
        {
            const int i = sepOrLeaf.inds[t];
            DEBUG_ONLY(
                if( i < 0 || i >= numSources )
                    LogicError("separator index was out of bounds");
            )
            const int q = RowToProcess( i, blocksize, commSize );
            ++recvRowSizes[q];
        }
    }
    for( int s=0; s<numDist; ++s )
    {
        const DistSeparator& sep = sepTree.distSeps[s];
        const DistSymmNodeInfo& node = info.distNodes[s+1];
        const Grid& grid = *node.grid;
        const int rowShift = grid.Col();
        const int rowStride = grid.Width();
        const int numInds = sep.inds.size();
        for( int t=rowShift; t<numInds; t+=rowStride )
        {
            const int i = sep.inds[t];
            DEBUG_ONLY(
                if( i < 0 || i >= numSources )
                    LogicError("separator index was out of bounds");
            )
            const int q = RowToProcess( i, blocksize, commSize );
            ++recvRowSizes[q];
        }
    }
    int numRecvRows=0;
    std::vector<int> recvRowOffs( commSize );
    for( int q=0; q<commSize; ++q )
    {
        recvRowOffs[q] = numRecvRows;
        numRecvRows += recvRowSizes[q];
    }
    std::vector<int> recvRows( numRecvRows );
    std::vector<int> offs = recvRowOffs;
    for( int s=0; s<numLocal; ++s )
    {
        const SepOrLeaf& sepOrLeaf = *sepTree.localSepsAndLeaves[s];
        const int numInds = sepOrLeaf.inds.size();
        for( int t=0; t<numInds; ++t )
        {
            const int i = sepOrLeaf.inds[t];
            DEBUG_ONLY(
                if( i < 0 || i >= numSources )
                    LogicError("separator index was out of bounds");
            )
            const int q = RowToProcess( i, blocksize, commSize );
            DEBUG_ONLY(
                if( offs[q] >= numRecvRows )
                    LogicError("offset got too large");
            )
            recvRows[offs[q]++] = i;
        }
    }
    for( int s=0; s<numDist; ++s )
    {
        const DistSeparator& sep = sepTree.distSeps[s];
        const DistSymmNodeInfo& node = info.distNodes[s+1];
        const Grid& grid = *node.grid;
        const int rowShift = grid.Col();
        const int rowStride = grid.Width();
        const int numInds = sep.inds.size();
        for( int t=rowShift; t<numInds; t+=rowStride )
        {
            const int i = sep.inds[t];
            DEBUG_ONLY(
                if( i < 0 || i >= numSources )
                    LogicError("separator index was out of bounds");
            )
            const int q = RowToProcess( i, blocksize, commSize );
            DEBUG_ONLY(
                if( offs[q] >= numRecvRows )
                    LogicError("offset got too large");
            )
            recvRows[offs[q]++] = i;
        }
    }

    // Retreive the list of rows that we must send to each process
    std::vector<int> sendRowSizes( commSize );
    mpi::AllToAll( &recvRowSizes[0], 1, &sendRowSizes[0], 1, comm );
    int numSendRows=0;
    std::vector<int> sendRowOffs( commSize );
    for( int q=0; q<commSize; ++q )
    {
        sendRowOffs[q] = numSendRows;
        numSendRows += sendRowSizes[q];
    }
    std::vector<int> sendRows( numSendRows );
    mpi::AllToAll
    ( &recvRows[0], &recvRowSizes[0], &recvRowOffs[0],
      &sendRows[0], &sendRowSizes[0], &sendRowOffs[0], comm );

    // Pack the number of nonzeros per row (and the nonzeros themselves)
    // TODO: Avoid sending upper-triangular data
    int numSendEntries=0;
    const int firstLocalRow = A.FirstLocalRow();
    std::vector<int> sendRowLengths( numSendRows ),
                     sendEntriesSizes( commSize, 0 ),
                     sendEntriesOffs( commSize );
    for( int q=0; q<commSize; ++q )
    {
        const int size = sendRowSizes[q];
        const int off = sendRowOffs[q];
        sendEntriesOffs[q] = numSendEntries;
        for( int s=0; s<size; ++s )
        {
            const int i = sendRows[s+off];
            const int iLocal = i - firstLocalRow;
            const int numConnections = A.NumConnections( iLocal );
            numSendEntries += numConnections;
            sendEntriesSizes[q] += numConnections;
            sendRowLengths[s+off] = numConnections;
        }
    }
    std::vector<T> sendEntries( numSendEntries );
    std::vector<int> sendTargets( numSendEntries );
    for( int q=0; q<commSize; ++q )
    {
        int index = sendEntriesOffs[q];
        const int size = sendRowSizes[q];
        const int off = sendRowOffs[q];
        for( int s=0; s<size; ++s )
        {
            const int i = sendRows[s+off];
            const int iLocal = i - firstLocalRow;
            const int numConnections = sendRowLengths[s+off];
            const int localEntryOff = A.LocalEntryOffset( iLocal );
            for( int t=0; t<numConnections; ++t )
            {
                const T value = A.Value( localEntryOff+t );
                const int col = A.Col( localEntryOff+t );
                const int targetOff = Find( targets, col );
                const int mappedTarget = mappedTargets[targetOff];
                DEBUG_ONLY(
                    if( index >= numSendEntries )
                        LogicError("send entry index got too big");
                )
                sendEntries[index] = (conjugate ? Conj(value) : value);
                sendTargets[index] = mappedTarget;
                ++index;
            }
        }
        DEBUG_ONLY(
            if( index != sendEntriesOffs[q]+sendEntriesSizes[q] )
                LogicError("index was not the correct value");
        )
    }

    // Send back the number of nonzeros per row and the nonzeros themselves
    std::vector<int> recvRowLengths( numRecvRows );
    mpi::AllToAll
    ( &sendRowLengths[0], &sendRowSizes[0], &sendRowOffs[0],
      &recvRowLengths[0], &recvRowSizes[0], &recvRowOffs[0], comm );
    int numRecvEntries=0;
    std::vector<int> recvEntriesSizes( commSize, 0 ),
                     recvEntriesOffs( commSize );
    for( int q=0; q<commSize; ++q )
    {
        const int size = recvRowSizes[q];
        const int off = recvRowOffs[q];
        for( int s=0; s<size; ++s )
            recvEntriesSizes[q] += recvRowLengths[off+s];

        recvEntriesOffs[q] = numRecvEntries; 
        numRecvEntries += recvEntriesSizes[q];
    }
    std::vector<T> recvEntries( numRecvEntries );
    std::vector<int> recvTargets( numRecvEntries );
    mpi::AllToAll
    ( &sendEntries[0], &sendEntriesSizes[0], &sendEntriesOffs[0],
      &recvEntries[0], &recvEntriesSizes[0], &recvEntriesOffs[0], comm );
    mpi::AllToAll
    ( &sendTargets[0], &sendEntriesSizes[0], &sendEntriesOffs[0],
      &recvTargets[0], &recvEntriesSizes[0], &recvEntriesOffs[0], comm );

    // Unpack the received entries
    offs = recvRowOffs;
    std::vector<int> entryOffs = recvEntriesOffs;
    localFronts.resize( numLocal );
    for( int s=0; s<numLocal; ++s )
    {
        SymmFront<T>& front = localFronts[s];
        const SepOrLeaf& sepOrLeaf = *sepTree.localSepsAndLeaves[s];
        const SymmNodeInfo& node = info.localNodes[s];
        const std::vector<int>& origLowerStruct = node.origLowerStruct;

        const int size = node.size;
        const int off = node.off;
        const int lowerSize = node.lowerStruct.size();
        Zeros( front.frontL, size+lowerSize, size );
        DEBUG_ONLY(
            if( size != (int)sepOrLeaf.inds.size() )
                LogicError("Mismatch between separator and node size");
        )

        for( int t=0; t<size; ++t )
        {
            const int i = sepOrLeaf.inds[t];
            const int q = RowToProcess( i, blocksize, commSize );

            int& entryOff = entryOffs[q];
            const int numEntries = recvRowLengths[offs[q]++];

            for( int k=0; k<numEntries; ++k )
            {
                const T value = recvEntries[entryOff];
                const int target = recvTargets[entryOff];
                ++entryOff;

                if( target < off+t )
                    continue;
                else if( target < off+size )
                {
                    front.frontL.Set( target-off, t, value );
                }
                else
                {
                    const int origOff = Find( origLowerStruct, target );
                    DEBUG_ONLY(
                        if( origOff >= (int)node.origLowerRelInds.size() )
                            LogicError("origLowerRelInds too small");
                    )
                    const int row = node.origLowerRelInds[origOff];
                    DEBUG_ONLY(
                        if( row < t )
                            LogicError("Tried to touch upper triangle");
                    )
                    front.frontL.Set( row, t, value );
                }
            }
        }
    }

    distFronts.resize( numDist+1 );
    for( int s=0; s<numDist; ++s )
    {
        DistSymmFront<T>& front = distFronts[s+1];
        const DistSeparator& sep = sepTree.distSeps[s];
        const DistSymmNodeInfo& node = info.distNodes[s+1];
        const std::vector<int>& origLowerStruct = node.origLowerStruct;

        const Grid& grid = *node.grid;
        const int colShift = grid.Row();
        const int rowShift = grid.Col();
        const int colStride = grid.Height();
        const int rowStride = grid.Width();

        const int size = node.size;
        const int off = node.off;
        const int lowerSize = node.lowerStruct.size();
        front.front2dL.SetGrid( grid );
        Zeros( front.front2dL, size+lowerSize, size );
        DEBUG_ONLY(
            if( size != (int)sep.inds.size() )
                LogicError("Mismatch in separator and node sizes");
        )

        for( int t=rowShift; t<size; t+=rowStride )
        {
            const int i = sep.inds[t];
            const int q = RowToProcess( i, blocksize, commSize );
            const int localCol = (t-rowShift) / rowStride;

            int& entryOff = entryOffs[q];
            const int numEntries = recvRowLengths[offs[q]++];

            for( int k=0; k<numEntries; ++k )
            {
                const T value = recvEntries[entryOff];
                const int target = recvTargets[entryOff];
                ++entryOff;

                if( target < off+t )
                    continue;
                else if( target < off+size )
                {
                    if( (target-off) % colStride == colShift )
                    {
                        const int row = target-off;
                        const int localRow = (row-colShift) / colStride;
                        front.front2dL.SetLocal( localRow, localCol, value );
                    }
                }
                else 
                {
                    const int origOff = Find( origLowerStruct, target );
                    DEBUG_ONLY(
                        if( origOff >= (int)node.origLowerRelInds.size() )
                            LogicError("origLowerRelInds too small");
                    )
                    const int row = node.origLowerRelInds[origOff];
                    DEBUG_ONLY(
                        if( row < t )
                            LogicError("Tried to touch upper triangle");
                    )
                    if( row % colStride == colShift )
                    {
                        const int localRow = (row-colShift) / colStride;
                        front.front2dL.SetLocal( localRow, localCol, value );
                    }
                }
            }
        }
    }
    DEBUG_ONLY(
        for( int q=0; q<commSize; ++q )
            if( entryOffs[q] != recvEntriesOffs[q]+recvEntriesSizes[q] )
                LogicError("entryOffs were incorrect");
    )
    
    // Copy information from the local root to the dist leaf
    {
        const DistSymmNodeInfo& node = info.distNodes[0];
        Matrix<T>& topLocal = localFronts.back().frontL;
        DistMatrix<T>& bottomDist = distFronts[0].front2dL;
        bottomDist.LockedAttach
        ( topLocal.Height(), topLocal.Width(), *node.grid, 0, 0, topLocal );
    }
}

template<typename T>
DistSymmFrontTree<T>::DistSymmFrontTree
( const DistSparseMatrix<T>& A, 
  const DistMap& reordering,
  const DistSeparatorTree& sepTree, 
  const DistSymmInfo& info,
  bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("DistSymmFrontTree::DistSymmFrontTree"))
    Initialize( A, reordering, sepTree, info, conjugate );
}

template<typename T>
void DistSymmFrontTree<T>::MemoryInfo
( double& numLocalEntries, double& minLocalEntries, double& maxLocalEntries, 
  double& numGlobalEntries ) const
{
    DEBUG_ONLY(CallStackEntry cse("DistSymmFrontTree::MemoryInfo"))
    numLocalEntries = numGlobalEntries = 0;
    const int numLocalFronts = localFronts.size();
    const int numDistFronts = distFronts.size();
    const bool frontsAre1d = FrontsAre1d( frontType );
    const Grid& grid = ( frontsAre1d ? distFronts.back().front1dL.Grid() 
                                     : distFronts.back().front2dL.Grid() );
    mpi::Comm comm = grid.Comm();

    for( int s=0; s<numLocalFronts; ++s )
    {
        const SymmFront<T>& front = localFronts[s];
        numLocalEntries += front.frontL.MemorySize();
        numLocalEntries += front.diag.MemorySize();
        numLocalEntries += front.subdiag.MemorySize();
        numLocalEntries += front.piv.MemorySize();
        numLocalEntries += front.work.MemorySize();
    }
    for( int s=1; s<numDistFronts; ++s )
    {
        const DistSymmFront<T>& front = distFronts[s];
        numLocalEntries += front.front1dL.AllocatedMemory();
        numLocalEntries += front.front2dL.AllocatedMemory();
        numLocalEntries += front.diag1d.AllocatedMemory();
        numLocalEntries += front.subdiag1d.AllocatedMemory();
        numLocalEntries += front.piv.AllocatedMemory();
        numLocalEntries += front.work1d.AllocatedMemory();
        numLocalEntries += front.work2d.AllocatedMemory();
    }

    minLocalEntries = mpi::AllReduce( numLocalEntries, mpi::MIN, comm );
    maxLocalEntries = mpi::AllReduce( numLocalEntries, mpi::MAX, comm );
    numGlobalEntries = mpi::AllReduce( numLocalEntries, mpi::SUM, comm );
}

template<typename T>
void DistSymmFrontTree<T>::TopLeftMemoryInfo
( double& numLocalEntries, double& minLocalEntries, double& maxLocalEntries, 
  double& numGlobalEntries ) const
{
    DEBUG_ONLY(CallStackEntry cse("DistSymmFrontTree::TopLeftMemInfo"))
    numLocalEntries = numGlobalEntries = 0;
    const int numLocalFronts = localFronts.size();
    const int numDistFronts = distFronts.size();
    const bool frontsAre1d = FrontsAre1d( frontType );
    const Grid& grid = ( frontsAre1d ? distFronts.back().front1dL.Grid() 
                                     : distFronts.back().front2dL.Grid() );
    mpi::Comm comm = grid.Comm();

    for( int s=0; s<numLocalFronts; ++s )
    {
        const SymmFront<T>& front = localFronts[s];
        Matrix<T> FTL, FBL;
        LockedPartitionDown( front.frontL, FTL, FBL, front.frontL.Width() );
        numLocalEntries += FTL.Height()*FTL.Width();
        numLocalEntries += front.diag.MemorySize();
        numLocalEntries += front.subdiag.MemorySize();
        numLocalEntries += front.piv.MemorySize();
    }
    for( int s=1; s<numDistFronts; ++s )
    {
        const DistSymmFront<T>& front = distFronts[s];
        if( frontsAre1d )
        {
            DistMatrix<T,VC,STAR> FTL(grid), FBL(grid);
            LockedPartitionDown
            ( front.front1dL, FTL, FBL, front.front1dL.Width() );
            numLocalEntries += FTL.LocalHeight()*FTL.LocalWidth();
        }
        else
        {
            DistMatrix<T> FTL(grid), FBL(grid);
            LockedPartitionDown
            ( front.front2dL, FTL, FBL, front.front2dL.Width() );
            numLocalEntries += FTL.LocalHeight()*FTL.LocalWidth();
        }
        numLocalEntries += front.diag1d.AllocatedMemory();
        numLocalEntries += front.subdiag1d.AllocatedMemory();
        numLocalEntries += front.piv.AllocatedMemory();
    }

    minLocalEntries = mpi::AllReduce( numLocalEntries, mpi::MIN, comm );
    maxLocalEntries = mpi::AllReduce( numLocalEntries, mpi::MAX, comm );
    numGlobalEntries = mpi::AllReduce( numLocalEntries, mpi::SUM, comm );
}

template<typename T>
void DistSymmFrontTree<T>::BottomLeftMemoryInfo
( double& numLocalEntries, double& minLocalEntries, double& maxLocalEntries, 
  double& numGlobalEntries ) const
{
    DEBUG_ONLY(CallStackEntry cse("DistSymmFrontTree::BottomLeftMemInfo"))
    numLocalEntries = numGlobalEntries = 0;
    const int numLocalFronts = localFronts.size();
    const int numDistFronts = distFronts.size();
    const bool frontsAre1d = FrontsAre1d( frontType );
    const Grid& grid = ( frontsAre1d ? distFronts.back().front1dL.Grid() 
                                     : distFronts.back().front2dL.Grid() );
    mpi::Comm comm = grid.Comm();

    for( int s=0; s<numLocalFronts; ++s )
    {
        const SymmFront<T>& front = localFronts[s];
        Matrix<T> FTL, FBL;
        LockedPartitionDown( front.frontL, FTL, FBL, front.frontL.Width() );
        numLocalEntries += FBL.Height()*FBL.Width();
    }
    for( int s=1; s<numDistFronts; ++s )
    {
        const DistSymmFront<T>& front = distFronts[s];
        if( frontsAre1d )
        {
            DistMatrix<T,VC,STAR> FTL(grid), FBL(grid);
            LockedPartitionDown
            ( front.front1dL, FTL, FBL, front.front1dL.Width() );
            numLocalEntries += FBL.LocalHeight()*FBL.LocalWidth();
        }
        else
        {
            DistMatrix<T> FTL(grid), FBL(grid);
            LockedPartitionDown
            ( front.front2dL, FTL, FBL, front.front2dL.Width() );
            numLocalEntries += FBL.LocalHeight()*FBL.LocalWidth();
        }
    }

    minLocalEntries = mpi::AllReduce( numLocalEntries, mpi::MIN, comm );
    maxLocalEntries = mpi::AllReduce( numLocalEntries, mpi::MAX, comm );
    numGlobalEntries = mpi::AllReduce( numLocalEntries, mpi::SUM, comm );
}

template<typename T>
void DistSymmFrontTree<T>::FactorizationWork
( double& numLocalFlops, double& minLocalFlops, double& maxLocalFlops, 
  double& numGlobalFlops, bool selInv ) const
{
    DEBUG_ONLY(CallStackEntry cse("DistSymmFrontTree::FactorizationWork"))
    numLocalFlops = numGlobalFlops = 0;
    const int numLocalFronts = localFronts.size();
    const int numDistFronts = distFronts.size();
    const bool frontsAre1d = FrontsAre1d( frontType );
    const Grid& grid = ( frontsAre1d ? distFronts.back().front1dL.Grid() 
                                     : distFronts.back().front2dL.Grid() );
    mpi::Comm comm = grid.Comm();

    for( int s=0; s<numLocalFronts; ++s )
    {
        const SymmFront<T>& front = localFronts[s];
        const double m = front.frontL.Height();
        const double n = front.frontL.Width();
        numLocalFlops += (1./3.)*n*n*n; // n x n LDL
        numLocalFlops += (m-n)*n*n; // n x n trsv, m-n r.h.s.
        numLocalFlops += (m-n)*(m-n)*n; // (m-n) x (m-n), rank-n
    }
    for( int s=1; s<numDistFronts; ++s )
    {
        const DistSymmFront<T>& front = distFronts[s];
        const double m = 
          ( frontsAre1d ? front.front1dL.Height() : front.front2dL.Height() );
        const double n = 
          ( frontsAre1d ? front.front1dL.Width() : front.front2dL.Width() );
        const double pFront = 
          ( frontsAre1d ? front.front1dL.Grid().Size()
                        : front.front2dL.Grid().Size() );
        numLocalFlops += (1./3.)*n*n*n/pFront;
        numLocalFlops += (m-n)*n*n/pFront;
        numLocalFlops += (m-n)*(m-n)*n/pFront; 
        if( selInv )
            numLocalFlops += (1./3.)*n*n*n/pFront;
    }

    // Since there are equal numbers of multiplies and adds, and the former
    // takes 6 times as much work in standard complex arithmetic, while the 
    // later only takes twice, the average is 4x more work
    if( IsComplex<T>::val )
        numLocalFlops *= 4;

    minLocalFlops = mpi::AllReduce( numLocalFlops, mpi::MIN, comm );
    maxLocalFlops = mpi::AllReduce( numLocalFlops, mpi::MAX, comm );
    numGlobalFlops = mpi::AllReduce( numLocalFlops, mpi::SUM, comm );
}

template<typename T>
void DistSymmFrontTree<T>::SolveWork
( double& numLocalFlops, double& minLocalFlops, double& maxLocalFlops, 
  double& numGlobalFlops, int numRhs ) const
{
    DEBUG_ONLY(CallStackEntry cse("DistSymmFrontTree::SolveWork"))
    numLocalFlops = numGlobalFlops = 0;
    const int numLocalFronts = localFronts.size();
    const int numDistFronts = distFronts.size();
    const bool frontsAre1d = FrontsAre1d( frontType );
    const Grid& grid = ( frontsAre1d ? distFronts.back().front1dL.Grid() 
                                     : distFronts.back().front2dL.Grid() );
    mpi::Comm comm = grid.Comm();

    for( int s=0; s<numLocalFronts; ++s )
    {
        const SymmFront<T>& front = localFronts[s];
        const double m = front.frontL.Height();
        const double n = front.frontL.Width();
        numLocalFlops += n*n;
        numLocalFlops += 2*(m-n)*n;
    }
    for( int s=1; s<numDistFronts; ++s )
    {
        const DistSymmFront<T>& front = distFronts[s];
        const double m = 
          ( frontsAre1d ? front.front1dL.Height() : front.front2dL.Height() );
        const double n = 
          ( frontsAre1d ? front.front1dL.Width() : front.front2dL.Width() );
        const double pFront = 
          ( frontsAre1d ? front.front1dL.Grid().Size()
                        : front.front2dL.Grid().Size() );
        numLocalFlops += n*n/pFront;
        numLocalFlops += 2*(m-n)*n/pFront;
    }

    // Since there are equal numbers of multiplies and adds, and the former
    // takes 6 times as much work in standard complex arithmetic, while the 
    // later only takes twice, the average is 4x more work
    if( IsComplex<T>::val )
        numLocalFlops *= 4;

    numLocalFlops *= numRhs;
    minLocalFlops = mpi::AllReduce( numLocalFlops, mpi::MIN, comm );
    maxLocalFlops = mpi::AllReduce( numLocalFlops, mpi::MAX, comm );
    numGlobalFlops = mpi::AllReduce( numLocalFlops, mpi::SUM, comm );
}

#define PROTO(T) template class DistSymmFrontTree<T>;
#include "El/macros/Instantiate.h"

} // namespace El
