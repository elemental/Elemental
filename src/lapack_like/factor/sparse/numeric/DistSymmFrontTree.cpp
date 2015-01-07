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
    const Int commSize = mpi::Size( comm );
    const Int numLocal = sepTree.localSepsAndLeaves.size();
    const Int numDist = sepTree.distSeps.size();

    // Get the reordered indices of the targets of our portion of the 
    // distributed sparse matrix
    std::set<Int> targetSet( graph.targets_.begin(), graph.targets_.end() );
    std::vector<Int> targets( targetSet.size() );
    std::copy( targetSet.begin(), targetSet.end(), targets.begin() );
    std::vector<Int> mappedTargets = targets;
    reordering.Translate( mappedTargets );

    // Set up the indices for the rows we need from each process
    std::vector<int> recvRowSizes( commSize, 0 );
    for( Int s=0; s<numLocal; ++s )
    {
        const SepOrLeaf& sepOrLeaf = *sepTree.localSepsAndLeaves[s];
        const Int numInds = sepOrLeaf.inds.size();
        for( Int t=0; t<numInds; ++t )
            ++recvRowSizes[ A.RowOwner(sepOrLeaf.inds[t]) ];
    }
    for( Int s=0; s<numDist; ++s )
    {
        const DistSeparator& sep = sepTree.distSeps[s];
        const DistSymmNodeInfo& node = info.distNodes[s+1];
        const Grid& grid = *node.grid;
        const Int rowShift = grid.Col();
        const Int rowStride = grid.Width();
        const Int numInds = sep.inds.size();
        for( Int t=rowShift; t<numInds; t+=rowStride )
            ++recvRowSizes[ A.RowOwner(sep.inds[t]) ];
    }
    std::vector<int> recvRowOffs;
    const Int numRecvRows = Scan( recvRowSizes, recvRowOffs );

    std::vector<Int> recvRows( numRecvRows );
    auto offs = recvRowOffs;
    for( Int s=0; s<numLocal; ++s )
    {
        const SepOrLeaf& sepOrLeaf = *sepTree.localSepsAndLeaves[s];
        const Int numInds = sepOrLeaf.inds.size();
        for( Int t=0; t<numInds; ++t )
        {
            const Int i = sepOrLeaf.inds[t];
            const Int q = A.RowOwner(i);
            recvRows[offs[q]++] = i;
        }
    }
    for( Int s=0; s<numDist; ++s )
    {
        const DistSeparator& sep = sepTree.distSeps[s];
        const DistSymmNodeInfo& node = info.distNodes[s+1];
        const Grid& grid = *node.grid;
        const Int rowShift = grid.Col();
        const Int rowStride = grid.Width();
        const Int numInds = sep.inds.size();
        for( Int t=rowShift; t<numInds; t+=rowStride )
        {
            const Int i = sep.inds[t];
            const Int q = A.RowOwner(i);
            recvRows[offs[q]++] = i;
        }
    }

    // Retreive the list of rows that we must send to each process
    std::vector<int> sendRowSizes( commSize );
    mpi::AllToAll( &recvRowSizes[0], 1, &sendRowSizes[0], 1, comm );
    std::vector<int> sendRowOffs;
    const Int numSendRows = Scan( sendRowSizes, sendRowOffs );
    std::vector<Int> sendRows( numSendRows );
    mpi::AllToAll
    ( &recvRows[0], &recvRowSizes[0], &recvRowOffs[0],
      &sendRows[0], &sendRowSizes[0], &sendRowOffs[0], comm );

    // Pack the number of nonzeros per row (and the nonzeros themselves)
    // TODO: Avoid sending upper-triangular data
    Int numSendEntries=0;
    const Int firstLocalRow = A.FirstLocalRow();
    std::vector<Int> sendRowLengths( numSendRows );
    std::vector<int> sendEntriesSizes( commSize, 0 ), 
                     sendEntriesOffs( commSize );
    for( Int q=0; q<commSize; ++q )
    {
        const Int size = sendRowSizes[q];
        const Int off = sendRowOffs[q];
        sendEntriesOffs[q] = numSendEntries;
        for( Int s=0; s<size; ++s )
        {
            const Int i = sendRows[s+off];
            const Int iLocal = i - firstLocalRow;
            const Int numConnections = A.NumConnections( iLocal );
            numSendEntries += numConnections;
            sendEntriesSizes[q] += numConnections;
            sendRowLengths[s+off] = numConnections;
        }
    }
    std::vector<T> sendEntries( numSendEntries );
    std::vector<Int> sendTargets( numSendEntries );
    for( Int q=0; q<commSize; ++q )
    {
        Int index = sendEntriesOffs[q];
        const Int size = sendRowSizes[q];
        const Int off = sendRowOffs[q];
        for( Int s=0; s<size; ++s )
        {
            const Int i = sendRows[s+off];
            const Int iLocal = i - firstLocalRow;
            const Int numConnections = sendRowLengths[s+off];
            const Int localEntryOff = A.EntryOffset( iLocal );
            for( Int t=0; t<numConnections; ++t )
            {
                const T value = A.Value( localEntryOff+t );
                const Int col = A.Col( localEntryOff+t );
                const Int targetOff = Find( targets, col );
                const Int mappedTarget = mappedTargets[targetOff];
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
    std::vector<Int> recvRowLengths( numRecvRows );
    mpi::AllToAll
    ( sendRowLengths.data(), sendRowSizes.data(), sendRowOffs.data(),
      recvRowLengths.data(), recvRowSizes.data(), recvRowOffs.data(), comm );
    Int numRecvEntries=0;
    std::vector<int> recvEntriesSizes( commSize, 0 ),
                     recvEntriesOffs( commSize );
    for( Int q=0; q<commSize; ++q )
    {
        const Int size = recvRowSizes[q];
        const Int off = recvRowOffs[q];
        for( Int s=0; s<size; ++s )
            recvEntriesSizes[q] += recvRowLengths[off+s];

        recvEntriesOffs[q] = numRecvEntries; 
        numRecvEntries += recvEntriesSizes[q];
    }
    std::vector<T> recvEntries( numRecvEntries );
    std::vector<Int> recvTargets( numRecvEntries );
    mpi::AllToAll
    ( sendEntries.data(), sendEntriesSizes.data(), sendEntriesOffs.data(),
      recvEntries.data(), recvEntriesSizes.data(), recvEntriesOffs.data(), 
      comm );
    mpi::AllToAll
    ( sendTargets.data(), sendEntriesSizes.data(), sendEntriesOffs.data(),
      recvTargets.data(), recvEntriesSizes.data(), recvEntriesOffs.data(), 
      comm );

    // Unpack the received entries
    offs = recvRowOffs;
    auto entryOffs = recvEntriesOffs;
    localFronts.resize( numLocal );
    for( Int s=0; s<numLocal; ++s )
    {
        SymmFront<T>& front = localFronts[s];
        const SepOrLeaf& sepOrLeaf = *sepTree.localSepsAndLeaves[s];
        const SymmNodeInfo& node = info.localNodes[s];
        const std::vector<Int>& origLowerStruct = node.origLowerStruct;

        const Int size = node.size;
        const Int off = node.off;
        const Int lowerSize = node.lowerStruct.size();
        Zeros( front.frontL, size+lowerSize, size );
        DEBUG_ONLY(
            if( size != (int)sepOrLeaf.inds.size() )
                LogicError("Mismatch between separator and node size");
        )

        for( Int t=0; t<size; ++t )
        {
            const Int i = sepOrLeaf.inds[t];
            const Int q = A.RowOwner(i);

            int& entryOff = entryOffs[q];
            const Int numEntries = recvRowLengths[offs[q]++];

            for( Int k=0; k<numEntries; ++k )
            {
                const T value = recvEntries[entryOff];
                const Int target = recvTargets[entryOff];
                ++entryOff;

                if( target < off+t )
                    continue;
                else if( target < off+size )
                {
                    front.frontL.Set( target-off, t, value );
                }
                else
                {
                    const Int origOff = Find( origLowerStruct, target );
                    const Int row = node.origLowerRelInds[origOff];
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
    for( Int s=0; s<numDist; ++s )
    {
        DistSymmFront<T>& front = distFronts[s+1];
        const DistSeparator& sep = sepTree.distSeps[s];
        const DistSymmNodeInfo& node = info.distNodes[s+1];
        const std::vector<Int>& origLowerStruct = node.origLowerStruct;

        const Grid& grid = *node.grid;
        const Int colShift = grid.Row();
        const Int rowShift = grid.Col();
        const Int colStride = grid.Height();
        const Int rowStride = grid.Width();

        const Int size = node.size;
        const Int off = node.off;
        const Int lowerSize = node.lowerStruct.size();
        front.front2dL.SetGrid( grid );
        Zeros( front.front2dL, size+lowerSize, size );
        DEBUG_ONLY(
            if( size != (int)sep.inds.size() )
                LogicError("Mismatch in separator and node sizes");
        )

        for( Int t=rowShift; t<size; t+=rowStride )
        {
            const Int i = sep.inds[t];
            const Int q = A.RowOwner(i);
            const Int localCol = (t-rowShift) / rowStride;

            Int& entryOff = entryOffs[q];
            const Int numEntries = recvRowLengths[offs[q]++];

            for( Int k=0; k<numEntries; ++k )
            {
                const T value = recvEntries[entryOff];
                const Int target = recvTargets[entryOff];
                ++entryOff;

                if( target < off+t )
                    continue;
                else if( target < off+size )
                {
                    if( (target-off) % colStride == colShift )
                    {
                        const Int row = target-off;
                        const Int localRow = (row-colShift) / colStride;
                        front.front2dL.SetLocal( localRow, localCol, value );
                    }
                }
                else 
                {
                    const Int origOff = Find( origLowerStruct, target );
                    const Int row = node.origLowerRelInds[origOff];
                    DEBUG_ONLY(
                        if( row < t )
                            LogicError("Tried to touch upper triangle");
                    )
                    if( row % colStride == colShift )
                    {
                        const Int localRow = (row-colShift) / colStride;
                        front.front2dL.SetLocal( localRow, localCol, value );
                    }
                }
            }
        }
    }
    DEBUG_ONLY(
        for( Int q=0; q<commSize; ++q )
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
    const Int numLocalFronts = localFronts.size();
    const Int numDistFronts = distFronts.size();
    const bool frontsAre1d = FrontsAre1d( frontType );
    const Grid& grid = ( frontsAre1d ? distFronts.back().front1dL.Grid() 
                                     : distFronts.back().front2dL.Grid() );
    mpi::Comm comm = grid.Comm();

    for( Int s=0; s<numLocalFronts; ++s )
    {
        const SymmFront<T>& front = localFronts[s];
        numLocalEntries += front.frontL.MemorySize();
        numLocalEntries += front.diag.MemorySize();
        numLocalEntries += front.subdiag.MemorySize();
        numLocalEntries += front.piv.MemorySize();
        numLocalEntries += front.work.MemorySize();
    }
    for( Int s=1; s<numDistFronts; ++s )
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
    const Int numLocalFronts = localFronts.size();
    const Int numDistFronts = distFronts.size();
    const bool frontsAre1d = FrontsAre1d( frontType );
    const Grid& grid = ( frontsAre1d ? distFronts.back().front1dL.Grid() 
                                     : distFronts.back().front2dL.Grid() );
    mpi::Comm comm = grid.Comm();

    for( Int s=0; s<numLocalFronts; ++s )
    {
        const SymmFront<T>& front = localFronts[s];
        Matrix<T> FTL, FBL;
        LockedPartitionDown( front.frontL, FTL, FBL, front.frontL.Width() );
        numLocalEntries += FTL.Height()*FTL.Width();
        numLocalEntries += front.diag.MemorySize();
        numLocalEntries += front.subdiag.MemorySize();
        numLocalEntries += front.piv.MemorySize();
    }
    for( Int s=1; s<numDistFronts; ++s )
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
    const Int numLocalFronts = localFronts.size();
    const Int numDistFronts = distFronts.size();
    const bool frontsAre1d = FrontsAre1d( frontType );
    const Grid& grid = ( frontsAre1d ? distFronts.back().front1dL.Grid() 
                                     : distFronts.back().front2dL.Grid() );
    mpi::Comm comm = grid.Comm();

    for( Int s=0; s<numLocalFronts; ++s )
    {
        const SymmFront<T>& front = localFronts[s];
        Matrix<T> FTL, FBL;
        LockedPartitionDown( front.frontL, FTL, FBL, front.frontL.Width() );
        numLocalEntries += FBL.Height()*FBL.Width();
    }
    for( Int s=1; s<numDistFronts; ++s )
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
    const Int numLocalFronts = localFronts.size();
    const Int numDistFronts = distFronts.size();
    const bool frontsAre1d = FrontsAre1d( frontType );
    const Grid& grid = ( frontsAre1d ? distFronts.back().front1dL.Grid() 
                                     : distFronts.back().front2dL.Grid() );
    mpi::Comm comm = grid.Comm();

    for( Int s=0; s<numLocalFronts; ++s )
    {
        const SymmFront<T>& front = localFronts[s];
        const double m = front.frontL.Height();
        const double n = front.frontL.Width();
        numLocalFlops += (1./3.)*n*n*n; // n x n LDL
        numLocalFlops += (m-n)*n*n; // n x n trsv, m-n r.h.s.
        numLocalFlops += (m-n)*(m-n)*n; // (m-n) x (m-n), rank-n
    }
    for( Int s=1; s<numDistFronts; ++s )
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
  double& numGlobalFlops, Int numRhs ) const
{
    DEBUG_ONLY(CallStackEntry cse("DistSymmFrontTree::SolveWork"))
    numLocalFlops = numGlobalFlops = 0;
    const Int numLocalFronts = localFronts.size();
    const Int numDistFronts = distFronts.size();
    const bool frontsAre1d = FrontsAre1d( frontType );
    const Grid& grid = ( frontsAre1d ? distFronts.back().front1dL.Grid() 
                                     : distFronts.back().front2dL.Grid() );
    mpi::Comm comm = grid.Comm();

    for( Int s=0; s<numLocalFronts; ++s )
    {
        const SymmFront<T>& front = localFronts[s];
        const double m = front.frontL.Height();
        const double n = front.frontL.Width();
        numLocalFlops += n*n;
        numLocalFlops += 2*(m-n)*n;
    }
    for( Int s=1; s<numDistFronts; ++s )
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
