/*
   Copyright 2009-2011, Jack Poulson.
   All rights reserved.

   Copyright 2011-2012, Jack Poulson, Lexing Ying, and
   The University of Texas at Austin.
   All rights reserved.

   Copyright 2013, Jack Poulson, Lexing Ying, and Stanford University.
   All rights reserved.

   Copyright 2013-2014, Jack Poulson and The Georgia Institute of Technology.
   All rights reserved.

   Copyright 2014-2015, Jack Poulson and Stanford University.
   All rights reserved.

   Copyright 2016, Jack Poulson.
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {
namespace ldl {

template<typename Field>
DistFront<Field>::DistFront( DistFront<Field>* parentNode )
: parent(parentNode)
{
    if( parentNode != nullptr )
    {
        type = parentNode->type;
        isHermitian = parentNode->isHermitian;
    }
}

template<typename Field>
DistFront<Field>::DistFront
( const DistSparseMatrix<Field>& A,
  const DistMap& reordering,
  const DistSeparator& sep,
  const DistNodeInfo& info,
  bool conjugate )
{
    EL_DEBUG_CSE
    Pull( A, reordering, sep, info, conjugate );
}

template<typename Field>
DistFront<Field>::~DistFront() { }

template<typename Field>
void UnpackEntriesLocal
( const Separator& sep,
  const NodeInfo& node,
        Front<Field>& front,
  const DistSparseMatrix<Field>& A,
  const vector<Int>& rRowLengths,
  const vector<Field>& rEntries,
  const vector<Int>& rTargets,
        vector<int>& offs,
        vector<int>& entryOffs )
{
    EL_DEBUG_CSE

    // Delete any existing children
    SwapClear( front.children );

    const Int numChildren = sep.children.size();
    front.children.resize( numChildren );
    for( Int c=0; c<numChildren; ++c )
    {
        front.children[c].reset( new Front<Field>(&front) );
        UnpackEntriesLocal
        ( *sep.children[c], *node.children[c], *front.children[c],
          A, rRowLengths, rEntries, rTargets, offs, entryOffs );
    }
    // Mark this node as a sparse leaf if it does not have any children
    // and is not a duplicate of a dense distributed node
    if( numChildren == 0 && !front.duplicate )
        front.sparseLeaf = true;

    const Int size = node.size;
    const Int off = node.off;
    const Int lowerSize = node.lowerStruct.size();

    if( front.sparseLeaf )
    {
        front.workSparse.Empty();
        Zeros( front.workSparse, size, size );
        Zeros( front.LDense, lowerSize, size );

        Int numSparseEntries = 0;
        auto offsCopy = offs;
        auto entryOffsCopy = entryOffs;
        for( Int t=0; t<size; ++t )
        {
            const Int i = sep.inds[t];
            const Int q = A.RowOwner(i);

            int& entryOff = entryOffsCopy[q];
            const Int numEntries = rRowLengths[offsCopy[q]++];
            for( Int k=0; k<numEntries; ++k )
            {
                const Int target = rTargets[entryOff++];
                EL_DEBUG_ONLY(
                  if( target < off+t )
                      LogicError("Received entry from upper triangle");
                )
                if( target < off+size )
                    ++numSparseEntries;
            }
        }
        front.workSparse.Reserve( numSparseEntries );

        for( Int t=0; t<size; ++t )
        {
            const Int i = sep.inds[t];
            const Int q = A.RowOwner(i);

            int& entryOff = entryOffs[q];
            const Int numEntries = rRowLengths[offs[q]++];

            for( Int k=0; k<numEntries; ++k )
            {
                const Field value = rEntries[entryOff];
                const Int target = rTargets[entryOff];
                ++entryOff;

                EL_DEBUG_ONLY(
                  if( target < off+t )
                      LogicError("Received entry from upper triangle");
                )
                if( target < off+size )
                {
                    const Field transVal =
                      front.isHermitian ? Conj(value) : value;
                    front.workSparse.QueueUpdate( target-off, t, transVal );
                }
                else
                {
                    // TODO(poulson): Avoid this binary search?
                    Int origOff = Find( node.origLowerStruct, target );
                    const Int row = node.origLowerRelInds[origOff];
                    front.LDense(row-size,t) = value;
                }
            }
        }
        front.workSparse.ProcessQueues();
        MakeSymmetric( LOWER, front.workSparse, front.isHermitian );
    }
    else
    {
        Zeros( front.LDense, size+lowerSize, size );
        for( Int t=0; t<size; ++t )
        {
            const Int i = sep.inds[t];
            const Int q = A.RowOwner(i);

            int& entryOff = entryOffs[q];
            const Int numEntries = rRowLengths[offs[q]++];

            for( Int k=0; k<numEntries; ++k )
            {
                const Field value = rEntries[entryOff];
                const Int target = rTargets[entryOff];
                ++entryOff;

                EL_DEBUG_ONLY(
                  if( target < off+t )
                      LogicError("Received entry from upper triangle");
                )
                if( target < off+size )
                {
                    front.LDense(target-off,t) = value;
                }
                else
                {
                    // TODO(poulson): Avoid this binary search?
                    Int origOff = Find( node.origLowerStruct, target );
                    const Int row = node.origLowerRelInds[origOff];
                    front.LDense(row,t) = value;
                }
            }
        }
    }
}

template<typename Field>
void UnpackEntries
( const DistSeparator& sep,
  const DistNodeInfo& node,
        DistFront<Field>& front,
  const DistSparseMatrix<Field>& A,
  const vector<Int>& rRowLengths,
  const vector<Field>& rEntries,
  const vector<Int>& rTargets,
        vector<int>& offs,
        vector<int>& entryOffs )
{
    EL_DEBUG_CSE
    const Grid& grid = node.Grid();

    if( sep.child == nullptr )
    {
        front.duplicate.reset( new Front<Field>(&front) );
        UnpackEntriesLocal
        ( *sep.duplicate, *node.duplicate, *front.duplicate,
          A, rRowLengths, rEntries, rTargets, offs, entryOffs );

        front.L2D.Attach( grid, front.duplicate->LDense );

        return;
    }
    front.child.reset( new DistFront<Field>(&front) );
    UnpackEntries
    ( *sep.child, *node.child, *front.child,
      A, rRowLengths, rEntries, rTargets, offs, entryOffs );

    const Int size = node.size;
    const Int off = node.off;
    const Int lowerSize = node.lowerStruct.size();
    front.L2D.SetGrid( grid );
    Zeros( front.L2D, size+lowerSize, size );

    const Int localWidth = front.L2D.LocalWidth();
    for( Int tLoc=0; tLoc<localWidth; ++tLoc )
    {
        const Int t = front.L2D.GlobalCol(tLoc);
        const Int i = sep.inds[t];
        const Int q = A.RowOwner(i);

        int& entryOff = entryOffs[q];
        const Int numEntries = rRowLengths[offs[q]++];

        for( Int k=0; k<numEntries; ++k )
        {
            const Field value = rEntries[entryOff];
            const Int target = rTargets[entryOff];
            ++entryOff;

            EL_DEBUG_ONLY(
              if( target < off+t )
                  LogicError("Received entry from upper triangle");
            )
            if( target < off+size )
            {
                front.L2D.Set( target-off, t, value );
            }
            else
            {
                // TODO(poulson): Avoid this binary search?
                const Int origOff = Find( node.origLowerStruct, target );
                const Int row = node.origLowerRelInds[origOff];
                front.L2D.Set( row, t, value );
            }
        }
    }
}

// NOTE:
// The current implementation (conjugate-)transposes A into the frontal tree
template<typename Field>
void DistFront<Field>::Pull
( const DistSparseMatrix<Field>& A,
  const DistMap& reordering,
  const DistSeparator& rootSep,
  const DistNodeInfo& rootInfo,
  bool conjugate )
{
    EL_DEBUG_CSE
    vector<Int> mappedSources, mappedTargets, colOffs;
    Pull
    ( A, reordering, rootSep, rootInfo,
      mappedSources, mappedTargets, colOffs,
      conjugate );
}

template<typename Field>
void DistFront<Field>::Pull
( const DistSparseMatrix<Field>& A,
  const DistMap& reordering,
  const DistSeparator& rootSep,
  const DistNodeInfo& rootInfo,
        vector<Int>& mappedSources,
        vector<Int>& mappedTargets,
        vector<Int>& colOffs,
  bool conjugate )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( A.LocalHeight() != reordering.NumLocalSources() )
          LogicError("Local mapping was not the right size");
    )
    const bool time = false;
    const Grid& grid = A.Grid();
    const int commRank = grid.Rank();
    const int commSize = grid.Size();
    Timer timer;

    A.MappedSources( reordering, mappedSources );
    A.MappedTargets( reordering, mappedTargets, colOffs );

    // Set up the indices for the rows we need from each process
    if( time && commRank == 0 )
        timer.Start();
    vector<int> rRowSizes( commSize, 0 );
    function<void(const Separator&)> rRowLocalAccumulate =
      [&]( const Separator& sep )
      {
          for( const auto& child : sep.children )
             rRowLocalAccumulate( *child );
          for( const Int& i : sep.inds )
              ++rRowSizes[ A.RowOwner(i) ];
      };
    function<void(const DistSeparator&,const DistNodeInfo&)>
      rRowAccumulate =
      [&]( const DistSeparator& sep, const DistNodeInfo& node )
      {
          if( sep.child == nullptr )
          {
              rRowLocalAccumulate( *sep.duplicate );
              return;
          }
          rRowAccumulate( *sep.child, *node.child );

          const Grid& grid = node.Grid();
          const Int rowShift = grid.Col();
          const Int rowStride = grid.Width();
          const Int numInds = sep.inds.size();
          for( Int t=rowShift; t<numInds; t+=rowStride )
              ++rRowSizes[ A.RowOwner(sep.inds[t]) ];
      };
    rRowAccumulate( rootSep, rootInfo );
    vector<int> rRowOffs;
    const Int numRecvRows = Scan( rRowSizes, rRowOffs );
    if( time && commRank == 0 )
        Output("Row index setup: ",timer.Stop()," secs");

    if( time && commRank == 0 )
        timer.Start();
    vector<Int> rRows( numRecvRows );
    auto offs = rRowOffs;
    function<void(const Separator&)> rRowsLocalPack =
      [&]( const Separator& sep )
      {
          for( const auto& childSep : sep.children )
              rRowsLocalPack( *childSep );
          for( Int i : sep.inds )
              rRows[offs[A.RowOwner(i)]++] = i;
      };
    function<void(const DistSeparator&,const DistNodeInfo&)> rRowsPack =
      [&]( const DistSeparator& sep, const DistNodeInfo& node )
      {
          if( sep.child == nullptr )
          {
              rRowsLocalPack( *sep.duplicate );
              return;
          }
          rRowsPack( *sep.child, *node.child );

          const Grid& grid = node.Grid();
          const Int rowShift = grid.Col();
          const Int rowStride = grid.Width();
          const Int numInds = sep.inds.size();
          for( Int t=rowShift; t<numInds; t+=rowStride )
          {
              const Int i = sep.inds[t];
              rRows[offs[A.RowOwner(i)]++] = i;
          }
      };
    rRowsPack( rootSep, rootInfo );
    if( time && commRank == 0 )
        Output("Row index pack: ",timer.Stop()," secs");

    // Retreive the list of rows that we must send to each process
    if( time && commRank == 0 )
        timer.Start();
    vector<int> sRowSizes( commSize );
    mpi::AllToAll( rRowSizes.data(), 1, sRowSizes.data(), 1, grid.Comm() );
    vector<int> sRowOffs;
    const Int numSendRows = Scan( sRowSizes, sRowOffs );
    vector<Int> sRows( numSendRows );
    mpi::AllToAll
    ( rRows.data(), rRowSizes.data(), rRowOffs.data(),
      sRows.data(), sRowSizes.data(), sRowOffs.data(), grid.Comm() );
    if( time && commRank == 0 )
        Output("AllToAll: ",timer.Stop()," secs");

    // Pack the number of nonzeros per row (and the nonzeros themselves)
    if( time && commRank == 0 )
        timer.Start();
    const Int firstLocalRow = A.FirstLocalRow();
    vector<Int> sRowLengths( numSendRows );
    vector<int> sEntriesSizes(commSize,0);
    for( Int q=0; q<commSize; ++q )
    {
        const Int size = sRowSizes[q];
        const Int off = sRowOffs[q];
        for( Int s=0; s<size; ++s )
        {
            const Int i = sRows[s+off];
            const Int iLoc = i-firstLocalRow;
            const Int jReord = mappedSources[iLoc];
            const Int rowOff = A.RowOffset( iLoc );
            const Int numConnections = A.NumConnections( iLoc );
            sRowLengths[s+off] = 0;
            for( Int e=0; e<numConnections; ++e )
            {
                const Int iReord = mappedTargets[colOffs[rowOff+e]];
                if( iReord >= jReord )
                {
                    ++sEntriesSizes[q];
                    ++sRowLengths[s+off];
                }
            }
        }
    }
    vector<int> sEntriesOffs;
    const int numSendEntries = Scan( sEntriesSizes, sEntriesOffs );
    vector<Field> sEntries( numSendEntries );
    vector<Int> sTargets( numSendEntries );
    for( Int q=0; q<commSize; ++q )
    {
        Int index = sEntriesOffs[q];
        const Int size = sRowSizes[q];
        const Int off = sRowOffs[q];
        for( Int s=0; s<size; ++s )
        {
            const Int i = sRows[s+off];
            const Int iLoc = i-firstLocalRow;
            const Int jReord = mappedSources[iLoc];
            const Int rowOff = A.RowOffset( iLoc );
            const Int numConnections = A.NumConnections( iLoc );
            for( Int e=0; e<numConnections; ++e )
            {
                const Int iReord = mappedTargets[colOffs[rowOff+e]];
                if( iReord >= jReord )
                {
                    const Field value = A.Value( rowOff+e );
                    sEntries[index] = (conjugate ? Conj(value) : value);
                    sTargets[index] = iReord;
                    ++index;
                }
            }
        }
        EL_DEBUG_ONLY(
          if( index != sEntriesOffs[q]+sEntriesSizes[q] )
              LogicError("index was not the correct value");
        )
    }
    if( time && commRank == 0 )
        Output("Payload pack: ",timer.Stop()," secs");

    // Send back the number of nonzeros per row and the nonzeros themselves
    if( time && commRank == 0 )
        timer.Start();
    vector<Int> rRowLengths( numRecvRows );
    mpi::AllToAll
    ( sRowLengths.data(), sRowSizes.data(), sRowOffs.data(),
      rRowLengths.data(), rRowSizes.data(), rRowOffs.data(), grid.Comm() );
    vector<int> rEntriesSizes(commSize,0);
    for( Int q=0; q<commSize; ++q )
    {
        const Int size = rRowSizes[q];
        const Int off = rRowOffs[q];
        for( Int s=0; s<size; ++s )
            rEntriesSizes[q] += rRowLengths[off+s];
    }
    vector<int> rEntriesOffs;
    const int numRecvEntries = Scan( rEntriesSizes, rEntriesOffs );
    vector<Field> rEntries( numRecvEntries );
    vector<Int> rTargets( numRecvEntries );
    mpi::AllToAll
    ( sEntries.data(), sEntriesSizes.data(), sEntriesOffs.data(),
      rEntries.data(), rEntriesSizes.data(), rEntriesOffs.data(), grid.Comm() );
    mpi::AllToAll
    ( sTargets.data(), sEntriesSizes.data(), sEntriesOffs.data(),
      rTargets.data(), rEntriesSizes.data(), rEntriesOffs.data(), grid.Comm() );
    if( time && commRank == 0 )
        Output("AllToAll time: ",timer.Stop()," secs");

    // Unpack the received entries
    if( time && commRank == 0 )
        timer.Start();
    // TODO(poulson): Modify constructor of [Dist]Front to default to SYMM_2D?
    type = SYMM_2D;
    isHermitian = conjugate;
    UnpackEntries
    ( rootSep, rootInfo, *this,
      A, rRowLengths, rEntries, rTargets, rRowOffs, rEntriesOffs );
    if( time && commRank == 0 )
        Output("Unpack: ",timer.Stop()," secs");
}

template<typename Field>
void DistFront<Field>::PullUpdate
( const DistSparseMatrix<Field>& A,
  const DistMap& reordering,
  const DistSeparator& rootSep,
  const DistNodeInfo& rootInfo )
{
    EL_DEBUG_CSE
    vector<Int> mappedSources, mappedTargets, colOffs;
    PullUpdate
    ( A, reordering, rootSep, rootInfo, mappedSources, mappedTargets, colOffs );
}

// TODO(poulson): Begin removing lambdas in a similar manner as for Pull
template<typename Field>
void DistFront<Field>::PullUpdate
( const DistSparseMatrix<Field>& A,
  const DistMap& reordering,
  const DistSeparator& rootSep,
  const DistNodeInfo& rootInfo,
        vector<Int>& mappedSources,
        vector<Int>& mappedTargets,
        vector<Int>& colOffs )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( A.LocalHeight() != reordering.NumLocalSources() )
          LogicError("Local mapping was not the right size");
    )
    const bool time = false;
    const Grid& grid = A.Grid();
    const int commSize = grid.Size();
    const int commRank = grid.Rank();
    Timer timer;

    A.MappedSources( reordering, mappedSources );
    A.MappedTargets( reordering, mappedTargets, colOffs );

    // Set up the indices for the rows we need from each process
    if( time && commRank == 0 )
        timer.Start();
    vector<int> rRowSizes( commSize, 0 );
    function<void(const Separator&)> rRowLocalAccumulate =
      [&]( const Separator& sep )
      {
          for( const auto& child : sep.children )
             rRowLocalAccumulate( *child );
          for( const Int& i : sep.inds )
              ++rRowSizes[ A.RowOwner(i) ];
      };
    function<void(const DistSeparator&,const DistNodeInfo&)>
      rRowAccumulate =
      [&]( const DistSeparator& sep, const DistNodeInfo& node )
      {
          if( sep.child == nullptr )
          {
              rRowLocalAccumulate( *sep.duplicate );
              return;
          }
          rRowAccumulate( *sep.child, *node.child );

          const Grid& grid = node.Grid();
          const Int rowShift = grid.Col();
          const Int rowStride = grid.Width();
          const Int numInds = sep.inds.size();
          for( Int t=rowShift; t<numInds; t+=rowStride )
              ++rRowSizes[ A.RowOwner(sep.inds[t]) ];
      };
    rRowAccumulate( rootSep, rootInfo );
    vector<int> rRowOffs;
    const Int numRecvRows = Scan( rRowSizes, rRowOffs );
    if( time && commRank == 0 )
        Output("Row index setup: ",timer.Stop()," secs");

    if( time && commRank == 0 )
        timer.Start();
    vector<Int> rRows( numRecvRows );
    auto offs = rRowOffs;
    function<void(const Separator&)> rRowsLocalPack =
      [&]( const Separator& sep )
      {
          for( const auto& childSep : sep.children )
              rRowsLocalPack( *childSep );
          for( Int i : sep.inds )
              rRows[offs[A.RowOwner(i)]++] = i;
      };
    function<void(const DistSeparator&,const DistNodeInfo&)> rRowsPack =
      [&]( const DistSeparator& sep, const DistNodeInfo& node )
      {
          if( sep.child == nullptr )
          {
              rRowsLocalPack( *sep.duplicate );
              return;
          }
          rRowsPack( *sep.child, *node.child );

          const Grid& grid = node.Grid();
          const Int rowShift = grid.Col();
          const Int rowStride = grid.Width();
          const Int numInds = sep.inds.size();
          for( Int t=rowShift; t<numInds; t+=rowStride )
          {
              const Int i = sep.inds[t];
              rRows[offs[A.RowOwner(i)]++] = i;
          }
      };
    rRowsPack( rootSep, rootInfo );
    if( time && commRank == 0 )
        Output("Row index pack: ",timer.Stop()," secs");

    // Retreive the list of rows that we must send to each process
    if( time && commRank == 0 )
        timer.Start();
    vector<int> sRowSizes( commSize );
    mpi::AllToAll( rRowSizes.data(), 1, sRowSizes.data(), 1, grid.Comm() );
    vector<int> sRowOffs;
    const Int numSendRows = Scan( sRowSizes, sRowOffs );
    vector<Int> sRows( numSendRows );
    mpi::AllToAll
    ( rRows.data(), rRowSizes.data(), rRowOffs.data(),
      sRows.data(), sRowSizes.data(), sRowOffs.data(), grid.Comm() );
    if( time && commRank == 0 )
        Output("AllToAll: ",timer.Stop()," secs");

    // Pack the number of nonzeros per row (and the nonzeros themselves)
    if( time && commRank == 0 )
        timer.Start();
    const Int firstLocalRow = A.FirstLocalRow();
    vector<Int> sRowLengths( numSendRows );
    vector<int> sEntriesSizes(commSize,0);
    for( Int q=0; q<commSize; ++q )
    {
        const Int size = sRowSizes[q];
        const Int off = sRowOffs[q];
        for( Int s=0; s<size; ++s )
        {
            const Int i = sRows[s+off];
            const Int iLoc = i-firstLocalRow;
            const Int jReord = mappedSources[iLoc];
            const Int rowOff = A.RowOffset( iLoc );
            const Int numConnections = A.NumConnections( iLoc );
            sRowLengths[s+off] = 0;
            for( Int e=0; e<numConnections; ++e )
            {
                const Int iReord = mappedTargets[colOffs[rowOff+e]];
                if( iReord >= jReord )
                {
                    ++sEntriesSizes[q];
                    ++sRowLengths[s+off];
                }
            }
        }
    }
    vector<int> sEntriesOffs;
    const int numSendEntries = Scan( sEntriesSizes, sEntriesOffs );
    vector<Field> sEntries( numSendEntries );
    vector<Int> sTargets( numSendEntries );
    for( Int q=0; q<commSize; ++q )
    {
        Int index = sEntriesOffs[q];
        const Int size = sRowSizes[q];
        const Int off = sRowOffs[q];
        for( Int s=0; s<size; ++s )
        {
            const Int i = sRows[s+off];
            const Int iLoc = i-firstLocalRow;
            const Int jReord = mappedSources[iLoc];
            const Int rowOff = A.RowOffset( iLoc );
            const Int numConnections = A.NumConnections( iLoc );
            for( Int e=0; e<numConnections; ++e )
            {
                const Int iReord = mappedTargets[colOffs[rowOff+e]];
                if( iReord >= jReord )
                {
                    const Field value = A.Value( rowOff+e );
                    sEntries[index] = (isHermitian ? Conj(value) : value);
                    sTargets[index] = iReord;
                    ++index;
                }
            }
        }
        EL_DEBUG_ONLY(
          if( index != sEntriesOffs[q]+sEntriesSizes[q] )
              LogicError("index was not the correct value");
        )
    }
    if( time && commRank == 0 )
        Output("Payload pack: ",timer.Stop()," secs");

    // Send back the number of nonzeros per row and the nonzeros themselves
    if( time && commRank == 0 )
        timer.Start();
    vector<Int> rRowLengths( numRecvRows );
    mpi::AllToAll
    ( sRowLengths.data(), sRowSizes.data(), sRowOffs.data(),
      rRowLengths.data(), rRowSizes.data(), rRowOffs.data(), grid.Comm() );
    vector<int> rEntriesSizes(commSize,0);
    for( Int q=0; q<commSize; ++q )
    {
        const Int size = rRowSizes[q];
        const Int off = rRowOffs[q];
        for( Int s=0; s<size; ++s )
            rEntriesSizes[q] += rRowLengths[off+s];
    }
    vector<int> rEntriesOffs;
    const int numRecvEntries = Scan( rEntriesSizes, rEntriesOffs );
    vector<Field> rEntries( numRecvEntries );
    vector<Int> rTargets( numRecvEntries );
    mpi::AllToAll
    ( sEntries.data(), sEntriesSizes.data(), sEntriesOffs.data(),
      rEntries.data(), rEntriesSizes.data(), rEntriesOffs.data(), grid.Comm() );
    mpi::AllToAll
    ( sTargets.data(), sEntriesSizes.data(), sEntriesOffs.data(),
      rTargets.data(), rEntriesSizes.data(), rEntriesOffs.data(), grid.Comm() );
    if( time && commRank == 0 )
        Output("AllToAll time: ",timer.Stop()," secs");

    // Unpack the received updates
    if( time && commRank == 0 )
        timer.Start();
    offs = rRowOffs;
    auto entryOffs = rEntriesOffs;
    function<void(const Separator&,const NodeInfo&,Front<Field>&)>
      unpackEntriesLocal =
      [&]( const Separator& sep, const NodeInfo& node, Front<Field>& front )
      {
        const Int numChildren = sep.children.size();
        for( Int c=0; c<numChildren; ++c )
            unpackEntriesLocal
            ( *sep.children[c], *node.children[c], *front.children[c] );

        const Int size = node.size;
        const Int off = node.off;

        if( front.sparseLeaf )
        {
            LogicError("Sparse leaves not supported in DistFront::PullUpdated");
        }
        else
        {
            for( Int t=0; t<size; ++t )
            {
                const Int i = sep.inds[t];
                const Int q = A.RowOwner(i);

                int& entryOff = entryOffs[q];
                const Int numEntries = rRowLengths[offs[q]++];

                for( Int k=0; k<numEntries; ++k )
                {
                    const Field value = rEntries[entryOff];
                    const Int target = rTargets[entryOff];
                    ++entryOff;

                    EL_DEBUG_ONLY(
                      if( target < off+t )
                          LogicError("Received entry from upper triangle");
                    )
                    if( target < off+size )
                    {
                        front.LDense(target-off,t) += value;
                    }
                    else
                    {
                        // TODO: Avoid this binary search?
                        Int origOff = Find( node.origLowerStruct, target );
                        const Int row = node.origLowerRelInds[origOff];
                        front.LDense(row,t) += value;
                    }
                }
            }
        }
      };
    function<void(const DistSeparator&,
                  const DistNodeInfo&,
                        DistFront<Field>&)> unpackEntries =
      [&]( const DistSeparator& sep, const DistNodeInfo& node,
                 DistFront<Field>& front )
      {
        if( sep.child == nullptr )
        {
            unpackEntriesLocal
            ( *sep.duplicate, *node.duplicate, *front.duplicate );
            return;
        }
        unpackEntries( *sep.child, *node.child, *front.child );

        const Int size = node.size;
        const Int off = node.off;
        const Int localWidth = front.L2D.LocalWidth();
        for( Int tLoc=0; tLoc<localWidth; ++tLoc )
        {
            const Int t = front.L2D.GlobalCol(tLoc);
            const Int i = sep.inds[t];
            const Int q = A.RowOwner(i);

            int& entryOff = entryOffs[q];
            const Int numEntries = rRowLengths[offs[q]++];

            for( Int k=0; k<numEntries; ++k )
            {
                const Field value = rEntries[entryOff];
                const Int target = rTargets[entryOff];
                ++entryOff;

                EL_DEBUG_ONLY(
                  if( target < off+t )
                      LogicError("Received entry from upper triangle");
                )
                if( target < off+size )
                {
                    front.L2D.Update( target-off, t, value );
                }
                else
                {
                    // TODO: Avoid this binary search?
                    const Int origOff = Find( node.origLowerStruct, target );
                    const Int row = node.origLowerRelInds[origOff];
                    front.L2D.Update( row, t, value );
                }
            }
        }
      };
    unpackEntries( rootSep, rootInfo, *this );
    if( time && commRank == 0 )
        Output("Unpack: ",timer.Stop()," secs");
    EL_DEBUG_ONLY(
      for( Int q=0; q<commSize; ++q )
          if( entryOffs[q] != rEntriesOffs[q]+rEntriesSizes[q] )
              LogicError("entryOffs were incorrect");
    )
}


template<typename Field>
void DistFront<Field>::Push
( DistSparseMatrix<Field>& A,
  const DistMap& reordering,
  const DistSeparator& rootSep,
  const DistNodeInfo& rootInfo ) const
{
    EL_DEBUG_CSE
    LogicError("This routine needs to be written");
}

template<typename Field>
void DistFront<Field>::Unpack
( DistSparseMatrix<Field>& A,
  const DistSeparator& rootSep,
  const DistNodeInfo& rootInfo ) const
{
    EL_DEBUG_CSE
    A.SetGrid( rootInfo.Grid() );
    const Int n = rootInfo.off + rootInfo.size;
    Zeros( A, n, n );

    const Int numLocalEntries = NumLocalEntries();
    A.Reserve( numLocalEntries, numLocalEntries );

    // Queue the updates
    // =================
    function<void(const Separator&,
                  const NodeInfo&,
                  const Front<Field>&)> localPack =
      [&]( const Separator& sep, const NodeInfo& node,
           const Front<Field>& front )
      {
        const Int numChildren = sep.children.size();
        for( Int c=0; c<numChildren; ++c )
            localPack
            ( *sep.children[c], *node.children[c], *front.children[c] );

        const Int structSize = node.lowerStruct.size();
        if( front.sparseLeaf )
        {
            // Queue the diagonal block
            const Int numEntries = front.LSparse.NumEntries();
            if( numEntries == 0 )
            {
                // Since the matrix has not yet been factored, use the original
                // sparse matrix in front.workSparse
                const Int numWorkEntries = front.workSparse.NumEntries();
                for( Int e=0; e<numWorkEntries; ++e )
                {
                    const Field value = front.workSparse.Value(e);
                    if( value != Field(0) )
                        A.QueueUpdate
                        ( front.workSparse.Row(e)+node.off,
                          front.workSparse.Col(e)+node.off,
                          value );
                }
            }
            else
            {
                // The matrix has already been factored
                for( Int e=0; e<numEntries; ++e )
                {
                    const Field value = front.LSparse.Value(e);
                    if( value != Field(0) )
                        A.QueueUpdate
                        ( front.LSparse.Row(e)+node.off,
                          front.LSparse.Col(e)+node.off,
                          value );
                }
            }

            // Queue the lower conectivity
            for( Int s=0; s<structSize; ++s )
            {
                const Int i = node.lowerStruct[s];
                for( Int t=0; t<node.size; ++t )
                {
                    const Field value = front.LDense(s,t);
                    if( value != Field(0) )
                        A.QueueUpdate( i, t+node.off, value );
                }
            }
        }
        else
        {
            for( Int s=0; s<node.size; ++s )
            {
                const Int i = node.off + s;
                for( Int t=0; t<=s; ++t )
                {
                    const Field value = front.LDense(s,t);
                    if( value != Field(0) )
                        A.QueueUpdate( i, t+node.off, value );
                }
            }

            for( Int s=0; s<structSize; ++s )
            {
                const Int i = node.lowerStruct[s];
                for( Int t=0; t<node.size; ++t )
                {
                    const Field value = front.LDense(node.size+s,t);
                    if( value != Field(0) )
                        A.QueueUpdate( i, t+node.off, value );
                }
            }
        }
      };
    function<void(const DistSeparator&,
                  const DistNodeInfo&,
                  const DistFront<Field>&)> pack =
      [&]( const DistSeparator& sep, const DistNodeInfo& node,
           const DistFront<Field>& front )
      {
        if( sep.duplicate != nullptr )
        {
            localPack( *sep.duplicate, *node.duplicate, *front.duplicate );
            return;
        }
        pack( *sep.child, *node.child, *front.child );

        if( FrontIs1D(front.type) )
        {
            const Int frontHeight = front.L1D.Height();
            auto FTL = front.L1D( IR(0,node.size), IR(0,node.size) );
            auto FBL = front.L1D( IR(node.size,frontHeight), IR(0,node.size) );

            const Int localWidth = FTL.LocalWidth();
            const Int topLocalHeight = FTL.LocalHeight();
            const Int botLocalHeight = FBL.LocalHeight();

            for( Int sLoc=0; sLoc<topLocalHeight; ++sLoc )
            {
                const Int s = FTL.GlobalRow(sLoc);
                const Int i = node.off + s;
                for( Int tLoc=0; tLoc<localWidth; ++tLoc )
                {
                    const Int t = FTL.GlobalCol(tLoc);
                    if( t <= s )
                    {
                        const Field value = FTL.GetLocal(sLoc,tLoc);
                        if( value != Field(0) )
                            A.QueueUpdate( i, t+node.off, value );
                    }
                }
            }

            for( Int sLoc=0; sLoc<botLocalHeight; ++sLoc )
            {
                const Int s = FBL.GlobalRow(sLoc);
                const Int i = node.lowerStruct[s];
                for( Int tLoc=0; tLoc<localWidth; ++tLoc )
                {
                    const Int t = FBL.GlobalCol(tLoc);
                    const Field value = FBL.GetLocal(sLoc,tLoc);
                    if( value != Field(0) )
                        A.QueueUpdate( i, t+node.off, value );
                }
            }
        }
        else
        {
            const Int frontHeight = front.L2D.Height();
            auto FTL = front.L2D( IR(0,node.size), IR(0,node.size) );
            auto FBL = front.L2D( IR(node.size,frontHeight), IR(0,node.size) );

            const Int localWidth = FTL.LocalWidth();
            const Int topLocalHeight = FTL.LocalHeight();
            const Int botLocalHeight = FBL.LocalHeight();

            for( Int sLoc=0; sLoc<topLocalHeight; ++sLoc )
            {
                const Int s = FTL.GlobalRow(sLoc);
                const Int i = node.off + s;
                for( Int tLoc=0; tLoc<localWidth; ++tLoc )
                {
                    const Int t = FTL.GlobalCol(tLoc);
                    if( t <= s )
                    {
                        const Field value = FTL.GetLocal(sLoc,tLoc);
                        if( value != Field(0) )
                            A.QueueUpdate( i, t+node.off, value );
                    }
                }
            }

            for( Int sLoc=0; sLoc<botLocalHeight; ++sLoc )
            {
                const Int s = FBL.GlobalRow(sLoc);
                const Int i = node.lowerStruct[s];
                for( Int tLoc=0; tLoc<localWidth; ++tLoc )
                {
                    const Int t = FBL.GlobalCol(tLoc);
                    const Field value = FBL.GetLocal(sLoc,tLoc);
                    if( value != Field(0) )
                        A.QueueUpdate( i, t+node.off, value );
                }
            }
        }
      };
    pack( rootSep, rootInfo, *this );

    A.ProcessQueues();
}

template<typename Field>
const DistFront<Field>&
DistFront<Field>::operator=( const DistFront<Field>& front )
{
    EL_DEBUG_CSE
    isHermitian = front.isHermitian;
    type = front.type;
    if( front.child == nullptr )
    {
        child = nullptr;
        duplicate.reset( new Front<Field>(this) );
        *duplicate = *front.duplicate;
        const Grid& grid = front.L2D.Grid();
        L2D.Attach( grid, front.duplicate->LDense );
        diag.Attach( grid, front.duplicate->diag );
        subdiag.Attach( grid, front.duplicate->subdiag );
        p.SetGrid( grid );
        p = front.duplicate->p;
    }
    else
    {
        duplicate = nullptr;
        child.reset( new DistFront<Field>(this) );
        *child = *front.child;
        L1D = front.L1D;
        L2D = front.L2D;
        diag = front.diag;
        subdiag = front.subdiag;
        p = front.p;
        work = front.work;
    }
    return *this;
}

template<typename Field>
Int DistFront<Field>::NumLocalEntries() const
{
    EL_DEBUG_CSE
    Int numEntries = 0;
    function<void(const DistFront<Field>&)> count =
      [&]( const DistFront<Field>& front )
      {
        if( front.duplicate != nullptr )
        {
            numEntries += front.duplicate->NumEntries();
            return;
        }
        count( *front.child );

        // Add in L
        numEntries += front.L1D.LocalHeight() * front.L1D.LocalWidth();
        numEntries += front.L2D.LocalHeight() * front.L2D.LocalWidth();

        // Add in the workspace
        numEntries += front.work.LocalHeight() * front.work.LocalWidth();
      };
    count( *this );
    return numEntries;
}

template<typename Field>
Int DistFront<Field>::NumTopLeftLocalEntries() const
{
    EL_DEBUG_CSE
    Int numEntries = 0;
    function<void(const DistFront<Field>&)> count =
      [&]( const DistFront<Field>& front )
      {
        if( front.duplicate != nullptr )
        {
            numEntries += front.duplicate->NumTopLeftEntries();
            return;
        }
        count( *front.child );

        if( FrontIs1D(front.type) )
        {
            const Int n = front.L1D.Width();
            auto FTL = front.L1D( IR(0,n), IR(0,n) );
            numEntries += FTL.LocalHeight() * FTL.LocalWidth();
        }
        else
        {
            const Int n = front.L2D.Width();
            auto FTL = front.L2D( IR(0,n), IR(0,n) );
            numEntries += FTL.LocalHeight() * FTL.LocalWidth();
        }
      };
    count( *this );
    return numEntries;
}

template<typename Field>
Int DistFront<Field>::NumBottomLeftLocalEntries() const
{
    EL_DEBUG_CSE
    Int numEntries = 0;
    function<void(const DistFront<Field>&)> count =
      [&]( const DistFront<Field>& front )
      {
        if( front.duplicate != nullptr )
        {
            numEntries += front.duplicate->NumBottomLeftEntries();
            return;
        }
        count( *front.child );

        if( FrontIs1D(front.type) )
        {
            const Int m = front.L1D.Height();
            const Int n = front.L1D.Width();
            auto FBL = front.L2D( IR(n,m), IR(0,n) );
            numEntries += FBL.LocalHeight() * FBL.LocalWidth();
        }
        else
        {
            const Int m = front.L2D.Height();
            const Int n = front.L2D.Width();
            auto FBL = front.L2D( IR(n,m), IR(0,n) );
            numEntries += FBL.LocalHeight() * FBL.LocalWidth();
        }
      };
    count( *this );
    return numEntries;
}

template<typename Field>
double DistFront<Field>::LocalFactorGFlops( bool selInv ) const
{
    EL_DEBUG_CSE
    double gflops = 0.;
    function<void(const DistFront<Field>&)> count =
      [&]( const DistFront<Field>& front )
      {
        if( front.duplicate != nullptr )
        {
            gflops += front.duplicate->FactorGFlops();
            return;
        }
        count( *front.child );

        double m, n, p;
        if( FrontIs1D(front.type) )
        {
            m = front.L1D.Height();
            n = front.L1D.Width();
            p = front.L1D.DistSize();
        }
        else
        {
            m = front.L2D.Height();
            n = front.L2D.Width();
            p = front.L2D.DistSize();
        }
        double realFrontFlops =
          ( selInv ? (2*n*n*n/3) + (m-n)*n + (m-n)*(m-n)*n
                   : (1*n*n*n/3) + (m-n)*n + (m-n)*(m-n)*n ) / p;
        gflops += (IsComplex<Field>::value ? 4*realFrontFlops
                                           : realFrontFlops)/1.e9;
      };
    count( *this );
    return gflops;
}

template<typename Field>
double DistFront<Field>::LocalSolveGFlops( Int numRHS ) const
{
    EL_DEBUG_CSE
    double gflops = 0.;
    function<void(const DistFront<Field>&)> count =
      [&]( const DistFront<Field>& front )
      {
        if( front.duplicate != nullptr )
        {
            gflops += front.duplicate->SolveGFlops( numRHS );
            return;
        }
        count( *front.child );

        double m, n, p;
        if( FrontIs1D(front.type) )
        {
            m = front.L1D.Height();
            n = front.L1D.Width();
            p = front.L1D.DistSize();
        }
        else
        {
            m = front.L2D.Height();
            n = front.L2D.Width();
            p = front.L2D.DistSize();
        }
        double realFrontFlops = (m*n*numRHS) / p;
        gflops += (IsComplex<Field>::value ? 4*realFrontFlops
                                           : realFrontFlops)/1.e9;
      };
    count( *this );
    return gflops;
}

template<typename Field>
void DistFront<Field>::ComputeRecvInds( const DistNodeInfo& info ) const
{
    EL_DEBUG_CSE

    vector<int> gridHeights, gridWidths;
    info.GetChildGridDims( gridHeights, gridWidths );

    const Grid& grid = info.Grid();
    const Grid& childGrid = info.child->Grid();

    const int teamSize = grid.Size();
    const int teamRank = grid.Rank();
    const bool onLeft = info.child->onLeft;
    const int childTeamSize = childGrid.Size();
    const int childTeamRank = childGrid.Rank();
    const bool inFirstTeam = ( childTeamRank == teamRank );
    const bool leftIsFirst = ( onLeft==inFirstTeam );
    vector<int> teamSizes(2), teamOffs(2);
    teamSizes[0] = ( onLeft ? childTeamSize : teamSize-childTeamSize );
    teamSizes[1] = teamSize - teamSizes[0];
    teamOffs[0] = ( leftIsFirst ? 0            : teamSizes[1] );
    teamOffs[1] = ( leftIsFirst ? teamSizes[0] : 0            );
    EL_DEBUG_ONLY(
      if( teamSizes[0] != gridHeights[0]*gridWidths[0] )
          RuntimeError("Computed left grid incorrectly");
      if( teamSizes[1] != gridHeights[1]*gridWidths[1] )
          RuntimeError("Computed right grid incorrectly");
    )

    commMeta.childRecvInds.resize( teamSize );
    for( int q=0; q<teamSize; ++q )
        commMeta.childRecvInds[q].clear();
    for( Int c=0; c<2; ++c )
    {
        // Compute the recv indices of the child from each process
        const Int numInds = info.childRelInds[c].size();
        vector<Int> rowInds, colInds;
        for( Int iChild=0; iChild<numInds; ++iChild )
        {
            if( L2D.IsLocalRow( info.childRelInds[c][iChild] ) )
                rowInds.push_back( iChild );
            if( L2D.IsLocalCol( info.childRelInds[c][iChild] ) )
                colInds.push_back( iChild );
        }

        vector<Int>::const_iterator it;
        const Int numColInds = colInds.size();
        const Int numRowInds = rowInds.size();
        for( Int jPre=0; jPre<numColInds; ++jPre )
        {
            const Int jChild = colInds[jPre];
            const Int j = info.childRelInds[c][jChild];
            const Int jLoc = L2D.LocalCol(j);

            const int childCol = (jChild+info.childSizes[c]) % gridWidths[c];

            // Find the first iPre that maps to the lower triangle
            it = std::lower_bound( rowInds.begin(), rowInds.end(), jChild );
            const Int iPreStart = Int(it-rowInds.begin());
            for( Int iPre=iPreStart; iPre<numRowInds; ++iPre )
            {
                const Int iChild = rowInds[iPre];
                const Int i = info.childRelInds[c][iChild];
                EL_DEBUG_ONLY(
                  if( iChild < jChild )
                      LogicError("Invalid iChild");
                )
                const Int iLoc = L2D.LocalRow(i);

                const int childRow = (iChild+info.childSizes[c])%gridHeights[c];
                const int childRank = childRow + childCol*gridHeights[c];

                const int q = teamOffs[c]+ childRank;
                commMeta.childRecvInds[q].push_back(iLoc);
                commMeta.childRecvInds[q].push_back(jLoc);
            }
        }
    }
}

template<typename Field>
void DistFront<Field>::ComputeCommMeta
( const DistNodeInfo& info, bool computeRecvInds ) const
{
    EL_DEBUG_CSE
    commMeta.Empty();
    if( child == nullptr )
        return;

    const int teamSize = L2D.DistSize();
    commMeta.numChildSendInds.resize( teamSize );
    El::MemZero( commMeta.numChildSendInds.data(), teamSize );

    auto& childFront = *child;
    auto& childInfo = *info.child;
    const auto& childRelInds =
      ( childInfo.onLeft ? info.childRelInds[0] : info.childRelInds[1] );

    const auto& FBR = childFront.work;
    const Int localHeight = FBR.LocalHeight();
    const Int localWidth = FBR.LocalWidth();
    for( Int jChildLoc=0; jChildLoc<localWidth; ++jChildLoc )
    {
        const Int jChild = FBR.GlobalCol(jChildLoc);
        const int j = childRelInds[jChild];
        for( Int iChildLoc=0; iChildLoc<localHeight; ++iChildLoc )
        {
            const Int iChild = FBR.GlobalRow(iChildLoc);
            if( iChild >= jChild )
            {
                const Int i = childRelInds[iChild];
                const int q = L2D.Owner( i, j );
                ++commMeta.numChildSendInds[q];
            }
        }
    }

    // This is optional since it requires a nontrivial amount of storage.
    if( computeRecvInds )
        ComputeRecvInds( info );
}

#define PROTO(Field) template struct DistFront<Field>;
#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace ldl
} // namespace El
