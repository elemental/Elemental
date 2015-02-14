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
   
   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename F>
DistSymmFront<F>::DistSymmFront( DistSymmFront<F>* parentNode )
: parent(parentNode), child(nullptr), duplicate(nullptr)
{ }

template<typename F>
DistSymmFront<F>::DistSymmFront
( const DistSparseMatrix<F>& A, 
  const DistMap& reordering,
  const DistSeparator& sep, 
  const DistSymmNodeInfo& info,
  bool conjugate )
: parent(nullptr), child(nullptr), duplicate(nullptr)
{
    DEBUG_ONLY(CallStackEntry cse("DistSymmFront::DistSymmFront"))
    Pull( A, reordering, sep, info, conjugate );
}

template<typename F>
void DistSymmFront<F>::Pull
( const DistSparseMatrix<F>& A, 
  const DistMap& reordering,
  const DistSeparator& rootSep, 
  const DistSymmNodeInfo& rootInfo,
  bool conjugate )
{
    DEBUG_ONLY(
      CallStackEntry cse("DistSymmFront::Pull");
      if( A.LocalHeight() != reordering.NumLocalSources() )
          LogicError("Local mapping was not the right size");
    )
   
    mpi::Comm comm = A.Comm();
    const DistGraph& graph = A.LockedDistGraph();
    const Int commSize = mpi::Size( comm );

    // Get the reordered indices of the targets of our portion of the 
    // distributed sparse matrix
    set<Int> targetSet( graph.targets_.begin(), graph.targets_.end() );
    vector<Int> targets;
    CopySTL( targetSet, targets );
    auto mappedTargets = targets;
    reordering.Translate( mappedTargets );

    // Set up the indices for the rows we need from each process
    vector<int> rRowSizes( commSize, 0 );
    function<void(const Separator&)> rRowLocalAccumulate = 
      [&]( const Separator& sep )
      {
          for( const Separator* child : sep.children )
             rRowLocalAccumulate( *child );
          for( const Int& i : sep.inds )
              ++rRowSizes[ A.RowOwner(i) ];
      };
    function<void(const DistSeparator&,const DistSymmNodeInfo&)> 
      rRowAccumulate =
      [&]( const DistSeparator& sep, const DistSymmNodeInfo& node )
      {
          if( sep.child == nullptr )
          {
              rRowLocalAccumulate( *sep.duplicate ); 
              return;
          }
          rRowAccumulate( *sep.child, *node.child );
          
          const Grid& grid = *node.grid;
          const Int rowShift = grid.Col();
          const Int rowStride = grid.Width();
          const Int numInds = sep.inds.size();
          for( Int t=rowShift; t<numInds; t+=rowStride )
              ++rRowSizes[ A.RowOwner(sep.inds[t]) ];
      };
    rRowAccumulate( rootSep, rootInfo );
    vector<int> rRowOffs;
    const Int numRecvRows = Scan( rRowSizes, rRowOffs );

    vector<Int> rRows( numRecvRows );
    auto offs = rRowOffs;
    function<void(const Separator&)> rRowsLocalPack = 
      [&]( const Separator& sep )
      {
          for( const Separator* childSep : sep.children )
              rRowsLocalPack( *childSep );
          for( Int i : sep.inds )
              rRows[offs[A.RowOwner(i)]++] = i;
      };
    function<void(const DistSeparator&,const DistSymmNodeInfo&)> rRowsPack = 
      [&]( const DistSeparator& sep, const DistSymmNodeInfo& node )
      {
          if( sep.child == nullptr )
          {
              rRowsLocalPack( *sep.duplicate );
              return;
          }
          rRowsPack( *sep.child, *node.child );
          
          const Grid& grid = *node.grid;
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

    // Retreive the list of rows that we must send to each process
    vector<int> sRowSizes( commSize );
    mpi::AllToAll( rRowSizes.data(), 1, sRowSizes.data(), 1, comm );
    vector<int> sRowOffs;
    const Int numSendRows = Scan( sRowSizes, sRowOffs );
    vector<Int> sRows( numSendRows );
    mpi::AllToAll
    ( rRows.data(), rRowSizes.data(), rRowOffs.data(),
      sRows.data(), sRowSizes.data(), sRowOffs.data(), comm );

    // Pack the number of nonzeros per row (and the nonzeros themselves)
    // TODO: Avoid sending upper-triangular data
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
            const Int numConnections = A.NumConnections( i-firstLocalRow );
            sEntriesSizes[q] += numConnections;
            sRowLengths[s+off] = numConnections;
        }
    }
    vector<int> sEntriesOffs;
    const int numSendEntries = Scan( sEntriesSizes, sEntriesOffs );
    vector<F> sEntries( numSendEntries );
    vector<Int> sTargets( numSendEntries );
    for( Int q=0; q<commSize; ++q )
    {
        Int index = sEntriesOffs[q];
        const Int size = sRowSizes[q];
        const Int off = sRowOffs[q];
        for( Int s=0; s<size; ++s )
        {
            const Int i = sRows[s+off];
            const Int numConnections = sRowLengths[s+off];
            const Int localEntryOff = A.EntryOffset( i-firstLocalRow );
            for( Int t=0; t<numConnections; ++t )
            {
                const F value = A.Value( localEntryOff+t );
                const Int col = A.Col( localEntryOff+t );
                const Int targetOff = Find( targets, col );
                const Int mappedTarget = mappedTargets[targetOff];
                sEntries[index] = (conjugate ? Conj(value) : value);
                sTargets[index] = mappedTarget;
                ++index;
            }
        }
        DEBUG_ONLY(
          if( index != sEntriesOffs[q]+sEntriesSizes[q] )
              LogicError("index was not the correct value");
        )
    }

    // Send back the number of nonzeros per row and the nonzeros themselves
    vector<Int> rRowLengths( numRecvRows );
    mpi::AllToAll
    ( sRowLengths.data(), sRowSizes.data(), sRowOffs.data(),
      rRowLengths.data(), rRowSizes.data(), rRowOffs.data(), comm );
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
    vector<F> rEntries( numRecvEntries );
    vector<Int> rTargets( numRecvEntries );
    mpi::AllToAll
    ( sEntries.data(), sEntriesSizes.data(), sEntriesOffs.data(),
      rEntries.data(), rEntriesSizes.data(), rEntriesOffs.data(), comm );
    mpi::AllToAll
    ( sTargets.data(), sEntriesSizes.data(), sEntriesOffs.data(),
      rTargets.data(), rEntriesSizes.data(), rEntriesOffs.data(), comm );

    // Unpack the received entries
    offs = rRowOffs;
    auto entryOffs = rEntriesOffs;
    function<void(const Separator&,const SymmNodeInfo&,SymmFront<F>&)> 
      unpackEntriesLocal = 
      [&]( const Separator& sep, const SymmNodeInfo& node, SymmFront<F>& front )
      {
          front.type = SYMM_2D;
          front.isHermitian = conjugate;
 
          const Int numChildren = sep.children.size();
          front.children.resize( numChildren );
          for( Int c=0; c<numChildren; ++c )
          {
              front.children[c] = new SymmFront<F>(&front);
              unpackEntriesLocal
              ( *sep.children[c], *node.children[c], *front.children[c] );
          }

          const Int size = node.size;
          const Int off = node.off;
          const Int lowerSize = node.lowerStruct.size();
          Zeros( front.L, size+lowerSize, size );

          for( Int t=0; t<size; ++t )
          {
              const Int i = sep.inds[t];
              const Int q = A.RowOwner(i);

              int& entryOff = entryOffs[q];
              const Int numEntries = rRowLengths[offs[q]++];

              for( Int k=0; k<numEntries; ++k )
              {
                  const F value = rEntries[entryOff];
                  const Int target = rTargets[entryOff];
                  ++entryOff;
  
                  if( target < off+t )
                      continue;
                  else if( target < off+size )
                  {
                      front.L.Set( target-off, t, value );
                  }
                  else
                  {
                      const Int origOff = Find( node.origLowerStruct, target );
                      const Int row = node.origLowerRelInds[origOff];
                      front.L.Set( row, t, value );
                  }
              }
          }
      };
    function<void(const DistSeparator&,
                  const DistSymmNodeInfo&,
                        DistSymmFront<F>&)> unpackEntries = 
      [&]( const DistSeparator& sep, const DistSymmNodeInfo& node, 
                 DistSymmFront<F>& front )
      {
          front.type = SYMM_2D;
          front.isHermitian = conjugate;
          const Grid& grid = *node.grid;

          if( sep.child == nullptr )
          {
              front.duplicate = new SymmFront<F>(&front);
              unpackEntriesLocal
              ( *sep.duplicate, *node.duplicate, *front.duplicate );

              front.L2D.Attach( grid, front.duplicate->L );

              return;
          }
          front.child = new DistSymmFront<F>(&front);
          unpackEntries( *sep.child, *node.child, *front.child );

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
                  const F value = rEntries[entryOff];
                  const Int target = rTargets[entryOff];
                  ++entryOff;

                  if( target < off+t )
                      continue;
                  else if( target < off+size )
                  {
                      front.L2D.Set( target-off, t, value );
                  }
                  else 
                  {
                      const Int origOff = Find( node.origLowerStruct, target );
                      const Int row = node.origLowerRelInds[origOff];
                      front.L2D.Set( row, t, value );
                  }
              }
          }
      };
    unpackEntries( rootSep, rootInfo, *this );
    DEBUG_ONLY(
      for( Int q=0; q<commSize; ++q )
          if( entryOffs[q] != rEntriesOffs[q]+rEntriesSizes[q] )
              LogicError("entryOffs were incorrect");
    )
}

template<typename F>
void DistSymmFront<F>::Push
( DistSparseMatrix<F>& A, 
  const DistMap& reordering,
  const DistSeparator& rootSep, 
  const DistSymmNodeInfo& rootInfo ) const
{
    DEBUG_ONLY(CallStackEntry cse("DistSymmFront::Push"))
    LogicError("This routine needs to be written");
}

template<typename F>
void DistSymmFront<F>::Unpack
( DistSparseMatrix<F>& A, 
  const DistSeparator& rootSep, 
  const DistSymmNodeInfo& rootInfo ) const
{
    DEBUG_ONLY(CallStackEntry cse("DistSymmFront::Unpack"))
    mpi::Comm comm = rootInfo.grid->Comm();
    const int commSize = mpi::Size(comm);
    A.SetComm( comm );
    const Int n = rootInfo.off + rootInfo.size;
    Zeros( A, n, n );
    
    // Compute the metadata for sending the entries of the frontal tree
    // TODO: Avoid metadata redundancy within each row of a front
    vector<int> sendSizes(commSize,0);
    function<void(const Separator&,
                  const SymmNodeInfo&,
                  const SymmFront<F>&)> localCount =
      [&]( const Separator& sep, const SymmNodeInfo& node, 
           const SymmFront<F>& front )
      {
        const Int numChildren = sep.children.size();
        for( Int c=0; c<numChildren; ++c )
            localCount
            ( *sep.children[c], *node.children[c], *front.children[c] );

        for( Int s=0; s<node.size; ++s )
        {
            const int q = A.RowOwner(node.off+s);
            for( Int t=0; t<=s; ++t ) 
                sendSizes[q]++;
        }

        const Int structSize = node.lowerStruct.size();
        for( Int s=0; s<structSize; ++s ) 
        {
            const int q = A.RowOwner(node.lowerStruct[s]);
            for( Int t=0; t<node.size; ++t )
                sendSizes[q]++;
        }
      };
    function<void(const DistSeparator&,
                  const DistSymmNodeInfo&,
                  const DistSymmFront<F>&)> count =  
      [&]( const DistSeparator& sep, const DistSymmNodeInfo& node,
           const DistSymmFront<F>& front )
      {
        if( sep.duplicate != nullptr )
        {
            localCount( *sep.duplicate, *node.duplicate, *front.duplicate );
            return;
        }
        count( *sep.child, *node.child, *front.child );

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
                const int q = A.RowOwner(i);
                for( Int tLoc=0; tLoc<localWidth; ++tLoc )
                    if( FTL.GlobalCol(tLoc) <= s )
                        ++sendSizes[q];
            }

            for( Int sLoc=0; sLoc<botLocalHeight; ++sLoc )
            {
                const Int s = FBL.GlobalRow(sLoc);
                const Int i = node.lowerStruct[s];
                const int q = A.RowOwner(i);
                for( Int tLoc=0; tLoc<localWidth; ++tLoc )
                    ++sendSizes[q];
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
                const int q = A.RowOwner(node.off+s);
                for( Int tLoc=0; tLoc<localWidth; ++tLoc )
                    if( FTL.GlobalCol(tLoc) <= s )
                        ++sendSizes[q];
            }

            for( Int sLoc=0; sLoc<botLocalHeight; ++sLoc )
            {
                const Int s = FBL.GlobalRow(sLoc);
                const Int i = node.lowerStruct[s];
                const int q = A.RowOwner(i);
                for( Int jLoc=0; jLoc<localWidth; ++jLoc )
                    ++sendSizes[q];
            }
        }
      };
    count( rootSep, rootInfo, *this );
    vector<int> recvSizes(commSize);
    mpi::AllToAll( sendSizes.data(), 1, recvSizes.data(), 1, comm );
    vector<int> sendOffs, recvOffs;
    const int totalSend = Scan( sendSizes, sendOffs );
    const int totalRecv = Scan( recvSizes, recvOffs );

    // Pack the data
    vector<Int> iSendBuf(totalSend), jSendBuf(totalSend);
    vector<F> vSendBuf(totalSend);
    auto offs = sendOffs;
    function<void(const Separator&,
                  const SymmNodeInfo&,
                  const SymmFront<F>&)> localPack =
      [&]( const Separator& sep, const SymmNodeInfo& node, 
           const SymmFront<F>& front )
      {
        const Int numChildren = sep.children.size();
        for( Int c=0; c<numChildren; ++c )
            localPack
            ( *sep.children[c], *node.children[c], *front.children[c] );

        for( Int s=0; s<node.size; ++s )
        {
            const Int i = node.off + s;
            const int q = A.RowOwner(i);
            for( Int t=0; t<=s; ++t ) 
            {
                iSendBuf[offs[q]] = i;
                jSendBuf[offs[q]] = node.off+t;
                vSendBuf[offs[q]] = front.L.Get(s,t);
                ++offs[q];
            }
        }

        const Int structSize = node.lowerStruct.size();
        for( Int s=0; s<structSize; ++s ) 
        {
            const Int i = node.lowerStruct[s];
            const int q = A.RowOwner(i);
            for( Int t=0; t<node.size; ++t )
            {
                iSendBuf[offs[q]] = i;
                jSendBuf[offs[q]] = node.off+t;
                vSendBuf[offs[q]] = front.L.Get(node.size+s,t);
                ++offs[q];
            }
        }
      };
    function<void(const DistSeparator&,
                  const DistSymmNodeInfo&,
                  const DistSymmFront<F>&)> pack =  
      [&]( const DistSeparator& sep, const DistSymmNodeInfo& node,
           const DistSymmFront<F>& front )
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
                const int q = A.RowOwner(i);
                for( Int tLoc=0; tLoc<localWidth; ++tLoc )
                {
                    const Int t = FTL.GlobalCol(tLoc);
                    if( t <= s )
                    {
                        iSendBuf[offs[q]] = i;
                        jSendBuf[offs[q]] = node.off+t;
                        vSendBuf[offs[q]] = FTL.GetLocal(sLoc,tLoc);
                        ++offs[q];
                    }
                }
            }

            for( Int sLoc=0; sLoc<botLocalHeight; ++sLoc )
            {
                const Int s = FBL.GlobalRow(sLoc);
                const Int i = node.lowerStruct[s];
                const int q = A.RowOwner(i);
                for( Int tLoc=0; tLoc<localWidth; ++tLoc )
                {
                    const Int t = FBL.GlobalCol(tLoc);
                    iSendBuf[offs[q]] = i;
                    jSendBuf[offs[q]] = node.off+t;
                    vSendBuf[offs[q]] = FBL.GetLocal(sLoc,tLoc);
                    ++offs[q];
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
                const int q = A.RowOwner(i);
                for( Int tLoc=0; tLoc<localWidth; ++tLoc )
                {
                    const Int t = FTL.GlobalCol(tLoc);
                    if( t <= s )
                    {
                        iSendBuf[offs[q]] = i;
                        jSendBuf[offs[q]] = node.off + t;
                        vSendBuf[offs[q]] = FTL.GetLocal(sLoc,tLoc);
                        ++offs[q];
                    }
                }
            }

            for( Int sLoc=0; sLoc<botLocalHeight; ++sLoc )
            {
                const Int s = FBL.GlobalRow(sLoc);
                const Int i = node.lowerStruct[s];
                const int q = A.RowOwner(i);
                for( Int tLoc=0; tLoc<localWidth; ++tLoc )
                {
                    const Int t = FBL.GlobalCol(tLoc);
                    iSendBuf[offs[q]] = i;
                    jSendBuf[offs[q]] = node.off + t;
                    vSendBuf[offs[q]] = FBL.GetLocal(sLoc,tLoc);
                    ++offs[q];
                }
            }
        }
      };
    pack( rootSep, rootInfo, *this );

    // Exchange the data
    vector<Int> iRecvBuf(totalRecv), jRecvBuf(totalRecv);
    vector<F> vRecvBuf(totalRecv);
    mpi::AllToAll
    ( iSendBuf.data(), sendSizes.data(), sendOffs.data(),
      iRecvBuf.data(), recvSizes.data(), recvOffs.data(), comm );
    mpi::AllToAll
    ( jSendBuf.data(), sendSizes.data(), sendOffs.data(),
      jRecvBuf.data(), recvSizes.data(), recvOffs.data(), comm );
    mpi::AllToAll
    ( vSendBuf.data(), sendSizes.data(), sendOffs.data(),
      vRecvBuf.data(), recvSizes.data(), recvOffs.data(), comm );
    
    // Unpack the data
    A.Reserve( totalRecv );
    const Int firstLocalRow = A.FirstLocalRow();
    for( Int e=0; e<totalRecv; ++e )
        A.QueueLocalUpdate
        ( iRecvBuf[e]-firstLocalRow, jRecvBuf[e], vRecvBuf[e] );
    A.MakeConsistent();

    MakeSymmetric( LOWER, A, isHermitian );
}

template<typename F>
double DistSymmFront<F>::NumLocalEntries() const
{
    DEBUG_ONLY(CallStackEntry cse("DistSymmFront::NumLocalEntries"))
    LogicError("This routine needs to be rewritten");
    return 0.;
}

template<typename F>
double DistSymmFront<F>::NumTopLeftLocalEntries() const
{
    DEBUG_ONLY(CallStackEntry cse("DistSymmFront::NumTopLeftLocalEntries"))
    LogicError("This routine needs to be rewritten");
    return 0.;
}

template<typename F>
double DistSymmFront<F>::NumBottomLeftLocalEntries() const
{
    DEBUG_ONLY(CallStackEntry cse("DistSymmFront::NumBottomLeftLocalEntries"))
    LogicError("This routine needs to be rewritten");
    return 0.;
}

template<typename F>
double DistSymmFront<F>::LocalFactorGFlops( bool selInv ) const
{
    DEBUG_ONLY(CallStackEntry cse("DistSymmFront::LocalFactorGFlops"))
    LogicError("This routine needs to be rewritten");
    return 0.;
}

template<typename F>
double DistSymmFront<F>::LocalSolveGFlops( Int numRHS ) const
{
    DEBUG_ONLY(CallStackEntry cse("DistSymmFront::LocalSolveGFlops"))
    LogicError("This routine needs to be rewritten");
    return 0.;
}

#define PROTO(F) template class DistSymmFront<F>;
#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
