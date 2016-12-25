/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   Copyright (c) 2012 Jack Poulson, Lexing Ying, and
   The University of Texas at Austin.
   All rights reserved.

   Copyright (c) 2013 Jack Poulson, Lexing Ying, and Stanford University.
   All rights reserved.

   Copyright (c) 2014 Jack Poulson and The Georgia Institute of Technology.
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {
namespace ldl {

Separator::Separator( Separator* parentNode )
: parent(parentNode)
{ }

Separator::Separator( DistSeparator* dupNode )
: duplicate(dupNode)
{
    EL_DEBUG_CSE
    off = duplicate->off;
    inds = duplicate->inds;
}

Separator::~Separator() { }

void Separator::BuildMap( vector<Int>& map ) const
{
    EL_DEBUG_CSE
    const Int numSources = off + inds.size();
    map.resize( numSources );

    function<void(const Separator&)> buildMap =
      [&]( const Separator& sep )
      {
        for( const auto& child : sep.children )
            buildMap( *child );
        for( size_t t=0; t<sep.inds.size(); ++t )
            map[sep.inds[t]] = sep.off + t;
      };
    buildMap( *this );
}

DistSeparator::DistSeparator( DistSeparator* parentNode )
: parent(parentNode)
{ }

DistSeparator::~DistSeparator() { }

void DistSeparator::BuildMap
( const DistNodeInfo& rootInfo, DistMap& map ) const
{
    EL_DEBUG_CSE
    const Int numSources = off + inds.size();
    const Grid& rootGrid = rootInfo.Grid();
    const int commSize = rootGrid.Size();

    map.SetGrid( rootGrid );
    map.Resize( numSources );

    vector<int> sendSizes( commSize, 0 );
    function<void(const NodeInfo&,const Separator&)> sendSizeLocalAccumulate =
      [&]( const NodeInfo& info, const Separator& sep )
      {
        const int numChildren = info.children.size();
        for( Int index=0; index<numChildren; ++index )
            sendSizeLocalAccumulate
            ( *info.children[index], *sep.children[index] );
        for( Int i : sep.inds )
            ++sendSizes[ map.RowOwner(i) ];
      };
    function<void(const DistNodeInfo&,const DistSeparator&)>
      sendSizeAccumulate =
      [&]( const DistNodeInfo& info, const DistSeparator& sep )
      {
          if( sep.child == nullptr )
          {
              const NodeInfo& infoDup = *info.duplicate;
              const Separator& sepDup = *sep.duplicate;
              const int numChildren = sepDup.children.size();
              for( Int index=0; index<numChildren; ++index )
                  sendSizeLocalAccumulate
                  ( *infoDup.children[index], *sepDup.children[index] );
          }
          else
              sendSizeAccumulate( *info.child, *sep.child );

          const Int numInds = sep.inds.size();
          const Grid& grid = info.Grid();
          const int teamSize = grid.Size();
          const int teamRank = grid.Rank();
          const Int numLocalInds = Length( numInds, teamRank, teamSize );
          for( Int tLocal=0; tLocal<numLocalInds; ++tLocal )
          {
              const Int t = teamRank + tLocal*teamSize;
              ++sendSizes[ map.RowOwner(sep.inds[t]) ];
          }
      };
    sendSizeAccumulate( rootInfo, *this );

    // Use a single-entry AllToAll to coordinate how many indices will be
    // exchanges
    vector<int> recvSizes( commSize );
    mpi::AllToAll( sendSizes.data(), 1, recvSizes.data(), 1, rootGrid.Comm() );

    // Pack the reordered indices
    vector<int> sendOffs;
    const int numSends = Scan( sendSizes, sendOffs );
    vector<Int> sendInds(numSends), sendOrigInds(numSends);
    auto offs = sendOffs;
    function<void(const NodeInfo&,const Separator&)> packRowsLocal =
      [&]( const NodeInfo& info, const Separator& sep )
      {
          const int numChildren = info.children.size();
          for( int index=0; index<numChildren; ++index )
              packRowsLocal( *info.children[index], *sep.children[index] );

          const Int numInds = sep.inds.size();
          for( Int t=0; t<numInds; ++t )
          {
              const Int i = sep.inds[t];
              const Int iMap = sep.off + t;
              const int q = map.RowOwner(i);
              sendOrigInds[offs[q]] = i;
              sendInds[offs[q]] = iMap;
              ++offs[q];
          }
      };
    function<void(const DistNodeInfo&,const DistSeparator&)> packRows =
      [&]( const DistNodeInfo& info, const DistSeparator& sep )
      {
          if( sep.child == nullptr )
          {
              const NodeInfo& infoDup = *info.duplicate;
              const Separator& sepDup = *sep.duplicate;
              const int numChildren = infoDup.children.size();
              for( int index=0; index<numChildren; ++index )
                  packRowsLocal
                  ( *infoDup.children[index], *sepDup.children[index] );
          }
          else
              packRows( *info.child, *sep.child );

          const Int numInds = sep.inds.size();
          const Grid& grid = info.Grid();
          const int teamSize = grid.Size();
          const int teamRank = grid.Rank();
          const Int numLocalInds = Length( numInds, teamRank, teamSize );
          for( Int tLocal=0; tLocal<numLocalInds; ++tLocal )
          {
              const Int t = teamRank + tLocal*teamSize;
              const Int i = sep.inds[t];
              const Int iMap = sep.off + t;
              const int q = map.RowOwner(i);
              sendOrigInds[offs[q]] = i;
              sendInds[offs[q]] = iMap;
              ++offs[q];
          }
      };
    packRows( rootInfo, *this );

    // Perform an AllToAll to exchange the reordered indices
    vector<int> recvOffs;
    const int numRecvs = Scan( recvSizes, recvOffs );
    EL_DEBUG_ONLY(
      const Int numLocalSources = map.NumLocalSources();
      if( numRecvs != numLocalSources )
          LogicError("incorrect number of recv indices");
    )
    vector<Int> recvInds( numRecvs );
    mpi::AllToAll
    ( sendInds.data(), sendSizes.data(), sendOffs.data(),
      recvInds.data(), recvSizes.data(), recvOffs.data(), rootGrid.Comm() );

    // Perform an AllToAll to exchange the original indices
    vector<Int> recvOrigInds( numRecvs );
    mpi::AllToAll
    ( sendOrigInds.data(), sendSizes.data(), sendOffs.data(),
      recvOrigInds.data(), recvSizes.data(), recvOffs.data(), rootGrid.Comm() );

    // Unpack the indices
    const Int firstLocalSource = map.FirstLocalSource();
    auto& mapLoc = map.Map();
    for( Int s=0; s<numRecvs; ++s )
        mapLoc[recvOrigInds[s]-firstLocalSource] = recvInds[s];
}

} // namespace ldl
} // namespace El
