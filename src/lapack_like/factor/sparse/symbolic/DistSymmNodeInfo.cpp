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

DistSymmNodeInfo::DistSymmNodeInfo( DistSymmNodeInfo* parentNode )
: parent(parentNode), child(nullptr), duplicate(nullptr)
{ }

DistSymmNodeInfo::~DistSymmNodeInfo()
{
    if( uncaught_exception() )
    {
        cerr << "Uncaught exception" << endl;
        DEBUG_ONLY(DumpCallStack())
        return;
    }

    delete child;
    delete duplicate;

    delete grid;
    mpi::Free( comm );
}

// TODO: Simplify and/or obsolete this horrendous implementation
void ComputeFactRecvInds( const DistSymmNodeInfo& node )
{
    DEBUG_ONLY(CallStackEntry cse("ComputeFactRecvInds"))
    // Communicate to get the grid sizes
    int childGridDims[4];
    GetChildGridDims( node, childGridDims );

    const int teamSize = mpi::Size( node.comm );
    const int teamRank = mpi::Rank( node.comm );
    const bool onLeft = node.child->onLeft;
    const int childTeamSize = mpi::Size( node.child->comm );
    vector<int> teamSizes(2), gridHeights(2), gridWidths(2);
    teamSizes[0] = ( onLeft ? childTeamSize : teamSize-childTeamSize );
    teamSizes[1] = teamSize - teamSizes[0];
    gridHeights[0] = childGridDims[0];
    gridWidths[0] = childGridDims[1];
    gridHeights[1] = childGridDims[2];
    gridWidths[1] = childGridDims[3];
    DEBUG_ONLY(
        if( teamSizes[0] != gridHeights[0]*gridWidths[0] )
            RuntimeError("Computed left grid incorrectly");
        if( teamSizes[1] != gridHeights[1]*gridWidths[1] )
            RuntimeError("Computed right grid incorrectly");
    )

    const FactorCommMeta& commMeta = node.factorMeta;
    commMeta.childRecvInds.resize( teamSize );

    const int gridHeight = node.grid->Height();
    const int gridWidth = node.grid->Width();
    const int gridRow = node.grid->Row();
    const int gridCol = node.grid->Col();

    const int childTeamRank = mpi::Rank( node.child->comm );
    const bool inFirstTeam = ( childTeamRank == teamRank );
    const bool leftIsFirst = ( onLeft==inFirstTeam );
    vector<int> teamOffs(2);
    teamOffs[0] = ( leftIsFirst ? 0            : teamSizes[1] );
    teamOffs[1] = ( leftIsFirst ? teamSizes[0] : 0            );

    // TODO: Loop over {0,1} instead of explicitly handling left and right
    for( Int c=0; c<2; ++c )
    {
        // Compute the recv indices of the child from each process 
        const Int numInds = node.childRelInds[c].size();
        vector<Int> rowInds, colInds;
        for( Int i=0; i<numInds; ++i )
            if( node.childRelInds[c][i] % gridHeight == gridRow )
                rowInds.push_back( i );
        for( Int i=0; i<numInds; ++i )
            if( node.childRelInds[c][i] % gridWidth == gridCol )
                colInds.push_back( i );

        vector<Int>::const_iterator it;
        const Int numColInds = colInds.size();
        const Int numRowInds = rowInds.size();
        for( Int jPre=0; jPre<numColInds; ++jPre )
        {
            const Int jChild = colInds[jPre];
            const Int j = node.childRelInds[c][jChild];
            DEBUG_ONLY(
                if( (j-gridCol) % gridWidth != 0 )
                    LogicError("Invalid j");
            )
            const Int jLoc = (j-gridCol) / gridWidth;
            const int childCol = (jChild+node.childSizes[c]) % gridWidths[c];

            // Find the first iPre that maps to the lower triangle
            it = std::lower_bound( rowInds.begin(), rowInds.end(), jChild );
            const Int iPreStart = Int(it-rowInds.begin());
            for( Int iPre=iPreStart; iPre<numRowInds; ++iPre )
            {
                const Int iChild = rowInds[iPre];
                const Int i = node.childRelInds[c][iChild];
                DEBUG_ONLY(
                    if( iChild < jChild )
                        LogicError("Invalid iChild");
                    if( (i-gridRow) % gridHeight != 0 )
                        LogicError("Invalid  i");
                )
                const Int iLoc = (i-gridRow) / gridHeight;

                const int childRow = (iChild+node.childSizes[c])%gridHeights[c];
                const int childRank = childRow+childCol*gridHeights[c];

                const int q = teamOffs[c]+ childRank;
                commMeta.childRecvInds[q].push_back(iLoc);
                commMeta.childRecvInds[q].push_back(jLoc);
            }
        }
    }
}

void GetChildGridDims( const DistSymmNodeInfo& node, int* childGridDims )
{
    const bool onLeft = node.child->onLeft;
    const int childTeamRank = mpi::Rank( node.child->comm );
    El::MemZero( childGridDims, 4 );
    if( onLeft && childTeamRank == 0 )
    {
        childGridDims[0] = node.child->grid->Height();
        childGridDims[1] = node.child->grid->Width();
    }
    else if( !onLeft && childTeamRank == 0 )
    {
        childGridDims[2] = node.child->grid->Height();
        childGridDims[3] = node.child->grid->Width();
    }
    mpi::AllReduce( childGridDims, 4, mpi::SUM, node.comm );
}

} // namespace El
