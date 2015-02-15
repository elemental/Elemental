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

template<typename T>
DistMatrixNode<T>::DistMatrixNode( DistMatrixNode<T>* parentNode )
: parent(parentNode), child(nullptr), duplicate(nullptr)
{ }

template<typename T>
DistMatrixNode<T>::DistMatrixNode
( const DistMap& invMap, const DistSymmNodeInfo& info,
  const DistMultiVec<T>& X )
: parent(nullptr), child(nullptr), duplicate(nullptr)
{
    DEBUG_ONLY(CallStackEntry cse("DistMatrixNode::DistMatrixNode"))
    Pull( invMap, info, X );
}

template<typename T>
DistMatrixNode<T>::DistMatrixNode( const DistMultiVecNode<T>& X )
: parent(nullptr), child(nullptr), duplicate(nullptr)
{
    DEBUG_ONLY(CallStackEntry cse("DistMatrixNode::DistMatrixNode"))
    *this = X;
}

template<typename T>
const DistMatrixNode<T>&
DistMatrixNode<T>::operator=( const DistMultiVecNode<T>& X )
{
    DEBUG_ONLY(CallStackEntry cse("DistMatrixNode::operator="))

    if( X.child == nullptr )
    {
        delete duplicate;
        duplicate = new MatrixNode<T>(this);
        *duplicate = *X.duplicate;

        matrix.Attach( X.matrix.Grid(), duplicate->matrix );

        return *this;
    }

    commMeta.Empty();
    matrix.SetGrid( X.matrix.Grid() );
    matrix = X.matrix;

    delete child;
    child = new DistMatrixNode<T>(this);
    *child = *X.child;

    return *this;
}

template<typename T>
void DistMatrixNode<T>::Pull
( const DistMap& invMap, const DistSymmNodeInfo& info,
  const DistMultiVec<T>& X )
{
    DEBUG_ONLY(CallStackEntry cse("DistMatrixNode::Pull"))
    DistMultiVecNode<T> XMultiVec( invMap, info, X );
    *this = XMultiVec;
}

template<typename T>
void DistMatrixNode<T>::Push
( const DistMap& invMap, const DistSymmNodeInfo& info,
        DistMultiVec<T>& X ) const
{
    DEBUG_ONLY(CallStackEntry cse("DistMatrixNode::Push"))
    DistMultiVecNode<T> XMultiVec( *this );
    XMultiVec.Push( invMap, info, X );
}

// NOTE: It is assumed that the child 'work' matrix is already filled
template<typename T>
void DistMatrixNode<T>::ComputeCommMeta( const DistSymmNodeInfo& info ) const
{
    DEBUG_ONLY(CallStackEntry cse("DistMatrixNode::ComputeCommMetas"))
    if( commMeta.numChildSendInds.size() != 0 )
        return;
    if( child == nullptr )
        return;
   
    const int teamSize = mpi::Size( info.comm );

    const Int numRHS = matrix.Width();
    const Int childSize = info.child->size;
    const Int updateSize = info.child->lowerStruct.size();
    const Int workSize = childSize + updateSize;
    
    auto childWT = child->work( IR(0,childSize),        IR(0,numRHS) );
    auto childWB = child->work( IR(childSize,workSize), IR(0,numRHS) );

    // Fill numChildSendInds
    // =====================
    commMeta.Empty();
    commMeta.numChildSendInds.resize( teamSize );
    MemZero( commMeta.numChildSendInds.data(), teamSize );
    const auto& childRelInds =
      ( info.child->onLeft ? info.childRelInds[0] : info.childRelInds[1] );
    const Int localHeight = childWB.LocalHeight();
    const Int localWidth = childWB.LocalWidth();
    for( Int iChildLoc=0; iChildLoc<localHeight; ++iChildLoc )
    {
        const Int iChild = childWB.GlobalRow(iChildLoc);
        const Int iParent = childRelInds[iChild];
        for( Int jChildLoc=0; jChildLoc<localWidth; ++jChildLoc )
        {
            const Int j = childWB.GlobalCol(jChildLoc);
            const int q = matrix.Owner(iParent,j);
            ++commMeta.numChildSendInds[q];
        }
    }

    // Compute the solve recv indices
    // ==============================
    const int teamRank = mpi::Rank( info.comm );
    const int childTeamSize = mpi::Size( info.child->comm );
    const int childTeamRank = mpi::Rank( info.child->comm );
    const bool inFirstTeam = ( childTeamRank == teamRank );
    const bool leftIsFirst = ( info.child->onLeft==inFirstTeam );
    vector<int> teamSizes(2), teamOffs(2);
    teamSizes[0] = info.child->onLeft ? childTeamSize : teamSize-childTeamSize;
    teamSizes[1] = teamSize - teamSizes[0];
    teamOffs[0] = ( leftIsFirst ? 0            : teamSizes[1] );
    teamOffs[1] = ( leftIsFirst ? teamSizes[0] : 0            );

    // Get the child grid dimensions
    vector<int> gridHeights, gridWidths;
    GetChildGridDims( info, gridHeights, gridWidths );

    commMeta.childRecvInds.resize( teamSize );
    for( int q=0; q<teamSize; ++q )
        commMeta.childRecvInds[q].clear();
    for( Int c=0; c<2; ++c )
    {
        const Int numInds = info.childRelInds[c].size();
        vector<Int> rowInds;
        for( Int i=0; i<numInds; ++i )
            if( matrix.IsLocalRow( info.childRelInds[c][i] ) )
                rowInds.push_back( i );

        const Int numRowInds = rowInds.size();
        for( Int iPre=0; iPre<numRowInds; ++iPre )
        {
            const Int iChild = rowInds[iPre];
            const Int i = info.childRelInds[c][iChild];
            const Int iLoc = matrix.LocalRow(i);
            const int childRow = (info.childSizes[c]+iChild) % gridHeights[c];
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const Int j = matrix.GlobalCol(jLoc);
                const int childCol = j % gridWidths[c];
                const int childRank = childRow + childCol*gridHeights[c];
                const int q = teamOffs[c] + childRank;
                commMeta.childRecvInds[q].push_back(iLoc);
                commMeta.childRecvInds[q].push_back(jLoc);
            }
        }
    }
}

#define PROTO(T) template class DistMatrixNode<T>;
#include "El/macros/Instantiate.h"

} // namespace El
