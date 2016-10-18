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
#include <El.hpp>

namespace El {
namespace ldl {

template<typename T>
DistMatrixNode<T>::DistMatrixNode( DistMatrixNode<T>* parentNode )
: parent(parentNode), child(nullptr), duplicate(nullptr)
{ }

template<typename T>
DistMatrixNode<T>::DistMatrixNode
( const DistMap& invMap,
  const DistNodeInfo& info,
  const DistMultiVec<T>& X )
: parent(nullptr), child(nullptr), duplicate(nullptr)
{
    DEBUG_CSE
    Pull( invMap, info, X );
}

template<typename T>
DistMatrixNode<T>::DistMatrixNode( const DistMultiVecNode<T>& X )
: parent(nullptr), child(nullptr), duplicate(nullptr)
{
    DEBUG_CSE
    *this = X;
}

template<typename T>
DistMatrixNode<T>::~DistMatrixNode()
{
    delete child;
    delete duplicate;
}

template<typename T>
const DistMatrixNode<T>&
DistMatrixNode<T>::operator=( const DistMultiVecNode<T>& X )
{
    DEBUG_CSE

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
( const DistMap& invMap,
  const DistNodeInfo& info,
  const DistMultiVec<T>& X )
{
    DEBUG_CSE
    DistMultiVecNode<T> XMultiVec( invMap, info, X );
    *this = XMultiVec;
}

template<typename T>
void DistMatrixNode<T>::Push
( const DistMap& invMap,
  const DistNodeInfo& info,
        DistMultiVec<T>& X ) const
{
    DEBUG_CSE
    DistMultiVecNode<T> XMultiVec( *this );
    XMultiVec.Push( invMap, info, X );
}

// NOTE: It is assumed that the child 'work' matrix is already filled
template<typename T>
void DistMatrixNode<T>::ComputeCommMeta( const DistNodeInfo& info ) const
{
    DEBUG_CSE
    if( commMeta.numChildSendInds.size() != 0 )
        return;
    if( child == nullptr )
        return;
    
    const Int numRHS = matrix.Width();
    const Int childSize = info.child->size;
    const Int updateSize = info.child->lowerStruct.size();
    const Int workSize = childSize + updateSize;
    
    auto childWT = child->work( IR(0,childSize),        IR(0,numRHS) );
    auto childWB = child->work( IR(childSize,workSize), IR(0,numRHS) );
    DEBUG_ONLY(
      if( matrix.ColAlign() != 0 || matrix.RowAlign() != 0 )
          LogicError
          ("matrix was not zero aligned: ",
           matrix.ColAlign()," ",matrix.RowAlign());
      if( childWT.ColAlign() != 0 || childWT.RowAlign() != 0 )
          LogicError
          ("childWT was not zero aligned: ",
           childWT.ColAlign()," ",childWT.RowAlign());
      if( child->work.Height() != workSize || child->work.Width() != numRHS )
          LogicError
          ("childWork was not the expected size: ",
           child->work.Height()," x ",child->work.Width());
    )

    // Fill numChildSendInds
    // =====================
    const int teamSize = mpi::Size( info.comm );
    commMeta.numChildSendInds.resize( teamSize );
    MemZero( commMeta.numChildSendInds.data(), teamSize );
    const Int myChild = ( info.child->onLeft ? 0 : 1 );
    const Int childLocalHeight = childWB.LocalHeight();
    const Int childLocalWidth = childWB.LocalWidth();
    for( Int iChildLoc=0; iChildLoc<childLocalHeight; ++iChildLoc )
    {
        const Int iChild = childWB.GlobalRow(iChildLoc);
        const Int iParent = info.childRelInds[myChild][iChild];
        for( Int jChildLoc=0; jChildLoc<childLocalWidth; ++jChildLoc )
        {
            const Int j = childWB.GlobalCol(jChildLoc);
            const int q = matrix.Owner(iParent,j);
            ++commMeta.numChildSendInds[q];
        }
    }

    // Compute the solve recv indices
    // ==============================
    vector<int> gridHeights, gridWidths;
    GetChildGridDims( info, gridHeights, gridWidths );

    const int teamRank = mpi::Rank( info.comm );
    const bool onLeft = info.child->onLeft;
    const int childTeamSize = mpi::Size( info.child->comm );
    const int childTeamRank = mpi::Rank( info.child->comm );
    const bool inFirstTeam = ( childTeamRank == teamRank );
    const bool leftIsFirst = ( onLeft==inFirstTeam );
    vector<int> teamSizes(2), teamOffs(2);
    teamSizes[0] = ( onLeft ? childTeamSize : teamSize-childTeamSize );
    teamSizes[1] = teamSize - teamSizes[0];
    teamOffs[0] = ( leftIsFirst ? 0            : teamSizes[1] );
    teamOffs[1] = ( leftIsFirst ? teamSizes[0] : 0            );
    DEBUG_ONLY(
      if( teamSizes[0] != gridHeights[0]*gridWidths[0] )
          RuntimeError("Computed left grid incorrectly");
      if( teamSizes[1] != gridHeights[1]*gridWidths[1] )
          RuntimeError("Computed right grid incorrectly");
    )

    const Int localWidth = matrix.LocalWidth();
    commMeta.childRecvInds.resize( teamSize );
    for( int q=0; q<teamSize; ++q )
        commMeta.childRecvInds[q].clear();
    for( Int c=0; c<2; ++c )
    {
        const Int numInds = info.childRelInds[c].size();
        vector<Int> rowInds;
        for( Int iChild=0; iChild<numInds; ++iChild )
            if( matrix.IsLocalRow( info.childRelInds[c][iChild] ) )
                rowInds.push_back( iChild );

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

#define PROTO(T) template struct DistMatrixNode<T>;
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace ldl
} // namespace El
