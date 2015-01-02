/*
   Copyright (c) 2009-2015, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, Stanford University, and the
   Georgia Insitute of Technology.
   All rights reserved.
 
   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_SPARSEDIRECT_NUMERIC_LOWERMULTIPLY_LOCAL_HPP
#define EL_SPARSEDIRECT_NUMERIC_LOWERMULTIPLY_LOCAL_HPP

namespace El {

template<typename T> 
inline void LocalLowerMultiplyNormal
( int diagOff, const DistSymmInfo& info, 
  const DistSymmFrontTree<T>& L, DistNodalMultiVec<T>& X )
{
    DEBUG_ONLY(CallStackEntry cse("LocalLowerMultiplyNormal"))
    const int numLocalNodes = info.localNodes.size();
    const int width = X.Width();
    for( int s=0; s<numLocalNodes; ++s )
    {
        const SymmNodeInfo& node = info.localNodes[s];
        const Matrix<T>& frontL = L.localFronts[s].frontL;
        Matrix<T>& W = L.localFronts[s].work;

        // Set up a workspace
        W.Resize( frontL.Height(), width );
        Matrix<T> WT, WB;
        PartitionDown( W, WT, WB, node.size );
        WT = X.localNodes[s];
        Zero( WB );

        // Multiply this block column of L against this node's portion of the
        // right-hand side and set W equal to the result
        FrontLowerMultiply( NORMAL, diagOff, frontL, W );

        // Update using the children (if they exist)
        const int numChildren = node.children.size();
        if( numChildren == 2 )
        {
            const int leftInd = node.children[0];
            const int rightInd = node.children[1];
            Matrix<T>& leftWork = L.localFronts[leftInd].work;
            Matrix<T>& rightWork = L.localFronts[rightInd].work;
            const int leftNodeSize = info.localNodes[leftInd].size;
            const int rightNodeSize = info.localNodes[rightInd].size;
            const int leftUpdateSize = leftWork.Height()-leftNodeSize;
            const int rightUpdateSize = rightWork.Height()-rightNodeSize;

            // Add the left child's update onto ours
            auto leftUpdate =
                LockedView( leftWork, leftNodeSize, 0, leftUpdateSize, width );
            for( int iChild=0; iChild<leftUpdateSize; ++iChild )
            {
                const int iFront = node.leftRelInds[iChild];
                for( int j=0; j<width; ++j )
                    W.Update( iFront, j, leftUpdate.Get(iChild,j) );
            }
            leftWork.Empty();

            // Add the right child's update onto ours
            auto rightUpdate =
                LockedView
                ( rightWork, rightNodeSize, 0, rightUpdateSize, width );
            for( int iChild=0; iChild<rightUpdateSize; ++iChild )
            {
                const int iFront = node.rightRelInds[iChild];
                for( int j=0; j<width; ++j )
                    W.Update( iFront, j, rightUpdate.Get(iChild,j) );
            }
            rightWork.Empty();
        }
        // else numChildren == 0 

        // Store this node's portion of the result
        X.localNodes[s] = WT;
    }
}

template<typename T> 
inline void LocalLowerMultiplyTranspose
( int diagOff, const DistSymmInfo& info, 
  const DistSymmFrontTree<T>& L, DistNodalMultiVec<T>& X,
  bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("LocalLowerMultiplyTranspose"))
    const int numLocalNodes = info.localNodes.size();
    const int width = X.Width();
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );

    for( int s=numLocalNodes-2; s>=0; --s )
    {
        const SymmNodeInfo& node = info.localNodes[s];
        const Matrix<T>& frontL = L.localFronts[s].frontL;
        Matrix<T>& W = L.localFronts[s].work;

        // Set up a workspace
        W.Resize( frontL.Height(), width );
        Matrix<T> WT, WB;
        PartitionDown( W, WT, WB, node.size );
        WT = X.localNodes[s];

        // Update using the parent's portion of the RHS
        const int parent = node.parent;
        const SymmNodeInfo& parentNode = info.localNodes[parent];
        Matrix<T>& parentWork = L.localFronts[parent].work;
        const int currentUpdateSize = WB.Height();
        const std::vector<int>& parentRelInds = 
            ( node.onLeft ? parentNode.leftRelInds : parentNode.rightRelInds );
        for( int iCurrent=0; iCurrent<currentUpdateSize; ++iCurrent )
        {
            const int iParent = parentRelInds[iCurrent]; 
            for( int j=0; j<width; ++j )
                WB.Set( iCurrent, j, parentWork.Get(iParent,j) );
        }

        // The left child is numbered lower than the right child, so we can 
        // safely free the parent's work if this node is the left child
        if( node.onLeft )
        {
            parentWork.Empty();
            if( parent == numLocalNodes-1 )
                L.distFronts[0].work1d.Empty();
        }

        // Make a copy of the unmodified RHS
        Matrix<T> XNode = W;

        // Multiply the (conjugate-)transpose of this block column of L against
        // this node's portion of the right-hand side.
        FrontLowerMultiply( orientation, diagOff, frontL, XNode );

        // Store this node's portion of the result
        Matrix<T> XNodeT, XNodeB;
        PartitionDown( XNode, XNodeT, XNodeB, node.size );
        X.localNodes[s] = XNodeT;
        XNode.Empty();
    }
    L.distFronts[0].work1d.Empty();
    L.localFronts.front().work.Empty();
}

} // namespace El

#endif // ifndef EL_SPARSEDIRECT_NUMERIC_LOWERMULTIPLY_LOCAL_HPP
