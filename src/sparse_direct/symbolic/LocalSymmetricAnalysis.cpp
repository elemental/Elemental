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

void LocalSymmetricAnalysis( const DistSymmElimTree& eTree, DistSymmInfo& info )
{
    DEBUG_ONLY(CallStackEntry cse("LocalSymmetricAnalysis"))
    const Int numNodes = eTree.localNodes.size();
    info.localNodes.resize( numNodes );

    // Perform the symbolic factorization
    Int myOff= 0;
    for( Int s=0; s<numNodes; ++s )
    {
        const SymmNode& node = *eTree.localNodes[s];
        SymmNodeInfo& nodeInfo = info.localNodes[s];
        nodeInfo.size = node.size;
        nodeInfo.off= node.off;
        nodeInfo.myOff= myOff;
        nodeInfo.parent = node.parent;
        nodeInfo.children = node.children;
        nodeInfo.origLowerStruct = node.lowerStruct;

        const Int numChildren = node.children.size();
        DEBUG_ONLY(
            if( numChildren != 0 && numChildren != 2 )
                LogicError("Tree must be built from bisections");
        )
        if( numChildren == 2 )
        {
            const Int left = node.children[0];
            const Int right = node.children[1];
            SymmNodeInfo& leftChild = info.localNodes[left];
            SymmNodeInfo& rightChild = info.localNodes[right];
            leftChild.onLeft = true;
            rightChild.onLeft = false;
            DEBUG_ONLY(
                if( !IsStrictlySorted(leftChild.lowerStruct) )
                {
                    if( IsSorted(leftChild.lowerStruct) )
                        LogicError("Repeat in left lower struct");
                    else
                        LogicError("Left lower struct not sorted");
                }
                if( !IsStrictlySorted(rightChild.lowerStruct) )
                {
                    if( IsSorted(rightChild.lowerStruct) )
                        LogicError("Repeat in right lower struct");
                    else
                        LogicError("Right lower struct not sorted");
                }
                if( !IsStrictlySorted(node.lowerStruct) )
                {
                    if( IsSorted(node.lowerStruct) )
                        LogicError("Repeat in original lower struct");
                    else
                        LogicError("Original lower struct not sorted");
                }
            )

            // Combine the structures of the children
            auto childrenStruct = 
                Union( leftChild.lowerStruct, rightChild.lowerStruct );

            // Now add in the original lower structure
            auto partialStruct = Union( node.lowerStruct, childrenStruct );

            // Now the node indices
            std::vector<Int> nodeInds( node.size );
            for( Int i=0; i<node.size; ++i )
                nodeInds[i] = node.off+ i;
            auto fullStruct = Union( partialStruct, nodeInds );

            // Construct the relative indices of the original lower structure
            nodeInfo.origLowerRelInds = 
                RelativeIndices( node.lowerStruct, fullStruct );

            // Construct the relative indices of the children
            nodeInfo.leftRelInds = 
                RelativeIndices( leftChild.lowerStruct, fullStruct );
            nodeInfo.rightRelInds =
                RelativeIndices( rightChild.lowerStruct, fullStruct );

            // Form lower struct of this node by removing node indices
            // (which take up the first node.size indices of fullStruct)
            const Int lowerStructSize = fullStruct.size()-node.size;
            nodeInfo.lowerStruct.resize( lowerStructSize );
            for( Int i=0; i<lowerStructSize; ++i )
                nodeInfo.lowerStruct[i] = fullStruct[node.size+i];
        }
        else // numChildren == 0, so this is a leaf node 
        {
            nodeInfo.lowerStruct = node.lowerStruct;
            
            // Construct the trivial relative indices of the original structure
            const Int numOrigLowerInds = node.lowerStruct.size();
            nodeInfo.origLowerRelInds.resize( numOrigLowerInds );
            for( Int i=0; i<numOrigLowerInds; ++i )
                nodeInfo.origLowerRelInds[i] = i + nodeInfo.size;
        }

        myOff += nodeInfo.size;
    }
}

} // namespace El
