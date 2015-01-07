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
#ifndef EL_SYMBOLIC_DISTSYMMELIMTREE_HPP
#define EL_SYMBOLIC_DISTSYMMELIMTREE_HPP

namespace El {

// 'Supernode' should perhaps be preferred to 'node', but since we will always
// use supernodes, the extra verbage is unnecessarily cumbersome.

struct SymmNode
{
    Int size, off; 
    Int parent; // -1 if root separator
    std::vector<Int> children;
    std::vector<Int> lowerStruct;
};

struct DistSymmNode
{
    bool onLeft; // irrelevant if root node
    mpi::Comm comm;
    Int size, off;
    std::vector<Int> lowerStruct;
};

struct DistSymmElimTree
{
    // NOTE: This is an array of pointers, as we will not know how many 
    //       are needed during construction
    std::vector<SymmNode*> localNodes;
    std::vector<DistSymmNode> distNodes;

    ~DistSymmElimTree()
    {
        if( std::uncaught_exception() )
        {
            std::cerr << "Uncaught exception in ~DistSymmElimTree" << std::endl;
            DEBUG_ONLY(DumpCallStack())
            return;
        }

        const Int numLocal = localNodes.size();
        for( Int i=0; i<numLocal; ++i )
            delete localNodes[i];

        const Int numDist = distNodes.size();
        for( Int i=0; i<numDist; ++i )
            mpi::Free( distNodes[i].comm );
    }
};

} // namespace El

#endif // ifndef EL_SYMBOLIC_DISTSYMMELIMTREE_HPP
