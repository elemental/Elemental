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
#ifndef EL_SYMBOLIC_DISTSEPARATORTREE_HPP
#define EL_SYMBOLIC_DISTSEPARATORTREE_HPP

namespace El {

// TODO: Rename to SeparatorOrLeaf?
struct SepOrLeaf
{
    int parent; // -1 if local root
    int off;
    std::vector<int> inds;
};

struct DistSeparator
{
    mpi::Comm comm;
    int off;
    std::vector<int> inds;
};

struct DistSeparatorTree
{
    // Full local binary tree
    //
    // NOTE: This is an array of pointers, as we will not know during 
    //       construction how many will need to be created
    std::vector<SepOrLeaf*> localSepsAndLeaves;

    // One path through top of binary tree 
    //
    // NOTE: does not include the single-process separator/leaf
    std::vector<DistSeparator> distSeps;

    ~DistSeparatorTree()
    {
        if( std::uncaught_exception() )
        {
            std::cerr << "Uncaught exception in ~DistSepTree" << std::endl;
            DEBUG_ONLY(DumpCallStack())
            return;
        }

        const int numLocal = localSepsAndLeaves.size();
        for( int i=0; i<numLocal; ++i )
            delete localSepsAndLeaves[i];

        const int numDist = distSeps.size();
        for( int i=0; i<numDist; ++i )
            mpi::Free( distSeps[i].comm );
    }
};

} // namespace El

#endif // ifndef EL_SYMBOLIC_DISTSEPARATORTREE_HPP
