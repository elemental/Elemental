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
#pragma once
#ifndef EL_FACTOR_LDL_SPARSE_SYMBOLIC_NODE_HPP
#define EL_FACTOR_LDL_SPARSE_SYMBOLIC_NODE_HPP

namespace El {
namespace ldl {

// 'Supernode' should perhaps be preferred to 'node', but since we will always
// use supernodes, the extra verbage is unnecessarily cumbersome.

struct DistNode;

struct Node
{
    Int size, off;
    vector<Int> lowerStruct;

    Node* parent; 
    vector<Node*> children;
    DistNode* duplicate;

    Node( Node* parentNode=nullptr )
    : parent(parentNode), duplicate(nullptr)
    { }

    Node( DistNode* duplicateNode );

    ~Node()
    {
        if( uncaught_exception() )
        {
            cerr << "Uncaught exception" << endl;
            DEBUG_ONLY(DumpCallStack())
            return;
        }

        for( const Node* child : children )
            delete child;
    }
};

struct DistNode
{
    bool onLeft; // irrelevant if root node
    mpi::Comm comm;
    Int size, off;
    vector<Int> lowerStruct;

    DistNode* parent;
    DistNode* child;
    Node* duplicate;

    DistNode( DistNode* parentNode=nullptr )
    : comm(mpi::COMM_WORLD), 
      parent(parentNode), child(nullptr), duplicate(nullptr)
    { }

    ~DistNode()
    {
        if( uncaught_exception() )
        {
            cerr << "Uncaught exception" << endl;
            DEBUG_ONLY(DumpCallStack())
            return;
        }

        delete child;
        delete duplicate;

        if( comm != mpi::COMM_WORLD )
            mpi::Free( comm );
    }
};

inline Node::Node( DistNode* duplicateNode )
: parent(nullptr), duplicate(duplicateNode)
{
    size = duplicate->size;
    off = duplicate->off;
    lowerStruct = duplicate->lowerStruct;
}

} // namespace ldl
} // namespace El

#endif // ifndef EL_FACTOR_LDL_SPARSE_SYMBOLIC_NODE_HPP
