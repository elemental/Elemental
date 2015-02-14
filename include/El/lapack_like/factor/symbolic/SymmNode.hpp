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
#ifndef EL_SYMBOLIC_SYMMNODE_HPP
#define EL_SYMBOLIC_SYMMNODE_HPP

namespace El {

// 'Supernode' should perhaps be preferred to 'node', but since we will always
// use supernodes, the extra verbage is unnecessarily cumbersome.

struct DistSymmNode;

struct SymmNode
{
    Int size, off; 
    vector<Int> lowerStruct;

    SymmNode* parent; 
    vector<SymmNode*> children;
    DistSymmNode* duplicate; 

    SymmNode( SymmNode* parentNode=nullptr );
    SymmNode( DistSymmNode* duplicateNode );
    ~SymmNode();
};

struct DistSymmNode
{
    bool onLeft; // irrelevant if root node
    mpi::Comm comm;
    Int size, off;
    vector<Int> lowerStruct;

    DistSymmNode* parent; 
    DistSymmNode* child; 
    SymmNode* duplicate;

    DistSymmNode( DistSymmNode* parentNode=nullptr );
    ~DistSymmNode();
};

} // namespace El

#endif // ifndef EL_SYMBOLIC_SYMMNODE_HPP
