/*
   Copyright (c) 2009-2016, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, Stanford University, and the
   Georgia Insitute of Technology.
   All rights reserved.
 
   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_FACTOR_LDL_SPARSE_SYMBOLIC_NODEINFO_HPP
#define EL_FACTOR_LDL_SPARSE_SYMBOLIC_NODEINFO_HPP

namespace El {
namespace ldl {

struct DistNodeInfo;

struct NodeInfo
{
    // Known before analysis
    // ---------------------
    Int size, off; 
    vector<Int> origLowerStruct;

    NodeInfo* parent;
    vector<NodeInfo*> children;
    DistNodeInfo* duplicate; 

    // Known after analysis
    // --------------------
    Int myOff;
    vector<Int> lowerStruct;
    vector<Int> origLowerRelInds;
    // (maps from the child update indices to our frontal indices).
    vector<vector<Int>> childRelInds;

    // Symbolic analysis for modification of SuiteSparse LDL
    // -----------------------------------------------------
    // NOTE: These are only used within leaf nodes
    vector<Int> LOffsets;
    vector<Int> LParents;

    NodeInfo( NodeInfo* parentNode=nullptr )
    : parent(parentNode), duplicate(nullptr)
    { }

    NodeInfo( DistNodeInfo* duplicateNode );

    ~NodeInfo()
    {
        if( uncaught_exception() )
        {
            cerr << "Uncaught exception" << endl;
            DEBUG_ONLY(DumpCallStack())
            return;
        }

        for( const NodeInfo* child : children )
            delete child;
    }
};

struct DistNodeInfo
{
    // Known before analysis
    // ---------------------
    Int size, off;
    vector<Int> origLowerStruct;
    bool onLeft;
    mpi::Comm comm;

    DistNodeInfo* parent; 
    DistNodeInfo* child; 
    NodeInfo* duplicate;

    // Known after analysis
    // --------------------
    Grid* grid;
    Int myOff;
    vector<Int> lowerStruct;
    vector<Int> origLowerRelInds;

    vector<Int> childSizes;
    // The relative indices of our children
    // (maps from the child update indices to our frontal indices).
    // These could be replaced with just the relative indices of our local 
    // submatrices of the child updates.
    vector<vector<Int>> childRelInds;

    DistNodeInfo( DistNodeInfo* parentNode=nullptr )
    : comm(mpi::COMM_WORLD), 
      parent(parentNode), child(nullptr), duplicate(nullptr),
      grid(nullptr)
    { }

    ~DistNodeInfo()
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
        if( comm != mpi::COMM_WORLD )
            mpi::Free( comm );
    }   
};

inline NodeInfo::NodeInfo( DistNodeInfo* duplicateNode )
: parent(nullptr), duplicate(duplicateNode)
{
    size = duplicate->size;
    off = duplicate->off;
    origLowerStruct = duplicate->origLowerStruct;

    myOff = duplicate->myOff;
    lowerStruct = duplicate->lowerStruct;
    origLowerRelInds = duplicate->origLowerRelInds;
}

} // namespace ldl
} // namespace El

#endif // ifndef EL_FACTOR_LDL_SPARSE_SYMBOLIC_NODEINFO_HPP
