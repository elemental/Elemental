/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   Copyright (c) 2012 Jack Poulson, Lexing Ying, and
   The University of Texas at Austin.
   All rights reserved.

   Copyright (c) 2013 Jack Poulson, Lexing Ying, and Stanford University.
   All rights reserved.

   Copyright (c) 2014 Jack Poulson and The Georgia Institute of Technology.
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

    NodeInfo* parent=nullptr;
    vector<NodeInfo*> children;
    DistNodeInfo* duplicate=nullptr;

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

    NodeInfo( NodeInfo* parentNode=nullptr );
    NodeInfo( DistNodeInfo* duplicateNode );
    ~NodeInfo();
};

struct DistNodeInfo
{
    // Known before analysis
    // ---------------------
    Int size, off;
    vector<Int> origLowerStruct;
    bool onLeft;

    DistNodeInfo* parent=nullptr;
    DistNodeInfo* child=nullptr;
    NodeInfo* duplicate=nullptr;

    const Grid* grid=nullptr;

    // Known after analysis
    // --------------------
    Int myOff;
    vector<Int> lowerStruct;
    vector<Int> origLowerRelInds;

    vector<Int> childSizes;
    // The relative indices of our children
    // (maps from the child update indices to our frontal indices).
    // These could be replaced with just the relative indices of our local
    // submatrices of the child updates.
    vector<vector<Int>> childRelInds;

    DistNodeInfo( DistNodeInfo* parentNode=nullptr );
    ~DistNodeInfo();

    void GetChildGridDims
    ( vector<int>& gridHeights, vector<int>& gridWidths ) const;
};

} // namespace ldl
} // namespace El

#endif // ifndef EL_FACTOR_LDL_SPARSE_SYMBOLIC_NODEINFO_HPP
