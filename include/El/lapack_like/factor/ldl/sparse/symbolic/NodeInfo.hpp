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

    // This observing pointer is to the parent node (should one exist).
    NodeInfo* parent=nullptr;

    // This observing pointer is to the equivalent distributed node and is only
    // used for the highest node in the sequential subtree. The duplicate will
    // always be 'distributed' over a single process.
    DistNodeInfo* duplicate=nullptr;

    // Unique pointers to any children of this node.
    vector<unique_ptr<NodeInfo>> children;

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

    // This observing pointer is to the parent node (should one exist).
    DistNodeInfo* parent=nullptr;

    // A unique pointer to the sequential equivalent node (should one exist).
    // Such a pointer should exist when this node is distributed over a single
    // process.
    unique_ptr<NodeInfo> duplicate;

    // A unique pointer to the distributed child shared by this process
    // (should one exist).
    unique_ptr<DistNodeInfo> child;

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

    // For constructing the root of the tree.
    explicit DistNodeInfo( const El::Grid& grid );

    // For constructing descendents in the tree.
    DistNodeInfo( DistNodeInfo* parentNode );

    ~DistNodeInfo();

    void SetRootGrid( const El::Grid& grid );
    void AssignGrid( unique_ptr<El::Grid>& grid );
    const El::Grid& Grid() const;

    void GetChildGridDims
    ( vector<int>& gridHeights, vector<int>& gridWidths ) const;

private:
    // If we are the root node, there is no need to construct a grid, so we
    // make use of an observing pointer.
    const El::Grid* rootGrid_=nullptr;

    // If we are not the root node, the following will point to a constructed
    // grid.
    // TODO(poulson): Accept a function for determining the grid dimensions.
    unique_ptr<El::Grid> grid_;
};

} // namespace ldl
} // namespace El

#endif // ifndef EL_FACTOR_LDL_SPARSE_SYMBOLIC_NODEINFO_HPP
