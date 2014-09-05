/*
   Copyright (c) 2009-2014, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, Stanford University, and the
   Georgia Insitute of Technology.
   All rights reserved.
 
   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_SYMBOLIC_DISTSYMMINFO_HPP
#define EL_SYMBOLIC_DISTSYMMINFO_HPP

namespace El {

struct SymmNodeInfo
{
    //
    // This is known before analysis
    //
    int size, off; 
    int parent; // -1 if root separator
    std::vector<int> children;
    std::vector<int> origLowerStruct;

    //
    // The following is computed during analysis
    //
    bool onLeft;
    int myOff;
    std::vector<int> lowerStruct;
    std::vector<int> origLowerRelInds;
    // (maps from the child update indices to our frontal indices).
    std::vector<int> leftRelInds, rightRelInds;
};

struct FactorCommMeta
{
    std::vector<int> numChildSendInds;
    // This information does not necessarily have to be kept and can be
    // computed from the above information (albeit somewhat expensively).
    mutable std::vector<std::vector<int>> childRecvInds;

    void EmptyChildRecvIndices() const
    { SwapClear(childRecvInds); }

    void Empty()
    {
        SwapClear( numChildSendInds );
        EmptyChildRecvIndices();
    }
};

struct MultiVecCommMeta
{
    int localOff, localSize;
    std::vector<int> numChildSendInds;
    std::vector<std::vector<int>> childRecvInds;

    void Empty()
    {
        SwapClear( numChildSendInds );
        SwapClear( childRecvInds );
    }
};

struct MatrixCommMeta
{
    std::vector<int> numChildSendInds;
    std::vector<std::vector<int>> childRecvInds;

    void Empty()
    {
        SwapClear( numChildSendInds );
        SwapClear( childRecvInds );
    }
};

struct DistSymmNodeInfo
{
    //
    // This is known before analysis
    //
    int size, off;
    std::vector<int> origLowerStruct;
    bool onLeft;
    mpi::Comm comm;

    //
    // The following is computed during analysis
    //
    Grid* grid;
    int myOff, leftSize, rightSize;
    std::vector<int> lowerStruct;
    std::vector<int> origLowerRelInds;

    // The relative indices of our child
    // (maps from the child update indices to our frontal indices).
    // These could be replaced with just the relative indices of our local 
    // submatrices of the child updates.
    std::vector<int> leftRelInds, rightRelInds;

    FactorCommMeta factorMeta;
    MultiVecCommMeta multiVecMeta;
};

struct DistSymmInfo
{
    std::vector<SymmNodeInfo> localNodes;
    std::vector<DistSymmNodeInfo> distNodes;
    ~DistSymmInfo();
};

// Utilities
void ComputeFactRecvInds
( const DistSymmNodeInfo& node, const DistSymmNodeInfo& childNode );
void GetChildGridDims
( const DistSymmNodeInfo& node, const DistSymmNodeInfo& childNode,
  int* childGridDims );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

inline
DistSymmInfo::~DistSymmInfo()
{
    const int numDist = distNodes.size();
    for( int s=0; s<numDist; ++s )
    {
        delete distNodes[s].grid;
        mpi::Free( distNodes[s].comm );
    }
}

} // namespace El

#endif // ifndef EL_SYMBOLIC_DISTSYMMINFO_HPP
