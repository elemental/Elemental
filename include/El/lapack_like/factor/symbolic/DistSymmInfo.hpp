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
#ifndef EL_SYMBOLIC_DISTSYMMINFO_HPP
#define EL_SYMBOLIC_DISTSYMMINFO_HPP

namespace El {

struct SymmNodeInfo
{
    //
    // This is known before analysis
    //
    Int size, off; 
    Int parent; // -1 if root separator
    std::vector<Int> children;
    std::vector<Int> origLowerStruct;

    //
    // The following is computed during analysis
    //
    bool onLeft;
    Int myOff;
    std::vector<Int> lowerStruct;
    std::vector<Int> origLowerRelInds;
    // (maps from the child update indices to our frontal indices).
    std::vector<Int> leftRelInds, rightRelInds;
};

struct FactorCommMeta
{
    std::vector<int> numChildSendInds;
    // This information does not necessarily have to be kept and can be
    // computed from the above information (albeit somewhat expensively).
    mutable std::vector<std::vector<Int>> childRecvInds;

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
    Int localOff, localSize;
    std::vector<int> numChildSendInds;
    std::vector<std::vector<Int>> childRecvInds;

    void Empty()
    {
        SwapClear( numChildSendInds );
        SwapClear( childRecvInds );
    }
};

struct MatrixCommMeta
{
    std::vector<int> numChildSendInds;
    std::vector<std::vector<Int>> childRecvInds;

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
    Int size, off;
    std::vector<Int> origLowerStruct;
    bool onLeft;
    mpi::Comm comm;

    //
    // The following is computed during analysis
    //
    Grid* grid;
    Int myOff, leftSize, rightSize;
    std::vector<Int> lowerStruct;
    std::vector<Int> origLowerRelInds;

    // The relative indices of our child
    // (maps from the child update indices to our frontal indices).
    // These could be replaced with just the relative indices of our local 
    // submatrices of the child updates.
    std::vector<Int> leftRelInds, rightRelInds;

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

} // namespace El

#endif // ifndef EL_SYMBOLIC_DISTSYMMINFO_HPP
