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
#ifndef EL_SYMBOLIC_SYMMNODEINFO_HPP
#define EL_SYMBOLIC_SYMMNODEINFO_HPP

namespace El {

struct DistSymmNodeInfo;

struct SymmNodeInfo
{
    // Known before analysis
    // ---------------------
    Int size, off; 
    vector<Int> origLowerStruct;

    SymmNodeInfo* parent;
    vector<SymmNodeInfo*> children;
    DistSymmNodeInfo* duplicate; 

    // Known after analysis
    // --------------------
    Int myOff;
    vector<Int> lowerStruct;
    vector<Int> origLowerRelInds;
    // (maps from the child update indices to our frontal indices).
    vector<vector<Int>> childRelInds;

    SymmNodeInfo( SymmNodeInfo* parentNode=nullptr );
    SymmNodeInfo( DistSymmNodeInfo* duplicateNode );
    ~SymmNodeInfo();
};

struct FactorCommMeta
{
    vector<int> numChildSendInds;
    // This information does not necessarily have to be kept and can be
    // computed from the above information (albeit somewhat expensively).
    mutable vector<vector<Int>> childRecvInds;

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
    vector<int> numChildSendInds;
    vector<vector<Int>> childRecvInds;

    void Empty()
    {
        SwapClear( numChildSendInds );
        SwapClear( childRecvInds );
    }
};

struct MatrixCommMeta
{
    vector<int> numChildSendInds;
    vector<vector<Int>> childRecvInds;

    void Empty()
    {
        SwapClear( numChildSendInds );
        SwapClear( childRecvInds );
    }
};

struct DistSymmNodeInfo
{
    // Known before analysis
    // ---------------------
    Int size, off;
    vector<Int> origLowerStruct;
    bool onLeft;
    mpi::Comm comm;

    DistSymmNodeInfo* parent; 
    DistSymmNodeInfo* child; 
    SymmNodeInfo* duplicate;

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

    mutable FactorCommMeta factorMeta;
    mutable MultiVecCommMeta multiVecMeta;

    DistSymmNodeInfo( DistSymmNodeInfo* parentNode=nullptr );
    ~DistSymmNodeInfo();
};

// Utilities
void ComputeFactRecvInds( const DistSymmNodeInfo& info );
void GetChildGridDims( const DistSymmNodeInfo& info, int* childGridDims );

} // namespace El

#endif // ifndef EL_SYMBOLIC_SYMMNODEINFO_HPP
