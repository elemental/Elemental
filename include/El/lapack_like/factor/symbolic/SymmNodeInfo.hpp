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

    DistSymmNodeInfo( DistSymmNodeInfo* parentNode=nullptr );
    ~DistSymmNodeInfo();
};

void GetChildGridDims
( const DistSymmNodeInfo& info, 
  vector<int>& gridHeights, vector<int>& gridWidths );

} // namespace El

#endif // ifndef EL_SYMBOLIC_SYMMNODEINFO_HPP
