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
#ifndef EL_SYMBOLIC_NATURALNESTEDDISSECTION_HPP
#define EL_SYMBOLIC_NATURALNESTEDDISSECTION_HPP

namespace El {

void NaturalNestedDissection
(       int nx, 
        int ny, 
        int nz,
  const DistGraph& graph, 
        DistMap& map,
        DistSeparatorTree& sepTree, 
        DistSymmInfo& info,
        int cutoff=128,
        bool storeFactRecvInds=false );

int NaturalBisect
(       int nx, 
        int ny, 
        int nz, 
  const Graph& graph, 
        int& nxLeft, 
        int& nyLeft, 
        int& nzLeft,
        Graph& leftChild, 
        int& nxRight, 
        int& nyRight, 
        int& nzRight,
        Graph& rightChild, 
        std::vector<int>& perm );

// NOTE: for two or more processes
int NaturalBisect
(       int nx, 
        int ny, 
        int nz,
  const DistGraph& graph, 
        int& nxChild, 
        int& nyChild, 
        int& nzChild,
        DistGraph& child, 
        DistMap& perm,
        bool& onLeft );

} // namespace El

#endif // ifndef EL_SYMBOLIC_NATURALNESTEDDISSECTION_HPP
