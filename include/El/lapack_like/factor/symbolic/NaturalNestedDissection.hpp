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
#ifndef EL_SYMBOLIC_NATURALNESTEDDISSECTION_HPP
#define EL_SYMBOLIC_NATURALNESTEDDISSECTION_HPP

namespace El {

void NaturalNestedDissection
(       Int nx, 
        Int ny, 
        Int nz,
  const DistGraph& graph, 
        DistMap& map,
        DistSeparatorTree& sepTree, 
        DistSymmInfo& info,
        Int cutoff=128,
        bool storeFactRecvInds=false );

Int NaturalBisect
(       Int nx, 
        Int ny, 
        Int nz, 
  const Graph& graph, 
        Int& nxLeft, 
        Int& nyLeft, 
        Int& nzLeft,
        Graph& leftChild, 
        Int& nxRight, 
        Int& nyRight, 
        Int& nzRight,
        Graph& rightChild, 
        std::vector<Int>& perm );

// NOTE: for two or more processes
Int NaturalBisect
(       Int nx, 
        Int ny, 
        Int nz,
  const DistGraph& graph, 
        Int& nxChild, 
        Int& nyChild, 
        Int& nzChild,
        DistGraph& child, 
        DistMap& perm,
        bool& onLeft );

} // namespace El

#endif // ifndef EL_SYMBOLIC_NATURALNESTEDDISSECTION_HPP
