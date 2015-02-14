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
#ifndef EL_FACTOR_SYMBOLIC_HPP
#define EL_FACTOR_SYMBOLIC_HPP

#include "./symbolic/Separator.hpp"
#include "./symbolic/SymmNode.hpp"
#include "./symbolic/SymmNodeInfo.hpp"
#include "./symbolic/NestedDissection.hpp"

namespace El {

Int SymmetricAnalysis
( const SymmNode& rootNode, SymmNodeInfo& rootInfo, Int myOff=0 );
void SymmetricAnalysis
( const DistSymmNode& rootNode, DistSymmNodeInfo& rootInfo,
  bool storeFactRecvInds=true );

void NaturalNestedDissection
(       Int nx,
        Int ny,
        Int nz,
  const Graph& graph,
        vector<Int>& map,
        Separator& rootSep,
        SymmNodeInfo& info,
        Int cutoff=128 );
void NaturalNestedDissection
(       Int nx,
        Int ny,
        Int nz,
  const DistGraph& graph,
        DistMap& map,
        DistSeparator& rootSep,
        DistSymmNodeInfo& info,
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
        vector<Int>& perm );

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

#endif // ifndef EL_FACTOR_SYMBOLIC_HPP
