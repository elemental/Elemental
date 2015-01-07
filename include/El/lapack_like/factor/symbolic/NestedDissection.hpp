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
#ifndef EL_SYMBOLIC_NESTEDDISSECTION_HPP
#define EL_SYMBOLIC_NESTEDDISSECTION_HPP

#ifdef EL_HAVE_PARMETIS

#include "parmetis.h"
extern "C" {
void ElBisect
( idx_t* nvtxs, idx_t* xAdj, idx_t* adjacency, idx_t* nseps, real_t* imbalance,
  idx_t* perm, idx_t* sizes );
void ElParallelBisect
( idx_t* vtxDist, idx_t* xAdj, idx_t* adjacency, 
  idx_t* nparseps, idx_t* nseqseps, real_t* imbalance, idx_t* options, 
  idx_t* perm, idx_t* sizes, MPI_Comm* comm );
} // extern "C"

#elif defined(EL_HAVE_METIS)

#include "metis.h"
extern "C" {
void ElBisect
( idx_t* nvtxs, idx_t* xAdj, idx_t* adjacency, idx_t* nseps, real_t* imbalance,
  idx_t* perm, idx_t* sizes );
} // extern "C"

#endif

namespace El {

struct BisectCtrl
{
    bool sequential;
    Int numDistSeps;
    Int numSeqSeps;
    Int cutoff;
    bool storeFactRecvInds;

    BisectCtrl()
    : sequential(true), numDistSeps(1), numSeqSeps(1), cutoff(128),
      storeFactRecvInds(false) 
    { }
};

void NestedDissection
( const DistGraph& graph, 
        DistMap& map,
        DistSeparatorTree& sepTree, 
        DistSymmInfo& info,
  const BisectCtrl& ctrl=BisectCtrl() );

Int Bisect
( const Graph& graph, 
        Graph& leftChild, 
        Graph& rightChild, 
        std::vector<Int>& perm, 
  const BisectCtrl& ctrl=BisectCtrl() );

// NOTE: for two or more processes
Int Bisect
( const DistGraph& graph, 
        DistGraph& child, 
        DistMap& perm,
        bool& onLeft,
  const BisectCtrl& ctrl=BisectCtrl() );

int DistributedDepth( mpi::Comm comm );
void EnsurePermutation( const std::vector<Int>& map );
void EnsurePermutation( const DistMap& map );
void ReverseOrder( DistSeparatorTree& sepTree, DistSymmElimTree& eTree );

void BuildChildrenFromPerm
( const Graph& graph, const std::vector<Int>& perm, 
  Int leftChildSize, Graph& leftChild,
  Int rightChildSize, Graph& rightChild );
void BuildChildFromPerm
( const DistGraph& graph, const DistMap& perm,
  Int leftChildSize, Int rightChildSize,
  bool& onLeft, DistGraph& child );

void BuildMap
( const DistGraph& graph, 
  const DistSeparatorTree& sepTree, 
        DistMap& map );

} // namespace El

#endif // ifndef EL_SYMBOLIC_NESTEDDISSECTION_HPP
