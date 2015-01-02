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
#ifndef EL_SYMBOLIC_SYMMETRICANALYSIS_HPP
#define EL_SYMBOLIC_SYMMETRICANALYSIS_HPP

namespace El {

void SymmetricAnalysis
( const DistSymmElimTree& eTree, DistSymmInfo& info, 
  bool storeFactRecvInds=true );

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

void LocalSymmetricAnalysis
( const DistSymmElimTree& eTree, DistSymmInfo& info );
void DistSymmetricAnalysis
( const DistSymmElimTree& eTree, DistSymmInfo& info, 
  bool storeFactRecvInds=true );

inline void SymmetricAnalysis
( const DistSymmElimTree& eTree, DistSymmInfo& info, bool storeFactRecvInds )
{
    DEBUG_ONLY(CallStackEntry cse("SymmetricAnalysis"))
    LocalSymmetricAnalysis( eTree, info );
    DistSymmetricAnalysis( eTree, info, storeFactRecvInds );
}

} // namespace El

#endif // ifndef EL_SYMBOLIC_SYMMETRICANALYSIS_HPP
