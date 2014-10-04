/*
   Copyright (c) 2009-2014, Jack Poulson, Lexing Ying,
   The University of Texas at Austin, Stanford University, and the
   Georgia Insitute of Technology.
   All rights reserved.
 
   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// This routine could be modified later so that it uses much less memory
// by replacing the '=' redistributions with piece-by-piece redistributions.
template<typename F>
void ChangeFrontType( DistSymmFrontTree<F>& L, SymmFrontType frontType )
{
    DEBUG_ONLY(CallStackEntry cse("ChangeFrontType"))
    // Check if this call can be a no-op
    if( frontType == L.frontType ) 
        return;
    const int numDistNodes = L.distFronts.size();    
    const SymmFrontType oldFrontType = L.frontType;
    DistSymmFront<F>& leafFront = L.distFronts[0];

    if( frontType == ConvertTo2d(oldFrontType) )
    {
        // 1d -> 2d
        if( leafFront.front1dL.Locked() )
            leafFront.front2dL.LockedAttach
            ( leafFront.front1dL.LocalHeight(), 
              leafFront.front1dL.LocalWidth(), 
              leafFront.front1dL.Grid(), 0, 0,
              leafFront.front1dL.LockedMatrix() );
        else
            leafFront.front2dL.Attach
            ( leafFront.front1dL.Height(), 
              leafFront.front1dL.Width(), 
              leafFront.front1dL.Grid(), 0, 0, 
              leafFront.front1dL.Matrix() );
        for( int s=1; s<numDistNodes; ++s )
        {
            DistSymmFront<F>& front = L.distFronts[s];
            front.front2dL.SetGrid( front.front1dL.Grid() );
            front.front2dL = front.front1dL;
            front.front1dL.Empty();
        }
    }
    else if( frontType == ConvertTo1d(oldFrontType) )
    {
        // 2d -> 1d
        if( leafFront.front2dL.Locked() )
            leafFront.front1dL.LockedAttach
            ( leafFront.front2dL.Height(),
              leafFront.front2dL.Width(), 
              leafFront.front2dL.Grid(), 0, 0, 
              leafFront.front2dL.LockedMatrix() );
        else
            leafFront.front1dL.Attach
            ( leafFront.front2dL.Height(), 
              leafFront.front2dL.Width(), 
              leafFront.front2dL.Grid(), 0, 0, 
              leafFront.front2dL.Matrix() );
        for( int s=1; s<numDistNodes; ++s )
        {
            DistSymmFront<F>& front = L.distFronts[s];
            front.front1dL.SetGrid( front.front2dL.Grid() );
            front.front1dL = front.front2dL;
            front.front2dL.Empty();
        }
    }
    else if( SelInvFactorization(frontType) && 
             ConvertTo2d(frontType) == ConvertTo2d(AppendSelInv(oldFrontType)) )
    {
        // We must perform selective inversion with a 2D distribution
        if( FrontsAre1d(oldFrontType) )
            ChangeFrontType( L, ConvertTo2d(oldFrontType) );
        // Perform selective inversion
        for( int s=1; s<numDistNodes; ++s )
        {
            // Invert the unit-diagonal lower triangle
            DistSymmFront<F>& front = L.distFronts[s];
            const int snSize = front.front2dL.Width();
            auto LT = front.front2dL( IR(0,snSize), IR(0,snSize) );
            TriangularInverse( LOWER, UNIT, LT );
        }
        // Convert to 1D if necessary
        if( FrontsAre1d(frontType) )
        {
            L.frontType = ConvertTo2d(frontType);
            ChangeFrontType( L, frontType );
        }
    }
    else
        LogicError("Unavailable front type change");
    L.frontType = frontType;
}

#define PROTO(F) \
  template void ChangeFrontType \
  ( DistSymmFrontTree<F>& L, SymmFrontType frontType );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
