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

template<typename T>
inline void ConvertFrom1dTo2d( DistSymmFrontTree<T>& L )
{
    DEBUG_ONLY(CallStackEntry cse("ConvertFrom1dTo2d"))
    if( !FrontsAre1d(L.frontType) )
        LogicError("Treated a 2D front type as 1D");

    DistSymmFront<T>& leafFront = L.distFronts[0];
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
    const Int numDistNodes = L.distFronts.size();    
    for( Int s=1; s<numDistNodes; ++s )
    {
        DistSymmFront<T>& front = L.distFronts[s];
        front.front2dL.SetGrid( front.front1dL.Grid() );
        front.front2dL = front.front1dL;
        front.front1dL.Empty();
    }
}

template<typename T>
inline void ConvertFrom2dTo1d( DistSymmFrontTree<T>& L )
{
    DEBUG_ONLY(CallStackEntry cse("ConvertFrom2dTo1d"))
    if( FrontsAre1d(L.frontType) )
        LogicError("Treated a 1D front type as 2D");

    DistSymmFront<T>& leafFront = L.distFronts[0];
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
    const Int numDistNodes = L.distFronts.size();    
    for( Int s=1; s<numDistNodes; ++s )
    {
        DistSymmFront<T>& front = L.distFronts[s];
        front.front1dL.SetGrid( front.front2dL.Grid() );
        front.front1dL = front.front2dL;
        front.front2dL.Empty();
    }
}

// This routine could be modified later so that it uses much less memory
// by replacing the '=' redistributions with piece-by-piece redistributions.
template<typename F>
void ChangeFrontType( DistSymmFrontTree<F>& L, SymmFrontType frontType )
{
    DEBUG_ONLY(CallStackEntry cse("ChangeFrontType"))
    // Check if this call can be a no-op
    if( frontType == L.frontType ) 
        return;

    if( frontType == SYMM_1D && FrontsAre1d(L.frontType) )
    {
        // No action is required
    }
    else if( frontType == SYMM_1D && !FrontsAre1d(L.frontType) )
    {
        ConvertFrom2dTo1d( L );    
    }
    else if( frontType == SYMM_2D && !FrontsAre1d(L.frontType) ) 
    {
        // No action is required
    }
    else if( frontType == SYMM_2D && FrontsAre1d(L.frontType) )
    {
        ConvertFrom1dTo2d( L );
    }
    else if( frontType == ConvertTo2d(L.frontType) )
    {
        ConvertFrom1dTo2d( L );
    }
    else if( frontType == ConvertTo1d(L.frontType) )
    {
        ConvertFrom2dTo1d( L );
    }
    else if( SelInvFactorization(frontType) && 
             ConvertTo2d(frontType) == ConvertTo2d(AppendSelInv(L.frontType)) )
    {
        // We must perform selective inversion with a 2D distribution
        if( FrontsAre1d(L.frontType) )
            ChangeFrontType( L, ConvertTo2d(L.frontType) );
        // Perform selective inversion
        const Int numDistNodes = L.distFronts.size();    
        for( Int s=1; s<numDistNodes; ++s )
        {
            // Invert the unit-diagonal lower triangle
            DistSymmFront<F>& front = L.distFronts[s];
            const Int snSize = front.front2dL.Width();
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
