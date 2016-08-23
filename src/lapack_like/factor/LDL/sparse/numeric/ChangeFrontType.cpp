/*
   Copyright (c) 2009-2012, Jack Poulson, Lexing Ying, and 
   The University of Texas at Austin.
   All rights reserved.

   Copyright (c) 2013, Jack Poulson, Lexing Ying, and Stanford University.
   All rights reserved.

   Copyright (c) 2013-2014, Jack Poulson and 
   The Georgia Institute of Technology.
   All rights reserved.

   Copyright (c) 2014-2015, Jack Poulson and Stanford University.
   All rights reserved.
   
   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {
namespace ldl {

template<typename F>
void ChangeFrontType( Front<F>& front, LDLFrontType type, bool recurse )
{
    DEBUG_CSE

    if( type == SYMM_1D || type == SYMM_2D || 
        type == ConvertTo1D(front.type) || type == ConvertTo2D(front.type) )
    {
        // No-op
    }
    else
        LogicError("Unavailable front type change");

    front.type = type;
    if( recurse )
        for( auto* child : front.children )
            ChangeFrontType( *child, type, recurse );
}

template<typename F>
void ChangeFrontType( DistFront<F>& front, LDLFrontType type, bool recurse )
{
    DEBUG_CSE

    if( type == SYMM_1D || type == ConvertTo1D(front.type) )
    {
        if( !FrontIs1D(front.type) )
        {
            if( front.duplicate != nullptr )
            {
                front.L1D.Attach( front.L2D.Grid(), front.duplicate->LDense ); 
            }
            else
            {
                front.L1D.SetGrid( front.L2D.Grid() );
                front.L1D = front.L2D;
            }
            front.L2D.Empty();
        }
    }
    else if( type == SYMM_2D || type == ConvertTo2D(front.type) ) 
    {
        if( FrontIs1D(front.type) )
        {
            if( front.duplicate != nullptr ) 
            {
                front.L2D.Attach( front.L1D.Grid(), front.duplicate->LDense );
            }
            else
            {
                front.L2D.SetGrid( front.L1D.Grid() );
                front.L2D = front.L1D;
            }
            front.L1D.Empty();
        }
    }
    else if( SelInvFactorization(type) && 
             ConvertTo2D(type) == ConvertTo2D(AppendSelInv(front.type)) )
    {
        // Switch to 2D as soon as possible
        if( FrontIs1D(front.type) && !FrontIs1D(type) )
        {
            if( front.duplicate != nullptr )
            {
                front.L2D.Attach( front.L1D.Grid(), front.duplicate->LDense );
            }
            else
            {
                front.L2D.SetGrid( front.L1D.Grid() );
                front.L2D = front.L1D;
            }
            front.L1D.Empty();
        }

        // Invert the unit-diagonal lower triangle if it is distributed
        if( front.child != nullptr )
        {
            if( FrontIs1D(front.type) && FrontIs1D(type) )
            {
                const Int snSize = front.L1D.Width();
                auto LT = front.L1D( IR(0,snSize), IR(0,snSize) );
                TriangularInverse( LOWER, UNIT, LT );
                FillDiagonal( LT, F(1) );
                MakeTrapezoidal( LOWER, LT );
            }
            else
            {
                const Int snSize = front.L2D.Width();
                auto LT = front.L2D( IR(0,snSize), IR(0,snSize) );
                TriangularInverse( LOWER, UNIT, LT );
                FillDiagonal( LT, F(1) );
                MakeTrapezoidal( LOWER, LT );
            }
        }

        // Switch to 1D as late as possible
        if( !FrontIs1D(front.type) && FrontIs1D(type) )
        {
            if( front.duplicate != nullptr )
            {
                front.L1D.Attach( front.L2D.Grid(), front.duplicate->LDense );
            }
            else
            {
                front.L1D.SetGrid( front.L2D.Grid() );
                front.L1D = front.L2D;
            }
            front.L2D.Empty();
        }
    }
    else
        LogicError("Unavailable front type change");

    if( front.child == nullptr )
    {
        if( SelInvFactorization(type) )
            front.type = RemoveSelInv(type);
        else
            front.type = type;
    }
    else
    {
        front.type = type;
        if( recurse )
            ChangeFrontType( *front.child, type, recurse );
    }
}

#define PROTO(F) \
  template void ChangeFrontType \
  ( Front<F>& front, LDLFrontType type, bool recurse ); \
  template void ChangeFrontType \
  ( DistFront<F>& front, LDLFrontType type, bool recurse );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace ldl
} // namespace El
