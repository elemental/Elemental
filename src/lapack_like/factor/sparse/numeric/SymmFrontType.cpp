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
#include "El.hpp"

namespace El {

bool Unfactored( SymmFrontType type )
{ return type == SYMM_1D || type == SYMM_2D; }

bool FrontIs1D( SymmFrontType type )
{
    return type == SYMM_1D                ||
           type == LDL_1D                 ||
           type == LDL_SELINV_1D          ||
           type == LDL_INTRAPIV_1D        ||
           type == LDL_INTRAPIV_SELINV_1D ||
           type == BLOCK_LDL_1D           ||
           type == BLOCK_LDL_INTRAPIV_1D;
}

bool BlockFactorization( SymmFrontType type )
{
    return type == BLOCK_LDL_1D ||
           type == BLOCK_LDL_2D ||
           type == BLOCK_LDL_INTRAPIV_1D ||
           type == BLOCK_LDL_INTRAPIV_2D;
}

bool SelInvFactorization( SymmFrontType type )
{
    return type == LDL_SELINV_1D ||
           type == LDL_SELINV_2D ||
           type == LDL_INTRAPIV_SELINV_1D ||
           type == LDL_INTRAPIV_SELINV_2D;
}

bool PivotedFactorization( SymmFrontType type )
{
    return type == LDL_INTRAPIV_1D ||
           type == LDL_INTRAPIV_2D ||
           type == LDL_INTRAPIV_SELINV_1D ||
           type == LDL_INTRAPIV_SELINV_2D ||
           type == BLOCK_LDL_INTRAPIV_1D  ||
           type == BLOCK_LDL_INTRAPIV_2D;
}

SymmFrontType ConvertTo2D( SymmFrontType type )
{
    DEBUG_ONLY(CallStackEntry cse("ConvertTo2D"))
    SymmFrontType newType;
    switch( type )
    {
    case SYMM_1D:
    case SYMM_2D:                newType = SYMM_2D;                break;
    case LDL_1D:
    case LDL_2D:                 newType = LDL_2D;                 break;
    case LDL_SELINV_1D:
    case LDL_SELINV_2D:          newType = LDL_SELINV_2D;          break;
    case LDL_INTRAPIV_1D:
    case LDL_INTRAPIV_2D:        newType = LDL_INTRAPIV_2D;        break;
    case LDL_INTRAPIV_SELINV_1D:
    case LDL_INTRAPIV_SELINV_2D: newType = LDL_INTRAPIV_SELINV_2D; break;
    case BLOCK_LDL_1D:
    case BLOCK_LDL_2D:           newType = BLOCK_LDL_2D;           break;
    case BLOCK_LDL_INTRAPIV_1D:
    case BLOCK_LDL_INTRAPIV_2D:  newType = BLOCK_LDL_INTRAPIV_2D;  break;
    default: LogicError("Invalid front type");
    }
    return newType;
}

SymmFrontType ConvertTo1D( SymmFrontType type )
{
    DEBUG_ONLY(CallStackEntry cse("ConvertTo1D"))
    SymmFrontType newType;
    switch( type )
    {
    case SYMM_1D:
    case SYMM_2D:                newType = SYMM_1D;                break;
    case LDL_1D:
    case LDL_2D:                 newType = LDL_1D;                 break;
    case LDL_SELINV_1D:
    case LDL_SELINV_2D:          newType = LDL_SELINV_1D;          break;
    case LDL_INTRAPIV_1D:
    case LDL_INTRAPIV_2D:        newType = LDL_INTRAPIV_1D;        break;
    case LDL_INTRAPIV_SELINV_1D:
    case LDL_INTRAPIV_SELINV_2D: newType = LDL_INTRAPIV_SELINV_1D; break;
    case BLOCK_LDL_1D:
    case BLOCK_LDL_2D:           newType = BLOCK_LDL_1D;           break;
    case BLOCK_LDL_INTRAPIV_1D:
    case BLOCK_LDL_INTRAPIV_2D:  newType = BLOCK_LDL_INTRAPIV_1D;  break;
    default: LogicError("Invalid front type");
    }
    return newType;
}

SymmFrontType AppendSelInv( SymmFrontType type )
{
    DEBUG_ONLY(CallStackEntry cse("AppendSelInv"))
    SymmFrontType newType;
    switch( type )
    {
    case LDL_1D:          newType = LDL_SELINV_1D; break;
    case LDL_2D:          newType = LDL_SELINV_2D; break;
    case LDL_INTRAPIV_1D: newType = LDL_INTRAPIV_SELINV_1D; break;
    case LDL_INTRAPIV_2D: newType = LDL_INTRAPIV_SELINV_2D; break;
    default: LogicError("Sel-inv does not make sense for this type");
    }
    return newType;
}

SymmFrontType RemoveSelInv( SymmFrontType type )
{
    DEBUG_ONLY(CallStackEntry cse("RemoveSelInv"))
    SymmFrontType newType;
    switch( type )
    {
    case LDL_SELINV_1D: newType = LDL_1D; break;
    case LDL_SELINV_2D: newType = LDL_2D; break;
    case LDL_INTRAPIV_SELINV_1D: newType = LDL_INTRAPIV_1D; break;
    case LDL_INTRAPIV_SELINV_2D: newType = LDL_INTRAPIV_2D; break;
    default: LogicError("This type did not involve selective inversion");
    }
    return newType;
}

SymmFrontType InitialFactorType( SymmFrontType type )
{
    if( Unfactored(type) )
        LogicError("Front type does not require factorization");
    if( BlockFactorization(type) )
        return ConvertTo2D(type);
    else if( PivotedFactorization(type) )
        return LDL_INTRAPIV_2D;
    else
        return LDL_2D;
}

} // namespace El
