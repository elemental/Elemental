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

bool Unfactored( LDLFrontType type )
{ return type == SYMM_1D || type == SYMM_2D; }

bool FrontIs1D( LDLFrontType type )
{
    return type == SYMM_1D                ||
           type == LDL_1D                 ||
           type == LDL_SELINV_1D          ||
           type == LDL_INTRAPIV_1D        ||
           type == LDL_INTRAPIV_SELINV_1D ||
           type == BLOCK_LDL_1D           ||
           type == BLOCK_LDL_INTRAPIV_1D;
}

bool BlockFactorization( LDLFrontType type )
{
    return type == BLOCK_LDL_1D ||
           type == BLOCK_LDL_2D ||
           type == BLOCK_LDL_INTRAPIV_1D ||
           type == BLOCK_LDL_INTRAPIV_2D;
}

bool SelInvFactorization( LDLFrontType type )
{
    return type == LDL_SELINV_1D ||
           type == LDL_SELINV_2D ||
           type == LDL_INTRAPIV_SELINV_1D ||
           type == LDL_INTRAPIV_SELINV_2D;
}

bool PivotedFactorization( LDLFrontType type )
{
    return type == LDL_INTRAPIV_1D ||
           type == LDL_INTRAPIV_2D ||
           type == LDL_INTRAPIV_SELINV_1D ||
           type == LDL_INTRAPIV_SELINV_2D ||
           type == BLOCK_LDL_INTRAPIV_1D  ||
           type == BLOCK_LDL_INTRAPIV_2D;
}

LDLFrontType ConvertTo2D( LDLFrontType type )
{
    EL_DEBUG_CSE
    LDLFrontType newType;
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

LDLFrontType ConvertTo1D( LDLFrontType type )
{
    EL_DEBUG_CSE
    LDLFrontType newType;
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

LDLFrontType AppendSelInv( LDLFrontType type )
{
    EL_DEBUG_CSE
    LDLFrontType newType=LDL_SELINV_2D;
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

LDLFrontType RemoveSelInv( LDLFrontType type )
{
    EL_DEBUG_CSE
    LDLFrontType newType=LDL_2D;
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

LDLFrontType InitialFactorType( LDLFrontType type )
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
