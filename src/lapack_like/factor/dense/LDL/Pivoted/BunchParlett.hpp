/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_LDL_PIVOTED_BUNCHPARLETT_HPP
#define EL_LDL_PIVOTED_BUNCHPARLETT_HPP

namespace El {
namespace ldl {
namespace pivot {

template<typename F>
inline LDLPivot
BunchParlett( const Matrix<F>& A, Base<F> gamma )
{
    DEBUG_ONLY(CallStackEntry cse("ldl::pivot::BunchParlett"))
    typedef Base<F> Real;
    if( gamma == Real(0) )
        gamma = LDLPivotConstant<Real>( BUNCH_PARLETT );

    const ValueInt<Real> diagMax = VectorMaxAbs( GetDiagonal(A) );
    const ValueIntPair<Real> offDiagMax = SymmetricMaxAbs( LOWER, A );

    LDLPivot pivot;
    if( diagMax.value >= gamma*offDiagMax.value )
    {
        pivot.nb = 1;
        pivot.from[0] = diagMax.index;
        return pivot;
    }
    else
    {
        pivot.nb = 2;
        pivot.from[0] = offDiagMax.indices[0];
        pivot.from[1] = offDiagMax.indices[1];
        return pivot;
    }
}

template<typename F>
inline LDLPivot
BunchParlett( const DistMatrix<F>& A, Base<F> gamma )
{
    DEBUG_ONLY(CallStackEntry cse("ldl::pivot::BunchParlett"))
    typedef Base<F> Real;
    if( gamma == Real(0) )
        gamma = LDLPivotConstant<Real>( BUNCH_PARLETT );

    const ValueInt<Real> diagMax = VectorMaxAbs( GetDiagonal(A) );
    const ValueIntPair<Real> offDiagMax = SymmetricMaxAbs( LOWER, A );

    LDLPivot pivot;
    if( diagMax.value >= gamma*offDiagMax.value )
    {
        pivot.nb = 1;
        pivot.from[0] = diagMax.index;
        return pivot;
    }
    else
    {
        pivot.nb = 2; 
        pivot.from[0] = offDiagMax.indices[0];
        pivot.from[1] = offDiagMax.indices[1];
        return pivot;
    }
}

} // namespace pivot
} // namespace ldl
} // namespace El

#endif // ifndef EL_LDL_PIVOTED_BUNCHPARLETT_HPP
