/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_LDL_PIVOTED_BUNCHPARLETT_HPP
#define EL_LDL_PIVOTED_BUNCHPARLETT_HPP

namespace El {
namespace ldl {
namespace pivot {

template<typename F>
inline LDLPivot
BunchParlett( const Matrix<F>& A, Base<F> gamma )
{
    EL_DEBUG_CSE
    typedef Base<F> Real;
    if( gamma == Real(0) )
        gamma = LDLPivotConstant<Real>( BUNCH_PARLETT );

    const ValueInt<Real> diagMax = VectorMaxAbsLoc( GetDiagonal(A) );
    const Entry<Real> offDiagMax = SymmetricMaxAbsLoc( LOWER, A );

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
        pivot.from[0] = offDiagMax.i;
        pivot.from[1] = offDiagMax.j;
        return pivot;
    }
}

template<typename F>
inline LDLPivot
BunchParlett( const DistMatrix<F>& A, Base<F> gamma )
{
    EL_DEBUG_CSE
    typedef Base<F> Real;
    if( gamma == Real(0) )
        gamma = LDLPivotConstant<Real>( BUNCH_PARLETT );

    const ValueInt<Real> diagMax = VectorMaxAbsLoc( GetDiagonal(A) );
    const Entry<Real> offDiagMax = SymmetricMaxAbsLoc( LOWER, A );

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
        pivot.from[0] = offDiagMax.i;
        pivot.from[1] = offDiagMax.j;
        return pivot;
    }
}

} // namespace pivot
} // namespace ldl
} // namespace El

#endif // ifndef EL_LDL_PIVOTED_BUNCHPARLETT_HPP
