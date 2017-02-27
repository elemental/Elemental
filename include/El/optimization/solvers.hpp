/*
   Copyright (c) 2009-2017, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_OPTIMIZATION_SOLVERS_HPP
#define EL_OPTIMIZATION_SOLVERS_HPP

namespace El {

template<typename Real>
Real RelativeComplementarityGap
( const Real& primalObj, const Real& dualObj, const Real& gap )
{
    EL_DEBUG_CSE
    Real relCompGap;
    if( primalObj < Real(0) )
        relCompGap = gap / -primalObj;
    else if( dualObj > Real(0) )
        relCompGap = gap / dualObj;
    else
        relCompGap = 2; // 200% error if the signs differ inadmissibly.
    return relCompGap;
}

template<typename Real>
Real RelativeObjectiveGap( const Real& primalObj, const Real& dualObj )
{
    EL_DEBUG_CSE
    const Real relObjGap =
      Abs(primalObj-dualObj) / (Max(Abs(primalObj),Abs(dualObj))+1);
    return relObjGap;
}

} // namespace El

#include <El/optimization/solvers/LP.hpp>
#include <El/optimization/solvers/QP.hpp>
#include <El/optimization/solvers/SOCP.hpp>

#endif // ifndef EL_OPTIMIZATION_SOLVERS_HPP
