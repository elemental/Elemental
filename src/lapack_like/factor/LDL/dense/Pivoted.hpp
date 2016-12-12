/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_LDL_PIVOTED_HPP
#define EL_LDL_PIVOTED_HPP

#include "./Pivoted/BunchKaufmanA.hpp"
// TODO: Bunch-Kaufman C
#include "./Pivoted/BunchKaufmanD.hpp"
#include "./Pivoted/BunchParlett.hpp"

#include "./Pivoted/Unblocked.hpp"
#include "./Pivoted/Panel.hpp"
#include "./Pivoted/Blocked.hpp"

namespace El {
namespace ldl {

template<typename F>
inline void
Pivoted
( Matrix<F>& A,
  Matrix<F>& dSub,
  Permutation& P,
  bool conjugate,
  const LDLPivotCtrl<Base<F>>& ctrl )
{
    EL_DEBUG_CSE
    switch( ctrl.pivotType )
    {
    case BUNCH_KAUFMAN_A:
    case BUNCH_KAUFMAN_C:
    case BUNCH_KAUFMAN_D:
        pivot::Blocked( A, dSub, P, conjugate, ctrl.pivotType, ctrl.gamma );
        break;
    default:
        pivot::Unblocked( A, dSub, P, conjugate, ctrl.pivotType, ctrl.gamma );
    }
}

template<typename F>
inline void
Pivoted
( AbstractDistMatrix<F>& A,
  AbstractDistMatrix<F>& dSub,
  DistPermutation& P,
  bool conjugate,
  const LDLPivotCtrl<Base<F>>& ctrl )
{
    EL_DEBUG_CSE
    switch( ctrl.pivotType )
    {
    case BUNCH_KAUFMAN_A:
    case BUNCH_KAUFMAN_C:
    case BUNCH_KAUFMAN_D:
        pivot::Blocked( A, dSub, P, conjugate, ctrl.pivotType, ctrl.gamma );
        break;
    default:
        pivot::Unblocked( A, dSub, P, conjugate, ctrl.pivotType, ctrl.gamma );
    }
}

} // namespace ldl
} // namespace El

#endif // ifndef EL_LDL_PIVOTED_HPP
