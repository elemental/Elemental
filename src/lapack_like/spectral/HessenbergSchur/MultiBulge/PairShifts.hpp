/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SCHUR_HESS_MULTIBULGE_PAIR_SHIFTS_HPP
#define EL_SCHUR_HESS_MULTIBULGE_PAIR_SHIFTS_HPP

namespace El {
namespace hess_schur {
namespace multibulge {

// NOTE: This should only be called for real matrices, where one can assume
// conjugate pairs of shifts
template<typename Real>
void PairShifts( Matrix<Complex<Real>>& shifts )
{
    EL_DEBUG_CSE
    const Int numShifts = shifts.Height();

    Complex<Real> tmp;
    for( Int i=numShifts-1; i>=2; i-=2 )
    {
        if( shifts(i).imag() != -shifts(i-1).imag() )
        {
            tmp = shifts(i);
            shifts(i) = shifts(i-1);
            shifts(i-1) = shifts(i-2);
            shifts(i-2) = tmp;
        }
    }

    EL_DEBUG_ONLY(
      for( Int i=numShifts-1; i>=2; i-=2 )
      {
          if( shifts(i).imag() != -shifts(i-1).imag() )
          {
              Print( shifts, "shifts" );
              RuntimeError("Shifts were not properly paired");
          }
      }
    )
}

} // namespace multibulge
} // namespace hess_schur
} // namespace El

#endif // ifndef EL_SCHUR_HESS_MULTIBULGE_PAIR_SHIFTS_HPP
