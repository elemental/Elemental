/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_HESS_SCHUR_MULTIBULGE_HANDLE_TWOBYTWO_HPP
#define EL_HESS_SCHUR_MULTIBULGE_HANDLE_TWOBYTWO_HPP

#include "../Simple.hpp"

namespace El {
namespace hess_schur {
namespace multibulge {

template<typename Real>
void HandleTwoByTwo
(       Matrix<Real>& H,
        Matrix<Complex<Real>>& w,
        Matrix<Real>& Z,
        Int iterBeg,
  const HessenbergSchurCtrl& ctrl )
{
    const Int n = H.Height();
    const Int nZ = Z.Height();
    Real c, s;
    schur::TwoByTwo
    ( H(iterBeg,  iterBeg), H(iterBeg,  iterBeg+1),
      H(iterBeg+1,iterBeg), H(iterBeg+1,iterBeg+1),
      w(iterBeg), w(iterBeg+1),
      c, s );
    if( ctrl.fullTriangle )
    {
        if( n > iterBeg+2 )
            blas::Rot
            ( n-(iterBeg+2),
              &H(iterBeg,  iterBeg+2), H.LDim(),
              &H(iterBeg+1,iterBeg+2), H.LDim(),
              c, s );
        blas::Rot
        ( iterBeg, &H(0,iterBeg), 1, &H(0,iterBeg+1), 1, c, s );
    }
    if( ctrl.wantSchurVecs )
        blas::Rot( nZ, &Z(0,iterBeg), 1, &Z(0,iterBeg+1), 1, c, s );
}

template<typename Real>
void HandleTwoByTwo
(       Matrix<Complex<Real>>& H,
        Matrix<Complex<Real>>& w,
        Matrix<Complex<Real>>& Z,
        Int iterBeg,
  const HessenbergSchurCtrl& ctrl )
{
    auto ctrlSub( ctrl );
    ctrlSub.winBeg = iterBeg;
    ctrlSub.winEnd = iterBeg + 2;
    Simple( H, w, Z, ctrlSub );
}

} // namespace multibulge
} // namespace hess_schur
} // namespace El

#endif // ifndef EL_HESS_SCHUR_MULTIBULGE_HANDLE_TWOBYTWO_HPP
