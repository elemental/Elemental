/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_HESS_SCHUR_SIMPLE_HPP
#define EL_HESS_SCHUR_SIMPLE_HPP

#include "./Simple/SingleShift.hpp"
#include "./Simple/DoubleShift.hpp"

namespace El {
namespace hess_schur {

template<typename Real>
HessenbergSchurInfo
Simple
( Matrix<Real>& H,
  Matrix<Complex<Real>>& w,
  Matrix<Real>& Z,
  const HessenbergSchurCtrl& ctrl )
{ return DoubleShift( H, w, Z, ctrl ); }

template<typename Real>
HessenbergSchurInfo
Simple
( Matrix<Complex<Real>>& H,
  Matrix<Complex<Real>>& w,
  Matrix<Complex<Real>>& Z,
  const HessenbergSchurCtrl& ctrl )
{ return SingleShift( H, w, Z, ctrl ); }

} // namespace hess_schur
} // namespace El

#endif // ifndef EL_HESS_SCHUR_SIMPLE_HPP
