/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_HESS_SCHUR_UTIL_MAKE_SUBDIAGONAL_REAL_HPP
#define EL_HESS_SCHUR_UTIL_MAKE_SUBDIAGONAL_REAL_HPP

#include "./Gather.hpp"

namespace El {
namespace hess_schur {
namespace util {

// Consider applying a unitary diagonal similarity transformation to make
// the subdiagonal entries of the upper-left block of the following complex
// matrix real:
//
//       ~ ~ ~ ~ ~ ~
//    -------------------------------
//   | x x x x x x x | x x x x x x x |
// ~ | x x x x x x x | x x x x x x x |
// ~ |   x x x x x x | x x x x x x x |
// ~ |     x x x x x | x x x x x x x |
// ~ |       x x x x | x x x x x x x |
// ~ |         x x x | x x x x x x x |
// ~ |           x x | x x x x x x x |
//   |---------------|---------------|
//   |             x | x x x x x x x |
//   |               | x x x x x x x |
//   |               |   x x x x x x |
//   |               |     x x x x x |
//   |               |       x x x x |
//   |               |         x x x |
//   |               |           x x |
//    -------------------------------
//
// As can be easily deduced from the diagram, the maximum row index that
// need be effected is the first past the upper-left window.

template<typename Real>
void MakeSubdiagonalReal
( Matrix<Complex<Real>>& H,
  Matrix<Complex<Real>>& Z,
  const HessenbergSchurCtrl& ctrl )
{
    EL_DEBUG_CSE
    const Real zero(0);
    const Int n = H.Height();
    const Int nZ = Z.Height();
    const Int winBeg = ( ctrl.winBeg==END ? n : ctrl.winBeg );
    const Int winEnd = ( ctrl.winEnd==END ? n : ctrl.winEnd );
    const Int scaleBeg = ( ctrl.fullTriangle ? 0 : winBeg );
    const Int scaleEnd = ( ctrl.fullTriangle ? n : winEnd );
    // TODO(poulson): Separate the phase formation and similarity application?
    for( Int i=winBeg+1; i<winEnd; ++i )
    {
        Complex<Real>& eta = H(i,i-1);
        if( ImagPart(eta) != zero )
        {
            Complex<Real> phase = eta / OneAbs(eta);
            phase = Conj(phase) / Abs(phase);
            eta = Abs(eta);
            blas::Scal( scaleEnd-i, phase, &H(i,i), H.LDim() );
            blas::Scal
            ( Min(scaleEnd,i+2)-scaleBeg, Conj(phase), &H(scaleBeg,i), 1 );
            if( ctrl.wantSchurVecs )
                blas::Scal( nZ, Conj(phase), &Z(0,i), 1 );
        }
    }
}

template<typename Real>
void MakeSubdiagonalReal
( DistMatrix<Complex<Real>,MC,MR,BLOCK>& H,
  DistMatrix<Complex<Real>,MC,MR,BLOCK>& Z,
  const HessenbergSchurCtrl& ctrl )
{
    EL_DEBUG_CSE
    const Real zero(0);
    const Int n = H.Height();
    const Int winBeg = ( ctrl.winBeg==END ? n : ctrl.winBeg );
    const Int winEnd = ( ctrl.winEnd==END ? n : ctrl.winEnd );
    const Int winSize = winEnd - winBeg;
    const auto winInd = IR(winBeg,winEnd);
    const Int scaleBeg = ( ctrl.fullTriangle ? 0 : winBeg );
    const Int scaleEnd = ( ctrl.fullTriangle ? n : winEnd );

    // Compute the phase of the unitary diagonal matrix whose similarity
    // transformations makes H have a real subdiagonal
    DistMatrix<Complex<Real>,STAR,STAR> hSubWin(H.Grid());
    util::GatherSubdiagonal( H, winInd, hSubWin );
    DistMatrix<Complex<Real>,STAR,STAR> phase(H.Grid());
    Zeros( phase, winSize-1, 1 );
    for( Int i=0; i<winSize-1; ++i )
    {
        Complex<Real> eta = hSubWin.GetLocal(i,0);
        if( ImagPart(eta) == zero )
        {
            phase.SetLocal(i,0,Real(1));
        }
        else
        {
            eta /= OneAbs(eta);
            phase.SetLocal(i,0,Conj(eta)/Abs(eta));
        }
    }

    // TODO(poulson): Provide support for DiagonalScaleTrapezoid for BlockMatrix
    auto HLeft = H( IR(winBeg+1,winEnd), IR(winBeg,scaleEnd) );
    auto HRight = H( IR(scaleBeg,Min(scaleEnd,winEnd+1)), IR(winBeg+1,winEnd) );
    DiagonalScale( LEFT, NORMAL, phase, HLeft );
    DiagonalScale( RIGHT, ADJOINT, phase, HRight );
    if( ctrl.wantSchurVecs )
    {
        auto ZWin = Z( ALL, IR(winBeg+1,winEnd) );
        DiagonalScale( RIGHT, ADJOINT, phase, ZWin );
    }
}

} // namespace util
} // namespace hess_schur
} // namespace El

#endif // ifndef EL_HESS_SCHUR_UTIL_MAKE_SUBDIAGONAL_REAL_HPP
