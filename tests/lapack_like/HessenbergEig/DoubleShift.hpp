/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SCHUR_HESSQR_DOUBLE_SHIFT_HPP
#define EL_SCHUR_HESSQR_DOUBLE_SHIFT_HPP

#include "./SingleShift.hpp"
#include "./DoubleShift/Sweep.hpp"

namespace El {
namespace schur {
namespace hess_qr {

template<typename Real,typename=EnableIf<IsReal<Real>>>
HessenbergQRInfo
DoubleShift
( Matrix<Real>& H,
  Matrix<Real>& wReal,
  Matrix<Real>& wImag,
  Matrix<Real>& Z,
  const HessenbergQRCtrl& ctrl )
{
    DEBUG_CSE
    const Real zero(0);
    const Int maxIter=30;    
    // Cf. LAPACK for these somewhat arbitrary constants
    const Real exceptScale0=Real(3)/Real(4),
               exceptScale1=Real(-4375)/Real(10000);
    const Int n = H.Height();
    const Int nZ = Z.Height();
    Int winBeg = ( ctrl.winBeg==END ? n : ctrl.winBeg );
    Int winEnd = ( ctrl.winEnd==END ? n : ctrl.winEnd );
    const Int windowSize = winEnd - winBeg;
    HessenbergQRInfo info;

    wReal.Resize( n, 1 );
    wImag.Resize( n, 1 );

    if( windowSize == 0 )
    {
        return info;
    }
    if( windowSize == 1 )
    {
        wReal(winBeg) = H(winBeg,winBeg);
        wImag(winBeg) = zero;
        return info;
    }

    // Follow LAPACK's suit and clear the two diagonals below the subdiagonal
    for( Int j=winBeg; j<winEnd-3; ++j ) 
    {
        H(j+2,j) = zero;
        H(j+3,j) = zero;
    }
    if( winBeg <= winEnd-3 )
        H(winEnd-1,winEnd-3) = zero;
    
    // Attempt to converge the eigenvalues one or two at a time
    auto ctrlSweep( ctrl );
    while( winBeg < winEnd )
    {
        Int iterBeg = winBeg;
        Int iter;
        for( iter=0; iter<maxIter; ++iter )
        {
            {
                auto winInd = IR(iterBeg,winEnd);
                iterBeg += DetectSmallSubdiagonal( H(winInd,winInd) );
            }
            if( iterBeg > winBeg )
            {
                H(iterBeg,iterBeg-1) = zero;
            }
            if( iterBeg == winEnd-1 )
            {
                wReal(iterBeg) = H(iterBeg,iterBeg);
                wImag(iterBeg) = zero;
                --winEnd;
                break;
            }
            else if( iterBeg == winEnd-2 )
            {
                Real c, s;
                lapack::TwoByTwoSchur
                ( H(winEnd-2,winEnd-2), H(winEnd-2,winEnd-1),
                  H(winEnd-1,winEnd-2), H(winEnd-1,winEnd-1),
                  c, s,
                  wReal(iterBeg),   wImag(iterBeg), 
                  wReal(iterBeg+1), wImag(iterBeg+1) );
                if( ctrl.fullTriangle )
                {
                    if( n > winEnd )
                        blas::Rot
                        ( n-winEnd,
                          &H(winEnd-2,winEnd), H.LDim(),
                          &H(winEnd-1,winEnd), H.LDim(),
                          c, s );
                    blas::Rot
                    ( winEnd-2,
                      &H(0,winEnd-2), 1,
                      &H(0,winEnd-1), 1,
                      c, s );
                }
                if( ctrl.wantSchurVecs )
                {
                    blas::Rot
                    ( nZ,
                      &Z(0,winEnd-2), 1,
                      &Z(0,winEnd-1), 1,
                      c, s );
                }
                winEnd -= 2;
                break;
            }

            // Pick either the Francis shifts or exceptional shifts
            Real eta00, eta01, eta10, eta11;
            if( iter == maxIter/3 )
            {
                const Real scale =
                  Abs(H(iterBeg+1,iterBeg)) + Abs(H(iterBeg+2,iterBeg+1));
                eta00 = exceptScale0*scale + H(iterBeg,iterBeg);
                eta01 = exceptScale1*scale;
                eta10 = scale;
                eta11 = eta00; 
            } 
            else if( iter == 2*maxIter/3 )
            {
                const Real scale = 
                  Abs(H(winEnd-1,winEnd-2)) + Abs(H(winEnd-2,winEnd-3));
                eta00 = exceptScale0*scale + H(winEnd-1,winEnd-1);
                eta01 = exceptScale1*scale;
                eta10 = scale;
                eta11 = eta00;
            }
            else
            {
                eta00 = H(winEnd-2,winEnd-2);
                eta01 = H(winEnd-2,winEnd-1);
                eta10 = H(winEnd-1,winEnd-2);
                eta11 = H(winEnd-1,winEnd-1);
            }
            Real shift0Real, shift0Imag, shift1Real, shift1Imag;
            double_shift::PrepareShifts
            ( eta00, eta01,
              eta10, eta11,
              shift0Real, shift0Imag,
              shift1Real, shift1Imag );

            ctrlSweep.winBeg = iterBeg;
            ctrlSweep.winEnd = winEnd;
            double_shift::Sweep
            ( H, shift0Real, shift0Imag, shift1Real, shift1Imag, Z, ctrlSweep );
            ++info.numIterations;
        }
        if( iter == maxIter )
        {
            if( ctrl.demandConverged )
                RuntimeError("QR iteration did not converge");
            else
                break;
        }
    }
    info.numUnconverged = winEnd-winBeg;
    return info;
}

} // namespace hess_qr
} // namespace schur
} // namespace El

#endif // ifndef EL_SCHUR_HESSQR_DOUBLE_SHIFT_HPP
