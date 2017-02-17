/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_HESS_SCHUR_MULTIBULGE_TWOBYTWO_HPP
#define EL_HESS_SCHUR_MULTIBULGE_TWOBYTWO_HPP

#include "../Simple.hpp"

namespace El {
namespace hess_schur {
namespace multibulge {

template<typename Real>
void TwoByTwo
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
        {
            // Apply the Givens rotation from the left to the region right of
            // the 2x2 diagonal block
            blas::Rot
            ( n-(iterBeg+2),
              &H(iterBeg,  iterBeg+2), H.LDim(),
              &H(iterBeg+1,iterBeg+2), H.LDim(),
              c, s );
        }
        // Apply the Givens rotation from the right to the region above the
        // 2x2 diagonal block
        blas::Rot
        ( iterBeg, &H(0,iterBeg), 1, &H(0,iterBeg+1), 1, c, s );
    }
    if( ctrl.wantSchurVecs )
    {
        // Apply the Givens rotation from the right
        blas::Rot( nZ, &Z(0,iterBeg), 1, &Z(0,iterBeg+1), 1, c, s );
    }
}

template<typename Real>
void TwoByTwo
(       Matrix<Complex<Real>>& H,
        Matrix<Complex<Real>>& w,
        Matrix<Complex<Real>>& Z,
        Int offset,
  const HessenbergSchurCtrl& ctrl )
{
    EL_DEBUG_CSE
    const Int n = H.Height();

    // Compute the 2x2 Schur decomposition HSub = ZSub TSub ZSub',
    // where HSub is overwritten by TSub, and wSub by diag(TSub).
    //
    // It might be slight overkill to call HessenbergSchur rather than a
    // specialized routine that returns the (c,s) pair of a Givens rotation
    // to define the Schur vectors, but there is not currently complex support
    // for such a routine. However, it is easy to observe that, given any 2x2
    // complex unitary matrix of Schur vectors, the fact that the columns can
    // have arbitrary phase implies that the top-left and bottom-right entries
    // can be rescaled to be real (and both equal to 'c').
    auto HSub = H( IR(offset,offset+2), IR(offset,offset+2) );
    auto wSub = w( IR(offset,offset+2), ALL );
    Matrix<Complex<Real>> ZSub;
    HessenbergSchur( HSub, wSub, ZSub );

    if( ctrl.fullTriangle )
    {
        if( n > offset+2 )
        {
            // Overwrite H((offset,offset+1),offset+2:end) *= ZSub'
            // (applied from the left)
            auto HRight = H( ALL, IR(offset+2,END) );
            Matrix<Complex<Real>> ZSubAdj;
            Adjoint( ZSub, ZSubAdj );
            Transform2x2Rows( ZSubAdj, HRight, offset, offset+1 );
        }
        // Overwrite H(0:offset,(offset,offset+1)) *= ZSub
        auto HTop = H( IR(0,offset), ALL );
        Transform2x2Cols( ZSub, HTop, offset, offset+1 );
    }
    if( ctrl.wantSchurVecs )
    {
        // Overwrite Z(:,(offset,offset+1)) *= ZSub
        Transform2x2Cols( ZSub, Z, offset, offset+1 );
    }
}

template<typename Field>
void TwoByTwo
(       DistMatrix<Field,MC,MR,BLOCK>& H,
        Field eta00,
        Field eta01,
        Field eta10,
        Field eta11,
        DistMatrix<Complex<Base<Field>>,STAR,STAR>& w,
        DistMatrix<Field,MC,MR,BLOCK>& Z,
        Int offset,
  const HessenbergSchurCtrl& ctrl )
{
    EL_DEBUG_CSE
    const Int n = H.Height();
    auto& wLoc = w.Matrix();

    // Compute the 2x2 Schur decomposition HSub = ZSub TSub ZSub',
    // where HSub is overwritten by TSub, and wSub by diag(TSub).
    //
    // It might be slight overkill to call HessenbergSchur rather than a
    // specialized routine that returns the (c,s) pair of a Givens rotation
    // to define the Schur vectors, but there is not currently complex support
    // for such a routine. However, it is easy to observe that, given any 2x2
    // complex unitary matrix of Schur vectors, the fact that the columns can
    // have arbitrary phase implies that the top-left and bottom-right entries
    // can be rescaled to be real (and both equal to 'c').
    Matrix<Field> HSub(2,2);
    HSub(0,0) = eta00;
    HSub(0,1) = eta01;
    HSub(1,0) = eta10;
    HSub(1,1) = eta11;
    auto wSub = wLoc( IR(offset,offset+2), ALL );
    Matrix<Field> ZSub;
    HessenbergSchur( HSub, wSub, ZSub );

    H.Set( offset,   offset,   HSub(0,0) );
    H.Set( offset,   offset+1, HSub(0,1) );
    H.Set( offset+1, offset,   HSub(1,0) );
    H.Set( offset+1, offset+1, HSub(1,1) );

    if( ctrl.fullTriangle )
    {
        if( n > offset+2 )
        {
            // Overwrite H((offset,offset+1),offset+2:end) *= ZSub'
            // (applied from the left)
            auto HRight = H( ALL, IR(offset+2,END) );
            Matrix<Field> ZSubAdj;
            Adjoint( ZSub, ZSubAdj );
            Transform2x2Rows( ZSubAdj, HRight, offset, offset+1 );
        }
        // Overwrite H(0:offset,(offset,offset+1)) *= ZSub
        auto HTop = H( IR(0,offset), ALL );
        Transform2x2Cols( ZSub, HTop, offset, offset+1 );
    }
    if( ctrl.wantSchurVecs )
    {
        // Overwrite Z(:,(offset,offset+1)) *= ZSub
        Transform2x2Cols( ZSub, Z, offset, offset+1 );
    }
}

} // namespace multibulge
} // namespace hess_schur
} // namespace El

#endif // ifndef EL_HESS_SCHUR_MULTIBULGE_TWOBYTWO_HPP
