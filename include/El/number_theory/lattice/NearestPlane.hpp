/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_LATTICE_NEAREST_PLANE_HPP
#define EL_LATTICE_NEAREST_PLANE_HPP

namespace El {

// TODO: Add correctness tests for this routine
template<typename F>
void NearestPlane
( const Matrix<F>& B,
  const Matrix<F>& QR,
  const Matrix<F>& t,
  const Matrix<Base<F>>& d,
  const Matrix<F>& T,
        Matrix<F>& Y,
  const LLLCtrl<Base<F>>& ctrl )
{
    EL_DEBUG_CSE
    const Int m = B.Height();
    const Int n = B.Width();
    const Int numRHS = T.Width();

    // Compute the components of T in the directions of the columns of Q
    Y = T;
    Matrix<F> RY( Y );
    qr::ApplyQ( LEFT, ADJOINT, QR, t, d, RY );

    // Run Babai's nearest plane algorithm on each column of Y.
    // The result is an approximation of the difference between each column
    // of Y and the nearest lattice point.
    // TODO(poulson): Extract this into a separate routine shared with LLL
    vector<F> xBuf(n);
    for( Int j=0; j<numRHS; ++j )
    {
        F* yBuf = &Y(0,j);
        F* rYBuf = &RY(0,j);

        Int numNonzero=0;
        for( Int i=n-1; i>=0; --i )
        {
            F chi = rYBuf[i] / QR(i,i);
            if( Abs(RealPart(chi)) > ctrl.eta ||
                Abs(ImagPart(chi)) > ctrl.eta )
            {
                chi = Round(chi);
                blas::Axpy
                ( i+1, -chi,
                  &QR(0,i), 1,
                  rYBuf,    1 );
                ++numNonzero;
            }
            else
                chi = 0;
            xBuf[i] = chi;
        }
        const float nonzeroRatio = float(numNonzero)/float(n);
        if( nonzeroRatio >= ctrl.blockingThresh )
        {
            blas::Gemv
            ( 'N', m, n,
              F(-1), B.Buffer(), B.LDim(),
                     &xBuf[0],   1,
              F(+1), yBuf,       1 );
        }
        else
        {
            for( Int i=n-1; i>=0; --i )
            {
                const F chi = xBuf[i];
                if( chi == F(0) )
                    continue;
                blas::Axpy
                ( m, -chi,
                  &B(0,i), 1,
                  yBuf,    1 );
            }
        }
    }

    // Y := T - Y
    Y *= -1;
    Y += T;
}

template<typename F>
void NearestPlane
( const Matrix<F>& B,
  const Matrix<F>& T,
        Matrix<F>& Y,
  const LLLCtrl<Base<F>>& ctrl )
{
    EL_DEBUG_CSE

    // LLL-reduce B
    Matrix<F> BRed( B ), QR, t;
    Matrix<Base<F>> d;
    auto info = LLLWithQ( BRed, QR, t, d, ctrl );

    auto BRedLeft = BRed( ALL, IR(0,info.rank) );
    auto QRLeft = QR( ALL, IR(0,info.rank) );
    auto tLeft = t( IR(0,info.rank), ALL );
    auto dLeft = d( IR(0,info.rank), ALL );
    NearestPlane( BRedLeft, QRLeft, tLeft, dLeft, T, Y, ctrl );
}

} // namespace El

#endif // ifndef EL_LATTICE_NEAREST_PLANE_HPP
