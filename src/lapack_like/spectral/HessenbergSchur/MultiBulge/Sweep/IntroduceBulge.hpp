/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SCHUR_HESS_MULTIBULGE_SWEEP_INTRODUCE_BULGE_HPP
#define EL_SCHUR_HESS_MULTIBULGE_SWEEP_INTRODUCE_BULGE_HPP

namespace El {
namespace hess_schur {
namespace multibulge {

template<typename Real>
void ImplicitQQuadraticSeed
( const Matrix<Real>& H,
  const Complex<Real>& shift0,
  const Complex<Real>& shift1,
        Real* v )
{
    EL_DEBUG_CSE
    const Real zero(0);
    const Int n = H.Height();
    EL_DEBUG_ONLY(
      if( n != 2 && n != 3 )
          LogicError("Expected n to be 2 or 3");
      const bool bothReal = ( shift0.imag() == zero && shift1.imag() == zero );
      const bool conjugate = ( shift0.imag() == -shift1.imag() );
      if( !bothReal && !conjugate )
          LogicError("Assumed shifts were either both real or conjugates");
    )
    if( n == 2 )
    {
        const Real& eta00 = H(0,0);
        const Real& eta01 = H(0,1);
        const Real& eta10 = H(1,0);
        const Real& eta11 = H(1,1);

        // It seems arbitrary whether the scale is computed relative
        // to shift0 or shift1, but we follow LAPACK's convention.
        // (While the choice is irrelevant for conjugate shifts, it is not for
        //  real shifts)
        const Real scale = OneAbs(eta00-shift1) + Abs(eta10);
        if( scale == zero )
        {
            v[0] = v[1] = zero;
        }
        else
        {
            // Normalize the first column by the scale
            Real eta10Scale = eta10 / scale;
            v[0] = eta10Scale*eta01 +
                   (eta00-shift0.real())*((eta00-shift1.real())/scale) -
                   shift0.imag()*(shift1.imag()/scale);
            v[1] = eta10Scale*(eta00+eta11-shift0.real()-shift1.real());
        }
    }
    else
    {
        const Real& eta00 = H(0,0);
        const Real& eta01 = H(0,1);
        const Real& eta02 = H(0,2);
        const Real& eta10 = H(1,0);
        const Real& eta11 = H(1,1);
        const Real& eta12 = H(1,2);
        const Real& eta20 = H(2,0);
        const Real& eta21 = H(2,1);
        const Real& eta22 = H(2,2);

        const Real scale = OneAbs(eta00-shift1) + Abs(eta10) + Abs(eta20);
        if( scale == zero )
        {
            v[0] = v[1] = v[2] = 0;
        }
        else
        {
            // Normalize the first column by the scale
            const Real eta10Scale = eta10 / scale;
            const Real eta20Scale = eta20 / scale;
            v[0] = (eta00-shift0.real())*((eta00-shift1.real())/scale) -
                   shift0.imag()*(shift1.imag()/scale) + eta01*eta10Scale +
                   eta02*eta20Scale;
            v[1] = eta10Scale*(eta00+eta11-shift0.real()-shift1.real()) +
                   eta12*eta20Scale;
            v[2] = eta20Scale*(eta00+eta22-shift0.real()-shift1.real()) +
                   eta21*eta10Scale;
        }
    }
}

template<typename Field>
void ImplicitQQuadraticSeed
( const Matrix<Field>& H,
  const Field& shift0,
  const Field& shift1,
        Field* v )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Real zero(0);
    const Int n = H.Height();
    EL_DEBUG_ONLY(
      if( n != 2 && n != 3 )
          LogicError("Expected n to be 2 or 3");
    )
    if( n == 2 )
    {
        const Field& eta00 = H(0,0);
        const Field& eta01 = H(0,1);
        const Field& eta10 = H(1,0);
        const Field& eta11 = H(1,1);

        // It seems arbitrary whether the scale is computed relative
        // to shift0 or shift1, but we follow LAPACK's convention.
        // (While the choice is irrelevant for conjugate shifts, it is not for
        //  real shifts)
        const Real scale = OneAbs(eta00-shift1) + OneAbs(eta10);
        if( scale == zero )
        {
            v[0] = v[1] = zero;
        }
        else
        {
            // Normalize the first column by the scale
            Field eta10Scale = eta10 / scale;
            v[0] = eta10Scale*eta01 + (eta00-shift0)*((eta00-shift1)/scale);
            v[1] = eta10Scale*(eta00+eta11-shift0-shift1);
        }
    }
    else
    {
        const Field& eta00 = H(0,0);
        const Field& eta01 = H(0,1);
        const Field& eta02 = H(0,2);
        const Field& eta10 = H(1,0);
        const Field& eta11 = H(1,1);
        const Field& eta12 = H(1,2);
        const Field& eta20 = H(2,0);
        const Field& eta21 = H(2,1);
        const Field& eta22 = H(2,2);

        const Real scale = OneAbs(eta00-shift1) + OneAbs(eta10) + OneAbs(eta20);
        if( scale == zero )
        {
            v[0] = v[1] = v[2] = 0;
        }
        else
        {
            // Normalize the first column by the scale
            const Field eta10Scale = eta10 / scale;
            const Field eta20Scale = eta20 / scale;
            v[0] = (eta00-shift0)*((eta00-shift1)/scale) +
                   eta01*eta10Scale + eta02*eta20Scale;
            v[1] = eta10Scale*(eta00+eta11-shift0-shift1) + eta12*eta20Scale;
            v[2] = eta20Scale*(eta00+eta22-shift0-shift1) + eta21*eta10Scale;
        }
    }
}

template<typename Field>
void IntroduceBulge
( const Matrix<Field>& H,
  const Complex<Base<Field>>& shift0,
  const Complex<Base<Field>>& shift1,
        Field* v )
{
    EL_DEBUG_CSE
    const Int n = H.Height();
    ImplicitQQuadraticSeed( H, shift0, shift1, v );
    Field beta = v[0];
    v[0] = lapack::Reflector( n, beta, &v[1], 1 );
}

} // namespace multibulge
} // namespace hess_schur
} // namespace El

#endif // ifndef EL_SCHUR_HESS_MULTIBULGE_SWEEP_INTRODUCE_BULGE_HPP
