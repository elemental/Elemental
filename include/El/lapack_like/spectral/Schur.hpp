/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SPECTRAL_SCHUR_HPP
#define EL_SPECTRAL_SCHUR_HPP

namespace El {
namespace schur {

// Put a two-by-two nonsymmetric real matrix into standard form
// ============================================================
// Compute the Schur factorization of a real 2x2 nonsymmetric matrix A
// in a manner similar to xLANV2, returning the cosine and sine terms as well
// as the real and imaginary parts of the two eigenvalues.
//
// Either A is overwritten with its real Schur factor (if it exists), or 
// it is put into the form 
//
//   | alpha00, alpha01 | = | c -s | | beta00 beta01 | |  c s |,
//   | alpha10, alpha11 |   | s  c | | beta10 beta11 | | -s c |
//
// where beta00 = beta11 and beta10*beta01 < 0, so that the two eigenvalues 
// are beta00 +- sqrt(beta10*beta01).
//

template<typename Real,typename>
void TwoByTwo
( Real& alpha00, Real& alpha01,
  Real& alpha10, Real& alpha11,
  Complex<Real>& lambda0,
  Complex<Real>& lambda1,
  Real& c, Real& s )
{
    const Real zero(0), one(1);
    const Real multiple(4);
    const Real epsilon = limits::Epsilon<Real>();

    if( alpha10 == zero )
    {
        c = one;
        s = zero;
    }
    else if( alpha01 == zero )
    {
        c = zero;
        s = one;
        Real tmp = alpha11;
        alpha11 = alpha00;
        alpha00 = tmp;
        alpha01 = -alpha10;
        alpha10 = zero;
    }
    else if( (alpha00-alpha11) == zero && Sgn(alpha01) != Sgn(alpha10) )
    {
        c = one;
        s = zero;
    }
    else
    {
        Real tmp = alpha00-alpha11;
        Real p = tmp/2;
        Real offDiagMax = Max( Abs(alpha01), Abs(alpha10) );
        Real offDiagMin = Min( Abs(alpha01), Abs(alpha10) );
        Real offDiagMinSigned = offDiagMin*Sgn(alpha01)*Sgn(alpha10);
        Real scale = Max( Abs(p), offDiagMax );
        Real z = (p/scale)*p + (offDiagMax/scale)*offDiagMinSigned;
        if( z >= multiple*epsilon )
        {
            // Compute the real eigenvalues
            z = p + Sqrt(scale)*Sqrt(z)*Sgn(p);
            alpha00 = alpha11 + z;
            alpha11 -= (offDiagMax/z)*offDiagMinSigned;

            // Compute the rotation matrix
            Real tau = SafeNorm( alpha10, z );
            c = z / tau;
            s = alpha10 / tau;
            alpha01 -= alpha10;
            alpha10 = zero;
        }
        else
        {
            // We have complex or (almost) equal real eigenvalues, so force
            // alpha00 and alpha11 to be equal 
            Real sigma = alpha01 + alpha10;
            Real tau = SafeNorm( sigma, tmp );
            c = Sqrt( (one + Abs(sigma)/tau)/2 );
            s = -(p/(tau*c))*Sgn(sigma);

            // B := A [c, -s; s, c]
            Real beta00 =  alpha00*c + alpha01*s;
            Real beta01 = -alpha00*s + alpha01*c;
            Real beta10 =  alpha10*c + alpha11*s;
            Real beta11 = -alpha10*s + alpha11*c;

            // A := [c, s; -s, c] B
            alpha00 =  c*beta00 + s*beta10;
            alpha01 =  c*beta01 + s*beta11;
            alpha10 = -s*beta00 + c*beta10;
            alpha11 = -s*beta01 + c*beta11;

            tmp = (alpha00+alpha11)/2;
            alpha00 = alpha11 = tmp;

            if( alpha10 != zero )
            {
                if( alpha01 != zero )
                {
                    if( Sgn(alpha01) == Sgn(alpha10) )
                    {
                        // We can reduce to (real) upper-triangular form
                        Real alpha01Sqrt = Sqrt(Abs(alpha01));
                        Real alpha10Sqrt = Sqrt(Abs(alpha10));
                        Real p = alpha01Sqrt*alpha10Sqrt*Sgn(alpha10);
                        tau = one / Sqrt(Abs(alpha01+alpha10));
                        alpha00 = tmp + p;
                        alpha11 = tmp - p;
                        alpha01 -= alpha10;
                        alpha10 = zero;
                        Real c1 = alpha01Sqrt*tau;
                        Real s1 = alpha10Sqrt*tau;
                        tmp = c*c1 - s*s1;
                        s = c*s1 + s*c1;
                        c = tmp;
                    }
                }
                else
                {
                    alpha01 = -alpha10;
                    alpha10 = zero;
                    tmp = c;
                    c = -s;
                    s = tmp;
                }
            }
        }
    }

    // Explicitly compute the eigenvalues
    lambda0 = alpha00;
    lambda1 = alpha11;
    if( alpha10 != zero )
    {
        lambda0.imag( Sqrt(Abs(alpha01))*Sqrt(Abs(alpha10)) );
        lambda1.imag( -lambda0.imag() );
    }
}

template<typename Real,typename>
void TwoByTwo
( Real& alpha00, Real& alpha01,
  Real& alpha10, Real& alpha11,
  Complex<Real>& lambda0,
  Complex<Real>& lambda1 )
{
    Real c, s;
    TwoByTwo
    ( alpha00, alpha01,
      alpha10, alpha11,
      lambda0, lambda1,
      c, s );
}

template<typename Real>
void TwoByTwo
( Complex<Real>& alpha00, Complex<Real>& alpha01,
  Complex<Real>& alpha10, Complex<Real>& alpha11,
  Complex<Real>& lambda0,
  Complex<Real>& lambda1 )
{
    typedef Complex<Real> F;

    const Real scale = OneAbs(alpha00) + OneAbs(alpha01) +
                       OneAbs(alpha10) + OneAbs(alpha11);
    alpha00 /= scale;
    alpha01 /= scale;
    alpha10 /= scale;
    alpha11 /= scale;
    const F halfTrace = (alpha00+alpha11) / Real(2);
    const F det =
      (alpha00-halfTrace)*(alpha11-halfTrace) - alpha01*alpha10;
    const F discrim = Sqrt( -det );
    lambda0 = (halfTrace+discrim)*scale;
    lambda1 = (halfTrace-discrim)*scale;
}

} // namespace schur
} // namespace El

#endif // ifndef EL_SPECTRAL_SCHUR_HPP
