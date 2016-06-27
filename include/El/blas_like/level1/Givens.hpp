/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_GIVENS_HPP
#define EL_BLAS_GIVENS_HPP

namespace El {

// CITATION
//
// Please see LAPACK Working Note 148,
//
//   D. Bindel, J. Demmel, W. Kahan, and O. Marques,
//   "On computing Givens rotations reliably and efficiently",
//   http://www.netlib.org/lapack/lawnspdf/lawn148.pdf
//
// which resulted in the LAPACK routines {s,d,c,z}lartg. But note
// that the LAPACK implementations slightly differ from said working
// note in that they round z^2 to the nearest radix rather than z
// (which results in a different result with float, but the same with
// double).
//
template<typename Real>
Real Givens( const Real& chi0, const Real& chi1, Real& c, Real& s )
{
    const Real zero(0), one(1);
    if( chi1 == zero ) 
    {
        c = 1;
        s = 0; 
        return chi0;
    }
    else if( chi0 == zero )
    {       
        c = 0;
        s = 1;
        return chi1;
    }
    
    const Real safeMin = limits::SafeMin(one);
    const Real safeMinToSquare = limits::SafeMinToSquare(one);
    const Real safeMaxToSquare = one / safeMinToSquare;

    Real scale = Max( Abs(chi0), Abs(chi1) );
    Real chi0Scaled = chi0;
    Real chi1Scaled = chi1;
    Int rescaleCounter = 0;
    if( scale >= safeMaxToSquare )
    {
        while( scale >= safeMaxToSquare )
        {
            chi0Scaled *= safeMinToSquare;
            chi1Scaled *= safeMinToSquare;
            scale *= safeMinToSquare;
            ++rescaleCounter; 
        }
    }
    else if( scale <= safeMinToSquare )
    {
        if( chi1 == zero || !limits::IsFinite(Abs(chi1)) )
        {
            c = 1;
            s = 0;
            return chi0;
        }
        while( scale <= safeMinToSquare )
        {
            chi0Scaled *= safeMaxToSquare;
            chi1Scaled *= safeMaxToSquare;
            scale *= safeMaxToSquare;
            --rescaleCounter;
        }
    }
    Real rho = Sqrt( chi0Scaled*chi0Scaled + chi1Scaled*chi1Scaled );
    c = chi0Scaled/rho;
    s = chi1Scaled/rho;
    if( rescaleCounter > 0 )
    {
        for( Int i=0; i<rescaleCounter; ++i )
            rho *= safeMaxToSquare;
    }
    else if( rescaleCounter < 0 )
    {
        for( Int i=0; i<-rescaleCounter; ++i )
            rho *= safeMinToSquare;
    }
    if( Abs(chi0) > Abs(chi1) && c < zero )
    {
        c = -c;
        s = -s;
        rho = -rho;
    }
    return rho;
}

template<typename Real>
Complex<Real> Givens
( const Complex<Real>& chi0,
  const Complex<Real>& chi1,
  Real& c,
  Complex<Real>& s )
{
    typedef Complex<Real> F;
    const Real one(1);
    const F zeroF(0);

    const Real safeMin = limits::SafeMin(one);
    const Real safeMinToSquare = limits::SafeMinToSquare(one);
    const Real safeMaxToSquare = one / safeMinToSquare;

    Real scale = Max( MaxAbs(chi0), MaxAbs(chi1) );
    F chi0Scaled = chi0;
    F chi1Scaled = chi1;
    Int rescaleCounter = 0;
    if( scale >= safeMaxToSquare )
    {
        while( scale >= safeMaxToSquare )
        {
            chi0Scaled *= safeMinToSquare;
            chi1Scaled *= safeMinToSquare;
            scale *= safeMinToSquare;
            ++rescaleCounter;
        }
    }
    else if( scale <= safeMinToSquare )
    {
        if( chi1 == zeroF || !limits::IsFinite(Abs(chi1)) )
        {
            c = 1;
            s = 0;
            return chi0;
        }
        while( scale <= safeMinToSquare )
        {
            chi0Scaled *= safeMaxToSquare;
            chi1Scaled *= safeMaxToSquare;
            scale *= safeMaxToSquare;
            --rescaleCounter;
        }
    }
    const Real chi0AbsSquare =
      chi0.real()*chi0.real() + chi0.imag()*chi0.imag();
    const Real chi1AbsSquare =
      chi1.real()*chi1.real() + chi1.imag()*chi1.imag();
    if( chi0AbsSquare <= Max(chi1AbsSquare,one)*safeMin )
    {
        // Handle the exceptional case where chi0 is very small
        if( chi0 == zeroF )
        {
            c = 0;
            s = Conj(chi1Scaled) / SafeAbs(chi0Scaled);
            return SafeAbs(chi1);
        }
        const Real chi0ScaledAbsSquare = SafeAbs(chi0Scaled);
        // By the nature of our initial rescaling, chi1 must be substantially
        // larger than chi0 and is in a safe range
        const Real chi1ScaledAbsSquare = Sqrt(chi1AbsSquare);
        c = chi0ScaledAbsSquare / chi1ScaledAbsSquare;
        F phase;
        if( MaxAbs(chi0) > one )
        {
            phase = chi0 / SafeAbs(chi0);
        }
        else
        {
            Complex<Real> delta = safeMaxToSquare*chi0;
            phase = delta / SafeAbs(delta);
        }
        s = phase*(Conj(chi1Scaled)/chi0ScaledAbsSquare);
        return c*chi0 + s*chi1;
    }
    else
    {
        // Typical branch where chi0 and chi1 are in a safe range
        const Real tau = Sqrt(one+(chi1AbsSquare/chi0AbsSquare));
        F rho = tau*chi0Scaled;
        c = one / tau;
        const Real delta = chi0AbsSquare + chi1AbsSquare;
        s = (rho/delta)*Conj(chi1Scaled);
        if( rescaleCounter > 0 )
        {
            for( Int i=0; i<rescaleCounter; ++i )
                rho *= safeMaxToSquare;
        }
        else if( rescaleCounter < 0 )
        {
            for( Int i=0; i<-rescaleCounter; ++i )
                rho *= safeMinToSquare;
        }
        return rho;
    }
}

} // namespace El

#endif // ifndef EL_BLAS_GIVENS_HPP
