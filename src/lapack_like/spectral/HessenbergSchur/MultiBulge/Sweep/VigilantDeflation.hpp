/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SCHUR_HESS_MULTIBULGE_SWEEP_VIGILANT_DEFLATION_HPP
#define EL_SCHUR_HESS_MULTIBULGE_SWEEP_VIGILANT_DEFLATION_HPP

namespace El {
namespace hess_schur {
namespace multibulge {

// To be performed during the application of the reflections from the right
template<typename Field>
void VigilantDeflation
( Matrix<Field>& H,
  Int winBeg,
  Int winEnd,
  Int packetBeg,
  Int firstBulge,
  Int numBulges,
  bool progress )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Real zero(0);
    const Field complexZero(0);
    const Real ulp = limits::Precision<Real>();
    const Real safeMin = limits::SafeMin<Real>();
    const Real smallNum = safeMin*(Real(H.Height())/ulp);

    for( Int bulge=firstBulge+numBulges-1; bulge>=firstBulge; --bulge )
    {
        const Int k = Min( winEnd-2, packetBeg+3*bulge );
        Field& eta00 = H(k  ,k  );
        Field& eta01 = H(k  ,k+1);
        Field& eta10 = H(k+1,k  );
        Field& eta11 = H(k+1,k+1);
        const Real eta00Abs = OneAbs(eta00);
        const Real eta10Abs = OneAbs(eta10);
        const Real eta11Abs = OneAbs(eta11);

        if( eta10 == complexZero )
            continue;
        Real localScale = eta00Abs + eta11Abs;
        if( localScale == zero )
        {
            if( k >= winBeg+3 )
                localScale += OneAbs(H(k,k-3));
            if( k >= winBeg+2 )
                localScale += OneAbs(H(k,k-2));
            if( k >= winBeg+1 )
                localScale += OneAbs(H(k,k-1));
            if( k < winEnd-2 )
                localScale += OneAbs(H(k+2,k+1));
            if( k < winEnd-3 )
                localScale += OneAbs(H(k+3,k+1));
            if( k < winEnd-4 )
                localScale += OneAbs(H(k+4,k+1));
        }
        if( eta10Abs <= Max(smallNum,ulp*localScale) )
        {
            const Real eta01Abs = OneAbs(eta01);
            const Real diagDiffAbs = OneAbs(eta00-eta11);
            Real offMax = Max( eta10Abs, eta01Abs );
            Real offMin = Min( eta10Abs, eta01Abs );
            Real diagMax = Max( eta11Abs, diagDiffAbs );
            Real diagMin = Min( eta11Abs, diagDiffAbs );
            Real scale = diagMax + offMax;
            Real localScale2 = diagMin*(diagMax/scale);
            if( localScale2 == zero ||
                offMin*(offMax/scale) <= Max(smallNum,ulp*localScale2) )
            {
                if( progress )
                    Output("Vigilant deflation of H(",k+1,",",k,")");
                eta10 = zero;
            }
        }
    }
}

} // namespace multibulge
} // namespace hess_schur
} // namespace El

#endif // ifndef EL_SCHUR_HESS_MULTIBULGE_SWEEP_VIGILANT_DEFLATION_HPP
