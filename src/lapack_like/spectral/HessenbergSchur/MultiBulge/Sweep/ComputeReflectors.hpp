/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SCHUR_HESS_MULTIBULGE_SWEEP_COMPUTE_REFLECTORS_HPP
#define EL_SCHUR_HESS_MULTIBULGE_SWEEP_COMPUTE_REFLECTORS_HPP

#include "./IntroduceBulge.hpp"

namespace El {
namespace hess_schur {
namespace multibulge {

template<typename F>
void ComputeReflectors
(       Matrix<F>& H,
        Int winBeg,
        Int winEnd,
  const Matrix<Complex<Base<F>>>& shifts,
        Matrix<F>& W,
        Int packetBeg,
        Int firstBulge,
        Int numBulges,
        bool progress )
{
    DEBUG_CSE
    typedef Base<F> Real;
    const Real realZero(0);
    const F zero(0);
    const Real ulp = limits::Precision<Real>();

    const Int lastBulge = firstBulge + numBulges - 1;
    const Int lastBulgeBeg = packetBeg + 3*lastBulge;
    const bool haveSmallBulge = ( lastBulgeBeg == winEnd-3 );
    DEBUG_ONLY(
      if( lastBulgeBeg > winEnd-3 )
          LogicError("Last bulge starts too late");
    )
    const Int numFullBulges = ( haveSmallBulge ? numBulges-1 : numBulges );

    // Set aside space for a Householder vector for a candidate reinflation of
    // a deflated bulge
    vector<F> wCand(3);

    if( haveSmallBulge )
    {
        const Int bulge = firstBulge + numFullBulges;
        const Complex<Real> shift0 = shifts(2*bulge);
        const Complex<Real> shift1 = shifts(2*bulge+1);
        const Int bulgeBeg = packetBeg + 3*bulge;
        F* w = &W(0,bulge);

        if( bulgeBeg == winBeg-1 )
        {
            auto ind1 = IR(bulgeBeg+1,bulgeBeg+3);
            auto H11BR = H(ind1,ind1);
            IntroduceBulge( H11BR, shift0, shift1, w );
        }
        else
        {
            // Find the reflection for chasing the 3x3 bulge
            F beta = H( bulgeBeg+1, bulgeBeg );
            w[1] = H( bulgeBeg+2, bulgeBeg );
            w[0] = lapack::Reflector( 2, beta, &w[1], 1 );
            H( bulgeBeg+1, bulgeBeg ) = beta;
            H( bulgeBeg+2, bulgeBeg ) = realZero;
        }
    }
    for( Int bulge=firstBulge+numFullBulges-1; bulge>=firstBulge; --bulge )
    {
        const Complex<Real> shift0 = shifts(2*bulge);
        const Complex<Real> shift1 = shifts(2*bulge+1);
        const Int bulgeBeg = packetBeg + 3*bulge;
        F* w = &W(0,bulge);
        auto ind1 = IR(bulgeBeg+1,bulgeBeg+4);
        auto H11BR = H(ind1,ind1);

        if( bulgeBeg == winBeg-1 )
        {
            IntroduceBulge( H11BR, shift0, shift1, w );
        }
        else
        {
            // Prepare to chase the bulge down a step
            F& eta10 = H(bulgeBeg+1,bulgeBeg);
            F& eta20 = H(bulgeBeg+2,bulgeBeg);
            F& eta30 = H(bulgeBeg+3,bulgeBeg);

            F beta = eta10;
            w[1] = eta20;
            w[2] = eta30;
            w[0] = lapack::Reflector( 3, beta, &w[1], 1 );
                    
            // "Vigilantly" search for a deflation
            // (deflation within the interior is exceedingly rare,
            //  but this check should be essentially free)
            const F& eta31 = H(bulgeBeg+3,bulgeBeg+1);
            const F& eta32 = H(bulgeBeg+3,bulgeBeg+2);
            if( eta30 != zero || eta31 != zero || eta32 == zero )
            {
                // Bulge has not collapsed
                eta10 = beta;
                eta20 = realZero;
                eta30 = realZero;
            }
            else
            {
                // Bulge has collapsed
                IntroduceBulge( H11BR, shift0, shift1, wCand.data() );
                F innerProd = wCand[0]*(eta10+Conj(wCand[1])*eta20);

                const F& eta00 = H(bulgeBeg,  bulgeBeg  );
                const F& eta11 = H(bulgeBeg+1,bulgeBeg+1);
                const F& eta22 = H(bulgeBeg+2,bulgeBeg+2);
                if( OneAbs(eta20-wCand[1]*innerProd) +
                    OneAbs(wCand[2]*innerProd) >=
                    ulp*(OneAbs(eta00)+OneAbs(eta11)+OneAbs(eta22)) )
                {
                    // The proposed bulge was unacceptable;
                    // continue using the collapsed one with regret
                    if( progress )
                        Output("Unacceptable replacement bulge at ",bulgeBeg);
                    eta10 = beta;
                    eta20 = realZero;
                    eta30 = realZero;
                }
                else
                {
                    // Accept the proposed replacement bulge
                    if( progress )
                        Output("Accepted replacement bulge at ",bulgeBeg);
                    eta10 -= innerProd;
                    eta20 = realZero; 
                    eta30 = realZero;
                    w[0] = wCand[0];
                    w[1] = wCand[1];
                    w[2] = wCand[2];
                }
            }
        }
    }
}

} // namespace multibulge
} // namespace hess_schur
} // namespace El

#endif // ifndef EL_SCHUR_HESS_MULTIBULGE_SWEEP_COMPUTE_REFLECTORS_HPP
