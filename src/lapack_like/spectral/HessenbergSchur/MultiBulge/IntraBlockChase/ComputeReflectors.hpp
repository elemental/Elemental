/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SCHUR_HESS_MULTIBULGE_INTRABLOCK_COMPUTE_REFLECTORS_HPP
#define EL_SCHUR_HESS_MULTIBULGE_INTRABLOCK_COMPUTE_REFLECTORS_HPP

namespace El {
namespace hess_schur {
namespace multibulge {
namespace intrablock {

// Compute the Householder reflections needed to chase a set of tightly-packed
// full 4x4 bulges one step downward from the top-left corner of H down to the
// bottom-right corner. See Figure 3 from
// 
//   R. Granat, Bo Kagstrom, and D. Kressner, "LAPACK Working Note #216:
//   A novel parallel QR algorithm for hybrid distributed memory HPC systems",
//
// [CITATION] for a diagram.
//
template<typename F>
void ComputeDiagonalBlockReflectors
(       Int step,
        Int numBulges, 
        Matrix<F>& H,
  const Matrix<Complex<Base<F>>>& shifts,
        Matrix<F>& W,
        bool progress )
{
    DEBUG_CSE
    typedef Base<F> Real;
    const Int n = H.Height();
    DEBUG_ONLY(
      if( step < 0 )
          LogicError("Intrablock chases do not introduce bulges");
    )

    // The top-left 4x4 bulge should be moved from position 0 to position 
    // (n-1) - 3*numBulges.
    const Int numSteps = (n-1) - 3*numBulges;
    if( numSteps <= step )
    {
        W.Resize( 3, 0 ); 
        return;
    }
    Zeros( W, 3, numBulges );

    const Real realZero(0);
    const F zero(0);
    const Real ulp = limits::Precision<Real>();

    // Set aside space for a Householder vector for a candidate reinflation of
    // a deflated bulge
    vector<F> wCand(3);

    for( Int bulge=numBulges-1; bulge>=0; --bulge )
    {
        const Complex<Real> shift0 = shifts(2*bulge);
        const Complex<Real> shift1 = shifts(2*bulge+1);
        const Int bulgeBeg = step + 3*bulge;
        F* w = &W(0,bulge);
        auto ind1 = IR(bulgeBeg+1,bulgeBeg+4);
        auto H11BR = H(ind1,ind1);

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
                    Output
                    ("Unacceptable intrablock replacement bulge at ",bulgeBeg);
                eta10 = beta;
                eta20 = realZero;
                eta30 = realZero;
            }
            else
            {
                // Accept the proposed replacement bulge
                if( progress )
                    Output
                    ("Accepted intrablock replacement bulge at ",bulgeBeg);
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

} // namespace intrablock
} // namespace multibulge
} // namespace hess_schur
} // namespace El

#endif // ifndef EL_SCHUR_HESS_MULTIBULGE_INTRABLOCK_COMPUTE_REFLECTORS_HPP
