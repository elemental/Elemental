/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SCHUR_HESS_MULTIBULGE_INTRABLOCK_APPLY_REFLECTORS_HPP
#define EL_SCHUR_HESS_MULTIBULGE_INTRABLOCK_APPLY_REFLECTORS_HPP

#include "../Sweep/ApplyReflector.hpp"
#include "../Sweep/VigilantDeflation.hpp"

namespace El {
namespace hess_schur {
namespace multibulge {
namespace intrablock {

template<typename F>
void ApplyReflectorsToDiagonalBlock
( Int step,
  Int numBulges,
  Matrix<F>& H,
  Matrix<F>& U,
  Matrix<F>& W,
  bool progress )
{
    DEBUG_CSE
    const Int n = H.Height();
    DEBUG_ONLY(
      if( step < 0 )
          LogicError("Intrablock chasing does not create bulges");
      if( step + 3*numBulges >= n )
          LogicError("Attempted to chase bulge outside of diagonal block");
    )

    // Apply from the left
    // ===================
    for( Int j=0; j<n; ++j )
    {
        // Avoid re-applying freshly generated reflections
        // (which is only relevant when (j-step) % 3 = 0)
        const Int applyBulgeEnd = Min( numBulges, (j-step+2)/3 );
        for( Int bulge=0; bulge<applyBulgeEnd; ++bulge )
        {
            const Int bulgeBeg = step + 3*bulge;
            const F* w = &W(0,bulge);
            ApplyLeftReflector
            ( H(bulgeBeg+1,j), H(bulgeBeg+2,j), H(bulgeBeg+3,j), w );
        }
    }

    // Apply from the right (excluding the fourth row to support vig. deflation)
    // =========================================================================
    for( Int bulge=numBulges-1; bulge>=0; --bulge )
    {
        const Int bulgeBeg = step + 3*bulge;
        const F* w = &W(0,bulge);
        for( Int i=0; i<n; ++i )
            ApplyRightReflector
            ( H(i,bulgeBeg+1), H(i,bulgeBeg+2), H(i,bulgeBeg+3), w );
        for( Int i=0; i<n-1; ++i ) 
            ApplyRightReflector
            ( U(i,bulgeBeg+1), U(i,bulgeBeg+2), U(i,bulgeBeg+3), w );
    }

    // Vigilant deflation check using Ahues/Tisseur criteria
    // =====================================================
    VigilantDeflation
    ( H, 0 /*winBeg*/, n /*winEnd*/, 0 /*packetBeg*/,
      0 /*firstVigBulge*/, numBulges /*numVigBulges*/, progress );

    // Form the last row of the result of applying from the right:
    //
    // | X X X X X |     | X X X X X |
    // | X X X X X |     | X X X X X |
    // | X X X X X | |-> |   X X X X |
    // | X X X X X |     |   X X X X |
    // |       X X |     |   X X X X |
    //
    // The last row is introduced from the transformation
    //
    //  H(k+4,k+1:k+3) -= conj(tau) H(k+4,k+1:k+3) [1; nu1; nu2] [1; nu1; nu2]'.
    //
    // For convenience, since H(k+4,k+1:k+2)=0, we start by computing 
    //
    //   innerProd = conj(tau) H(k+4,k+1:k+3) [1; nu1; nu2]
    //             = conj(tau) H(k+4,k+3) nu2
    //
    // TODO(poulson): Avoid memory allocations for several of the temporaries
    // below in the case where we have heap scalars (e.g., BigFloat).
    const Int lastRowBulgeEnd = Min( numBulges, (n-step-2)/3 );
    for( Int bulge=lastRowBulgeEnd-1; bulge>=0; --bulge )
    {
        const Int k = step + 3*bulge;
        const F* w = &W(0,bulge);
        const F& tau = w[0];
        const F& nu1 = w[1];
        const F& nu2 = w[2];

        const F innerProd = Conj(tau)*H(k+4,k+3)*nu2;
        H(k+4,k+1) = -innerProd;
        H(k+4,k+2) = -innerProd*Conj(nu1);
        H(k+4,k+3) -= innerProd*Conj(nu2);
    }
}

// Unfortunately, it seems to be the case that it is noticeably faster
// for this routine to manually inline the data access than to use the 
// (presumably inlined) Matrix::operator()(int,int) calls.
template<typename F>
void ApplyReflectorsToDiagonalBlockOpt
( Int step,
  Int numBulges,
  Matrix<F>& H,
  Matrix<F>& U,
  Matrix<F>& W,
  bool progress )
{
    DEBUG_CSE
    const Int n = H.Height();
    DEBUG_ONLY(
      if( step < 0 )
          LogicError("Intrablock chasing does not create bulges");
      if( step + 3*numBulges >= n )
          LogicError("Attempted to chase bulge outside of diagonal block");
    )
    F* HBuf = H.Buffer();
    F* UBuf = U.Buffer();
    F* WBuf = W.Buffer();
    const Int HLDim = H.LDim();
    const Int ULDim = U.LDim();
    const Int WLDim = W.LDim();

    // Apply from the left
    // ===================
    for( Int j=0; j<n; ++j )
    {
        // Avoid re-applying freshly generated reflections
        // (which is only relevant when (j-step) % 3 = 0)
        const Int applyBulgeEnd = Min( numBulges, (j-step+2)/3 );
        for( Int bulge=0; bulge<applyBulgeEnd; ++bulge )
        {
            const Int bulgeBeg = step + 3*bulge;
            const F* w = &WBuf[bulge*WLDim];
            ApplyLeftReflector
            ( HBuf[(bulgeBeg+1)+j*HLDim],
              HBuf[(bulgeBeg+2)+j*HLDim],
              HBuf[(bulgeBeg+3)+j*HLDim], w );
        }
    }

    // Apply from the right (excluding the fourth row to support vig. deflation)
    // =========================================================================
    for( Int bulge=numBulges-1; bulge>=0; --bulge )
    {
        const Int bulgeBeg = step + 3*bulge;
        const F* w = &WBuf[bulge*WLDim];
        for( Int i=0; i<n; ++i )
            ApplyRightReflector
            ( HBuf[i+(bulgeBeg+1)*HLDim],
              HBuf[i+(bulgeBeg+2)*HLDim],
              HBuf[i+(bulgeBeg+3)*HLDim], w );
        for( Int i=0; i<n-1; ++i ) 
            ApplyRightReflector
            ( UBuf[i+(bulgeBeg+1)*ULDim],
              UBuf[i+(bulgeBeg+2)*ULDim],
              UBuf[i+(bulgeBeg+3)*ULDim], w );
    }

    // Vigilant deflation check using Ahues/Tisseur criteria
    // =====================================================
    VigilantDeflation
    ( H, 0 /*winBeg*/, n /*winEnd*/, 0 /*packetBeg*/,
      0 /*firstVigBulge*/, numBulges /*numVigBulges*/, progress );

    // Form the last row of the result of applying from the right:
    //
    // | X X X X X |     | X X X X X |
    // | X X X X X |     | X X X X X |
    // | X X X X X | |-> |   X X X X |
    // | X X X X X |     |   X X X X |
    // |       X X |     |   X X X X |
    //
    // The last row is introduced from the transformation
    //
    //  H(k+4,k+1:k+3) -= conj(tau) H(k+4,k+1:k+3) [1; nu1; nu2] [1; nu1; nu2]'.
    //
    // For convenience, since H(k+4,k+1:k+2)=0, we start by computing 
    //
    //   innerProd = conj(tau) H(k+4,k+1:k+3) [1; nu1; nu2]
    //             = conj(tau) H(k+4,k+3) nu2
    //
    // TODO(poulson): Avoid memory allocations for several of the temporaries
    // below in the case where we have heap scalars (e.g., BigFloat).
    const Int lastRowBulgeEnd = Min( numBulges, (n-step-2)/3 );
    for( Int bulge=lastRowBulgeEnd-1; bulge>=0; --bulge )
    {
        const Int k = step + 3*bulge;
        const F* w = &W(0,bulge);
        const F& tau = w[0];
        const F& nu1 = w[1];
        const F& nu2 = w[2];

        const F innerProd = Conj(tau)*HBuf[(k+4)+(k+3)*HLDim]*nu2;
        HBuf[(k+4)+(k+1)*HLDim] = -innerProd;
        HBuf[(k+4)+(k+2)*HLDim] = -innerProd*Conj(nu1);
        HBuf[(k+4)+(k+3)*HLDim] -= innerProd*Conj(nu2);
    }
}

} // namespace intrablock
} // namespace multibulge
} // namespace hess_schur
} // namespace El

#endif // ifndef EL_SCHUR_HESS_MULTIBULGE_INTRABLOCK_APPLY_REFLECTORS_HPP
