/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

#include <El/blas_like/level2.hpp>

namespace El {
namespace trdtrmm {

template<typename F>
void LUnblocked( Matrix<F>& L, bool conjugate=false )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( L.Height() != L.Width() )
          LogicError("L must be square");
    )
    const Int n = L.Height();
    const Int ldim = L.LDim();

    Matrix<F> s10;

    for( Int k=0; k<n; ++k )
    {
        auto L00 = L( IR(0,k), IR(0,k) );
        auto l10 = L( IR(k),   IR(0,k) );

        // S10 := L10
        s10 = l10;

        // L10 := L10 / delta11
        const F deltaInv = F(1)/L.Get(k,k);
        l10 *= deltaInv;

        // L00 := L00 + l10' s10
        const F* l10Buf = l10.LockedBuffer();
        if( conjugate )
        {
            for( Int j=0; j<k; ++j )
            {
                F* L00Col = L00.Buffer(0,j);
                const F gamma = s10.Get(0,j);
                for( Int i=j; i<k; ++i )
                    L00Col[i] += Conj(l10Buf[i*ldim])*gamma;
            }
        }
        else
        {
            for( Int j=0; j<k; ++j )
            {
                F* L00Col = L00.Buffer(0,j);
                const F gamma = s10.Get(0,j);
                for( Int i=j; i<k; ++i )
                    L00Col[i] += l10Buf[i*ldim]*gamma;
            }
        }

        // L11 := 1 / delta11
        L.Set( k, k, deltaInv );
    }
}

template<typename F>
void LUnblocked( Matrix<F>& L, const Matrix<F>& dSub, bool conjugate=false )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( L.Height() != L.Width() )
          LogicError("L must be square");
    )
    const Int n = L.Height();
    const Int ldim = L.LDim();
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );

    Matrix<F> s10, S10, D11(2,2), D11Inv(2,2);

    Int k=0;
    while( k < n )
    {
        const Int nb = ( k<n-1 && dSub.Get(k,0) != F(0) ? 2 : 1 );

        if( nb == 1 )
        {
            auto L00 = L( IR(0,k),    IR(0,k) );
            auto l10 = L( IR(k,k+nb), IR(0,k) );

            // S10 := L10
            s10 = l10;

            // L10 := L10 / delta11
            const F deltaInv = F(1)/L.Get(k,k);
            l10 *= deltaInv;

            // L00 := L00 + l10' s10
            // TODO: Extend Trr for this case and then switch
            const F* l10Buf = l10.LockedBuffer();
            if( conjugate )
            {
                for( Int j=0; j<k; ++j )
                {
                    F* L00Col = L00.Buffer(0,j);
                    const F gamma = s10.Get(0,j);
                    for( Int i=j; i<k; ++i )
                        L00Col[i] += Conj(l10Buf[i*ldim])*gamma;
                }
            }
            else
            {
                for( Int j=0; j<k; ++j )
                {
                    F* L00Col = L00.Buffer(0,j);
                    const F gamma = s10.Get(0,j);
                    for( Int i=j; i<k; ++i )
                        L00Col[i] += l10Buf[i*ldim]*gamma;
                }
            }

            // lambda11 := 1 / delta11
            L.Set( k, k, deltaInv );
        }
        else
        {
            auto L00 = L( IR(0,k),    IR(0,k)    );
            auto L10 = L( IR(k,k+nb), IR(0,k)    );
            auto L11 = L( IR(k,k+nb), IR(k,k+nb) );

            // S10 := L10
            S10 = L10;

            // L10 := inv(D11) L10 
            D11.Set( 0, 0, L11.Get(0,0) );
            D11.Set( 1, 1, L11.Get(1,1) );
            D11.Set( 1, 0, dSub.Get(k,0) );

            D11Inv = D11;
            Symmetric2x2Inv( LOWER, D11Inv, conjugate );
            MakeSymmetric( LOWER, D11Inv, conjugate );
            Transform2x2Rows( D11Inv, L10, 0, 1 );

            // L00 := L00 + L10' S10
            // TODO: Custom rank-2 update 
            Trrk( LOWER, orientation, NORMAL, F(1), L10, S10, F(1), L00 );

            // L11 := inv(D11)
            L11.Set( 0, 0, D11Inv.Get(0,0) );
            L11.Set( 1, 0, D11Inv.Get(1,0) );
            L11.Set( 1, 1, D11Inv.Get(1,1) );
        }

        k += nb;
    }
}

template<typename F>
void UUnblocked( Matrix<F>& U, bool conjugate=false )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( U.Height() != U.Width() )
          LogicError("U must be square");
    )
    const Int n = U.Height();

    Matrix<F> s01;

    for( Int k=0; k<n; ++k )
    {
        auto U00 = U( IR(0,k), IR(0,k) );
        auto u01 = U( IR(0,k), IR(k)   );

        s01 = u01;

        // u01 := u01 / delta11
        const F deltaInv = F(1)/U.Get(k,k);
        u01 *= deltaInv;

        // U00 := U00 + s01 u01'
        Trr( UPPER, F(1), s01, u01, U00, conjugate );

        // lambda11 := 1 / delta11
        U.Set( k, k, deltaInv );
    }
}

} // namespace trdtrmm
} // namespace El
