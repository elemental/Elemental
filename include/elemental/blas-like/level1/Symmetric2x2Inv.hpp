/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_SYMMETRIC2X2INV_HPP
#define ELEM_SYMMETRIC2X2INV_HPP

namespace elem {

template<typename F>
inline void
Symmetric2x2Inv( UpperOrLower uplo, Matrix<F>& D, bool conjugate=false )
{
    DEBUG_ONLY(CallStackEntry cse("Symmetric2x2Inv"))
    typedef Base<F> Real;
    if( uplo == LOWER )
    {
        if( conjugate )
        {
            const Real delta11 = D.GetRealPart(0,0);
            const F delta21 = D.Get(1,0);
            const Real delta22 = D.GetRealPart(1,1);
            const Real delta21Abs = SafeAbs( delta21 );
            const Real phi21To11 = delta22 / delta21Abs;
            const Real phi21To22 = delta11 / delta21Abs;
            const F phi21 = delta21 / delta21Abs;
            const Real xi = (Real(1)/(phi21To11*phi21To22-Real(1)))/delta21Abs;

            D.SetRealPart( 0, 0,  xi*phi21To11 );
            D.Set(         1, 0, -xi*phi21     );
            D.SetRealPart( 1, 1,  xi*phi21To22 );
        }
        else
        {
            const F delta11 = D.Get(0,0);
            const F delta21 = D.Get(1,0);
            const F delta22 = D.Get(1,1);
            const F chi21To11 = -delta22 / delta21;
            const F chi21To22 = -delta11 / delta21;
            const F chi21 = (F(1)/(F(1)-chi21To11*chi21To22))/delta21;

            D.Set( 0, 0, chi21*chi21To11 );
            D.Set( 1, 0, chi21           );
            D.Set( 1, 1, chi21*chi21To22 );
        }
    }
    else
        LogicError("This option not yet supported");
}

} // namespace elem

#endif // ifndef ELEM_SYMMETRIC2X2INV_HPP
