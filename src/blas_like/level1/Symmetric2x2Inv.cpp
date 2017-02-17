/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>
#include <El/blas_like/level1.hpp>

namespace El {

template<typename Field>
void Symmetric2x2Inv( UpperOrLower uplo, Matrix<Field>& D, bool conjugate )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    if( uplo == LOWER )
    {
        if( conjugate )
        {
            const Real delta11 = RealPart(D(0,0));
            const Field delta21 = D(1,0);
            const Real delta22 = RealPart(D(1,1));
            const Real delta21Abs = SafeAbs( delta21 );
            const Real phi21To11 = delta22 / delta21Abs;
            const Real phi21To22 = delta11 / delta21Abs;
            const Field phi21 = delta21 / delta21Abs;
            const Real xi = (Real(1)/(phi21To11*phi21To22-Real(1)))/delta21Abs;

            D.SetRealPart( 0, 0,  xi*phi21To11 );
            D(1,0) = -xi*phi21;
            D.SetRealPart( 1, 1,  xi*phi21To22 );
        }
        else
        {
            const Field delta11 = D(0,0);
            const Field delta21 = D(1,0);
            const Field delta22 = D(1,1);
            const Field chi21To11 = -delta22 / delta21;
            const Field chi21To22 = -delta11 / delta21;
            const Field chi21 =
              (Field(1)/(Field(1)-chi21To11*chi21To22))/delta21;

            D(0,0) = chi21*chi21To11;
            D(1,0) = chi21;
            D(1,1) = chi21*chi21To22;
        }
    }
    else
        LogicError("This option not yet supported");
}

#define PROTO(Field) \
  template void Symmetric2x2Inv \
  ( UpperOrLower uplo, Matrix<Field>& A, bool conjugate );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
