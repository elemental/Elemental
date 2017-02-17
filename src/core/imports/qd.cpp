/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>

#ifdef EL_HAVE_QD

namespace {

unsigned oldControlWord=0;

}

namespace El {

void InitializeQD()
{
    fpu_fix_start( &::oldControlWord );
}

void FinalizeQD()
{
    fpu_fix_end( &::oldControlWord );
}

// TODO: Use a single templated implementation?
DoubleDouble::operator long double() const
{
    EL_DEBUG_CSE
    long double alpha = x[0];
    alpha += x[1];
    return alpha;
}

#ifdef EL_HAVE_QUAD
DoubleDouble::operator Quad() const
{
    EL_DEBUG_CSE
    Quad alpha = x[0];
    alpha += x[1];
    return alpha;
}
#endif

// TODO: Use a single templated implementation?
QuadDouble::operator long double() const
{
    EL_DEBUG_CSE
    long double alpha = x[0];
    for( Int j=1; j<4; ++j )
        alpha += x[j];
    return alpha;
}

#ifdef EL_HAVE_QUAD
QuadDouble::operator Quad() const
{
    EL_DEBUG_CSE
    Quad alpha = x[0];
    for( Int j=1; j<4; ++j )
        alpha += x[j];
    return alpha;
}
#endif

#ifdef EL_HAVE_MPC
DoubleDouble::operator BigFloat() const
{
    EL_DEBUG_CSE
    BigFloat alpha( x[0] );
    alpha += x[1];
    return alpha;
}

QuadDouble::operator BigFloat() const
{
    EL_DEBUG_CSE
    BigFloat alpha( x[0] );
    for( Int j=1; j<4; ++j )
        alpha += x[j];
    return alpha;
}
#endif

DoubleDouble::operator long long int() const
{
    // NOTE: This should be verified for correctness and is somewhat ad-hoc
    long long int result = 0;
    DoubleDouble alpha = *this;
    for( Int i=0; i<2; ++i )
    {
        double beta = double(alpha);
        long long int betaInt = static_cast<long long int>(beta);
        result += betaInt;
        alpha -= betaInt;
    }
    return result;
}

QuadDouble::operator long long int() const
{
    // NOTE: This should be verified for correctness and is somewhat ad-hoc
    long long int result = 0;
    QuadDouble alpha = *this;
    for( Int i=0; i<4; ++i )
    {
        double beta = double(alpha);
        long long int betaInt = static_cast<long long int>(beta);
        result += betaInt;
        alpha -= betaInt;
    }
    return result;
}

} // namespace El
#endif // ifdef EL_HAVE_QD
