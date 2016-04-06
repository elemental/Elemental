/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#ifdef EL_HAVE_QD

namespace El {

#ifdef EL_HAVE_QUAD
DoubleDouble::DoubleDouble( const Quad& a )
{
    DEBUG_ONLY(CSE cse("DoubleDouble::DoubleDouble [Quad]"))
    x[0] = a;
    x[1] = (a-x[0]);
}

QuadDouble::QuadDouble( const Quad& a )
{
    DEBUG_ONLY(CSE cse("QuadDouble::QuadDouble [Quad]"))
    Quad b = a;
    x[0] = b;
    for( Int j=1; j<4; ++j )
    {
        b -= x[j-1];
        x[j] = b;
    }
}

DoubleDouble::operator Quad() const
{
    DEBUG_ONLY(CSE cse("DoubleDouble::operator Quad"))
    return Quad(x[0]) + x[1];
}

QuadDouble::operator Quad() const
{
    DEBUG_ONLY(CSE cse("QuadDouble::operator Quad"))
    Quad alpha = x[0];
    for( Int j=1; j<4; ++j )
        alpha += x[j];
    return alpha;
}
#endif

#ifdef EL_HAVE_MPC
DoubleDouble::operator BigFloat() const
{
    DEBUG_ONLY(CSE cse("DoubleDouble::operator BigFloat"))
    BigFloat alpha( x[0] );
    alpha += x[1];
    return alpha;
}

QuadDouble::operator BigFloat() const
{
    DEBUG_ONLY(CSE cse("QuadDouble::operator BigFloat"))
    BigFloat alpha( x[0] );
    alpha += x[1];
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
