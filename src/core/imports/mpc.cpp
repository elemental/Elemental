/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
#ifdef EL_HAVE_MPC

namespace {

// To ease the initial integration of MPFR and MPC into Elemental, it will
// initially be required that all MPFR/MPC types use a globally-constant
// precision
mpfr_prec_t mpfrPrec=1024;

} // anonymous namespace

namespace El {
namespace mpc {

mpfr_prec_t Precision()
{ return ::mpfrPrec; }

void SetPrecision( mpfr_prec_t precision )
{ ::mpfrPrec = precision; }

mpfr_rnd_t RoundingMode()
{ return mpfr_get_default_rounding_mode(); }

byte* Serialize( Int n, const BigFloat* x, byte* buf )
{
    DEBUG_ONLY(CSE cse("mpc::Serialize"))
    for( Int j=0; j<n; ++j )
        buf = x[j].Serialize( buf );
    return buf;
}

// We assume that the BigFloat's in x are already constructed
byte* Deserialize( Int n, byte* buf, BigFloat* x )
{
    DEBUG_ONLY(CSE cse("mpc::Deserialize"))
    for( Int j=0; j<n; ++j )
        buf = x[j].Deserialize( buf );
    return buf;
}
const byte* Deserialize( Int n, const byte* buf, BigFloat* x )
{
    DEBUG_ONLY(CSE cse("mpc::Deserialize"))
    for( Int j=0; j<n; ++j )
        buf = x[j].Deserialize( buf );
    return buf;
}

} // namespace mpc
} // namespace El

#endif // ifdef EL_HAVE_MPC
