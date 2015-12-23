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

size_t numLimbs;

} // anonymous namespace

namespace El {
namespace mpc {

mpfr_prec_t Precision()
{ return mpfr_get_default_prec(); }

size_t NumLimbs()
{ return ::numLimbs; }

void SetPrecision( mpfr_prec_t prec )
{
    static bool previouslySet = false;
    if( previouslySet && prec == mpfr_get_default_prec() )
        return;

    if( previouslySet )
        FreeMPI();    
    mpfr_set_default_prec( prec ); 
    RegisterMPI();
    ::numLimbs = (prec-1) / GMP_NUMB_BITS + 1;
    previouslySet = true;
}

mpfr_rnd_t RoundingMode()
{ return mpfr_get_default_rounding_mode(); }

} // namespace mpc

byte* Serialize( Int n, const BigFloat* x, byte* buf )
{
    DEBUG_ONLY(CSE cse("mpc::Serialize"))
    for( Int j=0; j<n; ++j )
        buf = x[j].Serialize( buf );
    return buf;
}

void Serialize( Int n, const BigFloat* x, std::vector<byte>& buf )
{
    DEBUG_ONLY(CSE cse("mpc::Serialize"))
    if( n == 0 )
    {
        buf.resize(0);
        return;
    }
    const auto packedSize = x[0].SerializedSize();
    buf.resize( n*packedSize );
    Serialize( n, x, buf.data() );
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
void Deserialize( Int n, const std::vector<byte>& buf, BigFloat* x )
{ Deserialize( n, buf.data(), x ); }

} // namespace El

#endif // ifdef EL_HAVE_MPC
