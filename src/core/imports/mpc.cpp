/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
#ifdef EL_HAVE_MPC

namespace {

size_t numLimbs;
int numIntLimbs;

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
        mpi::DestroyBigFloatFamily();
    mpfr_set_default_prec( prec ); 
    ::numLimbs = (prec-1) / GMP_NUMB_BITS + 1;
    mpi::CreateBigFloatFamily();
    previouslySet = true;
}

void SetMinIntBits( int numBits )
{ 
    static bool previouslySet = false;
    const int numIntLimbsNew = (numBits+(GMP_NUMB_BITS-1)) / GMP_NUMB_BITS;
    if( previouslySet && ::numIntLimbs == numIntLimbsNew ) 
        return;

    if( previouslySet )
        mpi::DestroyBigIntFamily();
    ::numIntLimbs = numIntLimbsNew;
    mpi::CreateBigIntFamily();
    previouslySet = true;
}

int NumIntBits()
{ return ::numIntLimbs*GMP_NUMB_BITS; }

int NumIntLimbs()
{ return ::numIntLimbs; }

mpfr_rnd_t RoundingMode()
{ return mpfr_get_default_rounding_mode(); }

Int BinaryToDecimalPrecision( mpfr_prec_t prec )
{ return Int(Floor(prec*std::log10(2.))); }

} // namespace mpc
} // namespace El

#endif // ifdef EL_HAVE_MPC
