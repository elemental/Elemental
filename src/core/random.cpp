/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>

namespace {

// A common Mersenne twister configuration
std::mt19937 generator;

#ifdef EL_HAVE_MPC
gmp_randstate_t gmpRandState;
#endif

}

namespace El {

void InitializeRandom( bool deterministic )
{
    const unsigned rank = mpi::Rank( mpi::COMM_WORLD );
    const long secs = ( deterministic ? 21 : time(NULL) );
    const long seed = (secs<<16) | (rank & 0xFFFF);

    ::generator.seed( seed );

    srand( seed );

#ifdef EL_HAVE_MPC
    mpfr::SetMinIntBits( 256 );
    mpfr::SetPrecision( 256 );
    gmp_randinit_default( ::gmpRandState );
    gmp_randseed_ui( ::gmpRandState, seed );
#endif
}

void FinalizeRandom()
{
#ifdef EL_HAVE_MPC
    gmp_randclear( ::gmpRandState );
#endif
}

std::mt19937& Generator()
{ return ::generator; }

#ifdef EL_HAVE_MPC
namespace mpfr {

void RandomState( gmp_randstate_t randState )
{
    // It is surprisingly tedious to return the state...

    randState->_mp_seed->_mp_alloc = ::gmpRandState->_mp_seed->_mp_alloc;
    randState->_mp_seed->_mp_size = ::gmpRandState->_mp_seed->_mp_size;
    randState->_mp_seed->_mp_d = ::gmpRandState->_mp_seed->_mp_d;

    randState->_mp_alg = ::gmpRandState->_mp_alg;
    randState->_mp_algdata._mp_lc = ::gmpRandState->_mp_algdata._mp_lc;
}

} // namespace mpfr
#endif

bool BooleanCoinFlip()
{ return SampleUniform<double>(0,1) >= 0.5; }

Int CoinFlip() { return ( BooleanCoinFlip() ? 1 : -1 ); }

#ifdef EL_HAVE_QUAD
template<>
Quad SampleUniform( const Quad& a, const Quad& b )
{
    // Quad is typically a first-class scalar for STL (at least for GCC),
    // but not in this case
    return SampleUniformNaive( a, b );
}
#endif

#ifdef EL_HAVE_MPC
template<>
BigFloat SampleUniform( const BigFloat& a, const BigFloat& b )
{
    BigFloat sample;
    gmp_randstate_t randState;
    mpfr::RandomState( randState );

    while( 1 )
    {
        const int ret = mpfr_urandomb( sample.Pointer(), randState );
        if( ret == 0 )
            break;
    }
    return a + sample*(b-a);
}

template<>
BigInt SampleUniform( const BigInt& a, const BigInt& b )
{
    BigInt sample;
    gmp_randstate_t randState;
    mpfr::RandomState( randState );

    mpz_urandomb( sample.Pointer(), randState, b.NumBits() );
    return a+Mod(sample,b-a);
}

#endif // ifdef EL_HAVE_MPC

template<>
Int SampleUniform( const Int& a, const Int& b )
{
#ifdef EL_HAVE_CXX11RANDOM
    std::mt19937& gen = Generator();
    std::uniform_int_distribution<Int> intDist(a,b-1);
    return intDist(gen);
#else
    return a + (rand() % (b-a));
#endif
}
// TODO: BigInt version?

#ifdef EL_HAVE_QUAD
template<>
Quad SampleNormal( const Quad& mean, const Quad& stddev )
{
    // Quad is typically a first-class scalar for STL (at least for GCC),
    // but not in this case
    return SampleNormalMarsiglia( mean, stddev );
}

template<>
Complex<Quad> SampleNormal( const Complex<Quad>& mean, const Quad& stddev )
{
    // Quad is typically a first-class scalar for STL (at least for GCC),
    // but not in this case
    return SampleNormalMarsiglia( mean, stddev );
}
#endif

} // namespace El
