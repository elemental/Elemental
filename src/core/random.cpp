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

#ifdef EL_HAVE_QD
template<>
DoubleDouble SampleUniform( const DoubleDouble& a, const DoubleDouble& b )
{
    DoubleDouble sample;
    // TODO: Use a better random number generator
#ifdef EL_HAVE_CXX11RANDOM
    std::mt19937& gen = Generator();
    std::uniform_real_distribution<double> uni((double)a,(double)b);
    sample = uni(gen);
#else
    sample = (DoubleDouble(rand())/(DoubleDouble(RAND_MAX)+1))*(b-a) + a;
#endif
    return sample;
}

template<>
QuadDouble SampleUniform( const QuadDouble& a, const QuadDouble& b )
{
    QuadDouble sample;
    // TODO: Use a better random number generator
#ifdef EL_HAVE_CXX11RANDOM
    std::mt19937& gen = Generator();
    std::uniform_real_distribution<double> uni((double)a,(double)b);
    sample = uni(gen);
#else
    sample = (QuadDouble(rand())/(QuadDouble(RAND_MAX)+1))*(b-a) + a;
#endif
    return sample;
}

#endif

#ifdef EL_HAVE_QUAD
template<>
Quad SampleUniform( const Quad& a, const Quad& b )
{
    Quad sample;
#ifdef EL_HAVE_CXX11RANDOM
    std::mt19937& gen = Generator();
    std::uniform_real_distribution<long double> 
      uni((long double)a,(long double)b);
    sample = uni(gen);
#else
    sample = (Quad(rand())/(Quad(RAND_MAX)+1))*(b-a) + a;
#endif
    return sample;
}

template<>
Complex<Quad> SampleUniform( const Complex<Quad>& a, const Complex<Quad>& b )
{
    Complex<Quad> sample;
    sample.real( SampleUniform( a.real(), b.real() ) );
    sample.imag( SampleUniform( a.imag(), b.imag() ) );
    return sample;
}
#endif // ifdef EL_HAVE_QUAD

#ifdef EL_HAVE_MPC
template<>
BigInt SampleUniform( const BigInt& a, const BigInt& b )
{
    BigInt sample; 
    gmp_randstate_t randState;
    mpfr::RandomState( randState );
    
    mpz_urandomb( sample.Pointer(), randState, b.NumBits() );
    return a+Mod(sample,b-a);
}

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
Complex<BigFloat> SampleUniform
( const Complex<BigFloat>& a, const Complex<BigFloat>& b )
{
    Complex<BigFloat> sample;
    sample.real( SampleUniform( a.real(), b.real() ) );
    sample.imag( SampleUniform( a.imag(), b.imag() ) );
    return sample;
}
#endif // ifdef EL_HAVE_MPC

template<>
Int SampleUniform<Int>( const Int& a, const Int& b )
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

#ifdef EL_HAVE_QD
template<>
DoubleDouble
SampleNormal( const DoubleDouble& mean, const DoubleDouble& stddev )
{
    DoubleDouble sample;

#ifdef EL_HAVE_CXX11RANDOM
    // TODO: Use a better RNG
    std::mt19937& gen = Generator();
    std::normal_distribution<double> normal((double)mean,(double)stddev);
    sample = normal(gen);
#else
    // Run Marsiglia's polar method
    // ============================
    // NOTE: Half of the generated samples are thrown away in the case that
    //       F is real.
    while( true )
    {
        const DoubleDouble U = SampleUniform<DoubleDouble>(-1,1);
        const DoubleDouble V = SampleUniform<DoubleDouble>(-1,1);
        const DoubleDouble S = Sqrt(U*U+V*V);
        if( S > DoubleDouble(0) && S < DoubleDouble(1) )
        {
            const DoubleDouble W = Sqrt(-2*Log(S)/S);
            sample = mean + stddev*U*W;
            break;
        }
    }
#endif

    return sample;
}

template<>
QuadDouble SampleNormal( const QuadDouble& mean, const QuadDouble& stddev )
{
    QuadDouble sample;

#ifdef EL_HAVE_CXX11RANDOM
    // TODO: Use a better RNG
    std::mt19937& gen = Generator();
    std::normal_distribution<double> normal((double)mean,(double)stddev);
    sample = normal(gen);
#else
    // Run Marsiglia's polar method
    // ============================
    // NOTE: Half of the generated samples are thrown away in the case that
    //       F is real.
    while( true )
    {
        const QuadDouble U = SampleUniform<QuadDouble>(-1,1);
        const QuadDouble V = SampleUniform<QuadDouble>(-1,1);
        const QuadDouble S = Sqrt(U*U+V*V);
        if( S > QuadDouble(0) && S < QuadDouble(1) )
        {
            const QuadDouble W = Sqrt(-2*Log(S)/S);
            sample = mean + stddev*U*W;
            break;
        }
    }
#endif

    return sample;
}
#endif

#ifdef EL_HAVE_QUAD
template<>
Quad SampleNormal( const Quad& mean, const Quad& stddev )
{
    Quad sample;

#ifdef EL_HAVE_CXX11RANDOM
    std::mt19937& gen = Generator();
    std::normal_distribution<long double> 
      normal( (long double)mean, (long double)stddev );
    sample = normal(gen);
#else
    // Run Marsiglia's polar method
    // ============================
    // NOTE: Half of the generated samples are thrown away in the case that
    //       F is real.
    while( true )
    {
        const Quad U = SampleUniform<Quad>(-1,1);
        const Quad V = SampleUniform<Quad>(-1,1);
        const Quad S = Sqrt(U*U+V*V);
        if( S > Quad(0) && S < Quad(1) )
        {
            const Quad W = Sqrt(-2*Log(S)/S);
            sample = mean + stddev*U*W;
            break;
        }
    }
#endif

    return sample;
}

template<>
Complex<Quad> SampleNormal( const Complex<Quad>& mean, const Quad& stddev )
{
    Complex<Quad> sample;
    Quad stddevAdj = stddev;
    stddevAdj /= Sqrt(Quad(2));

    sample.real( SampleNormal( mean.real(), stddevAdj ) );
    sample.imag( SampleNormal( mean.imag(), stddevAdj ) );
    return sample;
}
#endif // ifdef EL_HAVE_QUAD

#ifdef EL_HAVE_MPC
template<>
BigFloat SampleNormal( const BigFloat& mean, const BigFloat& stddev )
{
    BigFloat sample;

    // Run Marsiglia's polar method
    // ============================
    // NOTE: Half of the generated samples are thrown away in the case that
    //       F is real.
    while( true )
    {
        const BigFloat U = SampleUniform<BigFloat>(-1,1);
        const BigFloat V = SampleUniform<BigFloat>(-1,1);
        const BigFloat S = Sqrt(U*U+V*V);
        if( S > BigFloat(0) && S < BigFloat(1) )
        {
            const BigFloat W = Sqrt(-2*Log(S)/S);
            sample = mean + stddev*U*W;
            break;
        }
    }

    return sample;
}

template<>
Complex<BigFloat> SampleNormal
( const Complex<BigFloat>& mean, const BigFloat& stddev )
{
    // TODO: Consider saving factor of two from Marsiglia's
    Complex<BigFloat> sample;
    BigFloat stddevAdj = stddev;
    stddevAdj /= Sqrt(BigFloat(2));

    sample.real( SampleNormal( mean.real(), stddevAdj ) );
    sample.imag( SampleNormal( mean.imag(), stddevAdj ) );
    return sample;
}
#endif // ifdef EL_HAVE_MPC

} // namespace El
