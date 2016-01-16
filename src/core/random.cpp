/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

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

#ifdef EL_HAVE_CXX11RANDOM
    std::mt19937& gen = Generator();
    std::uniform_real_distribution<long double> 
      realUni((long double)RealPart(a),(long double)RealPart(b));
    SetRealPart( sample, Quad(realUni(gen)) ); 

    std::uniform_real_distribution<long double> 
      imagUni((long double)ImagPart(a),(long double)ImagPart(b));
    SetImagPart( sample, Quad(imagUni(gen)) );
#else
    Quad aReal = RealPart(a);
    Quad aImag = ImagPart(a);
    Quad bReal = RealPart(b);
    Quad bImag = ImagPart(b);
    Quad realPart = (Quad(rand())/(Quad(RAND_MAX)+1))*(bReal-aReal) + aReal;
    SetRealPart( sample, realPart );

    Quad imagPart = (Quad(rand())/(Quad(RAND_MAX)+1))*(bImag-aImag) + aImag;
    SetImagPart( sample, imagPart );
#endif

    return sample;
}
#endif // ifdef EL_HAVE_QUAD

#ifdef EL_HAVE_MPC
template<>
BigInt SampleUniform( const BigInt& a, const BigInt& b )
{
    BigInt sample; 
    gmp_randstate_t randState;
    mpc::RandomState( randState );
    
    mpz_urandomb( sample.Pointer(), randState, b.NumBits() );
    return a+Mod(sample,b-a);
}

template<>
BigFloat SampleUniform( const BigFloat& a, const BigFloat& b )
{
    BigFloat sample; 
    gmp_randstate_t randState;
    mpc::RandomState( randState );

    while( 1 )
    {
        const int ret = mpfr_urandomb( sample.Pointer(), randState );
        if( ret == 0 )
            break;
    }
    return a + sample*(b-a);
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

#ifdef EL_HAVE_CXX11RANDOM
    std::mt19937& gen = Generator();

    std::normal_distribution<long double> 
      realNormal( (long double)RealPart(mean), (long double)stddevAdj );
    SetRealPart( sample, Quad(realNormal(gen)) );

    std::normal_distribution<long double>
      imagNormal( (long double)ImagPart(mean), (long double)stddevAdj );
    SetImagPart( sample, Quad(imagNormal(gen)) );
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
            SetRealPart( sample, RealPart(mean) + stddevAdj*U*W );
            SetImagPart( sample, ImagPart(mean) + stddevAdj*V*W );
            break;
        }
    }
#endif

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
#endif // ifdef EL_HAVE_MPC

} // namespace El
