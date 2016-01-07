/*
   Copyright (c) 2009-2015, Jack Poulson
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

#ifdef EL_HAVE_QUAD
template<>
Quad SampleUniform( Quad a, Quad b )
{
    Quad sample;

#ifdef EL_HAVE_CXX11RANDOM
    std::mt19937& gen = Generator();
    std::uniform_real_distribution<long double> 
      uni((long double)a,(long double)b);
    sample = uni(gen);
#else
    sample = (Real(rand())/(Real(RAND_MAX)+1))*(b-a) + a;
#endif

    return sample;
}

template<>
Complex<Quad> SampleUniform( Complex<Quad> a, Complex<Quad> b )
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
    Real aReal = RealPart(a);
    Real aImag = ImagPart(a);
    Real bReal = RealPart(b);
    Real bImag = ImagPart(b);
    Real realPart = (Real(rand())/(Real(RAND_MAX)+1))*(bReal-aReal) + aReal;
    SetRealPart( sample, realPart );

    Real imagPart = (Real(rand())/(Real(RAND_MAX)+1))*(bImag-aImag) + aImag;
    SetImagPart( sample, imagPart );
#endif

    return sample;
}
#endif // ifdef EL_HAVE_QUAD

#ifdef EL_HAVE_MPC
template<>
BigFloat SampleUniform( BigFloat a, BigFloat b )
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
Int SampleUniform<Int>( Int a, Int b )
{
#ifdef EL_HAVE_CXX11RANDOM
    std::mt19937& gen = Generator();
    std::uniform_int_distribution<Int> intDist(a,b-1); 
    return intDist(gen);
#else
    return a + (rand() % (b-a));
#endif
}

#ifdef EL_HAVE_QUAD
template<>
Quad SampleNormal( Quad mean, Quad stddev )
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
Complex<Quad> SampleNormal( Complex<Quad> mean, Quad stddev )
{
    Complex<Quad> sample;
    stddev = stddev / Sqrt(Quad(2));

#ifdef EL_HAVE_CXX11RANDOM
    std::mt19937& gen = Generator();

    std::normal_distribution<long double> 
      realNormal( (long double)RealPart(mean), (long double)stddev );
    SetRealPart( sample, Quad(realNormal(gen)) );

    std::normal_distribution<long double>
      imagNormal( (long double)ImagPart(mean), (long double)stddev );
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
            SetRealPart( sample, RealPart(mean) + stddev*U*W );
            SetImagPart( sample, ImagPart(mean) + stddev*V*W );
            break;
        }
    }
#endif

    return sample;
}
#endif // ifdef EL_HAVE_QUAD

#ifdef EL_HAVE_MPC
template<>
BigFloat SampleNormal( BigFloat mean, BigFloat stddev )
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

template<>
float SampleBall<float>( float center, float radius )
{ return SampleUniform<float>(center-radius,center+radius); }

template<>
double SampleBall<double>( double center, double radius )
{ return SampleUniform<double>(center-radius,center+radius); }

template<>
Complex<float> SampleBall<Complex<float>>( Complex<float> center, float radius )
{
    const float r = SampleUniform<float>(0,radius);
    const float angle = SampleUniform<float>(0.f,float(2*Pi<float>()));
    return center + Complex<float>(r*cos(angle),r*sin(angle));
}

template<>
Complex<double>
SampleBall<Complex<double>>( Complex<double> center, double radius )
{
    const double r = SampleUniform<double>(0,radius);
    const double angle = SampleUniform<double>(0.,2*Pi<double>());
    return center + Complex<double>(r*cos(angle),r*sin(angle));
}

#ifdef EL_HAVE_QUAD
template<>
Quad SampleBall<Quad>( Quad center, Quad radius )
{ return SampleUniform<Quad>(center-radius,center+radius); }

template<>
Complex<Quad> SampleBall<Complex<Quad>>( Complex<Quad> center, Quad radius )
{
    const Quad r = SampleUniform<Quad>(0,radius);
    const Quad angle = SampleUniform<Quad>(0.f,Quad(2*Pi<Quad>()));
    return center + Complex<Quad>(r*Cos(angle),r*Sin(angle));
}
#endif

#ifdef EL_HAVE_MPC
template<>
BigFloat SampleBall<BigFloat>( BigFloat center, BigFloat radius )
{ return SampleUniform<BigFloat>(center-radius,center+radius); }
#endif

// I'm not certain if there is any good way to define this
template<>
Int SampleBall<Int>( Int center, Int radius )
{
    const double u = SampleBall<double>( center, radius );
    return std::lround(u);
}

} // namespace El
