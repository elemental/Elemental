/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_LIMITS_HPP
#define EL_LIMITS_HPP

namespace El {
namespace limits {

// TODO: Extend this list of functions as needed

// NOTE: These functions take arguments so that types such as BigFloat,
//       which have configurable properties, can be passed a particular 
//       instance of the type to extract its limits

template<typename Real,typename=EnableIf<IsReal<Real>>>
struct ExponentTypeHelper
{ typedef Real type; };
#ifdef EL_HAVE_MPC
template<> struct ExponentTypeHelper<BigFloat>
{ typedef mpfr_exp_t type; };
#endif

template<typename Real,typename=EnableIf<IsReal<Real>>>
using ExponentType = typename ExponentTypeHelper<Real>::type;

template<typename Real,typename=EnableIf<IsReal<Real>>>
inline ExponentType<Real> MaxExponent( const Real& alpha=Real(1) )
{ return std::numeric_limits<Real>::max_exponent(); }

template<typename Real,typename=EnableIf<IsReal<Real>>>
inline ExponentType<Real> MinExponent( const Real& alpha=Real(1) )
{ return std::numeric_limits<Real>::min_exponent(); }

template<typename Real,typename=EnableIf<IsReal<Real>>>
inline Real Max( const Real& alpha=Real(1) )
{ return std::numeric_limits<Real>::max(); }
template<typename Real,typename=EnableIf<IsReal<Real>>>
inline Real Min( const Real& alpha=Real(1) )
{ return std::numeric_limits<Real>::min(); }
template<typename Real,typename=EnableIf<IsReal<Real>>>
inline Real Lowest( const Real& alpha=Real(1) )
{ return std::numeric_limits<Real>::lowest(); }

// NOTE: There is disagreement of a factor of two between
//       Demmel/LAPACK-style 'epsilon' and Higham/STL-style 'epsilon'.
//       We call the former 'epsilon' and the latter 'precision' for
//       consistency with LAPACK.
template<typename Real,typename=EnableIf<IsReal<Real>>>
inline Real Precision( const Real& alpha=Real(1) )
{ return std::numeric_limits<Real>::epsilon(); }
template<typename Real,typename=EnableIf<IsReal<Real>>>
inline Real Epsilon( const Real& alpha=Real(1) )
{ return Precision<Real>()*std::numeric_limits<Real>::round_error(); }

template<typename Real,typename=EnableIf<IsReal<Real>>>
inline Real SafeMin( const Real& alpha=Real(1) )
{
    const Real one = Real(1);
    const Real eps = Epsilon<Real>();
    const Real minVal = Min<Real>();
    const Real invMax = one/Max<Real>();
    return ( invMax>minVal ? invMax*(one+eps) : minVal );
}

template<typename Real,typename=EnableIf<IsReal<Real>>>
inline Real Infinity( const Real& alpha=Real(1) )
{ return std::numeric_limits<Real>::infinity(); }

template<typename Real,typename=EnableIf<IsReal<Real>>>
inline bool IsFinite( const Real& alpha )
{ return std::isfinite(alpha); }

#ifdef EL_HAVE_QUAD
template<>
inline ExponentType<Quad> MaxExponent<Quad>( const Quad& alpha )
{ return FLT128_MAX_EXP; }
template<>
inline ExponentType<Quad> MinExponent<Quad>( const Quad& alpha )
{ return FLT128_MIN_EXP; }

// NOTE: The following few macros require support for the -std=gnu++11
//       literal 'Q', which is *NOT* provided by GCC with -std=c++11, but 
//       *IS* provided by Intel with -std=c++11.
template<> inline Quad Max<Quad>( const Quad& alpha )
{ return FLT128_MAX; }
template<> inline Quad Min<Quad>( const Quad& alpha )
{ return FLT128_MIN; }
template<> inline Quad Lowest<Quad>( const Quad& alpha )
{ return -FLT128_MAX; }

template<> inline Quad Precision( const Quad& alpha )
{ return FLT128_EPSILON; }
template<> inline Quad Epsilon( const Quad& alpha )
{ return Precision<Quad>()/Quad(2); }

template<> inline Quad Infinity<Quad>( const Quad& alpha )
{
    // libquadmath does not document how to return infinity, so, for now,
    // we will instead return the maximum floating-point number
    return Max<Quad>(alpha);
}

template<>
inline bool IsFinite( const Quad& alpha )
{ return finiteq(alpha) != 0; }

#endif // ifdef EL_HAVE_QUAD

#ifdef EL_HAVE_MPC
template<>
inline ExponentType<BigFloat>
MaxExponent<BigFloat>( const BigFloat& alpha )
{ return mpfr_get_emax(); }

template<>
inline ExponentType<BigFloat>
MinExponent<BigFloat>( const BigFloat& alpha )
{ return mpfr_get_emin(); }

template<>
inline BigFloat Precision<BigFloat>( const BigFloat& alpha )
{ 
    // NOTE: This 'precision' is the number of bits in the mantissa
    auto p = alpha.Precision();

    // precision = b^{-(p-1)} = 2^{-(p-1)}
    BigFloat one(1,p);
    return one >> (p-1);
}

template<>
inline BigFloat Epsilon<BigFloat>( const BigFloat& alpha )
{ 
    // NOTE: This 'precision' is the number of bits in the mantissa
    auto p = alpha.Precision();

    // epsilon = b^{-(p-1)} / 2 = 2^{-p}
    BigFloat one(1,p);
    return one >> p;
}

template<>
inline BigFloat Max<BigFloat>( const BigFloat& alpha )
{
    BigFloat one(1,alpha.Precision());
    return (one-Epsilon(one)) << MaxExponent(one);
}

template<>
inline BigFloat Min<BigFloat>( const BigFloat& alpha )
{
    BigFloat one(1,alpha.Precision());
    return one << (MinExponent(one)-1);
}

template<>
inline BigFloat Lowest<BigFloat>( const BigFloat& alpha )
{ return -Max<BigFloat>(alpha); }

template<>
inline BigFloat Infinity<BigFloat>( const BigFloat& alpha )
{
    BigFloat inf(alpha.Precision());
    mpfr_set_inf( inf.Pointer(), 1 );
    return inf;
}

template<>
inline bool IsFinite( const BigFloat& alpha )
{ return mpfr_number_p( alpha.LockedPointer() ) != 0; }

#endif // ifdef EL_HAVE_MPC

} // namespace limits
} // namespace El

#endif // ifndef EL_LIMITS_HPP
