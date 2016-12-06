/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_LIMITS_HPP
#define EL_LIMITS_HPP

#include <climits> // for INT_MIN et al. due to the Intel C++11 limitations

namespace El {

template<typename Real>
struct IsFixedPrecision
{ static const bool value=true; };
#ifdef EL_HAVE_MPC
template<>
struct IsFixedPrecision<BigFloat>
{ static const bool value=false; };
#endif

template<typename Real>
struct MantissaBits;

template<> struct MantissaBits<unsigned>
{ static const unsigned value = 8*sizeof(unsigned); };
template<> struct MantissaBits<int>
{ static const unsigned value = 8*sizeof(int)-1; };

template<> struct MantissaBits<unsigned long>
{ static const unsigned value = 8*sizeof(unsigned long); };
template<> struct MantissaBits<long int>
{ static const unsigned value = 8*sizeof(long int)-1; };

template<> struct MantissaBits<unsigned long long>
{ static const unsigned value = 8*sizeof(unsigned long long); };
template<> struct MantissaBits<long long int>
{ static const unsigned value = 8*sizeof(long long int)-1; };

template<> struct MantissaBits<float>
{ static const unsigned value = 24; };
template<> struct MantissaBits<double>
{ static const unsigned value = 53; };

#ifdef EL_HAVE_QD
template<> struct MantissaBits<DoubleDouble>
{ static const unsigned value = 106; };
template<> struct MantissaBits<QuadDouble>
{ static const unsigned value = 212; };
#endif

#ifdef EL_HAVE_QUAD
template<> struct MantissaBits<Quad>
{ static const unsigned value = 113; };
#endif

// NOTE: The 'Num' is only prepended to avoid a symbol conflict
template<typename T>
Int NumMantissaBits( const T& alpha=T() )
{ return MantissaBits<T>::value; }
#ifdef EL_HAVE_MPC
template<>
inline Int NumMantissaBits<BigFloat>( const BigFloat& alpha )
{ return alpha.Precision(); }
template<>
inline Int NumMantissaBits<BigInt>( const BigInt& alpha )
{ return alpha.NumBits(); }
#endif

template<typename Real1,typename Real2>
struct MantissaIsLonger
{ static const bool value =
    MantissaBits<Real1>::value >
    MantissaBits<Real2>::value; };
#ifdef EL_HAVE_MPC
// While this isn't necessarily always true, it would be a capitally bad
// idea to use MPFR without using higher than the available fixed precision
template<typename Real2> struct MantissaIsLonger<BigFloat,Real2>
{ static const bool value = true; };
template<> struct MantissaIsLonger<BigFloat,BigFloat>
{ static const bool value = false; };
#endif

namespace limits {

// TODO: Extend this list of functions as needed

// NOTE: These functions take arguments so that types such as BigFloat,
//       which have configurable properties, can be passed a particular
//       instance of the type to extract its limits

template<typename Real,
         typename=EnableIf<IsReal<Real>>>
inline Real Base()
{ return std::numeric_limits<Real>::radix; }
#ifdef EL_HAVE_QD
template<>
inline DoubleDouble Base<DoubleDouble>()
{ return Base<double>(); }
template<>
inline QuadDouble Base<QuadDouble>()
{ return Base<double>(); }
#endif
#ifdef EL_HAVE_QUAD
// std::numeric_limits<__float128>::radix does not seem to be reliably defined
template<>
inline Quad Base<Quad>()
{ return Quad(2); }
#endif
#ifdef EL_HAVE_MPC
template<>
inline BigFloat Base<BigFloat>()
{ return BigFloat(2); }
#endif

// NOTE: There is disagreement by a multiple of the base between
//       Demmel/LAPACK-style 'epsilon' and Higham/STL-style 'epsilon'.
//       We call the former 'epsilon' and the latter 'precision' for
//       consistency with LAPACK.
template<typename Real,
         typename=EnableIf<IsReal<Real>>>
inline Real Precision( const Real& alpha=Real(1) )
{ return std::numeric_limits<Real>::epsilon(); }
template<typename Real,
         typename=EnableIf<IsReal<Real>>>
inline Real Epsilon( const Real& alpha=Real(1) )
{ return Precision<Real>()*std::numeric_limits<Real>::round_error(); }
#ifdef EL_HAVE_QD
template<> inline DoubleDouble Precision( const DoubleDouble& alpha )
{ return dd_real::_eps; }
template<> inline DoubleDouble Epsilon( const DoubleDouble& alpha )
{ return dd_real::_eps/Base<double>(); }

template<> inline QuadDouble Precision( const QuadDouble& alpha )
{ return qd_real::_eps; }
template<> inline QuadDouble Epsilon( const QuadDouble& alpha )
{ return qd_real::_eps/Base<double>(); }
#endif
#ifdef EL_HAVE_QUAD
// NOTE: The FLT128_* macros require support for the -std=gnu++11
//       literal 'Q', which is *NOT* provided by GCC with -std=c++11, but
//       *IS* provided by Intel with -std=c++11.
template<> inline Quad Precision( const Quad& alpha )
{ return FLT128_EPSILON; }
template<> inline Quad Epsilon( const Quad& alpha )
{ return Precision<Quad>()/Base<Quad>(); }
#endif
#ifdef EL_HAVE_MPC
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
#endif

template<typename Real,
         typename=EnableIf<IsReal<Real>>>
struct ExponentTypeHelper
{ typedef Real type; };
#ifdef EL_HAVE_MPC
template<> struct ExponentTypeHelper<BigFloat>
{ typedef mpfr_exp_t type; };
#endif

template<typename Real,
         typename=EnableIf<IsReal<Real>>>
using ExponentType = typename ExponentTypeHelper<Real>::type;

template<typename Real,
         typename=EnableIf<IsReal<Real>>>
inline ExponentType<Real> MaxExponent()
{ return std::numeric_limits<Real>::max_exponent; }
template<typename Real,
         typename=EnableIf<IsReal<Real>>>
inline ExponentType<Real> MinExponent()
{ return std::numeric_limits<Real>::min_exponent; }
#ifdef EL_HAVE_QD
template<>
inline ExponentType<DoubleDouble> MaxExponent<DoubleDouble>()
{ return MaxExponent<double>(); }
template<>
inline ExponentType<DoubleDouble> MinExponent<DoubleDouble>()
{ return MinExponent<double>(); }

template<>
inline ExponentType<QuadDouble> MaxExponent<QuadDouble>()
{ return MaxExponent<double>(); }
template<>
inline ExponentType<QuadDouble> MinExponent<QuadDouble>()
{ return MinExponent<double>(); }
#endif
#ifdef EL_HAVE_QUAD
template<>
inline ExponentType<Quad> MaxExponent<Quad>()
{ return FLT128_MAX_EXP; }
template<>
inline ExponentType<Quad> MinExponent<Quad>()
{ return FLT128_MIN_EXP; }
#endif
#ifdef EL_HAVE_MPC
// NOTE: The minimum and maximum exponents are *not* dependent upon MPFR's
// dynamic precision and are typically fixed at 2^30 - 1
template<>
inline ExponentType<BigFloat>
MaxExponent<BigFloat>()
{ return mpfr_get_emax(); }
template<>
inline ExponentType<BigFloat>
MinExponent<BigFloat>()
{ return mpfr_get_emin(); }
#endif

template<typename Real,
         typename=EnableIf<IsReal<Real>>>
inline Real Max( const Real& alpha=Real(1) )
{ return std::numeric_limits<Real>::max(); }
template<typename Real,
         typename=EnableIf<IsReal<Real>>>
inline Real Min( const Real& alpha=Real(1) )
{ return std::numeric_limits<Real>::min(); }
// NOTE: These implementations do not use std::numeric_limits<Real>::lowest()
//       because said routine does not seem to be provided by recent Intel
//       implementations of C++11
template<typename Real,
         typename=EnableIf<IsReal<Real>>>
inline Real Lowest( const Real& alpha=Real(1) )
{ return -std::numeric_limits<Real>::max(); }
#ifdef EL_USE_64BIT_INTS
template<> inline long long Lowest<long long>( const long long& alpha )
{ return LLONG_MIN; }
#else
template<> inline int Lowest<int>( const int& alpha )
{ return INT_MIN; }
#endif
#ifdef EL_HAVE_QD
template<> inline DoubleDouble Max<DoubleDouble>( const DoubleDouble& alpha )
{ return dd_real::_max; }
template<> inline DoubleDouble Min<DoubleDouble>( const DoubleDouble& alpha )
{ return dd_real::_min_normalized; }
template<> inline DoubleDouble Lowest<DoubleDouble>( const DoubleDouble& alpha )
{ return -dd_real::_max; }

template<> inline QuadDouble Max<QuadDouble>( const QuadDouble& alpha )
{ return qd_real::_max; }
template<> inline QuadDouble Min<QuadDouble>( const QuadDouble& alpha )
{ return qd_real::_min_normalized; }
template<> inline QuadDouble Lowest<QuadDouble>( const QuadDouble& alpha )
{ return -qd_real::_max; }
#endif
#ifdef EL_HAVE_QUAD
template<> inline Quad Max<Quad>( const Quad& alpha )
{ return FLT128_MAX; }
template<> inline Quad Min<Quad>( const Quad& alpha )
{ return FLT128_MIN; }
template<> inline Quad Lowest<Quad>( const Quad& alpha )
{ return -FLT128_MAX; }
#endif
#ifdef EL_HAVE_MPC
template<>
inline BigFloat Max<BigFloat>( const BigFloat& alpha )
{
    BigFloat one(1,alpha.Precision());
    return (one-Epsilon(one)) << MaxExponent<BigFloat>();
}
template<>
inline BigFloat Min<BigFloat>( const BigFloat& alpha )
{
    BigFloat one(1,alpha.Precision());
    return one << (MinExponent<BigFloat>()-1);
}
template<>
inline BigFloat Lowest<BigFloat>( const BigFloat& alpha )
{ return -Max<BigFloat>(alpha); }
#endif

template<typename Real,
         typename=EnableIf<IsReal<Real>>>
inline Real SafeMin( const Real& alpha=Real(1) )
{
    const Real one = Real(1);
    const Real eps = Epsilon<Real>();
    const Real minVal = Min<Real>();
    const Real invMax = one/Max<Real>();
    static const Real safeMin =
      ( invMax>minVal ? invMax*(one+eps) : minVal );
    return safeMin;
}
#ifdef EL_HAVE_MPC
template<>
inline BigFloat SafeMin( const BigFloat& alpha )
{
    // TODO: Decide how to only recompute this when the precision changes
    const BigFloat one = BigFloat(1,alpha.Precision());
    const BigFloat eps = Epsilon(alpha);
    const BigFloat minVal = Min(alpha);
    const BigFloat invMax = one/Max(alpha);
    return ( invMax>minVal ? invMax*(one+eps) : minVal );
}
#endif

// For now, this routine is only used within Givens in order to provide
// a comfortable upper bound on positive numbers that can be safely squared
// (and then safely added and substracted)
//
// CITATION
//
// Please see LAPACK Working Note 148,
//
//   D. Bindel, J. Demmel, W. Kahan, and O. Marques,
//   "On computing Givens rotations reliably and efficiently",
//   http://www.netlib.org/lapack/lawnspdf/lawn148.pdf
//
// which resulted in the LAPACK routines {s,d,c,z}lartg. But note
// that the LAPACK implementations slightly differ from said working
// note in that they round z^2 to the nearest radix rather than z
// (which results in a different result with float, but the same with
// double).
//
template<typename Real,
         typename=EnableIf<IsReal<Real>>>
inline Real SafeMinToSquare( const Real& alpha=Real(1) )
{
    const Real base = Base<Real>();
    const Real eps = Epsilon<Real>();
    const Real safeMin = SafeMin<Real>();
    static const Real safeMinToSquare =
      Pow( base, Round((Log(safeMin/eps)/Log(base))/Real(2)) );
    return safeMinToSquare;
}
#ifdef EL_HAVE_MPC
template<>
inline BigFloat SafeMinToSquare( const BigFloat& alpha )
{
    // TODO: Decide how to only recompute this when the precision changes
    const BigFloat two(2);
    const BigFloat safeMinToSquare =
      Pow( two, Round((Log2(SafeMin(alpha)/Epsilon(alpha)))/two) );
    return safeMinToSquare;
}
#endif

// At the moment, this is only used within the subroutine CubicSecular within
// the secular equation solver for the bidiagonal SVD. See LAPACK's
// {s,d}laed6 [CITATION] for the motivation for these constants.
template<typename Real,
         typename=EnableIf<IsReal<Real>>>
inline Real SafeMinToCube( const Real& alpha=Real(1) )
{
    const Real base = Base<Real>();
    const Real eps = Epsilon<Real>();
    const Real safeMin = SafeMin<Real>();
    static const Real safeMinToCube =
      Pow( base, Round((Log(safeMin/eps)/Log(base))/Real(3)) );
    return safeMinToCube;
}
#ifdef EL_HAVE_MPC
template<>
inline BigFloat SafeMinToCube( const BigFloat& alpha )
{
    // TODO: Decide how to only recompute this when the precision changes
    const BigFloat base(2);
    const BigFloat safeMinToSquare =
      Pow( base, Round((Log2(SafeMin(alpha)/Epsilon(alpha)))/3) );
    return safeMinToSquare;
}
#endif

template<typename Real,
         typename=EnableIf<IsReal<Real>>>
inline Real Infinity( const Real& alpha=Real(1) )
{ return std::numeric_limits<Real>::infinity(); }
#ifdef EL_HAVE_QD
template<> inline DoubleDouble Infinity<DoubleDouble>
( const DoubleDouble& alpha )
{ return dd_real::_inf; }
template<> inline QuadDouble Infinity<QuadDouble>
( const QuadDouble& alpha )
{ return qd_real::_inf; }
#endif
#ifdef EL_HAVE_QUAD
template<> inline Quad Infinity<Quad>( const Quad& alpha )
{
    // libquadmath does not document how to return infinity, so, for now,
    // we will instead return twice the maximum floating-point number
    return 2*Max<Quad>(alpha);
}
#endif
#ifdef EL_HAVE_MPC
template<>
inline BigFloat Infinity<BigFloat>( const BigFloat& alpha )
{
    BigFloat inf(0,alpha.Precision());
    mpfr_set_inf( inf.Pointer(), 1 );
    return inf;
}
#endif

template<typename Real,
         typename=EnableIf<IsReal<Real>>>
inline bool IsFinite( const Real& alpha )
{ return std::isfinite(alpha); }
#ifdef EL_HAVE_QD
template<>
inline bool IsFinite( const DoubleDouble& alpha )
{ return alpha.isfinite(); }
template<>
inline bool IsFinite( const QuadDouble& alpha )
{ return alpha.isfinite(); }
#endif
#ifdef EL_HAVE_QUAD
template<>
inline bool IsFinite( const Quad& alpha )
{ return finiteq(alpha) != 0; }
#endif
#ifdef EL_HAVE_MPC
template<>
inline bool IsFinite( const BigFloat& alpha )
{ return mpfr_number_p( alpha.LockedPointer() ) != 0; }
#endif

} // namespace limits

inline Int BinaryToDecimalPrecision( Int prec )
{ return Int(Floor(prec*std::log10(2.))); }

} // namespace El

#endif // ifndef EL_LIMITS_HPP
