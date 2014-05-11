/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_TYPES_DECL_HPP
#define ELEM_TYPES_DECL_HPP

namespace elem {

typedef unsigned char byte;

// If these are changes, you must make sure that they have 
// existing MPI datatypes. This is only sometimes true for 'long long'
#ifdef ELEM_USE_64BIT_INTS
typedef long long int Int;
typedef long long unsigned Unsigned;
#else
typedef int Int;
typedef unsigned Unsigned;
#endif
 
typedef Complex<float>  scomplex; 
typedef Complex<double> dcomplex;

template<typename Real>
struct ValueInt
{
    Real value;
    Int index;

    static bool Lesser( const ValueInt<Real>& a, const ValueInt<Real>& b )
    { return a.value < b.value; }
    static bool Greater( const ValueInt<Real>& a, const ValueInt<Real>& b )
    { return a.value > b.value; }
};

template<typename Real>
struct ValueInt<Complex<Real>>
{
    Complex<Real> value;
    Int index;

    static bool Lesser( const ValueInt<Real>& a, const ValueInt<Real>& b )
    { return Abs(a.value) < Abs(b.value); }
    static bool Greater( const ValueInt<Real>& a, const ValueInt<Real>& b )
    { return Abs(a.value) > Abs(b.value); }
};

template<typename Real>
struct ValueIntPair
{
    Real value;
    Int indices[2];
    
    static bool Lesser( const ValueInt<Real>& a, const ValueInt<Real>& b )
    { return a.value < b.value; }
    static bool Greater( const ValueInt<Real>& a, const ValueInt<Real>& b )
    { return a.value > b.value; }
};

template<typename Real>
struct ValueIntPair<Complex<Real>>
{
    Complex<Real> value;
    Int indices[2];
    
    static bool Lesser
    ( const ValueIntPair<Real>& a, const ValueIntPair<Real>& b )
    { return Abs(a.value) < Abs(b.value); }
    static bool Greater
    ( const ValueIntPair<Real>& a, const ValueIntPair<Real>& b )
    { return Abs(a.value) > Abs(b.value); }
};

// For the safe computation of products. The result is given by 
//   product = rho * exp(kappa*n)
// where rho lies in (usually on) the unit circle and kappa is real-valued.
template<typename F>
struct SafeProduct
{
    F rho;
    Base<F> kappa;
    Int n;

    SafeProduct( Int numEntries );
};

// The basic eigenvalue structure of a Hermitian matrix
struct InertiaType
{
    Int numPositive, numNegative, numZero;
};

// At some point this enum will be a member of each matrix class and will
// be used to simplify the syntax of various operations
namespace MatrixClassNS {
enum MatrixClass
{
    UNSPECIFIED,
    GENERAL,
    HERMITIAN,
    HERMITIAN_LOWER,
    HERMITIAN_UPPER, 
    SYMMETRIC,
    SYMMETRIC_LOWER,
    SYMMETRIC_UPPER,
    SKEW_SYMMETRIC,
    SKEW_SYMMETRIC_LOWER,
    SKEW_SYMMETRIC_UPPER,
    UNITARY,
    TRIANGULAR_LOWER,
    TRIANGULAR_LOWER_UNIT,
    TRIANGULAR_UPPER,
    TRIANGULAR_UPPER_UNIT,
    HESSENBERG_LOWER,
    HESSENBERG_UPPER,
    PERMUTATION,
    PERMUTATION_VECTOR,
    PIVOT_SEQUENCE,
    // Packed factorizations/decompositions
    LU_PACKED,
    QR_PACKED,
    RQ_PACKED,
    LQ_PACKED,
    QL_PACKED,
    BIDIAG_PACKED,
    TRIDIAG_LOWER_PACKED,
    TRIDIAG_UPPER_PACKED,
    HESSENBERG_LOWER_PACKED,
    HESSENBERG_UPPER_PACKED
};
}
using namespace MatrixClassNS;

namespace ConjugationNS {
enum Conjugation
{
    UNCONJUGATED,
    CONJUGATED
};
}
using namespace ConjugationNS;

namespace DistNS {
enum Dist
{
    MC,   // Col of a matrix distribution
    MD,   // Diagonal of a matrix distribution
    MR,   // Row of a matrix distribution
    VC,   // Col-major vector distribution
    VR,   // Row-major vector distribution
    STAR, // Give to every process
    CIRC  // Give to a single process
};
std::string DistToString( Dist distribution );
Dist StringToDist( std::string s );
}
using namespace DistNS;
typedef Dist Distribution;

template<Dist U,Dist V>
constexpr Dist DiagColDist() { return ( U==STAR ? V : U ); }
template<Dist U,Dist V>
constexpr Dist DiagRowDist() { return ( U==STAR ? U : V ); }

template<> constexpr Dist DiagColDist<MC,MR>() { return MD; }
template<> constexpr Dist DiagRowDist<MC,MR>() { return STAR; }
template<> constexpr Dist DiagColDist<MR,MC>() { return MD; }
template<> constexpr Dist DiagRowDist<MR,MC>() { return STAR; }

template<Dist U,Dist V>
constexpr Dist DiagInvColDist() { return ( U==STAR ? V : U ); }
template<Dist U,Dist V>
constexpr Dist DiagInvRowDist() { return ( U==STAR ? U : V ); }

template<> constexpr Dist DiagInvColDist<MD,STAR>() { return MC; }
template<> constexpr Dist DiagInvRowDist<MD,STAR>() { return MR; }
template<> constexpr Dist DiagInvColDist<STAR,MD>() { return MC; }
template<> constexpr Dist DiagInvRowDist<STAR,MD>() { return MR; }

template<Dist U> 
constexpr Dist GatheredDist() { return ( U==CIRC ? CIRC : STAR ); }

template<Dist U> constexpr Dist PartialDist() { return U; }
template<> constexpr Dist PartialDist<VC>() { return MC; }
template<> constexpr Dist PartialDist<VR>() { return MR; }

template<Dist U,Dist V> constexpr Dist ScatteredRowDist() { return V; }
template<> constexpr Dist ScatteredRowDist<VC,STAR>() { return MR; }
template<> constexpr Dist ScatteredRowDist<VR,STAR>() { return MC; }

template<Dist U,Dist V> constexpr Dist ScatteredColDist() 
{ return ScatteredRowDist<V,U>(); }

namespace ViewTypeNS {
enum ViewType
{
    OWNER = 0x0,
    VIEW = 0x1,
    OWNER_FIXED = 0x2,
    VIEW_FIXED = 0x3,
    LOCKED_OWNER = 0x4, // unused
    LOCKED_VIEW = 0x5,
    LOCKED_OWNER_FIXED = 0x6, // unused
    LOCKED_VIEW_FIXED = 0x7
};
static inline bool IsViewing( ViewType v )
{ return ( v & VIEW ) != 0; }
static inline bool IsFixedSize( ViewType v )
{ return ( v & OWNER_FIXED ) != 0; }
static inline bool IsLocked( ViewType v )
{ return ( v & LOCKED_OWNER ) != 0; }
}
using namespace ViewTypeNS;

namespace ForwardOrBackwardNS {
enum ForwardOrBackward
{
    FORWARD,
    BACKWARD
};
}
using namespace ForwardOrBackwardNS;

namespace GridOrderNS {
enum GridOrder
{
    ROW_MAJOR,
    COLUMN_MAJOR
};
}
using namespace GridOrderNS;

namespace LeftOrRightNS {
enum LeftOrRight
{
    LEFT,
    RIGHT
};
char LeftOrRightToChar( LeftOrRight side );
LeftOrRight CharToLeftOrRight( char c );
}
using namespace LeftOrRightNS;

namespace SortTypeNS {
enum SortType
{
    UNSORTED,
    DESCENDING,
    ASCENDING
};
}
using namespace SortTypeNS;

namespace NormTypeNS {
enum NormType
{
    ONE_NORM,           // Operator one norm
    INFINITY_NORM,      // Operator infinity norm
    ENTRYWISE_ONE_NORM, // One-norm of vectorized matrix
    MAX_NORM,           // Maximum entry-wise magnitude
    NUCLEAR_NORM,       // One-norm of the singular values
    FROBENIUS_NORM,     // Two-norm of the singular values
    TWO_NORM            // Infinity-norm of the singular values
};
}
using namespace NormTypeNS;

namespace OrientationNS {
enum Orientation
{
    NORMAL,
    TRANSPOSE,
    ADJOINT
};
char OrientationToChar( Orientation orientation );
Orientation CharToOrientation( char c );
}
using namespace OrientationNS;

namespace UnitOrNonUnitNS {
enum UnitOrNonUnit
{
    NON_UNIT,
    UNIT
};
char UnitOrNonUnitToChar( UnitOrNonUnit diag );
UnitOrNonUnit CharToUnitOrNonUnit( char c );
}
using namespace UnitOrNonUnitNS;

namespace UpperOrLowerNS {
enum UpperOrLower
{
    LOWER,
    UPPER
};
char UpperOrLowerToChar( UpperOrLower uplo );
UpperOrLower CharToUpperOrLower( char c );
}
using namespace UpperOrLowerNS;

namespace VerticalOrHorizontalNS {
enum VerticalOrHorizontal
{
    VERTICAL,
    HORIZONTAL
};
}
using namespace VerticalOrHorizontalNS;

} // namespace elem

#endif // ifndef ELEM_TYPES_DECL_HPP
