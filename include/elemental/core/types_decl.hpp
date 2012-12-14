/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {

typedef unsigned char byte;
 
typedef Complex<float>  scomplex; 
typedef Complex<double> dcomplex;

// For the safe computation of products. The result is given by 
//   product = rho * exp(kappa*n)
// where rho lies in (usually on) the unit circle and kappa is real-valued.
template<typename F,typename Int=int>
struct SafeProduct
{
    F rho;
    typename Base<F>::type kappa;
    Int n;

    SafeProduct( Int numEntries );
};

namespace conjugation_wrapper {
enum Conjugation
{
    UNCONJUGATED,
    CONJUGATED
};
}
using namespace conjugation_wrapper;

namespace distribution_wrapper {
enum Distribution
{
    MC,  // Col of a matrix distribution
    MD,  // Diagonal of a matrix distribution
    MR,  // Row of a matrix distribution
    VC,  // Col-major vector distribution
    VR,  // Row-major vector distribution
    STAR // Do not distribute
};
std::string DistToString( Distribution distribution );
Distribution StringToDist( std::string s );
}
using namespace distribution_wrapper;

namespace forward_or_backward_wrapper {
enum ForwardOrBackward
{
    FORWARD,
    BACKWARD
};
}
using namespace forward_or_backward_wrapper;

namespace grid_order_wrapper {
enum GridOrder
{
    ROW_MAJOR,
    COLUMN_MAJOR
};
}
using namespace grid_order_wrapper;

namespace left_or_right_wrapper {
enum LeftOrRight
{
    LEFT,
    RIGHT
};
char LeftOrRightToChar( LeftOrRight side );
LeftOrRight CharToLeftOrRight( char c );
}
using namespace left_or_right_wrapper;

namespace norm_type_wrapper {
enum NormType
{
    ONE_NORM,       // Operator one norm
    INFINITY_NORM,  // Operator infinity norm
    MAX_NORM,       // Maximum entry-wise magnitude
    NUCLEAR_NORM,   // One-norm of the singular values
    FROBENIUS_NORM, // Two-norm of the singular values
    TWO_NORM        // Infinity-norm of the singular values
};
}
using namespace norm_type_wrapper;

namespace orientation_wrapper {
enum Orientation
{
    NORMAL,
    TRANSPOSE,
    ADJOINT
};
char OrientationToChar( Orientation orientation );
Orientation CharToOrientation( char c );
}
using namespace orientation_wrapper;

namespace unit_or_non_unit_wrapper {
enum UnitOrNonUnit
{
    NON_UNIT,
    UNIT
};
char UnitOrNonUnitToChar( UnitOrNonUnit diag );
UnitOrNonUnit CharToUnitOrNonUnit( char c );
}
using namespace unit_or_non_unit_wrapper;

namespace upper_or_lower_wrapper {
enum UpperOrLower
{
    LOWER,
    UPPER
};
char UpperOrLowerToChar( UpperOrLower uplo );
UpperOrLower CharToUpperOrLower( char c );
}
using namespace upper_or_lower_wrapper;

namespace vertical_or_horizontal_wrapper {
enum VerticalOrHorizontal
{
    VERTICAL,
    HORIZONTAL
};
}
using namespace vertical_or_horizontal_wrapper;

} // namespace elem
