/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
using namespace El;

template<typename Real>
void QueryLimits( const std::string& title )
{
    const auto safeInv = limits::SafeMin<Real>()/limits::Epsilon<Real>();
    Output( title );
    Output("  base:      ",limits::Base<Real>());
    Output("  epsilon:   ",limits::Epsilon<Real>());
    Output("  precision: ",limits::Precision<Real>());
    Output("  safeMin:   ",limits::SafeMin<Real>());
    Output("  safeInv:   ",safeInv);
    Output("  min:       ",limits::Min<Real>());
    Output("  min/2:     ",limits::Min<Real>()/Real(2));
    Output("  max:       ",limits::Max<Real>());
    Output("  max*2:     ",limits::Max<Real>()*Real(2));
    Output("  lowest:    ",limits::Lowest<Real>());
    Output("  infinity:  ",limits::Infinity<Real>());
    Output("");
}

int 
main( int argc, char* argv[] )
{
    Environment env( argc, argv );
    const Int commRank = mpi::Rank();

    if( commRank == 0 )
    {
        QueryLimits<float>( "Single-precision:" );
        QueryLimits<double>( "Double-precision:" );
#ifdef EL_HAVE_QD
        QueryLimits<DoubleDouble>( "DoubleDouble:" );
        QueryLimits<QuadDouble>( "QuadDouble:" );
#endif
#ifdef EL_HAVE_QUAD
        QueryLimits<Quad>( "Quad-precision:" );
#endif
#ifdef EL_HAVE_MPC
        QueryLimits<BigFloat>( "BigFloat (Default):" );
        mpfr::SetPrecision( 64 );
        QueryLimits<BigFloat>( "BigFloat (64):" );
        mpfr::SetPrecision( 128 );
        QueryLimits<BigFloat>( "BigFloat (128):" );
        mpfr::SetPrecision( 256 );
        QueryLimits<BigFloat>( "BigFloat (256):" );
        mpfr::SetPrecision( 512 );
        QueryLimits<BigFloat>( "BigFloat (512):" );
        mpfr::SetPrecision( 1024 );
        QueryLimits<BigFloat>( "BigFloat (1024):" );
#endif
    }

    return 0;
}
