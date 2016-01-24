/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

template<typename Real>
void QueryLimits( const std::string& title )
{
        const auto safeInv = limits::SafeMin<Real>()/limits::Epsilon<Real>();
        Output( title );
        Output("  safeMin:   ",limits::SafeMin<Real>());
        Output("  epsilon:   ",limits::Epsilon<Real>());
        Output("  precision: ",limits::Precision<Real>());
        Output("  safeInv:   ",safeInv);
        Output("  min:       ",limits::Min<Real>());
        Output("  max:       ",limits::Max<Real>());
        Output("  lowest:    ",limits::Lowest<Real>());
        Output("  infinity:  ",limits::Infinity<Real>());
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
        mpc::SetPrecision( 64 );
        QueryLimits<BigFloat>( "BigFloat (64):" );
        mpc::SetPrecision( 128 );
        QueryLimits<BigFloat>( "BigFloat (128):" );
        mpc::SetPrecision( 256 );
        QueryLimits<BigFloat>( "BigFloat (256):" );
        mpc::SetPrecision( 512 );
        QueryLimits<BigFloat>( "BigFloat (512):" );
        mpc::SetPrecision( 1024 );
        QueryLimits<BigFloat>( "BigFloat (1024):" );
#endif
    }

    return 0;
}
