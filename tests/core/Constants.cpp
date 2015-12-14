/*
   Copyright (c) 2009-2015, Jack Poulson
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
        const auto safeInv = SafeMin<Real>()/Epsilon<Real>();
        Output( title );
        Output("  safeMin:   ",SafeMin<Real>());
        Output("  epsilon:   ",Epsilon<Real>());
        Output("  precision: ",Precision<Real>());
        Output("  safeInv:   ",safeInv);
        Output("  min:       ",Min<Real>());
        Output("  max:       ",Max<Real>());
        Output("  lowest:    ",Lowest<Real>());
        Output("  infinity:  ",Infinity<Real>());
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
#ifdef EL_HAVE_QUAD
        QueryLimits<Quad>( "Quad-precision:" );
#endif
#ifdef EL_HAVE_MPC
        QueryLimits<BigFloat>( "BigFloat:" );
#endif
    }

    return 0;
}
