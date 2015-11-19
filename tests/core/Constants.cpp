/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
#include <float.h>
using namespace El;

int 
main( int argc, char* argv[] )
{
    Environment env( argc, argv );
    const Int commRank = mpi::Rank();

    if( commRank == 0 )
    {
        const float safeMinFloat = lapack::MachineSafeMin<float>();
        const float epsilonFloat = lapack::MachineEpsilon<float>();
        const float safeInvFloat = safeMinFloat/epsilonFloat;
        const double safeMinDouble = lapack::MachineSafeMin<double>();
        const double epsilonDouble = lapack::MachineEpsilon<double>();
        const double safeInvDouble = safeMinDouble/epsilonDouble;
        std::cout << "Single precision:\n"
                  << "  safeMin: " << safeMinFloat << "\n"
                  << "  epsilon: " << epsilonFloat << "\n"
                  << "  safeInv: " << safeInvFloat << "\n"
                  << "  FLT_MIN: " << FLT_MIN << "\n"
                  << "  FLT_MAX: " << FLT_MAX << "\n"
                  << "\n"
                  << "Double precision:\n"
                  << "  safeMin: " << safeMinDouble << "\n"
                  << "  epsilon: " << epsilonDouble << "\n"
                  << "  safeInv: " << safeInvDouble << "\n"
                  << "  DBL_MIN: " << DBL_MIN << "\n"
                  << "  DBL_MAX: " << DBL_MAX << "\n"
                  << std::endl;
    }

    return 0;
}
