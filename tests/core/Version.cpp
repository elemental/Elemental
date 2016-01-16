/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace El;

int 
main( int argc, char* argv[] )
{
    Environment env( argc, argv );
    const Int commRank = mpi::Rank();

    if( commRank == 0 )
    {
        PrintVersion();
        PrintCCompilerInfo();
        PrintCxxCompilerInfo();
    }

    return 0;
}
