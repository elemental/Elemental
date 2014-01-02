/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental-lite.hpp"
using namespace elem;

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );

    try
    { LogicError("Phone number: ",4+4,3*2,7,"-",5,3,0,3*3); }
    catch( std::exception& e ) { ReportException(e); }

    Finalize();
    return 0;
}
