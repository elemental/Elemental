/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"
#include "El-C.h"
using namespace El;

#define RCG(gridHandle) reinterpret_cast<Grid*>(gridHandle)
#define RCG_const(gridHandle) reinterpret_cast<const Grid*>(gridHandle)

#define CATCH catch( std::exception& e ) { ReportException(e); }

extern "C" {

int ElGridRow( const ElGrid* gridHandle )
{
    int row = -1;
    try { row = RCG_const(gridHandle)->Row(); }
    CATCH
    return row;
}

int ElGridCol( const ElGrid* gridHandle )
{
    int col = -1;
    try { col = RCG_const(gridHandle)->Col(); }
    CATCH
    return col;
}

} // extern "C"
