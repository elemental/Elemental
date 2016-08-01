/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>
#include <iomanip>

namespace {

// A (per-process) output file for logging
std::ofstream logFile;

}

namespace El {

void OpenLog( const char* filename )
{
    if( ::logFile.is_open() )
        CloseLog();
    ::logFile.open( filename );
}

std::ostream& LogOS()
{
    if( !::logFile.is_open() )
    {
        std::ostringstream fileOS;
        fileOS << "El-Proc" << std::setfill('0') << std::setw(3)
               << mpi::Rank() << ".log";
        ::logFile.open( fileOS.str().c_str() );
    }
    return ::logFile; 
}

void CloseLog() { ::logFile.close(); }

} // namespace El
