/*
   Copyright (c) 2009-2011, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/
#include "elemental/environment.hpp"
#include "elemental/exports/c.h"

namespace {
// Since we can't return the class to C, instead return an integer ID and
// maintain a vector of all the created grids. This will need to be cleared
// by ElementalFinalize()
std::vector<elemental::Grid*> gridList;
}

extern "C" {

void ElementalClearGridList()
{
    for( int j=0; j< ::gridList.size(); ++j )
        delete ::gridList[j];
    ::gridList.clear();
}

int ElementalDefaultGrid( MPI_Comm comm )
{
    int gridHandle = ::gridList.size();
    ::gridList.push_back( new elemental::Grid( comm ) );
    return gridHandle;
}

int ElementalGrid( MPI_Comm comm, int r, int c )
{
    int gridHandle = ::gridList.size();
    ::gridList.push_back( new elemental::Grid( comm, r, c ) );
    return gridHandle;
}

// Don't yet export the viewingComm/owningGroup constructors.

int ElementalGridHeight( int gridHandle )
{
    return gridList[gridHandle]->Height();
}

int ElementalGridWidth( int gridHandle )
{
    return gridList[gridHandle]->Width();
}

int ElementalGridSize( int gridHandle )
{
    return gridList[gridHandle]->Size();
}

int ElementalInGrid( int gridHandle )
{
    return gridList[gridHandle]->InGrid();
}

int ElementalGridVCRank( int gridHandle )
{
    return gridList[gridHandle]->VCRank();
}

int ElementalGridVRRank( int gridHandle )
{
    return gridList[gridHandle]->VRRank();
}

int ElementalGridMCRank( int gridHandle )
{
    return gridList[gridHandle]->MCRank();
}

int ElementalGridMRRank( int gridHandle )
{
    return gridList[gridHandle]->MRRank();
}

} // extern "C"

