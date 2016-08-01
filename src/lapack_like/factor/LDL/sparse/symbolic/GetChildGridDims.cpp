/*
   Copyright 2009-2011, Jack Poulson.
   All rights reserved.

   Copyright 2011-2012, Jack Poulson, Lexing Ying, and 
   The University of Texas at Austin.
   All rights reserved.

   Copyright 2013, Jack Poulson, Lexing Ying, and Stanford University.
   All rights reserved.

   Copyright 2013-2014, Jack Poulson and The Georgia Institute of Technology.
   All rights reserved.

   Copyright 2014-2015, Jack Poulson and Stanford University.
   All rights reserved.
   
   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {
namespace ldl {

void GetChildGridDims
( const DistNodeInfo& info, vector<int>& gridHeights, vector<int>& gridWidths )
{
    gridHeights.resize( 2 );
    gridWidths.resize( 2 );
    gridHeights[0] = gridHeights[1] = 0;
    gridWidths[0] = gridWidths[1] = 0;

    const bool onLeft = info.child->onLeft;
    const int childTeamRank = mpi::Rank( info.child->comm );
    if( onLeft && childTeamRank == 0 )
    {
        gridHeights[0] = info.child->grid->Height();
        gridWidths[0] = info.child->grid->Width();
    }
    else if( !onLeft && childTeamRank == 0 )
    {
        gridHeights[1] = info.child->grid->Height();
        gridWidths[1] = info.child->grid->Width();
    }
    mpi::AllReduce( gridHeights.data(), 2, info.comm );
    mpi::AllReduce( gridWidths.data(), 2, info.comm );
}

} // namespace ldl
} // namespace El
