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

   Copyright 2016, Jack Poulson.
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {
namespace ldl {

NodeInfo::NodeInfo( NodeInfo* parentNode )
: parent(parentNode)
{ }

NodeInfo::NodeInfo( DistNodeInfo* duplicateNode )
: duplicate(duplicateNode)
{
    EL_DEBUG_CSE

    size = duplicate->size;
    off = duplicate->off;
    origLowerStruct = duplicate->origLowerStruct;

    myOff = duplicate->myOff;
    lowerStruct = duplicate->lowerStruct;
    origLowerRelInds = duplicate->origLowerRelInds;
}

NodeInfo::~NodeInfo() { }

DistNodeInfo::DistNodeInfo( const El::Grid& grid )
: rootGrid_(&grid)
{ }

DistNodeInfo::DistNodeInfo( DistNodeInfo* parentNode )
: parent(parentNode)
{ }

DistNodeInfo::~DistNodeInfo() { }

void DistNodeInfo::SetRootGrid( const El::Grid& grid )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( parent != nullptr )
          LogicError("This is not a root node");
    )
    rootGrid_ = &grid;
}

void DistNodeInfo::AssignGrid( unique_ptr<El::Grid>& grid )
{
    EL_DEBUG_CSE
    grid_ = std::move(grid);
}

const El::Grid& DistNodeInfo::Grid() const
{
    EL_DEBUG_CSE
    if( parent == nullptr )
    {
        return *rootGrid_;
    }
    else
    {
        const El::Grid* gridPtr = grid_.get();
        if( gridPtr == nullptr )
            LogicError("Tried to return pointer to non-existent Grid");
        return *gridPtr;
    }
}

void DistNodeInfo::GetChildGridDims
( vector<int>& gridHeights, vector<int>& gridWidths ) const
{
    EL_DEBUG_CSE

    // We will pack this with (gridHeight0,gridWidth0,gridHeight1,gridWidth1).
    vector<int> gridDims(4,0);
    if( child->onLeft )
    {
        if( child->Grid().Rank() == 0 )
        {
            gridDims[0] = child->Grid().Height();
            gridDims[1] = child->Grid().Width();
        }
    }
    else
    {
        if( child->Grid().Rank() == 0 )
        {
            gridDims[2] = child->Grid().Height();
            gridDims[3] = child->Grid().Width();
        }
    }
    mpi::AllReduce( gridDims.data(), 4, Grid().Comm() );

    gridHeights.resize( 2 );
    gridWidths.resize( 2 );
    gridHeights[0] = gridDims[0];
    gridWidths[0] = gridDims[1];
    gridHeights[1] = gridDims[2];
    gridWidths[1] = gridDims[3];
}

} // namespace ldl
} // namespace El
