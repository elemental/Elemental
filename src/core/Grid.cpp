/*
   Copyright (c) 2009-2010, Jack Poulson
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
using namespace elemental;
using namespace std;

elemental::Grid::Grid
( MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("Grid::Grid");
#endif
    _inGrid = true; // this is true by assumption for this constructor

    // Extract our rank, the underlying group, and the number of processes
    MPI_Comm_dup( comm, &_viewingComm );
    MPI_Comm_group( _viewingComm, &_viewingGroup );
    MPI_Comm_rank( _viewingComm, &_viewingRank );
    MPI_Comm_size( _viewingComm, &_p );

    // All processes own the grid, so we have to trivially split _viewingGroup
    _owningGroup = _viewingGroup;
    _notOwningGroup = MPI_GROUP_EMPTY;
    _owningRank = _viewingRank;

    // Set up the map from the owning to viewing ranks.
    _owningToViewingMap.resize(_p);
    for( int i=0; i<_p; ++i )
        _owningToViewingMap[i] = i;

    // Factor p
    _r = static_cast<int>(sqrt(static_cast<double>(_p)));
    while( _p % _r != 0 )
        ++_r;
    _c = _p / _r;

    SetUpGrid();

#ifndef RELEASE
    PopCallStack();
#endif
}

elemental::Grid::Grid
( MPI_Comm comm, int r, int c )
{
#ifndef RELEASE
    PushCallStack("Grid::Grid");
#endif
    _inGrid = true; // this is true by assumption for this constructor

    // Extract our rank, the underlying group, and the number of processes
    MPI_Comm_dup( comm, &_viewingComm );
    MPI_Comm_group( _viewingComm, &_viewingGroup );
    MPI_Comm_rank( _viewingComm, &_viewingRank );
    MPI_Comm_size( _viewingComm, &_p );

    // All processes own the grid, so we have to trivially split _viewingGroup
    _owningGroup = _viewingGroup;
    _notOwningGroup = MPI_GROUP_EMPTY;
    _owningRank = _viewingRank;

    // Set up the map from the owning to viewing ranks.
    _owningToViewingMap.resize(_p);
    for( int i=0; i<_p; ++i )
        _owningToViewingMap[i] = i;

    _r = r;
    _c = c;
    if( _r < 0 || _c < 0 )
        throw logic_error( "Process grid dimensions must be non-negative." );

    SetUpGrid();

#ifndef RELEASE
    PopCallStack();
#endif
}

elemental::Grid::Grid
( MPI_Comm viewingComm, MPI_Group owningGroup )
{
#ifndef RELEASE
    PushCallStack("Grid::Grid");
#endif

    // Extract our rank and the underlying group from the viewing comm
    MPI_Comm_dup( viewingComm, &_viewingComm );
    MPI_Comm_group( _viewingComm, &_viewingGroup );
    MPI_Comm_rank( _viewingComm, &_viewingRank );

    // Extract our rank and the number of processes from the owning group
    _owningGroup = owningGroup;
    MPI_Group_size( _owningGroup, &_p );
    MPI_Group_rank( _owningGroup, &_owningRank );
    _inGrid = ( _owningRank != MPI_UNDEFINED );

    // Create the complement of the owning group
    MPI_Group_difference( _viewingGroup, _owningGroup, &_notOwningGroup );

    // Set up the map from the owning to viewing ranks
    std::vector<int> ranks(_p);
    for( int i=0; i<_p; ++i )
        ranks[i] = i;
    _owningToViewingMap.resize(_p);
    MPI_Group_translate_ranks
    ( _owningGroup, _p, &ranks[0], _viewingGroup, &_owningToViewingMap[0] );

    // Factor p
    _r = static_cast<int>(sqrt(static_cast<double>(_p)));
    while( _p % _r != 0 )
        ++_r;
    _c = _p / _r;

    SetUpGrid();

#ifndef RELEASE
    PopCallStack();
#endif
}

// Currently forces a columnMajor absolute rank on the grid
elemental::Grid::Grid
( MPI_Comm viewingComm, MPI_Group owningGroup, int r, int c )
{
#ifndef RELEASE
    PushCallStack("Grid::Grid");
#endif

    // Extract our rank and the underlying group from the viewing comm
    MPI_Comm_dup( viewingComm, &_viewingComm );
    MPI_Comm_group( _viewingComm, &_viewingGroup );
    MPI_Comm_rank( _viewingComm, &_viewingRank );

    // Extract our rank and the number of processes from the owning group
    _owningGroup = owningGroup;
    MPI_Group_size( _owningGroup, &_p );
    MPI_Group_rank( _owningGroup, &_owningRank );
    _inGrid = ( _owningRank != MPI_UNDEFINED );

    // Create the complement of the owning group
    MPI_Group_difference( _viewingGroup, _owningGroup, &_notOwningGroup );

    // Set up the map from the owning to viewing ranks
    std::vector<int> ranks(_p);
    for( int i=0; i<_p; ++i )
        ranks[i] = i;
    _owningToViewingMap.resize(_p);
    MPI_Group_translate_ranks
    ( _owningGroup, _p, &ranks[0], _viewingGroup, &_owningToViewingMap[0] );

    _r = r;
    _c = c;

    if( _r < 0 || _c < 0 )
        throw logic_error( "Process grid dimensions must be non-negative." );
    
    SetUpGrid();

#ifndef RELEASE
    PopCallStack();
#endif
}

void
elemental::Grid::SetUpGrid()
{
#ifndef RELEASE
    PushCallStack("Grid::SetUpGrid");
#endif
    if( _p != _r*_c )
    {
        ostringstream msg;
        msg << "Number of processes must match grid size:\n"
            << "  p=" << _p << ", (r,c)=(" << _r << "," << _c << ")";
        throw logic_error( msg.str() );
    }

    _gcd = utilities::GCD( _r, _c );
    int lcm = _p / _gcd;

#ifndef RELEASE
    if( _owningRank == 0 )
    {
        cout << "Building process grid with:\n"
             << "  p=" << _p << ", (r,c)=(" << _r << "," << _c << ")\n"
             << "  gcd=" << _gcd << endl;
    }
#endif

    // Split the viewing comm into the owning and not owning subsets
    if( _inGrid )
        MPI_Comm_create( _viewingComm, _owningGroup, &_owningComm );    
    else
        MPI_Comm_create( _viewingComm, _notOwningGroup, &_notOwningComm );

    if( _inGrid )
    {
        // Create a cartesian communicator
        int dimensions[2] = { _c, _r };
        int periods[2] = { true, true };
        int reorder = false;
        MPI_Cart_create
        ( _owningComm, 2, dimensions, periods, reorder, &_cartComm );

        // Set up the MatrixCol and MatrixRow communicators
        int remainingDimensions[2];
        remainingDimensions[0] = false;
        remainingDimensions[1] = true;
        MPI_Cart_sub( _cartComm, remainingDimensions, &_matrixColComm );
        remainingDimensions[0] = true;
        remainingDimensions[1] = false;
        MPI_Cart_sub( _cartComm, remainingDimensions, &_matrixRowComm );
        MPI_Comm_rank( _matrixColComm, &_matrixColRank );
        MPI_Comm_rank( _matrixRowComm, &_matrixRowRank );

        // Set up the VectorCol and VectorRow communicators
        _vectorColRank = _matrixColRank + _r*_matrixRowRank;
        _vectorRowRank = _matrixRowRank + _c*_matrixColRank;
        MPI_Comm_split( _cartComm, 0, _vectorColRank, &_vectorColComm );
        MPI_Comm_split( _cartComm, 0, _vectorRowRank, &_vectorRowComm );

        // Compute which diagonal 'path' we're in, and what our rank is, then
        // perform AllGather world to store everyone's info
        _diagPathsAndRanks.resize(2*_p);
        vector<int> myDiagPathAndRank(2);
        myDiagPathAndRank[0] = (_matrixRowRank+_r-_matrixColRank) % _gcd;
        int diagPathRank = 0;
        int row = 0;
        int col = myDiagPathAndRank[0];
        for( int j=0; j<lcm; ++j )
        {
            if( row == _matrixColRank && col == _matrixRowRank )
            {
                myDiagPathAndRank[1] = diagPathRank;
                break;
            }
            else
            {
                row = (row + 1) % _r;
                col = (col + 1) % _c;
                ++diagPathRank;
            }
        }
        wrappers::mpi::AllGather
        ( &myDiagPathAndRank[0], 2, &_diagPathsAndRanks[0], 2, _vectorColComm );

#ifndef RELEASE
        MPI_Errhandler_set( _matrixColComm, MPI_ERRORS_RETURN );
        MPI_Errhandler_set( _matrixRowComm, MPI_ERRORS_RETURN );
        MPI_Errhandler_set( _vectorColComm, MPI_ERRORS_RETURN );
        MPI_Errhandler_set( _vectorRowComm, MPI_ERRORS_RETURN );
#endif
    }
    else
    {
        // diag paths and ranks are implicitly set to undefined
        _matrixColRank = MPI_UNDEFINED;
        _matrixRowRank = MPI_UNDEFINED;
        _vectorColRank = MPI_UNDEFINED;
        _vectorRowRank = MPI_UNDEFINED;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

elemental::Grid::~Grid()
{
    int finalized;
    MPI_Finalized( &finalized );
    if( !finalized )
    {
        if( _inGrid )
        {
            MPI_Comm_free( &_matrixColComm );
            MPI_Comm_free( &_matrixRowComm );
            MPI_Comm_free( &_vectorColComm );
            MPI_Comm_free( &_vectorRowComm );
            MPI_Comm_free( &_cartComm );
        }

        if( _inGrid )
            MPI_Comm_free( &_owningComm );
        else
            MPI_Comm_free( &_notOwningComm );

        if( _notOwningGroup != MPI_GROUP_EMPTY )
            MPI_Group_free( &_notOwningGroup );

        MPI_Comm_free( &_viewingComm );
    }
}


