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
#ifndef ELEMENTAL_PROCESSGRID_HPP
#define ELEMENTAL_PROCESSGRID_HPP 1

namespace elemental {

class Grid
{
    int _p;
    int _c;
    int _r;
    int _gcd;
    int _col;
    int _row;
    int _matrixColRank;
    int _matrixRowRank;
    int _vectorColRank;
    int _vectorRowRank;
    std::vector<int> _diagPathsAndRanks;

    mpi::Comm _viewingComm; // all processes that create the grid
    mpi::Group _viewingGroup;
    int _viewingRank; // our rank in the viewing communicator

    mpi::Group _owningGroup; // the processes that can own data
    mpi::Group _notOwningGroup; // contains the remaining processes

    std::vector<int> _vectorColToViewingMap;

    // Keep track of whether or not our process is in the grid. This is 
    // necessary to avoid calls like MPI_Comm_size when we're not in the
    // communicator's group. Note that we _can_ call MPI_Group_rank when not 
    // in the group and that the result is MPI_UNDEFINED.
    bool _inGrid;

    // Create a communicator for the processes that are in the process grid
    mpi::Comm _owningComm;
    mpi::Comm _notOwningComm; // necessary complimentary communicator
    int _owningRank;

    // These will only be valid if we are in the grid
    mpi::Comm _cartComm;  // the processes that are in the grid
    mpi::Comm _matrixColComm;
    mpi::Comm _matrixRowComm;
    mpi::Comm _vectorColComm;
    mpi::Comm _vectorRowComm;

    void SetUpGrid();

    // Disable copying this class due to MPI_Comm/MPI_Group ownership issues
    // and potential performance loss from duplicating MPI communicators, e.g.,
    // on Blue Gene/P there is supposedly a performance loss
    void operator=( Grid& );
    Grid( const Grid& );

    public:

    // For constructing grids where every process is a member
    Grid( mpi::Comm comm );
    Grid( mpi::Comm comm, int r, int c );
    
    // For constructing grids where only the 'owningGroup' processes are in the
    // grid. viewingComm must be valid for all processes creating the 
    // grid, not just those in the owning group.
    Grid( mpi::Comm viewingComm, mpi::Group owningGroup );
    Grid( mpi::Comm viewingComm, mpi::Group owningGroup, 
          int r, int c );

    ~Grid();

    bool InGrid() const;
        
    int Size() const;
    int Height() const;
    int Width() const;
    int GCD() const;
    int LCM() const;
    int DiagPath() const;
    int DiagPath( int vectorColRank ) const;
    int DiagPathRank() const;
    int DiagPathRank( int vectorColRank ) const;
    int MCRank() const;
    int MRRank() const;
    int VCRank() const;
    int VRRank() const;
    int OwningRank() const;
    int ViewingRank() const;
    int VCToViewingMap( int VCRank ) const;
    mpi::Group OwningGroup() const;
    mpi::Comm OwningComm() const;
    mpi::Comm ViewingComm() const;
    mpi::Comm MCComm() const;
    mpi::Comm MRComm() const;
    mpi::Comm VCComm() const;
    mpi::Comm VRComm() const;
};

bool operator== ( const Grid& A, const Grid& B );
bool operator!= ( const Grid& A, const Grid& B );

} // elemental

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

inline bool
elemental::Grid::InGrid() const
{ return _inGrid; }

inline int
elemental::Grid::Size() const
{ return _p; }

inline int
elemental::Grid::Height() const
{ return _r; }

inline int
elemental::Grid::Width() const
{ return _c; }

inline int
elemental::Grid::GCD() const
{ return _gcd; }

inline int
elemental::Grid::LCM() const
{ return _p/_gcd; }

inline int
elemental::Grid::DiagPath() const
{ 
    if( _inGrid )
        return _diagPathsAndRanks[2*_vectorColRank]; 
    else
        return mpi::UNDEFINED;
}

inline int
elemental::Grid::DiagPath( int vectorColRank ) const
{ 
    if( vectorColRank != mpi::UNDEFINED )
        return _diagPathsAndRanks[2*vectorColRank]; 
    else
        return mpi::UNDEFINED;
}

inline int
elemental::Grid::DiagPathRank() const
{ 
    if( _inGrid )
        return _diagPathsAndRanks[2*_vectorColRank+1];
    else
        return mpi::UNDEFINED;
}

inline int
elemental::Grid::DiagPathRank( int vectorColRank ) const
{ 
    if( vectorColRank != mpi::UNDEFINED )
        return _diagPathsAndRanks[2*vectorColRank+1]; 
    else
        return mpi::UNDEFINED;
}

inline int
elemental::Grid::MCRank() const
{ return _matrixColRank; }

inline int
elemental::Grid::MRRank() const
{ return _matrixRowRank; }

inline int
elemental::Grid::VCRank() const
{ return _vectorColRank; }

inline int
elemental::Grid::VRRank() const
{ return _vectorRowRank; }

inline int
elemental::Grid::OwningRank() const
{ return _owningRank; }

inline int
elemental::Grid::ViewingRank() const
{ return _viewingRank; }

inline int
elemental::Grid::VCToViewingMap( int VCRank ) const
{ return _vectorColToViewingMap[VCRank]; }

inline elemental::mpi::Group
elemental::Grid::OwningGroup() const
{ return _owningGroup; }

inline elemental::mpi::Comm
elemental::Grid::OwningComm() const
{ return _owningComm; }

inline elemental::mpi::Comm
elemental::Grid::ViewingComm() const
{ return _viewingComm; }

inline elemental::mpi::Comm
elemental::Grid::MCComm() const
{ return _matrixColComm; }

inline elemental::mpi::Comm
elemental::Grid::MRComm() const
{ return _matrixRowComm; }

inline elemental::mpi::Comm
elemental::Grid::VCComm() const
{ return _vectorColComm; }

inline elemental::mpi::Comm
elemental::Grid::VRComm() const
{ return _vectorRowComm; }

inline bool 
elemental::operator== ( const Grid& A, const Grid& B )
{ return &A == &B; }

inline bool 
elemental::operator!= ( const Grid& A, const Grid& B )
{ return &A != &B; }

inline 
elemental::Grid::Grid
( mpi::Comm comm )
{
#ifndef RELEASE
    PushCallStack("Grid::Grid");
#endif
    _inGrid = true; // this is true by assumption for this constructor

    // Extract our rank, the underlying group, and the number of processes
    mpi::CommDup( comm, _viewingComm );
    mpi::CommGroup( _viewingComm, _viewingGroup );
    _viewingRank = mpi::CommRank( _viewingComm );
    _p = mpi::CommSize( _viewingComm );

    // All processes own the grid, so we have to trivially split _viewingGroup
    _owningGroup = _viewingGroup;
    _notOwningGroup = mpi::GROUP_EMPTY;
    _owningRank = _viewingRank;

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

inline
elemental::Grid::Grid
( mpi::Comm comm, int r, int c )
{
#ifndef RELEASE
    PushCallStack("Grid::Grid");
#endif
    _inGrid = true; // this is true by assumption for this constructor

    // Extract our rank, the underlying group, and the number of processes
    mpi::CommDup( comm, _viewingComm );
    mpi::CommGroup( _viewingComm, _viewingGroup );
    _viewingRank = mpi::CommRank( _viewingComm );
    _p = mpi::CommSize( _viewingComm );

    // All processes own the grid, so we have to trivially split _viewingGroup
    _owningGroup = _viewingGroup;
    _notOwningGroup = mpi::GROUP_EMPTY;
    _owningRank = _viewingRank;

    _r = r;
    _c = c;
    if( _r < 0 || _c < 0 )
        throw std::logic_error
        ("Process grid dimensions must be non-negative");

    SetUpGrid();

#ifndef RELEASE
    PopCallStack();
#endif
}

inline
elemental::Grid::Grid
( mpi::Comm viewingComm, mpi::Group owningGroup )
{
#ifndef RELEASE
    PushCallStack("Grid::Grid");
#endif

    // Extract our rank and the underlying group from the viewing comm
    mpi::CommDup( viewingComm, _viewingComm );
    mpi::CommGroup( _viewingComm, _viewingGroup );
    _viewingRank = mpi::CommRank( _viewingComm );

    // Extract our rank and the number of processes from the owning group
    _owningGroup = owningGroup;
    _p = mpi::GroupSize( _owningGroup );
    _owningRank = mpi::GroupRank( _owningGroup );
    _inGrid = ( _owningRank != mpi::UNDEFINED );

    // Create the complement of the owning group
    mpi::GroupDifference
    ( _viewingGroup, _owningGroup, _notOwningGroup );

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
inline
elemental::Grid::Grid
( mpi::Comm viewingComm, mpi::Group owningGroup, 
  int r, int c )
{
#ifndef RELEASE
    PushCallStack("Grid::Grid");
#endif

    // Extract our rank and the underlying group from the viewing comm
    mpi::CommDup( viewingComm, _viewingComm );
    mpi::CommGroup( _viewingComm, _viewingGroup );
    _viewingRank = mpi::CommRank( _viewingComm );

    // Extract our rank and the number of processes from the owning group
    _owningGroup = owningGroup;
    _p = mpi::GroupSize( _owningGroup );
    _owningRank = mpi::GroupRank( _owningGroup );
    _inGrid = ( _owningRank != mpi::UNDEFINED );

    // Create the complement of the owning group
    mpi::GroupDifference
    ( _viewingGroup, _owningGroup, _notOwningGroup );

    _r = r;
    _c = c;

    if( _r < 0 || _c < 0 )
        throw std::logic_error
        ("Process grid dimensions must be non-negative");

    SetUpGrid();

#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
elemental::Grid::SetUpGrid()
{
#ifndef RELEASE
    PushCallStack("Grid::SetUpGrid");
#endif
    if( _p != _r*_c )
    {
        std::ostringstream msg;
        msg << "Number of processes must match grid size:\n"
            << "  p=" << _p << ", (r,c)=(" << _r << "," << _c << ")";
        throw std::logic_error( msg.str().c_str() );
    }

    _gcd = elemental::GCD( _r, _c );
    int lcm = _p / _gcd;

#ifndef RELEASE
    if( _owningRank == 0 )
    {
        std::cout << "Building process grid with:\n"
                  << "  p=" << _p << ", (r,c)=(" << _r << "," << _c << ")\n"
                  << "  gcd=" << _gcd << std::endl;
    }
#endif

    // Split the viewing comm into the owning and not owning subsets
    if( _inGrid )
        mpi::CommSplit( _viewingComm, true, _owningRank, _owningComm );
    else
        mpi::CommSplit( _viewingComm, false, 0, _notOwningComm );

    if( _inGrid )
    {
        // Create a cartesian communicator
        int dimensions[2] = { _c, _r };
        int periods[2] = { true, true };
        int reorder = false;
        mpi::CartCreate
        ( _owningComm, 2, dimensions, periods, reorder, _cartComm );

        // Set up the MatrixCol and MatrixRow communicators
        int remainingDimensions[2];
        remainingDimensions[0] = false;
        remainingDimensions[1] = true;
        mpi::CartSub( _cartComm, remainingDimensions, _matrixColComm );
        remainingDimensions[0] = true;
        remainingDimensions[1] = false;
        mpi::CartSub( _cartComm, remainingDimensions, _matrixRowComm );
        _matrixColRank = mpi::CommRank( _matrixColComm );
        _matrixRowRank = mpi::CommRank( _matrixRowComm );

        // Set up the VectorCol and VectorRow communicators
        _vectorColRank = _matrixColRank + _r*_matrixRowRank;
        _vectorRowRank = _matrixRowRank + _c*_matrixColRank;
        mpi::CommSplit( _cartComm, 0, _vectorColRank, _vectorColComm );
        mpi::CommSplit( _cartComm, 0, _vectorRowRank, _vectorRowComm );

        // Compute which diagonal 'path' we're in, and what our rank is, then
        // perform AllGather world to store everyone's info
        _diagPathsAndRanks.resize(2*_p);
        std::vector<int> myDiagPathAndRank(2);
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
        mpi::AllGather
        ( &myDiagPathAndRank[0], 2, &_diagPathsAndRanks[0], 2, _vectorColComm );

#ifndef RELEASE
        mpi::ErrorHandlerSet
        ( _matrixColComm, mpi::ERRORS_RETURN );
        mpi::ErrorHandlerSet
        ( _matrixRowComm, mpi::ERRORS_RETURN );
        mpi::ErrorHandlerSet
        ( _vectorColComm, mpi::ERRORS_RETURN );
        mpi::ErrorHandlerSet
        ( _vectorRowComm, mpi::ERRORS_RETURN );
#endif
    }
    else
    {
        // diag paths and ranks are implicitly set to undefined
        _matrixColRank = mpi::UNDEFINED;
        _matrixRowRank = mpi::UNDEFINED;
        _vectorColRank = mpi::UNDEFINED;
        _vectorRowRank = mpi::UNDEFINED;
    }

    // Set up the map from the VC group to the _viewingGroup ranks.
    // Since the VC communicator preserves the ordering of the _owningGroup
    // ranks, we can simply translate from _owningGroup.
    std::vector<int> ranks(_p);
    for( int i=0; i<_p; ++i )
        ranks[i] = i;
    _vectorColToViewingMap.resize(_p);
    mpi::GroupTranslateRanks
    ( _owningGroup, _p, &ranks[0], _viewingGroup, &_vectorColToViewingMap[0] );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline
elemental::Grid::~Grid()
{
    if( !mpi::Finalized() )
    {
        if( _inGrid )
        {
            mpi::CommFree( _matrixColComm );
            mpi::CommFree( _matrixRowComm );
            mpi::CommFree( _vectorColComm );
            mpi::CommFree( _vectorRowComm );
            mpi::CommFree( _cartComm );
        }

        if( _inGrid )
            mpi::CommFree( _owningComm );
        else
            mpi::CommFree( _notOwningComm );

        if( _notOwningGroup != mpi::GROUP_EMPTY )
            mpi::GroupFree( _notOwningGroup );

        mpi::CommFree( _viewingComm );
    }
}

#endif /* ELEMENTAL_PROCESSGRID_HPP */

