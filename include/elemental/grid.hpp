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

    imports::mpi::Comm _viewingComm; // all processes that create the grid
    imports::mpi::Group _viewingGroup;
    int _viewingRank; // our rank in the viewing communicator

    imports::mpi::Group _owningGroup; // the processes that can own data
    imports::mpi::Group _notOwningGroup; // contains the remaining processes

    std::vector<int> _vectorColToViewingMap;

    // Keep track of whether or not our process is in the grid. This is 
    // necessary to avoid calls like MPI_Comm_size when we're not in the
    // communicator's group. Note that we _can_ call MPI_Group_rank when not 
    // in the group and that the result is MPI_UNDEFINED.
    bool _inGrid;

    // Create a communicator for the processes that are in the process grid
    imports::mpi::Comm _owningComm;
    imports::mpi::Comm _notOwningComm; // necessary complimentary communicator
    int _owningRank;

    // These will only be valid if we are in the grid
    imports::mpi::Comm _cartComm;  // the processes that are in the grid
    imports::mpi::Comm _matrixColComm;
    imports::mpi::Comm _matrixRowComm;
    imports::mpi::Comm _vectorColComm;
    imports::mpi::Comm _vectorRowComm;

    void SetUpGrid();

    // Disable copying this class due to MPI_Comm/MPI_Group ownership issues
    // and potential performance loss from duplicating MPI communicators, e.g.,
    // on Blue Gene/P there is supposedly a performance loss
    void operator=( Grid& );
    Grid( const Grid& );

    public:

    // For constructing grids where every process is a member
    Grid( imports::mpi::Comm comm );
    Grid( imports::mpi::Comm comm, int r, int c );
    
    // For constructing grids where only the 'owningGroup' processes are in the
    // grid. viewingComm must be valid for all processes creating the 
    // grid, not just those in the owning group.
    Grid( imports::mpi::Comm viewingComm, imports::mpi::Group owningGroup );
    Grid( imports::mpi::Comm viewingComm, imports::mpi::Group owningGroup, 
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
    imports::mpi::Group OwningGroup() const;
    imports::mpi::Comm OwningComm() const;
    imports::mpi::Comm ViewingComm() const;
    imports::mpi::Comm MCComm() const;
    imports::mpi::Comm MRComm() const;
    imports::mpi::Comm VCComm() const;
    imports::mpi::Comm VRComm() const;
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
        return imports::mpi::UNDEFINED;
}

inline int
elemental::Grid::DiagPath( int vectorColRank ) const
{ 
    if( vectorColRank != imports::mpi::UNDEFINED )
        return _diagPathsAndRanks[2*vectorColRank]; 
    else
        return imports::mpi::UNDEFINED;
}

inline int
elemental::Grid::DiagPathRank() const
{ 
    if( _inGrid )
        return _diagPathsAndRanks[2*_vectorColRank+1];
    else
        return imports::mpi::UNDEFINED;
}

inline int
elemental::Grid::DiagPathRank( int vectorColRank ) const
{ 
    if( vectorColRank != imports::mpi::UNDEFINED )
        return _diagPathsAndRanks[2*vectorColRank+1]; 
    else
        return imports::mpi::UNDEFINED;
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

inline elemental::imports::mpi::Group
elemental::Grid::OwningGroup() const
{ return _owningGroup; }

inline elemental::imports::mpi::Comm
elemental::Grid::OwningComm() const
{ return _owningComm; }

inline elemental::imports::mpi::Comm
elemental::Grid::ViewingComm() const
{ return _viewingComm; }

inline elemental::imports::mpi::Comm
elemental::Grid::MCComm() const
{ return _matrixColComm; }

inline elemental::imports::mpi::Comm
elemental::Grid::MRComm() const
{ return _matrixRowComm; }

inline elemental::imports::mpi::Comm
elemental::Grid::VCComm() const
{ return _vectorColComm; }

inline elemental::imports::mpi::Comm
elemental::Grid::VRComm() const
{ return _vectorRowComm; }

inline bool 
elemental::operator== ( const Grid& A, const Grid& B )
{ return &A == &B; }

inline bool 
elemental::operator!= ( const Grid& A, const Grid& B )
{ return &A != &B; }

#endif /* ELEMENTAL_PROCESSGRID_HPP */

