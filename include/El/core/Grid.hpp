/*
   Copyright (c) 2009-2015, Jack Poulson
                      2013, Jed Brown 
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_GRID_HPP
#define EL_GRID_HPP

namespace El {

class Grid
{
public:
    explicit Grid
    ( mpi::Comm comm=mpi::COMM_WORLD, GridOrder order=COLUMN_MAJOR );
    explicit Grid( mpi::Comm comm, int height, GridOrder order=COLUMN_MAJOR );
    ~Grid();

    // Simple interface (simpler version of distributed-based interface)
    int Row() const; // MCRank()
    int Col() const; // MRRank()
    int Height() const EL_NOEXCEPT; // MCSize()
    int Width() const EL_NOEXCEPT;  // MRSize()
    int Size() const EL_NOEXCEPT;   // VCSize() and VRSize()
    int Rank() const;          // same as OwningRank()
    GridOrder Order() const;   // either COLUMN_MAJOR or ROW_MAJOR
    mpi::Comm ColComm() const EL_NOEXCEPT; // MCComm()
    mpi::Comm RowComm() const EL_NOEXCEPT; // MRComm()
    // VCComm (VRComm) if COLUMN_MAJOR (ROW_MAJOR)
    mpi::Comm Comm() const EL_NOEXCEPT;

    // Distribution-based interface
    int MCRank() const;
    int MRRank() const;
    int VCRank() const;
    int VRRank() const;
    int MCSize() const EL_NOEXCEPT;
    int MRSize() const EL_NOEXCEPT;
    int VCSize() const EL_NOEXCEPT;
    int VRSize() const EL_NOEXCEPT;
    mpi::Comm MCComm() const EL_NOEXCEPT;
    mpi::Comm MRComm() const EL_NOEXCEPT;
    mpi::Comm VCComm() const EL_NOEXCEPT;
    mpi::Comm VRComm() const EL_NOEXCEPT;
    mpi::Comm MDComm() const EL_NOEXCEPT;
    mpi::Comm MDPerpComm() const EL_NOEXCEPT;

    // Advanced routines
    explicit Grid
    ( mpi::Comm viewers, mpi::Group owners, int height, 
      GridOrder order=COLUMN_MAJOR );
    // greatest common denominator of grid height and width
    int GCD() const EL_NOEXCEPT;
    // lowest common multiple of grid height and width
    int LCM() const EL_NOEXCEPT;
    bool InGrid() const;
    bool HaveViewers() const EL_NOEXCEPT;
    int OwningRank() const;
    int ViewingRank() const;

    mpi::Group OwningGroup() const EL_NOEXCEPT;
    mpi::Comm OwningComm() const EL_NOEXCEPT;
    mpi::Comm ViewingComm() const EL_NOEXCEPT;
    int Diag() const;
    int Diag( int vcRank ) const EL_NOEXCEPT;
    int DiagRank() const;
    int DiagRank( int vcRank ) const EL_NOEXCEPT;

    int VCToVR( int vcRank ) const EL_NOEXCEPT;
    int VRToVC( int vrRank ) const EL_NOEXCEPT;
    int CoordsToVC
    ( Dist colDist, Dist rowDist, 
      int distRank, int crossRank=0, int redundant=0 ) const;
    int VCToViewing( int VCRank ) const EL_NOEXCEPT;

    static int FindFactor( int p ) EL_NOEXCEPT;

private:
    bool haveViewers_;
    int height_, size_, gcd_;
    GridOrder order_;
    vector<int> diagsAndRanks_;

    mpi::Comm viewingComm_; // all processes that create the grid
    mpi::Group viewingGroup_;
    vector<int> vcToViewing_;

    // Create a communicator for our owning team
    mpi::Comm owningComm_;
    mpi::Group owningGroup_;

    // These will only be valid if we are in the grid
    mpi::Comm cartComm_,  // the processes that are in the grid
              mcComm_, mrComm_,
              mdComm_, mdPerpComm_,
              vcComm_, vrComm_;

    void SetUpGrid();

    // Disable copying this class due to MPI_Comm/MPI_Group ownership issues
    // and potential performance loss from duplicating MPI communicators, e.g.,
    // on Blue Gene/P there is supposedly a performance loss
    const Grid& operator=( Grid& );
    Grid( const Grid& );
};

bool operator== ( const Grid& A, const Grid& B );
bool operator!= ( const Grid& A, const Grid& B );

// Return a grid constructed using mpi::COMM_WORLD.
const Grid& DefaultGrid();

inline void AssertSameGrids( const Grid& g1 ) { }

inline void AssertSameGrids( const Grid& g1, const Grid& g2 )
{
    if( g1 != g2 )
        LogicError("Grids did not match");
}

template<typename... Args>
inline void AssertSameGrids( const Grid& g1, const Grid& g2, Args&... args )
{
    if( g1 != g2 )
        LogicError("Grids did not match");
    AssertSameGrids( g2, args... );
}

} // namespace El

#endif // ifndef EL_GRID_HPP
