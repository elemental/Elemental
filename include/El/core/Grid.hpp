/*
   Copyright (c) 2009-2014, Jack Poulson
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
    int Row() const;           // MCRank()
    int Col() const;           // MRRank()
    int Rank() const;          // VCRank (VRRank) if COLUMN_MAJOR (ROW_MAJOR)
    int Height() const;        // MCSize()
    int Width() const;         // MRSize()
    int Size() const;          // VCSize() and VRSize()
    GridOrder Order() const;   // either COLUMN_MAJOR or ROW_MAJOR
    mpi::Comm ColComm() const; // MCComm()
    mpi::Comm RowComm() const; // MRComm()
    mpi::Comm Comm() const;    // VCComm (VRComm) if COLUMN_MAJOR (ROW_MAJOR)

    // Distribution-based interface
    int MCRank() const;
    int MRRank() const;
    int VCRank() const;
    int VRRank() const;
    int MCSize() const;
    int MRSize() const;
    int VCSize() const;
    int VRSize() const;
    mpi::Comm MCComm() const;
    mpi::Comm MRComm() const;
    mpi::Comm VCComm() const;
    mpi::Comm VRComm() const;
    mpi::Comm MDComm() const;
    mpi::Comm MDPerpComm() const;

    // Advanced routines
    explicit Grid
    ( mpi::Comm viewers, mpi::Group owners, int height, 
      GridOrder order=COLUMN_MAJOR );
    int GCD() const; // greatest common denominator of grid height and width
    int LCM() const; // lowest common multiple of grid height and width
    bool InGrid() const;
    bool HaveViewers() const;
    int OwningRank() const;
    int ViewingRank() const;
    int VCToViewingMap( int VCRank ) const;
    mpi::Group OwningGroup() const;
    mpi::Comm OwningComm() const;
    mpi::Comm ViewingComm() const;
    int DiagPath() const;
    int DiagPath( int vectorColRank ) const;
    int DiagPathRank() const;
    int DiagPathRank( int vectorColRank ) const;
    int FirstVCRank( int diagPath ) const;

    static int FindFactor( int p );

private:
    bool haveViewers_;
    int height_, size_, gcd_;
    GridOrder order_;
    std::vector<int> diagPathsAndRanks_;

    mpi::Comm viewingComm_; // all processes that create the grid
    mpi::Group viewingGroup_;
    std::vector<int> vectorColToViewingMap_;

    // Create a communicator for our owning team
    mpi::Comm owningComm_;
    mpi::Group owningGroup_;

    // These will only be valid if we are in the grid
    mpi::Comm cartComm_,  // the processes that are in the grid
              matrixColComm_, matrixRowComm_,
              matrixDiagComm_, matrixDiagPerpComm_,
              vectorColComm_, vectorRowComm_;

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

} // namespace El

#endif // ifndef EL_GRID_HPP
