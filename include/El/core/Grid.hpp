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
    int Row() const EL_NO_RELEASE_EXCEPT; // MCRank()
    int Col() const EL_NO_RELEASE_EXCEPT; // MRRank()
    int Height() const EL_NO_EXCEPT;       // MCSize()
    int Width() const EL_NO_EXCEPT;        // MRSize()
    int Size() const EL_NO_EXCEPT;         // VCSize() and VRSize()
    int Rank() const EL_NO_RELEASE_EXCEPT; // same as OwningRank()
    GridOrder Order() const EL_NO_EXCEPT;  // either COLUMN_MAJOR or ROW_MAJOR
    mpi::Comm ColComm() const EL_NO_EXCEPT; // MCComm()
    mpi::Comm RowComm() const EL_NO_EXCEPT; // MRComm()
    // VCComm (VRComm) if COLUMN_MAJOR (ROW_MAJOR)
    mpi::Comm Comm() const EL_NO_EXCEPT;

    // Distribution-based interface
    int MCRank() const EL_NO_RELEASE_EXCEPT;
    int MRRank() const EL_NO_RELEASE_EXCEPT;
    int VCRank() const EL_NO_RELEASE_EXCEPT;
    int VRRank() const EL_NO_RELEASE_EXCEPT;
    int MCSize() const EL_NO_EXCEPT;
    int MRSize() const EL_NO_EXCEPT;
    int VCSize() const EL_NO_EXCEPT;
    int VRSize() const EL_NO_EXCEPT;
    mpi::Comm MCComm() const EL_NO_EXCEPT;
    mpi::Comm MRComm() const EL_NO_EXCEPT;
    mpi::Comm VCComm() const EL_NO_EXCEPT;
    mpi::Comm VRComm() const EL_NO_EXCEPT;
    mpi::Comm MDComm() const EL_NO_EXCEPT;
    mpi::Comm MDPerpComm() const EL_NO_EXCEPT;

    // Advanced routines
    explicit Grid
    ( mpi::Comm viewers, mpi::Group owners, int height, 
      GridOrder order=COLUMN_MAJOR );
    // greatest common denominator of grid height and width
    int GCD() const EL_NO_EXCEPT;
    // lowest common multiple of grid height and width
    int LCM() const EL_NO_EXCEPT;
    bool InGrid() const EL_NO_RELEASE_EXCEPT;
    bool HaveViewers() const EL_NO_EXCEPT;
    int OwningRank() const EL_NO_RELEASE_EXCEPT;
    int ViewingRank() const EL_NO_RELEASE_EXCEPT;

    mpi::Group OwningGroup() const EL_NO_EXCEPT;
    mpi::Comm OwningComm() const EL_NO_EXCEPT;
    mpi::Comm ViewingComm() const EL_NO_EXCEPT;
    int Diag() const EL_NO_RELEASE_EXCEPT;
    int Diag( int vcRank ) const EL_NO_EXCEPT;
    int DiagRank() const EL_NO_RELEASE_EXCEPT;
    int DiagRank( int vcRank ) const EL_NO_EXCEPT;

    int VCToVR( int vcRank ) const EL_NO_EXCEPT;
    int VRToVC( int vrRank ) const EL_NO_EXCEPT;
    int CoordsToVC
    ( Dist colDist, Dist rowDist, 
      int distRank, int crossRank=0, int redundant=0 ) const
    EL_NO_RELEASE_EXCEPT;
    int VCToViewing( int VCRank ) const EL_NO_EXCEPT;

    static int FindFactor( int p ) EL_NO_EXCEPT;

    friend bool GridCompare( const Grid & g1, const Grid & g2 );

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

bool operator==( const Grid& A, const Grid& B ) EL_NO_EXCEPT;
bool operator!=( const Grid& A, const Grid& B ) EL_NO_EXCEPT;

// Return a grid constructed using mpi::COMM_WORLD.
const Grid& DefaultGrid() EL_NO_RELEASE_EXCEPT;

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

inline bool GridCompare
( const Grid & g1, const Grid & g2 )
{
    // haveViewers_
    if( g1.HaveViewers() != g2.HaveViewers() )
        return false;

    // height_
    if( g1.Height() != g2.Height() )
        return false;

    // size_
    if( g1.Size() != g2.Size() )
        return false;

    // gcd_
    if( g1.GCD() != g2.GCD() )
        return false;

    // order_
    if( g1.Order() != g2.Order() )
        return false;

    // diagsAndRanks_
    if( g1.diagsAndRanks_.size() != g2.diagsAndRanks_.size() )
        return false;
    size_t sizei = g1.diagsAndRanks_.size();
    for( size_t i=0; i<sizei; ++i )
        if( g1.diagsAndRanks_[i] != g2.diagsAndRanks_[i] )
            return false;

    // viewingComm_
    if( Congruent( g1.ViewingComm(), g2.ViewingComm() ) )
        return false;

    // viewingGroup_
    if( Congruent( g1.viewingGroup_, g2.viewingGroup_ ) )
        return false;

    // vcToViewing_
    if( g1.vcToViewing_.size() != g2.vcToViewing_.size() )
        return false;
    sizei = g1.vcToViewing_.size();
    for( size_t i=0; i<sizei; ++i )
        if( g1.vcToViewing_[i] != g2.vcToViewing_[i] )
            return false;

    // owningComm_
    if( Congruent( g1.OwningComm(), g2.OwningComm() ) )
        return false;

    // owningGroup_
    if( Congruent( g1.OwningGroup(), g2.OwningGroup() ) )
        return false;

    // cartComm_
    if( Congruent( g1.cartComm_, g2.cartComm_ ) )
        return false;

    // mcComm_
    if( Congruent( g1.MCComm(), g2.MCComm() ) )
        return false;

    // mrComm_
    if( Congruent( g1.MRComm(), g2.MRComm() ) )
        return false;

    // mdComm_
    if( Congruent( g1.MDComm(), g2.MDComm() ) )
        return false;

    // mdPerpComm_
    if( Congruent( g1.MDPerpComm(), g2.MDPerpComm() ) )
        return false;

    // vcComm_
    if( Congruent( g1.VCComm(), g2.VCComm() ) )
        return false;

    // vrComm_
    if( Congruent( g1.VRComm(), g2.VRComm() ) )
        return false;

    return true;
}

} // namespace El

#endif // ifndef EL_GRID_HPP
