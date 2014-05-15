/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"
#include "El-C.h"
using namespace El;

#define RCG(gridHandle) reinterpret_cast<Grid*>(gridHandle)
#define RCG_const(gridHandle) reinterpret_cast<const Grid*>(gridHandle)

#define CATCH catch( std::exception& e ) { ReportException(e); }

extern "C" {

const ElGrid* ElDefaultGrid()
{ return reinterpret_cast<const ElGrid*>(&DefaultGrid()); }

// Grid::Grid( MPI_Comm comm, GridOrder order )
// --------------------------------------------
ElGrid* ElGridCreate( MPI_Comm comm, ElGridOrderType orderC )
{
    ElGrid* gridHandle = 0;
    GridOrder order = static_cast<GridOrder>(orderC);
    try { gridHandle = reinterpret_cast<ElGrid*>(new Grid(comm,order)); }
    CATCH
    return gridHandle;
}

// Grid::Grid( MPI_Comm comm, int height, GridOrder order )
// --------------------------------------------------------
ElGrid* ElGridCreateSpecific
( MPI_Comm comm, int height, ElGridOrderType orderC )
{
    ElGrid* gridHandle = 0;
    GridOrder order = static_cast<GridOrder>(orderC);
    try { gridHandle = reinterpret_cast<ElGrid*>(new Grid(comm,height,order)); }
    CATCH
    return gridHandle;
}

// Grid::~Grid()
// -------------
void ElGridDestroy( const ElGrid* gridHandle )
{ delete RCG_const(gridHandle); }

// int Grid::Row() const
// ---------------------
int ElGridRow( const ElGrid* gridHandle )
{ return ElGridMCRank(gridHandle); }

// int Grid::Col() const
// ---------------------
int ElGridCol( const ElGrid* gridHandle )
{ return ElGridMRRank(gridHandle); }

// int Grid::Rank() const
// ----------------------
int ElGridRank( const ElGrid* gridHandle )
{ return ElGridVCRank(gridHandle); }

// int Grid::Height() const
// ------------------------
int ElGridHeight( const ElGrid* gridHandle )
{ return RCG_const(gridHandle)->Height(); }

// int Grid::Width() const
// -----------------------
int ElGridWidth( const ElGrid* gridHandle )
{ return RCG_const(gridHandle)->Width(); }

// int Grid::Size() const
// ----------------------
int ElGridSize( const ElGrid* gridHandle )
{ return RCG_const(gridHandle)->Size(); }

// GridOrder Grid::Order() const
// -----------------------------
ElGridOrderType ElGridOrder( const ElGrid* gridHandle )
{ return static_cast<ElGridOrderType>(RCG_const(gridHandle)->Order()); }

// mpi::Comm Grid::ColComm() const
// -------------------------------
MPI_Comm ElGridColComm( const ElGrid* gridHandle )
{ return RCG_const(gridHandle)->ColComm().comm; }

// mpi::Comm Grid::RowComm() const
// -------------------------------
MPI_Comm ElGridRowComm( const ElGrid* gridHandle )
{ return RCG_const(gridHandle)->RowComm().comm; }

// mpi::Comm Grid::Comm() const
// ----------------------------
MPI_Comm ElGridComm( const ElGrid* gridHandle )
{ return RCG_const(gridHandle)->Comm().comm; }

// int Grid::MCRank() const
// ------------------------
int ElGridMCRank( const ElGrid* gridHandle )
{
    int rank = -1;
    try { rank = RCG_const(gridHandle)->MCRank(); }
    CATCH
    return rank;
}

// int Grid::MRRank() const
// ------------------------
int ElGridMRRank( const ElGrid* gridHandle )
{
    int rank = -1;
    try { rank = RCG_const(gridHandle)->MRRank(); }
    CATCH
    return rank;
}

// int Grid::VCRank() const
// ------------------------
int ElGridVCRank( const ElGrid* gridHandle )
{ 
    int rank = -1;
    try { rank = RCG_const(gridHandle)->VCRank(); }
    CATCH
    return rank;
}

// int Grid::VRRank() const
// ------------------------
int ElGridVRRank( const ElGrid* gridHandle )
{ 
    int rank = -1;
    try { rank = RCG_const(gridHandle)->VRRank(); }
    CATCH
    return rank;
}

// int Grid::MCSize() const
// ------------------------
int ElGridMCSize( const ElGrid* gridHandle )
{ return RCG_const(gridHandle)->MCSize(); }

// int Grid::MRSize() const
// ------------------------
int ElGridMRSize( const ElGrid* gridHandle )
{ return RCG_const(gridHandle)->MRSize(); }

// int Grid::VCSize( const ElGrid* gridHandle ) const
// --------------------------------------------------
int ElGridVCSize( const ElGrid* gridHandle )
{ return RCG_const(gridHandle)->VCSize(); }

// int Grid::VRSize( const ElGrid* gridHandle ) const
// --------------------------------------------------
int ElGridVRSize( const ElGrid* gridHandle )
{ return RCG_const(gridHandle)->VRSize(); }

// mpi::Comm Grid::MCComm() const
// ------------------------------
MPI_Comm ElGridMCComm( const ElGrid* gridHandle )
{ return RCG_const(gridHandle)->MCComm().comm; }

// mpi::Comm Grid::MRComm() const
// ------------------------------
MPI_Comm ElGridMRComm( const ElGrid* gridHandle )
{ return RCG_const(gridHandle)->MRComm().comm; }

// mpi::Comm Grid::VCComm() const
// ------------------------------
MPI_Comm ElGridVCComm( const ElGrid* gridHandle )
{ return RCG_const(gridHandle)->VCComm().comm; }

// mpi::Comm Grid::VRComm() const
// ------------------------------
MPI_Comm ElGridVRComm( const ElGrid* gridHandle )
{ return RCG_const(gridHandle)->VRComm().comm; }

// mpi::Comm Grid::MDComm() const
// ------------------------------
MPI_Comm ElGridMDComm( const ElGrid* gridHandle )
{ return RCG_const(gridHandle)->MDComm().comm; }

// mpi::Comm Grid::MDPerpComm() const
// ----------------------------------
MPI_Comm ElGridMDPerpComm( const ElGrid* gridHandle )
{ return RCG_const(gridHandle)->MDPerpComm().comm; }

// Grid::Grid( mpi::Comm comm, mpi::Group owners, int height, GridOrder order )
// ----------------------------------------------------------------------------
ElGrid* ElGridCreateAdvanced
( MPI_Comm comm, MPI_Group owners, int height, ElGridOrderType orderC )
{
    ElGrid* gridHandle = 0;
    GridOrder order = static_cast<GridOrder>(orderC);
    try { gridHandle = reinterpret_cast<ElGrid*>
                       (new Grid(comm,owners,height,order)); }
    CATCH
    return gridHandle;
}

// int Grid::GCD() const
// ---------------------
int ElGridGCD( const ElGrid* gridHandle )
{ return RCG_const(gridHandle)->GCD(); }

// int Grid::LCM() const
// ---------------------
int ElGridLCM( const ElGrid* gridHandle )
{ return RCG_const(gridHandle)->LCM(); }

// bool Grid::InGrid() const
// -------------------------
bool ElGridInGrid( const ElGrid* gridHandle )
{ return RCG_const(gridHandle)->InGrid(); }

// bool Grid::HaveViewers() const
// ------------------------------
bool ElGridHaveViewers( const ElGrid* gridHandle )
{ return RCG_const(gridHandle)->HaveViewers(); }

// int Grid::OwningRank() const
// ----------------------------
int ElGridOwningRank( const ElGrid* gridHandle )
{ return RCG_const(gridHandle)->OwningRank(); }

// int Grid::ViewingRank() const
// -----------------------------
int ElGridViewingRank( const ElGrid* gridHandle )
{ return RCG_const(gridHandle)->ViewingRank(); }

// int Grid::VCToViewingMap( const ElGrid* grid, int vcRank ) const
// ----------------------------------------------------------------
int ElGridVCToViewingMap( const ElGrid* gridHandle, int vcRank )
{ return RCG_const(gridHandle)->VCToViewingMap(vcRank); }

// mpi::Group Grid::OwningGroup() const
// ------------------------------------
MPI_Group ElGridOwningGroup( const ElGrid* gridHandle )
{ return RCG_const(gridHandle)->OwningGroup().group; }

// mpi::Comm Grid::OwningComm() const
// ----------------------------------
MPI_Comm ElGridOwningComm( const ElGrid* gridHandle )
{ return RCG_const(gridHandle)->OwningComm().comm; }

// mpi::Comm Grid::ViewingComm() const
// -----------------------------------
MPI_Comm ElGridViewingComm( const ElGrid* gridHandle )
{ return RCG_const(gridHandle)->ViewingComm().comm; }

// int Grid::DiagPath( int vcRank ) const
// --------------------------------------
int ElGridDiagPath( const ElGrid* gridHandle, int vcRank )
{ return RCG_const(gridHandle)->DiagPath(vcRank); }

// int Grid::DiagPathRank( int vcRank ) const
// ------------------------------------------
int ElGridDiagPathRank( const ElGrid* gridHandle, int vcRank )
{ return RCG_const(gridHandle)->DiagPathRank(vcRank); }

// int Grid::FirstVCRank( int diagPath ) const
// -------------------------------------------
int ElGridFirstVCRank( const ElGrid* gridHandle, int vcRank )
{ return RCG_const(gridHandle)->FirstVCRank(vcRank); }

// static int Grid::FindFactor( int p )
// ------------------------------------
int ElGridFindFactor( int p )
{ return Grid::FindFactor(p); }

} // extern "C"
