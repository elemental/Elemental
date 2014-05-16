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

ElConstGrid ElDefaultGrid()
{ return (ElConstGrid)
         reinterpret_cast<const struct ElGridDummy*>(&DefaultGrid()); }

// Grid::Grid( MPI_Comm comm, GridOrder order )
// --------------------------------------------
ElGrid ElGridCreate( MPI_Comm comm, ElGridOrderType orderC )
{
    ElGrid gridHandle = 0;
    GridOrder order = static_cast<GridOrder>(orderC);
    try { gridHandle = reinterpret_cast<ElGrid>(new Grid(comm,order)); }
    CATCH
    return gridHandle;
}

// Grid::Grid( MPI_Comm comm, int height, GridOrder order )
// --------------------------------------------------------
ElGrid ElGridCreateSpecific
( MPI_Comm comm, int height, ElGridOrderType orderC )
{
    ElGrid gridHandle = 0;
    GridOrder order = static_cast<GridOrder>(orderC);
    try { gridHandle = reinterpret_cast<ElGrid>(new Grid(comm,height,order)); }
    CATCH
    return gridHandle;
}

// Grid::~Grid()
// -------------
void ElGridDestroy( ElConstGrid gridHandle )
{ delete RCG_const(gridHandle); }

// int Grid::Row() const
// ---------------------
int ElGridRow( ElConstGrid gridHandle )
{ return ElGridMCRank(gridHandle); }

// int Grid::Col() const
// ---------------------
int ElGridCol( ElConstGrid gridHandle )
{ return ElGridMRRank(gridHandle); }

// int Grid::Rank() const
// ----------------------
int ElGridRank( ElConstGrid gridHandle )
{ return ElGridVCRank(gridHandle); }

// int Grid::Height() const
// ------------------------
int ElGridHeight( ElConstGrid gridHandle )
{ return RCG_const(gridHandle)->Height(); }

// int Grid::Width() const
// -----------------------
int ElGridWidth( ElConstGrid gridHandle )
{ return RCG_const(gridHandle)->Width(); }

// int Grid::Size() const
// ----------------------
int ElGridSize( ElConstGrid gridHandle )
{ return RCG_const(gridHandle)->Size(); }

// GridOrder Grid::Order() const
// -----------------------------
ElGridOrderType ElGridOrder( ElConstGrid gridHandle )
{ return static_cast<ElGridOrderType>(RCG_const(gridHandle)->Order()); }

// mpi::Comm Grid::ColComm() const
// -------------------------------
MPI_Comm ElGridColComm( ElConstGrid gridHandle )
{ return RCG_const(gridHandle)->ColComm().comm; }

// mpi::Comm Grid::RowComm() const
// -------------------------------
MPI_Comm ElGridRowComm( ElConstGrid gridHandle )
{ return RCG_const(gridHandle)->RowComm().comm; }

// mpi::Comm Grid::Comm() const
// ----------------------------
MPI_Comm ElGridComm( ElConstGrid gridHandle )
{ return RCG_const(gridHandle)->Comm().comm; }

// int Grid::MCRank() const
// ------------------------
int ElGridMCRank( ElConstGrid gridHandle )
{
    int rank = -1;
    try { rank = RCG_const(gridHandle)->MCRank(); }
    CATCH
    return rank;
}

// int Grid::MRRank() const
// ------------------------
int ElGridMRRank( ElConstGrid gridHandle )
{
    int rank = -1;
    try { rank = RCG_const(gridHandle)->MRRank(); }
    CATCH
    return rank;
}

// int Grid::VCRank() const
// ------------------------
int ElGridVCRank( ElConstGrid gridHandle )
{ 
    int rank = -1;
    try { rank = RCG_const(gridHandle)->VCRank(); }
    CATCH
    return rank;
}

// int Grid::VRRank() const
// ------------------------
int ElGridVRRank( ElConstGrid gridHandle )
{ 
    int rank = -1;
    try { rank = RCG_const(gridHandle)->VRRank(); }
    CATCH
    return rank;
}

// int Grid::MCSize() const
// ------------------------
int ElGridMCSize( ElConstGrid gridHandle )
{ return RCG_const(gridHandle)->MCSize(); }

// int Grid::MRSize() const
// ------------------------
int ElGridMRSize( ElConstGrid gridHandle )
{ return RCG_const(gridHandle)->MRSize(); }

// int Grid::VCSize() const
// ------------------------
int ElGridVCSize( ElConstGrid gridHandle )
{ return RCG_const(gridHandle)->VCSize(); }

// int Grid::VRSize() const
// ------------------------
int ElGridVRSize( ElConstGrid gridHandle )
{ return RCG_const(gridHandle)->VRSize(); }

// mpi::Comm Grid::MCComm() const
// ------------------------------
MPI_Comm ElGridMCComm( ElConstGrid gridHandle )
{ return RCG_const(gridHandle)->MCComm().comm; }

// mpi::Comm Grid::MRComm() const
// ------------------------------
MPI_Comm ElGridMRComm( ElConstGrid gridHandle )
{ return RCG_const(gridHandle)->MRComm().comm; }

// mpi::Comm Grid::VCComm() const
// ------------------------------
MPI_Comm ElGridVCComm( ElConstGrid gridHandle )
{ return RCG_const(gridHandle)->VCComm().comm; }

// mpi::Comm Grid::VRComm() const
// ------------------------------
MPI_Comm ElGridVRComm( ElConstGrid gridHandle )
{ return RCG_const(gridHandle)->VRComm().comm; }

// mpi::Comm Grid::MDComm() const
// ------------------------------
MPI_Comm ElGridMDComm( ElConstGrid gridHandle )
{ return RCG_const(gridHandle)->MDComm().comm; }

// mpi::Comm Grid::MDPerpComm() const
// ----------------------------------
MPI_Comm ElGridMDPerpComm( ElConstGrid gridHandle )
{ return RCG_const(gridHandle)->MDPerpComm().comm; }

// Grid::Grid( mpi::Comm comm, mpi::Group owners, int height, GridOrder order )
// ----------------------------------------------------------------------------
ElGrid ElGridCreateAdvanced
( MPI_Comm comm, MPI_Group owners, int height, ElGridOrderType orderC )
{
    ElGrid gridHandle = 0;
    GridOrder order = static_cast<GridOrder>(orderC);
    try { gridHandle = reinterpret_cast<ElGrid>
                       (new Grid(comm,owners,height,order)); }
    CATCH
    return gridHandle;
}

// int Grid::GCD() const
// ---------------------
int ElGridGCD( ElConstGrid gridHandle )
{ return RCG_const(gridHandle)->GCD(); }

// int Grid::LCM() const
// ---------------------
int ElGridLCM( ElConstGrid gridHandle )
{ return RCG_const(gridHandle)->LCM(); }

// bool Grid::InGrid() const
// -------------------------
bool ElGridInGrid( ElConstGrid gridHandle )
{ return RCG_const(gridHandle)->InGrid(); }

// bool Grid::HaveViewers() const
// ------------------------------
bool ElGridHaveViewers( ElConstGrid gridHandle )
{ return RCG_const(gridHandle)->HaveViewers(); }

// int Grid::OwningRank() const
// ----------------------------
int ElGridOwningRank( ElConstGrid gridHandle )
{ return RCG_const(gridHandle)->OwningRank(); }

// int Grid::ViewingRank() const
// -----------------------------
int ElGridViewingRank( ElConstGrid gridHandle )
{ return RCG_const(gridHandle)->ViewingRank(); }

// int Grid::VCToViewingMap( int vcRank ) const
// --------------------------------------------
int ElGridVCToViewingMap( ElConstGrid gridHandle, int vcRank )
{ return RCG_const(gridHandle)->VCToViewingMap(vcRank); }

// mpi::Group Grid::OwningGroup() const
// ------------------------------------
MPI_Group ElGridOwningGroup( ElConstGrid gridHandle )
{ return RCG_const(gridHandle)->OwningGroup().group; }

// mpi::Comm Grid::OwningComm() const
// ----------------------------------
MPI_Comm ElGridOwningComm( ElConstGrid gridHandle )
{ return RCG_const(gridHandle)->OwningComm().comm; }

// mpi::Comm Grid::ViewingComm() const
// -----------------------------------
MPI_Comm ElGridViewingComm( ElConstGrid gridHandle )
{ return RCG_const(gridHandle)->ViewingComm().comm; }

// int Grid::DiagPath( int vcRank ) const
// --------------------------------------
int ElGridDiagPath( ElConstGrid gridHandle, int vcRank )
{ return RCG_const(gridHandle)->DiagPath(vcRank); }

// int Grid::DiagPathRank( int vcRank ) const
// ------------------------------------------
int ElGridDiagPathRank( ElConstGrid gridHandle, int vcRank )
{ return RCG_const(gridHandle)->DiagPathRank(vcRank); }

// int Grid::FirstVCRank( int diagPath ) const
// -------------------------------------------
int ElGridFirstVCRank( ElConstGrid gridHandle, int vcRank )
{ return RCG_const(gridHandle)->FirstVCRank(vcRank); }

// static int Grid::FindFactor( int p )
// ------------------------------------
int ElGridFindFactor( int p )
{ return Grid::FindFactor(p); }

} // extern "C"
