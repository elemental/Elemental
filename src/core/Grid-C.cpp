/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>
#include <El-lite.h>
using namespace El;

extern "C" {

ElError ElDefaultGrid( ElConstGrid* grid )
{ EL_TRY( *grid = CReflect(&Grid::Default()) ) }

ElError ElTrivialGrid( ElConstGrid* grid )
{ EL_TRY( *grid = CReflect(&Grid::Trivial()) ) }

// Grid::Grid( MPI_Comm comm, GridOrder order )
// --------------------------------------------
ElError ElGridCreate
( MPI_Comm comm, ElGridOrderType order, ElGrid* grid )
{ EL_TRY( *grid = CReflect(new Grid(comm,CReflect(order))) ) }

// Grid::Grid( MPI_Comm comm, int height, GridOrder order )
// --------------------------------------------------------
ElError ElGridCreateSpecific
( MPI_Comm comm, int height, ElGridOrderType order, ElGrid* grid )
{ EL_TRY( *grid = CReflect(new Grid(comm,height,CReflect(order))) ) }

// Grid::~Grid()
// -------------
ElError ElGridDestroy( ElConstGrid grid )
{ EL_TRY( delete CReflect(grid) ) }

// int Grid::Row() const
// ---------------------
ElError ElGridRow( ElConstGrid grid, int* row )
{ EL_TRY( ElGridMCRank(grid,row) ) }

// int Grid::Col() const
// ---------------------
ElError ElGridCol( ElConstGrid grid, int* col )
{ EL_TRY( ElGridMRRank(grid,col) ) }

// int Grid::Rank() const
// ----------------------
ElError ElGridRank( ElConstGrid grid, int* rank )
{ EL_TRY( ElGridVCRank(grid,rank) ) }

// int Grid::Height() const
// ------------------------
ElError ElGridHeight( ElConstGrid grid, int* height )
{ EL_TRY( *height = CReflect(grid)->Height() ) }

// int Grid::Width() const
// -----------------------
ElError ElGridWidth( ElConstGrid grid, int* width )
{ EL_TRY( *width = CReflect(grid)->Width() ) }

// int Grid::Size() const
// ----------------------
ElError ElGridSize( ElConstGrid grid, int* size )
{ EL_TRY( *size = CReflect(grid)->Size() ) }

// GridOrder Grid::Order() const
// -----------------------------
ElError ElGridOrder( ElConstGrid grid, ElGridOrderType* order )
{ EL_TRY( *order = CReflect(CReflect(grid)->Order()) ) }

// mpi::Comm Grid::ColComm() const
// -------------------------------
ElError ElGridColComm( ElConstGrid grid, MPI_Comm* comm )
{ EL_TRY( *comm = CReflect(grid)->ColComm().comm ) }

// mpi::Comm Grid::RowComm() const
// -------------------------------
ElError ElGridRowComm( ElConstGrid grid, MPI_Comm* comm )
{ EL_TRY( *comm = CReflect(grid)->RowComm().comm ) }

// mpi::Comm Grid::Comm() const
// ----------------------------
ElError ElGridComm( ElConstGrid grid, MPI_Comm* comm )
{ EL_TRY( *comm = CReflect(grid)->Comm().comm ) }

// int Grid::MCRank() const
// ------------------------
ElError ElGridMCRank( ElConstGrid grid, int* mcRank )
{ EL_TRY( *mcRank = CReflect(grid)->MCRank() ) }

// int Grid::MRRank() const
// ------------------------
ElError ElGridMRRank( ElConstGrid grid, int* mrRank )
{ EL_TRY( *mrRank = CReflect(grid)->MRRank() ) }

// int Grid::VCRank() const
// ------------------------
ElError ElGridVCRank( ElConstGrid grid, int* vcRank )
{ EL_TRY( *vcRank = CReflect(grid)->VCRank() ) }

// int Grid::VRRank() const
// ------------------------
ElError ElGridVRRank( ElConstGrid grid, int* vrRank )
{ EL_TRY( *vrRank = CReflect(grid)->VRRank() ) }

// int Grid::MCSize() const
// ------------------------
ElError ElGridMCSize( ElConstGrid grid, int* mcSize )
{ EL_TRY( *mcSize = CReflect(grid)->MCSize() ) }

// int Grid::MRSize() const
// ------------------------
ElError ElGridMRSize( ElConstGrid grid, int* mrSize )
{ EL_TRY( *mrSize = CReflect(grid)->MRSize() ) }

// int Grid::VCSize() const
// ------------------------
ElError ElGridVCSize( ElConstGrid grid, int* vcSize )
{ EL_TRY( *vcSize = CReflect(grid)->VCSize() ) }

// int Grid::VRSize() const
// ------------------------
ElError ElGridVRSize( ElConstGrid grid, int* vrSize )
{ EL_TRY( *vrSize = CReflect(grid)->VRSize() ) }

// mpi::Comm Grid::MCComm() const
// ------------------------------
ElError ElGridMCComm( ElConstGrid grid, MPI_Comm* comm )
{ EL_TRY( *comm = CReflect(grid)->MCComm().comm ) }

// mpi::Comm Grid::MRComm() const
// ------------------------------
ElError ElGridMRComm( ElConstGrid grid, MPI_Comm* comm )
{ EL_TRY( *comm = CReflect(grid)->MRComm().comm ) }

// mpi::Comm Grid::VCComm() const
// ------------------------------
ElError ElGridVCComm( ElConstGrid grid, MPI_Comm* comm )
{ EL_TRY( *comm = CReflect(grid)->VCComm().comm ) }

// mpi::Comm Grid::VRComm() const
// ------------------------------
ElError ElGridVRComm( ElConstGrid grid, MPI_Comm* comm )
{ EL_TRY( *comm = CReflect(grid)->VRComm().comm ) }

// mpi::Comm Grid::MDComm() const
// ------------------------------
ElError ElGridMDComm( ElConstGrid grid, MPI_Comm* comm )
{ EL_TRY( *comm = CReflect(grid)->MDComm().comm ) }

// mpi::Comm Grid::MDPerpComm() const
// ----------------------------------
ElError ElGridMDPerpComm( ElConstGrid grid, MPI_Comm* comm )
{ EL_TRY( *comm = CReflect(grid)->MDPerpComm().comm ) }

// Grid::Grid( mpi::Comm comm, mpi::Group owners, int height, GridOrder order )
// ----------------------------------------------------------------------------
ElError ElGridCreateAdvanced
( MPI_Comm comm, MPI_Group owners, int height, ElGridOrderType order,
  ElGrid* grid )
{ EL_TRY( *grid = CReflect(new Grid(comm,owners,height,CReflect(order))) ) }

// int Grid::GCD() const
// ---------------------
ElError ElGridGCD( ElConstGrid grid, int* gcd )
{ EL_TRY( *gcd = CReflect(grid)->GCD() ) }

// int Grid::LCM() const
// ---------------------
ElError ElGridLCM( ElConstGrid grid, int* lcm )
{ EL_TRY( *lcm = CReflect(grid)->LCM() ) }

// bool Grid::InGrid() const
// -------------------------
ElError ElGridInGrid( ElConstGrid grid, bool* inGrid )
{ EL_TRY( *inGrid = CReflect(grid)->InGrid() ) }

// bool Grid::HaveViewers() const
// ------------------------------
ElError ElGridHaveViewers( ElConstGrid grid, bool* haveViewers )
{ EL_TRY( *haveViewers = CReflect(grid)->HaveViewers() ) }

// int Grid::OwningRank() const
// ----------------------------
ElError ElGridOwningRank( ElConstGrid grid, int* owningRank )
{ EL_TRY( *owningRank = CReflect(grid)->OwningRank() ) }

// int Grid::ViewingRank() const
// -----------------------------
ElError ElGridViewingRank( ElConstGrid grid, int* viewingRank )
{ EL_TRY( *viewingRank = CReflect(grid)->ViewingRank() ) }

// int Grid::VCToViewing( int vcRank ) const
// -----------------------------------------
ElError ElGridVCToViewing( ElConstGrid grid, int vcRank, int* viewingRank )
{ EL_TRY( *viewingRank = CReflect(grid)->VCToViewing(vcRank) ) }

// mpi::Group Grid::OwningGroup() const
// ------------------------------------
ElError ElGridOwningGroup( ElConstGrid grid, MPI_Group* group )
{ EL_TRY( *group = CReflect(grid)->OwningGroup().group ) }

// mpi::Comm Grid::OwningComm() const
// ----------------------------------
ElError ElGridOwningComm( ElConstGrid grid, MPI_Comm* comm )
{ EL_TRY( *comm = CReflect(grid)->OwningComm().comm ) }

// mpi::Comm Grid::ViewingComm() const
// -----------------------------------
ElError ElGridViewingComm( ElConstGrid grid, MPI_Comm* comm )
{ EL_TRY( *comm = CReflect(grid)->ViewingComm().comm ) }

// int Grid::Diag( int vcRank ) const
// ----------------------------------
ElError ElGridDiag( ElConstGrid grid, int vcRank, int* diag )
{ EL_TRY( *diag = CReflect(grid)->Diag(vcRank) ) }

// int Grid::DiagRank( int vcRank ) const
// --------------------------------------
ElError ElGridDiagRank( ElConstGrid grid, int vcRank, int* diagRank )
{ EL_TRY( *diagRank = CReflect(grid)->DiagRank(vcRank) ) }

// static int Grid::DefaultHeight( int gridSize )
// ----------------------------------------------
ElError ElGridDefaultHeight( int gridSize, int* gridHeight )
{ EL_TRY( *gridHeight = Grid::DefaultHeight(gridSize) ) }

} // extern "C"
