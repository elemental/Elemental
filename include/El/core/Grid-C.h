/*
   Copyright (c) 2009-2014, Jack Poulson
                      2013, Jed Brown 
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_GRID_C_H
#define EL_GRID_C_H

#ifdef __cplusplus
extern "C" {
#endif

typedef       struct ElGridDummy* ElGrid;
typedef const struct ElGridDummy* ElConstGrid;

/* const Grid& DefaultGrid() 
   ------------------------- */
ElError ElDefaultGrid( ElConstGrid* grid );

/* Grid::Grid( mpi::Comm comm, GridOrder order ) 
   --------------------------------------------- */
ElError ElGridCreate( MPI_Comm comm, ElGridOrderType order, ElGrid* grid );

/* Grid::Grid( mpi::Comm comm, int height, GridOrder order )
   --------------------------------------------------------- */
ElError ElGridCreateSpecific
( MPI_Comm comm, int height, ElGridOrderType order, ElGrid* grid );

/* Grid::~Grid()
   ------------- */
ElError ElGridDestroy( ElConstGrid grid );

/* int Grid::Row() const 
   --------------------- */
ElError ElGridRow( ElConstGrid grid, int* row );

/* int Grid::Col() const
   --------------------- */
ElError ElGridCol( ElConstGrid grid, int* col );

/* int Grid::Rank() const
   ---------------------- */
ElError ElGridRank( ElConstGrid grid, int* rank );

/* int Grid::Height() const
   ------------------------ */
ElError ElGridHeight( ElConstGrid grid, int* height );

/* int Grid::Width() const
   ----------------------- */
ElError ElGridWidth( ElConstGrid grid, int* width );

/* int Grid::Size() const
   ---------------------- */
ElError ElGridSize( ElConstGrid grid, int* size );

/* GridOrder Grid::Order() const
   ----------------------------- */
ElError ElGridOrder( ElConstGrid grid, ElGridOrderType* order );

/* mpi::Comm Grid::ColComm() const 
   ------------------------------- */
ElError ElGridColComm( ElConstGrid grid, MPI_Comm* comm );

/* mpi::Comm Grid::RowComm() const
   ------------------------------- */
ElError ElGridRowComm( ElConstGrid grid, MPI_Comm* comm );

/* mpi::Comm Grid::Comm() const
   ---------------------------- */
ElError ElGridComm( ElConstGrid grid, MPI_Comm* comm );

/* int Grid::MCRank() const
   ------------------------ */
ElError ElGridMCRank( ElConstGrid grid, int* mcRank );

/* int Grid::MRRank() const
   ------------------------ */
ElError ElGridMRRank( ElConstGrid grid, int* mrRank );

/* int Grid::VCRank() const
   ------------------------ */
ElError ElGridVCRank( ElConstGrid grid, int* vcRank );

/* int Grid::VRRank() const
   ------------------------ */
ElError ElGridVRRank( ElConstGrid grid, int* vrRank );

/* int Grid::MCSize() const
   ------------------------ */
ElError ElGridMCSize( ElConstGrid grid, int* mcSize );

/* int Grid::MRSize() const
   ------------------------ */
ElError ElGridMRSize( ElConstGrid grid, int* mrSize );

/* int Grid::VCSize() const
   ------------------------ */
ElError ElGridVCSize( ElConstGrid grid, int* vcSize );

/* int Grid::VRSize() const
   ------------------------ */
ElError ElGridVRSize( ElConstGrid grid, int* vrSize );

/* mpi::Comm Grid::MCComm() const
   ------------------------------ */
ElError ElGridMCComm( ElConstGrid grid, MPI_Comm* comm );

/* mpi::Comm Grid::MRComm() const
   ------------------------------ */
ElError ElGridMRComm( ElConstGrid grid, MPI_Comm* comm );

/* mpi::Comm Grid::VCComm() const
   ------------------------------ */
ElError ElGridVCComm( ElConstGrid grid, MPI_Comm* comm );

/* mpi::Comm Grid::VRComm() const
   ------------------------------ */
ElError ElGridVRComm( ElConstGrid grid, MPI_Comm* comm );

/* mpi::Comm Grid::MDComm() const
   ------------------------------ */
ElError ElGridMDComm( ElConstGrid grid, MPI_Comm* comm );

/* mpi::Comm Grid::MDPerpComm() const 
   ---------------------------------- */
ElError ElGridMDPerpComm( ElConstGrid grid, MPI_Comm* comm );

/* Grid::Grid( mpi::Comm comm, mpi::Group owners, int height, GridOrder order )
   ----------------------------------------------------------------------------
*/
ElError ElGridCreateAdvanced
( MPI_Comm comm, MPI_Group owners, int height, ElGridOrderType order, 
  ElGrid* grid );

/* int Grid::LCD() const
   --------------------- */
ElError ElGridGCD( ElConstGrid grid, int* gcd );

/* int Grid::LCM() const
   --------------------- */
ElError ElGridLCM( ElConstGrid grid, int* lcm );

/* bool Grid::InGrid() const
   ------------------------- */
ElError ElGridInGrid( ElConstGrid grid, bool* inGrid );

/* bool Grid::HaveViewers() const
   ------------------------------ */
ElError ElGridHaveViewers( ElConstGrid grid, bool* haveViewers );

/* int Grid::OwningRank() const
   ---------------------------- */
ElError ElGridOwningRank( ElConstGrid grid, int* owningRank );

/* int Grid::ViewingRank() const
   ----------------------------- */
ElError ElGridViewingRank( ElConstGrid grid, int* viewingRank );

/* int Grid::VCToViewingMap( int vcRank ) const 
   -------------------------------------------- */
ElError ElGridVCToViewingMap( ElConstGrid grid, int VCRank, int* viewingRank );

/* mpi::Group Grid::OwningGroup() const
   ------------------------------------ */
ElError ElGridOwningGroup( ElConstGrid grid, MPI_Group* group );

/* mpi::Comm Grid::OwningComm() const
   ---------------------------------- */
ElError ElGridOwningComm( ElConstGrid grid, MPI_Comm* comm );

/* mpi::Comm Grid::ViewingComm() const
   ----------------------------------- */
ElError ElGridViewingComm( ElConstGrid grid, MPI_Comm* comm );

/* int Grid::DiagPath( int vcRank ) const
   -------------------------------------- */
ElError ElGridDiagPath( ElConstGrid grid, int VCRank, int* diagPath );

/* int Grid::DiagPathRank( int vcRank ) const
   ------------------------------------------ */
ElError ElGridDiagPathRank( ElConstGrid grid, int VCRank, int* diagPathRank );

/* int Grid::FirstVCRank( int diagPath ) const
   ------------------------------------------- */
ElError ElGridFirstVCRank( ElConstGrid grid, int diagPath, int* firstVCRank );

/* static int Grid::FindFactor( int p )
   ------------------------------------ */
ElError ElGridFindFactor( int p, int* factor );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_GRID_C_H */
