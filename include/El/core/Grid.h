/*
   Copyright (c) 2009-2016, Jack Poulson
                      2013, Jed Brown 
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_GRID_C_H
#define EL_GRID_C_H

#ifdef __cplusplus
extern "C" {
#endif

typedef       struct ElGridDummy* ElGrid;
typedef const struct ElGridDummy* ElConstGrid;

/* const Grid& DefaultGrid() 
   ------------------------- */
EL_EXPORT ElError ElDefaultGrid( ElConstGrid* grid );

/* Grid::Grid( mpi::Comm comm, GridOrder order ) 
   --------------------------------------------- */
EL_EXPORT ElError ElGridCreate
( MPI_Comm comm, ElGridOrderType order, ElGrid* grid );

/* Grid::Grid( mpi::Comm comm, int height, GridOrder order )
   --------------------------------------------------------- */
EL_EXPORT ElError ElGridCreateSpecific
( MPI_Comm comm, int height, ElGridOrderType order, ElGrid* grid );

/* Grid::~Grid()
   ------------- */
EL_EXPORT ElError ElGridDestroy( ElConstGrid grid );

/* int Grid::Row() const 
   --------------------- */
EL_EXPORT ElError ElGridRow( ElConstGrid grid, int* row );

/* int Grid::Col() const
   --------------------- */
EL_EXPORT ElError ElGridCol( ElConstGrid grid, int* col );

/* int Grid::Rank() const
   ---------------------- */
EL_EXPORT ElError ElGridRank( ElConstGrid grid, int* rank );

/* int Grid::Height() const
   ------------------------ */
EL_EXPORT ElError ElGridHeight( ElConstGrid grid, int* height );

/* int Grid::Width() const
   ----------------------- */
EL_EXPORT ElError ElGridWidth( ElConstGrid grid, int* width );

/* int Grid::Size() const
   ---------------------- */
EL_EXPORT ElError ElGridSize( ElConstGrid grid, int* size );

/* GridOrder Grid::Order() const
   ----------------------------- */
EL_EXPORT ElError ElGridOrder( ElConstGrid grid, ElGridOrderType* order );

/* mpi::Comm Grid::ColComm() const 
   ------------------------------- */
EL_EXPORT ElError ElGridColComm( ElConstGrid grid, MPI_Comm* comm );

/* mpi::Comm Grid::RowComm() const
   ------------------------------- */
EL_EXPORT ElError ElGridRowComm( ElConstGrid grid, MPI_Comm* comm );

/* mpi::Comm Grid::Comm() const
   ---------------------------- */
EL_EXPORT ElError ElGridComm( ElConstGrid grid, MPI_Comm* comm );

/* int Grid::MCRank() const
   ------------------------ */
EL_EXPORT ElError ElGridMCRank( ElConstGrid grid, int* mcRank );

/* int Grid::MRRank() const
   ------------------------ */
EL_EXPORT ElError ElGridMRRank( ElConstGrid grid, int* mrRank );

/* int Grid::VCRank() const
   ------------------------ */
EL_EXPORT ElError ElGridVCRank( ElConstGrid grid, int* vcRank );

/* int Grid::VRRank() const
   ------------------------ */
EL_EXPORT ElError ElGridVRRank( ElConstGrid grid, int* vrRank );

/* int Grid::MCSize() const
   ------------------------ */
EL_EXPORT ElError ElGridMCSize( ElConstGrid grid, int* mcSize );

/* int Grid::MRSize() const
   ------------------------ */
EL_EXPORT ElError ElGridMRSize( ElConstGrid grid, int* mrSize );

/* int Grid::VCSize() const
   ------------------------ */
EL_EXPORT ElError ElGridVCSize( ElConstGrid grid, int* vcSize );

/* int Grid::VRSize() const
   ------------------------ */
EL_EXPORT ElError ElGridVRSize( ElConstGrid grid, int* vrSize );

/* mpi::Comm Grid::MCComm() const
   ------------------------------ */
EL_EXPORT ElError ElGridMCComm( ElConstGrid grid, MPI_Comm* comm );

/* mpi::Comm Grid::MRComm() const
   ------------------------------ */
EL_EXPORT ElError ElGridMRComm( ElConstGrid grid, MPI_Comm* comm );

/* mpi::Comm Grid::VCComm() const
   ------------------------------ */
EL_EXPORT ElError ElGridVCComm( ElConstGrid grid, MPI_Comm* comm );

/* mpi::Comm Grid::VRComm() const
   ------------------------------ */
EL_EXPORT ElError ElGridVRComm( ElConstGrid grid, MPI_Comm* comm );

/* mpi::Comm Grid::MDComm() const
   ------------------------------ */
EL_EXPORT ElError ElGridMDComm( ElConstGrid grid, MPI_Comm* comm );

/* mpi::Comm Grid::MDPerpComm() const 
   ---------------------------------- */
EL_EXPORT ElError ElGridMDPerpComm( ElConstGrid grid, MPI_Comm* comm );

/* Grid::Grid( mpi::Comm comm, mpi::Group owners, int height, GridOrder order )
   ----------------------------------------------------------------------------
*/
EL_EXPORT ElError ElGridCreateAdvanced
( MPI_Comm comm, MPI_Group owners, int height, ElGridOrderType order, 
  ElGrid* grid );

/* int Grid::GCD() const
   --------------------- */
EL_EXPORT ElError ElGridGCD( ElConstGrid grid, int* gcd );

/* int Grid::LCM() const
   --------------------- */
EL_EXPORT ElError ElGridLCM( ElConstGrid grid, int* lcm );

/* bool Grid::InGrid() const
   ------------------------- */
EL_EXPORT ElError ElGridInGrid( ElConstGrid grid, bool* inGrid );

/* bool Grid::HaveViewers() const
   ------------------------------ */
EL_EXPORT ElError ElGridHaveViewers( ElConstGrid grid, bool* haveViewers );

/* int Grid::OwningRank() const
   ---------------------------- */
EL_EXPORT ElError ElGridOwningRank( ElConstGrid grid, int* owningRank );

/* int Grid::ViewingRank() const
   ----------------------------- */
EL_EXPORT ElError ElGridViewingRank( ElConstGrid grid, int* viewingRank );

/* mpi::Group Grid::OwningGroup() const
   ------------------------------------ */
EL_EXPORT ElError ElGridOwningGroup( ElConstGrid grid, MPI_Group* group );

/* mpi::Comm Grid::OwningComm() const
   ---------------------------------- */
EL_EXPORT ElError ElGridOwningComm( ElConstGrid grid, MPI_Comm* comm );

/* mpi::Comm Grid::ViewingComm() const
   ----------------------------------- */
EL_EXPORT ElError ElGridViewingComm( ElConstGrid grid, MPI_Comm* comm );

/* int Grid::Diag( int vcRank ) const
   ---------------------------------- */
EL_EXPORT ElError ElGridDiag( ElConstGrid grid, int VCRank, int* diag );

/* int Grid::DiagRank( int vcRank ) const
   -------------------------------------- */
EL_EXPORT ElError ElGridDiagRank( ElConstGrid grid, int VCRank, int* diagRank );

/* TODO: VCToVR, VRToVC, CoordsToVC */
/* int Grid::VCToViewing( int vcRank ) const 
   ----------------------------------------- */
EL_EXPORT ElError ElGridVCToViewing
( ElConstGrid grid, int VCRank, int* viewingRank );

/* static int Grid::FindFactor( int p )
   ------------------------------------ */
EL_EXPORT ElError ElGridFindFactor( int p, int* factor );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_GRID_C_H */
