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

#define CATCH \
  catch( std::bad_alloc& e ) \
  { ReportException(e); return EL_ALLOC_ERROR; } \
  catch( ArgException& e ) \
  { ReportException(e); return EL_ARG_ERROR; } \
  catch( std::logic_error& e ) \
  { ReportException(e); return EL_LOGIC_ERROR; } \
  catch( std::runtime_error& e ) \
  { ReportException(e); return EL_RUNTIME_ERROR; } \
  catch( std::exception& e ) \
  { ReportException(e); return EL_ERROR; }

extern "C" {

ElError ElDefaultGrid
( ElConstGrid* gridHandle )
{ 
    *gridHandle = (ElConstGrid)
      reinterpret_cast<const struct ElGridDummy*>(&DefaultGrid()); 
    return EL_SUCCESS;
}

// Grid::Grid( MPI_Comm comm, GridOrder order )
// --------------------------------------------
ElError ElGridCreate
( MPI_Comm comm, ElGridOrderType orderC, ElGrid* gridHandle )
{
    GridOrder order = static_cast<GridOrder>(orderC);
    try { *gridHandle = reinterpret_cast<ElGrid>(new Grid(comm,order)); }
    CATCH
    return EL_SUCCESS;
}

// Grid::Grid( MPI_Comm comm, int height, GridOrder order )
// --------------------------------------------------------
ElError ElGridCreateSpecific
( MPI_Comm comm, int height, ElGridOrderType orderC, ElGrid* gridHandle )
{
    GridOrder order = static_cast<GridOrder>(orderC);
    try { *gridHandle = reinterpret_cast<ElGrid>(new Grid(comm,height,order)); }
    CATCH
    return EL_SUCCESS;
}

// Grid::~Grid()
// -------------
ElError ElGridDestroy( ElConstGrid gridHandle )
{ 
    delete RCG_const(gridHandle); 
    return EL_SUCCESS;
}

// int Grid::Row() const
// ---------------------
ElError ElGridRow( ElConstGrid gridHandle, int* row )
{ return ElGridMCRank(gridHandle,row); }

// int Grid::Col() const
// ---------------------
ElError ElGridCol( ElConstGrid gridHandle, int* col )
{ return ElGridMRRank(gridHandle,col); }

// int Grid::Rank() const
// ----------------------
ElError ElGridRank( ElConstGrid gridHandle, int* rank )
{ return ElGridVCRank(gridHandle,rank); }

// int Grid::Height() const
// ------------------------
ElError ElGridHeight( ElConstGrid gridHandle, int* height )
{ 
    *height = RCG_const(gridHandle)->Height(); 
    return EL_SUCCESS;
}

// int Grid::Width() const
// -----------------------
ElError ElGridWidth( ElConstGrid gridHandle, int* width )
{ 
    *width = RCG_const(gridHandle)->Width(); 
    return EL_SUCCESS;
}

// int Grid::Size() const
// ----------------------
ElError ElGridSize( ElConstGrid gridHandle, int* size )
{ 
    *size = RCG_const(gridHandle)->Size(); 
    return EL_SUCCESS;
}

// GridOrder Grid::Order() const
// -----------------------------
ElError ElGridOrder( ElConstGrid gridHandle, ElGridOrderType* order )
{ 
    *order = static_cast<ElGridOrderType>(RCG_const(gridHandle)->Order()); 
    return EL_SUCCESS;
}

// mpi::Comm Grid::ColComm() const
// -------------------------------
ElError ElGridColComm( ElConstGrid gridHandle, MPI_Comm* comm )
{ 
    *comm = RCG_const(gridHandle)->ColComm().comm; 
    return EL_SUCCESS;
}

// mpi::Comm Grid::RowComm() const
// -------------------------------
ElError ElGridRowComm( ElConstGrid gridHandle, MPI_Comm* comm )
{ 
    *comm = RCG_const(gridHandle)->RowComm().comm; 
    return EL_SUCCESS;
}

// mpi::Comm Grid::Comm() const
// ----------------------------
ElError ElGridComm( ElConstGrid gridHandle, MPI_Comm* comm )
{ 
    *comm = RCG_const(gridHandle)->Comm().comm; 
    return EL_SUCCESS;
}

// int Grid::MCRank() const
// ------------------------
ElError ElGridMCRank( ElConstGrid gridHandle, int* mcRank )
{
    try { *mcRank = RCG_const(gridHandle)->MCRank(); }
    CATCH
    return EL_SUCCESS;
}

// int Grid::MRRank() const
// ------------------------
ElError ElGridMRRank( ElConstGrid gridHandle, int* mrRank )
{
    try { *mrRank = RCG_const(gridHandle)->MRRank(); }
    CATCH
    return EL_SUCCESS;
}

// int Grid::VCRank() const
// ------------------------
ElError ElGridVCRank( ElConstGrid gridHandle, int* vcRank )
{ 
    try { *vcRank = RCG_const(gridHandle)->VCRank(); }
    CATCH
    return EL_SUCCESS;
}

// int Grid::VRRank() const
// ------------------------
ElError ElGridVRRank( ElConstGrid gridHandle, int* vrRank )
{ 
    try { *vrRank = RCG_const(gridHandle)->VRRank(); }
    CATCH
    return EL_SUCCESS;
}

// int Grid::MCSize() const
// ------------------------
ElError ElGridMCSize( ElConstGrid gridHandle, int* mcSize )
{ 
    *mcSize = RCG_const(gridHandle)->MCSize(); 
    return EL_SUCCESS;
}

// int Grid::MRSize() const
// ------------------------
ElError ElGridMRSize( ElConstGrid gridHandle, int* mrSize )
{ 
    *mrSize = RCG_const(gridHandle)->MRSize(); 
    return EL_SUCCESS;
}

// int Grid::VCSize() const
// ------------------------
ElError ElGridVCSize( ElConstGrid gridHandle, int* vcSize )
{ 
    *vcSize = RCG_const(gridHandle)->VCSize(); 
    return EL_SUCCESS;
}

// int Grid::VRSize() const
// ------------------------
ElError ElGridVRSize( ElConstGrid gridHandle, int* vrSize )
{ 
    *vrSize = RCG_const(gridHandle)->VRSize(); 
    return EL_SUCCESS;
}

// mpi::Comm Grid::MCComm() const
// ------------------------------
ElError ElGridMCComm( ElConstGrid gridHandle, MPI_Comm* comm )
{ 
    *comm = RCG_const(gridHandle)->MCComm().comm; 
    return EL_SUCCESS;
}

// mpi::Comm Grid::MRComm() const
// ------------------------------
ElError ElGridMRComm( ElConstGrid gridHandle, MPI_Comm* comm )
{ 
    *comm = RCG_const(gridHandle)->MRComm().comm; 
    return EL_SUCCESS;
}

// mpi::Comm Grid::VCComm() const
// ------------------------------
ElError ElGridVCComm( ElConstGrid gridHandle, MPI_Comm* comm )
{ 
    *comm = RCG_const(gridHandle)->VCComm().comm; 
    return EL_SUCCESS;
}

// mpi::Comm Grid::VRComm() const
// ------------------------------
ElError ElGridVRComm( ElConstGrid gridHandle, MPI_Comm* comm )
{ 
    *comm = RCG_const(gridHandle)->VRComm().comm; 
    return EL_SUCCESS;
}

// mpi::Comm Grid::MDComm() const
// ------------------------------
ElError ElGridMDComm( ElConstGrid gridHandle, MPI_Comm* comm )
{ 
    *comm = RCG_const(gridHandle)->MDComm().comm; 
    return EL_SUCCESS;
}

// mpi::Comm Grid::MDPerpComm() const
// ----------------------------------
ElError ElGridMDPerpComm( ElConstGrid gridHandle, MPI_Comm* comm )
{ 
    *comm = RCG_const(gridHandle)->MDPerpComm().comm; 
    return EL_SUCCESS;
}

// Grid::Grid( mpi::Comm comm, mpi::Group owners, int height, GridOrder order )
// ----------------------------------------------------------------------------
ElError ElGridCreateAdvanced
( MPI_Comm comm, MPI_Group owners, int height, ElGridOrderType orderC,
  ElGrid* gridHandle )
{
    GridOrder order = static_cast<GridOrder>(orderC);
    try { *gridHandle = reinterpret_cast<ElGrid>
                        (new Grid(comm,owners,height,order)); }
    CATCH
    return EL_SUCCESS;
}

// int Grid::GCD() const
// ---------------------
ElError ElGridGCD( ElConstGrid gridHandle, int* gcd )
{ 
    *gcd = RCG_const(gridHandle)->GCD(); 
    return EL_SUCCESS;
}

// int Grid::LCM() const
// ---------------------
ElError ElGridLCM( ElConstGrid gridHandle, int* lcm )
{ 
    *lcm = RCG_const(gridHandle)->LCM(); 
    return EL_SUCCESS;
}

// bool Grid::InGrid() const
// -------------------------
ElError ElGridInGrid( ElConstGrid gridHandle, bool* inGrid )
{ 
    *inGrid = RCG_const(gridHandle)->InGrid(); 
    return EL_SUCCESS;
}

// bool Grid::HaveViewers() const
// ------------------------------
ElError ElGridHaveViewers( ElConstGrid gridHandle, bool* haveViewers )
{ 
    *haveViewers = RCG_const(gridHandle)->HaveViewers(); 
    return EL_SUCCESS;
}

// int Grid::OwningRank() const
// ----------------------------
ElError ElGridOwningRank( ElConstGrid gridHandle, int* owningRank )
{ 
    *owningRank = RCG_const(gridHandle)->OwningRank(); 
    return EL_SUCCESS;
}

// int Grid::ViewingRank() const
// -----------------------------
ElError ElGridViewingRank( ElConstGrid gridHandle, int* viewingRank )
{ 
    *viewingRank = RCG_const(gridHandle)->ViewingRank(); 
    return EL_SUCCESS;
}

// int Grid::VCToViewingMap( int vcRank ) const
// --------------------------------------------
ElError ElGridVCToViewingMap
( ElConstGrid gridHandle, int vcRank, int* viewingRank )
{ 
    *viewingRank = RCG_const(gridHandle)->VCToViewingMap(vcRank); 
    return EL_SUCCESS;
}

// mpi::Group Grid::OwningGroup() const
// ------------------------------------
ElError ElGridOwningGroup( ElConstGrid gridHandle, MPI_Group* group )
{ 
    *group = RCG_const(gridHandle)->OwningGroup().group; 
    return EL_SUCCESS;
}

// mpi::Comm Grid::OwningComm() const
// ----------------------------------
ElError ElGridOwningComm( ElConstGrid gridHandle, MPI_Comm* comm )
{ 
    *comm = RCG_const(gridHandle)->OwningComm().comm; 
    return EL_SUCCESS;
}

// mpi::Comm Grid::ViewingComm() const
// -----------------------------------
ElError ElGridViewingComm( ElConstGrid gridHandle, MPI_Comm* comm )
{ 
    *comm = RCG_const(gridHandle)->ViewingComm().comm; 
    return EL_SUCCESS;
}

// int Grid::DiagPath( int vcRank ) const
// --------------------------------------
ElError ElGridDiagPath( ElConstGrid gridHandle, int vcRank, int* diagPath )
{ 
    *diagPath = RCG_const(gridHandle)->DiagPath(vcRank); 
    return EL_SUCCESS;
}

// int Grid::DiagPathRank( int vcRank ) const
// ------------------------------------------
ElError ElGridDiagPathRank
( ElConstGrid gridHandle, int vcRank, int* diagPathRank )
{ 
    *diagPathRank = RCG_const(gridHandle)->DiagPathRank(vcRank); 
    return EL_SUCCESS;
}

// int Grid::FirstVCRank( int diagPath ) const
// -------------------------------------------
ElError ElGridFirstVCRank
( ElConstGrid gridHandle, int vcRank, int* firstVCRank )
{ 
    *firstVCRank = RCG_const(gridHandle)->FirstVCRank(vcRank); 
    return EL_SUCCESS;
}

// static int Grid::FindFactor( int p )
// ------------------------------------
ElError ElGridFindFactor( int p, int* factor )
{ 
    *factor = Grid::FindFactor(p); 
    return EL_SUCCESS;
}

} // extern "C"
