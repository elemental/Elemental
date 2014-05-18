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

ElError ElDefaultGrid( ElConstGrid* gridHandle )
{ 
    *gridHandle = Reinterpret(&DefaultGrid()); 
    return EL_SUCCESS;
}

// Grid::Grid( MPI_Comm comm, GridOrder order )
// --------------------------------------------
ElError ElGridCreate
( MPI_Comm comm, ElGridOrderType orderC, ElGrid* gridHandle )
{
    GridOrder order = Reinterpret(orderC);
    try { *gridHandle = Reinterpret(new Grid(comm,order)); }
    CATCH
    return EL_SUCCESS;
}

// Grid::Grid( MPI_Comm comm, int height, GridOrder order )
// --------------------------------------------------------
ElError ElGridCreateSpecific
( MPI_Comm comm, int height, ElGridOrderType orderC, ElGrid* gridHandle )
{
    GridOrder order = Reinterpret(orderC);
    try { *gridHandle = Reinterpret(new Grid(comm,height,order)); }
    CATCH
    return EL_SUCCESS;
}

// Grid::~Grid()
// -------------
ElError ElGridDestroy( ElConstGrid gridHandle )
{ 
    delete Reinterpret(gridHandle); 
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
    *height = Reinterpret(gridHandle)->Height(); 
    return EL_SUCCESS;
}

// int Grid::Width() const
// -----------------------
ElError ElGridWidth( ElConstGrid gridHandle, int* width )
{ 
    *width = Reinterpret(gridHandle)->Width(); 
    return EL_SUCCESS;
}

// int Grid::Size() const
// ----------------------
ElError ElGridSize( ElConstGrid gridHandle, int* size )
{ 
    *size = Reinterpret(gridHandle)->Size(); 
    return EL_SUCCESS;
}

// GridOrder Grid::Order() const
// -----------------------------
ElError ElGridOrder( ElConstGrid gridHandle, ElGridOrderType* order )
{ 
    *order = Reinterpret(Reinterpret(gridHandle)->Order()); 
    return EL_SUCCESS;
}

// mpi::Comm Grid::ColComm() const
// -------------------------------
ElError ElGridColComm( ElConstGrid gridHandle, MPI_Comm* comm )
{ 
    *comm = Reinterpret(gridHandle)->ColComm().comm; 
    return EL_SUCCESS;
}

// mpi::Comm Grid::RowComm() const
// -------------------------------
ElError ElGridRowComm( ElConstGrid gridHandle, MPI_Comm* comm )
{ 
    *comm = Reinterpret(gridHandle)->RowComm().comm; 
    return EL_SUCCESS;
}

// mpi::Comm Grid::Comm() const
// ----------------------------
ElError ElGridComm( ElConstGrid gridHandle, MPI_Comm* comm )
{ 
    *comm = Reinterpret(gridHandle)->Comm().comm; 
    return EL_SUCCESS;
}

// int Grid::MCRank() const
// ------------------------
ElError ElGridMCRank( ElConstGrid gridHandle, int* mcRank )
{
    try { *mcRank = Reinterpret(gridHandle)->MCRank(); }
    CATCH
    return EL_SUCCESS;
}

// int Grid::MRRank() const
// ------------------------
ElError ElGridMRRank( ElConstGrid gridHandle, int* mrRank )
{
    try { *mrRank = Reinterpret(gridHandle)->MRRank(); }
    CATCH
    return EL_SUCCESS;
}

// int Grid::VCRank() const
// ------------------------
ElError ElGridVCRank( ElConstGrid gridHandle, int* vcRank )
{ 
    try { *vcRank = Reinterpret(gridHandle)->VCRank(); }
    CATCH
    return EL_SUCCESS;
}

// int Grid::VRRank() const
// ------------------------
ElError ElGridVRRank( ElConstGrid gridHandle, int* vrRank )
{ 
    try { *vrRank = Reinterpret(gridHandle)->VRRank(); }
    CATCH
    return EL_SUCCESS;
}

// int Grid::MCSize() const
// ------------------------
ElError ElGridMCSize( ElConstGrid gridHandle, int* mcSize )
{ 
    *mcSize = Reinterpret(gridHandle)->MCSize(); 
    return EL_SUCCESS;
}

// int Grid::MRSize() const
// ------------------------
ElError ElGridMRSize( ElConstGrid gridHandle, int* mrSize )
{ 
    *mrSize = Reinterpret(gridHandle)->MRSize(); 
    return EL_SUCCESS;
}

// int Grid::VCSize() const
// ------------------------
ElError ElGridVCSize( ElConstGrid gridHandle, int* vcSize )
{ 
    *vcSize = Reinterpret(gridHandle)->VCSize(); 
    return EL_SUCCESS;
}

// int Grid::VRSize() const
// ------------------------
ElError ElGridVRSize( ElConstGrid gridHandle, int* vrSize )
{ 
    *vrSize = Reinterpret(gridHandle)->VRSize(); 
    return EL_SUCCESS;
}

// mpi::Comm Grid::MCComm() const
// ------------------------------
ElError ElGridMCComm( ElConstGrid gridHandle, MPI_Comm* comm )
{ 
    *comm = Reinterpret(gridHandle)->MCComm().comm; 
    return EL_SUCCESS;
}

// mpi::Comm Grid::MRComm() const
// ------------------------------
ElError ElGridMRComm( ElConstGrid gridHandle, MPI_Comm* comm )
{ 
    *comm = Reinterpret(gridHandle)->MRComm().comm; 
    return EL_SUCCESS;
}

// mpi::Comm Grid::VCComm() const
// ------------------------------
ElError ElGridVCComm( ElConstGrid gridHandle, MPI_Comm* comm )
{ 
    *comm = Reinterpret(gridHandle)->VCComm().comm; 
    return EL_SUCCESS;
}

// mpi::Comm Grid::VRComm() const
// ------------------------------
ElError ElGridVRComm( ElConstGrid gridHandle, MPI_Comm* comm )
{ 
    *comm = Reinterpret(gridHandle)->VRComm().comm; 
    return EL_SUCCESS;
}

// mpi::Comm Grid::MDComm() const
// ------------------------------
ElError ElGridMDComm( ElConstGrid gridHandle, MPI_Comm* comm )
{ 
    *comm = Reinterpret(gridHandle)->MDComm().comm; 
    return EL_SUCCESS;
}

// mpi::Comm Grid::MDPerpComm() const
// ----------------------------------
ElError ElGridMDPerpComm( ElConstGrid gridHandle, MPI_Comm* comm )
{ 
    *comm = Reinterpret(gridHandle)->MDPerpComm().comm; 
    return EL_SUCCESS;
}

// Grid::Grid( mpi::Comm comm, mpi::Group owners, int height, GridOrder order )
// ----------------------------------------------------------------------------
ElError ElGridCreateAdvanced
( MPI_Comm comm, MPI_Group owners, int height, ElGridOrderType orderC,
  ElGrid* gridHandle )
{
    GridOrder order = Reinterpret(orderC);
    try { *gridHandle = Reinterpret(new Grid(comm,owners,height,order)); }
    CATCH
    return EL_SUCCESS;
}

// int Grid::GCD() const
// ---------------------
ElError ElGridGCD( ElConstGrid gridHandle, int* gcd )
{ 
    *gcd = Reinterpret(gridHandle)->GCD(); 
    return EL_SUCCESS;
}

// int Grid::LCM() const
// ---------------------
ElError ElGridLCM( ElConstGrid gridHandle, int* lcm )
{ 
    *lcm = Reinterpret(gridHandle)->LCM(); 
    return EL_SUCCESS;
}

// bool Grid::InGrid() const
// -------------------------
ElError ElGridInGrid( ElConstGrid gridHandle, bool* inGrid )
{ 
    *inGrid = Reinterpret(gridHandle)->InGrid(); 
    return EL_SUCCESS;
}

// bool Grid::HaveViewers() const
// ------------------------------
ElError ElGridHaveViewers( ElConstGrid gridHandle, bool* haveViewers )
{ 
    *haveViewers = Reinterpret(gridHandle)->HaveViewers(); 
    return EL_SUCCESS;
}

// int Grid::OwningRank() const
// ----------------------------
ElError ElGridOwningRank( ElConstGrid gridHandle, int* owningRank )
{ 
    *owningRank = Reinterpret(gridHandle)->OwningRank(); 
    return EL_SUCCESS;
}

// int Grid::ViewingRank() const
// -----------------------------
ElError ElGridViewingRank( ElConstGrid gridHandle, int* viewingRank )
{ 
    *viewingRank = Reinterpret(gridHandle)->ViewingRank(); 
    return EL_SUCCESS;
}

// int Grid::VCToViewingMap( int vcRank ) const
// --------------------------------------------
ElError ElGridVCToViewingMap
( ElConstGrid gridHandle, int vcRank, int* viewingRank )
{ 
    *viewingRank = Reinterpret(gridHandle)->VCToViewingMap(vcRank); 
    return EL_SUCCESS;
}

// mpi::Group Grid::OwningGroup() const
// ------------------------------------
ElError ElGridOwningGroup( ElConstGrid gridHandle, MPI_Group* group )
{ 
    *group = Reinterpret(gridHandle)->OwningGroup().group; 
    return EL_SUCCESS;
}

// mpi::Comm Grid::OwningComm() const
// ----------------------------------
ElError ElGridOwningComm( ElConstGrid gridHandle, MPI_Comm* comm )
{ 
    *comm = Reinterpret(gridHandle)->OwningComm().comm; 
    return EL_SUCCESS;
}

// mpi::Comm Grid::ViewingComm() const
// -----------------------------------
ElError ElGridViewingComm( ElConstGrid gridHandle, MPI_Comm* comm )
{ 
    *comm = Reinterpret(gridHandle)->ViewingComm().comm; 
    return EL_SUCCESS;
}

// int Grid::DiagPath( int vcRank ) const
// --------------------------------------
ElError ElGridDiagPath( ElConstGrid gridHandle, int vcRank, int* diagPath )
{ 
    *diagPath = Reinterpret(gridHandle)->DiagPath(vcRank); 
    return EL_SUCCESS;
}

// int Grid::DiagPathRank( int vcRank ) const
// ------------------------------------------
ElError ElGridDiagPathRank
( ElConstGrid gridHandle, int vcRank, int* diagPathRank )
{ 
    *diagPathRank = Reinterpret(gridHandle)->DiagPathRank(vcRank); 
    return EL_SUCCESS;
}

// int Grid::FirstVCRank( int diagPath ) const
// -------------------------------------------
ElError ElGridFirstVCRank
( ElConstGrid gridHandle, int vcRank, int* firstVCRank )
{ 
    *firstVCRank = Reinterpret(gridHandle)->FirstVCRank(vcRank); 
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
