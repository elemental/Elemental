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

struct ElGrid; typedef struct ElGrid ElGrid;

const ElGrid* ElDefaultGrid();

ElGrid* ElGridCreate( MPI_Comm comm, ElGridOrderType order );
ElGrid* ElGridCreateSpecific
( MPI_Comm comm, int height, ElGridOrderType order );
void ElGridDestroy( const ElGrid* grid );

int ElGridRow( const ElGrid* grid );
int ElGridCol( const ElGrid* grid );
int ElGridRank( const ElGrid* grid );
int ElGridHeight( const ElGrid* grid );
int ElGridWidth( const ElGrid* grid );
int ElGridSize( const ElGrid* grid );

ElGridOrderType ElGridOrder( const ElGrid* grid );
MPI_Comm ElGridColComm( const ElGrid* grid );
MPI_Comm ElGridRowComm( const ElGrid* grid );
MPI_Comm ElGridComm( const ElGrid* grid );

int ElGridMCRank( const ElGrid* grid );
int ElGridMRRank( const ElGrid* grid );
int ElGridVCRank( const ElGrid* grid );
int ElGridVRRank( const ElGrid* grid );
int ElGridMCSize( const ElGrid* grid );
int ElGridMRSize( const ElGrid* grid );
int ElGridVCSize( const ElGrid* grid );
int ElGridVRSize( const ElGrid* grid );

MPI_Comm ElGridMCComm( const ElGrid* grid );
MPI_Comm ElGridMRComm( const ElGrid* grid );
MPI_Comm ElGridVCComm( const ElGrid* grid );
MPI_Comm ElGridVRComm( const ElGrid* grid );
MPI_Comm ElGridMDComm( const ElGrid* grid );
MPI_Comm ElGridMDPerpComm( const ElGrid* grid );

ElGrid* ElGridCreateAdvanced
( MPI_Comm comm, MPI_Group owners, int height, ElGridOrderType order );

int  ElGridGCD( const ElGrid* grid );
int  ElGridLCM( const ElGrid* grid );
bool ElGridInGrid( const ElGrid* grid );
bool ElGridHaveViewers( const ElGrid* grid );

int ElGridOwningRank( const ElGrid* grid );
int ElGridViewingRank( const ElGrid* grid );
int ElGridVCToViewingMap( const ElGrid* grid, int VCRank );

MPI_Group ElGridOwningGroup( const ElGrid* grid );
MPI_Comm  ElGridOwningComm( const ElGrid* grid );
MPI_Comm  ElGridViewingComm( const ElGrid* grid );

int ElGridDiagPath( const ElGrid* grid, int VCRank );
int ElGridDiagPathRank( const ElGrid* grid, int VCRank );
int ElGridFirstVCRank( const ElGrid* grid, int diagPath );

int ElGridFindFactor( int p );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_GRID_C_H */
