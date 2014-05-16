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

ElConstGrid ElDefaultGrid();

ElGrid ElGridCreate( MPI_Comm comm, ElGridOrderType order );
ElGrid ElGridCreateSpecific
( MPI_Comm comm, int height, ElGridOrderType order );
void ElGridDestroy( ElConstGrid grid );

int ElGridRow( ElConstGrid grid );
int ElGridCol( ElConstGrid grid );
int ElGridRank( ElConstGrid grid );
int ElGridHeight( ElConstGrid grid );
int ElGridWidth( ElConstGrid grid );
int ElGridSize( ElConstGrid grid );

ElGridOrderType ElGridOrder( ElConstGrid grid );
MPI_Comm ElGridColComm( ElConstGrid grid );
MPI_Comm ElGridRowComm( ElConstGrid grid );
MPI_Comm ElGridComm( ElConstGrid grid );

int ElGridMCRank( ElConstGrid grid );
int ElGridMRRank( ElConstGrid grid );
int ElGridVCRank( ElConstGrid grid );
int ElGridVRRank( ElConstGrid grid );
int ElGridMCSize( ElConstGrid grid );
int ElGridMRSize( ElConstGrid grid );
int ElGridVCSize( ElConstGrid grid );
int ElGridVRSize( ElConstGrid grid );

MPI_Comm ElGridMCComm( ElConstGrid grid );
MPI_Comm ElGridMRComm( ElConstGrid grid );
MPI_Comm ElGridVCComm( ElConstGrid grid );
MPI_Comm ElGridVRComm( ElConstGrid grid );
MPI_Comm ElGridMDComm( ElConstGrid grid );
MPI_Comm ElGridMDPerpComm( ElConstGrid grid );

ElGrid ElGridCreateAdvanced
( MPI_Comm comm, MPI_Group owners, int height, ElGridOrderType order );

int  ElGridGCD( ElConstGrid grid );
int  ElGridLCM( ElConstGrid grid );
bool ElGridInGrid( ElConstGrid grid );
bool ElGridHaveViewers( ElConstGrid grid );

int ElGridOwningRank( ElConstGrid grid );
int ElGridViewingRank( ElConstGrid grid );
int ElGridVCToViewingMap( ElConstGrid grid, int VCRank );

MPI_Group ElGridOwningGroup( ElConstGrid grid );
MPI_Comm  ElGridOwningComm( ElConstGrid grid );
MPI_Comm  ElGridViewingComm( ElConstGrid grid );

int ElGridDiagPath( ElConstGrid grid, int VCRank );
int ElGridDiagPathRank( ElConstGrid grid, int VCRank );
int ElGridFirstVCRank( ElConstGrid grid, int diagPath );

int ElGridFindFactor( int p );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_GRID_C_H */
