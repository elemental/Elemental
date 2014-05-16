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

ElError ElDefaultGrid( ElConstGrid* grid );

ElError ElGridCreate( MPI_Comm comm, ElGridOrderType order, ElGrid* grid );
ElError ElGridCreateSpecific
( MPI_Comm comm, int height, ElGridOrderType order, ElGrid* grid );
ElError ElGridDestroy( ElConstGrid grid );

ElError ElGridRow( ElConstGrid grid, int* row );
ElError ElGridCol( ElConstGrid grid, int* col );
ElError ElGridRank( ElConstGrid grid, int* rank );
ElError ElGridHeight( ElConstGrid grid, int* height );
ElError ElGridWidth( ElConstGrid grid, int* width );
ElError ElGridSize( ElConstGrid grid, int* size );

ElError ElGridOrder( ElConstGrid grid, ElGridOrderType* order );
ElError ElGridColComm( ElConstGrid grid, MPI_Comm* comm );
ElError ElGridRowComm( ElConstGrid grid, MPI_Comm* comm );
ElError ElGridComm( ElConstGrid grid, MPI_Comm* comm );

ElError ElGridMCRank( ElConstGrid grid, int* mcRank );
ElError ElGridMRRank( ElConstGrid grid, int* mrRank );
ElError ElGridVCRank( ElConstGrid grid, int* vcRank );
ElError ElGridVRRank( ElConstGrid grid, int* vrRank );
ElError ElGridMCSize( ElConstGrid grid, int* mcSize );
ElError ElGridMRSize( ElConstGrid grid, int* mrSize );
ElError ElGridVCSize( ElConstGrid grid, int* vcSize );
ElError ElGridVRSize( ElConstGrid grid, int* vrSize );

ElError ElGridMCComm( ElConstGrid grid, MPI_Comm* comm );
ElError ElGridMRComm( ElConstGrid grid, MPI_Comm* comm );
ElError ElGridVCComm( ElConstGrid grid, MPI_Comm* comm );
ElError ElGridVRComm( ElConstGrid grid, MPI_Comm* comm );
ElError ElGridMDComm( ElConstGrid grid, MPI_Comm* comm );
ElError ElGridMDPerpComm( ElConstGrid grid, MPI_Comm* comm );

ElError ElGridCreateAdvanced
( MPI_Comm comm, MPI_Group owners, int height, ElGridOrderType order, 
  ElGrid* grid );

ElError ElGridGCD( ElConstGrid grid, int* gcd );
ElError ElGridLCM( ElConstGrid grid, int* lcm );
ElError ElGridInGrid( ElConstGrid grid, bool* inGrid );
ElError ElGridHaveViewers( ElConstGrid grid, bool* haveViewers );

ElError ElGridOwningRank( ElConstGrid grid, int* owningRank );
ElError ElGridViewingRank( ElConstGrid grid, int* viewingRank );
ElError ElGridVCToViewingMap( ElConstGrid grid, int VCRank, int* viewingRank );

ElError ElGridOwningGroup( ElConstGrid grid, MPI_Group* group );
ElError ElGridOwningComm( ElConstGrid grid, MPI_Comm* comm );
ElError ElGridViewingComm( ElConstGrid grid, MPI_Comm* comm );

ElError ElGridDiagPath( ElConstGrid grid, int VCRank, int* diagPath );
ElError ElGridDiagPathRank( ElConstGrid grid, int VCRank, int* diagPathRank );
ElError ElGridFirstVCRank( ElConstGrid grid, int diagPath, int* firstVCRank );

ElError ElGridFindFactor( int p, int* factor );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_GRID_C_H */
