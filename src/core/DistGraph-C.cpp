/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
#include "El.h"
using namespace El;

extern "C" {

ElError ElDistGraphCreate( ElDistGraph* graph, MPI_Comm comm )
{ EL_TRY( *graph = CReflect( new DistGraph(mpi::Comm(comm)) ) ) }

ElError ElDistGraphDestroy( ElConstDistGraph graph )
{ EL_TRY( delete CReflect(graph) ) }

ElError ElDistGraphEmpty( ElDistGraph graph )
{ EL_TRY( CReflect(graph)->Empty() ) }

ElError ElDistGraphResize
( ElDistGraph graph, ElInt numSources, ElInt numTargets )
{ EL_TRY( CReflect(graph)->Resize( numSources, numTargets ) ) }

ElError ElDistGraphSetComm( ElDistGraph graph, MPI_Comm comm )
{ EL_TRY( CReflect(graph)->SetComm(comm) ) }

ElError ElDistGraphReserve
( ElDistGraph graph, ElInt numLocalEdges, ElInt numRemoteEdges )
{ EL_TRY( CReflect(graph)->Reserve(numLocalEdges,numRemoteEdges) ) }

ElError 
ElDistGraphConnect( ElDistGraph graph, ElInt row, ElInt col )
{ EL_TRY( CReflect(graph)->Connect( row, col ) ) }

ElError ElDistGraphConnectLocal( ElDistGraph graph, ElInt localRow, ElInt col )
{ EL_TRY( CReflect(graph)->ConnectLocal( localRow, col ) ) }

ElError 
ElDistGraphDisconnect( ElDistGraph graph, ElInt row, ElInt col )
{ EL_TRY( CReflect(graph)->Disconnect( row, col ) ) }

ElError ElDistGraphDisconnectLocal
( ElDistGraph graph, ElInt localRow, ElInt col )
{ EL_TRY( CReflect(graph)->DisconnectLocal( localRow, col ) ) }

ElError ElDistGraphQueueConnection
( ElDistGraph graph, ElInt row, ElInt col, bool passive )
{ EL_TRY( CReflect(graph)->QueueConnection( row, col, passive ) ) }

ElError ElDistGraphQueueLocalConnection
( ElDistGraph graph, ElInt localRow, ElInt col )
{ EL_TRY( CReflect(graph)->QueueLocalConnection( localRow, col ) ) }

ElError ElDistGraphQueueDisconnection
( ElDistGraph graph, ElInt row, ElInt col, bool passive )
{ EL_TRY( CReflect(graph)->QueueDisconnection( row, col, passive ) ) }

ElError ElDistGraphQueueLocalDisconnection
( ElDistGraph graph, ElInt localRow, ElInt col )
{ EL_TRY( CReflect(graph)->QueueLocalDisconnection( localRow, col ) ) }

ElError ElDistGraphProcessQueues( ElDistGraph graph )
{ EL_TRY( CReflect(graph)->ProcessQueues() ) }

ElError ElDistGraphProcessLocalQueues( ElDistGraph graph )
{ EL_TRY( CReflect(graph)->ProcessLocalQueues() ) }

ElError ElDistGraphNumSources( ElConstDistGraph graph, ElInt* numSources )
{ EL_TRY( *numSources = CReflect(graph)->NumSources() ) } 

ElError ElDistGraphNumTargets( ElConstDistGraph graph, ElInt* numTargets )
{ EL_TRY( *numTargets = CReflect(graph)->NumTargets() ) } 

ElError ElDistGraphFirstLocalSource
( ElConstDistGraph graph, ElInt* firstLocalSource )
{ EL_TRY( *firstLocalSource = CReflect(graph)->FirstLocalSource() ) }

ElError ElDistGraphNumLocalSources
( ElConstDistGraph graph, ElInt* numLocalSources )
{ EL_TRY( *numLocalSources = CReflect(graph)->NumLocalSources() ) }

ElError ElDistGraphNumLocalEdges( ElConstDistGraph graph, ElInt* numLocalEdges )
{ EL_TRY( *numLocalEdges = CReflect(graph)->NumLocalEdges() ) } 

ElError ElDistGraphCapacity( ElConstDistGraph graph, ElInt* capacity )
{ EL_TRY( *capacity = CReflect(graph)->Capacity() ) } 

ElError ElDistGraphLocallyConsistent( ElConstDistGraph graph, bool* consistent )
{ EL_TRY( *consistent = CReflect(graph)->LocallyConsistent() ) }

ElError ElDistGraphComm( ElConstDistGraph graph, MPI_Comm* comm )
{ EL_TRY( *comm = CReflect(graph)->Comm().comm ) }

ElError ElDistGraphBlocksize( ElConstDistGraph graph, ElInt* blocksize )
{ EL_TRY( *blocksize = CReflect(graph)->Blocksize() ) }

ElError ElDistGraphSource
( ElConstDistGraph graph, ElInt localEdge, ElInt* source )
{ EL_TRY( *source = CReflect(graph)->Source(localEdge) ) }

ElError ElDistGraphTarget
( ElConstDistGraph graph, ElInt localEdge, ElInt* target )
{ EL_TRY( *target = CReflect(graph)->Target(localEdge) ) }

ElError ElDistGraphSourceOffset
( ElConstDistGraph graph, ElInt localSource, ElInt* localSourceOffset )
{ EL_TRY( *localSourceOffset = CReflect(graph)->SourceOffset(localSource) ) }

ElError ElDistGraphOffset
( ElConstDistGraph graph, ElInt localSource, ElInt target, ElInt* localOffset )
{ EL_TRY( *localOffset = CReflect(graph)->Offset(localSource,target) ) }

ElError ElDistGraphNumConnections
( ElConstDistGraph graph, ElInt localSource, ElInt* numConnections )
{ EL_TRY( *numConnections = CReflect(graph)->NumConnections(localSource) ) }

ElError ElDistGraphImbalance( ElConstDistGraph graph, double* imbalance )
{ EL_TRY( *imbalance = CReflect(graph)->Imbalance() ) }

ElError ElDistGraphSourceBuffer( ElDistGraph graph, ElInt** sourceBuffer )
{ EL_TRY( *sourceBuffer = CReflect(graph)->SourceBuffer() ) }

ElError ElDistGraphLockedSourceBuffer
( ElConstDistGraph graph, const ElInt** sourceBuffer )
{ EL_TRY( *sourceBuffer = CReflect(graph)->LockedSourceBuffer() ) }

ElError ElDistGraphTargetBuffer( ElDistGraph graph, ElInt** targetBuffer )
{ EL_TRY( *targetBuffer = CReflect(graph)->TargetBuffer() ) }

ElError ElDistGraphLockedTargetBuffer
( ElConstDistGraph graph, const ElInt** targetBuffer )
{ EL_TRY( *targetBuffer = CReflect(graph)->LockedTargetBuffer() ) }

} // extern "C"
