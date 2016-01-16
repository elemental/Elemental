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

ElError ElGraphCreate( ElGraph* graph )
{ EL_TRY( *graph = CReflect( new Graph() ) ) }

ElError ElGraphDestroy( ElConstGraph graph )
{ EL_TRY( delete CReflect(graph) ) }

ElError ElGraphEmpty( ElGraph graph )
{ EL_TRY( CReflect(graph)->Empty() ) }

ElError ElGraphResize( ElGraph graph, ElInt numSources, ElInt numTargets )
{ EL_TRY( CReflect(graph)->Resize( numSources, numTargets ) ) }

ElError ElGraphReserve( ElGraph graph, ElInt numEdges )
{ EL_TRY( CReflect(graph)->Reserve(numEdges) ) }

ElError ElGraphConnect( ElGraph graph, ElInt source, ElInt target )
{ EL_TRY( CReflect(graph)->Connect( source, target ) ) }

ElError ElGraphDisconnect( ElGraph graph, ElInt source, ElInt target )
{ EL_TRY( CReflect(graph)->Disconnect( source, target ) ) }

ElError ElGraphQueueConnection( ElGraph graph, ElInt source, ElInt target )
{ EL_TRY( CReflect(graph)->QueueConnection( source, target ) ) }

ElError ElGraphQueueDisconnection( ElGraph graph, ElInt source, ElInt target )
{ EL_TRY( CReflect(graph)->QueueDisconnection( source, target ) ) }

ElError ElGraphProcessQueues( ElGraph graph )
{ EL_TRY( CReflect(graph)->ProcessQueues() ) }

ElError ElGraphNumSources( ElConstGraph graph, ElInt* numSources )
{ EL_TRY( *numSources = CReflect(graph)->NumSources() ) } 

ElError ElGraphNumTargets( ElConstGraph graph, ElInt* numTargets )
{ EL_TRY( *numTargets = CReflect(graph)->NumTargets() ) } 

ElError ElGraphNumEdges( ElConstGraph graph, ElInt* numEdges )
{ EL_TRY( *numEdges = CReflect(graph)->NumEdges() ) } 

ElError ElGraphCapacity( ElConstGraph graph, ElInt* capacity )
{ EL_TRY( *capacity = CReflect(graph)->Capacity() ) } 

ElError ElGraphConsistent( ElConstGraph graph, bool* consistent )
{ EL_TRY( *consistent = CReflect(graph)->Consistent() ) }

ElError ElGraphSource( ElConstGraph graph, ElInt edge, ElInt* source )
{ EL_TRY( *source = CReflect(graph)->Source(edge) ) }

ElError ElGraphTarget( ElConstGraph graph, ElInt edge, ElInt* target )
{ EL_TRY( *target = CReflect(graph)->Target(edge) ) }

ElError ElGraphSourceOffset
( ElConstGraph graph, ElInt source, ElInt* sourceOffset )
{ EL_TRY( *sourceOffset = CReflect(graph)->SourceOffset(source) ) }

ElError ElGraphOffset
( ElConstGraph graph, ElInt source, ElInt target, ElInt* offset )
{ EL_TRY( *offset = CReflect(graph)->Offset(source,target) ) }

ElError ElGraphNumConnections
( ElConstGraph graph, ElInt source, ElInt* numConnections )
{ EL_TRY( *numConnections = CReflect(graph)->NumConnections(source) ) }

ElError ElGraphSourceBuffer( ElGraph graph, ElInt** sourceBuffer )
{ EL_TRY( *sourceBuffer = CReflect(graph)->SourceBuffer() ) }

ElError ElGraphLockedSourceBuffer
( ElConstGraph graph, const ElInt** sourceBuffer )
{ EL_TRY( *sourceBuffer = CReflect(graph)->LockedSourceBuffer() ) }

ElError ElGraphTargetBuffer( ElGraph graph, ElInt** targetBuffer )
{ EL_TRY( *targetBuffer = CReflect(graph)->TargetBuffer() ) }

ElError ElGraphLockedTargetBuffer
( ElConstGraph graph, const ElInt** targetBuffer )
{ EL_TRY( *targetBuffer = CReflect(graph)->LockedTargetBuffer() ) }

} // extern "C"
