/*
   Copyright (c) 2009-2014, Jack Poulson
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

ElError ElGraphMakeConsistent( ElGraph graph )
{ EL_TRY( CReflect(graph)->MakeConsistent() ) }

ElError ElGraphReserve( ElGraph graph, ElInt numEdges )
{ EL_TRY( CReflect(graph)->Reserve(numEdges) ) }

ElError ElGraphInsert( ElGraph graph, ElInt row, ElInt col )
{ EL_TRY( CReflect(graph)->Insert( row, col ) ) }

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

ElError ElGraphEdgeOffset( ElConstGraph graph, ElInt source, ElInt* edgeOffset )
{ EL_TRY( *edgeOffset = CReflect(graph)->EdgeOffset(source) ) }

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
