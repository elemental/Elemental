/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_GRAPH_C_H
#define EL_GRAPH_C_H

#ifdef __cplusplus
extern "C" {
#endif

/* An anonymous struct meant as a placeholder for Graph
   ---------------------------------------------------- */
typedef struct ElGraph_Dummy* ElGraph;
typedef const struct ElGraph_Dummy* ElConstGraph;

/* Constructors and destructors
   ============================ */

/* Graph::Graph()
   -------------- */
EL_EXPORT ElError ElGraphCreate( ElGraph* graph );

/* Graph::~Graph()
   --------------- */
EL_EXPORT ElError ElGraphDestroy( ElConstGraph graph );

/* Assignment and reconfiguration
   ============================== */

/* void Graph::Empty()
   ------------------- */
EL_EXPORT ElError ElGraphEmpty( ElGraph graph );

/* void Graph::Resize( Int numSources, Int numTargets )
   ---------------------------------------------------- */
EL_EXPORT ElError ElGraphResize
( ElGraph graph, ElInt numSources, ElInt numTargets );

/* void Graph::Reserve( Int numEdges )
   ----------------------------------- */
EL_EXPORT ElError ElGraphReserve( ElGraph graph, ElInt numEdges );

/* void Graph::Connect( Int source, Int target )
   --------------------------------------------- */
EL_EXPORT ElError ElGraphConnect( ElGraph graph, ElInt source, ElInt target );

/* void Graph::Disconnect( Int source, Int target )
   ------------------------------------------------ */
EL_EXPORT ElError ElGraphDisconnect
( ElGraph graph, ElInt source, ElInt target );

/* void Graph::QueueConnection( Int source, Int target )
   ----------------------------------------------------- */
EL_EXPORT ElError ElGraphQueueConnection
( ElGraph graph, ElInt source, ElInt target );

/* void Graph::QueueDisconnection( Int source, Int target )
   -------------------------------------------------------- */
EL_EXPORT ElError ElGraphQueueDisconnection
( ElGraph graph, ElInt source, ElInt target );

/* void Graph::ProcessQueues()
   --------------------------- */ 
EL_EXPORT ElError ElGraphProcessQueues( ElGraph graph );

/* Queries
   ======= */

/* Int Graph::NumSources() const
   ----------------------------- */
EL_EXPORT ElError ElGraphNumSources( ElConstGraph graph, ElInt* numSources );

/* Int Graph::NumTargets() const
   ----------------------------- */
EL_EXPORT ElError ElGraphNumTargets( ElConstGraph graph, ElInt* numTargets );

/* Int Graph::NumEdges() const
   --------------------------- */
EL_EXPORT ElError ElGraphNumEdges( ElConstGraph graph, ElInt* numEdges );

/* Int Graph::Capacity() const
   --------------------------- */
EL_EXPORT ElError ElGraphCapacity( ElConstGraph graph, ElInt* capacity );

/* bool Graph::Consistent() const
   ------------------------------ */
EL_EXPORT ElError ElGraphConsistent( ElConstGraph graph, bool* consistent );

/* Int Graph::Source( Int edge ) const
   ----------------------------------- */
EL_EXPORT ElError ElGraphSource
( ElConstGraph graph, ElInt edge, ElInt* source );

/* Int Graph::Target( Int edge ) const
   ------------------------------------ */
EL_EXPORT ElError ElGraphTarget
( ElConstGraph graph, ElInt edge, ElInt* target );

/* Int Graph::SourceOffset( Int source ) const
   ------------------------------------------- */
EL_EXPORT ElError ElGraphSourceOffset
( ElConstGraph graph, ElInt source, ElInt* sourceOffset );

/* Int Graph::Offset( Int source, Int target ) const
   ------------------------------------------------- */
EL_EXPORT ElError ElGraphOffset
( ElConstGraph graph, ElInt source, ElInt target, ElInt* offset );

/* Int Graph::NumConnections( Int source ) const
   --------------------------------------------- */
EL_EXPORT ElError ElGraphNumConnections
( ElConstGraph graph, ElInt source, ElInt* numConnections );

/* Int* Graph::SourceBuffer()
   -------------------------- */
EL_EXPORT ElError ElGraphSourceBuffer( ElGraph graph, ElInt** sourceBuffer );

/* Int* Graph::LockedSourceBuffer() const
   -------------------------------------- */
EL_EXPORT ElError ElGraphLockedSourceBuffer
( ElConstGraph graph, const ElInt** sourceBuffer );

/* Int* Graph::TargetBuffer()
   -------------------------- */
EL_EXPORT ElError ElGraphTargetBuffer( ElGraph graph, ElInt** targetBuffer );

/* Int* Graph::LockedTargetBuffer() const
   -------------------------------------- */
EL_EXPORT ElError ElGraphLockedTargetBuffer
( ElConstGraph graph, const ElInt** targetBuffer );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_GRAPH_C_H */
