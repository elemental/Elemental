/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_DISTGRAPH_C_H
#define EL_DISTGRAPH_C_H

#ifdef __cplusplus
extern "C" {
#endif

/* An anonymous struct meant as a placeholder for DistGraph
   ---------------------------------------------------- */
typedef struct ElDistGraph_Dummy* ElDistGraph;
typedef const struct ElDistGraph_Dummy* ElConstDistGraph;

/* Constructors and destructors
   ============================ */

/* DistGraph::DistGraph( mpi::Comm comm )
   -------------------------------------- */
EL_EXPORT ElError ElDistGraphCreate( ElDistGraph* graph, MPI_Comm comm );

/* DistGraph::~DistGraph()
   ----------------------- */
EL_EXPORT ElError ElDistGraphDestroy( ElConstDistGraph graph );

/* Assignment and reconfiguration
   ============================== */

/* void DistGraph::Empty()
   ----------------------- */
EL_EXPORT ElError ElDistGraphEmpty( ElDistGraph graph );

/* void DistGraph::Resize( Int numSources, Int numTargets )
   -------------------------------------------------------- */
EL_EXPORT ElError ElDistGraphResize
( ElDistGraph graph, ElInt numSources, ElInt numTargets );

/* void DistGraph::SetComm( mpi::Comm comm )
   ----------------------------------------- */
EL_EXPORT ElError ElDistGraphSetComm( ElDistGraph graph, MPI_Comm comm );

/* void DistGraph::Reserve( Int numLocalEdges, Int numRemoteEdges )
   ---------------------------------------------------------------- */
EL_EXPORT ElError ElDistGraphReserve
( ElDistGraph graph, ElInt numLocalEdges, ElInt numRemoteEdges );

/* void DistGraph::Connect( Int row, Int col )
   ------------------------------------------- */
EL_EXPORT ElError 
ElDistGraphConnect( ElDistGraph graph, ElInt row, ElInt col );

/* void DistGraph::ConnectLocal( Int localRow, Int col )
   ----------------------------------------------------- */
EL_EXPORT ElError ElDistGraphConnectLocal
( ElDistGraph graph, ElInt localRow, ElInt col );

/* void DistGraph::Disconnect( Int row, Int col )
   ---------------------------------------------- */
EL_EXPORT ElError ElDistGraphDisconnect
( ElDistGraph graph, ElInt row, ElInt col );

/* void DistGraph::DisconnectLocal( Int localRow, Int col )
   -------------------------------------------------------- */
EL_EXPORT ElError ElDistGraphDisconnectLocal
( ElDistGraph graph, ElInt localRow, ElInt col );

/* void DistGraph::QueueConnection( Int row, Int col, bool passive )
   ----------------------------------------------------------------- */
EL_EXPORT ElError ElDistGraphQueueConnection
( ElDistGraph graph, ElInt row, ElInt col, bool passive );

/* void DistGraph::QueueLocalConnection( Int localRow, Int col )
   ------------------------------------------------------------- */
EL_EXPORT ElError ElDistGraphQueueLocalConnection
( ElDistGraph graph, ElInt localRow, ElInt col );

/* void DistGraph::QueueDisconnection( Int row, Int col, bool passive )
   -------------------------------------------------------------------- */
EL_EXPORT ElError ElDistGraphQueueDisconnection
( ElDistGraph graph, ElInt row, ElInt col, bool passive );

/* void DistGraph::QueueLocalDisconnection( Int localRow, Int col )
   ---------------------------------------------------------------- */
EL_EXPORT ElError ElDistGraphQueueLocalDisconnection
( ElDistGraph graph, ElInt localRow, ElInt col );

/* void DistGraph::ProcessQueues()
   ------------------------------- */ 
EL_EXPORT ElError ElDistGraphProcessQueues( ElDistGraph graph );

/* void DistGraph::ProcessLocalQueues()
   ------------------------------------ */
EL_EXPORT ElError ElDistGraphProcessLocalQueues( ElDistGraph graph );

/* Queries
   ======= */

/* Int DistGraph::NumSources() const
   --------------------------------- */
EL_EXPORT ElError ElDistGraphNumSources
( ElConstDistGraph graph, ElInt* numSources );

/* Int DistGraph::NumTargets() const
   --------------------------------- */
EL_EXPORT ElError ElDistGraphNumTargets
( ElConstDistGraph graph, ElInt* numTargets );

/* Int DistGraph::FirstLocalSource() const
   --------------------------------------- */
EL_EXPORT ElError ElDistGraphFirstLocalSource
( ElConstDistGraph graph, ElInt* firstLocalSource );

/* Int DistGraph::NumLocalSources() const
   -------------------------------------- */
EL_EXPORT ElError ElDistGraphNumLocalSources
( ElConstDistGraph graph, ElInt* numLocalSources );

/* Int DistGraph::NumLocalEdges() const
   ------------------------------------ */
EL_EXPORT ElError ElDistGraphNumLocalEdges
( ElConstDistGraph graph, ElInt* numLocalEdges );

/* Int DistGraph::Capacity() const
   ------------------------------- */
EL_EXPORT ElError ElDistGraphCapacity
( ElConstDistGraph graph, ElInt* capacity );

/* bool DistGraph::Consistent() const
   ---------------------------------- */
EL_EXPORT ElError ElDistGraphConsistent
( ElConstDistGraph graph, bool* consistent );

/* mpi::Comm DistGraph::Comm() const
   --------------------------------- */
EL_EXPORT ElError ElDistGraphComm( ElConstDistGraph graph, MPI_Comm* comm );

/* Int DistGraph::Blocksize() const
   -------------------------------- */
EL_EXPORT ElError ElDistGraphBlocksize
( ElConstDistGraph graph, ElInt* blocksize );

/* Int DistGraph::Source( Int localEdge ) const
   -------------------------------------------- */
EL_EXPORT ElError ElDistGraphSource
( ElConstDistGraph graph, ElInt localEdge, ElInt* source );

/* Int DistGraph::Target( Int localEdge ) const
   -------------------------------------------- */
EL_EXPORT ElError ElDistGraphTarget
( ElConstDistGraph graph, ElInt localEdge, ElInt* target );

/* Int DistGraph::SourceOffset( Int localSource ) const
   ---------------------------------------------------- */
EL_EXPORT ElError ElDistGraphSourceOffset
( ElConstDistGraph graph, ElInt localSource, ElInt* localSourceOffset );

/* Int DistGraph::Offset( Int localSource, Int target ) const
   ---------------------------------------------------------- */
EL_EXPORT ElError ElDistGraphOffset
( ElConstDistGraph graph, ElInt localSource, ElInt target, ElInt* localOffset );

/* Int DistGraph::NumConnections( Int localSource ) const
   ------------------------------------------------------ */
EL_EXPORT ElError ElDistGraphNumConnections
( ElConstDistGraph graph, ElInt localSource, ElInt* numConnections );

/* double DistGraph::Imbalance() const
   ----------------------------------- */
EL_EXPORT ElError ElDistGraphImbalance( ElConstDistGraph A, double* imbalance );

/* Int* DistGraph::SourceBuffer()
   ------------------------------ */
EL_EXPORT ElError ElDistGraphSourceBuffer
( ElDistGraph graph, ElInt** sourceBuffer );

/* Int* DistGraph::LockedSourceBuffer() const
   ------------------------------------------ */
EL_EXPORT ElError ElDistGraphLockedSourceBuffer
( ElConstDistGraph graph, const ElInt** sourceBuffer );

/* Int* DistGraph::TargetBuffer()
   ------------------------------ */
EL_EXPORT ElError ElDistGraphTargetBuffer
( ElDistGraph graph, ElInt** targetBuffer );

/* Int* DistGraph::LockedTargetBuffer() const
   ------------------------------------------ */
EL_EXPORT ElError ElDistGraphLockedTargetBuffer
( ElConstDistGraph graph, const ElInt** targetBuffer );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_DISTGRAPH_C_H */
