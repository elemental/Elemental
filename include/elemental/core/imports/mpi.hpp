/*
   Copyright (c) 2009-2013, Jack Poulson
                      2013, Jeff Hammond
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef CORE_MPI_HPP
#define CORE_MPI_HPP

namespace elem {
namespace mpi {

#if defined(HAVE_MPI3_NONBLOCKING_COLLECTIVES) || \
    defined(HAVE_MPIX_NONBLOCKING_COLLECTIVES)
#define HAVE_NONBLOCKING 1
#else
#define HAVE_NONBLOCKING 0
#endif

#ifdef HAVE_NONBLOCKING_COLLECTIVES
#ifdef HAVE_MPI3_NONBLOCKING_COLLECTIVES
#define NONBLOCKING_COLL(name) MPI_ ## name
#else
#define NONBLOCKING_COLL(name) MPIX_ ## name
#endif
#endif

// Datatype definitions
typedef MPI_Comm Comm;
typedef MPI_Datatype Datatype;
typedef MPI_Errhandler ErrorHandler;
typedef MPI_Group Group;
typedef MPI_Op Op;
typedef MPI_Request Request;
typedef MPI_Status Status;
typedef MPI_User_function UserFunction;

// Standard constants
const int ANY_SOURCE = MPI_ANY_SOURCE;
const int ANY_TAG = MPI_ANY_TAG;
const int THREAD_SINGLE = MPI_THREAD_SINGLE;
const int THREAD_FUNNELED = MPI_THREAD_FUNNELED;
const int THREAD_SERIALIZED = MPI_THREAD_SERIALIZED;
const int THREAD_MULTIPLE = MPI_THREAD_MULTIPLE;
const int UNDEFINED = MPI_UNDEFINED;
const Comm COMM_SELF = MPI_COMM_SELF;
const Comm COMM_WORLD = MPI_COMM_WORLD;
const ErrorHandler ERRORS_RETURN = MPI_ERRORS_RETURN;
const ErrorHandler ERRORS_ARE_FATAL = MPI_ERRORS_ARE_FATAL;
const Group GROUP_EMPTY = MPI_GROUP_EMPTY;
const Request REQUEST_NULL = MPI_REQUEST_NULL;
const Op MAX = MPI_MAX;
const Op MIN = MPI_MIN;
const Op PROD = MPI_PROD;
const Op SUM = MPI_SUM;
const Op LOGICAL_AND = MPI_LAND;
const Op LOGICAL_OR = MPI_LOR;
const Op LOGICAL_XOR = MPI_LXOR;
const Op BINARY_AND = MPI_BAND;
const Op BINARY_OR = MPI_BOR;
const Op BINARY_XOR = MPI_BXOR;

// Added constant(s)
const int MIN_COLL_MSG = 1; // minimum message size for collectives

//----------------------------------------------------------------------------//
// Routines                                                                   //
//----------------------------------------------------------------------------//

// Environment routines
void Initialize( int& argc, char**& argv );
int InitializeThread( int& argc, char**& argv, int required );
void Finalize();
bool Initialized();
bool Finalized();
int QueryThread();
double Time();
void OpCreate( UserFunction* func, bool commutes, Op& op );
void OpFree( Op& op );

// Communicator manipulation
int WorldRank();
int CommRank( Comm comm );
int CommSize( Comm comm );
void CommCreate( Comm parentComm, Group subsetGroup, Comm& subsetComm );
void CommDup( Comm original, Comm& duplicate );
void CommSplit( Comm comm, int color, int key, Comm& newComm );
void CommFree( Comm& comm );
bool CongruentComms( Comm comm1, Comm comm2 );
void ErrorHandlerSet( Comm comm, ErrorHandler errorHandler );

// Cartesian communicator routines
void CartCreate
( Comm comm, int numDims, const int* dimensions, const int* periods, 
  bool reorder, Comm& cartComm );
void CartSub
( Comm comm, const int* remainingDims, Comm& subComm );

// Group manipulation
int GroupRank( Group group );
int GroupSize( Group group );
void CommGroup( Comm comm, Group& group );
void GroupIncl( Group group, int n, const int* ranks, Group& subGroup );
void GroupDifference( Group parent, Group subset, Group& complement );
void GroupFree( Group& group );
void GroupTranslateRanks
( Group origGroup, int size, const int* origRanks, 
  Group newGroup,                  int* newRanks );

// Utilities
void Barrier( Comm comm );
void Wait( Request& request );
void Wait( Request& request, Status& status );
void WaitAll( int numRequests, Request* requests, Status* statuses );
bool Test( Request& request );
bool IProbe( int source, int tag, Comm comm, Status& status );

template<typename T>
int GetCount( Status& status );

// Point-to-point communication
template<typename R>
void Send( const R* buf, int count, int to, int tag, Comm comm );
template<typename R>
void Send( const Complex<R>* buf, int count, int to, int tag, Comm comm );

template<typename R>
void ISend
( const R* buf, int count, int to, int tag, Comm comm, Request& request );
template<typename R>
void ISend
( const Complex<R>* buf, int count, int to, int tag, Comm comm, 
  Request& request );

template<typename R>
void ISSend
( const R* buf, int count, int to, int tag, Comm comm, Request& request );
template<typename R>
void ISSend
( const Complex<R>* buf, int count, int to, int tag, Comm comm, 
  Request& request );

template<typename R>
void Recv( R* buf, int count, int from, int tag, Comm comm );
template<typename R>
void Recv( Complex<R>* buf, int count, int from, int tag, Comm comm );

template<typename R>
void IRecv
( R* buf, int count, int from, int tag, Comm comm, Request& request );
template<typename R>
void IRecv
( Complex<R>* buf, int count, int from, int tag, Comm comm, Request& request );


template<typename R>
void SendRecv
( const R* sbuf, int sc, int to,   int stag,
        R* rbuf, int rc, int from, int rtag, Comm comm );
template<typename R>
void SendRecv
( const Complex<R>* sbuf, int sc, int to,   int stag,
        Complex<R>* rbuf, int rc, int from, int rtag, Comm comm );

// Collective communication

template<typename R>
void Broadcast( R* buf, int count, int root, Comm comm );
template<typename R>
void Broadcast( Complex<R>* buf, int count, int root, Comm comm );

#ifdef HAVE_NONBLOCKING_COLLECTIVES
template<typename R>
void IBroadcast
( R* buf, int count, int root, Comm comm, Request& request );
template<typename R>
void IBroadcast
( Complex<R>* buf, int count, int root, Comm comm, Request& request );
#endif

template<typename R>
void Gather
( const R* sbuf, int sc,
        R* rbuf, int rc, int root, Comm comm );
template<typename R>
void  Gather
( const Complex<R>* sbuf, int sc,
        Complex<R>* rbuf, int rc, int root, Comm comm );

#ifdef HAVE_NONBLOCKING_COLLECTIVES
template<typename R>
void IGather
( const R* sbuf, int sc,
        R* rbuf, int rc, int root, Comm comm, Request& request );
template<typename R>
void IGather
( const Complex<R>* sbuf, int sc,
        Complex<R>* rbuf, int rc, int root, Comm comm, Request& request );
#endif

template<typename R>
void Gather
( const R* sbuf, int sc,
        R* rbuf, const int* rcs, const int* rds, int root, Comm comm );
template<typename R>
void Gather
( const Complex<R>* sbuf, int sc,
        Complex<R>* rbuf, const int* rcs, const int* rds, 
  int root, Comm comm );

template<typename R>
void AllGather
( const R* sbuf, int sc,
        R* rbuf, int rc, Comm comm );
template<typename R>
void AllGather
( const Complex<R>* sbuf, int sc,
        Complex<R>* rbuf, int rc, Comm comm );

template<typename R>
void AllGather
( const R* sbuf, int sc,
        R* rbuf, const int* rcs, const int* rds, Comm comm );
template<typename R>
void AllGather
( const Complex<R>* sbuf, int sc,
        Complex<R>* rbuf, const int* rcs, const int* rds, Comm comm );

template<typename R>
void Scatter
( const R* sbuf, int sc,
        R* rbuf, int rc, int root, Comm comm );
template<typename R>
void Scatter
( const Complex<R>* sbuf, int sc,
        Complex<R>* rbuf, int rc, int root, Comm comm );

// In-place option
template<typename R>
void Scatter( R* buf, int sc, int rc, int root, Comm comm );
template<typename R>
void Scatter( Complex<R>* buf, int sc, int rc, int root, Comm comm );

template<typename R>
void AllToAll
( const R* sbuf, int sc,
        R* rbuf, int rc, Comm comm );
template<typename R>
void AllToAll
( const Complex<R>* sbuf, int sc,
        Complex<R>* rbuf, int rc, Comm comm );

template<typename R>
void AllToAll
( const R* sbuf, const int* scs, const int* sds,
        R* rbuf, const int* rcs, const int* rds, Comm comm );
template<typename R>
void AllToAll
( const Complex<R>* sbuf, const int* scs, const int* sds,
        Complex<R>* rbuf, const int* rcs, const int* rds, Comm comm );

template<typename R>
void Reduce
( const R* sbuf, R* rbuf, int count, Op op, int root, Comm comm );
template<typename R>
void Reduce
( const Complex<R>* sbuf, Complex<R>* rbuf, int count, Op op, 
  int root, Comm comm );

// In-place option
template<typename R>
void Reduce( R* buf, int count, Op op, int root, Comm comm );
template<typename R>
void Reduce( Complex<R>* buf, int count, Op op, int root, Comm comm );

template<typename R>
void AllReduce
( const R* sbuf, R* rbuf, int count, Op op, Comm comm );
template<typename R>
void AllReduce
( const Complex<R>* sbuf, Complex<R>* rbuf, int count, Op op, Comm comm );

// In-place option
template<typename R>
void AllReduce( R* buf, int count, Op op, Comm comm );
template<typename R>
void AllReduce( Complex<R>* buf, int count, Op op, Comm comm );

template<typename R>
void ReduceScatter
( R* sbuf, R* rbuf, int rc, Op op, Comm comm );
template<typename R>
void ReduceScatter
( Complex<R>* sbuf, Complex<R>* rbuf, int rc, Op op, Comm comm );

// In-place option
template<typename R>
void ReduceScatter( R* buf, int rc, Op op, Comm comm );
template<typename R>
void ReduceScatter( Complex<R>* buf, int rc, Op op, Comm comm );

template<typename R>
void ReduceScatter
( const R* sbuf, R* rbuf, const int* rcs, Op op, Comm comm );
template<typename R>
void ReduceScatter
( const Complex<R>* sbuf, Complex<R>* rbuf, const int* rcs, Op op, Comm comm );

} // mpi
} // elem

#endif // ifndef CORE_MPI_HPP
