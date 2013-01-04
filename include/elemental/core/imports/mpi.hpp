/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

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
template<> int GetCount<byte>( Status& status );
template<> int GetCount<int>( Status& status );
template<> int GetCount<float>( Status& status );
template<> int GetCount<double>( Status& status );
template<> int GetCount<scomplex>( Status& status );
template<> int GetCount<dcomplex>( Status& status );

// Point-to-point communication

void Send( const byte* buf, int count, int to, int tag, Comm comm );
void Send( const int* buf, int count, int to, int tag, Comm comm );
void Send( const float* buf, int count, int to, int tag, Comm comm );
void Send( const double* buf, int count, int to, int tag, Comm comm );
void Send( const scomplex* buf, int count, int to, int tag, Comm comm );
void Send( const dcomplex* buf, int count, int to, int tag, Comm comm );

void ISend
( const byte* buf, int count, int to, int tag, Comm comm, Request& request );
void ISend
( const int* buf, int count, int to, int tag, Comm comm, Request& request );
void ISend
( const float* buf, int count, int to, int tag, Comm comm, Request& request );
void ISend
( const double* buf, int count, int to, int tag, Comm comm, Request& request );
void ISend
( const scomplex* buf, int count, int to, int tag, Comm comm, 
  Request& request );
void ISend
( const dcomplex* buf, int count, int to, int tag, Comm comm, 
  Request& request );

void ISSend
( const byte* buf, int count, int to, int tag, Comm comm, Request& request );
void ISSend
( const int* buf, int count, int to, int tag, Comm comm, Request& request );
void ISSend
( const float* buf, int count, int to, int tag, Comm comm, Request& request );
void ISSend
( const double* buf, int count, int to, int tag, Comm comm, Request& request );
void ISSend
( const scomplex* buf, int count, int to, int tag, Comm comm, 
  Request& request );
void ISSend
( const dcomplex* buf, int count, int to, int tag, Comm comm, 
 Request& request );
 
void Recv( byte* buf, int count, int from, int tag, Comm comm );
void Recv( int* buf, int count, int from, int tag, Comm comm );
void Recv( float* buf, int count, int from, int tag, Comm comm );
void Recv( double* buf, int count, int from, int tag, Comm comm );
void Recv( scomplex* buf, int count, int from, int tag, Comm comm );
void Recv( dcomplex* buf, int count, int from, int tag, Comm comm );

void IRecv
( byte* buf, int count, int from, int tag, Comm comm, Request& request );
void IRecv
( int* buf, int count, int from, int tag, Comm comm, Request& request );
void IRecv
( float* buf, int count, int from, int tag, Comm comm, Request& request );
void IRecv
( double* buf, int count, int from, int tag, Comm comm, Request& request );
void IRecv
( scomplex* buf, int count, int from, int tag, Comm comm, Request& request );
void IRecv
( dcomplex* buf, int count, int from, int tag, Comm comm, Request& request );

void SendRecv
( const byte* sbuf, int sc, int to,   int stag,
        byte* rbuf, int rc, int from, int rtag, Comm comm );
void SendRecv
( const int* sbuf, int sc, int to,   int stag,
        int* rbuf, int rc, int from, int rtag, Comm comm );
void SendRecv
( const float* sbuf, int sc, int to,   int stag,
        float* rbuf, int rc, int from, int rtag, Comm comm );
void SendRecv
( const double* sbuf, int sc, int to,   int stag,
        double* rbuf, int rc, int from, int rtag, Comm comm );
void SendRecv
( const scomplex* sbuf, int sc, int to,   int stag,
        scomplex* rbuf, int rc, int from, int rtag, Comm comm );
void SendRecv
( const dcomplex* sbuf, int sc, int to,   int stag,
        dcomplex* rbuf, int rc, int from, int rtag, Comm comm );

// Collective communication

void Broadcast( byte* buf, int count, int root, Comm comm );
void Broadcast( int* buf, int count, int root, Comm comm );
void Broadcast( float* buf, int count, int root, Comm comm );
void Broadcast( double* buf, int count, int root, Comm comm );
void Broadcast( scomplex* buf, int count, int root, Comm comm );
void Broadcast( dcomplex* buf, int count, int root, Comm comm );

#ifdef HAVE_NONBLOCKING_COLLECTIVES
void IBroadcast
( byte* buf, int count, int root, Comm comm, Request& request );
void IBroadcast
( int* buf, int count, int root, Comm comm, Request& request );
void IBroadcast
( float* buf, int count, int root, Comm comm, Request& request );
void IBroadcast
( double* buf, int count, int root, Comm comm, Request& request );
void IBroadcast
( scomplex* buf, int count, int root, Comm comm, Request& request );
void IBroadcast
( dcomplex* buf, int count, int root, Comm comm, Request& request );
#endif

void Gather
( const byte* sbuf, int sc,
        byte* rbuf, int rc, int root, Comm comm );
void Gather
( const int* sbuf, int sc,
        int* rbuf, int rc, int root, Comm comm );
void Gather
( const float* sbuf, int sc,
        float* rbuf, int rc, int root, Comm comm );
void Gather
( const double* sbuf, int sc,
        double* rbuf, int rc, int root, Comm comm );
void  Gather
( const scomplex* sbuf, int sc,
        scomplex* rbuf, int rc, int root, Comm comm );
void Gather
( const dcomplex* sbuf, int sc,
        dcomplex* rbuf, int rc, int root, Comm comm );

#ifdef HAVE_NONBLOCKING_COLLECTIVES
void IGather
( const byte* sbuf, int sc,
        byte* rbuf, int rc, int root, Comm comm, Request& request );
void IGather
( const int* sbuf, int sc,
        int* rbuf, int rc, int root, Comm comm, Request& request );
void IGather
( const float* sbuf, int sc,
        float* rbuf, int rc, int root, Comm comm, Request& request );
void IGather
( const double* sbuf, int sc,
        double* rbuf, int rc, int root, Comm comm, Request& request );
void IGather
( const scomplex* sbuf, int sc,
        scomplex* rbuf, int rc, int root, Comm comm, Request& request );
void IGather
( const dcomplex* sbuf, int sc,
        dcomplex* rbuf, int rc, int root, Comm comm, Request& request );
#endif

void Gather
( const byte* sbuf, int sc,
        byte* rbuf, const int* rcs, const int* rds, int root, Comm comm );
void Gather
( const int* sbuf, int sc,
        int* rbuf, const int* rcs, const int* rds, int root, Comm comm );
void Gather
( const float* sbuf, int sc,
        float* rbuf, const int* rcs, const int* rds, int root, Comm comm );
void Gather
( const double* sbuf, int sc,
        double* rbuf, const int* rcs, const int* rds, int root, Comm comm );
void Gather
( const scomplex* sbuf, int sc,
        scomplex* rbuf, const int* rcs, const int* rds, int root, Comm comm );
void Gather
( const dcomplex* sbuf, int sc,
        dcomplex* rbuf, const int* rcs, const int* rds, int root, Comm comm );

#ifdef HAVE_NONBLOCKING_COLLECTIVES
void IGather
( const byte* sbuf, int sc,
        byte* rbuf, const int* rcs, const int* rds, 
        int root, Comm comm, Request& request );
void IGather
( const int* sbuf, int sc,
        int* rbuf, const int* rcs, const int* rds, 
        int root, Comm comm, Request& request );
void IGather
( const float* sbuf, int sc,
        float* rbuf, const int* rcs, const int* rds, 
        int root, Comm comm, Request& request );
void IGather
( const double* sbuf, int sc,
        double* rbuf, const int* rcs, const int* rds, 
        int root, Comm comm, Request& request );
void IGather
( const scomplex* sbuf, int sc,
        scomplex* rbuf, const int* rcs, const int* rds, 
        int root, Comm comm, Request& request );
void IGather
( const dcomplex* sbuf, int sc,
        dcomplex* rbuf, const int* rcs, const int* rds, 
        int root, Comm comm, Request& request );
#endif

void AllGather
( const byte* sbuf, int sc,
        byte* rbuf, int rc, Comm comm );
void AllGather
( const int* sbuf, int sc,
        int* rbuf, int rc, Comm comm );
void AllGather
( const float* sbuf, int sc,
        float* rbuf, int rc, Comm comm );
void AllGather
( const double* sbuf, int sc,
        double* rbuf, int rc, Comm comm );
void AllGather
( const scomplex* sbuf, int sc,
        scomplex* rbuf, int rc, Comm comm );
void AllGather
( const dcomplex* sbuf, int sc,
        dcomplex* rbuf, int rc, Comm comm );

void AllGather
( const byte* sbuf, int sc,
        byte* rbuf, const int* rcs, const int* rds, Comm comm );
void AllGather
( const int* sbuf, int sc,
        int* rbuf, const int* rcs, const int* rds, Comm comm );
void AllGather
( const float* sbuf, int sc,
        float* rbuf, const int* rcs, const int* rds, Comm comm );
void AllGather
( const double* sbuf, int sc,
        double* rbuf, const int* rcs, const int* rds, Comm comm );
void AllGather
( const scomplex* sbuf, int sc,
        scomplex* rbuf, const int* rcs, const int* rds, Comm comm );
void AllGather
( const dcomplex* sbuf, int sc,
        dcomplex* rbuf, const int* rcs, const int* rds, Comm comm );

#ifdef HAVE_NONBLOCKING_COLLECTIVES
void IAllGather
( const byte* sbuf, int sc,
        byte* rbuf, int rc, Comm comm, Request& request );
void IAllGather
( const int* sbuf, int sc,
        int* rbuf, int rc, Comm comm, Request& request );
void IAllGather
( const float* sbuf, int sc,
        float* rbuf, int rc, Comm comm, Request& request );
void IAllGather
( const double* sbuf, int sc,
        double* rbuf, int rc, Comm comm, Request& request );
void IAllGather
( const scomplex* sbuf, int sc,
        scomplex* rbuf, int rc, Comm comm, Request& request );
void IAllGather
( const dcomplex* sbuf, int sc,
        dcomplex* rbuf, int rc, Comm comm, Request& request );
#endif

void Scatter
( const byte* sbuf, int sc,
        byte* rbuf, int rc, int root, Comm comm );
void Scatter
( const int* sbuf, int sc,
        int* rbuf, int rc, int root, Comm comm );
void Scatter
( const float* sbuf, int sc,
        float* rbuf, int rc, int root, Comm comm );
void Scatter
( const double* sbuf, int sc,
        double* rbuf, int rc, int root, Comm comm );
void Scatter
( const scomplex* sbuf, int sc,
        scomplex* rbuf, int rc, int root, Comm comm );
void Scatter
( const dcomplex* sbuf, int sc,
        dcomplex* rbuf, int rc, int root, Comm comm );

#ifdef HAVE_NONBLOCKING_COLLECTIVES
void IScatter
( const byte* sbuf, int sc,
        byte* rbuf, int rc, int root, Comm comm, Request& request );
void IScatter
( const int* sbuf, int sc,
        int* rbuf, int rc, int root, Comm comm, Request& request );
void IScatter
( const float* sbuf, int sc,
        float* rbuf, int rc, int root, Comm comm, Request& request );
void IScatter
( const double* sbuf, int sc,
        double* rbuf, int rc, int root, Comm comm, Request& request );
void IScatter
( const scomplex* sbuf, int sc,
        scomplex* rbuf, int rc, int root, Comm comm, Request& request );
void IScatter
( const dcomplex* sbuf, int sc,
        dcomplex* rbuf, int rc, int root, Comm comm, Request& request );
#endif

// In-place option
void Scatter( byte* buf, int sc, int rc, int root, Comm comm );
void Scatter( int* buf, int sc, int rc, int root, Comm comm );
void Scatter( float* buf, int sc, int rc, int root, Comm comm );
void Scatter( double* buf, int sc, int rc, int root, Comm comm );
void Scatter( scomplex* buf, int sc, int rc, int root, Comm comm );
void Scatter( dcomplex* buf, int sc, int rc, int root, Comm comm );

#ifdef HAVE_NONBLOCKING_COLLECTIVES
// TODO: Nonblocking in-place Scatter?
#endif
 
void AllToAll
( const byte* sbuf, int sc,
        byte* rbuf, int rc, Comm comm );
void AllToAll
( const int* sbuf, int sc,
        int* rbuf, int rc, Comm comm );
void AllToAll
( const float* sbuf, int sc,
        float* rbuf, int rc, Comm comm );
void AllToAll
( const double* sbuf, int sc,
        double* rbuf, int rc, Comm comm );
void AllToAll
( const scomplex* sbuf, int sc,
        scomplex* rbuf, int rc, Comm comm );
void AllToAll
( const dcomplex* sbuf, int sc,
        dcomplex* rbuf, int rc, Comm comm );

#ifdef HAVE_NONBLOCKING_COLLECTIVES
void IAllToAll
( const byte* sbuf, int sc,
        byte* rbuf, int rc, Comm comm, Request& request );
void IAllToAll
( const int* sbuf, int sc,
        int* rbuf, int rc, Comm comm, Request& request );
void IAllToAll
( const float* sbuf, int sc,
        float* rbuf, int rc, Comm comm, Request& request );
void IAllToAll
( const double* sbuf, int sc,
        double* rbuf, int rc, Comm comm, Request& request );
void IAllToAll
( const scomplex* sbuf, int sc,
        scomplex* rbuf, int rc, Comm comm, Request& request );
void IAllToAll
( const dcomplex* sbuf, int sc,
        dcomplex* rbuf, int rc, Comm comm, Request& request );
#endif

void AllToAll
( const byte* sbuf, const int* scs, const int* sds,
        byte* rbuf, const int* rcs, const int* rds, Comm comm );
void AllToAll
( const int* sbuf, const int* scs, const int* sds,
        int* rbuf, const int* rcs, const int* rds, Comm comm );
void AllToAll
( const float* sbuf, const int* scs, const int* sds,
        float* rbuf, const int* rcs, const int* rds, Comm comm );
void AllToAll
( const double* sbuf, const int* scs, const int* sds,
        double* rbuf, const int* rcs, const int* rds, Comm comm );
void AllToAll
( const scomplex* sbuf, const int* scs, const int* sds,
        scomplex* rbuf, const int* rcs, const int* rds, Comm comm );
void AllToAll
( const dcomplex* sbuf, const int* scs, const int* sds,
        dcomplex* rbuf, const int* rcs, const int* rds, Comm comm );

#ifdef HAVE_NONBLOCKING_COLLECTIVES
void IAllToAll
( const byte* sbuf, const int* scs, const int* sds,
        byte* rbuf, const int* rcs, const int* rds, 
        Comm comm, Request& request );
void IAllToAll
( const int* sbuf, const int* scs, const int* sds,
        int* rbuf, const int* rcs, const int* rds, 
        Comm comm, Request& request );
void IAllToAll
( const float* sbuf, const int* scs, const int* sds,
        float* rbuf, const int* rcs, const int* rds, 
        Comm comm, Request& request );
void IAllToAll
( const double* sbuf, const int* scs, const int* sds,
        double* rbuf, const int* rcs, const int* rds, 
        Comm comm, Request& request );
void IAllToAll
( const scomplex* sbuf, const int* scs, const int* sds,
        scomplex* rbuf, const int* rcs, const int* rds, 
        Comm comm, Request& request );
void IAllToAll
( const dcomplex* sbuf, const int* scs, const int* sds,
        dcomplex* rbuf, const int* rcs, const int* rds, 
        Comm comm, Request& request );
#endif

void Reduce
( const byte* sbuf, byte* rbuf, int count, Op op, int root, Comm comm );
void Reduce
( const int* sbuf, int* rbuf, int count, Op op, int root, Comm comm );
void Reduce
( const float* sbuf, float* rbuf, int count, Op op, int root, Comm comm );
void Reduce
( const double* sbuf, double* rbuf, int count, Op op, int root, Comm comm );
void Reduce
( const scomplex* sbuf, scomplex* rbuf, int count, Op op, int root, Comm comm );
void Reduce
( const dcomplex* sbuf, dcomplex* rbuf, int count, Op op, int root, Comm comm );

#ifdef HAVE_NONBLOCKING_COLLECTIVES
void IReduce
( const byte* sbuf, byte* rbuf, int count, 
  Op op, int root, Comm comm, Request& request );
void IReduce
( const int* sbuf, int* rbuf, int count, 
  Op op, int root, Comm comm, Request& request );
void IReduce
( const float* sbuf, float* rbuf, int count, 
  Op op, int root, Comm comm, Request& request );
void IReduce
( const double* sbuf, double* rbuf, int count, 
  Op op, int root, Comm comm, Request& request );
void IReduce
( const scomplex* sbuf, scomplex* rbuf, int count, 
  Op op, int root, Comm comm, Request& request );
void IReduce
( const dcomplex* sbuf, dcomplex* rbuf, int count, 
  Op op, int root, Comm comm, Request& request );
#endif

// In-place option
void Reduce( byte* buf, int count, Op op, int root, Comm comm );
void Reduce( int* buf, int count, Op op, int root, Comm comm );
void Reduce( float* buf, int count, Op op, int root, Comm comm );
void Reduce( double* buf, int count, Op op, int root, Comm comm );
void Reduce( scomplex* buf, int count, Op op, int root, Comm comm );
void Reduce( dcomplex* buf, int count, Op op, int root, Comm comm );

#ifdef HAVE_NONBLOCKING_COLLECTIVES
// TODO: Add in-place nonblocking Reduce?
#endif
    
void AllReduce
( const byte* sbuf, byte* rbuf, int count, Op op, Comm comm );
void AllReduce
( const int* sbuf, int* rbuf, int count, Op op, Comm comm );
void AllReduce
( const float* sbuf, float* rbuf, int count, Op op, Comm comm );
void AllReduce
( const double* sbuf, double* rbuf, int count, Op op, Comm comm );
void AllReduce
( const scomplex* sbuf, scomplex* rbuf, int count, Op op, Comm comm );
void AllReduce
( const dcomplex* sbuf, dcomplex* rbuf, int count, Op op, Comm comm );

#ifdef HAVE_NONBLOCKING_COLLECTIVES
void IAllReduce
( const byte* sbuf, byte* rbuf, int count, 
  Op op, Comm comm, Request& request );
void IAllReduce
( const int* sbuf, int* rbuf, int count, 
  Op op, Comm comm, Request& request );
void IAllReduce
( const float* sbuf, float* rbuf, int count, 
  Op op, Comm comm, Request& request );
void IAllReduce
( const double* sbuf, double* rbuf, int count, 
  Op op, Comm comm, Request& request );
void IAllReduce
( const scomplex* sbuf, scomplex* rbuf, int count, 
  Op op, Comm comm, Request& request );
void IAllReduce
( const dcomplex* sbuf, dcomplex* rbuf, int count, 
  Op op, Comm comm, Request& request );
#endif

// In-place option
void AllReduce( byte* buf, int count, Op op, Comm comm );
void AllReduce( int* buf, int count, Op op, Comm comm );
void AllReduce( float* buf, int count, Op op, Comm comm );
void AllReduce( double* buf, int count, Op op, Comm comm );
void AllReduce( scomplex* buf, int count, Op op, Comm comm );
void AllReduce( dcomplex* buf, int count, Op op, Comm comm );

#ifdef HAVE_NONBLOCKING_COLLECTIVES
// TODO: Add in-place nonblocking AllReduce?
#endif

void ReduceScatter( byte* sbuf, byte* rbuf, int rc, Op op, Comm comm );
void ReduceScatter( int* sbuf, int* rbuf, int rc, Op op, Comm comm );
void ReduceScatter( float* sbuf, float* rbuf, int rc, Op op, Comm comm );
void ReduceScatter( double* sbuf, double* rbuf, int rc, Op op, Comm comm );
void ReduceScatter( scomplex* sbuf, scomplex* rbuf, int rc, Op op, Comm comm );
void ReduceScatter( dcomplex* sbuf, dcomplex* rbuf, int rc, Op op, Comm comm );

#ifdef HAVE_NONBLOCKING_COLLECTIVES
void IReduceScatter
( byte* sbuf, byte* rbuf, int rc, Op op, Comm comm, Request& request );
void IReduceScatter
( int* sbuf, int* rbuf, int rc, Op op, Comm comm, Request& request );
void IReduceScatter
( float* sbuf, float* rbuf, int rc, Op op, Comm comm, Request& request );
void IReduceScatter
( double* sbuf, double* rbuf, int rc, Op op, Comm comm, Request& request );
void IReduceScatter
( scomplex* sbuf, scomplex* rbuf, int rc, Op op, Comm comm, Request& request );
void IReduceScatter
( dcomplex* sbuf, dcomplex* rbuf, int rc, Op op, Comm comm, Request& request );
#endif

// In-place option
void ReduceScatter( byte* buf, int rc, Op op, Comm comm );
void ReduceScatter( int* buf, int rc, Op op, Comm comm );
void ReduceScatter( float* buf, int rc, Op op, Comm comm );
void ReduceScatter( double* buf, int rc, Op op, Comm comm );
void ReduceScatter( scomplex* buf, int rc, Op op, Comm comm );
void ReduceScatter( dcomplex* buf, int rc, Op op, Comm comm );

#ifdef HAVE_NONBLOCKING_COLLECTIVES
// TODO: Add in-place nonblocking ReduceScatter?
#endif

void ReduceScatter
( const byte* sbuf, byte* rbuf, const int* rcs, Op op, Comm comm );
void ReduceScatter
( const int* sbuf, int* rbuf, const int* rcs, Op op, Comm comm );
void ReduceScatter
( const float* sbuf, float* rbuf, const int* rcs, Op op, Comm comm );
void ReduceScatter
( const double* sbuf, double* rbuf, const int* rcs, Op op, Comm comm );
void ReduceScatter
( const scomplex* sbuf, scomplex* rbuf, const int* rcs, Op op, Comm comm );
void ReduceScatter
( const dcomplex* sbuf, dcomplex* rbuf, const int* rcs, Op op, Comm comm );

#ifdef HAVE_NONBLOCKING_COLLECTIVES
void IReduceScatter
( const byte* sbuf, byte* rbuf, const int* rcs, 
  Op op, Comm comm, Request& request );
void IReduceScatter
( const int* sbuf, int* rbuf, const int* rcs, 
  Op op, Comm comm, Request& request );
void IReduceScatter
( const float* sbuf, float* rbuf, const int* rcs, 
  Op op, Comm comm, Request& request );
void IReduceScatter
( const double* sbuf, double* rbuf, const int* rcs, 
  Op op, Comm comm, Request& request );
void IReduceScatter
( const scomplex* sbuf, scomplex* rbuf, const int* rcs, 
  Op op, Comm comm, Request& request );
void IReduceScatter
( const dcomplex* sbuf, dcomplex* rbuf, const int* rcs, 
  Op op, Comm comm, Request& request );
#endif

} // mpi
} // elem
