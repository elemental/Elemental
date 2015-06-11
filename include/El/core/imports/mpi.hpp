/*
   Copyright (c) 2009-2015, Jack Poulson
                      2013, Jeff Hammond
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_IMPORTS_MPI_HPP
#define EL_IMPORTS_MPI_HPP

namespace El {

namespace mpi {

#if defined(EL_HAVE_MPI3_NONBLOCKING_COLLECTIVES) || \
    defined(EL_HAVE_MPIX_NONBLOCKING_COLLECTIVES)
#define EL_HAVE_NONBLOCKING 1
#else
#define EL_HAVE_NONBLOCKING 0
#endif

#ifdef EL_HAVE_NONBLOCKING_COLLECTIVES
#ifdef EL_HAVE_MPI3_NONBLOCKING_COLLECTIVES
#define EL_NONBLOCKING_COLL(name) MPI_ ## name
#else
#define EL_NONBLOCKING_COLL(name) MPIX_ ## name
#endif
#endif

struct Comm
{
    MPI_Comm comm;
    Comm( MPI_Comm mpiComm=MPI_COMM_NULL ) : comm(mpiComm) { }
};
inline bool operator==( const Comm& a, const Comm& b )
{ return a.comm == b.comm; }
inline bool operator!=( const Comm& a, const Comm& b )
{ return a.comm != b.comm; }

struct Group
{
    MPI_Group group;
    Group( MPI_Group mpiGroup=MPI_GROUP_NULL ) : group(mpiGroup) { }
};
inline bool operator==( const Group& a, const Group& b )
{ return a.group == b.group; }
inline bool operator!=( const Group& a, const Group& b )
{ return a.group != b.group; }

struct Op
{
    MPI_Op op;
    Op( MPI_Op mpiOp=MPI_OP_NULL ) : op(mpiOp) { }
};
inline bool operator==( const Op& a, const Op& b )
{ return a.op == b.op; }
inline bool operator!=( const Op& a, const Op& b )
{ return a.op != b.op; }

// Datatype definitions
// TODO: Convert these to structs/classes
typedef MPI_Aint Aint;
typedef MPI_Datatype Datatype;
typedef MPI_Errhandler ErrorHandler;
typedef MPI_Request Request;
typedef MPI_Status Status;
typedef MPI_User_function UserFunction;

// Standard constants
const int ANY_SOURCE = MPI_ANY_SOURCE;
const int ANY_TAG = MPI_ANY_TAG;
#ifdef EL_HAVE_MPI_QUERY_THREAD
const int THREAD_SINGLE = MPI_THREAD_SINGLE;
const int THREAD_FUNNELED = MPI_THREAD_FUNNELED;
const int THREAD_SERIALIZED = MPI_THREAD_SERIALIZED;
const int THREAD_MULTIPLE = MPI_THREAD_MULTIPLE;
#else
const int THREAD_SINGLE = 0;
const int THREAD_FUNNELED = 1;
const int THREAD_SERIALIZED = 2;
const int THREAD_MULTIPLE = 3;
#endif
const int UNDEFINED = MPI_UNDEFINED;
const Group GROUP_NULL = MPI_GROUP_NULL;
const Comm COMM_NULL = MPI_COMM_NULL;
const Comm COMM_SELF = MPI_COMM_SELF;
const Comm COMM_WORLD = MPI_COMM_WORLD;
const ErrorHandler ERRORS_RETURN = MPI_ERRORS_RETURN;
const ErrorHandler ERRORS_ARE_FATAL = MPI_ERRORS_ARE_FATAL;
const Group GROUP_EMPTY = MPI_GROUP_EMPTY;
const Request REQUEST_NULL = MPI_REQUEST_NULL;
const Op MAX = MPI_MAX;
const Op MIN = MPI_MIN;
const Op MAXLOC = MPI_MAXLOC;
const Op MINLOC = MPI_MINLOC;
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
inline int Pad( int count ) { return std::max(count,MIN_COLL_MSG); }

bool CommSameSizeAsInteger();
bool GroupSameSizeAsInteger();

// Environment routines
void Initialize( int& argc, char**& argv );
int InitializeThread( int& argc, char**& argv, int required );
void Finalize();
bool Initialized();
bool Finalized();
int QueryThread();
void Abort( Comm comm, int errCode );
double Time();
void Create( UserFunction* func, bool commutes, Op& op );
void Free( Op& op );
void Free( Datatype& type );

// Communicator manipulation
int WorldRank();
int WorldSize();
int Rank( Comm comm );
int Size( Comm comm );
void Create( Comm parentComm, Group subsetGroup, Comm& subsetComm );
void Dup( Comm original, Comm& duplicate );
void Split( Comm comm, int color, int key, Comm& newComm );
void Free( Comm& comm );
bool Congruent( Comm comm1, Comm comm2 );
void ErrorHandlerSet( Comm comm, ErrorHandler errorHandler );

// Cartesian communicator routines
void CartCreate
( Comm comm, int numDims, const int* dimensions, const int* periods, 
  bool reorder, Comm& cartComm );
void CartSub
( Comm comm, const int* remainingDims, Comm& subComm );

// Group manipulation
int Rank( Group group );
int Size( Group group );
void CommGroup( Comm comm, Group& group );
void Dup( Group group, Group& newGroup );
void Union( Group groupA, Group groupB, Group& newGroup );
void Incl( Group group, int n, const int* ranks, Group& subGroup );
void Excl( Group group, int n, const int* ranks, Group& subGroup );
void Difference( Group parent, Group subset, Group& complement );
void Free( Group& group );
int Translate( Group origGroup, int origRank, Group newGroup );
int Translate( Comm  origComm,  int origRank, Group newGroup );
int Translate( Group origGroup, int origRank, Comm  newComm  );
int Translate( Comm  origComm,  int origRank, Comm  newComm  );
void Translate
( Group origGroup, int size, const int* origRanks, 
  Group newGroup,                  int* newRanks );
void Translate
( Comm origComm,  int size, const int* origRanks, 
  Group newGroup,                 int* newRanks );
void Translate
( Group origGroup, int size, const int* origRanks, 
  Comm newComm,                    int* newRanks );
void Translate
( Comm origComm, int size, const int* origRanks, 
  Comm newComm,                  int* newRanks );

// Utilities
void Barrier( Comm comm );
void Wait( Request& request );
void Wait( Request& request, Status& status );
void WaitAll( int numRequests, Request* requests );
void WaitAll( int numRequests, Request* requests, Status* statuses );
bool Test( Request& request );
bool IProbe( int source, int tag, Comm comm, Status& status );

template<typename T>
int GetCount( Status& status );

// Point-to-point communication
// ============================

// Send
// ----
template<typename Real>
void TaggedSend
( const Real* buf, int count, int to, int tag, Comm comm );
template<typename Real>
void TaggedSend
( const Complex<Real>* buf, int count, int to, int tag, Comm comm );

// If the tag is irrelevant
template<typename T>
void Send( const T* buf, int count, int to, Comm comm );

// If the send-count is one
template<typename T>
void TaggedSend( T b, int to, int tag, Comm comm );

// If the send-count is one and the tag is irrelevant
template<typename T>
void Send( T b, int to, Comm comm );

// Non-blocking send
// -----------------
template<typename Real>
void TaggedISend
( const Real* buf, int count, int to, int tag, Comm comm, Request& request );
template<typename Real>
void TaggedISend
( const Complex<Real>* buf, int count, int to, int tag, Comm comm, 
  Request& request );

// If the tag is irrelevant
template<typename T>
void ISend( const T* buf, int count, int to, Comm comm, Request& request );

// If the send count is one
template<typename T>
void TaggedISend( T b, int to, int tag, Comm comm, Request& request );

// If the send count is one and the tag is irrelevant
template<typename T>
void ISend( T b, int to, Comm comm, Request& request );

// Non-blocking synchronous Send
// -----------------------------
template<typename Real>
void TaggedISSend
( const Real* buf, int count, int to, int tag, Comm comm, Request& request );
template<typename Real>
void TaggedISSend
( const Complex<Real>* buf, int count, int to, int tag, Comm comm, 
  Request& request );

// If the tag is irrelevant
template<typename T>
void ISSend( const T* buf, int count, int to, Comm comm, Request& request );

// If the send count is one
template<typename T>
void TaggedISSend( T b, int to, int tag, Comm comm, Request& request );

// If the send count is one and the tag is irrelevant
template<typename T>
void ISSend( T b, int to, Comm comm, Request& request );

// Recv
// ----
template<typename Real>
void TaggedRecv( Real* buf, int count, int from, int tag, Comm comm );
template<typename Real>
void TaggedRecv( Complex<Real>* buf, int count, int from, int tag, Comm comm );

// If the tag is irrelevant
template<typename T>
void Recv( T* buf, int count, int from, Comm comm );

// If the recv count is one
template<typename T>
T TaggedRecv( int from, int tag, Comm comm );

// If the recv count is one and the tag is irrelevant
template<typename T>
T Recv( int from, Comm comm );

// Non-blocking recv
// -----------------
template<typename Real>
void TaggedIRecv
( Real* buf, int count, int from, int tag, Comm comm, 
  Request& request );
template<typename Real>
void TaggedIRecv
( Complex<Real>* buf, int count, int from, int tag, Comm comm, 
  Request& request );

// If the tag is irrelevant
template<typename T>
void IRecv( T* buf, int count, int from, Comm comm, Request& request );

// If the recv count is one
template<typename T>
T TaggedIRecv( int from, int tag, Comm comm, Request& request );

// If the recv count is one and the tag is irrelevant
template<typename T>
T IRecv( int from, Comm comm, Request& request );

// SendRecv
// --------
template<typename Real>
void TaggedSendRecv
( const Real* sbuf, int sc, int to,   int stag,
        Real* rbuf, int rc, int from, int rtag, Comm comm );
template<typename Real>
void TaggedSendRecv
( const Complex<Real>* sbuf, int sc, int to,   int stag,
        Complex<Real>* rbuf, int rc, int from, int rtag, Comm comm );

// If the tags are irrelevant
template<typename T>
void SendRecv
( const T* sbuf, int sc, int to,
        T* rbuf, int rc, int from, Comm comm );

// If the send and recv counts are one
template<typename T>
T TaggedSendRecv( T sb, int to, int stag, int from, int rtag, Comm comm );

// If the send and recv counts are one and the tags don't matter
template<typename T>
T SendRecv( T sb, int to, int from, Comm comm );

// Single-buffer SendRecv
// ----------------------
template<typename Real>
void TaggedSendRecv
( Real* buf, int count, int to, int stag, int from, int rtag, Comm comm );
template<typename Real>
void TaggedSendRecv
( Complex<Real>* buf, int count, int to, int stag, int from, int rtag, 
  Comm comm );

// If the tags don't matter
template<typename T>
void SendRecv( T* buf, int count, int to, int from, Comm comm );

// Collective communication
// ========================

// Broadcast
// ---------
template<typename Real>
void Broadcast( Real* buf, int count, int root, Comm comm );
template<typename Real>
void Broadcast( Complex<Real>* buf, int count, int root, Comm comm );

// If the message length is one
template<typename T>
void Broadcast( T& b, int root, Comm comm );

// Non-blocking broadcast
// ----------------------
template<typename Real>
void IBroadcast
( Real* buf, int count, int root, Comm comm, Request& request );
template<typename Real>
void IBroadcast
( Complex<Real>* buf, int count, int root, Comm comm, Request& request );

// If the message length is one
template<typename T>
void IBroadcast( T& b, int root, Comm comm, Request& request );

// Gather
// ------
template<typename Real>
void Gather
( const Real* sbuf, int sc,
        Real* rbuf, int rc, int root, Comm comm );
template<typename Real>
void  Gather
( const Complex<Real>* sbuf, int sc,
        Complex<Real>* rbuf, int rc, int root, Comm comm );

// Non-blocking gather
// -------------------
template<typename Real>
void IGather
( const Real* sbuf, int sc,
        Real* rbuf, int rc, int root, Comm comm, Request& request );
template<typename Real>
void IGather
( const Complex<Real>* sbuf, int sc,
        Complex<Real>* rbuf, int rc, int root, Comm comm, Request& request );

// Gather with variable recv sizes
// -------------------------------
template<typename Real>
void Gather
( const Real* sbuf, int sc,
        Real* rbuf, const int* rcs, const int* rds, int root, Comm comm );
template<typename Real>
void Gather
( const Complex<Real>* sbuf, int sc,
        Complex<Real>* rbuf, const int* rcs, const int* rds, 
  int root, Comm comm );

// AllGather
// ---------
template<typename Real>
void AllGather
( const Real* sbuf, int sc,
        Real* rbuf, int rc, Comm comm );
template<typename Real>
void AllGather
( const Complex<Real>* sbuf, int sc,
        Complex<Real>* rbuf, int rc, Comm comm );

// AllGather with variable recv sizes
// ----------------------------------
template<typename Real>
void AllGather
( const Real* sbuf, int sc,
        Real* rbuf, const int* rcs, const int* rds, Comm comm );
template<typename Real>
void AllGather
( const Complex<Real>* sbuf, int sc,
        Complex<Real>* rbuf, const int* rcs, const int* rds, Comm comm );

// Scatter
// -------
template<typename Real>
void Scatter
( const Real* sbuf, int sc,
        Real* rbuf, int rc, int root, Comm comm );
template<typename Real>
void Scatter
( const Complex<Real>* sbuf, int sc,
        Complex<Real>* rbuf, int rc, int root, Comm comm );
// In-place option
template<typename Real>
void Scatter( Real* buf, int sc, int rc, int root, Comm comm );
template<typename Real>
void Scatter( Complex<Real>* buf, int sc, int rc, int root, Comm comm );

// AllToAll
// --------
template<typename Real>
void AllToAll
( const Real* sbuf, int sc,
        Real* rbuf, int rc, Comm comm );
template<typename Real>
void AllToAll
( const Complex<Real>* sbuf, int sc,
        Complex<Real>* rbuf, int rc, Comm comm );

// AllToAll with non-uniform send/recv sizes
// -----------------------------------------
template<typename Real>
void AllToAll
( const Real* sbuf, const int* scs, const int* sds,
        Real* rbuf, const int* rcs, const int* rds, Comm comm );
template<typename Real>
void AllToAll
( const Complex<Real>* sbuf, const int* scs, const int* sds,
        Complex<Real>* rbuf, const int* rcs, const int* rds, Comm comm );

template<typename T>
std::vector<T> AllToAll
( const std::vector<T>& sendBuf, 
  const std::vector<int>& sendCounts, 
  const std::vector<int>& sendDispls,
  mpi::Comm comm );

// Reduce
// ------
template<typename T>
void Reduce
( const T* sbuf, T* rbuf, int count, Op op, int root, Comm comm );
template<typename Real>
void Reduce
( const Complex<Real>* sbuf, Complex<Real>* rbuf, int count, Op op, 
  int root, Comm comm );

// Default to SUM
template<typename T>
void Reduce( const T* sbuf, T* rbuf, int count, int root, Comm comm );

// With a message-size of one
template<typename T>
T Reduce( T sb, Op op, int root, Comm comm );

// With a message-size of one and default to SUM
template<typename T>
T Reduce( T sb, int root, Comm comm );

// Single-buffer reduce
// --------------------
template<typename T>
void Reduce( T* buf, int count, Op op, int root, Comm comm );
template<typename Real>
void Reduce( Complex<Real>* buf, int count, Op op, int root, Comm comm );

// Default to SUM
template<typename T>
void Reduce( T* buf, int count, int root, Comm comm );

// AllReduce
// ---------
template<typename T>
void AllReduce( const T* sbuf, T* rbuf, int count, Op op, Comm comm );
template<typename Real>
void AllReduce
( const Complex<Real>* sbuf, Complex<Real>* rbuf, int count, Op op, Comm comm );

// Default to SUM
template<typename T>
void AllReduce( const T* sbuf, T* rbuf, int count, Comm comm );

// If the message-length is one
template<typename T>
T AllReduce( T sb, Op op, Comm comm );

// If the message-length is one (and default to SUM)
template<typename T>
T AllReduce( T sb, Comm comm );

// Single-buffer AllReduce
// -----------------------
template<typename T>
void AllReduce( T* buf, int count, Op op, Comm comm );
template<typename Real>
void AllReduce( Complex<Real>* buf, int count, Op op, Comm comm );

// Default to SUM
template<typename T>
void AllReduce( T* buf, int count, Comm comm );

// ReduceScatter
// -------------
template<typename Real>
void ReduceScatter
( Real* sbuf, Real* rbuf, int rc, Op op, Comm comm );
template<typename Real>
void ReduceScatter
( Complex<Real>* sbuf, Complex<Real>* rbuf, int rc, Op op, Comm comm );

// Default to SUM
template<typename T>
void ReduceScatter( T* sbuf, T* rbuf, int rc, Comm comm );

// Single-buffer ReduceScatter
// ---------------------------
template<typename Real>
void ReduceScatter( Real* buf, int rc, Op op, Comm comm );
template<typename Real>
void ReduceScatter( Complex<Real>* buf, int rc, Op op, Comm comm );

// Default to SUM
template<typename T>
void ReduceScatter( T* buf, int rc, Comm comm );

// Variable-length ReduceScatter
// -----------------------------
template<typename Real>
void ReduceScatter
( const Real* sbuf, Real* rbuf, const int* rcs, Op op, Comm comm );
template<typename Real>
void ReduceScatter
( const Complex<Real>* sbuf, Complex<Real>* rbuf, const int* rcs, Op op, 
  Comm comm );

// Default to SUM
template<typename T>
void ReduceScatter( const T* sbuf, T* rbuf, const int* rcs, Comm comm );

// Scan
// ----
template<typename T>
void Scan( const T* sbuf, T* rbuf, int count, Op op, Comm comm );
template<typename Real>
void Scan
( const Complex<Real>* sbuf, Complex<Real>* rbuf, int count, Op op, Comm comm );

// Default to SUM
template<typename T>
void Scan( const T* sbuf, T* rbuf, int count, Comm comm );

// With a message-size of one
template<typename T>
T Scan( T sb, Op op, Comm comm );

// With a message-size of one and default to SUM
template<typename T>
T Scan( T sb, Comm comm );

// Single-buffer scan
// ------------------
template<typename T>
void Scan( T* buf, int count, Op op, Comm comm );
template<typename Real>
void Scan( Complex<Real>* buf, int count, Op op, Comm comm );

// Default to SUM
template<typename T>
void Scan( T* buf, int count, Comm comm );

template<typename T>
void SparseAllToAll
( const std::vector<T>& sendBuffer,
  const std::vector<int>& sendCounts, 
  const std::vector<int>& sendOffs,
        std::vector<T>& recvBuffer,
  const std::vector<int>& recvCounts, 
  const std::vector<int>& recvOffs,
        Comm comm );

void VerifySendsAndRecvs
( const std::vector<int>& sendCounts,
  const std::vector<int>& recvCounts, Comm comm );

void CreateCustom();
void DestroyCustom();

template<typename T> Datatype TypeMap();
template<> Datatype TypeMap<byte>();
template<> Datatype TypeMap<int>();
template<> Datatype TypeMap<unsigned>();
template<> Datatype TypeMap<long int>();
template<> Datatype TypeMap<long unsigned>();
template<> Datatype TypeMap<long long int>();
template<> Datatype TypeMap<unsigned long long>();
template<> Datatype TypeMap<float>();
template<> Datatype TypeMap<double>();
template<> Datatype TypeMap<Complex<float>>();
template<> Datatype TypeMap<Complex<double>>();
#ifdef EL_HAVE_QUAD
template<> Datatype TypeMap<Complex<Quad>>();
#endif
template<> Datatype TypeMap<ValueInt<Int>>();
template<> Datatype TypeMap<ValueInt<float>>();
template<> Datatype TypeMap<ValueInt<double>>();
#ifdef EL_HAVE_QUAD
template<> Datatype TypeMap<ValueInt<Quad>>();
#endif
template<> Datatype TypeMap<ValueInt<Complex<float>>>();
template<> Datatype TypeMap<ValueInt<Complex<double>>>();
#ifdef EL_HAVE_QUAD
template<> Datatype TypeMap<ValueInt<Complex<Quad>>>();
#endif

template<> Datatype TypeMap<Entry<Int>>();
template<> Datatype TypeMap<Entry<float>>();
template<> Datatype TypeMap<Entry<double>>();
#ifdef EL_HAVE_QUAD
template<> Datatype TypeMap<Entry<Quad>>();
#endif
template<> Datatype TypeMap<Entry<Complex<float>>>();
template<> Datatype TypeMap<Entry<Complex<double>>>();
#ifdef EL_HAVE_QUAD
template<> Datatype TypeMap<Entry<Complex<Quad>>>();
#endif

template<typename Real> inline Op MaxOp() { return MAX; }
template<typename Real> inline Op MinOp() { return MIN; }
#ifdef EL_HAVE_QUAD
template<> Op MaxOp<Quad>();
template<> Op MinOp<Quad>();
#endif

template<typename T> inline Op SumOp() { return SUM; }
#ifdef EL_HAVE_QUAD
template<> Op SumOp<Quad>();
template<> Op SumOp<Complex<Quad>>();
#endif

template<typename Real> Op MaxLocOp();
template<> Op MaxLocOp<Int>();
template<> Op MaxLocOp<float>();
template<> Op MaxLocOp<double>();
#ifdef EL_HAVE_QUAD
template<> Op MaxLocOp<Quad>();
#endif

template<typename Real> Op MaxLocPairOp();
template<> Op MaxLocPairOp<Int>();
template<> Op MaxLocPairOp<float>();
template<> Op MaxLocPairOp<double>();
#ifdef EL_HAVE_QUAD
template<> Op MaxLocPairOp<Quad>();
#endif

template<typename Real> Op MinLocOp();
template<> Op MinLocOp<Int>();
template<> Op MinLocOp<float>();
template<> Op MinLocOp<double>();
#ifdef EL_HAVE_QUAD
template<> Op MinLocOp<Quad>();
#endif

template<typename Real> Op MinLocPairOp();
template<> Op MinLocPairOp<Int>();
template<> Op MinLocPairOp<float>();
template<> Op MinLocPairOp<double>();
#ifdef EL_HAVE_QUAD
template<> Op MinLocPairOp<Quad>();
#endif

} // mpi
} // elem

#endif // ifndef EL_IMPORTS_MPI_HPP
