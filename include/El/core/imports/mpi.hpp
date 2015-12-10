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

using std::function;
using std::vector;

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
    Comm( MPI_Comm mpiComm=MPI_COMM_WORLD ) EL_NO_EXCEPT : comm(mpiComm) { }

    inline int Rank() const EL_NO_RELEASE_EXCEPT;
    inline int Size() const EL_NO_RELEASE_EXCEPT;
};
inline bool operator==( const Comm& a, const Comm& b ) EL_NO_EXCEPT
{ return a.comm == b.comm; }
inline bool operator!=( const Comm& a, const Comm& b ) EL_NO_EXCEPT
{ return a.comm != b.comm; }
// Hopefully, despite the fact that MPI_Comm is opaque, the following will
// reliably hold (otherwise it must be extended). Typically, MPI_Comm is
// either 'int' or 'void*'.
inline bool operator<( const Comm& a, const Comm& b ) EL_NO_EXCEPT
{ return a.comm < b.comm; }

struct Group
{
    MPI_Group group;
    Group( MPI_Group mpiGroup=MPI_GROUP_NULL ) EL_NO_EXCEPT 
    : group(mpiGroup) { }

    inline int Rank() const EL_NO_RELEASE_EXCEPT;
    inline int Size() const EL_NO_RELEASE_EXCEPT;
};
inline bool operator==( const Group& a, const Group& b ) EL_NO_EXCEPT
{ return a.group == b.group; }
inline bool operator!=( const Group& a, const Group& b ) EL_NO_EXCEPT
{ return a.group != b.group; }

struct Op
{
    MPI_Op op;
    Op( MPI_Op mpiOp=MPI_SUM ) EL_NO_EXCEPT : op(mpiOp) { }
};
inline bool operator==( const Op& a, const Op& b ) EL_NO_EXCEPT
{ return a.op == b.op; }
inline bool operator!=( const Op& a, const Op& b ) EL_NO_EXCEPT
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

template<typename T> Op UserOp();
template<typename T> Op UserCommOp();

template<typename Real> inline Op MaxOp() EL_NO_EXCEPT { return MAX; }
template<typename Real> inline Op MinOp() EL_NO_EXCEPT { return MIN; }
#ifdef EL_HAVE_QUAD
template<> Op MaxOp<Quad>() EL_NO_EXCEPT;
template<> Op MinOp<Quad>() EL_NO_EXCEPT;
#endif

template<typename T> inline Op SumOp() EL_NO_EXCEPT { return SUM; }
#ifdef EL_HAVE_QUAD
template<> Op SumOp<Quad>() EL_NO_EXCEPT;
template<> Op SumOp<Complex<Quad>>() EL_NO_EXCEPT;
#endif

template<typename Real> Op MaxLocOp() EL_NO_EXCEPT;
template<> Op MaxLocOp<Int>() EL_NO_EXCEPT;
template<> Op MaxLocOp<float>() EL_NO_EXCEPT;
template<> Op MaxLocOp<double>() EL_NO_EXCEPT;
#ifdef EL_HAVE_QUAD
template<> Op MaxLocOp<Quad>() EL_NO_EXCEPT;
#endif

template<typename Real> Op MaxLocPairOp() EL_NO_EXCEPT;
template<> Op MaxLocPairOp<Int>() EL_NO_EXCEPT;
template<> Op MaxLocPairOp<float>() EL_NO_EXCEPT;
template<> Op MaxLocPairOp<double>() EL_NO_EXCEPT;
#ifdef EL_HAVE_QUAD
template<> Op MaxLocPairOp<Quad>() EL_NO_EXCEPT;
#endif

template<typename Real> Op MinLocOp() EL_NO_EXCEPT;
template<> Op MinLocOp<Int>() EL_NO_EXCEPT;
template<> Op MinLocOp<float>() EL_NO_EXCEPT;
template<> Op MinLocOp<double>() EL_NO_EXCEPT;
#ifdef EL_HAVE_QUAD
template<> Op MinLocOp<Quad>() EL_NO_EXCEPT;
#endif

template<typename Real> Op MinLocPairOp() EL_NO_EXCEPT;
template<> Op MinLocPairOp<Int>() EL_NO_EXCEPT;
template<> Op MinLocPairOp<float>() EL_NO_EXCEPT;
template<> Op MinLocPairOp<double>() EL_NO_EXCEPT;
#ifdef EL_HAVE_QUAD
template<> Op MinLocPairOp<Quad>() EL_NO_EXCEPT;
#endif

// Added constant(s)
const int MIN_COLL_MSG = 1; // minimum message size for collectives
inline int Pad( int count ) EL_NO_EXCEPT
{ return std::max(count,MIN_COLL_MSG); }

bool CommSameSizeAsInteger() EL_NO_EXCEPT;
bool GroupSameSizeAsInteger() EL_NO_EXCEPT;

// Environment routines
void Initialize( int& argc, char**& argv ) EL_NO_EXCEPT;
int InitializeThread( int& argc, char**& argv, int required ) EL_NO_EXCEPT;
void Finalize() EL_NO_EXCEPT;
bool Initialized() EL_NO_EXCEPT;
bool Finalized() EL_NO_EXCEPT;
int QueryThread() EL_NO_EXCEPT;
void Abort( Comm comm, int errCode ) EL_NO_EXCEPT;
double Time() EL_NO_EXCEPT;
void Create( UserFunction* func, bool commutes, Op& op ) EL_NO_RELEASE_EXCEPT;
void Free( Op& op ) EL_NO_RELEASE_EXCEPT;
void Free( Datatype& type ) EL_NO_RELEASE_EXCEPT;

// Communicator manipulation
int Rank( Comm comm=COMM_WORLD ) EL_NO_RELEASE_EXCEPT;
int Size( Comm comm=COMM_WORLD ) EL_NO_RELEASE_EXCEPT;
void Create
( Comm parentComm, Group subsetGroup, Comm& subsetComm ) EL_NO_RELEASE_EXCEPT;
void Dup( Comm original, Comm& duplicate ) EL_NO_RELEASE_EXCEPT;
void Split( Comm comm, int color, int key, Comm& newComm ) EL_NO_RELEASE_EXCEPT;
void Free( Comm& comm ) EL_NO_RELEASE_EXCEPT;
bool Congruent( Comm comm1, Comm comm2 ) EL_NO_RELEASE_EXCEPT;
void ErrorHandlerSet
( Comm comm, ErrorHandler errorHandler ) EL_NO_RELEASE_EXCEPT;

// Cartesian communicator routines
void CartCreate
( Comm comm, int numDims, const int* dimensions, const int* periods, 
  bool reorder, Comm& cartComm ) EL_NO_RELEASE_EXCEPT;
void CartSub
( Comm comm, const int* remainingDims, Comm& subComm ) EL_NO_RELEASE_EXCEPT;

// Group manipulation
int Rank( Group group ) EL_NO_RELEASE_EXCEPT;
int Size( Group group ) EL_NO_RELEASE_EXCEPT;
void CommGroup( Comm comm, Group& group ) EL_NO_RELEASE_EXCEPT;
void Dup( Group group, Group& newGroup ) EL_NO_RELEASE_EXCEPT;
void Union( Group groupA, Group groupB, Group& newGroup ) EL_NO_RELEASE_EXCEPT;
void Incl
( Group group, int n, const int* ranks, Group& subGroup ) EL_NO_RELEASE_EXCEPT;
void Excl
( Group group, int n, const int* ranks, Group& subGroup ) EL_NO_RELEASE_EXCEPT;
void Difference
( Group parent, Group subset, Group& complement ) EL_NO_RELEASE_EXCEPT;
void Free( Group& group ) EL_NO_RELEASE_EXCEPT;
bool Congruent( Group group1, Group group2 ) EL_NO_RELEASE_EXCEPT;
int Translate
( Group origGroup, int origRank, Group newGroup ) EL_NO_RELEASE_EXCEPT;
int Translate
( Comm  origComm,  int origRank, Group newGroup ) EL_NO_RELEASE_EXCEPT;
int Translate
( Group origGroup, int origRank, Comm  newComm  ) EL_NO_RELEASE_EXCEPT;
int Translate
( Comm  origComm,  int origRank, Comm  newComm  ) EL_NO_RELEASE_EXCEPT;
void Translate
( Group origGroup, int size, const int* origRanks, 
  Group newGroup,                  int* newRanks ) EL_NO_RELEASE_EXCEPT;
void Translate
( Comm origComm,  int size, const int* origRanks, 
  Group newGroup,                 int* newRanks ) EL_NO_RELEASE_EXCEPT;
void Translate
( Group origGroup, int size, const int* origRanks, 
  Comm newComm,                    int* newRanks ) EL_NO_RELEASE_EXCEPT;
void Translate
( Comm origComm, int size, const int* origRanks, 
  Comm newComm,                  int* newRanks ) EL_NO_RELEASE_EXCEPT;

// Utilities
void Barrier( Comm comm=COMM_WORLD ) EL_NO_RELEASE_EXCEPT;
void Wait( Request& request ) EL_NO_RELEASE_EXCEPT;
void Wait( Request& request, Status& status ) EL_NO_RELEASE_EXCEPT;
void WaitAll( int numRequests, Request* requests ) EL_NO_RELEASE_EXCEPT;
void WaitAll
( int numRequests, Request* requests, Status* statuses ) EL_NO_RELEASE_EXCEPT;
bool Test( Request& request ) EL_NO_RELEASE_EXCEPT;
bool IProbe
( int source, int tag, Comm comm, Status& status ) EL_NO_RELEASE_EXCEPT;

template<typename T>
int GetCount( Status& status ) EL_NO_RELEASE_EXCEPT;

// NOTE: This is instantiated for the standard datatypes
template<typename T>
void SetUserReduceFunc
( function<T(const T&,const T&)> func, bool commutative=true );

// Point-to-point communication
// ============================

// Send
// ----
template<typename Real>
void TaggedSend
( const Real* buf, int count, int to, int tag, Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT;
template<typename Real>
void TaggedSend
( const Complex<Real>* buf, int count, int to, int tag, Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT;

// If the tag is irrelevant
template<typename T>
void Send( const T* buf, int count, int to, Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT;

// If the send-count is one
template<typename T>
void TaggedSend( T b, int to, int tag, Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT;

// If the send-count is one and the tag is irrelevant
template<typename T>
void Send( T b, int to, Comm comm=COMM_WORLD ) EL_NO_RELEASE_EXCEPT;

// Non-blocking send
// -----------------
template<typename Real>
void TaggedISend
( const Real* buf, int count, int to, int tag, Comm comm, Request& request )
EL_NO_RELEASE_EXCEPT;
template<typename Real>
void TaggedISend
( const Complex<Real>* buf, int count, int to, int tag, Comm comm, 
  Request& request ) EL_NO_RELEASE_EXCEPT;

// If the tag is irrelevant
template<typename T>
void ISend( const T* buf, int count, int to, Comm comm, Request& request )
EL_NO_RELEASE_EXCEPT;

// If the send count is one
template<typename T>
void TaggedISend( T b, int to, int tag, Comm comm, Request& request )
EL_NO_RELEASE_EXCEPT;

// If the send count is one and the tag is irrelevant
template<typename T>
void ISend( T b, int to, Comm comm, Request& request ) EL_NO_RELEASE_EXCEPT;

// Non-blocking synchronous Send
// -----------------------------
template<typename Real>
void TaggedISSend
( const Real* buf, int count, int to, int tag, Comm comm, Request& request )
EL_NO_RELEASE_EXCEPT;
template<typename Real>
void TaggedISSend
( const Complex<Real>* buf, int count, int to, int tag, Comm comm, 
  Request& request ) EL_NO_RELEASE_EXCEPT;

// If the tag is irrelevant
template<typename T>
void ISSend( const T* buf, int count, int to, Comm comm, Request& request )
EL_NO_RELEASE_EXCEPT;

// If the send count is one
template<typename T>
void TaggedISSend( T b, int to, int tag, Comm comm, Request& request )
EL_NO_RELEASE_EXCEPT;

// If the send count is one and the tag is irrelevant
template<typename T>
void ISSend( T b, int to, Comm comm, Request& request )
EL_NO_RELEASE_EXCEPT;

// Recv
// ----
template<typename Real>
void TaggedRecv
( Real* buf, int count, int from, int tag, Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT;
template<typename Real>
void TaggedRecv
( Complex<Real>* buf, int count, int from, int tag, Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT;

// If the tag is irrelevant
template<typename T>
void Recv( T* buf, int count, int from, Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT;

// If the recv count is one
template<typename T>
T TaggedRecv( int from, int tag, Comm comm=COMM_WORLD ) EL_NO_RELEASE_EXCEPT;

// If the recv count is one and the tag is irrelevant
template<typename T>
T Recv( int from, Comm comm=COMM_WORLD ) EL_NO_RELEASE_EXCEPT;

// Non-blocking recv
// -----------------
template<typename Real>
void TaggedIRecv
( Real* buf, int count, int from, int tag, Comm comm, 
  Request& request ) EL_NO_RELEASE_EXCEPT;
template<typename Real>
void TaggedIRecv
( Complex<Real>* buf, int count, int from, int tag, Comm comm, 
  Request& request ) EL_NO_RELEASE_EXCEPT;

// If the tag is irrelevant
template<typename T>
void IRecv( T* buf, int count, int from, Comm comm, Request& request )
EL_NO_RELEASE_EXCEPT;

// If the recv count is one
template<typename T>
T TaggedIRecv( int from, int tag, Comm comm, Request& request )
EL_NO_RELEASE_EXCEPT;

// If the recv count is one and the tag is irrelevant
template<typename T>
T IRecv( int from, Comm comm, Request& request ) EL_NO_RELEASE_EXCEPT;

// SendRecv
// --------
template<typename Real>
void TaggedSendRecv
( const Real* sbuf, int sc, int to,   int stag,
        Real* rbuf, int rc, int from, int rtag, Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT;
template<typename Real>
void TaggedSendRecv
( const Complex<Real>* sbuf, int sc, int to,   int stag,
        Complex<Real>* rbuf, int rc, int from, int rtag, Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT;

// If the tags are irrelevant
template<typename T>
void SendRecv
( const T* sbuf, int sc, int to,
        T* rbuf, int rc, int from, Comm comm=COMM_WORLD ) EL_NO_RELEASE_EXCEPT;

// If the send and recv counts are one
template<typename T>
T TaggedSendRecv
( T sb, int to, int stag, int from, int rtag, Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT;

// If the send and recv counts are one and the tags don't matter
template<typename T>
T SendRecv( T sb, int to, int from, Comm comm=COMM_WORLD ) EL_NO_RELEASE_EXCEPT;

// Single-buffer SendRecv
// ----------------------
template<typename Real>
void TaggedSendRecv
( Real* buf, int count, int to, int stag, int from, int rtag,
  Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT;
template<typename Real>
void TaggedSendRecv
( Complex<Real>* buf, int count, int to, int stag, int from, int rtag, 
  Comm comm=COMM_WORLD ) EL_NO_RELEASE_EXCEPT;

// If the tags don't matter
template<typename T>
void SendRecv( T* buf, int count, int to, int from, Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT;

// Collective communication
// ========================

// Broadcast
// ---------
template<typename Real>
void Broadcast( Real* buf, int count, int root, Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT;
template<typename Real>
void Broadcast( Complex<Real>* buf, int count, int root, Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT;

// If the message length is one
template<typename T>
void Broadcast( T& b, int root, Comm comm=COMM_WORLD ) EL_NO_RELEASE_EXCEPT;

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
// Even though EL_AVOID_COMPLEX_MPI being defined implies that an std::vector
// copy of the input data will be created, and the memory allocation can clearly
// fail and throw an exception, said exception is not necessarily thrown on
// Linux platforms due to the "optimistic" allocation policy. Therefore we will
// go ahead and allow std::terminate to be called should such an std::bad_alloc
// exception occur in a Release build
template<typename Real>
void Gather
( const Real* sbuf, int sc,
        Real* rbuf, int rc, int root,
  Comm comm=COMM_WORLD ) EL_NO_RELEASE_EXCEPT;
template<typename Real>
void  Gather
( const Complex<Real>* sbuf, int sc,
        Complex<Real>* rbuf, int rc, int root,
  Comm comm=COMM_WORLD ) EL_NO_RELEASE_EXCEPT;

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
        Real* rbuf, const int* rcs, const int* rds, int root,
  Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT;
template<typename Real>
void Gather
( const Complex<Real>* sbuf, int sc,
        Complex<Real>* rbuf, const int* rcs, const int* rds, 
  int root, Comm comm=COMM_WORLD ) EL_NO_RELEASE_EXCEPT;

// AllGather
// ---------
// NOTE: See the corresponding note for Gather on std::bad_alloc exceptions
template<typename Real>
void AllGather
( const Real* sbuf, int sc,
        Real* rbuf, int rc, Comm comm=COMM_WORLD ) EL_NO_RELEASE_EXCEPT;
template<typename Real>
void AllGather
( const Complex<Real>* sbuf, int sc,
        Complex<Real>* rbuf, int rc, Comm comm=COMM_WORLD )
  EL_NO_RELEASE_EXCEPT;

// AllGather with variable recv sizes
// ----------------------------------
template<typename Real>
void AllGather
( const Real* sbuf, int sc,
        Real* rbuf, const int* rcs, const int* rds, Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT;
template<typename Real>
void AllGather
( const Complex<Real>* sbuf, int sc,
        Complex<Real>* rbuf, const int* rcs, const int* rds,
  Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT;

// Scatter
// -------
template<typename Real>
void Scatter
( const Real* sbuf, int sc,
        Real* rbuf, int rc, int root, Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT;
template<typename Real>
void Scatter
( const Complex<Real>* sbuf, int sc,
        Complex<Real>* rbuf, int rc, int root, Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT;
// In-place option
template<typename Real>
void Scatter
( Real* buf, int sc, int rc, int root, Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT;
template<typename Real>
void Scatter
( Complex<Real>* buf, int sc, int rc, int root, Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT;

// AllToAll
// --------
// NOTE: See the corresponding note on std::bad_alloc for Gather
template<typename Real>
void AllToAll
( const Real* sbuf, int sc,
        Real* rbuf, int rc, Comm comm=COMM_WORLD ) EL_NO_RELEASE_EXCEPT;
template<typename Real>
void AllToAll
( const Complex<Real>* sbuf, int sc,
        Complex<Real>* rbuf, int rc, Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT;

// AllToAll with non-uniform send/recv sizes
// -----------------------------------------
template<typename Real>
void AllToAll
( const Real* sbuf, const int* scs, const int* sds,
        Real* rbuf, const int* rcs, const int* rds, Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT;
template<typename Real>
void AllToAll
( const Complex<Real>* sbuf, const int* scs, const int* sds,
        Complex<Real>* rbuf, const int* rcs, const int* rds,
  Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT;

template<typename T>
vector<T> AllToAll
( const vector<T>& sendBuf, 
  const vector<int>& sendCounts, 
  const vector<int>& sendDispls,
  Comm comm=COMM_WORLD ) EL_NO_RELEASE_EXCEPT;

// Reduce
// ------
template<typename T>
void Reduce
( const T* sbuf, T* rbuf, int count, Op op, int root, Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT;

template<typename Real>
void Reduce
( const Complex<Real>* sbuf, Complex<Real>* rbuf, int count, Op op, 
  int root, Comm comm=COMM_WORLD ) EL_NO_RELEASE_EXCEPT;

template<typename T,class OpClass,typename=DisableIf<IsData<OpClass>>>
inline void Reduce
( const T* sb, T* rb, int count, OpClass op, bool commutative,
  int root, Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT
{
    SetUserReduceFunc( function<T(const T&,const T&)>(op), commutative );
    if( commutative )
        Reduce( sb, rb, count, UserCommOp<T>(), root, comm ); 
    else
        Reduce( sb, rb, count, UserOp<T>(), root, comm ); 
}

// Default to SUM
template<typename T>
void Reduce( const T* sbuf, T* rbuf, int count, int root, Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT;

// With a message-size of one
template<typename T>
T Reduce( T sb, Op op, int root, Comm comm=COMM_WORLD ) EL_NO_RELEASE_EXCEPT;

template<typename T,class OpClass,typename=DisableIf<IsData<OpClass>>>
inline T Reduce
( T sb, OpClass op, bool commutative, int root, Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT
{
    SetUserReduceFunc( function<T(const T&,const T&)>(op), commutative );
    if( commutative )
        return Reduce( sb, UserCommOp<T>(), root, comm ); 
    else
        return Reduce( sb, UserOp<T>(), root, comm ); 
}

// With a message-size of one and default to SUM
template<typename T>
T Reduce( T sb, int root, Comm comm=COMM_WORLD ) EL_NO_RELEASE_EXCEPT;

// Single-buffer reduce
// --------------------
template<typename T>
void Reduce( T* buf, int count, Op op, int root, Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT;

template<typename Real>
void Reduce
( Complex<Real>* buf, int count, Op op, int root, Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT;

template<typename T,class OpClass,typename=DisableIf<IsData<OpClass>>>
inline void Reduce
( T* buf, int count, OpClass op, bool commutative, int root,
  Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT
{
    SetUserReduceFunc( function<T(const T&,const T&)>(op), commutative );
    if( commutative )
        Reduce( buf, count, UserCommOp<T>(), root, comm ); 
    else
        Reduce( buf, count, UserOp<T>(), root, comm ); 
}

// Default to SUM
template<typename T>
void Reduce( T* buf, int count, int root, Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT;

// AllReduce
// ---------
template<typename T>
void AllReduce( const T* sbuf, T* rbuf, int count, Op op, Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT;

template<typename Real>
void AllReduce
( const Complex<Real>* sbuf, Complex<Real>* rbuf, int count, Op op,
  Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT;

template<typename T,class OpClass,typename=DisableIf<IsData<OpClass>>>
inline void AllReduce
( const T* sb, T* rb, int count, OpClass op, bool commutative,
  Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT
{
    SetUserReduceFunc( function<T(const T&,const T&)>(op), commutative );
    if( commutative )
        AllReduce( sb, rb, count, UserCommOp<T>(), comm ); 
    else
        AllReduce( sb, rb, count, UserOp<T>(), comm ); 
}

// Default to SUM
template<typename T>
void AllReduce( const T* sbuf, T* rbuf, int count, Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT;

// If the message-length is one
template<typename T>
T AllReduce( T sb, Op op, Comm comm=COMM_WORLD ) EL_NO_RELEASE_EXCEPT;

template<typename T,class OpClass,typename=DisableIf<IsData<OpClass>>>
inline T AllReduce
( T sb, OpClass op, bool commutative, Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT
{
    SetUserReduceFunc( function<T(const T&,const T&)>(op), commutative );
    if( commutative )
        return AllReduce( sb, UserCommOp<T>(), comm ); 
    else
        return AllReduce( sb, UserOp<T>(), comm ); 
}

// If the message-length is one (and default to SUM)
template<typename T>
T AllReduce( T sb, Comm comm=COMM_WORLD ) EL_NO_RELEASE_EXCEPT;

// Single-buffer AllReduce
// -----------------------
template<typename T>
void AllReduce( T* buf, int count, Op op, Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT;

template<typename Real>
void AllReduce( Complex<Real>* buf, int count, Op op, Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT;

template<typename T,class OpClass,typename=DisableIf<IsData<OpClass>>>
inline void AllReduce
( T* buf, int count, OpClass op, bool commutative, Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT
{
    SetUserReduceFunc( function<T(const T&,const T&)>(op), commutative );
    if( commutative )
        AllReduce( buf, count, UserCommOp<T>(), comm ); 
    else
        AllReduce( buf, count, UserOp<T>(), comm ); 
}

// Default to SUM
template<typename T>
void AllReduce( T* buf, int count, Comm comm=COMM_WORLD ) EL_NO_RELEASE_EXCEPT;

// ReduceScatter
// -------------
template<typename Real>
void ReduceScatter
( Real* sbuf, Real* rbuf, int rc, Op op, Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT;

template<typename Real>
void ReduceScatter
( Complex<Real>* sbuf, Complex<Real>* rbuf, int rc, Op op,
  Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT;

template<typename T,class OpClass,typename=DisableIf<IsData<OpClass>>>
inline void ReduceScatter
( const T* sb, T* rb, int count, OpClass op, bool commutative,
  Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT
{
    SetUserReduceFunc( function<T(const T&,const T&)>(op), commutative );
    if( commutative )
        ReduceScatter( sb, rb, count, UserCommOp<T>(), comm ); 
    else
        ReduceScatter( sb, rb, count, UserOp<T>(), comm ); 
}

// Default to SUM
template<typename T>
void ReduceScatter( T* sbuf, T* rbuf, int rc, Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT;

// Single-buffer ReduceScatter
// ---------------------------
template<typename Real>
void ReduceScatter( Real* buf, int rc, Op op, Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT;

template<typename Real>
void ReduceScatter( Complex<Real>* buf, int rc, Op op, Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT;

template<typename T,class OpClass,typename=DisableIf<IsData<OpClass>>>
inline void ReduceScatter
( T* buf, int count, OpClass op, bool commutative, Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT
{
    SetUserReduceFunc( function<T(const T&,const T&)>(op), commutative );
    if( commutative )
        ReduceScatter( buf, count, UserCommOp<T>(), comm ); 
    else
        ReduceScatter( buf, count, UserOp<T>(), comm ); 
}

// Default to SUM
template<typename T>
void ReduceScatter( T* buf, int rc, Comm comm=COMM_WORLD ) EL_NO_RELEASE_EXCEPT;

// Variable-length ReduceScatter
// -----------------------------
template<typename Real>
void ReduceScatter
( const Real* sbuf, Real* rbuf, const int* rcs, Op op, Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT;

template<typename Real>
void ReduceScatter
( const Complex<Real>* sbuf, Complex<Real>* rbuf, const int* rcs, Op op, 
  Comm comm=COMM_WORLD ) EL_NO_RELEASE_EXCEPT;

template<typename T,class OpClass,typename=DisableIf<IsData<OpClass>>>
inline void ReduceScatter
( const T* sb, T* rb, const int* rcs, OpClass op, bool commutative,
  Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT
{
    SetUserReduceFunc( function<T(const T&,const T&)>(op), commutative );
    if( commutative )
        ReduceScatter( sb, rb, rcs, UserCommOp<T>(), comm ); 
    else
        ReduceScatter( sb, rb, rcs, UserOp<T>(), comm ); 
}

// Default to SUM
template<typename T>
void ReduceScatter
( const T* sbuf, T* rbuf, const int* rcs, Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT;

// Scan
// ----
template<typename T>
void Scan( const T* sbuf, T* rbuf, int count, Op op, Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT;

template<typename Real>
void Scan
( const Complex<Real>* sbuf, Complex<Real>* rbuf, int count, Op op,
  Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT;

template<typename T,class OpClass,typename=DisableIf<IsData<OpClass>>>
inline void Scan
( const T* sb, T* rb, int count, OpClass op, bool commutative,
  int root, Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT
{
    SetUserReduceFunc( function<T(const T&,const T&)>(op), commutative );
    if( commutative )
        Scan( sb, rb, count, UserCommOp<T>(), root, comm ); 
    else
        Scan( sb, rb, count, UserOp<T>(), root, comm ); 
}

// Default to SUM
template<typename T>
void Scan
( const T* sbuf, T* rbuf, int count, Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT;

// With a message-size of one
template<typename T>
T Scan( T sb, Op op, Comm comm=COMM_WORLD ) EL_NO_RELEASE_EXCEPT;

template<typename T,class OpClass,typename=DisableIf<IsData<OpClass>>>
inline T Scan
( T sb, OpClass op, bool commutative, int root, Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT
{
    SetUserReduceFunc( function<T(const T&,const T&)>(op), commutative );
    if( commutative )
        return Scan( sb, UserCommOp<T>(), root, comm ); 
    else
        return Scan( sb, UserOp<T>(), root, comm ); 
}

// With a message-size of one and default to SUM
template<typename T>
T Scan( T sb, Comm comm=COMM_WORLD ) EL_NO_RELEASE_EXCEPT;

// Single-buffer scan
// ------------------
template<typename T>
void Scan( T* buf, int count, Op op, Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT;

template<typename Real>
void Scan( Complex<Real>* buf, int count, Op op, Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT;

template<typename T,class OpClass,typename=DisableIf<IsData<OpClass>>>
inline void Scan
( T* buf, int count, OpClass op, bool commutative, int root,
  Comm comm=COMM_WORLD )
EL_NO_RELEASE_EXCEPT
{
    SetUserReduceFunc( function<T(const T&,const T&)>(op), commutative );
    if( commutative )
        Scan( buf, count, UserCommOp<T>(), root, comm ); 
    else
        Scan( buf, count, UserOp<T>(), root, comm ); 
}

// Default to SUM
template<typename T>
void Scan( T* buf, int count, Comm comm=COMM_WORLD ) EL_NO_RELEASE_EXCEPT;

template<typename T>
void SparseAllToAll
( const vector<T>& sendBuffer,
  const vector<int>& sendCounts, 
  const vector<int>& sendOffs,
        vector<T>& recvBuffer,
  const vector<int>& recvCounts, 
  const vector<int>& recvOffs,
        Comm comm=COMM_WORLD ) EL_NO_RELEASE_EXCEPT;

void VerifySendsAndRecvs
( const vector<int>& sendCounts,
  const vector<int>& recvCounts, Comm comm=COMM_WORLD );

void CreateCustom() EL_NO_RELEASE_EXCEPT;
void DestroyCustom() EL_NO_RELEASE_EXCEPT;

template<typename T> Datatype TypeMap() EL_NO_EXCEPT;
template<> Datatype TypeMap<byte>() EL_NO_EXCEPT;
template<> Datatype TypeMap<short>() EL_NO_EXCEPT;
template<> Datatype TypeMap<int>() EL_NO_EXCEPT;
template<> Datatype TypeMap<unsigned>() EL_NO_EXCEPT;
template<> Datatype TypeMap<long int>() EL_NO_EXCEPT;
template<> Datatype TypeMap<long unsigned>() EL_NO_EXCEPT;
#ifdef EL_HAVE_MPI_LONG_LONG
template<> Datatype TypeMap<long long int>() EL_NO_EXCEPT;
template<> Datatype TypeMap<unsigned long long>() EL_NO_EXCEPT;
#endif
template<> Datatype TypeMap<float>() EL_NO_EXCEPT;
template<> Datatype TypeMap<double>() EL_NO_EXCEPT;
template<> Datatype TypeMap<Complex<float>>() EL_NO_EXCEPT;
template<> Datatype TypeMap<Complex<double>>() EL_NO_EXCEPT;
#ifdef EL_HAVE_QUAD
template<> Datatype TypeMap<Complex<Quad>>() EL_NO_EXCEPT;
#endif
template<> Datatype TypeMap<ValueInt<Int>>() EL_NO_EXCEPT;
template<> Datatype TypeMap<ValueInt<float>>() EL_NO_EXCEPT;
template<> Datatype TypeMap<ValueInt<double>>() EL_NO_EXCEPT;
#ifdef EL_HAVE_QUAD
template<> Datatype TypeMap<ValueInt<Quad>>() EL_NO_EXCEPT;
#endif
template<> Datatype TypeMap<ValueInt<Complex<float>>>() EL_NO_EXCEPT;
template<> Datatype TypeMap<ValueInt<Complex<double>>>() EL_NO_EXCEPT;
#ifdef EL_HAVE_QUAD
template<> Datatype TypeMap<ValueInt<Complex<Quad>>>() EL_NO_EXCEPT;
#endif

template<> Datatype TypeMap<Entry<Int>>() EL_NO_EXCEPT;
template<> Datatype TypeMap<Entry<float>>() EL_NO_EXCEPT;
template<> Datatype TypeMap<Entry<double>>() EL_NO_EXCEPT;
#ifdef EL_HAVE_QUAD
template<> Datatype TypeMap<Entry<Quad>>() EL_NO_EXCEPT;
#endif
template<> Datatype TypeMap<Entry<Complex<float>>>() EL_NO_EXCEPT;
template<> Datatype TypeMap<Entry<Complex<double>>>() EL_NO_EXCEPT;
#ifdef EL_HAVE_QUAD
template<> Datatype TypeMap<Entry<Complex<Quad>>>() EL_NO_EXCEPT;
#endif

// Convenience functions which might not be very useful
int Comm::Rank() const EL_NO_RELEASE_EXCEPT { return mpi::Rank(*this); }
int Comm::Size() const EL_NO_RELEASE_EXCEPT { return mpi::Size(*this); }
int Group::Rank() const EL_NO_RELEASE_EXCEPT { return mpi::Rank(*this); }
int Group::Size() const EL_NO_RELEASE_EXCEPT { return mpi::Size(*this); }

} // mpi
} // elem

#endif // ifndef EL_IMPORTS_MPI_HPP
