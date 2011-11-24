/*
   Copyright (c) 2009-2011, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/
#ifndef ELEMENTAL_IMPORTS_MPI_HPP
#define ELEMENTAL_IMPORTS_MPI_HPP 1

namespace elemental {
namespace mpi {

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
const Comm COMM_WORLD = MPI_COMM_WORLD;
const ErrorHandler ERRORS_RETURN = MPI_ERRORS_RETURN;
const ErrorHandler ERRORS_ARE_FATAL = MPI_ERRORS_ARE_FATAL;
const Group GROUP_EMPTY = MPI_GROUP_EMPTY;
const Request REQUEST_NULL = MPI_REQUEST_NULL;
const Op MAX = MPI_MAX;
const Op PROD = MPI_PROD;
const Op SUM = MPI_SUM;

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
bool Test( Request& request );
bool IProbe( int source, int tag, Comm comm, Status& status );
template<typename T>
int GetCount( Status& status );
template<> int GetCount<byte>( Status& status );
template<> int GetCount<int>( Status& status );
template<> int GetCount<float>( Status& status );
template<> int GetCount<double>( Status& status );
#ifndef WITHOUT_COMPLEX
template<> int GetCount<scomplex>( Status& status );
template<> int GetCount<dcomplex>( Status& status );
#endif

// Point-to-point communication

void Send( const byte* buf, int count, int to, int tag, Comm comm );
void Send( const int* buf, int count, int to, int tag, Comm comm );
void Send( const float* buf, int count, int to, int tag, Comm comm );
void Send( const double* buf, int count, int to, int tag, Comm comm );
#ifndef WITHOUT_COMPLEX
void Send( const scomplex* buf, int count, int to, int tag, Comm comm );
void Send( const dcomplex* buf, int count, int to, int tag, Comm comm );
#endif

void ISend
( const byte* buf, int count, int to, int tag, Comm comm, Request& request );
void ISend
( const int* buf, int count, int to, int tag, Comm comm, Request& request );
void ISend
( const float* buf, int count, int to, int tag, Comm comm, Request& request );
void ISend
( const double* buf, int count, int to, int tag, Comm comm, Request& request );
#ifndef WITHOUT_COMPLEX
void ISend
( const scomplex* buf, int count, int to, int tag, Comm comm, 
  Request& request );
void ISend
( const dcomplex* buf, int count, int to, int tag, Comm comm, 
  Request& request );
#endif

void ISSend
( const byte* buf, int count, int to, int tag, Comm comm, Request& request );
void ISSend
( const int* buf, int count, int to, int tag, Comm comm, Request& request );
void ISSend
( const float* buf, int count, int to, int tag, Comm comm, Request& request );
void ISSend
( const double* buf, int count, int to, int tag, Comm comm, Request& request );
#ifndef WITHOUT_COMPLEX
void ISSend
( const scomplex* buf, int count, int to, int tag, Comm comm, 
  Request& request );
void ISSend
( const dcomplex* buf, int count, int to, int tag, Comm comm, 
 Request& request );
#endif
 
void Recv( byte* buf, int count, int from, int tag, Comm comm );
void Recv( int* buf, int count, int from, int tag, Comm comm );
void Recv( float* buf, int count, int from, int tag, Comm comm );
void Recv( double* buf, int count, int from, int tag, Comm comm );
#ifndef WITHOUT_COMPLEX
void Recv( scomplex* buf, int count, int from, int tag, Comm comm );
void Recv( dcomplex* buf, int count, int from, int tag, Comm comm );
#endif

void IRecv
( byte* buf, int count, int from, int tag, Comm comm, Request& request );
void IRecv
( int* buf, int count, int from, int tag, Comm comm, Request& request );
void IRecv
( float* buf, int count, int from, int tag, Comm comm, Request& request );
void IRecv
( double* buf, int count, int from, int tag, Comm comm, Request& request );
#ifndef WITHOUT_COMPLEX
void IRecv
( scomplex* buf, int count, int from, int tag, Comm comm, Request& request );
void IRecv
( dcomplex* buf, int count, int from, int tag, Comm comm, Request& request );
#endif

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
#ifndef WITHOUT_COMPLEX
void SendRecv
( const scomplex* sbuf, int sc, int to,   int stag,
        scomplex* rbuf, int rc, int from, int rtag, Comm comm );
void SendRecv
( const dcomplex* sbuf, int sc, int to,   int stag,
        dcomplex* rbuf, int rc, int from, int rtag, Comm comm );
#endif

// Collective communication

void Broadcast( byte* buf, int count, int root, Comm comm );
void Broadcast( int* buf, int count, int root, Comm comm );
void Broadcast( float* buf, int count, int root, Comm comm );
void Broadcast( double* buf, int count, int root, Comm comm );
#ifndef WITHOUT_COMPLEX
void Broadcast( scomplex* buf, int count, int root, Comm comm );
void Broadcast( dcomplex* buf, int count, int root, Comm comm );
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
#ifndef WITHOUT_COMPLEX
void  Gather
( const scomplex* sbuf, int sc,
        scomplex* rbuf, int rc, int root, Comm comm );
void Gather
( const dcomplex* sbuf, int sc,
        dcomplex* rbuf, int rc, int root, Comm comm );
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
#ifndef WITHOUT_COMPLEX
void AllGather
( const scomplex* sbuf, int sc,
        scomplex* rbuf, int rc, Comm comm );
void AllGather
( const dcomplex* sbuf, int sc,
        dcomplex* rbuf, int rc, Comm comm );
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
#ifndef WITHOUT_COMPLEX
void Scatter
( const scomplex* sbuf, int sc,
        scomplex* rbuf, int rc, int root, Comm comm );
void Scatter
( const dcomplex* sbuf, int sc,
        dcomplex* rbuf, int rc, int root, Comm comm );
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
#ifndef WITHOUT_COMPLEX
void AllToAll
( const scomplex* sbuf, int sc,
        scomplex* rbuf, int rc, Comm comm );
void AllToAll
( const dcomplex* sbuf, int sc,
        dcomplex* rbuf, int rc, Comm comm );
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
#ifndef WITHOUT_COMPLEX
void AllToAll
( const scomplex* sbuf, const int* scs, const int* sds,
        scomplex* rbuf, const int* rcs, const int* rds, Comm comm );
void AllToAll
( const dcomplex* sbuf, const int* scs, const int* sds,
        dcomplex* rbuf, const int* rcs, const int* rds, Comm comm );
#endif

void Reduce
( const byte* sbuf, byte* rbuf, int count, Op op, int root, Comm comm );
void Reduce
( const int* sbuf, int* rbuf, int count, Op op, int root, Comm comm );
void Reduce
( const float* sbuf, float* rbuf, int count, Op op, int root, Comm comm );
void Reduce
( const double* sbuf, double* rbuf, int count, Op op, int root, Comm comm );
#ifndef WITHOUT_COMPLEX
void Reduce
( const scomplex* sbuf, scomplex* rbuf, int count, Op op, int root, Comm comm );
void Reduce
( const dcomplex* sbuf, dcomplex* rbuf, int count, Op op, int root, Comm comm );
#endif

// In-place option
void Reduce( byte* buf, int count, Op op, int root, Comm comm );
void Reduce( int* buf, int count, Op op, int root, Comm comm );
void Reduce( float* buf, int count, Op op, int root, Comm comm );
void Reduce( double* buf, int count, Op op, int root, Comm comm );
#ifndef WITHOUT_COMPLEX
void Reduce( scomplex* buf, int count, Op op, int root, Comm comm );
void Reduce( dcomplex* buf, int count, Op op, int root, Comm comm );
#endif
    
void AllReduce
( const byte* sbuf, byte* rbuf, int count, Op op, Comm comm );
void AllReduce
( const int* sbuf, int* rbuf, int count, Op op, Comm comm );
void AllReduce
( const float* sbuf, float* rbuf, int count, Op op, Comm comm );
void AllReduce
( const double* sbuf, double* rbuf, int count, Op op, Comm comm );
#ifndef WITHOUT_COMPLEX
void AllReduce
( const scomplex* sbuf, scomplex* rbuf, int count, Op op, Comm comm );
void AllReduce
( const dcomplex* sbuf, dcomplex* rbuf, int count, Op op, Comm comm );
#endif

// In-place option
void AllReduce( byte* buf, int count, Op op, Comm comm );
void AllReduce( int* buf, int count, Op op, Comm comm );
void AllReduce( float* buf, int count, Op op, Comm comm );
void AllReduce( double* buf, int count, Op op, Comm comm );
#ifndef WITHOUT_COMPLEX
void AllReduce( scomplex* buf, int count, Op op, Comm comm );
void AllReduce( dcomplex* buf, int count, Op op, Comm comm );
#endif
        
void ReduceScatter
( const byte* sbuf, byte* rbuf, const int* rcs, Op op, Comm comm );
void ReduceScatter
( const int* sbuf, int* rbuf, const int* rcs, Op op, Comm comm );
void ReduceScatter
( const float* sbuf, float* rbuf, const int* rcs, Op op, Comm comm );
void ReduceScatter
( const double* sbuf, double* rbuf, const int* rcs, Op op, Comm comm );
#ifndef WITHOUT_COMPLEX
void ReduceScatter
( const scomplex* sbuf, scomplex* rbuf, const int* rcs, Op op, Comm comm );
void ReduceScatter
( const dcomplex* sbuf, dcomplex* rbuf, const int* rcs, Op op, Comm comm );
#endif

} // mpi
} // elemental

#endif /* ELEMENTAL_IMPORTS_MPI_HPP */

