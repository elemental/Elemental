/*
   Copyright (c) 2009-2013, Jack Poulson
                      2013, Jeff Hammond
                      2013, Jed Brown
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental-lite.hpp"

namespace {

inline void 
SafeMpi( int mpiError )
{
#ifndef RELEASE
    if( mpiError != MPI_SUCCESS )    
    {
        char errorString[200];
        int lengthOfErrorString;
        MPI_Error_string( mpiError, errorString, &lengthOfErrorString );
        throw std::logic_error( errorString );
    }
#endif
}

} // anonymous namespace

namespace elem {
namespace mpi {

// NOTE: This data structure is inspired by Justin Holewinski's blog post at
//       http://jholewinski.org/blog/the-beauty-of-c-templates/
template<typename T>
struct MpiMap
{
    Datatype type;
    MpiMap();
};

template<>
MpiMap<byte>::MpiMap() : type(MPI_UNSIGNED_CHAR) { }

template<>
MpiMap<int>::MpiMap() : type(MPI_INT) { }

template<>
MpiMap<float>::MpiMap() : type(MPI_FLOAT) { }

template<>
MpiMap<double>::MpiMap() : type(MPI_DOUBLE) { }

template<>
MpiMap<Complex<float> >::MpiMap() : type(MPI_COMPLEX) { }

template<>
MpiMap<Complex<double> >::MpiMap() : type(MPI_DOUBLE_COMPLEX) { }

//----------------------------//
// MPI environmental routines //
//----------------------------//

void Initialize( int& argc, char**& argv )
{ MPI_Init( &argc, &argv ); }

int InitializeThread( int& argc, char**& argv, int required )
{ 
    int provided; 
    MPI_Init_thread( &argc, &argv, required, &provided ); 
    return provided;
}

void Finalize()
{ MPI_Finalize(); }

bool Initialized()
{ 
    int initialized;
    MPI_Initialized( &initialized );
    return initialized;
}

bool Finalized()
{
    int finalized;
    MPI_Finalized( &finalized );
    return finalized;
}

int QueryThread()
{
    int provided;
    MPI_Query_thread( &provided );
    return provided;
}

double Time()
{ return MPI_Wtime(); }

void OpCreate( UserFunction* func, bool commutes, Op& op )
{
#ifndef RELEASE
    PushCallStack("mpi::OpCreate");
#endif
    SafeMpi( MPI_Op_create( func, commutes, &op ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

void OpFree( Op& op )
{
#ifndef RELEASE
    PushCallStack("mpi::OpFree");
#endif
    SafeMpi( MPI_Op_free( &op ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

//---------------------------//
// Communicator manipulation //
//---------------------------//

int WorldRank()
{
#ifndef RELEASE
    PushCallStack("mpi::WorldRank");
#endif
    int rank;
    SafeMpi( MPI_Comm_rank( COMM_WORLD, &rank ) );
#ifndef RELEASE
    PopCallStack();
#endif
    return rank;
}

int CommRank( Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::CommRank");
#endif
    int rank;
    SafeMpi( MPI_Comm_rank( comm, &rank ) );
#ifndef RELEASE
    PopCallStack();
#endif
    return rank;
}

int CommSize( Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::CommSize");
#endif
    int size;
    SafeMpi( MPI_Comm_size( comm, &size ) );
#ifndef RELEASE
    PopCallStack();
#endif
    return size;
}

void CommCreate( Comm parentComm, Group subsetGroup, Comm& subsetComm )
{
#ifndef RELEASE
    PushCallStack("mpi::CommCreate");
#endif
    SafeMpi( MPI_Comm_create( parentComm, subsetGroup, &subsetComm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

void CommDup( Comm original, Comm& duplicate )
{
#ifndef RELEASE
    PushCallStack("mpi::CommDup");
#endif
    SafeMpi( MPI_Comm_dup( original, &duplicate ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

void CommSplit( Comm comm, int color, int key, Comm& newComm )
{
#ifndef RELEASE
    PushCallStack("mpi::CommSplit");
#endif
    SafeMpi( MPI_Comm_split( comm, color, key, &newComm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

void CommFree( Comm& comm )
{
#ifndef RELEASE
    PushCallStack("mpi::CommFree");
#endif
    SafeMpi( MPI_Comm_free( &comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

bool CongruentComms( Comm comm1, Comm comm2 )
{
#ifndef RELEASE
    PushCallStack("mpi::CongruentComms");
#endif
    int result;
    SafeMpi( MPI_Comm_compare( comm1, comm2, &result ) );
#ifndef RELEASE
    PopCallStack();
#endif
    return ( result == MPI_IDENT || result == MPI_CONGRUENT );
}

void ErrorHandlerSet( Comm comm, ErrorHandler errorHandler )
{
#ifndef RELEASE
    PushCallStack("mpi::ErrorHandlerSet");
#endif
#ifdef HAVE_MPI_COMM_SET_ERRHANDLER
    SafeMpi( MPI_Comm_set_errhandler( comm, errorHandler ) );
#else
    SafeMpi( MPI_Errhandler_set( comm, errorHandler ) );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

//---------------------------------//
// Cartesian communicator routines //
//---------------------------------//

void CartCreate
( Comm comm, int numDims, const int* dimensions, const int* periods, 
  bool reorder, Comm& cartComm )
{
#ifndef RELEASE
    PushCallStack("mpi::CartCreate");
#endif
    SafeMpi( 
        MPI_Cart_create
        ( comm, numDims, const_cast<int*>(dimensions), 
          const_cast<int*>(periods), reorder, &cartComm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void CartSub( Comm comm, const int* remainingDims, Comm& subComm )
{
#ifndef RELEASE
    PushCallStack("mpi::CartSub");
#endif
    SafeMpi( MPI_Cart_sub( comm, const_cast<int*>(remainingDims), &subComm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

//--------------------//
// Group manipulation //
//--------------------//

int GroupRank( Group group )
{
#ifndef RELEASE
    PushCallStack("mpi::GroupRank");
#endif
    int rank;
    SafeMpi( MPI_Group_rank( group, &rank ) );
#ifndef RELEASE
    PopCallStack();
#endif
    return rank;
}

int GroupSize( Group group )
{
#ifndef RELEASE
    PushCallStack("mpi::GroupSize");
#endif
    int size;
    SafeMpi( MPI_Group_size( group, &size ) );
#ifndef RELEASE
    PopCallStack();
#endif
    return size;
}

void CommGroup( Comm comm, Group& group )
{
#ifndef RELEASE
    PushCallStack("mpi::CommGroup");
#endif
    SafeMpi( MPI_Comm_group( comm, &group ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

void GroupIncl( Group group, int n, const int* ranks, Group& subGroup )
{
#ifndef RELEASE
    PushCallStack("mpi::GroupIncl");
#endif
    SafeMpi( MPI_Group_incl( group, n, const_cast<int*>(ranks), &subGroup ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

void GroupDifference( Group parent, Group subset, Group& complement )
{
#ifndef RELEASE
    PushCallStack("mpi::GroupDifference");
#endif
    SafeMpi( MPI_Group_difference( parent, subset, &complement ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

void GroupFree( Group& group )
{
#ifndef RELEASE
    PushCallStack("mpi::GroupFree");
#endif
    SafeMpi( MPI_Group_free( &group ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

void GroupTranslateRanks
( Group origGroup, int size, const int* origRanks, 
  Group newGroup,                  int* newRanks )
{
#ifndef RELEASE
    PushCallStack("mpi::GroupTranslateRanks");
#endif
    SafeMpi( 
        MPI_Group_translate_ranks
        ( origGroup, size, const_cast<int*>(origRanks), newGroup, newRanks ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

// Wait until every process in comm reaches this statement
void Barrier( Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Barrier");
#endif
    SafeMpi( MPI_Barrier( comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

// Test for completion
bool Test( Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::Test");
#endif
    Status status;
    int flag;
    SafeMpi( MPI_Test( &request, &flag, &status ) );
#ifndef RELEASE
    PopCallStack();
#endif
    return flag;
}

// Ensure that the request finishes before continuing
void Wait( Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::Wait");
#endif
    Status status;
    SafeMpi( MPI_Wait( &request, &status ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

// Ensure that the request finishes before continuing
void Wait( Request& request, Status& status )
{
#ifndef RELEASE
    PushCallStack("mpi::Wait");
#endif
    SafeMpi( MPI_Wait( &request, &status ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

// Ensure that several requests finish before continuing
void WaitAll( int numRequests, Request* requests, Status* statuses )
{
#ifndef RELEASE
    PushCallStack("mpi::WaitAll");
#endif
    SafeMpi( MPI_Waitall( numRequests, requests, statuses ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

// Nonblocking test for message completion
bool IProbe( int source, int tag, Comm comm, Status& status )
{
#ifndef RELEASE
    PushCallStack("mpi::IProbe");
#endif
    int flag;
    SafeMpi( MPI_Iprobe( source, tag, comm, &flag, &status ) );
#ifndef RELEASE
    PopCallStack();
#endif
    return flag;
}

template<typename T>
int GetCount( Status& status )
{
#ifndef RELEASE
    PushCallStack("mpi::GetCount");
#endif
    int count;
    MpiMap<T> map;
    SafeMpi( MPI_Get_count( &status, map.type, &count ) );
#ifndef RELEASE
    PopCallStack();
#endif
    return count;
}
template int GetCount<byte>( Status& status );
template int GetCount<int>( Status& status );
template int GetCount<float>( Status& status );
template int GetCount<double>( Status& status );
template int GetCount<Complex<float> >( Status& status );
template int GetCount<Complex<double> >( Status& status );

template<typename R>
void Send( const R* buf, int count, int to, int tag, Comm comm )
{ 
#ifndef RELEASE
    PushCallStack("mpi::Send");
#endif
    MpiMap<R> map;
    SafeMpi( MPI_Send( const_cast<R*>(buf), count, map.type, to, tag, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void Send( const Complex<R>* buf, int count, int to, int tag, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Send");
#endif
#ifdef AVOID_COMPLEX_MPI
    MpiMap<R> map;
    SafeMpi(
        MPI_Send
        ( const_cast<Complex<R>*>(buf), 2*count, map.type, to, tag, comm )
    );
#else
    MpiMap<Complex<R> > map;
    SafeMpi( 
        MPI_Send
        ( const_cast<Complex<R>*>(buf), count, map.type, to, tag, comm )
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

template void Send( const byte* buf, int count, int to, int tag, Comm comm );
template void Send( const int* buf, int count, int to, int tag, Comm comm );
template void Send( const float* buf, int count, int to, int tag, Comm comm );
template void Send( const double* buf, int count, int to, int tag, Comm comm );
template void Send( const Complex<float>* buf, int count, int to, int tag, Comm comm );
template void Send( const Complex<double>* buf, int count, int to, int tag, Comm comm );

template<typename R>
void ISend
( const R* buf, int count, int to, int tag, Comm comm, Request& request )
{ 
#ifndef RELEASE
    PushCallStack("mpi::ISend");
#endif
    MpiMap<R> map;
    SafeMpi( 
        MPI_Isend
        ( const_cast<R*>(buf), count, map.type, to, tag, comm, &request )
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void ISend
( const Complex<R>* buf, int count, int to, int tag, Comm comm, 
  Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::ISend");
#endif
#ifdef AVOID_COMPLEX_MPI
    MpiMap<R> map;
    SafeMpi(
        MPI_Isend
        ( const_cast<Complex<R>*>(buf), 2*count, map.type, to, tag, comm,
          &request )
    );
#else
    MpiMap<Complex<R> > map;
    SafeMpi( 
        MPI_Isend
        ( const_cast<Complex<R>*>(buf), count, map.type, to, tag, comm,
          &request )
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

template void ISend( const byte* buf, int count, int to, int tag, Comm comm, Request& request );
template void ISend( const int* buf, int count, int to, int tag, Comm comm, Request& request );
template void ISend( const float* buf, int count, int to, int tag, Comm comm, Request& request );
template void ISend( const double* buf, int count, int to, int tag, Comm comm, Request& request );
template void ISend( const Complex<float>* buf, int count, int to, int tag, Comm comm, Request& request );
template void ISend( const Complex<double>* buf, int count, int to, int tag, Comm comm, Request& request );

template<typename R>
void ISSend
( const R* buf, int count, int to, int tag, Comm comm, Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::ISSend");
#endif
    MpiMap<R> map;
    SafeMpi(
        MPI_Issend
        ( const_cast<R*>(buf), count, map.type, to, tag, comm, &request )
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void ISSend
( const Complex<R>* buf, int count, int to, int tag, Comm comm, 
  Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::ISSend");
#endif
#ifdef AVOID_COMPLEX_MPI
    MpiMap<R> map;
    SafeMpi(
        MPI_Issend
        ( const_cast<Complex<R>*>(buf), 2*count, map.type, to, tag, comm,
          &request )
    );
#else
    MpiMap<Complex<R> > map;
    SafeMpi(
        MPI_Issend
        ( const_cast<Complex<R>*>(buf), count, map.type, to, tag, comm,
          &request )
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

template void ISSend( const byte* buf, int count, int to, int tag, Comm comm, Request& request );
template void ISSend( const int* buf, int count, int to, int tag, Comm comm, Request& request );
template void ISSend( const float* buf, int count, int to, int tag, Comm comm, Request& request );
template void ISSend( const double* buf, int count, int to, int tag, Comm comm, Request& request );
template void ISSend( const Complex<float>* buf, int count, int to, int tag, Comm comm, Request& request );
template void ISSend( const Complex<double>* buf, int count, int to, int tag, Comm comm, Request& request );

template<typename R>
void Recv( R* buf, int count, int from, int tag, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Recv");
#endif
    MpiMap<R> map;
    Status status;
    SafeMpi( MPI_Recv( buf, count, map.type, from, tag, comm, &status ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void Recv( Complex<R>* buf, int count, int from, int tag, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Recv");
#endif
    Status status;
#ifdef AVOID_COMPLEX_MPI
    MpiMap<R> map;
    SafeMpi( MPI_Recv( buf, 2*count, map.type, from, tag, comm, &status ) );
#else
    MpiMap<Complex<R> > map;
    SafeMpi( MPI_Recv( buf, count, map.type, from, tag, comm, &status ) );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

template void Recv( byte* buf, int count, int from, int tag, Comm comm );
template void Recv( int* buf, int count, int from, int tag, Comm comm );
template void Recv( float* buf, int count, int from, int tag, Comm comm );
template void Recv( double* buf, int count, int from, int tag, Comm comm );
template void Recv( Complex<float>* buf, int count, int from, int tag, Comm comm );
template void Recv( Complex<double>* buf, int count, int from, int tag, Comm comm );

template<typename R>
void IRecv( R* buf, int count, int from, int tag, Comm comm, Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::IRecv");
#endif
    MpiMap<R> map;
    SafeMpi( MPI_Irecv( buf, count, map.type, from, tag, comm, &request ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void IRecv
( Complex<R>* buf, int count, int from, int tag, Comm comm, Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::IRecv");
#endif
#ifdef AVOID_COMPLEX_MPI
    MpiMap<R> map;
    SafeMpi( MPI_Irecv( buf, 2*count, map.type, from, tag, comm, &request ) );
#else
    MpiMap<Complex<R> > map;
    SafeMpi( MPI_Irecv( buf, count, map.type, from, tag, comm, &request ) );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

template void IRecv( byte* buf, int count, int from, int tag, Comm comm, Request& request );
template void IRecv( int* buf, int count, int from, int tag, Comm comm, Request& request );
template void IRecv( float* buf, int count, int from, int tag, Comm comm, Request& request );
template void IRecv( double* buf, int count, int from, int tag, Comm comm, Request& request );
template void IRecv( Complex<float>* buf, int count, int from, int tag, Comm comm, Request& request );
template void IRecv( Complex<double>* buf, int count, int from, int tag, Comm comm, Request& request );

template<typename R>
void SendRecv
( const R* sbuf, int sc, int to,   int stag,
        R* rbuf, int rc, int from, int rtag, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::SendRecv");
#endif
    Status status;
    MpiMap<R> map;
    SafeMpi( 
        MPI_Sendrecv
        ( const_cast<R*>(sbuf), sc, map.type, to,   stag,
          rbuf,                 rc, map.type, from, rtag, comm, &status )
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void SendRecv
( const Complex<R>* sbuf, int sc, int to,   int stag,
        Complex<R>* rbuf, int rc, int from, int rtag, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::SendRecv");
#endif
    Status status;
#ifdef AVOID_COMPLEX_MPI
    MpiMap<R> map;
    SafeMpi(
        MPI_Sendrecv
        ( const_cast<Complex<R>*>(sbuf), 2*sc, map.type, to,   stag,
          rbuf,                          2*rc, map.type, from, rtag, 
          comm, &status )
    );
#else
    MpiMap<Complex<R> > map;
    SafeMpi( 
        MPI_Sendrecv
        ( const_cast<Complex<R>*>(sbuf), sc, map.type, to,   stag,
          rbuf,                          rc, map.type, from, rtag, 
          comm, &status )
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

template void SendRecv
( const byte* sbuf, int sc, int to, int stag, 
        byte* rbuf, int rc, int from, int rtag, Comm comm );
template void SendRecv
( const int* sbuf, int sc, int to, int stag, 
        int* rbuf, int rc, int from, int rtag, Comm comm );
template void SendRecv
( const float* sbuf, int sc, int to, int stag, 
        float* rbuf, int rc, int from, int rtag, Comm comm );
template void SendRecv
( const double* sbuf, int sc, int to, int stag, 
        double* rbuf, int rc, int from, int rtag, Comm comm );
template void SendRecv
( const Complex<float>* sbuf, int sc, int to, int stag, 
        Complex<float>* rbuf, int rc, int from, int rtag, Comm comm );
template void SendRecv
( const Complex<double>* sbuf, int sc, int to, int stag, 
        Complex<double>* rbuf, int rc, int from, int rtag, Comm comm );

template<typename R>
void Broadcast( R* buf, int count, int root, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Broadcast");
#endif
    MpiMap<R> map;
    SafeMpi( MPI_Bcast( buf, count, map.type, root, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void Broadcast( Complex<R>* buf, int count, int root, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Broadcast");
#endif
#ifdef AVOID_COMPLEX_MPI
    MpiMap<R> map;
    SafeMpi( MPI_Bcast( buf, 2*count, map.type, root, comm ) );
#else
    MpiMap<Complex<R> > map;
    SafeMpi( MPI_Bcast( buf, count, map.type, root, comm ) );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

template void Broadcast( byte* buf, int count, int root, Comm comm );
template void Broadcast( int* buf, int count, int root, Comm comm );
template void Broadcast( float* buf, int count, int root, Comm comm );
template void Broadcast( double* buf, int count, int root, Comm comm );
template void Broadcast( Complex<float>* buf, int count, int root, Comm comm );
template void Broadcast( Complex<double>* buf, int count, int root, Comm comm );

#ifdef HAVE_NONBLOCKING_COLLECTIVES
template<typename R>
void IBroadcast( R* buf, int count, int root, Comm comm, Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::IBroadcast");
#endif
    MpiMap<R> map;
    SafeMpi( MPI_Ibcast( buf, count, map.type, root, comm, &request ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void IBroadcast
( Complex<R>* buf, int count, int root, Comm comm, Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::IBroadcast");
#endif
#ifdef AVOID_COMPLEX_MPI
    MpiMap<R> map;
    SafeMpi( MPI_Ibcast( buf, 2*count, map.type, root, comm, &request ) );
#else
    MpiMap<Complex<R> > map;
    SafeMpi( MPI_Ibcast( buf, count, map.type, root, comm, &request ) );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

template void IBroadcast( byte* buf, int count, int root, Comm comm, Request& request );
template void IBroadcast( int* buf, int count, int root, Comm comm, Request& request );
template void IBroadcast( float* buf, int count, int root, Comm comm, Request& request );
template void IBroadcast( double* buf, int count, int root, Comm comm, Request& request );
template void IBroadcast( Complex<float>* buf, int count, int root, Comm comm, Request& request );
template void IBroadcast( Complex<double>* buf, int count, int root, Comm comm, Request& request );
#endif // ifdef HAVE_NONBLOCKING_COLLECTIVES

template<typename R>
void Gather
( const R* sbuf, int sc,
        R* rbuf, int rc, int root, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Gather");
#endif
    MpiMap<R> map;
    SafeMpi( 
        MPI_Gather
        ( const_cast<R*>(sbuf), sc, map.type,
          rbuf,                 rc, map.type, root, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void Gather
( const Complex<R>* sbuf, int sc,
        Complex<R>* rbuf, int rc, int root, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Gather");
#endif
#ifdef AVOID_COMPLEX_MPI
    MpiMap<R> map;
    SafeMpi(
        MPI_Gather
        ( const_cast<Complex<R>*>(sbuf), 2*sc, map.type,
          rbuf,                          2*rc, map.type, root, comm )
    );
#else
    MpiMap<Complex<R> > map;
    SafeMpi( 
        MPI_Gather
        ( const_cast<Complex<R>*>(sbuf), sc, map.type,
          rbuf,                          rc, map.type, root, comm ) 
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

template void Gather( const byte* sbuf, int sc, byte* rbuf, int rc, int root, Comm comm );
template void Gather( const int* sbuf, int sc, int* rbuf, int rc, int root, Comm comm );
template void Gather( const float* sbuf, int sc, float* rbuf, int rc, int root, Comm comm );
template void Gather( const double* sbuf, int sc, double* rbuf, int rc, int root, Comm comm );
template void Gather( const Complex<float>* sbuf, int sc, Complex<float>* rbuf, int rc, int root, Comm comm );
template void Gather( const Complex<double>* sbuf, int sc, Complex<double>* rbuf, int rc, int root, Comm comm );

#ifdef HAVE_NONBLOCKING_COLLECTIVES
template<typename R>
void IGather
( const R* sbuf, int sc,
        R* rbuf, int rc, int root, Comm comm, Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::IGather");
#endif
    MpiMap<R> map;
    SafeMpi( 
        MPI_Igather
        ( const_cast<R*>(sbuf), sc, map.type,
          rbuf,                 rc, map.type, root, comm, &request )
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void IGather
( const Complex<R>* sbuf, int sc,
        Complex<R>* rbuf, int rc, int root, Comm comm, Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::IGather");
#endif
#ifdef AVOID_COMPLEX_MPI
    MpiMap<R> map;
    SafeMpi(
        MPI_Igather
        ( const_cast<Complex<R>*>(sbuf), 2*sc, map.type,
          rbuf,                          2*rc, map.type, 
          root, comm, &request )
    );
#else
    MpiMap<Complex<R> > map;
    SafeMpi( 
        MPI_Igather
        ( const_cast<Complex<R>*>(sbuf), sc, map.type,
          rbuf,                          rc, map.type, 
          root, comm, &request ) 
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

template void IGather
( const byte* sbuf, int sc, 
        byte* rbuf, int rc, int root, Comm comm, Request& request );
template void IGather
( const int* sbuf, int sc, 
        int* rbuf, int rc, int root, Comm comm, Request& request );
template void IGather
( const float* sbuf, int sc, 
        float* rbuf, int rc, int root, Comm comm, Request& request );
template void IGather
( const double* sbuf, int sc, 
        double* rbuf, int rc, int root, Comm comm, Request& request );
template void IGather
( const Complex<float>* sbuf, int sc, 
        Complex<float>* rbuf, int rc, int root, Comm comm, Request& request );
template void IGather
( const Complex<double>* sbuf, int sc, 
        Complex<double>* rbuf, int rc, int root, Comm comm, Request& request );
#endif // ifdef HAVE_NONBLOCKING_COLLECTIVES

template<typename R>
void Gather
( const R* sbuf, int sc,
        R* rbuf, const int* rcs, const int* rds, int root, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Gather");
#endif
    MpiMap<R> map;
    SafeMpi( 
        MPI_Gatherv
        ( const_cast<R*>(sbuf), 
          sc,       
          map.type,
          rbuf,                    
          const_cast<int*>(rcs), 
          const_cast<int*>(rds), 
          map.type,
          root, 
          comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void Gather
( const Complex<R>* sbuf, int sc,
        Complex<R>* rbuf, const int* rcs, const int* rds, int root, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Gather");
#endif
#ifdef AVOID_COMPLEX_MPI
    MpiMap<R> map;
    const int commRank = CommRank( comm );
    const int commSize = CommSize( comm );
    std::vector<int> rcsDouble, rdsDouble;
    if( commRank == root )
    {
        rcsDouble.resize( commSize );
        rdsDouble.resize( commSize );
        for( int i=0; i<commSize; ++i )
        {
            rcsDouble[i] = 2*rcs[i];
            rdsDouble[i] = 2*rds[i];
        }
    }
    SafeMpi(
        MPI_Gatherv
        ( const_cast<Complex<R>*>(sbuf), 2*sc,                         map.type,
          rbuf,                          &rcsDouble[0], &rdsDouble[0], map.type,
          root, comm )
    );
#else
    MpiMap<Complex<R> > map;
    SafeMpi( 
        MPI_Gatherv
        ( const_cast<Complex<R>*>(sbuf), 
          sc,       
          map.type,
          rbuf,  
          const_cast<int*>(rcs), 
          const_cast<int*>(rds), 
          map.type,
          root, 
          comm ) 
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

template void Gather
( const byte* sbuf, int sc, 
        byte* rbuf, const int* rcs, const int* rds, int root, Comm comm );
template void Gather
( const int* sbuf, int sc, 
        int* rbuf, const int* rcs, const int* rds, int root, Comm comm );
template void Gather
( const float* sbuf, int sc, 
        float* rbuf, const int* rcs, const int* rds, int root, Comm comm );
template void Gather
( const double* sbuf, int sc, 
        double* rbuf, const int* rcs, const int* rds, int root, Comm comm );
template void Gather
( const Complex<float>* sbuf, int sc, 
        Complex<float>* rbuf, const int* rcs, const int* rds, int root, Comm comm );
template void Gather
( const Complex<double>* sbuf, int sc, 
        Complex<double>* rbuf, const int* rcs, const int* rds, int root, Comm comm );

template<typename R>
void AllGather
( const R* sbuf, int sc,
        R* rbuf, int rc, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllGather");
#endif
#ifdef USE_BYTE_ALLGATHERS
    SafeMpi( 
        MPI_Allgather
        ( const_cast<R*>(sbuf), sizeof(R)*sc, MPI_UNSIGNED_CHAR, 
          rbuf,                 sizeof(R)*rc, MPI_UNSIGNED_CHAR, comm ) 
    );
#else
    MpiMap<R> map;
    SafeMpi( 
        MPI_Allgather
        ( const_cast<R*>(sbuf), sc, map.type, 
          rbuf,                 rc, map.type, comm ) 
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void AllGather
( const Complex<R>* sbuf, int sc,
        Complex<R>* rbuf, int rc, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllGather");
#endif
#ifdef USE_BYTE_ALLGATHERS
    SafeMpi( 
        MPI_Allgather
        ( const_cast<Complex<R>*>(sbuf), 2*sizeof(R)*sc, MPI_UNSIGNED_CHAR, 
          rbuf,                          2*sizeof(R)*rc, MPI_UNSIGNED_CHAR, 
          comm ) 
    );
#else
 #ifdef AVOID_COMPLEX_MPI
    MpiMap<R> map;
    SafeMpi(
        MPI_Allgather
        ( const_cast<Complex<R>*>(sbuf), 2*sc, map.type,
          rbuf,                          2*rc, map.type, comm )
    );
 #else
    MpiMap<Complex<R> > map;
    SafeMpi( 
        MPI_Allgather
        ( const_cast<Complex<R>*>(sbuf), sc, map.type,
          rbuf,                          rc, map.type, comm ) 
    );
 #endif
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

template void AllGather( const byte* sbuf, int sc, byte* rbuf, int rc, Comm comm );
template void AllGather( const int* sbuf, int sc, int* rbuf, int rc, Comm comm );
template void AllGather( const float* sbuf, int sc, float* rbuf, int rc, Comm comm );
template void AllGather( const double* sbuf, int sc, double* rbuf, int rc, Comm comm );
template void AllGather( const Complex<float>* sbuf, int sc, Complex<float>* rbuf, int rc, Comm comm );
template void AllGather( const Complex<double>* sbuf, int sc, Complex<double>* rbuf, int rc, Comm comm );

template<typename R>
void AllGather
( const R* sbuf, int sc,
        R* rbuf, const int* rcs, const int* rds, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllGather");
#endif
#ifdef USE_BYTE_ALLGATHERS
    const int commSize = CommSize( comm );
    std::vector<int> byteRcs( commSize ), byteRds( commSize );
    for( int i=0; i<commSize; ++i )
    {
        byteRcs[i] = sizeof(R)*rcs[i];
        byteRds[i] = sizeof(R)*rds[i];
    }
    SafeMpi( 
        MPI_Allgatherv
        ( const_cast<R*>(sbuf), sizeof(R)*sc, MPI_UNSIGNED_CHAR, 
          rbuf, &byteRcs[0], &byteRds[0],     MPI_UNSIGNED_CHAR, comm ) 
    );
#else
    MpiMap<R> map;
    SafeMpi( 
        MPI_Allgatherv
        ( const_cast<R*>(sbuf), 
          sc, 
          map.type, 
          rbuf,   
          const_cast<int*>(rcs), 
          const_cast<int*>(rds), 
          map.type, 
          comm ) 
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void AllGather
( const Complex<R>* sbuf, int sc,
        Complex<R>* rbuf, const int* rcs, const int* rds, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllGather");
#endif
#ifdef USE_BYTE_ALLGATHERS
    const int commSize = CommSize( comm );
    std::vector<int> byteRcs( commSize ), byteRds( commSize );
    for( int i=0; i<commSize; ++i )
    {
        byteRcs[i] = 2*sizeof(R)*rcs[i];
        byteRds[i] = 2*sizeof(R)*rds[i];
    }
    SafeMpi( 
        MPI_Allgatherv
        ( const_cast<Complex<R>*>(sbuf), 2*sizeof(R)*sc, MPI_UNSIGNED_CHAR, 
          rbuf, &byteRcs[0], &byteRds[0], MPI_UNSIGNED_CHAR, comm )
    );
#else
 #ifdef AVOID_COMPLEX_MPI
    MpiMap<R> map;
    const int commSize = CommSize( comm );
    std::vector<int> realRcs( commSize ), realRds( commSize );
    for( int i=0; i<commSize; ++i )
    {
        realRcs[i] = 2*rcs[i];
        realRds[i] = 2*rds[i];
    }
    SafeMpi(
        MPI_Allgatherv
        ( const_cast<Complex<R>*>(sbuf), 2*sc, map.type,
          rbuf, &realRcs[0], &realRds[0],      map.type, comm )
    );
 #else
    MpiMap<Complex<R> > map;
    SafeMpi( 
        MPI_Allgatherv
        ( const_cast<Complex<R>*>(sbuf), 
          sc, 
          map.type,
          rbuf, 
          const_cast<int*>(rcs), 
          const_cast<int*>(rds), 
          map.type,
          comm ) 
    );
 #endif
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

template void AllGather
( const byte* sbuf, int sc, 
        byte* rbuf, const int* rcs, const int* rds, Comm comm );
template void AllGather
( const int* sbuf, int sc, 
        int* rbuf, const int* rcs, const int* rds, Comm comm );
template void AllGather
( const float* sbuf, int sc, 
        float* rbuf, const int* rcs, const int* rds, Comm comm );
template void AllGather
( const double* sbuf, int sc, 
        double* rbuf, const int* rcs, const int* rds, Comm comm );
template void AllGather
( const Complex<float>* sbuf, int sc, 
        Complex<float>* rbuf, const int* rcs, const int* rds, Comm comm );
template void AllGather
( const Complex<double>* sbuf, int sc, 
        Complex<double>* rbuf, const int* rcs, const int* rds, Comm comm );

template<typename R>
void Scatter
( const R* sbuf, int sc,
        R* rbuf, int rc, int root, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Scatter");
#endif
    MpiMap<R> map;
    SafeMpi( 
        MPI_Scatter
        ( const_cast<R*>(sbuf), sc, map.type,
          rbuf,                 rc, map.type, root, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void Scatter
( const Complex<R>* sbuf, int sc,
        Complex<R>* rbuf, int rc, int root, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Scatter");
#endif
#ifdef AVOID_COMPLEX_MPI
    MpiMap<R> map;
    SafeMpi(
        MPI_Scatter
        ( const_cast<Complex<R>*>(sbuf), 2*sc, map.type,
          rbuf,                          2*rc, map.type, root, comm )
    );
#else
    MpiMap<Complex<R> > map;
    SafeMpi( 
        MPI_Scatter
        ( const_cast<Complex<R>*>(sbuf), sc, map.type,
          rbuf,                          rc, map.type, root, comm ) 
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

template void Scatter
( const byte* sbuf, int sc, 
        byte* rbuf, int rc, int root, Comm comm );
template void Scatter
( const int* sbuf, int sc, 
        int* rbuf, int rc, int root, Comm comm );
template void Scatter
( const float* sbuf, int sc, 
        float* rbuf, int rc, int root, Comm comm );
template void Scatter
( const double* sbuf, int sc, 
        double* rbuf, int rc, int root, Comm comm );
template void Scatter
( const Complex<float>* sbuf, int sc, 
        Complex<float>* rbuf, int rc, int root, Comm comm );
template void Scatter
( const Complex<double>* sbuf, int sc, 
        Complex<double>* rbuf, int rc, int root, Comm comm );


template<typename R>
void Scatter( R* buf, int sc, int rc, int root, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Scatter");
#endif
    MpiMap<R> map;
    const int commRank = CommRank( comm );
    if( commRank == root )
    {
#ifdef HAVE_MPI_IN_PLACE
        SafeMpi(
            MPI_Scatter
            ( buf,          sc, map.type, 
              MPI_IN_PLACE, rc, map.type, root, comm )
        );
#else
        const int commSize = CommSize( comm );
        std::vector<R> sendBuf( sc*commSize );
        MemCopy( &sendBuf[0], buf, sc*commSize );
        SafeMpi(
            MPI_Scatter
            ( &sendBuf[0], sc, map.type, 
              buf,         rc, map.type, root, comm )
        );
#endif
    }
    else
    {
        SafeMpi(
            MPI_Scatter
            ( 0,   sc, map.type, 
              buf, rc, map.type, root, comm )
        );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void Scatter( Complex<R>* buf, int sc, int rc, int root, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Scatter");
#endif
    const int commRank = CommRank( comm );
    if( commRank == root )
    {
#ifdef AVOID_COMPLEX_MPI
        MpiMap<R> map;
# ifdef HAVE_MPI_IN_PLACE
        SafeMpi(
            MPI_Scatter
            ( buf,          2*sc, map.type, 
              MPI_IN_PLACE, 2*rc, map.type, root, comm )
        );
# else
        const int commSize = CommSize( comm );
        std::vector<Complex<R> > sendBuf( sc*commSize );
        MemCopy( &sendBuf[0], buf, sc*commSize );
        SafeMpi(
            MPI_Scatter
            ( &sendBuf[0], 2*sc, map.type,          
              buf,         2*rc, map.type, root, comm )
        );
# endif
#else
# ifdef HAVE_MPI_IN_PLACE
        MpiMap<Complex<R> > map;
        SafeMpi(
            MPI_Scatter
            ( buf,          sc, map.type, 
              MPI_IN_PLACE, rc, map.type, root, comm )
        );
# else
        const int commSize = CommSize( comm );
        std::vector<Complex<R> > sendBuf( sc*commSize );
        MemCopy( &sendBuf[0], buf, sc*commSize );
        SafeMpi(
            MPI_Scatter
            ( &sendBuf[0], sc, map.type,
              buf,         rc, map.type, root, comm )
        );
# endif
#endif
    }
    else
    {
#ifdef AVOID_COMPLEX_MPI
        MpiMap<R> map;
        SafeMpi(
            MPI_Scatter
            ( 0,   2*sc, map.type, 
              buf, 2*rc, map.type, root, comm )
        );
#else
        MpiMap<Complex<R> > map;
        SafeMpi(
            MPI_Scatter
            ( 0,   sc, map.type, 
              buf, rc, map.type, root, comm )
        );
#endif
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void Scatter( byte* buf, int sc, int rc, int root, Comm comm );
template void Scatter( int* buf, int sc, int rc, int root, Comm comm );
template void Scatter( float* buf, int sc, int rc, int root, Comm comm );
template void Scatter( double* buf, int sc, int rc, int root, Comm comm );
template void Scatter( Complex<float>* buf, int sc, int rc, int root, Comm comm );
template void Scatter( Complex<double>* buf, int sc, int rc, int root, Comm comm );

template<typename R>
void AllToAll
( const R* sbuf, int sc,
        R* rbuf, int rc, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllToAll");
#endif
    MpiMap<R> map;
    SafeMpi( 
        MPI_Alltoall
        ( const_cast<R*>(sbuf), sc, map.type,
          rbuf,                 rc, map.type, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void AllToAll
( const Complex<R>* sbuf, int sc,
        Complex<R>* rbuf, int rc, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllToAll");
#endif
#ifdef AVOID_COMPLEX_MPI
    MpiMap<R> map;
    SafeMpi(
        MPI_Alltoall
        ( const_cast<Complex<R>*>(sbuf), 2*sc, map.type,
          rbuf,                          2*rc, map.type, comm )
    );
#else
    MpiMap<Complex<R> > map;
    SafeMpi( 
        MPI_Alltoall
        ( const_cast<Complex<R>*>(sbuf), sc, map.type,
          rbuf,                          rc, map.type, comm ) 
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

template void AllToAll
( const byte* sbuf, int sc, 
        byte* rbuf, int rc, Comm comm );
template void AllToAll
( const int* sbuf, int sc, 
        int* rbuf, int rc, Comm comm );
template void AllToAll
( const float* sbuf, int sc, 
        float* rbuf, int rc, Comm comm );
template void AllToAll
( const double* sbuf, int sc, 
        double* rbuf, int rc, Comm comm );
template void AllToAll
( const Complex<float>* sbuf, int sc, 
        Complex<float>* rbuf, int rc, Comm comm );
template void AllToAll
( const Complex<double>* sbuf, int sc, 
        Complex<double>* rbuf, int rc, Comm comm );

template<typename R>
void AllToAll
( const R* sbuf, const int* scs, const int* sds, 
        R* rbuf, const int* rcs, const int* rds, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllToAll");
#endif
    MpiMap<R> map;
    SafeMpi( 
        MPI_Alltoallv
        ( const_cast<R*>(sbuf), 
          const_cast<int*>(scs), 
          const_cast<int*>(sds), 
          map.type,
          rbuf, 
          const_cast<int*>(rcs), 
          const_cast<int*>(rds), 
          map.type,
          comm ) 
    ); 
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void AllToAll
( const Complex<R>* sbuf, const int* scs, const int* sds,
        Complex<R>* rbuf, const int* rcs, const int* rds, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllToAll");
#endif
#ifdef AVOID_COMPLEX_MPI
    MpiMap<R> map;
    int p;
    MPI_Comm_size( comm, &p );
    std::vector<int> scsDoubled(p);
    std::vector<int> sdsDoubled(p);
    std::vector<int> rcsDoubled(p);
    std::vector<int> rdsDoubled(p);
    for( int i=0; i<p; ++i )
        scsDoubled[i] = 2*scs[i];
    for( int i=0; i<p; ++i )
        sdsDoubled[i] = 2*sds[i];
    for( int i=0; i<p; ++i )
        rcsDoubled[i] = 2*rcs[i];
    for( int i=0; i<p; ++i )
        rdsDoubled[i] = 2*rds[i];
    SafeMpi(
        MPI_Alltoallv
        ( const_cast<Complex<R>*>(sbuf),
                &scsDoubled[0], &sdsDoubled[0], map.type,
          rbuf, &rcsDoubled[0], &rdsDoubled[0], map.type, comm )
    );
#else
    MpiMap<Complex<R> > map;
    SafeMpi( 
        MPI_Alltoallv
        ( const_cast<Complex<R>*>(sbuf), 
          const_cast<int*>(scs), 
          const_cast<int*>(sds), 
          map.type,
          rbuf, 
          const_cast<int*>(rcs), 
          const_cast<int*>(rds), 
          map.type,
          comm )
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

template void AllToAll
( const byte* sbuf, const int* scs, const int* sds,
        byte* rbuf, const int* rcs, const int* rds, Comm comm );
template void AllToAll
( const int* sbuf, const int* scs, const int* sds,
        int* rbuf, const int* rcs, const int* rds, Comm comm );
template void AllToAll
( const float* sbuf, const int* scs, const int* sds,
        float* rbuf, const int* rcs, const int* rds, Comm comm );
template void AllToAll
( const double* sbuf, const int* scs, const int* sds,
        double* rbuf, const int* rcs, const int* rds, Comm comm );
template void AllToAll
( const Complex<float>* sbuf, const int* scs, const int* sds,
        Complex<float>* rbuf, const int* rcs, const int* rds, Comm comm );
template void AllToAll
( const Complex<double>* sbuf, const int* scs, const int* sds,
        Complex<double>* rbuf, const int* rcs, const int* rds, Comm comm );

template<typename R>
void Reduce
( const R* sbuf, R* rbuf, int count, Op op, int root, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Reduce");
#endif
    MpiMap<R> map;
    if( count != 0 )
    {
        SafeMpi( 
            MPI_Reduce
            ( const_cast<R*>(sbuf), rbuf, count, map.type, op, root, comm )
        );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void Reduce
( const Complex<R>* sbuf, 
        Complex<R>* rbuf, int count, Op op, int root, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Reduce");
#endif
    if( count != 0 )
    {
#ifdef AVOID_COMPLEX_MPI
        if( op == SUM )
        {
            MpiMap<R> map;
            SafeMpi(
                MPI_Reduce
                ( const_cast<Complex<R>*>(sbuf),
                  rbuf, 2*count, map.type, op, root, comm )
            );
        }
        else
        {
            MpiMap<Complex<R> > map;
            SafeMpi(
                MPI_Reduce
                ( const_cast<Complex<R>*>(sbuf),
                  rbuf, count, map.type, op, root, comm )
            );
        }
#else
        MpiMap<Complex<R> > map;
        SafeMpi( 
            MPI_Reduce
            ( const_cast<Complex<R>*>(sbuf), 
              rbuf, count, map.type, op, root, comm ) 
        );
#endif
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void Reduce( const byte* sbuf, byte* rbuf, int count, Op op, int root, Comm comm );
template void Reduce( const int* sbuf, int* rbuf, int count, Op op, int root, Comm comm );
template void Reduce( const float* sbuf, float* rbuf, int count, Op op, int root, Comm comm );
template void Reduce( const double* sbuf, double* rbuf, int count, Op op, int root, Comm comm );
template void Reduce( const Complex<float>* sbuf, Complex<float>* rbuf, int count, Op op, int root, Comm comm );
template void Reduce( const Complex<double>* sbuf, Complex<double>* rbuf, int count, Op op, int root, Comm comm );

template<typename R>
void Reduce( R* buf, int count, Op op, int root, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Reduce");
#endif
    MpiMap<R> map;
    if( count != 0 )
    {
        const int commRank = CommRank( comm );
        if( commRank == root )
        {
#ifdef HAVE_MPI_IN_PLACE
            SafeMpi( 
                MPI_Reduce( MPI_IN_PLACE, buf, count, map.type, op, root, comm )
            );
#else
            std::vector<R> sendBuf( count );
            MemCopy( &sendBuf[0], buf, count );
            SafeMpi(
                MPI_Reduce( &sendBuf[0], buf, count, map.type, op, root, comm )
            );
#endif
        }
        else
            SafeMpi( MPI_Reduce( buf, 0, count, map.type, op, root, comm ) );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void Reduce( Complex<R>* buf, int count, Op op, int root, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Reduce");
#endif
    if( count != 0 )
    {
        const int commRank = CommRank( comm );
#ifdef AVOID_COMPLEX_MPI
        if( op == SUM )
        {
            MpiMap<R> map;
            if( commRank == root )
            {
# ifdef HAVE_MPI_IN_PLACE
                SafeMpi(
                    MPI_Reduce
                    ( MPI_IN_PLACE, buf, 2*count, map.type, op, root, comm )
                );
# else
                std::vector<Complex<R> > sendBuf( count );
                MemCopy( &sendBuf[0], buf, count );
                SafeMpi(
                    MPI_Reduce
                    ( &sendBuf[0], buf, 2*count, map.type, op, root, comm )
                );
# endif
            }
            else
                SafeMpi(
                    MPI_Reduce( buf, 0, 2*count, map.type, op, root, comm )
                );
        }
        else
        {
            MpiMap<Complex<R> > map;
            if( commRank == root )
            {
# ifdef HAVE_MPI_IN_PLACE
                SafeMpi(
                    MPI_Reduce
                    ( MPI_IN_PLACE, buf, count, map.type, op, root, comm )
                );
# else
                std::vector<Complex<R> > sendBuf( count );
                MemCopy( &sendBuf[0], buf, count );
                SafeMpi(
                    MPI_Reduce
                    ( &sendBuf[0], buf, count, map.type, op, root, comm )
                );
# endif
            }
            else
                SafeMpi(
                    MPI_Reduce( buf, 0, count, map.type, op, root, comm )
                );
        }
#else
        MpiMap<Complex<R> > map;
        if( commRank == root )
        {
# ifdef HAVE_MPI_IN_PLACE
            SafeMpi( 
                MPI_Reduce( MPI_IN_PLACE, buf, count, map.type, op, root, comm )
            );
# else
            std::vector<Complex<R> > sendBuf( count );
            MemCopy( &sendBuf[0], buf, count );
            SafeMpi(
                MPI_Reduce( &sendBuf[0], buf, count, map.type, op, root, comm )
            );
# endif
        }
        else
            SafeMpi( MPI_Reduce( buf, 0, count, map.type, op, root, comm ) );
#endif
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void Reduce( byte* buf, int count, Op op, int root, Comm comm );
template void Reduce( int* buf, int count, Op op, int root, Comm comm );
template void Reduce( float* buf, int count, Op op, int root, Comm comm );
template void Reduce( double* buf, int count, Op op, int root, Comm comm );
template void Reduce( Complex<float>* buf, int count, Op op, int root, Comm comm );
template void Reduce( Complex<double>* buf, int count, Op op, int root, Comm comm );

template<typename R>
void AllReduce( const R* sbuf, R* rbuf, int count, Op op, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllReduce");
#endif
    MpiMap<R> map;
    if( count != 0 )
    {
        SafeMpi( 
            MPI_Allreduce
            ( const_cast<R*>(sbuf), rbuf, count, map.type, op, comm )
        );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void AllReduce
( const Complex<R>* sbuf, Complex<R>* rbuf, int count, Op op, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllReduce");
#endif
    if( count != 0 )
    {
#ifdef AVOID_COMPLEX_MPI
        if( op == SUM )
        {
            MpiMap<R> map;
            SafeMpi(
                MPI_Allreduce
                ( const_cast<Complex<R>*>(sbuf),
                  rbuf, 2*count, map.type, op, comm )
            );
        }
        else
        {
            MpiMap<Complex<R> > map;
            SafeMpi(
                MPI_Allreduce
                ( const_cast<Complex<R>*>(sbuf),
                  rbuf, count, map.type, op, comm )
            );
        }
#else
        MpiMap<Complex<R> > map;
        SafeMpi( 
            MPI_Allreduce
            ( const_cast<Complex<R>*>(sbuf), 
              rbuf, count, map.type, op, comm ) 
        );
#endif
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void AllReduce( const byte* sbuf, byte* rbuf, int count, Op op, Comm comm );
template void AllReduce( const int* sbuf, int* rbuf, int count, Op op, Comm comm );
template void AllReduce( const float* sbuf, float* rbuf, int count, Op op, Comm comm );
template void AllReduce( const double* sbuf, double* rbuf, int count, Op op, Comm comm );
template void AllReduce( const Complex<float>* sbuf, Complex<float>* rbuf, int count, Op op, Comm comm );
template void AllReduce( const Complex<double>* sbuf, Complex<double>* rbuf, int count, Op op, Comm comm );

template<typename R>
void AllReduce( R* buf, int count, Op op, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllReduce");
#endif
    MpiMap<R> map;
    if( count != 0 )
    {
#ifdef HAVE_MPI_IN_PLACE
        SafeMpi( 
            MPI_Allreduce( MPI_IN_PLACE, buf, count, map.type, op, comm )
        );
#else
        std::vector<R> sendBuf( count );
        MemCopy( &sendBuf[0], buf, count );
        SafeMpi(
            MPI_Allreduce( &sendBuf[0], buf, count, map.type, op, comm )
        );
#endif
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void AllReduce( Complex<R>* buf, int count, Op op, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllReduce");
#endif
    if( count != 0 )
    {
#ifdef AVOID_COMPLEX_MPI
        if( op == SUM )
        {
            MpiMap<R> map;
# ifdef HAVE_MPI_IN_PLACE
            SafeMpi(
                MPI_Allreduce( MPI_IN_PLACE, buf, 2*count, map.type, op, comm )
            );
# else
            std::vector<Complex<R> > sendBuf( count );
            MemCopy( &sendBuf[0], buf, count );
            SafeMpi(
                MPI_Allreduce( &sendBuf[0], buf, 2*count, map.type, op, comm )
            );
# endif
        }
        else
        {
            MpiMap<Complex<R> > map;
# ifdef HAVE_MPI_IN_PLACE
            SafeMpi(
                MPI_Allreduce( MPI_IN_PLACE, buf, count, map.type, op, comm )
            );
# else
            std::vector<Complex<R> > sendBuf( count );
            MemCopy( &sendBuf[0], buf, count );
            SafeMpi(
                MPI_Allreduce( &sendBuf[0], buf, count, map.type, op, comm )
            );
# endif
        }
#else
        MpiMap<Complex<R> > map;
# ifdef HAVE_MPI_IN_PLACE
        SafeMpi( 
            MPI_Allreduce( MPI_IN_PLACE, buf, count, map.type, op, comm )
        );
# else
        std::vector<Complex<R> > sendBuf( count );
        MemCopy( &sendBuf[0], buf, count );
        SafeMpi( 
            MPI_Allreduce( &sendBuf[0], buf, count, map.type, op, comm )
        );
# endif
#endif
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void AllReduce( byte* buf, int count, Op op, Comm comm );
template void AllReduce( int* buf, int count, Op op, Comm comm );
template void AllReduce( float* buf, int count, Op op, Comm comm );
template void AllReduce( double* buf, int count, Op op, Comm comm );
template void AllReduce( Complex<float>* buf, int count, Op op, Comm comm );
template void AllReduce( Complex<double>* buf, int count, Op op, Comm comm );

template<typename R>
void ReduceScatter( R* sbuf, R* rbuf, int rc, Op op, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::ReduceScatter");
#endif
#ifdef REDUCE_SCATTER_BLOCK_VIA_ALLREDUCE
    const int commSize = CommSize( comm );
    const int commRank = CommRank( comm );
    AllReduce( sbuf, rc*commSize, op, comm );
    MemCopy( rbuf, &sbuf[commRank*rc], rc );
#elif defined(HAVE_MPI_REDUCE_SCATTER_BLOCK)
    MpiMap<R> map;
    SafeMpi( MPI_Reduce_scatter_block( sbuf, rbuf, rc, map.type, op, comm ) );
#else
    const int commSize = CommSize( comm );
    Reduce( sbuf, rc*commSize, op, 0, comm );
    Scatter( sbuf, rc, rbuf, rc, 0, comm );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void ReduceScatter
( Complex<R>* sbuf, Complex<R>* rbuf, int rc, Op op, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::ReduceScatter");
#endif
#ifdef REDUCE_SCATTER_BLOCK_VIA_ALLREDUCE
    const int commSize = CommSize( comm );
    const int commRank = CommRank( comm );
    AllReduce( sbuf, rc*commSize, op, comm );
    MemCopy( rbuf, &sbuf[commRank*rc], rc );
#elif defined(HAVE_MPI_REDUCE_SCATTER_BLOCK)
# ifdef AVOID_COMPLEX_MPI
    MpiMap<R> map;
    SafeMpi( MPI_Reduce_scatter_block( sbuf, rbuf, 2*rc, map.type, op, comm ) );
# else
    MpiMap<Complex<R> > map;
    SafeMpi( MPI_Reduce_scatter_block( sbuf, rbuf, rc, map.type, op, comm ) );
# endif
#else
    const int commSize = CommSize( comm );
    Reduce( sbuf, rc*commSize, op, 0, comm );
    Scatter( sbuf, rc, rbuf, rc, 0, comm );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

template void ReduceScatter( byte* sbuf, byte* rbuf, int rc, Op op, Comm comm );
template void ReduceScatter( int* sbuf, int* rbuf, int rc, Op op, Comm comm );
template void ReduceScatter( float* sbuf, float* rbuf, int rc, Op op, Comm comm );
template void ReduceScatter( double* sbuf, double* rbuf, int rc, Op op, Comm comm );
template void ReduceScatter( Complex<float>* sbuf, Complex<float>* rbuf, int rc, Op op, Comm comm );
template void ReduceScatter( Complex<double>* sbuf, Complex<double>* rbuf, int rc, Op op, Comm comm );

template<typename R>
void ReduceScatter( R* buf, int rc, Op op, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::ReduceScatter");
#endif
#ifdef REDUCE_SCATTER_BLOCK_VIA_ALLREDUCE
    const int commSize = CommSize( comm );
    const int commRank = CommRank( comm );
    AllReduce( buf, rc*commSize, op, comm );
    if( commRank != 0 )
        MemCopy( buf, &buf[commRank*rc], rc );
#elif defined(HAVE_MPI_REDUCE_SCATTER_BLOCK)
    MpiMap<R> map;
# ifdef HAVE_MPI_IN_PLACE
    SafeMpi( 
        MPI_Reduce_scatter_block( MPI_IN_PLACE, buf, rc, map.type, op, comm )
    );
# else
    const int commSize = CommSize( comm );
    std::vector<R> sendBuf( rc*commSize );
    MemCopy( &sendBuf[0], buf, rc*commSize );
    SafeMpi( 
        MPI_Reduce_scatter_block( &sendBuf[0], buf, rc, map.type, op, comm )
    );
# endif
#else
    const int commSize = CommSize( comm );
    Reduce( buf, rc*commSize, op, 0, comm );
    Scatter( buf, rc, rc, 0, comm );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

// TODO: Handle case where op is not summation
template<typename R>
void ReduceScatter( Complex<R>* buf, int rc, Op op, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::ReduceScatter");
#endif
#ifdef REDUCE_SCATTER_BLOCK_VIA_ALLREDUCE
    const int commSize = CommSize( comm );
    const int commRank = CommRank( comm );
    AllReduce( buf, rc*commSize, op, comm );
    if( commRank != 0 )
        MemCopy( buf, &buf[commRank*rc], rc );
#elif defined(HAVE_MPI_REDUCE_SCATTER_BLOCK)
# ifdef AVOID_COMPLEX_MPI
    MpiMap<R> map;
#  ifdef HAVE_MPI_IN_PLACE
    SafeMpi(
        MPI_Reduce_scatter_block( MPI_IN_PLACE, buf, 2*rc, map.type, op, comm )
    );
#  else 
    const int commSize = CommSize( comm );
    std::vector<Complex<R> > sendBuf( rc*commSize );
    MemCopy( &sendBuf[0], buf, rc*commSize );
    SafeMpi(
        MPI_Reduce_scatter_block( &sendBuf[0], buf, 2*rc, map.type, op, comm )
    );
#  endif
# else
    MpiMap<Complex<R> > map;
#  ifdef HAVE_MPI_IN_PLACE
    SafeMpi( 
        MPI_Reduce_scatter_block( MPI_IN_PLACE, buf, rc, map.type, op, comm )
    );
#  else
    const int commSize = CommSize( comm );
    std::vector<Complex<R> > sendBuf( rc*commSize );
    MemCopy( &sendBuf[0], buf, rc*commSize );
    SafeMpi( 
        MPI_Reduce_scatter_block( &sendBuf[0], buf, rc, map.type, op, comm )
    );
#  endif
# endif
#else
    const int commSize = CommSize( comm );
    Reduce( buf, rc*commSize, op, 0, comm );
    Scatter( buf, rc, rc, 0, comm );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

template void ReduceScatter( byte* buf, int rc, Op op, Comm comm );
template void ReduceScatter( int* buf, int rc, Op op, Comm comm );
template void ReduceScatter( float* buf, int rc, Op op, Comm comm );
template void ReduceScatter( double* buf, int rc, Op op, Comm comm );
template void ReduceScatter( Complex<float>* buf, int rc, Op op, Comm comm );
template void ReduceScatter( Complex<double>* buf, int rc, Op op, Comm comm );

template<typename R>
void ReduceScatter
( const R* sbuf, R* rbuf, const int* rcs, Op op, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::ReduceScatter");
#endif
    MpiMap<R> map;
    SafeMpi( 
        MPI_Reduce_scatter
        ( const_cast<R*>(sbuf), 
          rbuf, const_cast<int*>(rcs), map.type, op, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void ReduceScatter
( const Complex<R>* sbuf, Complex<R>* rbuf, const int* rcs, Op op, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::ReduceScatter");
#endif
#ifdef AVOID_COMPLEX_MPI
    if( op == SUM )
    {
        MpiMap<R> map;
        int p;
        MPI_Comm_size( comm, &p );
        std::vector<int> rcsDoubled(p);
        for( int i=0; i<p; ++i )
            rcsDoubled[i] = 2*rcs[i];
        SafeMpi(
            MPI_Reduce_scatter
            ( const_cast<Complex<R>*>(sbuf),
              rbuf, &rcsDoubled[0], map.type, op, comm )
        );
    }
    else
    {
        MpiMap<Complex<R> > map;
        SafeMpi(
            MPI_Reduce_scatter
            ( const_cast<Complex<R>*>(sbuf),
              rbuf, const_cast<int*>(rcs), map.type, op, comm )
        );
    }
#else
    MpiMap<Complex<R> > map;
    SafeMpi( 
        MPI_Reduce_scatter
        ( const_cast<Complex<R>*>(sbuf), 
          rbuf, const_cast<int*>(rcs), map.type, op, comm ) 
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

template void ReduceScatter( const byte* sbuf, byte* rbuf, const int* rcs, Op op, Comm comm );
template void ReduceScatter( const int* sbuf, int* rbuf, const int* rcs, Op op, Comm comm );
template void ReduceScatter( const float* sbuf, float* rbuf, const int* rcs, Op op, Comm comm );
template void ReduceScatter( const double* sbuf, double* rbuf, const int* rcs, Op op, Comm comm );
template void ReduceScatter( const Complex<float>* sbuf, Complex<float>* rbuf, const int* rcs, Op op, Comm comm );
template void ReduceScatter( const Complex<double>* sbuf, Complex<double>* rbuf, const int* rcs, Op op, Comm comm );

} // namespace mpi
} // namespace elem
