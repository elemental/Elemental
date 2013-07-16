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
//       and a later suggestion by Wolfgang Bangerth on SciComp that the type 
//       be made static const.
template<typename T>
struct MpiMap { static const Datatype type; };
template<>
const Datatype MpiMap<byte>::type              = MPI_UNSIGNED_CHAR;
template<>
const Datatype MpiMap<int>::type               = MPI_INT;
template<>
const Datatype MpiMap<float>::type             = MPI_FLOAT;
template<>
const Datatype MpiMap<double>::type            = MPI_DOUBLE;
template<>
const Datatype MpiMap<Complex<float> >::type   = MPI_COMPLEX;
template<>
const Datatype MpiMap<Complex<double> >::type  = MPI_DOUBLE_COMPLEX;
template<>
const Datatype MpiMap<ValueInt<int> >::type    = MPI_2INT;
template<>
const Datatype MpiMap<ValueInt<float> >::type  = MPI_FLOAT_INT;
template<>
const Datatype MpiMap<ValueInt<double> >::type = MPI_DOUBLE_INT;

//----------------------------//
// MPI environmental routines //
//----------------------------//

void Initialize( int& argc, char**& argv )
{ MPI_Init( &argc, &argv ); }

int InitializeThread( int& argc, char**& argv, int required )
{ 
    int provided; 
#ifdef HAVE_MPI_INIT_THREAD
    MPI_Init_thread( &argc, &argv, required, &provided ); 
#else
    MPI_Init( &argc, &argv );
    provided = 0; // equivalent to MPI_THREAD_SINGLE
#endif
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
#ifdef HAVE_MPI_QUERY_THREAD
    MPI_Query_thread( &provided );
#else
    provided = 0; // equivalent to MPI_THREAD_SINGLE
#endif
    return provided;
}

double Time()
{ return MPI_Wtime(); }

void OpCreate( UserFunction* func, bool commutes, Op& op )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::OpCreate");
#endif
    SafeMpi( MPI_Op_create( func, commutes, &op ) );
}

void OpFree( Op& op )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::OpFree");
#endif
    SafeMpi( MPI_Op_free( &op ) );
}

//---------------------------//
// Communicator manipulation //
//---------------------------//

int WorldRank()
{
#ifndef RELEASE
    CallStackEntry entry("mpi::WorldRank");
#endif
    int rank;
    SafeMpi( MPI_Comm_rank( COMM_WORLD, &rank ) );
    return rank;
}

int CommRank( Comm comm )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::CommRank");
#endif
    int rank;
    SafeMpi( MPI_Comm_rank( comm, &rank ) );
    return rank;
}

int CommSize( Comm comm )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::CommSize");
#endif
    int size;
    SafeMpi( MPI_Comm_size( comm, &size ) );
    return size;
}

void CommCreate( Comm parentComm, Group subsetGroup, Comm& subsetComm )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::CommCreate");
#endif
    SafeMpi( MPI_Comm_create( parentComm, subsetGroup, &subsetComm ) );
}

void CommDup( Comm original, Comm& duplicate )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::CommDup");
#endif
    SafeMpi( MPI_Comm_dup( original, &duplicate ) );
}

void CommSplit( Comm comm, int color, int key, Comm& newComm )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::CommSplit");
#endif
    SafeMpi( MPI_Comm_split( comm, color, key, &newComm ) );
}

void CommFree( Comm& comm )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::CommFree");
#endif
    SafeMpi( MPI_Comm_free( &comm ) );
}

bool CongruentComms( Comm comm1, Comm comm2 )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::CongruentComms");
#endif
    int result;
    SafeMpi( MPI_Comm_compare( comm1, comm2, &result ) );
    return ( result == MPI_IDENT || result == MPI_CONGRUENT );
}

void ErrorHandlerSet( Comm comm, ErrorHandler errorHandler )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::ErrorHandlerSet");
#endif
#ifdef HAVE_MPI_COMM_SET_ERRHANDLER
    SafeMpi( MPI_Comm_set_errhandler( comm, errorHandler ) );
#else
    SafeMpi( MPI_Errhandler_set( comm, errorHandler ) );
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
    CallStackEntry entry("mpi::CartCreate");
#endif
    SafeMpi
    ( MPI_Cart_create
      ( comm, numDims, const_cast<int*>(dimensions), 
        const_cast<int*>(periods), reorder, &cartComm ) );
}

void CartSub( Comm comm, const int* remainingDims, Comm& subComm )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::CartSub");
#endif
    SafeMpi( MPI_Cart_sub( comm, const_cast<int*>(remainingDims), &subComm ) );
}

//--------------------//
// Group manipulation //
//--------------------//

int GroupRank( Group group )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::GroupRank");
#endif
    int rank;
    SafeMpi( MPI_Group_rank( group, &rank ) );
    return rank;
}

int GroupSize( Group group )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::GroupSize");
#endif
    int size;
    SafeMpi( MPI_Group_size( group, &size ) );
    return size;
}

void CommGroup( Comm comm, Group& group )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::CommGroup");
#endif
    SafeMpi( MPI_Comm_group( comm, &group ) );
}

void GroupIncl( Group group, int n, const int* ranks, Group& subGroup )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::GroupIncl");
#endif
    SafeMpi( MPI_Group_incl( group, n, const_cast<int*>(ranks), &subGroup ) );
}

void GroupDifference( Group parent, Group subset, Group& complement )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::GroupDifference");
#endif
    SafeMpi( MPI_Group_difference( parent, subset, &complement ) );
}

void GroupFree( Group& group )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::GroupFree");
#endif
    SafeMpi( MPI_Group_free( &group ) );
}

void GroupTranslateRanks
( Group origGroup, int size, const int* origRanks, 
  Group newGroup,                  int* newRanks )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::GroupTranslateRanks");
#endif
    SafeMpi
    ( MPI_Group_translate_ranks
      ( origGroup, size, const_cast<int*>(origRanks), newGroup, newRanks ) );
}

// Wait until every process in comm reaches this statement
void Barrier( Comm comm )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::Barrier");
#endif
    SafeMpi( MPI_Barrier( comm ) );
}

// Test for completion
bool Test( Request& request )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::Test");
#endif
    Status status;
    int flag;
    SafeMpi( MPI_Test( &request, &flag, &status ) );
    return flag;
}

// Ensure that the request finishes before continuing
void Wait( Request& request )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::Wait");
#endif
    Status status;
    SafeMpi( MPI_Wait( &request, &status ) );
}

// Ensure that the request finishes before continuing
void Wait( Request& request, Status& status )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::Wait");
#endif
    SafeMpi( MPI_Wait( &request, &status ) );
}

// Ensure that several requests finish before continuing
void WaitAll( int numRequests, Request* requests )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::WaitAll");
#endif
    std::vector<Status> statuses( numRequests );
    SafeMpi( MPI_Waitall( numRequests, requests, &statuses[0] ) );
}

// Ensure that several requests finish before continuing
void WaitAll( int numRequests, Request* requests, Status* statuses )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::WaitAll");
#endif
    SafeMpi( MPI_Waitall( numRequests, requests, statuses ) );
}

// Nonblocking test for message completion
bool IProbe( int source, int tag, Comm comm, Status& status )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::IProbe");
#endif
    int flag;
    SafeMpi( MPI_Iprobe( source, tag, comm, &flag, &status ) );
    return flag;
}

template<typename T>
int GetCount( Status& status )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::GetCount");
#endif
    int count;
    SafeMpi( MPI_Get_count( &status, MpiMap<T>::type, &count ) );
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
    CallStackEntry entry("mpi::Send");
#endif
    SafeMpi
    ( MPI_Send( const_cast<R*>(buf), count, MpiMap<R>::type, to, tag, comm ) );
}

template<typename R>
void Send( const Complex<R>* buf, int count, int to, int tag, Comm comm )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::Send");
#endif
#ifdef AVOID_COMPLEX_MPI
    SafeMpi
    ( MPI_Send
      ( const_cast<Complex<R>*>(buf), 2*count, MpiMap<R>::type, to, 
        tag, comm ) );
#else
    SafeMpi
    ( MPI_Send
      ( const_cast<Complex<R>*>(buf), count, 
        MpiMap<Complex<R> >::type, to, tag, comm ) );
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
    CallStackEntry entry("mpi::ISend");
#endif
    SafeMpi
    ( MPI_Isend
      ( const_cast<R*>(buf), count, MpiMap<R>::type, to, 
        tag, comm, &request ) );
}

template<typename R>
void ISend
( const Complex<R>* buf, int count, int to, int tag, Comm comm, 
  Request& request )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::ISend");
#endif
#ifdef AVOID_COMPLEX_MPI
    SafeMpi
    ( MPI_Isend
      ( const_cast<Complex<R>*>(buf), 2*count, MpiMap<R>::type, to, tag, comm,
        &request ) );
#else
    SafeMpi
    ( MPI_Isend
      ( const_cast<Complex<R>*>(buf), count, 
        MpiMap<Complex<R> >::type, to, tag, comm, &request ) );
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
    CallStackEntry entry("mpi::ISSend");
#endif
    SafeMpi
    ( MPI_Issend
      ( const_cast<R*>(buf), count, MpiMap<R>::type, to, 
        tag, comm, &request ) );
}

template<typename R>
void ISSend
( const Complex<R>* buf, int count, int to, int tag, Comm comm, 
  Request& request )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::ISSend");
#endif
#ifdef AVOID_COMPLEX_MPI
    SafeMpi
    ( MPI_Issend
      ( const_cast<Complex<R>*>(buf), 2*count, MpiMap<R>::type, to, tag, comm,
        &request ) );
#else
    SafeMpi
    ( MPI_Issend
      ( const_cast<Complex<R>*>(buf), count, 
        MpiMap<Complex<R> >::type, to, tag, comm, &request ) );
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
    CallStackEntry entry("mpi::Recv");
#endif
    Status status;
    SafeMpi
    ( MPI_Recv( buf, count, MpiMap<R>::type, from, tag, comm, &status ) );
}

template<typename R>
void Recv( Complex<R>* buf, int count, int from, int tag, Comm comm )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::Recv");
#endif
    Status status;
#ifdef AVOID_COMPLEX_MPI
    SafeMpi
    ( MPI_Recv( buf, 2*count, MpiMap<R>::type, from, tag, comm, &status ) );
#else
    SafeMpi
    ( MPI_Recv
      ( buf, count, MpiMap<Complex<R> >::type, from, tag, comm, &status ) );
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
    CallStackEntry entry("mpi::IRecv");
#endif
    SafeMpi
    ( MPI_Irecv( buf, count, MpiMap<R>::type, from, tag, comm, &request ) );
}

template<typename R>
void IRecv
( Complex<R>* buf, int count, int from, int tag, Comm comm, Request& request )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::IRecv");
#endif
#ifdef AVOID_COMPLEX_MPI
    SafeMpi
    ( MPI_Irecv( buf, 2*count, MpiMap<R>::type, from, tag, comm, &request ) );
#else
    SafeMpi
    ( MPI_Irecv
      ( buf, count, MpiMap<Complex<R> >::type, from, tag, comm, &request ) );
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
    CallStackEntry entry("mpi::SendRecv");
#endif
    Status status;
    SafeMpi
    ( MPI_Sendrecv
      ( const_cast<R*>(sbuf), sc, MpiMap<R>::type, to,   stag,
        rbuf,                 rc, MpiMap<R>::type, from, rtag, 
        comm, &status ) );
}

template<typename R>
void SendRecv
( const Complex<R>* sbuf, int sc, int to,   int stag,
        Complex<R>* rbuf, int rc, int from, int rtag, Comm comm )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::SendRecv");
#endif
    Status status;
#ifdef AVOID_COMPLEX_MPI
    SafeMpi
    ( MPI_Sendrecv
      ( const_cast<Complex<R>*>(sbuf), 2*sc, MpiMap<R>::type, to,   stag,
        rbuf,                          2*rc, MpiMap<R>::type, from, rtag, 
        comm, &status ) );
#else
    SafeMpi
    ( MPI_Sendrecv
      ( const_cast<Complex<R>*>(sbuf), 
        sc, MpiMap<Complex<R> >::type, to,   stag,
        rbuf,                          
        rc, MpiMap<Complex<R> >::type, from, rtag, comm, &status ) );
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
void SendRecv
( R* buf, int count, int to, int stag, int from, int rtag, Comm comm )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::SendRecv");
#endif
    Status status;
    SafeMpi
    ( MPI_Sendrecv_replace
      ( buf, count, MpiMap<R>::type, to, stag, from, rtag, comm, &status ) );
}

template<typename R>
void SendRecv
( Complex<R>* buf, int count, int to, int stag, int from, int rtag, Comm comm )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::SendRecv");
#endif
    Status status;
#ifdef AVOID_COMPLEX_MPI
    SafeMpi
    ( MPI_Sendrecv_replace
      ( buf, 2*count, MpiMap<R>::type, to, stag, from, rtag, comm, &status ) );
#else
    SafeMpi
    ( MPI_Sendrecv_replace
      ( buf, count, MpiMap<Complex<R> >::type, 
        to, stag, from, rtag, comm, &status ) );
#endif
}

template void SendRecv
( byte* buf, int count, int to, int stag, int from, int rtag, Comm comm );
template void SendRecv
( int* buf, int count, int to, int stag, int from, int rtag, Comm comm );
template void SendRecv
( float* buf, int count, int to, int stag, int from, int rtag, Comm comm );
template void SendRecv
( double* buf, int count, int to, int stag, int from, int rtag, Comm comm );
template void SendRecv
( Complex<float>* buf, int count, int to, int stag, 
  int from, int rtag, Comm comm );
template void SendRecv
( Complex<double>* buf, int sc, int to, int stag, 
  int from, int rtag, Comm comm );

template<typename R>
void Broadcast( R* buf, int count, int root, Comm comm )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::Broadcast");
#endif
    SafeMpi( MPI_Bcast( buf, count, MpiMap<R>::type, root, comm ) );
}

template<typename R>
void Broadcast( Complex<R>* buf, int count, int root, Comm comm )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::Broadcast");
#endif
#ifdef AVOID_COMPLEX_MPI
    SafeMpi( MPI_Bcast( buf, 2*count, MpiMap<R>::type, root, comm ) );
#else
    SafeMpi( MPI_Bcast( buf, count, MpiMap<Complex<R> >::type, root, comm ) );
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
    CallStackEntry entry("mpi::IBroadcast");
#endif
    SafeMpi( MPI_Ibcast( buf, count, MpiMap<R>::type, root, comm, &request ) );
}

template<typename R>
void IBroadcast
( Complex<R>* buf, int count, int root, Comm comm, Request& request )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::IBroadcast");
#endif
#ifdef AVOID_COMPLEX_MPI
    SafeMpi
    ( MPI_Ibcast( buf, 2*count, MpiMap<R>::type, root, comm, &request ) );
#else
    SafeMpi
    ( MPI_Ibcast
      ( buf, count, MpiMap<Complex<R> >::type, root, comm, &request ) );
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
    CallStackEntry entry("mpi::Gather");
#endif
    SafeMpi
    ( MPI_Gather
      ( const_cast<R*>(sbuf), sc, MpiMap<R>::type,
        rbuf,                 rc, MpiMap<R>::type, root, comm ) );
}

template<typename R>
void Gather
( const Complex<R>* sbuf, int sc,
        Complex<R>* rbuf, int rc, int root, Comm comm )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::Gather");
#endif
#ifdef AVOID_COMPLEX_MPI
    SafeMpi
    ( MPI_Gather
      ( const_cast<Complex<R>*>(sbuf), 2*sc, MpiMap<R>::type,
        rbuf,                          2*rc, MpiMap<R>::type, root, comm ) );
#else
    SafeMpi
    ( MPI_Gather
      ( const_cast<Complex<R>*>(sbuf), sc, MpiMap<Complex<R> >::type,
        rbuf,                          rc, MpiMap<Complex<R> >::type, 
        root, comm ) );
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
    CallStackEntry entry("mpi::IGather");
#endif
    SafeMpi
    ( MPI_Igather
      ( const_cast<R*>(sbuf), sc, MpiMap<R>::type,
        rbuf,                 rc, MpiMap<R>::type, root, comm, &request ) );
}

template<typename R>
void IGather
( const Complex<R>* sbuf, int sc,
        Complex<R>* rbuf, int rc, int root, Comm comm, Request& request )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::IGather");
#endif
#ifdef AVOID_COMPLEX_MPI
    SafeMpi
    ( MPI_Igather
      ( const_cast<Complex<R>*>(sbuf), 2*sc, MpiMap<R>::type,
        rbuf,                          2*rc, MpiMap<R>::type, 
        root, comm, &request ) );
#else
    SafeMpi
    ( MPI_Igather
      ( const_cast<Complex<R>*>(sbuf), sc, MpiMap<Complex<R> >::type,
        rbuf,                          rc, MpiMap<Complex<R> >::type, 
        root, comm, &request ) );
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
    CallStackEntry entry("mpi::Gather");
#endif
    SafeMpi
    ( MPI_Gatherv
      ( const_cast<R*>(sbuf), 
        sc,       
        MpiMap<R>::type,
        rbuf,                    
        const_cast<int*>(rcs), 
        const_cast<int*>(rds), 
        MpiMap<R>::type,
        root, 
        comm ) );
}

template<typename R>
void Gather
( const Complex<R>* sbuf, int sc,
        Complex<R>* rbuf, const int* rcs, const int* rds, int root, Comm comm )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::Gather");
#endif
#ifdef AVOID_COMPLEX_MPI
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
    SafeMpi
    ( MPI_Gatherv
      ( const_cast<Complex<R>*>(sbuf), 2*sc, MpiMap<R>::type,
        rbuf, &rcsDouble[0], &rdsDouble[0], MpiMap<R>::type,
        root, comm ) );
#else
    SafeMpi
    ( MPI_Gatherv
      ( const_cast<Complex<R>*>(sbuf), 
        sc,       
        MpiMap<Complex<R> >::type,
        rbuf,  
        const_cast<int*>(rcs), 
        const_cast<int*>(rds), 
        MpiMap<Complex<R> >::type,
        root, 
        comm ) );
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
        Complex<float>* rbuf, const int* rcs, const int* rds, 
  int root, Comm comm );
template void Gather
( const Complex<double>* sbuf, int sc, 
        Complex<double>* rbuf, const int* rcs, const int* rds, 
  int root, Comm comm );

template<typename R>
void AllGather
( const R* sbuf, int sc,
        R* rbuf, int rc, Comm comm )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::AllGather");
#endif
#ifdef USE_BYTE_ALLGATHERS
    SafeMpi
    ( MPI_Allgather
      ( const_cast<R*>(sbuf), sizeof(R)*sc, MPI_UNSIGNED_CHAR, 
        rbuf,                 sizeof(R)*rc, MPI_UNSIGNED_CHAR, comm ) );
#else
    SafeMpi
    ( MPI_Allgather
      ( const_cast<R*>(sbuf), sc, MpiMap<R>::type, 
        rbuf,                 rc, MpiMap<R>::type, comm ) );
#endif
}

template<typename R>
void AllGather
( const Complex<R>* sbuf, int sc,
        Complex<R>* rbuf, int rc, Comm comm )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::AllGather");
#endif
#ifdef USE_BYTE_ALLGATHERS
    SafeMpi
    ( MPI_Allgather
      ( const_cast<Complex<R>*>(sbuf), 2*sizeof(R)*sc, MPI_UNSIGNED_CHAR, 
        rbuf,                          2*sizeof(R)*rc, MPI_UNSIGNED_CHAR, 
        comm ) );
#else
 #ifdef AVOID_COMPLEX_MPI
    SafeMpi
    ( MPI_Allgather
      ( const_cast<Complex<R>*>(sbuf), 2*sc, MpiMap<R>::type,
        rbuf,                          2*rc, MpiMap<R>::type, comm ) );
 #else
    SafeMpi
    ( MPI_Allgather
      ( const_cast<Complex<R>*>(sbuf), sc, MpiMap<Complex<R> >::type,
        rbuf,                          rc, MpiMap<Complex<R> >::type, comm ) );
 #endif
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
    CallStackEntry entry("mpi::AllGather");
#endif
#ifdef USE_BYTE_ALLGATHERS
    const int commSize = CommSize( comm );
    std::vector<int> byteRcs( commSize ), byteRds( commSize );
    for( int i=0; i<commSize; ++i )
    {
        byteRcs[i] = sizeof(R)*rcs[i];
        byteRds[i] = sizeof(R)*rds[i];
    }
    SafeMpi
    ( MPI_Allgatherv
      ( const_cast<R*>(sbuf), sizeof(R)*sc, MPI_UNSIGNED_CHAR, 
        rbuf, &byteRcs[0], &byteRds[0],     MPI_UNSIGNED_CHAR, comm ) );
#else
    SafeMpi
    ( MPI_Allgatherv
      ( const_cast<R*>(sbuf), 
        sc, 
        MpiMap<R>::type, 
        rbuf,   
        const_cast<int*>(rcs), 
        const_cast<int*>(rds), 
        MpiMap<R>::type, 
        comm ) );
#endif
}

template<typename R>
void AllGather
( const Complex<R>* sbuf, int sc,
        Complex<R>* rbuf, const int* rcs, const int* rds, Comm comm )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::AllGather");
#endif
#ifdef USE_BYTE_ALLGATHERS
    const int commSize = CommSize( comm );
    std::vector<int> byteRcs( commSize ), byteRds( commSize );
    for( int i=0; i<commSize; ++i )
    {
        byteRcs[i] = 2*sizeof(R)*rcs[i];
        byteRds[i] = 2*sizeof(R)*rds[i];
    }
    SafeMpi
    ( MPI_Allgatherv
      ( const_cast<Complex<R>*>(sbuf), 2*sizeof(R)*sc, MPI_UNSIGNED_CHAR, 
        rbuf, &byteRcs[0], &byteRds[0], MPI_UNSIGNED_CHAR, comm ) );
#else
 #ifdef AVOID_COMPLEX_MPI
    const int commSize = CommSize( comm );
    std::vector<int> realRcs( commSize ), realRds( commSize );
    for( int i=0; i<commSize; ++i )
    {
        realRcs[i] = 2*rcs[i];
        realRds[i] = 2*rds[i];
    }
    SafeMpi
    ( MPI_Allgatherv
      ( const_cast<Complex<R>*>(sbuf), 2*sc, MpiMap<R>::type,
        rbuf, &realRcs[0], &realRds[0],      MpiMap<R>::type, comm ) );
 #else
    SafeMpi
    ( MPI_Allgatherv
      ( const_cast<Complex<R>*>(sbuf), 
        sc, 
        MpiMap<Complex<R> >::type,
        rbuf, 
        const_cast<int*>(rcs), 
        const_cast<int*>(rds), 
        MpiMap<Complex<R> >::type,
        comm ) );
 #endif
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
    CallStackEntry entry("mpi::Scatter");
#endif
    SafeMpi
    ( MPI_Scatter
      ( const_cast<R*>(sbuf), sc, MpiMap<R>::type,
        rbuf,                 rc, MpiMap<R>::type, root, comm ) );
}

template<typename R>
void Scatter
( const Complex<R>* sbuf, int sc,
        Complex<R>* rbuf, int rc, int root, Comm comm )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::Scatter");
#endif
#ifdef AVOID_COMPLEX_MPI
    SafeMpi
    ( MPI_Scatter
      ( const_cast<Complex<R>*>(sbuf), 2*sc, MpiMap<R>::type,
        rbuf,                          2*rc, MpiMap<R>::type, root, comm ) );
#else
    SafeMpi
    ( MPI_Scatter
      ( const_cast<Complex<R>*>(sbuf), sc, MpiMap<Complex<R> >::type,
        rbuf,                          rc, MpiMap<Complex<R> >::type, 
        root, comm ) );
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
    CallStackEntry entry("mpi::Scatter");
#endif
    const int commRank = CommRank( comm );
    if( commRank == root )
    {
#ifdef HAVE_MPI_IN_PLACE
        SafeMpi
        ( MPI_Scatter
          ( buf,          sc, MpiMap<R>::type, 
            MPI_IN_PLACE, rc, MpiMap<R>::type, root, comm ) );
#else
        const int commSize = CommSize( comm );
        std::vector<R> sendBuf( sc*commSize );
        MemCopy( &sendBuf[0], buf, sc*commSize );
        SafeMpi
        ( MPI_Scatter
          ( &sendBuf[0], sc, MpiMap<R>::type, 
            buf,         rc, MpiMap<R>::type, root, comm ) );
#endif
    }
    else
    {
        SafeMpi
        ( MPI_Scatter
          ( 0,   sc, MpiMap<R>::type, 
            buf, rc, MpiMap<R>::type, root, comm ) );
    }
}

template<typename R>
void Scatter( Complex<R>* buf, int sc, int rc, int root, Comm comm )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::Scatter");
#endif
    const int commRank = CommRank( comm );
    if( commRank == root )
    {
#ifdef AVOID_COMPLEX_MPI
# ifdef HAVE_MPI_IN_PLACE
        SafeMpi
        ( MPI_Scatter
          ( buf,          2*sc, MpiMap<R>::type, 
            MPI_IN_PLACE, 2*rc, MpiMap<R>::type, root, comm ) );
# else
        const int commSize = CommSize( comm );
        std::vector<Complex<R> > sendBuf( sc*commSize );
        MemCopy( &sendBuf[0], buf, sc*commSize );
        SafeMpi
        ( MPI_Scatter
          ( &sendBuf[0], 2*sc, MpiMap<R>::type,          
            buf,         2*rc, MpiMap<R>::type, root, comm ) );
# endif
#else
# ifdef HAVE_MPI_IN_PLACE
        SafeMpi
        ( MPI_Scatter
          ( buf,          sc, MpiMap<Complex<R> >::type, 
            MPI_IN_PLACE, rc, MpiMap<Complex<R> >::type, root, comm ) );
# else
        const int commSize = CommSize( comm );
        std::vector<Complex<R> > sendBuf( sc*commSize );
        MemCopy( &sendBuf[0], buf, sc*commSize );
        SafeMpi
        ( MPI_Scatter
          ( &sendBuf[0], sc, MpiMap<Complex<R> >::type,
            buf,         rc, MpiMap<Complex<R> >::type, root, comm ) );
# endif
#endif
    }
    else
    {
#ifdef AVOID_COMPLEX_MPI
        SafeMpi
        ( MPI_Scatter
          ( 0,   2*sc, MpiMap<R>::type, 
            buf, 2*rc, MpiMap<R>::type, root, comm ) );
#else
        SafeMpi
        ( MPI_Scatter
          ( 0,   sc, MpiMap<Complex<R> >::type, 
            buf, rc, MpiMap<Complex<R> >::type, root, comm ) );
#endif
    }
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
    CallStackEntry entry("mpi::AllToAll");
#endif
    SafeMpi
    ( MPI_Alltoall
      ( const_cast<R*>(sbuf), sc, MpiMap<R>::type,
        rbuf,                 rc, MpiMap<R>::type, comm ) );
}

template<typename R>
void AllToAll
( const Complex<R>* sbuf, int sc,
        Complex<R>* rbuf, int rc, Comm comm )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::AllToAll");
#endif
#ifdef AVOID_COMPLEX_MPI
    SafeMpi
    ( MPI_Alltoall
      ( const_cast<Complex<R>*>(sbuf), 2*sc, MpiMap<R>::type,
        rbuf,                          2*rc, MpiMap<R>::type, comm ) );
#else
    SafeMpi
    ( MPI_Alltoall
      ( const_cast<Complex<R>*>(sbuf), sc, MpiMap<Complex<R> >::type,
        rbuf,                          rc, MpiMap<Complex<R> >::type, comm ) );
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
    CallStackEntry entry("mpi::AllToAll");
#endif
    SafeMpi
    ( MPI_Alltoallv
      ( const_cast<R*>(sbuf), 
        const_cast<int*>(scs), 
        const_cast<int*>(sds), 
        MpiMap<R>::type,
        rbuf, 
        const_cast<int*>(rcs), 
        const_cast<int*>(rds), 
        MpiMap<R>::type,
        comm ) );
}

template<typename R>
void AllToAll
( const Complex<R>* sbuf, const int* scs, const int* sds,
        Complex<R>* rbuf, const int* rcs, const int* rds, Comm comm )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::AllToAll");
#endif
#ifdef AVOID_COMPLEX_MPI
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
    SafeMpi
    ( MPI_Alltoallv
      ( const_cast<Complex<R>*>(sbuf),
              &scsDoubled[0], &sdsDoubled[0], MpiMap<R>::type,
        rbuf, &rcsDoubled[0], &rdsDoubled[0], MpiMap<R>::type, comm ) );
#else
    SafeMpi
    ( MPI_Alltoallv
      ( const_cast<Complex<R>*>(sbuf), 
        const_cast<int*>(scs), 
        const_cast<int*>(sds), 
        MpiMap<Complex<R> >::type,
        rbuf, 
        const_cast<int*>(rcs), 
        const_cast<int*>(rds), 
        MpiMap<Complex<R> >::type,
        comm ) );
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

template<typename T>
void Reduce
( const T* sbuf, T* rbuf, int count, Op op, int root, Comm comm )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::Reduce");
#endif
    if( count != 0 )
    {
        SafeMpi
        ( MPI_Reduce
          ( const_cast<T*>(sbuf), rbuf, count, MpiMap<T>::type, op, 
            root, comm ) );
    }
}

template<typename R>
void Reduce
( const Complex<R>* sbuf, 
        Complex<R>* rbuf, int count, Op op, int root, Comm comm )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::Reduce");
#endif
    if( count != 0 )
    {
#ifdef AVOID_COMPLEX_MPI
        if( op == SUM )
        {
            SafeMpi
            ( MPI_Reduce
              ( const_cast<Complex<R>*>(sbuf),
                rbuf, 2*count, MpiMap<R>::type, op, root, comm ) );
        }
        else
        {
            SafeMpi
            ( MPI_Reduce
              ( const_cast<Complex<R>*>(sbuf),
                rbuf, count, MpiMap<Complex<R> >::type, op, root, comm ) );
        }
#else
        SafeMpi
        ( MPI_Reduce
          ( const_cast<Complex<R>*>(sbuf), 
            rbuf, count, MpiMap<Complex<R> >::type, op, root, comm ) );
#endif
    }
}

template void Reduce( const byte* sbuf, byte* rbuf, int count, Op op, int root, Comm comm );
template void Reduce( const int* sbuf, int* rbuf, int count, Op op, int root, Comm comm );
template void Reduce( const float* sbuf, float* rbuf, int count, Op op, int root, Comm comm );
template void Reduce( const double* sbuf, double* rbuf, int count, Op op, int root, Comm comm );
template void Reduce( const Complex<float>* sbuf, Complex<float>* rbuf, int count, Op op, int root, Comm comm );
template void Reduce( const Complex<double>* sbuf, Complex<double>* rbuf, int count, Op op, int root, Comm comm );
template void Reduce( const ValueInt<int>* sbuf, ValueInt<int>* rbuf, int count, Op op, int root, Comm comm );
template void Reduce( const ValueInt<float>* sbuf, ValueInt<float>* rbuf, int count, Op op, int root, Comm comm );
template void Reduce( const ValueInt<double>* sbuf, ValueInt<double>* rbuf, int count, Op op, int root, Comm comm );

template<typename T>
void Reduce( T* buf, int count, Op op, int root, Comm comm )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::Reduce");
#endif
    if( count != 0 )
    {
        const int commRank = CommRank( comm );
        if( commRank == root )
        {
#ifdef HAVE_MPI_IN_PLACE
            SafeMpi
            ( MPI_Reduce
              ( MPI_IN_PLACE, buf, count, MpiMap<T>::type, op, root, comm ) );
#else
            std::vector<T> sendBuf( count );
            MemCopy( &sendBuf[0], buf, count );
            SafeMpi
            ( MPI_Reduce
              ( &sendBuf[0], buf, count, MpiMap<T>::type, op, root, comm ) );
#endif
        }
        else
            SafeMpi
            ( MPI_Reduce( buf, 0, count, MpiMap<T>::type, op, root, comm ) );
    }
}

template<typename R>
void Reduce( Complex<R>* buf, int count, Op op, int root, Comm comm )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::Reduce");
#endif
    if( count != 0 )
    {
        const int commRank = CommRank( comm );
#ifdef AVOID_COMPLEX_MPI
        if( op == SUM )
        {
            if( commRank == root )
            {
# ifdef HAVE_MPI_IN_PLACE
                SafeMpi
                ( MPI_Reduce
                  ( MPI_IN_PLACE, buf, 2*count, MpiMap<R>::type, op, 
                    root, comm ) );
# else
                std::vector<Complex<R> > sendBuf( count );
                MemCopy( &sendBuf[0], buf, count );
                SafeMpi
                ( MPI_Reduce
                  ( &sendBuf[0], buf, 2*count, MpiMap<R>::type, op, 
                    root, comm ) );
# endif
            }
            else
                SafeMpi
                ( MPI_Reduce
                  ( buf, 0, 2*count, MpiMap<R>::type, op, root, comm ) );
        }
        else
        {
            if( commRank == root )
            {
# ifdef HAVE_MPI_IN_PLACE
                SafeMpi
                ( MPI_Reduce
                  ( MPI_IN_PLACE, buf, count, MpiMap<Complex<R> >::type, op, 
                    root, comm ) );
# else
                std::vector<Complex<R> > sendBuf( count );
                MemCopy( &sendBuf[0], buf, count );
                SafeMpi
                ( MPI_Reduce
                  ( &sendBuf[0], buf, count, MpiMap<Complex<R> >::type, op, 
                    root, comm ) );
# endif
            }
            else
                SafeMpi
                ( MPI_Reduce
                  ( buf, 0, count, MpiMap<Complex<R> >::type, op, 
                    root, comm ) );
        }
#else
        if( commRank == root )
        {
# ifdef HAVE_MPI_IN_PLACE
            SafeMpi
            ( MPI_Reduce
              ( MPI_IN_PLACE, buf, count, MpiMap<Complex<R> >::type, op, 
                root, comm ) );
# else
            std::vector<Complex<R> > sendBuf( count );
            MemCopy( &sendBuf[0], buf, count );
            SafeMpi
            ( MPI_Reduce
              ( &sendBuf[0], buf, count, MpiMap<Complex<R> >::type, op, 
                root, comm ) );
# endif
        }
        else
            SafeMpi
            ( MPI_Reduce
              ( buf, 0, count, MpiMap<Complex<R> >::type, op, root, comm ) );
#endif
    }
}

template void Reduce( byte* buf, int count, Op op, int root, Comm comm );
template void Reduce( int* buf, int count, Op op, int root, Comm comm );
template void Reduce( float* buf, int count, Op op, int root, Comm comm );
template void Reduce( double* buf, int count, Op op, int root, Comm comm );
template void Reduce( Complex<float>* buf, int count, Op op, int root, Comm comm );
template void Reduce( Complex<double>* buf, int count, Op op, int root, Comm comm );
template void Reduce( ValueInt<int>* buf, int count, Op op, int root, Comm comm );
template void Reduce( ValueInt<float>* buf, int count, Op op, int root, Comm comm );
template void Reduce( ValueInt<double>* buf, int count, Op op, int root, Comm comm );

template<typename T>
void AllReduce( const T* sbuf, T* rbuf, int count, Op op, Comm comm )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::AllReduce");
#endif
    if( count != 0 )
    {
        SafeMpi
        ( MPI_Allreduce
          ( const_cast<T*>(sbuf), rbuf, count, MpiMap<T>::type, op, comm ) );
    }
}

template<typename R>
void AllReduce
( const Complex<R>* sbuf, Complex<R>* rbuf, int count, Op op, Comm comm )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::AllReduce");
#endif
    if( count != 0 )
    {
#ifdef AVOID_COMPLEX_MPI
        if( op == SUM )
        {
            SafeMpi
            ( MPI_Allreduce
                ( const_cast<Complex<R>*>(sbuf),
                  rbuf, 2*count, MpiMap<R>::type, op, comm ) );
        }
        else
        {
            SafeMpi
            ( MPI_Allreduce
              ( const_cast<Complex<R>*>(sbuf),
                rbuf, count, MpiMap<Complex<R> >::type, op, comm ) );
        }
#else
        SafeMpi
        ( MPI_Allreduce
          ( const_cast<Complex<R>*>(sbuf), 
            rbuf, count, MpiMap<Complex<R> >::type, op, comm ) );
#endif
    }
}

template void AllReduce( const byte* sbuf, byte* rbuf, int count, Op op, Comm comm );
template void AllReduce( const int* sbuf, int* rbuf, int count, Op op, Comm comm );
template void AllReduce( const float* sbuf, float* rbuf, int count, Op op, Comm comm );
template void AllReduce( const double* sbuf, double* rbuf, int count, Op op, Comm comm );
template void AllReduce( const Complex<float>* sbuf, Complex<float>* rbuf, int count, Op op, Comm comm );
template void AllReduce( const Complex<double>* sbuf, Complex<double>* rbuf, int count, Op op, Comm comm );
template void AllReduce( const ValueInt<int>* sbuf, ValueInt<int>* rbuf, int count, Op op, Comm comm );
template void AllReduce( const ValueInt<float>* sbuf, ValueInt<float>* rbuf, int count, Op op, Comm comm );
template void AllReduce( const ValueInt<double>* sbuf, ValueInt<double>* rbuf, int count, Op op, Comm comm );

template<typename T>
void AllReduce( T* buf, int count, Op op, Comm comm )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::AllReduce");
#endif
    if( count != 0 )
    {
#ifdef HAVE_MPI_IN_PLACE
        SafeMpi
        ( MPI_Allreduce
          ( MPI_IN_PLACE, buf, count, MpiMap<T>::type, op, comm ) );
#else
        std::vector<T> sendBuf( count );
        MemCopy( &sendBuf[0], buf, count );
        SafeMpi
        ( MPI_Allreduce( &sendBuf[0], buf, count, MpiMap<T>::type, op, comm ) );
#endif
    }
}

template<typename R>
void AllReduce( Complex<R>* buf, int count, Op op, Comm comm )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::AllReduce");
#endif
    if( count != 0 )
    {
#ifdef AVOID_COMPLEX_MPI
        if( op == SUM )
        {
# ifdef HAVE_MPI_IN_PLACE
            SafeMpi
            ( MPI_Allreduce
              ( MPI_IN_PLACE, buf, 2*count, MpiMap<R>::type, op, comm ) );
# else
            std::vector<Complex<R> > sendBuf( count );
            MemCopy( &sendBuf[0], buf, count );
            SafeMpi
            ( MPI_Allreduce
              ( &sendBuf[0], buf, 2*count, MpiMap<R>::type, op, comm ) );
# endif
        }
        else
        {
            MpiMap<Complex<R> > map;
# ifdef HAVE_MPI_IN_PLACE
            SafeMpi
            ( MPI_Allreduce
              ( MPI_IN_PLACE, buf, count, MpiMap<Complex<R> >::type, 
                op, comm ) );
# else
            std::vector<Complex<R> > sendBuf( count );
            MemCopy( &sendBuf[0], buf, count );
            SafeMpi
            ( MPI_Allreduce
              ( &sendBuf[0], buf, count, MpiMap<Complex<R> >::type, 
                op, comm ) );
# endif
        }
#else
# ifdef HAVE_MPI_IN_PLACE
        SafeMpi
        ( MPI_Allreduce
          ( MPI_IN_PLACE, buf, count, MpiMap<Complex<R> >::type, op, comm ) );
# else
        std::vector<Complex<R> > sendBuf( count );
        MemCopy( &sendBuf[0], buf, count );
        SafeMpi
        ( MPI_Allreduce
          ( &sendBuf[0], buf, count, MpiMap<Complex<R> >::type, op, comm ) );
# endif
#endif
    }
}

template void AllReduce( byte* buf, int count, Op op, Comm comm );
template void AllReduce( int* buf, int count, Op op, Comm comm );
template void AllReduce( float* buf, int count, Op op, Comm comm );
template void AllReduce( double* buf, int count, Op op, Comm comm );
template void AllReduce( Complex<float>* buf, int count, Op op, Comm comm );
template void AllReduce( Complex<double>* buf, int count, Op op, Comm comm );
template void AllReduce( ValueInt<int>* buf, int count, Op op, Comm comm );
template void AllReduce( ValueInt<float>* buf, int count, Op op, Comm comm );
template void AllReduce( ValueInt<double>* buf, int count, Op op, Comm comm );

template<typename R>
void ReduceScatter( R* sbuf, R* rbuf, int rc, Op op, Comm comm )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::ReduceScatter");
#endif
#ifdef REDUCE_SCATTER_BLOCK_VIA_ALLREDUCE
    const int commSize = CommSize( comm );
    const int commRank = CommRank( comm );
    AllReduce( sbuf, rc*commSize, op, comm );
    MemCopy( rbuf, &sbuf[commRank*rc], rc );
#elif defined(HAVE_MPI_REDUCE_SCATTER_BLOCK)
    SafeMpi
    ( MPI_Reduce_scatter_block( sbuf, rbuf, rc, MpiMap<R>::type, op, comm ) );
#else
    const int commSize = CommSize( comm );
    Reduce( sbuf, rc*commSize, op, 0, comm );
    Scatter( sbuf, rc, rbuf, rc, 0, comm );
#endif
}

template<typename R>
void ReduceScatter
( Complex<R>* sbuf, Complex<R>* rbuf, int rc, Op op, Comm comm )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::ReduceScatter");
#endif
#ifdef REDUCE_SCATTER_BLOCK_VIA_ALLREDUCE
    const int commSize = CommSize( comm );
    const int commRank = CommRank( comm );
    AllReduce( sbuf, rc*commSize, op, comm );
    MemCopy( rbuf, &sbuf[commRank*rc], rc );
#elif defined(HAVE_MPI_REDUCE_SCATTER_BLOCK)
# ifdef AVOID_COMPLEX_MPI
    SafeMpi
    ( MPI_Reduce_scatter_block( sbuf, rbuf, 2*rc, MpiMap<R>::type, op, comm ) );
# else
    SafeMpi
    ( MPI_Reduce_scatter_block
      ( sbuf, rbuf, rc, MpiMap<Complex<R> >::type, op, comm ) );
# endif
#else
    const int commSize = CommSize( comm );
    Reduce( sbuf, rc*commSize, op, 0, comm );
    Scatter( sbuf, rc, rbuf, rc, 0, comm );
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
    CallStackEntry entry("mpi::ReduceScatter");
#endif
#ifdef REDUCE_SCATTER_BLOCK_VIA_ALLREDUCE
    const int commSize = CommSize( comm );
    const int commRank = CommRank( comm );
    AllReduce( buf, rc*commSize, op, comm );
    if( commRank != 0 )
        MemCopy( buf, &buf[commRank*rc], rc );
#elif defined(HAVE_MPI_REDUCE_SCATTER_BLOCK)
# ifdef HAVE_MPI_IN_PLACE
    SafeMpi
    ( MPI_Reduce_scatter_block
      ( MPI_IN_PLACE, buf, rc, MpiMap<R>::type, op, comm ) );
# else
    const int commSize = CommSize( comm );
    std::vector<R> sendBuf( rc*commSize );
    MemCopy( &sendBuf[0], buf, rc*commSize );
    SafeMpi
    ( MPI_Reduce_scatter_block
      ( &sendBuf[0], buf, rc, MpiMap<R>::type, op, comm ) );
# endif
#else
    const int commSize = CommSize( comm );
    Reduce( buf, rc*commSize, op, 0, comm );
    Scatter( buf, rc, rc, 0, comm );
#endif
}

// TODO: Handle case where op is not summation
template<typename R>
void ReduceScatter( Complex<R>* buf, int rc, Op op, Comm comm )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::ReduceScatter");
#endif
#ifdef REDUCE_SCATTER_BLOCK_VIA_ALLREDUCE
    const int commSize = CommSize( comm );
    const int commRank = CommRank( comm );
    AllReduce( buf, rc*commSize, op, comm );
    if( commRank != 0 )
        MemCopy( buf, &buf[commRank*rc], rc );
#elif defined(HAVE_MPI_REDUCE_SCATTER_BLOCK)
# ifdef AVOID_COMPLEX_MPI
#  ifdef HAVE_MPI_IN_PLACE
    SafeMpi
    ( MPI_Reduce_scatter_block
      ( MPI_IN_PLACE, buf, 2*rc, MpiMap<R>::type, op, comm ) );
#  else 
    const int commSize = CommSize( comm );
    std::vector<Complex<R> > sendBuf( rc*commSize );
    MemCopy( &sendBuf[0], buf, rc*commSize );
    SafeMpi
    ( MPI_Reduce_scatter_block
      ( &sendBuf[0], buf, 2*rc, MpiMap<R>::type, op, comm ) );
#  endif
# else
#  ifdef HAVE_MPI_IN_PLACE
    SafeMpi
    ( MPI_Reduce_scatter_block
      ( MPI_IN_PLACE, buf, rc, MpiMap<Complex<R> >::type, op, comm ) );
#  else
    const int commSize = CommSize( comm );
    std::vector<Complex<R> > sendBuf( rc*commSize );
    MemCopy( &sendBuf[0], buf, rc*commSize );
    SafeMpi
    ( MPI_Reduce_scatter_block
      ( &sendBuf[0], buf, rc, MpiMap<Complex<R> >::type, op, comm ) );
#  endif
# endif
#else
    const int commSize = CommSize( comm );
    Reduce( buf, rc*commSize, op, 0, comm );
    Scatter( buf, rc, rc, 0, comm );
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
    CallStackEntry entry("mpi::ReduceScatter");
#endif
    SafeMpi
    ( MPI_Reduce_scatter
      ( const_cast<R*>(sbuf), 
        rbuf, const_cast<int*>(rcs), MpiMap<R>::type, op, comm ) );
}

template<typename R>
void ReduceScatter
( const Complex<R>* sbuf, Complex<R>* rbuf, const int* rcs, Op op, Comm comm )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::ReduceScatter");
#endif
#ifdef AVOID_COMPLEX_MPI
    if( op == SUM )
    {
        int p;
        MPI_Comm_size( comm, &p );
        std::vector<int> rcsDoubled(p);
        for( int i=0; i<p; ++i )
            rcsDoubled[i] = 2*rcs[i];
        SafeMpi
        ( MPI_Reduce_scatter
          ( const_cast<Complex<R>*>(sbuf),
            rbuf, &rcsDoubled[0], MpiMap<R>::type, op, comm ) );
    }
    else
    {
        SafeMpi
        ( MPI_Reduce_scatter
          ( const_cast<Complex<R>*>(sbuf),
            rbuf, const_cast<int*>(rcs), MpiMap<Complex<R> >::type, 
            op, comm ) );
    }
#else
    SafeMpi
    ( MPI_Reduce_scatter
      ( const_cast<Complex<R>*>(sbuf), 
        rbuf, const_cast<int*>(rcs), MpiMap<Complex<R> >::type, op, comm ) );
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
