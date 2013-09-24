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
        elem::RuntimeError( std::string(errorString) );
    }
#endif
}

} // anonymous namespace

namespace elem {
namespace mpi {

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
    SafeMpi( MPI_Waitall( numRequests, requests, statuses.data() ) );
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
bool IProbe( int source, Comm comm, Status& status )
{ return IProbe( source, 0, comm, status ); }

template<typename T>
int GetCount( Status& status )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::GetCount");
#endif
    int count;
    SafeMpi( MPI_Get_count( &status, TypeMap<T>(), &count ) );
    return count;
}
template int GetCount<byte>( Status& status );
template int GetCount<int>( Status& status );
template int GetCount<unsigned>( Status& status );
template int GetCount<long int>( Status& status );
template int GetCount<unsigned long>( Status& status );
#ifdef HAVE_MPI_LONG_LONG
template int GetCount<long long int>( Status& status );
template int GetCount<unsigned long long>( Status& status );
#endif
template int GetCount<float>( Status& status );
template int GetCount<double>( Status& status );
template int GetCount<Complex<float>>( Status& status );
template int GetCount<Complex<double>>( Status& status );

template<typename R>
void TaggedSend( const R* buf, int count, int to, int tag, Comm comm )
{ 
#ifndef RELEASE
    CallStackEntry entry("mpi::Send");
#endif
    SafeMpi
    ( MPI_Send( const_cast<R*>(buf), count, TypeMap<R>(), to, tag, comm ) );
}

template<typename R>
void TaggedSend( const Complex<R>* buf, int count, int to, int tag, Comm comm )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::Send");
#endif
#ifdef AVOID_COMPLEX_MPI
    SafeMpi
    ( MPI_Send
      ( const_cast<Complex<R>*>(buf), 2*count, TypeMap<R>(), to, 
        tag, comm ) );
#else
    SafeMpi
    ( MPI_Send
      ( const_cast<Complex<R>*>(buf), count, 
        TypeMap<Complex<R>>(), to, tag, comm ) );
#endif
}

template void TaggedSend( const byte* buf, int count, int to, int tag, Comm comm );
template void TaggedSend( const int* buf, int count, int to, int tag, Comm comm );
template void TaggedSend( const unsigned* buf, int count, int to, int tag, Comm comm  );
template void TaggedSend( const long int* buf, int count, int to, int tag, Comm comm );
template void TaggedSend( const unsigned long* buf, int count, int to, int tag, Comm comm  );
#ifdef HAVE_MPI_LONG_LONG
template void TaggedSend( const long long int* buf, int count, int to, int tag, Comm comm );
template void TaggedSend( const unsigned long long* buf, int count, int to, int tag, Comm comm  );
#endif
template void TaggedSend( const float* buf, int count, int to, int tag, Comm comm );
template void TaggedSend( const double* buf, int count, int to, int tag, Comm comm );
template void TaggedSend( const Complex<float>* buf, int count, int to, int tag, Comm comm );
template void TaggedSend( const Complex<double>* buf, int count, int to, int tag, Comm comm );

template<typename T>
void Send( const T* buf, int count, int to, Comm comm )
{ TaggedSend( buf, count, to, 0, comm ); }

template void Send( const byte* buf, int count, int to, Comm comm );
template void Send( const int* buf, int count, int to, Comm comm );
template void Send( const unsigned* buf, int count, int to, Comm comm );
template void Send( const long int* buf, int count, int to, Comm comm );
template void Send( const unsigned long* buf, int count, int to, Comm comm );
#ifdef HAVE_MPI_LONG_LONG
template void Send( const long long int* buf, int count, int to, Comm comm );
template void Send( const unsigned long long* buf, int count, int to, Comm comm );
#endif
template void Send( const float* buf, int count, int to, Comm comm );
template void Send( const double* buf, int count, int to, Comm comm );
template void Send( const Complex<float>* buf, int count, int to, Comm comm );
template void Send( const Complex<double>* buf, int count, int to, Comm comm );

template<typename T>
void TaggedSend( T b, int to, int tag, Comm comm )
{ TaggedSend( &b, 1, to, tag, comm ); }

template void TaggedSend( byte b, int to, int tag, Comm comm );
template void TaggedSend( int b, int to, int tag, Comm comm );
template void TaggedSend( unsigned b, int to, int tag, Comm comm );
template void TaggedSend( long int b, int to, int tag, Comm comm );
template void TaggedSend( unsigned long b, int to, int tag, Comm comm );
#ifdef HAVE_MPI_LONG_LONG
template void TaggedSend( long long int b, int to, int tag, Comm comm );
template void TaggedSend( unsigned long long b, int to, int tag, Comm comm );
#endif
template void TaggedSend( float b, int to, int tag, Comm comm );
template void TaggedSend( double b, int to, int tag, Comm comm );
template void TaggedSend( Complex<float> b, int to, int tag, Comm comm );
template void TaggedSend( Complex<double> b, int to, int tag, Comm comm );

template<typename T>
void Send( T b, int to, Comm comm )
{ TaggedSend( b, to, 0, comm ); }

template void Send( byte b, int to, Comm comm );
template void Send( int b, int to, Comm comm );
template void Send( unsigned b, int to, Comm comm );
template void Send( long int b, int to, Comm comm );
template void Send( unsigned long b, int to, Comm comm );
#ifdef HAVE_MPI_LONG_LONG
template void Send( long long int b, int to, Comm comm );
template void Send( unsigned long long b, int to, Comm comm );
#endif
template void Send( float b, int to, Comm comm );
template void Send( double b, int to, Comm comm );
template void Send( Complex<float> b, int to, Comm comm );
template void Send( Complex<double> b, int to, Comm comm );

template<typename R>
void TaggedISend
( const R* buf, int count, int to, int tag, Comm comm, Request& request )
{ 
#ifndef RELEASE
    CallStackEntry entry("mpi::ISend");
#endif
    SafeMpi
    ( MPI_Isend
      ( const_cast<R*>(buf), count, TypeMap<R>(), to, 
        tag, comm, &request ) );
}

template<typename R>
void TaggedISend
( const Complex<R>* buf, int count, int to, int tag, Comm comm, 
  Request& request )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::ISend");
#endif
#ifdef AVOID_COMPLEX_MPI
    SafeMpi
    ( MPI_Isend
      ( const_cast<Complex<R>*>(buf), 2*count, TypeMap<R>(), to, tag, comm,
        &request ) );
#else
    SafeMpi
    ( MPI_Isend
      ( const_cast<Complex<R>*>(buf), count, 
        TypeMap<Complex<R>>(), to, tag, comm, &request ) );
#endif
}

template void TaggedISend( const byte* buf, int count, int to, int tag, Comm comm, Request& request );
template void TaggedISend( const int* buf, int count, int to, int tag, Comm comm, Request& request );
template void TaggedISend( const unsigned* buf, int count, int to, int tag, Comm comm, Request& request );
template void TaggedISend( const long int* buf, int count, int to, int tag, Comm comm, Request& request );
template void TaggedISend( const unsigned long* buf, int count, int to, int tag, Comm comm, Request& request );
#ifdef HAVE_MPI_LONG_LONG
template void TaggedISend( const long long int* buf, int count, int to, int tag, Comm comm, Request& request );
template void TaggedISend( const unsigned long long* buf, int count, int to, int tag, Comm comm, Request& request );
#endif
template void TaggedISend( const float* buf, int count, int to, int tag, Comm comm, Request& request );
template void TaggedISend( const double* buf, int count, int to, int tag, Comm comm, Request& request );
template void TaggedISend( const Complex<float>* buf, int count, int to, int tag, Comm comm, Request& request );
template void TaggedISend( const Complex<double>* buf, int count, int to, int tag, Comm comm, Request& request );

template<typename T>
void ISend
( const T* buf, int count, int to, Comm comm, Request& request )
{ TaggedISend( buf, count, to, 0, comm, request ); } 

template void ISend( const byte* buf, int count, int to, Comm comm, Request& request );
template void ISend( const int* buf, int count, int to, Comm comm, Request& request );
template void ISend( const unsigned* buf, int count, int to, Comm comm, Request& request );
template void ISend( const long int* buf, int count, int to, Comm comm, Request& request );
template void ISend( const unsigned long* buf, int count, int to, Comm comm, Request& request );
#ifdef HAVE_MPI_LONG_LONG
template void ISend( const long long int* buf, int count, int to, Comm comm, Request& request );
template void ISend( const unsigned long long* buf, int count, int to, Comm comm, Request& request );
#endif
template void ISend( const float* buf, int count, int to, Comm comm, Request& request );
template void ISend( const double* buf, int count, int to, Comm comm, Request& request );
template void ISend( const Complex<float>* buf, int count, int to, Comm comm, Request& request );
template void ISend( const Complex<double>* buf, int count, int to, Comm comm, Request& request );

template<typename T>
void TaggedISend( T b, int to, int tag, Comm comm, Request& request )
{ TaggedISend( &b, 1, to, tag, comm, request ); }

template void TaggedISend( byte buf, int to, int tag, Comm comm, Request& request );
template void TaggedISend( int buf, int to, int tag, Comm comm, Request& request );
template void TaggedISend( unsigned buf, int to, int tag, Comm comm, Request& request );
template void TaggedISend( long int buf, int to, int tag, Comm comm, Request& request );
template void TaggedISend( unsigned long buf, int to, int tag, Comm comm, Request& request );
#ifdef HAVE_MPI_LONG_LONG
template void TaggedISend( long long int buf, int to, int tag, Comm comm, Request& request );
template void TaggedISend( unsigned long long buf, int to, int tag, Comm comm, Request& request );
#endif
template void TaggedISend( float buf, int to, int tag, Comm comm, Request& request );
template void TaggedISend( double buf, int to, int tag, Comm comm, Request& request );
template void TaggedISend( Complex<float> buf, int to, int tag, Comm comm, Request& request );
template void TaggedISend( Complex<double> buf, int to, int tag, Comm comm, Request& request );

template<typename T>
void ISend( T b, int to, Comm comm, Request& request )
{ TaggedISend( b, to, 0, comm, request ); }

template void ISend( byte buf, int to, Comm comm, Request& request );
template void ISend( int buf, int to, Comm comm, Request& request );
template void ISend( unsigned buf, int to, Comm comm, Request& request );
template void ISend( long int buf, int to, Comm comm, Request& request );
template void ISend( unsigned long buf, int to, Comm comm, Request& request );
#ifdef HAVE_MPI_LONG_LONG
template void ISend( long long int buf, int to, Comm comm, Request& request );
template void ISend( unsigned long long buf, int to, Comm comm, Request& request );
#endif
template void ISend( float buf, int to, Comm comm, Request& request );
template void ISend( double buf, int to, Comm comm, Request& request );
template void ISend( Complex<float> buf, int to, Comm comm, Request& request );
template void ISend( Complex<double> buf, int to, Comm comm, Request& request );

template<typename R>
void TaggedISSend
( const R* buf, int count, int to, int tag, Comm comm, Request& request )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::ISSend");
#endif
    SafeMpi
    ( MPI_Issend
      ( const_cast<R*>(buf), count, TypeMap<R>(), to, 
        tag, comm, &request ) );
}

template<typename R>
void TaggedISSend
( const Complex<R>* buf, int count, int to, int tag, Comm comm, 
  Request& request )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::ISSend");
#endif
#ifdef AVOID_COMPLEX_MPI
    SafeMpi
    ( MPI_Issend
      ( const_cast<Complex<R>*>(buf), 2*count, TypeMap<R>(), to, tag, comm,
        &request ) );
#else
    SafeMpi
    ( MPI_Issend
      ( const_cast<Complex<R>*>(buf), count, 
        TypeMap<Complex<R>>(), to, tag, comm, &request ) );
#endif
}

template void TaggedISSend( const byte* buf, int count, int to, int tag, Comm comm, Request& request );
template void TaggedISSend( const int* buf, int count, int to, int tag, Comm comm, Request& request );
template void TaggedISSend( const unsigned* buf, int count, int to, int tag, Comm comm, Request& request );
template void TaggedISSend( const long int* buf, int count, int to, int tag, Comm comm, Request& request );
template void TaggedISSend( const unsigned long* buf, int count, int to, int tag, Comm comm, Request& request );
#ifdef HAVE_MPI_LONG_LONG
template void TaggedISSend( const long long int* buf, int count, int to, int tag, Comm comm, Request& request );
template void TaggedISSend( const unsigned long long* buf, int count, int to, int tag, Comm comm, Request& request );
#endif
template void TaggedISSend( const float* buf, int count, int to, int tag, Comm comm, Request& request );
template void TaggedISSend( const double* buf, int count, int to, int tag, Comm comm, Request& request );
template void TaggedISSend( const Complex<float>* buf, int count, int to, int tag, Comm comm, Request& request );
template void TaggedISSend( const Complex<double>* buf, int count, int to, int tag, Comm comm, Request& request );

template<typename T>
void ISSend( const T* buf, int count, int to, Comm comm, Request& request )
{ TaggedISSend( buf, count, to, 0, comm, request ); }

template void ISSend( const byte* buf, int count, int to, Comm comm, Request& request );
template void ISSend( const int* buf, int count, int to, Comm comm, Request& request );
template void ISSend( const unsigned* buf, int count, int to, Comm comm, Request& request );
template void ISSend( const long int* buf, int count, int to, Comm comm, Request& request );
template void ISSend( const unsigned long* buf, int count, int to, Comm comm, Request& request );
#ifdef HAVE_MPI_LONG_LONG
template void ISSend( const long long int* buf, int count, int to, Comm comm, Request& request );
template void ISSend( const unsigned long long* buf, int count, int to, Comm comm, Request& request );
#endif
template void ISSend( const float* buf, int count, int to, Comm comm, Request& request );
template void ISSend( const double* buf, int count, int to, Comm comm, Request& request );
template void ISSend( const Complex<float>* buf, int count, int to, Comm comm, Request& request );
template void ISSend( const Complex<double>* buf, int count, int to, Comm comm, Request& request );

template<typename T>
void TaggedISSend( T b, int to, int tag, Comm comm, Request& request )
{ TaggedISSend( &b, 1, to, tag, comm, request ); }

template void TaggedISSend( byte b, int to, int tag, Comm comm, Request& request );
template void TaggedISSend( int b, int to, int tag, Comm comm, Request& request );
template void TaggedISSend( unsigned b, int to, int tag, Comm comm, Request& request );
template void TaggedISSend( long int b, int to, int tag, Comm comm, Request& request );
template void TaggedISSend( unsigned long b, int to, int tag, Comm comm, Request& request );
#ifdef HAVE_MPI_LONG_LONG
template void TaggedISSend( long long int b, int to, int tag, Comm comm, Request& request );
template void TaggedISSend( unsigned long long b, int to, int tag, Comm comm, Request& request );
#endif
template void TaggedISSend( float b, int to, int tag, Comm comm, Request& request );
template void TaggedISSend( double b, int to, int tag, Comm comm, Request& request );
template void TaggedISSend( Complex<float> b, int to, int tag, Comm comm, Request& request );
template void TaggedISSend( Complex<double> b, int to, int tag, Comm comm, Request& request );

template<typename R>
void TaggedRecv( R* buf, int count, int from, int tag, Comm comm )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::Recv");
#endif
    Status status;
    SafeMpi
    ( MPI_Recv( buf, count, TypeMap<R>(), from, tag, comm, &status ) );
}

template<typename R>
void TaggedRecv( Complex<R>* buf, int count, int from, int tag, Comm comm )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::Recv");
#endif
    Status status;
#ifdef AVOID_COMPLEX_MPI
    SafeMpi
    ( MPI_Recv( buf, 2*count, TypeMap<R>(), from, tag, comm, &status ) );
#else
    SafeMpi
    ( MPI_Recv
      ( buf, count, TypeMap<Complex<R>>(), from, tag, comm, &status ) );
#endif
}

template void TaggedRecv( byte* buf, int count, int from, int tag, Comm comm );
template void TaggedRecv( int* buf, int count, int from, int tag, Comm comm );
template void TaggedRecv( unsigned* buf, int count, int from, int tag, Comm comm );
template void TaggedRecv( long int* buf, int count, int from, int tag, Comm comm );
template void TaggedRecv( unsigned long* buf, int count, int from, int tag, Comm comm );
#ifdef HAVE_MPI_LONG_LONG
template void TaggedRecv( long long int* buf, int count, int from, int tag, Comm comm );
template void TaggedRecv( unsigned long long* buf, int count, int from, int tag, Comm comm );
#endif
template void TaggedRecv( float* buf, int count, int from, int tag, Comm comm );
template void TaggedRecv( double* buf, int count, int from, int tag, Comm comm );
template void TaggedRecv( Complex<float>* buf, int count, int from, int tag, Comm comm );
template void TaggedRecv( Complex<double>* buf, int count, int from, int tag, Comm comm );

template<typename T>
void Recv( T* buf, int count, int from, Comm comm )
{ TaggedRecv( buf, count, from, mpi::ANY_TAG, comm ); }

template void Recv( byte* buf, int count, int from, Comm comm );
template void Recv( int* buf, int count, int from, Comm comm );
template void Recv( unsigned* buf, int count, int from, Comm comm );
template void Recv( long int* buf, int count, int from, Comm comm );
template void Recv( unsigned long* buf, int count, int from, Comm comm );
#ifdef HAVE_MPI_LONG_LONG
template void Recv( long long int* buf, int count, int from, Comm comm );
template void Recv( unsigned long long* buf, int count, int from, Comm comm );
#endif
template void Recv( float* buf, int count, int from, Comm comm );
template void Recv( double* buf, int count, int from, Comm comm );
template void Recv( Complex<float>* buf, int count, int from, Comm comm );
template void Recv( Complex<double>* buf, int count, int from, Comm comm );

template<typename T>
T TaggedRecv( int from, int tag, Comm comm )
{ T b; TaggedRecv( &b, 1, from, tag, comm ); return b; }

template byte TaggedRecv( int from, int tag, Comm comm );
template int TaggedRecv( int from, int tag, Comm comm );
template unsigned TaggedRecv( int from, int tag, Comm comm );
template long int TaggedRecv( int from, int tag, Comm comm );
template unsigned long TaggedRecv( int from, int tag, Comm comm );
#ifdef HAVE_MPI_LONG_LONG
template long long int TaggedRecv( int from, int tag, Comm comm );
template unsigned long long TaggedRecv( int from, int tag, Comm comm );
#endif
template float TaggedRecv( int from, int tag, Comm comm );
template double TaggedRecv( int from, int tag, Comm comm );
template Complex<float> TaggedRecv( int from, int tag, Comm comm );
template Complex<double> TaggedRecv( int from, int tag, Comm comm );

template<typename T>
T Recv( int from, Comm comm )
{ return TaggedRecv<T>( from, mpi::ANY_TAG, comm ); }

template byte Recv( int from, Comm comm );
template int Recv( int from, Comm comm );
template unsigned Recv( int from, Comm comm );
template long int Recv( int from, Comm comm );
template unsigned long Recv( int from, Comm comm );
#ifdef HAVE_MPI_LONG_LONG
template long long int Recv( int from, Comm comm );
template unsigned long long Recv( int from, Comm comm );
#endif
template float Recv( int from, Comm comm );
template double Recv( int from, Comm comm );
template Complex<float> Recv( int from, Comm comm );
template Complex<double> Recv( int from, Comm comm );

template<typename R>
void TaggedIRecv
( R* buf, int count, int from, int tag, Comm comm, Request& request )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::IRecv");
#endif
    SafeMpi
    ( MPI_Irecv( buf, count, TypeMap<R>(), from, tag, comm, &request ) );
}

template<typename R>
void TaggedIRecv
( Complex<R>* buf, int count, int from, int tag, Comm comm, Request& request )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::IRecv");
#endif
#ifdef AVOID_COMPLEX_MPI
    SafeMpi
    ( MPI_Irecv( buf, 2*count, TypeMap<R>(), from, tag, comm, &request ) );
#else
    SafeMpi
    ( MPI_Irecv
      ( buf, count, TypeMap<Complex<R>>(), from, tag, comm, &request ) );
#endif
}

template void TaggedIRecv( byte* buf, int count, int from, int tag, Comm comm, Request& request );
template void TaggedIRecv( int* buf, int count, int from, int tag, Comm comm, Request& request );
template void TaggedIRecv( unsigned* buf, int count, int from, int tag, Comm comm, Request& request );
template void TaggedIRecv( long int* buf, int count, int from, int tag, Comm comm, Request& request );
template void TaggedIRecv( unsigned long* buf, int count, int from, int tag, Comm comm, Request& request );
#ifdef HAVE_MPI_LONG_LONG
template void TaggedIRecv( long long int* buf, int count, int from, int tag, Comm comm, Request& request );
template void TaggedIRecv( unsigned long long* buf, int count, int from, int tag, Comm comm, Request& request );
#endif
template void TaggedIRecv( float* buf, int count, int from, int tag, Comm comm, Request& request );
template void TaggedIRecv( double* buf, int count, int from, int tag, Comm comm, Request& request );
template void TaggedIRecv( Complex<float>* buf, int count, int from, int tag, Comm comm, Request& request );
template void TaggedIRecv( Complex<double>* buf, int count, int from, int tag, Comm comm, Request& request );

template<typename T>
void IRecv( T* buf, int count, int from, Comm comm, Request& request )
{ TaggedIRecv( buf, count, from, mpi::ANY_TAG, comm, request ); }

template void IRecv( byte* buf, int count, int from, Comm comm, Request& request );
template void IRecv( int* buf, int count, int from, Comm comm, Request& request );
template void IRecv( unsigned* buf, int count, int from, Comm comm, Request& request );
template void IRecv( long int* buf, int count, int from, Comm comm, Request& request );
template void IRecv( unsigned long* buf, int count, int from, Comm comm, Request& request );
#ifdef HAVE_MPI_LONG_LONG
template void IRecv( long long int* buf, int count, int from, Comm comm, Request& request );
template void IRecv( unsigned long long* buf, int count, int from, Comm comm, Request& request );
#endif
template void IRecv( float* buf, int count, int from, Comm comm, Request& request );
template void IRecv( double* buf, int count, int from, Comm comm, Request& request );
template void IRecv( Complex<float>* buf, int count, int from, Comm comm, Request& request );
template void IRecv( Complex<double>* buf, int count, int from, Comm comm, Request& request );

template<typename T>
T TaggedIRecv( int from, int tag, Comm comm, Request& request )
{ T b; TaggedIRecv( &b, 1, from, tag, comm, request ); return b; }

template byte TaggedIRecv( int from, int tag, Comm comm, Request& request );
template int TaggedIRecv( int from, int tag, Comm comm, Request& request );
template unsigned TaggedIRecv( int from, int tag, Comm comm, Request& request );
template long int TaggedIRecv( int from, int tag, Comm comm, Request& request );
template unsigned long TaggedIRecv( int from, int tag, Comm comm, Request& request );
#ifdef HAVE_MPI_LONG_LONG
template long long int TaggedIRecv( int from, int tag, Comm comm, Request& request );
template unsigned long long TaggedIRecv( int from, int tag, Comm comm, Request& request );
#endif
template float TaggedIRecv( int from, int tag, Comm comm, Request& request );
template double TaggedIRecv( int from, int tag, Comm comm, Request& request );
template Complex<float> TaggedIRecv( int from, int tag, Comm comm, Request& request );
template Complex<double> TaggedIRecv( int from, int tag, Comm comm, Request& request );

template<typename T>
T IRecv( int from, Comm comm, Request& request )
{ return TaggedIRecv<T>( from, mpi::ANY_TAG, comm, request ); }

template byte IRecv( int from, Comm comm, Request& request );
template int IRecv( int from, Comm comm, Request& request );
template unsigned IRecv( int from, Comm comm, Request& request );
template long int IRecv( int from, Comm comm, Request& request );
template unsigned long IRecv( int from, Comm comm, Request& request );
#ifdef HAVE_MPI_LONG_LONG
template long long int IRecv( int from, Comm comm, Request& request );
template unsigned long long IRecv( int from, Comm comm, Request& request );
#endif
template float IRecv( int from, Comm comm, Request& request );
template double IRecv( int from, Comm comm, Request& request );
template Complex<float> IRecv( int from, Comm comm, Request& request );
template Complex<double> IRecv( int from, Comm comm, Request& request );

template<typename R>
void TaggedSendRecv
( const R* sbuf, int sc, int to,   int stag,
        R* rbuf, int rc, int from, int rtag, Comm comm )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::SendRecv");
#endif
    Status status;
    SafeMpi
    ( MPI_Sendrecv
      ( const_cast<R*>(sbuf), sc, TypeMap<R>(), to,   stag,
        rbuf,                 rc, TypeMap<R>(), from, rtag, 
        comm, &status ) );
}

template<typename R>
void TaggedSendRecv
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
      ( const_cast<Complex<R>*>(sbuf), 2*sc, TypeMap<R>(), to,   stag,
        rbuf,                          2*rc, TypeMap<R>(), from, rtag, 
        comm, &status ) );
#else
    SafeMpi
    ( MPI_Sendrecv
      ( const_cast<Complex<R>*>(sbuf), 
        sc, TypeMap<Complex<R>>(), to,   stag,
        rbuf,                          
        rc, TypeMap<Complex<R>>(), from, rtag, comm, &status ) );
#endif
}

template void TaggedSendRecv
( const byte* sbuf, int sc, int to, int stag, 
        byte* rbuf, int rc, int from, int rtag, Comm comm );
template void TaggedSendRecv
( const int* sbuf, int sc, int to, int stag, 
        int* rbuf, int rc, int from, int rtag, Comm comm );
template void TaggedSendRecv
( const unsigned* sbuf, int sc, int to, int stag, 
        unsigned* rbuf, int rc, int from, int rtag, Comm comm );
template void TaggedSendRecv
( const long int* sbuf, int sc, int to, int stag, 
        long int* rbuf, int rc, int from, int rtag, Comm comm );
template void TaggedSendRecv
( const unsigned long* sbuf, int sc, int to, int stag, 
        unsigned long* rbuf, int rc, int from, int rtag, Comm comm );
#ifdef HAVE_MPI_LONG_LONG
template void TaggedSendRecv
( const long long int* sbuf, int sc, int to, int stag, 
        long long int* rbuf, int rc, int from, int rtag, Comm comm );
template void TaggedSendRecv
( const unsigned long long* sbuf, int sc, int to, int stag, 
        unsigned long long* rbuf, int rc, int from, int rtag, Comm comm );
#endif
template void TaggedSendRecv
( const float* sbuf, int sc, int to, int stag, 
        float* rbuf, int rc, int from, int rtag, Comm comm );
template void TaggedSendRecv
( const double* sbuf, int sc, int to, int stag, 
        double* rbuf, int rc, int from, int rtag, Comm comm );
template void TaggedSendRecv
( const Complex<float>* sbuf, int sc, int to, int stag, 
        Complex<float>* rbuf, int rc, int from, int rtag, Comm comm );
template void TaggedSendRecv
( const Complex<double>* sbuf, int sc, int to, int stag, 
        Complex<double>* rbuf, int rc, int from, int rtag, Comm comm );

template<typename T>
void SendRecv
( const T* sbuf, int sc, int to, 
        T* rbuf, int rc, int from, Comm comm )
{ TaggedSendRecv( sbuf, sc, to, 0, rbuf, rc, from, mpi::ANY_TAG, comm ); }

template void SendRecv
( const byte* sbuf, int sc, int to, 
        byte* rbuf, int rc, int from, Comm comm );
template void SendRecv
( const int* sbuf, int sc, int to, 
        int* rbuf, int rc, int from, Comm comm );
template void SendRecv
( const unsigned* sbuf, int sc, int to, 
        unsigned* rbuf, int rc, int from, Comm comm );
template void SendRecv
( const long int* sbuf, int sc, int to, 
        long int* rbuf, int rc, int from, Comm comm );
template void SendRecv
( const unsigned long* sbuf, int sc, int to, 
        unsigned long* rbuf, int rc, int from, Comm comm );
#ifdef HAVE_MPI_LONG_LONG
template void SendRecv
( const long long int* sbuf, int sc, int to, 
        long long int* rbuf, int rc, int from, Comm comm );
template void SendRecv
( const unsigned long long* sbuf, int sc, int to, 
        unsigned long long* rbuf, int rc, int from, Comm comm );
#endif
template void SendRecv
( const float* sbuf, int sc, int to,
        float* rbuf, int rc, int from, Comm comm );
template void SendRecv
( const double* sbuf, int sc, int to,
        double* rbuf, int rc, int from, Comm comm );
template void SendRecv
( const Complex<float>* sbuf, int sc, int to, 
        Complex<float>* rbuf, int rc, int from, Comm comm );
template void SendRecv
( const Complex<double>* sbuf, int sc, int to, 
        Complex<double>* rbuf, int rc, int from, Comm comm );

template<typename T>
T TaggedSendRecv( T sb, int to, int stag, int from, int rtag, Comm comm )
{ 
    T rb; 
    TaggedSendRecv( &sb, 1, to, stag, &rb, 1, from, rtag, comm ); 
    return rb; 
}

template byte TaggedSendRecv
( byte sb, int to, int stag, int from, int rtag, Comm comm );
template int TaggedSendRecv
( int sb, int to, int stag, int from, int rtag, Comm comm );
template unsigned TaggedSendRecv
( unsigned sb, int to, int stag, int from, int rtag, Comm comm );
template long int TaggedSendRecv
( long int sb, int to, int stag, int from, int rtag, Comm comm );
template unsigned long TaggedSendRecv
( unsigned long sb, int to, int stag, int from, int rtag, Comm comm );
#ifdef HAVE_MPI_LONG_LONG
template long long int TaggedSendRecv
( long long int sb, int to, int stag, int from, int rtag, Comm comm );
template unsigned long long TaggedSendRecv
( unsigned long long sb, int to, int stag, int from, int rtag, Comm comm );
#endif
template float TaggedSendRecv
( float sb, int to, int stag, int from, int rtag, Comm comm );
template double TaggedSendRecv
( double sb, int to, int stag, int from, int rtag, Comm comm );
template Complex<float> TaggedSendRecv
( Complex<float> sb, int to, int stag, int from, int rtag, Comm comm );
template Complex<double> TaggedSendRecv
( Complex<double> sb, int to, int stag, int from, int rtag, Comm comm );

template<typename T>
T SendRecv( T sb, int to, int from, Comm comm )
{ return TaggedSendRecv( sb, to, 0, from, mpi::ANY_TAG, comm ); }

template byte SendRecv( byte sb, int to, int from, Comm comm );
template int SendRecv( int sb, int to, int from, Comm comm );
template unsigned SendRecv( unsigned sb, int to, int from, Comm comm );
template long int SendRecv( long int sb, int to, int from, Comm comm );
template unsigned long SendRecv( unsigned long sb, int to, int from, Comm comm );
#ifdef HAVE_MPI_LONG_LONG
template long long int SendRecv( long long int sb, int to, int from, Comm comm );
template unsigned long long SendRecv( unsigned long long sb, int to, int from, Comm comm );
#endif
template float SendRecv( float sb, int to, int from, Comm comm );
template double SendRecv( double sb, int to, int from, Comm comm );
template Complex<float> SendRecv
( Complex<float> sb, int to, int from, Comm comm );
template Complex<double> SendRecv
( Complex<double> sb, int to, int from, Comm comm );

template<typename R>
void TaggedSendRecv
( R* buf, int count, int to, int stag, int from, int rtag, Comm comm )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::SendRecv");
#endif
    Status status;
    SafeMpi
    ( MPI_Sendrecv_replace
      ( buf, count, TypeMap<R>(), to, stag, from, rtag, comm, &status ) );
}

template<typename R>
void TaggedSendRecv
( Complex<R>* buf, int count, int to, int stag, int from, int rtag, Comm comm )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::SendRecv");
#endif
    Status status;
#ifdef AVOID_COMPLEX_MPI
    SafeMpi
    ( MPI_Sendrecv_replace
      ( buf, 2*count, TypeMap<R>(), to, stag, from, rtag, comm, &status ) );
#else
    SafeMpi
    ( MPI_Sendrecv_replace
      ( buf, count, TypeMap<Complex<R>>(), 
        to, stag, from, rtag, comm, &status ) );
#endif
}

template void TaggedSendRecv
( byte* buf, int count, int to, int stag, int from, int rtag, Comm comm );
template void TaggedSendRecv
( int* buf, int count, int to, int stag, int from, int rtag, Comm comm );
template void TaggedSendRecv
( unsigned* buf, int count, int to, int stag, int from, int rtag, Comm comm );
template void TaggedSendRecv
( long int* buf, int count, int to, int stag, int from, int rtag, Comm comm );
template void TaggedSendRecv
( unsigned long* buf, int count, int to, int stag, int from, int rtag, Comm comm );
#ifdef HAVE_MPI_LONG_LONG
template void TaggedSendRecv
( long long int* buf, int count, int to, int stag, int from, int rtag, Comm comm );
template void TaggedSendRecv
( unsigned long long* buf, int count, int to, int stag, int from, int rtag, Comm comm );
#endif
template void TaggedSendRecv
( float* buf, int count, int to, int stag, int from, int rtag, Comm comm );
template void TaggedSendRecv
( double* buf, int count, int to, int stag, int from, int rtag, Comm comm );
template void TaggedSendRecv
( Complex<float>* buf, int count, int to, int stag, 
  int from, int rtag, Comm comm );
template void TaggedSendRecv
( Complex<double>* buf, int count, int to, int stag, 
  int from, int rtag, Comm comm );

template<typename T>
void SendRecv( T* buf, int count, int to, int from, Comm comm )
{ TaggedSendRecv( buf, count, to, 0, from, mpi::ANY_TAG, comm ); }

template void SendRecv
( byte* buf, int count, int to, int from, Comm comm );
template void SendRecv
( int* buf, int count, int to, int from, Comm comm );
template void SendRecv
( unsigned* buf, int count, int to, int from, Comm comm );
template void SendRecv
( long int* buf, int count, int to, int from, Comm comm );
template void SendRecv
( unsigned long* buf, int count, int to, int from, Comm comm );
#ifdef HAVE_MPI_LONG_LONG
template void SendRecv
( long long int* buf, int count, int to, int from, Comm comm );
template void SendRecv
( unsigned long long* buf, int count, int to, int from, Comm comm );
#endif
template void SendRecv
( float* buf, int count, int to, int from, Comm comm );
template void SendRecv
( double* buf, int count, int to, int from, Comm comm );
template void SendRecv
( Complex<float>* buf, int count, int to, int from, Comm comm );
template void SendRecv
( Complex<double>* buf, int count, int to, int from, Comm comm );

template<typename R>
void Broadcast( R* buf, int count, int root, Comm comm )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::Broadcast");
#endif
    SafeMpi( MPI_Bcast( buf, count, TypeMap<R>(), root, comm ) );
}

template<typename R>
void Broadcast( Complex<R>* buf, int count, int root, Comm comm )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::Broadcast");
#endif
#ifdef AVOID_COMPLEX_MPI
    SafeMpi( MPI_Bcast( buf, 2*count, TypeMap<R>(), root, comm ) );
#else
    SafeMpi( MPI_Bcast( buf, count, TypeMap<Complex<R>>(), root, comm ) );
#endif
}

template void Broadcast( byte* buf, int count, int root, Comm comm );
template void Broadcast( int* buf, int count, int root, Comm comm );
template void Broadcast( unsigned* buf, int count, int root, Comm comm );
template void Broadcast( long int* buf, int count, int root, Comm comm );
template void Broadcast( unsigned long* buf, int count, int root, Comm comm );
#ifdef HAVE_MPI_LONG_LONG
template void Broadcast( long long int* buf, int count, int root, Comm comm );
template void Broadcast( unsigned long long* buf, int count, int root, Comm comm );
#endif
template void Broadcast( float* buf, int count, int root, Comm comm );
template void Broadcast( double* buf, int count, int root, Comm comm );
template void Broadcast( Complex<float>* buf, int count, int root, Comm comm );
template void Broadcast( Complex<double>* buf, int count, int root, Comm comm );

template<typename T>
void Broadcast( T& b, int root, Comm comm )
{ Broadcast( &b, 1, root, comm ); }

template void Broadcast( byte& b, int root, Comm comm );
template void Broadcast( int& b, int root, Comm comm );
template void Broadcast( unsigned& b, int root, Comm comm );
template void Broadcast( long int& b, int root, Comm comm );
template void Broadcast( unsigned long& b, int root, Comm comm );
#ifdef HAVE_MPI_LONG_LONG
template void Broadcast( long long int& b, int root, Comm comm );
template void Broadcast( unsigned long long& b, int root, Comm comm );
#endif
template void Broadcast( float& b, int root, Comm comm );
template void Broadcast( double& b, int root, Comm comm );
template void Broadcast( Complex<float>& b, int root, Comm comm );
template void Broadcast( Complex<double>& b, int root, Comm comm );

#ifdef HAVE_NONBLOCKING_COLLECTIVES
template<typename R>
void IBroadcast( R* buf, int count, int root, Comm comm, Request& request )
{
#ifndef RELEASE
    CallStackEntry entry("mpi::IBroadcast");
#endif
    SafeMpi( MPI_Ibcast( buf, count, TypeMap<R>(), root, comm, &request ) );
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
    ( MPI_Ibcast( buf, 2*count, TypeMap<R>(), root, comm, &request ) );
#else
    SafeMpi
    ( MPI_Ibcast
      ( buf, count, TypeMap<Complex<R>>(), root, comm, &request ) );
#endif
}

template void IBroadcast( byte* buf, int count, int root, Comm comm, Request& request );
template void IBroadcast( int* buf, int count, int root, Comm comm, Request& request );
template void IBroadcast( unsigned* buf, int count, int root, Comm comm, Request& request );
template void IBroadcast( long int* buf, int count, int root, Comm comm, Request& request );
template void IBroadcast( unsigned long* buf, int count, int root, Comm comm, Request& request );
#ifdef HAVE_MPI_LONG_LONG
template void IBroadcast( long long int* buf, int count, int root, Comm comm, Request& request );
template void IBroadcast( unsigned long long* buf, int count, int root, Comm comm, Request& request );
#endif
template void IBroadcast( float* buf, int count, int root, Comm comm, Request& request );
template void IBroadcast( double* buf, int count, int root, Comm comm, Request& request );
template void IBroadcast( Complex<float>* buf, int count, int root, Comm comm, Request& request );
template void IBroadcast( Complex<double>* buf, int count, int root, Comm comm, Request& request );

template<typename T>
void IBroadcast( T& b, int root, Comm comm, Request& request )
{ IBroadcast( &b, 1, root, comm, request ); }

template void IBroadcast( byte& b, int root, Comm comm, Request& request );
template void IBroadcast( int& b, int root, Comm comm, Request& request );
template void IBroadcast( unsigned& b, int root, Comm comm, Request& request );
template void IBroadcast( long int& b, int root, Comm comm, Request& request );
template void IBroadcast( unsigned long& b, int root, Comm comm, Request& request );
#ifdef HAVE_MPI_LONG_LONG
template void IBroadcast( long long int& b, int root, Comm comm, Request& request );
template void IBroadcast( unsigned long long& b, int root, Comm comm, Request& request );
#endif
template void IBroadcast( float& b, int root, Comm comm, Request& request );
template void IBroadcast( double& b, int root, Comm comm, Request& request );
template void IBroadcast( Complex<float>& b, int root, Comm comm, Request& request );
template void IBroadcast( Complex<double>& b, int root, Comm comm, Request& request );
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
      ( const_cast<R*>(sbuf), sc, TypeMap<R>(),
        rbuf,                 rc, TypeMap<R>(), root, comm ) );
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
      ( const_cast<Complex<R>*>(sbuf), 2*sc, TypeMap<R>(),
        rbuf,                          2*rc, TypeMap<R>(), root, comm ) );
#else
    SafeMpi
    ( MPI_Gather
      ( const_cast<Complex<R>*>(sbuf), sc, TypeMap<Complex<R>>(),
        rbuf,                          rc, TypeMap<Complex<R>>(), 
        root, comm ) );
#endif
}

template void Gather( const byte* sbuf, int sc, byte* rbuf, int rc, int root, Comm comm );
template void Gather( const int* sbuf, int sc, int* rbuf, int rc, int root, Comm comm );
template void Gather( const unsigned* sbuf, int sc, unsigned* rbuf, int rc, int root, Comm comm );
template void Gather( const long int* sbuf, int sc, long int* rbuf, int rc, int root, Comm comm );
template void Gather( const unsigned long* sbuf, int sc, unsigned long* rbuf, int rc, int root, Comm comm );
#ifdef HAVE_MPI_LONG_LONG
template void Gather( const long long int* sbuf, int sc, long long int* rbuf, int rc, int root, Comm comm );
template void Gather( const unsigned long long* sbuf, int sc, unsigned long long* rbuf, int rc, int root, Comm comm );
#endif
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
      ( const_cast<R*>(sbuf), sc, TypeMap<R>(),
        rbuf,                 rc, TypeMap<R>(), root, comm, &request ) );
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
      ( const_cast<Complex<R>*>(sbuf), 2*sc, TypeMap<R>(),
        rbuf,                          2*rc, TypeMap<R>(), 
        root, comm, &request ) );
#else
    SafeMpi
    ( MPI_Igather
      ( const_cast<Complex<R>*>(sbuf), sc, TypeMap<Complex<R>>(),
        rbuf,                          rc, TypeMap<Complex<R>>(), 
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
( const unsigned* sbuf, int sc, 
        unsigned* rbuf, int rc, int root, Comm comm, Request& request );
template void IGather
( const long int* sbuf, int sc, 
        long int* rbuf, int rc, int root, Comm comm, Request& request );
template void IGather
( const unsigned long* sbuf, int sc, 
        unsigned long* rbuf, int rc, int root, Comm comm, Request& request );
#ifdef HAVE_MPI_LONG_LONG
template void IGather
( const long long int* sbuf, int sc,
        long long int* rbuf, int rc, int root, Comm comm, Request& request );
template void IGather
( const unsigned long long* sbuf, int sc,
        unsigned long long* rbuf, int rc, int root, Comm comm, Request& request );
#endif
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
        TypeMap<R>(),
        rbuf,                    
        const_cast<int*>(rcs), 
        const_cast<int*>(rds), 
        TypeMap<R>(),
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
      ( const_cast<Complex<R>*>(sbuf), 2*sc, TypeMap<R>(),
        rbuf, rcsDouble.data(), rdsDouble.data(), TypeMap<R>(),
        root, comm ) );
#else
    SafeMpi
    ( MPI_Gatherv
      ( const_cast<Complex<R>*>(sbuf), 
        sc,       
        TypeMap<Complex<R>>(),
        rbuf,  
        const_cast<int*>(rcs), 
        const_cast<int*>(rds), 
        TypeMap<Complex<R>>(),
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
( const unsigned* sbuf, int sc, 
        unsigned* rbuf, const int* rcs, const int* rds, int root, Comm comm );
template void Gather
( const long int* sbuf, int sc, 
        long int* rbuf, const int* rcs, const int* rds, int root, Comm comm );
template void Gather
( const unsigned long* sbuf, int sc, 
        unsigned long* rbuf, const int* rcs, const int* rds, int root, Comm comm );
#ifdef HAVE_MPI_LONG_LONG
template void Gather
( const long long int* sbuf, int sc, 
        long long int* rbuf, const int* rcs, const int* rds, int root, Comm comm );
template void Gather
( const unsigned long long* sbuf, int sc, 
        unsigned long long* rbuf, const int* rcs, const int* rds, int root, Comm comm );
#endif
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
      ( const_cast<R*>(sbuf), sc, TypeMap<R>(), 
        rbuf,                 rc, TypeMap<R>(), comm ) );
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
      ( const_cast<Complex<R>*>(sbuf), 2*sc, TypeMap<R>(),
        rbuf,                          2*rc, TypeMap<R>(), comm ) );
 #else
    SafeMpi
    ( MPI_Allgather
      ( const_cast<Complex<R>*>(sbuf), sc, TypeMap<Complex<R>>(),
        rbuf,                          rc, TypeMap<Complex<R>>(), comm ) );
 #endif
#endif
}

template void AllGather( const byte* sbuf, int sc, byte* rbuf, int rc, Comm comm );
template void AllGather( const int* sbuf, int sc, int* rbuf, int rc, Comm comm );
template void AllGather( const unsigned* sbuf, int sc, unsigned* rbuf, int rc, Comm comm );
template void AllGather( const long int* sbuf, int sc, long int* rbuf, int rc, Comm comm );
template void AllGather( const unsigned long* sbuf, int sc, unsigned long* rbuf, int rc, Comm comm );
#ifdef HAVE_MPI_LONG_LONG
template void AllGather( const long long int* sbuf, int sc, long long int* rbuf, int rc, Comm comm );
template void AllGather( const unsigned long long* sbuf, int sc, unsigned long long* rbuf, int rc, Comm comm );
#endif
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
      ( const_cast<R*>(sbuf), sizeof(R)*sc,   MPI_UNSIGNED_CHAR, 
        rbuf, byteRcs.data(), byteRds.data(), MPI_UNSIGNED_CHAR, comm ) );
#else
    SafeMpi
    ( MPI_Allgatherv
      ( const_cast<R*>(sbuf), 
        sc, 
        TypeMap<R>(), 
        rbuf,   
        const_cast<int*>(rcs), 
        const_cast<int*>(rds), 
        TypeMap<R>(), 
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
        rbuf, byteRcs.data(), byteRds.data(), MPI_UNSIGNED_CHAR, comm ) );
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
      ( const_cast<Complex<R>*>(sbuf), 2*sc, TypeMap<R>(),
        rbuf, realRcs.data(), realRds.data(), TypeMap<R>(), comm ) );
 #else
    SafeMpi
    ( MPI_Allgatherv
      ( const_cast<Complex<R>*>(sbuf), 
        sc, 
        TypeMap<Complex<R>>(),
        rbuf, 
        const_cast<int*>(rcs), 
        const_cast<int*>(rds), 
        TypeMap<Complex<R>>(),
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
( const unsigned* sbuf, int sc, 
        unsigned* rbuf, const int* rcs, const int* rds, Comm comm );
template void AllGather
( const long int* sbuf, int sc, 
        long int* rbuf, const int* rcs, const int* rds, Comm comm );
template void AllGather
( const unsigned long* sbuf, int sc, 
        unsigned long* rbuf, const int* rcs, const int* rds, Comm comm );
#ifdef HAVE_MPI_LONG_LONG
template void AllGather
( const long long int* sbuf, int sc, 
        long long int* rbuf, const int* rcs, const int* rds, Comm comm );
template void AllGather
( const unsigned long long* sbuf, int sc, 
        unsigned long long* rbuf, const int* rcs, const int* rds, Comm comm );
#endif
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
      ( const_cast<R*>(sbuf), sc, TypeMap<R>(),
        rbuf,                 rc, TypeMap<R>(), root, comm ) );
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
      ( const_cast<Complex<R>*>(sbuf), 2*sc, TypeMap<R>(),
        rbuf,                          2*rc, TypeMap<R>(), root, comm ) );
#else
    SafeMpi
    ( MPI_Scatter
      ( const_cast<Complex<R>*>(sbuf), sc, TypeMap<Complex<R>>(),
        rbuf,                          rc, TypeMap<Complex<R>>(), 
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
( const unsigned* sbuf, int sc, 
        unsigned* rbuf, int rc, int root, Comm comm );
template void Scatter
( const long int* sbuf, int sc, 
        long int* rbuf, int rc, int root, Comm comm );
template void Scatter
( const unsigned long* sbuf, int sc, 
        unsigned long* rbuf, int rc, int root, Comm comm );
#ifdef HAVE_MPI_LONG_LONG
template void Scatter
( const long long int* sbuf, int sc, 
        long long int* rbuf, int rc, int root, Comm comm );
template void Scatter
( const unsigned long long* sbuf, int sc, 
        unsigned long long* rbuf, int rc, int root, Comm comm );
#endif
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
          ( buf,          sc, TypeMap<R>(), 
            MPI_IN_PLACE, rc, TypeMap<R>(), root, comm ) );
#else
        const int commSize = CommSize( comm );
        std::vector<R> sendBuf( sc*commSize );
        MemCopy( sendBuf.data(), buf, sc*commSize );
        SafeMpi
        ( MPI_Scatter
          ( sendBuf.data(), sc, TypeMap<R>(), 
            buf,            rc, TypeMap<R>(), root, comm ) );
#endif
    }
    else
    {
        SafeMpi
        ( MPI_Scatter
          ( 0,   sc, TypeMap<R>(), 
            buf, rc, TypeMap<R>(), root, comm ) );
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
          ( buf,          2*sc, TypeMap<R>(), 
            MPI_IN_PLACE, 2*rc, TypeMap<R>(), root, comm ) );
# else
        const int commSize = CommSize( comm );
        std::vector<Complex<R>> sendBuf( sc*commSize );
        MemCopy( sendBuf.data(), buf, sc*commSize );
        SafeMpi
        ( MPI_Scatter
          ( sendBuf.data(), 2*sc, TypeMap<R>(),          
            buf,            2*rc, TypeMap<R>(), root, comm ) );
# endif
#else
# ifdef HAVE_MPI_IN_PLACE
        SafeMpi
        ( MPI_Scatter
          ( buf,          sc, TypeMap<Complex<R>>(), 
            MPI_IN_PLACE, rc, TypeMap<Complex<R>>(), root, comm ) );
# else
        const int commSize = CommSize( comm );
        std::vector<Complex<R>> sendBuf( sc*commSize );
        MemCopy( sendBuf.data(), buf, sc*commSize );
        SafeMpi
        ( MPI_Scatter
          ( sendBuf.data(), sc, TypeMap<Complex<R>>(),
            buf,            rc, TypeMap<Complex<R>>(), root, comm ) );
# endif
#endif
    }
    else
    {
#ifdef AVOID_COMPLEX_MPI
        SafeMpi
        ( MPI_Scatter
          ( 0,   2*sc, TypeMap<R>(), 
            buf, 2*rc, TypeMap<R>(), root, comm ) );
#else
        SafeMpi
        ( MPI_Scatter
          ( 0,   sc, TypeMap<Complex<R>>(), 
            buf, rc, TypeMap<Complex<R>>(), root, comm ) );
#endif
    }
}

template void Scatter( byte* buf, int sc, int rc, int root, Comm comm );
template void Scatter( int* buf, int sc, int rc, int root, Comm comm );
template void Scatter( unsigned* buf, int sc, int rc, int root, Comm comm );
template void Scatter( long int* buf, int sc, int rc, int root, Comm comm );
template void Scatter( unsigned long* buf, int sc, int rc, int root, Comm comm );
#ifdef HAVE_MPI_LONG_LONG
template void Scatter( long long int* buf, int sc, int rc, int root, Comm comm );
template void Scatter( unsigned long long* buf, int sc, int rc, int root, Comm comm );
#endif
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
      ( const_cast<R*>(sbuf), sc, TypeMap<R>(),
        rbuf,                 rc, TypeMap<R>(), comm ) );
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
      ( const_cast<Complex<R>*>(sbuf), 2*sc, TypeMap<R>(),
        rbuf,                          2*rc, TypeMap<R>(), comm ) );
#else
    SafeMpi
    ( MPI_Alltoall
      ( const_cast<Complex<R>*>(sbuf), sc, TypeMap<Complex<R>>(),
        rbuf,                          rc, TypeMap<Complex<R>>(), comm ) );
#endif
}

template void AllToAll
( const byte* sbuf, int sc, 
        byte* rbuf, int rc, Comm comm );
template void AllToAll
( const int* sbuf, int sc, 
        int* rbuf, int rc, Comm comm );
template void AllToAll
( const unsigned* sbuf, int sc, 
        unsigned* rbuf, int rc, Comm comm );
template void AllToAll
( const long int* sbuf, int sc, 
        long int* rbuf, int rc, Comm comm );
template void AllToAll
( const unsigned long* sbuf, int sc, 
        unsigned long* rbuf, int rc, Comm comm );
#ifdef HAVE_MPI_LONG_LONG
template void AllToAll
( const long long int* sbuf, int sc, 
        long long int* rbuf, int rc, Comm comm );
template void AllToAll
( const unsigned long long* sbuf, int sc, 
        unsigned long long* rbuf, int rc, Comm comm );
#endif
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
        TypeMap<R>(),
        rbuf, 
        const_cast<int*>(rcs), 
        const_cast<int*>(rds), 
        TypeMap<R>(),
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
              scsDoubled.data(), sdsDoubled.data(), TypeMap<R>(),
        rbuf, rcsDoubled.data(), rdsDoubled.data(), TypeMap<R>(), comm ) );
#else
    SafeMpi
    ( MPI_Alltoallv
      ( const_cast<Complex<R>*>(sbuf), 
        const_cast<int*>(scs), 
        const_cast<int*>(sds), 
        TypeMap<Complex<R>>(),
        rbuf, 
        const_cast<int*>(rcs), 
        const_cast<int*>(rds), 
        TypeMap<Complex<R>>(),
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
( const unsigned* sbuf, const int* scs, const int* sds,
        unsigned* rbuf, const int* rcs, const int* rds, Comm comm );
template void AllToAll
( const long int* sbuf, const int* scs, const int* sds,
        long int* rbuf, const int* rcs, const int* rds, Comm comm );
template void AllToAll
( const unsigned long* sbuf, const int* scs, const int* sds,
        unsigned long* rbuf, const int* rcs, const int* rds, Comm comm );
#ifdef HAVE_MPI_LONG_LONG
template void AllToAll
( const long long int* sbuf, const int* scs, const int* sds,
        long long int* rbuf, const int* rcs, const int* rds, Comm comm );
template void AllToAll
( const unsigned long long* sbuf, const int* scs, const int* sds,
        unsigned long long* rbuf, const int* rcs, const int* rds, Comm comm );
#endif
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
          ( const_cast<T*>(sbuf), rbuf, count, TypeMap<T>(),
            op, root, comm ) );
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
                rbuf, 2*count, TypeMap<R>(), op, root, comm ) );
        }
        else
        {
            SafeMpi
            ( MPI_Reduce
              ( const_cast<Complex<R>*>(sbuf),
                rbuf, count, TypeMap<Complex<R>>(), op, root, comm ) );
        }
#else
        SafeMpi
        ( MPI_Reduce
          ( const_cast<Complex<R>*>(sbuf), 
            rbuf, count, TypeMap<Complex<R>>(), op, root, comm ) );
#endif
    }
}

template void Reduce( const byte* sbuf, byte* rbuf, int count, Op op, int root, Comm comm );
template void Reduce( const int* sbuf, int* rbuf, int count, Op op, int root, Comm comm );
template void Reduce( const unsigned* sbuf, unsigned* rbuf, int count, Op op, int root, Comm comm );
template void Reduce( const long int* sbuf, long int* rbuf, int count, Op op, int root, Comm comm );
template void Reduce( const unsigned long* sbuf, unsigned long* rbuf, int count, Op op, int root, Comm comm );
#ifdef HAVE_MPI_LONG_LONG
template void Reduce( const long long int* sbuf, long long int* rbuf, int count, Op op, int root, Comm comm );
template void Reduce( const unsigned long long* sbuf, unsigned long long* rbuf, int count, Op op, int root, Comm comm );
#endif
template void Reduce( const float* sbuf, float* rbuf, int count, Op op, int root, Comm comm );
template void Reduce( const double* sbuf, double* rbuf, int count, Op op, int root, Comm comm );
template void Reduce( const Complex<float>* sbuf, Complex<float>* rbuf, int count, Op op, int root, Comm comm );
template void Reduce( const Complex<double>* sbuf, Complex<double>* rbuf, int count, Op op, int root, Comm comm );
template void Reduce( const ValueInt<Int>* sbuf, ValueInt<Int>* rbuf, int count, Op op, int root, Comm comm );
template void Reduce( const ValueInt<float>* sbuf, ValueInt<float>* rbuf, int count, Op op, int root, Comm comm );
template void Reduce( const ValueInt<double>* sbuf, ValueInt<double>* rbuf, int count, Op op, int root, Comm comm );
template void Reduce( const ValueIntPair<Int>* sbuf, ValueIntPair<Int>* rbuf, int count, Op op, int root, Comm comm );
template void Reduce( const ValueIntPair<float>* sbuf, ValueIntPair<float>* rbuf, int count, Op op, int root, Comm comm );
template void Reduce( const ValueIntPair<double>* sbuf, ValueIntPair<double>* rbuf, int count, Op op, int root, Comm comm );

template<typename T>
void Reduce( const T* sbuf, T* rbuf, int count, int root, Comm comm )
{ Reduce( sbuf, rbuf, count, mpi::SUM, root, comm ); }

template void Reduce( const byte* sbuf, byte* rbuf, int count, int root, Comm comm );
template void Reduce( const int* sbuf, int* rbuf, int count, int root, Comm comm );
template void Reduce( const unsigned* sbuf, unsigned* rbuf, int count, int root, Comm comm );
template void Reduce( const long int* sbuf, long int* rbuf, int count, int root, Comm comm );
template void Reduce( const unsigned long* sbuf, unsigned long* rbuf, int count, int root, Comm comm );
#ifdef HAVE_MPI_LONG_LONG
template void Reduce( const long long int* sbuf, long long int* rbuf, int count, int root, Comm comm );
template void Reduce( const unsigned long long* sbuf, unsigned long long* rbuf, int count, int root, Comm comm );
#endif
template void Reduce( const float* sbuf, float* rbuf, int count, int root, Comm comm );
template void Reduce( const double* sbuf, double* rbuf, int count, int root, Comm comm );
template void Reduce( const Complex<float>* sbuf, Complex<float>* rbuf, int count, int root, Comm comm );
template void Reduce( const Complex<double>* sbuf, Complex<double>* rbuf, int count, int root, Comm comm );
template void Reduce( const ValueInt<Int>* sbuf, ValueInt<Int>* rbuf, int count, int root, Comm comm );
template void Reduce( const ValueInt<float>* sbuf, ValueInt<float>* rbuf, int count, int root, Comm comm );
template void Reduce( const ValueInt<double>* sbuf, ValueInt<double>* rbuf, int count, int root, Comm comm );
template void Reduce( const ValueIntPair<Int>* sbuf, ValueIntPair<Int>* rbuf, int count, int root, Comm comm );
template void Reduce( const ValueIntPair<float>* sbuf, ValueIntPair<float>* rbuf, int count, int root, Comm comm );
template void Reduce( const ValueIntPair<double>* sbuf, ValueIntPair<double>* rbuf, int count, int root, Comm comm );

template<typename T>
T Reduce( T sb, Op op, int root, Comm comm )
{ 
    T rb;
    Reduce( &sb, &rb, 1, op, root, comm );
    return rb;
}

template byte Reduce( byte sb, Op op, int root, Comm comm );
template int Reduce( int sb, Op op, int root, Comm comm );
template unsigned Reduce( unsigned sb, Op op, int root, Comm comm );
template long int Reduce( long int sb, Op op, int root, Comm comm );
template unsigned long Reduce( unsigned long sb, Op op, int root, Comm comm );
#ifdef HAVE_MPI_LONG_LONG
template long long int Reduce( long long int sb, Op op, int root, Comm comm );
template unsigned long long Reduce( unsigned long long sb, Op op, int root, Comm comm );
#endif
template float Reduce( float sb, Op op, int root, Comm comm );
template double Reduce( double sb, Op op, int root, Comm comm );
template Complex<float> Reduce( Complex<float> sb, Op op, int root, Comm comm );
template Complex<double> Reduce( Complex<double> sb, Op op, int root, Comm comm );
template ValueInt<Int> Reduce( ValueInt<Int> sb, Op op, int root, Comm comm );
template ValueInt<float> Reduce( ValueInt<float> sb, Op op, int root, Comm comm );
template ValueInt<double> Reduce( ValueInt<double> sb, Op op, int root, Comm comm );
template ValueIntPair<Int> Reduce( ValueIntPair<Int> sb, Op op, int root, Comm comm );
template ValueIntPair<float> Reduce( ValueIntPair<float> sb, Op op, int root, Comm comm );
template ValueIntPair<double> Reduce( ValueIntPair<double> sb, Op op, int root, Comm comm );

template<typename T>
T Reduce( T sb, int root, Comm comm )
{ 
    T rb;
    Reduce( &sb, &rb, 1, mpi::SUM, root, comm );
    return rb;
}

template byte Reduce( byte sb, int root, Comm comm );
template int Reduce( int sb, int root, Comm comm );
template unsigned Reduce( unsigned sb, int root, Comm comm );
template long int Reduce( long int sb, int root, Comm comm );
template unsigned long Reduce( unsigned long sb, int root, Comm comm );
#ifdef HAVE_MPI_LONG_LONG
template long long int Reduce( long long int sb, int root, Comm comm );
template unsigned long long Reduce( unsigned long long sb, int root, Comm comm );
#endif
template float Reduce( float sb, int root, Comm comm );
template double Reduce( double sb, int root, Comm comm );
template Complex<float> Reduce( Complex<float> sb, int root, Comm comm );
template Complex<double> Reduce( Complex<double> sb, int root, Comm comm );
template ValueInt<Int> Reduce( ValueInt<Int> sb, int root, Comm comm );
template ValueInt<float> Reduce( ValueInt<float> sb, int root, Comm comm );
template ValueInt<double> Reduce( ValueInt<double> sb, int root, Comm comm );
template ValueIntPair<Int> Reduce( ValueIntPair<Int> sb, int root, Comm comm );
template ValueIntPair<float> Reduce( ValueIntPair<float> sb, int root, Comm comm );
template ValueIntPair<double> Reduce( ValueIntPair<double> sb, int root, Comm comm );

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
              ( MPI_IN_PLACE, buf, count, TypeMap<T>(), op, root, comm ) );
#else
            std::vector<T> sendBuf( count );
            MemCopy( sendBuf.data(), buf, count );
            SafeMpi
            ( MPI_Reduce
              ( sendBuf.data(), buf, count, TypeMap<T>(), op, root, comm ) );
#endif
        }
        else
            SafeMpi
            ( MPI_Reduce( buf, 0, count, TypeMap<T>(), op, root, comm ) );
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
                  ( MPI_IN_PLACE, buf, 2*count, TypeMap<R>(), op, 
                    root, comm ) );
# else
                std::vector<Complex<R>> sendBuf( count );
                MemCopy( sendBuf.data(), buf, count );
                SafeMpi
                ( MPI_Reduce
                  ( sendBuf.data(), buf, 2*count, TypeMap<R>(), op, 
                    root, comm ) );
# endif
            }
            else
                SafeMpi
                ( MPI_Reduce
                  ( buf, 0, 2*count, TypeMap<R>(), op, root, comm ) );
        }
        else
        {
            if( commRank == root )
            {
# ifdef HAVE_MPI_IN_PLACE
                SafeMpi
                ( MPI_Reduce
                  ( MPI_IN_PLACE, buf, count, TypeMap<Complex<R>>(), op, 
                    root, comm ) );
# else
                std::vector<Complex<R>> sendBuf( count );
                MemCopy( sendBuf.data(), buf, count );
                SafeMpi
                ( MPI_Reduce
                  ( sendBuf.data(), buf, count, TypeMap<Complex<R>>(), op, 
                    root, comm ) );
# endif
            }
            else
                SafeMpi
                ( MPI_Reduce
                  ( buf, 0, count, TypeMap<Complex<R>>(), op, 
                    root, comm ) );
        }
#else
        if( commRank == root )
        {
# ifdef HAVE_MPI_IN_PLACE
            SafeMpi
            ( MPI_Reduce
              ( MPI_IN_PLACE, buf, count, TypeMap<Complex<R>>(), op, 
                root, comm ) );
# else
            std::vector<Complex<R>> sendBuf( count );
            MemCopy( sendBuf.data(), buf, count );
            SafeMpi
            ( MPI_Reduce
              ( sendBuf.data(), buf, count, TypeMap<Complex<R>>(), op, 
                root, comm ) );
# endif
        }
        else
            SafeMpi
            ( MPI_Reduce
              ( buf, 0, count, TypeMap<Complex<R>>(), op, root, comm ) );
#endif
    }
}

template void Reduce( byte* buf, int count, Op op, int root, Comm comm );
template void Reduce( int* buf, int count, Op op, int root, Comm comm );
template void Reduce( unsigned* buf, int count, Op op, int root, Comm comm );
template void Reduce( long int* buf, int count, Op op, int root, Comm comm );
template void Reduce( unsigned long* buf, int count, Op op, int root, Comm comm );
#ifdef HAVE_MPI_LONG_LONG
template void Reduce( long long int* buf, int count, Op op, int root, Comm comm );
template void Reduce( unsigned long long* buf, int count, Op op, int root, Comm comm );
#endif
template void Reduce( float* buf, int count, Op op, int root, Comm comm );
template void Reduce( double* buf, int count, Op op, int root, Comm comm );
template void Reduce( Complex<float>* buf, int count, Op op, int root, Comm comm );
template void Reduce( Complex<double>* buf, int count, Op op, int root, Comm comm );
template void Reduce( ValueInt<Int>* buf, int count, Op op, int root, Comm comm );
template void Reduce( ValueInt<float>* buf, int count, Op op, int root, Comm comm );
template void Reduce( ValueInt<double>* buf, int count, Op op, int root, Comm comm );
template void Reduce( ValueIntPair<Int>* buf, int count, Op op, int root, Comm comm );
template void Reduce( ValueIntPair<float>* buf, int count, Op op, int root, Comm comm );
template void Reduce( ValueIntPair<double>* buf, int count, Op op, int root, Comm comm );

template<typename T>
void Reduce( T* buf, int count, int root, Comm comm )
{ Reduce( buf, count, mpi::SUM, root, comm ); }

template void Reduce( byte* buf, int count, int root, Comm comm );
template void Reduce( int* buf, int count, int root, Comm comm );
template void Reduce( unsigned* buf, int count, int root, Comm comm );
template void Reduce( long int* buf, int count, int root, Comm comm );
template void Reduce( unsigned long* buf, int count, int root, Comm comm );
#ifdef HAVE_MPI_LONG_LONG
template void Reduce( long long int* buf, int count, int root, Comm comm );
template void Reduce( unsigned long long* buf, int count, int root, Comm comm );
#endif
template void Reduce( float* buf, int count, int root, Comm comm );
template void Reduce( double* buf, int count, int root, Comm comm );
template void Reduce( Complex<float>* buf, int count, int root, Comm comm );
template void Reduce( Complex<double>* buf, int count, int root, Comm comm );
template void Reduce( ValueInt<Int>* buf, int count, int root, Comm comm );
template void Reduce( ValueInt<float>* buf, int count, int root, Comm comm );
template void Reduce( ValueInt<double>* buf, int count, int root, Comm comm );
template void Reduce( ValueIntPair<Int>* buf, int count, int root, Comm comm );
template void Reduce( ValueIntPair<float>* buf, int count, int root, Comm comm );
template void Reduce( ValueIntPair<double>* buf, int count, int root, Comm comm );

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
          ( const_cast<T*>(sbuf), rbuf, count, TypeMap<T>(), op, comm ) );
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
                  rbuf, 2*count, TypeMap<R>(), op, comm ) );
        }
        else
        {
            SafeMpi
            ( MPI_Allreduce
              ( const_cast<Complex<R>*>(sbuf),
                rbuf, count, TypeMap<Complex<R>>(), op, comm ) );
        }
#else
        SafeMpi
        ( MPI_Allreduce
          ( const_cast<Complex<R>*>(sbuf), 
            rbuf, count, TypeMap<Complex<R>>(), op, comm ) );
#endif
    }
}

template void AllReduce( const byte* sbuf, byte* rbuf, int count, Op op, Comm comm );
template void AllReduce( const int* sbuf, int* rbuf, int count, Op op, Comm comm );
template void AllReduce( const unsigned* sbuf, unsigned* rbuf, int count, Op op, Comm comm );
template void AllReduce( const long int* sbuf, long int* rbuf, int count, Op op, Comm comm );
template void AllReduce( const unsigned long* sbuf, unsigned long* rbuf, int count, Op op, Comm comm );
#ifdef HAVE_MPI_LONG_LONG
template void AllReduce( const long long int* sbuf, long long int* rbuf, int count, Op op, Comm comm );
template void AllReduce( const unsigned long long* sbuf, unsigned long long* rbuf, int count, Op op, Comm comm );
#endif
template void AllReduce( const float* sbuf, float* rbuf, int count, Op op, Comm comm );
template void AllReduce( const double* sbuf, double* rbuf, int count, Op op, Comm comm );
template void AllReduce( const Complex<float>* sbuf, Complex<float>* rbuf, int count, Op op, Comm comm );
template void AllReduce( const Complex<double>* sbuf, Complex<double>* rbuf, int count, Op op, Comm comm );
template void AllReduce( const ValueInt<Int>* sbuf, ValueInt<Int>* rbuf, int count, Op op, Comm comm );
template void AllReduce( const ValueInt<float>* sbuf, ValueInt<float>* rbuf, int count, Op op, Comm comm );
template void AllReduce( const ValueInt<double>* sbuf, ValueInt<double>* rbuf, int count, Op op, Comm comm );
template void AllReduce( const ValueIntPair<Int>* sbuf, ValueIntPair<Int>* rbuf, int count, Op op, Comm comm );
template void AllReduce( const ValueIntPair<float>* sbuf, ValueIntPair<float>* rbuf, int count, Op op, Comm comm );
template void AllReduce( const ValueIntPair<double>* sbuf, ValueIntPair<double>* rbuf, int count, Op op, Comm comm );

template<typename T>
void AllReduce( const T* sbuf, T* rbuf, int count, Comm comm )
{ AllReduce( sbuf, rbuf, count, mpi::SUM, comm ); }

template void AllReduce( const byte* sbuf, byte* rbuf, int count, Comm comm );
template void AllReduce( const int* sbuf, int* rbuf, int count, Comm comm );
template void AllReduce( const unsigned* sbuf, unsigned* rbuf, int count, Comm comm );
template void AllReduce( const long int* sbuf, long int* rbuf, int count, Comm comm );
template void AllReduce( const unsigned long* sbuf, unsigned long* rbuf, int count, Comm comm );
#ifdef HAVE_MPI_LONG_LONG
template void AllReduce( const long long int* sbuf, long long int* rbuf, int count, Comm comm );
template void AllReduce( const unsigned long long* sbuf, unsigned long long* rbuf, int count, Comm comm );
#endif
template void AllReduce( const float* sbuf, float* rbuf, int count, Comm comm );
template void AllReduce( const double* sbuf, double* rbuf, int count, Comm comm );
template void AllReduce( const Complex<float>* sbuf, Complex<float>* rbuf, int count, Comm comm );
template void AllReduce( const Complex<double>* sbuf, Complex<double>* rbuf, int count, Comm comm );
template void AllReduce( const ValueInt<Int>* sbuf, ValueInt<Int>* rbuf, int count, Comm comm );
template void AllReduce( const ValueInt<float>* sbuf, ValueInt<float>* rbuf, int count, Comm comm );
template void AllReduce( const ValueInt<double>* sbuf, ValueInt<double>* rbuf, int count, Comm comm );
template void AllReduce( const ValueIntPair<Int>* sbuf, ValueIntPair<Int>* rbuf, int count, Comm comm );
template void AllReduce( const ValueIntPair<float>* sbuf, ValueIntPair<float>* rbuf, int count, Comm comm );
template void AllReduce( const ValueIntPair<double>* sbuf, ValueIntPair<double>* rbuf, int count, Comm comm );

template<typename T>
T AllReduce( T sb, Op op, Comm comm )
{ T rb; AllReduce( &sb, &rb, 1, op, comm ); return rb; }

template byte AllReduce( byte sb, Op op, Comm comm );
template int AllReduce( int sb, Op op, Comm comm );
template unsigned AllReduce( unsigned sb, Op op, Comm comm );
template long int AllReduce( long int sb, Op op, Comm comm );
template unsigned long AllReduce( unsigned long sb, Op op, Comm comm );
#ifdef HAVE_MPI_LONG_LONG
template long long int AllReduce( long long int sb, Op op, Comm comm );
template unsigned long long AllReduce( unsigned long long sb, Op op, Comm comm );
#endif
template float AllReduce( float sb, Op op, Comm comm );
template double AllReduce( double sb, Op op, Comm comm );
template Complex<float> AllReduce( Complex<float> sb, Op op, Comm comm );
template Complex<double> AllReduce( Complex<double> sb, Op op, Comm comm );
template ValueInt<Int> AllReduce( ValueInt<Int> sb, Op op, Comm comm );
template ValueInt<float> AllReduce( ValueInt<float> sb, Op op, Comm comm );
template ValueInt<double> AllReduce( ValueInt<double> sb, Op op, Comm comm );
template ValueIntPair<Int> AllReduce( ValueIntPair<Int> sb, Op op, Comm comm );
template ValueIntPair<float> AllReduce( ValueIntPair<float> sb, Op op, Comm comm );
template ValueIntPair<double> AllReduce( ValueIntPair<double> sb, Op op, Comm comm );

template<typename T>
T AllReduce( T sb, Comm comm )
{ return AllReduce( sb, mpi::SUM, comm ); }

template byte AllReduce( byte sb, Comm comm );
template int AllReduce( int sb, Comm comm );
template unsigned AllReduce( unsigned sb, Comm comm );
template long int AllReduce( long int sb, Comm comm );
template unsigned long AllReduce( unsigned long sb, Comm comm );
#ifdef HAVE_MPI_LONG_LONG
template long long int AllReduce( long long int sb, Comm comm );
template unsigned long long AllReduce( unsigned long long sb, Comm comm );
#endif
template float AllReduce( float sb, Comm comm );
template double AllReduce( double sb, Comm comm );
template Complex<float> AllReduce( Complex<float> sb, Comm comm );
template Complex<double> AllReduce( Complex<double> sb, Comm comm );
template ValueInt<Int> AllReduce( ValueInt<Int> sb, Comm comm );
template ValueInt<float> AllReduce( ValueInt<float> sb, Comm comm );
template ValueInt<double> AllReduce( ValueInt<double> sb, Comm comm );
template ValueIntPair<Int> AllReduce( ValueIntPair<Int> sb, Comm comm );
template ValueIntPair<float> AllReduce( ValueIntPair<float> sb, Comm comm );
template ValueIntPair<double> AllReduce( ValueIntPair<double> sb, Comm comm );

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
          ( MPI_IN_PLACE, buf, count, TypeMap<T>(), op, comm ) );
#else
        std::vector<T> sendBuf( count );
        MemCopy( sendBuf.data(), buf, count );
        SafeMpi
        ( MPI_Allreduce( sendBuf.data(), buf, count, TypeMap<T>(), op, comm ) );
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
              ( MPI_IN_PLACE, buf, 2*count, TypeMap<R>(), op, comm ) );
# else
            std::vector<Complex<R>> sendBuf( count );
            MemCopy( sendBuf.data(), buf, count );
            SafeMpi
            ( MPI_Allreduce
              ( sendBuf.data(), buf, 2*count, TypeMap<R>(), op, comm ) );
# endif
        }
        else
        {
# ifdef HAVE_MPI_IN_PLACE
            SafeMpi
            ( MPI_Allreduce
              ( MPI_IN_PLACE, buf, count, TypeMap<Complex<R>>(), 
                op, comm ) );
# else
            std::vector<Complex<R>> sendBuf( count );
            MemCopy( sendBuf.data(), buf, count );
            SafeMpi
            ( MPI_Allreduce
              ( sendBuf.data(), buf, count, TypeMap<Complex<R>>(), 
                op, comm ) );
# endif
        }
#else
# ifdef HAVE_MPI_IN_PLACE
        SafeMpi
        ( MPI_Allreduce
          ( MPI_IN_PLACE, buf, count, TypeMap<Complex<R>>(), op, comm ) );
# else
        std::vector<Complex<R>> sendBuf( count );
        MemCopy( sendBuf.data(), buf, count );
        SafeMpi
        ( MPI_Allreduce
          ( sendBuf.data(), buf, count, TypeMap<Complex<R>>(), op, comm ) );
# endif
#endif
    }
}

template void AllReduce( byte* buf, int count, Op op, Comm comm );
template void AllReduce( int* buf, int count, Op op, Comm comm );
template void AllReduce( unsigned* buf, int count, Op op, Comm comm );
template void AllReduce( long int* buf, int count, Op op, Comm comm );
template void AllReduce( unsigned long* buf, int count, Op op, Comm comm );
#ifdef HAVE_MPI_LONG_LONG
template void AllReduce( long long int* buf, int count, Op op, Comm comm );
template void AllReduce( unsigned long long* buf, int count, Op op, Comm comm );
#endif
template void AllReduce( float* buf, int count, Op op, Comm comm );
template void AllReduce( double* buf, int count, Op op, Comm comm );
template void AllReduce( Complex<float>* buf, int count, Op op, Comm comm );
template void AllReduce( Complex<double>* buf, int count, Op op, Comm comm );
template void AllReduce( ValueInt<Int>* buf, int count, Op op, Comm comm );
template void AllReduce( ValueInt<float>* buf, int count, Op op, Comm comm );
template void AllReduce( ValueInt<double>* buf, int count, Op op, Comm comm );
template void AllReduce( ValueIntPair<Int>* buf, int count, Op op, Comm comm );
template void AllReduce( ValueIntPair<float>* buf, int count, Op op, Comm comm );
template void AllReduce( ValueIntPair<double>* buf, int count, Op op, Comm comm );

template<typename T>
void AllReduce( T* buf, int count, Comm comm )
{ AllReduce( buf, count, mpi::SUM, comm ); }

template void AllReduce( byte* buf, int count, Comm comm );
template void AllReduce( int* buf, int count, Comm comm );
template void AllReduce( unsigned* buf, int count, Comm comm );
template void AllReduce( long int* buf, int count, Comm comm );
template void AllReduce( unsigned long* buf, int count, Comm comm );
#ifdef HAVE_MPI_LONG_LONG
template void AllReduce( long long int* buf, int count, Comm comm );
template void AllReduce( unsigned long long* buf, int count, Comm comm );
#endif
template void AllReduce( float* buf, int count, Comm comm );
template void AllReduce( double* buf, int count, Comm comm );
template void AllReduce( Complex<float>* buf, int count, Comm comm );
template void AllReduce( Complex<double>* buf, int count, Comm comm );
template void AllReduce( ValueInt<Int>* buf, int count, Comm comm );
template void AllReduce( ValueInt<float>* buf, int count, Comm comm );
template void AllReduce( ValueInt<double>* buf, int count, Comm comm );
template void AllReduce( ValueIntPair<Int>* buf, int count, Comm comm );
template void AllReduce( ValueIntPair<float>* buf, int count, Comm comm );
template void AllReduce( ValueIntPair<double>* buf, int count, Comm comm );

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
    ( MPI_Reduce_scatter_block( sbuf, rbuf, rc, TypeMap<R>(), op, comm ) );
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
    ( MPI_Reduce_scatter_block( sbuf, rbuf, 2*rc, TypeMap<R>(), op, comm ) );
# else
    SafeMpi
    ( MPI_Reduce_scatter_block
      ( sbuf, rbuf, rc, TypeMap<Complex<R>>(), op, comm ) );
# endif
#else
    const int commSize = CommSize( comm );
    Reduce( sbuf, rc*commSize, op, 0, comm );
    Scatter( sbuf, rc, rbuf, rc, 0, comm );
#endif
}

template void ReduceScatter( byte* sbuf, byte* rbuf, int rc, Op op, Comm comm );
template void ReduceScatter( int* sbuf, int* rbuf, int rc, Op op, Comm comm );
template void ReduceScatter( unsigned* sbuf, unsigned* rbuf, int rc, Op op, Comm comm );
template void ReduceScatter( long int* sbuf, long int* rbuf, int rc, Op op, Comm comm );
template void ReduceScatter( unsigned long* sbuf, unsigned long* rbuf, int rc, Op op, Comm comm );
#ifdef HAVE_MPI_LONG_LONG
template void ReduceScatter( long long int* sbuf, long long int* rbuf, int rc, Op op, Comm comm );
template void ReduceScatter( unsigned long long* sbuf, unsigned long long* rbuf, int rc, Op op, Comm comm );
#endif
template void ReduceScatter( float* sbuf, float* rbuf, int rc, Op op, Comm comm );
template void ReduceScatter( double* sbuf, double* rbuf, int rc, Op op, Comm comm );
template void ReduceScatter( Complex<float>* sbuf, Complex<float>* rbuf, int rc, Op op, Comm comm );
template void ReduceScatter( Complex<double>* sbuf, Complex<double>* rbuf, int rc, Op op, Comm comm );

template<typename T>
void ReduceScatter( T* sbuf, T* rbuf, int rc, Comm comm )
{ ReduceScatter( sbuf, rbuf, rc, mpi::SUM, comm ); }

template void ReduceScatter( byte* sbuf, byte* rbuf, int rc, Comm comm );
template void ReduceScatter( int* sbuf, int* rbuf, int rc, Comm comm );
template void ReduceScatter( unsigned* sbuf, unsigned* rbuf, int rc, Comm comm );
template void ReduceScatter( long int* sbuf, long int* rbuf, int rc, Comm comm );
template void ReduceScatter( unsigned long* sbuf, unsigned long* rbuf, int rc, Comm comm );
#ifdef HAVE_MPI_LONG_LONG
template void ReduceScatter( long long int* sbuf, long long int* rbuf, int rc, Comm comm );
template void ReduceScatter( unsigned long long* sbuf, unsigned long long* rbuf, int rc, Comm comm );
#endif
template void ReduceScatter( float* sbuf, float* rbuf, int rc, Comm comm );
template void ReduceScatter( double* sbuf, double* rbuf, int rc, Comm comm );
template void ReduceScatter( Complex<float>* sbuf, Complex<float>* rbuf, int rc, Comm comm );
template void ReduceScatter( Complex<double>* sbuf, Complex<double>* rbuf, int rc, Comm comm );

template<typename T>
T ReduceScatter( T sb, Op op, Comm comm )
{ T rb; ReduceScatter( &sb, &rb, 1, op, comm ); return rb; }

template byte ReduceScatter( byte sb, Op op, Comm comm );
template int ReduceScatter( int sb, Op op, Comm comm );
template unsigned ReduceScatter( unsigned sb, Op op, Comm comm );
template long int ReduceScatter( long int sb, Op op, Comm comm );
template unsigned long ReduceScatter( unsigned long sb, Op op, Comm comm );
#ifdef HAVE_MPI_LONG_LONG
template long long int ReduceScatter( long long int sb, Op op, Comm comm );
template unsigned long long ReduceScatter( unsigned long long sb, Op op, Comm comm );
#endif
template float ReduceScatter( float sb, Op op, Comm comm );
template double ReduceScatter( double sb, Op op, Comm comm );
template Complex<float> ReduceScatter( Complex<float> sb, Op op, Comm comm );
template Complex<double> ReduceScatter( Complex<double> sb, Op op, Comm comm );

template<typename T>
T ReduceScatter( T sb, Comm comm )
{ return ReduceScatter( sb, mpi::SUM, comm ); }

template byte ReduceScatter( byte sb, Comm comm );
template int ReduceScatter( int sb, Comm comm );
template unsigned ReduceScatter( unsigned sb, Comm comm );
template long int ReduceScatter( long int sb, Comm comm );
template unsigned long ReduceScatter( unsigned long sb, Comm comm );
#ifdef HAVE_MPI_LONG_LONG
template long long int ReduceScatter( long long int sb, Comm comm );
template unsigned long long ReduceScatter( unsigned long long sb, Comm comm );
#endif
template float ReduceScatter( float sb, Comm comm );
template double ReduceScatter( double sb, Comm comm );
template Complex<float> ReduceScatter( Complex<float> sb, Comm comm );
template Complex<double> ReduceScatter( Complex<double> sb, Comm comm );

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
      ( MPI_IN_PLACE, buf, rc, TypeMap<R>(), op, comm ) );
# else
    const int commSize = CommSize( comm );
    std::vector<R> sendBuf( rc*commSize );
    MemCopy( sendBuf.data(), buf, rc*commSize );
    SafeMpi
    ( MPI_Reduce_scatter_block
      ( sendBuf.data(), buf, rc, TypeMap<R>(), op, comm ) );
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
      ( MPI_IN_PLACE, buf, 2*rc, TypeMap<R>(), op, comm ) );
#  else 
    const int commSize = CommSize( comm );
    std::vector<Complex<R>> sendBuf( rc*commSize );
    MemCopy( sendBuf.data(), buf, rc*commSize );
    SafeMpi
    ( MPI_Reduce_scatter_block
      ( sendBuf.data(), buf, 2*rc, TypeMap<R>(), op, comm ) );
#  endif
# else
#  ifdef HAVE_MPI_IN_PLACE
    SafeMpi
    ( MPI_Reduce_scatter_block
      ( MPI_IN_PLACE, buf, rc, TypeMap<Complex<R>>(), op, comm ) );
#  else
    const int commSize = CommSize( comm );
    std::vector<Complex<R>> sendBuf( rc*commSize );
    MemCopy( sendBuf.data(), buf, rc*commSize );
    SafeMpi
    ( MPI_Reduce_scatter_block
      ( sendBuf.data(), buf, rc, TypeMap<Complex<R>>(), op, comm ) );
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
template void ReduceScatter( unsigned* buf, int rc, Op op, Comm comm );
template void ReduceScatter( long int* buf, int rc, Op op, Comm comm );
template void ReduceScatter( unsigned long* buf, int rc, Op op, Comm comm );
#ifdef HAVE_MPI_LONG_LONG
template void ReduceScatter( long long int* buf, int rc, Op op, Comm comm );
template void ReduceScatter( unsigned long long* buf, int rc, Op op, Comm comm );
#endif
template void ReduceScatter( float* buf, int rc, Op op, Comm comm );
template void ReduceScatter( double* buf, int rc, Op op, Comm comm );
template void ReduceScatter( Complex<float>* buf, int rc, Op op, Comm comm );
template void ReduceScatter( Complex<double>* buf, int rc, Op op, Comm comm );

template<typename T>
void ReduceScatter( T* buf, int rc, Comm comm )
{ ReduceScatter( buf, rc, mpi::SUM, comm ); }

template void ReduceScatter( byte* buf, int rc, Comm comm );
template void ReduceScatter( int* buf, int rc, Comm comm );
template void ReduceScatter( unsigned* buf, int rc, Comm comm );
template void ReduceScatter( long int* buf, int rc, Comm comm );
template void ReduceScatter( unsigned long* buf, int rc, Comm comm );
#ifdef HAVE_MPI_LONG_LONG
template void ReduceScatter( long long int* buf, int rc, Comm comm );
template void ReduceScatter( unsigned long long* buf, int rc, Comm comm );
#endif
template void ReduceScatter( float* buf, int rc, Comm comm );
template void ReduceScatter( double* buf, int rc, Comm comm );
template void ReduceScatter( Complex<float>* buf, int rc, Comm comm );
template void ReduceScatter( Complex<double>* buf, int rc, Comm comm );

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
        rbuf, const_cast<int*>(rcs), TypeMap<R>(), op, comm ) );
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
            rbuf, rcsDoubled.data(), TypeMap<R>(), op, comm ) );
    }
    else
    {
        SafeMpi
        ( MPI_Reduce_scatter
          ( const_cast<Complex<R>*>(sbuf),
            rbuf, const_cast<int*>(rcs), TypeMap<Complex<R>>(), 
            op, comm ) );
    }
#else
    SafeMpi
    ( MPI_Reduce_scatter
      ( const_cast<Complex<R>*>(sbuf), 
        rbuf, const_cast<int*>(rcs), TypeMap<Complex<R>>(), op, comm ) );
#endif
}

template void ReduceScatter( const byte* sbuf, byte* rbuf, const int* rcs, Op op, Comm comm );
template void ReduceScatter( const int* sbuf, int* rbuf, const int* rcs, Op op, Comm comm );
template void ReduceScatter( const unsigned* sbuf, unsigned* rbuf, const int* rcs, Op op, Comm comm );
template void ReduceScatter( const long int* sbuf, long int* rbuf, const int* rcs, Op op, Comm comm );
template void ReduceScatter( const unsigned long* sbuf, unsigned long* rbuf, const int* rcs, Op op, Comm comm );
#ifdef HAVE_MPI_LONG_LONG
template void ReduceScatter( const long long int* sbuf, long long int* rbuf, const int* rcs, Op op, Comm comm );
template void ReduceScatter( const unsigned long long* sbuf, unsigned long long* rbuf, const int* rcs, Op op, Comm comm );
#endif
template void ReduceScatter( const float* sbuf, float* rbuf, const int* rcs, Op op, Comm comm );
template void ReduceScatter( const double* sbuf, double* rbuf, const int* rcs, Op op, Comm comm );
template void ReduceScatter( const Complex<float>* sbuf, Complex<float>* rbuf, const int* rcs, Op op, Comm comm );
template void ReduceScatter( const Complex<double>* sbuf, Complex<double>* rbuf, const int* rcs, Op op, Comm comm );

template<typename T>
void ReduceScatter( const T* sbuf, T* rbuf, const int* rcs, Comm comm )
{ ReduceScatter( sbuf, rbuf, rcs, mpi::SUM, comm ); }

template void ReduceScatter( const byte* sbuf, byte* rbuf, const int* rcs, Comm comm );
template void ReduceScatter( const int* sbuf, int* rbuf, const int* rcs, Comm comm );
template void ReduceScatter( const unsigned* sbuf, unsigned* rbuf, const int* rcs, Comm comm );
template void ReduceScatter( const long int* sbuf, long int* rbuf, const int* rcs, Comm comm );
template void ReduceScatter( const unsigned long* sbuf, unsigned long* rbuf, const int* rcs, Comm comm );
#ifdef HAVE_MPI_LONG_LONG
template void ReduceScatter( const long long int* sbuf, long long int* rbuf, const int* rcs, Comm comm );
template void ReduceScatter( const unsigned long long* sbuf, unsigned long long* rbuf, const int* rcs, Comm comm );
#endif
template void ReduceScatter( const float* sbuf, float* rbuf, const int* rcs, Comm comm );
template void ReduceScatter( const double* sbuf, double* rbuf, const int* rcs, Comm comm );
template void ReduceScatter( const Complex<float>* sbuf, Complex<float>* rbuf, const int* rcs, Comm comm );
template void ReduceScatter( const Complex<double>* sbuf, Complex<double>* rbuf, const int* rcs, Comm comm );

} // namespace mpi
} // namespace elem
