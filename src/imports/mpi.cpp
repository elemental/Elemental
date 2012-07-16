/*
   Copyright (c) 2009-2012, Jack Poulson
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
#include "elemental/core/environment.hpp"

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
    SafeMpi( MPI_Comm_rank( mpi::COMM_WORLD, &rank ) );
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
    SafeMpi( MPI_Errhandler_set( comm, errorHandler ) );
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

template<>
int GetCount<byte>( Status& status )
{
#ifndef RELEASE
    PushCallStack("mpi::GetCount");
#endif
    int count;
    SafeMpi( MPI_Get_count( &status, MPI_UNSIGNED_CHAR, &count ) );
#ifndef RELEASE
    PopCallStack();
#endif
    return count;
}

template<>
int GetCount<int>( Status& status )
{
#ifndef RELEASE
    PushCallStack("mpi::GetCount");
#endif
    int count;
    SafeMpi( MPI_Get_count( &status, MPI_INT, &count ) );
#ifndef RELEASE
    PopCallStack();
#endif
    return count;
}

template<>
int GetCount<float>( Status& status )
{
#ifndef RELEASE
    PushCallStack("mpi::GetCount");
#endif
    int count;
    SafeMpi( MPI_Get_count( &status, MPI_FLOAT, &count ) );
#ifndef RELEASE
    PopCallStack();
#endif
    return count;
}

template<>
int GetCount<double>( Status& status )
{
#ifndef RELEASE
    PushCallStack("mpi::GetCount");
#endif
    int count;
    SafeMpi( MPI_Get_count( &status, MPI_DOUBLE, &count ) );
#ifndef RELEASE
    PopCallStack();
#endif
    return count;
}

template<>
int GetCount<scomplex>( Status& status )
{
#ifndef RELEASE
    PushCallStack("mpi::GetCount");
#endif
    int count;
    SafeMpi( MPI_Get_count( &status, MPI_COMPLEX, &count ) );
#ifndef RELEASE
    PopCallStack();
#endif
    return count;
}

template<>
int GetCount<dcomplex>( Status& status )
{
#ifndef RELEASE
    PushCallStack("mpi::GetCount");
#endif
    int count;
    SafeMpi( MPI_Get_count( &status, MPI_DOUBLE_COMPLEX, &count ) );
#ifndef RELEASE
    PopCallStack();
#endif
    return count;
}

void Send( const byte* buf, int count, int to, int tag, Comm comm )
{ 
#ifndef RELEASE
    PushCallStack("mpi::Send");
#endif
    SafeMpi( 
        MPI_Send
        ( const_cast<byte*>(buf), count, MPI_UNSIGNED_CHAR, to, tag, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void ISend
( const byte* buf, int count, int to, int tag, Comm comm, Request& request )
{ 
#ifndef RELEASE
    PushCallStack("mpi::ISend");
#endif
    SafeMpi( 
        MPI_Isend
        ( const_cast<byte*>(buf), count, MPI_UNSIGNED_CHAR, to, tag, comm,
          &request ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void ISSend
( const byte* buf, int count, int to, int tag, Comm comm, Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::ISSend");
#endif
    SafeMpi(
        MPI_Issend
        ( const_cast<byte*>(buf), count, MPI_UNSIGNED_CHAR, to, tag, comm, 
          &request )
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void Send( const int* buf, int count, int to, int tag, Comm comm )
{ 
#ifndef RELEASE
    PushCallStack("mpi::Send");
#endif
    SafeMpi( 
        MPI_Send( const_cast<int*>(buf), count, MPI_INT, to, tag, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void ISend
( const int* buf, int count, int to, int tag, Comm comm, Request& request )
{ 
#ifndef RELEASE
    PushCallStack("mpi::ISend");
#endif
    SafeMpi( 
        MPI_Isend
        ( const_cast<int*>(buf), count, MPI_INT, to, tag, comm, &request ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void ISSend
( const int* buf, int count, int to, int tag, Comm comm, Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::ISSend");
#endif
    SafeMpi(
        MPI_Issend
        ( const_cast<int*>(buf), count, MPI_INT, to, tag, comm, &request )
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void Send( const float* buf, int count, int to, int tag, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Send");
#endif
    SafeMpi( 
        MPI_Send( const_cast<float*>(buf), count, MPI_FLOAT, to, tag, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void ISend
( const float* buf, int count, int to, int tag, Comm comm, Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::ISend");
#endif
    SafeMpi( 
        MPI_Isend
        ( const_cast<float*>(buf), count, MPI_FLOAT, to, tag, comm, &request ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void ISSend
( const float* buf, int count, int to, int tag, Comm comm, Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::ISSend");
#endif
    SafeMpi(
        MPI_Issend
        ( const_cast<float*>(buf), count, MPI_FLOAT, to, tag, comm, &request )
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void Send( const double* buf, int count, int to, int tag, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Send");
#endif
    SafeMpi( 
        MPI_Send( const_cast<double*>(buf), count, MPI_DOUBLE, to, tag, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void ISend
( const double* buf, int count, int to, int tag, Comm comm, Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::ISend");
#endif
    SafeMpi( 
        MPI_Isend
        ( const_cast<double*>(buf), count, MPI_DOUBLE, to, tag, comm, &request )
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void ISSend
( const double* buf, int count, int to, int tag, Comm comm, Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::ISSend");
#endif
    SafeMpi(
        MPI_Issend
        ( const_cast<double*>(buf), count, MPI_DOUBLE, to, tag, comm, &request )
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void Send( const scomplex* buf, int count, int to, int tag, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Send");
#endif
#ifdef AVOID_COMPLEX_MPI
    SafeMpi(
        MPI_Send
        ( const_cast<scomplex*>(buf), 2*count, MPI_FLOAT, to, tag, comm )
    );
#else
    SafeMpi( 
        MPI_Send
        ( const_cast<scomplex*>(buf), count, MPI_COMPLEX, to, tag, comm )
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void ISend
( const scomplex* buf, int count, int to, int tag, Comm comm, Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::ISend");
#endif
#ifdef AVOID_COMPLEX_MPI
    SafeMpi(
        MPI_Isend
        ( const_cast<scomplex*>(buf), 2*count, MPI_FLOAT, to, tag, comm,
          &request )
    );
#else
    SafeMpi( 
        MPI_Isend
        ( const_cast<scomplex*>(buf), count, MPI_COMPLEX, to, tag, comm,
          &request )
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void ISSend
( const scomplex* buf, int count, int to, int tag, Comm comm, Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::ISSend");
#endif
#ifdef AVOID_COMPLEX_MPI
    SafeMpi(
        MPI_Issend
        ( const_cast<scomplex*>(buf), 2*count, MPI_FLOAT, to, tag, comm,
          &request )
    );
#else
    SafeMpi(
        MPI_Issend
        ( const_cast<scomplex*>(buf), count, MPI_COMPLEX, to, tag, comm,
          &request )
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void Send( const dcomplex* buf, int count, int to, int tag, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Send");
#endif
#ifdef AVOID_COMPLEX_MPI
    SafeMpi(
        MPI_Send
        ( const_cast<dcomplex*>(buf), 2*count, MPI_DOUBLE, to, tag, comm )
    );
#else
    SafeMpi( 
        MPI_Send
        ( const_cast<dcomplex*>(buf), count, MPI_DOUBLE_COMPLEX, to, tag, comm )
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void ISend
( const dcomplex* buf, int count, int to, int tag, Comm comm, Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::ISend");
#endif
#ifdef AVOID_COMPLEX_MPI
    SafeMpi(
        MPI_Isend
        ( const_cast<dcomplex*>(buf), 2*count, MPI_DOUBLE, to, tag, comm,
          &request )
    );
#else
    SafeMpi( 
        MPI_Isend
        ( const_cast<dcomplex*>(buf), count, MPI_DOUBLE_COMPLEX, to, tag, comm,
          &request )
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void ISSend
( const dcomplex* buf, int count, int to, int tag, Comm comm, Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::ISSend");
#endif
#ifdef AVOID_COMPLEX_MPI
    SafeMpi(
        MPI_Issend
        ( const_cast<dcomplex*>(buf), 2*count, MPI_DOUBLE, to, tag, comm,
          &request )
    );
#else
    SafeMpi(
        MPI_Issend
        ( const_cast<dcomplex*>(buf), count, MPI_DOUBLE_COMPLEX, to, tag, comm,
          &request )
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void Recv( byte* buf, int count, int from, int tag, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Recv");
#endif
    Status status;
    SafeMpi( 
        MPI_Recv( buf, count, MPI_UNSIGNED_CHAR, from, tag, comm, &status ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void IRecv
( byte* buf, int count, int from, int tag, Comm comm, Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::IRecv");
#endif
    SafeMpi( 
        MPI_Irecv( buf, count, MPI_UNSIGNED_CHAR, from, tag, comm, &request ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void Recv( int* buf, int count, int from, int tag, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Recv");
#endif
    Status status;
    SafeMpi( MPI_Recv( buf, count, MPI_INT, from, tag, comm, &status ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

void IRecv
( int* buf, int count, int from, int tag, Comm comm, Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::IRecv");
#endif
    SafeMpi( MPI_Irecv( buf, count, MPI_INT, from, tag, comm, &request ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

void Recv( float* buf, int count, int from, int tag, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Recv");
#endif
    Status status;
    SafeMpi( MPI_Recv( buf, count, MPI_FLOAT, from, tag, comm, &status ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

void IRecv
( float* buf, int count, int from, int tag, Comm comm, Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::IRecv");
#endif
    SafeMpi( MPI_Irecv( buf, count, MPI_FLOAT, from, tag, comm, &request ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

void Recv( double* buf, int count, int from, int tag, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Recv");
#endif
    Status status;
    SafeMpi( MPI_Recv( buf, count, MPI_DOUBLE, from, tag, comm, &status ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

void IRecv
( double* buf, int count, int from, int tag, Comm comm, Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::IRecv");
#endif
    SafeMpi( MPI_Irecv( buf, count, MPI_DOUBLE, from, tag, comm, &request ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

void Recv( scomplex* buf, int count, int from, int tag, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Recv");
#endif
    Status status;
#ifdef AVOID_COMPLEX_MPI
    SafeMpi( MPI_Recv( buf, 2*count, MPI_FLOAT, from, tag, comm, &status ) );
#else
    SafeMpi( MPI_Recv( buf, count, MPI_COMPLEX, from, tag, comm, &status ) );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void IRecv
( scomplex* buf, int count, int from, int tag, Comm comm, Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::IRecv");
#endif
#ifdef AVOID_COMPLEX_MPI
    SafeMpi( MPI_Irecv( buf, 2*count, MPI_FLOAT, from, tag, comm, &request ) );
#else
    SafeMpi( MPI_Irecv( buf, count, MPI_COMPLEX, from, tag, comm, &request ) );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void Recv( dcomplex* buf, int count, int from, int tag, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Recv");
#endif
    Status status;
#ifdef AVOID_COMPLEX_MPI
    SafeMpi(
        MPI_Recv( buf, 2*count, MPI_DOUBLE, from, tag, comm, &status )
    );
#else
    SafeMpi( 
        MPI_Recv( buf, count, MPI_DOUBLE_COMPLEX, from, tag, comm, &status )
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void IRecv
( dcomplex* buf, int count, int from, int tag, Comm comm, Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::IRecv");
#endif
#ifdef AVOID_COMPLEX_MPI
    SafeMpi(
        MPI_Irecv( buf, 2*count, MPI_DOUBLE, from, tag, comm, &request )
    );
#else
    SafeMpi( 
        MPI_Irecv( buf, count, MPI_DOUBLE_COMPLEX, from, tag, comm, &request )
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void SendRecv
( const byte* sbuf, int sc, int to,   int stag,
        byte* rbuf, int rc, int from, int rtag, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::SendRecv");
#endif
    Status status;
    SafeMpi( 
        MPI_Sendrecv
        ( const_cast<byte*>(sbuf), sc, MPI_UNSIGNED_CHAR, to, stag,
          rbuf, rc, MPI_UNSIGNED_CHAR, from, rtag, comm, &status ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void SendRecv
( const int* sbuf, int sc, int to,   int stag,
        int* rbuf, int rc, int from, int rtag, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::SendRecv");
#endif
    Status status;
    SafeMpi( 
        MPI_Sendrecv
        ( const_cast<int*>(sbuf), sc, MPI_INT, to, stag,
          rbuf, rc, MPI_INT, from, rtag, comm, &status ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void SendRecv
( const float* sbuf, int sc, int to,   int stag,
        float* rbuf, int rc, int from, int rtag, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::SendRecv");
#endif
    Status status;
    SafeMpi( 
        MPI_Sendrecv
        ( const_cast<float*>(sbuf), sc, MPI_FLOAT, to, stag,
          rbuf, rc, MPI_FLOAT, from, rtag, comm, &status ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void SendRecv
( const double* sbuf, int sc, int to,   int stag,
        double* rbuf, int rc, int from, int rtag, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::SendRecv");
#endif
    Status status;
    SafeMpi( 
        MPI_Sendrecv
        ( const_cast<double*>(sbuf), sc, MPI_DOUBLE, to, stag,
          rbuf, rc, MPI_DOUBLE, from, rtag, comm, &status ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void SendRecv
( const scomplex* sbuf, int sc, int to,   int stag,
        scomplex* rbuf, int rc, int from, int rtag, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::SendRecv");
#endif
    Status status;
#ifdef AVOID_COMPLEX_MPI
    SafeMpi(
        MPI_Sendrecv
        ( const_cast<scomplex*>(sbuf), 2*sc, MPI_FLOAT, to,   stag,
          rbuf,                        2*rc, MPI_FLOAT, from, rtag, 
          comm, &status )
    );
#else
    SafeMpi( 
        MPI_Sendrecv
        ( const_cast<scomplex*>(sbuf), sc, MPI_COMPLEX, to,   stag,
          rbuf,                        rc, MPI_COMPLEX, from, rtag, 
          comm, &status )
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void SendRecv
( const dcomplex* sbuf, int sc, int to,   int stag,
        dcomplex* rbuf, int rc, int from, int rtag, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::SendRecv");
#endif
    Status status;
#ifdef AVOID_COMPLEX_MPI
    SafeMpi(
        MPI_Sendrecv
        ( const_cast<dcomplex*>(sbuf), 2*sc, MPI_DOUBLE, to,   stag,
          rbuf,                        2*rc, MPI_DOUBLE, from, rtag, 
          comm, &status )
    );
#else
    SafeMpi( 
        MPI_Sendrecv
        ( const_cast<dcomplex*>(sbuf), sc, MPI_DOUBLE_COMPLEX, to,   stag,
          rbuf,                        rc, MPI_DOUBLE_COMPLEX, from, rtag, 
          comm, &status )
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void Broadcast( byte* buf, int count, int root, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Broadcast");
#endif
    SafeMpi( MPI_Bcast( buf, count, MPI_UNSIGNED_CHAR, root, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

void Broadcast( int* buf, int count, int root, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Broadcast");
#endif
    SafeMpi( MPI_Bcast( buf, count, MPI_INT, root, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

void Broadcast( float* buf, int count, int root, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Broadcast");
#endif
    SafeMpi( MPI_Bcast( buf, count, MPI_FLOAT, root, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

void Broadcast( double* buf, int count, int root, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Broadcast");
#endif
    SafeMpi( MPI_Bcast( buf, count, MPI_DOUBLE, root, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

void Broadcast( scomplex* buf, int count, int root, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Broadcast");
#endif
#ifdef AVOID_COMPLEX_MPI
    SafeMpi( MPI_Bcast( buf, 2*count, MPI_FLOAT, root, comm ) );
#else
    SafeMpi( MPI_Bcast( buf, count, MPI_COMPLEX, root, comm ) );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void Broadcast( dcomplex* buf, int count, int root, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Broadcast");
#endif
#ifdef AVOID_COMPLEX_MPI
    SafeMpi( MPI_Bcast( buf, 2*count, MPI_DOUBLE, root, comm ) );
#else
    SafeMpi( MPI_Bcast( buf, count, MPI_DOUBLE_COMPLEX, root, comm ) );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifdef HAVE_NONBLOCKING_COLLECTIVES
void IBroadcast( byte* buf, int count, int root, Comm comm, Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::IBroadcast");
#endif
    SafeMpi( 
        MPI_Ibcast( buf, count, MPI_UNSIGNED_CHAR, root, comm, &request ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void IBroadcast( int* buf, int count, int root, Comm comm, Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::IBroadcast");
#endif
    SafeMpi( MPI_Ibcast( buf, count, MPI_INT, root, comm, &request ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

void IBroadcast( float* buf, int count, int root, Comm comm, Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::IBroadcast");
#endif
    SafeMpi( MPI_Ibcast( buf, count, MPI_FLOAT, root, comm, &request ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

void IBroadcast( double* buf, int count, int root, Comm comm, Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::IBroadcast");
#endif
    SafeMpi( MPI_Ibcast( buf, count, MPI_DOUBLE, root, comm, &request ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

void IBroadcast
( scomplex* buf, int count, int root, Comm comm, Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::IBroadcast");
#endif
#ifdef AVOID_COMPLEX_MPI
    SafeMpi( MPI_Ibcast( buf, 2*count, MPI_FLOAT, root, comm, &request ) );
#else
    SafeMpi( MPI_Ibcast( buf, count, MPI_COMPLEX, root, comm, &request ) );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void IBroadcast
( dcomplex* buf, int count, int root, Comm comm, Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::IBroadcast");
#endif
#ifdef AVOID_COMPLEX_MPI
    SafeMpi( MPI_Ibcast( buf, 2*count, MPI_DOUBLE, root, comm, &request ) );
#else
    SafeMpi( 
        MPI_Ibcast( buf, count, MPI_DOUBLE_COMPLEX, root, comm, &request ) 
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // ifdef HAVE_NONBLOCKING_COLLECTIVES

void Gather
( const byte* sbuf, int sc,
        byte* rbuf, int rc, int root, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Gather");
#endif
    SafeMpi( 
        MPI_Gather
        ( const_cast<byte*>(sbuf), sc, MPI_UNSIGNED_CHAR,
          rbuf,                    rc, MPI_UNSIGNED_CHAR, root, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void Gather
( const int* sbuf, int sc,
        int* rbuf, int rc, int root, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Gather");
#endif
    SafeMpi( 
        MPI_Gather
        ( const_cast<int*>(sbuf), sc, MPI_INT,
          rbuf,                   rc, MPI_INT, root, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void Gather
( const float* sbuf, int sc,
        float* rbuf, int rc, int root, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Gather");
#endif
    SafeMpi( 
        MPI_Gather
        ( const_cast<float*>(sbuf), sc, MPI_FLOAT,
          rbuf,                     rc, MPI_FLOAT, root, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void Gather
( const double* sbuf, int sc,
        double* rbuf, int rc, int root, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Gather");
#endif
    SafeMpi( 
        MPI_Gather
        ( const_cast<double*>(sbuf), sc, MPI_DOUBLE,
          rbuf,                      rc, MPI_DOUBLE, root, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void Gather
( const scomplex* sbuf, int sc,
        scomplex* rbuf, int rc, int root, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Gather");
#endif
#ifdef AVOID_COMPLEX_MPI
    SafeMpi(
        MPI_Gather
        ( const_cast<scomplex*>(sbuf), 2*sc, MPI_FLOAT,
          rbuf,                        2*rc, MPI_FLOAT, root, comm )
    );
#else
    SafeMpi( 
        MPI_Gather
        ( const_cast<scomplex*>(sbuf), sc, MPI_COMPLEX,
          rbuf,                        rc, MPI_COMPLEX, root, comm ) 
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void Gather
( const dcomplex* sbuf, int sc,
        dcomplex* rbuf, int rc, int root, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Gather");
#endif
#ifdef AVOID_COMPLEX_MPI
    SafeMpi(
        MPI_Gather
        ( const_cast<dcomplex*>(sbuf), 2*sc, MPI_DOUBLE,
          rbuf,                        2*rc, MPI_DOUBLE, root, comm )
    );
#else
    SafeMpi( 
        MPI_Gather
        ( const_cast<dcomplex*>(sbuf), sc, MPI_DOUBLE_COMPLEX,
          rbuf,                        rc, MPI_DOUBLE_COMPLEX, root, comm )
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifdef HAVE_NONBLOCKING_COLLECTIVES
void IGather
( const byte* sbuf, int sc,
        byte* rbuf, int rc, int root, Comm comm, Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::IGather");
#endif
    SafeMpi( 
        MPI_Igather
        ( const_cast<byte*>(sbuf), sc, MPI_UNSIGNED_CHAR,
          rbuf,                    rc, MPI_UNSIGNED_CHAR, 
          root, comm, &request ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void IGather
( const int* sbuf, int sc,
        int* rbuf, int rc, int root, Comm comm, Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::IGather");
#endif
    SafeMpi( 
        MPI_Igather
        ( const_cast<int*>(sbuf), sc, MPI_INT,
          rbuf,                   rc, MPI_INT, 
          root, comm, &request ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void Igather
( const float* sbuf, int sc,
        float* rbuf, int rc, int root, Comm comm, Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::IGather");
#endif
    SafeMpi( 
        MPI_Igather
        ( const_cast<float*>(sbuf), sc, MPI_FLOAT,
          rbuf,                     rc, MPI_FLOAT, 
          root, comm, &request ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void IGather
( const double* sbuf, int sc,
        double* rbuf, int rc, int root, Comm comm, Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::IGather");
#endif
    SafeMpi( 
        MPI_Igather
        ( const_cast<double*>(sbuf), sc, MPI_DOUBLE,
          rbuf,                      rc, MPI_DOUBLE, 
          root, comm, &request ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void IGather
( const scomplex* sbuf, int sc,
        scomplex* rbuf, int rc, int root, Comm comm, Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::IGather");
#endif
#ifdef AVOID_COMPLEX_MPI
    SafeMpi(
        MPI_Igather
        ( const_cast<scomplex*>(sbuf), 2*sc, MPI_FLOAT,
          rbuf,                        2*rc, MPI_FLOAT, 
          root, comm, &request )
    );
#else
    SafeMpi( 
        MPI_Igather
        ( const_cast<scomplex*>(sbuf), sc, MPI_COMPLEX,
          rbuf,                        rc, MPI_COMPLEX, 
          root, comm, &request ) 
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void IGather
( const dcomplex* sbuf, int sc,
        dcomplex* rbuf, int rc, int root, Comm comm, Request& request )
{
#ifndef RELEASE
    PushCallStack("mpi::IGather");
#endif
#ifdef AVOID_COMPLEX_MPI
    SafeMpi(
        MPI_Igather
        ( const_cast<dcomplex*>(sbuf), 2*sc, MPI_DOUBLE,
          rbuf,                        2*rc, MPI_DOUBLE, 
          root, comm, &request )
    );
#else
    SafeMpi( 
        MPI_Igather
        ( const_cast<dcomplex*>(sbuf), sc, MPI_DOUBLE_COMPLEX,
          rbuf,                        rc, MPI_DOUBLE_COMPLEX, 
          root, comm, &request )
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // ifdef HAVE_NONBLOCKING_COLLECTIVES

void Gather
( const byte* sbuf, int sc,
        byte* rbuf, const int* rcs, const int* rds, int root, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Gather");
#endif
    SafeMpi( 
        MPI_Gatherv
        ( const_cast<byte*>(sbuf), 
          sc,       
          MPI_UNSIGNED_CHAR,
          rbuf,                    
          const_cast<int*>(rcs), 
          const_cast<int*>(rds), 
          MPI_UNSIGNED_CHAR, 
          root, 
          comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void Gather
( const int* sbuf, int sc,
        int* rbuf, const int* rcs, const int* rds, int root, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Gather");
#endif
    SafeMpi( 
        MPI_Gatherv
        ( const_cast<int*>(sbuf), 
          sc,       
          MPI_INT,
          rbuf, 
          const_cast<int*>(rcs), 
          const_cast<int*>(rds), 
          MPI_INT, 
          root, 
          comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void Gather
( const float* sbuf, int sc,
        float* rbuf, const int* rcs, const int* rds, int root, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Gather");
#endif
    SafeMpi( 
        MPI_Gatherv
        ( const_cast<float*>(sbuf), 
          sc,       
          MPI_FLOAT,
          rbuf, 
          const_cast<int*>(rcs), 
          const_cast<int*>(rds), 
          MPI_FLOAT, 
          root, 
          comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void Gather
( const double* sbuf, int sc,
        double* rbuf, const int* rcs, const int* rds, int root, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Gather");
#endif
    SafeMpi( 
        MPI_Gatherv
        ( const_cast<double*>(sbuf), 
          sc,       
          MPI_DOUBLE,
          rbuf, 
          const_cast<int*>(rcs), 
          const_cast<int*>(rds), 
          MPI_DOUBLE, 
          root, 
          comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void Gather
( const scomplex* sbuf, int sc,
        scomplex* rbuf, const int* rcs, const int* rds, int root, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Gather");
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
    SafeMpi(
        MPI_Gatherv
        ( const_cast<scomplex*>(sbuf), 2*sc,                         MPI_FLOAT,
          rbuf,                        &rcsDouble[0], &rdsDouble[0], MPI_FLOAT, 
          root, comm )
    );
#else
    SafeMpi( 
        MPI_Gatherv
        ( const_cast<scomplex*>(sbuf), 
          sc,       
          MPI_COMPLEX,
          rbuf,  
          const_cast<int*>(rcs), 
          const_cast<int*>(rds), 
          MPI_COMPLEX, 
          root, 
          comm ) 
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void Gather
( const dcomplex* sbuf, int sc,
        dcomplex* rbuf, const int* rcs, const int* rds, int root, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Gather");
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
    SafeMpi(
        MPI_Gatherv
        ( const_cast<dcomplex*>(sbuf), 2*sc,                         MPI_DOUBLE,
          rbuf,                        &rcsDouble[0], &rdsDouble[0], MPI_DOUBLE,
          root, comm )
    );
#else
    SafeMpi( 
        MPI_Gatherv
        ( const_cast<dcomplex*>(sbuf), 
          sc,       
          MPI_DOUBLE_COMPLEX,
          rbuf,  
          const_cast<int*>(rcs), 
          const_cast<int*>(rds), 
          MPI_DOUBLE_COMPLEX, 
          root, 
          comm )
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void AllGather
( const byte* sbuf, int sc,
        byte* rbuf, int rc, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllGather");
#endif
    SafeMpi( 
        MPI_Allgather
        ( const_cast<byte*>(sbuf), sc, MPI_UNSIGNED_CHAR, 
          rbuf,                    rc, MPI_UNSIGNED_CHAR, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void AllGather
( const int* sbuf, int sc,
        int* rbuf, int rc, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllGather");
#endif
#ifdef USE_BYTE_ALLGATHERS
    SafeMpi( 
        MPI_Allgather
        ( const_cast<int*>(sbuf), sizeof(int)*sc, MPI_UNSIGNED_CHAR, 
          rbuf,                   sizeof(int)*rc, MPI_UNSIGNED_CHAR, comm ) 
    );
#else
    SafeMpi( 
        MPI_Allgather
        ( const_cast<int*>(sbuf), sc, MPI_INT, 
          rbuf,                   rc, MPI_INT, comm ) 
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void AllGather
( const float* sbuf, int sc,
        float* rbuf, int rc, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllGather");
#endif
#ifdef USE_BYTE_ALLGATHERS
    SafeMpi( 
        MPI_Allgather
        ( const_cast<float*>(sbuf), sizeof(float)*sc, MPI_UNSIGNED_CHAR, 
          rbuf,                     sizeof(float)*rc, MPI_UNSIGNED_CHAR, comm ) 
    );
#else
    SafeMpi( 
        MPI_Allgather
        ( const_cast<float*>(sbuf), sc, MPI_FLOAT, 
          rbuf,                     rc, MPI_FLOAT, comm ) 
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void AllGather
( const double* sbuf, int sc,
        double* rbuf, int rc, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllGather");
#endif
#ifdef USE_BYTE_ALLGATHERS
    SafeMpi( 
        MPI_Allgather
        ( const_cast<double*>(sbuf), sizeof(double)*sc, MPI_UNSIGNED_CHAR, 
          rbuf,                      sizeof(double)*rc, MPI_UNSIGNED_CHAR, 
          comm ) 
    );
#else
    SafeMpi( 
        MPI_Allgather
        ( const_cast<double*>(sbuf), sc, MPI_DOUBLE,
          rbuf,                      rc, MPI_DOUBLE, comm ) 
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void AllGather
( const scomplex* sbuf, int sc,
        scomplex* rbuf, int rc, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllGather");
#endif
#ifdef USE_BYTE_ALLGATHERS
    SafeMpi( 
        MPI_Allgather
        ( const_cast<scomplex*>(sbuf), 2*sizeof(float)*sc, MPI_UNSIGNED_CHAR, 
          rbuf,                        2*sizeof(float)*rc, MPI_UNSIGNED_CHAR, 
          comm ) 
    );
#else
 #ifdef AVOID_COMPLEX_MPI
    SafeMpi(
        MPI_Allgather
        ( const_cast<scomplex*>(sbuf), 2*sc, MPI_FLOAT,
          rbuf,                        2*rc, MPI_FLOAT, comm )
    );
 #else
    SafeMpi( 
        MPI_Allgather
        ( const_cast<scomplex*>(sbuf), sc, MPI_COMPLEX,
          rbuf,                        rc, MPI_COMPLEX, comm ) 
    );
 #endif
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void AllGather
( const dcomplex* sbuf, int sc,
        dcomplex* rbuf, int rc, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllGather");
#endif
#ifdef USE_BYTE_ALLGATHERS
    SafeMpi( 
        MPI_Allgather
        ( const_cast<dcomplex*>(sbuf), 2*sizeof(double)*sc, MPI_UNSIGNED_CHAR, 
          rbuf,                        2*sizeof(double)*rc, MPI_UNSIGNED_CHAR, 
          comm ) 
    );
#else
 #ifdef AVOID_COMPLEX_MPI
    SafeMpi(
        MPI_Allgather
        ( const_cast<dcomplex*>(sbuf), 2*sc, MPI_DOUBLE,
          rbuf,                        2*rc, MPI_DOUBLE, comm )
    );
 #else
    SafeMpi( 
        MPI_Allgather
        ( const_cast<dcomplex*>(sbuf), sc, MPI_DOUBLE_COMPLEX,
          rbuf,                        rc, MPI_DOUBLE_COMPLEX, comm ) 
    );
 #endif
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void AllGather
( const byte* sbuf, int sc,
        byte* rbuf, const int* rcs, const int* rds, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllGather");
#endif
    SafeMpi( 
        MPI_Allgatherv
        ( const_cast<byte*>(sbuf), 
          sc,       
          MPI_UNSIGNED_CHAR, 
          rbuf, 
          const_cast<int*>(rcs), 
          const_cast<int*>(rds), 
          MPI_UNSIGNED_CHAR, 
          comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void AllGather
( const int* sbuf, int sc,
        int* rbuf, const int* rcs, const int* rds, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllGather");
#endif
#ifdef USE_BYTE_ALLGATHERS
    const int commSize = mpi::CommSize( comm );
    std::vector<int> byteRcs( commSize ), byteRds( commSize );
    for( int i=0; i<commSize; ++i )
    {
        byteRcs[i] = sizeof(int)*rcs[i];
        byteRds[i] = sizeof(int)*rds[i];
    }
    SafeMpi( 
        MPI_Allgatherv
        ( const_cast<int*>(sbuf), sizeof(int)*sc,   MPI_UNSIGNED_CHAR, 
          rbuf, &byteRcs[0], &byteRds[0], MPI_UNSIGNED_CHAR, comm ) 
    );
#else
    SafeMpi( 
        MPI_Allgatherv
        ( const_cast<int*>(sbuf), 
          sc, 
          MPI_INT, 
          rbuf,   
          const_cast<int*>(rcs), 
          const_cast<int*>(rds), 
          MPI_INT, 
          comm ) 
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void AllGather
( const float* sbuf, int sc,
        float* rbuf, const int* rcs, const int* rds, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllGather");
#endif
#ifdef USE_BYTE_ALLGATHERS
    const int commSize = mpi::CommSize( comm );
    std::vector<int> byteRcs( commSize ), byteRds( commSize );
    for( int i=0; i<commSize; ++i )
    {
        byteRcs[i] = sizeof(float)*rcs[i];
        byteRds[i] = sizeof(float)*rds[i];
    }
    SafeMpi( 
        MPI_Allgatherv
        ( const_cast<float*>(sbuf), sizeof(float)*sc, MPI_UNSIGNED_CHAR, 
          rbuf, &byteRcs[0], &byteRds[0], MPI_UNSIGNED_CHAR, comm ) 
    );
#else
    SafeMpi( 
        MPI_Allgatherv
        ( const_cast<float*>(sbuf), 
          sc, 
          MPI_FLOAT, 
          rbuf, 
          const_cast<int*>(rcs), 
          const_cast<int*>(rds), 
          MPI_FLOAT, 
          comm ) 
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void AllGather
( const double* sbuf, int sc,
        double* rbuf, const int* rcs, const int* rds, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllGather");
#endif
#ifdef USE_BYTE_ALLGATHERS
    const int commSize = mpi::CommSize( comm );
    std::vector<int> byteRcs( commSize ), byteRds( commSize );
    for( int i=0; i<commSize; ++i )
    {
        byteRcs[i] = sizeof(double)*rcs[i];
        byteRds[i] = sizeof(double)*rds[i];
    }
    SafeMpi( 
        MPI_Allgatherv
        ( const_cast<double*>(sbuf), sizeof(double)*sc, MPI_UNSIGNED_CHAR, 
          rbuf, &byteRcs[0], &byteRds[0], MPI_UNSIGNED_CHAR, comm )
    );
#else
    SafeMpi( 
        MPI_Allgatherv
        ( const_cast<double*>(sbuf), 
          sc, 
          MPI_DOUBLE,
          rbuf, 
          const_cast<int*>(rcs), 
          const_cast<int*>(rds),
          MPI_DOUBLE, 
          comm ) 
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void AllGather
( const scomplex* sbuf, int sc,
        scomplex* rbuf, const int* rcs, const int* rds, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllGather");
#endif
#ifdef USE_BYTE_ALLGATHERS
    const int commSize = mpi::CommSize( comm );
    std::vector<int> byteRcs( commSize ), byteRds( commSize );
    for( int i=0; i<commSize; ++i )
    {
        byteRcs[i] = 2*sizeof(float)*rcs[i];
        byteRds[i] = 2*sizeof(float)*rds[i];
    }
    SafeMpi( 
        MPI_Allgatherv
        ( const_cast<scomplex*>(sbuf), 2*sizeof(float)*sc, MPI_UNSIGNED_CHAR, 
          rbuf, &byteRcs[0], &byteRds[0], MPI_UNSIGNED_CHAR, comm )
    );
#else
 #ifdef AVOID_COMPLEX_MPI
    const int commSize = mpi::CommSize( comm );
    std::vector<int> realRcs( commSize ), realRds( commSize );
    for( int i=0; i<commSize; ++i )
    {
        realRcs[i] = 2*rcs[i];
        realRds[i] = 2*rds[i];
    }
    SafeMpi(
        MPI_Allgatherv
        ( const_cast<scomplex*>(sbuf), 2*sc, MPI_FLOAT,
          rbuf, &realRcs[0], &realRds[0], MPI_FLOAT, comm )
    );
 #else
    SafeMpi( 
        MPI_Allgatherv
        ( const_cast<scomplex*>(sbuf), 
          sc, 
          MPI_COMPLEX,
          rbuf, 
          const_cast<int*>(rcs), 
          const_cast<int*>(rds), 
          MPI_COMPLEX, 
          comm ) 
    );
 #endif
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void AllGather
( const dcomplex* sbuf, int sc,
        dcomplex* rbuf, const int* rcs, const int* rds, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllGather");
#endif
#ifdef USE_BYTE_ALLGATHERS
    const int commSize = mpi::CommSize( comm );
    std::vector<int> byteRcs( commSize ), byteRds( commSize );
    for( int i=0; i<commSize; ++i )
    {
        byteRcs[i] = 2*sizeof(double)*rcs[i];
        byteRds[i] = 2*sizeof(double)*rds[i];
    }
    SafeMpi( 
        MPI_Allgatherv
        ( const_cast<dcomplex*>(sbuf), 2*sizeof(double)*sc, MPI_UNSIGNED_CHAR, 
          rbuf, &byteRcs[0], &byteRds[0], MPI_UNSIGNED_CHAR, comm )
    );
#else
 #ifdef AVOID_COMPLEX_MPI
    const int commSize = mpi::CommSize( comm );
    std::vector<int> realRcs( commSize ), realRds( commSize );
    for( int i=0; i<commSize; ++i )
    {
        realRcs[i] = 2*rcs[i];
        realRds[i] = 2*rds[i];
    }
    SafeMpi(
        MPI_Allgatherv
        ( const_cast<dcomplex*>(sbuf), 2*sc, MPI_DOUBLE,
          rbuf, &realRcs[0], &realRds[0], MPI_DOUBLE, comm )
    );
 #else
    SafeMpi( 
        MPI_Allgatherv
        ( const_cast<dcomplex*>(sbuf), 
          sc, 
          MPI_DOUBLE_COMPLEX,
          rbuf, 
          const_cast<int*>(rcs), 
          const_cast<int*>(rds), 
          MPI_DOUBLE_COMPLEX, 
          comm ) 
    );
 #endif
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void Scatter
( const byte* sbuf, int sc,
        byte* rbuf, int rc, int root, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Scatter");
#endif
    SafeMpi( 
        MPI_Scatter
        ( const_cast<byte*>(sbuf), sc, MPI_UNSIGNED_CHAR,
          rbuf,                    rc, MPI_UNSIGNED_CHAR, root, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void Scatter
( const int* sbuf, int sc,
        int* rbuf, int rc, int root, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Scatter");
#endif
    SafeMpi( 
        MPI_Scatter
        ( const_cast<int*>(sbuf), sc, MPI_INT,
          rbuf,                   rc, MPI_INT, root, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void Scatter
( const float* sbuf, int sc,
        float* rbuf, int rc, int root, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Scatter");
#endif
    SafeMpi( 
        MPI_Scatter
        ( const_cast<float*>(sbuf), sc, MPI_FLOAT,
          rbuf,                     rc, MPI_FLOAT, root, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void Scatter
( const double* sbuf, int sc,
        double* rbuf, int rc, int root, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Scatter");
#endif
    SafeMpi( 
        MPI_Scatter
        ( const_cast<double*>(sbuf), sc, MPI_DOUBLE,
          rbuf,                      rc, MPI_DOUBLE, root, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void Scatter
( const scomplex* sbuf, int sc,
        scomplex* rbuf, int rc, int root, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Scatter");
#endif
#ifdef AVOID_COMPLEX_MPI
    SafeMpi(
        MPI_Scatter
        ( const_cast<scomplex*>(sbuf), 2*sc, MPI_FLOAT,
          rbuf,                        2*rc, MPI_FLOAT, root, comm )
    );
#else
    SafeMpi( 
        MPI_Scatter
        ( const_cast<scomplex*>(sbuf), sc, MPI_COMPLEX,
          rbuf,                        rc, MPI_COMPLEX, root, comm ) 
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void Scatter
( const dcomplex* sbuf, int sc,
        dcomplex* rbuf, int rc, int root, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Scatter");
#endif
#ifdef AVOID_COMPLEX_MPI
    SafeMpi(
        MPI_Scatter
        ( const_cast<dcomplex*>(sbuf), 2*sc, MPI_DOUBLE,
          rbuf,                        2*rc, MPI_DOUBLE, root, comm )
    );
#else
    SafeMpi( 
        MPI_Scatter
        ( const_cast<dcomplex*>(sbuf), sc, MPI_DOUBLE_COMPLEX,
          rbuf,                        rc, MPI_DOUBLE_COMPLEX, root, comm ) 
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void Scatter( byte* buf, int sc, int rc, int root, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Scatter");
#endif
    const int commRank = CommRank( comm );
    if( commRank == root )
    {
        SafeMpi(
            MPI_Scatter
            ( buf,          sc, MPI_UNSIGNED_CHAR, 
              MPI_IN_PLACE, rc, MPI_UNSIGNED_CHAR, root, comm )
        );
    }
    else
    {
        SafeMpi(
            MPI_Scatter
            ( 0,   sc, MPI_UNSIGNED_CHAR, 
              buf, rc, MPI_UNSIGNED_CHAR, root, comm )
        );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void Scatter( int* buf, int sc, int rc, int root, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Scatter");
#endif
    const int commRank = CommRank( comm );
    if( commRank == root )
    {
        SafeMpi(
            MPI_Scatter
            ( buf,          sc, MPI_INT, 
              MPI_IN_PLACE, rc, MPI_INT, root, comm )
        );
    }
    else
    {
        SafeMpi(
            MPI_Scatter
            ( 0,   sc, MPI_INT, 
              buf, rc, MPI_INT, root, comm )
        );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void Scatter( float* buf, int sc, int rc, int root, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Scatter");
#endif
    const int commRank = CommRank( comm );
    if( commRank == root )
    {
        SafeMpi(
            MPI_Scatter
            ( buf,          sc, MPI_FLOAT, 
              MPI_IN_PLACE, rc, MPI_FLOAT, root, comm )
        );
    }
    else
    {
        SafeMpi(
            MPI_Scatter
            ( 0,   sc, MPI_FLOAT, 
              buf, rc, MPI_FLOAT, root, comm )
        );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void Scatter( double* buf, int sc, int rc, int root, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Scatter");
#endif
    const int commRank = CommRank( comm );
    if( commRank == root )
    {
        SafeMpi(
            MPI_Scatter
            ( buf,          sc, MPI_DOUBLE, 
              MPI_IN_PLACE, rc, MPI_DOUBLE, root, comm )
        );
    }
    else
    {
        SafeMpi(
            MPI_Scatter
            ( 0,   sc, MPI_DOUBLE, 
              buf, rc, MPI_DOUBLE, root, comm )
        );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void Scatter( scomplex* buf, int sc, int rc, int root, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Scatter");
#endif
    const int commRank = CommRank( comm );
    if( commRank == root )
    {
#ifdef AVOID_COMPLEX_MPI
        SafeMpi(
            MPI_Scatter
            ( buf,          2*sc, MPI_FLOAT, 
              MPI_IN_PLACE, 2*rc, MPI_FLOAT, root, comm )
        );
#else
        SafeMpi(
            MPI_Scatter
            ( buf,          sc, MPI_COMPLEX, 
              MPI_IN_PLACE, rc, MPI_COMPLEX, root, comm )
        );
#endif
    }
    else
    {
#ifdef AVOID_COMPLEX_MPI
        SafeMpi(
            MPI_Scatter
            ( 0,   2*sc, MPI_FLOAT, 
              buf, 2*rc, MPI_FLOAT, root, comm )
        );
#else
        SafeMpi(
            MPI_Scatter
            ( 0,   sc, MPI_COMPLEX, 
              buf, rc, MPI_COMPLEX, root, comm )
        );
#endif
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void Scatter( dcomplex* buf, int sc, int rc, int root, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Scatter");
#endif
    const int commRank = CommRank( comm );
    if( commRank == root )
    {
#ifdef AVOID_COMPLEX_MPI
        SafeMpi(
            MPI_Scatter
            ( buf,          2*sc, MPI_DOUBLE, 
              MPI_IN_PLACE, 2*rc, MPI_DOUBLE, root, comm )
        );
#else
        SafeMpi(
            MPI_Scatter
            ( buf,          sc, MPI_DOUBLE_COMPLEX, 
              MPI_IN_PLACE, rc, MPI_DOUBLE_COMPLEX, root, comm )
        );
#endif
    }
    else
    {
#ifdef AVOID_COMPLEX_MPI
        SafeMpi(
            MPI_Scatter
            ( 0,   2*sc, MPI_DOUBLE, 
              buf, 2*rc, MPI_DOUBLE, root, comm )
        );
#else
        SafeMpi(
            MPI_Scatter
            ( 0,   sc, MPI_DOUBLE_COMPLEX, 
              buf, rc, MPI_DOUBLE_COMPLEX, root, comm )
        );
#endif
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void AllToAll
( const byte* sbuf, int sc,
        byte* rbuf, int rc, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllToAll");
#endif
    SafeMpi( 
        MPI_Alltoall
        ( const_cast<byte*>(sbuf), sc, MPI_UNSIGNED_CHAR,
          rbuf,                    rc, MPI_UNSIGNED_CHAR, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void AllToAll
( const int* sbuf, int sc,
        int* rbuf, int rc, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllToAll");
#endif
    SafeMpi( 
        MPI_Alltoall
        ( const_cast<int*>(sbuf), sc, MPI_INT,
          rbuf,                   rc, MPI_INT, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void AllToAll
( const float* sbuf, int sc,
        float* rbuf, int rc, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllToAll");
#endif
    SafeMpi( 
        MPI_Alltoall
        ( const_cast<float*>(sbuf), sc, MPI_FLOAT,
          rbuf,                     rc, MPI_FLOAT, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void AllToAll
( const double* sbuf, int sc,
        double* rbuf, int rc, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllToAll");
#endif
    SafeMpi( 
        MPI_Alltoall
        ( const_cast<double*>(sbuf), sc, MPI_DOUBLE,
          rbuf,                      rc, MPI_DOUBLE, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void AllToAll
( const scomplex* sbuf, int sc,
        scomplex* rbuf, int rc, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllToAll");
#endif
#ifdef AVOID_COMPLEX_MPI
    SafeMpi(
        MPI_Alltoall
        ( const_cast<scomplex*>(sbuf), 2*sc, MPI_FLOAT,
          rbuf,                        2*rc, MPI_FLOAT, comm )
    );
#else
    SafeMpi( 
        MPI_Alltoall
        ( const_cast<scomplex*>(sbuf), sc, MPI_COMPLEX,
          rbuf,                        rc, MPI_COMPLEX, comm ) 
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void AllToAll
( const dcomplex* sbuf, int sc,
        dcomplex* rbuf, int rc, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllToAll");
#endif
#ifdef AVOID_COMPLEX_MPI
    SafeMpi(
        MPI_Alltoall
        ( const_cast<dcomplex*>(sbuf), 2*sc, MPI_DOUBLE,
          rbuf,                        2*rc, MPI_DOUBLE, comm )
    );
#else
    SafeMpi( 
        MPI_Alltoall
        ( const_cast<dcomplex*>(sbuf), sc, MPI_DOUBLE_COMPLEX,
          rbuf,                        rc, MPI_DOUBLE_COMPLEX, comm ) 
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void AllToAll
( const byte* sbuf, const int* scs, const int* sds, 
        byte* rbuf, const int* rcs, const int* rds, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllToAll");
#endif
    SafeMpi( 
        MPI_Alltoallv
        ( const_cast<byte*>(sbuf), 
          const_cast<int*>(scs), 
          const_cast<int*>(sds), 
          MPI_UNSIGNED_CHAR,
          rbuf, 
          const_cast<int*>(rcs), 
          const_cast<int*>(rds), 
          MPI_UNSIGNED_CHAR, 
          comm ) 
    ); 
#ifndef RELEASE
    PopCallStack();
#endif
}

void AllToAll
( const int* sbuf, const int* scs, const int* sds, 
        int* rbuf, const int* rcs, const int* rds, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllToAll");
#endif
    SafeMpi( 
        MPI_Alltoallv
        ( const_cast<int*>(sbuf), 
          const_cast<int*>(scs), 
          const_cast<int*>(sds), 
          MPI_INT,
          rbuf, 
          const_cast<int*>(rcs), 
          const_cast<int*>(rds), 
          MPI_INT, 
          comm ) 
    ); 
#ifndef RELEASE
    PopCallStack();
#endif
}

void AllToAll
( const float* sbuf, const int* scs, const int* sds,
        float* rbuf, const int* rcs, const int* rds, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllToAll");
#endif
    SafeMpi( 
        MPI_Alltoallv
        ( const_cast<float*>(sbuf), 
          const_cast<int*>(scs), 
          const_cast<int*>(sds), 
          MPI_FLOAT,
          rbuf, 
          const_cast<int*>(rcs), 
          const_cast<int*>(rds), 
          MPI_FLOAT, 
          comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void AllToAll
( const double* sbuf, const int* scs, const int* sds,
        double* rbuf, const int* rcs, const int* rds, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllToAll");
#endif
    SafeMpi( 
        MPI_Alltoallv
        ( const_cast<double*>(sbuf), 
          const_cast<int*>(scs), 
          const_cast<int*>(sds), 
          MPI_DOUBLE,
          rbuf, 
          const_cast<int*>(rcs), 
          const_cast<int*>(rds), 
          MPI_DOUBLE, 
          comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void AllToAll
( const scomplex* sbuf, const int* scs, const int* sds,
        scomplex* rbuf, const int* rcs, const int* rds, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllToAll");
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
    SafeMpi(
        MPI_Alltoallv
        ( const_cast<scomplex*>(sbuf),
                &scsDoubled[0], &sdsDoubled[0], MPI_FLOAT,
          rbuf, &rcsDoubled[0], &rdsDoubled[0], MPI_FLOAT, comm )
    );
#else
    SafeMpi( 
        MPI_Alltoallv
        ( const_cast<scomplex*>(sbuf), 
          const_cast<int*>(scs), 
          const_cast<int*>(sds), 
          MPI_COMPLEX,
          rbuf, 
          const_cast<int*>(rcs), 
          const_cast<int*>(rds), 
          MPI_COMPLEX, 
          comm )
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void AllToAll
( const dcomplex* sbuf, const int* scs, const int* sds,
        dcomplex* rbuf, const int* rcs, const int* rds, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllToAll");
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
    SafeMpi(
        MPI_Alltoallv
        ( const_cast<dcomplex*>(sbuf),
                &scsDoubled[0], &sdsDoubled[0], MPI_DOUBLE,
          rbuf, &rcsDoubled[0], &rdsDoubled[0], MPI_DOUBLE, comm )
    );
#else
    SafeMpi( 
        MPI_Alltoallv
        ( const_cast<dcomplex*>(sbuf), 
          const_cast<int*>(scs), 
          const_cast<int*>(sds), 
          MPI_DOUBLE_COMPLEX,
          rbuf, 
          const_cast<int*>(rcs), 
          const_cast<int*>(rds), 
          MPI_DOUBLE_COMPLEX, 
          comm ) 
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void Reduce
( const byte* sbuf, byte* rbuf, int count, Op op, int root, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Reduce");
#endif
    if( count != 0 )
    {
        SafeMpi( 
            MPI_Reduce
            ( const_cast<byte*>(sbuf), rbuf, count, MPI_UNSIGNED_CHAR, op, root,
              comm ) 
        );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void Reduce( byte* buf, int count, Op op, int root, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Reduce");
#endif
    if( count != 0 )
    {
        const int commRank = CommRank( comm );
        if( commRank == root )
            SafeMpi( 
                MPI_Reduce
                ( MPI_IN_PLACE, buf, count, MPI_UNSIGNED_CHAR, op, root, comm )
            );
        else
            SafeMpi(
                MPI_Reduce
                ( buf, 0, count, MPI_UNSIGNED_CHAR, op, root, comm )
            );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void Reduce
( const int* sbuf, int* rbuf, int count, Op op, int root, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Reduce");
#endif
    if( count != 0 )
    {
        SafeMpi( 
            MPI_Reduce
            ( const_cast<int*>(sbuf), rbuf, count, MPI_INT, op, root, comm ) 
        );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void Reduce( int* buf, int count, Op op, int root, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Reduce");
#endif
    if( count != 0 )
    {
        const int commRank = CommRank( comm );
        if( commRank == root )
            SafeMpi( 
                MPI_Reduce
                ( MPI_IN_PLACE, buf, count, MPI_INT, op, root, comm )
            );
        else
            SafeMpi(
                MPI_Reduce
                ( buf, 0, count, MPI_INT, op, root, comm )
            );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void Reduce
( const float* sbuf, float* rbuf, int count, Op op, int root, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Reduce");
#endif
    if( count != 0 )
    {
        SafeMpi( 
            MPI_Reduce
            ( const_cast<float*>(sbuf), rbuf, count, MPI_FLOAT, op, root, comm ) 
        );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void Reduce( float* buf, int count, Op op, int root, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Reduce");
#endif
    if( count != 0 )
    {
        const int commRank = CommRank( comm );
        if( commRank == root )
            SafeMpi( 
                MPI_Reduce
                ( MPI_IN_PLACE, buf, count, MPI_FLOAT, op, root, comm )
            );
        else
            SafeMpi(
                MPI_Reduce
                ( buf, 0, count, MPI_FLOAT, op, root, comm )
            );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void Reduce
( const double* sbuf, double* rbuf, int count, Op op, int root, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Reduce");
#endif
    if( count != 0 )
    {
        SafeMpi( 
            MPI_Reduce
            ( const_cast<double*>(sbuf), rbuf, count, MPI_DOUBLE, op, root, 
              comm ) 
        );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void Reduce( double* buf, int count, Op op, int root, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Reduce");
#endif
    if( count != 0 )
    {
        const int commRank = CommRank( comm );
        if( commRank == root )
            SafeMpi( 
                MPI_Reduce
                ( MPI_IN_PLACE, buf, count, MPI_DOUBLE, op, root, comm )
            );
        else
            SafeMpi(
                MPI_Reduce
                ( buf, 0, count, MPI_DOUBLE, op, root, comm )
            );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void Reduce
( const scomplex* sbuf, scomplex* rbuf, int count, Op op, int root, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Reduce");
#endif
    if( count != 0 )
    {
#ifdef AVOID_COMPLEX_MPI
        if( op == SUM )
        {
            SafeMpi(
                MPI_Reduce
                ( const_cast<scomplex*>(sbuf),
                  rbuf, 2*count, MPI_FLOAT, op, root, comm )
            );
        }
        else
        {
            SafeMpi(
                MPI_Reduce
                ( const_cast<scomplex*>(sbuf),
                  rbuf, count, MPI_COMPLEX, op, root, comm )
            );
        }
#else
        SafeMpi( 
            MPI_Reduce
            ( const_cast<scomplex*>(sbuf), 
              rbuf, count, MPI_COMPLEX, op, root, comm ) 
        );
#endif
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void Reduce( scomplex* buf, int count, Op op, int root, Comm comm )
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
            if( commRank == root )
                SafeMpi(
                    MPI_Reduce
                    ( MPI_IN_PLACE, buf, 2*count, MPI_FLOAT, op, root, comm )
                );
            else
                SafeMpi(
                    MPI_Reduce
                    ( buf, 0, 2*count, MPI_FLOAT, op, root, comm )
                );
        }
        else
        {
            if( commRank == root )
                SafeMpi(
                    MPI_Reduce
                    ( MPI_IN_PLACE, buf, count, MPI_COMPLEX, op, root, comm )
                );
            else
                SafeMpi(
                    MPI_Reduce
                    ( buf, 0, count, MPI_COMPLEX, op, root, comm )
                );
        }
#else
        if( commRank == root )
            SafeMpi( 
                MPI_Reduce
                ( MPI_IN_PLACE, buf, count, MPI_COMPLEX, op, root, comm )
            );
        else
            SafeMpi(
                MPI_Reduce
                ( buf, 0, count, MPI_COMPLEX, op, root, comm )
            );
#endif
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void Reduce
( const dcomplex* sbuf, dcomplex* rbuf, int count, Op op, int root, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::Reduce");
#endif
    if( count != 0 )
    {
#ifdef AVOID_COMPLEX_MPI
        if( op == SUM )
        {
            SafeMpi(
                MPI_Reduce
                ( const_cast<dcomplex*>(sbuf),
                  rbuf, 2*count, MPI_DOUBLE, op, root, comm )
            );
        }
        else
        {
            SafeMpi(
                MPI_Reduce
                ( const_cast<dcomplex*>(sbuf),
                  rbuf, count, MPI_DOUBLE_COMPLEX, op, root, comm )
            );
        }
#else
        SafeMpi( 
            MPI_Reduce
            ( const_cast<dcomplex*>(sbuf), 
              rbuf, count, MPI_DOUBLE_COMPLEX, op, root, comm ) 
        );
#endif
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void Reduce( dcomplex* buf, int count, Op op, int root, Comm comm )
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
            if( commRank == root )
                SafeMpi(
                    MPI_Reduce
                    ( MPI_IN_PLACE, buf, 2*count, MPI_DOUBLE, op, root, comm )
                );
            else
                SafeMpi(
                    MPI_Reduce
                    ( buf, 0, 2*count, MPI_DOUBLE, op, root, comm )
                );
        }
        else
        {
            if( commRank == root )
                SafeMpi(
                    MPI_Reduce
                    ( MPI_IN_PLACE, buf, count, MPI_DOUBLE_COMPLEX, 
                      op, root, comm )
                );
            else
                SafeMpi(
                    MPI_Reduce
                    ( buf, 0, count, MPI_DOUBLE_COMPLEX, op, root, comm )
                );
        }
#else
        if( commRank == root )
            SafeMpi( 
                MPI_Reduce
                ( MPI_IN_PLACE, buf, count, MPI_DOUBLE_COMPLEX, op, root, comm )
            );
        else
            SafeMpi(
                MPI_Reduce
                ( buf, 0, count, MPI_DOUBLE_COMPLEX, op, root, comm )
            );
#endif
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void AllReduce( const byte* sbuf, byte* rbuf, int count, Op op, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllReduce");
#endif
    if( count != 0 )
    {
        SafeMpi( 
            MPI_Allreduce
            ( const_cast<byte*>(sbuf), rbuf, count, MPI_UNSIGNED_CHAR, op, 
              comm ) 
        );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void AllReduce( byte* buf, int count, Op op, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllReduce");
#endif
    if( count != 0 )
    {
        SafeMpi( 
            MPI_Allreduce
            ( MPI_IN_PLACE, buf, count, MPI_UNSIGNED_CHAR, op, comm )
        );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void AllReduce( const int* sbuf, int* rbuf, int count, Op op, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllReduce");
#endif
    if( count != 0 )
    {
        SafeMpi( 
            MPI_Allreduce
            ( const_cast<int*>(sbuf), rbuf, count, MPI_INT, op, comm ) 
        );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void AllReduce( int* buf, int count, Op op, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllReduce");
#endif
    if( count != 0 )
    {
        SafeMpi( 
            MPI_Allreduce( MPI_IN_PLACE, buf, count, MPI_INT, op, comm )
        );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void AllReduce( const float* sbuf, float* rbuf, int count, Op op, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllReduce");
#endif
    if( count != 0 )
    {
        SafeMpi( 
            MPI_Allreduce
            ( const_cast<float*>(sbuf), rbuf, count, MPI_FLOAT, op, comm ) 
        );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void AllReduce( float* buf, int count, Op op, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllReduce");
#endif
    if( count != 0 )
    {
        SafeMpi( 
            MPI_Allreduce( MPI_IN_PLACE, buf, count, MPI_FLOAT, op, comm )
        );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void AllReduce( const double* sbuf, double* rbuf, int count, Op op, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllReduce");
#endif
    if( count != 0 )
    {
        SafeMpi( 
            MPI_Allreduce
            ( const_cast<double*>(sbuf), rbuf, count, MPI_DOUBLE, op, comm ) 
        );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void AllReduce( double* buf, int count, Op op, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllReduce");
#endif
    if( count != 0 )
    {
        SafeMpi( 
            MPI_Allreduce( MPI_IN_PLACE, buf, count, MPI_DOUBLE, op, comm )
        );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void AllReduce
( const scomplex* sbuf, scomplex* rbuf, int count, Op op, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllReduce");
#endif
    if( count != 0 )
    {
#ifdef AVOID_COMPLEX_MPI
        if( op == SUM )
        {
            SafeMpi(
                MPI_Allreduce
                ( const_cast<scomplex*>(sbuf),
                  rbuf, 2*count, MPI_FLOAT, op, comm )
            );
        }
        else
        {
            SafeMpi(
                MPI_Allreduce
                ( const_cast<scomplex*>(sbuf),
                  rbuf, count, MPI_COMPLEX, op, comm )
            );
        }
#else
        SafeMpi( 
            MPI_Allreduce
            ( const_cast<scomplex*>(sbuf), 
              rbuf, count, MPI_COMPLEX, op, comm ) 
        );
#endif
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void AllReduce( scomplex* buf, int count, Op op, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllReduce");
#endif
    if( count != 0 )
    {
#ifdef AVOID_COMPLEX_MPI
        if( op == SUM )
        {
            SafeMpi(
                MPI_Allreduce( MPI_IN_PLACE, buf, 2*count, MPI_FLOAT, op, comm )
            );
        }
        else
        {
            SafeMpi(
                MPI_Allreduce( MPI_IN_PLACE, buf, count, MPI_COMPLEX, op, comm )
            );
        }
#else
        SafeMpi( 
            MPI_Allreduce( MPI_IN_PLACE, buf, count, MPI_COMPLEX, op, comm )
        );
#endif
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void AllReduce
( const dcomplex* sbuf, dcomplex* rbuf, int count, Op op, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllReduce");
#endif
    if( count != 0 )
    {
#ifdef AVOID_COMPLEX_MPI
        if( op == SUM )
        {
            SafeMpi(
                MPI_Allreduce
                ( const_cast<dcomplex*>(sbuf),
                  rbuf, 2*count, MPI_DOUBLE, op, comm )
            );
        }
        else
        {
            SafeMpi(
                MPI_Allreduce
                ( const_cast<dcomplex*>(sbuf),
                  rbuf, count, MPI_DOUBLE_COMPLEX, op, comm )
            );
        }
#else
        SafeMpi( 
            MPI_Allreduce
            ( const_cast<dcomplex*>(sbuf), 
              rbuf, count, MPI_DOUBLE_COMPLEX, op, comm )
        );
#endif
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void AllReduce( dcomplex* buf, int count, Op op, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::AllReduce");
#endif
    if( count != 0 )
    {
#ifdef AVOID_COMPLEX_MPI
        if( op == SUM )
        {
            SafeMpi(
                MPI_Allreduce
                ( MPI_IN_PLACE, buf, 2*count, MPI_DOUBLE, op, comm )
            );
        }
        else
        {
            SafeMpi(
                MPI_Allreduce
                ( MPI_IN_PLACE, buf, count, MPI_DOUBLE_COMPLEX, op, comm )
            );
        }
#else
        SafeMpi( 
            MPI_Allreduce
            ( MPI_IN_PLACE, buf, count, MPI_DOUBLE_COMPLEX, op, comm )
        );
#endif
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

void ReduceScatter( byte* sbuf, byte* rbuf, int rc, Op op, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::ReduceScatter");
#endif
#ifdef HAVE_REDUCE_SCATTER_BLOCK
    SafeMpi( 
        MPI_Reduce_scatter_block( sbuf, rbuf, rc, MPI_UNSIGNED_CHAR, op, comm )
    );
#elif REDUCE_SCATTER_BLOCK_VIA_ALLREDUCE
    const int commSize = CommSize( comm );
    const int commRank = CommRank( comm );
    AllReduce( sbuf, rc*commSize, op, comm );
    std::memcpy( rbuf, &sbuf[commRank*rc], rc );
#else
    const int commSize = CommSize( comm );
    Reduce( sbuf, rc*commSize, op, 0, comm );
    Scatter( sbuf, rc, rbuf, rc, 0, comm );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void ReduceScatter( int* sbuf, int* rbuf, int rc, Op op, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::ReduceScatter");
#endif
#ifdef HAVE_REDUCE_SCATTER_BLOCK
    SafeMpi( 
        MPI_Reduce_scatter_block( sbuf, rbuf, rc, MPI_INT, op, comm )
    );
#elif REDUCE_SCATTER_BLOCK_VIA_ALLREDUCE
    const int commSize = CommSize( comm );
    const int commRank = CommRank( comm );
    AllReduce( sbuf, rc*commSize, op, comm );
    std::memcpy( rbuf, &sbuf[commRank*rc], rc*sizeof(int) );
#else
    const int commSize = CommSize( comm );
    Reduce( sbuf, rc*commSize, op, 0, comm );
    Scatter( sbuf, rc, rbuf, rc, 0, comm );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void ReduceScatter( float* sbuf, float* rbuf, int rc, Op op, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::ReduceScatter");
#endif
#ifdef HAVE_REDUCE_SCATTER_BLOCK
    SafeMpi( 
        MPI_Reduce_scatter_block( sbuf, rbuf, rc, MPI_FLOAT, op, comm )
    );
#elif REDUCE_SCATTER_BLOCK_VIA_ALLREDUCE
    const int commSize = CommSize( comm );
    const int commRank = CommRank( comm );
    AllReduce( sbuf, rc*commSize, op, comm );
    std::memcpy( rbuf, &sbuf[commRank*rc], rc*sizeof(float) );
#else
    const int commSize = CommSize( comm );
    Reduce( sbuf, rc*commSize, op, 0, comm );
    Scatter( sbuf, rc, rbuf, rc, 0, comm );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void ReduceScatter( double* sbuf, double* rbuf, int rc, Op op, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::ReduceScatter");
#endif
#ifdef HAVE_REDUCE_SCATTER_BLOCK
    SafeMpi( 
        MPI_Reduce_scatter_block( sbuf, rbuf, rc, MPI_DOUBLE, op, comm )
    );
#elif REDUCE_SCATTER_BLOCK_VIA_ALLREDUCE
    const int commSize = CommSize( comm );
    const int commRank = CommRank( comm );
    AllReduce( sbuf, rc*commSize, op, comm );
    std::memcpy( rbuf, &sbuf[commRank*rc], rc*sizeof(double) );
#else
    const int commSize = CommSize( comm );
    Reduce( sbuf, rc*commSize, op, 0, comm );
    Scatter( sbuf, rc, rbuf, rc, 0, comm );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void ReduceScatter( scomplex* sbuf, scomplex* rbuf, int rc, Op op, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::ReduceScatter");
#endif
#ifdef HAVE_REDUCE_SCATTER_BLOCK
# ifdef AVOID_COMPLEX_MPI
    SafeMpi(
        MPI_Reduce_scatter_block
        ( sbuf, rbuf, 2*rc, MPI_FLOAT, op, comm )
    );
# else
    SafeMpi( 
        MPI_Reduce_scatter_block
        ( sbuf, rbuf, rc, MPI_COMPLEX, op, comm )
    );
# endif
#elif REDUCE_SCATTER_BLOCK_VIA_ALLREDUCE
    const int commSize = CommSize( comm );
    const int commRank = CommRank( comm );
    AllReduce( sbuf, rc*commSize, op, comm );
    std::memcpy( rbuf, &sbuf[commRank*rc], rc*sizeof(scomplex) );
#else
    const int commSize = CommSize( comm );
    Reduce( sbuf, rc*commSize, op, 0, comm );
    Scatter( sbuf, rc, rbuf, rc, 0, comm );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void ReduceScatter( dcomplex* sbuf, dcomplex* rbuf, int rc, Op op, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::ReduceScatter");
#endif
#ifdef HAVE_REDUCE_SCATTER_BLOCK
# ifdef AVOID_COMPLEX_MPI
    SafeMpi(
        MPI_Reduce_scatter_block
        ( sbuf, rbuf, 2*rc, MPI_DOUBLE, op, comm )
    );
# else
    SafeMpi( 
        MPI_Reduce_scatter_block
        ( sbuf, rbuf, rc, MPI_DOUBLE_COMPLEX, op, comm )
    );
# endif
#elif REDUCE_SCATTER_BLOCK_VIA_ALLREDUCE
    const int commSize = CommSize( comm );
    const int commRank = CommRank( comm );
    AllReduce( sbuf, rc*commSize, op, comm );
    std::memcpy( rbuf, &sbuf[commRank*rc], rc*sizeof(dcomplex) );
#else
    const int commSize = CommSize( comm );
    Reduce( sbuf, rc*commSize, op, 0, comm );
    Scatter( sbuf, rc, rbuf, rc, 0, comm );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void ReduceScatter( byte* buf, int rc, Op op, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::ReduceScatter");
#endif
#ifdef HAVE_REDUCE_SCATTER_BLOCK
    SafeMpi( 
        MPI_Reduce_scatter_block
        ( MPI_IN_PLACE, buf, rc, MPI_UNSIGNED_CHAR, op, comm )
    );
#elif REDUCE_SCATTER_BLOCK_VIA_ALLREDUCE
    const int commSize = CommSize( comm );
    const int commRank = CommRank( comm );
    AllReduce( buf, rc*commSize, op, comm );
    if( commRank != 0 )
        std::memcpy( buf, &buf[commRank*rc], rc );
#else
    const int commSize = CommSize( comm );
    Reduce( buf, rc*commSize, op, 0, comm );
    Scatter( buf, rc, rc, 0, comm );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void ReduceScatter( int* buf, int rc, Op op, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::ReduceScatter");
#endif
#ifdef HAVE_REDUCE_SCATTER_BLOCK
    SafeMpi( 
        MPI_Reduce_scatter_block( MPI_IN_PLACE, buf, rc, MPI_INT, op, comm )
    );
#elif REDUCE_SCATTER_BLOCK_VIA_ALLREDUCE
    const int commSize = CommSize( comm );
    const int commRank = CommRank( comm );
    AllReduce( buf, rc*commSize, op, comm );
    if( commRank != 0 )
        std::memcpy( buf, &buf[commRank*rc], rc*sizeof(int) );
#else
    const int commSize = CommSize( comm );
    Reduce( buf, rc*commSize, op, 0, comm );
    Scatter( buf, rc, rc, 0, comm );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void ReduceScatter( float* buf, int rc, Op op, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::ReduceScatter");
#endif
#ifdef HAVE_REDUCE_SCATTER_BLOCK
    SafeMpi( 
        MPI_Reduce_scatter_block( MPI_IN_PLACE, buf, rc, MPI_FLOAT, op, comm )
    );
#elif REDUCE_SCATTER_BLOCK_VIA_ALLREDUCE
    const int commSize = CommSize( comm );
    const int commRank = CommRank( comm );
    AllReduce( buf, rc*commSize, op, comm );
    if( commRank != 0 )
        std::memcpy( buf, &buf[commRank*rc], rc*sizeof(float) );
#else
    const int commSize = CommSize( comm );
    Reduce( buf, rc*commSize, op, 0, comm );
    Scatter( buf, rc, rc, 0, comm );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void ReduceScatter( double* buf, int rc, Op op, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::ReduceScatter");
#endif
#ifdef HAVE_REDUCE_SCATTER_BLOCK
    SafeMpi( 
        MPI_Reduce_scatter_block( MPI_IN_PLACE, buf, rc, MPI_DOUBLE, op, comm )
    );
#elif REDUCE_SCATTER_BLOCK_VIA_ALLREDUCE
    const int commSize = CommSize( comm );
    const int commRank = CommRank( comm );
    AllReduce( buf, rc*commSize, op, comm );
    if( commRank != 0 )
        std::memcpy( buf, &buf[commRank*rc], rc*sizeof(double) );
#else
    const int commSize = CommSize( comm );
    Reduce( buf, rc*commSize, op, 0, comm );
    Scatter( buf, rc, rc, 0, comm );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void ReduceScatter( scomplex* buf, int rc, Op op, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::ReduceScatter");
#endif
#ifdef HAVE_REDUCE_SCATTER_BLOCK
# ifdef AVOID_COMPLEX_MPI
    SafeMpi(
        MPI_Reduce_scatter_block
        ( MPI_IN_PLACE, buf, 2*rc, MPI_FLOAT, op, comm )
    );
# else
    SafeMpi( 
        MPI_Reduce_scatter_block
        ( MPI_IN_PLACE, buf, rc, MPI_COMPLEX, op, comm )
    );
# endif
#elif REDUCE_SCATTER_BLOCK_VIA_ALLREDUCE
    const int commSize = CommSize( comm );
    const int commRank = CommRank( comm );
    AllReduce( buf, rc*commSize, op, comm );
    if( commRank != 0 )
        std::memcpy( buf, &buf[commRank*rc], rc*sizeof(scomplex) );
#else
    const int commSize = CommSize( comm );
    Reduce( buf, rc*commSize, op, 0, comm );
    Scatter( buf, rc, rc, 0, comm );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void ReduceScatter( dcomplex* buf, int rc, Op op, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::ReduceScatter");
#endif
#ifdef HAVE_REDUCE_SCATTER_BLOCK
# ifdef AVOID_COMPLEX_MPI
    SafeMpi(
        MPI_Reduce_scatter_block
        ( MPI_IN_PLACE, buf, 2*rc, MPI_DOUBLE, op, comm )
    );
# else
    SafeMpi( 
        MPI_Reduce_scatter_block
        ( MPI_IN_PLACE, buf, rc, MPI_DOUBLE_COMPLEX, op, comm )
    );
# endif
#elif REDUCE_SCATTER_BLOCK_VIA_ALLREDUCE
    const int commSize = CommSize( comm );
    const int commRank = CommRank( comm );
    AllReduce( buf, rc*commSize, op, comm );
    if( commRank != 0 )
        std::memcpy( buf, &buf[commRank*rc], rc*sizeof(dcomplex) );
#else
    const int commSize = CommSize( comm );
    Reduce( buf, rc*commSize, op, 0, comm );
    Scatter( buf, rc, rc, 0, comm );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void ReduceScatter
( const byte* sbuf, byte* rbuf, const int* rcs, Op op, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::ReduceScatter");
#endif
    SafeMpi( 
        MPI_Reduce_scatter
        ( const_cast<byte*>(sbuf), 
          rbuf, const_cast<int*>(rcs), MPI_UNSIGNED_CHAR, op, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void ReduceScatter
( const int* sbuf, int* rbuf, const int* rcs, Op op, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::ReduceScatter");
#endif
    SafeMpi( 
        MPI_Reduce_scatter
        ( const_cast<int*>(sbuf), 
          rbuf, const_cast<int*>(rcs), MPI_INT, op, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void ReduceScatter
( const float* sbuf, float* rbuf, const int* rcs, Op op, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::ReduceScatter");
#endif
    SafeMpi( 
        MPI_Reduce_scatter
        ( const_cast<float*>(sbuf), 
          rbuf, const_cast<int*>(rcs), MPI_FLOAT, op, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void ReduceScatter
( const double* sbuf, double* rbuf, const int* rcs, Op op, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::ReduceScatter");
#endif
    SafeMpi( 
        MPI_Reduce_scatter
        ( const_cast<double*>(sbuf), 
          rbuf, const_cast<int*>(rcs), MPI_DOUBLE, op, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

void ReduceScatter
( const scomplex* sbuf, scomplex* rbuf, const int* rcs, Op op, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::ReduceScatter");
#endif
#ifdef AVOID_COMPLEX_MPI
    if( op == SUM )
    {
        int p;
        MPI_Comm_size( comm, &p );
        std::vector<int> rcsDoubled(p);
        for( int i=0; i<p; ++i )
            rcsDoubled[i] = 2*rcs[i];
        SafeMpi(
            MPI_Reduce_scatter
            ( const_cast<scomplex*>(sbuf),
              rbuf, &rcsDoubled[0], MPI_FLOAT, op, comm )
        );
    }
    else
    {
        SafeMpi(
            MPI_Reduce_scatter
            ( const_cast<scomplex*>(sbuf),
              rbuf, const_cast<int*>(rcs), MPI_COMPLEX, op, comm )
        );
    }
#else
    SafeMpi( 
        MPI_Reduce_scatter
        ( const_cast<scomplex*>(sbuf), 
          rbuf, const_cast<int*>(rcs), MPI_COMPLEX, op, comm ) 
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

void ReduceScatter
( const dcomplex* sbuf, dcomplex* rbuf, const int* rcs, Op op, Comm comm )
{
#ifndef RELEASE
    PushCallStack("mpi::ReduceScatter");
#endif
#ifdef AVOID_COMPLEX_MPI
    if( op == SUM )
    {
        int p;
        MPI_Comm_size( comm, &p );
        std::vector<int> rcsDoubled(p);
        for( int i=0; i<p; ++i )
            rcsDoubled[i] = 2*rcs[i];
        SafeMpi(
            MPI_Reduce_scatter
            ( const_cast<dcomplex*>(sbuf),
              rbuf, &rcsDoubled[0], MPI_DOUBLE, op, comm )
        );
    }
    else
    {
        SafeMpi(
            MPI_Reduce_scatter
            ( const_cast<dcomplex*>(sbuf),
              rbuf, const_cast<int*>(rcs), MPI_DOUBLE_COMPLEX, op, comm )
        );
    }
#else
    SafeMpi( 
        MPI_Reduce_scatter
        ( const_cast<dcomplex*>(sbuf), 
          rbuf, const_cast<int*>(rcs), MPI_DOUBLE_COMPLEX, op, comm ) 
    );
#endif
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace mpi
} // namespace elem
