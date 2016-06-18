/*
   Copyright (c) 2009-2016, Jack Poulson
                      2013, Jeff Hammond
                      2013, Jed Brown
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>

// TODO: Introduce macros to shorten the explicit instantiation code

typedef unsigned char* UCP;

namespace {

inline void 
SafeMpi( int mpiError ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(
      if( mpiError != MPI_SUCCESS )    
      {
          char errorString[MPI_MAX_ERROR_STRING];
          int lengthOfErrorString;
          MPI_Error_string( mpiError, errorString, &lengthOfErrorString );
          El::RuntimeError( std::string(errorString) );
      }
    )
}

} // anonymous namespace

namespace El {
namespace mpi {

bool CommSameSizeAsInteger() EL_NO_EXCEPT
{ return sizeof(MPI_Comm) == sizeof(int); }

bool GroupSameSizeAsInteger() EL_NO_EXCEPT
{ return sizeof(MPI_Group) == sizeof(int); }

// MPI environmental routines
// ==========================

void Initialize( int& argc, char**& argv ) EL_NO_EXCEPT
{ MPI_Init( &argc, &argv ); }

int InitializeThread( int& argc, char**& argv, int required ) EL_NO_EXCEPT
{ 
    int provided; 
#ifdef EL_HAVE_MPI_INIT_THREAD
    MPI_Init_thread( &argc, &argv, required, &provided ); 
#else
    MPI_Init( &argc, &argv );
    provided = 0; // equivalent to MPI_THREAD_SINGLE
#endif
    return provided;
}

void Finalize() EL_NO_EXCEPT
{ MPI_Finalize(); }

bool Initialized() EL_NO_EXCEPT
{ 
    int initialized;
    MPI_Initialized( &initialized );
    return initialized;
}

bool Finalized() EL_NO_EXCEPT
{
    int finalized;
    MPI_Finalized( &finalized );
    return finalized;
}

int QueryThread() EL_NO_EXCEPT
{
    int provided;
#ifdef EL_HAVE_MPI_QUERY_THREAD
    MPI_Query_thread( &provided );
#else
    provided = 0; // equivalent to MPI_THREAD_SINGLE
#endif
    return provided;
}

void Abort( Comm comm, int errCode ) EL_NO_EXCEPT
{ MPI_Abort( comm.comm, errCode ); }

double Time() EL_NO_EXCEPT { return MPI_Wtime(); }

void Create( UserFunction* func, bool commutes, Op& op ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    SafeMpi( MPI_Op_create( func, commutes, &op.op ) );
}

void Free( Op& op ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    SafeMpi( MPI_Op_free( &op.op ) );
}

void Free( Datatype& type ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    SafeMpi( MPI_Type_free( &type ) );
}

// Communicator manipulation 
// =========================
int Rank( Comm comm ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    if( comm != COMM_NULL )
    {
        int rank;
        SafeMpi( MPI_Comm_rank( comm.comm, &rank ) );
        return rank;
    }
    else return UNDEFINED;
}

int Size( Comm comm ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    if( comm != COMM_NULL )
    {
        int size;
        SafeMpi( MPI_Comm_size( comm.comm, &size ) );
        return size;
    } 
    else return UNDEFINED;
}

void Create( Comm parentComm, Group subsetGroup, Comm& subsetComm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    SafeMpi( 
        MPI_Comm_create( parentComm.comm, subsetGroup.group, &subsetComm.comm ) 
    );
}

void Dup( Comm original, Comm& duplicate ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    SafeMpi( MPI_Comm_dup( original.comm, &duplicate.comm ) );
}

void Split( Comm comm, int color, int key, Comm& newComm ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    SafeMpi( MPI_Comm_split( comm.comm, color, key, &newComm.comm ) );
}

void Free( Comm& comm ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    SafeMpi( MPI_Comm_free( &comm.comm ) );
}

bool Congruent( Comm comm1, Comm comm2 ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    int result;
    SafeMpi( MPI_Comm_compare( comm1.comm, comm2.comm, &result ) );
    return ( result == MPI_IDENT || result == MPI_CONGRUENT );
}

void ErrorHandlerSet( Comm comm, ErrorHandler errorHandler )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
#ifdef EL_HAVE_MPI_COMM_SET_ERRHANDLER
    SafeMpi( MPI_Comm_set_errhandler( comm.comm, errorHandler ) );
#else
    SafeMpi( MPI_Errhandler_set( comm.comm, errorHandler ) );
#endif
}

// Cartesian communicator routines 
// ===============================

void CartCreate
( Comm comm, int numDims, const int* dimensions, const int* periods, 
  bool reorder, Comm& cartComm ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    SafeMpi
    ( MPI_Cart_create
      ( comm.comm, numDims, const_cast<int*>(dimensions), 
        const_cast<int*>(periods), reorder, &cartComm.comm ) );
}

void CartSub( Comm comm, const int* remainingDims, Comm& subComm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    SafeMpi( 
      MPI_Cart_sub
      ( comm.comm, const_cast<int*>(remainingDims), &subComm.comm ) 
    );
}

// Group manipulation 
// ==================

int Rank( Group group ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    int rank;
    SafeMpi( MPI_Group_rank( group.group, &rank ) );
    return rank;
}

int Size( Group group ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    int size;
    SafeMpi( MPI_Group_size( group.group, &size ) );
    return size;
}

void CommGroup( Comm comm, Group& group ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    SafeMpi( MPI_Comm_group( comm.comm, &group.group ) );
}

void Dup( Group group, Group& newGroup ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    // For some reason, MPI_Group_dup does not exist
    Excl( group, 0, 0, newGroup ); 
}

void Union( Group groupA, Group groupB, Group& newGroup ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    SafeMpi( MPI_Group_union( groupA.group, groupB.group, &newGroup.group ) );
}

void Incl( Group group, int n, const int* ranks, Group& subGroup )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    SafeMpi( 
      MPI_Group_incl
      ( group.group, n, const_cast<int*>(ranks), &subGroup.group ) 
    );
}

void Excl( Group group, int n, const int* ranks, Group& subGroup )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    SafeMpi(
      MPI_Group_excl
      ( group.group, n, const_cast<int*>(ranks), &subGroup.group )
    );
}

void Difference( Group parent, Group subset, Group& complement )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    SafeMpi( 
      MPI_Group_difference( parent.group, subset.group, &complement.group ) 
    );
}

void Free( Group& group ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    SafeMpi( MPI_Group_free( &group.group ) );
}

bool Congruent( Group group1, Group group2 ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    int result;
    SafeMpi( MPI_Group_compare( group1.group, group2.group, &result ) );
    return ( result == MPI_IDENT );
}

// Rank translations
// =================

int Translate( Group origGroup, int origRank, Group newGroup )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    int newRank;
    Translate( origGroup, 1, &origRank, newGroup, &newRank );
    return newRank;
}

int Translate( Comm origComm, int origRank, Group newGroup )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    int newRank;
    Translate( origComm, 1, &origRank, newGroup, &newRank );
    return newRank;
}

int Translate( Group origGroup, int origRank, Comm newComm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    int newRank;
    Translate( origGroup, 1, &origRank, newComm, &newRank );
    return newRank;
}

int Translate( Comm origComm, int origRank, Comm newComm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    int newRank;
    Translate( origComm, 1, &origRank, newComm, &newRank );
    return newRank;
}

void Translate
( Group origGroup, int size, const int* origRanks, 
  Group newGroup,                  int* newRanks ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    SafeMpi
    ( MPI_Group_translate_ranks
      ( origGroup.group, size, const_cast<int*>(origRanks), 
        newGroup.group, newRanks ) );
}

void Translate
( Comm origComm,  int size, const int* origRanks, 
  Group newGroup,                 int* newRanks ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    Group origGroup;
    CommGroup( origComm, origGroup );
    Translate( origGroup, size, origRanks, newGroup, newRanks );
    Free( origGroup );
}

void Translate
( Group origGroup,  int size, const int* origRanks, 
  Comm newComm,                     int* newRanks ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    Group newGroup;
    CommGroup( newComm,  newGroup  );
    Translate( origGroup, size, origRanks, newGroup, newRanks );
    Free( newGroup  );
}

void Translate
( Comm origComm,  int size, const int* origRanks, 
  Comm newComm,                   int* newRanks ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    Group origGroup, newGroup;
    CommGroup( origComm, origGroup );
    CommGroup( newComm,  newGroup  );
    Translate( origGroup, size, origRanks, newGroup, newRanks );
    Free( origGroup );
    Free( newGroup  );
}

// Various utilities
// =================

// Wait until every process in comm reaches this statement
void Barrier( Comm comm ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    SafeMpi( MPI_Barrier( comm.comm ) );
}

// Test for completion
template<typename T>
bool Test( Request<T>& request ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    Status status;
    int flag;
    SafeMpi( MPI_Test( &request.backend, &flag, &status ) );
    return flag;
}

// Ensure that the request finishes before continuing
template<typename T>
void Wait( Request<T>& request ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    Status status;
    Wait( request, status );
}

// Ensure that the request finishes before continuing
template<typename T,typename>
void Wait( Request<T>& request, Status& status ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    SafeMpi( MPI_Wait( &request.backend, &status ) );
}

// Ensure that several requests finish before continuing
template<typename T>
void WaitAll( int numRequests, Request<T>* requests ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    vector<Status> statuses( numRequests );
    WaitAll( numRequests, requests, statuses.data() );
}

// Ensure that several requests finish before continuing
template<typename T,typename>
void WaitAll( int numRequests, Request<T>* requests, Status* statuses )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
#ifndef EL_MPI_REQUEST_IS_NOT_POINTER
    // Both MPICH and OpenMPI define MPI_Request to be a pointer to a structure,
    // which implies that the following code is legal. AFAIK, there are not
    // any popular MPI implementations which should break this logic, but 
    // the alternative #ifdef logic is provided in case a breakage is observed.
    vector<MPI_Request> backends( numRequests );
    for( Int j=0; j<numRequests; ++j )
        backends[j] = requests[j].backend;
    SafeMpi( MPI_Waitall( numRequests, backends.data(), statuses ) );
    // NOTE: This write back will almost always be superfluous, but it ensures
    //       that any changes to the pointer are propagated
    for( Int j=0; j<numRequests; ++j )
        requests[j].backend = backends[j];
#else
    for( Int j=0; j<numRequests; ++j )
    {
        Status status;
        MPI_Wait( &requests[j].backend, &status );
    }
#endif
}

template<typename T,typename,typename>
void Wait( Request<T>& request, Status& status ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    SafeMpi( MPI_Wait( &request.backend, &status ) );
    if( request.receivingPacked )
    {
        Deserialize
        ( request.recvCount, request.buffer.data(), request.unpackedRecvBuf );
        request.receivingPacked = false;
    }
    request.buffer.clear();
}

template<typename T,typename,typename>
void WaitAll( int numRequests, Request<T>* requests, Status* statuses )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
#ifndef EL_MPI_REQUEST_IS_NOT_POINTER
    // Both MPICH and OpenMPI define MPI_Request to be a pointer to a structure,
    // which implies that the following code is legal. AFAIK, there are not
    // any popular MPI implementations which should break this logic, but 
    // the alternative #ifdef logic is provided in case a breakage is observed.
    vector<MPI_Request> backends( numRequests );
    for( Int j=0; j<numRequests; ++j )
        backends[j] = requests[j].backend;
    SafeMpi( MPI_Waitall( numRequests, backends.data(), statuses ) );
    // NOTE: This write back will almost always be superfluous, but it ensures
    //       that any changes to the pointer are propagated
    for( Int j=0; j<numRequests; ++j )
        requests[j].backend = backends[j];
#else
    for( Int j=0; j<numRequests; ++j )
    {
        Status status;
        MPI_Wait( &requests[j].backend, &status );
    }
#endif
    for( Int j=0; j<numRequests; ++j )
    {
        if( requests[j].receivingPacked )
        {
            Deserialize
            ( requests[j].recvCount,
              requests[j].buffer.data(),
              requests[j].unpackedRecvBuf );
            requests[j].receivingPacked = false;
        }
        requests[j].buffer.clear();
    }
}

// Nonblocking test for message completion
bool IProbe( int source, int tag, Comm comm, Status& status )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    int flag;
    SafeMpi( MPI_Iprobe( source, tag, comm.comm, &flag, &status ) );
    return flag;
}
bool IProbe( int source, Comm comm, Status& status ) EL_NO_RELEASE_EXCEPT
{ return IProbe( source, 0, comm, status ); }

template<typename T>
int GetCount( Status& status ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    int count;
    SafeMpi( MPI_Get_count( &status, TypeMap<T>(), &count ) );
    return count;
}

template<typename Real,typename>
void TaggedSend( const Real* buf, int count, int to, int tag, Comm comm )
EL_NO_RELEASE_EXCEPT
{ 
    DEBUG_CSE
    SafeMpi
    ( MPI_Send
      ( const_cast<Real*>(buf), count, TypeMap<Real>(), to, tag, comm.comm ) );
}

template<typename Real,typename>
void TaggedSend
( const Complex<Real>* buf, int count, int to, int tag, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
#ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi
    ( MPI_Send
      ( const_cast<Complex<Real>*>(buf), 2*count, TypeMap<Real>(), to, 
        tag, comm.comm ) );
#else
    SafeMpi
    ( MPI_Send
      ( const_cast<Complex<Real>*>(buf), count, 
        TypeMap<Complex<Real>>(), to, tag, comm.comm ) );
#endif
}

template<typename T,typename,typename>
void TaggedSend( const T* buf, int count, int to, int tag, Comm comm )
{
    DEBUG_CSE
    std::vector<byte> packedBuf;
    Serialize( count, buf, packedBuf );
    SafeMpi
    ( MPI_Send( packedBuf.data(), count, TypeMap<T>(), to, tag, comm.comm ) );
}

template<typename T>
void Send( const T* buf, int count, int to, Comm comm ) EL_NO_RELEASE_EXCEPT
{ TaggedSend( buf, count, to, 0, comm ); }

template<typename T>
void TaggedSend( T b, int to, int tag, Comm comm ) EL_NO_RELEASE_EXCEPT
{ TaggedSend( &b, 1, to, tag, comm ); }

template<typename T>
void Send( T b, int to, Comm comm ) EL_NO_RELEASE_EXCEPT
{ TaggedSend( b, to, 0, comm ); }

template<typename Real,typename>
void TaggedISend
( const Real* buf, int count, int to, int tag, Comm comm,
  Request<Real>& request )
EL_NO_RELEASE_EXCEPT
{ 
    DEBUG_CSE
    SafeMpi
    ( MPI_Isend
      ( const_cast<Real*>(buf), count, TypeMap<Real>(), to, 
        tag, comm.comm, &request.backend ) );
}

template<typename Real,typename>
void TaggedISend
( const Complex<Real>* buf, int count, int to, int tag, Comm comm, 
  Request<Complex<Real>>& request ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
#ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi
    ( MPI_Isend
      ( const_cast<Complex<Real>*>(buf), 2*count, 
        TypeMap<Real>(), to, tag, comm.comm, &request.backend ) );
#else
    SafeMpi
    ( MPI_Isend
      ( const_cast<Complex<Real>*>(buf), count, 
        TypeMap<Complex<Real>>(), to, tag, comm.comm, &request.backend ) );
#endif
}

template<typename T,typename,typename>
void TaggedISend
( const T* buf, int count, int to, int tag, Comm comm, Request<T>& request )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    Serialize( count, buf, request.buffer );
    SafeMpi
    ( MPI_Isend
      ( request.buffer.data(), count, TypeMap<T>(), to, tag, comm.comm,
        &request.backend ) );
}

template<typename T>
void ISend
( const T* buf, int count, int to, Comm comm, Request<T>& request )
EL_NO_RELEASE_EXCEPT
{ TaggedISend( buf, count, to, 0, comm, request ); } 

template<typename T>
void TaggedISend( T b, int to, int tag, Comm comm, Request<T>& request )
EL_NO_RELEASE_EXCEPT
{ TaggedISend( &b, 1, to, tag, comm, request ); }

template<typename T>
void ISend( T b, int to, Comm comm, Request<T>& request )
EL_NO_RELEASE_EXCEPT
{ TaggedISend( b, to, 0, comm, request ); }

template<typename Real,typename>
void TaggedIRSend
( const Real* buf, int count, int to, int tag, Comm comm,
  Request<Real>& request )
EL_NO_RELEASE_EXCEPT
{ 
    DEBUG_CSE
    SafeMpi
    ( MPI_Irsend
      ( const_cast<Real*>(buf), count, TypeMap<Real>(), to, 
        tag, comm.comm, &request.backend ) );
}

template<typename Real,typename>
void TaggedIRSend
( const Complex<Real>* buf, int count, int to, int tag, Comm comm, 
  Request<Complex<Real>>& request ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
#ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi
    ( MPI_Irsend
      ( const_cast<Complex<Real>*>(buf), 2*count, 
        TypeMap<Real>(), to, tag, comm.comm, &request.backend ) );
#else
    SafeMpi
    ( MPI_Irsend
      ( const_cast<Complex<Real>*>(buf), count, 
        TypeMap<Complex<Real>>(), to, tag, comm.comm, &request.backend ) );
#endif
}

template<typename T,typename,typename>
void TaggedIRSend
( const T* buf, int count, int to, int tag, Comm comm,
  Request<T>& request )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    Serialize( count, buf, request.buffer );
    SafeMpi
    ( MPI_Irsend
      ( request.buffer.data(), count, TypeMap<T>(), to, 
        tag, comm.comm, &request.backend ) );
}

template<typename T>
void IRSend
( const T* buf, int count, int to, Comm comm, Request<T>& request )
EL_NO_RELEASE_EXCEPT
{ TaggedIRSend( buf, count, to, 0, comm, request ); } 

template<typename T>
void TaggedIRSend( T b, int to, int tag, Comm comm, Request<T>& request )
EL_NO_RELEASE_EXCEPT
{ TaggedIRSend( &b, 1, to, tag, comm, request ); }

template<typename T>
void IRSend( T b, int to, Comm comm, Request<T>& request )
EL_NO_RELEASE_EXCEPT
{ TaggedIRSend( b, to, 0, comm, request ); }

template<typename Real,typename>
void TaggedISSend
( const Real* buf, int count, int to, int tag, Comm comm,
  Request<Real>& request ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    SafeMpi
    ( MPI_Issend
      ( const_cast<Real*>(buf), count, TypeMap<Real>(), to, 
        tag, comm.comm, &request.backend ) );
}

template<typename Real,typename>
void TaggedISSend
( const Complex<Real>* buf, int count, int to, int tag, Comm comm, 
  Request<Complex<Real>>& request ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
#ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi
    ( MPI_Issend
      ( const_cast<Complex<Real>*>(buf), 2*count, 
        TypeMap<Real>(), to, tag, comm.comm, &request.backend ) );
#else
    SafeMpi
    ( MPI_Issend
      ( const_cast<Complex<Real>*>(buf), count, 
        TypeMap<Complex<Real>>(), to, tag, comm.comm, &request.backend ) );
#endif
}

template<typename T,typename,typename>
void TaggedISSend
( const T* buf, int count, int to, int tag, Comm comm, Request<T>& request )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    Serialize( count, buf, request.buffer );
    SafeMpi
    ( MPI_Issend
      ( request.buffer.data(), count, TypeMap<T>(), to, 
        tag, comm.comm, &request.backend ) );
}

template<typename T>
void ISSend( const T* buf, int count, int to, Comm comm, Request<T>& request )
EL_NO_RELEASE_EXCEPT
{ TaggedISSend( buf, count, to, 0, comm, request ); }

template<typename T>
void TaggedISSend( T b, int to, int tag, Comm comm, Request<T>& request )
EL_NO_RELEASE_EXCEPT
{ TaggedISSend( &b, 1, to, tag, comm, request ); }

template<typename Real,typename>
void TaggedRecv( Real* buf, int count, int from, int tag, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    Status status;
    SafeMpi
    ( MPI_Recv( buf, count, TypeMap<Real>(), from, tag, comm.comm, &status ) );
}

template<typename Real,typename>
void TaggedRecv( Complex<Real>* buf, int count, int from, int tag, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    Status status;
#ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi
    ( MPI_Recv( buf, 2*count, TypeMap<Real>(), from, tag, comm.comm, &status ) );
#else
    SafeMpi
    ( MPI_Recv
      ( buf, count, TypeMap<Complex<Real>>(), from, tag, comm.comm, &status ) );
#endif
}

template<typename T,typename,typename>
void TaggedRecv( T* buf, int count, int from, int tag, Comm comm )
{
    DEBUG_CSE
    std::vector<byte> packedBuf;
    ReserveSerialized( count, buf, packedBuf );
    Status status;
    SafeMpi
    ( MPI_Recv
      ( packedBuf.data(), count, TypeMap<T>(), from, tag,
        comm.comm, &status ) );
    Deserialize( count, packedBuf, buf );
}

template<typename T>
void Recv( T* buf, int count, int from, Comm comm )
EL_NO_RELEASE_EXCEPT
{ TaggedRecv( buf, count, from, ANY_TAG, comm ); }

template<typename T>
T TaggedRecv( int from, int tag, Comm comm )
EL_NO_RELEASE_EXCEPT
{ T b; TaggedRecv( &b, 1, from, tag, comm ); return b; }

template<typename T>
T Recv( int from, Comm comm )
EL_NO_RELEASE_EXCEPT
{ return TaggedRecv<T>( from, ANY_TAG, comm ); }

template<typename Real,typename>
void TaggedIRecv
( Real* buf, int count, int from, int tag, Comm comm, Request<Real>& request )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    SafeMpi
    ( MPI_Irecv
      ( buf, count, TypeMap<Real>(), from, tag, comm.comm, &request.backend ) );
}

template<typename Real,typename>
void TaggedIRecv
( Complex<Real>* buf, int count, int from, int tag, Comm comm,
  Request<Complex<Real>>& request )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
#ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi
    ( MPI_Irecv
      ( buf, 2*count, TypeMap<Real>(), from, tag, comm.comm,
        &request.backend ) );
#else
    SafeMpi
    ( MPI_Irecv
      ( buf, count, TypeMap<Complex<Real>>(), from, tag, comm.comm,
        &request.backend ) );
#endif
}

template<typename T,typename,typename>
void TaggedIRecv
( T* buf, int count, int from, int tag, Comm comm, Request<T>& request )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    request.receivingPacked = true;
    request.recvCount = count;
    request.unpackedRecvBuf = buf;
    ReserveSerialized( count, buf, request.buffer );
    SafeMpi
    ( MPI_Irecv
      ( request.buffer.data(), count, TypeMap<T>(), from, tag, comm.comm,
        &request.backend ) );
}

template<typename T>
void IRecv( T* buf, int count, int from, Comm comm, Request<T>& request )
EL_NO_RELEASE_EXCEPT
{ TaggedIRecv( buf, count, from, ANY_TAG, comm, request ); }

template<typename T>
T TaggedIRecv( int from, int tag, Comm comm, Request<T>& request )
EL_NO_RELEASE_EXCEPT
{ T b; TaggedIRecv( &b, 1, from, tag, comm, request ); return b; }

template<typename T>
T IRecv( int from, Comm comm, Request<T>& request )
EL_NO_RELEASE_EXCEPT
{ return TaggedIRecv<T>( from, ANY_TAG, comm, request ); }

template<typename Real,typename>
void TaggedSendRecv
( const Real* sbuf, int sc, int to,   int stag,
        Real* rbuf, int rc, int from, int rtag, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    Status status;
    SafeMpi
    ( MPI_Sendrecv
      ( const_cast<Real*>(sbuf), sc, TypeMap<Real>(), to,   stag,
        rbuf,                    rc, TypeMap<Real>(), from, rtag, 
        comm.comm, &status ) );
}

template<typename Real,typename>
void TaggedSendRecv
( const Complex<Real>* sbuf, int sc, int to,   int stag,
        Complex<Real>* rbuf, int rc, int from, int rtag, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    Status status;
#ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi
    ( MPI_Sendrecv
      ( const_cast<Complex<Real>*>(sbuf), 2*sc, TypeMap<Real>(), to,   stag,
        rbuf,                             2*rc, TypeMap<Real>(), from, rtag, 
        comm.comm, &status ) );
#else
    SafeMpi
    ( MPI_Sendrecv
      ( const_cast<Complex<Real>*>(sbuf), 
        sc, TypeMap<Complex<Real>>(), to,   stag,
        rbuf,                          
        rc, TypeMap<Complex<Real>>(), from, rtag, comm.comm, &status ) );
#endif
}

template<typename T,typename,typename>
void TaggedSendRecv
( const T* sbuf, int sc, int to,   int stag,
        T* rbuf, int rc, int from, int rtag, Comm comm )
{
    DEBUG_CSE
    Status status;
    std::vector<byte> packedSend, packedRecv;
    Serialize( sc, sbuf, packedSend );
    ReserveSerialized( rc, rbuf, packedRecv );
    SafeMpi
    ( MPI_Sendrecv
      ( packedSend.data(), sc, TypeMap<T>(), to,   stag,
        packedRecv.data(), rc, TypeMap<T>(), from, rtag, 
        comm.comm, &status ) );
    Deserialize( rc, packedRecv, rbuf );
}

template<typename T>
void SendRecv
( const T* sbuf, int sc, int to, 
        T* rbuf, int rc, int from, Comm comm )
EL_NO_RELEASE_EXCEPT
{ TaggedSendRecv( sbuf, sc, to, 0, rbuf, rc, from, ANY_TAG, comm ); }

template<typename T>
T TaggedSendRecv( T sb, int to, int stag, int from, int rtag, Comm comm )
EL_NO_RELEASE_EXCEPT
{ 
    T rb; 
    TaggedSendRecv( &sb, 1, to, stag, &rb, 1, from, rtag, comm ); 
    return rb; 
}

template<typename T>
T SendRecv( T sb, int to, int from, Comm comm )
EL_NO_RELEASE_EXCEPT
{ return TaggedSendRecv( sb, to, 0, from, ANY_TAG, comm ); }

template<typename Real,typename>
void TaggedSendRecv
( Real* buf, int count, int to, int stag, int from, int rtag, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    Status status;
    SafeMpi
    ( MPI_Sendrecv_replace
      ( buf, count, TypeMap<Real>(), to, stag, from, rtag, comm.comm,
        &status ) );
}

template<typename Real,typename>
void TaggedSendRecv
( Complex<Real>* buf, int count, int to, int stag, int from, int rtag,
  Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    Status status;
#ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi
    ( MPI_Sendrecv_replace
      ( buf, 2*count, TypeMap<Real>(), to, stag, from, rtag, comm.comm, 
        &status ) );
#else
    SafeMpi
    ( MPI_Sendrecv_replace
      ( buf, count, TypeMap<Complex<Real>>(), 
        to, stag, from, rtag, comm.comm, &status ) );
#endif
}

template<typename T,typename,typename>
void TaggedSendRecv
( T* buf, int count, int to, int stag, int from, int rtag, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    std::vector<byte> packedBuf;
    ReserveSerialized( count, buf, packedBuf );
    Serialize( count, buf, packedBuf );
    Status status;
    SafeMpi
    ( MPI_Sendrecv_replace
      ( packedBuf.data(), count, TypeMap<T>(), to, stag, from, rtag,
        comm.comm, &status ) );
    Deserialize( count, packedBuf, buf );
}

template<typename T>
void SendRecv( T* buf, int count, int to, int from, Comm comm )
EL_NO_RELEASE_EXCEPT
{ TaggedSendRecv( buf, count, to, 0, from, ANY_TAG, comm ); }

template<typename Real,typename>
void Broadcast( Real* buf, int count, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    if( Size(comm) == 1 || count == 0 )
        return;
    SafeMpi( MPI_Bcast( buf, count, TypeMap<Real>(), root, comm.comm ) );
}

template<typename Real,typename>
void Broadcast( Complex<Real>* buf, int count, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    if( Size(comm) == 1 )
        return;
#ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi( MPI_Bcast( buf, 2*count, TypeMap<Real>(), root, comm.comm ) );
#else
    SafeMpi( MPI_Bcast( buf, count, TypeMap<Complex<Real>>(), root, comm.comm ) );
#endif
}

template<typename T,typename,typename>
void Broadcast( T* buf, int count, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    if( Size(comm) == 1 || count == 0 )
        return;
    std::vector<byte> packedBuf;
    Serialize( count, buf, packedBuf );
    SafeMpi(
      MPI_Bcast( packedBuf.data(), count, TypeMap<T>(), root, comm.comm )
    );
    Deserialize( count, packedBuf, buf );
}

template<typename T>
void Broadcast( T& b, int root, Comm comm ) EL_NO_RELEASE_EXCEPT
{ Broadcast( &b, 1, root, comm ); }

template<typename Real,typename>
void IBroadcast
( Real* buf, int count, int root, Comm comm, Request<Real>& request )
{
    DEBUG_CSE
#ifdef EL_HAVE_NONBLOCKING_COLLECTIVES
    SafeMpi
    ( MPI_Ibcast
      ( buf, count, TypeMap<Real>(), root, comm.comm, &request.backend ) );
#else
    LogicError("Elemental was not configured with non-blocking support");
#endif
}

template<typename Real,typename>
void IBroadcast
( Complex<Real>* buf, int count, int root, Comm comm,
  Request<Complex<Real>>& request )
{
    DEBUG_CSE
#ifdef EL_HAVE_NONBLOCKING_COLLECTIVES
#ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi
    ( MPI_Ibcast
      ( buf, 2*count, TypeMap<Real>(), root, comm.comm, &request.backend ) );
#else
    SafeMpi
    ( MPI_Ibcast
      ( buf, count, TypeMap<Complex<Real>>(), root, comm.comm,
        &request.backend ) );
#endif
#else
    LogicError("Elemental was not configured with non-blocking support");
#endif
}

template<typename T,typename,typename>
void IBroadcast
( T* buf, int count, int root, Comm comm, Request<T>& request )
{
    DEBUG_CSE
#ifdef EL_HAVE_NONBLOCKING_COLLECTIVES
    request.receivingPacked = true;
    request.recvCount = count;
    request.unpackedRecvBuf = buf;
    ReserveSerialized( count, buf, request.buffer );
    SafeMpi
    ( MPI_Ibcast
      ( request.buffer.data(), count, TypeMap<Real>(), root, comm.comm,
        &request.backend ) );
#else
    LogicError("Elemental was not configured with non-blocking support");
#endif
}

template<typename T>
void IBroadcast( T& b, int root, Comm comm, Request<T>& request )
{ IBroadcast( &b, 1, root, comm, request ); }

template<typename Real,typename>
void Gather
( const Real* sbuf, int sc,
        Real* rbuf, int rc, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    SafeMpi
    ( MPI_Gather
      ( const_cast<Real*>(sbuf), sc, TypeMap<Real>(),
        rbuf,                    rc, TypeMap<Real>(), root, comm.comm ) );
}

template<typename Real,typename>
void Gather
( const Complex<Real>* sbuf, int sc,
        Complex<Real>* rbuf, int rc, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
#ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi
    ( MPI_Gather
      ( const_cast<Complex<Real>*>(sbuf), 2*sc, TypeMap<Real>(),
        rbuf,                             2*rc, TypeMap<Real>(),
        root, comm.comm ) );
#else
    SafeMpi
    ( MPI_Gather
      ( const_cast<Complex<Real>*>(sbuf), sc, TypeMap<Complex<Real>>(),
        rbuf,                             rc, TypeMap<Complex<Real>>(), 
        root, comm.comm ) );
#endif
}

template<typename T,typename,typename>
void Gather
( const T* sbuf, int sc,
        T* rbuf, int rc, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    const int commSize = mpi::Size(comm);
    const int commRank = mpi::Rank(comm);
    const int totalRecv = rc*commSize;

    std::vector<byte> packedSend, packedRecv;
    Serialize( sc, sbuf, packedSend );

    if( commRank == root )
        ReserveSerialized( totalRecv, rbuf, packedRecv );
    SafeMpi
    ( MPI_Gather
      ( packedSend.data(), sc, TypeMap<T>(),
        packedRecv.data(), rc, TypeMap<T>(), root, comm.comm ) );
    if( commRank == root )
        Deserialize( totalRecv, packedRecv, rbuf );
}

template<typename Real,typename>
void IGather
( const Real* sbuf, int sc,
        Real* rbuf, int rc,
  int root, Comm comm,
  Request<Real>& request )
{
    DEBUG_CSE
#ifdef EL_HAVE_NONBLOCKING_COLLECTIVES
    SafeMpi
    ( MPI_Igather
      ( const_cast<Real*>(sbuf), sc, TypeMap<Real>(),
        rbuf,                    rc, TypeMap<Real>(), root, comm.comm,
        &request.backend ) );
#else
    LogicError("Elemental was not configured with non-blocking support");
#endif
}

template<typename Real,typename>
void IGather
( const Complex<Real>* sbuf, int sc,
        Complex<Real>* rbuf, int rc,
  int root, Comm comm,
  Request<Complex<Real>>& request )
{
    DEBUG_CSE
#ifdef EL_HAVE_NONBLOCKING_COLLECTIVES
#ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi
    ( MPI_Igather
      ( const_cast<Complex<Real>*>(sbuf), 2*sc, TypeMap<Real>(),
        rbuf,                             2*rc, TypeMap<Real>(), 
        root, comm.comm, &request.backend ) );
#else
    SafeMpi
    ( MPI_Igather
      ( const_cast<Complex<Real>*>(sbuf), sc, TypeMap<Complex<Real>>(),
        rbuf,                             rc, TypeMap<Complex<Real>>(), 
        root, comm.comm, &request.backend ) );
#endif
#else
    LogicError("Elemental was not configured with non-blocking support");
#endif
}

template<typename T,typename,typename>
void IGather
( const T* sbuf, int sc,
        T* rbuf, int rc,
  int root, Comm comm,
  Request<T>& request )
{
    DEBUG_CSE
#ifdef EL_HAVE_NONBLOCKING_COLLECTIVES
    if( mpi::Rank(comm) == root )
    {
        const int commSize = mpi::Size(comm);
        request.receivingPacked = true;
        request.recvCount = rc*commSize;
        request.unpackedRecvBuf = rbuf;
        ReserveSerialized( rc*commSize, rbuf, request.buffer );
    }
    SafeMpi
    ( MPI_Igather
      ( request.buffer.data(), sc, TypeMap<Real>(),
        rbuf,                  rc, TypeMap<Real>(), root, comm.comm,
        &request.backend ) );
#else
    LogicError("Elemental was not configured with non-blocking support");
#endif
}

template<typename Real,typename>
void Gather
( const Real* sbuf, int sc,
        Real* rbuf, const int* rcs, const int* rds,
  int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    SafeMpi
    ( MPI_Gatherv
      ( const_cast<Real*>(sbuf), 
        sc,       
        TypeMap<Real>(),
        rbuf,                    
        const_cast<int*>(rcs), 
        const_cast<int*>(rds), 
        TypeMap<Real>(),
        root, 
        comm.comm ) );
}

template<typename Real,typename>
void Gather
( const Complex<Real>* sbuf, int sc,
        Complex<Real>* rbuf, const int* rcs, const int* rds, int root,
  Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
#ifdef EL_AVOID_COMPLEX_MPI
    const int commRank = Rank( comm );
    const int commSize = Size( comm );
    vector<int> rcsDouble, rdsDouble;
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
      ( const_cast<Complex<Real>*>(sbuf), 2*sc, TypeMap<Real>(),
        rbuf, rcsDouble.data(), rdsDouble.data(), TypeMap<Real>(),
        root, comm.comm ) );
#else
    SafeMpi
    ( MPI_Gatherv
      ( const_cast<Complex<Real>*>(sbuf), 
        sc,       
        TypeMap<Complex<Real>>(),
        rbuf,  
        const_cast<int*>(rcs), 
        const_cast<int*>(rds), 
        TypeMap<Complex<Real>>(),
        root, 
        comm.comm ) );
#endif
}

template<typename T,typename,typename>
void Gather
( const T* sbuf, int sc,
        T* rbuf, const int* rcs, const int* rds,
  int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    const int commSize = mpi::Size(comm);
    const int commRank = mpi::Rank(comm);
    int totalRecv=0;
    if( commRank == root )
        totalRecv = rcs[commSize-1]+rds[commSize-1];

    std::vector<byte> packedSend, packedRecv;
    Serialize( sc, sbuf, packedSend );

    if( commRank == root )
        ReserveSerialized( totalRecv, rbuf, packedRecv );
    SafeMpi
    ( MPI_Gatherv
      ( packedSend.data(),
        sc,
        TypeMap<T>(),
        packedRecv.data(),
        const_cast<int*>(rcs),
        const_cast<int*>(rds),
        TypeMap<T>(),
        root,
        comm.comm ) );
    if( commRank == root )
        Deserialize( totalRecv, packedRecv, rbuf );
}

template<typename Real,typename>
void AllGather
( const Real* sbuf, int sc,
        Real* rbuf, int rc, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
#ifdef EL_USE_BYTE_ALLGATHERS
    SafeMpi
    ( MPI_Allgather
      ( (UCP)const_cast<Real*>(sbuf), sizeof(Real)*sc, MPI_UNSIGNED_CHAR, 
        (UCP)rbuf,                    sizeof(Real)*rc, MPI_UNSIGNED_CHAR, 
        comm.comm ) );
#else
    SafeMpi
    ( MPI_Allgather
      ( const_cast<Real*>(sbuf), sc, TypeMap<Real>(), 
        rbuf,                    rc, TypeMap<Real>(), comm.comm ) );
#endif
}

template<typename Real,typename>
void AllGather
( const Complex<Real>* sbuf, int sc,
        Complex<Real>* rbuf, int rc, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
#ifdef EL_USE_BYTE_ALLGATHERS
    SafeMpi
    ( MPI_Allgather
      ( (UCP)const_cast<Complex<Real>*>(sbuf),
        2*sizeof(Real)*sc, MPI_UNSIGNED_CHAR, 
        (UCP)rbuf,
        2*sizeof(Real)*rc, MPI_UNSIGNED_CHAR, 
        comm.comm ) );
#else
 #ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi
    ( MPI_Allgather
      ( const_cast<Complex<Real>*>(sbuf), 2*sc, TypeMap<Real>(),
        rbuf,                             2*rc, TypeMap<Real>(),
        comm.comm ) );
 #else
    SafeMpi
    ( MPI_Allgather
      ( const_cast<Complex<Real>*>(sbuf), sc, TypeMap<Complex<Real>>(),
        rbuf,                             rc, TypeMap<Complex<Real>>(),
        comm.comm ) );
 #endif
#endif
}

template<typename T,typename,typename>
void AllGather
( const T* sbuf, int sc,
        T* rbuf, int rc, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    const int commSize = mpi::Size(comm);
    const int totalRecv = rc*commSize;

    std::vector<byte> packedSend, packedRecv;
    Serialize( sc, sbuf, packedSend );

    ReserveSerialized( totalRecv, rbuf, packedRecv );
    SafeMpi
    ( MPI_Allgather
      ( packedSend.data(), sc, TypeMap<T>(),
        packedRecv.data(), rc, TypeMap<T>(), comm.comm ) );
    Deserialize( totalRecv, packedRecv, rbuf );
}

template<typename Real,typename>
void AllGather
( const Real* sbuf, int sc,
        Real* rbuf, const int* rcs, const int* rds, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
#ifdef EL_USE_BYTE_ALLGATHERS
    const int commSize = Size( comm );
    vector<int> byteRcs( commSize ), byteRds( commSize );
    for( int i=0; i<commSize; ++i )
    {
        byteRcs[i] = sizeof(Real)*rcs[i];
        byteRds[i] = sizeof(Real)*rds[i];
    }
    SafeMpi
    ( MPI_Allgatherv
      ( (UCP)const_cast<Real*>(sbuf), sizeof(Real)*sc,   MPI_UNSIGNED_CHAR, 
        (UCP)rbuf, byteRcs.data(), byteRds.data(), MPI_UNSIGNED_CHAR, 
        comm.comm ) );
#else
    SafeMpi
    ( MPI_Allgatherv
      ( const_cast<Real*>(sbuf), 
        sc, 
        TypeMap<Real>(), 
        rbuf,   
        const_cast<int*>(rcs), 
        const_cast<int*>(rds), 
        TypeMap<Real>(), 
        comm.comm ) );
#endif
}

template<typename Real,typename>
void AllGather
( const Complex<Real>* sbuf, int sc,
        Complex<Real>* rbuf, const int* rcs, const int* rds, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
#ifdef EL_USE_BYTE_ALLGATHERS
    const int commSize = Size( comm );
    vector<int> byteRcs( commSize ), byteRds( commSize );
    for( int i=0; i<commSize; ++i )
    {
        byteRcs[i] = 2*sizeof(Real)*rcs[i];
        byteRds[i] = 2*sizeof(Real)*rds[i];
    }
    SafeMpi
    ( MPI_Allgatherv
      ( (UCP)const_cast<Complex<Real>*>(sbuf),
        2*sizeof(Real)*sc, MPI_UNSIGNED_CHAR, 
        (UCP)rbuf, byteRcs.data(), byteRds.data(),
        MPI_UNSIGNED_CHAR, 
        comm.comm ) );
#else
 #ifdef EL_AVOID_COMPLEX_MPI
    const int commSize = Size( comm );
    vector<int> realRcs( commSize ), realRds( commSize );
    for( int i=0; i<commSize; ++i )
    {
        realRcs[i] = 2*rcs[i];
        realRds[i] = 2*rds[i];
    }
    SafeMpi
    ( MPI_Allgatherv
      ( const_cast<Complex<Real>*>(sbuf), 2*sc, TypeMap<Real>(),
        rbuf, realRcs.data(), realRds.data(), TypeMap<Real>(), comm.comm ) );
 #else
    SafeMpi
    ( MPI_Allgatherv
      ( const_cast<Complex<Real>*>(sbuf), 
        sc, 
        TypeMap<Complex<Real>>(),
        rbuf, 
        const_cast<int*>(rcs), 
        const_cast<int*>(rds), 
        TypeMap<Complex<Real>>(),
        comm.comm ) );
 #endif
#endif
}

template<typename T,typename,typename>
void AllGather
( const T* sbuf, int sc,
        T* rbuf, const int* rcs, const int* rds, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    const int commSize = mpi::Size(comm);
    const int totalRecv = rcs[commSize-1]+rds[commSize-1];

    std::vector<byte> packedSend, packedRecv;
    Serialize( sc, sbuf, packedSend );

    ReserveSerialized( totalRecv, rbuf, packedRecv );
    SafeMpi
    ( MPI_Allgatherv
      ( packedSend.data(),
        sc,
        TypeMap<T>(),
        packedRecv.data(),
        const_cast<int*>(rcs),
        const_cast<int*>(rds),
        TypeMap<T>(),
        comm.comm ) );
    Deserialize( totalRecv, packedRecv, rbuf );
}

template<typename Real,typename>
void Scatter
( const Real* sbuf, int sc,
        Real* rbuf, int rc, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    SafeMpi
    ( MPI_Scatter
      ( const_cast<Real*>(sbuf), sc, TypeMap<Real>(),
        rbuf,                    rc, TypeMap<Real>(), root, comm.comm ) );
}

template<typename Real,typename>
void Scatter
( const Complex<Real>* sbuf, int sc,
        Complex<Real>* rbuf, int rc, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
#ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi
    ( MPI_Scatter
      ( const_cast<Complex<Real>*>(sbuf), 2*sc, TypeMap<Real>(),
        rbuf,                             2*rc, TypeMap<Real>(), root,
        comm.comm ) );
#else
    SafeMpi
    ( MPI_Scatter
      ( const_cast<Complex<Real>*>(sbuf), sc, TypeMap<Complex<Real>>(),
        rbuf,                             rc, TypeMap<Complex<Real>>(), 
        root, comm.comm ) );
#endif
}

template<typename T,typename,typename>
void Scatter
( const T* sbuf, int sc,
        T* rbuf, int rc, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    const int commSize = mpi::Size(comm);
    const int commRank = mpi::Rank(comm);
    const int totalSend = sc*commSize;

    std::vector<byte> packedSend, packedRecv;
    if( commRank == root )
        Serialize( totalSend, sbuf, packedSend );

    ReserveSerialized( rc, rbuf, packedRecv );
    SafeMpi
    ( MPI_Scatter
      ( packedSend.data(), sc, TypeMap<T>(),
        packedRecv.data(), rc, TypeMap<T>(), root, comm.comm ) );
    Deserialize( rc, packedRecv, rbuf );
}

template<typename Real,typename>
void Scatter( Real* buf, int sc, int rc, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    const int commRank = Rank( comm );
    if( commRank == root )
    {
        SafeMpi
        ( MPI_Scatter
          ( buf,          sc, TypeMap<Real>(), 
            MPI_IN_PLACE, rc, TypeMap<Real>(), root, comm.comm ) );
    }
    else
    {
        SafeMpi
        ( MPI_Scatter
          ( 0,   sc, TypeMap<Real>(), 
            buf, rc, TypeMap<Real>(), root, comm.comm ) );
    }
}

template<typename Real,typename>
void Scatter( Complex<Real>* buf, int sc, int rc, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    const int commRank = Rank( comm );
    if( commRank == root )
    {
#ifdef EL_AVOID_COMPLEX_MPI
        SafeMpi
        ( MPI_Scatter
          ( buf,          2*sc, TypeMap<Real>(), 
            MPI_IN_PLACE, 2*rc, TypeMap<Real>(), root, comm.comm ) );
#else
        SafeMpi
        ( MPI_Scatter
          ( buf,          sc, TypeMap<Complex<Real>>(), 
            MPI_IN_PLACE, rc, TypeMap<Complex<Real>>(), root, comm.comm ) );
#endif
    }
    else
    {
#ifdef EL_AVOID_COMPLEX_MPI
        SafeMpi
        ( MPI_Scatter
          ( 0,   2*sc, TypeMap<Real>(), 
            buf, 2*rc, TypeMap<Real>(), root, comm.comm ) );
#else
        SafeMpi
        ( MPI_Scatter
          ( 0,   sc, TypeMap<Complex<Real>>(), 
            buf, rc, TypeMap<Complex<Real>>(), root, comm.comm ) );
#endif
    }
}

template<typename T,typename,typename>
void Scatter( T* buf, int sc, int rc, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    const int commSize = mpi::Size(comm);
    const int commRank = mpi::Rank(comm);
    const int totalSend = sc*commSize;

    // TODO: Use in-place option?

    std::vector<byte> packedSend, packedRecv;
    if( commRank == root )
        Serialize( totalSend, buf, packedSend );

    ReserveSerialized( rc, buf, packedRecv );
    SafeMpi
    ( MPI_Scatter
      ( packedSend.data(), sc, TypeMap<T>(),
        packedRecv.data(), rc, TypeMap<T>(), root, comm.comm ) );
    Deserialize( rc, packedRecv, buf );
}

template<typename Real,typename>
void AllToAll
( const Real* sbuf, int sc,
        Real* rbuf, int rc, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    SafeMpi
    ( MPI_Alltoall
      ( const_cast<Real*>(sbuf), sc, TypeMap<Real>(),
        rbuf,                    rc, TypeMap<Real>(), comm.comm ) );
}

template<typename Real,typename>
void AllToAll
( const Complex<Real>* sbuf, int sc,
        Complex<Real>* rbuf, int rc, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
#ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi
    ( MPI_Alltoall
      ( const_cast<Complex<Real>*>(sbuf),
        2*sc, TypeMap<Real>(),
        rbuf,
        2*rc, TypeMap<Real>(), comm.comm ) );
#else
    SafeMpi
    ( MPI_Alltoall
      ( const_cast<Complex<Real>*>(sbuf),
        sc, TypeMap<Complex<Real>>(),
        rbuf,
        rc, TypeMap<Complex<Real>>(), comm.comm ) );
#endif
}

template<typename T,typename,typename>
void AllToAll
( const T* sbuf, int sc,
        T* rbuf, int rc, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    const int commSize = mpi::Size( comm );
    const int totalSend = sc*commSize;
    const int totalRecv = rc*commSize;

    std::vector<byte> packedSend, packedRecv;
    Serialize( totalSend, sbuf, packedSend );
    ReserveSerialized( totalRecv, rbuf, packedRecv );
    SafeMpi
    ( MPI_Alltoall
      ( packedSend.data(), sc, TypeMap<T>(),
        packedRecv.data(), rc, TypeMap<T>(), comm.comm ) );
    Deserialize( totalRecv, packedRecv, rbuf );
}

template<typename Real,typename>
void AllToAll
( const Real* sbuf, const int* scs, const int* sds, 
        Real* rbuf, const int* rcs, const int* rds, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    SafeMpi
    ( MPI_Alltoallv
      ( const_cast<Real*>(sbuf), 
        const_cast<int*>(scs), 
        const_cast<int*>(sds), 
        TypeMap<Real>(),
        rbuf, 
        const_cast<int*>(rcs), 
        const_cast<int*>(rds), 
        TypeMap<Real>(),
        comm.comm ) );
}

template<typename Real,typename>
void AllToAll
( const Complex<Real>* sbuf, const int* scs, const int* sds,
        Complex<Real>* rbuf, const int* rcs, const int* rds, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
#ifdef EL_AVOID_COMPLEX_MPI
    int p;
    MPI_Comm_size( comm.comm, &p );
    vector<int> scsDoubled(p), sdsDoubled(p),
                rcsDoubled(p), rdsDoubled(p);
    for( int i=0; i<p; ++i )
    {
        scsDoubled[i] = 2*scs[i];
        sdsDoubled[i] = 2*sds[i];
        rcsDoubled[i] = 2*rcs[i];
        rdsDoubled[i] = 2*rds[i];
    }
    SafeMpi
    ( MPI_Alltoallv
      ( const_cast<Complex<Real>*>(sbuf),
              scsDoubled.data(), sdsDoubled.data(), TypeMap<Real>(),
        rbuf, rcsDoubled.data(), rdsDoubled.data(), TypeMap<Real>(), comm.comm ) );
#else
    SafeMpi
    ( MPI_Alltoallv
      ( const_cast<Complex<Real>*>(sbuf), 
        const_cast<int*>(scs), 
        const_cast<int*>(sds), 
        TypeMap<Complex<Real>>(),
        rbuf, 
        const_cast<int*>(rcs), 
        const_cast<int*>(rds), 
        TypeMap<Complex<Real>>(),
        comm.comm ) );
#endif
}

template<typename T,typename,typename>
void AllToAll
( const T* sbuf, const int* scs, const int* sds,
        T* rbuf, const int* rcs, const int* rds, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    const int commSize = mpi::Size( comm );
    const int totalSend = scs[commSize-1]+sds[commSize-1];
    const int totalRecv = rcs[commSize-1]+rds[commSize-1];

    std::vector<byte> packedSend, packedRecv;
    Serialize( totalSend, sbuf, packedSend );
    ReserveSerialized( totalRecv, rbuf, packedRecv );
    SafeMpi
    ( MPI_Alltoallv
      ( packedSend.data(),
        const_cast<int*>(scs), const_cast<int*>(sds), TypeMap<T>(),
        packedRecv.data(),
        const_cast<int*>(rcs), const_cast<int*>(rds), TypeMap<T>(),
        comm.comm ) );
    Deserialize( totalRecv, packedRecv, rbuf );
}

template<typename T>
vector<T> AllToAll
( const vector<T>& sendBuf,
  const vector<int>& sendCounts, 
  const vector<int>& sendOffs,
  Comm comm )
EL_NO_RELEASE_EXCEPT
{
    const int commSize = Size( comm ); 
    vector<int> recvCounts(commSize);
    AllToAll( sendCounts.data(), 1, recvCounts.data(), 1, comm ); 
    vector<int> recvOffs;
    const int totalRecv = El::Scan( recvCounts, recvOffs );
    vector<T> recvBuf(totalRecv);
    AllToAll
    ( sendBuf.data(), sendCounts.data(), sendOffs.data(),
      recvBuf.data(), recvCounts.data(), recvOffs.data(), comm );
    return recvBuf;
}

template<typename Real,typename>
void Reduce
( const Real* sbuf, Real* rbuf, int count, Op op, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    if( count == 0 )
        return;

    MPI_Op opC;
    if( op == SUM )
        opC = SumOp<Real>().op; 
    else if( op == MAX )
        opC = MaxOp<Real>().op;
    else if( op == MIN )
        opC = MinOp<Real>().op;
    else
        opC = op.op;

    SafeMpi
    ( MPI_Reduce
      ( const_cast<Real*>(sbuf), rbuf, count, TypeMap<Real>(),
        opC, root, comm.comm ) );
}

template<typename Real,typename>
void Reduce
( const Complex<Real>* sbuf, 
        Complex<Real>* rbuf, int count, Op op, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    if( count == 0 )
        return;

#ifdef EL_AVOID_COMPLEX_MPI
    if( op == SUM )
    {
        MPI_Op opC = SumOp<Real>().op;
        SafeMpi
        ( MPI_Reduce
          ( const_cast<Complex<Real>*>(sbuf),
            rbuf, 2*count, TypeMap<Real>(), opC, 
            root, comm.comm ) );
    }
    else
    {
        MPI_Op opC = op.op;
        SafeMpi
        ( MPI_Reduce
          ( const_cast<Complex<Real>*>(sbuf),
            rbuf, count, TypeMap<Complex<Real>>(), opC, root, comm.comm ) );
    }
#else
    MPI_Op opC = ( op==SUM ? SumOp<Complex<Real>>().op : op.op );
    SafeMpi
    ( MPI_Reduce
      ( const_cast<Complex<Real>*>(sbuf), 
        rbuf, count, TypeMap<Complex<Real>>(), opC, root, comm.comm ) );
#endif
}

template<typename T,typename,typename>
void Reduce
( const T* sbuf, T* rbuf, int count, Op op, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    if( count == 0 )
        return;

    MPI_Op opC;
    if( op == SUM )
        opC = SumOp<T>().op; 
    else if( op == MAX )
        opC = MaxOp<T>().op;
    else if( op == MIN )
        opC = MinOp<T>().op;
    else
        opC = op.op;

    const int commRank = mpi::Rank(comm);
    std::vector<byte> packedSend, packedRecv;
    Serialize( count, sbuf, packedSend );

    if( commRank == root )
        ReserveSerialized( count, rbuf, packedRecv );
    SafeMpi
    ( MPI_Reduce
      ( packedSend.data(), packedRecv.data(), count, TypeMap<T>(),
        opC, root, comm.comm ) );
    if( commRank == root )
        Deserialize( count, packedRecv, rbuf );
}

template<typename T>
void Reduce( const T* sbuf, T* rbuf, int count, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{ Reduce( sbuf, rbuf, count, SUM, root, comm ); }

template<typename T>
T Reduce( T sb, Op op, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{ 
    T rb;
    Reduce( &sb, &rb, 1, op, root, comm );
    return rb;
}

template<typename T>
T Reduce( T sb, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{ 
    T rb;
    Reduce( &sb, &rb, 1, SUM, root, comm );
    return rb;
}

template<typename Real,typename>
void Reduce( Real* buf, int count, Op op, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    if( count == 0 || Size(comm) == 1 )
        return;

    MPI_Op opC;
    if( op == SUM )
        opC = SumOp<Real>().op; 
    else if( op == MAX )
        opC = MaxOp<Real>().op;
    else if( op == MIN )
        opC = MinOp<Real>().op;
    else
        opC = op.op;

    const int commRank = Rank( comm );
    if( commRank == root )
    {
        SafeMpi
        ( MPI_Reduce
          ( MPI_IN_PLACE, buf, count, TypeMap<Real>(), opC, root, 
            comm.comm ) );
    }
    else
        SafeMpi
        ( MPI_Reduce
          ( buf, 0, count, TypeMap<Real>(), opC, root, comm.comm ) );
}

template<typename Real,typename>
void Reduce( Complex<Real>* buf, int count, Op op, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    if( Size(comm) == 1 )
        return;
    if( count != 0 )
    {
        const int commRank = Rank( comm );
#ifdef EL_AVOID_COMPLEX_MPI
        if( op == SUM )
        {
            MPI_Op opC = SumOp<Real>().op;
            if( commRank == root )
            {
                SafeMpi
                ( MPI_Reduce
                  ( MPI_IN_PLACE, buf, 2*count, TypeMap<Real>(), opC, 
                    root, comm.comm ) );
            }
            else
                SafeMpi
                ( MPI_Reduce
                  ( buf, 0, 2*count, TypeMap<Real>(), opC, root, comm.comm ) );
        }
        else
        {
            MPI_Op opC = op.op;
            if( commRank == root )
            {
                SafeMpi
                ( MPI_Reduce
                  ( MPI_IN_PLACE, buf, count, TypeMap<Complex<Real>>(), opC, 
                    root, comm.comm ) );
            }
            else
                SafeMpi
                ( MPI_Reduce
                  ( buf, 0, count, TypeMap<Complex<Real>>(), opC, 
                    root, comm.comm ) );
        }
#else
        MPI_Op opC = ( op==SUM ? SumOp<Complex<Real>>().op : op.op );
        if( commRank == root )
        {
            SafeMpi
            ( MPI_Reduce
              ( MPI_IN_PLACE, buf, count, TypeMap<Complex<Real>>(), opC, 
                root, comm.comm ) );
        }
        else
            SafeMpi
            ( MPI_Reduce
              ( buf, 0, count, TypeMap<Complex<Real>>(), opC, root, 
                comm.comm ) );
#endif
    }
}

template<typename T,typename,typename>
void Reduce( T* buf, int count, Op op, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    if( count == 0 )
        return;

    MPI_Op opC;
    if( op == SUM )
        opC = SumOp<T>().op; 
    else if( op == MAX )
        opC = MaxOp<T>().op;
    else if( op == MIN )
        opC = MinOp<T>().op;
    else
        opC = op.op;

    // TODO: Use in-place option?

    const int commRank = mpi::Rank(comm);
    std::vector<byte> packedSend, packedRecv;
    Serialize( count, buf, packedSend );

    if( commRank == root )
        ReserveSerialized( count, buf, packedRecv );
    SafeMpi
    ( MPI_Reduce
      ( packedSend.data(), packedRecv.data(), count, TypeMap<T>(),
        opC, root, comm.comm ) );
    if( commRank == root )
        Deserialize( count, packedRecv, buf );
}

template<typename T>
void Reduce( T* buf, int count, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{ Reduce( buf, count, SUM, root, comm ); }

template<typename Real,typename>
void AllReduce( const Real* sbuf, Real* rbuf, int count, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    if( count != 0 )
    {
        MPI_Op opC;
        if( op == SUM )
            opC = SumOp<Real>().op; 
        else if( op == MAX )
            opC = MaxOp<Real>().op;
        else if( op == MIN )
            opC = MinOp<Real>().op;
        else
            opC = op.op;

        SafeMpi
        ( MPI_Allreduce
          ( const_cast<Real*>(sbuf), rbuf, count, TypeMap<Real>(), opC, 
            comm.comm ) );
    }
}

template<typename Real,typename>
void AllReduce
( const Complex<Real>* sbuf, Complex<Real>* rbuf, int count, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    if( count != 0 )
    {
#ifdef EL_AVOID_COMPLEX_MPI
        if( op == SUM )
        {
            MPI_Op opC = SumOp<Real>().op; 
            SafeMpi
            ( MPI_Allreduce
                ( const_cast<Complex<Real>*>(sbuf),
                  rbuf, 2*count, TypeMap<Real>(), opC, comm.comm ) );
        }
        else
        {
            MPI_Op opC = op.op;
            SafeMpi
            ( MPI_Allreduce
              ( const_cast<Complex<Real>*>(sbuf),
                rbuf, count, TypeMap<Complex<Real>>(), opC, comm.comm ) );
        }
#else
        MPI_Op opC = ( op==SUM ? SumOp<Complex<Real>>().op : op.op );
        SafeMpi
        ( MPI_Allreduce
          ( const_cast<Complex<Real>*>(sbuf), 
            rbuf, count, TypeMap<Complex<Real>>(), opC, comm.comm ) );
#endif
    }
}

template<typename T,typename,typename>
void AllReduce
( const T* sbuf, T* rbuf, int count, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    if( count == 0 )
        return;

    MPI_Op opC;
    if( op == SUM )
        opC = SumOp<T>().op; 
    else if( op == MAX )
        opC = MaxOp<T>().op;
    else if( op == MIN )
        opC = MinOp<T>().op;
    else
        opC = op.op;

    std::vector<byte> packedSend, packedRecv;
    Serialize( count, sbuf, packedSend );

    ReserveSerialized( count, rbuf, packedRecv );
    SafeMpi
    ( MPI_Allreduce
      ( packedSend.data(), packedRecv.data(), count, TypeMap<T>(),
        opC, comm.comm ) );
    Deserialize( count, packedRecv, rbuf );
}

template<typename T>
void AllReduce( const T* sbuf, T* rbuf, int count, Comm comm )
EL_NO_RELEASE_EXCEPT
{ AllReduce( sbuf, rbuf, count, SUM, comm ); }

template<typename T>
T AllReduce( T sb, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{ T rb; AllReduce( &sb, &rb, 1, op, comm ); return rb; }

template<typename T>
T AllReduce( T sb, Comm comm )
EL_NO_RELEASE_EXCEPT
{ return AllReduce( sb, SUM, comm ); }

template<typename Real,typename>
void AllReduce( Real* buf, int count, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    if( count == 0 || Size(comm) == 1 )
        return;

    MPI_Op opC;
    if( op == SUM )
        opC = SumOp<Real>().op; 
    else if( op == MAX )
        opC = MaxOp<Real>().op;
    else if( op == MIN )
        opC = MinOp<Real>().op;
    else
        opC = op.op;

    SafeMpi
    ( MPI_Allreduce
      ( MPI_IN_PLACE, buf, count, TypeMap<Real>(), opC, comm.comm ) );
}

template<typename Real,typename>
void AllReduce( Complex<Real>* buf, int count, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    if( count == 0 || Size(comm) == 1 )
        return;

#ifdef EL_AVOID_COMPLEX_MPI
    if( op == SUM )
    {
        MPI_Op opC = SumOp<Real>().op;
        SafeMpi
        ( MPI_Allreduce
          ( MPI_IN_PLACE, buf, 2*count, TypeMap<Real>(), opC, comm.comm ) );
    }
    else
    {
        MPI_Op opC = op.op;
        SafeMpi
        ( MPI_Allreduce
          ( MPI_IN_PLACE, buf, count, TypeMap<Complex<Real>>(), 
            opC, comm.comm ) );
    }
#else
    MPI_Op opC = ( op==SUM ? SumOp<Complex<Real>>().op : op.op );
    SafeMpi
    ( MPI_Allreduce
      ( MPI_IN_PLACE, buf, count, TypeMap<Complex<Real>>(), opC, 
        comm.comm ) );
#endif
}

template<typename T,typename,typename>
void AllReduce( T* buf, int count, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    if( count == 0 )
        return;

    MPI_Op opC;
    if( op == SUM )
        opC = SumOp<T>().op; 
    else if( op == MAX )
        opC = MaxOp<T>().op;
    else if( op == MIN )
        opC = MinOp<T>().op;
    else
        opC = op.op;

    std::vector<byte> packedSend, packedRecv;
    Serialize( count, buf, packedSend );

    ReserveSerialized( count, buf, packedRecv );
    SafeMpi
    ( MPI_Allreduce
      ( packedSend.data(), packedRecv.data(), count, TypeMap<T>(),
        opC, comm.comm ) );
    Deserialize( count, packedRecv, buf );
}

template<typename T>
void AllReduce( T* buf, int count, Comm comm )
EL_NO_RELEASE_EXCEPT
{ AllReduce( buf, count, SUM, comm ); }

template<typename Real,typename>
void ReduceScatter( Real* sbuf, Real* rbuf, int rc, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    if( rc == 0 )
        return;
#ifdef EL_REDUCE_SCATTER_BLOCK_VIA_ALLREDUCE
    const int commSize = Size( comm );
    const int commRank = Rank( comm );
    AllReduce( sbuf, rc*commSize, op, comm );
    MemCopy( rbuf, &sbuf[commRank*rc], rc );
#elif defined(EL_HAVE_MPI_REDUCE_SCATTER_BLOCK)
    MPI_Op opC;
    if( op == SUM )
        opC = SumOp<Real>().op; 
    else if( op == MAX )
        opC = MaxOp<Real>().op;
    else if( op == MIN )
        opC = MinOp<Real>().op;
    else
        opC = op.op;
    SafeMpi
    ( MPI_Reduce_scatter_block
      ( sbuf, rbuf, rc, TypeMap<Real>(), opC, comm.comm ) );
#else
    const int commSize = Size( comm );
    Reduce( sbuf, rc*commSize, op, 0, comm );
    Scatter( sbuf, rc, rbuf, rc, 0, comm );
#endif
}

template<typename Real,typename>
void ReduceScatter
( Complex<Real>* sbuf, Complex<Real>* rbuf, int rc, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    if( rc == 0 )
        return;

#ifdef EL_REDUCE_SCATTER_BLOCK_VIA_ALLREDUCE
    const int commSize = Size( comm );
    const int commRank = Rank( comm );
    AllReduce( sbuf, rc*commSize, op, comm );
    MemCopy( rbuf, &sbuf[commRank*rc], rc );
#elif defined(EL_HAVE_MPI_REDUCE_SCATTER_BLOCK)
# ifdef EL_AVOID_COMPLEX_MPI
    MPI_Op opC = ( op==SUM ? SumOp<Real>().op : op.op );
    SafeMpi
    ( MPI_Reduce_scatter_block
      ( sbuf, rbuf, 2*rc, TypeMap<Real>(), opC, comm.comm ) );
# else
    MPI_Op opC = ( op==SUM ? SumOp<Complex<Real>>().op : op.op );
    SafeMpi
    ( MPI_Reduce_scatter_block
      ( sbuf, rbuf, rc, TypeMap<Complex<Real>>(), opC, comm.comm ) );
# endif
#else
    const int commSize = Size( comm );
    Reduce( sbuf, rc*commSize, op, 0, comm );
    Scatter( sbuf, rc, rbuf, rc, 0, comm );
#endif
}

template<typename T,typename,typename>
void ReduceScatter( T* sbuf, T* rbuf, int rc, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    if( rc == 0 )
        return;
    const int commSize = mpi::Size(comm);
    const int totalSend = rc*commSize;
    const int totalRecv = rc;

    // TODO: Add AllReduce approach via EL_REDUCE_SCATTER_BLOCK_VIA_ALLREDUCE
#if defined(EL_HAVE_MPI_REDUCE_SCATTER_BLOCK)
    MPI_Op opC;
    if( op == SUM )
        opC = SumOp<T>().op; 
    else if( op == MAX )
        opC = MaxOp<T>().op;
    else if( op == MIN )
        opC = MinOp<T>().op;
    else
        opC = op.op;

    std::vector<byte> packedSend, packedRecv;
    Serialize( totalSend, sbuf, packedSend );

    ReserveSerialized( totalRecv, rbuf, packedRecv );
    SafeMpi
    ( MPI_Reduce_scatter_block
      ( packedSend.data(), packedRecv.data(), rc, TypeMap<T>(),
        opC, comm.comm ) );

    Deserialize( totalRecv, packedRecv, rbuf );
#else
    Reduce( sbuf, totalSend, op, 0, comm );
    Scatter( sbuf, rc, rbuf, rc, 0, comm );
#endif
}

template<typename T>
void ReduceScatter( T* sbuf, T* rbuf, int rc, Comm comm )
EL_NO_RELEASE_EXCEPT
{ ReduceScatter( sbuf, rbuf, rc, SUM, comm ); }

template<typename T>
T ReduceScatter( T sb, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{ T rb; ReduceScatter( &sb, &rb, 1, op, comm ); return rb; }

template<typename T>
T ReduceScatter( T sb, Comm comm )
EL_NO_RELEASE_EXCEPT
{ return ReduceScatter( sb, SUM, comm ); }

template<typename Real,typename>
void ReduceScatter( Real* buf, int rc, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    if( rc == 0 || Size(comm) == 1 )
        return;

#ifdef EL_REDUCE_SCATTER_BLOCK_VIA_ALLREDUCE
    const int commSize = Size( comm );
    const int commRank = Rank( comm );
    AllReduce( buf, rc*commSize, op, comm );
    if( commRank != 0 )
        MemCopy( buf, &buf[commRank*rc], rc );
#elif defined(EL_HAVE_MPI_REDUCE_SCATTER_BLOCK)
    MPI_Op opC;
    if( op == SUM )
        opC = SumOp<Real>().op; 
    else if( op == MAX )
        opC = MaxOp<Real>().op;
    else if( op == MIN )
        opC = MinOp<Real>().op;
    else
        opC = op.op;
    SafeMpi
    ( MPI_Reduce_scatter_block
      ( MPI_IN_PLACE, buf, rc, TypeMap<Real>(), opC, comm.comm ) );
#else
    const int commSize = Size( comm );
    Reduce( buf, rc*commSize, op, 0, comm );
    Scatter( buf, rc, rc, 0, comm );
#endif
}

// TODO: Handle case where op is not summation
template<typename Real,typename>
void ReduceScatter( Complex<Real>* buf, int rc, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    if( rc == 0 || Size(comm) == 1 )
        return;

#ifdef EL_REDUCE_SCATTER_BLOCK_VIA_ALLREDUCE
    const int commSize = Size( comm );
    const int commRank = Rank( comm );
    AllReduce( buf, rc*commSize, op, comm );
    if( commRank != 0 )
        MemCopy( buf, &buf[commRank*rc], rc );
#elif defined(EL_HAVE_MPI_REDUCE_SCATTER_BLOCK)
# ifdef EL_AVOID_COMPLEX_MPI
    MPI_Op opC = ( op==SUM ? SumOp<Real>().op : op.op );
    SafeMpi
    ( MPI_Reduce_scatter_block
      ( MPI_IN_PLACE, buf, 2*rc, TypeMap<Real>(), opC, comm.comm ) );
# else
    MPI_Op opC = ( op==SUM ? SumOp<Complex<Real>>().op : op.op );
    SafeMpi
    ( MPI_Reduce_scatter_block
      ( MPI_IN_PLACE, buf, rc, TypeMap<Complex<Real>>(), opC, comm.comm ) );
# endif
#else
    const int commSize = Size( comm );
    Reduce( buf, rc*commSize, op, 0, comm );
    Scatter( buf, rc, rc, 0, comm );
#endif
}

template<typename T,typename,typename>
void ReduceScatter( T* buf, int rc, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    if( rc == 0 )
        return;
    const int commSize = mpi::Size(comm);
    const int totalSend = rc*commSize;
    const int totalRecv = rc;

    // TODO: Add AllReduce approach via EL_REDUCE_SCATTER_BLOCK_VIA_ALLREDUCE
#if defined(EL_HAVE_MPI_REDUCE_SCATTER_BLOCK)
    MPI_Op opC;
    if( op == SUM )
        opC = SumOp<T>().op; 
    else if( op == MAX )
        opC = MaxOp<T>().op;
    else if( op == MIN )
        opC = MinOp<T>().op;
    else
        opC = op.op;

    std::vector<byte> packedSend, packedRecv;
    Serialize( totalSend, buf, packedSend );

    ReserveSerialized( totalRecv, buf, packedRecv );
    SafeMpi
    ( MPI_Reduce_scatter_block
      ( packedSend.data(), packedRecv.data(), rc, TypeMap<T>(),
        opC, comm.comm ) );

    Deserialize( totalRecv, packedRecv, buf );
#else
    Reduce( buf, totalSend, op, 0, comm );
    Scatter( buf, rc, rc, 0, comm );
#endif
}

template<typename T>
void ReduceScatter( T* buf, int rc, Comm comm )
EL_NO_RELEASE_EXCEPT
{ ReduceScatter( buf, rc, SUM, comm ); }

template<typename Real,typename>
void ReduceScatter
( const Real* sbuf, Real* rbuf, const int* rcs, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    MPI_Op opC;
    if( op == SUM )
        opC = SumOp<Real>().op; 
    else if( op == MAX )
        opC = MaxOp<Real>().op;
    else if( op == MIN )
        opC = MinOp<Real>().op;
    else
        opC = op.op;

    SafeMpi
    ( MPI_Reduce_scatter
      ( const_cast<Real*>(sbuf), 
        rbuf, const_cast<int*>(rcs), TypeMap<Real>(), opC, comm.comm ) );
}

template<typename Real,typename>
void ReduceScatter
( const Complex<Real>* sbuf, Complex<Real>* rbuf, const int* rcs,
  Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
#ifdef EL_AVOID_COMPLEX_MPI
    if( op == SUM )
    {
        MPI_Op opC = SumOp<Real>().op;
        int p;
        MPI_Comm_size( comm.comm, &p );
        vector<int> rcsDoubled(p);
        for( int i=0; i<p; ++i )
            rcsDoubled[i] = 2*rcs[i];
        SafeMpi
        ( MPI_Reduce_scatter
          ( const_cast<Complex<Real>*>(sbuf),
            rbuf, rcsDoubled.data(), TypeMap<Real>(), opC, comm.comm ) );
    }
    else
    {
        MPI_Op opC = op.op;
        SafeMpi
        ( MPI_Reduce_scatter
          ( const_cast<Complex<Real>*>(sbuf),
            rbuf, const_cast<int*>(rcs), TypeMap<Complex<Real>>(), 
            opC, comm.comm ) );
    }
#else
    MPI_Op opC = ( op==SUM ? SumOp<Complex<Real>>().op : op.op );
    SafeMpi
    ( MPI_Reduce_scatter
      ( const_cast<Complex<Real>*>(sbuf), 
        rbuf, const_cast<int*>(rcs), TypeMap<Complex<Real>>(), opC, 
        comm.comm ) );
#endif
}

template<typename T,typename,typename>
void ReduceScatter
( const T* sbuf, T* rbuf, const int* rcs, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    const int commRank = mpi::Rank(comm);
    const int commSize = mpi::Size(comm);
    int totalSend=0;
    for( int q=0; q<commSize; ++q )
        totalSend += rcs[q];
    const int totalRecv = rcs[commRank];

    MPI_Op opC;
    if( op == SUM )
        opC = SumOp<T>().op; 
    else if( op == MAX )
        opC = MaxOp<T>().op;
    else if( op == MIN )
        opC = MinOp<T>().op;
    else
        opC = op.op;

    std::vector<byte> packedSend, packedRecv;
    Serialize( totalSend, sbuf, packedSend );
    ReserveSerialized( totalRecv, rbuf, packedRecv );
    SafeMpi
    ( MPI_Reduce_scatter
      ( packedSend.data(), packedRecv.data(), const_cast<int*>(rcs),
        TypeMap<T>(), opC, comm.comm ) );
    Deserialize( totalRecv, packedRecv, rbuf );
}

template<typename T>
void ReduceScatter( const T* sbuf, T* rbuf, const int* rcs, Comm comm )
EL_NO_RELEASE_EXCEPT
{ ReduceScatter( sbuf, rbuf, rcs, SUM, comm ); }

void VerifySendsAndRecvs
( const vector<int>& sendCounts,
  const vector<int>& recvCounts, Comm comm )
{
    DEBUG_CSE
    const int commSize = Size( comm );
    vector<int> actualRecvCounts(commSize);
    AllToAll
    ( sendCounts.data(),       1,
      actualRecvCounts.data(), 1, comm );
    for( int q=0; q<commSize; ++q )
        if( actualRecvCounts[q] != recvCounts[q] )
            LogicError
            ("Expected recv count of ",recvCounts[q],
             " but recv'd ",actualRecvCounts[q]," from process ",q);
}

template<typename Real,typename>
void Scan( const Real* sbuf, Real* rbuf, int count, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    if( count != 0 )
    {
        MPI_Op opC;
        if( op == SUM )
            opC = SumOp<Real>().op; 
        else if( op == MAX )
            opC = MaxOp<Real>().op;
        else if( op == MIN )
            opC = MinOp<Real>().op;
        else
            opC = op.op;

        SafeMpi
        ( MPI_Scan
          ( const_cast<Real*>(sbuf), rbuf, count, TypeMap<Real>(),
            opC, comm.comm ) );
    }
}

template<typename Real,typename>
void Scan
( const Complex<Real>* sbuf, 
        Complex<Real>* rbuf, int count, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    if( count != 0 )
    {
#ifdef EL_AVOID_COMPLEX_MPI
        if( op == SUM )
        {
            MPI_Op opC = SumOp<Real>().op;
            SafeMpi
            ( MPI_Scan
              ( const_cast<Complex<Real>*>(sbuf),
                rbuf, 2*count, TypeMap<Real>(), opC, comm.comm ) );
        }
        else
        {
            MPI_Op opC = op.op;
            SafeMpi
            ( MPI_Scan
              ( const_cast<Complex<Real>*>(sbuf),
                rbuf, count, TypeMap<Complex<Real>>(), opC, comm.comm ) );
        }
#else
        MPI_Op opC = ( op==SUM ? SumOp<Complex<Real>>().op : op.op );
        SafeMpi
        ( MPI_Scan
          ( const_cast<Complex<Real>*>(sbuf), 
            rbuf, count, TypeMap<Complex<Real>>(), opC, comm.comm ) );
#endif
    }
}

template<typename T,typename,typename>
void Scan( const T* sbuf, T* rbuf, int count, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    if( count == 0 )
        return;

    MPI_Op opC;
    if( op == SUM )
        opC = SumOp<T>().op; 
    else if( op == MAX )
        opC = MaxOp<T>().op;
    else if( op == MIN )
        opC = MinOp<T>().op;
    else
        opC = op.op;

    std::vector<byte> packedSend, packedRecv;
    Serialize( count, sbuf, packedSend );
    ReserveSerialized( count, rbuf, packedRecv );
    SafeMpi
    ( MPI_Scan
      ( packedSend.data(), packedRecv.data(), count, TypeMap<T>(),
        opC, comm.comm ) );
    Deserialize( count, packedRecv, rbuf );
}

template<typename T>
void Scan( const T* sbuf, T* rbuf, int count, Comm comm )
EL_NO_RELEASE_EXCEPT
{ Scan( sbuf, rbuf, count, SUM, comm ); }

template<typename T>
T Scan( T sb, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{ 
    T rb;
    Scan( &sb, &rb, 1, op, comm );
    return rb;
}

template<typename T>
T Scan( T sb, Comm comm )
EL_NO_RELEASE_EXCEPT
{ 
    T rb;
    Scan( &sb, &rb, 1, SUM, comm );
    return rb;
}

template<typename Real,typename>
void Scan( Real* buf, int count, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    if( count != 0 )
    {
        MPI_Op opC;
        if( op == SUM )
            opC = SumOp<Real>().op; 
        else if( op == MAX )
            opC = MaxOp<Real>().op;
        else if( op == MIN )
            opC = MinOp<Real>().op;
        else
            opC = op.op;

        SafeMpi
        ( MPI_Scan
          ( MPI_IN_PLACE, buf, count, TypeMap<Real>(), opC, comm.comm ) );
    }
}

template<typename Real,typename>
void Scan( Complex<Real>* buf, int count, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    if( count != 0 )
    {
#ifdef EL_AVOID_COMPLEX_MPI
        if( op == SUM )
        {
            MPI_Op opC = SumOp<Real>().op;
            SafeMpi
            ( MPI_Scan
              ( MPI_IN_PLACE, buf, 2*count, TypeMap<Real>(), opC, comm.comm ) );
        }
        else
        {
            MPI_Op opC = op.op;
            SafeMpi
            ( MPI_Scan
              ( MPI_IN_PLACE, buf, count, TypeMap<Complex<Real>>(), opC, 
                comm.comm ) );
        }
#else
        MPI_Op opC = ( op==SUM ? SumOp<Complex<Real>>().op : op.op );
        SafeMpi
        ( MPI_Scan
          ( MPI_IN_PLACE, buf, count, TypeMap<Complex<Real>>(), opC, 
            comm.comm ) );
#endif
    }
}

template<typename T,typename,typename>
void Scan( T* buf, int count, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    if( count == 0 )
        return;

    MPI_Op opC;
    if( op == SUM )
        opC = SumOp<T>().op; 
    else if( op == MAX )
        opC = MaxOp<T>().op;
    else if( op == MIN )
        opC = MinOp<T>().op;
    else
        opC = op.op;

    std::vector<byte> packedSend, packedRecv;
    Serialize( count, buf, packedSend );
    ReserveSerialized( count, buf, packedRecv );
    SafeMpi
    ( MPI_Scan
      ( packedSend.data(), packedRecv.data(), count, TypeMap<T>(),
        opC, comm.comm ) );
    Deserialize( count, packedRecv, buf );
}

template<typename T>
void Scan( T* buf, int count, Comm comm )
EL_NO_RELEASE_EXCEPT
{ Scan( buf, count, SUM, comm ); }

template<typename T>
void SparseAllToAll
( const vector<T>& sendBuffer,
  const vector<int>& sendCounts,
  const vector<int>& sendDispls,
        vector<T>& recvBuffer,
  const vector<int>& recvCounts,
  const vector<int>& recvDispls,
        Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_CSE
    DEBUG_ONLY(VerifySendsAndRecvs( sendCounts, recvCounts, comm ))
#ifdef EL_USE_CUSTOM_ALLTOALLV
    const int commSize = Size( comm );
    int numSends=0,numRecvs=0;
    for( int q=0; q<commSize; ++q )
    {
        if( sendCounts[q] != 0 )
            ++numSends;
        if( recvCounts[q] != 0 )
            ++numRecvs;
    }
    vector<Status> statuses(numSends+numRecvs);
    vector<Request<T>> requests(numSends+numRecvs);
    int rCount=0;
    for( int q=0; q<commSize; ++q )
    {
        int count = recvCounts[q];
        int displ = recvDispls[q];
        if( count != 0 )
            IRecv( &recvBuffer[displ], count, q, comm, requests[rCount++] );
    }

    // Ensure that recvs are posted before the sends
    // (Invalid MPI_Irecv's have been observed otherwise)
    Barrier( comm );

    for( int q=0; q<commSize; ++q )
    {
        int count = sendCounts[q];
        int displ = sendDispls[q];
        if( count != 0 )
            IRSend( &sendBuffer[displ], count, q, comm, requests[rCount++] );
    }
    WaitAll( numSends+numRecvs, requests.data(), statuses.data() );
#else
    AllToAll
    ( sendBuffer.data(), sendCounts.data(), sendDispls.data(),
      recvBuffer.data(), recvCounts.data(), recvDispls.data(), comm );
#endif
}

#define MPI_PROTO(T) \
  template bool Test( Request<T>& request ) EL_NO_RELEASE_EXCEPT; \
  template void Wait( Request<T>& request ) EL_NO_RELEASE_EXCEPT; \
  template void Wait( Request<T>& request, Status& status ) \
  EL_NO_RELEASE_EXCEPT; \
  template void WaitAll( int numRequests, Request<T>* requests ) \
  EL_NO_RELEASE_EXCEPT; \
  template void WaitAll \
  ( int numRequests, Request<T>* requests, Status* statuses ) \
  EL_NO_RELEASE_EXCEPT; \
  template int GetCount<T>( Status& status ) EL_NO_RELEASE_EXCEPT; \
  template void TaggedSend \
  ( const T* buf, int count, int to, int tag, Comm comm ) \
  EL_NO_RELEASE_EXCEPT; \
  template void Send( const T* buf, int count, int to, Comm comm ) \
  EL_NO_RELEASE_EXCEPT; \
  template void TaggedSend( T b, int to, int tag, Comm comm ) \
  EL_NO_RELEASE_EXCEPT; \
  template void Send( T b, int to, Comm comm ) \
  EL_NO_RELEASE_EXCEPT; \
  template void TaggedISend \
  ( const T* buf, int count, int to, int tag, Comm comm, Request<T>& request ) \
  EL_NO_RELEASE_EXCEPT; \
  template void ISend \
  ( const T* buf, int count, int to, Comm comm, Request<T>& request ) \
  EL_NO_RELEASE_EXCEPT; \
  template void TaggedISend \
  ( T buf, int to, int tag, Comm comm, Request<T>& request ) \
  EL_NO_RELEASE_EXCEPT; \
  template void ISend( T buf, int to, Comm comm, Request<T>& request ) \
  EL_NO_RELEASE_EXCEPT; \
  template void TaggedISSend \
  ( const T* buf, int count, int to, int tag, Comm comm, Request<T>& request ) \
  EL_NO_RELEASE_EXCEPT; \
  template void ISSend \
  ( const T* buf, int count, int to, Comm comm, Request<T>& request ) \
  EL_NO_RELEASE_EXCEPT; \
  template void TaggedISSend \
  ( T b, int to, int tag, Comm comm, Request<T>& request ) \
  EL_NO_RELEASE_EXCEPT; \
  template void TaggedRecv \
  ( T* buf, int count, int from, int tag, Comm comm ) EL_NO_RELEASE_EXCEPT; \
  template void Recv( T* buf, int count, int from, Comm comm ) \
  EL_NO_RELEASE_EXCEPT; \
  template T TaggedRecv<T>( int from, int tag, Comm comm ) \
  EL_NO_RELEASE_EXCEPT; \
  template T Recv( int from, Comm comm ) EL_NO_RELEASE_EXCEPT; \
  template void TaggedIRecv \
  ( T* buf, int count, int from, int tag, Comm comm, Request<T>& request ) \
  EL_NO_RELEASE_EXCEPT; \
  template void IRecv \
  ( T* buf, int count, int from, Comm comm, Request<T>& request ) \
  EL_NO_RELEASE_EXCEPT; \
  template T TaggedIRecv<T> \
  ( int from, int tag, Comm comm, Request<T>& request ) EL_NO_RELEASE_EXCEPT; \
  template T IRecv<T>( int from, Comm comm, Request<T>& request ) \
  EL_NO_RELEASE_EXCEPT; \
  template void TaggedSendRecv \
  ( const T* sbuf, int sc, int to,   int stag, \
          T* rbuf, int rc, int from, int rtag, Comm comm ) \
  EL_NO_RELEASE_EXCEPT; \
  template void SendRecv \
  ( const T* sbuf, int sc, int to, \
          T* rbuf, int rc, int from, Comm comm ) EL_NO_RELEASE_EXCEPT; \
  template T TaggedSendRecv \
  ( T sb, int to, int stag, int from, int rtag, Comm comm ); \
  template T SendRecv( T sb, int to, int from, Comm comm ) \
  EL_NO_RELEASE_EXCEPT; \
  template void TaggedSendRecv \
  ( T* buf, int count, int to, int stag, int from, int rtag, Comm comm ) \
  EL_NO_RELEASE_EXCEPT; \
  template void SendRecv \
  ( T* buf, int count, int to, int from, Comm comm ) \
  EL_NO_RELEASE_EXCEPT; \
  template void Broadcast( T* buf, int count, int root, Comm comm ) \
  EL_NO_RELEASE_EXCEPT; \
  template void Broadcast( T& b, int root, Comm comm ) \
  EL_NO_RELEASE_EXCEPT; \
  template void IBroadcast \
  ( T* buf, int count, int root, Comm comm, Request<T>& request ); \
  template void IBroadcast \
  ( T& b, int root, Comm comm, Request<T>& request ); \
  template void Gather \
  ( const T* sbuf, int sc, T* rbuf, int rc, int root, Comm comm ) \
  EL_NO_RELEASE_EXCEPT; \
  template void IGather \
  ( const T* sbuf, int sc, \
          T* rbuf, int rc, \
    int root, Comm comm, Request<T>& request ); \
  template void Gather \
  ( const T* sbuf, int sc, \
          T* rbuf, const int* rcs, const int* rds, int root, Comm comm ) \
  EL_NO_RELEASE_EXCEPT; \
  template void AllGather( const T* sbuf, int sc, T* rbuf, int rc, Comm comm ) \
  EL_NO_RELEASE_EXCEPT; \
  template void AllGather \
  ( const T* sbuf, int sc, \
          T* rbuf, const int* rcs, const int* rds, Comm comm ) \
  EL_NO_RELEASE_EXCEPT; \
  template void Scatter \
  ( const T* sbuf, int sc, \
          T* rbuf, int rc, int root, Comm comm ) \
  EL_NO_RELEASE_EXCEPT; \
  template void Scatter( T* buf, int sc, int rc, int root, Comm comm ) \
  EL_NO_RELEASE_EXCEPT; \
  template void AllToAll \
  ( const T* sbuf, int sc, \
          T* rbuf, int rc, Comm comm ) \
  EL_NO_RELEASE_EXCEPT; \
  template void AllToAll \
  ( const T* sbuf, const int* scs, const int* sds, \
          T* rbuf, const int* rcs, const int* rds, Comm comm ) \
  EL_NO_RELEASE_EXCEPT; \
  template vector<T> AllToAll \
  ( const vector<T>& sendBuf, \
    const vector<int>& sendCounts, \
    const vector<int>& sendOffs, \
    Comm comm ) \
  EL_NO_RELEASE_EXCEPT; \
  template void Reduce \
  ( const T* sbuf, T* rbuf, int count, Op op, int root, Comm comm ) \
  EL_NO_RELEASE_EXCEPT; \
  template void Reduce \
  ( const T* sbuf, T* rbuf, int count, int root, Comm comm ) \
  EL_NO_RELEASE_EXCEPT; \
  template T Reduce( T sb, Op op, int root, Comm comm ) \
  EL_NO_RELEASE_EXCEPT; \
  template T Reduce( T sb, int root, Comm comm ) \
  EL_NO_RELEASE_EXCEPT; \
  template void Reduce( T* buf, int count, Op op, int root, Comm comm ) \
  EL_NO_RELEASE_EXCEPT; \
  template void Reduce( T* buf, int count, int root, Comm comm ) \
  EL_NO_RELEASE_EXCEPT; \
  template void AllReduce \
  ( const T* sbuf, T* rbuf, int count, Op op, Comm comm ) \
  EL_NO_RELEASE_EXCEPT; \
  template void AllReduce( const T* sbuf, T* rbuf, int count, Comm comm ) \
  EL_NO_RELEASE_EXCEPT; \
  template T AllReduce( T sb, Op op, Comm comm ) \
  EL_NO_RELEASE_EXCEPT; \
  template T AllReduce( T sb, Comm comm ) \
  EL_NO_RELEASE_EXCEPT; \
  template void AllReduce( T* buf, int count, Op op, Comm comm ) \
  EL_NO_RELEASE_EXCEPT; \
  template void AllReduce( T* buf, int count, Comm comm ) \
  EL_NO_RELEASE_EXCEPT; \
  template void ReduceScatter( T* sbuf, T* rbuf, int rc, Op op, Comm comm ) \
  EL_NO_RELEASE_EXCEPT; \
  template void ReduceScatter( T* sbuf, T* rbuf, int rc, Comm comm ) \
  EL_NO_RELEASE_EXCEPT; \
  template T ReduceScatter( T sb, Op op, Comm comm ) \
  EL_NO_RELEASE_EXCEPT; \
  template T ReduceScatter( T sb, Comm comm ) \
  EL_NO_RELEASE_EXCEPT; \
  template void ReduceScatter( T* buf, int rc, Op op, Comm comm ) \
  EL_NO_RELEASE_EXCEPT; \
  template void ReduceScatter( T* buf, int rc, Comm comm ) \
  EL_NO_RELEASE_EXCEPT; \
  template void ReduceScatter \
  ( const T* sbuf, T* rbuf, const int* rcs, Op op, Comm comm ) \
  EL_NO_RELEASE_EXCEPT; \
  template void ReduceScatter \
  ( const T* sbuf, T* rbuf, const int* rcs, Comm comm ) \
  EL_NO_RELEASE_EXCEPT; \
  template void Scan( const T* sbuf, T* rbuf, int count, Op op, Comm comm ) \
  EL_NO_RELEASE_EXCEPT; \
  template void Scan( const T* sbuf, T* rbuf, int count, Comm comm ) \
  EL_NO_RELEASE_EXCEPT; \
  template T Scan( T sb, Op op, Comm comm ) \
  EL_NO_RELEASE_EXCEPT; \
  template T Scan( T sb, Comm comm ) \
  EL_NO_RELEASE_EXCEPT; \
  template void Scan( T* buf, int count, Op op, Comm comm ) \
  EL_NO_RELEASE_EXCEPT; \
  template void Scan( T* buf, int count, Comm comm ) \
  EL_NO_RELEASE_EXCEPT;

MPI_PROTO(byte)
MPI_PROTO(int)
MPI_PROTO(unsigned)
MPI_PROTO(long int)
MPI_PROTO(unsigned long)
#ifdef EL_HAVE_MPI_LONG_LONG
MPI_PROTO(long long int)
MPI_PROTO(unsigned long long)
#endif
MPI_PROTO(ValueInt<Int>)
MPI_PROTO(Entry<Int>)
MPI_PROTO(float)
MPI_PROTO(Complex<float>)
MPI_PROTO(ValueInt<float>)
MPI_PROTO(ValueInt<Complex<float>>)
MPI_PROTO(Entry<float>)
MPI_PROTO(Entry<Complex<float>>)
MPI_PROTO(double)
MPI_PROTO(Complex<double>)
MPI_PROTO(ValueInt<double>)
MPI_PROTO(ValueInt<Complex<double>>)
MPI_PROTO(Entry<double>)
MPI_PROTO(Entry<Complex<double>>)
#ifdef EL_HAVE_QD
MPI_PROTO(DoubleDouble)
MPI_PROTO(QuadDouble)
MPI_PROTO(Complex<DoubleDouble>)
MPI_PROTO(Complex<QuadDouble>)
MPI_PROTO(ValueInt<DoubleDouble>)
MPI_PROTO(ValueInt<QuadDouble>)
MPI_PROTO(ValueInt<Complex<DoubleDouble>>)
MPI_PROTO(ValueInt<Complex<QuadDouble>>)
MPI_PROTO(Entry<DoubleDouble>)
MPI_PROTO(Entry<QuadDouble>)
MPI_PROTO(Entry<Complex<DoubleDouble>>)
MPI_PROTO(Entry<Complex<QuadDouble>>)
#endif
#ifdef EL_HAVE_QUAD
MPI_PROTO(Quad)
MPI_PROTO(Complex<Quad>)
MPI_PROTO(ValueInt<Quad>)
MPI_PROTO(ValueInt<Complex<Quad>>)
MPI_PROTO(Entry<Quad>)
MPI_PROTO(Entry<Complex<Quad>>)
#endif
#ifdef EL_HAVE_MPC
MPI_PROTO(BigInt)
MPI_PROTO(BigFloat)
MPI_PROTO(Complex<BigFloat>)
MPI_PROTO(ValueInt<BigInt>)
MPI_PROTO(ValueInt<BigFloat>)
MPI_PROTO(ValueInt<Complex<BigFloat>>)
MPI_PROTO(Entry<BigInt>)
MPI_PROTO(Entry<BigFloat>)
MPI_PROTO(Entry<Complex<BigFloat>>)
#endif

#define PROTO(T) \
  template void SparseAllToAll \
  ( const vector<T>& sendBuffer, \
    const vector<int>& sendCounts, \
    const vector<int>& sendDispls, \
          vector<T>& recvBuffer, \
    const vector<int>& recvCounts, \
    const vector<int>& recvDispls, \
          Comm comm ) EL_NO_RELEASE_EXCEPT;

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace mpi
} // namespace El
