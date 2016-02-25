/*
   Copyright (c) 2009-2016, Jack Poulson
                      2013, Jeff Hammond
                      2013, Jed Brown
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

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
    DEBUG_ONLY(CSE cse("mpi::Create"))
    SafeMpi( MPI_Op_create( func, commutes, &op.op ) );
}

void Free( Op& op ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Free"))
    SafeMpi( MPI_Op_free( &op.op ) );
}

void Free( Datatype& type ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Free"))
    SafeMpi( MPI_Type_free( &type ) );
}

// Communicator manipulation 
// =========================
int Rank( Comm comm ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Rank"))
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
    DEBUG_ONLY(CSE cse("mpi::Size"))
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
    DEBUG_ONLY(CSE cse("mpi::Create"))
    SafeMpi( 
        MPI_Comm_create( parentComm.comm, subsetGroup.group, &subsetComm.comm ) 
    );
}

void Dup( Comm original, Comm& duplicate ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Dup"))
    SafeMpi( MPI_Comm_dup( original.comm, &duplicate.comm ) );
}

void Split( Comm comm, int color, int key, Comm& newComm ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Split"))
    SafeMpi( MPI_Comm_split( comm.comm, color, key, &newComm.comm ) );
}

void Free( Comm& comm ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Free"))
    SafeMpi( MPI_Comm_free( &comm.comm ) );
}

bool Congruent( Comm comm1, Comm comm2 ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Congruent"))
    int result;
    SafeMpi( MPI_Comm_compare( comm1.comm, comm2.comm, &result ) );
    return ( result == MPI_IDENT || result == MPI_CONGRUENT );
}

void ErrorHandlerSet( Comm comm, ErrorHandler errorHandler )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::ErrorHandlerSet"))
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
    DEBUG_ONLY(CSE cse("mpi::CartCreate"))
    SafeMpi
    ( MPI_Cart_create
      ( comm.comm, numDims, const_cast<int*>(dimensions), 
        const_cast<int*>(periods), reorder, &cartComm.comm ) );
}

void CartSub( Comm comm, const int* remainingDims, Comm& subComm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::CartSub"))
    SafeMpi( 
      MPI_Cart_sub
      ( comm.comm, const_cast<int*>(remainingDims), &subComm.comm ) 
    );
}

// Group manipulation 
// ==================

int Rank( Group group ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Rank"))
    int rank;
    SafeMpi( MPI_Group_rank( group.group, &rank ) );
    return rank;
}

int Size( Group group ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Size"))
    int size;
    SafeMpi( MPI_Group_size( group.group, &size ) );
    return size;
}

void CommGroup( Comm comm, Group& group ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::CommGroup"))
    SafeMpi( MPI_Comm_group( comm.comm, &group.group ) );
}

void Dup( Group group, Group& newGroup ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Dup"))
    // For some reason, MPI_Group_dup does not exist
    Excl( group, 0, 0, newGroup ); 
}

void Union( Group groupA, Group groupB, Group& newGroup ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Union"))
    SafeMpi( MPI_Group_union( groupA.group, groupB.group, &newGroup.group ) );
}

void Incl( Group group, int n, const int* ranks, Group& subGroup )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Incl"))
    SafeMpi( 
      MPI_Group_incl
      ( group.group, n, const_cast<int*>(ranks), &subGroup.group ) 
    );
}

void Excl( Group group, int n, const int* ranks, Group& subGroup )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Excl"))
    SafeMpi(
      MPI_Group_excl
      ( group.group, n, const_cast<int*>(ranks), &subGroup.group )
    );
}

void Difference( Group parent, Group subset, Group& complement )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Difference"))
    SafeMpi( 
      MPI_Group_difference( parent.group, subset.group, &complement.group ) 
    );
}

void Free( Group& group ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Free"))
    SafeMpi( MPI_Group_free( &group.group ) );
}

bool Congruent( Group group1, Group group2 ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Congruent"))
    int result;
    SafeMpi( MPI_Group_compare( group1.group, group2.group, &result ) );
    return ( result == MPI_IDENT );
}

// Rank translations
// =================

int Translate( Group origGroup, int origRank, Group newGroup )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Translate"))
    int newRank;
    Translate( origGroup, 1, &origRank, newGroup, &newRank );
    return newRank;
}

int Translate( Comm origComm, int origRank, Group newGroup )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Translate"))
    int newRank;
    Translate( origComm, 1, &origRank, newGroup, &newRank );
    return newRank;
}

int Translate( Group origGroup, int origRank, Comm newComm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Translate"))
    int newRank;
    Translate( origGroup, 1, &origRank, newComm, &newRank );
    return newRank;
}

int Translate( Comm origComm, int origRank, Comm newComm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Translate"))
    int newRank;
    Translate( origComm, 1, &origRank, newComm, &newRank );
    return newRank;
}

void Translate
( Group origGroup, int size, const int* origRanks, 
  Group newGroup,                  int* newRanks ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Translate"))
    SafeMpi
    ( MPI_Group_translate_ranks
      ( origGroup.group, size, const_cast<int*>(origRanks), 
        newGroup.group, newRanks ) );
}

void Translate
( Comm origComm,  int size, const int* origRanks, 
  Group newGroup,                 int* newRanks ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Translate"))
    Group origGroup;
    CommGroup( origComm, origGroup );
    Translate( origGroup, size, origRanks, newGroup, newRanks );
    Free( origGroup );
}

void Translate
( Group origGroup,  int size, const int* origRanks, 
  Comm newComm,                     int* newRanks ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Translate"))
    Group newGroup;
    CommGroup( newComm,  newGroup  );
    Translate( origGroup, size, origRanks, newGroup, newRanks );
    Free( newGroup  );
}

void Translate
( Comm origComm,  int size, const int* origRanks, 
  Comm newComm,                   int* newRanks ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Translate"))
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
    DEBUG_ONLY(CSE cse("mpi::Barrier"))
    SafeMpi( MPI_Barrier( comm.comm ) );
}

// Test for completion
template<typename T>
bool Test( Request<T>& request ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Test"))
    Status status;
    int flag;
    SafeMpi( MPI_Test( &request.backend, &flag, &status ) );
    return flag;
}

// Ensure that the request finishes before continuing
template<typename T>
void Wait( Request<T>& request ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Wait"))
    Status status;
    Wait( request, status );
}

// Ensure that the request finishes before continuing
template<typename T>
void Wait( Request<T>& request, Status& status ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Wait"))
    SafeMpi( MPI_Wait( &request.backend, &status ) );
}

// Ensure that several requests finish before continuing
template<typename T>
void WaitAll( int numRequests, Request<T>* requests ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::WaitAll"))
    vector<Status> statuses( numRequests );
    WaitAll( numRequests, requests, statuses.data() );
}

// Ensure that several requests finish before continuing
template<typename T>
void WaitAll( int numRequests, Request<T>* requests, Status* statuses )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::WaitAll"))
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

#ifdef EL_HAVE_MPC
template<typename T>
void PackedWait( Request<T>& request, Status& status ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::PackedWait"))
    SafeMpi( MPI_Wait( &request.backend, &status ) );
    if( request.receivingPacked )
    {
        Deserialize
        ( request.recvCount, request.buffer.data(), request.unpackedRecvBuf );
        request.receivingPacked = false;
    }
    request.buffer.clear();
}

template<>
void Wait( Request<BigInt>& request, Status& status ) EL_NO_RELEASE_EXCEPT
{ PackedWait( request, status ); }
template<>
void Wait( Request<ValueInt<BigInt>>& request, Status& status )
EL_NO_RELEASE_EXCEPT
{ PackedWait( request, status ); }
template<>
void Wait( Request<Entry<BigInt>>& request, Status& status )
EL_NO_RELEASE_EXCEPT
{ PackedWait( request, status ); }

template<>
void Wait( Request<BigFloat>& request, Status& status ) EL_NO_RELEASE_EXCEPT
{ PackedWait( request, status ); }
template<>
void Wait( Request<ValueInt<BigFloat>>& request, Status& status )
EL_NO_RELEASE_EXCEPT
{ PackedWait( request, status ); }
template<>
void Wait( Request<Entry<BigFloat>>& request, Status& status )
EL_NO_RELEASE_EXCEPT
{ PackedWait( request, status ); }

template<typename T>
void PackedWaitAll( int numRequests, Request<T>* requests, Status* statuses )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::PackedWaitAll"))
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

template<>
void WaitAll
( int numRequests, Request<BigInt>* requests, Status* statuses )
{ PackedWaitAll( numRequests, requests, statuses ); }
template<>
void WaitAll
( int numRequests, Request<ValueInt<BigInt>>* requests, Status* statuses )
{ PackedWaitAll( numRequests, requests, statuses ); }
template<>
void WaitAll
( int numRequests, Request<Entry<BigInt>>* requests, Status* statuses )
{ PackedWaitAll( numRequests, requests, statuses ); }

template<>
void WaitAll
( int numRequests, Request<BigFloat>* requests, Status* statuses )
{ PackedWaitAll( numRequests, requests, statuses ); }
template<>
void WaitAll
( int numRequests, Request<ValueInt<BigFloat>>* requests, Status* statuses )
{ PackedWaitAll( numRequests, requests, statuses ); }
template<>
void WaitAll
( int numRequests, Request<Entry<BigFloat>>* requests, Status* statuses )
{ PackedWaitAll( numRequests, requests, statuses ); }
#endif

// Nonblocking test for message completion
bool IProbe( int source, int tag, Comm comm, Status& status )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::IProbe"))
    int flag;
    SafeMpi( MPI_Iprobe( source, tag, comm.comm, &flag, &status ) );
    return flag;
}
bool IProbe( int source, Comm comm, Status& status ) EL_NO_RELEASE_EXCEPT
{ return IProbe( source, 0, comm, status ); }

template<typename T>
int GetCount( Status& status ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::GetCount"))
    int count;
    SafeMpi( MPI_Get_count( &status, TypeMap<T>(), &count ) );
    return count;
}

template<typename Real>
void TaggedSend( const Real* buf, int count, int to, int tag, Comm comm )
EL_NO_RELEASE_EXCEPT
{ 
    DEBUG_ONLY(CSE cse("mpi::TaggedSend"))
    SafeMpi
    ( MPI_Send
      ( const_cast<Real*>(buf), count, TypeMap<Real>(), to, tag, comm.comm ) );
}
#ifdef EL_HAVE_MPC
template<typename T>
void PackedTaggedSend( const T* buf, int count, int to, int tag, Comm comm )
{
    DEBUG_ONLY(CSE cse("mpi::PackedTaggedSend"))
    std::vector<byte> packedBuf;
    Serialize( count, buf, packedBuf );
    SafeMpi
    ( MPI_Send( packedBuf.data(), count, TypeMap<T>(), to, tag, comm.comm ) );
}
template<>
void TaggedSend( const BigInt* buf, int count, int to, int tag, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::TaggedSend [BigInt]"))
    PackedTaggedSend( buf, count, to, tag, comm );
}
template<>
void TaggedSend
( const ValueInt<BigInt>* buf, int count, int to, int tag, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::TaggedSend [ValueInt<BigInt>]"))
    PackedTaggedSend( buf, count, to, tag, comm );
}
template<>
void TaggedSend
( const Entry<BigInt>* buf, int count, int to, int tag, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::TaggedSend [Entry<BigInt>]"))
    PackedTaggedSend( buf, count, to, tag, comm );
}

template<>
void TaggedSend( const BigFloat* buf, int count, int to, int tag, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::TaggedSend [BigFloat]"))
    PackedTaggedSend( buf, count, to, tag, comm );
}
template<>
void TaggedSend
( const ValueInt<BigFloat>* buf, int count, int to, int tag, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::TaggedSend [ValueInt<BigFloat>]"))
    PackedTaggedSend( buf, count, to, tag, comm );
}
template<>
void TaggedSend
( const Entry<BigFloat>* buf, int count, int to, int tag, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::TaggedSend [Entry<BigFloat>]"))
    PackedTaggedSend( buf, count, to, tag, comm );
}
#endif

template<typename Real>
void TaggedSend
( const Complex<Real>* buf, int count, int to, int tag, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::TaggedSend"))
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

template<typename T>
void Send( const T* buf, int count, int to, Comm comm ) EL_NO_RELEASE_EXCEPT
{ TaggedSend( buf, count, to, 0, comm ); }

template<typename T>
void TaggedSend( T b, int to, int tag, Comm comm ) EL_NO_RELEASE_EXCEPT
{ TaggedSend( &b, 1, to, tag, comm ); }

template<typename T>
void Send( T b, int to, Comm comm ) EL_NO_RELEASE_EXCEPT
{ TaggedSend( b, to, 0, comm ); }

template<typename Real>
void TaggedISend
( const Real* buf, int count, int to, int tag, Comm comm,
  Request<Real>& request )
EL_NO_RELEASE_EXCEPT
{ 
    DEBUG_ONLY(CSE cse("mpi::ISend"))
    SafeMpi
    ( MPI_Isend
      ( const_cast<Real*>(buf), count, TypeMap<Real>(), to, 
        tag, comm.comm, &request.backend ) );
}
#ifdef EL_HAVE_MPC
template<typename T>
void PackedTaggedISend
( const T* buf, int count, int to, int tag, Comm comm, Request<T>& request )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::PackedTaggedISend"))
    Serialize( count, buf, request.buffer );
    SafeMpi
    ( MPI_Isend
      ( request.buffer.data(), count, TypeMap<T>(), to, tag, comm.comm,
        &request.backend ) );
}

template<>
void TaggedISend
( const BigInt* buf, int count, int to, int tag, Comm comm,
  Request<BigInt>& request )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::TaggedISend [BigInt]"))
    PackedTaggedISend( buf, count, to, tag, comm, request );
}
template<>
void TaggedISend
( const ValueInt<BigInt>* buf, int count, int to, int tag, Comm comm,
  Request<ValueInt<BigInt>>& request )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::TaggedISend [ValueInt<BigInt>]"))
    PackedTaggedISend( buf, count, to, tag, comm, request );
}
template<>
void TaggedISend
( const Entry<BigInt>* buf, int count, int to, int tag, Comm comm,
  Request<Entry<BigInt>>& request )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::TaggedISend [Entry<BigInt>]"))
    PackedTaggedISend( buf, count, to, tag, comm, request );
}

template<>
void TaggedISend
( const BigFloat* buf, int count, int to, int tag, Comm comm,
  Request<BigFloat>& request )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::TaggedISend [BigFloat]"))
    PackedTaggedISend( buf, count, to, tag, comm, request );
}
template<>
void TaggedISend
( const ValueInt<BigFloat>* buf, int count, int to, int tag, Comm comm,
  Request<ValueInt<BigFloat>>& request )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::TaggedISend [ValueInt<BigFloat>]"))
    PackedTaggedISend( buf, count, to, tag, comm, request );
}
template<>
void TaggedISend
( const Entry<BigFloat>* buf, int count, int to, int tag, Comm comm,
  Request<Entry<BigFloat>>& request )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::TaggedISend [Entry<BigFloat>]"))
    PackedTaggedISend( buf, count, to, tag, comm, request );
}
#endif

template<typename Real>
void TaggedISend
( const Complex<Real>* buf, int count, int to, int tag, Comm comm, 
  Request<Complex<Real>>& request ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::ISend"))
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

template<typename Real>
void TaggedIRSend
( const Real* buf, int count, int to, int tag, Comm comm,
  Request<Real>& request )
EL_NO_RELEASE_EXCEPT
{ 
    DEBUG_ONLY(CSE cse("mpi::IRSend"))
    SafeMpi
    ( MPI_Irsend
      ( const_cast<Real*>(buf), count, TypeMap<Real>(), to, 
        tag, comm.comm, &request.backend ) );
}
#ifdef EL_HAVE_MPC
template<typename T>
void PackedTaggedIRSend
( const T* buf, int count, int to, int tag, Comm comm,
  Request<T>& request )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::PackedTaggedIRSend"))
    Serialize( count, buf, request.buffer );
    SafeMpi
    ( MPI_Irsend
      ( request.buffer.data(), count, TypeMap<T>(), to, 
        tag, comm.comm, &request.backend ) );
}

template<>
void TaggedIRSend
( const BigInt* buf, int count, int to, int tag, Comm comm,
  Request<BigInt>& request )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::TaggedIRSend [BigInt]"))
    PackedTaggedIRSend( buf, count, to, tag, comm, request );
}
template<>
void TaggedIRSend
( const ValueInt<BigInt>* buf, int count, int to, int tag, Comm comm,
  Request<ValueInt<BigInt>>& request )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::TaggedIRSend [ValueInt<BigInt>]"))
    PackedTaggedIRSend( buf, count, to, tag, comm, request );
}
template<>
void TaggedIRSend
( const Entry<BigInt>* buf, int count, int to, int tag, Comm comm,
  Request<Entry<BigInt>>& request )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::TaggedIRSend [Entry<BigInt>]"))
    PackedTaggedIRSend( buf, count, to, tag, comm, request );
}

template<>
void TaggedIRSend
( const BigFloat* buf, int count, int to, int tag, Comm comm,
  Request<BigFloat>& request )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::TaggedIRSend [BigFloat]"))
    PackedTaggedIRSend( buf, count, to, tag, comm, request );
}
template<>
void TaggedIRSend
( const ValueInt<BigFloat>* buf, int count, int to, int tag, Comm comm,
  Request<ValueInt<BigFloat>>& request )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::TaggedIRSend [ValueInt<BigFloat>]"))
    PackedTaggedIRSend( buf, count, to, tag, comm, request );
}
template<>
void TaggedIRSend
( const Entry<BigFloat>* buf, int count, int to, int tag, Comm comm,
  Request<Entry<BigFloat>>& request )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::TaggedIRSend [Entry<BigFloat>]"))
    PackedTaggedIRSend( buf, count, to, tag, comm, request );
}
#endif

template<typename Real>
void TaggedIRSend
( const Complex<Real>* buf, int count, int to, int tag, Comm comm, 
  Request<Complex<Real>>& request ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::IRSend"))
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

template<typename Real>
void TaggedISSend
( const Real* buf, int count, int to, int tag, Comm comm,
  Request<Real>& request ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::ISSend"))
    SafeMpi
    ( MPI_Issend
      ( const_cast<Real*>(buf), count, TypeMap<Real>(), to, 
        tag, comm.comm, &request.backend ) );
}
#ifdef EL_HAVE_MPC
template<typename T>
void PackedTaggedISSend
( const T* buf, int count, int to, int tag, Comm comm, Request<T>& request )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::PackedISSend"))
    Serialize( count, buf, request.buffer );
    SafeMpi
    ( MPI_Issend
      ( request.buffer.data(), count, TypeMap<T>(), to, 
        tag, comm.comm, &request.backend ) );
}

template<>
void TaggedISSend
( const BigInt* buf, int count, int to, int tag, Comm comm,
  Request<BigInt>& request ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::ISSend [BigInt]"))
    PackedTaggedISSend( buf, count, to, tag, comm, request );
}
template<>
void TaggedISSend
( const ValueInt<BigInt>* buf, int count, int to, int tag, Comm comm,
  Request<ValueInt<BigInt>>& request ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::ISSend [ValueInt<BigInt>]"))
    PackedTaggedISSend( buf, count, to, tag, comm, request );
}
template<>
void TaggedISSend
( const Entry<BigInt>* buf, int count, int to, int tag, Comm comm,
  Request<Entry<BigInt>>& request ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::ISSend [Entry<BigInt>]"))
    PackedTaggedISSend( buf, count, to, tag, comm, request );
}

template<>
void TaggedISSend
( const BigFloat* buf, int count, int to, int tag, Comm comm,
  Request<BigFloat>& request ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::ISSend [BigFloat]"))
    PackedTaggedISSend( buf, count, to, tag, comm, request );
}
template<>
void TaggedISSend
( const ValueInt<BigFloat>* buf, int count, int to, int tag, Comm comm,
  Request<ValueInt<BigFloat>>& request )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::ISSend [ValueInt<BigFloat>]"))
    PackedTaggedISSend( buf, count, to, tag, comm, request );
}
template<>
void TaggedISSend
( const Entry<BigFloat>* buf, int count, int to, int tag, Comm comm,
  Request<Entry<BigFloat>>& request ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::ISSend [Entry<BigFloat>]"))
    PackedTaggedISSend( buf, count, to, tag, comm, request );
}
#endif

template<typename Real>
void TaggedISSend
( const Complex<Real>* buf, int count, int to, int tag, Comm comm, 
  Request<Complex<Real>>& request ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::ISSend"))
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

template<typename T>
void ISSend( const T* buf, int count, int to, Comm comm, Request<T>& request )
EL_NO_RELEASE_EXCEPT
{ TaggedISSend( buf, count, to, 0, comm, request ); }

template<typename T>
void TaggedISSend( T b, int to, int tag, Comm comm, Request<T>& request )
EL_NO_RELEASE_EXCEPT
{ TaggedISSend( &b, 1, to, tag, comm, request ); }

template<typename Real>
void TaggedRecv( Real* buf, int count, int from, int tag, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::TaggedRecv"))
    Status status;
    SafeMpi
    ( MPI_Recv( buf, count, TypeMap<Real>(), from, tag, comm.comm, &status ) );
}
#ifdef EL_HAVE_MPC
template<typename T>
void PackedTaggedRecv( T* buf, int count, int from, int tag, Comm comm )
{
    DEBUG_ONLY(CSE cse("mpi::PackedTaggedRecv"))
    std::vector<byte> packedBuf;
    ReserveSerialized( count, buf, packedBuf );
    Status status;
    SafeMpi
    ( MPI_Recv
      ( packedBuf.data(), count, TypeMap<T>(), from, tag,
        comm.comm, &status ) );
    Deserialize( count, packedBuf, buf );
}

template<>
void TaggedRecv( BigInt* buf, int count, int from, int tag, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::TaggedRecv [BigInt]"))
    PackedTaggedRecv( buf, count, from, tag, comm );
}
template<>
void TaggedRecv
( ValueInt<BigInt>* buf, int count, int from, int tag, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::TaggedRecv [ValueInt<BigInt>]"))
    PackedTaggedRecv( buf, count, from, tag, comm );
}
template<>
void TaggedRecv( Entry<BigInt>* buf, int count, int from, int tag, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::TaggedRecv [Entry<BigInt>]"))
    PackedTaggedRecv( buf, count, from, tag, comm );
}

template<>
void TaggedRecv( BigFloat* buf, int count, int from, int tag, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::TaggedRecv [BigFloat]"))
    PackedTaggedRecv( buf, count, from, tag, comm );
}
template<>
void TaggedRecv
( ValueInt<BigFloat>* buf, int count, int from, int tag, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::TaggedRecv [ValueInt<BigFloat>]"))
    PackedTaggedRecv( buf, count, from, tag, comm );
}
template<>
void TaggedRecv( Entry<BigFloat>* buf, int count, int from, int tag, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::TaggedRecv [Entry<BigFloat>]"))
    PackedTaggedRecv( buf, count, from, tag, comm );
}
#endif

template<typename Real>
void TaggedRecv( Complex<Real>* buf, int count, int from, int tag, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Recv"))
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

template<typename Real>
void TaggedIRecv
( Real* buf, int count, int from, int tag, Comm comm, Request<Real>& request )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::TaggedIRecv"))
    SafeMpi
    ( MPI_Irecv
      ( buf, count, TypeMap<Real>(), from, tag, comm.comm, &request.backend ) );
}
#ifdef EL_HAVE_MPC
template<typename T>
void PackedTaggedIRecv
( T* buf, int count, int from, int tag, Comm comm, Request<T>& request )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::PackedTaggedIRecv"))
    request.receivingPacked = true;
    request.recvCount = count;
    request.unpackedRecvBuf = buf;
    ReserveSerialized( count, buf, request.buffer );
    SafeMpi
    ( MPI_Irecv
      ( request.buffer.data(), count, TypeMap<T>(), from, tag, comm.comm,
        &request.backend ) );
}

template<>
void TaggedIRecv
( BigInt* buf, int count, int from, int tag, Comm comm,
  Request<BigInt>& request ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::TaggedIRecv [BigInt]"))
    PackedTaggedIRecv( buf, count, from, tag, comm, request );
}
template<>
void TaggedIRecv
( ValueInt<BigInt>* buf, int count, int from, int tag, Comm comm,
  Request<ValueInt<BigInt>>& request )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::TaggedIRecv [ValueInt<BigInt>]"))
    PackedTaggedIRecv( buf, count, from, tag, comm, request );
}
template<>
void TaggedIRecv
( Entry<BigInt>* buf, int count, int from, int tag, Comm comm,
  Request<Entry<BigInt>>& request )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::TaggedIRecv [Entry<BigInt>]"))
    PackedTaggedIRecv( buf, count, from, tag, comm, request );
}

template<>
void TaggedIRecv
( BigFloat* buf, int count, int from, int tag, Comm comm,
  Request<BigFloat>& request ) EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::TaggedIRecv [BigFloat]"))
    PackedTaggedIRecv( buf, count, from, tag, comm, request );
}
template<>
void TaggedIRecv
( ValueInt<BigFloat>* buf, int count, int from, int tag, Comm comm,
  Request<ValueInt<BigFloat>>& request )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::TaggedIRecv [ValueInt<BigFloat>]"))
    PackedTaggedIRecv( buf, count, from, tag, comm, request );
}
template<>
void TaggedIRecv
( Entry<BigFloat>* buf, int count, int from, int tag, Comm comm,
  Request<Entry<BigFloat>>& request )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::TaggedIRecv [Entry<BigFloat>]"))
    PackedTaggedIRecv( buf, count, from, tag, comm, request );
}
#endif

template<typename Real>
void TaggedIRecv
( Complex<Real>* buf, int count, int from, int tag, Comm comm,
  Request<Complex<Real>>& request )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::IRecv"))
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

template<typename Real>
void TaggedSendRecv
( const Real* sbuf, int sc, int to,   int stag,
        Real* rbuf, int rc, int from, int rtag, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::TaggedSendRecv"))
    Status status;
    SafeMpi
    ( MPI_Sendrecv
      ( const_cast<Real*>(sbuf), sc, TypeMap<Real>(), to,   stag,
        rbuf,                    rc, TypeMap<Real>(), from, rtag, 
        comm.comm, &status ) );
}
#ifdef EL_HAVE_MPC
template<typename T>
void PackedTaggedSendRecv
( const T* sbuf, int sc, int to,   int stag,
        T* rbuf, int rc, int from, int rtag, Comm comm )
{
    DEBUG_ONLY(CSE cse("mpi::PackedTaggedSendRecv"))
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

template<>
void TaggedSendRecv
( const BigInt* sbuf, int sc, int to,   int stag,
        BigInt* rbuf, int rc, int from, int rtag, Comm comm )
{
    DEBUG_ONLY(CSE cse("mpi::TaggedSendRecv [BigInt]"))
    PackedTaggedSendRecv( sbuf, sc, to, stag, rbuf, rc, from, rtag, comm );
}
template<>
void TaggedSendRecv
( const ValueInt<BigInt>* sbuf, int sc, int to,   int stag,
        ValueInt<BigInt>* rbuf, int rc, int from, int rtag, Comm comm )
{
    DEBUG_ONLY(CSE cse("mpi::TaggedSendRecv [ValueInt<BigInt>]"))
    PackedTaggedSendRecv( sbuf, sc, to, stag, rbuf, rc, from, rtag, comm );
}
template<>
void TaggedSendRecv
( const Entry<BigInt>* sbuf, int sc, int to,   int stag,
        Entry<BigInt>* rbuf, int rc, int from, int rtag, Comm comm )
{
    DEBUG_ONLY(CSE cse("mpi::TaggedSendRecv [Entry<BigInt>]"))
    PackedTaggedSendRecv( sbuf, sc, to, stag, rbuf, rc, from, rtag, comm );
}

template<>
void TaggedSendRecv
( const BigFloat* sbuf, int sc, int to,   int stag,
        BigFloat* rbuf, int rc, int from, int rtag, Comm comm )
{
    DEBUG_ONLY(CSE cse("mpi::TaggedSendRecv [BigFloat]"))
    PackedTaggedSendRecv( sbuf, sc, to, stag, rbuf, rc, from, rtag, comm );
}
template<>
void TaggedSendRecv
( const ValueInt<BigFloat>* sbuf, int sc, int to,   int stag,
        ValueInt<BigFloat>* rbuf, int rc, int from, int rtag, Comm comm )
{
    DEBUG_ONLY(CSE cse("mpi::TaggedSendRecv [ValueInt<BigFloat>]"))
    PackedTaggedSendRecv( sbuf, sc, to, stag, rbuf, rc, from, rtag, comm );
}
template<>
void TaggedSendRecv
( const Entry<BigFloat>* sbuf, int sc, int to,   int stag,
        Entry<BigFloat>* rbuf, int rc, int from, int rtag, Comm comm )
{
    DEBUG_ONLY(CSE cse("mpi::TaggedSendRecv [Entry<BigFloat>]"))
    PackedTaggedSendRecv( sbuf, sc, to, stag, rbuf, rc, from, rtag, comm );
}
#endif

template<typename Real>
void TaggedSendRecv
( const Complex<Real>* sbuf, int sc, int to,   int stag,
        Complex<Real>* rbuf, int rc, int from, int rtag, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::SendRecv"))
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

template<typename Real>
void TaggedSendRecv
( Real* buf, int count, int to, int stag, int from, int rtag, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::SendRecv"))
    Status status;
    SafeMpi
    ( MPI_Sendrecv_replace
      ( buf, count, TypeMap<Real>(), to, stag, from, rtag, comm.comm,
        &status ) );
}
#ifdef EL_HAVE_MPC
template<typename T>
void PackedTaggedSendRecv
( T* buf, int count, int to, int stag, int from, int rtag, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::PackedTaggedSendRecv"))
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

template<>
void TaggedSendRecv
( BigInt* buf, int count, int to, int stag, int from, int rtag, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::TaggedSendRecv [BigInt]"))
    PackedTaggedSendRecv( buf, count, to, stag, from, rtag, comm );
}
template<>
void TaggedSendRecv
( ValueInt<BigInt>* buf, int count, int to, int stag, int from, int rtag,
  Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::TaggedSendRecv [ValueInt<BigInt>]"))
    PackedTaggedSendRecv( buf, count, to, stag, from, rtag, comm );
}
template<>
void TaggedSendRecv
( Entry<BigInt>* buf, int count, int to, int stag, int from, int rtag,
  Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::TaggedSendRecv [Entry<BigInt>]"))
    PackedTaggedSendRecv( buf, count, to, stag, from, rtag, comm );
}

template<>
void TaggedSendRecv
( BigFloat* buf, int count, int to, int stag, int from, int rtag, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::TaggedSendRecv [BigFloat]"))
    PackedTaggedSendRecv( buf, count, to, stag, from, rtag, comm );
}
template<>
void TaggedSendRecv
( ValueInt<BigFloat>* buf, int count, int to, int stag, int from, int rtag,
  Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::TaggedSendRecv [ValueInt<BigFloat>]"))
    PackedTaggedSendRecv( buf, count, to, stag, from, rtag, comm );
}
template<>
void TaggedSendRecv
( Entry<BigFloat>* buf, int count, int to, int stag, int from, int rtag,
  Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::TaggedSendRecv [Entry<BigFloat>]"))
    PackedTaggedSendRecv( buf, count, to, stag, from, rtag, comm );
}
#endif

template<typename Real>
void TaggedSendRecv
( Complex<Real>* buf, int count, int to, int stag, int from, int rtag,
  Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::SendRecv"))
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

template<typename T>
void SendRecv( T* buf, int count, int to, int from, Comm comm )
EL_NO_RELEASE_EXCEPT
{ TaggedSendRecv( buf, count, to, 0, from, ANY_TAG, comm ); }

template<typename Real>
void Broadcast( Real* buf, int count, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Broadcast"))
    if( Size(comm) == 1 || count == 0 )
        return;
    SafeMpi( MPI_Bcast( buf, count, TypeMap<Real>(), root, comm.comm ) );
}
#ifdef EL_HAVE_MPC
template<typename T>
void PackedBroadcast( T* buf, int count, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::PackedBroadcast"))
    if( Size(comm) == 1 || count == 0 )
        return;
    std::vector<byte> packedBuf;
    Serialize( count, buf, packedBuf );
    SafeMpi(
      MPI_Bcast( packedBuf.data(), count, TypeMap<T>(), root, comm.comm )
    );
    Deserialize( count, packedBuf, buf );
}

template<>
void Broadcast( BigInt* buf, int count, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Broadcast [BigInt]"))
    PackedBroadcast( buf, count, root, comm );
}
template<>
void Broadcast( ValueInt<BigInt>* buf, int count, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Broadcast [ValueInt<BigInt>]"))
    PackedBroadcast( buf, count, root, comm );
}
template<>
void Broadcast( Entry<BigInt>* buf, int count, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Broadcast [Entry<BigInt>]"))
    PackedBroadcast( buf, count, root, comm );
}

template<>
void Broadcast( BigFloat* buf, int count, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Broadcast [BigFloat]"))
    PackedBroadcast( buf, count, root, comm );
}
template<>
void Broadcast( ValueInt<BigFloat>* buf, int count, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Broadcast [ValueInt<BigFloat>]"))
    PackedBroadcast( buf, count, root, comm );
}
template<>
void Broadcast( Entry<BigFloat>* buf, int count, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Broadcast [Entry<BigFloat>]"))
    PackedBroadcast( buf, count, root, comm );
}
#endif

template<typename Real>
void Broadcast( Complex<Real>* buf, int count, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Broadcast"))
    if( Size(comm) == 1 )
        return;
#ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi( MPI_Bcast( buf, 2*count, TypeMap<Real>(), root, comm.comm ) );
#else
    SafeMpi( MPI_Bcast( buf, count, TypeMap<Complex<Real>>(), root, comm.comm ) );
#endif
}

template<typename T>
void Broadcast( T& b, int root, Comm comm ) EL_NO_RELEASE_EXCEPT
{ Broadcast( &b, 1, root, comm ); }

template<typename Real>
void IBroadcast
( Real* buf, int count, int root, Comm comm, Request<Real>& request )
{
    DEBUG_ONLY(CSE cse("mpi::IBroadcast"))
#ifdef EL_HAVE_NONBLOCKING_COLLECTIVES
    SafeMpi
    ( MPI_Ibcast
      ( buf, count, TypeMap<Real>(), root, comm.comm, &request.backend ) );
#else
    LogicError("Elemental was not configured with non-blocking support");
#endif
}
#ifdef EL_HAVE_MPC
template<typename T>
void PackedIBroadcast
( T* buf, int count, int root, Comm comm, Request<T>& request )
{
    DEBUG_ONLY(CSE cse("mpi::PackedIBroadcast"))
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

template<>
void IBroadcast
( BigInt* buf, int count, int root, Comm comm, Request<BigInt>& request )
{
    DEBUG_ONLY(CSE cse("mpi::IBroadcast [BigInt]"))
    PackedIBroadcast( buf, count, root, comm, request );
}
template<>
void IBroadcast
( ValueInt<BigInt>* buf, int count, int root, Comm comm,
  Request<ValueInt<BigInt>>& request )
{
    DEBUG_ONLY(CSE cse("mpi::IBroadcast [ValueInt<BigInt>]"))
    PackedIBroadcast( buf, count, root, comm, request );
}
template<>
void IBroadcast
( Entry<BigInt>* buf, int count, int root, Comm comm,
  Request<Entry<BigInt>>& request )
{
    DEBUG_ONLY(CSE cse("mpi::IBroadcast [Entry<BigInt>]"))
    PackedIBroadcast( buf, count, root, comm, request );
}

template<>
void IBroadcast
( BigFloat* buf, int count, int root, Comm comm, Request<BigFloat>& request )
{
    DEBUG_ONLY(CSE cse("mpi::IBroadcast [BigFloat]"))
    PackedIBroadcast( buf, count, root, comm, request );
}
template<>
void IBroadcast
( ValueInt<BigFloat>* buf, int count, int root, Comm comm,
  Request<ValueInt<BigFloat>>& request )
{
    DEBUG_ONLY(CSE cse("mpi::IBroadcast [ValueInt<BigFloat>]"))
    PackedIBroadcast( buf, count, root, comm, request );
}
template<>
void IBroadcast
( Entry<BigFloat>* buf, int count, int root, Comm comm,
  Request<Entry<BigFloat>>& request )
{
    DEBUG_ONLY(CSE cse("mpi::IBroadcast [Entry<BigFloat>]"))
    PackedIBroadcast( buf, count, root, comm, request );
}
#endif

template<typename Real>
void IBroadcast
( Complex<Real>* buf, int count, int root, Comm comm,
  Request<Complex<Real>>& request )
{
    DEBUG_ONLY(CSE cse("mpi::IBroadcast"))
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

template<typename T>
void IBroadcast( T& b, int root, Comm comm, Request<T>& request )
{ IBroadcast( &b, 1, root, comm, request ); }

template<typename Real>
void Gather
( const Real* sbuf, int sc,
        Real* rbuf, int rc, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Gather"))
    SafeMpi
    ( MPI_Gather
      ( const_cast<Real*>(sbuf), sc, TypeMap<Real>(),
        rbuf,                    rc, TypeMap<Real>(), root, comm.comm ) );
}
#ifdef EL_HAVE_MPC
template<typename T>
void PackedGather
( const T* sbuf, int sc,
        T* rbuf, int rc, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::PackedGather"))
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

template<>
void Gather
( const BigInt* sbuf, int sc,
        BigInt* rbuf, int rc, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Gather [BigInt]"))
    PackedGather( sbuf, sc, rbuf, rc, root, comm );
}
template<>
void Gather
( const ValueInt<BigInt>* sbuf, int sc,
        ValueInt<BigInt>* rbuf, int rc, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Gather [ValueInt<BigInt>]"))
    PackedGather( sbuf, sc, rbuf, rc, root, comm );
}
template<>
void Gather
( const Entry<BigInt>* sbuf, int sc,
        Entry<BigInt>* rbuf, int rc, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Gather [Entry<BigInt>]"))
    PackedGather( sbuf, sc, rbuf, rc, root, comm );
}

template<>
void Gather
( const BigFloat* sbuf, int sc,
        BigFloat* rbuf, int rc, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Gather [BigFloat]"))
    PackedGather( sbuf, sc, rbuf, rc, root, comm );
}
template<>
void Gather
( const ValueInt<BigFloat>* sbuf, int sc,
        ValueInt<BigFloat>* rbuf, int rc, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Gather [ValueInt<BigFloat>]"))
    PackedGather( sbuf, sc, rbuf, rc, root, comm );
}
template<>
void Gather
( const Entry<BigFloat>* sbuf, int sc,
        Entry<BigFloat>* rbuf, int rc, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Gather [Entry<BigFloat>]"))
    PackedGather( sbuf, sc, rbuf, rc, root, comm );
}
#endif

template<typename Real>
void Gather
( const Complex<Real>* sbuf, int sc,
        Complex<Real>* rbuf, int rc, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Gather"))
#ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi
    ( MPI_Gather
      ( const_cast<Complex<Real>*>(sbuf), 2*sc, TypeMap<Real>(),
        rbuf,                             2*rc, TypeMap<Real>(), root, comm.comm ) );
#else
    SafeMpi
    ( MPI_Gather
      ( const_cast<Complex<Real>*>(sbuf), sc, TypeMap<Complex<Real>>(),
        rbuf,                             rc, TypeMap<Complex<Real>>(), 
        root, comm.comm ) );
#endif
}

template<typename Real>
void IGather
( const Real* sbuf, int sc,
        Real* rbuf, int rc,
  int root, Comm comm,
  Request<Real>& request )
{
    DEBUG_ONLY(CSE cse("mpi::IGather"))
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
#ifdef EL_HAVE_MPC
template<typename T>
void PackedIGather
( const T* sbuf, int sc,
        T* rbuf, int rc,
  int root, Comm comm,
  Request<T>& request )
{
    DEBUG_ONLY(CSE cse("mpi::PackedIGather"))
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

template<>
void IGather
( const BigInt* sbuf, int sc,
        BigInt* rbuf, int rc,
  int root, Comm comm,
  Request<BigInt>& request )
{
    DEBUG_ONLY(CSE cse("mpi::IGather [BigInt]"))
    PackedIGather( sbuf, sc, rbuf, rc, root, comm, request );
}
template<>
void IGather
( const ValueInt<BigInt>* sbuf, int sc,
        ValueInt<BigInt>* rbuf, int rc,
  int root, Comm comm,
  Request<ValueInt<BigInt>>& request )
{
    DEBUG_ONLY(CSE cse("mpi::IGather [ValueInt<BigInt>]"))
    PackedIGather( sbuf, sc, rbuf, rc, root, comm, request );
}
template<>
void IGather
( const Entry<BigInt>* sbuf, int sc,
        Entry<BigInt>* rbuf, int rc,
  int root, Comm comm, Request<Entry<BigInt>>& request )
{
    DEBUG_ONLY(CSE cse("mpi::IGather [Entry<BigInt>]"))
    PackedIGather( sbuf, sc, rbuf, rc, root, comm, request );
}

template<>
void IGather
( const BigFloat* sbuf, int sc,
        BigFloat* rbuf, int rc,
  int root, Comm comm,
  Request<BigFloat>& request )
{
    DEBUG_ONLY(CSE cse("mpi::IGather [BigFloat]"))
    PackedIGather( sbuf, sc, rbuf, rc, root, comm, request );
}
template<>
void IGather
( const ValueInt<BigFloat>* sbuf, int sc,
        ValueInt<BigFloat>* rbuf, int rc,
  int root, Comm comm,
  Request<ValueInt<BigFloat>>& request )
{
    DEBUG_ONLY(CSE cse("mpi::IGather [ValueInt<BigFloat>]"))
    PackedIGather( sbuf, sc, rbuf, rc, root, comm, request );
}
template<>
void IGather
( const Entry<BigFloat>* sbuf, int sc,
        Entry<BigFloat>* rbuf, int rc,
  int root, Comm comm,
  Request<Entry<BigFloat>>& request )
{
    DEBUG_ONLY(CSE cse("mpi::IGather [Entry<BigFloat>]"))
    PackedIGather( sbuf, sc, rbuf, rc, root, comm, request );
}
#endif

template<typename Real>
void IGather
( const Complex<Real>* sbuf, int sc,
        Complex<Real>* rbuf, int rc,
  int root, Comm comm,
  Request<Complex<Real>>& request )
{
    DEBUG_ONLY(CSE cse("mpi::IGather"))
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

template<typename Real>
void Gather
( const Real* sbuf, int sc,
        Real* rbuf, const int* rcs, const int* rds,
  int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Gather"))
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
#ifdef EL_HAVE_MPC
template<typename T>
void PackedGather
( const T* sbuf, int sc,
        T* rbuf, const int* rcs, const int* rds,
  int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::PackedGather"))
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

template<>
void Gather
( const BigInt* sbuf, int sc,
        BigInt* rbuf, const int* rcs, const int* rds, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Gather [BigInt]"))
    PackedGather( sbuf, sc, rbuf, rcs, rds, root, comm );
}
template<>
void Gather
( const ValueInt<BigInt>* sbuf, int sc,
        ValueInt<BigInt>* rbuf, const int* rcs, const int* rds,
  int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Gather [ValueInt<BigInt>]"))
    PackedGather( sbuf, sc, rbuf, rcs, rds, root, comm );
}
template<>
void Gather
( const Entry<BigInt>* sbuf, int sc,
        Entry<BigInt>* rbuf, const int* rcs, const int* rds,
  int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Gather [Entry<BigInt>]"))
    PackedGather( sbuf, sc, rbuf, rcs, rds, root, comm );
}

template<>
void Gather
( const BigFloat* sbuf, int sc,
        BigFloat* rbuf, const int* rcs, const int* rds, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Gather [BigFloat]"))
    PackedGather( sbuf, sc, rbuf, rcs, rds, root, comm );
}
template<>
void Gather
( const ValueInt<BigFloat>* sbuf, int sc,
        ValueInt<BigFloat>* rbuf, const int* rcs, const int* rds,
  int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Gather [ValueInt<BigFloat>]"))
    PackedGather( sbuf, sc, rbuf, rcs, rds, root, comm );
}
template<>
void Gather
( const Entry<BigFloat>* sbuf, int sc,
        Entry<BigFloat>* rbuf, const int* rcs, const int* rds,
  int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Gather [Entry<BigFloat>]"))
    PackedGather( sbuf, sc, rbuf, rcs, rds, root, comm );
}
#endif

template<typename Real>
void Gather
( const Complex<Real>* sbuf, int sc,
        Complex<Real>* rbuf, const int* rcs, const int* rds, int root,
  Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Gather"))
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

template<typename Real>
void AllGather
( const Real* sbuf, int sc,
        Real* rbuf, int rc, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::AllGather"))
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
#ifdef EL_HAVE_MPC
template<typename T>
void PackedAllGather
( const T* sbuf, int sc,
        T* rbuf, int rc, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::PackedAllGather"))
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

template<>
void AllGather
( const BigInt* sbuf, int sc,
        BigInt* rbuf, int rc, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::AllGather [BigInt]"))
    PackedAllGather( sbuf, sc, rbuf, rc, comm );
}
template<>
void AllGather
( const ValueInt<BigInt>* sbuf, int sc,
        ValueInt<BigInt>* rbuf, int rc, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::AllGather [ValueInt<BigInt>]"))
    PackedAllGather( sbuf, sc, rbuf, rc, comm );
}
template<>
void AllGather
( const Entry<BigInt>* sbuf, int sc,
        Entry<BigInt>* rbuf, int rc, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::AllGather [Entry<BigInt>]"))
    PackedAllGather( sbuf, sc, rbuf, rc, comm );
}

template<>
void AllGather
( const BigFloat* sbuf, int sc,
        BigFloat* rbuf, int rc, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::AllGather [BigFloat]"))
    PackedAllGather( sbuf, sc, rbuf, rc, comm );
}
template<>
void AllGather
( const ValueInt<BigFloat>* sbuf, int sc,
        ValueInt<BigFloat>* rbuf, int rc, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::AllGather [ValueInt<BigFloat>]"))
    PackedAllGather( sbuf, sc, rbuf, rc, comm );
}
template<>
void AllGather
( const Entry<BigFloat>* sbuf, int sc,
        Entry<BigFloat>* rbuf, int rc, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::AllGather [Entry<BigFloat>]"))
    PackedAllGather( sbuf, sc, rbuf, rc, comm );
}
#endif

template<typename Real>
void AllGather
( const Complex<Real>* sbuf, int sc,
        Complex<Real>* rbuf, int rc, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::AllGather"))
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
        rbuf,                             2*rc, TypeMap<Real>(), comm.comm ) );
 #else
    SafeMpi
    ( MPI_Allgather
      ( const_cast<Complex<Real>*>(sbuf), sc, TypeMap<Complex<Real>>(),
        rbuf,                             rc, TypeMap<Complex<Real>>(),
        comm.comm ) );
 #endif
#endif
}

template<typename Real>
void AllGather
( const Real* sbuf, int sc,
        Real* rbuf, const int* rcs, const int* rds, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::AllGather"))
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
#ifdef EL_HAVE_MPC
template<typename T>
void PackedAllGather
( const T* sbuf, int sc,
        T* rbuf, const int* rcs, const int* rds, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::PackedAllgather"))
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

template<>
void AllGather
( const BigInt* sbuf, int sc,
        BigInt* rbuf, const int* rcs, const int* rds, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Allgather [BigInt]"))
    PackedAllGather( sbuf, sc, rbuf, rcs, rds, comm );
}
template<>
void AllGather
( const ValueInt<BigInt>* sbuf, int sc,
        ValueInt<BigInt>* rbuf, const int* rcs, const int* rds, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Allgather [ValueInt<BigInt>]"))
    PackedAllGather( sbuf, sc, rbuf, rcs, rds, comm );
}
template<>
void AllGather
( const Entry<BigInt>* sbuf, int sc,
        Entry<BigInt>* rbuf, const int* rcs, const int* rds, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Allgather [Entry<BigInt>]"))
    PackedAllGather( sbuf, sc, rbuf, rcs, rds, comm );
}

template<>
void AllGather
( const BigFloat* sbuf, int sc,
        BigFloat* rbuf, const int* rcs, const int* rds, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Allgather [BigFloat]"))
    PackedAllGather( sbuf, sc, rbuf, rcs, rds, comm );
}
template<>
void AllGather
( const ValueInt<BigFloat>* sbuf, int sc,
        ValueInt<BigFloat>* rbuf, const int* rcs, const int* rds, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Allgather [ValueInt<BigFloat>]"))
    PackedAllGather( sbuf, sc, rbuf, rcs, rds, comm );
}
template<>
void AllGather
( const Entry<BigFloat>* sbuf, int sc,
        Entry<BigFloat>* rbuf, const int* rcs, const int* rds, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Allgather [Entry<BigFloat>]"))
    PackedAllGather( sbuf, sc, rbuf, rcs, rds, comm );
}
#endif

template<typename Real>
void AllGather
( const Complex<Real>* sbuf, int sc,
        Complex<Real>* rbuf, const int* rcs, const int* rds, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::AllGather"))
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

template<typename Real>
void Scatter
( const Real* sbuf, int sc,
        Real* rbuf, int rc, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Scatter"))
    SafeMpi
    ( MPI_Scatter
      ( const_cast<Real*>(sbuf), sc, TypeMap<Real>(),
        rbuf,                    rc, TypeMap<Real>(), root, comm.comm ) );
}
#ifdef EL_HAVE_MPC
template<typename T>
void PackedScatter
( const T* sbuf, int sc,
        T* rbuf, int rc, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::PackedScatter"))
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

template<>
void Scatter
( const BigInt* sbuf, int sc,
        BigInt* rbuf, int rc, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Scatter [BigInt]"))
    PackedScatter( sbuf, sc, rbuf, rc, root, comm );
}
template<>
void Scatter
( const ValueInt<BigInt>* sbuf, int sc,
        ValueInt<BigInt>* rbuf, int rc, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Scatter [ValueInt<BigInt>]"))
    PackedScatter( sbuf, sc, rbuf, rc, root, comm );
}
template<>
void Scatter
( const Entry<BigInt>* sbuf, int sc,
        Entry<BigInt>* rbuf, int rc, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Scatter [Entry<BigInt>]"))
    PackedScatter( sbuf, sc, rbuf, rc, root, comm );
}

template<>
void Scatter
( const BigFloat* sbuf, int sc,
        BigFloat* rbuf, int rc, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Scatter [BigFloat]"))
    PackedScatter( sbuf, sc, rbuf, rc, root, comm );
}
template<>
void Scatter
( const ValueInt<BigFloat>* sbuf, int sc,
        ValueInt<BigFloat>* rbuf, int rc, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Scatter [ValueInt<BigFloat>]"))
    PackedScatter( sbuf, sc, rbuf, rc, root, comm );
}
template<>
void Scatter
( const Entry<BigFloat>* sbuf, int sc,
        Entry<BigFloat>* rbuf, int rc, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Scatter [Entry<BigFloat>]"))
    PackedScatter( sbuf, sc, rbuf, rc, root, comm );
}
#endif

template<typename Real>
void Scatter
( const Complex<Real>* sbuf, int sc,
        Complex<Real>* rbuf, int rc, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Scatter"))
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

template<typename Real>
void Scatter( Real* buf, int sc, int rc, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Scatter"))
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
#ifdef EL_HAVE_MPC
template<typename T>
void PackedScatter( T* buf, int sc, int rc, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::PackedScatter"))
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

template<>
void Scatter( BigInt* buf, int sc, int rc, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Scatter [BigInt]"))
    PackedScatter( buf, sc, rc, root, comm );
}
template<>
void Scatter( ValueInt<BigInt>* buf, int sc, int rc, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Scatter [ValueInt<BigInt>]"))
    PackedScatter( buf, sc, rc, root, comm );
}
template<>
void Scatter( Entry<BigInt>* buf, int sc, int rc, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Scatter [Entry<BigInt>]"))
    PackedScatter( buf, sc, rc, root, comm );
}

template<>
void Scatter( BigFloat* buf, int sc, int rc, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Scatter [BigFloat]"))
    PackedScatter( buf, sc, rc, root, comm );
}
template<>
void Scatter( ValueInt<BigFloat>* buf, int sc, int rc, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Scatter [ValueInt<BigFloat>]"))
    PackedScatter( buf, sc, rc, root, comm );
}
template<>
void Scatter( Entry<BigFloat>* buf, int sc, int rc, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Scatter [Entry<BigFloat>]"))
    PackedScatter( buf, sc, rc, root, comm );
}
#endif

template<typename Real>
void Scatter( Complex<Real>* buf, int sc, int rc, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Scatter"))
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

template<typename Real>
void AllToAll
( const Real* sbuf, int sc,
        Real* rbuf, int rc, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::AllToAll"))
    SafeMpi
    ( MPI_Alltoall
      ( const_cast<Real*>(sbuf), sc, TypeMap<Real>(),
        rbuf,                    rc, TypeMap<Real>(), comm.comm ) );
}
#ifdef EL_HAVE_MPC
template<typename T>
void PackedAllToAll
( const T* sbuf, int sc,
        T* rbuf, int rc, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::PackedAllToAll"))
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

template<>
void AllToAll
( const BigInt* sbuf, int sc,
        BigInt* rbuf, int rc, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::AllToAll [BigInt]"))
    PackedAllToAll( sbuf, sc, rbuf, rc, comm );
}
template<>
void AllToAll
( const ValueInt<BigInt>* sbuf, int sc,
        ValueInt<BigInt>* rbuf, int rc, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::AllToAll [ValueInt<BigInt>]"))
    PackedAllToAll( sbuf, sc, rbuf, rc, comm );
}
template<>
void AllToAll
( const Entry<BigInt>* sbuf, int sc,
        Entry<BigInt>* rbuf, int rc, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::AllToAll [Entry<BigInt>]"))
    PackedAllToAll( sbuf, sc, rbuf, rc, comm );
}

template<>
void AllToAll
( const BigFloat* sbuf, int sc,
        BigFloat* rbuf, int rc, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::AllToAll [BigFloat]"))
    PackedAllToAll( sbuf, sc, rbuf, rc, comm );
}
template<>
void AllToAll
( const ValueInt<BigFloat>* sbuf, int sc,
        ValueInt<BigFloat>* rbuf, int rc, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::AllToAll [ValueInt<BigFloat>]"))
    PackedAllToAll( sbuf, sc, rbuf, rc, comm );
}
template<>
void AllToAll
( const Entry<BigFloat>* sbuf, int sc,
        Entry<BigFloat>* rbuf, int rc, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::AllToAll [Entry<BigFloat>]"))
    PackedAllToAll( sbuf, sc, rbuf, rc, comm );
}
#endif

template<typename Real>
void AllToAll
( const Complex<Real>* sbuf, int sc,
        Complex<Real>* rbuf, int rc, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::AllToAll"))
#ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi
    ( MPI_Alltoall
      ( const_cast<Complex<Real>*>(sbuf), 2*sc, TypeMap<Real>(),
        rbuf,                             2*rc, TypeMap<Real>(), comm.comm ) );
#else
    SafeMpi
    ( MPI_Alltoall
      ( const_cast<Complex<Real>*>(sbuf), sc, TypeMap<Complex<Real>>(),
        rbuf,                             rc, TypeMap<Complex<Real>>(), comm.comm ) );
#endif
}

template<typename Real>
void AllToAll
( const Real* sbuf, const int* scs, const int* sds, 
        Real* rbuf, const int* rcs, const int* rds, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::AllToAll"))
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

#ifdef EL_HAVE_MPC
template<typename T>
void PackedAllToAll
( const T* sbuf, const int* scs, const int* sds,
        T* rbuf, const int* rcs, const int* rds, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::PackedAllToAll"))
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

template<>
void AllToAll
( const BigInt* sbuf, const int* scs, const int* sds,
        BigInt* rbuf, const int* rcs, const int* rds, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::AllToAll [BigInt]"))
    PackedAllToAll( sbuf, scs, sds, rbuf, rcs, rds, comm );
}
template<>
void AllToAll
( const ValueInt<BigInt>* sbuf, const int* scs, const int* sds,
        ValueInt<BigInt>* rbuf, const int* rcs, const int* rds, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::AllToAll [ValueInt<BigInt>]"))
    PackedAllToAll( sbuf, scs, sds, rbuf, rcs, rds, comm );
}
template<>
void AllToAll
( const Entry<BigInt>* sbuf, const int* scs, const int* sds,
        Entry<BigInt>* rbuf, const int* rcs, const int* rds, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::AllToAll [Entry<BigInt>]"))
    PackedAllToAll( sbuf, scs, sds, rbuf, rcs, rds, comm );
}

template<>
void AllToAll
( const BigFloat* sbuf, const int* scs, const int* sds,
        BigFloat* rbuf, const int* rcs, const int* rds, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::AllToAll [BigFloat]"))
    PackedAllToAll( sbuf, scs, sds, rbuf, rcs, rds, comm );
}
template<>
void AllToAll
( const ValueInt<BigFloat>* sbuf, const int* scs, const int* sds,
        ValueInt<BigFloat>* rbuf, const int* rcs, const int* rds, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::AllToAll [ValueInt<BigFloat>]"))
    PackedAllToAll( sbuf, scs, sds, rbuf, rcs, rds, comm );
}
template<>
void AllToAll
( const Entry<BigFloat>* sbuf, const int* scs, const int* sds,
        Entry<BigFloat>* rbuf, const int* rcs, const int* rds, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::AllToAll [Entry<BigFloat>]"))
    PackedAllToAll( sbuf, scs, sds, rbuf, rcs, rds, comm );
}
#endif

template<typename Real>
void AllToAll
( const Complex<Real>* sbuf, const int* scs, const int* sds,
        Complex<Real>* rbuf, const int* rcs, const int* rds, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::AllToAll"))
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

template<typename Real>
void Reduce
( const Real* sbuf, Real* rbuf, int count, Op op, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Reduce"))
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

#ifdef EL_HAVE_MPC
template<typename T>
void PackedReduce
( const T* sbuf, T* rbuf, int count, Op op, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::PackedReduce"))
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

template<>
void Reduce
( const BigInt* sbuf, BigInt* rbuf, int count, Op op, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Reduce [BigInt]"))
    PackedReduce( sbuf, rbuf, count, op, root, comm );
}
template<>
void Reduce
( const ValueInt<BigInt>* sbuf,
        ValueInt<BigInt>* rbuf, int count, Op op, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Reduce [ValueInt<BigInt>]"))
    PackedReduce( sbuf, rbuf, count, op, root, comm );
}
template<>
void Reduce
( const Entry<BigInt>* sbuf,
        Entry<BigInt>* rbuf, int count, Op op, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Reduce [Entry<BigInt>]"))
    PackedReduce( sbuf, rbuf, count, op, root, comm );
}

template<>
void Reduce
( const BigFloat* sbuf, BigFloat* rbuf, int count, Op op, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Reduce [BigFloat]"))
    PackedReduce( sbuf, rbuf, count, op, root, comm );
}
template<>
void Reduce
( const ValueInt<BigFloat>* sbuf,
        ValueInt<BigFloat>* rbuf, int count, Op op, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Reduce [ValueInt<BigFloat>]"))
    PackedReduce( sbuf, rbuf, count, op, root, comm );
}
template<>
void Reduce
( const Entry<BigFloat>* sbuf,
        Entry<BigFloat>* rbuf, int count, Op op, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Reduce [Entry<BigFloat>]"))
    PackedReduce( sbuf, rbuf, count, op, root, comm );
}
#endif

template<typename Real>
void Reduce
( const Complex<Real>* sbuf, 
        Complex<Real>* rbuf, int count, Op op, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Reduce"))
    if( count == 0 )
        return;

    MPI_Op opC;
    if( op == SUM )
        opC = SumOp<Complex<Real>>().op; 
    else
        opC = op.op;

#ifdef EL_AVOID_COMPLEX_MPI
    if( op == SUM )
    {
        SafeMpi
        ( MPI_Reduce
          ( const_cast<Complex<Real>*>(sbuf),
            rbuf, 2*count, TypeMap<Real>(), opC, 
            root, comm.comm ) );
    }
    else
    {
        SafeMpi
        ( MPI_Reduce
          ( const_cast<Complex<Real>*>(sbuf),
            rbuf, count, TypeMap<Complex<Real>>(), opC, root, comm.comm ) );
    }
#else
    SafeMpi
    ( MPI_Reduce
      ( const_cast<Complex<Real>*>(sbuf), 
        rbuf, count, TypeMap<Complex<Real>>(), opC, root, comm.comm ) );
#endif
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

template<typename Real>
void Reduce( Real* buf, int count, Op op, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Reduce"))
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

#ifdef EL_HAVE_MPC
template<typename T>
void PackedReduce( T* buf, int count, Op op, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::PackedReduce"))
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

template<>
void Reduce( BigInt* buf, int count, Op op, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Reduce [BigInt]"))
    PackedReduce( buf, count, op, root, comm );
}
template<>
void Reduce( ValueInt<BigInt>* buf, int count, Op op, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Reduce [ValueInt<BigInt>]"))
    PackedReduce( buf, count, op, root, comm );
}
template<>
void Reduce( Entry<BigInt>* buf, int count, Op op, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Reduce [Entry<BigInt>]"))
    PackedReduce( buf, count, op, root, comm );
}

template<>
void Reduce( BigFloat* buf, int count, Op op, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Reduce [BigFloat]"))
    PackedReduce( buf, count, op, root, comm );
}
template<>
void Reduce( ValueInt<BigFloat>* buf, int count, Op op, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Reduce [ValueInt<BigFloat>]"))
    PackedReduce( buf, count, op, root, comm );
}
template<>
void Reduce( Entry<BigFloat>* buf, int count, Op op, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Reduce [Entry<BigFloat>]"))
    PackedReduce( buf, count, op, root, comm );
}
#endif

template<typename Real>
void Reduce( Complex<Real>* buf, int count, Op op, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Reduce"))
    if( Size(comm) == 1 )
        return;
    if( count != 0 )
    {
        MPI_Op opC;
        if( op == SUM )
            opC = SumOp<Complex<Real>>().op; 
        else
            opC = op.op;

        const int commRank = Rank( comm );
#ifdef EL_AVOID_COMPLEX_MPI
        if( op == SUM )
        {
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

template<typename T>
void Reduce( T* buf, int count, int root, Comm comm )
EL_NO_RELEASE_EXCEPT
{ Reduce( buf, count, SUM, root, comm ); }

template<typename Real>
void AllReduce( const Real* sbuf, Real* rbuf, int count, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::AllReduce"))
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

#ifdef EL_HAVE_MPC
template<typename T>
void PackedAllReduce
( const T* sbuf, T* rbuf, int count, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::PackedAllReduce"))
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

template<>
void AllReduce
( const BigInt* sbuf, BigInt* rbuf, int count, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::AllReduce [BigInt]"))
    PackedAllReduce( sbuf, rbuf, count, op, comm );
}
template<>
void AllReduce
( const ValueInt<BigInt>* sbuf,
        ValueInt<BigInt>* rbuf, int count, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::AllReduce [ValueInt<BigInt>]"))
    PackedAllReduce( sbuf, rbuf, count, op, comm );
}
template<>
void AllReduce
( const Entry<BigInt>* sbuf,
        Entry<BigInt>* rbuf, int count, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::AllReduce [Entry<BigInt>]"))
    PackedAllReduce( sbuf, rbuf, count, op, comm );
}

template<>
void AllReduce
( const BigFloat* sbuf, BigFloat* rbuf, int count, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::AllReduce [BigFloat]"))
    PackedAllReduce( sbuf, rbuf, count, op, comm );
}
template<>
void AllReduce
( const ValueInt<BigFloat>* sbuf,
        ValueInt<BigFloat>* rbuf, int count, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::AllReduce [ValueInt<BigFloat>]"))
    PackedAllReduce( sbuf, rbuf, count, op, comm );
}
template<>
void AllReduce
( const Entry<BigFloat>* sbuf,
        Entry<BigFloat>* rbuf, int count, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::AllReduce [Entry<BigFloat>]"))
    PackedAllReduce( sbuf, rbuf, count, op, comm );
}
#endif

template<typename Real>
void AllReduce
( const Complex<Real>* sbuf, Complex<Real>* rbuf, int count, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::AllReduce"))
    if( count != 0 )
    {
        MPI_Op opC;
        if( op == SUM )
            opC = SumOp<Complex<Real>>().op; 
        else
            opC = op.op;

#ifdef EL_AVOID_COMPLEX_MPI
        if( op == SUM )
        {
            SafeMpi
            ( MPI_Allreduce
                ( const_cast<Complex<Real>*>(sbuf),
                  rbuf, 2*count, TypeMap<Real>(), opC, comm.comm ) );
        }
        else
        {
            SafeMpi
            ( MPI_Allreduce
              ( const_cast<Complex<Real>*>(sbuf),
                rbuf, count, TypeMap<Complex<Real>>(), opC, comm.comm ) );
        }
#else
        SafeMpi
        ( MPI_Allreduce
          ( const_cast<Complex<Real>*>(sbuf), 
            rbuf, count, TypeMap<Complex<Real>>(), opC, comm.comm ) );
#endif
    }
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

template<typename Real>
void AllReduce( Real* buf, int count, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::AllReduce"))
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

#ifdef EL_HAVE_MPC
template<typename T>
void PackedAllReduce( T* buf, int count, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::PackedAllReduce"))
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

template<>
void AllReduce( BigInt* buf, int count, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::AllReduce [BigInt]"))
    PackedAllReduce( buf, count, op, comm );
}
template<>
void AllReduce( ValueInt<BigInt>* buf, int count, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::AllReduce [ValueInt<BigInt>]"))
    PackedAllReduce( buf, count, op, comm );
}
template<>
void AllReduce( Entry<BigInt>* buf, int count, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::AllReduce [Entry<BigInt>]"))
    PackedAllReduce( buf, count, op, comm );
}

template<>
void AllReduce( BigFloat* buf, int count, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::AllReduce [BigFloat]"))
    PackedAllReduce( buf, count, op, comm );
}
template<>
void AllReduce( ValueInt<BigFloat>* buf, int count, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::AllReduce [ValueInt<BigFloat>]"))
    PackedAllReduce( buf, count, op, comm );
}
template<>
void AllReduce( Entry<BigFloat>* buf, int count, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::AllReduce [Entry<BigFloat>]"))
    PackedAllReduce( buf, count, op, comm );
}
#endif

template<typename Real>
void AllReduce( Complex<Real>* buf, int count, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::AllReduce"))
    if( count == 0 || Size(comm) == 1 )
        return;

    MPI_Op opC;
    if( op == SUM )
        opC = SumOp<Complex<Real>>().op; 
    else
        opC = op.op;

#ifdef EL_AVOID_COMPLEX_MPI
    if( op == SUM )
    {
        SafeMpi
        ( MPI_Allreduce
          ( MPI_IN_PLACE, buf, 2*count, TypeMap<Real>(), opC, comm.comm ) );
    }
    else
    {
        SafeMpi
        ( MPI_Allreduce
          ( MPI_IN_PLACE, buf, count, TypeMap<Complex<Real>>(), 
            opC, comm.comm ) );
    }
#else
    SafeMpi
    ( MPI_Allreduce
      ( MPI_IN_PLACE, buf, count, TypeMap<Complex<Real>>(), opC, 
        comm.comm ) );
#endif
}

template<typename T>
void AllReduce( T* buf, int count, Comm comm )
EL_NO_RELEASE_EXCEPT
{ AllReduce( buf, count, SUM, comm ); }

template<typename Real>
void ReduceScatter( Real* sbuf, Real* rbuf, int rc, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::ReduceScatter"))
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

#ifdef EL_HAVE_MPC
template<typename T>
void PackedReduceScatter( T* sbuf, T* rbuf, int rc, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::PackedReduceScatter"))
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

template<>
void ReduceScatter( BigInt* sbuf, BigInt* rbuf, int rc, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::ReduceScatter [BigInt]"))
    PackedReduceScatter( sbuf, rbuf, rc, op, comm );
}
template<>
void ReduceScatter
( ValueInt<BigInt>* sbuf,
  ValueInt<BigInt>* rbuf, int rc, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::ReduceScatter [ValueInt<BigInt>]"))
    PackedReduceScatter( sbuf, rbuf, rc, op, comm );
}
template<>
void ReduceScatter
( Entry<BigInt>* sbuf,
  Entry<BigInt>* rbuf, int rc, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::ReduceScatter [Entry<BigInt>]"))
    PackedReduceScatter( sbuf, rbuf, rc, op, comm );
}

template<>
void ReduceScatter( BigFloat* sbuf, BigFloat* rbuf, int rc, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::ReduceScatter [BigFloat]"))
    PackedReduceScatter( sbuf, rbuf, rc, op, comm );
}
template<>
void ReduceScatter
( ValueInt<BigFloat>* sbuf,
  ValueInt<BigFloat>* rbuf, int rc, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::ReduceScatter [ValueInt<BigFloat>]"))
    PackedReduceScatter( sbuf, rbuf, rc, op, comm );
}
template<>
void ReduceScatter
( Entry<BigFloat>* sbuf,
  Entry<BigFloat>* rbuf, int rc, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::ReduceScatter [Entry<BigFloat>]"))
    PackedReduceScatter( sbuf, rbuf, rc, op, comm );
}
#endif

template<typename Real>
void ReduceScatter
( Complex<Real>* sbuf, Complex<Real>* rbuf, int rc, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::ReduceScatter"))
    if( rc == 0 )
        return;
    MPI_Op opC;
    if( op == SUM )
        opC = SumOp<Complex<Real>>().op; 
    else
        opC = op.op;

#ifdef EL_REDUCE_SCATTER_BLOCK_VIA_ALLREDUCE
    const int commSize = Size( comm );
    const int commRank = Rank( comm );
    AllReduce( sbuf, rc*commSize, opC, comm );
    MemCopy( rbuf, &sbuf[commRank*rc], rc );
#elif defined(EL_HAVE_MPI_REDUCE_SCATTER_BLOCK)
# ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi
    ( MPI_Reduce_scatter_block
      ( sbuf, rbuf, 2*rc, TypeMap<Real>(), opC, comm.comm ) );
# else
    SafeMpi
    ( MPI_Reduce_scatter_block
      ( sbuf, rbuf, rc, TypeMap<Complex<Real>>(), opC, comm.comm ) );
# endif
#else
    const int commSize = Size( comm );
    Reduce( sbuf, rc*commSize, opC, 0, comm );
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

template<typename Real>
void ReduceScatter( Real* buf, int rc, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::ReduceScatter"))
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

#ifdef EL_HAVE_MPC
template<typename T>
void PackedReduceScatter( T* buf, int rc, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::PackedReduceScatter"))
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

template<>
void ReduceScatter( BigInt* buf, int rc, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::ReduceScatter [BigInt]"))
    PackedReduceScatter( buf, rc, op, comm );
}
template<>
void ReduceScatter( ValueInt<BigInt>* buf, int rc, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::ReduceScatter [ValueInt<BigInt>]"))
    PackedReduceScatter( buf, rc, op, comm );
}
template<>
void ReduceScatter( Entry<BigInt>* buf, int rc, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::ReduceScatter [Entry<BigInt>]"))
    PackedReduceScatter( buf, rc, op, comm );
}

template<>
void ReduceScatter( BigFloat* buf, int rc, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::ReduceScatter [BigFloat]"))
    PackedReduceScatter( buf, rc, op, comm );
}
template<>
void ReduceScatter( ValueInt<BigFloat>* buf, int rc, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::ReduceScatter [ValueInt<BigFloat>]"))
    PackedReduceScatter( buf, rc, op, comm );
}
template<>
void ReduceScatter( Entry<BigFloat>* buf, int rc, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::ReduceScatter [Entry<BigFloat>]"))
    PackedReduceScatter( buf, rc, op, comm );
}
#endif

// TODO: Handle case where op is not summation
template<typename Real>
void ReduceScatter( Complex<Real>* buf, int rc, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::ReduceScatter"))
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
        opC = SumOp<Complex<Real>>().op; 
    else
        opC = op.op;
# ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi
    ( MPI_Reduce_scatter_block
      ( MPI_IN_PLACE, buf, 2*rc, TypeMap<Real>(), opC, comm.comm ) );
# else
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

template<typename T>
void ReduceScatter( T* buf, int rc, Comm comm )
EL_NO_RELEASE_EXCEPT
{ ReduceScatter( buf, rc, SUM, comm ); }

template<typename Real>
void ReduceScatter
( const Real* sbuf, Real* rbuf, const int* rcs, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::ReduceScatter"))
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

#ifdef EL_HAVE_MPC
template<typename T>
void PackedReduceScatter
( const T* sbuf, T* rbuf, const int* rcs, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::PackedReduceScatter"))
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

template<>
void ReduceScatter
( const BigInt* sbuf, BigInt* rbuf, const int* rcs, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::ReduceScatter [BigInt]"))
    PackedReduceScatter( sbuf, rbuf, rcs, op, comm );
}
template<>
void ReduceScatter
( const ValueInt<BigInt>* sbuf,
        ValueInt<BigInt>* rbuf, const int* rcs, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::ReduceScatter [ValueInt<BigInt>]"))
    PackedReduceScatter( sbuf, rbuf, rcs, op, comm );
}
template<>
void ReduceScatter
( const Entry<BigInt>* sbuf,
        Entry<BigInt>* rbuf, const int* rcs, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::ReduceScatter [Entry<BigInt>]"))
    PackedReduceScatter( sbuf, rbuf, rcs, op, comm );
}

template<>
void ReduceScatter
( const BigFloat* sbuf, BigFloat* rbuf, const int* rcs, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::ReduceScatter [BigFloat]"))
    PackedReduceScatter( sbuf, rbuf, rcs, op, comm );
}
template<>
void ReduceScatter
( const ValueInt<BigFloat>* sbuf,
        ValueInt<BigFloat>* rbuf, const int* rcs, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::ReduceScatter [ValueInt<BigFloat>]"))
    PackedReduceScatter( sbuf, rbuf, rcs, op, comm );
}
template<>
void ReduceScatter
( const Entry<BigFloat>* sbuf,
        Entry<BigFloat>* rbuf, const int* rcs, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::ReduceScatter [Entry<BigFloat>]"))
    PackedReduceScatter( sbuf, rbuf, rcs, op, comm );
}
#endif

template<typename Real>
void ReduceScatter
( const Complex<Real>* sbuf, Complex<Real>* rbuf, const int* rcs,
  Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::ReduceScatter"))
    MPI_Op opC;
    if( op == SUM )
        opC = SumOp<Complex<Real>>().op; 
    else
        opC = op.op;

#ifdef EL_AVOID_COMPLEX_MPI
    if( op == SUM )
    {
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
        SafeMpi
        ( MPI_Reduce_scatter
          ( const_cast<Complex<Real>*>(sbuf),
            rbuf, const_cast<int*>(rcs), TypeMap<Complex<Real>>(), 
            opC, comm.comm ) );
    }
#else
    SafeMpi
    ( MPI_Reduce_scatter
      ( const_cast<Complex<Real>*>(sbuf), 
        rbuf, const_cast<int*>(rcs), TypeMap<Complex<Real>>(), opC, 
        comm.comm ) );
#endif
}

template<typename T>
void ReduceScatter( const T* sbuf, T* rbuf, const int* rcs, Comm comm )
EL_NO_RELEASE_EXCEPT
{ ReduceScatter( sbuf, rbuf, rcs, SUM, comm ); }

void VerifySendsAndRecvs
( const vector<int>& sendCounts,
  const vector<int>& recvCounts, Comm comm )
{
    DEBUG_ONLY(CSE cse("mpi::VerifySendsAndRecvs"))
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

template<typename Real>
void Scan( const Real* sbuf, Real* rbuf, int count, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Scan"))
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

#ifdef EL_HAVE_MPC
template<typename T>
void PackedScan( const T* sbuf, T* rbuf, int count, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::PackedScan"))
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

template<>
void Scan
( const BigInt* sbuf, BigInt* rbuf, int count, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Scan [BigInt]"))
    PackedScan( sbuf, rbuf, count, op, comm );
}
template<>
void Scan
( const ValueInt<BigInt>* sbuf,
        ValueInt<BigInt>* rbuf, int count, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Scan [ValueInt<BigInt>]"))
    PackedScan( sbuf, rbuf, count, op, comm );
}
template<>
void Scan
( const Entry<BigInt>* sbuf,
        Entry<BigInt>* rbuf, int count, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Scan [Entry<BigInt>]"))
    PackedScan( sbuf, rbuf, count, op, comm );
}

template<>
void Scan
( const BigFloat* sbuf, BigFloat* rbuf, int count, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Scan [BigFloat]"))
    PackedScan( sbuf, rbuf, count, op, comm );
}
template<>
void Scan
( const ValueInt<BigFloat>* sbuf,
        ValueInt<BigFloat>* rbuf, int count, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Scan [ValueInt<BigFloat>]"))
    PackedScan( sbuf, rbuf, count, op, comm );
}
template<>
void Scan
( const Entry<BigFloat>* sbuf,
        Entry<BigFloat>* rbuf, int count, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Scan [Entry<BigFloat>]"))
    PackedScan( sbuf, rbuf, count, op, comm );
}
#endif

template<typename Real>
void Scan
( const Complex<Real>* sbuf, 
        Complex<Real>* rbuf, int count, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Scan"))
    if( count != 0 )
    {
        MPI_Op opC;
        if( op == SUM )
            opC = SumOp<Complex<Real>>().op; 
        else
            opC = op.op;

#ifdef EL_AVOID_COMPLEX_MPI
        if( op == SUM )
        {
            SafeMpi
            ( MPI_Scan
              ( const_cast<Complex<Real>*>(sbuf),
                rbuf, 2*count, TypeMap<Real>(), opC, comm.comm ) );
        }
        else
        {
            SafeMpi
            ( MPI_Scan
              ( const_cast<Complex<Real>*>(sbuf),
                rbuf, count, TypeMap<Complex<Real>>(), opC, comm.comm ) );
        }
#else
        SafeMpi
        ( MPI_Scan
          ( const_cast<Complex<Real>*>(sbuf), 
            rbuf, count, TypeMap<Complex<Real>>(), opC, comm.comm ) );
#endif
    }
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

template<typename Real>
void Scan( Real* buf, int count, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Scan"))
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

#ifdef EL_HAVE_MPC
template<typename T>
void PackedScan( T* buf, int count, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::PackedScan"))
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

template<>
void Scan( BigInt* buf, int count, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Scan [BigInt]"))
    PackedScan( buf, count, op, comm );
}
template<>
void Scan( ValueInt<BigInt>* buf, int count, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Scan [ValueInt<BigInt>]"))
    PackedScan( buf, count, op, comm );
}
template<>
void Scan( Entry<BigInt>* buf, int count, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Scan [Entry<BigInt>]"))
    PackedScan( buf, count, op, comm );
}

template<>
void Scan( BigFloat* buf, int count, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Scan [BigFloat]"))
    PackedScan( buf, count, op, comm );
}
template<>
void Scan( ValueInt<BigFloat>* buf, int count, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Scan [ValueInt<BigFloat>]"))
    PackedScan( buf, count, op, comm );
}
template<>
void Scan( Entry<BigFloat>* buf, int count, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Scan [Entry<BigFloat>]"))
    PackedScan( buf, count, op, comm );
}
#endif

template<typename Real>
void Scan( Complex<Real>* buf, int count, Op op, Comm comm )
EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(CSE cse("mpi::Scan"))
    if( count != 0 )
    {
        MPI_Op opC;
        if( op == SUM )
            opC = SumOp<Complex<Real>>().op; 
        else
            opC = op.op;

#ifdef EL_AVOID_COMPLEX_MPI
        if( op == SUM )
        {
            SafeMpi
            ( MPI_Scan
              ( MPI_IN_PLACE, buf, 2*count, TypeMap<Real>(), opC, comm.comm ) );
        }
        else
        {
            SafeMpi
            ( MPI_Scan
              ( MPI_IN_PLACE, buf, count, TypeMap<Complex<Real>>(), opC, 
                comm.comm ) );
        }
#else
        SafeMpi
        ( MPI_Scan
          ( MPI_IN_PLACE, buf, count, TypeMap<Complex<Real>>(), opC, 
            comm.comm ) );
#endif
    }
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
    DEBUG_ONLY(
      CSE cse("mpi::SparseAllToAll");
      VerifySendsAndRecvs( sendCounts, recvCounts, comm );
    )
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
MPI_PROTO(ValueInt<DoubleDouble>)
MPI_PROTO(ValueInt<QuadDouble>)
MPI_PROTO(Entry<DoubleDouble>)
MPI_PROTO(Entry<QuadDouble>)
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
// TODO: MPI_PROTO(Complex<BigFloat>)
MPI_PROTO(ValueInt<BigInt>)
MPI_PROTO(ValueInt<BigFloat>)
// TODO: MPI_PROTO(ValueInt<Complex<BigFloat>>)
MPI_PROTO(Entry<BigInt>)
MPI_PROTO(Entry<BigFloat>)
// TODO: MPI_PROTO(Entry<Complex<BigFloat>>)
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
#include "El/macros/Instantiate.h"

} // namespace mpi
} // namespace El
