/*
   Copyright (c) 2009-2015, Jack Poulson
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
SafeMpi( int mpiError )
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

bool CommSameSizeAsInteger()
{ return sizeof(MPI_Comm) == sizeof(int); }

bool GroupSameSizeAsInteger()
{ return sizeof(MPI_Group) == sizeof(int); }

// MPI environmental routines
// ==========================

void Initialize( int& argc, char**& argv )
{ MPI_Init( &argc, &argv ); }

int InitializeThread( int& argc, char**& argv, int required )
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
#ifdef EL_HAVE_MPI_QUERY_THREAD
    MPI_Query_thread( &provided );
#else
    provided = 0; // equivalent to MPI_THREAD_SINGLE
#endif
    return provided;
}

void Abort( Comm comm, int errCode )
{ MPI_Abort( comm.comm, errCode ); }

double Time()
{ return MPI_Wtime(); }

void Create( UserFunction* func, bool commutes, Op& op )
{
    DEBUG_ONLY(CSE cse("mpi::Create"))
    SafeMpi( MPI_Op_create( func, commutes, &op.op ) );
}

void Free( Op& op )
{
    DEBUG_ONLY(CSE cse("mpi::Free"))
    SafeMpi( MPI_Op_free( &op.op ) );
}

void Free( Datatype& type )
{
    DEBUG_ONLY(CSE cse("mpi::Free"))
    SafeMpi( MPI_Type_free( &type ) );
}

// Communicator manipulation 
// =========================

int WorldRank()
{
    DEBUG_ONLY(CSE cse("mpi::WorldRank"))
    return Rank( mpi::COMM_WORLD ); 
}

int WorldSize()
{
    DEBUG_ONLY(CSE cse("mpi::WorldSize"))
    return Size( mpi::COMM_WORLD );
}

int Rank( Comm comm )
{
    DEBUG_ONLY(CSE cse("mpi::Rank"))
    if( comm != COMM_NULL )
    {
        int rank;
        SafeMpi( MPI_Comm_rank( comm.comm, &rank ) );
        return rank;
    }
    else return mpi::UNDEFINED;
}

int Size( Comm comm )
{
    DEBUG_ONLY(CSE cse("mpi::Size"))
    if( comm != COMM_NULL )
    {
        int size;
        SafeMpi( MPI_Comm_size( comm.comm, &size ) );
        return size;
    } 
    else return mpi::UNDEFINED;
}

void Create( Comm parentComm, Group subsetGroup, Comm& subsetComm )
{
    DEBUG_ONLY(CSE cse("mpi::Create"))
    SafeMpi( 
        MPI_Comm_create( parentComm.comm, subsetGroup.group, &subsetComm.comm ) 
    );
}

void Dup( Comm original, Comm& duplicate )
{
    DEBUG_ONLY(CSE cse("mpi::Dup"))
    SafeMpi( MPI_Comm_dup( original.comm, &duplicate.comm ) );
}

void Split( Comm comm, int color, int key, Comm& newComm )
{
    DEBUG_ONLY(CSE cse("mpi::Split"))
    SafeMpi( MPI_Comm_split( comm.comm, color, key, &newComm.comm ) );
}

void Free( Comm& comm )
{
    DEBUG_ONLY(CSE cse("mpi::Free"))
    SafeMpi( MPI_Comm_free( &comm.comm ) );
}

bool Congruent( Comm comm1, Comm comm2 )
{
    DEBUG_ONLY(CSE cse("mpi::Congruent"))
    int result;
    SafeMpi( MPI_Comm_compare( comm1.comm, comm2.comm, &result ) );
    return ( result == MPI_IDENT || result == MPI_CONGRUENT );
}

void ErrorHandlerSet( Comm comm, ErrorHandler errorHandler )
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
  bool reorder, Comm& cartComm )
{
    DEBUG_ONLY(CSE cse("mpi::CartCreate"))
    SafeMpi
    ( MPI_Cart_create
      ( comm.comm, numDims, const_cast<int*>(dimensions), 
        const_cast<int*>(periods), reorder, &cartComm.comm ) );
}

void CartSub( Comm comm, const int* remainingDims, Comm& subComm )
{
    DEBUG_ONLY(CSE cse("mpi::CartSub"))
    SafeMpi( 
        MPI_Cart_sub
        ( comm.comm, const_cast<int*>(remainingDims), &subComm.comm ) 
    );
}

// Group manipulation 
// ==================

int Rank( Group group )
{
    DEBUG_ONLY(CSE cse("mpi::Rank"))
    int rank;
    SafeMpi( MPI_Group_rank( group.group, &rank ) );
    return rank;
}

int Size( Group group )
{
    DEBUG_ONLY(CSE cse("mpi::Size"))
    int size;
    SafeMpi( MPI_Group_size( group.group, &size ) );
    return size;
}

void CommGroup( Comm comm, Group& group )
{
    DEBUG_ONLY(CSE cse("mpi::CommGroup"))
    SafeMpi( MPI_Comm_group( comm.comm, &group.group ) );
}

void Dup( Group group, Group& newGroup )
{
    DEBUG_ONLY(CSE cse("mpi::Dup"))
    // For some reason, MPI_Group_dup does not exist
    Excl( group, 0, 0, newGroup ); 
}

void Union( Group groupA, Group groupB, Group& newGroup )
{
    DEBUG_ONLY(CSE cse("mpi::Union"))
    SafeMpi( MPI_Group_union( groupA.group, groupB.group, &newGroup.group ) );
}

void Incl( Group group, int n, const int* ranks, Group& subGroup )
{
    DEBUG_ONLY(CSE cse("mpi::Incl"))
    SafeMpi( 
        MPI_Group_incl
        ( group.group, n, const_cast<int*>(ranks), &subGroup.group ) 
    );
}

void Excl( Group group, int n, const int* ranks, Group& subGroup )
{
    DEBUG_ONLY(CSE cse("mpi::Excl"))
    SafeMpi(
        MPI_Group_excl
        ( group.group, n, const_cast<int*>(ranks), &subGroup.group )
    );
}

void Difference( Group parent, Group subset, Group& complement )
{
    DEBUG_ONLY(CSE cse("mpi::Difference"))
    SafeMpi( 
        MPI_Group_difference( parent.group, subset.group, &complement.group ) 
    );
}

void Free( Group& group )
{
    DEBUG_ONLY(CSE cse("mpi::Free"))
    SafeMpi( MPI_Group_free( &group.group ) );
}

// Rank translations
// =================

int Translate( Group origGroup, int origRank, Group newGroup )
{
    DEBUG_ONLY(CSE cse("mpi::Translate"))
    int newRank;
    Translate( origGroup, 1, &origRank, newGroup, &newRank );
    return newRank;
}

int Translate( Comm origComm, int origRank, Group newGroup )
{
    DEBUG_ONLY(CSE cse("mpi::Translate"))
    int newRank;
    Translate( origComm, 1, &origRank, newGroup, &newRank );
    return newRank;
}

int Translate( Group origGroup, int origRank, Comm newComm )
{
    DEBUG_ONLY(CSE cse("mpi::Translate"))
    int newRank;
    Translate( origGroup, 1, &origRank, newComm, &newRank );
    return newRank;
}

int Translate( Comm origComm, int origRank, Comm newComm )
{
    DEBUG_ONLY(CSE cse("mpi::Translate"))
    int newRank;
    Translate( origComm, 1, &origRank, newComm, &newRank );
    return newRank;
}

void Translate
( Group origGroup, int size, const int* origRanks, 
  Group newGroup,                  int* newRanks )
{
    DEBUG_ONLY(CSE cse("mpi::Translate"))
    SafeMpi
    ( MPI_Group_translate_ranks
      ( origGroup.group, size, const_cast<int*>(origRanks), 
        newGroup.group, newRanks ) );
}

void Translate
( Comm origComm,  int size, const int* origRanks, 
  Group newGroup,                 int* newRanks )
{
    DEBUG_ONLY(CSE cse("mpi::Translate"))
    Group origGroup;
    CommGroup( origComm, origGroup );
    Translate( origGroup, size, origRanks, newGroup, newRanks );
    Free( origGroup );
}

void Translate
( Group origGroup,  int size, const int* origRanks, 
  Comm newComm,                     int* newRanks )
{
    DEBUG_ONLY(CSE cse("mpi::Translate"))
    Group newGroup;
    CommGroup( newComm,  newGroup  );
    Translate( origGroup, size, origRanks, newGroup, newRanks );
    Free( newGroup  );
}

void Translate
( Comm origComm,  int size, const int* origRanks, 
  Comm newComm,                   int* newRanks )
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
void Barrier( Comm comm )
{
    DEBUG_ONLY(CSE cse("mpi::Barrier"))
    SafeMpi( MPI_Barrier( comm.comm ) );
}

// Test for completion
bool Test( Request& request )
{
    DEBUG_ONLY(CSE cse("mpi::Test"))
    Status status;
    int flag;
    SafeMpi( MPI_Test( &request, &flag, &status ) );
    return flag;
}

// Ensure that the request finishes before continuing
void Wait( Request& request )
{
    DEBUG_ONLY(CSE cse("mpi::Wait"))
    Status status;
    SafeMpi( MPI_Wait( &request, &status ) );
}

// Ensure that the request finishes before continuing
void Wait( Request& request, Status& status )
{
    DEBUG_ONLY(CSE cse("mpi::Wait"))
    SafeMpi( MPI_Wait( &request, &status ) );
}

// Ensure that several requests finish before continuing
void WaitAll( int numRequests, Request* requests )
{
    DEBUG_ONLY(CSE cse("mpi::WaitAll"))
    vector<Status> statuses( numRequests );
    SafeMpi( MPI_Waitall( numRequests, requests, statuses.data() ) );
}

// Ensure that several requests finish before continuing
void WaitAll( int numRequests, Request* requests, Status* statuses )
{
    DEBUG_ONLY(CSE cse("mpi::WaitAll"))
    SafeMpi( MPI_Waitall( numRequests, requests, statuses ) );
}

// Nonblocking test for message completion
bool IProbe( int source, int tag, Comm comm, Status& status )
{
    DEBUG_ONLY(CSE cse("mpi::IProbe"))
    int flag;
    SafeMpi( MPI_Iprobe( source, tag, comm.comm, &flag, &status ) );
    return flag;
}
bool IProbe( int source, Comm comm, Status& status )
{ return IProbe( source, 0, comm, status ); }

template<typename T>
int GetCount( Status& status )
{
    DEBUG_ONLY(CSE cse("mpi::GetCount"))
    int count;
    SafeMpi( MPI_Get_count( &status, TypeMap<T>(), &count ) );
    return count;
}

template<typename Real>
void TaggedSend( const Real* buf, int count, int to, int tag, Comm comm )
{ 
    DEBUG_ONLY(CSE cse("mpi::Send"))
    SafeMpi( 
        MPI_Send
        ( const_cast<Real*>(buf), count, TypeMap<Real>(), to, tag, comm.comm )
    );
}

template<typename Real>
void TaggedSend
( const Complex<Real>* buf, int count, int to, int tag, Comm comm )
{
    DEBUG_ONLY(CSE cse("mpi::Send"))
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
void Send( const T* buf, int count, int to, Comm comm )
{ TaggedSend( buf, count, to, 0, comm ); }

template<typename T>
void TaggedSend( T b, int to, int tag, Comm comm )
{ TaggedSend( &b, 1, to, tag, comm ); }

template<typename T>
void Send( T b, int to, Comm comm )
{ TaggedSend( b, to, 0, comm ); }

template<typename Real>
void TaggedISend
( const Real* buf, int count, int to, int tag, Comm comm, Request& request )
{ 
    DEBUG_ONLY(CSE cse("mpi::ISend"))
    SafeMpi
    ( MPI_Isend
      ( const_cast<Real*>(buf), count, TypeMap<Real>(), to, 
        tag, comm.comm, &request ) );
}

template<typename Real>
void TaggedISend
( const Complex<Real>* buf, int count, int to, int tag, Comm comm, 
  Request& request )
{
    DEBUG_ONLY(CSE cse("mpi::ISend"))
#ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi
    ( MPI_Isend
      ( const_cast<Complex<Real>*>(buf), 2*count, 
        TypeMap<Real>(), to, tag, comm.comm, &request ) );
#else
    SafeMpi
    ( MPI_Isend
      ( const_cast<Complex<Real>*>(buf), count, 
        TypeMap<Complex<Real>>(), to, tag, comm.comm, &request ) );
#endif
}

template<typename T>
void ISend
( const T* buf, int count, int to, Comm comm, Request& request )
{ TaggedISend( buf, count, to, 0, comm, request ); } 

template<typename T>
void TaggedISend( T b, int to, int tag, Comm comm, Request& request )
{ TaggedISend( &b, 1, to, tag, comm, request ); }

template<typename T>
void ISend( T b, int to, Comm comm, Request& request )
{ TaggedISend( b, to, 0, comm, request ); }

template<typename Real>
void TaggedISSend
( const Real* buf, int count, int to, int tag, Comm comm, Request& request )
{
    DEBUG_ONLY(CSE cse("mpi::ISSend"))
    SafeMpi
    ( MPI_Issend
      ( const_cast<Real*>(buf), count, TypeMap<Real>(), to, 
        tag, comm.comm, &request ) );
}

template<typename Real>
void TaggedISSend
( const Complex<Real>* buf, int count, int to, int tag, Comm comm, 
  Request& request )
{
    DEBUG_ONLY(CSE cse("mpi::ISSend"))
#ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi
    ( MPI_Issend
      ( const_cast<Complex<Real>*>(buf), 2*count, 
        TypeMap<Real>(), to, tag, comm.comm, &request ) );
#else
    SafeMpi
    ( MPI_Issend
      ( const_cast<Complex<Real>*>(buf), count, 
        TypeMap<Complex<Real>>(), to, tag, comm.comm, &request ) );
#endif
}

template<typename T>
void ISSend( const T* buf, int count, int to, Comm comm, Request& request )
{ TaggedISSend( buf, count, to, 0, comm, request ); }

template<typename T>
void TaggedISSend( T b, int to, int tag, Comm comm, Request& request )
{ TaggedISSend( &b, 1, to, tag, comm, request ); }

template<typename Real>
void TaggedRecv( Real* buf, int count, int from, int tag, Comm comm )
{
    DEBUG_ONLY(CSE cse("mpi::Recv"))
    Status status;
    SafeMpi
    ( MPI_Recv( buf, count, TypeMap<Real>(), from, tag, comm.comm, &status ) );
}

template<typename Real>
void TaggedRecv( Complex<Real>* buf, int count, int from, int tag, Comm comm )
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
{ TaggedRecv( buf, count, from, mpi::ANY_TAG, comm ); }

template<typename T>
T TaggedRecv( int from, int tag, Comm comm )
{ T b; TaggedRecv( &b, 1, from, tag, comm ); return b; }

template<typename T>
T Recv( int from, Comm comm )
{ return TaggedRecv<T>( from, mpi::ANY_TAG, comm ); }

template<typename Real>
void TaggedIRecv
( Real* buf, int count, int from, int tag, Comm comm, Request& request )
{
    DEBUG_ONLY(CSE cse("mpi::IRecv"))
    SafeMpi
    ( MPI_Irecv
      ( buf, count, TypeMap<Real>(), from, tag, comm.comm, &request ) );
}

template<typename Real>
void TaggedIRecv
( Complex<Real>* buf, int count, int from, int tag, 
  Comm comm, Request& request )
{
    DEBUG_ONLY(CSE cse("mpi::IRecv"))
#ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi
    ( MPI_Irecv( buf, 2*count, TypeMap<Real>(), from, tag, comm.comm, &request ) );
#else
    SafeMpi
    ( MPI_Irecv
      ( buf, count, TypeMap<Complex<Real>>(), from, tag, comm.comm, &request ) );
#endif
}

template<typename T>
void IRecv( T* buf, int count, int from, Comm comm, Request& request )
{ TaggedIRecv( buf, count, from, mpi::ANY_TAG, comm, request ); }

template<typename T>
T TaggedIRecv( int from, int tag, Comm comm, Request& request )
{ T b; TaggedIRecv( &b, 1, from, tag, comm, request ); return b; }

template<typename T>
T IRecv( int from, Comm comm, Request& request )
{ return TaggedIRecv<T>( from, mpi::ANY_TAG, comm, request ); }

template<typename Real>
void TaggedSendRecv
( const Real* sbuf, int sc, int to,   int stag,
        Real* rbuf, int rc, int from, int rtag, Comm comm )
{
    DEBUG_ONLY(CSE cse("mpi::SendRecv"))
    Status status;
    SafeMpi
    ( MPI_Sendrecv
      ( const_cast<Real*>(sbuf), sc, TypeMap<Real>(), to,   stag,
        rbuf,                    rc, TypeMap<Real>(), from, rtag, 
        comm.comm, &status ) );
}

template<typename Real>
void TaggedSendRecv
( const Complex<Real>* sbuf, int sc, int to,   int stag,
        Complex<Real>* rbuf, int rc, int from, int rtag, Comm comm )
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
{ TaggedSendRecv( sbuf, sc, to, 0, rbuf, rc, from, mpi::ANY_TAG, comm ); }

template<typename T>
T TaggedSendRecv( T sb, int to, int stag, int from, int rtag, Comm comm )
{ 
    T rb; 
    TaggedSendRecv( &sb, 1, to, stag, &rb, 1, from, rtag, comm ); 
    return rb; 
}

template<typename T>
T SendRecv( T sb, int to, int from, Comm comm )
{ return TaggedSendRecv( sb, to, 0, from, mpi::ANY_TAG, comm ); }

template<typename Real>
void TaggedSendRecv
( Real* buf, int count, int to, int stag, int from, int rtag, Comm comm )
{
    DEBUG_ONLY(CSE cse("mpi::SendRecv"))
    Status status;
    SafeMpi
    ( MPI_Sendrecv_replace
      ( buf, count, TypeMap<Real>(), to, stag, from, rtag, comm.comm, &status ) );
}

template<typename Real>
void TaggedSendRecv
( Complex<Real>* buf, int count, int to, int stag, int from, int rtag, Comm comm )
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
{ TaggedSendRecv( buf, count, to, 0, from, mpi::ANY_TAG, comm ); }

template<typename Real>
void Broadcast( Real* buf, int count, int root, Comm comm )
{
    DEBUG_ONLY(CSE cse("mpi::Broadcast"))
    SafeMpi( MPI_Bcast( buf, count, TypeMap<Real>(), root, comm.comm ) );
}

template<typename Real>
void Broadcast( Complex<Real>* buf, int count, int root, Comm comm )
{
    DEBUG_ONLY(CSE cse("mpi::Broadcast"))
#ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi( MPI_Bcast( buf, 2*count, TypeMap<Real>(), root, comm.comm ) );
#else
    SafeMpi( MPI_Bcast( buf, count, TypeMap<Complex<Real>>(), root, comm.comm ) );
#endif
}

template<typename T>
void Broadcast( T& b, int root, Comm comm )
{ Broadcast( &b, 1, root, comm ); }

template<typename Real>
void IBroadcast( Real* buf, int count, int root, Comm comm, Request& request )
{
    DEBUG_ONLY(CSE cse("mpi::IBroadcast"))
#ifdef EL_HAVE_NONBLOCKING_COLLECTIVES
    SafeMpi
    ( MPI_Ibcast( buf, count, TypeMap<Real>(), root, comm.comm, &request ) );
#else
    LogicError("Elemental was not configured with non-blocking support");
#endif
}

template<typename Real>
void IBroadcast
( Complex<Real>* buf, int count, int root, Comm comm, Request& request )
{
    DEBUG_ONLY(CSE cse("mpi::IBroadcast"))
#ifdef EL_HAVE_NONBLOCKING_COLLECTIVES
#ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi
    ( MPI_Ibcast( buf, 2*count, TypeMap<Real>(), root, comm.comm, &request ) );
#else
    SafeMpi
    ( MPI_Ibcast
      ( buf, count, TypeMap<Complex<Real>>(), root, comm.comm, &request ) );
#endif
#else
    LogicError("Elemental was not configured with non-blocking support");
#endif
}

template<typename T>
void IBroadcast( T& b, int root, Comm comm, Request& request )
{ IBroadcast( &b, 1, root, comm, request ); }

template<typename Real>
void Gather
( const Real* sbuf, int sc,
        Real* rbuf, int rc, int root, Comm comm )
{
    DEBUG_ONLY(CSE cse("mpi::Gather"))
    SafeMpi
    ( MPI_Gather
      ( const_cast<Real*>(sbuf), sc, TypeMap<Real>(),
        rbuf,                    rc, TypeMap<Real>(), root, comm.comm ) );
}

template<typename Real>
void Gather
( const Complex<Real>* sbuf, int sc,
        Complex<Real>* rbuf, int rc, int root, Comm comm )
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
        Real* rbuf, int rc, int root, Comm comm, Request& request )
{
    DEBUG_ONLY(CSE cse("mpi::IGather"))
#ifdef EL_HAVE_NONBLOCKING_COLLECTIVES
    SafeMpi
    ( MPI_Igather
      ( const_cast<Real*>(sbuf), sc, TypeMap<Real>(),
        rbuf,                    rc, TypeMap<Real>(), root, comm.comm, &request ) );
#else
    LogicError("Elemental was not configured with non-blocking support");
#endif
}

template<typename Real>
void IGather
( const Complex<Real>* sbuf, int sc,
        Complex<Real>* rbuf, int rc, int root, Comm comm, Request& request )
{
    DEBUG_ONLY(CSE cse("mpi::IGather"))
#ifdef EL_HAVE_NONBLOCKING_COLLECTIVES
#ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi
    ( MPI_Igather
      ( const_cast<Complex<Real>*>(sbuf), 2*sc, TypeMap<Real>(),
        rbuf,                             2*rc, TypeMap<Real>(), 
        root, comm.comm, &request ) );
#else
    SafeMpi
    ( MPI_Igather
      ( const_cast<Complex<Real>*>(sbuf), sc, TypeMap<Complex<Real>>(),
        rbuf,                             rc, TypeMap<Complex<Real>>(), 
        root, comm.comm, &request ) );
#endif
#else
    LogicError("Elemental was not configured with non-blocking support");
#endif
}

template<typename Real>
void Gather
( const Real* sbuf, int sc,
        Real* rbuf, const int* rcs, const int* rds, int root, Comm comm )
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

template<typename Real>
void Gather
( const Complex<Real>* sbuf, int sc,
        Complex<Real>* rbuf, const int* rcs, const int* rds, int root, Comm comm )
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

template<typename Real>
void AllGather
( const Complex<Real>* sbuf, int sc,
        Complex<Real>* rbuf, int rc, Comm comm )
{
    DEBUG_ONLY(CSE cse("mpi::AllGather"))
#ifdef EL_USE_BYTE_ALLGATHERS
    SafeMpi
    ( MPI_Allgather
      ( (UCP)const_cast<Complex<Real>*>(sbuf), 2*sizeof(Real)*sc, MPI_UNSIGNED_CHAR, 
        (UCP)rbuf,                             2*sizeof(Real)*rc, MPI_UNSIGNED_CHAR, 
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
        rbuf,                             rc, TypeMap<Complex<Real>>(), comm.comm ) );
 #endif
#endif
}

template<typename Real>
void AllGather
( const Real* sbuf, int sc,
        Real* rbuf, const int* rcs, const int* rds, Comm comm )
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

template<typename Real>
void AllGather
( const Complex<Real>* sbuf, int sc,
        Complex<Real>* rbuf, const int* rcs, const int* rds, Comm comm )
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
      ( (UCP)const_cast<Complex<Real>*>(sbuf), 2*sizeof(Real)*sc, MPI_UNSIGNED_CHAR, 
        (UCP)rbuf, byteRcs.data(), byteRds.data(),                MPI_UNSIGNED_CHAR, 
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
{
    DEBUG_ONLY(CSE cse("mpi::Scatter"))
    SafeMpi
    ( MPI_Scatter
      ( const_cast<Real*>(sbuf), sc, TypeMap<Real>(),
        rbuf,                    rc, TypeMap<Real>(), root, comm.comm ) );
}

template<typename Real>
void Scatter
( const Complex<Real>* sbuf, int sc,
        Complex<Real>* rbuf, int rc, int root, Comm comm )
{
    DEBUG_ONLY(CSE cse("mpi::Scatter"))
#ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi
    ( MPI_Scatter
      ( const_cast<Complex<Real>*>(sbuf), 2*sc, TypeMap<Real>(),
        rbuf,                             2*rc, TypeMap<Real>(), root, comm.comm ) );
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
{
    DEBUG_ONLY(CSE cse("mpi::Scatter"))
    const int commRank = Rank( comm );
    if( commRank == root )
    {
#ifdef EL_HAVE_MPI_IN_PLACE
        SafeMpi
        ( MPI_Scatter
          ( buf,          sc, TypeMap<Real>(), 
            MPI_IN_PLACE, rc, TypeMap<Real>(), root, comm.comm ) );
#else
        const int commSize = Size( comm );
        vector<Real> sendBuf( sc*commSize );
        MemCopy( sendBuf.data(), buf, sc*commSize );
        SafeMpi
        ( MPI_Scatter
          ( sendBuf.data(), sc, TypeMap<Real>(), 
            buf,            rc, TypeMap<Real>(), root, comm.comm ) );
#endif
    }
    else
    {
        SafeMpi
        ( MPI_Scatter
          ( 0,   sc, TypeMap<Real>(), 
            buf, rc, TypeMap<Real>(), root, comm.comm ) );
    }
}

template<typename Real>
void Scatter( Complex<Real>* buf, int sc, int rc, int root, Comm comm )
{
    DEBUG_ONLY(CSE cse("mpi::Scatter"))
    const int commRank = Rank( comm );
    if( commRank == root )
    {
#ifdef EL_AVOID_COMPLEX_MPI
# ifdef EL_HAVE_MPI_IN_PLACE
        SafeMpi
        ( MPI_Scatter
          ( buf,          2*sc, TypeMap<Real>(), 
            MPI_IN_PLACE, 2*rc, TypeMap<Real>(), root, comm.comm ) );
# else
        const int commSize = Size( comm );
        vector<Complex<Real>> sendBuf( sc*commSize );
        MemCopy( sendBuf.data(), buf, sc*commSize );
        SafeMpi
        ( MPI_Scatter
          ( sendBuf.data(), 2*sc, TypeMap<Real>(),          
            buf,            2*rc, TypeMap<Real>(), root, comm.comm ) );
# endif
#else
# ifdef EL_HAVE_MPI_IN_PLACE
        SafeMpi
        ( MPI_Scatter
          ( buf,          sc, TypeMap<Complex<Real>>(), 
            MPI_IN_PLACE, rc, TypeMap<Complex<Real>>(), root, comm.comm ) );
# else
        const int commSize = Size( comm );
        vector<Complex<Real>> sendBuf( sc*commSize );
        MemCopy( sendBuf.data(), buf, sc*commSize );
        SafeMpi
        ( MPI_Scatter
          ( sendBuf.data(), sc, TypeMap<Complex<Real>>(),
            buf,            rc, TypeMap<Complex<Real>>(), root, comm.comm ) );
# endif
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
{
    DEBUG_ONLY(CSE cse("mpi::AllToAll"))
    SafeMpi
    ( MPI_Alltoall
      ( const_cast<Real*>(sbuf), sc, TypeMap<Real>(),
        rbuf,                    rc, TypeMap<Real>(), comm.comm ) );
}

template<typename Real>
void AllToAll
( const Complex<Real>* sbuf, int sc,
        Complex<Real>* rbuf, int rc, Comm comm )
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

template<typename Real>
void AllToAll
( const Complex<Real>* sbuf, const int* scs, const int* sds,
        Complex<Real>* rbuf, const int* rcs, const int* rds, Comm comm )
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
std::vector<T> AllToAll
( const std::vector<T>& sendBuf,
  const std::vector<int>& sendCounts, 
  const std::vector<int>& sendOffs,
  mpi::Comm comm )
{
    const int commSize = mpi::Size( comm ); 
    std::vector<int> recvCounts(commSize);
    mpi::AllToAll( sendCounts.data(), 1, recvCounts.data(), 1, comm ); 
    std::vector<int> recvOffs;
    const int totalRecv = El::Scan( recvCounts, recvOffs );
    std::vector<T> recvBuf(totalRecv);
    mpi::AllToAll
    ( sendBuf.data(), sendCounts.data(), sendOffs.data(),
      recvBuf.data(), recvCounts.data(), recvOffs.data(), comm );
    return recvBuf;
}

template<typename Real>
void Reduce
( const Real* sbuf, Real* rbuf, int count, Op op, int root, Comm comm )
{
    DEBUG_ONLY(CSE cse("mpi::Reduce"))
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
        ( MPI_Reduce
          ( const_cast<Real*>(sbuf), rbuf, count, TypeMap<Real>(),
            opC, root, comm.comm ) );
    }
}

template<typename Real>
void Reduce
( const Complex<Real>* sbuf, 
        Complex<Real>* rbuf, int count, Op op, int root, Comm comm )
{
    DEBUG_ONLY(CSE cse("mpi::Reduce"))
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
}

template<typename T>
void Reduce( const T* sbuf, T* rbuf, int count, int root, Comm comm )
{ Reduce( sbuf, rbuf, count, mpi::SUM, root, comm ); }

template<typename T>
T Reduce( T sb, Op op, int root, Comm comm )
{ 
    T rb;
    Reduce( &sb, &rb, 1, op, root, comm );
    return rb;
}

template<typename T>
T Reduce( T sb, int root, Comm comm )
{ 
    T rb;
    Reduce( &sb, &rb, 1, mpi::SUM, root, comm );
    return rb;
}

template<typename Real>
void Reduce( Real* buf, int count, Op op, int root, Comm comm )
{
    DEBUG_ONLY(CSE cse("mpi::Reduce"))
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

        const int commRank = Rank( comm );
        if( commRank == root )
        {
#ifdef EL_HAVE_MPI_IN_PLACE
            SafeMpi
            ( MPI_Reduce
              ( MPI_IN_PLACE, buf, count, TypeMap<Real>(), opC, root, 
                comm.comm ) );
#else
            vector<Real> sendBuf( count );
            MemCopy( sendBuf.data(), buf, count );
            SafeMpi
            ( MPI_Reduce
              ( sendBuf.data(), buf, count, TypeMap<Real>(), opC, root, 
                comm.comm ) );
#endif
        }
        else
            SafeMpi
            ( MPI_Reduce
              ( buf, 0, count, TypeMap<Real>(), opC, root, comm.comm ) );
    }
}

template<typename Real>
void Reduce( Complex<Real>* buf, int count, Op op, int root, Comm comm )
{
    DEBUG_ONLY(CSE cse("mpi::Reduce"))
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
# ifdef EL_HAVE_MPI_IN_PLACE
                SafeMpi
                ( MPI_Reduce
                  ( MPI_IN_PLACE, buf, 2*count, TypeMap<Real>(), opC, 
                    root, comm.comm ) );
# else
                vector<Complex<Real>> sendBuf( count );
                MemCopy( sendBuf.data(), buf, count );
                SafeMpi
                ( MPI_Reduce
                  ( sendBuf.data(), buf, 2*count, TypeMap<Real>(), opC, 
                    root, comm.comm ) );
# endif
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
# ifdef EL_HAVE_MPI_IN_PLACE
                SafeMpi
                ( MPI_Reduce
                  ( MPI_IN_PLACE, buf, count, TypeMap<Complex<Real>>(), opC, 
                    root, comm.comm ) );
# else
                vector<Complex<Real>> sendBuf( count );
                MemCopy( sendBuf.data(), buf, count );
                SafeMpi
                ( MPI_Reduce
                  ( sendBuf.data(), buf, count, TypeMap<Complex<Real>>(), opC, 
                    root, comm.comm ) );
# endif
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
# ifdef EL_HAVE_MPI_IN_PLACE
            SafeMpi
            ( MPI_Reduce
              ( MPI_IN_PLACE, buf, count, TypeMap<Complex<Real>>(), opC, 
                root, comm.comm ) );
# else
            vector<Complex<Real>> sendBuf( count );
            MemCopy( sendBuf.data(), buf, count );
            SafeMpi
            ( MPI_Reduce
              ( sendBuf.data(), buf, count, TypeMap<Complex<Real>>(), opC, 
                root, comm.comm ) );
# endif
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
{ Reduce( buf, count, mpi::SUM, root, comm ); }

template<typename Real>
void AllReduce( const Real* sbuf, Real* rbuf, int count, Op op, Comm comm )
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

template<typename Real>
void AllReduce
( const Complex<Real>* sbuf, Complex<Real>* rbuf, int count, Op op, Comm comm )
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
{ AllReduce( sbuf, rbuf, count, mpi::SUM, comm ); }

template<typename T>
T AllReduce( T sb, Op op, Comm comm )
{ T rb; AllReduce( &sb, &rb, 1, op, comm ); return rb; }

template<typename T>
T AllReduce( T sb, Comm comm )
{ return AllReduce( sb, mpi::SUM, comm ); }

template<typename Real>
void AllReduce( Real* buf, int count, Op op, Comm comm )
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

#ifdef EL_HAVE_MPI_IN_PLACE
        SafeMpi
        ( MPI_Allreduce
          ( MPI_IN_PLACE, buf, count, TypeMap<Real>(), opC, comm.comm ) );
#else
        vector<Real> sendBuf( count );
        MemCopy( sendBuf.data(), buf, count );
        SafeMpi
        ( MPI_Allreduce
          ( sendBuf.data(), buf, count, TypeMap<Real>(), opC, comm.comm ) );
#endif
    }
}

template<typename Real>
void AllReduce( Complex<Real>* buf, int count, Op op, Comm comm )
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
# ifdef EL_HAVE_MPI_IN_PLACE
            SafeMpi
            ( MPI_Allreduce
              ( MPI_IN_PLACE, buf, 2*count, TypeMap<Real>(), opC, comm.comm ) );
# else
            vector<Complex<Real>> sendBuf( count );
            MemCopy( sendBuf.data(), buf, count );
            SafeMpi
            ( MPI_Allreduce
              ( sendBuf.data(), buf, 2*count, TypeMap<Real>(), opC, 
                comm.comm ) );
# endif
        }
        else
        {
# ifdef EL_HAVE_MPI_IN_PLACE
            SafeMpi
            ( MPI_Allreduce
              ( MPI_IN_PLACE, buf, count, TypeMap<Complex<Real>>(), 
                opC, comm.comm ) );
# else
            vector<Complex<Real>> sendBuf( count );
            MemCopy( sendBuf.data(), buf, count );
            SafeMpi
            ( MPI_Allreduce
              ( sendBuf.data(), buf, count, TypeMap<Complex<Real>>(), 
                opC, comm.comm ) );
# endif
        }
#else
# ifdef EL_HAVE_MPI_IN_PLACE
        SafeMpi
        ( MPI_Allreduce
          ( MPI_IN_PLACE, buf, count, TypeMap<Complex<Real>>(), opC, 
            comm.comm ) );
# else
        vector<Complex<Real>> sendBuf( count );
        MemCopy( sendBuf.data(), buf, count );
        SafeMpi
        ( MPI_Allreduce
          ( sendBuf.data(), buf, count, TypeMap<Complex<Real>>(), opC, 
            comm.comm ) );
# endif
#endif
    }
}

template<typename T>
void AllReduce( T* buf, int count, Comm comm )
{ AllReduce( buf, count, mpi::SUM, comm ); }

template<typename Real>
void ReduceScatter( Real* sbuf, Real* rbuf, int rc, Op op, Comm comm )
{
    DEBUG_ONLY(CSE cse("mpi::ReduceScatter"))
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

template<typename Real>
void ReduceScatter
( Complex<Real>* sbuf, Complex<Real>* rbuf, int rc, Op op, Comm comm )
{
    DEBUG_ONLY(CSE cse("mpi::ReduceScatter"))
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
{ ReduceScatter( sbuf, rbuf, rc, mpi::SUM, comm ); }

template<typename T>
T ReduceScatter( T sb, Op op, Comm comm )
{ T rb; ReduceScatter( &sb, &rb, 1, op, comm ); return rb; }

template<typename T>
T ReduceScatter( T sb, Comm comm )
{ return ReduceScatter( sb, mpi::SUM, comm ); }

template<typename Real>
void ReduceScatter( Real* buf, int rc, Op op, Comm comm )
{
    DEBUG_ONLY(CSE cse("mpi::ReduceScatter"))
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
# ifdef EL_HAVE_MPI_IN_PLACE
    SafeMpi
    ( MPI_Reduce_scatter_block
      ( MPI_IN_PLACE, buf, rc, TypeMap<Real>(), opC, comm.comm ) );
# else
    const int commSize = Size( comm );
    vector<Real> sendBuf( rc*commSize );
    MemCopy( sendBuf.data(), buf, rc*commSize );
    SafeMpi
    ( MPI_Reduce_scatter_block
      ( sendBuf.data(), buf, rc, TypeMap<Real>(), opC, comm.comm ) );
# endif
#else
    const int commSize = Size( comm );
    Reduce( buf, rc*commSize, op, 0, comm );
    Scatter( buf, rc, rc, 0, comm );
#endif
}

// TODO: Handle case where op is not summation
template<typename Real>
void ReduceScatter( Complex<Real>* buf, int rc, Op op, Comm comm )
{
    DEBUG_ONLY(CSE cse("mpi::ReduceScatter"))
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
#  ifdef EL_HAVE_MPI_IN_PLACE
    SafeMpi
    ( MPI_Reduce_scatter_block
      ( MPI_IN_PLACE, buf, 2*rc, TypeMap<Real>(), opC, comm.comm ) );
#  else 
    const int commSize = Size( comm );
    vector<Complex<Real>> sendBuf( rc*commSize );
    MemCopy( sendBuf.data(), buf, rc*commSize );
    SafeMpi
    ( MPI_Reduce_scatter_block
      ( sendBuf.data(), buf, 2*rc, TypeMap<Real>(), opC, comm.comm ) );
#  endif
# else
#  ifdef EL_HAVE_MPI_IN_PLACE
    SafeMpi
    ( MPI_Reduce_scatter_block
      ( MPI_IN_PLACE, buf, rc, TypeMap<Complex<Real>>(), opC, comm.comm ) );
#  else
    const int commSize = Size( comm );
    vector<Complex<Real>> sendBuf( rc*commSize );
    MemCopy( sendBuf.data(), buf, rc*commSize );
    SafeMpi
    ( MPI_Reduce_scatter_block
      ( sendBuf.data(), buf, rc, TypeMap<Complex<Real>>(), opC, comm.comm ) );
#  endif
# endif
#else
    const int commSize = Size( comm );
    Reduce( buf, rc*commSize, op, 0, comm );
    Scatter( buf, rc, rc, 0, comm );
#endif
}

template<typename T>
void ReduceScatter( T* buf, int rc, Comm comm )
{ ReduceScatter( buf, rc, mpi::SUM, comm ); }

template<typename Real>
void ReduceScatter
( const Real* sbuf, Real* rbuf, const int* rcs, Op op, Comm comm )
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

template<typename Real>
void ReduceScatter
( const Complex<Real>* sbuf, Complex<Real>* rbuf, const int* rcs, Op op, Comm comm )
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
{ ReduceScatter( sbuf, rbuf, rcs, mpi::SUM, comm ); }

void VerifySendsAndRecvs
( const vector<int>& sendCounts,
  const vector<int>& recvCounts, mpi::Comm comm )
{
    DEBUG_ONLY(CSE cse("mpi::VerifySendsAndRecvs"))
    const int commSize = mpi::Size( comm );
    vector<int> actualRecvCounts(commSize);
    mpi::AllToAll
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

template<typename Real>
void Scan
( const Complex<Real>* sbuf, 
        Complex<Real>* rbuf, int count, Op op, Comm comm )
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
{ Scan( sbuf, rbuf, count, mpi::SUM, comm ); }

template<typename T>
T Scan( T sb, Op op, Comm comm )
{ 
    T rb;
    Scan( &sb, &rb, 1, op, comm );
    return rb;
}

template<typename T>
T Scan( T sb, Comm comm )
{ 
    T rb;
    Scan( &sb, &rb, 1, mpi::SUM, comm );
    return rb;
}

template<typename Real>
void Scan( Real* buf, int count, Op op, Comm comm )
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

#ifdef EL_HAVE_MPI_IN_PLACE
        SafeMpi
        ( MPI_Scan
          ( MPI_IN_PLACE, buf, count, TypeMap<Real>(), opC, comm.comm ) );
#else
        vector<Real> sendBuf( count );
        MemCopy( sendBuf.data(), buf, count );
        SafeMpi
        ( MPI_Scan
          ( sendBuf.data(), buf, count, TypeMap<Real>(), opC, comm.comm ) );
#endif
    }
}

template<typename Real>
void Scan( Complex<Real>* buf, int count, Op op, Comm comm )
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
# ifdef EL_HAVE_MPI_IN_PLACE
            SafeMpi
            ( MPI_Scan
              ( MPI_IN_PLACE, buf, 2*count, TypeMap<Real>(), opC, comm.comm ) );
# else
            vector<Complex<Real>> sendBuf( count );
            MemCopy( sendBuf.data(), buf, count );
            SafeMpi
            ( MPI_Scan
              ( sendBuf.data(), buf, 2*count, TypeMap<Real>(), opC, 
                comm.comm ) );
# endif
        }
        else
        {
# ifdef EL_HAVE_MPI_IN_PLACE
            SafeMpi
            ( MPI_Scan
              ( MPI_IN_PLACE, buf, count, TypeMap<Complex<Real>>(), opC, 
                comm.comm ) );
# else
            vector<Complex<Real>> sendBuf( count );
            MemCopy( sendBuf.data(), buf, count );
            SafeMpi
            ( MPI_Scan
              ( sendBuf.data(), buf, count, TypeMap<Complex<Real>>(), opC, 
                comm.comm ) );
# endif
        }
#else
# ifdef EL_HAVE_MPI_IN_PLACE
        SafeMpi
        ( MPI_Scan
          ( MPI_IN_PLACE, buf, count, TypeMap<Complex<Real>>(), opC, 
            comm.comm ) );
# else
        vector<Complex<Real>> sendBuf( count );
        MemCopy( sendBuf.data(), buf, count );
        SafeMpi
        ( MPI_Scan
          ( sendBuf.data(), buf, count, TypeMap<Complex<Real>>(), opC, 
            comm.comm ) );
# endif
#endif
    }
}

template<typename T>
void Scan( T* buf, int count, Comm comm )
{ Scan( buf, count, mpi::SUM, comm ); }

template<typename T>
void SparseAllToAll
( const vector<T>& sendBuffer,
  const vector<int>& sendCounts, const vector<int>& sendDispls,
        vector<T>& recvBuffer,
  const vector<int>& recvCounts, const vector<int>& recvDispls,
        mpi::Comm comm )
{
#ifdef EL_USE_CUSTOM_ALLTOALLV
    const int commSize = mpi::Size( comm );
    int numSends=0,numRecvs=0;
    for( int q=0; q<commSize; ++q )
    {
        if( sendCounts[q] != 0 )
            ++numSends;
        if( recvCounts[q] != 0 )
            ++numRecvs;
    }
    vector<mpi::Status> statuses(numSends+numRecvs);
    vector<mpi::Request> requests(numSends+numRecvs);
    int rCount=0;
    for( int q=0; q<commSize; ++q )
    {
        int count = recvCounts[q];
        int displ = recvDispls[q];
        if( count != 0 )
            mpi::IRecv
            ( &recvBuffer[displ], count, q, comm, requests[rCount++] );
    }
#ifdef EL_BARRIER_IN_ALLTOALLV
    // This should help ensure that recvs are posted before the sends
    mpi::Barrier( comm );
#endif
    for( int q=0; q<commSize; ++q )
    {
        int count = sendCounts[q];
        int displ = sendDispls[q];
        if( count != 0 )
            mpi::ISend
            ( &sendBuffer[displ], count, q, comm, requests[rCount++] );
    }
    mpi::WaitAll( numSends+numRecvs, requests.data(), statuses.data() );
#else
    mpi::AllToAll
    ( sendBuffer.data(), sendCounts.data(), sendDispls.data(),
      recvBuffer.data(), recvCounts.data(), recvDispls.data(), comm );
#endif
}

#define MPI_PROTO(T) \
  template int GetCount<T>( Status& status ); \
  template void TaggedSend( const T* buf, int count, int to, int tag, Comm comm ); \
  template void Send( const T* buf, int count, int to, Comm comm ); \
  template void TaggedSend( T b, int to, int tag, Comm comm ); \
  template void Send( T b, int to, Comm comm ); \
  template void TaggedISend( const T* buf, int count, int to, int tag, Comm comm, Request& request ); \
  template void ISend( const T* buf, int count, int to, Comm comm, Request& request ); \
  template void TaggedISend( T buf, int to, int tag, Comm comm, Request& request ); \
  template void ISend( T buf, int to, Comm comm, Request& request ); \
  template void TaggedISSend( const T* buf, int count, int to, int tag, Comm comm, Request& request ); \
  template void ISSend( const T* buf, int count, int to, Comm comm, Request& request ); \
  template void TaggedISSend( T b, int to, int tag, Comm comm, Request& request ); \
  template void TaggedRecv( T* buf, int count, int from, int tag, Comm comm ); \
  template void Recv( T* buf, int count, int from, Comm comm ); \
  template T TaggedRecv<T>( int from, int tag, Comm comm ); \
  template T Recv( int from, Comm comm ); \
  template void TaggedIRecv( T* buf, int count, int from, int tag, Comm comm, Request& request ); \
  template void IRecv( T* buf, int count, int from, Comm comm, Request& request ); \
  template T TaggedIRecv<T>( int from, int tag, Comm comm, Request& request ); \
  template T IRecv<T>( int from, Comm comm, Request& request ); \
  template void TaggedSendRecv \
  ( const T* sbuf, int sc, int to,   int stag, \
          T* rbuf, int rc, int from, int rtag, Comm comm ); \
  template void SendRecv \
  ( const T* sbuf, int sc, int to, \
          T* rbuf, int rc, int from, Comm comm ); \
  template T TaggedSendRecv \
  ( T sb, int to, int stag, int from, int rtag, Comm comm ); \
  template T SendRecv( T sb, int to, int from, Comm comm ); \
  template void TaggedSendRecv \
  ( T* buf, int count, int to, int stag, int from, int rtag, Comm comm ); \
  template void SendRecv \
  ( T* buf, int count, int to, int from, Comm comm ); \
  template void Broadcast( T* buf, int count, int root, Comm comm ); \
  template void Broadcast( T& b, int root, Comm comm ); \
  template void IBroadcast( T* buf, int count, int root, Comm comm, Request& request ); \
  template void IBroadcast( T& b, int root, Comm comm, Request& request ); \
  template void Gather( const T* sbuf, int sc, T* rbuf, int rc, int root, Comm comm ); \
  template void IGather \
  ( const T* sbuf, int sc, \
          T* rbuf, int rc, int root, Comm comm, Request& request ); \
  template void Gather \
  ( const T* sbuf, int sc, \
          T* rbuf, const int* rcs, const int* rds, int root, Comm comm ); \
  template void AllGather( const T* sbuf, int sc, T* rbuf, int rc, Comm comm ); \
  template void AllGather \
  ( const T* sbuf, int sc, \
          T* rbuf, const int* rcs, const int* rds, Comm comm ); \
  template void Scatter \
  ( const T* sbuf, int sc, \
          T* rbuf, int rc, int root, Comm comm ); \
  template void Scatter( T* buf, int sc, int rc, int root, Comm comm ); \
  template void AllToAll \
  ( const T* sbuf, int sc, \
          T* rbuf, int rc, Comm comm ); \
  template void AllToAll \
  ( const T* sbuf, const int* scs, const int* sds, \
          T* rbuf, const int* rcs, const int* rds, Comm comm ); \
  template std::vector<T> AllToAll \
  ( const std::vector<T>& sendBuf, \
    const std::vector<int>& sendCounts, \
    const std::vector<int>& sendOffs, \
    mpi::Comm comm ); \
  template void Reduce( const T* sbuf, T* rbuf, int count, Op op, int root, Comm comm ); \
  template void Reduce( const T* sbuf, T* rbuf, int count, int root, Comm comm ); \
  template T Reduce( T sb, Op op, int root, Comm comm ); \
  template T Reduce( T sb, int root, Comm comm ); \
  template void Reduce( T* buf, int count, Op op, int root, Comm comm ); \
  template void Reduce( T* buf, int count, int root, Comm comm ); \
  template void AllReduce( const T* sbuf, T* rbuf, int count, Op op, Comm comm ); \
  template void AllReduce( const T* sbuf, T* rbuf, int count, Comm comm ); \
  template T AllReduce( T sb, Op op, Comm comm ); \
  template T AllReduce( T sb, Comm comm ); \
  template void AllReduce( T* buf, int count, Op op, Comm comm ); \
  template void AllReduce( T* buf, int count, Comm comm ); \
  template void ReduceScatter( T* sbuf, T* rbuf, int rc, Op op, Comm comm ); \
  template void ReduceScatter( T* sbuf, T* rbuf, int rc, Comm comm ); \
  template T ReduceScatter( T sb, Op op, Comm comm ); \
  template T ReduceScatter( T sb, Comm comm ); \
  template void ReduceScatter( T* buf, int rc, Op op, Comm comm ); \
  template void ReduceScatter( T* buf, int rc, Comm comm ); \
  template void ReduceScatter( const T* sbuf, T* rbuf, const int* rcs, Op op, Comm comm ); \
  template void ReduceScatter( const T* sbuf, T* rbuf, const int* rcs, Comm comm ); \
  template void Scan( const T* sbuf, T* rbuf, int count, Op op, Comm comm ); \
  template void Scan( const T* sbuf, T* rbuf, int count, Comm comm ); \
  template T Scan( T sb, Op op, Comm comm ); \
  template T Scan( T sb, Comm comm ); \
  template void Scan( T* buf, int count, Op op, Comm comm ); \
  template void Scan( T* buf, int count, Comm comm );

MPI_PROTO(byte)
MPI_PROTO(int)
MPI_PROTO(unsigned)
MPI_PROTO(long int)
MPI_PROTO(unsigned long)
#ifdef EL_HAVE_MPI_LONG_LONG
MPI_PROTO(long long int)
MPI_PROTO(unsigned long long)
#endif
MPI_PROTO(float)
MPI_PROTO(double)
#ifdef EL_HAVE_QUAD
MPI_PROTO(Quad)
#endif
MPI_PROTO(Complex<float>)
MPI_PROTO(Complex<double>)
#ifdef EL_HAVE_QUAD
MPI_PROTO(Complex<Quad>)
#endif

MPI_PROTO(ValueInt<Int>)
MPI_PROTO(ValueInt<float>)
MPI_PROTO(ValueInt<double>)
#ifdef EL_HAVE_QUAD
MPI_PROTO(ValueInt<Quad>)
#endif
MPI_PROTO(ValueInt<Complex<float>>)
MPI_PROTO(ValueInt<Complex<double>>)
#ifdef EL_HAVE_QUAD
MPI_PROTO(ValueInt<Complex<Quad>>)
#endif
MPI_PROTO(Entry<Int>)
MPI_PROTO(Entry<float>)
MPI_PROTO(Entry<double>)
#ifdef EL_HAVE_QUAD
MPI_PROTO(Entry<Quad>)
#endif
MPI_PROTO(Entry<Complex<float>>)
MPI_PROTO(Entry<Complex<double>>)
#ifdef EL_HAVE_QUAD
MPI_PROTO(Entry<Complex<Quad>>)
#endif

#define PROTO(T) \
  template void SparseAllToAll \
  ( const vector<T>& sendBuffer, \
    const vector<int>& sendCounts, const vector<int>& sendDispls, \
          vector<T>& recvBuffer, \
    const vector<int>& recvCounts, const vector<int>& recvDispls, \
          mpi::Comm comm );
#include "El/macros/Instantiate.h"

} // namespace mpi
} // namespace El
