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
    DEBUG_ONLY(CallStackEntry cse("mpi::Create"))
    SafeMpi( MPI_Op_create( func, commutes, &op.op ) );
}

void Free( Op& op )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::Free"))
    SafeMpi( MPI_Op_free( &op.op ) );
}

void Free( Datatype& type )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::Free"))
    SafeMpi( MPI_Type_free( &type ) );
}

// Communicator manipulation 
// =========================

int WorldRank()
{
    DEBUG_ONLY(CallStackEntry cse("mpi::WorldRank"))
    return Rank( mpi::COMM_WORLD ); 
}

int Rank( Comm comm )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::Rank"))
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
    DEBUG_ONLY(CallStackEntry cse("mpi::Size"))
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
    DEBUG_ONLY(CallStackEntry cse("mpi::Create"))
    SafeMpi( 
        MPI_Comm_create( parentComm.comm, subsetGroup.group, &subsetComm.comm ) 
    );
}

void Dup( Comm original, Comm& duplicate )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::Dup"))
    SafeMpi( MPI_Comm_dup( original.comm, &duplicate.comm ) );
}

void Split( Comm comm, int color, int key, Comm& newComm )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::Split"))
    SafeMpi( MPI_Comm_split( comm.comm, color, key, &newComm.comm ) );
}

void Free( Comm& comm )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::Free"))
    SafeMpi( MPI_Comm_free( &comm.comm ) );
}

bool Congruent( Comm comm1, Comm comm2 )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::Congruent"))
    int result;
    SafeMpi( MPI_Comm_compare( comm1.comm, comm2.comm, &result ) );
    return ( result == MPI_IDENT || result == MPI_CONGRUENT );
}

void ErrorHandlerSet( Comm comm, ErrorHandler errorHandler )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::ErrorHandlerSet"))
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
    DEBUG_ONLY(CallStackEntry cse("mpi::CartCreate"))
    SafeMpi
    ( MPI_Cart_create
      ( comm.comm, numDims, const_cast<int*>(dimensions), 
        const_cast<int*>(periods), reorder, &cartComm.comm ) );
}

void CartSub( Comm comm, const int* remainingDims, Comm& subComm )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::CartSub"))
    SafeMpi( 
        MPI_Cart_sub
        ( comm.comm, const_cast<int*>(remainingDims), &subComm.comm ) 
    );
}

// Group manipulation 
// ==================

int Rank( Group group )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::Rank"))
    int rank;
    SafeMpi( MPI_Group_rank( group.group, &rank ) );
    return rank;
}

int Size( Group group )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::Size"))
    int size;
    SafeMpi( MPI_Group_size( group.group, &size ) );
    return size;
}

void CommGroup( Comm comm, Group& group )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::CommGroup"))
    SafeMpi( MPI_Comm_group( comm.comm, &group.group ) );
}

void Dup( Group group, Group& newGroup )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::Dup"))
    // For some reason, MPI_Group_dup does not exist
    Excl( group, 0, 0, newGroup ); 
}

void Union( Group groupA, Group groupB, Group& newGroup )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::Union"))
    SafeMpi( MPI_Group_union( groupA.group, groupB.group, &newGroup.group ) );
}

void Incl( Group group, int n, const int* ranks, Group& subGroup )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::Incl"))
    SafeMpi( 
        MPI_Group_incl
        ( group.group, n, const_cast<int*>(ranks), &subGroup.group ) 
    );
}

void Excl( Group group, int n, const int* ranks, Group& subGroup )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::Excl"))
    SafeMpi(
        MPI_Group_excl
        ( group.group, n, const_cast<int*>(ranks), &subGroup.group )
    );
}

void Difference( Group parent, Group subset, Group& complement )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::Difference"))
    SafeMpi( 
        MPI_Group_difference( parent.group, subset.group, &complement.group ) 
    );
}

void Free( Group& group )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::Free"))
    SafeMpi( MPI_Group_free( &group.group ) );
}

// Rank translations
// =================

int Translate( Group origGroup, int origRank, Group newGroup )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::Translate"))
    int newRank;
    Translate( origGroup, 1, &origRank, newGroup, &newRank );
    return newRank;
}

int Translate( Comm origComm, int origRank, Group newGroup )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::Translate"))
    int newRank;
    Translate( origComm, 1, &origRank, newGroup, &newRank );
    return newRank;
}

int Translate( Group origGroup, int origRank, Comm newComm )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::Translate"))
    int newRank;
    Translate( origGroup, 1, &origRank, newComm, &newRank );
    return newRank;
}

int Translate( Comm origComm, int origRank, Comm newComm )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::Translate"))
    int newRank;
    Translate( origComm, 1, &origRank, newComm, &newRank );
    return newRank;
}

void Translate
( Group origGroup, int size, const int* origRanks, 
  Group newGroup,                  int* newRanks )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::Translate"))
    SafeMpi
    ( MPI_Group_translate_ranks
      ( origGroup.group, size, const_cast<int*>(origRanks), 
        newGroup.group, newRanks ) );
}

void Translate
( Comm origComm,  int size, const int* origRanks, 
  Group newGroup,                 int* newRanks )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::Translate"))
    Group origGroup;
    CommGroup( origComm, origGroup );
    Translate( origGroup, size, origRanks, newGroup, newRanks );
    Free( origGroup );
}

void Translate
( Group origGroup,  int size, const int* origRanks, 
  Comm newComm,                     int* newRanks )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::Translate"))
    Group newGroup;
    CommGroup( newComm,  newGroup  );
    Translate( origGroup, size, origRanks, newGroup, newRanks );
    Free( newGroup  );
}

void Translate
( Comm origComm,  int size, const int* origRanks, 
  Comm newComm,                   int* newRanks )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::Translate"))
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
    DEBUG_ONLY(CallStackEntry cse("mpi::Barrier"))
    SafeMpi( MPI_Barrier( comm.comm ) );
}

// Test for completion
bool Test( Request& request )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::Test"))
    Status status;
    int flag;
    SafeMpi( MPI_Test( &request, &flag, &status ) );
    return flag;
}

// Ensure that the request finishes before continuing
void Wait( Request& request )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::Wait"))
    Status status;
    SafeMpi( MPI_Wait( &request, &status ) );
}

// Ensure that the request finishes before continuing
void Wait( Request& request, Status& status )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::Wait"))
    SafeMpi( MPI_Wait( &request, &status ) );
}

// Ensure that several requests finish before continuing
void WaitAll( int numRequests, Request* requests )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::WaitAll"))
    vector<Status> statuses( numRequests );
    SafeMpi( MPI_Waitall( numRequests, requests, statuses.data() ) );
}

// Ensure that several requests finish before continuing
void WaitAll( int numRequests, Request* requests, Status* statuses )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::WaitAll"))
    SafeMpi( MPI_Waitall( numRequests, requests, statuses ) );
}

// Nonblocking test for message completion
bool IProbe( int source, int tag, Comm comm, Status& status )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::IProbe"))
    int flag;
    SafeMpi( MPI_Iprobe( source, tag, comm.comm, &flag, &status ) );
    return flag;
}
bool IProbe( int source, Comm comm, Status& status )
{ return IProbe( source, 0, comm, status ); }

template<typename T>
int GetCount( Status& status )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::GetCount"))
    int count;
    SafeMpi( MPI_Get_count( &status, TypeMap<T>(), &count ) );
    return count;
}
template int GetCount<byte>( Status& status );
template int GetCount<int>( Status& status );
template int GetCount<unsigned>( Status& status );
template int GetCount<long int>( Status& status );
template int GetCount<unsigned long>( Status& status );
#ifdef EL_HAVE_MPI_LONG_LONG
template int GetCount<long long int>( Status& status );
template int GetCount<unsigned long long>( Status& status );
#endif
template int GetCount<float>( Status& status );
template int GetCount<double>( Status& status );
#ifdef EL_HAVE_QUAD
template int GetCount<Quad>( Status& status );
#endif
template int GetCount<Complex<float>>( Status& status );
template int GetCount<Complex<double>>( Status& status );
#ifdef EL_HAVE_QUAD
template int GetCount<Complex<Quad>>( Status& status );
#endif

template<typename Real>
void TaggedSend( const Real* buf, int count, int to, int tag, Comm comm )
{ 
    DEBUG_ONLY(CallStackEntry cse("mpi::Send"))
    SafeMpi( 
        MPI_Send
        ( const_cast<Real*>(buf), count, TypeMap<Real>(), to, tag, comm.comm )
    );
}

template<typename Real>
void TaggedSend
( const Complex<Real>* buf, int count, int to, int tag, Comm comm )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::Send"))
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

template void TaggedSend( const byte* buf, int count, int to, int tag, Comm comm );
template void TaggedSend( const int* buf, int count, int to, int tag, Comm comm );
template void TaggedSend( const unsigned* buf, int count, int to, int tag, Comm comm  );
template void TaggedSend( const long int* buf, int count, int to, int tag, Comm comm );
template void TaggedSend( const unsigned long* buf, int count, int to, int tag, Comm comm  );
#ifdef EL_HAVE_MPI_LONG_LONG
template void TaggedSend( const long long int* buf, int count, int to, int tag, Comm comm );
template void TaggedSend( const unsigned long long* buf, int count, int to, int tag, Comm comm  );
#endif
template void TaggedSend( const float* buf, int count, int to, int tag, Comm comm );
template void TaggedSend( const double* buf, int count, int to, int tag, Comm comm );
#ifdef EL_HAVE_QUAD
template void TaggedSend( const Quad* buf, int count, int to, int tag, Comm comm );
#endif
template void TaggedSend( const Complex<float>* buf, int count, int to, int tag, Comm comm );
template void TaggedSend( const Complex<double>* buf, int count, int to, int tag, Comm comm );
#ifdef EL_HAVE_QUAD
template void TaggedSend( const Complex<Quad>* buf, int count, int to, int tag, Comm comm );
#endif

template<typename T>
void Send( const T* buf, int count, int to, Comm comm )
{ TaggedSend( buf, count, to, 0, comm ); }

template void Send( const byte* buf, int count, int to, Comm comm );
template void Send( const int* buf, int count, int to, Comm comm );
template void Send( const unsigned* buf, int count, int to, Comm comm );
template void Send( const long int* buf, int count, int to, Comm comm );
template void Send( const unsigned long* buf, int count, int to, Comm comm );
#ifdef EL_HAVE_MPI_LONG_LONG
template void Send( const long long int* buf, int count, int to, Comm comm );
template void Send( const unsigned long long* buf, int count, int to, Comm comm );
#endif
template void Send( const float* buf, int count, int to, Comm comm );
template void Send( const double* buf, int count, int to, Comm comm );
#ifdef EL_HAVE_QUAD
template void Send( const Quad* buf, int count, int to, Comm comm );
#endif
template void Send( const Complex<float>* buf, int count, int to, Comm comm );
template void Send( const Complex<double>* buf, int count, int to, Comm comm );
#ifdef EL_HAVE_QUAD
template void Send( const Complex<Quad>* buf, int count, int to, Comm comm );
#endif

template<typename T>
void TaggedSend( T b, int to, int tag, Comm comm )
{ TaggedSend( &b, 1, to, tag, comm ); }

template void TaggedSend( byte b, int to, int tag, Comm comm );
template void TaggedSend( int b, int to, int tag, Comm comm );
template void TaggedSend( unsigned b, int to, int tag, Comm comm );
template void TaggedSend( long int b, int to, int tag, Comm comm );
template void TaggedSend( unsigned long b, int to, int tag, Comm comm );
#ifdef EL_HAVE_MPI_LONG_LONG
template void TaggedSend( long long int b, int to, int tag, Comm comm );
template void TaggedSend( unsigned long long b, int to, int tag, Comm comm );
#endif
template void TaggedSend( float b, int to, int tag, Comm comm );
template void TaggedSend( double b, int to, int tag, Comm comm );
#ifdef EL_HAVE_QUAD
template void TaggedSend( Quad b, int to, int tag, Comm comm );
#endif
template void TaggedSend( Complex<float> b, int to, int tag, Comm comm );
template void TaggedSend( Complex<double> b, int to, int tag, Comm comm );
#ifdef EL_HAVE_QUAD
template void TaggedSend( Complex<Quad> b, int to, int tag, Comm comm );
#endif

template<typename T>
void Send( T b, int to, Comm comm )
{ TaggedSend( b, to, 0, comm ); }

template void Send( byte b, int to, Comm comm );
template void Send( int b, int to, Comm comm );
template void Send( unsigned b, int to, Comm comm );
template void Send( long int b, int to, Comm comm );
template void Send( unsigned long b, int to, Comm comm );
#ifdef EL_HAVE_MPI_LONG_LONG
template void Send( long long int b, int to, Comm comm );
template void Send( unsigned long long b, int to, Comm comm );
#endif
template void Send( float b, int to, Comm comm );
template void Send( double b, int to, Comm comm );
#ifdef EL_HAVE_QUAD
template void Send( Quad b, int to, Comm comm );
#endif
template void Send( Complex<float> b, int to, Comm comm );
template void Send( Complex<double> b, int to, Comm comm );
#ifdef EL_HAVE_QUAD
template void Send( Complex<Quad> b, int to, Comm comm );
#endif

template<typename Real>
void TaggedISend
( const Real* buf, int count, int to, int tag, Comm comm, Request& request )
{ 
    DEBUG_ONLY(CallStackEntry cse("mpi::ISend"))
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
    DEBUG_ONLY(CallStackEntry cse("mpi::ISend"))
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

template void TaggedISend( const byte* buf, int count, int to, int tag, Comm comm, Request& request );
template void TaggedISend( const int* buf, int count, int to, int tag, Comm comm, Request& request );
template void TaggedISend( const unsigned* buf, int count, int to, int tag, Comm comm, Request& request );
template void TaggedISend( const long int* buf, int count, int to, int tag, Comm comm, Request& request );
template void TaggedISend( const unsigned long* buf, int count, int to, int tag, Comm comm, Request& request );
#ifdef EL_HAVE_MPI_LONG_LONG
template void TaggedISend( const long long int* buf, int count, int to, int tag, Comm comm, Request& request );
template void TaggedISend( const unsigned long long* buf, int count, int to, int tag, Comm comm, Request& request );
#endif
template void TaggedISend( const float* buf, int count, int to, int tag, Comm comm, Request& request );
template void TaggedISend( const double* buf, int count, int to, int tag, Comm comm, Request& request );
#ifdef EL_HAVE_QUAD
template void TaggedISend( const Quad* buf, int count, int to, int tag, Comm comm, Request& request );
#endif
template void TaggedISend( const Complex<float>* buf, int count, int to, int tag, Comm comm, Request& request );
template void TaggedISend( const Complex<double>* buf, int count, int to, int tag, Comm comm, Request& request );
#ifdef EL_HAVE_QUAD
template void TaggedISend( const Complex<Quad>* buf, int count, int to, int tag, Comm comm, Request& request );
#endif

template<typename T>
void ISend
( const T* buf, int count, int to, Comm comm, Request& request )
{ TaggedISend( buf, count, to, 0, comm, request ); } 

template void ISend( const byte* buf, int count, int to, Comm comm, Request& request );
template void ISend( const int* buf, int count, int to, Comm comm, Request& request );
template void ISend( const unsigned* buf, int count, int to, Comm comm, Request& request );
template void ISend( const long int* buf, int count, int to, Comm comm, Request& request );
template void ISend( const unsigned long* buf, int count, int to, Comm comm, Request& request );
#ifdef EL_HAVE_MPI_LONG_LONG
template void ISend( const long long int* buf, int count, int to, Comm comm, Request& request );
template void ISend( const unsigned long long* buf, int count, int to, Comm comm, Request& request );
#endif
template void ISend( const float* buf, int count, int to, Comm comm, Request& request );
template void ISend( const double* buf, int count, int to, Comm comm, Request& request );
#ifdef EL_HAVE_QUAD
template void ISend( const Quad* buf, int count, int to, Comm comm, Request& request );
#endif
template void ISend( const Complex<float>* buf, int count, int to, Comm comm, Request& request );
template void ISend( const Complex<double>* buf, int count, int to, Comm comm, Request& request );
#ifdef EL_HAVE_QUAD
template void ISend( const Complex<Quad>* buf, int count, int to, Comm comm, Request& request );
#endif

template<typename T>
void TaggedISend( T b, int to, int tag, Comm comm, Request& request )
{ TaggedISend( &b, 1, to, tag, comm, request ); }

template void TaggedISend( byte buf, int to, int tag, Comm comm, Request& request );
template void TaggedISend( int buf, int to, int tag, Comm comm, Request& request );
template void TaggedISend( unsigned buf, int to, int tag, Comm comm, Request& request );
template void TaggedISend( long int buf, int to, int tag, Comm comm, Request& request );
template void TaggedISend( unsigned long buf, int to, int tag, Comm comm, Request& request );
#ifdef EL_HAVE_MPI_LONG_LONG
template void TaggedISend( long long int buf, int to, int tag, Comm comm, Request& request );
template void TaggedISend( unsigned long long buf, int to, int tag, Comm comm, Request& request );
#endif
template void TaggedISend( float buf, int to, int tag, Comm comm, Request& request );
template void TaggedISend( double buf, int to, int tag, Comm comm, Request& request );
#ifdef EL_HAVE_QUAD
template void TaggedISend( Quad buf, int to, int tag, Comm comm, Request& request );
#endif
template void TaggedISend( Complex<float> buf, int to, int tag, Comm comm, Request& request );
template void TaggedISend( Complex<double> buf, int to, int tag, Comm comm, Request& request );
#ifdef EL_HAVE_QUAD
template void TaggedISend( Complex<Quad> buf, int to, int tag, Comm comm, Request& request );
#endif

template<typename T>
void ISend( T b, int to, Comm comm, Request& request )
{ TaggedISend( b, to, 0, comm, request ); }

template void ISend( byte buf, int to, Comm comm, Request& request );
template void ISend( int buf, int to, Comm comm, Request& request );
template void ISend( unsigned buf, int to, Comm comm, Request& request );
template void ISend( long int buf, int to, Comm comm, Request& request );
template void ISend( unsigned long buf, int to, Comm comm, Request& request );
#ifdef EL_HAVE_MPI_LONG_LONG
template void ISend( long long int buf, int to, Comm comm, Request& request );
template void ISend( unsigned long long buf, int to, Comm comm, Request& request );
#endif
template void ISend( float buf, int to, Comm comm, Request& request );
template void ISend( double buf, int to, Comm comm, Request& request );
#ifdef EL_HAVE_QUAD
template void ISend( Quad buf, int to, Comm comm, Request& request );
#endif
template void ISend( Complex<float> buf, int to, Comm comm, Request& request );
template void ISend( Complex<double> buf, int to, Comm comm, Request& request );
#ifdef EL_HAVE_QUAD
template void ISend( Complex<Quad> buf, int to, Comm comm, Request& request );
#endif

template<typename Real>
void TaggedISSend
( const Real* buf, int count, int to, int tag, Comm comm, Request& request )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::ISSend"))
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
    DEBUG_ONLY(CallStackEntry cse("mpi::ISSend"))
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

template void TaggedISSend( const byte* buf, int count, int to, int tag, Comm comm, Request& request );
template void TaggedISSend( const int* buf, int count, int to, int tag, Comm comm, Request& request );
template void TaggedISSend( const unsigned* buf, int count, int to, int tag, Comm comm, Request& request );
template void TaggedISSend( const long int* buf, int count, int to, int tag, Comm comm, Request& request );
template void TaggedISSend( const unsigned long* buf, int count, int to, int tag, Comm comm, Request& request );
#ifdef EL_HAVE_MPI_LONG_LONG
template void TaggedISSend( const long long int* buf, int count, int to, int tag, Comm comm, Request& request );
template void TaggedISSend( const unsigned long long* buf, int count, int to, int tag, Comm comm, Request& request );
#endif
template void TaggedISSend( const float* buf, int count, int to, int tag, Comm comm, Request& request );
template void TaggedISSend( const double* buf, int count, int to, int tag, Comm comm, Request& request );
#ifdef EL_HAVE_QUAD
template void TaggedISSend( const Quad* buf, int count, int to, int tag, Comm comm, Request& request );
#endif
template void TaggedISSend( const Complex<float>* buf, int count, int to, int tag, Comm comm, Request& request );
template void TaggedISSend( const Complex<double>* buf, int count, int to, int tag, Comm comm, Request& request );
#ifdef EL_HAVE_QUAD
template void TaggedISSend( const Complex<Quad>* buf, int count, int to, int tag, Comm comm, Request& request );
#endif

template<typename T>
void ISSend( const T* buf, int count, int to, Comm comm, Request& request )
{ TaggedISSend( buf, count, to, 0, comm, request ); }

template void ISSend( const byte* buf, int count, int to, Comm comm, Request& request );
template void ISSend( const int* buf, int count, int to, Comm comm, Request& request );
template void ISSend( const unsigned* buf, int count, int to, Comm comm, Request& request );
template void ISSend( const long int* buf, int count, int to, Comm comm, Request& request );
template void ISSend( const unsigned long* buf, int count, int to, Comm comm, Request& request );
#ifdef EL_HAVE_MPI_LONG_LONG
template void ISSend( const long long int* buf, int count, int to, Comm comm, Request& request );
template void ISSend( const unsigned long long* buf, int count, int to, Comm comm, Request& request );
#endif
template void ISSend( const float* buf, int count, int to, Comm comm, Request& request );
template void ISSend( const double* buf, int count, int to, Comm comm, Request& request );
#ifdef EL_HAVE_QUAD
template void ISSend( const Quad* buf, int count, int to, Comm comm, Request& request );
#endif
template void ISSend( const Complex<float>* buf, int count, int to, Comm comm, Request& request );
template void ISSend( const Complex<double>* buf, int count, int to, Comm comm, Request& request );
#ifdef EL_HAVE_QUAD
template void ISSend( const Complex<Quad>* buf, int count, int to, Comm comm, Request& request );
#endif

template<typename T>
void TaggedISSend( T b, int to, int tag, Comm comm, Request& request )
{ TaggedISSend( &b, 1, to, tag, comm, request ); }

template void TaggedISSend( byte b, int to, int tag, Comm comm, Request& request );
template void TaggedISSend( int b, int to, int tag, Comm comm, Request& request );
template void TaggedISSend( unsigned b, int to, int tag, Comm comm, Request& request );
template void TaggedISSend( long int b, int to, int tag, Comm comm, Request& request );
template void TaggedISSend( unsigned long b, int to, int tag, Comm comm, Request& request );
#ifdef EL_HAVE_MPI_LONG_LONG
template void TaggedISSend( long long int b, int to, int tag, Comm comm, Request& request );
template void TaggedISSend( unsigned long long b, int to, int tag, Comm comm, Request& request );
#endif
template void TaggedISSend( float b, int to, int tag, Comm comm, Request& request );
template void TaggedISSend( double b, int to, int tag, Comm comm, Request& request );
#ifdef EL_HAVE_QUAD
template void TaggedISSend( Quad b, int to, int tag, Comm comm, Request& request );
#endif
template void TaggedISSend( Complex<float> b, int to, int tag, Comm comm, Request& request );
template void TaggedISSend( Complex<double> b, int to, int tag, Comm comm, Request& request );
#ifdef EL_HAVE_QUAD
template void TaggedISSend( Complex<Quad> b, int to, int tag, Comm comm, Request& request );
#endif

template<typename Real>
void TaggedRecv( Real* buf, int count, int from, int tag, Comm comm )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::Recv"))
    Status status;
    SafeMpi
    ( MPI_Recv( buf, count, TypeMap<Real>(), from, tag, comm.comm, &status ) );
}

template<typename Real>
void TaggedRecv( Complex<Real>* buf, int count, int from, int tag, Comm comm )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::Recv"))
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

template void TaggedRecv( byte* buf, int count, int from, int tag, Comm comm );
template void TaggedRecv( int* buf, int count, int from, int tag, Comm comm );
template void TaggedRecv( unsigned* buf, int count, int from, int tag, Comm comm );
template void TaggedRecv( long int* buf, int count, int from, int tag, Comm comm );
template void TaggedRecv( unsigned long* buf, int count, int from, int tag, Comm comm );
#ifdef EL_HAVE_MPI_LONG_LONG
template void TaggedRecv( long long int* buf, int count, int from, int tag, Comm comm );
template void TaggedRecv( unsigned long long* buf, int count, int from, int tag, Comm comm );
#endif
template void TaggedRecv( float* buf, int count, int from, int tag, Comm comm );
template void TaggedRecv( double* buf, int count, int from, int tag, Comm comm );
#ifdef EL_HAVE_QUAD
template void TaggedRecv( Quad* buf, int count, int from, int tag, Comm comm );
#endif
template void TaggedRecv( Complex<float>* buf, int count, int from, int tag, Comm comm );
template void TaggedRecv( Complex<double>* buf, int count, int from, int tag, Comm comm );
#ifdef EL_HAVE_QUAD
template void TaggedRecv( Complex<Quad>* buf, int count, int from, int tag, Comm comm );
#endif

template<typename T>
void Recv( T* buf, int count, int from, Comm comm )
{ TaggedRecv( buf, count, from, mpi::ANY_TAG, comm ); }

template void Recv( byte* buf, int count, int from, Comm comm );
template void Recv( int* buf, int count, int from, Comm comm );
template void Recv( unsigned* buf, int count, int from, Comm comm );
template void Recv( long int* buf, int count, int from, Comm comm );
template void Recv( unsigned long* buf, int count, int from, Comm comm );
#ifdef EL_HAVE_MPI_LONG_LONG
template void Recv( long long int* buf, int count, int from, Comm comm );
template void Recv( unsigned long long* buf, int count, int from, Comm comm );
#endif
template void Recv( float* buf, int count, int from, Comm comm );
template void Recv( double* buf, int count, int from, Comm comm );
#ifdef EL_HAVE_QUAD
template void Recv( Quad* buf, int count, int from, Comm comm );
#endif
template void Recv( Complex<float>* buf, int count, int from, Comm comm );
template void Recv( Complex<double>* buf, int count, int from, Comm comm );
#ifdef EL_HAVE_QUAD
template void Recv( Complex<Quad>* buf, int count, int from, Comm comm );
#endif

template<typename T>
T TaggedRecv( int from, int tag, Comm comm )
{ T b; TaggedRecv( &b, 1, from, tag, comm ); return b; }

template byte TaggedRecv( int from, int tag, Comm comm );
template int TaggedRecv( int from, int tag, Comm comm );
template unsigned TaggedRecv( int from, int tag, Comm comm );
template long int TaggedRecv( int from, int tag, Comm comm );
template unsigned long TaggedRecv( int from, int tag, Comm comm );
#ifdef EL_HAVE_MPI_LONG_LONG
template long long int TaggedRecv( int from, int tag, Comm comm );
template unsigned long long TaggedRecv( int from, int tag, Comm comm );
#endif
template float TaggedRecv( int from, int tag, Comm comm );
template double TaggedRecv( int from, int tag, Comm comm );
#ifdef EL_HAVE_QUAD
template Quad TaggedRecv( int from, int tag, Comm comm );
#endif
template Complex<float> TaggedRecv( int from, int tag, Comm comm );
template Complex<double> TaggedRecv( int from, int tag, Comm comm );
#ifdef EL_HAVE_QUAD
template Complex<Quad> TaggedRecv( int from, int tag, Comm comm );
#endif

template<typename T>
T Recv( int from, Comm comm )
{ return TaggedRecv<T>( from, mpi::ANY_TAG, comm ); }

template byte Recv( int from, Comm comm );
template int Recv( int from, Comm comm );
template unsigned Recv( int from, Comm comm );
template long int Recv( int from, Comm comm );
template unsigned long Recv( int from, Comm comm );
#ifdef EL_HAVE_MPI_LONG_LONG
template long long int Recv( int from, Comm comm );
template unsigned long long Recv( int from, Comm comm );
#endif
template float Recv( int from, Comm comm );
template double Recv( int from, Comm comm );
#ifdef EL_HAVE_QUAD
template Quad Recv( int from, Comm comm );
#endif
template Complex<float> Recv( int from, Comm comm );
template Complex<double> Recv( int from, Comm comm );
#ifdef EL_HAVE_QUAD
template Complex<Quad> Recv( int from, Comm comm );
#endif

template<typename Real>
void TaggedIRecv
( Real* buf, int count, int from, int tag, Comm comm, Request& request )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::IRecv"))
    SafeMpi
    ( MPI_Irecv
      ( buf, count, TypeMap<Real>(), from, tag, comm.comm, &request ) );
}

template<typename Real>
void TaggedIRecv
( Complex<Real>* buf, int count, int from, int tag, 
  Comm comm, Request& request )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::IRecv"))
#ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi
    ( MPI_Irecv( buf, 2*count, TypeMap<Real>(), from, tag, comm.comm, &request ) );
#else
    SafeMpi
    ( MPI_Irecv
      ( buf, count, TypeMap<Complex<Real>>(), from, tag, comm.comm, &request ) );
#endif
}

template void TaggedIRecv( byte* buf, int count, int from, int tag, Comm comm, Request& request );
template void TaggedIRecv( int* buf, int count, int from, int tag, Comm comm, Request& request );
template void TaggedIRecv( unsigned* buf, int count, int from, int tag, Comm comm, Request& request );
template void TaggedIRecv( long int* buf, int count, int from, int tag, Comm comm, Request& request );
template void TaggedIRecv( unsigned long* buf, int count, int from, int tag, Comm comm, Request& request );
#ifdef EL_HAVE_MPI_LONG_LONG
template void TaggedIRecv( long long int* buf, int count, int from, int tag, Comm comm, Request& request );
template void TaggedIRecv( unsigned long long* buf, int count, int from, int tag, Comm comm, Request& request );
#endif
template void TaggedIRecv( float* buf, int count, int from, int tag, Comm comm, Request& request );
template void TaggedIRecv( double* buf, int count, int from, int tag, Comm comm, Request& request );
#ifdef EL_HAVE_QUAD
template void TaggedIRecv( Quad* buf, int count, int from, int tag, Comm comm, Request& request );
#endif
template void TaggedIRecv( Complex<float>* buf, int count, int from, int tag, Comm comm, Request& request );
template void TaggedIRecv( Complex<double>* buf, int count, int from, int tag, Comm comm, Request& request );
#ifdef EL_HAVE_QUAD
template void TaggedIRecv( Complex<Quad>* buf, int count, int from, int tag, Comm comm, Request& request );
#endif

template<typename T>
void IRecv( T* buf, int count, int from, Comm comm, Request& request )
{ TaggedIRecv( buf, count, from, mpi::ANY_TAG, comm, request ); }

template void IRecv( byte* buf, int count, int from, Comm comm, Request& request );
template void IRecv( int* buf, int count, int from, Comm comm, Request& request );
template void IRecv( unsigned* buf, int count, int from, Comm comm, Request& request );
template void IRecv( long int* buf, int count, int from, Comm comm, Request& request );
template void IRecv( unsigned long* buf, int count, int from, Comm comm, Request& request );
#ifdef EL_HAVE_MPI_LONG_LONG
template void IRecv( long long int* buf, int count, int from, Comm comm, Request& request );
template void IRecv( unsigned long long* buf, int count, int from, Comm comm, Request& request );
#endif
template void IRecv( float* buf, int count, int from, Comm comm, Request& request );
template void IRecv( double* buf, int count, int from, Comm comm, Request& request );
#ifdef EL_HAVE_QUAD
template void IRecv( Quad* buf, int count, int from, Comm comm, Request& request );
#endif
template void IRecv( Complex<float>* buf, int count, int from, Comm comm, Request& request );
template void IRecv( Complex<double>* buf, int count, int from, Comm comm, Request& request );
#ifdef EL_HAVE_QUAD
template void IRecv( Complex<Quad>* buf, int count, int from, Comm comm, Request& request );
#endif

template<typename T>
T TaggedIRecv( int from, int tag, Comm comm, Request& request )
{ T b; TaggedIRecv( &b, 1, from, tag, comm, request ); return b; }

template byte TaggedIRecv( int from, int tag, Comm comm, Request& request );
template int TaggedIRecv( int from, int tag, Comm comm, Request& request );
template unsigned TaggedIRecv( int from, int tag, Comm comm, Request& request );
template long int TaggedIRecv( int from, int tag, Comm comm, Request& request );
template unsigned long TaggedIRecv( int from, int tag, Comm comm, Request& request );
#ifdef EL_HAVE_MPI_LONG_LONG
template long long int TaggedIRecv( int from, int tag, Comm comm, Request& request );
template unsigned long long TaggedIRecv( int from, int tag, Comm comm, Request& request );
#endif
template float TaggedIRecv( int from, int tag, Comm comm, Request& request );
template double TaggedIRecv( int from, int tag, Comm comm, Request& request );
#ifdef EL_HAVE_QUAD
template Quad TaggedIRecv( int from, int tag, Comm comm, Request& request );
#endif
template Complex<float> TaggedIRecv( int from, int tag, Comm comm, Request& request );
template Complex<double> TaggedIRecv( int from, int tag, Comm comm, Request& request );
#ifdef EL_HAVE_QUAD
template Complex<Quad> TaggedIRecv( int from, int tag, Comm comm, Request& request );
#endif

template<typename T>
T IRecv( int from, Comm comm, Request& request )
{ return TaggedIRecv<T>( from, mpi::ANY_TAG, comm, request ); }

template byte IRecv( int from, Comm comm, Request& request );
template int IRecv( int from, Comm comm, Request& request );
template unsigned IRecv( int from, Comm comm, Request& request );
template long int IRecv( int from, Comm comm, Request& request );
template unsigned long IRecv( int from, Comm comm, Request& request );
#ifdef EL_HAVE_MPI_LONG_LONG
template long long int IRecv( int from, Comm comm, Request& request );
template unsigned long long IRecv( int from, Comm comm, Request& request );
#endif
template float IRecv( int from, Comm comm, Request& request );
template double IRecv( int from, Comm comm, Request& request );
#ifdef EL_HAVE_QUAD
template Quad IRecv( int from, Comm comm, Request& request );
#endif
template Complex<float> IRecv( int from, Comm comm, Request& request );
template Complex<double> IRecv( int from, Comm comm, Request& request );
#ifdef EL_HAVE_QUAD
template Complex<Quad> IRecv( int from, Comm comm, Request& request );
#endif

template<typename Real>
void TaggedSendRecv
( const Real* sbuf, int sc, int to,   int stag,
        Real* rbuf, int rc, int from, int rtag, Comm comm )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::SendRecv"))
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
    DEBUG_ONLY(CallStackEntry cse("mpi::SendRecv"))
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

template void TaggedSendRecv
( const byte* sbuf, int sc, int to,   int stag, 
        byte* rbuf, int rc, int from, int rtag, Comm comm );
template void TaggedSendRecv
( const int* sbuf, int sc, int to,   int stag, 
        int* rbuf, int rc, int from, int rtag, Comm comm );
template void TaggedSendRecv
( const unsigned* sbuf, int sc, int to,   int stag, 
        unsigned* rbuf, int rc, int from, int rtag, Comm comm );
template void TaggedSendRecv
( const long int* sbuf, int sc, int to,   int stag, 
        long int* rbuf, int rc, int from, int rtag, Comm comm );
template void TaggedSendRecv
( const unsigned long* sbuf, int sc, int to,   int stag, 
        unsigned long* rbuf, int rc, int from, int rtag, Comm comm );
#ifdef EL_HAVE_MPI_LONG_LONG
template void TaggedSendRecv
( const long long int* sbuf, int sc, int to,   int stag, 
        long long int* rbuf, int rc, int from, int rtag, Comm comm );
template void TaggedSendRecv
( const unsigned long long* sbuf, int sc, int to,   int stag, 
        unsigned long long* rbuf, int rc, int from, int rtag, Comm comm );
#endif
template void TaggedSendRecv
( const float* sbuf, int sc, int to,   int stag, 
        float* rbuf, int rc, int from, int rtag, Comm comm );
template void TaggedSendRecv
( const double* sbuf, int sc, int to,   int stag, 
        double* rbuf, int rc, int from, int rtag, Comm comm );
#ifdef EL_HAVE_QUAD
template void TaggedSendRecv
( const Quad* sbuf, int sc, int to,   int stag, 
        Quad* rbuf, int rc, int from, int rtag, Comm comm );
#endif
template void TaggedSendRecv
( const Complex<float>* sbuf, int sc, int to,   int stag, 
        Complex<float>* rbuf, int rc, int from, int rtag, Comm comm );
template void TaggedSendRecv
( const Complex<double>* sbuf, int sc, int to, int stag, 
        Complex<double>* rbuf, int rc, int from, int rtag, Comm comm );
#ifdef EL_HAVE_QUAD
template void TaggedSendRecv
( const Complex<Quad>* sbuf, int sc, int to,   int stag, 
        Complex<Quad>* rbuf, int rc, int from, int rtag, Comm comm );
#endif

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
#ifdef EL_HAVE_MPI_LONG_LONG
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
#ifdef EL_HAVE_QUAD
template void SendRecv
( const Quad* sbuf, int sc, int to,
        Quad* rbuf, int rc, int from, Comm comm );
#endif
template void SendRecv
( const Complex<float>* sbuf, int sc, int to, 
        Complex<float>* rbuf, int rc, int from, Comm comm );
template void SendRecv
( const Complex<double>* sbuf, int sc, int to, 
        Complex<double>* rbuf, int rc, int from, Comm comm );
#ifdef EL_HAVE_QUAD
template void SendRecv
( const Complex<Quad>* sbuf, int sc, int to,
        Complex<Quad>* rbuf, int rc, int from, Comm comm );
#endif

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
#ifdef EL_HAVE_MPI_LONG_LONG
template long long int TaggedSendRecv
( long long int sb, int to, int stag, int from, int rtag, Comm comm );
template unsigned long long TaggedSendRecv
( unsigned long long sb, int to, int stag, int from, int rtag, Comm comm );
#endif
template float TaggedSendRecv
( float sb, int to, int stag, int from, int rtag, Comm comm );
template double TaggedSendRecv
( double sb, int to, int stag, int from, int rtag, Comm comm );
#ifdef EL_HAVE_QUAD
template Quad TaggedSendRecv
( Quad sb, int to, int stag, int from, int rtag, Comm comm );
#endif
template Complex<float> TaggedSendRecv
( Complex<float> sb, int to, int stag, int from, int rtag, Comm comm );
template Complex<double> TaggedSendRecv
( Complex<double> sb, int to, int stag, int from, int rtag, Comm comm );
#ifdef EL_HAVE_QUAD
template Complex<Quad> TaggedSendRecv
( Complex<Quad> sb, int to, int stag, int from, int rtag, Comm comm );
#endif

template<typename T>
T SendRecv( T sb, int to, int from, Comm comm )
{ return TaggedSendRecv( sb, to, 0, from, mpi::ANY_TAG, comm ); }

template byte SendRecv( byte sb, int to, int from, Comm comm );
template int SendRecv( int sb, int to, int from, Comm comm );
template unsigned SendRecv( unsigned sb, int to, int from, Comm comm );
template long int SendRecv( long int sb, int to, int from, Comm comm );
template unsigned long SendRecv( unsigned long sb, int to, int from, Comm comm );
#ifdef EL_HAVE_MPI_LONG_LONG
template long long int SendRecv( long long int sb, int to, int from, Comm comm );
template unsigned long long SendRecv( unsigned long long sb, int to, int from, Comm comm );
#endif
template float SendRecv( float sb, int to, int from, Comm comm );
template double SendRecv( double sb, int to, int from, Comm comm );
#ifdef EL_HAVE_QUAD
template Quad SendRecv( Quad sb, int to, int from, Comm comm );
#endif
template Complex<float> SendRecv
( Complex<float> sb, int to, int from, Comm comm );
template Complex<double> SendRecv
( Complex<double> sb, int to, int from, Comm comm );
#ifdef EL_HAVE_QUAD
template Complex<Quad> SendRecv( Complex<Quad> sb, int to, int from, Comm comm );
#endif

template<typename Real>
void TaggedSendRecv
( Real* buf, int count, int to, int stag, int from, int rtag, Comm comm )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::SendRecv"))
    Status status;
    SafeMpi
    ( MPI_Sendrecv_replace
      ( buf, count, TypeMap<Real>(), to, stag, from, rtag, comm.comm, &status ) );
}

template<typename Real>
void TaggedSendRecv
( Complex<Real>* buf, int count, int to, int stag, int from, int rtag, Comm comm )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::SendRecv"))
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
#ifdef EL_HAVE_MPI_LONG_LONG
template void TaggedSendRecv
( long long int* buf, int count, int to, int stag, int from, int rtag, Comm comm );
template void TaggedSendRecv
( unsigned long long* buf, int count, int to, int stag, int from, int rtag, Comm comm );
#endif
template void TaggedSendRecv
( float* buf, int count, int to, int stag, int from, int rtag, Comm comm );
template void TaggedSendRecv
( double* buf, int count, int to, int stag, int from, int rtag, Comm comm );
#ifdef EL_HAVE_QUAD
template void TaggedSendRecv
( Quad* buf, int count, int to, int stag, int from, int rtag, Comm comm );
#endif
template void TaggedSendRecv
( Complex<float>* buf, int count, int to, int stag, 
  int from, int rtag, Comm comm );
template void TaggedSendRecv
( Complex<double>* buf, int count, int to, int stag, 
  int from, int rtag, Comm comm );
#ifdef EL_HAVE_QUAD
template void TaggedSendRecv
( Complex<Quad>* buf, int count, int to, int stag, 
  int from, int rtag, Comm comm );
#endif

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
#ifdef EL_HAVE_MPI_LONG_LONG
template void SendRecv
( long long int* buf, int count, int to, int from, Comm comm );
template void SendRecv
( unsigned long long* buf, int count, int to, int from, Comm comm );
#endif
template void SendRecv
( float* buf, int count, int to, int from, Comm comm );
template void SendRecv
( double* buf, int count, int to, int from, Comm comm );
#ifdef EL_HAVE_QUAD
template void SendRecv
( Quad* buf, int count, int to, int from, Comm comm );
#endif
template void SendRecv
( Complex<float>* buf, int count, int to, int from, Comm comm );
template void SendRecv
( Complex<double>* buf, int count, int to, int from, Comm comm );
#ifdef EL_HAVE_QUAD
template void SendRecv
( Complex<Quad>* buf, int count, int to, int from, Comm comm );
#endif

template<typename Real>
void Broadcast( Real* buf, int count, int root, Comm comm )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::Broadcast"))
    SafeMpi( MPI_Bcast( buf, count, TypeMap<Real>(), root, comm.comm ) );
}

template<typename Real>
void Broadcast( Complex<Real>* buf, int count, int root, Comm comm )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::Broadcast"))
#ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi( MPI_Bcast( buf, 2*count, TypeMap<Real>(), root, comm.comm ) );
#else
    SafeMpi( MPI_Bcast( buf, count, TypeMap<Complex<Real>>(), root, comm.comm ) );
#endif
}

template void Broadcast( byte* buf, int count, int root, Comm comm );
template void Broadcast( int* buf, int count, int root, Comm comm );
template void Broadcast( unsigned* buf, int count, int root, Comm comm );
template void Broadcast( long int* buf, int count, int root, Comm comm );
template void Broadcast( unsigned long* buf, int count, int root, Comm comm );
#ifdef EL_HAVE_MPI_LONG_LONG
template void Broadcast( long long int* buf, int count, int root, Comm comm );
template void Broadcast( unsigned long long* buf, int count, int root, Comm comm );
#endif
template void Broadcast( float* buf, int count, int root, Comm comm );
template void Broadcast( double* buf, int count, int root, Comm comm );
#ifdef EL_HAVE_QUAD
template void Broadcast( Quad* buf, int count, int root, Comm comm );
#endif
template void Broadcast( Complex<float>* buf, int count, int root, Comm comm );
template void Broadcast( Complex<double>* buf, int count, int root, Comm comm );
#ifdef EL_HAVE_QUAD
template void Broadcast( Complex<Quad>* buf, int count, int root, Comm comm );
#endif

template<typename T>
void Broadcast( T& b, int root, Comm comm )
{ Broadcast( &b, 1, root, comm ); }

template void Broadcast( byte& b, int root, Comm comm );
template void Broadcast( int& b, int root, Comm comm );
template void Broadcast( unsigned& b, int root, Comm comm );
template void Broadcast( long int& b, int root, Comm comm );
template void Broadcast( unsigned long& b, int root, Comm comm );
#ifdef EL_HAVE_MPI_LONG_LONG
template void Broadcast( long long int& b, int root, Comm comm );
template void Broadcast( unsigned long long& b, int root, Comm comm );
#endif
template void Broadcast( float& b, int root, Comm comm );
template void Broadcast( double& b, int root, Comm comm );
#ifdef EL_HAVE_QUAD
template void Broadcast( Quad& b, int root, Comm comm );
#endif
template void Broadcast( Complex<float>& b, int root, Comm comm );
template void Broadcast( Complex<double>& b, int root, Comm comm );
#ifdef EL_HAVE_QUAD
template void Broadcast( Complex<Quad>& b, int root, Comm comm );
#endif

#ifdef EL_HAVE_NONBLOCKING_COLLECTIVES
template<typename Real>
void IBroadcast( Real* buf, int count, int root, Comm comm, Request& request )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::IBroadcast"))
    SafeMpi
    ( MPI_Ibcast( buf, count, TypeMap<Real>(), root, comm.comm, &request ) );
}

template<typename Real>
void IBroadcast
( Complex<Real>* buf, int count, int root, Comm comm, Request& request )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::IBroadcast"))
#ifdef EL_AVOID_COMPLEX_MPI
    SafeMpi
    ( MPI_Ibcast( buf, 2*count, TypeMap<Real>(), root, comm.comm, &request ) );
#else
    SafeMpi
    ( MPI_Ibcast
      ( buf, count, TypeMap<Complex<Real>>(), root, comm.comm, &request ) );
#endif
}

template void IBroadcast( byte* buf, int count, int root, Comm comm, Request& request );
template void IBroadcast( int* buf, int count, int root, Comm comm, Request& request );
template void IBroadcast( unsigned* buf, int count, int root, Comm comm, Request& request );
template void IBroadcast( long int* buf, int count, int root, Comm comm, Request& request );
template void IBroadcast( unsigned long* buf, int count, int root, Comm comm, Request& request );
#ifdef EL_HAVE_MPI_LONG_LONG
template void IBroadcast( long long int* buf, int count, int root, Comm comm, Request& request );
template void IBroadcast( unsigned long long* buf, int count, int root, Comm comm, Request& request );
#endif
template void IBroadcast( float* buf, int count, int root, Comm comm, Request& request );
template void IBroadcast( double* buf, int count, int root, Comm comm, Request& request );
#ifdef EL_HAVE_QUAD
template void IBroadcast( Quad* buf, int count, int root, Comm comm, Request& request );
#endif
template void IBroadcast( Complex<float>* buf, int count, int root, Comm comm, Request& request );
template void IBroadcast( Complex<double>* buf, int count, int root, Comm comm, Request& request );
#ifdef EL_HAVE_QUAD
template void IBroadcast( Complex<Quad>* buf, int count, int root, Comm comm, Request& request );
#endif

template<typename T>
void IBroadcast( T& b, int root, Comm comm, Request& request )
{ IBroadcast( &b, 1, root, comm, request ); }

template void IBroadcast( byte& b, int root, Comm comm, Request& request );
template void IBroadcast( int& b, int root, Comm comm, Request& request );
template void IBroadcast( unsigned& b, int root, Comm comm, Request& request );
template void IBroadcast( long int& b, int root, Comm comm, Request& request );
template void IBroadcast( unsigned long& b, int root, Comm comm, Request& request );
#ifdef EL_HAVE_MPI_LONG_LONG
template void IBroadcast( long long int& b, int root, Comm comm, Request& request );
template void IBroadcast( unsigned long long& b, int root, Comm comm, Request& request );
#endif
template void IBroadcast( float& b, int root, Comm comm, Request& request );
template void IBroadcast( double& b, int root, Comm comm, Request& request );
#ifdef EL_HAVE_QUAD
template void IBroadcast( Quad& b, int root, Comm comm, Request& request );
#endif
template void IBroadcast( Complex<float>& b, int root, Comm comm, Request& request );
template void IBroadcast( Complex<double>& b, int root, Comm comm, Request& request );
#ifdef EL_HAVE_QUAD
template void IBroadcast( Complex<Quad>& b, int root, Comm comm, Request& request );
#endif
#endif // ifdef EL_HAVE_NONBLOCKING_COLLECTIVES

template<typename Real>
void Gather
( const Real* sbuf, int sc,
        Real* rbuf, int rc, int root, Comm comm )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::Gather"))
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
    DEBUG_ONLY(CallStackEntry cse("mpi::Gather"))
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

template void Gather( const byte* sbuf, int sc, byte* rbuf, int rc, int root, Comm comm );
template void Gather( const int* sbuf, int sc, int* rbuf, int rc, int root, Comm comm );
template void Gather( const unsigned* sbuf, int sc, unsigned* rbuf, int rc, int root, Comm comm );
template void Gather( const long int* sbuf, int sc, long int* rbuf, int rc, int root, Comm comm );
template void Gather( const unsigned long* sbuf, int sc, unsigned long* rbuf, int rc, int root, Comm comm );
#ifdef EL_HAVE_MPI_LONG_LONG
template void Gather( const long long int* sbuf, int sc, long long int* rbuf, int rc, int root, Comm comm );
template void Gather( const unsigned long long* sbuf, int sc, unsigned long long* rbuf, int rc, int root, Comm comm );
#endif
template void Gather( const float* sbuf, int sc, float* rbuf, int rc, int root, Comm comm );
template void Gather( const double* sbuf, int sc, double* rbuf, int rc, int root, Comm comm );
#ifdef EL_HAVE_QUAD
template void Gather( const Quad* sbuf, int sc, Quad* rbuf, int rc, int root, Comm comm );
#endif
template void Gather( const Complex<float>* sbuf, int sc, Complex<float>* rbuf, int rc, int root, Comm comm );
template void Gather( const Complex<double>* sbuf, int sc, Complex<double>* rbuf, int rc, int root, Comm comm );
#ifdef EL_HAVE_QUAD
template void Gather( const Complex<Quad>* sbuf, int sc, Complex<Quad>* rbuf, int rc, int root, Comm comm );
#endif

#ifdef EL_HAVE_NONBLOCKING_COLLECTIVES
template<typename Real>
void IGather
( const Real* sbuf, int sc,
        Real* rbuf, int rc, int root, Comm comm, Request& request )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::IGather"))
    SafeMpi
    ( MPI_Igather
      ( const_cast<Real*>(sbuf), sc, TypeMap<Real>(),
        rbuf,                    rc, TypeMap<Real>(), root, comm.comm, &request ) );
}

template<typename Real>
void IGather
( const Complex<Real>* sbuf, int sc,
        Complex<Real>* rbuf, int rc, int root, Comm comm, Request& request )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::IGather"))
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
#ifdef EL_HAVE_MPI_LONG_LONG
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
#ifdef EL_HAVE_QUAD
template void IGather
( const Quad* sbuf, int sc, 
        Quad* rbuf, int rc, int root, Comm comm, Request& request );
#endif
template void IGather
( const Complex<float>* sbuf, int sc, 
        Complex<float>* rbuf, int rc, int root, Comm comm, Request& request );
template void IGather
( const Complex<double>* sbuf, int sc, 
        Complex<double>* rbuf, int rc, int root, Comm comm, Request& request );
#ifdef EL_HAVE_QUAD
template void IGather
( const Complex<Quad>* sbuf, int sc, 
        Complex<Quad>* rbuf, int rc, int root, Comm comm, Request& request );
#endif
#endif // ifdef EL_HAVE_NONBLOCKING_COLLECTIVES

template<typename Real>
void Gather
( const Real* sbuf, int sc,
        Real* rbuf, const int* rcs, const int* rds, int root, Comm comm )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::Gather"))
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
    DEBUG_ONLY(CallStackEntry cse("mpi::Gather"))
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
#ifdef EL_HAVE_MPI_LONG_LONG
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
#ifdef EL_HAVE_QUAD
template void Gather
( const Quad* sbuf, int sc, 
        Quad* rbuf, const int* rcs, const int* rds, int root, Comm comm );
#endif
template void Gather
( const Complex<float>* sbuf, int sc, 
        Complex<float>* rbuf, const int* rcs, const int* rds, 
  int root, Comm comm );
template void Gather
( const Complex<double>* sbuf, int sc, 
        Complex<double>* rbuf, const int* rcs, const int* rds, 
  int root, Comm comm );
#ifdef EL_HAVE_QUAD
template void Gather
( const Complex<Quad>* sbuf, int sc, 
        Complex<Quad>* rbuf, const int* rcs, const int* rds, 
  int root, Comm comm );
#endif

template<typename Real>
void AllGather
( const Real* sbuf, int sc,
        Real* rbuf, int rc, Comm comm )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::AllGather"))
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
    DEBUG_ONLY(CallStackEntry cse("mpi::AllGather"))
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

template void AllGather( const byte* sbuf, int sc, byte* rbuf, int rc, Comm comm );
template void AllGather( const int* sbuf, int sc, int* rbuf, int rc, Comm comm );
template void AllGather( const unsigned* sbuf, int sc, unsigned* rbuf, int rc, Comm comm );
template void AllGather( const long int* sbuf, int sc, long int* rbuf, int rc, Comm comm );
template void AllGather( const unsigned long* sbuf, int sc, unsigned long* rbuf, int rc, Comm comm );
#ifdef EL_HAVE_MPI_LONG_LONG
template void AllGather( const long long int* sbuf, int sc, long long int* rbuf, int rc, Comm comm );
template void AllGather( const unsigned long long* sbuf, int sc, unsigned long long* rbuf, int rc, Comm comm );
#endif
template void AllGather( const float* sbuf, int sc, float* rbuf, int rc, Comm comm );
template void AllGather( const double* sbuf, int sc, double* rbuf, int rc, Comm comm );
#ifdef EL_HAVE_QUAD
template void AllGather( const Quad* sbuf, int sc, Quad* rbuf, int rc, Comm comm );
#endif
template void AllGather( const Complex<float>* sbuf, int sc, Complex<float>* rbuf, int rc, Comm comm );
template void AllGather( const Complex<double>* sbuf, int sc, Complex<double>* rbuf, int rc, Comm comm );
#ifdef EL_HAVE_QUAD
template void AllGather( const Complex<Quad>* sbuf, int sc, Complex<Quad>* rbuf, int rc, Comm comm );
#endif

template<typename Real>
void AllGather
( const Real* sbuf, int sc,
        Real* rbuf, const int* rcs, const int* rds, Comm comm )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::AllGather"))
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
    DEBUG_ONLY(CallStackEntry cse("mpi::AllGather"))
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
#ifdef EL_HAVE_MPI_LONG_LONG
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
#ifdef EL_HAVE_QUAD
template void AllGather
( const Quad* sbuf, int sc, 
        Quad* rbuf, const int* rcs, const int* rds, Comm comm );
#endif
template void AllGather
( const Complex<float>* sbuf, int sc, 
        Complex<float>* rbuf, const int* rcs, const int* rds, Comm comm );
template void AllGather
( const Complex<double>* sbuf, int sc, 
        Complex<double>* rbuf, const int* rcs, const int* rds, Comm comm );
#ifdef EL_HAVE_QUAD
template void AllGather
( const Complex<Quad>* sbuf, int sc, 
        Complex<Quad>* rbuf, const int* rcs, const int* rds, Comm comm );
#endif

template<typename Real>
void Scatter
( const Real* sbuf, int sc,
        Real* rbuf, int rc, int root, Comm comm )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::Scatter"))
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
    DEBUG_ONLY(CallStackEntry cse("mpi::Scatter"))
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
#ifdef EL_HAVE_MPI_LONG_LONG
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
#ifdef EL_HAVE_QUAD
template void Scatter
( const Quad* sbuf, int sc, 
        Quad* rbuf, int rc, int root, Comm comm );
#endif
template void Scatter
( const Complex<float>* sbuf, int sc, 
        Complex<float>* rbuf, int rc, int root, Comm comm );
template void Scatter
( const Complex<double>* sbuf, int sc, 
        Complex<double>* rbuf, int rc, int root, Comm comm );
#ifdef EL_HAVE_QUAD
template void Scatter
( const Complex<Quad>* sbuf, int sc, 
        Complex<Quad>* rbuf, int rc, int root, Comm comm );
#endif

template<typename Real>
void Scatter( Real* buf, int sc, int rc, int root, Comm comm )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::Scatter"))
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
    DEBUG_ONLY(CallStackEntry cse("mpi::Scatter"))
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

template void Scatter( byte* buf, int sc, int rc, int root, Comm comm );
template void Scatter( int* buf, int sc, int rc, int root, Comm comm );
template void Scatter( unsigned* buf, int sc, int rc, int root, Comm comm );
template void Scatter( long int* buf, int sc, int rc, int root, Comm comm );
template void Scatter( unsigned long* buf, int sc, int rc, int root, Comm comm );
#ifdef EL_HAVE_MPI_LONG_LONG
template void Scatter( long long int* buf, int sc, int rc, int root, Comm comm );
template void Scatter( unsigned long long* buf, int sc, int rc, int root, Comm comm );
#endif
template void Scatter( float* buf, int sc, int rc, int root, Comm comm );
template void Scatter( double* buf, int sc, int rc, int root, Comm comm );
#ifdef EL_HAVE_QUAD
template void Scatter( Quad* buf, int sc, int rc, int root, Comm comm );
#endif
template void Scatter( Complex<float>* buf, int sc, int rc, int root, Comm comm );
template void Scatter( Complex<double>* buf, int sc, int rc, int root, Comm comm );
#ifdef EL_HAVE_QUAD
template void Scatter( Complex<Quad>* buf, int sc, int rc, int root, Comm comm );
#endif

template<typename Real>
void AllToAll
( const Real* sbuf, int sc,
        Real* rbuf, int rc, Comm comm )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::AllToAll"))
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
    DEBUG_ONLY(CallStackEntry cse("mpi::AllToAll"))
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
#ifdef EL_HAVE_MPI_LONG_LONG
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
#ifdef EL_HAVE_QUAD
template void AllToAll
( const Quad* sbuf, int sc, 
        Quad* rbuf, int rc, Comm comm );
#endif
template void AllToAll
( const Complex<float>* sbuf, int sc, 
        Complex<float>* rbuf, int rc, Comm comm );
template void AllToAll
( const Complex<double>* sbuf, int sc, 
        Complex<double>* rbuf, int rc, Comm comm );
#ifdef EL_HAVE_QUAD
template void AllToAll
( const Complex<Quad>* sbuf, int sc, 
        Complex<Quad>* rbuf, int rc, Comm comm );
#endif

template<typename Real>
void AllToAll
( const Real* sbuf, const int* scs, const int* sds, 
        Real* rbuf, const int* rcs, const int* rds, Comm comm )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::AllToAll"))
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
    DEBUG_ONLY(CallStackEntry cse("mpi::AllToAll"))
#ifdef EL_AVOID_COMPLEX_MPI
    int p;
    MPI_Comm_size( comm.comm, &p );
    vector<int> scsDoubled(p);
    vector<int> sdsDoubled(p);
    vector<int> rcsDoubled(p);
    vector<int> rdsDoubled(p);
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
#ifdef EL_HAVE_MPI_LONG_LONG
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
#ifdef EL_HAVE_QUAD
template void AllToAll
( const Quad* sbuf, const int* scs, const int* sds,
        Quad* rbuf, const int* rcs, const int* rds, Comm comm );
#endif
template void AllToAll
( const Complex<float>* sbuf, const int* scs, const int* sds,
        Complex<float>* rbuf, const int* rcs, const int* rds, Comm comm );
template void AllToAll
( const Complex<double>* sbuf, const int* scs, const int* sds,
        Complex<double>* rbuf, const int* rcs, const int* rds, Comm comm );
#ifdef EL_HAVE_QUAD
template void AllToAll
( const Complex<Quad>* sbuf, const int* scs, const int* sds,
        Complex<Quad>* rbuf, const int* rcs, const int* rds, Comm comm );
#endif

template<typename Real>
void Reduce
( const Real* sbuf, Real* rbuf, int count, Op op, int root, Comm comm )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::Reduce"))
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
    DEBUG_ONLY(CallStackEntry cse("mpi::Reduce"))
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

template void Reduce( const byte* sbuf, byte* rbuf, int count, Op op, int root, Comm comm );
template void Reduce( const int* sbuf, int* rbuf, int count, Op op, int root, Comm comm );
template void Reduce( const unsigned* sbuf, unsigned* rbuf, int count, Op op, int root, Comm comm );
template void Reduce( const long int* sbuf, long int* rbuf, int count, Op op, int root, Comm comm );
template void Reduce( const unsigned long* sbuf, unsigned long* rbuf, int count, Op op, int root, Comm comm );
#ifdef EL_HAVE_MPI_LONG_LONG
template void Reduce( const long long int* sbuf, long long int* rbuf, int count, Op op, int root, Comm comm );
template void Reduce( const unsigned long long* sbuf, unsigned long long* rbuf, int count, Op op, int root, Comm comm );
#endif
template void Reduce( const float* sbuf, float* rbuf, int count, Op op, int root, Comm comm );
template void Reduce( const double* sbuf, double* rbuf, int count, Op op, int root, Comm comm );
#ifdef EL_HAVE_QUAD
template void Reduce( const Quad* sbuf, Quad* rbuf, int count, Op op, int root, Comm comm );
#endif
template void Reduce( const Complex<float>* sbuf, Complex<float>* rbuf, int count, Op op, int root, Comm comm );
template void Reduce( const Complex<double>* sbuf, Complex<double>* rbuf, int count, Op op, int root, Comm comm );
#ifdef EL_HAVE_QUAD
template void Reduce( const Complex<Quad>* sbuf, Complex<Quad>* rbuf, int count, Op op, int root, Comm comm );
#endif
template void Reduce( const ValueInt<Int>* sbuf, ValueInt<Int>* rbuf, int count, Op op, int root, Comm comm );
template void Reduce( const ValueInt<float>* sbuf, ValueInt<float>* rbuf, int count, Op op, int root, Comm comm );
template void Reduce( const ValueInt<double>* sbuf, ValueInt<double>* rbuf, int count, Op op, int root, Comm comm );
#ifdef EL_HAVE_QUAD
template void Reduce( const ValueInt<Quad>* sbuf, ValueInt<Quad>* rbuf, int count, Op op, int root, Comm comm );
#endif
template void Reduce( const ValueIntPair<Int>* sbuf, ValueIntPair<Int>* rbuf, int count, Op op, int root, Comm comm );
template void Reduce( const ValueIntPair<float>* sbuf, ValueIntPair<float>* rbuf, int count, Op op, int root, Comm comm );
template void Reduce( const ValueIntPair<double>* sbuf, ValueIntPair<double>* rbuf, int count, Op op, int root, Comm comm );
#ifdef EL_HAVE_QUAD
template void Reduce( const ValueIntPair<Quad>* sbuf, ValueIntPair<Quad>* rbuf, int count, Op op, int root, Comm comm );
#endif

template<typename T>
void Reduce( const T* sbuf, T* rbuf, int count, int root, Comm comm )
{ Reduce( sbuf, rbuf, count, mpi::SUM, root, comm ); }

template void Reduce( const byte* sbuf, byte* rbuf, int count, int root, Comm comm );
template void Reduce( const int* sbuf, int* rbuf, int count, int root, Comm comm );
template void Reduce( const unsigned* sbuf, unsigned* rbuf, int count, int root, Comm comm );
template void Reduce( const long int* sbuf, long int* rbuf, int count, int root, Comm comm );
template void Reduce( const unsigned long* sbuf, unsigned long* rbuf, int count, int root, Comm comm );
#ifdef EL_HAVE_MPI_LONG_LONG
template void Reduce( const long long int* sbuf, long long int* rbuf, int count, int root, Comm comm );
template void Reduce( const unsigned long long* sbuf, unsigned long long* rbuf, int count, int root, Comm comm );
#endif
template void Reduce( const float* sbuf, float* rbuf, int count, int root, Comm comm );
template void Reduce( const double* sbuf, double* rbuf, int count, int root, Comm comm );
#ifdef EL_HAVE_QUAD
template void Reduce( const Quad* sbuf, Quad* rbuf, int count, int root, Comm comm );
#endif
template void Reduce( const Complex<float>* sbuf, Complex<float>* rbuf, int count, int root, Comm comm );
template void Reduce( const Complex<double>* sbuf, Complex<double>* rbuf, int count, int root, Comm comm );
template void Reduce( const ValueInt<Int>* sbuf, ValueInt<Int>* rbuf, int count, int root, Comm comm );
template void Reduce( const ValueInt<float>* sbuf, ValueInt<float>* rbuf, int count, int root, Comm comm );
template void Reduce( const ValueInt<double>* sbuf, ValueInt<double>* rbuf, int count, int root, Comm comm );
#ifdef EL_HAVE_QUAD
template void Reduce( const ValueInt<Quad>* sbuf, ValueInt<Quad>* rbuf, int count, int root, Comm comm );
#endif
template void Reduce( const ValueIntPair<Int>* sbuf, ValueIntPair<Int>* rbuf, int count, int root, Comm comm );
template void Reduce( const ValueIntPair<float>* sbuf, ValueIntPair<float>* rbuf, int count, int root, Comm comm );
template void Reduce( const ValueIntPair<double>* sbuf, ValueIntPair<double>* rbuf, int count, int root, Comm comm );
#ifdef EL_HAVE_QUAD
template void Reduce( const ValueIntPair<Quad>* sbuf, ValueIntPair<Quad>* rbuf, int count, int root, Comm comm );
#endif

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
#ifdef EL_HAVE_MPI_LONG_LONG
template long long int Reduce( long long int sb, Op op, int root, Comm comm );
template unsigned long long Reduce( unsigned long long sb, Op op, int root, Comm comm );
#endif
template float Reduce( float sb, Op op, int root, Comm comm );
template double Reduce( double sb, Op op, int root, Comm comm );
#ifdef EL_HAVE_QUAD
template Quad Reduce( Quad sb, Op op, int root, Comm comm );
#endif
template Complex<float> Reduce( Complex<float> sb, Op op, int root, Comm comm );
template Complex<double> Reduce( Complex<double> sb, Op op, int root, Comm comm );
#ifdef EL_HAVE_QUAD
template Complex<Quad> Reduce( Complex<Quad> sb, Op op, int root, Comm comm );
#endif
template ValueInt<Int> Reduce( ValueInt<Int> sb, Op op, int root, Comm comm );
template ValueInt<float> Reduce( ValueInt<float> sb, Op op, int root, Comm comm );
template ValueInt<double> Reduce( ValueInt<double> sb, Op op, int root, Comm comm );
#ifdef EL_HAVE_QUAD
template ValueInt<Quad> Reduce( ValueInt<Quad> sb, Op op, int root, Comm comm );
#endif
template ValueIntPair<Int> Reduce( ValueIntPair<Int> sb, Op op, int root, Comm comm );
template ValueIntPair<float> Reduce( ValueIntPair<float> sb, Op op, int root, Comm comm );
template ValueIntPair<double> Reduce( ValueIntPair<double> sb, Op op, int root, Comm comm );
#ifdef EL_HAVE_QUAD
template ValueIntPair<Quad> Reduce( ValueIntPair<Quad> sb, Op op, int root, Comm comm );
#endif

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
#ifdef EL_HAVE_MPI_LONG_LONG
template long long int Reduce( long long int sb, int root, Comm comm );
template unsigned long long Reduce( unsigned long long sb, int root, Comm comm );
#endif
template float Reduce( float sb, int root, Comm comm );
template double Reduce( double sb, int root, Comm comm );
#ifdef EL_HAVE_QUAD
template Quad Reduce( Quad sb, int root, Comm comm );
#endif
template Complex<float> Reduce( Complex<float> sb, int root, Comm comm );
template Complex<double> Reduce( Complex<double> sb, int root, Comm comm );
#ifdef EL_HAVE_QUAD
template Complex<Quad> Reduce( Complex<Quad> sb, int root, Comm comm );
#endif
template ValueInt<Int> Reduce( ValueInt<Int> sb, int root, Comm comm );
template ValueInt<float> Reduce( ValueInt<float> sb, int root, Comm comm );
template ValueInt<double> Reduce( ValueInt<double> sb, int root, Comm comm );
#ifdef EL_HAVE_QUAD
template ValueInt<Quad> Reduce( ValueInt<Quad> sb, int root, Comm comm );
#endif
template ValueIntPair<Int> Reduce( ValueIntPair<Int> sb, int root, Comm comm );
template ValueIntPair<float> Reduce( ValueIntPair<float> sb, int root, Comm comm );
template ValueIntPair<double> Reduce( ValueIntPair<double> sb, int root, Comm comm );
#ifdef EL_HAVE_QUAD
template ValueIntPair<Quad> Reduce( ValueIntPair<Quad> sb, int root, Comm comm );
#endif

template<typename Real>
void Reduce( Real* buf, int count, Op op, int root, Comm comm )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::Reduce"))
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
    DEBUG_ONLY(CallStackEntry cse("mpi::Reduce"))
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

template void Reduce( byte* buf, int count, Op op, int root, Comm comm );
template void Reduce( int* buf, int count, Op op, int root, Comm comm );
template void Reduce( unsigned* buf, int count, Op op, int root, Comm comm );
template void Reduce( long int* buf, int count, Op op, int root, Comm comm );
template void Reduce( unsigned long* buf, int count, Op op, int root, Comm comm );
#ifdef EL_HAVE_MPI_LONG_LONG
template void Reduce( long long int* buf, int count, Op op, int root, Comm comm );
template void Reduce( unsigned long long* buf, int count, Op op, int root, Comm comm );
#endif
template void Reduce( float* buf, int count, Op op, int root, Comm comm );
template void Reduce( double* buf, int count, Op op, int root, Comm comm );
#ifdef EL_HAVE_QUAD
template void Reduce( Quad* buf, int count, Op op, int root, Comm comm );
#endif
template void Reduce( Complex<float>* buf, int count, Op op, int root, Comm comm );
template void Reduce( Complex<double>* buf, int count, Op op, int root, Comm comm );
#ifdef EL_HAVE_QUAD
template void Reduce( Complex<Quad>* buf, int count, Op op, int root, Comm comm );
#endif
template void Reduce( ValueInt<Int>* buf, int count, Op op, int root, Comm comm );
template void Reduce( ValueInt<float>* buf, int count, Op op, int root, Comm comm );
template void Reduce( ValueInt<double>* buf, int count, Op op, int root, Comm comm );
#ifdef EL_HAVE_QUAD
template void Reduce( ValueInt<Quad>* buf, int count, Op op, int root, Comm comm );
#endif
template void Reduce( ValueIntPair<Int>* buf, int count, Op op, int root, Comm comm );
template void Reduce( ValueIntPair<float>* buf, int count, Op op, int root, Comm comm );
template void Reduce( ValueIntPair<double>* buf, int count, Op op, int root, Comm comm );
#ifdef EL_HAVE_QUAD
template void Reduce( ValueIntPair<Quad>* buf, int count, Op op, int root, Comm comm );
#endif

template<typename T>
void Reduce( T* buf, int count, int root, Comm comm )
{ Reduce( buf, count, mpi::SUM, root, comm ); }

template void Reduce( byte* buf, int count, int root, Comm comm );
template void Reduce( int* buf, int count, int root, Comm comm );
template void Reduce( unsigned* buf, int count, int root, Comm comm );
template void Reduce( long int* buf, int count, int root, Comm comm );
template void Reduce( unsigned long* buf, int count, int root, Comm comm );
#ifdef EL_HAVE_MPI_LONG_LONG
template void Reduce( long long int* buf, int count, int root, Comm comm );
template void Reduce( unsigned long long* buf, int count, int root, Comm comm );
#endif
template void Reduce( float* buf, int count, int root, Comm comm );
template void Reduce( double* buf, int count, int root, Comm comm );
#ifdef EL_HAVE_QUAD
template void Reduce( Quad* buf, int count, int root, Comm comm );
#endif
template void Reduce( Complex<float>* buf, int count, int root, Comm comm );
template void Reduce( Complex<double>* buf, int count, int root, Comm comm );
#ifdef EL_HAVE_QUAD
template void Reduce( Complex<Quad>* buf, int count, int root, Comm comm );
#endif
template void Reduce( ValueInt<Int>* buf, int count, int root, Comm comm );
template void Reduce( ValueInt<float>* buf, int count, int root, Comm comm );
template void Reduce( ValueInt<double>* buf, int count, int root, Comm comm );
#ifdef EL_HAVE_QUAD
template void Reduce( ValueInt<Quad>* buf, int count, int root, Comm comm );
#endif
template void Reduce( ValueIntPair<Int>* buf, int count, int root, Comm comm );
template void Reduce( ValueIntPair<float>* buf, int count, int root, Comm comm );
template void Reduce( ValueIntPair<double>* buf, int count, int root, Comm comm );
#ifdef EL_HAVE_QUAD
template void Reduce( ValueIntPair<Quad>* buf, int count, int root, Comm comm );
#endif

template<typename Real>
void AllReduce( const Real* sbuf, Real* rbuf, int count, Op op, Comm comm )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::AllReduce"))
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
    DEBUG_ONLY(CallStackEntry cse("mpi::AllReduce"))
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

template void AllReduce( const byte* sbuf, byte* rbuf, int count, Op op, Comm comm );
template void AllReduce( const int* sbuf, int* rbuf, int count, Op op, Comm comm );
template void AllReduce( const unsigned* sbuf, unsigned* rbuf, int count, Op op, Comm comm );
template void AllReduce( const long int* sbuf, long int* rbuf, int count, Op op, Comm comm );
template void AllReduce( const unsigned long* sbuf, unsigned long* rbuf, int count, Op op, Comm comm );
#ifdef EL_HAVE_MPI_LONG_LONG
template void AllReduce( const long long int* sbuf, long long int* rbuf, int count, Op op, Comm comm );
template void AllReduce( const unsigned long long* sbuf, unsigned long long* rbuf, int count, Op op, Comm comm );
#endif
template void AllReduce( const float* sbuf, float* rbuf, int count, Op op, Comm comm );
template void AllReduce( const double* sbuf, double* rbuf, int count, Op op, Comm comm );
#ifdef EL_HAVE_QUAD
template void AllReduce( const Quad* sbuf, Quad* rbuf, int count, Op op, Comm comm );
#endif
template void AllReduce( const Complex<float>* sbuf, Complex<float>* rbuf, int count, Op op, Comm comm );
template void AllReduce( const Complex<double>* sbuf, Complex<double>* rbuf, int count, Op op, Comm comm );
#ifdef EL_HAVE_QUAD
template void AllReduce<Quad>( const Complex<Quad>* sbuf, Complex<Quad>* rbuf, int count, Op op, Comm comm );
#endif
template void AllReduce( const ValueInt<Int>* sbuf, ValueInt<Int>* rbuf, int count, Op op, Comm comm );
template void AllReduce( const ValueInt<float>* sbuf, ValueInt<float>* rbuf, int count, Op op, Comm comm );
template void AllReduce( const ValueInt<double>* sbuf, ValueInt<double>* rbuf, int count, Op op, Comm comm );
#ifdef EL_HAVE_QUAD
template void AllReduce( const ValueInt<Quad>* sbuf, ValueInt<Quad>* rbuf, int count, Op op, Comm comm );
#endif
template void AllReduce( const ValueIntPair<Int>* sbuf, ValueIntPair<Int>* rbuf, int count, Op op, Comm comm );
template void AllReduce( const ValueIntPair<float>* sbuf, ValueIntPair<float>* rbuf, int count, Op op, Comm comm );
template void AllReduce( const ValueIntPair<double>* sbuf, ValueIntPair<double>* rbuf, int count, Op op, Comm comm );
#ifdef EL_HAVE_QUAD
template void AllReduce( const ValueIntPair<Quad>* sbuf, ValueIntPair<Quad>* rbuf, int count, Op op, Comm comm );
#endif

template<typename T>
void AllReduce( const T* sbuf, T* rbuf, int count, Comm comm )
{ AllReduce( sbuf, rbuf, count, mpi::SUM, comm ); }

template void AllReduce( const byte* sbuf, byte* rbuf, int count, Comm comm );
template void AllReduce( const int* sbuf, int* rbuf, int count, Comm comm );
template void AllReduce( const unsigned* sbuf, unsigned* rbuf, int count, Comm comm );
template void AllReduce( const long int* sbuf, long int* rbuf, int count, Comm comm );
template void AllReduce( const unsigned long* sbuf, unsigned long* rbuf, int count, Comm comm );
#ifdef EL_HAVE_MPI_LONG_LONG
template void AllReduce( const long long int* sbuf, long long int* rbuf, int count, Comm comm );
template void AllReduce( const unsigned long long* sbuf, unsigned long long* rbuf, int count, Comm comm );
#endif
template void AllReduce( const float* sbuf, float* rbuf, int count, Comm comm );
template void AllReduce( const double* sbuf, double* rbuf, int count, Comm comm );
#ifdef EL_HAVE_QUAD
template void AllReduce( const Quad* sbuf, Quad* rbuf, int count, Comm comm );
#endif
template void AllReduce( const Complex<float>* sbuf, Complex<float>* rbuf, int count, Comm comm );
template void AllReduce( const Complex<double>* sbuf, Complex<double>* rbuf, int count, Comm comm );
#ifdef EL_HAVE_QUAD
template void AllReduce( const Complex<Quad>* sbuf, Complex<Quad>* rbuf, int count, Comm comm );
#endif
template void AllReduce( const ValueInt<Int>* sbuf, ValueInt<Int>* rbuf, int count, Comm comm );
template void AllReduce( const ValueInt<float>* sbuf, ValueInt<float>* rbuf, int count, Comm comm );
template void AllReduce( const ValueInt<double>* sbuf, ValueInt<double>* rbuf, int count, Comm comm );
#ifdef EL_HAVE_QUAD
template void AllReduce( const ValueInt<Quad>* sbuf, ValueInt<Quad>* rbuf, int count, Comm comm );
#endif
template void AllReduce( const ValueIntPair<Int>* sbuf, ValueIntPair<Int>* rbuf, int count, Comm comm );
template void AllReduce( const ValueIntPair<float>* sbuf, ValueIntPair<float>* rbuf, int count, Comm comm );
template void AllReduce( const ValueIntPair<double>* sbuf, ValueIntPair<double>* rbuf, int count, Comm comm );
#ifdef EL_HAVE_QUAD
template void AllReduce( const ValueIntPair<Quad>* sbuf, ValueIntPair<Quad>* rbuf, int count, Comm comm );
#endif

template<typename T>
T AllReduce( T sb, Op op, Comm comm )
{ T rb; AllReduce( &sb, &rb, 1, op, comm ); return rb; }

template byte AllReduce( byte sb, Op op, Comm comm );
template int AllReduce( int sb, Op op, Comm comm );
template unsigned AllReduce( unsigned sb, Op op, Comm comm );
template long int AllReduce( long int sb, Op op, Comm comm );
template unsigned long AllReduce( unsigned long sb, Op op, Comm comm );
#ifdef EL_HAVE_MPI_LONG_LONG
template long long int AllReduce( long long int sb, Op op, Comm comm );
template unsigned long long AllReduce( unsigned long long sb, Op op, Comm comm );
#endif
template float AllReduce( float sb, Op op, Comm comm );
template double AllReduce( double sb, Op op, Comm comm );
#ifdef EL_HAVE_QUAD
template Quad AllReduce( Quad sb, Op op, Comm comm );
#endif
template Complex<float> AllReduce( Complex<float> sb, Op op, Comm comm );
template Complex<double> AllReduce( Complex<double> sb, Op op, Comm comm );
#ifdef EL_HAVE_QUAD
template Complex<Quad> AllReduce( Complex<Quad> sb, Op op, Comm comm );
#endif
template ValueInt<Int> AllReduce( ValueInt<Int> sb, Op op, Comm comm );
template ValueInt<float> AllReduce( ValueInt<float> sb, Op op, Comm comm );
template ValueInt<double> AllReduce( ValueInt<double> sb, Op op, Comm comm );
#ifdef EL_HAVE_QUAD
template ValueInt<Quad> AllReduce( ValueInt<Quad> sb, Op op, Comm comm );
#endif
template ValueIntPair<Int> AllReduce( ValueIntPair<Int> sb, Op op, Comm comm );
template ValueIntPair<float> AllReduce( ValueIntPair<float> sb, Op op, Comm comm );
template ValueIntPair<double> AllReduce( ValueIntPair<double> sb, Op op, Comm comm );
#ifdef EL_HAVE_QUAD
template ValueIntPair<Quad> AllReduce( ValueIntPair<Quad> sb, Op op, Comm comm );
#endif

template<typename T>
T AllReduce( T sb, Comm comm )
{ return AllReduce( sb, mpi::SUM, comm ); }

template byte AllReduce( byte sb, Comm comm );
template int AllReduce( int sb, Comm comm );
template unsigned AllReduce( unsigned sb, Comm comm );
template long int AllReduce( long int sb, Comm comm );
template unsigned long AllReduce( unsigned long sb, Comm comm );
#ifdef EL_HAVE_MPI_LONG_LONG
template long long int AllReduce( long long int sb, Comm comm );
template unsigned long long AllReduce( unsigned long long sb, Comm comm );
#endif
template float AllReduce( float sb, Comm comm );
template double AllReduce( double sb, Comm comm );
#ifdef EL_HAVE_QUAD
template Quad AllReduce( Quad sb, Comm comm );
#endif
template Complex<float> AllReduce( Complex<float> sb, Comm comm );
template Complex<double> AllReduce( Complex<double> sb, Comm comm );
#ifdef EL_HAVE_QUAD
template Complex<Quad> AllReduce( Complex<Quad> sb, Comm comm );
#endif
template ValueInt<Int> AllReduce( ValueInt<Int> sb, Comm comm );
template ValueInt<float> AllReduce( ValueInt<float> sb, Comm comm );
template ValueInt<double> AllReduce( ValueInt<double> sb, Comm comm );
#ifdef EL_HAVE_QUAD
template ValueInt<Quad> AllReduce( ValueInt<Quad> sb, Comm comm );
#endif
template ValueIntPair<Int> AllReduce( ValueIntPair<Int> sb, Comm comm );
template ValueIntPair<float> AllReduce( ValueIntPair<float> sb, Comm comm );
template ValueIntPair<double> AllReduce( ValueIntPair<double> sb, Comm comm );
#ifdef EL_HAVE_QUAD
template ValueIntPair<Quad> AllReduce( ValueIntPair<Quad> sb, Comm comm );
#endif

template<typename Real>
void AllReduce( Real* buf, int count, Op op, Comm comm )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::AllReduce"))
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
    DEBUG_ONLY(CallStackEntry cse("mpi::AllReduce"))
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

template void AllReduce( byte* buf, int count, Op op, Comm comm );
template void AllReduce( int* buf, int count, Op op, Comm comm );
template void AllReduce( unsigned* buf, int count, Op op, Comm comm );
template void AllReduce( long int* buf, int count, Op op, Comm comm );
template void AllReduce( unsigned long* buf, int count, Op op, Comm comm );
#ifdef EL_HAVE_MPI_LONG_LONG
template void AllReduce( long long int* buf, int count, Op op, Comm comm );
template void AllReduce( unsigned long long* buf, int count, Op op, Comm comm );
#endif
template void AllReduce( float* buf, int count, Op op, Comm comm );
template void AllReduce( double* buf, int count, Op op, Comm comm );
#ifdef EL_HAVE_QUAD
template void AllReduce( Quad* buf, int count, Op op, Comm comm );
#endif
template void AllReduce( Complex<float>* buf, int count, Op op, Comm comm );
template void AllReduce( Complex<double>* buf, int count, Op op, Comm comm );
#ifdef EL_HAVE_QUAD
template void AllReduce( Complex<Quad>* buf, int count, Op op, Comm comm );
#endif
template void AllReduce( ValueInt<Int>* buf, int count, Op op, Comm comm );
template void AllReduce( ValueInt<float>* buf, int count, Op op, Comm comm );
template void AllReduce( ValueInt<double>* buf, int count, Op op, Comm comm );
#ifdef EL_HAVE_QUAD
template void AllReduce( ValueInt<Quad>* buf, int count, Op op, Comm comm );
#endif
template void AllReduce( ValueIntPair<Int>* buf, int count, Op op, Comm comm );
template void AllReduce( ValueIntPair<float>* buf, int count, Op op, Comm comm );
template void AllReduce( ValueIntPair<double>* buf, int count, Op op, Comm comm );
#ifdef EL_HAVE_QUAD
template void AllReduce( ValueIntPair<Quad>* buf, int count, Op op, Comm comm );
#endif

template<typename T>
void AllReduce( T* buf, int count, Comm comm )
{ AllReduce( buf, count, mpi::SUM, comm ); }

template void AllReduce( byte* buf, int count, Comm comm );
template void AllReduce( int* buf, int count, Comm comm );
template void AllReduce( unsigned* buf, int count, Comm comm );
template void AllReduce( long int* buf, int count, Comm comm );
template void AllReduce( unsigned long* buf, int count, Comm comm );
#ifdef EL_HAVE_MPI_LONG_LONG
template void AllReduce( long long int* buf, int count, Comm comm );
template void AllReduce( unsigned long long* buf, int count, Comm comm );
#endif
template void AllReduce( float* buf, int count, Comm comm );
template void AllReduce( double* buf, int count, Comm comm );
#ifdef EL_HAVE_QUAD
template void AllReduce( Quad* buf, int count, Comm comm );
#endif
template void AllReduce( Complex<float>* buf, int count, Comm comm );
template void AllReduce( Complex<double>* buf, int count, Comm comm );
#ifdef EL_HAVE_QUAD
template void AllReduce( Complex<Quad>* buf, int count, Comm comm );
#endif
template void AllReduce( ValueInt<Int>* buf, int count, Comm comm );
template void AllReduce( ValueInt<float>* buf, int count, Comm comm );
template void AllReduce( ValueInt<double>* buf, int count, Comm comm );
#ifdef EL_HAVE_QUAD
template void AllReduce( ValueInt<Quad>* buf, int count, Comm comm );
#endif
template void AllReduce( ValueIntPair<Int>* buf, int count, Comm comm );
template void AllReduce( ValueIntPair<float>* buf, int count, Comm comm );
template void AllReduce( ValueIntPair<double>* buf, int count, Comm comm );
#ifdef EL_HAVE_QUAD
template void AllReduce( ValueIntPair<Quad>* buf, int count, Comm comm );
#endif

template<typename Real>
void ReduceScatter( Real* sbuf, Real* rbuf, int rc, Op op, Comm comm )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::ReduceScatter"))
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
    DEBUG_ONLY(CallStackEntry cse("mpi::ReduceScatter"))
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

template void ReduceScatter( byte* sbuf, byte* rbuf, int rc, Op op, Comm comm );
template void ReduceScatter( int* sbuf, int* rbuf, int rc, Op op, Comm comm );
template void ReduceScatter( unsigned* sbuf, unsigned* rbuf, int rc, Op op, Comm comm );
template void ReduceScatter( long int* sbuf, long int* rbuf, int rc, Op op, Comm comm );
template void ReduceScatter( unsigned long* sbuf, unsigned long* rbuf, int rc, Op op, Comm comm );
#ifdef EL_HAVE_MPI_LONG_LONG
template void ReduceScatter( long long int* sbuf, long long int* rbuf, int rc, Op op, Comm comm );
template void ReduceScatter( unsigned long long* sbuf, unsigned long long* rbuf, int rc, Op op, Comm comm );
#endif
template void ReduceScatter( float* sbuf, float* rbuf, int rc, Op op, Comm comm );
template void ReduceScatter( double* sbuf, double* rbuf, int rc, Op op, Comm comm );
#ifdef EL_HAVE_QUAD
template void ReduceScatter( Quad* sbuf, Quad* rbuf, int rc, Op op, Comm comm );
#endif
template void ReduceScatter( Complex<float>* sbuf, Complex<float>* rbuf, int rc, Op op, Comm comm );
template void ReduceScatter( Complex<double>* sbuf, Complex<double>* rbuf, int rc, Op op, Comm comm );
#ifdef EL_HAVE_QUAD
template void ReduceScatter( Complex<Quad>* sbuf, Complex<Quad>* rbuf, int rc, Op op, Comm comm );
#endif

template<typename T>
void ReduceScatter( T* sbuf, T* rbuf, int rc, Comm comm )
{ ReduceScatter( sbuf, rbuf, rc, mpi::SUM, comm ); }

template void ReduceScatter( byte* sbuf, byte* rbuf, int rc, Comm comm );
template void ReduceScatter( int* sbuf, int* rbuf, int rc, Comm comm );
template void ReduceScatter( unsigned* sbuf, unsigned* rbuf, int rc, Comm comm );
template void ReduceScatter( long int* sbuf, long int* rbuf, int rc, Comm comm );
template void ReduceScatter( unsigned long* sbuf, unsigned long* rbuf, int rc, Comm comm );
#ifdef EL_HAVE_MPI_LONG_LONG
template void ReduceScatter( long long int* sbuf, long long int* rbuf, int rc, Comm comm );
template void ReduceScatter( unsigned long long* sbuf, unsigned long long* rbuf, int rc, Comm comm );
#endif
template void ReduceScatter( float* sbuf, float* rbuf, int rc, Comm comm );
template void ReduceScatter( double* sbuf, double* rbuf, int rc, Comm comm );
#ifdef EL_HAVE_QUAD
template void ReduceScatter( Quad* sbuf, Quad* rbuf, int rc, Comm comm );
#endif
template void ReduceScatter( Complex<float>* sbuf, Complex<float>* rbuf, int rc, Comm comm );
template void ReduceScatter( Complex<double>* sbuf, Complex<double>* rbuf, int rc, Comm comm );
#ifdef EL_HAVE_QUAD
template void ReduceScatter( Complex<Quad>* sbuf, Complex<Quad>* rbuf, int rc, Comm comm );
#endif

template<typename T>
T ReduceScatter( T sb, Op op, Comm comm )
{ T rb; ReduceScatter( &sb, &rb, 1, op, comm ); return rb; }

template byte ReduceScatter( byte sb, Op op, Comm comm );
template int ReduceScatter( int sb, Op op, Comm comm );
template unsigned ReduceScatter( unsigned sb, Op op, Comm comm );
template long int ReduceScatter( long int sb, Op op, Comm comm );
template unsigned long ReduceScatter( unsigned long sb, Op op, Comm comm );
#ifdef EL_HAVE_MPI_LONG_LONG
template long long int ReduceScatter( long long int sb, Op op, Comm comm );
template unsigned long long ReduceScatter( unsigned long long sb, Op op, Comm comm );
#endif
template float ReduceScatter( float sb, Op op, Comm comm );
template double ReduceScatter( double sb, Op op, Comm comm );
#ifdef EL_HAVE_QUAD
template Quad ReduceScatter( Quad sb, Op op, Comm comm );
#endif
template Complex<float> ReduceScatter( Complex<float> sb, Op op, Comm comm );
template Complex<double> ReduceScatter( Complex<double> sb, Op op, Comm comm );
#ifdef EL_HAVE_QUAD
template Complex<Quad> ReduceScatter( Complex<Quad> sb, Op op, Comm comm );
#endif

template<typename T>
T ReduceScatter( T sb, Comm comm )
{ return ReduceScatter( sb, mpi::SUM, comm ); }

template byte ReduceScatter( byte sb, Comm comm );
template int ReduceScatter( int sb, Comm comm );
template unsigned ReduceScatter( unsigned sb, Comm comm );
template long int ReduceScatter( long int sb, Comm comm );
template unsigned long ReduceScatter( unsigned long sb, Comm comm );
#ifdef EL_HAVE_MPI_LONG_LONG
template long long int ReduceScatter( long long int sb, Comm comm );
template unsigned long long ReduceScatter( unsigned long long sb, Comm comm );
#endif
template float ReduceScatter( float sb, Comm comm );
template double ReduceScatter( double sb, Comm comm );
#ifdef EL_HAVE_QUAD
template Quad ReduceScatter( Quad sb, Comm comm );
#endif
template Complex<float> ReduceScatter( Complex<float> sb, Comm comm );
template Complex<double> ReduceScatter( Complex<double> sb, Comm comm );
#ifdef EL_HAVE_QUAD
template Complex<Quad> ReduceScatter( Complex<Quad> sb, Comm comm );
#endif

template<typename Real>
void ReduceScatter( Real* buf, int rc, Op op, Comm comm )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::ReduceScatter"))
    MPI_Op opC;
    if( op == SUM )
        opC = SumOp<Real>().op; 
    else if( op == MAX )
        opC = MaxOp<Real>().op;
    else if( op == MIN )
        opC = MinOp<Real>().op;
    else
        opC = op.op;

#ifdef EL_REDUCE_SCATTER_BLOCK_VIA_ALLREDUCE
    const int commSize = Size( comm );
    const int commRank = Rank( comm );
    AllReduce( buf, rc*commSize, op, comm );
    if( commRank != 0 )
        MemCopy( buf, &buf[commRank*rc], rc );
#elif defined(EL_HAVE_MPI_REDUCE_SCATTER_BLOCK)
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
    DEBUG_ONLY(CallStackEntry cse("mpi::ReduceScatter"))
    MPI_Op opC;
    if( op == SUM )
        opC = SumOp<Complex<Real>>().op; 
    else
        opC = op.op;

#ifdef EL_REDUCE_SCATTER_BLOCK_VIA_ALLREDUCE
    const int commSize = Size( comm );
    const int commRank = Rank( comm );
    AllReduce( buf, rc*commSize, op, comm );
    if( commRank != 0 )
        MemCopy( buf, &buf[commRank*rc], rc );
#elif defined(EL_HAVE_MPI_REDUCE_SCATTER_BLOCK)
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

template void ReduceScatter( byte* buf, int rc, Op op, Comm comm );
template void ReduceScatter( int* buf, int rc, Op op, Comm comm );
template void ReduceScatter( unsigned* buf, int rc, Op op, Comm comm );
template void ReduceScatter( long int* buf, int rc, Op op, Comm comm );
template void ReduceScatter( unsigned long* buf, int rc, Op op, Comm comm );
#ifdef EL_HAVE_MPI_LONG_LONG
template void ReduceScatter( long long int* buf, int rc, Op op, Comm comm );
template void ReduceScatter( unsigned long long* buf, int rc, Op op, Comm comm );
#endif
template void ReduceScatter( float* buf, int rc, Op op, Comm comm );
template void ReduceScatter( double* buf, int rc, Op op, Comm comm );
#ifdef EL_HAVE_QUAD
template void ReduceScatter( Quad* buf, int rc, Op op, Comm comm );
#endif
template void ReduceScatter( Complex<float>* buf, int rc, Op op, Comm comm );
template void ReduceScatter( Complex<double>* buf, int rc, Op op, Comm comm );
#ifdef EL_HAVE_QUAD
template void ReduceScatter( Complex<Quad>* buf, int rc, Op op, Comm comm );
#endif

template<typename T>
void ReduceScatter( T* buf, int rc, Comm comm )
{ ReduceScatter( buf, rc, mpi::SUM, comm ); }

template void ReduceScatter( byte* buf, int rc, Comm comm );
template void ReduceScatter( int* buf, int rc, Comm comm );
template void ReduceScatter( unsigned* buf, int rc, Comm comm );
template void ReduceScatter( long int* buf, int rc, Comm comm );
template void ReduceScatter( unsigned long* buf, int rc, Comm comm );
#ifdef EL_HAVE_MPI_LONG_LONG
template void ReduceScatter( long long int* buf, int rc, Comm comm );
template void ReduceScatter( unsigned long long* buf, int rc, Comm comm );
#endif
template void ReduceScatter( float* buf, int rc, Comm comm );
template void ReduceScatter( double* buf, int rc, Comm comm );
#ifdef EL_HAVE_QUAD
template void ReduceScatter( Quad* buf, int rc, Comm comm );
#endif
template void ReduceScatter( Complex<float>* buf, int rc, Comm comm );
template void ReduceScatter( Complex<double>* buf, int rc, Comm comm );
#ifdef EL_HAVE_QUAD
template void ReduceScatter( Complex<Quad>* buf, int rc, Comm comm );
#endif

template<typename Real>
void ReduceScatter
( const Real* sbuf, Real* rbuf, const int* rcs, Op op, Comm comm )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::ReduceScatter"))
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
    DEBUG_ONLY(CallStackEntry cse("mpi::ReduceScatter"))
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

template void ReduceScatter( const byte* sbuf, byte* rbuf, const int* rcs, Op op, Comm comm );
template void ReduceScatter( const int* sbuf, int* rbuf, const int* rcs, Op op, Comm comm );
template void ReduceScatter( const unsigned* sbuf, unsigned* rbuf, const int* rcs, Op op, Comm comm );
template void ReduceScatter( const long int* sbuf, long int* rbuf, const int* rcs, Op op, Comm comm );
template void ReduceScatter( const unsigned long* sbuf, unsigned long* rbuf, const int* rcs, Op op, Comm comm );
#ifdef EL_HAVE_MPI_LONG_LONG
template void ReduceScatter( const long long int* sbuf, long long int* rbuf, const int* rcs, Op op, Comm comm );
template void ReduceScatter( const unsigned long long* sbuf, unsigned long long* rbuf, const int* rcs, Op op, Comm comm );
#endif
template void ReduceScatter( const float* sbuf, float* rbuf, const int* rcs, Op op, Comm comm );
template void ReduceScatter( const double* sbuf, double* rbuf, const int* rcs, Op op, Comm comm );
#ifdef EL_HAVE_QUAD
template void ReduceScatter( const Quad* sbuf, Quad* rbuf, const int* rcs, Op op, Comm comm );
#endif
template void ReduceScatter( const Complex<float>* sbuf, Complex<float>* rbuf, const int* rcs, Op op, Comm comm );
template void ReduceScatter( const Complex<double>* sbuf, Complex<double>* rbuf, const int* rcs, Op op, Comm comm );
#ifdef EL_HAVE_QUAD
template void ReduceScatter( const Complex<Quad>* sbuf, Complex<Quad>* rbuf, const int* rcs, Op op, Comm comm );
#endif

template<typename T>
void ReduceScatter( const T* sbuf, T* rbuf, const int* rcs, Comm comm )
{ ReduceScatter( sbuf, rbuf, rcs, mpi::SUM, comm ); }

template void ReduceScatter( const byte* sbuf, byte* rbuf, const int* rcs, Comm comm );
template void ReduceScatter( const int* sbuf, int* rbuf, const int* rcs, Comm comm );
template void ReduceScatter( const unsigned* sbuf, unsigned* rbuf, const int* rcs, Comm comm );
template void ReduceScatter( const long int* sbuf, long int* rbuf, const int* rcs, Comm comm );
template void ReduceScatter( const unsigned long* sbuf, unsigned long* rbuf, const int* rcs, Comm comm );
#ifdef EL_HAVE_MPI_LONG_LONG
template void ReduceScatter( const long long int* sbuf, long long int* rbuf, const int* rcs, Comm comm );
template void ReduceScatter( const unsigned long long* sbuf, unsigned long long* rbuf, const int* rcs, Comm comm );
#endif
template void ReduceScatter( const float* sbuf, float* rbuf, const int* rcs, Comm comm );
template void ReduceScatter( const double* sbuf, double* rbuf, const int* rcs, Comm comm );
#ifdef EL_HAVE_QUAD
template void ReduceScatter( const Quad* sbuf, Quad* rbuf, const int* rcs, Comm comm );
#endif
template void ReduceScatter( const Complex<float>* sbuf, Complex<float>* rbuf, const int* rcs, Comm comm );
template void ReduceScatter( const Complex<double>* sbuf, Complex<double>* rbuf, const int* rcs, Comm comm );
#ifdef EL_HAVE_QUAD
template void ReduceScatter( const Complex<Quad>* sbuf, Complex<Quad>* rbuf, const int* rcs, Comm comm );
#endif

void VerifySendsAndRecvs
( const vector<int>& sendCounts,
  const vector<int>& recvCounts, mpi::Comm comm )
{
    DEBUG_ONLY(CallStackEntry cse("mpi::VerifySendsAndRecvs"))
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
