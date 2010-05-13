/*
   Copyright 2009-2010 Jack Poulson

   This file is part of Elemental.

   Elemental is free software: you can redistribute it and/or modify it under
   the terms of the GNU Lesser General Public License as published by the
   Free Software Foundation; either version 3 of the License, or 
   (at your option) any later version.

   Elemental is distributed in the hope that it will be useful, but 
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with Elemental. If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef ELEMENTAL_WRAPPERS_MPI_HPP
#define ELEMENTAL_WRAPPERS_MPI_HPP 1

#include "Elemental/Environment.hpp"

namespace Elemental {
namespace wrappers {
namespace MPI {

// Minimum number of elements that each process must contribute
// in collective communications
const int MinCollectContrib = 1;

double 
Time();

void 
Barrier( MPI_Comm comm );

void
Send
( int* buf, int count, int to, int tag, MPI_Comm comm );

void
Send
( float* buf, int count, int to, int tag, MPI_Comm comm );

void
Send
( double* buf, int count, int to, int tag, MPI_Comm comm );

#ifndef WITHOUT_COMPLEX
void
Send
( scomplex* buf, int count, int to, int tag, MPI_Comm comm );

void
Send
( dcomplex* buf, int count, int to, int tag, MPI_Comm comm );
#endif
        
void
Recv
( int* buf, int count, int from, int tag, MPI_Comm comm );

void
Recv
( float* buf, int count, int from, int tag, MPI_Comm comm );

void
Recv
( double* buf, int count, int from, int tag, MPI_Comm comm );

#ifndef WITHOUT_COMPLEX
void
Recv
( scomplex* buf, int count, int from, int tag, MPI_Comm comm );

void
Recv
( dcomplex* buf, int count, int from, int tag, MPI_Comm comm );
#endif

void
SendRecv
( int* sbuf, int sc, int to,   int stag,
  int* rbuf, int rc, int from, int rtag, MPI_Comm comm );

void
SendRecv
( float* sbuf, int sc, int to,   int stag,
  float* rbuf, int rc, int from, int rtag, MPI_Comm comm );

void
SendRecv
( double* sbuf, int sc, int to,   int stag,
  double* rbuf, int rc, int from, int rtag, MPI_Comm comm );

#ifndef WITHOUT_COMPLEX
void
SendRecv
( scomplex* sbuf, int sc, int to,   int stag,
  scomplex* rbuf, int rc, int from, int rtag, MPI_Comm comm );

void
SendRecv
( dcomplex* sbuf, int sc, int to,   int stag,
  dcomplex* rbuf, int rc, int from, int rtag, MPI_Comm comm );
#endif

void
Broadcast( int* buf, int count, int root, MPI_Comm comm );

void
Broadcast( float* buf, int count, int root, MPI_Comm comm );

void
Broadcast( double* buf, int count, int root, MPI_Comm comm );

#ifndef WITHOUT_COMPLEX
void
Broadcast
( scomplex* buf, int count, int root, MPI_Comm comm );

void
Broadcast
( dcomplex* buf, int count, int root, MPI_Comm comm );
#endif

void
Gather
( int* sbuf, int sc,
  int* rbuf, int rc, int root, MPI_Comm comm );

void
Gather
( float* sbuf, int sc,
  float* rbuf, int rc, int root, MPI_Comm comm );

void
Gather
( double* sbuf, int sc,
  double* rbuf, int rc, int root, MPI_Comm comm );

#ifndef WITHOUT_COMPLEX
void 
Gather
( scomplex* sbuf, int sc,
  scomplex* rbuf, int rc, int root, MPI_Comm comm );

void
Gather
( dcomplex* sbuf, int sc,
  dcomplex* rbuf, int rc, int root, MPI_Comm comm );
#endif
    
void
AllGather
( int* sbuf, int sc,
  int* rbuf, int rc, MPI_Comm comm );

void
AllGather
( float* sbuf, int sc,
  float* rbuf, int rc, MPI_Comm comm );

void
AllGather
( double* sbuf, int sc,
  double* rbuf, int rc, MPI_Comm comm );

#ifndef WITHOUT_COMPLEX
void
AllGather
( scomplex* sbuf, int sc,
  scomplex* rbuf, int rc, MPI_Comm comm );

void
AllGather
( dcomplex* sbuf, int sc,
  dcomplex* rbuf, int rc, MPI_Comm comm );
#endif

void
Scatter
( int* sbuf, int sc,
  int* rbuf, int rc, int root, MPI_Comm comm );

void
Scatter
( float* sbuf, int sc,
  float* rbuf, int rc, int root, MPI_Comm comm );

void
Scatter
( double* sbuf, int sc,
  double* rbuf, int rc, int root, MPI_Comm comm );

#ifndef WITHOUT_COMPLEX
void
Scatter
( scomplex* sbuf, int sc,
  scomplex* rbuf, int rc, int root, MPI_Comm comm );

void
Scatter
( dcomplex* sbuf, int sc,
  dcomplex* rbuf, int rc, int root, MPI_Comm comm );
#endif
    
void
AllToAll
( int* sbuf, int sc,
  int* rbuf, int rc, MPI_Comm comm );

void
AllToAll
( float* sbuf, int sc,
  float* rbuf, int rc, MPI_Comm comm );

void
AllToAll
( double* sbuf, int sc,
  double* rbuf, int rc, MPI_Comm comm );

#ifndef WITHOUT_COMPLEX
void
AllToAll
( scomplex* sbuf, int sc,
  scomplex* rbuf, int rc, MPI_Comm comm );

void
AllToAll
( dcomplex* sbuf, int sc,
  dcomplex* rbuf, int rc, MPI_Comm comm );
#endif

void
AllToAll
( int* sbuf, int* scs, int* sds,
  int* rbuf, int* rcs, int* rds, MPI_Comm comm );

void
AllToAll
( float* sbuf, int* scs, int* sds,
  float* rbuf, int* rcs, int* rds, MPI_Comm comm );

void
AllToAll
( double* sbuf, int* scs, int* sds,
  double* rbuf, int* rcs, int* rds, MPI_Comm comm );

#ifndef WITHOUT_COMPLEX
void
AllToAll
( scomplex* sbuf, int* scs, int* sds,
  scomplex* rbuf, int* rcs, int* rds, MPI_Comm comm );

void
AllToAll
( dcomplex* sbuf, int* scs, int* sds,
  dcomplex* rbuf, int* rcs, int* rds, MPI_Comm comm );
#endif

void
Reduce
( int* sbuf, int* rbuf, int count, MPI_Op op, int root, 
  MPI_Comm comm );

void
Reduce
( float* sbuf, float* rbuf, int count, MPI_Op op, int root,
  MPI_Comm comm );

void
Reduce
( double* sbuf, double* rbuf, int count, MPI_Op op, int root,
  MPI_Comm comm );

#ifndef WITHOUT_COMPLEX
void
Reduce
( scomplex* sbuf, scomplex* rbuf, int count, MPI_Op op, int root, 
  MPI_Comm comm );

void
Reduce
( dcomplex* sbuf, dcomplex* rbuf, int count, MPI_Op op, int root, 
  MPI_Comm comm );
#endif
    
void
AllReduce
( char* sbuf, char* rbuf, int count, MPI_Op op, MPI_Comm comm );

void
AllReduce
( int* sbuf, int* rbuf, int count, MPI_Op op, MPI_Comm comm );

void
AllReduce
( float* sbuf, float* rbuf, int count, MPI_Op op, MPI_Comm comm );

void
AllReduce
( double* sbuf, double* rbuf, int count, MPI_Op op, MPI_Comm comm );
#ifndef WITHOUT_COMPLEX
void
AllReduce
( scomplex* sbuf, scomplex* rbuf, int count, MPI_Op op, MPI_Comm comm );

void
AllReduce
( dcomplex* sbuf, dcomplex* rbuf, int count, MPI_Op op, MPI_Comm comm );
#endif
        
void
ReduceScatter
( int* sbuf, int* rbuf, int* rcs, MPI_Op op, MPI_Comm comm );

void
ReduceScatter
( float* sbuf, float* rbuf, int* rcs, MPI_Op op, MPI_Comm comm );

void
ReduceScatter
( double* sbuf, double* rbuf, int* rcs, MPI_Op op, MPI_Comm comm );
#ifndef WITHOUT_COMPLEX
void
ReduceScatter
( scomplex* sbuf, scomplex* rbuf, int* rcs, MPI_Op op, MPI_Comm com );

void
ReduceScatter
( dcomplex* sbuf, dcomplex* rbuf, int* rcs, MPI_Op op, MPI_Comm comm );
#endif

} // MPI
} // wrapppers
} // Elemental

//----------------------------------------------------------------------------//
// Implementation begins here                                                 //
//----------------------------------------------------------------------------//

#ifndef RELEASE
#define SAFE_MPI( routine ) \
{ \
    int retval = routine; \
    if( retval != MPI_SUCCESS ) \
    { \
        std::cout << "MPI return value = " << retval << std::endl; \
        char errorString[200]; \
        int lengthOfErrorString; \
        MPI_Error_string( retval, errorString, &lengthOfErrorString ); \
        std::cout << "Error msg: " << errorString << std::endl; \
        DumpCallStack(); \
        throw std::exception(); \
    } \
}
#else
#define SAFE_MPI( routine ) routine
#endif

inline double
Elemental::wrappers::MPI::Time()
{ return MPI_Wtime(); }

inline void
Elemental::wrappers::MPI::Barrier( MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::Barrier");
#endif
    SAFE_MPI( MPI_Barrier( comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
Elemental::wrappers::MPI::Send
( int* buf, int count, int to, int tag, MPI_Comm comm )
{ 
#ifndef RELEASE
    PushCallStack("wrappers::MPI::Send");
#endif
    SAFE_MPI( MPI_Send( buf, count, MPI_INT, to, tag, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
Elemental::wrappers::MPI::Send
( float* buf, int count, int to, int tag, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::Send");
#endif
    SAFE_MPI( MPI_Send( buf, count, MPI_FLOAT, to, tag, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
Elemental::wrappers::MPI::Send
( double* buf, int count, int to, int tag, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::Send");
#endif
    SAFE_MPI( MPI_Send( buf, count, MPI_DOUBLE, to, tag, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
inline void
Elemental::wrappers::MPI::Send
( scomplex* buf, int count, int to, int tag, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::Send");
#endif
    SAFE_MPI( MPI_Send( buf, count, MPI_COMPLEX, to, tag, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
Elemental::wrappers::MPI::Send
( dcomplex* buf, int count, int to, int tag, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::Send");
#endif
    SAFE_MPI( MPI_Send( buf, count, MPI_DOUBLE_COMPLEX, to, tag, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif

inline void
Elemental::wrappers::MPI::Recv
( int* buf, int count, int from, int tag, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::Recv");
#endif
    MPI_Status status;
    SAFE_MPI( MPI_Recv( buf, count, MPI_INT, from, tag, comm, &status ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
Elemental::wrappers::MPI::Recv
( float* buf, int count, int from, int tag, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::Recv");
#endif
    MPI_Status status;
    SAFE_MPI( MPI_Recv( buf, count, MPI_FLOAT, from, tag, comm, &status ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
Elemental::wrappers::MPI::Recv
( double* buf, int count, int from, int tag, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::Recv");
#endif
    MPI_Status status;
    SAFE_MPI( MPI_Recv( buf, count, MPI_DOUBLE, from, tag, comm, &status ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
inline void
Elemental::wrappers::MPI::Recv
( scomplex* buf, int count, int from, int tag, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::Recv");
#endif
    MPI_Status status;
    SAFE_MPI( MPI_Recv( buf, count, MPI_COMPLEX, from, tag, comm, &status ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
Elemental::wrappers::MPI::Recv
( dcomplex* buf, int count, int from, int tag, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::Recv");
#endif
    MPI_Status status;
    SAFE_MPI( MPI_Recv( buf, count, MPI_DOUBLE_COMPLEX, from, tag, comm, 
                        &status ) );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif

inline void
Elemental::wrappers::MPI::SendRecv
( int* sbuf, int sc, int to,   int stag,
  int* rbuf, int rc, int from, int rtag, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::SendRecv");
#endif
    MPI_Status status;
    SAFE_MPI( MPI_Sendrecv( sbuf, sc, MPI_INT, to,   stag,
                            rbuf, rc, MPI_INT, from, rtag, comm, &status ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
Elemental::wrappers::MPI::SendRecv
( float* sbuf, int sc, int to,   int stag,
  float* rbuf, int rc, int from, int rtag, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::SendRecv");
#endif
    MPI_Status status;
    SAFE_MPI( MPI_Sendrecv( sbuf, sc, MPI_FLOAT, to,   stag,
                            rbuf, rc, MPI_FLOAT, from, rtag, comm, &status ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
Elemental::wrappers::MPI::SendRecv
( double* sbuf, int sc, int to,   int stag,
  double* rbuf, int rc, int from, int rtag, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::SendRecv");
#endif
    MPI_Status status;
    SAFE_MPI( MPI_Sendrecv( sbuf, sc, MPI_DOUBLE, to,   stag,
                            rbuf, rc, MPI_DOUBLE, from, rtag, comm, &status ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
inline void
Elemental::wrappers::MPI::SendRecv
( scomplex* sbuf, int sc, int to,   int stag,
  scomplex* rbuf, int rc, int from, int rtag, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::SendRecv");
#endif
    MPI_Status status;
    SAFE_MPI( MPI_Sendrecv( sbuf, sc, MPI_COMPLEX, to,   stag,
                            rbuf, rc, MPI_COMPLEX, from, rtag, comm, 
                            &status ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
Elemental::wrappers::MPI::SendRecv
( dcomplex* sbuf, int sc, int to,   int stag,
  dcomplex* rbuf, int rc, int from, int rtag, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::SendRecv");
#endif
    MPI_Status status;
    SAFE_MPI( MPI_Sendrecv( sbuf, sc, MPI_DOUBLE_COMPLEX, to,   stag,
                            rbuf, rc, MPI_DOUBLE_COMPLEX, from, rtag, comm, 
                            &status ) );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif

inline void
Elemental::wrappers::MPI::Broadcast
( int* buf, int count, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::Broadcast");
#endif
    SAFE_MPI( MPI_Bcast( buf, count, MPI_INT, root, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
Elemental::wrappers::MPI::Broadcast
( float* buf, int count, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::Broadcast");
#endif
    SAFE_MPI( MPI_Bcast( buf, count, MPI_FLOAT, root, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
Elemental::wrappers::MPI::Broadcast
( double* buf, int count, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::Broadcast");
#endif
    SAFE_MPI( MPI_Bcast( buf, count, MPI_DOUBLE, root, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
inline void
Elemental::wrappers::MPI::Broadcast
( scomplex* buf, int count, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::Broadcast");
#endif
    SAFE_MPI( MPI_Bcast( buf, count, MPI_COMPLEX, root, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
Elemental::wrappers::MPI::Broadcast
( dcomplex* buf, int count, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::Broadcast");
#endif
    SAFE_MPI( MPI_Bcast( buf, count, MPI_DOUBLE_COMPLEX, root, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif

inline void
Elemental::wrappers::MPI::Gather
( int* sbuf, int sc,
  int* rbuf, int rc, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::Gather");
#endif
    SAFE_MPI( MPI_Gather( sbuf, sc, MPI_INT,
                          rbuf, rc, MPI_INT, root, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
Elemental::wrappers::MPI::Gather
( float* sbuf, int sc,
  float* rbuf, int rc, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::Gather");
#endif
    SAFE_MPI( MPI_Gather( sbuf, sc, MPI_FLOAT,
                          rbuf, rc, MPI_FLOAT, root, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
Elemental::wrappers::MPI::Gather
( double* sbuf, int sc,
  double* rbuf, int rc, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::Gather");
#endif
    SAFE_MPI( MPI_Gather( sbuf, sc, MPI_DOUBLE,
                          rbuf, rc, MPI_DOUBLE, root, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
inline void
Elemental::wrappers::MPI::Gather
( scomplex* sbuf, int sc,
  scomplex* rbuf, int rc, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::Gather");
#endif
    SAFE_MPI( MPI_Gather( sbuf, sc, MPI_COMPLEX,
                          rbuf, rc, MPI_COMPLEX, root, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
Elemental::wrappers::MPI::Gather
( dcomplex* sbuf, int sc,
  dcomplex* rbuf, int rc, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::Gather");
#endif
    SAFE_MPI( MPI_Gather( sbuf, sc, MPI_DOUBLE_COMPLEX,
                          rbuf, rc, MPI_DOUBLE_COMPLEX, root, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif

inline void
Elemental::wrappers::MPI::AllGather
( int* sbuf, int sc,
  int* rbuf, int rc, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::AllGather");
#endif
    SAFE_MPI( MPI_Allgather( sbuf, sc, MPI_INT,
                             rbuf, rc, MPI_INT, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
Elemental::wrappers::MPI::AllGather
( float* sbuf, int sc,
  float* rbuf, int rc, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::AllGather");
#endif
    SAFE_MPI( MPI_Allgather( sbuf, sc, MPI_FLOAT,
                             rbuf, rc, MPI_FLOAT, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
Elemental::wrappers::MPI::AllGather
( double* sbuf, int sc,
  double* rbuf, int rc, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::AllGather");
#endif
    SAFE_MPI( MPI_Allgather( sbuf, sc, MPI_DOUBLE,
                             rbuf, rc, MPI_DOUBLE, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
inline void
Elemental::wrappers::MPI::AllGather
( scomplex* sbuf, int sc,
  scomplex* rbuf, int rc, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::AllGather");
#endif
    SAFE_MPI( MPI_Allgather( sbuf, sc, MPI_COMPLEX,
                             rbuf, rc, MPI_COMPLEX, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
Elemental::wrappers::MPI::AllGather
( dcomplex* sbuf, int sc,
  dcomplex* rbuf, int rc, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::AllGather");
#endif
    SAFE_MPI( MPI_Allgather( sbuf, sc, MPI_DOUBLE_COMPLEX,
                             rbuf, rc, MPI_DOUBLE_COMPLEX, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif

inline void
Elemental::wrappers::MPI::Scatter
( int* sbuf, int sc,
  int* rbuf, int rc, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::Scatter");
#endif
    SAFE_MPI( MPI_Scatter( sbuf, sc, MPI_INT,
                           rbuf, rc, MPI_INT, root, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
Elemental::wrappers::MPI::Scatter
( float* sbuf, int sc,
  float* rbuf, int rc, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::Scatter");
#endif
    SAFE_MPI( MPI_Scatter( sbuf, sc, MPI_FLOAT,
                           rbuf, rc, MPI_FLOAT, root, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
Elemental::wrappers::MPI::Scatter
( double* sbuf, int sc,
  double* rbuf, int rc, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::Scatter");
#endif
    SAFE_MPI( MPI_Scatter( sbuf, sc, MPI_DOUBLE,
                           rbuf, rc, MPI_DOUBLE, root, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
inline void
Elemental::wrappers::MPI::Scatter
( scomplex* sbuf, int sc,
  scomplex* rbuf, int rc, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::Scatter");
#endif
    SAFE_MPI( MPI_Scatter( sbuf, sc, MPI_COMPLEX,
                           rbuf, rc, MPI_COMPLEX, root, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
Elemental::wrappers::MPI::Scatter
( dcomplex* sbuf, int sc,
  dcomplex* rbuf, int rc, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::Scatter");
#endif
    SAFE_MPI( MPI_Scatter( sbuf, sc, MPI_DOUBLE_COMPLEX,
                           rbuf, rc, MPI_DOUBLE_COMPLEX, root, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif

inline void
Elemental::wrappers::MPI::AllToAll
( int* sbuf, int sc,
  int* rbuf, int rc, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::AllToAll");
#endif
    SAFE_MPI( MPI_Alltoall( sbuf, sc, MPI_INT,
                            rbuf, rc, MPI_INT, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
Elemental::wrappers::MPI::AllToAll
( float* sbuf, int sc,
  float* rbuf, int rc, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::AllToAll");
#endif
    SAFE_MPI( MPI_Alltoall( sbuf, sc, MPI_FLOAT,
                            rbuf, rc, MPI_FLOAT, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
Elemental::wrappers::MPI::AllToAll
( double* sbuf, int sc,
  double* rbuf, int rc, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::AllToAll");
#endif
    SAFE_MPI( MPI_Alltoall( sbuf, sc, MPI_DOUBLE,
                            rbuf, rc, MPI_DOUBLE, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
inline void
Elemental::wrappers::MPI::AllToAll
( scomplex* sbuf, int sc,
  scomplex* rbuf, int rc, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::AllToAll");
#endif
    SAFE_MPI( MPI_Alltoall( sbuf, sc, MPI_COMPLEX,
                            rbuf, rc, MPI_COMPLEX, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
Elemental::wrappers::MPI::AllToAll
( dcomplex* sbuf, int sc,
  dcomplex* rbuf, int rc, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::AllToAll");
#endif
    SAFE_MPI( MPI_Alltoall( sbuf, sc, MPI_DOUBLE_COMPLEX,
                            rbuf, rc, MPI_DOUBLE_COMPLEX, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif

inline void
Elemental::wrappers::MPI::AllToAll
( int* sbuf, int* scs, int* sds, 
  int* rbuf, int* rcs, int* rds, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::AllToAll");
#endif
    SAFE_MPI( MPI_Alltoallv( sbuf, scs, sds, MPI_INT,
                             rbuf, rcs, rds, MPI_INT, comm ) ); 
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
Elemental::wrappers::MPI::AllToAll
( float* sbuf, int* scs, int* sds,
  float* rbuf, int* rcs, int* rds, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::AllToAll");
#endif
    SAFE_MPI( MPI_Alltoallv( sbuf, scs, sds, MPI_FLOAT,
                             rbuf, rcs, rds, MPI_FLOAT, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
Elemental::wrappers::MPI::AllToAll
( double* sbuf, int* scs, int* sds,
  double* rbuf, int* rcs, int* rds, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::AllToAll");
#endif
    SAFE_MPI( MPI_Alltoallv( sbuf, scs, sds, MPI_DOUBLE,
                             rbuf, rcs, rds, MPI_DOUBLE, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
inline void
Elemental::wrappers::MPI::AllToAll
( scomplex* sbuf, int* scs, int* sds,
  scomplex* rbuf, int* rcs, int* rds, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::AllToAll");
#endif
    SAFE_MPI( MPI_Alltoallv( sbuf, scs, sds, MPI_COMPLEX,
                             rbuf, rcs, rds, MPI_COMPLEX, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
Elemental::wrappers::MPI::AllToAll
( dcomplex* sbuf, int* scs, int* sds,
  dcomplex* rbuf, int* rcs, int* rds, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::AllToAll");
#endif
    SAFE_MPI( MPI_Alltoallv( sbuf, scs, sds, MPI_DOUBLE_COMPLEX,
                             rbuf, rcs, rds, MPI_DOUBLE_COMPLEX, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif

inline void
Elemental::wrappers::MPI::Reduce
( int* sbuf, int* rbuf, int count, MPI_Op op, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::Reduce");
#endif
    SAFE_MPI( MPI_Reduce( sbuf, rbuf, count, MPI_INT, op, root, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
Elemental::wrappers::MPI::Reduce
( float* sbuf, float* rbuf, int count, MPI_Op op, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::Reduce");
#endif
    SAFE_MPI( MPI_Reduce( sbuf, rbuf, count, MPI_FLOAT, op, root, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
Elemental::wrappers::MPI::Reduce
( double* sbuf, double* rbuf, int count, MPI_Op op, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::Reduce");
#endif
    SAFE_MPI( MPI_Reduce( sbuf, rbuf, count, MPI_DOUBLE, op, root, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
inline void
Elemental::wrappers::MPI::Reduce
( scomplex* sbuf, scomplex* rbuf, int count, MPI_Op op, int root, 
  MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::Reduce");
#endif
    SAFE_MPI( MPI_Reduce( sbuf, rbuf, count, MPI_COMPLEX, op, root, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
Elemental::wrappers::MPI::Reduce
( dcomplex* sbuf, dcomplex* rbuf, int count, MPI_Op op, int root, 
  MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::Reduce");
#endif
    SAFE_MPI( MPI_Reduce( sbuf, rbuf, count, MPI_DOUBLE_COMPLEX, op, root, 
              comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif

inline void
Elemental::wrappers::MPI::AllReduce
( char* sbuf, char* rbuf, int count, MPI_Op op, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::AllReduce");
#endif
    SAFE_MPI( MPI_Allreduce( sbuf, rbuf, count, MPI_CHAR, op, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
Elemental::wrappers::MPI::AllReduce
( int* sbuf, int* rbuf, int count, MPI_Op op, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::AllReduce");
#endif
    SAFE_MPI( MPI_Allreduce( sbuf, rbuf, count, MPI_INT, op, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
Elemental::wrappers::MPI::AllReduce
( float* sbuf, float* rbuf, int count, MPI_Op op, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::AllReduce");
#endif
    SAFE_MPI( MPI_Allreduce( sbuf, rbuf, count, MPI_FLOAT, op, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
Elemental::wrappers::MPI::AllReduce
( double* sbuf, double* rbuf, int count, MPI_Op op, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::AllReduce");
#endif
    SAFE_MPI( MPI_Allreduce( sbuf, rbuf, count, MPI_DOUBLE, op, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
inline void
Elemental::wrappers::MPI::AllReduce
( scomplex* sbuf, scomplex* rbuf, int count, 
  MPI_Op op, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::AllReduce");
#endif
    SAFE_MPI( MPI_Allreduce( sbuf, rbuf, count, MPI_COMPLEX, op, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
Elemental::wrappers::MPI::AllReduce
( dcomplex* sbuf, dcomplex* rbuf, int count, MPI_Op op, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::AllReduce");
#endif
    SAFE_MPI( MPI_Allreduce( sbuf, rbuf, count, MPI_DOUBLE_COMPLEX, op, 
                             comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif

inline void
Elemental::wrappers::MPI::ReduceScatter
( int* sbuf, int* rbuf, int* rcs, MPI_Op op, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::ReduceScatter");
#endif
    SAFE_MPI( MPI_Reduce_scatter( sbuf, rbuf, rcs, MPI_INT, op, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
Elemental::wrappers::MPI::ReduceScatter
( float* sbuf, float* rbuf, int* rcs, MPI_Op op, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::ReduceScatter");
#endif
    SAFE_MPI( MPI_Reduce_scatter( sbuf, rbuf, rcs, MPI_FLOAT, op, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
Elemental::wrappers::MPI::ReduceScatter
( double* sbuf, double* rbuf, int* rcs, MPI_Op op, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::ReduceScatter");
#endif
    SAFE_MPI( MPI_Reduce_scatter( sbuf, rbuf, rcs, MPI_DOUBLE, op, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
inline void
Elemental::wrappers::MPI::ReduceScatter
( scomplex* sbuf, scomplex* rbuf, int* rcs, MPI_Op op, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::ReduceScatter");
#endif
    SAFE_MPI( MPI_Reduce_scatter( sbuf, rbuf, rcs, MPI_COMPLEX, op, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
Elemental::wrappers::MPI::ReduceScatter
( dcomplex* sbuf, dcomplex* rbuf, int* rcs, MPI_Op op, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::MPI::ReduceScatter");
#endif
    SAFE_MPI( MPI_Reduce_scatter( sbuf, rbuf, rcs, MPI_DOUBLE_COMPLEX, op, 
                                  comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif

#undef SAFE_MPI

#endif /* ELEMENTAL_WRAPPERS_MPI_HPP */

