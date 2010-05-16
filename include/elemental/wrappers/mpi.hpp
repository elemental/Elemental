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

#include "elemental/environment.hpp"

namespace elemental {
namespace wrappers {
namespace mpi {

// Minimum number of elements that each process must contribute
// in collective communications
const int MinCollectContrib = 1;

double 
Time();

void 
Barrier( MPI_Comm comm );

void
Send
( const int* buf, int count, int to, int tag, MPI_Comm comm );

void
Send
( const float* buf, int count, int to, int tag, MPI_Comm comm );

void
Send
( const double* buf, int count, int to, int tag, MPI_Comm comm );

#ifndef WITHOUT_COMPLEX
void
Send
( const scomplex* buf, int count, int to, int tag, MPI_Comm comm );

void
Send
( const dcomplex* buf, int count, int to, int tag, MPI_Comm comm );
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
( const int* sbuf, int sc, int to,   int stag,
        int* rbuf, int rc, int from, int rtag, MPI_Comm comm );

void
SendRecv
( const float* sbuf, int sc, int to,   int stag,
        float* rbuf, int rc, int from, int rtag, MPI_Comm comm );

void
SendRecv
( const double* sbuf, int sc, int to,   int stag,
        double* rbuf, int rc, int from, int rtag, MPI_Comm comm );

#ifndef WITHOUT_COMPLEX
void
SendRecv
( const scomplex* sbuf, int sc, int to,   int stag,
        scomplex* rbuf, int rc, int from, int rtag, MPI_Comm comm );

void
SendRecv
( const dcomplex* sbuf, int sc, int to,   int stag,
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
( const int* sbuf, int sc,
        int* rbuf, int rc, int root, MPI_Comm comm );

void
Gather
( const float* sbuf, int sc,
        float* rbuf, int rc, int root, MPI_Comm comm );

void
Gather
( const double* sbuf, int sc,
        double* rbuf, int rc, int root, MPI_Comm comm );

#ifndef WITHOUT_COMPLEX
void 
Gather
( const scomplex* sbuf, int sc,
        scomplex* rbuf, int rc, int root, MPI_Comm comm );

void
Gather
( const dcomplex* sbuf, int sc,
        dcomplex* rbuf, int rc, int root, MPI_Comm comm );
#endif
    
void
AllGather
( const int* sbuf, int sc,
        int* rbuf, int rc, MPI_Comm comm );

void
AllGather
( const float* sbuf, int sc,
        float* rbuf, int rc, MPI_Comm comm );

void
AllGather
( const double* sbuf, int sc,
        double* rbuf, int rc, MPI_Comm comm );

#ifndef WITHOUT_COMPLEX
void
AllGather
( const scomplex* sbuf, int sc,
        scomplex* rbuf, int rc, MPI_Comm comm );

void
AllGather
( const dcomplex* sbuf, int sc,
        dcomplex* rbuf, int rc, MPI_Comm comm );
#endif

void
Scatter
( const int* sbuf, int sc,
        int* rbuf, int rc, int root, MPI_Comm comm );

void
Scatter
( const float* sbuf, int sc,
        float* rbuf, int rc, int root, MPI_Comm comm );

void
Scatter
( const double* sbuf, int sc,
        double* rbuf, int rc, int root, MPI_Comm comm );

#ifndef WITHOUT_COMPLEX
void
Scatter
( const scomplex* sbuf, int sc,
        scomplex* rbuf, int rc, int root, MPI_Comm comm );

void
Scatter
( const dcomplex* sbuf, int sc,
        dcomplex* rbuf, int rc, int root, MPI_Comm comm );
#endif
    
void
AllToAll
( const int* sbuf, int sc,
        int* rbuf, int rc, MPI_Comm comm );

void
AllToAll
( const float* sbuf, int sc,
        float* rbuf, int rc, MPI_Comm comm );

void
AllToAll
( const double* sbuf, int sc,
        double* rbuf, int rc, MPI_Comm comm );

#ifndef WITHOUT_COMPLEX
void
AllToAll
( const scomplex* sbuf, int sc,
        scomplex* rbuf, int rc, MPI_Comm comm );

void
AllToAll
( const dcomplex* sbuf, int sc,
        dcomplex* rbuf, int rc, MPI_Comm comm );
#endif

void
AllToAll
( const int* sbuf, const int* scs, const int* sds,
        int* rbuf, const int* rcs, const int* rds, MPI_Comm comm );

void
AllToAll
( const float* sbuf, const int* scs, const int* sds,
        float* rbuf, const int* rcs, const int* rds, MPI_Comm comm );

void
AllToAll
( const double* sbuf, const int* scs, const int* sds,
        double* rbuf, const int* rcs, const int* rds, MPI_Comm comm );

#ifndef WITHOUT_COMPLEX
void
AllToAll
( const scomplex* sbuf, const int* scs, const int* sds,
        scomplex* rbuf, const int* rcs, const int* rds, MPI_Comm comm );

void
AllToAll
( const dcomplex* sbuf, const int* scs, const int* sds,
        dcomplex* rbuf, const int* rcs, const int* rds, MPI_Comm comm );
#endif

void
Reduce
( const int* sbuf, int* rbuf, int count, MPI_Op op, int root, MPI_Comm comm );

void
Reduce
( const float* sbuf, float* rbuf, int count, MPI_Op op, int root,
  MPI_Comm comm );

void
Reduce
( const double* sbuf, double* rbuf, int count, MPI_Op op, int root,
  MPI_Comm comm );

#ifndef WITHOUT_COMPLEX
void
Reduce
( const scomplex* sbuf, scomplex* rbuf, int count, MPI_Op op, int root, 
  MPI_Comm comm );

void
Reduce
( const dcomplex* sbuf, dcomplex* rbuf, int count, MPI_Op op, int root, 
  MPI_Comm comm );
#endif
    
void
AllReduce
( const char* sbuf, char* rbuf, int count, MPI_Op op, MPI_Comm comm );

void
AllReduce
( const int* sbuf, int* rbuf, int count, MPI_Op op, MPI_Comm comm );

void
AllReduce
( const float* sbuf, float* rbuf, int count, MPI_Op op, MPI_Comm comm );

void
AllReduce
( const double* sbuf, double* rbuf, int count, MPI_Op op, MPI_Comm comm );
#ifndef WITHOUT_COMPLEX
void
AllReduce
( const scomplex* sbuf, scomplex* rbuf, int count, MPI_Op op, MPI_Comm comm );

void
AllReduce
( const dcomplex* sbuf, dcomplex* rbuf, int count, MPI_Op op, MPI_Comm comm );
#endif
        
void
ReduceScatter
( const int* sbuf, int* rbuf, const int* rcs, MPI_Op op, MPI_Comm comm );

void
ReduceScatter
( const float* sbuf, float* rbuf, const int* rcs, MPI_Op op, MPI_Comm comm );

void
ReduceScatter
( const double* sbuf, double* rbuf, const int* rcs, MPI_Op op, MPI_Comm comm );
#ifndef WITHOUT_COMPLEX
void
ReduceScatter
( const scomplex* sbuf, scomplex* rbuf, const int* rcs, MPI_Op op, 
  MPI_Comm com );

void
ReduceScatter
( const dcomplex* sbuf, dcomplex* rbuf, const int* rcs, MPI_Op op, 
  MPI_Comm comm );
#endif

} // mpi
} // wrapppers
} // elemental

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
        throw errorString; \
    } \
}
#else
#define SAFE_MPI( routine ) routine
#endif

inline double
elemental::wrappers::mpi::Time()
{ return MPI_Wtime(); }

inline void
elemental::wrappers::mpi::Barrier( MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::Barrier");
#endif
    SAFE_MPI( MPI_Barrier( comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
elemental::wrappers::mpi::Send
( const int* buf, int count, int to, int tag, MPI_Comm comm )
{ 
#ifndef RELEASE
    PushCallStack("wrappers::mpi::Send");
#endif
    SAFE_MPI( 
        MPI_Send( const_cast<int*>(buf), count, MPI_INT, to, tag, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
elemental::wrappers::mpi::Send
( const float* buf, int count, int to, int tag, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::Send");
#endif
    SAFE_MPI( 
        MPI_Send( const_cast<float*>(buf), count, MPI_FLOAT, to, tag, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
elemental::wrappers::mpi::Send
( const double* buf, int count, int to, int tag, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::Send");
#endif
    SAFE_MPI( 
        MPI_Send( const_cast<double*>(buf), count, MPI_DOUBLE, to, tag, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
inline void
elemental::wrappers::mpi::Send
( const scomplex* buf, int count, int to, int tag, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::Send");
#endif
    SAFE_MPI( 
        MPI_Send( const_cast<scomplex*>(buf), count, MPI_COMPLEX, to, tag, 
                  comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
elemental::wrappers::mpi::Send
( const dcomplex* buf, int count, int to, int tag, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::Send");
#endif
    SAFE_MPI( 
        MPI_Send( const_cast<dcomplex*>(buf), count, MPI_DOUBLE_COMPLEX, to, 
                  tag, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif

inline void
elemental::wrappers::mpi::Recv
( int* buf, int count, int from, int tag, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::Recv");
#endif
    MPI_Status status;
    SAFE_MPI( MPI_Recv( buf, count, MPI_INT, from, tag, comm, &status ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
elemental::wrappers::mpi::Recv
( float* buf, int count, int from, int tag, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::Recv");
#endif
    MPI_Status status;
    SAFE_MPI( MPI_Recv( buf, count, MPI_FLOAT, from, tag, comm, &status ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
elemental::wrappers::mpi::Recv
( double* buf, int count, int from, int tag, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::Recv");
#endif
    MPI_Status status;
    SAFE_MPI( MPI_Recv( buf, count, MPI_DOUBLE, from, tag, comm, &status ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
inline void
elemental::wrappers::mpi::Recv
( scomplex* buf, int count, int from, int tag, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::Recv");
#endif
    MPI_Status status;
    SAFE_MPI( MPI_Recv( buf, count, MPI_COMPLEX, from, tag, comm, &status ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
elemental::wrappers::mpi::Recv
( dcomplex* buf, int count, int from, int tag, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::Recv");
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
elemental::wrappers::mpi::SendRecv
( const int* sbuf, int sc, int to,   int stag,
        int* rbuf, int rc, int from, int rtag, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::SendRecv");
#endif
    MPI_Status status;
    SAFE_MPI( 
        MPI_Sendrecv( const_cast<int*>(sbuf), sc, MPI_INT, to, stag,
                      rbuf, rc, MPI_INT, from, rtag, comm, &status ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
elemental::wrappers::mpi::SendRecv
( const float* sbuf, int sc, int to,   int stag,
        float* rbuf, int rc, int from, int rtag, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::SendRecv");
#endif
    MPI_Status status;
    SAFE_MPI( 
        MPI_Sendrecv( const_cast<float*>(sbuf), sc, MPI_FLOAT, to, stag,
                      rbuf, rc, MPI_FLOAT, from, rtag, comm, &status ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
elemental::wrappers::mpi::SendRecv
( const double* sbuf, int sc, int to,   int stag,
        double* rbuf, int rc, int from, int rtag, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::SendRecv");
#endif
    MPI_Status status;
    SAFE_MPI( 
        MPI_Sendrecv( const_cast<double*>(sbuf), sc, MPI_DOUBLE, to, stag,
                      rbuf, rc, MPI_DOUBLE, from, rtag, comm, &status ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
inline void
elemental::wrappers::mpi::SendRecv
( const scomplex* sbuf, int sc, int to,   int stag,
        scomplex* rbuf, int rc, int from, int rtag, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::SendRecv");
#endif
    MPI_Status status;
    SAFE_MPI( 
        MPI_Sendrecv( const_cast<scomplex*>(sbuf), sc, MPI_COMPLEX, to, stag,
                      rbuf, rc, MPI_COMPLEX, from, rtag, comm, &status )
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
elemental::wrappers::mpi::SendRecv
( const dcomplex* sbuf, int sc, int to,   int stag,
        dcomplex* rbuf, int rc, int from, int rtag, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::SendRecv");
#endif
    MPI_Status status;
    SAFE_MPI( 
        MPI_Sendrecv( const_cast<dcomplex*>(sbuf), 
                            sc, MPI_DOUBLE_COMPLEX, to,   stag,
                      rbuf, rc, MPI_DOUBLE_COMPLEX, from, rtag, comm, &status )
    );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif

inline void
elemental::wrappers::mpi::Broadcast
( int* buf, int count, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::Broadcast");
#endif
    SAFE_MPI( MPI_Bcast( buf, count, MPI_INT, root, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
elemental::wrappers::mpi::Broadcast
( float* buf, int count, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::Broadcast");
#endif
    SAFE_MPI( MPI_Bcast( buf, count, MPI_FLOAT, root, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
elemental::wrappers::mpi::Broadcast
( double* buf, int count, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::Broadcast");
#endif
    SAFE_MPI( MPI_Bcast( buf, count, MPI_DOUBLE, root, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
inline void
elemental::wrappers::mpi::Broadcast
( scomplex* buf, int count, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::Broadcast");
#endif
    SAFE_MPI( MPI_Bcast( buf, count, MPI_COMPLEX, root, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
elemental::wrappers::mpi::Broadcast
( dcomplex* buf, int count, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::Broadcast");
#endif
    SAFE_MPI( MPI_Bcast( buf, count, MPI_DOUBLE_COMPLEX, root, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif

inline void
elemental::wrappers::mpi::Gather
( const int* sbuf, int sc,
        int* rbuf, int rc, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::Gather");
#endif
    SAFE_MPI( MPI_Gather( const_cast<int*>(sbuf), sc, MPI_INT,
                          rbuf, rc, MPI_INT, root, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
elemental::wrappers::mpi::Gather
( const float* sbuf, int sc,
        float* rbuf, int rc, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::Gather");
#endif
    SAFE_MPI( MPI_Gather( const_cast<float*>(sbuf), sc, MPI_FLOAT,
                          rbuf, rc, MPI_FLOAT, root, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
elemental::wrappers::mpi::Gather
( const double* sbuf, int sc,
        double* rbuf, int rc, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::Gather");
#endif
    SAFE_MPI( MPI_Gather( const_cast<double*>(sbuf), sc, MPI_DOUBLE,
                          rbuf, rc, MPI_DOUBLE, root, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
inline void
elemental::wrappers::mpi::Gather
( const scomplex* sbuf, int sc,
        scomplex* rbuf, int rc, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::Gather");
#endif
    SAFE_MPI( MPI_Gather( const_cast<scomplex*>(sbuf), sc, MPI_COMPLEX,
                          rbuf, rc, MPI_COMPLEX, root, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
elemental::wrappers::mpi::Gather
( const dcomplex* sbuf, int sc,
        dcomplex* rbuf, int rc, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::Gather");
#endif
    SAFE_MPI( MPI_Gather( const_cast<dcomplex*>(sbuf), sc, MPI_DOUBLE_COMPLEX,
                          rbuf, rc, MPI_DOUBLE_COMPLEX, root, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif

inline void
elemental::wrappers::mpi::AllGather
( const int* sbuf, int sc,
        int* rbuf, int rc, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::AllGather");
#endif
    SAFE_MPI( MPI_Allgather( const_cast<int*>(sbuf), sc, MPI_INT,
                             rbuf, rc, MPI_INT, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
elemental::wrappers::mpi::AllGather
( const float* sbuf, int sc,
        float* rbuf, int rc, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::AllGather");
#endif
    SAFE_MPI( MPI_Allgather( const_cast<float*>(sbuf), sc, MPI_FLOAT,
                             rbuf, rc, MPI_FLOAT, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
elemental::wrappers::mpi::AllGather
( const double* sbuf, int sc,
        double* rbuf, int rc, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::AllGather");
#endif
    SAFE_MPI( MPI_Allgather( const_cast<double*>(sbuf), sc, MPI_DOUBLE,
                             rbuf, rc, MPI_DOUBLE, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
inline void
elemental::wrappers::mpi::AllGather
( const scomplex* sbuf, int sc,
        scomplex* rbuf, int rc, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::AllGather");
#endif
    SAFE_MPI( MPI_Allgather( const_cast<scomplex*>(sbuf), sc, MPI_COMPLEX,
                             rbuf, rc, MPI_COMPLEX, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
elemental::wrappers::mpi::AllGather
( const dcomplex* sbuf, int sc,
        dcomplex* rbuf, int rc, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::AllGather");
#endif
    SAFE_MPI( MPI_Allgather( const_cast<dcomplex*>(sbuf), 
                                   sc, MPI_DOUBLE_COMPLEX,
                             rbuf, rc, MPI_DOUBLE_COMPLEX, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif

inline void
elemental::wrappers::mpi::Scatter
( const int* sbuf, int sc,
        int* rbuf, int rc, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::Scatter");
#endif
    SAFE_MPI( MPI_Scatter( const_cast<int*>(sbuf), sc, MPI_INT,
                           rbuf, rc, MPI_INT, root, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
elemental::wrappers::mpi::Scatter
( const float* sbuf, int sc,
        float* rbuf, int rc, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::Scatter");
#endif
    SAFE_MPI( MPI_Scatter( const_cast<float*>(sbuf), sc, MPI_FLOAT,
                           rbuf, rc, MPI_FLOAT, root, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
elemental::wrappers::mpi::Scatter
( const double* sbuf, int sc,
        double* rbuf, int rc, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::Scatter");
#endif
    SAFE_MPI( MPI_Scatter( const_cast<double*>(sbuf), sc, MPI_DOUBLE,
                           rbuf, rc, MPI_DOUBLE, root, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
inline void
elemental::wrappers::mpi::Scatter
( const scomplex* sbuf, int sc,
        scomplex* rbuf, int rc, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::Scatter");
#endif
    SAFE_MPI( MPI_Scatter( const_cast<scomplex*>(sbuf), sc, MPI_COMPLEX,
                           rbuf, rc, MPI_COMPLEX, root, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
elemental::wrappers::mpi::Scatter
( const dcomplex* sbuf, int sc,
        dcomplex* rbuf, int rc, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::Scatter");
#endif
    SAFE_MPI( MPI_Scatter( const_cast<dcomplex*>(sbuf), sc, MPI_DOUBLE_COMPLEX,
                           rbuf, rc, MPI_DOUBLE_COMPLEX, root, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif

inline void
elemental::wrappers::mpi::AllToAll
( const int* sbuf, int sc,
        int* rbuf, int rc, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::AllToAll");
#endif
    SAFE_MPI( MPI_Alltoall( const_cast<int*>(sbuf), sc, MPI_INT,
                            rbuf, rc, MPI_INT, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
elemental::wrappers::mpi::AllToAll
( const float* sbuf, int sc,
        float* rbuf, int rc, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::AllToAll");
#endif
    SAFE_MPI( MPI_Alltoall( const_cast<float*>(sbuf), sc, MPI_FLOAT,
                            rbuf, rc, MPI_FLOAT, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
elemental::wrappers::mpi::AllToAll
( const double* sbuf, int sc,
        double* rbuf, int rc, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::AllToAll");
#endif
    SAFE_MPI( MPI_Alltoall( const_cast<double*>(sbuf), sc, MPI_DOUBLE,
                            rbuf, rc, MPI_DOUBLE, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
inline void
elemental::wrappers::mpi::AllToAll
( const scomplex* sbuf, int sc,
        scomplex* rbuf, int rc, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::AllToAll");
#endif
    SAFE_MPI( MPI_Alltoall( const_cast<scomplex*>(sbuf), sc, MPI_COMPLEX,
                            rbuf, rc, MPI_COMPLEX, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
elemental::wrappers::mpi::AllToAll
( const dcomplex* sbuf, int sc,
        dcomplex* rbuf, int rc, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::AllToAll");
#endif
    SAFE_MPI( MPI_Alltoall( const_cast<dcomplex*>(sbuf), sc, MPI_DOUBLE_COMPLEX,
                            rbuf, rc, MPI_DOUBLE_COMPLEX, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif

inline void
elemental::wrappers::mpi::AllToAll
( const int* sbuf, const int* scs, const int* sds, 
        int* rbuf, const int* rcs, const int* rds, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::AllToAll");
#endif
    SAFE_MPI( 
        MPI_Alltoallv( const_cast<int*>(sbuf), 
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

inline void
elemental::wrappers::mpi::AllToAll
( const float* sbuf, const int* scs, const int* sds,
        float* rbuf, const int* rcs, const int* rds, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::AllToAll");
#endif
    SAFE_MPI( 
        MPI_Alltoallv( const_cast<float*>(sbuf), 
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

inline void
elemental::wrappers::mpi::AllToAll
( const double* sbuf, const int* scs, const int* sds,
        double* rbuf, const int* rcs, const int* rds, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::AllToAll");
#endif
    SAFE_MPI( 
        MPI_Alltoallv( const_cast<double*>(sbuf), 
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

#ifndef WITHOUT_COMPLEX
inline void
elemental::wrappers::mpi::AllToAll
( const scomplex* sbuf, const int* scs, const int* sds,
        scomplex* rbuf, const int* rcs, const int* rds, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::AllToAll");
#endif
    SAFE_MPI( 
        MPI_Alltoallv( const_cast<scomplex*>(sbuf), 
                       const_cast<int*>(scs), 
                       const_cast<int*>(sds), 
                       MPI_COMPLEX,
                       rbuf, 
                       const_cast<int*>(rcs), 
                       const_cast<int*>(rds), 
                       MPI_COMPLEX, 
                       comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
elemental::wrappers::mpi::AllToAll
( const dcomplex* sbuf, const int* scs, const int* sds,
        dcomplex* rbuf, const int* rcs, const int* rds, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::AllToAll");
#endif
    SAFE_MPI( 
        MPI_Alltoallv( const_cast<dcomplex*>(sbuf), 
                       const_cast<int*>(scs), 
                       const_cast<int*>(sds), 
                       MPI_DOUBLE_COMPLEX,
                       rbuf, 
                       const_cast<int*>(rcs), 
                       const_cast<int*>(rds), 
                       MPI_DOUBLE_COMPLEX, 
                       comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif

inline void
elemental::wrappers::mpi::Reduce
( const int* sbuf, int* rbuf, int count, MPI_Op op, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::Reduce");
#endif
    SAFE_MPI( 
        MPI_Reduce( const_cast<int*>(sbuf), 
                    rbuf, count, MPI_INT, op, root, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
elemental::wrappers::mpi::Reduce
( const float* sbuf, float* rbuf, int count, MPI_Op op, int root, 
  MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::Reduce");
#endif
    SAFE_MPI( 
        MPI_Reduce( const_cast<float*>(sbuf), 
                    rbuf, count, MPI_FLOAT, op, root, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
elemental::wrappers::mpi::Reduce
( const double* sbuf, double* rbuf, int count, MPI_Op op, int root, 
  MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::Reduce");
#endif
    SAFE_MPI( 
        MPI_Reduce( const_cast<double*>(sbuf), 
                    rbuf, count, MPI_DOUBLE, op, root, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
inline void
elemental::wrappers::mpi::Reduce
( const scomplex* sbuf, scomplex* rbuf, int count, MPI_Op op, int root, 
  MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::Reduce");
#endif
    SAFE_MPI( 
        MPI_Reduce( const_cast<scomplex*>(sbuf), 
                    rbuf, count, MPI_COMPLEX, op, root, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
elemental::wrappers::mpi::Reduce
( const dcomplex* sbuf, dcomplex* rbuf, int count, MPI_Op op, int root, 
  MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::Reduce");
#endif
    SAFE_MPI( 
        MPI_Reduce( const_cast<dcomplex*>(sbuf), 
                    rbuf, count, MPI_DOUBLE_COMPLEX, op, root, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif

inline void
elemental::wrappers::mpi::AllReduce
( const char* sbuf, char* rbuf, int count, MPI_Op op, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::AllReduce");
#endif
    SAFE_MPI( 
        MPI_Allreduce( const_cast<char*>(sbuf), 
                       rbuf, count, MPI_CHAR, op, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
elemental::wrappers::mpi::AllReduce
( const int* sbuf, int* rbuf, int count, MPI_Op op, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::AllReduce");
#endif
    SAFE_MPI( 
        MPI_Allreduce( const_cast<int*>(sbuf), 
                       rbuf, count, MPI_INT, op, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
elemental::wrappers::mpi::AllReduce
( const float* sbuf, float* rbuf, int count, MPI_Op op, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::AllReduce");
#endif
    SAFE_MPI( 
        MPI_Allreduce( const_cast<float*>(sbuf), 
                       rbuf, count, MPI_FLOAT, op, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
elemental::wrappers::mpi::AllReduce
( const double* sbuf, double* rbuf, int count, MPI_Op op, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::AllReduce");
#endif
    SAFE_MPI( 
        MPI_Allreduce( const_cast<double*>(sbuf), 
                       rbuf, count, MPI_DOUBLE, op, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
inline void
elemental::wrappers::mpi::AllReduce
( const scomplex* sbuf, scomplex* rbuf, int count, MPI_Op op, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::AllReduce");
#endif
    SAFE_MPI( 
        MPI_Allreduce( const_cast<scomplex*>(sbuf), 
                       rbuf, count, MPI_COMPLEX, op, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
elemental::wrappers::mpi::AllReduce
( const dcomplex* sbuf, dcomplex* rbuf, int count, MPI_Op op, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::AllReduce");
#endif
    SAFE_MPI( 
        MPI_Allreduce( const_cast<dcomplex*>(sbuf), 
                       rbuf, count, MPI_DOUBLE_COMPLEX, op, comm )
    );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif

inline void
elemental::wrappers::mpi::ReduceScatter
( const int* sbuf, int* rbuf, const int* rcs, MPI_Op op, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::ReduceScatter");
#endif
    SAFE_MPI( 
        MPI_Reduce_scatter( const_cast<int*>(sbuf), 
                            rbuf, 
                            const_cast<int*>(rcs), MPI_INT, op, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
elemental::wrappers::mpi::ReduceScatter
( const float* sbuf, float* rbuf, const int* rcs, MPI_Op op, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::ReduceScatter");
#endif
    SAFE_MPI( 
        MPI_Reduce_scatter( const_cast<float*>(sbuf), 
                            rbuf, 
                            const_cast<int*>(rcs), MPI_FLOAT, op, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
elemental::wrappers::mpi::ReduceScatter
( const double* sbuf, double* rbuf, const int* rcs, MPI_Op op, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::ReduceScatter");
#endif
    SAFE_MPI( 
        MPI_Reduce_scatter( const_cast<double*>(sbuf), 
                            rbuf, 
                            const_cast<int*>(rcs), MPI_DOUBLE, op, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
inline void
elemental::wrappers::mpi::ReduceScatter
( const scomplex* sbuf, scomplex* rbuf, const int* rcs, MPI_Op op, 
  MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::ReduceScatter");
#endif
    SAFE_MPI( 
        MPI_Reduce_scatter( const_cast<scomplex*>(sbuf), 
                            rbuf, 
                            const_cast<int*>(rcs), MPI_COMPLEX, op, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
elemental::wrappers::mpi::ReduceScatter
( const dcomplex* sbuf, dcomplex* rbuf, const int* rcs, MPI_Op op, 
  MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::ReduceScatter");
#endif
    SAFE_MPI( 
        MPI_Reduce_scatter( const_cast<dcomplex*>(sbuf), 
                            rbuf, 
                            const_cast<int*>(rcs), 
                            MPI_DOUBLE_COMPLEX, op, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif

#undef SAFE_MPI

#endif /* ELEMENTAL_WRAPPERS_MPI_HPP */

