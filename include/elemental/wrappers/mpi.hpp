/*
   This file is part of elemental, a library for distributed-memory dense 
   linear algebra.

   Copyright (C) 2009-2010 Jack Poulson <jack.poulson@gmail.com>

   This program is released under the terms of the license contained in the 
   file LICENSE.
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
        throw errorString;
    }
#endif
}

} // anonymous namespace

inline double
elemental::wrappers::mpi::Time()
{ return MPI_Wtime(); }

inline void
elemental::wrappers::mpi::Barrier( MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::Barrier");
#endif
    SafeMpi( MPI_Barrier( comm ) );
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
    SafeMpi( 
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
    SafeMpi( 
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
    SafeMpi( 
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
    SafeMpi( 
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
    SafeMpi( 
        MPI_Send( const_cast<dcomplex*>(buf), count, MPI_DOUBLE_COMPLEX, to, 
                  tag, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

inline void
elemental::wrappers::mpi::Recv
( int* buf, int count, int from, int tag, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::Recv");
#endif
    MPI_Status status;
    SafeMpi( MPI_Recv( buf, count, MPI_INT, from, tag, comm, &status ) );
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
    SafeMpi( MPI_Recv( buf, count, MPI_FLOAT, from, tag, comm, &status ) );
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
    SafeMpi( MPI_Recv( buf, count, MPI_DOUBLE, from, tag, comm, &status ) );
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
    SafeMpi( MPI_Recv( buf, count, MPI_COMPLEX, from, tag, comm, &status ) );
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
    SafeMpi( 
        MPI_Recv( buf, count, MPI_DOUBLE_COMPLEX, from, tag, comm, &status );
    );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

inline void
elemental::wrappers::mpi::SendRecv
( const int* sbuf, int sc, int to,   int stag,
        int* rbuf, int rc, int from, int rtag, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::SendRecv");
#endif
    MPI_Status status;
    SafeMpi( 
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
    SafeMpi( 
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
    SafeMpi( 
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
    SafeMpi( 
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
    SafeMpi( 
        MPI_Sendrecv( const_cast<dcomplex*>(sbuf), 
                            sc, MPI_DOUBLE_COMPLEX, to,   stag,
                      rbuf, rc, MPI_DOUBLE_COMPLEX, from, rtag, comm, &status )
    );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

inline void
elemental::wrappers::mpi::Broadcast
( int* buf, int count, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::Broadcast");
#endif
    SafeMpi( MPI_Bcast( buf, count, MPI_INT, root, comm ) );
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
    SafeMpi( MPI_Bcast( buf, count, MPI_FLOAT, root, comm ) );
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
    SafeMpi( MPI_Bcast( buf, count, MPI_DOUBLE, root, comm ) );
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
    SafeMpi( MPI_Bcast( buf, count, MPI_COMPLEX, root, comm ) );
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
    SafeMpi( MPI_Bcast( buf, count, MPI_DOUBLE_COMPLEX, root, comm ) );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

inline void
elemental::wrappers::mpi::Gather
( const int* sbuf, int sc,
        int* rbuf, int rc, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::Gather");
#endif
    SafeMpi( 
        MPI_Gather( const_cast<int*>(sbuf), sc, MPI_INT,
                    rbuf, rc, MPI_INT, root, comm ) 
    );
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
    SafeMpi( 
        MPI_Gather( const_cast<float*>(sbuf), sc, MPI_FLOAT,
                    rbuf, rc, MPI_FLOAT, root, comm ) 
    );
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
    SafeMpi( 
        MPI_Gather( const_cast<double*>(sbuf), sc, MPI_DOUBLE,
                    rbuf, rc, MPI_DOUBLE, root, comm ) 
    );
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
    SafeMpi( 
        MPI_Gather( const_cast<scomplex*>(sbuf), sc, MPI_COMPLEX,
                    rbuf, rc, MPI_COMPLEX, root, comm ) 
    );
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
    SafeMpi( 
        MPI_Gather( const_cast<dcomplex*>(sbuf), sc, MPI_DOUBLE_COMPLEX,
                    rbuf, rc, MPI_DOUBLE_COMPLEX, root, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

inline void
elemental::wrappers::mpi::AllGather
( const int* sbuf, int sc,
        int* rbuf, int rc, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::AllGather");
#endif
    SafeMpi( 
        MPI_Allgather( const_cast<int*>(sbuf), sc, MPI_INT,
                       rbuf, rc, MPI_INT, comm ) 
    );
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
    SafeMpi( 
        MPI_Allgather( const_cast<float*>(sbuf), sc, MPI_FLOAT,
                       rbuf, rc, MPI_FLOAT, comm ) 
    );
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
    SafeMpi( 
        MPI_Allgather( const_cast<double*>(sbuf), sc, MPI_DOUBLE,
                       rbuf, rc, MPI_DOUBLE, comm ) 
    );
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
    SafeMpi( 
        MPI_Allgather( const_cast<scomplex*>(sbuf), sc, MPI_COMPLEX,
                       rbuf, rc, MPI_COMPLEX, comm ) 
    );
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
    SafeMpi( 
        MPI_Allgather( const_cast<dcomplex*>(sbuf), 
                             sc, MPI_DOUBLE_COMPLEX,
                       rbuf, rc, MPI_DOUBLE_COMPLEX, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

inline void
elemental::wrappers::mpi::Scatter
( const int* sbuf, int sc,
        int* rbuf, int rc, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::Scatter");
#endif
    SafeMpi( 
        MPI_Scatter( const_cast<int*>(sbuf), sc, MPI_INT,
                     rbuf, rc, MPI_INT, root, comm ) 
    );
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
    SafeMpi( 
        MPI_Scatter( const_cast<float*>(sbuf), sc, MPI_FLOAT,
                     rbuf, rc, MPI_FLOAT, root, comm ) 
    );
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
    SafeMpi( 
        MPI_Scatter( const_cast<double*>(sbuf), sc, MPI_DOUBLE,
                     rbuf, rc, MPI_DOUBLE, root, comm ) 
    );
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
    SafeMpi( 
        MPI_Scatter( const_cast<scomplex*>(sbuf), sc, MPI_COMPLEX,
                     rbuf, rc, MPI_COMPLEX, root, comm ) 
    );
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
    SafeMpi( 
        MPI_Scatter( const_cast<dcomplex*>(sbuf), sc, MPI_DOUBLE_COMPLEX,
                     rbuf, rc, MPI_DOUBLE_COMPLEX, root, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

inline void
elemental::wrappers::mpi::AllToAll
( const int* sbuf, int sc,
        int* rbuf, int rc, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::AllToAll");
#endif
    SafeMpi( 
        MPI_Alltoall( const_cast<int*>(sbuf), sc, MPI_INT,
                      rbuf, rc, MPI_INT, comm ) 
    );
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
    SafeMpi( 
        MPI_Alltoall( const_cast<float*>(sbuf), sc, MPI_FLOAT,
                      rbuf, rc, MPI_FLOAT, comm ) 
    );
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
    SafeMpi( 
        MPI_Alltoall( const_cast<double*>(sbuf), sc, MPI_DOUBLE,
                      rbuf, rc, MPI_DOUBLE, comm ) 
    );
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
    SafeMpi( 
        MPI_Alltoall( const_cast<scomplex*>(sbuf), sc, MPI_COMPLEX,
                      rbuf, rc, MPI_COMPLEX, comm ) 
    );
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
    SafeMpi( 
        MPI_Alltoall( const_cast<dcomplex*>(sbuf), sc, MPI_DOUBLE_COMPLEX,
                      rbuf, rc, MPI_DOUBLE_COMPLEX, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

inline void
elemental::wrappers::mpi::AllToAll
( const int* sbuf, const int* scs, const int* sds, 
        int* rbuf, const int* rcs, const int* rds, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::AllToAll");
#endif
    SafeMpi( 
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
    SafeMpi( 
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
    SafeMpi( 
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
    SafeMpi( 
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
    SafeMpi( 
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
#endif // WITHOUT_COMPLEX

inline void
elemental::wrappers::mpi::Reduce
( const int* sbuf, int* rbuf, int count, MPI_Op op, int root, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::Reduce");
#endif
    SafeMpi( 
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
    SafeMpi( 
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
    SafeMpi( 
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
    SafeMpi( 
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
    SafeMpi( 
        MPI_Reduce( const_cast<dcomplex*>(sbuf), 
                    rbuf, count, MPI_DOUBLE_COMPLEX, op, root, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

inline void
elemental::wrappers::mpi::AllReduce
( const char* sbuf, char* rbuf, int count, MPI_Op op, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::AllReduce");
#endif
    SafeMpi( 
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
    SafeMpi( 
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
    SafeMpi( 
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
    SafeMpi( 
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
    SafeMpi( 
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
    SafeMpi( 
        MPI_Allreduce( const_cast<dcomplex*>(sbuf), 
                       rbuf, count, MPI_DOUBLE_COMPLEX, op, comm )
    );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

inline void
elemental::wrappers::mpi::ReduceScatter
( const int* sbuf, int* rbuf, const int* rcs, MPI_Op op, MPI_Comm comm )
{
#ifndef RELEASE
    PushCallStack("wrappers::mpi::ReduceScatter");
#endif
    SafeMpi( 
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
    SafeMpi( 
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
    SafeMpi( 
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
    SafeMpi( 
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
    SafeMpi( 
        MPI_Reduce_scatter( const_cast<dcomplex*>(sbuf), 
                            rbuf, 
                            const_cast<int*>(rcs), 
                            MPI_DOUBLE_COMPLEX, op, comm ) 
    );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

#endif /* ELEMENTAL_WRAPPERS_MPI_HPP */

