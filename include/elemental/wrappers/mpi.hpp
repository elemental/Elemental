/*
   Copyright (c) 2009-2011, Jack Poulson
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
#ifndef ELEMENTAL_WRAPPERS_MPI_HPP
#define ELEMENTAL_WRAPPERS_MPI_HPP 1

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
Wait( MPI_Request& request );

bool
CongruentComms( MPI_Comm comm1, MPI_Comm comm2 );

void
Send
( const int* buf, int count, int to, int tag, MPI_Comm comm );
void
ISend
( const int* buf, int count, int to, int tag, MPI_Comm comm, 
  MPI_Request& request );

void
Send
( const float* buf, int count, int to, int tag, MPI_Comm comm );
void
ISend
( const float* buf, int count, int to, int tag, MPI_Comm comm,
  MPI_Request& request );

void
Send
( const double* buf, int count, int to, int tag, MPI_Comm comm );
void
ISend
( const double* buf, int count, int to, int tag, MPI_Comm comm,
  MPI_Request& request );

#ifndef WITHOUT_COMPLEX
void
Send
( const scomplex* buf, int count, int to, int tag, MPI_Comm comm );
void
ISend
( const scomplex* buf, int count, int to, int tag, MPI_Comm comm,
  MPI_Request& request );

void
Send
( const dcomplex* buf, int count, int to, int tag, MPI_Comm comm );
void
ISend
( const dcomplex* buf, int count, int to, int tag, MPI_Comm comm,
  MPI_Request& request );
#endif
        
void
Recv
( int* buf, int count, int from, int tag, MPI_Comm comm );
void
IRecv
( int* buf, int count, int from, int tag, MPI_Comm comm,
  MPI_Request& request );

void
Recv
( float* buf, int count, int from, int tag, MPI_Comm comm );
void
IRecv
( float* buf, int count, int from, int tag, MPI_Comm comm,
  MPI_Request& request );

void
Recv
( double* buf, int count, int from, int tag, MPI_Comm comm );
void
IRecv
( double* buf, int count, int from, int tag, MPI_Comm comm,
  MPI_Request& request );

#ifndef WITHOUT_COMPLEX
void
Recv
( scomplex* buf, int count, int from, int tag, MPI_Comm comm );
void
IRecv
( scomplex* buf, int count, int from, int tag, MPI_Comm comm,
  MPI_Request& request );

void
Recv
( dcomplex* buf, int count, int from, int tag, MPI_Comm comm );
void
IRecv
( dcomplex* buf, int count, int from, int tag, MPI_Comm comm,
  MPI_Request& request );
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

#endif /* ELEMENTAL_WRAPPERS_MPI_HPP */

