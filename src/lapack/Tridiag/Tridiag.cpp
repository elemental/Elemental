/*
   Copyright (c) 2009-2010, Jack Poulson
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
#include "elemental/lapack_internal.hpp"
using namespace elemental;
using namespace std;

template<typename R>
void
elemental::lapack::Tridiag
( Shape shape, DistMatrix<R,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("lapack::Tridiag");
#endif
    // The old approach was to use a general (i.e., nonsquare) process grid,
    // but it is usually significantly faster to redistribute to a (smaller)
    // square process grid, tridiagonalize, and then redistribute back to the
    // full grid.
    if( shape == Lower )
        lapack::internal::TridiagL( A );
    else
        lapack::internal::TridiagU( A );

    // This will be enabled as soon as the redistributions are written
    /*
    const Grid& g = A.Grid();
    int p = g.Size();
    int pSqrt = static_cast<int>(sqrt(static_cast<double>(p)));

    if( g.InGrid() )
    {
        MPI_Comm owningComm = g.OwningComm();
        MPI_Group owningGroup;
        MPI_Comm_group( owningComm, &owningGroup );
        
        std::vector<int> ranks(pSqrt*pSqrt);
        for( int i=0; i<ranks.size(); ++i )
            ranks[i] = i;
        MPI_Group squareGroup;
        MPI_Group_incl( owningGroup, pSqrt*pSqrt, &ranks[0], &squareGroup );

        const Grid squareGrid( owningComm, squareGroup, pSqrt, pSqrt );
        bool inSquareGroup = squareGrid.InGrid();

        // Determine padding size

        DistMatrix<R,MC,MR> ASquare(m,m,squareGrid);

        // Redistribute to from A to ASquare

        if( shape == Lower )        
            lapack::internal::TridiagLSquare( ASquare, paddingSize );
        else
            lapack::internal::TridiagUSquare( ASquare, paddingSize );

        // Redistribute from ASquare to A
    }
    */

#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename R>
void
elemental::lapack::Tridiag
( Shape shape, 
  DistMatrix<complex<R>,MC,  MR  >& A,
  DistMatrix<complex<R>,Star,Star>& t )
{
#ifndef RELEASE
    PushCallStack("lapack::Tridiag");
#endif
    // The old approach was to use a general (i.e., nonsquare) process grid,
    // but it is usually significantly faster to redistribute to a (smaller)
    // square process grid, tridiagonalize, and then redistribute back to the
    // full grid.
    if( shape == Lower )
        lapack::internal::TridiagL( A, t );
    else
        lapack::internal::TridiagU( A, t );

    // This will be enabled as soon as the redistributions are written
    /*
    const Grid& g = A.Grid();
    int p = g.Size();
    int pSqrt = static_cast<int>(sqrt(static_cast<double>(p)));

    if( g.InGrid() )
    {
        MPI_Comm owningComm = g.OwningComm();
        MPI_Group owningGroup;
        MPI_Comm_group( owningComm, &owningGroup );
        
        std::vector<int> ranks(pSqrt*pSqrt);
        for( int i=0; i<ranks.size(); ++i )
            ranks[i] = i;
        MPI_Group squareGroup;
        MPI_Group_incl( owningGroup, pSqrt*pSqrt, &ranks[0], &squareGroup );

        const Grid squareGrid( owningComm, squareGroup, pSqrt, pSqrt );
        bool inSquareGroup = squareGrid.InGrid();

        // Determine padding size

        DistMatrix<complex<R>,MC,MR> ASquare(m,m,squareGrid);
        DistMatrix<complex<R>,Star,Star> tSquare(squareGrid);

        // Redistribute to from A to ASquare

        if( shape == Lower )        
            lapack::internal::TridiagLSquare( ASquare, tSquare, paddingSize );
        else
            lapack::internal::TridiagUSquare( ASquare, tSquare, paddingSize );

        // Redistribute from ASquare to A
        // Redistribute from tSquare to t
    }
    */
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

template void elemental::lapack::Tridiag
( Shape shape, 
  DistMatrix<float,MC,MR>& A );

template void elemental::lapack::Tridiag
( Shape shape, 
  DistMatrix<double,MC,MR>& A );

#ifndef WITHOUT_COMPLEX
template void elemental::lapack::Tridiag
( Shape shape,
  DistMatrix<scomplex,MC,  MR  >& A,
  DistMatrix<scomplex,Star,Star>& t );

template void elemental::lapack::Tridiag
( Shape shape,
  DistMatrix<dcomplex,MC,  MR  >& A,
  DistMatrix<dcomplex,Star,Star>& t );
#endif

