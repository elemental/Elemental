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

namespace elem {

template<typename R>
inline void
HouseholderSolve
( Orientation orientation, DistMatrix<R,MC,MR>& A, DistMatrix<R,MC,MR>& B )
{
#ifndef RELEASE
    PushCallStack("HouseholderSolve");
    if( A.Grid() != B.Grid() )
        throw std::logic_error("Grids do not match");
#endif
    // TODO: Add scaling
    const int m = A.Height();
    const int n = A.Width();
    if( orientation == NORMAL )
    {
        if( m >= n )
        {
            if( m != B.Height() )
                throw std::logic_error("A and B do not conform");

            QR( A );
            // Apply Q' to B
            ApplyPackedReflectors
            ( LEFT, LOWER, VERTICAL, FORWARD, 0, A, B );
            // Shrink B to its new height
            B.ResizeTo( n, B.Width() );
            // Solve against R (checking for singularities)
            DistMatrix<R,MC,MR> AT;
            AT.LockedView( A, 0, 0, n, n );
            Trsm( LEFT, UPPER, NORMAL, NON_UNIT, (R)1, AT, B, true );
        }
        else
        {
            if( n != B.Height() )
                throw std::logic_error
                ("B should be passed in with padding for the solution");

            DistMatrix<R,MC,MR> BT,
                                BB;
            PartitionDown( B, BT,
                              BB, m );

            LQ( A );
            // Solve against L (checking for singularities)
            DistMatrix<R,MC,MR> AL;
            AL.LockedView( A, 0, 0, m, m );
            Trsm( LEFT, LOWER, NORMAL, NON_UNIT, (R)1, AL, BT, true );
            // Apply Q' to B (explicitly zero the bottom of B first)
            Zero( BB );
            ApplyPackedReflectors( LEFT, UPPER, HORIZONTAL, FORWARD, 0, A, B );
        }
    }
    else // orientation == ADJOINT
    {
        if( m >= n )
        {
            if( m != B.Height() )
                throw std::logic_error
                ("B should be passed in with padding for the solution");

            DistMatrix<R,MC,MR> BT,
                                BB;
            PartitionDown( B, BT,
                              BB, n );

            QR( A );
            // Solve against R' (checking for singularities)
            DistMatrix<R,MC,MR> AT;
            AT.LockedView( A, 0, 0, n, n );
            Trsm( LEFT, UPPER, ADJOINT, NON_UNIT, (R)1, AT, BT, true );
            // Apply Q to B (explicitly zero the bottom of B first)
            Zero( BB );
            ApplyPackedReflectors( LEFT, LOWER, VERTICAL, BACKWARD, 0, A, B );
        }
        else
        {
            if( n != B.Height() )
                throw std::logic_error("A and B do not conform");

            LQ( A );
            // Apply Q to B
            ApplyPackedReflectors( LEFT, UPPER, HORIZONTAL, BACKWARD, 0, A, B );
            B.ResizeTo( m, B.Width() );
            // Solve against L' (check for singularities)
            DistMatrix<R,MC,MR> AL;
            AL.LockedView( A, 0, 0, m, m );
            Trsm( LEFT, LOWER, ADJOINT, NON_UNIT, (R)1, AL, B, true );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R> 
inline void
HouseholderSolve
( Orientation orientation, 
  DistMatrix<Complex<R>,MC,MR>& A, 
  DistMatrix<Complex<R>,MC,MR>& B )
{
#ifndef RELEASE
    PushCallStack("HouseholderSolve");
    if( A.Grid() != B.Grid() )
        throw std::logic_error("Grids do not match");
    if( orientation == TRANSPOSE )
        throw std::logic_error("Invalid orientation");
#endif
    // TODO: Add scaling
    typedef Complex<R> C;
    const int m = A.Height();
    const int n = A.Width();
    DistMatrix<C,MD,STAR> t;
    if( orientation == NORMAL )
    {
        if( m >= n )
        {
            if( m != B.Height() )
                throw std::logic_error("A and B do not conform");

            QR( A, t );
            // Apply Q' to B
            ApplyPackedReflectors
            ( LEFT, LOWER, VERTICAL, FORWARD, CONJUGATED, 0, A, t, B );
            // Shrink B to its new height
            B.ResizeTo( n, B.Width() );
            // Solve against R (checking for singularities)
            DistMatrix<C,MC,MR> AT;
            AT.LockedView( A, 0, 0, n, n );
            Trsm( LEFT, UPPER, NORMAL, NON_UNIT, (C)1, AT, B, true );
        }
        else
        {
            if( n != B.Height() )
                throw std::logic_error
                ("B should be passed in with padding for the solution");

            DistMatrix<C,MC,MR> BT,
                                BB;
            PartitionDown( B, BT,
                              BB, m );

            LQ( A, t );
            // Solve against L (checking for singularities)
            DistMatrix<C,MC,MR> AL;
            AL.LockedView( A, 0, 0, m, m );
            Trsm( LEFT, LOWER, NORMAL, NON_UNIT, (C)1, AL, BT, true );
            // Apply Q' to B (explicitly zero the bottom of B first)
            Zero( BB );
            ApplyPackedReflectors
            ( LEFT, UPPER, HORIZONTAL, FORWARD, CONJUGATED, 0, A, t, B );
        }
    }
    else // orientation == ADJOINT
    {
        if( m >= n )
        {
            if( m != B.Height() )
                throw std::logic_error
                ("B should be passed in with padding for the solution");

            DistMatrix<C,MC,MR> BT,
                                BB;
            PartitionDown( B, BT,
                              BB, n );

            QR( A, t );
            // Solve against R' (checking for singularities)
            DistMatrix<C,MC,MR> AT;
            AT.LockedView( A, 0, 0, n, n );
            Trsm( LEFT, UPPER, ADJOINT, NON_UNIT, (C)1, AT, BT, true );
            // Apply Q to B (explicitly zero the bottom of B first)
            Zero( BB );
            ApplyPackedReflectors
            ( LEFT, LOWER, VERTICAL, BACKWARD, CONJUGATED, 0, A, t, B );
        }
        else
        {
            if( n != B.Height() )
                throw std::logic_error("A and B do not conform");

            LQ( A, t );
            // Apply Q to B
            ApplyPackedReflectors
            ( LEFT, UPPER, HORIZONTAL, BACKWARD, CONJUGATED, 0, A, t, B );
            B.ResizeTo( m, B.Width() );
            // Solve against L' (check for singularities)
            DistMatrix<C,MC,MR> AL;
            AL.LockedView( A, 0, 0, m, m );
            Trsm( LEFT, LOWER, ADJOINT, NON_UNIT, (C)1, AL, B, true );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
