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

template<typename R> // representation of a real number
inline void
elemental::advanced::HouseholderSolve
( Orientation orientation, DistMatrix<R,MC,MR>& A, DistMatrix<R,MC,MR>& X )
{
#ifndef RELEASE
    PushCallStack("advanced::HouseholderSolve");
    if( A.Grid() != X.Grid() )
        throw std::logic_error("Grids do not match");
#endif
    // TODO: Add scaling
    const int m = A.Height();
    const int n = A.Width();
    if( orientation == NORMAL )
    {
        if( m >= n )
        {
            if( m != X.Height() )
                throw std::logic_error("A and X do not conform");

            advanced::QR( A );
            // Apply Q' to X
            advanced::ApplyPackedReflectors
            ( LEFT, LOWER, VERTICAL, FORWARD, 0, A, X );
            // Shrink X to its new height
            X.Resize( n, X.Width() );
            // Solve against R (checking for singularities)
            DistMatrix<R,MC,MR> AT;
            AT.LockedView( A, 0, 0, n, n );
            basic::Trsm( LEFT, UPPER, NORMAL, NON_UNIT, (R)1, AT, X, true );
        }
        else
        {
            if( n != X.Height() )
                throw std::logic_error
                ("X should be passed in with padding for the solution");

            DistMatrix<R,MC,MR> XT,
                                XB;
            PartitionDown( X, XT,
                              XB, m );

            advanced::LQ( A );
            // Solve against L (checking for singularities)
            DistMatrix<R,MC,MR> AL;
            AL.LockedView( A, 0, 0, m, m );
            basic::Trsm( LEFT, LOWER, NORMAL, NON_UNIT, (R)1, AL, XT, true );
            // Apply Q' to X (explicitly zero the bottom of X first)
            XB.SetToZero();
            advanced::ApplyPackedReflectors
            ( LEFT, UPPER, HORIZONTAL, FORWARD, 0, A, X );
        }
    }
    else // orientation == ADJOINT
    {
        if( m >= n )
        {
            if( m != X.Height() )
                throw std::logic_error
                ("X should be passed in with padding for the solution");

            DistMatrix<R,MC,MR> XT,
                                XB;
            PartitionDown( X, XT,
                              XB, n );

            advanced::QR( A );
            // Solve against R' (checking for singularities)
            DistMatrix<R,MC,MR> AT;
            AT.LockedView( A, 0, 0, n, n );
            basic::Trsm
            ( LEFT, UPPER, ADJOINT, NON_UNIT, (R)1, AT, XT, true );
            // Apply Q to X (explicitly zero the bottom of X first)
            XB.SetToZero();
            advanced::ApplyPackedReflectors
            ( LEFT, LOWER, VERTICAL, BACKWARD, 0, A, X );
        }
        else
        {
            if( n != X.Height() )
                throw std::logic_error("A and X do not conform");

            advanced::LQ( A );
            // Apply Q to X
            advanced::ApplyPackedReflectors
            ( LEFT, UPPER, HORIZONTAL, BACKWARD, 0, A, X );
            X.Resize( m, X.Width() );
            // Solve against L' (check for singularities)
            DistMatrix<R,MC,MR> AL;
            AL.LockedView( A, 0, 0, m, m );
            basic::Trsm
            ( LEFT, LOWER, ADJOINT, NON_UNIT, (R)1, AL, X, true );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename R> // representation of a real number
inline void
elemental::advanced::HouseholderSolve
( Orientation orientation, 
  DistMatrix<std::complex<R>,MC,MR>& A, 
  DistMatrix<std::complex<R>,MC,MR>& X )
{
#ifndef RELEASE
    PushCallStack("advanced::HouseholderSolve");
    if( A.Grid() != X.Grid() )
        throw std::logic_error("Grids do not match");
    if( orientation == TRANSPOSE )
        throw std::logic_error("Invalid orientation");
#endif
    // TODO: Add scaling
    typedef std::complex<R> C;
    const int m = A.Height();
    const int n = A.Width();
    DistMatrix<C,MD,STAR> t;
    if( orientation == NORMAL )
    {
        if( m >= n )
        {
            if( m != X.Height() )
                throw std::logic_error("A and X do not conform");

            advanced::QR( A, t );
            // Apply Q' to X
            advanced::ApplyPackedReflectors
            ( LEFT, LOWER, VERTICAL, FORWARD, CONJUGATED, 0, A, t, X );
            // Shrink X to its new height
            X.Resize( n, X.Width() );
            // Solve against R (checking for singularities)
            DistMatrix<C,MC,MR> AT;
            AT.LockedView( A, 0, 0, n, n );
            basic::Trsm( LEFT, UPPER, NORMAL, NON_UNIT, (C)1, AT, X, true );
        }
        else
        {
            if( n != X.Height() )
                throw std::logic_error
                ("X should be passed in with padding for the solution");

            DistMatrix<C,MC,MR> XT,
                                XB;
            PartitionDown( X, XT,
                              XB, m );

            advanced::LQ( A, t );
            // Solve against L (checking for singularities)
            DistMatrix<C,MC,MR> AL;
            AL.LockedView( A, 0, 0, m, m );
            basic::Trsm( LEFT, LOWER, NORMAL, NON_UNIT, (C)1, AL, XT, true );
            // Apply Q' to X (explicitly zero the bottom of X first)
            XB.SetToZero();
            advanced::ApplyPackedReflectors
            ( LEFT, UPPER, HORIZONTAL, FORWARD, CONJUGATED, 0, A, t, X );
        }
    }
    else // orientation == ADJOINT
    {
        if( m >= n )
        {
            if( m != X.Height() )
                throw std::logic_error
                ("X should be passed in with padding for the solution");

            DistMatrix<C,MC,MR> XT,
                                XB;
            PartitionDown( X, XT,
                              XB, n );

            advanced::QR( A, t );
            // Solve against R' (checking for singularities)
            DistMatrix<C,MC,MR> AT;
            AT.LockedView( A, 0, 0, n, n );
            basic::Trsm
            ( LEFT, UPPER, ADJOINT, NON_UNIT, (C)1, AT, XT, true );
            // Apply Q to X (explicitly zero the bottom of X first)
            XB.SetToZero();
            advanced::ApplyPackedReflectors
            ( LEFT, LOWER, VERTICAL, BACKWARD, CONJUGATED, 0, A, t, X );
        }
        else
        {
            if( n != X.Height() )
                throw std::logic_error("A and X do not conform");

            advanced::LQ( A, t );
            // Apply Q to X
            advanced::ApplyPackedReflectors
            ( LEFT, UPPER, HORIZONTAL, BACKWARD, CONJUGATED, 0, A, t, X );
            X.Resize( m, X.Width() );
            // Solve against L' (check for singularities)
            DistMatrix<C,MC,MR> AL;
            AL.LockedView( A, 0, 0, m, m );
            basic::Trsm
            ( LEFT, LOWER, ADJOINT, NON_UNIT, (C)1, AL, X, true );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX
