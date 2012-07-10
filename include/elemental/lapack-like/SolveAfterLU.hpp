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

template<typename F> 
inline void
SolveAfterLU( Orientation orientation, const Matrix<F>& A, Matrix<F>& B )
{
#ifndef RELEASE
    PushCallStack("SolveAfterLU");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( A.Height() != B.Height() )
        throw std::logic_error("A and B must be the same height");
#endif
    if( B.Width() == 1 )
    {
        if( orientation == NORMAL )
        {
            Trsv( LOWER, NORMAL, UNIT, A, B );
            Trsv( UPPER, NORMAL, NON_UNIT, A, B );
        }
        else 
        {
            Trsv( UPPER, orientation, NON_UNIT, A, B );
            Trsv( LOWER, orientation, UNIT, A, B );
        }
    }
    else
    {
        if( orientation == NORMAL )
        {
            Trsm( LEFT, LOWER, NORMAL, UNIT, (F)1, A, B );
            Trsm( LEFT, UPPER, NORMAL, NON_UNIT, (F)1, A, B );
        }
        else
        {
            Trsm( LEFT, UPPER, orientation, NON_UNIT, (F)1, A, B );
            Trsm( LEFT, LOWER, orientation, UNIT, (F)1, A, B );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> 
inline void
SolveAfterLU( Orientation orientation, const DistMatrix<F>& A, DistMatrix<F>& B )
{
#ifndef RELEASE
    PushCallStack("SolveAfterLU");
    if( A.Grid() != B.Grid() )
        throw std::logic_error("{A,B} must be distributed over the same grid");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( A.Height() != B.Height() )
        throw std::logic_error("A and B must be the same height");
#endif
    if( B.Width() == 1 )
    {
        if( orientation == NORMAL )
        {
            Trsv( LOWER, NORMAL, UNIT, A, B );
            Trsv( UPPER, NORMAL, NON_UNIT, A, B );
        }
        else 
        {
            Trsv( UPPER, orientation, NON_UNIT, A, B );
            Trsv( LOWER, orientation, UNIT, A, B );
        }
    }
    else
    {
        if( orientation == NORMAL )
        {
            Trsm( LEFT, LOWER, NORMAL, UNIT, (F)1, A, B );
            Trsm( LEFT, UPPER, NORMAL, NON_UNIT, (F)1, A, B );
        }
        else
        {
            Trsm( LEFT, UPPER, orientation, NON_UNIT, (F)1, A, B );
            Trsm( LEFT, LOWER, orientation, UNIT, (F)1, A, B );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> 
inline void
SolveAfterLU
( Orientation orientation, 
  const Matrix<F>& A, const Matrix<int>& p, Matrix<F>& B )
{
#ifndef RELEASE
    PushCallStack("SolveAfterLU");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( A.Height() != B.Height() )
        throw std::logic_error("A and B must be the same height");
    if( p.Height() != A.Height() )
        throw std::logic_error("A and p must be the same height");
#endif
    if( B.Width() == 1 )
    {
        if( orientation == NORMAL )
        {
            ApplyRowPivots( B, p );
            Trsv( LOWER, NORMAL, UNIT, A, B );
            Trsv( UPPER, NORMAL, NON_UNIT, A, B );
        }
        else 
        {
            Trsv( UPPER, orientation, NON_UNIT, A, B );
            Trsv( LOWER, orientation, UNIT, A, B );
            ApplyInverseRowPivots( B, p );
        }
    }
    else
    {
        if( orientation == NORMAL )
        {
            ApplyRowPivots( B, p );
            Trsm( LEFT, LOWER, NORMAL, UNIT, (F)1, A, B );
            Trsm( LEFT, UPPER, NORMAL, NON_UNIT, (F)1, A, B );
        }
        else
        {
            Trsm( LEFT, UPPER, orientation, NON_UNIT, (F)1, A, B );
            Trsm( LEFT, LOWER, orientation, UNIT, (F)1, A, B );
            ApplyInverseRowPivots( B, p );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> 
inline void
SolveAfterLU
( Orientation orientation, 
  const DistMatrix<F>& A, const DistMatrix<int,VC,STAR>& p, DistMatrix<F>& B )
{
#ifndef RELEASE
    PushCallStack("SolveAfterLU");
    if( A.Grid() != B.Grid() || A.Grid() != p.Grid() )
        throw std::logic_error("{A,B} must be distributed over the same grid");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( A.Height() != B.Height() )
        throw std::logic_error("A and B must be the same height");
    if( A.Height() != p.Height() )
        throw std::logic_error("A and p must be the same height");
#endif
    if( B.Width() == 1 )
    {
        if( orientation == NORMAL )
        {
            ApplyRowPivots( B, p );
            Trsv( LOWER, NORMAL, UNIT, A, B );
            Trsv( UPPER, NORMAL, NON_UNIT, A, B );
        }
        else
        {
            Trsv( UPPER, orientation, NON_UNIT, A, B );
            Trsv( LOWER, orientation, UNIT, A, B );
            ApplyInverseRowPivots( B, p );
        }
    }
    else
    {
        if( orientation == NORMAL )
        {
            ApplyRowPivots( B, p );
            Trsm( LEFT, LOWER, NORMAL, UNIT, (F)1, A, B );
            Trsm( LEFT, UPPER, NORMAL, NON_UNIT, (F)1, A, B );
        }
        else
        {
            Trsm( LEFT, UPPER, orientation, NON_UNIT, (F)1, A, B );
            Trsm( LEFT, LOWER, orientation, UNIT, (F)1, A, B );
            ApplyInverseRowPivots( B, p );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
