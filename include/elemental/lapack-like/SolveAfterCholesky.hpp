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
SolveAfterCholesky
( UpperOrLower uplo, Orientation orientation, 
  const Matrix<F>& A, Matrix<F>& B )
{
#ifndef RELEASE
    PushCallStack("SolveAfterCholesky");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( A.Height() != B.Height() )
        throw std::logic_error("A and B must be the same height");
#endif
    if( B.Width() == 1 )
    {
        if( uplo == LOWER )
        {
            if( orientation == TRANSPOSE )
                Conj( B );
            Trsv( LOWER, NORMAL, NON_UNIT, A, B );
            Trsv( LOWER, ADJOINT, NON_UNIT, A, B );
            if( orientation == TRANSPOSE )
                Conj( B );
        }
        else
        {
            if( orientation == TRANSPOSE )
                Conj( B );
            Trsv( UPPER, ADJOINT, NON_UNIT, A, B );
            Trsv( UPPER, NORMAL, NON_UNIT, A, B );
            if( orientation == TRANSPOSE )
                Conj( B );
        }
    }
    else
    {
        if( uplo == LOWER )
        {
            if( orientation == TRANSPOSE )
                Conj( B );
            Trsm( LEFT, LOWER, NORMAL, NON_UNIT, (F)1, A, B );
            Trsm( LEFT, LOWER, ADJOINT, NON_UNIT, (F)1, A, B );
            if( orientation == TRANSPOSE )
                Conj( B );
        }
        else
        {
            if( orientation == TRANSPOSE )
                Conj( B );
            Trsm( LEFT, UPPER, ADJOINT, NON_UNIT, (F)1, A, B );
            Trsm( LEFT, UPPER, NORMAL, NON_UNIT, (F)1, A, B );
            if( orientation == TRANSPOSE )
                Conj( B );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> 
inline void
SolveAfterCholesky
( UpperOrLower uplo, Orientation orientation, 
  const DistMatrix<F>& A, DistMatrix<F>& B )
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
        if( uplo == LOWER )
        {
            if( orientation == TRANSPOSE )
                Conj( B );
            Trsv( LOWER, NORMAL, NON_UNIT, A, B );
            Trsv( LOWER, ADJOINT, NON_UNIT, A, B );
            if( orientation == TRANSPOSE )
                Conj( B );
        }
        else
        {
            if( orientation == TRANSPOSE )
                Conj( B );
            Trsv( UPPER, ADJOINT, NON_UNIT, A, B );
            Trsv( UPPER, NORMAL, NON_UNIT, A, B );
            if( orientation == TRANSPOSE )
                Conj( B );
        }
    }
    else
    {
        if( uplo == LOWER )
        {
            if( orientation == TRANSPOSE )
                Conj( B );
            Trsm( LEFT, LOWER, NORMAL, NON_UNIT, (F)1, A, B );
            Trsm( LEFT, LOWER, ADJOINT, NON_UNIT, (F)1, A, B );
            if( orientation == TRANSPOSE )
                Conj( B );
        }
        else
        {
            if( orientation == TRANSPOSE )
                Conj( B );
            Trsm( LEFT, UPPER, ADJOINT, NON_UNIT, (F)1, A, B );
            Trsm( LEFT, UPPER, NORMAL, NON_UNIT, (F)1, A, B );
            if( orientation == TRANSPOSE )
                Conj( B );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
