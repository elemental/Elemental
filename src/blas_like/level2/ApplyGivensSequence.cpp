/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>
#include <El/blas_like/level2.hpp>

namespace El {

// [CITATION] LAPACK's {s,d,c,z}lasr

// TODO: Optimized versions which avoid temporaries and/or directly work on the
// underlying raw data buffers

template<typename F,typename=DisableIf<IsReal<F>>>
void ApplyVariableLeft
( Int i,
  const Base<F>& c,
  const F& s,
        Matrix<F>& A,
        F& tmp,
  const Base<F>& one=Base<F>(1),
  const F& zeroF=F(0) )
{
    if( c == one && s == zeroF )
        return;
    const Int n = A.Width();
    const F sConj = Conj(s);
    for( Int j=0; j<n; ++j )
    {
        tmp = A(i+1,j);
        A(i+1,j) = c*tmp - sConj*A(i,j);
        A(i,  j) = s*tmp +     c*A(i,j);
    }
}

template<typename F>
void ApplyVariableLeft
( Int i,
  const Base<F>& c,
  const Base<F>& s,
        Matrix<F>& A,
        F& tmp,
  const Base<F>& one=Base<F>(1),
  const Base<F>& zero=Base<F>(0) )
{
    if( c == one && s == zero )
        return;
    const Int n = A.Width();
    for( Int j=0; j<n; ++j )
    {
        tmp = A(i+1,j);
        A(i+1,j) = c*tmp - s*A(i,j);
        A(i,  j) = s*tmp + c*A(i,j);
    }
}

template<typename F,typename=DisableIf<IsReal<F>>>
void ApplyVariableRight
( Int j,
  const Base<F>& c,
  const F& s,
        Matrix<F>& A,
        F& tmp,
  const Base<F>& one=Base<F>(1),
  const F& zeroF=F(0) )
{
    if( c == one && s == zeroF )
        return;
    const Int m = A.Height();
    const F sConj = Conj(s);
    for( Int i=0; i<m; ++i )
    {
        tmp = A(i,j+1);
        A(i,j+1) = c*tmp - sConj*A(i,j);
        A(i,j  ) = s*tmp +     c*A(i,j);
    }
}

template<typename F>
void ApplyVariableRight
( Int j,
  const Base<F>& c,
  const Base<F>& s,
        Matrix<F>& A,
        F& tmp,
  const Base<F>& one=Base<F>(1),
  const Base<F>& zero=Base<F>(0) )
{
    if( c == one && s == zero )
        return;
    const Int m = A.Height();
    for( Int i=0; i<m; ++i )
    {
        tmp = A(i,j+1);
        A(i,j+1) = c*tmp - s*A(i,j);
        A(i,j  ) = s*tmp + c*A(i,j);
    }
}

template<typename F,typename=DisableIf<IsReal<F>>>
void ApplyTopLeft
( Int i,
  const Base<F>& c,
  const F& s,
        Matrix<F>& A,
        F& tmp,
  const Base<F>& one=Base<F>(1),
  const F& zeroF=F(0) )
{
    if( c == one && s == zeroF )
        return;
    const Int n = A.Width();
    const F sConj = Conj(s);
    for( Int j=0; j<n; ++j )
    {
        tmp = A(i,j);
        A(i,j) = c*tmp - sConj*A(0,j);
        A(0,j) = s*tmp +     c*A(0,j);
    }
}

template<typename F>
void ApplyTopLeft
( Int i,
  const Base<F>& c,
  const Base<F>& s,
        Matrix<F>& A,
        F& tmp,
  const Base<F>& one=Base<F>(1),
  const Base<F>& zero=Base<F>(0) )
{
    if( c == one && s == zero )
        return;
    const Int n = A.Width();
    for( Int j=0; j<n; ++j )
    {
        tmp = A(i,j);
        A(i,j) = c*tmp - s*A(0,j);
        A(0,j) = s*tmp + c*A(0,j);
    }
}

template<typename F,typename=DisableIf<IsReal<F>>>
void ApplyTopRight
( Int j,
  const Base<F>& c,
  const F& s,
        Matrix<F>& A,
        F& tmp,
  const Base<F>& one=Base<F>(1),
  const F& zeroF=F(0) )
{
    if( c == one && s == zeroF )
        return;
    const Int m = A.Height();
    const F sConj = Conj(s);
    for( Int i=0; i<m; ++i )
    {
        tmp = A(i,j);
        A(i,j) = c*tmp - sConj*A(i,0);
        A(i,0) = s*tmp +     c*A(i,0);
    }
}

template<typename F>
void ApplyTopRight
( Int j,
  const Base<F>& c,
  const Base<F>& s,
        Matrix<F>& A,
        F& tmp,
  const Base<F>& one=Base<F>(1),
  const Base<F>& zero=Base<F>(0) )
{
    if( c == one && s == zero )
        return;
    const Int m = A.Height();
    for( Int i=0; i<m; ++i )
    {
        tmp = A(i,j);
        A(i,j) = c*tmp - s*A(i,0);
        A(i,0) = s*tmp + c*A(i,0);
    }
}

template<typename F,typename=DisableIf<IsReal<F>>>
void ApplyBottomLeft
( Int i,
  const Base<F>& c,
  const F& s,
        Matrix<F>& A,
        F& tmp,
  const Base<F>& one=Base<F>(1),
  const F& zeroF=F(0) )
{
    if( c == one && s == zeroF )
        return;
    const Int m = A.Height();
    const Int n = A.Width();
    const F sConj = Conj(s);
    for( Int j=0; j<n; ++j )
    {
        tmp = A(i,j);
        A(i,  j) = s*A(m-1,j) +     c*tmp;
        A(m-1,j) = c*A(m-1,j) - sConj*tmp;
    }
}

template<typename F>
void ApplyBottomLeft
( Int i,
  const Base<F>& c,
  const Base<F>& s,
        Matrix<F>& A,
        F& tmp,
  const Base<F>& one=Base<F>(1),
  const Base<F>& zero=Base<F>(0) )
{
    if( c == one && s == zero )
        return;
    const Int m = A.Height();
    const Int n = A.Width();
    for( Int j=0; j<n; ++j )
    {
        tmp = A(i,j);
        A(i,  j) = s*A(m-1,j) + c*tmp;
        A(m-1,j) = c*A(m-1,j) - s*tmp;
    }
}

template<typename F,typename=DisableIf<IsReal<F>>>
void ApplyBottomRight
( Int j,
  const Base<F>& c,
  const F& s,
        Matrix<F>& A,
        F& tmp,
  const Base<F>& one=Base<F>(1),
  const F& zeroF=F(0) )
{
    if( c == one && s == zeroF )
        return;
    const Int m = A.Height();
    const Int n = A.Width();
    const F sConj = Conj(s);
    for( Int i=0; i<m; ++i )
    {
        tmp = A(i,j);
        A(i,j  ) = s*A(i,n-1) +     c*tmp;
        A(i,n-1) = c*A(i,n-1) - sConj*tmp;
    }
}

template<typename F>
void ApplyBottomRight
( Int j,
  const Base<F>& c,
  const Base<F>& s,
        Matrix<F>& A,
        F& tmp,
  const Base<F>& one=Base<F>(1),
  const Base<F>& zero=Base<F>(0) )
{
    if( c == one && s == zero )
        return;
    const Int m = A.Height();
    const Int n = A.Width();
    for( Int i=0; i<m; ++i )
    {
        tmp = A(i,j);
        A(i,j  ) = s*A(i,n-1) + c*tmp;
        A(i,n-1) = c*A(i,n-1) - s*tmp;
    }
}

template<typename F,typename>
void ApplyGivensSequence
( LeftOrRight side, GivensSequenceType seqType, ForwardOrBackward direction,
  const Matrix<Base<F>>& cList,
  const Matrix<F>& sList,
  Matrix<F>& A )
{
    EL_DEBUG_CSE
    // TODO(poulson): Assert the correct lengths of cList and sList
    typedef Base<F> Real;
    const Real one(1);
    const F zeroF(0);
    const Int m = A.Height();    
    const Int n = A.Width();
    if( m == 0 || n == 0 )
        return;

    F tmp;
    if( side == LEFT )
    {
        if( seqType == VARIABLE_GIVENS_SEQUENCE )
        {
            if( direction == FORWARD )
            {
                for( Int i=0; i<m-1; ++i )
                {
                    ApplyVariableLeft
                    ( i, cList(i), sList(i), A, tmp, one, zeroF );
                }
            }
            else
            {
                for( Int i=m-2; i>=0; --i )
                {
                    ApplyVariableLeft
                    ( i, cList(i), sList(i), A, tmp, one, zeroF );
                }
            }
        }
        else if( seqType == TOP_GIVENS_SEQUENCE )
        {
            if( direction == FORWARD )
            {
                for( Int i=1; i<m; ++i )
                {
                    ApplyTopLeft
                    ( i, cList(i-1), sList(i-1), A, tmp, one, zeroF );
                }
            }
            else
            {
                for( Int i=m-1; i>=1; --i )
                {
                    ApplyTopLeft
                    ( i, cList(i-1), sList(i-1), A, tmp, one, zeroF );
                }
            }
        }
        else
        {
            if( direction == FORWARD )
            {
                for( Int i=0; i<m-1; ++i )
                {
                    ApplyBottomLeft
                    ( i, cList(i), sList(i), A, tmp, one, zeroF );
                }
            }
            else
            {
                for( Int i=m-2; i>=0; --i )
                {
                    ApplyBottomLeft
                    ( i, cList(i), sList(i), A, tmp, one, zeroF );
                }
            }
        }
    }
    else
    {
        if( seqType == VARIABLE_GIVENS_SEQUENCE )
        {
            if( direction == FORWARD )
            {
                for( Int j=0; j<n-1; ++j )
                {
                    ApplyVariableRight
                    ( j, cList(j), sList(j), A, tmp, one, zeroF );
                }
            }
            else
            {
                for( Int j=n-2; j>=0; --j )
                {
                    ApplyVariableRight
                    ( j, cList(j), sList(j), A, tmp, one, zeroF );
                }
            }
        }
        else if( seqType == TOP_GIVENS_SEQUENCE )
        {
            if( direction == FORWARD )
            {
                for( Int j=1; j<n; ++j )
                {
                    ApplyTopRight
                    ( j, cList(j-1), sList(j-1), A, tmp, one, zeroF );
                }
            }
            else
            {
                for( Int j=n-1; j>=1; --j )
                {
                    ApplyTopRight
                    ( j, cList(j-1), sList(j-1), A, tmp, one, zeroF );
                }
            }
        }
        else
        {
            if( direction == FORWARD )
            {
                for( Int j=0; j<n-1; ++j )
                {
                    ApplyBottomRight
                    ( j, cList(j), sList(j), A, tmp, one, zeroF );
                }
            }
            else
            {
                for( Int j=n-2; j>=0; --j )
                {
                    ApplyBottomRight
                    ( j, cList(j), sList(j), A, tmp, one, zeroF );
                }
            }
        }
    }
}

template<typename F>
void ApplyGivensSequence
( LeftOrRight side, GivensSequenceType seqType, ForwardOrBackward direction,
  const Matrix<Base<F>>& cList,
  const Matrix<Base<F>>& sList,
  Matrix<F>& A )
{
    EL_DEBUG_CSE
    // TODO(poulson): Assert the correct lengths of cList and sList
    typedef Base<F> Real;
    const Real one(1), zero(0);
    const Int m = A.Height();    
    const Int n = A.Width();
    if( m == 0 || n == 0 )
        return;

    F tmp;
    if( side == LEFT )
    {
        if( seqType == VARIABLE_GIVENS_SEQUENCE )
        {
            if( direction == FORWARD )
            {
                for( Int i=0; i<m-1; ++i )
                {
                    ApplyVariableLeft
                    ( i, cList(i), sList(i), A, tmp, one, zero );
                }
            }
            else
            {
                for( Int i=m-2; i>=0; --i )
                {
                    ApplyVariableLeft
                    ( i, cList(i), sList(i), A, tmp, one, zero );
                }
            }
        }
        else if( seqType == TOP_GIVENS_SEQUENCE )
        {
            if( direction == FORWARD )
            {
                for( Int i=1; i<m; ++i )
                {
                    ApplyTopLeft
                    ( i, cList(i-1), sList(i-1), A, tmp, one, zero );
                }
            }
            else
            {
                for( Int i=m-1; i>=1; --i )
                {
                    ApplyTopLeft
                    ( i, cList(i-1), sList(i-1), A, tmp, one, zero );
                }
            }
        }
        else
        {
            if( direction == FORWARD )
            {
                for( Int i=0; i<m-1; ++i )
                {
                    ApplyBottomLeft
                    ( i, cList(i), sList(i), A, tmp, one, zero );
                }
            }
            else
            {
                for( Int i=m-2; i>=0; --i )
                {
                    ApplyBottomLeft
                    ( i, cList(i), sList(i), A, tmp, one, zero );
                }
            }
        }
    }
    else
    {
        if( seqType == VARIABLE_GIVENS_SEQUENCE )
        {
            if( direction == FORWARD )
            {
                for( Int j=0; j<n-1; ++j )
                {
                    ApplyVariableRight
                    ( j, cList(j), sList(j), A, tmp, one, zero );
                }
            }
            else
            {
                for( Int j=n-2; j>=0; --j )
                {
                    ApplyVariableRight
                    ( j, cList(j), sList(j), A, tmp, one, zero );
                }
            }
        }
        else if( seqType == TOP_GIVENS_SEQUENCE )
        {
            if( direction == FORWARD )
            {
                for( Int j=1; j<n; ++j )
                {
                    ApplyTopRight
                    ( j, cList(j-1), sList(j-1), A, tmp, one, zero );
                }
            }
            else
            {
                for( Int j=n-1; j>=1; --j )
                {
                    ApplyTopRight
                    ( j, cList(j-1), sList(j-1), A, tmp, one, zero );
                }
            }
        }
        else
        {
            if( direction == FORWARD )
            {
                for( Int j=0; j<n-1; ++j )
                {
                    ApplyBottomRight
                    ( j, cList(j), sList(j), A, tmp, one, zero );
                }
            }
            else
            {
                for( Int j=n-2; j>=0; --j )
                {
                    ApplyBottomRight
                    ( j, cList(j), sList(j), A, tmp, one, zero );
                }
            }
        }
    }
}

#define PROTO_REAL(F) \
  template void ApplyGivensSequence \
  ( LeftOrRight side, GivensSequenceType seqType, ForwardOrBackward direction, \
    const Matrix<Base<F>>& cList, \
    const Matrix<Base<F>>& sList, \
    Matrix<F>& A );

#define PROTO(F) \
  PROTO_REAL(F) \
  template void ApplyGivensSequence \
  ( LeftOrRight side, GivensSequenceType seqType, ForwardOrBackward direction, \
    const Matrix<Base<F>>& cList, \
    const Matrix<F>& sList, \
    Matrix<F>& A );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
