/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>
#include <El/blas_like/level1.hpp>

namespace El {

template<typename T>
void Swap( Orientation orientation, Matrix<T>& X, Matrix<T>& Y )
{
    DEBUG_CSE
    const Int mX = X.Height();
    const Int nX = X.Width();

    if( orientation == NORMAL )
    {
        DEBUG_ONLY(
          if( Y.Height() != mX || Y.Width() != nX )
              LogicError("Invalid submatrix sizes");
        )
        // TODO: Optimize memory access patterns
        if( mX > nX )
        {
            for( Int j=0; j<nX; ++j )
                blas::Swap( mX, X.Buffer(0,j), 1, Y.Buffer(0,j), 1 );
        }
        else
        {
            for( Int i=0; i<mX; ++i )
                blas::Swap
                ( nX, X.Buffer(i,0), X.LDim(), Y.Buffer(i,0), Y.LDim() );
        }
    }
    else
    {
        const bool conjugate = ( orientation==ADJOINT );
        DEBUG_ONLY(
          if( Y.Width() != mX || Y.Height() != nX )
              LogicError("Invalid submatrix sizes");
        )
        // TODO: Optimize memory access patterns
        for( Int j=0; j<nX; ++j )
        {
            if( conjugate )
            {
                for( Int i=0; i<mX; ++i )
                {
                    const T alpha = X(i,j);
                    X(i,j) = Conj(Y(j,i));
                    Y(j,i) = Conj(alpha);
                }
            }
            else
            {
                blas::Swap( mX, X.Buffer(0,j), 1, Y.Buffer(j,0), Y.LDim() );
            }
        }
    }
}

template<typename T>
void Swap
( Orientation orientation,
  AbstractDistMatrix<T>& X,
  AbstractDistMatrix<T>& Y )
{
    DEBUG_CSE
    if( orientation == NORMAL )
    {
        DEBUG_ONLY(
          if( Y.Height() != X.Height() || Y.Width() != X.Width() )
              LogicError("Invalid submatrix sizes");
        )
        // TODO: Optimize communication
        unique_ptr<AbstractDistMatrix<T>> XCopy( X.Copy() );
        Copy( Y, X );
        Copy( *XCopy, Y );
    }
    else
    {
        const bool conjugate = ( orientation==ADJOINT );
        DEBUG_ONLY(
          if( Y.Width() != X.Height() || Y.Height() != X.Width() )
              LogicError("Invalid submatrix sizes");
        )
        // TODO: Optimize communication
        unique_ptr<AbstractDistMatrix<T>> XCopy( X.Copy() );       
        Transpose( Y, X, conjugate );
        Transpose( *XCopy, Y, conjugate );
    }
}

template<typename T>
void RowSwap( Matrix<T>& A, Int to, Int from )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( to < 0 || to >= A.Height() || from < 0 || from >= A.Height() )
          LogicError
          ("Attempted invalid row swap, (",to,",",from,") of matrix of height ",
           A.Height());
    )
    if( to == from )
        return;
    const Int n = A.Width();
    T* ABuf = A.Buffer();
    const Int ALDim = A.LDim();
    blas::Swap( n, &ABuf[to], ALDim, &ABuf[from], ALDim );
}

template<typename T>
void RowSwap( AbstractDistMatrix<T>& A, Int to, Int from )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( to < 0 || to >= A.Height() || from < 0 || from >= A.Height() )
          LogicError
          ("Attempted invalid row swap, (",to,",",from,") of matrix of height ",
           A.Height());
    )

    if( to == from )
        return;
    if( !A.Participating() )
        return;
    const Int nLocal = A.LocalWidth();
    const Int colAlign = A.ColAlign();
    const Int colShift = A.ColShift();
    const Int colStride = A.ColStride();
    const Int toMod = Mod(to,colStride);
    const Int fromMod = Mod(from,colStride);
    T* ABuf = A.Buffer();
    const Int ALDim = A.LDim();

    if( toMod == fromMod )
    {
        if( toMod == colShift )
        {
            const Int iLocTo = (to-colShift) / colStride;
            const Int iLocFrom = (from-colShift) / colStride;
            blas::Swap( nLocal, &ABuf[iLocTo], ALDim, &ABuf[iLocFrom], ALDim ); 
        }
    }
    else if( toMod == colShift )
    {
        const Int iLocTo = (to-colShift) / colStride;
        const int fromOwner = Mod(from+colAlign,colStride);
        vector<T> buf;
        FastResize( buf, nLocal );
        for( Int jLoc=0; jLoc<nLocal; ++jLoc )
            buf[jLoc] = ABuf[iLocTo+jLoc*ALDim];
        mpi::SendRecv( buf.data(), nLocal, fromOwner, fromOwner, A.ColComm() );
        for( Int jLoc=0; jLoc<nLocal; ++jLoc )
            ABuf[iLocTo+jLoc*ALDim] = buf[jLoc];
    }
    else if( fromMod == colShift )
    {
        const Int iLocFrom = (from-colShift) / colStride;
        const int toOwner = Mod(to+colAlign,colStride);
        vector<T> buf;
        FastResize( buf, nLocal );
        for( Int jLoc=0; jLoc<nLocal; ++jLoc )
            buf[jLoc] = ABuf[iLocFrom+jLoc*ALDim];
        mpi::SendRecv( buf.data(), nLocal, toOwner, toOwner, A.ColComm() );
        for( Int jLoc=0; jLoc<nLocal; ++jLoc )
            ABuf[iLocFrom+jLoc*ALDim] = buf[jLoc];
    }
}

template<typename T>
void ColSwap( Matrix<T>& A, Int to, Int from )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( to < 0 || to >= A.Width() || from < 0 || from >= A.Width() )
          LogicError
          ("Attempted invalid col swap, (",to,",",from,") of matrix of width ",
           A.Width());
    )

    if( to == from )
        return;
    const Int m = A.Height();
    T* ABuf = A.Buffer();
    const Int ALDim = A.LDim();
    blas::Swap( m, &ABuf[to*ALDim], 1, &ABuf[from*ALDim], 1 );
}

template<typename T>
void ColSwap( AbstractDistMatrix<T>& A, Int to, Int from )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( to < 0 || to >= A.Width() || from < 0 || from >= A.Width() )
          LogicError
          ("Attempted invalid col swap, (",to,",",from,") of matrix of width ",
           A.Width());
    )
    if( to == from )
        return;
    if( !A.Participating() )
        return;
    const Int mLocal = A.LocalHeight();
    const Int rowAlign = A.RowAlign();
    const Int rowShift = A.RowShift();
    const Int rowStride = A.RowStride();
    const Int toMod = Mod(to,rowStride);
    const Int fromMod = Mod(from,rowStride);
    T* ABuf = A.Buffer();
    const Int ALDim = A.LDim();

    if( toMod == fromMod )
    {
        const Int jLocTo = (to-rowShift) / rowStride;
        const Int jLocFrom = (from-rowShift) / rowStride;
        if( toMod == rowShift )
            blas::Swap
            ( mLocal, &ABuf[jLocTo*ALDim], 1, &ABuf[jLocFrom*ALDim], 1 );
    }
    else if( toMod == rowShift )
    {
        const Int jLocTo = (to-rowShift) / rowStride;
        const int fromOwner = Mod(from+rowAlign,rowStride);
        mpi::SendRecv
        ( &ABuf[jLocTo*ALDim], mLocal, fromOwner, fromOwner, A.RowComm() );
    }
    else if( fromMod == rowShift )
    {
        const Int jLocFrom = (from-rowShift) / rowStride;
        const int toOwner = Mod(to+rowAlign,rowStride);
        mpi::SendRecv
        ( &ABuf[jLocFrom*ALDim], mLocal, toOwner, toOwner, A.RowComm() );
    }
}

template<typename T>
void SymmetricSwap
( UpperOrLower uplo, Matrix<T>& A, Int to, Int from, bool conjugate )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( A.Height() != A.Width() )
          LogicError("A must be square");
      if( to < 0 || to >= A.Height() || from < 0 || from >= A.Height() )
          LogicError
          ("Attempted invalid symmetric swap, (",to,",",from,
           ") of matrix of size ",A.Height());
    )

    if( to == from )
    {
        if( conjugate )
            A.MakeReal( to, to );
        return;
    }
    const Int n = A.Height();
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );
    if( to > from )
        std::swap( to, from ); 
    if( uplo == LOWER )
    { 
        // Bottom swap
        if( from+1 < n )
        {
            auto ABot = A( IR(from+1,n), IR(0,n) );
            ColSwap( ABot, to, from );
        }
        // Inner swap
        if( to+1 < from )
        {
            auto aToInner = A( IR(to+1,from), IR(to,to+1) );
            auto aFromInner = A( IR(from,from+1), IR(to+1,from) );
            Swap( orientation, aToInner, aFromInner );
        }
        // Corner swap
        if( conjugate )
            A.Conjugate( from, to );
        // Diagonal swap
        {
            const T value = A(from,from);
            A(from,from) = A(to,to);
            A(to,to) = value;
            if( conjugate )
            {
                A.MakeReal( to, to );
                A.MakeReal( from, from );
            }
        }
        // Left swap
        if( to > 0 )
        {
            auto ALeft = A( IR(0,n), IR(0,to) );
            RowSwap( ALeft, to, from ); 
        }
    }
    else
    {
        // Right swap
        if( from+1 < n )
        {
            auto ARight = A( IR(0,n), IR(from+1,n) );
            RowSwap( ARight, to, from );
        }
        // Inner swap
        if( to+1 < from )
        {
            auto aToInner = A( IR(to,to+1), IR(to+1,from) );
            auto aFromInner = A( IR(to+1,from), IR(from,from+1) );
            Swap( orientation, aToInner, aFromInner );
        }
        // Corner swap
        if( conjugate )
            A.Conjugate( to, from );
        // Diagonal swap
        {
            const T value = A(from,from);
            A(from,from) = A(to,to);
            A(to,to) = value;
            if( conjugate )
            {
                A.MakeReal( to, to );
                A.MakeReal( from, from );
            }
        }
        // Top swap
        if( to > 0 )
        {
            auto ATop = A( IR(0,to), IR(0,n) );
            ColSwap( ATop, to, from ); 
        }
    }
}

template<typename T>
void SymmetricSwap
( UpperOrLower uplo, AbstractDistMatrix<T>& A, 
  Int to, Int from, bool conjugate )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( A.Height() != A.Width() )
          LogicError("A must be square");
      if( to < 0 || to >= A.Height() || from < 0 || from >= A.Height() )
          LogicError
          ("Attempted invalid symmetric swap, (",to,",",from,
           ") of matrix of size ",A.Height());
    )
    typedef unique_ptr<AbstractDistMatrix<T>> ADMPtr;

    if( to == from )
    {
        if( conjugate )
            A.MakeReal( to, to );
        return;
    }
    const Int n = A.Height();
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );
    if( to > from )
        std::swap( to, from );
    if( uplo == LOWER )
    { 
        // Bottom swap
        if( from+1 < n )
        {
            ADMPtr ABot( A.Construct( A.Grid(), A.Root() ) );
            View( *ABot, A, IR(from+1,n), IR(0,n) );
            ColSwap( *ABot, to, from );
        }
        // Inner swap
        if( to+1 < from )
        {
            ADMPtr aToInner( A.Construct( A.Grid(), A.Root() ) );
            ADMPtr aFromInner( A.Construct( A.Grid(), A.Root() ) );
            View( *aToInner, A, IR(to+1,from), IR(to,to+1) );
            View( *aFromInner, A, IR(from,from+1), IR(to+1,from) );
            Swap( orientation, *aToInner, *aFromInner );
        }
        // Corner swap
        if( conjugate )
            A.Conjugate( from, to );
        // Diagonal swap
        {
            const T value = A.Get(from,from);
            A.Set( from, from, A.Get(to,to) );
            A.Set( to,   to,   value        );
            if( conjugate )
            {
                A.MakeReal( to, to );
                A.MakeReal( from, from );
            }
        }
        // Left swap
        if( to > 0 )
        {
            ADMPtr ALeft( A.Construct( A.Grid(), A.Root() ) );
            View( *ALeft, A, IR(0,n), IR(0,to) );
            RowSwap( *ALeft, to, from ); 
        }
    }
    else
    {
        // Right swap
        if( from+1 < n )
        {
            ADMPtr ARight( A.Construct( A.Grid(), A.Root() ) );
            View( *ARight, A, IR(0,n), IR(from+1,n) );
            RowSwap( *ARight, to, from );
        }
        // Inner swap
        if( to+1 < from )
        {
            ADMPtr aToInner( A.Construct( A.Grid(), A.Root() ) );
            ADMPtr aFromInner( A.Construct( A.Grid(), A.Root() ) );
            View( *aToInner, A, IR(to,to+1), IR(to+1,from) );
            View( *aFromInner, A, IR(to+1,from), IR(from,from+1) );
            Swap( orientation, *aToInner, *aFromInner );
        }
        // Corner swap
        if( conjugate )
            A.Conjugate( to, from );
        // Diagonal swap
        {
            const T value = A.Get(from,from);
            A.Set( from, from, A.Get(to,to) );
            A.Set( to,   to,   value        );
            if( conjugate )
            {
                A.MakeReal( to, to );
                A.MakeReal( from, from );
            }
        }
        // Top swap
        if( to > 0 )
        {
            ADMPtr ATop( A.Construct( A.Grid(), A.Root() ) );
            View( *ATop, A, IR(0,to), IR(0,n) );
            ColSwap( *ATop, to, from ); 
        }
    }
}

template<typename T>
void HermitianSwap( UpperOrLower uplo, Matrix<T>& A, Int to, Int from )
{
    DEBUG_CSE
    SymmetricSwap( uplo, A, to, from, true );
}

template<typename T>
void HermitianSwap
( UpperOrLower uplo, AbstractDistMatrix<T>& A, Int to, Int from )
{
    DEBUG_CSE
    SymmetricSwap( uplo, A, to, from, true );
}

#define PROTO(T) \
  template void Swap( Orientation orientation, Matrix<T>& X, Matrix<T>& Y ); \
  template void Swap \
  ( Orientation orientation, \
    AbstractDistMatrix<T>& X, AbstractDistMatrix<T>& Y ); \
  template void RowSwap( Matrix<T>& A, Int to, Int from ); \
  template void RowSwap( AbstractDistMatrix<T>& A, Int to, Int from ); \
  template void ColSwap( Matrix<T>& A, Int to, Int from ); \
  template void ColSwap( AbstractDistMatrix<T>& A, Int to, Int from ); \
  template void SymmetricSwap \
  ( UpperOrLower uplo, Matrix<T>& A, Int to, Int from, bool conjugate ); \
  template void SymmetricSwap \
  ( UpperOrLower uplo, AbstractDistMatrix<T>& A, Int to, Int from, \
    bool conjugate ); \
  template void HermitianSwap \
  ( UpperOrLower uplo, Matrix<T>& A, Int to, Int from ); \
  template void HermitianSwap \
  ( UpperOrLower uplo, AbstractDistMatrix<T>& A, Int to, Int from );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
