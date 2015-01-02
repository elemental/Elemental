/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename T>
void Swap( Orientation orientation, Matrix<T>& X, Matrix<T>& Y )
{
    DEBUG_ONLY(CallStackEntry cse("Swap"))
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
                    const T alpha = X.Get(i,j);
                    X.Set( i, j, Conj(Y.Get(j,i)) );
                    Y.Set( j, i, Conj(alpha)      );
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
( Orientation orientation, AbstractDistMatrix<T>& X, AbstractDistMatrix<T>& Y )
{
    DEBUG_ONLY(CallStackEntry cse("Swap"))
    const Grid& g = X.Grid();
    if( orientation == NORMAL )
    {
        DEBUG_ONLY(
            if( Y.Height() != X.Height() || Y.Width() != X.Width() )
                LogicError("Invalid submatrix sizes");
        )
        // TODO: Optimize communication

        std::unique_ptr<AbstractDistMatrix<T>> 
          YLikeX( X.Construct(g,X.Root()) );
        YLikeX->AlignWith( X.DistData() );
        Copy( Y, *YLikeX );

        std::unique_ptr<AbstractDistMatrix<T>> 
          XLikeY( Y.Construct(g,Y.Root()) );
        XLikeY->AlignWith( Y.DistData() );
        Copy( X, *XLikeY );

        Copy( XLikeY->Matrix(), Y.Matrix() );
        Copy( YLikeX->Matrix(), X.Matrix() );
    }
    else
    {
        const bool conjugate = ( orientation==ADJOINT );
        DEBUG_ONLY(
            if( Y.Width() != X.Height() || Y.Height() != X.Width() )
                LogicError("Invalid submatrix sizes");
        )

        // TODO: Optimize communication
        std::unique_ptr<AbstractDistMatrix<T>> 
          YTransLikeX( X.Construct(g,X.Root()) );
        YTransLikeX->AlignWith( X.DistData() );
        Transpose( Y, *YTransLikeX, conjugate );

        std::unique_ptr<AbstractDistMatrix<T>> 
          XTransLikeY( Y.Construct(g,Y.Root()) );
        XTransLikeY->AlignWith( Y.DistData() );
        Transpose( X, *XTransLikeY, conjugate );

        Copy( XTransLikeY->Matrix(), Y.Matrix() );
        Copy( YTransLikeX->Matrix(), X.Matrix() );
    }
}

template<typename T>
void RowSwap( Matrix<T>& A, Int to, Int from )
{
    DEBUG_ONLY(CallStackEntry cse("RowSwap"))
    if( to == from )
        return;
    const Int n = A.Width();
    auto aToRow = A( IR(to,to+1), IR(0,n) );
    auto aFromRow = A( IR(from,from+1), IR(0,n) );
    Swap( NORMAL, aToRow, aFromRow );
}

template<typename T>
void RowSwap( AbstractDistMatrix<T>& A, Int to, Int from )
{
    DEBUG_ONLY(CallStackEntry cse("RowSwap"))
    if( to == from )
        return;
    if( !A.Participating() )
        return;
    const Int n = A.Width();
    const Int nLocal = A.LocalWidth();
    std::unique_ptr<AbstractDistMatrix<T>> 
      aToRow( A.Construct(A.Grid(),A.Root()) );
    std::unique_ptr<AbstractDistMatrix<T>> 
      aFromRow( A.Construct(A.Grid(),A.Root()) );
    View( *aToRow, A, IR(to,to+1), IR(0,n) );
    View( *aFromRow, A, IR(from,from+1), IR(0,n) );
    if( aToRow->ColAlign() == aFromRow->ColAlign() )
    {
        if( aToRow->ColShift() == 0 )
            Swap( NORMAL, aToRow->Matrix(), aFromRow->Matrix() );
    }
    else if( aToRow->ColShift() == 0 )
    {
        std::vector<T> buf( nLocal );
        for( Int jLoc=0; jLoc<nLocal; ++jLoc )
            buf[jLoc] = aToRow->GetLocal(0,jLoc);
        mpi::SendRecv
        ( buf.data(), nLocal, 
          aFromRow->ColAlign(), aFromRow->ColAlign(), A.ColComm() );
        for( Int jLoc=0; jLoc<nLocal; ++jLoc )
            aToRow->SetLocal(0,jLoc,buf[jLoc]);
    }
    else if( aFromRow->ColShift() == 0 )
    {
        std::vector<T> buf( nLocal );
        for( Int jLoc=0; jLoc<nLocal; ++jLoc )
            buf[jLoc] = aFromRow->GetLocal(0,jLoc);
        mpi::SendRecv
        ( buf.data(), nLocal, 
          aToRow->ColAlign(), aToRow->ColAlign(), A.ColComm() );
        for( Int jLoc=0; jLoc<nLocal; ++jLoc )
            aFromRow->SetLocal(0,jLoc,buf[jLoc]);
    }
}

template<typename T>
void ColSwap( Matrix<T>& A, Int to, Int from )
{
    DEBUG_ONLY(CallStackEntry cse("ColSwap"))
    if( to == from )
        return;
    const Int m = A.Height();
    auto aToCol = A( IR(0,m), IR(to,to+1) );
    auto aFromCol = A( IR(0,m), IR(from,from+1) );
    Swap( NORMAL, aToCol, aFromCol );
}

template<typename T>
void ColSwap( AbstractDistMatrix<T>& A, Int to, Int from )
{
    DEBUG_ONLY(CallStackEntry cse("ColSwap"))
    if( to == from )
        return;
    if( !A.Participating() )
        return;
    const Int m = A.Height();
    const Int mLocal = A.LocalHeight();
    std::unique_ptr<AbstractDistMatrix<T>> 
      aToCol( A.Construct(A.Grid(),A.Root()) );
    std::unique_ptr<AbstractDistMatrix<T>> 
      aFromCol( A.Construct(A.Grid(),A.Root()) );
    View( *aToCol, A, IR(0,m), IR(to,to+1) );
    View( *aFromCol, A, IR(0,m), IR(from,from+1) );
    if( aToCol->RowAlign() == aFromCol->RowAlign() )
    {
        if( aToCol->RowShift() == 0 )
            Swap( NORMAL, aToCol->Matrix(), aFromCol->Matrix() );
    }
    else if( aToCol->RowShift() == 0 )
    {
        mpi::SendRecv
        ( aToCol->Buffer(), mLocal, 
          aFromCol->RowAlign(), aFromCol->RowAlign(), A.RowComm() );
    }
    else if( aFromCol->RowShift() == 0 )
    {
        mpi::SendRecv
        ( aFromCol->Buffer(), mLocal, 
          aToCol->RowAlign(), aToCol->RowAlign(), A.RowComm() );
    }
}

template<typename T>
void SymmetricSwap
( UpperOrLower uplo, Matrix<T>& A, Int to, Int from, bool conjugate )
{
    DEBUG_ONLY(
        CallStackEntry cse("SymmetricSwap");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
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
        std::swap( to, from ); if( uplo == LOWER )
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
    DEBUG_ONLY(
        CallStackEntry cse("SymmetricSwap");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
    )
    typedef std::unique_ptr<AbstractDistMatrix<T>> ADMPtr;

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
    DEBUG_ONLY(CallStackEntry cse("HermitianSwap"))
    SymmetricSwap( uplo, A, to, from, true );
}

template<typename T>
void HermitianSwap
( UpperOrLower uplo, AbstractDistMatrix<T>& A, Int to, Int from )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianSwap"))
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

#include "El/macros/Instantiate.h"

} // namespace El
