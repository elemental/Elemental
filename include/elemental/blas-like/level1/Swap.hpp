/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_BLAS_SWAP_HPP
#define ELEM_BLAS_SWAP_HPP

namespace elem {

template<typename F>
inline void Swap( Orientation orientation, Matrix<F>& X, Matrix<F>& Y )
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
                    const F alpha = X.Get(i,j);
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

template<typename F,Dist U1,Dist V1,
                    Dist U2,Dist V2>
inline void Swap
( Orientation orientation, DistMatrix<F,U1,V1>& X, DistMatrix<F,U2,V2>& Y )
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

        DistMatrix<F,U1,V1> YLikeX(g);
        YLikeX.AlignWith( X );
        YLikeX = Y;

        DistMatrix<F,U2,V2> XLikeY(g);
        XLikeY.AlignWith( Y );
        XLikeY = X; 

        Y = XLikeY;
        X = YLikeX; 
    }
    else
    {
        const bool conjugate = ( orientation==ADJOINT );
        DEBUG_ONLY(
            if( Y.Width() != X.Height() || Y.Height() != X.Width() )
                LogicError("Invalid submatrix sizes");
        )

        // TODO: Optimize communication

        DistMatrix<F,U1,V1> YTransLikeX(g);
        YTransLikeX.AlignWith( X );
        Transpose( Y, YTransLikeX, conjugate );

        DistMatrix<F,U2,V2> XTransLikeY(g);
        XTransLikeY.AlignWith( Y );
        Transpose( X, XTransLikeY, conjugate );

        Y = XTransLikeY;
        X = YTransLikeX;
    }
}

template<typename F>
inline void RowSwap( Matrix<F>& A, Int to, Int from )
{
    DEBUG_ONLY(CallStackEntry cse("RowSwap"))
    if( to == from )
        return;
    const Int n = A.Width();
    auto aToRow   = ViewRange( A, to,   0, to+1,   n );
    auto aFromRow = ViewRange( A, from, 0, from+1, n );
    Swap( NORMAL, aToRow, aFromRow );
}

template<typename F,Dist U,Dist V>
inline void RowSwap( DistMatrix<F,U,V>& A, Int to, Int from )
{
    DEBUG_ONLY(CallStackEntry cse("RowSwap"))
    if( to == from )
        return;
    if( !A.Participating() )
        return;
    const Int n = A.Width();
    const Int nLocal = A.LocalWidth();
    auto aToRow   = ViewRange( A, to,   0, to+1,   n );
    auto aFromRow = ViewRange( A, from, 0, from+1, n );
    if( aToRow.ColAlign() == aFromRow.ColAlign() )
    {
        if( aToRow.ColShift() == 0 )
            Swap( NORMAL, aToRow.Matrix(), aFromRow.Matrix() );
    }
    else if( aToRow.ColShift() == 0 )
    {
        std::vector<F> buf( nLocal );
        for( Int jLoc=0; jLoc<nLocal; ++jLoc )
            buf[jLoc] = aToRow.GetLocal(0,jLoc);
        mpi::SendRecv
        ( buf.data(), nLocal, 
          aFromRow.ColAlign(), aFromRow.ColAlign(), A.ColComm() );
        for( Int jLoc=0; jLoc<nLocal; ++jLoc )
            aToRow.SetLocal(0,jLoc,buf[jLoc]);
    }
    else if( aFromRow.ColShift() == 0 )
    {
        std::vector<F> buf( nLocal );
        for( Int jLoc=0; jLoc<nLocal; ++jLoc )
            buf[jLoc] = aFromRow.GetLocal(0,jLoc);
        mpi::SendRecv
        ( buf.data(), nLocal, 
          aToRow.ColAlign(), aToRow.ColAlign(), A.ColComm() );
        for( Int jLoc=0; jLoc<nLocal; ++jLoc )
            aFromRow.SetLocal(0,jLoc,buf[jLoc]);
    }
}

template<typename F>
inline void ColSwap( Matrix<F>& A, Int to, Int from )
{
    DEBUG_ONLY(CallStackEntry cse("ColSwap"))
    if( to == from )
        return;
    const Int m = A.Height();
    auto aToCol   = ViewRange( A, 0, to,   m, to+1   );
    auto aFromCol = ViewRange( A, 0, from, m, from+1 );
    Swap( NORMAL, aToCol, aFromCol );
}

template<typename F,Dist U,Dist V>
inline void ColSwap( DistMatrix<F,U,V>& A, Int to, Int from )
{
    DEBUG_ONLY(CallStackEntry cse("ColSwap"))
    if( to == from )
        return;
    if( !A.Participating() )
        return;
    const Int m = A.Height();
    const Int mLocal = A.LocalHeight();
    auto aToCol   = ViewRange( A, 0, to,   m, to+1   );
    auto aFromCol = ViewRange( A, 0, from, m, from+1 );
    if( aToCol.RowAlign() == aFromCol.RowAlign() )
    {
        if( aToCol.RowShift() == 0 )
            Swap( NORMAL, aToCol.Matrix(), aFromCol.Matrix() );
    }
    else if( aToCol.RowShift() == 0 )
    {
        mpi::SendRecv
        ( aToCol.Buffer(), mLocal, 
          aFromCol.RowAlign(), aFromCol.RowAlign(), A.RowComm() );
    }
    else if( aFromCol.RowShift() == 0 )
    {
        mpi::SendRecv
        ( aFromCol.Buffer(), mLocal, 
          aToCol.RowAlign(), aToCol.RowAlign(), A.RowComm() );
    }
}

template<typename F>
inline void SymmetricSwap
( UpperOrLower uplo, Matrix<F>& A, int to, int from, bool conjugate=false )
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
        std::swap( to, from );
    if( uplo == LOWER )
    { 
        // Bottom swap
        if( from+1 < n )
        {
            auto ABot = ViewRange( A, from+1, 0, n, n );
            ColSwap( ABot, to, from );
        }
        // Inner swap
        if( to+1 < from )
        {
            auto aToInner   = ViewRange( A, to+1, to,   from,   to+1 );
            auto aFromInner = ViewRange( A, from, to+1, from+1, from );
            Swap( orientation, aToInner, aFromInner );
        }
        // Corner swap
        if( conjugate )
            A.Conjugate( from, to );
        // Diagonal swap
        {
            const F value = A.Get(from,from);
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
            auto ALeft = ViewRange( A, 0, 0, n, to );
            RowSwap( ALeft, to, from ); 
        }
    }
    else
    {
        // Right swap
        if( from+1 < n )
        {
            auto ARight = ViewRange( A, 0, from+1, n, n );
            RowSwap( ARight, to, from );
        }
        // Inner swap
        if( to+1 < from )
        {
            auto aToInner   = ViewRange( A, to,   to+1, to+1, from   );
            auto aFromInner = ViewRange( A, to+1, from, from, from+1 );
            Swap( orientation, aToInner, aFromInner );
        }
        // Corner swap
        if( conjugate )
            A.Conjugate( to, from );
        // Diagonal swap
        {
            const F value = A.Get(from,from);
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
            auto ATop = ViewRange( A, 0, 0, to, n );
            ColSwap( ATop, to, from ); 
        }
    }
}

template<typename F,Dist U,Dist V>
inline void SymmetricSwap
( UpperOrLower uplo, DistMatrix<F,U,V>& A, 
  int to, int from, bool conjugate=false )
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
        std::swap( to, from );
    if( uplo == LOWER )
    { 
        // Bottom swap
        if( from+1 < n )
        {
            auto ABot = ViewRange( A, from+1, 0, n, n );
            ColSwap( ABot, to, from );
        }
        // Inner swap
        if( to+1 < from )
        {
            auto aToInner   = ViewRange( A, to+1, to,   from,   to+1 );
            auto aFromInner = ViewRange( A, from, to+1, from+1, from );
            Swap( orientation, aToInner, aFromInner );
        }
        // Corner swap
        if( conjugate )
            A.Conjugate( from, to );
        // Diagonal swap
        {
            const F value = A.Get(from,from);
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
            auto ALeft = ViewRange( A, 0, 0, n, to );
            RowSwap( ALeft, to, from ); 
        }
    }
    else
    {
        // Right swap
        if( from+1 < n )
        {
            auto ARight = ViewRange( A, 0, from+1, n, n );
            RowSwap( ARight, to, from );
        }
        // Inner swap
        if( to+1 < from )
        {
            auto aToInner   = ViewRange( A, to,   to+1, to+1, from   );
            auto aFromInner = ViewRange( A, to+1, from, from, from+1 );
            Swap( orientation, aToInner, aFromInner );
        }
        // Corner swap
        if( conjugate )
            A.Conjugate( to, from );
        // Diagonal swap
        {
            const F value = A.Get(from,from);
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
            auto ATop = ViewRange( A, 0, 0, to, n );
            ColSwap( ATop, to, from ); 
        }
    }
}

template<typename F>
inline void HermitianSwap( UpperOrLower uplo, Matrix<F>& A, int to, int from )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianSwap"))
    SymmetricSwap( uplo, A, to, from, true );
}

template<typename F,Dist U,Dist V>
inline void HermitianSwap
( UpperOrLower uplo, DistMatrix<F,U,V>& A, int to, int from )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianSwap"))
    SymmetricSwap( uplo, A, to, from, true );
}

} // namespace elem

#endif // ifndef ELEM_BLAS_SWAP_HPP
