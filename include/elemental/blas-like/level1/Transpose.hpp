/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_TRANSPOSE_HPP
#define ELEM_TRANSPOSE_HPP

namespace elem {

template<typename T>
inline void
Transpose( const Matrix<T>& A, Matrix<T>& B, bool conjugate=false )
{
    DEBUG_ONLY(CallStackEntry cse("Transpose"))
    const Int m = A.Height();
    const Int n = A.Width();
    B.Resize( n, m );
    if( conjugate )
    {
        for( Int j=0; j<n; ++j )
            for( Int i=0; i<m; ++i )
                B.Set(j,i,Conj(A.Get(i,j)));
    }
    else
    {
        for( Int j=0; j<n; ++j )
            for( Int i=0; i<m; ++i )
                B.Set(j,i,A.Get(i,j));
    }
}

template<typename T,Dist U,Dist V,
                    Dist W,Dist Z>
inline void
Transpose
( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B, bool conjugate=false )
{
    DEBUG_ONLY(CallStackEntry cse("Transpose"))
    if( U == Z && V == W && 
        (A.ColAlign()==B.RowAlign() || !B.RowConstrained()) &&
        (A.RowAlign()==B.ColAlign() || !B.ColConstrained()) )
    {
        B.Align( A.RowAlign(), A.ColAlign() );
        B.Resize( A.Width(), A.Height() );
        Transpose( A.LockedMatrix(), B.Matrix(), conjugate );
    }
    else
    {
        DistMatrix<T,Z,W> C( B.Grid() );
        C.AlignRowsWith( B );
        C.AlignColsWith( B );
        C = A;
        B.Resize( A.Width(), A.Height() );
        Transpose( C.LockedMatrix(), B.Matrix(), conjugate );
    }
}

template<typename T,Dist U,Dist V,
                    Dist W,Dist Z>
inline void
Transpose
( const BlockDistMatrix<T,U,V>& A, BlockDistMatrix<T,W,Z>& B, 
  bool conjugate=false )
{
    DEBUG_ONLY(CallStackEntry cse("Transpose"))
    // TODO: Member function to quickly query this information?!?
    //       Perhaps A.ColsAlignedWithRows( B )?
    //               A.ColsAlignedWithCols( B )?
    const bool transposeColAligned = 
        U               == Z              &&
        A.BlockHeight() == B.BlockWidth() &&
        A.ColAlign()    == B.RowAlign()   &&
        A.ColCut()      == B.RowCut();
    const bool transposeRowAligned =
        V              == W               &&
        A.BlockWidth() == B.BlockHeight() &&
        A.RowAlign()   == B.ColAlign()    &&
        A.RowCut()     == B.ColCut();
    if( (transposeColAligned || !B.RowConstrained()) &&
        (transposeRowAligned || !B.ColConstrained()) )
    {
        B.Align
        ( A.BlockWidth(), A.BlockHeight(), 
          A.RowAlign(), A.ColAlign(), A.RowCut(), A.ColCut() );
        B.Resize( A.Width(), A.Height() );
        Transpose( A.LockedMatrix(), B.Matrix(), conjugate );
    }
    else
    {
        BlockDistMatrix<T,Z,W> C( B.Grid() );
        C.AlignRowsWith( B );
        C.AlignColsWith( B );
        C = A;
        B.Resize( A.Width(), A.Height() );
        Transpose( C.LockedMatrix(), B.Matrix(), conjugate );
    }
}

} // namespace elem

#endif // ifndef ELEM_TRANSPOSE_HPP
