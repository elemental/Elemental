/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_COPY_HPP
#define ELEM_COPY_HPP

namespace elem {

template<typename T>
inline void
Copy( const Matrix<T>& A, Matrix<T>& B )
{
    DEBUG_ONLY(CallStackEntry cse("Copy"))
    B = A;
}

template<typename Real>
inline void
Copy( const Matrix<Real>& A, Matrix<Complex<Real>>& B )
{
    DEBUG_ONLY(CallStackEntry cse("Copy"))
    const Int m = A.Height();
    const Int n = A.Width();
    B.Resize( m, n );
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            B.Set( i, j, A.Get(i,j) );
}

template<typename T,Dist U,Dist V,Dist W,Dist Z>
inline void
Copy( const DistMatrix<T,U,V>& A, DistMatrix<T,W,Z>& B )
{
    DEBUG_ONLY(CallStackEntry cse("Copy"))
    B = A;
}

template<typename Real,Dist U,Dist V,Dist W,Dist Z>
inline void
Copy( const DistMatrix<Real,U,V>& A, DistMatrix<Complex<Real>,W,Z>& B )
{
    DEBUG_ONLY(CallStackEntry cse("Copy"))

    if( U == W && V == Z )
    {
        if( !B.ColConstrained() )
            B.AlignCols( A.ColAlign() );
        if( !B.RowConstrained() )
            B.AlignRows( A.RowAlign() );
        if( A.ColAlign() == B.ColAlign() && A.RowAlign() == B.RowAlign() )
        {
            B.Resize( A.Height(), A.Width() );
            Copy( A.LockedMatrix(), B.Matrix() );
            return;
        }
    }

    DistMatrix<Real,W,Z> BReal(A.Grid());
    BReal.AlignWith( B );
    BReal = A;
    B.Resize( A.Height(), A.Width() );
    Copy( BReal.LockedMatrix(), B.Matrix() );
}

template<typename T,Dist U,Dist V,Dist W,Dist Z>
inline void
Copy( const BlockDistMatrix<T,U,V>& A, BlockDistMatrix<T,W,Z>& B )
{
    DEBUG_ONLY(CallStackEntry cse("Copy"))
    B = A;
}

template<typename Real,Dist U,Dist V,Dist W,Dist Z>
inline void
Copy
( const BlockDistMatrix<Real,U,V>& A, BlockDistMatrix<Complex<Real>,W,Z>& B )
{
    DEBUG_ONLY(CallStackEntry cse("Copy"))

    if( U == W && V == Z )
    {
        if( !B.ColConstrained() )
            B.AlignCols( A.ColAlign() );
        if( !B.RowConstrained() )
            B.AlignRows( A.RowAlign() );
        if( A.ColAlign() == B.ColAlign() && 
            A.RowAlign() == B.RowAlign() && 
            A.ColCut() == B.ColCut() &&
            A.RowCut() == B.RowCut() )
        {
            B.Resize( A.Height(), A.Width() );
            Copy( A.LockedMatrix(), B.Matrix() );
            return;
        }
    }

    BlockDistMatrix<Real,W,Z> BReal(A.Grid());
    BReal.AlignWith( B );
    BReal = A;
    B.Resize( A.Height(), A.Width() );
    Copy( BReal.LockedMatrix(), B.Matrix() );
}

} // namespace elem

#endif // ifndef ELEM_COPY_HPP
