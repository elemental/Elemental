/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_COPY_HPP
#define EL_COPY_HPP

namespace El {

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

template<typename T>
inline void
Copy( const DynamicDistMatrix<T>& ADyn, DynamicDistMatrix<T>& BDyn )
{
    DEBUG_ONLY(CallStackEntry cse("Copy"))
    #define INNER_IF_CONVERT_AND_COPY(A,BDYN,CDIST,RDIST) \
      if( BDYN.U == CDIST && BDYN.V == RDIST ) \
      { \
          auto B = dynamic_cast<DistMatrix<T,CDIST,RDIST>*>(BDYN.ADM); \
          if( B == nullptr ) \
              RuntimeError("Dynamic cast failed"); \
          *B = *A; \
      }
    #define INNER_ELSEIF_CONVERT_AND_COPY(A,BDYN,CDIST,RDIST) \
      else INNER_IF_CONVERT_AND_COPY(A,BDYN,CDIST,RDIST)
    #define IF_CONVERT_AND_COPY(ADYN,BDYN,CDIST,RDIST) \
      if( ADYN.U == CDIST && ADYN.V == RDIST ) \
      { \
          auto A = dynamic_cast<const DistMatrix<T,CDIST,RDIST>*>(ADYN.ADM); \
          if( A == nullptr ) \
              RuntimeError("Dynamic cast failed"); \
          INNER_IF_CONVERT_AND_COPY(    A,BDYN,CIRC,CIRC) \
          INNER_ELSEIF_CONVERT_AND_COPY(A,BDYN,MC,  MR  ) \
          INNER_ELSEIF_CONVERT_AND_COPY(A,BDYN,MD,  STAR) \
          INNER_ELSEIF_CONVERT_AND_COPY(A,BDYN,MR,  MC  ) \
          INNER_ELSEIF_CONVERT_AND_COPY(A,BDYN,MR,  STAR) \
          INNER_ELSEIF_CONVERT_AND_COPY(A,BDYN,STAR,MC  ) \
          INNER_ELSEIF_CONVERT_AND_COPY(A,BDYN,STAR,MD  ) \
          INNER_ELSEIF_CONVERT_AND_COPY(A,BDYN,STAR,MR  ) \
          INNER_ELSEIF_CONVERT_AND_COPY(A,BDYN,STAR,STAR) \
          INNER_ELSEIF_CONVERT_AND_COPY(A,BDYN,STAR,VC  ) \
          INNER_ELSEIF_CONVERT_AND_COPY(A,BDYN,STAR,VR  ) \
          INNER_ELSEIF_CONVERT_AND_COPY(A,BDYN,VC,  STAR) \
          INNER_ELSEIF_CONVERT_AND_COPY(A,BDYN,VR,  STAR) \
      }
    #define ELSEIF_CONVERT_AND_COPY(ADYN,BDYN,CDIST,RDIST) \
      else IF_CONVERT_AND_COPY(ADYN,BDYN,CDIST,RDIST)

    IF_CONVERT_AND_COPY(    ADyn,BDyn,CIRC,CIRC)
    ELSEIF_CONVERT_AND_COPY(ADyn,BDyn,MC,  MR  )
    ELSEIF_CONVERT_AND_COPY(ADyn,BDyn,MC,  STAR)
    ELSEIF_CONVERT_AND_COPY(ADyn,BDyn,MD,  STAR)
    ELSEIF_CONVERT_AND_COPY(ADyn,BDyn,MR,  MC  )
    ELSEIF_CONVERT_AND_COPY(ADyn,BDyn,MR,  STAR)
    ELSEIF_CONVERT_AND_COPY(ADyn,BDyn,STAR,MD  )
    ELSEIF_CONVERT_AND_COPY(ADyn,BDyn,STAR,MR  )
    ELSEIF_CONVERT_AND_COPY(ADyn,BDyn,STAR,STAR)
    ELSEIF_CONVERT_AND_COPY(ADyn,BDyn,STAR,VC  )
    ELSEIF_CONVERT_AND_COPY(ADyn,BDyn,STAR,VR  )
    ELSEIF_CONVERT_AND_COPY(ADyn,BDyn,VC,  STAR)
    ELSEIF_CONVERT_AND_COPY(ADyn,BDyn,VR,  STAR)

    #undef ELSEIF_CONVERT_AND_COPY
    #undef IF_CONVERT_AND_COPY
    #undef INNER_ELSEIF_CONVERT_AND_COPY
    #undef INNER_IF_CONVERT_AND_COPY
}

} // namespace El

#endif // ifndef EL_COPY_HPP
