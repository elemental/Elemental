/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_DIAGONALSCALETRAPEZOID_HPP
#define EL_BLAS_DIAGONALSCALETRAPEZOID_HPP

namespace El {

template<typename TDiag,typename T>
void DiagonalScaleTrapezoid
( LeftOrRight side,
  UpperOrLower uplo,
  Orientation orientation, 
  const Matrix<TDiag>& d,
        Matrix<T>& A,
  Int offset )
{
    DEBUG_ONLY(
      CSE cse("DiagonalScaleTrapezoid");
      if( side==LEFT && (d.Height()!=A.Height() || d.Width()!=1) )
          LogicError("d should have been a vector of the height of A");
      if( side==RIGHT && (d.Height()!=A.Width() || d.Width()!=1) )
          LogicError("d should have been a vector of the width of A");
    )
    const Int m = A.Height();
    const Int n = A.Width();
    const Int diagLength = A.DiagonalLength(offset);
    const Int ldim = A.LDim();
    T* ABuf = A.Buffer();
    const bool conjugate = ( orientation==ADJOINT );

    const Int iOff = ( offset>=0 ? 0      : -offset );
    const Int jOff = ( offset>=0 ? offset : 0       );

    if( uplo == LOWER && side == LEFT )
    {
        // Scale from the left up to the diagonal
        for( Int i=iOff; i<m; ++i )
        {
            const Int k = i-iOff;
            const Int j = k+jOff;
            const TDiag alpha = ( conjugate ? Conj(d.Get(i,0)) : d.Get(i,0) );
            blas::Scal( Min(j+1,n), alpha, &ABuf[i], ldim );
        }
    }
    else if( uplo == UPPER && side == LEFT )
    {
        // Scale from the diagonal to the right
        for( Int i=0; i<iOff+diagLength; ++i )
        {
            const Int k = i-iOff;
            const Int j = k+jOff;
            const Int jLeft = Max(j,0);
            const TDiag alpha = ( conjugate ? Conj(d.Get(i,0)) : d.Get(i,0) );
            blas::Scal( n-jLeft, alpha, &ABuf[i+jLeft*ldim], ldim );
        }
    }
    else if( uplo == LOWER && side == RIGHT )
    {
        // Scale from the diagonal downwards
        for( Int j=0; j<jOff+diagLength; ++j )
        {
            const Int k = j-jOff;
            const Int i = k+iOff;
            const Int iTop = Max(i,0);
            const TDiag alpha = ( conjugate ? Conj(d.Get(j,0)) : d.Get(j,0) );
            blas::Scal( m-iTop, alpha, &ABuf[iTop+j*ldim], 1 );
        }
    }
    else /* uplo == UPPER && side == RIGHT */
    {
        // Scale downward to the diagonal
        for( Int j=jOff; j<n; ++j )
        {
            const Int k = j-jOff;
            const Int i = k+iOff;
            const TDiag alpha = ( conjugate ? Conj(d.Get(j,0)) : d.Get(j,0) );
            blas::Scal( Min(i+1,m), alpha, &ABuf[j*ldim], 1 );
        }
    }
}

template<typename TDiag,typename T,Dist U,Dist V>
void DiagonalScaleTrapezoid
( LeftOrRight side,
  UpperOrLower uplo,
  Orientation orientation,
  const ElementalMatrix<TDiag>& dPre,
        DistMatrix<T,U,V>& A,
  Int offset )
{
    DEBUG_ONLY(CSE cse("DiagonalScaleTrapezoid"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int mLoc = A.LocalHeight();
    const Int nLoc = A.LocalWidth();
    const bool conjugate = ( orientation==ADJOINT );

    const Int diagLength = A.DiagonalLength(offset);
    const Int ldim = A.LDim();
    T* ABuf = A.Buffer();

    const Int iOff = ( offset>=0 ? 0      : -offset );
    const Int jOff = ( offset>=0 ? offset : 0       );

    if( side == LEFT )
    {
        ElementalProxyCtrl ctrl;
        ctrl.rootConstrain = true;
        ctrl.colConstrain = true;
        ctrl.root = A.Root();
        ctrl.colAlign = A.ColAlign();

        DistMatrixReadProxy<TDiag,TDiag,U,Collect<V>()> dProx( dPre, ctrl );
        auto& d = dProx.GetLocked();

        if( uplo == LOWER )
        {
            // Scale from the left up to the diagonal
            for( Int iLoc=0; iLoc<mLoc; ++iLoc )            
            {
                const Int i = A.GlobalRow(iLoc);
                if( i >= iOff )
                {
                    const Int k = i-iOff;
                    const Int j = k+jOff;
                    const Int width = Min(j+1,n);
                    const Int localWidth = A.LocalColOffset(width);
                    const TDiag alpha = 
                        ( conjugate ? Conj(d.GetLocal(iLoc,0))
                                    :      d.GetLocal(iLoc,0) );
                    blas::Scal( localWidth, alpha, &ABuf[iLoc], ldim );
                }
            }
        }
        else
        {
            // Scale from the diagonal to the right
            for( Int iLoc=0; iLoc<mLoc; ++iLoc )
            {
                const Int i = A.GlobalRow(iLoc);
                if( i < iOff+diagLength )
                {
                    const Int k = i-iOff;
                    const Int j = k+jOff;
                    const Int jLeft = Max(j,0);
                    const Int jLeftLoc = A.LocalColOffset(jLeft);
                    const TDiag alpha = 
                        ( conjugate ? Conj(d.GetLocal(iLoc,0))
                                    :      d.GetLocal(iLoc,0) );
                    blas::Scal
                    ( nLoc-jLeftLoc, alpha, &ABuf[iLoc+jLeftLoc*ldim], ldim );
                }
            }
        }    
    }
    else
    {
        ElementalProxyCtrl ctrl;
        ctrl.rootConstrain = true;
        ctrl.colConstrain = true;
        ctrl.root = A.Root();
        ctrl.colAlign = A.RowAlign();
        
        DistMatrixReadProxy<TDiag,TDiag,V,Collect<U>()> dProx( dPre, ctrl );
        auto& d = dProx.GetLocked();

        if( uplo == LOWER )
        {
            // Scale from the diagonal downwards
            for( Int jLoc=0; jLoc<nLoc; ++jLoc )
            {
                const Int j = A.GlobalCol(jLoc);
                if( j < jOff+diagLength )
                {
                    const Int k = j-jOff;
                    const Int i = k+iOff;
                    const Int iTop = Max(i,0);
                    const Int iTopLoc = A.LocalRowOffset(iTop);
                    const TDiag alpha = 
                        ( conjugate ? Conj(d.GetLocal(jLoc,0))
                                    :      d.GetLocal(jLoc,0) );
                    blas::Scal
                    ( mLoc-iTopLoc, alpha, &ABuf[iTopLoc+jLoc*ldim], 1 );
                }
            }
        }
        else 
        {
            // Scale downward to the diagonal
            for( Int jLoc=0; jLoc<nLoc; ++jLoc )
            {
                const Int j = A.GlobalCol(jLoc);
                if( j >= jOff )
                {
                    const Int k = j-jOff;
                    const Int i = k+iOff;
                    const Int height = Min(i+1,m);
                    const Int localHeight = A.LocalRowOffset(height);
                    const TDiag alpha = 
                        ( conjugate ? Conj(d.GetLocal(jLoc,0))
                                    :      d.GetLocal(jLoc,0) );
                    blas::Scal( localHeight, alpha, &ABuf[jLoc*ldim], 1 );
                }
            }
        }
    }
}

template<typename TDiag,typename T>
void DiagonalScaleTrapezoid
( LeftOrRight side,
  UpperOrLower uplo,
  Orientation orientation,
  const ElementalMatrix<TDiag>& d,
        ElementalMatrix<T>& A,
  Int offset )
{
    DEBUG_ONLY(CSE cse("DiagonalScale"))
    #define GUARD(CDIST,RDIST) A.ColDist() == CDIST && A.RowDist() == RDIST
    #define PAYLOAD(CDIST,RDIST) \
        auto& ACast = static_cast<DistMatrix<T,CDIST,RDIST>&>(A); \
        DiagonalScaleTrapezoid( side, uplo, orientation, d, ACast, offset );
    #include "El/macros/GuardAndPayload.h"
}

template<typename TDiag,typename T>
void DiagonalScaleTrapezoid
( LeftOrRight side,
  UpperOrLower uplo,
  Orientation orientation,
  const Matrix<TDiag>& d,
        SparseMatrix<T>& A,
  Int offset )
{
    DEBUG_ONLY(CSE cse("DiagonalScaleTrapezoid"))
    LogicError("This routine is not yet written");
}

template<typename TDiag,typename T>
void DiagonalScaleTrapezoid
( LeftOrRight side,
  UpperOrLower uplo,
  Orientation orientation,
  const DistMultiVec<TDiag>& d,
        DistSparseMatrix<T>& A,
  Int offset )
{
    DEBUG_ONLY(CSE cse("DiagonalScaleTrapezoid"))
    LogicError("This routine is not yet written");
}

#ifdef EL_INSTANTIATE_BLAS_LEVEL1
# define EL_EXTERN
#else
# define EL_EXTERN extern
#endif

#define PROTO(T) \
  EL_EXTERN template void DiagonalScaleTrapezoid \
  ( LeftOrRight side, \
    UpperOrLower uplo, \
    Orientation orientation, \
    const Matrix<T>& d, \
          Matrix<T>& A, \
    Int offset ); \
  EL_EXTERN template void DiagonalScaleTrapezoid \
  ( LeftOrRight side, \
    UpperOrLower uplo, \
    Orientation orientation, \
    const ElementalMatrix<T>& d, \
          ElementalMatrix<T>& A, \
    Int offset ); \
  EL_EXTERN template void DiagonalScaleTrapezoid \
  ( LeftOrRight side, \
    UpperOrLower uplo, \
    Orientation orientation, \
    const Matrix<T>& d, \
          SparseMatrix<T>& A, \
    Int offset ); \
  EL_EXTERN template void DiagonalScaleTrapezoid \
  ( LeftOrRight side, \
    UpperOrLower uplo, \
    Orientation orientation, \
    const DistMultiVec<T>& d, \
          DistSparseMatrix<T>& A, \
    Int offset );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

#undef EL_EXTERN

} // namespace El

#endif // ifndef EL_BLAS_DIAGONALSCALETRAPEZOID_HPP
