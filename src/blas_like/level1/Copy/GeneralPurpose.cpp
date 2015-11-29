/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#include "El/blas_like/level1/copy_internal.hpp"

namespace El {
namespace copy {

template<typename T>
void GeneralPurpose
( const AbstractDistMatrix<T>& A,
        AbstractDistMatrix<T>& B ) 
{
    DEBUG_ONLY(CSE cse("copy::GeneralPurpose"))
    const Int height = A.Height();
    const Int width = A.Width();
    if( A.Grid().Size() == 1 )
    {
        B.Resize( height, width );
        Copy( A.LockedMatrix(), B.Matrix() );
    }
    else
    {
        const Int localHeight = A.LocalHeight();
        const Int localWidth = A.LocalWidth();
        Zeros( B, height, width );

        // TODO: Break into smaller pieces to avoid excessive memory usage?
        if( A.RedundantRank() == 0 )
        {
            vector<Int> globalRows(localHeight);
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                globalRows[iLoc] = A.GlobalRow(iLoc);

            B.Reserve( localHeight*localWidth );
            const T* ABuf = A.LockedBuffer();
            const Int ALDim = A.LDim();
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const Int j = A.GlobalCol(jLoc);
                for( Int iLoc=0; iLoc<localHeight; ++iLoc ) 
                {
                    const Int i = globalRows[iLoc];
                    B.QueueUpdate( i, j, ABuf[iLoc+jLoc*ALDim] );
                }
            }    
        }

        const bool includeViewers = (A.Grid() != B.Grid());
        B.ProcessQueues( includeViewers );
    }
}

template<typename S,typename T,typename>
void GeneralPurpose
( const AbstractDistMatrix<S>& A,
        AbstractDistMatrix<T>& B ) 
{
    DEBUG_ONLY(CSE cse("copy::GeneralPurpose"))
    const Int height = A.Height();
    const Int width = A.Width();
    if( A.Grid().Size() == 1 )
    {
        B.Resize( height, width );
        Copy( A.LockedMatrix(), B.Matrix() );
    }
    else
    {
        const Int localHeight = A.LocalHeight();
        const Int localWidth = A.LocalWidth();
        Zeros( B, height, width );

        // TODO: Break into smaller pieces to avoid excessive memory usage?
        if( A.RedundantRank() == 0 )
        {
            vector<Int> globalRows(localHeight);
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                globalRows[iLoc] = A.GlobalRow(iLoc);

            B.Reserve( localHeight*localWidth );
            const S* ABuf = A.LockedBuffer();
            const Int ALDim = A.LDim();
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const Int j = A.GlobalCol(jLoc);
                for( Int iLoc=0; iLoc<localHeight; ++iLoc ) 
                {
                    const Int i = globalRows[iLoc];
                    B.QueueUpdate( i, j, T(ABuf[iLoc+jLoc*ALDim]) );
                }
            }    
        }

        const bool includeViewers = (A.Grid() != B.Grid());
        B.ProcessQueues( includeViewers );
    }
}

#define CONVERT(S,T) \
  template void GeneralPurpose \
  ( const AbstractDistMatrix<S>& A, \
          AbstractDistMatrix<T>& B );

#define SAME(T) CONVERT(T,T)

#define PROTO_INT(T) SAME(T) 

#define PROTO_REAL(Real) \
  SAME(Real) \
  CONVERT(Int,Real) \
  CONVERT(Real,Complex<Real>)

#define PROTO_COMPLEX(C) \
  SAME(C) \
  CONVERT(Int,C)

#ifdef EL_HAVE_QUAD

#define PROTO_FLOAT \
  PROTO_REAL(float) \
  CONVERT(float,double) \
  CONVERT(float,Quad) \
  CONVERT(float,Complex<double>) \
  CONVERT(float,Complex<Quad>)

#define PROTO_DOUBLE \
  PROTO_REAL(double) \
  CONVERT(double,float) \
  CONVERT(double,Quad) \
  CONVERT(double,Complex<float>) \
  CONVERT(double,Complex<Quad>)

#define PROTO_QUAD \
  PROTO_REAL(Quad) \
  CONVERT(Quad,float) \
  CONVERT(Quad,double) \
  CONVERT(Quad,Complex<float>) \
  CONVERT(Quad,Complex<double>)

#define PROTO_COMPLEX_FLOAT \
  PROTO_COMPLEX(Complex<float>) \
  CONVERT(Complex<float>,Complex<double>) \
  CONVERT(Complex<float>,Complex<Quad>)

#define PROTO_COMPLEX_DOUBLE \
  PROTO_COMPLEX(Complex<double>) \
  CONVERT(Complex<double>,Complex<float>) \
  CONVERT(Complex<double>,Complex<Quad>)

#define PROTO_COMPLEX_QUAD \
  PROTO_COMPLEX(Complex<Quad>) \
  CONVERT(Complex<Quad>,Complex<float>) \
  CONVERT(Complex<Quad>,Complex<double>)

#else

#define PROTO_FLOAT \
  PROTO_REAL(float) \
  CONVERT(float,double) \
  CONVERT(float,Complex<double>)

#define PROTO_DOUBLE \
  PROTO_REAL(double) \
  CONVERT(double,float) \
  CONVERT(double,Complex<float>)

#define PROTO_COMPLEX_FLOAT \
  PROTO_COMPLEX(Complex<float>) \
  CONVERT(Complex<float>,Complex<double>)

#define PROTO_COMPLEX_DOUBLE \
  PROTO_COMPLEX(Complex<double>) \
  CONVERT(Complex<double>,Complex<float>)

#endif

#define EL_ENABLE_QUAD
#include "El/macros/Instantiate.h"

} // namespace copy
} // namespace El
