/*
   Copyright 2009-2010 Jack Poulson

   This file is part of Elemental.

   Elemental is free software: you can redistribute it and/or modify it under
   the terms of the GNU Lesser General Public License as published by the
   Free Software Foundation; either version 3 of the License, or 
   (at your option) any later version.

   Elemental is distributed in the hope that it will be useful, but 
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with Elemental. If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef ELEMENTAL_BLAS_INTERNAL_HPP
#define ELEMENTAL_BLAS_INTERNAL_HPP 1

#include "elemental/blas.hpp"

namespace elemental {
namespace blas {
namespace internal {

//----------------------------------------------------------------//
// Distributed BLAS: Level 1                                      //
//----------------------------------------------------------------//
            
// Pseudo-partial-specializations of BLAS::Dot
template<typename T, Distribution U, Distribution V>
T
Dot
( const DistMatrix<T,U,V>& x, const DistMatrix<T,MC,MR>& y );

template<typename T, Distribution U, Distribution V>
T
Dot
( const DistMatrix<T,U,V>& x, const DistMatrix<T,MC,Star>& y );

template<typename T, Distribution U, Distribution V>
T
Dot
( const DistMatrix<T,U,V>& x, const DistMatrix<T,Star,MR>& y );
            
template<typename T, Distribution U, Distribution V>
T
Dot
( const DistMatrix<T,U,V>& x, const DistMatrix<T,MR,MC>& y );
            
template<typename T, Distribution U, Distribution V>
T
Dot
( const DistMatrix<T,U,V>& x, const DistMatrix<T,MR,Star>& y );

template<typename T, Distribution U, Distribution V>
T
Dot
( const DistMatrix<T,U,V>& x, const DistMatrix<T,Star,MC>& y );

template<typename T, Distribution U, Distribution V>
T
Dot
( const DistMatrix<T,U,V>& x, const DistMatrix<T,VC,Star>& y );

template<typename T, Distribution U, Distribution V>
T
Dot
( const DistMatrix<T,U,V>& x, const DistMatrix<T,Star,VC>& y );

template<typename T, Distribution U, Distribution V>
T
Dot
( const DistMatrix<T,U,V>& x, const DistMatrix<T,VR,Star>& y );

template<typename T, Distribution U, Distribution V>
T
Dot
( const DistMatrix<T,U,V>& x, const DistMatrix<T,Star,VR>& y );
            
template<typename T, Distribution U, Distribution V>
T
Dot
( const DistMatrix<T,U,V>& x, const DistMatrix<T,Star,Star>& y );

// Pseudo-partial-specializations of BLAS::Dotu
template<typename T, Distribution U, Distribution V>
T
Dotu
( const DistMatrix<T,U,V>& x, const DistMatrix<T,MC,MR>& y );

template<typename T, Distribution U, Distribution V>
T
Dotu
( const DistMatrix<T,U,V>& x, const DistMatrix<T,MC,Star>& y );

template<typename T, Distribution U, Distribution V>
T
Dotu
( const DistMatrix<T,U,V>& x, const DistMatrix<T,Star,MR>& y );
            
template<typename T, Distribution U, Distribution V>
T
Dotu
( const DistMatrix<T,U,V>& x, const DistMatrix<T,MR,MC>& y );
            
template<typename T, Distribution U, Distribution V>
T
Dotu
( const DistMatrix<T,U,V>& x, const DistMatrix<T,MR,Star>& y );

template<typename T, Distribution U, Distribution V>
T
Dotu
( const DistMatrix<T,U,V>& x, const DistMatrix<T,Star,MC>& y );

template<typename T, Distribution U, Distribution V>
T
Dotu
( const DistMatrix<T,U,V>& x, const DistMatrix<T,VC,Star>& y );

template<typename T, Distribution U, Distribution V>
T
Dotu
( const DistMatrix<T,U,V>& x, const DistMatrix<T,Star,VC>& y );

template<typename T, Distribution U, Distribution V>
T
Dotu
( const DistMatrix<T,U,V>& x, const DistMatrix<T,VR,Star>& y );

template<typename T, Distribution U, Distribution V>
T
Dotu
( const DistMatrix<T,U,V>& x, const DistMatrix<T,Star,VR>& y );
            
template<typename T, Distribution U, Distribution V>
T
Dotu
( const DistMatrix<T,U,V>& x, const DistMatrix<T,Star,Star>& y );

//----------------------------------------------------------------------------//
// Distributed BLAS: Level 2                                                  //
//----------------------------------------------------------------------------//
            
// Gemv where A is not transposed
template<typename T>
void
GemvN
( T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& x,
  T beta,        DistMatrix<T,MC,MR>& y );

// Gemv where A is (conjugate) transposed
template<typename T>
void
GemvT
( Orientation orientation,
  T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& x,
  T beta,        DistMatrix<T,MC,MR>& y );

// This is for the case where x is a column vector.
//
// Returns the unreduced components z[MC,* ] and z[MR,* ]:
//     z[MC,* ] := alpha tril(A)[MC,MR] x[MR,* ]
//     z[MR,* ] := alpha (trils(A)[MC,MR])^H x[MC,* ]
template<typename T>
void
HemvColAccumulate
( Shape shape,
  T alpha, 
  const DistMatrix<T,MC,MR  >& A,
  const DistMatrix<T,MC,Star>& x_MC_Star,
  const DistMatrix<T,MR,Star>& x_MR_Star,
        DistMatrix<T,MC,Star>& z_MC_Star,
        DistMatrix<T,MR,Star>& z_MR_Star
);

// This is for the case where x is a column vector and A is lower.
//
// Returns the unreduced components z[MC,* ] and z[MR,* ]:
//     z[MC,* ] := alpha tril(A)[MC,MR] x[MR,* ]
//     z[MR,* ] := alpha (trils(A)[MC,MR])^H x[MC,* ]
template<typename T>
void
HemvColAccumulateL
( T alpha, 
  const DistMatrix<T,MC,MR  >& A,
  const DistMatrix<T,MC,Star>& x_MC_Star,
  const DistMatrix<T,MR,Star>& x_MR_Star,
        DistMatrix<T,MC,Star>& z_MC_Star,
        DistMatrix<T,MR,Star>& z_MR_Star
);

// This is for the case where x is a column vector and A is upper.
//
// Returns the unreduced components z[MC,* ] and z[MR,* ]:
//     z[MC,* ] := alpha triu(A)[MC,MR] x[MR,* ]
//     z[MR,* ] := alpha (trius(A)[MC,MR])^H x[MC,* ]
template<typename T>
void
HemvColAccumulateU
( T alpha, 
  const DistMatrix<T,MC,MR  >& A,
  const DistMatrix<T,MC,Star>& x_MC_Star,
  const DistMatrix<T,MR,Star>& x_MR_Star,
        DistMatrix<T,MC,Star>& z_MC_Star,
        DistMatrix<T,MR,Star>& z_MR_Star
);

// This is for the case where x is a row vector.
//
// Returns the unreduced components z[MC,* ] and z[MR,* ]:
//     z[MC,* ] := alpha tril(A)[MC,MR] (x[* ,MR])^T
//     z[MR,* ] := alpha (trils(A)[MC,MR])^H (x[* ,MC])^T
template<typename T>
void
HemvRowAccumulate
( Shape shape,
  T alpha, 
  const DistMatrix<T,MC,  MR  >& A,
  const DistMatrix<T,Star,MC  >& x_Star_MC,
  const DistMatrix<T,Star,MR  >& x_Star_MR,
        DistMatrix<T,MC,  Star>& z_MC_Star,
        DistMatrix<T,MR,  Star>& z_MR_Star
);

// This is for the case where x is a row vector and A is lower.
//
// Returns the unreduced components z[MC,* ] and z[MR,* ]:
//     z[MC,* ] := alpha tril(A)[MC,MR] (x[* ,MR])^T
//     z[MR,* ] := alpha (trils(A)[MC,MR])^H (x[* ,MC])^T
template<typename T>
void
HemvRowAccumulateL
( T alpha, 
  const DistMatrix<T,MC,  MR  >& A,
  const DistMatrix<T,Star,MC  >& x_Star_MC,
  const DistMatrix<T,Star,MR  >& x_Star_MR,
        DistMatrix<T,MC,  Star>& z_MC_Star,
        DistMatrix<T,MR,  Star>& z_MR_Star
);

// This is for the case where x is a row vector and A is upper.
//
// Returns the unreduced components z[MC,* ] and z[MR,* ]:
//     z[MC,* ] := alpha triu(A)[MC,MR] (x[* ,MR])^T
//     z[MR,* ] := alpha (trius(A)[MC,MR])^H (x[* ,MC])^T
template<typename T>
void
HemvRowAccumulateU
( T alpha, 
  const DistMatrix<T,MC,  MR  >& A,
  const DistMatrix<T,Star,MC  >& x_Star_MC,
  const DistMatrix<T,Star,MR  >& x_Star_MR,
        DistMatrix<T,MC,  Star>& z_MC_Star,
        DistMatrix<T,MR,  Star>& z_MR_Star
);


// This is for the case where x is a column vector.
//
// Returns the unreduced components z[MC,* ] and z[MR,* ]:
//     z[MC,* ] := alpha tril(A)[MC,MR] x[MR,* ]
//     z[MR,* ] := alpha (trils(A)[MC,MR])^T x[MC,* ]
template<typename T>
void
SymvColAccumulate
( Shape shape,
  T alpha, 
  const DistMatrix<T,MC,MR  >& A,
  const DistMatrix<T,MC,Star>& x_MC_Star,
  const DistMatrix<T,MR,Star>& x_MR_Star,
        DistMatrix<T,MC,Star>& z_MC_Star,
        DistMatrix<T,MR,Star>& z_MR_Star
);

// This is for the case where x is a column vector and A is lower.
//
// Returns the unreduced components z[MC,* ] and z[MR,* ]:
//     z[MC,* ] := alpha tril(A)[MC,MR] x[MR,* ]
//     z[MR,* ] := alpha (trils(A)[MC,MR])^T x[MC,* ]
template<typename T>
void
SymvColAccumulateL
( T alpha, 
  const DistMatrix<T,MC,MR  >& A,
  const DistMatrix<T,MC,Star>& x_MC_Star,
  const DistMatrix<T,MR,Star>& x_MR_Star,
        DistMatrix<T,MC,Star>& z_MC_Star,
        DistMatrix<T,MR,Star>& z_MR_Star
);

// This is for the case where x is a column vector and A is upper.
//
// Returns the unreduced components z[MC,* ] and z[MR,* ]:
//     z[MC,* ] := alpha triu(A)[MC,MR] x[MR,* ]
//     z[MR,* ] := alpha (trius(A)[MC,MR])^T x[MC,* ]
template<typename T>
void
SymvColAccumulateU
( T alpha, 
  const DistMatrix<T,MC,MR  >& A,
  const DistMatrix<T,MC,Star>& x_MC_Star,
  const DistMatrix<T,MR,Star>& x_MR_Star,
        DistMatrix<T,MC,Star>& z_MC_Star,
        DistMatrix<T,MR,Star>& z_MR_Star
);

// This is for the case where x is a row vector.
//
// Returns the unreduced components z[MC,* ] and z[MR,* ]:
//     z[MC,* ] := alpha tril(A)[MC,MR] (x[* ,MR])^T
//     z[MR,* ] := alpha (trils(A)[MC,MR])^T (x[* ,MC])^T
template<typename T>
void
SymvRowAccumulate
( Shape shape,
  T alpha, 
  const DistMatrix<T,MC,  MR>& A,
  const DistMatrix<T,Star,MC>& x_Star_MC,
  const DistMatrix<T,Star,MR>& x_Star_MR,
        DistMatrix<T,Star,MC>& z_Star_MC,
        DistMatrix<T,Star,MR>& z_Star_MR
);

// This is for the case where x is a row vector and A is lower.
//
// Returns the unreduced components z[MC,* ] and z[MR,* ]:
//     z[MC,* ] := alpha tril(A)[MC,MR] (x[* ,MR])^T
//     z[MR,* ] := alpha (trils(A)[MC,MR])^T (x[* ,MC])^T
template<typename T>
void
SymvRowAccumulateL
( T alpha, 
  const DistMatrix<T,MC,  MR>& A,
  const DistMatrix<T,Star,MC>& x_Star_MC,
  const DistMatrix<T,Star,MR>& x_Star_MR,
        DistMatrix<T,Star,MC>& z_Star_MC,
        DistMatrix<T,Star,MR>& z_Star_MR
);

// This is for the case where x is a row vector and A is upper.
//
// Returns the unreduced components z[MC,* ] and z[MR,* ]:
//     z[MC,* ] := alpha triu(A)[MC,MR] (x[* ,MR])^T
//     z[MR,* ] := alpha (trius(A)[MC,MR])^T (x[* ,MC])^T
template<typename T>
void
SymvRowAccumulateU
( T alpha, 
  const DistMatrix<T,MC,  MR>& A,
  const DistMatrix<T,Star,MC>& x_Star_MC,
  const DistMatrix<T,Star,MR>& x_Star_MR,
        DistMatrix<T,Star,MC>& z_Star_MC,
        DistMatrix<T,Star,MR>& z_Star_MR
);

template<typename T>
void
TrsvLN
( Diagonal diagonal, const DistMatrix<T,MC,MR>& L, DistMatrix<T,MC,MR>& x );

template<typename T>
void
TrsvLT
( Orientation orientation, Diagonal diagonal,
  const DistMatrix<T,MC,MR>& L, DistMatrix<T,MC,MR>& x );

template<typename T>
void
TrsvUN
( Diagonal diagonal, const DistMatrix<T,MC,MR>& U, DistMatrix<T,MC,MR>& x );

template<typename T>
void
TrsvUT
( Orientation orientation, Diagonal diagonal,
  const DistMatrix<T,MC,MR>& U, DistMatrix<T,MC,MR>& x );

//----------------------------------------------------------------------------//
// Distributed BLAS: Level 3                                                  //
//----------------------------------------------------------------------------//

// Gemm where we avoid redistributing A.
template<typename T>
void
GemmA
( Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

// Gemm where we avoid redistributing B.
template<typename T>
void
GemmB
( Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

// Gemm where we avoid redistributing C.
template<typename T>
void
GemmC
( Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

// Gemm for panel-panel dot products.
template<typename T>
void
GemmDot
( Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

// Normal Normal Gemm.
template<typename T>
void
GemmNN
( T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

// Normal Normal Gemm where we avoid redistributing A.
template<typename T>
void
GemmNNA
( T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

// Normal Normal Gemm where we avoid redistributing B.
template<typename T>
void
GemmNNB
( T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

// Normal Normal Gemm where we avoid redistributing C.
template<typename T>
void
GemmNNC
( T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

// Normal Normal Gemm for panel-panel dot product
template<typename T>
void
GemmNNDot
( T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

// Normal (Conjugate)Transpose Gemm.
template<typename T>
void
GemmNT
( Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

// Normal (Conjugate)Transpose Gemm where we avoid redistributing A.
template<typename T>
void
GemmNTA
( Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

// Normal (Conjugate)Transpose Gemm where we avoid redistributing B.
template<typename T>
void
GemmNTB
( Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

// Normal (Conjugate)Transpose Gemm where we avoid redistributing C.
template<typename T>
void
GemmNTC
( Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

// Normal (Conjugate)Transpose Gemm for panel-panel dot product
template<typename T>
void
GemmNTDot
( Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );
 
// (Conjugate)Transpose Normal Gemm.
template<typename T>
void
GemmTN
( Orientation orientationOfA,
  T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

// (Conjugate)Transpose Normal Gemm where we avoid redistributing A.
template<typename T>
void
GemmTNA
( Orientation orientationOfA,
  T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

// (Conjugate)Transpose Normal Gemm where we avoid redistributing B.
template<typename T>
void
GemmTNB
( Orientation orientationOfA,
  T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

// (Conjugate)Transpose Normal Gemm where we avoid redistributing C.
template<typename T>
void
GemmTNC
( Orientation orientationOfA,
  T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

// (Conjugate)Transpose Normal Gemm for panel-panel dot product
template<typename T>
void
GemmTNDot
( Orientation orientationOfA,
  T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

// (Conjugate)Transpose (Conjugate)Transpose Gemm.
template<typename T>
void
GemmTT
( Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

// (Conjugate)Transpose (Conjugate)Transpose Gemm where we avoid 
// redistributing A.
template<typename T>
void
GemmTTA
( Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

// (Conjugate)Transpose (Conjugate)Transpose Gemm where we avoid 
// redistributing B.
template<typename T>
void
GemmTTB
( Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

// (Conjugate)Transpose (Conjugate)Transpose Gemm where we avoid 
// redistributing C.
template<typename T>
void
GemmTTC
( Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

// (Conjugate)Transpose (Conjugate)Transpose Gemm for panel-panel
// dot product
template<typename T>
void
GemmTTDot
( Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

// Hemm
template<typename T>
void
Hemm
( Side side, Shape shape,
  T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

// Left Lower Hemm
template<typename T>
void
HemmLL
( T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

// Left Lower Hemm where we avoid redistributing C
template<typename T>
void
HemmLLC
( T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

// Left Upper Hemm
template<typename T>
void
HemmLU
( T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

// Left Upper Hemm where we avoid redistributing C
template<typename T>
void
HemmLUC
( T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

// Right Lower Hemm
template<typename T>
void
HemmRL
( T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

// Right Lower Hemm where we avoid redistributing C
template<typename T>
void
HemmRLC
( T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

// Right Upper Hemm
template<typename T>
void
HemmRU
( T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

// Right Upper Hemm where we avoid redistributing C
template<typename T>
void
HemmRUC
( T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

// Lower, Normal Her2k
template<typename T>
void
Her2kLN
( T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

// Lower, ConjugateTranspose Her2k
template<typename T>
void
Her2kLC
( T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

// Upper, Normal Her2k
template<typename T>
void
Her2kUN
( T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

// Upper, ConjugateTranspose Her2k
template<typename T>
void
Her2kUC
( T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

// Lower, Normal Herk
template<typename T>
void
HerkLN
( T alpha, const DistMatrix<T,MC,MR>& A, T beta, DistMatrix<T,MC,MR>& C );

// Lower, ConjugateTranspose Herk
template<typename T>
void
HerkLC
( T alpha, const DistMatrix<T,MC,MR>& A, T beta, DistMatrix<T,MC,MR>& C );

// Upper, Normal Herk
template<typename T>
void
HerkUN
( T alpha, const DistMatrix<T,MC,MR>& A, T beta, DistMatrix<T,MC,MR>& C );

// Upper, ConjugateTranspose Herk
template<typename T>
void
HerkUC
( T alpha, const DistMatrix<T,MC,MR>& A, T beta, DistMatrix<T,MC,MR>& C );

// Symm
template<typename T>
void
Symm
( Side side, Shape shape,
  T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

// Left Lower Symm
template<typename T>
void
SymmLL
( T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

// Left Lower Symm where we avoid redistributing C
template<typename T>
void
SymmLLC
( T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

// Left Upper Symm
template<typename T>
void
SymmLU
( T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

// Left Upper Symm where we avoid redistributing C
template<typename T>
void
SymmLUC
( T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

// Right Lower Symm
template<typename T>
void
SymmRL
( T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

// Right Lower Symm where we avoid redistributing C
template<typename T>
void
SymmRLC
( T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

// Right Upper Symm
template<typename T>
void
SymmRU
( T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

// Right Upper Symm where we avoid redistributing C
template<typename T>
void
SymmRUC
( T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

// Lower, Normal Syr2k
template<typename T>
void
Syr2kLN
( T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

// Lower, Transpose Syr2k
template<typename T>
void
Syr2kLT
( T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

// Upper, Normal Syr2k
template<typename T>
void
Syr2kUN
( T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

// Upper, Transpose Syr2k
template<typename T>
void
Syr2kUT
( T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C );

// Lower, Normal Syrk
template<typename T>
void
SyrkLN
( T alpha, const DistMatrix<T,MC,MR>& A, T beta, DistMatrix<T,MC,MR>& C );

// Lower, Transpose Syrk
template<typename T>
void
SyrkLT
( T alpha, const DistMatrix<T,MC,MR>& A, T beta, DistMatrix<T,MC,MR>& C );

// Upper, Normal Syrk
template<typename T>
void
SyrkUN
( T alpha, const DistMatrix<T,MC,MR>& A, T beta, DistMatrix<T,MC,MR>& C );

// Upper, Transpose Syrk
template<typename T>
void
SyrkUT
( T alpha, const DistMatrix<T,MC,MR>& A, T beta, DistMatrix<T,MC,MR>& C );

// Triangular Rank-K Update:
// tril(C) := alpha tril( A*B ) + beta tril(C)
//   or 
// triu(C) := alpha triu( A*B ) + beta triu(C)
template<typename T>
void
TriangularRankK
( Shape shape,
  T alpha, const DistMatrix<T,MC,Star>& A, const DistMatrix<T,Star,MR>& B,
  T beta,        DistMatrix<T,MC,MR  >& C );

// Triangular Rank-K Update:
// tril(C) := alpha tril( A^T*B ) + beta tril(C)
//   or 
// triu(C) := alpha triu( A^T*B ) + beta triu(C)
template<typename T>
void
TriangularRankK
( Shape shape,
  T alpha, const DistMatrix<T,Star,MC>& A, const DistMatrix<T,Star,MR>& B,
  T beta,        DistMatrix<T,MC,  MR>& C );

// Triangular Rank-K Update:
// tril(C) := alpha tril( A^T*B^T ) + beta tril(C)
//   or 
// triu(C) := alpha triu( A^T*B^T ) + beta triu(C)
template<typename T>
void
TriangularRankK
( Shape shape,
  T alpha, const DistMatrix<T,Star,MC>& A, const DistMatrix<T,MR,Star>& B,
  T beta,        DistMatrix<T,MC,  MR>& C );

// Triangular Rank-2K Update:
// tril(C) := alpha tril( A1*B1^T + A2*B2 ) + beta tril(C)
//   or
// triu(C) := alpha triu( A1*B1^T + A2*B2 ) + beta triu(C)
template<typename T>
void
TriangularRank2K
( Shape shape,
  T alpha, const DistMatrix<T,MC,Star>& A1, const DistMatrix<T,MC,Star>& A2, 
           const DistMatrix<T,MR,Star>& B1, const DistMatrix<T,Star,MR>& B2,
  T beta,        DistMatrix<T,MC,MR  >& C );

// Triangular Rank-2K Update:
// tril(C) := alpha tril( A1*B1 + A2*B2 ) + beta tril(C)
//   or
// triu(C) := alpha triu( A1*B1 + A2*B2 ) + beta triu(C)
template<typename T>
void
TriangularRank2K
( Shape shape,
  T alpha, const DistMatrix<T,MC,Star>& A1, const DistMatrix<T,MC,Star>& A2, 
           const DistMatrix<T,Star,MR>& B1, const DistMatrix<T,Star,MR>& B2,
  T beta,        DistMatrix<T,MC,  MR>& C );

// Triangular Rank-2K Update:
// tril(C) := alpha tril( A1^T*B1 + A2^T*B2 ) + beta tril(C)
//   or
// triu(C) := alpha triu( A1^T*B1 + A2^T*B2 ) + beta triu(C)
template<typename T>
void
TriangularRank2K
( Shape shape,
  T alpha, const DistMatrix<T,Star,MC>& A1, const DistMatrix<T,Star,MC>& A2, 
           const DistMatrix<T,Star,MR>& B1, const DistMatrix<T,Star,MR>& B2,
  T beta,        DistMatrix<T,MC,  MR>& C );

// Triangular Rank-2K Update:
// tril(C) := alpha tril( A1^T*B1^T + A2^T*B2^T ) + beta tril(C)
//   or
// triu(C) := alpha triu( A1^T*B1^T + A2^T*B2^T ) + beta triu(C)
template<typename T>
void
TriangularRank2K
( Shape shape,
  T alpha, const DistMatrix<T,Star,MC>& A1, const DistMatrix<T,Star,MC>& A2, 
           const DistMatrix<T,MR,Star>& B1, const DistMatrix<T,MR,Star>& B2,
  T beta,        DistMatrix<T,MC,  MR>& C );

// Left, Lower, Normal Trmm
template<typename T>
void
TrmmLLN
( Diagonal diagonal,
  T alpha, const DistMatrix<T,MC,MR>& L, DistMatrix<T,MC,MR>& X );

// Left, Lower, (Conjugate)Transpose Trmm
template<typename T>
void
TrmmLLT
( Orientation orientation, Diagonal diagonal,
  T alpha, const DistMatrix<T,MC,MR>& L, DistMatrix<T,MC,MR>& X );

// Left, Upper, Normal Trmm
template<typename T>
void
TrmmLUN
( Diagonal diagonal,
  T alpha, const DistMatrix<T,MC,MR>& U, DistMatrix<T,MC,MR>& X );

// Left, Upper, (Conjugate)Transpose Trmm
template<typename T>
void
TrmmLUT
( Orientation orientation, Diagonal diagonal,
  T alpha, const DistMatrix<T,MC,MR>& U, DistMatrix<T,MC,MR>& X );

// Right, Lower, Normal Trmm
template<typename T>
void
TrmmRLN
( Diagonal diagonal,
  T alpha, const DistMatrix<T,MC,MR>& L, DistMatrix<T,MC,MR>& X );

// Right, Lower, (Conjugate)Transpose Trmm
template<typename T>
void
TrmmRLT
( Orientation orientation, Diagonal diagonal,
  T alpha, const DistMatrix<T,MC,MR>& L, DistMatrix<T,MC,MR>& X );

// Right, Upper, Normal Trmm
template<typename T>
void
TrmmRUN
( Diagonal diagonal,
  T alpha, const DistMatrix<T,MC,MR>& U, DistMatrix<T,MC,MR>& X );

// Right, Upper, (Conjugate)Transpose Trmm
template<typename T>
void
TrmmRUT
( Orientation orientation, Diagonal diagonal,
  T alpha, const DistMatrix<T,MC,MR>& U, DistMatrix<T,MC,MR>& X );

// Left, Lower, Normal Trsm
template<typename T>
void
TrsmLLN
( Diagonal diagonal,
  T alpha, const DistMatrix<T,MC,MR>& L, DistMatrix<T,MC,MR>& X );

// Left, Lower, (Conjugate)Transpose Trsm
template<typename T>
void
TrsmLLT
( Orientation orientation, Diagonal diagonal,
  T alpha, const DistMatrix<T,MC,MR>& L, DistMatrix<T,MC,MR>& X );

// Left, Upper, Normal Trsm
template<typename T>
void
TrsmLUN
( Diagonal diagonal,
  T alpha, const DistMatrix<T,MC,MR>& U, DistMatrix<T,MC,MR>& X );

// Left, Upper, (Conjugate)Transpose Trsm
template<typename T>
void
TrsmLUT
( Orientation orientation, Diagonal diagonal,
  T alpha, const DistMatrix<T,MC,MR>& U, DistMatrix<T,MC,MR>& X );

// Right, Lower, Normal Trsm
template<typename T>
void
TrsmRLN
( Diagonal diagonal,
  T alpha, const DistMatrix<T,MC,MR>& L, DistMatrix<T,MC,MR>& X );

// Right, Lower, (Conjugate)Transpose Trsm
template<typename T>
void
TrsmRLT
( Orientation orientation, Diagonal diagonal,
  T alpha, const DistMatrix<T,MC,MR>& L, DistMatrix<T,MC,MR>& X );

// Right, Upper, Normal Trsm
template<typename T>
void
TrsmRUN
( Diagonal diagonal,
  T alpha, const DistMatrix<T,MC,MR>& U, DistMatrix<T,MC,MR>& X );

// Right, Upper, (Conjugate)Transpose Trsm
template<typename T>
void
TrsmRUT
( Orientation orientation, Diagonal diagonal,
  T alpha, const DistMatrix<T,MC,MR>& U, DistMatrix<T,MC,MR>& X );

//----------------------------------------------------------------------------//
// Level 2 BLAS Utility Functions                                             //
//----------------------------------------------------------------------------//
template<typename T>
double
SymvGFlops
( int m, double seconds );
 
//----------------------------------------------------------------------------//
// Level 3 BLAS Utility Functions                                             //
//----------------------------------------------------------------------------//
template<typename T>
double 
GemmGFlops
( int m, int n, int k, double seconds );

template<typename T>
double
HemmGFlops
( Side side, int m, int n, double seconds );

template<typename T>
double
Her2kGFlops
( int m, int k, double seconds );

template<typename T>
double
HerkGFlops
( int m, int k, double seconds );
            
template<typename T>
double
SymmGFlops
( Side side, int m, int n, double seconds );
            
template<typename T>
double
Syr2kGFlops
( int m, int k, double seconds );

template<typename T>
double
SyrkGFlops
( int m, int k, double seconds );
 
template<typename T>
double
TrmmGFlops
( Side side, int m, int n, double seconds );
            
template<typename T>
double
TrsmGFlops
( Side side, int m, int n, double seconds );

} // internal
} // blas
} // elemental

//----------------------------------------------------------------------------//
// Implementations begin here                                                 //
//----------------------------------------------------------------------------//

namespace elemental {
namespace blas {
namespace internal {

// Level 2 Utility functions
template<>
inline double
SymvGFlops<float>
( int m, double seconds )
{ return (1.*m*m)/(1.e9*seconds); }

template<>
inline double
SymvGFlops<double>
( int m, double seconds )
{ return SymvGFlops<float>(m,seconds); }

#ifndef WITHOUT_COMPLEX
template<>
inline double
SymvGFlops<scomplex>
( int m, double seconds )
{ return 4.*SymvGFlops<float>(m,seconds); }

template<>
inline double
SymvGFlops<dcomplex>
( int m, double seconds )
{ return 4.*SymvGFlops<float>(m,seconds); }
#endif

// Level 3 Utility functions
template<>
inline double
GemmGFlops<float>
( int m, int n, int k, double seconds )
{ return (2.*m*n*k)/(1.e9*seconds); }

template<>
inline double
GemmGFlops<double>
( int m, int n, int k, double seconds )
{ return GemmGFlops<float>(m,n,k,seconds); }

#ifndef WITHOUT_COMPLEX
template<>
inline double
GemmGFlops<scomplex>
( int m, int n, int k, double seconds )
{ return 4.*GemmGFlops<float>(m,n,k,seconds); }

template<>
inline double
GemmGFlops<dcomplex>
( int m, int n, int k, double seconds )
{ return 4.*GemmGFlops<float>(m,n,k,seconds); }
#endif

template<>
inline double
HemmGFlops<float>
( Side side, int m, int n, double seconds )
{
    if( side == Left )
        return (2.*m*m*n)/(1.e9*seconds);
    else
        return (2.*m*n*n)/(1.e9*seconds);
}

template<>
inline double
HemmGFlops<double>
( Side side, int m, int n, double seconds )
{ return HemmGFlops<float>(side,m,n,seconds); }

#ifndef WITHOUT_COMPLEX
template<>
inline double
HemmGFlops<scomplex>
( Side side, int m, int n, double seconds )
{ return 4.*HemmGFlops<float>(side,m,n,seconds); }

template<>
inline double
HemmGFlops<dcomplex>
( Side side, int m, int n, double seconds )
{ return 4.*HemmGFlops<float>(side,m,n,seconds); }
#endif

template<>
inline double
Her2kGFlops<float>
( int m, int k, double seconds )
{ return (2.*m*m*k)/(1.e9*seconds); }

template<>
inline double
Her2kGFlops<double>
( int m, int k, double seconds )
{ return Her2kGFlops<float>(m,k,seconds); }

#ifndef WITHOUT_COMPLEX
template<>
inline double
Her2kGFlops<scomplex>
( int m, int k, double seconds )
{ return 4.*Her2kGFlops<float>(m,k,seconds); }

template<>
inline double
Her2kGFlops<dcomplex>
( int m, int k, double seconds )
{ return 4.*Her2kGFlops<float>(m,k,seconds); }
#endif

template<>
inline double
HerkGFlops<float>
( int m, int k, double seconds )
{ return (1.*m*m*k)/(1.e9*seconds); }

template<>
inline double
HerkGFlops<double>
( int m, int k, double seconds )
{ return HerkGFlops<float>(m,k,seconds); }

#ifndef WITHOUT_COMPLEX
template<>
inline double
HerkGFlops<scomplex>
( int m, int k, double seconds )
{ return 4.*HerkGFlops<float>(m,k,seconds); }

template<>
inline double
HerkGFlops<dcomplex>
( int m, int k, double seconds )
{ return 4.*HerkGFlops<float>(m,k,seconds); }
#endif
            
template<>
inline double
SymmGFlops<float>
( Side side, int m, int n, double seconds )
{
    if( side == Left )
        return (2.*m*m*n)/(1.e9*seconds);
    else
        return (2.*m*n*n)/(1.e9*seconds);
}

template<>
inline double
SymmGFlops<double>
( Side side, int m, int n, double seconds ) 
{ return SymmGFlops<float>(side,m,n,seconds); }
            
#ifndef WITHOUT_COMPLEX
template<>
inline double
SymmGFlops<scomplex>
( Side side, int m, int n, double seconds )
{ return 4.*SymmGFlops<float>(side,m,n,seconds); }

template<>
inline double
SymmGFlops<dcomplex>
( Side side, int m, int n, double seconds )
{ return 4.*SymmGFlops<float>(side,m,n,seconds); }
#endif
            
template<>
inline double
Syr2kGFlops<float>
( int m, int k, double seconds )
{ return (2.*m*m*k)/(1.e9*seconds); }

template<>
inline double
Syr2kGFlops<double>
( int m, int k, double seconds )
{ return Syr2kGFlops<float>(m,k,seconds); }
            
#ifndef WITHOUT_COMPLEX
template<>
inline double
Syr2kGFlops<scomplex>
( int m, int k, double seconds )
{ return 4.*Syr2kGFlops<float>(m,k,seconds); }

template<>
inline double
Syr2kGFlops<dcomplex>
( int m, int k, double seconds )
{ return 4.*Syr2kGFlops<float>(m,k,seconds); }
#endif
            
template<>
inline double
SyrkGFlops<float>
( int m, int k, double seconds )
{ return (1.*m*m*k)/(1.e9*seconds); }

template<>
inline double
SyrkGFlops<double>
( int m, int k, double seconds )
{ return SyrkGFlops<float>(m,k,seconds); }
            
#ifndef WITHOUT_COMPLEX
template<>
inline double
SyrkGFlops<scomplex>
( int m, int k, double seconds )
{ return 4.*SyrkGFlops<float>(m,k,seconds); }
            
template<>
inline double
SyrkGFlops<dcomplex>
( int m, int k, double seconds )
{ return 4.*SyrkGFlops<float>(m,k,seconds); }
#endif
            
template<>
inline double
TrmmGFlops<float>
( Side side, int m, int n, double seconds )
{
    if( side == Left )
        return (1.*m*m*n)/(1.e9*seconds);
    else
        return (1.*m*n*n)/(1.e9*seconds);
}

template<>
inline double
TrmmGFlops<double>
( Side side, int m, int n, double seconds )
{ return TrmmGFlops<float>(side,m,n,seconds); }

#ifndef WITHOUT_COMPLEX
template<>
inline double
TrmmGFlops<scomplex>
( Side side, int m, int n, double seconds )
{ return 4.*TrmmGFlops<float>(side,m,n,seconds); }

template<>
inline double
TrmmGFlops<dcomplex>
( Side side, int m, int n, double seconds )
{ return 4.*TrmmGFlops<float>(side,m,n,seconds); }
#endif
            
template<>
inline double
TrsmGFlops<float>
( Side side, int m, int n, double seconds )
{
    if( side == Left )
        return (1.*m*m*n)/(1.e9*seconds);
    else
        return (1.*m*n*n)/(1.e9*seconds);
}

template<>
inline double
TrsmGFlops<double>
( Side side, int m, int n, double seconds )
{ return TrsmGFlops<float>(side,m,n,seconds); }
            
#ifndef WITHOUT_COMPLEX
template<>
inline double
TrsmGFlops<scomplex>
( Side side, int m, int n, double seconds )
{ return 4.*TrsmGFlops<float>(side,m,n,seconds); }

template<>
inline double
TrsmGFlops<dcomplex>
( Side side, int m, int n, double seconds )
{ return 4.*TrsmGFlops<float>(side,m,n,seconds); }
#endif
            
} // internal
} // blas
} // elemental

#endif /* ELEMENTAL_BLAS_INTERNAL_HPP */

