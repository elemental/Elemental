/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>
#include <El/blas_like/level1.hpp>
#include <El-lite.h>
#include <El/blas_like/level1.h>
using namespace El;

extern "C" {

ElError ElCopyGraph( ElConstGraph A, ElGraph B )
{ EL_TRY( Copy( *CReflect(A), *CReflect(B) ) ) }
ElError ElCopyDistGraph( ElConstDistGraph A, ElDistGraph B )
{ EL_TRY( Copy( *CReflect(A), *CReflect(B) ) ) }

ElError ElCopyGraphFromRoot( ElConstDistGraph GDist, ElGraph G )
{ EL_TRY( CopyFromRoot( *CReflect(GDist), *CReflect(G) ) ) }
ElError ElCopyGraphFromNonRoot( ElConstDistGraph GDist, int root )
{ EL_TRY( CopyFromNonRoot( *CReflect(GDist), root ) ) }

ElError ElGetSubgraph
( ElConstGraph graph, ElRange_i I, ElRange_i J, ElGraph subgraph )
{ EL_TRY( GetSubgraph
          ( *CReflect(graph), CReflect(I), CReflect(J), 
            *CReflect(subgraph) ) ) }
ElError ElGetSubgraphDist
( ElConstDistGraph graph, ElRange_i I, ElRange_i J, ElDistGraph subgraph )
{ EL_TRY( GetSubgraph
          ( *CReflect(graph), CReflect(I), CReflect(J), 
            *CReflect(subgraph) ) ) }

#define C_PROTO_BASE(SIG,SIGBASE,T) \
  /* Y := alpha X + Y */ \
  ElError ElAxpy_ ## SIG \
  ( CREFLECT(T) alpha, ElConstMatrix_ ## SIG X, ElMatrix_ ## SIG Y ) \
  { EL_TRY( Axpy( CReflect(alpha), *CReflect(X), *CReflect(Y) ) ) } \
  ElError ElAxpyDist_ ## SIG \
  ( CREFLECT(T) alpha, ElConstDistMatrix_ ## SIG X, ElDistMatrix_ ## SIG Y ) \
  { EL_TRY( Axpy( CReflect(alpha), *CReflect(X), *CReflect(Y) ) ) } \
  ElError ElAxpySparse_ ## SIG \
  ( CREFLECT(T) alpha, \
    ElConstSparseMatrix_ ## SIG X, ElSparseMatrix_ ## SIG Y ) \
  { EL_TRY( Axpy( CReflect(alpha), *CReflect(X), *CReflect(Y) ) ) } \
  ElError ElAxpyDistSparse_ ## SIG \
  ( CREFLECT(T) alpha, \
    ElConstDistSparseMatrix_ ## SIG X, ElDistSparseMatrix_ ## SIG Y ) \
  { EL_TRY( Axpy( CReflect(alpha), *CReflect(X), *CReflect(Y) ) ) } \
  ElError ElAxpyDistMultiVec_ ## SIG \
  ( CREFLECT(T) alpha, \
    ElConstDistMultiVec_ ## SIG X, ElDistMultiVec_ ## SIG Y ) \
  { EL_TRY( Axpy( CReflect(alpha), *CReflect(X), *CReflect(Y) ) ) } \
  /* tri(Y) := tri(alpha X + Y) */ \
  ElError ElAxpyTrapezoid_ ## SIG \
  ( ElUpperOrLower uplo, CREFLECT(T) alpha, \
    ElConstMatrix_ ## SIG X, ElMatrix_ ## SIG Y, ElInt offset ) \
  { EL_TRY \
    ( AxpyTrapezoid \
      ( CReflect(uplo), CReflect(alpha), \
        *CReflect(X), *CReflect(Y), offset ) ) } \
  ElError ElAxpyTrapezoidDist_ ## SIG \
  ( ElUpperOrLower uplo, CREFLECT(T) alpha, \
    ElConstDistMatrix_ ## SIG X, ElDistMatrix_ ## SIG Y, ElInt offset ) \
  { EL_TRY \
    ( AxpyTrapezoid \
      ( CReflect(uplo), CReflect(alpha), \
        *CReflect(X), *CReflect(Y), offset ) ) } \
  ElError ElAxpyTrapezoidSparse_ ## SIG \
  ( ElUpperOrLower uplo, CREFLECT(T) alpha, \
    ElConstSparseMatrix_ ## SIG X, ElSparseMatrix_ ## SIG Y, ElInt offset ) \
  { EL_TRY \
    ( AxpyTrapezoid \
      ( CReflect(uplo), CReflect(alpha), \
        *CReflect(X), *CReflect(Y), offset ) ) } \
  ElError ElAxpyTrapezoidDistSparse_ ## SIG \
  ( ElUpperOrLower uplo, CREFLECT(T) alpha, \
    ElConstDistSparseMatrix_ ## SIG X, ElDistSparseMatrix_ ## SIG Y, \
    ElInt offset ) \
  { EL_TRY \
    ( AxpyTrapezoid \
      ( CReflect(uplo), CReflect(alpha), \
        *CReflect(X), *CReflect(Y), offset ) ) } \
  /* Horizontal concatenation */ \
  ElError ElHCat_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG B, \
    ElMatrix_ ## SIG C ) \
  { EL_TRY( HCat( *CReflect(A), *CReflect(B), *CReflect(C) ) ) } \
  ElError ElHCatDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG B, \
    ElDistMatrix_ ## SIG C ) \
  { EL_TRY( HCat( *CReflect(A), *CReflect(B), *CReflect(C) ) ) } \
  ElError ElHCatSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElConstSparseMatrix_ ## SIG B, \
    ElSparseMatrix_ ## SIG C ) \
  { EL_TRY( HCat( *CReflect(A), *CReflect(B), *CReflect(C) ) ) } \
  ElError ElHCatDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElConstDistSparseMatrix_ ## SIG B, \
    ElDistSparseMatrix_ ## SIG C ) \
  { EL_TRY( HCat( *CReflect(A), *CReflect(B), *CReflect(C) ) ) } \
  ElError ElHCatDistMultiVec_ ## SIG \
  ( ElConstDistMultiVec_ ## SIG A, ElConstDistMultiVec_ ## SIG B, \
    ElDistMultiVec_ ## SIG C ) \
  { EL_TRY( HCat( *CReflect(A), *CReflect(B), *CReflect(C) ) ) } \
  /* Vertical concatenation */ \
  ElError ElVCat_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG B, \
    ElMatrix_ ## SIG C ) \
  { EL_TRY( VCat( *CReflect(A), *CReflect(B), *CReflect(C) ) ) } \
  ElError ElVCatDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG B, \
    ElDistMatrix_ ## SIG C ) \
  { EL_TRY( VCat( *CReflect(A), *CReflect(B), *CReflect(C) ) ) } \
  ElError ElVCatSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElConstSparseMatrix_ ## SIG B, \
    ElSparseMatrix_ ## SIG C ) \
  { EL_TRY( VCat( *CReflect(A), *CReflect(B), *CReflect(C) ) ) } \
  ElError ElVCatDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElConstDistSparseMatrix_ ## SIG B, \
    ElDistSparseMatrix_ ## SIG C ) \
  { EL_TRY( VCat( *CReflect(A), *CReflect(B), *CReflect(C) ) ) } \
  ElError ElVCatDistMultiVec_ ## SIG \
  ( ElConstDistMultiVec_ ## SIG A, ElConstDistMultiVec_ ## SIG B, \
    ElDistMultiVec_ ## SIG C ) \
  { EL_TRY( VCat( *CReflect(A), *CReflect(B), *CReflect(C) ) ) } \
  /* B = A */ \
  ElError ElCopy_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElMatrix_ ## SIG B ) \
  { EL_TRY( Copy( *CReflect(A), *CReflect(B) ) ) } \
  ElError ElCopyDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG B ) \
  { EL_TRY( Copy( *CReflect(A), *CReflect(B) ) ) } \
  ElError ElCopySparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElSparseMatrix_ ## SIG B ) \
  { EL_TRY( Copy( *CReflect(A), *CReflect(B) ) ) } \
  ElError ElCopySparseToDense_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElMatrix_ ## SIG B ) \
  { EL_TRY( Copy( *CReflect(A), *CReflect(B) ) ) } \
  ElError ElCopyDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElDistSparseMatrix_ ## SIG B ) \
  { EL_TRY( Copy( *CReflect(A), *CReflect(B) ) ) } \
  ElError ElCopyDistSparseToDense_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElDistMatrix_ ## SIG B ) \
  { EL_TRY( Copy( *CReflect(A), *CReflect(B) ) ) } \
  ElError ElCopySparseMatrixFromRoot_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG ADist, ElSparseMatrix_ ## SIG A ) \
  { EL_TRY( CopyFromRoot( *CReflect(ADist), *CReflect(A) ) ) } \
  ElError ElCopySparseMatrixFromNonRoot_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG ADist, int root ) \
  { EL_TRY( CopyFromNonRoot( *CReflect(ADist), root ) ) } \
  ElError ElCopyDistMultiVec_ ## SIG \
  ( ElConstDistMultiVec_ ## SIG A, ElDistMultiVec_ ## SIG B ) \
  { EL_TRY( Copy( *CReflect(A), *CReflect(B) ) ) } \
  ElError ElCopyMultiVecFromRoot_ ## SIG \
  ( ElConstDistMultiVec_ ## SIG A, ElMatrix_ ## SIG B ) \
  { EL_TRY( CopyFromRoot( *CReflect(A), *CReflect(B) ) ) } \
  ElError ElCopyMultiVecFromNonRoot_ ## SIG \
  ( ElConstDistMultiVec_ ## SIG A, int root ) \
  { EL_TRY( CopyFromNonRoot( *CReflect(A), root ) ) } \
  /* Dot product (<A,B>=vec(A)^H vec(B)) */ \
  ElError ElDot_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG B, CREFLECT(T)* prod ) \
  { EL_TRY( *prod = CReflect(Dot(*CReflect(A),*CReflect(B))) ) } \
  ElError ElDotDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG B, \
    CREFLECT(T)* prod ) \
  { EL_TRY( *prod = CReflect(Dot(*CReflect(A),*CReflect(B))) ) } \
  ElError ElDotDistMultiVec_ ## SIG \
  ( ElConstDistMultiVec_ ## SIG A, ElConstDistMultiVec_ ## SIG B, \
    CREFLECT(T)* prod ) \
  { EL_TRY( *prod = CReflect(Dot(*CReflect(A),*CReflect(B))) ) } \
  /* Unconjugated dot product, vec(A)^T vec(B) */ \
  ElError ElDotu_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG B, CREFLECT(T)* prod ) \
  { EL_TRY( *prod = CReflect(Dotu(*CReflect(A),*CReflect(B))) ) } \
  ElError ElDotuDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG B, \
    CREFLECT(T)* prod ) \
  { EL_TRY( *prod = CReflect(Dotu(*CReflect(A),*CReflect(B))) ) } \
  ElError ElDotuDistMultiVec_ ## SIG \
  ( ElConstDistMultiVec_ ## SIG A, ElConstDistMultiVec_ ## SIG B, \
    CREFLECT(T)* prod ) \
  { EL_TRY( *prod = CReflect(Dotu(*CReflect(A),*CReflect(B))) ) } \
  /* EntrywiseFill */ \
  ElError ElEntrywiseFill_ ## SIG \
  ( ElMatrix_ ## SIG A, CREFLECT(T) (*fill)() ) \
  { try { \
      auto newFill = [&]() { return CReflect(fill()); }; \
      EntrywiseFill( *CReflect(A), function<T(void)>(newFill) ); \
    } EL_CATCH; return EL_SUCCESS; } \
  ElError ElEntrywiseFillDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, CREFLECT(T) (*fill)() ) \
  { try { \
      auto newFill = [&]() { return CReflect(fill()); }; \
      EntrywiseFill( *CReflect(A), function<T(void)>(newFill) ); \
    } EL_CATCH; return EL_SUCCESS; } \
  ElError ElEntrywiseFillDistMultiVec_ ## SIG \
  ( ElDistMultiVec_ ## SIG A, CREFLECT(T) (*fill)() ) \
  { try { \
      auto newFill = [&]() { return CReflect(fill()); }; \
      EntrywiseFill( *CReflect(A), function<T(void)>(newFill) ); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* EntrywiseMap */ \
  ElError ElEntrywiseMap_ ## SIG \
  ( ElMatrix_ ## SIG A, CREFLECT(T) (*func)(CREFLECT(T)) ) \
  { try { \
      auto newMap = [&]( T alpha ) \
        { return CReflect(func(CReflect(alpha))); }; \
      EntrywiseMap( *CReflect(A), function<T(const T&)>(newMap) ); \
    } EL_CATCH; return EL_SUCCESS; } \
  ElError ElEntrywiseMapDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, CREFLECT(T) (*func)(CREFLECT(T)) ) \
  { try { \
      auto newMap = [&]( T alpha ) \
        { return CReflect(func(CReflect(alpha))); }; \
      EntrywiseMap( *CReflect(A), function<T(const T&)>(newMap) ); \
    } EL_CATCH; return EL_SUCCESS; } \
  ElError ElEntrywiseMapSparse_ ## SIG \
  ( ElSparseMatrix_ ## SIG A, CREFLECT(T) (*func)(CREFLECT(T)) ) \
  { try { \
      auto newMap = [&]( T alpha ) \
        { return CReflect(func(CReflect(alpha))); }; \
      EntrywiseMap( *CReflect(A), function<T(const T&)>(newMap) ); \
    } EL_CATCH; return EL_SUCCESS; } \
  ElError ElEntrywiseMapDistSparse_ ## SIG \
  ( ElDistSparseMatrix_ ## SIG A, CREFLECT(T) (*func)(CREFLECT(T)) ) \
  { try { \
      auto newMap = [&]( T alpha ) \
        { return CReflect(func(CReflect(alpha))); }; \
      EntrywiseMap( *CReflect(A), function<T(const T&)>(newMap) ); \
    } EL_CATCH; return EL_SUCCESS; } \
  ElError ElEntrywiseMapDistMultiVec_ ## SIG \
  ( ElDistMultiVec_ ## SIG A, CREFLECT(T) (*func)(CREFLECT(T)) ) \
  { try { \
      auto newMap = [&]( T alpha ) \
        { return CReflect(func(CReflect(alpha))); }; \
      EntrywiseMap( *CReflect(A), function<T(const T&)>(newMap) ); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* Fill */ \
  ElError ElFill_ ## SIG \
  ( ElMatrix_ ## SIG A, CREFLECT(T) alpha ) \
  { EL_TRY( Fill( *CReflect(A), CReflect(alpha) ) ) } \
  ElError ElFillDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, CREFLECT(T) alpha ) \
  { EL_TRY( Fill( *CReflect(A), CReflect(alpha) ) ) } \
  ElError ElFillDistMultiVec_ ## SIG \
  ( ElDistMultiVec_ ## SIG A, CREFLECT(T) alpha ) \
  { EL_TRY( Fill( *CReflect(A), CReflect(alpha) ) ) } \
  ElError ElFillSparse_ ## SIG \
  ( ElSparseMatrix_ ## SIG A, CREFLECT(T) alpha ) \
  { EL_TRY( Fill( *CReflect(A), CReflect(alpha) ) ) } \
  ElError ElFillDistSparse_ ## SIG \
  ( ElDistSparseMatrix_ ## SIG A, CREFLECT(T) alpha ) \
  { EL_TRY( Fill( *CReflect(A), CReflect(alpha) ) ) } \
  /* FillDiagonal */ \
  ElError ElFillDiagonal_ ## SIG \
  ( ElMatrix_ ## SIG A, CREFLECT(T) alpha, ElInt offset ) \
  { EL_TRY( FillDiagonal( *CReflect(A), CReflect(alpha), offset ) ) } \
  ElError ElFillDiagonalDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, CREFLECT(T) alpha, ElInt offset ) \
  { EL_TRY( FillDiagonal( *CReflect(A), CReflect(alpha), offset ) ) } \
  /* Full */ \
  ElError ElFull_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElMatrix_ ## SIG B ) \
  { EL_TRY( Copy( *CReflect(A), *CReflect(B) ) ) } \
  ElError ElFullDist_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElDistMatrix_ ## SIG B ) \
  { EL_TRY( Copy( *CReflect(A), *CReflect(B) ) ) } \
  /* GetSubmatrix */ \
  ElError ElGetContigSubmatrix_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElRange_i I, ElRange_i J, \
    ElMatrix_ ## SIG ASub ) \
  { EL_TRY( GetSubmatrix \
            ( *CReflect(A), CReflect(I), CReflect(J), *CReflect(ASub) ) ) } \
  ElError ElGetContigSubmatrixDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElRange_i I, ElRange_i J, \
    ElDistMatrix_ ## SIG ASub ) \
  { EL_TRY( GetSubmatrix \
            ( *CReflect(A), CReflect(I), CReflect(J), *CReflect(ASub) ) ) } \
  ElError ElGetContigSubmatrixSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElRange_i I, ElRange_i J, \
    ElSparseMatrix_ ## SIG ASub ) \
  { EL_TRY( GetSubmatrix \
            ( *CReflect(A), CReflect(I), CReflect(J), *CReflect(ASub) ) ) } \
  ElError ElGetContigSubmatrixDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElRange_i I, ElRange_i J, \
    ElDistSparseMatrix_ ## SIG ASub ) \
  { EL_TRY( GetSubmatrix \
            ( *CReflect(A), CReflect(I), CReflect(J), *CReflect(ASub) ) ) } \
  ElError ElGetContigSubmatrixDistMultiVec_ ## SIG \
  ( ElConstDistMultiVec_ ## SIG A, ElRange_i I, ElRange_i J, \
    ElDistMultiVec_ ## SIG ASub ) \
  { EL_TRY( GetSubmatrix \
            ( *CReflect(A), CReflect(I), CReflect(J), *CReflect(ASub) ) ) } \
  /* Hadamard */ \
  ElError ElHadamard_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG B, ElMatrix_ ## SIG C ) \
  { EL_TRY( Hadamard(*CReflect(A),*CReflect(B),*CReflect(C)) ) } \
  ElError ElHadamardDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG B, \
    ElDistMatrix_ ## SIG C ) \
  { EL_TRY( Hadamard(*CReflect(A),*CReflect(B),*CReflect(C)) ) } \
  ElError ElHadamardDistMultiVec_ ## SIG \
  ( ElConstDistMultiVec_ ## SIG A, ElConstDistMultiVec_ ## SIG B, \
    ElDistMultiVec_ ## SIG C ) \
  { EL_TRY( Hadamard(*CReflect(A),*CReflect(B),*CReflect(C)) ) } \
  /* Hilbert-Schmidt inner product (same as Dot) */ \
  ElError ElHilbertSchmidt_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG B, CREFLECT(T)* prod ) \
  { EL_TRY( *prod = \
      CReflect(HilbertSchmidt(*CReflect(A),*CReflect(B))) ) } \
  ElError ElHilbertSchmidtDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElConstDistMatrix_ ## SIG B, \
    CREFLECT(T)* prod ) \
  { EL_TRY( *prod = \
      CReflect(HilbertSchmidt(*CReflect(A),*CReflect(B))) ) } \
  ElError ElHilbertSchmidtDistMultiVec_ ## SIG \
  ( ElConstDistMultiVec_ ## SIG A, ElConstDistMultiVec_ ## SIG B, \
    CREFLECT(T)* prod ) \
  { EL_TRY( *prod = CReflect(HilbertSchmidt(*CReflect(A),*CReflect(B))) ) } \
  /* Imaginary part */ \
  ElError ElImagPart_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElMatrix_ ## SIGBASE AImag ) \
  { EL_TRY( ImagPart( *CReflect(A), *CReflect(AImag) ) ) } \
  ElError ElImagPartDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElDistMatrix_ ## SIGBASE AImag ) \
  { EL_TRY( ImagPart( *CReflect(A), *CReflect(AImag) ) ) } \
  /* IndexDependentFill */ \
  ElError ElIndexDependentFill_ ## SIG \
  ( ElMatrix_ ## SIG A, CREFLECT(T) (*fill)(ElInt,ElInt) ) \
  { try { \
      auto newFill = [&]( Int i, Int j ) { return CReflect(fill(i,j)); }; \
      IndexDependentFill \
      ( *CReflect(A), function<T(Int,Int)>(newFill) ); \
    } EL_CATCH; return EL_SUCCESS; } \
  ElError ElIndexDependentFillDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, CREFLECT(T) (*fill)(ElInt,ElInt) ) \
  { try { \
      auto newFill = [&]( Int i, Int j ) { return CReflect(fill(i,j)); }; \
      IndexDependentFill \
      ( *CReflect(A), function<T(Int,Int)>(newFill) ); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* IndexDependentMap */ \
  ElError ElIndexDependentMap_ ## SIG \
  ( ElMatrix_ ## SIG A, CREFLECT(T) (*func)(ElInt,ElInt,CREFLECT(T)) ) \
  { try { \
      auto newMap = [&]( Int i, Int j, T alpha ) \
        { return CReflect(func(i,j,CReflect(alpha))); }; \
      IndexDependentMap \
      ( *CReflect(A), function<T(Int,Int,const T&)>(newMap) ); \
    } EL_CATCH; return EL_SUCCESS; } \
  ElError ElIndexDependentMapDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, CREFLECT(T) (*func)(ElInt,ElInt,CREFLECT(T)) ) \
  { try { \
      auto newMap = [&]( Int i, Int j, T alpha ) \
        { return CReflect(func(i,j,CReflect(alpha))); }; \
      IndexDependentMap \
      ( *CReflect(A), function<T(Int,Int,const T&)>(newMap) ); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* Kronecker */ \
  ElError ElKronecker_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG B, ElMatrix_ ## SIG C ) \
  { EL_TRY( Kronecker( *CReflect(A), *CReflect(B), *CReflect(C) ) ) } \
  ElError ElKroneckerDist_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElConstMatrix_ ## SIG B, ElDistMatrix_ ## SIG C ) \
  { EL_TRY( Kronecker( *CReflect(A), *CReflect(B), *CReflect(C) ) ) } \
  ElError ElKroneckerSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElConstSparseMatrix_ ## SIG B, \
    ElSparseMatrix_ ## SIG C ) \
  { EL_TRY( Kronecker( *CReflect(A), *CReflect(B), *CReflect(C) ) ) } \
  ElError ElKroneckerDistSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElConstSparseMatrix_ ## SIG B, \
    ElDistSparseMatrix_ ## SIG C ) \
  { EL_TRY( Kronecker( *CReflect(A), *CReflect(B), *CReflect(C) ) ) } \
  /* MakeSymmetric */ \
  ElError ElMakeSymmetric_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A ) \
  { EL_TRY( MakeSymmetric( CReflect(uplo), *CReflect(A) ) ) } \
  ElError ElMakeSymmetricDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A ) \
  { EL_TRY( MakeSymmetric( CReflect(uplo), *CReflect(A) ) ) } \
  ElError ElMakeSymmetricSparse_ ## SIG \
  ( ElUpperOrLower uplo, ElSparseMatrix_ ## SIG A ) \
  { EL_TRY( MakeSymmetric( CReflect(uplo), *CReflect(A) ) ) } \
  ElError ElMakeSymmetricDistSparse_ ## SIG \
  ( ElUpperOrLower uplo, ElDistSparseMatrix_ ## SIG A ) \
  { EL_TRY( MakeSymmetric( CReflect(uplo), *CReflect(A) ) ) } \
  /* MakeTrapezoidal */ \
  ElError ElMakeTrapezoidal_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A, ElInt offset ) \
  { EL_TRY( MakeTrapezoidal( CReflect(uplo), *CReflect(A), offset ) ) } \
  ElError ElMakeTrapezoidalDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A, ElInt offset ) \
  { EL_TRY( MakeTrapezoidal( CReflect(uplo), *CReflect(A), offset ) ) } \
  ElError ElMakeTrapezoidalSparse_ ## SIG \
  ( ElUpperOrLower uplo, ElSparseMatrix_ ## SIG A, ElInt offset ) \
  { EL_TRY( MakeTrapezoidal( CReflect(uplo), *CReflect(A), offset ) ) } \
  ElError ElMakeTrapezoidalDistSparse_ ## SIG \
  ( ElUpperOrLower uplo, ElDistSparseMatrix_ ## SIG A, ElInt offset ) \
  { EL_TRY( MakeTrapezoidal( CReflect(uplo), *CReflect(A), offset ) ) } \
  /* MaxAbs */ \
  /* TODO */ \
  /* MaxAbsLoc */ \
  ElError ElMaxAbsLoc_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElEntry_ ## SIGBASE *entry ) \
  { EL_TRY( *entry = CReflect(MaxAbsLoc(*CReflect(A))) ) } \
  ElError ElMaxAbsLocDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElEntry_ ## SIGBASE *entry ) \
  { EL_TRY( *entry = CReflect(MaxAbsLoc(*CReflect(A))) ) } \
  ElError ElSymmetricMaxAbsLoc_ ## SIG \
  ( ElUpperOrLower uplo, ElConstMatrix_ ## SIG A, \
    ElEntry_ ## SIGBASE *entry ) \
  { EL_TRY( *entry = \
      CReflect(SymmetricMaxAbsLoc(CReflect(uplo),*CReflect(A))) ) } \
  ElError ElSymmetricMaxAbsLocDist_ ## SIG \
  ( ElUpperOrLower uplo, ElConstDistMatrix_ ## SIG A, \
    ElEntry_ ## SIGBASE *entry ) \
  { EL_TRY( *entry = \
      CReflect(SymmetricMaxAbsLoc(CReflect(uplo),*CReflect(A))) ) } \
  ElError ElVectorMaxAbsLoc_ ## SIG \
  ( ElConstMatrix_ ## SIG x, ElValueInt_ ## SIGBASE *entry ) \
  { EL_TRY( *entry = CReflect(VectorMaxAbsLoc(*CReflect(x))) ) } \
  ElError ElVectorMaxAbsLocDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG x, ElValueInt_ ## SIGBASE *entry ) \
  { EL_TRY( *entry = CReflect(VectorMaxAbsLoc(*CReflect(x))) ) } \
  /* MinAbs */ \
  /* TODO */ \
  /* MinAbsLoc */ \
  ElError ElMinAbsLoc_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElEntry_ ## SIGBASE *entry ) \
  { EL_TRY( *entry = CReflect(MinAbsLoc(*CReflect(A))) ) } \
  ElError ElMinAbsLocDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElEntry_ ## SIGBASE *entry ) \
  { EL_TRY( *entry = CReflect(MinAbsLoc(*CReflect(A))) ) } \
  ElError ElSymmetricMinAbsLoc_ ## SIG \
  ( ElUpperOrLower uplo, ElConstMatrix_ ## SIG A, \
    ElEntry_ ## SIGBASE *entry ) \
  { EL_TRY( *entry = \
      CReflect(SymmetricMinAbsLoc(CReflect(uplo),*CReflect(A))) ) } \
  ElError ElSymmetricMinAbsLocDist_ ## SIG \
  ( ElUpperOrLower uplo, ElConstDistMatrix_ ## SIG A, \
    ElEntry_ ## SIGBASE *entry ) \
  { EL_TRY( *entry = \
      CReflect(SymmetricMinAbsLoc(CReflect(uplo),*CReflect(A))) ) } \
  ElError ElVectorMinAbsLoc_ ## SIG \
  ( ElConstMatrix_ ## SIG x, ElValueInt_ ## SIGBASE *entry ) \
  { EL_TRY( *entry = CReflect(VectorMinAbsLoc(*CReflect(x))) ) } \
  ElError ElVectorMinAbsLocDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG x, ElValueInt_ ## SIGBASE *entry ) \
  { EL_TRY( *entry = CReflect(VectorMinAbsLoc(*CReflect(x))) ) } \
  /* TODO: QuasiDiagonalScale */ \
  /* Round */ \
  ElError ElRound_ ## SIG ( ElMatrix_ ## SIG A ) \
  { EL_TRY( Round( *CReflect(A) ) ) } \
  ElError ElRoundDist_ ## SIG ( ElDistMatrix_ ## SIG A ) \
  { EL_TRY( Round( *CReflect(A) ) ) } \
  ElError ElRoundDistMultiVec_ ## SIG ( ElDistMultiVec_ ## SIG A ) \
  { EL_TRY( Round( *CReflect(A) ) ) } \
  /* Scale */ \
  ElError ElScale_ ## SIG \
  ( CREFLECT(T) alpha, ElMatrix_ ## SIG A ) \
  { EL_TRY( Scale( CReflect(alpha), *CReflect(A) ) ) } \
  ElError ElScaleDist_ ## SIG \
  ( CREFLECT(T) alpha, ElDistMatrix_ ## SIG A ) \
  { EL_TRY( Scale( CReflect(alpha), *CReflect(A) ) ) } \
  ElError ElScaleSparse_ ## SIG \
  ( CREFLECT(T) alpha, ElSparseMatrix_ ## SIG A ) \
  { EL_TRY( Scale( CReflect(alpha), *CReflect(A) ) ) } \
  ElError ElScaleDistSparse_ ## SIG \
  ( CREFLECT(T) alpha, ElDistSparseMatrix_ ## SIG A ) \
  { EL_TRY( Scale( CReflect(alpha), *CReflect(A) ) ) } \
  ElError ElScaleDistMultiVec_ ## SIG \
  ( CREFLECT(T) alpha, ElDistMultiVec_ ## SIG A ) \
  { EL_TRY( Scale( CReflect(alpha), *CReflect(A) ) ) } \
  /* Real part */ \
  ElError ElRealPart_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElMatrix_ ## SIGBASE AReal ) \
  { EL_TRY( RealPart( *CReflect(A), *CReflect(AReal) ) ) } \
  ElError ElRealPartDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElDistMatrix_ ## SIGBASE AReal ) \
  { EL_TRY( RealPart( *CReflect(A), *CReflect(AReal) ) ) } \
  /* Reshape */ \
  ElError ElReshape_ ## SIG \
  ( ElInt m, ElInt n, ElConstMatrix_ ## SIG A, ElMatrix_ ## SIG B ) \
  { EL_TRY( Reshape( m, n, *CReflect(A), *CReflect(B) ) ) } \
  ElError ElReshapeDist_ ## SIG \
  ( ElInt m, ElInt n, ElConstDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG B ) \
  { EL_TRY( Reshape( m, n, *CReflect(A), *CReflect(B) ) ) } \
  ElError ElReshapeSparse_ ## SIG \
  ( ElInt m, ElInt n, \
    ElConstSparseMatrix_ ## SIG A, ElSparseMatrix_ ## SIG B ) \
  { EL_TRY( Reshape( m, n, *CReflect(A), *CReflect(B) ) ) } \
  ElError ElReshapeDistSparse_ ## SIG \
  ( ElInt m, ElInt n, \
    ElConstDistSparseMatrix_ ## SIG A, ElDistSparseMatrix_ ## SIG B ) \
  { EL_TRY( Reshape( m, n, *CReflect(A), *CReflect(B) ) ) } \
  /* ScaleTrapezoid */ \
  ElError ElScaleTrapezoid_ ## SIG \
  ( CREFLECT(T) alpha, ElUpperOrLower uplo, ElMatrix_ ## SIG A, ElInt offset ) \
  { EL_TRY( ScaleTrapezoid \
      ( CReflect(alpha), CReflect(uplo), *CReflect(A), offset ) ) } \
  ElError ElScaleTrapezoidDist_ ## SIG \
  ( CREFLECT(T) alpha, ElUpperOrLower uplo, ElDistMatrix_ ## SIG A, \
    ElInt offset ) \
  { EL_TRY( ScaleTrapezoid \
      ( CReflect(alpha), CReflect(uplo), *CReflect(A), offset ) ) } \
  ElError ElScaleTrapezoidSparse_ ## SIG \
  ( CREFLECT(T) alpha, ElUpperOrLower uplo, ElSparseMatrix_ ## SIG A, \
    ElInt offset ) \
  { EL_TRY( ScaleTrapezoid \
      ( CReflect(alpha), CReflect(uplo), *CReflect(A), offset ) ) } \
  ElError ElScaleTrapezoidDistSparse_ ## SIG \
  ( CREFLECT(T) alpha, ElUpperOrLower uplo, ElDistSparseMatrix_ ## SIG A, \
    ElInt offset ) \
  { EL_TRY( ScaleTrapezoid \
      ( CReflect(alpha), CReflect(uplo), *CReflect(A), offset ) ) } \
  /* SetDiagonal */ \
  /* TODO */ \
  /* Shift */ \
  ElError ElShift_ ## SIG \
  ( ElMatrix_ ## SIG A, CREFLECT(T) alpha ) \
  { EL_TRY( Shift( *CReflect(A), CReflect(alpha) ) ) } \
  ElError ElShiftDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, CREFLECT(T) alpha ) \
  { EL_TRY( Shift( *CReflect(A), CReflect(alpha) ) ) } \
  ElError ElShiftDistMultiVec_ ## SIG \
  ( ElDistMultiVec_ ## SIG A, CREFLECT(T) alpha ) \
  { EL_TRY( Shift( *CReflect(A), CReflect(alpha) ) ) } \
  /* ShiftDiagonal */ \
  ElError ElShiftDiagonal_ ## SIG \
  ( ElMatrix_ ## SIG A, CREFLECT(T) alpha, ElInt offset ) \
  { EL_TRY( ShiftDiagonal( *CReflect(A), CReflect(alpha), offset ) ) } \
  ElError ElShiftDiagonalDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, CREFLECT(T) alpha, ElInt offset ) \
  { EL_TRY( ShiftDiagonal( *CReflect(A), CReflect(alpha), offset ) ) } \
  ElError ElShiftDiagonalSparse_ ## SIG \
  ( ElSparseMatrix_ ## SIG A, CREFLECT(T) alpha, ElInt offset ) \
  { EL_TRY( ShiftDiagonal( *CReflect(A), CReflect(alpha), offset ) ) } \
  ElError ElShiftDiagonalDistSparse_ ## SIG \
  ( ElDistSparseMatrix_ ## SIG A, CREFLECT(T) alpha, ElInt offset ) \
  { EL_TRY( ShiftDiagonal( *CReflect(A), CReflect(alpha), offset ) ) } \
  /* Swap */ \
  ElError ElSwap_ ## SIG \
  ( ElOrientation orientation, ElMatrix_ ## SIG X, ElMatrix_ ## SIG Y ) \
  { EL_TRY( \
      Swap( CReflect(orientation), *CReflect(X), *CReflect(Y) ) ) } \
  ElError ElSwapDist_ ## SIG \
  ( ElOrientation orientation, \
    ElDistMatrix_ ## SIG X, ElDistMatrix_ ## SIG Y ) \
  { EL_TRY( \
      Swap( CReflect(orientation), *CReflect(X), *CReflect(Y) ) ) } \
  ElError ElRowSwap_ ## SIG \
  ( ElMatrix_ ## SIG A, ElInt to, ElInt from ) \
  { EL_TRY( RowSwap( *CReflect(A), to, from ) ) } \
  ElError ElRowSwapDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElInt to, ElInt from ) \
  { EL_TRY( RowSwap( *CReflect(A), to, from ) ) } \
  ElError ElColSwap_ ## SIG \
  ( ElMatrix_ ## SIG A, ElInt to, ElInt from ) \
  { EL_TRY( ColSwap( *CReflect(A), to, from ) ) } \
  ElError ElColSwapDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElInt to, ElInt from ) \
  { EL_TRY( ColSwap( *CReflect(A), to, from ) ) } \
  ElError ElSymmetricSwap_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A, ElInt to, ElInt from ) \
  { EL_TRY( SymmetricSwap( CReflect(uplo), *CReflect(A), to, from ) ) } \
  ElError ElSymmetricSwapDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A, ElInt to, ElInt from ) \
  { EL_TRY( SymmetricSwap( CReflect(uplo), *CReflect(A), to, from ) ) } \
  /* TODO: Transform2x2 */ \
  /* B = A^T */ \
  ElError ElTranspose_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElMatrix_ ## SIG B ) \
  { EL_TRY( Transpose(*CReflect(A),*CReflect(B),false) ) } \
  ElError ElTransposeDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG B ) \
  { EL_TRY( Transpose(*CReflect(A),*CReflect(B),false) ) } \
  ElError ElTransposeSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElSparseMatrix_ ## SIG B ) \
  { EL_TRY( Transpose(*CReflect(A),*CReflect(B),false) ) } \
  ElError ElTransposeDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElDistSparseMatrix_ ## SIG B ) \
  { EL_TRY( Transpose(*CReflect(A),*CReflect(B),false) ) } \
  /* UpdateDiagonal */ \
  /* TODO */ \
  /* Zero */ \
  ElError ElZero_ ## SIG ( ElMatrix_ ## SIG A ) \
  { EL_TRY( Zero( *CReflect(A) ) ) } \
  ElError ElZeroDist_ ## SIG ( ElDistMatrix_ ## SIG A ) \
  { EL_TRY( Zero( *CReflect(A) ) ) } \
  ElError ElZeroSparse_ ## SIG ( ElSparseMatrix_ ## SIG A ) \
  { EL_TRY( Zero( *CReflect(A) ) ) } \
  ElError ElZeroDistSparse_ ## SIG ( ElDistSparseMatrix_ ## SIG A ) \
  { EL_TRY( Zero( *CReflect(A) ) ) } \
  ElError ElZeroDistMultiVec_ ## SIG ( ElDistMultiVec_ ## SIG A ) \
  { EL_TRY( Zero( *CReflect(A) ) ) }

#define C_PROTO_REALONLY(SIG,Real) \
  /* DiagonalSolve */ \
  ElError ElDiagonalSolve_ ## SIG \
  ( ElLeftOrRight side, \
    ElConstMatrix_ ## SIG d, ElMatrix_ ## SIG A ) \
  { EL_TRY( \
      DiagonalSolve( CReflect(side), NORMAL, *CReflect(d), *CReflect(A) ) ) } \
  ElError ElDiagonalSolveDist_ ## SIG \
  ( ElLeftOrRight side, \
    ElConstDistMatrix_ ## SIG d, ElDistMatrix_ ## SIG A ) \
  { EL_TRY( \
      DiagonalSolve( CReflect(side), NORMAL, *CReflect(d), *CReflect(A) ) ) } \
  ElError ElDiagonalSolveSparse_ ## SIG \
  ( ElLeftOrRight side, \
    ElConstMatrix_ ## SIG d, ElSparseMatrix_ ## SIG A ) \
  { EL_TRY( \
      DiagonalSolve( CReflect(side), NORMAL, *CReflect(d), *CReflect(A) ) ) } \
  ElError ElDiagonalSolveDistSparse_ ## SIG \
  ( ElLeftOrRight side, \
    ElConstDistMultiVec_ ## SIG d, ElDistSparseMatrix_ ## SIG A ) \
  { EL_TRY( \
      DiagonalSolve( CReflect(side), NORMAL, *CReflect(d), *CReflect(A) ) ) } \
  ElError ElDiagonalSolveDistMultiVec_ ## SIG \
  ( ElLeftOrRight side, \
    ElConstDistMultiVec_ ## SIG d, ElDistMultiVec_ ## SIG A ) \
  { EL_TRY( \
      DiagonalSolve( CReflect(side), NORMAL, *CReflect(d), *CReflect(A) ) ) }

#define C_PROTO_NOCOMPLEX(SIG,T) \
  /* DiagonalScale */ \
  ElError ElDiagonalScale_ ## SIG \
  ( ElLeftOrRight side, \
    ElConstMatrix_ ## SIG d, ElMatrix_ ## SIG A ) \
  { EL_TRY( \
      DiagonalScale( CReflect(side), NORMAL, *CReflect(d), *CReflect(A) ) ) } \
  ElError ElDiagonalScaleDist_ ## SIG \
  ( ElLeftOrRight side, \
    ElConstDistMatrix_ ## SIG d, ElDistMatrix_ ## SIG A ) \
  { EL_TRY( \
      DiagonalScale( CReflect(side), NORMAL, *CReflect(d), *CReflect(A) ) ) } \
  ElError ElDiagonalScaleSparse_ ## SIG \
  ( ElLeftOrRight side, \
    ElConstMatrix_ ## SIG d, ElSparseMatrix_ ## SIG A ) \
  { EL_TRY( \
      DiagonalScale( CReflect(side), NORMAL, *CReflect(d), *CReflect(A) ) ) } \
  ElError ElDiagonalScaleDistSparse_ ## SIG \
  ( ElLeftOrRight side, \
    ElConstDistMultiVec_ ## SIG d, ElDistSparseMatrix_ ## SIG A ) \
  { EL_TRY( \
      DiagonalScale( CReflect(side), NORMAL, *CReflect(d), *CReflect(A) ) ) } \
  ElError ElDiagonalScaleDistMultiVec_ ## SIG \
  ( ElLeftOrRight side, \
    ElConstDistMultiVec_ ## SIG d, ElDistMultiVec_ ## SIG A ) \
  { EL_TRY( \
      DiagonalScale( CReflect(side), NORMAL, *CReflect(d), *CReflect(A) ) ) } \
  /* DiagonalScaleTrapezoid */ \
  ElError ElDiagonalScaleTrapezoid_ ## SIG \
  ( ElLeftOrRight side, ElUpperOrLower uplo, \
    ElConstMatrix_ ## SIG d, ElMatrix_ ## SIG A, ElInt offset ) \
  { EL_TRY( \
      DiagonalScaleTrapezoid \
      ( CReflect(side), CReflect(uplo), NORMAL, \
        *CReflect(d), *CReflect(A), offset ) ) } \
  ElError ElDiagonalScaleTrapezoidDist_ ## SIG \
  ( ElLeftOrRight side, ElUpperOrLower uplo, \
    ElConstDistMatrix_ ## SIG d, ElDistMatrix_ ## SIG A, ElInt offset ) \
  { EL_TRY( \
      DiagonalScaleTrapezoid \
      ( CReflect(side), CReflect(uplo), NORMAL, \
        *CReflect(d), *CReflect(A), offset ) ) } \
  ElError ElDiagonalScaleTrapezoidSparse_ ## SIG \
  ( ElLeftOrRight side, ElUpperOrLower uplo, \
    ElConstMatrix_ ## SIG d, ElSparseMatrix_ ## SIG A, ElInt offset ) \
  { EL_TRY( \
      DiagonalScaleTrapezoid \
      ( CReflect(side), CReflect(uplo), NORMAL, \
        *CReflect(d), *CReflect(A), offset ) ) } \
  ElError ElDiagonalScaleTrapezoidDistSparse_ ## SIG \
  ( ElLeftOrRight side, ElUpperOrLower uplo, \
    ElConstDistMultiVec_ ## SIG d, ElDistSparseMatrix_ ## SIG A, \
    ElInt offset ) \
  { EL_TRY( \
      DiagonalScaleTrapezoid \
      ( CReflect(side), CReflect(uplo), NORMAL, \
        *CReflect(d), *CReflect(A), offset ) ) } \
  /* Max */ \
  ElError ElMax_ ## SIG \
  ( ElConstMatrix_ ## SIG A, CREFLECT(T)* value ) \
  { EL_TRY( *value = CReflect(Max(*CReflect(A))) ) } \
  ElError ElMaxDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, CREFLECT(T)* value ) \
  { EL_TRY( *value = CReflect(Max(*CReflect(A))) ) } \
  ElError ElSymmetricMax_ ## SIG \
  ( ElUpperOrLower uplo, ElConstMatrix_ ## SIG A, CREFLECT(T)* value ) \
  { EL_TRY( *value = \
    CReflect(SymmetricMax(CReflect(uplo),*CReflect(A))) ) } \
  ElError ElSymmetricMaxDist_ ## SIG \
  ( ElUpperOrLower uplo, ElConstDistMatrix_ ## SIG A, CREFLECT(T)* value ) \
  { EL_TRY( *value = \
    CReflect(SymmetricMax(CReflect(uplo),*CReflect(A))) ) } \
  /* MaxLoc */ \
  ElError ElMaxLoc_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElEntry_ ## SIG *entry ) \
  { EL_TRY( *entry = CReflect(MaxLoc(*CReflect(A))) ) } \
  ElError ElMaxLocDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElEntry_ ## SIG *entry ) \
  { EL_TRY( *entry = CReflect(MaxLoc(*CReflect(A))) ) } \
  ElError ElSymmetricMaxLoc_ ## SIG \
  ( ElUpperOrLower uplo, ElConstMatrix_ ## SIG A, ElEntry_ ## SIG *entry ) \
  { EL_TRY( *entry = \
    CReflect(SymmetricMaxLoc(CReflect(uplo),*CReflect(A))) ) } \
  ElError ElSymmetricMaxLocDist_ ## SIG \
  ( ElUpperOrLower uplo, ElConstDistMatrix_ ## SIG A, ElEntry_ ## SIG *entry ) \
  { EL_TRY( *entry = \
    CReflect(SymmetricMaxLoc(CReflect(uplo),*CReflect(A))) ) } \
  ElError ElVectorMaxLoc_ ## SIG \
  ( ElConstMatrix_ ## SIG x, ElValueInt_ ## SIG *entry ) \
  { EL_TRY( *entry = CReflect(VectorMaxLoc(*CReflect(x))) ) } \
  ElError ElVectorMaxLocDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG x, ElValueInt_ ## SIG *entry ) \
  { EL_TRY( *entry = CReflect(VectorMaxLoc(*CReflect(x))) ) } \
  /* Min */ \
  ElError ElMin_ ## SIG \
  ( ElConstMatrix_ ## SIG A, CREFLECT(T)* value ) \
  { EL_TRY( *value = CReflect(Min(*CReflect(A))) ) } \
  ElError ElMinDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, CREFLECT(T)* value ) \
  { EL_TRY( *value = CReflect(Min(*CReflect(A))) ) } \
  ElError ElSymmetricMin_ ## SIG \
  ( ElUpperOrLower uplo, ElConstMatrix_ ## SIG A, CREFLECT(T)* value ) \
  { EL_TRY( *value = \
    CReflect(SymmetricMin(CReflect(uplo),*CReflect(A))) ) } \
  ElError ElSymmetricMinDist_ ## SIG \
  ( ElUpperOrLower uplo, ElConstDistMatrix_ ## SIG A, CREFLECT(T)* value ) \
  { EL_TRY( *value = \
    CReflect(SymmetricMin(CReflect(uplo),*CReflect(A))) ) } \
  /* MinLoc */ \
  ElError ElMinLoc_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElEntry_ ## SIG *entry ) \
  { EL_TRY( *entry = CReflect(MinLoc(*CReflect(A))) ) } \
  ElError ElMinLocDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElEntry_ ## SIG *entry ) \
  { EL_TRY( *entry = CReflect(MinLoc(*CReflect(A))) ) } \
  ElError ElSymmetricMinLoc_ ## SIG \
  ( ElUpperOrLower uplo, ElConstMatrix_ ## SIG A, ElEntry_ ## SIG *entry ) \
  { EL_TRY( *entry = \
    CReflect(SymmetricMinLoc(CReflect(uplo),*CReflect(A))) ) } \
  ElError ElSymmetricMinLocDist_ ## SIG \
  ( ElUpperOrLower uplo, ElConstDistMatrix_ ## SIG A, ElEntry_ ## SIG *entry ) \
  { EL_TRY( *entry = \
    CReflect(SymmetricMinLoc(CReflect(uplo),*CReflect(A))) ) } \
  ElError ElVectorMinLoc_ ## SIG \
  ( ElConstMatrix_ ## SIG x, ElValueInt_ ## SIG *entry ) \
  { EL_TRY( *entry = CReflect(VectorMinLoc(*CReflect(x))) ) } \
  ElError ElVectorMinLocDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG x, ElValueInt_ ## SIG *entry ) \
  { EL_TRY( *entry = CReflect(VectorMinLoc(*CReflect(x))) ) }

#define C_PROTO_FIELD(SIG,SIGBASE,F) \
  /* Column norms */ \
  ElError ElColumnTwoNormsDistMultiVec_ ## SIG \
  ( ElConstDistMultiVec_ ## SIG A, ElMatrix_ ## SIGBASE norms ) \
  { EL_TRY( ColumnTwoNorms( *CReflect(A), *CReflect(norms) ) ) } \
  /* Nrm2 (same as FrobeniusNorm) */ \
  ElError ElNrm2_ ## SIG \
  ( ElConstMatrix_ ## SIG A, Base<F>* norm ) \
  { EL_TRY( *norm = Nrm2(*CReflect(A)) ) } \
  ElError ElNrm2Dist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, Base<F>* norm ) \
  { EL_TRY( *norm = Nrm2(*CReflect(A)) ) } \
  ElError ElNrm2DistMultiVec_ ## SIG \
  ( ElConstDistMultiVec_ ## SIG A, Base<F>* norm ) \
  { EL_TRY( *norm = Nrm2( *CReflect(A) ) ) } \
  /* TODO: QuasiDiagonalSolve */ \
  /* TODO: Symmetric2x2Inv */

#define C_PROTO_INT(SIG,T) \
  C_PROTO_BASE(SIG,SIG,T) \
  C_PROTO_NOCOMPLEX(SIG,T)

#define C_PROTO_REAL(SIG,Real) \
  C_PROTO_BASE(SIG,SIG,Real) \
  C_PROTO_NOCOMPLEX(SIG,Real) \
  C_PROTO_REALONLY(SIG,Real) \
  C_PROTO_FIELD(SIG,SIG,Real)

#define C_PROTO_COMPLEX(SIG,SIGBASE,T) \
  C_PROTO_BASE(SIG,SIGBASE,T) \
  C_PROTO_FIELD(SIG,SIGBASE,T) \
  /* B = A^H */ \
  ElError ElAdjoint_ ## SIG \
  ( ElConstMatrix_ ## SIG A, ElMatrix_ ## SIG B ) \
  { EL_TRY( Adjoint( *CReflect(A), *CReflect(B) ) ) } \
  ElError ElAdjointDist_ ## SIG \
  ( ElConstDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG B ) \
  { EL_TRY( Adjoint( *CReflect(A), *CReflect(B) ) ) } \
  ElError ElAdjointSparse_ ## SIG \
  ( ElConstSparseMatrix_ ## SIG A, ElSparseMatrix_ ## SIG B ) \
  { EL_TRY( Adjoint( *CReflect(A), *CReflect(B) ) ) } \
  ElError ElAdjointDistSparse_ ## SIG \
  ( ElConstDistSparseMatrix_ ## SIG A, ElDistSparseMatrix_ ## SIG B ) \
  { EL_TRY( Adjoint( *CReflect(A), *CReflect(B) ) ) } \
  /* Conjugate */ \
  ElError ElConjugate_ ## SIG( ElMatrix_ ## SIG A ) \
  { EL_TRY( Conjugate( *CReflect(A) ) ) } \
  ElError ElConjugateDist_ ## SIG( ElDistMatrix_ ## SIG A ) \
  { EL_TRY( Conjugate( *CReflect(A) ) ) } \
  /* DiagonalScale */ \
  ElError ElDiagonalScale_ ## SIG \
  ( ElLeftOrRight side, ElOrientation orientation, \
    ElConstMatrix_ ## SIG d, ElMatrix_ ## SIG A ) \
  { EL_TRY( \
      DiagonalScale \
      ( CReflect(side), CReflect(orientation), \
        *CReflect(d), *CReflect(A) ) ) } \
  ElError ElDiagonalScaleDist_ ## SIG \
  ( ElLeftOrRight side, ElOrientation orientation, \
    ElConstDistMatrix_ ## SIG d, ElDistMatrix_ ## SIG A ) \
  { EL_TRY( \
      DiagonalScale \
      ( CReflect(side), CReflect(orientation), \
        *CReflect(d), *CReflect(A) ) ) } \
  ElError ElDiagonalScaleSparse_ ## SIG \
  ( ElLeftOrRight side, ElOrientation orientation, \
    ElConstMatrix_ ## SIG d, ElSparseMatrix_ ## SIG A ) \
  { EL_TRY( \
      DiagonalScale \
      ( CReflect(side), CReflect(orientation), \
        *CReflect(d), *CReflect(A) ) ) } \
  ElError ElDiagonalScaleDistSparse_ ## SIG \
  ( ElLeftOrRight side, ElOrientation orientation, \
    ElConstDistMultiVec_ ## SIG d, ElDistSparseMatrix_ ## SIG A ) \
  { EL_TRY( \
      DiagonalScale \
      ( CReflect(side), CReflect(orientation), \
        *CReflect(d), *CReflect(A) ) ) } \
  ElError ElDiagonalScaleDistMultiVec_ ## SIG \
  ( ElLeftOrRight side, ElOrientation orientation, \
    ElConstDistMultiVec_ ## SIG d, ElDistMultiVec_ ## SIG A ) \
  { EL_TRY( \
      DiagonalScale \
      ( CReflect(side), CReflect(orientation), \
        *CReflect(d), *CReflect(A) ) ) } \
  /* DiagonalScaleTrapezoid */ \
  ElError ElDiagonalScaleTrapezoid_ ## SIG \
  ( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, \
    ElConstMatrix_ ## SIG d, ElMatrix_ ## SIG A, ElInt offset ) \
  { EL_TRY( \
      DiagonalScaleTrapezoid \
      ( CReflect(side), CReflect(uplo), CReflect(orientation), \
        *CReflect(d), *CReflect(A), offset ) ) } \
  ElError ElDiagonalScaleTrapezoidDist_ ## SIG \
  ( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, \
    ElConstDistMatrix_ ## SIG d, ElDistMatrix_ ## SIG A, ElInt offset ) \
  { EL_TRY( \
      DiagonalScaleTrapezoid \
      ( CReflect(side), CReflect(uplo), CReflect(orientation), \
        *CReflect(d), *CReflect(A), offset ) ) } \
  ElError ElDiagonalScaleTrapezoidSparse_ ## SIG \
  ( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, \
    ElConstMatrix_ ## SIG d, ElSparseMatrix_ ## SIG A, ElInt offset ) \
  { EL_TRY( \
      DiagonalScaleTrapezoid \
      ( CReflect(side), CReflect(uplo), CReflect(orientation), \
        *CReflect(d), *CReflect(A), offset ) ) } \
  ElError ElDiagonalScaleTrapezoidDistSparse_ ## SIG \
  ( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation, \
    ElConstDistMultiVec_ ## SIG d, ElDistSparseMatrix_ ## SIG A, \
    ElInt offset ) \
  { EL_TRY( \
      DiagonalScaleTrapezoid \
      ( CReflect(side), CReflect(uplo), CReflect(orientation), \
        *CReflect(d), *CReflect(A), offset ) ) } \
  /* DiagonalSolve */ \
  ElError ElDiagonalSolve_ ## SIG \
  ( ElLeftOrRight side, ElOrientation orientation, \
    ElConstMatrix_ ## SIG d, ElMatrix_ ## SIG A ) \
  { EL_TRY( \
      DiagonalSolve \
      ( CReflect(side), CReflect(orientation), \
        *CReflect(d), *CReflect(A) ) ) } \
  ElError ElDiagonalSolveDist_ ## SIG \
  ( ElLeftOrRight side, ElOrientation orientation, \
    ElConstDistMatrix_ ## SIG d, ElDistMatrix_ ## SIG A ) \
  { EL_TRY( \
      DiagonalSolve \
      ( CReflect(side), CReflect(orientation), \
        *CReflect(d), *CReflect(A) ) ) } \
  ElError ElDiagonalSolveSparse_ ## SIG \
  ( ElLeftOrRight side, ElOrientation orientation, \
    ElConstMatrix_ ## SIG d, ElSparseMatrix_ ## SIG A ) \
  { EL_TRY( \
      DiagonalSolve \
      ( CReflect(side), CReflect(orientation), \
        *CReflect(d), *CReflect(A) ) ) } \
  ElError ElDiagonalSolveDistSparse_ ## SIG \
  ( ElLeftOrRight side, ElOrientation orientation, \
    ElConstDistMultiVec_ ## SIG d, ElDistSparseMatrix_ ## SIG A ) \
  { EL_TRY( \
      DiagonalSolve \
      ( CReflect(side), CReflect(orientation), \
        *CReflect(d), *CReflect(A) ) ) } \
  ElError ElDiagonalSolveDistMultiVec_ ## SIG \
  ( ElLeftOrRight side, ElOrientation orientation, \
    ElConstDistMultiVec_ ## SIG d, ElDistMultiVec_ ## SIG A ) \
  { EL_TRY( \
      DiagonalSolve \
      ( CReflect(side), CReflect(orientation), \
        *CReflect(d), *CReflect(A) ) ) } \
  /* HermitianSwap */ \
  ElError ElHermitianSwap_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A, ElInt to, ElInt from ) \
  { EL_TRY( HermitianSwap( CReflect(uplo), *CReflect(A), to, from ) ) } \
  ElError ElHermitianSwapDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A, ElInt to, ElInt from ) \
  { EL_TRY( HermitianSwap( CReflect(uplo), *CReflect(A), to, from ) ) } \
  /* MakeHermitian */ \
  ElError ElMakeHermitian_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A ) \
  { EL_TRY( MakeHermitian( CReflect(uplo), *CReflect(A) ) ) } \
  ElError ElMakeHermitianDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A ) \
  { EL_TRY( MakeHermitian( CReflect(uplo), *CReflect(A) ) ) } \
  ElError ElMakeHermitianSparse_ ## SIG \
  ( ElUpperOrLower uplo, ElSparseMatrix_ ## SIG A ) \
  { EL_TRY( MakeHermitian( CReflect(uplo), *CReflect(A) ) ) } \
  ElError ElMakeHermitianDistSparse_ ## SIG \
  ( ElUpperOrLower uplo, ElDistSparseMatrix_ ## SIG A ) \
  { EL_TRY( MakeHermitian( CReflect(uplo), *CReflect(A) ) ) } \
  /* MakeReal */ \
  ElError ElMakeReal_ ## SIG \
  ( ElMatrix_ ## SIG A ) { EL_TRY( MakeReal( *CReflect(A) ) ) } \
  ElError ElMakeRealDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A ) { EL_TRY( MakeReal( *CReflect(A) ) ) }

#include <El/macros/CInstantiate.h>

} // extern "C"
