/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>
#include <El/matrices.hpp>
#include <El-lite.h>
#include <El/matrices.h>
using namespace El;

extern "C" {

#define C_PROTO_BASE(SIG,SIGBASE,T) \
  /* Bernoulli */ \
  ElError ElBernoulli_ ## SIG \
  ( ElMatrix_ ## SIG A, ElInt m, ElInt n, double p ) \
  { EL_TRY( Bernoulli( *CReflect(A), m, n, p ) ) } \
  ElError ElBernoulliDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElInt m, ElInt n, double p ) \
  { EL_TRY( Bernoulli( *CReflect(A), m, n, p ) ) } \
  /* Circulant */ \
  ElError ElCirculant_ ## SIG \
  ( ElMatrix_ ## SIG A, ElInt aSize, CREFLECT(T)* aBuf ) \
  { try { \
      vector<T> a( CReflect(aBuf), CReflect(aBuf)+aSize ); \
      Circulant( *CReflect(A), a ); \
    } EL_CATCH; return EL_SUCCESS; } \
  ElError ElCirculantDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElInt aSize, CREFLECT(T)* aBuf ) \
  { try { \
      vector<T> a( CReflect(aBuf), CReflect(aBuf)+aSize ); \
      Circulant( *CReflect(A), a ); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* Diagonal */ \
  ElError ElDiagonal_ ## SIG \
  ( ElMatrix_ ## SIG A, ElMatrix_ ## SIG d ) \
  { EL_TRY( Diagonal( *CReflect(A), *CReflect(d) ) ) } \
  ElError ElDiagonalDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElDistMatrix_ ## SIG d ) \
  { EL_TRY( Diagonal( *CReflect(A), *CReflect(d) ) ) } \
  ElError ElDiagonalSparse_ ## SIG \
  ( ElSparseMatrix_ ## SIG A, ElMatrix_ ## SIG d ) \
  { EL_TRY( Diagonal( *CReflect(A), *CReflect(d) ) ) } \
  ElError ElDiagonalDistSparse_ ## SIG \
  ( ElDistSparseMatrix_ ## SIG A, ElDistMultiVec_ ## SIG d ) \
  { EL_TRY( Diagonal( *CReflect(A), *CReflect(d) ) ) } \
  /* Dynamic regularization counter-example */ \
  ElError ElDynamicRegCounter_ ## SIG ( ElMatrix_ ## SIG A, ElInt n ) \
  { EL_TRY( DynamicRegCounter( *CReflect(A), n ) ) } \
  ElError ElDynamicRegCounterDist_ ## SIG ( ElDistMatrix_ ## SIG A, ElInt n ) \
  { EL_TRY( DynamicRegCounter( *CReflect(A), n ) ) } \
  ElError ElDynamicRegCounterSparse_ ## SIG \
  ( ElSparseMatrix_ ## SIG A, ElInt n ) \
  { EL_TRY( DynamicRegCounter( *CReflect(A), n ) ) } \
  ElError ElDynamicRegCounterDistSparse_ ## SIG \
  ( ElDistSparseMatrix_ ## SIG A, ElInt n ) \
  { EL_TRY( DynamicRegCounter( *CReflect(A), n ) ) } \
  /* Forsythe */ \
  ElError ElForsythe_ ## SIG \
  ( ElMatrix_ ## SIG J, ElInt n, CREFLECT(T) alpha, CREFLECT(T) lambda ) \
  { EL_TRY( \
      Forsythe( \
        *CReflect(J), n, CReflect(alpha), CReflect(lambda) ) ) } \
  ElError ElForsytheDist_ ## SIG \
  ( ElDistMatrix_ ## SIG J, ElInt n, CREFLECT(T) alpha, CREFLECT(T) lambda ) \
  { EL_TRY( \
      Forsythe( \
        *CReflect(J), n, CReflect(alpha), CReflect(lambda) ) ) } \
  /* GCD matrix */ \
  ElError ElGCDMatrix_ ## SIG ( ElMatrix_ ## SIG G, ElInt m, ElInt n ) \
  { EL_TRY( GCDMatrix( *CReflect(G), m, n ) ) } \
  ElError ElGCDMatrixDist_ ## SIG ( ElDistMatrix_ ## SIG G, ElInt m, ElInt n ) \
  { EL_TRY( GCDMatrix( *CReflect(G), m, n ) ) } \
  /* Gear */ \
  ElError ElGear_ ## SIG \
  ( ElMatrix_ ## SIG G, ElInt n, ElInt s, ElInt t ) \
  { EL_TRY( Gear( *CReflect(G), n, s, t ) ) } \
  ElError ElGearDist_ ## SIG \
  ( ElDistMatrix_ ## SIG G, ElInt n, ElInt s, ElInt t ) \
  { EL_TRY( Gear( *CReflect(G), n, s, t ) ) } \
  /* GEPP Growth */ \
  ElError ElGEPPGrowth_ ## SIG ( ElMatrix_ ## SIG A, ElInt n ) \
  { EL_TRY( GEPPGrowth( *CReflect(A), n ) ) } \
  ElError ElGEPPGrowthDist_ ## SIG ( ElDistMatrix_ ## SIG A, ElInt n ) \
  { EL_TRY( GEPPGrowth( *CReflect(A), n ) ) } \
  /* Grcar */ \
  ElError ElGrcar_ ## SIG ( ElMatrix_ ## SIG A, ElInt n, ElInt k ) \
  { EL_TRY( Grcar( *CReflect(A), n, k ) ) } \
  ElError ElGrcarDist_ ## SIG ( ElDistMatrix_ ## SIG A, ElInt n, ElInt k ) \
  { EL_TRY( Grcar( *CReflect(A), n, k ) ) } \
  /* Hankel */ \
  ElError ElHankel_ ## SIG \
  ( ElMatrix_ ## SIG A, ElInt m, ElInt n, ElInt aSize, CREFLECT(T)* aBuf ) \
  { try { \
      vector<T> a( CReflect(aBuf), CReflect(aBuf)+aSize ); \
      Hankel( *CReflect(A), m, n, a ); \
    } EL_CATCH; return EL_SUCCESS; } \
  ElError ElHankelDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElInt m, ElInt n, ElInt aSize, CREFLECT(T)* aBuf ) \
  { try { \
      vector<T> a( CReflect(aBuf), CReflect(aBuf)+aSize ); \
      Hankel( *CReflect(A), m, n, a ); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* Hanowa */ \
  ElError ElHanowa_ ## SIG \
  ( ElMatrix_ ## SIG A, ElInt n, CREFLECT(T) mu ) \
  { EL_TRY( Hanowa( *CReflect(A), n, CReflect(mu) ) ) } \
  ElError ElHanowaDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElInt n, CREFLECT(T) mu ) \
  { EL_TRY( Hanowa( *CReflect(A), n, CReflect(mu) ) ) } \
  /* Identity */ \
  ElError ElIdentity_ ## SIG ( ElMatrix_ ## SIG A, ElInt m, ElInt n ) \
  { EL_TRY( Identity( *CReflect(A), m, n ) ) } \
  ElError ElIdentityDist_ ## SIG ( ElDistMatrix_ ## SIG A, ElInt m, ElInt n ) \
  { EL_TRY( Identity( *CReflect(A), m, n ) ) } \
  ElError ElIdentitySparse_ ## SIG \
  ( ElSparseMatrix_ ## SIG A, ElInt m, ElInt n ) \
  { EL_TRY( Identity( *CReflect(A), m, n ) ) } \
  ElError ElIdentityDistSparse_ ## SIG \
  ( ElDistSparseMatrix_ ## SIG A, ElInt m, ElInt n ) \
  { EL_TRY( Identity( *CReflect(A), m, n ) ) } \
  /* Jordan */ \
  ElError ElJordan_ ## SIG \
  ( ElMatrix_ ## SIG A, ElInt n, CREFLECT(T) lambda ) \
  { EL_TRY( Jordan( *CReflect(A), n, CReflect(lambda) ) ) } \
  ElError ElJordanDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElInt n, CREFLECT(T) lambda ) \
  { EL_TRY( Jordan( *CReflect(A), n, CReflect(lambda) ) ) } \
  /* Jordan-Cholesky */ \
  ElError ElJordanCholesky_ ## SIG ( ElMatrix_ ## SIG A, ElInt n ) \
  { EL_TRY( JordanCholesky( *CReflect(A), n ) ) } \
  ElError ElJordanCholeskyDist_ ## SIG ( ElDistMatrix_ ## SIG A, ElInt n ) \
  { EL_TRY( JordanCholesky( *CReflect(A), n ) ) } \
  ElError ElJordanCholeskySparse_ ## SIG \
  ( ElSparseMatrix_ ## SIG A, ElInt n ) \
  { EL_TRY( JordanCholesky( *CReflect(A), n ) ) } \
  ElError ElJordanCholeskyDistSparse_ ## SIG \
  ( ElDistSparseMatrix_ ## SIG A, ElInt n ) \
  { EL_TRY( JordanCholesky( *CReflect(A), n ) ) } \
  /* KMS */ \
  ElError ElKMS_ ## SIG \
  ( ElMatrix_ ## SIG K, ElInt n, CREFLECT(T) rho ) \
  { EL_TRY( KMS( *CReflect(K), n, CReflect(rho) ) ) } \
  ElError ElKMSDist_ ## SIG \
  ( ElDistMatrix_ ## SIG K, ElInt n, CREFLECT(T) rho ) \
  { EL_TRY( KMS( *CReflect(K), n, CReflect(rho) ) ) } \
  /* Lauchli */ \
  ElError ElLauchli_ ## SIG \
  ( ElMatrix_ ## SIG A, ElInt n, CREFLECT(T) mu ) \
  { EL_TRY( Lauchli( *CReflect(A), n, CReflect(mu) ) ) } \
  ElError ElLauchliDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElInt n, CREFLECT(T) mu ) \
  { EL_TRY( Lauchli( *CReflect(A), n, CReflect(mu) ) ) } \
  /* MinIJ */ \
  ElError ElMinIJ_ ## SIG ( ElMatrix_ ## SIG A, ElInt n ) \
  { EL_TRY( MinIJ( *CReflect(A), n ) ) } \
  ElError ElMinIJDist_ ## SIG ( ElDistMatrix_ ## SIG A, ElInt n ) \
  { EL_TRY( MinIJ( *CReflect(A), n ) ) } \
  /* Ones */ \
  ElError ElOnes_ ## SIG ( ElMatrix_ ## SIG A, ElInt m, ElInt n ) \
  { EL_TRY( Ones( *CReflect(A), m, n ) ) } \
  ElError ElOnesDist_ ## SIG ( ElDistMatrix_ ## SIG A, ElInt m, ElInt n ) \
  { EL_TRY( Ones( *CReflect(A), m, n ) ) } \
  ElError ElOnesDistMultiVec_ ## SIG \
  ( ElDistMultiVec_ ## SIG A, ElInt m, ElInt n ) \
  { EL_TRY( Ones( *CReflect(A), m, n ) ) } \
  ElError ElOnesSparse_ ## SIG \
  ( ElSparseMatrix_ ## SIG A, ElInt m, ElInt n ) \
  { EL_TRY( Ones( *CReflect(A), m, n ) ) } \
  ElError ElOnesDistSparse_ ## SIG \
  ( ElDistSparseMatrix_ ## SIG A, ElInt m, ElInt n ) \
  { EL_TRY( Ones( *CReflect(A), m, n ) ) } \
  /* 1-2-1 */ \
  ElError ElOneTwoOne_ ## SIG ( ElMatrix_ ## SIG A, ElInt n ) \
  { EL_TRY( OneTwoOne( *CReflect(A), n ) ) } \
  ElError ElOneTwoOneDist_ ## SIG ( ElDistMatrix_ ## SIG A, ElInt n ) \
  { EL_TRY( OneTwoOne( *CReflect(A), n ) ) } \
  /* Redheffer */ \
  ElError ElRedheffer_ ## SIG ( ElMatrix_ ## SIG A, ElInt n ) \
  { EL_TRY( Redheffer( *CReflect(A), n ) ) } \
  ElError ElRedhefferDist_ ## SIG ( ElDistMatrix_ ## SIG A, ElInt n ) \
  { EL_TRY( Redheffer( *CReflect(A), n ) ) } \
  /* Rademacher */ \
  ElError ElRademacher_ ## SIG ( ElMatrix_ ## SIG A, ElInt m, ElInt n ) \
  { EL_TRY( Rademacher( *CReflect(A), m, n ) ) } \
  ElError ElRademacherDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElInt m, ElInt n ) \
  { EL_TRY( Rademacher( *CReflect(A), m, n ) ) } \
  /* Three-valued */ \
  ElError ElThreeValued_ ## SIG \
  ( ElMatrix_ ## SIG A, ElInt m, ElInt n, double p ) \
  { EL_TRY( ThreeValued( *CReflect(A), m, n, p ) ) } \
  ElError ElThreeValuedDist_ ## SIG  \
  ( ElDistMatrix_ ## SIG A, ElInt m, ElInt n, double p ) \
  { EL_TRY( ThreeValued( *CReflect(A), m, n, p ) ) } \
  /* Toeplitz */ \
  ElError ElToeplitz_ ## SIG \
  ( ElMatrix_ ## SIG A, ElInt m, ElInt n, ElInt aSize, CREFLECT(T)* aBuf ) \
  { try { \
      vector<T> a( CReflect(aBuf), CReflect(aBuf)+aSize ); \
      Toeplitz( *CReflect(A), m, n, a ); \
    } EL_CATCH; return EL_SUCCESS; } \
  ElError ElToeplitzDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElInt m, ElInt n, ElInt aSize, CREFLECT(T)* aBuf ) \
  { try { \
      vector<T> a( CReflect(aBuf), CReflect(aBuf)+aSize ); \
      Toeplitz( *CReflect(A), m, n, a ); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* TriW */ \
  ElError ElTriW_ ## SIG \
  ( ElMatrix_ ## SIG A, ElInt n, CREFLECT(T) alpha, ElInt k ) \
  { EL_TRY( TriW( *CReflect(A), n, CReflect(alpha), k ) ) } \
  ElError ElTriWDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElInt n, CREFLECT(T) alpha, ElInt k ) \
  { EL_TRY( TriW( *CReflect(A), n, CReflect(alpha), k ) ) } \
  /* Uniform */ \
  ElError ElUniform_ ## SIG \
  ( ElMatrix_ ## SIG A, ElInt m, ElInt n, \
    CREFLECT(T) center, Base<T> radius ) \
  { EL_TRY( Uniform( *CReflect(A), m, n, CReflect(center), radius ) ) } \
  ElError ElUniformDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElInt m, ElInt n, \
    CREFLECT(T) center, Base<T> radius ) \
  { EL_TRY( Uniform( *CReflect(A), m, n, CReflect(center), radius ) ) } \
  ElError ElUniformDistMultiVec_ ## SIG \
  ( ElDistMultiVec_ ## SIG A, ElInt m, ElInt n, \
    CREFLECT(T) center, Base<T> radius ) \
  { EL_TRY( Uniform( *CReflect(A), m, n, CReflect(center), radius ) ) } \
  /* Walsh */ \
  ElError ElWalsh_ ## SIG ( ElMatrix_ ## SIG A, ElInt k, bool binary ) \
  { EL_TRY( Walsh( *CReflect(A), k, binary ) ) } \
  ElError ElWalshDist_ ## SIG ( ElDistMatrix_ ## SIG A, ElInt k, bool binary ) \
  { EL_TRY( Walsh( *CReflect(A), k, binary ) ) } \
  /* Wilkinson */ \
  ElError ElWilkinson_ ## SIG ( ElMatrix_ ## SIG A, ElInt k ) \
  { EL_TRY( Wilkinson( *CReflect(A), k ) ) } \
  ElError ElWilkinsonDist_ ## SIG ( ElDistMatrix_ ## SIG A, ElInt k ) \
  { EL_TRY( Wilkinson( *CReflect(A), k ) ) } \
  /* Zeros */ \
  ElError ElZeros_ ## SIG ( ElMatrix_ ## SIG A, ElInt m, ElInt n ) \
  { EL_TRY( Zeros( *CReflect(A), m, n ) ) } \
  ElError ElZerosDist_ ## SIG ( ElDistMatrix_ ## SIG A, ElInt m, ElInt n ) \
  { EL_TRY( Zeros( *CReflect(A), m, n ) ) } \
  ElError ElZerosSparse_ ## SIG ( ElSparseMatrix_ ## SIG A, ElInt m, ElInt n ) \
  { EL_TRY( Zeros( *CReflect(A), m, n ) ) } \
  ElError ElZerosDistSparse_ ## SIG \
  ( ElDistSparseMatrix_ ## SIG A, ElInt m, ElInt n ) \
  { EL_TRY( Zeros( *CReflect(A), m, n ) ) } \
  ElError ElZerosDistMultiVec_ ## SIG \
  ( ElDistMultiVec_ ## SIG A, ElInt m, ElInt n ) \
  { EL_TRY( Zeros( *CReflect(A), m, n ) ) }

#define C_PROTO_FIELD(SIG,SIGBASE,T) \
  /* Cauchy */ \
  ElError ElCauchy_ ## SIG \
  ( ElMatrix_ ## SIG A, ElInt xSize, CREFLECT(T)* xBuf, \
                        ElInt ySize, CREFLECT(T)* yBuf ) \
  { try { \
      vector<T> x( CReflect(xBuf), CReflect(xBuf)+xSize ); \
      vector<T> y( CReflect(yBuf), CReflect(yBuf)+ySize ); \
      Cauchy( *CReflect(A), x, y ); \
    } EL_CATCH; return EL_SUCCESS; } \
  ElError ElCauchyDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElInt xSize, CREFLECT(T)* xBuf, \
                            ElInt ySize, CREFLECT(T)* yBuf ) \
  { try { \
      vector<T> x( CReflect(xBuf), CReflect(xBuf)+xSize ); \
      vector<T> y( CReflect(yBuf), CReflect(yBuf)+ySize ); \
      Cauchy( *CReflect(A), x, y ); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* Cauchy-like */ \
  ElError ElCauchyLike_ ## SIG \
  ( ElMatrix_ ## SIG A, ElInt rSize, CREFLECT(T)* rBuf, \
                        ElInt sSize, CREFLECT(T)* sBuf, \
                        ElInt xSize, CREFLECT(T)* xBuf, \
                        ElInt ySize, CREFLECT(T)* yBuf ) \
  { try { \
      vector<T> r( CReflect(rBuf), CReflect(rBuf)+rSize ); \
      vector<T> s( CReflect(sBuf), CReflect(sBuf)+sSize ); \
      vector<T> x( CReflect(xBuf), CReflect(xBuf)+xSize ); \
      vector<T> y( CReflect(yBuf), CReflect(yBuf)+ySize ); \
      CauchyLike( *CReflect(A), r, s, x, y ); \
    } EL_CATCH; return EL_SUCCESS; } \
  ElError ElCauchyLikeDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElInt rSize, CREFLECT(T)* rBuf, \
                            ElInt sSize, CREFLECT(T)* sBuf, \
                            ElInt xSize, CREFLECT(T)* xBuf, \
                            ElInt ySize, CREFLECT(T)* yBuf ) \
  { try { \
      vector<T> r( CReflect(rBuf), CReflect(rBuf)+rSize ); \
      vector<T> s( CReflect(sBuf), CReflect(sBuf)+sSize ); \
      vector<T> x( CReflect(xBuf), CReflect(xBuf)+xSize ); \
      vector<T> y( CReflect(yBuf), CReflect(yBuf)+ySize ); \
      CauchyLike( *CReflect(A), r, s, x, y ); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* Demmel */ \
  ElError ElDemmel_ ## SIG ( ElMatrix_ ## SIG A, ElInt n ) \
  { EL_TRY( Demmel( *CReflect(A), n ) ) } \
  ElError ElDemmelDist_ ## SIG ( ElDistMatrix_ ## SIG A, ElInt n ) \
  { EL_TRY( Demmel( *CReflect(A), n ) ) } \
  /* Druinsky-Toledo */ \
  ElError ElDruinskyToledo_ ## SIG ( ElMatrix_ ## SIG A, ElInt n ) \
  { EL_TRY( DruinskyToledo( *CReflect(A), n ) ) } \
  ElError ElDruinskyToledoDist_ ## SIG ( ElDistMatrix_ ## SIG A, ElInt n ) \
  { EL_TRY( DruinskyToledo( *CReflect(A), n ) ) } \
  /* Ehrenfest */ \
  ElError ElEhrenfest_ ## SIG ( ElMatrix_ ## SIG P, ElInt n ) \
  { EL_TRY( Ehrenfest( *CReflect(P), n ) ) } \
  ElError ElEhrenfestDist_ ## SIG ( ElDistMatrix_ ## SIG P, ElInt n ) \
  { EL_TRY( Ehrenfest( *CReflect(P), n ) ) } \
  ElError ElEhrenfestStationary_ ## SIG \
  ( ElMatrix_ ## SIG PInf, ElInt n ) \
  { EL_TRY( EhrenfestStationary( *CReflect(PInf), n ) ) } \
  ElError ElEhrenfestStationaryDist_ ## SIG \
  ( ElDistMatrix_ ## SIG PInf, ElInt n ) \
  { EL_TRY( EhrenfestStationary( *CReflect(PInf), n ) ) } \
  ElError ElEhrenfestDecay_ ## SIG ( ElMatrix_ ## SIG A, ElInt n ) \
  { EL_TRY( EhrenfestDecay( *CReflect(A), n ) ) } \
  ElError ElEhrenfestDecayDist_ ## SIG ( ElDistMatrix_ ## SIG A, ElInt n ) \
  { EL_TRY( EhrenfestDecay( *CReflect(A), n ) ) } \
  /* ExtendedKahan */ \
  ElError ElExtendedKahan_ ## SIG \
  ( ElMatrix_ ## SIG A, ElInt k, Base<T> phi, Base<T> mu ) \
  { EL_TRY( ExtendedKahan( *CReflect(A), k, phi, mu ) ) } \
  ElError ElExtendedKahanDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElInt k, Base<T> phi, Base<T> mu ) \
  { EL_TRY( ExtendedKahan( *CReflect(A), k, phi, mu ) ) } \
  /* Fiedler */ \
  ElError ElFiedler_ ## SIG \
  ( ElMatrix_ ## SIG A, ElInt cSize, CREFLECT(T)* cBuf ) \
  { try { \
      vector<T> c( CReflect(cBuf), CReflect(cBuf)+cSize ); \
      Fiedler( *CReflect(A), c ); \
    } EL_CATCH; return EL_SUCCESS; } \
  ElError ElFiedlerDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElInt cSize, CREFLECT(T)* cBuf ) \
  { try { \
      vector<T> c( CReflect(cBuf), CReflect(cBuf)+cSize ); \
      Fiedler( *CReflect(A), c ); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* Gaussian */ \
  ElError ElGaussian_ ## SIG \
  ( ElMatrix_ ## SIG A, ElInt m, ElInt n, \
    CREFLECT(T) mean, Base<T> stddev ) \
  { EL_TRY( Gaussian( *CReflect(A), m, n, CReflect(mean), stddev ) ) } \
  ElError ElGaussianDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElInt m, ElInt n, \
    CREFLECT(T) mean, Base<T> stddev ) \
  { EL_TRY( Gaussian( *CReflect(A), m, n, CReflect(mean), stddev ) ) } \
  ElError ElGaussianDistMultiVec_ ## SIG \
  ( ElDistMultiVec_ ## SIG A, ElInt m, ElInt n, \
    CREFLECT(T) mean, Base<T> stddev ) \
  { EL_TRY( Gaussian( *CReflect(A), m, n, CReflect(mean), stddev ) ) } \
  /* Golub Klema Stewart */ \
  ElError ElGKS_ ## SIG ( ElMatrix_ ## SIG A, ElInt n ) \
  { EL_TRY( GKS( *CReflect(A), n ) ) } \
  ElError ElGKSDist_ ## SIG ( ElDistMatrix_ ## SIG A, ElInt n ) \
  { EL_TRY( GKS( *CReflect(A), n ) ) } \
  /* Haar */ \
  ElError ElHaar_ ## SIG ( ElMatrix_ ## SIG A, ElInt n ) \
  { EL_TRY( Haar( *CReflect(A), n ) ) } \
  ElError ElHaarDist_ ## SIG ( ElDistMatrix_ ## SIG A, ElInt n ) \
  { EL_TRY( Haar( *CReflect(A), n ) ) } \
  ElError ElImplicitHaar_ ## SIG \
  ( ElMatrix_ ## SIG A, \
    ElMatrix_ ## SIG t, ElMatrix_ ## SIGBASE d, ElInt n ) \
  { EL_TRY( \
      ImplicitHaar( *CReflect(A), *CReflect(t), *CReflect(d), n ) ) } \
  ElError ElImplicitHaarDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, \
    ElDistMatrix_ ## SIG t, ElDistMatrix_ ## SIGBASE d, ElInt n ) \
  { EL_TRY( \
      ImplicitHaar( *CReflect(A), *CReflect(t), *CReflect(d), n ) ) } \
  /* Hatano-Nelson */ \
  ElError ElHatanoNelson_ ## SIG \
  ( ElMatrix_ ## SIG A, ElInt n, CREFLECT(T) center, Base<T> radius, \
    CREFLECT(T) g, bool periodic ) \
  { EL_TRY( \
      HatanoNelson( \
        *CReflect(A), n, CReflect(center), radius, CReflect(g), \
        periodic ) ) } \
  ElError ElHatanoNelsonDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElInt n, CREFLECT(T) center, Base<T> radius, \
    CREFLECT(T) g, bool periodic ) \
  { EL_TRY( \
      HatanoNelson( \
        *CReflect(A), n, CReflect(center), radius, CReflect(g), \
        periodic ) ) } \
  /* Helmholtz */ \
  ElError ElHelmholtz1D_ ## SIG \
  ( ElMatrix_ ## SIG H, ElInt nx, CREFLECT(T) shift ) \
  { EL_TRY( Helmholtz( *CReflect(H), nx, CReflect(shift) ) ) } \
  ElError ElHelmholtz1DDist_ ## SIG \
  ( ElDistMatrix_ ## SIG H, ElInt nx, CREFLECT(T) shift ) \
  { EL_TRY( Helmholtz( *CReflect(H), nx, CReflect(shift) ) ) } \
  ElError ElHelmholtz1DSparse_ ## SIG \
  ( ElSparseMatrix_ ## SIG H, ElInt nx, CREFLECT(T) shift ) \
  { EL_TRY( Helmholtz( *CReflect(H), nx, CReflect(shift) ) ) } \
  ElError ElHelmholtz1DDistSparse_ ## SIG \
  ( ElDistSparseMatrix_ ## SIG H, ElInt nx, CREFLECT(T) shift ) \
  { EL_TRY( Helmholtz( *CReflect(H), nx, CReflect(shift) ) ) } \
  ElError ElHelmholtz2D_ ## SIG \
  ( ElMatrix_ ## SIG H, ElInt nx, ElInt ny, CREFLECT(T) shift ) \
  { EL_TRY( Helmholtz( *CReflect(H), nx, ny, CReflect(shift) ) ) } \
  ElError ElHelmholtz2DDist_ ## SIG \
  ( ElDistMatrix_ ## SIG H, ElInt nx, ElInt ny, CREFLECT(T) shift ) \
  { EL_TRY( Helmholtz( *CReflect(H), nx, ny, CReflect(shift) ) ) } \
  ElError ElHelmholtz2DSparse_ ## SIG \
  ( ElSparseMatrix_ ## SIG H, ElInt nx, ElInt ny, CREFLECT(T) shift ) \
  { EL_TRY( Helmholtz( *CReflect(H), nx, ny, CReflect(shift) ) ) } \
  ElError ElHelmholtz2DDistSparse_ ## SIG \
  ( ElDistSparseMatrix_ ## SIG H, ElInt nx, ElInt ny, CREFLECT(T) shift ) \
  { EL_TRY( Helmholtz( *CReflect(H), nx, ny, CReflect(shift) ) ) } \
  ElError ElHelmholtz3D_ ## SIG \
  ( ElMatrix_ ## SIG H, ElInt nx, ElInt ny, ElInt nz, CREFLECT(T) shift ) \
  { EL_TRY( Helmholtz( *CReflect(H), nx, ny, nz, CReflect(shift) ) ) } \
  ElError ElHelmholtz3DDist_ ## SIG \
  ( ElDistMatrix_ ## SIG H, ElInt nx, ElInt ny, ElInt nz, CREFLECT(T) shift ) \
  { EL_TRY( Helmholtz( *CReflect(H), nx, ny, nz, CReflect(shift) ) ) } \
  ElError ElHelmholtz3DSparse_ ## SIG \
  ( ElSparseMatrix_ ## SIG H, ElInt nx, ElInt ny, ElInt nz, \
    CREFLECT(T) shift ) \
  { EL_TRY( Helmholtz( *CReflect(H), nx, ny, nz, CReflect(shift) ) ) } \
  ElError ElHelmholtz3DDistSparse_ ## SIG \
  ( ElDistSparseMatrix_ ## SIG H, ElInt nx, ElInt ny, ElInt nz, \
    CREFLECT(T) shift ) \
  { EL_TRY( Helmholtz( *CReflect(H), nx, ny, nz, CReflect(shift) ) ) } \
  /* Hermitian uniform spectrum */ \
  ElError ElHermitianUniformSpectrum_ ## SIG \
  ( ElMatrix_ ## SIG A, ElInt n, Base<T> lower, Base<T> upper ) \
  { EL_TRY( HermitianUniformSpectrum( *CReflect(A), n, lower, upper ) ) } \
  ElError ElHermitianUniformSpectrumDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElInt n, Base<T> lower, Base<T> upper ) \
  { EL_TRY( HermitianUniformSpectrum( *CReflect(A), n, lower, upper ) ) } \
  /* Hilbert */ \
  ElError ElHilbert_ ## SIG ( ElMatrix_ ## SIG A, ElInt n ) \
  { EL_TRY( Hilbert( *CReflect(A), n ) ) } \
  ElError ElHilbertDist_ ## SIG ( ElDistMatrix_ ## SIG A, ElInt n ) \
  { EL_TRY( Hilbert( *CReflect(A), n ) ) } \
  /* Kahan */ \
  ElError ElKahan_ ## SIG \
  ( ElMatrix_ ## SIG A, ElInt n, CREFLECT(T) phi ) \
  { EL_TRY( Kahan( *CReflect(A), n, CReflect(phi) ) ) } \
  ElError ElKahanDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElInt n, CREFLECT(T) phi ) \
  { EL_TRY( Kahan( *CReflect(A), n, CReflect(phi) ) ) } \
  /* Laplacian */ \
  ElError ElLaplacian1D_ ## SIG ( ElMatrix_ ## SIG L, ElInt nx ) \
  { EL_TRY( Laplacian( *CReflect(L), nx ) ) } \
  ElError ElLaplacian1DDist_ ## SIG ( ElDistMatrix_ ## SIG L, ElInt nx ) \
  { EL_TRY( Laplacian( *CReflect(L), nx ) ) } \
  ElError ElLaplacian1DSparse_ ## SIG ( ElSparseMatrix_ ## SIG L, ElInt nx ) \
  { EL_TRY( Laplacian( *CReflect(L), nx ) ) } \
  ElError ElLaplacian1DDistSparse_ ## SIG \
  ( ElDistSparseMatrix_ ## SIG L, ElInt nx ) \
  { EL_TRY( Laplacian( *CReflect(L), nx ) ) } \
  ElError ElLaplacian2D_ ## SIG \
  ( ElMatrix_ ## SIG L, ElInt nx, ElInt ny ) \
  { EL_TRY( Laplacian( *CReflect(L), nx, ny ) ) } \
  ElError ElLaplacian2DDist_ ## SIG \
  ( ElDistMatrix_ ## SIG L, ElInt nx, ElInt ny ) \
  { EL_TRY( Laplacian( *CReflect(L), nx, ny ) ) } \
  ElError ElLaplacian2DSparse_ ## SIG \
  ( ElSparseMatrix_ ## SIG L, ElInt nx, ElInt ny ) \
  { EL_TRY( Laplacian( *CReflect(L), nx, ny ) ) } \
  ElError ElLaplacian2DDistSparse_ ## SIG \
  ( ElDistSparseMatrix_ ## SIG L, ElInt nx, ElInt ny ) \
  { EL_TRY( Laplacian( *CReflect(L), nx, ny ) ) } \
  ElError ElLaplacian3D_ ## SIG \
  ( ElMatrix_ ## SIG L, ElInt nx, ElInt ny, ElInt nz ) \
  { EL_TRY( Laplacian( *CReflect(L), nx, ny, nz ) ) } \
  ElError ElLaplacian3DDist_ ## SIG \
  ( ElDistMatrix_ ## SIG L, ElInt nx, ElInt ny, ElInt nz ) \
  { EL_TRY( Laplacian( *CReflect(L), nx, ny, nz ) ) } \
  ElError ElLaplacian3DSparse_ ## SIG \
  ( ElSparseMatrix_ ## SIG L, ElInt nx, ElInt ny, ElInt nz ) \
  { EL_TRY( Laplacian( *CReflect(L), nx, ny, nz ) ) } \
  ElError ElLaplacian3DDistSparse_ ## SIG \
  ( ElDistSparseMatrix_ ## SIG L, ElInt nx, ElInt ny, ElInt nz ) \
  { EL_TRY( Laplacian( *CReflect(L), nx, ny, nz ) ) } \
  /* Legendre */ \
  ElError ElLegendre_ ## SIG ( ElMatrix_ ## SIG A, ElInt n ) \
  { EL_TRY( Legendre( *CReflect(A), n ) ) } \
  ElError ElLegendreDist_ ## SIG ( ElDistMatrix_ ## SIG A, ElInt n ) \
  { EL_TRY( Legendre( *CReflect(A), n ) ) } \
  /* Lehmer */ \
  ElError ElLehmer_ ## SIG ( ElMatrix_ ## SIG L, ElInt n ) \
  { EL_TRY( Lehmer( *CReflect(L), n ) ) } \
  ElError ElLehmerDist_ ## SIG ( ElDistMatrix_ ## SIG L, ElInt n ) \
  { EL_TRY( Lehmer( *CReflect(L), n ) ) } \
  /* Lotkin */ \
  ElError ElLotkin_ ## SIG ( ElMatrix_ ## SIG A, ElInt n ) \
  { EL_TRY( Lotkin( *CReflect(A), n ) ) } \
  ElError ElLotkinDist_ ## SIG ( ElDistMatrix_ ## SIG A, ElInt n ) \
  { EL_TRY( Lotkin( *CReflect(A), n ) ) }  \
  /* Parter */ \
  ElError ElParter_ ## SIG ( ElMatrix_ ## SIG A, ElInt n ) \
  { EL_TRY( Parter( *CReflect(A), n ) ) } \
  ElError ElParterDist_ ## SIG ( ElDistMatrix_ ## SIG A, ElInt n ) \
  { EL_TRY( Parter( *CReflect(A), n ) ) } \
  /* Pei */ \
  ElError ElPei_ ## SIG \
  ( ElMatrix_ ## SIG A, ElInt n, CREFLECT(T) alpha ) \
  { EL_TRY( Pei( *CReflect(A), n, CReflect(alpha) ) ) } \
  ElError ElPeiDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElInt n, CREFLECT(T) alpha ) \
  { EL_TRY( Pei( *CReflect(A), n, CReflect(alpha) ) ) } \
  /* Riffle */ \
  ElError ElRiffle_ ## SIG ( ElMatrix_ ## SIG P, ElInt n ) \
  { EL_TRY( Riffle( *CReflect(P), n ) ) } \
  ElError ElRiffleDist_ ## SIG ( ElDistMatrix_ ## SIG P, ElInt n ) \
  { EL_TRY( Riffle( *CReflect(P), n ) ) } \
  ElError ElRiffleStationary_ ## SIG \
  ( ElMatrix_ ## SIG PInf, ElInt n ) \
  { EL_TRY( RiffleStationary( *CReflect(PInf), n ) ) } \
  ElError ElRiffleStationaryDist_ ## SIG \
  ( ElDistMatrix_ ## SIG PInf, ElInt n ) \
  { EL_TRY( RiffleStationary( *CReflect(PInf), n ) ) } \
  ElError ElRiffleDecay_ ## SIG ( ElMatrix_ ## SIG A, ElInt n ) \
  { EL_TRY( RiffleDecay( *CReflect(A), n ) ) } \
  ElError ElRiffleDecayDist_ ## SIG ( ElDistMatrix_ ## SIG A, ElInt n ) \
  { EL_TRY( RiffleDecay( *CReflect(A), n ) ) } \
  /* Ris */ \
  ElError ElRis_ ## SIG ( ElMatrix_ ## SIG A, ElInt n ) \
  { EL_TRY( Ris( *CReflect(A), n ) ) } \
  ElError ElRisDist_ ## SIG ( ElDistMatrix_ ## SIG A, ElInt n ) \
  { EL_TRY( Ris( *CReflect(A), n ) ) } \
  /* Triangle */ \
  ElError ElTriangle_ ## SIG ( ElMatrix_ ## SIG A, ElInt n ) \
  { EL_TRY( Triangle( *CReflect(A), n ) ) } \
  ElError ElTriangleDist_ ## SIG ( ElDistMatrix_ ## SIG A, ElInt n ) \
  { EL_TRY( Triangle( *CReflect(A), n ) ) } \
  /* Wigner */ \
  ElError ElWigner_ ## SIG \
  ( ElMatrix_ ## SIG A, ElInt n, CREFLECT(T) mean, Base<T> stddev ) \
  { EL_TRY( Wigner( *CReflect(A), n, CReflect(mean), stddev ) ) } \
  ElError ElWignerDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElInt n, CREFLECT(T) mean, Base<T> stddev ) \
  { EL_TRY( Wigner( *CReflect(A), n, CReflect(mean), stddev ) ) }

#define C_PROTO_INT(SIG,T) \
  C_PROTO_BASE(SIG,SIG,T)

#define C_PROTO_REAL(SIG,T) \
  C_PROTO_BASE(SIG,SIG,T) \
  C_PROTO_FIELD(SIG,SIG,T)

#define C_PROTO_COMPLEX(SIG,SIGBASE,T) \
  C_PROTO_BASE(SIG,SIGBASE,T) \
  C_PROTO_FIELD(SIG,SIGBASE,T) \
  /* Bull's head */ \
  ElError ElBullsHead_ ## SIG ( ElMatrix_ ## SIG A, ElInt n ) \
  { EL_TRY( BullsHead( *CReflect(A), n ) ) } \
  ElError ElBullsHeadDist_ ## SIG ( ElDistMatrix_ ## SIG A, ElInt n ) \
  { EL_TRY( BullsHead( *CReflect(A), n ) ) } \
  /* Egorov */ \
  ElError ElEgorov_ ## SIG \
  ( ElMatrix_ ## SIG A, Base<T> (*phase)(ElInt,ElInt), ElInt n ) \
  { try { \
      function<Base<T>(Int,Int)> phaseFunc(phase); \
      Egorov( *CReflect(A), phaseFunc, n ); \
    } EL_CATCH; return EL_SUCCESS; } \
  ElError ElEgorovDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, Base<T> (*phase)(ElInt,ElInt), ElInt n ) \
  { try { \
      function<Base<T>(Int,Int)> phaseFunc(phase); \
      Egorov( *CReflect(A), phaseFunc, n ); \
    } EL_CATCH; return EL_SUCCESS; } \
  /* Fox-Li */ \
  ElError ElFoxLi_ ## SIG \
  ( ElMatrix_ ## SIG A, ElInt n, Base<T> omega ) \
  { EL_TRY( FoxLi( *CReflect(A), n, omega ) ) } \
  ElError ElFoxLiDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElInt n, Base<T> omega ) \
  { EL_TRY( FoxLi( *CReflect(A), n, omega ) ) } \
  /* Fourier */ \
  ElError ElFourier_ ## SIG ( ElMatrix_ ## SIG A, ElInt n ) \
  { EL_TRY( Fourier( *CReflect(A), n ) ) } \
  ElError ElFourierDist_ ## SIG ( ElDistMatrix_ ## SIG A, ElInt n ) \
  { EL_TRY( Fourier( *CReflect(A), n ) ) } \
  /* Helmholtz with PML */ \
  ElError ElHelmholtzPML1D_ ## SIG \
  ( ElMatrix_ ## SIG H, ElInt nx, CREFLECT(T) omega, \
    ElInt numPmlPoints, Base<T> sigma, Base<T> pmlExp ) \
  { EL_TRY( \
      HelmholtzPML( \
        *CReflect(H), nx, CReflect(omega), \
        numPmlPoints, sigma, pmlExp ) ) } \
  ElError ElHelmholtzPML1DDist_ ## SIG \
  ( ElDistMatrix_ ## SIG H, ElInt nx, CREFLECT(T) omega, \
    ElInt numPmlPoints, Base<T> sigma, Base<T> pmlExp ) \
  { EL_TRY( \
      HelmholtzPML( \
        *CReflect(H), nx, CReflect(omega), \
        numPmlPoints, sigma, pmlExp ) ) } \
  ElError ElHelmholtzPML1DSparse_ ## SIG \
  ( ElSparseMatrix_ ## SIG H, ElInt nx, CREFLECT(T) omega, \
    ElInt numPmlPoints, Base<T> sigma, Base<T> pmlExp ) \
  { EL_TRY( \
      HelmholtzPML( \
        *CReflect(H), nx, CReflect(omega), \
        numPmlPoints, sigma, pmlExp ) ) } \
  ElError ElHelmholtzPML1DDistSparse_ ## SIG \
  ( ElDistSparseMatrix_ ## SIG H, ElInt nx, CREFLECT(T) omega, \
    ElInt numPmlPoints, Base<T> sigma, Base<T> pmlExp ) \
  { EL_TRY( \
      HelmholtzPML( \
        *CReflect(H), nx, CReflect(omega), \
        numPmlPoints, sigma, pmlExp ) ) } \
  ElError ElHelmholtzPML2D_ ## SIG \
  ( ElMatrix_ ## SIG H, ElInt nx, ElInt ny, CREFLECT(T) omega, \
    ElInt numPmlPoints, Base<T> sigma, Base<T> pmlExp ) \
  { EL_TRY( \
      HelmholtzPML( \
        *CReflect(H), nx, ny, CReflect(omega), \
        numPmlPoints, sigma, pmlExp ) ) } \
  ElError ElHelmholtzPML2DDist_ ## SIG \
  ( ElDistMatrix_ ## SIG H, ElInt nx, ElInt ny, CREFLECT(T) omega, \
    ElInt numPmlPoints, Base<T> sigma, Base<T> pmlExp ) \
  { EL_TRY( \
      HelmholtzPML( \
        *CReflect(H), nx, ny, CReflect(omega), \
        numPmlPoints, sigma, pmlExp ) ) } \
  ElError ElHelmholtzPML2DSparse_ ## SIG \
  ( ElSparseMatrix_ ## SIG H, ElInt nx, ElInt ny, CREFLECT(T) omega, \
    ElInt numPmlPoints, Base<T> sigma, Base<T> pmlExp ) \
  { EL_TRY( \
      HelmholtzPML( \
        *CReflect(H), nx, ny, CReflect(omega), \
        numPmlPoints, sigma, pmlExp ) ) } \
  ElError ElHelmholtzPML2DDistSparse_ ## SIG \
  ( ElDistSparseMatrix_ ## SIG H, ElInt nx, ElInt ny, CREFLECT(T) omega, \
    ElInt numPmlPoints, Base<T> sigma, Base<T> pmlExp ) \
  { EL_TRY( \
      HelmholtzPML( \
        *CReflect(H), nx, ny, CReflect(omega), \
        numPmlPoints, sigma, pmlExp ) ) } \
  ElError ElHelmholtzPML3D_ ## SIG \
  ( ElMatrix_ ## SIG H, ElInt nx, ElInt ny, ElInt nz, CREFLECT(T) omega, \
    ElInt numPmlPoints, Base<T> sigma, Base<T> pmlExp ) \
  { EL_TRY( \
      HelmholtzPML( \
        *CReflect(H), nx, ny, nz, CReflect(omega), \
        numPmlPoints, sigma, pmlExp ) ) } \
  ElError ElHelmholtzPML3DDist_ ## SIG \
  ( ElDistMatrix_ ## SIG H, ElInt nx, ElInt ny, ElInt nz, CREFLECT(T) omega, \
    ElInt numPmlPoints, Base<T> sigma, Base<T> pmlExp ) \
  { EL_TRY( \
      HelmholtzPML( \
        *CReflect(H), nx, ny, nz, CReflect(omega), \
        numPmlPoints, sigma, pmlExp ) ) } \
  ElError ElHelmholtzPML3DSparse_ ## SIG \
  ( ElSparseMatrix_ ## SIG H, ElInt nx, ElInt ny, ElInt nz, CREFLECT(T) omega, \
    ElInt numPmlPoints, Base<T> sigma, Base<T> pmlExp ) \
  { EL_TRY( \
      HelmholtzPML( \
        *CReflect(H), nx, ny, nz, CReflect(omega), \
        numPmlPoints, sigma, pmlExp ) ) } \
  ElError ElHelmholtzPML3DDistSparse_ ## SIG \
  ( ElDistSparseMatrix_ ## SIG H, ElInt nx, ElInt ny, ElInt nz, \
    CREFLECT(T) omega, ElInt numPmlPoints, Base<T> sigma, Base<T> pmlExp ) \
  { EL_TRY( \
      HelmholtzPML( \
        *CReflect(H), nx, ny, nz, CReflect(omega), \
        numPmlPoints, sigma, pmlExp ) ) } \
  /* NormalUniformSpectrum */ \
  ElError ElNormalUniformSpectrum_ ## SIG \
  ( ElMatrix_ ## SIG A, ElInt n, CREFLECT(T) center, Base<T> radius ) \
  { EL_TRY( \
      NormalUniformSpectrum( \
        *CReflect(A), n, CReflect(center), radius ) ) } \
  ElError ElNormalUniformSpectrumDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElInt n, CREFLECT(T) center, Base<T> radius ) \
  { EL_TRY( \
      NormalUniformSpectrum( \
        *CReflect(A), n, CReflect(center), radius ) ) } \
  /* Trefethen-Embree */ \
  ElError ElTrefethenEmbree_ ## SIG ( ElMatrix_ ## SIG A, ElInt n ) \
  { EL_TRY( TrefethenEmbree( *CReflect(A), n ) ) } \
  ElError ElTrefethenEmbreeDist_ ## SIG ( ElDistMatrix_ ## SIG A, ElInt n ) \
  { EL_TRY( TrefethenEmbree( *CReflect(A), n ) ) } \
  /* Uniform Helmholtz Green's */ \
  ElError ElUniformHelmholtzGreens_ ## SIG \
  ( ElMatrix_ ## SIG A, ElInt n, Base<T> lambda ) \
  { EL_TRY( UniformHelmholtzGreens( *CReflect(A), n, lambda ) ) } \
  ElError ElUniformHelmholtzGreensDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElInt n, Base<T> lambda ) \
  { EL_TRY( UniformHelmholtzGreens( *CReflect(A), n, lambda ) ) } \
  /* Whale */ \
  ElError ElWhale_ ## SIG ( ElMatrix_ ## SIG A, ElInt n ) \
  { EL_TRY( Whale( *CReflect(A), n ) ) } \
  ElError ElWhaleDist_ ## SIG ( ElDistMatrix_ ## SIG A, ElInt n ) \
  { EL_TRY( Whale( *CReflect(A), n ) ) }

#include <El/macros/CInstantiate.h>

} // extern "C"
