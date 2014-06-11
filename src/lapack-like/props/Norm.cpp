/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

#include "./Norm/Entrywise.hpp"
#include "./Norm/EntrywiseOne.hpp"
#include "./Norm/Frobenius.hpp"
#include "./Norm/Infinity.hpp"
#include "./Norm/KyFan.hpp"
#include "./Norm/Max.hpp"
#include "./Norm/One.hpp"

#include "./Norm/Nuclear.hpp"
#include "./Norm/Schatten.hpp"
#include "./Norm/Two.hpp"

#include "./Norm/TwoEstimate.hpp"

#include "./Norm/Zero.hpp"

namespace El {

template<typename F>
Base<F> Norm( const Matrix<F>& A, NormType type )
{
    DEBUG_ONLY(CallStackEntry cse("Norm"))
    Base<F> norm = 0;
    switch( type )
    {
    // The following norms are rather cheap to compute
    case ENTRYWISE_ONE_NORM:
        norm = EntrywiseOneNorm( A );
        break;
    case FROBENIUS_NORM: 
        norm = FrobeniusNorm( A );
        break;
    case INFINITY_NORM:
        norm = InfinityNorm( A );
        break;
    case MAX_NORM:
        norm = MaxNorm( A );
        break;
    case ONE_NORM:
        norm = OneNorm( A );
        break;
    // The following two norms make use of an SVD
    case NUCLEAR_NORM:
        norm = NuclearNorm( A );
        break;
    case TWO_NORM:
        norm = TwoNorm( A );
        break;
    }
    return norm;
}

template<typename F>
Base<F> SymmetricNorm( UpperOrLower uplo, const Matrix<F>& A, NormType type )
{
    DEBUG_ONLY(CallStackEntry cse("SymmetricNorm"))
    Base<F> norm = 0;
    switch( type )
    {
    // The following norms are rather cheap to compute
    case FROBENIUS_NORM:
        norm = SymmetricFrobeniusNorm( uplo, A );
        break;
    case ENTRYWISE_ONE_NORM:
        norm = SymmetricEntrywiseOneNorm( uplo, A );
        break;
    case INFINITY_NORM:
        norm = SymmetricInfinityNorm( uplo, A );
        break;
    case MAX_NORM:
        norm = SymmetricMaxNorm( uplo, A );
        break;
    case ONE_NORM:
        norm = SymmetricOneNorm( uplo, A );
        break;
    // The following norms make use of an SVD
    case NUCLEAR_NORM:
        norm = SymmetricNuclearNorm( uplo, A );
        break;
    case TWO_NORM:
        norm = SymmetricTwoNorm( uplo, A );
        break;
    }
    return norm;
}

template<typename F>
Base<F> HermitianNorm( UpperOrLower uplo, const Matrix<F>& A, NormType type )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianNorm"))
    Base<F> norm = 0;
    switch( type )
    {
    // The following norms are rather cheap to compute
    case FROBENIUS_NORM:
        norm = HermitianFrobeniusNorm( uplo, A );
        break;
    case ENTRYWISE_ONE_NORM:
        norm = HermitianEntrywiseOneNorm( uplo, A );
        break;
    case INFINITY_NORM:
        norm = HermitianInfinityNorm( uplo, A );
        break;
    case MAX_NORM:
        norm = HermitianMaxNorm( uplo, A );
        break;
    case ONE_NORM:
        norm = HermitianOneNorm( uplo, A );
        break;
    // The following norms make use of an SVD
    case NUCLEAR_NORM:
        norm = HermitianNuclearNorm( uplo, A );
        break;
    case TWO_NORM:
        norm = HermitianTwoNorm( uplo, A );
        break;
    }
    return norm;
}

template<typename F,Dist U,Dist V> 
Base<F> Norm( const DistMatrix<F,U,V>& A, NormType type )
{
    DEBUG_ONLY(CallStackEntry cse("Norm"))
    Base<F> norm = 0;
    switch( type )
    {
    // The following norms are rather cheap to compute
    case FROBENIUS_NORM: 
        norm = FrobeniusNorm( A );
        break;
    case ENTRYWISE_ONE_NORM:
        norm = EntrywiseOneNorm( A );
        break;
    case INFINITY_NORM:
        norm = InfinityNorm( A );
        break;
    case MAX_NORM:
        norm = MaxNorm( A );
        break;
    case ONE_NORM:
        norm = OneNorm( A );
        break;
    // The following norms make use of an SVD
    case NUCLEAR_NORM:
        norm = NuclearNorm( A );
        break;
    case TWO_NORM:
        norm = TwoNorm( A );
        break;
    }
    return norm;
}

template<typename F,Dist U,Dist V>
Base<F> SymmetricNorm
( UpperOrLower uplo, const DistMatrix<F,U,V>& A, NormType type )
{
    DEBUG_ONLY(CallStackEntry cse("SymmetricNorm"))
    Base<F> norm = 0;
    switch( type )
    {
    // The following norms are rather cheap to compute
    case FROBENIUS_NORM:
        norm = SymmetricFrobeniusNorm( uplo, A );
        break;
    case ENTRYWISE_ONE_NORM:
        norm = SymmetricEntrywiseOneNorm( uplo, A );
        break;
    case INFINITY_NORM:
        norm = SymmetricInfinityNorm( uplo, A );
        break;
    case MAX_NORM:
        norm = SymmetricMaxNorm( uplo, A );
        break;
    case ONE_NORM:
        norm = SymmetricOneNorm( uplo, A );
        break;
    // The following norms make use of an SVD
    case NUCLEAR_NORM:
        norm = SymmetricNuclearNorm( uplo, A );
        break;
    case TWO_NORM:
        norm = SymmetricTwoNorm( uplo, A );
        break;
    }
    return norm;
}

template<typename F,Dist U,Dist V>
Base<F> HermitianNorm
( UpperOrLower uplo, const DistMatrix<F,U,V>& A, NormType type )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianNorm"))
    Base<F> norm = 0;
    switch( type )
    {
    // The following norms are rather cheap to compute
    case FROBENIUS_NORM:
        norm = HermitianFrobeniusNorm( uplo, A );
        break;
    case ENTRYWISE_ONE_NORM:
        norm = HermitianEntrywiseOneNorm( uplo, A );
        break;
    case INFINITY_NORM:
        norm = HermitianInfinityNorm( uplo, A );
        break;
    case MAX_NORM:
        norm = HermitianMaxNorm( uplo, A );
        break;
    case ONE_NORM:
        norm = HermitianOneNorm( uplo, A );
        break;
    // The following norms make use of an SVD
    case NUCLEAR_NORM:
        norm = HermitianNuclearNorm( uplo, A );
        break;
    case TWO_NORM:
        norm = HermitianTwoNorm( uplo, A );
        break;
    }
    return norm;
}

#define PROTO_DIST(F,U,V) \
  template Base<F> Norm( const DistMatrix<F,U,V>& A, NormType type ); \
  template Base<F> HermitianNorm \
  ( UpperOrLower uplo, const DistMatrix<F,U,V>& A, NormType type ); \
  template Base<F> SymmetricNorm \
  ( UpperOrLower uplo, const DistMatrix<F,U,V>& A, NormType type ); \
  template Base<F> KyFanNorm( const DistMatrix<F,U,V>& A, Int k ); \
  template Base<F> HermitianKyFanNorm \
  ( UpperOrLower uplo, const DistMatrix<F,U,V>& A, Int k ); \
  template Base<F> SymmetricKyFanNorm \
  ( UpperOrLower uplo, const DistMatrix<F,U,V>& A, Int k ); \
  template Base<F> NuclearNorm( const DistMatrix<F,U,V>& A ); \
  template Base<F> HermitianNuclearNorm \
  ( UpperOrLower uplo, const DistMatrix<F,U,V>& A ); \
  template Base<F> SymmetricNuclearNorm \
  ( UpperOrLower uplo, const DistMatrix<F,U,V>& A ); \
  template Base<F> SchattenNorm( const DistMatrix<F,U,V>& A, Base<F> p ); \
  template Base<F> HermitianSchattenNorm \
  ( UpperOrLower uplo, const DistMatrix<F,U,V>& A, Base<F> p ); \
  template Base<F> SymmetricSchattenNorm \
  ( UpperOrLower uplo, const DistMatrix<F,U,V>& A, Base<F> p ); \
  template Base<F> TwoNorm( const DistMatrix<F,U,V>& A ); \
  template Base<F> HermitianTwoNorm \
  ( UpperOrLower uplo, const DistMatrix<F,U,V>& A ); \
  template Base<F> SymmetricTwoNorm \
  ( UpperOrLower uplo, const DistMatrix<F,U,V>& A );

#define PROTO(F) \
  template Base<F> Norm( const Matrix<F>& A, NormType type ); \
  template Base<F> HermitianNorm \
  ( UpperOrLower uplo, const Matrix<F>& A, NormType type ); \
  template Base<F> SymmetricNorm \
  ( UpperOrLower uplo, const Matrix<F>& A, NormType type ); \
  template Base<F> EntrywiseNorm \
  ( const Matrix<F>& A, Base<F> p ); \
  template Base<F> EntrywiseNorm \
  ( const AbstractDistMatrix<F>& A, Base<F> p ); \
  template Base<F> HermitianEntrywiseNorm \
  ( UpperOrLower uplo, const Matrix<F>& A, Base<F> p ); \
  template Base<F> HermitianEntrywiseNorm \
  ( UpperOrLower uplo, const AbstractDistMatrix<F>& A, Base<F> p ); \
  template Base<F> SymmetricEntrywiseNorm \
  ( UpperOrLower uplo, const Matrix<F>& A, Base<F> p ); \
  template Base<F> SymmetricEntrywiseNorm \
  ( UpperOrLower uplo, const AbstractDistMatrix<F>& A, Base<F> p ); \
  template Base<F> EntrywiseOneNorm( const Matrix<F>& A ); \
  template Base<F> EntrywiseOneNorm ( const AbstractDistMatrix<F>& A ); \
  template Base<F> HermitianEntrywiseOneNorm \
  ( UpperOrLower uplo, const Matrix<F>& A ); \
  template Base<F> HermitianEntrywiseOneNorm \
  ( UpperOrLower uplo, const AbstractDistMatrix<F>& A ); \
  template Base<F> SymmetricEntrywiseOneNorm \
  ( UpperOrLower uplo, const Matrix<F>& A ); \
  template Base<F> SymmetricEntrywiseOneNorm \
  ( UpperOrLower uplo, const AbstractDistMatrix<F>& A ); \
  template Base<F> FrobeniusNorm( const Matrix<F>& A ); \
  template Base<F> FrobeniusNorm ( const AbstractDistMatrix<F>& A ); \
  template Base<F> HermitianFrobeniusNorm \
  ( UpperOrLower uplo, const Matrix<F>& A ); \
  template Base<F> HermitianFrobeniusNorm \
  ( UpperOrLower uplo, const AbstractDistMatrix<F>& A ); \
  template Base<F> SymmetricFrobeniusNorm \
  ( UpperOrLower uplo, const Matrix<F>& A ); \
  template Base<F> SymmetricFrobeniusNorm \
  ( UpperOrLower uplo, const AbstractDistMatrix<F>& A ); \
  template Base<F> InfinityNorm( const Matrix<F>& A ); \
  template Base<F> InfinityNorm ( const AbstractDistMatrix<F>& A ); \
  template Base<F> HermitianInfinityNorm \
  ( UpperOrLower uplo, const Matrix<F>& A ); \
  template Base<F> HermitianInfinityNorm \
  ( UpperOrLower uplo, const AbstractDistMatrix<F>& A ); \
  template Base<F> SymmetricInfinityNorm \
  ( UpperOrLower uplo, const Matrix<F>& A ); \
  template Base<F> SymmetricInfinityNorm \
  ( UpperOrLower uplo, const AbstractDistMatrix<F>& A ); \
  template Base<F> KyFanNorm( const Matrix<F>& A, Int k ); \
  template Base<F> HermitianKyFanNorm \
  ( UpperOrLower uplo, const Matrix<F>& A, Int k ); \
  template Base<F> SymmetricKyFanNorm \
  ( UpperOrLower uplo, const Matrix<F>& A, Int k ); \
  template Base<F> MaxNorm( const Matrix<F>& A ); \
  template Base<F> MaxNorm ( const AbstractDistMatrix<F>& A ); \
  template Base<F> HermitianMaxNorm \
  ( UpperOrLower uplo, const Matrix<F>& A ); \
  template Base<F> HermitianMaxNorm \
  ( UpperOrLower uplo, const AbstractDistMatrix<F>& A ); \
  template Base<F> SymmetricMaxNorm \
  ( UpperOrLower uplo, const Matrix<F>& A ); \
  template Base<F> SymmetricMaxNorm \
  ( UpperOrLower uplo, const AbstractDistMatrix<F>& A ); \
  template Base<F> NuclearNorm( const Matrix<F>& A ); \
  template Base<F> HermitianNuclearNorm \
  ( UpperOrLower uplo, const Matrix<F>& A ); \
  template Base<F> SymmetricNuclearNorm \
  ( UpperOrLower uplo, const Matrix<F>& A ); \
  template Base<F> OneNorm( const Matrix<F>& A ); \
  template Base<F> OneNorm ( const AbstractDistMatrix<F>& A ); \
  template Base<F> HermitianOneNorm \
  ( UpperOrLower uplo, const Matrix<F>& A ); \
  template Base<F> HermitianOneNorm \
  ( UpperOrLower uplo, const AbstractDistMatrix<F>& A ); \
  template Base<F> SymmetricOneNorm \
  ( UpperOrLower uplo, const Matrix<F>& A ); \
  template Base<F> SymmetricOneNorm \
  ( UpperOrLower uplo, const AbstractDistMatrix<F>& A ); \
  template Base<F> SchattenNorm( const Matrix<F>& A, Base<F> p ); \
  template Base<F> HermitianSchattenNorm \
  ( UpperOrLower uplo, const Matrix<F>& A, Base<F> p ); \
  template Base<F> SymmetricSchattenNorm \
  ( UpperOrLower uplo, const Matrix<F>& A, Base<F> p ); \
  template Base<F> TwoNormEstimate \
  ( const Matrix<F>& A, Base<F> tol, Int maxIts ); \
  template Base<F> TwoNormEstimate \
  ( const DistMatrix<F>& A, Base<F> tol, Int maxIts ); \
  template Base<F> HermitianTwoNormEstimate \
  ( UpperOrLower uplo, const Matrix<F>& A, Base<F> tol, Int maxIts ); \
  template Base<F> HermitianTwoNormEstimate \
  ( UpperOrLower uplo, const DistMatrix<F>& A, Base<F> tol, Int maxIts ); \
  template Base<F> SymmetricTwoNormEstimate \
  ( UpperOrLower uplo, const Matrix<F>& A, Base<F> tol, Int maxIts ); \
  template Base<F> SymmetricTwoNormEstimate \
  ( UpperOrLower uplo, const DistMatrix<F>& A, Base<F> tol, Int maxIts ); \
  template Int ZeroNorm( const Matrix<F>& A ); \
  template Int ZeroNorm( const AbstractDistMatrix<F>& A ); \
  PROTO_DIST(F,CIRC,CIRC) \
  PROTO_DIST(F,MC,  MR  ) \
  PROTO_DIST(F,MC,  STAR) \
  PROTO_DIST(F,MD,  STAR) \
  PROTO_DIST(F,MR,  MC  ) \
  PROTO_DIST(F,MR,  STAR) \
  PROTO_DIST(F,STAR,MC  ) \
  PROTO_DIST(F,STAR,MD  ) \
  PROTO_DIST(F,STAR,MR  ) \
  PROTO_DIST(F,STAR,STAR) \
  PROTO_DIST(F,STAR,VC  ) \
  PROTO_DIST(F,STAR,VR  ) \
  PROTO_DIST(F,VC,  STAR) \
  PROTO_DIST(F,VR,  STAR)

PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

template Int ZeroNorm( const Matrix<Int>& A );
template Int ZeroNorm( const AbstractDistMatrix<Int>& A );
template Int MaxNorm( const Matrix<Int>& A );
template Int MaxNorm ( const AbstractDistMatrix<Int>& A );
template Int HermitianMaxNorm
( UpperOrLower uplo, const Matrix<Int>& A );
template Int HermitianMaxNorm
( UpperOrLower uplo, const AbstractDistMatrix<Int>& A );
template Int SymmetricMaxNorm
( UpperOrLower uplo, const Matrix<Int>& A );
template Int SymmetricMaxNorm
( UpperOrLower uplo, const AbstractDistMatrix<Int>& A );

} // namespace El
