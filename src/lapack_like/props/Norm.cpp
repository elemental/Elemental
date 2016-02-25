/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename F>
Base<F> Norm( const Matrix<F>& A, NormType type )
{
    DEBUG_ONLY(CSE cse("Norm"))
    Base<F> norm = 0;
    switch( type )
    {
    // The following norms are rather cheap to compute
    case ENTRYWISE_ONE_NORM:
        norm = EntrywiseNorm( A, Base<F>(1) );
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
    DEBUG_ONLY(CSE cse("SymmetricNorm"))
    Base<F> norm = 0;
    switch( type )
    {
    // The following norms are rather cheap to compute
    case FROBENIUS_NORM:
        norm = SymmetricFrobeniusNorm( uplo, A );
        break;
    case ENTRYWISE_ONE_NORM:
        norm = SymmetricEntrywiseNorm( uplo, A, Base<F>(1) );
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
    DEBUG_ONLY(CSE cse("HermitianNorm"))
    Base<F> norm = 0;
    switch( type )
    {
    // The following norms are rather cheap to compute
    case FROBENIUS_NORM:
        norm = HermitianFrobeniusNorm( uplo, A );
        break;
    case ENTRYWISE_ONE_NORM:
        norm = HermitianEntrywiseNorm( uplo, A, Base<F>(1) );
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

template<typename F> 
Base<F> Norm( const ElementalMatrix<F>& A, NormType type )
{
    DEBUG_ONLY(CSE cse("Norm"))
    Base<F> norm = 0;
    switch( type )
    {
    // The following norms are rather cheap to compute
    case FROBENIUS_NORM: 
        norm = FrobeniusNorm( A );
        break;
    case ENTRYWISE_ONE_NORM:
        norm = EntrywiseNorm( A, Base<F>(1) );
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

template<typename F>
Base<F> SymmetricNorm
( UpperOrLower uplo, const ElementalMatrix<F>& A, NormType type )
{
    DEBUG_ONLY(CSE cse("SymmetricNorm"))
    Base<F> norm = 0;
    switch( type )
    {
    // The following norms are rather cheap to compute
    case FROBENIUS_NORM:
        norm = SymmetricFrobeniusNorm( uplo, A );
        break;
    case ENTRYWISE_ONE_NORM:
        norm = SymmetricEntrywiseNorm( uplo, A, Base<F>(1) );
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
Base<F> HermitianNorm
( UpperOrLower uplo, const ElementalMatrix<F>& A, NormType type )
{
    DEBUG_ONLY(CSE cse("HermitianNorm"))
    Base<F> norm = 0;
    switch( type )
    {
    // The following norms are rather cheap to compute
    case FROBENIUS_NORM:
        norm = HermitianFrobeniusNorm( uplo, A );
        break;
    case ENTRYWISE_ONE_NORM:
        norm = HermitianEntrywiseNorm( uplo, A, Base<F>(1) );
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

#define PROTO(F) \
  template Base<F> Norm( const Matrix<F>& A, NormType type ); \
  template Base<F> Norm( const ElementalMatrix<F>& A, NormType type ); \
  template Base<F> HermitianNorm \
  ( UpperOrLower uplo, const Matrix<F>& A, NormType type ); \
  template Base<F> HermitianNorm \
  ( UpperOrLower uplo, const ElementalMatrix<F>& A, NormType type ); \
  template Base<F> SymmetricNorm \
  ( UpperOrLower uplo, const Matrix<F>& A, NormType type ); \
  template Base<F> SymmetricNorm \
  ( UpperOrLower uplo, const ElementalMatrix<F>& A, NormType type );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
