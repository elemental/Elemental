/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {

// Grab the full SVD of the general matrix A, A = U diag(s) V^H
// ============================================================

// NOTE: A is ovewritten with its packed reduction to tridiagonal form
template<typename Field>
void HermitianSVD
( UpperOrLower uplo,
  Matrix<Field>& A,
  Matrix<Field>& U,
  Matrix<Base<Field>>& s,
  Matrix<Field>& V,
  bool overwrite )
{
    EL_DEBUG_CSE
    if( !overwrite )
    {
        auto ACopy( A );
        HermitianSVD( uplo, ACopy, U, s, V, true );
        return;
    }

    // Grab an eigenvalue decomposition of A
    HermitianEig( uplo, A, s, V );

    const Int n = A.Height();
    U.Resize( n, n );
    Matrix<Base<Field>> sSgn(n,1);
    for( Int j=0; j<n; ++j )
    {
        const Base<Field> sigma = s(j);
        if( sigma >= 0 )
        {
            sSgn(j) = Base<Field>(1);
        }
        else
        {
            sSgn(j) = Base<Field>(-1);
            s(j) = -sigma;
        }
    }

    auto pairs = TaggedSort( s, DESCENDING );
    Matrix<Field> VPerm(n,n);
    Matrix<Base<Field>> sSgnPerm(n,1);
    for( Int j=0; j<n; ++j )
    {
        MemCopy( VPerm.Buffer(0,j), V.Buffer(0,pairs[j].index), n );
        s(j) = pairs[j].value;
        sSgnPerm(j) = sSgn(pairs[j].index);
    }
    V = VPerm;
    U = V;
    DiagonalScale( RIGHT, NORMAL, sSgnPerm, U );
}

template<typename Field>
void HermitianSVD
( UpperOrLower uplo,
  const Matrix<Field>& A,
        Matrix<Field>& U,
        Matrix<Base<Field>>& s,
        Matrix<Field>& V )
{
    EL_DEBUG_CSE
    auto ACopy( A );
    HermitianSVD( uplo, ACopy, U, s, V, true );
}

// TODO(poulson): Add support for HermitianEigSubset and HermitianEigCtrl
template<typename Field>
void HermitianSVD
( UpperOrLower uplo,
  AbstractDistMatrix<Field>& A,
  AbstractDistMatrix<Field>& U,
  AbstractDistMatrix<Base<Field>>& s,
  AbstractDistMatrix<Field>& V,
  bool overwrite )
{
    EL_DEBUG_CSE
    const Int n = A.Height();
    if( !overwrite )
    {
        DistMatrix<Field> ACopy( A );
        HermitianSVD( uplo, ACopy, U, s, V, true );
        return;
    }

    // Grab an eigenvalue decomposition of A
    HermitianEig( uplo, A, s, V );

    typedef Base<Field> Real;
    DistMatrix<Real,STAR,STAR> sSgn( s );
    auto sgnLambda = []( const Real& sigma ) { return Sgn(sigma,false); };
    EntrywiseMap( sSgn, MakeFunction(sgnLambda) );

    // Set the singular values to the absolute value of the eigenvalues
    auto absLambda = []( const Real& sigma ) { return Abs(sigma); };
    EntrywiseMap( s, MakeFunction(absLambda) );

    auto pairs = TaggedSort( s, DESCENDING );
    DistMatrix<Field,VC,STAR> V_VC_STAR( V );
    DistMatrix<Field,VC,STAR> VPerm_VC_STAR(A.Grid());
    DistMatrix<Real,STAR,STAR> sSgnPerm(n,1,A.Grid());
    VPerm_VC_STAR.AlignWith( V_VC_STAR );
    VPerm_VC_STAR.Resize( n, n );
    const Int nLocal = V_VC_STAR.LocalHeight();
    for( Int j=0; j<n; ++j )
    {
        MemCopy
        ( VPerm_VC_STAR.Buffer(0,j),
          V_VC_STAR.Buffer(0,pairs[j].index), nLocal );
        s.Set( j, 0, pairs[j].value );
        sSgnPerm.Set( j, 0, sSgn.Get(pairs[j].index,0) );
    }
    Copy( VPerm_VC_STAR, V );
    Copy( V, U );
    DiagonalScale( RIGHT, NORMAL, sSgnPerm, U );
}

template<typename Field>
void HermitianSVD
( UpperOrLower uplo,
  const AbstractDistMatrix<Field>& A,
        AbstractDistMatrix<Field>& U,
        AbstractDistMatrix<Base<Field>>& s,
        AbstractDistMatrix<Field>& V )
{
    EL_DEBUG_CSE
    DistMatrix<Field> ACopy( A );
    HermitianSVD( uplo, ACopy, U, s, V, true );
}

// Return the singular values
// ==========================
// NOTE: A is ovewritten with its packed reduction to tridiagonal form

template<typename Field>
void HermitianSVD
( UpperOrLower uplo,
  Matrix<Field>& A,
  Matrix<Base<Field>>& s,
  bool overwrite )
{
    EL_DEBUG_CSE
    if( !overwrite )
    {
        Matrix<Field> ACopy( A );
        HermitianSVD( uplo, ACopy, s, true );
        return;
    }

    typedef Base<Field> Real;

    // Grab the eigenvalues of A
    HermitianEig( uplo, A, s );

    // Set the singular values to the absolute value of the eigenvalues
    auto absLambda = []( const Real& sigma ) { return Abs(sigma); };
    EntrywiseMap( s, MakeFunction(absLambda) );

    Sort( s, DESCENDING );
}

template<typename Field>
void HermitianSVD
( UpperOrLower uplo, const Matrix<Field>& A, Matrix<Base<Field>>& s )
{
    EL_DEBUG_CSE
    Matrix<Field> ACopy( A );
    HermitianSVD( uplo, ACopy, s, true );
}

template<typename Field>
void HermitianSVD
( UpperOrLower uplo,
  AbstractDistMatrix<Field>& A,
  AbstractDistMatrix<Base<Field>>& s,
  bool overwrite )
{
    EL_DEBUG_CSE
    if( !overwrite )
    {
        DistMatrix<Field> ACopy( A );
        HermitianSVD( uplo, ACopy, s, true );
        return;
    }

    // Grab the eigenvalues of A
    HermitianEig( uplo, A, s );

    // Set the singular values to the absolute value of the eigenvalues
    typedef Base<Field> Real;
    auto absLambda = []( const Real& sigma ) { return Abs(sigma); };
    EntrywiseMap( s, MakeFunction(absLambda) );

    Sort( s, DESCENDING );
}

template<typename Field>
void HermitianSVD
( UpperOrLower uplo,
  const AbstractDistMatrix<Field>& A,
        AbstractDistMatrix<Base<Field>>& s )
{
    EL_DEBUG_CSE
    DistMatrix<Field> ACopy( A );
    HermitianSVD( uplo, ACopy, s, true );
}

#define PROTO(Field) \
  template void HermitianSVD \
  ( UpperOrLower uplo, \
    Matrix<Field>& A, \
    Matrix<Base<Field>>& s, \
    bool overwrite ); \
  template void HermitianSVD \
  ( UpperOrLower uplo, \
    const Matrix<Field>& A, \
          Matrix<Base<Field>>& s ); \
  template void HermitianSVD \
  ( UpperOrLower uplo, \
    AbstractDistMatrix<Field>& A, \
    AbstractDistMatrix<Base<Field>>& s, \
    bool overwrite ); \
  template void HermitianSVD \
  ( UpperOrLower uplo, \
    const AbstractDistMatrix<Field>& A, \
          AbstractDistMatrix<Base<Field>>& s ); \
  template void HermitianSVD \
  ( UpperOrLower uplo, \
    Matrix<Field>& A, \
    Matrix<Field>& U, \
    Matrix<Base<Field>>& s, \
    Matrix<Field>& V, \
    bool overwrite ); \
  template void HermitianSVD \
  ( UpperOrLower uplo, \
    const Matrix<Field>& A, \
          Matrix<Field>& U, \
          Matrix<Base<Field>>& s, \
          Matrix<Field>& V ); \
  template void HermitianSVD \
  ( UpperOrLower uplo, \
    AbstractDistMatrix<Field>& A, \
    AbstractDistMatrix<Field>& U, \
    AbstractDistMatrix<Base<Field>>& s, \
    AbstractDistMatrix<Field>& V, \
    bool overwrite ); \
  template void HermitianSVD \
  ( UpperOrLower uplo, \
    const AbstractDistMatrix<Field>& A, \
          AbstractDistMatrix<Field>& U, \
          AbstractDistMatrix<Base<Field>>& s, \
          AbstractDistMatrix<Field>& V );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
