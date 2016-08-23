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
template<typename F>
void HermitianSVD
( UpperOrLower uplo,
  Matrix<F>& A,
  Matrix<F>& U,
  Matrix<Base<F>>& s,
  Matrix<F>& V,
  bool overwrite )
{
    DEBUG_CSE
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
    Matrix<Base<F>> sSgn(n,1);
    for( Int j=0; j<n; ++j )
    {
        const Base<F> sigma = s(j);
        if( sigma >= 0 )
        {
            sSgn(j) = Base<F>(1);
        }
        else
        {
            sSgn(j) = Base<F>(-1);
            s(j) = -sigma;
        }
    }

    auto pairs = TaggedSort( s, DESCENDING );
    Matrix<F> VPerm(n,n);
    Matrix<Base<F>> sSgnPerm(n,1);
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

template<typename F>
void HermitianSVD
( UpperOrLower uplo,
  const Matrix<F>& A,
        Matrix<F>& U,
        Matrix<Base<F>>& s,
        Matrix<F>& V )
{
    DEBUG_CSE
    auto ACopy( A );
    HermitianSVD( uplo, ACopy, U, s, V, true );
}

// TODO: Add support for HermitianEigSubset and HermitianEigCtrl
template<typename F>
void HermitianSVD
( UpperOrLower uplo,
  ElementalMatrix<F>& A,
  ElementalMatrix<F>& U,
  ElementalMatrix<Base<F>>& s, 
  ElementalMatrix<F>& V,
  bool overwrite )
{
    DEBUG_CSE
    const Int n = A.Height();
    if( !overwrite )
    {
        DistMatrix<F> ACopy( A );
        HermitianSVD( uplo, ACopy, U, s, V, true );
        return;
    }

    // Grab an eigenvalue decomposition of A
    HermitianEig( uplo, A, s, V );

    typedef Base<F> Real;
    DistMatrix<Real,STAR,STAR> sSgn( s );
    auto sgnLambda = []( Real sigma ) { return Sgn(sigma,false); };
    EntrywiseMap( sSgn, function<Real(Real)>(sgnLambda) );

    // Set the singular values to the absolute value of the eigenvalues
    auto absLambda = []( Real sigma ) { return Abs(sigma); };
    EntrywiseMap( s, function<Real(Real)>(absLambda) );

    auto pairs = TaggedSort( s, DESCENDING );
    DistMatrix<F,VC,STAR> V_VC_STAR( V );
    DistMatrix<F,VC,STAR> VPerm_VC_STAR(A.Grid());
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
    V = VPerm_VC_STAR;
    U = V;
    DiagonalScale( RIGHT, NORMAL, sSgnPerm, U );
}

template<typename F>
void HermitianSVD
( UpperOrLower uplo,
  const ElementalMatrix<F>& A,
        ElementalMatrix<F>& U,
        ElementalMatrix<Base<F>>& s, 
        ElementalMatrix<F>& V )
{
    DEBUG_CSE
    DistMatrix<F> ACopy( A );
    HermitianSVD( uplo, ACopy, U, s, V, true );
}

// Return the singular values
// ==========================
// NOTE: A is ovewritten with its packed reduction to tridiagonal form

template<typename F>
void HermitianSVD
( UpperOrLower uplo,
  Matrix<F>& A,
  Matrix<Base<F>>& s,
  bool overwrite )
{
    DEBUG_CSE
    if( !overwrite )
    {
        Matrix<F> ACopy( A );
        HermitianSVD( uplo, ACopy, s, true );
        return;
    }

    typedef Base<F> Real;

    // Grab the eigenvalues of A
    HermitianEig( uplo, A, s );

    // Set the singular values to the absolute value of the eigenvalues
    auto absLambda = []( Real sigma ) { return Abs(sigma); };
    EntrywiseMap( s, function<Real(Real)>(absLambda) );

    Sort( s, DESCENDING );
}

template<typename F>
void HermitianSVD
( UpperOrLower uplo, const Matrix<F>& A, Matrix<Base<F>>& s )
{
    DEBUG_CSE
    Matrix<F> ACopy( A );
    HermitianSVD( uplo, ACopy, s, true );
}

template<typename F>
void HermitianSVD
( UpperOrLower uplo,
  ElementalMatrix<F>& A,
  ElementalMatrix<Base<F>>& s,
  bool overwrite )
{
    DEBUG_CSE
    if( !overwrite )
    {
        DistMatrix<F> ACopy( A );
        HermitianSVD( uplo, ACopy, s, true );
        return;
    }

    // Grab the eigenvalues of A
    HermitianEig( uplo, A, s );

    // Set the singular values to the absolute value of the eigenvalues
    typedef Base<F> Real;
    auto absLambda = []( Real sigma ) { return Abs(sigma); };
    EntrywiseMap( s, function<Real(Real)>(absLambda) );

    Sort( s, DESCENDING );
}

template<typename F>
void HermitianSVD
( UpperOrLower uplo,
  const ElementalMatrix<F>& A,
        ElementalMatrix<Base<F>>& s )
{
    DEBUG_CSE
    DistMatrix<F> ACopy( A );
    HermitianSVD( uplo, ACopy, s, true );
}

#define PROTO(F) \
  template void HermitianSVD \
  ( UpperOrLower uplo, \
    Matrix<F>& A, \
    Matrix<Base<F>>& s, \
    bool overwrite ); \
  template void HermitianSVD \
  ( UpperOrLower uplo, \
    const Matrix<F>& A, \
          Matrix<Base<F>>& s ); \
  template void HermitianSVD \
  ( UpperOrLower uplo, \
    ElementalMatrix<F>& A, \
    ElementalMatrix<Base<F>>& s, \
    bool overwrite ); \
  template void HermitianSVD \
  ( UpperOrLower uplo, \
    const ElementalMatrix<F>& A, \
          ElementalMatrix<Base<F>>& s ); \
  template void HermitianSVD \
  ( UpperOrLower uplo, \
    Matrix<F>& A, \
    Matrix<F>& U, \
    Matrix<Base<F>>& s, \
    Matrix<F>& V, \
    bool overwrite ); \
  template void HermitianSVD \
  ( UpperOrLower uplo, \
    const Matrix<F>& A, \
          Matrix<F>& U, \
          Matrix<Base<F>>& s, \
          Matrix<F>& V ); \
  template void HermitianSVD \
  ( UpperOrLower uplo, \
    ElementalMatrix<F>& A, \
    ElementalMatrix<F>& U, \
    ElementalMatrix<Base<F>>& s, \
    ElementalMatrix<F>& V, \
    bool overwrite ); \
  template void HermitianSVD \
  ( UpperOrLower uplo, \
    const ElementalMatrix<F>& A, \
          ElementalMatrix<F>& U, \
          ElementalMatrix<Base<F>>& s, \
          ElementalMatrix<F>& V );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
