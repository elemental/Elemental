/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

extern "C" {

void EL_BLAS(strsm)
( const char* side, const char* uplo, const char* transA, const char* diag,
  const BlasInt* m, const BlasInt* n,
  const float* alpha, const float* A, const BlasInt* ALDim,
                            float* B, const BlasInt* BLDim );
void EL_BLAS(dtrsm)
( const char* side, const char* uplo, const char* transA, const char* diag,
  const BlasInt* m, const BlasInt* n,
  const double* alpha, const double* A, const BlasInt* ALDim,
                             double* B, const BlasInt* BLDim );
void EL_BLAS(ctrsm)
( const char* side, const char* uplo, const char* transA, const char* diag,
  const BlasInt* m, const BlasInt* n,
  const scomplex* alpha,
  const scomplex* A, const BlasInt* ALDim,
        scomplex* B, const BlasInt* BLDim );
void EL_BLAS(ztrsm)
( const char* side, const char* uplo, const char* transA, const char* diag,
  const BlasInt* m, const BlasInt* n,
  const dcomplex* alpha,
  const dcomplex* A, const BlasInt* ALDim,
        dcomplex* B, const BlasInt* BLDim );

} // extern "C"

namespace El {
namespace blas {

template<typename F>
void Trsm
( char side, char uplo, char trans, char unit,
  BlasInt m, BlasInt n,
  const F& alpha,
  const F* A, BlasInt ALDim,
        F* B, BlasInt BLDim )
{
    // NOTE: Temporaries are avoided since constructing a BigInt/BigFloat
    //       involves a memory allocation
    const bool onLeft = ( std::toupper(side) == 'L' );
    const bool lower = ( std::toupper(uplo) == 'L' );
    const bool conjugate = ( std::toupper(trans) == 'C' );
    const bool unitDiag = ( std::toupper(unit) == 'U' );

    // Scale B
    for( BlasInt j=0; j<n; ++j )
        for( BlasInt i=0; i<m; ++i )
            B[i+j*BLDim] *= alpha;

    F alpha11, alpha11Conj;
    if( onLeft )
    {
        if( std::toupper(trans) == 'N' )
        {
            if( lower )
            {
                for( BlasInt k=0; k<m; ++k )     
                {
                    if( !unitDiag )
                    {
                        alpha11 = A[k+k*ALDim];
                        for( BlasInt j=0; j<n; ++j )      
                            B[k+j*BLDim] /= alpha11;
                    }
                    Geru
                    ( m-(k+1), n,
                      F(-1), &A[(k+1)+k*ALDim], 1,
                             &B[k],             BLDim,
                             &B[k+1],           BLDim );
                }
            }
            else
            {
                for( BlasInt k=m-1; k>=0; --k )     
                {
                    if( !unitDiag )
                    {
                        alpha11 = A[k+k*ALDim];
                        for( BlasInt j=0; j<n; ++j )      
                            B[k+j*BLDim] /= alpha11;
                    }
                    Geru
                    ( k, n,
                      F(-1), &A[k*ALDim], 1,
                             &B[k],       BLDim,
                              B,          BLDim );
                }
            }
        }
        else if( conjugate )
        {
            if( lower )
            {
               std::vector<F> aRow(m);
               for( BlasInt k=m-1; k>=0; --k )
                {
                    if( !unitDiag )
                    {
                        Conj( A[k+k*ALDim], alpha11Conj );
                        for( BlasInt j=0; j<n; ++j )
                            B[k+j*BLDim] /= alpha11Conj;
                    }
                    for( Int s=0; s<k; ++s )
                        Conj( A[k+s*ALDim], aRow[s] );
                    Geru
                    ( k, n,
                      F(-1), aRow.data(), 1,
                             &B[k],       BLDim,
                              B,          BLDim );
                }
            }
            else
            {
                for( BlasInt k=0; k<m; ++k )
                {
                    if( !unitDiag )
                    {
                        Conj( A[k+k*ALDim], alpha11Conj );
                        for( BlasInt j=0; j<n; ++j )
                            B[k+j*BLDim] /= alpha11Conj;
                    }
                    std::vector<F> aRow(m);
                    for( Int s=0; s<m-(k+1); ++s )
                        Conj( A[k+(k+1+s)*ALDim], aRow[s] );
                    Geru
                    ( m-(k+1), n,
                      F(-1), aRow.data(), 1,
                             &B[k],       BLDim,
                             &B[k+1],     BLDim );
                }
            }
        }
        else
        {
            if( lower )
            {
                for( BlasInt k=m-1; k>=0; --k )
                {
                    if( !unitDiag )
                    {
                        alpha11 = A[k+k*ALDim];
                        for( BlasInt j=0; j<n; ++j )
                            B[k+j*BLDim] /= alpha11;
                    }
                    Geru
                    ( k, n,
                      F(-1), &A[k], ALDim,
                             &B[k], BLDim,
                              B,    BLDim );
                }
            }
            else
            {
                for( BlasInt k=0; k<m; ++k )
                {
                    if( !unitDiag )
                    {
                        alpha11 = A[k+k*ALDim];
                        for( BlasInt j=0; j<n; ++j )
                            B[k+j*BLDim] /= alpha11;
                    }
                    Geru
                    ( m-(k+1), n,
                      F(-1), &A[k+(k+1)*ALDim], ALDim,
                             &B[k],             BLDim,
                             &B[k+1],           BLDim );
                }
            }
        }
    }
    else
    {
        if( std::toupper(trans) == 'N' )
        {
            if( lower )
            {
                for( BlasInt k=n-1; k>=0; --k )
                {
                    if( !unitDiag )
                    {
                        alpha11 = A[k+k*ALDim];
                        for( BlasInt i=0; i<m; ++i )
                            B[i+k*BLDim] /= alpha11;
                    }
                    blas::Geru
                    ( m, k,
                      F(-1), &B[k*BLDim], 1,
                             &A[k],       ALDim,
                              B,          BLDim );
                }
            }
            else
            {
                for( BlasInt k=0; k<n; ++k )
                {
                    if( !unitDiag )
                    {
                        alpha11 = A[k+k*ALDim];
                        for( BlasInt i=0; i<m; ++i )
                            B[i+k*BLDim] /= alpha11;
                    }
                    blas::Geru
                    ( m, n-(k+1),
                      F(-1), &B[k*BLDim],       1,
                             &A[k+(k+1)*ALDim], ALDim,
                             &B[(k+1)*BLDim],   BLDim );
                }
            }
        }
        else if( conjugate )
        {
            if( lower )
            {
                for( BlasInt k=0; k<n; ++k )
                {
                    if( !unitDiag )
                    {
                        Conj( A[k+k*ALDim], alpha11Conj );
                        for( BlasInt i=0; i<m; ++i )
                            B[i+k*BLDim] /= alpha11Conj;
                    }
                    blas::Ger
                    ( m, n-(k+1),
                      F(-1), &B[k*BLDim],       1,
                             &A[(k+1)+k*ALDim], 1,
                             &B[(k+1)*BLDim],   BLDim );
                }
            }
            else
            {
                for( BlasInt k=n-1; k>=0; --k )
                {
                    if( !unitDiag )
                    {
                        Conj( A[k+k*ALDim], alpha11Conj );
                        for( BlasInt i=0; i<m; ++i )
                            B[i+k*BLDim] /= alpha11Conj;
                    }
                    blas::Ger
                    ( m, k,
                      F(-1), &B[k*BLDim], 1,
                             &A[k*ALDim], 1,
                              B,          BLDim );
                }
            }
        }
        else
        {
            if( lower )
            {
                for( BlasInt k=0; k<n; ++k )
                {
                    if( !unitDiag )
                    {
                        alpha11 = A[k+k*ALDim];
                        for( BlasInt i=0; i<m; ++i )
                            B[i+k*BLDim] /= alpha11;
                    }
                    blas::Geru
                    ( m, n-(k+1),
                      F(-1), &B[k*BLDim],       1,
                             &A[(k+1)+k*ALDim], 1,
                             &B[(k+1)*BLDim],   BLDim );
                }
            }
            else
            {
                for( BlasInt k=n-1; k>=0; --k )
                {
                    if( !unitDiag )
                    {
                        alpha11 = A[k+k*ALDim];
                        for( BlasInt i=0; i<m; ++i )
                            B[i+k*BLDim] /= alpha11;
                    }
                    blas::Geru
                    ( m, k,
                      F(-1), &B[k*BLDim], 1,
                             &A[k*ALDim], 1,
                              B,          BLDim );
                }
            }
        }
    }
}
#ifdef EL_HAVE_QD
template void Trsm
( char side, char uplo, char trans, char unit,
  BlasInt m, BlasInt n,
  const DoubleDouble& alpha,
  const DoubleDouble* A, BlasInt ALDim,
        DoubleDouble* B, BlasInt BLDim );
template void Trsm
( char side, char uplo, char trans, char unit,
  BlasInt m, BlasInt n,
  const QuadDouble& alpha,
  const QuadDouble* A, BlasInt ALDim,
        QuadDouble* B, BlasInt BLDim );
template void Trsm
( char side, char uplo, char trans, char unit,
  BlasInt m, BlasInt n,
  const Complex<DoubleDouble>& alpha,
  const Complex<DoubleDouble>* A, BlasInt ALDim,
        Complex<DoubleDouble>* B, BlasInt BLDim );
template void Trsm
( char side, char uplo, char trans, char unit,
  BlasInt m, BlasInt n,
  const Complex<QuadDouble>& alpha,
  const Complex<QuadDouble>* A, BlasInt ALDim,
        Complex<QuadDouble>* B, BlasInt BLDim );
#endif
#ifdef EL_HAVE_QUAD
template void Trsm
( char side, char uplo, char trans, char unit,
  BlasInt m, BlasInt n,
  const Quad& alpha,
  const Quad* A, BlasInt ALDim,
        Quad* B, BlasInt BLDim );
template void Trsm
( char side, char uplo, char trans, char unit,
  BlasInt m, BlasInt n,
  const Complex<Quad>& alpha,
  const Complex<Quad>* A, BlasInt ALDim,
        Complex<Quad>* B, BlasInt BLDim );
#endif
#ifdef EL_HAVE_MPC
template void Trsm
( char side, char uplo, char trans, char unit,
  BlasInt m, BlasInt n,
  const BigFloat& alpha,
  const BigFloat* A, BlasInt ALDim,
        BigFloat* B, BlasInt BLDim );
template void Trsm
( char side, char uplo, char trans, char unit,
  BlasInt m, BlasInt n,
  const Complex<BigFloat>& alpha,
  const Complex<BigFloat>* A, BlasInt ALDim,
        Complex<BigFloat>* B, BlasInt BLDim );
#endif

void Trsm
( char side, char uplo, char trans, char unit,
  BlasInt m, BlasInt n,
  const float& alpha,
  const float* A, BlasInt ALDim,
        float* B, BlasInt BLDim )
{
    const char fixedTrans = ( std::toupper(trans) == 'C' ? 'T' : trans );
    EL_BLAS(strsm)
    ( &side, &uplo, &fixedTrans, &unit, &m, &n, &alpha, A, &ALDim, B, &BLDim );
} 

void Trsm
( char side, char uplo, char trans, char unit,
  BlasInt m, BlasInt n,
  const double& alpha,
  const double* A, BlasInt ALDim,
        double* B, BlasInt BLDim )
{
    const char fixedTrans = ( std::toupper(trans) == 'C' ? 'T' : trans );
    EL_BLAS(dtrsm)
    ( &side, &uplo, &fixedTrans, &unit, &m, &n, &alpha, A, &ALDim, B, &BLDim );
} 

void Trsm
( char side, char uplo, char trans, char unit,
  BlasInt m, BlasInt n,
  const scomplex& alpha,
  const scomplex* A, BlasInt ALDim,
        scomplex* B, BlasInt BLDim )
{
    EL_BLAS(ctrsm)
    ( &side, &uplo, &trans, &unit, &m, &n, &alpha, A, &ALDim, B, &BLDim );
} 

void Trsm
( char side, char uplo, char trans, char unit,
  BlasInt m, BlasInt n,
  const dcomplex& alpha,
  const dcomplex* A, BlasInt ALDim,
        dcomplex* B, BlasInt BLDim )
{
    EL_BLAS(ztrsm)
    ( &side, &uplo, &trans, &unit, &m, &n, &alpha, A, &ALDim, B, &BLDim );
} 

} // namespace blas
} // namespace El
