/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

extern "C" {

void EL_BLAS(strmm)
( const char* side, const char* uplo, const char* trans, const char* diag,
  const BlasInt* m, const BlasInt* n,
  const float* alpha, const float* A, const BlasInt* ALDim,
                            float* B, const BlasInt* BLDim );
void EL_BLAS(dtrmm)
( const char* side, const char* uplo, const char* trans, const char* diag,
  const BlasInt* m, const BlasInt* n,
  const double* alpha, const double* A, const BlasInt* ALDim,
                             double* B, const BlasInt* BLDim );
void EL_BLAS(ctrmm)
( const char* side, const char* uplo, const char* trans, const char* diag,
  const BlasInt* m, const BlasInt* n,
  const scomplex* alpha,
  const scomplex* A, const BlasInt* ALDim,
        scomplex* B, const BlasInt* BLDim );
void EL_BLAS(ztrmm)
( const char* side, const char* uplo, const char* trans, const char* diag,
  const BlasInt* m, const BlasInt* n,
  const dcomplex* alpha,
  const dcomplex* A, const BlasInt* ALDim,
        dcomplex* B, const BlasInt* BLDim );

} // extern "C"

namespace El {
namespace blas {

template<typename T>
void Trmm
( char side, char uplo, char trans, char unit,
  BlasInt m, BlasInt n,
  const T& alpha,
  const T* A, BlasInt ALDim,
        T* B, BlasInt BLDim )
{
    const bool onLeft = ( std::toupper(side) == 'L' );
    const bool conjugate = ( std::toupper(trans) == 'C' );

    // Scale B
    for( BlasInt j=0; j<n; ++j )
        for( BlasInt i=0; i<m; ++i )
            B[i+j*BLDim] *= alpha;

    // TODO: Legitimate blocked implementations...it seems offensive to
    //       repeatedly stream all of the triangular matrix through memory
    //       for each row/column of B

    if( onLeft )
    {
        for( BlasInt j=0; j<n; ++j )
            Trmv( uplo, trans, unit, m, A, ALDim, &B[j*BLDim], 1 );
    }
    else
    {
        char newTrans = ( std::toupper(trans) == 'N' ? 'T' : 'N' );
        for( BlasInt i=0; i<m; ++i )
        {
            if( conjugate )
                for( Int j=0; j<n; ++j )
                    B[i+j*BLDim] = Conj(B[i+j*BLDim]);
            Trmv( uplo, newTrans, unit, n, A, ALDim, &B[i], BLDim );
            // TODO: Avoid temporaries here
            if( conjugate )
                for( Int j=0; j<n; ++j )
                    B[i+j*BLDim] = Conj(B[i+j*BLDim]);
        }
    }
}
template void Trmm
( char side, char uplo, char trans, char unit,
  BlasInt m, BlasInt n,
  const Int& alpha,
  const Int* A, BlasInt ALDim,
        Int* B, BlasInt BLDim );
#ifdef EL_HAVE_QD
template void Trmm
( char side, char uplo, char trans, char unit,
  BlasInt m, BlasInt n,
  const DoubleDouble& alpha,
  const DoubleDouble* A, BlasInt ALDim,
        DoubleDouble* B, BlasInt BLDim );
template void Trmm
( char side, char uplo, char trans, char unit,
  BlasInt m, BlasInt n,
  const QuadDouble& alpha,
  const QuadDouble* A, BlasInt ALDim,
        QuadDouble* B, BlasInt BLDim );
template void Trmm
( char side, char uplo, char trans, char unit,
  BlasInt m, BlasInt n,
  const Complex<DoubleDouble>& alpha,
  const Complex<DoubleDouble>* A, BlasInt ALDim,
        Complex<DoubleDouble>* B, BlasInt BLDim );
template void Trmm
( char side, char uplo, char trans, char unit,
  BlasInt m, BlasInt n,
  const Complex<QuadDouble>& alpha,
  const Complex<QuadDouble>* A, BlasInt ALDim,
        Complex<QuadDouble>* B, BlasInt BLDim );
#endif
#ifdef EL_HAVE_QUAD
template void Trmm
( char side, char uplo, char trans, char unit,
  BlasInt m, BlasInt n,
  const Quad& alpha,
  const Quad* A, BlasInt ALDim,
        Quad* B, BlasInt BLDim );
template void Trmm
( char side, char uplo, char trans, char unit,
  BlasInt m, BlasInt n,
  const Complex<Quad>& alpha,
  const Complex<Quad>* A, BlasInt ALDim,
        Complex<Quad>* B, BlasInt BLDim );
#endif
#ifdef EL_HAVE_MPC
template void Trmm
( char side, char uplo, char trans, char unit,
  BlasInt m, BlasInt n,
  const BigInt& alpha,
  const BigInt* A, BlasInt ALDim,
        BigInt* B, BlasInt BLDim );
template void Trmm
( char side, char uplo, char trans, char unit,
  BlasInt m, BlasInt n,
  const BigFloat& alpha,
  const BigFloat* A, BlasInt ALDim,
        BigFloat* B, BlasInt BLDim );
template void Trmm
( char side, char uplo, char trans, char unit,
  BlasInt m, BlasInt n,
  const Complex<BigFloat>& alpha,
  const Complex<BigFloat>* A, BlasInt ALDim,
        Complex<BigFloat>* B, BlasInt BLDim );
#endif

void Trmm
( char side, char uplo, char trans, char unit, BlasInt m, BlasInt n,
  const float& alpha,
  const float* A, BlasInt ALDim,
        float* B, BlasInt BLDim )
{
    const char fixedTrans = ( std::toupper(trans) == 'C' ? 'T' : trans );    
    EL_BLAS(strmm)
    ( &side, &uplo, &fixedTrans, &unit, &m, &n, &alpha, A, &ALDim, B, &BLDim );
}

void Trmm
( char side, char uplo, char trans, char unit, BlasInt m, BlasInt n,
  const double& alpha,
  const double* A, BlasInt ALDim,
        double* B, BlasInt BLDim )
{
    const char fixedTrans = ( std::toupper(trans) == 'C' ? 'T' : trans );    
    EL_BLAS(dtrmm)
    ( &side, &uplo, &fixedTrans, &unit, &m, &n, &alpha, A, &ALDim, B, &BLDim );
}

void Trmm
( char side, char uplo, char trans, char unit, BlasInt m, BlasInt n,
  const scomplex& alpha,
  const scomplex* A, BlasInt ALDim,
        scomplex* B, BlasInt BLDim )
{
    EL_BLAS(ctrmm)
    ( &side, &uplo, &trans, &unit, &m, &n, &alpha, A, &ALDim, B, &BLDim );
}

void Trmm
( char side, char uplo, char trans, char unit, BlasInt m, BlasInt n,
  const dcomplex& alpha,
  const dcomplex* A, BlasInt ALDim,
        dcomplex* B, BlasInt BLDim )
{
    EL_BLAS(ztrmm)
    ( &side, &uplo, &trans, &unit, &m, &n, &alpha, A, &ALDim, B, &BLDim );
}

} // namespace blas
} // namespace El
