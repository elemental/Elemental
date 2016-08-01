/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SCHUR_QUASITRIANGEIG_HPP
#define EL_SCHUR_QUASITRIANGEIG_HPP

namespace El {
namespace schur {

template<typename F>
void QuasiTriangEig
( const Matrix<F>& dMain,
  const Matrix<F>& dSub,
  const Matrix<F>& dSup,
  Matrix<Complex<Base<F>>>& w )
{
    DEBUG_CSE
    const Int n = dMain.Height();
    Matrix<F> H11(2,2);
    w.Resize( n, 1 );

    Int j=0;
    while( j < n )
    {
        if( j == n-1 || dSub(j) == F(0) )
        {
            w(j) = dMain(j);
            ++j;
        }
        else
        {
            H11(0,0) = dMain(j);
            H11(1,0) = dSub(j);
            H11(0,1) = dSup(j);
            H11(1,1) = dMain(j+1);
            lapack::HessenbergEig( 2, H11.Buffer(), H11.LDim(), w.Buffer(j,0) );
            j += 2;
        }
    }
}

template<typename F>
void QuasiTriangEig( const Matrix<F>& U, Matrix<Complex<Base<F>>>& w )
{
    DEBUG_CSE
    auto dMain = GetDiagonal(U);
    auto dSub = GetDiagonal(U,-1);
    auto dSup = GetDiagonal(U,+1);
    QuasiTriangEig( dMain, dSub, dSup, w );
}

template<typename F>
Matrix<Complex<Base<F>>> QuasiTriangEig( const Matrix<F>& U )
{
    DEBUG_CSE
    Matrix<Complex<Base<F>>> w;
    QuasiTriangEig( U, w );
    return w;
}

template<typename F>
void QuasiTriangEig
( const ElementalMatrix<F>& UPre,
        ElementalMatrix<Complex<Base<F>>>& w )
{
    DEBUG_CSE

    DistMatrixReadProxy<F,F,MC,MR> UProx( UPre );
    auto& U = UProx.GetLocked();

    const Grid& g = U.Grid();
    DistMatrix<F,STAR,STAR> dMain(g), dSub(g), dSup(g);
    DistMatrix<Complex<Base<F>>,STAR,STAR> w_STAR_STAR(g);
    dMain = GetDiagonal(U);
    dSub = GetDiagonal(U,-1);
    dSup = GetDiagonal(U,+1);
    w_STAR_STAR.Resize( U.Height(), 1 );
    QuasiTriangEig
    ( dMain.Matrix(), dSub.Matrix(), dSup.Matrix(), w_STAR_STAR.Matrix() );
    Copy( w_STAR_STAR, w );
}

template<typename F>
DistMatrix<Complex<Base<F>>,VR,STAR> 
QuasiTriangEig( const ElementalMatrix<F>& U )
{
    DEBUG_CSE
    DistMatrix<Complex<Base<F>>,VR,STAR> w(U.Grid());
    QuasiTriangEig( U, w );
    return w;
}

} // namespace schur
} // namespace El

#endif // ifndef EL_SCHUR_QUASITRIANGEIG_HPP
