/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_SCHUR_QUASITRIANGEIG_HPP
#define EL_SCHUR_QUASITRIANGEIG_HPP

namespace El {
namespace schur {

template<typename F>
void QuasiTriangEig
( const Matrix<F>& dMain, const Matrix<F>& dSub, const Matrix<F>& dSup,
  Matrix<Complex<Base<F>>>& w )
{
    DEBUG_ONLY(CallStackEntry cse("schur::QuasiTriangEig"))
    const Int n = dMain.Height();
    Matrix<F> H11(2,2);
    w.Resize( n, 1 );

    Int j=0;
    while( j < n )
    {
        if( j == n-1 || dSub.Get(j,0) == F(0) )
        {
            w.Set( j, 0, dMain.Get(j,0) );
            ++j;
        }
        else
        {
            H11.Set( 0, 0, dMain.Get(j,0) );
            H11.Set( 1, 0, dSub.Get(j,0) );
            H11.Set( 0, 1, dSup.Get(j,0) );
            H11.Set( 1, 1, dMain.Get(j+1,0) );
            lapack::HessenbergEig( 2, H11.Buffer(), H11.LDim(), w.Buffer(j,0) );
            j += 2;
        }
    }
}

template<typename F>
void QuasiTriangEig( const Matrix<F>& U, Matrix<Complex<Base<F>>>& w )
{
    DEBUG_ONLY(CallStackEntry cse("schur::QuasiTriangEig"))
    auto dMain = U.GetDiagonal();
    auto dSub = U.GetDiagonal(-1);
    auto dSup = U.GetDiagonal(+1);
    QuasiTriangEig( dMain, dSub, dSup, w );
}

template<typename F>
Matrix<Complex<Base<F>>> QuasiTriangEig( const Matrix<F>& U )
{
    DEBUG_ONLY(CallStackEntry cse("schur::QuasiTriangEig"))
    Matrix<Complex<Base<F>>> w;
    QuasiTriangEig( U, w );
    return w;
}

template<typename F,Dist colDist,Dist rowDist>
void QuasiTriangEig
( const DistMatrix<F>& U, DistMatrix<Complex<Base<F>>,colDist,rowDist>& w )
{
    DEBUG_ONLY(CallStackEntry cse("schur::QuasiTriangEig"))
    const Grid& g = U.Grid();
    DistMatrix<F,STAR,STAR> dMain(g), dSub(g), dSup(g);
    DistMatrix<Complex<Base<F>>,STAR,STAR> w_STAR_STAR(g);
    dMain = U.GetDiagonal();
    dSub = U.GetDiagonal(-1);
    dSup = U.GetDiagonal(+1);
    w_STAR_STAR.Resize( U.Height(), 1 );
    QuasiTriangEig
    ( dMain.Matrix(), dSub.Matrix(), dSup.Matrix(), w_STAR_STAR.Matrix() );
    w = w_STAR_STAR;
}

template<typename F>
DistMatrix<Complex<Base<F>>,VR,STAR> QuasiTriangEig( const DistMatrix<F>& U )
{
    DEBUG_ONLY(CallStackEntry cse("schur::QuasiTriangEig"))
    DistMatrix<Complex<Base<F>>,VR,STAR> w(U.Grid());
    QuasiTriangEig( U, w );
    return w;
}

} // namespace schur
} // namespace El

#endif // ifndef EL_SCHUR_QUASITRIANGEIG_HPP
