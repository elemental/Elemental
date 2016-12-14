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
    EL_DEBUG_CSE
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
            auto w1 = w( IR(j,j+2), ALL );
            HessenbergSchur( H11, w1 );
            j += 2;
        }
    }
}

template<typename F>
void QuasiTriangEig( const Matrix<F>& U, Matrix<Complex<Base<F>>>& w )
{
    EL_DEBUG_CSE
    auto dMain = GetDiagonal(U);
    auto dSub = GetDiagonal(U,-1);
    auto dSup = GetDiagonal(U,+1);
    QuasiTriangEig( dMain, dSub, dSup, w );
}

template<typename F>
Matrix<Complex<Base<F>>> QuasiTriangEig( const Matrix<F>& U )
{
    EL_DEBUG_CSE
    Matrix<Complex<Base<F>>> w;
    QuasiTriangEig( U, w );
    return w;
}

template<typename F>
void QuasiTriangEig
( const AbstractDistMatrix<F>& UPre,
        AbstractDistMatrix<Complex<Base<F>>>& wPre )
{
    EL_DEBUG_CSE
    typedef Base<F> Real;

    DistMatrixReadProxy<F,F,MC,MR> UProx( UPre );
    auto& U = UProx.GetLocked();

    DistMatrixWriteProxy<Complex<Real>,Complex<Real>,STAR,STAR> wProx( wPre );
    auto& w = wProx.Get();

    const Grid& g = U.Grid();
    DistMatrix<F,STAR,STAR> dMain(g), dSub(g), dSup(g);
    dMain = GetDiagonal(U);
    dSub = GetDiagonal(U,-1);
    dSup = GetDiagonal(U,+1);
    w.Resize( U.Height(), 1 );
    QuasiTriangEig( dMain.Matrix(), dSub.Matrix(), dSup.Matrix(), w.Matrix() );
}

template<typename F>
DistMatrix<Complex<Base<F>>,VR,STAR> 
QuasiTriangEig( const AbstractDistMatrix<F>& U )
{
    EL_DEBUG_CSE
    DistMatrix<Complex<Base<F>>,VR,STAR> w(U.Grid());
    QuasiTriangEig( U, w );
    return w;
}

} // namespace schur
} // namespace El

#endif // ifndef EL_SCHUR_QUASITRIANGEIG_HPP
