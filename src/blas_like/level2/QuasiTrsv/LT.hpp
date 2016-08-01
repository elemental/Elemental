/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace El {
namespace quasitrsv {

template<typename F>
void LTUnb
( Orientation orientation, const Matrix<F>& L, Matrix<F>& x, 
  bool checkIfSingular=false )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( L.Height() != L.Width() )
          LogicError("L must be square");
      if( x.Width() != 1 && x.Height() != 1 )
          LogicError("x must be a vector");
      const Int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
      if( L.Width() != xLength )
          LogicError("Nonconformal");
      if( orientation == NORMAL )
          LogicError("Invalid orientation");
    )
    typedef Base<F> Real;
    const bool conjugate = ( orientation == ADJOINT );
    if( conjugate )
        Conjugate( x );

    F* xBuf = x.Buffer();
    const F* LBuf = L.LockedBuffer();
    const Int incx = ( x.Width()==1 ? 1 : x.LDim() );
    const Int ldl = L.LDim();
    const Int m = L.Height();
    Int k=m-1;
    while( k >= 0 )
    {
        const bool in2x2 = ( k>0 && LBuf[(k-1)+k*ldl] != F(0) );
        if( in2x2 )
        {
            --k;
            // Solve the 2x2 linear system via a 2x2 LQ decomposition produced
            // by the Givens rotation
            //    | L(k,k) L(k,k+1) | | c -conj(s) | = | gamma11 0 |
            //                        | s    c     |
            // and by also forming the bottom two entries of the 2x2 resulting
            // lower-triangular matrix, say gamma21 and gamma22
            //
            // Extract the 2x2 diagonal block, D
            const F delta11 = LBuf[   k +   k *ldl];
            const F delta12 = LBuf[   k +(k+1)*ldl];
            const F delta21 = LBuf[(k+1)+   k *ldl];
            const F delta22 = LBuf[(k+1)+(k+1)*ldl];
            // Decompose D = L Q
            Real c; F s;
            const F gamma11 = Givens( delta11, delta12, c, s );
            const F gamma21 =        c*delta21 + s*delta22;
            const F gamma22 = -Conj(s)*delta21 + c*delta22;
            if( checkIfSingular )
            {
                // TODO: Instead check if values are too small in magnitude
                if( gamma11 == F(0) || gamma22 == F(0) )
                    LogicError("Singular diagonal block detected");
            }
            // Solve against Q^T
            const F chi1 = xBuf[ k   *incx];
            const F chi2 = xBuf[(k+1)*incx];
            xBuf[ k   *incx] =        c*chi1 + s*chi2;
            xBuf[(k+1)*incx] = -Conj(s)*chi1 + c*chi2;
            // Solve against R^T
            xBuf[(k+1)*incx] /= gamma22;
            xBuf[ k   *incx] -= gamma21*xBuf[(k+1)*incx];
            xBuf[ k   *incx] /= gamma11;

            // Update x0 := x0 - L10^T x1
            blas::Axpy( k, -xBuf[ k   *incx], &LBuf[k  ], ldl, xBuf, incx );
            blas::Axpy( k, -xBuf[(k+1)*incx], &LBuf[k+1], ldl, xBuf, incx );
        }
        else
        {
            if( checkIfSingular )
                if( LBuf[k+k*ldl] == F(0) )
                    LogicError("Singular diagonal entry detected");
            // Solve the 1x1 linear system
            xBuf[k*incx] /= LBuf[k+k*ldl];

            // Update x0 := x0 - l10^T chi_1
            blas::Axpy( k, -xBuf[k*incx], &LBuf[k], ldl, xBuf, incx );
        }
        --k;
    }
    if( conjugate )
        Conjugate( x );
}

template<typename F>
void LT
( Orientation orientation, const Matrix<F>& L, Matrix<F>& x, 
  bool checkIfSingular=false )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( L.Height() != L.Width() )
          LogicError("L must be square");
      if( x.Width() != 1 && x.Height() != 1 )
          LogicError("x must be a vector");
      const Int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
      if( L.Width() != xLength )
          LogicError("Nonconformal");
      if( orientation == NORMAL )
          LogicError("Invalid orientation");
    )
    const bool vert = ( x.Width()==1 );
    const bool conjugate = ( orientation==ADJOINT );
    if( conjugate )
        Conjugate( x );

    Matrix<F> x0, x1;
    const Int m = L.Height();
    const Int bsize = Blocksize();
    const Int kLast = LastOffset( m, bsize );
    Int k=kLast, kOld=m;
    while( true )
    {
        const bool in2x2 = ( k>0 && L.Get(k-1,k) != F(0) );
        if( in2x2 )
            --k;
        const IR ind0( 0, k ), ind1( k, kOld );

        auto L10 = LockedView( L, ind1, ind0 );
        auto L11 = LockedView( L, ind1, ind1 );

        if( vert )
        {
            View( x0, x, ind0, ALL );
            View( x1, x, ind1, ALL );
        }
        else
        {
            View( x0, x, ALL, ind0 );
            View( x1, x, ALL, ind1 );
        }

        quasitrsv::LTUnb( TRANSPOSE, L11, x1, checkIfSingular );
        Gemv( TRANSPOSE, F(-1), L10, x1, F(1), x0 );

        if( k == 0 )
            break;
        kOld = k;
        k -= Min(bsize,k);
    }
    if( conjugate )
        Conjugate( x );
}

template<typename F>
void LT
( Orientation orientation, 
  const AbstractDistMatrix<F>& LPre,
        AbstractDistMatrix<F>& xPre,
  bool checkIfSingular=false )
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertSameGrids( LPre, xPre );
      if( LPre.Height() != LPre.Width() )
          LogicError("L must be square");
      if( xPre.Width() != 1 && xPre.Height() != 1 )
          LogicError("x must be a vector");
      const Int xLength = 
          ( xPre.Width() == 1 ? xPre.Height() : xPre.Width() );
      if( LPre.Width() != xLength )
          LogicError("Nonconformal");
      if( orientation == NORMAL )
          LogicError("Invalid orientation");
    )
    const Int m = LPre.Height();
    const Int bsize = Blocksize();
    const Int kLast = LastOffset( m, bsize );
    const Grid& g = LPre.Grid();
    const bool conjugate = ( orientation==ADJOINT );
    if( conjugate )
        Conjugate( xPre );

    DistMatrixReadProxy<F,F,MC,MR> LProx( LPre );
    DistMatrixReadWriteProxy<F,F,MC,MR> xProx( xPre );
    auto& L = LProx.GetLocked();
    auto& x = xProx.Get();

    // Matrix views 
    DistMatrix<F> L10(g), L11(g), x1(g);

    // Temporary distributions
    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g), x1_STAR_STAR(g);

    if( x.Width() == 1 )
    {
        DistMatrix<F,MR,STAR> x1_MR_STAR(g);
        DistMatrix<F,MC,STAR> z_MC_STAR(g);

        // Views of z[MC,* ], which will store updates to x
        DistMatrix<F,MC,STAR> z0_MC_STAR(g), z1_MC_STAR(g);

        z_MC_STAR.AlignWith( L );
        z_MC_STAR.Resize( m, 1 );
        Zero( z_MC_STAR );

        Int k=kLast, kOld=m;
        while( true )
        {
            const bool in2x2 = ( k>0 && L.Get(k-1,k) != F(0) );
            if( in2x2 )
                --k;
            const IR ind0( 0, k ), ind1( k, kOld );

            LockedView( L10, L, ind1, ind0 );
            LockedView( L11, L, ind1, ind1 );

            View( x1, x, ind1, ALL );

            View( z0_MC_STAR, z_MC_STAR, ind0, ALL );
            View( z1_MC_STAR, z_MC_STAR, ind1, ALL );

            if( kOld != m )
                AxpyContract( F(1), z1_MC_STAR, x1 );

            x1_STAR_STAR = x1;
            L11_STAR_STAR = L11;
            quasitrsv::LT
            ( TRANSPOSE, L11_STAR_STAR.LockedMatrix(), x1_STAR_STAR.Matrix(),
              checkIfSingular );
            x1 = x1_STAR_STAR;

            x1_MR_STAR.AlignWith( L10 );
            x1_MR_STAR = x1_STAR_STAR;
            LocalGemv( TRANSPOSE, F(-1), L10, x1_MR_STAR, F(1), z0_MC_STAR );

            if( k == 0 )
                break;
            kOld = k;
            k -= Min(bsize,k);
        }
    }
    else
    {
        DistMatrix<F,STAR,MR  > x1_STAR_MR(g);
        DistMatrix<F,MC,  MR  > z1(g);
        DistMatrix<F,MR,  MC  > z1_MR_MC(g);
        DistMatrix<F,STAR,MC  > z_STAR_MC(g);

        // Views of z[* ,MC]
        DistMatrix<F,STAR,MC> z0_STAR_MC(g), z1_STAR_MC(g);

        z_STAR_MC.AlignWith( L );
        z_STAR_MC.Resize( 1, m );
        Zero( z_STAR_MC );

        Int k=kLast, kOld=m;
        while( true )
        {
            const bool in2x2 = ( k>0 && L.Get(k-1,k) != F(0) );
            if( in2x2 )
                --k;
            const IR ind0( 0, k ), ind1( k, kOld );

            LockedView( L10, L, ind1, ind0 );
            LockedView( L11, L, ind1, ind1 );

            View( x1, x, ALL, ind1 );

            View( z0_STAR_MC, z_STAR_MC, ALL, ind0 );
            View( z1_STAR_MC, z_STAR_MC, ALL, ind1 );

            if( kOld != m )
            {
                Contract( z1_STAR_MC, z1_MR_MC );
                x1 += z1_MR_MC;
            }

            x1_STAR_STAR = x1;
            L11_STAR_STAR = L11;
            quasitrsv::LT
            ( TRANSPOSE, L11_STAR_STAR.LockedMatrix(), x1_STAR_STAR.Matrix(),
              checkIfSingular );
            x1 = x1_STAR_STAR;

            x1_STAR_MR.AlignWith( L10 );
            x1_STAR_MR = x1_STAR_STAR;
            LocalGemv( TRANSPOSE, F(-1), L10, x1_STAR_MR, F(1), z0_STAR_MC );

            if( k == 0 )
                break;
            kOld = k;
            k -= Min(bsize,k);
        }
    }
    if( conjugate )
        Conjugate( x );
}

} // namespace quasitrsv
} // namespace El
