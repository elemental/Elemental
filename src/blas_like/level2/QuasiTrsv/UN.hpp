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
void UNUnb( const Matrix<F>& U, Matrix<F>& x, bool checkIfSingular=false )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( U.Height() != U.Width() )
          LogicError("U must be square");
      if( x.Width() != 1 && x.Height() != 1 )
          LogicError("x must be a vector");
      const Int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
      if( U.Width() != xLength )
          LogicError("Nonconformal");
    )
    typedef Base<F> Real;

    F* xBuf = x.Buffer();
    const F* UBuf = U.LockedBuffer();
    const Int incx = ( x.Width()==1 ? 1 : x.LDim() );
    const Int ldu = U.LDim();
    const Int m = U.Height();
    Int k=m-1;
    while( k >= 0 )
    {
        const bool in2x2 = ( k>0 && UBuf[k+(k-1)*ldu] != F(0) );
        if( in2x2 )
        {
            --k;
            // Solve the 2x2 linear system via a 2x2 QR decomposition produced
            // by the Givens rotation
            //    | c        s | | U(k,  k) | = | gamma11 | 
            //    | -conj(s) c | | U(k+1,k) |   | 0       |
            //
            // and by also forming the right two entries of the 2x2 resulting
            // upper-triangular matrix, say gamma12 and gamma22
            //
            // Extract the 2x2 diagonal block, D
            const F delta11 = UBuf[ k   + k   *ldu];
            const F delta12 = UBuf[ k   +(k+1)*ldu];
            const F delta21 = UBuf[(k+1)+ k   *ldu];
            const F delta22 = UBuf[(k+1)+(k+1)*ldu];
            // Decompose D = Q R
            Real c; F s;
            const F gamma11 = Givens( delta11, delta21, c, s );
            const F gamma12 =        c*delta12 + s*delta22;
            const F gamma22 = -Conj(s)*delta12 + c*delta22;
            if( checkIfSingular )
            {
                // TODO: Instead check if values are too small in magnitude
                if( gamma11 == F(0) || gamma22 == F(0) )
                    LogicError("Singular diagonal block detected");
            }
            // Solve against Q
            const F chi1 = xBuf[ k   *incx];
            const F chi2 = xBuf[(k+1)*incx];
            xBuf[ k   *incx] =        c*chi1 + s*chi2;
            xBuf[(k+1)*incx] = -Conj(s)*chi1 + c*chi2;
            // Solve against R
            xBuf[(k+1)*incx] /= gamma22;
            xBuf[ k   *incx] -= gamma12*xBuf[(k+1)*incx];
            xBuf[ k   *incx] /= gamma11;

            // Update x0 := x0 - U01 x1
            blas::Axpy( k, -xBuf[ k   *incx], &UBuf[ k   *ldu], 1, xBuf, incx );
            blas::Axpy( k, -xBuf[(k+1)*incx], &UBuf[(k+1)*ldu], 1, xBuf, incx );
        }
        else
        {
            if( checkIfSingular )
                if( UBuf[k+k*ldu] == F(0) )
                    LogicError("Singular diagonal entry detected");
            // Solve the 1x1 linear system
            xBuf[k*incx] /= UBuf[k+k*ldu];

            // Update x0 := x0 - u01 chi_1
            blas::Axpy( k, -xBuf[k*incx], &UBuf[k*ldu], 1, xBuf, incx );
        }
        --k;
    }
}

template<typename F>
void UN( const Matrix<F>& U, Matrix<F>& x, bool checkIfSingular=false )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( U.Height() != U.Width() )
          LogicError("U must be square");
      if( x.Width() != 1 && x.Height() != 1 )
          LogicError("x must be a vector");
      const Int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
      if( U.Width() != xLength )
          LogicError("Nonconformal");
    )
    const bool vert = ( x.Width()==1 );

    Matrix<F> x0, x1;
    const Int m = U.Height();
    const Int bsize = Blocksize();
    const Int kLast = LastOffset( m, bsize );
    Int k=kLast, kOld=m;
    while( true )
    {
        const bool in2x2 = ( k>0 && U.Get(k,k-1) != F(0) );
        if( in2x2 )
            --k;
        const IR ind0( 0, k ), ind1( k, kOld );

        auto U01 = LockedView( U, ind0, ind1 );
        auto U11 = LockedView( U, ind1, ind1 );

        if( vert )
        {
            View( x0, x, 0, 0, k,    1 );
            View( x1, x, k, 0, kOld, 1 );
        }
        else
        {
            View( x0, x, 0, 0, 1, k    );
            View( x1, x, 0, k, 1, kOld );
        }

        quasitrsv::UNUnb( U11, x1, checkIfSingular );
        Gemv( NORMAL, F(-1), U01, x1, F(1), x0 );

        if( k == 0 )
            break;
        kOld = k;
        k -= Min(bsize,k);
    }
}

template<typename F>
void UN
( const AbstractDistMatrix<F>& UPre,
        AbstractDistMatrix<F>& xPre,
  bool checkIfSingular=false )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      AssertSameGrids( UPre, xPre );
      if( UPre.Height() != UPre.Width() )
          LogicError("U must be square");
      if( xPre.Width() != 1 && xPre.Height() != 1 )
          LogicError("x must be a vector");
      const Int xLength = 
          ( xPre.Width() == 1 ? xPre.Height() : xPre.Width() );
      if( UPre.Width() != xLength )
          LogicError("Nonconformal");
    )
    const Int m = UPre.Height();
    const Int bsize = Blocksize();
    const Int kLast = LastOffset( m, bsize );
    const Grid& g = UPre.Grid();

    DistMatrixReadProxy<F,F,MC,MR> UProx( UPre );
    DistMatrixReadWriteProxy<F,F,MC,MR> xProx( xPre );
    auto& U = UProx.GetLocked();
    auto& x = xProx.Get();

    // Matrix views 
    DistMatrix<F> U01(g), U11(g), x1(g);

    // Temporary distributions
    DistMatrix<F,STAR,STAR> U11_STAR_STAR(g), x1_STAR_STAR(g);

    if( x.Width() == 1 )
    {
        DistMatrix<F,MR,STAR> x1_MR_STAR(g);
        DistMatrix<F,MC,STAR> z_MC_STAR(g);

        // Views of z[MC,* ], which will store updates to x
        DistMatrix<F,MC,STAR> z0_MC_STAR(g), z1_MC_STAR(g);

        z_MC_STAR.AlignWith( U );
        z_MC_STAR.Resize( m, 1 );
        Zero( z_MC_STAR );

        Int k=kLast, kOld=m;
        while( true )
        {
            const bool in2x2 = ( k>0 && U.Get(k,k-1) != F(0) );
            if( in2x2 )
                --k;
            const IR ind0( 0, k ), ind1( k, kOld );

            LockedView( U01, U, ind0, ind1 );
            LockedView( U11, U, ind1, ind1 );

            View( x1, x, ind1, ALL );

            View( z0_MC_STAR, z_MC_STAR, ind0, ALL );
            View( z1_MC_STAR, z_MC_STAR, ind1, ALL );

            if( kOld != m )
                AxpyContract( F(1), z1_MC_STAR, x1 );

            x1_STAR_STAR = x1;
            U11_STAR_STAR = U11;
            quasitrsv::UN
            ( U11_STAR_STAR.LockedMatrix(), x1_STAR_STAR.Matrix(), 
              checkIfSingular );
            x1 = x1_STAR_STAR;

            x1_MR_STAR.AlignWith( U01 );
            x1_MR_STAR = x1_STAR_STAR;
            LocalGemv( NORMAL, F(-1), U01, x1_MR_STAR, F(1), z0_MC_STAR );

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

        z_STAR_MC.AlignWith( U );
        z_STAR_MC.Resize( 1, m );
        Zero( z_STAR_MC );

        Int k=kLast, kOld=m;
        while( true )
        {
            const bool in2x2 = ( k>0 && U.Get(k,k-1) != F(0) );
            if( in2x2 )
                --k;
            const IR ind0( 0, k ), ind1( k, kOld );

            LockedView( U01, U, ind0, ind1 );
            LockedView( U11, U, ind1, ind1 );

            View( x1, x, ALL, ind1 );

            View( z0_STAR_MC, z_STAR_MC, ALL, ind0 );
            View( z1_STAR_MC, z_STAR_MC, ALL, ind1 );

            if( kOld != m )
            {
                Contract( z1_STAR_MC, z1_MR_MC );
                x1 += z1_MR_MC;
            }

            x1_STAR_STAR = x1;
            U11_STAR_STAR = U11;
            quasitrsv::UN
            ( U11_STAR_STAR.LockedMatrix(), x1_STAR_STAR.Matrix(),
              checkIfSingular );
            x1 = x1_STAR_STAR;

            x1_STAR_MR.AlignWith( U01 );
            x1_STAR_MR = x1_STAR_STAR;
            LocalGemv( NORMAL, F(-1), U01, x1_STAR_MR, F(1), z0_STAR_MC );

            if( k == 0 )
                break;
            kOld = k;
            k -= Min(bsize,k);
        }
    }
}

} // namespace quasitrsv
} // namespace El
