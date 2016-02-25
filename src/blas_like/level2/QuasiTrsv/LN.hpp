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
inline void
LNUnb( const Matrix<F>& L, Matrix<F>& x, bool checkIfSingular=false )
{
    DEBUG_ONLY(
      CSE cse("quasitrsv::LNUnb");
      if( L.Height() != L.Width() )
          LogicError("L must be square");
      if( x.Width() != 1 && x.Height() != 1 )
          LogicError("x must be a vector");
      const Int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
      if( L.Width() != xLength )
          LogicError("Nonconformal");
    )
    typedef Base<F> Real;

    F* xBuf = x.Buffer();
    const F* LBuf = L.LockedBuffer();
    const Int incx = ( x.Width()==1 ? 1 : x.LDim() );
    const Int ldl = L.LDim();
    const Int m = L.Height();
    for( Int k=0; k<m; ++k )
    {
        const bool in2x2 = ( k+1<m && LBuf[k+(k+1)*ldl] != F(0) );
        if( in2x2 )
        {
            // Solve the 2x2 linear system via a 2x2 LQ decomposition produced
            // by the Givens rotation
            //    | L(k,k) L(k,k+1) | | c -conj(s) | = | gamma11 0 |
            //                        | s    c     |
            // and by also forming the bottom two entries of the 2x2 resulting
            // lower-triangular matrix, say gamma21 and gamma22
            //
            // Extract the 2x2 diagonal block, D
            const F delta11 = LBuf[ k   + k   *ldl];
            const F delta12 = LBuf[ k   +(k+1)*ldl];
            const F delta21 = LBuf[(k+1)+ k   *ldl];
            const F delta22 = LBuf[(k+1)+(k+1)*ldl];
            // Decompose D = L Q
            Real c; F s;
            const F gamma11 = lapack::Givens( delta11, delta12, &c, &s );
            const F gamma21 =        c*delta21 + s*delta22;
            const F gamma22 = -Conj(s)*delta21 + c*delta22; 
            if( checkIfSingular )
            {
                // TODO: Instead check if values are too small in magnitude
                if( gamma11 == F(0) || gamma22 == F(0) ) 
                    LogicError("Singular diagonal block detected");
            }
            // Solve against L
            xBuf[ k   *incx] /= gamma11;
            xBuf[(k+1)*incx] -= gamma21*xBuf[k*incx];
            xBuf[(k+1)*incx] /= gamma22;
            // Solve against Q
            const F chi1 = xBuf[ k   *incx];
            const F chi2 = xBuf[(k+1)*incx];
            xBuf[ k   *incx] = c*chi1 - Conj(s)*chi2;
            xBuf[(k+1)*incx] = s*chi1 +       c*chi2;

            // Update x2 := x2 - L21 x1
            blas::Axpy
            ( m-(k+2), -xBuf[ k   *incx], 
              &LBuf[(k+2)+ k   *ldl], 1, &xBuf[(k+2)*incx], incx );
            blas::Axpy
            ( m-(k+2), -xBuf[(k+1)*incx], 
              &LBuf[(k+2)+(k+1)*ldl], 1, &xBuf[(k+2)*incx], incx );

            k += 2;
        }
        else
        {
            if( checkIfSingular )
                if( LBuf[k+k*ldl] == F(0) )
                    LogicError("Singular diagonal entry detected");
            // Solve the 1x1 linear system
            xBuf[k] /= LBuf[k+k*ldl];

            // Update x2 := x2 - l21 chi_1
            blas::Axpy
            ( m-(k+1), -xBuf[k*incx], 
              &LBuf[(k+1)+k*ldl], 1, &xBuf[(k+1)*incx], incx );

            k += 1;
        }
    }
}

template<typename F>
inline void
LN( const Matrix<F>& L, Matrix<F>& x, bool checkIfSingular=false )
{
    DEBUG_ONLY(
      CSE cse("quasitrsv::LN");
      if( L.Height() != L.Width() )
          LogicError("L must be square");
      if( x.Width() != 1 && x.Height() != 1 )
          LogicError("x must be a vector");
      const Int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
      if( L.Width() != xLength )
          LogicError("Nonconformal");
    )
    const bool vert = ( x.Width()==1 );

    Matrix<F> x1, x2;
    const Int m = L.Height();
    const Int bsize = Blocksize();
    Int k=0;
    while( k < m )
    {
        const Int nbProp = Min(bsize,m-k);
        const bool in2x2 = ( k+nbProp<m && L.Get(k+nbProp-1,k+nbProp) != F(0) );
        const Int nb = ( in2x2 ? nbProp+1 : nbProp );
        const IR ind1( k, k+nb ), ind2( k+nb, END );

        auto L11 = LockedView( L, ind1, ind1 );
        auto L21 = LockedView( L, ind2, ind1 );

        if( vert )
        {
            View( x1, x, ind1, ALL );
            View( x2, x, ind2, ALL );
        }
        else
        {
            View( x1, x, ALL, ind1 );
            View( x2, x, ALL, ind2 );
        }

        quasitrsv::LNUnb( L11, x1, checkIfSingular );
        Gemv( NORMAL, F(-1), L21, x1, F(1), x2 );

        k += nb;
    }
}

template<typename F>
inline void
LN
( const ElementalMatrix<F>& LPre,
        ElementalMatrix<F>& xPre, 
  bool checkIfSingular=false )
{
    DEBUG_ONLY(
      CSE cse("quasitrsv::LN");
      AssertSameGrids( LPre, xPre );
      if( LPre.Height() != LPre.Width() )
          LogicError("L must be square");
      if( xPre.Width() != 1 && xPre.Height() != 1 )
          LogicError("x must be a vector");
      const Int xLength = 
          ( xPre.Width() == 1 ? xPre.Height() : xPre.Width() );
      if( LPre.Width() != xLength )
          LogicError("Nonconformal");
    )
    const Int m = LPre.Height();
    const Int bsize = Blocksize();
    const Grid& g = LPre.Grid();

    DistMatrixReadProxy<F,F,MC,MR> LProx( LPre );
    DistMatrixReadWriteProxy<F,F,MC,MR> xProx( xPre );
    auto& L = LProx.GetLocked();
    auto& x = xProx.Get();

    // Matrix views 
    DistMatrix<F> L11(g), L21(g), x1(g);

    // Temporary distributions
    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g), x1_STAR_STAR(g);

    if( x.Width() == 1 )
    {
        DistMatrix<F,MR,STAR> x1_MR_STAR(g);
        DistMatrix<F,MC,STAR> z_MC_STAR(g);

        // Views of z[MC,* ], which will store updates to x
        DistMatrix<F,MC,STAR> z1_MC_STAR(g), z2_MC_STAR(g);

        z_MC_STAR.AlignWith( L );
        Zeros( z_MC_STAR, m, 1 );

        Int k=0;
        while( k < m )
        {
            const Int nbProp = Min(bsize,m-k);
            const bool in2x2 = 
                ( k+nbProp<m && L.Get(k+nbProp-1,k+nbProp) != F(0) );
            const Int nb = ( in2x2 ? nbProp+1 : nbProp );

            LockedView( L11, L, k,    k, k+nb, k+nb );
            LockedView( L21, L, k+nb, k, m,    k+nb );

            View( x1, x, k,    0, k+nb, 1 );

            View( z1_MC_STAR, z_MC_STAR, k,    0, k+nb, 1 );
            View( z2_MC_STAR, z_MC_STAR, k+nb, 0, m,    1 ); 

            if( k != 0 )
                AxpyContract( F(1), z1_MC_STAR, x1 );

            x1_STAR_STAR = x1;
            L11_STAR_STAR = L11;
            quasitrsv::LN
            ( L11_STAR_STAR.LockedMatrix(), x1_STAR_STAR.Matrix(), 
              checkIfSingular );
            x1 = x1_STAR_STAR;

            x1_MR_STAR.AlignWith( L21 );
            x1_MR_STAR = x1_STAR_STAR;
            LocalGemv( NORMAL, F(-1), L21, x1_MR_STAR, F(1), z2_MC_STAR );

            k += nb;
        }
    }
    else
    {
        DistMatrix<F,STAR,MR> x1_STAR_MR(g);
        DistMatrix<F,MR,  MC> z1_MR_MC(g);
        DistMatrix<F,STAR,MC> z_STAR_MC(g);

        // Views of z[* ,MC]
        DistMatrix<F,STAR,MC> z1_STAR_MC(g), z2_STAR_MC(g);

        z_STAR_MC.AlignWith( L );
        Zeros( z_STAR_MC, 1, m );

        Int k=0;
        while( k < m )
        {
            const Int nbProp = Min(bsize,m-k);
            const bool in2x2 = 
                ( k+nbProp<m && L.Get(k+nbProp-1,k+nbProp) != F(0) );
            const Int nb = ( in2x2 ? nbProp+1 : nbProp );
            const IR ind1( k, k+nb ), ind2( k+nb, END );

            LockedView( L11, L, ind1, ind1 );
            LockedView( L21, L, ind2, ind1 );

            View( x1, x, ALL, ind1 );

            View( z1_STAR_MC, z_STAR_MC, ALL, ind1 );
            View( z2_STAR_MC, z_STAR_MC, ALL, ind2 );

            if( k != 0 )
            {
                Contract( z1_STAR_MC, z1_MR_MC );
                x1 += z1_MR_MC;
            }

            x1_STAR_STAR = x1;
            L11_STAR_STAR = L11;
            quasitrsv::LN
            ( L11_STAR_STAR.LockedMatrix(), x1_STAR_STAR.Matrix(), 
              checkIfSingular );
            x1 = x1_STAR_STAR;

            x1_STAR_MR.AlignWith( L21 );
            x1_STAR_MR = x1_STAR_STAR;
            LocalGemv( NORMAL, F(-1), L21, x1_STAR_MR, F(1), z2_STAR_MC );

            k += nb;
        }
    }
}

} // namespace quasitrsv
} // namespace El
