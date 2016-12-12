/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace El {
namespace symv {

// s += alpha A  q
// t += alpha A' r
template<typename T>
void FusedRowPanelGemvs
( bool conjugate, T alpha, 
  const Matrix<T>& A, const Matrix<T>& q, const Matrix<T>& r, 
                            Matrix<T>& s,       Matrix<T>& t,
  Int bsize )
{
    const Int m = A.Height();
    const Int n = A.Width();
    const char transChar = ( conjugate ? 'C' : 'T' );
    const T* ABuf = A.LockedBuffer();
    const T* qBuf = q.LockedBuffer();
    const T* rBuf = r.LockedBuffer();
          T* sBuf = s.Buffer();
          T* tBuf = t.Buffer();
    const Int ALDim = A.LDim();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(n-k,bsize);
        blas::Gemv
        ( 'N', m, nb, alpha, 
          &ABuf[k*ALDim], ALDim, &qBuf[k], 1, T(1), sBuf, 1 );
        blas::Gemv
        ( transChar, m, nb, alpha,
          &ABuf[k*ALDim], ALDim, rBuf,     1, T(1), &tBuf[k], 1 );
    }
}

template<typename T>
void LocalColAccumulateUGeneral
( T alpha, 
  const DistMatrix<T>& A,
  const DistMatrix<T,MC,STAR>& x_MC_STAR,
  const DistMatrix<T,MR,STAR>& x_MR_STAR,
        DistMatrix<T,MC,STAR>& z_MC_STAR,
        DistMatrix<T,MR,STAR>& z_MR_STAR,
  bool conjugate, const SymvCtrl<T>& ctrl )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      AssertSameGrids( A, x_MC_STAR, x_MR_STAR, z_MC_STAR, z_MR_STAR );
      if( x_MC_STAR.Width() != 1 || x_MR_STAR.Width() != 1 ||
          z_MC_STAR.Width() != 1 || z_MR_STAR.Width() != 1 )
          LogicError("Expected x and z to be column vectors");
      if( A.Height() != A.Width() || 
          A.Height() != x_MC_STAR.Height() ||
          A.Height() != x_MR_STAR.Height() ||
          A.Height() != z_MC_STAR.Height() ||
          A.Height() != z_MR_STAR.Height() )
          LogicError
          ("Nonconformal: \n",
           "  A ~ ",A.Height()," x ",A.Width(),"\n",
           "  x[MC,* ] ~ ",x_MC_STAR.Height()," x ",x_MC_STAR.Width(),"\n",
           "  x[MR,* ] ~ ",x_MR_STAR.Height()," x ",x_MR_STAR.Width(),"\n",
           "  z[MC,* ] ~ ",z_MC_STAR.Height()," x ",z_MC_STAR.Width(),"\n",
           "  z[MR,* ] ~ ",z_MR_STAR.Height()," x ",z_MR_STAR.Width(),"\n");
      if( x_MC_STAR.ColAlign() != A.ColAlign() ||
          x_MR_STAR.ColAlign() != A.RowAlign() ||
          z_MC_STAR.ColAlign() != A.ColAlign() ||
          z_MR_STAR.ColAlign() != A.RowAlign() )
          LogicError("Partial matrix distributions are misaligned");
    )
    const Grid& g = A.Grid();
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );

    DistMatrix<T> D11(g);

    // We want our local gemvs to be of width blocksize, so we will 
    // temporarily change to max(r,c) times the current blocksize
    const Int n = A.Height();
    const Int bsize = Max(g.Height(),g.Width())*ctrl.bsize;
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k); 
        const Range<Int> ind1( k, k+nb ), ind2( k+nb, n );
 
        auto A11 = A( ind1, ind1 );
        auto A12 = A( ind1, ind2 );
        auto x1_MC_STAR = x_MC_STAR( ind1, ALL );
        auto x1_MR_STAR = x_MR_STAR( ind1, ALL );
        auto x2_MR_STAR = x_MR_STAR( ind2, ALL );
        auto z1_MC_STAR = z_MC_STAR( ind1, ALL );
        auto z1_MR_STAR = z_MR_STAR( ind1, ALL );
        auto z2_MR_STAR = z_MR_STAR( ind2, ALL );

        D11.AlignWith( A11 );
        // TODO: These diagonal block updates can be greatly improved
        D11 = A11;
        MakeTrapezoidal( UPPER, D11 );
        LocalGemv( NORMAL,      alpha, D11, x1_MR_STAR, T(1), z1_MC_STAR );
        FillDiagonal( D11, T(0) );
        LocalGemv( orientation, alpha, D11, x1_MC_STAR, T(1), z1_MR_STAR );
        
        // TODO: Expose the fusion blocksize as another parameter
        FusedRowPanelGemvs
        ( conjugate, alpha, A12.LockedMatrix(),
          x2_MR_STAR.LockedMatrix(), x1_MC_STAR.LockedMatrix(),
          z1_MC_STAR.Matrix(),       z2_MR_STAR.Matrix(), ctrl.bsize );
        // NOTE: The following are mathematically equivalent but should 
        //       be slower in practice for cache reasons
        //LocalGemv( NORMAL,      alpha, A12, x2_MR_STAR, T(1), z1_MC_STAR );
        //LocalGemv( orientation, alpha, A12, x1_MC_STAR, T(1), z2_MR_STAR );
    }
}

template<typename T>
void LocalColAccumulateUSquareTwoTrmv
( T alpha, 
  const DistMatrix<T>& A,
  const DistMatrix<T,MC,STAR>& x_MC_STAR,
  const DistMatrix<T,MR,STAR>& x_MR_STAR,
        DistMatrix<T,MC,STAR>& z_MC_STAR,
        DistMatrix<T,MR,STAR>& z_MR_STAR,
  bool conjugate )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      AssertSameGrids( A, x_MC_STAR, x_MR_STAR, z_MC_STAR, z_MR_STAR );
      if( x_MC_STAR.Width() != 1 || x_MR_STAR.Width() != 1 ||
          z_MC_STAR.Width() != 1 || z_MR_STAR.Width() != 1 )
          LogicError("Expected x and z to be column vectors");
      if( A.Height() != A.Width() || 
          A.Height() != x_MC_STAR.Height() ||
          A.Height() != x_MR_STAR.Height() ||
          A.Height() != z_MC_STAR.Height() ||
          A.Height() != z_MR_STAR.Height() )
          LogicError
          ("Nonconformal: \n",
           "  A ~ ",A.Height()," x ",A.Width(),"\n",
           "  x[MC,* ] ~ ",x_MC_STAR.Height()," x ",x_MC_STAR.Width(),"\n",
           "  x[MR,* ] ~ ",x_MR_STAR.Height()," x ",x_MR_STAR.Width(),"\n",
           "  z[MC,* ] ~ ",z_MC_STAR.Height()," x ",z_MC_STAR.Width(),"\n",
           "  z[MR,* ] ~ ",z_MR_STAR.Height()," x ",z_MR_STAR.Width(),"\n");
      if( x_MC_STAR.ColAlign() != A.ColAlign() ||
          x_MR_STAR.ColAlign() != A.RowAlign() ||
          z_MC_STAR.ColAlign() != A.ColAlign() ||
          z_MR_STAR.ColAlign() != A.RowAlign() )
          LogicError("Partial matrix distributions are misaligned");
      if( A.Grid().Height() != A.Grid().Width() )
          LogicError("Grid must be square");
    )

    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    if( A.ColShift() < A.RowShift() )
    {
        // We are above the diagonal, so we can multiply without an 
        // offset for triu(A)[MC,MR] and triu(A,1)'[MR,MC]
        if( localHeight != 0 )
        {
            // Our local portion of z[MC,* ] might be one entry longer
            // than A.LocalWidth(), so go ahead and set the last entry
            // to 0.
            z_MC_STAR.SetLocal(localHeight-1,0,T(0));
        }
        MemCopy
        ( z_MC_STAR.Buffer(),
          x_MR_STAR.LockedBuffer(), localWidth );
        MemCopy
        ( z_MR_STAR.Buffer(),
          x_MC_STAR.LockedBuffer(), localWidth );

        blas::Trmv
        ( 'U', 'N', 'N', localWidth, 
          A.LockedBuffer(),   A.LDim(),
          z_MC_STAR.Buffer(), 1 );
        blas::Trmv
        ( 'U', 'C', 'N', localWidth, 
          A.LockedBuffer(),   A.LDim(),
          z_MR_STAR.Buffer(), 1 );
    }
    else if( A.ColShift() > A.RowShift() )
    {
        // We are below the diagonal, so we need to use an offset 
        // for both triu(A)[MC,MR] and triu(A,+1)'[MR,MC]
        if( localHeight != 0 )
        {
            // The last entry of z[MC,* ] will be zero due to the forced
            // offset
            z_MC_STAR.SetLocal(localHeight-1,0,T(0));
        }
        if( localWidth != 0 )
        {
            MemCopy
            ( z_MC_STAR.Buffer(),
              x_MR_STAR.LockedBuffer(1,0), localWidth-1 );
            // The first entry of z[MR,* ] will be zero due to the forced
            // offset
            z_MR_STAR.SetLocal(0,0,T(0));
            MemCopy
            ( z_MR_STAR.Buffer(1,0),
              x_MC_STAR.LockedBuffer(), localWidth-1 );

            blas::Trmv
            ( 'U', 'N', 'N', localWidth-1,
              A.LockedBuffer(0,1), A.LDim(),
              z_MC_STAR.Buffer(),  1 );
            blas::Trmv
            ( 'U', 'C', 'N', localWidth-1,
              A.LockedBuffer(0,1),   A.LDim(),
              z_MR_STAR.Buffer(1,0), 1 );
        }
    }
    else
    {
        // We are on the diagonal, so we only need an offset for
        // triu(A,+1)'[MR,MC]
        if( localWidth != 0 )
        {
            MemCopy
            ( z_MC_STAR.Buffer(),
              x_MR_STAR.LockedBuffer(), localHeight );
            // The first entry of z[MR,* ] must be zero due to the forced
            // offset
            z_MR_STAR.SetLocal(0,0,T(0));
            MemCopy
            ( z_MR_STAR.Buffer(1,0),
              x_MC_STAR.LockedBuffer(), localWidth-1 );

            blas::Trmv
            ( 'U', 'N', 'N', localHeight, 
              A.LockedBuffer(),   A.LDim(),
              z_MC_STAR.Buffer(), 1 );
            blas::Trmv
            ( 'U', 'C', 'N', localWidth-1,
              A.LockedBuffer(0,1),   A.LDim(),
              z_MR_STAR.Buffer(1,0), 1 );
        }
    }
    z_MC_STAR *= alpha;
    z_MR_STAR *= alpha;
}

template<typename T>
void LocalColAccumulateU
( T alpha, 
  const DistMatrix<T>& A,
  const DistMatrix<T,MC,STAR>& x_MC_STAR,
  const DistMatrix<T,MR,STAR>& x_MR_STAR,
        DistMatrix<T,MC,STAR>& z_MC_STAR,
        DistMatrix<T,MR,STAR>& z_MR_STAR,
  bool conjugate, const SymvCtrl<T>& ctrl )
{
    EL_DEBUG_CSE
    if( ctrl.avoidTrmvBasedLocalSymv || A.Grid().Height() != A.Grid().Width() )
        LocalColAccumulateUGeneral
        ( alpha, A, x_MC_STAR, x_MR_STAR, z_MC_STAR, z_MR_STAR, conjugate,
          ctrl );
    else
        LocalColAccumulateUSquareTwoTrmv
        ( alpha, A, x_MC_STAR, x_MR_STAR, z_MC_STAR, z_MR_STAR, conjugate );
}

template<typename T>
void LocalRowAccumulateU
( T alpha, 
  const DistMatrix<T>& A,
  const DistMatrix<T,STAR,MC>& x_STAR_MC,
  const DistMatrix<T,STAR,MR>& x_STAR_MR,
        DistMatrix<T,STAR,MC>& z_STAR_MC,
        DistMatrix<T,STAR,MR>& z_STAR_MR,
  bool conjugate, const SymvCtrl<T>& ctrl )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      AssertSameGrids( A, x_STAR_MC, x_STAR_MR, z_STAR_MC, z_STAR_MR );
      if( x_STAR_MC.Height() != 1 || x_STAR_MR.Height() != 1 ||
          z_STAR_MC.Height() != 1 || z_STAR_MR.Height() != 1 )
          LogicError("Expected x and z to be row vectors");
      if( A.Height() != A.Width() || 
          A.Height() != x_STAR_MC.Width() ||
          A.Height() != x_STAR_MR.Width() ||
          A.Height() != z_STAR_MC.Width() ||
          A.Height() != z_STAR_MR.Width() )
          LogicError
          ("Nonconformal: \n",
           "  A ~ ",A.Height()," x ",A.Width(),"\n",
           "  x[* ,MC] ~ ",x_STAR_MC.Height()," x ",x_STAR_MC.Width(),"\n",
           "  x[* ,MR] ~ ",x_STAR_MR.Height()," x ",x_STAR_MR.Width(),"\n",
           "  z[* ,MC] ~ ",z_STAR_MC.Height()," x ",z_STAR_MC.Width(),"\n",
           "  z[* ,MR] ~ ",z_STAR_MR.Height()," x ",z_STAR_MR.Width(),"\n");
      if( x_STAR_MC.RowAlign() != A.ColAlign() ||
          x_STAR_MR.RowAlign() != A.RowAlign() ||
          z_STAR_MC.RowAlign() != A.ColAlign() ||
          z_STAR_MR.RowAlign() != A.RowAlign() )
          LogicError("Partial matrix distributions are misaligned");
    )
    const Grid& g = A.Grid();
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );

    DistMatrix<T> D11(g);

    // We want our local gemvs to be of width blocksize, so we will 
    // temporarily change to max(r,c) times the current blocksize
    const Int n = A.Height();
    const Int bsize = Max(g.Height(),g.Width())*ctrl.bsize;
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);
        const Range<Int> ind1( k, k+nb ), ind2( k+nb, n );

        auto A11 = A( ind1, ind1 );
        auto A12 = A( ind1, ind2 );
        auto x1_STAR_MC = x_STAR_MC( ALL, ind1 );
        auto x1_STAR_MR = x_STAR_MR( ALL, ind1 );
        auto x2_STAR_MR = x_STAR_MR( ALL, ind2 );
        auto z1_STAR_MC = z_STAR_MC( ALL, ind1 );
        auto z1_STAR_MR = z_STAR_MR( ALL, ind1 );
        auto z2_STAR_MR = z_STAR_MR( ALL, ind2 );

        D11.AlignWith( A11 );
        // TODO: These diagonal block updates can be greatly improved
        D11 = A11;
        MakeTrapezoidal( UPPER, D11 );
        LocalGemv( NORMAL,      alpha, D11, x1_STAR_MR, T(1), z1_STAR_MC );
        FillDiagonal( D11, T(0) );
        LocalGemv( orientation, alpha, D11, x1_STAR_MC, T(1), z1_STAR_MR );

        LocalGemv( NORMAL,      alpha, A12, x2_STAR_MR, T(1), z1_STAR_MC );
        LocalGemv( orientation, alpha, A12, x1_STAR_MC, T(1), z2_STAR_MR );
    }
}

} // namespace symv
} // namespace El
