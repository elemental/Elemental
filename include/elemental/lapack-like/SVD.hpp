/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/

#ifdef HAVE_FLA_BSVD
// TODO: Move this into a file at include/elemental/imports/flame.hpp
extern "C" {
typedef int FLA_Error;
FLA_Error FLA_Bsvd_v_opd_var1( int       min_m_n,
                               int       m_U,
                               int       m_V,
                               int       n_GH,
                               int       n_iter_max,
                               double*   buff_d, int inc_d,
                               double*   buff_e, int inc_e,
                               elem::dcomplex* buff_G, int rs_G, int cs_G,
                               elem::dcomplex* buff_H, int rs_H, int cs_H,
                               double*   buff_U, int rs_U, int cs_U,
                               double*   buff_V, int rs_V, int cs_V,
                               int       b_alg );
FLA_Error FLA_Bsvd_v_opz_var1( int       min_m_n,
                               int       m_U,
                               int       m_V,
                               int       n_GH,
                               int       n_iter_max,
                               double*   buff_d, int inc_d,
                               double*   buff_e, int inc_e,
                               elem::dcomplex* buff_G, int rs_G, int cs_G,
                               elem::dcomplex* buff_H, int rs_H, int cs_H,
                               elem::dcomplex* buff_U, int rs_U, int cs_U,
                               elem::dcomplex* buff_V, int rs_V, int cs_V,
                               int       b_alg );
}
#endif // HAVE_FLA_BSVD

namespace elem {

namespace svd {

template<typename R>
inline void
SimpleSVD
( DistMatrix<R>& A,
  DistMatrix<R,VR,STAR>& s,
  DistMatrix<R>& V )
{
#ifndef RELEASE
    PushCallStack("svd::SimpleSVD");
#endif
    const int m = A.Height();
    const int n = A.Width();
    const int k = std::min( m, n );
    const int offdiagonal = ( m>=n ? 1 : -1 );
    const char uplo = ( m>=n ? 'U' : 'L' );
    const Grid& grid = A.Grid();

    // Bidiagonalize A
    Bidiag( A );

    // Grab copies of the diagonal and sub/super-diagonal of A
    DistMatrix<R,MD,STAR> d_MD_STAR( grid ), 
                          e_MD_STAR( grid );
    A.GetDiagonal( d_MD_STAR );
    A.GetDiagonal( e_MD_STAR, offdiagonal );

    // In order to use serial QR kernels, we need the full bidiagonal matrix
    // on each process
    DistMatrix<R,STAR,STAR> d_STAR_STAR( d_MD_STAR ),
                            e_STAR_STAR( e_MD_STAR );

    // Initialize U and VTrans to the appropriate identity matrices.
    DistMatrix<R,VC,STAR> U_VC_STAR( grid );
    DistMatrix<R,STAR,VC> VTrans_STAR_VC( grid );
    U_VC_STAR.AlignWith( A );
    VTrans_STAR_VC.AlignWith( V );
    Identity( m, k, U_VC_STAR );
    Identity( k, n, VTrans_STAR_VC );

    // Compute the SVD of the bidiagonal matrix and accumulate the Givens
    // rotations into our local portion of U and VTrans
    Matrix<R>& ULocal = U_VC_STAR.LocalMatrix();
    Matrix<R>& VTransLocal = VTrans_STAR_VC.LocalMatrix();
    lapack::BidiagQRAlg
    ( uplo, k, VTransLocal.Width(), ULocal.Height(),
      d_STAR_STAR.LocalBuffer(), e_STAR_STAR.LocalBuffer(), 
      VTransLocal.Buffer(), VTransLocal.LDim(), 
      ULocal.Buffer(), ULocal.LDim() );

    // Make a copy of A (for the Householder vectors) and pull the necessary 
    // portions of U and VTrans into a standard matrix dist.
    DistMatrix<R> B( A );
    if( m >= n )
    {
        DistMatrix<R> AT( grid ),
                      AB( grid );
        DistMatrix<R,VC,STAR> UT_VC_STAR( grid ), 
                              UB_VC_STAR( grid );
        PartitionDown( A, AT,
                          AB, n );
        PartitionDown( U_VC_STAR, UT_VC_STAR,
                                  UB_VC_STAR, n );
        AT = UT_VC_STAR;
        MakeZeros( AB );
        Transpose( VTrans_STAR_VC, V );
    }
    else
    {
        DistMatrix<R> VT( grid ), 
                      VB( grid );
        DistMatrix<R,STAR,VC> VTransL_STAR_VC( grid ), VTransR_STAR_VC( grid );
        PartitionDown( V, VT, 
                          VB, m );
        PartitionRight( VTrans_STAR_VC, VTransL_STAR_VC, VTransR_STAR_VC, m );
        Transpose( VTransL_STAR_VC, VT );
        MakeZeros( VB );
    }

    // Backtransform U and V
    if( m >= n )
    {
        ApplyPackedReflectors
        ( LEFT, LOWER, VERTICAL, BACKWARD, 0, B, A );
        ApplyPackedReflectors
        ( LEFT, UPPER, HORIZONTAL, BACKWARD, 1, B, V );
    }
    else
    {
        ApplyPackedReflectors
        ( LEFT, LOWER, VERTICAL, BACKWARD, -1, B, A );
        ApplyPackedReflectors
        ( LEFT, UPPER, HORIZONTAL, BACKWARD, 0, B, V );
    }

    // Copy out the appropriate subset of the singular values
    s = d_STAR_STAR;
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifdef HAVE_FLA_BSVD
template<>
inline void
SimpleSVD
( DistMatrix<double>& A,
  DistMatrix<double,VR,STAR>& s,
  DistMatrix<double>& V )
{
#ifndef RELEASE
    PushCallStack("svd::SimpleSVD");
#endif
    typedef double R;
    typedef Complex<R> C;

    const int m = A.Height();
    const int n = A.Width();
    const int k = std::min( m, n );
    const int offdiagonal = ( m>=n ? 1 : -1 );
    const char uplo = ( m>=n ? 'U' : 'L' );
    const Grid& grid = A.Grid();

    // Bidiagonalize A
    Bidiag( A );

    // Grab copies of the diagonal and sub/super-diagonal of A
    DistMatrix<R,MD,STAR> d_MD_STAR( grid ), 
                          e_MD_STAR( grid );
    A.GetDiagonal( d_MD_STAR );
    A.GetDiagonal( e_MD_STAR, offdiagonal );

    // In order to use serial QR kernels, we need the full bidiagonal matrix
    // on each process
    DistMatrix<R,STAR,STAR> d_STAR_STAR( d_MD_STAR ),
                            e_STAR_STAR( e_MD_STAR );

    // Initialize U and V to the appropriate identity matrices.
    DistMatrix<R,VC,STAR> U_VC_STAR( grid );
    DistMatrix<R,VC,STAR> V_VC_STAR( grid );
    U_VC_STAR.AlignWith( A );
    V_VC_STAR.AlignWith( V );
    Identity( m, k, U_VC_STAR );
    Identity( n, k, V_VC_STAR );

    // Compute the SVD of the bidiagonal matrix and accumulate the Givens
    // rotations into our local portion of U and V
    // NOTE: This _only_ works in the case where m >= n
    const int numAccum = 32;
    const int maxNumIts = 30;
    const int bAlg = 512;
    std::vector<C> GBuffer( (k-1)*numAccum ),
                   HBuffer( (k-1)*numAccum );
    FLA_Bsvd_v_opd_var1
    ( k, U_VC_STAR.LocalHeight(), V_VC_STAR.LocalHeight(), 
      numAccum, maxNumIts,
      d_STAR_STAR.LocalBuffer(), 1,
      e_STAR_STAR.LocalBuffer(), 1,
      &GBuffer[0], 1, k-1,
      &HBuffer[0], 1, k-1,
      U_VC_STAR.LocalBuffer(), 1, U_VC_STAR.LocalLDim(),
      V_VC_STAR.LocalBuffer(), 1, V_VC_STAR.LocalLDim(),
      bAlg );

    // Make a copy of A (for the Householder vectors) and pull the necessary 
    // portions of U and V into a standard matrix dist.
    DistMatrix<R> B( A );
    if( m >= n )
    {
        DistMatrix<R> AT( grid ),
                      AB( grid );
        DistMatrix<R,VC,STAR> UT_VC_STAR( grid ), 
                              UB_VC_STAR( grid );
        PartitionDown( A, AT,
                          AB, n );
        PartitionDown( U_VC_STAR, UT_VC_STAR,
                                  UB_VC_STAR, n );
        AT = UT_VC_STAR;
        MakeZeros( AB );
        V = V_VC_STAR;
    }
    else
    {
        DistMatrix<R> VT( grid ), 
                      VB( grid );
        DistMatrix<R,VC,STAR> VT_VC_STAR( grid ), 
                              VB_VC_STAR( grid );
        PartitionDown( V, VT, 
                          VB, m );
        PartitionDown
        ( V_VC_STAR, VT_VC_STAR, 
                     VB_VC_STAR, m );
        VT = VT_VC_STAR;
        MakeZeros( VB );
    }

    // Backtransform U and V
    if( m >= n )
    {
        ApplyPackedReflectors
        ( LEFT, LOWER, VERTICAL, BACKWARD, 0, B, A );
        ApplyPackedReflectors
        ( LEFT, UPPER, HORIZONTAL, BACKWARD, 1, B, V );
    }
    else
    {
        ApplyPackedReflectors
        ( LEFT, LOWER, VERTICAL, BACKWARD, -1, B, A );
        ApplyPackedReflectors
        ( LEFT, UPPER, HORIZONTAL, BACKWARD, 0, B, V );
    }

    // Copy out the appropriate subset of the singular values
    s = d_STAR_STAR;
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // HAVE_FLA_BSVD

template<typename R>
inline void
SimpleSingularValues
( DistMatrix<R>& A,
  DistMatrix<R,VR,STAR>& s )
{
#ifndef RELEASE
    PushCallStack("svd::SimpleSingularValues");
#endif
    const int m = A.Height();
    const int n = A.Width();
    const int k = std::min( m, n );
    const int offdiagonal = ( m>=n ? 1 : -1 );
    const char uplo = ( m>=n ? 'U' : 'L' );
    const Grid& grid = A.Grid();

    // Bidiagonalize A
    Bidiag( A );

    // Grab copies of the diagonal and sub/super-diagonal of A
    DistMatrix<R,MD,STAR> d_MD_STAR( grid ), 
                          e_MD_STAR( grid );
    A.GetDiagonal( d_MD_STAR );
    A.GetDiagonal( e_MD_STAR, offdiagonal );

    // In order to use serial QR kernels, we need the full bidiagonal matrix
    // on each process
    DistMatrix<R,STAR,STAR> d_STAR_STAR( d_MD_STAR ),
                            e_STAR_STAR( e_MD_STAR );

    // Compute the singular values of the bidiagonal matrix
    lapack::BidiagQRAlg
    ( uplo, k, 0, 0,
      d_STAR_STAR.LocalBuffer(), e_STAR_STAR.LocalBuffer(), 
      0, 1, 0, 1 );

    // Copy out the appropriate subset of the singular values
    s = d_STAR_STAR;
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline void
SimpleSVD
( DistMatrix<Complex<R> >& A,
  DistMatrix<R,VR,STAR>& s,
  DistMatrix<Complex<R> >& V )
{
#ifndef RELEASE
    PushCallStack("svd::SimpleSVD");
#endif
    typedef Complex<R> C;
    const int m = A.Height();
    const int n = A.Width();
    const int k = std::min( m, n );
    const int offdiagonal = ( m>=n ? 1 : -1 );
    const char uplo = ( m>=n ? 'U' : 'L' );
    const Grid& grid = A.Grid();

    // Bidiagonalize A
    DistMatrix<C,STAR,STAR> tP( grid ), tQ( grid );
    Bidiag( A, tP, tQ );

    // Grab copies of the diagonal and sub/super-diagonal of A
    DistMatrix<R,MD,STAR> d_MD_STAR( grid ),
                          e_MD_STAR( grid );
    A.GetRealDiagonal( d_MD_STAR );
    A.GetRealDiagonal( e_MD_STAR, offdiagonal );

    // on each process
    DistMatrix<R,STAR,STAR> d_STAR_STAR( d_MD_STAR ),
                            e_STAR_STAR( e_MD_STAR );

    // Initialize U and VAdj to the appropriate identity matrices
    DistMatrix<C,VC,STAR> U_VC_STAR( grid );
    DistMatrix<C,STAR,VC> VAdj_STAR_VC( grid );
    U_VC_STAR.AlignWith( A );
    VAdj_STAR_VC.AlignWith( V );
    Identity( m, k, U_VC_STAR );
    Identity( k, n, VAdj_STAR_VC );

    // Compute the SVD of the bidiagonal matrix and accumulate the Givens
    // rotations into our local portion of U and VAdj
    Matrix<C>& ULocal = U_VC_STAR.LocalMatrix();
    Matrix<C>& VAdjLocal = VAdj_STAR_VC.LocalMatrix();
    lapack::BidiagQRAlg
    ( uplo, k, VAdjLocal.Width(), ULocal.Height(),
      d_STAR_STAR.LocalBuffer(), e_STAR_STAR.LocalBuffer(), 
      VAdjLocal.Buffer(), VAdjLocal.LDim(), 
      ULocal.Buffer(), ULocal.LDim() );

    // Make a copy of A (for the Householder vectors) and pull the necessary 
    // portions of U and VAdj into a standard matrix dist.
    DistMatrix<C> B( A );
    if( m >= n )
    {
        DistMatrix<C> AT( grid ),
                      AB( grid );
        DistMatrix<C,VC,STAR> UT_VC_STAR( grid ),
                              UB_VC_STAR( grid );
        PartitionDown( A, AT,
                          AB, n );
        PartitionDown( U_VC_STAR, UT_VC_STAR,
                                  UB_VC_STAR, n );
        AT = UT_VC_STAR;
        MakeZeros( AB );
        Adjoint( VAdj_STAR_VC, V );
    }
    else
    {
        DistMatrix<C> VT( grid ), 
                      VB( grid );
        DistMatrix<C,STAR,VC> VAdjL_STAR_VC( grid ), VAdjR_STAR_VC( grid );
        PartitionDown( V, VT, 
                          VB, m );
        PartitionRight( VAdj_STAR_VC, VAdjL_STAR_VC, VAdjR_STAR_VC, m );
        Adjoint( VAdjL_STAR_VC, VT );
        MakeZeros( VB );
    }

    // Backtransform U and V
    if( m >= n )
    {
        ApplyPackedReflectors
        ( LEFT, LOWER, VERTICAL, BACKWARD, UNCONJUGATED, 0, B, tQ, A );
        ApplyPackedReflectors
        ( LEFT, UPPER, HORIZONTAL, BACKWARD, UNCONJUGATED, 1, B, tP, V );
    }
    else
    {
        ApplyPackedReflectors
        ( LEFT, LOWER, VERTICAL, BACKWARD, UNCONJUGATED, -1, B, tQ, A );
        ApplyPackedReflectors
        ( LEFT, UPPER, HORIZONTAL, BACKWARD, UNCONJUGATED, 0, B, tP, V );
    }

    // Copy out the appropriate subset of the singular values
    s = d_STAR_STAR;
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifdef HAVE_FLA_BSVD
template<>
inline void
SimpleSVD
( DistMatrix<Complex<double> >& A,
  DistMatrix<double,VR,STAR>& s,
  DistMatrix<Complex<double> >& V )
{
#ifndef RELEASE
    PushCallStack("svd::SimpleSVD");
#endif
    typedef double R;
    typedef Complex<R> C;

    const int m = A.Height();
    const int n = A.Width();
    const int k = std::min( m, n );
    const int offdiagonal = ( m>=n ? 1 : -1 );
    const char uplo = ( m>=n ? 'U' : 'L' );
    const Grid& grid = A.Grid();

    // Bidiagonalize A
    DistMatrix<C,STAR,STAR> tP( grid ), tQ( grid );
    Bidiag( A, tP, tQ );

    // Grab copies of the diagonal and sub/super-diagonal of A
    DistMatrix<R,MD,STAR> d_MD_STAR( grid ),
                          e_MD_STAR( grid );
    A.GetRealDiagonal( d_MD_STAR );
    A.GetRealDiagonal( e_MD_STAR, offdiagonal );

    // on each process
    DistMatrix<R,STAR,STAR> d_STAR_STAR( d_MD_STAR ),
                            e_STAR_STAR( e_MD_STAR );

    // Initialize U and VAdj to the appropriate identity matrices
    DistMatrix<C,VC,STAR> U_VC_STAR( grid );
    DistMatrix<C,VC,STAR> V_VC_STAR( grid );
    U_VC_STAR.AlignWith( A );
    V_VC_STAR.AlignWith( V );
    Identity( m, k, U_VC_STAR );
    Identity( n, k, V_VC_STAR );

    // Compute the SVD of the bidiagonal matrix and accumulate the Givens
    // rotations into our local portion of U and V
    // NOTE: This _only_ works in the case where m >= n
    const int numAccum = 32;
    const int maxNumIts = 30;
    const int bAlg = 512;
    std::vector<C> GBuffer( (k-1)*numAccum ),
                   HBuffer( (k-1)*numAccum );
    FLA_Bsvd_v_opz_var1
    ( k, U_VC_STAR.LocalHeight(), V_VC_STAR.LocalHeight(), 
      numAccum, maxNumIts,
      d_STAR_STAR.LocalBuffer(), 1,
      e_STAR_STAR.LocalBuffer(), 1,
      &GBuffer[0], 1, k-1,
      &HBuffer[0], 1, k-1,
      U_VC_STAR.LocalBuffer(), 1, U_VC_STAR.LocalLDim(),
      V_VC_STAR.LocalBuffer(), 1, V_VC_STAR.LocalLDim(),
      bAlg );

    // Make a copy of A (for the Householder vectors) and pull the necessary 
    // portions of U and V into a standard matrix dist.
    DistMatrix<C> B( A );
    if( m >= n )
    {
        DistMatrix<C> AT( grid ),
                      AB( grid );
        DistMatrix<C,VC,STAR> UT_VC_STAR( grid ), 
                              UB_VC_STAR( grid );
        PartitionDown( A, AT,
                          AB, n );
        PartitionDown( U_VC_STAR, UT_VC_STAR,
                                  UB_VC_STAR, n );
        AT = UT_VC_STAR;
        MakeZeros( AB );
        V = V_VC_STAR;
    }
    else
    {
        DistMatrix<C> VT( grid ), 
                      VB( grid );
        DistMatrix<C,VC,STAR> VT_VC_STAR( grid ), 
                              VB_VC_STAR( grid );
        PartitionDown( V, VT, 
                          VB, m );
        PartitionDown
        ( V_VC_STAR, VT_VC_STAR, 
                     VB_VC_STAR, m );
        VT = VT_VC_STAR;
        MakeZeros( VB );
    }

    // Backtransform U and V
    if( m >= n )
    {
        ApplyPackedReflectors
        ( LEFT, LOWER, VERTICAL, BACKWARD, UNCONJUGATED, 0, B, tQ, A );
        ApplyPackedReflectors
        ( LEFT, UPPER, HORIZONTAL, BACKWARD, UNCONJUGATED, 1, B, tP, V );
    }
    else
    {
        ApplyPackedReflectors
        ( LEFT, LOWER, VERTICAL, BACKWARD, UNCONJUGATED, -1, B, tQ, A );
        ApplyPackedReflectors
        ( LEFT, UPPER, HORIZONTAL, BACKWARD, UNCONJUGATED, 0, B, tP, V );
    }

    // Copy out the appropriate subset of the singular values
    s = d_STAR_STAR;
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // HAVE_FLA_BSVD

template<typename F>
inline void
CheckScale
( DistMatrix<F>& A, bool& needRescaling, typename Base<F>::type& scale )
{
    typedef typename Base<F>::type R;

    scale = 1;
    needRescaling = false;
    const R oneNormOfA = Norm( A, ONE_NORM );
    const R safeMin = lapack::MachineSafeMin<R>();
    const R precision = lapack::MachinePrecision<R>();
    const R smallNumber = safeMin/precision;
    const R bigNumber = 1/smallNumber;
    const R rhoMin = Sqrt(smallNumber);
    const R rhoMax = std::min( Sqrt(bigNumber), 1/Sqrt(Sqrt(safeMin)) );

    if( oneNormOfA > 0 && oneNormOfA < rhoMin )
    {
        needRescaling = true;
        scale = rhoMin/oneNormOfA;
    }
    else if( oneNormOfA > rhoMax )
    {
        needRescaling = true;
        scale = rhoMax/oneNormOfA;
    }
}

template<typename F>
inline void
DivideAndConquerSVD
( Matrix<F>& A, Matrix<typename Base<F>::type>& s, Matrix<F>& V )
{
#ifndef RELEASE
    PushCallStack("svd::DivideAndConquerSVD");
#endif
    typedef typename Base<F>::type R;

    const int m = A.Height();
    const int n = A.Width();
    const int k = std::min(m,n);
    s.ResizeTo( k, 1 );
    Matrix<F> U( m, k );
    Matrix<F> VAdj( k, n );
    lapack::DivideAndConquerSVD
    ( m, n, A.Buffer(), A.LDim(), s.Buffer(), U.Buffer(), U.LDim(),
      VAdj.Buffer(), VAdj.LDim() );

    A = U;
    Adjoint( VAdj, V );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
QRSVD( Matrix<F>& A, Matrix<typename Base<F>::type>& s, Matrix<F>& V )
{
#ifndef RELEASE
    PushCallStack("svd::QRSVD");
#endif
    typedef typename Base<F>::type R;

    const int m = A.Height();
    const int n = A.Width();
    const int k = std::min(m,n);
    s.ResizeTo( k, 1 );
    Matrix<F> U( m, k );
    Matrix<F> VAdj( k, n );
    lapack::QRSVD
    ( m, n, A.Buffer(), A.LDim(), s.Buffer(), U.Buffer(), U.LDim(),
      VAdj.Buffer(), VAdj.LDim() );

    A = U;
    Adjoint( VAdj, V );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace svd

//----------------------------------------------------------------------------//
// Grab the full SVD of the general matrix A, A = U diag(s) V^H.              //
// On exit, A is overwritten with U.                                          //
//----------------------------------------------------------------------------//

template<typename F>
inline void
SVD( Matrix<F>& A, Matrix<typename Base<F>::type>& s, Matrix<F>& V, bool useQR )
{
#ifndef RELEASE
    PushCallStack("SVD");
#endif
    if( useQR )
        svd::QRSVD( A, s, V );
    else
        svd::DivideAndConquerSVD( A, s, V );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
SVD
( DistMatrix<F>& A,
  DistMatrix<typename Base<F>::type,VR,STAR>& s,
  DistMatrix<F>& V )
{
#ifndef RELEASE
    PushCallStack("SVD");
#endif
    typedef typename Base<F>::type R;

    // Check if we need to rescale the matrix, and do so if necessary
    bool needRescaling;
    R scale;
    svd::CheckScale( A, needRescaling, scale );
    if( needRescaling )
        Scale( scale, A );

    // TODO: Switch between different algorithms. For instance, starting 
    //       with a QR decomposition of tall-skinny matrices.
    if( A.Height() >= A.Width() )
    {
        svd::SimpleSVD( A, s, V );
    }
    else
    {
        // Lower bidiagonalization is not yet supported, so we instead play a 
        // trick to get the SVD of A.
        Adjoint( A, V );
        svd::SimpleSVD( V, s, A );
    }

    // Rescale the singular values if necessary
    if( needRescaling )
        Scal( 1/scale, s );
#ifndef RELEASE
    PopCallStack();
#endif
}

//----------------------------------------------------------------------------//
// Grab the singular values of the general matrix A.                          //
//----------------------------------------------------------------------------//
template<typename F>
inline void
SingularValues( Matrix<F>& A, Matrix<typename Base<F>::type>& s )
{
#ifndef RELEASE
    PushCallStack("SingularValues");
#endif
    typedef typename Base<F>::type R;

    const int m = A.Height();
    const int n = A.Width();
    const int k = std::min(m,n);
    s.ResizeTo( k, 1 );
    lapack::SingularValues( m, n, A.Buffer(), A.LDim(), s.Buffer() );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F>
inline void
SingularValues
( DistMatrix<F>& A,
  DistMatrix<typename Base<F>::type,VR,STAR>& s )
{
#ifndef RELEASE
    PushCallStack("SingularValues");
#endif
    typedef typename Base<F>::type R;

    // Check if we need to rescale the matrix, and do so if necessary
    bool needRescaling;
    R scale;
    svd::CheckScale( A, needRescaling, scale );
    if( needRescaling )
        Scale( scale, A );

    // TODO: Switch between different algorithms. For instance, starting 
    //       with a QR decomposition of tall-skinny matrices.
    if( A.Height() >= A.Width() )
    {
        svd::SimpleSingularValues( A, s );
    }
    else
    {
        // Lower bidiagonalization is not yet supported, so we instead play a 
        // trick to get the SVD of A.
        DistMatrix<F> AAdj( A.Grid() );
        Adjoint( A, AAdj );
        svd::SimpleSingularValues( AAdj, s );
    }

    // Rescale the singular values if necessary
    if( needRescaling )
        Scal( 1/scale, s );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
