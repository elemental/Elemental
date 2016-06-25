/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SCHUR_QR_HPP
#define EL_SCHUR_QR_HPP

namespace El {
namespace schur {

template<typename F>
void QR
( Matrix<F>& A,
  Matrix<Complex<Base<F>>>& w,
  bool fullTriangle,
  bool time=false )
{
    DEBUG_CSE
    const Int n = A.Height();
    Timer timer;

    Matrix<F> phase;
    if( time )
        timer.Start();
    Hessenberg( UPPER, A, phase );
    if( time )
        Output("  Hessenberg reduction: ",timer.Stop()," seconds");
    MakeTrapezoidal( UPPER, A, -1 );

    HessenbergSchurCtrl hessSchurCtrl;
    hessSchurCtrl.fullTriangle = fullTriangle;
    hessSchurCtrl.wantSchurVecs = false;
    hessSchurCtrl.demandConverged = true;
    hessSchurCtrl.useAED = true;
    hessSchurCtrl.recursiveAED = true;
    if( time )
        timer.Start();
    HessenbergSchur( A, w, hessSchurCtrl ); 
    if( time )
        Output("  HessenbergSchur: ",timer.Stop()," seconds");

    if( IsComplex<F>::value )
        MakeTrapezoidal( UPPER, A );
    else
    {
        MakeTrapezoidal( UPPER, A, -1 );
        DEBUG_ONLY(CheckRealSchur(A))
    }
}

template<typename F>
void QR
( Matrix<F>& A,
  Matrix<Complex<Base<F>>>& w,
  Matrix<F>& Q, 
  bool fullTriangle,
  bool time=false )
{
    DEBUG_CSE
    const Int n = A.Height();
    Timer timer;

    Matrix<F> phase;
    if( time )
        timer.Start();
    Hessenberg( UPPER, A, phase );
    if( time )
        Output("  Hessenberg reduction: ",timer.Stop()," seconds");
    hessenberg::FormQ( UPPER, A, phase, Q );
    MakeTrapezoidal( UPPER, A, -1 );

    HessenbergSchurCtrl hessSchurCtrl;
    hessSchurCtrl.fullTriangle = fullTriangle;
    hessSchurCtrl.wantSchurVecs = true;
    hessSchurCtrl.demandConverged = true;
    hessSchurCtrl.useAED = true;
    hessSchurCtrl.recursiveAED = true;
    if( time )
        timer.Start();
    HessenbergSchur( A, w, Q, hessSchurCtrl ); 
    if( time )
        Output("  HessenbergSchur: ",timer.Stop()," seconds");

    if( IsComplex<F>::value )
        MakeTrapezoidal( UPPER, A );
    else
    {
        MakeTrapezoidal( UPPER, A, -1 );
        DEBUG_ONLY(CheckRealSchur(A))
    }
}

template<typename F>
void QR
( DistMatrix<F,MC,MR,BLOCK>& A,
  ElementalMatrix<Complex<Base<F>>>& w,
  bool fullTriangle,
  const HessQRCtrl& ctrl,
  bool time=false )
{
    DEBUG_CSE
    AssertScaLAPACKSupport();
#ifdef EL_HAVE_SCALAPACK
    const Int n = A.Height();
    const int bHandle = blacs::Handle( A );
    const int context = blacs::GridInit( bHandle, A );
    auto descA = FillDesc( A, context );

    Timer timer;
    const int gridRank = A.Grid().Rank();

    // Reduce the matrix to upper-Hessenberg form in an elemental form
    DistMatrix<F> AElem( A );
    DistMatrix<F,STAR,STAR> t( A.Grid() );
    if( time && gridRank == 0 )
        timer.Start();
    Hessenberg( UPPER, AElem, t );
    if( time && gridRank == 0 )
        Output("  ScaLAPACK Hessenberg: ",timer.Stop()," seconds"); 
    MakeTrapezoidal( UPPER, AElem, -1 );
    if( time && gridRank == 0 )
        timer.Start();
    A = AElem;
    if( time && gridRank == 0 )
        Output
        ("  Redist of A from elemental to block: ",timer.Stop()," seconds");

    // Run the QR algorithm in block form
    DistMatrix<Complex<Base<F>>,STAR,STAR> w_STAR_STAR( n, 1, A.Grid() );

    if( time && gridRank == 0 )
        timer.Start();
#define FORCE_WANTZ_TRUE 1
#if FORCE_WANTZ_TRUE
    DistMatrix<F,MC,MR,BLOCK> Z(n,n,A.Grid(),A.BlockHeight(),A.BlockWidth());
    Identity( Z, n, n );
    bool multiplyZ=true;
    scalapack::HessenbergSchur
    ( n,
      A.Buffer(), descA.data(),
      w_STAR_STAR.Buffer(),
      Z.Buffer(), descA.data(),
      fullTriangle, multiplyZ, ctrl.distAED );
#else
    scalapack::HessenbergSchur
    ( n,
      A.Buffer(), descA.data(),
      w_STAR_STAR.Buffer(),
      fullTriangle, ctrl.distAED );
#endif
    if( time && gridRank == 0 )
        Output("  ScaLAPACK HessenbergSchur: ",timer.Stop()," seconds");
    Copy( w_STAR_STAR, w );

    // TODO: Cache context, handle, and exit BLACS during El::Finalize()
    blacs::FreeGrid( context );
    blacs::FreeHandle( bHandle );
#endif
    if( IsComplex<F>::value )
        MakeTrapezoidal( UPPER, A );
    else
    {
        MakeTrapezoidal( UPPER, A, -1 );
        // NOTE: This routine is not yet implemented
        //DEBUG_ONLY(CheckRealSchur(A))
    }
}

template<typename F>
void QR
( DistMatrix<F,MC,MR,BLOCK>& A,
  ElementalMatrix<Complex<Base<F>>>& w,
  DistMatrix<F,MC,MR,BLOCK>& Q,
  bool fullTriangle,
  const HessQRCtrl& ctrl,
  bool time=false )
{
    DEBUG_CSE
    AssertScaLAPACKSupport();
#ifdef EL_HAVE_SCALAPACK
    const Int n = A.Height();
    const int bHandle = blacs::Handle( A );
    const int context = blacs::GridInit( bHandle, A );
    Q.AlignWith( A );
    Q.Resize( n, n, A.LDim() );
    auto descA = FillDesc( A, context );
    auto descQ = FillDesc( Q, context );

    Timer timer;
    const int gridRank = A.Grid().Rank();

    // Reduce A to upper-Hessenberg form in an element-wise distribution
    // and form the explicit reflector matrix
    DistMatrix<F> AElem( A ), QElem( A.Grid() );
    DistMatrix<F,STAR,STAR> phase( A.Grid() );
    if( time && gridRank == 0 )
        timer.Start();
    Hessenberg( UPPER, AElem, phase );
    if( time && gridRank == 0 )
        Output("  Hessenberg: ",timer.Stop()," seconds");
    // There is not yet a 'form Q'
    if( time && gridRank == 0 )
        timer.Start();
    Identity( QElem, n, n ); 
    hessenberg::FormQ( UPPER, AElem, phase, QElem );
    if( time && gridRank == 0 )
        Output("  hessenberg::FormQ: ",timer.Stop()," seconds");
    MakeTrapezoidal( UPPER, AElem, -1 );
    if( time && gridRank == 0 )
        timer.Start();
    A = AElem;
    Q = QElem;
    if( time && gridRank == 0 )
        Output
        ("  Redist of (A,Q) from elemental to block: ",timer.Stop()," seconds");
    
    // Compute the Schur decomposition in block form, multiplying the 
    // accumulated Householder reflectors from the right
    DistMatrix<Complex<Base<F>>,STAR,STAR> w_STAR_STAR( n, 1, A.Grid() );
    const bool multiplyQ = true;
    if( time && gridRank == 0 )
        timer.Start();
    scalapack::HessenbergSchur
    ( n,
      A.Buffer(), descA.data(),
      w_STAR_STAR.Buffer(), 
      Q.Buffer(), descQ.data(),
      fullTriangle, multiplyQ, ctrl.distAED );
    if( time && gridRank == 0 )
        Output("  ScaLAPACK HessenbergSchur: ",timer.Stop()," seconds");
    Copy( w_STAR_STAR, w );

    // TODO: Cache context, handle, and exit BLACS during El::Finalize()
    blacs::FreeGrid( context );
    blacs::FreeHandle( bHandle );
#endif
    if( IsComplex<F>::value )
        MakeTrapezoidal( UPPER, A );
    else
    {
        MakeTrapezoidal( UPPER, A, -1 );
        // NOTE: This routine is not yet implemented
        //DEBUG_ONLY(CheckRealSchur(A))
    }
}

template<typename F>
void QR
( ElementalMatrix<F>& APre,
  ElementalMatrix<Complex<Base<F>>>& w, 
  bool fullTriangle,
  const HessQRCtrl& ctrl,
  bool time=false )
{
    DEBUG_CSE
    AssertScaLAPACKSupport();

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.Get();

#ifdef EL_HAVE_SCALAPACK
    Timer timer;
    const int gridRank = A.Grid().Rank();

    // Reduce the matrix to upper-Hessenberg form in an elemental form
    DistMatrix<F,STAR,STAR> t( A.Grid() );
    if( time && gridRank == 0 )
        timer.Start();
    Hessenberg( UPPER, A, t );
    if( time && gridRank == 0 )
        Output("  Hessenberg: ",timer.Stop()," seconds");
    MakeTrapezoidal( UPPER, A, -1 );

    // Run the QR algorithm in block form
    // TODO: Create schur::HessenbergQR
    const Int n = A.Height(); 
    const Int mb = ctrl.blockHeight;
    const Int nb = ctrl.blockWidth;
    DistMatrix<F,MC,MR,BLOCK> ABlock( n, n, A.Grid(), mb, nb );
    if( time && gridRank == 0 )
        timer.Start();
    ABlock = A;
    if( time && gridRank == 0 )
        Output("  Redist. from elemental to block: ",timer.Stop()," seconds");
    const int bHandle = blacs::Handle( ABlock );
    const int context = blacs::GridInit( bHandle, ABlock );
    blacs::Desc descA = FillDesc( ABlock, context );
    DistMatrix<Complex<Base<F>>,STAR,STAR> w_STAR_STAR( n, 1, A.Grid() );

    if( time && gridRank == 0 )    
        timer.Start();
#define FORCE_WANTZ_TRUE 1
#if FORCE_WANTZ_TRUE
    DistMatrix<F,MC,MR,BLOCK> Z(n,n,A.Grid(),mb,nb);
    Identity( Z, n, n );
    blacs::Desc descZ = FillDesc( Z, context );
    bool multiplyZ=true;
    scalapack::HessenbergSchur
    ( n,
      ABlock.Buffer(), descA.data(),
      w_STAR_STAR.Buffer(),
      Z.Buffer(), descZ.data(),
      fullTriangle, multiplyZ, ctrl.distAED );
#else
    scalapack::HessenbergSchur
    ( n, 
      ABlock.Buffer(), descA.data(),
      w_STAR_STAR.Buffer(), 
      fullTriangle, ctrl.distAED );
#endif
    if( time && gridRank == 0 )    
        Output("  scalapack::HessenbergSchur: ",timer.Stop()," seconds");

    if( time && gridRank == 0 )
        timer.Start();
    A = ABlock;
    if( time && gridRank == 0 )
        Output("  Redist. from block to elemental: ",timer.Stop()," seconds");
    Copy( w_STAR_STAR, w );

    // TODO: Cache context, handle, and exit BLACS during El::Finalize()
    blacs::FreeGrid( context );
    blacs::FreeHandle( bHandle );
#endif
    if( IsComplex<F>::value )
        MakeTrapezoidal( UPPER, A );
    else
    {
        MakeTrapezoidal( UPPER, A, -1 );
        DEBUG_ONLY(CheckRealSchur(A))
    }
}

template<typename F>
void QR
( ElementalMatrix<F>& APre,
  ElementalMatrix<Complex<Base<F>>>& w, 
  ElementalMatrix<F>& QPre,
  bool fullTriangle,
  const HessQRCtrl& ctrl,
  bool time=false )
{
    DEBUG_CSE
    AssertScaLAPACKSupport();

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    DistMatrixWriteProxy<F,F,MC,MR> QProx( QPre );
    auto& A = AProx.Get();
    auto& Q = QProx.Get();

#ifdef EL_HAVE_SCALAPACK
    Timer timer;
    const int gridRank = A.Grid().Rank();

    const Int n = A.Height();
    // Reduce A to upper-Hessenberg form in an element-wise distribution
    // and form the explicit reflector matrix
    DistMatrix<F,STAR,STAR> phase( A.Grid() );
    if( time && gridRank == 0 )
        timer.Start();
    Hessenberg( UPPER, A, phase );
    if( time && gridRank == 0 )
        Output("  Hessenberg: ",timer.Stop()," seconds");
    // There is not yet a 'form Q'
    Identity( Q, n, n ); 
    if( time && gridRank == 0 )
        timer.Start();
    hessenberg::FormQ( UPPER, A, phase, Q );
    if( time && gridRank == 0 )
        Output("  hessenberg::FormQ: ",timer.Stop()," seconds");
    MakeTrapezoidal( UPPER, A, -1 );

    // Run the Hessenberg QR algorithm in block form
    const Int mb = ctrl.blockHeight;
    const Int nb = ctrl.blockWidth;
    DistMatrix<F,MC,MR,BLOCK>
      ABlock( n, n, A.Grid(), mb, nb ), 
      QBlock( n, n, A.Grid(), mb, nb );
    if( time && gridRank == 0 )
        timer.Start();
    ABlock = A;
    QBlock = Q;
    if( time && gridRank == 0 )
        Output
        ("  Redist of (A,Q) from elemental to block: ",timer.Stop()," seconds");
    const int bHandle = blacs::Handle( ABlock );
    const int context = blacs::GridInit( bHandle, ABlock );
    auto descA = FillDesc( ABlock, context );
    auto descQ = FillDesc( QBlock, context );

    // Compute the Schur decomposition in block form, multiplying the 
    // accumulated Householder reflectors from the right
    DistMatrix<Complex<Base<F>>,STAR,STAR> w_STAR_STAR( n, 1, A.Grid() );
    const bool multiplyQ = true;
    if( time && gridRank == 0 )
        timer.Start();
    scalapack::HessenbergSchur
    ( n,
      ABlock.Buffer(), descA.data(),
      w_STAR_STAR.Buffer(), 
      QBlock.Buffer(), descQ.data(),
      fullTriangle, multiplyQ, ctrl.distAED );
    if( time && gridRank == 0 )
        Output("  scalapack::HessenbergSchur: ",timer.Stop()," seconds");
    if( time && gridRank == 0 )
        timer.Start();
    A = ABlock;
    Q = QBlock;
    if( time && gridRank == 0 )
        Output
        ("  Redist of (A,Q) from block to elemental: ",timer.Stop()," seconds");
    Copy( w_STAR_STAR, w );

    // TODO: Cache context, handle, and exit BLACS during El::Finalize()
    blacs::FreeGrid( context );
    blacs::FreeHandle( bHandle );
#endif
    if( IsComplex<F>::value )
        MakeTrapezoidal( UPPER, A );
    else
    {
        MakeTrapezoidal( UPPER, A, -1 );
        DEBUG_ONLY(CheckRealSchur(A))
    }
}

} // namespace schur
} // namespace El

#endif // ifndef EL_SCHUR_QR_HPP
