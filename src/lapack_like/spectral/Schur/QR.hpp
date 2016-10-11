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
  const SchurCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    Timer timer;

    Matrix<F> householderScalars;
    if( ctrl.time )
        timer.Start();
    Hessenberg( UPPER, A, householderScalars );
    if( ctrl.time )
        Output("  Hessenberg reduction: ",timer.Stop()," seconds");
    MakeTrapezoidal( UPPER, A, -1 );

    if( ctrl.time )
        timer.Start();
    HessenbergSchur( A, w, ctrl.hessSchurCtrl ); 
    if( ctrl.time )
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
  const SchurCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    Timer timer;

    Matrix<F> householderScalars;
    if( ctrl.time )
        timer.Start();
    Hessenberg( UPPER, A, householderScalars );
    if( ctrl.time )
        Output("  Hessenberg reduction: ",timer.Stop()," seconds");
    hessenberg::FormQ( UPPER, A, householderScalars, Q );
    MakeTrapezoidal( UPPER, A, -1 );

    auto hessSchurCtrl( ctrl.hessSchurCtrl );
    hessSchurCtrl.accumulateSchurVecs = true;
    if( ctrl.time )
        timer.Start();
    HessenbergSchur( A, w, Q, hessSchurCtrl ); 
    if( ctrl.time )
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
  const SchurCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    const Int n = A.Height();
    const Grid& grid = A.Grid();
    Timer timer;

    // Reduce the matrix to upper-Hessenberg form in an elemental form
    DistMatrix<F> AElem( A );
    DistMatrix<F,STAR,STAR> t(grid);
    if( ctrl.time && grid.Rank() == 0 )
        timer.Start();
    Hessenberg( UPPER, AElem, t );
    if( ctrl.time && grid.Rank() == 0 )
        Output("  ScaLAPACK Hessenberg: ",timer.Stop()," seconds"); 
    MakeTrapezoidal( UPPER, AElem, -1 );
    if( ctrl.time && grid.Rank() == 0 )
        timer.Start();
    A = AElem;
    if( ctrl.time && grid.Rank() == 0 )
        Output
        ("  Redist of A from elemental to block: ",timer.Stop()," seconds");

    // Run the QR algorithm in block form
    DistMatrix<Complex<Base<F>>,STAR,STAR> w_STAR_STAR( n, 1, grid );

    if( ctrl.time && grid.Rank() == 0 )
        timer.Start();
    HessenbergSchur( A, w_STAR_STAR, ctrl.hessSchurCtrl );
    /*
#ifdef EL_HAVE_SCALAPACK
    const int bHandle = blacs::Handle( A );
    const int context = blacs::GridInit( bHandle, A );
    auto descA = FillDesc( A, context );
    // TODO(poulson): Move the following into HessenbergSchur
#define FORCE_WANTZ_TRUE 1
#if FORCE_WANTZ_TRUE
    DistMatrix<F,MC,MR,BLOCK> Z(n,n,grid,A.BlockHeight(),A.BlockWidth());
    Identity( Z, n, n );
    bool accumulateSchurVecs=true;
    scalapack::HessenbergSchur
    ( n,
      A.Buffer(), descA.data(),
      w_STAR_STAR.Buffer(),
      Z.Buffer(), descA.data(),
      ctrl.hessSchurCtrl.fullTriangle,
      accumulateSchurVecs,
      ctrl.hessSchurCtrl.scalapackAED );
#else
    scalapack::HessenbergSchur
    ( n,
      A.Buffer(), descA.data(),
      w_STAR_STAR.Buffer(),
      ctrl.hessSchurCtrl.fullTriangle,
      ctrl.hessSchurCtrl.scalapackAED );
#endif
    if( ctrl.time && gridRank == 0 )
        Output("  ScaLAPACK HessenbergSchur: ",timer.Stop()," seconds");

    // TODO: Cache context, handle, and exit BLACS during El::Finalize()
    blacs::FreeGrid( context );
    blacs::FreeHandle( bHandle );
#endif
    */
    Copy( w_STAR_STAR, w );

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
  const SchurCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    const Int n = A.Height();
    const Grid& grid = A.Grid();
    Timer timer;

    // Reduce A to upper-Hessenberg form in an element-wise distribution
    // and form the explicit reflector matrix
    DistMatrix<F> AElem( A ), QElem( A.Grid() );
    DistMatrix<F,STAR,STAR> householderScalars( A.Grid() );
    if( ctrl.time && grid.Rank() == 0 )
        timer.Start();
    Hessenberg( UPPER, AElem, householderScalars );
    if( ctrl.time && grid.Rank() == 0 )
        Output("  Hessenberg: ",timer.Stop()," seconds");
    // There is not yet a 'form Q'
    if( ctrl.time && grid.Rank() == 0 )
        timer.Start();
    Identity( QElem, n, n ); 
    hessenberg::FormQ( UPPER, AElem, householderScalars, QElem );
    if( ctrl.time && grid.Rank() == 0 )
        Output("  hessenberg::FormQ: ",timer.Stop()," seconds");
    MakeTrapezoidal( UPPER, AElem, -1 );
    if( ctrl.time && grid.Rank() == 0 )
        timer.Start();
    A = AElem;
    Q.AlignWith( A );
    Q.Resize( n, n, A.LDim() );
    Q = QElem;
    if( ctrl.time && grid.Rank() == 0 )
        Output
        ("  Redist of (A,Q) from elemental to block: ",timer.Stop()," seconds");
    
    // Compute the Schur decomposition in block form, multiplying the 
    // accumulated Householder reflectors from the right
    DistMatrix<Complex<Base<F>>,STAR,STAR> w_STAR_STAR( n, 1, A.Grid() );
    auto hessSchurCtrl = ctrl.hessSchurCtrl;
    hessSchurCtrl.accumulateSchurVecs = true;
    HessenbergSchur( A, w_STAR_STAR, Q, hessSchurCtrl );
    /*
#ifdef EL_HAVE_SCALAPACK
    const int bHandle = blacs::Handle( A );
    const int context = blacs::GridInit( bHandle, A );
    auto descA = FillDesc( A, context );
    auto descQ = FillDesc( Q, context );

    const bool accumulateSchurVecs = true;
    if( ctrl.time && grid.Rank() == 0 )
        timer.Start();
    scalapack::HessenbergSchur
    ( n,
      A.Buffer(), descA.data(),
      w_STAR_STAR.Buffer(), 
      Q.Buffer(), descQ.data(),
      ctrl.hessSchurCtrl.fullTriangle,
      accumulateSchurVecs,
      ctrl.hessSchurCtrl.scalapackAED );
    if( ctrl.time && grid.Rank() == 0 )
        Output("  ScaLAPACK HessenbergSchur: ",timer.Stop()," seconds");

    // TODO: Cache context, handle, and exit BLACS during El::Finalize()
    blacs::FreeGrid( context );
    blacs::FreeHandle( bHandle );
#endif
    */
    Copy( w_STAR_STAR, w );

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
  const SchurCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    Timer timer;

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.Get();

    const Int n = A.Height(); 
    const Grid& grid = A.Grid();

    // Reduce the matrix to upper-Hessenberg form in an elemental form
    DistMatrix<F,STAR,STAR> t( A.Grid() );
    if( ctrl.time && grid.Rank() == 0 )
        timer.Start();
    Hessenberg( UPPER, A, t );
    if( ctrl.time && grid.Rank() == 0 )
        Output("  Hessenberg: ",timer.Stop()," seconds");
    MakeTrapezoidal( UPPER, A, -1 );

    DistMatrix<Complex<Base<F>>,STAR,STAR> w_STAR_STAR( n, 1, A.Grid() );
    HessenbergSchur( A, w_STAR_STAR, ctrl.hessSchurCtrl );
    /*
#ifdef EL_HAVE_SCALAPACK
    // Run the QR algorithm in block form
    // TODO: Create schur::HessenbergQR
    const Int nb = ctrl.hessSchurCtrl.blockHeight;
    DistMatrix<F,MC,MR,BLOCK> ABlock( n, n, A.Grid(), nb, nb );
    if( ctrl.time && grid.Rank() == 0 )
        timer.Start();
    ABlock = A;
    if( ctrl.time && grid.Rank() == 0 )
        Output("  Redist. from elemental to block: ",timer.Stop()," seconds");
    const int bHandle = blacs::Handle( ABlock );
    const int context = blacs::GridInit( bHandle, ABlock );
    blacs::Desc descA = FillDesc( ABlock, context );
    if( ctrl.time && grid.Rank() == 0 )    
        timer.Start();
#define FORCE_WANTZ_TRUE 1
#if FORCE_WANTZ_TRUE
    DistMatrix<F,MC,MR,BLOCK> Z(n,n,A.Grid(),nb,nb);
    Identity( Z, n, n );
    blacs::Desc descZ = FillDesc( Z, context );
    bool accumulateSchurVecs=true;
    scalapack::HessenbergSchur
    ( n,
      ABlock.Buffer(), descA.data(),
      w_STAR_STAR.Buffer(),
      Z.Buffer(), descZ.data(),
      ctrl.hessSchurCtrl.fullTriangle,
      accumulateSchurVecs,
      ctrl.hessSchurCtrl.scalapackAED );
#else
    scalapack::HessenbergSchur
    ( n, 
      ABlock.Buffer(), descA.data(),
      w_STAR_STAR.Buffer(), 
      ctrl.hessSchurCtrl.fullTriangle,
      ctrl.hessSchurCtrl.scalapackAED );
#endif
    // TODO: Cache context, handle, and exit BLACS during El::Finalize()
    blacs::FreeGrid( context );
    blacs::FreeHandle( bHandle );

    if( ctrl.time && grid.Rank() == 0 )    
        Output("  scalapack::HessenbergSchur: ",timer.Stop()," seconds");

    if( ctrl.time && grid.Rank() == 0 )
        timer.Start();
    A = ABlock;
    if( ctrl.time && grid.Rank() == 0 )
        Output("  Redist. from block to elemental: ",timer.Stop()," seconds");
#endif
    */
    Copy( w_STAR_STAR, w );

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
  const SchurCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    DistMatrixWriteProxy<F,F,MC,MR> QProx( QPre );
    auto& A = AProx.Get();
    auto& Q = QProx.Get();

    const Int n = A.Height();
    const Grid& grid = A.Grid();
    Timer timer;

    // Reduce A to upper-Hessenberg form in an element-wise distribution
    // and form the explicit reflector matrix
    DistMatrix<F,STAR,STAR> householderScalars( A.Grid() );
    if( ctrl.time && grid.Rank() == 0 )
        timer.Start();
    Hessenberg( UPPER, A, householderScalars );
    if( ctrl.time && grid.Rank() == 0 )
        Output("  Hessenberg: ",timer.Stop()," seconds");
    // There is not yet a 'form Q'
    Identity( Q, n, n ); 
    if( ctrl.time && grid.Rank() == 0 )
        timer.Start();
    hessenberg::FormQ( UPPER, A, householderScalars, Q );
    if( ctrl.time && grid.Rank() == 0 )
        Output("  hessenberg::FormQ: ",timer.Stop()," seconds");
    MakeTrapezoidal( UPPER, A, -1 );

    DistMatrix<Complex<Base<F>>,STAR,STAR> w_STAR_STAR( n, 1, A.Grid() );
    auto hessSchurCtrl = ctrl.hessSchurCtrl;
    hessSchurCtrl.accumulateSchurVecs = true;
    HessenbergSchur( A, w_STAR_STAR, Q, hessSchurCtrl );
    /*
#ifdef EL_HAVE_SCALAPACK
    // Run the Hessenberg QR algorithm in block form
    const Int nb = ctrl.hessSchurCtrl.blockHeight;
    DistMatrix<F,MC,MR,BLOCK>
      ABlock( n, n, A.Grid(), nb, nb ), 
      QBlock( n, n, A.Grid(), nb, nb );
    if( ctrl.time && grid.Rank() == 0 )
        timer.Start();
    ABlock = A;
    QBlock = Q;
    if( ctrl.time && grid.Rank() == 0 )
        Output
        ("  Redist of (A,Q) from elemental to block: ",timer.Stop()," seconds");
    const int bHandle = blacs::Handle( ABlock );
    const int context = blacs::GridInit( bHandle, ABlock );
    auto descA = FillDesc( ABlock, context );
    auto descQ = FillDesc( QBlock, context );

    // Compute the Schur decomposition in block form, multiplying the 
    // accumulated Householder reflectors from the right
    const bool accumulateSchurVecs = true;
    if( ctrl.time && grid.Rank() == 0 )
        timer.Start();
    scalapack::HessenbergSchur
    ( n,
      ABlock.Buffer(), descA.data(),
      w_STAR_STAR.Buffer(), 
      QBlock.Buffer(), descQ.data(),
      ctrl.hessSchurCtrl.fullTriangle,
      accumulateSchurVecs,
      ctrl.hessSchurCtrl.scalapackAED );
    if( ctrl.time && grid.Rank() == 0 )
        Output("  scalapack::HessenbergSchur: ",timer.Stop()," seconds");
    if( ctrl.time && grid.Rank() == 0 )
        timer.Start();
    A = ABlock;
    Q = QBlock;
    if( ctrl.time && grid.Rank() == 0 )
        Output
        ("  Redist of (A,Q) from block to elemental: ",timer.Stop()," seconds");

    // TODO: Cache context, handle, and exit BLACS during El::Finalize()
    blacs::FreeGrid( context );
    blacs::FreeHandle( bHandle );
#endif
    */
    Copy( w_STAR_STAR, w );

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
