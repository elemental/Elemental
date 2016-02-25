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
inline void
QR
( Matrix<F>& A,
  Matrix<Complex<Base<F>>>& w,
  bool fullTriangle )
{
    DEBUG_ONLY(CSE cse("schur::QR"))
    const Int n = A.Height();
    w.Resize( n, 1 );
    lapack::Schur( n, A.Buffer(), A.LDim(), w.Buffer(), fullTriangle );
    if( IsComplex<F>::value )
        MakeTrapezoidal( UPPER, A );
    else
    {
        MakeTrapezoidal( UPPER, A, -1 );
        DEBUG_ONLY(CheckRealSchur(A))
    }
}

template<typename F>
inline void
QR
( Matrix<F>& A,
  Matrix<Complex<Base<F>>>& w,
  Matrix<F>& Q, 
  bool fullTriangle )
{
    DEBUG_ONLY(CSE cse("schur::QR"))
    const Int n = A.Height();
    Q.Resize( n, n );
    w.Resize( n, 1 );
    lapack::Schur
    ( n, A.Buffer(), A.LDim(), w.Buffer(), Q.Buffer(), Q.LDim(), fullTriangle );
    if( IsComplex<F>::value )
        MakeTrapezoidal( UPPER, A );
    else
    {
        MakeTrapezoidal( UPPER, A, -1 );
        DEBUG_ONLY(CheckRealSchur(A))
    }
}

template<typename F>
inline void
QR
( DistMatrix<F,MC,MR,BLOCK>& A,
  ElementalMatrix<Complex<Base<F>>>& w,
  bool fullTriangle,
  const HessQRCtrl& ctrl )
{
    DEBUG_ONLY(CSE cse("schur::QR"))
    AssertScaLAPACKSupport();
#ifdef EL_HAVE_SCALAPACK
    const Int n = A.Height();
    const int bHandle = blacs::Handle( A );
    const int context = blacs::GridInit( bHandle, A );
    auto descA = FillDesc( A, context );

    // Reduce the matrix to upper-Hessenberg form in an elemental form
    DistMatrix<F> AElem( A );
    DistMatrix<F,STAR,STAR> t( A.Grid() );
    Hessenberg( UPPER, AElem, t );
    MakeTrapezoidal( UPPER, AElem, -1 );
    A = AElem;

    // Run the QR algorithm in block form
    DistMatrix<Complex<Base<F>>,STAR,STAR> w_STAR_STAR( n, 1, A.Grid() );

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
inline void
QR
( DistMatrix<F,MC,MR,BLOCK>& A,
  ElementalMatrix<Complex<Base<F>>>& w,
  DistMatrix<F,MC,MR,BLOCK>& Q,
  bool fullTriangle, const HessQRCtrl& ctrl )
{
    DEBUG_ONLY(CSE cse("schur::QR"))
    AssertScaLAPACKSupport();
#ifdef EL_HAVE_SCALAPACK
    const Int n = A.Height();
    const int bHandle = blacs::Handle( A );
    const int context = blacs::GridInit( bHandle, A );
    Q.AlignWith( A );
    Q.Resize( n, n, A.LDim() );
    auto descA = FillDesc( A, context );
    auto descQ = FillDesc( Q, context );

    // Reduce A to upper-Hessenberg form in an element-wise distribution
    // and form the explicit reflector matrix
    DistMatrix<F> AElem( A ), QElem( A.Grid() );
    DistMatrix<F,STAR,STAR> t( A.Grid() );
    Hessenberg( UPPER, AElem, t );
    // There is not yet a 'form Q'
    Identity( QElem, n, n ); 
    hessenberg::ApplyQ( LEFT, UPPER, NORMAL, AElem, t, QElem );
    MakeTrapezoidal( UPPER, AElem, -1 );
    A = AElem;
    Q = QElem;
    
    // Compute the Schur decomposition in block form, multiplying the 
    // accumulated Householder reflectors from the right
    DistMatrix<Complex<Base<F>>,STAR,STAR> w_STAR_STAR( n, 1, A.Grid() );
    const bool multiplyQ = true;
    scalapack::HessenbergSchur
    ( n,
      A.Buffer(), descA.data(),
      w_STAR_STAR.Buffer(), 
      Q.Buffer(), descQ.data(),
      fullTriangle, multiplyQ, ctrl.distAED );
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
inline void
QR
( ElementalMatrix<F>& APre,
  ElementalMatrix<Complex<Base<F>>>& w, 
  bool fullTriangle, const HessQRCtrl& ctrl )
{
    DEBUG_ONLY(CSE cse("schur::QR"))
    AssertScaLAPACKSupport();

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.Get();

#ifdef EL_HAVE_SCALAPACK
    // Reduce the matrix to upper-Hessenberg form in an elemental form
    DistMatrix<F,STAR,STAR> t( A.Grid() );
    Hessenberg( UPPER, A, t );
    MakeTrapezoidal( UPPER, A, -1 );

    // Run the QR algorithm in block form
    // TODO: Create schur::HessenbergQR
    const Int n = A.Height(); 
    const Int mb = ctrl.blockHeight;
    const Int nb = ctrl.blockWidth;
    DistMatrix<F,MC,MR,BLOCK> ABlock( n, n, A.Grid(), mb, nb );
    ABlock = A;
    const int bHandle = blacs::Handle( ABlock );
    const int context = blacs::GridInit( bHandle, ABlock );
    blacs::Desc descA = FillDesc( ABlock, context );
    DistMatrix<Complex<Base<F>>,STAR,STAR> w_STAR_STAR( n, 1, A.Grid() );

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

    A = ABlock;
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
inline void
QR
( ElementalMatrix<F>& APre,
  ElementalMatrix<Complex<Base<F>>>& w, 
  ElementalMatrix<F>& QPre,
  bool fullTriangle,
  const HessQRCtrl& ctrl )
{
    DEBUG_ONLY(CSE cse("schur::QR"))
    AssertScaLAPACKSupport();

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    DistMatrixWriteProxy<F,F,MC,MR> QProx( QPre );
    auto& A = AProx.Get();
    auto& Q = QProx.Get();

#ifdef EL_HAVE_SCALAPACK
    const Int n = A.Height();
    // Reduce A to upper-Hessenberg form in an element-wise distribution
    // and form the explicit reflector matrix
    DistMatrix<F,STAR,STAR> t( A.Grid() );
    Hessenberg( UPPER, A, t );
    // There is not yet a 'form Q'
    Identity( Q, n, n ); 
    hessenberg::ApplyQ( LEFT, UPPER, NORMAL, A, t, Q );
    MakeTrapezoidal( UPPER, A, -1 );

    // Run the Hessenberg QR algorithm in block form
    const Int mb = ctrl.blockHeight;
    const Int nb = ctrl.blockWidth;
    DistMatrix<F,MC,MR,BLOCK>
      ABlock( n, n, A.Grid(), mb, nb ), 
      QBlock( n, n, A.Grid(), mb, nb );
    ABlock = A;
    QBlock = Q;
    const int bHandle = blacs::Handle( ABlock );
    const int context = blacs::GridInit( bHandle, ABlock );
    auto descA = FillDesc( ABlock, context );
    auto descQ = FillDesc( QBlock, context );

    // Compute the Schur decomposition in block form, multiplying the 
    // accumulated Householder reflectors from the right
    DistMatrix<Complex<Base<F>>,STAR,STAR> w_STAR_STAR( n, 1, A.Grid() );
    const bool multiplyQ = true;
    scalapack::HessenbergSchur
    ( n,
      ABlock.Buffer(), descA.data(),
      w_STAR_STAR.Buffer(), 
      QBlock.Buffer(), descQ.data(),
      fullTriangle, multiplyQ, ctrl.distAED );
    A = ABlock;
    Q = QBlock;
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
