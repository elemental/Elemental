/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

// See Chapter 5 of Nicholas J. Higham's "Functions of Matrices: Theory and
// Computation", which is currently available at:
// http://www.siam.org/books/ot104/OT104HighamChapter5.pdf

namespace El {

namespace sign {

template<typename Field>
void
NewtonStep
( const Matrix<Field>& X,
        Matrix<Field>& XNew,
  SignScaling scaling=SIGN_SCALE_FROB )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;

    // Calculate mu while forming XNew := inv(X)
    Real mu=1;
    Permutation P;
    XNew = X;
    LU( XNew, P );
    if( scaling == SIGN_SCALE_DET )
    {
        SafeProduct<Field> det = det::AfterLUPartialPiv( XNew, P );
        mu = Real(1)/Exp(det.kappa);
    }
    inverse::AfterLUPartialPiv( XNew, P );
    if( scaling == SIGN_SCALE_FROB )
        mu = Sqrt( FrobeniusNorm(XNew)/FrobeniusNorm(X) );

    // Overwrite XNew with the new iterate
    const Real halfMu = mu/Real(2);
    const Real halfMuInv = Real(1)/(2*mu);
    XNew *= halfMuInv;
    Axpy( halfMu, X, XNew );
}

template<typename Field>
void
NewtonStep
( const DistMatrix<Field>& X,
        DistMatrix<Field>& XNew,
  SignScaling scaling=SIGN_SCALE_FROB )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;

    // Calculate mu while forming B := inv(X)
    Real mu=1;
    DistPermutation P( X.Grid() );
    XNew = X;
    LU( XNew, P );
    if( scaling == SIGN_SCALE_DET )
    {
        SafeProduct<Field> det = det::AfterLUPartialPiv( XNew, P );
        mu = Real(1)/Exp(det.kappa);
    }
    inverse::AfterLUPartialPiv( XNew, P );
    if( scaling == SIGN_SCALE_FROB )
        mu = Sqrt( FrobeniusNorm(XNew)/FrobeniusNorm(X) );

    // Overwrite XNew with the new iterate
    const Real halfMu = mu/Real(2);
    const Real halfMuInv = Real(1)/(2*mu);
    XNew *= halfMuInv;
    Axpy( halfMu, X, XNew );
}

template<typename Field>
void
NewtonSchulzStep
( const Matrix<Field>& X,
        Matrix<Field>& XTmp,
        Matrix<Field>& XNew )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Int n = X.Height();

    // XTmp := 3I - X^2
    Identity( XTmp, n, n );
    Gemm( NORMAL, NORMAL, Real(-1), X, X, Real(3), XTmp );

    // XNew := 1/2 X XTmp
    Gemm( NORMAL, NORMAL, Real(1)/Real(2), X, XTmp, XNew );
}

template<typename Field>
void
NewtonSchulzStep
( const DistMatrix<Field>& X,
        DistMatrix<Field>& XTmp,
        DistMatrix<Field>& XNew )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Int n = X.Height();

    // XTmp := 3I - X^2
    Identity( XTmp, n, n );
    Gemm( NORMAL, NORMAL, Real(-1), X, X, Real(3), XTmp );

    // XNew := 1/2 X XTmp
    Gemm( NORMAL, NORMAL, Real(1)/Real(2), X, XTmp, XNew );
}

// Please see Chapter 5 of Higham's
// "Functions of Matrices: Theory and Computation" for motivation behind
// the different choices of p, which are usually in {0,1,2}
template<typename Field>
Int
Newton( Matrix<Field>& A, const SignCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    Real tol = ctrl.tol;
    if( tol == Real(0) )
        tol = A.Height()*limits::Epsilon<Real>();

    Int numIts=0;
    Matrix<Field> B;
    Matrix<Field> *X=&A, *XNew=&B;
    while( numIts < ctrl.maxIts )
    {
        // Overwrite XNew with the new iterate
        NewtonStep( *X, *XNew, ctrl.scaling );

        // Use the difference in the iterates to test for convergence
        Axpy( Real(-1), *XNew, *X );
        const Real oneDiff = OneNorm( *X );
        const Real oneNew = OneNorm( *XNew );

        // Ensure that X holds the current iterate and break if possible
        ++numIts;
        std::swap( X, XNew );
        if( ctrl.progress )
            cout << "after " << numIts << " Newton iter's: "
                 << "oneDiff=" << oneDiff << ", oneNew=" << oneNew
                 << ", oneDiff/oneNew=" << oneDiff/oneNew << ", tol="
                 << tol << endl;
        if( oneDiff/oneNew <= Pow(oneNew,ctrl.power)*tol )
            break;
    }
    if( X != &A )
        A = *X;
    return numIts;
}

template<typename Field>
Int
Newton( DistMatrix<Field>& A, const SignCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    Real tol = ctrl.tol;
    if( tol == Real(0) )
        tol = A.Height()*limits::Epsilon<Real>();

    Int numIts=0;
    DistMatrix<Field> B( A.Grid() );
    DistMatrix<Field> *X=&A, *XNew=&B;
    while( numIts < ctrl.maxIts )
    {
        // Overwrite XNew with the new iterate
        NewtonStep( *X, *XNew, ctrl.scaling );

        // Use the difference in the iterates to test for convergence
        Axpy( Real(-1), *XNew, *X );
        const Real oneDiff = OneNorm( *X );
        const Real oneNew = OneNorm( *XNew );

        // Ensure that X holds the current iterate and break if possible
        ++numIts;
        std::swap( X, XNew );
        if( ctrl.progress && A.Grid().Rank() == 0 )
            cout << "after " << numIts << " Newton iter's: "
                 << "oneDiff=" << oneDiff << ", oneNew=" << oneNew
                 << ", oneDiff/oneNew=" << oneDiff/oneNew << ", tol="
                 << tol << endl;
        if( oneDiff/oneNew <= Pow(oneNew,ctrl.power)*tol )
            break;
    }
    if( X != &A )
        A = *X;
    return numIts;
}

// TODO: NewtonSchulzHybrid which estimates when || X^2 - I ||_2 < 1

} // namespace sign

template<typename Field>
void Sign( Matrix<Field>& A, const SignCtrl<Base<Field>> ctrl )
{
    EL_DEBUG_CSE
    sign::Newton( A, ctrl );
}

template<typename Field>
void Sign
( Matrix<Field>& A, Matrix<Field>& N, const SignCtrl<Base<Field>> ctrl )
{
    EL_DEBUG_CSE
    Matrix<Field> ACopy( A );
    sign::Newton( A, ctrl );
    Gemm( NORMAL, NORMAL, Field(1), A, ACopy, N );
}

template<typename Field>
void Sign( AbstractDistMatrix<Field>& APre, const SignCtrl<Base<Field>> ctrl )
{
    EL_DEBUG_CSE

    DistMatrixReadWriteProxy<Field,Field,MC,MR> AProx( APre );
    auto& A = AProx.Get();

    sign::Newton( A, ctrl );
}

template<typename Field>
void Sign
( AbstractDistMatrix<Field>& APre,
  AbstractDistMatrix<Field>& NPre,
  const SignCtrl<Base<Field>> ctrl )
{
    EL_DEBUG_CSE

    DistMatrixReadWriteProxy<Field,Field,MC,MR> AProx( APre );
    DistMatrixWriteProxy<Field,Field,MC,MR> NProx( NPre );
    auto& A = AProx.Get();
    auto& N = NProx.Get();

    DistMatrix<Field> ACopy( A );
    sign::Newton( A, ctrl );
    Gemm( NORMAL, NORMAL, Field(1), A, ACopy, N );
}

// The Hermitian sign decomposition is equivalent to the Hermitian polar
// decomposition... A = (U sgn(Lambda) U') (U sgn(Lambda)Lambda U')
//                    = (U sgn(Lambda) U') (U |Lambda| U')

// Even though sgn(lambda) isn't well-defined when lambda=0, we will extend it
// from the right so that the sign decomposition of a singular Hermitian matrix
// is a polar decomposition (which always exists).

template<typename Field>
void HermitianSign
( UpperOrLower uplo, Matrix<Field>& A, const HermitianEigCtrl<Field>& ctrl )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;

    // Get the EVD of A
    Matrix<Real> w;
    Matrix<Field> Q;
    auto ctrlMod( ctrl );
    ctrlMod.tridiagEigCtrl.sort = UNSORTED;
    HermitianEig( uplo, A, w, Q, ctrlMod );

    const Int n = A.Height();
    for( Int i=0; i<n; ++i )
    {
        const Real omega = w(i);
        if( omega >= 0 )
            w(i) = Real(1);
        else
            w(i) = Real(-1);
    }

    // Reform the Hermitian matrix with the modified eigenvalues
    HermitianFromEVD( uplo, A, w, Q );
}

template<typename Field>
void HermitianSign
( UpperOrLower uplo,
  Matrix<Field>& A,
  Matrix<Field>& N,
  const HermitianEigCtrl<Field>& ctrl )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;

    // Get the EVD of A
    Matrix<Real> w;
    Matrix<Field> Q;
    auto ctrlMod( ctrl );
    ctrlMod.tridiagEigCtrl.sort = UNSORTED;
    HermitianEig( uplo, A, w, Q, ctrlMod );

    const Int n = A.Height();
    Matrix<Real> wSgn( n, 1 ), wAbs( n, 1 );
    for( Int i=0; i<n; ++i )
    {
        const Real omega = w(i);
        if( omega >= 0 )
        {
            wSgn(i) = Real(1);
            wAbs(i) = omega;
        }
        else
        {
            wSgn(i) = Real(-1);
            wAbs(i) = -omega;
        }
    }

    // Form the Hermitian matrices with modified eigenvalues
    HermitianFromEVD( uplo, A, wSgn, Q );
    HermitianFromEVD( uplo, N, wAbs, Q );
}

template<typename Field>
void HermitianSign
( UpperOrLower uplo,
  AbstractDistMatrix<Field>& APre,
  const HermitianEigCtrl<Field>& ctrl )
{
    EL_DEBUG_CSE

    DistMatrixReadWriteProxy<Field,Field,MC,MR> AProx( APre );
    auto& A = AProx.Get();

    // Get the EVD of A
    typedef Base<Field> Real;
    const Grid& g = A.Grid();
    DistMatrix<Real,VR,STAR> w(g);
    DistMatrix<Field> Q(g);
    auto ctrlMod( ctrl );
    ctrlMod.tridiagEigCtrl.sort = UNSORTED;
    HermitianEig( uplo, A, w, Q, ctrlMod );

    const Int numLocalEigs = w.LocalHeight();
    for( Int iLoc=0; iLoc<numLocalEigs; ++iLoc )
    {
        const Real omega = w.GetLocal(iLoc,0);
        if( omega >= 0 )
            w.SetLocal(iLoc,0,Real(1));
        else
            w.SetLocal(iLoc,0,Real(-1));
    }

    // Reform the Hermitian matrix with the modified eigenvalues
    HermitianFromEVD( uplo, A, w, Q );
}

template<typename Field>
void HermitianSign
( UpperOrLower uplo,
  AbstractDistMatrix<Field>& APre,
  AbstractDistMatrix<Field>& NPre,
  const HermitianEigCtrl<Field>& ctrl )
{
    EL_DEBUG_CSE

    DistMatrixReadWriteProxy<Field,Field,MC,MR> AProx( APre );
    DistMatrixWriteProxy<Field,Field,MC,MR> NProx( NPre );
    auto& A = AProx.Get();
    auto& N = NProx.Get();

    // Get the EVD of A
    typedef Base<Field> Real;
    const Grid& g = A.Grid();
    DistMatrix<Real,VR,STAR> w(g);
    DistMatrix<Field> Q(g);
    auto ctrlMod( ctrl );
    ctrlMod.tridiagEigCtrl.sort = UNSORTED;
    HermitianEig( uplo, A, w, Q, ctrlMod );

    const Int n = A.Height();
    const Int numLocalEigs = w.LocalHeight();
    DistMatrix<Real,VR,STAR> wSgn(g), wAbs(g);
    wSgn.AlignWith( w );
    wAbs.AlignWith( w );
    wSgn.Resize( n, 1 );
    wAbs.Resize( n, 1 );
    for( Int iLoc=0; iLoc<numLocalEigs; ++iLoc )
    {
        const Real omega = w.GetLocal(iLoc,0);
        if( omega >= 0 )
        {
            wSgn.SetLocal(iLoc,0,Real(1));
            wAbs.SetLocal(iLoc,0,omega);
        }
        else
        {
            wSgn.SetLocal(iLoc,0,Real(-1));
            wAbs.SetLocal(iLoc,0,-omega);
        }
    }

    // Form the Hermitian matrix with the modified eigenvalues
    HermitianFromEVD( uplo, A, wSgn, Q );
    HermitianFromEVD( uplo, N, wAbs, Q );
}

#define PROTO(Field) \
  template void Sign \
  ( Matrix<Field>& A, const SignCtrl<Base<Field>> ctrl ); \
  template void Sign \
  ( AbstractDistMatrix<Field>& A, const SignCtrl<Base<Field>> ctrl ); \
  template void Sign \
  ( Matrix<Field>& A, Matrix<Field>& N, const SignCtrl<Base<Field>> ctrl ); \
  template void Sign \
  ( AbstractDistMatrix<Field>& A, AbstractDistMatrix<Field>& N, \
    const SignCtrl<Base<Field>> ctrl ); \
  template void HermitianSign \
  ( UpperOrLower uplo, Matrix<Field>& A, \
    const HermitianEigCtrl<Field>& ctrl ); \
  template void HermitianSign \
  ( UpperOrLower uplo, AbstractDistMatrix<Field>& A, \
    const HermitianEigCtrl<Field>& ctrl ); \
  template void HermitianSign \
  ( UpperOrLower uplo, Matrix<Field>& A, Matrix<Field>& N, \
    const HermitianEigCtrl<Field>& ctrl ); \
  template void HermitianSign \
  ( UpperOrLower uplo, AbstractDistMatrix<Field>& A, \
    AbstractDistMatrix<Field>& N, \
    const HermitianEigCtrl<Field>& ctrl );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
