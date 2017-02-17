/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_HERM_TRIDIAG_EIG_QR_HPP
#define EL_HERM_TRIDIAG_EIG_QR_HPP
namespace El {
namespace herm_tridiag_eig {
namespace qr {

template<typename Real,
         typename=EnableIf<IsReal<Real>>>
Real WilkinsonShift
( const Real& alpha00, const Real& alpha01, const Real& alpha11 )
{
    EL_DEBUG_CSE
    // Following Parlett's "The Symmetric Eigenvalue Problem" [CITATION],
    // the eigenvalue of [alpha00, alpha01; alpha01, alpha11] closest to
    // alpha00 should not be naively evaluated as
    //
    //   omega = (alpha00+alpha11)/2 - sgn(delta)*sqrt(delta^2 + alpha01^2),
    //
    // where delta = (alpha11-alpha00)/2. Instead, one can recognize that
    //
    //   (alpha00+alpha11)/2 = alpha00 + delta = alpha00 + sgn(delta)*|delta|,
    //
    // so that
    //
    //   omega = alpha00 - sgn(delta)*(sqrt(delta^2 + alpha01^2) - |delta|).
    //
    // And since
    //
    //   alpha01^2                         |delta|-sqrt(delta^2+alpha01^2)
    //   ------------------------------- * ------------------------------- =
    //   |delta|+sqrt(delta^2+alpha01^2)   |delta|-sqrt(delta^2+alpha01^2)
    //
    //   alpha01^2*(|delta|-sqrt(delta^2+alpha01^2))
    //   ------------------------------------------- =
    //         |delta|^2-(|delta|^2+alpha01^2)
    //
    //   alpha01^2*(|delta|-sqrt(delta^2+alpha01^2))
    // - ------------------------------------------- =
    //                  alpha01^2
    //
    //         sqrt(delta^2+alpha01^2) - |delta|,
    //
    //  omega = alpha00-sgn(delta)*alpha01^2/(|delta|+sqrt(delta^2+alpha01^2)).
    //
    // This latter formula is recommended by Parlett (see subsection 8.9,
    // "Residual bounds using Wilkinson's shift"), but LAPACK's {s,l}steqr
    // [CITATION] goes a step further and computes sqrt(delta^2+alpha01^2) as
    //
    //   sqrt(delta^2 + alpha01^2) = |alpha01|*sqrt((delta/alpha01)^2 + 1),
    //
    // which, defining gamma=delta/alpha01, implies that
    //
    //  omega = alpha00 - alpha01/(gamma+sgn(gamma)*sqrt(gamma^2+1)).
    //
    const Real gamma = (alpha11-alpha00)/(2*alpha01);
    const Real rho = SafeNorm( gamma, Real(1) );
    // Following Parlett, demand that sgn(0) = 1.
    const Real sgnGamma = Sgn(gamma,false);
    const Real omega = alpha00 - alpha01 / (gamma + sgnGamma*rho);
    return omega;
}

// Cf. EISPACK's imtql1 [CITATION], which is based upon
//
//  Augustin (Austin) Dubrulle,
//  "A short note on the implicit QL algorithm for symmetric tridiagonal
//   matrices", 1969,
//  [CITATION]
//
// but computes the square-root of the sum of squares more carefully. A nice
// overview of the history of Dubrulle's contribution is given in the
// "Cleve's Corner" article
//
//  Cleve Moler,
//  "Dubrulle creates a faster tridiagonal QR algorithm", 2015.
//  [CITATION]
//
// We follow the suit of LAPACK's {s,d}steqr by representing the 'q' and 'g'
// variables from Dubrulle's algorithm with the single variable 'g' and using
// 't' for both the safe computation of 'e_i' and for 't_i'.
//
// TODO(poulson): Introduce [winBeg,winEnd) to avoid parent allocation
template<typename Field>
void QLSweep
( Matrix<Base<Field>>& d,
  Matrix<Base<Field>>& e,
  Matrix<Base<Field>>& cList,
  Matrix<Base<Field>>& sList,
  Matrix<Field>& Q,
  const Base<Field>& shift,
  bool wantEigVecs )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Int n = d.Height();
    const Real zero(0), one(1), two(2);
    if( wantEigVecs )
    {
        cList.Resize( n-1, 1 );
        sList.Resize( n-1, 1 );
    }

    Real s(1), c(1);
    Real t, h, p;
    Real g = d(n-1) - shift;
    Real u = zero;

    for( Int j=n-2; j>=0; --j )
    {
        h = c*e(j);
        p = s*e(j);
        t = Givens( g, p, c, s );
        if( j != n-2 )
            e(j+1) = t;
        g = d(j+1) - u;
        t = (d(j)-g)*s + two*c*h;
        u = s*t;
        d(j+1) = g + u;
        g = c*t - h;

        if( wantEigVecs )
        {
            cList(j) = c;
            sList(j) = -s;
        }
    }
    d(0) -= u;
    e(0) = g;
    if( wantEigVecs )
    {
        ApplyGivensSequence
        ( RIGHT, VARIABLE_GIVENS_SEQUENCE, BACKWARD, cList, sList, Q );
    }
}

template<typename Field>
void QRSweep
( Matrix<Base<Field>>& d,
  Matrix<Base<Field>>& e,
  Matrix<Base<Field>>& cList,
  Matrix<Base<Field>>& sList,
  Matrix<Field>& Q,
  const Base<Field>& shift,
  bool wantEigVecs )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Int n = d.Height();
    const Real zero(0), one(1), two(2);
    cList.Resize( n-1, 1 );
    sList.Resize( n-1, 1 );

    Real s(1), c(1);
    Real t, h, p;
    Real g = d(0) - shift;
    Real u = zero;

    for( Int j=0; j<n-1; ++j )
    {
        h = c*e(j);
        p = s*e(j);
        t = Givens( g, p, c, s );
        if( j != 0 )
            e(j-1) = t;
        g = d(j) - u;
        t = (d(j+1)-g)*s + two*c*h;
        u = s*t;
        d(j) = g + u;
        g = c*t - h;

        if( wantEigVecs )
        {
            cList(j) = c;
            sList(j) = s;
        }
    }
    d(n-1) -= u;
    e(n-2) = g;
    if( wantEigVecs )
    {
        ApplyGivensSequence
        ( RIGHT, VARIABLE_GIVENS_SEQUENCE, FORWARD, cList, sList, Q );
    }
}

// TODO(poulson): Support for what Parlett calls "ultimate shifts" in
// subsubsection 8.14.5 of "The symmetric eigenvalue problem" [CITATION].
//
// TODO(poulson): Support for what Parlett calls "Saad's shifts" in
// subsubsection 8.14.5 of "The symmetric eigenvalue problem" [CITATION].
//
template<typename Field>
herm_tridiag_eig::QRInfo
Helper
( Matrix<Base<Field>>& d,
  Matrix<Base<Field>>& e,
  Matrix<Field>& Q,
  const HermitianTridiagEigCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Int n = d.Height();
    const Int mQ = Q.Height();
    herm_tridiag_eig::QRInfo info;

    if( n <= 1 )
        return info;

    const Real zero(0), one(1);
    const Real eps = limits::Epsilon<Real>();
    const Real epsSquared = eps*eps;
    const Real safeMin = limits::SafeMin<Real>();
    const Real safeMax = one / safeMin;

    // [CITATION] Cf. LAPACK's {s,d}steqr for the choice of these bounds
    const Real normMin = Sqrt(safeMin) / epsSquared;
    const Real normMax = Sqrt(safeMax) / Real(3);

    Matrix<Real> cList(n-1,1), sList(n-1,1);
    Matrix<Real> dSub, eSub;
    Matrix<Field> QSub;

    const Int maxIter = n*ctrl.qrCtrl.maxIterPerEig;
    Int winBeg = 0;
    while( winBeg < n )
    {
        if( winBeg > 0 )
            e(winBeg-1) = zero;
        Int subWinBeg = winBeg;
        Int subWinEnd = n;
        // End the subwindow at the first sufficiently small off-diagonal
        // (if one exists)
        for( Int j=subWinBeg; j<n-1; ++j )
        {
            const Real etaAbs = Abs(e(j));
            if( etaAbs == zero )
            {
                if( ctrl.progress )
                    Output("  Detected zero at e(",j,")");
                subWinEnd = j+1;
                break;
            }
            if( etaAbs <= (Sqrt(Abs(d(j)))*Sqrt(Abs(d(j+1))))*eps )
            {
                if( ctrl.progress )
                    Output("  Deflating e(",j,")");
                e(j) = zero;
                subWinEnd = j+1;
                break;
            }
        }

        // Once we finish the current subwindow, our window will begin at the
        // end of the subwindow
        winBeg = subWinEnd;

        if( subWinEnd == subWinBeg+1 )
            continue;

        View( dSub, d, IR(subWinBeg,subWinEnd), ALL );
        View( eSub, e, IR(subWinBeg,Min(subWinEnd,n-1)), ALL );
        const Real infNormSub = HermitianTridiagInfinityNorm( dSub, eSub );
        bool scaledDown=false, scaledUp=false;
        if( infNormSub == zero )
            continue;
        if( infNormSub > normMax )
        {
            if( ctrl.progress )
                Output("Scaling down from ",infNormSub," to ",normMax);
            scaledDown = true;
            SafeScaleHermitianTridiag( infNormSub, normMax, dSub, eSub );
        }
        else if( infNormSub < normMin )
        {
            if( ctrl.progress )
                Output("Scaling up from ",infNormSub," to ",normMin);
            scaledUp = true;
            SafeScaleHermitianTridiag( infNormSub, normMin, dSub, eSub );
        }

        // Follow LAPACK's suit and use the simplest possible test for grading
        bool useQR = ( Abs(d(subWinEnd-1)) < Abs(d(subWinBeg)) );
        if( !useQR )
        {
            // QL iteration
            if( ctrl.progress )
            {
                Output
                ("QL iteration at iter ",info.numIterations,
                 " over [",subWinEnd,",",subWinEnd,")");
            }
            while( true )
            {
                Int iterEnd = subWinEnd;
                if( subWinBeg != subWinEnd-1 )
                {
                    for( Int j=subWinBeg; j<subWinEnd-1; ++j )
                    {
                        const Real etaAbs = Abs(e(j));
                        const Real etaAbsSquared = etaAbs*etaAbs;
                        if( etaAbsSquared <=
                            (epsSquared*Abs(d(j)))*Abs(d(j+1))+safeMin )
                        {
                            if( ctrl.progress )
                                Output("  Deflating e(",j,")=",e(j));
                            e(j) = zero;
                            iterEnd = j+1;
                            break;
                        }
                    }
                }

                if( subWinBeg+1 == iterEnd )
                {
                    // Deflate the isolated eigenvalue
                    if( ctrl.progress )
                        Output("  Deflating eigenvalue at ",subWinBeg);
                    subWinBeg += 1;
                    if( subWinBeg < subWinEnd )
                        continue;
                    else
                        break;
                }
                else if( subWinBeg+2 == iterEnd )
                {
                    // Deflate the isolated pair of eigenvalues
                    if( ctrl.progress )
                        Output("  Deflating eigenvalue pair at ",subWinBeg);
                    Real lambda0, lambda1;
                    if( ctrl.wantEigVecs )
                    {
                        Real c, s;
                        herm_eig::TwoByTwo
                        ( d(subWinBeg), e(subWinBeg), d(subWinBeg+1),
                          lambda0, lambda1, c, s,
                          ctrl.qrCtrl.fullAccuracyTwoByTwo );
                        // Apply the Givens rotation from the right to Q
                        blas::Rot
                        ( mQ, &Q(0,subWinBeg), 1, &Q(0,subWinBeg+1), 1, c, s );
                    }
                    else
                    {
                        herm_eig::TwoByTwo
                        ( d(subWinBeg), e(subWinBeg), d(subWinBeg+1),
                          lambda0, lambda1, ctrl.qrCtrl.fullAccuracyTwoByTwo );
                    }
                    d(subWinBeg) = lambda0;
                    d(subWinBeg+1) = lambda1;
                    e(subWinBeg) = zero;
                    subWinBeg += 2;
                    if( subWinBeg < subWinEnd )
                        continue;
                    else
                        break;
                }

                if( info.numIterations == maxIter )
                {
                    break;
                }
                ++info.numIterations;

                // TODO(poulson): Decide if it is worthwhile to avoid the cost
                // of these views
                View( dSub, d, IR(subWinBeg,iterEnd), ALL );
                View( eSub, e, IR(subWinBeg,Min(iterEnd,n-1)), ALL );
                if( ctrl.wantEigVecs )
                {
                    View( QSub, Q, ALL, IR(subWinBeg,iterEnd) );
                }

                Real shift = WilkinsonShift( dSub(0), eSub(0), dSub(1) );
                QLSweep
                ( dSub, eSub, cList, sList, QSub, shift, ctrl.wantEigVecs );
            }
        }
        else
        {
            // QR iteration
            if( ctrl.progress )
            {
                Output
                ("QR iteration at iter ",info.numIterations,
                 " over [",subWinBeg,",",subWinEnd,")");
            }
            while( true )
            {
                Int iterBeg = subWinBeg;
                if( subWinBeg != subWinEnd-1 )
                {
                    for( Int j=subWinEnd-1; j>=subWinBeg+1; --j )
                    {
                        const Real etaAbs = Abs(e(j-1));
                        const Real etaAbsSquared = etaAbs*etaAbs;
                        if( etaAbsSquared <=
                            (epsSquared*Abs(d(j)))*Abs(d(j-1))+safeMin )
                        {
                            if( ctrl.progress )
                                Output("  Deflating e(",j-1,")=",e(j-1));
                            e(j-1) = zero;
                            iterBeg = j;
                            break;
                        }
                    }
                }

                if( iterBeg+1 == subWinEnd )
                {
                    // Deflate the isolated eigenvalue
                    if( ctrl.progress )
                        Output("  Deflating eigenvalue at ",iterBeg);
                    subWinEnd -= 1;
                    if( subWinEnd > subWinBeg )
                        continue;
                    else
                        break;
                }
                else if( iterBeg+2 == subWinEnd )
                {
                    // Deflate the isolated pair of eigenvalues
                    if( ctrl.progress )
                        Output("  Deflating eigenvalue pair at ",iterBeg);
                    Real lambda0, lambda1;
                    if( ctrl.wantEigVecs )
                    {
                        Real c, s;
                        herm_eig::TwoByTwo
                        ( d(subWinEnd-2), e(subWinEnd-2), d(subWinEnd-1),
                          lambda0, lambda1, c, s,
                          ctrl.qrCtrl.fullAccuracyTwoByTwo );
                        // Apply the Givens rotation from the right to Q
                        blas::Rot
                        ( mQ, &Q(0,subWinEnd-2), 1, &Q(0,subWinEnd-1), 1,
                          c, s );
                    }
                    else
                    {
                        herm_eig::TwoByTwo
                        ( d(subWinEnd-2), e(subWinEnd-2), d(subWinEnd-1),
                          lambda0, lambda1, ctrl.qrCtrl.fullAccuracyTwoByTwo );
                    }
                    d(subWinEnd-2) = lambda0;
                    d(subWinEnd-1) = lambda1;
                    e(subWinEnd-2) = zero;
                    subWinEnd -= 2;
                    if( subWinEnd > subWinBeg )
                        continue;
                    else
                        break;
                }

                if( info.numIterations == maxIter )
                {
                    break;
                }
                ++info.numIterations;

                // TODO(poulson): Decide if it is worthwhile to avoid the cost
                // of these views
                View( dSub, d, IR(iterBeg,subWinEnd), ALL );
                View( eSub, e, IR(iterBeg,Min(subWinEnd,n-1)), ALL );
                if( ctrl.wantEigVecs )
                {
                    View( QSub, Q, ALL, IR(iterBeg,subWinEnd) );
                }

                Real shift =
                  WilkinsonShift
                  ( d(subWinEnd-1), e(subWinEnd-2), d(subWinEnd-2) );
                QRSweep
                ( dSub, eSub, cList, sList, QSub, shift, ctrl.wantEigVecs );
            }
        }

        if( scaledDown )
        {
            SafeScaleHermitianTridiag( normMax, infNormSub, dSub, eSub );
        }
        else if( scaledUp )
        {
            SafeScaleHermitianTridiag( normMin, infNormSub, dSub, eSub );
        }
        if( info.numIterations >= maxIter )
        {
            for( Int i=0; i<n-1; ++i )
                if( e(i) != zero )
                    ++info.numUnconverged;
            if( ctrl.qrCtrl.demandConverged )
                RuntimeError
                (info.numUnconverged," eigenvalues did not converge");
            return info;
        }
    }

    return info;
}

} // namespace qr

// TODO(poulson): Lift these routines up somewhere else?
void AllocatePackedQRInfo( vector<Int>& packedQRInfo )
{
    EL_DEBUG_CSE
    packedQRInfo.resize( 2 );
}

void PackQRInfo
( const herm_tridiag_eig::QRInfo& qrInfo, vector<Int>& packedQRInfo )
{
    EL_DEBUG_CSE
    if( packedQRInfo.size() != 2 )
        LogicError("Expected packedQRInfo to be of size 2");
    Int offset = 0;
    packedQRInfo[offset++] = qrInfo.numUnconverged;
    packedQRInfo[offset++] = qrInfo.numIterations;
}

void UnpackQRInfo
( const vector<Int>& packedQRInfo, herm_tridiag_eig::QRInfo& qrInfo )
{
    EL_DEBUG_CSE
    if( packedQRInfo.size() != 2 )
        LogicError("Expected packedQRInfo to be of size 2");
    Int offset = 0;
    qrInfo.numUnconverged = packedQRInfo[offset++];
    qrInfo.numIterations = packedQRInfo[offset++];
}

template<typename Real,typename=EnableIf<IsReal<Real>>>
herm_tridiag_eig::QRInfo
QRAlg
( Matrix<Real>& mainDiag,
  Matrix<Real>& subDiag,
  const HermitianTridiagEigCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    auto ctrlMod( ctrl );
    ctrlMod.wantEigVecs = false;
    Matrix<Real> Q;
    return qr::Helper( mainDiag, subDiag, Q, ctrlMod );
}

template<typename Real,typename=EnableIf<IsReal<Real>>>
herm_tridiag_eig::QRInfo
QRAlg
( AbstractDistMatrix<Real>& mainDiagPre,
  AbstractDistMatrix<Real>& subDiagPre,
  const HermitianTridiagEigCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE

    DistMatrixReadWriteProxy<Real,Real,STAR,STAR> mainDiagProx( mainDiagPre );
    DistMatrixReadProxy<Real,Real,STAR,STAR> subDiagProx( subDiagPre );
    auto& mainDiag = mainDiagProx.Get();
    auto& subDiag = subDiagProx.Get();

    auto ctrlMod( ctrl );
    ctrlMod.wantEigVecs = false;
    Matrix<Real> QLoc;

    herm_tridiag_eig::QRInfo info;
    if( ctrl.qrCtrl.broadcast )
    {
        const Grid& grid = mainDiag.Grid();
        vector<Int> packedQRInfo;
        AllocatePackedQRInfo( packedQRInfo );
        if( grid.VCRank() == 0 )
        {
            info =
              qr::Helper( mainDiag.Matrix(), subDiag.Matrix(), QLoc, ctrlMod );
            PackQRInfo( info, packedQRInfo );
        }
        El::Broadcast( mainDiag.Matrix(), grid.VCComm(), 0 );
        mpi::Broadcast
        ( packedQRInfo.data(), packedQRInfo.size(), 0, grid.VCComm() );
        UnpackQRInfo( packedQRInfo, info );
    }
    else
    {
        // Let's cross our fingers and ignore the forward instability
        info = qr::Helper( mainDiag.Matrix(), subDiag.Matrix(), QLoc, ctrlMod );
    }
    return info;
}

template<typename Real>
herm_tridiag_eig::QRInfo
QRAlg
( Matrix<Real>& mainDiag,
  Matrix<Real>& subDiag,
  Matrix<Real>& Q,
  const HermitianTridiagEigCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Int n = mainDiag.Height();
    auto ctrlMod( ctrl );
    ctrlMod.wantEigVecs = true;
    if( ctrl.accumulateEigVecs )
    {
        if( Q.Width() != n )
            LogicError("Q was an invalid size");
    }
    else
    {
        // Allow for computing a subset of rows of Q
        Identity( Q, n, n );
    }
    return qr::Helper( mainDiag, subDiag, Q, ctrlMod );
}

template<typename Real>
herm_tridiag_eig::QRInfo
QRAlg
( Matrix<Real>& mainDiag,
  Matrix<Real>& subDiag,
  Matrix<Complex<Real>>& Q,
  const HermitianTridiagEigCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Int n = mainDiag.Height();
    auto ctrlMod( ctrl );
    ctrlMod.wantEigVecs = true;
    if( ctrl.accumulateEigVecs )
    {
        if( Q.Width() != n )
            LogicError("Q was an invalid size");
        return qr::Helper( mainDiag, subDiag, Q, ctrlMod );
    }
    else
    {
        // Allow for computing a subset of rows of Q
        Matrix<Real> QReal;
        Identity( QReal, n, n );
        auto info = qr::Helper( mainDiag, subDiag, QReal, ctrlMod );
        Q = QReal;
        return info;
    }
}

template<typename Real>
herm_tridiag_eig::QRInfo
QRAlg
( AbstractDistMatrix<Real>& mainDiagPre,
  AbstractDistMatrix<Real>& subDiagPre,
  AbstractDistMatrix<Real>& QPre,
  const HermitianTridiagEigCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Int n = mainDiagPre.Height();

    DistMatrixReadWriteProxy<Real,Real,STAR,STAR> mainDiagProx( mainDiagPre );
    DistMatrixReadProxy<Real,Real,STAR,STAR> subDiagProx( subDiagPre );
    auto& mainDiag = mainDiagProx.Get();
    auto& subDiag = subDiagProx.Get();

    if( ctrl.accumulateEigVecs )
    {
        DistMatrixReadWriteProxy<Real,Real,VC,STAR> QProx( QPre );
        auto& Q = QProx.Get();
        if( Q.Width() != n )
            LogicError("Q was an invalid size");

        // WARNING: Forward instability can easily yield non-determinism and
        // lead to the following trivial parallelization yielding nonsensical
        // results. However, such an issue appears to be quite rare.
        return
          qr::Helper( mainDiag.Matrix(), subDiag.Matrix(), Q.Matrix(), ctrl );
    }
    else
    {
        DistMatrixWriteProxy<Real,Real,VC,STAR> QProx(QPre);
        auto& Q = QProx.Get();
        Identity( Q, n, n );

        // WARNING: Forward instability can easily yield non-determinism and
        // lead to the following trivial parallelization yielding nonsensical
        // results. However, such an issue appears to be quite rare.
        return
          qr::Helper( mainDiag.Matrix(), subDiag.Matrix(), Q.Matrix(), ctrl );
    }
}

template<typename Real>
herm_tridiag_eig::QRInfo
QRAlg
( AbstractDistMatrix<Real>& mainDiagPre,
  AbstractDistMatrix<Real>& subDiagPre,
  AbstractDistMatrix<Complex<Real>>& QPre,
  const HermitianTridiagEigCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Int n = mainDiagPre.Height();
    typedef Complex<Real> Field;

    DistMatrixReadWriteProxy<Real,Real,STAR,STAR> mainDiagProx( mainDiagPre );
    DistMatrixReadProxy<Real,Real,STAR,STAR> subDiagProx( subDiagPre );
    auto& mainDiag = mainDiagProx.Get();
    auto& subDiag = subDiagProx.Get();

    if( ctrl.accumulateEigVecs )
    {
        DistMatrixReadWriteProxy<Field,Field,VC,STAR> QProx( QPre );
        auto& Q = QProx.Get();
        if( Q.Width() != n )
            LogicError("Q was an invalid size");

        // WARNING: Forward instability can easily yield non-determinism and
        // lead to the following trivial parallelization yielding nonsensical
        // results. However, such an issue appears to be quite rare.
        return
          qr::Helper( mainDiag.Matrix(), subDiag.Matrix(), Q.Matrix(), ctrl );
    }
    else
    {
        DistMatrix<Real,VC,STAR> QReal(QPre.Grid());
        Identity( QReal, n, n );

        // WARNING: Forward instability can easily yield non-determinism and
        // lead to the following trivial parallelization yielding nonsensical
        // results. However, such an issue appears to be quite rare.
        auto info =
          qr::Helper
          ( mainDiag.Matrix(), subDiag.Matrix(), QReal.Matrix(), ctrl );

        Copy( QReal, QPre );
        return info;
    }
}

} // namespace herm_tridiag_eig
} // namespace El

#endif // ifndef EL_HERM_TRIDIAG_EIG_QR_HPP
