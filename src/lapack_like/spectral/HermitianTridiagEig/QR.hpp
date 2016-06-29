/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SPECTRAL_HERM_TRIDIAG_EIG_QR_HPP
#define EL_SPECTRAL_HERM_TRIDIAG_EIG_QR_HPP
namespace El {
namespace herm_eig {

template<typename Real,typename=EnableIf<IsReal<Real>>>
Real WilkinsonShift
( const Real& alpha00, const Real& alpha01, const Real& alpha11 )
{
    DEBUG_CSE
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
    bool symmetric = false;
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
template<typename Real,typename=EnableIf<IsReal<Real>>>
void QLSweep
( Matrix<Real>& d,
  Matrix<Real>& e,
  Matrix<Real>& cList,
  Matrix<Real>& sList,
  Matrix<Real>& Q,
  const Real& shift,
  bool wantEigVecs )
{
    DEBUG_CSE
    const Int n = d.Height();
    const Real zero(0), one(1), two(2);
    cList.Resize( n-1, 1 );
    sList.Resize( n-1, 1 );

    Real s(1), c(1);
    Real t, h, p;
    Real g = d(n-1) - shift;
    Real u = zero;

    for( Int j=n-2; j>=0; --j )
    {
        h = c*e(j);
        p = s*e(j);
        t = SafeNorm( g, p, c, s );
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
            sList(j) = s;
        }
    }
    d(0) -= u;
    e(0) = g;
    if( wantEigVecs )
    {
        ApplyGivensSequence
        ( RIGHT, VARIABLE_GIVENS_SEQUENCE, BACKWARD,
          Q.Height(), n, cList, sList, Q );
    }
}

template<typename Real,typename=EnableIf<IsReal<Real>>>
void QRSweep
( Matrix<Real>& d,
  Matrix<Real>& e,
  Matrix<Real>& cList,
  Matrix<Real>& sList,
  Matrix<Real>& Q,
  const Real& shift,
  bool wantEigVecs )
{
    DEBUG_CSE
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
        ( RIGHT, VARIABLE_GIVENS_SEQUENCE, FORWARD,
          Q.Height(), n, cList, sList, Q );
    }
}

// TODO(poulson): Support for what Parlett calls "ultimate shifts" in 
// subsubsection 8.14.5 of "The symmetric eigenvalue problem" [CITATION].
//
// TODO(poulson): Support for what Parlett calls "Saad's shifts" in 
// subsubsection 8.14.5 of "The symmetric eigenvalue problem" [CITATION].
//
template<typename Real,typename>
HermitianTridiagQRInfo
Helper
( Matrix<Real>& d,
  Matrix<Real>& e, 
  Matrix<Real>& w,
  Matrix<Real>& Q,
  const HermitianTridiagQRCtrl& ctrl )
{
    DEBUG_CSE
    const Int n = d.Height();
    HermitianTridiagQRInfo info;
    bool fullAccuracy = true; // TODO: Make configurable

    w.Resize( n, 1 );
    if( n == 0 )
        return;
    if( n == 1 )
    {
        w(0) = d(0);
        if( ctrl.wantEigVecs )
            Q(0,0) = 1;
        return;
    }

    const Real zero(0), one(1);
    const Real eps = limits::Epsilon<Real>();
    const Real epsSquared = eps*eps;
    const Real safeMin = limits::SafeMin<Real>();
    const Real safeMax = one / safeMin;

    // [CITATION] Cf. LAPACK's {s,d}steqr for the choice of these bounds
    const Real normMin = Sqrt(safeMin) / epsSquared;
    const Real normMax = Sqrt(safeMax) / Real(3);

    Matrix<Real> cList(n-1,1), sList(n-1,1);
    Matrix<Real> dSub, eSub, QSub;

    const Int maxIter = n*ctrl.maxIterPerEig; 
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
                subWinEnd = j+1;
                break;
            }
            if( etaAbs <= (Sqrt(Abs(d(j)))*Sqrt(Abs(d(j+1))))*eps )
            {
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
        const Real infNormSub = HermitianTridiagonalInfinityNorm( dSub, eSub );
        bool scaledDown=false, scaledUp=false;
        if( infNormSub == zero )
            continue;
        if( infNormSub > normMax )
        {
            scaledDown = true;
            SafeScaleHermitianTridiagonal( infNormSub, normMax, dSub, eSub );
        }
        else if( infNormSub < normMin )
        {
            scaledUp = true;
            SafeScaleHermitianTridiagonal( infNormSub, normMin, dSub, eSub );
        }

        // Follow LAPACK's suit and use the simplest possible test for grading
        bool useQR = ( Abs(d(subWinEnd-1)) < Abs(d(subWinBeg)) );
        if( subWinEnd > subWinBeg )
        {
            // QL iteration
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
                            e(j) = zero;
                            iterEnd = j+1;
                            break;
                        }
                    }
                }

                if( subWinBeg+1 == iterEnd )
                {
                    // Deflate the isolated eigenvalue
                    subWinBeg += 1;
                    if( subWinBeg < subWinEnd )
                        continue;
                    else
                        break;
                }
                else if( subWinBeg+2 == iterEnd )
                {
                    // Deflate the isolated pair of eigenvalues
                    Real lambda0, lambda1;
                    if( ctrl.wantEigVecs )
                    {
                        Real c, s;
                        TwoByTwo
                        ( d(subWinBeg), e(subWinBeg), d(subWinBeg+1),
                          lambda0, lambda1, c, s, fullAccuracy );
                        // Apply the Givens rotation from the right to Q
                        blas::Rot
                        ( n, &Q(0,subWinBeg), 1, &Q(0,subWinBeg+1), 1, c, s );
                    }
                    else
                    {
                        TwoByTwo
                        ( d(subWinBeg), e(subWinBeg), d(subWinBeg+1),
                          lambda0, lambda1, fullAccuracy );
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
                View( dSub, d, IR(subWinBeg,subWinEnd), ALL );
                View( eSub, e, IR(subWinBeg,Min(subWinEnd,n-1)), ALL );
                if( ctrl.wantEigVecs )
                {
                    View( QSub, Q, ALL, IR(subWinBeg,subWinEnd) );
                }

                Real shift = WilkinsonShift( dSub(0), eSub(0), dSub(1) );
                QLSweep
                ( dSub, eSub, cList, sList, QSub, shift, ctrl.wantEigVecs );
            }
        }
        else
        {
            // QR iteration
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
                            e(j-1) = zero;
                            iterBeg = j;
                            break;
                        }
                    }
                }

                if( iterBeg+1 == subWinEnd )
                {
                    // Deflate the isolated eigenvalue
                    subWinEnd -= 1;
                    if( subWinEnd > subWinBeg )
                        continue;
                    else
                        break;
                }
                else if( iterBeg+2 == subWinEnd )
                {
                    // Deflate the isolated pair of eigenvalues
                    Real lambda0, lambda1;
                    if( ctrl.wantEigVecs )
                    {
                        Real c, s;
                        TwoByTwo
                        ( d(subWinEnd-2), e(subWinEnd-2), d(subWinEnd-1),
                          lambda0, lambda1, c, s, fullAccuracy ); 
                        // Apply the Givens rotation from the right to Q
                        blas::Rot
                        ( n, &Q(0,subWinEnd-2), 1, &Q(0,subWinEnd-1), 1, c, s );
                    }
                    else
                    {
                        TwoByTwo
                        ( d(subWinEnd-2), e(subWinEnd-2), d(subWinEnd-1),
                          lambda0, lambda1, fullAccuracy );
                    }
                    d(subWinEnd-2) = lambda0;
                    d(subWinEnd-1) = lambda1;
                    e(subWinBeg) = zero;
                    subWinEnd -= 2;
                    if( subWinEnd < subWinBeg )
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
                View( dSub, d, IR(subWinBeg,subWinEnd), ALL );
                View( eSub, e, IR(subWinBeg,Min(subWinEnd,n-1)), ALL );
                if( ctrl.wantEigVecs )
                {
                    View( QSub, Q, ALL, IR(subWinBeg,subWinEnd) );
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
            if( ctrl.demandConverged )
                RuntimeError
                (info.numUnconverged," eigenvalues did not converge");
            return info;
        }
    }

    if( ctrl.wantEigVecs )
    {
        Sort( w, Q, ASCENDING );
    }
    else
    {
        Sort( w, ASCENDING );
    }

    return info;
}

template<typename Real,typename>
HermitianTridiagQRInfo
TridiagQR
( Matrix<Real>& mainDiag,
  Matrix<Real>& subDiag, 
  Matrix<Real>& w,
  const HermitianTridiagQRCtrl& ctrl )
{
    DEBUG_CSE
    auto ctrlMod( ctrl );
    ctrlMod.wantEigVecs = false;
    Matrix<Real> Q;
    return Helper( mainDiag, subDiag, w, Q, ctrlMod );
}

template<typename Real,typename>
HermitianTridiagQRInfo
TridiagQR
( Matrix<Real>& mainDiag,
  Matrix<Real>& subDiag,
  Matrix<Real>& w,
  Matrix<Real>& Q,
  const HermitianTridiagQRCtrl& ctrl )
{
    DEBUG_CSE
    const Int n = mainDiag.Height();
    auto ctrlMod( ctrl );
    ctrlMod.wantEigVecs = true;
    // TODO: Allow for computing a subset of eigenvectors
    Q.Resize( n, n );
    return Helper( mainDiag, subDiag, w, Q, ctrlMod );
}

} // namespace herm_eig
} // namespace El

#endif // ifndef EL_SPECTRAL_HERM_TRIDIAG_EIG_QR_HPP
