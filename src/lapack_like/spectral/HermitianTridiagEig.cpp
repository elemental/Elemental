/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

#include "./HermitianTridiagEig/QR.hpp"

// NOTE: dSubReal and QReal could be packed into their complex counterparts

namespace El {

// Return eigenvalues
// ==================

namespace herm_eig {

template<typename Real>
void SortAndFilter
( Matrix<Real>& w, const HermitianTridiagEigCtrl<Real>& ctrl )
{
    const Int n = w.Height();
    if( ctrl.subset.indexSubset )
    {
        Sort( w, ctrl.sort );
        auto wCopy = w;
        w = wCopy(IR(ctrl.subset.lowerIndex,ctrl.subset.upperIndex+1),ALL);
    }
    else if( ctrl.subset.rangeSubset )
    {
        Int numValid = 0;
        for( Int j=0; j<n; ++j )
            if( w(j) > ctrl.subset.lowerBound &&
                w(j) <= ctrl.subset.upperBound )
                ++numValid;

        Matrix<Real> wFilter(numValid,1);
        numValid = 0;
        for( Int j=0; j<n; ++j )
           if( w(j) > ctrl.subset.lowerBound &&
               w(j) <= ctrl.subset.upperBound )
               wFilter(numValid++) = w(j);

        w = wFilter;
        Sort( w, ctrl.sort );
    }
    else
    {
        Sort( w, ctrl.sort );
    }
}

template<typename F>
void SortAndFilter
( Matrix<Base<F>>& w,
  Matrix<F>& Q,
  const HermitianTridiagEigCtrl<Base<F>>& ctrl )
{
    const Int n = w.Height();
    if( ctrl.subset.indexSubset )
    {
        auto sortPairs = TaggedSort( w, ctrl.sort );
        for( Int j=0; j<n; ++j )
            w(j) = sortPairs[j].value;
        ApplyTaggedSortToEachRow( sortPairs, Q );

        auto wCopy = w;
        auto QCopy = Q;
        w = wCopy(IR(ctrl.subset.lowerIndex,ctrl.subset.upperIndex+1),ALL);
        Q = QCopy(ALL,IR(ctrl.subset.lowerIndex,ctrl.subset.upperIndex+1));
    }
    else if( ctrl.subset.rangeSubset )
    {
        Int numValid = 0;
        for( Int j=0; j<n; ++j )
            if( w(j) > ctrl.subset.lowerBound &&
                w(j) <= ctrl.subset.upperBound )
                ++numValid;

        // TODO(poulson): Avoid unnecessary extra matrix for Q
        Matrix<Base<F>> wFilter(numValid,1);
        Matrix<F> QFilter(n,numValid);
        numValid = 0;
        for( Int j=0; j<n; ++j )
        {
           if( w(j) > ctrl.subset.lowerBound &&
               w(j) <= ctrl.subset.upperBound )
           {
               wFilter(numValid) = w(j);
               auto qFilterCol = QFilter(ALL,IR(numValid));
               qFilterCol = Q(ALL,IR(j));
               ++numValid;
           }
        }

        w = wFilter;
        Q = QFilter;

        auto sortPairs = TaggedSort( w, ctrl.sort );
        for( Int j=0; j<numValid; ++j )
            w(j) = sortPairs[j].value;
        ApplyTaggedSortToEachRow( sortPairs, Q );
    }
    else
    {
        auto sortPairs = TaggedSort( w, ctrl.sort );
        for( Int j=0; j<n; ++j )
            w(j) = sortPairs[j].value;
        ApplyTaggedSortToEachRow( sortPairs, Q );
    }
}

} // namespace herm_eig

namespace herm_tridiag_eig {

template<typename Real,typename=EnableIf<IsBlasScalar<Real>>>
HermitianTridiagEigInfo
Helper
(       Matrix<Real>& d,
        Matrix<Real>& dSub,
        Matrix<Real>& w,
  const HermitianTridiagEigCtrl<Real>& ctrl )
{
    const Int n = d.Height();
    HermitianTridiagEigInfo info;

    if( ctrl.useQR )
    {
        w = d;    
        info.qrInfo = QRAlg( w, dSub, ctrl );
        herm_eig::SortAndFilter( w, ctrl );
    }
    else
    {
        w.Resize( n, 1 );
        if( ctrl.subset.rangeSubset )
        {
             const Int k = lapack::SymmetricTridiagEig
              ( BlasInt(n), d.Buffer(), dSub.Buffer(), w.Buffer(), 
                ctrl.subset.lowerBound, ctrl.subset.upperBound );
             w.Resize( k, 1 );
        }
        else if( ctrl.subset.indexSubset )
        {
            const Int numEig = ctrl.subset.upperIndex-ctrl.subset.lowerIndex+1;
            lapack::SymmetricTridiagEig
            ( BlasInt(n), d.Buffer(), dSub.Buffer(), w.Buffer(), 
              BlasInt(ctrl.subset.lowerIndex),
              BlasInt(ctrl.subset.upperIndex) );
            w.Resize( numEig, 1 );
        }
        else
            lapack::SymmetricTridiagEig
            ( BlasInt(n), d.Buffer(), dSub.Buffer(), w.Buffer() );
        Sort( w, ctrl.sort );
    }

    return info;
}

template<typename Real,typename=DisableIf<IsBlasScalar<Real>>,typename=void>
HermitianTridiagEigInfo
Helper
(       Matrix<Real>& d,
        Matrix<Real>& dSub,
        Matrix<Real>& w,
  const HermitianTridiagEigCtrl<Real>& ctrl )
{
    const Int n = d.Height();
    HermitianTridiagEigInfo info;

    w = d;    
    info.qrInfo = QRAlg( w, dSub, ctrl );
    herm_eig::SortAndFilter( w, ctrl );

    return info;
}

template<typename Real>
HermitianTridiagEigInfo
Helper
(       Matrix<Real>& d,
        Matrix<Complex<Real>>& dSub,
        Matrix<Real>& w,
  const HermitianTridiagEigCtrl<Real>& ctrl )
{
    typedef Complex<Real> C;
    const Int n = d.Height();
    Matrix<Real> dSubReal( n-1, 1 );
    C yLast = 1;
    for( Int j=0; j<n-1; ++j )
    {
        const C psi = dSub(j);
        const Real psiAbs = Abs(psi);
        if( psiAbs == Real(0) )
            yLast = 1;
        else
            yLast = ComplexFromPolar(Real(1),Arg(psi*yLast));
        dSubReal(j) = psiAbs;
    }
    return HermitianTridiagEig( d, dSubReal, w, ctrl );
}

} // namespace herm_tridiag_eig

template<typename F>
HermitianTridiagEigInfo
HermitianTridiagEig
(       Matrix<Base<F>>& d,
        Matrix<F>& dSub,
        Matrix<Base<F>>& w,
  const HermitianTridiagEigCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    return herm_tridiag_eig::Helper( d, dSub, w, ctrl );
}

namespace herm_tridiag_eig {

template<typename Real,typename=EnableIf<IsBlasScalar<Real>>>
HermitianTridiagEigInfo
Helper
( const AbstractDistMatrix<Real>& d,
  const AbstractDistMatrix<Real>& dSub,
        AbstractDistMatrix<Real>& wPre,
  const HermitianTridiagEigCtrl<Real>& ctrl )
{
    DEBUG_CSE
    if( ctrl.useQR )
    {
        LogicError("Distributed QR not yet supported");
    }
    HermitianTridiagEigInfo info;

    ElementalProxyCtrl wCtrl;
    wCtrl.colConstrain = true;
    wCtrl.colAlign = 0;

    DistMatrixWriteProxy<Real,Real,VR,STAR> wProx( wPre, wCtrl );
    auto& w = wProx.Get();

    // Force the computation to take place with double-precision since PMRRR
    // currently only supports this case
    const Int n = d.Height();
    const Grid& g = d.Grid();
    DistMatrix<double,STAR,STAR> d_STAR_STAR(g), dSub_STAR_STAR(g);
    Copy( d, d_STAR_STAR );
    dSub_STAR_STAR.Resize( n-1, 1, n );
    Copy( dSub, dSub_STAR_STAR );

    vector<double> wVector(n);
    herm_tridiag_eig::Info rangeInfo;
    if( ctrl.subset.rangeSubset )
        rangeInfo = herm_tridiag_eig::Eig
          ( int(n), d_STAR_STAR.Buffer(), dSub_STAR_STAR.Buffer(), 
            wVector.data(), w.ColComm(), 
            ctrl.subset.lowerBound, ctrl.subset.upperBound );
    else if( ctrl.subset.indexSubset )
        rangeInfo = herm_tridiag_eig::Eig
          ( int(n), d_STAR_STAR.Buffer(), dSub_STAR_STAR.Buffer(), 
            wVector.data(), w.ColComm(), 
            int(ctrl.subset.lowerIndex), int(ctrl.subset.upperIndex) );
    else
        rangeInfo = herm_tridiag_eig::Eig
          ( int(n), d_STAR_STAR.Buffer(), dSub_STAR_STAR.Buffer(), 
            wVector.data(), w.ColComm() );
    w.Resize( rangeInfo.numGlobalEigenvalues, 1 );
    for( Int iLoc=0; iLoc<w.LocalHeight(); ++iLoc )
        w.SetLocal( iLoc, 0, Real(wVector[iLoc]) );
    Sort( w, ctrl.sort );

    return info;
}

template<typename Real,typename=DisableIf<IsBlasScalar<Real>>,typename=void>
HermitianTridiagEigInfo
Helper
( const AbstractDistMatrix<Real>& d,
  const AbstractDistMatrix<Real>& dSub,
        AbstractDistMatrix<Real>& wPre,
  const HermitianTridiagEigCtrl<Real>& ctrl )
{
    DEBUG_CSE
    LogicError
    ("Distributed non-standard HermitianTridiagEig not yet supported");
    HermitianTridiagEigInfo info;
    return info;
}

template<typename Real,typename=EnableIf<IsBlasScalar<Real>>>
HermitianTridiagEigInfo
Helper
( const AbstractDistMatrix<Real         >& d,
  const AbstractDistMatrix<Complex<Real>>& dSub,
        AbstractDistMatrix<Real         >& wPre, 
  const HermitianTridiagEigCtrl<Real>& ctrl )
{
    DEBUG_CSE
    if( ctrl.useQR )
    {
        LogicError("Distributed QR not yet supported");
    }
    HermitianTridiagEigInfo info;

    ElementalProxyCtrl wCtrl;
    wCtrl.colConstrain = true;
    wCtrl.colAlign = 0;

    DistMatrixWriteProxy<Real,Real,VR,STAR> wProx( wPre, wCtrl );
    auto& w = wProx.Get();

    // Force the computation to take place with double-precision since PMRRR
    // currently only supports this case
    const Int n = d.Height();
    const Grid& g = d.Grid();
    DistMatrix<double,STAR,STAR> d_STAR_STAR(g);
    DistMatrix<Complex<double>,STAR,STAR> dSub_STAR_STAR(g);
    Copy( d, d_STAR_STAR );
    dSub_STAR_STAR.Resize( n-1, 1, n );
    Copy( dSub, dSub_STAR_STAR );
    auto& dSubLoc = dSub_STAR_STAR.Matrix();

    DistMatrix<double,STAR,STAR> dSubReal(g);
    dSubReal.Resize( n-1, 1, n );
    auto& dSubRealLoc = dSubReal.Matrix();

    Complex<double> yLast = 1;
    for( Int j=0; j<n-1; ++j )
    {
        const Complex<double> psi = dSubLoc(j);
        const double psiAbs = Abs(psi);
        if( psiAbs == double(0) )
            yLast = 1;
        else
            yLast = ComplexFromPolar(double(1),Arg(psi*yLast));
        dSubRealLoc(j) = psiAbs;
    }

    herm_tridiag_eig::Info rangeInfo;
    vector<double> wVector(n);
    if( ctrl.subset.rangeSubset )
    {
        rangeInfo = herm_tridiag_eig::Eig
          ( int(n), d_STAR_STAR.Buffer(), dSubReal.Buffer(),
            wVector.data(), w.ColComm(),
            ctrl.subset.lowerBound, ctrl.subset.upperBound );
    }
    else if( ctrl.subset.indexSubset )
    {
        rangeInfo = herm_tridiag_eig::Eig
          ( int(n), d_STAR_STAR.Buffer(), dSubReal.Buffer(),
            wVector.data(), w.ColComm(),
            int(ctrl.subset.lowerIndex), int(ctrl.subset.upperIndex) );
    }
    else
    {
        rangeInfo = herm_tridiag_eig::Eig
          ( int(n), d_STAR_STAR.Buffer(), dSubReal.Buffer(),
            wVector.data(), w.ColComm() );
    }
    w.Resize( rangeInfo.numGlobalEigenvalues, 1 );
    auto& wLoc = w.Matrix();
    for( Int iLoc=0; iLoc<w.LocalHeight(); ++iLoc )
        wLoc(iLoc) = Real(wVector[iLoc]);

    Sort( w, ctrl.sort );

    return info;
}

template<typename Real,typename=DisableIf<IsBlasScalar<Real>>,typename=void>
HermitianTridiagEigInfo
Helper
( const AbstractDistMatrix<Real         >& d,
  const AbstractDistMatrix<Complex<Real>>& dSub,
        AbstractDistMatrix<Real         >& wPre, 
  const HermitianTridiagEigCtrl<Real>& ctrl )
{
    DEBUG_CSE
    LogicError
    ("Distributed non-standard HermitianTridiagEig not yet supported");
    HermitianTridiagEigInfo info;
    return info;
}

} // namespace herm_tridiag_eig

template<typename F>
HermitianTridiagEigInfo
HermitianTridiagEig
( const AbstractDistMatrix<Base<F>>& d,
  const AbstractDistMatrix<F      >& dSub,
        AbstractDistMatrix<Base<F>>& w, 
  const HermitianTridiagEigCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    return herm_tridiag_eig::Helper( d, dSub, w, ctrl );
}

// Return eigenpairs
// =================

namespace herm_tridiag_eig {

template<typename Real,typename=EnableIf<IsBlasScalar<Real>>>
HermitianTridiagEigInfo
Helper
(       Matrix<Real>& d,
        Matrix<Real>& dSub,
        Matrix<Real>& w,
        Matrix<Real>& Q,
  const HermitianTridiagEigCtrl<Real>& ctrl )
{
    const Int n = d.Height();
    HermitianTridiagEigInfo info;

    if( ctrl.useQR )
    {
        w = d;    
        info.qrInfo = QRAlg( w, dSub, Q, ctrl );
        herm_eig::SortAndFilter( w, Q, ctrl );
    }
    else
    {
        w.Resize( n, 1 );
        if( ctrl.subset.rangeSubset )
        {
             Q.Resize( n, n );
             const Int k = lapack::SymmetricTridiagEig
              ( BlasInt(n), d.Buffer(), dSub.Buffer(), w.Buffer(), 
                Q.Buffer(), BlasInt(Q.LDim()),
                ctrl.subset.lowerBound, ctrl.subset.upperBound );
             w.Resize( k, 1 );
             Q.Resize( n, k );
        }
        else if( ctrl.subset.indexSubset )
        {
            const Int numEig = ctrl.subset.upperIndex-ctrl.subset.lowerIndex+1;
            Q.Resize( n, numEig );
            lapack::SymmetricTridiagEig
            ( BlasInt(n), d.Buffer(), dSub.Buffer(), w.Buffer(), 
              Q.Buffer(), BlasInt(Q.LDim()),
              BlasInt(ctrl.subset.lowerIndex),
              BlasInt(ctrl.subset.upperIndex) );
            w.Resize( numEig, 1 );
        }
        else
        {
            Q.Resize( n, n );
            lapack::SymmetricTridiagEig
            ( BlasInt(n), d.Buffer(), dSub.Buffer(), w.Buffer(), 
              Q.Buffer(), BlasInt(Q.LDim()) );
        }
        auto sortPairs = TaggedSort( w, ctrl.sort );
        for( Int j=0; j<n; ++j )
            w(j) = sortPairs[j].value;
        ApplyTaggedSortToEachRow( sortPairs, Q );
    }

    return info;
}

template<typename Real,typename=DisableIf<IsBlasScalar<Real>>,typename=void>
HermitianTridiagEigInfo
Helper
(       Matrix<Real>& d,
        Matrix<Real>& dSub,
        Matrix<Real>& w,
        Matrix<Real>& Q,
  const HermitianTridiagEigCtrl<Real>& ctrl )
{
    const Int n = d.Height();
    HermitianTridiagEigInfo info;

    w = d;    
    info.qrInfo = QRAlg( w, dSub, Q, ctrl );
    herm_eig::SortAndFilter( w, Q, ctrl );

    return info;
}

// (Y^H T Y) QHat = QHat Lambda
// T (Y QHat) = (Y QHat) Lambda
template<typename Real>
HermitianTridiagEigInfo
Helper
(       Matrix<Real>& d,
        Matrix<Complex<Real>>& dSub,
        Matrix<Real>& w, 
        Matrix<Complex<Real>>& Q,
  const HermitianTridiagEigCtrl<Real>& ctrl )
{
    typedef Complex<Real> C;
    const Int n = d.Height();
    HermitianTridiagEigInfo info;

    Matrix<Real> dSubReal( n-1, 1 );
    Matrix<C> y( n, 1 );
    y(0) = 1;
    for( Int j=0; j<n-1; ++j )
    {
        const C psi = dSub(j);
        const Real psiAbs = Abs(psi);
        if( psiAbs == Real(0) )
            y(j+1) = 1;
        else
            y(j+1) = ComplexFromPolar(Real(1),Arg(psi*y(j)));
        dSubReal(j) = psiAbs;
    }
    Matrix<Real> QReal;
    HermitianTridiagEig( d, dSubReal, w, QReal, ctrl );
    Q.Resize( n, QReal.Width() );
    for( Int j=0; j<QReal.Width(); ++j )
        for( Int i=0; i<n; ++i )
            Q(i,j) = y(i)*QReal(i,j);

    return info;
}

} // namespace herm_tridiag_eig

template<typename F>
HermitianTridiagEigInfo
HermitianTridiagEig
(       Matrix<Base<F>>& d,
        Matrix<F>& dSub,
        Matrix<Base<F>>& w,
        Matrix<F>& Q, 
  const HermitianTridiagEigCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    return herm_tridiag_eig::Helper( d, dSub, w, Q, ctrl );
}

namespace herm_tridiag_eig {

template<typename Real,typename=EnableIf<IsBlasScalar<Real>>>
HermitianTridiagEigInfo
Helper
( const AbstractDistMatrix<Real>& d,
  const AbstractDistMatrix<Real>& dSub,
        AbstractDistMatrix<Real>& wPre, 
        AbstractDistMatrix<Real>& QPre, 
  const HermitianTridiagEigCtrl<Real>& ctrl )
{
    // NOTE: The computation forces double-precision due to PMRRR limitations
    const Int n = d.Height();
    const Grid& g = d.Grid();
    HermitianTridiagEigInfo info;
    if( ctrl.useQR )
    {
        LogicError("Distributed QR not yet supported");
    }

    ElementalProxyCtrl wCtrl, QCtrl;
    wCtrl.colConstrain = true;
    wCtrl.colAlign = 0;
    QCtrl.rowConstrain = true;
    QCtrl.rowAlign = 0;

    DistMatrixWriteProxy<Real,Real,VR,STAR> wProx( wPre, wCtrl );
    DistMatrixWriteProxy<Real,double,STAR,VR> QProx( QPre, QCtrl );
    auto& w = wProx.Get();
    auto& Q = QProx.Get();

    DistMatrix<double,STAR,STAR> d_STAR_STAR(g), dSub_STAR_STAR(g);
    Copy( d, d_STAR_STAR );
    dSub_STAR_STAR.Resize( n-1, 1, n );
    Copy( dSub, dSub_STAR_STAR );

    Int k;
    if( ctrl.subset.rangeSubset )
    {
        vector<double> dVector(n), dSubVector(n), wVector(n);
        MemCopy( dVector.data(), d_STAR_STAR.Buffer(), n );
        MemCopy( dSubVector.data(), dSub_STAR_STAR.Buffer(), n-1 );
        auto estimate = herm_tridiag_eig::EigEstimate
          ( int(n), dVector.data(), dSubVector.data(),
            wVector.data(), w.ColComm(),
            ctrl.subset.lowerBound, ctrl.subset.upperBound );
        SwapClear( dVector );
        SwapClear( dSubVector );
        k = estimate.numGlobalEigenvalues;
    }
    else if( ctrl.subset.indexSubset )
        k = ( n==0 ? 0 : ctrl.subset.upperIndex-ctrl.subset.lowerIndex+1 );
    else
        k = n;
    Q.Resize( n, k );

    herm_tridiag_eig::Info rangeInfo;
    vector<double> wVector(n);
    if( ctrl.subset.rangeSubset )
        rangeInfo = herm_tridiag_eig::Eig
          ( int(n), d_STAR_STAR.Buffer(), dSub_STAR_STAR.Buffer(), 
            wVector.data(), Q.Buffer(), Q.LDim(), w.ColComm(),
            ctrl.subset.lowerBound, ctrl.subset.upperBound );
    else if( ctrl.subset.indexSubset )
        rangeInfo = herm_tridiag_eig::Eig
          ( int(n), d_STAR_STAR.Buffer(), dSub_STAR_STAR.Buffer(), 
            wVector.data(), Q.Buffer(), Q.LDim(), w.ColComm(),
            int(ctrl.subset.lowerIndex), int(ctrl.subset.upperIndex) );
    else
        rangeInfo = herm_tridiag_eig::Eig
          ( int(n), d_STAR_STAR.Buffer(), dSub_STAR_STAR.Buffer(), 
            wVector.data(), Q.Buffer(), Q.LDim(), w.ColComm() );
    w.Resize( rangeInfo.numGlobalEigenvalues, 1 );
    Q.Resize( n, rangeInfo.numGlobalEigenvalues );
    for( Int iLoc=0; iLoc<w.LocalHeight(); ++iLoc )
        w.SetLocal( iLoc, 0, Real(wVector[iLoc]) );

    auto sortPairs = TaggedSort( w, ctrl.sort );
    for( Int j=0; j<n; ++j )
        w.Set( j, 0, sortPairs[j].value );
    ApplyTaggedSortToEachRow( sortPairs, Q );

    return info;
}

template<typename Real,typename=DisableIf<IsBlasScalar<Real>>,typename=void>
HermitianTridiagEigInfo
Helper
( const AbstractDistMatrix<Real>& d,
  const AbstractDistMatrix<Real>& dSub,
        AbstractDistMatrix<Real>& wPre, 
        AbstractDistMatrix<Real>& QPre, 
  const HermitianTridiagEigCtrl<Real>& ctrl )
{
    HermitianTridiagEigInfo info;
    LogicError
    ("Distributed non-standard HermitianTridiagEig not yet supported");
    return info;
}

template<typename Real,typename=EnableIf<IsBlasScalar<Real>>>
HermitianTridiagEigInfo
Helper
( const AbstractDistMatrix<Real         >& d,
  const AbstractDistMatrix<Complex<Real>>& dSub,
        AbstractDistMatrix<Real         >& wPre, 
        AbstractDistMatrix<Complex<Real>>& QPre, 
  const HermitianTridiagEigCtrl<Real>& ctrl )
{
    // NOTE: The computation forces double-precision due to PMRRR limitations
    const Int n = d.Height();
    const Grid& g = d.Grid();
    typedef Complex<Real> C;
    HermitianTridiagEigInfo info;
    if( ctrl.useQR )
    {
        LogicError("Distributed QR not yet supported");
    }

    DistMatrix<double,STAR,STAR> d_STAR_STAR(g);
    DistMatrix<Complex<double>,STAR,STAR> dSub_STAR_STAR(g);
    Copy( d, d_STAR_STAR );
    dSub_STAR_STAR.Resize( n-1, 1, n );
    Copy( dSub, dSub_STAR_STAR );
    auto& dSubLoc = dSub_STAR_STAR.Matrix();

    DistMatrix<Complex<double>,STAR,STAR> y(n,1,g);
    auto& yLoc = y.Matrix();

    DistMatrix<double,STAR,STAR> dSubReal(g);
    dSubReal.Resize( n-1, 1, n );
    auto dSubRealLoc = dSubReal.Matrix();

    yLoc(0) = 1;
    for( Int j=0; j<n-1; ++j )
    {
        const auto psi = dSubLoc(j);
        const double psiAbs = Abs(psi);
        if( psiAbs == double(0) )
            yLoc(j+1) = 1;
        else
            yLoc(j+1) = ComplexFromPolar(double(1),Arg(psi*yLoc(j)));
        dSubRealLoc(j) = psiAbs;
    }

    ElementalProxyCtrl wCtrl, QCtrl;
    wCtrl.colConstrain = true;
    wCtrl.colAlign = 0;
    QCtrl.rowConstrain = true;
    QCtrl.rowAlign = 0;

    DistMatrixWriteProxy<Real,Real,VR,STAR> wProx( wPre, wCtrl );
    DistMatrixWriteProxy<C,C,STAR,VR> QProx( QPre, QCtrl );
    auto& w = wProx.Get();
    auto& Q = QProx.Get();

    Int k;
    if( ctrl.subset.rangeSubset )
    {
        vector<double> dVector(n), dSubVector(n), wVector(n);
        MemCopy( dVector.data(), d_STAR_STAR.Buffer(), n );
        MemCopy( dSubVector.data(), dSubReal.Buffer(), n-1 );
        auto estimate = herm_tridiag_eig::EigEstimate
          ( int(n), dVector.data(), dSubVector.data(),
            wVector.data(), w.ColComm(),
            ctrl.subset.lowerBound, ctrl.subset.upperBound );
        SwapClear( dVector );
        SwapClear( dSubVector );
        k = estimate.numGlobalEigenvalues;
    }
    else if( ctrl.subset.indexSubset )
        k = ( n==0 ? 0 : ctrl.subset.upperIndex-ctrl.subset.lowerIndex+1 );
    else
        k = n;
    DistMatrix<double,STAR,VR> QReal(g);
    QReal.Resize( n, k );

    herm_tridiag_eig::Info rangeInfo;
    vector<double> wVector(n);
    if( ctrl.subset.rangeSubset )
        rangeInfo = herm_tridiag_eig::Eig
          ( int(n), d_STAR_STAR.Buffer(), dSubReal.Buffer(), 
            wVector.data(), QReal.Buffer(), QReal.LDim(), w.ColComm(),
            ctrl.subset.lowerBound, ctrl.subset.upperBound );
    else if( ctrl.subset.indexSubset )
        rangeInfo = herm_tridiag_eig::Eig
          ( int(n), d_STAR_STAR.Buffer(), dSubReal.Buffer(), 
            wVector.data(), QReal.Buffer(), QReal.LDim(), w.ColComm(),
            int(ctrl.subset.lowerIndex), int(ctrl.subset.upperIndex) );
    else
        rangeInfo = herm_tridiag_eig::Eig
          ( int(n), d_STAR_STAR.Buffer(), dSubReal.Buffer(), 
            wVector.data(), QReal.Buffer(), QReal.LDim(), w.ColComm() );

    w.Resize( rangeInfo.numGlobalEigenvalues, 1 );
    auto& wLoc = w.Matrix();
    for( Int iLoc=0; iLoc<w.LocalHeight(); ++iLoc )
        wLoc(iLoc) = wVector[iLoc];

    QReal.Resize( n, rangeInfo.numGlobalEigenvalues );

    auto sortPairs = TaggedSort( w, ctrl.sort );
    for( Int j=0; j<n; ++j )
        w.Set( j, 0, sortPairs[j].value );
    ApplyTaggedSortToEachRow( sortPairs, QReal );

    Q.Resize( n, rangeInfo.numGlobalEigenvalues );
    auto& QLoc = Q.Matrix();
    auto& QRealLoc = QReal.Matrix();
    for( Int jLoc=0; jLoc<Q.LocalWidth(); ++jLoc )
        for( Int i=0; i<n; ++i )
            QLoc(i,jLoc) = C(yLoc(i)*QRealLoc(i,jLoc));

    return info;
}

template<typename Real,typename=DisableIf<IsBlasScalar<Real>>,typename=void>
HermitianTridiagEigInfo
Helper
( const AbstractDistMatrix<Real         >& d,
  const AbstractDistMatrix<Complex<Real>>& dSub,
        AbstractDistMatrix<Real         >& wPre, 
        AbstractDistMatrix<Complex<Real>>& QPre, 
  const HermitianTridiagEigCtrl<Real>& ctrl )
{
    HermitianTridiagEigInfo info;
    LogicError
    ("Distributed non-standard HermitianTridiagEig not yet supported");
    return info;
}

} // namespace herm_tridiag_eig

template<typename F>
HermitianTridiagEigInfo
HermitianTridiagEig
( const AbstractDistMatrix<Base<F>>& d,
  const AbstractDistMatrix<F>& dSub,
        AbstractDistMatrix<Base<F>>& w,
        AbstractDistMatrix<F>& Q, 
  const HermitianTridiagEigCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    return herm_tridiag_eig::Helper( d, dSub, w, Q, ctrl );
}

namespace herm_tridiag_eig {

template<typename Real,typename=EnableIf<IsBlasScalar<Real>>>
Int MRRREstimateHelper
( const AbstractDistMatrix<Real>& d,
  const AbstractDistMatrix<Real>& dSub,
        mpi::Comm wColComm,
        Real vl,
        Real vu )
{
    DEBUG_CSE
    const Int n = d.Height();
    DistMatrix<double,STAR,STAR> d_STAR_STAR( d.Grid() );
    DistMatrix<double,STAR,STAR> dSub_STAR_STAR( d.Grid() );
    Copy( d, d_STAR_STAR );
    dSub_STAR_STAR.Resize( n-1, 1, n );
    Copy( dSub, dSub_STAR_STAR );
    vector<double> dVector(n), dSubVector(n), wVector(n);
    MemCopy( dVector.data(), d_STAR_STAR.Buffer(), n );
    MemCopy( dSubVector.data(), dSub_STAR_STAR.Buffer(), n-1 );
    auto estimate = herm_tridiag_eig::EigEstimate
    ( int(n), dVector.data(), dSubVector.data(), wVector.data(), wColComm,
      vl, vu );
    return estimate.numGlobalEigenvalues;
}

template<typename Real,typename=DisableIf<IsBlasScalar<Real>>,typename=void>
Int MRRREstimateHelper
( const AbstractDistMatrix<Real>& d,
  const AbstractDistMatrix<Real>& dSub,
        mpi::Comm wColComm,
        Real vl,
        Real vu )
{
    DEBUG_CSE
    LogicError
    ("HermitianTridiagEigEstimate not yet supported for nonstandard datatypes");
    return 0;
}

// Q is assumed to be sufficiently large and properly aligned
template<typename Real,typename=EnableIf<IsBlasScalar<Real>>>
HermitianTridiagEigInfo
MRRRPostEstimateHelper
( const AbstractDistMatrix<Real>& d,
  const AbstractDistMatrix<Real>& dSub,
        AbstractDistMatrix<Real>& wPre,
        AbstractDistMatrix<Real>& QPre,
        SortType sort,
        Real vl,
        Real vu )
{
    DEBUG_CSE
    HermitianTridiagEigInfo info;

    ElementalProxyCtrl wCtrl, QCtrl;
    wCtrl.colConstrain = true;
    wCtrl.colAlign = 0;
    QCtrl.rowConstrain = true;
    QCtrl.rowAlign = 0;

    DistMatrixWriteProxy<Real,Real,VR,STAR> wProx( wPre, wCtrl );
    DistMatrixWriteProxy<Real,double,STAR,VR> QProx( QPre, QCtrl );
    auto& w = wProx.Get();
    auto& Q = QProx.Get();

    const Int n = d.Height();
    DistMatrix<double,STAR,STAR> d_STAR_STAR( d.Grid() ),
                                 dSub_STAR_STAR( d.Grid() );
    Copy( d, d_STAR_STAR );
    dSub_STAR_STAR.Resize( n-1, 1, n );
    Copy( dSub, dSub_STAR_STAR );

    vector<double> wVector(n);
    auto rangeInfo = herm_tridiag_eig::Eig
    ( int(n), d_STAR_STAR.Buffer(), dSub_STAR_STAR.Buffer(), wVector.data(),
      Q.Buffer(), Q.LDim(), w.ColComm(), vl, vu );
    const Int k = rangeInfo.numGlobalEigenvalues;

    w.Resize( k, 1 );
    for( Int iLoc=0; iLoc<w.LocalHeight(); ++iLoc )
        w.SetLocal( iLoc, 0, Real(wVector[iLoc]) );

    // Shrink Q
    Q.Resize( n, k );

    auto sortPairs = TaggedSort( w, sort );
    for( Int j=0; j<n; ++j )
        w.Set( j, 0, sortPairs[j].value );
    ApplyTaggedSortToEachRow( sortPairs, Q );

    return info;
}

template<typename Real,typename=DisableIf<IsBlasScalar<Real>>,typename=void>
HermitianTridiagEigInfo
MRRRPostEstimateHelper
( const AbstractDistMatrix<Real>& d,
  const AbstractDistMatrix<Real>& dSub,
        AbstractDistMatrix<Real>& wPre,
        AbstractDistMatrix<Real>& QPre,
        SortType sort,
        Real vl,
        Real vu )
{
    DEBUG_CSE
    LogicError
    ("HermitianTridiagEigEstimate not yet supported for nonstandard datatypes");
    HermitianTridiagEigInfo info;
    return info;
}

template<typename Real>
Int MRRREstimate
( const AbstractDistMatrix<Real>& d,
  const AbstractDistMatrix<Real>& dSub,
        mpi::Comm wColComm,
        Real vl,
        Real vu )
{
    DEBUG_CSE
    return MRRREstimateHelper( d, dSub, wColComm, vl, vu );
}

// Q is assumed to be sufficiently large and properly aligned
template<typename Real>
HermitianTridiagEigInfo
MRRRPostEstimate
( const AbstractDistMatrix<Real>& d,
  const AbstractDistMatrix<Real>& dSub,
        AbstractDistMatrix<Real>& w,
        AbstractDistMatrix<Real>& Q,
        SortType sort,
        Real vl,
        Real vu )
{
    DEBUG_CSE
    return MRRRPostEstimateHelper( d, dSub, w, Q, sort, vl, vu );
}

} // namespace herm_tridiag_eig

#define PROTO(F) \
  template void herm_eig::SortAndFilter \
  ( Matrix<Base<F>>& w, \
    Matrix<F>& Q, \
    const HermitianTridiagEigCtrl<Base<F>>& ctrl ); \
  template HermitianTridiagEigInfo HermitianTridiagEig \
  (       Matrix<Base<F>>& d, \
          Matrix<F>& dSub, \
          Matrix<Base<F>>& w, \
    const HermitianTridiagEigCtrl<Base<F>>& ctrl ); \
  template HermitianTridiagEigInfo HermitianTridiagEig \
  (       Matrix<Base<F>>& d, \
          Matrix<F>& dSub, \
          Matrix<Base<F>>& w, \
          Matrix<F>& Q, \
    const HermitianTridiagEigCtrl<Base<F>>& ctrl ); \
  template HermitianTridiagEigInfo HermitianTridiagEig \
  ( const AbstractDistMatrix<Base<F>>& d, \
    const AbstractDistMatrix<F>& dSub, \
          AbstractDistMatrix<Base<F>>& w, \
    const HermitianTridiagEigCtrl<Base<F>>& ctrl ); \
  template HermitianTridiagEigInfo HermitianTridiagEig \
  ( const AbstractDistMatrix<Base<F>>& d, \
    const AbstractDistMatrix<F>& dSub, \
          AbstractDistMatrix<Base<F>>& w, \
          AbstractDistMatrix<F>& Q, \
    const HermitianTridiagEigCtrl<Base<F>>& ctrl );

#define PROTO_REAL(Real) \
  PROTO(Real) \
  template void herm_eig::SortAndFilter \
  ( Matrix<Real>& w, \
    const HermitianTridiagEigCtrl<Real>& ctrl ); \
  template Int herm_tridiag_eig::MRRREstimate \
  ( const AbstractDistMatrix<Real>& d, \
    const AbstractDistMatrix<Real>& dSub, \
          mpi::Comm wColComm, Real vl, Real vu ); \
  template HermitianTridiagEigInfo herm_tridiag_eig::MRRRPostEstimate \
  ( const AbstractDistMatrix<Real>& d, \
    const AbstractDistMatrix<Real>& dSub, \
          AbstractDistMatrix<Real>& w, \
          AbstractDistMatrix<Real>& Q, \
          SortType sort, \
          Real vl, \
          Real vu );

#define EL_NO_INT_PROTO
#define EL_ENABLE_QUAD
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
