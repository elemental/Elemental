/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

#include "./HermitianTridiagEig/QR.hpp"
#include "./HermitianTridiagEig/DivideAndConquer.hpp"

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

template<typename Real>
void SortAndFilter
( AbstractDistMatrix<Real>& wPre, const HermitianTridiagEigCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE

    DistMatrixReadWriteProxy<Real,Real,STAR,STAR> wProx( wPre );
    auto& w = wProx.Get();

    const Int n = w.Height();
    const Grid& g = w.Grid();

    if( ctrl.subset.indexSubset )
    {
        Sort( w, ctrl.sort );
        DistMatrix<Real,STAR,STAR> wCopy( w );
        w = wCopy(IR(ctrl.subset.lowerIndex,ctrl.subset.upperIndex+1),ALL);
    }
    else if( ctrl.subset.rangeSubset )
    {
        Int numValid = 0;

        for( Int j=0; j<n; ++j )
            if( w.GetLocal(j,0) > ctrl.subset.lowerBound &&
                w.GetLocal(j,0) <= ctrl.subset.upperBound )
                ++numValid;

        DistMatrix<Real,STAR,STAR> wFilter(numValid,1,g);
        numValid = 0;
        for( Int j=0; j<n; ++j )
           if( w.GetLocal(j,0) > ctrl.subset.lowerBound &&
               w.GetLocal(j,0) <= ctrl.subset.upperBound )
               wFilter.SetLocal( numValid++, 0, w.GetLocal(j,0) );

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

template<typename F>
void SortAndFilter
( AbstractDistMatrix<Base<F>>& wPre,
  AbstractDistMatrix<F>& QPre,
  const HermitianTridiagEigCtrl<Base<F>>& ctrl )
{
    EL_DEBUG_CSE
    typedef Base<F> Real;

    DistMatrixReadWriteProxy<Real,Real,STAR,STAR> wProx( wPre );
    DistMatrixReadWriteProxy<F,F,VC,STAR> QProx( QPre );
    auto& w = wProx.Get();
    auto& Q = QProx.Get();

    const Int n = w.Height();
    const Grid& g = w.Grid();

    if( ctrl.subset.indexSubset )
    {
        auto sortPairs = TaggedSort( w, ctrl.sort );
        for( Int j=0; j<n; ++j )
            w.SetLocal( j, 0, sortPairs[j].value );
        ApplyTaggedSortToEachRow( sortPairs, Q );

        DistMatrix<Real,STAR,STAR> wCopy( w );
        DistMatrix<F,VC,STAR> QCopy( Q );
        w = wCopy(IR(ctrl.subset.lowerIndex,ctrl.subset.upperIndex+1),ALL);
        Q = QCopy(ALL,IR(ctrl.subset.lowerIndex,ctrl.subset.upperIndex+1));
    }
    else if( ctrl.subset.rangeSubset )
    {
        Int numValid = 0;

        for( Int j=0; j<n; ++j )
            if( w.GetLocal(j,0) > ctrl.subset.lowerBound &&
                w.GetLocal(j,0) <= ctrl.subset.upperBound )
                ++numValid;

        DistMatrix<Real,STAR,STAR> wFilter(numValid,1,g);
        DistMatrix<F,VC,STAR> QFilter(n,numValid,g);
        numValid = 0;
        for( Int j=0; j<n; ++j )
        {
           if( w.GetLocal(j,0) > ctrl.subset.lowerBound &&
               w.GetLocal(j,0) <= ctrl.subset.upperBound )
           {
               wFilter.SetLocal( numValid, 0, w.GetLocal(j,0) );
               auto qFilterCol = QFilter(ALL,IR(numValid));
               qFilterCol = Q(ALL,IR(j));
               ++numValid;
           }
        }

        w = wFilter;
        Q = QFilter;

        auto sortPairs = TaggedSort( w, ctrl.sort );
        for( Int j=0; j<numValid; ++j )
            w.SetLocal( j, 0, sortPairs[j].value );
        ApplyTaggedSortToEachRow( sortPairs, Q );
    }
    else
    {
        auto sortPairs = TaggedSort( w, ctrl.sort );
        for( Int j=0; j<n; ++j )
            w.SetLocal( j, 0, sortPairs[j].value );
        ApplyTaggedSortToEachRow( sortPairs, Q );
    }
}

} // namespace herm_eig

namespace herm_tridiag_eig {

template<typename Real>
void RemovePhase
( const Matrix<Complex<Real>>& dSub,
        Matrix<Real>& dSubReal )
{
    EL_DEBUG_CSE
    const Int n = dSub.Height() + 1;
    dSubReal.Resize( n-1, 1 );
    Complex<Real> phaseLast(1);
    for( Int j=0; j<n-1; ++j )
    {
        const auto psi = dSub(j);
        const auto psiAbs = Abs(psi);
        if( psiAbs == Real(0) )
            phaseLast = 1;
        else
            phaseLast = ComplexFromPolar(Real(1),Arg(psi*phaseLast));
        dSubReal(j) = psiAbs;
    }
}

template<typename Real>
void RemovePhase
( const DistMatrix<Complex<Real>,STAR,STAR>& dSub,
        DistMatrix<Real,STAR,STAR>& dSubReal )
{
    EL_DEBUG_CSE
    const Int n = dSub.Height() + 1;
    const Grid& g = dSub.Grid();
    dSubReal.SetGrid( g );
    dSubReal.Resize( n-1, 1 );
    RemovePhase( dSub.LockedMatrix(), dSubReal.Matrix() );
}

template<typename Real>
void RemovePhase
( const Matrix<Complex<Real>>& dSub,
        Matrix<Real>& dSubReal,
        Matrix<Complex<Real>>& phase )
{
    EL_DEBUG_CSE
    const Int n = dSub.Height() + 1;
    dSubReal.Resize( n-1, 1 );
    phase.Resize( n, 1 );
    phase(0) = 1;
    for( Int j=0; j<n-1; ++j )
    {
        const auto psi = dSub(j);
        const auto psiAbs = Abs(psi);
        if( psiAbs == Real(0) )
            phase(j+1) = 1;
        else
            phase(j+1) = ComplexFromPolar(Real(1),Arg(psi*phase(j)));
        dSubReal(j) = psiAbs;
    }
}

template<typename Real>
void RemovePhase
( const DistMatrix<Complex<Real>,STAR,STAR>& dSub,
        DistMatrix<Real,STAR,STAR>& dSubReal,
        DistMatrix<Complex<Real>,STAR,STAR>& phase )
{
    EL_DEBUG_CSE
    const Int n = dSub.Height() + 1;
    const Grid& g = dSub.Grid();
    dSubReal.SetGrid( g );
    dSubReal.Resize( n-1, 1 );
    phase.SetGrid( g );
    phase.Resize( n, 1 );
    RemovePhase( dSub.LockedMatrix(), dSubReal.Matrix(), phase.Matrix() );
}

template<typename Real>
HermitianTridiagEigInfo
QRHelper
( const Matrix<Real>& d,
        Matrix<Real>& dSub,
        Matrix<Real>& w,
  const HermitianTridiagEigCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    HermitianTridiagEigInfo info;
    w = d;    
    info.qrInfo = QRAlg( w, dSub, ctrl );
    herm_eig::SortAndFilter( w, ctrl );
    return info;
}

template<typename Real>
HermitianTridiagEigInfo
DCHelper
(       Matrix<Real>& d,
        Matrix<Real>& dSub,
        Matrix<Real>& w,
  const HermitianTridiagEigCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    HermitianTridiagEigInfo info;
    auto ctrlMod( ctrl );
    ctrlMod.wantEigVecs = false;
    Matrix<Real> Q;
    info.dcInfo = DivideAndConquer( d, dSub, w, Q, ctrlMod );
    herm_eig::SortAndFilter( w, ctrl );
    return info;
}

template<typename Real,typename=EnableIf<IsBlasScalar<Real>>>
HermitianTridiagEigInfo
LAPACKHelper
(       Matrix<Real>& d,
        Matrix<Real>& dSub,
        Matrix<Real>& w,
  const HermitianTridiagEigCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Int n = d.Height();
    HermitianTridiagEigInfo info;

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

    return info;
}

template<typename Real,typename=EnableIf<IsBlasScalar<Real>>>
HermitianTridiagEigInfo
Helper
( const Matrix<Real>& d,
  const Matrix<Real>& dSub,
        Matrix<Real>& w,
  const HermitianTridiagEigCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    if( ctrl.alg == HERM_TRIDIAG_EIG_QR )
    {
        // Only dSub needs to be modifiable
        auto dSubMod( dSub );
        return QRHelper( d, dSubMod, w, ctrl );
    }
    else if( ctrl.alg == HERM_TRIDIAG_EIG_DC )
    {
        auto dMod( d );
        auto dSubMod( dSub );
        return DCHelper( dMod, dSubMod, w, ctrl );
    }
    // Both d and dSub need to be modifiable
    auto dMod( d );
    auto dSubMod( dSub );
    return LAPACKHelper( dMod, dSubMod, w, ctrl );
}

template<typename Real,typename=DisableIf<IsBlasScalar<Real>>,typename=void>
HermitianTridiagEigInfo
Helper
( const Matrix<Real>& d,
  const Matrix<Real>& dSub,
        Matrix<Real>& w,
  const HermitianTridiagEigCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    if( ctrl.alg == HERM_TRIDIAG_EIG_QR )
    {
        // Only dSub needs to be modifiable
        auto dSubMod( dSub );
        return QRHelper( d, dSubMod, w, ctrl );
    }
    else
    {
        auto dMod( d );
        auto dSubMod( dSub );
        return DCHelper( dMod, dSubMod, w, ctrl );
    }
}

template<typename Real>
HermitianTridiagEigInfo
Helper
( const Matrix<Real>& d,
  const Matrix<Complex<Real>>& dSub,
        Matrix<Real>& w,
  const HermitianTridiagEigCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    // TODO(poulson): Avoid making another copy of dSubReal later
    Matrix<Real> dSubReal;
    RemovePhase( dSub, dSubReal );
    return HermitianTridiagEig( d, dSubReal, w, ctrl );
}

} // namespace herm_tridiag_eig

template<typename F>
HermitianTridiagEigInfo
HermitianTridiagEig
( const Matrix<Base<F>>& d,
  const Matrix<F>& dSub,
        Matrix<Base<F>>& w,
  const HermitianTridiagEigCtrl<Base<F>>& ctrl )
{
    EL_DEBUG_CSE
    return herm_tridiag_eig::Helper( d, dSub, w, ctrl );
}

namespace herm_tridiag_eig {

template<typename Real>
HermitianTridiagEigInfo
QRHelper
( const AbstractDistMatrix<Real>& d,
  const AbstractDistMatrix<Real>& dSub,
        AbstractDistMatrix<Real>& w, 
  const HermitianTridiagEigCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    HermitianTridiagEigInfo info;
    DistMatrix<Real,STAR,STAR> d_STAR_STAR(d), dSub_STAR_STAR(dSub);
    info.qrInfo = QRAlg( d_STAR_STAR, dSub_STAR_STAR, ctrl );
    herm_eig::SortAndFilter( d_STAR_STAR, ctrl );
    Copy( d_STAR_STAR, w );
    return info;
}

template<typename Real>
HermitianTridiagEigInfo
QRHelper
( const AbstractDistMatrix<Real         >& d,
  const AbstractDistMatrix<Complex<Real>>& dSub,
        AbstractDistMatrix<Real         >& w,
  const HermitianTridiagEigCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    HermitianTridiagEigInfo info;
    const Grid& g = d.Grid();

    DistMatrix<Real,STAR,STAR> d_STAR_STAR( d );
    DistMatrix<Complex<Real>,STAR,STAR> dSub_STAR_STAR( dSub );

    DistMatrix<Real,STAR,STAR> dSubReal(g);
    RemovePhase( dSub_STAR_STAR, dSubReal );

    info.qrInfo = QRAlg( d_STAR_STAR, dSubReal, ctrl );
    herm_eig::SortAndFilter( d_STAR_STAR, ctrl );
    Copy( d_STAR_STAR, w );

    return info;
}

template<typename Real>
HermitianTridiagEigInfo
DCHelper
( const AbstractDistMatrix<Real>& d,
  const AbstractDistMatrix<Real>& dSub,
        AbstractDistMatrix<Real>& wPre,
  const HermitianTridiagEigCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    HermitianTridiagEigInfo info;
    DistMatrix<Real,STAR,STAR> d_STAR_STAR(d), dSub_STAR_STAR(dSub);

    DistMatrixWriteProxy<Real,Real,STAR,STAR> wProx( wPre );
    auto& w = wProx.Get();

    auto ctrlMod( ctrl );
    ctrlMod.wantEigVecs = false;
    DistMatrix<Real> Q(w.Grid());
    info.dcInfo =
      DivideAndConquer
      ( d_STAR_STAR.Matrix(), dSub_STAR_STAR.Matrix(), w, Q, ctrlMod );
    herm_eig::SortAndFilter( w, ctrl );
    return info;
}

template<typename Real>
HermitianTridiagEigInfo
DCHelper
( const AbstractDistMatrix<Real         >& d,
  const AbstractDistMatrix<Complex<Real>>& dSub,
        AbstractDistMatrix<Real         >& wPre,
  const HermitianTridiagEigCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    HermitianTridiagEigInfo info;
    const Grid& g = d.Grid();

    DistMatrix<Real,STAR,STAR> d_STAR_STAR( d );
    DistMatrix<Complex<Real>,STAR,STAR> dSub_STAR_STAR( dSub );

    DistMatrix<Real,STAR,STAR> dSubReal(g);
    RemovePhase( dSub_STAR_STAR, dSubReal );

    DistMatrixWriteProxy<Real,Real,STAR,STAR> wProx( wPre );
    auto& w = wProx.Get();

    auto ctrlMod( ctrl );
    ctrlMod.wantEigVecs = false;
    DistMatrix<Real> Q(w.Grid());
    info.dcInfo =
      DivideAndConquer( d_STAR_STAR.Matrix(), dSubReal.Matrix(), w, Q, ctrl );
    herm_eig::SortAndFilter( w, ctrl );

    return info;
}

template<typename Real,typename=EnableIf<IsBlasScalar<Real>>>
HermitianTridiagEigInfo
MRRRHelper
( const AbstractDistMatrix<Real>& d,
  const AbstractDistMatrix<Real>& dSub,
        AbstractDistMatrix<Real>& wPre,
  const HermitianTridiagEigCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
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

template<typename Real,typename=EnableIf<IsBlasScalar<Real>>>
HermitianTridiagEigInfo
Helper
( const AbstractDistMatrix<Real>& d,
  const AbstractDistMatrix<Real>& dSub,
        AbstractDistMatrix<Real>& wPre,
  const HermitianTridiagEigCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    if( ctrl.alg == HERM_TRIDIAG_EIG_QR )
    {
        return QRHelper( d, dSub, wPre, ctrl );
    }
    else if( ctrl.alg == HERM_TRIDIAG_EIG_DC )
    {
        return DCHelper( d, dSub, wPre, ctrl );
    }
    else
    {
        return MRRRHelper( d, dSub, wPre, ctrl );
    }
}

template<typename Real,typename=DisableIf<IsBlasScalar<Real>>,typename=void>
HermitianTridiagEigInfo
Helper
( const AbstractDistMatrix<Real>& d,
  const AbstractDistMatrix<Real>& dSub,
        AbstractDistMatrix<Real>& w, 
  const HermitianTridiagEigCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    if( ctrl.alg == HERM_TRIDIAG_EIG_QR )
    {
        return QRHelper( d, dSub, w, ctrl );
    }
    else
    {
        return DCHelper( d, dSub, w, ctrl );
    }
}

template<typename Real,typename=EnableIf<IsBlasScalar<Real>>>
HermitianTridiagEigInfo
MRRRHelper
( const AbstractDistMatrix<Real         >& d,
  const AbstractDistMatrix<Complex<Real>>& dSub,
        AbstractDistMatrix<Real         >& wPre, 
  const HermitianTridiagEigCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
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

    DistMatrix<double,STAR,STAR> dSubReal(g);
    RemovePhase( dSub_STAR_STAR, dSubReal );

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

template<typename Real,typename=EnableIf<IsBlasScalar<Real>>>
HermitianTridiagEigInfo
Helper
( const AbstractDistMatrix<Real         >& d,
  const AbstractDistMatrix<Complex<Real>>& dSub,
        AbstractDistMatrix<Real         >& wPre, 
  const HermitianTridiagEigCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    if( ctrl.alg == HERM_TRIDIAG_EIG_QR )
    {
        return QRHelper( d, dSub, wPre, ctrl );
    }
    else if( ctrl.alg == HERM_TRIDIAG_EIG_DC )
    {
        return DCHelper( d, dSub, wPre, ctrl );
    }
    else
    {
        return MRRRHelper( d, dSub, wPre, ctrl );
    }
}

template<typename Real,typename=DisableIf<IsBlasScalar<Real>>,typename=void>
HermitianTridiagEigInfo
Helper
( const AbstractDistMatrix<Real         >& d,
  const AbstractDistMatrix<Complex<Real>>& dSub,
        AbstractDistMatrix<Real         >& w, 
  const HermitianTridiagEigCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    if( ctrl.alg == HERM_TRIDIAG_EIG_QR )
    {
        return QRHelper( d, dSub, w, ctrl );
    }
    else
    {
        return DCHelper( d, dSub, w, ctrl );
    }
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
    EL_DEBUG_CSE
    return herm_tridiag_eig::Helper( d, dSub, w, ctrl );
}

// Return eigenpairs
// =================

namespace herm_tridiag_eig {

template<typename Real>
HermitianTridiagEigInfo
QRHelper
( const Matrix<Real>& d,
        Matrix<Real>& dSub,
        Matrix<Real>& w,
        Matrix<Real>& Q,
  const HermitianTridiagEigCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    HermitianTridiagEigInfo info;
    w = d;    
    info.qrInfo = QRAlg( w, dSub, Q, ctrl );
    herm_eig::SortAndFilter( w, Q, ctrl );
    return info;
}

template<typename Real>
HermitianTridiagEigInfo
DCHelper
(       Matrix<Real>& d,
        Matrix<Real>& dSub,
        Matrix<Real>& w,
        Matrix<Real>& Q,
  const HermitianTridiagEigCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    HermitianTridiagEigInfo info;
    if( ctrl.accumulateEigVecs )
    {
        LogicError("Accumulating eigenvectors not yet supported for D&C");
    }
    else
    {
        info.dcInfo = DivideAndConquer( d, dSub, w, Q, ctrl );
        herm_eig::SortAndFilter( w, Q, ctrl );
    }
    return info;
}

template<typename Real,typename=EnableIf<IsBlasScalar<Real>>>
HermitianTridiagEigInfo
LAPACKHelper
(       Matrix<Real>& d,
        Matrix<Real>& dSub,
        Matrix<Real>& w,
        Matrix<Real>& Q,
  const HermitianTridiagEigCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    if( ctrl.accumulateEigVecs )
        LogicError("Accumulation not yet supported for LAPACK wrapper");
    const Int n = d.Height();
    HermitianTridiagEigInfo info;
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

    return info;
}

template<typename Real,typename=EnableIf<IsBlasScalar<Real>>>
HermitianTridiagEigInfo
Helper
( const Matrix<Real>& d,
  const Matrix<Real>& dSub,
        Matrix<Real>& w,
        Matrix<Real>& Q,
  const HermitianTridiagEigCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    if( ctrl.alg == HERM_TRIDIAG_EIG_QR )
    {
        // Only dSub needs to be modified
        auto dSubMod( dSub );
        return QRHelper( d, dSubMod, w, Q, ctrl );
    }
    else if( ctrl.alg == HERM_TRIDIAG_EIG_DC )
    {
        auto dMod( d );
        auto dSubMod( dSub );
        return DCHelper( dMod, dSubMod, w, Q, ctrl );
    }
    else
    {
        // Both d and dSub need to be modified
        auto dMod( d );
        auto dSubMod( dSub );
        return LAPACKHelper( dMod, dSubMod, w, Q, ctrl );
    }
}

template<typename Real,typename=DisableIf<IsBlasScalar<Real>>,typename=void>
HermitianTridiagEigInfo
Helper
( const Matrix<Real>& d,
  const Matrix<Real>& dSub,
        Matrix<Real>& w,
        Matrix<Real>& Q,
  const HermitianTridiagEigCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    if( ctrl.alg == HERM_TRIDIAG_EIG_QR )
    {
        // Only dSub needs to be modified
        auto dSubMod( dSub );
        return QRHelper( d, dSubMod, w, Q, ctrl );
    }
    else
    {
        auto dMod( d );
        auto dSubMod( dSub );
        return DCHelper( dMod, dSubMod, w, Q, ctrl );
    }
}

// (Y^H T Y) QHat = QHat Lambda
// T (Y QHat) = (Y QHat) Lambda
template<typename Real>
HermitianTridiagEigInfo
Helper
( const Matrix<Real>& d,
  const Matrix<Complex<Real>>& dSub,
        Matrix<Real>& w, 
        Matrix<Complex<Real>>& Q,
  const HermitianTridiagEigCtrl<Real>& ctrl )
{
    HermitianTridiagEigInfo info;

    // TODO(poulson): Avoid another copy of dSubReal down the road
    Matrix<Real> dSubReal;
    Matrix<Complex<Real>> phase;
    RemovePhase( dSub, dSubReal, phase );

    Matrix<Real> QReal;
    HermitianTridiagEig( d, dSubReal, w, QReal, ctrl );

    Copy( QReal, Q );
    DiagonalScale( LEFT, NORMAL, phase, Q );

    return info;
}

} // namespace herm_tridiag_eig

template<typename F>
HermitianTridiagEigInfo
HermitianTridiagEig
( const Matrix<Base<F>>& d,
  const Matrix<F>& dSub,
        Matrix<Base<F>>& w,
        Matrix<F>& Q, 
  const HermitianTridiagEigCtrl<Base<F>>& ctrl )
{
    EL_DEBUG_CSE
    return herm_tridiag_eig::Helper( d, dSub, w, Q, ctrl );
}

namespace herm_tridiag_eig {

template<typename Real>
HermitianTridiagEigInfo
QRHelper
( const AbstractDistMatrix<Real>& d,
  const AbstractDistMatrix<Real>& dSub,
        AbstractDistMatrix<Real>& w, 
        AbstractDistMatrix<Real>& QPre, 
  const HermitianTridiagEigCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    HermitianTridiagEigInfo info;
    DistMatrix<Real,STAR,STAR> d_STAR_STAR(d), dSub_STAR_STAR(dSub);

    if( ctrl.accumulateEigVecs )
    {
        DistMatrixReadWriteProxy<Real,Real,VC,STAR> QProx( QPre );
        auto& Q = QProx.Get();
        info.qrInfo = QRAlg( d_STAR_STAR, dSub_STAR_STAR, Q, ctrl );
        herm_eig::SortAndFilter( d_STAR_STAR, Q, ctrl );
    }
    else
    {
        DistMatrixWriteProxy<Real,Real,VC,STAR> QProx( QPre );
        auto& Q = QProx.Get();

        info.qrInfo = QRAlg( d_STAR_STAR, dSub_STAR_STAR, Q, ctrl );
        herm_eig::SortAndFilter( d_STAR_STAR, Q, ctrl );
    }
    Copy( d_STAR_STAR, w );

    return info;
}

template<typename Real>
HermitianTridiagEigInfo
QRHelper
( const AbstractDistMatrix<Real         >& d,
  const AbstractDistMatrix<Complex<Real>>& dSub,
        AbstractDistMatrix<Real         >& w, 
        AbstractDistMatrix<Complex<Real>>& QPre, 
  const HermitianTridiagEigCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    typedef Complex<Real> F;
    const Grid& g = d.Grid();
    HermitianTridiagEigInfo info;

    DistMatrix<Real,STAR,STAR> d_STAR_STAR( d );
    DistMatrix<F,STAR,STAR> dSub_STAR_STAR( dSub );

    DistMatrix<Real,STAR,STAR> dSubReal(g);
    DistMatrix<F,STAR,STAR> phase(g);
    RemovePhase( dSub_STAR_STAR, dSubReal, phase );

    if( ctrl.accumulateEigVecs )
    {
        DistMatrixReadWriteProxy<F,F,VC,STAR> QProx(QPre);
        auto& Q = QProx.Get();

        info.qrInfo = QRAlg( d_STAR_STAR, dSubReal, Q, ctrl );
        herm_eig::SortAndFilter( d_STAR_STAR, Q, ctrl );

        DiagonalScale( LEFT, NORMAL, phase, Q );
    }
    else
    {
        const Int n = d.Height();
        DistMatrix<Real,VC,STAR> QReal(g);
        Identity( QReal, n, n );

        info.qrInfo = QRAlg( d_STAR_STAR, dSubReal, QReal, ctrl );
        herm_eig::SortAndFilter( d_STAR_STAR, QReal, ctrl );

        Copy( QReal, QPre );
        DiagonalScale( LEFT, NORMAL, phase, QPre );
    }
    Copy( d_STAR_STAR, w );

    return info;
}

template<typename Real>
HermitianTridiagEigInfo
DCHelper
( const AbstractDistMatrix<Real>& d,
  const AbstractDistMatrix<Real>& dSub,
        AbstractDistMatrix<Real>& wPre, 
        AbstractDistMatrix<Real>& QPre, 
  const HermitianTridiagEigCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    HermitianTridiagEigInfo info;
    DistMatrix<Real,STAR,STAR> d_STAR_STAR(d), dSub_STAR_STAR(dSub);

    if( ctrl.accumulateEigVecs )
    {
        LogicError("Accumulation not yet supported for D&C");
    }
    else
    {
        DistMatrixWriteProxy<Real,Real,STAR,STAR> wProx( wPre );
        auto& w = wProx.Get();
        DistMatrixWriteProxy<Real,Real,MC,MR> QProx( QPre );
        auto& Q = QProx.Get();
        info.dcInfo =
          DivideAndConquer
          ( d_STAR_STAR.Matrix(), dSub_STAR_STAR.Matrix(), w, Q, ctrl );
        herm_eig::SortAndFilter( w, Q, ctrl );
    }

    return info;
}

template<typename Real>
HermitianTridiagEigInfo
DCHelper
( const AbstractDistMatrix<Real         >& d,
  const AbstractDistMatrix<Complex<Real>>& dSub,
        AbstractDistMatrix<Real         >& wPre,
        AbstractDistMatrix<Complex<Real>>& Q,
  const HermitianTridiagEigCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    typedef Complex<Real> F;
    const Grid& g = d.Grid();
    HermitianTridiagEigInfo info;

    DistMatrix<Real,STAR,STAR> d_STAR_STAR( d );
    DistMatrix<F,STAR,STAR> dSub_STAR_STAR( dSub );

    DistMatrix<Real,STAR,STAR> dSubReal(g);
    DistMatrix<F,STAR,STAR> phase(g);
    RemovePhase( dSub_STAR_STAR, dSubReal, phase );

    if( ctrl.accumulateEigVecs )
    {
        LogicError("Accumulation not yet supported for D&C");
    }
    else
    {
        DistMatrix<Real,MC,MR> QReal(g);

        DistMatrixWriteProxy<Real,Real,STAR,STAR> wProx( wPre );
        auto& w = wProx.Get();

        info.dcInfo =
          DivideAndConquer
          ( d_STAR_STAR.Matrix(), dSubReal.Matrix(), w, QReal, ctrl );
        herm_eig::SortAndFilter( w, QReal, ctrl );

        Copy( QReal, Q );
        DiagonalScale( LEFT, NORMAL, phase, Q );
    }

    return info;
}

template<typename Real,typename=EnableIf<IsBlasScalar<Real>>>
HermitianTridiagEigInfo
MRRRHelper
( const AbstractDistMatrix<Real>& d,
  const AbstractDistMatrix<Real>& dSub,
        AbstractDistMatrix<Real>& wPre, 
        AbstractDistMatrix<Real>& QPre, 
  const HermitianTridiagEigCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    // NOTE: The computation forces double-precision due to PMRRR limitations
    const Int n = d.Height();
    const Grid& g = d.Grid();
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

template<typename Real,typename=EnableIf<IsBlasScalar<Real>>>
HermitianTridiagEigInfo
MRRRHelper
( const AbstractDistMatrix<Real         >& d,
  const AbstractDistMatrix<Complex<Real>>& dSub,
        AbstractDistMatrix<Real         >& wPre, 
        AbstractDistMatrix<Complex<Real>>& QPre, 
  const HermitianTridiagEigCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    // NOTE: The computation forces double-precision due to PMRRR limitations
    const Int n = d.Height();
    const Grid& g = d.Grid();
    typedef Complex<Real> C;
    HermitianTridiagEigInfo info;

    DistMatrix<double,STAR,STAR> d_STAR_STAR(g);
    DistMatrix<Complex<double>,STAR,STAR> dSub_STAR_STAR(g);
    Copy( d, d_STAR_STAR );
    dSub_STAR_STAR.Resize( n-1, 1, n );
    Copy( dSub, dSub_STAR_STAR );

    DistMatrix<double,STAR,STAR> dSubReal(g);
    DistMatrix<Complex<double>,STAR,STAR> phase(g);
    RemovePhase( dSub_STAR_STAR, dSubReal, phase );

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

    Copy( QReal, Q );
    DiagonalScale( LEFT, NORMAL, phase, Q ); 

    return info;
}

template<typename Real,typename=EnableIf<IsBlasScalar<Real>>>
HermitianTridiagEigInfo
Helper
( const AbstractDistMatrix<Real>& d,
  const AbstractDistMatrix<Real>& dSub,
        AbstractDistMatrix<Real>& w, 
        AbstractDistMatrix<Real>& Q, 
  const HermitianTridiagEigCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    if( ctrl.alg == HERM_TRIDIAG_EIG_QR )
    {
        return QRHelper( d, dSub, w, Q, ctrl );
    }
    else if( ctrl.alg == HERM_TRIDIAG_EIG_DC )
    {
        return DCHelper( d, dSub, w, Q, ctrl );
    }
    else
    {
        return MRRRHelper( d, dSub, w, Q, ctrl );
    }
}

template<typename Real,typename=DisableIf<IsBlasScalar<Real>>,typename=void>
HermitianTridiagEigInfo
Helper
( const AbstractDistMatrix<Real>& d,
  const AbstractDistMatrix<Real>& dSub,
        AbstractDistMatrix<Real>& w, 
        AbstractDistMatrix<Real>& Q, 
  const HermitianTridiagEigCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    if( ctrl.alg == HERM_TRIDIAG_EIG_QR )
    {
        return QRHelper( d, dSub, w, Q, ctrl );
    }
    else
    {
        return DCHelper( d, dSub, w, Q, ctrl );
    }
}

template<typename Real,typename=EnableIf<IsBlasScalar<Real>>>
HermitianTridiagEigInfo
Helper
( const AbstractDistMatrix<Real         >& d,
  const AbstractDistMatrix<Complex<Real>>& dSub,
        AbstractDistMatrix<Real         >& w, 
        AbstractDistMatrix<Complex<Real>>& Q, 
  const HermitianTridiagEigCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    if( ctrl.alg == HERM_TRIDIAG_EIG_QR )
    {
        return QRHelper( d, dSub, w, Q, ctrl );
    }
    else if( ctrl.alg == HERM_TRIDIAG_EIG_DC )
    {
        return DCHelper( d, dSub, w, Q, ctrl );
    }
    else
    {
        return MRRRHelper( d, dSub, w, Q, ctrl );
    }
}

template<typename Real,typename=DisableIf<IsBlasScalar<Real>>,typename=void>
HermitianTridiagEigInfo
Helper
( const AbstractDistMatrix<Real         >& d,
  const AbstractDistMatrix<Complex<Real>>& dSub,
        AbstractDistMatrix<Real         >& w, 
        AbstractDistMatrix<Complex<Real>>& QPre, 
  const HermitianTridiagEigCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    if( ctrl.alg == HERM_TRIDIAG_EIG_QR )
    {
        return QRHelper( d, dSub, w, QPre, ctrl );
    }
    else
    {
        return DCHelper( d, dSub, w, QPre, ctrl );
    }
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
    EL_DEBUG_CSE
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
    EL_DEBUG_CSE
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
    EL_DEBUG_CSE
    LogicError("MRRREstimate not yet supported for nonstandard datatypes");
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
    EL_DEBUG_CSE
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
    EL_DEBUG_CSE
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
    EL_DEBUG_CSE
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
    EL_DEBUG_CSE
    return MRRRPostEstimateHelper( d, dSub, w, Q, sort, vl, vu );
}

} // namespace herm_tridiag_eig

#define PROTO(F) \
  template void herm_eig::SortAndFilter \
  ( Matrix<Base<F>>& w, \
    Matrix<F>& Q, \
    const HermitianTridiagEigCtrl<Base<F>>& ctrl ); \
  template void herm_eig::SortAndFilter \
  ( AbstractDistMatrix<Base<F>>& w, \
    AbstractDistMatrix<F>& Q, \
    const HermitianTridiagEigCtrl<Base<F>>& ctrl ); \
  template HermitianTridiagEigInfo HermitianTridiagEig \
  ( const Matrix<Base<F>>& d, \
    const Matrix<F>& dSub, \
          Matrix<Base<F>>& w, \
    const HermitianTridiagEigCtrl<Base<F>>& ctrl ); \
  template HermitianTridiagEigInfo HermitianTridiagEig \
  ( const Matrix<Base<F>>& d, \
    const Matrix<F>& dSub, \
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
  template void herm_eig::SortAndFilter \
  ( AbstractDistMatrix<Real>& w, \
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
