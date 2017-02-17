/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

#include "./Pseudospectra/Util.hpp"
#include "./Pseudospectra/Power.hpp"
#include "./Pseudospectra/Lanczos.hpp"
#include "./Pseudospectra/IRA.hpp"
#include "./Pseudospectra/IRL.hpp"
#include "./Pseudospectra/Analytic.hpp"

// For one-norm pseudospectra. An adaptation of the more robust algorithm of
// Higham and Tisseur will hopefully be implemented soon.
#include "./Pseudospectra/HagerHigham.hpp"

namespace El {

template<typename Field>
Matrix<Int> TriangularSpectralCloud
( const Matrix<Field>& UPre,
  const Matrix<Complex<Base<Field>>>& shifts,
        Matrix<Base<Field>>& invNorms,
        PseudospecCtrl<Base<Field>> psCtrl )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    typedef Complex<Real> C;

    // Force U to be complex in as cheap of a manner as possible
    MatrixReadProxy<Field,C> UProx( UPre );
    auto& U = UProx.GetLocked();

    // Check if the off-diagonal is sufficiently small; if so, compute the
    // pseudospectra analytically from the eigenvalues. This also takes care
    // of the case where the matrix is a constant multiple of the identity
    // matrix, which, after shifting, can lead to the zero matrix, which would
    // cause problems for the Lanczos convergence criteria.
    if( pspec::TriangIsNormal( U, psCtrl.tol ) )
    {
        Matrix<Int> itCounts;
        if( psCtrl.progress )
            Output("Matrix was numerically normal");
        auto w = GetDiagonal(U);
        if( psCtrl.norm == PS_TWO_NORM )
            pspec::Analytic( w, shifts, invNorms, psCtrl.snapCtrl );
        else
            LogicError("Analytic one-norm pseudospectra not yet supported");
        Zeros( itCounts, shifts.Height(), 1 );
        return itCounts;
    }

    psCtrl.schur = true;
    if( psCtrl.norm == PS_TWO_NORM )
    {
        if( psCtrl.arnoldi )
        {
            if( psCtrl.basisSize > 1 )
                return pspec::IRA( U, shifts, invNorms, psCtrl );
            else
                return pspec::Lanczos( U, shifts, invNorms, psCtrl );
        }
        else
            return pspec::Power( U, shifts, invNorms, psCtrl );
    }
    else
        return pspec::HagerHigham( U, shifts, invNorms, psCtrl );
        // Q is assumed to be the identity
}

template<typename Field>
Matrix<Int> TriangularSpectralCloud
( const Matrix<Field>& UPre,
  const Matrix<Field>& QPre,
  const Matrix<Complex<Base<Field>>>& shifts,
        Matrix<Base<Field>>& invNorms,
        PseudospecCtrl<Base<Field>> psCtrl )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    typedef Complex<Real> C;

    // Force U to be complex as cheaply as possible
    MatrixReadProxy<Field,C> UProx( UPre );
    auto& U = UProx.GetLocked();

    // Check if the off-diagonal is sufficiently small; if so, compute the
    // pseudospectrum analytically from the eigenvalues. This also takes care
    // of the case where the matrix is a constant multiple of the identity
    // matrix, which, after shifting, can lead to the zero matrix, which would
    // cause problems for the Lanczos convergence criteria.
    if( pspec::TriangIsNormal( U, psCtrl.tol ) )
    {
        Matrix<Int> itCounts;
        if( psCtrl.progress )
            Output("Matrix was numerically normal");
        auto w = GetDiagonal(U);
        if( psCtrl.norm == PS_TWO_NORM )
            pspec::Analytic( w, shifts, invNorms, psCtrl.snapCtrl );
        else
            LogicError("Analytic one-norm pseudospectra not yet supported");
        Zeros( itCounts, shifts.Height(), 1 );
        return itCounts;
    }

    psCtrl.schur = true;
    if( psCtrl.norm == PS_TWO_NORM )
    {
        if( psCtrl.arnoldi )
        {
            if( psCtrl.basisSize > 1 )
                return pspec::IRA( U, shifts, invNorms, psCtrl );
            else
                return pspec::Lanczos( U, shifts, invNorms, psCtrl );
        }
        else
            return pspec::Power( U, shifts, invNorms, psCtrl );
    }
    else
    {
        // Force Q to be complex as cheaply as possible
        MatrixReadProxy<Field,C> QProx( QPre );
        auto& Q = QProx.GetLocked();
        return pspec::HagerHigham( U, Q, shifts, invNorms, psCtrl );
    }
}

template<typename Real>
Matrix<Int> QuasiTriangularSpectralCloud
( const Matrix<Real>& U,
  const Matrix<Complex<Real>>& shifts,
        Matrix<Real>& invNorms,
        PseudospecCtrl<Real> psCtrl )
{
    EL_DEBUG_CSE

    // Check if the off-diagonal is sufficiently small; if so, compute the
    // pseudospectrum analytically from the eigenvalues. This also takes care
    // of the case where the matrix is a constant multiple of the identity
    // matrix, which, after shifting, can lead to the zero matrix, which would
    // cause problems for the Lanczos convergence criteria.
    if( pspec::QuasiTriangIsNormal( U, psCtrl.tol ) )
    {
        Matrix<Int> itCounts;
        if( psCtrl.progress )
            Output("Matrix was numerically normal");
        const auto w = schur::QuasiTriangEig( U );
        if( psCtrl.norm == PS_TWO_NORM )
            pspec::Analytic( w, shifts, invNorms, psCtrl.snapCtrl );
        else
            LogicError("Analytic one-norm pseudospectra not yet supported");
        Zeros( itCounts, shifts.Height(), 1 );
        return itCounts;
    }

    psCtrl.schur = true;
    if( psCtrl.norm == PS_ONE_NORM )
        LogicError("This option is not yet written");
    return pspec::IRA( U, shifts, invNorms, psCtrl );
}

template<typename Real>
Matrix<Int> QuasiTriangularSpectralCloud
( const Matrix<Real>& U,
  const Matrix<Real>& Q,
  const Matrix<Complex<Real>>& shifts,
        Matrix<Real>& invNorms,
        PseudospecCtrl<Real> psCtrl )
{
    EL_DEBUG_CSE

    // Check if the off-diagonal is sufficiently small; if so, compute the
    // pseudospectrum analytically from the eigenvalues. This also takes care
    // of the case where the matrix is a constant multiple of the identity
    // matrix, which, after shifting, can lead to the zero matrix, which would
    // cause problems for the Lanczos convergence criteria.
    if( pspec::QuasiTriangIsNormal( U, psCtrl.tol ) )
    {
        Matrix<Int> itCounts;
        if( psCtrl.progress )
            Output("Matrix was numerically normal");
        auto w = schur::QuasiTriangEig( U );
        if( psCtrl.norm == PS_TWO_NORM )
            pspec::Analytic( w, shifts, invNorms, psCtrl.snapCtrl );
        else
            LogicError("Analytic one-norm pseudospectra not yet supported");
        Zeros( itCounts, shifts.Height(), 1 );
        return itCounts;
    }

    psCtrl.schur = true;
    if( psCtrl.norm == PS_ONE_NORM )
        LogicError("This option is not yet written");
    return pspec::IRA( U, shifts, invNorms, psCtrl );
}

template<typename Field>
Matrix<Int> HessenbergSpectralCloud
( const Matrix<Field>& HPre,
  const Matrix<Complex<Base<Field>>>& shifts,
        Matrix<Base<Field>>& invNorms,
        PseudospecCtrl<Base<Field>> psCtrl )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    typedef Complex<Real> C;

    // Force H to be complex as cheaply as possible
    MatrixReadProxy<Field,C> HProx( HPre );
    auto& H = HProx.GetLocked();

    // TODO: Check if the subdiagonal is numerically zero, and, if so, revert to
    //       triangular version of SpectralCloud?
    psCtrl.schur = false;
    if( psCtrl.norm == PS_TWO_NORM )
    {
        if( psCtrl.arnoldi )
        {
            if( psCtrl.basisSize > 1 )
                return pspec::IRA( H, shifts, invNorms, psCtrl );
            else
                return pspec::Lanczos( H, shifts, invNorms, psCtrl );
        }
        else
            return pspec::Power( H, shifts, invNorms, psCtrl );
    }
    else
        return pspec::HagerHigham( H, shifts, invNorms, psCtrl );
        // Q is assumed to be the identity
}

template<typename Field>
Matrix<Int> HessenbergSpectralCloud
( const Matrix<Field>& HPre,
  const Matrix<Field>& QPre,
  const Matrix<Complex<Base<Field>>>& shifts,
        Matrix<Base<Field>>& invNorms,
        PseudospecCtrl<Base<Field>> psCtrl )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    typedef Complex<Real> C;

    // Force H to be complex as cheaply as possible
    MatrixReadProxy<Field,C> HProx( HPre );
    auto& H = HProx.GetLocked();

    // TODO: Check if the subdiagonal is numerically zero, and, if so, revert to
    //       triangular version of SpectralCloud?
    psCtrl.schur = false;
    if( psCtrl.norm == PS_TWO_NORM )
    {
        if( psCtrl.arnoldi )
        {
            if( psCtrl.basisSize > 1 )
                return pspec::IRA( H, shifts, invNorms, psCtrl );
            else
                return pspec::Lanczos( H, shifts, invNorms, psCtrl );
        }
        else
            return pspec::Power( H, shifts, invNorms, psCtrl );
    }
    else
    {
        // Force Q to be complex as cheaply as possible
        MatrixReadProxy<Field,C> QProx( QPre );
        auto& Q = QProx.GetLocked();
        return pspec::HagerHigham( H, Q, shifts, invNorms, psCtrl );
    }
}

template<typename Field>
DistMatrix<Int,VR,STAR> TriangularSpectralCloud
( const AbstractDistMatrix<Field>& UPre,
  const AbstractDistMatrix<Complex<Base<Field>>>& shiftsPre,
        AbstractDistMatrix<Base<Field>>& invNorms,
        PseudospecCtrl<Base<Field>> psCtrl )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    typedef Complex<Real> C;
    const Grid& g = UPre.Grid();

    // Force 'U' to be complex and in a [MC,MR] distribution
    DistMatrixReadProxy<Field,C,MC,MR> UProx( UPre );
    auto& U = UProx.GetLocked();

    // Force 'shifts' to be in a [VR,STAR] distribution
    DistMatrixReadProxy<C,C,VR,STAR> shiftsProx( shiftsPre );
    auto& shifts = shiftsProx.GetLocked();

    // Check if the off-diagonal is sufficiently small; if so, compute the
    // pseudospectrum analytically from the eigenvalues. This also takes care
    // of the case where the matrix is a constant multiple of the identity
    // matrix, which, after shifting, can lead to the zero matrix, which would
    // cause problems for the Lanczos convergence criteria.
    if( pspec::TriangIsNormal( U, psCtrl.tol ) )
    {
        DistMatrix<Int,VR,STAR> itCounts(g);
        if( psCtrl.progress && g.Rank() == 0 )
            Output("Matrix was numerically normal");
        auto w = GetDiagonal(U);
        if( psCtrl.norm == PS_TWO_NORM )
            pspec::Analytic( w, shifts, invNorms, psCtrl.snapCtrl );
        else
            LogicError("Analytic one-norm pseudospectra not yet supported");
        itCounts.AlignWith( shifts );
        Zeros( itCounts, shifts.Height(), 1 );
        return itCounts;
    }

    psCtrl.schur = true;
    if( psCtrl.norm == PS_TWO_NORM )
    {
        if( psCtrl.arnoldi )
        {
            if( psCtrl.basisSize > 1 )
                return pspec::IRA( U, shifts, invNorms, psCtrl );
            else
                return pspec::Lanczos( U, shifts, invNorms, psCtrl );
        }
        else
            return pspec::Power( U, shifts, invNorms, psCtrl );
    }
    else
        return pspec::HagerHigham( U, shifts, invNorms, psCtrl );
}

template<typename Field>
DistMatrix<Int,VR,STAR> TriangularSpectralCloud
( const AbstractDistMatrix<Field>& UPre,
  const AbstractDistMatrix<Field>& QPre,
  const AbstractDistMatrix<Complex<Base<Field>>>& shiftsPre,
        AbstractDistMatrix<Base<Field>>& invNorms,
        PseudospecCtrl<Base<Field>> psCtrl )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    typedef Complex<Real> C;
    const Grid& g = UPre.Grid();

    // Force 'U' to be complex and in a [MC,MR] distribution
    DistMatrixReadProxy<Field,C,MC,MR> UProx( UPre );
    auto& U = UProx.GetLocked();

    // Force 'shifts' to be in a [VR,STAR] distribution
    DistMatrixReadProxy<C,C,VR,STAR> shiftsProx( shiftsPre );
    auto& shifts = shiftsProx.GetLocked();

    // Check if the off-diagonal is sufficiently small; if so, compute the
    // pseudospectrum analytically from the eigenvalues. This also takes care
    // of the case where the matrix is a constant multiple of the identity
    // matrix, which, after shifting, can lead to the zero matrix, which would
    // cause problems for the Lanczos convergence criteria.
    if( pspec::TriangIsNormal( U, psCtrl.tol ) )
    {
        DistMatrix<Int,VR,STAR> itCounts(g);
        if( psCtrl.progress && g.Rank() == 0 )
            Output("Matrix was numerically normal");
        auto w = GetDiagonal(U);
        if( psCtrl.norm == PS_TWO_NORM )
            pspec::Analytic( w, shifts, invNorms, psCtrl.snapCtrl );
        else
            LogicError("Analytic one-norm pseudospectra not yet supported");
        itCounts.AlignWith( shifts );
        Zeros( itCounts, shifts.Height(), 1 );
        return itCounts;
    }

    psCtrl.schur = true;
    if( psCtrl.norm == PS_TWO_NORM )
    {
        if( psCtrl.arnoldi )
        {
            if( psCtrl.basisSize > 1 )
                return pspec::IRA( U, shifts, invNorms, psCtrl );
            else
                return pspec::Lanczos( U, shifts, invNorms, psCtrl );
        }
        else
            return pspec::Power( U, shifts, invNorms, psCtrl );
    }
    else
    {
        // Force 'Q' to be complex and in a [MC,MR] distribution as cheaply
        // as possible
        DistMatrixReadProxy<Field,C,MC,MR> QProx( QPre );
        auto& Q = QProx.GetLocked();
        return pspec::HagerHigham( U, Q, shifts, invNorms, psCtrl );
    }
}

template<typename Real>
DistMatrix<Int,VR,STAR> QuasiTriangularSpectralCloud
( const AbstractDistMatrix<Real>& UPre,
  const AbstractDistMatrix<Complex<Real>>& shiftsPre,
        AbstractDistMatrix<Real>& invNorms,
        PseudospecCtrl<Real> psCtrl )
{
    EL_DEBUG_CSE
    typedef Complex<Real> C;
    const Grid& g = UPre.Grid();

    // Force 'U' to be in a [MC,MR] distribution
    DistMatrixReadProxy<Real,Real,MC,MR> UProx( UPre );
    auto& U = UProx.GetLocked();

    // Force 'shifts' to be in a [VR,STAR] distribution
    DistMatrixReadProxy<C,C,VR,STAR> shiftsProx( shiftsPre );
    auto& shifts = shiftsProx.GetLocked();

    // Check if the off-diagonal is sufficiently small; if so, compute the
    // pseudospectrum analytically from the eigenvalues. This also takes care
    // of the case where the matrix is a constant multiple of the identity
    // matrix, which, after shifting, can lead to the zero matrix, which would
    // cause problems for the Lanczos convergence criteria.
    if( pspec::QuasiTriangIsNormal( U, psCtrl.tol ) )
    {
        DistMatrix<Int,VR,STAR> itCounts(g);
        if( psCtrl.progress && g.Rank() == 0 )
            Output("Matrix was numerically normal");
        auto w = schur::QuasiTriangEig( U );
        if( psCtrl.norm == PS_TWO_NORM )
            pspec::Analytic( w, shifts, invNorms, psCtrl.snapCtrl );
        else
            LogicError("Analytic one-norm pseudospectra not yet supported");
        itCounts.AlignWith( shifts );
        Zeros( itCounts, shifts.Height(), 1 );
        return itCounts;
    }

    psCtrl.schur = true;
    if( psCtrl.norm == PS_ONE_NORM )
        LogicError("This option is not yet written");
    return pspec::IRA( U, shifts, invNorms, psCtrl );
}

template<typename Real>
DistMatrix<Int,VR,STAR> QuasiTriangularSpectralCloud
( const AbstractDistMatrix<Real>& UPre,
  const AbstractDistMatrix<Real>& QPre,
  const AbstractDistMatrix<Complex<Real>>& shiftsPre,
        AbstractDistMatrix<Real>& invNorms,
        PseudospecCtrl<Real> psCtrl )
{
    EL_DEBUG_CSE
    typedef Complex<Real> C;
    const Grid& g = UPre.Grid();

    // Force 'U' to be in a [MC,MR] distribution
    DistMatrixReadProxy<Real,Real,MC,MR> UProx( UPre );
    auto& U = UProx.GetLocked();

    // Force 'shifts' to be in a [VR,STAR] distribution
    DistMatrixReadProxy<C,C,VR,STAR> shiftsProx( shiftsPre );
    auto& shifts = shiftsProx.GetLocked();

    // Check if the off-diagonal is sufficiently small; if so, compute the
    // pseudospectrum analytically from the eigenvalues. This also takes care
    // of the case where the matrix is a constant multiple of the identity
    // matrix, which, after shifting, can lead to the zero matrix, which would
    // cause problems for the Lanczos convergence criteria.
    if( pspec::QuasiTriangIsNormal( U, psCtrl.tol ) )
    {
        DistMatrix<Int,VR,STAR> itCounts(g);
        if( psCtrl.progress && g.Rank() == 0 )
            Output("Matrix was numerically normal");
        auto w = schur::QuasiTriangEig( U );
        if( psCtrl.norm == PS_TWO_NORM )
            pspec::Analytic( w, shifts, invNorms, psCtrl.snapCtrl );
        else
            LogicError("Analytic one-norm pseudospectra not yet supported");
        itCounts.AlignWith( shifts );
        Zeros( itCounts, shifts.Height(), 1 );
        return itCounts;
    }

    psCtrl.schur = true;
    if( psCtrl.norm == PS_ONE_NORM )
        LogicError("This option is not yet written");
    return pspec::IRA( U, shifts, invNorms, psCtrl );
}

template<typename Field>
DistMatrix<Int,VR,STAR> HessenbergSpectralCloud
( const AbstractDistMatrix<Field>& HPre,
  const AbstractDistMatrix<Complex<Base<Field>>>& shiftsPre,
        AbstractDistMatrix<Base<Field>>& invNorms,
        PseudospecCtrl<Base<Field>> psCtrl )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    typedef Complex<Real> C;

    // Force 'H' to be complex in a [MC,MR] distribution
    DistMatrixReadProxy<Field,C,MC,MR> HProx( HPre );
    auto& H = HProx.GetLocked();

    // Force 'shifts' to be in a [VR,STAR] distribution
    DistMatrixReadProxy<C,C,VR,STAR> shiftsProx( shiftsPre );
    auto& shifts = shiftsProx.GetLocked();

    // TODO: Check if the subdiagonal is sufficiently small, and, if so, revert
    //       to TriangularSpectralCloud
    psCtrl.schur = false;
    if( psCtrl.norm == PS_TWO_NORM )
    {
        if( psCtrl.arnoldi )
        {
            if( psCtrl.basisSize > 1 )
                return pspec::IRA( H, shifts, invNorms, psCtrl );
            else
                return pspec::Lanczos( H, shifts, invNorms, psCtrl );
        }
        else
            return pspec::Power( H, shifts, invNorms, psCtrl );
    }
    else
        return pspec::HagerHigham( H, shifts, invNorms, psCtrl );
}

template<typename Field>
DistMatrix<Int,VR,STAR> HessenbergSpectralCloud
( const AbstractDistMatrix<Field>& HPre,
  const AbstractDistMatrix<Field>& QPre,
  const AbstractDistMatrix<Complex<Base<Field>>>& shiftsPre,
        AbstractDistMatrix<Base<Field>>& invNorms,
        PseudospecCtrl<Base<Field>> psCtrl )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    typedef Complex<Real> C;

    // Force 'H' to be complex and in a [MC,MR] distribution
    DistMatrixReadProxy<Field,C,MC,MR> HProx( HPre );
    auto& H = HProx.GetLocked();

    // Force 'shifts' to be in a [VR,STAR] distribution
    DistMatrixReadProxy<C,C,VR,STAR> shiftsProx( shiftsPre );
    auto& shifts = shiftsProx.GetLocked();

    // TODO: Check if the subdiagonal is sufficiently small, and, if so, revert
    //       to TriangularSpectralCloud
    psCtrl.schur = false;
    if( psCtrl.norm == PS_TWO_NORM )
    {
        if( psCtrl.arnoldi )
        {
            if( psCtrl.basisSize > 1 )
                return pspec::IRA( H, shifts, invNorms, psCtrl );
            else
                return pspec::Lanczos( H, shifts, invNorms, psCtrl );
        }
        else
            return pspec::Power( H, shifts, invNorms, psCtrl );
    }
    else
    {
        // Force 'Q' to be complex and in a [MC,MR] distribution
        DistMatrixReadProxy<Field,C,MC,MR> QProx( QPre );
        auto& Q = QProx.GetLocked();
        return pspec::HagerHigham( H, Q, shifts, invNorms, psCtrl );
    }
}

namespace pspec {

template<typename Real>
Matrix<Int> Helper
( const Matrix<Real>& A,
  const Matrix<Complex<Real>>& shifts,
        Matrix<Real>& invNorms,
        PseudospecCtrl<Real> psCtrl )
{
    EL_DEBUG_CSE
    typedef Complex<Real> C;

    if( psCtrl.forceComplexSchur )
    {
        Matrix<C> ACpx;
        Copy( A, ACpx );
        return SpectralCloud( ACpx, shifts, invNorms, psCtrl );
    }

    if( !psCtrl.schur )
        LogicError("Real Hessenberg algorithm not yet supported");
    Matrix<Real> U( A );
    Matrix<C> w;
    auto schurCtrl( psCtrl.schurCtrl );
    schurCtrl.hessSchurCtrl.fullTriangle = true;
    if( psCtrl.norm == PS_TWO_NORM )
    {

        Schur( U, w, schurCtrl );
        if( psCtrl.forceComplexPs )
        {
            Matrix<C> UCpx;
            schur::RealToComplex( U, UCpx );
            return TriangularSpectralCloud( UCpx, shifts, invNorms, psCtrl );
        }
        return QuasiTriangularSpectralCloud( U, shifts, invNorms, psCtrl );
    }
    else
    {
        Matrix<Real> Q;
        Schur( U, w, Q, schurCtrl );
        if( psCtrl.forceComplexPs )
        {
            Matrix<C> UCpx, QCpx;
            schur::RealToComplex( U, Q, UCpx, QCpx );
            return TriangularSpectralCloud
                   ( UCpx, QCpx, shifts, invNorms, psCtrl );
        }
        return QuasiTriangularSpectralCloud( U, Q, shifts, invNorms, psCtrl );
    }
}

template<typename Real>
DistMatrix<Int,VR,STAR> Helper
( const AbstractDistMatrix<Real>& A,
  const AbstractDistMatrix<Complex<Real>>& shifts,
        AbstractDistMatrix<Real>& invNorms,
        PseudospecCtrl<Real> psCtrl )
{
    EL_DEBUG_CSE
    typedef Complex<Real> C;
    const Grid& g = A.Grid();

    if( psCtrl.forceComplexSchur )
    {
        DistMatrix<C> ACpx(g);
        Copy( A, ACpx );
        return SpectralCloud( ACpx, shifts, invNorms, psCtrl );
    }

    if( !psCtrl.schur )
        LogicError("Real Hessenberg algorithm not yet supported");
    DistMatrix<Real> U( A );
    DistMatrix<C,VR,STAR> w(g);
    auto schurCtrl( psCtrl.schurCtrl );
    schurCtrl.hessSchurCtrl.fullTriangle = true;
    if( psCtrl.norm == PS_TWO_NORM )
    {
        Schur( U, w, schurCtrl );
        if( psCtrl.forceComplexPs )
        {
            DistMatrix<C> UCpx(g);
            schur::RealToComplex( U, UCpx );
            return TriangularSpectralCloud( UCpx, shifts, invNorms, psCtrl );
        }
        return QuasiTriangularSpectralCloud( U, shifts, invNorms, psCtrl );
    }
    else
    {
        DistMatrix<Real> Q(g);
        Schur( U, w, Q, schurCtrl );
        if( psCtrl.forceComplexPs )
        {
            DistMatrix<C> UCpx(g), QCpx(g);
            schur::RealToComplex( U, Q, UCpx, QCpx );
            return TriangularSpectralCloud
                   ( UCpx, QCpx, shifts, invNorms, psCtrl );
        }
        return QuasiTriangularSpectralCloud( U, Q, shifts, invNorms, psCtrl );
    }
}

template<typename Real>
Matrix<Int> Helper
( const Matrix<Complex<Real>>& A,
  const Matrix<Complex<Real>>& shifts,
        Matrix<Real>& invNorms,
        PseudospecCtrl<Real> psCtrl )
{
    EL_DEBUG_CSE
    typedef Complex<Real> C;

    Matrix<C> U( A );
    if( psCtrl.norm == PS_TWO_NORM )
    {
        if( psCtrl.schur )
        {
            Matrix<C> w;
            auto schurCtrl( psCtrl.schurCtrl );
            schurCtrl.hessSchurCtrl.fullTriangle = true;
            Schur( U, w, schurCtrl );
            return TriangularSpectralCloud( U, shifts, invNorms, psCtrl );
        }
        else
        {
            hessenberg::ExplicitCondensed( UPPER, U );
            return HessenbergSpectralCloud( U, shifts, invNorms, psCtrl );
        }
    }
    else
    {
        Matrix<C> Q;
        if( psCtrl.schur )
        {
            Matrix<C> w;
            auto schurCtrl( psCtrl.schurCtrl );
            schurCtrl.hessSchurCtrl.fullTriangle = true;
            Schur( U, w, Q, schurCtrl );
            return TriangularSpectralCloud( U, Q, shifts, invNorms, psCtrl );
        }
        else
        {
            Matrix<C> t;
            Hessenberg( UPPER, U, t );
            Identity( Q, A.Height(), A.Height() );
            hessenberg::ApplyQ( LEFT, UPPER, NORMAL, U, t, Q );
            return HessenbergSpectralCloud( U, Q, shifts, invNorms, psCtrl );
        }
    }
}

template<typename Real>
DistMatrix<Int,VR,STAR> Helper
( const AbstractDistMatrix<Complex<Real>>& A,
  const AbstractDistMatrix<Complex<Real>>& shifts,
        AbstractDistMatrix<Real>& invNorms,
        PseudospecCtrl<Real> psCtrl )
{
    EL_DEBUG_CSE
    typedef Complex<Real> C;

    const Grid& g = A.Grid();
    DistMatrix<C> U( A );

    if( psCtrl.norm == PS_TWO_NORM )
    {
        if( psCtrl.schur )
        {
            DistMatrix<C,VR,STAR> w(g);
            auto schurCtrl( psCtrl.schurCtrl );
            schurCtrl.hessSchurCtrl.fullTriangle = true;
            Schur( U, w, schurCtrl );
            return TriangularSpectralCloud( U, shifts, invNorms, psCtrl );
        }
        else
        {
            hessenberg::ExplicitCondensed( UPPER, U );
            return HessenbergSpectralCloud( U, shifts, invNorms, psCtrl );
        }
    }
    else
    {
        DistMatrix<C> Q(g);
        if( psCtrl.schur )
        {
            DistMatrix<C,VR,STAR> w(g);
            auto schurCtrl( psCtrl.schurCtrl );
            schurCtrl.hessSchurCtrl.fullTriangle = true;
            Schur( U, w, Q, schurCtrl );
            return TriangularSpectralCloud( U, Q, shifts, invNorms, psCtrl );
        }
        else
        {
            DistMatrix<C,STAR,STAR> t(g);
            Hessenberg( UPPER, U, t );
            Identity( Q, U.Height(), U.Height() );
            hessenberg::ApplyQ( LEFT, UPPER, NORMAL, U, t, Q );
            return HessenbergSpectralCloud( U, Q, shifts, invNorms, psCtrl );
        }
    }
}

} // namespace pspec

template<typename Field>
Matrix<Int> SpectralCloud
( const Matrix<Field>& A,
  const Matrix<Complex<Base<Field>>>& shifts,
        Matrix<Base<Field>>& invNorms,
        PseudospecCtrl<Base<Field>> psCtrl )
{
    EL_DEBUG_CSE
    return pspec::Helper( A, shifts, invNorms, psCtrl );
}

template<typename Field>
DistMatrix<Int,VR,STAR> SpectralCloud
( const AbstractDistMatrix<Field>& A,
  const AbstractDistMatrix<Complex<Base<Field>>>& shifts,
        AbstractDistMatrix<Base<Field>>& invNorms,
        PseudospecCtrl<Base<Field>> psCtrl )
{
    EL_DEBUG_CSE
    return pspec::Helper( A, shifts, invNorms, psCtrl );
}

// Treat each pixel as being located a cell center and tesselate a box with
// said square cells
template<typename Field>
Matrix<Int> TriangularSpectralWindow
( const Matrix<Field>& U,
        Matrix<Base<Field>>& invNormMap,
  Complex<Base<Field>> center,
  Base<Field> realWidth,
  Base<Field> imagWidth,
  Int realSize,
  Int imagSize,
  PseudospecCtrl<Base<Field>> psCtrl )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    typedef Complex<Real> C;

    psCtrl.snapCtrl.realSize = realSize;
    psCtrl.snapCtrl.imagSize = imagSize;

    psCtrl.center = center;
    psCtrl.realWidth = realWidth;
    psCtrl.imagWidth = imagWidth;

    const Real realStep = realWidth/realSize;
    const Real imagStep = imagWidth/imagSize;
    const C corner = center + C(-realWidth/2,imagWidth/2);
    Matrix<C> shifts( realSize*imagSize, 1 );
    for( Int j=0; j<realSize*imagSize; ++j )
    {
        const Int x = j / imagSize;
        const Int y = j % imagSize;
        shifts.Set( j, 0, corner+C((x+0.5)*realStep,-(y+0.5)*imagStep) );
    }

    // Form the vector of invNorms
    Matrix<Real> invNorms;
    auto itCounts = TriangularSpectralCloud( U, shifts, invNorms, psCtrl );

    // Rearrange the vectors into grids
    Matrix<Int> itCountMap;
    pspec::ReshapeIntoGrid( realSize, imagSize, invNorms, invNormMap );
    pspec::ReshapeIntoGrid( realSize, imagSize, itCounts, itCountMap );
    return itCountMap;
}

template<typename Field>
Matrix<Int> TriangularSpectralWindow
( const Matrix<Field>& U,
  const Matrix<Field>& Q,
        Matrix<Base<Field>>& invNormMap,
  Complex<Base<Field>> center,
  Base<Field> realWidth,
  Base<Field> imagWidth,
  Int realSize,
  Int imagSize,
  PseudospecCtrl<Base<Field>> psCtrl )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    typedef Complex<Real> C;

    psCtrl.snapCtrl.realSize = realSize;
    psCtrl.snapCtrl.imagSize = imagSize;

    psCtrl.center = center;
    psCtrl.realWidth = realWidth;
    psCtrl.imagWidth = imagWidth;

    const Real realStep = realWidth/realSize;
    const Real imagStep = imagWidth/imagSize;
    const C corner = center + C(-realWidth/2,imagWidth/2);
    Matrix<C> shifts( realSize*imagSize, 1 );
    for( Int j=0; j<realSize*imagSize; ++j )
    {
        const Int x = j / imagSize;
        const Int y = j % imagSize;
        shifts.Set( j, 0, corner+C((x+0.5)*realStep,-(y+0.5)*imagStep) );
    }

    // Form the vector of invNorms
    Matrix<Real> invNorms;
    auto itCounts = TriangularSpectralCloud( U, Q, shifts, invNorms, psCtrl );

    // Rearrange the vectors into grids
    Matrix<Int> itCountMap;
    pspec::ReshapeIntoGrid( realSize, imagSize, invNorms, invNormMap );
    pspec::ReshapeIntoGrid( realSize, imagSize, itCounts, itCountMap );
    return itCountMap;
}

template<typename Real>
Matrix<Int> QuasiTriangularSpectralWindow
( const Matrix<Real>& U,
        Matrix<Real>& invNormMap,
  Complex<Real> center,
  Real realWidth,
  Real imagWidth,
  Int realSize,
  Int imagSize,
  PseudospecCtrl<Real> psCtrl )
{
    EL_DEBUG_CSE
    typedef Complex<Real> C;

    psCtrl.snapCtrl.realSize = realSize;
    psCtrl.snapCtrl.imagSize = imagSize;

    psCtrl.center = center;
    psCtrl.realWidth = realWidth;
    psCtrl.imagWidth = imagWidth;

    const Real realStep = realWidth/realSize;
    const Real imagStep = imagWidth/imagSize;
    const C corner = center + C(-realWidth/2,imagWidth/2);
    Matrix<C> shifts( realSize*imagSize, 1 );
    for( Int j=0; j<realSize*imagSize; ++j )
    {
        const Int x = j / imagSize;
        const Int y = j % imagSize;
        shifts.Set( j, 0, corner+C((x+0.5)*realStep,-(y+0.5)*imagStep) );
    }

    // Form the vector of invNorms
    Matrix<Real> invNorms;
    auto itCounts = QuasiTriangularSpectralCloud( U, shifts, invNorms, psCtrl );

    // Rearrange the vectors into grids
    Matrix<Int> itCountMap;
    pspec::ReshapeIntoGrid( realSize, imagSize, invNorms, invNormMap );
    pspec::ReshapeIntoGrid( realSize, imagSize, itCounts, itCountMap );
    return itCountMap;
}

template<typename Real>
Matrix<Int> QuasiTriangularSpectralWindow
( const Matrix<Real>& U,
  const Matrix<Real>& Q,
        Matrix<Real>& invNormMap,
  Complex<Real> center,
  Real realWidth,
  Real imagWidth,
  Int realSize,
  Int imagSize,
  PseudospecCtrl<Real> psCtrl )
{
    EL_DEBUG_CSE
    typedef Complex<Real> C;

    psCtrl.snapCtrl.realSize = realSize;
    psCtrl.snapCtrl.imagSize = imagSize;

    psCtrl.center = center;
    psCtrl.realWidth = realWidth;
    psCtrl.imagWidth = imagWidth;

    const Real realStep = realWidth/realSize;
    const Real imagStep = imagWidth/imagSize;
    const C corner = center + C(-realWidth/2,imagWidth/2);
    Matrix<C> shifts( realSize*imagSize, 1 );
    for( Int j=0; j<realSize*imagSize; ++j )
    {
        const Int x = j / imagSize;
        const Int y = j % imagSize;
        shifts.Set( j, 0, corner+C((x+0.5)*realStep,-(y+0.5)*imagStep) );
    }

    // Form the vector of invNorms
    Matrix<Real> invNorms;
    auto itCounts =
        QuasiTriangularSpectralCloud( U, Q, shifts, invNorms, psCtrl );

    // Rearrange the vectors into grids
    Matrix<Int> itCountMap;
    pspec::ReshapeIntoGrid( realSize, imagSize, invNorms, invNormMap );
    pspec::ReshapeIntoGrid( realSize, imagSize, itCounts, itCountMap );
    return itCountMap;
}

template<typename Field>
Matrix<Int> HessenbergSpectralWindow
( const Matrix<Field>& H,
        Matrix<Base<Field>>& invNormMap,
  Complex<Base<Field>> center,
  Base<Field> realWidth,
  Base<Field> imagWidth,
  Int realSize,
  Int imagSize,
  PseudospecCtrl<Base<Field>> psCtrl )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    typedef Complex<Real> C;

    psCtrl.snapCtrl.realSize = realSize;
    psCtrl.snapCtrl.imagSize = imagSize;

    psCtrl.center = center;
    psCtrl.realWidth = realWidth;
    psCtrl.imagWidth = imagWidth;

    const Real realStep = realWidth/realSize;
    const Real imagStep = imagWidth/imagSize;
    const C corner = center + C(-realWidth/2,imagWidth/2);
    Matrix<C> shifts( realSize*imagSize, 1 );
    for( Int j=0; j<realSize*imagSize; ++j )
    {
        const Int x = j / imagSize;
        const Int y = j % imagSize;
        shifts.Set( j, 0, corner+C((x+0.5)*realStep,-(y+0.5)*imagStep) );
    }

    // Form the vector of invNorms
    Matrix<Real> invNorms;
    auto itCounts = HessenbergSpectralCloud( H, shifts, invNorms, psCtrl );

    // Rearrange the vectors into grids
    Matrix<Int> itCountMap;
    pspec::ReshapeIntoGrid( realSize, imagSize, invNorms, invNormMap );
    pspec::ReshapeIntoGrid( realSize, imagSize, itCounts, itCountMap );
    return itCountMap;
}

template<typename Field>
Matrix<Int> HessenbergSpectralWindow
( const Matrix<Field>& H,
  const Matrix<Field>& Q,
        Matrix<Base<Field>>& invNormMap,
  Complex<Base<Field>> center,
  Base<Field> realWidth,
  Base<Field> imagWidth,
  Int realSize,
  Int imagSize,
  PseudospecCtrl<Base<Field>> psCtrl )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    typedef Complex<Real> C;

    psCtrl.snapCtrl.realSize = realSize;
    psCtrl.snapCtrl.imagSize = imagSize;

    psCtrl.center = center;
    psCtrl.realWidth = realWidth;
    psCtrl.imagWidth = imagWidth;

    const Real realStep = realWidth/realSize;
    const Real imagStep = imagWidth/imagSize;
    const C corner = center + C(-realWidth/2,imagWidth/2);
    Matrix<C> shifts( realSize*imagSize, 1 );
    for( Int j=0; j<realSize*imagSize; ++j )
    {
        const Int x = j / imagSize;
        const Int y = j % imagSize;
        shifts.Set( j, 0, corner+C((x+0.5)*realStep,-(y+0.5)*imagStep) );
    }

    // Form the vector of invNorms
    Matrix<Real> invNorms;
    auto itCounts = HessenbergSpectralCloud( H, Q, shifts, invNorms, psCtrl );

    // Rearrange the vectors into grids
    Matrix<Int> itCountMap;
    pspec::ReshapeIntoGrid( realSize, imagSize, invNorms, invNormMap );
    pspec::ReshapeIntoGrid( realSize, imagSize, itCounts, itCountMap );
    return itCountMap;
}

template<typename Field>
DistMatrix<Int> TriangularSpectralWindow
( const AbstractDistMatrix<Field>& U,
        AbstractDistMatrix<Base<Field>>& invNormMap,
  Complex<Base<Field>> center,
  Base<Field> realWidth,
  Base<Field> imagWidth,
  Int realSize,
  Int imagSize,
  PseudospecCtrl<Base<Field>> psCtrl )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    typedef Complex<Real> C;
    const Grid& g = U.Grid();

    psCtrl.snapCtrl.realSize = realSize;
    psCtrl.snapCtrl.imagSize = imagSize;

    psCtrl.center = center;
    psCtrl.realWidth = realWidth;
    psCtrl.imagWidth = imagWidth;

    const Real realStep = realWidth/realSize;
    const Real imagStep = imagWidth/imagSize;
    const C corner = center + C(-realWidth/2,imagWidth/2);
    DistMatrix<C,VR,STAR> shifts( realSize*imagSize, 1, g );
    const Int numLocShifts = shifts.LocalHeight();
    for( Int iLoc=0; iLoc<numLocShifts; ++iLoc )
    {
        const Int i = shifts.GlobalRow(iLoc);
        const Int x = i / imagSize;
        const Int y = i % imagSize;
        shifts.SetLocal
        ( iLoc, 0, corner+C((x+0.5)*realStep,-(y+0.5)*imagStep) );
    }

    // Form the vector of invNorms
    DistMatrix<Real,VR,STAR> invNorms(g);
    auto itCounts = TriangularSpectralCloud( U, shifts, invNorms, psCtrl );

    // Rearrange the vectors into grids
    DistMatrix<Int> itCountMap(g);
    pspec::ReshapeIntoGrid( realSize, imagSize, invNorms, invNormMap );
    pspec::ReshapeIntoGrid( realSize, imagSize, itCounts, itCountMap );
    return itCountMap;
}

template<typename Field>
DistMatrix<Int> TriangularSpectralWindow
( const AbstractDistMatrix<Field>& U,
  const AbstractDistMatrix<Field>& Q,
        AbstractDistMatrix<Base<Field>>& invNormMap,
  Complex<Base<Field>> center,
  Base<Field> realWidth,
  Base<Field> imagWidth,
  Int realSize,
  Int imagSize,
  PseudospecCtrl<Base<Field>> psCtrl )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    typedef Complex<Real> C;
    const Grid& g = U.Grid();

    psCtrl.snapCtrl.realSize = realSize;
    psCtrl.snapCtrl.imagSize = imagSize;

    psCtrl.center = center;
    psCtrl.realWidth = realWidth;
    psCtrl.imagWidth = imagWidth;

    const Real realStep = realWidth/realSize;
    const Real imagStep = imagWidth/imagSize;
    const C corner = center + C(-realWidth/2,imagWidth/2);
    DistMatrix<C,VR,STAR> shifts( realSize*imagSize, 1, g );
    const Int numLocShifts = shifts.LocalHeight();
    for( Int iLoc=0; iLoc<numLocShifts; ++iLoc )
    {
        const Int i = shifts.GlobalRow(iLoc);
        const Int x = i / imagSize;
        const Int y = i % imagSize;
        shifts.SetLocal
        ( iLoc, 0, corner+C((x+0.5)*realStep,-(y+0.5)*imagStep) );
    }

    // Form the vector of invNorms
    DistMatrix<Real,VR,STAR> invNorms(g);
    auto itCounts = TriangularSpectralCloud( U, Q, shifts, invNorms, psCtrl );

    // Rearrange the vectors into grids
    DistMatrix<Int> itCountMap(g);
    pspec::ReshapeIntoGrid( realSize, imagSize, invNorms, invNormMap );
    pspec::ReshapeIntoGrid( realSize, imagSize, itCounts, itCountMap );
    return itCountMap;
}

template<typename Real>
DistMatrix<Int> QuasiTriangularSpectralWindow
( const AbstractDistMatrix<Real>& U,
        AbstractDistMatrix<Real>& invNormMap,
  Complex<Real> center,
  Real realWidth,
  Real imagWidth,
  Int realSize,
  Int imagSize,
  PseudospecCtrl<Real> psCtrl )
{
    EL_DEBUG_CSE
    typedef Complex<Real> C;
    const Grid& g = U.Grid();

    psCtrl.snapCtrl.realSize = realSize;
    psCtrl.snapCtrl.imagSize = imagSize;

    psCtrl.center = center;
    psCtrl.realWidth = realWidth;
    psCtrl.imagWidth = imagWidth;

    const Real realStep = realWidth/realSize;
    const Real imagStep = imagWidth/imagSize;
    const C corner = center + C(-realWidth/2,imagWidth/2);
    DistMatrix<C,VR,STAR> shifts( realSize*imagSize, 1, g );
    const Int numLocShifts = shifts.LocalHeight();
    for( Int iLoc=0; iLoc<numLocShifts; ++iLoc )
    {
        const Int i = shifts.GlobalRow(iLoc);
        const Int x = i / imagSize;
        const Int y = i % imagSize;
        shifts.SetLocal
        ( iLoc, 0, corner+C((x+0.5)*realStep,-(y+0.5)*imagStep) );
    }

    // Form the vector of invNorms
    DistMatrix<Real,VR,STAR> invNorms(g);
    auto itCounts = QuasiTriangularSpectralCloud( U, shifts, invNorms, psCtrl );

    // Rearrange the vectors into grids
    DistMatrix<Int> itCountMap(g);
    pspec::ReshapeIntoGrid( realSize, imagSize, invNorms, invNormMap );
    pspec::ReshapeIntoGrid( realSize, imagSize, itCounts, itCountMap );
    return itCountMap;
}

template<typename Real>
DistMatrix<Int> QuasiTriangularSpectralWindow
( const AbstractDistMatrix<Real>& U,
  const AbstractDistMatrix<Real>& Q,
        AbstractDistMatrix<Real>& invNormMap,
  Complex<Real> center,
  Real realWidth,
  Real imagWidth,
  Int realSize,
  Int imagSize,
  PseudospecCtrl<Real> psCtrl )
{
    EL_DEBUG_CSE
    typedef Complex<Real> C;
    const Grid& g = U.Grid();

    psCtrl.snapCtrl.realSize = realSize;
    psCtrl.snapCtrl.imagSize = imagSize;

    psCtrl.center = center;
    psCtrl.realWidth = realWidth;
    psCtrl.imagWidth = imagWidth;

    const Real realStep = realWidth/realSize;
    const Real imagStep = imagWidth/imagSize;
    const C corner = center + C(-realWidth/2,imagWidth/2);
    DistMatrix<C,VR,STAR> shifts( realSize*imagSize, 1, g );
    const Int numLocShifts = shifts.LocalHeight();
    for( Int iLoc=0; iLoc<numLocShifts; ++iLoc )
    {
        const Int i = shifts.GlobalRow(iLoc);
        const Int x = i / imagSize;
        const Int y = i % imagSize;
        shifts.SetLocal
        ( iLoc, 0, corner+C((x+0.5)*realStep,-(y+0.5)*imagStep) );
    }

    // Form the vector of invNorms
    DistMatrix<Real,VR,STAR> invNorms(g);
    auto itCounts =
        QuasiTriangularSpectralCloud( U, Q, shifts, invNorms, psCtrl );

    // Rearrange the vectors into grids
    DistMatrix<Int> itCountMap(g);
    pspec::ReshapeIntoGrid( realSize, imagSize, invNorms, invNormMap );
    pspec::ReshapeIntoGrid( realSize, imagSize, itCounts, itCountMap );
    return itCountMap;
}

template<typename Field>
DistMatrix<Int> HessenbergSpectralWindow
( const AbstractDistMatrix<Field>& H,
        AbstractDistMatrix<Base<Field>>& invNormMap,
  Complex<Base<Field>> center,
  Base<Field> realWidth,
  Base<Field> imagWidth,
  Int realSize,
  Int imagSize,
  PseudospecCtrl<Base<Field>> psCtrl )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    typedef Complex<Real> C;
    const Grid& g = H.Grid();

    psCtrl.snapCtrl.realSize = realSize;
    psCtrl.snapCtrl.imagSize = imagSize;

    psCtrl.center = center;
    psCtrl.realWidth = realWidth;
    psCtrl.imagWidth = imagWidth;

    const Real realStep = realWidth/realSize;
    const Real imagStep = imagWidth/imagSize;
    const C corner = center + C(-realWidth/2,imagWidth/2);
    DistMatrix<C,VR,STAR> shifts( realSize*imagSize, 1, g );
    const Int numLocShifts = shifts.LocalHeight();
    for( Int iLoc=0; iLoc<numLocShifts; ++iLoc )
    {
        const Int i = shifts.GlobalRow(iLoc);
        const Int x = i / imagSize;
        const Int y = i % imagSize;
        shifts.SetLocal
        ( iLoc, 0, corner+C((x+0.5)*realStep,-(y+0.5)*imagStep) );
    }

    // Form the vector of invNorms
    DistMatrix<Real,VR,STAR> invNorms(g);
    auto itCounts = HessenbergSpectralCloud( H, shifts, invNorms, psCtrl );

    // Rearrange the vectors into grids
    DistMatrix<Int> itCountMap(g);
    pspec::ReshapeIntoGrid( realSize, imagSize, invNorms, invNormMap );
    pspec::ReshapeIntoGrid( realSize, imagSize, itCounts, itCountMap );
    return itCountMap;
}

template<typename Field>
DistMatrix<Int> HessenbergSpectralWindow
( const AbstractDistMatrix<Field>& H,
  const AbstractDistMatrix<Field>& Q,
        AbstractDistMatrix<Base<Field>>& invNormMap,
  Complex<Base<Field>> center,
  Base<Field> realWidth,
  Base<Field> imagWidth,
  Int realSize,
  Int imagSize,
  PseudospecCtrl<Base<Field>> psCtrl )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    typedef Complex<Real> C;
    const Grid& g = H.Grid();

    psCtrl.snapCtrl.realSize = realSize;
    psCtrl.snapCtrl.imagSize = imagSize;

    psCtrl.center = center;
    psCtrl.realWidth = realWidth;
    psCtrl.imagWidth = imagWidth;

    const Real realStep = realWidth/realSize;
    const Real imagStep = imagWidth/imagSize;
    const C corner = center + C(-realWidth/2,imagWidth/2);
    DistMatrix<C,VR,STAR> shifts( realSize*imagSize, 1, g );
    const Int numLocShifts = shifts.LocalHeight();
    for( Int iLoc=0; iLoc<numLocShifts; ++iLoc )
    {
        const Int i = shifts.GlobalRow(iLoc);
        const Int x = i / imagSize;
        const Int y = i % imagSize;
        shifts.SetLocal
        ( iLoc, 0, corner+C((x+0.5)*realStep,-(y+0.5)*imagStep) );
    }

    // Form the vector of invNorms
    DistMatrix<Real,VR,STAR> invNorms(g);
    auto itCounts = HessenbergSpectralCloud( H, Q, shifts, invNorms, psCtrl );

    // Rearrange the vectors into grids
    DistMatrix<Int> itCountMap(g);
    pspec::ReshapeIntoGrid( realSize, imagSize, invNorms, invNormMap );
    pspec::ReshapeIntoGrid( realSize, imagSize, itCounts, itCountMap );
    return itCountMap;
}

template<typename Field>
Matrix<Int> SpectralWindow
( const Matrix<Field>& A,
        Matrix<Base<Field>>& invNormMap,
  Complex<Base<Field>> center,
  Base<Field> realWidth,
  Base<Field> imagWidth,
  Int realSize,
  Int imagSize,
  PseudospecCtrl<Base<Field>> psCtrl )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    typedef Complex<Real> C;

    psCtrl.snapCtrl.realSize = realSize;
    psCtrl.snapCtrl.imagSize = imagSize;

    psCtrl.center = center;
    psCtrl.realWidth = realWidth;
    psCtrl.imagWidth = imagWidth;

    const Real realStep = realWidth/realSize;
    const Real imagStep = imagWidth/imagSize;
    const C corner = center + C(-realWidth/2,imagWidth/2);
    Matrix<C> shifts( realSize*imagSize, 1 );
    for( Int j=0; j<realSize*imagSize; ++j )
    {
        const Int x = j / imagSize;
        const Int y = j % imagSize;
        shifts.Set( j, 0, corner+C((x+0.5)*realStep,-(y+0.5)*imagStep) );
    }

    // Form the vector of invNorms
    Matrix<Real> invNorms;
    auto itCounts = SpectralCloud( A, shifts, invNorms, psCtrl );

    // Rearrange the vectors into grids
    Matrix<Int> itCountMap;
    pspec::ReshapeIntoGrid( realSize, imagSize, invNorms, invNormMap );
    pspec::ReshapeIntoGrid( realSize, imagSize, itCounts, itCountMap );
    return itCountMap;
}

template<typename Field>
DistMatrix<Int> SpectralWindow
( const AbstractDistMatrix<Field>& A,
        AbstractDistMatrix<Base<Field>>& invNormMap,
  Complex<Base<Field>> center,
  Base<Field> realWidth,
  Base<Field> imagWidth,
  Int realSize,
  Int imagSize,
  PseudospecCtrl<Base<Field>> psCtrl )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    typedef Complex<Real> C;
    const Grid& g = A.Grid();

    psCtrl.snapCtrl.realSize = realSize;
    psCtrl.snapCtrl.imagSize = imagSize;

    psCtrl.center = center;
    psCtrl.realWidth = realWidth;
    psCtrl.imagWidth = imagWidth;

    const Real realStep = realWidth/realSize;
    const Real imagStep = imagWidth/imagSize;
    const C corner = center + C(-realWidth/2,imagWidth/2);
    DistMatrix<C,VR,STAR> shifts( realSize*imagSize, 1, g );
    const Int numLocShifts = shifts.LocalHeight();
    for( Int iLoc=0; iLoc<numLocShifts; ++iLoc )
    {
        const Int i = shifts.GlobalRow(iLoc);
        const Int x = i / imagSize;
        const Int y = i % imagSize;
        shifts.SetLocal
        ( iLoc, 0, corner+C((x+0.5)*realStep,-(y+0.5)*imagStep) );
    }

    // Form the vector of invNorms
    DistMatrix<Real,VR,STAR> invNorms(g);
    auto itCounts = SpectralCloud( A, shifts, invNorms, psCtrl );

    // Rearrange the vectors into grids
    DistMatrix<Int> itCountMap(g);
    pspec::ReshapeIntoGrid( realSize, imagSize, invNorms, invNormMap );
    pspec::ReshapeIntoGrid( realSize, imagSize, itCounts, itCountMap );
    return itCountMap;
}

template<typename Field>
Matrix<Int> TriangularSpectralPortrait
( const Matrix<Field>& U,
        Matrix<Base<Field>>& invNormMap,
  Int realSize,
  Int imagSize,
  SpectralBox<Base<Field>>& box,
  PseudospecCtrl<Base<Field>> psCtrl )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Real radius = MaxNorm( GetDiagonal(U) );
    const Real oneNorm = OneNorm( U );

    // Essentially three cases are handled here:
    // 1) The zero matrix (force the pseudospectrum width to 1)
    // 2) Typical matrices (use a small multiple of the spectral radius)
    // 3) Highly non-normal matrices (e.g., triangular with zero main diagonal)
    Real width;
    if( oneNorm == Real(0) )
    {
        width = 1;
        if( psCtrl.progress )
            Output("Setting width to 1 to handle zero matrix");
    }
    else if( radius >= 0.2*oneNorm )
    {
        width = 2.5*radius;
        if( psCtrl.progress )
            Output
            ("Setting width to ",width," based on spectral radius, ",radius);
    }
    else
    {
        width = 0.8*oneNorm;
        if( psCtrl.progress )
            Output("Setting width to ",width," based on one norm, ",oneNorm);
    }
    Complex<Real> center(0,0);

    box.center = center;
    box.realWidth = width;
    box.imagWidth = width;

    return TriangularSpectralWindow
           ( U, invNormMap, center, width, width, realSize, imagSize, psCtrl );
}

template<typename Field>
Matrix<Int> TriangularSpectralPortrait
( const Matrix<Field>& U,
  const Matrix<Field>& Q,
        Matrix<Base<Field>>& invNormMap,
  Int realSize,
  Int imagSize,
  SpectralBox<Base<Field>>& box,
  PseudospecCtrl<Base<Field>> psCtrl )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Real radius = MaxNorm( GetDiagonal(U) );
    const Real oneNorm = OneNorm( U );

    // Essentially three cases are handled here:
    // 1) The zero matrix (force the pseudospectrum width to 1)
    // 2) Typical matrices (use a small multiple of the spectral radius)
    // 3) Highly non-normal matrices (e.g., triangular with zero main diagonal)
    Real width;
    if( oneNorm == Real(0) )
    {
        width = 1;
        if( psCtrl.progress )
            Output("Setting width to 1 to handle zero matrix");
    }
    else if( radius >= oneNorm/5 )
    {
        width = 2.5*radius;
        if( psCtrl.progress )
            Output
            ("Setting width to ",width," based on spectral radius, ",radius);
    }
    else
    {
        width = 0.8*oneNorm;
        if( psCtrl.progress )
            Output("Setting width to ",width," based on one norm, ",oneNorm);
    }
    Complex<Real> center(0,0);

    box.center = center;
    box.realWidth = width;
    box.imagWidth = width;

    return TriangularSpectralWindow
           ( U, Q, invNormMap, center, width, width, realSize, imagSize,
             psCtrl );
}

template<typename Real>
Matrix<Int> QuasiTriangularSpectralPortrait
( const Matrix<Real>& U,
        Matrix<Real>& invNormMap,
  Int realSize,
  Int imagSize,
  SpectralBox<Real>& box,
  PseudospecCtrl<Real> psCtrl )
{
    EL_DEBUG_CSE

    const auto w = schur::QuasiTriangEig( U );
    const Real radius = MaxNorm( w );
    const Real oneNorm = OneNorm( U );

    // Essentially three cases are handled here:
    // 1) The zero matrix (force the pseudospectrum width to 1)
    // 2) Typical matrices (use a small multiple of the spectral radius)
    // 3) Highly non-normal matrices (e.g., triangular with zero main diagonal)
    Real width;
    if( oneNorm == Real(0) )
    {
        width = 1;
        if( psCtrl.progress )
            Output("Setting width to 1 to handle zero matrix");
    }
    else if( radius >= 0.2*oneNorm )
    {
        width = 2.5*radius;
        if( psCtrl.progress )
            Output
            ("Setting width to ",width," based on spectral radius, ",radius);
    }
    else
    {
        width = 0.8*oneNorm;
        if( psCtrl.progress )
            Output("Setting width to ",width," based on one norm, ",oneNorm);
    }
    Complex<Real> center(0,0);

    box.center = center;
    box.realWidth = width;
    box.imagWidth = width;

    return QuasiTriangularSpectralWindow
           ( U, invNormMap, center, width, width, realSize, imagSize, psCtrl );
}

template<typename Real>
Matrix<Int> QuasiTriangularSpectralPortrait
( const Matrix<Real>& U,
  const Matrix<Real>& Q,
        Matrix<Real>& invNormMap,
  Int realSize,
  Int imagSize,
  SpectralBox<Real>& box,
  PseudospecCtrl<Real> psCtrl )
{
    EL_DEBUG_CSE

    const auto w = schur::QuasiTriangEig( U );
    const Real radius = MaxNorm( w );
    const Real oneNorm = OneNorm( U );

    // Essentially three cases are handled here:
    // 1) The zero matrix (force the pseudospectrum width to 1)
    // 2) Typical matrices (use a small multiple of the spectral radius)
    // 3) Highly non-normal matrices (e.g., triangular with zero main diagonal)
    Real width;
    if( oneNorm == Real(0) )
    {
        width = 1;
        if( psCtrl.progress )
            Output("Setting width to 1 to handle zero matrix");
    }
    else if( radius >= 0.2*oneNorm )
    {
        width = 2.5*radius;
        if( psCtrl.progress )
            Output
            ("Setting width to ",width," based on spectral radius, ",radius);
    }
    else
    {
        width = 0.8*oneNorm;
        if( psCtrl.progress )
            Output
            ("Setting width to ",width," based on one-norm, ",oneNorm);
    }
    Complex<Real> center(0,0);

    box.center = center;
    box.realWidth = width;
    box.imagWidth = width;

    return QuasiTriangularSpectralWindow
           ( U, Q, invNormMap, center, width, width, realSize, imagSize,
             psCtrl );
}

template<typename Field>
Matrix<Int> HessenbergSpectralPortrait
( const Matrix<Field>& H,
        Matrix<Base<Field>>& invNormMap,
  Int realSize,
  Int imagSize,
  SpectralBox<Base<Field>>& box,
  PseudospecCtrl<Base<Field>> psCtrl )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Real infNorm = InfinityNorm( H );
    const Real oneNorm = OneNorm( H );

    Real width;
    if( oneNorm == Real(0) )
    {
        width = 1;
        if( psCtrl.progress )
            Output("Setting width to 1 to handle zero matrix");
    }
    else
    {
        width = 0.8*Max(oneNorm,infNorm);
        if( psCtrl.progress )
            Output
            ("Setting width to ",width," based on the one norm, ",oneNorm,
             ", and infinity norm, ",infNorm);
    }
    Complex<Real> center(0,0);

    box.center = center;
    box.realWidth = width;
    box.imagWidth = width;

    return HessenbergSpectralWindow
           ( H, invNormMap, center, width, width, realSize, imagSize, psCtrl );
}

template<typename Field>
Matrix<Int> HessenbergSpectralPortrait
( const Matrix<Field>& H,
  const Matrix<Field>& Q,
        Matrix<Base<Field>>& invNormMap,
  Int realSize,
  Int imagSize,
  SpectralBox<Base<Field>>& box,
  PseudospecCtrl<Base<Field>> psCtrl )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Real infNorm = InfinityNorm( H );
    const Real oneNorm = OneNorm( H );

    Real width;
    if( oneNorm == Real(0) )
    {
        width = 1;
        if( psCtrl.progress )
            Output("Setting width to 1 to handle zero matrix");
    }
    else
    {
        width = 0.8*Max(oneNorm,infNorm);
        if( psCtrl.progress )
            Output
            ("Setting width to ",width," based on the one norm, ",oneNorm,
             ", and infinity norm, ",infNorm);
    }
    Complex<Real> center(0,0);

    box.center = center;
    box.realWidth = width;
    box.imagWidth = width;

    return HessenbergSpectralWindow
           ( H, Q, invNormMap, center, width, width, realSize, imagSize,
             psCtrl );
}

template<typename Field>
DistMatrix<Int> TriangularSpectralPortrait
( const AbstractDistMatrix<Field>& UPre,
        AbstractDistMatrix<Base<Field>>& invNormMap,
  Int realSize,
  Int imagSize,
  SpectralBox<Base<Field>>& box,
  PseudospecCtrl<Base<Field>> psCtrl )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Grid& g = UPre.Grid();

    // Force 'U' to be in a [MC,MR] distribution so that we can get its diagonal
    DistMatrixReadProxy<Field,Field,MC,MR> UProx( UPre );
    auto& U = UProx.GetLocked();

    const Real radius = MaxNorm( GetDiagonal(U) );
    const Real oneNorm = OneNorm( U );

    // Essentially three cases are handled here:
    // 1) The zero matrix (force the pseudospectrum width to 1)
    // 2) Typical matrices (use a small multiple of the spectral radius)
    // 3) Highly non-normal matrices (e.g., triangular with zero main diagonal)
    Real width;
    if( oneNorm == Real(0) )
    {
        width = 1;
        if( psCtrl.progress && g.Rank() == 0 )
            Output("Setting width to 1 to handle zero matrix");
    }
    else if( radius >= 0.2*oneNorm )
    {
        width = 2.5*radius;
        if( psCtrl.progress && g.Rank() == 0 )
            Output
            ("Setting width to ",width,
             " based on the spectral radius, ",radius);
    }
    else
    {
        width = 0.8*oneNorm;
        if( psCtrl.progress && g.Rank() == 0 )
            Output
            ("Setting width to ",width," based on the one norm, ",oneNorm);
    }
    Complex<Real> center(0,0);

    box.center = center;
    box.realWidth = width;
    box.imagWidth = width;

    return TriangularSpectralWindow
           ( U, invNormMap, center, width, width, realSize, imagSize, psCtrl );
}

template<typename Field>
DistMatrix<Int> TriangularSpectralPortrait
( const AbstractDistMatrix<Field>& UPre,
  const AbstractDistMatrix<Field>& Q,
        AbstractDistMatrix<Base<Field>>& invNormMap,
  Int realSize,
  Int imagSize,
  SpectralBox<Base<Field>>& box,
  PseudospecCtrl<Base<Field>> psCtrl )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Grid& g = UPre.Grid();

    // Force 'U' to be in a [MC,MR] distribution so that we can get its diagonal
    DistMatrixReadProxy<Field,Field,MC,MR> UProx( UPre );
    auto& U = UProx.GetLocked();

    const Real radius = MaxNorm( GetDiagonal(U) );
    const Real oneNorm = OneNorm( U );

    // Essentially three cases are handled here:
    // 1) The zero matrix (force the pseudospectrum width to 1)
    // 2) Typical matrices (use a small multiple of the spectral radius)
    // 3) Highly non-normal matrices (e.g., triangular with zero main diagonal)
    Real width;
    if( oneNorm == Real(0) )
    {
        width = 1;
        if( psCtrl.progress && g.Rank() == 0 )
            Output("Setting width to 1 to handle zero matrix");
    }
    else if( radius >= 0.2*oneNorm )
    {
        width = 2.5*radius;
        if( psCtrl.progress && g.Rank() == 0 )
            Output
            ("Setting width to ",width,
             " based on the spectral radius, ",radius);
    }
    else
    {
        width = 0.8*oneNorm;
        if( psCtrl.progress && g.Rank() == 0 )
            Output
            ("Setting width to ",width," based on the one norm, ",oneNorm);
    }
    Complex<Real> center(0,0);

    box.center = center;
    box.realWidth = width;
    box.imagWidth = width;

    return TriangularSpectralWindow
           ( U, Q, invNormMap, center, width, width, realSize, imagSize,
             psCtrl );
}

template<typename Real>
DistMatrix<Int> QuasiTriangularSpectralPortrait
( const AbstractDistMatrix<Real>& UPre,
        AbstractDistMatrix<Real>& invNormMap,
  Int realSize,
  Int imagSize,
  SpectralBox<Real>& box,
  PseudospecCtrl<Real> psCtrl )
{
    EL_DEBUG_CSE
    const Grid& g = UPre.Grid();

    // Force 'U' to be in a [MC,MR] distribution to get its eigenvalues
    DistMatrixReadProxy<Real,Real,MC,MR> UProx( UPre );
    auto& U = UProx.GetLocked();
    const auto w = schur::QuasiTriangEig( U );

    const Real radius = MaxNorm( w );
    const Real oneNorm = OneNorm( U );

    // Essentially three cases are handled here:
    // 1) The zero matrix (force the pseudospectrum width to 1)
    // 2) Typical matrices (use a small multiple of the spectral radius)
    // 3) Highly non-normal matrices (e.g., triangular with zero main diagonal)
    Real width;
    if( oneNorm == Real(0) )
    {
        width = 1;
        if( psCtrl.progress && g.Rank() == 0 )
            Output("Setting width to 1 to handle zero matrix");
    }
    else if( radius >= 0.2*oneNorm )
    {
        width = 2.5*radius;
        if( psCtrl.progress && g.Rank() == 0 )
            Output
            ("Setting width to ",width,
             " based on the spectral radius, ",radius);
    }
    else
    {
        width = 0.8*oneNorm;
        if( psCtrl.progress && g.Rank() == 0 )
            Output
            ("Setting width to ",width," based on the one norm, ",oneNorm);
    }
    Complex<Real> center(0,0);

    box.center = center;
    box.realWidth = width;
    box.imagWidth = width;

    return QuasiTriangularSpectralWindow
           ( U, invNormMap, center, width, width, realSize, imagSize, psCtrl );
}

template<typename Real>
DistMatrix<Int> QuasiTriangularSpectralPortrait
( const AbstractDistMatrix<Real>& UPre,
  const AbstractDistMatrix<Real>& Q,
        AbstractDistMatrix<Real>& invNormMap,
  Int realSize,
  Int imagSize,
  SpectralBox<Real>& box,
  PseudospecCtrl<Real> psCtrl )
{
    EL_DEBUG_CSE
    const Grid& g = UPre.Grid();

    // Force 'U' to be in a [MC,MR] distribution to get its eigenvalues
    DistMatrixReadProxy<Real,Real,MC,MR> UProx( UPre );
    auto& U = UProx.GetLocked();
    const auto w = schur::QuasiTriangEig( U );

    const Real radius = MaxNorm( w );
    const Real oneNorm = OneNorm( U );

    // Essentially three cases are handled here:
    // 1) The zero matrix (force the pseudospectrum width to 1)
    // 2) Typical matrices (use a small multiple of the spectral radius)
    // 3) Highly non-normal matrices (e.g., triangular with zero main diagonal)
    Real width;
    if( oneNorm == Real(0) )
    {
        width = 1;
        if( psCtrl.progress && g.Rank() == 0 )
            Output("Setting width to 1 to handle zero matrix");
    }
    else if( radius >= 0.2*oneNorm )
    {
        width = 2.5*radius;
        if( psCtrl.progress && g.Rank() == 0 )
            Output
            ("Setting width to ",width,
             " based on the spectral radius, ",radius);
    }
    else
    {
        width = 0.8*oneNorm;
        if( psCtrl.progress && g.Rank() == 0 )
            Output
            ("Setting width to ",width," based on the one norm, ",oneNorm);
    }
    Complex<Real> center(0,0);

    box.center = center;
    box.realWidth = width;
    box.imagWidth = width;

    return QuasiTriangularSpectralWindow
           ( U, Q, invNormMap, center, width, width, realSize, imagSize,
             psCtrl );
}

template<typename Field>
DistMatrix<Int> HessenbergSpectralPortrait
( const AbstractDistMatrix<Field>& H,
        AbstractDistMatrix<Base<Field>>& invNormMap,
  Int realSize,
  Int imagSize,
  SpectralBox<Base<Field>>& box,
  PseudospecCtrl<Base<Field>> psCtrl )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Real oneNorm = OneNorm( H );
    const Real infNorm = InfinityNorm( H );

    Real width;
    if( oneNorm == Real(0) )
    {
        width = 1;
        if( psCtrl.progress && H.Grid().Rank() == 0 )
            Output("Setting width to 1 to handle zero matrix");
    }
    else
    {
        width = 0.8*Max(oneNorm,infNorm);
        if( psCtrl.progress && H.Grid().Rank() == 0 )
            Output
            ("Setting width to ",width," based on the one norm, ",oneNorm,
             ", and infinity norm, ",infNorm);
    }
    Complex<Real> center(0,0);

    box.center = center;
    box.realWidth = width;
    box.imagWidth = width;

    return HessenbergSpectralWindow
           ( H, invNormMap, center, width, width, realSize, imagSize, psCtrl );
}

template<typename Field>
DistMatrix<Int> HessenbergSpectralPortrait
( const AbstractDistMatrix<Field>& H,
  const AbstractDistMatrix<Field>& Q,
        AbstractDistMatrix<Base<Field>>& invNormMap,
  Int realSize,
  Int imagSize,
  SpectralBox<Base<Field>>& box,
  PseudospecCtrl<Base<Field>> psCtrl )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Real oneNorm = OneNorm( H );
    const Real infNorm = InfinityNorm( H );

    Real width;
    if( oneNorm == Real(0) )
    {
        width = 1;
        if( psCtrl.progress && H.Grid().Rank() == 0 )
            Output("Setting width to 1 to handle zero matrix");
    }
    else
    {
        width = 0.8*Max(oneNorm,infNorm);
        if( psCtrl.progress && H.Grid().Rank() == 0 )
            Output
            ("Setting width to ",width," based on the one norm, ",oneNorm,
             ", and infinity norm, ",infNorm);
    }
    Complex<Real> center(0,0);

    box.center = center;
    box.realWidth = width;
    box.imagWidth = width;

    return HessenbergSpectralWindow
           ( H, Q, invNormMap, center, width, width, realSize, imagSize,
             psCtrl );
}

namespace pspec {

template<typename Real>
Matrix<Int> Helper
( const Matrix<Complex<Real>>& A,
        Matrix<Real>& invNormMap,
  Int realSize,
  Int imagSize,
  SpectralBox<Real>& box,
  PseudospecCtrl<Real> psCtrl )
{
    EL_DEBUG_CSE
    typedef Complex<Real> C;

    Matrix<C> B( A );
    if( psCtrl.norm == PS_TWO_NORM )
    {
        if( psCtrl.schur )
        {
            Matrix<C> w;
            auto schurCtrl( psCtrl.schurCtrl );
            schurCtrl.hessSchurCtrl.fullTriangle = true;
            Schur( B, w, schurCtrl );
            return TriangularSpectralPortrait
                   ( B, invNormMap, realSize, imagSize, box, psCtrl );
        }
        else
        {
            hessenberg::ExplicitCondensed( UPPER, B );
            return HessenbergSpectralPortrait
                   ( B, invNormMap, realSize, imagSize, box, psCtrl );
        }
    }
    else
    {
        Matrix<C> Q;
        if( psCtrl.schur )
        {
            Matrix<C> w;
            auto schurCtrl( psCtrl.schurCtrl );
            schurCtrl.hessSchurCtrl.fullTriangle = true;
            Schur( B, w, Q, schurCtrl );
            return TriangularSpectralPortrait
                   ( B, Q, invNormMap, realSize, imagSize, box, psCtrl );
        }
        else
        {
            Matrix<C> t;
            Hessenberg( UPPER, B, t );
            Identity( Q, B.Height(), B.Height() );
            hessenberg::ApplyQ( LEFT, UPPER, NORMAL, B, t, Q );
            return HessenbergSpectralPortrait
                   ( B, Q, invNormMap, realSize, imagSize, box, psCtrl );
        }
    }
}

template<typename Real>
DistMatrix<Int> Helper
( const AbstractDistMatrix<Complex<Real>>& A,
        AbstractDistMatrix<Real>& invNormMap,
  Int realSize,
  Int imagSize,
  SpectralBox<Real>& box,
  PseudospecCtrl<Real> psCtrl )
{
    EL_DEBUG_CSE
    typedef Complex<Real> C;
    const Grid& g = A.Grid();
    DistMatrix<C> B( A );

    if( psCtrl.norm == PS_TWO_NORM )
    {
        if( psCtrl.schur )
        {
            DistMatrix<C,VR,STAR> w(g);
            auto schurCtrl( psCtrl.schurCtrl );
            schurCtrl.hessSchurCtrl.fullTriangle = true;
            Schur( B, w, schurCtrl );
            return TriangularSpectralPortrait
                   ( B, invNormMap, realSize, imagSize, box, psCtrl );
        }
        else
        {
            hessenberg::ExplicitCondensed( UPPER, B );
            return HessenbergSpectralPortrait
                   ( B, invNormMap, realSize, imagSize, box, psCtrl );
        }
    }
    else
    {
        DistMatrix<C> Q(g);
        if( psCtrl.schur )
        {
            DistMatrix<C,VR,STAR> w(g);
            auto schurCtrl( psCtrl.schurCtrl );
            schurCtrl.hessSchurCtrl.fullTriangle = true;
            Schur( B, w, Q, schurCtrl );
            return TriangularSpectralPortrait
                   ( B, Q, invNormMap, realSize, imagSize, box, psCtrl );
        }
        else
        {
            DistMatrix<C,STAR,STAR> t(g);
            Hessenberg( UPPER, B, t );
            Identity( Q, B.Height(), B.Height() );
            hessenberg::ApplyQ( LEFT, UPPER, NORMAL, B, t, Q );
            return HessenbergSpectralPortrait
                   ( B, Q, invNormMap, realSize, imagSize, box, psCtrl );
        }
    }
}

template<typename Real>
Matrix<Int> Helper
( const Matrix<Real>& A,
        Matrix<Real>& invNormMap,
  Int realSize,
  Int imagSize,
  SpectralBox<Real>& box,
  PseudospecCtrl<Real> psCtrl )
{
    EL_DEBUG_CSE
    typedef Complex<Real> C;

    if( psCtrl.forceComplexSchur )
    {
        Matrix<C> ACpx;
        Copy( A, ACpx );
        return Helper( ACpx, invNormMap, realSize, imagSize, box, psCtrl );
    }

    if( !psCtrl.schur )
        LogicError("Real Hessenberg algorithm not yet supported");
    Matrix<Real> B( A );
    Matrix<C> w;
    auto schurCtrl( psCtrl.schurCtrl );
    schurCtrl.hessSchurCtrl.fullTriangle = true;
    if( psCtrl.norm == PS_TWO_NORM )
    {
        Schur( B, w, schurCtrl );
        if( psCtrl.forceComplexPs )
        {
            Matrix<C> BCpx;
            schur::RealToComplex( B, BCpx );
            return TriangularSpectralPortrait
                   ( BCpx, invNormMap, realSize, imagSize, box, psCtrl );
        }
        return QuasiTriangularSpectralPortrait
               ( B, invNormMap, realSize, imagSize, box, psCtrl );
    }
    else
    {
        Matrix<Real> Q;
        Schur( B, w, Q, schurCtrl );
        if( psCtrl.forceComplexPs )
        {
            LogicError("Real to complex full Schur not yet supported");
            Matrix<C> BCpx, QCpx;
            schur::RealToComplex( B, Q, BCpx, QCpx );
            return TriangularSpectralPortrait
                   ( BCpx, invNormMap, realSize, imagSize, box, psCtrl );
        }
        return QuasiTriangularSpectralPortrait
               ( B, Q, invNormMap, realSize, imagSize, box, psCtrl );
    }
}

template<typename Real>
DistMatrix<Int> Helper
( const AbstractDistMatrix<Real>& A,
        AbstractDistMatrix<Real>& invNormMap,
  Int realSize,
  Int imagSize,
  SpectralBox<Real>& box,
  PseudospecCtrl<Real> psCtrl )
{
    EL_DEBUG_CSE
    typedef Complex<Real> C;
    const Grid& g = A.Grid();

    if( psCtrl.forceComplexSchur )
    {
        DistMatrix<C> ACpx(g);
        Copy( A, ACpx );
        return SpectralPortrait
               ( ACpx, invNormMap, realSize, imagSize, box, psCtrl );
    }

    if( !psCtrl.schur )
        LogicError("Real Hessenberg algorithm not yet supported");
    DistMatrix<Real> B( A );

    DistMatrix<C,VR,STAR> w(g);
    auto schurCtrl( psCtrl.schurCtrl );
    schurCtrl.hessSchurCtrl.fullTriangle = true;
    if( psCtrl.norm == PS_TWO_NORM )
    {
        Schur( B, w, schurCtrl );
        if( psCtrl.forceComplexPs )
        {
            DistMatrix<C> BCpx(g);
            schur::RealToComplex( B, BCpx );
            return TriangularSpectralPortrait
                   ( BCpx, invNormMap, realSize, imagSize, box, psCtrl );
        }
        return QuasiTriangularSpectralPortrait
               ( B, invNormMap, realSize, imagSize, box, psCtrl );
    }
    else
    {
        DistMatrix<Real> Q(g);
        Schur( B, w, Q, schurCtrl );
        if( psCtrl.forceComplexPs )
        {
            DistMatrix<C> BCpx(g), QCpx(g);
            schur::RealToComplex( B, Q, BCpx, QCpx );
            return TriangularSpectralPortrait
                   ( BCpx, QCpx, invNormMap, realSize, imagSize, box, psCtrl );
        }
        return QuasiTriangularSpectralPortrait
               ( B, Q, invNormMap, realSize, imagSize, box, psCtrl );
    }
}

} // namespace pspec

template<typename Field>
Matrix<Int> SpectralPortrait
( const Matrix<Field>& A,
        Matrix<Base<Field>>& invNormMap,
  Int realSize,
  Int imagSize,
  SpectralBox<Base<Field>>& box,
  PseudospecCtrl<Base<Field>> psCtrl )
{
    EL_DEBUG_CSE
    return pspec::Helper( A, invNormMap, realSize, imagSize, box, psCtrl );
}

template<typename Field>
DistMatrix<Int> SpectralPortrait
( const AbstractDistMatrix<Field>& A,
        AbstractDistMatrix<Base<Field>>& invNormMap,
  Int realSize,
  Int imagSize,
  SpectralBox<Base<Field>>& box,
  PseudospecCtrl<Base<Field>> psCtrl )
{
    EL_DEBUG_CSE
    return pspec::Helper( A, invNormMap, realSize, imagSize, box, psCtrl );
}

#define PROTO(Field) \
  template Matrix<Int> SpectralCloud \
  ( const Matrix<Field>& A, \
    const Matrix<Complex<Base<Field>>>& shifts, \
          Matrix<Base<Field>>& invNorms, \
          PseudospecCtrl<Base<Field>> psCtrl ); \
  template DistMatrix<Int,VR,STAR> SpectralCloud \
  ( const AbstractDistMatrix<Field>& A, \
    const AbstractDistMatrix<Complex<Base<Field>>>& shifts, \
          AbstractDistMatrix<Base<Field>>& invNorms, \
          PseudospecCtrl<Base<Field>> psCtrl ); \
  template Matrix<Int> TriangularSpectralCloud \
  ( const Matrix<Field>& U, \
    const Matrix<Complex<Base<Field>>>& shifts, \
          Matrix<Base<Field>>& invNorms, \
          PseudospecCtrl<Base<Field>> psCtrl ); \
  template DistMatrix<Int,VR,STAR> TriangularSpectralCloud \
  ( const AbstractDistMatrix<Field>& U, \
    const AbstractDistMatrix<Complex<Base<Field>>>& shifts, \
          AbstractDistMatrix<Base<Field>>& invNorms, \
          PseudospecCtrl<Base<Field>> psCtrl ); \
  template Matrix<Int> TriangularSpectralCloud \
  ( const Matrix<Field>& U, \
    const Matrix<Field>& Q, \
    const Matrix<Complex<Base<Field>>>& shifts, \
          Matrix<Base<Field>>& invNorms, \
          PseudospecCtrl<Base<Field>> psCtrl ); \
  template DistMatrix<Int,VR,STAR> TriangularSpectralCloud \
  ( const AbstractDistMatrix<Field>& U, \
    const AbstractDistMatrix<Field>& Q, \
    const AbstractDistMatrix<Complex<Base<Field>>>& shifts, \
          AbstractDistMatrix<Base<Field>>& invNorms, \
          PseudospecCtrl<Base<Field>> psCtrl ); \
  template Matrix<Int> HessenbergSpectralCloud \
  ( const Matrix<Field>& H, \
    const Matrix<Complex<Base<Field>>>& shifts, \
          Matrix<Base<Field>>& invNorms, \
          PseudospecCtrl<Base<Field>> psCtrl ); \
  template DistMatrix<Int,VR,STAR> HessenbergSpectralCloud \
  ( const AbstractDistMatrix<Field>& H, \
    const AbstractDistMatrix<Complex<Base<Field>>>& shifts, \
          AbstractDistMatrix<Base<Field>>& invNorms, \
          PseudospecCtrl<Base<Field>> psCtrl ); \
  template Matrix<Int> HessenbergSpectralCloud \
  ( const Matrix<Field>& H, \
    const Matrix<Field>& Q, \
    const Matrix<Complex<Base<Field>>>& shifts, \
          Matrix<Base<Field>>& invNorms, \
          PseudospecCtrl<Base<Field>> psCtrl ); \
  template DistMatrix<Int,VR,STAR> HessenbergSpectralCloud \
  ( const AbstractDistMatrix<Field>& H, \
    const AbstractDistMatrix<Field>& Q, \
    const AbstractDistMatrix<Complex<Base<Field>>>& shifts, \
          AbstractDistMatrix<Base<Field>>& invNorms, \
          PseudospecCtrl<Base<Field>> psCtrl ); \
  template Matrix<Int> TriangularSpectralWindow \
  ( const Matrix<Field>& U, \
          Matrix<Base<Field>>& invNormMap, \
          Complex<Base<Field>> center, \
          Base<Field> realWidth, \
          Base<Field> imagWidth, \
          Int realSize, \
          Int imagSize, \
          PseudospecCtrl<Base<Field>> psCtrl ); \
  template DistMatrix<Int> TriangularSpectralWindow \
  ( const AbstractDistMatrix<Field>& U, \
          AbstractDistMatrix<Base<Field>>& invNormMap, \
          Complex<Base<Field>> center, \
          Base<Field> realWidth, \
          Base<Field> imagWidth, \
          Int realSize, \
          Int imagSize, \
          PseudospecCtrl<Base<Field>> psCtrl ); \
  template Matrix<Int> TriangularSpectralWindow \
  ( const Matrix<Field>& U, \
    const Matrix<Field>& Q, \
          Matrix<Base<Field>>& invNormMap, \
          Complex<Base<Field>> center, \
          Base<Field> realWidth, \
          Base<Field> imagWidth, \
          Int realSize, \
          Int imagSize, \
          PseudospecCtrl<Base<Field>> psCtrl ); \
  template DistMatrix<Int> TriangularSpectralWindow \
  ( const AbstractDistMatrix<Field>& U, \
    const AbstractDistMatrix<Field>& Q, \
          AbstractDistMatrix<Base<Field>>& invNormMap, \
          Complex<Base<Field>> center, \
          Base<Field> realWidth, \
          Base<Field> imagWidth, \
          Int realSize, \
          Int imagSize, \
          PseudospecCtrl<Base<Field>> psCtrl ); \
  template Matrix<Int> HessenbergSpectralWindow \
  ( const Matrix<Field>& H, \
          Matrix<Base<Field>>& invNormMap, \
          Complex<Base<Field>> center, \
          Base<Field> realWidth, \
          Base<Field> imagWidth, \
          Int realSize, \
          Int imagSize, \
          PseudospecCtrl<Base<Field>> psCtrl ); \
  template DistMatrix<Int> HessenbergSpectralWindow \
  ( const AbstractDistMatrix<Field>& H, \
          AbstractDistMatrix<Base<Field>>& invNormMap, \
          Complex<Base<Field>> center, \
          Base<Field> realWidth, \
          Base<Field> imagWidth, \
          Int realSize, \
          Int imagSize, \
          PseudospecCtrl<Base<Field>> psCtrl ); \
  template Matrix<Int> HessenbergSpectralWindow \
  ( const Matrix<Field>& H, \
    const Matrix<Field>& Q, \
          Matrix<Base<Field>>& invNormMap, \
          Complex<Base<Field>> center, \
          Base<Field> realWidth, \
          Base<Field> imagWidth, \
          Int realSize, \
          Int imagSize, \
          PseudospecCtrl<Base<Field>> psCtrl ); \
  template DistMatrix<Int> HessenbergSpectralWindow \
  ( const AbstractDistMatrix<Field>& H, \
    const AbstractDistMatrix<Field>& Q, \
          AbstractDistMatrix<Base<Field>>& invNormMap, \
    Complex<Base<Field>> center, \
    Base<Field> realWidth, \
    Base<Field> imagWidth, \
    Int realSize, Int imagSize, PseudospecCtrl<Base<Field>> psCtrl ); \
  template Matrix<Int> SpectralWindow \
  ( const Matrix<Field>& A, Matrix<Base<Field>>& invNormMap, \
    Complex<Base<Field>> center, \
    Base<Field> realWidth, \
    Base<Field> imagWidth, \
    Int realSize, Int imagSize, PseudospecCtrl<Base<Field>> psCtrl ); \
  template DistMatrix<Int> SpectralWindow \
  ( const AbstractDistMatrix<Field>& A, \
          AbstractDistMatrix<Base<Field>>& invNormMap, \
    Complex<Base<Field>> center, Base<Field> realWidth, Base<Field> imagWidth, \
    Int realSize, Int imagSize, PseudospecCtrl<Base<Field>> psCtrl ); \
  template Matrix<Int> SpectralPortrait \
  ( const Matrix<Field>& A, Matrix<Base<Field>>& invNormMap, \
    Int realSize, Int imagSize, SpectralBox<Base<Field>>& box, \
    PseudospecCtrl<Base<Field>> psCtrl ); \
  template DistMatrix<Int> SpectralPortrait \
  ( const AbstractDistMatrix<Field>& A, \
          AbstractDistMatrix<Base<Field>>& invNormMap, \
    Int realSize, Int imagSize, SpectralBox<Base<Field>>& box, \
    PseudospecCtrl<Base<Field>> psCtrl ); \
  template Matrix<Int> TriangularSpectralPortrait \
  ( const Matrix<Field>& U, Matrix<Base<Field>>& invNormMap, \
    Int realSize, Int imagSize, SpectralBox<Base<Field>>& box, \
    PseudospecCtrl<Base<Field>> psCtrl ); \
  template DistMatrix<Int> TriangularSpectralPortrait \
  ( const AbstractDistMatrix<Field>& U, AbstractDistMatrix<Base<Field>>& invNormMap, \
    Int realSize, Int imagSize, SpectralBox<Base<Field>>& box, \
    PseudospecCtrl<Base<Field>> psCtrl ); \
  template Matrix<Int> TriangularSpectralPortrait \
  ( const Matrix<Field>& U, \
    const Matrix<Field>& Q, \
          Matrix<Base<Field>>& invNormMap, \
    Int realSize, Int imagSize, SpectralBox<Base<Field>>& box, \
    PseudospecCtrl<Base<Field>> psCtrl ); \
  template DistMatrix<Int> TriangularSpectralPortrait \
  ( const AbstractDistMatrix<Field>& U, const AbstractDistMatrix<Field>& Q, \
    AbstractDistMatrix<Base<Field>>& invNormMap, \
    Int realSize, Int imagSize, SpectralBox<Base<Field>>& box, \
    PseudospecCtrl<Base<Field>> psCtrl ); \
  template Matrix<Int> HessenbergSpectralPortrait \
  ( const Matrix<Field>& H, Matrix<Base<Field>>& invNormMap, \
    Int realSize, Int imagSize, SpectralBox<Base<Field>>& box, \
    PseudospecCtrl<Base<Field>> psCtrl ); \
  template DistMatrix<Int> HessenbergSpectralPortrait \
  ( const AbstractDistMatrix<Field>& H, AbstractDistMatrix<Base<Field>>& invNormMap, \
    Int realSize, Int imagSize, SpectralBox<Base<Field>>& box, \
    PseudospecCtrl<Base<Field>> psCtrl ); \
  template Matrix<Int> HessenbergSpectralPortrait \
  ( const Matrix<Field>& H, \
    const Matrix<Field>& Q, \
          Matrix<Base<Field>>& invNormMap, \
    Int realSize, Int imagSize, SpectralBox<Base<Field>>& box, \
    PseudospecCtrl<Base<Field>> psCtrl ); \
  template DistMatrix<Int> HessenbergSpectralPortrait \
  ( const AbstractDistMatrix<Field>& H, \
    const AbstractDistMatrix<Field>& Q, \
          AbstractDistMatrix<Base<Field>>& invNormMap, \
          Int realSize, \
          Int imagSize, \
          SpectralBox<Base<Field>>& box, \
          PseudospecCtrl<Base<Field>> psCtrl );

#define PROTO_REAL(Real) \
  PROTO(Real) \
  template Matrix<Int> QuasiTriangularSpectralCloud \
  ( const Matrix<Real>& U, \
    const Matrix<Complex<Real>>& shifts, \
          Matrix<Real>& invNorms, \
          PseudospecCtrl<Real> psCtrl ); \
  template DistMatrix<Int,VR,STAR> QuasiTriangularSpectralCloud \
  ( const AbstractDistMatrix<Real>& U, \
    const AbstractDistMatrix<Complex<Real>>& shifts, \
          AbstractDistMatrix<Real>& invNorms, \
          PseudospecCtrl<Real> psCtrl ); \
  template Matrix<Int> QuasiTriangularSpectralCloud \
  ( const Matrix<Real>& U, \
    const Matrix<Real>& Q, \
    const Matrix<Complex<Real>>& shifts, \
          Matrix<Real>& invNorms, \
          PseudospecCtrl<Real> psCtrl ); \
  template DistMatrix<Int,VR,STAR> QuasiTriangularSpectralCloud \
  ( const AbstractDistMatrix<Real>& U, \
    const AbstractDistMatrix<Real>& Q, \
    const AbstractDistMatrix<Complex<Real>>& shifts, \
          AbstractDistMatrix<Real>& invNorms, \
          PseudospecCtrl<Real> psCtrl ); \
  template Matrix<Int> QuasiTriangularSpectralWindow \
  ( const Matrix<Real>& U, \
          Matrix<Real>& invNormMap, \
          Complex<Real> center, \
          Real realWidth, \
          Real imagWidth, \
          Int realSize, \
          Int imagSize, \
          PseudospecCtrl<Real> psCtrl ); \
  template DistMatrix<Int> QuasiTriangularSpectralWindow \
  ( const AbstractDistMatrix<Real>& U, \
          AbstractDistMatrix<Real>& invNormMap, \
          Complex<Real> center, \
          Real realWidth, \
          Real imagWidth, \
          Int realSize, \
          Int imagSize, \
          PseudospecCtrl<Real> psCtrl ); \
  template Matrix<Int> QuasiTriangularSpectralWindow \
  ( const Matrix<Real>& U, \
    const Matrix<Real>& Q, \
          Matrix<Real>& invNormMap, \
          Complex<Real> center, \
          Real realWidth, \
          Real imagWidth, \
          Int realSize, \
          Int imagSize, \
          PseudospecCtrl<Real> psCtrl ); \
  template DistMatrix<Int> QuasiTriangularSpectralWindow \
  ( const AbstractDistMatrix<Real>& U, \
    const AbstractDistMatrix<Real>& Q, \
          AbstractDistMatrix<Real>& invNormMap, \
          Complex<Real> center, \
          Real realWidth, \
          Real imagWidth, \
          Int realSize, \
          Int imagSize, \
          PseudospecCtrl<Real> psCtrl ); \
  template Matrix<Int> QuasiTriangularSpectralPortrait \
  ( const Matrix<Real>& U, \
          Matrix<Real>& invNormMap, \
          Int realSize, \
          Int imagSize, \
          SpectralBox<Real>& box, \
          PseudospecCtrl<Real> psCtrl ); \
  template DistMatrix<Int> QuasiTriangularSpectralPortrait \
  ( const AbstractDistMatrix<Real>& U, \
          AbstractDistMatrix<Real>& invNormMap, \
          Int realSize, \
          Int imagSize, \
          SpectralBox<Real>& box, \
          PseudospecCtrl<Real> psCtrl ); \
  template Matrix<Int> QuasiTriangularSpectralPortrait \
  ( const Matrix<Real>& U, \
    const Matrix<Real>& Q, \
          Matrix<Real>& invNormMap, \
          Int realSize, \
          Int imagSize, \
          SpectralBox<Real>& box, \
          PseudospecCtrl<Real> psCtrl ); \
  template DistMatrix<Int> QuasiTriangularSpectralPortrait \
  ( const AbstractDistMatrix<Real>& U, \
    const AbstractDistMatrix<Real>& Q, \
          AbstractDistMatrix<Real>& invNormMap, \
          Int realSize, \
          Int imagSize, \
          SpectralBox<Real>& box, \
          PseudospecCtrl<Real> psCtrl );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
