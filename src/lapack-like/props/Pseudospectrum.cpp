/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

#include EL_IDENTITY_INC

#include "./Pseudospectrum/Util.hpp"
#include "./Pseudospectrum/Power.hpp"
#include "./Pseudospectrum/Lanczos.hpp"
#include "./Pseudospectrum/IRA.hpp"
#include "./Pseudospectrum/IRL.hpp"
#include "./Pseudospectrum/Analytic.hpp"

// For one-norm pseudospectra. An adaptation of the more robust algorithm of
// Higham and Tisseur will hopefully be implemented soon.
#include "./Pseudospectrum/HagerHigham.hpp"

namespace El {

template<typename Real>
Matrix<Int> TriangularPseudospectrum
( const Matrix<Complex<Real>>& U, const Matrix<Complex<Real>>& shifts, 
  Matrix<Real>& invNorms, PseudospecCtrl<Real> psCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("TriangularPseudospectrum"))

    // Check if the off-diagonal is sufficiently small; if so, compute the 
    // pseudospectrum analytically from the eigenvalues. This also takes care
    // of the case where the matrix is a constant multiple of the identity 
    // matrix, which, after shifting, can lead to the zero matrix, which would 
    // cause problems for the Lanczos convergence criteria.
    if( pspec::TriangIsNormal( U, psCtrl.tol ) )
    {
        Matrix<Int> itCounts;
        if( psCtrl.progress )
            std::cout << "Matrix was numerically normal" << std::endl;
        auto w = U.GetDiagonal();
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

template<typename Real>
Matrix<Int> TriangularPseudospectrum
( const Matrix<Complex<Real>>& U, const Matrix<Complex<Real>>& Q, 
  const Matrix<Complex<Real>>& shifts, 
  Matrix<Real>& invNorms, PseudospecCtrl<Real> psCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("TriangularPseudospectrum"))

    // Check if the off-diagonal is sufficiently small; if so, compute the 
    // pseudospectrum analytically from the eigenvalues. This also takes care
    // of the case where the matrix is a constant multiple of the identity 
    // matrix, which, after shifting, can lead to the zero matrix, which would 
    // cause problems for the Lanczos convergence criteria.
    if( pspec::TriangIsNormal( U, psCtrl.tol ) )
    {
        Matrix<Int> itCounts;
        if( psCtrl.progress )
            std::cout << "Matrix was numerically normal" << std::endl;
        auto w = U.GetDiagonal();
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
        return pspec::HagerHigham( U, Q, shifts, invNorms, psCtrl );
}

template<typename Real>
Matrix<Int> TriangularPseudospectrum
( const Matrix<Real>& U, const Matrix<Complex<Real>>& shifts, 
  Matrix<Real>& invNorms, PseudospecCtrl<Real> psCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("TriangularPseudospectrum"))
    Matrix<Complex<Real>> UCpx;
    Copy( U, UCpx );

    // TODO: Use a real multi-shift TRSM instead
    return TriangularPseudospectrum( UCpx, shifts, invNorms, psCtrl );
}

template<typename Real>
Matrix<Int> TriangularPseudospectrum
( const Matrix<Real>& U, const Matrix<Real>& Q, 
  const Matrix<Complex<Real>>& shifts, Matrix<Real>& invNorms, 
  PseudospecCtrl<Real> psCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("TriangularPseudospectrum"))
    Matrix<Complex<Real>> UCpx, QCpx;
    Copy( U, UCpx );
    Copy( Q, QCpx );

    // TODO: Use a real multi-shift TRSM instead
    return TriangularPseudospectrum( UCpx, QCpx, shifts, invNorms, psCtrl );
}

template<typename Real>
Matrix<Int> QuasiTriangularPseudospectrum
( const Matrix<Real>& U, 
  const Matrix<Complex<Real>>& shifts, 
  Matrix<Real>& invNorms, PseudospecCtrl<Real> psCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("QuasiTriangularPseudospectrum"))

    // Check if the off-diagonal is sufficiently small; if so, compute the 
    // pseudospectrum analytically from the eigenvalues. This also takes care
    // of the case where the matrix is a constant multiple of the identity 
    // matrix, which, after shifting, can lead to the zero matrix, which would 
    // cause problems for the Lanczos convergence criteria.
    if( pspec::QuasiTriangIsNormal( U, psCtrl.tol ) )
    {
        Matrix<Int> itCounts;
        if( psCtrl.progress )
            std::cout << "Matrix was numerically normal" << std::endl;
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
Matrix<Int> QuasiTriangularPseudospectrum
( const Matrix<Real>& U,  const Matrix<Real>& Q,
  const Matrix<Complex<Real>>& shifts,
  Matrix<Real>& invNorms, PseudospecCtrl<Real> psCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("QuasiTriangularPseudospectrum"))

    // Check if the off-diagonal is sufficiently small; if so, compute the 
    // pseudospectrum analytically from the eigenvalues. This also takes care
    // of the case where the matrix is a constant multiple of the identity 
    // matrix, which, after shifting, can lead to the zero matrix, which would 
    // cause problems for the Lanczos convergence criteria.
    if( pspec::QuasiTriangIsNormal( U, psCtrl.tol ) )
    {
        Matrix<Int> itCounts;
        if( psCtrl.progress )
            std::cout << "Matrix was numerically normal" << std::endl;
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

template<typename Real>
Matrix<Int> HessenbergPseudospectrum
( const Matrix<Complex<Real>>& H, const Matrix<Complex<Real>>& shifts, 
  Matrix<Real>& invNorms, PseudospecCtrl<Real> psCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("HessenbergPseudospectrum"))

    // TODO: Check if the subdiagonal is numerically zero, and, if so, revert to
    //       triangular version of Pseudospectrum?
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

template<typename Real>
Matrix<Int> HessenbergPseudospectrum
( const Matrix<Complex<Real>>& H, 
  const Matrix<Complex<Real>>& Q, 
  const Matrix<Complex<Real>>& shifts, 
  Matrix<Real>& invNorms, PseudospecCtrl<Real> psCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("HessenbergPseudospectrum"))

    // TODO: Check if the subdiagonal is numerically zero, and, if so, revert to
    //       triangular version of Pseudospectrum?
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
        return pspec::HagerHigham( H, Q, shifts, invNorms, psCtrl );
}

template<typename Real>
Matrix<Int> HessenbergPseudospectrum
( const Matrix<Real>& H, const Matrix<Complex<Real>>& shifts, 
  Matrix<Real>& invNorms, PseudospecCtrl<Real> psCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("HessenbergPseudospectrum"))
    Matrix<Complex<Real>> HCpx;
    Copy( H, HCpx );

    // TODO: Use a real multi-shift Hess. solve instead
    return HessenbergPseudospectrum( HCpx, shifts, invNorms, psCtrl );
}

template<typename Real>
Matrix<Int> HessenbergPseudospectrum
( const Matrix<Real>& H, 
  const Matrix<Real>& Q, 
  const Matrix<Complex<Real>>& shifts, 
  Matrix<Real>& invNorms, PseudospecCtrl<Real> psCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("HessenbergPseudospectrum"))
    Matrix<Complex<Real>> HCpx, QCpx;
    Copy( H, HCpx );
    Copy( Q, QCpx );

    // TODO: Use a real multi-shift Hess. solve instead
    return HessenbergPseudospectrum( HCpx, QCpx, shifts, invNorms, psCtrl );
}

template<typename Real>
DistMatrix<Int,VR,STAR> TriangularPseudospectrum
( const DistMatrix<Complex<Real>>& U, 
  const DistMatrix<Complex<Real>,VR,STAR>& shifts,
  DistMatrix<Real,VR,STAR>& invNorms, PseudospecCtrl<Real> psCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("TriangularPseudospectrum"))
    const Grid& g = U.Grid();

    // Check if the off-diagonal is sufficiently small; if so, compute the 
    // pseudospectrum analytically from the eigenvalues. This also takes care
    // of the case where the matrix is a constant multiple of the identity 
    // matrix, which, after shifting, can lead to the zero matrix, which would 
    // cause problems for the Lanczos convergence criteria.
    if( pspec::TriangIsNormal( U, psCtrl.tol ) )
    {
        DistMatrix<Int,VR,STAR> itCounts(g);
        if( psCtrl.progress && g.Rank() == 0 )
            std::cout << "Matrix was numerically normal" << std::endl;
        auto w = U.GetDiagonal();
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

template<typename Real>
DistMatrix<Int,VR,STAR> TriangularPseudospectrum
( const DistMatrix<Complex<Real>>& U, 
  const DistMatrix<Complex<Real>>& Q, 
  const DistMatrix<Complex<Real>,VR,STAR>& shifts,
  DistMatrix<Real,VR,STAR>& invNorms, PseudospecCtrl<Real> psCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("TriangularPseudospectrum"))
    const Grid& g = U.Grid();

    // Check if the off-diagonal is sufficiently small; if so, compute the 
    // pseudospectrum analytically from the eigenvalues. This also takes care
    // of the case where the matrix is a constant multiple of the identity 
    // matrix, which, after shifting, can lead to the zero matrix, which would 
    // cause problems for the Lanczos convergence criteria.
    if( pspec::TriangIsNormal( U, psCtrl.tol ) )
    {
        DistMatrix<Int,VR,STAR> itCounts(g);
        if( psCtrl.progress && g.Rank() == 0 )
            std::cout << "Matrix was numerically normal" << std::endl;
        auto w = U.GetDiagonal();
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
        return pspec::HagerHigham( U, Q, shifts, invNorms, psCtrl );
}

template<typename Real>
DistMatrix<Int,VR,STAR> TriangularPseudospectrum
( const DistMatrix<Real>& U, const DistMatrix<Complex<Real>,VR,STAR>& shifts,
  DistMatrix<Real,VR,STAR>& invNorms, PseudospecCtrl<Real> psCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("TriangularPseudospectrum"))
    DistMatrix<Complex<Real>> UCpx(U.Grid());
    Copy( U, UCpx );

    // TODO: Use a real multi-shift TRSM instead
    return TriangularPseudospectrum( UCpx, shifts, invNorms, psCtrl );
}

template<typename Real>
DistMatrix<Int,VR,STAR> TriangularPseudospectrum
( const DistMatrix<Real>& U, 
  const DistMatrix<Real>& Q,
  const DistMatrix<Complex<Real>,VR,STAR>& shifts,
  DistMatrix<Real,VR,STAR>& invNorms, PseudospecCtrl<Real> psCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("TriangularPseudospectrum"))
    DistMatrix<Complex<Real>> UCpx(U.Grid()), QCpx(U.Grid());
    Copy( U, UCpx );
    Copy( Q, QCpx );

    // TODO: Use a real multi-shift TRSM instead
    return TriangularPseudospectrum( UCpx, QCpx, shifts, invNorms, psCtrl );
}

template<typename Real>
DistMatrix<Int,VR,STAR> QuasiTriangularPseudospectrum
( const DistMatrix<Real>& U, 
  const DistMatrix<Complex<Real>,VR,STAR>& shifts,
  DistMatrix<Real,VR,STAR>& invNorms, PseudospecCtrl<Real> psCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("QuasiTriangularPseudospectrum"))
    const Grid& g = U.Grid();

    // Check if the off-diagonal is sufficiently small; if so, compute the 
    // pseudospectrum analytically from the eigenvalues. This also takes care
    // of the case where the matrix is a constant multiple of the identity 
    // matrix, which, after shifting, can lead to the zero matrix, which would 
    // cause problems for the Lanczos convergence criteria.
    if( pspec::QuasiTriangIsNormal( U, psCtrl.tol ) )
    {
        DistMatrix<Int,VR,STAR> itCounts(g);
        if( psCtrl.progress && g.Rank() == 0 )
            std::cout << "Matrix was numerically normal" << std::endl;
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
DistMatrix<Int,VR,STAR> QuasiTriangularPseudospectrum
( const DistMatrix<Real>& U, 
  const DistMatrix<Real>& Q,
  const DistMatrix<Complex<Real>,VR,STAR>& shifts,
  DistMatrix<Real,VR,STAR>& invNorms, PseudospecCtrl<Real> psCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("QuasiTriangularPseudospectrum"))
    const Grid& g = U.Grid();

    // Check if the off-diagonal is sufficiently small; if so, compute the 
    // pseudospectrum analytically from the eigenvalues. This also takes care
    // of the case where the matrix is a constant multiple of the identity 
    // matrix, which, after shifting, can lead to the zero matrix, which would 
    // cause problems for the Lanczos convergence criteria.
    if( pspec::QuasiTriangIsNormal( U, psCtrl.tol ) )
    {
        DistMatrix<Int,VR,STAR> itCounts(g);
        if( psCtrl.progress && g.Rank() == 0 )
            std::cout << "Matrix was numerically normal" << std::endl;
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
DistMatrix<Int,VR,STAR> HessenbergPseudospectrum
( const DistMatrix<Complex<Real>>& H, 
  const DistMatrix<Complex<Real>,VR,STAR>& shifts,
  DistMatrix<Real,VR,STAR>& invNorms, PseudospecCtrl<Real> psCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("HessenbergPseudospectrum"))

    // TODO: Check if the subdiagonal is sufficiently small, and, if so, revert
    //       to TriangularPseudospectrum
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

template<typename Real>
DistMatrix<Int,VR,STAR> HessenbergPseudospectrum
( const DistMatrix<Complex<Real>>& H, 
  const DistMatrix<Complex<Real>>& Q,
  const DistMatrix<Complex<Real>,VR,STAR>& shifts,
  DistMatrix<Real,VR,STAR>& invNorms, PseudospecCtrl<Real> psCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("HessenbergPseudospectrum"))

    // TODO: Check if the subdiagonal is sufficiently small, and, if so, revert
    //       to TriangularPseudospectrum
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
        return pspec::HagerHigham( H, Q, shifts, invNorms, psCtrl );
}

template<typename Real>
DistMatrix<Int,VR,STAR> HessenbergPseudospectrum
( const DistMatrix<Real>& H, const DistMatrix<Complex<Real>,VR,STAR>& shifts,
  DistMatrix<Real,VR,STAR>& invNorms, PseudospecCtrl<Real> psCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("HessenbergPseudospectrum"))
    DistMatrix<Complex<Real>> HCpx(H.Grid());
    Copy( H, HCpx );

    return HessenbergPseudospectrum( HCpx, shifts, invNorms, psCtrl );
}

template<typename Real>
DistMatrix<Int,VR,STAR> HessenbergPseudospectrum
( const DistMatrix<Real>& H, 
  const DistMatrix<Real>& Q,
  const DistMatrix<Complex<Real>,VR,STAR>& shifts,
  DistMatrix<Real,VR,STAR>& invNorms, PseudospecCtrl<Real> psCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("HessenbergPseudospectrum"))
    DistMatrix<Complex<Real>> HCpx(H.Grid()), QCpx(H.Grid());
    Copy( H, HCpx );
    Copy( Q, QCpx );

    return HessenbergPseudospectrum( HCpx, QCpx, shifts, invNorms, psCtrl );
}

template<typename Real>
Matrix<Int> Pseudospectrum
( const Matrix<Real>& A, const Matrix<Complex<Real>>& shifts, 
  Matrix<Real>& invNorms, PseudospecCtrl<Real> psCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("Pseudospectrum"))
    typedef Complex<Real> C;

    if( psCtrl.forceComplexSchur )
    {
        Matrix<C> ACpx;
        Copy( A, ACpx );
        return Pseudospectrum( ACpx, shifts, invNorms, psCtrl );
    }

    if( !psCtrl.schur )
        LogicError("Real Hessenberg algorithm not yet supported");
    Matrix<Real> U( A );
    Matrix<C> w;
    const bool fullTriangle = true;
    if( psCtrl.norm == PS_TWO_NORM )
    {
        Schur( U, w, fullTriangle, psCtrl.schurCtrl );
        if( psCtrl.forceComplexPs )
        {
            Matrix<C> UCpx;
            schur::RealToComplex( U, UCpx );
            return TriangularPseudospectrum( UCpx, shifts, invNorms, psCtrl );
        }
        return QuasiTriangularPseudospectrum( U, shifts, invNorms, psCtrl );
    }
    else
    {
        Matrix<Real> Q;
        Schur( U, w, Q, fullTriangle, psCtrl.schurCtrl );
        if( psCtrl.forceComplexPs )
        {
            LogicError("Real to complex full Schur not yet supported");
            /*
            Matrix<C> UCpx, QCpx;
            schur::RealToComplex( U, Q, UCpx, QCpx );
            return TriangularPseudospectrum
                   ( UCpx, QCpx, shifts, invNorms, psCtrl );
            */
        }
        return QuasiTriangularPseudospectrum( U, Q, shifts, invNorms, psCtrl );
    }
}

template<typename Real>
Matrix<Int> Pseudospectrum
( const Matrix<Complex<Real>>& A, const Matrix<Complex<Real>>& shifts, 
  Matrix<Real>& invNorms, PseudospecCtrl<Real> psCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("Pseudospectrum"))
    typedef Complex<Real> C;

    Matrix<C> U( A );
    if( psCtrl.norm == PS_TWO_NORM )
    {
        if( psCtrl.schur )
        {
            Matrix<C> w;
            const bool fullTriangle = true;
            Schur( U, w, fullTriangle, psCtrl.schurCtrl );
            return TriangularPseudospectrum( U, shifts, invNorms, psCtrl );
        }
        else
        {
            Hessenberg( UPPER, U );
            return HessenbergPseudospectrum( U, shifts, invNorms, psCtrl );
        }
    }
    else
    {
        Matrix<C> Q;
        if( psCtrl.schur )
        {
            Matrix<C> w;
            const bool fullTriangle = true;
            Schur( U, w, Q, fullTriangle, psCtrl.schurCtrl );
            return TriangularPseudospectrum( U, Q, shifts, invNorms, psCtrl );
        }
        else
        {
            Matrix<C> t;
            Hessenberg( UPPER, U, t );
            Identity( Q, A.Height(), A.Height() );
            hessenberg::ApplyQ( UPPER, LEFT, NORMAL, U, t, Q );
            return HessenbergPseudospectrum( U, Q, shifts, invNorms, psCtrl );
        }
    }
}

template<typename Real>
DistMatrix<Int,VR,STAR> Pseudospectrum
( const DistMatrix<Real>& A, const DistMatrix<Complex<Real>,VR,STAR>& shifts,
  DistMatrix<Real,VR,STAR>& invNorms, PseudospecCtrl<Real> psCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("Pseudospectrum"))
    typedef Complex<Real> C;
    const Grid& g = A.Grid();

    if( psCtrl.forceComplexSchur )
    {
        DistMatrix<C> ACpx(g);
        Copy( A, ACpx );
        return Pseudospectrum( ACpx, shifts, invNorms, psCtrl );
    }

    if( !psCtrl.schur )
        LogicError("Real Hessenberg algorithm not yet supported");
    DistMatrix<Real> U( A );
    DistMatrix<C,VR,STAR> w(g);
    const bool fullTriangle = true;
    if( psCtrl.norm == PS_TWO_NORM )
    {
        Schur( U, w, fullTriangle, psCtrl.schurCtrl );
        if( psCtrl.forceComplexPs )
        {
            DistMatrix<C> UCpx(g);
            schur::RealToComplex( U, UCpx );
            return TriangularPseudospectrum( UCpx, shifts, invNorms, psCtrl );
        }
        return QuasiTriangularPseudospectrum( U, shifts, invNorms, psCtrl );
    }
    else
    {
        DistMatrix<Real> Q(g);
        Schur( U, w, Q, fullTriangle, psCtrl.schurCtrl );
        if( psCtrl.forceComplexPs )
        {
            LogicError("Real to complex full Schur not yet supported");
            /*
            DistMatrix<C> UCpx(g), QCpx(g);
            schur::RealToComplex( U, Q, UCpx, QCpx );
            return TriangularPseudospectrum
                   ( UCpx, QCpx, shifts, invNorms, psCtrl );
            */
        }
        return QuasiTriangularPseudospectrum( U, Q, shifts, invNorms, psCtrl );
    }
}

template<typename Real>
DistMatrix<Int,VR,STAR> Pseudospectrum
( const DistMatrix<Complex<Real>>& A, 
  const DistMatrix<Complex<Real>,VR,STAR>& shifts,
  DistMatrix<Real,VR,STAR>& invNorms, PseudospecCtrl<Real> psCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("Pseudospectrum"))
    typedef Complex<Real> C;

    const Grid& g = A.Grid();
    DistMatrix<C> U( A );

    if( psCtrl.norm == PS_TWO_NORM )
    {
        if( psCtrl.schur )
        {
            DistMatrix<C,VR,STAR> w(g);
            const bool fullTriangle = true;
            Schur( U, w, fullTriangle, psCtrl.schurCtrl );
            return TriangularPseudospectrum( U, shifts, invNorms, psCtrl );
        }
        else
        {
            Hessenberg( UPPER, U );
            return HessenbergPseudospectrum( U, shifts, invNorms, psCtrl );
        }
    }
    else
    {
        DistMatrix<C> Q(g);
        if( psCtrl.schur )
        {
            DistMatrix<C,VR,STAR> w(g);
            const bool fullTriangle = true;
            Schur( U, w, Q, fullTriangle, psCtrl.schurCtrl );
            return TriangularPseudospectrum( U, Q, shifts, invNorms, psCtrl );
        }
        else
        {
            DistMatrix<C,STAR,STAR> t(g);
            Hessenberg( UPPER, U, t );
            Identity( Q, U.Height(), U.Height() );
            hessenberg::ApplyQ( UPPER, LEFT, NORMAL, U, t, Q );
            return HessenbergPseudospectrum( U, Q, shifts, invNorms, psCtrl );
        }
    }
}

// Treat each pixel as being located a cell center and tesselate a box with
// said square cells
template<typename F>
Matrix<Int> TriangularPseudospectrum
( const Matrix<F>& U, Matrix<Base<F>>& invNormMap, 
  Complex<Base<F>> center, Base<F> realWidth, Base<F> imagWidth,
  Int realSize, Int imagSize, PseudospecCtrl<Base<F>> psCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("TriangularPseudospectrum"))
    typedef Base<F> Real;
    typedef Complex<Real> C;

    psCtrl.snapCtrl.realSize = realSize;
    psCtrl.snapCtrl.imagSize = imagSize;

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
        TriangularPseudospectrum( U, shifts, invNorms, psCtrl );

    // Rearrange the vectors into grids
    Matrix<Int> itCountMap; 
    pspec::ReshapeIntoGrid( realSize, imagSize, invNorms, invNormMap );
    pspec::ReshapeIntoGrid( realSize, imagSize, itCounts, itCountMap );
    return itCountMap;
}

template<typename F>
Matrix<Int> TriangularPseudospectrum
( const Matrix<F>& U, const Matrix<F>& Q, Matrix<Base<F>>& invNormMap, 
  Complex<Base<F>> center, Base<F> realWidth, Base<F> imagWidth,
  Int realSize, Int imagSize, PseudospecCtrl<Base<F>> psCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("TriangularPseudospectrum"))
    typedef Base<F> Real;
    typedef Complex<Real> C;

    psCtrl.snapCtrl.realSize = realSize;
    psCtrl.snapCtrl.imagSize = imagSize;

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
        TriangularPseudospectrum( U, Q, shifts, invNorms, psCtrl );

    // Rearrange the vectors into grids
    Matrix<Int> itCountMap; 
    pspec::ReshapeIntoGrid( realSize, imagSize, invNorms, invNormMap );
    pspec::ReshapeIntoGrid( realSize, imagSize, itCounts, itCountMap );
    return itCountMap;
}

template<typename Real>
Matrix<Int> QuasiTriangularPseudospectrum
( const Matrix<Real>& U, 
  Matrix<Real>& invNormMap, 
  Complex<Real> center, Real realWidth, Real imagWidth,
  Int realSize, Int imagSize, PseudospecCtrl<Real> psCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("QuasiTriangularPseudospectrum"))
    typedef Complex<Real> C;

    psCtrl.snapCtrl.realSize = realSize;
    psCtrl.snapCtrl.imagSize = imagSize;

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
        QuasiTriangularPseudospectrum( U, shifts, invNorms, psCtrl );

    // Rearrange the vectors into grids
    Matrix<Int> itCountMap; 
    pspec::ReshapeIntoGrid( realSize, imagSize, invNorms, invNormMap );
    pspec::ReshapeIntoGrid( realSize, imagSize, itCounts, itCountMap );
    return itCountMap;
}

template<typename Real>
Matrix<Int> QuasiTriangularPseudospectrum
( const Matrix<Real>& U, 
  const Matrix<Real>& Q,
  Matrix<Real>& invNormMap, 
  Complex<Real> center, Real realWidth, Real imagWidth,
  Int realSize, Int imagSize, PseudospecCtrl<Real> psCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("QuasiTriangularPseudospectrum"))
    typedef Complex<Real> C;

    psCtrl.snapCtrl.realSize = realSize;
    psCtrl.snapCtrl.imagSize = imagSize;

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
        QuasiTriangularPseudospectrum( U, Q, shifts, invNorms, psCtrl );

    // Rearrange the vectors into grids
    Matrix<Int> itCountMap; 
    pspec::ReshapeIntoGrid( realSize, imagSize, invNorms, invNormMap );
    pspec::ReshapeIntoGrid( realSize, imagSize, itCounts, itCountMap );
    return itCountMap;
}

template<typename F>
Matrix<Int> HessenbergPseudospectrum
( const Matrix<F>& H, Matrix<Base<F>>& invNormMap, 
  Complex<Base<F>> center, Base<F> realWidth, Base<F> imagWidth,
  Int realSize, Int imagSize, PseudospecCtrl<Base<F>> psCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("HessenbergPseudospectrum"))
    typedef Base<F> Real;
    typedef Complex<Real> C;

    psCtrl.snapCtrl.realSize = realSize;
    psCtrl.snapCtrl.imagSize = imagSize;

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
    auto itCounts = HessenbergPseudospectrum( H, shifts, invNorms, psCtrl );

    // Rearrange the vectors into grids
    Matrix<Int> itCountMap; 
    pspec::ReshapeIntoGrid( realSize, imagSize, invNorms, invNormMap );
    pspec::ReshapeIntoGrid( realSize, imagSize, itCounts, itCountMap );
    return itCountMap;
}

template<typename F>
Matrix<Int> HessenbergPseudospectrum
( const Matrix<F>& H, const Matrix<F>& Q, Matrix<Base<F>>& invNormMap, 
  Complex<Base<F>> center, Base<F> realWidth, Base<F> imagWidth,
  Int realSize, Int imagSize, PseudospecCtrl<Base<F>> psCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("HessenbergPseudospectrum"))
    typedef Base<F> Real;
    typedef Complex<Real> C;

    psCtrl.snapCtrl.realSize = realSize;
    psCtrl.snapCtrl.imagSize = imagSize;

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
    auto itCounts = HessenbergPseudospectrum( H, Q, shifts, invNorms, psCtrl );

    // Rearrange the vectors into grids
    Matrix<Int> itCountMap; 
    pspec::ReshapeIntoGrid( realSize, imagSize, invNorms, invNormMap );
    pspec::ReshapeIntoGrid( realSize, imagSize, itCounts, itCountMap );
    return itCountMap;
}

template<typename F>
DistMatrix<Int> TriangularPseudospectrum
( const DistMatrix<F>& U, DistMatrix<Base<F>>& invNormMap, 
  Complex<Base<F>> center, Base<F> realWidth, Base<F> imagWidth, 
  Int realSize, Int imagSize, PseudospecCtrl<Base<F>> psCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("TriangularPseudospectrum"))
    typedef Base<F> Real;
    typedef Complex<Real> C;
    const Grid& g = U.Grid();

    psCtrl.snapCtrl.realSize = realSize;
    psCtrl.snapCtrl.imagSize = imagSize;

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
    auto itCounts = TriangularPseudospectrum( U, shifts, invNorms, psCtrl );

    // Rearrange the vectors into grids
    DistMatrix<Int> itCountMap(g); 
    pspec::ReshapeIntoGrid( realSize, imagSize, invNorms, invNormMap );
    pspec::ReshapeIntoGrid( realSize, imagSize, itCounts, itCountMap );
    return itCountMap;
}

template<typename F>
DistMatrix<Int> TriangularPseudospectrum
( const DistMatrix<F>& U, const DistMatrix<F>& Q, 
  DistMatrix<Base<F>>& invNormMap, 
  Complex<Base<F>> center, Base<F> realWidth, Base<F> imagWidth, 
  Int realSize, Int imagSize, PseudospecCtrl<Base<F>> psCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("TriangularPseudospectrum"))
    typedef Base<F> Real;
    typedef Complex<Real> C;
    const Grid& g = U.Grid();

    psCtrl.snapCtrl.realSize = realSize;
    psCtrl.snapCtrl.imagSize = imagSize;

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
    auto itCounts = TriangularPseudospectrum( U, Q, shifts, invNorms, psCtrl );

    // Rearrange the vectors into grids
    DistMatrix<Int> itCountMap(g); 
    pspec::ReshapeIntoGrid( realSize, imagSize, invNorms, invNormMap );
    pspec::ReshapeIntoGrid( realSize, imagSize, itCounts, itCountMap );
    return itCountMap;
}

template<typename Real>
DistMatrix<Int> QuasiTriangularPseudospectrum
( const DistMatrix<Real>& U, 
  DistMatrix<Real>& invNormMap, 
  Complex<Real> center, Real realWidth, Real imagWidth, 
  Int realSize, Int imagSize, PseudospecCtrl<Real> psCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("QuasiTriangularPseudospectrum"))
    typedef Complex<Real> C;
    const Grid& g = U.Grid();

    psCtrl.snapCtrl.realSize = realSize;
    psCtrl.snapCtrl.imagSize = imagSize;

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
        QuasiTriangularPseudospectrum( U, shifts, invNorms, psCtrl );

    // Rearrange the vectors into grids
    DistMatrix<Int> itCountMap(g); 
    pspec::ReshapeIntoGrid( realSize, imagSize, invNorms, invNormMap );
    pspec::ReshapeIntoGrid( realSize, imagSize, itCounts, itCountMap );
    return itCountMap;
}

template<typename Real>
DistMatrix<Int> QuasiTriangularPseudospectrum
( const DistMatrix<Real>& U, 
  const DistMatrix<Real>& Q,
  DistMatrix<Real>& invNormMap, 
  Complex<Real> center, Real realWidth, Real imagWidth, 
  Int realSize, Int imagSize, PseudospecCtrl<Real> psCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("QuasiTriangularPseudospectrum"))
    typedef Complex<Real> C;
    const Grid& g = U.Grid();

    psCtrl.snapCtrl.realSize = realSize;
    psCtrl.snapCtrl.imagSize = imagSize;

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
        QuasiTriangularPseudospectrum( U, Q, shifts, invNorms, psCtrl );

    // Rearrange the vectors into grids
    DistMatrix<Int> itCountMap(g); 
    pspec::ReshapeIntoGrid( realSize, imagSize, invNorms, invNormMap );
    pspec::ReshapeIntoGrid( realSize, imagSize, itCounts, itCountMap );
    return itCountMap;
}

template<typename F>
DistMatrix<Int> HessenbergPseudospectrum
( const DistMatrix<F>& H, DistMatrix<Base<F>>& invNormMap, 
  Complex<Base<F>> center, Base<F> realWidth, Base<F> imagWidth, 
  Int realSize, Int imagSize, PseudospecCtrl<Base<F>> psCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("HessenbergPseudospectrum"))
    typedef Base<F> Real;
    typedef Complex<Real> C;
    const Grid& g = H.Grid();

    psCtrl.snapCtrl.realSize = realSize;
    psCtrl.snapCtrl.imagSize = imagSize;

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
    auto itCounts = HessenbergPseudospectrum( H, shifts, invNorms, psCtrl );

    // Rearrange the vectors into grids
    DistMatrix<Int> itCountMap(g); 
    pspec::ReshapeIntoGrid( realSize, imagSize, invNorms, invNormMap );
    pspec::ReshapeIntoGrid( realSize, imagSize, itCounts, itCountMap );
    return itCountMap;
}

template<typename F>
DistMatrix<Int> HessenbergPseudospectrum
( const DistMatrix<F>& H, 
  const DistMatrix<F>& Q,
  DistMatrix<Base<F>>& invNormMap, 
  Complex<Base<F>> center, Base<F> realWidth, Base<F> imagWidth, 
  Int realSize, Int imagSize, PseudospecCtrl<Base<F>> psCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("HessenbergPseudospectrum"))
    typedef Base<F> Real;
    typedef Complex<Real> C;
    const Grid& g = H.Grid();

    psCtrl.snapCtrl.realSize = realSize;
    psCtrl.snapCtrl.imagSize = imagSize;

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
    auto itCounts = HessenbergPseudospectrum( H, Q, shifts, invNorms, psCtrl );

    // Rearrange the vectors into grids
    DistMatrix<Int> itCountMap(g); 
    pspec::ReshapeIntoGrid( realSize, imagSize, invNorms, invNormMap );
    pspec::ReshapeIntoGrid( realSize, imagSize, itCounts, itCountMap );
    return itCountMap;
}

template<typename F>
Matrix<Int> Pseudospectrum
( const Matrix<F>& A, Matrix<Base<F>>& invNormMap, 
  Complex<Base<F>> center, Base<F> realWidth, Base<F> imagWidth,
  Int realSize, Int imagSize, PseudospecCtrl<Base<F>> psCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("Pseudospectrum"))
    typedef Base<F> Real;
    typedef Complex<Real> C;

    psCtrl.snapCtrl.realSize = realSize;
    psCtrl.snapCtrl.imagSize = imagSize;

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
    auto itCounts = Pseudospectrum( A, shifts, invNorms, psCtrl );

    // Rearrange the vectors into grids
    Matrix<Int> itCountMap; 
    pspec::ReshapeIntoGrid( realSize, imagSize, invNorms, invNormMap );
    pspec::ReshapeIntoGrid( realSize, imagSize, itCounts, itCountMap );
    return itCountMap;
}

template<typename F>
DistMatrix<Int> Pseudospectrum
( const DistMatrix<F>& A, DistMatrix<Base<F>>& invNormMap, 
  Complex<Base<F>> center, Base<F> realWidth, Base<F> imagWidth, 
  Int realSize, Int imagSize, PseudospecCtrl<Base<F>> psCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("Pseudospectrum"))
    typedef Base<F> Real;
    typedef Complex<Real> C;
    const Grid& g = A.Grid();

    psCtrl.snapCtrl.realSize = realSize;
    psCtrl.snapCtrl.imagSize = imagSize;

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
    auto itCounts = Pseudospectrum( A, shifts, invNorms, psCtrl );

    // Rearrange the vectors into grids
    DistMatrix<Int> itCountMap(g); 
    pspec::ReshapeIntoGrid( realSize, imagSize, invNorms, invNormMap );
    pspec::ReshapeIntoGrid( realSize, imagSize, itCounts, itCountMap );
    return itCountMap;
}

template<typename F>
Matrix<Int> TriangularPseudospectrum
( const Matrix<F>& U, Matrix<Base<F>>& invNormMap, 
  Complex<Base<F>> center, Int realSize, Int imagSize, 
  PseudospecCtrl<Base<F>> psCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("TriangularPseudospectrum"))

    auto diag = U.GetDiagonal();
    const Base<F> radius = MaxNorm( diag );
    const Base<F> oneNorm = OneNorm( U );

    // Essentially three cases are handled here:
    // 1) The zero matrix (force the pseudospectrum width to 1)
    // 2) Typical matrices (use a small multiple of the spectral radius)
    // 3) Highly non-normal matrices (e.g., triangular with zero main diagonal)
    Base<F> width;
    if( oneNorm == Base<F>(0) )
    {
        width = 1;
        if( psCtrl.progress )
            std::cout << "Setting width to 1 to handle zero matrix" 
                      << std::endl;
    }
    else if( radius >= 0.2*oneNorm )
    {
        width = 2.5*radius;
        if( psCtrl.progress )
            std::cout << "Setting width to " << width 
                      << " based on the spectral radius, " << radius 
                      << std::endl;
    }
    else
    {
        width = 0.8*oneNorm;
        if( psCtrl.progress )
            std::cout << "Setting width to " << width 
                      << " based on the one norm, " << oneNorm << std::endl;
    }

    return TriangularPseudospectrum
           ( U, invNormMap, center, width, width, realSize, imagSize, psCtrl );
}

template<typename F>
Matrix<Int> TriangularPseudospectrum
( const Matrix<F>& U, const Matrix<F>& Q, Matrix<Base<F>>& invNormMap, 
  Complex<Base<F>> center, Int realSize, Int imagSize, 
  PseudospecCtrl<Base<F>> psCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("TriangularPseudospectrum"))

    auto diag = U.GetDiagonal();
    const Base<F> radius = MaxNorm( diag );
    const Base<F> oneNorm = OneNorm( U );

    // Essentially three cases are handled here:
    // 1) The zero matrix (force the pseudospectrum width to 1)
    // 2) Typical matrices (use a small multiple of the spectral radius)
    // 3) Highly non-normal matrices (e.g., triangular with zero main diagonal)
    Base<F> width;
    if( oneNorm == Base<F>(0) )
    {
        width = 1;
        if( psCtrl.progress )
            std::cout << "Setting width to 1 to handle zero matrix" 
                      << std::endl;
    }
    else if( radius >= 0.2*oneNorm )
    {
        width = 2.5*radius;
        if( psCtrl.progress )
            std::cout << "Setting width to " << width 
                      << " based on the spectral radius, " << radius 
                      << std::endl;
    }
    else
    {
        width = 0.8*oneNorm;
        if( psCtrl.progress )
            std::cout << "Setting width to " << width 
                      << " based on the one norm, " << oneNorm << std::endl;
    }

    return TriangularPseudospectrum
           ( U, Q, invNormMap, center, width, width, realSize, imagSize, 
             psCtrl );
}

template<typename Real>
Matrix<Int> QuasiTriangularPseudospectrum
( const Matrix<Real>& U, 
  Matrix<Real>& invNormMap, 
  Complex<Real> center, Int realSize, Int imagSize, 
  PseudospecCtrl<Real> psCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("QuasiTriangularPseudospectrum"))

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
            std::cout << "Setting width to 1 to handle zero matrix" 
                      << std::endl;
    }
    else if( radius >= 0.2*oneNorm )
    {
        width = 2.5*radius;
        if( psCtrl.progress )
            std::cout << "Setting width to " << width 
                      << " based on the spectral radius, " << radius 
                      << std::endl;
    }
    else
    {
        width = 0.8*oneNorm;
        if( psCtrl.progress )
            std::cout << "Setting width to " << width 
                      << " based on the one norm, " << oneNorm << std::endl;
    }

    return QuasiTriangularPseudospectrum
           ( U, invNormMap, center, width, width, realSize, imagSize, psCtrl );
}

template<typename Real>
Matrix<Int> QuasiTriangularPseudospectrum
( const Matrix<Real>& U, 
  const Matrix<Real>& Q,
  Matrix<Real>& invNormMap, 
  Complex<Real> center, Int realSize, Int imagSize, 
  PseudospecCtrl<Real> psCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("QuasiTriangularPseudospectrum"))

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
            std::cout << "Setting width to 1 to handle zero matrix" 
                      << std::endl;
    }
    else if( radius >= 0.2*oneNorm )
    {
        width = 2.5*radius;
        if( psCtrl.progress )
            std::cout << "Setting width to " << width 
                      << " based on the spectral radius, " << radius 
                      << std::endl;
    }
    else
    {
        width = 0.8*oneNorm;
        if( psCtrl.progress )
            std::cout << "Setting width to " << width 
                      << " based on the one norm, " << oneNorm << std::endl;
    }

    return QuasiTriangularPseudospectrum
           ( U, Q, invNormMap, center, width, width, realSize, imagSize, 
             psCtrl );
}

template<typename F>
Matrix<Int> HessenbergPseudospectrum
( const Matrix<F>& H, Matrix<Base<F>>& invNormMap, 
  Complex<Base<F>> center, Int realSize, Int imagSize, 
  PseudospecCtrl<Base<F>> psCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("HessenbergPseudospectrum"))

    const Base<F> infNorm = InfinityNorm( H );
    const Base<F> oneNorm = OneNorm( H );
    Base<F> width;
    if( oneNorm == Base<F>(0) )
    {
        width = 1;
        if( psCtrl.progress )
            std::cout << "Setting width to 1 to handle zero matrix" 
                      << std::endl;
    }
    else
    {
        width = 0.8*Max(oneNorm,infNorm);
        if( psCtrl.progress )
            std::cout << "Setting width to " << width 
                      << " based on the one norm, " << oneNorm 
                      << ", and infinity norm, " << infNorm << std::endl;
    }

    return HessenbergPseudospectrum
           ( H, invNormMap, center, width, width, realSize, imagSize, psCtrl );
}

template<typename F>
Matrix<Int> HessenbergPseudospectrum
( const Matrix<F>& H, const Matrix<F>& Q, Matrix<Base<F>>& invNormMap, 
  Complex<Base<F>> center, Int realSize, Int imagSize, 
  PseudospecCtrl<Base<F>> psCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("HessenbergPseudospectrum"))

    const Base<F> infNorm = InfinityNorm( H );
    const Base<F> oneNorm = OneNorm( H );
    Base<F> width;
    if( oneNorm == Base<F>(0) )
    {
        width = 1;
        if( psCtrl.progress )
            std::cout << "Setting width to 1 to handle zero matrix" 
                      << std::endl;
    }
    else
    {
        width = 0.8*Max(oneNorm,infNorm);
        if( psCtrl.progress )
            std::cout << "Setting width to " << width 
                      << " based on the one norm, " << oneNorm 
                      << ", and infinity norm, " << infNorm << std::endl;
    }

    return HessenbergPseudospectrum
           ( H, Q, invNormMap, center, width, width, realSize, imagSize, 
             psCtrl );
}

template<typename F>
DistMatrix<Int> TriangularPseudospectrum
( const DistMatrix<F>& U, DistMatrix<Base<F>>& invNormMap, 
  Complex<Base<F>> center, Int realSize, Int imagSize,
  PseudospecCtrl<Base<F>> psCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("TriangularPseudospectrum"))

    auto diag = U.GetDiagonal();
    const Base<F> radius = MaxNorm( diag );
    const Base<F> oneNorm = OneNorm( U );

    // Essentially three cases are handled here:
    // 1) The zero matrix (force the pseudospectrum width to 1)
    // 2) Typical matrices (use a small multiple of the spectral radius)
    // 3) Highly non-normal matrices (e.g., triangular with zero main diagonal)
    Base<F> width;
    if( oneNorm == Base<F>(0) )
    {
        width = 1;
        if( psCtrl.progress && U.Grid().Rank() == 0 )
            std::cout << "Setting width to 1 to handle zero matrix"
                      << std::endl;
    }
    else if( radius >= 0.2*oneNorm )
    {
        width = 2.5*radius;
        if( psCtrl.progress && U.Grid().Rank() == 0 )
            std::cout << "Setting width to " << width 
                      << " based on the spectral radius, " << radius 
                      << std::endl;
    }
    else
    {
        width = 0.8*oneNorm;
        if( psCtrl.progress && U.Grid().Rank() == 0 )
            std::cout << "Setting width to " << width
                      << " based on the one norm, " << oneNorm << std::endl;
    }

    return TriangularPseudospectrum
           ( U, invNormMap, center, width, width, realSize, imagSize, psCtrl );
}

template<typename F>
DistMatrix<Int> TriangularPseudospectrum
( const DistMatrix<F>& U, 
  const DistMatrix<F>& Q, DistMatrix<Base<F>>& invNormMap, 
  Complex<Base<F>> center, Int realSize, Int imagSize,
  PseudospecCtrl<Base<F>> psCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("TriangularPseudospectrum"))

    auto diag = U.GetDiagonal();
    const Base<F> radius = MaxNorm( diag );
    const Base<F> oneNorm = OneNorm( U );

    // Essentially three cases are handled here:
    // 1) The zero matrix (force the pseudospectrum width to 1)
    // 2) Typical matrices (use a small multiple of the spectral radius)
    // 3) Highly non-normal matrices (e.g., triangular with zero main diagonal)
    Base<F> width;
    if( oneNorm == Base<F>(0) )
    {
        width = 1;
        if( psCtrl.progress && U.Grid().Rank() == 0 )
            std::cout << "Setting width to 1 to handle zero matrix"
                      << std::endl;
    }
    else if( radius >= 0.2*oneNorm )
    {
        width = 2.5*radius;
        if( psCtrl.progress && U.Grid().Rank() == 0 )
            std::cout << "Setting width to " << width 
                      << " based on the spectral radius, " << radius 
                      << std::endl;
    }
    else
    {
        width = 0.8*oneNorm;
        if( psCtrl.progress && U.Grid().Rank() == 0 )
            std::cout << "Setting width to " << width
                      << " based on the one norm, " << oneNorm << std::endl;
    }

    return TriangularPseudospectrum
           ( U, Q, invNormMap, center, width, width, realSize, imagSize, 
             psCtrl );
}

template<typename Real>
DistMatrix<Int> QuasiTriangularPseudospectrum
( const DistMatrix<Real>& U, 
  DistMatrix<Real>& invNormMap, 
  Complex<Real> center, Int realSize, Int imagSize,
  PseudospecCtrl<Real> psCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("QuasiTriangularPseudospectrum"))

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
        if( psCtrl.progress && U.Grid().Rank() == 0 )
            std::cout << "Setting width to 1 to handle zero matrix"
                      << std::endl;
    }
    else if( radius >= 0.2*oneNorm )
    {
        width = 2.5*radius;
        if( psCtrl.progress && U.Grid().Rank() == 0 )
            std::cout << "Setting width to " << width 
                      << " based on the spectral radius, " << radius 
                      << std::endl;
    }
    else
    {
        width = 0.8*oneNorm;
        if( psCtrl.progress && U.Grid().Rank() == 0 )
            std::cout << "Setting width to " << width
                      << " based on the one norm, " << oneNorm << std::endl;
    }

    return QuasiTriangularPseudospectrum
           ( U, invNormMap, center, width, width, realSize, imagSize, psCtrl );
}

template<typename Real>
DistMatrix<Int> QuasiTriangularPseudospectrum
( const DistMatrix<Real>& U, 
  const DistMatrix<Real>& Q,
  DistMatrix<Real>& invNormMap, 
  Complex<Real> center, Int realSize, Int imagSize,
  PseudospecCtrl<Real> psCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("QuasiTriangularPseudospectrum"))

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
        if( psCtrl.progress && U.Grid().Rank() == 0 )
            std::cout << "Setting width to 1 to handle zero matrix"
                      << std::endl;
    }
    else if( radius >= 0.2*oneNorm )
    {
        width = 2.5*radius;
        if( psCtrl.progress && U.Grid().Rank() == 0 )
            std::cout << "Setting width to " << width 
                      << " based on the spectral radius, " << radius 
                      << std::endl;
    }
    else
    {
        width = 0.8*oneNorm;
        if( psCtrl.progress && U.Grid().Rank() == 0 )
            std::cout << "Setting width to " << width
                      << " based on the one norm, " << oneNorm << std::endl;
    }

    return QuasiTriangularPseudospectrum
           ( U, Q, invNormMap, center, width, width, realSize, imagSize, 
             psCtrl );
}

template<typename F>
DistMatrix<Int> HessenbergPseudospectrum
( const DistMatrix<F>& H, DistMatrix<Base<F>>& invNormMap, 
  Complex<Base<F>> center, Int realSize, Int imagSize,
  PseudospecCtrl<Base<F>> psCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("HessenbergPseudospectrum"))

    const Base<F> oneNorm = OneNorm( H );
    const Base<F> infNorm = InfinityNorm( H );
    Base<F> width;
    if( oneNorm == Base<F>(0) )
    {
        width = 1;
        if( psCtrl.progress && H.Grid().Rank() == 0 )
            std::cout << "Setting width to 1 to handle zero matrix"
                      << std::endl;
    }
    else
    {
        width = 0.8*Max(oneNorm,infNorm);
        if( psCtrl.progress && H.Grid().Rank() == 0 )
            std::cout << "Setting width to " << width 
                      << " based on the one norm, " << oneNorm 
                      << ", and infinity norm, " << infNorm << std::endl;
    }

    return HessenbergPseudospectrum
           ( H, invNormMap, center, width, width, realSize, imagSize, psCtrl );
}

template<typename F>
DistMatrix<Int> HessenbergPseudospectrum
( const DistMatrix<F>& H, 
  const DistMatrix<F>& Q, DistMatrix<Base<F>>& invNormMap, 
  Complex<Base<F>> center, Int realSize, Int imagSize,
  PseudospecCtrl<Base<F>> psCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("HessenbergPseudospectrum"))

    const Base<F> oneNorm = OneNorm( H );
    const Base<F> infNorm = InfinityNorm( H );
    Base<F> width;
    if( oneNorm == Base<F>(0) )
    {
        width = 1;
        if( psCtrl.progress && H.Grid().Rank() == 0 )
            std::cout << "Setting width to 1 to handle zero matrix"
                      << std::endl;
    }
    else
    {
        width = 0.8*Max(oneNorm,infNorm);
        if( psCtrl.progress && H.Grid().Rank() == 0 )
            std::cout << "Setting width to " << width 
                      << " based on the one norm, " << oneNorm 
                      << ", and infinity norm, " << infNorm << std::endl;
    }

    return HessenbergPseudospectrum
           ( H, Q, invNormMap, center, width, width, realSize, imagSize, 
             psCtrl );
}

template<typename Real>
Matrix<Int> Pseudospectrum
( const Matrix<Complex<Real>>& A, Matrix<Real>& invNormMap, 
  Complex<Real> center, Int realSize, Int imagSize, 
  PseudospecCtrl<Real> psCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("Pseudospectrum"))
    typedef Complex<Real> C;

    Matrix<C> B( A );
    if( psCtrl.norm == PS_TWO_NORM )
    {
        if( psCtrl.schur )
        {
            Matrix<C> w;
            const bool fullTriangle = true;
            Schur( B, w, fullTriangle, psCtrl.schurCtrl );
            return TriangularPseudospectrum
                   ( B, invNormMap, center, realSize, imagSize, psCtrl );
        }
        else
        {
            Hessenberg( UPPER, B );
            return HessenbergPseudospectrum
                   ( B, invNormMap, center, realSize, imagSize, psCtrl );
        }
    }
    else
    {
        Matrix<C> Q;
        if( psCtrl.schur )
        {
            Matrix<C> w;
            const bool fullTriangle = true;
            Schur( B, w, Q, fullTriangle, psCtrl.schurCtrl );
            return TriangularPseudospectrum
                   ( B, Q, invNormMap, center, realSize, imagSize, psCtrl );
        }
        else
        {
            Matrix<C> t;
            Hessenberg( UPPER, B, t );
            Identity( Q, B.Height(), B.Height() );
            hessenberg::ApplyQ( UPPER, LEFT, NORMAL, B, t, Q );
            return HessenbergPseudospectrum
                   ( B, Q, invNormMap, center, realSize, imagSize, psCtrl );
        }
    }
}

template<typename Real>
Matrix<Int> Pseudospectrum
( const Matrix<Real>& A, Matrix<Real>& invNormMap, 
  Complex<Real> center, Int realSize, Int imagSize, 
  PseudospecCtrl<Real> psCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("Pseudospectrum"))
    typedef Complex<Real> C;

    if( psCtrl.forceComplexSchur )
    {
        Matrix<C> ACpx;
        Copy( A, ACpx );
        return Pseudospectrum
               ( ACpx, invNormMap, center, realSize, imagSize, psCtrl );
    }

    if( !psCtrl.schur )
        LogicError("Real Hessenberg algorithm not yet supported");
    Matrix<Real> B( A );
    Matrix<C> w;
    const bool fullTriangle = true;
    if( psCtrl.norm == PS_TWO_NORM )
    {
        Schur( B, w, fullTriangle, psCtrl.schurCtrl );
        if( psCtrl.forceComplexPs )
        {
            Matrix<C> BCpx;
            schur::RealToComplex( B, BCpx );
            return TriangularPseudospectrum
                   ( BCpx, invNormMap, center, realSize, imagSize, psCtrl );
        }
        return QuasiTriangularPseudospectrum
               ( B, invNormMap, center, realSize, imagSize, psCtrl );
    }
    else
    {
        Matrix<Real> Q;
        Schur( B, w, Q, fullTriangle, psCtrl.schurCtrl );
        if( psCtrl.forceComplexPs )
        {
            LogicError("Real to complex full Schur not yet supported");
            /*
            Matrix<C> BCpx, QCpx;
            schur::RealToComplex( B, Q, BCpx, QCpx );
            return TriangularPseudospectrum
                   ( BCpx, invNormMap, center, realSize, imagSize, psCtrl );
            */
        }
        return QuasiTriangularPseudospectrum
               ( B, Q, invNormMap, center, realSize, imagSize, psCtrl );
    }
}

template<typename Real>
DistMatrix<Int> Pseudospectrum
( const DistMatrix<Complex<Real>>& A, DistMatrix<Real>& invNormMap, 
  Complex<Real> center, Int realSize, Int imagSize,
  PseudospecCtrl<Real> psCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("Pseudospectrum"))
    typedef Complex<Real> C;

    const Grid& g = A.Grid();
    DistMatrix<C> B( A );

    if( psCtrl.norm == PS_TWO_NORM )
    {
        if( psCtrl.schur )
        {
            DistMatrix<C,VR,STAR> w(g);
            const bool fullTriangle = true;
            Schur( B, w, fullTriangle, psCtrl.schurCtrl );
            return TriangularPseudospectrum
                   ( B, invNormMap, center, realSize, imagSize, psCtrl );
        }
        else
        {
            Hessenberg( UPPER, B );
            return HessenbergPseudospectrum
                   ( B, invNormMap, center, realSize, imagSize, psCtrl );
        }
    }
    else
    {
        DistMatrix<C> Q(g);
        if( psCtrl.schur )
        {
            DistMatrix<C,VR,STAR> w(g);
            const bool fullTriangle = true;
            Schur( B, w, Q, fullTriangle, psCtrl.schurCtrl );
            return TriangularPseudospectrum
                   ( B, Q, invNormMap, center, realSize, imagSize, psCtrl );
        }
        else
        {
            DistMatrix<C,STAR,STAR> t(g);
            Hessenberg( UPPER, B, t );
            Identity( Q, B.Height(), B.Height() );
            hessenberg::ApplyQ( UPPER, LEFT, NORMAL, B, t, Q );
            return HessenbergPseudospectrum
                   ( B, Q, invNormMap, center, realSize, imagSize, psCtrl );
        }
    }
}

template<typename Real>
DistMatrix<Int> Pseudospectrum
( const DistMatrix<Real>& A, DistMatrix<Real>& invNormMap, 
  Complex<Real> center, Int realSize, Int imagSize,
  PseudospecCtrl<Real> psCtrl )
{
    DEBUG_ONLY(CallStackEntry cse("Pseudospectrum"))
    typedef Complex<Real> C;
    const Grid& g = A.Grid();

    if( psCtrl.forceComplexSchur )
    {
        DistMatrix<C> ACpx(g);
        Copy( A, ACpx );
        return Pseudospectrum
               ( ACpx, invNormMap, center, realSize, imagSize, psCtrl );
    }

    if( !psCtrl.schur )
        LogicError("Real Hessenberg algorithm not yet supported");
    DistMatrix<Real> B( A );

    DistMatrix<C,VR,STAR> w(g);
    const bool fullTriangle = true;
    if( psCtrl.norm == PS_TWO_NORM )
    {
        Schur( B, w, fullTriangle, psCtrl.schurCtrl );
        if( psCtrl.forceComplexPs ) 
        {
            DistMatrix<C> BCpx(g);
            schur::RealToComplex( B, BCpx );
            return TriangularPseudospectrum
                   ( BCpx, invNormMap, center, realSize, imagSize, psCtrl );
        }
        return QuasiTriangularPseudospectrum
               ( B, invNormMap, center, realSize, imagSize, psCtrl );
    }
    else
    {
        DistMatrix<Real> Q(g); 
        Schur( B, w, Q, fullTriangle, psCtrl.schurCtrl );
        if( psCtrl.forceComplexPs ) 
        {
            LogicError("Real to complex full Schur not yet supported");
            /*
            DistMatrix<C> BCpx(g), QCpx(g);
            schur::RealToComplex( B, Q, BCpx, QCpx );
            return TriangularPseudospectrum
                   ( BCpx, QCpx, invNormMap, center, realSize, imagSize, 
                     psCtrl );
            */
        }
        return QuasiTriangularPseudospectrum
               ( B, Q, invNormMap, center, realSize, imagSize, psCtrl );
    }
}

#define PROTO_REAL(Real) \
  template Matrix<Int> TriangularPseudospectrum \
  ( const Matrix<Complex<Real>>& U, const Matrix<Complex<Real>>& shifts, \
    Matrix<Real>& invNorms, \
    PseudospecCtrl<Real> psCtrl ); \
  template DistMatrix<Int,VR,STAR> TriangularPseudospectrum \
  ( const DistMatrix<Complex<Real>>& U, \
    const DistMatrix<Complex<Real>,VR,STAR>& shifts, \
    DistMatrix<Real,VR,STAR>& invNorms, \
    PseudospecCtrl<Real> psCtrl ); \
  template Matrix<Int> TriangularPseudospectrum \
  ( const Matrix<Complex<Real>>& U, const Matrix<Complex<Real>>& Q, \
    const Matrix<Complex<Real>>& shifts, \
    Matrix<Real>& invNorms, \
    PseudospecCtrl<Real> psCtrl ); \
  template DistMatrix<Int,VR,STAR> TriangularPseudospectrum \
  ( const DistMatrix<Complex<Real>>& U, \
    const DistMatrix<Complex<Real>>& Q, \
    const DistMatrix<Complex<Real>,VR,STAR>& shifts, \
    DistMatrix<Real,VR,STAR>& invNorms, \
    PseudospecCtrl<Real> psCtrl ); \
  template Matrix<Int> TriangularPseudospectrum \
  ( const Matrix<Real>& U, const Matrix<Complex<Real>>& shifts, \
    Matrix<Real>& invNorms, \
    PseudospecCtrl<Real> psCtrl ); \
  template DistMatrix<Int,VR,STAR> TriangularPseudospectrum \
  ( const DistMatrix<Real>& U, \
    const DistMatrix<Complex<Real>,VR,STAR>& shifts, \
    DistMatrix<Real,VR,STAR>& invNorms, \
    PseudospecCtrl<Real> psCtrl ); \
  template Matrix<Int> TriangularPseudospectrum \
  ( const Matrix<Real>& U, const Matrix<Real>& Q, \
    const Matrix<Complex<Real>>& shifts, Matrix<Real>& invNorms, \
    PseudospecCtrl<Real> psCtrl ); \
  template DistMatrix<Int,VR,STAR> TriangularPseudospectrum \
  ( const DistMatrix<Real>& U, \
    const DistMatrix<Real>& Q, \
    const DistMatrix<Complex<Real>,VR,STAR>& shifts, \
    DistMatrix<Real,VR,STAR>& invNorms, \
    PseudospecCtrl<Real> psCtrl ); \
  template Matrix<Int> QuasiTriangularPseudospectrum \
  ( const Matrix<Real>& U, \
    const Matrix<Complex<Real>>& shifts, \
    Matrix<Real>& invNorms, \
    PseudospecCtrl<Real> psCtrl ); \
  template DistMatrix<Int,VR,STAR> QuasiTriangularPseudospectrum \
  ( const DistMatrix<Real>& U, \
    const DistMatrix<Complex<Real>,VR,STAR>& shifts, \
    DistMatrix<Real,VR,STAR>& invNorms, \
    PseudospecCtrl<Real> psCtrl ); \
  template Matrix<Int> QuasiTriangularPseudospectrum \
  ( const Matrix<Real>& U, const Matrix<Real>& Q, \
    const Matrix<Complex<Real>>& shifts, \
    Matrix<Real>& invNorms, \
    PseudospecCtrl<Real> psCtrl ); \
  template DistMatrix<Int,VR,STAR> QuasiTriangularPseudospectrum \
  ( const DistMatrix<Real>& U, \
    const DistMatrix<Real>& Q, \
    const DistMatrix<Complex<Real>,VR,STAR>& shifts, \
    DistMatrix<Real,VR,STAR>& invNorms, \
    PseudospecCtrl<Real> psCtrl ); \
  template Matrix<Int> HessenbergPseudospectrum \
  ( const Matrix<Complex<Real>>& H, const Matrix<Complex<Real>>& shifts, \
    Matrix<Real>& invNorms, \
    PseudospecCtrl<Real> psCtrl ); \
  template DistMatrix<Int,VR,STAR> HessenbergPseudospectrum \
  ( const DistMatrix<Complex<Real>>& H, \
    const DistMatrix<Complex<Real>,VR,STAR>& shifts, \
    DistMatrix<Real,VR,STAR>& invNorms, \
    PseudospecCtrl<Real> psCtrl ); \
  template Matrix<Int> HessenbergPseudospectrum \
  ( const Matrix<Complex<Real>>& H, \
    const Matrix<Complex<Real>>& Q, \
    const Matrix<Complex<Real>>& shifts, \
    Matrix<Real>& invNorms, \
    PseudospecCtrl<Real> psCtrl ); \
  template DistMatrix<Int,VR,STAR> HessenbergPseudospectrum \
  ( const DistMatrix<Complex<Real>>& H, \
    const DistMatrix<Complex<Real>>& Q, \
    const DistMatrix<Complex<Real>,VR,STAR>& shifts, \
    DistMatrix<Real,VR,STAR>& invNorms, \
    PseudospecCtrl<Real> psCtrl ); \
  template Matrix<Int> HessenbergPseudospectrum \
  ( const Matrix<Real>& H, const Matrix<Complex<Real>>& shifts, \
    Matrix<Real>& invNorms, \
    PseudospecCtrl<Real> psCtrl ); \
  template DistMatrix<Int,VR,STAR> HessenbergPseudospectrum \
  ( const DistMatrix<Real>& H, \
    const DistMatrix<Complex<Real>,VR,STAR>& shifts, \
    DistMatrix<Real,VR,STAR>& invNorms, \
    PseudospecCtrl<Real> psCtrl ); \
  template Matrix<Int> HessenbergPseudospectrum \
  ( const Matrix<Real>& H, \
    const Matrix<Real>& Q, \
    const Matrix<Complex<Real>>& shifts, \
    Matrix<Real>& invNorms, \
    PseudospecCtrl<Real> psCtrl ); \
  template DistMatrix<Int,VR,STAR> HessenbergPseudospectrum \
  ( const DistMatrix<Real>& H, \
    const DistMatrix<Real>& Q, \
    const DistMatrix<Complex<Real>,VR,STAR>& shifts, \
    DistMatrix<Real,VR,STAR>& invNorms, \
    PseudospecCtrl<Real> psCtrl ); \
  template Matrix<Int> Pseudospectrum \
  ( const Matrix<Real>& A, const Matrix<Complex<Real>>& shifts, \
    Matrix<Real>& invNorms, \
    PseudospecCtrl<Real> psCtrl ); \
  template DistMatrix<Int,VR,STAR> Pseudospectrum \
  ( const DistMatrix<Real>& A, \
    const DistMatrix<Complex<Real>,VR,STAR>& shifts, \
    DistMatrix<Real,VR,STAR>& invNorms, \
    PseudospecCtrl<Real> psCtrl ); \
  template Matrix<Int> Pseudospectrum \
  ( const Matrix<Complex<Real>>& A, const Matrix<Complex<Real>>& shifts, \
    Matrix<Real>& invNorms, PseudospecCtrl<Real> psCtrl ); \
  template DistMatrix<Int,VR,STAR> Pseudospectrum \
  ( const DistMatrix<Complex<Real>>& A, \
    const DistMatrix<Complex<Real>,VR,STAR>& shifts, \
    DistMatrix<Real,VR,STAR>& invNorms, \
    PseudospecCtrl<Real> psCtrl ); \
  template Matrix<Int> QuasiTriangularPseudospectrum \
  ( const Matrix<Real>& U, \
    Matrix<Real>& invNormMap, \
    Complex<Real> center, Real realWidth, Real imagWidth, \
    Int realSize, Int imagSize, \
    PseudospecCtrl<Real> psCtrl ); \
  template DistMatrix<Int> QuasiTriangularPseudospectrum \
  ( const DistMatrix<Real>& U, \
    DistMatrix<Real>& invNormMap, \
    Complex<Real> center, Real realWidth, Real imagWidth, \
    Int realSize, Int imagSize, \
    PseudospecCtrl<Real> psCtrl ); \
  template Matrix<Int> QuasiTriangularPseudospectrum \
  ( const Matrix<Real>& U, \
    const Matrix<Real>& Q, \
    Matrix<Real>& invNormMap, \
    Complex<Real> center, Real realWidth, Real imagWidth, \
    Int realSize, Int imagSize, \
    PseudospecCtrl<Real> psCtrl ); \
  template DistMatrix<Int> QuasiTriangularPseudospectrum \
  ( const DistMatrix<Real>& U, \
    const DistMatrix<Real>& Q, \
    DistMatrix<Real>& invNormMap, \
    Complex<Real> center, Real realWidth, Real imagWidth, \
    Int realSize, Int imagSize, \
    PseudospecCtrl<Real> psCtrl ); \
  template Matrix<Int> QuasiTriangularPseudospectrum \
  ( const Matrix<Real>& U, \
    Matrix<Real>& invNormMap, \
    Complex<Real> center, Int realSize, Int imagSize, \
    PseudospecCtrl<Real> psCtrl ); \
  template DistMatrix<Int> QuasiTriangularPseudospectrum \
  ( const DistMatrix<Real>& U, \
    DistMatrix<Real>& invNormMap, \
    Complex<Real> center, Int realSize, Int imagSize, \
    PseudospecCtrl<Real> psCtrl ); \
  template Matrix<Int> QuasiTriangularPseudospectrum \
  ( const Matrix<Real>& U, const Matrix<Real>& Q, \
    Matrix<Real>& invNormMap, \
    Complex<Real> center, Int realSize, Int imagSize, \
    PseudospecCtrl<Real> psCtrl ); \
  template DistMatrix<Int> QuasiTriangularPseudospectrum \
  ( const DistMatrix<Real>& U, \
    const DistMatrix<Real>& Q, \
    DistMatrix<Real>& invNormMap, \
    Complex<Real> center, Int realSize, Int imagSize, \
    PseudospecCtrl<Real> psCtrl ); \
  template Matrix<Int> Pseudospectrum \
  ( const Matrix<Complex<Real>>& A, Matrix<Real>& invNormMap, \
    Complex<Real> center, Int realSize, Int imagSize, \
    PseudospecCtrl<Real> psCtrl ); \
  template DistMatrix<Int> Pseudospectrum \
  ( const DistMatrix<Complex<Real>>& A, DistMatrix<Real>& invNormMap, \
    Complex<Real> center, Int realSize, Int imagSize, \
    PseudospecCtrl<Real> psCtrl ); \
  template Matrix<Int> Pseudospectrum \
  ( const Matrix<Real>& A, Matrix<Real>& invNormMap, \
    Complex<Real> center, Int realSize, Int imagSize, \
    PseudospecCtrl<Real> psCtrl ); \
  template DistMatrix<Int> Pseudospectrum \
  ( const DistMatrix<Real>& A, DistMatrix<Real>& invNormMap, \
    Complex<Real> center, Int realSize, Int imagSize, \
    PseudospecCtrl<Real> psCtrl );

#define PROTO(F) \
  template Matrix<Int> TriangularPseudospectrum \
  ( const Matrix<F>& U, Matrix<Base<F>>& invNormMap, \
    Complex<Base<F>> center, Base<F> realWidth, Base<F> imagWidth, \
    Int realSize, Int imagSize, \
    PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() ); \
  template DistMatrix<Int> TriangularPseudospectrum \
  ( const DistMatrix<F>& U, DistMatrix<Base<F>>& invNormMap, \
    Complex<Base<F>> center, Base<F> realWidth, Base<F> imagWidth, \
    Int realSize, Int imagSize, \
    PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() ); \
  template Matrix<Int> TriangularPseudospectrum \
  ( const Matrix<F>& U, const Matrix<F>& Q, Matrix<Base<F>>& invNormMap, \
    Complex<Base<F>> center, Base<F> realWidth, Base<F> imagWidth, \
    Int realSize, Int imagSize, \
    PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() ); \
  template DistMatrix<Int> TriangularPseudospectrum \
  ( const DistMatrix<F>& U, const DistMatrix<F>& Q, \
    DistMatrix<Base<F>>& invNormMap, \
    Complex<Base<F>> center, Base<F> realWidth, Base<F> imagWidth, \
    Int realSize, Int imagSize, \
    PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() ); \
  template Matrix<Int> HessenbergPseudospectrum \
  ( const Matrix<F>& H, Matrix<Base<F>>& invNormMap, \
    Complex<Base<F>> center, Base<F> realWidth, Base<F> imagWidth, \
    Int realSize, Int imagSize, \
    PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() ); \
  template DistMatrix<Int> HessenbergPseudospectrum \
  ( const DistMatrix<F>& H, DistMatrix<Base<F>>& invNormMap, \
    Complex<Base<F>> center, Base<F> realWidth, Base<F> imagWidth, \
    Int realSize, Int imagSize, \
    PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() ); \
  template Matrix<Int> HessenbergPseudospectrum \
  ( const Matrix<F>& H, const Matrix<F>& Q, Matrix<Base<F>>& invNormMap, \
    Complex<Base<F>> center, Base<F> realWidth, Base<F> imagWidth, \
    Int realSize, Int imagSize, \
    PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() ); \
  template DistMatrix<Int> HessenbergPseudospectrum \
  ( const DistMatrix<F>& H, \
    const DistMatrix<F>& Q, \
    DistMatrix<Base<F>>& invNormMap, \
    Complex<Base<F>> center, Base<F> realWidth, Base<F> imagWidth, \
    Int realSize, Int imagSize, \
    PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() ); \
  template Matrix<Int> Pseudospectrum \
  ( const Matrix<F>& A, Matrix<Base<F>>& invNormMap, \
    Complex<Base<F>> center, Base<F> realWidth, Base<F> imagWidth, \
    Int realSize, Int imagSize, \
    PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() ); \
  template DistMatrix<Int> Pseudospectrum \
  ( const DistMatrix<F>& A, DistMatrix<Base<F>>& invNormMap, \
    Complex<Base<F>> center, Base<F> realWidth, Base<F> imagWidth, \
    Int realSize, Int imagSize, \
    PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() ); \
  template Matrix<Int> TriangularPseudospectrum \
  ( const Matrix<F>& U, Matrix<Base<F>>& invNormMap, \
    Complex<Base<F>> center, Int realSize, Int imagSize, \
    PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() ); \
  template DistMatrix<Int> TriangularPseudospectrum \
  ( const DistMatrix<F>& U, DistMatrix<Base<F>>& invNormMap, \
    Complex<Base<F>> center, Int realSize, Int imagSize, \
    PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() ); \
  template Matrix<Int> TriangularPseudospectrum \
  ( const Matrix<F>& U, const Matrix<F>& Q, Matrix<Base<F>>& invNormMap, \
    Complex<Base<F>> center, Int realSize, Int imagSize, \
    PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() ); \
  template DistMatrix<Int> TriangularPseudospectrum \
  ( const DistMatrix<F>& U, const DistMatrix<F>& Q, \
    DistMatrix<Base<F>>& invNormMap, \
    Complex<Base<F>> center, Int realSize, Int imagSize, \
    PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() ); \
  template Matrix<Int> HessenbergPseudospectrum \
  ( const Matrix<F>& H, Matrix<Base<F>>& invNormMap, \
    Complex<Base<F>> center, Int realSize, Int imagSize, \
    PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() ); \
  template DistMatrix<Int> HessenbergPseudospectrum \
  ( const DistMatrix<F>& H, DistMatrix<Base<F>>& invNormMap, \
    Complex<Base<F>> center, Int realSize, Int imagSize, \
    PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() ); \
  template Matrix<Int> HessenbergPseudospectrum \
  ( const Matrix<F>& H, const Matrix<F>& Q, Matrix<Base<F>>& invNormMap, \
    Complex<Base<F>> center, Int realSize, Int imagSize, \
    PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() ); \
  template DistMatrix<Int> HessenbergPseudospectrum \
  ( const DistMatrix<F>& H, \
    const DistMatrix<F>& Q, DistMatrix<Base<F>>& invNormMap, \
    Complex<Base<F>> center, Int realSize, Int imagSize, \
    PseudospecCtrl<Base<F>> psCtrl=PseudospecCtrl<Base<F>>() );

PROTO_REAL(float)
PROTO_REAL(double)

PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
