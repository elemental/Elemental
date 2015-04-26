/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {
namespace reg_qsd_ldl {

// TODO: Switch to returning the relative residual of the refined solution

template<typename F>
inline Int RegularizedSolveAfterNoPromote
( const SparseMatrix<F>& A,    const Matrix<Base<F>>& reg,
  const vector<Int>& invMap,   const SymmNodeInfo& info,
  const SymmFront<F>& front,         Matrix<F>& b,
  Base<F> relTol,                    Int maxRefineIts,
  bool progress )
{
    DEBUG_ONLY(
      CallStackEntry cse("reg_qsd_ldl::RegularizedSolveAfterNoPromote")
    )
    auto bOrig = b;
    const Base<F> bNorm = Nrm2( b );

    // Compute the initial guess
    // =========================
    Matrix<F> x;
    MatrixNode<F> xNodal( invMap, info, b );
    ldl::SolveAfter( info, front, xNodal );
    xNodal.Push( invMap, info, x );

    Int refineIt = 0;
    if( maxRefineIts > 0 )
    {
        Matrix<F> dx, y;
        y = x;
        DiagonalScale( LEFT, NORMAL, reg, y );
        Multiply( NORMAL, F(1), A, x, F(1), y );
        Axpy( F(-1), y, b );
        Base<F> errorNorm = Nrm2( b );
        if( progress )
            cout << "    original rel error: " << errorNorm/bNorm << endl;
 
        while( true )
        {
            if( errorNorm/bNorm <= relTol )
            {
                if( progress )
                    cout << "    " << errorNorm/bNorm << " <= " << relTol
                         << endl;
                break;
            }

            // Compute the proposed update to the solution
            // -------------------------------------------
            xNodal.Pull( invMap, info, b );
            ldl::SolveAfter( info, front, xNodal );
            xNodal.Push( invMap, info, dx );
            Axpy( F(1), dx, x );

            // Check the new residual
            // ----------------------
            b = bOrig;
            y = x;
            DiagonalScale( LEFT, NORMAL, reg, y );
            Multiply( NORMAL, F(1), A, x, F(1), y );
            Axpy( F(-1), y, b );
            auto newErrorNorm = Nrm2( b );
            if( progress )
                cout << "    refined rel error: " << newErrorNorm/bNorm << endl;

            errorNorm = newErrorNorm;
            ++refineIt;
            if( refineIt >= maxRefineIts )
                RuntimeError("Iterative refinement did not converge in time"); 
        }
    }
    // Store the final result
    // ======================
    b = x;
    return refineIt;
}

// A will be the original sparse matrix for this system, 
// but the applied regularization is of the form diag(dR)*diag(reg)*diag(dC),
// and the sparse-direct factorization is of inv(diag(dR)) A inv(diag(dC)).
template<typename F>
inline Int RegularizedSolveAfterNoPromote
( const SparseMatrix<F>& A,    const Matrix<Base<F>>& reg,
  const Matrix<Base<F>>& d, 
  const vector<Int>& invMap,   const SymmNodeInfo& info,
  const SymmFront<F>& front,         Matrix<F>& b,
  Base<F> relTol,                    Int maxRefineIts,
  bool progress )
{
    DEBUG_ONLY(
      CallStackEntry cse("reg_qsd_ldl::RegularizedSolveAfterNoPromote")
    )
    auto bOrig = b;
    const Base<F> bNorm = Nrm2( b );

    // Compute the initial guess
    // =========================
    Matrix<F> x;
    DiagonalSolve( LEFT, NORMAL, d, b );
    MatrixNode<F> xNodal( invMap, info, b );
    ldl::SolveAfter( info, front, xNodal );
    xNodal.Push( invMap, info, x );
    DiagonalSolve( LEFT, NORMAL, d, x );

    Int refineIt = 0;
    if( maxRefineIts > 0 )
    {
        Matrix<F> dx, y;
        y = x;
        DiagonalScale( LEFT, NORMAL, d, y );
        DiagonalScale( LEFT, NORMAL, reg, y );
        DiagonalScale( LEFT, NORMAL, d, y );
        Multiply( NORMAL, F(1), A, x, F(1), y );
        Axpy( F(-1), y, b );
        Base<F> errorNorm = Nrm2( b );
        if( progress )
            cout << "    original rel error: " << errorNorm/bNorm << endl;
 
        while( true )
        {
            if( errorNorm/bNorm <= relTol )
            {
                if( progress )
                    cout << "    " << errorNorm/bNorm << " <= " << relTol
                         << endl;
                break;
            }

            // Compute the proposed update to the solution
            // -------------------------------------------
            DiagonalSolve( LEFT, NORMAL, d, b );
            xNodal.Pull( invMap, info, b );
            ldl::SolveAfter( info, front, xNodal );
            xNodal.Push( invMap, info, dx );
            DiagonalSolve( LEFT, NORMAL, d, dx );
            Axpy( F(1), dx, x );

            // Check the new residual
            // ----------------------
            b = bOrig;
            y = x;
            DiagonalScale( LEFT, NORMAL, d, y );
            DiagonalScale( LEFT, NORMAL, reg, y );
            DiagonalScale( LEFT, NORMAL, d, y );
            Multiply( NORMAL, F(1), A, x, F(1), y );
            Axpy( F(-1), y, b );
            auto newErrorNorm = Nrm2( b );
            if( progress )
                cout << "    refined rel error: " << newErrorNorm/bNorm << endl;

            errorNorm = newErrorNorm;
            ++refineIt;
            if( refineIt >= maxRefineIts )
                RuntimeError("Iterative refinement did not converge in time"); 
        }
    }
    // Store the final result
    // ======================
    b = x;
    return refineIt;
}

template<typename F>
inline Int RegularizedSolveAfterPromote
( const SparseMatrix<F>& A,    const Matrix<Base<F>>& reg,
  const vector<Int>& invMap,   const SymmNodeInfo& info,
  const SymmFront<F>& front,         Matrix<F>& b,
  Base<F> relTol,                    Int maxRefineIts,
  bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("reg_qsd_ldl::RegularizedSolveAfterPromote"))
    typedef Base<F> Real;
    typedef Promote<Real> RealProm;
    typedef Promote<F> FProm;

    Matrix<FProm> bProm, bOrigProm;
    Copy( b, bProm );
    Copy( b, bOrigProm );
    const RealProm bNorm = Nrm2( bOrigProm );

    Matrix<RealProm> regProm;
    Copy( reg, regProm );

    // Compute the initial guess
    // =========================
    MatrixNode<F> xNodal( invMap, info, b );
    ldl::SolveAfter( info, front, xNodal );
    xNodal.Push( invMap, info, b );
    Matrix<FProm> xProm;
    Copy( b, xProm );

    SparseMatrix<FProm> AProm;
    Copy( A, AProm );

    Int refineIt = 0;
    if( maxRefineIts > 0 )
    {
        Matrix<FProm> dxProm, yProm;
        yProm = xProm;
        DiagonalScale( LEFT, NORMAL, regProm, yProm );
        Multiply( NORMAL, FProm(1), AProm, xProm, FProm(1), yProm );
        Axpy( FProm(-1), yProm, bProm );
        auto errorNorm = Nrm2( bProm );
        if( progress )
            cout << "    original rel error: " << errorNorm/bNorm << endl;
 
        while( true )
        {
            if( errorNorm/bNorm <= relTol )
            {
                if( progress )
                    cout << "    " << errorNorm/bNorm << " <= " << relTol
                         << endl;
                break;
            }

            // Compute the proposed update to the solution
            // -------------------------------------------
            Copy( bProm, b );
            xNodal.Pull( invMap, info, b );
            ldl::SolveAfter( info, front, xNodal );
            xNodal.Push( invMap, info, b );
            Copy( b, dxProm );
            Axpy( FProm(1), dxProm, xProm );

            // Check the new residual
            // ----------------------
            bProm = bOrigProm;
            yProm = xProm;
            DiagonalScale( LEFT, NORMAL, regProm, yProm );
            Multiply( NORMAL, FProm(1), AProm, xProm, FProm(1), yProm );
            Axpy( FProm(-1), yProm, bProm );
            auto newErrorNorm = Nrm2( bProm );
            if( progress )
                cout << "    refined rel error: " << newErrorNorm/bNorm << endl;

            errorNorm = newErrorNorm;
            ++refineIt;
            if( refineIt >= maxRefineIts )
                RuntimeError("Iterative refinement did not converge in time"); 
        }
    }
    // Store the final result
    // ======================
    Copy( xProm, b );
    return refineIt;
}

template<typename F>
inline Int RegularizedSolveAfterPromote
( const SparseMatrix<F>& A,    const Matrix<Base<F>>& reg,
  const Matrix<Base<F>>& d, 
  const vector<Int>& invMap,   const SymmNodeInfo& info,
  const SymmFront<F>& front,         Matrix<F>& b,
  Base<F> relTol,                    Int maxRefineIts,
  bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("reg_qsd_ldl::RegularizedSolveAfterPromote"))
    typedef Base<F> Real;
    typedef Promote<Real> RealProm;
    typedef Promote<F> FProm;

    Matrix<FProm> bProm, bOrigProm;
    Copy( b, bProm );
    Copy( b, bOrigProm );
    const RealProm bNorm = Nrm2( bOrigProm );

    Matrix<RealProm> dProm;
    Copy( d, dProm );

    Matrix<RealProm> regProm;
    Copy( reg, regProm );

    // Compute the initial guess
    // =========================
    DiagonalSolve( LEFT, NORMAL, d, b );
    MatrixNode<F> xNodal( invMap, info, b );
    ldl::SolveAfter( info, front, xNodal );
    xNodal.Push( invMap, info, b );
    DiagonalSolve( LEFT, NORMAL, d, b );

    Matrix<FProm> xProm;
    Copy( b, xProm );

    SparseMatrix<FProm> AProm;
    Copy( A, AProm );

    Int refineIt = 0;
    if( maxRefineIts > 0 )
    {
        Matrix<FProm> dxProm, yProm;
        yProm = xProm;
        DiagonalScale( LEFT, NORMAL, dProm, yProm );
        DiagonalScale( LEFT, NORMAL, regProm, yProm );
        DiagonalScale( LEFT, NORMAL, dProm, yProm );
        Multiply( NORMAL, FProm(1), AProm, xProm, FProm(1), yProm );
        Axpy( FProm(-1), yProm, bProm );
        auto errorNorm = Nrm2( bProm );
        if( progress )
            cout << "    original rel error: " << errorNorm/bNorm << endl;
 
        while( true )
        {
            if( errorNorm/bNorm <= relTol )
            {
                if( progress )
                    cout << "    " << errorNorm/bNorm << " <= " << relTol
                         << endl;
                break;
            }

            // Compute the proposed update to the solution
            // -------------------------------------------
            Copy( bProm, b );
            DiagonalSolve( LEFT, NORMAL, d, b );
            xNodal.Pull( invMap, info, b );
            ldl::SolveAfter( info, front, xNodal );
            xNodal.Push( invMap, info, b );
            DiagonalSolve( LEFT, NORMAL, d, b );
            Copy( b, dxProm );
            Axpy( FProm(1), dxProm, xProm );

            // Check the new residual
            // ----------------------
            bProm = bOrigProm;
            yProm = xProm;
            DiagonalScale( LEFT, NORMAL, dProm, yProm );
            DiagonalScale( LEFT, NORMAL, regProm, yProm );
            DiagonalScale( LEFT, NORMAL, dProm, yProm );
            Multiply( NORMAL, FProm(1), AProm, xProm, FProm(1), yProm );
            Axpy( FProm(-1), yProm, bProm );
            auto newErrorNorm = Nrm2( bProm );
            if( progress )
                cout << "    refined rel error: " << newErrorNorm/bNorm << endl;

            errorNorm = newErrorNorm;
            ++refineIt;
            if( refineIt >= maxRefineIts )
                RuntimeError("Iterative refinement did not converge in time"); 
        }
    }
    // Store the final result
    // ======================
    Copy( xProm, b );
    return refineIt;
}

template<typename F>
Int RegularizedSolveAfter
( const SparseMatrix<F>& A,    const Matrix<Base<F>>& reg,
  const vector<Int>& invMap,   const SymmNodeInfo& info,
  const SymmFront<F>& front,         Matrix<F>& b,
  Base<F> relTol,                    Int maxRefineIts,
  bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("reg_qsd_ldl::RegularizedSolveAfter"))
    return RegularizedSolveAfterPromote
           ( A, reg, invMap, info, front, b, relTol, maxRefineIts, progress );
}

template<typename F>
Int RegularizedSolveAfter
( const SparseMatrix<F>& A,    const Matrix<Base<F>>& reg,
  const Matrix<Base<F>>& d, 
  const vector<Int>& invMap,   const SymmNodeInfo& info,
  const SymmFront<F>& front,         Matrix<F>& b,
  Base<F> relTol,                    Int maxRefineIts,
  bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("reg_qsd_ldl::RegularizedSolveAfter"))
    return RegularizedSolveAfterPromote
           ( A, reg, d, invMap, info, front, 
             b, relTol, maxRefineIts, progress );
}

template<typename F>
inline Int RegularizedSolveAfterNoPromote
( const DistSparseMatrix<F>& A,      const DistMultiVec<Base<F>>& reg,
  const DistMap& invMap,             const DistSymmNodeInfo& info,
  const DistSymmFront<F>& front,           DistMultiVec<F>& b,
  Base<F> relTol,                          Int maxRefineIts,
  bool progress )
{
    DEBUG_ONLY(
      CallStackEntry cse("reg_qsd_ldl::RegularizedSolveAfterNoPromote")
    )
    mpi::Comm comm = b.Comm();
    const Int commRank = mpi::Rank(comm);

    DistMultiVec<F> bOrig(comm);
    bOrig = b;
    const Base<F> bNorm = Nrm2( b );

    // Compute the initial guess
    // =========================
    DistMultiVec<F> x(comm);
    DistMultiVecNode<F> xNodal( invMap, info, b );
    ldl::SolveAfter( info, front, xNodal );
    xNodal.Push( invMap, info, x );

    Int refineIt = 0;
    if( maxRefineIts > 0 )
    {
        DistMultiVec<F> dx(comm), y(comm);
        y = x;
        DiagonalScale( LEFT, NORMAL, reg, y );
        Multiply( NORMAL, F(1), A, x, F(1), y );
        Axpy( F(-1), y, b );
        Base<F> errorNorm = Nrm2( b );
        if( progress && commRank == 0 )
            cout << "    original rel error: " << errorNorm/bNorm << endl;
        while( true )
        {
            if( errorNorm/bNorm <= relTol )
            {
                if( progress && commRank == 0 )
                    cout << "    " << errorNorm/bNorm << " <= " << relTol
                         << endl;
                break;
            }

            // Compute the proposed update to the solution
            // -------------------------------------------
            xNodal.Pull( invMap, info, b );
            ldl::SolveAfter( info, front, xNodal );
            xNodal.Push( invMap, info, dx );
            Axpy( F(1), dx, x );

            // Compute the new residual
            // ------------------------
            b = bOrig;
            y = x;
            DiagonalScale( LEFT, NORMAL, reg, y );
            Multiply( NORMAL, F(1), A, x, F(1), y );
            Axpy( F(-1), y, b );
            Base<F> newErrorNorm = Nrm2( b );
            if( progress && commRank == 0 )
                cout << "    refined rel error: " << newErrorNorm/bNorm << endl;

            errorNorm = newErrorNorm;
            ++refineIt;
            if( refineIt >= maxRefineIts )
                RuntimeError("Iterative refinement did not converge in time"); 
        }
    }
    // Store the final result
    // ======================
    b = x;
    return refineIt;
}

template<typename F>
inline Int RegularizedSolveAfterNoPromote
( const DistSparseMatrix<F>& A,      const DistMultiVec<Base<F>>& reg,
  const DistMultiVec<Base<F>>& d, 
  const DistMap& invMap,             const DistSymmNodeInfo& info,
  const DistSymmFront<F>& front,           DistMultiVec<F>& b,
  Base<F> relTol,                          Int maxRefineIts,
  bool progress )
{
    DEBUG_ONLY(
      CallStackEntry cse("reg_qsd_ldl::RegularizedSolveAfterNoPromote")
    )
    mpi::Comm comm = b.Comm();
    const Int commRank = mpi::Rank(comm);

    DistMultiVec<F> bOrig(comm);
    bOrig = b;
    const Base<F> bNorm = Nrm2( b );

    // Compute the initial guess
    // =========================
    DistMultiVec<F> x(comm);
    DiagonalSolve( LEFT, NORMAL, d, b );
    DistMultiVecNode<F> xNodal( invMap, info, b );
    ldl::SolveAfter( info, front, xNodal );
    xNodal.Push( invMap, info, x );
    DiagonalSolve( LEFT, NORMAL, d, x );

    Int refineIt = 0;
    if( maxRefineIts > 0 )
    {
        DistMultiVec<F> dx(comm), y(comm);
        y = x;
        DiagonalScale( LEFT, NORMAL, d, y );
        DiagonalScale( LEFT, NORMAL, reg, y );
        DiagonalScale( LEFT, NORMAL, d, y );
        Multiply( NORMAL, F(1), A, x, F(1), y );
        Axpy( F(-1), y, b );
        Base<F> errorNorm = Nrm2( b );
        if( progress && commRank == 0 )
            cout << "    original rel error: " << errorNorm/bNorm << endl;
        while( true )
        {
            if( errorNorm/bNorm <= relTol )
            {
                if( progress && commRank == 0 )
                    cout << "    " << errorNorm/bNorm << " <= " << relTol
                         << endl;
                break;
            }

            // Compute the proposed update to the solution
            // -------------------------------------------
            DiagonalSolve( LEFT, NORMAL, d, b );
            xNodal.Pull( invMap, info, b );
            ldl::SolveAfter( info, front, xNodal );
            xNodal.Push( invMap, info, dx );
            DiagonalSolve( LEFT, NORMAL, d, dx );
            Axpy( F(1), dx, x );

            // Compute the new residual
            // ------------------------
            b = bOrig;
            y = x;
            DiagonalScale( LEFT, NORMAL, d, y );
            DiagonalScale( LEFT, NORMAL, reg, y );
            DiagonalScale( LEFT, NORMAL, d, y );
            Multiply( NORMAL, F(1), A, x, F(1), y );
            Axpy( F(-1), y, b );
            Base<F> newErrorNorm = Nrm2( b );
            if( progress && commRank == 0 )
                cout << "    refined rel error: " << newErrorNorm/bNorm << endl;

            errorNorm = newErrorNorm;
            ++refineIt;
            if( refineIt >= maxRefineIts )
                RuntimeError("Iterative refinement did not converge in time"); 
        }
    }
    // Store the final result
    // ======================
    b = x;
    return refineIt;
}

template<typename F>
inline Int RegularizedSolveAfterPromote
( const DistSparseMatrix<F>& A,      const DistMultiVec<Base<F>>& reg,
  const DistMap& invMap,             const DistSymmNodeInfo& info,
  const DistSymmFront<F>& front,           DistMultiVec<F>& b,
  Base<F> relTol,                          Int maxRefineIts,
  bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("reg_qsd_ldl::RegularizedSolveAfterPromote"))
    typedef Base<F> Real;
    typedef Promote<Real> RealProm;
    typedef Promote<F> FProm;

    mpi::Comm comm = b.Comm();
    const Int commRank = mpi::Rank(comm);

    DistMultiVec<FProm> bProm(comm), bOrigProm(comm);
    Copy( b, bProm ); 
    Copy( b, bOrigProm );
    const auto bNorm = Nrm2( bProm );

    DistMultiVec<RealProm> regProm(comm);
    Copy( reg, regProm );

    // Compute the initial guess
    // =========================
    DistMultiVecNode<F> xNodal( invMap, info, b );
    ldl::SolveAfter( info, front, xNodal );
    xNodal.Push( invMap, info, b );
    DistMultiVec<FProm> xProm(comm);
    Copy( b, xProm );

    DistSparseMatrix<FProm> AProm(comm);
    Copy( A, AProm );

    Int refineIt = 0;
    if( maxRefineIts > 0 )
    {
        DistMultiVec<FProm> dxProm(comm), yProm(comm);
        yProm = xProm;
        DiagonalScale( LEFT, NORMAL, regProm, yProm );
        Multiply( NORMAL, FProm(1), AProm, xProm, FProm(1), yProm );
        Axpy( FProm(-1), yProm, bProm );
        auto errorNorm = Nrm2( bProm );
        if( progress && commRank == 0 )
            cout << "    original rel error: " << errorNorm/bNorm << endl;
        while( true )
        {
            if( errorNorm/bNorm <= relTol )
            {
                if( progress && commRank == 0 )
                    cout << "    " << errorNorm/bNorm << " <= " << relTol
                         << endl;
                break;
            }

            // Compute the proposed update to the solution
            // -------------------------------------------
            Copy( bProm, b );
            xNodal.Pull( invMap, info, b );
            ldl::SolveAfter( info, front, xNodal );
            xNodal.Push( invMap, info, b );
            Copy( b, dxProm );
            Axpy( FProm(1), dxProm, xProm );

            // Check the new residual
            // ----------------------
            bProm = bOrigProm;
            yProm = xProm;
            DiagonalScale( LEFT, NORMAL, regProm, yProm );
            Multiply( NORMAL, FProm(1), AProm, xProm, FProm(1), yProm );
            Axpy( FProm(-1), yProm, bProm );
            auto newErrorNorm = Nrm2( bProm );
            if( progress && commRank == 0 )
                cout << "    refined rel error: " << newErrorNorm/bNorm << endl;

            errorNorm = newErrorNorm;
            ++refineIt;
            if( refineIt >= maxRefineIts )
                RuntimeError("Iterative refinement did not converge in time"); 
        }
    }
    // Store the final result
    // ======================
    Copy( xProm, b );
    return refineIt;
}

template<typename F>
inline Int RegularizedSolveAfterPromote
( const DistSparseMatrix<F>& A,      const DistMultiVec<Base<F>>& reg,
  const DistMultiVec<Base<F>>& d, 
  const DistMap& invMap,             const DistSymmNodeInfo& info,
  const DistSymmFront<F>& front,           DistMultiVec<F>& b,
  Base<F> relTol,                          Int maxRefineIts,
  bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("reg_qsd_ldl::RegularizedSolveAfterPromote"))
    typedef Base<F> Real;
    typedef Promote<Real> RealProm;
    typedef Promote<F> FProm;

    mpi::Comm comm = b.Comm();
    const Int commRank = mpi::Rank(comm);

    DistMultiVec<FProm> bProm(comm), bOrigProm(comm);
    Copy( b, bProm ); 
    Copy( b, bOrigProm );
    const auto bNorm = Nrm2( bProm );

    DistMultiVec<RealProm> dProm(comm);
    Copy( d, dProm );

    DistMultiVec<RealProm> regProm(comm);
    Copy( reg, regProm );

    // Compute the initial guess
    // =========================
    DiagonalSolve( LEFT, NORMAL, d, b );
    DistMultiVecNode<F> xNodal( invMap, info, b );
    ldl::SolveAfter( info, front, xNodal );
    xNodal.Push( invMap, info, b );
    DiagonalSolve( LEFT, NORMAL, d, b );

    DistMultiVec<FProm> xProm(comm);
    Copy( b, xProm );

    DistSparseMatrix<FProm> AProm(comm);
    Copy( A, AProm );

    Int refineIt = 0;
    if( maxRefineIts > 0 )
    {
        DistMultiVec<FProm> dxProm(comm), yProm(comm);
        yProm = xProm;
        DiagonalScale( LEFT, NORMAL, dProm, yProm );
        DiagonalScale( LEFT, NORMAL, regProm, yProm );
        DiagonalScale( LEFT, NORMAL, dProm, yProm );
        Multiply( NORMAL, FProm(1), AProm, xProm, FProm(1), yProm );
        Axpy( FProm(-1), yProm, bProm );
        auto errorNorm = Nrm2( bProm );
        if( progress && commRank == 0 )
            cout << "    original rel error: " << errorNorm/bNorm << endl;
        while( true )
        {
            if( errorNorm/bNorm <= relTol )
            {
                if( progress && commRank == 0 )
                    cout << "    " << errorNorm/bNorm << " <= " << relTol
                         << endl;
                break;
            }

            // Compute the proposed update to the solution
            // -------------------------------------------
            Copy( bProm, b );
            DiagonalSolve( LEFT, NORMAL, d, b );
            xNodal.Pull( invMap, info, b );
            ldl::SolveAfter( info, front, xNodal );
            xNodal.Push( invMap, info, b );
            DiagonalSolve( LEFT, NORMAL, d, b );
            Copy( b, dxProm );
            Axpy( FProm(1), dxProm, xProm );

            // Check the new residual
            // ----------------------
            bProm = bOrigProm;
            yProm = xProm;
            DiagonalScale( LEFT, NORMAL, dProm, yProm );
            DiagonalScale( LEFT, NORMAL, regProm, yProm );
            DiagonalScale( LEFT, NORMAL, dProm, yProm );
            Multiply( NORMAL, FProm(1), AProm, xProm, FProm(1), yProm );
            Axpy( FProm(-1), yProm, bProm );
            auto newErrorNorm = Nrm2( bProm );
            if( progress && commRank == 0 )
                cout << "    refined rel error: " << newErrorNorm/bNorm << endl;

            errorNorm = newErrorNorm;
            ++refineIt;
            if( refineIt >= maxRefineIts )
                RuntimeError("Iterative refinement did not converge in time"); 
        }
    }
    // Store the final result
    // ======================
    Copy( xProm, b );
    return refineIt;
}

template<typename F>
Int RegularizedSolveAfter
( const DistSparseMatrix<F>& A,      const DistMultiVec<Base<F>>& reg,
  const DistMap& invMap,             const DistSymmNodeInfo& info,
  const DistSymmFront<F>& front,           DistMultiVec<F>& b,
  Base<F> relTol,                          Int maxRefineIts,
  bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("reg_qsd_ldl::RegularizedSolveAfter"))
    return RegularizedSolveAfterPromote
    ( A, reg, invMap, info, front, b, relTol, maxRefineIts, progress );
}

template<typename F>
Int RegularizedSolveAfter
( const DistSparseMatrix<F>& A,    const DistMultiVec<Base<F>>& reg,
  const DistMultiVec<Base<F>>& d, 
  const DistMap& invMap,           const DistSymmNodeInfo& info,
  const DistSymmFront<F>& front,         DistMultiVec<F>& b,
  Base<F> relTol,                        Int maxRefineIts,
  bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("reg_qsd_ldl::RegularizedSolveAfter"))
    return RegularizedSolveAfterPromote
    ( A, reg, d, invMap, info, front, b, relTol, maxRefineIts, progress );
}

template<typename F>
Int IRSolveAfter
( const SparseMatrix<F>& A,  const Matrix<Base<F>>& reg,
  const vector<Int>& invMap, const SymmNodeInfo& info,
  const SymmFront<F>& front,       Matrix<F>& b,
  Base<F> relTol,                  Int maxRefineIts,
  bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("reg_qsd_ldl::IRSolveAfter"))
    auto bOrig = b;
    const Base<F> bNorm = Nrm2( b );

    // Compute the initial guess
    // =========================
    auto x = b;
    RegularizedSolveAfter
    ( A, reg, invMap, info, front, x, relTol, maxRefineIts, progress );

    Int refineIt = 0;
    if( maxRefineIts > 0 )
    {
        Matrix<F> dx;
        Multiply( NORMAL, F(-1), A, x, F(1), b );
        Base<F> errorNorm = Nrm2( b );
        if( progress )
            cout << "    original rel error: " << errorNorm/bNorm << endl;
        while( true )
        {
            if( errorNorm/bNorm <= relTol )
            {
                if( progress )
                    cout << "    " << errorNorm/bNorm << " <= " << relTol
                         << endl;
                break;
            }

            // Compute the proposed update to the solution
            // -------------------------------------------
            dx = b;
            RegularizedSolveAfter
            ( A, reg, invMap, info, front, dx, relTol, maxRefineIts, progress );
            Axpy( F(1), dx, x );

            // Compute the new residual
            // ------------------------
            b = bOrig;
            Multiply( NORMAL, F(-1), A, x, F(1), b );
            Base<F> newErrorNorm = Nrm2( b );
            if( progress )
                cout << "    refined rel error: " << newErrorNorm/bNorm << endl;

            errorNorm = newErrorNorm;
            ++refineIt;
            if( refineIt >= maxRefineIts )
                RuntimeError("Iterative refinement did not converge in time"); 
        }
    }
    // Store the final result
    // ======================
    b = x;
    return refineIt;
}

template<typename F>
Int IRSolveAfter
( const SparseMatrix<F>& A,  const Matrix<Base<F>>& reg,
  const Matrix<Base<F>>& d,
  const vector<Int>& invMap, const SymmNodeInfo& info,
  const SymmFront<F>& front,       Matrix<F>& b,
  Base<F> relTol,                  Int maxRefineIts,
  bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("reg_qsd_ldl::IRSolveAfter"))
    auto bOrig = b;
    const Base<F> bNorm = Nrm2( b );

    // Compute the initial guess
    // =========================
    auto x = b;
    RegularizedSolveAfter
    ( A, reg, d, invMap, info, front, x, relTol, maxRefineIts, progress );

    Int refineIt = 0;
    if( maxRefineIts > 0 )
    {
        Matrix<F> dx;
        Multiply( NORMAL, F(-1), A, x, F(1), b );
        Base<F> errorNorm = Nrm2( b );
        if( progress )
            cout << "    original rel error: " << errorNorm/bNorm << endl;
        while( true )
        {
            if( errorNorm/bNorm <= relTol )
            {
                if( progress )
                    cout << "    " << errorNorm/bNorm << " <= " << relTol
                         << endl;
                break;
            }

            // Compute the proposed update to the solution
            // -------------------------------------------
            dx = b;
            RegularizedSolveAfter
            ( A, reg, d, invMap, info, front, 
              dx, relTol, maxRefineIts, progress );
            Axpy( F(1), dx, x );

            // Compute the new residual
            // ------------------------
            b = bOrig;
            Multiply( NORMAL, F(-1), A, x, F(1), b );
            Base<F> newErrorNorm = Nrm2( b );
            if( progress )
                cout << "    refined rel error: " << newErrorNorm/bNorm << endl;

            errorNorm = newErrorNorm;
            ++refineIt;
            if( refineIt >= maxRefineIts )
                RuntimeError("Iterative refinement did not converge in time"); 
        }
    }
    // Store the final result
    // ======================
    b = x;
    return refineIt;
}

template<typename F>
Int IRSolveAfter
( const DistSparseMatrix<F>& A,  const DistMultiVec<Base<F>>& reg,
  const DistMap& invMap,         const DistSymmNodeInfo& info,
  const DistSymmFront<F>& front,       DistMultiVec<F>& b,
  Base<F> relTol,                      Int maxRefineIts,
  bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("reg_qsd_ldl::IRSolveAfter"))
    mpi::Comm comm = b.Comm();
    const Int commRank = mpi::Rank(comm);

    DistMultiVec<F> bOrig(comm);
    bOrig = b;
    const Base<F> bNorm = Nrm2( b );

    // Compute the initial guess
    // =========================
    DistMultiVec<F> x(comm);
    x = b;
    RegularizedSolveAfter
    ( A, reg, invMap, info, front, x, relTol, maxRefineIts, progress );

    Int refineIt = 0;
    if( maxRefineIts > 0 )
    {
        DistMultiVec<F> dx(comm);
        Multiply( NORMAL, F(-1), A, x, F(1), b );
        Base<F> errorNorm = Nrm2( b );
        if( progress && commRank == 0 )
            cout << "    original rel error: " << errorNorm/bNorm << endl;
        while( true )
        {
            if( errorNorm/bNorm <= relTol )
            {
                if( progress && commRank == 0 )
                    cout << "    " << errorNorm/bNorm << " <= " << relTol
                         << endl;
                break;
            }

            // Compute the proposed update to the solution
            // -------------------------------------------
            dx = b;
            RegularizedSolveAfter
            ( A, reg, invMap, info, front, dx, relTol, maxRefineIts, progress );
            Axpy( F(1), dx, x );

            // If the proposed update lowers the residual, accept it
            // -----------------------------------------------------
            b = bOrig;
            Multiply( NORMAL, F(-1), A, x, F(1), b );
            Base<F> newErrorNorm = Nrm2( b );
            if( progress && commRank == 0 )
                cout << "    refined rel error: " << newErrorNorm/bNorm << endl;

            errorNorm = newErrorNorm;
            ++refineIt;
            if( refineIt >= maxRefineIts )
                RuntimeError("Refinement did not converge in time");
        }
    }
    // Store the final result
    // ======================
    b = x;
    return refineIt;
}

template<typename F>
Int IRSolveAfter
( const DistSparseMatrix<F>& A,    const DistMultiVec<Base<F>>& reg,
  const DistMultiVec<Base<F>>& d,
  const DistMap& invMap,           const DistSymmNodeInfo& info,
  const DistSymmFront<F>& front,         DistMultiVec<F>& b,
  Base<F> relTol,                        Int maxRefineIts,
  bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("reg_qsd_ldl::IRSolveAfter"))
    mpi::Comm comm = b.Comm();
    const Int commRank = mpi::Rank(comm);

    DistMultiVec<F> bOrig(comm);
    bOrig = b;
    const Base<F> bNorm = Nrm2( b );

    // Compute the initial guess
    // =========================
    DistMultiVec<F> x(comm);
    x = b;
    RegularizedSolveAfter
    ( A, reg, d, invMap, info, front, x, relTol, maxRefineIts, progress );

    Int refineIt = 0;
    if( maxRefineIts > 0 )
    {
        DistMultiVec<F> dx(comm);
        Multiply( NORMAL, F(-1), A, x, F(1), b );
        Base<F> errorNorm = Nrm2( b );
        if( progress && commRank == 0 )
            cout << "    original rel error: " << errorNorm/bNorm << endl;
        while( true )
        {
            if( errorNorm/bNorm <= relTol )
            {
                if( progress && commRank == 0 )
                    cout << "    " << errorNorm/bNorm << " <= " << relTol
                         << endl;
                break;
            }

            // Compute the proposed update to the solution
            // -------------------------------------------
            dx = b;
            RegularizedSolveAfter
            ( A, reg, d, invMap, info, front, 
              dx, relTol, maxRefineIts, progress );
            Axpy( F(1), dx, x );

            // If the proposed update lowers the residual, accept it
            // -----------------------------------------------------
            b = bOrig;
            Multiply( NORMAL, F(-1), A, x, F(1), b );
            Base<F> newErrorNorm = Nrm2( b );
            if( progress && commRank == 0 )
                cout << "    refined rel error: " << newErrorNorm/bNorm << endl;

            errorNorm = newErrorNorm;
            ++refineIt;
            if( refineIt >= maxRefineIts )
                RuntimeError("Refinement did not converge in time");
        }
    }
    // Store the final result
    // ======================
    b = x;
    return refineIt;
}

template<typename F>
Int LGMRESSolveAfter
( const SparseMatrix<F>& A,   const Matrix<Base<F>>& reg,
  const vector<Int>& invMap,  const SymmNodeInfo& info,
  const SymmFront<F>& front,        Matrix<F>& b,
  Base<F> relTol,                   Int k,
  Base<F> relTolRefine,             Int maxRefineIts,
  bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("reg_qsd_ldl::LGMRESSolveAfter"))
    typedef Base<F> Real;
    const Int n = A.Height();

    // x := 0
    // ======
    Matrix<F> x;
    Zeros( x, n, 1 );

    // w := b (= b - A x_0)
    // ====================
    Matrix<F> w; 
    w = b;
    const Real origResidNorm = Nrm2( w );
    if( origResidNorm == Real(0) )
        return 0;

    Int iter=0;
    Int maxLargeRefines=0;
    bool converged = false;
    Matrix<Real> cs;
    Matrix<F> sn, H, t;
    Matrix<F> x0, q, V;
    while( !converged )
    {
        if( progress )
            cout << "  Starting GMRES iteration " << iter << endl;
        Zeros( cs, k, 1 );
        Zeros( sn, k, 1 );
        Zeros( H,  k, k );
        Zeros( V, n, k );
        
        // x0 := x
        // =======
        x0 = x;

        // w := inv(M) w
        // =============
        Int refineIts = RegularizedSolveAfter
        ( A, reg, invMap, info, front, w, 
          relTolRefine, maxRefineIts, progress );
        maxLargeRefines = Max( refineIts, maxLargeRefines );

        // beta := || w ||_2
        // =================
        const Real beta = Nrm2( w );

        // v0 := w / beta
        // ==============
        auto v0 = V( ALL, IR(0) );
        v0 = w;
        Scale( Real(1)/beta, v0 ); 

        // t := beta e_0
        // =============
        Zeros( t, k+1, 1 );
        t.Set( 0, 0, beta );

        // Run one round of GMRES(k)
        // =========================
        for( Int j=0; j<k; ++j )
        {
            // w := A v_j
            // ----------
            Multiply( NORMAL, F(1), A, V(ALL,IR(j)), F(0), w );

            // w := inv(M) w
            // -------------
            Int refineIts = RegularizedSolveAfter
            ( A, reg, invMap, info, front, w, 
              relTolRefine, maxRefineIts, progress );
            maxLargeRefines = Max( refineIts, maxLargeRefines );

            // Run the j'th step of Arnoldi
            // ----------------------------
            for( Int i=0; i<=j; ++i )
            {
                // H(i,j) := v_i' w
                // ^^^^^^^^^^^^^^^^
                auto vi = V( ALL, IR(i) );
                H.Set( i, j, Dot(vi,w) ); 
              
                // w := w - H(i,j) v_i
                // ^^^^^^^^^^^^^^^^^^^
                Axpy( -H.Get(i,j), vi, w );
            }
            const Real delta = Nrm2( w );
            if( std::isnan(delta) )
                RuntimeError("Arnoldi step produced a NaN");
            if( delta == Real(0) )
                RuntimeError("Lucky breakdowns not yet handled");
            if( j+1 != k )
            {
                // v_{j+1} := w / delta
                // ^^^^^^^^^^^^^^^^^^^^^^^^^^
                auto vjp1 = V( ALL, IR(j+1) );
                vjp1 = w;
                Scale( Real(1)/delta, vjp1 );
            }

            // Apply existing rotations to the new column of H
            // -----------------------------------------------
            for( Int i=0; i<j; ++i )
            {
                const Real c = cs.Get(i,0);
                const F s = sn.Get(i,0);
                const F sConj = Conj(s);
                const F eta_i_j = H.Get(i,j);
                const F eta_ip1_j = H.Get(i+1,j);
                H.Set( i,   j,  c    *eta_i_j + s*eta_ip1_j );
                H.Set( i+1, j, -sConj*eta_i_j + c*eta_ip1_j );
            }

            // Generate and apply a new rotation to both H and the rotated
            // beta*e_0 vector, t, then solve the minimum residual problem
            // -----------------------------------------------------------
            const F eta_j_j = H.Get(j,j);
            const F eta_jp1_j = delta;
            if( std::isnan(RealPart(eta_j_j))   || 
                std::isnan(ImagPart(eta_j_j))   ||
                std::isnan(RealPart(eta_jp1_j)) || 
                std::isnan(ImagPart(eta_jp1_j)) )
                RuntimeError("Either H(j,j) or H(j+1,j) was NaN");
            Real c;
            F s;
            F rho = lapack::Givens( eta_j_j, eta_jp1_j, &c, &s );
            if( std::isnan(c) || 
                std::isnan(RealPart(s)) || std::isnan(ImagPart(s)) || 
                std::isnan(RealPart(rho)) || std::isnan(ImagPart(rho)) )
                RuntimeError("Givens rotation produced a NaN");
            H.Set( j, j, rho );
            cs.Set( j, 0, c );
            sn.Set( j, 0, s );
            // Apply the rotation to the rotated beta*e_0 vector
            // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            const F sConj = Conj(s); 
            const F tau_j = t.Get(j,0);
            const F tau_jp1 = t.Get(j+1,0); 
            t.Set( j,   0,  c    *tau_j + s*tau_jp1 );
            t.Set( j+1, 0, -sConj*tau_j + c*tau_jp1 );
            // Minimize the residual
            // ^^^^^^^^^^^^^^^^^^^^^
            auto tT = t( IR(0,j+1), ALL );
            auto HTL = H( IR(0,j+1), IR(0,j+1) );
            auto y = tT;
            Trsv( UPPER, NORMAL, NON_UNIT, HTL, y );
            // x := x0 + Vj y
            // ^^^^^^^^^^^^^^
            x = x0;
            for( Int i=0; i<=j; ++i )
            {
                const F eta_i = y.Get(i,0);
                Axpy( eta_i, V( ALL, IR(i) ), x );
            }

            // w := b - A x
            // ------------
            w = b;
            Multiply( NORMAL, F(-1), A, x, F(1), w );

            // Residual checks
            // ---------------
            const Real residNorm = Nrm2( w );
            if( std::isnan(residNorm) )
                RuntimeError("Residual norm was NaN");
            const Real relResidNorm = residNorm/origResidNorm; 
            if( relResidNorm < relTol )
            {
                if( progress )
                    cout << "  converged with relative tolerance: "
                         << relResidNorm << endl;
                converged = true;
                ++iter;
                break;
            }
            else
            {
                if( progress )
                    cout << "  finished iteration " << iter << " with "
                         << "relResidNorm=" << relResidNorm << endl;
            }
            ++iter;
        }
    }
    b = x;
    return maxLargeRefines;
}

template<typename F>
Int LGMRESSolveAfter
( const SparseMatrix<F>& A,   const Matrix<Base<F>>& reg,
  const Matrix<Base<F>>& d,
  const vector<Int>& invMap,  const SymmNodeInfo& info,
  const SymmFront<F>& front,        Matrix<F>& b,
  Base<F> relTol,                   Int k,
  Base<F> relTolRefine,             Int maxRefineIts,
  bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("reg_qsd_ldl::LGMRESSolveAfter"))
    typedef Base<F> Real;
    const Int n = A.Height();

    // x := 0
    // ======
    Matrix<F> x;
    Zeros( x, n, 1 );

    // w := b (= b - A x_0)
    // ====================
    Matrix<F> w; 
    w = b;
    const Real origResidNorm = Nrm2( w );
    if( origResidNorm == Real(0) )
        return 0;

    Int iter=0;
    Int maxLargeRefines=0;
    bool converged = false;
    Matrix<Real> cs;
    Matrix<F> sn, H, t;
    Matrix<F> x0, q, V;
    while( !converged )
    {
        if( progress )
            cout << "  Starting GMRES iteration " << iter << endl;
        Zeros( cs, k, 1 );
        Zeros( sn, k, 1 );
        Zeros( H,  k, k );
        Zeros( V, n, k );
        
        // x0 := x
        // =======
        x0 = x;

        // w := inv(M) w
        // =============
        Int refineIts = RegularizedSolveAfter
        ( A, reg, d, invMap, info, front, w, 
          relTolRefine, maxRefineIts, progress );
        maxLargeRefines = Max( refineIts, maxLargeRefines );

        // beta := || w ||_2
        // =================
        const Real beta = Nrm2( w );

        // v0 := w / beta
        // ==============
        auto v0 = V( ALL, IR(0) );
        v0 = w;
        Scale( Real(1)/beta, v0 ); 

        // t := beta e_0
        // =============
        Zeros( t, k+1, 1 );
        t.Set( 0, 0, beta );

        // Run one round of GMRES(k)
        // =========================
        for( Int j=0; j<k; ++j )
        {
            // w := A v_j
            // ----------
            Multiply( NORMAL, F(1), A, V(ALL,IR(j)), F(0), w );

            // w := inv(M) w
            // -------------
            Int refineIts = RegularizedSolveAfter
            ( A, reg, d, invMap, info, front, w, 
              relTolRefine, maxRefineIts, progress );
            maxLargeRefines = Max( refineIts, maxLargeRefines );

            // Run the j'th step of Arnoldi
            // ----------------------------
            for( Int i=0; i<=j; ++i )
            {
                // H(i,j) := v_i' w
                // ^^^^^^^^^^^^^^^^
                auto vi = V( ALL, IR(i) );
                H.Set( i, j, Dot(vi,w) ); 
              
                // w := w - H(i,j) v_i
                // ^^^^^^^^^^^^^^^^^^^
                Axpy( -H.Get(i,j), vi, w );
            }
            const Real delta = Nrm2( w );
            if( std::isnan(delta) )
                RuntimeError("Arnoldi step produced a NaN");
            if( delta == Real(0) )
                RuntimeError("Lucky breakdowns not yet handled");
            if( j+1 != k )
            {
                // v_{j+1} := w / delta
                // ^^^^^^^^^^^^^^^^^^^^^^^^^^
                auto vjp1 = V( ALL, IR(j+1) );
                vjp1 = w;
                Scale( Real(1)/delta, vjp1 );
            }

            // Apply existing rotations to the new column of H
            // -----------------------------------------------
            for( Int i=0; i<j; ++i )
            {
                const Real c = cs.Get(i,0);
                const F s = sn.Get(i,0);
                const F sConj = Conj(s);
                const F eta_i_j = H.Get(i,j);
                const F eta_ip1_j = H.Get(i+1,j);
                H.Set( i,   j,  c    *eta_i_j + s*eta_ip1_j );
                H.Set( i+1, j, -sConj*eta_i_j + c*eta_ip1_j );
            }

            // Generate and apply a new rotation to both H and the rotated
            // beta*e_0 vector, t, then solve the minimum residual problem
            // -----------------------------------------------------------
            const F eta_j_j = H.Get(j,j);
            const F eta_jp1_j = delta;
            if( std::isnan(RealPart(eta_j_j))   || 
                std::isnan(ImagPart(eta_j_j))   ||
                std::isnan(RealPart(eta_jp1_j)) || 
                std::isnan(ImagPart(eta_jp1_j)) )
                RuntimeError("Either H(j,j) or H(j+1,j) was NaN");
            Real c;
            F s;
            F rho = lapack::Givens( eta_j_j, eta_jp1_j, &c, &s );
            if( std::isnan(c) || 
                std::isnan(RealPart(s)) || std::isnan(ImagPart(s)) || 
                std::isnan(RealPart(rho)) || std::isnan(ImagPart(rho)) )
                RuntimeError("Givens rotation produced a NaN");
            H.Set( j, j, rho );
            cs.Set( j, 0, c );
            sn.Set( j, 0, s );
            // Apply the rotation to the rotated beta*e_0 vector
            // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            const F sConj = Conj(s); 
            const F tau_j = t.Get(j,0);
            const F tau_jp1 = t.Get(j+1,0); 
            t.Set( j,   0,  c    *tau_j + s*tau_jp1 );
            t.Set( j+1, 0, -sConj*tau_j + c*tau_jp1 );
            // Minimize the residual
            // ^^^^^^^^^^^^^^^^^^^^^
            auto tT = t( IR(0,j+1), ALL );
            auto HTL = H( IR(0,j+1), IR(0,j+1) );
            auto y = tT;
            Trsv( UPPER, NORMAL, NON_UNIT, HTL, y );
            // x := x0 + Vj y
            // ^^^^^^^^^^^^^^
            x = x0;
            for( Int i=0; i<=j; ++i )
            {
                const F eta_i = y.Get(i,0);
                Axpy( eta_i, V( ALL, IR(i) ), x );
            }

            // w := b - A x
            // ------------
            w = b;
            Multiply( NORMAL, F(-1), A, x, F(1), w );

            // Residual checks
            // ---------------
            const Real residNorm = Nrm2( w );
            if( std::isnan(residNorm) )
                RuntimeError("Residual norm was NaN");
            const Real relResidNorm = residNorm/origResidNorm; 
            if( relResidNorm < relTol )
            {
                if( progress )
                    cout << "  converged with relative tolerance: "
                         << relResidNorm << endl;
                converged = true;
                ++iter;
                break;
            }
            else
            {
                if( progress )
                    cout << "  finished iteration " << iter << " with "
                         << "relResidNorm=" << relResidNorm << endl;
            }
            ++iter;
        }
    }
    b = x;
    return maxLargeRefines;
}

template<typename F>
Int LGMRESSolveAfter
( const DistSparseMatrix<F>& A,  const DistMultiVec<Base<F>>& reg,
  const DistMap& invMap,         const DistSymmNodeInfo& info,
  const DistSymmFront<F>& front,       DistMultiVec<F>& b,
  Base<F> relTol,                      Int k,
  Base<F> relTolRefine,                Int maxRefineIts,
  bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("reg_qsd_ldl::LGMRESSolveAfter"))
    typedef Base<F> Real;
    const Int n = A.Height();
    mpi::Comm comm = A.Comm();
    const Int commRank = mpi::Rank(comm);

    // x := 0
    // ======
    DistMultiVec<F> x(comm);
    Zeros( x, n, 1 );

    // w := b (= b - A x_0)
    // ====================
    DistMultiVec<F> w(comm); 
    w = b;
    const Real origResidNorm = Nrm2( w );
    if( origResidNorm == Real(0) )
        return 0;

    Int iter=0;
    Int maxLargeRefines=0;
    bool converged = false;
    Matrix<Real> cs;
    Matrix<F> sn, H, t;
    DistMultiVec<F> x0(comm), q(comm), V(comm);
    while( !converged )
    {
        if( progress && commRank == 0 )
            cout << "  Starting GMRES iteration " << iter << endl;
        Zeros( cs, k, 1 );
        Zeros( sn, k, 1 );
        Zeros( H,  k, k );
        Zeros( V, n, k );
        // TODO: Extend DistMultiVec so that it can be directly manipulated
        //       rather than requiring access to the local Matrix and staging
        //       through the temporary vector q
        auto& VLoc = V.Matrix();
        Zeros( q, n, 1 );
        
        // x0 := x
        // =======
        x0 = x;

        // w := inv(M) w
        // =============
        Int refineIts = RegularizedSolveAfter
        ( A, reg, invMap, info, front, w, 
          relTolRefine, maxRefineIts, progress );
        maxLargeRefines = Max( refineIts, maxLargeRefines );

        // beta := || w ||_2
        // =================
        const Real beta = Nrm2( w );

        // v0 := w / beta
        // ==============
        auto v0Loc = VLoc( ALL, IR(0) );
        v0Loc = w.Matrix();
        Scale( Real(1)/beta, v0Loc ); 

        // t := beta e_0
        // =============
        Zeros( t, k+1, 1 );
        t.Set( 0, 0, beta );

        // Run one round of GMRES(k)
        // =========================
        for( Int j=0; j<k; ++j )
        {
            // w := A v_j
            // ----------
            q.Matrix() = VLoc( ALL, IR(j) );
            Multiply( NORMAL, F(1), A, q, F(0), w );

            // w := inv(M) w
            // -------------
            Int refineIts = RegularizedSolveAfter
            ( A, reg, invMap, info, front, w, 
              relTolRefine, maxRefineIts, progress );
            maxLargeRefines = Max( refineIts, maxLargeRefines );

            // Run the j'th step of Arnoldi
            // ----------------------------
            for( Int i=0; i<=j; ++i )
            {
                // H(i,j) := v_i' w
                // ^^^^^^^^^^^^^^^^
                q.Matrix() = VLoc( ALL, IR(i) );
                H.Set( i, j, Dot(q,w) ); 
              
                // w := w - H(i,j) v_i
                // ^^^^^^^^^^^^^^^^^^^
                Axpy( -H.Get(i,j), q, w );
            }
            const Real delta = Nrm2( w );
            if( std::isnan(delta) )
                RuntimeError("Arnoldi step produced a NaN");
            if( delta == Real(0) )
                RuntimeError("Lucky breakdowns not yet handled");
            if( j+1 != k )
            {
                // v_{j+1} := w / delta
                // ^^^^^^^^^^^^^^^^^^^^^^^^^^
                auto v_jp1Loc = VLoc( ALL, IR(j+1) );
                v_jp1Loc = w.Matrix();
                Scale( Real(1)/delta, v_jp1Loc );
            }

            // Apply existing rotations to the new column of H
            // -----------------------------------------------
            for( Int i=0; i<j; ++i )
            {
                const Real c = cs.Get(i,0);
                const F s = sn.Get(i,0);
                const F sConj = Conj(s);
                const F eta_i_j = H.Get(i,j);
                const F eta_ip1_j = H.Get(i+1,j);
                H.Set( i,   j,  c    *eta_i_j + s*eta_ip1_j );
                H.Set( i+1, j, -sConj*eta_i_j + c*eta_ip1_j );
            }

            // Generate and apply a new rotation to both H and the rotated
            // beta*e_0 vector, t, then solve the minimum residual problem
            // -----------------------------------------------------------
            const F eta_j_j = H.Get(j,j);
            const F eta_jp1_j = delta;
            if( std::isnan(RealPart(eta_j_j))   || 
                std::isnan(ImagPart(eta_j_j))   ||
                std::isnan(RealPart(eta_jp1_j)) || 
                std::isnan(ImagPart(eta_jp1_j)) )
                RuntimeError("Either H(j,j) or H(j+1,j) was NaN");
            Real c;
            F s;
            F rho = lapack::Givens( eta_j_j, eta_jp1_j, &c, &s );
            if( std::isnan(c) || 
                std::isnan(RealPart(s)) || std::isnan(ImagPart(s)) || 
                std::isnan(RealPart(rho)) || std::isnan(ImagPart(rho)) )
                RuntimeError("Givens rotation produced a NaN");
            H.Set( j, j, rho );
            cs.Set( j, 0, c );
            sn.Set( j, 0, s );
            // Apply the rotation to the rotated beta*e_0 vector
            // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            const F sConj = Conj(s); 
            const F tau_j = t.Get(j,0);
            const F tau_jp1 = t.Get(j+1,0); 
            t.Set( j,   0,  c    *tau_j + s*tau_jp1 );
            t.Set( j+1, 0, -sConj*tau_j + c*tau_jp1 );
            // Minimize the residual
            // ^^^^^^^^^^^^^^^^^^^^^
            auto tT = t( IR(0,j+1), ALL );
            auto HTL = H( IR(0,j+1), IR(0,j+1) );
            auto y = tT;
            Trsv( UPPER, NORMAL, NON_UNIT, HTL, y );
            // x := x0 + Vj y
            // ^^^^^^^^^^^^^^
            x = x0;
            for( Int i=0; i<=j; ++i )
            {
                const F eta_i = y.Get(i,0);
                Axpy( eta_i, VLoc( ALL, IR(i) ), x.Matrix() );
            }

            // w := b - A x
            // ------------
            w = b;
            Multiply( NORMAL, F(-1), A, x, F(1), w );

            // Residual checks
            // ---------------
            const Real residNorm = Nrm2( w );
            if( std::isnan(residNorm) )
                RuntimeError("Residual norm was NaN");
            const Real relResidNorm = residNorm/origResidNorm; 
            if( relResidNorm < relTol )
            {
                if( progress && commRank == 0 )
                    cout << "  converged with relative tolerance: "
                         << relResidNorm << endl;
                converged = true;
                ++iter;
                break;
            }
            else
            {
                if( progress && commRank == 0 )
                    cout << "  finished iteration " << iter << " with "
                         << "relResidNorm=" << relResidNorm << endl;
            }
            ++iter;
        }
    }
    b = x;
    return maxLargeRefines;
}

template<typename F>
Int LGMRESSolveAfter
( const DistSparseMatrix<F>& A,    const DistMultiVec<Base<F>>& reg,
  const DistMultiVec<Base<F>>& d,
  const DistMap& invMap,           const DistSymmNodeInfo& info,
  const DistSymmFront<F>& front,         DistMultiVec<F>& b,
  Base<F> relTol,                        Int k,
  Base<F> relTolRefine,                  Int maxRefineIts,
  bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("reg_qsd_ldl::LGMRESSolveAfter"))
    typedef Base<F> Real;
    const Int n = A.Height();
    mpi::Comm comm = A.Comm();
    const Int commRank = mpi::Rank(comm);

    // x := 0
    // ======
    DistMultiVec<F> x(comm);
    Zeros( x, n, 1 );

    // w := b (= b - A x_0)
    // ====================
    DistMultiVec<F> w(comm); 
    w = b;
    const Real origResidNorm = Nrm2( w );
    if( origResidNorm == Real(0) )
        return 0;

    Int iter=0;
    Int maxLargeRefines=0;
    bool converged = false;
    Matrix<Real> cs;
    Matrix<F> sn, H, t;
    DistMultiVec<F> x0(comm), q(comm), V(comm);
    while( !converged )
    {
        if( progress && commRank == 0 )
            cout << "  Starting GMRES iteration " << iter << endl;
        Zeros( cs, k, 1 );
        Zeros( sn, k, 1 );
        Zeros( H,  k, k );
        Zeros( V, n, k );
        // TODO: Extend DistMultiVec so that it can be directly manipulated
        //       rather than requiring access to the local Matrix and staging
        //       through the temporary vector q
        auto& VLoc = V.Matrix();
        Zeros( q, n, 1 );
        
        // x0 := x
        // =======
        x0 = x;

        // w := inv(M) w
        // =============
        Int refineIts = RegularizedSolveAfter
        ( A, reg, d, invMap, info, front, w, 
          relTolRefine, maxRefineIts, progress );
        maxLargeRefines = Max( refineIts, maxLargeRefines );

        // beta := || w ||_2
        // =================
        const Real beta = Nrm2( w );

        // v0 := w / beta
        // ==============
        auto v0Loc = VLoc( ALL, IR(0) );
        v0Loc = w.Matrix();
        Scale( Real(1)/beta, v0Loc ); 

        // t := beta e_0
        // =============
        Zeros( t, k+1, 1 );
        t.Set( 0, 0, beta );

        // Run one round of GMRES(k)
        // =========================
        for( Int j=0; j<k; ++j )
        {
            // w := A v_j
            // ----------
            q.Matrix() = VLoc( ALL, IR(j) );
            Multiply( NORMAL, F(1), A, q, F(0), w );

            // w := inv(M) w
            // -------------
            Int refineIts = RegularizedSolveAfter
            ( A, reg, d, invMap, info, front, w, 
              relTolRefine, maxRefineIts, progress );
            maxLargeRefines = Max( refineIts, maxLargeRefines );

            // Run the j'th step of Arnoldi
            // ----------------------------
            for( Int i=0; i<=j; ++i )
            {
                // H(i,j) := v_i' w
                // ^^^^^^^^^^^^^^^^
                q.Matrix() = VLoc( ALL, IR(i) );
                H.Set( i, j, Dot(q,w) ); 
              
                // w := w - H(i,j) v_i
                // ^^^^^^^^^^^^^^^^^^^
                Axpy( -H.Get(i,j), q, w );
            }
            const Real delta = Nrm2( w );
            if( std::isnan(delta) )
                RuntimeError("Arnoldi step produced a NaN");
            if( delta == Real(0) )
                RuntimeError("Lucky breakdowns not yet handled");
            if( j+1 != k )
            {
                // v_{j+1} := w / delta
                // ^^^^^^^^^^^^^^^^^^^^^^^^^^
                auto v_jp1Loc = VLoc( ALL, IR(j+1) );
                v_jp1Loc = w.Matrix();
                Scale( Real(1)/delta, v_jp1Loc );
            }

            // Apply existing rotations to the new column of H
            // -----------------------------------------------
            for( Int i=0; i<j; ++i )
            {
                const Real c = cs.Get(i,0);
                const F s = sn.Get(i,0);
                const F sConj = Conj(s);
                const F eta_i_j = H.Get(i,j);
                const F eta_ip1_j = H.Get(i+1,j);
                H.Set( i,   j,  c    *eta_i_j + s*eta_ip1_j );
                H.Set( i+1, j, -sConj*eta_i_j + c*eta_ip1_j );
            }

            // Generate and apply a new rotation to both H and the rotated
            // beta*e_0 vector, t, then solve the minimum residual problem
            // -----------------------------------------------------------
            const F eta_j_j = H.Get(j,j);
            const F eta_jp1_j = delta;
            if( std::isnan(RealPart(eta_j_j))   || 
                std::isnan(ImagPart(eta_j_j))   ||
                std::isnan(RealPart(eta_jp1_j)) || 
                std::isnan(ImagPart(eta_jp1_j)) )
                RuntimeError("Either H(j,j) or H(j+1,j) was NaN");
            Real c;
            F s;
            F rho = lapack::Givens( eta_j_j, eta_jp1_j, &c, &s );
            if( std::isnan(c) || 
                std::isnan(RealPart(s)) || std::isnan(ImagPart(s)) || 
                std::isnan(RealPart(rho)) || std::isnan(ImagPart(rho)) )
                RuntimeError("Givens rotation produced a NaN");
            H.Set( j, j, rho );
            cs.Set( j, 0, c );
            sn.Set( j, 0, s );
            // Apply the rotation to the rotated beta*e_0 vector
            // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            const F sConj = Conj(s); 
            const F tau_j = t.Get(j,0);
            const F tau_jp1 = t.Get(j+1,0); 
            t.Set( j,   0,  c    *tau_j + s*tau_jp1 );
            t.Set( j+1, 0, -sConj*tau_j + c*tau_jp1 );
            // Minimize the residual
            // ^^^^^^^^^^^^^^^^^^^^^
            auto tT = t( IR(0,j+1), ALL );
            auto HTL = H( IR(0,j+1), IR(0,j+1) );
            auto y = tT;
            Trsv( UPPER, NORMAL, NON_UNIT, HTL, y );
            // x := x0 + Vj y
            // ^^^^^^^^^^^^^^
            x = x0;
            for( Int i=0; i<=j; ++i )
            {
                const F eta_i = y.Get(i,0);
                Axpy( eta_i, VLoc( ALL, IR(i) ), x.Matrix() );
            }

            // w := b - A x
            // ------------
            w = b;
            Multiply( NORMAL, F(-1), A, x, F(1), w );

            // Residual checks
            // ---------------
            const Real residNorm = Nrm2( w );
            if( std::isnan(residNorm) )
                RuntimeError("Residual norm was NaN");
            const Real relResidNorm = residNorm/origResidNorm; 
            if( relResidNorm < relTol )
            {
                if( progress && commRank == 0 )
                    cout << "  converged with relative tolerance: "
                         << relResidNorm << endl;
                converged = true;
                ++iter;
                break;
            }
            else
            {
                if( progress && commRank == 0 )
                    cout << "  finished iteration " << iter << " with "
                         << "relResidNorm=" << relResidNorm << endl;
            }
            ++iter;
        }
    }
    b = x;
    return maxLargeRefines;
}

// The pseudocode for Flexible GMRES can be found in "Algorithm 2.2" in
//   Youcef Saad
//   "A flexible inner-outer preconditioned GMRES algorithm"
//   SIAM J. Sci. Comput., Vol. 14, No. 2, pp. 461--469, 1993.

template<typename F>
Int FGMRESSolveAfter
( const SparseMatrix<F>& A,   const Matrix<Base<F>>& reg,
  const vector<Int>& invMap,  const SymmNodeInfo& info,
  const SymmFront<F>& front,        Matrix<F>& b,
  Base<F> relTol,                   Int k,
  Base<F> relTolRefine,             Int maxRefineIts,
  bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("reg_qsd_ldl::FGMRESSolveAfter"))
    typedef Base<F> Real;
    const Int n = A.Height();

    // x := 0
    // ======
    Matrix<F> x;
    Zeros( x, n, 1 );

    // w := b (= b - A x_0)
    // ====================
    auto w = b;
    const Real origResidNorm = Nrm2( w );
    if( progress )
        cout << "origResidNorm: " << origResidNorm << endl;
    if( origResidNorm == Real(0) )
        return 0;

    Int iter=0;
    Int maxLargeRefines=0;
    bool converged = false;
    Matrix<Real> cs;
    Matrix<F> sn, H, t;
    Matrix<F> x0, V, Z;
    while( !converged )
    {
        if( progress )
            cout << "  Starting FGMRES iteration " << iter << endl;
        Zeros( cs, k, 1 );
        Zeros( sn, k, 1 );
        Zeros( H,  k, k );
        Zeros( V, n, k );
        Zeros( Z, n, k );
        
        // x0 := x
        // =======
        x0 = x;

        // NOTE: w = b - A x already

        // beta := || w ||_2
        // =================
        const Real beta = Nrm2( w );

        // v0 := w / beta
        // ==============
        auto v0 = V( ALL, IR(0) );
        v0 = w;
        Scale( Real(1)/beta, v0 ); 

        // t := beta e_0
        // =============
        Zeros( t, k+1, 1 );
        t.Set( 0, 0, beta );

        // Run one round of GMRES(k)
        // =========================
        for( Int j=0; j<k; ++j )
        {
            // z_j := inv(M) v_j
            // =================
            auto vj = V( ALL, IR(j) );
            auto zj = Z( ALL, IR(j) );
            zj = vj;
            Int refineIts = RegularizedSolveAfter
            ( A, reg, invMap, info, front, zj, 
              relTolRefine, maxRefineIts, progress );
            maxLargeRefines = Max( refineIts, maxLargeRefines );

            // w := A z_j
            // ----------
            Multiply( NORMAL, F(1), A, zj, F(0), w );

            // Run the j'th step of Arnoldi
            // ----------------------------
            for( Int i=0; i<=j; ++i )
            {
                // H(i,j) := v_i' w
                // ^^^^^^^^^^^^^^^^
                auto vi = V( ALL, IR(i) );
                H.Set( i, j, Dot(vi,w) ); 
              
                // w := w - H(i,j) v_i
                // ^^^^^^^^^^^^^^^^^^^
                Axpy( -H.Get(i,j), vi, w );
            }
            const Real delta = Nrm2( w );
            if( std::isnan(delta) )
                RuntimeError("Arnoldi step produced a NaN");
            if( delta == Real(0) )
                RuntimeError("Lucky breakdowns not yet handled");
            if( j+1 != k )
            {
                // v_{j+1} := w / delta
                // ^^^^^^^^^^^^^^^^^^^^^^^^^^
                auto vjp1 = V( ALL, IR(j+1) );
                vjp1 = w;
                Scale( Real(1)/delta, vjp1 );
            }

            // Apply existing rotations to the new column of H
            // -----------------------------------------------
            for( Int i=0; i<j; ++i )
            {
                const Real c = cs.Get(i,0);
                const F s = sn.Get(i,0);
                const F sConj = Conj(s);
                const F eta_i_j = H.Get(i,j);
                const F eta_ip1_j = H.Get(i+1,j);
                H.Set( i,   j,  c    *eta_i_j + s*eta_ip1_j );
                H.Set( i+1, j, -sConj*eta_i_j + c*eta_ip1_j );
            }

            // Generate and apply a new rotation to both H and the rotated
            // beta*e_0 vector, t, then solve the minimum residual problem
            // -----------------------------------------------------------
            const F eta_j_j = H.Get(j,j);
            const F eta_jp1_j = delta;
            if( std::isnan(RealPart(eta_j_j))   || 
                std::isnan(ImagPart(eta_j_j))   ||
                std::isnan(RealPart(eta_jp1_j)) || 
                std::isnan(ImagPart(eta_jp1_j)) )
                RuntimeError("Either H(j,j) or H(j+1,j) was NaN");
            Real c;
            F s;
            F rho = lapack::Givens( eta_j_j, eta_jp1_j, &c, &s );
            if( std::isnan(c) || 
                std::isnan(RealPart(s)) || std::isnan(ImagPart(s)) || 
                std::isnan(RealPart(rho)) || std::isnan(ImagPart(rho)) )
                RuntimeError("Givens rotation produced a NaN");
            H.Set( j, j, rho );
            cs.Set( j, 0, c );
            sn.Set( j, 0, s );
            // Apply the rotation to the rotated beta*e_0 vector
            // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            const F sConj = Conj(s); 
            const F tau_j = t.Get(j,0);
            const F tau_jp1 = t.Get(j+1,0); 
            t.Set( j,   0,  c    *tau_j + s*tau_jp1 );
            t.Set( j+1, 0, -sConj*tau_j + c*tau_jp1 );
            // Minimize the residual
            // ^^^^^^^^^^^^^^^^^^^^^
            auto tT = t( IR(0,j+1), ALL );
            auto HTL = H( IR(0,j+1), IR(0,j+1) );
            auto y = tT;
            Trsv( UPPER, NORMAL, NON_UNIT, HTL, y );
            // x := x0 + Zj y
            // ^^^^^^^^^^^^^^
            x = x0;
            for( Int i=0; i<=j; ++i )
            {
                const F eta_i = y.Get(i,0);
                Axpy( eta_i, Z( ALL, IR(i) ), x );
            }

            // w := b - A x
            // ------------
            w = b;
            Multiply( NORMAL, F(-1), A, x, F(1), w );

            // Residual checks
            // ---------------
            const Real residNorm = Nrm2( w );
            if( std::isnan(residNorm) )
                RuntimeError("Residual norm was NaN");
            const Real relResidNorm = residNorm/origResidNorm; 
            if( relResidNorm < relTol )
            {
                if( progress )
                    cout << "  converged with relative tolerance: "
                         << relResidNorm << endl;
                converged = true;
                ++iter;
                break;
            }
            else
            {
                if( progress )
                    cout << "  finished iteration " << iter << " with "
                         << "relResidNorm=" << relResidNorm << endl;
            }
            ++iter;
        }
    }
    b = x;
    return maxLargeRefines;
}

template<typename F>
Int FGMRESSolveAfter
( const SparseMatrix<F>& A,   const Matrix<Base<F>>& reg,
  const Matrix<Base<F>>& d,
  const vector<Int>& invMap,  const SymmNodeInfo& info,
  const SymmFront<F>& front,        Matrix<F>& b,
  Base<F> relTol,                   Int k,
  Base<F> relTolRefine,             Int maxRefineIts,
  bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("reg_qsd_ldl::FGMRESSolveAfter"))
    typedef Base<F> Real;
    const Int n = A.Height();

    // x := 0
    // ======
    Matrix<F> x;
    Zeros( x, n, 1 );

    // w := b (= b - A x_0)
    // ====================
    auto w = b;
    const Real origResidNorm = Nrm2( w );
    if( progress )
        cout << "origResidNorm: " << origResidNorm << endl;
    if( origResidNorm == Real(0) )
        return 0;

    Int iter=0;
    Int maxLargeRefines=0;
    bool converged = false;
    Matrix<Real> cs;
    Matrix<F> sn, H, t;
    Matrix<F> x0, V, Z;
    while( !converged )
    {
        if( progress )
            cout << "  Starting FGMRES iteration " << iter << endl;
        Zeros( cs, k, 1 );
        Zeros( sn, k, 1 );
        Zeros( H,  k, k );
        Zeros( V, n, k );
        Zeros( Z, n, k );
        
        // x0 := x
        // =======
        x0 = x;

        // NOTE: w = b - A x already

        // beta := || w ||_2
        // =================
        const Real beta = Nrm2( w );

        // v0 := w / beta
        // ==============
        auto v0 = V( ALL, IR(0) );
        v0 = w;
        Scale( Real(1)/beta, v0 ); 

        // t := beta e_0
        // =============
        Zeros( t, k+1, 1 );
        t.Set( 0, 0, beta );

        // Run one round of GMRES(k)
        // =========================
        for( Int j=0; j<k; ++j )
        {
            // z_j := inv(M) v_j
            // =================
            auto vj = V( ALL, IR(j) );
            auto zj = Z( ALL, IR(j) );
            zj = vj;
            Int refineIts = RegularizedSolveAfter
            ( A, reg, d, invMap, info, front, zj, 
              relTolRefine, maxRefineIts, progress );
            maxLargeRefines = Max( refineIts, maxLargeRefines );

            // w := A z_j
            // ----------
            Multiply( NORMAL, F(1), A, zj, F(0), w );

            // Run the j'th step of Arnoldi
            // ----------------------------
            for( Int i=0; i<=j; ++i )
            {
                // H(i,j) := v_i' w
                // ^^^^^^^^^^^^^^^^
                auto vi = V( ALL, IR(i) );
                H.Set( i, j, Dot(vi,w) ); 
              
                // w := w - H(i,j) v_i
                // ^^^^^^^^^^^^^^^^^^^
                Axpy( -H.Get(i,j), vi, w );
            }
            const Real delta = Nrm2( w );
            if( std::isnan(delta) )
                RuntimeError("Arnoldi step produced a NaN");
            if( delta == Real(0) )
                RuntimeError("Lucky breakdowns not yet handled");
            if( j+1 != k )
            {
                // v_{j+1} := w / delta
                // ^^^^^^^^^^^^^^^^^^^^^^^^^^
                auto vjp1 = V( ALL, IR(j+1) );
                vjp1 = w;
                Scale( Real(1)/delta, vjp1 );
            }

            // Apply existing rotations to the new column of H
            // -----------------------------------------------
            for( Int i=0; i<j; ++i )
            {
                const Real c = cs.Get(i,0);
                const F s = sn.Get(i,0);
                const F sConj = Conj(s);
                const F eta_i_j = H.Get(i,j);
                const F eta_ip1_j = H.Get(i+1,j);
                H.Set( i,   j,  c    *eta_i_j + s*eta_ip1_j );
                H.Set( i+1, j, -sConj*eta_i_j + c*eta_ip1_j );
            }

            // Generate and apply a new rotation to both H and the rotated
            // beta*e_0 vector, t, then solve the minimum residual problem
            // -----------------------------------------------------------
            const F eta_j_j = H.Get(j,j);
            const F eta_jp1_j = delta;
            if( std::isnan(RealPart(eta_j_j))   || 
                std::isnan(ImagPart(eta_j_j))   ||
                std::isnan(RealPart(eta_jp1_j)) || 
                std::isnan(ImagPart(eta_jp1_j)) )
                RuntimeError("Either H(j,j) or H(j+1,j) was NaN");
            Real c;
            F s;
            F rho = lapack::Givens( eta_j_j, eta_jp1_j, &c, &s );
            if( std::isnan(c) || 
                std::isnan(RealPart(s)) || std::isnan(ImagPart(s)) || 
                std::isnan(RealPart(rho)) || std::isnan(ImagPart(rho)) )
                RuntimeError("Givens rotation produced a NaN");
            H.Set( j, j, rho );
            cs.Set( j, 0, c );
            sn.Set( j, 0, s );
            // Apply the rotation to the rotated beta*e_0 vector
            // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            const F sConj = Conj(s); 
            const F tau_j = t.Get(j,0);
            const F tau_jp1 = t.Get(j+1,0); 
            t.Set( j,   0,  c    *tau_j + s*tau_jp1 );
            t.Set( j+1, 0, -sConj*tau_j + c*tau_jp1 );
            // Minimize the residual
            // ^^^^^^^^^^^^^^^^^^^^^
            auto tT = t( IR(0,j+1), ALL );
            auto HTL = H( IR(0,j+1), IR(0,j+1) );
            auto y = tT;
            Trsv( UPPER, NORMAL, NON_UNIT, HTL, y );
            // x := x0 + Zj y
            // ^^^^^^^^^^^^^^
            x = x0;
            for( Int i=0; i<=j; ++i )
            {
                const F eta_i = y.Get(i,0);
                Axpy( eta_i, Z( ALL, IR(i) ), x );
            }

            // w := b - A x
            // ------------
            w = b;
            Multiply( NORMAL, F(-1), A, x, F(1), w );

            // Residual checks
            // ---------------
            const Real residNorm = Nrm2( w );
            if( std::isnan(residNorm) )
                RuntimeError("Residual norm was NaN");
            const Real relResidNorm = residNorm/origResidNorm; 
            if( relResidNorm < relTol )
            {
                if( progress )
                    cout << "  converged with relative tolerance: "
                         << relResidNorm << endl;
                converged = true;
                ++iter;
                break;
            }
            else
            {
                if( progress )
                    cout << "  finished iteration " << iter << " with "
                         << "relResidNorm=" << relResidNorm << endl;
            }
            ++iter;
        }
    }
    b = x;
    return maxLargeRefines;
}

template<typename F>
Int FGMRESSolveAfter
( const DistSparseMatrix<F>& A,  const DistMultiVec<Base<F>>& reg,
  const DistMap& invMap,         const DistSymmNodeInfo& info,
  const DistSymmFront<F>& front,       DistMultiVec<F>& b,
  Base<F> relTol,                      Int k,
  Base<F> relTolRefine,                Int maxRefineIts,
  bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("reg_qsd_ldl::FGMRESSolveAfter"))
    typedef Base<F> Real;
    const Int n = A.Height();
    mpi::Comm comm = A.Comm();
    const Int commRank = mpi::Rank(comm);

    // x := 0
    // ======
    DistMultiVec<F> x(comm);
    Zeros( x, n, 1 );

    // w := b (= b - A x_0)
    // ====================
    DistMultiVec<F> w(comm); 
    w = b;
    const Real origResidNorm = Nrm2( w );
    if( progress && commRank == 0 )
        cout << "origResidNorm: " << origResidNorm << endl;
    if( origResidNorm == Real(0) )
        return 0;

    Int iter=0;
    Int maxLargeRefines=0;
    bool converged = false;
    Matrix<Real> cs;
    Matrix<F> sn, H, t;
    DistMultiVec<F> x0(comm), q(comm), V(comm), Z(comm);
    while( !converged )
    {
        if( progress && commRank == 0 )
            cout << "  Starting FGMRES iteration " << iter << endl;
        Zeros( cs, k, 1 );
        Zeros( sn, k, 1 );
        Zeros( H,  k, k );
        Zeros( V, n, k );
        Zeros( Z, n, k );
        // TODO: Extend DistMultiVec so that it can be directly manipulated
        //       rather than requiring access to the local Matrix and staging
        //       through the temporary vector q
        auto& VLoc = V.Matrix();
        auto& ZLoc = Z.Matrix();
        Zeros( q, n, 1 );
        
        // x0 := x
        // =======
        x0 = x;

        // NOTE: w = b - A x already

        // beta := || w ||_2
        // =================
        const Real beta = Nrm2( w );

        // v0 := w / beta
        // ==============
        auto v0Loc = VLoc( ALL, IR(0) );
        v0Loc = w.Matrix();
        Scale( Real(1)/beta, v0Loc ); 

        // t := beta e_0
        // =============
        Zeros( t, k+1, 1 );
        t.Set( 0, 0, beta );

        // Run one round of GMRES(k)
        // =========================
        for( Int j=0; j<k; ++j )
        {
            // z_j := inv(M) v_j
            // =================
            auto vjLoc = VLoc( ALL, IR(j) );
            auto zjLoc = ZLoc( ALL, IR(j) );
            q.Matrix() = vjLoc;
            Int refineIts = RegularizedSolveAfter
            ( A, reg, invMap, info, front, q, 
              relTolRefine, maxRefineIts, progress );
            maxLargeRefines = Max( refineIts, maxLargeRefines );
            zjLoc = q.Matrix();

            // w := A z_j
            // ----------
            // NOTE: q currently contains z_j
            Multiply( NORMAL, F(1), A, q, F(0), w );

            // Run the j'th step of Arnoldi
            // ----------------------------
            for( Int i=0; i<=j; ++i )
            {
                // H(i,j) := v_i' w
                // ^^^^^^^^^^^^^^^^
                q.Matrix() = VLoc( ALL, IR(i) );
                H.Set( i, j, Dot(q,w) ); 
              
                // w := w - H(i,j) v_i
                // ^^^^^^^^^^^^^^^^^^^
                Axpy( -H.Get(i,j), q, w );
            }
            const Real delta = Nrm2( w );
            if( std::isnan(delta) )
                RuntimeError("Arnoldi step produced a NaN");
            if( delta == Real(0) )
                RuntimeError("Lucky breakdowns not yet handled");
            if( j+1 != k )
            {
                // v_{j+1} := w / delta
                // ^^^^^^^^^^^^^^^^^^^^^^^^^^
                auto v_jp1Loc = VLoc( ALL, IR(j+1) );
                v_jp1Loc = w.Matrix();
                Scale( Real(1)/delta, v_jp1Loc );
            }

            // Apply existing rotations to the new column of H
            // -----------------------------------------------
            for( Int i=0; i<j; ++i )
            {
                const Real c = cs.Get(i,0);
                const F s = sn.Get(i,0);
                const F sConj = Conj(s);
                const F eta_i_j = H.Get(i,j);
                const F eta_ip1_j = H.Get(i+1,j);
                H.Set( i,   j,  c    *eta_i_j + s*eta_ip1_j );
                H.Set( i+1, j, -sConj*eta_i_j + c*eta_ip1_j );
            }

            // Generate and apply a new rotation to both H and the rotated
            // beta*e_0 vector, t, then solve the minimum residual problem
            // -----------------------------------------------------------
            const F eta_j_j = H.Get(j,j);
            const F eta_jp1_j = delta;
            if( std::isnan(RealPart(eta_j_j))   || 
                std::isnan(ImagPart(eta_j_j))   ||
                std::isnan(RealPart(eta_jp1_j)) || 
                std::isnan(ImagPart(eta_jp1_j)) )
                RuntimeError("Either H(j,j) or H(j+1,j) was NaN");
            Real c;
            F s;
            F rho = lapack::Givens( eta_j_j, eta_jp1_j, &c, &s );
            if( std::isnan(c) || 
                std::isnan(RealPart(s)) || std::isnan(ImagPart(s)) || 
                std::isnan(RealPart(rho)) || std::isnan(ImagPart(rho)) )
                RuntimeError("Givens rotation produced a NaN");
            H.Set( j, j, rho );
            cs.Set( j, 0, c );
            sn.Set( j, 0, s );
            // Apply the rotation to the rotated beta*e_0 vector
            // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            const F sConj = Conj(s); 
            const F tau_j = t.Get(j,0);
            const F tau_jp1 = t.Get(j+1,0); 
            t.Set( j,   0,  c    *tau_j + s*tau_jp1 );
            t.Set( j+1, 0, -sConj*tau_j + c*tau_jp1 );
            // Minimize the residual
            // ^^^^^^^^^^^^^^^^^^^^^
            auto tT = t( IR(0,j+1), ALL );
            auto HTL = H( IR(0,j+1), IR(0,j+1) );
            auto y = tT;
            Trsv( UPPER, NORMAL, NON_UNIT, HTL, y );
            // x := x0 + Zj y
            // ^^^^^^^^^^^^^^
            x = x0;
            for( Int i=0; i<=j; ++i )
            {
                const F eta_i = y.Get(i,0);
                Axpy( eta_i, ZLoc( ALL, IR(i) ), x.Matrix() );
            }

            // w := b - A x
            // ------------
            w = b;
            Multiply( NORMAL, F(-1), A, x, F(1), w );

            // Residual checks
            // ---------------
            const Real residNorm = Nrm2( w );
            if( std::isnan(residNorm) )
                RuntimeError("Residual norm was NaN");
            const Real relResidNorm = residNorm/origResidNorm; 
            if( relResidNorm < relTol )
            {
                if( progress && commRank == 0 )
                    cout << "  converged with relative tolerance: "
                         << relResidNorm << endl;
                converged = true;
                ++iter;
                break;
            }
            else
            {
                if( progress && commRank == 0 )
                    cout << "  finished iteration " << iter << " with "
                         << "relResidNorm=" << relResidNorm << endl;
            }
            ++iter;
        }
    }
    b = x;
    return maxLargeRefines;
}

template<typename F>
Int FGMRESSolveAfter
( const DistSparseMatrix<F>& A,    const DistMultiVec<Base<F>>& reg,
  const DistMultiVec<Base<F>>& d,
  const DistMap& invMap,           const DistSymmNodeInfo& info,
  const DistSymmFront<F>& front,         DistMultiVec<F>& b,
  Base<F> relTol,                        Int k,
  Base<F> relTolRefine,                  Int maxRefineIts,
  bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("reg_qsd_ldl::FGMRESSolveAfter"))
    typedef Base<F> Real;
    const Int n = A.Height();
    mpi::Comm comm = A.Comm();
    const Int commRank = mpi::Rank(comm);

    // x := 0
    // ======
    DistMultiVec<F> x(comm);
    Zeros( x, n, 1 );

    // w := b (= b - A x_0)
    // ====================
    DistMultiVec<F> w(comm); 
    w = b;
    const Real origResidNorm = Nrm2( w );
    if( progress && commRank == 0 )
        cout << "origResidNorm: " << origResidNorm << endl;
    if( origResidNorm == Real(0) )
        return 0;

    Int iter=0;
    Int maxLargeRefines=0;
    bool converged = false;
    Matrix<Real> cs;
    Matrix<F> sn, H, t;
    DistMultiVec<F> x0(comm), q(comm), V(comm), Z(comm);
    while( !converged )
    {
        if( progress && commRank == 0 )
            cout << "  Starting FGMRES iteration " << iter << endl;
        Zeros( cs, k, 1 );
        Zeros( sn, k, 1 );
        Zeros( H,  k, k );
        Zeros( V, n, k );
        Zeros( Z, n, k );
        // TODO: Extend DistMultiVec so that it can be directly manipulated
        //       rather than requiring access to the local Matrix and staging
        //       through the temporary vector q
        auto& VLoc = V.Matrix();
        auto& ZLoc = Z.Matrix();
        Zeros( q, n, 1 );
        
        // x0 := x
        // =======
        x0 = x;

        // NOTE: w = b - A x already

        // beta := || w ||_2
        // =================
        const Real beta = Nrm2( w );

        // v0 := w / beta
        // ==============
        auto v0Loc = VLoc( ALL, IR(0) );
        v0Loc = w.Matrix();
        Scale( Real(1)/beta, v0Loc ); 

        // t := beta e_0
        // =============
        Zeros( t, k+1, 1 );
        t.Set( 0, 0, beta );

        // Run one round of GMRES(k)
        // =========================
        for( Int j=0; j<k; ++j )
        {
            // z_j := inv(M) v_j
            // =================
            auto vjLoc = VLoc( ALL, IR(j) );
            auto zjLoc = ZLoc( ALL, IR(j) );
            q.Matrix() = vjLoc;
            Int refineIts = RegularizedSolveAfter
            ( A, reg, d, invMap, info, front, q, 
              relTolRefine, maxRefineIts, progress );
            maxLargeRefines = Max( refineIts, maxLargeRefines );
            zjLoc = q.Matrix();

            // w := A z_j
            // ----------
            // NOTE: q currently contains z_j
            Multiply( NORMAL, F(1), A, q, F(0), w );

            // Run the j'th step of Arnoldi
            // ----------------------------
            for( Int i=0; i<=j; ++i )
            {
                // H(i,j) := v_i' w
                // ^^^^^^^^^^^^^^^^
                q.Matrix() = VLoc( ALL, IR(i) );
                H.Set( i, j, Dot(q,w) ); 
              
                // w := w - H(i,j) v_i
                // ^^^^^^^^^^^^^^^^^^^
                Axpy( -H.Get(i,j), q, w );
            }
            const Real delta = Nrm2( w );
            if( std::isnan(delta) )
                RuntimeError("Arnoldi step produced a NaN");
            if( delta == Real(0) )
                RuntimeError("Lucky breakdowns not yet handled");
            if( j+1 != k )
            {
                // v_{j+1} := w / delta
                // ^^^^^^^^^^^^^^^^^^^^^^^^^^
                auto v_jp1Loc = VLoc( ALL, IR(j+1) );
                v_jp1Loc = w.Matrix();
                Scale( Real(1)/delta, v_jp1Loc );
            }

            // Apply existing rotations to the new column of H
            // -----------------------------------------------
            for( Int i=0; i<j; ++i )
            {
                const Real c = cs.Get(i,0);
                const F s = sn.Get(i,0);
                const F sConj = Conj(s);
                const F eta_i_j = H.Get(i,j);
                const F eta_ip1_j = H.Get(i+1,j);
                H.Set( i,   j,  c    *eta_i_j + s*eta_ip1_j );
                H.Set( i+1, j, -sConj*eta_i_j + c*eta_ip1_j );
            }

            // Generate and apply a new rotation to both H and the rotated
            // beta*e_0 vector, t, then solve the minimum residual problem
            // -----------------------------------------------------------
            const F eta_j_j = H.Get(j,j);
            const F eta_jp1_j = delta;
            if( std::isnan(RealPart(eta_j_j))   || 
                std::isnan(ImagPart(eta_j_j))   ||
                std::isnan(RealPart(eta_jp1_j)) || 
                std::isnan(ImagPart(eta_jp1_j)) )
                RuntimeError("Either H(j,j) or H(j+1,j) was NaN");
            Real c;
            F s;
            F rho = lapack::Givens( eta_j_j, eta_jp1_j, &c, &s );
            if( std::isnan(c) || 
                std::isnan(RealPart(s)) || std::isnan(ImagPart(s)) || 
                std::isnan(RealPart(rho)) || std::isnan(ImagPart(rho)) )
                RuntimeError("Givens rotation produced a NaN");
            H.Set( j, j, rho );
            cs.Set( j, 0, c );
            sn.Set( j, 0, s );
            // Apply the rotation to the rotated beta*e_0 vector
            // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            const F sConj = Conj(s); 
            const F tau_j = t.Get(j,0);
            const F tau_jp1 = t.Get(j+1,0); 
            t.Set( j,   0,  c    *tau_j + s*tau_jp1 );
            t.Set( j+1, 0, -sConj*tau_j + c*tau_jp1 );
            // Minimize the residual
            // ^^^^^^^^^^^^^^^^^^^^^
            auto tT = t( IR(0,j+1), ALL );
            auto HTL = H( IR(0,j+1), IR(0,j+1) );
            auto y = tT;
            Trsv( UPPER, NORMAL, NON_UNIT, HTL, y );
            // x := x0 + Zj y
            // ^^^^^^^^^^^^^^
            x = x0;
            for( Int i=0; i<=j; ++i )
            {
                const F eta_i = y.Get(i,0);
                Axpy( eta_i, ZLoc( ALL, IR(i) ), x.Matrix() );
            }

            // w := b - A x
            // ------------
            w = b;
            Multiply( NORMAL, F(-1), A, x, F(1), w );

            // Residual checks
            // ---------------
            const Real residNorm = Nrm2( w );
            if( std::isnan(residNorm) )
                RuntimeError("Residual norm was NaN");
            const Real relResidNorm = residNorm/origResidNorm; 
            if( relResidNorm < relTol )
            {
                if( progress && commRank == 0 )
                    cout << "  converged with relative tolerance: "
                         << relResidNorm << endl;
                converged = true;
                ++iter;
                break;
            }
            else
            {
                if( progress && commRank == 0 )
                    cout << "  finished iteration " << iter << " with "
                         << "relResidNorm=" << relResidNorm << endl;
            }
            ++iter;
        }
    }
    b = x;
    return maxLargeRefines;
}

// TODO: Add RGMRES

template<typename F>
Int SolveAfter
( const SparseMatrix<F>& A,   const Matrix<Base<F>>& reg,
  const vector<Int>& invMap,  const SymmNodeInfo& info,
  const SymmFront<F>& front,        Matrix<F>& b,
  const RegQSDCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("reg_qsd_ldl::SolveAfter"))
    switch( ctrl.alg )
    {
    case REG_REFINE_FGMRES:
        return FGMRESSolveAfter
        ( A, reg, invMap, info, front, b, 
          ctrl.relTol, ctrl.restart, ctrl.relTolRefine, ctrl.maxRefineIts, 
          ctrl.progress );
    case REG_REFINE_LGMRES:
        return LGMRESSolveAfter
        ( A, reg, invMap, info, front, b, 
          ctrl.relTol, ctrl.restart, ctrl.relTolRefine, ctrl.maxRefineIts, 
          ctrl.progress );
    case REG_REFINE_IR:
        return IRSolveAfter
        ( A, reg, invMap, info, front, b, 
          ctrl.relTolRefine, ctrl.maxRefineIts, ctrl.progress );
    case REG_REFINE_IR_MOD:    
        return RegularizedSolveAfter
        ( A, reg, invMap, info, front, b, 
          ctrl.relTolRefine, ctrl.maxRefineIts, ctrl.progress );
    default:
        LogicError("Invalid refinement algorithm");
    }
}

template<typename F>
Int SolveAfter
( const SparseMatrix<F>& A,   const Matrix<Base<F>>& reg,
  const Matrix<Base<F>>& d,
  const vector<Int>& invMap,  const SymmNodeInfo& info,
  const SymmFront<F>& front,        Matrix<F>& b,
  const RegQSDCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("reg_qsd_ldl::SolveAfter"))
    switch( ctrl.alg )
    {
    case REG_REFINE_FGMRES:
        return FGMRESSolveAfter
        ( A, reg, d, invMap, info, front, b, 
          ctrl.relTol, ctrl.restart, ctrl.relTolRefine, ctrl.maxRefineIts, 
          ctrl.progress );
    case REG_REFINE_LGMRES:
        return LGMRESSolveAfter
        ( A, reg, d, invMap, info, front, b, 
          ctrl.relTol, ctrl.restart, ctrl.relTolRefine, ctrl.maxRefineIts, 
          ctrl.progress );
    case REG_REFINE_IR:
        return IRSolveAfter
        ( A, reg, d, invMap, info, front, b, 
          ctrl.relTolRefine, ctrl.maxRefineIts, ctrl.progress );
    case REG_REFINE_IR_MOD:    
        return RegularizedSolveAfter
        ( A, reg, d, invMap, info, front, b, 
          ctrl.relTolRefine, ctrl.maxRefineIts, ctrl.progress );
    default:
        LogicError("Invalid refinement algorithm");
    }
}

template<typename F>
Int SolveAfter
( const DistSparseMatrix<F>& A,   const DistMultiVec<Base<F>>& reg,
  const DistMap& invMap,          const DistSymmNodeInfo& info,
  const DistSymmFront<F>& front,        DistMultiVec<F>& b,
  const RegQSDCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("reg_qsd_ldl::SolveAfter"))
    switch( ctrl.alg )
    {
    case REG_REFINE_FGMRES:
        return FGMRESSolveAfter
        ( A, reg, invMap, info, front, b, 
          ctrl.relTol, ctrl.restart, ctrl.relTolRefine, ctrl.maxRefineIts, 
          ctrl.progress );
    case REG_REFINE_LGMRES:
        return LGMRESSolveAfter
        ( A, reg, invMap, info, front, b, 
          ctrl.relTol, ctrl.restart, ctrl.relTolRefine, ctrl.maxRefineIts, 
          ctrl.progress );
    case REG_REFINE_IR:
        return IRSolveAfter
        ( A, reg, invMap, info, front, b, 
          ctrl.relTolRefine, ctrl.maxRefineIts, ctrl.progress );
    case REG_REFINE_IR_MOD:    
        return RegularizedSolveAfter
        ( A, reg, invMap, info, front, b, 
          ctrl.relTolRefine, ctrl.maxRefineIts, ctrl.progress );
    default:
        LogicError("Invalid refinement algorithm");
    }
}

template<typename F>
Int SolveAfter
( const DistSparseMatrix<F>& A,    const DistMultiVec<Base<F>>& reg,
  const DistMultiVec<Base<F>>& d,
  const DistMap& invMap,           const DistSymmNodeInfo& info,
  const DistSymmFront<F>& front,         DistMultiVec<F>& b,
  const RegQSDCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("reg_qsd_ldl::SolveAfter"))
    switch( ctrl.alg )
    {
    case REG_REFINE_FGMRES:
        return FGMRESSolveAfter
        ( A, reg, d, invMap, info, front, b, 
          ctrl.relTol, ctrl.restart, ctrl.relTolRefine, ctrl.maxRefineIts, 
          ctrl.progress );
    case REG_REFINE_LGMRES:
        return LGMRESSolveAfter
        ( A, reg, d, invMap, info, front, b, 
          ctrl.relTol, ctrl.restart, ctrl.relTolRefine, ctrl.maxRefineIts, 
          ctrl.progress );
    case REG_REFINE_IR:
        return IRSolveAfter
        ( A, reg, d, invMap, info, front, b, 
          ctrl.relTolRefine, ctrl.maxRefineIts, ctrl.progress );
    case REG_REFINE_IR_MOD:    
        return RegularizedSolveAfter
        ( A, reg, d, invMap, info, front, b, 
          ctrl.relTolRefine, ctrl.maxRefineIts, ctrl.progress );
    default:
        LogicError("Invalid refinement algorithm");
    }
}

#define PROTO(F) \
  template Int RegularizedSolveAfter \
  ( const SparseMatrix<F>& A,   const Matrix<Base<F>>& reg, \
    const vector<Int>& invMap,  const SymmNodeInfo& info, \
    const SymmFront<F>& front,        Matrix<F>& b, \
    Base<F> relTol,             Int maxRefineIts, \
    bool progress ); \
  template Int RegularizedSolveAfter \
  ( const SparseMatrix<F>& A,   const Matrix<Base<F>>& reg, \
    const Matrix<Base<F>>& d, \
    const vector<Int>& invMap,  const SymmNodeInfo& info, \
    const SymmFront<F>& front,        Matrix<F>& b, \
    Base<F> relTol,             Int maxRefineIts, \
    bool progress ); \
  template Int RegularizedSolveAfter \
  ( const DistSparseMatrix<F>& A,  const DistMultiVec<Base<F>>& reg, \
    const DistMap& invMap,         const DistSymmNodeInfo& info, \
    const DistSymmFront<F>& front,       DistMultiVec<F>& b, \
    Base<F> relTol,                Int maxRefineIts, \
    bool progress ); \
  template Int RegularizedSolveAfter \
  ( const DistSparseMatrix<F>& A,    const DistMultiVec<Base<F>>& reg, \
    const DistMultiVec<Base<F>>& d, \
    const DistMap& invMap,           const DistSymmNodeInfo& info, \
    const DistSymmFront<F>& front,         DistMultiVec<F>& b, \
    Base<F> relTol,                  Int maxRefineIts, \
    bool progress ); \
  template Int SolveAfter \
  ( const SparseMatrix<F>& A,  const Matrix<Base<F>>& reg, \
    const vector<Int>& invMap, const SymmNodeInfo& info, \
    const SymmFront<F>& front,       Matrix<F>& b, \
    const RegQSDCtrl<Base<F>>& ctrl ); \
  template Int SolveAfter \
  ( const SparseMatrix<F>& A,  const Matrix<Base<F>>& reg, \
    const Matrix<Base<F>>& d, \
    const vector<Int>& invMap, const SymmNodeInfo& info, \
    const SymmFront<F>& front,       Matrix<F>& b, \
    const RegQSDCtrl<Base<F>>& ctrl ); \
  template Int SolveAfter \
  ( const DistSparseMatrix<F>& A,      const DistMultiVec<Base<F>>& reg, \
    const DistMap& invMap,             const DistSymmNodeInfo& info, \
    const DistSymmFront<F>& front,           DistMultiVec<F>& b, \
    const RegQSDCtrl<Base<F>>& ctrl ); \
  template Int SolveAfter \
  ( const DistSparseMatrix<F>& A,      const DistMultiVec<Base<F>>& reg, \
    const DistMultiVec<Base<F>>& d, \
    const DistMap& invMap,             const DistSymmNodeInfo& info, \
    const DistSymmFront<F>& front,           DistMultiVec<F>& b, \
    const RegQSDCtrl<Base<F>>& ctrl );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace reg_qsd_ldl
} // namespace El
