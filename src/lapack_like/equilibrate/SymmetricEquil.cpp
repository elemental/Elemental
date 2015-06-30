/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename F>
void SymmetricEquil
( Matrix<F>& A, Matrix<Base<F>>& d, 
  bool geomEquil, bool diagEquil,
  bool scaleTwoNorm, Int basisSize, 
  bool progress )
{
    DEBUG_ONLY(CSE cse("SymmetricEquil"))
    // TODO: Ensure A is square
    const Int n = A.Height();
    Ones( d, n, 1 );
    if( geomEquil )
        SymmetricGeomEquil( A, d, progress );
    else if( diagEquil )
    {
        if( progress )
            cout << "    Diagonal equilibration not yet enabled for dense "
                    "matrices" << endl;
    }
       
    if( scaleTwoNorm )
    {
        if( progress )
            cout << "    Two-norm estimation with basisSize=" << basisSize
                 << " not yet enabled for dense matrices" << endl;
    }
}

template<typename F>
void SymmetricEquil
( AbstractDistMatrix<F>& APre, AbstractDistMatrix<Base<F>>& dPre, 
  bool geomEquil, bool diagEquil,
  bool scaleTwoNorm, Int basisSize, bool progress )
{
    DEBUG_ONLY(CSE cse("SymmetricEquil"))
    ProxyCtrl control;
    control.colConstrain = true;
    control.rowConstrain = true;
    control.colAlign = 0;
    control.rowAlign = 0;
    auto APtr = ReadWriteProxy<F,MC,MR>(&APre,control);
    auto& A = *APtr;
    const Int n = A.Height();

    Ones( dPre, n, 1 );
    if( geomEquil )
        SymmetricGeomEquil( A, dPre, progress );
    else if( diagEquil )
    {
        if( progress )
            cout << "    Diagonal equilibration not yet enabled for dense "
                    "matrices" << endl;
    }

    if( scaleTwoNorm )
    {
        if( progress && mpi::Rank(A.DistComm()) )
            cout << "    Two-norm estimation with basisSize=" << basisSize
                 << " not yet enabled for dense matrices" << endl;
    }
}

template<typename F>
void SymmetricEquil
( SparseMatrix<F>& A, Matrix<Base<F>>& d, 
  bool geomEquil, bool diagEquil,
  bool scaleTwoNorm, Int basisSize, bool progress )
{
    DEBUG_ONLY(CSE cse("SymmetricEquil"))
    typedef Base<F> Real;
    const Int n = A.Height();

    if( geomEquil )
    {
        SymmetricGeomEquil( A, d, progress );
    }
    else if( diagEquil )
    {
        auto maxSqrtLambda = []( F delta ) 
                             { return Sqrt(Max(Abs(delta),Real(1))); };
        function<Real(F)> maxSqrt( maxSqrtLambda );
        GetMappedDiagonal( A, d, maxSqrt );
        if( progress )
        {
            const Real maxNorm = MaxNorm( d ); 
            cout << "    || d ||_max = " << maxNorm << endl;
        }
        DiagonalSolve( LEFT, NORMAL, d, A );
        DiagonalSolve( RIGHT, NORMAL, d, A );
    }
    else
        Ones( d, n, 1 );

    if( scaleTwoNorm )
    {
        Real twoNormEst = HermitianTwoNormEstimate( A, basisSize );
        if( progress )
            cout << "    Estimated two-norm as " << twoNormEst << endl;
        Scale( Real(1)/twoNormEst, A );
        Scale( Sqrt(twoNormEst), d );
    } 
}

template<typename F>
void SymmetricEquil
( DistSparseMatrix<F>& A, DistMultiVec<Base<F>>& d, 
  bool geomEquil, bool diagEquil,
  bool scaleTwoNorm, Int basisSize, bool progress, bool time )
{
    DEBUG_ONLY(CSE cse("SymmetricEquil"))
    typedef Base<F> Real;
    mpi::Comm comm = A.Comm();
    const int commRank = mpi::Rank(comm);
    const Int n = A.Height();
    Timer timer;

    d.SetComm( comm );
    if( geomEquil )
    {
        if( commRank == 0 && time )
            timer.Start();
        SymmetricGeomEquil( A, d, progress );
        if( commRank == 0 && time )
            cout << "  SymmetricGeomEquil time: " << timer.Stop() << " secs" 
                 << endl;
    }
    else if( diagEquil )
    {

        auto maxSqrtLambda = []( F delta )
                             { return Sqrt(Max(Abs(delta),Real(1))); };
        function<Real(F)> maxSqrt( maxSqrtLambda );
        if( commRank == 0 && time )
            timer.Start();
        GetMappedDiagonal( A, d, maxSqrt );
        if( commRank == 0 && time )
            cout << "    Get mapped diag time: " << timer.Stop() << endl;
        if( commRank == 0 && time )
            timer.Start();
        DiagonalSolve( LEFT, NORMAL, d, A );
        if( commRank == 0 && time )
            cout << "    Left diag solve time: " << timer.Stop() << endl;
        if( commRank == 0 && time )
            timer.Start();
        DiagonalSolve( RIGHT, NORMAL, d, A );
        if( commRank == 0 && time )
            cout << "    Right diag solve time: " << timer.Stop() << endl;
        if( progress )
        {
            const Real maxNorm = MaxNorm( d );
            if( commRank == 0 ) 
                cout << "    Diagonally equilibrated with || d ||_max = " 
                     << maxNorm << endl;
        }
    }
    else
        Ones( d, n, 1 );

    if( scaleTwoNorm )
    {
        if( commRank == 0 && time )
            timer.Start();
        Real twoNormEst = HermitianTwoNormEstimate( A, basisSize );
        if( progress && commRank == 0 )
            cout << "    Estimated two-norm as " << twoNormEst << endl;
        Scale( Real(1)/twoNormEst, A );
        Scale( Sqrt(twoNormEst), d );
        if( commRank == 0 && time )
            cout << "    Two-norm estimation time: " << timer.Stop() 
                 << " secs" << endl;
    } 
}

#define PROTO(F) \
  template void SymmetricEquil \
  ( Matrix<F>& A, Matrix<Base<F>>& d, \
    bool geomEquil, bool diagEquil, \
    bool scaleTwoNorm, Int basisSize, bool progress ); \
  template void SymmetricEquil \
  ( AbstractDistMatrix<F>& A,  AbstractDistMatrix<Base<F>>& d, \
    bool geomEquil, bool diagEquil, \
    bool scaleTwoNorm, Int basisSize, bool progress ); \
  template void SymmetricEquil \
  ( SparseMatrix<F>& A, Matrix<Base<F>>& d, \
    bool geomEquil, bool diagEquil, \
    bool scaleTwoNorm, Int basisSize, bool progress ); \
  template void SymmetricEquil \
  ( DistSparseMatrix<F>& A, DistMultiVec<Base<F>>& d, \
    bool geomEquil, bool diagEquil, \
    bool scaleTwoNorm, Int basisSize, bool progress, bool time );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
