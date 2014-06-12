/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

namespace El {
namespace kmeans {

template<typename F>
void AssignClusters
( const Matrix<F>& D, Matrix<Int>& cluster, Matrix<Base<F>>& dist )
{
    DEBUG_ONLY(CallStackEntry cse("kmeans::AssignClusters"))
    typedef Base<F> Real;
    const Int numPoints = D.Height();
    const Int numClusters = D.Width();
    cluster.Resize( numPoints, 1 );
    dist.Resize( numPoints, 1 );

    const Base<F> maxNorm = MaxNorm( D );
    std::vector<ValueInt<Real>> choices(numPoints);
    for( Int i=0; i<numPoints; ++i )
    {
        ValueInt<Real> choice;
        choice.value = 2*maxNorm;
        choice.index = -1;

        for( Int j=0; j<numClusters; ++j )
        {
            const Real rho = RealPart(D.Get(i,j));
            if( rho < choice.value ) 
            {
                choice.value = rho;
                choice.index = j;
            }
        }
        cluster.Set( i, 0, choice.index );
        dist.Set( i, 0, choice.value );
    }
}

template<typename F>
void AssignClusters
( const DistMatrix<F>& D, 
        DistMatrix<Int,    MC,STAR>& cluster, 
        DistMatrix<Base<F>,MC,STAR>& dist )
{
    DEBUG_ONLY(CallStackEntry cse("kmeans::AssignClusters"))
    typedef Base<F> Real;
    const Int numPoints = D.Height();
    const Int numClusters = D.Width();

    cluster.AlignWith( D );
    dist.AlignWith( D );
    cluster.Resize( numPoints, 1 );
    dist.Resize( numPoints, 1 );

    const Int numLocalPoints = D.LocalHeight();
    const Int numLocalClusters = D.LocalWidth();
    const Base<F> maxNorm = MaxNorm( D );
    std::vector<ValueInt<Real>> localChoices(numLocalPoints);
    for( Int iLoc=0; iLoc<numLocalPoints; ++iLoc )
    {
        localChoices[iLoc].value = 2*maxNorm;
        localChoices[iLoc].index = -1;

        for( Int jLoc=0; jLoc<numLocalClusters; ++jLoc )
        {
            const Int j = D.GlobalCol(jLoc);
            const Real rho = RealPart(D.GetLocal(iLoc,jLoc));
            if( rho < localChoices[iLoc].value ) 
            {
                localChoices[iLoc].value = rho;
                localChoices[iLoc].index = j;
            }
        }
    }
    mpi::AllReduce
    ( localChoices.data(), numLocalPoints, mpi::MinLocOp<Real>(), 
      D.RowComm() );
    for( Int iLoc=0; iLoc<numLocalPoints; ++iLoc )
    {
        cluster.SetLocal( iLoc, 0, localChoices[iLoc].index );
        dist.SetLocal( iLoc, 0, localChoices[iLoc].value );
    }
}

#define PROTO(F) \
  template void AssignClusters \
  ( const Matrix<F>& D, Matrix<Int>& cluster, Matrix<Base<F>>& dist ); \
  template void AssignClusters \
  ( const DistMatrix<F>& D,  \
          DistMatrix<Int,    MC,STAR>& cluster, \
          DistMatrix<Base<F>,MC,STAR>& dist );

PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace kmeans
} // namespace El
