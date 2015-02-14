/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {
namespace kmeans {

// NOTE: The result will be real, even if F is complex

template<typename F>
void Distance( const Matrix<F>& X, const Matrix<F>& C, Matrix<F>& D )
{
    DEBUG_ONLY(CallStackEntry cse("kmeans::Distance"))
    const Int numPoints = X.Width();
    const Int numClusters = C.Width();

    // Compute the cross term, 2 RealPart( X^H C ) 
    Gemm( ADJOINT, NORMAL, F(-2), X, C, D );
    MakeReal( D );
    
    Matrix<Base<F>> xNorms, cNorms;
    ColumnNorms( X, xNorms );
    ColumnNorms( C, cNorms );

    auto squareMap = []( Base<F> alpha ) { return alpha*alpha; };
    EntrywiseMap( xNorms, function<Base<F>(Base<F>)>(squareMap) );
    EntrywiseMap( cNorms, function<Base<F>(Base<F>)>(squareMap) );

    for( Int j=0; j<numClusters; ++j )
        for( Int i=0; i<numPoints; ++i )
            D.Update( i, j, xNorms.Get(i,0)+cNorms.Get(j,0) );
}

template<typename F>
void Distance
( const DistMatrix<F>& X, const DistMatrix<F>& C, DistMatrix<F>& D )
{
    DEBUG_ONLY(CallStackEntry cse("kmeans::Distance"))

    // Compute the cross term, 2 RealPart( X^H C ) 
    Gemm( ADJOINT, NORMAL, F(-2), X, C, D );
    MakeReal( D );
    
    DistMatrix<Base<F>,MR,STAR> xNorms_MR_STAR(X.Grid()), 
                                cNorms_MR_STAR(X.Grid());
    cNorms_MR_STAR.AlignWith( D ); 
    ColumnNorms( X, xNorms_MR_STAR );
    ColumnNorms( C, cNorms_MR_STAR );

    auto squareMap = []( Base<F> alpha ) { return alpha*alpha; };
    EntrywiseMap( xNorms_MR_STAR, function<Base<F>(Base<F>)>(squareMap) );
    EntrywiseMap( cNorms_MR_STAR, function<Base<F>(Base<F>)>(squareMap) );

    DistMatrix<Base<F>,MC,STAR> xNorms_MC_STAR(X.Grid());
    xNorms_MC_STAR.AlignWith( D );
    xNorms_MC_STAR = xNorms_MR_STAR;

    for( Int jLoc=0; jLoc<D.LocalWidth(); ++jLoc )
        for( Int iLoc=0; iLoc<D.LocalHeight(); ++iLoc )
            D.UpdateLocal( iLoc, jLoc, xNorms_MC_STAR.GetLocal(iLoc,0)+
                                       cNorms_MR_STAR.GetLocal(jLoc,0) );
}

#define PROTO(F) \
  template void Distance \
  ( const Matrix<F>& X, const Matrix<F>& C, Matrix<F>& D ); \
  template void Distance \
  ( const DistMatrix<F>& X, const DistMatrix<F>& C, DistMatrix<F>& D );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace kmeans
} // namespace El
