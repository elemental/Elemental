/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SCHUR_CONDENSE_HPP
#define EL_SCHUR_CONDENSE_HPP

namespace El {
namespace schur {

template<typename F>
void Condense
( Matrix<F>& A,
  Matrix<Complex<Base<F>>>& w,
  const SchurCtrl<Base<F>>& ctrl )
{
    EL_DEBUG_CSE
    Timer timer;

    Matrix<F> householderScalars;
    if( ctrl.time )
        timer.Start();
    Hessenberg( UPPER, A, householderScalars );
    if( ctrl.time )
        Output("  Hessenberg reduction: ",timer.Stop()," seconds");
    MakeTrapezoidal( UPPER, A, -1 );

    if( ctrl.time )
        timer.Start();
    HessenbergSchur( A, w, ctrl.hessSchurCtrl ); 
    if( ctrl.time )
        Output("  HessenbergSchur: ",timer.Stop()," seconds");

    if( IsComplex<F>::value )
        MakeTrapezoidal( UPPER, A );
    else
    {
        MakeTrapezoidal( UPPER, A, -1 );
        EL_DEBUG_ONLY(CheckRealSchur(A))
    }
}

template<typename F>
void Condense
( Matrix<F>& A,
  Matrix<Complex<Base<F>>>& w,
  Matrix<F>& Q, 
  const SchurCtrl<Base<F>>& ctrl )
{
    EL_DEBUG_CSE
    Timer timer;

    Matrix<F> householderScalars;
    if( ctrl.time )
        timer.Start();
    Hessenberg( UPPER, A, householderScalars );
    if( ctrl.time )
        Output("  Hessenberg reduction: ",timer.Stop()," seconds");
    hessenberg::FormQ( UPPER, A, householderScalars, Q );
    MakeTrapezoidal( UPPER, A, -1 );

    auto hessSchurCtrl( ctrl.hessSchurCtrl );
    hessSchurCtrl.accumulateSchurVecs = true;
    if( ctrl.time )
        timer.Start();
    HessenbergSchur( A, w, Q, hessSchurCtrl ); 
    if( ctrl.time )
        Output("  HessenbergSchur: ",timer.Stop()," seconds");

    if( IsComplex<F>::value )
        MakeTrapezoidal( UPPER, A );
    else
    {
        MakeTrapezoidal( UPPER, A, -1 );
        EL_DEBUG_ONLY(CheckRealSchur(A))
    }
}

template<typename F>
void Condense
( AbstractDistMatrix<F>& A,
  AbstractDistMatrix<Complex<Base<F>>>& w,
  const SchurCtrl<Base<F>>& ctrl )
{
    EL_DEBUG_CSE
    const Grid& grid = A.Grid();
    Timer timer;

    // Reduce the matrix to upper-Hessenberg form in an elemental form
    if( ctrl.time && grid.Rank() == 0 )
        timer.Start();
    hessenberg::ExplicitCondensed( UPPER, A );
    if( ctrl.time && grid.Rank() == 0 )
        Output("  Hessenberg reduction: ",timer.Stop()," seconds"); 

    // Run the black-box Hessenberg Schur decomposition
    if( ctrl.time && grid.Rank() == 0 )
        timer.Start();
    HessenbergSchur( A, w, ctrl.hessSchurCtrl );
    if( ctrl.time && grid.Rank() == 0 )
        Output("  HessenbergSchur: ",timer.Stop()," seconds");

    if( IsComplex<F>::value )
        MakeTrapezoidal( UPPER, A );
    else
    {
        MakeTrapezoidal( UPPER, A, -1 );
        // NOTE: This routine is not yet implemented
        //DEBUG_ONLY(CheckRealSchur(A))
    }
}

template<typename F>
void Condense
( AbstractDistMatrix<F>& A,
  AbstractDistMatrix<Complex<Base<F>>>& w,
  AbstractDistMatrix<F>& Q,
  const SchurCtrl<Base<F>>& ctrl )
{
    EL_DEBUG_CSE
    const Grid& grid = A.Grid();
    Timer timer;

    // Reduce A to upper-Hessenberg form
    DistMatrix<F,STAR,STAR> householderScalars( A.Grid() );
    if( ctrl.time && grid.Rank() == 0 )
        timer.Start();
    Hessenberg( UPPER, A, householderScalars );
    if( ctrl.time && grid.Rank() == 0 )
        Output("  Hessenberg reduction: ",timer.Stop()," seconds");

    // Explicitly accumulate the Householder transformations into Q
    if( ctrl.time && grid.Rank() == 0 )
        timer.Start();
    hessenberg::FormQ( UPPER, A, householderScalars, Q );
    if( ctrl.time && grid.Rank() == 0 )
        Output("  hessenberg::FormQ: ",timer.Stop()," seconds");
    MakeTrapezoidal( UPPER, A, -1 );
    
    // Call the black-box HessenbergSchur decomposition
    if( ctrl.time && grid.Rank() == 0 )
        timer.Start();
    auto hessSchurCtrl = ctrl.hessSchurCtrl;
    hessSchurCtrl.accumulateSchurVecs = true;
    HessenbergSchur( A, w, Q, hessSchurCtrl );
    if( ctrl.time && grid.Rank() == 0 )
        Output("  HessenbergSchur: ",timer.Stop()," seconds");

    if( IsComplex<F>::value )
        MakeTrapezoidal( UPPER, A );
    else
    {
        MakeTrapezoidal( UPPER, A, -1 );
        // NOTE: This routine is not yet implemented
        //DEBUG_ONLY(CheckRealSchur(A))
    }
}

} // namespace schur
} // namespace El

#endif // ifndef EL_SCHUR_CONDENSE_HPP
