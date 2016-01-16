/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
#include "El/io/SpyWindow.hpp"

#ifdef EL_HAVE_QT5
# include <QApplication>
#endif

namespace El {

template<typename T>
void Spy( const Matrix<T>& A, string title, Base<T> tol )
{
    DEBUG_ONLY(CSE cse("Spy"))
#ifdef EL_HAVE_QT5
    if( GuiDisabled() )
        LogicError("GUI was disabled");

    const Int m = A.Height();
    const Int n = A.Width();
    Matrix<Int>* ASpy = new Matrix<Int>( m, n );
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            ASpy->Set( i, j, ( Abs(A.Get(i,j))>tol ? 1 : 0 ) );

    QString qTitle = QString::fromStdString( title );
    SpyWindow* spyWindow = new SpyWindow;
    spyWindow->Spy( ASpy, qTitle );
    spyWindow->show();

    // Spend at most 200 milliseconds rendering
    ProcessEvents( 200 );
#else
    LogicError("Qt5 not available");
#endif // ifdef EL_HAVE_QT5
}

template<typename T>
void Spy( const AbstractDistMatrix<T>& A, string title, Base<T> tol )
{
    DEBUG_ONLY(CSE cse("Spy"))
#ifdef EL_HAVE_QT5
    if( GuiDisabled() )
        LogicError("GUI was disabled");
    if( A.ColStride() == 1 && A.RowStride() == 1 )
    {
        if( A.CrossRank() == A.Root() && A.RedundantRank() == 0 )
            Spy( A.LockedMatrix(), title, tol );
    }
    else
    {
        DistMatrix<T,CIRC,CIRC> A_CIRC_CIRC( A );
        if( A_CIRC_CIRC.CrossRank() == A_CIRC_CIRC.Root() )
            Spy( A_CIRC_CIRC.Matrix(), title, tol );
    }
#else
    LogicError("Qt5 not available");
#endif // ifdef EL_HAVE_QT5
}

#define PROTO(T) \
  template void Spy ( const Matrix<T>& A, string title, Base<T> tol ); \
  template void Spy \
  ( const AbstractDistMatrix<T>& A, string title, Base<T> tol );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
