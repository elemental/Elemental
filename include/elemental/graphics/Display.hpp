/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef GRAPHICS_DISPLAY_HPP
#define GRAPHICS_DISPLAY_HPP

#ifdef HAVE_QT5

#include "elemental/graphics/display_window_decl.hpp"
#include "elemental/graphics/complex_display_window_decl.hpp"

#include "elemental/graphics/display_widget_impl.hpp"

namespace elem {

template<typename T>
inline void
Display( const Matrix<T>& A, std::string title="" )
{
#ifndef RELEASE
    CallStackEntry entry("Display");
#endif
    // Convert A to double-precision since Qt's MOC does not support templates
    const int m = A.Height();
    const int n = A.Width();
    Matrix<double>* ADouble = new Matrix<double>( m, n );
    for( int j=0; j<n; ++j )
        for( int i=0; i<m; ++i )
            ADouble->Set( i, j, A.Get(i,j) );

    QString qTitle = QString::fromStdString( title );
    DisplayWindow* displayWindow = new DisplayWindow;
    displayWindow->Display( ADouble, qTitle );
    displayWindow->show();
    RegisterDisplayWindow( displayWindow );
}

template<typename T>
inline void
Display( const Matrix<Complex<T> >& A, std::string title="" )
{
#ifndef RELEASE
    CallStackEntry entry("Display");
#endif
    // Convert A to double-precision since Qt's MOC does not support templates
    const int m = A.Height();
    const int n = A.Width();
    Matrix<Complex<double> >* ADouble = new Matrix<Complex<double> >( m, n );
    for( int j=0; j<n; ++j )
        for( int i=0; i<m; ++i )
            ADouble->Set( i, j, A.Get(i,j) );

    QString qTitle = QString::fromStdString( title );
    ComplexDisplayWindow* displayWindow = new ComplexDisplayWindow;
    displayWindow->Display( ADouble, qTitle );
    displayWindow->show();
    RegisterComplexDisplayWindow( displayWindow );
}

template<typename T,Distribution U,Distribution V>
inline void
Display( const DistMatrix<T,U,V>& A, std::string title="" )
{
#ifndef RELEASE
    CallStackEntry entry("Display");
#endif
    // TODO: Avoid giving every process a full copy and think about avoiding
    //       the extra copy needed when the underlying datatype is not 'double'
    DistMatrix<T,STAR,STAR> A_STAR_STAR( A );
    if( A.Grid().Rank() == 0 )
        Display( A_STAR_STAR.Matrix(), title );
}

} // namespace elem

#endif // ifdef HAVE_QT5

#endif // ifndef GRAPHICS_DISPLAY_HPP
