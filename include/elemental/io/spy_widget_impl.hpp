/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_IO_SPYWIDGET_IMPL_HPP
#define ELEM_IO_SPYWIDGET_IMPL_HPP
#ifdef HAVE_QT5

#include <QPainter>
#include <QPixmap>
#include <QStylePainter>

#include "elemental/io/ColorMap.hpp"

#include "elemental/io/display_widget_decl.hpp"

namespace elem {

inline
SpyWidget::SpyWidget( QWidget* parent )
: QWidget(parent)
{ }

inline
SpyWidget::~SpyWidget()
{ }

inline void 
SpyWidget::paintEvent( QPaintEvent* event )
{
#ifndef RELEASE
    CallStackEntry entry("SpyWidget::paintEvent");
#endif
    QStylePainter painter( this );
    painter.drawPixmap( 0, 0, pixmap_ );
}

inline void 
SpyWidget::Spy( const Matrix<int>* A )
{
#ifndef RELEASE
    CallStackEntry entry("SpyWidget::Spy");
#endif
    const int m = A->Height();
    const int n = A->Width();

    // TODO: Parameterize these instead
    const int mPix = std::max( 500, m );
    const int nPix = std::max( 500, n );
    const double mRatio = double(m) / double(mPix);
    const double nRatio = double(n) / double(nPix);
    pixmap_ = QPixmap( nPix, mPix );
    resize( nPix, mPix );

    // Paint the matrix
    QPainter painter( &pixmap_ );
    painter.initFrom( this );
    for( int jPix=0; jPix<nPix; ++jPix ) 
    {
        const int j = nRatio*jPix;
        for( int iPix=0; iPix<mPix; ++iPix )
        {
            const int i = mRatio*iPix;
            const int z = A->Get(i,j);
            if( z == 0 )
                painter.setPen( qRgba(0,0,0,255) );
            else 
                painter.setPen( qRgba(255,255,255,255) );
            painter.drawPoint( jPix, iPix );
        }
    }

    // Refresh the widget
    update();
}

} // namespace elem

#endif // ifdef HAVE_QT5
#endif // ifndef ELEM_IO_SPYWIDGET_IMPL_HPP
