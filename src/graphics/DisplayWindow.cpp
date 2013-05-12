/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental-lite.hpp"
#include "elemental/graphics.hpp"

#ifdef HAVE_QT5

#include <QBoxLayout>

namespace elem {

DisplayWindow::DisplayWindow( QWidget* parent )
: QWidget(parent)
{
#ifndef RELEASE
    CallStackEntry entry("DisplayWindow::DisplayWindow");
#endif
    matrix_ = 0;

    // For the real matrix
    QHBoxLayout* matrixLayout = new QHBoxLayout();
    display_ = new DisplayWidget<double>();
    scroll_ = new QScrollArea();
    scroll_->setWidget( display_ );
    matrixLayout->addWidget( scroll_ );

    // Add two buttons at the bottom
    QHBoxLayout* buttonLayout = new QHBoxLayout();
    QPushButton* localButton = new QPushButton("Local");
    QPushButton* globalButton = new QPushButton("Global");
    buttonLayout->addWidget( localButton );
    buttonLayout->addWidget( globalButton );

    // Stack the matrix on top of the buttons
    QVBoxLayout* mainLayout = new QVBoxLayout();
    mainLayout->addLayout( matrixLayout );
    mainLayout->addLayout( buttonLayout );
    setLayout( mainLayout );

    connect( localButton, SIGNAL(clicked()), this, SLOT(UseLocalScale()) );
    connect( globalButton, SIGNAL(clicked()), this, SLOT(UseGlobalScale()) );

    // Elemental needs to know if a window was opened for cleanup purposes
    OpenedWindow();
}

DisplayWindow::~DisplayWindow()
{ delete matrix_; }

void
DisplayWindow::Display( const Matrix<double>* matrix, QString title )
{
#ifndef RELEASE
    CallStackEntry entry("DisplayWindow::Display");
#endif
    if( matrix_ != 0 )
        delete matrix_;
    matrix_ = matrix;

    setWindowTitle( title );
    display_->DisplayReal( matrix );
}

void
DisplayWindow::Display
( const Matrix<double>* matrix, double minVal, double maxVal, QString title )
{
#ifndef RELEASE
    CallStackEntry entry("DisplayWindow::Display");
#endif
    if( matrix_ != 0 )
        delete matrix_;
    matrix_ = matrix;

    setWindowTitle( title );
    display_->DisplayReal( matrix, minVal, maxVal );
}

void
DisplayWindow::UseLocalScale()
{
#ifndef RELEASE
    CallStackEntry entry("DisplayWindow::UseLocalScale");
#endif
    display_->DisplayReal( matrix_ );
}

void
DisplayWindow::UseGlobalScale()
{
#ifndef RELEASE
    CallStackEntry entry("DisplayWindow::UseGlobalScale");
#endif
    const double minVal = MinRealWindowVal();
    const double maxVal = MaxRealWindowVal();
    display_->DisplayReal( matrix_, minVal, maxVal );
}

} // namespace elem

#endif // ifdef HAVE_QT5
