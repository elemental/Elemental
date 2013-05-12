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

ComplexDisplayWindow::ComplexDisplayWindow( QWidget* parent )
: QWidget(parent)
{
#ifndef RELEASE
    CallStackEntry entry("ComplexDisplayWindow::ComplexDisplayWindow");
#endif
    matrix_ = 0;

    QHBoxLayout* matrixLayout = new QHBoxLayout();

    // Display the real data
    realDisplay_ = new DisplayWidget<Complex<double> >();
    realScroll_ = new QScrollArea();
    realScroll_->setWidget( realDisplay_ );
    matrixLayout->addWidget( realScroll_ );

    // Display the imaginary data
    imagDisplay_ = new DisplayWidget<Complex<double> >();
    imagScroll_ = new QScrollArea();
    imagScroll_->setWidget( imagDisplay_ );
    matrixLayout->addWidget( imagScroll_ );

    // Add two buttons underneath the two matrices
    QHBoxLayout* buttonLayout = new QHBoxLayout();
    QPushButton* localButton = new QPushButton("Local");
    QPushButton* globalButton = new QPushButton("Global");
    buttonLayout->addWidget( localButton );
    buttonLayout->addWidget( globalButton );

    QVBoxLayout* mainLayout = new QVBoxLayout();
    mainLayout->addLayout( matrixLayout );
    mainLayout->addLayout( buttonLayout );
    setLayout( mainLayout );

    connect( localButton, SIGNAL(clicked()), this, SLOT(UseLocalScale()) );
    connect( globalButton, SIGNAL(clicked()), this, SLOT(UseGlobalScale()) );
    setAttribute( Qt::WA_DeleteOnClose );

    // Elemental needs to know if a window was opened for cleanup purposes
    OpenedWindow();
}

ComplexDisplayWindow::~ComplexDisplayWindow()
{ delete matrix_; } 

void 
ComplexDisplayWindow::Display
( const Matrix<Complex<double> >* matrix, QString title )
{
#ifndef RELEASE
    CallStackEntry entry("ComplexDisplayWindow::Display");
#endif
    if( matrix_ != 0 )
        delete matrix_; 
    matrix_ = matrix;

    setWindowTitle( title );
    realDisplay_->DisplayReal( matrix );
    imagDisplay_->DisplayImag( matrix );
}

void 
ComplexDisplayWindow::Display
( const Matrix<Complex<double> >* matrix, 
  double minRealVal, double maxRealVal, 
  double minImagVal, double maxImagVal,
  QString title )
{
#ifndef RELEASE
    CallStackEntry entry("ComplexDisplayWindow::Display");
#endif
    if( matrix_ != 0 )
        delete matrix_;
    matrix_ = matrix;

    setWindowTitle( title );
    realDisplay_->DisplayReal( matrix, minRealVal, maxRealVal );
    imagDisplay_->DisplayImag( matrix, minImagVal, maxImagVal );
}

void
ComplexDisplayWindow::UseLocalScale()
{
#ifndef RELEASE
    CallStackEntry entry("ComplexDisplayWindow::UseLocalScale");
#endif
    realDisplay_->DisplayReal( matrix_ );
    imagDisplay_->DisplayImag( matrix_ );
}

void 
ComplexDisplayWindow::UseGlobalScale()
{
#ifndef RELEASE
    CallStackEntry entry("ComplexDisplayWindow::UseGlobalScale");
#endif
    const double minRealVal = MinRealWindowVal();
    const double maxRealVal = MaxRealWindowVal();
    const double minImagVal = MinImagWindowVal();
    const double maxImagVal = MaxImagWindowVal();
    realDisplay_->DisplayReal( matrix_, minRealVal, maxRealVal );
    imagDisplay_->DisplayImag( matrix_, minImagVal, maxImagVal );
}

} // namespace elem

#endif // ifdef HAVE_QT5
