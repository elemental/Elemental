/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_IMPORTS_QT5_HPP
#define EL_IMPORTS_QT5_HPP

#ifdef EL_HAVE_QT5

namespace El {

void InitializeQt5( int& argc, char**& argv );
void FinalizeQt5();

// Return true if Qt5 was detected, but no GUI is allowed.
// This is useful if images are to be written to file using Qt5's wrappers.
bool GuiDisabled();

// When Elemental is finalized, if no window was opened, then it must call 
// app.exit() instead
void OpenedWindow();

// For keeping track of the extreme values visualized so far
double MinRealWindowVal();
double MaxRealWindowVal();
double MinImagWindowVal();
double MaxImagWindowVal();
void UpdateMinRealWindowVal( double minVal );
void UpdateMaxRealWindowVal( double maxVal );
void UpdateMinImagWindowVal( double minVal );
void UpdateMaxImagWindowVal( double maxVal );

} // namespace El

#endif // ifdef EL_HAVE_QT5

#endif // ifndef EL_IMPORTS_QT5_HPP
