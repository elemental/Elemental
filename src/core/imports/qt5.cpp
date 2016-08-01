/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>

#ifdef EL_HAVE_QT5
 #include <QApplication>
#endif

namespace {
using namespace El;

#ifdef EL_HAVE_QT5
// The command-line arguments should be passed into Qt5 in a manner which
// ensures that they do not fall out of scope until the last Qt5 call.
// The best way to do so is to make a copy and pass in the copy.
int argcSave;
char** argvSave;

bool guiDisabled;
bool elemInitializedQt = false;
bool elemOpenedWindow = false;
QCoreApplication* coreApp;
bool haveMinRealWindowVal=false, haveMaxRealWindowVal=false,
     haveMinImagWindowVal=false, haveMaxImagWindowVal=false;
double minRealWindowVal, maxRealWindowVal,
       minImagWindowVal, maxImagWindowVal;
#endif
}

namespace El {

#ifdef EL_HAVE_QT5
bool GuiDisabled()
{ return ::guiDisabled; }

void OpenedWindow()
{ ::elemOpenedWindow = true; }

double MinRealWindowVal()
{
    if( ::haveMinRealWindowVal )
        return ::minRealWindowVal;
    else
        return 0;
}

double MaxRealWindowVal()
{
    if( ::haveMaxRealWindowVal )
        return ::maxRealWindowVal;
    else
        return 0;
}

double MinImagWindowVal()
{
    if( ::haveMinImagWindowVal )
        return ::minImagWindowVal;
    else
        return 0;
}

double MaxImagWindowVal()
{
    if( ::haveMaxImagWindowVal )
        return ::maxImagWindowVal;
    else
        return 0;
}

void UpdateMinRealWindowVal( double minVal )
{
    if( ::haveMinRealWindowVal )
        ::minRealWindowVal = Min( ::minRealWindowVal, minVal );
    else
        ::minRealWindowVal = minVal;
    ::haveMinRealWindowVal = true;
}

void UpdateMaxRealWindowVal( double maxVal )
{
    if( ::haveMaxRealWindowVal )
        ::maxRealWindowVal = Max( ::maxRealWindowVal, maxVal );
    else
        ::maxRealWindowVal = maxVal;
    ::haveMaxRealWindowVal = true;
}

void UpdateMinImagWindowVal( double minVal )
{
    if( ::haveMinImagWindowVal )
        ::minImagWindowVal = Min( ::minImagWindowVal, minVal );
    else
        ::minImagWindowVal = minVal;
    ::haveMinImagWindowVal = true;
}

void UpdateMaxImagWindowVal( double maxVal )
{
    if( ::haveMaxImagWindowVal )
        ::maxImagWindowVal = Max( ::maxImagWindowVal, maxVal );
    else
        ::maxImagWindowVal = maxVal;
    ::haveMaxImagWindowVal = true;
}

void InitializeQt5( int& argc, char**& argv )
{
    ::coreApp = QCoreApplication::instance();
    if( ::coreApp == 0 )
    {
        // Test for whether the GUI should be disabled
        ::guiDisabled = false;
        for( int i=1; i<argc; ++i )
            if( !qstrcmp(argv[i],"-no-gui") )
                ::guiDisabled = true;

        ::argcSave = argc;
        ::argvSave = new char*[argc+1];
        for( int i=0; i<argc; ++i )
        {
            ::argvSave[i] = new char[std::strlen(argv[i])+1];
            std::strcpy( ::argvSave[i], argv[i] );
        }
       ::argvSave[argc] = nullptr;

        if( ::guiDisabled )
            ::coreApp = new QCoreApplication( ::argcSave, ::argvSave );
        else
            ::coreApp = new QApplication( ::argcSave, ::argvSave );
        ::elemInitializedQt = true;
    }
}

void FinalizeQt5()
{
    if( ::elemInitializedQt )
    {
        if( ::elemOpenedWindow )
            ::coreApp->exec();
        else
            ::coreApp->exit();
        delete ::coreApp;

        // Delete the copies of argc and argv
        for( int i=0; i< ::argcSave; ++i )
            delete[] ::argvSave[i]; 
        delete[] ::argvSave;
    }
}
#endif // ifdef EL_HAVE_QT5

} // namespace El
