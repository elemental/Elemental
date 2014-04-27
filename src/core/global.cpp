/*
   Copyright (c) 2009-2014, Jack Poulson
                      2013, Jeff Hammond
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "elemental-lite.hpp"
#ifdef ELEM_HAVE_QT5
 #include <QApplication>
#endif

namespace {
using namespace elem;

Int numElemInits = 0;
bool elemInitializedMpi = false;

std::stack<Int> blocksizeStack;
Grid* defaultGrid = 0;
Args* args = 0;

// Default blocksizes for BlockDistMatrix
Int blockHeight=32, blockWidth=32;

// A common Mersenne twister configuration
std::mt19937 generator;

// Debugging
DEBUG_ONLY(std::stack<std::string> callStack)

// Tuning parameters for basic routines
Int localSymvFloatBlocksize = 64;
Int localSymvDoubleBlocksize = 64;
Int localSymvComplexFloatBlocksize = 64;
Int localSymvComplexDoubleBlocksize = 64;

Int localTrr2kFloatBlocksize = 64;
Int localTrr2kDoubleBlocksize = 64;
Int localTrr2kComplexFloatBlocksize = 64;
Int localTrr2kComplexDoubleBlocksize = 64;

Int localTrrkFloatBlocksize = 64;
Int localTrrkDoubleBlocksize = 64;
Int localTrrkComplexFloatBlocksize = 64;
Int localTrrkComplexDoubleBlocksize = 64;

// Qt5
ColorMap colorMap=RED_BLACK_GREEN;
Int numDiscreteColors = 15;
#ifdef ELEM_HAVE_QT5
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

namespace elem {

void PrintVersion( std::ostream& os )
{
    os << "Elemental version information:\n"
       << "  Git revision: " << ELEM_GIT_SHA1 << "\n"
       << "  Version:      " << Elemental_VERSION_MAJOR << "."
                             << Elemental_VERSION_MINOR << "\n"
       << "  Build type:   " << ELEM_CMAKE_BUILD_TYPE << "\n"
       << std::endl;
}

void PrintConfig( std::ostream& os )
{
    os << "Elemental configuration:\n"
       << "  Math libraries:               " << ELEM_MATH_LIBS << "\n"
       << "  Have FLAME bidiagonal SVD:    " 
#ifdef ELEM_HAVE_FLA_BSVD
       << "YES\n"
#else
       << "NO\n"
#endif
       << "  Have OpenMP:                  "
#ifdef ELEM_HAVE_OPENMP
       << "YES\n"
#else
       << "NO\n"
#endif
       << "  Have Qt5:                     "
#ifdef ELEM_HAVE_QT5
       << "YES\n"
#else
       << "NO\n"
#endif
       << "  Have F90 interface:           "
#ifdef ELEM_HAVE_F90_INTERFACE
       << "YES\n"
#else
       << "NO\n"
#endif
       << "  Avoiding complex MPI:         "
#ifdef ELEM_AVOID_COMPLEX_MPI
       << "YES\n"
#else
       << "NO\n"
#endif
       << "  Have MPI_Reducescatter_block: "
#ifdef ELEM_HAVE_MPI_REDUCE_SCATTER_BLOCK
       << "YES\n"
#else
       << "NO\n"
#endif
       << "  Have MPI_IN_PLACE:            "
#ifdef ELEM_HAVE_MPI_IN_PLACE
       << "YES\n"
#else
       << "NO\n"
#endif
       << "  AllReduce ReduceScatterBlock: "
#ifdef ELEM_REDUCE_SCATTER_BLOCK_VIA_ALLREDUCE
       << "YES\n"
#else
       << "NO\n"
#endif
       << "  Use byte Allgathers:          "
#ifdef ELEM_USE_BYTE_ALLGATHERS
       << "YES\n"
#else
       << "NO\n"
#endif
       << std::endl;
}

void PrintCCompilerInfo( std::ostream& os )
{
    os << "Elemental's C compiler info:\n"
       << "  ELEM_CMAKE_C_COMPILER:    " << ELEM_CMAKE_C_COMPILER << "\n"
       << "  ELEM_MPI_C_COMPILER:      " << ELEM_MPI_C_COMPILER << "\n"
       << "  ELEM_MPI_C_INCLUDE_PATH:  " << ELEM_MPI_C_INCLUDE_PATH << "\n"
       << "  ELEM_MPI_C_COMPILE_FLAGS: " << ELEM_MPI_C_COMPILE_FLAGS << "\n"
       << "  ELEM_MPI_C_LINK_FLAGS:    " << ELEM_MPI_C_LINK_FLAGS << "\n"
       << "  ELEM_MPI_C_LIBRARIES:     " << ELEM_MPI_C_LIBRARIES << "\n"
       << std::endl;
}

void PrintCxxCompilerInfo( std::ostream& os )
{
    os << "Elemental's C++ compiler info:\n"
       << "  ELEM_CMAKE_CXX_COMPILER:    " << ELEM_CMAKE_CXX_COMPILER << "\n"
       << "  ELEM_CXX_FLAGS:             " << ELEM_CXX_FLAGS << "\n"
       << "  ELEM_MPI_CXX_COMPILER:      " << ELEM_MPI_CXX_COMPILER << "\n"
       << "  ELEM_MPI_CXX_INCLUDE_PATH:  " << ELEM_MPI_CXX_INCLUDE_PATH << "\n"
       << "  ELEM_MPI_CXX_COMPILE_FLAGS: " << ELEM_MPI_CXX_COMPILE_FLAGS << "\n"
       << "  ELEM_MPI_CXX_LINK_FLAGS:    " << ELEM_MPI_CXX_LINK_FLAGS << "\n"
       << "  ELEM_MPI_CXX_LIBRARIES:     " << ELEM_MPI_CXX_LIBRARIES << "\n"
       << std::endl;
}

void SetColorMap( ColorMap map )
{ ::colorMap = map; }

ColorMap GetColorMap()
{ return ::colorMap; }

void SetNumDiscreteColors( Int numChunks )
{ ::numDiscreteColors = numChunks; }

Int NumDiscreteColors()
{ return ::numDiscreteColors; }

#ifdef ELEM_HAVE_QT5
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
        ::minRealWindowVal = std::min( ::minRealWindowVal, minVal );
    else
        ::minRealWindowVal = minVal;
    ::haveMinRealWindowVal = true;
}

void UpdateMaxRealWindowVal( double maxVal )
{
    if( ::haveMaxRealWindowVal )
        ::maxRealWindowVal = std::max( ::maxRealWindowVal, maxVal );
    else
        ::maxRealWindowVal = maxVal;
    ::haveMaxRealWindowVal = true;
}

void UpdateMinImagWindowVal( double minVal )
{
    if( ::haveMinImagWindowVal )
        ::minImagWindowVal = std::min( ::minImagWindowVal, minVal );
    else
        ::minImagWindowVal = minVal;
    ::haveMinImagWindowVal = true;
}

void UpdateMaxImagWindowVal( double maxVal )
{
    if( ::haveMaxImagWindowVal )
        ::maxImagWindowVal = std::max( ::maxImagWindowVal, maxVal );
    else
        ::maxImagWindowVal = maxVal;
    ::haveMaxImagWindowVal = true;
}
#endif // ifdef ELEM_HAVE_QT5

bool Initialized()
{ return ::numElemInits > 0; }

void Initialize( int& argc, char**& argv )
{
    if( ::numElemInits > 0 )
    {
        ++::numElemInits;
        return;
    }

    ::args = new Args( argc, argv );

    ::numElemInits = 1;
    if( !mpi::Initialized() )
    {
        if( mpi::Finalized() )
        {
            LogicError
            ("Cannot initialize elemental after finalizing MPI");
        }
#ifdef ELEM_HAVE_OPENMP
        const Int provided = 
            mpi::InitializeThread
            ( argc, argv, mpi::THREAD_MULTIPLE );
        const Int commRank = mpi::Rank( mpi::COMM_WORLD );
        if( provided != mpi::THREAD_MULTIPLE && commRank == 0 )
        {
            std::cerr << "WARNING: Could not achieve THREAD_MULTIPLE support."
                      << std::endl;
        }
#else
        mpi::Initialize( argc, argv );
#endif
        ::elemInitializedMpi = true;
    }
    else
    {
#ifdef ELEM_HAVE_OPENMP
        const Int provided = mpi::QueryThread();
        if( provided != mpi::THREAD_MULTIPLE )
        {
            throw std::runtime_error
            ("MPI initialized with inadequate thread support for Elemental");
        }
#endif
    }

#ifdef ELEM_HAVE_QT5
    ::coreApp = QCoreApplication::instance();
    if( ::coreApp == 0 )
    {
        // Test for whether the GUI should be disabled
        ::guiDisabled = false;
        for( int i=1; i<argc; ++i )
            if( !qstrcmp(argv[i],"-no-gui") )
                ::guiDisabled = true;
        if( ::guiDisabled )
            ::coreApp = new QCoreApplication( argc, argv );
        else
            ::coreApp = new QApplication( argc, argv );        
        ::elemInitializedQt = true;
    }
#endif

    // Queue a default algorithmic blocksize
    while( ! ::blocksizeStack.empty() )
        ::blocksizeStack.pop();
    ::blocksizeStack.push( 128 );

    // Build the default grid
    defaultGrid = new Grid( mpi::COMM_WORLD );

    // Create the types and ops needed for ValueInt
    mpi::CreateValueIntType<Int>();
    mpi::CreateValueIntType<float>();
    mpi::CreateValueIntType<double>();
    mpi::CreateMaxLocOp<Int>();
    mpi::CreateMaxLocOp<float>();
    mpi::CreateMaxLocOp<double>();
    mpi::CreateMinLocOp<Int>();
    mpi::CreateMinLocOp<float>();
    mpi::CreateMinLocOp<double>();

    // Do the same for ValueIntPair
    mpi::CreateValueIntPairType<Int>();
    mpi::CreateValueIntPairType<float>();
    mpi::CreateValueIntPairType<double>();
    mpi::CreateMaxLocPairOp<Int>();
    mpi::CreateMaxLocPairOp<float>();
    mpi::CreateMaxLocPairOp<double>();
    mpi::CreateMinLocPairOp<Int>();
    mpi::CreateMinLocPairOp<float>();
    mpi::CreateMinLocPairOp<double>();

    const unsigned rank = mpi::Rank( mpi::COMM_WORLD );
    // TODO: Allow for switching on/off reproducibility?
    //const long secs = time(NULL);
    const long secs = 21;
    const long seed = (secs<<16) | (rank & 0xFFFF);
    ::generator.seed( seed );
    srand( seed );
}

void Finalize()
{
    DEBUG_ONLY(CallStackEntry cse("Finalize"))
    if( ::numElemInits <= 0 )
        LogicError("Finalized Elemental more than initialized");
    --::numElemInits;

    if( mpi::Finalized() )
    {
        std::cerr << "Warning: MPI was finalized before Elemental." 
                  << std::endl;
    }
    if( ::numElemInits == 0 )
    {
        delete ::args;
        ::args = 0;

        if( ::elemInitializedMpi )
        {
            // Destroy the types and ops needed for ValueInt
            mpi::DestroyValueIntType<Int>();
            mpi::DestroyValueIntType<float>();
            mpi::DestroyValueIntType<double>();
            mpi::DestroyMaxLocOp<Int>();
            mpi::DestroyMaxLocOp<float>();
            mpi::DestroyMaxLocOp<double>();
            mpi::DestroyMinLocOp<Int>();
            mpi::DestroyMinLocOp<float>();
            mpi::DestroyMinLocOp<double>();

            // Do the same for ValueIntPair
            mpi::DestroyValueIntPairType<Int>();
            mpi::DestroyValueIntPairType<float>();
            mpi::DestroyValueIntPairType<double>();
            mpi::DestroyMaxLocPairOp<Int>();
            mpi::DestroyMaxLocPairOp<float>();
            mpi::DestroyMaxLocPairOp<double>();
            mpi::DestroyMinLocPairOp<Int>();
            mpi::DestroyMinLocPairOp<float>();
            mpi::DestroyMinLocPairOp<double>();

            // Delete the default grid
            delete ::defaultGrid;
            ::defaultGrid = 0;

            mpi::Finalize();
        }

#ifdef ELEM_HAVE_QT5
        if( ::elemInitializedQt )
        {
            if( ::elemOpenedWindow )
                ::coreApp->exec();
            else
                ::coreApp->exit();
            delete ::coreApp;
        }
#endif

        delete ::defaultGrid;
        ::defaultGrid = 0;
        while( ! ::blocksizeStack.empty() )
            ::blocksizeStack.pop();
    }
}

Args& GetArgs()
{ 
    if( args == 0 )
        throw std::runtime_error("No available instance of Args");
    return *::args; 
}

Int Blocksize()
{ return ::blocksizeStack.top(); }

void SetBlocksize( Int blocksize )
{ ::blocksizeStack.top() = blocksize; }

void PushBlocksizeStack( Int blocksize )
{ ::blocksizeStack.push( blocksize ); }

void PopBlocksizeStack()
{ ::blocksizeStack.pop(); }

const Grid& DefaultGrid()
{
    DEBUG_ONLY(
        CallStackEntry cse("DefaultGrid");
        if( ::defaultGrid == 0 )
            LogicError
            ("Attempted to return a non-existant default grid. Please ensure "
             "that Elemental is initialized before creating a DistMatrix.");
    )
    return *::defaultGrid;
}

Int DefaultBlockHeight()
{ return ::blockHeight; }

Int DefaultBlockWidth()
{ return ::blockWidth; }

void SetDefaultBlockHeight( Int mb )
{ ::blockHeight = mb; }

void SetDefaultBlockWidth( Int nb )
{ ::blockWidth = nb; }

std::mt19937& Generator()
{ return ::generator; }

// If we are not in RELEASE mode, then implement wrappers for a CallStack
DEBUG_ONLY(

    void PushCallStack( std::string s )
    { 
#ifdef ELEM_HAVE_OPENMP
        if( omp_get_thread_num() != 0 )
            return;
#endif // ELEM_HAVE_OPENMP
        ::callStack.push(s); 
    }

    void PopCallStack()
    { 
#ifdef ELEM_HAVE_OPENMP
        if( omp_get_thread_num() != 0 )
            return;
#endif // ELEM_HAVE_OPENMP
        ::callStack.pop(); 
    }

    void DumpCallStack( std::ostream& os )
    {
        std::ostringstream msg;
        while( ! ::callStack.empty() )
        {
            msg << "[" << ::callStack.size() << "]: " << ::callStack.top() 
                << "\n";
            ::callStack.pop();
        }
        os << msg.str();
        os.flush();
    }

) // DEBUG_ONLY

template<>
void SetLocalSymvBlocksize<float>( Int blocksize )
{ ::localSymvFloatBlocksize = blocksize; }

template<>
void SetLocalSymvBlocksize<double>( Int blocksize )
{ ::localSymvDoubleBlocksize = blocksize; }

template<>
void SetLocalSymvBlocksize<Complex<float>>( Int blocksize )
{ ::localSymvComplexFloatBlocksize = blocksize; }

template<>
void SetLocalSymvBlocksize<Complex<double>>( Int blocksize )
{ ::localSymvComplexDoubleBlocksize = blocksize; }

template<>
Int LocalSymvBlocksize<float>()
{ return ::localSymvFloatBlocksize; }

template<>
Int LocalSymvBlocksize<double>()
{ return ::localSymvDoubleBlocksize; }

template<>
Int LocalSymvBlocksize<Complex<float>>()
{ return ::localSymvComplexFloatBlocksize; }

template<>
Int LocalSymvBlocksize<Complex<double>>()
{ return ::localSymvComplexDoubleBlocksize; }

template<>
void SetLocalTrr2kBlocksize<float>( Int blocksize )
{ ::localTrr2kFloatBlocksize = blocksize; }

template<>
void SetLocalTrr2kBlocksize<double>( Int blocksize )
{ ::localTrr2kDoubleBlocksize = blocksize; }

template<>
void SetLocalTrr2kBlocksize<Complex<float>>( Int blocksize )
{ ::localTrr2kComplexFloatBlocksize = blocksize; }

template<>
void SetLocalTrr2kBlocksize<Complex<double>>( Int blocksize )
{ ::localTrr2kComplexDoubleBlocksize = blocksize; }

template<>
Int LocalTrr2kBlocksize<float>()
{ return ::localTrr2kFloatBlocksize; }

template<>
Int LocalTrr2kBlocksize<double>()
{ return ::localTrr2kDoubleBlocksize; }

template<>
Int LocalTrr2kBlocksize<Complex<float>>()
{ return ::localTrr2kComplexFloatBlocksize; }

template<>
Int LocalTrr2kBlocksize<Complex<double>>()
{ return ::localTrr2kComplexDoubleBlocksize; }

template<>
void SetLocalTrrkBlocksize<float>( Int blocksize )
{ ::localTrrkFloatBlocksize = blocksize; }

template<>
void SetLocalTrrkBlocksize<double>( Int blocksize )
{ ::localTrrkDoubleBlocksize = blocksize; }

template<>
void SetLocalTrrkBlocksize<Complex<float>>( Int blocksize )
{ ::localTrrkComplexFloatBlocksize = blocksize; }

template<>
void SetLocalTrrkBlocksize<Complex<double>>( Int blocksize )
{ ::localTrrkComplexDoubleBlocksize = blocksize; }

template<>
Int LocalTrrkBlocksize<float>()
{ return ::localTrrkFloatBlocksize; }

template<>
Int LocalTrrkBlocksize<double>()
{ return ::localTrrkDoubleBlocksize; }

template<>
Int LocalTrrkBlocksize<Complex<float>>()
{ return ::localTrrkComplexFloatBlocksize; }

template<>
Int LocalTrrkBlocksize<Complex<double>>()
{ return ::localTrrkComplexDoubleBlocksize; }

} // namespace elem
