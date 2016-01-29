/*
   Copyright (c) 2009-2016, Jack Poulson
                      2013, Jeff Hammond
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#include <algorithm>
#include <iomanip>
#include <set>
#include <stack>

#ifdef EL_HAVE_QT5
 #include <QApplication>
#endif

namespace {
using namespace El;

Int numElemInits = 0;
bool elemInitializedMpi = false;
#ifdef EL_HAVE_QD
unsigned oldControlWord=0;
#endif

std::stack<Int> blocksizeStack;
Grid* defaultGrid = 0;
Args* args = 0;

// Default blocksizes for BlockMatrix
Int blockHeight=32, blockWidth=32;

// A common Mersenne twister configuration
std::mt19937 generator;

#ifdef EL_HAVE_MPC
gmp_randstate_t gmpRandState;
#endif

// Debugging
DEBUG_ONLY(
  std::stack<string> callStack;
  bool tracingEnabled = false;
)

// A (per-process) output file for logging
std::ofstream logFile;

// Output/logging
Int indentLevel=0;
Int spacesPerIndent=2;

// Tuning parameters for basic routines
Int localSymvIntBlocksize = 64;
Int localSymvFloatBlocksize = 64;
Int localSymvDoubleBlocksize = 64;
Int localSymvComplexFloatBlocksize = 64;
Int localSymvComplexDoubleBlocksize = 64;
#ifdef EL_HAVE_QD
Int localSymvDoubleDoubleBlocksize = 64;
Int localSymvQuadDoubleBlocksize = 64;
#endif
#ifdef EL_HAVE_QUAD
Int localSymvQuadBlocksize = 64;
Int localSymvComplexQuadBlocksize = 64;
#endif
#ifdef EL_HAVE_MPC
Int localSymvBigIntBlocksize = 64;
Int localSymvBigFloatBlocksize = 64;
#endif

Int localTrr2kIntBlocksize = 64;
Int localTrr2kFloatBlocksize = 64;
Int localTrr2kDoubleBlocksize = 64;
Int localTrr2kComplexFloatBlocksize = 64;
Int localTrr2kComplexDoubleBlocksize = 64;
#ifdef EL_HAVE_QD
Int localTrr2kDoubleDoubleBlocksize = 64;
Int localTrr2kQuadDoubleBlocksize = 64;
#endif
#ifdef EL_HAVE_QUAD
Int localTrr2kQuadBlocksize = 64;
Int localTrr2kComplexQuadBlocksize = 64;
#endif
#ifdef EL_HAVE_MPC
Int localTrr2kBigIntBlocksize = 64;
Int localTrr2kBigFloatBlocksize = 64;
#endif

Int localTrrkIntBlocksize = 64;
Int localTrrkFloatBlocksize = 64;
Int localTrrkDoubleBlocksize = 64;
Int localTrrkComplexFloatBlocksize = 64;
Int localTrrkComplexDoubleBlocksize = 64;
#ifdef EL_HAVE_QD
Int localTrrkDoubleDoubleBlocksize = 64;
Int localTrrkQuadDoubleBlocksize = 64;
#endif
#ifdef EL_HAVE_QUAD
Int localTrrkQuadBlocksize = 64;
Int localTrrkComplexQuadBlocksize = 64;
#endif
#ifdef EL_HAVE_MPC
Int localTrrkBigIntBlocksize = 64;
Int localTrrkBigFloatBlocksize = 64;
#endif

// Qt5
// ===
// TODO: Move this into its own file?
ColorMap colorMap=RED_BLACK_GREEN;
Int numDiscreteColors = 15;
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

void PrintVersion( ostream& os )
{
    os << "Elemental version information:\n"
       << "  Git revision: " << EL_GIT_SHA1 << "\n"
       << "  Version:      " << EL_VERSION_MAJOR << "."
                             << EL_VERSION_MINOR << "\n"
       << "  Build type:   " << EL_CMAKE_BUILD_TYPE << "\n"
       << endl;
}

void PrintConfig( ostream& os )
{
    os << 
      "Elemental configuration:\n" <<
      "  Math libraries:               " << EL_MATH_LIBS << "\n"
#ifdef EL_HAVE_FLA_BSVD
      "  Have FLAME bidiagonal SVD:    YES\n"
#else
      "  Have FLAME bidiagonal SVD:    NO\n"
#endif
#ifdef EL_HYBRID
      "  Hybrid mode:                  YES\n"
#else
      "  Hybrid mode:                  NO\n"
#endif
#ifdef EL_HAVE_QT5
      "  Have Qt5:                     YES\n"
#else
      "  Have Qt5:                     NO\n"
#endif
#ifdef EL_AVOID_COMPLEX_MPI
      "  Avoiding complex MPI:         YES\n"
#else
      "  Avoiding complex MPI:         NO\n"
#endif
#ifdef EL_HAVE_MPI_REDUCE_SCATTER_BLOCK
      "  Have MPI_Reducescatter_block: YES\n"
#else
      "  Have MPI_Reducescatter_block: NO\n"
#endif
#ifdef EL_REDUCE_SCATTER_BLOCK_VIA_ALLREDUCE
      "  AllReduce ReduceScatterBlock: YES\n"
#else
      "  AllReduce ReduceScatterBlock: NO\n"
#endif
#ifdef EL_USE_BYTE_ALLGATHERS
      "  Use byte AllGathers:          YES\n"
#else
      "  Use byte AllGathers:          NO\n"
#endif
       << endl;
}

void PrintCCompilerInfo( ostream& os )
{
    os << "Elemental's C compiler info:\n"
       << "  EL_CMAKE_C_COMPILER:    " << EL_CMAKE_C_COMPILER << "\n"
       << "  EL_MPI_C_COMPILER:      " << EL_MPI_C_COMPILER << "\n"
       << "  EL_MPI_C_INCLUDE_PATH:  " << EL_MPI_C_INCLUDE_PATH << "\n"
       << "  EL_MPI_C_COMPILE_FLAGS: " << EL_MPI_C_COMPILE_FLAGS << "\n"
       << "  EL_MPI_LINK_FLAGS:      " << EL_MPI_LINK_FLAGS << "\n"
       << "  EL_MPI_C_LIBRARIES:     " << EL_MPI_C_LIBRARIES << "\n"
       << endl;
}

void PrintCxxCompilerInfo( ostream& os )
{
    os << "Elemental's C++ compiler info:\n"
       << "  EL_CMAKE_CXX_COMPILER:    " << EL_CMAKE_CXX_COMPILER << "\n"
       << "  EL_CXX_FLAGS:             " << EL_CXX_FLAGS << "\n"
       << "  EL_MPI_CXX_COMPILER:      " << EL_MPI_CXX_COMPILER << "\n"
       << "  EL_MPI_CXX_INCLUDE_PATH:  " << EL_MPI_CXX_INCLUDE_PATH << "\n"
       << "  EL_MPI_CXX_COMPILE_FLAGS: " << EL_MPI_CXX_COMPILE_FLAGS << "\n"
       << "  EL_MPI_LINK_FLAGS:        " << EL_MPI_LINK_FLAGS << "\n"
       << "  EL_MPI_CXX_LIBRARIES:     " << EL_MPI_CXX_LIBRARIES << "\n"
       << endl;
}

bool Using64BitInt()
{
#ifdef EL_USE_64BIT_INTS
    return true;
#else
    return false;
#endif
}

bool Using64BitBlasInt()
{
#ifdef EL_USE_64BIT_BLAS_INTS
    return true;
#else
    return false;
#endif
}

void SetColorMap( ColorMap map )
{ ::colorMap = map; }

ColorMap GetColorMap()
{ return ::colorMap; }

void SetNumDiscreteColors( Int numChunks )
{ ::numDiscreteColors = numChunks; }

Int NumDiscreteColors()
{ return ::numDiscreteColors; }

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
#endif // ifdef EL_HAVE_QT5

bool Initialized()
{ return ::numElemInits > 0; }

void Initialize()
{
    int argc=0;
    char** argv=NULL;
    Initialize( argc, argv );
}

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
#ifdef EL_HYBRID
        const Int provided = 
            mpi::InitializeThread
            ( argc, argv, mpi::THREAD_MULTIPLE );
        const int commRank = mpi::Rank( mpi::COMM_WORLD );
        if( provided != mpi::THREAD_MULTIPLE && commRank == 0 )
        {
            cerr << "WARNING: Could not achieve THREAD_MULTIPLE support."
                 << endl;
        }
#else
        mpi::Initialize( argc, argv );
#endif
        ::elemInitializedMpi = true;
    }
    else
    {
#ifdef EL_HYBRID
        const Int provided = mpi::QueryThread();
        if( provided != mpi::THREAD_MULTIPLE )
        {
            throw std::runtime_error
            ("MPI initialized with inadequate thread support for Elemental");
        }
#endif
    }

#ifdef EL_HAVE_QT5
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
#endif

    // Queue a default algorithmic blocksize
    while( ! ::blocksizeStack.empty() )
        ::blocksizeStack.pop();
    ::blocksizeStack.push( 128 );

    // Build the default grid
    defaultGrid = new Grid( mpi::COMM_WORLD );

#ifdef EL_HAVE_QD
    // TODO: Enable and disable when entering routines using QD
    fpu_fix_start( &::oldControlWord );
#endif

    const unsigned rank = mpi::Rank( mpi::COMM_WORLD );
    // TODO: Allow for switching on/off reproducibility?
    //const long secs = time(NULL);
    const long secs = 21;
    const long seed = (secs<<16) | (rank & 0xFFFF);
    ::generator.seed( seed );
    srand( seed );
#ifdef EL_HAVE_MPC
    mpc::SetMinIntBits( 256 );
    mpc::SetPrecision( 256 );
    gmp_randinit_default( ::gmpRandState );
    gmp_randseed_ui( ::gmpRandState, seed );
#endif

    // Create the types and ops
    // NOTE: mpc::SetPrecision created the BigFloat types
    mpi::CreateCustom();
}

void Finalize()
{
    DEBUG_ONLY(CSE cse("Finalize"))
    if( ::numElemInits <= 0 )
    { 
        cerr << "Finalized Elemental more times than initialized" << endl;
        return;
    }
    --::numElemInits;

    if( mpi::Finalized() )
        cerr << "Warning: MPI was finalized before Elemental." << endl;
    if( ::numElemInits == 0 )
    {
        delete ::args;
        ::args = 0;
       
        // Destroy the types and ops
        mpi::DestroyCustom();

        // Delete the default grid
        delete ::defaultGrid;
        ::defaultGrid = 0;

#ifdef EL_HAVE_QT5
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
#endif
        if( ::elemInitializedMpi )
            mpi::Finalize();


        while( ! ::blocksizeStack.empty() )
            ::blocksizeStack.pop();

#ifdef EL_HAVE_QD
        // TODO: Enable and disable when entering routines using QD
        fpu_fix_end( &::oldControlWord );
#endif

#ifdef EL_HAVE_MPC
        gmp_randclear( ::gmpRandState );
#endif
    }

    DEBUG_ONLY( CloseLog() )
}

Args& GetArgs()
{ 
    if( args == 0 )
        throw std::runtime_error("No available instance of Args");
    return *::args; 
}

Int Blocksize()
{ 
    DEBUG_ONLY(
      if( ::blocksizeStack.empty() )
          LogicError("Attempted to extract blocksize from empty stack");
    )
    return ::blocksizeStack.top(); 
}

void SetBlocksize( Int blocksize )
{ 
    DEBUG_ONLY(
      if( ::blocksizeStack.empty() )
          LogicError("Attempted to set blocksize at top of empty stack");
    )
    ::blocksizeStack.top() = blocksize; 
}

void PushBlocksizeStack( Int blocksize )
{ ::blocksizeStack.push( blocksize ); }

void PopBlocksizeStack()
{
    DEBUG_ONLY(
      if( ::blocksizeStack.empty() )
          LogicError("Attempted to pop an empty blocksize stack");
    )
    ::blocksizeStack.pop();
}

const Grid& DefaultGrid() EL_NO_RELEASE_EXCEPT
{
    DEBUG_ONLY(
      CSE cse("DefaultGrid");
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

#ifdef EL_HAVE_MPC
namespace mpc {

void RandomState( gmp_randstate_t randState )
{
    // It is surprisingly tedious to return the state...

    randState->_mp_seed->_mp_alloc = ::gmpRandState->_mp_seed->_mp_alloc;
    randState->_mp_seed->_mp_size = ::gmpRandState->_mp_seed->_mp_size;
    randState->_mp_seed->_mp_d = ::gmpRandState->_mp_seed->_mp_d;

    randState->_mp_alg = ::gmpRandState->_mp_alg;
    randState->_mp_algdata._mp_lc = ::gmpRandState->_mp_algdata._mp_lc;
}

} // namespace mpc
#endif

void Args::HandleVersion( ostream& os ) const
{
    string version = "--version";
    char** arg = std::find( argv_, argv_+argc_, version );
    const bool foundVersion = ( arg != argv_+argc_ );
    if( foundVersion )
    {
        if( mpi::Rank() == 0 )
            PrintVersion();
        throw ArgException();
    }
}

void Args::HandleBuild( ostream& os ) const
{
    string build = "--build";
    char** arg = std::find( argv_, argv_+argc_, build );
    const bool foundBuild = ( arg != argv_+argc_ );
    if( foundBuild )
    {
        if( mpi::Rank() == 0 )
        {
            PrintVersion();
            PrintConfig();
            PrintCCompilerInfo();
            PrintCxxCompilerInfo();
        }
        throw ArgException();
    }
}

void ReportException( const exception& e, ostream& os )
{
    try
    {
        const ArgException& argExcept = dynamic_cast<const ArgException&>(e);
        if( string(argExcept.what()) != "" )
            os << argExcept.what() << endl;
        DEBUG_ONLY(DumpCallStack(os))
    }
    catch( UnrecoverableException& recovExcept )
    {
        if( string(e.what()) != "" )
        {
            os << "Process " << mpi::Rank() 
               << " caught an unrecoverable exception with message:\n"
               << e.what() << endl;
        }
        DEBUG_ONLY(DumpCallStack(os))
        mpi::Abort( mpi::COMM_WORLD, 1 );
    }
    catch( exception& castExcept ) 
    { 
        if( string(e.what()) != "" )
        {
            os << "Process " << mpi::Rank() << " caught error message:\n"
               << e.what() << endl;
        }
        DEBUG_ONLY(DumpCallStack(os))
    }
}

void ComplainIfDebug()
{
    DEBUG_ONLY(
        if( mpi::Rank() == 0 )
        {
            Output("=======================================================");
            Output(" In debug mode! Do not expect competitive performance! ");
            Output("=======================================================");
        }
    )
}

// If we are not in RELEASE mode, then implement wrappers for a CallStack
DEBUG_ONLY(

    void EnableTracing() { ::tracingEnabled = true; }
    void DisableTracing() { ::tracingEnabled = false; }

    void PushCallStack( string s )
    { 
#ifdef EL_HYBRID
        if( omp_get_thread_num() != 0 )
            return;
#endif
        const size_t maxStackSize = 300;
        if( ::callStack.size() > maxStackSize )
        {
            DumpCallStack();
            return;
        }
        ::callStack.push(s); 
        if( ::tracingEnabled )
        {
            const int stackSize = ::callStack.size();
            ostringstream os;
            for( int j=0; j<stackSize; ++j )
                os << " "; 
            os << s << endl;
            cout <<  os.str();
        }
    }

    void PopCallStack()
    { 
#ifdef EL_HYBRID
        if( omp_get_thread_num() != 0 )
            return;
#endif
        if( ::callStack.empty() )
            LogicError("Attempted to pop an empty call stack");
        ::callStack.pop(); 
    }

    void DumpCallStack( ostream& os )
    {
        ostringstream msg;
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

void OpenLog( const char* filename )
{
    if( ::logFile.is_open() )
        CloseLog();
    ::logFile.open( filename );
}

std::ostream& LogOS()
{
    if( !::logFile.is_open() )
    {
        ostringstream fileOS;
        fileOS << "El-Proc" << std::setfill('0') << std::setw(3)
               << mpi::Rank() << ".log";
        ::logFile.open( fileOS.str().c_str() );
    }
    return ::logFile; 
}

void CloseLog() { ::logFile.close(); }

Int PushIndent() { return ::indentLevel++; }
Int PopIndent() { return ::indentLevel--; }
void SetIndent( Int indent ) { ::indentLevel = indent; }
void ClearIndent() { ::indentLevel = 0; }
Int IndentLevel() { return ::indentLevel; }

string Indent()
{
    string ind;
    for( Int i=0; i < ::spacesPerIndent * ::indentLevel; ++i )
        ind = ind + " ";
    return ind;
}

template<>
void SetLocalSymvBlocksize<Int>( Int blocksize )
{ ::localSymvIntBlocksize = blocksize; }

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

#ifdef EL_HAVE_QD
template<>
void SetLocalSymvBlocksize<DoubleDouble>( Int blocksize )
{ ::localSymvDoubleDoubleBlocksize = blocksize; }
template<>
void SetLocalSymvBlocksize<QuadDouble>( Int blocksize )
{ ::localSymvQuadDoubleBlocksize = blocksize; }
#endif

#ifdef EL_HAVE_QUAD
template<>
void SetLocalSymvBlocksize<Quad>( Int blocksize )
{ ::localSymvQuadBlocksize = blocksize; }

template<>
void SetLocalSymvBlocksize<Complex<Quad>>( Int blocksize )
{ ::localSymvComplexQuadBlocksize = blocksize; }
#endif

#ifdef EL_HAVE_MPC
template<>
void SetLocalSymvBlocksize<BigInt>( Int blocksize )
{ ::localSymvBigIntBlocksize = blocksize; }

template<>
void SetLocalSymvBlocksize<BigFloat>( Int blocksize )
{ ::localSymvBigFloatBlocksize = blocksize; }
#endif

template<>
Int LocalSymvBlocksize<Int>()
{ return ::localSymvIntBlocksize; }

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

#ifdef EL_HAVE_QD
template<>
Int LocalSymvBlocksize<DoubleDouble>()
{ return ::localSymvDoubleDoubleBlocksize; }
template<>
Int LocalSymvBlocksize<QuadDouble>()
{ return ::localSymvQuadDoubleBlocksize; }
#endif

#ifdef EL_HAVE_QUAD
template<>
Int LocalSymvBlocksize<Quad>()
{ return ::localSymvQuadBlocksize; }

template<>
Int LocalSymvBlocksize<Complex<Quad>>()
{ return ::localSymvComplexQuadBlocksize; }
#endif

#ifdef EL_HAVE_MPC
template<>
Int LocalSymvBlocksize<BigInt>()
{ return ::localSymvBigIntBlocksize; }
template<>
Int LocalSymvBlocksize<BigFloat>()
{ return ::localSymvBigFloatBlocksize; }
#endif

template<>
void SetLocalTrr2kBlocksize<Int>( Int blocksize )
{ ::localTrr2kIntBlocksize = blocksize; }

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

#ifdef EL_HAVE_QD
template<>
void SetLocalTrr2kBlocksize<DoubleDouble>( Int blocksize )
{ ::localTrr2kDoubleDoubleBlocksize = blocksize; }

template<>
void SetLocalTrr2kBlocksize<QuadDouble>( Int blocksize )
{ ::localTrr2kQuadDoubleBlocksize = blocksize; }
#endif

#ifdef EL_HAVE_QUAD
template<>
void SetLocalTrr2kBlocksize<Quad>( Int blocksize )
{ ::localTrr2kQuadBlocksize = blocksize; }

template<>
void SetLocalTrr2kBlocksize<Complex<Quad>>( Int blocksize )
{ ::localTrr2kComplexQuadBlocksize = blocksize; }
#endif

#ifdef EL_HAVE_MPC
template<>
void SetLocalTrr2kBlocksize<BigInt>( Int blocksize )
{ ::localTrr2kBigIntBlocksize = blocksize; }
template<>
void SetLocalTrr2kBlocksize<BigFloat>( Int blocksize )
{ ::localTrr2kBigFloatBlocksize = blocksize; }
#endif

template<>
Int LocalTrr2kBlocksize<Int>()
{ return ::localTrr2kIntBlocksize; }

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

#ifdef EL_HAVE_QD
template<>
Int LocalTrr2kBlocksize<DoubleDouble>()
{ return ::localTrr2kDoubleDoubleBlocksize; }

template<>
Int LocalTrr2kBlocksize<QuadDouble>()
{ return ::localTrr2kQuadDoubleBlocksize; }
#endif

#ifdef EL_HAVE_QUAD
template<>
Int LocalTrr2kBlocksize<Quad>()
{ return ::localTrr2kQuadBlocksize; }

template<>
Int LocalTrr2kBlocksize<Complex<Quad>>()
{ return ::localTrr2kComplexQuadBlocksize; }
#endif

#ifdef EL_HAVE_MPC
template<>
Int LocalTrr2kBlocksize<BigInt>()
{ return ::localTrr2kBigIntBlocksize; }
template<>
Int LocalTrr2kBlocksize<BigFloat>()
{ return ::localTrr2kBigFloatBlocksize; }
#endif

template<>
void SetLocalTrrkBlocksize<Int>( Int blocksize )
{ ::localTrrkIntBlocksize = blocksize; }

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

#ifdef EL_HAVE_QD
template<>
void SetLocalTrrkBlocksize<DoubleDouble>( Int blocksize )
{ ::localTrrkDoubleDoubleBlocksize = blocksize; }

template<>
void SetLocalTrrkBlocksize<QuadDouble>( Int blocksize )
{ ::localTrrkQuadDoubleBlocksize = blocksize; }
#endif

#ifdef EL_HAVE_QUAD
template<>
void SetLocalTrrkBlocksize<Quad>( Int blocksize )
{ ::localTrrkQuadBlocksize = blocksize; }

template<>
void SetLocalTrrkBlocksize<Complex<Quad>>( Int blocksize )
{ ::localTrrkComplexQuadBlocksize = blocksize; }
#endif

#ifdef EL_HAVE_MPC
template<>
void SetLocalTrrkBlocksize<BigInt>( Int blocksize )
{ ::localTrrkBigIntBlocksize = blocksize; }
template<>
void SetLocalTrrkBlocksize<BigFloat>( Int blocksize )
{ ::localTrrkBigFloatBlocksize = blocksize; }
#endif

template<>
Int LocalTrrkBlocksize<Int>()
{ return ::localTrrkIntBlocksize; }

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

#ifdef EL_HAVE_QD
template<>
Int LocalTrrkBlocksize<DoubleDouble>()
{ return ::localTrrkDoubleDoubleBlocksize; }

template<>
Int LocalTrrkBlocksize<QuadDouble>()
{ return ::localTrrkQuadDoubleBlocksize; }
#endif

#ifdef EL_HAVE_QUAD
template<>
Int LocalTrrkBlocksize<Quad>()
{ return ::localTrrkQuadBlocksize; }

template<>
Int LocalTrrkBlocksize<Complex<Quad>>()
{ return ::localTrrkComplexQuadBlocksize; }
#endif

#ifdef EL_HAVE_MPC
template<>
Int LocalTrrkBlocksize<BigInt>()
{ return ::localTrrkBigIntBlocksize; }
template<>
Int LocalTrrkBlocksize<BigFloat>()
{ return ::localTrrkBigFloatBlocksize; }
#endif

template<typename T>
bool IsSorted( const vector<T>& x )
{
    const Int vecLength = x.size();
    for( Int i=1; i<vecLength; ++i )
    {
        if( x[i] < x[i-1] )
            return false;
    }
    return true;
}

// While is_strictly_sorted exists in Boost, it does not exist in the STL (yet)
template<typename T>
bool IsStrictlySorted( const vector<T>& x )
{
    const Int vecLength = x.size();
    for( Int i=1; i<vecLength; ++i )
    {
        if( x[i] <= x[i-1] )
            return false;
    }
    return true;
}

void Union
( vector<Int>& both, const vector<Int>& first, const vector<Int>& second )
{
    both.resize( first.size()+second.size() );
    auto it = std::set_union
      ( first.cbegin(),  first.cend(), 
        second.cbegin(), second.cend(),
        both.begin() );
    both.resize( Int(it-both.begin()) );
}

vector<Int>
Union( const vector<Int>& first, const vector<Int>& second )
{
    vector<Int> both;
    Union( both, first, second );
    return both;
}

void RelativeIndices
( vector<Int>& relInds, const vector<Int>& sub, const vector<Int>& full )
{
    const Int numSub = sub.size();
    relInds.resize( numSub );
    auto it = full.cbegin();
    for( Int i=0; i<numSub; ++i )
    {
        const Int index = sub[i];
        it = std::lower_bound( it, full.cend(), index );
        DEBUG_ONLY(
          if( it == full.cend() )
              LogicError("Index was not found");
        )
        relInds[i] = Int(it-full.cbegin());
    }
}

vector<Int> RelativeIndices( const vector<Int>& sub, const vector<Int>& full )
{
    vector<Int> relInds;
    RelativeIndices( relInds, sub, full );
    return relInds;
}

Int Find( const vector<Int>& sortedInds, Int index )
{
    DEBUG_ONLY(CSE cse("Find"))
    auto it = std::lower_bound( sortedInds.cbegin(), sortedInds.cend(), index );
    DEBUG_ONLY(
      if( it == sortedInds.cend() )
          LogicError("All indices were smaller");
      if( *it != index )
          LogicError("Could not find index");
    )
    return it - sortedInds.cbegin();
}

#define EL_NO_COMPLEX_PROTO
#define PROTO(T) \
  template bool IsSorted( const vector<T>& x ); \
  template bool IsStrictlySorted( const vector<T>& x );
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
