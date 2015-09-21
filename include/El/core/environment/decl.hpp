/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_ENVIRONMENT_DECL_HPP
#define EL_ENVIRONMENT_DECL_HPP

namespace El {

using std::size_t;

using std::array;
using std::complex;
using std::deque;
using std::function;
using std::pair;
using std::set;
using std::vector;

using std::make_shared;
using std::shared_ptr;
using std::unique_ptr;

using std::move;

using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::ostream;
using std::ostringstream;

using std::exception;
using std::uncaught_exception;

void PrintVersion( ostream& os=cout );
void PrintConfig( ostream& os=cout );
void PrintCCompilerInfo( ostream& os=cout );
void PrintCxxCompilerInfo( ostream& os=cout );
bool Using64BitInt();
bool Using64BitBlasInt();

// For initializing and finalizing Elemental
void Initialize( int& argc, char**& argv );
void Finalize();
bool Initialized();

// For getting the MPI argument instance (for internal usage)
class Args : public choice::MpiArgs
{
public:
    Args
    ( int argc, char** argv,
      mpi::Comm comm=mpi::COMM_WORLD, ostream& error=cerr )
    : choice::MpiArgs(argc,argv,comm,error)
    { }
    virtual ~Args() { }
protected:
    virtual void HandleVersion( ostream& os=cout ) const;
    virtual void HandleBuild( ostream& os=cout ) const;
};
Args& GetArgs();

// For processing command-line arguments
template<typename T>
T Input( string name, string desc );
template<typename T>
T Input( string name, string desc, T defaultVal );
void ProcessInput();
void PrintInputReport();

// For getting and setting the algorithmic blocksize
Int Blocksize();
void SetBlocksize( Int blocksize );

// For manipulating the algorithmic blocksize as a stack
void PushBlocksizeStack( Int blocksize );
void PopBlocksizeStack();

Int DefaultBlockHeight();
Int DefaultBlockWidth();
void SetDefaultBlockHeight( Int blockHeight );
void SetDefaultBlockWidth( Int blockWidth );

std::mt19937& Generator();

template<typename T>
inline const T& Max( const T& m, const T& n ) EL_NO_EXCEPT
{ return std::max(m,n); }

inline const Int& Max( const Int& m, const Int& n ) EL_NO_EXCEPT
{ return std::max(m,n); }

template<typename T>
inline const T& Min( const T& m, const T& n ) EL_NO_EXCEPT
{ return std::min(m,n); }

inline const Int& Min( const Int& m, const Int& n ) EL_NO_EXCEPT
{ return std::min(m,n); }

// Replacement for std::memcpy, which is known to often be suboptimal.
// Notice the sizeof(T) is no longer required.
template<typename T>
void MemCopy( T* dest, const T* source, size_t numEntries );

template<typename T>
void MemSwap( T* a, T* b, T* temp, size_t numEntries );

// Generalization of std::memcpy so that unit strides are not required
template<typename T>
void StridedMemCopy
(       T* dest,   Int destStride,
  const T* source, Int sourceStride, Int numEntries );

template<typename S,typename T>
inline void CopySTL( const S& a, T& b )
{
    b.resize( a.size() ); 
    std::copy( a.begin(), a.end(), b.begin() );
}

// Replacement for std::memset, which is likely suboptimal and hard to extend
// to non-POD datatypes. Notice that sizeof(T) is no longer required.
template<typename T>
void MemZero( T* buffer, size_t numEntries );

// Clear the contents of x by swapping with an empty object of the same type
template<typename T>
void SwapClear( T& x );

inline void BuildStream( ostringstream& os ) { }

template<typename T,typename... Args>
inline void BuildStream( ostringstream& os, T item, Args... args )
{
    os << item;
    BuildStream( os, args... );
}

template<typename... Args>
inline void LogicError( Args... args )
{
    ostringstream os;
    BuildStream( os, args... );
    os << endl;
    throw std::logic_error( os.str().c_str() );
}

template<typename... Args>
inline void RuntimeError( Args... args )
{
    ostringstream os;
    BuildStream( os, args... );
    os << endl;
    throw std::logic_error( os.str().c_str() );
}

// This is the only place that Elemental is currently using duck-typing.
// I'm not sure if it's a good idea to use it more often.
template<class MatType>
inline string 
DimsString( const MatType& A, string label="Matrix" )
{ 
    ostringstream os;
    os << label << " ~ " << A.Height() << " x " << A.Width();
    return os.str();
}

// This is defined in choice.hpp
class ArgException;

// An exception which signifies that a matrix was unexpectedly singular.
class SingularMatrixException : public std::runtime_error 
{
public:
    SingularMatrixException( const char* msg="Matrix was singular" ) 
    : std::runtime_error( msg ) { }
};

// An exception which signifies a zero pivot was chosen, though the matrix
// may not actually be singular
class ZeroPivotException : public std::runtime_error
{
public:
    ZeroPivotException( const char* msg="Zero pivot was chosen" )
    : std::runtime_error( msg ) { }
};

// An exception which signifies that a matrix was unexpectedly non-HPD
class NonHPDMatrixException  : public std::runtime_error
{
public:
    NonHPDMatrixException( const char* msg="Matrix was not HPD" )
    : std::runtime_error( msg ) { }
};

// An exception which signifies that a matrix was unexpectedly non-HPSD
class NonHPSDMatrixException  : public std::runtime_error
{
public:
    NonHPSDMatrixException( const char* msg="Matrix was not HPSD" )
    : std::runtime_error( msg ) { }
};

DEBUG_ONLY(
    void PushCallStack( string s );
    void PopCallStack();
    void DumpCallStack( ostream& os=cerr );

    class CallStackEntry 
    {
    public:
        CallStackEntry( string s ) 
        { 
            if( !uncaught_exception() )
                PushCallStack(s); 
        }
        ~CallStackEntry() 
        { 
            if( !uncaught_exception() )
                PopCallStack(); 
        }
    };
    typedef CallStackEntry CSE;
)

#ifndef EL_RELEASE
    void LogFileOpen( char* filename );

    std::ostream & LogFileOS();

    void LogFileCoutStr( std::string str );

    template<typename... Args>
    inline void LogFileCout( Args... args )
    {
        std::ostringstream str;
        BuildStream( str, args... );
        LogFileCoutStr( str.str() );
    }

    void LogFileClose();
#endif

void ReportException( const exception& e, ostream& os=cerr );

void ComplainIfDebug();

Int PushIndent();
Int PopIndent();
void SetIndent( Int level );
void ClearIndent();
Int IndentLevel();
std::string Indent();

template<typename... Args>
inline void Output( Args... args )
{
    ostringstream os;
    os << Indent();
    BuildStream( os, args... );
    os << endl;
    cout << os.str();
}

// TODO: OutputRoot?

template<typename T>
void EnsureConsistent( T alpha, mpi::Comm comm, string name="" );

// This will be guaranteed by C++14 via std::make_unique
template<typename T, typename ...Args>
inline unique_ptr<T> MakeUnique( Args&& ...args )
{ return unique_ptr<T>( new T( std::forward<Args>(args)... ) ); }

template<typename T>
T Scan( const vector<T>& counts, vector<T>& offsets );

template<typename T>
bool IsSorted( const vector<T>& x );
// While is_strictly_sorted exists in Boost, it does not exist in the STL (yet)
template<typename T>
bool IsStrictlySorted( const vector<T>& x );

void Union
( vector<Int>& both,
  const vector<Int>& first, const vector<Int>& second );
vector<Int>
Union( const vector<Int>& first, const vector<Int>& second );

void RelativeIndices
( vector<Int>& relInds, const vector<Int>& sub, const vector<Int>& full );
vector<Int> RelativeIndices( const vector<Int>& sub, const vector<Int>& full );

Int Find( const vector<Int>& sortedInds, Int index );
Int Find( const vector<Int>& sortedInds, Int index, const string& msg );

template<typename F>
void UpdateScaledSquare
( F alpha, Base<F>& scale, Base<F>& scaledSquare ) EL_NO_EXCEPT;

} // namespace El

#endif // ifndef EL_ENVIRONMENT_DECL_HPP
