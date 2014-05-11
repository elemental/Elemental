/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_ENVIRONMENT_DECL_HPP
#define ELEM_ENVIRONMENT_DECL_HPP

namespace elem {

void PrintVersion( std::ostream& os=std::cout );
void PrintConfig( std::ostream& os=std::cout );
void PrintCCompilerInfo( std::ostream& os=std::cout );
void PrintCxxCompilerInfo( std::ostream& os=std::cout );

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
      mpi::Comm comm=mpi::COMM_WORLD, std::ostream& error=std::cerr )
    : choice::MpiArgs(argc,argv,comm,error)
    { }
    virtual ~Args() { }
protected:
    virtual void HandleVersion( std::ostream& os=std::cout ) const;
    virtual void HandleBuild( std::ostream& os=std::cout ) const;
};
Args& GetArgs();

// For processing command-line arguments
template<typename T>
T Input( std::string name, std::string desc );
template<typename T>
T Input( std::string name, std::string desc, T defaultVal );
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
inline T Max( T m, T n )
{ return std::max(m,n); }

inline Int Max( Int m, Int n )
{ return std::max(m,n); }

template<typename T>
inline T Min( T m, T n )
{ return std::min(m,n); }

inline Int Min( Int m, Int n )
{ return std::min(m,n); }

// Replacement for std::memcpy, which is known to often be suboptimal.
// Notice the sizeof(T) is no longer required.
template<typename T>
void MemCopy( T* dest, const T* source, std::size_t numEntries );

template<typename T>
void MemSwap( T* a, T* b, T* temp, std::size_t numEntries );

// Generalization of std::memcpy so that unit strides are not required
template<typename T>
void StridedMemCopy
(       T* dest,   std::size_t destStride,
  const T* source, std::size_t sourceStride, std::size_t numEntries );

// Replacement for std::memset, which is likely suboptimal and hard to extend
// to non-POD datatypes. Notice that sizeof(T) is no longer required.
template<typename T>
void MemZero( T* buffer, std::size_t numEntries );

// Clear the contents of x by swapping with an empty object of the same type
template<typename T>
void SwapClear( T& x );

inline void BuildStream( std::ostringstream& os ) { }

template<typename T,typename... Args>
inline void BuildStream( std::ostringstream& os, T item, Args... args )
{
    os << item;
    BuildStream( os, args... );
}

template<typename... Args>
inline void LogicError( Args... args )
{
    std::ostringstream os;
    BuildStream( os, args... );
    os << std::endl;
    throw std::logic_error( os.str().c_str() );
}

template<typename... Args>
inline void RuntimeError( Args... args )
{
    std::ostringstream os;
    BuildStream( os, args... );
    os << std::endl;
    throw std::logic_error( os.str().c_str() );
}

// This is the only place that Elemental is currently using duck-typing.
// I'm not sure if it's a good idea to use it more often.
template<class MatType>
inline std::string 
DimsString( const MatType& A, std::string label="Matrix" )
{ 
    std::ostringstream os;
    os << label << " ~ " << A.Height() << " x " << A.Width();
    return os.str();
}

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
    void PushCallStack( std::string s );
    void PopCallStack();
    void DumpCallStack( std::ostream& os=std::cerr );

    class CallStackEntry 
    {
    public:
        CallStackEntry( std::string s ) 
        { 
            if( !std::uncaught_exception() )
                PushCallStack(s); 
        }
        ~CallStackEntry() 
        { 
            if( !std::uncaught_exception() )
                PopCallStack(); 
        }
    };
)

void ReportException( const std::exception& e, std::ostream& os=std::cerr );
class ArgException;

void ComplainIfDebug();

template<typename T>
void EnsureConsistent( T alpha, mpi::Comm comm, std::string name="" );

} // namespace elem

#endif // ifndef ELEM_ENVIRONMENT_DECL_HPP
