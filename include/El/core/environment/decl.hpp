/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_ENVIRONMENT_DECL_HPP
#define EL_ENVIRONMENT_DECL_HPP

namespace El {

using std::size_t;

using std::array;
using std::function;
using std::pair;
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
using std::istream;
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

// For manually initializing and finalizing Elemental; their direct usage
// in C++ programs now deprecated.
void Initialize();
void Initialize( int& argc, char**& argv );
void Finalize();
bool Initialized();

// For initializing/finalizing Elemental using RAII
class Environment
{
public:
    Environment() { Initialize(); }
    Environment( int& argc, char**& argv ) { Initialize( argc, argv ); }
    ~Environment() { Finalize(); }
};

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
void EmptyBlocksizeStack();

template<typename T,
         typename=EnableIf<IsScalar<T>>>
const T& Max( const T& m, const T& n ) EL_NO_EXCEPT;
inline const Int& Max( const Int& m, const Int& n ) EL_NO_EXCEPT;

template<typename T,
         typename=EnableIf<IsScalar<T>>>
const T& Min( const T& m, const T& n ) EL_NO_EXCEPT;
inline const Int& Min( const Int& m, const Int& n ) EL_NO_EXCEPT;

// Replacement for std::memcpy, which is known to often be suboptimal.
// Notice the sizeof(T) is no longer required.
template<typename T,
         typename=EnableIf<IsPacked<T>>>
void MemCopy
(       T* dest,
  const T* source,
        size_t numEntries );
template<typename T,
         typename=DisableIf<IsPacked<T>>,
         typename=void>
void MemCopy
(       T* dest,
  const T* source,
        size_t numEntries );

template<typename T,
         typename=EnableIf<IsPacked<T>>>
void MemSwap
( T* a,
  T* b,
  T* temp,
  size_t numEntries );
template<typename T,
         typename=DisableIf<IsPacked<T>>,
         typename=void>
void MemSwap
( T* a,
  T* b,
  T* temp,
  size_t numEntries );

// Generalization of std::memcpy so that unit strides are not required
template<typename T,
         typename=EnableIf<IsPacked<T>>>
void StridedMemCopy
(       T* dest,   Int destStride,
  const T* source, Int sourceStride, Int numEntries );
template<typename T,
         typename=DisableIf<IsPacked<T>>,
         typename=void>
void StridedMemCopy
(       T* dest,   Int destStride,
  const T* source, Int sourceStride, Int numEntries );

// A thin wrapper around std::copy
template<typename S,typename T>
void CopySTL( const S& a, T& b );

// Replacement for std::memset, which is likely suboptimal and hard to extend
// to non-POD datatypes. Notice that sizeof(T) is no longer required.
template<typename T,
         typename=EnableIf<IsPacked<T>>>
void MemZero( T* buffer, size_t numEntries );
template<typename T,
         typename=DisableIf<IsPacked<T>>,
         typename=void>
void MemZero( T* buffer, size_t numEntries );

// Clear the contents of x by swapping with an empty object of the same type
template<typename T>
void SwapClear( T& x );

// Reserve memory in a vector without zero-initializing the variables unless
// valgrind is currently running or the datatype *requires* construction.
template<typename T,
         typename=EnableIf<IsPacked<T>>>
void FastResize( vector<T>& v, Int numEntries );
template<typename T,
         typename=DisableIf<IsPacked<T>>,
         typename=void>
void FastResize( vector<T>& v, Int numEntries );

inline void BuildStream( ostringstream& os ) { }

template<typename T,typename... ArgPack>
void BuildStream( ostringstream& os, const T& item, const ArgPack& ... args );
template<typename... ArgPack>
string BuildString( const ArgPack& ... args );

class UnrecoverableException : public std::runtime_error
{
public:
    UnrecoverableException( const char* msg="Unrecoverable exception" )
    : std::runtime_error( msg ) { }
};

template<typename... ArgPack>
void UnrecoverableError( const ArgPack& ... args );
template<typename... ArgPack>
void LogicError( const ArgPack& ... args );
template<typename... ArgPack>
void RuntimeError( const ArgPack& ... args );

// This is the only place that Elemental is currently using duck-typing.
// I'm not sure if it's a good idea to use it more often.
template<class MatType>
string DimsString( const MatType& A, string label="Matrix" );

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

EL_DEBUG_ONLY(
    void EnableTracing();
    void DisableTracing();

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

void OpenLog( const char* filename );

std::ostream & LogOS();

template<typename... ArgPack>
void Log( const ArgPack& ... args );

void CloseLog();

void ReportException( const exception& e, ostream& os=cerr );

void ComplainIfDebug();

Int PushIndent();
Int PopIndent();
void SetIndent( Int level );
void ClearIndent();
Int IndentLevel();
std::string Indent();

template<typename... ArgPack>
void Output( const ArgPack& ... args );

template<typename... ArgPack>
void OutputFromRoot( mpi::Comm comm, const ArgPack& ... args );

template<typename T>
void EnsureConsistent( T alpha, mpi::Comm comm, string name="" );

// This will be guaranteed by C++14 via std::make_unique
template<typename T,typename ...ArgPack>
unique_ptr<T> MakeUnique( ArgPack&& ...args )
{ return unique_ptr<T>( new T( std::forward<ArgPack>(args)... ) ); }

// MakeFunction allows for the automatic conversion of a lambda to an
// std::function and is similar to http://stackoverflow.com/a/24068396/1119818.

/*
namespace make_function {

template<typename T>
struct Helper { using type = void; };
template<typename Ret,typename Class,typename... Args>
struct Helper<Ret(Class::*)(Args...) const>
{ using type = std::function<Ret(Args...)>; };

} // namespace make_function

template<typename Function>
typename make_function::Helper<decltype(&Function::operator())>::type
MakeFunction(Function const& func) { return func; }
*/

// Handles generic types that are functors, delegate to its 'operator()'
template<typename T>
struct function_traits : public function_traits<decltype(&T::operator())>
{};

// Handles pointers to member functions
template<typename ClassType,typename ReturnType,typename... Args>
struct function_traits<ReturnType(ClassType::*)(Args...) const>
{
    enum { arity = sizeof...(Args) };
    typedef function<ReturnType (Args...)> f_type;
};

// Handles pointers to member functions
template<typename ClassType,typename ReturnType,typename... Args>
struct function_traits<ReturnType(ClassType::*)(Args...)>
{
    enum { arity = sizeof...(Args) };
    typedef function<ReturnType (Args...)> f_type;
};

// Handles function pointers
template<typename ReturnType,typename... Args>
struct function_traits<ReturnType (*)(Args...)>
{
    enum { arity = sizeof...(Args) };
    typedef function<ReturnType (Args...)> f_type;
};

template<typename L>
static typename function_traits<L>::f_type MakeFunction(L l)
{ return (typename function_traits<L>::f_type)(l); }

// Handles std::bind & multiple function call operator()'s
template<typename ReturnType,typename... Args,class T>
auto MakeFunction(T&& t)
  -> std::function<decltype(ReturnType(t(std::declval<Args>()...)))(Args...)>
{ return {std::forward<T>(t)}; }

// For explicit overloads
template<typename ReturnType,typename... Args>
auto MakeFunction(ReturnType(*p)(Args...))
    -> std::function<ReturnType(Args...)>
{ return {p}; }

// For explicit overloads
template<typename ReturnType,typename... Args,typename ClassType>
auto MakeFunction(ReturnType(ClassType::*p)(Args...))
    -> std::function<ReturnType(Args...)>
{ return {p}; }

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

// Insists that the index can be found
Int Find( const vector<Int>& sortedInds, Int index );

#ifdef EL_HAVE_PRETTY_FUNCTION
# define EL_FUNCTION __PRETTY_FUNCTION__
#else
# define EL_FUNCTION __func__
#endif

#define EL_LOGIC_ERROR(...) \
 El::LogicError(EL_FUNCTION," in ",__FILE__,"@",__LINE__,": ",__VA_ARGS__);
#define EL_RUNTIME_ERROR(...) \
 El::RuntimeError(EL_FUNCTION," in ",__FILE__,"@",__LINE__,": ",__VA_ARGS__);
#define EL_DEBUG_CSE EL_DEBUG_ONLY(El::CSE cse(EL_FUNCTION))

} // namespace El

#endif // ifndef EL_ENVIRONMENT_DECL_HPP
