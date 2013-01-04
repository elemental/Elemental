/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {

// For initializing and finalizing Elemental
void Initialize( int& argc, char**& argv );
void Finalize();
bool Initialized();

// For getting the MPI argument instance (for internal usage)
class MpiArgs;
MpiArgs& GetArgs();

// For processing command-line arguments
template<typename T>
T Input( std::string name, std::string desc );
template<typename T>
T Input( std::string name, std::string desc, T defaultVal );
void ProcessInput();
void PrintInputReport();

// For getting and setting the algorithmic blocksize
int Blocksize();
void SetBlocksize( int blocksize );

// For manipulating the algorithmic blocksize as a stack
void PushBlocksizeStack( int blocksize );
void PopBlocksizeStack();

// Replacement for std::memcpy, which is known to often be suboptimal.
// Notice the sizeof(T) is no longer required.
template<typename T>
void MemCopy( T* dest, const T* source, std::size_t numEntries );

// Generalization of std::memcpy so that unit strides are not required
template<typename T>
void StridedMemCopy
(       T* dest,   std::size_t destStride,
  const T* source, std::size_t sourceStride, std::size_t numEntries );

// Replacement for std::memset, which is likely suboptimal and hard to extend
// to non-POD datatypes. Notice that sizeof(T) is no longer required.
template<typename T>
void MemZero( T* buffer, std::size_t numEntries );

// Euclidean (l_2) magnitudes
template<typename R>
R Abs( const R& alpha );
template<typename R>
R Abs( const Complex<R>& alpha );

// Square-root free (l_1) magnitudes
template<typename R>
R FastAbs( const R& alpha );
template<typename R>
R FastAbs( const Complex<R>& alpha );

// Return the real part of a real or complex number
template<typename R>
R RealPart( const R& alpha );
template<typename R>
R RealPart( const Complex<R>& alpha );

// Return the imaginary part of a real or complex number
template<typename R>
R ImagPart( const R& alpha );
template<typename R>
R ImagPart( const Complex<R>& alpha );

// Conjugation
template<typename R>
R Conj( const R& alpha );
template<typename R>
Complex<R> Conj( const Complex<R>& alpha );

// Square root
template<typename R>
R Sqrt( const R& alpha );
template<typename R>
Complex<R> Sqrt( const Complex<R>& alpha );

// Cosine
template<typename R>
R Cos( const R& alpha );
template<typename R>
Complex<R> Cos( const Complex<R>& alpha );

// Sine
template<typename R>
R Sin( const R& alpha );
template<typename R>
Complex<R> Sin( const Complex<R>& alpha );

// Tangent
template<typename R>
R Tan( const R& alpha );
template<typename R>
Complex<R> Tan( const Complex<R>& alpha );

// Hyperbolic cosine
template<typename R>
R Cosh( const R& alpha );
template<typename R>
Complex<R> Cosh( const Complex<R>& alpha );

// Hyperbolic sine
template<typename R>
R Sinh( const R& alpha );
template<typename R>
Complex<R> Sinh( const Complex<R>& alpha );

// Inverse cosine
template<typename R>
R Acos( const R& alpha );
// TODO
/*
template<typename R>
Complex<R> Acos( const Complex<R>& alpha );
*/

// Inverse sine
template<typename R>
R Asin( const R& alpha );
// TODO
/*
template<typename R>
Complex<R> Asin( const Complex<R>& alpha );
*/

// Inverse tangent
template<typename R>
R Atan( const R& alpha );
// TODO
/*
template<typename R>
Complex<R> Atan( const Complex<R>& alpha );
*/

// Coordinate-based inverse tangent
template<typename R>
R Atan2( const R& y, const R& x );

// Inverse hyperbolic cosine
template<typename R>
R Acosh( const R& alpha );
// TODO
/*
template<typename R>
Complex<R> Acosh( const Complex<R>& alpha );
*/

// Inverse hyperbolic sine
template<typename R>
R Asinh( const R& alpha );
// TODO
/*
template<typename R>
Complex<R> Asinh( const Complex<R>& alpha );
*/

// Inverse hyperbolic tangent
template<typename R>
R Atanh( const R& alpha );
// TODO
/*
template<typename R>
Complex<R> Atanh( const Complex<R>& alpha );
*/

// Complex argument
template<typename R>
R Arg( const R& alpha );
template<typename R>
R Arg( const Complex<R>& alpha );

// Convert polar coordinates to the complex number
template<typename R>
Complex<R> Polar( const R& r, const R& theta=0 ); 

// Exponential
template<typename R>
R Exp( const R& alpha );
template<typename R>
Complex<R> Exp( const Complex<R>& alpha );

// Power, return alpha^beta
// TODO: Mixed versions, such as a real number to an integer power, 
//       or a complex number to a real power
template<typename R>
R Pow( const R& alpha, const R& beta );
template<typename R>
Complex<R> Pow( const Complex<R>& alpha, const Complex<R>& beta );

// Logarithm
template<typename R>
R Log( const R& alpha );
template<typename R>
Complex<R> Log( const Complex<R>& alpha );

// An exception which signifies that a matrix was unexpectedly singular.
class SingularMatrixException : public std::runtime_error 
{
public:
    SingularMatrixException( const char* msg="Matrix was singular" ) 
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

#ifndef RELEASE
void PushCallStack( std::string s );
void PopCallStack();
void DumpCallStack();
#endif // ifndef RELEASE

// We define an output stream that does nothing. This is done so that the 
// root process can be used to print data to a file's ostream while all other 
// processes use a null ostream. This is used within the DistMatrix class's
// 'Write' functions.
struct NullStream : std::ostream
{
    struct NullStreamBuffer : std::streambuf
    {
        int overflow( int c ) { return traits_type::not_eof(c); }
    } nullStreamBuffer_;

    NullStream() 
    : std::ios(&nullStreamBuffer_), std::ostream(&nullStreamBuffer_) 
    { }
};

} // namespace elem
