/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

template<typename T,
         typename=El::EnableIf<El::IsReal<T>>>
void LambdaTest()
{
    if( El::mpi::Rank() == 0 )
        El::Output("Testing with ",El::TypeName<T>());

    El::mpi::Comm comm = El::mpi::COMM_WORLD;
    const int rank = El::mpi::Rank( comm );

    const T sumOfRanks =
      El::mpi::AllReduce
      ( T(rank), []( const T& x, const T& y ) { return x+y; },
        true, comm );
    if( rank == 0 )
        El::Output("sum of ranks=",sumOfRanks);

    const T minRank =
      El::mpi::AllReduce
      ( T(rank), []( const T& x, const T& y ) { return El::Min(x,y); },
        true, comm );
    if( rank == 0 )
        El::Output("minimum rank=",minRank);

    const T maxRank =
      El::mpi::AllReduce
      ( T(rank), []( const T& x, const T& y ) { return El::Max(x,y); },
        true, comm );
    if( rank == 0 )
        El::Output("maximum rank=",maxRank);

    const T prodRank =
      El::mpi::AllReduce
      ( T(rank+1), []( const T& x, const T& y ) { return x*y; }, true, comm );
    if( rank == 0 )
        El::Output("product of ranks (plus one)=",prodRank);
}

template<typename T,
         typename=El::DisableIf<El::IsReal<T>>,
         typename=void>
void LambdaTest()
{
    if( El::mpi::Rank() == 0 )
        El::Output("Testing with ",El::TypeName<T>());

    El::mpi::Comm comm = El::mpi::COMM_WORLD;
    const int rank = El::mpi::Rank( comm );

    const T sumOfRanks =
      El::mpi::AllReduce
      ( T(rank), []( const T& x, const T& y ) { return x+y; },
        true, comm );
    if( rank == 0 )
        El::Output("sum of ranks=",sumOfRanks);

    const T prodRank =
      El::mpi::AllReduce
      ( T(rank+1), []( const T& x, const T& y ) { return x*y; }, true, comm );
    if( rank == 0 )
        El::Output("product of ranks (plus one)=",prodRank);
}

int main( int argc, char* argv[] )
{
    El::Environment env( argc, argv );

    try
    {
        LambdaTest<El::Int>();

        LambdaTest<float>();
        LambdaTest<El::Complex<float>>();

        LambdaTest<double>();
        LambdaTest<El::Complex<double>>();

#ifdef EL_HAVE_QD
        LambdaTest<El::DoubleDouble>();
        LambdaTest<El::Complex<El::DoubleDouble>>();
        LambdaTest<El::QuadDouble>();
        LambdaTest<El::Complex<El::QuadDouble>>();
#endif

#ifdef EL_HAVE_QUAD
        LambdaTest<El::Quad>();
        LambdaTest<El::Complex<El::Quad>>();
#endif

#ifdef EL_HAVE_MPC
        LambdaTest<El::BigInt>();
        LambdaTest<El::BigFloat>();
        LambdaTest<El::Complex<El::BigFloat>>();
#endif
    }
    catch( std::exception& e ) { El::ReportException(e); }

    return 0;
}
