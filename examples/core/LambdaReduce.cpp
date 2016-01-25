#include "El.hpp"
using namespace El;

template<typename T,typename=EnableIf<IsReal<T>>>
void LambdaTest()
{
    if( mpi::Rank() == 0 )
        Output("Testing with ",TypeName<T>());

    mpi::Comm comm = mpi::COMM_WORLD;
    const int rank = mpi::Rank( comm );

    const T sumOfRanks =
      mpi::AllReduce
      ( T(rank), []( const T& x, const T& y ) { return x+y; },
        true, comm );
    if( rank == 0 )
        Output("sum of ranks=",sumOfRanks);

    const T minRank = 
      mpi::AllReduce
      ( T(rank), []( const T& x, const T& y ) { return Min(x,y); },
        true, comm );
    if( rank == 0 )
        Output("minimum rank=",minRank);

    const T maxRank = 
      mpi::AllReduce
      ( T(rank), []( const T& x, const T& y ) { return Max(x,y); },
        true, comm );
    if( rank == 0 )
        Output("maximum rank=",maxRank);

    const T prodRank = 
      mpi::AllReduce
      ( T(rank+1), []( const T& x, const T& y ) { return x*y; }, true, comm );
    if( rank == 0 )
        Output("product of ranks (plus one)=",prodRank);
}

template<typename T,typename=DisableIf<IsReal<T>>,typename=void>
void LambdaTest()
{
    if( mpi::Rank() == 0 )
        Output("Testing with ",TypeName<T>());

    mpi::Comm comm = mpi::COMM_WORLD;
    const int rank = mpi::Rank( comm );

    const T sumOfRanks =
      mpi::AllReduce
      ( T(rank), []( const T& x, const T& y ) { return x+y; },
        true, comm );
    if( rank == 0 )
        Output("sum of ranks=",sumOfRanks);

    const T prodRank = 
      mpi::AllReduce
      ( T(rank+1), []( const T& x, const T& y ) { return x*y; }, true, comm );
    if( rank == 0 )
        Output("product of ranks (plus one)=",prodRank);
}

int main( int argc, char* argv[] )
{
    Environment env( argc, argv );

    try
    {
        LambdaTest<Int>();

        LambdaTest<float>();
        LambdaTest<Complex<float>>();

        LambdaTest<double>();
        LambdaTest<Complex<double>>();

#ifdef EL_HAVE_QD
        LambdaTest<DoubleDouble>();
        LambdaTest<QuadDouble>();
#endif

#ifdef EL_HAVE_QUAD
        LambdaTest<Quad>();
        LambdaTest<Complex<Quad>>();
#endif

#ifdef EL_HAVE_MPC
        LambdaTest<BigInt>();
        LambdaTest<BigFloat>();
#endif
    }
    catch( std::exception& e ) { ReportException(e); }

    return 0;
}
