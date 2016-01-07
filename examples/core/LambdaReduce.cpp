#include "El.hpp"
using namespace El;

int main( int argc, char* argv[] )
{
    Environment env( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const Int rank = mpi::Rank(comm);

    try
    {
        const Int sumOfRanks =
          mpi::AllReduce
          ( rank,
            []( const Int& x, const Int& y ) { return x+y; },
            true, comm );
        if( rank == 0 )
            Output("sum of ranks=",sumOfRanks);

        const Int minRank = 
          mpi::AllReduce
          ( rank,
            []( const Int& x, const Int& y ) { return Min(x,y); },
            true, comm );
        if( rank == 0 )
            Output("minimum rank=",minRank);

        const Int maxRank = 
          mpi::AllReduce
          ( rank,
            []( const Int& x, const Int& y ) { return Max(x,y); },
            true, comm );
        if( rank == 0 )
            Output("maximum rank=",maxRank);

        const Int prodRank = 
          mpi::AllReduce
          ( rank+1,
            []( const Int& x, const Int& y ) { return x*y; },
            true, comm );
        if( rank == 0 )
            Output("product of ranks (plus one)=",prodRank);
    }
    catch( std::exception& e ) { ReportException(e); }

    return 0;
}
