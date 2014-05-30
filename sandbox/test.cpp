#include "El.hpp"
using namespace El;

int main( int argc, char* argv[] )
{
    Initialize( argc, argv );

    const Int m = Input("--m","matrix height",100);
    const Int n = Input("--n","matrix width",100);

    DistMatrix<double> A;
    Uniform( A, m, n );
    Display( A, "A" );

    Finalize();
    return 0;
}
