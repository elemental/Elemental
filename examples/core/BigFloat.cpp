#include <El.hpp>
using namespace El;

int main( int argc, char* argv[] )
{
    Environment env( argc, argv );

    try
    {
        const Int m = Input("--m","matrix height",100);
        const Int n = Input("--n","matrix width",100);
        ProcessInput();

        std::cout.precision(20);

#ifdef EL_HAVE_MPC
        Matrix<BigFloat> A;
        Uniform( A, m, n );
        Print( A, "A" );
#else
        Output("Elemental was not built with MPC support");
#endif
    }
    catch( std::exception& e ) { ReportException(e); }

    return 0;
}
