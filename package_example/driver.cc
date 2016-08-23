#include <El.hpp>

int main( int argc, char* argv[] )
{
    El::Environment env( argc, argv );
    El::Matrix<double> A;
    El::Identity( A, 8, 8 );
    El::Print( A, "A" );
    return 0;
}
