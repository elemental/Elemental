#include "El.hpp"
using namespace El;

int main( int argc, char* argv[] )
{
    Initialize( argc, argv );

    const Int m = Input("--m","matrix height",100);
    const Int n = Input("--n","matrix width",20);

    DistMatrix<double> A;
    Uniform( A, m, n );

    DistMatrix<Int,VR,STAR> p;
    DistMatrix<double,STAR,VR> Z;
    QRCtrl<double> ctrl;
    ctrl.boundRank = true;
    ctrl.maxRank = 20;
    ID( A, p, Z, ctrl );
    if( mpi::WorldRank() == 0 )
        std::cout << Z.Height() << ", " << Z.Width() << std::endl;

    Finalize();
    return 0;
}
