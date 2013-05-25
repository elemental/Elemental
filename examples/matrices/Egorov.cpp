/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
// NOTE: It is possible to simply include "elemental.hpp" instead
#include "elemental-lite.hpp"
#include "elemental/matrices/Egorov.hpp"
#include "elemental/graphics.hpp"
using namespace elem;

template<typename R>
class FourierPhase {
public:
    FourierPhase( int n ) : n_(n), pi_(4*Atan(R(1))) { }
    R operator()( int i, int j ) const { return (-2*pi_*i*j)/n_; }
private:
    int n_;
    R pi_;
};

template<typename R>
class Phase {
public:
    Phase( int n ) : n_(n), pi_(4*Atan(R(1))) { }
    R operator()( int i, int j ) const 
    { return (-2*pi_*i*j)/n_ + Sqrt(R(i)*R(i) + R(j)*R(j)); }
private:
    int n_;
    R pi_;
};

int 
main( int argc, char* argv[] )
{
    Initialize( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::CommRank( comm );

    try
    {
        const int n = Input("--size","size of matrix",10);
        const bool print = Input("--print","print matrix?",true);
#ifdef HAVE_QT5
        const bool display = Input("--display","display matrix?",true);
#endif
        ProcessInput();
        PrintInputReport();

        FourierPhase<double> fourier( n );
        Phase<double> phase( n );

        DistMatrix<Complex<double> > F, G;
        Egorov( F, fourier, n );
        Egorov( G, phase, n ); 

        if( print )
        {
            F.Print("Egorov with Fourier phase:");
            G.Print("Egorov with more general phase:");
        }
#ifdef HAVE_QT5
        if( display )
        {
            Display( F, "Egorov with Fourier phase" );
            Display( G, "Egorov with more general phase" );
        }
#endif
    }
    catch( ArgException& e )
    {
        // There is nothing to do
    }
    catch( std::exception& e )
    {
        std::ostringstream os;
        os << "Process " << commRank << " caught error message:\n"
           << e.what() << std::endl;
        std::cerr << os.str();
#ifndef RELEASE
        DumpCallStack();
#endif
    }

    Finalize();
    return 0;
}
