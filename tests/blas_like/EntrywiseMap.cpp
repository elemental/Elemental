/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
using namespace El;

// Convert a random variable drawn uniformly from the unit ball to a
// power-law random variable
template<typename T>
T PowerLawVariable( T x, Base<T> alpha )
{
    typedef Base<T> Real;

    // Generate uniform random variable in [0,0.999]
    Real u;
    if( sizeof(Real) == sizeof(T) )
    {
        u = RealPart(T(0.5)*(x+T(1)));
    }
    else
    {
        u = ImagPart( Log(x) ) / (2*M_PI);
        if( u < Real(0) )
            u += Real(1);
    }
    u = Max(Min(u, Real(1)), Real(0));
      
    // Generate random number from power-law distribution
    int iters = int(Pow((Real(1)-u), -Real(1)/(alpha -Real(1))));

    // Increment variable until it is equal to power-law variable
    T y = T(0);
    for(int i=0; i<iters; ++i)
      y += T(1);

    return y;
}

template<typename T>
void TestEntrywiseMap
( Int m,
  Int n,
  std::function<T(const T&)> func,
  Int numThreads,
  const Grid& g,
  bool print )
{
    OutputFromRoot(g.Comm(),"Testing with ",TypeName<T>());
    PushIndent();

    double runTime, opsPerSec;
    Timer timer;

#ifdef EL_HYBRID
    if( numThreads > 0 )
    {
        omp_set_num_threads(numThreads);
    }
#endif

    // Generate random matrix
    DistMatrix<T> A(g);
    Uniform( A, m, n );
    if( print )
    {
        Print( A, "A" );
    }

    // Apply entrywise map
    if( g.Rank() == 0 )
    {
        Output("  Starting EntrywiseMap");
    }
    timer.Start();
    EntrywiseMap( A, func );
    mpi::Barrier( g.Comm() );
    runTime = timer.Stop();

    // Print results
    opsPerSec = double(m)*double(n) / runTime;
    OutputFromRoot(g.Comm(),"Finished in ",runTime," seconds (",opsPerSec," ops/s)");
    if( print )
    {
        Print( A, "func(A)" );
    }

    PopIndent();
}

int
main( int argc, char* argv[] )
{
    Environment env( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;

    try
    {
        // Get command-line arguments
        int gridHeight = Input("--gridHeight","height of process grid",0);
        const Int m = Input("--m","height of matrix",100);
        const Int n = Input("--n","width of matrix",120);
        const Int func = Input("--func","function (0=exp,1=reLU,2=unbalanced function)",0);
        const Int numThreads = Input("--numThreads","number of OpenMP threads (-1=default)",-1);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const bool colMajor = Input("--colMajor","column-major ordering?",true);
        const bool print = Input("--print","print matrices?",false);
        ProcessInput();
        PrintInputReport();

        // Set Elemental parameters
        if( gridHeight == 0 )
            gridHeight = Grid::DefaultHeight( mpi::Size(comm) );
        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( comm, gridHeight, order );
        SetBlocksize( nb );
        ComplainIfDebug();

        // Message
        OutputFromRoot(comm,"Testing EntrywiseMap");

        // Choose function
        std::function<float(const float&)> funcFloat;
        std::function<Complex<float>(const Complex<float>&)> funcComplexFloat;
        std::function<double(const double&)> funcDouble;
        std::function<Complex<double>(const Complex<double>&)> funcComplexDouble;
        switch(func)
        {
        case 0:
            // Exponential function
            funcFloat = Exp<float>;
            funcComplexFloat = Exp<Complex<float>>;
            funcDouble = Exp<double>;
            funcComplexDouble = Exp<Complex<double>>;
            break;
        case 1:
            // Rectified linear function
            funcFloat = [](float x) { return Max(x, 0.f); };
            funcComplexFloat = [](Complex<float> x) { return RealPart(x) > 0.f ? x : 0.f; };
            funcDouble = [](double x) { return Max(x, 0.); };
            funcComplexDouble = [](Complex<double> x) { return RealPart(x) > 0. ? x : 0.; };
            break;
        case 2:
            // Unbalanced function (number of iterations is power-law)
            funcFloat
              = ([](const float &x) { return PowerLawVariable<float>(x, 2.f); } );
            funcComplexFloat
              = ([](const Complex<float> &x)
                 { return PowerLawVariable<Complex<float>>(x, 2.f); });
            funcDouble
              = ( [](const double &x) { return PowerLawVariable<double>(x, 2.); } );
            funcComplexDouble
              = ([](const Complex<double> &x)
                 { return PowerLawVariable<Complex<double>>(x, 2.); });
            break;
        default:
            std::cerr << "Invalid function\n";
            return EXIT_FAILURE;
        }

        // Run tests
        TestEntrywiseMap<float>( m, n, funcFloat, numThreads, g, print );
        TestEntrywiseMap<Complex<float>>( m, n, funcComplexFloat, numThreads, g, print );
        TestEntrywiseMap<double>( m, n, funcDouble, numThreads, g, print );
        TestEntrywiseMap<Complex<double>>( m, n, funcComplexDouble, numThreads, g, print );
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
