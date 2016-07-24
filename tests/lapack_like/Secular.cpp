/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
using namespace El;

extern "C" {

void EL_LAPACK(dlasd4)
( const BlasInt* n,
  const BlasInt* i,
  const double* d,
  const double* z,
  double* dMinusShift,
  const double* rho,
  double* sigma,
  double* dPlusShift,
  BlasInt* info );

} // extern "C"

namespace El {
namespace bidiag_svd {
namespace dc {

template<typename Real,typename=EnableIf<IsReal<Real>>>
struct SecularInfo
{
    Real singularValue;
    Int numIterations = 0;
    Int numAlternations = 0;
    Int numCubicFailures = 0;
};

template<typename Real>
struct SecularCtrl
{
    const Int maxIterations = 400; // This is LAPACK's limit
    const Real sufficientDecay = Real(1)/Real(10);
    FlipOrClip negativeFix;
    const bool progress = true;
};

template<typename Real,typename=EnableIf<IsReal<Real>>>
SecularInfo<Real>
Secular
( Int whichSingularValue,
  const Matrix<Real>& d,
  const Real& rho,
  const Matrix<Real>& u,
        Matrix<Real>& dMinusShift,
        Matrix<Real>& dPlusShift,
  const SecularCtrl<Real>& ctrl=SecularCtrl<Real>() );

} // namespace dc
} // namespace bidiag_svd
} // namespace El

int main( int argc, char* argv[] )
{
    Environment env( argc, argv );

    try
    {
        const Int n = Input("--n","matrix size",100);
        ProcessInput();

        Matrix<double> d, u;
        Uniform( d, n, 1, double(2.), double(2.) );
        Sort( d );
        Gaussian( u, n, 1 );
        u *= double(1) / FrobeniusNorm( u );

        const double rho = SampleUniform( double(1), double(0.5) );
         
        Print( d, "d" );
        Output( "rho=", rho );
        Print( u, "u" );

        Matrix<double> wSecular(n,1);
        Matrix<double> dMinusShift(n,1), dPlusShift(n,1);
        for( Int i=0; i<n; ++i )
        {
            BlasInt infoLAPACK;
            double sigmaLAPACK;
            const BlasInt ip1 = i+1;
            EL_LAPACK(dlasd4)
            ( &n, &ip1, d.Buffer(), u.Buffer(),
              dMinusShift.Buffer(), &rho, &sigmaLAPACK,
              dPlusShift.Buffer(), &infoLAPACK );
            Output("sigmaLAPACK=",sigmaLAPACK);
            Output("");
            Output("");

            auto info =
              bidiag_svd::dc::Secular( i, d, rho, u, dMinusShift, dPlusShift );
            wSecular(i) = info.singularValue*info.singularValue;

            Output
            ("i=", i, ", root=", info.singularValue, ", numIterations=",
             info.numIterations, ", numAlternations=", info.numAlternations, 
             ", numCubicFailures=", info.numCubicFailures);
            Output("");
            Output("");
        }

        Matrix<double> A, w;
        Matrix<double> dSquared;
        Hadamard( d, d, dSquared );
        Diagonal( A, dSquared );
        Syrk( LOWER, NORMAL, rho, u, double(1), A );
        auto hermEigInfo = HermitianEig( LOWER, A, w );
        Print( w, "w" );

        auto wDiff( w );
        wDiff -= wSecular;
        const double diffNorm = FrobeniusNorm( wDiff );
        Output("|| w - wSecular ||_F = ", diffNorm);
    }
    catch( std::exception& e ) { ReportException(e); }

    return 0;
}
