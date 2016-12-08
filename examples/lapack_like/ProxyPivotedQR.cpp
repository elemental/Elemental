/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

// NOTE: This definition is nearly identical to StandardProxy but is
//       meant to demonstrate how to manually build and use a pivot proxy
template<typename Field>
class Proxy
{
private:
    El::Int numPower_, numOversample_;

public:
    Proxy( El::Int numPower=1, El::Int numOversample=10 )
    : numPower_(numPower), numOversample_(numOversample)
    { }

    void operator()
    ( const El::Matrix<Field>& A,
            El::Permutation& Omega,
            El::Int numPivots,
            bool smallestFirst=false ) const
    {
        const El::Int m = A.Height();

        // Generate a Gaussian random matrix
        El::Matrix<Field> G;
        El::Gaussian( G, numPivots+numOversample_, m );

        // Form G (A A^H)^q A = G A (A^H A)^2
        El::Matrix<Field> Y, Z;
        El::Gemm( El::NORMAL, El::NORMAL, Field(1), G, A, Y );
        for( El::Int powerIter=0; powerIter<numPower_; ++powerIter )
        {
            El::Gemm( El::NORMAL, El::ADJOINT, Field(1), Y, A, Z );
            El::Gemm( El::NORMAL, El::NORMAL, Field(1), Z, A, Y );
        }

        El::QRCtrl<El::Base<Field>> ctrl;
        ctrl.boundRank = true;
        ctrl.maxRank = numPivots;
        ctrl.smallestFirst = smallestFirst;
        El::Matrix<Field> householderScalars;
        El::Matrix<El::Base<Field>> signature;
        El::QR( Y, householderScalars, signature, Omega, ctrl );
    }

    void operator()
    ( const El::AbstractDistMatrix<Field>& APre,
            El::DistPermutation& Omega,
            El::Int numPivots,
            bool smallestFirst=false ) const
    {
        const El::Int m = APre.Height();
        const El::Grid& grid = APre.Grid();

        El::DistMatrixReadProxy<Field,Field,El::MC,El::MR> AProxy( APre );
        auto& A = AProxy.GetLocked();

        // Generate a Gaussian random matrix
        El::DistMatrix<Field> G(grid);
        El::Gaussian( G, numPivots+numOversample_, m );

        // Form G (A A^H)^q A = G A (A^H A)^2
        El::DistMatrix<Field> Y(grid), Z(grid);
        El::Gemm( El::NORMAL, El::NORMAL, Field(1), G, A, Y );

        for( El::Int powerIter=0; powerIter<numPower_; ++powerIter )
        {
            El::Gemm( El::NORMAL, El::ADJOINT, Field(1), Y, A, Z );
            El::Gemm( El::NORMAL, El::NORMAL, Field(1), Z, A, Y );
        }

        El::QRCtrl<El::Base<Field>> ctrl;
        ctrl.boundRank = true;
        ctrl.maxRank = numPivots;
        ctrl.smallestFirst = smallestFirst;
        El::DistMatrix<Field,El::MD,El::STAR> householderScalars(grid);
        El::DistMatrix<El::Base<Field>,El::MD,El::STAR> signature(grid);
        El::QR( Y, householderScalars, signature, Omega, ctrl );
    }
};

int main( int argc, char* argv[] )
{
    El::Environment env( argc, argv );
    El::mpi::Comm comm = El::mpi::COMM_WORLD;
    const int commRank = El::mpi::Rank(comm);

    try
    {
        const El::Int n = El::Input("--n","matrix size",1000);
        const El::Int nb = El::Input("--nb","blocksize",64);
        //const double phi = El::Input("--phi","Kahan parameter",0.5);
        const bool panelPiv = El::Input("--panelPiv","panel pivoting?",false);
        const El::Int oversample =
          El::Input("--oversample","oversample factor",10);
        const El::Int numPower =
          El::Input("--numPower","# of power iterations",1);
        const bool smallestFirst =
          El::Input("--smallestFirst","smallest norms first?",false);
        const bool print = El::Input("--print","print matrices?",false);
        El::ProcessInput();
        El::PrintInputReport();

        El::SetBlocksize( nb );

        const El::Grid grid(comm);
        El::DistMatrix<double> A(grid);
        //El::Kahan( A, n, phi );
        El::Uniform( A, n, n );
        auto ACopy = A;
        if( print )
            El::Print( A, "A" );
        El::Timer timer;

        if( commRank == 0 )
            timer.Start();
        El::DistMatrix<double> householderScalars(grid), signature(grid);
        El::DistPermutation Omega(grid);
        Proxy<double> prox(numPower,oversample);
        El::qr::ProxyHouseholder
        ( A, householderScalars, signature, Omega, prox, panelPiv,
          smallestFirst );
        if( commRank == 0 )
            El::Output("Proxy QR time: ",timer.Stop()," seconds");
        if( print )
        {
            El::Print( A, "QR" );
            El::Print( householderScalars, "householderScalars" );
            El::Print( signature, "signature" );

            El::DistMatrix<El::Int> OmegaFull(grid);
            Omega.ExplicitMatrix( OmegaFull );
            El::Print( OmegaFull, "Omega" );
        }
        El::DistMatrix<double> diagR(grid);
        El::GetDiagonal( A, diagR );
        El::Print( diagR, "diag(R)" );

        A = ACopy;
        if( commRank == 0 )
            timer.Start();
        El::QRCtrl<double> ctrl;
        ctrl.smallestFirst = smallestFirst;
        El::QR( A, householderScalars, signature, Omega, ctrl );
        if( commRank == 0 )
            El::Output("Businger-Golub time: ",timer.Stop()," seconds");
        if( print )
        {
            El::Print( A, "QR" );
            El::Print( householderScalars, "householderScalars" );
            El::Print( signature, "signature" );

            El::DistMatrix<El::Int> OmegaFull;
            Omega.ExplicitMatrix( OmegaFull );
            El::Print( OmegaFull, "Omega" );
        }
        El::GetDiagonal( A, diagR );
        El::Print( diagR, "diag(R)" );

        A = ACopy;
        if( commRank == 0 )
            timer.Start();
        El::QR( A, householderScalars, signature );
        if( commRank == 0 )
            El::Output("Standard QR time: ",timer.Stop()," seconds");
        if( print )
        {
            El::Print( A, "QR" );
            El::Print( householderScalars, "householderScalars" );
            El::Print( signature, "signature" );
        }
        El::GetDiagonal( A, diagR );
        El::Print( diagR, "diag(R)" );
    }
    catch( std::exception& e ) { El::ReportException(e); }

    return 0;
}
