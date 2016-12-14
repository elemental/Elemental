/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

int
main( int argc, char* argv[] )
{
    El::Environment env( argc, argv );

    try
    {
        const El::Int n = El::Input("--size","width of matrix",100);
        const bool display = El::Input("--display","display matrices?",false);
        const bool print = El::Input("--print","print matrices?",false);
        const bool smallestFirst =
          El::Input("--smallestFirst","smallest norm first?",false);
        El::ProcessInput();
        El::PrintInputReport();
        const El::Int m = n;

        El::QRCtrl<double> ctrl;
        ctrl.smallestFirst = smallestFirst;

        El::DistMatrix<El::Complex<double>> A;
        El::GKS( A, n );
        const double frobA = El::FrobeniusNorm( A );
        if( display )
            El::Display( A, "A" );
        if( print )
            El::Print( A, "A" );

        // Compute the pivoted QR decomposition of A, but do not overwrite A
        auto QRPiv( A );
        El::DistMatrix<El::Complex<double>> householderScalarsPiv;
        El::DistMatrix<double> signaturePiv;
        El::DistPermutation Omega;
        El::Timer PQRtimer;
        if( El::mpi::Rank() == 0 )
            PQRtimer.Start();
        El::QR( QRPiv, householderScalarsPiv, signaturePiv, Omega );
        if( El::mpi::Rank() == 0 )
            PQRtimer.Stop();
        if( display )
        {
            El::Display( QRPiv, "QRPiv" );
            El::Display( householderScalarsPiv, "householderScalarsPiv" );
            El::Display( signaturePiv, "signaturePiv" );

            El::DistMatrix<El::Int> OmegaFull;
            Omega.ExplicitMatrix( OmegaFull );
            El::Display( OmegaFull, "Omega" );
        }
        if( print )
        {
            El::Print( QRPiv, "QRPiv" );
            El::Print( householderScalarsPiv, "householderScalarsPiv" );
            El::Print( signaturePiv, "signaturePiv" );

            El::DistMatrix<El::Int> OmegaFull;
            Omega.ExplicitMatrix( OmegaFull );
            El::Print( OmegaFull, "Omega" );
        }

        // Compute the standard QR decomposition of A
        auto QRNoPiv( A );
        El::DistMatrix<El::Complex<double>> householderScalarsNoPiv;
        El::DistMatrix<double> signatureNoPiv;
        El::Timer QRtimer;
        if( El::mpi::Rank() == 0 )
            QRtimer.Start();
        El::QR( QRNoPiv, householderScalarsNoPiv, signatureNoPiv );
        if( El::mpi::Rank() == 0 )
            QRtimer.Start();
        if( display )
        {
            El::Display( QRNoPiv, "QRNoPiv" );
            El::Display( householderScalarsNoPiv, "householderScalarsNoPiv" );
            El::Display( signatureNoPiv, "signatureNoPiv" );
        }
        if( print )
        {
            El::Print( QRNoPiv, "QRNoPiv" );
            El::Print( householderScalarsNoPiv, "householderScalarsNoPiv" );
            El::Print( signatureNoPiv, "signatureNoPiv" );
        }

        // Check the error in the pivoted QR factorization,
        // || A Omega^T - Q R ||_F / || A ||_F
        auto E( QRPiv );
        El::MakeTrapezoidal( El::UPPER, E );
        El::qr::ApplyQ
        ( El::LEFT, El::NORMAL, QRPiv, householderScalarsPiv, signaturePiv, E );
        Omega.InversePermuteCols( E );
        E -= A;
        const double frobQRPiv = El::FrobeniusNorm( E );
        if( display )
            El::Display( E, "A P - Q R" );
        if( print )
            El::Print( E, "A P - Q R" );

        // Check the error in the standard QR factorization,
        // || A - Q R ||_F / || A ||_F
        E = QRNoPiv;
        El::MakeTrapezoidal( El::UPPER, E );
        El::qr::ApplyQ
        ( El::LEFT, El::NORMAL,
          QRNoPiv, householderScalarsNoPiv, signatureNoPiv, E );
        E -= A;
        const double frobQRNoPiv = El::FrobeniusNorm( E );
        if( display )
            El::Display( E, "A - Q R" );
        if( print )
            El::Print( E, "A - Q R" );

        // Check orthogonality of pivoted Q, || I - Q^H Q ||_F / || A ||_F
        El::Identity( E, m, n );
        El::qr::ApplyQ
        ( El::LEFT, El::NORMAL,
          QRPiv, householderScalarsPiv, signaturePiv, E );
        El::qr::ApplyQ
        ( El::LEFT, El::ADJOINT,
          QRPiv, householderScalarsPiv, signaturePiv, E );
        const El::Int k = El::Min(m,n);
        auto EUpper = E( El::IR(0,k), El::IR(0,k) );
        El::ShiftDiagonal( EUpper, El::Complex<double>(-1) );
        const double frobOrthogPiv = El::FrobeniusNorm( EUpper );
        if( display )
            El::Display( E, "pivoted I - Q^H Q" );
        if( print )
            El::Print( E, "pivoted I - Q^H Q" );

        // Check orthogonality of unpivoted Q, || I - Q^H Q ||_F / || A ||_F
        El::Identity( E, m, n );
        El::qr::ApplyQ
        ( El::LEFT, El::NORMAL,
          QRPiv, householderScalarsPiv, signaturePiv, E );
        El::qr::ApplyQ
        ( El::LEFT, El::ADJOINT,
          QRPiv, householderScalarsPiv, signaturePiv, E );
        EUpper = E( El::IR(0,k), El::IR(0,k) );
        El::ShiftDiagonal( EUpper, El::Complex<double>(-1) );
        const double frobOrthogNoPiv = El::FrobeniusNorm( EUpper );
        if( display )
            El::Display( E, "unpivoted I - Q^H Q" );
        if( print )
            El::Print( E, "unpivoted I - Q^H Q" );

        if( El::mpi::Rank() == 0 )
        {
            El::Output("Pivot QR time: ",PQRtimer.Total()," secs\n",
                       "      QR time: ",QRtimer.Total()," secs");
            El::Output
            ("|| A ||_F = ",frobA,"\n\n",
             "With pivoting: \n",
             "    || A P - Q R ||_F / || A ||_F = ",frobQRPiv/frobA,"\n",
             "    || I - Q^H Q ||_F / || A ||_F = ",frobOrthogPiv/frobA,"\n\n",
             "Without pivoting: \n",
             "    || A - Q R ||_F / || A ||_F = ",frobQRNoPiv/frobA,"\n",
             "    || I - Q^H Q ||_F / || A ||_F = ",frobOrthogNoPiv/frobA,"\n");
        }
    }
    catch( std::exception& e ) { El::ReportException(e); }

    return 0;
}
