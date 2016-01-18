/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
using namespace std;
using namespace El;

template<typename F>
void TestCorrectness
( bool print,
  Pencil pencil, UpperOrLower uplo,
  const ElementalMatrix<F>& AOrig,
  const ElementalMatrix<F>& BOrig,
  const ElementalMatrix<F>& A,
  const ElementalMatrix<F>& B,
  const ElementalMatrix<Base<F>>& w,
  const ElementalMatrix<F>& X )
{
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    const Int n = X.Height();
    const Int k = X.Width();

    DistMatrix<Real,MR,STAR> w_MR_STAR(true,X.RowAlign(),g); 
    w_MR_STAR = w;

    if( pencil == AXBX )
    {
        // Set Y := BXW, where W is the diagonal eigenvalue matrix
        DistMatrix<F> Y( g );
        Y.AlignWith( X );
        Zeros( Y, n, k );
        Hemm( LEFT, uplo, F(1), BOrig, X, F(0), Y );
        DiagonalScale( RIGHT, NORMAL, w_MR_STAR, Y );
        // Y := Y - AX = BXW - AX
        Hemm( LEFT, uplo, F(-1), AOrig, X, F(1), Y );

        const Real frobNormA = HermitianFrobeniusNorm( uplo, AOrig );
        const Real frobNormB = HermitianFrobeniusNorm( uplo, BOrig );
        Real frobNormE = FrobeniusNorm( Y );
        if( g.Rank() == 0 )
            cout << "    ||A X - B X W||_F / max(||A||_F,||B||_F) = " 
                 << frobNormE/Max(frobNormA,frobNormB) << std::endl;
        DistMatrix<F> Z(g);
        Z = X;
        if( uplo == LOWER )
            Trmm( LEFT, LOWER, ADJOINT, NON_UNIT, F(1), B, Z );
        else
            Trmm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), B, Z );
        Identity( Y, k, k );
        Herk( uplo, ADJOINT, Real(-1), Z, Real(1), Y );
        frobNormE = FrobeniusNorm( Y );
        if( g.Rank() == 0 )
            cout << "    ||X^H B X - I||_F  / max(||A||_F,||B||_F) = " 
                 << frobNormE/Max(frobNormA,frobNormB) << endl;
    }
    else if( pencil == ABX )
    {
        // Set Y := BX
        DistMatrix<F> Y( g );
        Y.AlignWith( X );
        Zeros( Y, n, k );
        Hemm( LEFT, uplo, F(1), BOrig, X, F(0), Y );
        // Set Z := AY = ABX
        DistMatrix<F> Z( n, k, g );
        Hemm( LEFT, uplo, F(1), AOrig, Y, F(0), Z );
        // Set Z := Z - XW = ABX - XW
        for( Int jLoc=0; jLoc<Z.LocalWidth(); ++jLoc )
        {
            const Real omega = w_MR_STAR.GetLocal(jLoc,0); 
            for( Int iLoc=0; iLoc<Z.LocalHeight(); ++iLoc )
            {
                const F chi = X.GetLocal(iLoc,jLoc);
                const F zeta = Z.GetLocal(iLoc,jLoc);
                Z.SetLocal(iLoc,jLoc,zeta-omega*chi);
            }
        }
        // Find the infinity norms of A, B, X, and ABX-XW
        const Real frobNormA = HermitianFrobeniusNorm( uplo, AOrig );
        const Real frobNormB = HermitianFrobeniusNorm( uplo, BOrig );
        Real frobNormE = FrobeniusNorm( Z );
        if( g.Rank() == 0 )
            cout << "    ||A B X - X W||_F  / max(||A||_F,||B||_F) = " 
                 << frobNormE/Max(frobNormA,frobNormB) << std::endl;
        Z = X;
        if( uplo == LOWER )
            Trmm( LEFT, LOWER, ADJOINT, NON_UNIT, F(1), B, Z );
        else
            Trmm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), B, Z );
        Identity( Y, k, k );
        Herk( uplo, ADJOINT, Real(-1), Z, Real(1), Y );
        frobNormE = FrobeniusNorm( Y );
        if( g.Rank() == 0 )
            cout << "    ||X^H B X - I||_F / Max(||A||_F,||B||_F) = " 
                 << frobNormE/Max(frobNormA,frobNormB) << endl;
    }
    else /* pencil == BAX */
    {
        // Set Y := AX
        DistMatrix<F> Y( g );
        Y.AlignWith( X );
        Zeros( Y, n, k );
        Hemm( LEFT, uplo, F(1), AOrig, X, F(0), Y );
        // Set Z := BY = BAX
        DistMatrix<F> Z( n, k, g );
        Hemm( LEFT, uplo, F(1), BOrig, Y, F(0), Z );
        // Set Z := Z - XW = BAX-XW
        for( Int jLoc=0; jLoc<Z.LocalWidth(); ++jLoc )
        {
            const Real omega = w_MR_STAR.GetLocal(jLoc,0); 
            for( Int iLoc=0; iLoc<Z.LocalHeight(); ++iLoc )
            {
                const F chi = X.GetLocal(iLoc,jLoc);
                const F zeta = Z.GetLocal(iLoc,jLoc);
                Z.SetLocal(iLoc,jLoc,zeta-omega*chi);
            }
        }
        // Find the infinity norms of A, B, X, and BAX-XW
        const Real frobNormA = HermitianFrobeniusNorm( uplo, AOrig );
        const Real frobNormB = HermitianFrobeniusNorm( uplo, BOrig );
        Real frobNormE = FrobeniusNorm( Z );
        if( g.Rank() == 0 )
            cout << "    ||B A X - X W||_F / Max(||A||_F,||B||_F) = " 
                 << frobNormE/Max(frobNormA,frobNormB) << std::endl;
        Z = X;
        if( uplo == LOWER )
            Trsm( LEFT, LOWER, NORMAL, NON_UNIT, F(1), B, Z );
        else
            Trsm( LEFT, UPPER, ADJOINT, NON_UNIT, F(1), B, Z );
        Identity( Y, k, k );
        Herk( uplo, ADJOINT, Real(-1), Z, Real(1), Y );
        frobNormE = FrobeniusNorm( Y );
        if( g.Rank() == 0 )
            cout << "    ||X^H B^-1 X - I||_F  / Max(||A||_F,||B||_F) = " 
                 << frobNormE/Max(frobNormA,frobNormB) << endl;
    }
}

template<typename F,Dist U=MC,Dist V=MR,Dist S=VR>
void TestHermitianGenDefEig
( bool testCorrectness,
  bool print,
  Pencil pencil,
  bool onlyEigvals,
  UpperOrLower uplo, 
  Int m,
  SortType sort,
  const Grid& g, 
  const HermitianEigSubset<Base<F>> subset, 
  const HermitianEigCtrl<F>& ctrl )
{
    typedef Base<F> Real;
    DistMatrix<F,U,V> A(g), B(g), AOrig(g), BOrig(g);
    DistMatrix<Real,S,STAR> w(g);
    DistMatrix<F> X(g);
    if( g.Rank() == 0 )
        Output("Testing with ",TypeName<F>());

    HermitianUniformSpectrum( A, m, 1, 10 );
    if( pencil == BAX )
    {
        // Because we will multiply by L three times, generate HPD B more 
        // carefully than just adding m to its diagonal entries.
        Zeros( B, m, m );
        DistMatrix<F> C(g);
        Uniform( C, m, m );
        Herk( uplo, ADJOINT, Real(1), C, Real(0), B );
    }
    else
        HermitianUniformSpectrum( B, m, 1, 10 );

    if( testCorrectness && !onlyEigvals )
    {
        AOrig = A;
        BOrig = B;
    }
    if( print )
    {
        Print( A, "A" );
        Print( B, "B" );
    }

    if( g.Rank() == 0 )
        cout << "  Starting Hermitian Generalized-Definite Eigensolver..."
             << std::endl;
    mpi::Barrier( g.Comm() );
    const double startTime = mpi::Time();
    if( onlyEigvals )
        HermitianGenDefEig( pencil, uplo, A, B, w, sort, subset, ctrl );
    else
        HermitianGenDefEig( pencil, uplo, A, B, w, X, sort, subset, ctrl );
    mpi::Barrier( g.Comm() );
    const double runTime = mpi::Time() - startTime;
    if( g.Rank() == 0 )
        cout << "  Time = " << runTime << " seconds." << endl;
    if( print )
    {
        Print( w, "eigenvalues:" );
        if( !onlyEigvals )
            Print( X, "eigenvectors:" );
    }
    if( testCorrectness && !onlyEigvals )
        TestCorrectness( print, pencil, uplo, AOrig, BOrig, A, B, w, X );
}

int 
main( int argc, char* argv[] )
{
    Environment env( argc, argv );
    mpi::Comm comm = mpi::COMM_WORLD;
    const int commRank = mpi::Rank( comm );
    const int commSize = mpi::Size( comm );

    try
    {
        Int r = Input("--gridHeight","height of process grid",0);
        const bool colMajor = Input("--colMajor","column-major ordering?",true);
        const Int pencilInt = Input("--pencil",
             "1 is A x = lambda B x, "
             "2 is A B x = lambda x, "
             "3 is B A x = lambda x",1);
        const bool onlyEigvals = Input 
            ("--onlyEigvals","only compute eigenvalues?",false);
        const char range = Input
            ("--range",
             "range of eigenpairs: 'A' for all, 'I' for index range, "
             "'V' for value range",'A');
        const Int il = Input("--il","lower bound of index range",0);
        const Int iu = Input("--iu","upper bound of index range",100);
        const double vl = Input("--vl","lower bound of value range",0.);
        const double vu = Input("--vu","upper bound of value range",100.);
        const Int sortInt = Input("--sort","sort type",0);
        const char uploChar = Input("--uplo","upper or lower storage: L/U",'L');
        const Int m = Input("--height","height of matrix",100);
        const Int nb = Input("--nb","algorithmic blocksize",96);
        const Int nbLocal = Input("--nbLocal","local blocksize",32);
        const bool avoidTrmv = 
            Input("--avoidTrmv","avoid Trmv based Symv",true);
        const bool testCorrectness = Input
            ("--correctness","test correctness?",true);
        const bool print = Input("--print","print matrices?",false);
        const bool testReal = Input("--testReal","test real matrices?",true);
        const bool testCpx = Input("--testCpx","test complex matrices?",true);
        const bool timeStages = Input("--timeStages","time stages?",true);
        ProcessInput();
        PrintInputReport();

        if( r == 0 )
            r = Grid::FindFactor( commSize );
        const GridOrder order = ( colMajor ? COLUMN_MAJOR : ROW_MAJOR );
        const Grid g( comm, r, order );
        const UpperOrLower uplo = CharToUpperOrLower( uploChar );
        SetBlocksize( nb );
        if( range != 'A' && range != 'I' && range != 'V' )
            throw runtime_error("'range' must be 'A', 'I', or 'V'");
        const SortType sort = static_cast<SortType>(sortInt);
        if( testCorrectness && onlyEigvals && commRank==0 )
            cout << "Cannot test correctness with only eigenvalues." << endl;
        Pencil pencil;
        string pencilString;
        if( pencilInt == 1 )
        {
            pencil = AXBX;
            pencilString = "AXBX";
        }
        else if( pencilInt == 2 )
        {
            pencil = ABX;
            pencilString = "ABX";
        }
        else if( pencilInt == 3 )
        {
            pencil = BAX;
            pencilString = "BAX";
        }
        else
            LogicError("Invalid pencil integer");
        ComplainIfDebug();
        if( commRank == 0 )
            cout << "Will test " 
                 << ( uplo==LOWER ? "lower" : "upper" )
                 << " " << pencilString << " HermitianGenDefEig." << endl;

        HermitianEigSubset<double> subset;
        if( range == 'I' )
        {
            subset.indexSubset = true;
            subset.lowerIndex = il;
            subset.upperIndex = iu;
        }
        else if( range == 'V' )
        {
            subset.rangeSubset = true;
            subset.lowerBound = vl;
            subset.upperBound = vu;
        }

        HermitianEigCtrl<double> ctrl_d;
        HermitianEigCtrl<Complex<double>> ctrl_z;
        ctrl_d.timeStages = timeStages;
        ctrl_z.timeStages = timeStages;
        ctrl_d.tridiagCtrl.symvCtrl.bsize = nbLocal;
        ctrl_z.tridiagCtrl.symvCtrl.bsize = nbLocal;
        ctrl_d.tridiagCtrl.symvCtrl.avoidTrmvBasedLocalSymv = avoidTrmv;
        ctrl_z.tridiagCtrl.symvCtrl.avoidTrmvBasedLocalSymv = avoidTrmv;

        if( commRank == 0 )
            cout << "Normal tridiag algorithms:" << endl;
        ctrl_d.tridiagCtrl.approach = HERMITIAN_TRIDIAG_NORMAL;
        ctrl_z.tridiagCtrl.approach = HERMITIAN_TRIDIAG_NORMAL;
        if( testReal )
            TestHermitianGenDefEig<double>
            ( testCorrectness, print, 
              pencil, onlyEigvals, uplo, m, sort, g, subset, ctrl_d );
        if( testCpx )
            TestHermitianGenDefEig<Complex<double>>
            ( testCorrectness, print, 
              pencil, onlyEigvals, uplo, m, sort, g, subset, ctrl_z );

        if( commRank == 0 )
            cout << "Square row-major algorithms:" << endl;
        ctrl_d.tridiagCtrl.approach = HERMITIAN_TRIDIAG_SQUARE;
        ctrl_z.tridiagCtrl.approach = HERMITIAN_TRIDIAG_SQUARE;
        ctrl_d.tridiagCtrl.order = ROW_MAJOR;
        ctrl_z.tridiagCtrl.order = ROW_MAJOR;
        if( testReal )
            TestHermitianGenDefEig<double>
            ( testCorrectness, print, 
              pencil, onlyEigvals, uplo, m, sort, g, subset, ctrl_d );
        if( testCpx )
            TestHermitianGenDefEig<Complex<double>>
            ( testCorrectness, print, 
              pencil, onlyEigvals, uplo, m, sort, g, subset, ctrl_z );

        if( commRank == 0 )
            cout << "Square column-major algorithms:" << endl;
        ctrl_d.tridiagCtrl.approach = HERMITIAN_TRIDIAG_SQUARE;
        ctrl_z.tridiagCtrl.approach = HERMITIAN_TRIDIAG_SQUARE;
        ctrl_d.tridiagCtrl.order = COLUMN_MAJOR;
        ctrl_z.tridiagCtrl.order = COLUMN_MAJOR;
        if( testReal )
            TestHermitianGenDefEig<double>
            ( testCorrectness, print, 
              pencil, onlyEigvals, uplo, m, sort, g, subset, ctrl_d );
        if( testCpx )
            TestHermitianGenDefEig<Complex<double>>
            ( testCorrectness, print, 
              pencil, onlyEigvals, uplo, m, sort, g, subset, ctrl_z );

        // Also test with non-standard distributions
        if( commRank == 0 )
            cout << "Nonstandard distributions:" << endl;
        if( testReal )
            TestHermitianGenDefEig<double,MR,MC,MC>
            ( testCorrectness, print, 
              pencil, onlyEigvals, uplo, m, sort, g, subset, ctrl_d );
        if( testCpx )
            TestHermitianGenDefEig<Complex<double>,MR,MC,MC>
            ( testCorrectness, print, 
              pencil, onlyEigvals, uplo, m, sort, g, subset, ctrl_z );
    }
    catch( exception& e ) { ReportException(e); }

    return 0;
}
