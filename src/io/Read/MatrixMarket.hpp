/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_READ_MATRIXMARKET_HPP
#define EL_READ_MATRIXMARKET_HPP

namespace El {
namespace read {

template<typename T>
void MatrixMarket( Matrix<T>& A, const string filename )
{
    DEBUG_CSE
    typedef Base<T> Real;
    std::ifstream file( filename.c_str() );
    if( !file.is_open() )
        RuntimeError("Could not open ",filename);

    // Read the header
    // ===============
    // Attempt to pull in the various header components
    // ------------------------------------------------
    string line, stamp, object, format, field, symmetry;
    if( !std::getline( file, line ) )
        RuntimeError("Could not extract header line");
    {
        std::stringstream lineStream( line );
        lineStream >> stamp; 
        if( stamp != string("%%MatrixMarket") )
            RuntimeError("Invalid Matrix Market stamp: ",stamp);
        if( !(lineStream >> object) ) 
            RuntimeError("Missing Matrix Market object");
        if( !(lineStream >> format) )
            RuntimeError("Missing Matrix Market format");
        if( !(lineStream >> field) )
            RuntimeError("Missing Matrix Market field");
        if( !(lineStream >> symmetry) )
            RuntimeError("Missing Matrix Market symmetry");
    }
    // Ensure that the header components are individually valid
    // --------------------------------------------------------
    const bool isMatrix = ( object == string("matrix") );
    const bool isArray = ( format == string("array") );
    const bool isComplex = ( field == string("complex") );
    const bool isPattern = ( field == string("pattern") );
    const bool isGeneral = ( symmetry == string("general") );
    const bool isSymmetric = ( symmetry == string("symmetric") );
    const bool isSkewSymmetric = ( symmetry == string("skew-symmetric") );
    const bool isHermitian = ( symmetry == string("hermitian") );
    if( !isMatrix && object != string("vector") )
        RuntimeError("Invalid Matrix Market object: ",object);
    if( !isArray && format != string("coordinate") )
        RuntimeError("Invalid Matrix Market format: ",format);
    if( !isComplex && !isPattern && 
        field != string("real") && 
        field != string("double") &&
        field != string("integer") )
        RuntimeError("Invalid Matrix Market field: ",field);
    if( !isGeneral && !isSymmetric && !isSkewSymmetric && !isHermitian )
        RuntimeError("Invalid Matrix Market symmetry: ",symmetry);
    // Ensure that the components are consistent
    // -----------------------------------------
    if( isArray && isPattern )
        RuntimeError("Pattern field requires coordinate format");
    // NOTE: This constraint is only enforced because of the note located at
    //       http://people.sc.fsu.edu/~jburkardt/data/mm/mm.html
    if( isSkewSymmetric && isPattern )
        RuntimeError("Pattern field incompatible with skew-symmetry");
    if( isHermitian && !isComplex )
        RuntimeError("Hermitian symmetry requires complex data");

    // Skip the comment lines
    // ======================
    while( file.peek() == '%' ) 
        std::getline( file, line );
  
    int m, n;
    if( !std::getline( file, line ) )
        RuntimeError("Could not extract the size line");
    if( isArray )
    {
        // Read in the matrix dimensions
        // =============================
        if( isMatrix )
        {
            std::stringstream lineStream( line );
            if( !(lineStream >> m) )
                RuntimeError("Missing matrix height: ",line);
            if( !(lineStream >> n) )
                RuntimeError("Missing matrix width: ",line);
        }
        else
        {
            std::stringstream lineStream( line );
            if( !(lineStream >> m) )
                RuntimeError("Missing vector height: ",line);
            n = 1;
        }

        // Resize the matrix
        // =================
        Zeros( A, m, n );

        // Now read in the data
        // ====================
        Real realPart, imagPart;
        for( Int j=0; j<n; ++j )
        {
            for( Int i=0; i<m; ++i )
            {
                if( !std::getline( file, line ) )
                    RuntimeError("Could not get entry (",i,",",j,")");
                std::stringstream lineStream( line );
                if( !(lineStream >> realPart) )
                    RuntimeError
                    ("Could not extract real part of entry (",i,",",j,")");
                A.SetRealPart( i, j, realPart );
                if( isComplex )
                {
                    if( !(lineStream >> imagPart) )
                        RuntimeError
                        ("Could not extract imag part of entry (",i,",",j,")");
                    A.SetImagPart( i, j, imagPart );
                }
            }
        }
    }
    else
    {
        // Read in the matrix dimensions and number of nonzeros
        // ====================================================
        int numNonzero;
        if( isMatrix )
        {
            std::stringstream lineStream( line );
            if( !(lineStream >> m) )
                RuntimeError("Missing matrix height: ",line);
            if( !(lineStream >> n) )
                RuntimeError("Missing matrix width: ",line);
            if( !(lineStream >> numNonzero) )
                RuntimeError("Missing nonzeros entry: ",line);
        }
        else
        {
            std::stringstream lineStream( line );
            if( !(lineStream >> m) )
                RuntimeError("Missing vector height: ",line);
            n = 1;
            if( !(lineStream >> numNonzero) )
                RuntimeError("Missing nonzeros entry: ",line);
        }

        // Create a matrix of zeros
        // ========================
        Zeros( A, m, n );

        // Fill in the nonzero entries
        // ===========================
        int i, j;
        Real realPart, imagPart;
        for( Int k=0; k<numNonzero; ++k )
        {
            if( !std::getline( file, line ) )
                RuntimeError("Could not get nonzero ",k);
            std::stringstream lineStream( line );
            if( !(lineStream >> i) )
                RuntimeError("Could not extract row coordinate of nonzero ",k);
            --i; // convert from Fortran to C indexing
            if( isMatrix )
            {
                if( !(lineStream >> j) )
                    RuntimeError
                    ("Could not extract col coordinate of nonzero ",k);
                --j;
            }
            else
                j = 0;

            if( isPattern )
            {
                A.Set( i, j, T(1) );
            }
            else
            {
                if( !(lineStream >> realPart) )
                    RuntimeError
                    ("Could not extract real part of entry (",i,",",j,")");
                A.UpdateRealPart( i, j, realPart );
                if( isComplex )
                {
                    if( !(lineStream >> imagPart) )
                        RuntimeError
                        ("Could not extract imag part of entry (",i,",",j,")");
                    A.UpdateImagPart( i, j, imagPart );
                }
            }
        }
    }

    if( isSymmetric )
        MakeSymmetric( LOWER, A );
    if( isHermitian )
        MakeHermitian( LOWER, A );
    // I'm not certain of what the MM standard is for complex skew-symmetry,
    // so I'll default to assuming no conjugation
    const bool conjugateSkew = false;
    if( isSkewSymmetric )
    {
        MakeSymmetric( LOWER, A, conjugateSkew );
        ScaleTrapezoid( T(-1), UPPER, A, 1 );
    }
}

template<typename T>
void MatrixMarket( AbstractDistMatrix<T>& A, const string filename )
{
    DEBUG_CSE
    // TODO: Use a WriteProxy instead
    DistMatrix<T,CIRC,CIRC> A_CIRC_CIRC( A.Grid() );
    if( A_CIRC_CIRC.CrossRank() == A_CIRC_CIRC.Root() )
    {
        MatrixMarket( A_CIRC_CIRC.Matrix(), filename );
        A_CIRC_CIRC.Resize
        ( A_CIRC_CIRC.Matrix().Height(), A_CIRC_CIRC.Matrix().Width() );
    }
    A_CIRC_CIRC.MakeSizeConsistent();
    Copy( A_CIRC_CIRC, A );
}

template<typename T>
void MatrixMarket( SparseMatrix<T>& A, const string filename )
{
    DEBUG_CSE
    typedef Base<T> Real;
    std::ifstream file( filename.c_str() );
    if( !file.is_open() )
        RuntimeError("Could not open ",filename);

    // Read the header
    // ===============
    // Attempt to pull in the various header components
    // ------------------------------------------------
    string line, stamp, object, format, field, symmetry;
    if( !std::getline( file, line ) )
        RuntimeError("Could not extract header line");
    {
        std::stringstream lineStream( line );
        lineStream >> stamp; 
        if( stamp != string("%%MatrixMarket") )
            RuntimeError("Invalid Matrix Market stamp: ",stamp);
        if( !(lineStream >> object) ) 
            RuntimeError("Missing Matrix Market object");
        if( !(lineStream >> format) )
            RuntimeError("Missing Matrix Market format");
        if( !(lineStream >> field) )
            RuntimeError("Missing Matrix Market field");
        if( !(lineStream >> symmetry) )
            RuntimeError("Missing Matrix Market symmetry");
    }
    // Ensure that the header components are individually valid
    // --------------------------------------------------------
    const bool isMatrix = ( object == string("matrix") );
    const bool isArray = ( format == string("array") );
    const bool isComplex = ( field == string("complex") );
    const bool isPattern = ( field == string("pattern") );
    const bool isGeneral = ( symmetry == string("general") );
    const bool isSymmetric = ( symmetry == string("symmetric") );
    const bool isSkewSymmetric = ( symmetry == string("skew-symmetric") );
    const bool isHermitian = ( symmetry == string("hermitian") );
    if( !isMatrix && object != string("vector") )
        RuntimeError("Invalid Matrix Market object: ",object);
    if( !isArray && format != string("coordinate") )
        RuntimeError("Invalid Matrix Market format: ",format);
    if( !isComplex && !isPattern && 
        field != string("real") && 
        field != string("double") &&
        field != string("integer") )
        RuntimeError("Invalid Matrix Market field: ",field);
    if( !isGeneral && !isSymmetric && !isSkewSymmetric && !isHermitian )
        RuntimeError("Invalid Matrix Market symmetry: ",symmetry);
    // Ensure that the components are consistent
    // -----------------------------------------
    if( isArray && isPattern )
        RuntimeError("Pattern field requires coordinate format");
    // NOTE: This constraint is only enforced because of the note located at
    //       http://people.sc.fsu.edu/~jburkardt/data/mm/mm.html
    if( isSkewSymmetric && isPattern )
        RuntimeError("Pattern field incompatible with skew-symmetry");
    if( isHermitian && !isComplex )
        RuntimeError("Hermitian symmetry requires complex data");

    if( isArray )
    {
        LogicError
        ("Attempted to load dense MatrixMarket format into SparseMatrix");
    }

    // Skip the comment lines
    // ======================
    while( file.peek() == '%' ) 
        std::getline( file, line );
  
    int m, n;
    if( !std::getline( file, line ) )
        RuntimeError("Could not extract the size line");

    // Read in the matrix dimensions and number of nonzeros
    // ====================================================
    int numNonzero;
    if( isMatrix )
    {
        std::stringstream lineStream( line );
        if( !(lineStream >> m) )
            RuntimeError("Missing matrix height: ",line);
        if( !(lineStream >> n) )
            RuntimeError("Missing matrix width: ",line);
        if( !(lineStream >> numNonzero) )
            RuntimeError("Missing nonzeros entry: ",line);
    }
    else
    {
        std::stringstream lineStream( line );
        if( !(lineStream >> m) )
            RuntimeError("Missing vector height: ",line);
        n = 1;
        if( !(lineStream >> numNonzero) )
            RuntimeError("Missing nonzeros entry: ",line);
    }

    // Create a matrix of zeros
    // ========================
    Zeros( A, m, n );

    // Fill in the nonzero entries
    // ===========================
    int i, j;
    Real realPart, imagPart;
    T value;
    A.Reserve( numNonzero );
    for( Int k=0; k<numNonzero; ++k )
    {
        if( !std::getline( file, line ) )
            RuntimeError("Could not get nonzero ",k);
        std::stringstream lineStream( line );
        if( !(lineStream >> i) )
            RuntimeError("Could not extract row coordinate of nonzero ",k);
        --i; // convert from Fortran to C indexing
        if( isMatrix )
        {
            if( !(lineStream >> j) )
                RuntimeError
                ("Could not extract col coordinate of nonzero ",k);
            --j;
        }
        else
            j = 0;

        if( isPattern )
        {
            A.QueueUpdate( i, j, T(1) );
        }
        else
        {
            if( isComplex )
            {
                if( !(lineStream >> realPart) )
                    RuntimeError
                    ("Could not extract real part of entry (",i,",",j,")");
                if( !(lineStream >> imagPart) )
                    RuntimeError
                    ("Could not extract imag part of entry (",i,",",j,")");
                SetRealPart( value, realPart );
                SetImagPart( value, imagPart );
                A.QueueUpdate( i, j, value );
            }
            else
            {
                if( !(lineStream >> realPart) )
                    RuntimeError
                    ("Could not extract real entry (",i,",",j,")");
                A.QueueUpdate( i, j, T(realPart) );
            }
        }
    }
    A.ProcessQueues();

    if( isSymmetric )
        MakeSymmetric( LOWER, A );
    if( isHermitian )
        MakeHermitian( LOWER, A );
    // I'm not certain of what the MM standard is for complex skew-symmetry,
    // so I'll default to assuming no conjugation
    const bool conjugateSkew = false;
    if( isSkewSymmetric )
    {
        MakeSymmetric( LOWER, A, conjugateSkew );
        ScaleTrapezoid( T(-1), UPPER, A, 1 );
    }
}

template<typename T>
void MatrixMarket( DistSparseMatrix<T>& A, const string filename )
{
    DEBUG_CSE
    typedef Base<T> Real;
    std::ifstream file( filename.c_str() );
    if( !file.is_open() )
        RuntimeError("Could not open ",filename);

    // Read the header
    // ===============
    // Attempt to pull in the various header components
    // ------------------------------------------------
    string line, stamp, object, format, field, symmetry;
    if( !std::getline( file, line ) )
        RuntimeError("Could not extract header line");
    {
        std::stringstream lineStream( line );
        lineStream >> stamp; 
        if( stamp != string("%%MatrixMarket") )
            RuntimeError("Invalid Matrix Market stamp: ",stamp);
        if( !(lineStream >> object) ) 
            RuntimeError("Missing Matrix Market object");
        if( !(lineStream >> format) )
            RuntimeError("Missing Matrix Market format");
        if( !(lineStream >> field) )
            RuntimeError("Missing Matrix Market field");
        if( !(lineStream >> symmetry) )
            RuntimeError("Missing Matrix Market symmetry");
    }
    // Ensure that the header components are individually valid
    // --------------------------------------------------------
    const bool isMatrix = ( object == string("matrix") );
    const bool isArray = ( format == string("array") );
    const bool isComplex = ( field == string("complex") );
    const bool isPattern = ( field == string("pattern") );
    const bool isGeneral = ( symmetry == string("general") );
    const bool isSymmetric = ( symmetry == string("symmetric") );
    const bool isSkewSymmetric = ( symmetry == string("skew-symmetric") );
    const bool isHermitian = ( symmetry == string("hermitian") );
    if( !isMatrix && object != string("vector") )
        RuntimeError("Invalid Matrix Market object: ",object);
    if( !isArray && format != string("coordinate") )
        RuntimeError("Invalid Matrix Market format: ",format);
    if( !isComplex && !isPattern && 
        field != string("real") && 
        field != string("double") &&
        field != string("integer") )
        RuntimeError("Invalid Matrix Market field: ",field);
    if( !isGeneral && !isSymmetric && !isSkewSymmetric && !isHermitian )
        RuntimeError("Invalid Matrix Market symmetry: ",symmetry);
    // Ensure that the components are consistent
    // -----------------------------------------
    if( isArray && isPattern )
        RuntimeError("Pattern field requires coordinate format");
    // NOTE: This constraint is only enforced because of the note located at
    //       http://people.sc.fsu.edu/~jburkardt/data/mm/mm.html
    if( isSkewSymmetric && isPattern )
        RuntimeError("Pattern field incompatible with skew-symmetry");
    if( isHermitian && !isComplex )
        RuntimeError("Hermitian symmetry requires complex data");

    if( isArray )
    {
        LogicError
        ("Attempted to load dense MatrixMarket format into SparseMatrix");
    }

    // Skip the comment lines
    // ======================
    while( file.peek() == '%' ) 
        std::getline( file, line );
  
    int m, n;
    if( !std::getline( file, line ) )
        RuntimeError("Could not extract the size line");

    // Read in the matrix dimensions and number of nonzeros
    // ====================================================
    int numNonzero;
    if( isMatrix )
    {
        std::stringstream lineStream( line );
        if( !(lineStream >> m) )
            RuntimeError("Missing matrix height: ",line);
        if( !(lineStream >> n) )
            RuntimeError("Missing matrix width: ",line);
        if( !(lineStream >> numNonzero) )
            RuntimeError("Missing nonzeros entry: ",line);
    }
    else
    {
        std::stringstream lineStream( line );
        if( !(lineStream >> m) )
            RuntimeError("Missing vector height: ",line);
        n = 1;
        if( !(lineStream >> numNonzero) )
            RuntimeError("Missing nonzeros entry: ",line);
    }

    // Create a matrix of zeros
    // ========================
    Zeros( A, m, n );

    // Fill in the local nonzero entries by passively calling QueueUpdate
    // ==================================================================
    int i, j;
    Real realPart, imagPart;
    T value;
    const int commSize = mpi::Size(A.Comm());
    A.Reserve( numNonzero/commSize ); // Assume an even nonzero distribution
    bool passive = true;
    for( Int k=0; k<numNonzero; ++k )
    {
        if( !std::getline( file, line ) )
            RuntimeError("Could not get nonzero ",k);
        std::stringstream lineStream( line );
        if( !(lineStream >> i) )
            RuntimeError("Could not extract row coordinate of nonzero ",k);
        --i; // convert from Fortran to C indexing
        if( isMatrix )
        {
            if( !(lineStream >> j) )
                RuntimeError
                ("Could not extract col coordinate of nonzero ",k);
            --j;
        }
        else
            j = 0;

        if( isPattern )
        {
            A.QueueUpdate( i, j, T(1), passive );
        }
        else
        {
            if( isComplex )
            {
                if( !(lineStream >> realPart) )
                    RuntimeError
                    ("Could not extract real part of entry (",i,",",j,")");
                if( !(lineStream >> imagPart) )
                    RuntimeError
                    ("Could not extract imag part of entry (",i,",",j,")");
                SetRealPart( value, realPart );
                SetImagPart( value, imagPart );
                A.QueueUpdate( i, j, value, passive );
            }
            else
            {
                if( !(lineStream >> realPart) )
                    RuntimeError
                    ("Could not extract real entry (",i,",",j,")");
                A.QueueUpdate( i, j, T(realPart), passive );
            }
        }
    }
    A.ProcessLocalQueues();

    if( isSymmetric )
    {
        MakeSymmetric( LOWER, A );
    }
    if( isHermitian )
    {
        MakeHermitian( LOWER, A );
    }

    // I'm not certain of what the MM standard is for complex skew-symmetry,
    // so I'll default to assuming no conjugation
    const bool conjugateSkew = false;
    if( isSkewSymmetric )
    {
        MakeSymmetric( LOWER, A, conjugateSkew );
        ScaleTrapezoid( T(-1), UPPER, A, 1 );
    }
}

} // namespace read
} // namespace El

#endif // ifndef EL_READ_MATRIXMARKET_HPP
