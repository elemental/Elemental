/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_READ_MATRIXMARKET_HPP
#define EL_READ_MATRIXMARKET_HPP

#include EL_ZEROS_INC

namespace El {
namespace read {

template<typename T>
inline void
MatrixMarket( Matrix<T>& A, const std::string filename )
{
    DEBUG_ONLY(CallStackEntry cse("read::MatrixMarket"))
    typedef Base<T> Real;
    std::ifstream file( filename.c_str() );
    if( !file.is_open() )
        RuntimeError("Could not open ",filename);

    // Read the header
    // ===============
    // Attempt to pull in the various header components
    // ------------------------------------------------
    std::string line, stamp, object, format, field, symmetry;
    if( !std::getline( file, line ) )
        RuntimeError("Could not extract header line");
    {
        std::stringstream lineStream( line );
        lineStream >> stamp; 
        if( stamp != std::string("%%MatrixMarket") )
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
    const bool isMatrix = ( object == std::string("matrix") );
    const bool isArray = ( format == std::string("array") );
    const bool isComplex = ( field == std::string("complex") );
    const bool isPattern = ( field == std::string("pattern") );
    const bool isGeneral = ( symmetry == std::string("general") );
    const bool isSymmetric = ( symmetry == std::string("symmetric") );
    const bool isSkewSymmetric = ( symmetry == std::string("skew-symmetric") );
    const bool isHermitian = ( symmetry == std::string("hermitian") );
    if( !isMatrix && object != std::string("vector") )
        RuntimeError("Invalid Matrix Market object: ",object);
    if( !isArray && format != std::string("coordinate") )
        RuntimeError("Invalid Matrix Market format: ",format);
    if( !isComplex && !isPattern && 
        field != std::string("real") && 
        field != std::string("double") &&
        field != std::string("integer") )
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

template<typename T,Dist U,Dist V>
inline void
MatrixMarket( DistMatrix<T,U,V>& A, const std::string filename )
{
    DEBUG_ONLY(CallStackEntry cse("read::MatrixMarket"))
    DistMatrix<T,CIRC,CIRC> A_CIRC_CIRC( A.Grid() );
    if( A_CIRC_CIRC.CrossRank() == A_CIRC_CIRC.Root() )
    {
        MatrixMarket( A_CIRC_CIRC.Matrix(), filename );
        A_CIRC_CIRC.Resize
        ( A_CIRC_CIRC.Matrix().Height(), A_CIRC_CIRC.Matrix().Width() );
    }
    A_CIRC_CIRC.MakeSizeConsistent();
    A = A_CIRC_CIRC;
}

template<typename T,Dist U,Dist V>
inline void
MatrixMarket( BlockDistMatrix<T,U,V>& A, const std::string filename )
{
    DEBUG_ONLY(CallStackEntry cse("read::MatrixMarket"))
    BlockDistMatrix<T,CIRC,CIRC> A_CIRC_CIRC( A.Grid() );
    if( A_CIRC_CIRC.CrossRank() == A_CIRC_CIRC.Root() )
    {
        MatrixMarket( A_CIRC_CIRC.Matrix(), filename );
        A_CIRC_CIRC.Resize
        ( A_CIRC_CIRC.Matrix().Height(), A_CIRC_CIRC.Matrix().Width() );
    }
    A_CIRC_CIRC.MakeSizeConsistent();
    A = A_CIRC_CIRC;
}

} // namespace read
} // namespace El

#endif // ifndef EL_READ_MATRIXMARKET_HPP
