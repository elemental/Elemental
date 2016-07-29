/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_WRITE_MATRIXMARKET_HPP
#define EL_WRITE_MATRIXMARKET_HPP

namespace El {
namespace write {

template<typename T>
void MatrixMarket( const Matrix<T>& A, string basename="matrix" )
{
    DEBUG_CSE
    
    string filename = basename + "." + FileExtension(MATRIX_MARKET);
    ofstream file( filename.c_str(), std::ios::binary );
    if( !file.is_open() )
        RuntimeError("Could not open ",filename);

    // Write the header
    // ================
    {
        ostringstream os;
        os << "%%MatrixMarket matrix array ";
        if( IsComplex<T>::value )
            os << "complex "; 
        else
            os << "real ";
        os << "general\n";
        file << os.str();
    }
    
    // Write the size line
    // ===================
    const Int m = A.Height();
    const Int n = A.Width();
    file << BuildString(m," ",n,"\n");
    
    // Write the entries
    // =================
    for( Int j=0; j<n; ++j )
    {
        for( Int i=0; i<m; ++i )
        {
            ostringstream os;
            os << A.GetRealPart(i,j);
            if( IsComplex<T>::value )
                os << " " << A.GetImagPart(i,j);
            os << "\n";
            file << os.str();
        }
    }
}

template<typename T>
void MatrixMarket( const SparseMatrix<T>& A, string basename="matrix" )
{
    DEBUG_CSE
    
    string filename = basename + "." + FileExtension(MATRIX_MARKET);
    ofstream file( filename.c_str(), std::ios::binary );
    if( !file.is_open() )
        RuntimeError("Could not open ",filename);

    // Write the header
    // ================
    {
        ostringstream os;
        os << "%%MatrixMarket matrix coordinate ";
        if( IsComplex<T>::value )
            os << "complex "; 
        else
            os << "real ";
        os << "general\n";
        file << os.str();
    }
    
    // Write the size line
    // ===================
    const Int m = A.Height();
    const Int n = A.Width();
    const Int numNonzeros = A.NumEntries();
    file << BuildString(m," ",n," ",numNonzeros,"\n");
    
    // Write the entries
    // =================
    for( Int e=0; e<numNonzeros; ++e )
    {
        const Int i = A.Row(e);
        const Int j = A.Col(e);
        const T value = A.Value(e);
        if( IsComplex<T>::value )
        {
            file << 
              BuildString(i," ",j," ",RealPart(value)," ",ImagPart(value),"\n");
        }
        else
        {
            file << BuildString(i," ",j," ",RealPart(value),"\n");
        }
    }
}

} // namespace write
} // namespace El

#endif // ifndef EL_WRITE_MATRIXMARKET_HPP
