/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_WRITE_MATRIXMARKET_HPP
#define ELEM_WRITE_MATRIXMARKET_HPP

namespace elem {
namespace write {

template<typename T>
inline void
MatrixMarket( const Matrix<T>& A, std::string basename="matrix" )
{
    DEBUG_ONLY(CallStackEntry cse("write::MatrixMarket"))
    
    std::string filename = basename + "." + FileExtension(MATRIX_MARKET);
    std::ofstream file( filename.c_str(), std::ios::binary );
    if( !file.is_open() )
        RuntimeError("Could not open ",filename);

    // Write the header
    // ================
    {
        std::ostringstream os;
        os << "%%MatrixMarket matrix array ";
        if( IsComplex<T>::val )
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
    {
        std::ostringstream os; 
        os << m << " " << n << "\n";
        file << os.str();
    }
    
    // Write the entries
    // =================
    for( Int j=0; j<n; ++j )
    {
        for( Int i=0; i<m; ++i )
        {
            std::ostringstream os;
            os << A.GetRealPart(i,j);
            if( IsComplex<T>::val )
                os << " " << A.GetImagPart(i,j);
            os << "\n";
            file << os.str();
        }
    }
}

} // namespace write
} // namespace elem

#endif // ifndef ELEM_WRITE_MATRIXMARKET_HPP
