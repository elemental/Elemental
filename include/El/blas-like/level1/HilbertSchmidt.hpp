/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_HILBERTSCHMIDT_HPP
#define EL_HILBERTSCHMIDT_HPP

namespace El {

template<typename F> 
F HilbertSchmidt( const Matrix<F>& A, const Matrix<F>& B );

template<typename F>
F HilbertSchmidt
( const AbstractDistMatrix<F>& A, const AbstractDistMatrix<F>& C );

} // namespace El

#endif // ifndef EL_HILBERTSCHMIDT_HPP
