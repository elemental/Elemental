
/*
   Copyright (c) 2009-2016, Ryan H. Lewis
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_GRAPH_ALGORITHMS_SEIDEL_HPP
#define EL_GRAPH_ALGORITHMS_SEIDEL_HPP
namespace El {
/***
 * This algorithm is based off of the algorithm APD in the paper
 * ` On the all-pairs-shortest-path problem` by Raimund Seidel
 *  Found: http://dl.acm.org/citation.cfm?id=129784
 * 	   http://duch.mimuw.edu.pl/~mucha/teaching/alp2006/seidel92.pdf
 * TODO: Make this work with DistGraph input. 
 * TODO: Write Unit Tests
 * TODO: Verify Correctness.
 */
template< typename T>
DistMatrix<T> Seidel( const DistMatrix<T>& A){
	DistMatrix<T> Z( A.Height(), A.Width(), A.Grid());
	Multiply( El::NORMAL, A, A, T(0), Z);
	DistMatrix<T> B( A);
	bool all_ones = true;
	for(Int i = 0; i < B.LocalHeight(); ++i){
		for(Int j = i+1; j < B.LocalWidth(); ++j){
			if( Z.Get(B.GlobalRow( i), B.GlobalColumn( j)) > 0){ 
				B.SetLocal( i,j, T(1)); 
			} else {
			   all_ones = false;
			}
		}
	}
	auto comm = B.Comm();
	auto rank = mpi::Rank( comm);
	all_ones = mpi::AllReduce( all_ones, [](bool& x, bool& y){ return x & y; }, true, comm);
	//Set b_ij = 1 iff i neq j and (a_ij = 1 or Z_ij > 0)	
	if( all_ones){
		//B = 2*B - A;
		Axpy(T(1), B, B);
		Axpy(T(-1), A, B);
		return B;
	}
	Z = Seidel( B);
	Multiply( El::NORMAL, Z, A, T(0), B);
	DistMatrix<T> Deg(A.Height(), 1, A.Grid());
	DistMatrix<T> O(A.Width(), 1, A.Grid());
	for(Int i = 0; i < O.Width(); ++i){  O.Set(i,0, 1); }
 	Gemv( El::NORMAL, A, O, T(0), Deg);
	for( Int i = 0; i < B.LocalHeight(); ++i){
	  for( Int j = 0; j < B.LocalWidth(); ++j){
		auto zij = Z.Get(B.GlobalRow( i), B.GlobalColumn( j));
	  	if( B.GetLocal(i,j) >= zij*O.Get(B.GlobalColumn(j),0)){
			B.SetLocal(i,j,T(2)*zij);
		} else {
			B.SetLocal(i,j,T(2)*zij - T(1));
		}
	   }
	}
	return B;
}

}
#endif // ifndef EL_GRAPH_ALGORITHMS_SEIDEL_HPP
