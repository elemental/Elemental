/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/

#include <El.hpp>
using namespace El;

template <typename T, DistWrap W>
void TestColumnTwoNorms(Int m, Int n, const Grid& g, bool print) {
  // Generate random matrix to test.
  DistMatrix<T, MC, MR, W> A(g);
  Uniform(A, m, n);
  if (print) {
    Print(A, "A");
  }
  DistMatrix<T, MR, STAR, W> norms(g);
  ColumnTwoNorms(A, norms);
  if (print) {
    Print(norms, "norms");
  }
  for (Int j = 0; j < A.LocalWidth(); ++j) {
    T got = norms.GetLocal(j, 0);
    T expected = 0;
    for (Int i = 0; i < A.LocalHeight(); ++i) {
      T val = A.GetLocal(i, j);
      expected += val * val;
    }
    expected = mpi::AllReduce(expected, g.ColComm());
    expected = Sqrt(expected);
    if (Abs(got - expected) > 1e-5) {
      Output("Results do not match, norms(", j, ")=", got,
             " instead of ", expected);
      RuntimeError("got != expected");
    }
  }
}

template <typename T, DistWrap W>
void TestColumnMaxNorms(Int m, Int n, const Grid& g, bool print) {
  // Generate random matrix to test.
  DistMatrix<T, MC, MR, W> A(g);
  Uniform(A, m, n);
  if (print) {
    Print(A, "A");
  }
  DistMatrix<T, MR, STAR, W> norms(g);
  ColumnMaxNorms(A, norms);
  if (print) {
    Print(norms, "norms");
  }
  for (Int j = 0; j < A.LocalWidth(); ++j) {
    T got = norms.GetLocal(j, 0);
    T expected = 0;
    for (Int i = 0; i < A.LocalHeight(); ++i) {
      expected = Max(expected, Abs(A.GetLocal(i, j)));
    }
    T r;
    mpi::AllReduce(&expected, &r, 1, mpi::MAX, g.ColComm());
    expected = r;
    if (got != expected) {
      Output("Results do not match, norms(", j, ")=", got,
             " instead of ", expected);
      RuntimeError("got != expected");
    }
  }
}

int main(int argc, char** argv) {
  Environment env(argc, argv);
  mpi::Comm comm = mpi::COMM_WORLD;
  try {
    const Int m = Input("--m", "height", 100);
    const Int n = Input("--n", "width", 100);
    const bool print = Input("--print", "print matrices?", false);
    ProcessInput();
    PrintInputReport();

    const Grid g(comm);
    OutputFromRoot(comm, "Testing ColumnTwoNorms");
    TestColumnTwoNorms<float, ELEMENT>(m, n, g, print);
    TestColumnTwoNorms<float, BLOCK>(m, n, g, print);
    TestColumnTwoNorms<double, ELEMENT>(m, n, g, print);
    TestColumnTwoNorms<double, BLOCK>(m, n, g, print);
    OutputFromRoot(comm, "Testing ColumnMaxNorms");
    TestColumnMaxNorms<float, ELEMENT>(m, n, g, print);
    TestColumnMaxNorms<float, BLOCK>(m, n, g, print);
    TestColumnMaxNorms<double, ELEMENT>(m, n, g, print);
    TestColumnMaxNorms<double, BLOCK>(m, n, g, print);
  } catch (exception& e) {
    ReportException(e);
  }
}
