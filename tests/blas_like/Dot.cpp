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
void TestDot(Int m, Int n, const Grid& g, bool print) {
  // Generate random matrices to test.
  DistMatrix<T, MC, MR, W> A(g);
  Uniform(A, m, n);
  DistMatrix<T, MC, MR, W> B(g);
  Uniform(B, m, n);
  if (print) {
    Print(A, "A");
    Print(B, "B");
  }
  // Do Dot.
  T got = Dot(A, B);
  if (print) {
    OutputFromRoot(g.Comm(), "result=", got);
  }
  // Manually check results.
  T expected = 0;
  for (Int j = 0; j < A.LocalWidth(); ++j) {
    for (Int i = 0; i < B.LocalHeight(); ++i) {
      expected += A.GetLocal(i, j) * B.GetLocal(i, j);
    }
  }
  expected = mpi::AllReduce(expected, g.Comm());
  if (Abs(got - expected) > 1e-6) {
    Output("Results do not match, got=", got,
           " instead of ", expected);
    RuntimeError("got != expected");
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
    OutputFromRoot(comm, "Testing Dot");
    TestDot<float, ELEMENT>(m, n, g, print);
    TestDot<float, BLOCK>(m, n, g, print);
    TestDot<double, ELEMENT>(m, n, g, print);
    TestDot<double, BLOCK>(m, n, g, print);
  } catch (exception& e) {
    ReportException(e);
  }
}
