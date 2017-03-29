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
void TestDot(Int m, Int n, const Grid& g, bool print)
{
  // Generate random matrices to test.
  DistMatrix<T, MC, MR, W> A(g);
  Uniform(A, m, n);
  DistMatrix<T, MC, MR, W> B(g);
  Uniform(B, m, n);
  if (print)
  {
    Print(A, "A");
    Print(B, "B");
  }
  // Do Dot.
  T got = Dot(A, B);
  if (print)
    OutputFromRoot(g.Comm(), "result=", got);
  // Manually check results.
  T expected = 0;
  for (Int j = 0; j < A.LocalWidth(); ++j)
    for (Int i = 0; i < B.LocalHeight(); ++i)
      expected += Conj(A.GetLocal(i, j)) * B.GetLocal(i, j);
  expected = mpi::AllReduce(expected, g.Comm());
  // The constant here is large because this is not an especially stable way
  // to compute the dot product, but it provides a dumb implementation baseline.
  if (Abs(got - expected) > 700 * limits::Epsilon<El::Base<T>>())
  {
    Output("Results do not match, got=", got,
           " instead of ", expected);
    RuntimeError("got != expected");
  }
}

int main(int argc, char** argv)
{
  Environment env(argc, argv);
  mpi::Comm comm = mpi::COMM_WORLD;
  try
  {
    const Int m = Input("--m", "height", 100);
    const Int n = Input("--n", "width", 100);
    const bool print = Input("--print", "print matrices?", false);
    ProcessInput();
    PrintInputReport();

    const Grid g(comm);
    OutputFromRoot(comm, "Testing Dot");
    TestDot<float, ELEMENT>(m, n, g, print);
    TestDot<float, BLOCK>(m, n, g, print);
    TestDot<Complex<float>, ELEMENT>(m, n, g, print);
    TestDot<Complex<float>, BLOCK>(m, n, g, print);
    TestDot<double, ELEMENT>(m, n, g, print);
    TestDot<double, BLOCK>(m, n, g, print);
    TestDot<Complex<double>, ELEMENT>(m, n, g, print);
    TestDot<Complex<double>, BLOCK>(m, n, g, print);
#if defined(EL_HAVE_QD)
    TestDot<DoubleDouble, ELEMENT>(m, n, g, print);
    TestDot<DoubleDouble, BLOCK>(m, n, g, print);
    TestDot<Complex<DoubleDouble>, ELEMENT>(m, n, g, print);
    TestDot<Complex<DoubleDouble>, BLOCK>(m, n, g, print);
    TestDot<QuadDouble, ELEMENT>(m, n, g, print);
    TestDot<QuadDouble, BLOCK>(m, n, g, print);
    TestDot<Complex<QuadDouble>, ELEMENT>(m, n, g, print);
    TestDot<Complex<QuadDouble>, BLOCK>(m, n, g, print);
#endif
#if defined(EL_HAVE_QUAD)
    TestDot<Quad, ELEMENT>(m, n, g, print);
    TestDot<Quad, BLOCK>(m, n, g, print);
    TestDot<Complex<Quad>, ELEMENT>(m, n, g, print);
    TestDot<Complex<Quad>, BLOCK>(m, n, g, print);
#endif
#if defined(EL_HAVE_MPC)
    TestDot<BigFloat, ELEMENT>(m, n, g, print);
    TestDot<BigFloat, BLOCK>(m, n, g, print);
    TestDot<Complex<BigFloat>, ELEMENT>(m, n, g, print);
    TestDot<Complex<BigFloat>, BLOCK>(m, n, g, print);
#endif
  }
  catch (exception& e)
  {
    ReportException(e);
  }
}
