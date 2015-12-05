import El

n = 1000
output = True

A = El.Matrix(El.zTag)
El.Uniform( A, n, n, 0., 10.)

T = El.Matrix(El.zTag)
El.Copy(A,T)
w, Q = El.Schur(T,fullTriangle=True,vectors=True)

X = El.TriangEig(T)

if output:
  El.Print( A, "A" )
  El.Print( w, "w" )
  El.Print( Q, "Q" )
  El.Print( T, "T" )
  El.Print( X, "X" )
