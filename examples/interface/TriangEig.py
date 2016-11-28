import El

n = 1000
outputMatrices = False
outputVector = True

A = El.Matrix(El.zTag)
El.Uniform( A, n, n, 0., 10.)

T = El.Matrix(El.zTag)
El.Copy(A,T)
w, Q = El.Schur(T,fullTriangle=True,vectors=True)

X = El.TriangEig(T)

if outputMatrices:
  El.Print( A, "A" )
  El.Print( Q, "Q" )
  El.Print( T, "T" )
  El.Print( X, "X" )
if outputVector:
  El.Print( w, "w" )

El.Finalize()
