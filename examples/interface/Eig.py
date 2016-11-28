import El

n = 1000
output = True

A = El.Matrix(El.dTag)
El.Uniform( A, n, n, 0., 10.)
if output:
  El.Print( A, "AOrig" )

w, X = El.Eig( A )

if output:
  El.Print( A, "A" )
  El.Print( w, "w" )
  El.Print( X, "X" )

El.Finalize()
