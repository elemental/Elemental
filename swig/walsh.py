from elem import *

# Create an eighth-order Walsh matrix and Display it
A = DistMatrix_d()
Walsh( A, 8 )
Display( A, "Walsh matrix" )
ProcessEvents( 800 )

# Overwrite it with the L from its LDL^H factorization
LDLH( A )
MakeTriangular( LOWER, A )
SetDiagonal( A, 1 )
Display( A, "Lower-triangular factor" )
