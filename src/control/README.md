### `src/control/`

A few matrix sign function based solvers for control theory:

-  `Lyapunov.hpp`: Solves A X + X A' = C for X when A has its eigenvalues
   in the open right-half plane
-  `Riccati.hpp`: Solves X K X - A' X - X A = L for X when K and L are 
   Hermitian.
-  `Sylvester.hpp`: Solves A X + X B = C for X when A and B both have all of 
   their eigenvalues in the open right-half plane

#### TODO

Implement algorithms from Benner, Quintana-Orti, and Quintana-Orti's 
"Solving Stable Sylvester Equations via Rational Iterative Schemes".
