%
% Copyright (c) 2009-2014, Jack Poulson
% All rights reserved.
%
% This file is part of Elemental and is under the BSD 2-Clause License, 
% which can be found in the LICENSE file in the root directory, or at 
% http://opensource.org/licenses/BSD-2-Clause
% 
n=100;

A=randn(n,n);
[L,U,P]=lu(A);
p=P*(1:n)';

x=randn(n,1);
y=A*x;
xMan=U\(L\y(p));
norm(x-xMan)/norm(x)

u=randn(n,1);
v=randn(n,1);
tau=0.2;
[LNew,UNew,pNew]=LUModSimple(L,U,p,u,v,tau);

x=randn(n,1);
y=(A+u*v')*x;
xMan=UNew\(LNew\y(pNew));
norm(x-xMan)/norm(x)
