%
% Copyright (c) 2009-2014, Jack Poulson
% All rights reserved.
%
% This file is part of Elemental and is under the BSD 2-Clause License, 
% which can be found in the LICENSE file in the root directory, or at 
% http://opensource.org/licenses/BSD-2-Clause
% 
n=3000;

A=randn(n,n);
tic;
[L,U,P]=lu(A);
toc;
APack=tril(L,-1)+U;
p=P*(1:n)';

x=randn(n,1);
y=A*x;
xMan=U\(L\y(p));
origRelResid=norm(x-xMan)/norm(x);
fprintf('||eOrig||_2 / ||xOrig||_2 = %e\n',origRelResid);

u=randn(n,1);
v=randn(n,1);
tau=0.1;
tic;
[APackNew,pNew]=LUMod(APack,p,u,v,tau);
toc;
LNew=tril(APackNew,-1)+eye(n,n);
UNew=triu(APackNew);

x=randn(n,1);
y=(A+u*v')*x;
xMan=UNew\(LNew\y(pNew));
relResid=norm(x-xMan)/norm(x);
fprintf('||e||_2 / ||x||_2 = %e\n',relResid);
