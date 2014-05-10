function [A,p] = LUMod( A, p, u, v, tau )
%
% Copyright (c) 2009-2014, Jack Poulson
% All rights reserved.
%
% This file is part of Elemental and is under the BSD 2-Clause License, 
% which can be found in the LICENSE file in the root directory, or at 
% http://opensource.org/licenses/BSD-2-Clause
% 
% Begin with an LU factorization with partial pivoting, 
%     A = P^T L U,
% and turn it into a partially-pivoted LU factorization of
%     A + u v',
% say
%     (A + u v') = P^T L ( U + w v'),
% w = inv(L) P u.
%
% Please see subsection 2.1 from 
%     Peter Stange, Andreas Griewank, and Matthias Bollhofer,
%     "On the efficient update of rectangular LU factorizations subject to
%      low rank modifications"
% which discusses the technique of Schwetlick and Kielbasinski described in
% "Numerische Lineare Algebra".
%
[m,n] = size(A);
minDim = min(m,n);
if minDim ~= m,
  error('It is assumed that height(A) <= width(A)');
end

% w := inv(L) P u
w=u(p);
L=tril(A,-1)+diag(ones(minDim,1));
w=L\w;

% Maintain an external vector for the temporary subdiagonal of U
uSub=zeros(minDim-1,1);

% Reduce w to a multiple of e0
for i=minDim-1:-1:1,
  % Decide if we should pivot the i'th and i+1'th rows of w
  rightTerm = abs(A(i+1,i)*w(i)+w(i+1));
  if abs(w(i)) < tau*rightTerm,
    % P := P_i P
    p([i,i+1])=p([i+1,i]);

    % Simultaneously perform 
    %   U := P_i U and
    %   L := P_i L P_i^T
    %
    % Then update
    %     L := L T_{i,L}^{-1},
    %     U := T_{i,L} U, 
    %     w := T_{i,L} P_i w,
    % where T_{i,L} is the Gauss transform which zeros (P_i w)_{i+1}.
    % 
    % More succinctly,
    %     gamma    := w(i) / w(i+1),
    %     w(i)     := w(i+1), 
    %     w(i+1)   := 0,
    %     L(:,i)   += gamma L(:,i+1),
    %     U(i+1,:) -= gamma U(i,:).
    gamma = w(i) / w(i+1);
    uSub(i) = A(i,i);
    lambdaSub = A(i+1,i);
    lambda_ii = 1 + gamma*lambdaSub;
    w(i) = w(i+1);
    A(i,  i) = gamma;
    A(i+1,i) = 0;

    lBiCopy = A(i+2:end,i);
    A(i+2:end,[i,i+1]) = A(i+2:end,[i+1,i]);
    A(i+2:end,i) = A(i+2:end,i) + gamma*lBiCopy;

    uip1RCopy = A(i+1,i+1:end);
    A([i,i+1],:)=A([i+1,i],:);
    A(i+1,i+1:end) = A(i+1,i+1:end) - gamma*uip1RCopy;

    % Force L back to *unit* lower-triangular form via the transform
    %     L := L T_{i,U}^{-1} D^{-1}, 
    % where D is diagonal and responsible for forcing L(i,i) and 
    % L(i+1,i+1) back to 1. The effect on L is:
    %     eta       := L(i,i+1)/L(i,i),
    %     L(:,i+1)  -= eta L(:,i),
    %     delta_i   := L(i,i),
    %     delta_ip1 := L(i+1,i+1),
    %     L(:,i)   /= delta_i,
    %     L(:,i+1) /= delta_ip1,
    % while the effect on U is
    %     U(i,:)   += eta U(i+1,:)
    %     U(i,:)   *= delta_i,
    %     U(i+1,:) *= delta_{i+1},
    % and the effect on w is
    %     w(i) *= delta_i.
    eta = lambdaSub/lambda_ii;
    delta_i = lambda_ii;
    delta_ip1 = 1 - eta*gamma;
    A(i+2:end,i+1) = (A(i+2:end,i+1) - eta*A(i+2:end,i)) / delta_ip1;
    A(i+1:end,i  ) = A(i+1:end,i  ) / delta_i;

    A(i,i      ) = (               eta*uSub(i)       ) * delta_i;
    A(i,i+1:end) = (A(i,i+1:end) + eta*A(i+1,i+1:end)) * delta_i;
    uSub(i)        = uSub(i)       *delta_ip1;
    A(i+1,i+1:end) = A(i+1,i+1:end)*delta_ip1;

    w(i) = w(i)*delta_i;
  else
    % Update
    %     L := L T_{i,L}^{-1},
    %     U := T_{i,L} U, 
    %     w := T_{i,L} w,
    % where T_{i,L} is the Gauss transform which zeros w_{i+1}.
    %
    % More succinctly,
    %     gamma    := w(i+1) / w(i),
    %     L(:,i)   += gamma L(:,i+1),
    %     U(i+1,:) -= gamma U(i,:),
    %     w(i+1)   := 0.
    gamma = w(i+1) / w(i);
    A(i+1,    i) = A(i+1,    i) + gamma;
    A(i+2:end,i) = A(i+2:end,i) + gamma*A(i+2:end,i+1);
    uSub(i) = -gamma*A(i,i);
    A(i+1,i+1:end) = A(i+1,i+1:end) - gamma*A(i,i+1:end);
  end
end

% Add the modified w v' into U
A(1,:) = A(1,:) + w(1)*v';

% Transform U from upper-Hessenberg to upper-triangular form
for i=1:minDim-1,
  % Decide if we should pivot the i'th and i+1'th rows U
  rightTerm = abs(A(i+1,i)*A(i,i)+uSub(i));
  if abs(A(i,i)) < tau*rightTerm,
    % P := P_i P
    p([i,i+1])=p([i+1,i]);

    % Simultaneously perform 
    %   U := P_i U and
    %   L := P_i L P_i^T.
    %
    % Then update
    %     L := L T_{i,L}^{-1},
    %     U := T_{i,L} U, 
    % where T_{i,L} is the Gauss transform which zeros U(i+1,i).
    % 
    % More succinctly,
    %     gamma    := U(i+1,i) / U(i,i),
    %     L(:,i)   += gamma L(:,i+1),
    %     U(i+1,:) -= gamma U(i,:),
    gamma = A(i,i) / uSub(i);
    lambdaSub = A(i+1,i);
    lambda_ii = 1 + gamma*lambdaSub;
    A(i+1,i) = uSub(i);
    A(i,i) = gamma;

    lBiCopy = A(i+2:end,i);
    A(i+2:end,[i,i+1]) = A(i+2:end,[i+1,i]);
    A(i+2:end,i) = A(i+2:end,i) + gamma*lBiCopy;

    uip1RCopy = A(i+1,i+1:end);
    A([i,i+1],:)=A([i+1,i],:);
    A(i+1,i+1:end) = A(i+1,i+1:end) - gamma*uip1RCopy;

    % Force L back to *unit* lower-triangular form via the transform
    %     L := L T_{i,U}^{-1} D^{-1}, 
    % where D is diagonal and responsible for forcing L(i,i) and 
    % L(i+1,i+1) back to 1. The effect on L is:
    %     eta       := L(i,i+1)/L(i,i),
    %     L(:,i+1)  -= eta L(:,i),
    %     delta_i   := L(i,i),
    %     delta_ip1 := L(i+1,i+1),
    %     L(:,i)   /= delta_i,
    %     L(:,i+1) /= delta_ip1,
    % while the effect on U is
    %     U(i,:)   += eta U(i+1,:)
    %     U(i,:)   *= delta(i),
    %     U(i+1,:) *= delta(i+1).
    eta = lambdaSub/lambda_ii;
    delta_i = lambda_ii;
    delta_ip1 = 1 - eta*gamma;
    A(i+2:end,i+1) = (A(i+2:end,i+1) - eta*A(i+2:end,i)) / delta_ip1;
    A(i+1:end,i  ) =  A(i+1:end,i  )                     / delta_i;

    A(i,i) = A(i,i)*delta_i;
    A(i,  i+1:end) = (A(i,  i+1:end) + eta*A(i+1,i+1:end))*delta_i;
    A(i+1,i+1:end) =  A(i+1,i+1:end)                      *delta_ip1;
  else
    % Update
    %     L := L T_{i,L}^{-1},
    %     U := T_{i,L} U, 
    % where T_{i,L} is the Gauss transform which zeros U(i+1,i).
    %
    % More succinctly,
    %     gamma    := U(i+1,i) / U(i,i),
    %     L(:,i)   += gamma L(:,i+1),
    %     U(i+1,:) -= gamma U(i,:).
    gamma = uSub(i) / A(i,i);
    A(i+1,    i      ) = A(i+1,    i      ) + gamma;
    A(i+2:end,i      ) = A(i+2:end,i      ) + gamma*A(i+2:end,i+1    );
    A(i+1,    i+1:end) = A(i+1,    i+1:end) - gamma*A(i,      i+1:end);
  end
end
