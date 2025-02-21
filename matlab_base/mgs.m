function [Q, R] = mgs(A)

%   [Q R] = mgs(A)
%
%    This function computes the QR factorization of a general input
%  matrix A using the modified Gram-Schmidt algorithm. Q is an
%  isometry and R is upper-triangular. R is built up one row at
%  a time from left to right.
%
%    Copyright 1995 by Carlos F. Borges. All rights reserved.


% Initialize
[n m] = size(A);
Q = zeros(n,m);
R = zeros(m,m);

qhat = A(:,1);
R(1,1) = norm(qhat);
qhat = qhat/R(1,1);
Q(:,1) = qhat;

% Loop across the columns.
for j=2:m

  % Compute the projection of the current vector onto the remaining
  % columns of A. 
  % for k=j:m
  %   R(j-1,k) = qhat'*A(:,k);
  % end
  % To do this without a loop you use the following line instead.
  R(j-1,j:m) = qhat'*A(:,j:m);

  % Subtract out the projections so that the remaining columns of A
  % will be orthogonal to all of the columns of Q.
  % for k=j:m
  %   A(:,k) = A(:,k) - R(j-1,k)*qhat;
  % end
  % To do this without a loop you use the following line instead.
  A(:,j:m) = A(:,j:m) - qhat*R(j-1,j:m);

  % Get the next vector and add it to Q.
  qhat = A(:,j);
  R(j,j) = norm(qhat);
  qhat = qhat/R(j,j);
  Q(:,j) = qhat;

end
