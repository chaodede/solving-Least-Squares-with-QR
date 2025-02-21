function [Q, R] = cgs(A)

%   [Q R] = cgs(A)
%
%    This function computes the QR factorization of a general input
%  matrix A using the classical Gram-Schmidt algorithm. Q is an
%  isometry and R is upper-triangular. R is built up one column at
%  a time from left to right.
%
%    Copyright 1995 by Carlos F. Borges. All rights reserved.


% Initialize
[n, m] = size(A);
Q = zeros(n,m);
R = zeros(m,m);

% Loop across the columns.
for j=1:m

  % Grab the next column of A.
  qhat = A(:,j);

  % Compute the projection of the current vector onto the previous
  % columns of Q. 
  % for i=1:j-1
  %   R(i,j) = qhat'*Q(:,i);
  % end
  % To do this without a loop you use the following line instead.
  R(:,j) = Q'*qhat;

  % Subtract out the projections so that qhat will be orthogonal to
  % the previous vectors.
  % for i=1:j-1
  %   qhat = qhat - R(i,j)*Q(:,i);
  % end
  % To do this without a loop you use the following line instead.
  qhat = qhat - Q*R(:,j);
  % NOTICE: Although this single line is mathematically equivalent
  % to the loop above, it is not computationally equivalent. Why?
  % You may notice that using this single line reduces the size of
  % Q*R - A but does not noticeably change the size of Q'Q - I.
  % Can you explain this phenomenon?

  % Normalize qhat and append it to Q.
  R(j,j) = norm(qhat);
  Q(:,j) = qhat/R(j,j);

end
