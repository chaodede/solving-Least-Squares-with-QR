function [x] = upper_solve (U, b)

%    x = upper_solve (U, b)
%
%      This function solves linear systems, Ux = b, where the matrix U is
%    upper triangular (U(i,j) = 0 if j < i).  The code fails if U is
%    singular (there will be a divide by zero).
%
%    Copyright 1994 by Carlos F. Borges. All rights reserved.



% Initialize

[n, m] = size(U);

% Check that U is square.

if (n ~= m)
  error('U is not square.');
end


% Check that b is an nx1 column vector.

if (norm(size(b) - [n 1]) ~= 0)
  error('U and b are not conformably sized');
end

% We have a valid linear system. Begin solution process.

% Allocate the solution vector.

x = zeros(n,1);

for k=n:-1:2

  x(k) = b(k)/U(k,k);
  b(1:k-1) = b(1:k-1) - x(k) * U(1:k-1,k);

end

x(1) = b(1)/U(1,1);
