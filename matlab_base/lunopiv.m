function [A] = lunopiv(A)

%    [A] = lunopiv (A)
%
%
%
%      This function computes the LU factorization with no pivoting.
%
%    Copyright 1994 by Carlos F. Borges. All rights reserved.



% Initialize

[n m] = size(A);

if (n > m)
  dim = m-1;
else
  dim = n-1;
end;

for i=1:dim

  A(i+1:n,i) =  A(i+1:n,i)/A(i,i);
  A(i+1:n,i+1:m) = A(i+1:n,i+1:m) - A(i+1:n,i) * A(i,i+1:m);

end

if (n > m)
  A(m+1:n,m) = A(m+1:n,m) / A(m,m);
end 



