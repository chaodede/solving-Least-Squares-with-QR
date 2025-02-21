function [v] = get_house(x,i,j)
%
%    This function computes a Householder vector that zeros out the
%  i+1'st to j'th elements of a vector x and puts that energy into the
%  i'th element.  It prenormalizes the vector so that the form of a
%  Householder reflection can be reduced to I - v*v'.
%
%  Copyright 1995 by Carlos F. Borges. All rights reserved.

% Initialization.
n = length(x);
v = zeros(n,1);

% Copy that part of x to be worked on to the corresponding positions in v.
v(i:j) = x(i:j);

% Compute the proper Householder vector.
v(i) = v(i) + mysign(x(i)) * norm(x(i:j));

% Normalize the result so that H = I - v*v'. Includes an error check for 
% the trivial reflection.

if ((v'*v) > 0)
  v = v*sqrt(2/(v'*v));
end

