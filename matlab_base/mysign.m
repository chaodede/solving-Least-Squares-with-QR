function [s] = mysign(x)
%
%  Returns the sign of a number +1 if x >= 0 and -1 if x < 0.

if x >=0
  s = 1;
else
  s= -1;
end

