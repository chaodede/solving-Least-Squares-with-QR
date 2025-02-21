%
% myldlt - computes LDL^T factorization of
%          symmetric matrix A
%
function [L,D]=myldlt(A)
[m,n]=size(A);
L=zeros(n);
v=zeros(n,1);
d=zeros(1,n);
for j=1:n
    Lj=A(j:n,j);
    for k=1:j-1
        v(k)=d(k)*L(j,k);
        Lj=Lj-L(j:n,k)*v(k);
    end
    d(j)=Lj(1);
    L(j:n,j)=Lj/d(j);
end
D=diag(d);
