function [x,D]=cheb(n,a,b)
%CHEB Ascending Chebyshev-Lobatto nodes and physical first derivative.
validateattributes(n,{'numeric'},{'scalar','integer','>=',2});
z=-cos(pi*(0:n)'/n);w=(-1).^(0:n)';w([1,end])=w([1,end])/2;
delta=z-z';D=(w'./w)./(delta+eye(n+1));
D(1:n+2:end)=0;D=D-diag(sum(D,2));
x=(a+b)/2+(b-a)*z/2;D=D*(2/(b-a));
x([1,end])=[a;b];
end
