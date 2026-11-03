function q=sourcePolynomial(coeff,x,y,R)
% Physical q represented in normalized Cartesian coordinates; no extension.
x=x/R;y=y/R;q=zeros(size(x));
for b=size(coeff,2):-1:1
    px=zeros(size(x));for a=size(coeff,1):-1:1,px=px.*x+coeff(a,b);end
    q=q.*y+px;
end
end
