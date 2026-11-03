function result=greenPolynomial(coeff,x,y,R)
%GREEN_POLYNOMIAL Exact half-plane Green field of a finite-half-disk source.
% coeff(a+1,b+1) multiplies (x/R)^a*(y/R)^b for 0<=a<=7, 0<=b<=3.
% Source is this polynomial in 0<=y, x^2+y^2<=R^2 and zero outside.
% Evaluation is on the CLOSED upper half disk.  The solution satisfies
% -Delta(psi)=q, psi=0 on y=0, and matches the decaying exterior Green field.
% result contains psi, u1=-psi_y, u2=psi_x, q. No spatial quadrature is used.
validateattributes(coeff,{'numeric'},{'2d','real'});
validateattributes(R,{'numeric'},{'scalar','real','finite','positive'});
assert(size(coeff,1)<=8&&size(coeff,2)<=4,'greenPolynomial:Degree','Maximum degrees are (7,3).');
z=complex(x/R,y/R);xi=real(z);eta=imag(z);
assert(all(isfinite(z),'all')&&all(eta>=-64*eps,'all')&& ...
    all(abs(z)<=1+128*eps,'all'),'greenPolynomial:Domain','Evaluate on the closed upper half disk.');
% Roundoff outside the requested closed disk is projected only radially.
outside=abs(z)>1;z(outside)=z(outside)./abs(z(outside));xi=real(z);eta=imag(z);
P=zeros(size(z));Px=P;Py=P;q=P;H=complex(P);Hp=H;
for aa=1:size(coeff,1)
    for bb=1:size(coeff,2)
        amplitude=coeff(aa,bb);if amplitude==0,continue;end
        a=aa-1;b=bb-1;d=a+b;m=d+2;
        q=q+amplitude*xi.^a.*eta.^b;
        c=-1/((b+1)*(b+2));particularTopMode=0;
        for j=0:floor(a/2)
            ax=a-2*j;by=b+2+2*j;
            if j>0,c=-c*(ax+2)*(ax+1)/(by*(by-1));end
            P=P+amplitude*c*xi.^ax.*eta.^by;
            if ax>0,Px=Px+amplitude*c*ax*xi.^(ax-1).*eta.^by;end
            Py=Py+amplitude*c*by*xi.^ax.*eta.^(by-1);
            if mod(b,2)==1,particularTopMode=particularTopMode-2*imag(c/(2^m*1i^by));end
        end
        angular=angular_coefficients(a,b);
        if mod(b,2)==1
            % P has a free degree-m sine harmonic, although q has degree
            % d=m-2. Remove it before applying the exterior DtN correction.
            H=H-amplitude*particularTopMode*z.^m;
            Hp=Hp-amplitude*particularTopMode*m*z.^(m-1);
            for k=1:d
                beta=-2*imag(angular(k+d+1));if beta==0,continue;end
                c=-amplitude*beta/(2*k*(k-m));
                H=H+c*z.^k;Hp=Hp+c*k*z.^(k-1);
            end
        else
            parity=1+mod(d,2); % odd n for d even, even n for d odd.
            for k=0:d
                c=real(angular(k+d+1));if k>0,c=2*c;end
                if c==0,continue;end
                [F,Fp]=even_source_sum(z,k,m,parity);
                c=-2*amplitude*c/pi;H=H+c*F;Hp=Hp+c*Fp;
            end
        end
    end
end
result=struct('psi',R^2*(P+imag(H)),'u1',-R*(Py+real(Hp)), ...
    'u2',R*(Px+imag(Hp)),'q',q);
end

function f=angular_coefficients(a,b)
% cos(theta)^a sin(theta)^b = sum_{k=-d}^d f_k exp(i k theta).
d=a+b;f=complex(zeros(1,2*d+1));
for j=0:a
    for l=0:b
        k=d-2*(j+l);
        f(k+d+1)=f(k+d+1)+nchoosek(a,j)*nchoosek(b,l)*(-1)^l/(2^d*1i^b);
    end
end
end

function [F,Fp]=even_source_sum(z,k,m,n0)
% Sum z^n/((n^2-k^2)*(n-m)); n=n0,n0+2,... . Its pole indices
% have the opposite parity, so there are no resonant denominators.
F=complex(zeros(size(z)));Fp=F;small=abs(z)<=.65;
if any(small,'all')
    n=n0:2:129;v=1./((n.^2-k*k).*(n-m));w=z(small);
    [p,dp]=ascending_polynomial(v,w.^2);
    F(small)=w.^n0.*p;Fp(small)=w.^(n0-1).*(n0*p+2*w.^2.*dp);
end
if ~any(~small,'all'),return;end
w=z(~small);t=(1-w).*(1+w);A1=complex(zeros(size(w)));
interior=w~=1&w~=-1;A1(interior)=t(interior).*atanh(w(interior));
if k==0
    % F=(z^m-1)*atanh(z)/m^2 - oddLi2(z)/m + N_m/m^2.
    % H' uses the already-factored logarithm coefficient, finite at +/-1.
    [Q,~]=ascending_polynomial(ones(1,m/2),w.^2);
    T=.5*(ipm_cr.poisson.fem.complexDilog(w)-ipm_cr.poisson.fem.complexDilog(-w));
    [N,Np]=S_correction(m,n0,w);
    F(~small)=-Q.*A1/(m*m)-T/m+N/(m*m);
    Fp(~small)=-Q.*A1./(m*w)-Q/(m*m)+Np/(m*m);
else
    A=1/(2*k*(k-m));B=1/(2*k*(k+m));C=1/(m*m-k*k);
    % z^k*M(z)=A*z^(2k)+B+C*z^(m+k) has a double root at z^2=1.
    p=zeros(1,(m+k)/2+1);p(1)=B;p(k+1)=p(k+1)+A;p(end)=p(end)+C;
    p=divide_one_minus_w(divide_one_minus_w(p));
    [Q,Qp]=ascending_polynomial(p,w.^2);A2=t.*A1;
    [N1,D1]=S_correction(k,n0,w);[N2,D2]=S_correction(-k,n0,w);[N3,D3]=S_correction(m,n0,w);
    N=A*N1+B*N2+C*N3;Np=A*D1+B*D2+C*D3;
    F(~small)=w.^(-k).*Q.*A2+N;
    Fp(~small)=w.^(-k).*((-k./w).*Q.*A2-4*w.*Q.*A1+2*w.*Qp.*A2+t.*Q)+Np;
end
end

function q=divide_one_minus_w(p)
% Ascending coefficients. Exact algebra gives zero discarded remainder.
% Its floating-point evaluation is qualified by independent Green tests.
q=cumsum(p(1:end-1));
end

function [N,Np]=S_correction(a,n0,z)
% S_a=z^a atanh(z)+N_a, summing n=n0,n0+2,... .
N=complex(zeros(size(z)));Np=N;
if a>=n0+1
    for n=n0:2:a-1
        c=1/(n-a);N=N+c*z.^n;Np=Np+c*n*z.^(n-1);
    end
else
    for j=1:2:n0-a-2
        n=a+j;c=-1/j;N=N+c*z.^n;
        if n~=0,Np=Np+c*n*z.^(n-1);end
    end
end
end

function [p,dp]=ascending_polynomial(c,z)
p=complex(zeros(size(z)))+c(end);dp=complex(zeros(size(z)));
for j=numel(c)-1:-1:1,dp=dp.*z+p;p=p.*z+c(j);end
end
