function value=complexDilog(z)
%COMPLEX_DILOG Principal Li_2 on the closed unit disk, without toolboxes.
% Transformations reduce every power series to |z| <= 3/4.  In particular,
% values near +/-1 do not rely on a slowly converging boundary power series.
assert(all(isfinite(z),'all')&&all(abs(z)<=1+64*eps,'all'), ...
    'greenPolynomial:DilogDomain','Li2 input must lie in the closed unit disk.');
value=zeros(size(z),'like',complex(z));
one=z==1;minus=z==-1;value(one)=pi^2/6;value(minus)=-pi^2/12;
remaining=~(one|minus);
small=remaining&abs(z)<=.75;
if any(small,'all'),value(small)=power_series(z(small));end
negative=remaining&~small&real(z)<=0;
if any(negative,'all')
    w=z(negative);value(negative)=-ipm_cr.poisson.fem.complexDilog(w./(w-1))-.5*log(1-w).^2;
end
right=remaining&~small&~negative&real(z)>.5;
if any(right,'all')
    w=z(right);value(right)=pi^2/6-log(w).*log(1-w)-ipm_cr.poisson.fem.complexDilog(1-w);
end
middle=remaining&~small&~negative&~right;
if any(middle,'all')
    w=z(middle);value(middle)=.5*ipm_cr.poisson.fem.complexDilog(w.^2)-ipm_cr.poisson.fem.complexDilog(-w);
end
end

function value=power_series(z)
value=zeros(size(z),'like',z);correction=value;power=z;
for k=1:128
    term=power/(k*k)-correction;
    next=value+term;correction=(next-value)-term;value=next;power=power.*z;
end
end
