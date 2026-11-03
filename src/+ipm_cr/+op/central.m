function [plus,minus]=central(field,spacing,dissipation)
%CENTRAL Fourth-order centered derivative with explicit sixth-order damping.
% Each row is a logical coordinate line. With LF speed alpha, the existing
% pair contract gives -v*D4(field)-alpha*Q6(field). The positive periodic Q6
% has symbol dissipation*sin(k*h/2)^6/h and vanishes through degree five.
% Physical derivatives follow by dividing by the actual mapping Jacobian.
% Wall/exterior closures use the same one-sided quintic extension in both
% directions; their stability must be checked in the coupled problem.
validateattributes(field,{'double'},{'2d','real'});
validateattributes(spacing,{'numeric'},{'scalar','real','finite','positive'});
if nargin<3,dissipation=.2;end
validateattributes(dissipation,{'numeric'},{'scalar','real','finite','>=',0,'<=',1});
n=size(field,2);assert(n>=7,'ipm_cr:CentralWidth','Seven samples are required.');
weights=[6,-15,20,-15,6,-1;21,-70,105,-84,35,-6;56,-210,336,-280,120,-21]';
left=field(:,1);right=field(:,end);
lo=(field(:,1:6)-left)*weights;
hi=(field(:,end:-1:end-5)-right)*weights;
dm1=shift(-1);dm2=shift(-2);dp1=shift(1);dp2=shift(2);
D=(dm2-8*dm1+8*dp1-dp2)/(12*spacing);
if dissipation==0,Q=zeros(size(D));
else,Q=dissipation*(-15*(dm1+dp1)+6*(dm2+dp2)-shift(-3)-shift(3))/(64*spacing);end
plus=D+Q;minus=D-Q;

    function value=shift(k)
        % Form every difference about its native center. Exterior offsets
        % are not added to a large boundary value before subtraction.
        m=abs(k);
        if k<0
            value=[left-field(:,1:m)+lo(:,m:-1:1),field(:,1:n-m)-field(:,1+m:n)];
        else
            value=[field(:,1+m:n)-field(:,1:n-m),right-field(:,n-m+1:n)+hi(:,1:m)];
        end
    end
end
