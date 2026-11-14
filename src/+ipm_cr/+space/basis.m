classdef basis
    %BASIS Quintic spline axes and spatial derivatives on native sites.
    methods (Static)
        function a=build(nodes)
            order=6;
            nodes=nodes(:);n=numel(nodes);p=5;
            assert(n>p&&all(diff(nodes)>0),'ipm_cr:SplineNodes','Use at least six increasing native nodes.');
            % Averaged knots satisfy Schoenberg-Whitney at the original sample sites.
            inner=zeros(n-p-1,1);
            for k=1:p,inner=inner+nodes((1:numel(inner))+k)/p;end
            t=[repmat(nodes(1),p+1,1);inner;repmat(nodes(end),p+1,1)];
            breaks=[nodes(1);inner;nodes(end)];
            k=(1:order-1)';off=k./sqrt(4*k.^2-1);
            [V,L]=eig(diag(off,1)+diag(off,-1));
            [x,idx]=sort(diag(L));w=2*V(1,idx).^2;
            h=diff(breaks);qx=reshape((breaks(1:end-1)+h.*((x'+1)/2)).',[],1);
            qw=reshape((h*w/2).',[],1);
            a=struct('nodes',nodes,'breaks',breaks,'knots',t,'degree',p, ...
                'elementDegrees',p*ones(numel(h),1),'partitionConditioning',struct(), ...
                'boundaryDegreePolicy',struct('strategy','clamped_quintic_spline'), ...
                'qx',qx,'qw',qw,'quadratureOrder',order);
            [a.Bq,a.Dq]=ipm_cr.space.basis.evaluate(a,qx);
            [a.E,a.D]=ipm_cr.space.basis.evaluate(a,nodes);
            a.M=a.Bq'*(qw.*a.Bq);a.K=a.Dq'*(qw.*a.Dq);
            a.M=(a.M+a.M')/2;a.K=(a.K+a.K')/2;
            a.interpolation=decomposition(a.E,'lu');
            a.mass=decomposition(a.M,'chol');
            a.convection=a.Bq'*(qw.*a.Dq);
            % Penalize the jumps of the fifth derivative. Lower derivatives remain
            % continuous, including at the open boundaries. The penalty is PSD.
            % The medium-grid transport spectrum requires weaker jumps for RK4 at
            % the advective CFL; this coefficient scales the same high-order form.
            nb=numel(inner);rows=(1:nb)';left=p+rows;
            J=sparse([rows;rows],[left+1;left],[ones(nb,1);-ones(nb,1)],nb,numel(t)-1);
            for j=1:p,J=J*derivative(t,j);end
            spans=diff([nodes(1);inner;nodes(end)]);
            a.penalty=J'*(min(spans(1:end-1),spans(2:end)).^(2*p).*J)/65536;
        end

        function [B,D]=evaluate(a,queries)
            q=queries(:);t=a.knots;p=a.degree;n=numel(t);
            assert(all(q>=t(1)&q<=t(end)),'ipm_cr:SplineQuery','Query is outside its spline chart.');
            [knots,start]=unique(t,'last');span=start(discretize(q,knots));
            span(q==t(end))=start(end-1);
            i=(1:numel(q))';
            B=sparse(i,span,1,numel(q),n-1);
            for k=1:p
                if k==p,D=B*derivative(t,k);end
                [row,j,value]=find(B);
                row=row(:);j=j(:);value=value(:);
                ll=t(j+k)-t(j);
                keep=j<=n-k-1 & ll>0;
                A=sparse(row(keep),j(keep),value(keep).*(q(row(keep))-t(j(keep)))./ll(keep),numel(q),n-k-1);
                % Contribution of old basis j to the new basis j-1.
                keep=j>1 & (t(j+k)-t(j))>0;
                B=A+sparse(row(keep),j(keep)-1,value(keep).*(t(j(keep)+k)-q(row(keep)))./ll(keep), ...
                    numel(q),n-k-1);
            end
        end

        function [center,penalty]=differentiate(a,field,weak)
            % Columns are independent lines. Retain constants exactly by translation.
            f=field-field(1,:);
            coeff=a.interpolation\f;
            penalty=[];
            if weak
                Z=a.mass\[a.convection*coeff,a.penalty*coeff];
                Z=a.E*Z;k=size(field,2);center=Z(:,1:k);penalty=Z(:,k+1:end);
            else
                center=a.D*coeff;
            end
        end

        function d=gradient(rho,g,spline)
            %DIFF Shared spline derivative, annihilating constants first.
            d=struct();
            D=@(axis,field)ipm_cr.space.basis.differentiate(axis,field,false);
            d.core=struct('x',D(spline.x,rho.core.').','y',D(spline.y,rho.core));
            fr=D(spline.t,rho.polar)./g.polar.r;ft=D(spline.theta,rho.polar.').';
            d.polar=struct('x',fr.*g.polar.cosTheta-ft.*g.polar.sinTheta./g.polar.r, ...
                    'y',fr.*g.polar.sinTheta+ft.*g.polar.cosTheta./g.polar.r);
        end
    end
end

function G=derivative(t,p)
j=(1:numel(t)-p-1)';
l=t(j+p)-t(j);r=t(j+p+1)-t(j+1);
L=l>0;R=r>0;
G=sparse([j(L);j(R)+1],[j(L);j(R)],[p./l(L);-p./r(R)],numel(j)+1,numel(j));
end
