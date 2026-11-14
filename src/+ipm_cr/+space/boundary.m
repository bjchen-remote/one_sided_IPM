classdef boundary
    %BOUNDARY Angular DtN modes, their radial factors and optional verification load.
    methods (Static)
        function p=prepare(p)
            p.M=p.thetaAxis.M(p.jj,p.jj);
            [p.V,p.mu,p.eigenDefect]=mass_eig(p.thetaAxis.K(p.jj,p.jj),p.M);
            p.radialFree=2:p.nr;
            m=numel(p.jj);
            p.operators=cell(1,m);p.factors=cell(1,m);p.harmonic=zeros(p.nr,m);
            for j=1:m
                A=p.tAxis.K+p.mu(j)*p.tAxis.M;
                A(end,end)=A(end,end)+sqrt(p.mu(j));
                p.operators{j}=A;p.factors{j}=decomposition(A(p.radialFree,p.radialFree),'chol');
                h=zeros(p.nr,1);h(1)=1;
                h(p.radialFree)=p.factors{j}\(-A(p.radialFree,1));p.harmonic(:,j)=h;
            end
            p.A=kron(p.thetaAxis.M,p.tAxis.K)+kron(p.thetaAxis.K,p.tAxis.M);
            boundaryMatrix=zeros(p.nt);boundaryMatrix(p.jj,p.jj)=p.M*p.V*diag(sqrt(p.mu))*p.V.'*p.M;
            endRow=sparse(p.nr,p.nr);endRow(end,end)=1;
            p.A=p.A+kron(sparse(boundaryMatrix),endRow);
        end

        function load=load(p,boundary)
            % An empty load imposes homogeneous DtN. Independent tests may prescribe
            % a complete Robin datum at the outer circle.
            load=zeros(p.nt,1);
            if isempty(boundary),return;end
            assert(isstruct(boundary)&&isscalar(boundary)&&isfield(boundary,'value'), ...
                'ipm_cr:FEMExteriorData','Supply the complete outer Robin datum.');
            if isa(boundary.value,'function_handle')
                theta=p.thetaAxis.qx;R=p.r(end);
                q=sample_callback(boundary.value,R*cos(theta),R*sin(theta));
            else
                values=boundary.value(:);
                assert(numel(values)==p.nt&&isreal(values),'ipm_cr:FEMExteriorData');
                q=p.thetaAxis.Bq*(p.thetaAxis.interpolation\values);
            end
            load=p.thetaAxis.Bq'*(p.thetaAxis.qw.*q(:));
        end
    end
end

function [V,lambda,defect]=mass_eig(K,M)
K=full((K+K')/2);M=full((M+M')/2);L=chol(M,'lower');
H=L\K/L';H=(H+H')/2;[Q,D]=eig(H);[lambda,order]=sort(diag(D));V=L'\Q(:,order);
for j=1:size(V,2),[~,i]=max(abs(V(:,j)));if V(i,j)<0,V(:,j)=-V(:,j);end,end
defect=max(norm(V'*M*V-eye(numel(lambda)),inf), ...
    norm(K*V-M*V*diag(lambda),inf)/max(1,norm(K*V,inf)));
end

function q=sample_callback(callback,x,y)
% Evaluate a real-valued callback on a FEM query array.
q=callback(x,y);if isscalar(q),q=repmat(q,size(x));end
assert(isnumeric(q)&&isreal(q)&&isequal(size(q),size(x)), ...
    'ipm_cr:FEMSource','Callbacks must return real numeric values matching the supplied query shape.');
end
