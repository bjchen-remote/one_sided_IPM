classdef solver < handle
    % Experimental same-transport-node overlapping Poisson backend.
    % Seven-point physical-axis D2: >=fifth local consistency on regular
    % nonuniform families; complete psi AND velocity target is fourth order.
    properties (SetAccess=private)
        g; c; p; coupling; sourceMaps; auditMaps; buildSeconds; config
        Counters=struct('builds',1,'solves',0,'factorReuses',0, ...
            'accepted',0,'rejected',0,'unaudited',0,'directFallbacks',0)
    end
    methods
        function obj=solver(g,options)
            if nargin<2,options=struct();end
            timer=tic;
            assert(isstruct(g)&&isscalar(g)&&all(isfield(g,{'core','polar'})), ...
                'ipm_cr:FDGeometry','Use a complete core/annulus transport geometry.');
            assert(isstruct(options)&&isscalar(options)&& ...
                all(ismember(fieldnames(options),{'width','traceWidth','outerBoundary'})), ...
                'ipm_cr:FDOptions','Unknown same-mesh FD option.');
            obj.g=g;
            cfg=struct('width',7,'traceWidth',7,'outerBoundary','decay');
            names=fieldnames(options);
            for k=1:numel(names),cfg.(names{k})=options.(names{k});end
            for name={'width','traceWidth'}
                value=cfg.(name{1});
                assert(isnumeric(value)&&isreal(value)&&isscalar(value)&&isfinite(value)&&value==7, ...
                    'ipm_cr:FDOptions','The verified FD configuration uses width=traceWidth=7.');
            end
            assert(((ischar(cfg.outerBoundary)&&isrow(cfg.outerBoundary))|| ...
                (isstring(cfg.outerBoundary)&&isscalar(cfg.outerBoundary)))&& ...
                strcmp(cfg.outerBoundary,'decay'),'ipm_cr:FDOptions', ...
                'The same-mesh backend retains source-free exterior DtN (decay).');
            cfg.outerBoundary=char(cfg.outerBoundary);
            obj.config=cfg;w=cfg.width;
            c=g.core;c.Dx=local_fd(c.x,1,w);c.Dy=local_fd(c.y,1,w);
            c.Dxx=local_fd(c.x,2,w);c.Dyy=local_fd(c.y,2,w);
            c.ix=2:c.nx-1;c.iy=2:c.ny-1;
            mask=false(size(c.X));mask(:,[1,end])=true;mask(end,:)=true;mask(1,:)=false;
            c.traceIndex=find(mask);
            [c.Vx,c.lx,c.condX]=positive_eig(-c.Dxx(c.ix,c.ix));
            [c.Vy,c.ly,c.condY]=positive_eig(-c.Dyy(c.iy,c.iy));
            c.Vxi=c.Vx\eye(numel(c.ix));c.Vyi=c.Vy\eye(numel(c.iy));
            c.denominator=c.ly+c.lx.';
            p=g.polar;p.t=log(p.r);p.jj=2:p.nt-1;
            p.Dt=local_fd(p.t,1,w);p.Dtt=local_fd(p.t,2,w);
            p.Dtheta=local_fd(p.theta,1,w);p.Dtheta2=local_fd(p.theta,2,w);
            [p.V,p.mu,p.eigenCondition]=positive_eig(-p.Dtheta2(p.jj,p.jj));
            p.Vinv=p.V\eye(numel(p.jj));m=numel(p.jj);
            p.A=cell(1,m);p.factor=cell(1,m);p.scale=zeros(p.nr,m);p.harmonic=zeros(p.nr,m);
            for j=1:m
                A=-p.Dtt+p.mu(j)*speye(p.nr);A(1,:)=0;A(1,1)=1;
                if strcmp(cfg.outerBoundary,'decay')
                    A(end,:)=p.Dt(end,:);A(end,end)=A(end,end)+sqrt(p.mu(j));
                else
                    A(end,:)=0;A(end,end)=1;
                end
                sc=1./full(max(abs(A),[],2));
                p.A{j}=A;p.scale(:,j)=sc;p.factor{j}=decomposition(spdiags(sc,0,p.nr,p.nr)*A,'lu');
                b=zeros(p.nr,1);b(1)=1;p.harmonic(:,j)=refined_modal(p,b,j);
            end
            obj.c=c;obj.p=p;
            RCA=trace_map(g,'polar',c.X(c.traceIndex),c.Y(c.traceIndex),cfg.traceWidth);
            RAC=trace_map(g,'core',p.X(1,p.jj)',p.Y(1,p.jj)',cfg.traceWidth);
            response=zeros(numel(c.X),m);T=zeros(m);
            for j=1:m
                lam=zeros(1,m);lam(j)=1;
                pa=obj.polar_harmonic(lam);
                cc=obj.core_solve(zeros(size(c.X)),RCA*pa(:));
                response(:,j)=cc(:);T(:,j)=RAC*cc(:);
            end
            S=eye(m)-T;sv=svd(S);
            assert(sv(end)>1e-11,'same_mesh:Schur','Singular or poorly conditioned trace system.');
            obj.coupling=struct('RCA',RCA,'RAC',RAC,'S',S,'factor',decomposition(S,'lu'), ...
                'coreResponse',response,'sigmaMin',sv(end),'condition',sv(1)/sv(end), ...
                'mapNorm',max(norm(RCA,inf),norm(RAC,inf)), ...
                'responseInfinity',norm(response,inf));
            z=struct('core',zeros(size(c.X)),'polar',zeros(size(p.X)));
            [~,mc]=ipm_cr.mesh.sample(z,g,c.X,c.Y);[~,ma]=ipm_cr.mesh.sample(z,g,p.X,p.Y);
            obj.sourceMaps=struct('core',mc,'polar',ma);
            % Non-node query cloud independently reconstructs psi derivatives.
            hi=min([g.core.activeXBounds(2),g.core.activeYBounds(2)]);
            lo=g.polar.innerActiveRadius;
            assert(hi>lo,'same_mesh:Overlap','No independent active overlap query region.');
            rq=lo+(hi-lo)*[.17;.48;.81];th=((1:31)+.21)*pi/32.5;
            [tt,rr]=meshgrid(th,rq);x=rr.*cos(tt);y=rr.*sin(tt);
            obj.auditMaps=struct('x',x,'y',y,'core',query_map(c.x,c.y,x,y,'core',w), ...
                'polar',query_map(p.t,p.theta,x,y,'polar',w));
            obj.buildSeconds=toc(timer);
        end
        function [flow,info]=solve(obj,source,doAudit,outer,analyticSource)
            if nargin<3,doAudit=true;end
            if nargin<4,outer=[];end
            if nargin<5,analyticSource=[];end
            validateattributes(doAudit,{'logical'},{'scalar'});
            assert(isempty(analyticSource)||isa(analyticSource,'function_handle'), ...
                'ipm_cr:FDSource','An analytic source offset must be a function handle.');
            total=tic;phase=tic;c=obj.c;p=obj.p;cp=obj.coupling;
            if isa(source,'function_handle')
                qc=sample_source(source,c.X,c.Y);qa=sample_source(source,p.X,p.Y);mode='analytic';
            else
                validate_native(source,obj.g);
                qc=source_map(obj.sourceMaps.core,source);qa=source_map(obj.sourceMaps.polar,source);mode='native_unique_owner';
            end
            if ~isempty(analyticSource)
                qc=qc+sample_source(analyticSource,c.X,c.Y);qa=qa+sample_source(analyticSource,p.X,p.Y);
            end
            assert(all(isfinite(qc),'all')&&all(isfinite(qa),'all'), ...
                'ipm_cr:FDSource','The complete Poisson source must remain finite.');
            obj.Counters.solves=obj.Counters.solves+1;
            obj.Counters.factorReuses=obj.Counters.factorReuses+1;
            sourceSeconds=toc(phase);phase=tic;
            [pa0,modal0,rhs]=obj.polar_solve(qa,[],outer);
            cc0=obj.core_solve(qc,cp.RCA*pa0(:));
            b=cp.RAC*cc0(:);lambda=cp.factor\b;
            for j=1:2
                defect=b-cp.S*lambda;
                trial=lambda+cp.factor\defect;
                if norm(b-cp.S*trial,inf)>=norm(defect,inf),break;end
                lambda=trial;
            end
            cc=cc0+reshape(cp.coreResponse*lambda,size(cc0));
            lm=lambda.'*p.Vinv.';modal=modal0+p.harmonic.*lm;
            pa=zeros(size(p.X));pa(:,p.jj)=real_checked(modal*p.V.');
            solveSeconds=toc(phase);phase=tic;
            cy=apply_difference(c.Dy,cc);cx=apply_difference(c.Dx,cc')';
            pt=apply_difference(p.Dt,pa);pth=apply_difference(p.Dtheta,pa')';
            u1=-(p.sinTheta.*pt+p.cosTheta.*pth)./p.r;
            u2=(p.cosTheta.*pt-p.sinTheta.*pth)./p.r;
            cx(1,:)=0;u2(:,[1,end])=0;
            flow=struct('core',struct('psi',cc,'u1',-cy,'u2',cx), ...
                'polar',struct('psi',pa,'u1',u1,'u2',u2));
            returnSeconds=toc(phase);
            info=struct('sourceSeconds',sourceSeconds,'solveSeconds',solveSeconds, ...
                'returnSeconds',returnSeconds,'kernelSeconds',sourceSeconds+solveSeconds+returnSeconds, ...
                'sourceMode',mode,'audited',doAudit);
            % Both diagnostic kernels retain original-row algebraic checks.
            res=-c.Dyy*cc-cc*c.Dxx'-qc;
            den=abs(c.Dyy)*abs(cc)+abs(cc)*abs(c.Dxx')+abs(qc);
            backward=max(abs(res(c.iy,c.ix))./max(den(c.iy,c.ix),realmin),[],'all');
            rhs(1,:)=lm;
            for j=1:numel(p.mu)
                d=p.A{j}*modal(:,j)-rhs(:,j);
                denom=abs(p.A{j})*abs(modal(:,j))+abs(rhs(:,j));
                backward=max(backward,max(abs(d)./max(denom,realmin)));
            end
            initialBackward=backward;fallback=false;
            if backward>1e-10
                % Retain the same original equations and gate. Superposition
                % of tiny modes can lose relative digits: use direct patch
                % solves with the already obtained coupled Dirichlet trace.
                fallback=true;
                [pa,modal,rhs]=obj.polar_solve(qa,lambda.',outer);
                cc=obj.core_solve(qc,cp.RCA*pa(:));
                cy=apply_difference(c.Dy,cc);cx=apply_difference(c.Dx,cc')';
                pt=apply_difference(p.Dt,pa);pth=apply_difference(p.Dtheta,pa')';
                u1=-(p.sinTheta.*pt+p.cosTheta.*pth)./p.r;
                u2=(p.cosTheta.*pt-p.sinTheta.*pth)./p.r;cx(1,:)=0;u2(:,[1,end])=0;
                flow.core=struct('psi',cc,'u1',-cy,'u2',cx);flow.polar=struct('psi',pa,'u1',u1,'u2',u2);
                res=-c.Dyy*cc-cc*c.Dxx'-qc;den=abs(c.Dyy)*abs(cc)+abs(cc)*abs(c.Dxx')+abs(qc);
                backward=max(abs(res(c.iy,c.ix))./max(den(c.iy,c.ix),realmin),[],'all');
                for j=1:numel(p.mu)
                    d=p.A{j}*modal(:,j)-rhs(:,j);denom=abs(p.A{j})*abs(modal(:,j))+abs(rhs(:,j));
                    backward=max(backward,max(abs(d)./max(denom,realmin)));
                end
            end
            info.backwardError=backward;info.algebraicAccepted=backward<1e-10;
            info.directFallback=fallback;info.beforeFallbackBackward=initialBackward;
            assert(info.algebraicAccepted,'same_mesh:Backward','Original-row backward error %.4g exceeds 1e-10.',backward);
            if doAudit
                trace=max(norm(cp.RAC*cc(:)-lambda,inf),norm(cp.RCA*pa(:)-cc(c.traceIndex),inf));
                [ac,aa]=obj.overlap(flow);
                info.traceAbsolute=trace;
                info.overlapPsi=max(abs(ac.psi-aa.psi));
                info.overlapVelocity=max([abs(ac.u1-aa.u1);abs(ac.u2-aa.u2)]);
                info.schurCondition=cp.condition;info.schurSigmaMin=cp.sigmaMin;
                info.accepted=backward<1e-10&&trace<1e-9;
                info.auditDoesNotCertifyContinuum=true;
                info.traceGateAccepted=trace<1e-9;
                info.traceResidual=trace;info.overlapValueMismatch=info.overlapPsi;
                info.overlapVelocityMismatch=info.overlapVelocity;
                if info.accepted,obj.Counters.accepted=obj.Counters.accepted+1;
                else,obj.Counters.rejected=obj.Counters.rejected+1;end
            else
                info.accepted=NaN;info.traceGateAccepted=NaN;
                info.traceAbsolute=NaN;info.traceResidual=NaN;
                info.overlapPsi=NaN;info.overlapVelocity=NaN;
                info.overlapValueMismatch=NaN;info.overlapVelocityMismatch=NaN;
                obj.Counters.unaudited=obj.Counters.unaudited+1;
            end
            if fallback,obj.Counters.directFallbacks=obj.Counters.directFallbacks+1;end
            info.acceptanceScope='finite-grid algebraic backward error and, when audited, artificial Dirichlet trace only';
            info.spectralSpatialAuditPerformed=false;info.continuumAccuracyCertified=false;
            info.spatialAccuracyCertified=false;info.auditDoesNotCertifyContinuum=true;
            info.linearBackwardError=backward;info.linearSolveBackwardError=backward;
            info.linearSolveComponentwiseBackwardError=backward;
            info.relativeResidual=backward;info.algebraicResidual=backward;
            info.normalFluxMismatch=NaN;info.strongResidual=NaN;info.strongResidualAvailable=false;
            info.nativeReconstructionDefect=struct('available',false,'scaledMaximum',NaN, ...
                'independentContinuumCertificate',false,'independentSource',false);
            info.sourceInterpolationAccuracyCertified=false;info.outerSourceTailCertified=false;
            info.analyticBackgroundSource=~isempty(analyticSource);
            info.sourceDefinition='unique-owner already differentiated residual plus analytic source offset at the original physical nodes';
            info.pde='-Delta psi=source; u=(-psi_y,psi_x)';
            info.samePhysicalNodes=true;info.sameDiscreteStreamfunctionDerivatives=true;
            info.factorBuildCount=1;info.factorReuseCount=obj.Counters.factorReuses;
            info.buildSeconds=obj.buildSeconds;
            info.totalSeconds=toc(total);
        end
        function [ac,aa]=overlap(obj,flow)
            ac=evaluate_plan(obj.auditMaps.core,flow.core.psi);
            aa=evaluate_plan(obj.auditMaps.polar,flow.polar.psi);
        end
        function v=evaluate(obj,flow,x,y,chart)
            if strcmp(chart,'core'),plan=query_map(obj.c.x,obj.c.y,x,y,'core',obj.config.width);
            else,plan=query_map(obj.p.t,obj.p.theta,x,y,'polar',obj.config.width);end
            v=evaluate_plan(plan,flow.(chart).psi);
        end
        function U=core_solve(obj,q,trace)
            c=obj.c;U=zeros(size(c.X));U(c.traceIndex)=trace;
            lap=-c.Dyy*U-U*c.Dxx';b=q(c.iy,c.ix)-lap(c.iy,c.ix);
            U(c.iy,c.ix)=real_checked(c.Vy*((c.Vyi*b*c.Vxi.')./c.denominator)*c.Vx.');
            for j=1:2
                residual=q+c.Dyy*U+U*c.Dxx';
                denominator=abs(q)+abs(c.Dyy)*abs(U)+abs(U)*abs(c.Dxx');
                if max(abs(residual(c.iy,c.ix))./max(denominator(c.iy,c.ix),realmin),[],'all')<1e-12,break;end
                d=real_checked(c.Vy*((c.Vyi*residual(c.iy,c.ix)*c.Vxi.')./c.denominator)*c.Vx.');
                U(c.iy,c.ix)=U(c.iy,c.ix)+d;
            end
        end
        function [U,modal,rhs]=polar_solve(obj,q,inner,outer)
            p=obj.p;m=numel(p.mu);
            if isempty(inner),inner=zeros(1,m);end
            if isempty(outer),outer=zeros(1,m);end
            if numel(outer)==p.nt,outer=outer(p.jj);end
            inner=reshape(inner,1,[]);outer=reshape(outer,1,[]);
            assert(numel(inner)==m&&numel(outer)==m,'same_mesh:Boundary','Incorrect trace length.');
            rhs=(q(:,p.jj)*p.Vinv.').*p.r.^2;
            rhs(1,:)=inner*p.Vinv.';rhs(end,:)=outer*p.Vinv.';
            modal=zeros(p.nr,m);
            for j=1:m,modal(:,j)=refined_modal(p,rhs(:,j),j);end
            U=zeros(size(p.X));U(:,p.jj)=real_checked(modal*p.V.');
        end
        function U=polar_harmonic(obj,inner)
            p=obj.p;modal=p.harmonic.*(inner*p.Vinv.');
            U=zeros(size(p.X));U(:,p.jj)=real_checked(modal*p.V.');
        end
    end
end

function u=refined_modal(p,b,j)
A=p.A{j};u=p.factor{j}\(p.scale(:,j).*b);u(1)=b(1);
for k=1:3
    r=b-A*u;back=max(abs(r)./max(abs(A)*abs(u)+abs(b),realmin));
    if back<128*eps,break;end
    z=u+p.factor{j}\(p.scale(:,j).*r);z(1)=b(1);
    rb=b-A*z;next=max(abs(rb)./max(abs(A)*abs(z)+abs(b),realmin));
    if next>=back,break;end
    u=z;
end
end

function D=local_fd(z,m,w)
z=z(:);n=numel(z);assert(n>=w);ii=zeros(n*w,1);jj=ii;vv=ii;
for i=1:n
    first=min(max(i-floor(w/2),1),n-w+1);j=first:first+w-1;
    dz=z(j)-z(i);scale=max(abs(dz));a=(dz/scale)';
    M=a.^((0:w-1)');b=zeros(w,1);b(m+1)=factorial(m)/scale^m;
    v=M\b;v(j==i)=v(j==i)-sum(v);
    slots=(i-1)*w+(1:w);ii(slots)=i;jj(slots)=j;vv(slots)=v;
end
D=sparse(ii,jj,vv,n,n);
end

function [V,lambda,kappa]=positive_eig(A)
[V,L]=eig(full(A));lambda=diag(L);
% One-sided high-order closures need not have real spectra. Retain the full
% conjugate eigensystem; positive real part selects a decaying principal
% square root for the semidiscrete exterior, without dropping imaginary data.
[~,k]=sortrows([real(lambda),imag(lambda)],[1,2]);lambda=lambda(k);V=V(:,k);
assert(all(real(lambda)>0),'same_mesh:Spectrum','Spectrum outside the positive half-plane.');
kappa=cond(V);assert(kappa<1e9,'same_mesh:Spectrum','Ill-conditioned local eigenbasis.');
end

function r=real_checked(z)
defect=max(abs(imag(z)),[],'all');scale=max(1,max(abs(z),[],'all'));
assert(defect<1e-9*scale,'same_mesh:ComplexReturn','Conjugate modes failed to reconstruct a real field.');
r=real(z);
end

function E=trace_map(g,chart,x,y,w)
p=g.(chart);indices=p.activeIndices;
if strcmp(chart,'core')
    M=ipm_cr.mesh.tensorMap(p.y(p.activeRows),p.x(p.activeColumns),y,x,w);
else
    M=ipm_cr.mesh.tensorMap(p.r(p.activeRows),p.theta,hypot(x,y),atan2(y,x),w);
end
[i,j,v]=find(M);E=sparse(i,indices(j),v,numel(x),numel(p.X));
end

function q=source_map(m,s)
q=ipm_cr.mesh.sample(s,m);
end

function value=sample_source(callback,x,y)
value=callback(x,y);
if isscalar(value),value=repmat(value,size(x));end
assert(isnumeric(value)&&isreal(value)&&isequal(size(value),size(x))&&all(isfinite(value),'all'), ...
    'ipm_cr:FDSource','Analytic source values must be finite real arrays of the query shape.');
end

function validate_native(source,g)
assert(isstruct(source)&&isscalar(source)&&all(isfield(source,{'core','polar'}))&& ...
    isnumeric(source.core)&&isreal(source.core)&&isequal(size(source.core),size(g.core.X))&& ...
    isnumeric(source.polar)&&isreal(source.polar)&&isequal(size(source.polar),size(g.polar.X))&& ...
    all(isfinite(source.core),'all')&&all(isfinite(source.polar),'all'), ...
    'ipm_cr:FDSource','Supply complete finite real already differentiated core/polar source arrays.');
end

function result=apply_difference(D,U)
[row,col,weight]=find(D);M=sparse(row,1:numel(row),weight,size(D,1),numel(row));
result=M*(U(col,:)-U(row,:));
end

function plan=query_map(a,b,x,y,chart,w)
% Independent local degree-six reconstruction AND its analytical derivative.
shape=size(x);x=x(:);y=y(:);n=numel(x);
if strcmp(chart,'core'),qa=x;qb=y;else,qa=log(hypot(x,y));qb=atan2(y,x);end
a=a(:);b=b(:);na=numel(a);nb=numel(b);
I=zeros(n*w*w,1);J=I;v=I;va=I;vb=I;
for k=1:n
    [ja,wa,da]=weights(a,qa(k),w);[jb,wb,db]=weights(b,qb(k),w);
    [ib,ia]=ndgrid(jb,ja);
    if strcmp(chart,'core'),j=ib+(ia-1)*nb;else,j=ia+(ib-1)*na;end
    slots=(k-1)*w*w+(1:w*w);
    I(slots)=k;J(slots)=j(:);W=wb*wa';WA=wb*da';WB=db*wa';
    v(slots)=W(:);va(slots)=WA(:);vb(slots)=WB(:);
end
plan=struct('E',sparse(I,J,v,n,na*nb),'Da',sparse(I,J,va,n,na*nb), ...
    'Db',sparse(I,J,vb,n,na*nb),'x',x,'y',y,'chart',chart,'shape',shape);
end

function [j,v,d]=weights(z,q,w)
[~,i]=min(abs(z-q));first=min(max(i-floor(w/2),1),numel(z)-w+1);j=(first:first+w-1)';
dz=z(j)-q;scale=max(abs(dz));M=(dz'/scale).^((0:w-1)');
b=zeros(w,1);b(1)=1;v=M\b;b(:)=0;b(2)=1/scale;d=M\b;
end

function v=evaluate_plan(p,U)
psi=p.E*U(:);a=p.Da*U(:);b=p.Db*U(:);
if strcmp(p.chart,'core'),v=struct('psi',psi,'u1',-b,'u2',a);
else
    r=hypot(p.x,p.y);co=p.x./r;si=p.y./r;
    v=struct('psi',psi,'u1',-(si.*a+co.*b)./r,'u2',(co.*a-si.*b)./r);
end
end
