classdef poisson < handle
    %SOLVER Overlapping Galerkin FEM, explicit exterior data, same-psi u.
    % Core coefficients have y rows/x columns. Annulus coefficients have
    % log-radius rows/theta columns. No pointwise D2 enters the FEM solve.
    properties (SetAccess=private)
        g; c; p; coupling; buildSeconds; auditMaps; elementSummary
        Counters=struct('builds',1,'solves',0,'factorReuses',0)
    end
    properties (Access=private)
        returns
    end
    methods
        function obj=poisson(g)
            timer=tic;obj.g=g;
            c=struct('xAxis',ipm_cr.space.basis.build(g.core.x), ...
                'yAxis',ipm_cr.space.basis.build(g.core.y));
            c.x=g.core.x;c.y=g.core.y;c.X=g.core.X;c.Y=g.core.Y;
            c.nx=numel(c.x);c.ny=numel(c.y);[c.qX,c.qY]=meshgrid(c.xAxis.qx,c.yAxis.qx);
            c.ix=2:c.nx-1;c.iy=2:c.ny-1;
            free=false(c.ny,c.nx);free(c.iy,c.ix)=true;c.free=find(free);
            trace=false(c.ny,c.nx);trace(:,[1,end])=true;trace(end,:)=true;trace(1,:)=false;
            if isfield(g.core,'stair')
                c=assemble_staircase(c,g.core);
            else
                c.traceIndex=find(trace);c.A=kron(c.xAxis.M,c.yAxis.K)+kron(c.xAxis.K,c.yAxis.M);
            end
            c.factor=cholesky_factor(c.A(c.free,c.free));
            p=struct('tAxis',ipm_cr.space.basis.build(log(g.polar.r)), ...
                'thetaAxis',ipm_cr.space.basis.build(g.polar.theta));
            p.t=p.tAxis.nodes;
            p.r=g.polar.r;p.theta=g.polar.theta;p.X=g.polar.X;p.Y=g.polar.Y;
            p.nr=numel(p.t);p.nt=numel(p.theta);p.jj=2:p.nt-1;
            [th,tq]=meshgrid(p.thetaAxis.qx,p.tAxis.qx);p.qR=exp(tq);
            p.qX=p.qR.*cos(th);p.qY=p.qR.*sin(th);
            radial=p.tAxis;
            p.sourceMass=radial.Bq'*((radial.qw.*p.qR(:,1).^2).*radial.Bq);
            p=ipm_cr.space.boundary.prepare(p);
            m=numel(p.jj);
            free=false(p.nr,p.nt);free(p.radialFree,p.jj)=true;p.free=find(free);
            obj.elementSummary=degree_summary(c,p);
            if ~isfield(c,'traceLoad')
                c.traceX=c.X(c.traceIndex);c.traceY=c.Y(c.traceIndex);
            end
            mapCA=ipm_cr.space.velocity.plan(p.tAxis,p.thetaAxis,c.traceX,c.traceY,'polar');
            mapAC=ipm_cr.space.velocity.plan(c.xAxis,c.yAxis,p.X(1,p.jj)',p.Y(1,p.jj)','core');
            p.innerTraceEval=p.thetaAxis.E(p.jj,p.jj);
            p.innerTraceFit=decomposition(p.innerTraceEval,'lu');
            if isfield(c,'outsideSupport')
                c.splineOutside=ipm_cr.space.velocity.plan(p.tAxis,p.thetaAxis, ...
                    c.X(c.outsideSupport),c.Y(c.outsideSupport),'polar');
            end
            obj.c=c;obj.p=p;
            if isfield(c,'traceLoad')
                T=spline_trace_response(c,p,mapCA.E,mapAC.E);
            else
                T=zeros(m);
                for j=1:m
                    lambda=zeros(1,m);lambda(j)=1;pa=obj.polar_harmonic(lambda);
                    cc=obj.core_solve(zeros(size(c.X)),mapCA.E*pa(:));
                    T(:,j)=mapAC.E*cc(:);
                end
            end
            S=eye(m)-T;
            obj.coupling=struct('RCA',mapCA.E,'RAC',mapAC.E,'S',S, ...
                'factor',decomposition(S,'lu'),'condition',cond(S));
            obj.returns=struct('core',ipm_cr.space.velocity.native(c.xAxis,c.yAxis,g.core,'core'), ...
                'polar',ipm_cr.space.velocity.native(p.tAxis,p.thetaAxis,g.polar,'polar'));
            [x,y,coverage]=audit_cloud(g);
            obj.auditMaps=struct('x',x,'y',y, ...
                'core',ipm_cr.space.velocity.plan(c.xAxis,c.yAxis,x,y,'core'), ...
                'polar',ipm_cr.space.velocity.plan(p.tAxis,p.thetaAxis,x,y,'polar'),'coverage',coverage);
            obj.buildSeconds=toc(timer);
        end

        function [flow,info]=solve(obj,source,boundary)
            if nargin<3,boundary=[];end
            timer=tic;phase=tic;c=obj.c;p=obj.p;
            load=ipm_cr.space.boundary.load(p,boundary);
            if isa(source,'function_handle')
                qc=sample_callback(source,c.qX,c.qY);qp=sample_callback(source,p.qX,p.qY);
                if isfield(c,'stairQuadrature'),qc(~c.stairQuadrature)=0;end
                Fc=c.yAxis.Bq'*(c.yAxis.qw.*qc.*c.xAxis.qw')*c.xAxis.Bq;
                Fp=p.tAxis.Bq'*(p.tAxis.qw.*(p.qR.^2.*qp).*p.thetaAxis.qw')*p.thetaAxis.Bq;
                sourceMode='analytic';quadratureCount=numel(qc)+numel(qp);
            else
                validate_source(source,obj.g);
                source=ipm_cr.mesh.exchange(source,obj.g);
                core=spline_coefficients(c.yAxis,c.xAxis,source.core);
                if isfield(c,'sourceMass')
                    Fc=reshape(c.sourceMass*core(:),c.ny,c.nx);
                else
                    Fc=c.yAxis.M*core*c.xAxis.M.';
                end
                polar=spline_coefficients(p.tAxis,p.thetaAxis,source.polar);
                Fp=p.sourceMass*polar*p.thetaAxis.M.';
                sourceMode='cached_native_spline_load';quadratureCount=0;
            end
            Fp(end,:)=Fp(end,:)+load.';
            sourceSeconds=toc(phase);phase=tic;
            modalRhs=Fp(:,p.jj)*p.V;
            [pa0,~]=obj.polar_solve(Fp,zeros(1,numel(p.jj)),modalRhs);
            cc0=obj.core_solve(Fc,obj.coupling.RCA*pa0(:));
            b=obj.coupling.RAC*cc0(:);lambda=obj.coupling.factor\b;
            [pa,~]=obj.polar_solve(Fp,lambda.',modalRhs);
            coreTrace=obj.coupling.RCA*pa(:);
            [cc,Fc]=obj.core_solve(Fc,coreTrace);
            coreBackward=weak_backward(c.A,cc,Fc,c.free);
            polarBackward=weak_backward(p.A,pa,Fp,p.free);backward=max(coreBackward,polarBackward);
            inner=p.innerTraceEval*pa(1,p.jj).';
            trace=max(norm(obj.coupling.RAC*cc(:)-lambda,inf),norm(inner-lambda,inf));
            if ~isfield(c,'traceLoad')
                ccNodes=c.yAxis.E*cc*c.xAxis.E.';
                trace=max(trace,norm(coreTrace-ccNodes(c.traceIndex),inf));
            end
            essential=max([max(abs(cc(1,:))),max(abs(pa(:,[1,end])),[],'all'), ...
                norm(inner-lambda,inf)]);
            solveSeconds=toc(phase);phase=tic;
            flow=struct('core',ipm_cr.space.velocity.evaluate(obj.returns.core,cc), ...
                'polar',ipm_cr.space.velocity.evaluate(obj.returns.polar,pa), ...
                'coefficients',struct('core',cc,'polar',pa), ...
                'representation','quintic tensor B-spline Galerkin');
            if isfield(c,'splineOutside')
                outside=ipm_cr.space.velocity.evaluate(c.splineOutside,pa);
                for names={'psi','u1','u2'}
                    k=names{1};flow.core.(k)(c.outsideSupport)=outside.(k);
                end
            end
            returnSeconds=toc(phase);[oc,op]=obj.overlap(flow);
            obj.Counters.solves=obj.Counters.solves+1;obj.Counters.factorReuses=obj.Counters.factorReuses+1;
            info=struct('solverCompleted',true,'accuracyAssessed',false, ...
                'backwardError',backward,'weakBackwardError',backward, ...
                'coreWeakBackwardError',coreBackward,'polarWeakBackwardError',polarBackward, ...
                'traceAbsolute',trace,'couplingCondition',obj.coupling.condition, ...
                'essentialBoundaryAbsolute',essential,'sourceMode',sourceMode, ...
                'sourceSeconds',sourceSeconds,'solveSeconds',solveSeconds,'returnSeconds',returnSeconds, ...
                'kernelSeconds',sourceSeconds+solveSeconds+returnSeconds,'totalSeconds',toc(timer), ...
                'buildSeconds',obj.buildSeconds,'degree',c.xAxis.degree,'quadOrder',c.xAxis.quadratureOrder, ...
                'massLumped',false,'discretizationContract','overlap_quintic_bspline_galerkin_v1', ...
                'coefficientDof',numel(cc)+numel(pa),'freeWeakDof',numel(c.free)+numel(p.free), ...
                'nodalLayout','native','sameNativeGrid',true, ...
                'independentPoissonGridUsed',false,'noIndependentPoissonGrid',true, ...
                'nativeReturnRegriddingUsed',false, ...
                'actualDegreeRange',obj.elementSummary.range, ...
                'actualElementDegreeCounts',obj.elementSummary.counts, ...
                'outerElementDegreePolicy',obj.p.tAxis.boundaryDegreePolicy, ...
                'sourceQuadratureCount',quadratureCount,'sparseOperatorNnz',nnz(c.A)+nnz(p.A), ...
                'linearSolveComponentwiseBackwardError',backward,'linearSolveBackwardError',backward, ...
                'fixedWeakResidualCorrections',1, ...
                'linearSolveRule','cached factor solve followed by one unconditional original-weak-residual correction per patch solve', ...
                'diagnosticScope','assembled weak free rows, traces and boundary values measured only; no experiment acceptance threshold', ...
                'overlapValueMismatch',maximum_abs(oc.psi(:)-op.psi(:)), ...
                'overlapVelocityMismatch',maximum_abs([oc.u1(:)-op.u1(:);oc.u2(:)-op.u2(:)]), ...
                'overlapAuditCoverage',obj.auditMaps.coverage, ...
                'outerBoundary','forced_dtn', ...
                'outerSourceTailCertified',false,'strongResidualAvailable',false,'continuumAccuracyCertified',false, ...
                'sharedFaceDerivative','continuous spline derivative at simple interior knots', ...
                'linearSolverStorage',factor_storage(obj), ...
                'complexityScope','sparse core, banded radial modes, dense angular transforms and trace Schur solve');
            info.finiteOutput=finite_flow(flow);
            if isfield(c,'traceLoad')
                info.coreWeakTraceMismatch=norm(c.traceEval*cc(:)-coreTrace,inf);
            end
        end

        function v=evaluate(obj,flow,x,y,chart)
            if strcmp(chart,'core'),plan=ipm_cr.space.velocity.plan(obj.c.xAxis,obj.c.yAxis,x,y,'core');
            else,assert(strcmp(chart,'polar'),'ipm_cr:FEMChart');plan=ipm_cr.space.velocity.plan(obj.p.tAxis,obj.p.thetaAxis,x,y,'polar');end
            v=ipm_cr.space.velocity.evaluate(plan,flow.coefficients.(chart));
        end
        function [a,b]=overlap(obj,flow)
            a=ipm_cr.space.velocity.evaluate(obj.auditMaps.core,flow.coefficients.core);
            b=ipm_cr.space.velocity.evaluate(obj.auditMaps.polar,flow.coefficients.polar);
        end
        function [U,F]=core_solve(obj,F,trace)
            c=obj.c;U=zeros(c.ny,c.nx);
            if isfield(c,'traceLoad')
                F=F+reshape(c.traceLoad*trace,size(F));
            else
                U(c.traceIndex)=trace;
                top=U(end,:);
                U(:,1)=c.yAxis.interpolation\U(:,1);
                U(:,end)=c.yAxis.interpolation\U(:,end);
                U(end,:)=(c.xAxis.interpolation\top.').';
            end
            U=weak_solve(c.A,c.factor,U,F,c.free);
        end
        function [U,Z]=polar_solve(obj,F,inner,modalRhs)
            p=obj.p;
            inner=(p.innerTraceFit\inner(:)).';
            if nargin<4,modalRhs=F(:,p.jj)*p.V;end
            inner=inner*p.M*p.V;
            Z=zeros(p.nr,numel(p.mu));
            for j=1:numel(p.mu)
                A=p.operators{j};z=zeros(p.nr,1);z(1)=inner(j);
                b=modalRhs(:,j)-A*z;f=p.radialFree;z(f)=p.factors{j}\b(f);
                Z(:,j)=z;
            end
            U=zeros(p.nr,p.nt);U(:,p.jj)=Z*p.V.';
            [U,Z]=polar_correct(p,U,F,Z);
        end
        function U=polar_harmonic(obj,inner)
            p=obj.p;
            inner=(p.innerTraceFit\inner(:)).';
            Z=p.harmonic.*(inner*p.M*p.V);U=zeros(p.nr,p.nt);U(:,p.jj)=Z*p.V.';
        end
    end
end

function coeff=spline_coefficients(row,col,native)
% Preserve small far-field values before the radial r^2 integration.
coeff=row.interpolation\native;
coeff=(col.interpolation\coeff.').';
end
function result=degree_summary(c,p)
axes={c.xAxis,c.yAxis,p.tAxis,p.thetaAxis};names={'coreX','coreY','polarLogRadius','polarTheta'};
allDegrees=[];counts=struct();
for k=1:4
    values=axes{k}.elementDegrees;degrees=unique(values)';number=zeros(size(degrees));
    for j=1:numel(degrees),number(j)=nnz(values==degrees(j));end
    counts.(names{k})=struct('degree',degrees,'count',number,'elementCount',numel(values), ...
        'partitionConditioning',axes{k}.partitionConditioning,'boundaryDegreePolicy',axes{k}.boundaryDegreePolicy);
    allDegrees=[allDegrees;values]; %#ok<AGROW>
end
result=struct('range',[min(allDegrees),max(allDegrees)],'counts',counts);
end
function [x,y,coverage]=audit_cloud(g)
% Audit the actual common ACTIVE intersection, which may consist of near-
% wall sectors rather than a complete semicircle on a coarse/remeshed grid.
theta=[1e-5,((1:25)+.17)*pi/26.5,pi-1e-5];co=cos(theta);si=sin(theta);
xbound=g.core.activeXBounds(2)*ones(size(theta));
xbound(co<0)=-g.core.activeXBounds(1);
hi=min(xbound./abs(co),g.core.activeYBounds(2)./si);lo=g.polar.innerActiveRadius;
if isfield(g.core,'stair'),hi=ipm_cr.mesh.staircase('ray',g.core.stair,theta).';end
keep=hi>lo+128*eps(max([1,lo,max(hi)]));
kept=theta(keep);r=lo+[.2;.5;.8]*(hi(keep)-lo);th=repmat(kept,3,1);
x=r.*cos(th);y=r.*sin(th);
inside=all(x>=g.core.activeXBounds(1)&x<=g.core.activeXBounds(2)& ...
    y>=g.core.activeYBounds(1)&y<=g.core.activeYBounds(2)&r>lo,'all');
coverage=struct('kind','offgrid_common_active_radial_sectors','requestedAngles',numel(theta), ...
    'keptAngles',nnz(keep),'theta',kept,'queryCount',numel(x),'radialFractions',[.2,.5,.8], ...
    'innerActiveRadius',lo,'outerActiveRadiusByAngle',hi(keep),'allRequestedAnglesCovered',all(keep), ...
    'wholeSemicircleRequired',false,'available',~isempty(x),'actualActiveIntersectionVerified',inside);
end
function f=cholesky_factor(A)
% The mathematical free weak system is positive definite. Diagonal
% congruence balances its representation; diagnostics use the original A.
diagonal=diag(A);assert(all(isfinite(diagonal))&&all(diagonal>0), ...
    'ipm_cr:FEMCholesky','FEM free weak rows require positive finite diagonal entries.');
scale=1./sqrt(diagonal);D=spdiags(scale,0,numel(scale),numel(scale));
B=D*((A+A')/2)*D;
% Retain the sparse Cholesky solver, including its internal ordering and
% solve workspace. Repeated backslashes on an exported L rebuild that work.
f=struct('solver',decomposition(B,'chol','lower'),'scale',scale);
end
function T=spline_trace_response(c,p,toCore,toPolar)
% Batch harmonic responses. Only retain the trace Schur matrix.
mode=p.V.'*p.M*(p.innerTraceFit\eye(numel(p.jj)));
polarTrace=zeros(size(toCore,1),size(mode,1));
for k=1:numel(p.mu)
    response=zeros(p.nr,p.nt);response(:,p.jj)=p.harmonic(:,k)*p.V(:,k).';
    polarTrace(:,k)=toCore*response(:);
end
load=c.traceLoad(c.free,:)*polarTrace;
T=(toPolar(:,c.free)*cholesky_solve(c.factor,load))*mode;
end
function x=cholesky_solve(f,b)
x=f.scale.*(f.solver\(f.scale.*b));
end
function U=weak_solve(A,f,U,F,free)
r=F(:)-A*U(:);U(free)=cholesky_solve(f,r(free));
% A fixed second factor application is part of the linear solve, including
% setup harmonics. It never tests a tolerance or selects between answers.
r=F(:)-A*U(:);U(free)=U(free)+cholesky_solve(f,r(free));
end
function [U,Z]=polar_correct(p,U,F,Z)
% Correct the original physical weak rows, including angular roundoff,
% rather than refining only the separated radial reference systems.
r=F-reshape(p.A*U(:),size(U));rhs=r(:,p.jj)*p.V;delta=zeros(size(Z));
for j=1:numel(p.mu)
    delta(p.radialFree,j)=p.factors{j}\rhs(p.radialFree,j);
end
Z=Z+delta;U(:,p.jj)=Z*p.V.';
end
function s=factor_storage(obj)
radialMass=obj.p.sourceMass;
a=obj.p.operators;b=obj.p.factors;v=obj.p.V;h=obj.p.harmonic;
modal=whos('a','b','v','h');load=whos('radialMass');
s=struct('patchFactorNnz',NaN,'patchFactorBytes',NaN, ...
    'patchFactorStorageMeasured',false,'patchFactorRepresentation','cached_sparse_cholesky_decomposition', ...
    'patchFactorStorageScope','MATLAB does not expose internal factor storage; measure process memory separately', ...
    'traceDimension',size(obj.coupling.S,1),'dtNMatrixIncludedInSparseOperator',true, ...
    'polarModalPlanBytes',sum([modal.bytes]),'sourceLoadCacheBytes',load.bytes);
if isfield(obj.c,'sourceMass')
    stairMass=obj.c.sourceMass;
    entry=whos('stairMass');s.sourceLoadCacheBytes=s.sourceLoadCacheBytes+entry.bytes;
end
end
function value=weak_backward(A,U,F,free)
r=A*U(:)-F(:);den=abs(A)*abs(U(:))+abs(F(:));
value=max(abs(r(free))./max(den(free),realmin));
end
function validate_source(q,g)
assert(isstruct(q)&&isscalar(q)&&all(isfield(q,{'core','polar'}))&& ...
    isnumeric(q.core)&&isreal(q.core)&&isequal(size(q.core),size(g.core.X))&& ...
    isnumeric(q.polar)&&isreal(q.polar)&&isequal(size(q.polar),size(g.polar.X)), ...
    'ipm_cr:FEMSource','Native input is the already differentiated source on both complete transport arrays.');
end
function value=maximum_abs(values)
if isempty(values),value=NaN;else,value=max(abs(values(:)));end
end
function yes=finite_flow(flow)
yes=all(isfinite(flow.core.psi),'all')&&all(isfinite(flow.polar.psi),'all')&& ...
    all(isfinite(flow.core.u1),'all')&&all(isfinite(flow.core.u2),'all')&& ...
    all(isfinite(flow.polar.u1),'all')&&all(isfinite(flow.polar.u2),'all');
end

function q=sample_callback(callback,x,y)
% Evaluate a real-valued callback on a FEM query array.
q=callback(x,y);if isscalar(q),q=repmat(q,size(x));end
assert(isnumeric(q)&&isreal(q)&&isequal(size(q),size(x)), ...
    'ipm_cr:FEMSource','Callbacks must return real numeric values matching the supplied query shape.');
end

function c=assemble_staircase(c,core)
% Integrate whole spline cells and apply symmetric Nitsche traces.
x=c.xAxis;y=c.yAxis;nx=c.nx;ny=c.ny;ndof=nx*ny;
% Retain the legacy Poisson collar and cover it by whole knot cells.
segments=core.stair.elements;
lx=ipm_cr.mesh.staircase('partition',c.x);
ly=ipm_cr.mesh.staircase('partition',c.y);
mx=x.breaks(1:end-1)<lx.breaks(2:end).' & x.breaks(2:end)>lx.breaks(1:end-1).';
my=y.breaks(1:end-1)<ly.breaks(2:end).' & y.breaks(2:end)>ly.breaks(1:end-1).';
mask=spones(sparse(my)*sparse(segments)*sparse(mx'))>0;
rowsPerCell=x.quadratureOrder;
qmask=repelem(mask,rowsPerCell,rowsPerCell);
c.stairQuadrature=full(qmask);
ii=cell(size(mask,1),1);jj=ii;vv=ii;
mi=ii;mj=ii;mv=ii;
for j=1:size(mask,1)
    q=(j-1)*rowsPerCell+(1:rowsPerCell);
    b=y.Bq(q,:);d=y.Dq(q,:);w=y.qw(q);
    use=c.stairQuadrature(q(1),:);
    ex=x.Bq(use,:);dx=x.Dq(use,:);wx=x.qw(use);
    xMass=ex'*(wx.*ex);yMass=b'*(w.*b);
    [ii{j},jj{j},vv{j}]=find(kron(xMass,d'*(w.*d))+kron(dx'*(wx.*dx),yMass));
    [mi{j},mj{j},mv{j}]=find(kron(xMass,yMass));
end
c.A=sparse(vertcat(ii{:}),vertcat(jj{:}),vertcat(vv{:}),ndof,ndof);
c.sourceMass=sparse(vertcat(mi{:}),vertcat(mj{:}),vertcat(mv{:}),ndof,ndof);
% Edge quadrature is in the same knot partition as the volume space.
% With whole cells the inverse trace estimate is p^2/h in the normal direction.
pad=false(size(mask)+2);pad(2:end-1,2:end-1)=mask;
E={};normal={};weights={};penalty={};xs={};ys={};edgeCount=0;
for dir=[-1 1]
    for ax=1:2
        if ax==1
            edges=mask & ~pad(2:end-1,2+dir:end-1+dir);
        else
            edges=mask & ~pad(2+dir:end-1+dir,2:end-1);
            if dir==-1,edges(1,:)=false;end % strong homogeneous wall
        end
        [j,i]=find(edges);
        for n=1:numel(i)
            col=i(n);row=j(n);edge=(dir+1)/2;
            if ax==1
                q=(row-1)*rowsPerCell+(1:rowsPerCell);
                b=y.Bq(q,:);w=y.qw(q);v=x.breaks(col+edge);
                [e,d]=ipm_cr.space.basis.evaluate(x,v);
                B=kron(e,b);N=dir*kron(d,b);
                xx=v*ones(size(w));yy=y.qx(q);h=x.breaks(col+1)-x.breaks(col);
            else
                q=(col-1)*rowsPerCell+(1:rowsPerCell);
                b=x.Bq(q,:);w=x.qw(q);v=y.breaks(row+edge);
                [e,d]=ipm_cr.space.basis.evaluate(y,v);
                B=kron(b,e);N=dir*kron(b,d);
                xx=x.qx(q);yy=v*ones(size(w));h=y.breaks(row+1)-y.breaks(row);
            end
            edgeCount=edgeCount+1;
            E{edgeCount}=B;normal{edgeCount}=N;weights{edgeCount}=w; %#ok<AGROW>
            penalty{edgeCount}=2*(x.degree+1)^2/h*w; %#ok<AGROW>
            xs{edgeCount}=xx;ys{edgeCount}=yy; %#ok<AGROW>
        end
    end
end
c.traceEval=vertcat(E{:});N=vertcat(normal{:});weights=vertcat(weights{:});
wc=vertcat(penalty{:});
c.traceLoad=(wc.*c.traceEval-weights.*N)';
c.A=c.A+c.traceLoad*c.traceEval-c.traceEval'*(weights.*N);
c.traceX=vertcat(xs{:});c.traceY=vertcat(ys{:});
c.A=(c.A+c.A')/2;
free=reshape(full(diag(c.A))~=0,ny,nx);free(1,:)=false;c.free=find(free);
c.traceIndex=[];
ix=discretize(c.x,x.breaks);iy=discretize(c.y,y.breaks);
% At a knot either adjacent cell gives the same interior spline trace.
onX=c.x(:)==x.breaks(ix);onY=c.y(:)==y.breaks(iy);
ixLeft=max(1,ix(:)-onX);iyLeft=max(1,iy(:)-onY);
nodal=mask(iy,ix)|mask(iyLeft,ix)|mask(iy,ixLeft)|mask(iyLeft,ixLeft);
assert(all(nodal(core.active)),'ipm_cr:SplineCoverage', ...
    'The spline core must cover every active transport site.');
c.outsideSupport=find(~nodal);
end
