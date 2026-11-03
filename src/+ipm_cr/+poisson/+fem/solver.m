classdef solver < handle
    %SOLVER Overlapping Galerkin FEM, explicit exterior data, same-psi u.
    % Core coefficients have y rows/x columns. Annulus coefficients have
    % log-radius rows/theta columns. No pointwise D2 enters the FEM solve.
    properties (SetAccess=private)
        g; c; p; coupling; config; buildSeconds; auditMaps; elementSummary
        Counters=struct('builds',1,'solves',0,'factorReuses',0,'sourceMapBuilds',0)
    end
    properties (Access=private)
        sourceMaps=[]
        returns
    end
    methods
        function obj=solver(g,options)
            if nargin<2,options=struct();end
            timer=tic;cfg=struct('degree',5,'quadOrder',8,'stride',[], ...
                'nodalLayout','native','outerBoundary','green_dirichlet','linearSolver','sparse_direct');
            assert(isstruct(options)&&isscalar(options)&&all(ismember(fieldnames(options),fieldnames(cfg))), ...
                'ipm_cr:FEMOptions','Unknown FEM option.');
            for name=fieldnames(options)',cfg.(name{1})=options.(name{1});end
            validateattributes(cfg.degree,{'numeric'},{'scalar','integer','>=',2,'<=',6});
            validateattributes(cfg.quadOrder,{'numeric'},{'scalar','integer','>=',cfg.degree+1});
            cfg.nodalLayout=char(validatestring(cfg.nodalLayout,{'native','gll_reference'}));
            native=strcmp(cfg.nodalLayout,'native');
            if native
                assert(isempty(cfg.stride)||(isnumeric(cfg.stride)&&isscalar(cfg.stride)&& ...
                    isreal(cfg.stride)&&isfinite(cfg.stride)&&cfg.stride==cfg.degree), ...
                    'ipm_cr:FEMNativeGrid','Native FEM stride must be empty or equal its element degree.');
            else
                if isempty(cfg.stride),cfg.stride=4;end
                validateattributes(cfg.stride,{'numeric'},{'scalar','integer','positive'});
            end
            assert(any(strcmp(cfg.outerBoundary,{'green_dirichlet','forced_dtn'})), ...
                'ipm_cr:FEMOptions','Use green_dirichlet or forced_dtn.');
            assert(any(strcmp(cfg.linearSolver,{'sparse_direct','modal'})), ...
                'ipm_cr:FEMOptions','Use sparse_direct or modal.');
            cfg.outerBoundary=char(cfg.outerBoundary);cfg.linearSolver=char(cfg.linearSolver);
            sparseDirect=strcmp(cfg.linearSolver,'sparse_direct');obj.g=g;obj.config=cfg;
            c=struct('xAxis',make_axis(g.core.x,cfg,false),'yAxis',make_axis(g.core.y,cfg,false));
            if native
                c.x=g.core.x;c.y=g.core.y;c.X=g.core.X;c.Y=g.core.Y;
            else
                c.x=c.xAxis.nodes;c.y=c.yAxis.nodes;[c.X,c.Y]=meshgrid(c.x,c.y);
            end
            c.nx=numel(c.x);c.ny=numel(c.y);[c.qX,c.qY]=meshgrid(c.xAxis.qx,c.yAxis.qx);
            c.ix=2:c.nx-1;c.iy=2:c.ny-1;
            free=false(c.ny,c.nx);free(c.iy,c.ix)=true;c.free=find(free);
            trace=false(c.ny,c.nx);trace(:,[1,end])=true;trace(end,:)=true;trace(1,:)=false;
            c.traceIndex=find(trace);c.A=kron(c.xAxis.M,c.yAxis.K)+kron(c.xAxis.K,c.yAxis.M);
            if sparseDirect
                c.factor=cholesky_factor(c.A(c.free,c.free));
            else
                [c.Vx,c.lx,c.eigenDefectX]=mass_eig(c.xAxis.K(c.ix,c.ix),c.xAxis.M(c.ix,c.ix));
                [c.Vy,c.ly,c.eigenDefectY]=mass_eig(c.yAxis.K(c.iy,c.iy),c.yAxis.M(c.iy,c.iy));
                c.denominator=c.ly+c.lx.';
            end
            p=struct('tAxis',make_axis(log(g.polar.r),cfg,true), ...
                'thetaAxis',make_axis(g.polar.theta,cfg,false));
            p.t=p.tAxis.nodes;
            if native
                % exp(log(r)) is not bitwise r. Physical unknowns and
                % Cartesian point arrays are exactly the transport arrays.
                p.r=g.polar.r;p.theta=g.polar.theta;p.X=g.polar.X;p.Y=g.polar.Y;
            else
                p.r=exp(p.t);p.theta=p.thetaAxis.nodes.';
                [th,rr]=meshgrid(p.theta,p.r);p.X=rr.*cos(th);p.Y=rr.*sin(th);p.Y(:,[1,end])=0;
            end
            p.nr=numel(p.t);p.nt=numel(p.theta);p.jj=2:p.nt-1;
            [th,tq]=meshgrid(p.thetaAxis.qx,p.tAxis.qx);p.qR=exp(tq);
            p.qX=p.qR.*cos(th);p.qY=p.qR.*sin(th);
            p.M=p.thetaAxis.M(p.jj,p.jj);
            if ~sparseDirect||strcmp(cfg.outerBoundary,'forced_dtn')
                % The angular eigensystem is needed only during sparse
                % DtN setup; sparse warm solves never transform by it.
                [p.V,p.mu,p.eigenDefect]=mass_eig(p.thetaAxis.K(p.jj,p.jj),p.M);
            end
            p.radialFree=2:p.nr-1;
            if strcmp(cfg.outerBoundary,'forced_dtn'),p.radialFree=2:p.nr;end
            m=numel(p.jj);
            if ~sparseDirect
                p.operators=cell(1,m);p.factors=cell(1,m);p.harmonic=zeros(p.nr,m);
                for j=1:m
                    A=p.tAxis.K+p.mu(j)*p.tAxis.M;
                    if strcmp(cfg.outerBoundary,'forced_dtn'),A(end,end)=A(end,end)+sqrt(p.mu(j));end
                    p.operators{j}=A;p.factors{j}=decomposition(A(p.radialFree,p.radialFree),'chol');
                    h=zeros(p.nr,1);h(1)=1;
                    h(p.radialFree)=p.factors{j}\(-A(p.radialFree,1));p.harmonic(:,j)=h;
                end
            end
            p.A=kron(p.thetaAxis.M,p.tAxis.K)+kron(p.thetaAxis.K,p.tAxis.M);
            if strcmp(cfg.outerBoundary,'forced_dtn')
                boundaryMatrix=zeros(p.nt);boundaryMatrix(p.jj,p.jj)=p.M*p.V*diag(sqrt(p.mu))*p.V.'*p.M;
                endRow=sparse(p.nr,p.nr);endRow(end,end)=1;
                p.A=p.A+kron(sparse(boundaryMatrix),endRow);
            end
            free=false(p.nr,p.nt);free(p.radialFree,p.jj)=true;p.free=find(free);
            if sparseDirect,p.factor=cholesky_factor(p.A(p.free,p.free));end
            obj.c=c;obj.p=p;obj.elementSummary=degree_summary(c,p);
            mapCA=query_plan(p.tAxis,p.thetaAxis,c.X(c.traceIndex),c.Y(c.traceIndex),'polar','average');
            mapAC=query_plan(c.xAxis,c.yAxis,p.X(1,p.jj)',p.Y(1,p.jj)','core','average');
            T=zeros(m);response=[];
            if ~sparseDirect,response=zeros(numel(c.X),m);end
            for j=1:m
                lambda=zeros(1,m);lambda(j)=1;pa=obj.polar_harmonic(lambda);
                cc=obj.core_solve(zeros(size(c.X)),mapCA.E*pa(:));
                if ~sparseDirect,response(:,j)=cc(:);end
                T(:,j)=mapAC.E*cc(:);
            end
            S=eye(m)-T;
            obj.coupling=struct('RCA',mapCA.E,'RAC',mapAC.E,'S',S, ...
                'factor',decomposition(S,'lu'),'condition',cond(S));
            if ~sparseDirect,obj.coupling.coreResponse=response;end
            if native
                obj.returns=struct('core',native_plan(c.xAxis,c.yAxis,g.core,'core'), ...
                    'polar',native_plan(p.tAxis,p.thetaAxis,g.polar,'polar'));
            else
                obj.returns=struct('core',query_plan(c.xAxis,c.yAxis,g.core.X,g.core.Y,'core','average'), ...
                    'polar',query_plan(p.tAxis,p.thetaAxis,g.polar.X,g.polar.Y,'polar','average'));
            end
            [x,y,coverage]=audit_cloud(g);
            obj.auditMaps=struct('x',x,'y',y, ...
                'core',query_plan(c.xAxis,c.yAxis,x,y,'core','average'), ...
                'polar',query_plan(p.tAxis,p.thetaAxis,x,y,'polar','average'),'coverage',coverage);
            obj.buildSeconds=toc(timer);
        end

        function [flow,info]=solve(obj,source,boundary,analyticSource)
            if nargin<4,analyticSource=[];end
            if nargin<3,error('ipm_cr:FEMExteriorData','Explicit exterior data are mandatory.');end
            timer=tic;phase=tic;c=obj.c;p=obj.p;
            bd=obj.boundary_data(boundary);
            if isa(source,'function_handle')
                qc=sample(source,c.qX,c.qY);qp=sample(source,p.qX,p.qY);sourceMode='analytic';
            else
                validate_source(source,obj.g);
                % Only active q values are authority. Reconstruct collars
                % simultaneously before any Gauss-point source stencil.
                source=ipm_cr.mesh.fill(source,obj.g);
                if isempty(obj.sourceMaps)
                    [qc,mc]=ipm_cr.mesh.sample(source,obj.g,c.qX,c.qY);
                    [qp,mp]=ipm_cr.mesh.sample(source,obj.g,p.qX,p.qY);
                    assert(mc.allValid&&mp.allValid,'ipm_cr:FEMSourceCoverage');
                    obj.sourceMaps=struct('core',mc,'polar',mp);obj.Counters.sourceMapBuilds=obj.Counters.sourceMapBuilds+1;
                else
                    qc=apply_source(obj.sourceMaps.core,source);qp=apply_source(obj.sourceMaps.polar,source);
                end
                sourceMode='native_unique_owner';
            end
            if ~isempty(analyticSource)
                assert(isa(analyticSource,'function_handle'),'ipm_cr:FEMSource');
                qc=qc+sample(analyticSource,c.qX,c.qY);qp=qp+sample(analyticSource,p.qX,p.qY);
            end
            Fc=c.yAxis.Bq'*(c.yAxis.qw.*qc.*c.xAxis.qw')*c.xAxis.Bq;
            Fp=p.tAxis.Bq'*(p.tAxis.qw.*(p.qR.^2.*qp).*p.thetaAxis.qw')*p.thetaAxis.Bq;
            if strcmp(obj.config.outerBoundary,'forced_dtn'),Fp(end,:)=Fp(end,:)+bd.load.';end
            sourceSeconds=toc(phase);phase=tic;
            [pa0,~]=obj.polar_solve(Fp,zeros(1,numel(p.jj)),bd);
            cc0=obj.core_solve(Fc,obj.coupling.RCA*pa0(:));
            b=obj.coupling.RAC*cc0(:);lambda=obj.coupling.factor\b;
            if strcmp(obj.config.linearSolver,'sparse_direct')
                % Re-solve the patches at their final traces. This avoids
                % the Ncore-by-trace dense response cache and its warm
                % multiplication, and preserves the original weak rows.
                [pa,~]=obj.polar_solve(Fp,lambda.',bd);
                cc=obj.core_solve(Fc,obj.coupling.RCA*pa(:));
            else
                cc=cc0+reshape(obj.coupling.coreResponse*lambda,size(cc0));
                pa=pa0+obj.polar_harmonic(lambda.');
            end
            coreBackward=weak_backward(c.A,cc,Fc,c.free);
            polarBackward=weak_backward(p.A,pa,Fp,p.free);backward=max(coreBackward,polarBackward);
            trace=max(norm(obj.coupling.RAC*cc(:)-lambda,inf), ...
                norm(obj.coupling.RCA*pa(:)-cc(c.traceIndex),inf));
            essential=max([max(abs(cc(1,:))),max(abs(pa(:,[1,end])),[],'all'), ...
                norm(pa(1,p.jj)-lambda.',inf)]);
            outerEssential=NaN;
            if strcmp(obj.config.outerBoundary,'green_dirichlet')
                outerEssential=norm(pa(end,:)-bd.values,inf);essential=max(essential,outerEssential);
            end
            solveSeconds=toc(phase);phase=tic;
            flow=struct('core',evaluate_plan(obj.returns.core,cc), ...
                'polar',evaluate_plan(obj.returns.polar,pa), ...
                'coefficients',struct('core',cc,'polar',pa), ...
                'representation',sprintf('C0 %s tensor Galerkin streamfunction; actual element degrees %d-%d; analytical element derivatives', ...
                    obj.config.nodalLayout,obj.elementSummary.range));
            returnSeconds=toc(phase);[oc,op]=obj.overlap(flow);
            obj.Counters.solves=obj.Counters.solves+1;obj.Counters.factorReuses=obj.Counters.factorReuses+1;
            info=struct('solverCompleted',true,'accuracyAssessed',false, ...
                'backwardError',backward,'weakBackwardError',backward, ...
                'coreWeakBackwardError',coreBackward,'polarWeakBackwardError',polarBackward, ...
                'traceAbsolute',trace,'couplingCondition',obj.coupling.condition, ...
                'essentialBoundaryAbsolute',essential,'outerEssentialBoundaryAbsolute',outerEssential, ...
                'sourceMode',sourceMode,'analyticSourceOffset',~isempty(analyticSource), ...
                'sourceSeconds',sourceSeconds,'solveSeconds',solveSeconds,'returnSeconds',returnSeconds, ...
                'kernelSeconds',sourceSeconds+solveSeconds+returnSeconds,'totalSeconds',toc(timer), ...
                'buildSeconds',obj.buildSeconds,'degree',obj.config.degree,'quadOrder',obj.config.quadOrder, ...
                'massLumped',false,'discretizationContract','overlap_C0_Qp_galerkin_explicit_exterior_v1', ...
                'coefficientDof',numel(cc)+numel(pa),'freeWeakDof',numel(c.free)+numel(p.free), ...
                'nodalLayout',obj.config.nodalLayout,'sameNativeGrid',strcmp(obj.config.nodalLayout,'native'), ...
                'independentPoissonGridUsed',~strcmp(obj.config.nodalLayout,'native'), ...
                'noIndependentPoissonGrid',strcmp(obj.config.nodalLayout,'native'), ...
                'nativeReturnRegriddingUsed',~strcmp(obj.config.nodalLayout,'native'), ...
                'actualDegreeRange',obj.elementSummary.range, ...
                'actualElementDegreeCounts',obj.elementSummary.counts, ...
                'outerElementDegreePolicy',obj.p.tAxis.boundaryDegreePolicy, ...
                'sourceQuadratureCount',numel(qc)+numel(qp),'sparseOperatorNnz',nnz(c.A)+nnz(p.A), ...
                'linearSolveComponentwiseBackwardError',backward,'linearSolveBackwardError',backward, ...
                'fixedWeakResidualCorrections',1, ...
                'linearSolveRule','cached factor solve followed by one unconditional original-weak-residual correction per patch solve', ...
                'diagnosticScope','assembled weak free rows, traces and boundary values measured only; no experiment acceptance threshold', ...
                'overlapValueMismatch',maximum_abs(oc.psi(:)-op.psi(:)), ...
                'overlapVelocityMismatch',maximum_abs([oc.u1(:)-op.u1(:);oc.u2(:)-op.u2(:)]), ...
                'overlapAuditCoverage',obj.auditMaps.coverage, ...
                'exteriorStatus',boundary.exteriorStatus,'outerBoundary',obj.config.outerBoundary, ...
                'outerSourceTailCertified',false,'strongResidualAvailable',false,'continuumAccuracyCertified',false, ...
                'sharedFaceDerivative','8-ulp coordinate-roundtrip snap, then average analytical one-sided derivatives of the same FE field', ...
                'sourceInterpolationAccuracyCertified',false,'globalConformingMeshClaimed',false, ...
                'linearSolver',obj.config.linearSolver,'linearSolverStorage',factor_storage(obj), ...
                'denseModalWarmTransforms',strcmp(obj.config.linearSolver,'modal'), ...
                'complexityScope','warm patch solves scale with stored sparse factors; setup includes trace responses and dense DtN eigensystem; no universal complexity certificate');
            info.finiteOutput=finite_flow(flow);
            info.boundaryCornerValuesSupplied=bd.wallCornerValues;
        end

        function v=evaluate(obj,flow,x,y,chart,side)
            if nargin<6,side='average';end
            if strcmp(chart,'core'),plan=query_plan(obj.c.xAxis,obj.c.yAxis,x,y,'core',side);
            else,assert(strcmp(chart,'polar'),'ipm_cr:FEMChart');plan=query_plan(obj.p.tAxis,obj.p.thetaAxis,x,y,'polar',side);end
            v=evaluate_plan(plan,flow.coefficients.(chart));
        end
        function [a,b]=overlap(obj,flow)
            a=evaluate_plan(obj.auditMaps.core,flow.coefficients.core);
            b=evaluate_plan(obj.auditMaps.polar,flow.coefficients.polar);
        end
        function bd=boundary_data(obj,boundary)
            assert(isstruct(boundary)&&isscalar(boundary)&&all(isfield(boundary,{'kind','exteriorStatus','value'})), ...
                'ipm_cr:FEMExteriorData','Supply kind, exteriorStatus and complete boundary value.');
            assert(strcmp(boundary.kind,obj.config.outerBoundary)&&any(strcmp(boundary.exteriorStatus, ...
                {'known_complete_source','declared_compact_source','modeled_tail_with_bound','manufactured_boundary'})), ...
                'ipm_cr:FEMExteriorData','Unknown exterior information or a different boundary operator is not allowed.');
            p=obj.p;R=p.r(end);bd=struct('values',zeros(1,p.nt),'load',zeros(p.nt,1),'wallCornerValues',[NaN,NaN]);
            if strcmp(boundary.kind,'green_dirichlet')
                if isa(boundary.value,'function_handle'),values=sample(boundary.value,p.X(end,:),p.Y(end,:));
                else,values=reshape(boundary.value,1,[]);end
                assert(numel(values)==p.nt&&isreal(values), ...
                    'ipm_cr:FEMExteriorData','Dirichlet data must match the FEM angular nodes.');
                % Homogeneous wall data define the corner values. Supplied
                % exterior values are retained as diagnostics, not tested
                % against a run-time accuracy tolerance.
                bd.wallCornerValues=values([1,end]);
                values([1,end])=0;bd.values=values;
            else
                theta=p.thetaAxis.qx;
                if isa(boundary.value,'function_handle')
                    q=sample(boundary.value,R*cos(theta),R*sin(theta));
                else
                    values=boundary.value(:);assert(numel(values)==p.nt&&isreal(values), ...
                        'ipm_cr:FEMExteriorData','Nodal Robin data must match FEM theta.');
                    q=p.thetaAxis.Bq*values;
                end
                bd.load=p.thetaAxis.Bq'*(p.thetaAxis.qw.*q(:));
            end
        end

        function U=core_solve(obj,F,trace)
            c=obj.c;U=zeros(c.ny,c.nx);U(c.traceIndex)=trace;
            if strcmp(obj.config.linearSolver,'sparse_direct')
                U=weak_solve(c.A,c.factor,U,F,c.free);return
            end
            rhs=F-reshape(c.A*U(:),size(U));U(c.iy,c.ix)=core_inverse(c,rhs(c.iy,c.ix));
            rhs=F-reshape(c.A*U(:),size(U));
            U(c.iy,c.ix)=U(c.iy,c.ix)+core_inverse(c,rhs(c.iy,c.ix));
        end
        function [U,Z]=polar_solve(obj,F,inner,bd)
            p=obj.p;
            if strcmp(obj.config.linearSolver,'sparse_direct')
                U=zeros(p.nr,p.nt);U(1,p.jj)=inner;
                if strcmp(obj.config.outerBoundary,'green_dirichlet'),U(end,p.jj)=bd.values(p.jj);end
                U=weak_solve(p.A,p.factor,U,F,p.free);
                % Z denotes modal coefficients only in the reference path.
                Z=[];return
            end
            rhs=F(:,p.jj)*p.V;inner=inner*p.M*p.V;
            outer=bd.values(p.jj)*p.M*p.V;Z=zeros(p.nr,numel(p.mu));
            for j=1:numel(p.mu)
                A=p.operators{j};z=zeros(p.nr,1);z(1)=inner(j);
                if strcmp(obj.config.outerBoundary,'green_dirichlet'),z(end)=outer(j);end
                b=rhs(:,j)-A*z;f=p.radialFree;z(f)=p.factors{j}\b(f);
                Z(:,j)=z;
            end
            U=zeros(p.nr,p.nt);U(:,p.jj)=Z*p.V.';
            [U,Z]=polar_correct(p,U,F,Z);
        end
        function U=polar_harmonic(obj,inner)
            p=obj.p;
            if strcmp(obj.config.linearSolver,'sparse_direct')
                U=obj.polar_solve(zeros(p.nr,p.nt),inner,struct('values',zeros(1,p.nt)));return
            end
            Z=p.harmonic.*(inner*p.M*p.V);U=zeros(p.nr,p.nt);U(:,p.jj)=Z*p.V.';
            U=polar_correct(p,U,zeros(size(U)),Z);
        end
    end
end

function a=make_axis(nodes,cfg,outer)
input=nodes;
if strcmp(cfg.nodalLayout,'gll_reference')
    if outer,input=outer_breaks(nodes,cfg.stride);else,input=breaks(nodes,cfg.stride);end
end
a=ipm_cr.poisson.fem.axis(input,cfg.degree,cfg.quadOrder,cfg.nodalLayout,outer);
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
function plan=native_plan(a,b,g,chart)
% Exact native psi is the coefficient array. Only its FE derivatives need
% evaluation, by sparse one-dimensional maps on those same native nodes.
[~,Da]=ipm_cr.poisson.fem.assembly('evaluate',a,a.nodes);
[~,Db]=ipm_cr.poisson.fem.assembly('evaluate',b,b.nodes);
plan=struct('native',true,'chart',chart,'Da',Da,'Db',Db,'shape',size(g.X));
if strcmp(chart,'polar')
    plan.radius=g.r(:);plan.cosTheta=g.cosTheta;plan.sinTheta=g.sinTheta;
end
end
function [x,y,coverage]=audit_cloud(g)
% Audit the actual common ACTIVE intersection, which may consist of near-
% wall sectors rather than a complete semicircle on a coarse/remeshed grid.
theta=[1e-5,((1:25)+.17)*pi/26.5,pi-1e-5];co=cos(theta);si=sin(theta);
xbound=g.core.activeXBounds(2)*ones(size(theta));
xbound(co<0)=-g.core.activeXBounds(1);
hi=min(xbound./abs(co),g.core.activeYBounds(2)./si);lo=g.polar.innerActiveRadius;
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
function x=breaks(nodes,stride)
x=nodes(:);index=unique([1:stride:numel(x),numel(x)]);x=x(index);
end
function x=outer_breaks(nodes,stride)
% Fix the outer-element phase when cropping changes the inner annular row.
x=nodes(:);index=sort(unique([1,numel(x):-stride:1]));x=x(index);
end
function [V,lambda,defect]=mass_eig(K,M)
K=full((K+K')/2);M=full((M+M')/2);L=chol(M,'lower');
H=L\K/L';H=(H+H')/2;[Q,D]=eig(H);[lambda,order]=sort(diag(D));V=L'\Q(:,order);
for j=1:size(V,2),[~,i]=max(abs(V(:,j)));if V(i,j)<0,V(:,j)=-V(:,j);end,end
defect=max(norm(V'*M*V-eye(numel(lambda)),inf), ...
    norm(K*V-M*V*diag(lambda),inf)/max(1,norm(K*V,inf)));
end
function U=core_inverse(c,F)
U=c.Vy*((c.Vy'*F*c.Vx)./c.denominator)*c.Vx';
end
function f=cholesky_factor(A)
% The mathematical free weak system is positive definite. Diagonal
% congruence balances its representation; diagnostics use the original A.
diagonal=diag(A);assert(all(isfinite(diagonal))&&all(diagonal>0), ...
    'ipm_cr:FEMCholesky','FEM free weak rows require positive finite diagonal entries.');
scale=1./sqrt(diagonal);D=spdiags(scale,0,numel(scale),numel(scale));
B=D*((A+A')/2)*D;order=symamd(B);L=chol(B(order,order),'lower');
f=struct('L',L,'order',order(:),'scale',scale);
end
function x=cholesky_solve(f,b)
rhs=f.scale.*b;x=zeros(size(b));x(f.order)=f.L'\(f.L\rhs(f.order));x=f.scale.*x;
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
s=struct('patchFactorNnz',0,'patchFactorBytes',0,'coreResponseBytes',0, ...
    'traceDimension',size(obj.coupling.S,1),'dtNMatrixIncludedInSparseOperator', ...
    strcmp(obj.config.outerBoundary,'forced_dtn'));
if strcmp(obj.config.linearSolver,'sparse_direct')
    for name={'c','p'}
        f=obj.(name{1}).factor;entry=whos('f');
        s.patchFactorNnz=s.patchFactorNnz+nnz(f.L);
        s.patchFactorBytes=s.patchFactorBytes+entry.bytes;
    end
else
    s.coreResponseBytes=8*numel(obj.coupling.coreResponse);
end
end
function value=weak_backward(A,U,F,free)
r=A*U(:)-F(:);den=abs(A)*abs(U(:))+abs(F(:));
value=max(abs(r(free))./max(den(free),realmin));
end
function q=sample(callback,x,y)
q=callback(x,y);if isscalar(q),q=repmat(q,size(x));end
assert(isnumeric(q)&&isreal(q)&&isequal(size(q),size(x)), ...
    'ipm_cr:FEMSource','Callbacks must return real numeric values matching the supplied query shape.');
end
function validate_source(q,g)
assert(isstruct(q)&&isscalar(q)&&all(isfield(q,{'core','polar'}))&& ...
    isnumeric(q.core)&&isreal(q.core)&&isequal(size(q.core),size(g.core.X))&& ...
    isnumeric(q.polar)&&isreal(q.polar)&&isequal(size(q.polar),size(g.polar.X)), ...
    'ipm_cr:FEMSource','Native input is the already differentiated source on both complete transport arrays.');
end
function q=apply_source(m,s)
q=ipm_cr.mesh.sample(s,m);
end
function plan=query_plan(a,b,x,y,chart,side)
assert(isequal(size(x),size(y)),'ipm_cr:FEMQueryShape');shape=size(x);x=x(:);y=y(:);
if strcmp(chart,'core'),qa=x;qb=y;else,qa=log(hypot(x,y));qb=atan2(y,x);end
qa=roundoff_clamp(qa,a.breaks);qb=roundoff_clamp(qb,b.breaks);
[Ea,Da]=ipm_cr.poisson.fem.assembly('evaluate',a,qa,side);[Eb,Db]=ipm_cr.poisson.fem.assembly('evaluate',b,qb,side);
plan=struct('E',row_tensor(Ea,Eb,chart),'Da',row_tensor(Da,Eb,chart), ...
    'Db',row_tensor(Ea,Db,chart),'shape',shape,'x',x,'y',y,'chart',chart,'side',side);
end
function q=roundoff_clamp(q,interval)
lo=interval(1);hi=interval(end);tol=64*eps(max([1,abs(lo),abs(hi)]));
assert(all(q>=lo-tol&q<=hi+tol),'ipm_cr:FEMQueryOutside','Query is outside its FEM patch.');
q=max(lo,min(hi,q));
% Polar physical points make an r -> (x,y) -> hypot -> log round trip.
% Snap only a few coordinate ulps to an existing element face so its
% declared averaged derivative does not become an accidental one-sided one.
[distance,index]=min(abs(q-interval(:).'),[],2);
onFace=distance<=8*eps(max(1,abs(q)));q(onFace)=interval(index(onFace));
end
function M=row_tensor(A,B,chart)
if size(A,1)==0,M=sparse(0,size(A,2)*size(B,2));return;end
[ca,va]=row_lists(A);[cb,vb]=row_lists(B);n=size(A,1);na=size(A,2);nb=size(B,2);
ii={};jj={};vv={};count=0;
for a=1:size(ca,2),for b=1:size(cb,2)
    use=ca(:,a)>0&cb(:,b)>0&va(:,a)~=0&vb(:,b)~=0;if ~any(use),continue;end
    count=count+1;ii{count}=find(use); %#ok<AGROW>
    if strcmp(chart,'core'),jj{count}=cb(use,b)+(ca(use,a)-1)*nb;
    else,jj{count}=ca(use,a)+(cb(use,b)-1)*na;end %#ok<AGROW>
    vv{count}=va(use,a).*vb(use,b); %#ok<AGROW>
end,end
if count==0,M=sparse(n,na*nb);else,M=sparse(vertcat(ii{:}),vertcat(jj{:}),vertcat(vv{:}),n,na*nb);end
end
function [columns,values]=row_lists(A)
[c,r,v]=find(A');counts=accumarray(r,1,[size(A,1),1]);w=max([0;counts]);
columns=zeros(size(A,1),w);values=columns;counts(:)=0;
for j=1:numel(r),i=r(j);counts(i)=counts(i)+1;columns(i,counts(i))=c(j);values(i,counts(i))=v(j);end
end
function v=evaluate_plan(plan,U)
if isfield(plan,'native')&&plan.native
    if strcmp(plan.chart,'core')
        u1=-(plan.Db*U);u2=U*plan.Da';
    else
        dt=plan.Da*U;da=U*plan.Db';
        u1=-(plan.sinTheta.*dt+plan.cosTheta.*da)./plan.radius;
        u2=(plan.cosTheta.*dt-plan.sinTheta.*da)./plan.radius;
    end
    v=struct('psi',U,'u1',u1,'u2',u2);return
end
psi=plan.E*U(:);a=plan.Da*U(:);b=plan.Db*U(:);
if strcmp(plan.chart,'core'),u1=-b;u2=a;
else
    r=hypot(plan.x,plan.y);co=plan.x./r;si=plan.y./r;
    u1=-(si.*a+co.*b)./r;u2=(co.*a-si.*b)./r;
end
v=struct('psi',reshape(psi,plan.shape),'u1',reshape(u1,plan.shape),'u2',reshape(u2,plan.shape));
end
function value=maximum_abs(values)
if isempty(values),value=NaN;else,value=max(abs(values(:)));end
end
function yes=finite_flow(flow)
yes=all(isfinite(flow.core.psi),'all')&&all(isfinite(flow.polar.psi),'all')&& ...
    all(isfinite(flow.core.u1),'all')&&all(isfinite(flow.core.u2),'all')&& ...
    all(isfinite(flow.polar.u1),'all')&&all(isfinite(flow.polar.u2),'all');
end
