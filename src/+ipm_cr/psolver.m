classdef psolver < handle
    %SOLVER Own geometry, cached patch factors and overlap Schur complement.
    % A source can be a continuous function or struct('core',qC,'polar',qP).
    % Reconfiguration is transactional; old solutions cannot use new geometry.
    properties (SetAccess=private)
        Config
        Core
        Polar
        Coupling
        Counters=struct('builds',0,'solves',0,'factorReuses',0,'accepted',0,'rejected',0, ...
            'linearSuperpositions',0,'superpositionFallbacks',0)
        Id
        BuildSeconds
    end
    methods
        function obj=psolver(options)
            if nargin<1,options=struct();end
            obj.reconfigure(options);
        end
        function reconfigure(obj,options)
            if isempty(obj.Config),next=ipm_cr.pconfig(options);
            else
                next=obj.Config;names=fieldnames(options);
                for k=1:numel(names),next.(names{k})=options.(names{k});end
                next=ipm_cr.pconfig(next);
                if isequaln(next,obj.Config),return;end
            end
            timer=tic;
            core=ipm_cr.poisson.core.build(next);polar=ipm_cr.poisson.polar.build(next);
            nm=numel(polar.interior);B=zeros(numel(core.traceIndex),nm);T=zeros(nm);
            coreResponseMatrix=zeros(numel(core.X),nm);
            ix=polar.X(1,polar.interior)';iy=polar.Y(1,polar.interior)';
            for j=1:nm
                harmonic=struct('modal',polar.harmonicModes.*polar.Vinv(:,j).');
                v=ipm_cr.poisson.polar.evaluate(polar,harmonic,core.traceX,core.traceY);
                assert(all(isfinite(v.psi)),'ipm_cr:CouplingGeometry','Nonfinite polar trace map.');
                B(:,j)=v.psi;
                h=ipm_cr.poisson.core.solve(core,zeros(size(core.X)),v.psi);
                % This homogeneous response is already needed for T. Keep
                % its full field for subsequent Schur linear combinations.
                coreResponseMatrix(:,j)=h.psi(:);
                v=ipm_cr.poisson.core.evaluate(core,h,ix,iy);T(:,j)=v.psi;
            end
            S=eye(nm)-T;
            assert(rcond(S)>1e-12,'ipm_cr:InterfaceCondition','Ill-conditioned overlap Schur system.');
            coupling=struct('boundaryMap',B,'matrix',S,'factor',decomposition(S,'lu'), ...
                'condition',cond(S),'innerX',ix,'innerY',iy, ...
                'coreResponseMatrix',coreResponseMatrix);
            % Commit only after all construction and validation have succeeded.
            obj.Config=next;obj.Core=core;obj.Polar=polar;obj.Coupling=coupling;
            obj.Id=tempname;obj.BuildSeconds=toc(timer);
            obj.Counters.builds=obj.Counters.builds+1;
        end
        function [s,info]=solve(obj,source,outer)
            c=obj.Core;p=obj.Polar;g=obj.Config;timer=tic;
            if nargin<3,outer=zeros(1,p.nt);end
            if isa(outer,'function_handle'),outer=outer(p.X(end,:),p.Y(end,:));end
            if strcmp(g.outerBoundary,'dirichlet')&&nargin<3
                error('ipm_cr:OuterBoundary','Supply explicit outer Dirichlet data.');
            end
            [qc,qp,sourceFcn]=obj.sample(source);
            obj.Counters.solves=obj.Counters.solves+1;
            obj.Counters.factorReuses=obj.Counters.factorReuses+1;
            p0=ipm_cr.poisson.polar.solve(p,qp,[],outer);
            boundary=ipm_cr.poisson.polar.evaluate(p,p0,c.traceX,c.traceY);
            c0=ipm_cr.poisson.core.solve(c,qc,boundary.psi);
            trace=ipm_cr.poisson.core.evaluate(c,c0,obj.Coupling.innerX,obj.Coupling.innerY);
            lambda=obj.Coupling.factor\trace.psi;
            % Refine against the original small Schur equation.
            for k=1:2
                defect=trace.psi-obj.Coupling.matrix*lambda;
                candidate=lambda+obj.Coupling.factor\defect;
                if norm(trace.psi-obj.Coupling.matrix*candidate,inf)>=norm(defect,inf),break;end
                lambda=candidate;
            end
            inner=zeros(1,p.nt);inner(p.interior)=lambda;
            % Same discrete equations: source/outer particular solution
            % plus the cached homogeneous response to the Schur trace.
            ps=ipm_cr.poisson.polar.withInnerTrace(p,p0,inner);
            cs=struct('psi',c0.psi+reshape(obj.Coupling.coreResponseMatrix*lambda,size(c.X)));
            obj.Counters.linearSuperpositions=obj.Counters.linearSuperpositions+1;
            s=struct('core',cs,'polar',ps,'solverId',obj.Id);
            info=obj.audit(s,qc,qp,sourceFcn,trace.psi,lambda);
            candidateBackward=info.linearBackwardError;fallback=false;
            if ~isfinite(candidateBackward)||candidateBackward>g.linearTolerance
                % Superposition can lose digits through cancellation. Fall
                % back to the ORIGINAL direct solves and ORIGINAL audit;
                % no tolerance or spatial acceptance condition is weakened.
                fallback=true;
                obj.Counters.superpositionFallbacks=obj.Counters.superpositionFallbacks+1;
                ps=ipm_cr.poisson.polar.solve(p,qp,inner,outer);
                cs=ipm_cr.poisson.core.solve(c,qc,boundary.psi+obj.Coupling.boundaryMap*lambda,true);
                s=struct('core',cs,'polar',ps,'solverId',obj.Id);
                info=obj.audit(s,qc,qp,sourceFcn,trace.psi,lambda);
                info.coreOriginalPdeRefinement=cs.refinement;
            end
            info.linearSuperposition=struct('attempted',true,'used',~fallback, ...
                'directFallbackUsed',fallback,'candidateLinearBackwardError',candidateBackward, ...
                'sameDiscreteSystems',true,'cachedPolarHomogeneousModes',true, ...
                'cachedCoreResponseMatrix',true,'coreResponseCount',size(obj.Coupling.coreResponseMatrix,2), ...
                'polarDirectSolvesThisRhs',1+double(fallback), ...
                'coreDirectSolvesThisRhs',1+double(fallback), ...
                'accuracyGatesRelaxed',false);
            info.solveSeconds=toc(timer);info.buildSeconds=obj.BuildSeconds;
            info.factorBuildCount=obj.Counters.builds;
            info.factorReuseCount=obj.Counters.factorReuses;
            info.outerBoundary=g.outerBoundary;info.outerRadius=g.rOuter;
            info.outerSourceTailCertified=false;
            if info.accepted,obj.Counters.accepted=obj.Counters.accepted+1;
            else,obj.Counters.rejected=obj.Counters.rejected+1;end
            s.info=info;
            if g.requireAccuracy&&~info.accepted
                error('ipm_cr:AccuracyGate', ...
                    'Spatial audit failed: linear %.3g, trace %.3g, value %.3g, velocity %.3g, strong %.3g. Refine geometry.', ...
                    info.linearBackwardError,info.traceResidual,info.overlapValueMismatch, ...
                    info.overlapVelocityMismatch,info.strongResidual);
            end
        end
        function [s,info]=velocity(obj,rho,outer)
            %VELOCITY Physical Darcy velocity from density on both native grids.
            [rc,rp]=obj.sample(rho);c=obj.Core;p=obj.Polar;
            % Differentiate a centered density to preserve additive invariance.
            origin=ipm_cr.internal.bary(c.y,0)*rc*ipm_cr.internal.bary(c.x,0)';
            rc=rc-origin;rp=rp-origin;
            qc=rc*c.Dx';
            qp=((p.Dt*rp).*cos(p.theta)-(rp*p.Dtheta').*sin(p.theta))./p.r;
            if nargin<3,[s,info]=obj.solve(struct('core',qc,'polar',qp));
            else,[s,info]=obj.solve(struct('core',qc,'polar',qp),outer);end
            info.sourceOrigin='native spectral Cartesian derivative of centered density';
            info.sourceInterpolationAccuracyCertified=false;s.info=info;
        end
        function v=evaluate(obj,s,x,y)
            obj.validateSolution(s);
            assert(isequal(size(x),size(y))&&isreal(x)&&isreal(y), ...
                'ipm_cr:Queries','Use matching real query arrays.');
            shape=size(x);x=x(:);y=y(:);g=obj.Config;
            valid=isfinite(x)&isfinite(y)&y>=0&hypot(x,y)<=g.rOuter*(1+16*eps);
            % Ownership is unique: use the Cartesian interpolant in its core.
            core=valid&abs(x)<=g.coreHalfWidth&y<=g.coreHeight;
            polar=valid&~core;
            fields={'psi','u1','u2','laplacian'};
            for k=1:numel(fields),v.(fields{k})=nan(size(x));end
            if any(core)
                a=ipm_cr.poisson.core.evaluate(obj.Core,s.core,x(core),y(core));
                for k=1:numel(fields),v.(fields{k})(core)=a.(fields{k});end
            end
            if any(polar)
                a=ipm_cr.poisson.polar.evaluate(obj.Polar,s.polar,x(polar),y(polar));
                for k=1:numel(fields),v.(fields{k})(polar)=a.(fields{k});end
            end
            for k=1:numel(fields),v.(fields{k})=reshape(v.(fields{k}),shape);end
            v.owner=reshape(double(core)+2*double(polar),shape);
        end
    end
    methods (Access=private)
        function validateSolution(obj,s)
            assert(isstruct(s)&&isfield(s,'solverId')&&strcmp(s.solverId,obj.Id), ...
                'ipm_cr:StaleSolution','The solution belongs to different or reconfigured geometry.');
        end
        function [a,b,callback]=sample(obj,input)
            c=obj.Core;p=obj.Polar;callback=[];
            if isa(input,'function_handle')
                callback=input;a=input(c.X,c.Y);b=input(p.X,p.Y);
            else
                assert(isstruct(input)&&isfield(input,'core')&&isfield(input,'polar'), ...
                    'ipm_cr:Source','Use a continuous function or core/polar arrays.');
                a=input.core;b=input.polar;
            end
            if isscalar(a),a=repmat(a,size(c.X));end
            if isscalar(b),b=repmat(b,size(p.X));end
            assert(isequal(size(a),size(c.X))&&isequal(size(b),size(p.X))&& ...
                isreal(a)&&isreal(b)&&all(isfinite(a),'all')&&all(isfinite(b),'all'), ...
                'ipm_cr:Source','Native arrays must be finite, real and correctly sized.');
        end
        function info=audit(obj,s,qc,qp,callback,rhs,lambda)
            c=obj.Core;p=obj.Polar;g=obj.Config;
            scale=max([max(abs(qc),[],'all'),max(abs(qp),[],'all'),realmin]);
            amplitude=max([max(abs(s.core.psi),[],'all'),max(abs(s.polar.psi),[],'all'),realmin]);
            residual=-c.Dyy*s.core.psi-s.core.psi*c.Dxx'-qc;
            denom=abs(c.Dyy)*abs(s.core.psi)+abs(s.core.psi)*abs(c.Dxx')+abs(qc);
            backward=max(abs(residual(c.iy,c.ix))./max(denom(c.iy,c.ix),realmin),[],'all');
            for k=1:numel(p.mu)
                A=p.operators{k};u=s.polar.modal(:,k);b=s.polar.rhs(:,k);
                backward=max(backward,max(abs(A*u-b)./max(abs(A)*abs(u)+abs(b),realmin)));
            end
            algebraicTrace=norm(obj.Coupling.matrix*lambda-rhs,inf)/max(amplitude,realmin);
            actualInner=ipm_cr.poisson.core.evaluate(c,s.core,obj.Coupling.innerX,obj.Coupling.innerY);
            actualOuter=ipm_cr.poisson.polar.evaluate(p,s.polar,c.traceX,c.traceY);
            traceResidual=max([norm(actualInner.psi-lambda,inf), ...
                norm(actualOuter.psi-s.core.psi(c.traceIndex),inf)])/max(amplitude,realmin);
            traceResidual=max(traceResidual,algebraicTrace);
            % Independent non-nodal points in the overlap and on core edges.
            theta=pi*((0:g.auditPoints-1)+.37)/g.auditPoints;
            radius=linspace(g.rInner,.98*min(g.coreHalfWidth,g.coreHeight),5)';
            x=radius*cos(theta);y=radius*sin(theta);x=x(:);y=y(:);
            nx=x./hypot(x,y);ny=y./hypot(x,y);
            z=linspace(0,1,g.auditPoints+2)';z=z(2:end-1);
            x=[x;-g.coreHalfWidth*ones(size(z));g.coreHalfWidth*ones(size(z));g.coreHalfWidth*(2*z-1)];
            y=[y;g.coreHeight*z;g.coreHeight*z;g.coreHeight*ones(size(z))];
            nx=[nx;-ones(size(z));ones(size(z));zeros(size(z))];
            ny=[ny;zeros(size(z));zeros(size(z));ones(size(z))];
            keep=hypot(x,y)>=g.rInner;x=x(keep);y=y(keep);nx=nx(keep);ny=ny(keep);
            a=ipm_cr.poisson.core.evaluate(c,s.core,x,y);b=ipm_cr.poisson.polar.evaluate(p,s.polar,x,y);
            valueMismatch=max(abs(a.psi-b.psi));
            velocityMismatch=max([abs(a.u1-b.u1);abs(a.u2-b.u2)]);
            % Circle radial normals and actual left/right/top core normals.
            flux=max(abs(-ny.*(a.u1-b.u1)+nx.*(a.u2-b.u2)));
            velocityScale=max([abs(a.u1);abs(a.u2);abs(b.u1);abs(b.u2);realmin]);
            % A unit floor makes tolerances mixed absolute/relative and keeps
            % the zero solution well-defined; report physical mismatches too.
            valueScale=max(1,amplitude);vScale=max(1,velocityScale);qScale=max(1,scale);
            strong=NaN;
            if ~isempty(callback)
                source=ipm_cr.internal.auditSource(callback,x,y);
                strong=max([abs(-a.laplacian-source);abs(-b.laplacian-source)])/qScale;
                z=((0:24)+.29)/25;
                [xx,yy]=meshgrid(g.coreHalfWidth*(2*z-1),g.coreHeight*z);
                v=ipm_cr.poisson.core.evaluate(c,s.core,xx(:),yy(:));
                f=ipm_cr.internal.auditSource(callback,xx(:),yy(:));
                strong=max(strong,max(abs(-v.laplacian-f),[],'all')/qScale);
                % Audit every log-r panel, including source-bearing exterior.
                mid=(p.logEdges(1:end-1)+p.logEdges(2:end))/2;
                half=diff(p.logEdges)/2;
                ts=reshape(mid+[-sqrt(3/5);0;sqrt(3/5)]*half,1,[]);
                [tt,aa]=ndgrid(ts,theta);xx=exp(tt).*cos(aa);yy=exp(tt).*sin(aa);
                v=ipm_cr.poisson.polar.evaluate(p,s.polar,xx,yy);f=ipm_cr.internal.auditSource(callback,xx,yy);
                strong=max(strong,max(abs(-v.laplacian-f),[],'all')/qScale);
            end
            accepted=isfinite(backward)&&backward<=g.linearTolerance&& ...
                traceResidual<=g.traceTolerance&&valueMismatch/valueScale<=g.overlapValueTolerance&& ...
                velocityMismatch/vScale<=g.overlapFluxTolerance&& ...
                flux/vScale<=g.overlapFluxTolerance&& ...
                (isempty(callback)||(isfinite(strong)&&strong<=g.strongTolerance));
            info=struct('accepted',accepted,'linearBackwardError',backward, ...
                'traceResidual',traceResidual,'schurResidual',algebraicTrace, ...
                'overlapValueMismatch',valueMismatch, ...
                'overlapVelocityMismatch',velocityMismatch,'normalFluxMismatch',flux, ...
                'scaledOverlapValueMismatch',valueMismatch/valueScale, ...
                'scaledOverlapVelocityMismatch',velocityMismatch/vScale, ...
                'scaledNormalFluxMismatch',flux/vScale, ...
                'strongResidual',strong,'strongResidualAvailable',~isempty(callback), ...
                'sourceScale',scale,'interfaceCondition',obj.Coupling.condition, ...
                'pde','-Delta psi=source; u=(-psi_y,psi_x)', ...
                'accuracyCertifiedByManufacturedSolution',false, ...
                'coreEigenvectorCondition',c.eigenvectorCondition);
        end
    end
end
