classdef enriched < handle
    %ENRICHED Optional outer-corner Green enrichment of the unchanged FEM.
    % q=qP+(q-qP), where qP matches total-order-three jets at (+/-R,0).
    % The polynomial's exact finite-halfdisk Green field is global and shared
    % by both patches. Only the remainder is represented by FEM coefficients.
    % It is applied only to an explicitly declared compact source. Known
    % smooth exterior/Dirichlet problems pass through the unchanged raw FEM:
    % adding a finite-disk corner field there can create a singular remainder.
    properties (SetAccess=private)
        raw
        g
        c
        p
        coupling
        auditMaps
        config
        buildSeconds
        greenBuildSeconds
        greenResponseCacheBytes
    end
    properties (Dependent,SetAccess=private)
        Counters
    end
    properties (Access=private)
        responses
        selection
        coefficientIndex
        cornerFits=0
    end
    methods
        function obj=enriched(g,options)
            if nargin<2,options=struct();end
            timer=tic;obj.raw=ipm_cr.poisson.fem.solver(g,options);
            obj.g=g;obj.c=obj.raw.c;obj.p=obj.raw.p;
            obj.coupling=obj.raw.coupling;obj.auditMaps=obj.raw.auditMaps;
            obj.config=obj.raw.config;obj.config.cornerCorrection=true;
            phase=tic;c=obj.c;p=obj.p;nc=numel(g.core.X);np=numel(g.polar.X);
            x=[g.core.X(:);g.polar.X(:);c.X(c.traceIndex);p.X(1,p.jj)'; ...
                obj.auditMaps.x(:);p.X(end,:)'];
            y=[g.core.Y(:);g.polar.Y(:);c.Y(c.traceIndex);p.Y(1,p.jj)'; ...
                obj.auditMaps.y(:);p.Y(end,:)'];
            sizes=[nc,np,numel(c.traceIndex),numel(p.jj),numel(obj.auditMaps.x),p.nt];
            names={'core','polar','coreTrace','innerTrace','overlap','outerTrace'};
            first=0;for k=1:numel(names)
                obj.selection.(names{k})=first+(1:sizes(k));first=first+sizes(k);
            end
            index=[];for b=0:3,index=[index,(1:8-2*b)+8*b];end %#ok<AGROW>
            obj.coefficientIndex=index;
            response=struct('psi',zeros(numel(x),20),'u1',zeros(numel(x),20),'u2',zeros(numel(x),20));
            for j=1:20
                coefficient=zeros(8,4);coefficient(index(j))=1;
                v=ipm_cr.poisson.fem.greenPolynomial(coefficient,x,y,g.polar.r(end));
                for field={'psi','u1','u2'},name=field{1};response.(name)(:,j)=v.(name)(:);end
            end
            obj.responses=response;
            obj.greenResponseCacheBytes=8*3*numel(x)*20;
            obj.greenBuildSeconds=toc(phase);obj.buildSeconds=toc(timer);
        end

        function [flow,info]=solve(obj,source,boundary,analyticSource)
            if nargin<4,analyticSource=[];end
            assert(nargin>=3,'ipm_cr:FEMExteriorData','Explicit exterior data are mandatory.');
            timer=tic;originalBoundary=obj.raw.boundary_data(boundary);
            if ~strcmp(boundary.exteriorStatus,'declared_compact_source')
                [flow,info]=obj.raw.solve(source,boundary,analyticSource);
                info.cornerCorrection=false;info.cornerCorrectionSkipped=true;
                info.cornerCorrectionReason='finite-disk corner enrichment requires declared_compact_source; original source and exterior data retained';
                info.polynomialCoefficientCount=0;info.cornerFitSeconds=0;info.greenApplySeconds=0;
                info.buildSeconds=obj.buildSeconds;info.greenBuildSeconds=obj.greenBuildSeconds;
                info.greenResponseCacheBytes=obj.greenResponseCacheBytes;info.totalSeconds=toc(timer);
                return
            end
            phase=tic;
            [coefficient,fit]=ipm_cr.poisson.fem.corner(source,obj.g,analyticSource);
            obj.cornerFits=obj.cornerFits+1;R=obj.g.polar.r(end);
            zeroCorrection=~any(coefficient(:));
            if zeroCorrection
                offset=analyticSource;
            elseif isempty(analyticSource)
                offset=@(x,y)-ipm_cr.poisson.fem.sourcePolynomial(coefficient,x,y,R);
            else
                offset=@(x,y)analyticSource(x,y)-ipm_cr.poisson.fem.sourcePolynomial(coefficient,x,y,R);
            end
            if zeroCorrection
                zero=zeros(size(obj.responses.psi,1),1);
                green=struct('psi',zero,'u1',zero,'u2',zero);
            else
                values=coefficient(obj.coefficientIndex);
                green=struct();for field={'psi','u1','u2'}
                    name=field{1};green.(name)=obj.responses.(name)*values(:);
                end
            end
            remainderBoundary=boundary;
            if strcmp(boundary.kind,'green_dirichlet')&&~zeroCorrection
                remainderBoundary.value=originalBoundary.values-green.psi(obj.selection.outerTrace)';
            end
            fitSeconds=toc(phase);
            [flow,info]=obj.raw.solve(source,remainderBoundary,offset);
            phase=tic;
            if ~zeroCorrection
                for patch={'core','polar'}
                    key=patch{1};indices=obj.selection.(key);shape=size(obj.g.(key).X);
                    for field={'psi','u1','u2'}
                        name=field{1};flow.(key).(name)=flow.(key).(name)+reshape(green.(name)(indices),shape);
                    end
                end
            end
            flow.enrichment=struct('coefficients',coefficient,'radius',R, ...
                'source','normalized Cartesian polynomial, total-order-three outer-corner jets', ...
                'field','exact global finite-halfdisk Green field; zero exterior source');
            flow.representation=[flow.representation,' + exact global polynomial Green field'];
            cc=flow.coefficients.core;pa=flow.coefficients.polar;
            gCore=green.psi(obj.selection.coreTrace);gInner=green.psi(obj.selection.innerTrace);
            fullTrace=max(norm((cc(obj.c.traceIndex)+gCore)-(obj.coupling.RCA*pa(:)+gCore),inf), ...
                norm((obj.coupling.RAC*cc(:)+gInner)-(pa(1,obj.p.jj)'+gInner),inf));
            fullEssential=max([max(abs(flow.core.psi(1,:))),max(abs(flow.polar.psi(:,[1,end])),[],'all')]);
            outerEssential=NaN;
            if strcmp(boundary.kind,'green_dirichlet')
                outerEssential=norm(pa(end,:)+green.psi(obj.selection.outerTrace)'-originalBoundary.values,inf);
                fullEssential=max(fullEssential,outerEssential);
            end
            [a,b]=obj.raw.overlap(flow);indices=obj.selection.overlap;
            for field={'psi','u1','u2'}
                name=field{1};v=reshape(green.(name)(indices),size(obj.auditMaps.x));
                a.(name)=a.(name)+v;b.(name)=b.(name)+v;
            end
            info.remainderWeakBackwardError=info.weakBackwardError;
            info.remainderTraceAbsolute=info.traceAbsolute;
            info.remainderEssentialBoundaryAbsolute=info.essentialBoundaryAbsolute;
            info.traceAbsolute=max(info.traceAbsolute,fullTrace);
            info.fullSharedTraceAbsolute=fullTrace;
            info.essentialBoundaryAbsolute=max(info.essentialBoundaryAbsolute,fullEssential);
            info.outerEssentialBoundaryAbsolute=outerEssential;
            info.overlapValueMismatch=maximum_abs(a.psi(:)-b.psi(:));
            info.overlapVelocityMismatch=maximum_abs([a.u1(:)-b.u1(:);a.u2(:)-b.u2(:)]);
            info.cornerCorrection=true;info.cornerCorrectionSkipped=false;
            info.cornerCorrectionReason='declared compact finite-disk source; preserve nonzero outer-wall cutoff jets';
            info.cornerFit=fit;
            info.zeroCornerCorrection=zeroCorrection;
            info.polynomialCoefficientCount=20;
            info.originalAnalyticSourceOffset=~isempty(analyticSource);
            info.sourceCorrectionEvaluation='native source unchanged; original analytic source minus qP at Gauss points';
            if zeroCorrection
                info.sourceCorrectionEvaluation='all corner coefficients zero; original source and analytic callback passed directly to FEM';
            end
            info.discretizationContract='overlap_C0_Qp_galerkin_plus_global_polynomial_Green_v1';
            info.diagnosticScope='FEM remainder weak rows plus full shared traces and wall/outer boundary measured only; no experiment acceptance threshold';
            info.weakResidualScope='assembled FEM equations for q-qP only; full field also includes analytical Green(qP)';
            info.cornerFitSeconds=fitSeconds;info.greenApplySeconds=toc(phase);
            info.sourceSeconds=info.sourceSeconds+fitSeconds;
            info.returnSeconds=info.returnSeconds+info.greenApplySeconds;
            info.kernelSeconds=info.sourceSeconds+info.solveSeconds+info.returnSeconds;
            info.totalSeconds=toc(timer);info.buildSeconds=obj.buildSeconds;
            info.greenBuildSeconds=obj.greenBuildSeconds;
            info.greenResponseCacheBytes=obj.greenResponseCacheBytes;
            info.finiteOutput=all(isfinite(flow.core.psi),'all')&&all(isfinite(flow.polar.psi),'all')&& ...
                all(isfinite(flow.core.u1),'all')&&all(isfinite(flow.core.u2),'all')&& ...
                all(isfinite(flow.polar.u1),'all')&&all(isfinite(flow.polar.u2),'all');
        end

        function value=evaluate(obj,flow,x,y,chart,side)
            if nargin<6,side='average';end
            value=obj.raw.evaluate(flow,x,y,chart,side);
            if ~isfield(flow,'enrichment'),return;end
            g=ipm_cr.poisson.fem.greenPolynomial(flow.enrichment.coefficients,x,y,flow.enrichment.radius);
            for field={'psi','u1','u2'},name=field{1};value.(name)=value.(name)+g.(name);end
        end

        function [a,b]=overlap(obj,flow)
            a=obj.evaluate(flow,obj.auditMaps.x,obj.auditMaps.y,'core');
            b=obj.evaluate(flow,obj.auditMaps.x,obj.auditMaps.y,'polar');
        end

        function value=get.Counters(obj)
            value=obj.raw.Counters;value.greenResponseCacheBuilds=1;
            value.cornerFits=obj.cornerFits;
        end
    end
end
function value=maximum_abs(values)
if isempty(values),value=NaN;else,value=max(abs(values(:)));end
end
