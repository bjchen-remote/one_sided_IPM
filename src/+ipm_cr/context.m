classdef context < handle
    %CONTEXT One prepared simulation; mesh remains a value struct.
    % Fixed arrays are read-only here. Only owned solve counters change.
    properties (SetAccess=private)
        mesh
        wall
        opts
        id
        key
        wallR
        wallMap
        identity
    end
    properties (Access=private)
        pmap
        calls=0
        backend='spectral'
    end
    methods
        function obj=context(g,wall,opts)
            % Build every candidate before publishing the prepared context.
            if ~isfield(opts,'poissonBackend'),opts.poissonBackend='spectral';end
            if ~isfield(opts,'fdOptions'),opts.fdOptions=struct();end
            if ~isfield(opts,'femOptions'),opts.femOptions=struct();end
            assert(((ischar(opts.poissonBackend)&&isrow(opts.poissonBackend))|| ...
                (isstring(opts.poissonBackend)&&isscalar(opts.poissonBackend)))&& ...
                any(strcmp(opts.poissonBackend,{'spectral','same_mesh_fd','fem'})), ...
                'ipm_cr:PoissonBackend','Use spectral, same_mesh_fd, or the experimental fem backend.');
            opts.poissonBackend=char(opts.poissonBackend);
            assert(isstruct(opts.fdOptions)&&isscalar(opts.fdOptions), ...
                'ipm_cr:FDOptions','Use a scalar same-mesh FD options structure.');
            [opts.femOptions,femBoundary]=ipm_cr.poisson.fem.options(opts.femOptions);
            initial='data/initial.json';if isfield(opts,'initial'),initial=opts.initial;end
            provenance=struct('source',ipm_cr.io.identity('source'), ...
                'authority',ipm_cr.io.identity('authority',ipm_cr.root(),initial));
            if strcmp(opts.poissonBackend,'spectral')
                plan=ipm_cr.poisson.plan(g,opts.solverOptions);
            elseif strcmp(opts.poissonBackend,'same_mesh_fd')
                fd=ipm_cr.poisson.fd.solver(g,opts.fdOptions);
                plan=struct('solver',fd,'buildSeconds',fd.buildSeconds);
            else
                options=rmfield(opts.femOptions,{'exteriorModel','cornerCorrection'});
                if opts.femOptions.cornerCorrection
                    fem=ipm_cr.poisson.fem.enriched(g,options);
                else
                    fem=ipm_cr.poisson.fem.solver(g,options);
                end
                plan=struct('solver',fem,'buildSeconds',fem.buildSeconds,'boundary',femBoundary);
            end
            r=unique([0;abs(wall.x(:))]);X=[r,-r];Y=zeros(size(X));
            empty=struct('core',zeros(size(g.core.X)),'polar',zeros(size(g.polar.X)));
            [~,map]=ipm_cr.mesh.sample(empty,g,X,Y);
            assert(map.allValid,'ipm_cr:WallCoverage','Every original wall point needs an active owner.');
            obj.mesh=g;obj.wall=wall;obj.opts=opts;
            assert(isequaln(provenance.source,ipm_cr.io.identity('source'))&& ...
                isequaln(provenance.authority,ipm_cr.io.identity('authority',ipm_cr.root(),initial)), ...
                'ipm_cr:ContextIdentity','Source or initial datum changed while preparing the context.');
            obj.identity=provenance;
            obj.key=ipm_cr.internal.key(g,wall,opts,provenance);obj.id=tempname;
            obj.pmap=plan;obj.wallR=r;obj.wallMap=map;
            obj.backend=opts.poissonBackend;
        end
        function [flow,info]=velocity(obj,source,analyticSource,includePsi)
            if nargin<3,analyticSource=[];end
            if nargin<4,includePsi=true;end
            if strcmp(obj.backend,'spectral')
                [flow,info]=ipm_cr.poisson.apply(source,obj.mesh,obj.pmap, ...
                    analyticSource,includePsi,obj.calls==0);
            elseif strcmp(obj.backend,'same_mesh_fd')
                validateattributes(includePsi,{'logical'},{'scalar'});
                [flow,info]=obj.pmap.solver.solve(source,includePsi,[],analyticSource);
                if ~includePsi,flow.core.psi=[];flow.polar.psi=[];end
                info.poissonBackend='same_mesh_fd';
                info.definition='minus_laplacian_psi_equals_physical_source_on_unique_active_transport_owner';
                info.discretizationContract='experimental_same_mesh_fd7_trace_schur_v1';
                info.factorCacheRebuilt=obj.calls==0;
                info.sourceMapCacheReused=obj.calls>0;info.returnMapCacheReused=obj.calls>0;
                info.explicitGeometryBuildSeconds=obj.pmap.buildSeconds;
                info.sourceInterpolationSeconds=info.sourceSeconds;
                info.returnInterpolationSeconds=info.returnSeconds;
                info.returnOperation='same-node discrete streamfunction derivatives; no spectral return interpolation';
                info.totalAdapterSolveSeconds=info.totalSeconds;info.totalCallSeconds=info.totalSeconds;
                info.returnIncludesPsi=includePsi;
                info.maximumNativeReconstructionDefectScaledBySource=NaN;
                info.nativeReconstructionDefectSourceMapsCached=false;
                info.nativeReconstructionDefectIsIndependentContinuumResidual=false;
                info.sameSemiDiscreteOperatorClaimed=false;info.iterativeLinearSolverUsed=false;
                info.strongResidualIsAlgebraicResidual=false;info.strongResidualIsAccuracyCertificate=false;
                info.wallNormalVelocityInfinity=max([max(abs(flow.core.u2(obj.mesh.core.Y==0))), ...
                    max(abs(flow.polar.u2(:,[1,end])),[],'all')]);
            else
                validateattributes(includePsi,{'logical'},{'scalar'});
                sourceMapsBefore=obj.pmap.solver.Counters.sourceMapBuilds;
                [flow,info]=obj.pmap.solver.solve(source,obj.pmap.boundary,analyticSource);
                if ~includePsi,flow.core.psi=[];flow.polar.psi=[];end
                info.poissonBackend='fem';
                info.definition='minus_laplacian_psi_equals_native_residual_reconstruction_plus_analytic_source';
                info.exteriorModel=obj.opts.femOptions.exteriorModel;
                info.linearBackwardError=info.weakBackwardError;
                info.traceResidual=info.traceAbsolute;
                info.normalFluxMismatch=NaN;info.strongResidual=NaN;
                info.factorCacheRebuilt=obj.calls==0;
                info.sourceMapCacheReused=sourceMapsBefore>0&&~isa(source,'function_handle');
                info.returnMapCacheReused=true;
                info.explicitGeometryBuildSeconds=obj.pmap.buildSeconds;
                info.sourceInterpolationSeconds=info.sourceSeconds;
                info.returnInterpolationSeconds=info.returnSeconds;
                info.returnOperation='native-node FEM coefficients and analytical element derivatives; no independent solution grid';
                info.sameNativeGrid=true;info.independentPoissonGridUsed=false;
                info.quadraturePointsAreUnknowns=false;
                if obj.opts.femOptions.cornerCorrection
                    info.returnOperation='native-node FEM remainder plus one global finite-disk polynomial Green field';
                end
                info.totalAdapterSolveSeconds=info.totalSeconds;info.totalCallSeconds=info.totalSeconds;
                info.returnIncludesPsi=includePsi;
                info.sameSemiDiscreteOperatorClaimed=false;info.iterativeLinearSolverUsed=false;
                info.strongResidualIsAlgebraicResidual=false;info.strongResidualIsAccuracyCertificate=false;
                info.wallNormalVelocityInfinity=max([max(abs(flow.core.u2(obj.mesh.core.Y==0))), ...
                    max(abs(flow.polar.u2(:,[1,end])),[],'all')]);
            end
            obj.calls=obj.calls+1;
        end
        function value=stats(obj)
            value=struct('rhsPoissonCalls',obj.calls,'poisson',obj.pmap.solver.Counters, ...
                'buildSeconds',obj.pmap.buildSeconds);
            value.backend=obj.backend;
        end
    end
end
