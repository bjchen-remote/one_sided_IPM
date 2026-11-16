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
        wallProbe
        identity
        spline
    end
    properties (Access=private)
        pmap
        calls=0
    end
    methods
        function obj=context(g,wall,opts)
            provenance=struct('source',ipm_cr.io.identity('source'), ...
                'authority',ipm_cr.io.identity('authority',ipm_cr.io.identity('root'),opts.initial));
            fem=ipm_cr.space.poisson(g);
            obj.spline=struct('x',fem.c.xAxis,'y',fem.c.yAxis, ...
                't',fem.p.tAxis,'theta',fem.p.thetaAxis);
            r=unique([0;abs(wall.x(:))]);X=[r,-r];Y=zeros(size(X));
            empty=struct('core',zeros(size(g.core.X)),'polar',zeros(size(g.polar.X)));
            [~,map]=ipm_cr.mesh.interpolate.sample(empty,g,X,Y);
            assert(map.allValid,'ipm_cr:WallCoverage','Every original wall point needs an active owner.');
            obj.mesh=g;obj.wall=wall;obj.opts=opts;
            assert(isequaln(provenance.source,ipm_cr.io.identity('source'))&& ...
                isequaln(provenance.authority,ipm_cr.io.identity('authority',ipm_cr.io.identity('root'),opts.initial)), ...
                'ipm_cr:ContextIdentity','Source or initial datum changed while preparing the context.');
            obj.identity=provenance;
            obj.key=context_key(g,wall,opts,provenance);obj.id=tempname;
            obj.pmap=fem;obj.wallR=r;obj.wallMap=map;
            obj.wallProbe=ipm_cr.mesh.interpolate.wallPlan(g,map);
        end
        function [flow,info]=velocity(obj,source,includePsi)
            if nargin<3,includePsi=true;end
            [flow,info]=obj.pmap.solve(source);
            if ~includePsi,flow.core.psi=[];flow.polar.psi=[];end
            info.definition='minus_laplacian_psi_equals_native_complete_source_reconstruction';
            info.exteriorModel='finite_source_disk_zero_exterior';
            info.linearBackwardError=info.weakBackwardError;
            info.traceResidual=info.traceAbsolute;
            info.normalFluxMismatch=NaN;info.strongResidual=NaN;
            info.factorCacheRebuilt=obj.calls==0;
            info.returnIncludesPsi=includePsi;
            info.wallNormalVelocityInfinity=max([max(abs(flow.core.u2(obj.mesh.core.Y==0))), ...
                max(abs(flow.polar.u2(:,[1,end])),[],'all')]);
            obj.calls=obj.calls+1;
        end
        function value=stats(obj)
            value=struct('rhsPoissonCalls',obj.calls,'poisson',obj.pmap.Counters, ...
                'buildSeconds',obj.pmap.buildSeconds);
            value.backend='bspline4';
        end
    end
end

function key=context_key(g,wall,opts,identity)
%KEY Mathematical identity for an RHS; time-step policy is independent.
names={'cfl','maxDt'};
opts=rmfield(opts,names(isfield(opts,names)));
key=struct('geometry',g.geometryKey,'wall',wall,'method',opts,'identity',identity);
end
