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
    end
    methods
        function obj=context(g,wall,opts)
            % Build every candidate before publishing the prepared context.
            initial='data/initial.json';if isfield(opts,'initial'),initial=opts.initial;end
            provenance=struct('source',ipm_cr.io.identity('source'), ...
                'authority',ipm_cr.io.identity('authority',ipm_cr.root(),initial));
            plan=ipm_cr.poisson.plan(g,opts.solverOptions);
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
        end
        function [flow,info]=velocity(obj,source,analyticSource,includePsi)
            if nargin<3,analyticSource=[];end
            if nargin<4,includePsi=true;end
            [flow,info]=ipm_cr.poisson.apply(source,obj.mesh,obj.pmap, ...
                analyticSource,includePsi,obj.calls==0);
            obj.calls=obj.calls+1;
        end
        function value=stats(obj)
            value=struct('rhsPoissonCalls',obj.calls,'poisson',obj.pmap.solver.Counters, ...
                'buildSeconds',obj.pmap.buildSeconds);
        end
    end
end
