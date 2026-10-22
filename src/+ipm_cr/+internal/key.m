function key=key(g,wall,opts)
%KEY Mathematical identity for an RHS; time-step policy is independent.
names={'cfl','maxDt','minDt','timeIntegrator'};
opts=rmfield(opts,names(isfield(opts,names)));
key=struct('geometry',g.geometryKey,'wall',wall,'method',opts);
end
