function result=apply(prepared,vx,vy,request)
%APPLY Apply one frozen transport reconstruction using named requests.
% The default always evaluates the complete active-patch transport CFL.
% Only private gauge probes should explicitly request cfl=false.
if nargin<4,request=struct();end
assert(isa(prepared,'function_handle'),'ipm_cr:TransportPrepared', ...
    'Supply the current RHS reconstruction returned by op.prepare.');
result=prepared(vx,vy,request);
end
